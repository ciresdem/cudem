### uncertainties.py
##
## Copyright (c) 2020 - 2023 Regents of the University of Colorado
##
## Permission is hereby granted, free of charge, to any person obtaining a copy 
## of this software and associated documentation files (the "Software"), to deal 
## in the Software without restriction, including without limitation the rights 
## to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies 
## of the Software, and to permit persons to whom the Software is furnished to do so, 
## subject to the following conditions:
##
## The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
##
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
## INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR 
## PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE 
## FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
## ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
##
### Commentary:
##
## Generate uncertainty raster(s) from DEMs
## currently:
## - interpolation uncertainty
##
## TODO:
## source data uncertainty, slope (within cell)...
## generate a "Low Confidence Raster" <- uncertainty raster
## point density, land cover, classification (lidar)
## datum, resampling, spatial resolution, aspect, morph change.
## rank factors? 
## multi threading !!
##
### Code:

import os
import sys
import numpy as np
import math
import json

from osgeo import gdal

import cudem
from cudem import utils
from cudem import regions
from cudem import dlim
from cudem import demfun
from cudem import xyzfun
from cudem import waffles

## ==============================================
## Waffles Interpolation Uncertainty module
##
## Calculate the uncertainty of interpolated cells in a
## (waffles) raster
##
## ==============================================
class InterpolationUncertainty: #(waffles.Waffle):

    def __init__(self, dem=None, percentile=95, sims=None, chnk_lvl=None, max_sample=None):
        """calculate cell-level interpolation uncertainty

        Args:
          dem (Waffle): a waffles generated DEM (or constructed Waffle object)
          percentile (int): max percentile
          sims (int): number of split-sample simulations
          chnk_lvl (int): the 'chunk-level'
        """
        
        self.dem = dem
        self.percentile = percentile
        self.sims = sims
        self.max_sample = max_sample
        self.chnk_lvl = chnk_lvl        
        self._zones = ['LD0','LD1','LD2','MD0','MD1','MD2','HD0', 'HD1', 'HD2']
        #self._zones = ['low-dens','mid-dens','high-dens']
        self.prox = None
        self.slope = None

    def _mask_analysis(self, src_gdal, region=None):
        ds_config = demfun.gather_infos(src_gdal)
        if region is not None:
            srcwin = region.srcwin(ds_config['geoT'], ds_config['nx'], ds_config['ny'])
        else:
            srcwin = (0, 0, ds_config['nx'], ds_config['ny'])
            
        ds_arr = src_gdal.GetRasterBand(1).ReadAsArray(srcwin[0], srcwin[1], srcwin[2], srcwin[3])
        msk_sum = np.sum(ds_arr)
        msk_max = float(srcwin[2] * srcwin[3])
        msk_perc = float((msk_sum / msk_max) * 100.)
        dst_arr = None

        return(msk_sum, msk_max, msk_perc)

    def _prox_analysis(self, src_gdal, region=None):
        ds_config = demfun.gather_infos(src_gdal)
        if region is not None:
            srcwin = region.srcwin(ds_config['geoT'], ds_config['nx'], ds_config['ny'])
        else:
            srcwin = (0, 0, ds_config['nx'], ds_config['ny'])
            
        ds_arr = src_gdal.GetRasterBand(1).ReadAsArray(srcwin[0], srcwin[1], srcwin[2], srcwin[3])
        prox_perc = np.percentile(ds_arr, 95)
        dst_arr = None

        return(prox_perc)

    def _gen_prox(self):
        self.prox = '{}_prox.tif'.format(self.dem.mod)
        utils.echo_msg('generating proximity grid {}...'.format(self.prox))
        demfun.proximity(self.dem.mask_fn, self.prox)
        if self.dem.dst_srs is not None:
            demfun.set_srs(self.prox, self.dem.dst_srs)

    def _gen_slope(self):
        self.slope = '{}_slope.tif'.format(self.dem.mod)
        utils.echo_msg('generating slope grid {}...'.format(self.slope))
        demfun.slope(self.dem.fn, self.slope)
        if self.dem.dst_srs is not None:
            demfun.set_srs(self.slope, self.dem.dst_srs)

    def _regions_sort(self, trainers, t_num=25, verbose=False):
        """sort regions by distance; regions is a list of regions [xmin, xmax, ymin, ymax].

        returns the sorted region-list
        """

        train_sorted = []
        for z, train in enumerate(trainers):
            train_d = []
            np.random.shuffle(train)
            train_total = len(train)
            while True:
                if verbose:
                    utils.echo_msg_inline('sorting training tiles [{}]'.format(len(train)))
                    
                if len(train) == 0:
                    break
                
                this_center = train[0][0].center()
                train_d.append(train[0])
                train = train[1:]
                if len(train_d) > t_num or len(train) == 0:
                    break
                
                dsts = [utils.euc_dst(this_center, x[0].center()) for x in train]
                min_dst = np.percentile(dsts, 50)
                d_t = lambda t: utils.euc_dst(this_center, t[0].center()) > min_dst
                np.random.shuffle(train)
                train.sort(reverse=True, key=d_t)
            #if verbose:
            #    utils.echo_msg(' '.join([x[0].format('gmt') for x in train_d[:t_num]]))
                
            train_sorted.append(train_d)
            
        if verbose:
            utils.echo_msg_inline('sorting training tiles [OK]\n')
            
        return(train_sorted)

    def _gmt_select_split(self, o_xyz, sub_region, sub_bn, verbose = False):
        """split an xyz file into an inner and outer region.

        Args:
          o_xyz (str): a pathname to an xyz file
          sub_region (list): a region list [xmin, xmax, ymin, ymax]
          sub_bn (str): a basename for the selected data
          verbose (bool): increase verbosity

        Returns:
          list: [inner_region, outer_region]

        TODO: update for GDAL.
        """

        out_inner = None
        out_outer = None
        utils.run_cmd('gmt gmtset IO_COL_SEPARATOR space')
        gmt_s_inner = 'gmt gmtselect -V {} {} > {}_inner.xyz'.format(
            o_xyz, sub_region.format('gmt'), sub_bn
        )
        out, status = utils.run_cmd(gmt_s_inner, verbose = verbose)
        if status == 0:
            out_inner = '{}_inner.xyz'.format(sub_bn)
            
        gmt_s_outer = 'gmt gmtselect -V {} {} -Ir > {}_outer.xyz'.format(
            o_xyz, sub_region.format('gmt'), sub_bn
        )
        out, status = utils.run_cmd(gmt_s_outer, verbose=verbose)
        if status == 0:
            out_outer = '{}_outer.xyz'.format(sub_bn)
            
        return([out_inner, out_outer])

    def _err_fit_plot(self, xdata, ydata, out, fitfunc, bins_final, std_final, sampling_den, max_int_dist, dst_name='unc', xa='distance'):
        """plot a best fit plot with matplotlib

        Args:
          xdata (list): list of x-axis data
          ydata (list): list of y-axis data

        """

        #try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib.offsetbox import AnchoredText

        short_name="All Terrain"
        
        fig = plt.figure()
        ax = plt.subplot(111)

        plt_data=ax.scatter(bins_final, std_final, zorder=1, label='Error St. Dev.', marker="o", color="black", s=30)
        #plt_best_fit,=ax.plot(xdata,ydata, zorder=1, linewidth=2.0)
        plt_best_fit,=ax.plot(xdata, fitfunc(out, xdata), '-')
        
        box = ax.get_position()
        ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])

        plt.tick_params(
            axis='x',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            bottom='on',      # ticks along the bottom edge are off
            top='off',         # ticks along the top edge are off
            labelbottom='on') # labels along the bottom edge are off

        anchored_text = AnchoredText(short_name, loc=2)
        anchored_text.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
        ax.add_artist(anchored_text)

        anchored_text2 = AnchoredText(" $y = {%gx}^{%g}$ "%(out[1],out[2]), loc=1)
        #add r2 value using below
        #anchored_text2 = AnchoredText(" $y = {%gx}^{%g}$      $r^2=%g$ "%(coeff1,coeff2,rsquared), loc=1)
        anchored_text2.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
        ax.add_artist(anchored_text2)

        str_ss_samp_den="Sampling Density = " + str(sampling_den) + " %"
        anchored_text3 = AnchoredText(str_ss_samp_den, loc=4)
        anchored_text3.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
        ax.add_artist(anchored_text3)

        plt.legend([plt_data, plt_best_fit], ['Interpolation Error St Dev', 'Best-Fit Line'], loc='upper center', bbox_to_anchor=(0.5, 1.15), fancybox=True, shadow=True, ncol=2, fontsize=14)

        plt.tick_params(
            axis='x',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            bottom='on',      # ticks along the bottom edge are off
            top='off',         # ticks along the top edge are off
            labelbottom='on') # labels along the bottom edge are off
        
        plt.xlabel('Distance from Measurement (cells)', fontsize=14)
        plt.ylabel('Interpolation Error St Dev (m)', fontsize=14)
        plt.xlim(xmin=0)
        plt.xlim(xmax=int(max_int_dist)+1)
        plt.ylim(ymin=0)
        y_max=max(std_final)+(0.25*max(std_final))
        plt.ylim(ymax=y_max)

        # plt.plot(xdata, ydata, 'o')
        # plt.plot(xdata, fitfunc(out, xdata), '-')
        #plt.xlabel(xa)
        #plt.ylabel('Interpolation Error (m)')
        out_png = '{}_bf.png'.format(dst_name)
        plt.savefig(out_png)
        plt.close()

        #except: utils.echo_error_msg('you need to install matplotlib to run uncertainty plots...')

    def _err_scatter_plot(self, error_arr, dist_arr, mean, std, max_int_dist, bins_orig, sampling_den, dst_name='unc', xa='distance'):
        """plot a scatter plot with matplotlib

        Args:
          error_arr (array): an array of errors
          dist_arr (array): an array of distances

        """

        #try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib.offsetbox import AnchoredText

        short_name="All Terrain"

        fig = plt.figure()
        ax = plt.subplot(111)
        plt_data=ax.scatter(dist_arr, error_arr, zorder=1, label="Measurements", marker=".", color="black", s=20)
        plt_data_uncert=ax.errorbar(bins_orig, mean, yerr=std, fmt='r-', linewidth=3)
        box = ax.get_position()
        ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])

        plt.tick_params(
            axis='x',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            bottom='on',      # ticks along the bottom edge are off
            top='off',         # ticks along the top edge are off
            labelbottom='on') # labels along the bottom edge are off

        anchored_text = AnchoredText(short_name, loc=2)
        anchored_text.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
        ax.add_artist(anchored_text)

        str_ss_samp_den="Sampling Density = " + str(sampling_den) + " %"
        anchored_text3 = AnchoredText(str_ss_samp_den, loc=4)
        anchored_text3.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
        ax.add_artist(anchored_text3)

        plt.legend([plt_data, plt_data_uncert], ["Interpolation Error", "Mean +/- St. Deviation"], loc="upper center", bbox_to_anchor=(0.5, 1.15), fancybox=True, shadow=True, ncol=2, fontsize=14)
        plt.xlabel("Distance from Measurement (cells)", fontsize=14)
        plt.ylabel("Interpolation Error (m)", fontsize=14)
        plt.xlim(xmin=0)
        plt.xlim(xmax=int(max_int_dist)+1)

        #plt.xlabel(xa)
        #plt.ylabel('Interpolation Error (m)')
        out_png = "{}_scatter.png".format(dst_name)
        plt.savefig(out_png)
        plt.close()

        #xcept: utils.echo_error_msg('you need to install matplotlib to run uncertainty plots...')

    def _err2coeff(self, err_arr, sampling_den, coeff_guess=[0, 0.1, 0.2], dst_name='unc', xa='distance'):
        """calculate and plot the error coefficient given err_arr which is 
        a 2 col array with `err dist

        Args:
          error_arr (array): an array of errors and distances

        Returns:
          list: [coefficient-list]
        """

        from scipy import optimize
        
        error = err_arr[:,0]
        distance = err_arr[:,1]
        max_err = np.max(error)
        min_err = np.min(error)
        max_int_dist = int(np.max(distance))
        nbins = 10
        n, _ = np.histogram(distance, bins = nbins)
        while 0 in n:
            nbins -= 1
            n, _ = np.histogram(distance, bins=nbins)
            
        serror, _ = np.histogram(distance, bins=nbins, weights=error)
        serror2, _ = np.histogram(distance, bins=nbins, weights=error**2)

        mean = serror / n
        std = np.sqrt(serror2 / n - mean * mean)
        ydata = np.insert(std, 0, 0)
        bins_orig=(_[1:] + _[:-1]) / 2
        
        xdata = np.insert(bins_orig, 0, 0)
        xdata[xdata - 0 < 0.0001] = 0.0001
        while len(xdata) < 3:
            xdata = np.append(xdata, 0)
            ydata = np.append(ydata, 0)
            
        fitfunc = lambda p, x: p[0] + p[1] * (x ** p[2])
        errfunc = lambda p, x, y: y - fitfunc(p, x)
        out, cov, infodict, mesg, ier = optimize.leastsq(
            errfunc, coeff_guess, args=(xdata, ydata), full_output=True
        )
        #ry:

        self._err_fit_plot(xdata, ydata, out, fitfunc, bins_orig, std, sampling_den, max_int_dist, dst_name, xa)
        self._err_scatter_plot(error, distance, mean, std, max_int_dist, bins_orig, sampling_den, dst_name, xa)
        #except:
        #   utils.echo_error_msg('unable to generate error plots, please check configs.')
            
        return(out)

    def _sub_region_analysis(self, sub_regions):
        """sub-region analysis"""
        
        sub_zones = {}
        dem_ds = gdal.Open(self.dem.fn)
        msk_ds = gdal.Open(self.dem.mask_fn)
        prox_ds = gdal.Open(self.prox)
        slp_ds = gdal.Open(self.slope)
        with utils.CliProgress(total=len(sub_regions), message='analyzing {} sub-regions'.format(len(sub_regions))) as pbar:
            for sc, sub_region in enumerate(sub_regions):
                pbar.update()
                s_sum, s_g_max, s_perc = self._mask_analysis(msk_ds, region=sub_region)
                #s_sum, s_g_max, s_perc = self._mask_analysis(slp_ds, region=sub_region)
                s_dc = demfun.gather_infos(dem_ds, region=sub_region, scan=True)
                p_perc = self._prox_analysis(prox_ds, region=sub_region)
                slp_perc = self._prox_analysis(slp_ds, region=sub_region)
                #slp_perc = 0

                # if p_perc < self.prox_perc_33 or abs(p_perc - self.prox_perc_33) < 0.01:
                #     zone = self._zones[2]
                # elif p_perc < self.prox_perc_66 or abs(p_perc - self.prox_perc_66) < 0.01:
                #     zone = self._zones[1]
                # else:
                #     zone = self._zones[0]
                zone = None
                if p_perc < self.prox_perc_33 or abs(p_perc - self.prox_perc_33) < 0.01:
                    if slp_perc < self.slp_perc_33 or abs(slp_perc - self.slp_perc_33) < 0.01:
                        zone = self._zones[6]
                    elif slp_perc < self.slp_perc_66 or abs(slp_perc - self.slp_perc_66) < 0.01:
                        zone = self._zones[7]
                    else:
                        zone = None
                        #    zone = self._zones[8]

                elif p_perc < self.prox_perc_66 or abs(p_perc - self.prox_perc_66) < 0.01:
                    if slp_perc < self.slp_perc_33 or abs(slp_perc - self.slp_perc_33) < 0.01:
                        zone = self._zones[3]
                    elif slp_perc < self.slp_perc_66 or abs(slp_perc - self.slp_perc_66) < 0.01:
                        zone = self._zones[4]
                    else:
                        zone = None
                        #    zone = self._zones[5]
                else:
                    if slp_perc < self.slp_perc_33 or abs(slp_perc - self.slp_perc_33) < 0.01:
                        zone = self._zones[0]
                    elif slp_perc < self.slp_perc_66 or abs(slp_perc - self.slp_perc_66) < 0.01:
                        zone = self._zones[1]
                    else:
                        zone = None
                        #    zone = self._zones[2]

                if zone is not None:
                    sub_zones[sc + 1] = [sub_region, s_g_max, s_sum, s_perc, p_perc, zone]
                #sub_zones[sc + 1] = [sub_region, s_g_max, s_sum, s_perc, slp_perc, zone]
                #sub_zones[sc + 1] = [sub_region, s_g_max, s_sum, s_perc, p_perc, slp_perc, s_dc['zr'][0], s_dc['zr'][1], zone]
            
        dem_ds = msk_ds = prox_ds = slp_ds = None
        return(sub_zones)

    def _split_sample(self, trainers, perc, max_dist):
        """split-sample simulations and error calculations
        sims = max-simulations
        """
            
        #utils.echo_msg('simulation\terrors\tproximity-coeff\tp_diff\tslp-coeff\tslp_diff')
        #utils.echo_msg('simulation\terrors\tmin-error\tmax-error\tmean-error\tproximity-coeff\tp_diff')
        utils.echo_msg('simulation\terrors\tmean-error\tproximity-coeff\tp_diff')
        sim = 0
        status = 0
        last_ec_d = None
        #s_dp = sub_dp = None
        s_dp = []

        with utils.CliProgress(message='performing MAX {} SPLIT-SAMPLE simulations looking for MIN {} sample errors'.format(self.sims, self.max_sample)) as pbar:
            while True:
                status = 0
                sim += 1
                #utils.echo_msg('sorting training tiles by distance...')
                trains = self._regions_sort(trainers, verbose=True)
                tot_trains = len([x for s in trains for x in s])
                #utils.echo_msg('sorted sub-regions into {} training tiles.'.format(tot_trains))

                for z, train in enumerate(trains):
                    train_h = train[:25]
                    ss_samp = perc

                    ## ==============================================
                    ## perform split-sample analysis on each training region.
                    ## ==============================================
                    for n, sub_region in enumerate(train_h):
                        #print(sub_region)
                        #_prog.update()
                        pbar.update()
                        ss_samp = perc
                        this_region = sub_region[0].copy()
                        #if sub_region[3] < ss_samp:
                        #   ss_samp = sub_region[3]
                        ss_samp = sub_region[3]

                        ## ==============================================
                        ## extract the xyz data for the region from the DEM
                        ## ==============================================
                        o_xyz = '{}_{}.xyz'.format(self.dem.name, n)
                        ds = gdal.Open(self.dem.fn)
                        ds_config = demfun.gather_infos(ds)
                        b_region = this_region.copy()
                        b_region.buffer(pct=20, x_inc=self.dem.xinc, y_inc=self.dem.yinc)
                        srcwin = b_region.srcwin(ds_config['geoT'], ds_config['nx'], ds_config['ny'])
                        #cnt = 0

                        ## TODO: extract weights here as well...
                        with open(o_xyz, 'w') as o_fh:
                            for xyz in demfun.parse(ds, srcwin=srcwin, mask=self.dem.mask_fn):
                                #cnt+=1
                                xyz.dump(dst_port=o_fh)

                        ds = None

                        if os.stat(o_xyz).st_size != 0:
                            ## ==============================================
                            ## split the xyz data to inner/outer; outer is
                            ## the data buffer, inner will be randomly sampled
                            ## ==============================================

                            s_inner, s_outer = self._gmt_select_split(
                                o_xyz, this_region, 'sub_{}'.format(n), verbose=False
                            )
                            if os.stat(s_inner).st_size != 0:
                                sub_xyz = np.loadtxt(s_inner, ndmin=2, delimiter=' ')                        
                                ss_len = len(sub_xyz)
                                #sx_cnt = int(sub_region[2] * (ss_samp / 100.)) if ss_samp is not None else ss_len-1
                                #if ss_samp is not None:
                                #sx_cnt = int(sub_region[1] * (ss_samp / 100.)) + 1
                                sx_cnt = int(ss_len * (ss_samp / 100.))
                                #else:
                                #    #sx_cnt = ss_len-1
                                #    sx_cnt = 1

                                #if sx_cnt >= ss_len: sx_cnt = ss_len=-1
                                sx_cnt = 1 if sx_cnt < 1 or sx_cnt >= ss_len else sx_cnt

                                #print('sx_cnt', sx_cnt)
                                #print('ss_len', ss_len)
                                #print('sub_xyz_head', len(sub_xyz[sx_cnt:]))
                                #_prog.update('performing MAX {} SPLIT-SAMPLE simulations looking for MIN {} sample errors: {} {} {} {} {}'.format(self.sims, self.max_sample, sub_region[5], ss_samp, sx_cnt, ss_len, len(sub_xyz[sx_cnt:])))
                                sub_xyz_head = 'sub_{}_head_{}.xyz'.format(n, sx_cnt)
                                np.random.shuffle(sub_xyz)
                                np.savetxt(sub_xyz_head, sub_xyz[:sx_cnt], '%f', ' ')

                                ## ==============================================
                                ## generate the random-sample DEM
                                ## ==============================================
                                waff = waffles.WaffleFactory(
                                    mod=self.dem.mod,
                                    data=[s_outer, sub_xyz_head],
                                    src_region=b_region,
                                    xinc=self.dem.xinc,
                                    yinc=self.dem.yinc,
                                    name='sub_{}'.format(n),
                                    node=self.dem.node,
                                    fmt=self.dem.fmt,
                                    extend=self.dem.extend,
                                    extend_proc=self.dem.extend_proc,
                                    weights=self.dem.weights,
                                    sample=self.dem.sample,
                                    clip=self.dem.clip,
                                    dst_srs=self.dem.dst_srs,
                                    mask=True,
                                    verbose=True,
                                    clobber=True,
                                )
                                waff.mod_args = self.dem.mod_args
                                wf = waff.acquire().generate()

                                if wf.valid_p():
                                    ## ==============================================
                                    ## generate the random-sample data PROX and SLOPE
                                    ## ==============================================
                                    sub_prox = '{}_prox.tif'.format(wf.name)
                                    demfun.proximity('{}_m.tif'.format(wf.name), sub_prox)
                                    #sub_slp = '{}_slp.tif'.format(wf.name)
                                    #demfun.slope(wf.fn, sub_slp)

                                    ## ==============================================
                                    ## Calculate the random-sample errors
                                    ## ==============================================
                                    sub_xyd = demfun.query(sub_xyz[sx_cnt:], wf.fn, 'xyd')
                                    sub_dp = demfun.query(sub_xyd, sub_prox, 'xyzg')
                                    #sub_ds = demfun.query(sub_dp, self.slope, 'g')
                                    #if len(sub_dp) > 0:
                                    #   if sub_dp.shape[0] == sub_ds.shape[0]:
                                    #       sub_dp = np.append(sub_dp, sub_ds, 1)
                                    #   else: sub_dp = []
                                else:
                                    sub_dp = None

                                utils.remove_glob(sub_xyz_head)
                                if sub_dp is not None and len(sub_dp) > 0:
                                    try:
                                        s_dp = np.vstack((s_dp, sub_dp))
                                    except:
                                        s_dp = sub_dp

                        utils.remove_glob(o_xyz, 'sub_{}*'.format(n))

                if s_dp is not None and len(s_dp) > 0:
                    d_max = self.region_info[self.dem.name][4]
                    #s_max = self.region_info[self.dem.name][5]
                    s_dp = s_dp[s_dp[:,3] < max_dist,:]
                    s_dp = s_dp[s_dp[:,3] >= 1,:]
                    prox_err = s_dp[:,[2,3]]

                    if last_ec_d is None:
                        last_ec_d = [0, 0.1, 0.2]
                        last_ec_diff = 10
                    else:
                        last_ec_diff = abs(last_ec_d[2] - last_ec_d[1])

                    ec_d = self._err2coeff(prox_err[:50000000], perc, coeff_guess=last_ec_d, dst_name=self.dem.name + '_prox', xa='Distance to Nearest Measurement (cells)')
                    #ec_d = self._err2coeff(prox_err[:50000000], coeff_guess=[0, .1, .2], dst_name=self.dem.name + '_prox', xa='distance')
                    ec_diff = abs(ec_d[2] - ec_d[1])
                    ec_l_diff = abs(last_ec_diff - ec_diff)
                    #s_dp = s_dp[s_dp[:,4] < s_max,:]
                    #slp_err = s_dp[:,[2,4]]
                    #ec_s = self._err2coeff(slp_err[:50000000], coeff_guess=[0, 0.1, 0.2], dst_name = self.dem.name + '_slp', xa = 'slope')
                    # ec_s = [0, 1, 2]
                    #utils.echo_msg('{}\t{}\t{}\t{}\t{}\t{}'.format(sim, len(s_dp), ec_d, ec_d[2] - ec_d[1], ec_s, ec_s[2] - ec_s[1]))
                    utils.echo_msg('{}\t{}\t{}\t{}\t{}'.format(sim, len(s_dp), np.mean(prox_err, axis=0)[0], ec_d, ec_l_diff))

                    if len(s_dp) < self.max_sample:
                        last_ec_d = ec_d
                        continue

                    if ec_d[0] == 0 and ec_d[1] == 0.1 and ec_d[2] == 0.2:
                        last_ec_d = ec_d
                        continue

                    if sim >= int(self.sims):
                        break

                    if abs(last_ec_diff - ec_diff) < 0.01:
                        break

                    # if len(s_dp) >= int(self.region_info[self.dem.name][1] / 10):
                    #     break

                    last_ec_d = ec_d
                else:
                    utils.echo_msg('{}\t{}\t{}\t{}\t{}'.format(sim, len(s_dp), None, None, None, None))

        np.savetxt('prox_err.xyz', prox_err, '%f', ' ')                
        #_prog.end(status, 'performed {} SPLIT-SAMPLE simulations'.format(sim))
        return([ec_d])
        
    def run(self):
        s_dp = s_ds = None
        unc_out = {}
        zones = ['low-dens-low-slp','low-dens-mid-slp','low-dens-high-slp','mid-dens-low-slp','mid-dens-mid-slp','mid-dens-high-slp','high-dens-low-slp', 'high-dens-mid-slp', 'high-dens-low-slp']
        #zones = ['low-dens','mid-dens','high-dens','low-slp','mid-slp','high-slp']
        #zones = ['low-dens','mid-dens','high-dens']
        utils.echo_msg('running INTERPOLATION uncertainty module using {}...'.format(self.dem.mod))
        
        if self.prox is None:
            self._gen_prox()

        if self.slope is None:
            self._gen_slope()
            
        ## ==============================================
        ## region and der. analysis
        ## ==============================================
        self.region_info = {}
        msk_ds = gdal.Open(self.dem.mask_fn)
        num_sum, g_max, num_perc = self._mask_analysis(msk_ds)
        msk_ds = None

        self.prox_percentile = demfun.percentile(self.prox, self.percentile)
        self.prox_perc_33 = demfun.percentile(self.prox, 25)
        self.prox_perc_66 = demfun.percentile(self.prox, 75)
        self.prox_perc_100 = demfun.percentile(self.prox, 100)

        self.slp_percentile = demfun.percentile(self.slope, self.percentile)
        self.slp_perc_33 = demfun.percentile(self.slope, 25)
        self.slp_perc_66 = demfun.percentile(self.slope, 75)
        self.slp_perc_100 = demfun.percentile(self.slope, 100)

        #self.region_info[self.dem.name] = [self.dem.region, g_max, num_sum, num_perc, self.prox_percentile, self.slp_percentile]
        self.region_info[self.dem.name] = [self.dem.region, g_max, num_sum, num_perc, self.prox_percentile] 
        for x in self.region_info.keys():
            utils.echo_msg('region: {}: {}'.format(x, self.region_info[x]))

        ## ==============================================
        ## chunk region into sub regions
        ## ==============================================
        #print(self.region_info[self.dem.name])
        #print(num_sum, g_max, num_perc)
        chnk_inc = int((num_sum / math.sqrt(g_max)) / num_perc) * 2
        #chnk_inc = 24
        sub_regions = self.dem.region.chunk(self.dem.xinc, chnk_inc)
        utils.echo_msg('chunked region into {} sub-regions @ {}x{} cells.'.format(len(sub_regions), chnk_inc, chnk_inc))

        ## ==============================================
        ## sub-region analysis
        ## ==============================================
        sub_zones = self._sub_region_analysis(sub_regions)
        #sub_zones.append(self._sub_region_analysis(sub_regions, prox_or_slp='slp'))

        ## ==============================================
        ## sub-region density and percentiles
        ## ==============================================
        s_dens = np.array([sub_zones[x][3] for x in sub_zones.keys()])
        s_5perc = np.percentile(s_dens, 5)
        s_dens = None
        #utils.echo_msg('Sampling density for region is: {:.16f}'.format(s_5perc))
        utils.echo_msg('Sampling density for region is: {:.16f}'.format(num_perc))

        ## ==============================================
        ## zone analysis / generate training regions
        ## ==============================================
        trainers = []
        t_perc = 95
        s_perc = 50

        for z, this_zone in enumerate(self._zones):
            #print(sub_zones) #sub_zones[x][8] (with slope) [5] for prox only
            #print(this_zone)
            tile_set = [sub_zones[x] for x in sub_zones.keys() if sub_zones[x][5] == self._zones[z]]
            #print(tile_set)
            if len(tile_set) > 0:
                d_50perc = np.percentile(np.array([x[3] for x in tile_set]), 50)
            else:
                continue
            
            t_trainers = [x for x in tile_set if x[3] < d_50perc or abs(x[3] - d_50perc) < 0.01]
            utils.echo_msg(
                'possible {} training zones: {} @ MAX {}'.format(
                    self._zones[z].upper(), len(t_trainers), d_50perc
                )
            )
            trainers.append(t_trainers)
        #utils.echo_msg('sorting training tiles by distance...')
        # trains = self._regions_sort(trainers, verbose = False)
        # tot_trains = len([x for s in trains for x in s])
        # utils.echo_msg('sorted sub-regions into {} training tiles.'.format(tot_trains))
        utils.echo_msg('analyzed {} sub-regions.'.format(len(sub_regions)))

        ## ==============================================
        ## split-sample simulations and error calculations
        ## sims = max-simulations
        ## ==============================================
        if self.sims is None:
            #self.sims = int(len(sub_regions)/tot_trains)
            self.sims = 12

        if self.max_sample is None:
            self.max_sample = int((self.region_info[self.dem.name][1] - self.region_info[self.dem.name][2]) * .25)
            #self.max_sample = 1
            
        #ec_d, ec_s = self._split_sample(trains, s_5perc)
        #ec_d = self._split_sample(trains, s_5perc)[0]
        ec_d = self._split_sample(trainers, num_perc, chnk_inc/2)[0]

        ## ==============================================
        ## Save/Output results
        ## apply error coefficient to full proximity grid
        ## TODO: USE numpy/gdal instead!
        ## ==============================================
        utils.echo_msg('applying coefficient to PROXIMITY grid')
        #if self.dem.gc['GMT'] is None:
        utils.run_cmd('gdal_calc.py -A {} --outfile {}_prox_unc.tif --calc "{}+({}*(A**{}))"'.format(self.prox, self.dem.name, 0, ec_d[1], ec_d[2]), verbose = True)
        #else:
        #    math_cmd = 'gmt grdmath {} 0 AND ABS {} POW {} MUL {} ADD = {}_prox_unc.tif=gd+n-9999:GTiff\
        #'.format(self.prox, ec_d[2], ec_d[1], 0, self.dem.name)
        #utils.run_cmd(math_cmd, verbose = self.dem.verbose)
            
        if self.dem.dst_srs is not None:
            status = demfun.set_srs('{}_prox_unc.tif'.format(self.dem.name), src_srs=self.dem.dst_srs)
            
        utils.echo_msg('applied coefficient {} to PROXIMITY grid'.format(ec_d))

        # utils.echo_msg('applying coefficient to SLOPE grid')
        # #if self.dem.gc['GMT'] is None:
        # utils.run_cmd('gdal_calc.py -A {} --outfile {}_slp_unc.tif --calc "{}+({}*(A**{}))"'.format(self.slope, self.dem.name, 0, ec_s[1], ec_s[2]), verbose = True)
        # #else:
        # #    math_cmd = 'gmt grdmath {} 0 AND ABS {} POW {} MUL {} ADD = {}_slp_unc.tif=gd+n-9999:GTiff\
        # #    '.format(self.slope, ec_s[2], ec_s[1], 0, self.dem.name)
        # #    utils.run_cmd(math_cmd, verbose = self.dem.verbose)
        # if self.dem.dst_srs is not None: status = demfun.set_srs('{}_prox_unc.tif'.format(self.dem.name), src_srs=self.dem.dst_srs)
        # utils.echo_msg('applied coefficient {} to SLOPE grid'.format(ec_s))

        utils.remove_glob(self.prox)
        utils.remove_glob(self.slope)

        unc_out['prox_unc'] = ['{}_prox_unc.tif'.format(self.dem.name), 'raster']
        unc_out['prox_bf'] = ['{}_prox_bf.png'.format(self.dem.name), 'image']
        unc_out['prox_scatter'] = ['{}_prox_scatter.png'.format(self.dem.name), 'image']

        return(unc_out, 0)

## ==============================================
## Command-line Interface (CLI)
## $ uncertainties
##
## uncertainties cli
## ==============================================
_waffles_module_short_desc = lambda: ', '.join(
    ['{}'.format(key) for key in waffles.WaffleFactory()._modules])
_waffles_module_long_desc = lambda x: 'waffles modules:\n% waffles ... -M <mod>:key=val:key=val...\n\n  ' + '\n  '.join(
    ['\033[1m{:14}\033[0m{}\n'.format(key, x[key]['description']) for key in x]) + '\n'

uncertainties_cli_usage = """{cmd} ({wf_version}): Generate DEMs and derivatives.

usage: {cmd} [OPTIONS] WAFFLES_CONFIG

Options:
  -q, --quiet\t\tLower verbosity to a quiet.

  --help\t\tPrint the usage text
  --version\t\tPrint the version information

Generate a waffles config files with the `waffles` command.

CIRES DEM home page: <http://ciresgroups.colorado.edu/coastalDEM>
""".format(cmd=os.path.basename(sys.argv[0]),
           dl_formats=utils._cudem_module_name_short_desc(dlim.DatasetFactory.data_types),
           modules=_waffles_module_short_desc(),
           wf_version=cudem.__version__)

def uncertainties_cli(argv = sys.argv):
    """run uncertainties from command-line

    See `uncertainties_cli_usage` for full cli options.
    """
    
    wg_user = None
    status = 0
    i = 1
    wg = {}
    wg['verbose'] = True
    wg['clobber'] = False
    sims = None
    max_sample = None
    
    while i < len(argv):
        arg = argv[i]
        if arg == '--quiet' or arg == '-q': wg['verbose'] = False
        elif arg == '-sims' or arg == '--sims' or arg == '-S':
            sims = utils.int_or(argv[i + 1])
            i = i + 1
        elif arg == '-sample' or arg == '--sample' or arg == '-P':
            max_sample = utils.int_or(argv[i + 1])
            i = i + 1
        elif arg == '--help' or arg == '-h':
            sys.stderr.write(uncertainties_cli_usage)
            sys.exit(0)
        elif arg == '--version' or arg == '-v':
            sys.stdout.write('{}\n'.format(cudem.__version__))
            sys.exit(0)
        elif arg[0] == '-':
            sys.stdout.write(uncertainties_cli_usage)
            utils.echo_error_msg('{} is not a valid waffles cli switch'.format(arg))
            sys.exit(0)
        else: wg_user = arg
        i += 1

    ## ==============================================
    ## load the user wg json and run waffles with that.
    ## ==============================================
    if wg_user is not None:
        if os.path.exists(wg_user):
            try:
                with open(wg_user, 'r') as wgj:
                    wg = json.load(wgj)
                    if wg['src_region'] is not None:
                        wg['src_region'] = regions.Region().from_list(
                            wg['src_region']
                        )

                    this_waffle = waffles.WaffleFactory(**wg).acquire()
                    this_waffle.mask = True
                    this_waffle.clobber = False

                    if not this_waffle.valid_p():
                        this_waffle.generate()

                    ## run interpolation uncertainty
                    i = InterpolationUncertainty(dem=this_waffle, sims=sims, max_sample=max_sample).run()
                    utils.echo_msg(this_waffle)

                    sys.exit(0)
            except Exception as e:
                #utils.echo_error_msg(e)
                sys.exit(-1)
        else:
            utils.echo_error_msg(
                'specified waffles config file does not exist, {}'.format(wg_user)
            )
            sys.stderr.write(waffles_cli_usage)
            sys.exit(-1)
    else:
        utils.echo_error_msg(
            'you must supply a waffles config file; see waffles --help for more information.'
        )
        sys.exit(-1)

### End        
