### uncertainties.py
##
## Copyright (c) 2010 - 2021 CIRES Coastal DEM Team
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
### Code:

import os
import sys
from osgeo import gdal
import numpy as np
import math
import cudem
from cudem import utils
from cudem import regions
from cudem import dlim
from cudem import demfun
from cudem import xyzfun
from cudem import waffles

## ==============================================
## Waffles Interpolation Uncertainty module
## ==============================================
class InterpolationUncertainty: #(waffles.Waffle):

    def __init__(self, dem=None, percentile=95, sims=None, chnk_lvl=None):
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
        self.chnk_lvl = chnk_lvl        
        #self._zones = ['low-dens','mid-dens','high-dens','low-slp','mid-slp','high-slp']
        self._zones = ['low-dens','mid-dens','high-dens']

        self.prox = None
        self.slope = None

    def _mask_analysis(self, src_gdal, region=None):
        ds_config = demfun.gather_infos(src_gdal)
        if region is not None:
            srcwin = region.srcwin(ds_config['geoT'], ds_config['nx'], ds_config['ny'])
        else: srcwin = (0, 0, ds_config['nx'], ds_config['ny'])
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
        else: srcwin = (0, 0, ds_config['nx'], ds_config['ny'])
        ds_arr = src_gdal.GetRasterBand(1).ReadAsArray(srcwin[0], srcwin[1], srcwin[2], srcwin[3])

        prox_perc = np.percentile(ds_arr, 95)
        dst_arr = None

        return(prox_perc)

    def _gen_prox(self):
        self.prox = '{}_prox.tif'.format(self.dem.mod)
        utils.echo_msg('generating proximity grid {}...'.format(self.prox))
        demfun.proximity(self.dem.mask_fn, self.prox)
        if self.dem.epsg is not None: demfun.set_epsg(self.prox, self.dem.epsg)

    def _gen_slope(self):
        self.slope = '{}_slope.tif'.format(self.dem.mod)
        utils.echo_msg('generating proximity grid {}...'.format(self.slope))
        demfun.proximity(self.dem.fn, self.slope)
        if self.dem.epsg is not None: demfun.set_epsg(self.slope, self.dem.epsg)

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
            if verbose:
                utils.echo_msg(' '.join([x[0].format('gmt') for x in train_d[:t_num]]))
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
        gmt_s_inner = 'gmt gmtselect -V {} {} > {}_inner.xyz'.format(o_xyz, sub_region.format('gmt'), sub_bn)
        out, status = utils.run_cmd(gmt_s_inner, verbose = verbose)
        if status == 0:
            out_inner = '{}_inner.xyz'.format(sub_bn)
        gmt_s_outer = 'gmt gmtselect -V {} {} -Ir > {}_outer.xyz'.format(o_xyz, sub_region.format('gmt'), sub_bn)
        out, status = utils.run_cmd(gmt_s_outer, verbose=verbose)
        if status == 0:
            out_outer = '{}_outer.xyz'.format(sub_bn)
        return([out_inner, out_outer])

    def _err_fit_plot(self, xdata, ydata, out, fitfunc, dst_name='unc', xa='distance'):
        """plot a best fit plot with matplotlib

        Args:
          xdata (list): list of x-axis data
          ydata (list): list of y-axis data

        """

        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            from matplotlib.offsetbox import AnchoredText

            plt.plot(xdata, ydata, 'o')
            plt.plot(xdata, fitfunc(out, xdata), '-')
            plt.xlabel(xa)
            plt.ylabel('error (m)')
            out_png = '{}_bf.png'.format(dst_name)
            plt.savefig(out_png)
            plt.close()

        except: utils.echo_error_msg('you need to install matplotlib to run uncertainty plots...')

    def _err_scatter_plot(self, error_arr, dist_arr, dst_name='unc', xa='distance'):
        """plot a scatter plot with matplotlib

        Args:
          error_arr (array): an array of errors
          dist_arr (array): an array of distances

        """

        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            from matplotlib.offsetbox import AnchoredText

            plt.scatter(dist_arr, error_arr)
            #plt.title('Scatter')
            plt.xlabel(xa)
            plt.ylabel('error (m)')
            out_png = '{}_scatter.png'.format(dst_name)
            plt.savefig(out_png)
            plt.close()

        except: utils.echo_error_msg('you need to install matplotlib to run uncertainty plots...')

    def _err2coeff(self, err_arr, coeff_guess=[0, 0.1, 0.2], dst_name='unc', xa='distance'):
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

        max_int_dist = np.max(distance)
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
        fitfunc = lambda p, x: p[0] + p[1] * (x ** p[2])
        errfunc = lambda p, x, y: y - fitfunc(p, x)
        out, cov, infodict, mesg, ier = optimize.leastsq(errfunc, coeff_guess, args=(xdata, ydata), full_output=True)
        try:
            self._err_fit_plot(xdata, ydata, out, fitfunc, dst_name, xa)
            self._err_scatter_plot(error, distance, dst_name, xa)
        except: utils.echo_error_msg('unable to generate error plots, please check configs.')
        return(out)

    def _sub_region_analysis(self, sub_regions):
        """sub-region analysis"""
        
        utils.echo_msg('analyzing {} sub-regions...'.format(len(sub_regions)))
        sub_zones = {}
        dem_ds = gdal.Open(self.dem.fn)
        msk_ds = gdal.Open(self.dem.mask_fn)
        prox_ds = gdal.Open(self.prox)
        #slp_ds = gdal.Open(self.slope)
        _prog = utils.CliProgress('analyzing {} sub-regions.'.format(len(sub_regions)))
        for sc, sub_region in enumerate(sub_regions):
            _prog.update_perc((sc, len(sub_regions)))
            #utils.echo_msg_inline('analyzing sub-regions [{}]'.format(sc))
            s_sum, s_g_max, s_perc = self._mask_analysis(msk_ds, region=sub_region)
            p_perc = self._prox_analysis(prox_ds, region=sub_region)
            #slp_perc = self._prox_analysis(slp_ds, region=sub_region)
            #slp_perc = 0
            s_dc = demfun.gather_infos(dem_ds, region=sub_region, scan=True)
            if p_perc < self.prox_perc_33 or abs(p_perc - self.prox_perc_33) < 0.01:
                zone = self._zones[2]
            elif p_perc < self.prox_perc_66 or abs(p_perc - self.prox_perc_66) < 0.01:
                zone = self._zones[1]
            else:
                zone = self._zones[0]
            # if slp_perc < self.slp_perc_33 or abs(slp_perc - self.slp_perc_33) < 0.01:
            #     zone = self._zones[3]
            # elif slp_perc < self.slp_perc_66 or abs(slp_perc - self.slp_perc_66) < 0.01:
            #     zone = self._zones[4]
            # else:
            #     zone = self._zones[5]

            #sub_zones[sc + 1] = [sub_region, s_g_max, s_sum, s_perc, p_perc, slp_perc, s_dc['zr'][0], s_dc['zr'][1], zone]
            sub_zones[sc + 1] = [sub_region, s_g_max, s_sum, s_perc, p_perc, zone]
        dem_ds = msk_ds = prox_ds = slp_ds = None
        _prog.end(0, 'analyzed {} sub-regions.'.format(len(sub_regions)))
        #utils.echo_msg_inline('analyzing sub-regions [OK]\n')
        return(sub_zones)

    def _split_sample(self, trains, perc):
        """split-sample simulations and error calculations
        sims = max-simulations
        """
            
        #utils.echo_msg('performing MAX {} SPLIT-SAMPLE simulations...'.format(self.sims))
        _prog = utils.CliProgress('performing MAX {} SPLIT-SAMPLE simulations'.format(self.sims))
        #utils.echo_msg('simulation\terrors\tproximity-coeff\tp_diff\tslp-coeff\tslp_diff')
        utils.echo_msg('simulation\terrors\tproximity-coeff\tp_diff')
        
        sim = 0
        status = 0
        last_ec_d = None
        
        while True:
            status = 0
            sim += 1
            #trains = self._regions_sort(trainers, verbose=False)
            for z, train in enumerate(trains):
                train_h = train[:25]
                ss_samp = perc

                ## ==============================================
                ## perform split-sample analysis on each training region.
                ## ==============================================
                for n, sub_region in enumerate(train_h):
                    ss_samp = perc
                    #perc = int(float(n+(len(train_h) * z))/(len(train_h)*len(trains)) * 100)
                    #_prog.update_perc((int(float(n+(len(train_h) * z))), len(train_h)*len(trains)))
                    _prog.update()
                    this_region = sub_region[0]
                    if sub_region[3] < ss_samp:
                        ss_samp = None

                    ## ==============================================
                    ## extract the xyz data for the region from the DEM
                    ## ==============================================
                    o_xyz = '{}_{}.xyz'.format(self.dem.name, n)
                    ds = gdal.Open(self.dem.fn)
                    ds_config = demfun.gather_infos(ds)
                    b_region = this_region
                    b_region.buffer(20*self.dem.inc)
                    srcwin = b_region.srcwin(ds_config['geoT'], ds_config['nx'], ds_config['ny'])
                    with open(o_xyz, 'w') as o_fh:
                        for xyz in demfun.parse(ds, srcwin=srcwin, mask=self.dem.mask_fn):
                            xyz.dump(dst_port=o_fh)
                    ds = None

                    if os.stat(o_xyz).st_size != 0:
                        ## ==============================================
                        ## split the xyz data to inner/outer; outer is
                        ## the data buffer, inner will be randomly sampled
                        ## ==============================================
                        s_inner, s_outer = self._gmt_select_split(o_xyz, this_region, 'sub_{}'.format(n), verbose=False)
                        if os.stat(s_inner).st_size != 0:
                            sub_xyz = np.loadtxt(s_inner, ndmin=2, delimiter=' ')
                        else: sub_xyz = []
                        ss_len = len(sub_xyz)
                        if ss_samp is not None:
                            sx_cnt = int(sub_region[1] * (ss_samp / 100.)) + 1
                        else: sx_cnt = 1

                        sub_xyz_head = 'sub_{}_head.xyz'.format(n)
                        np.random.shuffle(sub_xyz)
                        np.savetxt(sub_xyz_head, sub_xyz[:sx_cnt], '%f', ' ')

                        ## ==============================================
                        ## generate the random-sample DEM
                        ## ==============================================
                        waff = waffles.WaffleFactory(mod=self.dem.mod, data=[s_outer, sub_xyz_head], src_region=this_region,
                                                   inc=self.dem.inc, name='sub_{}'.format(n), node=self.dem.node, fmt=self.dem.fmt,
                                                   extend=self.dem.extend, extend_proc=self.dem.extend_proc, weights=self.dem.weights,
                                                   sample=self.dem.sample, clip=self.dem.clip, chunk=None, epsg=self.dem.epsg,
                                                   mask=True, verbose=True, clobber=True)
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
                            #sub_dp = gdalfun.gdal_query(sub_xyd, sub_prox, 'zg')
                            sub_dp = demfun.query(sub_xyd, sub_prox, 'xyzg')
                            #sub_ds = demfun.query(sub_dp, self.slope, 'g')
                            
                            #if len(sub_dp) > 0:
                            #   if sub_dp.shape[0] == sub_ds.shape[0]:
                            #       sub_dp = np.append(sub_dp, sub_ds, 1)
                            #   else: sub_dp = []
                        else: sub_dp = None
                        utils.remove_glob(sub_xyz_head)

                        #if s_dp is not None: 
                        if sub_dp is not None and len(sub_dp) > 0:
                            try:
                                s_dp = np.concatenate((s_dp, sub_dp), axis = 0)
                            except: s_dp = sub_dp
                        #else: s_dp = sub_dp
                    utils.remove_glob(o_xyz, 'sub_{}*'.format(n))

            if len(s_dp) > 0:
                d_max = self.region_info[self.dem.name][4]
                #s_max = self.region_info[self.dem.name][5]
                s_dp = s_dp[s_dp[:,3] < d_max,:]
                s_dp = s_dp[s_dp[:,3] > 0,:]
                prox_err = s_dp[:,[2,3]]

                if last_ec_d is None:
                    last_ec_d = [0, 0.1, 0.2]
                    last_ec_diff = 10
                else: last_ec_diff = abs(last_ec_d[2] - last_ec_d[1])

                ec_d = self._err2coeff(prox_err[:50000000], coeff_guess=last_ec_d, dst_name=self.dem.name + '_prox', xa='distance')
                ec_diff = abs(ec_d[2] - ec_d[1])
                ec_l_diff = abs(last_ec_diff - ec_diff)

                # s_dp = s_dp[s_dp[:,4] < s_max,:]
                # slp_err = s_dp[:,[2,4]]
                # #print(slp_err)
                # #ec_s = self._err2coeff(slp_err[:50000000], coeff_guess=[0, 0.1, 0.2], dst_name = self.dem.name + '_slp', xa = 'slope')
                # ec_s = [0, 1, 2]
                # utils.echo_msg('{}\t{}\t{}\t{}\t{}\t{}'.format(sim, len(s_dp), ec_d, ec_d[2] - ec_d[1], ec_s, ec_s[2] - ec_s[1]))
                utils.echo_msg('{}\t{}\t{}\t{}'.format(sim, len(s_dp), ec_d, ec_l_diff))

                #if ec_d[2] < 0.0001: continue
                #if abs(ec_d[2] - ec_d[1]) > 2: continue
                if ec_d[0] == 0 and ec_d[1] == 0.1 and ec_d[2] == 0.2:
                    continue
                if sim >= int(self.sims):
                    break
                if abs(last_ec_diff - ec_diff) == 0:
                    break
                #if abs(last_ec_diff - ec_diff) < 0.001: break
                #if len(s_dp) >= int(self.region_info[self.dem.name][1] / 10): break
                last_ec_d = ec_d
                #else: utils.echo_msg('{}\t{}\t{}\t{}\t{}\t{}'.format(sim, len(s_dp), None, None, None, None))
            else:
                utils.echo_msg('{}\t{}\t{}\t{}'.format(sim, len(s_dp), None, None))
        _prog.end(status, 'performed {} SPLIT-SAMPLE simulations'.format(sim))
        return([ec_d])
        
    def run(self):

        s_dp = s_ds = None
        unc_out = {}
        zones = ['low-dens','mid-dens','high-dens','low-slp','mid-slp','high-slp']
        utils.echo_msg('running INTERPOLATION uncertainty module using {}...'.format(self.dem.mod))
        
        if self.prox is None:
            self._gen_prox()

        # if self.slope is None:
        #     self._gen_slope()
            
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

        # self.slp_percentile = demfun.percentile(self.slope, self.percentile)
        # self.slp_perc_33 = demfun.percentile(self.slope, 25)
        # self.slp_perc_66 = demfun.percentile(self.slope, 75)
        # self.slp_perc_100 = demfun.percentile(self.slope, 100)

        #self.region_info[self.dem.name] = [self.dem.region, g_max, num_sum, num_perc, self.prox_percentile, self.slp_percentile]
        self.region_info[self.dem.name] = [self.dem.region, g_max, num_sum, num_perc, self.prox_percentile] 
        for x in self.region_info.keys():
            utils.echo_msg('region: {}: {}'.format(x, self.region_info[x]))

        ## ==============================================
        ## chunk region into sub regions
        ## ==============================================
        chnk_inc = int((self.region_info[self.dem.name][1] / math.sqrt(g_max)) / self.region_info[self.dem.name][3])
        #chnk_inc = 250
        sub_regions = self.dem.region.chunk(self.dem.inc, chnk_inc)
        utils.echo_msg('chunked region into {} sub-regions @ {}x{} cells.'.format(len(sub_regions), chnk_inc, chnk_inc))

        ## ==============================================
        ## sub-region analysis
        ## ==============================================
        sub_zones = self._sub_region_analysis(sub_regions)

        ## ==============================================
        ## sub-region density and percentiles
        ## ==============================================
        s_dens = np.array([sub_zones[x][3] for x in sub_zones.keys()])
        s_5perc = np.percentile(s_dens, 5)
        s_dens = None
        utils.echo_msg('Sampling density for region is: {:.16f}'.format(s_5perc))

        ## ==============================================
        ## zone analysis / generate training regions
        ## ==============================================
        trainers = []
        t_perc = 95
        s_perc = 50

        for z, this_zone in enumerate(self._zones):
            #print(sub_zones) sub_zones[x][8] (with slope)
            tile_set = [sub_zones[x] for x in sub_zones.keys() if sub_zones[x][5] == self._zones[z]]
            if len(tile_set) > 0:
                d_50perc = np.percentile(np.array([x[3] for x in tile_set]), 50)
            else: continue
            t_trainers = [x for x in tile_set if x[3] < d_50perc or abs(x[3] - d_50perc) < 0.01]
            utils.echo_msg('possible {} training zones: {} @ MAX {}'.format(self._zones[z].upper(), len(t_trainers), d_50perc))
            trainers.append(t_trainers)

        utils.echo_msg('sorting training tiles by distance...')
        trains = self._regions_sort(trainers, verbose = False)
        tot_trains = len([x for s in trains for x in s])
        utils.echo_msg('sorted sub-regions into {} training tiles.'.format(tot_trains))
        utils.echo_msg('analyzed {} sub-regions.'.format(len(sub_regions)))

        ## ==============================================
        ## split-sample simulations and error calculations
        ## sims = max-simulations
        ## ==============================================
        if self.sims is None:
            self.sims = int(len(sub_regions)/tot_trains)
        
        ec_d = self._split_sample(trains, s_5perc)[0]

        ## ==============================================
        ## Save/Output results
        ## apply error coefficient to full proximity grid
        ## TODO: USE numpy/gdal instead!
        ## ==============================================
        utils.echo_msg('applying coefficient to PROXIMITY grid')
        if self.dem.gc['GMT'] is None:
            utils.run_cmd('gdal_calc.py -A {} --outfile {}_prox_unc.tif --calc "{}+({}*(A**{}))"'.format(self.prox, self.dem.name, 0, ec_d[1], ec_d[2]), verbose = True)
        else:
            math_cmd = 'gmt grdmath {} 0 AND ABS {} POW {} MUL {} ADD = {}_prox_unc.tif=gd+n-9999:GTiff\
            '.format(self.prox, ec_d[2], ec_d[1], 0, self.dem.name)
            utils.run_cmd(math_cmd, verbose = self.dem.verbose)
        if self.dem.epsg is not None: status = demfun.set_epsg('{}_prox_unc.tif'.format(self.dem.name), epsg=self.dem.epsg)
        utils.echo_msg('applied coefficient {} to PROXIMITY grid'.format(ec_d))

        # utils.echo_msg('applying coefficient to SLOPE grid')
        # if self.dem.gc['GMT'] is None:
        #     utils.run_cmd('gdal_calc.py -A {} --outfile {}_slp_unc.tif --calc "{}+({}*(A**{}))"'.format(self.slope, self.dem.name, 0, ec_s[1], ec_s[2]), verbose = True)
        # else:
        #     math_cmd = 'gmt grdmath {} 0 AND ABS {} POW {} MUL {} ADD = {}_slp_unc.tif=gd+n-9999:GTiff\
        #     '.format(self.slope, ec_s[2], ec_s[1], 0, self.dem.name)
        #     utils.run_cmd(math_cmd, verbose = self.dem.verbose)
        # if self.dem.epsg is not None: status = demfun.set_epsg('{}_prox_unc.tif'.format(self.dem.name), epsg=self.dem.epsg)
        # utils.echo_msg('applied coefficient {} to SLOPE grid'.format(ec_s))

        utils.remove_glob(self.prox)
        #utils.remove_glob(self.slope)

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

usage: {cmd} [OPTIONS] DATALIST

Options:
  -R, --region\t\tSpecifies the desired REGION;
\t\t\tThis can either be a GMT-style region ( -R xmin/xmax/ymin/ymax )
\t\t\tor an OGR-compatible vector file with regional polygons. 
\t\t\tIf a vector file is supplied it will search each region found therein.
\t\t\tIf omitted, use the region gathered from the data in DATALIST.
  -E, --increment\tGridding CELL-SIZE in native units or GMT-style increments.
\t\t\tappend :<inc> to resample the output to the given <inc>: -E.3333333s:.1111111s
  -F, --format\t\tOutput grid FORMAT. [GTiff]
  -M, --module\t\tDesired Waffles MODULE and options. (see available Modules below)
\t\t\tsyntax is -M module:mod_opt=mod_val:mod_opt1=mod_val1:...
  -O, --output-name\tBASENAME for all outputs.
  -P, --epsg\t\tHorizontal projection of data as EPSG code [4326]
  -X, --extend\t\tNumber of cells with which to EXTEND the REGION.
\t\t\tappend :<num> to extend the processing region: -X6:12
  -T, --filter\t\tFILTER the output DEM using one or multiple filters. <fltr:fltr_val:split_value=z>
\t\t\tAvailable filters:
\t\t\t1: perform a Gaussian filter at -T1:<factor>.
\t\t\t2: use a Cosine Arch Filter at -T2:<dist(km)> search distance.
\t\t\t3: Spike Filter at -T3:<stand-dev. threshhold>.
\t\t\tThe -T switch may be set multiple times to perform multiple filters.
\t\t\tAppend :split_value=<num> to only filter values below z-value <num>.
\t\t\te.g. -T1:10:split_value=0 to smooth bathymetry (z<0) using Gaussian filter
  -C, --clip\t\tCLIP the output to the clip polygon -C<clip_ply.shp:invert=False>
  -G, --wg-config\tA waffles config JSON file. If supplied, will overwrite all other options.
\t\t\tgenerate a waffles_config JSON file using the --config flag.

  -p, --prefix\t\tSet BASENAME (-O) to PREFIX (append inc/region/year info to output BASENAME).
  -r, --grid-node\tUse grid-node registration, default is pixel-node
  -w, --weights\t\tUse weights provided in the datalist to weight overlapping data.

  -c, --continue\tDon't clobber existing files.
  -q, --quiet\t\tLower verbosity to a quiet.

  --help\t\tPrint the usage text
  --config\t\tSave the waffles config JSON and major datalist
  --modules\t\tDisply the module descriptions and usage
  --version\t\tPrint the version information

Datalists and data formats:
  A datalist is a file that contains a number of datalist entries, 
  while an entry is a space-delineated line:
  `/path/to/data format weight data,meta,data`

Supported datalist formats: 
  {dl_formats}

Modules (see waffles --modules <module-name> for more info):
  {modules}

CIRES DEM home page: <http://ciresgroups.colorado.edu/coastalDEM>
""".format(cmd=os.path.basename(sys.argv[0]),
           dl_formats=dlim._datalist_fmts_short_desc(),
           modules=_waffles_module_short_desc(),
           wf_version=cudem.__version__)

def uncertainties_cli(argv = sys.argv):
    """run waffles from command-line

    See `waffles_cli_usage` for full cli options.
    """
    
    dls = []
    i_regions = []
    these_regions = []
    module = None
    wg_user = None
    want_prefix = False
    status = 0
    i = 1
    wg = {}
    wg['verbose'] = True
    wg['sample'] = None
    wg['fltr'] = []
    
    while i < len(argv):
        arg = argv[i]
        if arg == '--region' or arg == '-R':
            i_regions.append(str(argv[i + 1]))
            i += 1
        elif arg[:2] == '-R': i_regions.append(str(arg[2:]))
        elif arg == '--module' or arg == '-M':
            module = str(argv[i + 1])
            i += 1
        elif arg[:2] == '-M': module = str(arg[2:])
        elif arg == '--increment' or arg == '-E':
            incs = argv[i + 1].split(':')
            wg['inc'] = utils.str2inc(incs[0])
            if len(incs) > 1: wg['sample'] = utils.str2inc(incs[1])
            i = i + 1
        elif arg[:2] == '-E':
            incs = arg[2:].split(':')
            wg['inc'] = utils.str2inc(arg[2:].split(':')[0])
            if len(incs) > 1: wg['sample'] = utils.str2inc(incs[1])
        elif arg == '--outname' or arg == '-O':
            wg['name'] = argv[i + 1]
            i += 1
        elif arg[:2] == '-O': wg['name'] = arg[2:]
        elif arg == '--format' or arg == '-F':
            wg['fmt'] = argv[i + 1]
            i += 1
        elif arg[:2] == '-F': wg['fmt'] = arg[2:]
        elif arg == '--filter' or arg == '-T':
            wg['fltr'].append(argv[i + 1])
            i += 1
        elif arg[:2] == '-T': wg['fltr'].append(arg[2:])
        elif arg == '--extend' or arg == '-X':
            exts = argv[i + 1].split(':')
            wg['extend'] = utils.int_or(exts[0], 0)
            if len(exts) > 1: wg['extend_proc'] = utils.int_or(exts[1], 10)
            i += 1
        elif arg[:2] == '-X':
            exts = arg[2:].split(':')
            wg['extend'] = utils.int_or(exts[0], 0)
            if len(exts) > 1: wg['extend_proc'] = utils.int_or(exts[1], 10)
        elif arg == '--wg-config' or arg == '-G':
            wg_user = argv[i + 1]
            i += 1
        elif arg[:2] == '-G': wg_user = arg[2:]
        elif arg == '--clip' or arg == '-C':
            wg['clip'] = argv[i + 1]
            i = i + 1
        elif arg[:2] == '-C': wg['clip'] = arg[2:]
        elif arg == '--chunk' or arg == '-K':
            wg['chunk'] = utils.int_or(argv[i + 1], None)
            i = i + 1
        elif arg[:2] == '-K': wg['chunk'] = utils.int_or(arg[2:], None)
        elif arg == '--epsg' or arg == '-P':
            wg['epsg'] = utils.int_or(argv[i + 1], 4326)
            i = i + 1
        elif arg[:2] == '-P': wg['epsg'] = utils.int_or(arg[2:], 4326)
        
        elif arg == '-w' or arg == '--weights': wg['weights'] = 1
        elif arg == '-t' or arg == '--threads': want_threads = True
        elif arg == '-p' or arg == '--prefix': want_prefix = True
        elif arg == '-a' or arg == '--archive': wg['archive'] = True
        elif arg == '-m' or arg == '--mask': wg['mask'] = True
        elif arg == '-s' or arg == '--spat': wg['spat'] = True
        elif arg == '-c' or arg == '--continue': wg['clobber'] = False
        elif arg == '-r' or arg == '--grid-node': wg['node'] = 'grid'

        elif arg == '--quiet' or arg == '-q': wg['verbose'] = False
        elif arg == '--config': want_config = True
        elif arg == '--modules' or arg == '-m':
            try:
                if argv[i + 1] in waffles.WaffleFactory()._modules.keys():
                    sys.stderr.write(_waffles_module_long_desc({k: waffles.WaffleFactory()._modules[k] for k in (argv[i + 1],)}))
                else: sys.stderr.write(_waffles_module_long_desc(waffles.WaffleFactory()._modules))
            except: sys.stderr.write(_waffles_module_long_desc(waffles.WaffleFactory()._modules))
            sys.exit(0)
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
        else: dls.append(arg)
        i += 1

    if module is None:
        sys.stderr.write(uncertainties_cli_usage)
        utils.echo_error_msg('''must specify a waffles -M module.''')
        sys.exit(-1)
                
    ## ==============================================
    ## check the increment
    ## ==============================================
    if wg['inc'] is None:
        sys.stderr.write(uncertainties_cli_usage)
        utils.echo_error_msg('''must specify a gridding increment.''')
        sys.exit(-1)
    
    ## ==============================================
    ## set the datalists and names
    ## ==============================================
    wg['data'] = dls

    for i_region in i_regions:
        tmp_region = regions.Region().from_string(i_region)
        if tmp_region.valid_p():
            these_regions.append(tmp_region)
        else:
            tmp_region = regions.ogr_wkts(i_region)
            for i in tmp_region:
                if i.valid_p():
                    these_regions.append(i)
                    
    if len(these_regions) == 0:
        these_regions = []
        utils.echo_error_msg('failed to parse region(s), {}'.format(i_regions))
    else:
        if wg['verbose']: utils.echo_msg('parsed {} region(s)'.format(len(these_regions)))
    
    for i, this_region in enumerate(these_regions):
        wg['src_region'] = this_region
        if want_prefix or len(these_regions) > 1:
            wg['name'] = utils.append_fn(wg['name'], wg['src_region'], wg['sample'] if wg['sample'] is not None else wg['inc'])

        wg['mask'] = True
        print(module)
        wf = waffles.WaffleFactory(mod=module, **wg).acquire()
        if not wf.valid_p():
            wf.generate()

        #print(wf.mod)
        #print(wf.mod_args)
        i = InterpolationUncertainty(dem=wf).run()
        utils.echo_msg(wf)
    
# ## ==============================================
# ## testing
# ## ==============================================
# region = regions.Region().from_list([480000.01, 482999.99, 4395000.0, 4397999.99])
# epsg = 26913

# module = 'triangulate'
# module_args = ()

# inc = 1
# name='CO_SoPlatteRiver_MergeTriangulate1'

# datasets = [
#     'USGS_LPC_CO_SoPlatteRiver_Lot5_2013_13SDD480395_LAS_2015.xyz',
#     'USGS_LPC_CO_SoPlatteRiver_Lot5_2013_13SDD481395_LAS_2015.xyz',
#     'USGS_LPC_CO_SoPlatteRiver_Lot5_2013_13SDD480396_LAS_2015.xyz',
#     'USGS_LPC_CO_SoPlatteRiver_Lot5_2013_13SDD481396_LAS_2015.xyz',
# ]

# waffles_dem = waffles.WaffleFactory(
#     mod = module,
#     data = datasets,
#     src_region = region,
#     inc = inc,
#     name = name,
#     epsg = epsg,
# ).acquire(*module_args)

# waffles_dem.fn = 'CO_SoPlatteRiver_MergeTriangulate1.tif'
# waffles_dem.waffled = True

# print(waffles_dem.valid_p())

# i = InterpolationUncertainty(dem=waffles_dem, data_mask='CO_SoPlatteRiver_MergeTriangulate1_msk.tif')
# i.run()

# print(i.dem)
