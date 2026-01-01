### uncertainty.py
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
##
## uncertainty.py is part of CUDEM
##
## Permission is hereby granted, free of charge, to any person obtaining a copy 
## of this software and associated documentation files (the "Software"), to deal 
## in the Software without restriction, including without limitation the rights 
## to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies 
## of the Software, and to permit persons to whom the Software is furnished to do so, 
## subject to the following conditions:
##
## The above copyright notice and this permission notice shall be included in all
## copies or substantial portions of the Software.
##
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
## INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR 
## PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE 
## FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
## ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
## SOFTWARE.
##
### Commentary:
##
### Code:

import os
import math
import numpy as np
from osgeo import gdal

from cudem import utils
from cudem import gdalfun
from cudem.waffles.waffles import Waffle

class WafflesUncertainty(Waffle):
    """Calculate cell-level interpolation uncertainty

    -----------
    Parameters:

    waffles_module (str): waffles module string
    percentile (int): max percentile
    sims (int): number of split-sample simulations
    chnk_lvl (int): the 'chunk-level'
    max_sample (int): the maximum sample errors
    max_errors (int): the maximum accumulated errors
    accumulate (bool): accumulate errors
    """
    
    def __init__(
            self,
            waffles_module='IDW',
            percentile=95,
            sims=10,
            chnk_lvl=None,
            max_sample=None,
            max_errors=5000000,
            accumulate=False,
            **kwargs):

        ## parse the waffles module
        self.waffles_module_args = {}
        tmp_waffles = Waffle()
        for kpam, kval in kwargs.items():
            if kpam not in tmp_waffles.__dict__:
                self.waffles_module_args[kpam] = kval

        for kpam, kval in self.waffles_module_args.items():
            del kwargs[kpam]
            
        super().__init__(**kwargs)
        self.waffles_module = waffles_module
        self.percentile = utils.float_or(percentile, 95)
        self.sims = sims
        self.max_sample = max_sample
        self.chnk_lvl = chnk_lvl
        self.accumulate = accumulate
        self.max_errors = max_errors

        ## set up the accumulated errors file
        self.prox_errs = '{}_errs.dat.gz'.format(self.waffles_module.split(':')[0])
        self.prox_errs_local = self.prox_errs
        if not os.path.exists(self.prox_errs):
            if os.path.exists(os.path.join(utils.CUDEM_DATA, self.prox_errs)):
                self.prox_errs = os.path.join(utils.CUDEM_DATA, self.prox_errs)
            else:
                utils.touch(self.prox_errs)
                self.accumulate = True
        
        self._zones = ['LD0','LD1','LD2','MD0','MD1','MD2','HD0', 'HD1', 'HD2']
        self.prox = None
        self.slope = None

        
    def _mask_analysis(self, src_gdal, region=None):
        """scan the mask raster and gather some infos...

        returns the number of filled cells, the total number of cells 
        and the percent of total cells filled.
        """
       
        ds_config = gdalfun.gdal_infos(src_gdal)
        if region is not None:
            srcwin = region.srcwin(
                ds_config['geoT'], ds_config['nx'], ds_config['ny']
            )
        else:
            srcwin = (0, 0, ds_config['nx'], ds_config['ny'])
          
        ds_arr = src_gdal.GetRasterBand(1).ReadAsArray(
            srcwin[0], srcwin[1], srcwin[2], srcwin[3]
        )
        ds_arr[ds_arr == ds_config['ndv']] = np.nan
        ds_arr[~np.isnan(ds_arr)] = 1
        msk_sum = np.nansum(ds_arr)
        msk_max = float(srcwin[2] * srcwin[3])
        msk_perc = float((msk_sum / msk_max) * 100.)
        dst_arr = None

        return(msk_sum, msk_max, msk_perc)

    
    def _prox_analysis(self, src_gdal, region = None, band = 1):
        """scan the proximity raster and gather some infos...

        returns the percentile (self.percentile) of values in the srcwin
        """
        
        ds_config = gdalfun.gdal_infos(src_gdal)
        if region is not None:
            srcwin = region.srcwin(
                ds_config['geoT'], ds_config['nx'], ds_config['ny']
            )
        else:
            srcwin = (0, 0, ds_config['nx'], ds_config['ny'])
            
        ds_arr = src_gdal.GetRasterBand(band).ReadAsArray(*srcwin).astype(float)
        ds_arr[ds_arr == ds_config['ndv']] = np.nan
        prox_perc = np.nanpercentile(ds_arr, self.percentile)
        dst_arr = None

        return(prox_perc)

    
    def _generate_proximity_raster(self, out_prox = None):
        """
        generate a proximity grid from the data mask raster
        
        returns the output proximity grid's fn
        """
        
        if out_prox is None:
            out_prox = utils.make_temp_fn(f'{self.waffles_module}_prox.tif')
            # out_prox = utils.make_temp_fn(
            #     '{}_prox.tif'.format(self.params['mod_args']['waffles_module'])
            # )

        if self.verbose:
            utils.echo_msg(
                f'generating proximity grid {out_prox}...'
            )
            
        gdalfun.gdal_proximity(self.stack, out_prox, distunits='PIXEL')
        if self.dst_srs is not None:
            gdalfun.gdal_set_srs(out_prox, self.dst_srs)

        return(out_prox)

    
    def _generate_slope_raster(self, out_slope = None):
        """
        generate a slope grid from the elevation raster
        
        returns the output slope grid's fn
        """

        if out_slope is None:
            out_slope = utils.make_temp_fn(f'{self.waffles_module}_slope.tif')
            #     '{}_slope.tif'.format(
            #         self.params['mod_args']['waffles_module']
            #     )
            # )

        if self.verbose:
            utils.echo_msg(
                f'generating slope grid {out_slope}...'
            )
            
        gdalfun.gdal_slope(self.stack, out_slope)
        if self.dst_srs is not None:
            gdalfun.gdal_set_srs(out_slope, self.dst_srs)

        return(out_slope)

    
    def _generate_interpolated_src_raster(self):
        """
        generate an interpolated source uncertainty raster
        """

        from cudem.waffles.waffles import WaffleFactory
        
        #gdalfun.gdal_get_array(self.stack, band=4)
        src_unc_name = '{}_src_unc'.format(self.name)
        src_unc_surface = WaffleFactory(
            mod='nearest',
            data=['{},200:band_no=4,1'.format(self.stack)],
            src_region=self._proc_region(),
            xinc=self.xinc,
            yinc=self.yinc,
            name=src_unc_name,
            node='pixel',
            want_weight=False,
            want_uncertainty=False,
            dst_srs=self.dst_srs,
            srs_transform=self.srs_transform,
            clobber=True,
            verbose=self.verbose,
            keep_auxiliary=False
        )._acquire_module()
        src_unc_surface.initialize()
        src_unc_surface.generate()
        return(src_unc_surface.fn)

    
    def _regions_sort(self, trainers, t_num = 25):
        """sort regions (trainers is a list of regions) by 
        distance from one another; 

        a region is a list: [xmin, xmax, ymin, ymax].

        -----------
        Parameters:

        trainers (region-list): a list of regions to sort
        t_num (int): total output number of regions
        
        -----------
        Returns:

        the sorted region-list
        """

        train_sorted = []
        for z, train in enumerate(trainers):
            train_d = []
            np.random.shuffle(train)
            train_total = len(train)
            while True:
                # if self.verbose:
                #     utils.echo_msg_inline(
                #         'sorting training tiles [{}]'.format(len(train))
                #     )
                    
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
                
            ## uncomment to print out the sorted regions...
            #if self.verbose:
            #    utils.echo_msg(' '.join([x[0].format('gmt') for x in train_d[:t_num]]))
                
            train_sorted.append(train_d)
            
        #if self.verbose:
        #    utils.echo_msg_inline('sorting training tiles [OK]\n')
            
        return(train_sorted)

    
    def _select_split(self, o_xyz, sub_region, sub_bn):
        """split an xyz file into an inner and outer region.

        -----------
        Parameters:

        o_xyz (fn): input xyz file-name to split
        sub_region(regions.Region()): the region to split the xyz with.
        sub_bn (str): the basename of the output split xyz files.
        """
        
        out_inner = '{}_inner.xyz'.format(sub_bn)
        out_outer = '{}_outer.xyz'.format(sub_bn)
        xyz_ds = xyzfile.XYZFile(
            fn=o_xyz,
            data_format=168,
            src_region=sub_region
        ).initialize()
        with open(out_inner, 'w') as sub_inner:
            xyz_ds.dump_xyz_direct(dst_port=sub_inner)
            
        xyz_ds.invert_region = True
        with open(out_outer, 'w') as sub_outer:
            xyz_ds.dump_xyz_direct(dst_port=sub_outer)
            
        return([out_inner, out_outer])    

    
    def _sub_region_analysis(self, sub_regions):
        """analyze a list of sub-regions and assign them to 
        various zones (self._zones)

        -----------
        Parameters:

        sub_regions (list): a list of regions to analyze

        -----------
        Returns:

        the sub-zones extracted from the sub-regions
        """
        
        sub_zones = {}
        stack_ds = gdal.Open(self.stack)
        prox_ds = gdal.Open(self.prox)
        slp_ds = gdal.Open(self.slope)
        with utils.ccp(
                total=len(sub_regions),
                desc=f'analyzing {len(sub_regions)} sub-regions',
                leave=self.verbose
        ) as pbar:
            for sc, sub_region in enumerate(sub_regions):
                pbar.update()
                s_sum, s_g_max, s_perc = self._mask_analysis(
                    stack_ds, region=sub_region
                )
                if s_sum == 0:
                    continue

                p_perc = self._prox_analysis(prox_ds, region=sub_region)
                slp_perc = self._prox_analysis(slp_ds, region=sub_region)
                zone = None
                
                ## assign the region to the zone based on the density/slope
                if p_perc <= self.prox_perc_33:# or abs(p_perc - self.prox_perc_33) < 0.01:
                    if slp_perc <= self.slp_perc_33:# or abs(slp_perc - self.slp_perc_33) < 0.01:
                        zone = self._zones[6]
                    elif slp_perc <= self.slp_perc_66:# or abs(slp_perc - self.slp_perc_66) < 0.01:
                        zone = self._zones[7]
                    else:
                        zone = self._zones[8]
                elif p_perc <= self.prox_perc_66:# or abs(p_perc - self.prox_perc_66) < 0.01:
                    if slp_perc <= self.slp_perc_33:# or abs(slp_perc - self.slp_perc_33) < 0.01:
                        zone = self._zones[3]
                    elif slp_perc <= self.slp_perc_66:# or abs(slp_perc - self.slp_perc_66) < 0.01:
                        zone = self._zones[4]
                    else:
                        zone = self._zones[5]
                else:
                    if slp_perc <= self.slp_perc_33:# or abs(slp_perc - self.slp_perc_33) < 0.01:
                        zone = self._zones[0]
                    elif slp_perc <= self.slp_perc_66:# or abs(slp_perc - self.slp_perc_66) < 0.01:
                        zone = self._zones[1]
                    else:
                        zone = self._zones[2]

                if zone is not None:
                    sub_zones[sc + 1] = [sub_region, s_g_max, s_sum, s_perc, p_perc, zone]
            
        stack_ds = prox_ds = slp_ds = None
        return(sub_zones)

    
    def _split_sample(self, trainers, perc):
        """split-sample simulations and error calculations

        -----------
        Parameters:

        trainers (list): a list of training regions (sorted)
        perc (float): sampling density

        -----------
        Returns:

        distance error array
        """

        from cudem.waffles.waffles import WaffleFactory
        
        last_ec_d = None
        s_dp = []
        status = 0
        trains = self._regions_sort(trainers)
        all_trains = [x for s in trains for x in s[:5]]
        tot_trains = len(all_trains)
        
        with utils.ccp(
                desc='performing SPLIT-SAMPLE simulation',
                leave=False,
                total=tot_trains
        ) as pbar:
            for n, sub_region in enumerate(all_trains):
                ## perform split-sample analysis on each training region.
                pbar.update()
                ss_samp = perc
                this_region = sub_region[0].copy()
                if sub_region[3] < ss_samp:
                   ss_samp = sub_region[3]

                ## extract the xyz data for the region from the DEM
                o_xyz = utils.make_temp_fn('{}_{}.xyz'.format(self.name, n))
                with gdalfun.gdal_datasource(self.stack) as ds:
                   ds_config = gdalfun.gdal_infos(ds)
                   b_region = this_region.copy()
                   b_region.buffer(pct=20, x_inc=self.xinc, y_inc=self.yinc)
                   srcwin = b_region.srcwin(
                       ds_config['geoT'], ds_config['nx'], ds_config['ny']
                   )

                   ## TODO: extract weights here as well...
                   with open(o_xyz, 'w') as o_fh:
                       for xyz in gdalfun.gdal_parse(ds, srcwin=srcwin):
                           xyz.dump(dst_port=o_fh)

                if os.stat(o_xyz).st_size == 0:
                    continue

                ## split the xyz data to inner/outer; outer is
                ## the data buffer, inner will be randomly sampled
                s_inner, s_outer = self._select_split(
                    o_xyz, this_region, utils.make_temp_fn('sub_{}'.format(n))
                )
                if os.stat(s_inner).st_size == 0:
                    utils.echo_warning_msg('no inner points, cont.')
                    continue

                if os.stat(s_outer).st_size == 0:
                    utils.echo_warning_msg('no outer points, cont.')
                    continue

                sub_xyz = np.loadtxt(s_inner, ndmin=2, delimiter=' ')                        
                ss_len = len(sub_xyz)

                ## determine the sampling density
                #sx_cnt = int(sub_region[2] * (ss_samp / 100.)) if ss_samp is not None else ss_len-1
                sx_cnt = int(sub_region[1] * (ss_samp / 100.)) + 1
                ##sx_cnt = int(ss_len * (ss_samp / 100.))
                sx_cnt = 1 if sx_cnt < 1 or sx_cnt >= ss_len else sx_cnt
                sub_xyz_head = utils.make_temp_fn('sub_{}_head_{}.xyz'.format(n, sx_cnt))
                np.random.shuffle(sub_xyz)
                np.savetxt(sub_xyz_head, sub_xyz[sx_cnt:], '%f', ' ')

                ## generate the random-sample DEM
                this_mod = '{}:{}'.format(
                    self.waffles_module, factory.dict2args(self.waffles_module_args)
                )
                kwargs = self.params['kwargs']
                kwargs['name'] = utils.make_temp_fn('sub_{}'.format(n))
                kwargs['data'] = [s_outer, sub_xyz_head]
                kwargs['src_region'] = b_region
                kwargs['want_uncertainty'] = False
                kwargs['want_sm'] = False
                kwargs['verbose'] = False
                kwargs['clobber'] = True
                this_waffle = WaffleFactory(mod=this_mod, **kwargs)._acquire_module()
                this_waffle.initialize()
                wf = this_waffle.generate()
                if not WaffleDEM(
                        wf.fn,
                        cache_dir=self.cache_dir,
                        verbose=self.verbose
                ).initialize().valid_p():
                    continue

                ## generate the random-sample data PROX and SLOPE
                sub_prox = '{}_prox.tif'.format(wf.name)
                gdalfun.gdal_proximity(wf.stack, sub_prox, distunits='PIXEL')

                ## Calculate the random-sample errors
                ## todo: account for source uncertainty (rms with xyz?)
                sub_xyd = gdalfun.gdal_query(sub_xyz[:sx_cnt], wf.fn, 'xyd')
                sub_dp = gdalfun.gdal_query(sub_xyd, sub_prox, 'zg')
                utils.remove_glob('{}*'.format(sub_xyz_head))
                
                if sub_dp is not None and len(sub_dp) > 0:
                    try:
                        s_dp = np.vstack((s_dp, sub_dp))
                    except:
                        s_dp = sub_dp

                utils.remove_glob(
                    '{}*'.format(wf.stack),
                    '{}*'.format(o_xyz),
                    s_inner,
                    s_outer,
                    wf.fn
                )
                s_dp_m = []
                ## bin the error data
                if s_dp is not None and len(s_dp) > 0:
                    err_count = len(s_dp)
                    ds = np.unique(s_dp[:,1])
                    for d in ds:
                        arr=np.array(
                            [(True if x == d else False) for x in s_dp[:,1]]
                        )
                        if arr.any():
                            arr_count = np.count_nonzero(arr)
                            err_perc = (arr_count / err_count)
                            d_err_count = int(self.max_errors * err_perc)
                            err_sum = np.histogram(
                                s_dp[:,0][arr],
                                d_err_count,
                                weights=s_dp[:,0][arr]
                            )[0]
                            err_cnt = np.histogram(s_dp[:,0][arr], d_err_count)[0]
                            err_sum = err_sum[np.nonzero(err_cnt)]
                            err_cnt = err_cnt[np.nonzero(err_cnt)]
                            d_errs = err_sum/err_cnt
                            d_dist = np.full((d_errs.size, 1), d)
                            dist_errs = np.hstack(
                                (d_errs.reshape((d_errs.size, 1)), d_dist)
                            )
                            if len(s_dp_m) == 0:
                                s_dp_m = np.array(dist_errs)
                            else:
                                s_dp_m = np.vstack((s_dp_m, dist_errs))
        s_dp = np.array(s_dp_m)
        return(s_dp)

    
    def get_accumulated_coefficient(self):
        """load the distance/error points from the acuumulated file,
        calculate the error coefficient and return both.
        """
        
        s_dp = []        
        if os.path.exists(self.prox_errs):
            s_dp = np.loadtxt(self.prox_errs)
            
        utils.echo_msg(
            'loaded {} errors from {}'.format(
                len(s_dp), self.prox_errs
            )
        )
        pre_ec_d = [0, .1, .2]
        if len(s_dp) > 1:
            max_dist = np.nanpercentile(s_dp[:,1], 95)
            pre_ec_d = utils._err2coeff(
                s_dp[s_dp[:,1] <= max_dist],
                self.percentile,
                coeff_guess=pre_ec_d
            )

        return(pre_ec_d, s_dp)

    
    def apply_coefficient_(self, ec_d, want_multi_band=False):
        """apply the error coefficeint `ec_d` to the proximity 
        raster and add it to self.stack band 4 (uncertainty).
        """
        
        if self.prox is None:
            self.prox = self._generate_proximity_raster(
                f'{self.name}_u.tif'
            )

        src_unc_raster = self._generate_interpolated_src_raster()
            
        if self.verbose:
            utils.echo_msg(
                f'applying coefficient {ec_d} to PROXIMITY grid {self.prox}'
            )
            
        with gdalfun.gdal_datasource(self.prox, update=True) as prox_ds:
            prox_inf = gdalfun.gdal_infos(prox_ds)
            prox_band = prox_ds.GetRasterBand(1)
            prox_arr = prox_band.ReadAsArray().astype(float)
            
            with gdalfun.gdal_datasource(src_unc_raster) as src_unc_ds:
                src_unc_inf = gdalfun.gdal_infos(src_unc_ds)
                src_unc_band = src_unc_ds.GetRasterBand(1)
                src_unc_arr = src_unc_band.ReadAsArray().astype(float)
                src_unc_arr[src_unc_arr == src_unc_inf['ndv']] = np.nan
                out_arr = ec_d[0] + ec_d[1] * (prox_arr**ec_d[2])
                out_arr = np.sqrt(np.power(out_arr, 2), np.power(src_unc_arr, 2))
                
                #out_arr[~np.isnan(src_unc_arr)] = unc_arr[~np.isnan(src_unc_arr)]                
                #with gdalfun.gdal_datasource(self.stack) as stack_ds:
                #    unc_inf = gdalfun.gdal_infos(stack_ds, band=4)
                #    unc_band = stack_ds.GetRasterBand(4)
                #    unc_arr = unc_band.ReadAsArray()
                # unc_arr[unc_arr == unc_inf['ndv']] = np.nan
                # out_arr = ec_d[0] + ec_d[1] * (prox_arr**ec_d[2])
                # out_arr[~np.isnan(unc_arr)] = unc_arr[~np.isnan(unc_arr)]
                
        unc_out = gdalfun.gdal_write(
            out_arr, '{}.{}'.format(self.name, 'tif'), self.ds_config
        )[0]
        if self.dst_srs is not None:
            status = gdalfun.gdal_set_srs(self.prox, src_srs=self.dst_srs)

        if self.verbose:
            utils.echo_msg(
                f'applied coefficient {ec_d} to PROXIMITY grid'
            )
            
        return(unc_out)

    
    def apply_coefficient(self, ec_d, fmt = 'GTiff'):
        """apply the error coefficeint `ec_d` to the proximity 
        raster and add it to self.stack band 4 (uncertainty).

        output a multi-banded raster that includes the source uncertainty from
        self.stack, the proximity from data grid, the interpolation uncertainty
        and the final TVU.
        """

        unc_out = '{}.{}'.format(self.name, 'tif')
        gdt = gdal.GDT_Float32
        driver = gdal.GetDriverByName(fmt)
        xcount, ycount, dst_gt = self.region.geo_transform(
            x_inc=self.xinc, y_inc=self.yinc, node='grid'
        )
        if xcount <= 0 or ycount <=0:
            utils.echo_error_msg(
                (f'could not create grid of {xcount}x{ycount} '
                 f'cells with {self.xinc}/{self.yinc} increments '
                 f'on region: {self.region}')
            )
            sys.exit(-1)
        
        if os.path.exists(unc_out):
            status = driver.Delete(unc_out)
            if status != 0:
                utils.remove_glob('{}*'.format(unc_out))
        
        dst_ds = driver.Create(
            unc_out,
            xcount,
            ycount,
            4,
            gdt,
            options=['COMPRESS=LZW',
                     'PREDICTOR=2',
                     'TILED=YES',
                     'BIGTIFF=YES']
        )

        if dst_ds is None:
            utils.echo_error_msg(
                'failed to create uncertainty grid {} {} {} {} {}...'.format(
                    out_file, xcount, ycount, gdt, fmt
                )
            )
            sys.exit(-1)

        dst_ds.SetGeoTransform(dst_gt)
        unc_bands = {
            'tvu': dst_ds.GetRasterBand(1),
            'src_uncertainty': dst_ds.GetRasterBand(2),
            'interpolation_uncertainty': dst_ds.GetRasterBand(3),
            'proximity': dst_ds.GetRasterBand(4)
        }
        unc_data = {
            'tvu': None,
            'src_uncertainty': None,
            'interpolation_uncertainty': None,
            'proximity': None
        }
        
        for key in unc_bands.keys():
            unc_bands[key].SetNoDataValue(np.nan)
            unc_bands[key].SetDescription(key)

        #dst_inf = gdalfun.gdal_infos(dst_ds)
        if self.prox is None:
            self.prox = self._generate_proximity_raster(
                f'{self.name}_u.tif'
            )

        src_unc_raster = self._generate_interpolated_src_raster()
        if self.verbose:
            utils.echo_msg(
                f'applying coefficient {ec_d} to PROXIMITY grid {self.prox}'
            )
            
        with gdalfun.gdal_datasource(self.prox, update=True) as prox_ds:
            prox_inf = gdalfun.gdal_infos(prox_ds)
            srcwin = self.region.srcwin(
                prox_inf['geoT'],
                prox_ds.RasterXSize,
                prox_ds.RasterYSize,
                node='grid'
            )
            prox_band = prox_ds.GetRasterBand(1)
            prox_arr = prox_band.ReadAsArray().astype(float)
            unc_bands['proximity'].WriteArray(
                prox_arr[srcwin[0]:srcwin[0]+srcwin[2],
                         srcwin[1]:srcwin[1]+srcwin[3]]
            )
            with gdalfun.gdal_datasource(src_unc_raster) as src_unc_ds:
                src_unc_inf = gdalfun.gdal_infos(src_unc_ds)
                src_unc_band = src_unc_ds.GetRasterBand(1)
                src_unc_arr = src_unc_band.ReadAsArray().astype(float)
                src_unc_arr[src_unc_arr == src_unc_inf['ndv']] = np.nan

                unc_bands['src_uncertainty'].WriteArray(
                    src_unc_arr[srcwin[0]:srcwin[0]+srcwin[2],
                                srcwin[1]:srcwin[1]+srcwin[3]]
                )
                
                interp_arr = ec_d[0] + ec_d[1] * (prox_arr**ec_d[2])
                interp_arr[prox_arr == 0] = 0
                unc_bands['interpolation_uncertainty'].WriteArray(
                    interp_arr[srcwin[0]:srcwin[0]+srcwin[2],
                            srcwin[1]:srcwin[1]+srcwin[3]]
                )
                
                out_arr = np.sqrt(np.power(interp_arr, 2) + np.power(src_unc_arr, 2))
                unc_bands['tvu'].WriteArray(
                    out_arr[srcwin[0]:srcwin[0]+srcwin[2],
                            srcwin[1]:srcwin[1]+srcwin[3]]
                )

            # with gdalfun.gdal_datasource(self.stack) as stack_ds:
            #     unc_inf = gdalfun.gdal_infos(stack_ds, band=4)
            #     unc_stack_band = stack_ds.GetRasterBand(4)
            #     unc_arr = unc_stack_band.ReadAsArray()
            #     unc_arr[unc_arr == unc_inf['ndv']] = np.nan
            #     unc_arr[np.isnan(unc_arr)] = 0
            #     unc_bands['src_uncertainty'].WriteArray(
            #         unc_arr[srcwin[0]:srcwin[0]+srcwin[2],
            #                 srcwin[1]:srcwin[1]+srcwin[3]]
            #     )
            #     out_arr = ec_d[0] + ec_d[1] * (prox_arr**ec_d[2])
            #     out_arr[prox_arr == 0] = 0
            #     unc_bands['interpolation_uncertainty'].WriteArray(
            #         out_arr[srcwin[0]:srcwin[0]+srcwin[2],
            #                 srcwin[1]:srcwin[1]+srcwin[3]]
            #     )
            #     unc_arr[unc_arr == 0] = np.nan
            #     out_arr[~np.isnan(unc_arr)] = unc_arr[~np.isnan(unc_arr)]
            #     unc_bands['tvu'].WriteArray(
            #         out_arr[srcwin[0]:srcwin[0]+srcwin[2],
            #                 srcwin[1]:srcwin[1]+srcwin[3]]
            #     )
                
        #unc_out = gdalfun.gdal_write(out_arr, '{}.{}'.format(self.name, 'tif'), self.ds_config)[0]
        
        if self.dst_srs is not None:
            status = gdalfun.gdal_set_srs(
                self.prox, src_srs=self.dst_srs
            )

        if self.verbose:
            utils.echo_msg(
                f'applied coefficient {ec_d} to PROXIMITY grid'
            )

        dst_ds = None
        utils.remove_glob(f'{src_unc_raster}*')
        return(unc_out)

    
    def run(self):
        """run the waffles uncertainty module"""
        
        if self.waffles_module.split(':')[0] not in [
                'IDW', 'linear', 'cubic', 'nearest', 'gmt-surface',
                'gmt-triangulate', 'gmt-nearneighbor', 'mbgrid',
                'gdal-linear', 'gdal-nearest', 'gdal-average',
                'gdal-invdst', 'flatten', 'cudem'
        ]:
            utils.echo_warning_msg(
                'cannot perform interpolation uncertainty estimation with {}'.format(
                self.waffles_module.split(':')[0]
                )
            )
            if self.verbose:
                utils.echo_msg('extracting uncertainty band from stack...')

            # band 4 is the uncertainty band in stacks
            gdalfun.gdal_extract_band(self.stack, self.fn, band=4) 
            with gdalfun.gdal_datasource(self.fn, update=True) as unc_ds:
                unc_band = unc_ds.GetRasterBand(1)
                unc_band.SetDescription('tvu')
                
            return(self)
        
        s_dp = s_ds = None
        unc_out = {}
        if self.verbose:
            utils.echo_msg(
                f'running UNCERTAINTY module using {self.waffles_module}...'
            )
            utils.echo_msg(
                'using {}; accumulate is {}'.format(
                    self.prox_errs, self.accumulate
                )
            )
            
        if self.prox is None:
            self.prox = self._generate_proximity_raster()

        if self.slope is None:
            self.slope = self._generate_slope_raster()

        pre_ec_d, s_dp = self.get_accumulated_coefficient() 
        if len(s_dp) <= 1:
            self.accumulate = True

        if not self.accumulate:
            unc_out = self.apply_coefficient(pre_ec_d)
            utils.remove_glob(self.slope)
            utils.remove_glob(self.prox)
            return(unc_out, 0)                        
        else:
            ## region and der. analysis
            self.region_info = {}
            with gdalfun.gdal_datasource(self.stack) as tmp_ds:
                num_sum, g_max, num_perc = self._mask_analysis(tmp_ds)

            self.prox_percentile = gdalfun.gdal_percentile(self.prox, self.percentile)
            self.prox_perc_33 = gdalfun.gdal_percentile(self.prox, 25)
            self.prox_perc_66 = gdalfun.gdal_percentile(self.prox, 75)
            self.prox_perc_100 = gdalfun.gdal_percentile(self.prox, 100)

            self.slp_percentile = gdalfun.gdal_percentile(self.slope, self.percentile)
            self.slp_perc_33 = gdalfun.gdal_percentile(self.slope, 25)
            self.slp_perc_66 = gdalfun.gdal_percentile(self.slope, 75)
            self.slp_perc_100 = gdalfun.gdal_percentile(self.slope, 100)

            self.region_info[self.name] = [self.region,
                                           g_max,
                                           num_sum,
                                           num_perc,
                                           self.prox_percentile]
            for x in self.region_info.keys():
                utils.echo_msg(
                    f'region: {x}: {self.region_info[x]}'
                )

            ## chunk region into sub regions
            chnk_inc = int((num_sum / math.sqrt(g_max)) / num_perc) * 2
            chnk_inc = chnk_inc if chnk_inc > 10 else 10
            utils.echo_msg('chunk inc is: {}'.format(chnk_inc))

            sub_regions = self.region.chunk(self.xinc, chnk_inc)
            utils.echo_msg(
                'chunked region into {} sub-regions @ {}x{} cells.'.format(
                    len(sub_regions), chnk_inc, chnk_inc
                )
            )

            ## sub-region analysis
            sub_zones = self._sub_region_analysis(sub_regions)

            ## sub-region density and percentiles
            s_dens = np.array([sub_zones[x][3] for x in sub_zones.keys()])
            s_5perc = np.percentile(s_dens, 5)
            s_dens = None
            utils.echo_msg(
                'Sampling density for region is: {:.16f}'.format(num_perc)
            )

            ## zone analysis / generate training regions
            trainers = []
            t_perc = 95
            s_perc = 50
            for z, this_zone in enumerate(self._zones):
                tile_set = [
                    sub_zones[x] for x in sub_zones.keys() if sub_zones[x][5] == self._zones[z]
                ]
                if len(tile_set) > 0:
                    d_50perc = np.percentile(np.array([x[3] for x in tile_set]), 50)
                else:
                    continue

                t_trainers = [
                    x for x in tile_set if x[3] < d_50perc or abs(x[3] - d_50perc) < 0.01
                ]
                utils.echo_msg(
                    'possible {} training zones: {} @ MAX {}'.format(
                        self._zones[z].upper(), len(t_trainers), d_50perc
                    )
                )
                trainers.append(t_trainers)

            utils.echo_msg(
                'analyzed {} sub-regions.'.format(len(sub_regions))
            )

            ## split-sample simulations and error calculations
            ## sims = max-simulations
            if self.sims is None:
                self.sims = int(len(sub_regions)/tot_trains)

            if self.max_sample is None:
                self.max_sample = int(
                    (self.region_info[self.name][1] \
                     - self.region_info[self.name][2]) \
                    * .005
                )

            sim = 0
            max_dist = gdalfun.gdal_percentile(self.prox, 95)
            if self.verbose:
                utils.echo_msg(
                    'max sample is {}, max sims is {}'.format(
                        self.max_sample, self.sims
                    )
                )
                utils.echo_msg(
                    'pre ec_d is {}'.format(pre_ec_d)
                )
                utils.echo_msg(
                    'performing at least {} simulations, looking for {} errors'.format(
                        self.sims, self.max_sample
                    )
                )                
                utils.echo_msg(
                    'max distance is {}'.format(max_dist)
                )
                utils.echo_msg(
                    'simulation\terrors\tmean-error\tproximity-coeff'
                )
                
            while True:
                sim += 1                
                ## run the split-sample simulation(s)
                sample_dp = self._split_sample(trainers, num_perc)
                if len(s_dp) == 0:
                    s_dp = sample_dp
                else:
                    s_dp = np.vstack((s_dp, sample_dp))

                err_count = len(s_dp)
                if err_count == 0:
                    utils.echo_error_msg(
                        'did not gather any errors, check configuration'
                    )
                    break

                ## bin the error data
                ds = np.unique(s_dp[:,1])
                s_dp_m = None
                for d in ds:
                    arr=np.array([(True if x == d else False) for x in s_dp[:,1]])
                    if arr.any():
                        arr_count = np.count_nonzero(arr)
                        err_perc = (arr_count / err_count)
                        d_err_count = int(self.max_errors * err_perc)
                        err_sum = np.histogram(
                            s_dp[:,0][arr],
                            d_err_count,
                            weights=s_dp[:,0][arr]
                        )[0]
                        err_cnt = np.histogram(s_dp[:,0][arr], d_err_count)[0]
                        err_sum = err_sum[np.nonzero(err_cnt)]
                        err_cnt = err_cnt[np.nonzero(err_cnt)]
                        d_errs = err_sum/err_cnt
                        d_dist = np.full((d_errs.size, 1), d)
                        dist_errs = np.hstack((d_errs.reshape((d_errs.size, 1)), d_dist))

                        if s_dp_m is None:
                            s_dp_m = np.array(dist_errs)
                        else:
                            s_dp_m = np.vstack((s_dp_m, dist_errs))

                s_dp = np.array(s_dp_m)
                if self.accumulate:
                    np.savetxt(self.prox_errs_local, s_dp, '%f', ' ')

                max_dist = np.nanpercentile(s_dp[:,1], 95)
                # if self.verbose:
                #     utils.echo_msg('max distance is {}'.format(max(s_dp[:,1])))
                #     utils.echo_msg('max distance 95th percentile is {}'.format(max_dist))

                ec_d = utils._err2coeff(
                    s_dp[s_dp[:,1] <= max_dist],
                    num_perc,
                    coeff_guess=pre_ec_d,
                    plots=True
                )
                pre_ec_d = ec_d
                if self.verbose:
                    utils.echo_msg(
                        '{}\t{}\t{}\t{}'.format(
                            sim, len(s_dp), np.mean(s_dp, axis=0)[0], ec_d
                        )
                    )

                ## continue if we got back the default err coeff
                if ec_d[0] == 0 and ec_d[1] == 0.1 and ec_d[2] == 0.2:
                    continue

                ## continue if we haven't reached max_sample
                #if len(s_dp) < self.max_sample:
                #    continue

                ## break if we gathered enough simulation errors
                if sim >= int(self.sims): 
                    break

            ## Save/Output results, apply error coefficient to full proximity grid
            unc_out = self.apply_coefficient(ec_d)
            utils.remove_glob(self.slope, self.prox)

        return(self)


### End
