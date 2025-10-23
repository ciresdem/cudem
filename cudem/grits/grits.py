### grits.py - GRId filTerS
##
## Copyright (c) 2024 - 2025 Regents of the University of Colorado
##
## Permission is hereby granted, free of charge, to any person obtaining a copy 
## of this software and associated documentation files (the "Software"), to deal 
## in the Software without restriction, including without limitation the rights 
## to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies 
## of the Software, and to permit persons to whom the Software is furnished to do so, 
## subject to the following conditions:
##
## The above copyright notice and this permission notice shall be included in all copies or
## substantial portions of the Software.
##
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
## INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR 
## PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE 
## FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
## ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
## SOFTWARE.
##
###############################################################################
### Commentary:
##
## filter grids using various methods
##
## Grits modules are sub-classes of the Grits class.
## Define a new sub-class to create a new DEM filter.
##
### Code:

import os
import numpy as np
from osgeo import gdal

from cudem import utils
from cudem import gdalfun

class Grits:
    """DEM Filtering.

    Filter DEMs using various filtering methods. 

    Define a sub-class to make a new grits filter.
    """
    
    def __init__(
            self,
            src_dem: str = None,
            dst_dem: str = None,
            band: int = 1,
            min_z: float = None,
            max_z: float = None,
            min_weight: float = None,
            max_weight: float = None,
            count_mask: any = None,
            weight_mask: any = None,
            uncertainty_mask: any = None,
            cache_dir: str = './',
            verbose: bool = True,
            params: dict = {},
            **kwargs: any
    ):
        self.src_dem = src_dem
        self.dst_dem = dst_dem
        self.band = utils.int_or(band, 1)
        self.min_z = utils.float_or(min_z)
        self.max_z = utils.float_or(max_z)
        self.min_weight = utils.float_or(min_weight)
        self.max_weight = utils.float_or(max_weight)
        self.count_mask = count_mask
        self.weight_mask = weight_mask
        self.uncertainty_mask = uncertainty_mask
        self.cache_dir = cache_dir
        self.verbose = verbose
        self.params = params
        self.kwargs = kwargs

        if self.dst_dem is None:
            if self.src_dem is not None:
                self.dst_dem = utils.make_temp_fn('{}_filtered.{}'.format(
                    utils.fn_basename2(self.src_dem),
                    utils.fn_ext(self.src_dem)
                ), temp_dir=self.cache_dir)
            else:
                self.dst_dem = 'grits_filtered.tif'

                
    def __call__(self):
        return(self.generate())

    
    def init_ds(self, src_ds: any = None):
        self.ds_config = gdalfun.gdal_infos(src_ds)
        self.ds_band = src_ds.GetRasterBand(self.band)
        self.gt = self.ds_config['geoT']

        if self.ds_band.GetNoDataValue() is None:
            self.ds_band.SetNoDataValue(self.ds_config['ndv'])
        
        ## setup the associted mask data (uncertainty, weights, counts)
        self.weight_is_fn = False
        self.weight_is_band = False
        self.uncertainty_is_fn = False
        self.uncertainty_is_band = False
        self.count_is_fn = False
        self.count_is_band = False
        if self.weight_mask is not None:
            self.weight_is_band = False
            self.weight_is_fn = False
            if utils.int_or(self.weight_mask) is not None:
                self.weight_is_band = True
                self.weight_mask = utils.int_or(self.weight_mask)
            elif os.path.exists(self.weight_mask):
                self.weight_is_fn = True
            else:
                self.weight_mask = None        

        if self.uncertainty_mask is not None:
            self.unc_is_band = False
            self.unc_is_fn = False
            if utils.int_or(self.uncertainty_mask) is not None:
                self.unc_is_band = True
                self.uncertainty_mask = utils.int_or(self.uncertainty_mask)
            elif os.path.exists(self.uncertainty_mask):
                self.unc_is_fn = True
            else:
                self.uncertainty_mask = None
                
        if self.count_mask is not None:
            self.cnt_is_band = False
            self.cnt_is_fn = False
            if utils.int_or(self.count_mask) is not None:
                self.cnt_is_band = True
                self.count_mask = utils.int_or(self.count_mask)
            elif os.path.exists(self.count_mask):
                self.cnt_is_fn = True
            else:
                self.count_mask = None

                
    def generate(self):
        if self.verbose:
            utils.echo_msg(
                f'filtering {self.src_dem} using {self}'
            )
            
        self.run()
        self.split_by_z()
        self.split_by_weight()
        return(self)        

    
    def run(self):
        raise(NotImplementedError)

    
    def copy_src_dem(self):
        with gdalfun.gdal_datasource(
                self.src_dem, update=False
        ) as src_ds:
            if src_ds is not None:
                src_infos = gdalfun.gdal_infos(src_ds)
                driver = gdal.GetDriverByName(src_infos['fmt'])
                copy_ds = driver.CreateCopy(
                    self.dst_dem, src_ds, 1, options=['COMPRESS=DEFLATE']
                )
            else:
                copy_ds = None
                
        return(copy_ds)

    
    def _density(self, src_arr):
        nonzero = np.count_nonzero(~np.isnan(src_arr))
        dd = nonzero / src_arr.size
        
        return(dd)

    
    def split_by_z(self):
        """Split the filtered DEM by z-value"""

        if self.max_z is not None or self.min_z is not None:
            utils.echo_msg(
                f'split by z:{self.min_z} {self.max_z}'
            )
            with gdalfun.gdal_datasource(self.src_dem) as src_ds:
                if src_ds is not None:
                    self.init_ds(src_ds)
                    elev_array = self.ds_band.ReadAsArray()
                    mask_array = np.zeros(
                        (self.ds_config['ny'], self.ds_config['nx'])
                    )
                    mask_array[elev_array == self.ds_config['ndv']] = np.nan
                    if self.min_z is not None:
                        mask_array[elev_array > self.min_z] = 1
                        if self.max_z is not None:
                            mask_array[elev_array > self.max_z] = 0
                        
                    elif self.max_z is not None:
                        mask_array[elev_array < self.max_z] = 1
                        if self.min_z is not None:
                            mask_array[elev_array < self.min_z] = 0
                        
                    elev_array[mask_array == 1] = 0

                    ## todo: all bands
                    with gdalfun.gdal_datasource(self.dst_dem, update=True) as s_ds:
                        if s_ds is not None:
                            #for b in range(1, s_ds.RasterCount+1):
                            s_band = s_ds.GetRasterBand(1)
                            s_array = s_band.ReadAsArray()
                            s_array = s_array * mask_array
                            smoothed_array = s_array + elev_array
                            elev_array = None
                            s_band.WriteArray(smoothed_array)
                            
        return(self)

    
    def split_by_weight(self):
        """Split the filtered DEM by z-value"""

        if self.max_weight is not None \
           or self.min_weight is not None:
            if self.weight_mask is not None:
                utils.echo_msg(
                    f'split by weight: {self.min_weight} {self.max_weight}'
                )
                with gdalfun.gdal_datasource(self.src_dem) as src_ds:
                    if src_ds is not None:
                        self.init_ds(src_ds)
                        elev_array = self.ds_band.ReadAsArray()

                        # uncertainty ds
                        weight_band = None
                        if self.weight_is_fn:
                            weight_ds = gdal.Open(self.weight_mask)
                            weight_band = weight_ds.GetRasterBand(1)
                        elif self.weight_is_band:
                            weight_band = src_ds.GetRasterBand(self.weight_mask)

                        if weight_band is not None:
                            weight_array = weight_band.ReadAsArray()
                            weight_array[(weight_array == self.ds_config['ndv'])] = 0

                        mask_array = np.zeros(
                            (self.ds_config['ny'], self.ds_config['nx'])
                        )
                        mask_array[elev_array == self.ds_config['ndv']] = np.nan
                        mask_array[weight_array == self.ds_config['ndv']] = np.nan

                        if self.min_weight is not None:
                            mask_array[weight_array > self.min_weight] = 1
                            if self.max_weight is not None:
                                mask_array[weight_array > self.max_weight] = 0

                        elif self.max_weight is not None:
                            mask_array[weight_array < self.max_weight] = 1
                            if self.min_weight is not None:
                                mask_array[weight_array < self.min_weight] = 0

                        elev_array[mask_array == 1] = 0

                        ## todo: all bands
                        with gdalfun.gdal_datasource(
                                self.dst_dem, update=True
                        ) as s_ds:
                            if s_ds is not None:
                                s_band = s_ds.GetRasterBand(1)
                                s_array = s_band.ReadAsArray()
                                s_array = s_array * mask_array
                                smoothed_array = s_array + elev_array
                                elev_array = None
                                s_band.WriteArray(smoothed_array)
        return(self)

    
    def get_outliers(self, in_array: any, percentile: float = 75,
                     k: float = 1.5, verbose: bool = False):
        """get the outliers from in_array based on the percentile

        https://en.wikipedia.org/wiki/Interquartile_range
        """

        if verbose:
            utils.echo_msg(
                f'input percentile: {percentile}'
            )

        # if np.isnan(percentile):
        #     percentile = 75
            
        if percentile < 0:
            percentile = 0
            
        if percentile > 100:
            percentile = 100

        max_percentile = percentile
        #min_percentile = 100 - percentile
        min_percentile = percentile-50

        if min_percentile < 0:
            min_percentile = 0

        if verbose:
            utils.echo_msg(
                f'percentiles: {min_percentile}>>{max_percentile}'
            )

        if np.all(np.isnan(in_array)):
            upper_limit = np.nan
            lower_limit = np.nan
        else:
            perc_max = np.nanpercentile(in_array, max_percentile)
            perc_min = np.nanpercentile(in_array, min_percentile)
            iqr_p = (perc_max - perc_min) * k
            upper_limit = perc_max + iqr_p
            lower_limit = perc_min - iqr_p

        return(upper_limit, lower_limit)
            
### End
