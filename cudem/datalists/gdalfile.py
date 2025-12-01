### inf.py - DataLists IMproved
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
##
## inf.py is part of CUDEM
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
###############################################################################
### Commentary:
##
### Examples:
##
### TODO:
##
### Code:

import os
import re
import numpy as np
from tqdm import tqdm

from osgeo import gdal

from cudem import utils
from cudem import regions
from cudem import gdalfun
from cudem.datalists.dlim import ElevationDataset

class GDALFile(ElevationDataset):
    """providing a GDAL raster dataset parser.

    Process/Parse GDAL supported raster files.
    See GDAL for more information regarding supported formats.
    
    generate_inf - generate an inf file for the RASTER data
    yield_xyz - yield the RASTER data as xyz
    yield_array - yield the RASTER data as an array

    -----------
    Parameters:

    mask: raster dataset to use as MASK dataset OR a band number
    weight_mask: raster dataset to use as weights (per-cell) OR a band number
    otherwise will use single value (weight) from superclass.
    uncertainty_mask: raster dataset to use as uncertainty (per-cell) OR a band number
    uncertainty_mask_to_meter: conversion value for uncertainty data to meter
    otherwise will use a single value (uncertainty) from superclass.
    open_options: GDAL open_options for raster dataset
    sample: sample method to use in resamplinig
    check_path: check to make sure path exists
    super_grid: Force processing of a supergrid (BAG files) (True/False)
    band_no: the band number of the elevation data
    remove_flat: remove flattened data from the input
    node: force node registration of either 'pixel' or 'grid'
    """
    
    def __init__(self,
                 weight_mask=None,
                 uncertainty_mask=None,
                 x_band=None,
                 y_band=None,
                 uncertainty_mask_to_meter=1,
                 open_options=None,
                 sample=None,
                 check_path=True,
                 super_grid=False,
                 band_no=1,
                 remove_flat=False,
                 node=None,
                 resample_and_warp=True,
                 yield_chunk=True,
                 **kwargs):
        super().__init__(**kwargs)
        # associated raster file/band holding weight data
        self.weight_mask = weight_mask 
        # associated raster file/band holding uncertainty data
        self.uncertainty_mask = uncertainty_mask 
        # uncertainty data factor to meters
        self.uncertainty_mask_to_meter = uncertainty_mask_to_meter
        # GDAL open-options
        self.open_options = open_options
        # GDAL resampling method
        self.sample = sample
        # check for the path to the input data file
        self.check_path = check_path
        # input has super_grids (force their use)
        self.super_grid = super_grid
        # band number holding elevation data
        self.band_no = band_no
        # temporary elevation band
        self.tmp_elev_band = None
        # temporary uncertainty band
        self.tmp_unc_band = None
        # temporary weight band
        self.tmp_weight_band = None
        # the GDAL dataset object
        self.src_ds = None
        # remove flattened data from input
        self.remove_flat = remove_flat 
        self.flat_removed = False
        # band holding x values
        self.x_band = x_band
        # band holding y values
        self.y_band = y_band
        # input is 'pixel' or 'grid' registered (force)
        self.node = node 
        self.resample_and_warp = resample_and_warp
        self.yield_chunk = yield_chunk
        
        if self.fn.startswith('http') \
           or self.fn.startswith('/vsicurl/') \
           or self.fn.startswith('BAG'):
            self.check_path = False

        if self.valid_p() and self.src_srs is None:
            if self.infos.src_srs is None:
                self.src_srs = self.init_srs(self.fn)
            else:
                self.src_srs = self.infos.src_srs
                
    def destroy_ds(self):
        self.src_ds = None

        
    def init_srs(self, src_ds):
        """initialize the srs from the gdal file.

        try to split the horizontal and vertical and them 
        combine them...
        """
        
        if self.src_srs is None:
            return(gdalfun.gdal_get_srs(src_ds))
        else:
            return(self.src_srs)

        
    def generate_inf(self):
        with gdalfun.gdal_datasource(self.fn) as src_ds:
            if src_ds is not None:
                if self.src_srs is None:
                    self.infos.src_srs = self.init_srs(src_ds)
                else:
                    self.infos.src_srs = self.src_srs

                ds_infos = gdalfun.gdal_infos(src_ds)
                this_region = regions.Region(
                    src_srs=self.src_srs
                ).from_geo_transform(
                    geo_transform=ds_infos['geoT'],
                    x_count=ds_infos['nx'],
                    y_count=ds_infos['ny']
                )

                # if scan:
                #    zr = src_ds.GetRasterBand(
                #        utils.int_or(self.band_no, 1)
                #    ).ComputeRasterMinMax()

                ##gdalfun.gdal_polygonize(src_ds, dst_layer, verbose=True)
                #this_region.zmin, this_region.zmax = zr[0], zr[1]
                #self.infos.minmax = this_region.export_as_list(include_z=True)
                self.infos.minmax = this_region.export_as_list()
                self.infos.wkt = this_region.export_as_wkt()
                self.infos.numpts = ds_infos['nb']

        return(self.infos)

    
    def get_srcwin(self, gt, x_size, y_size, node='grid'):
        if self.region is not None:
            if self.invert_region:
                srcwin = (
                    0, 0, x_size, y_size, node
                )
            else:
                if self.transform['transformer'] is not None:
                    if self.transform['trans_region'] is not None \
                       and self.transform['trans_region'].valid_p(
                            check_xy = True
                    ):
                        srcwin = self.transform['trans_region'].srcwin(
                            gt, x_size, y_size, node
                        )
                    else:
                        srcwin = (
                            0, 0, x_size, y_size, node
                        )

                else:
                    srcwin = self.region.srcwin(
                        gt, x_size, y_size, node
                    )

        else:
            srcwin = (
                0, 0, x_size, y_size
            )
            
        return(srcwin)

    
    def yield_points(self):
        """initialize the raster dataset

        if x/y incs are set, will warp raster to that resolution.
        """

        if self.fn is None or (self.check_path and not os.path.exists(self.fn)):
            utils.echo_warning_msg(f'{self.fn} doesn\'t exist!')
            return(None)    

        ## apply any open_options that are specified
        try:
            self.open_options = self.open_options.split('/')
        except AttributeError:
            self.open_options = self.open_options
        except:
            self.open_options = None
            
        ## set up any transformations and other options
        #inf_region = regions.Region().from_string(self.infos.wkt)
        self.sample_alg = self.sample if self.sample is not None else self.sample_alg
        self.dem_infos = gdalfun.gdal_infos(self.fn)
        if self.node is None:
            self.node = gdalfun.gdal_get_node(self.fn, 'pixel')
            
        #self.resample_and_warp = True        
        if (self.x_inc is None and self.y_inc is None) or self.region is None:
            self.resample_and_warp = False
        else:
            ## only perform 'upsamples' if the input raster is of higher resolution than the
            ## target self.x_inc, etc. then treat the data as points and do any 'downsampling'
            ## in stacks.
            if self.transform['trans_region'] is not None:
                raster_is_higher_res = np.prod(
                    self.transform['trans_region'].geo_transform(
                        x_inc=self.dem_infos['geoT'][1]
                    )[:2]
                ) > np.prod(
                    self.transform['trans_region'].geo_transform(
                        x_inc=self.x_inc
                    )[:2]
                )

            else:
                raster_is_higher_res = False

            if raster_is_higher_res:
                self.resample_and_warp = False
                
        if self.node == 'grid':
            self.resample_and_warp = False

        ndv = utils.float_or(gdalfun.gdal_get_ndv(self.fn), -9999)
        if self.region is not None:
            self.warp_region = self.region.copy()
        else:
            if self.transform['transformer'] is not None:
                self.warp_region = self.transform['trans_region'].copy()
            else:
                self.warp_region = self.inf_region.copy()
            
        tmp_elev_fn = utils.make_temp_fn(self.fn, temp_dir=self.cache_dir)
        tmp_unc_fn = utils.make_temp_fn(self.fn, temp_dir=self.cache_dir)
        tmp_weight_fn = utils.make_temp_fn(self.fn, temp_dir=self.cache_dir)
        if self.remove_flat:
            grits_filter = grits.grits.GritsFactory(
                mod='flats',
                src_dem=self.fn,
                cache_dir=self.cache_dir,
                verbose=True
            )._acquire_module()
            if grits_filter is not None:
                grits_filter()                
                self.fn = grits_filter.dst_dem
                self.flat_removed = True
            
        ## resample/warp src gdal file to specified x/y inc/transformer respectively
        if self.resample_and_warp:
            if self.transform['transformer'] is not None:
                self.transform['transformer'] = None

            if self.sample_alg == 'auto':
                if self.stack_mode == 'min':
                    self.sample_alg = 'min'
                elif self.stack_mode == 'max':
                    self.sample_alg = 'max'
                elif not raster_is_higher_res:
                    self.sample_alg = 'bilinear'
                else:
                    self.sample_alg = 'average'

            tmp_ds = self.fn
            if self.open_options is not None:
                self.src_ds = gdal.OpenEx(
                    self.fn, open_options=self.open_options
                )
                if self.src_ds is None:
                    self.src_ds = gdal.Open(self.fn)
            else:
                self.src_ds = gdal.Open(self.fn)

            if self.src_ds is None:
                return(None)

            ## Sample/Warp
            ## resmaple and/or warp dataset based on target srs and x_inc/y_inc
            ## doing this in MEM has a bug, fix if able
            ## extract necessary bands before warping,
            ## gadlwarp does not work with multiband rasters!
            tmp_warp = utils.make_temp_fn(f'{tmp_ds}', temp_dir=self.cache_dir)
            in_bands = self.src_ds.RasterCount
            src_ds_config = gdalfun.gdal_infos(self.src_ds)
            src_gt = src_ds_config['geoT']

            if in_bands > 1:
                ## the srcwin for to extract data
                if self.transform['trans_region'] is not None:
                    srcwin_region = self.transform['trans_region'].copy()
                elif self.region is not None:
                    srcwin_region = self.region.copy()
                else:
                    srcwin_region = None

                if srcwin_region is not None:
                    srcwin = srcwin_region.srcwin(
                        src_gt,
                        src_ds_config['nx'],
                        src_ds_config['ny'],
                        node='pixel'
                    )
                else:
                    srcwin = None

                self.tmp_elev_band = tmp_elev_fn
                if self.verbose:
                    utils.echo_msg(
                        (f'extracting elevation data from {self.fn} '
                        f'to {self.tmp_elev_band}')
                    )
                    
                self.tmp_elev_band, status = gdalfun.gdal_extract_band(
                    self.src_ds, self.tmp_elev_band,
                    band=self.band_no,
                    exclude=[],
                    srcwin=srcwin,
                    inverse=False
                )
                tmp_ds = self.tmp_elev_band
                if tmp_ds is None:
                    return(None)
                
                ## band is now 1!
                self.band_no = 1
                if utils.int_or(self.uncertainty_mask) is not None:
                    self.tmp_unc_band = tmp_unc_fn
                    if self.verbose:
                        utils.echo_msg(
                            (f'extracting uncertainty mask from {self.fn} '
                             f'to {self.tmp_unc_band}')
                        )
                        
                    self.tmp_unc_band, status = gdalfun.gdal_extract_band(
                        self.src_ds, self.tmp_unc_band,
                        band=self.uncertainty_mask,
                        exclude=[],
                        srcwin=srcwin,
                        inverse=False
                    )
                    self.uncertainty_mask = self.tmp_unc_band

                if utils.int_or(self.weight_mask) is not None:                    
                    self.tmp_weight_band = tmp_weight_fn
                    if self.verbose:
                        utils.echo_msg(
                            (f'extracting weight mask from {self.fn} '
                             f'to {self.tmp_weight_band}')
                        )
                        
                    self.tmp_weight_band, status = gdalfun.gdal_extract_band(
                        self.src_ds, self.tmp_weight_band,
                        band=self.weight_mask,
                        exclude=[],
                        srcwin=srcwin,
                        inverse=False
                    )
                    self.weight_mask = self.tmp_weight_band

            warp_ = gdalfun.sample_warp(
                tmp_ds, tmp_warp, self.x_inc, self.y_inc,
                src_srs=self.transform['src_horz_crs'].to_proj4()
                if self.transform['src_horz_crs'] is not None \
                else None,
                dst_srs=self.transform['dst_horz_crs'].to_proj4() \
                if self.transform['dst_horz_crs'] is not None \
                else None,                
                src_region=self.warp_region,
                sample_alg=self.sample_alg,
                ndv=ndv,
                verbose=self.verbose,
                co=["COMPRESS=DEFLATE", "TILED=YES"]
            )[0]
            tmp_ds = warp_ = None
            
            ## the following seems to be redundant...
            warp_ds = gdal.Open(tmp_warp)
            if warp_ds is not None:
                ## clip wapred ds to warped srcwin
                warp_ds_config = gdalfun.gdal_infos(warp_ds)
                gt = warp_ds_config['geoT']
                srcwin = self.warp_region.srcwin(
                    gt, warp_ds.RasterXSize, warp_ds.RasterYSize, node='grid'
                )
                dst_gt = (gt[0] + (srcwin[0] * gt[1]),
                          gt[1],
                          0.,
                          gt[3] + (srcwin[1] * gt[5]),
                          0.,
                          gt[5])
                out_ds_config = gdalfun.gdal_set_infos(
                    srcwin[2], srcwin[3], srcwin[2] * srcwin[3], dst_gt,
                    warp_ds_config['proj'], warp_ds_config['dt'],
                    warp_ds_config['ndv'], warp_ds_config['fmt'],
                    None, None
                )

                in_bands = warp_ds.RasterCount
                self.src_ds = gdalfun.gdal_mem_ds(out_ds_config, bands=in_bands)
                if self.src_ds is not None:
                    for band in range(1, in_bands+1):
                        this_band = self.src_ds.GetRasterBand(band)
                        this_band.WriteArray(
                            warp_ds.GetRasterBand(band).ReadAsArray(*srcwin)
                        )
                        self.src_ds.FlushCache()

                warp_ds = None
                utils.remove_glob(tmp_warp)
        else:
            if self.open_options:
                self.src_ds = gdal.OpenEx(
                    self.fn, open_options=self.open_options
                )
                if self.src_ds is None:
                    if self.verbose:
                        utils.echo_warning_msg(
                            ('could not open file using open '
                             f'options {self.open_options}')
                        )
                        
                    self.src_ds = gdal.Open(self.fn)
            else:
                self.src_ds = gdal.Open(self.fn)
                
        self.src_dem_infos = gdalfun.gdal_infos(self.src_ds)
        if self.src_ds is None:
            if self.verbose:
                utils.echo_error_msg(
                    f'could not load raster file {self.fn}'
                )

        band = self.src_ds.GetRasterBand(utils.int_or(self.band_no, 1))
        gt = self.src_ds.GetGeoTransform()
        ndv = utils.float_or(band.GetNoDataValue())
        mask_band = None
        weight_band = None
        uncertainty_band = None
        nodata = ['{:g}'.format(-9999), 'nan', float('nan')]            
        if ndv is not None:
            nodata.append('{:g}'.format(ndv))

        band_data = count_data = weight_data = mask_data = None
                
        src_dem_x_inc = self.src_dem_infos['geoT'][1]
        src_dem_y_inc = -1*self.src_dem_infos['geoT'][5]
        src_dem_region = regions.Region().from_geo_transform(
            self.src_dem_infos['geoT'],
            self.src_dem_infos['nx'],
            self.src_dem_infos['ny']
        )

        #######################################################################
        ## todo: always warp these to src_ds
        ## weight mask, each cell should have the corresponding weight
        ## weight_mask can either be a seperate gdal file or a band number
        ## corresponding to the appropriate band in src_ds
        #######################################################################
        if self.weight_mask is not None:
            if utils.int_or(self.weight_mask) is not None:
                weight_band = self.src_ds.GetRasterBand(int(self.weight_mask))
            elif os.path.exists(self.weight_mask):
                # some numbers now return true here (file-descriptors),
                # check for int first!
                src_weight = gdalfun.sample_warp(
                    self.weight_mask, None, src_dem_x_inc, src_dem_y_inc,
                    src_srs=self.transform['src_horz_crs'].to_proj4() \
                    if self.transform['src_horz_crs'] is not None \
                    else None,
                    dst_srs=self.transform['dst_horz_crs'].to_proj4() \
                    if self.transform['dst_horz_crs'] is not None \
                    else None,
                    src_region=src_dem_region,
                    sample_alg=self.sample_alg,
                    ndv=ndv,
                    verbose=self.verbose,
                    co=["COMPRESS=DEFLATE", "TILED=YES"]
                )[0]
                weight_band = src_weight.GetRasterBand(1)

            else:
                utils.echo_warning_msg(
                    f'could not load weight mask {self.weight_mask}'
                )
                weight_band = None

        #######################################################################
        ## uncertainty mask, each cell should have the corresponding uncertainty
        ## uncertainty_mask can either be a seperate gdal file or a band number
        ## corresponding to the appropriate band in src_ds
        #######################################################################
        if self.uncertainty_mask is not None:
            if utils.int_or(self.uncertainty_mask):
                uncertainty_band = self.src_ds.GetRasterBand(
                    int(self.uncertainty_mask)
                )
            elif os.path.exists(self.uncertainty_mask):
                src_uncertainty = gdalfun.sample_warp(
                    self.uncertainty_mask, None, src_dem_x_inc, src_dem_y_inc,
                    src_srs=self.transform['src_horz_crs'].to_proj4() \
                    if self.transform['src_horz_crs'] is not None \
                    else None,
                    dst_srs=self.transform['dst_horz_crs'].to_proj4() \
                    if self.transform['dst_horz_crs'] is not None \
                    else None,
                    src_region=src_dem_region,
                    sample_alg='bilinear',
                    ndv=ndv,
                    verbose=self.verbose,
                    co=["COMPRESS=DEFLATE", "TILED=YES"]
                )[0]
                uncertainty_band = src_uncertainty.GetRasterBand(1)
            else:
                utils.echo_warning_msg(
                    f'could not load uncertainty mask {self.uncertainty_mask}'
                )
                uncertainty_band = None

        #######################################################################
        ## uncertainty from the vertical transformation
        #######################################################################
        if self.transform['trans_fn_unc'] is not None:
            trans_uncertainty = gdalfun.sample_warp(
                self.transform['trans_fn_unc'], None, src_dem_x_inc, src_dem_y_inc,
                src_srs='+proj=longlat +datum=WGS84 +ellps=WGS84',
                dst_srs=self.transform['dst_horz_crs'].to_proj4() \
                if self.transform['dst_horz_crs'] is not None \
                else None,
                src_region=src_dem_region,
                sample_alg='bilinear',
                ndv=ndv,
                verbose=self.verbose,
                co=["COMPRESS=DEFLATE", "TILED=YES"]
            )[0]

            if uncertainty_band is not None:
                trans_uncertainty_band = trans_uncertainty.GetRasterBand(1)
                trans_uncertainty_arr = trans_uncertainty_band.ReadAsArray()
                uncertainty_arr = uncertainty_band.ReadAsArray()
                uncertainty_arr *= self.uncertainty_mask_to_meter
                uncertainty_arr = np.sqrt(uncertainty_arr**2 + trans_uncertainty_arr**2)
                uncertainty_band.WriteArray(uncertainty_arr)
                trans_uncertainty_band = None
            else:
                uncertainty_band = trans_uncertainty.GetRasterBand(1)

        if self.yield_chunk:
            #######################################################################
            ## parse through the data / chunks
            #######################################################################
            for srcwin in utils.yield_srcwin(
                    n_size=(self.src_ds.RasterYSize, self.src_ds.RasterXSize),
                    n_chunk=4000,
                    msg=f'chunking {self.fn} @ {self.src_ds.RasterXSize}/{self.src_ds.RasterYSize}',
                    verbose=False
            ):                
                band_data = band.ReadAsArray(*srcwin).astype(float)
                if ndv is not None and not np.isnan(ndv):
                    band_data[band_data == ndv] = np.nan

                if np.all(np.isnan(band_data)):
                    continue

                ## weights
                if weight_band is not None:
                    weight_data = weight_band.ReadAsArray(*srcwin)
                    weight_ndv = float(weight_band.GetNoDataValue())
                    if not np.isnan(weight_ndv):
                        weight_data[weight_data==weight_ndv] = np.nan
                else:
                    weight_data = np.ones(band_data.shape)

                ## uncertainty
                if uncertainty_band is not None:
                    uncertainty_data = uncertainty_band.ReadAsArray(*srcwin)
                    uncertainty_ndv = float(uncertainty_band.GetNoDataValue())
                    if not np.isnan(uncertainty_ndv):
                        uncertainty_data[uncertainty_data==uncertainty_ndv] = np.nan

                    if self.transform['trans_fn_unc'] is None:
                        uncertainty_data *= self.uncertainty_mask_to_meter

                else:
                    uncertainty_data = np.zeros(band_data.shape)

                ## convert grid array to points
                if self.x_band is None and self.y_band is None:
                                        
                    x_precision = 12#len(str(gt[0]).split('.')[-1])
                    y_precision = 12#len(str(gt[3]).split('.')[-1])
                    #while True:
                    # 'self.node to 'pixel', this breaks if set to
                    # 'grid' even if 'grid-node'
                    geo_x_origin, geo_y_origin = utils._pixel2geo(
                        srcwin[0], srcwin[1], gt,
                        node='pixel',
                        x_precision=x_precision,
                        y_precision=y_precision
                    ) 
                    geo_x_end, geo_y_end = utils._pixel2geo(
                        srcwin[0] + srcwin[2], srcwin[1] + srcwin[3], gt,
                        node='grid',
                        x_precision=x_precision,
                        y_precision=y_precision
                    )
                    lon_array = np.arange(geo_x_origin, geo_x_end, gt[1])
                    lat_array = np.arange(geo_y_origin, geo_y_end, gt[5])

                    lon_data = np.tile(lon_array, (band_data.shape[0], 1))
                    lat_data = np.tile(lat_array[:,None], (1, band_data.shape[1]))

                    try:
                        assert lon_data.shape == lat_data.shape
                        assert lon_data.shape == band_data.shape
                        dataset = np.column_stack(
                            (lon_data.flatten(), lat_data.flatten(), band_data.flatten(),
                             weight_data.flatten(), uncertainty_data.flatten())
                        )
                    except Exception as e:
                        utils.echo_error_msg(e)
                        pass

                else:
                    lon_band = self.src_ds.GetRasterBand(self.x_band)
                    lon_array = lon_band.ReadAsArray(*srcwin).astype(float)
                    lon_array[np.isnan(band_data)] = np.nan

                    lat_band = self.src_ds.GetRasterBand(self.y_band)
                    lat_array = lat_band.ReadAsArray(*srcwin).astype(float)
                    lat_array[np.isnan(band_data)] = np.nan
                    dataset = np.column_stack(
                        (lon_array[0], lat_array[0], band_data[0],
                         weight_data[0], uncertainty_data[0])
                    )

                points = np.rec.fromrecords(dataset, names='x, y, z, w, u')
                points =  points[~np.isnan(points['z'])]
                dataset = band_data = weight_data = None
                uncertainty_data = lat_array = lon_array = None
                utils.remove_glob(tmp_elev_fn, tmp_unc_fn, tmp_weight_fn)
                yield(points)
        else:
            #######################################################################
            ## parse through the data / scanline
            ## todo: option to parse by chunk
            #######################################################################
            srcwin = self.get_srcwin(
                gt, self.src_ds.RasterXSize, self.src_ds.RasterYSize,
                node=self.node
            )
            for y in range(srcwin[1], (srcwin[1] + srcwin[3]), 1):
                band_data = band.ReadAsArray(
                    srcwin[0], y, srcwin[2], 1
                ).astype(float)
                if ndv is not None and not np.isnan(ndv):
                    band_data[band_data == ndv] = np.nan

                if np.all(np.isnan(band_data)):
                    continue

                ## weights
                if weight_band is not None:
                    weight_data = weight_band.ReadAsArray(
                        srcwin[0], y, srcwin[2], 1
                    )
                    weight_ndv = float(weight_band.GetNoDataValue())
                    if not np.isnan(weight_ndv):
                        weight_data[weight_data==weight_ndv] = np.nan
                else:
                    weight_data = np.ones(band_data.shape)

                ## uncertainty
                if uncertainty_band is not None:
                    uncertainty_data = uncertainty_band.ReadAsArray(
                        srcwin[0], y, srcwin[2], 1
                    )
                    uncertainty_ndv = float(uncertainty_band.GetNoDataValue())
                    if not np.isnan(uncertainty_ndv):
                        uncertainty_data[uncertainty_data==uncertainty_ndv] = np.nan

                    if self.transform['trans_fn_unc'] is None:
                        uncertainty_data *= self.uncertainty_mask_to_meter

                else:
                    uncertainty_data = np.zeros(band_data.shape)

                ## convert grid array to points
                if self.x_band is None and self.y_band is None:
                    x_precision = 12#len(str(gt[0]).split('.')[-1])
                    y_precision = 12#len(str(gt[3]).split('.')[-1])
                    while True:
                        # 'self.node to 'pixel', this breaks if set to
                        # 'grid' even if 'grid-node'
                        geo_x_origin, geo_y_origin = utils._pixel2geo(
                            srcwin[0], y, gt,
                            node='pixel',
                            x_precision=x_precision,
                            y_precision=y_precision
                        ) 
                        geo_x_end, geo_y_end = utils._pixel2geo(
                            srcwin[0] + srcwin[2], y, gt,
                            node='grid',
                            x_precision=x_precision,
                            y_precision=y_precision
                        )
                        lon_array = np.arange(geo_x_origin, geo_x_end, gt[1])

                        if lon_array.shape == band_data[0].shape:
                            break
                        else:
                            if x_precision < 0 or y_precision < 0:
                                break

                            x_precision -= 1
                            y_precision -= 1

                    lat_array = np.zeros((lon_array.shape))
                    lat_array[:] = geo_y_origin
                    # dataset = np.column_stack(
                    #     (lon_array, lat_array, band_data[0],
                    #      weight_data[0], uncertainty_data[0])
                    # )

                    try:
                        assert lon_array.shape == lat_array.shape
                        assert lon_array.shape == band_data[0].shape
                        dataset = np.column_stack(
                            (lon_array, lat_array, band_data[0],
                             weight_data[0], uncertainty_data[0])
                        )
                    except Exception as e:
                        utils.echo_error_msg(e)
                        pass
                else:
                    lon_band = self.src_ds.GetRasterBand(self.x_band)
                    lon_array = lon_band.ReadAsArray(
                        srcwin[0], y, srcwin[2], 1
                    ).astype(float)
                    lon_array[np.isnan(band_data)] = np.nan

                    lat_band = self.src_ds.GetRasterBand(self.y_band)
                    lat_array = lat_band.ReadAsArray(
                        srcwin[0], y, srcwin[2], 1
                    ).astype(float)
                    lat_array[np.isnan(band_data)] = np.nan
                    dataset = np.column_stack(
                        (lon_array[0], lat_array[0], band_data[0],
                         weight_data[0], uncertainty_data[0])
                    )

                points = np.rec.fromrecords(dataset, names='x, y, z, w, u')
                points =  points[~np.isnan(points['z'])]
                dataset = band_data = weight_data = None
                uncertainty_data = lat_array = lon_array = None
                utils.remove_glob(tmp_elev_fn, tmp_unc_fn, tmp_weight_fn)
                yield(points)
            
        src_uncertainty = src_weight = trans_uncertainty = self.src_ds = None
        # delete the filtered dem...
        if self.remove_flat and self.flat_removed:
            utils.remove_glob(self.fn)

            
class BAGFile(ElevationDataset):
    """providing a BAG raster dataset parser.

    Process supergrids at native resolution if they
    exist, otherwise process as normal grid.

    generate_inf - generate an inf file for the BAG data
    init_srs - discover the srs of the BAG file
    parse - parse the gdal datasets out of the BAG file

    -----------
    Parameters:

    explode (bool): Explode the BAG and process each super grid seperately.
    force_vr (bool): Force VR processing (if BAG file has bad header info)
    vr_strategy (str): VR strategy to use (MIN/MAX/AUTO)
    """

    def __init__(self,
                 explode = False,
                 force_vr = False,
                 vr_resampled_grid = False,
                 vr_interpolate = True,
                 vr_strategy = 'MIN',
                 min_weight = 0,
                 **kwargs):
        super().__init__(**kwargs)
        self.explode = explode
        self.min_weight = utils.float_or(min_weight, 0)
        self.force_vr = force_vr
        self.vr_resampled_grid = vr_resampled_grid
        self.vr_interpolate = vr_interpolate
        if self.vr_interpolate:
            self.vr_resampled_grid = False
            
        self.vr_strategy = vr_strategy
        if self.src_srs is None:
            self.src_srs = self.init_srs()


    def init_srs(self):
        if self.src_srs is None:
            src_horz, src_vert = gdalfun.split_srs(
                gdalfun.gdal_get_srs(self.fn), as_epsg=False
            )
            if src_horz is None and src_vert is None:
                return(None)
            
            if src_vert is None:
                src_vert = 5866
            # if 'MSL' in src_vert:
            #     src_vert = '5703'

            self.src_srs = gdalfun.combine_epsgs(
                src_horz, src_vert, name='BAG Combined'
            )
            return(self.src_srs)
        else:
            return(self.src_srs)

        
    def generate_inf(self):
        if self.src_srs is None:
            self.infos.src_srs = self.init_srs()
        else:
            self.infos.src_srs = self.src_srs

        with gdalfun.gdal_datasource(self.fn) as src_ds:
            if src_ds is not None:
                ds_infos = gdalfun.gdal_infos(src_ds)
                this_region = regions.Region(src_srs=self.src_srs).from_geo_transform(
                    geo_transform=ds_infos['geoT'],
                    x_count=ds_infos['nx'],
                    y_count=ds_infos['ny']
                )
                # bag band 1 is elevation
                zr = src_ds.GetRasterBand(1).ComputeRasterMinMax() 
                this_region.zmin, this_region.zmax = zr[0], zr[1] 
                self.infos.minmax = this_region.export_as_list(include_z=True)
                self.infos.wkt = this_region.export_as_wkt()
                self.infos.numpts = ds_infos['nb']
                
        return(self.infos)

    def parse(self, resample=True):
        from .dlim import DatasetFactory
        
        parse_entry_p = True
        if self.region is not None:            
            if not regions.regions_intersect_p(
                    self.inf_region,
                    self.region \
                    if self.transform['trans_region'] is None \
                    else self.transform['trans_region']
            ):
                parse_entry_p = False

        if parse_entry_p:
            mt = gdal.Info(self.fn, format='json')['metadata']['']
            ds_infos = gdalfun.gdal_infos(self.fn)
            x_res = ds_infos['geoT'][1]
            sub_weight = max((3 * (10 if x_res <=3 else 1))/x_res, self.min_weight)
            oo = []
            if self.data_region is not None and self.data_region.valid_p():
                oo.append('MINX={}'.format(self.data_region.xmin))
                oo.append('MAXX={}'.format(self.data_region.xmax))
                oo.append('MINY={}'.format(self.data_region.ymin))
                oo.append('MAXY={}'.format(self.data_region.ymax))

            if ('HAS_SUPERGRIDS' in mt.keys() \
                and mt['HAS_SUPERGRIDS'] == 'TRUE') \
               or self.force_vr \
               or 'MAX_RESOLUTION_X' in mt.keys() \
               or 'MAX_RESOLUTION_Y' in mt.keys():
                if self.explode:
                    min_res = float(mt['MIN_RESOLUTION_X'])
                    max_res = float(mt['MAX_RESOLUTION_X'])
                    oo.append("MODE=LIST_SUPERGRIDS")
                    src_ds = gdal.OpenEx(self.fn, open_options=oo)
                    sub_datasets = src_ds.GetSubDatasets()
                    src_ds = None

                    utils.echo_msg(f'VRBAG min resolution is {min_res}')
                    utils.echo_msg(f'VRBAG max resolution is {max_res}')
                    with tqdm(
                            total=len(sub_datasets),
                            desc=f'parsing {len(sub_datasets)} supergrids from BAG file {self.fn}',
                            leave=True
                    ) as pbar:                
                        for sub_dataset in sub_datasets:
                            pbar.update()
                            res = sub_dataset[-1].split(',')[-2:]
                            #utils.echo_msg(res)
                            res = [float(re.findall(r'\d+\.\d+|\d+', x)[0]) for x in res]
                            sub_weight = max((3 * (10 if res[0] <=3 else 1))/res[0], self.min_weight)
                            sub_ds = DatasetFactory(
                                **self._set_params(
                                    mod=sub_dataset[0],
                                    data_format=200,
                                    src_srs=self.src_srs,
                                    band_no=1,
                                    uncertainty_mask_to_meter=0.01,
                                    check_path=False,
                                    super_grid=True,
                                    #node='pixel',
                                    #weight_multiplier=sub_weight,
                                    weight=sub_weight*(self.weight if self.weight is not None else 1),
                                    uncertainty_mask=2,
                                )
                            )._acquire_module()
                            self.data_entries.append(sub_ds)
                            sub_ds.initialize()
                            for gdal_ds in sub_ds.parse():
                                yield(gdal_ds)

                                utils.remove_glob(f'{gdal_ds.fn}.inf')


                elif self.vr_resampled_grid or self.vr_interpolate:
                    if self.vr_resampled_grid:
                        oo.append("MODE=RESAMPLED_GRID")
                    elif self.vr_interpolate:
                        oo.append("MODE=INTERPOLATED")

                    oo.append("RES_STRATEGY={}".format(self.vr_strategy))
                    sub_ds = DatasetFactory(
                        **self._set_params(
                            mod=self.fn,
                            data_format=200,
                            band_no=1,
                            open_options=oo,
                            weight=sub_weight*(self.weight if self.weight is not None else 1),
                            uncertainty_mask=2,
                            uncertainty_mask_to_meter=0.01,
                        )
                    )._acquire_module()
                    self.data_entries.append(sub_ds)
                    sub_ds.initialize()
                    for gdal_ds in sub_ds.parse():
                        yield(gdal_ds)
                else: # use vrbag.py
                    tmp_bag_as_tif = utils.make_temp_fn(
                        '{}_tmp.tif'.format(utils.fn_basename2(self.fn))
                    )
                    # scale cellsize to meters,
                    # todo: check if input is degress/meters/feet
                    #self.x_inc * 111120 
                    sr_cell_size = None
                    vrbag.interpolate_vr_bag(
                        self.fn, tmp_bag_as_tif,
                        self.cache_dir, sr_cell_size=sr_cell_size,
                        use_blocks=True,
                        method='linear',
                        nodata=3.4028234663852886e+38
                    )
                    sub_ds = DatasetFactory(
                        **self._set_params(
                            mod=tmp_bag_as_tif,
                            data_format=200,
                            band_no=1,
                            weight=sub_weight*(self.weight if self.weight is not None else 1),
                            uncertainty_mask=2,
                            uncertainty_mask_to_meter=0.01,
                        )
                    )._acquire_module()
                    self.data_entries.append(sub_ds)
                    sub_ds.initialize()
                    for gdal_ds in sub_ds.parse():
                        yield(gdal_ds)

            else:
                ds_infos = gdalfun.gdal_infos(self.fn)
                x_res = ds_infos['geoT'][1]
                sub_weight = max((3 * (10 if x_res <=3 else 1))/x_res, self.min_weight)
                sub_ds = DatasetFactory(
                    **self._set_params(
                        mod=self.fn,
                        data_format=200,
                        band_no=1,
                        uncertainty_mask=2,
                        uncertainty_mask_to_meter=0.01,
                        weight=sub_weight*(self.weight if self.weight is not None else 1),
                    )
                )._acquire_module()
                self.data_entries.append(sub_ds)
                sub_ds.initialize()
                yield(sub_ds)
                #for gdal_ds in sub_ds.parse():
                #    yield(gdal_ds)

                
### End
