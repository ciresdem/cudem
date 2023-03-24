### waffles.py
##
## Copyright (c) 2010 - 2023 Regents of the University of Colorado
##
## waffles.py is part of CUDEM
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
## Generate DEMs from a variety of data sources.
## Use CLI command 'waffles'
##
## Supported input datatypes include:
## datalist, las/laz, gdal, bag, xyz, mbs, fetches
##
## Supported gridding modules include:
## surface (GMT), triangulate (GMT), nearneighbor (GMT), mbgrid (MB-System), IDW (CUDEM), num (CUDEM/GMT),
## coastline (CUDEM), cudem (CUDEM), stacks (CUDEM), inv-dst (GDAL), linear (GDAL), average (GDAL),
## nearest (GDAL), scipy (SCIPY)
##
## GMT, GDAL and MB-System are required for full functionality.
##
### Code:

import sys
import os
import math
import json
import time
from tqdm import tqdm
import traceback

import numpy as np
from scipy import interpolate
from scipy import spatial
from scipy import ndimage
import threading        
from osgeo import gdal
from osgeo import ogr
from osgeo import osr

import cudem

from cudem import dlim
from cudem import datasets
from cudem import regions
from cudem import utils
from cudem import xyzfun
from cudem import demfun
from cudem import metadata
from cudem import vdatumfun

## Data cache directory
## mostly stores fetches related data here
waffles_cache = utils.cudem_cache()

## TODO
## add upper/lower limits to waffle class and
## remove from individual waffle modules.
class Waffle:
    """Representing a WAFFLES DEM/MODULE.
    Specific Gridding modules are sub-classes of this class.
    See WaffleFactory for module specifications and generation.

    Procedures:
      yield_xyz() - yield the xyz data from self.data
      dump_xyz() - dump the xyz data from self.data to port
      run() - run the WAFFLES module (function set via module sub-class)
      generate() - run and process the WAFFLES module
    """
    
    def __init__(
            self,
            data=[],
            src_region=None,
            inc=None,
            xinc=None,
            yinc=None,
            xsize=None,
            ysize=None,
            name='waffles_dem',
            node='pixel',
            fmt='GTiff',
            extend=0,
            extend_proc=0,
            weights=None,
            fltr=[],
            sample='bilinear',
            xsample=None,
            ysample=None,
            clip=None,
            chunk=None,
            dst_srs=None,
            srs_transform=False,
            verbose=False,
            archive=False,
            mask=False,
            spat=False,
            clobber=True,
            ndv=-9999,
            block=False,
            cache_dir=waffles_cache
    ):

        self.data = data
        self.datalist = None
        self.region = src_region
        self.inc = inc
        self.xinc = xinc
        self.yinc = yinc
        self.sample = sample
        self.xsample = xsample
        self.ysample = ysample
        self.name = name
        self.node = node
        self.fmt = fmt
        self.extend = extend
        self.extend_proc = extend_proc
        self.weights = weights
        self.fltr = fltr
        self.clip = clip
        self.chunk = chunk
        self.dst_srs = dst_srs
        self.srs_transform = srs_transform
        self.archive = archive
        self.mask = mask
        self.clobber = clobber
        self.verbose = verbose
        self.cache_dir = cache_dir
        self.gc = utils.config_check()
        self.spat = spat
        self.ndv = ndv
        self.block = block
        self.block_t = None
        self.ogr_ds = None
        self._init_regions()
        self.data_ = data
        self.fn = '{}.tif'.format(self.name)
        self.mask_fn = '{}_m.tif'.format(self.name)
        self.waffled = False
        self.aux_dems = []

        ## setting set_incs to True will force dlim to process the data to the set increments (raster mainly)
        self._init_data(set_incs=True)
        #self._init_data()
            
    def _init_regions(self):

        if self.node == 'grid':
            self.region = self.region.buffer(x_bv=self.xinc*.5, y_bv=self.yinc*.5)
            
        self.d_region = self._dist_region()
        self.p_region = self._proc_region()
        self.c_region = self._coast_region()
        self.ps_region = self.p_region.copy()
        self.ps_region = self.ps_region.buffer(
            x_bv=self.xinc*-.5, y_bv=self.yinc*-.5, x_inc=self.xinc, y_inc=self.yinc
        )

        if self.verbose:
            utils.echo_msg('target region: {}'.format(self.region))
            utils.echo_msg('processing region: {}'.format(self.p_region))
        
    def _init_data(self, set_incs=False):
        """Initialize the data for processing
        parses data paths to dlim dataset objects.

        set `set_incs` to True to block/sample datasets to given increment
        this function sets `self.data` to a list of dataset objects.
        """

        self.data = [dlim.DatasetFactory(
            fn=" ".join(['-' if x == "" else x for x in dl.split(",")]),
            src_region=self.p_region,
            verbose=self.verbose,
            dst_srs=self.dst_srs if self.srs_transform else None,
            weight=self.weights,
            x_inc=self.xinc if set_incs else None,
            y_inc=self.yinc if set_incs else None,
            sample_alg=self.sample,
            cache_dir=self.cache_dir
        ).acquire() for dl in self.data_]
        self.data = [d for d in self.data if d is not None]

    def _coast_region(self):
        """processing region (extended by self.extend_proc)"""

        cr = self.d_region.copy()
        return(cr.buffer(pct=self.extend_proc, x_inc=self.xinc, y_inc=self.yinc))
    
    def _proc_region(self):
        """processing region (extended by percentage self.extend_proc)"""

        pr = self.d_region.copy()
        return(pr.buffer(pct=self.extend_proc, x_inc=self.xinc, y_inc=self.yinc))
    
    def _dist_region(self):
        """distribution region (extended by self.extend)."""
        
        #dr = regions.Region().from_region(self.region)
        dr = self.region.copy()
        if self.xsample is None and self.ysample is None:
            return(
                dr.buffer(
                    x_bv=(self.xinc*self.extend),
                    y_bv=(self.yinc*self.extend)
                )
            )
        else:
            return(
                dr.buffer(
                    x_bv=(self.xsample*self.extend),
                    y_bv=(self.ysample*self.extend)
                )
            )
    def copy(self):
        """return a copy of the waffle object"""
        
        return(Waffle(
            data=self.data_,
            src_region=self.region,
            inc=self.inc,
            xinc=self.xinc,
            yinc=self.yinc,
            name=self.name,
            node=self.node,
            fmt=self.fmt,
            extend=self.extend,
            extend_proc=self.extend_proc,
            weights=self.weights,
            fltr=self.fltr,
            sample=self.sample,
            xsample=self.xsample,
            ysample=self.ysample,
            clip=self.clip,
            chunk=self.chunk,
            dst_srs=self.dst_srs,
            srs_transform=self.srs_transform,
            verbose=self.verbose,
            archive=self.archive,
            mask=self.mask,
            spat=self.spat,
            clobber=self.clobber,
            ndv=self.ndv,
            block=self.block
        ))

    ## todo: add blend option
    def _stacks_array(
            self,
            supercede=False,
            out_name=None,
            method='weighted_mean'
    ):
        """stack incoming arrays together

        method is either 'weighted_mean' or 'supercede'
        """
        
        if not self.weights:
            self.weights = 1

        xcount, ycount, dst_gt = self.p_region.geo_transform(
            x_inc=self.xinc, y_inc=self.yinc, node='grid'
        )

        gdt = gdal.GDT_Float32
        c_gdt = gdal.GDT_Int32
        driver = gdal.GetDriverByName(self.fmt)

        z_ds = driver.Create(
            '{}_s.tif'.format(out_name), xcount, ycount, 1, gdt,
            options=['COMPRESS=LZW', 'PREDICTOR=2', 'TILED=YES']
        )
        z_ds.SetGeoTransform(dst_gt)
        z_band = z_ds.GetRasterBand(1)
        z_band.SetNoDataValue(self.ndv)
        
        w_ds = driver.Create(
            '{}_w.tif'.format(out_name), xcount, ycount, 1, gdt,
            options=['COMPRESS=LZW', 'PREDICTOR=2', 'TILED=YES']
        )
        w_ds.SetGeoTransform(dst_gt)
        w_band = w_ds.GetRasterBand(1)
        w_band.SetNoDataValue(self.ndv)
                
        c_ds = driver.Create(
            '{}_c.tif'.format(out_name), xcount, ycount, 1, gdt
        )
        c_ds.SetGeoTransform(dst_gt)
        c_band = c_ds.GetRasterBand(1)
        c_band.SetNoDataValue(self.ndv)

        if self.verbose:
            utils.echo_msg('stacking data to {}/{} grid using {} method to {}'.format(
                ycount, xcount, 'supercede' if supercede else 'weighted mean', out_name
            ))

        ## incoming arrays can be quite large...perhaps chunks these
        for arrs, srcwin, gt in self.yield_array():
            arr = arrs['z']
            w_arr = arrs['weight']
            c_arr = arrs['count']
            c_arr[np.isnan(arr)] = 0
            w_arr[np.isnan(arr)] = 0
            arr[np.isnan(arr)] = 0
            z_array = z_band.ReadAsArray(srcwin[0], srcwin[1], srcwin[2], srcwin[3])
            w_array = w_band.ReadAsArray(srcwin[0], srcwin[1], srcwin[2], srcwin[3])
            c_array = c_band.ReadAsArray(srcwin[0], srcwin[1], srcwin[2], srcwin[3])

            z_array[z_array == self.ndv] = 0
            w_array[w_array == self.ndv] = 0
            c_array[c_array == self.ndv] = 0

            c_array += c_arr
            if supercede:
                z_array[w_arr > w_array] = arr[w_arr > w_array]
                w_array[w_arr > w_array] = w_arr[w_arr > w_array]
            else:
                z_array += (arr * w_arr)
                w_array += w_arr

            c_array[c_array == 0] = self.ndv
            w_array[w_array == 0] = self.ndv
            z_array[np.isnan(w_array)] = self.ndv
            z_array[w_array == self.ndv] = self.ndv
            
            # write out results
            z_band.WriteArray(z_array, srcwin[0], srcwin[1])
            w_band.WriteArray(w_array, srcwin[0], srcwin[1])
            c_band.WriteArray(c_array, srcwin[0], srcwin[1])
            arr = w_arr = c_arr = z_array = w_array = c_array = None

        ## Finalize and close datasets
        if not supercede:
            srcwin = (0, 0, z_ds.RasterXSize, z_ds.RasterYSize)
            for y in range(
                    srcwin[1], srcwin[1] + srcwin[3], 1
            ):
                z_data = z_band.ReadAsArray(
                    srcwin[0], y, srcwin[2], 1
                )
                c_data = c_band.ReadAsArray(
                    srcwin[0], y, srcwin[2], 1
                )
                w_data = w_band.ReadAsArray(
                    srcwin[0], y, srcwin[2], 1
                )

                z_data[z_data == self.ndv] = np.nan
                w_data[w_data == self.ndv] = np.nan
                c_data[c_data == self.ndv] = np.nan
                z_data[np.isnan(w_data)] = np.nan
                
                w_data = w_data / c_data
                z_data = (z_data/w_data)/c_data

                z_data[np.isnan(z_data)] = self.ndv
                c_data[np.isnan(c_data)] = self.ndv
                w_data[np.isnan(w_data)] = self.ndv
                
                z_band.WriteArray(z_data, srcwin[0], y)
                c_band.WriteArray(c_data, srcwin[0], y)
                w_band.WriteArray(w_data, srcwin[0], y)

        #utils.echo_msg('maximum stacked value: {}'.format(z_data.max()))
        #utils.echo_msg('minimum stacked value: {}'.format(z_data.min()))
        z_ds = c_ds = w_ds = None
    
    def yield_array(self, **kwargs):
        """yield the arrays from the datalist for use in gridding"""
        
        if self.mask:
            xcount, ycount, dst_gt = self.region.geo_transform(x_inc=self.xinc, y_inc=self.yinc)
            mask_array = np.zeros((ycount, xcount))
            ds_config = demfun.set_infos(
                xcount, ycount, (xcount*ycount), dst_gt, utils.sr_wkt(self.dst_srs),
                gdal.GDT_Float32, self.ndv, 'GTiff'
        )

        with tqdm(desc='{}: parsing ARRAY data'.format(utils._command_name)) as pbar:
            for xdl in self.data:
                for array in xdl.yield_array():
                    if self.mask:
                        cnt_arr = array[0]['count']
                        srcwin = array[1]
                        gt = array[2]
                        cnt_arr[np.isnan(cnt_arr)] = 0
                        mask_array[srcwin[1]:srcwin[1]+srcwin[3],srcwin[0]:srcwin[0]+srcwin[2]] += cnt_arr
                        mask_array[mask_array > 0] = 1

                    yield(array)
                    pbar.update()
                
        if self.mask:
            utils.gdal_write(mask_array, self.mask_fn, ds_config, verbose=self.verbose)
            mask_array = None
                    
    def yield_xyz(self, region=None, **kwargs):
        """yields the xyz data"""

        if self.mask:
            xcount, ycount, dst_gt = self.region.geo_transform(x_inc=self.xinc, y_inc=self.yinc)
            msk_array = np.zeros((ycount, xcount))
            ds_config = demfun.set_infos(
                xcount, ycount, (xcount*ycount), dst_gt, utils.sr_wkt(self.dst_srs),
                gdal.GDT_Float32, self.ndv, 'GTiff'
        )

        with tqdm(desc='{}: parsing XYZ data'.format(utils._command_name)) as pbar:
            for xdl in self.data:
                if self.spat:
                    xyz_yield = metadata.SpatialMetadata(
                        data=[xdl.fn],
                        src_region=self.p_region,
                        inc=self.xinc,
                        extend=self.extend,
                        dst_srs=self.dst_srs if self.srs_transform else None,
                        node=self.node,
                        name=self.name,
                        verbose=self.verbose,
                        make_valid=True
                    ).yield_xyz()

                elif self.archive:
                    xyz_yield = xdl.archive_xyz()
                else:
                    xyz_yield = xdl.yield_xyz()

                if self.block:
                    #xyz_yield = self._xyz_block(xyz_yield, out_name=self.block) if utils.str_or(self.block) != 'False' else self._xyz_block(xyz_yield)
                    xyz_yield = xdl.block_xyz()

                for xyz in xyz_yield:
                    yield(xyz)
                    pbar.update()
                    if self.mask:
                        if regions.xyz_in_region_p(xyz, self.region):
                            xpos, ypos = utils._geo2pixel(
                                xyz.x, xyz.y, dst_gt, 'pixel'
                            )
                            try:
                                msk_array[ypos, xpos] = 1
                            except:
                                pass
        if self.mask:
            out, status = utils.gdal_write(msk_array, self.mask_fn, ds_config)    

        if self.spat:
            dst_layer = '{}_sm'.format(self.name)
            dst_vector = dst_layer + '.shp'
            utils.run_cmd(
                'ogr2ogr -clipsrc {} __tmp_clip.shp {} -overwrite -nlt POLYGON -skipfailures'.format(
                    self.d_region.format('ul_lr'), dst_vector
                ),
                verbose=True
            )
            utils.run_cmd(
                'ogr2ogr {} __tmp_clip.shp -overwrite'.format(dst_vector),
                verbose=True
            )
            utils.remove_glob('__tmp_clip.*')

    def dump_xyz(self, dst_port=sys.stdout, encode=False, **kwargs):
        """dump the xyz data to dst_port"""

        #with tqdm(desc='dumping xyz data...') as pbar:
        for xyz in self.yield_xyz(**kwargs):
            #pbar.update()
            xyz.dump(
                include_w = True if self.weights is not None else False,
                dst_port=dst_port,
                encode=encode,
                **kwargs
            )

    def _process(self, fn=None, filter_=False):
        """process the outpt WAFFLES DEM (filter, clip, cut, set)."""
        
        if fn is None:
            fn = self.fn

        ## SET NODATA
        demfun.set_nodata(fn, nodata=self.ndv, convert_array=True, verbose=self.verbose)

        ## FILTER
        if filter_:
            if len(self.fltr) > 0:
                for f in self.fltr:
                    fltr_val = None
                    split_val = None
                    fltr_opts = f.split(':')
                    fltr = fltr_opts[0]
                    if len(fltr_opts) > 1:
                        fltr_val = fltr_opts[1]
                        
                    if len(fltr_opts) > 2:
                        split_val= fltr_opts[2]

                    # fails if fltr_val in float
                    if demfun.filter_(
                            fn, '__tmp_fltr.tif', fltr=fltr, fltr_val=fltr_val, split_val=split_val,
                    ) == 0:
                        os.rename('__tmp_fltr.tif', fn)

        ## SAMPLE
        if self.xsample is not None or self.ysample is not None:
            if demfun.sample_warp(fn, '__tmp_sample.tif', self.xsample, self.ysample,
                             src_region=self.p_region, sample_alg=self.sample,
                             ndv=self.ndv, verbose=self.verbose)[1] == 0:
                os.rename('__tmp_sample.tif', fn)


        ## CLIP
        if self.clip is not None:
            clip_args = {}
            cp = self.clip.split(':')
            clip_args['src_ply'] = cp[0]
            clip_args['verbose'] = self.verbose
            clip_args = utils.args2dict(cp[1:], clip_args)
            if clip_args['src_ply'] == 'coastline':
                self.coast = WaffleFactory(
                    mod='coastline:polygonize=False',
                    data=self.data_,
                    src_region=self.p_region,
                    xinc=self.xsample if self.xsample is not None else self.xinc,
                    yinc=self.ysample if self.ysample is not None else self.yinc,
                    name='tmp_coast',
                    node=self.node,
                    weights=self.weights,
                    dst_srs=self.dst_srs,
                    srs_transform=self.srs_transform,
                    clobber=True,
                    verbose=self.verbose,
                ).acquire().generate()
                
                demfun.mask_(fn, self.coast.fn, '__tmp_clip__.tif', msk_value=1)
                os.rename('__tmp_clip__.tif', '{}'.format(fn))
                
            else:
                ## if clip vector is empty, will return surface with all upper_value :(
                ## maybe have it remove all data in such cases (if invert is true)
                if demfun.clip(fn, '__tmp_clip__.tif', **clip_args)[1] == 0:
                    os.rename('__tmp_clip__.tif', '{}'.format(fn))

        ## CUT
        #if demfun.cut(fn, self.d_region, '__tmp_cut__.tif', node='grid', mode=None)[1] == 0:
        _tmp_cut, cut_status = demfun.cut(fn, self.d_region, '__tmp_cut__.tif', node='grid', mode=None)
        if cut_status == 0:
            try:
                os.rename(_tmp_cut, fn)
            except Exception as e:
                utils.echo_error_msg('could not cut {}; {}'.format(fn, e))
        else:
            utils.echo_error_msg('could not cut {}'.format(fn))

        ## SET SRS/METADATA
        ## if set_srs fails, set_metadata will skip first entry...
        demfun.set_srs(fn, self.dst_srs, verbose=self.verbose)
        demfun.set_metadata(fn, node=self.node, cudem=True, verbose=self.verbose)

        ## REFORMAT
        if self.fmt != 'GTiff':
            out_dem = utils.gdal2gdal(fn, dst_fmt=self.fmt)
            if out_dem is not None:
                utils.remove_glob(fn)
                
        return(self)
    
    def generate(self):
        """run and process the WAFFLES module"""

        if os.path.exists(self.fn):
            if not self.clobber:
                utils.echo_warning_msg(
                    'DEM {} already exists, skipping...'.format(self.fn)
                )
                if self.mask:
                    if not os.path.exists(self.mask_fn):
                        [x for x in self.yield_xyz()]
                        self._process(fn=self.mask_fn, filter_=False)
                return(self)
        else:
            if not os.path.exists(os.path.dirname(self.fn)):
                try:
                    os.makedirs(os.path.dirname(self.fn))
                except: pass
            
        ## Generate in Chunks of self.chunk by self.chunk
        if self.chunk is not None:
            xcount, ycount, dst_gt = self.p_region.geo_transform(x_inc=self.xinc, y_inc=self.yinc, node='grid')
            count = 0
            chunks = []
            for srcwin in utils.yield_srcwin((ycount, xcount), self.chunk):
                count += 1
                #print(srcwin)
                this_geo_x_origin, this_geo_y_origin = utils._pixel2geo(srcwin[0], srcwin[1], dst_gt)
                this_geo_x_end, this_geo_y_end = utils._pixel2geo(srcwin[0]+srcwin[2], srcwin[1]+srcwin[3], dst_gt)
                this_gt = [this_geo_x_origin, float(dst_gt[1]), 0.0, this_geo_y_origin, 0.0, float(dst_gt[5])]
                #this_region = self.region.copy()
                this_region = regions.Region()
                this_region.from_geo_transform(geo_transform=this_gt, x_count=srcwin[2], y_count=srcwin[3])
                this_region.buffer(pct=10, x_inc = self.xinc, y_inc = self.yinc)
                this_waffle = WaffleFactory()._modules[self.mod]['class'](
                    data=self.data_,
                    src_region=this_region,
                    inc=self.inc,
                    xinc=self.xinc,
                    yinc=self.yinc,
                    name='{}_{}'.format(self.name, count),
                    node=self.node,
                    fmt=self.fmt,
                    extend=self.extend,
                    extend_proc=self.extend_proc+1,
                    weights=self.weights,
                    fltr=self.fltr,
                    sample=self.sample,
                    xsample=self.xsample,
                    ysample=self.ysample,
                    clip=self.clip,
                    chunk=None,
                    dst_srs=self.dst_srs,
                    srs_transform=self.srs_transform,
                    verbose=self.verbose,
                    archive=self.archive,
                    mask=self.mask,
                    spat=self.spat,
                    clobber=self.clobber,
                    block=self.block,
                    ndv=self.ndv,
                    cache_dir=self.cache_dir,
                    **self.mod_args
                )
                
                this_waffle.run()
                if this_waffle.valid_p():
                    this_waffle._process(filter_=True)
                    chunks.append(this_waffle.fn)
                
            if len(chunks) > 0:
                ## todo: use demfun.sample_warp
                g = gdal.Warp(self.fn, chunks, format='GTiff', resampleAlg='cubicspline',
                              options=["COMPRESS=LZW", "TILED=YES"])
                g = None
                
            utils.remove_glob(*chunks)
            # if self.valid_p():
            #     return(self._process(filter_=True))
            # else:
            #     return(self)
        else:
            self.run()
            if self.mask:
                if os.path.exists(self.mask_fn):
                    self._process(fn=self.mask_fn, filter_=False)

        if self.valid_p():
            self.waffled = True
            [self._process(fn=x, filter_=False) for x in self.aux_dems]
            return(self._process(filter_=True))
        else:
            return(self)

    def set_limits(self, upper_limit=None, lower_limit=None):
        upper_limit = utils.float_or(upper_limit)
        lower_limit = utils.float_or(lower_limit)
        
        if upper_limit is not None or lower_limit is not None:
            src_ds = gdal.Open(self.fn, 1)
            if src_ds is not None:
                src_band = src_ds.GetRasterBand(1)
                band_data = src_band.ReadAsArray()

                if upper_limit is not None:
                    utils.echo_msg('setting upper_limit to {}'.format(upper_limit))
                    band_data[band_data > upper_limit] = upper_limit

                if lower_limit is not None:
                    utils.echo_msg('setting lower_limit to {}'.format(lower_limit))
                    band_data[band_data < lower_limit] = lower_limit

                src_band.WriteArray(band_data)
                src_ds = None
                
        return(self)
            
    def valid_p(self):
        """check if the output WAFFLES DEM is valid"""

        if not os.path.exists(self.fn):
            return(False)
        
        if self.mask:
            if not os.path.exists(self.mask_fn):
                return(False)
            
        gdi = demfun.infos(self.fn, scan=True)        
        if gdi is not None:
            if np.isnan(gdi['zr'][0]):
                return(False)
        else:
            return(False)
        
        return(True)            

    def run(self):
        """run the WAFFLES module (set via sub-module class)."""
        
        raise(NotImplementedError)
    
## ==============================================
## GMT Surface
##
## TODO: update to use pygmt
## ==============================================
class GMTSurface(Waffle):
    """SPLINE DEM via GMT surface
    
Generate a DEM using GMT's surface command
see gmt surface --help for more info.

    ---
    Parameters:
    
    tension=[0-1] - spline tension.
    relaxation=[val] - spline relaxation factor.
    lower_limit=[val] - constrain interpolation to lower limit.
    upper_limit=[val] - constrain interpolation to upper limit.
    """
    
    def __init__(self, tension=.35, relaxation=1, max_radius=None,
                 lower_limit=None, upper_limit=None, aspect=None,
                 breakline=None, convergence=None, blockmean=True,
                 geographic=True, **kwargs):
        """generate a DEM with GMT surface"""

        self.mod = 'surface'
        self.mod_args = {
            'tension':tension,
            'relaxation':relaxation,
            'aspect':aspect,
            'max_radius':max_radius,
            'lower_limit':lower_limit,
            'upper_limit':upper_limit,
            'breakline':breakline,
            'convergence':convergence,
            'blockmean':blockmean,
            'geographic':geographic
        }
        super().__init__(**kwargs)
        if self.gc['GMT'] is None:
            utils.echo_error_msg(
                'GMT must be installed to use the SURFACE module'
            )
            return(None, -1)
        
        self.tension = tension
        self.convergence = utils.float_or(convergence)
        self.relaxation = relaxation
        #self.lower_limit = utils.float_or(lower_limit, 'd')
        #self.upper_limit = utils.float_or(upper_limit, 'd')
        self.upper_limit = upper_limit
        self.lower_limit = lower_limit
        self.breakline = breakline
        self.max_radius = max_radius
        self.aspect = aspect
        self.blockmean = blockmean
        self.geographic = geographic
        out, status = utils.run_cmd(
            'gmt gmtset IO_COL_SEPARATOR = SPACE',
            verbose = False
        )        
        
    def run(self):

        xsize, ysize, gt = self.p_region.geo_transform(x_inc=self.xinc, node='grid')
        #print(self.ps_region.xmax-self.ps_region.xmin)
        #print(self.xinc * xsize)
        
        dem_surf_cmd = ('')
        if self.blockmean:
            dem_surf_cmd = (
                'gmt blockmean {} -I{:.16f}/{:.16f}+e{}{} -V |'.format(
                    self.ps_region.format('gmt'),
                    self.xinc,
                    self.yinc,
                    ' -W' if self.weights else '',
                    ' -fg' if self.geographic else '',
                )
            )

        dem_surf_cmd += (
            'gmt surface -V {} -I{:.16f}/{:.16f}+e -G{}.tif=gd+n{}:GTiff -T{} -Z{} {}{}{}{}{}'.format(
                self.ps_region.format('gmt'),
                self.xinc,
                self.yinc,
                self.name,
                self.ndv,
                self.tension,
                self.relaxation,
                ' -D{}'.format(self.breakline) if self.breakline is not None else '',
                ' -M{}'.format(self.max_radius) if self.max_radius is not None else '',
                ' -C{}'.format(self.convergence) if self.convergence is not None else '',
                ' -A{}'.format(self.aspect) if self.aspect is not None else '',
                ' -fg' if self.geographic else '',
            )
        )
        
        # dem_surf_cmd = (
        #     'gmt blockmean {} -I{:.14f}/{:.14f}{} -V | gmt surface -V {} -I{:.14f}/{:.14f} -G{}.tif=gd+n{}:GTiff -T{} -Z{} {}{}{}{}'.format(
        #         self.ps_region.format('gmt'),
        #         self.xinc,
        #         self.yinc,
        #         ' -W' if self.weights else '',
        #         self.ps_region.format('gmt'),
        #         self.xinc,
        #         self.yinc,
        #         self.name,
        #         self.ndv,
        #         self.tension,
        #         self.relaxation,
        #         ' -D{}'.format(self.breakline) if self.breakline is not None else '',
        #         ' -M{}'.format(self.max_radius) if self.max_radius is not None else '',
        #         ' -C{}'.format(self.convergence) if self.convergence is not None else '',
        #         ' -A{}'.format(self.aspect) if self.aspect is not None else ''
        #     )
        # )

        ## -L* messes up the output solution! use self.set_limits()
        #' -Ll{}'.format(self.lower_limit) if self.lower_limit is not None else '',
        #' -Lu{}'.format(self.upper_limit) if self.upper_limit is not None else '',
        
        out, status = utils.run_cmd(
            dem_surf_cmd,
            verbose=self.verbose,
            data_fun=lambda p: self.dump_xyz(
                dst_port=p, encode=True
            )
        )

        return(self.set_limits(self.upper_limit, self.lower_limit))
               #return(self)

## ==============================================
## GMT Triangulate
##
## TODO: update to use pygmt
## ==============================================
class GMTTriangulate(Waffle):
    """TRIANGULATION DEM via GMT triangulate
    
Generate a DEM using GMT's triangulate command.
see gmt triangulate --help for more info.        
    """
    
    def __init__(self, **kwargs):
        """generate a DEM with GMT triangulate"""

        self.mod = 'triangulate'
        self.mod_args = {}
        super().__init__(**kwargs)        
        if self.gc['GMT'] is None:
            utils.echo_error_msg(
                'GMT must be installed to use the TRIANGULATE module'
            )
            return(None, -1)

        out, status = utils.run_cmd(
            'gmt gmtset IO_COL_SEPARATOR = SPACE',
            verbose = False
        )        
        
    def run(self):
        dem_tri_cmd = 'gmt blockmean {} -I{:.14f}/{:.14f}{} -V | gmt triangulate -V {} -I{:.14f}/{:.14f} -G{}.tif=gd:GTiff'.format(
            self.ps_region.format('gmt'),
            self.xinc,
            self.yinc,
            ' -W' if self.weights else '',
            self.ps_region.format('gmt'),
            self.xinc,
            self.yinc,
            self.name
        )
        out, status = utils.run_cmd(
            dem_tri_cmd,
            verbose = self.verbose,
            data_fun = lambda p: self.dump_xyz(
                dst_port=p, encode=True
            )
        )        
        return(self)

## ==============================================
## GMT Near Neighbor
##
## TODO: update to use pygmt
## ==============================================
class GMTNearNeighbor(Waffle):
    """NEARNEIGHBOR DEM via GMT nearneighbor
    
Generate a DEM using GMT's nearneighbor command.
see gmt nearneighbor --help for more info.

    ---
    Parameters:
    
    radius=[val]\t\tsearch radius
    sectors=[val]\t\tsector information
    """
    
    def __init__(self, radius=None, sectors=None, **kwargs):
        """generate a DEM with GMT nearneighbor"""

        self.mod = 'nearneighbor'
        self.mod_args = {
            'radius':radius,
            'sectors':sectors
        }
        super().__init__(**kwargs) 
        if self.gc['GMT'] is None:
            utils.echo_error_msg(
                'GMT must be installed to use the NEARNEIGHBOR module'
            )
            return(None, -1)

        out, status = utils.run_cmd(
            'gmt gmtset IO_COL_SEPARATOR = SPACE',
            verbose = False
        )        
        self.radius = radius
        self.sectors = sectors
        
    def run(self):
        dem_nn_cmd = 'gmt blockmean {} -I{:.14f}/{:.14f}{} -V | gmt nearneighbor -V {} -I{:.14f}/{:.14f} -G{}.tif=gd+n{}:GTiff{}{}{}'.format(
            self.ps_region.format('gmt'),
            self.xinc,
            self.yinc,
            ' -W' if self.weights else '',
            self.ps_region.format('gmt'),
            self.xinc,
            self.yinc,
            self.name,
            self.ndv,
            ' -W' if self.weights else '',
            ' -N{}'.format(self.sectors) if self.sectors is not None else '',
            ' -S{}'.format(self.radius) if self.radius is not None else ' -S{}'.format(self.xinc),
        )
        out, status = utils.run_cmd(
            dem_nn_cmd,
            verbose = self.verbose,
            data_fun = lambda p: self.dump_xyz(
                dst_port=p, encode=True
            )
        )        
        return(self)

## ==============================================
## MB-System mbgrid
## ==============================================
class WafflesMBGrid(Waffle):
    """SPLINE DEM via MB-System's mbgrid.
    
Generate a DEM using MB-System's mbgrid command.
By default, will use MB-Systems datalist processes.
set `use_datalist=True` to use CUDEM's dlim instead.
see mbgrid --help for more info

< mbgrid:dist='10/3':tension=35:use_datalists=False >

    ---
    Parameters:
    
    dist=[val] - the dist variable to use in mbgrid
    tension=[val] - the spline tension value (0-inf)
    use_datalist=[True/False] - use built-in datalists rather than mbdatalist
    """
    
    def __init__(self, dist='10/3', tension=35, use_datalists=False, nc=False, **kwargs):
        self.mod = 'mbgrid'
        self.mod_args = {
            'dist':dist,
            'tension':tension,
            'use_datalists':use_datalists,
            'nc':nc
        }
        super().__init__(**kwargs)
        if self.gc['MBGRID'] is None:
            utils.echo_error_msg(
                'MB-System must be installed to use the MBGRID module'
            )
            return(None, -1)

        if self.gc['GMT'] is None:
            utils.echo_error_msg(
                'GMT must be installed to use the MBGRID module'
            )
            return(None, -1)

        self.nc = nc
        self.dist = dist
        self.tension = tension
        self.use_datalists = use_datalists
        self.mod = 'mbgrid'
                
    def _gmt_num_msk(self, num_grd, dst_msk):
        """generate a num-msk from a NUM grid using GMT grdmath

        Args:
          num_grd (str): pathname to a source `num` grid file
          dst_msk (str): pathname to a destination `msk` grid file

        Returns:
          list: [cmd-output, cmd-return-code]
        """

        num_msk_cmd = 'gmt grdmath -V {} 0 MUL 1 ADD 0 AND = {}'.format(
            num_grd, dst_msk
        )
        return(
            utils.run_cmd(
                num_msk_cmd, verbose=self.verbose
            )
        )

    def _gmt_grd2gdal(self, src_grd, dst_fmt='GTiff'):
        """convert the grd file to tif using GMT

        Args:
          src_grd (str): a pathname to a grid file
          dst_fmt (str): the output GDAL format string

        Returns:
          str: the gdal file name or None
        """

        dst_gdal = '{}.{}'.format(
            os.path.basename(src_grd).split('.')[0], utils.gdal_fext(dst_fmt)
        )        
        grd2gdal_cmd = 'gmt grdconvert {} {}=gd+n{}:{} -V'.format(
            src_grd, dst_gdal, self.ndv, dst_fmt
        )
        out, status = utils.run_cmd(
            grd2gdal_cmd, verbose=self.verbose
        )
        if status == 0:
            return(dst_gdal)
        else:
            return(None)
        
    def _gmt_grdsample(self, src_grd, dst_fmt='GTiff'):
        """convert the grd file to tif using GMT

        Args:
          src_grd (str): a pathname to a grid file
          dst_fmt (str): the output GDAL format string

        Returns:
          str: the gdal file name or None
        """

        dst_gdal = '{}.{}'.format(
            os.path.basename(src_grd).split('.')[0], utils.gdal_fext(dst_fmt)
        )
        grdsample_cmd = 'gmt grdsample {} -T -G{}=gd+n{}:{} -V'.format(
            src_grd, dst_gdal, self.ndv, dst_fmt
        )        
        out, status = utils.run_cmd(
            grdsample_cmd, verbose=self.verbose
        )        
        if status == 0:
            return(dst_gdal)        
        else:
            return(None)
    
    def run(self):
        if self.use_datalists:
            archive = self.archive
            self.archive = True
            [x for x in self.yield_xyz()]
            self.archive = archive

        mb_region = self.p_region.copy()
        mb_region = mb_region.buffer(x_bv=self.xinc*-.5, y_bv=self.yinc*-.5)
        xsize, ysize, gt = self.p_region.geo_transform(x_inc=self.xinc, node='grid')
        print(self.data[0].archive_datalist)
        if self.data[0].archive_datalist is not None:
            dl_name = self.data[0].archive_datalist
        else:
            dl_name = self.data[0].fn
            
        mbgrid_cmd = 'mbgrid -I{} {} -D{}/{} -O{} -A2 -F1 -N -C{} -S0 -X0.1 -T{} {}'.format(
            dl_name,
            mb_region.format('gmt'),
            xsize,
            ysize,
            self.name,
            self.dist,
            self.tension,
            '-M' if self.mask else ''
        )
        for out in utils.yield_cmd(mbgrid_cmd, verbose=self.verbose):
            sys.stderr.write('{}'.format(out))
            
        utils.gdal2gdal('{}.grd'.format(self.name))
        utils.remove_glob('*.cmd', '*.mb-1', '{}.grd'.format(self.name))
        if self.use_datalists and not self.archive:
            utils.remove_glob('archive')

        if self.mask:
            num_grd = '{}_num.grd'.format(self.name)
            dst_msk = '{}_m.tif=gd+n{}:GTiff'.format(self.name, self.ndv)
            self.mask_fn = dst_msk
            out, status = self._gmt_num_msk(
                num_grd, dst_msk, verbose=self.verbose
            )
            utils.remove_glob(num_grd, '*_sd.grd')
            if not self.use_datalists:
                if self.spat or self.archive:
                    [x for x in self.yield_xyz()]
                    
        return(self)

## ==============================================
## Waffles 'num' - no interpolation
##
## just use stacks...
## ==============================================
class WafflesNum(Waffle):
    """Uninterpolated DEM populated by <mode>.
    
Generate an uninterpolated DEM using <mode> option.
Using mode of 'A<mode>' uses GMT's xyz2grd command, 
see gmt xyz2grd --help for more info.

mode keys: k (mask), m (mean), n (num), w (wet), A<mode> (gmt xyz2grd)

< num:mode=n >

    ---
    Parameters:
    
    mode=[key] - specify mode of grid population
    """
    
    def __init__(self, mode='n', min_count=None, **kwargs):
        """generate an uninterpolated Grid
        `mode` of `n` generates a num grid
        `mode` of `m` generates a mean grid
        `mode` of `k` generates a mask grid
        `mode` of `w` generates a wet/dry mask grid
        `mode` of `A` generates via GMT 'xyz2grd'
        """
        
        self.mod = 'num'
        self.mod_args = {
            'mode':mode,
            'min_count':min_count
        }
        super().__init__(**kwargs)        
        self.mode = mode
        self.min_count = min_count
        
    def _xyz_num(self):
        """Create a GDAL supported grid from xyz data """

        xcount, ycount, dst_gt = self.p_region.geo_transform(x_inc=self.xinc, y_inc=self.yinc)
        if self.verbose:
            utils.echo_msg(
                'generating uninterpolated NUM grid `{}` @ {}/{}'.format(
                    self.mode, ycount, xcount
                )
            )          
            
        if self.mode == 'm' or self.mode == 'w':
            sum_array = np.zeros((ycount, xcount))
            
        count_array = np.zeros((ycount, xcount))
        gdt = gdal.GDT_Float32
        ds_config = demfun.set_infos(
            xcount,
            ycount,
            xcount * ycount,
            dst_gt,
            utils.sr_wkt(self.dst_srs),
            gdt,
            self.ndv,
            'GTiff'
        )
        
        for this_xyz in self.yield_xyz():
            if regions.xyz_in_region_p(this_xyz, self.p_region):
                xpos, ypos = utils._geo2pixel(
                    this_xyz.x, this_xyz.y, dst_gt, 'pixel'
                )
                
                try:
                    if self.mode == 'm' or self.mode == 'w':
                        sum_array[ypos, xpos] += this_xyz.z
                    if self.mode == 'n' or self.mode == 'm' or self.mode == 'w':
                        count_array[ypos, xpos] += 1
                    else:
                        count_array[ypos, xpos] = 1
                        
                except Exception as e:
                    pass
                
        if self.mode == 'm' or self.mode == 'w':
            count_array[count_array == 0] = np.nan
            out_array = (sum_array/count_array)
            if self.mode == 'w':
                out_array[out_array >= 0] = 1
                out_array[out_array < 0] = 0
        else:
            out_array = count_array

        out_array[np.isnan(out_array)] = self.ndv            
        utils.gdal_write(out_array, self.fn, ds_config)
        
    def run(self):
        if self.mode.startswith('A'):
            if self.gc['GMT'] is None:
                utils.echo_error_msg(
                    'GMT must be installed to use the Mode `A` with the NUM module'
                )
                return(None, -1)

            out, status = utils.run_cmd('gmt gmtset IO_COL_SEPARATOR = SPACE', verbose = False)            
            dem_xyz2grd_cmd = 'gmt xyz2grd -{} -V {} -I{:.14f}/{:.14f}+e -G{}.tif=gd:GTiff'.format(
                self.mode,
                self.ps_region.format('gmt'),
                self.xinc,
                self.yinc,
                self.name
            )
            
            out, status = utils.run_cmd(
                dem_xyz2grd_cmd,
                verbose=self.verbose,
                data_fun=lambda p: self.dump_xyz(
                    dst_port=p, encode=True
                )
            )
        else:
            self._xyz_num()
                        
        return(self)

## ==============================================
## Waffles IDW
## ==============================================
class Invdisttree():
    """ inverse-distance-weighted interpolation using KDTree:
    @Denis via https://stackoverflow.com/questions/3104781/inverse-distance-weighted-idw-interpolation-with-python
    https://creativecommons.org/licenses/by-nc-sa/3.0/

invdisttree = Invdisttree( X, z )  -- data points, values
interpol = invdisttree( q, nnear=3, eps=0, p=1, weights=None, stat=0 )
    interpolates z from the 3 points nearest each query point q;
    For example, interpol[ a query point q ]
    finds the 3 data points nearest q, at distances d1 d2 d3
    and returns the IDW average of the values z1 z2 z3
        (z1/d1 + z2/d2 + z3/d3)
        / (1/d1 + 1/d2 + 1/d3)
        = .55 z1 + .27 z2 + .18 z3  for distances 1 2 3

    q may be one point, or a batch of points.
    eps: approximate nearest, dist <= (1 + eps) * true nearest
    p: use 1 / distance**p
    weights: optional multipliers for 1 / distance**p, of the same shape as q
    stat: accumulate wsum, wn for average weights

How many nearest neighbors should one take ?
a) start with 8 11 14 .. 28 in 2d 3d 4d .. 10d; see Wendel's formula
b) make 3 runs with nnear= e.g. 6 8 10, and look at the results --
    |interpol 6 - interpol 8| etc., or |f - interpol*| if you have f(q).
    I find that runtimes don't increase much at all with nnear -- ymmv.

p=1, p=2 ?
    p=2 weights nearer points more, farther points less.
    In 2d, the circles around query points have areas ~ distance**2,
    so p=2 is inverse-area weighting. For example,
        (z1/area1 + z2/area2 + z3/area3)
        / (1/area1 + 1/area2 + 1/area3)
        = .74 z1 + .18 z2 + .08 z3  for distances 1 2 3
    Similarly, in 3d, p=3 is inverse-volume weighting.

Scaling:
    if different X coordinates measure different things, Euclidean distance
    can be way off.  For example, if X0 is in the range 0 to 1
    but X1 0 to 1000, the X1 distances will swamp X0;
    rescale the data, i.e. make X0.std() ~= X1.std() .

A nice property of IDW is that it's scale-free around query points:
if I have values z1 z2 z3 from 3 points at distances d1 d2 d3,
the IDW average
    (z1/d1 + z2/d2 + z3/d3)
    / (1/d1 + 1/d2 + 1/d3)
is the same for distances 1 2 3, or 10 20 30 -- only the ratios matter.
In contrast, the commonly-used Gaussian kernel exp( - (distance/h)**2 )
is exceedingly sensitive to distance and to h.

quite heavy on memory when large grid-size...

    """
# anykernel( dj / av dj ) is also scale-free
# error analysis, |f(x) - idw(x)| ? todo: regular grid, nnear ndim+1, 2*ndim

    def __init__( self, X, z, leafsize=10, stat=0 ):
        assert len(X) == len(z), "len(X) %d != len(z) %d" % (len(X), len(z))
        self.tree = spatial.cKDTree( X, leafsize=leafsize )  # build the tree
        self.z = z
        self.stat = stat
        self.wn = 0
        self.wsum = None;

    def __call__(self, q, nnear=6, eps=0, p=1, dub=np.inf, weights=None):
        # nnear nearest neighbours of each query point --
        q = np.asarray(q)
        qdim = q.ndim
        if qdim == 1:
            q = np.array([q])
        if self.wsum is None:
            self.wsum = np.zeros(nnear)

        eps = utils.float_or(eps, .1)
        self.distances, self.ix = self.tree.query(q, k=nnear, eps=eps, distance_upper_bound=dub)
        interpol = np.zeros((len(self.distances),) + np.shape(self.z[0]))
        jinterpol = 0
        for dist, ix in zip( self.distances, self.ix ):
            if np.any(np.isinf(dist)):
                wz = np.nan
            elif nnear == 1:
                wz = self.z[ix]
            elif dist[0] < 1e-10:
                wz = self.z[ix[0]]
            else:  # weight z s by 1/dist --
                w = 1 / dist**p
                if weights is not None:
                    w *= weights[ix]  # >= 0
                    
                w /= np.sum(w)
                wz = np.dot( w, self.z[ix] )
                if self.stat:
                    self.wn += 1
                    self.wsum += w
                    
            interpol[jinterpol] = wz
            jinterpol += 1
        return interpol if qdim > 1  else interpol[0]

class WafflesIDW(Waffle):
    """INVERSE DISTANCE WEIGHTED DEM
    
Generate a DEM using an Inverse Distance Weighted algorithm.
If weights are used, will generate a UIDW DEM, using weight values as inverse uncertainty,
as described here: https://ir.library.oregonstate.edu/concern/graduate_projects/79407x932
and here: https://stackoverflow.com/questions/3104781/inverse-distance-weighted-idw-interpolation-with-python

< IDW:min_points=8:radius=inf:power=1:upper_limit=None:lower_limit=None:supercede=False:keep_auxiliary=False:chunk_size=None >

    ---
    Parameters:
    
    power=[val] - weight**power
    min_points=[val] - minimum neighbor points for IDW
    radius=[val] - search radius (in cells), only fill data cells within radius from data
    upper_limit=[val] - Restrict output DEM to values below val
    lower_limit=[val] - Restrict output DEM to values above val
    supercede=[True/False] - supercede higher weighted data,
    keep_auxiliary=[True/False] - retain auxiliary files
    chunk_size=[val] - size of chunks in pixels
    """
    
    def __init__(
            self,
            power=1,
            min_points=8,
            upper_limit=None,
            lower_limit=None,
            radius=None,
            supercede=False,
            keep_auxiliary=False,
            chunk_size=None,
            **kwargs
    ):
        self.mod = 'IDW'
        self.mod_args = {
            'power':power,
            'min_points':min_points,
            'upper_limit':upper_limit,
            'lower_limit':lower_limit,
            'radius':radius,
            'supercede':supercede,
            'keep_auxiliary':keep_auxiliary,
            'chunk_size':chunk_size,
        }
        super().__init__(**kwargs)
        self.power = utils.float_or(power)
        self.min_points = utils.int_or(min_points)
        self.radius = np.inf if radius is None else utils.str2inc(radius) 
        self.upper_limit = utils.float_or(upper_limit)
        self.lower_limit = utils.float_or(lower_limit)
        self.supercede = supercede
        self.keep_auxiliary = keep_auxiliary
        self.chunk_size = chunk_size
        self.chunk_step = None
        
    def run(self):
        xcount, ycount, dst_gt = self.p_region.geo_transform(x_inc=self.xinc, y_inc=self.yinc)
        ds_config = demfun.set_infos(
            xcount,
            ycount,
            xcount * ycount,
            dst_gt,
            utils.sr_wkt(self.dst_srs),
            gdal.GDT_Float32,
            self.ndv,
            self.fmt
        )

        if self.verbose:
            if self.min_points:
                utils.echo_msg(
                    'generating IDW grid @ {}/{} looking for at least {} neighbors within {} pixels'.format(
                        ycount, xcount, self.min_points, self.radius
                    )
                )
            else:
                utils.echo_msg(
                    'generating IDW grid @ {}/{}'.format(ycount, xcount)
                )
            i=0

        self._stacks_array(
            out_name='{}_idw_stack'.format(self.name),
            supercede=self.supercede
        )
        n = '{}_idw_stack_s.tif'.format(self.name)
        w = '{}_idw_stack_w.tif'.format(self.name)
        c = '{}_idw_stack_c.tif'.format(self.name)

        if self.chunk_size is None:
            n_chunk = int(ds_config['nx'] * .1)
            n_chunk = 10 if n_chunk < 10 else n_chunk
        else: n_chunk = self.chunk_size
        if self.chunk_step is None:
            n_step = int(n_chunk/2)
        else: n_step = self.chunk_step

        points_ds = gdal.Open(n)
        points_band = points_ds.GetRasterBand(1)
        points_no_data = points_band.GetNoDataValue()
        
        weights_ds = gdal.Open(w)
        weights_band = weights_ds.GetRasterBand(1)
        weights_no_data = weights_band.GetNoDataValue()

        try:
            interp_ds = points_ds.GetDriver().Create(
                self.fn, points_ds.RasterXSize, points_ds.RasterYSize, bands=1, eType=points_band.DataType,
                options=["BLOCKXSIZE=256", "BLOCKYSIZE=256", "TILED=YES", "COMPRESS=LZW", "BIGTIFF=YES"]
            )
            interp_ds.SetProjection(points_ds.GetProjection())
            interp_ds.SetGeoTransform(points_ds.GetGeoTransform())
            interp_band = interp_ds.GetRasterBand(1)
            interp_band.SetNoDataValue(np.nan)
        except:
            return(self)

        points_array = points_band.ReadAsArray()
        point_indices = np.nonzero(points_array != points_no_data)
        point_values = points_array[point_indices]
        points_ds = points_array = None

        if self.weights:
            weights_array = weights_band.ReadAsArray()
            weight_values = weights_array[point_indices]
            weights_ds = weights_array = None
        else:
            weight_values = None

        invdisttree = Invdisttree(np.transpose(point_indices), point_values, leafsize=10, stat=1)
        for srcwin in utils.yield_srcwin((ycount, xcount), n_chunk=n_chunk):                
            if len(point_indices[0]):
                xi, yi = np.mgrid[srcwin[0]:srcwin[0]+srcwin[2],
                                  srcwin[1]:srcwin[1]+srcwin[3]]
                interp_data = invdisttree(
                    np.vstack((yi.flatten(), xi.flatten())).T,
                    nnear=self.min_points,
                    eps=.1,
                    p=self.power,
                    dub=self.radius,
                    weights=weight_values
                )
                interp_data = np.reshape(interp_data, (srcwin[2], srcwin[3]))
                interp_band.WriteArray(interp_data.T, srcwin[0], srcwin[1])                    
        interp_ds = point_values = weight_values = None
        if not self.keep_auxiliary:
            utils.remove_glob('{}*'.format(n), '{}*'.format(w), '{}*'.format(c))
        else:
            os.rename(w, '{}_w.tif'.format(self.name))
            os.rename(c, '{}_c.tif'.format(self.name))
            os.rename(n, '{}_n.tif'.format(self.name))
            self.aux_dems.append('{}_w.tif'.format(self.name))
            self.aux_dems.append('{}_c.tif'.format(self.name))
            self.aux_dems.append('{}_n.tif'.format(self.name))
        
        return(self)    
    
## ==============================================
## Scipy gridding (linear, cubic, nearest)
##
## https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.griddata.html
## TODO: rename this in some way
## ==============================================
class WafflesSciPy(Waffle):
    """Generate DEM using Scipy gridding interpolation
    
Generate a DEM using Scipy's gridding interpolation
Optional gridding methods are 'linear', 'cubic' and 'nearest'
https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.griddata.html
            
< scipy:method=<method>:supercede=False:keep_auxiliary=False:chunk_size=None:chunk_buffer=40 >

---
Parameters:
    
method=[linear/cubic/nearest] - interpolation method to use
supercede=[True/False] - supercede higher weighted data,
keep_auxiliary=[True/False] - retain auxiliary files
chunk_size=[val] - size of chunks in pixels
chunk_buffer=[val] - size of the chunk buffer in pixels
    """
    
    def __init__(
            self,
            method='linear',
            supercede=False,
            keep_auxiliary=False,
            chunk_size=None,
            chunk_buffer=40,
            **kwargs):
        """generate a `scipy` dem"""
        
        self.mod = 'scipy'
        self.mod_args = {
            'method':method,
            'supercede':supercede,
            'keep_auxiliary':keep_auxiliary,
            'chunk_size':chunk_size,
            'chunk_buffer':chunk_buffer,
        }
        super().__init__(**kwargs)
        self.methods = ['linear', 'cubic', 'nearest']
        self.method = method
        self.supercede = supercede
        self.keep_auxiliary = keep_auxiliary
        self.chunk_size = utils.int_or(chunk_size)
        self.chunk_step = None
        self.chunk_buffer = utils.int_or(chunk_buffer)

    def run(self):
        if self.method not in self.methods:
            utils.echo_error_msg(
                '{} is not a valid interpolation method, options are {}'.format(
                    self.method, self.methods
                )
            )
            return(self)
        
        xcount, ycount, dst_gt = self.p_region.geo_transform(x_inc=self.xinc, y_inc=self.yinc)
        ds_config = demfun.set_infos(
            xcount,
            ycount,
            xcount * ycount,
            dst_gt,
            utils.sr_wkt(self.dst_srs),
            gdal.GDT_Float32,
            self.ndv,
            self.fmt
        )

        self._stacks_array(
            out_name='{}_scipy_stack'.format(self.name),
            supercede=self.supercede
        )
        n = '{}_scipy_stack_s.tif'.format(self.name)
        w = '{}_scipy_stack_w.tif'.format(self.name)
        c = '{}_scipy_stack_c.tif'.format(self.name)

        if self.chunk_size is None:
            n_chunk = int(ds_config['nx'] * .1)
            n_chunk = 10 if n_chunk < 10 else n_chunk
        else: n_chunk = self.chunk_size
        if self.chunk_step is None:
            n_step = int(n_chunk/2)
        else: n_step = self.chunk_step

        points_ds = gdal.Open(n)
        points_band = points_ds.GetRasterBand(1)
        points_no_data = points_band.GetNoDataValue()
        
        try:
            interp_ds = points_ds.GetDriver().Create(
                self.fn, points_ds.RasterXSize, points_ds.RasterYSize, bands=1, eType=points_band.DataType,
                options=["BLOCKXSIZE=256", "BLOCKYSIZE=256", "TILED=YES", "COMPRESS=LZW", "BIGTIFF=YES"]
            )
            interp_ds.SetProjection(points_ds.GetProjection())
            interp_ds.SetGeoTransform(points_ds.GetGeoTransform())
            interp_band = interp_ds.GetRasterBand(1)
            interp_band.SetNoDataValue(np.nan)
        except:
            return(self)
        
        if self.verbose:
            utils.echo_msg('buffering srcwin by {} pixels'.format(self.chunk_buffer))
       
        for srcwin in utils.yield_srcwin((ycount, xcount), n_chunk=n_chunk, verbose=self.verbose): #, step=n_step):                
            srcwin_buff = utils.buffer_srcwin(srcwin, (ycount, xcount), self.chunk_buffer)
            points_array = points_band.ReadAsArray(*srcwin_buff)
            point_indices = np.nonzero(points_array != points_no_data)
            if len(point_indices[0]):
                point_values = points_array[point_indices]
                xi, yi = np.mgrid[0:srcwin_buff[2],
                                  0:srcwin_buff[3]]
                
                try:
                    interp_data = interpolate.griddata(
                        np.transpose(point_indices), point_values,
                        (xi, yi), method=self.method
                    )
                    y_origin = srcwin[1]-srcwin_buff[1]
                    x_origin = srcwin[0]-srcwin_buff[0]
                    y_size = y_origin + srcwin[3]
                    x_size = x_origin + srcwin[2]
                    interp_data = interp_data[y_origin:y_size,x_origin:x_size]
                    interp_band.WriteArray(interp_data, srcwin[0], srcwin[1])
                except Exception as e:
                    #print(e)
                    continue
        interp_ds = points_ds = point_values = weight_values = None
        if not self.keep_auxiliary:
            utils.remove_glob('{}*'.format(n), '{}*'.format(w), '{}*'.format(c))
        else:
            os.rename(w, '{}_w.tif'.format(self.name))
            os.rename(c, '{}_c.tif'.format(self.name))
            os.rename(n, '{}_n.tif'.format(self.name))
            self.aux_dems.append('{}_w.tif'.format(self.name))
            self.aux_dems.append('{}_c.tif'.format(self.name))
            self.aux_dems.append('{}_n.tif'.format(self.name))
            
        return(self)
    
## ==============================================
## Waffles VDatum
## ==============================================
class WafflesVDatum(Waffle):
    """VDATUM transformation grid.
Generate a Vertical DATUM transformation grid.

< vdatum:vdatum_in=None:vdatum_out=None >

    ---
    Parameters:
    
    vdatum_in=[vdatum] - input vertical datum
    vdatum_out=[vdatum] - output vertical datum
    """

    def __init__(self, vdatum_in=None, vdatum_out=None, **kwargs):
        self.mod = 'vdatum'
        self.mod_args = {
            'vdatum_in':vdatum_in,
            'vdatum_out':vdatum_out
        }
        super().__init__(**kwargs)
        self.vdatum_in = vdatum_in
        self.vdatum_out = vdatum_out
        
        import cudem.vdatums

    def run(self):
        cudem.vdatums.VerticalTransform(
            self.p_region,
            self.xinc,
            self.yinc,
            self.vdatum_in,
            self.vdatum_out,
            cache_dir=waffles_cache
        ).run(outfile='{}.tif'.format(self.name))
        
        return(self)
    
## ==============================================
## GDAL gridding (invdst, linear, nearest, average)
## ==============================================
class WafflesGDALGrid(Waffle):
    """Waffles GDAL_GRID module.

    see gdal_grid for more info and gridding algorithms
    """
    
    def __init__(self, **kwargs):
        """run gdal grid using alg_str

        parse the data through xyzfun.xyz_block to get weighted mean before
        building the GDAL dataset to pass into gdal_grid
        
        Args: 
          alg_str (str): the gdal_grid algorithm string
        """
        
        super().__init__(**kwargs)
        self.alg_str = 'linear:radius=-1'
        self.mod = self.alg_str.split(':')[0]

    def _xyz_ds(self):
        """Make a point vector OGR DataSet Object from src_xyz

        for use in gdal gridding functions.
        """

        dst_ogr = '{}'.format(self.name)
        self.ogr_ds = gdal.GetDriverByName('Memory').Create(
            '', 0, 0, 0, gdal.GDT_Unknown
        )
        layer = self.ogr_ds.CreateLayer(
            dst_ogr,
            geom_type=ogr.wkbPoint25D
        )
        fd = ogr.FieldDefn('long', ogr.OFTReal)
        fd.SetWidth(10)
        fd.SetPrecision(8)
        layer.CreateField(fd)
        fd = ogr.FieldDefn('lat', ogr.OFTReal)
        fd.SetWidth(10)
        fd.SetPrecision(8)
        layer.CreateField(fd)
        fd = ogr.FieldDefn('elev', ogr.OFTReal)
        fd.SetWidth(12)
        fd.SetPrecision(12)
        layer.CreateField(fd)
        
        if self.weights:
            fd = ogr.FieldDefn('weight', ogr.OFTReal)
            fd.SetWidth(6)
            fd.SetPrecision(6)
            layer.CreateField(fd)
            
        f = ogr.Feature(feature_def=layer.GetLayerDefn())        
        for this_xyz in self.yield_xyz():
            f.SetField(0, this_xyz.x)
            f.SetField(1, this_xyz.y)
            f.SetField(2, float(this_xyz.z))
            if self.weights:
                f.SetField(3, this_xyz.w)
                
            wkt = this_xyz.export_as_wkt(include_z=True)
            g = ogr.CreateGeometryFromWkt(wkt)
            f.SetGeometryDirectly(g)
            layer.CreateFeature(f)
            
        return(self.ogr_ds)
        
    def run(self):
        xcount, ycount, dst_gt = self.p_region.geo_transform(x_inc=self.xinc, y_inc=self.yinc)
        utils.echo_msg(
            'running GDAL GRID {} algorithm @ {} and {}/{}...'.format(
                self.alg_str.split(':')[0], self.p_region.format('fn'), xcount, ycount
            )
        )
        #_prog_update = lambda x, y, z: _prog.update()
        ds = self._xyz_ds()
        if ds.GetLayer().GetFeatureCount() == 0:
            utils.echo_error_msg('no input data')
            
        gd_opts = gdal.GridOptions(
            outputType = gdal.GDT_Float32,
            noData = self.ndv,
            format = 'GTiff',
            width = xcount,
            height = ycount,
            algorithm = self.alg_str,
            #callback = _prog_update if self.verbose else None,
            outputBounds = [
                self.p_region.xmin,
                self.p_region.ymax,
                self.p_region.xmax,
                self.p_region.ymin
            ]
        )
        gdal.Grid(
            '{}.tif'.format(self.name), ds, options = gd_opts
        )
        demfun.set_nodata(
            '{}.tif'.format(self.name, nodata=self.ndv, convert_array=False)
        )
        # _prog.end(
        #     0,
        #     'ran GDAL GRID {} algorithm @ {}.'.format(
        #         self.alg_str.split(':')[0], self.p_region.format('fn')
        #     )
        # )        
        ds = None
        return(self)

class WafflesLinear(WafflesGDALGrid):
    """LINEAR DEM via gdal_grid
Generate a DEM using GDAL's gdal_grid command.
see gdal_grid --help for more info

< linear:radius=-1 >

    ---
    Parameters:
    
    radius=[val] - search radius
    """
    
    def __init__(self, radius=None, nodata=-9999, **kwargs):
        self.mod = 'linear'
        self.mod_args = {
            'radius':radius,
            'nodata':nodata
        }
        super().__init__(**kwargs)        
        radius = self.xinc * 4 if radius is None else utils.str2inc(radius)
        self.alg_str = 'linear:radius={}:nodata={}'.format(radius, nodata)
        
class WafflesInvDst(WafflesGDALGrid):
    """INVERSE DISTANCE DEM via gdal_grid
    
Generate a DEM using GDAL's gdal_grid command.
see gdal_grid --help for more info

< invdst:power=2.0:smoothing=0.0:radius1=0:radius2=0:angle=0:max_points=0:min_points=0:nodata=0 >

    ---
    Parameters:
    
    radius1=[val] - search radius 1
    radius2=[val] - search radius 2
    power=[val] - weight**power
    min_points=[val] - minimum points per IDW bucket (use to fill entire DEM)
    """
    def __init__(self, power = 2.0, smoothing = 0.0, radius1 = None, radius2 = None, angle = 0.0,
                   max_points = 0, min_points = 0, nodata = -9999, **kwargs):
        self.mod = 'invdst'
        self.mod_args = {
            'power':power,
            'smoothing':smoothing,
            'radius1':radius1,
            'radius2':radius2,
            'angle':angle,
            'max_points':max_points,
            'min_points':min_points,
            'nodata':nodata
        }
        super().__init__(**kwargs)
        radius1 = self.xinc * 2 if radius1 is None else utils.str2inc(radius1)
        radius2 = self.yinc * 2 if radius2 is None else utils.str2inc(radius2)
        self.alg_str = 'invdist:power={}:smoothing={}:radius1={}:radius2={}:angle={}:max_points={}:min_points={}:nodata={}'\
            .format(power, smoothing, radius1, radius2, angle, max_points, min_points, nodata)
                
class WafflesMovingAverage(WafflesGDALGrid):
    """Moving AVERAGE DEM via gdal_grid
    
Generate a DEM using GDAL's gdal_grid command.
see gdal_grid --help for more info

< average:radius1=0:radius2=0:angle=0:min_points=0:nodata=0 >

    ---
    Parameters:
    
    radius1=[val] - search radius 1
    radius2=[val] - search radius 2
    min_points=[val] - minimum points per bucket (use to fill entire DEM)
    """
    
    def __init__(self, radius1=None, radius2=None, angle=0.0, min_points=0, nodata=-9999, **kwargs):
        self.mod = 'average'
        self.mod_args = {
            'radius1':radius1,
            'radius2':radius2,
            'angle':angle,
            'min_points':min_points,
            'nodata':nodata
        }
        super().__init__(**kwargs)
        radius1 = self.xinc * 2 if radius1 is None else utils.str2inc(radius1)
        radius2 = self.yinc * 2 if radius2 is None else utils.str2inc(radius2)
        self.alg_str = 'average:radius1={}:radius2={}:angle={}:min_points={}:nodata={}'\
            .format(radius1, radius2, angle, min_points, nodata)
                
class WafflesNearest(WafflesGDALGrid):
    """NEAREST DEM via gdal_grid
    
Generate a DEM using GDAL's gdal_grid command.
see gdal_grid --help for more info

< nearest:radius1=0:radius2=0:angle=0:nodata=0 >

    ---
    Parameters:
    
    radius1=[val] - search radius 1
    radius2=[val] - search radius 2
    angle=[val] - angle
    nodata=[val] - nodata
    """
    
    def __init__(self, radius1=None, radius2=None, angle=0.0, nodata=-9999, **kwargs):
        self.mod = 'nearest'
        self.mod_args = {
            'radius1':radius1,
            'radius2':radius2,
            'angle':angle,
            'nodata':nodata
        }
        super().__init__(**kwargs)
        radius1 = self.xinc * 2 if radius1 is None else utils.str2inc(radius1)
        radius2 = self.yinc * 2 if radius2 is None else utils.str2inc(radius2)
        self.alg_str = 'nearest:radius1={}:radius2={}:angle={}:nodata={}'\
            .format(radius1, radius2, angle, nodata)

## ==============================================
## Waffles 'CUDEM' gridding
##
## combined gridding method (stacks/surface/coastline/IDW)
## ==============================================
class WafflesCUDEM(Waffle):
    """CUDEM integrated DEM generation.
    
Generate an topo/bathy integrated DEM using a variety of data sources.
Will iterate <pre_count> pre-surfaces at lower-resolutions.
Each pre-surface will be clipped to <landmask> if it exists and smoothed with <smoothing> factor.
Each pre-surface is used in subsequent pre-surface(s)/final DEM at each iterative weight.

generate a DEM with `pre_surface`s which are generated
at lower resolution and with various weight threshholds.

To generate a typical CUDEM tile, generate 1 pre-surface ('bathy_surface'), clipped to a coastline.
Use a min_weight that excludes low-resolution bathymetry data from being used as input in the final
DEM generation. 
    
    < cudem:landmask=None:min_weight=[75th percentile]:smoothing=None:pre_count=1:mode=surface:supercede=True >

    ---
    Parameters:
    
    landmask=[path] - path to coastline vector mask or set as `coastline` to auto-generate
    min_weight=[val] - the minumum weight to include in the final DEM
    smoothing=[val] - the Gaussian smoothing value to apply to pre-surface(s)
    pre_count=[val] - number of pre-surface iterations to perform
    poly_count=[True/False]
    keep_auxiliary=[True/False]
    mode=surface
    filter_outliers=[val]
    supercede=[True/False]
    """
    
    def __init__(
            self,
            min_weight=None,
            pre_count=1,
            smoothing=None,
            landmask=False,
            poly_count=True,
            keep_auxiliary=False,
            mode='surface',
            filter_outliers=None,
            supercede=True,
            **kwargs
    ):
        self.mod = 'cudem'
        self.mod_args = {
            'min_weight':min_weight,
            'pre_count':pre_count,
            'smoothing':smoothing,
            'landmask':landmask,
            'poly_count':poly_count,
            'keep_auxiliary':keep_auxiliary,
            'mode':mode,
            'supercede':supercede,
            'filter_outliers':filter_outliers,
        }
        try:
            super().__init__(**kwargs)
        except Exception as e:
            utils.echo_error_msg(e)
            sys.exit()

        self.min_weight = utils.float_or(min_weight)
        self.pre_count = utils.int_or(pre_count, 1)
        self.smoothing = utils.int_or(smoothing)
        self.landmask = landmask
        self.poly_count = poly_count
        self.keep_auxiliary = keep_auxiliary
        self.mode = mode
        self.supercede = supercede
        self.filter_outliers = utils.int_or(filter_outliers)
        if self.filter_outliers is not None:
            self.filter_outliers = 1 if self.filter_outliers > 9 or self.filter_outliers < 1 else self.filter_outliers
        
    def run(self):
        pre = self.pre_count
        pre_weight = 0
        final_region = self.d_region.copy()
        pre_region = self.p_region.copy()
        pre_region.wmin = None
        pre_clip = None 
        upper_limit = None
        coast = '{}_cst'.format(self.name)
        
        ## Block/Stack the data with weights
        self._stacks_array(
            out_name='{}_stack'.format(self.name), supercede=self.supercede
        )
        n = '{}_stack_s.tif'.format(self.name)
        w = '{}_stack_w.tif'.format(self.name)
        c = '{}_stack_c.tif'.format(self.name)
        
        ## Remove outliers from the stacked data
        if self.filter_outliers is not None:
            demfun.filter_outliers_slp(
                n, '_tmp_fltr.tif', agg_level=self.filter_outliers, replace=True
            )
            os.rename('_tmp_fltr.tif', n)
            demfun.mask_(w, n, '_tmp_w.tif')
            os.rename('_tmp_w.tif', w)
            demfun.mask_(c, n, '_tmp_c.tif')
            os.rename('_tmp_c.tif', c)

        if self.min_weight is None:
            self.min_weight = demfun.percentile(w, 75)
        
        if self.verbose:
            utils.echo_msg('cudem min weight is: {}'.format(self.min_weight))
            
        pre_data = ['{},200:weight_mask={}:sample=average,1'.format(n, w)]

        ## Generate Coastline
        self.coast = None
        if self.landmask:
            upper_limit = -0.1
            if not os.path.exists(utils.str_or(self.landmask)):
               if os.path.exists('{}.shp'.format(utils.append_fn(self.landmask, self.region, self.xinc))):
                   coastline = '{}.shp'.format(utils.append_fn(self.landmask, self.region, self.xinc))
               else:
                   self.coast = WaffleFactory(
                       mod='coastline:polygonize=False',
                       data=pre_data,
                       src_region=self.p_region,
                       xinc=self.xinc,
                       yinc=self.yinc,
                       name=coast,
                       node=self.node,
                       weights=self.weights,
                       dst_srs=self.dst_srs,
                       srs_transform=self.srs_transform,
                       clobber=True,
                       verbose=self.verbose
                   ).acquire().generate()
                   coastline = '{}.shp'.format(self.coast.name)
            else:
               coastline = self.landmask
            pre_clip = '{}'.format(coastline)

        ## Grid/Stack the data `pre` times concluding in full resolution @ min_weight
        while pre >= 0:
            pre_xinc = float(self.xinc * (3**pre))
            pre_yinc = float(self.yinc * (3**pre))
            xsample = self.xinc * (3**(pre - 1))
            ysample = self.yinc * (3**(pre - 1))
            if xsample == 0: xsample = self.xinc
            if ysample == 0: ysample = self.yinc
            pre_filter=['1:{}'.format(self.smoothing)] if self.smoothing is not None else []
            
            ## if not final output, setup the configuration for the pre-surface
            if pre != self.pre_count:
                pre_weight = self.min_weight/(pre + 1) if pre > 0 else self.min_weight
                if pre_weight == 0: pre_weight = 1-e20
                pre_data = [
                    '{},200:weight_mask={}:sample=average,1'.format(n, w),
                    '{}.tif,200,{}'.format(
                        utils.append_fn('_pre_surface', pre_region, pre+1), pre_weight
                    )
                ]
                pre_region.wmin = pre_weight

            if pre == 0:
                pre_region.zmax = None

            pre_name = utils.append_fn('_pre_surface', pre_region, pre)
            if self.verbose:
                utils.echo_msg('pre region: {}'.format(pre_region))

            ## Grid/Stack the iteration
            if pre == self.pre_count:
                waffles_mod_surface = 'surface:tension=.65:upper_limit={}:blockmean=False:geographic=True:aspect=m'.format(upper_limit)
                waffles_mod_surfstack = 'surface:tension=.65:upper_limit={}:blockmean=False:geographic=True:aspect=m'.format(upper_limit)
            elif pre != 0:
                waffles_mod_surface = 'stacks:supercede=True:upper_limit={}'.format(upper_limit if upper_limit is not None else None)
                waffles_mod_surfstack = 'stacks:supercede=True:upper_limit={}'.format(upper_limit if upper_limit is not None else None)
            else:
                waffles_mod_surface = 'IDW'
                waffles_mod_surfstack = 'stacks:supercede=True'
            
            waffles_mod = waffles_mod_surface if self.landmask else waffles_mod_surfstack
            pre_surface = WaffleFactory(
                mod=waffles_mod,
                data=pre_data,
                src_region=pre_region if pre !=0 else final_region,
                xinc=pre_xinc if pre !=0 else self.xinc,
                yinc=pre_yinc if pre !=0 else self.yinc,
                name=pre_name if pre !=0 else self.name,
                node=self.node,
                fltr=pre_filter if pre !=0 else [],
                weights=1,
                dst_srs=self.dst_srs,
                srs_transform=self.srs_transform,
                clobber=True,
                verbose=self.verbose,
            ).acquire().generate()

            if self.coast is not None:
                if pre !=0:
                    demfun.mask_(pre_surface.fn, self.coast.fn, '__tmp_coast_clip.tif', msk_value=1)
                    os.rename('__tmp_coast_clip.tif', pre_surface.fn)

            pre -= 1
                    
        if not self.keep_auxiliary:
            utils.remove_glob('*_pre_surface*')
            utils.remove_glob('{}*'.format(n), '{}*'.format(c), '{}.*'.format(coast))

        os.rename(w, '{}_w.tif'.format(self.name))
        self.aux_dems.append('{}_w.tif'.format(self.name))
            
        return(self)
        
## ==============================================
## Waffles Raster Stacking
## ==============================================
class WafflesStacks_ETOPO(Waffle):
    """Waffles STACKS gridding module.

    stack data to generate DEM. No interpolation
    occurs with this module. To guarantee a full DEM,
    use a background DEM with a low weight, such as GMRT or GEBCO,
    which will be stacked upon to create the final DEM.

    XYZ input data is converted to raster; currently uses `mbgrid` to
    grid the XYZ data before stacking.

    set `supercede` to True to have higher weighted data overwrite 
    lower weighted data in overlapping cells. Default performs a weighted
    mean of overlapping data.

    note: this class was used for ETOPO2022, updated function is _stacks_array()
    """
    def __init__(
            self,
            supercede=False,
            keep_weights=False,
            keep_count=False,
            min_count=None,
            upper_limit=None,
            lower_limit=None,
            **kwargs
    ):
        self.mod = 'stacks'
        self.mod_args = {
            'supercede':supercede,
            'keep_weights':keep_weights,
            'keep_count':keep_count,
            'min_count':min_count,
            'upper_limit':upper_limit,
            'lower_limit':lower_limit
        }
        try:
            super().__init__(**kwargs)
        except Exception as e:
            utils.echo_error_msg(e)
            sys.exit()
        self.supercede = supercede
        self.keep_weights = keep_weights
        self.keep_count = keep_count
        self.min_count = min_count
        self.upper_limit = utils.float_or(upper_limit)
        self.lower_limit = utils.float_or(lower_limit)
        self._init_data(set_incs=True)
        if self.weights is None:
            self.weights = 1
        
    def run(self):
        xcount, ycount, dst_gt = self.p_region.geo_transform(
            x_inc=self.xinc, y_inc=self.yinc
        )
        gdt = gdal.GDT_Float32
        z_array = np.zeros((ycount, xcount))
        count_array = np.zeros((ycount, xcount))
        weight_array = np.zeros((ycount, xcount))

        if self.verbose:
            utils.echo_msg('stacking data to {}/{} grid using {} method'.format(
                ycount, xcount, 'supercede' if self.supercede else 'weighted mean'
            ))
            
        for arr, srcwin, gt, w in self.yield_array():
            if utils.float_or(w) is not None:
                w_arr = np.zeros((srcwin[3], srcwin[2]))
                w_arr[~np.isnan(arr)] = w if self.weights else 1
            else:
                w_arr = w

            w_arr[np.isnan(arr)] = 0
            w_arr[np.isnan(w_arr)] = 0
            c_arr = np.zeros((srcwin[3], srcwin[2]))
            c_arr[~np.isnan(arr)] = 1
            count_array[srcwin[1]:srcwin[1]+srcwin[3],srcwin[0]:srcwin[0]+srcwin[2]] += c_arr
            arr[np.isnan(arr)] = 0
            if not self.supercede:
                z_array[srcwin[1]:srcwin[1]+srcwin[3],srcwin[0]:srcwin[0]+srcwin[2]] += (arr * w_arr)
                weight_array[srcwin[1]:srcwin[1]+srcwin[3],srcwin[0]:srcwin[0]+srcwin[2]] += w_arr
            else:
                tmp_z = z_array[srcwin[1]:srcwin[1]+srcwin[3],srcwin[0]:srcwin[0]+srcwin[2]]
                tmp_w = weight_array[srcwin[1]:srcwin[1]+srcwin[3],srcwin[0]:srcwin[0]+srcwin[2]]
                tmp_z[w_arr > tmp_w] = arr[w_arr > tmp_w]
                tmp_w[w_arr > tmp_w] = w_arr[w_arr > tmp_w]
                z_array[srcwin[1]:srcwin[1]+srcwin[3],srcwin[0]:srcwin[0]+srcwin[2]] = tmp_z
                weight_array[srcwin[1]:srcwin[1]+srcwin[3],srcwin[0]:srcwin[0]+srcwin[2]] = tmp_w
            
        count_array[count_array == 0] = np.nan
        weight_array[weight_array == 0] = np.nan
        z_array[np.isnan(weight_array)] = np.nan
        if not self.supercede:
            weight_array = weight_array/count_array
            z_array = (z_array/weight_array)/count_array

        if self.min_count is not None:
            z_array[count_array < min_count] = np.nan
            weight_array[count_array < min_count] = np.nan

        if self.upper_limit is not None:
            z_array[z_array > self.upper_limit] = self.upper_limit
            
        if self.lower_limit is not None:
            z_array[z_array < self.lower_limit] = self.lower_limit
            
        z_array[np.isnan(z_array)] = self.ndv
        weight_array[np.isnan(weight_array)] = self.ndv
        count_array[np.isnan(count_array)] = self.ndv
        ds_config = demfun.set_infos(
            xcount,
            ycount,
            xcount * ycount,
            dst_gt,
            self.dst_srs,
            gdal.GDT_Float32,
            self.ndv,
            'GTiff'
        )
        utils.gdal_write(z_array, '{}.tif'.format(self.name), ds_config, verbose=True)
        if self.keep_weights:
            utils.gdal_write(weight_array, '{}_w.tif'.format(self.name), ds_config, verbose=True)
            self.aux_dems.append('{}_w.tif'.format(self.name))
            
        if self.keep_count:
            utils.gdal_write(count_array, '{}_c.tif'.format(self.name), ds_config, verbose=True)
            self.aux_dems.append('{}_c.tif'.format(self.name))

        if self.verbose:
            utils.echo_msg('stacked data to {}/{} grid using {} method'.format(
                ycount, xcount, 'supercede' if self.supercede else 'weighted mean'
            ))
            
        return(self)

class WafflesStacks(Waffle):
    """STACK data into a DEM.
    
Generate a DEM using a raster STACKing method. 
By default, will calculate the [weighted]-mean where overlapping cells occur. 
Set supercede to True to overwrite overlapping cells with higher weighted data.

stack data to generate DEM. No interpolation
occurs with this module. To guarantee a full DEM,
use a background DEM with a low weight, such as GMRT or GEBCO,
which will be stacked upon to create the final DEM.
    
< stacks:supercede=False:keep_weights=False:keep_count=False:min_count=None >

    ---
    Parameters:
    
    supercede=[True/False] - superced data cells with higher weighted data
    keep_weights=[True/False] - retain weight raster
    keep_count=[True/False] - retain count raster
    min_count=[val] - only retain data cells if they contain `min_count` overlapping data
    """
    
    def __init__(
            self,
            supercede=False,
            keep_weights=False,
            keep_count=False,
            min_count=None,
            upper_limit=None,
            lower_limit=None,
            **kwargs
    ):
        self.mod = 'stacks'
        self.mod_args = {
            'supercede':supercede,
            'keep_weights':keep_weights,
            'keep_count':keep_count,
            'min_count':min_count,
            'upper_limit':upper_limit,
            'lower_limit':lower_limit
        }
        try:
            super().__init__(**kwargs)
        except Exception as e:
            utils.echo_error_msg(e)
            sys.exit()
        self.supercede = supercede
        self.keep_weights = keep_weights
        self.keep_count = keep_count
        self.min_count = min_count
        self.upper_limit = utils.float_or(upper_limit)
        self.lower_limit = utils.float_or(lower_limit)
        self._init_data(set_incs=True)
        
    def run(self):
        self._stacks_array(
            out_name=self.name,
            supercede=self.supercede,
        )
        os.rename('{}_s.tif'.format(self.name), '{}.tif'.format(self.name))
        
        if self.keep_count:
            self.aux_dems.append('{}_c.tif'.format(self.name))
        else:
            utils.remove_glob('{}_c.tif'.format(self.name))
            
        if self.keep_weights:
            self.aux_dems.append('{}_w.tif'.format(self.name))
        else:
            utils.remove_glob('{}_w.tif'.format(self.name))
        
        return(self)

## ==============================================
## Waffles Coastline/Landmask
## ==============================================
class WafflesCoastline(Waffle):
    """COASTLINE (land/etc-mask) generation
    
Generate a coastline (land/etc-mask) using a variety of sources. 
User data can be provided to provide source for further land masking. 
Output raster will mask land-areas as 1 and oceans/(lakes/buildings) as 0.
Output vector will polygonize land-areas.

< coastline:want_gmrt=False:want_nhd=True:want_lakes=False:want_buildings=False:invert=False:polygonize=True >

    ---
    Parameters:
    
    want_gmrt=[True/False] - use GMRT to fill background (will use Copernicus otherwise)
    want_copernicus=[True/False] - use COPERNICUS to fill background
    want_nhd=[True/False] - use high-resolution NHD to fill US coastal zones
    want_lakes=[True/False] - mask LAKES using HYDROLAKES
    want_buildings=[True/False] - mask BUILDINGS using OSM
    osm_tries=[val] - OSM max server attempts
    min_building_length=[val] - only use buildings larger than val
    want_wsf=[True/False] - mask BUILDINGS using WSF
    invert=[True/False] - invert the output results
    polygonize=[True/False] - polygonize the output
    """
    
    def __init__(
            self,
            want_nhd=True,
            want_copernicus=True,
            want_gmrt=False,
            want_lakes=False,
            want_buildings=False,
            min_building_length=None,
            want_osm_planet=False,
            invert=False,
            polygonize=True,
            osm_tries=5,
            want_wsf=False,
            **kwargs
    ):
        """Generate a landmask from various sources.

        sources include:
        GMRT
        Copernicus
        NHD
        HydroLakes
        OSM Buildings
        *user-data*
        
        set polyginize to False to skip polyginizing the landmask, set
        polygonize to an integer to control how many output polygons will
        be generated.

        output raster mask will have 0 over water and 1 over land.
        """

        self.mod = 'coastline'
        self.mod_args = {
            'want_nhd':want_nhd,
            'want_copernicus':want_copernicus,
            'want_gmrt':want_gmrt,
            'want_lakes':want_lakes,
            'want_buildings':want_buildings,
            'want_wsf':want_wsf,
            'min_building_length':min_building_length,
            'want_osm_planet':want_osm_planet,
            'osm_tries':osm_tries,
            'invert':invert,
            'polygonize':polygonize
        }
        super().__init__(**kwargs)

        self.want_nhd = want_nhd
        self.want_gmrt = want_gmrt
        self.want_copernicus = want_copernicus
        self.want_lakes = want_lakes
        self.want_buildings = want_buildings
        self.want_wsf = want_wsf
        self.min_building_length = min_building_length
        self.want_osm_planet = want_osm_planet
        self.invert = invert
        self.polygonize = polygonize
        self.osm_tries = utils.int_or(osm_tries, 5)

        self.coast_array = None
        self.ds_config = None

        import cudem.fetches.copernicus
        import cudem.fetches.tnm
        import cudem.fetches.gmrt
        import cudem.fetches.hydrolakes
        import cudem.fetches.osm
        import cudem.fetches.wsf
        import cudem.fetches.utils

        self.f_region = self.p_region.copy()
        self.f_region.buffer(pct=5, x_inc=self.xinc, y_inc=self.yinc)
        self.f_region.src_srs = self.dst_srs
        self.wgs_region = self.f_region.copy()
        if self.dst_srs is not None:
            self.wgs_region.warp('epsg:4326')
        else:
            self.dst_srs = 'epsg:4326'            

        horz_epsg, vert_epsg = utils.epsg_from_input(self.dst_srs)
        self.cst_srs = osr.SpatialReference()
        self.cst_srs.SetFromUserInput('epsg:{}'.format(horz_epsg))
            
    def run(self):
        self._load_coast_mask()

        if self.want_gmrt:
            self._load_gmrt()
            
        if self.want_copernicus:
            self._load_copernicus()

        if self.want_nhd:
            self._load_nhd()

        if self.want_lakes:
            self._load_lakes()

        if self.want_buildings:
            self._load_bldgs()

        if self.want_wsf:
            self._load_wsf()
            
        if len(self.data) > 0:
            self._load_data()

        if self.verbose:
            utils.echo_msg(
                'finanlizing array for region {} at {} {}...'.format(
                    self.p_region.format('gmt'), self.ds_config['nx'], self.ds_config['ny']
                )
            )
        self._finalize_array()
        self._write_coast_array()
        if self.polygonize:
            if utils.int_or(self.polygonize) is not None:
                self._write_coast_poly(poly_count=self.polygonize)
            else:
                self._write_coast_poly()
            
        return(self)

    def _finalize_array(self):
        self.coast_array[self.coast_array > 0] = 1
        self.coast_array[self.coast_array <= 0] = 0
        if self.invert:
            self.coast_array[self.coast_array == 0] = 2
            self.coast_array[self.coast_array == 1] = 0
            self.coast_array[self.coast_array == 2] = 1
            
    def _load_coast_mask(self):
        """create a nodata grid"""
        
        xcount, ycount, gt = self.p_region.geo_transform(x_inc=self.xinc, y_inc=self.yinc)
        self.ds_config = demfun.set_infos(
            xcount,
            ycount,
            xcount * ycount,
            gt,
            utils.sr_wkt(self.cst_srs),
            gdal.GDT_Float32,
            self.ndv,
            'GTiff'
        )        
        self.coast_array = np.zeros( (ycount, xcount) )

    def _load_gmrt(self):
        """GMRT - Global low-res.

        Used to fill un-set cells.
        """
        
        this_gmrt = cudem.fetches.gmrt.GMRT(
            src_region=self.wgs_region, weight=self.weights, verbose=self.verbose, layer='topo', outdir=self.cache_dir
        )
        #this_gmrt._outdir = self.cache_dir
        this_gmrt.run()

        fr = cudem.fetches.utils.fetch_results(this_gmrt, want_proc=False)
        fr.daemon = True
        fr.start()
        fr.join()

        gmrt_tif = this_gmrt.results[0]
        gmrt_ds = demfun.generate_mem_ds(self.ds_config, name='gmrt')
        gdal.Warp(gmrt_ds, gmrt_tif[1], dstSRS=self.cst_srs, resampleAlg=self.sample)
        gmrt_ds_arr = gmrt_ds.GetRasterBand(1).ReadAsArray()
        gmrt_ds_arr[gmrt_ds_arr > 0] = 1
        gmrt_ds_arr[gmrt_ds_arr < 0] = 0
        self.coast_array += gmrt_ds_arr
        gmrt_ds = gmrt_ds_arr = None
        
    def _load_copernicus(self):
        """copernicus"""

        this_cop = cudem.fetches.copernicus.CopernicusDEM(
            src_region=self.wgs_region, weight=self.weights, verbose=self.verbose, datatype='1', outdir=self.cache_dir
        )
        #this_cop._outdir = self.cache_dir
        this_cop.run()

        fr = cudem.fetches.utils.fetch_results(this_cop, want_proc=False, check_size=False)
        fr.daemon = True
        fr.start()
        fr.join()
        
        for i, cop_tif in enumerate(this_cop.results):
            demfun.set_nodata(cop_tif[1], 0, verbose=False)
            cop_ds = demfun.generate_mem_ds(self.ds_config, name='copernicus')
            gdal.Warp(
                cop_ds, cop_tif[1], dstSRS=self.cst_srs, resampleAlg=self.sample,
                callback=False, srcNodata=0
            )
            cop_ds_arr = cop_ds.GetRasterBand(1).ReadAsArray()
            cop_ds_arr[cop_ds_arr != 0] = 1
            self.coast_array += cop_ds_arr
            cop_ds = cop_ds_arr = None

    def _load_wsf(self):
        """wsf"""

        this_wsf = cudem.fetches.wsf.WSF(
            src_region=self.wgs_region, weight=self.weights, verbose=self.verbose, outdir=self.cache_dir
        )
        #this_wsf._outdir = self.cache_dir
        this_wsf.run()

        fr = cudem.fetches.utils.fetch_results(this_wsf, want_proc=False, check_size=False)
        fr.daemon = True
        fr.start()
        fr.join()
       
        for i, wsf_tif in enumerate(this_wsf.results):
            demfun.set_nodata(wsf_tif[1], 0, verbose=False)
            wsf_ds = demfun.generate_mem_ds(self.ds_config, name='wsf')

            gdal.Warp(
                wsf_ds, wsf_tif[1], dstSRS=self.cst_srs, resampleAlg='cubicspline',
                callback=gdal.TermProgress, srcNodata=0, outputType=gdal.GDT_Float32
            )
            wsf_ds_arr = wsf_ds.GetRasterBand(1).ReadAsArray()
            wsf_ds_arr[wsf_ds_arr != 0 ] = -1
            self.coast_array += wsf_ds_arr
            wsf_ds = wsf_ds_arr = None
            
    def _load_nhd(self):
        """USGS NHD (HIGH-RES U.S. Only)
        Fetch NHD (NHD High/Plus) data from TNM to fill in near-shore areas. 
        High resoultion data varies by location...
        """

        this_tnm = cudem.fetches.tnm.TheNationalMap(
            src_region=self.wgs_region,
            weight=self.weights,
            verbose=self.verbose,
            where="Name LIKE '%Hydro%'",
            extents='HU-8 Subbasin,HU-4 Subregion',
            outdir=self.cache_dir
        )
        this_tnm.run()
        fr = cudem.fetches.utils.fetch_results(this_tnm, want_proc=False)
        fr.daemon = True
        fr.start()
        fr.join()

        tnm_ds = demfun.generate_mem_ds(self.ds_config, name='nhd')
        if len(this_tnm.results) > 0:
            for i, tnm_zip in enumerate(this_tnm.results):
                tnm_zips = utils.unzip(tnm_zip[1], self.cache_dir)
                gdb = '.'.join(tnm_zip[1].split('.')[:-1]) + '.gdb'
                utils.run_cmd(
                    'ogr2ogr -update -append nhdArea_merge.shp {} NHDArea -where "FType=312 OR FType=336 OR FType=445 OR FType=460 OR FType=537" -clipdst {} 2>/dev/null'.format(
                        gdb, self.p_region.format('ul_lr')
                    ),
                    verbose=self.verbose
                )
                utils.run_cmd(
                    'ogr2ogr -update -append nhdArea_merge.shp {} NHDPlusBurnWaterBody -clipdst {} 2>/dev/null'.format(
                        gdb, self.p_region.format('ul_lr')
                    ),
                    verbose=self.verbose
                )
                utils.run_cmd(
                    'ogr2ogr -update -append nhdArea_merge.shp {} NHDWaterBody -where "FType=493 OR FType=466" -clipdst {} 2>/dev/null'.format(
                        gdb, self.p_region.format('ul_lr')
                    ),
                    verbose=self.verbose
                )

            utils.run_cmd(
                'gdal_rasterize -burn 1 nhdArea_merge.shp nhdArea_merge.tif -te {} -ts {} {} -ot Int32'.format(
                    self.p_region.format('te'),
                    self.ds_config['nx'],
                    self.ds_config['ny'],
                ),
                verbose=self.verbose
            )

            try:
                gdal.Warp(tnm_ds, 'nhdArea_merge.tif', dstSRS=self.cst_srs, resampleAlg=self.sample)
            except:
                tnm_ds = None
                
            #tnm_ds = gdal.Open('nhdArea_merge.tif')
            if tnm_ds is not None:
                tnm_ds_arr = tnm_ds.GetRasterBand(1).ReadAsArray()
                tnm_ds_arr[tnm_ds_arr < 1] = 0
                self.coast_array -= tnm_ds_arr
                tnm_ds = tnm_ds_arr = None
                
            utils.remove_glob('nhdArea_merge.*')

    def _load_lakes(self):
        """HydroLakes -- Global Lakes"""
        
        this_lakes = cudem.fetches.hydrolakes.HydroLakes(
            src_region=self.wgs_region, weight=self.weights, verbose=self.verbose, outdir=self.cache_dir
        )
        #this_lakes._outdir = self.cache_dir
        this_lakes.run()
        
        fr = cudem.fetches.utils.fetch_results(this_lakes, want_proc=False)
        fr.daemon = True
        fr.start()
        fr.join()
        
        lakes_shp = None
        lakes_zip = this_lakes.results[0][1]
        lakes_shps = utils.unzip(lakes_zip, self.cache_dir)
        for i in lakes_shps:
            if i.split('.')[-1] == 'shp':
                lakes_shp = i

        lakes_ds = demfun.generate_mem_ds(self.ds_config, name='lakes')
        lakes_warp_ds = demfun.generate_mem_ds(self.ds_config, name='lakes')
        lk_ds = ogr.Open(lakes_shp)
        if lk_ds is not None:
            lk_layer = lk_ds.GetLayer()
            lk_layer.SetSpatialFilter(self.f_region.export_as_geom())
            gdal.RasterizeLayer(lakes_ds, [1], lk_layer, burn_values=[-1])
            gdal.Warp(lakes_warp_ds, lakes_ds, dstSRS=self.cst_srs, resampleAlg=self.sample)
            lakes_ds_arr = lakes_warp_ds.GetRasterBand(1).ReadAsArray()
            self.coast_array[lakes_ds_arr == -1] = 0
            lakes_ds = lk_ds = lakes_warp_ds = None
        else:
            utils.echo_error_msg('could not open {}'.format(lakes_shp))

    def _load_bldgs(self):
        """load buildings from OSM
        
        OSM has a size limit, so will chunk the region and
        do multiple osm calls
        """

        # bldg_ds = demfun.generate_mem_ds(self.ds_config, name='bldg')
        # bldg_warp_ds = demfun.generate_mem_ds(self.ds_config, name='bldg')
        this_osm = cudem.fetches.osm.OpenStreetMap(
            src_region=self.wgs_region, weight=self.weights, verbose=self.verbose,
            planet=self.want_osm_planet, chunks=True, q='buildings', fmt='osm',
            min_length=self.min_building_length, outdir=self.cache_dir
        )
        #this_osm._outdir = self.cache_dir
        this_osm.run()

        # fr = cudem.fetches.utils.fetch_results(this_osm, want_proc=False, check_size=False)
        # fr.daemon = True
        # fr.start()
        # fr.join()
        
        os.environ["OGR_OSM_OPTIONS"] = "INTERLEAVED_READING=YES"
        os.environ["OGR_OSM_OPTIONS"] = "OGR_INTERLEAVED_READING=YES"

        with tqdm(total=len(this_osm.results), desc='processing OSM buildings') as pbar:
            for n, osm_result in enumerate(this_osm.results):
                pbar.update()
                if cudem.fetches.utils.Fetch(osm_result[0], verbose=True).fetch_file(osm_result[1], check_size=False, tries=self.osm_tries, read_timeout=3600) >= 0:
                    #if True:
                    if osm_result[-1] == 'bz2':
                        osm_planet = utils.unbz2(osm_result[1], self.cache_dir)
                        osm_file = utils.ogr_clip(osm_planet, self.wgs_region)
                        _clipped = True
                    elif osm_result[-1] == 'pbf':
                        osm_file = utils.ogr_clip(osm_result[1], self.wgs_region, 'multipolygons')
                        #osm_file = utils.ogr_clip(osm_result[1], self.wgs_region)
                        _clipped = True
                    else:
                        osm_file = osm_result[1]
                        _clipped = False

                    # utils.run_cmd(
                    #     'gdal_rasterize -burn -1 -l multipolygons {} bldg_osm.tif -where "building!=\'\'" -te {} -ts {} {} -ot Int32'.format(
                    #         osm_file,
                    #         self.p_region.format('te'),
                    #         self.ds_config['nx'],
                    #         self.ds_config['ny'],
                    #     ),
                    #     verbose=True
                    # )

                    # osm_ds = ogr.Open(osm_file)
                    # osm_layer = osm_ds.GetLayer('multipolygons')
                    # #print(osm_layer.GetFeatureCount())
                    # osm_layer_dfn = osm_layer.GetLayerDefn()
                    # print(osm_layer_dfn)
                    # #osm_fi = osm_layer_dfn.GetFieldIndex()
                    # #print(osm_fi)
                    # osm_fi = osm_layer_dfn.GetFieldIndex('building')
                    # print(osm_fi)
                    # osm_ds = None

                    if os.path.getsize(osm_file) == 366:
                        continue

                    out, status = utils.run_cmd(
                        'gdal_rasterize -burn -1 -l multipolygons {} bldg_osm.tif -te {} -ts {} {} -ot Int32 -q'.format(
                            osm_file,
                            self.p_region.format('te'),
                            self.ds_config['nx'],
                            self.ds_config['ny'],
                        ),
                        verbose=False
                    )

                    if status == 0:
                        bldg_ds = gdal.Open('bldg_osm.tif')
                        #bldg_ds = gdal.Warp('', 'bldg_osm.tif', format='MEM', dstSRS=self.cst_srs, resampleAlg=self.sample, callback=gdal.TermProgress)
                        if bldg_ds is not None:
                            bldg_ds_arr = bldg_ds.GetRasterBand(1).ReadAsArray()
                            self.coast_array[bldg_ds_arr == -1] = 0
                            bldg_ds = bldg_ds_arr = None

                        bldg_ds = None

                    utils.remove_glob('bldg_osm.tif*')

                    # osm_ds = ogr.Open(osm_file)
                    # if osm_ds is not None:
                    #     osm_layer = osm_ds.GetLayer('multipolygons')
                    #     #osm_layer.SetSpatialFilter(self.wgs_region.export_as_geom())
                    #     osm_layer.SetAttributeFilter("building!=''")
                    #     gdal.RasterizeLayer(bldg_ds, [1], osm_layer, burn_values=[-1])
                    #     gdal.Warp(bldg_warp_ds, bldg_ds, dstSRS=self.cst_srs, resampleAlg=self.sample)
                    #     bldg_arr = bldg_warp_ds.GetRasterBand(1).ReadAsArray()
                    #     self.coast_array[bldg_arr == -1] = 0
                    # else:
                    #     utils.echo_error_msg('could not open ogr dataset {}'.format(osm_file))
                    # osm_ds = None
                    #if _clipped:
                    #    utils.remove_glob(osm_file)
                #else:
                #/    utils.echo_error_msg('failed to fetch {}'.format(osm_result[0]))

        bldg_ds = bldg_warp_ds = None
        
    def _load_data(self):
        """load data from user datalist and amend coast_array"""

        for this_data in self.data:
            #print(this_data)
            for this_arr in this_data.yield_array():
                data_arr = this_arr[0]['z']
                srcwin = this_arr[1]
                data_arr[np.isnan(data_arr)] = 0
                data_arr[data_arr > 0] = 1
                data_arr[data_arr < 0] = -1
                self.coast_array[srcwin[1]:srcwin[1]+srcwin[3],
                                 srcwin[0]:srcwin[0]+srcwin[2]] += data_arr

    def _write_coast_array(self):
        """write coast_array to file"""

        if self.verbose:
            utils.echo_msg('writing array to {}.tif...'.format(self.name))
            
        utils.gdal_write(
            self.coast_array, '{}.tif'.format(self.name), self.ds_config,
        )

    def _write_coast_poly(self, poly_count=None):
        """convert to coast_array vector"""

        if self.verbose:
            utils.echo_msg('polygonizing array to {}.shp...'.format(self.name))
        
        poly_count = utils.int_or(poly_count)
        tmp_ds = ogr.GetDriverByName('ESRI Shapefile').CreateDataSource(
            'tmp_c_{}.shp'.format(self.name)
        )
        
        if tmp_ds is not None:
            tmp_layer = tmp_ds.CreateLayer('tmp_c_{}'.format(self.name), None, ogr.wkbMultiPolygon)
            tmp_layer.CreateField(ogr.FieldDefn('DN', ogr.OFTInteger))
            demfun.polygonize('{}.tif'.format(self.name), tmp_layer, verbose=self.verbose)
            tmp_ds = None
            
        utils.run_cmd(
            'ogr2ogr -dialect SQLITE -sql "SELECT * FROM tmp_c_{} WHERE DN=0 {}" {}.shp tmp_c_{}.shp'.format(
                self.name, 'order by ST_AREA(geometry) desc limit {}'.format(poly_count) if poly_count is not None else '', self.name, self.name),
            verbose=self.verbose
        )
        utils.remove_glob('tmp_c_{}.*'.format(self.name))
        utils.run_cmd(
            'ogrinfo -dialect SQLITE -sql "UPDATE {} SET geometry = ST_MakeValid(geometry)" {}.shp'.format(
                self.name, self.name),
            verbose=self.verbose
        )
        
        #utils.gdal_prj_file(self.name + '.prj', self.cst_srs)
        print(self.dst_srs)
        print(self.cst_srs.ExportToProj4())
        utils.gdal_prj_file(self.name + '.prj', self.dst_srs)

## ==============================================
## Waffles Lakes Bathymetry
## ==============================================
class WafflesLakes(Waffle):
    """Estimate lake bathymetry.
    
By default, will return lake bathymetry as depth values (positive down), 
to get elevations (positive up), set apply_elevations=True.

< lakes:apply_elevations=False:min_area=None:max_area=None:min_id=None:max_id=None:depth=globathy >

    ---
    Parameters:
    
    apply_elevations=[True/False] - use COPERNICUS to apply lake level elevations to output
    min_area=[val] - minimum lake area to consider
    max_area=[val] - maximum lake area to consider
    min_id=[val] - minimum lake ID to consider
    max_id=[val] - maximum lake ID to consider
    depth=[globathy/hydrolakes/val] - obtain the depth value from GloBathy, HydroLakes or constant value
    """
    
    def __init__(
            self,
            apply_elevations=False,
            min_area=None,
            max_area=None,
            min_id=None,
            max_id=None,
            depth='globathy',
            elevations='copernicus',
            **kwargs
    ):
        self.mod = 'lakes'
        self.mod_args = {}

        super().__init__(**kwargs)
        self._mask = None
        self.apply_elevations = apply_elevations
        self.min_area = min_area
        self.max_area = max_area
        self.min_id = min_id
        self.max_id = max_id
        self.depth = depth
        self.elevations = elevations
        self.ds_config = None
            
        self.wgs_region = self.p_region.copy()
        if self.dst_srs is not None:
            self.wgs_region.warp('epsg:4326')            
        else:
            self.dst_srs = 'epsg:4326'

    def _init_bathy(self):
        """create a nodata grid"""
        
        xcount, ycount, gt = self.p_region.geo_transform(x_inc=self.xinc, y_inc=self.yinc)
        self.ds_config = demfun.set_infos(
            xcount,
            ycount,
            xcount * ycount,
            gt,
            utils.sr_wkt(self.dst_srs),
            gdal.GDT_Float32,
            self.ndv,
            'GTiff'
        )
            
    def _fetch_lakes(self):
        """fetch hydrolakes polygons"""

        this_lakes = cudem.fetches.hydrolakes.HydroLakes(
            src_region=self.p_region, weight=self.weights, verbose=self.verbose, outdir=self.cache_dir
        )
        #this_lakes._outdir = self.cache_dir
        this_lakes.run()

        fr = cudem.fetches.utils.fetch_results(
            this_lakes, want_proc=False, check_size=False
        )
        fr.daemon = True
        fr.start()
        fr.join()

        lakes_shp = None
        lakes_zip = this_lakes.results[0][1]
        lakes_shps = utils.unzip(lakes_zip, self.cache_dir)
        for i in lakes_shps:
            if i.split('.')[-1] == 'shp':
                lakes_shp = i

        return(lakes_shp)

    def _fetch_globathy(self, ids=[]):
        """fetch globathy csv data and process into dict"""
        
        import csv
        import cudem.fetches.utils
        
        _globathy_url = 'https://springernature.figshare.com/ndownloader/files/28919991'
        globathy_zip = os.path.join(self.cache_dir, 'globathy_parameters.zip')
        cudem.fetches.utils.Fetch(
            _globathy_url, verbose=self.verbose
        ).fetch_file(
            globathy_zip, check_size=False
        )
        globathy_csvs = utils.unzip(globathy_zip, self.cache_dir)        
        globathy_csv = os.path.join(self.cache_dir, 'GLOBathy_basic_parameters/GLOBathy_basic_parameters(ALL_LAKES).csv')
        with open(globathy_csv, mode='r') as globc:
            reader = csv.reader(globc)
            next(reader)
            if len(ids) > 0:
                globd = {}
                for row in reader:
                    if int(row[0]) in ids:
                        globd[int(row[0])] = float(row[-1])
                        ids.remove(int(row[0]))
                        
                    if len(ids) == 0:
                        break
            else:
                globd = {int(row[0]):float(row[-1]) for row in reader}

        return(globd)

    def _fetch_gmrt(self, gmrt_region=None):
        """GMRT - Global low-res.
        """

        if gmrt_region is None:
            gmrt_region = self.p_region
        
        this_gmrt = cudem.fetches.gmrt.GMRT(
            src_region=gmrt_region, weight=self.weights, verbose=self.verbose, layer='topo', outdir=self.cache_dir
        )
        #this_gmrt._outdir = self.cache_dir
        this_gmrt.run()

        fr = cudem.fetches.utils.fetch_results(this_gmrt, want_proc=False)
        fr.daemon = True
        fr.start()
        fr.join()

        dst_srs = osr.SpatialReference()
        dst_srs.SetFromUserInput(self.dst_srs)
        
        gmrt_tif = this_gmrt.results[0]
        gmrt_ds = demfun.generate_mem_ds(self.ds_config, name='gmrt')
        gdal.Warp(gmrt_ds, gmrt_tif[1], dstSRS=dst_srs, resampleAlg=self.sample)
        return(gmrt_ds)
    
    def _fetch_copernicus(self, cop_region=None):
        """copernicus"""

        if cop_region is None:
            cop_region = self.p_region
            
        this_cop = cudem.fetches.copernicus.CopernicusDEM(
            src_region=cop_region, weight=self.weights, verbose=self.verbose, datatype='1', outdir=self.cache_dir
        )
        #this_cop._outdir = self.cache_dir
        this_cop.run()

        fr = cudem.fetches.utils.fetch_results(this_cop, want_proc=False, check_size=False)
        fr.daemon = True
        fr.start()
        fr.join()

        dst_srs = osr.SpatialReference()
        dst_srs.SetFromUserInput(self.dst_srs)
        cop_ds = demfun.generate_mem_ds(self.ds_config, name='copernicus')
        [gdal.Warp(cop_ds, cop_tif[1], dstSRS=dst_srs, resampleAlg=self.sample) for cop_tif in this_cop.results]
        
        return(cop_ds)
    
    def generate_mem_ogr(self, geom, srs):
        """Create temporary polygon vector layer from feature geometry 
        so that we can rasterize it (Rasterize needs a layer)
        """
        
        ds = ogr.GetDriverByName('MEMORY').CreateDataSource('')
        Layer = ds.CreateLayer('', geom_type=ogr.wkbPolygon, srs=srs)
        outfeature = ogr.Feature(Layer.GetLayerDefn())
        outfeature.SetGeometry(geom)
        Layer.SetFeature(outfeature)

        return(ds)

    ##
    ## Adapted from GLOBathy
    def apply_calculation(self, shore_distance_arr, lake_depths_arr, shore_arr=None):
        """
        Apply distance calculation, which is each pixel's distance to shore, multiplied
        by the maximum depth, all divided by the maximum distance to shore. This provides
        a smooth slope from shore to lake max depth.

        shore_distance - Input numpy array containing array of distances to shoreline.
            Must contain positive values for distance away from shore and 0 elsewhere.
        max_depth - Input value with maximum lake depth.
        NoDataVal - value to assign to all non-lake pixels (zero distance to shore).
        """

        labels, nfeatures = ndimage.label(shore_distance_arr)
        nfeatures = np.arange(1, nfeatures +1)
        maxes = ndimage.maximum(shore_distance_arr, labels, nfeatures)
        max_dist_arr = np.zeros(np.shape(shore_distance_arr))
        with tqdm(total=len(nfeatures), desc='applying labels...') as pbar:
            for n, x in enumerate(nfeatures):
                pbar.update()
                max_dist_arr[labels==x] = maxes[x-1]
        
        max_dist_arr[max_dist_arr == 0] = np.nan
        bathy_arr = (shore_distance_arr * lake_depths_arr) / max_dist_arr
        bathy_arr[bathy_arr == 0] = np.nan

        if shore_arr is None \
           or shore_arr.size == 0 \
           or shore_arr[~np.isnan(bathy_arr)].size == 0 \
           or shore_arr[~np.isnan(bathy_arr)].max() == 0:
            utils.echo_warning_msg('invalid shore array, using default shore value of zero')
            bathy_arr = 0 - bathy_arr
        else:
            bathy_arr = shore_arr - bathy_arr

        utils.echo_msg('applied shore elevations to lake depths')
        bathy_arr[np.isnan(bathy_arr)] = 0
        return(bathy_arr)    
    
    def run(self):
        #self.p_region.buffer(pct=2)
        
        lakes_shp = self._fetch_lakes()
        lk_ds = ogr.Open(lakes_shp, 1)
        lk_layer = lk_ds.GetLayer()

        ## filter layer to region
        filter_region = self.p_region.copy()
        lk_layer.SetSpatialFilter(filter_region.export_as_geom())

        ## filter by ID
        if self.max_id is not None:
            lk_layer.SetAttributeFilter('Hylak_id < {}'.format(self.max_id))
            
        if self.min_id is not None:
            lk_layer.SetAttributeFilter('Hylak_id > {}'.format(self.min_id))

        ## filter by Area
        if self.max_area is not None:
            lk_layer.SetAttributeFilter('Lake_area < {}'.format(self.max_area))
            
        if self.min_area is not None:
            lk_layer.SetAttributeFilter('Lake_area > {}'.format(self.min_area))
        
        lk_features = lk_layer.GetFeatureCount()
        lk_regions = self.p_region.copy()
        with tqdm(total=len(lk_layer), desc='processing {} lakes'.format(lk_features)) as pbar:
            for lk_f in lk_layer:
                pbar.update()
                this_region = regions.Region()
                lk_geom = lk_f.GetGeometryRef()
                lk_wkt = lk_geom.ExportToWkt()
                this_region.from_list(ogr.CreateGeometryFromWkt(lk_wkt).GetEnvelope())
                lk_regions = regions.regions_merge(lk_regions, this_region)

        while not regions.regions_within_ogr_p(self.p_region, lk_regions):
            utils.echo_msg('buffering region by 2 percent to gather all lake boundaries...')
            self.p_region.buffer(pct=2, x_inc=self.xinc, y_inc=self.yinc)

        self._init_bathy()
            
        ## fetch and initialize the copernicus data
        if self.elevations == 'copernicus':
            cop_ds = self._fetch_copernicus(cop_region=self.p_region)
            cop_band = cop_ds.GetRasterBand(1)
        elif self.elevations == 'gmrt':
            cop_ds = self._fetch_gmrt(gmrt_region=self.p_region)
            cop_band = cop_ds.GetRasterBand(1)
        elif utils.float_or(self.elevations) is not None:
            cop_band = None
            cop_arr = np.zeros((self.ds_config['nx'], self.ds_config['ny']))
            cop_arr[:] = self.elevations
        elif os.path.exists(self.elevations):
            elev_ds = gdal.Open(self.elevations)
            if elev_ds is not None:
                dst_srs = osr.SpatialReference()
                dst_srs.SetFromUserInput(self.dst_srs)
                cop_ds = demfun.generate_mem_ds(self.ds_config, name='cop')
                gdal.Warp(cop_ds, elev_ds, dstSRS=dst_srs, resampleAlg=self.sample)
                cop_band = cop_ds.GetRasterBand(1)
            else:
                cop_band = None
                cop_arr = None
        else:
            cop_band = None
            cop_arr = None
        
        ## initialize the tmp datasources
        prox_ds = demfun.generate_mem_ds(self.ds_config, name='prox')
        msk_ds = demfun.generate_mem_ds(self.ds_config, name='msk')
        msk_band = None
        globd = None
        
        ## get lake ids and globathy depths
        lk_ids = []
        [lk_ids.append(feat.GetField('Hylak_id')) for feat in lk_layer]

        utils.echo_msg('using Lake IDS: {}'.format(lk_ids))
        if self.depth == 'globathy':
            globd = self._fetch_globathy(ids=lk_ids[:])
            ## rasterize hydrolakes using id
            gdal.RasterizeLayer(msk_ds, [1], lk_layer, options=["ATTRIBUTE=Hylak_id"], callback=gdal.TermProgress)
            msk_ds.FlushCache()
            msk_band = msk_ds.GetRasterBand(1)
            msk_band.SetNoDataValue(self.ds_config['ndv'])

            ## assign max depth from globathy
            msk_arr = msk_band.ReadAsArray()
            with tqdm(total=len(lk_ids), desc='Assigning Globathy Depths to rasterized lakes...') as pbar:
                for n, this_id in enumerate(lk_ids):
                    depth = globd[this_id]
                    msk_arr[msk_arr == this_id] = depth
                    pbar.update()
            
        elif self.depth == 'hydrolakes':
            gdal.RasterizeLayer(msk_ds, [1], lk_layer, options=["ATTRIBUTE=Depth_avg"], callback=gdal.TermProgress)
            msk_ds.FlushCache()
            msk_band = msk_ds.GetRasterBand(1)
            msk_band.SetNoDataValue(self.ds_config['ndv'])
            msk_arr = msk_band.ReadAsArray()
            
        elif utils.float_or(self.depth) is not None:
            msk_arr = np.zeros((self.ds_config['nx'], self.ds_config['ny']))
            msk_arr[:] = self.depth
            
        else:            
            msk_arr = np.ones((self.ds_config['nx'], self.ds_config['ny']))
            
        ## calculate proximity of lake cells to shore
        if msk_band is None:
            gdal.RasterizeLayer(msk_ds, [1], lk_layer, options=["ATTRIBUTE=Hylak_id"], callback=gdal.TermProgress)
            msk_ds.FlushCache()
            msk_band = msk_ds.GetRasterBand(1)
            msk_band.SetNoDataValue(self.ds_config['ndv'])
            
        lk_ds = None
        prox_band = prox_ds.GetRasterBand(1)
        proximity_options = ["VALUES=0", "DISTUNITS=PIXEL"]
        gdal.ComputeProximity(msk_band, prox_band, options=proximity_options, callback=gdal.TermProgress)        
        prox_arr = prox_band.ReadAsArray()
        if cop_band is not None:
            cop_arr = cop_band.ReadAsArray()

        utils.echo_msg('Calculating simulated lake depths...')
        ## apply calculation from globathy
        bathy_arr = self.apply_calculation(
            prox_arr,
            msk_arr,
            shore_arr=cop_arr,
        )

        bathy_arr[bathy_arr == 0] = self.ndv        
        utils.gdal_write(
            bathy_arr, '{}.tif'.format(self.name), self.ds_config,
        )            

        prox_ds = msk_ds = cop_ds = None
        return(self)

    def apply_calculation2(self, shore_distance_arr, max_depth=0, shore_arr=None, shore_elev=0):
        """
        Apply distance calculation, which is each pixel's distance to shore, multiplied
        by the maximum depth, all divided by the maximum distance to shore. This provides
        a smooth slope from shore to lake max depth.

        shore_distance - Input numpy array containing array of distances to shoreline.
            Must contain positive values for distance away from shore and 0 elsewhere.
        max_depth - Input value with maximum lake depth.
        NoDataVal - value to assign to all non-lake pixels (zero distance to shore).
        """

        max_dist = float(1e-7) if shore_distance_arr.max() <= 0 else shore_distance_arr.max()
        bathy_arr = (shore_distance_arr * max_depth) / float(max_dist)
        bathy_arr[bathy_arr == 0] = np.nan

        if shore_arr is None \
           or shore_arr.size == 0 \
           or shore_arr[~np.isnan(bathy_arr)].size == 0 \
           or shore_arr[~np.isnan(bathy_arr)].max() == 0:
            bathy_arr = shore_elev - bathy_arr
        else:
            bathy_arr = shore_arr - bathy_arr
            
        bathy_arr[np.isnan(bathy_arr)] = 0
        return(bathy_arr)
    
    def run2(self):        
        self._load_bathy()        
        lakes_shp = self._fetch_lakes()
        cop_ds = self._fetch_copernicus()
        cop_band = cop_ds.GetRasterBand(1)
        prox_ds = demfun.generate_mem_ds(self.ds_config, name='prox')
        msk_ds = demfun.generate_mem_ds(self.ds_config, name='msk')

        lk_ds = ogr.Open(lakes_shp)
        lk_layer = lk_ds.GetLayer()
        lk_layer.SetSpatialFilter(self.p_region.export_as_geom())
        lk_features = lk_layer.GetFeatureCount()

        lk_ids = []
        for i, feat in enumerate(lk_layer):
            lk_ids.append(feat.GetField('Hylak_id'))

        globd = self._fetch_globathy(ids=lk_ids)

        with tqdm(total=len(lk_layer), desc='processing lakes') as pbar:
            for i, feat in enumerate(lk_layer):
                pbar.update()
                lk_id = feat.GetField('Hylak_id')
                lk_elevation = feat.GetField('Elevation')
                lk_depth_glb = globd[lk_id]

                lk_name = feat.GetField('Lake_name')
                lk_depth = feat.GetField('Depth_avg')
                lk_area = feat.GetField('Lake_area')
                utils.echo_msg('lake: {}'.format(lk_name))
                utils.echo_msg('lake elevation: {}'.format(lk_elevation))
                utils.echo_msg('lake_depth: {}'.format(lk_depth))
                utils.echo_msg('lake_area: {}'.format(lk_area))
                utils.echo_msg('cell_area: {}'.format(self.xinc*self.yinc))
                utils.echo_msg('lake_depth (globathy): {}'.format(lk_depth_glb))

                tmp_ds = self.generate_mem_ogr(feat.GetGeometryRef(), lk_layer.GetSpatialRef())
                tmp_layer = tmp_ds.GetLayer()            
                gdal.RasterizeLayer(msk_ds, [1], tmp_layer, burn_values=[1])            
                msk_band = msk_ds.GetRasterBand(1)
                msk_band.SetNoDataValue(self.ds_config['ndv'])

                prox_band = prox_ds.GetRasterBand(1)
                proximity_options = ["VALUES=0", "DISTUNITS=PIXEL"]
                gdal.ComputeProximity(msk_band, prox_band, options=proximity_options)

                prox_arr = prox_band.ReadAsArray()
                msk_arr = msk_band.ReadAsArray()

                self.bathy_arr += self.apply_calculation(
                    prox_band.ReadAsArray(),
                    max_depth=lk_depth_glb,
                    shore_arr=cop_band.ReadAsArray(),
                    shore_elev=lk_elevation
                )

                prox_arr[msk_arr == 1] = 0
                prox_ds.GetRasterBand(1).WriteArray(prox_arr)
                msk_arr[msk_arr == 1] = 0
                msk_ds.GetRasterBand(1).WriteArray(msk_arr)
                tmp_ds = None
            
        self.bathy_arr[self.bathy_arr == 0] = self.ndv        
        utils.gdal_write(
            self.bathy_arr, '{}.tif'.format(self.name), self.ds_config,
        )            
        prox_ds = msk_ds = None
        return(self)    
    
## ==============================================
## Waffles DEM update - testing
## ==============================================
class WafflesPatch(Waffle):
    """PATCH an existing DEM with new data.
    
Patch an existing DEM with data from the datalist.

< patch:min_weight=.5:dem=None >

    ---
    Parameters:
    
    dem=[path] - the path the the DEM to update
    min_weight=[val] - the minumum data weight to include in the patched DEM
    """
    
    def __init__(
            self,
            radius=None,
            min_weight=1,
            max_diff=.25,
            dem=None,
            **kwargs
    ):
        self.mod = 'update'
        self.mod_args = {
            'radius':radius,
            'min_weight':min_weight,
            'max_diff':max_diff,
            'dem':dem
        }
        try:
            super().__init__(**kwargs)
        except Exception as e:
            utils.echo_error_msg(e)
            sys.exit()

        self.radius = utils.str2inc(radius)
        self.min_weight = utils.float_or(min_weight)
        self.max_diff = utils.float_or(max_diff)

        if dem is not None:
            if os.path.exists(dem):
                self.dem = dem
        elif os.path.exists('{}.tif'.format(self.name)):
            self.dem = '{}.tif'.format(self.name)
            self.name = '{}_update'.format(self.name)
        else:
            utils.echo_error_msg('must specify DEM to patch (:dem=fn) to run the patch module.')
            return(None)

    def yield_diff(self, src_dem, max_diff=.25):
        '''query a gdal-compatible grid file with xyz data.
        out_form dictates return values

        yields out_form results'''

        def con_dec(x, dec):
            '''Return a float string with n decimals
            (used for ascii output).'''

            if x is None:
                utils.echo_error_msg('Attempting to convert a None value.')
                return
            fstr = "%." + str(dec) + "f"
            return fstr % x
      
        try:
            ds = gdal.Open(src_dem)
        except: ds = None
        if ds is not None:
            ds_config = demfun.gather_infos(ds)
            ds_band = ds.GetRasterBand(1)
            ds_gt = ds_config['geoT']
            ds_nd = ds_config['ndv']
            tgrid = ds_band.ReadAsArray()
            dsband = ds = None

            for xyz in self.yield_xyz():
                #if xyz.x > ds_gt[0] and xyz.y < float(ds_gt[3]):
                try: 
                    xpos, ypos = utils._geo2pixel(xyz.x, xyz.y, ds_gt, 'pixel')
                    g = tgrid[ypos, xpos]
                except: g = ds_nd
                
                if g != ds_nd:
                    d = xyz.z - g
                    s = math.fabs(d / (xyz.z + (g+0.00000001)))

                    if s < max_diff:
                        xyz.z = d
                        yield(xyz)
            ds = None
        
    def query_dump(self, dst_port=sys.stdout, encode=False,  max_diff=.25, **kwargs):
        for xyz in self.yield_diff(self.dem, max_diff):
            xyz.dump(
                include_w = True if self.weights is not None else False,
                dst_port=dst_port,
                encode=encode,
                **kwargs
            )
        
    def run(self):
        dem_infos = demfun.infos(self.dem)
        dem_region = regions.Region().from_geo_transform(geo_transform=dem_infos['geoT'], x_count=dem_infos['nx'], y_count=dem_infos['ny'])

        if not regions.regions_intersect_p(self.region, dem_region):
            utils.echo_error_msg('input region does not intersect with input DEM')


        # grid the difference to array using query_dump / num
        # polygonize the differences and add small buffer (1% or so)
        # make zero array, inverse clipped to buffered polygonized diffs
        # surface zero array and diffs...
        # add surfaced diffs to self.dem
                        
        # diff_cmd = 'gmt blockmedian {region} -I{xinc}/{yinc} | gmt surface {region} -I{xinc}/{yinc} -G_diff.tif=gd+n-9999:GTiff -T.1 -Z1.2 -V -rp -C.5 -Lud -Lld -M{radius}'.format(
        #     region=dem_region, xinc=dem_infos['geoT'][1], yinc=-1*dem_infos['geoT'][5], radius=self.radius
        # )
        diff_cmd = 'gmt blockmedian {region} -I{xinc}/{yinc} | gmt surface {region} -I{xinc}/{yinc} -G_diff.tif=gd+n{ndv}:GTiff -T.1 -Z1.2 -V -Lud -Lld'.format(
            region=self.region.format('gmt'), xinc=self.xinc, yinc=self.yinc, ndv=self.ndv, radius=self.radius
        )

        out, status = utils.run_cmd(
            diff_cmd,
            verbose=self.verbose,
            data_fun=lambda p: self.query_dump(
                dst_port=p, encode=True, max_diff=self.max_diff
            )
        )

        utils.echo_msg('smoothing diff grid...')
        smooth_dem, status = demfun.blur('_diff.tif', '_tmp_smooth.tif', 5)
        #out, status = demfun.grdfilter('_diff.tif', '_tmp_smooth.tif', dist=self.radius, node='pixel', verbose=self.verbose)

        if self.xinc != dem_infos['geoT']:
            utils.echo_msg('resampling diff grid...')
            if demfun.sample_warp('_diff.tif', '__tmp_sample.tif', dem_infos['geoT'][1], -1*dem_infos['geoT'][5],
                             src_region=dem_region)[1] == 0:
                os.rename('__tmp_sample.tif', '_tmp_smooth.tif')
            else:
                utils.echo_warning_msg('failed to resample diff grid')

        #if self.xinc != dem_infos['geoT']:
        utils.echo_msg('resampling diff grid...')
        diff_ds = demfun.sample_warp(
            '_diff.tif',
            None,
            dem_infos['geoT'][1],
            -1*dem_infos['geoT'][5],
            src_region=dem_region
        )[0]
        
        #utils.remove_glob('{}.tif'.format(self.name))

        utils.echo_msg('applying diff grid to dem')
        #diff_ds = gdal.Open('_tmp_smooth.tif')
        diff_band = diff_ds.GetRasterBand(1)
        diff_arr = diff_band.ReadAsArray()
        diff_arr[diff_arr == dem_infos['ndv']] = 0
        diff_arr[diff_arr == self.ndv] = 0
        
        dem_ds = gdal.Open(self.dem)
        dem_band = dem_ds.GetRasterBand(1)
        dem_arr = dem_band.ReadAsArray()

        update_dem_arr = diff_arr + dem_arr
        utils.gdal_write(update_dem_arr, '{}.tif'.format(self.name), dem_infos)
        
        diff_ds = dem_ds = None
        #utils.remove_glob('_tmp_smooth.tif', '_diff.tif')
        return(self)
    
class WaffleFactory():
    """Find and generate a WAFFLE object for DEM generation."""
    
    _modules = {
        'surface': {
            'name': 'surface',
            'datalist-p': True,
            'class':GMTSurface,
        },
        'triangulate': {
            'name': 'triangulate',
            'datalist-p': True,
            'class': GMTTriangulate,
        },
        'nearneighbor': {
            'name': 'nearneihbor',
            'datalist-p': True,
            'class': GMTNearNeighbor,
        },
        'num': {
            'name': 'num',
            'datalist-p': True,
            'class': WafflesNum,
        },
        'linear': {
            'name': 'linear',
            'datalist-p': True,
            'class': WafflesLinear,
        },
        'nearest': {
            'name': 'nearest',
            'datalist-p': True,
            'class': WafflesNearest,
        },
        'average': {
            'name': 'average',
            'datalist-p': True,
            'class': WafflesMovingAverage,
        },
        'invdst': {
            'name': 'invdst',
            'datalist-p': True,
            'class': WafflesInvDst,
        },
        'IDW': {
            'name': 'IDW',
            'datalist-p': True,
            'class': WafflesIDW,
        },
        'vdatum': {
            'name': 'vdatum',
            'datalist-p': False,
            'class': WafflesVDatum,
        },
        'mbgrid': {
            'name': 'mbgrid',
            'datalist-p': True,
            'class': WafflesMBGrid,
        },
        'coastline': {
            'name': 'coastline',
            'datalist-p': False,
            'class': WafflesCoastline,
        },
        'cudem': {
            'name': 'cudem',
            'datalist-p': True,
            'class': WafflesCUDEM,
        },
        'patch': {
            'name': 'patch',
            'datalist-p': True,
            'class': WafflesPatch,
        },
        'lakes': {
            'name': 'lakes',
            'datalist-p': False,
            'class': WafflesLakes,
        },
        'stacks': {
            'name': 'stacks',
            'datalist-p': True,
            'class': WafflesStacks,
        },
        'scipy': {
            'name': 'scipy',
            'datalist-p': True,
            'class': WafflesSciPy,
        },
    }

    def __init__(
            self,
            mod=None,
            data=[],
            src_region=None,
            inc=None,
            xinc=None,
            yinc=None,
            name='waffles_dem',
            node='pixel',
            fmt='GTiff',
            extend=0,
            extend_proc=0,
            weights=None,
            fltr=[],
            sample='bilinear',
            xsample=None,
            ysample=None,
            clip=None,
            chunk=None,
            dst_srs=None,
            srs_transform=False,
            verbose=False,
            archive=False,
            mask=False,
            spat=False,
            clobber=True,
            ndv=-9999,
            block=False,
            cache_dir=waffles_cache
    ):
        self.mod = mod
        self.mod_args = {}
        self.data = data
        self.region = src_region
        self.inc = inc
        self.xinc = xinc
        self.yinc = yinc
        self.sample = sample
        self.xsample = xsample
        self.ysample = ysample
        self.name = name
        self.node = node
        self.fmt = fmt
        self.extend = extend
        self.extend_proc = extend_proc
        self.weights = weights
        self.fltr = fltr
        self.clip = clip
        self.chunk = chunk
        self.dst_srs = dst_srs
        self.srs_transform = srs_transform
        self.archive = archive
        self.mask = mask
        self.spat = spat
        self.clobber = clobber
        self.verbose = verbose
        self.ndv = ndv
        self.block = block
        self.cache_dir = cache_dir

        if self.mod is not None:
            self._parse_mod(self.mod)

        #self._set_config()
        
    def _parse_mod(self, mod):
        opts = mod.split(':')
        if opts[0] in self._modules.keys():
            if len(opts) > 1:
                self.mod_args = utils.args2dict(list(opts[1:]), {})
            self.mod_name = opts[0]
        else:
            utils.echo_error_msg(
                'invalid module name `{}`'.format(opts[0])
            )
            return(None)
        return(self.mod_name, self.mod_args)
        
    def _export_config(self, parse_data=False):
        """export the waffles config info as a dictionary"""

        def _init_data():
            data = [dlim.DatasetFactory(
                fn=" ".join(['-' if x == "" else x for x in dl.split(",")]),
                src_region=self.region.copy().buffer(pct=.25),
                verbose=self.verbose,
                dst_srs=self.dst_srs,
                weight=self.weights,
                cache_dir=self.cache_dir
            ).acquire() for dl in self.data]
            
            return([d for d in data if d is not None])
        
        if parse_data:
            _data = _init_data()
            _data = ["{},{},{}".format(
                os.path.abspath(i.fn) if not i.remote else i.fn, i.data_format, i.weight
            ) for s in [[x for x in d.parse()] for d in _data] for i in s]
        else:
            _data = self.data

        self._config = {
            'mod': self.mod,
            'data': _data,
            'src_region': self.region.export_as_list(
                include_z=True, include_w=True
            ) if self.region is not None else None,
            'xinc': self.xinc,
            'yinc': self.yinc,
            'sample': self.sample,
            'xsample': self.xsample,
            'ysample': self.ysample,
            'name': self.name,
            'node': self.node,
            'fmt': self.fmt,
            'extend': self.extend,
            'extend_proc': self.extend_proc,
            'weights': self.weights,
            'fltr': self.fltr,
            'clip': self.clip,
            'chunk': self.chunk,
            'dst_srs': self.dst_srs,
            'srs_transform': self.srs_transform,
            'verbose': self.verbose,
            'archive': self.archive,
            'mask': self.mask,
            'spat': self.spat,
            'clobber': self.clobber,
            'ndv': self.ndv,
            'block': self.block,
            'cache_dir': self.cache_dir
        }
        return(self._config)

    def acquire(self):
        try:
            return(
                self._modules[self.mod_name]['class'](
                    data=self.data,
                    src_region=self.region,
                    inc=self.inc,
                    xinc=self.xinc,
                    yinc=self.yinc,
                    xsample=self.xsample,
                    ysample=self.ysample,
                    name=self.name,
                    node=self.node,
                    fmt=self.fmt,
                    extend=self.extend,
                    extend_proc=self.extend_proc,
                    weights=self.weights,
                    fltr=self.fltr,
                    sample=self.sample,
                    clip=self.clip,
                    chunk=self.chunk,
                    dst_srs=self.dst_srs,
                    srs_transform=self.srs_transform,
                    archive=self.archive,
                    mask=self.mask,
                    spat=self.spat,
                    clobber=self.clobber,
                    verbose=self.verbose,
                    cache_dir=self.cache_dir,
                    ndv=self.ndv,
                    block=self.block,
                    **self.mod_args
                )
            )
        except Exception as e:
           utils.echo_error_msg(e)
           return(None)
    
    def acquire_from_config(self, config):
        def waffles_config(wc):
            print(wc)
            if wc['src_region'] is not None:
                wc['src_region'] = regions.Region().from_list(
                    wc['src_region']
                )
            return(wc)
        
        mod_name = list(config.keys())[0]
        args = waffles_config(config[mod_name])
        return(self._modules[mod_name]['class'](args))
    
## ==============================================
## Command-line Interface (CLI)
## $ waffles
##
## waffles cli
## ==============================================
waffles_cli_usage = """{cmd} ({wf_version}): Generate DEMs and derivatives.

usage: {cmd} [OPTIONS] DATALIST

Options:
  -R, --region\t\tSpecifies the desired REGION;
\t\t\tWhere a REGION is xmin/xmax/ymin/ymax[/zmin/zmax[/wmin/wmax]]
\t\t\tUse '-' to indicate no bounding range; e.g. -R -/-/-/-/-10/10/1/-
\t\t\tOR an OGR-compatible vector file with regional polygons. 
\t\t\tWhere the REGION is /path/to/vector[:zmin/zmax[/wmin/wmax]].
\t\t\tIf a vector file is supplied, will use each region found therein.
  -E, --increment\tGridding INCREMENT and RESAMPLE-INCREMENT in native units.
\t\t\tWhere INCREMENT is x-inc[/y-inc][:sample-x-inc/sample-y-inc]
  -S, --sample_alg\tReSAMPLE algorithm to use (from gdalwarp)
\t\t\tSet as 'auto' to use 'average' when down-sampling and 'bilinear' when up-sampling
\t\t\tThis switch controls resampling of input raster datasets as well as resampling
\t\t\tthe final DEM if RESAMPLE-INCREMENT is set in -E
  -F, --format\t\tOutput grid FORMAT. [GTiff]
  -M, --module\t\tDesired Waffles MODULE and options. (see available Modules below)
\t\t\tWhere MODULE is module[:mod_opt=mod_val[:mod_opt1=mod_val1[:...]]]
  -O, --output-name\tBASENAME for all outputs.
  -P, --t_srs\t\tProjection of REGION and output DEM.
  -X, --extend\t\tNumber of cells with which to EXTEND the output DEM REGION and a 
\t\t\tpercentage to extend the processing REGION.
\t\t\tWhere EXTEND is dem-extend(cell-count)[:processing-extend(percentage)]
\t\t\te.g. -X6:10 to extend the DEM REGION by 6 cells and the processing region by 10 
\t\t\tpercent of the input REGION.
  -T, --filter\t\tFILTER the output DEM using one or multiple filters. 
\t\t\tWhere FILTER is fltr_id[:fltr_val[:split_value=z]]
\t\t\tAvailable FILTERS:
\t\t\t1: perform a Gaussian Filter at -T1:<factor>.
\t\t\t2: use a Cosine Arch Filter at -T2:<dist(km)> search distance.
\t\t\t3: perform an Outlier Filter at -T3:<aggression>.
\t\t\tThe -T switch may be set multiple times to perform multiple filters.
\t\t\tAppend :split_value=<num> to only filter values below z-value <num>.
\t\t\te.g. -T1:10:split_value=0 to smooth bathymetry (z<0) using Gaussian filter
  -C, --clip\t\tCLIP the output to the clip polygon -C<clip_ply.shp:invert=False>
  -K, --chunk\t\tGenerate the DEM in CHUNKs
  -G, --wg-config\tA waffles config JSON file. If supplied, will overwrite all other options.
\t\t\tGenerate a waffles_config JSON file using the --config flag.
  -D, --cache-dir\tCACHE Directory for storing temp data.
\t\t\tDefault Cache Directory is ~/.cudem_cache; cache will be cleared after a waffles session
\t\t\tto retain the data, use the --keep-cache flag
  -N, --nodata\t\tNODATA value of output DEM

  -f, --transform\tTransform all data to PROJECTION value set with --t_srs/-P where applicable.
  -p, --prefix\t\tSet BASENAME (-O) to PREFIX (append <RES>_nYYxYY_wXXxXX_<YEAR>v<VERSION> info to output BASENAME).
\t\t\tnote: Set Resolution, Year and Version by setting this to 'res=X:year=XXXX:version=X', 
\t\t\tleave blank for default of <INCREMENT>, <CURRENT_YEAR> and <1>, respectively.
  -r, --grid-node\tUse grid-node registration, default is pixel-node
  -w, --weights\t\tUse weights provided in the datalist to weight overlapping data.

  -a, --archive\t\tARCHIVE the datalist to the given region.
  -k, --keep-cache\tKEEP the cache data intact after run
  -m, --mask\t\tGenerate a data MASK raster.
  -s, --spat\t\tGenerate SPATIAL-METADATA.
  -c, --continue\tDon't clobber existing files.
  -q, --quiet\t\tLower verbosity to a quiet.

  --help\t\tPrint the usage text
  --config\t\tSave the waffles config JSON and major datalist
  --modules\t\tDisplay the module descriptions and usage
  --version\t\tPrint the version information

Datalists and data formats:
  A datalist is a file that contains a number of datalist entries, 
  while an entry is a space-delineated line:
  `path [format weight [name source date type resolution hdatum vdatum url]]`

Supported datalist formats: 
  {dl_formats}

Modules (see waffles --modules <module-name> for more info):
  {modules}

CIRES DEM home page: <http://ciresgroups.colorado.edu/coastalDEM>
""".format(cmd=os.path.basename(sys.argv[0]),
           dl_formats=utils._cudem_module_name_short_desc(dlim.DatasetFactory.data_types),
           modules=utils._cudem_module_short_desc(WaffleFactory._modules),
           wf_version=cudem.__version__)

def waffles_cli(argv = sys.argv):
    """run waffles from command-line

    See `waffles_cli_usage` for full cli options.
    """
    
    dls = []
    i_regions = []
    these_regions = []
    module = None
    wg_user = None
    want_prefix = False
    prefix_args = {}
    want_config = False
    keep_cache = False
    status = 0
    i = 1
    wg = {}
    wg['verbose'] = True
    wg['sample'] = 'bilinear'
    wg['xsample'] = None
    wg['ysample'] = None
    wg['dst_srs'] = None
    wg['srs_transform'] = False
    wg['fltr'] = []
    wg['name'] = 'waffles'
    wg['cache_dir'] = waffles_cache
    wg['ndv'] = -9999
    
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
            xy_inc = incs[0].split('/')
            wg['xinc'] = utils.str2inc(xy_inc[0])
            if len(xy_inc) > 1:
                wg['yinc'] = utils.str2inc(xy_inc[1])
            else:
                wg['yinc'] = utils.str2inc(xy_inc[0])
            if len(incs) > 1:
                xy_samples = incs[1].split('/')
                wg['xsample'] = utils.str2inc(xy_samples[0])
                if len(xy_samples) > 1:
                    wg['ysample'] = utils.str2inc(xy_samples[1])
                else:
                    wg['ysample'] = utils.str2inc(xy_samples[0])
            i = i + 1
        elif arg[:2] == '-E':
            incs = arg[2:].split(':')
            xy_inc = incs[0].split('/')
            wg['xinc'] = utils.str2inc(xy_inc[0])
            if len(xy_inc) > 1:
                wg['yinc'] = utils.str2inc(xy_inc[1])
            else:
                wg['yinc'] = utils.str2inc(xy_inc[0])
            if len(incs) > 1:
                xy_samples = incs[1].split('/')
                wg['xsample'] = utils.str2inc(xy_samples[0])
                if len(xy_samples) > 1:
                    wg['ysample'] = utils.str2inc(xy_samples[1])
                else:
                    wg['ysample'] = utils.str2inc(xy_samples[0])

        elif arg == '--sample_alg' or arg == '-S':
            wg['sample'] = argv[i + 1]
            i += 1
        elif arg[:2] == '-S': wg['sample'] = arg[2:]
                    
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
        elif arg == '--t_srs' or arg == '-P' or arg == '-t_srs':
            wg['dst_srs'] = utils.str_or(argv[i + 1], 'epsg:4326')
            i = i + 1
        elif arg[:2] == '-P': wg['dst_srs'] = utils.str_or(arg[2:], 'epsg:4326')
        ## update cache_dir to default to current utils.cache_dir or if an arg, dont add .cudem_cache!
        elif arg == '--cache-dir' or arg == '-D' or arg == '-cache-dir':
            wg['cache_dir'] = os.path.join(utils.str_or(argv[i + 1], os.path.expanduser('~')), '.cudem_cache')
            i = i + 1
        elif arg[:2] == '-D': wg['cache_dir'] = os.path.join(utils.str_or(argv[i + 1], os.path.expanduser('~')), '.cudem_cache')
        elif arg == '--nodata' or arg == '-N' or arg == '-ndv':
            wg['ndv'] = utils.float_or(argv[i + 1], -9999)
            i = i + 1
        elif arg[:2] == '-D': wg['ndv'] = utils.float_or(argv[i + 1], -9999)
        
        elif arg == '--transform' or arg == '-f' or arg == '-transform':
            wg['srs_transform'] = True
            if wg['dst_srs'] is None:
                wg['dst_srs'] = 'epsg:4326'
        elif arg == '-w' or arg == '--weights':
            if 'weights' not in wg.keys():
                wg['weights'] = 1
            else:
                wg['weights'] += 1
        elif arg == '-p' or arg == '--prefix':
            want_prefix = True
            try:
                prefix_opts = argv[i + 1].split(':')
                prefix_args = utils.args2dict(prefix_opts, prefix_args)
                if len(prefix_args) > 0:
                    i += 1
            except:
                pass

        elif arg == '-k' or arg == '--keep-cache': keep_cache = True
        elif arg == '-t' or arg == '--threads': want_threads = True
        elif arg == '-a' or arg == '--archive': wg['archive'] = True
        elif arg == '-m' or arg == '--mask': wg['mask'] = True
        elif arg == '-s' or arg == '--spat': wg['spat'] = True
        elif arg == '-c' or arg == '--continue': wg['clobber'] = False
        elif arg == '-r' or arg == '--grid-node': wg['node'] = 'grid'

        elif arg == '--quiet' or arg == '-q': wg['verbose'] = False
        elif arg == '--config': want_config = True
        elif arg == '--modules' or arg == '-m':
            utils.echo_modules(WaffleFactory._modules, None if i+1 >= len(argv) else sys.argv[i+1])
            sys.exit(0)
        elif arg == '--help' or arg == '-h':
            sys.stderr.write(waffles_cli_usage)
            sys.exit(0)
        elif arg == '--version' or arg == '-v':
            sys.stdout.write('{}\n'.format(cudem.__version__))
            sys.exit(0)
        elif arg[0] == '-':
            sys.stderr.write(waffles_cli_usage)
            utils.echo_error_msg('{} is not a valid waffles cli switch'.format(arg))
            sys.exit(0)
        else: dls.append(arg)
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

                    this_waffle = WaffleFactory(**wg).acquire()
                    this_waffle.generate()

                    sys.exit(0)
            except Exception as e:
                utils.echo_error_msg(e)
                traceback.print_exc()
                sys.exit(-1)
        else:
            utils.echo_error_msg(
                'specified waffles config file does not exist, {}'.format(wg_user)
            )
            sys.stderr.write(waffles_cli_usage)
            sys.exit(-1)

    ## ==============================================
    ## Otherwise run from cli options...
    ## set the dem module
    ## ==============================================        
    if module is None:
        sys.stderr.write(waffles_cli_usage)
        utils.echo_error_msg(
            '''must specify a waffles -M module.'''
        )
        sys.exit(-1)

    if module.split(':')[0] not in WaffleFactory()._modules.keys():
        utils.echo_error_msg(
            '''{} is not a valid waffles module, available modules are: {}'''.format(
                module.split(':')[0], _waffles_module_short_desc()
            )
        )
        sys.exit(-1)
        
    if WaffleFactory()._modules[module.split(':')[0]]['datalist-p']:
        if len(dls) == 0:
            sys.stderr.write(waffles_cli_usage)
            utils.echo_error_msg('''must specify a datalist/entry, try gmrt or srtm for global data.''')
            sys.exit(-1)

    ## ==============================================
    ## check the increment
    ## ==============================================
    if 'xinc' in wg.keys():
        if wg['xinc'] is None:
            sys.stderr.write(waffles_cli_usage)
            utils.echo_error_msg('''must specify a gridding increment.''')
            sys.exit(-1)
        else:
            if wg['yinc'] is None:
                wg['yinc'] = wg['xinc']
    else:
        sys.stderr.write(waffles_cli_usage)
        utils.echo_error_msg('''must specify a gridding increment.''')
        sys.exit(-1)      
    
    ## ==============================================
    ## set the datalists and names
    ## ==============================================
    wg['data'] = dls    
    these_regions = regions.parse_cli_region(i_regions, wg['verbose'])
    name = wg['name']
    for i, this_region in enumerate(these_regions):
        wg['src_region'] = this_region
        if want_prefix or len(these_regions) > 1:
            wg['name'] = utils.append_fn(
                name, wg['src_region'], wg['xsample'] if wg['xsample'] is not None else wg['xinc'], **prefix_args
            )
        if want_config:
            this_waffle = WaffleFactory(mod=module, **wg)
            this_wg = this_waffle._export_config(parse_data=True)
            utils.echo_msg(json.dumps(this_wg, indent=4, sort_keys=True))
            with open('{}.json'.format(this_wg['name']), 'w') as wg_json:
                utils.echo_msg('generating waffles config file: {}.json'.format(this_wg['name']))
                wg_json.write(json.dumps(this_wg, indent=4, sort_keys=True))
        else:
            this_waffle = WaffleFactory(mod=module, **wg)
            if this_waffle is not None:
                if wg['verbose']:
                    utils.echo_msg('------------------------------------------------ :config')
                    this_wg = this_waffle._export_config(parse_data=False)
                    utils.echo_msg(json.dumps(this_wg, sort_keys=True))
                    utils.echo_msg('------------------------------------------------ :config')
                this_waffle_module = this_waffle.acquire()
                if this_waffle_module is not None:
                    this_waffle_module.generate()
                else:
                    if wg['verbose']:
                        utils.echo_error_msg('could not acquire waffles module {}'.format(module))
                        
            #utils.echo_msg('generated DEM: {} @ {}/{}'.format(wf.fn, wf.))
        if not keep_cache:
            utils.remove_glob(wg['cache_dir'])
### End
