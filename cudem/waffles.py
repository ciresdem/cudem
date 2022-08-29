### waffles.py
##
## Copyright (c) 2010 - 2022 Regents of the University of Colorado
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
## nearest (GDAL)
##
## GMT, GDAL and MB-System are required for full functionality.
##
### Code:

import sys
import os
import math
import json
import time

import numpy as np
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
        self.ps_region = self.ps_region.buffer(x_bv=self.xinc*-.5, y_bv=self.yinc*-.5)

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
            sample_alg=self.sample
        ).acquire() for dl in self.data_]
        self.data = [d for d in self.data if d is not None]

    def _coast_region(self):
        """processing region (extended by self.extend_proc)"""
        
        cr = regions.Region().from_region(self.region)
        xcount, ycount, gt = cr.geo_transform(x_inc=self.xinc, y_inc=self.yinc)
        ewp = xcount * ((self.extend_proc+5) * .01)
        nsp = ycount * ((self.extend_proc+5) * .01)
        x_bv = int((ewp + nsp) / 2)
        y_bv = int((ewp + nsp) / 2)

        return(
            cr.buffer(
                x_bv=self.xinc*x_bv,
                y_bv=self.yinc*y_bv
            )
        )
    
    def _proc_region(self):
        """processing region (extended by percentage self.extend_proc)"""
        
        pr = regions.Region().from_region(self.d_region)
        xcount, ycount, gt = pr.geo_transform(x_inc=self.xinc, y_inc=self.yinc)
        ewp = xcount * (self.extend_proc * .01)
        nsp = ycount * (self.extend_proc * .01)
        x_bv = int((ewp + nsp) / 2)
        y_bv = int((ewp + nsp) / 2)

        return(
            pr.buffer(
                x_bv=self.xinc*x_bv,
                y_bv=self.yinc*y_bv
            )
        )
    
    def _dist_region(self):
        """distribution region (extended by self.extend)."""
        
        dr = regions.Region().from_region(self.region)
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

    ## TODO: move following functions to datasets
    ##
    def _xyz_ds(self, src_xyz):
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
        for this_xyz in self.yield_xyz():#src_xyz:
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

    def _xyz_block_array(self, src_xyz, min_count=None, out_name=None):
        """block the src_xyz data to the mean block value

        Yields:
          list: xyz data for each block with data
        """

        xcount, ycount, dst_gt = self.p_region.geo_transform(
            x_inc=self.xinc, y_inc=self.yinc
        )
        gdt = gdal.GDT_Float32
        z_array = np.zeros((ycount, xcount))
        count_array = np.zeros((ycount, xcount))
        if self.weights:
            weight_array = np.zeros((ycount, xcount))
            
        if self.verbose:
            utils.echo_msg(
                'blocking data to {}/{} grid'.format(ycount, xcount)
            )
        for this_xyz in src_xyz:
            xpos, ypos = utils._geo2pixel(
                this_xyz.x, this_xyz.y, dst_gt
            )
            try:
                z_array[ypos, xpos] += this_xyz.z * this_xyz.w
                count_array[ypos, xpos] += 1
                if self.weights:
                    weight_array[ypos, xpos] += this_xyz.w

            except: pass

        count_array[count_array == 0] = np.nan
        if self.weights:
            weight_array[weight_array == 0] = np.nan
            weight_array = (weight_array/count_array)
            z_array = (z_array/weight_array)/count_array
        else:
            z_array = (z_array/count_array)
            weight_array = np.ones((ycount, xcount))

        if min_count is not None:
            z_array[count_array < min_count] = np.nan
            weight_array[count_array < min_count] = np.nan
        
        if out_name is not None:
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
            utils.gdal_write(z_array, '{}_n.tif'.format(out_name), ds_config, verbose=True)
            utils.gdal_write(weight_array, '{}_w.tif'.format(out_name), ds_config, verbose=True)
            utils.gdal_write(count_array, '{}_c.tif'.format(out_name), ds_config, verbose=True)

        z_array[z_array == self.ndv] = np.nan
        weight_array[weight_array == self.ndv] = np.nan
        count_array[count_array == self.ndv] = np.nan
        return(z_array, weight_array, count_array, dst_gt)
    
    def _xyz_block(self, src_xyz, out_name=None):
        """block the src_xyz data to the mean block value

        Yields:
          list: xyz data for each block with data
        """

        z_array, weight_array, count_array, dst_gt = self._xyz_block_array(src_xyz, out_name=out_name)
        ycount, xcount = z_array.shape

        for y in range(0, ycount):
            for x in range(0, xcount):
                z = z_array[y,x]
                if not np.isnan(z):
                    geo_x, geo_y = utils._pixel2geo(x, y, dst_gt)
                    out_xyz = xyzfun.XYZPoint(
                        x=geo_x, y=geo_y, z=z, w=weight_array[y,x]
                    )
                    yield(out_xyz)

        z_array = weight_array = count_array = None

    def _xyz_block_t(self, src_xyz):
        """block the src_xyz data for fast lookup"""

        xcount, ycount, dst_gt = self.p_region.geo_transform(x_inc=self.xinc, y_inc=self.yinc)
        self.block_t = np.empty((ycount, xcount), dtype=object)
        for y in range(0, ycount):
            for x in range(0, xcount):
                self.block_t[y,x] = []

        if self.verbose:
            utils.echo_msg(
                'blocking data into {} buckets'.format(ycount*xcount)
            )
            
        for this_xyz in src_xyz:
            if regions.xyz_in_region_p(this_xyz, self.p_region):
                xpos, ypos = utils._geo2pixel(
                    this_xyz.x, this_xyz.y, dst_gt
                )
                if xpos < xcount and ypos < ycount:
                    self.block_t[ypos,xpos].append(this_xyz.copy())
        
    def _xyz_mask(self, src_xyz, dst_gdal, dst_format='GTiff'):
        """Create a num grid mask of xyz data. The output grid
        will contain 1 where data exists and 0 where no data exists.

        yields the xyz data
        """

        xcount, ycount, dst_gt = self.region.geo_transform(x_inc=self.xinc, y_inc=self.yinc)
        ptArray = np.zeros((ycount, xcount))
        ds_config = demfun.set_infos(
            xcount, ycount, (xcount*ycount), dst_gt, utils.sr_wkt(self.dst_srs),
            gdal.GDT_Float32, self.ndv, 'GTiff'
        )
        for this_xyz in src_xyz:
            yield(this_xyz)            
            if regions.xyz_in_region_p(this_xyz, self.region):
                xpos, ypos = utils._geo2pixel(
                    this_xyz.x, this_xyz.y, dst_gt
                )
                try:
                    ptArray[ypos, xpos] = 1
                except:
                    pass
                
        out, status = utils.gdal_write(ptArray, dst_gdal, ds_config)    

    def yield_array(self, **kwargs):
        for xdl in self.data:
            for array in xdl.yield_array():
                yield(array)
        
    ## TODO: allow spat-meta and archive at same time...
    def yield_xyz(self, region=None, **kwargs):
        """yields the xyz data"""

        # this_datalist = dlim.init_data(self.data_, self.p_region, None, self.dst_srs if self.srs_transform else None, (self.xinc, self.yinc) if block else (None, None), self.verbose)
        # if this_datalist is not None and this_datalist.valid_p(
        #     fmts=dlim.DatasetFactory.data_types[this_datalist.data_format]['fmts']
        # ):
            
        #     if self.archive:
        #         xyz_yield = this_datalist.archive_xyz()
        #     else:
        #         xyz_yield = this_datalist.xyz_yield      
        #     if self.mask:
        #         xyz_yield = self._xyz_mask(xyz_yield, self.mask_fn)

        #     for xyz in xyz_yield:
        #         yield(xyz)
                
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
                
            if self.mask:
                xyz_yield = self._xyz_mask(xyz_yield, self.mask_fn)
                
            if self.block:
                xyz_yield = self._xyz_block(xyz_yield, out_name=self.block) if utils.str_or(self.block) != 'False' else self._xyz_block(xyz_yield)
                
            for xyz in xyz_yield:
                yield(xyz)

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
        
        for xyz in self.yield_xyz(**kwargs):
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
        
        demfun.set_nodata(fn, nodata=self.ndv, convert_array=True)
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
            
        if self.xsample is not None or self.ysample is not None:
            if demfun.sample_warp(fn, '__tmp_sample.tif', self.xsample, self.ysample,
                             src_region=self.p_region, sample_alg=self.sample,
                             ndv=self.ndv, verbose=self.verbose)[1] == 0:
                os.rename('__tmp_sample.tif', fn)
            
        if self.clip is not None:
            clip_args = {}
            cp = self.clip.split(':')
            clip_args['src_ply'] = cp[0]
            clip_args = utils.args2dict(cp[1:], clip_args)
            if clip_args['src_ply'] == 'coastline':
                self.coast = WaffleFactory(
                    mod='coastline:invert=True',
                    data=[],
                    src_region=self.c_region,
                    xinc=self.xinc,
                    yinc=self.yinc,
                    name='tmp_coast',
                    node=self.node,
                    extend=self.extend+12,
                    weights=self.weights,
                    dst_srs=self.dst_srs,
                    srs_transform=self.srs_transform,
                    clobber=True,
                    verbose=self.verbose,
                ).acquire().generate()
                clip_args['src_ply'] = 'tmp_coast.shp'

            ## if clip vector is empty, will return surface with all upper_value :(
            ## maybe have it remove all data in such cases (if invert is true)
            if demfun.clip(fn, '__tmp_clip__.tif', **clip_args)[1] == 0:
                os.rename('__tmp_clip__.tif', '{}'.format(fn))

        if demfun.cut(fn, self.d_region, '__tmp_cut__.tif', node='grid' if self.mod == 'mbgrid' else 'pixel', mode=None)[1] == 0:
            try:
                os.rename('__tmp_cut__.tif', '{}'.format(fn))
            except Exception as e:
                utils.echo_error_msg('could not cut {}; {}'.format(fn, e))

        ## if set_srs fails, set_metadata will skip first entry...
        demfun.set_srs(fn, self.dst_srs)
        demfun.set_metadata(fn, node=self.node, cudem=True)
        if self.fmt != 'GTiff':
            out_dem = utils.gdal2gdal(fn, dst_fmt=self.fmt)
            if out_dem is not None:
                utils.remove_glob(fn)
                
        return(self)
    
    def generate(self):
        """run and process the WAFFLES module"""
        
        if os.path.exists(self.fn):
            if not self.clobber:
                utils.echo_msg(
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
            xcount, ycount, dst_gt = self.p_region.geo_transform(x_inc=self.xinc, y_inc=self.yinc)
            count = 0
            chunks = []
            for srcwin in utils.yield_srcwin((ycount, xcount), self.chunk):
                count += 1
                this_geo_x_origin, this_geo_y_origin = utils._pixel2geo(srcwin[0], srcwin[1], dst_gt)
                this_geo_x_end, this_geo_y_end = utils._pixel2geo(srcwin[0]+srcwin[2], srcwin[1]+srcwin[3], dst_gt)
                this_gt = [this_geo_x_origin, float(dst_gt[1]), 0.0, this_geo_y_origin, 0.0, float(dst_gt[5])]
                this_region = self.region.copy()
                this_region.from_geo_transform(geo_transform=this_gt, x_count=srcwin[2], y_count=srcwin[3])
                this_waffle = WaffleFactory()._modules[self.mod]['class'](
                    data=self.data_,
                    src_region=this_region,
                    inc=self.inc,
                    xinc=self.xinc,
                    yinc=self.yinc,
                    name='{}_{}'.format(self.name, count),
                    node=self.node,
                    fmt=self.fmt,
                    extend=self.extend+10,
                    extend_proc=self.extend_proc+20,
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
            if self.valid_p():
                return(self._process(filter_=True))
            else:
                return(self)
        else:
            self.run()
            
            if self.mask:
                if os.path.exists(self.mask_fn):
                    self._process(fn=self.mask_fn, filter_=False)
                    
            self.waffled = True
            
            [self._process(fn=x, filter_=False) for x in self.aux_dems]
            
            if self.valid_p():
                return(self._process(filter_=True))
            else:
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
## ==============================================
class GMTSurface(Waffle):
    """Waffles GMT surface module.
    Grid the data using GMT 'surface'.
    Data passes through GMT 'blockmean' using weighted mean value if self.weights is True
    """
    
    def __init__(self, tension=.35, relaxation=1.2, max_radius=None,
                 lower_limit='d', upper_limit='d',
                 breakline=None, convergence=None,
                 **kwargs):
        """generate a DEM with GMT surface"""

        self.mod = 'surface'
        self.mod_args = {
            'tension':tension,
            'relaxation':relaxation,
            'max_radius':max_radius,
            'lower_limit':lower_limit,
            'upper_limit':upper_limit,
            'breakline':breakline,
            'convergence':convergence
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
        self.lower_limit = utils.float_or(lower_limit, 'd')
        self.upper_limit = utils.float_or(upper_limit, 'd')
        self.breakline = breakline
        self.max_radius = max_radius
        out, status = utils.run_cmd(
            'gmt gmtset IO_COL_SEPARATOR = SPACE',
            verbose = False
        )        
        
    def run(self):
        dem_surf_cmd = (
            'gmt blockmean {} -I{:.10f}/{:.10f}{} -V | gmt surface -V {} -I{:.10f}/{:.10f} -G{}.tif=gd+n{}:GTiff -T{} -Z{} -Ll{} -Lu{}{}{}{}'.format(
                self.ps_region.format('gmt'),
                self.xinc,
                self.yinc,
                ' -W' if self.weights else '',
                self.ps_region.format('gmt'),
                self.xinc,
                self.yinc,
                self.name,
                self.ndv,
                self.tension,
                self.relaxation,
                self.lower_limit,
                self.upper_limit,
                ' -D{}'.format(self.breakline) if self.breakline is not None else '',
                ' -M{}'.format(self.max_radius) if self.max_radius is not None else '',
                ' -C{}'.format(self.convergence) if self.convergence is not None else ''
            )
        )
        
        out, status = utils.run_cmd(
            dem_surf_cmd,
            verbose=self.verbose,
            data_fun=lambda p: self.dump_xyz(
                dst_port=p, encode=True
            )
        )
        
        return(self)

## ==============================================
## GMT Triangulate
## ==============================================
class GMTTriangulate(Waffle):
    """Waffles GMT triangulate module
    Grid the data using GMT 'triangulate'.
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
        dem_tri_cmd = 'gmt blockmean {} -I{:.10f}/{:.10f}{} -V | gmt triangulate -V {} -I{:.10f}/{:.10f} -G{}.tif=gd:GTiff'.format(
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
## ==============================================
class GMTNearNeighbor(Waffle):
    """Waffles GMT nearneighbor module
    Grid the data using GMT 'nearneighbor'.
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
        dem_nn_cmd = 'gmt blockmean {} -I{:.10f}/{:.10f}{} -V | gmt nearneighbor -V {} -I{:.10f}/{:.10f} -G{}.tif=gd+n{}:GTiff{}{}{}'.format(
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
    """Waffles MBGrid module
    Does not use internal datalists processing,
    input datalist must be compatible with MB-SYSTEM.
    Use dlim --archive to convert a dlim datalist to MB-SYSTEM
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
        # if use_datalists:
        #     #datalist_archive(wg, arch_dir = '.mb_tmp_datalist', verbose = True)
        #     archive = wg['archive']
        #     wg['archive'] = True
        #     for xyz in waffles_yield_datalist(wg): pass
        #     wg['datalist'] = datalists.datalist_major(['archive/{}.datalist'.format(wg['name'])])
        #     wg['archive'] = archive

        mb_region = self.p_region.copy()
        mb_region = mb_region.buffer(x_bv=self.xinc*-.5, y_bv=self.yinc*-.5)
        xsize, ysize, gt = self.p_region.geo_transform(x_inc=self.xinc, node='grid')
        mbgrid_cmd = 'mbgrid -I{} {} -D{}/{} -O{} -A2 -F1 -N -C{} -S0 -X0.1 -T{} {}'.format(
            self.data[0].fn,
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
## ==============================================
class WafflesNum(Waffle):
    """Generate an Uninterpolated grid from XYZ data;
    num methods include: mask, mean, num, landmask and any gmt grd2xyz -A option.
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
        
    def _xyz_num(self, src_xyz):
        """Create a GDAL supported grid from xyz data """

        xcount, ycount, dst_gt = self.p_region.geo_transform(x_inc=self.xinc, y_inc=self.yinc)
        if self.verbose:
            progress = utils.CliProgress(
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
        
        for this_xyz in src_xyz:
            if regions.xyz_in_region_p(this_xyz, self.p_region):
                xpos, ypos = utils._geo2pixel(
                    this_xyz.x, this_xyz.y, dst_gt
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
        if self. verbose:
            progress.end(
                0,
                'generated uninterpolated num grid `{}` @ {}/{}'.format(
                    self.mode, ycount, xcount
                )
            )
            
        utils.gdal_write(out_array, self.fn, ds_config)
        
    def run(self):
        if self.mode.startswith('A'):
            if self.gc['GMT'] is None:
                utils.echo_error_msg(
                    'GMT must be installed to use the Mode `A` with the NUM module'
                )
                return(None, -1)

            out, status = utils.run_cmd('gmt gmtset IO_COL_SEPARATOR = SPACE', verbose = False)            
            dem_xyz2grd_cmd = 'gmt xyz2grd -{} -V {} -I{:.10f}/{:.10f} -G{}.tif=gd:GTiff'.format(
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
            self._xyz_block_array(self.yield_xyz(), min_count=self.min_count, out_name=self.name)
            if self.mode != 'm':
                utils.remove_glob('{}_n.tif'.format(self.name))
            else:
                os.rename('{}_n.tif'.format(self.name), '{}.tif'.format(self.name))
                
            if self.mode != 'n':
                utils.remove_glob('{}_c.tif'.format(self.name))

            if not self.weights:
                utils.remove_glob('{}_w.tif'.format(self.name))
            
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
    """Inverse Distance Weighted."""
    
    def __init__(self, power=1, min_points=8, upper_limit=None, lower_limit=None, radius=None, **kwargs):
        self.mod = 'IDW'
        self.mod_args = {
            'power':power,
            'min_points':min_points,
            'upper_limit':upper_limit,
            'lower_limit':lower_limit,
            'radius':radius
        }
        super().__init__(**kwargs)
        self.power = utils.float_or(power)
        self.min_points = utils.int_or(min_points)
        self.radius = np.inf if radius is None else utils.str2inc(radius) 
        self.upper_limit = utils.float_or(upper_limit)
        self.lower_limit = utils.float_or(lower_limit)

    def run(self):

        ## load and block data into rasters
        ## yield_srcwin through the dem generation
        ## combine the srcwin DEMs
        
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
        outArray = np.empty((ycount, xcount))
        outArray[:] = np.nan
        if self.verbose:
            if self.min_points:
                progress = utils.CliProgress(
                    'generating IDW grid @ {}/{} looking for at least {} neighbors'.format(
                        ycount, xcount, self.min_points
                    )
                )
            else:
                progress = utils.CliProgress(
                    'generating IDW grid @ {}/{}'.format(ycount, xcount)
                )
            i=0
            
        x, y, z, w = [], [], [], []
        for xyz in self.yield_xyz():
            x.append(xyz.x)
            y.append(xyz.y)
            z.append(xyz.z)
            w.append(xyz.w)

        x, y, z, w = np.array(x), np.array(y), np.array(z), np.array(w)        
        obs = np.vstack((x, y)).T
        x = y = None
        xi = np.linspace(self.p_region.xmin, self.p_region.xmax, xcount+1)
        yi = np.linspace(self.p_region.ymin, self.p_region.ymax, ycount+1)
        xi, yi = np.meshgrid(xi, yi)
        xi, yi = xi.flatten(), yi.flatten()
        ask = np.vstack((xi, yi)).T
        #xi = yi = None
        # if no data break now...
        if len(obs) == 0:
            if self.verbose:
                progress.end(
                    -1,
                    'generated IDW grid @ {}/{}'.format(
                        ycount, xcount
                    )
                )
            return(self)

        utils.echo_msg('loaded {} observations for processing...'.format(len(obs)))
        invdisttree = Invdisttree(obs, z, leafsize=10, stat=1)
        interpol = invdisttree(
            ask,
            nnear=self.min_points,
            eps=.1,
            p=self.power,
            dub=self.radius,
            weights=w if self.weights else None
        )
        z = obs = ask = None
        
        if self.upper_limit is not None:
            interpol[interpol > self.upper_limit] = self.upper_limit

        if self.lower_limit is not None:
            interpol[interpol < self.lower_limit] = self.lower_limit
        
        ## fill the grid, interpol.reshape((ycount, xcount)) is flipped...
        ## use 'grid' option in _geo2pixel (nodes from linespace are
        ## at ul corner of cells
        for n, this_z in enumerate(interpol):
            xpos, ypos = utils._geo2pixel(
                xi[n], yi[n], dst_gt, 'grid'
            )
            try:
                outArray[ypos, xpos] = this_z
            except: pass

        if self.verbose:
            progress.end(
                0,
                'generated IDW grid @ {}/{}'.format(
                    ycount, xcount
                )
            )
            
        outArray[np.isnan(outArray)] = self.ndv
        out, status = utils.gdal_write(
            outArray, '{}.tif'.format(self.name), ds_config
        )        
        return(self)    

## ==============================================
## old IDW class...slow, but low mem!
## use WafflesIDW instead.
## ==============================================
class WafflesUIDW(Waffle):
    """Uncertainty Weighted Inverse Distance Weighted.
    
    see: https://ir.library.oregonstate.edu/concern/graduate_projects/79407x932
    """
    
    def __init__(
            self, radius=None, power=2, min_points=None, **kwargs
    ):
        self.mod = 'IDW'
        self.mod_args = {
            'radius':radius,
            'power':power,
            'min_points':min_points
        }
        super().__init__(**kwargs)
        if radius is not None:
            self.radius = utils.str2inc(radius)
        else:
            self.radius = self.xinc
        
        self.power = utils.float_or(power)
        self.min_points = utils.int_or(min_points)
        
    def _distance(self, pnt0, pnt1):
        return(math.sqrt(sum([(a-b)**2 for a, b in zip(pnt0, pnt1)])))
            
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
        outArray = np.empty((ycount, xcount))
        outArray[:] = np.nan
        if self.verbose:
            if self.min_points:
                progress = utils.CliProgress(
                    'generating IDW grid @ {} and {}/{} looking for at least {} volunteers'.format(
                        self.radius, ycount, xcount, self.min_points
                    )
                )
            else:
                progress = utils.CliProgress(
                    'generating IDW grid @ {} and {}/{}'.format(
                        self.radius, ycount, xcount
                    )
                )
            i=0

        x, y, z, w = [], [], [], []
        for xyz in self.yield_xyz():
            x.append(xyz.x)
            y.append(xyz.y)
            z.append(xyz.z)
            w.append(xyz.w)

        x, y, z, w = np.array(x), np.array(y), np.array(z), np.array(w)        
        obs = np.vstack((x, y)).T
        tree = spatial.cKDTree(obs)
        for y_g in range(0, ycount):
            if self.verbose:
                i+=1
                progress.update_perc((i, ycount))
                
            for x_g in range(0, xcount):
                xg, yg = utils._pixel2geo(x_g, y_g, dst_gt)
                bucket = tree.query_ball_point([xg, yg], self.radius)
                if bucket:
                    zs, ds, ws = [], [], []                    
                    for obs_index in bucket:
                        d = self._distance(obs[obs_index], [xg, yg])
                        zs.append(z[obs_index])
                        ws.append(w[obs_index]**self.power)
                        if d > 0:
                            ds.append(1./d**self.power)
                        else: ds.append(0)

                    distances = np.transpose(ds)
                    weights = np.transpose(ws)
                    out_weights = distances*weights if self.weights else distances
                    sums = sum(out_weights)
                    if sums > 0:
                        outArray[y_g, x_g] = np.dot(zs, out_weights)/sums
                    else:
                        outArray[y_g, x_g] = sum(zs)/len(zs)
                            
        if self.verbose:
            progress.end(
                0,
                'generated IDW grid {}/{}'.format(
                    ycount, xcount
                )
            )

        ds = None            
        outArray[np.isnan(outArray)] = self.ndv
        out, status = utils.gdal_write(
            outArray, '{}.tif'.format(self.name), ds_config
        )        
        return(self)    

## ==============================================
## Waffles VDatum
## ==============================================
class WafflesVDatum(Waffle):
    """vertical datum transformation grid"""

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
            self.p_region, self.xinc, self.yinc, self.vdatum_in, self.vdatum_out, cache_dir=waffles_cache
        ).run(outfile='{}.tif'.format(self.name))
        
        return(self)

## ==============================================
## GDAL gridding
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
        
    def run(self):
        xcount, ycount, dst_gt = self.p_region.geo_transform(x_inc=self.xinc, y_inc=self.yinc)
        _prog = utils.CliProgress(
            'running GDAL GRID {} algorithm @ {} and {}/{}...'.format(
                self.alg_str.split(':')[0], self.p_region.format('fn'), xcount, ycount
            )
        )
        _prog_update = lambda x, y, z: _prog.update()
        ds = self._xyz_ds(self.yield_xyz())
        #ds = vectorize_xyz()
        if ds.GetLayer().GetFeatureCount() == 0:
            utils.echo_error_msg('no input data')
            
        gd_opts = gdal.GridOptions(
            outputType = gdal.GDT_Float32,
            noData = self.ndv,
            format = 'GTiff',
            width = xcount,
            height = ycount,
            algorithm = self.alg_str,
            callback = _prog_update if self.verbose else None,
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
        _prog.end(
            0,
            'ran GDAL GRID {} algorithm @ {}.'.format(
                self.alg_str.split(':')[0], self.p_region.format('fn')
            )
        )        
        ds = None
        return(self)

class WafflesLinear(WafflesGDALGrid):
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
## ...and some tests
## ==============================================
class WafflesCUDEM(Waffle):
    """Waffles CUDEM gridding module

    generate a DEM with `pre_surface`s which are generated
    at lower resolution and with various weight threshholds.

To generate a typical CUDEM tile, generate 1 pre-surface ('bathy_surface'), clipped to a coastline.
Use a min_weight that excludes low-resolution bathymetry data from being used as input in the final
DEM generation. 
    """
    
    def __init__(
            self,
            min_weight=None,
            pre_count=1,
            smoothing=None,
            landmask=False,
            poly_count=3,
            keep_auxilary=False,
            **kwargs
    ):
        self.mod = 'cudem'
        self.mod_args = {
            'min_weight':min_weight,
            'pre_count':pre_count,
            'smoothing':smoothing,
            'landmask':landmask,
            'poly_count':poly_count,
            'keep_auxilary':keep_auxilary
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
        self.keep_auxilary = keep_auxilary
        
    def run(self):
        pre = self.pre_count
        self.p_region.buffer(pct=self.pre_count)
        pre_weight = 0
        pre_region = self.p_region.copy()
        pre_region.wmin = None
        pre_clip = None 
        upper_limit = None
        coast = '{}_cst'.format(self.name)

        utils.echo_msg('target region: {}'.format(self.region))
        utils.echo_msg('processing region: {}'.format(self.p_region))
        
        self._xyz_block_array(self.yield_xyz(), out_name=self.name)
        n = '{}_n.tif'.format(self.name)
        w = '{}_w.tif'.format(self.name)
        c = '{}_c.tif'.format(self.name)
        if self.min_weight is None:
            self.min_weight = demfun.percentile(w, 75)
            
        utils.echo_msg('cudem min weight is: {}'.format(self.min_weight))
        pre_data = ['{},200:weight_mask={},1'.format(n, w)]
        
        if self.landmask:
            upper_limit = -0.1
            pre_region.zmax = 1
            if not os.path.exists(utils.str_or(self.landmask)):
                if os.path.exists('{}.shp'.format(utils.append_fn(self.landmask, self.region, self.xinc))):
                    coastline = '{}.shp'.format(utils.append_fn(self.landmask, self.region, self.xinc))
                else:
                    self.coast = WaffleFactory(
                        mod='coastline:polygonize={}'.format(self.poly_count),
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
                        verbose=self.verbose,
                    ).acquire().generate()
                    coastline = '{}.shp'.format(self.coast.name)
            else:
                coastline = self.landmask
            pre_clip = '{}:invert=True'.format(coastline)
        
        while pre >= 0:
            pre_xinc = self.xinc * (3**pre)
            pre_yinc = self.yinc * (3**pre)
            xsample = self.xinc * (3**(pre - 1))
            ysample = self.yinc * (3**(pre - 1))
            if xsample == 0: xsample = self.xinc
            if ysample == 0: ysample = self.yinc
            
            pre_filter=['1:{}'.format(self.smoothing)] if self.smoothing is not None else []
            if pre != self.pre_count:
                pre_weight = self.min_weight/(pre + 1) if pre > 0 else self.min_weight
                #pre_weight = self.min_weight/(3**pre) if pre > 0 else self.min_weight
                if pre_weight == 0: pre_weight = 1-e20
                pre_data = [
                    '{},200:weight_mask={},1'.format(n, w),
                    '{}.tif,200,{}'.format(
                        utils.append_fn('_pre_surface', pre_region, pre+1), pre_weight
                    )
                ]
                pre_region.wmin = pre_weight

            pre_region = self._proc_region()
            pre_region.buffer(pct=pre)

            if pre == 0:
                pre_region.zmax = None

            if pre != self.pre_count:
                pre_region.wmin = pre_weight
                
            pre_name = utils.append_fn('_pre_surface', pre_region, pre)
            utils.echo_msg('pre region: {}'.format(pre_region))
            
            pre_surface = WaffleFactory(
                mod='surface:tension=1:upper_limit={}'.format(upper_limit if pre !=0 else 'd') if pre !=0 else 'surface:tension=1',
                data=pre_data,
                src_region=pre_region,
                xinc=pre_xinc if pre !=0 else self.xinc,
                yinc=pre_yinc if pre != 0 else self.yinc,
                name=pre_name if pre !=0 else self.name,
                node=self.node,
                fltr=pre_filter if pre !=0 else [],
                weights=1,
                dst_srs=self.dst_srs,
                srs_transform=self.srs_transform,
                clobber=True,
                sample=self.sample,
                xsample=xsample if pre !=0 else None,
                ysample=ysample if pre !=0 else None,
                verbose=self.verbose,
                clip=pre_clip if pre !=0 else None,
            ).acquire().generate()                
            pre -= 1
            
        if not self.keep_auxilary:
            utils.remove_glob('*_pre_surface*')
            utils.remove_glob('{}*'.format(n), '{}*'.format(w), '{}*'.format(c), '{}.*'.format(coast))
            
        return(self)
    
class WafflesCUDEM_stacks(Waffle):
    """Waffles CUDEM gridding module

    generate a DEM with `pre_surface`s which are generated
    at lower resolution and with various weight threshholds.
    """
    
    def __init__(
            self,
            min_weight=None,
            pre_count=1,
            smoothing=None,
            landmask=False,
            poly_count=3,
            keep_auxilary=False,
            **kwargs
    ):
        self.mod = 'cudem'
        self.mod_args = {
            'min_weight':min_weight,
            'pre_count':pre_count,
            'smoothing':smoothing,
            'landmask':landmask,
            'poly_count':poly_count,
            'keep_auxilary':keep_auxilary
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
        self.keep_auxilary = keep_auxilary
        
    def run(self):
        pre = self.pre_count
        pre_weight = 0
        pre_region = self.p_region.copy()
        pre_region.wmin = None
        pre_clip = None 
        upper_limit = None
        coast = '{}_cst'.format(self.name)

        pre_stack = WaffleFactory(
            mod='stacks:supercede=True:keep_weights=True',
            data=self.data_,
            src_region=pre_region,
            xinc=self.xinc,
            yinc=self.yinc,
            name='{}_n'.format(self.name),
            node=self.node,
            weights=1,
            dst_srs=self.dst_srs,
            srs_transform=self.srs_transform,
            clobber=True,
            verbose=self.verbose,
        ).acquire().generate()
        
        n = pre_stack.fn
        w = '{}_w.tif'.format('.'.join(pre_stack.fn.split('.')[:-1]))

        if self.min_weight is None:
            self.min_weight = demfun.percentile(w, 75)
            
        utils.echo_msg('cudem min weight is: {}'.format(self.min_weight))
        pre_data = ['{},200:weight_mask={},1'.format(n, w)]
        
        if self.landmask:
            upper_limit = -0.1
            pre_region.zmax = 1
            if not os.path.exists(utils.str_or(self.landmask)):
                if os.path.exists('{}.shp'.format(utils.append_fn(self.landmask, self.region, self.xinc))):
                    coastline = '{}.shp'.format(utils.append_fn(self.landmask, self.region, self.xinc))
                else:
                    self.coast = WaffleFactory(
                        mod='coastline:polygonize={}'.format(self.poly_count),
                        data=pre_data,
                        src_region=pre_region,
                        xinc=self.xinc,
                        yinc=self.yinc,
                        name=coast,
                        node=self.node,
                        weights=self.weights,
                        dst_srs=self.dst_srs,
                        srs_transform=self.srs_transform,
                        clobber=True,
                        verbose=self.verbose,
                    ).acquire().generate()
                    coastline = '{}.shp'.format(self.coast.name)
            else:
                coastline = self.landmask
            pre_clip = '{}:invert=True'.format(coastline)
        
        while pre >= 0:
            pre_xinc = self.xinc * (3**pre)
            pre_yinc = self.yinc * (3**pre)
            xsample = self.xinc * (3**(pre - 1))
            ysample = self.yinc * (3**(pre - 1))
            if xsample == 0: xsample = self.xinc
            if ysample == 0: ysample = self.yinc

            pre_name = utils.append_fn('_pre_surface', pre_region, pre)
            pre_filter=['1:{}'.format(self.smoothing)] if self.smoothing is not None else []
            if pre != self.pre_count:
                pre_weight = self.min_weight/(pre + 1) if pre > 0 else self.min_weight
                if pre_weight == 0: pre_weight = 1-e20
                pre_data = [
                    '{},200:weight_mask={},1'.format(n, w),
                    '{}.tif,200,{}'.format(
                        utils.append_fn('_pre_surface', pre_region, pre+1), pre_weight
                    )
                ]
                pre_region.wmin = pre_weight

            if pre == 0:
                pre_region.zmax = None

            pre_surface = WaffleFactory(
                mod='stacks',
                data=pre_data,
                src_region=pre_region,
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
                clip=pre_clip if pre !=0 else None,
            ).acquire().generate()                
            pre -= 1
            
        if not self.keep_auxilary:
            utils.remove_glob('*_pre_surface*')
            utils.remove_glob('{}*'.format(n), '{}*'.format(w), '{}.*'.format(coast))
            
        return(self)
    
## ==============================================
## Waffles Raster Stacking
## ==============================================
class WafflesStacks(Waffle):
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
        
        return(self)
    
## ==============================================
## Waffles Coastline/Landmask
## ==============================================
class WafflesCoastline(Waffle):
    def __init__(
            self,
            want_nhd=True,
            want_copernicus=True,
            want_gmrt=False,
            want_lakes=False,
            want_buildings=False,
            invert=False,
            polygonize=True,
            osm_tries=3,
            **kwargs
    ):
        """Generate a landmask from various sources.
        
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
            'invert':invert,
            'polygonize':polygonize
        }
        super().__init__(**kwargs)

        self.want_nhd = want_nhd
        self.want_gmrt = want_gmrt
        self.want_copernicus = want_copernicus
        self.want_lakes = want_lakes
        self.want_buildings = want_buildings
        self.invert = invert
        self.polygonize = polygonize
        self.osm_tries = utils.int_or(osm_tries, 3)

        self.coast_array = None
        self.ds_config = None

        import cudem.fetches.copernicus
        import cudem.fetches.tnm
        import cudem.fetches.gmrt
        import cudem.fetches.hydrolakes
        import cudem.fetches.osm
        import cudem.fetches.utils

        self.f_region = self.p_region.copy()
        self.f_region.buffer(pct=5)
        self.f_region.src_srs = self.dst_srs
        self.wgs_region = self.f_region.copy()
        if self.dst_srs is not None:
            self.wgs_region.warp('epsg:4326')
        else:
            self.dst_srs = 'epsg:4326'
        
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
            
        if len(self.data) > 0:
            self._load_data()

        utils.echo_msg(
            'finanlizing array for region {} at {} {}...'.format(
                self.p_region.format('gmt'), self.ds_config['nx'], self.ds_config['ny']
            )
        )
        self._finalize_array()
        utils.echo_msg('writing array to {}.tif...'.format(self.name))        
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
            utils.sr_wkt(self.dst_srs),
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
            src_region=self.f_region, weight=self.weights, verbose=self.verbose, layer='topo'
        )
        this_gmrt._outdir = self.cache_dir
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
        gmrt_ds_arr = gmrt_ds.GetRasterBand(1).ReadAsArray()
        gmrt_ds_arr[gmrt_ds_arr > 0] = 1
        gmrt_ds_arr[gmrt_ds_arr < 0] = 0
        self.coast_array += gmrt_ds_arr
        gmrt_ds = gmrt_ds_arr = None
        
    def _load_copernicus(self):
        """copernicus"""

        this_cop = cudem.fetches.copernicus.CopernicusDEM(
            src_region=self.wgs_region, weight=self.weights, verbose=self.verbose, datatype='1'
        )
        this_cop._outdir = self.cache_dir
        this_cop.run()

        fr = cudem.fetches.utils.fetch_results(this_cop, want_proc=False)
        fr.daemon = True
        fr.start()
        fr.join()

        dst_srs = osr.SpatialReference()
        dst_srs.SetFromUserInput(self.dst_srs)
        dst_srs.AutoIdentifyEPSG()
        dst_auth = dst_srs.GetAttrValue('AUTHORITY', 1)
        dst_srs.SetFromUserInput('epsg:{}'.format(dst_auth))
       
        for i, cop_tif in enumerate(this_cop.results):
            demfun.set_nodata(cop_tif[1], 0)            
            cop_ds = demfun.generate_mem_ds(self.ds_config, name='copernicus')
            gdal.Warp(
                cop_ds, cop_tif[1], dstSRS=dst_srs, resampleAlg=self.sample,
                callback=gdal.TermProgress, srcNodata=0
            )
            cop_ds_arr = cop_ds.GetRasterBand(1).ReadAsArray()
            cop_ds_arr[cop_ds_arr != 0] = 1
            self.coast_array += cop_ds_arr
            cop_ds = cop_ds_arr = None
    
    def _load_nhd(self):
        """USGS NHD (HIGH-RES U.S. Only)
        Fetch NHD (NHD High/Plus) data from TNM to fill in near-shore areas. 
        High resoultion data varies by location...
        """

        this_tnm = cudem.fetches.tnm.TheNationalMap(
            src_region=self.wgs_region,
            weight=self.weights,
            verbose=self.verbose,
            where=["Name LIKE '%Hydro%'"],
            extents='HU-8 Subbasin,HU-4 Subregion'
        )
        this_tnm._outdir = self.cache_dir
        this_tnm.run()
        fr = cudem.fetches.utils.fetch_results(this_tnm, want_proc=False)
        fr.daemon = True
        fr.start()
        fr.join()
        if len(this_tnm.results) > 0:
            for i, tnm_zip in enumerate(this_tnm.results):
                tnm_zips = utils.unzip(tnm_zip[1], self.cache_dir)
                gdb = '.'.join(tnm_zip[1].split('.')[:-1]) + '.gdb'
                utils.run_cmd(
                    'ogr2ogr -update -append nhdArea_merge.shp {} NHDArea -where "FType=312 OR FType=336 OR FType=445 OR FType=460 OR FType=537" -clipdst {} 2>/dev/null'.format(
                        gdb, self.p_region.format('ul_lr')
                    ),
                    verbose=True
                )
                utils.run_cmd(
                    'ogr2ogr -update -append nhdArea_merge.shp {} NHDPlusBurnWaterBody -clipdst {} 2>/dev/null'.format(
                        gdb, self.p_region.format('ul_lr')
                    ),
                    verbose=True
                )
                utils.run_cmd(
                    'ogr2ogr -update -append nhdArea_merge.shp {} NHDWaterBody -where "FType=493 OR FType=466" -clipdst {} 2>/dev/null'.format(
                        gdb, self.p_region.format('ul_lr')
                    ),
                    verbose=True
                )
                #utils.remove_glob(tnm_zip[1], *tnm_zips, gdb)
                #utils.remove_glob(*tnm_zips, gdb)

            utils.run_cmd(
                'gdal_rasterize -burn 1 nhdArea_merge.shp nhdArea_merge.tif -te {} -ts {} {} -ot Int32'.format(
                    self.p_region.format('te'),
                    self.ds_config['nx'],
                    self.ds_config['ny'],
                ),
                verbose=True
            )
            tnm_ds = gdal.Open('nhdArea_merge.tif')
            if tnm_ds is not None:
                tnm_ds_arr = tnm_ds.GetRasterBand(1).ReadAsArray()
                tnm_ds_arr[tnm_ds_arr < 1] = 0
                self.coast_array -= tnm_ds_arr
                tnm_ds = tnm_ds_arr = None
                
            utils.remove_glob('nhdArea_merge.*')

    def _load_lakes(self):
        """HydroLakes -- Global Lakes"""
        
        this_lakes = cudem.fetches.hydrolakes.HydroLakes(
            src_region=self.f_region, weight=self.weights, verbose=self.verbose
        )
        this_lakes._outdir = self.cache_dir
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
        lk_ds = ogr.Open(lakes_shp)
        if lk_ds is not None:
            lk_layer = lk_ds.GetLayer()
            lk_layer.SetSpatialFilter(self.f_region.export_as_geom())
            gdal.RasterizeLayer(lakes_ds, [1], lk_layer, burn_values=[-1])
            lakes_ds_arr = lakes_ds.GetRasterBand(1).ReadAsArray()
            self.coast_array[lakes_ds_arr == -1] = 0
            lakes_ds = lk_ds = None
        else:
            utils.echo_error_msg('could not open {}'.format(lakes_shp))

    def _load_bldgs(self):
        """load buildings from OSM
        
        OSM has a size limit, so will chunk the region and
        do multiple osm calls
        """

        bldg_ds = demfun.generate_mem_ds(self.ds_config, name='bldg')
        this_osm = cudem.fetches.osm.OpenStreetMap(
            src_region=self.f_region, weight=self.weights, verbose=self.verbose
        )
        this_osm._outdir = self.cache_dir
        this_osm.run()
        for osm_result in this_osm.results:
            if cudem.fetches.utils.Fetch(osm_result[0], verbose=self.verbose).fetch_file(osm_result[1], check_size=False, tries=self.osm_tries) >= 0:
                osm_ds = ogr.Open(osm_result[1])
                if osm_ds is not None:
                    osm_layer = osm_ds.GetLayer('multipolygons')
                    osm_layer.SetAttributeFilter("building!=''")
                    gdal.RasterizeLayer(bldg_ds, [1], osm_layer, burn_values=[-1])
                    bldg_arr = bldg_ds.GetRasterBand(1).ReadAsArray()
                    self.coast_array[bldg_arr == -1] = 0
                else:
                    utils.echo_error_msg('could not open ogr dataset {}'.format(osm_result[1]))
            #else:
            #/    utils.echo_error_msg('failed to fetch {}'.format(osm_result[0]))
            
        bldg_ds = None
        
    def _load_data(self):
        """load data from user datalist and amend coast_array"""

        ##TODO: use stacks instead
        for this_data in self.data:
            for this_xyz in this_data.yield_xyz():
                xpos, ypos = utils._geo2pixel(
                    this_xyz.x, this_xyz.y, self.ds_config['geoT']
                )
                try:
                    if this_xyz.z >= 0:
                        self.coast_array[ypos, xpos] += 1
                    else:
                        self.coast_array[ypos, xpos] -= 1
                except: pass
        
    def _write_coast_array(self):
        """write coast_array to file"""
        
        utils.gdal_write(
            self.coast_array, '{}.tif'.format(self.name), self.ds_config,
        )

    def _write_coast_poly(self, poly_count=3):
        """convert to coast_array vector"""

        tmp_ds = ogr.GetDriverByName('ESRI Shapefile').CreateDataSource(
            'tmp_c_{}.shp'.format(self.name)
        )
        
        if tmp_ds is not None:
            tmp_layer = tmp_ds.CreateLayer('tmp_c_{}'.format(self.name), None, ogr.wkbMultiPolygon)
            tmp_layer.CreateField(ogr.FieldDefn('DN', ogr.OFTInteger))
            demfun.polygonize('{}.tif'.format(self.name), tmp_layer, verbose=self.verbose)
            tmp_ds = None
            
        # utils.run_cmd(
        #     'ogr2ogr -dialect SQLITE -sql "SELECT * FROM tmp_c_{} WHERE DN=0 order by ST_AREA(geometry) desc limit {}" {}.shp tmp_c_{}.shp'.format(
        #         self.name, poly_count, self.name, self.name),
        #     verbose=True
        # )
        utils.run_cmd(
            'ogr2ogr {}.shp tmp_c_{}.shp'.format(
                self.name, self.name, self.name),
            verbose=True
        )
        
        utils.remove_glob('tmp_c_{}.*'.format(self.name))
        utils.run_cmd(
            'ogrinfo -dialect SQLITE -sql "UPDATE {} SET geometry = ST_MakeValid(geometry)" {}.shp'.format(
                self.name, self.name)
        )

        utils.gdal_prj_file(self.name + '.prj', self.dst_srs)

## ==============================================
## Waffles Bathymetry - testing
## ==============================================
class WafflesLakes(Waffle):
    def __init__(self, apply_elevations=False, min_area=None, max_area=None, min_id=None, max_id=None, depth='globathy', elevations='copernicus', **kwargs):
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
            src_region=self.p_region, weight=self.weights, verbose=self.verbose
        )
        this_lakes._outdir = self.cache_dir
        this_lakes.run()

        fr = cudem.fetches.utils.fetch_results(this_lakes, want_proc=False, check_size=False)
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
        cudem.fetches.utils.Fetch(_globathy_url, verbose=self.verbose).fetch_file(globathy_zip, check_size=False)
        globathy_csvs = utils.unzip(globathy_zip, self.cache_dir)        
        globathy_csv = os.path.join(self.cache_dir, 'GLOBathy_basic_parameters/GLOBathy_basic_parameters(ALL_LAKES).csv')

        _prog = utils.CliProgress('parsing globathy parameters...')
        with open(globathy_csv, mode='r') as globc:
            reader = csv.reader(globc)
            next(reader)
            if len(ids) > 0:
                globd = {}
                for row in reader:
                    if int(row[0]) in ids:
                        _prog.update()
                        globd[int(row[0])] = float(row[-1])
                        ids.remove(int(row[0]))
                    if len(ids) == 0:
                        break
            else:
                globd = {int(row[0]):float(row[-1]) for row in reader}

        _prog.end(0, 'parsed {} globathy parameters'.format(len(globd)))
        return(globd)

    def _fetch_gmrt(self, gmrt_region=None):
        """GMRT - Global low-res.
        """

        if gmrt_region is None:
            gmrt_region = self.p_region
        
        this_gmrt = cudem.fetches.gmrt.GMRT(
            src_region=gmrt_region, weight=self.weights, verbose=self.verbose, layer='topo'
        )
        this_gmrt._outdir = self.cache_dir
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
            src_region=cop_region, weight=self.weights, verbose=self.verbose, datatype='1'
        )
        this_cop._outdir = self.cache_dir
        this_cop.run()

        fr = cudem.fetches.utils.fetch_results(this_cop, want_proc=False, check_size=False)
        fr.daemon = True
        fr.start()
        fr.join()

        dst_srs = osr.SpatialReference()
        dst_srs.SetFromUserInput(self.dst_srs)
        cop_ds = demfun.generate_mem_ds(self.ds_config, name='copernicus')

        _cop_prog = utils.CliProgress('combining Copernicus tiles')
        [gdal.Warp(cop_ds, cop_tif[1], dstSRS=dst_srs, resampleAlg=self.sample) for cop_tif in this_cop.results]
        #for i, cop_tif in enumerate(this_cop.results):
        #    _cop_prog.update_perc((i, len(this_cop.results)))
        #    gdal.Warp(cop_ds, cop_tif[1], dstSRS=dst_srs, resampleAlg=self.sample)
        _cop_prog.end(0, 'combined {} Copernicus tiles'.format(len(this_cop.results)))
        
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
        #utils.echo_msg('set labels')
        nfeatures = np.arange(1, nfeatures +1)
        maxes = ndimage.maximum(shore_distance_arr, labels, nfeatures)
        #utils.echo_msg('set maxes')
        max_dist_arr = np.zeros(np.shape(shore_distance_arr))
        #utils.echo_msg('set max_dist_arr')
        p = utils.CliProgress('applying labesl to array...')
        for n, x in enumerate(nfeatures):
            p.update_perc((n,len(nfeatures)))
            max_dist_arr[labels==x] = maxes[x-1]
        p.end(0, 'applied labels to max_dist_arr')
        
        max_dist_arr[max_dist_arr == 0] = np.nan
        bathy_arr = (shore_distance_arr * lake_depths_arr) / max_dist_arr
        bathy_arr[bathy_arr == 0] = np.nan

        #utils.echo_msg('applied calculation to array')
        
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
        self.p_region.buffer(pct=2)
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
        
        lk_prog = utils.CliProgress('processing {} lakes'.format(lk_features))

        lk_regions = self.p_region.copy()
        for lk_f in lk_layer:
            this_region = regions.Region()
            lk_geom = lk_f.GetGeometryRef()
            lk_wkt = lk_geom.ExportToWkt()
            this_region.from_list(ogr.CreateGeometryFromWkt(lk_wkt).GetEnvelope())
            lk_regions = regions.regions_merge(lk_regions, this_region)

        while not regions.regions_within_ogr_p(self.p_region, lk_regions):
            self.p_region.buffer(pct=2)

        self._init_bathy()
        
        ## fetch and initialize the copernicus data
        if self.elevations == 'copernicus':
            cop_ds = self._fetch_copernicus(cop_region=self.p_region)
            cop_band = cop_ds.GetRasterBand(1)
        elif self.elevations == 'gmrt':
            cop_ds = self._fetch_gmrt(gmrt_region=self.p_region)
            cop_band = cop_ds.GetRasterBand(1)
        
        ## initialize the tmp datasources
        prox_ds = demfun.generate_mem_ds(self.ds_config, name='prox')
        msk_ds = demfun.generate_mem_ds(self.ds_config, name='msk')
        
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
            _p = utils.CliProgress('Assigning Globathy Depths to rasterized lakes...')
            for n, this_id in enumerate(lk_ids):
                #if self.depth == 'globathy':
                depth = globd[this_id]
                # elif self.depth == 'hydrolakes':
                #     lk_layer.SetAttributeFilter('Hylak_id = {}'.format(this_id))
                #     for lk_f in lk_layer:
                #         depth = float(lk_f.GetField('Depth_avg'))
                # elif utils.float_or(self.depth) is not None:
                #     depth = utils.float_or(self.depth)
                # else:
                #     depth = 1

                #msk_arr[msk_arr == this_id] = globd[this_id]
                msk_arr[msk_arr == this_id] = depth
                _p.update_perc((n, len(lk_ids)))
            _p.end(0, 'Assigned {} Globathy Depths to rasterized lakes.'.format(len(lk_ids)))
            
        elif self.depth == 'hydrolakes':
            globd = None
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
            
        lk_ds = None
            
        ## calculate proximity of lake cells to shore
        prox_band = prox_ds.GetRasterBand(1)
        proximity_options = ["VALUES=0", "DISTUNITS=PIXEL"]
        gdal.ComputeProximity(msk_band, prox_band, options=proximity_options, callback=gdal.TermProgress)
        
        #utils.echo_msg('Assigning Globathy Depths to rasterized lakes...')
        #msk_arr[msk_arr[np.in1d(msk_arr, globd.keys())]] = globd.values()
        #print(zip(*globd.keys()))
        #print(list(globd.values()))
        #msk_arr[list(globd.keys()),:] = list(globd.values())
        #msk_arr[zip(*globd.keys())] = list(globd.values())
        #msk_arr[np.isin(msk_arr, list(globd.keys()))] = list(globd.values())
        
        prox_arr = prox_band.ReadAsArray()
        cop_arr = cop_band.ReadAsArray()

        utils.echo_msg('Calculating simulated lake depths...')
        ## apply calculation from globathy
        bathy_arr = self.apply_calculation(
            prox_arr,
            msk_arr,
            shore_arr=cop_arr,
        )

        lk_prog.end(0, 'processed {} lakes'.format(lk_features))
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
        lk_prog = utils.CliProgress('processing {} lakes'.format(lk_features))

        lk_ids = []
        for i, feat in enumerate(lk_layer):
            lk_ids.append(feat.GetField('Hylak_id'))

        globd = self._fetch_globathy(ids=lk_ids)
        
        for i, feat in enumerate(lk_layer):
            lk_prog.update_perc((i, lk_features))
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
            
        lk_prog.end(0, 'processed {} lakes'.format(lk_features))
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
            utils.echo_error_msg('must specify DEM to update (:dem=fn) to run the update module.')
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
                    xpos, ypos = utils._geo2pixel(xyz.x, xyz.y, ds_gt)
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


        ## grid the difference to array using query_dump / num
        ## polygonize the differences and add small buffer (1% or so)
        ## make zero array, inverse clipped to buffered polygonized diffs
        ## surface zero array and diffs...
        ## add surfaced diffs to self.dem
                        
        # # diff_cmd = 'gmt blockmedian {region} -I{xinc}/{yinc} | gmt surface {region} -I{xinc}/{yinc} -G_diff.tif=gd+n-9999:GTiff -T.1 -Z1.2 -V -rp -C.5 -Lud -Lld -M{radius}'.format(
        # #     region=dem_region, xinc=dem_infos['geoT'][1], yinc=-1*dem_infos['geoT'][5], radius=self.radius
        # # )
        # diff_cmd = 'gmt blockmedian {region} -I{xinc}/{yinc} | gmt surface {region} -I{xinc}/{yinc} -G_diff.tif=gd+n{ndv}:GTiff -T.1 -Z1.2 -V -Lud -Lld'.format(
        #     region=self.region, xinc=self.xinc, yinc=self.yinc, ndv=self.ndv, radius=self.radius
        # )

        # out, status = utils.run_cmd(
        #     diff_cmd,
        #     verbose=self.verbose,
        #     data_fun=lambda p: self.query_dump(
        #         dst_port=p, encode=True, max_diff=self.max_diff
        #     )
        # )

        #utils.echo_msg('smoothing diff grid...')
        #smooth_dem, status = demfun.blur('_diff.tif', '_tmp_smooth.tif', 5)

        #out, status = demfun.grdfilter('_diff.tif', '_tmp_smooth.tif', dist=self.radius, node='pixel', verbose=self.verbose)

        # if self.xinc != dem_infos['geoT']:
        #     utils.echo_msg('resampling diff grid...')
        #     if demfun.sample_warp('_diff.tif', '__tmp_sample.tif', dem_infos['geoT'][1], -1*dem_infos['geoT'][5],
        #                      src_region=dem_region)[1] == 0:
        #         os.rename('__tmp_sample.tif', '_tmp_smooth.tif')
        #     else:
        #         utils.echo_warning_msg('failed to resample diff grid')

        #if self.xinc != dem_infos['geoT']:
        utils.echo_msg('resampling diff grid...')
        diff_ds = demfun.sample_warp(
            '_diff.tif',
            None,
            dem_infos['geoT'][1],
            -1*dem_infos['geoT'][5],
            src_region=dem_region)[0]
        
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
            'description': '''SPLINE DEM via GMT surface\n
Generate a DEM using GMT's surface command
see gmt surface --help for more info.
        
< surface:tension=.35:relaxation=1.2:lower_limit=d:upper_limit=d >
 :tension=[0-1]\t\tspline tension.
 :relaxation=[val]\tspline relaxation factor.
 :lower_limit=[val]\tconstrain interpolation to lower limit.
 :upper_limit=[val]\tconstrain interpolation to upper limit.''',
        },
        'triangulate': {
            'name': 'triangulate',
            'datalist-p': True,
            'class': GMTTriangulate,
            'description': '''TRIANGULATION DEM via GMT triangulate\n
Generate a DEM using GMT's triangulate command.
see gmt triangulate --help for more info.        

< triangulate >''',
        },
        'nearneighbor': {
            'name': 'nearneihbor',
            'datalist-p': True,
            'class': GMTNearNeighbor,
            'description': """NEARNEIGHBOR DEM via GMT nearneighbor\n
Generate a DEM using GMT's nearneighbor command.
see gmt nearneighbor --help for more info.

< nearneighbor:radius=None:sector=None >
 :radius=[val]\t\tsearch radius
 :sectors=[val]\t\tsector information""",
        },
        'num': {
            'name': 'num',
            'datalist-p': True,
            'class': WafflesNum,
            'description': """Uninterpolated DEM populated by <mode>.\n
Generate an uninterpolated DEM using <mode> option.
Using mode of 'A<mode>' uses GMT's xyz2grd command, 
see gmt xyz2grd --help for more info.

mode keys: k (mask), m (mean), n (num), w (wet), A<mode> (gmt xyz2grd)

< num:mode=n >
 :mode=[key]\tspecify mode of grid population""",
        },
        'linear': {
            'name': 'linear',
            'datalist-p': True,
            'class': WafflesLinear,
            'description': """LINEAR DEM via gdal_grid\n
Generate a DEM using GDAL's gdal_grid command.
see gdal_grid --help for more info

< linear:radius=-1 >
 :radius=[val]\t\tsearch radius""",
        },
        'nearest': {
            'name': 'nearest',
            'datalist-p': True,
            'class': WafflesNearest,
            'description': """NEAREST DEM via gdal_grid\n
Generate a DEM using GDAL's gdal_grid command.
see gdal_grid --help for more info

< nearest:radius1=0:radius2=0:angle=0:nodata=0 >
 :radius1=[val]\t\tsearch radius 1
 :radius2=[val]\t\tsearch radius 2
 :angle=[val]\t\tangle
 :nodata=[val]\t\tnodata""",
        },
        'average': {
            'name': 'average',
            'datalist-p': True,
            'class': WafflesMovingAverage,
            'description': """Moving AVERAGE DEM via gdal_grid\n
Generate a DEM using GDAL's gdal_grid command.
see gdal_grid --help for more info

< average:radius1=0:radius2=0:angle=0:min_points=0:nodata=0 >
 :radius1=[val]\t\tsearch radius 1
 :radius2=[val]\t\tsearch radius 2
 :min_points=[val]\tminimum points per bucket (use to fill entire DEM)""",
        },
        'invdst': {
            'name': 'invdst',
            'datalist-p': True,
            'class': WafflesInvDst,
            'description': """INVERSE DISTANCE DEM via gdal_grid\n
Generate a DEM using GDAL's gdal_grid command.
see gdal_grid --help for more info

< invdst:power=2.0:smoothing=0.0:radius1=0:radius2=0:angle=0:max_points=0:min_points=0:nodata=0 >
 :radius1=[val]\t\tsearch radius 1
 :radius2=[val]\t\tsearch radius 2
 :power=[val]\t\tweight**power
 :min_points=[val]\tminimum points per IDW bucket (use to fill entire DEM)""",
        },
        'IDW': {
            'name': 'IDW',
            'datalist-p': True,
            'class': WafflesIDW,
            'description': """INVERSE DISTANCE WEIGHTED DEM\n
Generate a DEM using an Inverse Distance Weighted algorithm.
If weights are used, will generate a UIDW DEM, using weight values as inverse uncertainty,
as described here: https://ir.library.oregonstate.edu/concern/graduate_projects/79407x932
and here: https://stackoverflow.com/questions/3104781/inverse-distance-weighted-idw-interpolation-with-python

< IDW:min_points=8:radius=inf:power=1:block=False >
 :power=[val]\t\tweight**power
 :min_points=[val]\tminimum neighbor points for IDW
 :radius=[val]\t\tsearch radius, only fill data cells within radius from data
 :block=[True/False]\tblock the data before performing the IDW routine""",
        },
        'vdatum': {
            'name': 'vdatum',
            'datalist-p': False,
            'class': WafflesVDatum,
            'description': """VDATUM transformation grid.\n
Generate a Vertical DATUM transformation grid.

< vdatum:vdatum_in=None:vdatum_out=None >
 :vdatum_in=[vdatum]\tinput vertical datum
 :vdatum_out=[vdatum]\toutput vertical datum""",
        },
        'mbgrid': {
            'name': 'mbgrid',
            'datalist-p': True,
            'class': WafflesMBGrid,
            'description': """SPLINE DEM via MB-System's mbgrid.\n
Generate a DEM using MB-System's mbgrid command.
By default, will use MB-Systems datalist processes.
set `use_datalist=True` to use CUDEM's dlim instead.
see mbgrid --help for more info

< mbgrid:dist='10/3':tension=35:use_datalists=False >
 :dist=[val]\t\t\tthe dist variable to use in mbgrid
 :tension=[val]\t\t\tthe spline tension value (0-inf)
 :use_datalist=[True/False]\tuse built-in datalists rather than mbdatalist""",
        },
        'coastline': {
            'name': 'coastline',
            'datalist-p': False,
            'class': WafflesCoastline,
            'description': """COASTLINE (land/etc-mask) generation.\n
Generate a coastline (land/etc-mask) using a variety of sources. 
User data can be provided to provide source for further land masking. 
Output raster will mask land-areas as 1 and oceans/(lakes/buildings) as 0.
Output vector will polygonize land-areas.

< coastline:want_gmrt=False:want_nhd=True:want_lakes=False:want_buildings=False:invert=False:polygonize=True >
 :want_gmrt=[True/False]\tuse GMRT to fill background (will use Copernicus otherwise)
 :want_nhd=[True/False]\t\tuse high-resolution NHD to fill US coastal zones
 :want_lakes=[True/False\tmask LAKES using HYDROLAKES
 :want_buildings=[True/False]\tmask BUILDINGS using OSM
 :invert=[True/False]\t\tinvert the output results
 :polygonize=[True/False]\tpolygonize the output""",
        },
        'cudem': {
            'name': 'cudem',
            'datalist-p': True,
            'class': WafflesCUDEM,
            'description': """CUDEM integrated DEM generation. <beta>\n
Generate an topo/bathy integrated DEM using a variety of data sources.
Will iterate <pre_count> pre-surfaces at lower-resolutions.
Each pre-surface will be clipped to <landmask> if it exists and smoothed with <smoothing> factor.
Each pre-surface is used in subsequent pre-surface(s)/final DEM at each iterative weight.

< cudem:landmask=None:min_weight=[75th percentile]:smoothing=None:pre_count=1 >
 :landmask=[path]\tpath to coastline vector mask or set as `coastline` to auto-generate
 :min_weight=[val]\tthe minumum weight to include in the final DEM
 :smoothing=[val]\tthe Gaussian smoothing value to apply to pre-surface(s)
 :pre_count=[val]\tnumber of pre-surface iterations to perform""",
        },
        'patch': {
            'name': 'patch',
            'datalist-p': True,
            'class': WafflesPatch,
            'description': """PATCH an existing DEM with new data. <beta>\n
Patch an existing DEM with data from the datalist.

< patch:min_weight=.5:dem=None >
 :dem=[path]\t\tthe path the the DEM to update
 :min_weight=[val]\tthe minumum data weight to include in the patched DEM""",
        },
        'lakes': {
            'name': 'lakes',
            'datalist-p': False,
            'class': WafflesLakes,
            'description': """estimate lake bathymetry.\n
By default, will return lake bathymetry as depth values (positive down), 
to get elevations (positive up), set apply_elevations=True.

< lakes:apply_elevations=False:min_area=None:max_area=None:min_id=None:max_id=None:depth=globathy >
 :apply_elevations=[True/False]\t\tuse COPERNICUS to apply lake level elevations to output
 :min_area=[val]\t\t\tminimum lake area to consider
 :max_area=[val]\t\t\tmaximum lake area to consider
 :min_id=[val]\t\t\t\tminimum lake ID to consider
 :max_id=[val]\t\t\t\tmaximum lake ID to consider
 :depth=[globathy/hydrolakes/val]\tobtain the depth value from GloBathy, HydroLakes or constant value""",
        },
        'stacks': {
            'name': 'stacks',
            'datalist-p': True,
            'class': WafflesStacks,
            'description': """STACK data into a DEM. <beta>\n
Generate a DEM using a raster STACKing method. 
By default, will calculate the [weighted]-mean where overlapping cells occur. 
Set supercede to True to overwrite overlapping cells with higher weighted data.

< stacks:supercede=False:keep_weights=False:keep_count=False:min_count=None >
 :supercede=[True/False]\tsuperced data cells with higher weighted data
 :keep_weights=[True/False]\tretain weight raster
 :keep_count=[True/False]\tretain count raster
 :min_count=[val]\t\tonly retain data cells if they contain `min_count` overlapping data""",
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
                weight=self.weights
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
_waffles_module_short_desc = lambda: ', '.join(
    ['{}'.format(key) for key in WaffleFactory()._modules])
_waffles_module_long_desc = lambda x: 'waffles modules:\n% waffles ... -M <mod>:key=val:key=val...\n\n  ' + '\n  '.join(
    ['\033[1m{:14}\033[0m{}\n'.format(key, x[key]['description']) for key in x]) + '\n'

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
  --modules\t\tDisply the module descriptions and usage
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
           dl_formats=dlim._datalist_fmts_short_desc(),
           modules=_waffles_module_short_desc(),
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
            try:
                if argv[i + 1] in WaffleFactory()._modules.keys():
                    sys.stderr.write(
                        _waffles_module_long_desc({k: WaffleFactory()._modules[k] for k in (argv[i + 1],)})
                    )
                else: sys.stderr.write(
                        _waffles_module_long_desc(WaffleFactory()._modules)
                )
            except: sys.stderr.write(
                    _waffles_module_long_desc(WaffleFactory()._modules)
            )
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

    for i_region in i_regions:
        tmp_region = regions.Region().from_string(i_region)
        if tmp_region.valid_p(check_xy=True):
            these_regions.append(tmp_region)
        else:
            i_region_s = i_region.split(':')
            tmp_region = regions.ogr_wkts(i_region_s[0])
            for i in tmp_region:
                if i.valid_p():
                    if len(i_region_s) > 1:
                        these_regions.append(
                            regions.Region().from_string(
                                '/'.join([i.format('str'), i_region_s[1]])
                            )
                        )
                    else:
                        these_regions.append(i)
                    
    if len(these_regions) == 0:
        these_regions = []
        utils.echo_error_msg('failed to parse region(s), {}'.format(i_regions))
    else:
        if wg['verbose']: utils.echo_msg('parsed {} region(s)'.format(len(these_regions)))

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
            #this_wg = this_wg._config
            utils.echo_msg(json.dumps(this_wg, indent=4, sort_keys=True))
            with open('{}.json'.format(this_wg['name']), 'w') as wg_json:
                utils.echo_msg('generating waffles config file: {}.json'.format(this_wg['name']))
                wg_json.write(json.dumps(this_wg, indent=4, sort_keys=True))
        else:
            #this_waffle = WaffleFactory(mod=module, **wg).acquire()
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
