### waffles.py
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
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
## Generate DEMs from a variety of data sources.
## Use CLI command 'waffles'
##
## Supported input datatypes include:
## datalist, las/laz, gdal, bag, ogr, xyz, mbs, fetches
## see cudem.dlim for more information on supported datasets
##
## Supported gridding modules include:
## gmt-surface (GMT), gmt-triangulate (GMT), gmt-nearneighbor (GMT),
## mbgrid (MB-System), IDW (CUDEM), num (CUDEM/GMT), coastline (CUDEM),
## cudem (CUDEM), stacks (CUDEM), gdal-inv-dst (GDAL), gdal-linear (GDAL),
## gdal-average (GDAL), gdal-nearest (GDAL), linear (SCIPY), cubic (SCIPY),
## nearest (SCIPY), scratch (CUDEM)
##
## GMT, GDAL and MB-System are required for full functionality.
##
## Process all input data through cudem.dlim first, minimally, using
## cudem.dlim.init_data(list-of-data).
##
## Data will be processed through cudem.dlim._stacks prior to any interpolation.
## This will produce a weighted-mean (or weight-superceded)
## grid of the data to use for interpolation. See cudem.dlim for more information.
##
## The interpolated DEM will be post-processed through Waffle._process() where
## it will be filtered, resampled, clipped, set limits, cut to final output size,
## filled with metadata. Optionally, if want_uncertainty is set, _process() will
## also calculate the interpolation uncertainty.
##
### Code:


import sys
import os
import math
import json
import time
import glob
import traceback
from tqdm import tqdm
from tqdm import trange

import numpy as np
import h5py as h5
import netCDF4 as nc
import scipy
from scipy import interpolate
from scipy import spatial
from scipy import ndimage
import threading
import multiprocessing as mp
try:
   import Queue as queue
except: import queue as queue

from osgeo import gdal
from osgeo import ogr
from osgeo import osr

import cudem
from cudem import dlim
from cudem import regions
from cudem import utils
from cudem import xyzfun
from cudem import gdalfun
from cudem import vdatums
from cudem import factory
from cudem import fetches
from cudem import grits
from cudem import perspecto
from cudem import srsfun

## Data cache directory, hold temp data, fetch data, etc here.
waffles_cache = utils.cudem_cache()
gc = utils.config_check()
gdal.DontUseExceptions()
ogr.DontUseExceptions()
osr.DontUseExceptions()
gdal.SetConfigOption(
    'CPL_LOG', 'NUL' if gc['platform'] == 'win32' else '/dev/null'
)


class Waffle:
    """Representing a WAFFLES DEM/MODULE.
    Specific Gridding modules are sub-classes of this class.
    See WaffleFactory for module specifications and generation.

    -----------
    Procedures:
      yield_xyz() - yield the xyz data from self.data
      dump_xyz() - dump the xyz data from self.data to port
      run() - run the WAFFLES module (function set via module sub-class)
      generate() - run and process the WAFFLES module
    """
    
    def __init__(
            self,
            data: list = [],
            src_region: regions.Region = None,
            inc: str = None,
            xinc: str = None,
            yinc: str = None,
            xsize: int = None,
            ysize: int = None,
            name: str = 'waffles_dem',
            node: str = 'pixel',
            fmt: str = 'GTiff',
            extend: int = 0,
            extend_proc: float = 0,
            want_weight: bool = False,
            want_uncertainty: bool = False,
            fltr: list = [],
            sample: str = 'bilinear',
            xsample: str = None,
            ysample: str = None,
            clip: str = None,
            chunk: int = None,
            dst_srs: str = None,
            srs_transform: bool = False,
            verbose: bool = False,
            archive: bool = False,
            want_mask: bool = False,
            keep_auxiliary: bool = False,
            want_sm: bool = False,
            clobber: bool = True,
            ndv: float = -9999,
            block: bool = False,
            cache_dir: str = waffles_cache,
            stack_mode: str = 'mean',
            upper_limit: float = None,
            lower_limit: float = None,
            proximity_limit: int = None,
            size_limit: int = None,
            percentile_limit: float = None,
            expand_limit: float = None,
            count_limit: int = None,
            uncertainty_limit: int = None,
            flatten_nodata_values: bool = False,
            want_stack: bool = True,
            co: list = [],
            params: dict = {}
    ):
        self.params = params # the factory parameters
        self.data = data # list of data paths/fetches modules to grid
        self.datalist = None # the datalist which holds the processed datasets
        self.region = src_region # the region to grid
        self.src_region = src_region # the region to grid
        self.inc = inc # the gridding increments [xinc, yinc]
        self.xinc = xinc # the x/lon gridding increment
        self.yinc = yinc # the y/lat gridding increment
        self.sample = sample # the gdal sample algorithm to use when needed
        self.xsample = xsample # the x/lon increment to sample the output dem
        self.ysample = ysample # the y/lat increment to sample the output dem
        self.name = name # the output dem basename
        self.node = node # the grid node method, either 'grid' or 'pixel'
        self.fmt = fmt # the gdal-compatible output dem file format
        self.extend = extend # extend the dem region by this many pixels
        self.extend_proc = extend_proc # extend the dem processing region by this percentage
        self.want_weight = want_weight # use weights, either None or 1
        self.want_uncertainty = want_uncertainty # apply/calculate uncertainty
        self.fltr = fltr # a list of filters (see cudem.grits for options)
        self.clip = clip # ogr compatible vector file or keyword module to clip output dem
        self.chunk = chunk # process the dem in this many chunks
        self.dst_srs = dst_srs # the output dem projection
        self.srs_transform = srs_transform # transform data to the dst_srs
        self.archive = archive # archive the data used in this dem
        self.want_mask = want_mask # mask the incoming datalist
        #self.supercede = supercede # higher weighted data supercedes lower weighted data
        self.stack_mode = stack_mode
        self.upper_limit = utils.float_or(upper_limit) # upper limit of z values in output
        self.lower_limit = utils.float_or(lower_limit) # lower limit of z values in output
        self.proximity_limit = utils.int_or(proximity_limit) # proximity limit of interpolation
        self.size_limit = utils.int_or(size_limit) # size limit of interpolation
        self.percentile_limit = utils.float_or(percentile_limit) # percentile limit of interpolation
        self.expand_limit = utils.float_or(expand_limit) # percentile limit of interpolation
        self.count_limit = utils.int_or(count_limit) # limit by stack count per cell
        self.uncertainty_limit = utils.int_or(uncertainty_limit) # limit by stack upper uncertainty
        self.keep_auxiliary = keep_auxiliary # keep auxiliary outputs
        self.clobber = clobber # clobber the output dem file
        self.verbose = verbose # increase verbosity
        self.cache_dir = cache_dir # directory path to store cahced data
        self.ndv = ndv # no data value for the dem
        self.block = block # block the data (defunct)
        self.block_t = None # block the data (defunct)
        self.ogr_ds = None # datasets as an ogr object
        self.data_ = data # store data paths here, self.data gets processed to dlim datasets
        self.fn = '{}.tif'.format(self.name) # output dem filename
        self.want_sm = want_sm # generate spatial metadata
        self.aux_dems = [] # list of auxiliary dems fns
        self.want_stack = want_stack # generate the stacked rasters
        self.stack = None # multi-banded stacked raster from cudem.dlim
        self.stack_ds = None # the stacked raster as a dlim dataset object
        self.co = co # the gdal creation options
        self.flatten_nodata_values = flatten_nodata_values # flatten any remaining nodata values
        self.output_files = {}
        self.status = 0

        
    def __str__(self):
        return('<Waffles: {}>'.format(self.name))

    
    def __repr__(self):
        return('<Waffles: {}>'.format(self.name))

    
    def initialize(self):
        if self.verbose:
            utils.echo_msg(
                f'initializing waffles module < \033[1m{self.params["mod"]}\033[m >'
            )

        # output dem filename
        self.fn = '{}.{}'.format(self.name, gdalfun.gdal_fext(self.fmt))
        # cudem config file holding foriegn programs and versions
        self.gc = utils.config_check()
        # initialize regions
        self._init_regions()
        # initialize increments
        self._init_incs() 

        if isinstance(self.co, list):
            if len(self.co) == 0:
                self.co = ["COMPRESS=DEFLATE", "TILED=YES"]
                
        else:
            self.co = ["COMPRESS=DEFLATE", "TILED=YES"]

        ## initialize data, setting set_incs to True will force dlim to process the
        ## data to the set increments
        if self.want_stack:
            self._init_data(set_incs=True) 

        self.xcount, self.ycount, self.dst_gt \
            = self.p_region.geo_transform(
                x_inc=self.xinc, y_inc=self.yinc, node='grid'
            )
        self.ds_config = gdalfun.gdal_set_infos(
            self.xcount, self.ycount, (self.xcount*self.ycount),
            self.dst_gt, gdalfun.osr_wkt(self.dst_srs),
            gdal.GDT_Float32, self.ndv, self.fmt, None, None
        )
        if self.verbose:
            utils.echo_msg(
                f'output size: {self.ds_config["nx"]}/{self.ds_config["ny"]}'
            )
        
        self.status = self._init()
        return(self)            

    
    def _init(self):
        return(0)

    
    def __call__(self):
        self.initialize()
        if self.status == 0:
            return(self.generate())
        else:
            utils.echo_warning_msg(
                'failed to initialize from sub-module'
            )

            
    def _init_regions(self):
        """Initialize and set regions.
        
        regions set here include:
        d_region: Distribution Region
        p_region: Processing Region
        c_region: Coastline Region
        ps_region: Processing region for GMT/MB-System
        """
        
        if isinstance(self.region, list):
            self.region = regions.Region().from_list(self.region)
        elif not isinstance(self.region, regions.Region):
            raise ValueError(
                f'could not parse region: {self.region}'
            )
        
        if self.node == 'grid':
            self.region = self.region.buffer(
                x_bv=self.xinc*.5, y_bv=self.yinc*.5
            )
            
        self.d_region = self._dist_region()
        self.p_region = self._proc_region()
        self.c_region = self._coast_region()
        self.ps_region = self.p_region.copy()
        self.ps_region = self.ps_region.buffer(
            x_bv=self.xinc*.5, y_bv=self.yinc*.5,
            x_inc=self.xinc, y_inc=self.yinc
        )

        if self.verbose:
            utils.echo_msg(f'input region: {self.region}')
            utils.echo_msg(f'distribution region: {self.d_region}')
            utils.echo_msg(f'processing region: {self.p_region}')

            
    def _init_data(self, set_incs=False):
        """Initialize the data for processing
        parses data paths to dlim dataset objects.

        set `set_incs` to True to block/sample datasets to given 
        increment

        this function sets `self.data` to a list of dataset objects.
        """

        stack_fltr=[]
        point_fltr=[]
        if isinstance(self.fltr, list):
            for f in self.fltr:
                if f.split(':')[0] in grits.GritsFactory._modules.keys():
                    grits_filter = grits.GritsFactory(mod=f)._acquire_module()
                    if grits_filter is not None:
                        if 'stacks' in grits_filter.kwargs.keys():
                            if grits_filter.kwargs['stacks']:
                                stack_fltr.append(f)
                elif f.split(':')[0] in dlim.PointFilterFactory._modules.keys():
                    point_filter = dlim.PointFilterFactory(mod=f)._acquire_module()
                    if point_filter is not None:
                        point_fltr.append(f)

        self.data = dlim.init_data(
            self.data,
            region=self.p_region,
            src_srs=None,
            dst_srs=self.dst_srs,
            xy_inc=(self.xinc, self.yinc),
            sample_alg=self.sample,
            want_weight=self.want_weight,
            want_uncertainty=self.want_uncertainty,
            want_verbose=self.verbose,
            want_mask=self.want_mask,
            pnt_fltrs=point_fltr,
            stack_fltrs=stack_fltr,
            invert_region=False,
            cache_dir=self.cache_dir,
            stack_mode=self.stack_mode,
        )

        if self.data is not None:
            self.data.initialize()
            if not self.want_weight:
                self.data.weight = None
                
            if not self.want_uncertainty:
                self.data.uncertainty = None
        else:
            return(None)

        
    def _init_incs(self):
        """Initialize increments
        
        xinc/yinc are the DEM increments, in native units.
        xsample/ysample set the output DEM increments, 
        in native units.
        """
        
        self.xinc = utils.str2inc(self.xinc)
        self.yinc = utils.str2inc(self.yinc)
        self.xsample = utils.str2inc(self.xsample)
        self.ysample = utils.str2inc(self.ysample)

        if self.verbose:
            utils.echo_msg(
                f'gridding increments: {self.xinc}/{self.yinc}'
            )
            utils.echo_msg(
                'output increments: {}/{}'.format(
                    self.xsample \
                    if self.xsample is not None \
                    else self.xinc,
                    self.ysample \
                    if self.ysample is not None \
                    else self.yinc
                )
            )

            
    def _coast_region(self):
        """coastline region 
        (extended by percentage self.extend_proc)
        """

        cr = self.d_region.copy()
        return(
            cr.buffer(
                pct=self.extend_proc, x_inc=self.xinc, y_inc=self.yinc
            )
        )

    
    def _proc_region(self):
        """processing region 
        (extended by percentage self.extend_proc)
        """

        pr = self.d_region.copy()
        return(
            pr.buffer(
                pct=self.extend_proc, x_inc=self.xinc, y_inc=self.yinc
            )
        )

    
    def _dist_region(self):
        """distribution region 
        (extended by self.extend).
        """
        
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

        
    def dump_xyz(self, dst_port=sys.stdout, encode=False):
        """dump the stacked xyz data to dst_port

        use this to dump data into a foreign cli program, 
        such as GMT.
        """

        for xyz in self.stack_ds.yield_xyz():
            xyz.dump(
                include_w = self.want_weight,
                include_u = self.want_uncertainty,
                dst_port=dst_port,
                encode=encode,
            )

            
    def generate(self):
        """run and process the WAFFLES module.

        Generate the data 'stack' and use that to run the 
        waffles module and generate the DEM.
        """

        if self.data is None:
            return(self)

        #self.output_files = {}
        ## todo: move to init
        if os.path.exists(self.fn):
            if not self.clobber:
                utils.echo_warning_msg(
                    f'DEM {self.fn} already exists, skipping...'
                )
                return(self)
            else:
                utils.echo_warning_msg(
                    f'DEM {self.fn} exists and will be clobbered.'
                )
                status = gdal.GetDriverByName(self.fmt).Delete(self.fn)
                if status != 0:
                    utils.remove_glob(f'{self.fn}*')

        else:
            if not os.path.exists(os.path.dirname(self.fn)):
                try:
                    os.makedirs(os.path.dirname(self.fn))
                except: pass

        if self.chunk is not None:
            ## Generate in Chunks of self.chunk by self.chunk and
            ## merge chunks into final DEM
            chunks = []
            stack_chunks = []
            mask_chunks = []
            aux_chunks = []
            for srcwin in utils.yield_srcwin(
                    (self.ycount, self.xcount),
                    self.chunk,
                    verbose=self.verbose
            ):
                this_geo_x_origin, this_geo_y_origin \
                    = utils._pixel2geo(
                        srcwin[0], srcwin[1], self.dst_gt
                    )
                this_geo_x_end, this_geo_y_end \
                    = utils._pixel2geo(
                        srcwin[0]+srcwin[2], srcwin[1]+srcwin[3], self.dst_gt
                    )
                this_gt = [
                    this_geo_x_origin,
                    float(self.dst_gt[1]),
                    0.0,
                    this_geo_y_origin,
                    0.0,
                    float(self.dst_gt[5])
                ]
                this_region = regions.Region()
                this_region.from_geo_transform(
                    geo_transform=this_gt, x_count=srcwin[2], y_count=srcwin[3]
                )
                this_region.buffer(pct=10, x_inc = self.xinc, y_inc = self.yinc)
                this_params = self.params.copy()
                this_params['kwargs']['src_region'] = this_region
                this_params['kwargs']['chunk'] = None
                this_params['kwargs']['name'] \
                    = utils.append_fn(
                        '_chunk', this_region, self.xinc, high_res=True
                    )
                this_waffle = WaffleFactory().load_parameter_dict(this_params)
                this_waffle_module = this_waffle._acquire_module()
                this_waffle_module.initialize()
                this_waffle_module.generate()

                ## append the chunk to the chunk list if the DEM generated ok
                if WaffleDEM(
                        this_waffle_module.fn,
                        cache_dir=self.cache_dir,
                        verbose=self.verbose
                ).initialize().valid_p():
                    chunks.append(this_waffle_module.fn)

                if WaffleDEM(
                        this_waffle_module.stack,
                        cache_dir=self.cache_dir,
                        verbose=self.verbose
                ).initialize().valid_p():
                    stack_chunks.append(this_waffle_module.stack)
                    
                mask_name = f'{this_waffle_module.name}_msk'
                mask_fn = '{}.{}'.format(
                    os.path.join(this_waffle_module.cache_dir, mask_name),
                    gdalfun.gdal_fext(this_waffle_module.fmt)
                )

                mask_chunks.append(this_waffle_module.msk_fn)
                aux_chunks.append(this_waffle_module.aux_dems)

            ## combine DEM chunks
            if len(chunks) > 0:
                g = gdal.Warp(
                    self.fn, chunks,
                    format=self.fmt,
                    resampleAlg='cubicspline',
                    options=["COMPRESS=LZW", "TILED=YES"]
                )
                g = None

            ## combine STACK chunks
            ## multi-band
            if len(stack_chunks) > 0:
                g = gdal.Warp(
                    self.fn,
                    stack_chunks,
                    format=self.fmt,
                    resampleAlg='cubicspline',
                    options=["COMPRESS=LZW", "TILED=YES"]
                )
                g = None

            ## combine MASK chunks
            ## multi-band
            if len(mask_chunks) > 0:
                g = gdal.Warp(
                    mask_fn,
                    mask_chunks,
                    format=self.fmt,
                    resampleAlg='cubicspline',
                    options=["COMPRESS=LZW", "TILED=YES"]
                )
                g = None
                
            ## combine AUXILIARY chunks
            if len(aux_chunks) > 0:
                for ac, aux_dem in enumerate(aux_chunks):
                    g = gdal.Warp(
                        aux_dem,
                        aux_chunks[ac],
                        format=self.fmt,
                        resampleAlg='cubicspline',
                        options=["COMPRESS=LZW", "TILED=YES"]
                    )
                    g = None
                
            utils.remove_glob(*chunks)
            for aux_chunk in aux_chunks:
                utils.remove_glob(aux_chunk)
        else:
            ## stack the data and run the waffles module
            mask_fn = None
            if self.want_stack:
                stack_name = f'{os.path.basename(self.name)}_stack'
                mask_name = f'{stack_name}_msk'
                mask_fn = (f'{os.path.join(self.cache_dir, mask_name)}'
                           f'.{gdalfun.gdal_fext(self.fmt)}')

                ## Threaded...no worky
                # num_threads = 8
                # try:
                #     sk = dlim.stacks_ds(
                #         self.data,
                #         n_threads=num_threads,
                #         out_name=os.path.join(self.cache_dir, stack_name),
                #         supercede=self.supercede,
                #         want_mask=self.want_mask
                #     )
                #     sk.daemon = True
                
                #     sk.start()
                #     sk.join()                
                # except (KeyboardInterrupt, SystemExit):
                #     utils.echo_error_msg(
                #         'user breakage...please wait while fetches exits.'
                #     )
                #     stop_threads = True
                #     while not sk.arr_q.empty():
                #         try:
                #             sk.arr_q.get(False)
                #         except Empty:
                #             continue
                
                #         sk.arr_q.task_done()                        
                
                # self.stack = sk.out_file
                
                ## generate the stack
                stack_fn = os.path.join(
                    self.cache_dir,
                    f'{stack_name}.{gdalfun.gdal_fext("GTiff")}'
                )
                ## for stacks_h5
                # stack_fn = os.path.join(
                #     self.cache_dir, f'{stack_name}.csg'
                # )
                stack_bn = utils.fn_basename2(stack_fn)
                if not self.clobber and os.path.exists(stack_fn):
                    self.stack = stack_fn
                    if not WaffleDEM(
                            self.stack,
                            cache_dir=self.cache_dir,
                            verbose=self.verbose
                    ).initialize().valid_p():
                        self.stack = self.data._stacks(out_name=stack_bn)
                        #, mode=self.stack_mode)#supercede=self.supercede)
                else:
                    self.stack = self.data._stacks(out_name=stack_bn)
                    #, mode=self.stack_mode)#supercede=self.supercede)
                    
                self.stack_ds = dlim.GDALFile(
                    fn=self.stack,
                    band_no=1,
                    weight_mask=3,
                    uncertainty_mask=4,
                    data_format=200,
                    src_srs=self.dst_srs,
                    dst_srs=self.dst_srs,
                    x_inc=self.xinc,
                    y_inc=self.yinc,
                    src_region=self.p_region,
                    weight=1,
                    verbose=self.verbose
                ).initialize()

                ## rename the mask_fn outside of cachedir
                if self.want_mask:
                    mask_fn = f'{mask_name}.{gdalfun.gdal_fext(self.fmt)}'

                ## apply the count or uncertainty limit to the stack grid
                if self.count_limit is not None \
                   or self.uncertainty_limit is not None:
                    with gdalfun.gdal_datasource(self.stack, update=True) as stack_ds:
                        stack_infos = gdalfun.gdal_infos(stack_ds)

                        if self.count_limit:
                            count_band = stack_ds.GetRasterBand(2)
                            count_arr = count_band.ReadAsArray()
                            count_mask = count_arr <= self.count_limit
                            
                        if self.uncertainty_limit:
                            uncertainty_band = stack_ds.GetRasterBand(4)
                            uncertainty_arr = uncertainty_band.ReadAsArray()
                            uncertainty_mask = uncertainty_arr <= self.uncertainty_limit

                        for band in range(1, stack_ds.RasterCount+1):
                            this_band = stack_ds.GetRasterBand(band)
                            this_arr = this_band.ReadAsArray()
                            if self.count_limit:
                                this_arr[count_mask] = stack_infos['ndv']
                                
                            if self.uncertainty_limit:
                                this_arr[uncertainty_mask] = stack_infos['ndv']
                                
                            this_band.WriteArray(this_arr)
                            
                ## run the waffles module
                if WaffleDEM(
                        self.stack,
                        cache_dir=self.cache_dir,
                        verbose=self.verbose
                ).initialize().valid_p():
                    self.run()
                    
            else:
                self.run()    
                
            # if self.node == 'grid':
            #     self.region = self.region.buffer(
            #         x_bv=-self.xinc*.5, y_bv=-self.yinc*.5
            #     )
            ## calculate estimated uncertainty of the interpolation
            unc_fn = None
            if self.want_uncertainty:
                iu = WaffleFactory(
                    mod=('uncertainty:percentile=95:accumulate=False'
                         f':waffles_module={self.params["mod"]}'),
                    **self.params['kwargs']
                )._acquire_module()
                iu.name = f'{self.params["kwargs"]["name"]}_u'
                iu.want_uncertainty = False
                iu.want_mask = False
                iu.stack = self.stack
                iu.initialize()
                iu.run()

                unc_dem = WaffleDEM(
                    iu.fn,
                    cache_dir=self.cache_dir,
                    verbose=self.verbose,
                    want_scan=True,
                    co=self.co
                ).initialize()
                if unc_dem.valid_p():
                    unc_dem.process(
                        ndv=self.ndv,
                        xsample=self.xsample,
                        ysample=self.ysample,
                        region=self.d_region,
                        clip_str=self.clip,
                        node=self.node,
                        dst_srs=self.dst_srs,
                        dst_fmt=self.fmt,
                        set_metadata=False,
                        dst_dir=os.path.dirname(self.fn)
                    )
                    self.output_files['uncertainty'] = iu.fn
                unc_fn = iu.fn
            
            ## post-process the DEM(s)
            waffle_dem = WaffleDEM(
                self.fn,
                cache_dir=self.cache_dir,
                verbose=self.verbose,
                co=self.co
            ).initialize()
            if waffle_dem.valid_p():
                if mask_fn is not None and os.path.exists(mask_fn):
                    ## set the projection to the mask,
                    ## it gets lost in translation from mem to tif
                    if self.dst_srs is not None:
                        gdalfun.gdal_set_srs(mask_fn, src_srs=self.dst_srs)
                    
                waffle_dem.process(
                    ndv=self.ndv,
                    xsample=self.xsample,
                    ysample=self.ysample,
                    region=self.d_region,
                    clip_str=self.clip,
                    node=self.node,
                    upper_limit=self.upper_limit,
                    lower_limit=self.lower_limit,
                    size_limit=self.size_limit,
                    proximity_limit=self.proximity_limit,
                    percentile_limit=self.percentile_limit,
                    expand_limit=self.expand_limit,
                    dst_srs=self.dst_srs,
                    dst_fmt=self.fmt,
                    stack_fn=self.stack,
                    mask_fn=mask_fn,
                    unc_fn=unc_fn,
                    filter_=self.fltr,
                    flatten_nodata_values=self.flatten_nodata_values,
                    want_nc=self.keep_auxiliary,
                    want_h5=self.keep_auxiliary
                )
                self.output_files['DEM'] = self.fn
                
            ## post-process the mask, etc.
            if self.want_sm:
                if mask_fn is not None:
                    mask_dem = WaffleDEM(
                        mask_fn,
                        cache_dir=self.cache_dir,
                        verbose=self.verbose,
                        want_scan=True,
                        co=self.co
                    ).initialize()
                    if mask_dem.valid_p():
                        mask_dem.process(
                            ndv=0,
                            xsample=self.xsample,
                            ysample=self.ysample,
                            region=self.d_region,
                            clip_str=self.clip,
                            node=self.node,
                            dst_srs=self.dst_srs,
                            dst_fmt=self.fmt,
                            set_metadata=False,
                            dst_fn='{}_msk.{}'.format(
                                os.path.basename(self.name),
                                gdalfun.gdal_fext(self.fmt)
                            ),
                            dst_dir=os.path.dirname(self.fn)
                        )
                        self.output_files['mask'] = mask_dem.fn
                        if self.want_sm:
                            self.output_files['spatial-metadata'] = []
                            ## with gdal_footprint
                            ## gdal_footprint can stall on a large number of mask bands
                            ## has_gdal_footprint = utils.cmd_exists('gdal_footprint')
                            with gdalfun.gdal_datasource(mask_dem.fn) as msk_ds:
                                # if has_gdal_footprint:
                                #     sm_files, sm_fmt = dlim.ogr_mask_footprints(
                                #         msk_ds, verbose=True, mask_level=0
                                #     )
                                # else:
                                sm_layer, sm_fmt = dlim.polygonize_mask_multibands(
                                    msk_ds, verbose=True
                                )
                                sm_files = glob.glob(f'{sm_layer}.*')
                            
                            for f in sm_files:
                                out_sm = '{}_sm.{}'.format(
                                    self.name, f.split('.')[-1]
                                )
                                if os.path.exists(out_sm):
                                    utils.remove_glob(out_sm)

                                os.rename(f, out_sm)
                                self.output_files['spatial-metadata'].append(out_sm)
                else:
                    utils.echo_warning_msg(
                        'mask DEM is invalid...'
                    )
                    
            ## post-process any auxiliary rasters
            for aux_dem in self.aux_dems:
                aux_dem = WaffleDEM(
                    aux_dem,
                    cache_dir=self.cache_dir,
                    verbose=False,
                    co=self.co
                ).initialize()
                if aux_dem.valid_p():
                    aux_dem.process(
                        ndv=None,
                        xsample=self.xsample,
                        ysample=self.ysample,
                        region=self.d_region,
                        clip_str=self.clip,
                        node=self.node,
                        dst_srs=self.dst_srs,
                        dst_fmt=self.fmt,
                        dst_dir=os.path.dirname(self.fn),
                        set_metadata=False
                    )
                    self.output_files['stack'] = aux_dem.fn

            if self.keep_auxiliary:
                self.output_files['aux netcdf'] = f'{self.name}.nc'
            ## reset the self.stack to new post-processed fn and ds
            # if self.want_stack and self.keep_auxiliary:
            #     self.stack = os.path.join(
            #         os.path.dirname(self.fn), os.path.basename(self.stack)
            #     )
            #     self.stack_ds = dlim.GDALFile(
            #         fn=self.stack,
            #         band_no=1,
            #         weight_mask=3,
            #         uncertainty_mask=4,
            #         data_format=200,
            #         src_srs=self.dst_srs,
            #         dst_srs=self.dst_srs,
            #         x_inc=self.xinc,
            #         y_inc=self.yinc,
            #         src_region=self.p_region,
            #         weight=1,
            #         verbose=self.verbose
            #     ).initialize()

        if self.verbose:
            utils.echo_msg(
                f'output files: {self.output_files}'
            )
            
        return(self)                
        
    def run(self):
        """run the WAFFLES module (set via sub-module class)."""
        
        raise(NotImplementedError)

    
class WafflesScratch(Waffle):
    """SCRATCH Module. 
    
    Don't generate any DEMs, only auxiliary data, 
    including the raster stack.
    
    -----------
    Parameters:
    
    min_count=[val] - only retain data cells if they contain 
                      `min_count` overlapping data

    < scratch:min_count=None >
    """

    ## todo: add parameters for specifying outputs...
    def __init__(self, min_count=None, **kwargs):
        super().__init__(**kwargs)
        self.min_count = min_count

        
    def run(self):
        return(self)

    
class WafflesStacks(Waffle):
    """STACK data into a DEM. 
    
    Generate a DEM using a raster STACKing method. 
    By default, will calculate the [weighted]-mean where 
    overlapping cells occur. 

    Set `stack_mode` to 'supercede' to True to overwrite overlapping 
    cells with higher weighted data.

    stack data to generate DEM. No interpolation
    occurs with this module. To guarantee a full DEM,
    use a background DEM with a low weight, such as GMRT or GEBCO,
    which will be stacked upon to create the final DEM.
    
    -----------
    Parameters:
    
    min_count=[val] - only retain data cells if they contain `min_count` 
                      overlapping data
    
    < stacks:min_count=None >
    """

    ## todo: add parameters for specifying outputs...
    def __init__(self, min_count = None, **kwargs):
        super().__init__(**kwargs)
        self.min_count = min_count
        
    def run(self):
        z_ds = gdal.GetDriverByName(self.fmt).Create(
            '{}.{}'.format(self.name, gdalfun.gdal_fext(self.fmt)),
            self.xcount,
            self.ycount,
            1,
            self.ds_config['dt'],
            options=self.co
        )
        z_ds.SetGeoTransform(self.dst_gt)
        z_band = z_ds.GetRasterBand(1)
        z_band.SetNoDataValue(self.ndv)
        for arrs, srcwin, gt in self.stack_ds.yield_array():
            out_z = (arrs['z'] / arrs['weight']) / arrs['count']
            out_z[np.isnan(out_z)] = self.ndv
            z_band.WriteArray(out_z, srcwin[0], srcwin[1])
            
        z_ds = None
        if self.verbose:
            utils.echo_msg(
                f'stacked data to {self.fn}'
            )
            
        return(self)

    
class WafflesFlatten(Waffle):
    """Stack the data into a DEM and then hydro-flatten all the 
    void areas.

    specify 'size_threshold' to only flatten voids above threshold.

    -----------
    Parameters:
    
    min_count=[val] - only retain data cells if they contain `min_count` 
                      overlapping data
    size_threshold=[val] - the minimum size void to flatten (in cells)
    """
    
    ## todo: add parameters for specifying outputs...
    def __init__(self, min_count=None, size_threshold=1, **kwargs):
        super().__init__(**kwargs)
        self.min_count = min_count
        self.size_threshold = size_threshold

        
    def run(self):
        gdalfun.cudem_flatten_no_data_zones(
            self.stack,
            dst_dem=self.fn,
            band=1,
            size_threshold=self.size_threshold
        )
        
        return(self)

    
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

    def __call__(
            self, q,
            nnear=6,
            eps=0,
            p=1,
            dub=np.inf,
            weights=None,
            uncertainties=None
    ):
        # nnear nearest neighbours of each query point --
        q = np.asarray(q)
        qdim = q.ndim
        if qdim == 1:
            q = np.array([q])
            
        if self.wsum is None:
            self.wsum = np.zeros(nnear)

        eps = utils.float_or(eps, .1)
        self.distances, self.ix \
            = self.tree.query(q, k=nnear, eps=eps, distance_upper_bound=dub)
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

                # if uncertainties is not None:
                #     w *= (1/uncertainties[ix])
                    
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
    If weights are used, will generate a UIDW DEM, using weight 
    values as inverse uncertainty, as described here: 
    https://ir.library.oregonstate.edu/concern/graduate_projects/79407x932
    and here: 
    https://stackoverflow.com/questions/3104781/inverse-distance-weighted-idw-interpolation-with-python

    -----------
    Parameters:
    
    power=[val] - weight**power
    min_points=[val] - minimum neighbor points for IDW
    radius=[val] - search radius (in cells), only fill data cells 
                   within radius from data
    chunk_size=[val] - size of chunks in pixels

    < IDW:min_points=8:radius=inf:power=1:chunk_size=None >
    """
    
    def __init__(
            self,
            power=1,
            min_points=8,
            radius=None,
            chunk_size=None,
            **kwargs
    ):
        super().__init__(**kwargs)
        self.power = utils.float_or(power)
        self.min_points = utils.int_or(min_points)
        self.radius = np.inf if radius is None else utils.str2inc(radius) 
        self.chunk_size = chunk_size
        self.chunk_step = None

        
    def run(self):
        if self.verbose:
            if self.min_points:
                utils.echo_msg(
                    (f'generating IDW grid @ {self.ycount}/{self.xcount} '
                     f'looking for at least {self.min_points} neighbors '
                     f'within {self.radius} pixels')
                )
            else:
                utils.echo_msg(
                    f'generating IDW grid @ {self.ycount}/{self.xcount}'
                )
            i=0

        if self.chunk_size is None:
            n_chunk = int(self.ds_config['nx'] * .1)
            n_chunk = 10 if n_chunk < 10 else n_chunk
        else:
            n_chunk = self.chunk_size
            
        if self.chunk_step is None:
            n_step = int(n_chunk/2)
        else:
            n_step = self.chunk_step

        stack_ds = gdal.Open(self.stack)
        points_band = stack_ds.GetRasterBand(1)
        points_no_data = points_band.GetNoDataValue()
        weights_band = stack_ds.GetRasterBand(3)
        weights_no_data = weights_band.GetNoDataValue()
        uncertainty_band = stack_ds.GetRasterBand(4)
        uncertainty_no_data = uncertainty_band.GetNoDataValue()

        #try:
        interp_ds = stack_ds.GetDriver().Create(
            self.fn,
            stack_ds.RasterXSize,
            stack_ds.RasterYSize,
            bands=1,
            eType=points_band.DataType,
            options=["BLOCKXSIZE=256",
                     "BLOCKYSIZE=256",
                     "TILED=YES",
                     "COMPRESS=LZW",
                     "BIGTIFF=YES"]
        )
        interp_ds.SetProjection(stack_ds.GetProjection())
        interp_ds.SetGeoTransform(stack_ds.GetGeoTransform())
        interp_band = interp_ds.GetRasterBand(1)
        interp_band.SetNoDataValue(np.nan)
        #except:
        #    return(self)

        points_array = points_band.ReadAsArray()
        #points_array[points_array == points_no_data] = np.nan
        #point_indices = np.nonzero(~np.isnan(points_array))
        point_indices = np.nonzero(points_array != points_no_data)
        point_values = points_array[point_indices]
        points_array = None

        if self.want_weight:
            weights_array = weights_band.ReadAsArray()
            weight_values = weights_array[point_indices]
            weights_array = None
        else:
            weight_values = None

        if self.want_uncertainty:
            uncertainty_array = uncertainty_band.ReadAsArray()
            uncertainty_values = uncertainty_array[point_indices]
            uncertainty_array = None
        else:
            uncertainty_values = None

        stack_ds = None
        invdisttree = Invdisttree(
            np.transpose(point_indices),
            point_values,
            leafsize=10,
            stat=1
        )
        for srcwin in utils.yield_srcwin(
                (self.ycount, self.xcount),
                n_chunk=n_chunk,
                msg='Generating IDW DEM',
                end_msg='Generated IDW DEM',
                verbose=self.verbose
        ):
            # if np.count_nonzero(np.isnan(points_array)) == 0:
            #     #if not np.all(~np.nan(points_array)):
            #     # y_origin = srcwin[1]-srcwin_buff[1]
            #     # x_origin = srcwin[0]-srcwin_buff[0]
            #     # y_size = y_origin + srcwin[3]
            #     # x_size = x_origin + srcwin[2]
            #     # points_array = points_array[y_origin:y_size,x_origin:x_size]
            #     interp_band.WriteArray(points_array, srcwin[0], srcwin[1])

            if len(point_indices[0]):
                xi, yi = np.mgrid[srcwin[0]:srcwin[0]+srcwin[2],
                                  srcwin[1]:srcwin[1]+srcwin[3]]
                interp_data = invdisttree(
                    np.vstack((yi.flatten(), xi.flatten())).T,
                    nnear=self.min_points,
                    eps=.1,
                    p=self.power,
                    dub=self.radius,
                    weights=weight_values,
                    uncertainties=uncertainty_values
                )
                interp_data = np.reshape(interp_data, (srcwin[2], srcwin[3]))
                interp_band.WriteArray(interp_data.T, srcwin[0], srcwin[1])
                
        interp_ds = point_values = weight_values = None
        
        return(self)    

    
def write_array_queue(wq, q, m):
    while True:
        wq_args = wq.get()
        m.interp_band.WriteArray(wq_args[0], wq_args[1][0], wq_args[1][1])
        wq.task_done()

        
def scipy_queue(q, wq, m, p):
    while True:
        this_srcwin = q.get()
        p.update()
        try:
            interp_array = m.grid_srcwin(this_srcwin)
        except Exception as e:
            utils.echo_msg(e)
            utils.echo_warning_msg(
                f'failed to grid srcwin {this_srcwin}, placing back into queue'
            )
            q.put(this_srcwin)
            q.task_done()
            continue

        wq.put([interp_array, this_srcwin])
        q.task_done()

        
class grid_scipy(threading.Thread):
    """Scipy interpolate.griddata gridding 
    (linear, cubic, nearest)

    https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.griddata.html
    
    set 'num_threads' to int over 1 to generate in multiple threads...
    """
    
    def __init__(self, mod, n_threads=3):
        threading.Thread.__init__(self)
        self.mod = mod
        self.scipy_q = queue.Queue()
        self.grid_q = queue.Queue()
        self.n_threads = n_threads

        
    def run(self):
        for this_srcwin in utils.yield_srcwin(
                n_size=(self.mod.ycount, self.mod.xcount),
                n_chunk=self.mod.chunk_size,
                step=self.mod.chunk_size,
                verbose=True
        ):
            self.scipy_q.put(this_srcwin)

        self.pbar = tqdm(
            desc=f'gridding data to {self.mod}',
            total=self.scipy_q.qsize()
        )
        for _ in range(1):
            tg = threading.Thread(
                target=write_array_queue,
                args=(self.grid_q, self.scipy_q, self.mod)
            )
            tg.daemon = True
            tg.start()

        for _ in range(self.n_threads):
            t = threading.Thread(
                target=scipy_queue,
                args=(self.scipy_q, self.grid_q, self.mod, self.pbar)
            )
            t.daemon = True
            t.start()

        self.grid_q.join()
        self.scipy_q.join()
        self.pbar.close()

        
class WafflesSciPy(Waffle):
    """Generate DEM using Scipy gridding interpolation
    
    Generate a DEM using Scipy's gridding interpolation
    Optional gridding methods are 'linear', 'cubic' and 'nearest'
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.griddata.html

    -----------
    Parameters:
    
    method=[linear/cubic/nearest] - interpolation method to use
    chunk_size=[val] - size of chunks in pixels
    chunk_buffer=[val] - size of the chunk buffer in pixels
    chunk_step=[val] - size of the chunk step in pixels
    num_threads=[val] - number of threads to use in interpolation, 
                        use None to process as a single thread

    < scipy:method=<method>:chunk_size=None:chunk_buffer=40 >
    """
    
    def __init__(
            self,
            method='linear',
            chunk_size=None,
            chunk_buffer=20,
            chunk_step=None,
            num_threads=None,
            **kwargs
    ):
        """generate a `scipy` dem"""
        
        super().__init__(**kwargs)
        self.methods = ['linear', 'cubic', 'nearest']
        self.method = method
        self.chunk_size = utils.int_or(chunk_size)
        self.chunk_step = utils.int_or(chunk_step)
        self.chunk_buffer = utils.int_or(chunk_buffer)
        self.num_threads = utils.int_or(num_threads)

        
    ## this 'run' command runs the scipy module in a multiple threads
    def _run(self):
        if self.num_threads is None:
            return(self._run())
            
        self.open()
        try:
            gs = grid_scipy(self, n_threads=self.num_threads)
            gs.daemon = True
    
            gs.start()
            gs.join()
        except (KeyboardInterrupt, SystemExit):
            utils.echo_error_msg(
                'user breakage...please wait while fetches exits.'
            )
            stop_threads = True
            while not gs.scipy_q.empty():
                try:
                    gs.scipy_q.get(False)
                except Empty:
                    continue
                        
                gs.scipy_q.task_done()
            
        self.close()
        
        return(self)

    
    def open(self):
        if self.method not in self.methods:
            utils.echo_error_msg(
                (f'{self.method} is not a valid interpolation method, '
                 f'options are {self.methods}')
            )
            return(self)
        
        if self.chunk_size is None:
            self.chunk_size = int(self.ds_config['nx'] * .1)
            self.chunk_size = 10 \
                if self.chunk_size < 10 \
                   else self.chunk_size
            
        if self.chunk_step is None:
            self.chunk_step = int(self.chunk_size/2)

        self.stack_ds = gdal.Open(self.stack)
        self.points_band = self.stack_ds.GetRasterBand(1)
        self.points_no_data = self.points_band.GetNoDataValue()        
        self.interp_ds = self.stack_ds.GetDriver().Create(
            self.fn,
            self.stack_ds.RasterXSize,
            self.stack_ds.RasterYSize,
            bands=1,
            eType=self.points_band.DataType,
            options=["BLOCKXSIZE=256",
                     "BLOCKYSIZE=256",
                     "TILED=YES",
                     "COMPRESS=LZW",
                     "BIGTIFF=YES"]
        )
        if self.interp_ds is not None:
            self.interp_ds.SetProjection(self.stack_ds.GetProjection())
            self.interp_ds.SetGeoTransform(self.stack_ds.GetGeoTransform())
            self.interp_band = self.interp_ds.GetRasterBand(1)
            self.interp_band.SetNoDataValue(np.nan)
        else:
            utils.echo_error_msg(
                f'could not create {self.fn}...'
            )
            return(self)
        
        if self.verbose:
            utils.echo_msg(
                f'buffering srcwin by {self.chunk_buffer} pixels'
            )

        self.points_array = self.points_band.ReadAsArray()
        self.stack_ds = None

        
    def close(self):
        self.interp_ds = self.stack_ds = None            

        
    def grid_srcwin(self, srcwin):
        srcwin_buff = utils.buffer_srcwin(
            srcwin, (self.ycount, self.xcount), self.chunk_buffer
        )
        points_array = self.points_array[
            srcwin_buff[1]:srcwin_buff[1]+srcwin_buff[3],
            srcwin_buff[0]:srcwin_buff[0]+srcwin_buff[2]
        ]
        points_array[points_array == self.points_no_data] = np.nan
        point_indices = np.nonzero(~np.isnan(points_array))
        if np.count_nonzero(np.isnan(points_array)) == 0:
            y_origin = srcwin[1]-srcwin_buff[1]
            x_origin = srcwin[0]-srcwin_buff[0]
            y_size = y_origin + srcwin[3]
            x_size = x_origin + srcwin[2]
            points_array = points_array[y_origin:y_size,x_origin:x_size]
            return(points_array)

        elif len(point_indices[0]):
            point_values = points_array[point_indices]
            xi, yi = np.mgrid[0:srcwin_buff[2],
                              0:srcwin_buff[3]]
            #try:
            interp_data = interpolate.griddata(
                np.transpose(point_indices), point_values,
                (xi, yi), method=self.method
            )
            # while np.any(interp_data[np.isnan(interp_data)]):
            #     utils.echo_msg('nodata in {}'.format(srcwin))
            #     point_indices = np.nonzero(~np.isnan(interp_data))
            #     point_values = interp_data[point_indices]
            #     interp_data = interpolate.griddata(
            #         np.transpose(point_indices), point_values,
            #         (xi, yi), method=self.method
            #     )

            y_origin = srcwin[1]-srcwin_buff[1]
            x_origin = srcwin[0]-srcwin_buff[0]
            y_size = y_origin + srcwin[3]
            x_size = x_origin + srcwin[2]
            interp_data = interp_data[y_origin:y_size,x_origin:x_size]
            return(interp_data)
            #except Exception as e:
            #    return(points_array)
                
        return(None)

    
    ## this 'run' command runs the scipy module in a single thread
    def run(self):
        if self.method not in self.methods:
            utils.echo_error_msg(
                (f'{self.method} is not a valid interpolation method, '
                 f'options are {self.methods}')
            )
            return(self)
        
        if self.chunk_size is None:
            n_chunk = int(self.ds_config['nx'] * .1)
            n_chunk = 10 if n_chunk < 10 else n_chunk
        else:
            n_chunk = self.chunk_size
            
        if self.chunk_step is None:
            n_step = int(n_chunk/2)
            #n_step = n_chunk
        else:
            n_step = self.chunk_step

        stack_ds = gdal.Open(self.stack)
        points_band = stack_ds.GetRasterBand(1)
        points_no_data = points_band.GetNoDataValue()
        interp_ds = stack_ds.GetDriver().Create(
            self.fn,
            stack_ds.RasterXSize,
            stack_ds.RasterYSize,
            bands=1,
            eType=points_band.DataType,
            options=["BLOCKXSIZE=256",
                     "BLOCKYSIZE=256",
                     "TILED=YES",
                     "COMPRESS=LZW",
                     "BIGTIFF=YES"]
        )
        if interp_ds is not None:
            interp_ds.SetProjection(stack_ds.GetProjection())
            interp_ds.SetGeoTransform(stack_ds.GetGeoTransform())
            interp_band = interp_ds.GetRasterBand(1)
            interp_band.SetNoDataValue(np.nan)
        else:
            utils.echo_error_msg(
                f'could not create {self.fn}...'
            )
            return(self)
        
        if self.verbose:
            utils.echo_msg(
                f'buffering srcwin by {self.chunk_buffer} pixels'
            )
       
        for srcwin in utils.yield_srcwin(
                (self.ycount, self.xcount),
                n_chunk=n_chunk,
                msg='generating {} grid'.format(self.method),
                verbose=self.verbose,
                step=n_step
        ):
            chunk_buffer = self.chunk_buffer
            srcwin_buff = utils.buffer_srcwin(
                srcwin, (self.ycount, self.xcount), chunk_buffer
            )
            points_array = points_band.ReadAsArray(*srcwin_buff)
            points_array[points_array == points_no_data] = np.nan
            point_indices = np.nonzero(~np.isnan(points_array))
            if np.count_nonzero(np.isnan(points_array)) == 0:
                y_origin = srcwin[1]-srcwin_buff[1]
                x_origin = srcwin[0]-srcwin_buff[0]
                y_size = y_origin + srcwin[3]
                x_size = x_origin + srcwin[2]
                points_array = points_array[y_origin:y_size,x_origin:x_size]
                interp_band.WriteArray(points_array, srcwin[0], srcwin[1])

            elif len(point_indices[0]):
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
                    continue
                
        interp_ds = stack_ds = point_values = weight_values = None
        
        return(self)

    
class WafflesLinear(WafflesSciPy):
    """LINEAR (triangulated) DEM
    
    -----------
    Parameters:
    
    chunk_size (int): size of chunks in pixels
    chunk_step (iint):  size of chunks in pixels
    chunk_buffer (int):  size of the chunk buffer in pixels
    """
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.method = 'linear'

        
class WafflesCubic(WafflesSciPy):
    """CUBIC (triangulated) DEM
    
    -----------
    Parameters:
    
    chunk_size=[val] - size of chunks in pixels
    chunk_buffer=[val] - size of the chunk buffer in pixels
    """
        
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.method = 'cubic'

        
class WafflesNearest(WafflesSciPy):
    """NEAREST neighbor DEM
    
    -----------
    Parameters:
    
    chunk_size=[val] - size of chunks in pixels
    chunk_buffer=[val] - size of the chunk buffer in pixels
    """
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.method = 'nearest'

        
class GMTSurface(Waffle):
    """SPLINE DEM via GMT surface
    
    Generate a DEM using GMT's surface command
    see gmt surface --help for more info.

    -----------
    Parameters:
   
    tension=[0-1] - spline tension.
    relaxation=[val] - spline relaxation factor.
    aspect=[val/None] - gridding aspect
    breakline=[path/None] - use xyz dataset at `path` as a breakline
    convergence=[val/None] - gridding convergence
    blockmean=[True/False] - pipe the data through gmt blockmean before gridding
    geographic=[True/Faslse] - data/grid are geographic

    < gmt-surface:tension=.35:relaxation=1:max_radius=None:aspect=None:breakline=None:convergence=None:blockmean=False:geographic=True >
    """
    
    def __init__(self, tension=1, relaxation=1, max_radius=None,
                 aspect=None, breakline=None, convergence=None,
                 blockmean=False, geographic=True, **kwargs):
        super().__init__(**kwargs)
        if utils.float_or(tension) is not None:
            if utils.float_or(tension) > 1:
                self.tension = 1
            else:
                self.tension = tension
        else:
            self.tension = .35

        self.convergence = utils.float_or(convergence)
        self.relaxation = relaxation
        self.breakline = breakline
        self.max_radius = max_radius
        self.aspect = aspect
        self.blockmean = blockmean
        self.geographic = geographic

        
    def run(self):
        if self.gc['GMT'] is None:
            utils.echo_error_msg(
                'GMT must be installed to use the SURFACE module'
            )
            return(None, -1)

        out, status = utils.run_cmd(
            'gmt gmtset IO_COL_SEPARATOR = SPACE',
            verbose = False,
            cwd = self.cache_dir
        )        

        #self.gmt_region = self.ps_region.copy()
        dem_surf_cmd = ('')
        if self.blockmean:
            dem_surf_cmd = (
                'gmt blockmean {} -I{:.16f}/{:.16f}+e{}{} -V |'.format(
                    self.ps_region.format('gmt'), self.xinc, self.yinc,
                    ' -W' if self.want_weight else '', ' -fg' if self.geographic else ''
                )
            )

        ## mrl: removed -rp and switched p_region to ps_region
        ## (pre 6.5.0 will shift the grid otherwise)
        dem_surf_cmd += (
            'gmt surface -V {} -I{:.16f}/{:.16f}+e -G{}.tif=gd+n{}:GTiff -T{} -Z{} {}{}{}{}{}{}'.format(
                self.ps_region.format('gmt'), self.xinc, self.yinc,
                self.name, self.ndv, self.tension, self.relaxation,
                ' -Qr' if self.gc['GMT'] >= '6.5.0' else '',
                ' -D{}'.format(self.breakline) if self.breakline is not None else '',
                ' -M{}'.format(self.max_radius) if self.max_radius is not None else '',
                ' -C{}'.format(self.convergence) if self.convergence is not None else '',
                ' -A{}'.format(self.aspect) if self.aspect is not None else '',
                ' -fg' if self.geographic else '',
            )
        )

        # dem_surf_cmd += (
        #     'gmt surface -V {} -I{:.16f}/{:.16f} -G{}.tif=gd+n{}:GTiff -rp -T{} -Z{} {}{}{}{}{}{}'.format(
        #         self.p_region.format('gmt'), self.xinc, self.yinc,
        #         self.name, self.ndv, self.tension, self.relaxation,
        #         '' if self.gc['GMT'] >= '6.5.0' else '',
        #         ' -D{}'.format(self.breakline) if self.breakline is not None else '',
        #         ' -M{}'.format(self.max_radius) if self.max_radius is not None else '',
        #         ' -C{}'.format(self.convergence) if self.convergence is not None else '',
        #         ' -A{}'.format(self.aspect) if self.aspect is not None else '',
        #         ' -fg' if self.geographic else '',
        #     )
        # )

        out, status = utils.run_cmd(
            dem_surf_cmd,
            verbose=self.verbose,
            data_fun=lambda p: self.dump_xyz(
                dst_port=p, encode=True
            )
        )
            
        return(self)

    
class GMTTriangulate(Waffle):
    """TRIANGULATION DEM via GMT triangulate
    
    Generate a DEM using GMT's triangulate command.
    see gmt triangulate --help for more info.

    < gmt-triangulate >
    """
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)        
        
    def run(self):
        if self.gc['GMT'] is None:
            utils.echo_error_msg(
                'GMT must be installed to use the TRIANGULATE module'
            )
            return(None, -1)

        out, status = utils.run_cmd(
            'gmt gmtset IO_COL_SEPARATOR = SPACE',
            verbose = False
        )        
        
        dem_tri_cmd = 'gmt triangulate -V {} -I{:.14f}/{:.14f} -G{}.tif=gd:GTiff'.format(
            self.ps_region.format('gmt'), self.xinc, self.yinc, self.name
        )
        out, status = utils.run_cmd(
            dem_tri_cmd,
            verbose = self.verbose,
            data_fun = lambda p: self.dump_xyz(
                dst_port=p, encode=True
            )
        )
        
        return(self)

    
class GMTNearNeighbor(Waffle):
    """NEARNEIGHBOR DEM via GMT nearneighbor
    
    Generate a DEM using GMT's nearneighbor command.
    see gmt nearneighbor --help for more info.
    
    -----------
    Parameters:
    
    radius=[val] - search radius
    sectors=[val] - sector information
    
    < gmt-nearneighbor:radius=None:sectors=None >
    """
    
    def __init__(self, radius=None, sectors=None, **kwargs):
        super().__init__(**kwargs) 
        self.radius = radius
        self.sectors = sectors

        
    def run(self):
        if self.gc['GMT'] is None:
            utils.echo_error_msg(
                'GMT must be installed to use the NEARNEIGHBOR module'
            )
            return(None, -1)

        out, status = utils.run_cmd(
            'gmt gmtset IO_COL_SEPARATOR = SPACE',
            verbose = False
        )        

        dem_nn_cmd = 'gmt nearneighbor -V {} -I{:.14f}/{:.14f} -G{}.tif=gd+n{}:GTiff{}{}{}'.format(
            self.ps_region.format('gmt'),
            self.xinc,
            self.yinc,
            self.name,
            self.ndv,
            ' -W' \
            if self.want_weight \
            else '', ' -N{}'.format(self.sectors) \
            if self.sectors is not None \
            else '',
            ' -S{}'.format(self.radius) \
            if self.radius is not None \
            else ' -S{}'.format(self.xinc),
        )
        out, status = utils.run_cmd(
            dem_nn_cmd,
            verbose = self.verbose,
            data_fun = lambda p: self.dump_xyz(
                dst_port=p, encode=True
            )
        )
        
        return(self)

    
class WafflesMBGrid(Waffle):
    """SPLINE DEM via MB-System's mbgrid.
    
    Generate a DEM using MB-System's mbgrid command.
    By default, will use MB-Systems datalist processes.
    set `use_datalist=True` to use CUDEM's dlim instead.
    see mbgrid --help for more info

    -----------
    Parameters:
    
    dist=[val] - the dist variable to use in mbgrid
    tension=[val] - the spline tension value (0-inf)
    use_stack=[True/False] - use built-in datalists rather than 
                             mbdatalist
    
    < mbgrid:dist='10/3':tension=35:use_datalists=False >
    """
    
    def __init__(
            self,
            dist='10/3',
            tension=0,
            use_stack=True,
            nc=False,
            **kwargs
    ):
        super().__init__(**kwargs)
        self.nc = nc
        self.dist = dist
        self.tension = tension
        self.use_stack = use_stack

        
    def _gmt_num_msk(self, num_grd, dst_msk):
        """generate a num-msk from a NUM grid using GMT grdmath

        Args:
          num_grd (str): pathname to a source `num` grid file
          dst_msk (str): pathname to a destination `msk` grid file

        Returns:
          list: [cmd-output, cmd-return-code]
        """

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
            os.path.basename(src_grd).split('.')[0], gdalfun.gdal_fext(dst_fmt)
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
            os.path.basename(src_grd).split('.')[0], galfun.gdal_fext(dst_fmt)
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

        
    def stack2mbdatalist(self):
        mb_datalist_fn = os.path.join(self.cache_dir, '_tmp_mb.datalist')
        mb_stack_xyz = os.path.join(
            self.cache_dir, '{}.xyz'.format(utils.fn_basename2(self.stack))
        )
        with open(mb_datalist_fn, 'w') as dl:
            dl.write('{} 168 1'.format(mb_stack_xyz))
                
        with open(mb_stack_xyz, 'w') as stack_xyz:
            self.dump_xyz(stack_xyz)

        return(mb_datalist_fn)

    
    def run(self):
        mb_datalist = self.stack2mbdatalist() if self.use_stack else self.data.fn
        #self.mb_region = self.ps_region.copy()
        self.mb_region = self.ps_region.copy()
        #out_name = os.path.join(self.cache_dir, self.name)
        out_name = self.name
        
        mb_xcount, mb_ycount, mb_gt \
            = self.ps_region.geo_transform(
                x_inc=self.xinc, y_inc=self.yinc, node=self.node
            )
        mbgrid_cmd = 'mbgrid -I{} {} -D{}/{} -O{} -A2 -F1 -N -C{} -S0 -X0 -T{}'.format(
            mb_datalist,
            self.mb_region.format('gmt'),
            mb_xcount,
            mb_ycount,
            out_name,
            self.dist,
            self.tension
        )
        out, status = utils.run_cmd(mbgrid_cmd, verbose=self.verbose)
        if status == 0:
            out = gdalfun.gdal2gdal(
                '{}.grd'.format(out_name),
                dst_dem=self.fn,
                dst_fmt=self.fmt,
                co=True
            )
                    
        return(self)

    
class WafflesGDALGrid(Waffle):
    """Waffles GDAL_GRID module.

    see gdal_grid for more info and gridding algorithms
    """
    
    def __init__(self, **kwargs):
        """run gdal grid using alg_str

        parse the data through xyzfun.xyz_block to get weighted 
        mean before building the GDAL dataset to pass into gdal_grid
        
        Args: 
          alg_str (str): the gdal_grid algorithm string
        """
        
        super().__init__(**kwargs)
        self.alg_str = 'linear:radius=-1'
        self.mod = self.alg_str.split(':')[0]

        
    def _vectorize_stack(self):
        """Make a point vector OGR DataSet Object from src_xyz

        for use in gdal gridding functions.
        """

        dst_ogr = '{}'.format(self.name)
        ogr_ds = gdal.GetDriverByName('Memory').Create(
            '', 0, 0, 0, gdal.GDT_Unknown
        )
        layer = ogr_ds.CreateLayer(
            dst_ogr, geom_type=ogr.wkbPoint25D
        )
        for x in ['long', 'lat', 'elev', 'weight']:
            fd = ogr.FieldDefn(x, ogr.OFTReal)
            fd.SetWidth(12)
            fd.SetPrecision(8)
            layer.CreateField(fd)
            
        f = ogr.Feature(feature_def=layer.GetLayerDefn())        
        with tqdm(desc='vectorizing stack', leave=self.verbose) as pbar:
            for this_xyz in self.stack_ds.yield_xyz():
                pbar.update()
                f.SetField(0, this_xyz.x)
                f.SetField(1, this_xyz.y)
                f.SetField(2, float(this_xyz.z))
                if self.want_weight:
                    f.SetField(3, this_xyz.w)

                wkt = this_xyz.export_as_wkt(include_z=True)
                g = ogr.CreateGeometryFromWkt(wkt)
                f.SetGeometryDirectly(g)
                layer.CreateFeature(f)
            
        return(ogr_ds)

    
    def run(self):
        if self.verbose:
            utils.echo_msg(
                'running GDAL GRID {} algorithm @ {} and {}/{}...'.format(
                    self.alg_str.split(':')[0],
                    self.p_region.format('fn'),
                    self.xcount,
                    self.ycount
                )
            )
        _prog = tqdm(
            desc=f'running GDAL GRID {self.alg_str} algorithm',
            leave=self.verbose
        )
        _prog_update = lambda x, y, z: _prog.update()
        ogr_ds = self._vectorize_stack()
        if ogr_ds.GetLayer().GetFeatureCount() == 0:
            utils.echo_error_msg('no input data')
            
        gd_opts = gdal.GridOptions(
            outputType = gdal.GDT_Float32,
            noData = self.ndv,
            format = 'GTiff',
            width = self.xcount,
            height = self.ycount,
            algorithm = self.alg_str,
            callback = _prog_update if self.verbose else None,
            outputBounds = [
                self.p_region.xmin, self.p_region.ymax,
                self.p_region.xmax, self.p_region.ymin
            ]
        )
        gdal.Grid(
            '{}.tif'.format(self.name),
            ogr_ds,
            options=gd_opts
        )
        gdalfun.gdal_set_ndv(
            '{}.tif'.format(
                self.name, ndv=self.ndv, convert_array=False
            )
        )
        ogr_ds = None
        return(self)

    
class GDALLinear(WafflesGDALGrid):
    """LINEAR DEM via gdal_grid
    Generate a DEM using GDAL's gdal_grid command.
    see gdal_grid --help for more info

    -----------
    Parameters:
    
    radius=[val] - search radius
    
    < gdal-linear:radius=-1 >
    """
    
    def __init__(self, radius=None, nodata=-9999, **kwargs):
        super().__init__(**kwargs)        
        radius = self.xinc * 4 if radius is None else utils.str2inc(radius)
        self.alg_str = 'linear:radius={}:nodata={}'.format(radius, nodata)

        
class GDALInvDst(WafflesGDALGrid):
    """INVERSE DISTANCE DEM via gdal_grid
    
    Generate a DEM using GDAL's gdal_grid command.
    see gdal_grid --help for more info

    -----------
    Parameters:
    
    radius1=[val] - search radius 1
    radius2=[val] - search radius 2
    power=[val] - weight**power
    min_points=[val] - minimum points per IDW bucket (use to fill entire DEM)

    < gdal-invdst:power=2.0:smoothing=0.0:radius1=0:radius2=0:angle=0:max_points=0:min_points=0:nodata=0 >
    """
    
    def __init__(
            self,
            power=2.0,
            smoothing=0.0,
            radius1=None,
            radius2=None,
            angle=0.0,
            max_points=0,
            min_points=0,
            nodata=-9999,
            **kwargs
    ):
        super().__init__(**kwargs)
        radius1 = self.xinc * 2 if radius1 is None else utils.str2inc(radius1)
        radius2 = self.yinc * 2 if radius2 is None else utils.str2inc(radius2)
        self.alg_str = (f'invdist:power={power}'
                        f':smoothing={smoothing}'
                        f':radius1={radius1}'
                        f':radius2={radius2}'
                        f':angle={angle}'
                        f':max_points={max_points}'
                        f':min_points={min_points}'
                        f':nodata={nodata}')

        
class GDALMovingAverage(WafflesGDALGrid):
    """Moving AVERAGE DEM via gdal_grid
    
    Generate a DEM using GDAL's gdal_grid command.
    see gdal_grid --help for more info
    
    -----------
    Parameters:
    
    radius1=[val] - search radius 1
    radius2=[val] - search radius 2
    min_points=[val] - minimum points per bucket 
                       (use to fill entire DEM)

    < gdal-average:radius1=0:radius2=0:angle=0:min_points=0:nodata=0 >
    """
    
    def __init__(
            self,
            radius1=None,
            radius2=None,
            angle=0.0,
            min_points=0,
            nodata=-9999,
            **kwargs
    ):
        super().__init__(**kwargs)
        radius1 = self.xinc * 2 if radius1 is None else utils.str2inc(radius1)
        radius2 = self.yinc * 2 if radius2 is None else utils.str2inc(radius2)
        self.alg_str = (f'average:radius1={radius1}'
                        f':radius2={radius2}'
                        f':angle={angle}'
                        f':min_points={min_points}'
                        f':nodata={nodata}')

        
class GDALNearest(WafflesGDALGrid):
    """NEAREST DEM via gdal_grid
    
    Generate a DEM using GDAL's gdal_grid command.
    see gdal_grid --help for more info

    -----------
    Parameters:
    
    radius1=[val] - search radius 1
    radius2=[val] - search radius 2
    angle=[val] - angle
    nodata=[val] - nodata

    < gdal-nearest:radius1=0:radius2=0:angle=0:nodata=0 >
    """
    
    def __init__(
            self,
            radius1=None,
            radius2=None,
            angle=0.0,
            nodata=-9999,
            **kwargs
    ):
        super().__init__(**kwargs)
        radius1 = self.xinc * 2 if radius1 is None else utils.str2inc(radius1)
        radius2 = self.yinc * 2 if radius2 is None else utils.str2inc(radius2)
        self.alg_str = (f'nearest:radius1={radius1}'
                        f':radius2={radius2}'
                        f':angle={angle}'
                        f':nodata={nodata}')
        
    
class WafflesNum(Waffle):
    """Uninterpolated DEM populated by <mode>.
    
    Generate an uninterpolated DEM using <mode> option.
    Using mode of 'A<mode>' uses GMT's xyz2grd command, 
    see gmt xyz2grd --help for more info.
    
    mode keys: k (mask), m (mean), n (num), w (wet), A<mode> (gmt xyz2grd)

    This module is depreciated!

    -----------
    Parameters:
    
    mode=[key] - specify mode of grid population
    min_count=[val] - minimum number of points per cell

    < num:mode=n:min_count=None >
    """
    
    def __init__(self, mode='n', min_count=None, **kwargs):
        """generate an uninterpolated Grid
        `mode` of `n` generates a num grid
        `mode` of `m` generates a mean grid
        `mode` of `k` generates a mask grid
        `mode` of `w` generates a wet/dry mask grid
        `mode` of `A` generates via GMT 'xyz2grd'
        """
        
        super().__init__(**kwargs)        
        self.mode = mode
        self.min_count = min_count
        
    def _xyz_num(self):
        """Create a GDAL supported grid from xyz data """

        xcount, ycount, dst_gt = self.p_region.geo_transform(
            x_inc=self.xinc, y_inc=self.yinc
        )
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
        ds_config = gdalfun.gdal_set_infos(
            xcount,
            ycount,
            xcount * ycount,
            dst_gt,
            gdalfun.osr_wkt(self.dst_srs),
            gdt,
            self.ndv,
            'GTiff',
            None,
            None
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
        gdalfun.gdal_write(out_array, self.fn, ds_config)

        
    def run(self):
        if self.mode.startswith('A'):
            if self.gc['GMT'] is None:
                utils.echo_error_msg(
                    'GMT must be installed to use the Mode `A` with the NUM module'
                )
                return(None, -1)

            out, status = utils.run_cmd(
                'gmt gmtset IO_COL_SEPARATOR = SPACE', verbose = False
            ) 
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

    
class WafflesVDatum(Waffle):
    """VDATUM transformation grid.
    Generate a Vertical DATUM transformation grid.

    -----------
    Parameters:

    mode=[waffle-module] - the waffles module to use to interpolate/extrapolate 
    vdatum_in=[vdatum] - input vertical datum
    vdatum_out=[vdatum] - output vertical datum
    
    < mode=gmt-surface:vdatum:vdatum_in=None:vdatum_out=None >
    """

    def __init__(self, mode='IDW', vdatum_in=None, vdatum_out=None, **kwargs):
        self.valid_modes = [
            'gmt-surface',
            'IDW',
            'linear',
            'cubic',
            'nearest',
            'gmt-triangulate',
            'mbgrid'
        ]
        self.mode = mode
        self.mode_args = {}
        if self.mode not in self.valid_modes:
            self.mode = 'IDW'
            
        tmp_waffles = Waffle()
        tmp_waffles_mode = WaffleFactory(mod=self.mode)._acquire_module()
        for kpam, kval in kwargs.items():
            if kpam not in tmp_waffles.__dict__:
                if kpam in tmp_waffles_mode.__dict__:
                    self.mode_args[kpam] = kval

        for kpam, kval in self.mode_args.items():
            del kwargs[kpam]

        super().__init__(**kwargs)
        self.vdatum_in = vdatum_in
        self.vdatum_out = vdatum_out        

        
    def run(self):
        status = vdatums.VerticalTransform(
            self.mode,
            self.p_region,
            self.xinc,
            self.yinc,
            self.vdatum_in,
            self.vdatum_out,
            node=self.node,
            cache_dir=waffles_cache
        ).run(outfile='{}.tif'.format(self.name))

        utils.echo_msg(status)
        return(self)

    
class WafflesCoastline(Waffle):
    """COASTLINE (land/etc-mask) generation
    
    Generate a coastline (land/etc-mask) using a variety of sources. 
    User data can be provided to provide source for further land masking. 
    Output raster will mask land-areas as 1 and oceans/(lakes/buildings) as 0.
    Output vector will polygonize land-areas.
    
    -----------
    Parameters:
    
    want_osm_coast=[True/False] - use OSM to fill background
    want_enc_coast=[True/False] - use ENC charts to fill background
    want_gmrt=[True/False] - use GMRT to fill background (will use Copernicus otherwise)
    want_copernicus=[True/False] - use COPERNICUS to fill background
    want_cudem=[True/False] - use CUDEM to fill background
    want_nhd=[True/False] - use high-resolution NHD to fill US coastal zones
    want_lakes=[True/False] - mask LAKES using HYDROLAKES
    invert_lakes=[True/False] - invert the lake mask 
                                (invert=True to remove lakes from the waterbodies)
    want_osm_buildings=[True/False] - mask BUILDINGS using OSM
    want_bing_buildings=[True/False] - mask BUILDINGS using BING
    want_wsf_buildings=[True/False] - mask BUILDINGS using WSF
    osm_tries=[val] - OSM max server attempts
    min_building_length=[val] - only use buildings larger than val
    invert=[True/False] - invert the output results
    polygonize=[True/False] - polygonize the output
    min_weight=[val] - weight applied to fetched coastal data

    < coastline:want_gmrt=False:want_nhd=True:want_lakes=False:want_buildings=False:invert=False:polygonize=True >
    """
    
    def __init__(
            self,
            want_osm_coast=True,
            want_enc_coast=False,
            want_nhd=True,
            want_nhd_plus=False,
            want_copernicus=False,
            want_cudem=False,
            want_gmrt=False,
            want_lakes=False,
            invert_lakes=False,
            want_buildings=False,
            min_building_length=None,
            want_osm_planet=False,
            want_bing_buildings=False,
            want_osm_buildings=False,
            want_wsf_buildings=False,
            invert=False,
            polygonize=True, # float(True) is 1.0
            osm_tries=5,
            want_wsf=False,
            min_weight=1,
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

        super().__init__(**kwargs)
        self.want_osm_coast = want_osm_coast
        self.want_enc_coast = want_enc_coast
        self.want_nhd = want_nhd
        self.want_nhd_plus = want_nhd_plus
        self.want_gmrt = want_gmrt
        self.want_copernicus = want_copernicus
        self.want_cudem = want_cudem
        self.want_lakes = want_lakes
        self.invert_lakes = invert_lakes
        self.want_buildings = want_buildings
        self.want_wsf = want_wsf
        self.want_bing_buildings = want_bing_buildings
        self.want_osm_buildings = want_osm_buildings
        self.want_wsf_buildings = want_wsf_buildings
        self.min_building_length = min_building_length
        self.want_osm_planet = want_osm_planet
        self.invert = invert
        self.polygonize = polygonize
        self.osm_tries = utils.int_or(osm_tries, 5)
        self.coast_array = None
        self.ds_config = None
        self.min_weight = utils.float_or(min_weight, 1)

        
    def fetch_data(self, fetches_module, check_size=True):
        this_fetches = fetches.FetchesFactory(
            mod=fetches_module,
            src_region=self.wgs_region,
            verbose=self.verbose,
            outdir=self.cache_dir,
            callback=fetches.fetches_callback
        )._acquire_module()        
        this_fetches.run()
        fr = fetches.fetch_results(this_fetches, check_size=check_size)
        fr.daemon = True
        fr.start()
        fr.join()

        return(fr)

    
    def run(self):
        self.f_region = self.p_region.copy()
        self.f_region.buffer(pct=5, x_inc=self.xinc, y_inc=self.yinc)
        self.f_region.src_srs = self.dst_srs
        self.wgs_region = self.f_region.copy()
        self.wgs_srs = 'epsg:4326'
        if self.dst_srs is not None:
            self.wgs_region.warp(self.wgs_srs)
        else:
            self.dst_srs = self.wgs_srs
            
        horz_epsg, vert_epsg = gdalfun.epsg_from_input(self.dst_srs)
        self.cst_srs = horz_epsg

        self.mem_config = self.ds_config
        self.mem_config['dt'] = gdal.GDT_Byte
        
        self._load_coast_mask()

        self.want_nhd = True if self.want_nhd_plus else self.want_nhd

        if self.want_osm_coast:
            self._load_osm_coast()

        if self.want_enc_coast:
            self._load_enc_coast()
            #sys.exit()
            
        if self.want_gmrt:
            self._load_gmrt()
            
        if self.want_copernicus:
            self._load_copernicus()

        if self.want_cudem:
            self._load_cudem()

        if self.want_nhd:
            self._load_nhd()

        if self.want_lakes:
            self._load_lakes()

        if self.want_buildings or self.want_bing_buildings:
            self._load_bing_bldgs()

        if self.want_osm_buildings:
            self._load_osm_bldgs()

        if self.want_wsf or self.want_wsf_buildings:
            self._load_wsf_bldgs()
            
        if self.want_stack:
            self._load_data()
            
        if self.verbose:
            utils.echo_msg(
                'finanlizing array for region {} at {} {}...'.format(
                    self.p_region.format('gmt'),
                    self.ds_config['nx'],
                    self.ds_config['ny']
                )
            )
        self._finalize_array()
        self._write_coast_array()
        if self.polygonize:
            if utils.int_or(self.polygonize) is not None:
                self._write_coast_poly(poly_count=self.polygonize)
            else:
                self._write_coast_poly()

        if np.all(self.coast_array == 0) or np.all(self.coast_array == 1):
            return(None)
        else:
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
        
        xcount, ycount, gt = self.p_region.geo_transform(
            x_inc=self.xinc, y_inc=self.yinc
        )
        self.ds_config = gdalfun.gdal_set_infos(
            xcount,
            ycount,
            xcount * ycount,
            gt,
            self.cst_srs,
            gdal.GDT_Float32,
            self.ndv,
            'GTiff',
            None,
            None
        )        
        self.coast_array = np.zeros( (ycount, xcount) )

        
    def _apply_coast_mask(self, msk_arr, srcwin=None):
        if srcwin is None:
            srcwin = (0,0,self.ds_config['nx'],self.ds_config['ny'])            
            self.coast_array[
                srcwin[1]:srcwin[1]+srcwin[3],
                srcwin[0]:srcwin[0]+srcwin[2]
            ] += (msk_arr * self.min_weight)

            
    def _load_osm_coast(self):
        this_cst = self.fetch_data('osm:q=coastline', check_size=False)
        if this_cst is not None:
            with tqdm(
                    total=len(this_cst.results),
                    desc='processing osm coastline',
                    leave=True
            ) as pbar:
                for n, cst_result in enumerate(this_cst.results):                
                    if cst_result[-1] == 0:
                        pbar.update()
                        cst_osm = cst_result[1]
                        out = fetches.polygonize_osm_coastline(
                            cst_osm, utils.make_temp_fn(
                                utils.fn_basename2(cst_osm) + '_coast.shp',
                                temp_dir=self.cache_dir
                            ),
                            region=self.wgs_region,
                            include_landmask=False,
                            landmask_is_watermask=True,
                            line_buffer=0.0000001
                        )

                        cst_ds = gdalfun.gdal_mem_ds(
                            self.mem_config,
                            name='osm_coast',
                            src_srs=self.wgs_srs,
                            co=self.co
                        )
                        cst_ds_ = ogr.Open(out)
                        if cst_ds_ is not None:
                            cst_layer = cst_ds_.GetLayer()
                            gdal.RasterizeLayer(
                                cst_ds, [1], cst_layer, burn_values=[1]
                            )

                            for srcwin in utils.yield_srcwin(
                                    (self.ds_config['ny'], self.ds_config['nx']),
                                    1000,
                                    verbose=self.verbose
                            ):
                                cst_ds_arr = cst_ds.GetRasterBand(1).ReadAsArray(*srcwin)
                                self.coast_array[
                                    srcwin[1]:srcwin[1]+srcwin[3],
                                    srcwin[0]:srcwin[0]+srcwin[2]
                                ] += cst_ds_arr
                            cst_ds = cst_ds_ = None
                            
                        else:
                            utils.echo_error_msg(
                                f'could not open {out}'
                            )

                            
    def _load_enc_coast(self):
        this_cst = self.fetch_data('charts')
        if this_cst is not None:
            with tqdm(
                    total=len(this_cst.results),
                    desc='processing charts coastline',
                    leave=True
            ) as pbar:
                for n, cst_result in enumerate(this_cst.results):
                    if cst_result[-1] == 0:
                        pbar.update()
                        cst_enc = cst_result[1]
                        src_000s = utils.p_unzip(
                            os.path.join(cst_enc),
                            exts=['000'],
                            outdir=self.cache_dir,
                            verbose=self.verbose            
                        )
                        for src_000 in src_000s:
                            enc_level = utils.int_or(
                                os.path.basename(src_000)[2], 0
                            )
                            if enc_level <= 3:
                                continue

                            enc_ds = gdalfun.gdal_mem_ds(
                                self.mem_config,
                                name='enc',
                                src_srs=self.wgs_srs,
                                co=self.co
                            )
                            cst_ds_ = ogr.Open(src_000)
                            if cst_ds_ is not None:
                                cst_layer = cst_ds_.GetLayer('SEAARE')                            
                                _ = gdal.RasterizeLayer(
                                    enc_ds, [1], cst_layer, burn_values=[1]
                                )
                                cst_ds_ = cst_layer = None
                                
                            if enc_ds is not None:
                                for srcwin in utils.yield_srcwin(
                                        (self.ds_config['ny'], self.ds_config['nx']),
                                        1000,
                                        verbose=self.verbose
                                ):
                                    enc_ds_arr = enc_ds.GetRasterBand(1).ReadAsArray(*srcwin)
                                    self.coast_array[
                                        srcwin[1]:srcwin[1]+srcwin[3],
                                        srcwin[0]:srcwin[0]+srcwin[2]
                                    ] -= (enc_ds_arr * self.min_weight)
                                    
                                enc_ds = enc_ds_arr = None

                                
    def _load_gmrt(self):
        """GMRT - Global low-res.

        Used to fill un-set cells.
        """
        
        this_gmrt = self.fetch_data('gmrt:layer=topo')        
        gmrt_result = this_gmrt.results[0]
        if gmrt_result[-1] == 0:
            gmrt_tif = gmrt_result[1]
            gmrt_ds = gdalfun.gdal_mem_ds(
                self.ds_config,
                name='gmrt',
                src_srs=self.wgs_srs,
                co=self.co
            )
            gdal.Warp(
                gmrt_ds,
                gmrt_tif,
                dstSRS=self.cst_srs,
                resampleAlg=self.sample
            )
            gmrt_ds_arr = gmrt_ds.GetRasterBand(1).ReadAsArray()
            gmrt_ds_arr[gmrt_ds_arr > 0] = 1
            gmrt_ds_arr[gmrt_ds_arr < 0] = 0
            self.coast_array += (gmrt_ds_arr * self.min_weight)
            gmrt_ds = gmrt_ds_arr = None

            
    def _load_copernicus(self):
        """copernicus"""

        this_cop = self.fetch_data(
            'copernicus:datatype=1', check_size=False
        )
        for i, cop_result in enumerate(this_cop.results):
            if cop_result[-1] == 0:
                cop_tif = cop_result[1]
                gdalfun.gdal_set_ndv(cop_tif, 0, verbose=False)
                cop_ds = gdalfun.gdal_mem_ds(
                    self.ds_config,
                    name='copernicus',
                    src_srs=self.wgs_srs,
                    co=self.co
                )
                gdal.Warp(
                    cop_ds,
                    cop_tif,
                    dstSRS=self.cst_srs,
                    resampleAlg=self.sample,
                    callback=False,
                    srcNodata=0
                )
                for srcwin in utils.yield_srcwin(
                        (self.ds_config['ny'], self.ds_config['nx']),
                        1000,
                        verbose=self.verbose
                ):
                    cop_ds_arr = cop_ds.GetRasterBand(1).ReadAsArray(*srcwin)
                    cop_ds_arr[cop_ds_arr != 0] = 1
                    self.coast_array[
                        srcwin[1]:srcwin[1]+srcwin[3],
                        srcwin[0]:srcwin[0]+srcwin[2]
                    ] += (cop_ds_arr * self.min_weight)
                    
                cop_ds = cop_ds_arr = None

                
    def _load_cudem(self):
        """cudem"""

        this_cudem = self.fetch_data(
            'CUDEM:datatype=ninth', check_size=False
        )
        for i, cudem_result in enumerate(this_cudem.results):
            if cudem_result[-1] == 0:
                cudem_tif = cudem_result[1]
                cudem_ds = gdalfun.gdal_mem_ds(
                    self.mem_config,
                    name='cudem',
                    src_srs=self.wgs_srs,
                    co=self.co
                )
                gdal.Warp(
                    cudem_ds, cudem_tif, dstSRS=self.cst_srs, resampleAlg=self.sample,
                    callback=False, srcNodata=0
                )
                for srcwin in utils.yield_srcwin(
                        (self.ds_config['ny'], self.ds_config['nx']),
                        1000,
                        verbose=self.verbose
                ):
                    cudem_ds_arr = cudem_ds.GetRasterBand(1).ReadAsArray(*srcwin)
                    cudem_ds_arr[cudem_ds_arr != 0] = 1
                    self.coast_array[
                        srcwin[1]:srcwin[1]+srcwin[3],
                        srcwin[0]:srcwin[0]+srcwin[2]
                    ] += (cudem_ds_arr * self.min_weight)
                    
                cudem_ds = cudem_ds_arr = None

                
    def _load_nhd(self):
        """USGS NHD (HIGH-RES U.S. Only)
        Fetch NHD (NHD High/Plus) data from TNM to fill in 
        near-shore areas. 
        
        High resoultion data varies by location...
        """

        this_tnm = self.fetch_data(
            "tnm:datasets=13/14:extents='HU-8 Subbasin,HU-4 Subregion'"
        )
        if len(this_tnm.results) > 0:
            tnm_ds = gdalfun.gdal_mem_ds(
                self.ds_config,
                name='nhd',
                src_srs=self.wgs_srs,
                co=self.co
            )
            for i, tnm_result in enumerate(this_tnm.results):
                if tnm_result[-1] == 0:
                    #tnm_zip = os.path.join(this_tnm._outdir, tnm_result[1])
                    tnm_zip = tnm_result[1]
                    if not os.path.exists(tnm_zip):
                        break

                    tnm_zips = utils.unzip(tnm_zip, self.cache_dir, verbose=False)
                    gdb = '.'.join(tnm_zip.split('.')[:-1]) + '.gdb'
                    utils.run_cmd(
                        'ogr2ogr -update -append nhdArea_merge.shp {} NHDArea -where "FType=312 OR FType=336 OR FType=445 OR FType=460 OR FType=537" -clipdst {} 2>/dev/null'.format(
                            gdb, self.p_region.format('ul_lr')
                        ),
                        verbose=self.verbose
                    )
                    utils.run_cmd(
                        'ogr2ogr -update -append nhdArea_merge.shp {} NHDWaterbody -where "FType=493 OR FType=466" -clipdst {} 2>/dev/null'.format(
                            gdb, self.p_region.format('ul_lr')
                        ),
                        verbose=self.verbose
                    )

                    if self.want_nhd_plus:
                        utils.run_cmd(
                            'ogr2ogr -update -append nhdArea_merge.shp {} NHDPlusBurnWaterbody -clipdst {} 2>/dev/null'.format(
                                gdb, self.p_region.format('ul_lr')
                            ),
                            verbose=self.verbose
                        )

                try:
                    utils.run_cmd(
                        'gdal_rasterize -burn 1 nhdArea_merge.shp nhdArea_merge.tif -te {} -ts {} {} -ot Int32'.format(
                            self.p_region.format('te'),
                            self.ds_config['nx'],
                            self.ds_config['ny'],
                        ),
                        verbose=self.verbose
                    )

                    gdal.Warp(
                        tnm_ds,
                        'nhdArea_merge.tif',
                        dstSRS=self.cst_srs,
                        resampleAlg=self.sample
                    )
                except:
                    tnm_ds = None

                if tnm_ds is not None:
                    tnm_ds_arr = tnm_ds.GetRasterBand(1).ReadAsArray()
                    tnm_ds_arr[tnm_ds_arr < 1] = 0
                    self.coast_array -= (tnm_ds_arr * self.min_weight)
                    #self.coast_array[tnm_ds_arr < 1] = 0
                    tnm_ds = tnm_ds_arr = None

                utils.remove_glob('nhdArea_merge.*')

                
    def _load_lakes(self):
        """HydroLakes -- Global Lakes"""
        
        this_lakes = self.fetch_data('hydrolakes')
        lakes_shp = None
        if this_lakes.results[0][-1] == 0:
            lakes_zip = this_lakes.results[0][1]
            lakes_shps = utils.unzip(lakes_zip, self.cache_dir)
            for i in lakes_shps:
                if i.split('.')[-1] == 'shp':
                    lakes_shp = i

            lakes_ds = gdalfun.gdal_mem_ds(
                self.ds_config,
                name='lakes',
                src_srs=self.wgs_srs,
                co=self.co
            )
            lakes_warp_ds = gdalfun.gdal_mem_ds(
                self.ds_config,
                name='lakes',
                src_srs=self.wgs_srs,
                co=self.co
            )
            lk_ds = ogr.Open(lakes_shp)
            if lk_ds is not None:
                lk_layer = lk_ds.GetLayer()
                lk_layer.SetSpatialFilter(self.f_region.export_as_geom())
                gdal.RasterizeLayer(lakes_ds, [1], lk_layer, burn_values=[-1])
                gdal.Warp(
                    lakes_warp_ds,
                    lakes_ds,
                    dstSRS=self.cst_srs,
                    resampleAlg=self.sample
                )
                lakes_ds_arr = lakes_warp_ds.GetRasterBand(1).ReadAsArray()
                # update for weights
                self.coast_array[lakes_ds_arr == -1] = 0 if not self.invert_lakes else 1 
                lakes_ds = lk_ds = lakes_warp_ds = None
            else:
                utils.echo_error_msg(
                    f'could not open {lakes_shp}'
                )

                
    def _load_wsf_bldgs(self):
        """wsf"""

        this_wsf = self.fetch_data('wsf', check_size=False)        
        for i, wsf_result in enumerate(this_wsf.results):
            if wsf_result[-1] == 0:
                #wsf_tif = os.path.join(this_wsf._outdir, wsf_result[1])
                wsf_tif = wsf_result[1]
                gdalfun.gdal_set_ndv(wsf_tif, 0, verbose=False)
                wsf_ds = gdalfun.gdal_mem_ds(
                    self.ds_config,
                    name='wsf',
                    src_srs=self.wgs_srs,
                    co=self.co
                )
                gdal.Warp(
                    wsf_ds,
                    wsf_tif,
                    dstSRS=self.cst_srs,
                    resampleAlg='cubicspline',
                    callback=gdal.TermProgress,
                    srcNodata=0,
                    outputType=gdal.GDT_Float32
                )
                wsf_ds_arr = wsf_ds.GetRasterBand(1).ReadAsArray()
                wsf_ds_arr[wsf_ds_arr != 0 ] = -1
                self.coast_array += (wsf_ds_arr * self.min_weight)
                wsf_ds = wsf_ds_arr = None

                
    def _load_osm_bldgs(self):
        """load buildings from OSM
        
        OSM has a size limit, so will chunk the region and
        do multiple osm calls
        """

        this_osm = self.fetch_data(
            "osm:q=buildings:fmt=osm:chunks=True:min_length={}:planet={}".format(
                self.min_building_length, self.want_osm_planet
            )
        )
        os.environ["OGR_OSM_OPTIONS"] = "INTERLEAVED_READING=YES"
        os.environ["OGR_OSM_OPTIONS"] = "OGR_INTERLEAVED_READING=YES"
        with tqdm(
                total=len(this_osm.results),
                desc='processing OSM buildings',
                leave=self.verbose
        ) as pbar:
            for n, osm_result in enumerate(this_osm.results):
                if osm_result[-1] == 0:
                    pbar.update()
                    #osm_z = os.path.join(this_osm._outdir, osm_result[1])
                    osm_z = osm_result[1]
                    # if fetches.Fetch(osm_result[0], verbose=True).fetch_file(
                    #         osm_z, check_size=False, tries=self.osm_tries, read_timeout=3600
                    # ) >= 0:
                    #if True:
                    if osm_result[-1] == 'bz2':
                        osm_planet = utils.unbz2(osm_z, self.cache_dir)
                        osm_file = utils.ogr_clip(osm_planet, self.wgs_region)
                        _clipped = True
                    elif osm_result[-1] == 'pbf':
                        osm_file = utils.ogr_clip(osm_z, self.wgs_region, 'multipolygons')
                        _clipped = True
                    else:
                        osm_file = osm_z
                        _clipped = False

                    if os.path.getsize(osm_file) == 366:
                        continue

                    out, status = utils.run_cmd(
                        (f'gdal_rasterize -burn -1 -l multipolygons {osm_file} '
                         f'bldg_osm.tif -te {self.p_region.format("te")} '
                         f'-ts {self.ds_config["nx"]} '
                         f'{self.ds_config["ny"]} -ot Int32 -q'),
                        verbose=False
                    )

                    if status == 0:
                        bldg_ds = gdal.Open('bldg_osm.tif')
                        if bldg_ds is not None:
                            bldg_ds_arr = bldg_ds.GetRasterBand(1).ReadAsArray()
                            # update for weights
                            self.coast_array[bldg_ds_arr == -1] = 0 
                            bldg_ds = bldg_ds_arr = None

                        bldg_ds = None
                    else:
                        utils.echo_warning_msg(
                            'could not parse buildings!'
                        )

                    utils.remove_glob('bldg_osm.tif*')

        bldg_ds = bldg_warp_ds = None

        
    def _load_bing_bldgs(self):
        """load buildings from BING      
        """

        this_bing = self.fetch_data("bing_bfp", check_size=True)
        os.environ["OGR_OSM_OPTIONS"] = "INTERLEAVED_READING=YES"
        os.environ["OGR_OSM_OPTIONS"] = "OGR_INTERLEAVED_READING=YES"
        with tqdm(
                total=len(this_bing.results),
                desc='processing BING buildings',
                leave=self.verbose
        ) as pbar:
            for n, bing_result in enumerate(this_bing.results):                
                if bing_result[-1] == 0:
                    pbar.update()

                    bing_gz = bing_result[1]
                    try:
                        bing_gj = utils.gunzip(bing_gz, self.cache_dir)
                        os.rename(bing_gj, bing_gj + '.geojson')
                        bing_gj = bing_gj + '.geojson'

                        out, status = utils.run_cmd(
                            ('gdal_rasterize -burn -1 -l '
                             f'{os.path.basename(utils.fn_basename2(bing_gj))} '
                             f'{bing_gj} bldg_bing.tif '
                             f'-te {self.p_region.format("te")} '
                             f'-ts {self.ds_config["nx"]} '
                             f'{self.ds_config["ny"],} -ot Int32'),
                            verbose=True
                        )
                    except Exception as e:
                        utils.remove_glob(bing_gz)
                        utils.echo_error_msg(
                            f'could not process bing bfp, {e}'
                        )
                        status = -1
                    
                    if status == 0:
                        bldg_ds = gdal.Open('bldg_bing.tif')
                        if bldg_ds is not None:
                            bldg_ds_arr = bldg_ds.GetRasterBand(1).ReadAsArray()
                            # update for weights
                            self.coast_array[bldg_ds_arr == -1] = 0 
                            bldg_ds = bldg_ds_arr = None

                        bldg_ds = None
                    else:
                        utils.echo_warning_msg(
                            'could not parse buildings!'
                        )

                    utils.remove_glob('bldg_bing.tif*')

        bldg_ds = bldg_warp_ds = None

        
    def _load_data(self):
        """load data from user datalist and amend coast_array"""
        for this_arr in self.stack_ds.yield_array():
            data_arr = this_arr[0]['z']
            weight_arr = this_arr[0]['weight']
            weight_arr[np.isnan(weight_arr)] = 0
            srcwin = this_arr[1]
            data_arr[np.isnan(data_arr)] = 0
            data_arr[data_arr > 0] = 1
            data_arr[data_arr < 0] = -1
            self.coast_array[
                srcwin[1]:srcwin[1]+srcwin[3],
                srcwin[0]:srcwin[0]+srcwin[2]
            ] += (data_arr*weight_arr)
            
        # for this_arr in self.stack_ds.yield_array():
        #     data_arr = this_arr[0]['z']
        #     weight_arr = this_arr[0]['weight']
        #     weight_arr[np.isnan(weight_arr)] = 0
        #     srcwin = this_arr[1]
        #     data_arr[np.isnan(data_arr)] = 0
        #     data_mask = data_arr >= 0
        #     data_mask = scipy.ndimage.binary_fill_holes(data_mask)

        #     data_arr = np.ones(data_arr.shape)
        #     #data_arr[data_mask] = 1
        #     data_arr[~data_mask] = -1
        #     #data_arr[data_arr > 0] = 1
        #     #data_arr[data_arr < 0] = -1

        #     # #expanded_arr = scipy.ndimage.binary_dilation(data_arr >= 0, iterations=1, structure=np.ones((12,12)))
        #     # #contracted_arr = scipy.ndimage.binary_erosion(expanded_arr, iterations=12, structure=np.ones((12,12)))
        #     # expanded_arr = utils.expand_for(data_arr > 0, shiftx=12, shifty=12)
        #     # contracted_arr = np.invert(utils.expand_for(np.invert(expanded_arr), shiftx=12, shifty=12))
            
        #     # #data_arr[np.isnan(data_arr)] = 0
        #     # data_arr[contracted_arr] = 1
        #     # data_arr[~contracted_arr] = -1
        #     # weight_arr[contracted_arr & np.isnan(weight_arr)] = 10
        #     # weight_arr[contracted_arr] = 10

        #     # utils.echo_msg_bold(contracted_arr)
        #     # utils.echo_msg_bold(data_arr)
        #     self.coast_array[
        #         srcwin[1]:srcwin[1]+srcwin[3],
        #         srcwin[0]:srcwin[0]+srcwin[2]
        #     ] += (data_arr*weight_arr)

            
    def _write_coast_array(self):
        """write coast_array to file"""

        if self.verbose:
            utils.echo_msg(
                f'writing array to {self.name}.tif...'
            )
            
        gdalfun.gdal_write(
            self.coast_array,
            f'{self.name}.tif',
            self.ds_config,
        )

        
    def _write_coast_poly(self, poly_count=None):
        """convert to coast_array vector"""

        if self.verbose:
            utils.echo_msg(
                f'polygonizing {poly_count} features from array to {self.name}.shp...'
            )
        
        poly_count = utils.int_or(poly_count)
        tmp_ds = ogr.GetDriverByName('ESRI Shapefile').CreateDataSource(
            '{}_tmp_c.shp'.format(self.name)
        )
        
        if tmp_ds is not None:
            tmp_layer = tmp_ds.CreateLayer(
                '{}_tmp_c'.format(os.path.basename(self.name)),
                None,
                ogr.wkbMultiPolygon
            )
            tmp_layer.CreateField(ogr.FieldDefn('DN', ogr.OFTInteger))
            gdalfun.gdal_polygonize(
                '{}.tif'.format(self.name), tmp_layer, verbose=self.verbose
            )
            tmp_ds = None
            
        utils.run_cmd(
            'ogr2ogr -dialect SQLITE -sql "SELECT * FROM {}_tmp_c WHERE DN=0{}" {}.shp {}_tmp_c.shp'.format(
                os.path.basename(self.name),
                ' order by ST_AREA(geometry) desc limit {}'.format(poly_count) \
                if poly_count is not None \
                else '',
                self.name,
                self.name
            ),
            verbose=self.verbose
        )        
        utils.remove_glob(f'{self.name}_tmp_c.*')
        utils.run_cmd(
            'ogrinfo -dialect SQLITE -sql "UPDATE {} SET geometry = ST_MakeValid(geometry)" {}.shp'.format(
                os.path.basename(self.name), self.name),
            verbose=self.verbose
        )        
        gdalfun.osr_prj_file(self.name + '.prj', self.dst_srs)

        
class WafflesLakes(Waffle):
    """Estimate lake bathymetry.
    
    By default, will return lake bathymetry as depth values (positive down), 
    to get elevations (positive up), set apply_elevations=True.

    -----------
    Parameters:
    
    apply_elevations=[True/False] - use COPERNICUS to apply lake level elevations to output
    min_area=[val] - minimum lake area to consider
    max_area=[val] - maximum lake area to consider
    min_id=[val] - minimum lake ID to consider
    max_id=[val] - maximum lake ID to consider
    depth=[globathy/hydrolakes/val] - obtain the depth value from GloBathy, 
                                      HydroLakes or constant value
    
    < lakes:apply_elevations=False:min_area=None:max_area=None:min_id=None:max_id=None:depth=globathy:elevations=copernicus >
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

        #self.initialize()
        #self.init_lakes()

        
    def _fetch_lakes(self):
        """fetch hydrolakes polygons"""

        this_lakes = fetches.HydroLakes(
            src_region=self.p_region,
            verbose=self.verbose,
            outdir=self.cache_dir
        )
        this_lakes.run()
        fr = fetches.fetch_results(this_lakes, check_size=False)
        fr.daemon = True
        fr.start()
        fr.join()

        lakes_shp = None
        lakes_zip = os.path.join(this_lakes._outdir, this_lakes.results[0]['dst_fn'])
        lakes_shps = utils.unzip(lakes_zip, self.cache_dir)
        for i in lakes_shps:
            if i.split('.')[-1] == 'shp':
                lakes_shp = i

        return(lakes_shp)

    
    def _fetch_globathy(self, ids=[]):
        """fetch globathy csv data and process into dict"""
        
        import csv
        
        _globathy_url = 'https://springernature.figshare.com/ndownloader/files/28919991'
        globathy_zip = os.path.join(self.cache_dir, 'globathy_parameters.zip')
        fetches.Fetch(
            _globathy_url, verbose=self.verbose
        ).fetch_file(globathy_zip, check_size=False)
        globathy_csvs = utils.unzip(globathy_zip, self.cache_dir)        
        globathy_csv = os.path.join(
            self.cache_dir,
            'GLOBathy_basic_parameters/GLOBathy_basic_parameters(ALL_LAKES).csv'
        )
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
        
        this_gmrt = fetches.GMRT(
            src_region=gmrt_region,
            verbose=self.verbose,
            layer='topo',
            outdir=self.cache_dir
        )
        this_gmrt.run()

        fr = fetches.fetch_results(this_gmrt)
        fr.daemon = True
        fr.start()
        fr.join()

        dst_srs = osr.SpatialReference()
        dst_srs.SetFromUserInput(self.dst_srs)
        
        gmrt_tif = os.path.join(this_gmrt._outdir, this_gmrt.results[0]['dst_fn'])
        gmrt_ds = gdalfun.gdal_mem_ds(self.ds_config, name='gmrt', co=self.co)
        gdal.Warp(gmrt_ds, gmrt_tif, dstSRS=dst_srs, resampleAlg=self.sample)
        return(gmrt_ds)

    
    def _fetch_copernicus(self, cop_region=None):
        """copernicus"""

        if cop_region is None:
            cop_region = self.p_region
            
        this_cop = fetches.CopernicusDEM(
            src_region=cop_region,
            verbose=self.verbose,
            datatype='1',
            outdir=self.cache_dir
        )
        this_cop.run()

        fr = fetches.fetch_results(this_cop, check_size=False)
        fr.daemon = True
        fr.start()
        fr.join()

        dst_srs = osr.SpatialReference()
        dst_srs.SetFromUserInput(self.dst_srs)
        
        cop_ds = gdalfun.gdal_mem_ds(self.ds_config, name='copernicus', co=self.co)
        [gdal.Warp(
            cop_ds,
            os.path.join(this_cop._outdir, cop_result['dst_fn']),
            dstSRS=dst_srs,
            resampleAlg=self.sample
        ) for cop_result in this_cop.results]
        
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
        with tqdm(
                total=len(nfeatures),
                desc='applying labels.',
                leave=self.verbose
        ) as pbar:
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
            utils.echo_warning_msg(
                'invalid shore array, using default shore value of zero'
            )
            bathy_arr = 0 - bathy_arr
        else:
            bathy_arr = shore_arr - bathy_arr

        utils.echo_msg('applied shore elevations to lake depths')
        bathy_arr[np.isnan(bathy_arr)] = 0
        return(bathy_arr)    

    
    def _init(self):
        self.wgs_region = self.p_region.copy()
        if self.dst_srs is not None:
            self.wgs_region.warp('epsg:4326')            
        else:
            self.dst_srs = 'epsg:4326'
        
        #self.p_region.buffer(pct=2)
        
        lakes_shp = self._fetch_lakes()
        self.lk_ds = ogr.Open(lakes_shp, 1)
        self.lk_layer = self.lk_ds.GetLayer()

        ## filter layer to region
        filter_region = self.p_region.copy()
        self.lk_layer.SetSpatialFilter(filter_region.export_as_geom())

        ## filter by ID
        if self.max_id is not None:
            self.lk_layer.SetAttributeFilter(
                f'Hylak_id < {self.max_id}'
            )
            
        if self.min_id is not None:
            self.lk_layer.SetAttributeFilter(
                f'Hylak_id > {self.min_id}'
            )

        ## filter by Area
        if self.max_area is not None:
            self.lk_layer.SetAttributeFilter(
                f'Lake_area < {self.max_area}'
            )
            
        if self.min_area is not None:
            self.lk_layer.SetAttributeFilter(
                f'Lake_area > {self.min_area}'
            )
            
        lk_features = self.lk_layer.GetFeatureCount()
        if lk_features == 0:
            utils.echo_error_msg('no lakes found in region')
            return(-1)
            #return(self)

        ## get lake ids and globathy depths
        self.lk_ids = []
        [self.lk_ids.append(feat.GetField('Hylak_id')) for feat in self.lk_layer]
        utils.echo_msg('using Lake IDS: {}'.format(self.lk_ids))
        
        lk_regions = self.p_region.copy()
        with tqdm(
                total=len(self.lk_layer),
                desc=f'processing {lk_features} lakes',
                leave=self.verbose
        ) as pbar:            
            for lk_f in self.lk_layer:
                pbar.update()
                this_region = regions.Region()
                lk_geom = lk_f.GetGeometryRef()
                lk_wkt = lk_geom.ExportToWkt()
                this_region.from_list(
                    ogr.CreateGeometryFromWkt(lk_wkt).GetEnvelope()
                )
                lk_regions = regions.regions_merge(lk_regions, this_region)

        while not regions.regions_within_ogr_p(self.p_region, lk_regions):
            utils.echo_msg(
                ('buffering region by 2 percent to gather all '
                 f'lake boundaries...{self.p_region}')
            )
            self.p_region.buffer(pct=2, x_inc=self.xinc, y_inc=self.yinc)

        return(0)

    
    def run(self):
        ## fetch and initialize the shoreline data
        cop_band = None
        cop_arr = None
        ds_config = gdalfun.gdal_infos(self.stack)
        if self.elevations == 'copernicus':
            cop_ds = self._fetch_copernicus(cop_region=self.p_region)
            cop_band = cop_ds.GetRasterBand(1)
        elif self.elevations == 'gmrt':
            cop_ds = self._fetch_gmrt(gmrt_region=self.p_region)
            cop_band = cop_ds.GetRasterBand(1)
        elif utils.float_or(self.elevations) is not None: # single value
            cop_band = None
            cop_arr = np.zeros((ds_config['nx'], ds_config['ny']))
            cop_arr[:] = self.elevations
        elif self.elevations == 'self': # from stacks, interpolated
            cop_band = None
            tmp_arr = gdalfun.gdal_get_array(self.stack, band=1)[0]
            cop_arr = gdalfun.generate_mem_ds(
                gdalfun.gdal_infos(self.stack),
                band_data=tmp_arr,
                srcwin=None,
                return_array=True,
                interpolation='flats'
            )
        elif os.path.exists(self.elevations): # from input raster
            elev_ds = gdal.Open(self.elevations)
            if elev_ds is not None:
                dst_srs = osr.SpatialReference()
                dst_srs.SetFromUserInput(self.dst_srs)
                cop_ds = gdalfun.gdal_mem_ds(
                    ds_config,
                    name='cop',
                    co=self.co
                )
                gdal.Warp(
                    cop_ds,
                    elev_ds,
                    dstSRS=dst_srs,
                    resampleAlg=self.sample
                )
                cop_band = cop_ds.GetRasterBand(1)

        if cop_band is not None:
            cop_arr = cop_band.ReadAsArray()
                
        ## initialize the tmp datasources
        prox_ds = gdalfun.gdal_mem_ds(ds_config, name='prox', co=self.co)
        msk_ds = gdalfun.gdal_mem_ds(ds_config, name='msk', co=self.co)
        msk_band = None
        globd = None
        
        if len(self.lk_ids) == 0:
            return(self)
        
        if self.depth == 'globathy':
            globd = self._fetch_globathy(ids=self.lk_ids[:])
            ## rasterize hydrolakes using id
            gdal.RasterizeLayer(
                msk_ds,
                [1],
                self.lk_layer,
                options=["ATTRIBUTE=Hylak_id"],
                callback=gdal.TermProgress
            )
            msk_ds.FlushCache()
            msk_band = msk_ds.GetRasterBand(1)
            msk_band.SetNoDataValue(ds_config['ndv'])

            ## assign max depth from globathy
            msk_arr = msk_band.ReadAsArray()
            with tqdm(
                    total=len(self.lk_ids),
                    desc='Assigning Globathy Depths to rasterized lakes...',
                    leave=self.verbose
            ) as pbar:
                
                for n, this_id in enumerate(self.lk_ids):
                    depth = globd[this_id]
                    msk_arr[msk_arr == this_id] = depth
                    pbar.update()
            
        elif self.depth == 'hydrolakes':
            gdal.RasterizeLayer(
                msk_ds,
                [1],
                self.lk_layer,
                options=["ATTRIBUTE=Depth_avg"],
                callback=gdal.TermProgress
            )
            msk_ds.FlushCache()
            msk_band = msk_ds.GetRasterBand(1)
            msk_band.SetNoDataValue(ds_config['ndv'])
            msk_arr = msk_band.ReadAsArray()
            
        elif utils.float_or(self.depth) is not None:
            msk_arr = np.zeros((ds_config['nx'], ds_config['ny']))
            msk_arr[:] = self.depth            
        else:            
            msk_arr = np.ones((ds_config['nx'], ds_config['ny']))

        ## calculate proximity of lake cells to shore
        if msk_band is None:
            gdal.RasterizeLayer(
                msk_ds,
                [1],
                self.lk_layer,
                options=["ATTRIBUTE=Hylak_id"],
                callback=gdal.TermProgress
            )
            msk_ds.FlushCache()
            msk_band = msk_ds.GetRasterBand(1)
            msk_band.SetNoDataValue(ds_config['ndv'])
            
        self.lk_ds = None
        prox_band = prox_ds.GetRasterBand(1)
        proximity_options = ["VALUES=0", "DISTUNITS=PIXEL"]
        gdal.ComputeProximity(
            msk_band,
            prox_band,
            options=proximity_options,
            callback=gdal.TermProgress
        )
        prox_arr = prox_band.ReadAsArray()

        ## apply calculation from globathy
        utils.echo_msg('Calculating simulated lake depths...')

        if self.depth == 'flatten':
            cop_arr[prox_arr == 0] = 0
            bathy_arr = cop_arr
        else:
            bathy_arr = self.apply_calculation(
                prox_arr,
                msk_arr,
                shore_arr=cop_arr,
            )
            
        bathy_arr[bathy_arr == 0] = self.ndv            
        gdalfun.gdal_write(
            bathy_arr, '{}.tif'.format(self.name), ds_config,
        )            

        prox_ds = msk_ds = cop_ds = None
        return(self)

    
class WafflesCUDEM(Waffle):
    """CUDEM integrated DEM generation.
    
    Generate an topo/bathy integrated DEM using a variety of data sources.
    Will iterate <pre_count> pre-surfaces at lower-resolutions, 
    specified with <inc_levels>, e.g. .3s/1s

    Each pre-surface will be clipped to <landmask>, if it exists, and smoothed 
    with <pre_smoothing> factor

    Each pre-surface is used in subsequent pre-surface(s)/final DEM at each 
    iterative weight specified in <weight_levels>

    generate a DEM with `pre_surface`s which are generated at lower resolution(s) 
    and with various weight threshholds

    To generate a typical CUDEM tile, generate 1 pre-surface ('bathy_surface'), 
    clipped to a coastline

    Use a <weight_levels>, e.g. 1/.5 that excludes low-resolution bathymetry data 
    from being used as input in the final DEM generation

    e.g.
    cudem:pre_mode=gmt-surface:tension=.65:geographic=True:landmask=True:polygonize=2:weight_levels=1/.5:inc_levels=.3333333s/1s:pre_count=2:pre_upper_limit=-.1:pre_smoothing=2:flatten=95

    -----------
    Parameters:
    
    pre_mode (str) - the waffles module to perform the initial pre-surface
    pre_count (int): number of pre-surface iterations to perform
    weight_levels () - the minimum weights for each of the pre-surface iterations
    inc_levels () -  the increments for each of the pre-surface iterations
    pre_upper_limit (float) - the upper elevation limit of the pre-surfaces 
                              (used with landmask)
    pre_smoothing (float) - the smoothing (blur) factor to apply to the pre-surfaces
    pre_verbose (bool) - increase the verbosity of pre-surface generation
    landmask (bool): path to coastline vector mask or set as `coastline` to auto-generate
    invert_landmask (bool): invert the coastline mask, invert assumes mask over ocean, while
                            non invert assumes mask over land.
    want_supercede (bool) - supercede subsquent pre-surfaces
    flatten (float) - the nodata-size percentile above which to flatten
    exclude_lakes (bool) - exclude lakes from the landmask
    filter_outliers (float) - percentile at which to filter outliers
    <mode-opts> - options for the waffles module specified in 'mode'
    <coastline-opts> - options for the coastline module when generating the landmask
    """

    def __init__(
            self,
            pre_mode='gmt-surface',
            pre_count=2,
            pre_upper_limit=-0.1,
            pre_smoothing=None,
            weight_levels=None,
            inc_levels=None,
            landmask=False,
            invert_landmask=True,
            filter_outliers=None,
            want_supercede=False,
            flatten=None,
            exclude_lakes=False,
            mode=None,
            min_weight='1',
            pre_verbose=False,
            final_mode='IDW',
            keep_pre_surfaces=False,
            **kwargs
    ):
        self.valid_modes = [
            'gmt-surface', 'IDW', 'linear', 'cubic',
            'nearest', 'gmt-triangulate', 'mbgrid'
        ]
        self.coastline_args = {}
        tmp_waffles = Waffle()
        tmp_coastline = WafflesCoastline()
        for kpam, kval in kwargs.items():
            if kpam not in tmp_waffles.__dict__:
                if kpam in tmp_coastline.__dict__:
                    self.coastline_args[kpam] = kval

        for kpam, kval in self.coastline_args.items():
            del kwargs[kpam]

        if exclude_lakes:
            self.coastline_args['want_lakes'] = True
            self.coastline_args['invert_lakes'] = True

        if mode is not None:
            pre_mode = mode
            
        self.pre_mode = pre_mode
        self.pre_mode_args = {}
        if self.pre_mode not in self.valid_modes:
            utils.echo_warning_msg(
                (f'{self.pre_mode} is not a valid waffles module! '
                 'falling back to `gmt-surface`')
            )
            self.pre_mode = 'gmt-surface'

        tmp_waffles_mode = WaffleFactory(mod=self.pre_mode)._acquire_module()
        for kpam, kval in kwargs.items():
            if kpam not in tmp_waffles.__dict__:
                if kpam in tmp_waffles_mode.__dict__:
                    self.pre_mode_args[kpam] = kval

        for kpam, kval in self.pre_mode_args.items():
            del kwargs[kpam]

        self.final_mode = final_mode
        self.final_mode_args = {}
        if self.final_mode not in self.valid_modes:
            self.final_mode = 'IDW'
            
        super().__init__(**kwargs)
        self.pre_count = utils.int_or(pre_count, 1)
        self.min_weight = min_weight        
        if min_weight is not None:
            weight_levels = min_weight
        
        if weight_levels is not None:
            self.weight_levels = [utils.float_or(w) for w in weight_levels.split('/')]
            self.weight_levels = self.weight_levels[:self.pre_count]
        else:
            self.weight_levels = []
            
        if inc_levels is not None:
            self.inc_levels = [utils.str2inc(i) for i in inc_levels.split('/')]
            self.inc_levels = self.inc_levels[:self.pre_count]
        else:
            self.inc_levels = []
            
        self.landmask = landmask
        self.invert_landmask = invert_landmask
        self.pre_upper_limit = utils.float_or(pre_upper_limit, -0.1) \
            if landmask \
               else None
        
        self.filter_outliers = filter_outliers
        self.want_supercede = want_supercede
        self.pre_smoothing = utils.float_or(pre_smoothing)
        if self.pre_smoothing is not None:
            self.pre_smoothing = ['blur:blur_factor={}'.format(self.pre_smoothing)]
            
        self.flatten = utils.float_or(flatten)
        self.want_weight = True
        self.pre_verbose = pre_verbose
        self.keep_pre_surfaces = keep_pre_surfaces

        
    ## todo: remove coastline after processing...
    def generate_coastline(self, pre_data=None):
        cst_region = self.p_region.copy()
        cst_region.wmin = self.weight_levels[0]
        utils.echo_msg('coast region is: {}'.format(cst_region))
        cst_fn = '{}_cst'.format(
            os.path.join(self.cache_dir, os.path.basename(self.name))
        )
        this_coastline = 'coastline:{}'.format(
            factory.dict2args(self.coastline_args)
        )
        coastline = WaffleFactory(
            mod=this_coastline,
            data=pre_data,
            src_region=cst_region,
            want_weight=True,
            min_weight=self.weight_levels[0],
            xinc=self.xinc,
            yinc=self.yinc,
            name=cst_fn,
            node=self.node,
            dst_srs=self.dst_srs,
            srs_transform=self.srs_transform,
            clobber=True,
            cache_dir=self.cache_dir,
            verbose=True)._acquire_module()
        coastline.initialize()
        coastline.generate()

        if coastline is not None:
            return('{}.shp:invert={}'.format(coastline.name, self.invert_landmask))
        else:
            return(None)

        
    def run(self):
        ## set the weights if not already set correctly
        self.weight_levels.sort(reverse=True)
        if len(self.weight_levels) < self.pre_count+1:
            tmp_pre = self.pre_count - len(self.weight_levels)
            while tmp_pre >= 0:
                if len(self.weight_levels) == 0:
                    # self.stack isn't set yet!
                    tmp_weight = gdalfun.gdal_percentile(
                        self.stack, perc=99, band=3
                    ) 
                    tmp_weight = utils.float_or(tmp_weight, 1)
                else:
                    tmp_weight = self.weight_levels[-1]/(tmp_pre + 1)
                    if tmp_weight == 0:
                        tmp_weight = 1-e20

                if tmp_pre == 0:
                    tmp_weight = 0
                    
                self.weight_levels.append(tmp_weight)                
                tmp_pre -= 1

        ## set the increments if not already set correctly
        self.inc_levels.sort()
        if len(self.inc_levels) < self.pre_count+1:
            tmp_pre = self.pre_count - len(self.inc_levels)
            if len(self.inc_levels) == 0:
                tmp_xinc = self.xinc
                self.inc_levels.append(tmp_xinc)
                
            for t in range(1, tmp_pre+1):
                tmp_xinc = float(self.inc_levels[-1] * 3)
                tmp_yinc = float(self.inc_levels[-1] * 3)
                self.inc_levels.append(tmp_xinc)

        if self.inc_levels[0] != self.xinc or len(self.inc_levels) < self.pre_count+1:
            self.inc_levels.insert(0, self.xinc)
            self.inc_levels = self.inc_levels[:self.pre_count+1]
            self.inc_levels.sort()
        
        if self.verbose:
            utils.echo_msg_bold('==============================================')
            utils.echo_msg('')
            #utils.echo_msg_bold('cudem pre-mode: {}'.format(self.pre_mode))
                
            utils.echo_msg_bold(
                f'cudem generating {self.pre_count} pre-surface(s):'
            )
            for i in range(len(self.inc_levels)-1, -1, -1):
                utils.echo_msg_bold(
                    '* {} {} <{}> @ {} using data with a minimum weight of {}'.format(
                        'pre_surface' if i !=0 else 'final surface',
                        i,
                        self.pre_mode if i == self.pre_count \
                        else 'stacks' if i != 0 \
                        else self.final_mode,
                        self.inc_levels[i],
                        self.weight_levels[i],
                    )
                )
                
            utils.echo_msg('')
            utils.echo_msg_bold(f'cudem using stack file: {self.stack}')
            if self.landmask:
                if isinstance(self.landmask, str):
                    if os.path.exists(self.landmask.split(':')[0]):
                        utils.echo_msg_bold(
                            f'cudem using coastline: {self.landmask}'
                        )
                    else:
                        self.landmask = True
                        
                if isinstance(self.landmask, bool):
                    utils.echo_msg_bold(
                        'cudem using coastline: Waffles COASTLINE module'
                    )
                
            utils.echo_msg_bold(f'cudem flattening: {self.flatten}')
            utils.echo_msg_bold(f'cudem output DEM: {self.name}')
            utils.echo_msg('')
            utils.echo_msg_bold('==============================================')

        orig_stack = self.stack
        pre = self.pre_count
        pre_weight = 0 # initial run will use all data weights
        pre_region = self.p_region.copy()
        pre_region.wmin = None

        # _pre_name_minus = os.path.join(
        #     self.cache_dir, utils.append_fn('_pre_surface', pre_region, pre)
        # )
        # ## initial data to pass through surface (stack)
        stack_data_entry = (f'{self.stack},200:band_no=1:weight_mask=3:'
                            'uncertainty_mask=4:sample=average,1')
        pre_data = [stack_data_entry]
         
        ## generate coastline
        pre_clip = None
        if self.landmask:            
            if isinstance(self.landmask, str):
                if os.path.exists(self.landmask.split(':')[0]):
                    # todo: update to make 'invert' an option
                    pre_clip = '{}:invert={}'.format(self.landmask, self.invert_landmask) 

            if pre_clip is None:
                coast_data = [
                    (f'{self.stack},200:band_no=1:weight_mask=3:'
                     'uncertainty_mask=4:sample=cubicspline,1')]
                coastline = self.generate_coastline(pre_data=coast_data)
                #utils.echo_msg_bold(coastline)
                pre_clip = coastline

        ## Grid/Stack the data `pre` times concluding in full
        ## resolution with data > min_weight
        pre_surfaces = []
        with tqdm(
                total=self.pre_count+1,
                desc='generating CUDEM surfaces',
                leave=self.verbose
        ) as pbar:
            for pre in range(self.pre_count, -1, -1): # pre_count - 1
                pre_xinc = self.inc_levels[pre]
                pre_yinc = self.inc_levels[pre]
                xsample = self.inc_levels[pre-1] if pre != 0 else self.xinc
                ysample = self.inc_levels[pre-1] if pre != 0 else self.yinc

                ## initial data to pass through surface (stack)

                
                ## if not final or initial output, setup the configuration
                ## for the pre-surface
                if pre != self.pre_count:
                    pre_weight = self.weight_levels[pre]
                    _pre_name_plus = os.path.join(
                        self.cache_dir, utils.append_fn('_pre_surface', pre_region, pre+1)
                    )
                    _pre_unc_name = f'{_pre_name_plus}_u.tif' \
                        if self.want_uncertainty \
                           else None
                    pre_data_entry = (f'{_pre_name_plus}.tif,200'
                                      f':uncertainty_mask={_pre_unc_name}'
                                      f':sample=cubicspline:check_path=True'
                                      f',{pre_weight-.1}')

                    ## rq filter points to previous pre_surface
                    # _pre_name_minus = os.path.join(
                    #     self.cache_dir, utils.append_fn('_pre_surface', pre_region, pre+1)
                    # )
                    # #rq_threshold = max(1, 5+(pre*2))
                    # rq_threshold = 5
                    # stack_data_entry = (f'{self.stack},200:band_no=1:weight_mask=3:'
                    #                     'uncertainty_mask=4:sample=average'
                    #                     f':pnt_fltrs="rq:threshold={rq_threshold}:raster={_pre_name_minus}.tif",1')
                    # #pre_data = [stack_data_entry]
                    
                    pre_data = [stack_data_entry, pre_data_entry]
                    pre_region.wmin = None#pre_weight

                ## reset pre_region for final grid
                if pre == 0:
                    pre_region = self.p_region.copy()
                    pre_region.wmin = None#self.weight_levels[pre]

                _pre_name = os.path.join(
                    self.cache_dir, utils.append_fn('_pre_surface', pre_region, pre)
                )
                # if self.verbose:
                #     utils.echo_msg('pre region: {}'.format(pre_region))
                #     utils.echo_msg('pre data: {}'.format(pre_data))
                #     utils.echo_msg('pre weight: {}'.format(pre_weight))

                last_fltr = [
                    (f'weights:stacks=True:weight_threshold={pre_weight}'
                     ':buffer_cells=1:verbose=False')
                     ]
                waffles_mod = '{}:{}'.format(
                    self.pre_mode,
                    factory.dict2args(self.pre_mode_args)
                ) if pre==self.pre_count \
                else 'stacks'\
                     if pre != 0 \
                        else 'IDW'
                utils.echo_msg(
                    'cudem gridding surface {} @ {} {}/{} to {} using {}...'.format(
                        pre, pre_region, pre_xinc, pre_yinc, _pre_name, waffles_mod
                    )
                )
                pre_surface = WaffleFactory(
                    mod=waffles_mod,
                    data=pre_data,
                    src_region=pre_region,
                    xinc=pre_xinc,
                    yinc=pre_yinc,
                    xsample=xsample,
                    ysample=ysample,
                    name=_pre_name,
                    node='pixel',
                    want_weight=True,
                    want_uncertainty=self.want_uncertainty,
                    dst_srs=self.dst_srs,
                    srs_transform=self.srs_transform,
                    clobber=True,
                    verbose=self.pre_verbose,
                    clip=pre_clip if pre !=0 else None,
                    stack_mode='mixed:weight_threshold={}'.format(pre_weight),
                    #self.weight_levels[0]),# if pre == 0 else 'mean',
                    #stack_mode='supercede' if (pre == 0 and self.want_supercede)
                    #else self.stack_mode,
                    upper_limit=self.pre_upper_limit if pre != 0 else None,
                    keep_auxiliary=False,
                    fltr=self.pre_smoothing if pre != 0 else None,
                    #fltr=last_fltr,
                    percentile_limit=self.flatten if pre == 0 else None,
                    cache_dir=self.cache_dir,
                )._acquire_module()
                pre_surface.initialize()
                pre_surface.generate()
                pbar.update()
                pre_surfaces.append(pre_surface.name)

        ## todo add option to flatten here...or move flatten up
        gdalfun.cudem_flatten_no_data_zones(
            pre_surface.fn, dst_dem=self.fn, band=1, size_threshold=1
        )

        ## reset the stack for uncertainty
        self.stack = orig_stack
        if not self.keep_pre_surfaces:
            utils.remove_glob(*[f"{x}*" for x in pre_surfaces])
        
        return(self)

    
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
            if os.path.exists(os.path.join(utils.cudem_data, self.prox_errs)):
                self.prox_errs = os.path.join(utils.cudem_data, self.prox_errs)
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
            out_prox = utils.make_temp_fn(
                '{}_prox.tif'.format(self.params['mod_args']['waffles_module'])
            )

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
            out_slope = utils.make_temp_fn(
                '{}_slope.tif'.format(
                    self.params['mod_args']['waffles_module']
                )
            )

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
        xyz_ds = dlim.XYZFile(
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
        with tqdm(
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
            
        last_ec_d = None
        s_dp = []
        status = 0
        trains = self._regions_sort(trainers)
        all_trains = [x for s in trains for x in s[:5]]
        tot_trains = len(all_trains)
        
        with tqdm(
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
                'running UNCERTAINTY module using {}...'.format(
                    self.params['mod_args']['waffles_module']
                )
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

    
class WafflesCUBE(Waffle):
    """
    BathyCUBE - doesn't seem to work as expected, 
    likely doing something wrong here...

    https://github.com/noaa-ocs-hydrography/bathygrid
    
    """
    
    def __init__(
            self,
            chunk_size=None,
            chunk_buffer=40,
            **kwargs):
        """generate a `CUBE` dem"""
        
        super().__init__(**kwargs)
        self.chunk_size = utils.int_or(chunk_size)
        self.chunk_step = None
        self.chunk_buffer = chunk_buffer

        
    def run(self):

        try:
            import bathycube.cube as cube
        except:
            utils.echo_error_msg(
                'could not import bathycube, it may not be installed'
            )
            return(self)

        from scipy import optimize
        
        # _x = []
        # _y = []
        # _z = []
        # for this_xyz in self.yield_xyz():
        #     _x.append(this_xyz.x)
        #     _y.append(this_xyz.y)
        #     _z.append(this_xyz.z)

        # _x = np.array(_x)
        # _y = np.array(_y)
        # _z = np.array(_z)

        # print(_z.size)
        # print(_z.shape)
        # print(len(_z))
        
        if self.chunk_size is None:
            n_chunk = int(self.ds_config['nx'] * .1)
            n_chunk = 10 if n_chunk < 10 else n_chunk
        else: n_chunk = self.chunk_size
        if self.chunk_step is None:
            n_step = int(n_chunk/2)
        else: n_step = self.chunk_step

        points_ds = gdal.Open(self.stacked_rasters['z'])
        points_band = points_ds.GetRasterBand(1)
        points_no_data = points_band.GetNoDataValue()
        
        #try:
        interp_ds = points_ds.GetDriver().Create(
            self.fn,
            points_ds.RasterXSize,
            points_ds.RasterYSize,
            bands=1,
            eType=points_band.DataType,
            options=["BLOCKXSIZE=256",
                     "BLOCKYSIZE=256",
                     "TILED=YES",
                     "COMPRESS=LZW",
                     "BIGTIFF=YES"]
        )
        interp_ds.SetProjection(points_ds.GetProjection())
        interp_ds.SetGeoTransform(points_ds.GetGeoTransform())
        #interp_ds.SetGeoTransform(ds_config['geoT'])
        interp_band = interp_ds.GetRasterBand(1)
        interp_band.SetNoDataValue(np.nan)

        uncert_ds = points_ds.GetDriver().Create(
            '{}_unc.tif'.format(self.name),
            points_ds.RasterXSize,
            points_ds.RasterYSize,
            bands=1,
            eType=points_band.DataType,
            options=["BLOCKXSIZE=256",
                     "BLOCKYSIZE=256",
                     "TILED=YES",
                     "COMPRESS=LZW",
                     "BIGTIFF=YES"]
        )
        uncert_ds.SetProjection(points_ds.GetProjection())
        uncert_ds.SetGeoTransform(points_ds.GetGeoTransform())
        #uncert_ds.SetGeoTransform(ds_config['geoT'])
        uncert_band = uncert_ds.GetRasterBand(1)
        uncert_band.SetNoDataValue(np.nan)
            
        # except:
        #     return(self)
        
        if self.verbose:
            utils.echo_msg(
                'buffering srcwin by {} pixels'.format(self.chunk_buffer)
            )

        _x, _y = np.mgrid[0:ds_config['nx'], 0:ds_config['ny']]
        _x = _x.ravel()
        _y = _y.ravel()

        _z = points_band.ReadAsArray()
        _z = _z.T
        _z = _z.ravel()
        point_indices = np.nonzero(_z != ds_config['ndv'])
        point_values = _z[point_indices]        
        xi = _x[point_indices]
        yi = _y[point_indices]

        thu = np.ones(len(point_values))
        thu[:] = 2
        #thu = thu[point_indices]
        a = .25
        b = .0075
        tvu = np.sqrt((a**2 + (b * abs(point_values))**2))
        numrows, numcols = (ds_config['nx'], ds_config['ny'])
        res_x, res_y = ds_config['geoT'][1], ds_config['geoT'][5]*-1
        depth_grid, uncertainty_grid, ratio_grid, numhyp_grid = cube.run_cube_gridding(
            point_values, thu, tvu, xi, yi, self.ds_config['nx'], self.ds_config['ny'],
            min(_x), max(_y), 'local', 'order1a', 1, 1)
        print(depth_grid)
        depth_grid = np.flip(depth_grid)
        depth_grid = np.fliplr(depth_grid)
        interp_band.WriteArray(depth_grid)
        uncertainty_grid = np.flip(uncertainty_grid)
        uncertainty_grid = np.fliplr(uncertainty_grid)
        uncert_band.WriteArray(uncertainty_grid)
            
        return(self)

    
class WafflesBGrid(Waffle):
    """
    BathyCUBE - doesn't seem to work as expected, 
    likely doing something wrong here...

    https://github.com/noaa-ocs-hydrography/bathygrid
    
    """
    
    def __init__(
            self,
            chunk_size=None,
            chunk_buffer=40,
            **kwargs):
        """generate a `CUBE` dem"""
        
        super().__init__(**kwargs)
        self.chunk_size = utils.int_or(chunk_size)
        self.chunk_step = None
        self.chunk_buffer = chunk_buffer

        
    def run(self):

        try:
            import bathycube.cube as cube
        except:
            utils.echo_error_msg(
                'could not import bathycube, it may not be installed'
            )
            return(self)

        from scipy import optimize
        
        # self._stacks_array(
        #     out_name='{}_cube_stack'.format(self.name),
        #     supercede=self.supercede
        # )
        # n = '{}_cube_stack_s.tif'.format(self.name)
        # w = '{}_cube_stack_w.tif'.format(self.name)
        # c = '{}_cube_stack_c.tif'.format(self.name)
        _x = []
        _y = []
        _z = []
        for this_xyz in self.yield_xyz():
            _x.append(this_xyz.x)
            _y.append(this_xyz.y)
            _z.append(this_xyz.z)

        _x = np.array(_x)
        _y = np.array(_y)
        _z = np.array(_z)

        dtyp = [('x', np.float64),
                ('y', np.float64),
                ('z', np.float32),
                ('tvu', np.float32),
                ('thu', np.float32)]
        data = np.empty(len(_x), dtype=dtyp)
        data['x'] = _x
        data['y'] = _y
        data['z'] = _z
        
        # print(_z.size)
        # print(_z.shape)
        # print(len(_z))
        
        if self.chunk_size is None:
            n_chunk = int(ds_config['nx'] * .1)
            n_chunk = 10 if n_chunk < 10 else n_chunk
        else: n_chunk = self.chunk_size
        if self.chunk_step is None:
            n_step = int(n_chunk/2)
        else: n_step = self.chunk_step

        if self.verbose:
            utils.echo_msg(
                'buffering srcwin by {} pixels'.format(self.chunk_buffer)
            )

        thu = np.ones(len(_z))
        thu[:] = 2
        #thu = thu[point_indices]
        a = .25
        b = .0075
        tvu = np.sqrt((a**2 + (b * abs(_z))**2))

        data['tvu'] = tvu
        data['thu'] = thu
        numrows, numcols = (ds_config['nx'], ds_config['ny'])
        res_x, res_y = ds_config['geoT'][1], ds_config['geoT'][5]*-1


        from bathygrid.convenience import create_grid
        # a single resolution grid that is entirely within computer memory
        bg = create_grid(folder_path='./', grid_type='single_resolution')
        # add points from two multibeam lines, EPSG:26917 with vertical reference 'waterline'
        bg.add_points(data, 'test1', ['line1', 'line2'], 26965, 'waterline')
        print(bg.points_count)
        assert not bg.is_empty

        # grid by looking up the mean depth of each tile to determine resolution
        bg.grid()

        print(bg.resolutions)
        out_tif = os.path.join(bg.output_folder, self.fn)
        bg.export(out_tif, export_format='geotiff')
        
        # # check to see if the new bags are written
        # new_bag = os.path.join(bg.output_folder, 'outtiff_0.5.bag')
        # new_bag_two = os.path.join(bg.output_folder, 'outtiff_1.0.bag')
        # assert os.path.exists(new_bag)
        # assert os.path.exists(new_bag_two)

        # # Get the total number of cells in the variable resolution grid for each resolution
        # bg.cell_count
        # bg.coverage_area
        
        # depth_grid, uncertainty_grid, ratio_grid, numhyp_grid = cube.run_cube_gridding(
        #     point_values,
        #     thu,
        #     tvu,
        #     xi,
        #     yi,
        #     ds_config['nx'],
        #     ds_config['ny'],
        #     min(_x),
        #     max(_y),
        #     'local',
        #     'order1a',
        #     1,
        #     1,
        # )
        # print(depth_grid)
        # depth_grid = np.flip(depth_grid)
        # depth_grid = np.fliplr(depth_grid)
        # interp_band.WriteArray(depth_grid)
        # uncertainty_grid = np.flip(uncertainty_grid)
        # uncertainty_grid = np.fliplr(uncertainty_grid)
        # uncert_band.WriteArray(uncertainty_grid)
            
        return(self)

    
class WafflesPatch(Waffle):
    """PATCH an existing DEM with new data.
    
    Patch an existing DEM with data from the datalist.

    -----------
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
        super().__init__(**kwargs)
        self.radius = utils.str2inc(radius)
        self.min_weight = utils.float_or(min_weight)
        self.max_diff = utils.float_or(max_diff)
        self.dem = dem

        
    def yield_diff(self, src_dem, max_diff=.25):
        '''query a gdal-compatible grid file with xyz data.
        out_form dictates return values

        yields out_form results'''

        def con_dec(x, dec):
            '''Return a float string with n decimals
            (used for ascii output).'''

            if x is None:
                utils.echo_error_msg(
                    'Attempting to convert a None value.'
                )
                return
            fstr = "%." + str(dec) + "f"
            return fstr % x
      
        try:
            ds = gdal.Open(src_dem)
        except: ds = None
        if ds is not None:
            ds_config = gdalfun.gdal_infos(ds)
            ds_band = ds.GetRasterBand(1)
            ds_gt = ds_config['geoT']
            ds_nd = ds_config['ndv']
            tgrid = ds_band.ReadAsArray()
            dsband = ds = None

            for xyz in self.yield_xyz():
                #if xyz.x > ds_gt[0] and xyz.y < float(ds_gt[3]):
                try: 
                    xpos, ypos = utils._geo2pixel(
                        xyz.x, xyz.y, ds_gt, 'pixel'
                    )
                    g = tgrid[ypos, xpos]
                except: g = ds_nd
                
                if g != ds_nd:
                    d = xyz.z - g
                    s = math.fabs(d / (xyz.z + (g+0.00000001)))

                    if s < max_diff:
                        xyz.z = d
                        yield(xyz)
            ds = None

            
    def query_dump(
            self,
            dst_port=sys.stdout,
            encode=False,
            max_diff=.25,
            **kwargs
    ):
        for xyz in self.yield_diff(self.dem, max_diff):
            xyz.dump(
                include_w = self.want_weight,
                include_u = self.want_uncertainty,
                dst_port=dst_port,
                encode=encode,
                **kwargs
            )

            
    def run(self):
        if dem is not None:
            if os.path.exists(dem):
                self.dem = dem
        elif os.path.exists('{}.tif'.format(self.name)):
            self.dem = '{}.tif'.format(self.name)
            self.name = '{}_update'.format(self.name)
        else:
            utils.echo_error_msg(
                'must specify DEM to patch (:dem=fn) to run the patch module.'
            )
            return(None)

        dem_infos = gdalfun.gdal_infos(self.dem)
        dem_region = regions.Region().from_geo_transform(
            geo_transform=dem_infos['geoT'],
            x_count=dem_infos['nx'],
            y_count=dem_infos['ny']
        )

        if not regions.regions_intersect_p(self.region, dem_region):
            utils.echo_error_msg(
                'input region does not intersect with input DEM'
            )

        ## grid the difference to array using query_dump / num
        ## polygonize the differences and add small buffer (1% or so)
        ## make zero array, inverse clipped to buffered polygonized diffs
        ## surface zero array and diffs...
        ## add surfaced diffs to self.dem
        diff_cmd = 'gmt blockmedian {region} -I{xinc}/{yinc} | gmt surface {region} -I{xinc}/{yinc} -G_diff.tif=gd+n{ndv}:GTiff -T.1 -Z1.2 -V -Lud -Lld'.format(
            region=self.region.format('gmt'),
            xinc=self.xinc,
            yinc=self.yinc,
            ndv=self.ndv,
            radius=self.radius
        )

        out, status = utils.run_cmd(
            diff_cmd,
            verbose=self.verbose,
            data_fun=lambda p: self.query_dump(
                dst_port=p, encode=True, max_diff=self.max_diff
            )
        )

        utils.echo_msg('smoothing diff grid...')
        smooth_dem, status = gdalfun.gdal_blur('_diff.tif', '_tmp_smooth.tif', 5)

        if self.xinc != dem_infos['geoT']:
            utils.echo_msg('resampling diff grid...')
            if gdalfun.sample_warp(
                    '_diff.tif',
                    '__tmp_sample.tif',
                    dem_infos['geoT'][1],
                    -1*dem_infos['geoT'][5],
                    src_region=dem_region,
                    verbose=self.verbose
            )[1] == 0:
                os.replace('__tmp_sample.tif', '_tmp_smooth.tif')
            else:
                utils.echo_warning_msg('failed to resample diff grid')

        #if self.xinc != dem_infos['geoT']:
        if self.verbose:
            utils.echo_msg('resampling diff grid...')
            
        diff_ds = gdalfun.sample_warp(
            '_diff.tif',
            None,
            dem_infos['geoT'][1],
            -1*dem_infos['geoT'][5],
            src_region=dem_region,
            verbose=self.verbose
        )[0]
        
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
        gdalfun.gdal_write(
            update_dem_arr, '{}.tif'.format(self.name), dem_infos
        )
        
        diff_ds = dem_ds = None
        return(self)

# class WafflesMMI(Waffle):
#     """ Multiresolution interpolation function
#     [mapa, xmin,xmax, ymin,ymax, mapa0]=mrinterp(x,y,z, dxy, slev, xmin,xmax,ymin,ymax)
#     x, y, z: coordinates and value
#     dxy: spatial resolution of the raster (same for X and Y)
#     slev: first interpolation level (optional) and smoothing level (number of levels without exact interpolation)
#     mapa: interpolated raster matrix
#     xmin, xmax, ymin, ymax: extent
#     mapa0: raster with the projection of the original data points
#     """

#     def __init__(self, **kwargs):
#         super().__init__(**kwargs)
        
#     #utility functions
#     #sum2
#     def sum2(self, m):
#         return(np.sum(m))
    
#     #duplicate
#     def dupl2(self, a):
#         ha, wa = np.shape(a)
#         #ha = size(a,1);
#         #wa = size(a,2);
#         b = np.nan(2*a.size)
#         #b = nan(2*size(a),class(a))
#         b(1:2:2*ha,1:2:2*wa) = a
#         b(2:2:2*ha,1:2:2*wa) = a
#         b(1:2:2*ha,2:2:2*wa) = a
#         b(2:2:2*ha,2:2:2*wa) = a


#         #function [mapa, xmin,xmax, ymin,ymax, mapa0]=mrinterp(x,y,z, dxy, slev, xmin,xmax,ymin,ymax)

#     def next_power_of_2(self, x):  
#         return 1 if x == 0 else 2**(x - 1).bit_length()
        
#     def run(self):
        
#         #data type
#         dcls=class(z);

#         #compute map size
#         #if( ~exist('xmin') || ~exist('xmax') )
#         #xmin=min(x);
#         #xmax=max(x);

#         #dx=xmax-xmin;
#         dx = self.region.xmax - self.region.xmin
#         #wa=ceil(dx/dxy);
#         wa = np.ceil(dx / dxy)
#         #Nlx=nextpow2(wa);
#         Nlx = self.next_power_of_2(wa)
#         #w=pow2(Nlx);
#         w = Nlx ** 2
#         #xmin=xmin-(w-wa)*dxy/2;
#         xmin = xmin - (w - wa) * dxy / 2
#         #xmax=xmin+w*dxy;
#         xmax = xmin + w * dxy

#         #if( ~exist('ymin') || ~exist('ymax') )
#         #ymin=min(y);
#         #ymax=max(y);
#         #endif
#         #dy=ymax-ymin;
#         dy = self.region.ymax - self.region.ymin
#         #ha=ceil(dy/dxy);
#         ha = np.ceil(dy / dxy)
#         #Nly=nextpow2(ha);
#         Nly = self.next_power_of_2(ha)
#         #h=pow2(Nly);
#         h = Nly ** 2
#         #ymin=ymin-(h-ha)*dxy/2;
#         ymin = ymin - (h - ha) * dxy / 2
#         #ymax=ymin+h*dxy;
#         ymax = ymin + h * dxy

#         #m=(xmin<=x & x<=xmax & ymin<y & y<=ymax);        
#         #x=x(m);
#         #y=y(m);
#         #z=z(m);
#         #clear m
  
#   Nl=min([Nlx,Nly]);
#   bl=min([w,h]);

#   %raw projection
#   i=1+floor((ymax-y)/dxy);
#   j=1+floor((x-xmin)/dxy);
#   ij=i+h*(j-1);
  
#   %map number of points and sum of values
#   map0N=zeros(h,w, dcls);
#   map0S=zeros(h,w, dcls);
#   for n=1:length(z)
#     ijn=ij(n);
#     map0N(ijn)=map0N(ijn)+1;
#     map0S(ijn)=map0S(ijn)+z(n);
#   endfor

#   %OPTIONAL: prune scarcely populated pixels
#   if( false )
#     nmin=1;
#     while( sum(map0N(:)>nmin) > 0.99*sum(map0N(:)>0) )
#       nmin=nmin+1;
#     endwhile
#     map0S(map0N>=nmin)=map0S(map0N>=nmin)./map0N(map0N>=nmin);
#     map0S(map0N<nmin)=0;
#     map0N=(map0N>=nmin); %start afresh from equally weithted pixels (!)
#   endif
  
#   mapa0=nan(size(map0S), dcls);
#     mapa0(map0N>0) = map0S(map0N>0) ./ map0N(map0N>0);
 
#   clear x y z i j ij
 
#   %scaling levels
#   if( length(slev) == 1 )
#     slev0=0;
#     slev1=slev;
#   else
#     slev0=slev(1);
#     slev1=slev(2);
#   endif
    
#   %multiscale iteration
#   b=bl/2^slev0;
#   for l=0:Nl-slev0
#     ni=h/b;
#     nj=w/b;
    
#     %data sums and number
#     maplS=zeros(ni,nj, dcls);
#     maplN=zeros(ni,nj, dcls);
#     if( l <= Nl - slev0 - slev1 ) %smoothing level
#       if( b > 1 )
#         for i=1:ni
#           for j=1:nj
#             maplN(i,j)=sum2( map0N(1+(i-1)*b:i*b, 1+(j-1)*b:j*b) );
#             maplS(i,j)=sum2( map0S(1+(i-1)*b:i*b, 1+(j-1)*b:j*b) );
#           endfor
#         endfor
#       else
#         maplN=1.0*map0N;
#         maplS=map0S;
#       endif
#     endif
    
#     if( l>0 )
#       %mK=[0 1 0; 1 0 1; 0 1 0];
#       mK=[0.614 1 0.614; 1 0 1; 0.614 1 0.614];
#         mK=mK/sum(mK(:));

#       %compute low resolution mean height
#       ma=(mapaN>0); %where data exists at lower resolution
#       mapaz=zeros(size(mapaS), dcls);
#       mapaz(ma)=mapaS(ma)./mapaN(ma);
      
#       %pad low-resolution map with current resolution data
#       m=(maplN>0); %where data exists at current resolution
#       mapaN(m)=maplN(m);
#       %mapaS(m)=maplS(m);
#       mapaz(m)=maplS(m)./maplN(m);
      
#       %compute smoothed version of low resolution variables
#       cmapaM=conv2(ma,mK, 'same');
#       cmapaN=conv2(mapaN,mK, 'same');
#       %cmapaS=conv2(mapaS,mK, 'same');
#       cmapaZ=conv2(mapaz,mK, 'same');
      
#       %fill in missing data at pixels with no data at current resolution
#       m=~m & (cmapaM>0); %where data does not exist at current resolution but it does at lower 
#       maplN(m)=cmapaN(m) ./ cmapaM(m);
#       maplS(m)=(cmapaZ(m) .* cmapaN(m)) ./ power(cmapaM(m),2);

#       clear m ma
#       clear cmapaM cmapaN cmapaZ %cmapaS
#       clear mapaz mapaN mapaS
      
#     endif
    
#     %next resolution
#     if( l<Nl-slev0)
#       b=b/2;
#       %disp( [ 'PROGRESS: b=' num2str(b) ] );
#       mapaS=dupl2(maplS)/4;
#       mapaN=dupl2(maplN)/4;
#     endif

#   endfor

#   %return mapa
#   mapa=nan(size(maplS), dcls);
#     mapa(maplN>0) = maplS(maplN>0) ./ maplN(maplN>0);
  
#   %one last check
#   if( any(size(mapa0) ~= size(mapa)) )
#     warning "mapa0 and mapa have different sizes"
#   endif
  
# endfunction

    
def nodata_count_mask(src_dem, band=1):
    """Make a data mask of nodata zones.

    output array will have the number of cells per nodata zone
    """
    
    mn = None
    with gdalfun.gdal_datasource(src_dem) as src_ds:
        if src_ds is not None:
            src_band = src_ds.GetRasterBand(band)
            src_arr = src_band.ReadAsArray().astype(float)
            src_config = gdalfun.gdal_infos(src_ds)
            src_arr[src_arr == src_config['ndv']] = np.nan

            ## generate the mask array
            msk_arr = np.zeros((src_config['ny'], src_config['nx']))
            msk_arr[np.isnan(src_arr)] = 1

            ## group adjacent non-zero cells
            l, n = scipy.ndimage.label(msk_arr)

            ## get the total number of cells in each group
            mn = scipy.ndimage.sum_labels(msk_arr, labels=l, index=l)
    
    return(mn)


class WaffleDEM:
    """WaffleDEM which holds a gdal DEM to process
    WaffleDEM(fn='module_output.tif')
    """
    
    def __init__(
            self,
            fn: str = 'this_waffle.tif',
            ds_config: dict = {},
            cache_dir: str = waffles_cache,
            verbose: bool = True,
            waffle: Waffle = None,
            want_scan: bool = True,
            co: list = ["COMPRESS=DEFLATE", "TILED=YES"]
    ):
        self.fn = fn # the dem filename
        self.ds_config = ds_config # a dem config dictionary (see gdalfun.gdal_infos)
        self.cache_dir = cache_dir # cache dir for auxiliary data
        self.verbose = verbose # verbosity
        self.dem_region = None # the dem regions.Region()
        self.waffle = waffle # the waffles module that generated the DEM
        self.want_scan = want_scan # scan DEM for min/max values
        self.co = co

        
    def initialize(self):
        if os.path.exists(self.fn):
            if self.fn.endswith('csg'):
                dem_ds = h5.File(self.fn)
                if dem_ds is not None:
                    dst_gt = tuple(
                        [float(x) for x in dem_ds['crs'].attrs['GeoTransform'].split()]
                    )
                    self.ds_config = gdalfun.gdal_set_infos(
                        dem_ds['lon'].size, dem_ds['lat'].size,
                        dem_ds['lon'].size * dem_ds['lat'].size,
                        dst_gt, dem_ds['crs'].attrs['crs_wkt'],
                        'h5', np.nan, 'h5', None, None
                    )
                    self.dem_region = regions.Region().from_geo_transform(
                        dst_gt, self.ds_config['nx'], self.ds_config['ny']
                    )
                    self.ds_config['zr'] \
                        = [dem_ds['stack/z'][...].min(), dem_ds['stack/z'][...].max()]
                    
                    dem_ds = None
                else:
                    utils.echo_warning_msg(f'could not open dem: {self.fn}')
            else:            
                dem_ds = gdal.Open(self.fn, 1)
                if dem_ds is not None:
                    self.ds_config = gdalfun.gdal_infos(dem_ds, scan=self.want_scan)
                    self.dem_region = regions.Region().from_geo_transform(
                        self.ds_config['geoT'], self.ds_config['nx'], self.ds_config['ny']
                    )

                    dem_ds = None
                else:
                    utils.echo_warning_msg(f'could not open dem: {self.fn}')
            
        return(self)

    
    def valid_p(self):
        """check if the WAFFLES DEM appears to be valid"""

        if self.fn is None or not os.path.exists(self.fn):
            utils.echo_warning_msg(
                f'{self.fn} does not exist'
            )
            return(False)
        
        self.initialize()
        if self.ds_config is None:
            utils.echo_warning_msg(
                f'could not parse dem: {self.fn}'
            )
            return(False)
        
        if not 'zr' in self.ds_config:
            utils.echo_msg(
                f'dem {self.fn} has no z values?'
            )
            return(False)

        if self.ds_config['raster_count'] == 0:
            utils.echo_warning_msg(
                f'dem {self.fn} has no bands'
            )
            return(False)
        else:
            band_check = []
            for band_num in range(1, self.ds_config['raster_count']+1):
                #band_infos = gdalfun.gdal_infos(self.fn, scan=True, band=band_num)
                if np.isnan(self.ds_config['zr'][0]) \
                   or np.isnan(self.ds_config['zr'][1]):
                    band_check.append(0)
                else:
                    band_check.append(1)

            if not np.any(band_check):
                utils.echo_warning_msg(
                    f'dem {self.fn} is all nan'
                )
                return(False)
                
        return(True)

    
    def process(self, filter_ = None, ndv = None, xsample = None, ysample = None,
                region = None, node='pixel', clip_str = None, upper_limit = None,
                lower_limit = None, size_limit = None, proximity_limit = None,
                percentile_limit = None, expand_limit = None, dst_srs = None, dst_fmt = None,
                dst_fn = None, dst_dir = None, set_metadata = True, stack_fn = None, mask_fn = None,
                unc_fn = None, flatten_nodata_values = False, want_nc = False,
                want_h5 = False):
        """Process the DEM using various optional functions.

        set the nodata value, srs, metadata, limits; resample, filter, clip, cut, output.
        """

        if self.verbose:
            utils.echo_msg(
                f'post processing DEM {self.fn}...'
            )
            
        if self.ds_config is None:
            self.initialize()
            
        ## set nodata value
        if ndv is not None:
            self.set_nodata(ndv)
        else:
            if self.ds_config['ndv'] is not None:
                self.set_nodata(self.ds_config['ndv'])

        ## set interpolation limits
        self.set_interpolation_limits(
            stack_fn=stack_fn,
            size_limit=size_limit,
            proximity_limit=proximity_limit,
            percentile_limit=percentile_limit,
            expand_limit=expand_limit
        )

        ## filtering the DEM will change the weights/uncertainty
        if filter_ is not None:
            self.filter_(filter_)

        ## resamples all bands
        self.resample(xsample=xsample, ysample=ysample, ndv=ndv, region=region)

        ## setting limits will change the weights/uncertainty for flattened data
        self.set_limits(upper_limit=upper_limit, lower_limit=lower_limit)

        ## flatten all nodata values
        if flatten_nodata_values:
            flattened_fn = utils.make_temp_fn(
                '{}_flat.tif'.format(utils.fn_basename2(self.fn)),
                self.cache_dir
            )
            gdalfun.cudem_flatten_no_data_zones(
                self.fn, dst_dem=flattened_fn, band=1, size_threshold=1
            )
            os.rename(flattened_fn, self.fn)

        ## clip/cut
        self.clip(clip_str=clip_str)
        if region is not None:
            self.cut(region=region, node='grid')#'pixel' if node == 'grid' else 'grid')

        ## set projection, metadata, reformat and move to final location
        if dst_srs is not None:
            self.set_srs(dst_srs=dst_srs)

        if set_metadata:
            self.set_metadata(node=node)

        ## put everything in a hdf5
        if want_h5:
            if stack_fn is not None:
                self.stack2hd5(
                    stack_fn=stack_fn, mask_fn=mask_fn, unc_fn=unc_fn, node=node
                )

        ## put everything in a netcdf
        if want_nc:
            if stack_fn is not None:
                self.stack2nc(
                    stack_fn=stack_fn, mask_fn=mask_fn, unc_fn=unc_fn, node=node
                )
            
        ## reformat and move final output
        self.reformat(out_fmt=dst_fmt)
        self.move(out_fn=dst_fn, out_dir=dst_dir)

        if self.verbose:
            utils.echo_msg(
                f'Processed DEM: {self.fn}'
            )

            
    def set_nodata(self, ndv):
        if self.ds_config['ndv'] != ndv:
            gdalfun.gdal_set_ndv(
                self.fn, ndv=ndv, convert_array=True, verbose=self.verbose
            )
            self.ds_config['ndv'] = ndv

            if self.verbose:
                utils.echo_msg(
                    f'set nodata value to {ndv}.'
                )

                
    def filter_(self, fltr=[]):
        for f in fltr:
            if f.split(':')[0] in grits.GritsFactory._modules.keys():
                grits_filter = grits.GritsFactory(
                    mod=f,
                    src_dem=self.fn
                )._acquire_module()
                if grits_filter is not None:
                    if 'stacks' in grits_filter.kwargs.keys():
                        if grits_filter.kwargs['stacks']:
                            continue

                    grits_filter()
                    os.replace(grits_filter.dst_dem, self.fn)

                    
    def resample(
            self,
            region=None,
            xsample=None,
            ysample=None,
            ndv=-9999,
            sample_alg='cubicspline'
    ):
        if xsample is not None or ysample is not None:
            #warp_fn = os.path.join(self.cache_dir, '__tmp_sample.tif')
            warp_fn = utils.make_temp_fn('__tmp_sample.tif', temp_dir = self.cache_dir)
            if gdalfun.sample_warp(
                    self.fn,
                    warp_fn,
                    xsample,
                    ysample,
                    src_region=region,
                    sample_alg=sample_alg,
                    ndv=ndv,
                    verbose=self.verbose,
                    co=self.co
            )[1] == 0:
                os.replace(warp_fn, self.fn)
                self.initialize()

            if self.verbose:
                utils.echo_msg(
                    'resampled data to {}/{} using {}.'.format(
                        xsample, ysample, sample_alg
                    )
                )

                
    def clip(self, clip_str=None):
        ## todo: update for multi-band
        ## todo: re-add coastline option
        if clip_str is not None:
            clip_args = {}
            cp = clip_str.split(':')
            clip_args['src_ply'] = cp[0]
            clip_args['verbose'] = self.verbose
            clip_args = factory.args2dict(cp[1:], clip_args)

            if clip_args['src_ply'] == 'coastline':
                self.coast = WaffleFactory(
                    mod='coastline:polygonize=False',
                    data=self.data_,
                    src_region=self.p_region,
                    xinc=self.xsample if self.xsample is not None else self.xinc,
                    yinc=self.ysample if self.ysample is not None else self.yinc,
                    name='tmp_coast',
                    node=self.node,
                    want_weight=self.want_weight,
                    want_uncertainty=self.want_uncertainty,
                    dst_srs=self.dst_srs,
                    srs_transform=self.srs_transform,
                    clobber=True,
                    cache_dir=self.cache_dir,
                    verbose=self.verbose
                )._acquire_module()
                self.coast.initialize()
                self.coast.generate()
                gdalfun.gdal_mask(
                    fn,
                    self.coast.fn, '__tmp_clip__.tif',
                    msk_value=1,
                    verbose=self.verbose
                )
                os.replace('__tmp_clip__.tif', f'{fn}')

            if os.path.exists(clip_args['src_ply']):
                if ogr.Open(clip_args['src_ply']) is not None:
                    tmp_clip = utils.make_temp_fn(
                        '__tmp_clip__.tif', temp_dir=self.cache_dir
                    )
                    if gdalfun.gdal_clip(self.fn, tmp_clip, **clip_args)[1] == 0:
                        os.replace(tmp_clip, self.fn)
                        if self.verbose:
                            utils.echo_msg(
                                f'{tmp_clip} -> {self.fn}'
                            )
                        self.initialize()
                        if self.verbose:
                            utils.echo_msg(
                                f'clipped data with {clip_str}.'
                            )
                    else:
                        utils.echo_error_msg(
                            f'clip failed, {clip_str}'
                        )
                else:
                    utils.echo_error_msg(
                        f'could not read {clip_args["src_ply"]}'
                    )
                    
            else:
                utils.echo_warning_msg(
                    ('could not find clip ogr source/clip '
                     f'keyword {clip_args["src_ply"]}')
                )

                
    def cut(self, region=None, node='grid'):
        if region is not None:
            _tmp_cut, cut_status = gdalfun.gdal_cut_trans(
                self.fn,
                region,
                utils.make_temp_fn(
                    '__tmp_cut__.tif',
                    temp_dir=self.cache_dir
                ),
                node=node,
                co=self.co,
                verbose=self.verbose
            )
            if cut_status == 0:
                os.replace(_tmp_cut, self.fn)
                self.initialize()

                if self.verbose:
                    utils.echo_msg(
                        f'cut data to {region}...'
                    )
            else:
                utils.echo_error_msg(
                    f'failed to cut data to {region}...'
                )


    def expand_for_mask(self, stack_fn=None, cells=10):

        ## expand/contract the stack to determine mask
        if stack_fn is not None:
            with gdalfun.gdal_datasource(stack_fn) as stack_ds:
                stack_infos = gdalfun.gdal_infos(stack_ds)
                this_band = stack_ds.GetRasterBand(1)
                this_array = this_band.ReadAsArray()
                mask_array = this_array != stack_infos['ndv']
        
                mask_array = utils.expand_for(
                    mask_array, shiftx=cells, shifty=cells
                )

                mask_array = utils.expand_for(
                    mask_array, shiftx=cells*-1, shifty=cells*-1
                )

            with gdalfun.gdal_datasource(self.fn, update=True) as dem_ds:
                dem_infos = gdalfun.gdal_infos(dem_ds)
                this_band = dem_ds.GetRasterBand(1)
                this_array = this_band.ReadAsArray()
                this_array[~mask_array] = dem_infos['ndv']

                this_band.WriteArray(this_array)
                #this_band.FlushCache()                

                
    def set_interpolation_limits(
            self,
            stack_fn=None,
            size_limit=None,
            proximity_limit=None,
            percentile_limit=None,
            expand_limit=None
    ):
        """set interpolation limits"""

        if stack_fn is not None:
            ## optionally mask nodata zones based on proximity
            ## todo: adjust proxmimity grid to only limit areas
            ##       where all neighbors are below proximity limit...
            if proximity_limit is not None:
                tmp_prox = gdalfun.gdal_proximity(
                    stack_fn, utils.make_temp_fn('_prox.tif'), band=2
                )
                mn = gdalfun.gdal_get_array(tmp_prox)[0]
                if mn is not None:
                    with gdalfun.gdal_datasource(self.fn, update=True) as src_ds:
                        if src_ds is not None:
                            src_config = gdalfun.gdal_infos(src_ds)
                            src_band = src_ds.GetRasterBand(1)
                            src_array = src_band.ReadAsArray()
                            src_array[mn >= proximity_limit] = src_band.GetNoDataValue()
                            src_band.WriteArray(src_array)

                    utils.remove_glob(tmp_prox)

                    if self.verbose:
                        utils.echo_msg(
                            f'set proximity interpolation limit to {proximity_limit}.'
                        )
                else:
                    if self.verbose:
                        utils.echo_warning_msg(
                            'could not set proximity limit'
                        )
                        
            ## optionally mask nodata zones based on size_threshold
            if size_limit is not None:
                mn = nodata_count_mask(stack_fn, band=2)
                if mn is not None:
                    with gdalfun.gdal_datasource(self.fn, update=True) as src_ds:
                        if src_ds is not None:
                            src_config = gdalfun.gdal_infos(src_ds)
                            src_band = src_ds.GetRasterBand(1)
                            src_array = src_band.ReadAsArray()
                            src_array[mn >= size_limit] = src_band.GetNoDataValue()
                            src_band.WriteArray(src_array)

                    if self.verbose:
                        utils.echo_msg(
                            f'set size interpolation limit to {size_limit}.'
                        )
                else:
                    if self.verbose:
                        utils.echo_warning_msg('could not set size limit')
                        
            ## optionally mask nodata zones based on percentile_threshold
            if percentile_limit is not None:
                mn = nodata_count_mask(stack_fn, band=2)
                if mn is not None:
                    mn[mn == 0] = np.nan
                    mn_percentile = np.nanpercentile(mn, percentile_limit)
                    with gdalfun.gdal_datasource(self.fn, update=True) as src_ds:
                        if src_ds is not None:
                            src_config = gdalfun.gdal_infos(src_ds)
                            src_band = src_ds.GetRasterBand(1)
                            src_array = src_band.ReadAsArray()
                            src_array[mn >= mn_percentile] = src_band.GetNoDataValue()
                            src_band.WriteArray(src_array)

                    if self.verbose:
                        utils.echo_msg(
                            f'set percentile interpolation limit to {mn_percentile}.'
                        )
                else:
                    if self.verbose:
                        utils.echo_warning_msg('could not set percentile limit')

            ## optionally mask nodata zones based on size_threshold
            if expand_limit is not None:
                #mn = nodata_count_mask(stack_fn, band=2)
                #if mn is not None:
                mask_array = None
                with gdalfun.gdal_datasource(stack_fn) as stack_ds:
                    stack_infos = gdalfun.gdal_infos(stack_ds)
                    this_band = stack_ds.GetRasterBand(1)
                    this_array = this_band.ReadAsArray()
                    mask_array = this_array != stack_infos['ndv']

                    mask_array = utils.expand_for(mask_array, shiftx=expand_limit, shifty=expand_limit)
                    mask_array = np.invert(utils.expand_for(np.invert(mask_array), shiftx=expand_limit, shifty=expand_limit))
                    
                if mask_array is not None:
                    with gdalfun.gdal_datasource(self.fn, update=True) as src_ds:
                        if src_ds is not None:
                            src_config = gdalfun.gdal_infos(src_ds)
                            src_band = src_ds.GetRasterBand(1)
                            src_array = src_band.ReadAsArray()
                            src_array[~mask_array] = src_band.GetNoDataValue()
                            src_band.WriteArray(src_array)

                    if self.verbose:
                        utils.echo_msg(
                            f'set expandede interpolation limit to {expand_limit}.'
                        )
                else:
                    if self.verbose:
                        utils.echo_warning_msg('could not set expand limit')
                        
                        
    def set_limits(self, upper_limit=None, lower_limit=None, band=1):
        ## limit in other bands?? or chose band to limit??
        upper_limit = utils.float_or(upper_limit)
        lower_limit = utils.float_or(lower_limit)
        if upper_limit is not None or lower_limit is not None:
            dem_ds = gdal.Open(self.fn, 1)
            if dem_ds is not None:
                src_band = dem_ds.GetRasterBand(band)
                band_data = src_band.ReadAsArray()

                if upper_limit is not None:
                    band_data[band_data > upper_limit] = upper_limit
                    if self.verbose:
                        utils.echo_msg(
                            f'set upper limit to {upper_limit}.'
                        )

                if lower_limit is not None:
                    band_data[band_data < lower_limit] = lower_limit
                    if self.verbose:
                        utils.echo_msg(
                            f'set lower limit to {lower_limit}.'
                        )

                src_band.WriteArray(band_data)
                
                dem_ds = None
                self.initialize()

                
    def set_srs(self, dst_srs=None):
        if dst_srs is not None:
            gdalfun.gdal_set_srs(
                self.fn, src_srs=dst_srs, verbose=self.verbose
            )
            self.initialize()
            if self.verbose:
                utils.echo_msg(
                    f'set SRS to {dst_srs}...'
                )

                
    def reformat(self, out_fmt=None):
        if out_fmt is not None:
            if out_fmt != self.ds_config['fmt']:
                out_fn = '{}.{}'.format(
                    utils.fn_basename2(self.fn), gdalfun.gdal_fext(out_fmt)
                )
                out_ds = gdal.GetDriverByName(out_fmt).CreateCopy(
                    out_fn, 0, co=self.co
                )
                if out_ds is not None:
                    utils.remove_glob(self.fn)
                    self.fn = out_fn
                    self.initialize()

                    if self.verbose:
                        utils.echo_msg(
                            f'formatted data to {out_fmt}.'
                        )

                        
    def move(self, out_fn=None, out_dir=None):
        if out_fn is not None:
            _out_fn = os.path.join(os.path.dirname(self.fn), out_fn)
            os.replace(self.fn, _out_fn)                
            self.fn = _out_fn
            self.initialize()
            if self.verbose:
                utils.echo_msg(
                    'moved output DEM from {} to {}.'.format(
                        os.path.basename(self.fn), os.path.abspath(out_fn)
                    )
                )
            
        if out_dir is not None:
            _out_fn = os.path.join(out_dir, os.path.basename(self.fn))
            os.replace(self.fn, _out_fn)
            self.fn = _out_fn
            self.initialize()
            if self.verbose:
                utils.echo_msg(
                    'moved output DEM {} to {}.'.format(
                        os.path.basename(self.fn), os.path.abspath(out_dir)
                    )
                )

                
    def set_metadata(self, cudem=False, node='pixel'):
        """add metadata to the waffled raster

        Args: 
          cudem (bool): add CUDEM metadata
        """
        
        dem_ds = gdal.Open(self.fn, 1)
        if dem_ds is not None:
            md = self.ds_config['metadata']
            md['AREA_OR_POINT'] = 'Point'
            dem_ds.SetMetadata(md)
            dem_ds = None
            
        dem_ds = gdal.Open(self.fn, 1)
        if dem_ds is not None:
            md = self.ds_config['metadata']
            md['TIFFTAG_DATETIME'] = '{}'.format(utils.this_date())

            if node == 'pixel':
                md['AREA_OR_POINT'] = 'Area'
                md['NC_GLOBAL#node_offset'] = '1'
                md['tos#node_offset'] = '1'
            else:
                md['AREA_OR_POINT'] = 'Point'
                md['NC_GLOBAL#node_offset'] = '0'
                md['tos#node_offset'] = '0'

            if cudem:
                md['TIFFTAG_COPYRIGHT'] = ('DOC/NOAA/NESDIS/NCEI > National '
                                           'Centers for Environmental Information, '
                                           'NESDIS, NOAA, U.S. Department of Commerce')
                if self.ds_config['zr'][1] < 0:
                    tb = 'Bathymetry'
                elif self.ds_config['zr'][0] > 0:
                    tb = 'Topography'
                else:
                    tb = 'Topography-Bathymetry'

                srs=osr.SpatialReference(wkt=self.ds_config['proj'])
                vdatum=srs.GetAttrValue('vert_cs')
                md['TIFFTAG_IMAGEDESCRIPTION'] = '{}; {}'.format(
                    tb, '' if vdatum is None else vdatum
                )

            if dem_ds.SetMetadata(md) != 0:
                utils.echo_error_msg('failed to correctly set metadata')

            dem_ds = None
            if self.verbose:
                utils.echo_msg(
                    f'set DEM metadata: {md}.'
                )

    def stack2nc(self, stack_fn=None, mask_fn=None, unc_fn=None, node='pixel'):
        def generate_dims(fn, grp):
            dem_infos = gdalfun.gdal_infos(fn)
            nx = dem_infos['nx']
            ny = dem_infos['ny']
            geoT = dem_infos['geoT']
            proj = dem_infos['proj']
            cstype = srsfun.srs_get_cstype(proj)
            horz_epsg, vert_epsg = srsfun.epsg_from_input(proj)
            if cstype == 'GEOGCS':
                lat_name = 'lat'
                lat_long_name = 'latitude'
                lat_units = 'degrees_north'
                lon_name = 'lon'
                lon_long_name = 'longitude'
                lon_units = 'degrees_east'
            else:
                lat_name = 'y'
                lat_long_name = 'northing'
                lat_units = 'meters'
                lon_name = 'x'
                lon_long_name = 'easting'
                lon_units = 'meters'
            
            if node == 'pixel': # pixel-node
                lon_start = geoT[0] + (geoT[1] / 2)
                lat_start = geoT[3] + (geoT[5] / 2)
                grp.GDAL_AREA_OR_POINT = 'AREA'
                grp.node_offset = 1
            else: # grid-node
                lon_start = geoT[0]
                lat_start = geoT[3]
                grp.GDAL_AREA_OR_POINT = 'POINT'
                grp.node_offset = 0

            lon_end = geoT[0] + (geoT[1] * nx)
            lat_end = geoT[3] + (geoT[5] * ny)
            lon_inc = geoT[1]
            lat_inc = geoT[5]

            grp.Conventions = 'CF-1.5'
            #nc_ds.node_offset = 1

            # ************  Grid Mapping  ********
            crs_var = grp.createVariable('crs', 'S1', (), compression='zlib')
            crs_var.long_name = 'CRS definition'
            crs_var.GeoTransfrom = ' '.join([str(x) for x in geoT])
            crs_var.crs_wkt = proj
            crs_var.spatial_ref = proj
            crs_var.horz_crs_epsg = horz_epsg
            crs_var.vert_crs_epsg = 'EPSG:{}'.format(vert_epsg)

            # ************  LATITUDE  ************
            lat_array = np.arange(lat_start, lat_end, lat_inc)
            lat_dim = grp.createDimension(lat_name, ny) # lat dim
            lat_var = grp.createVariable(
                lat_name, 'd', (lat_name,), compression='zlib'
            ) # lat var
            lat_var.long_name = lat_long_name
            lat_var.standard_name = lat_long_name
            lat_var.units = lat_units
            lat_var.actual_range = [lat_array.min(), lat_array.max()]
            lat_var[:] = lat_array

            # ************  LONGITUDE  ***********
            lon_array = np.arange(lon_start, lon_end, lon_inc)
            lon_dim = grp.createDimension(lon_name, nx) # lon dim
            lon_var = grp.createVariable(
                lon_name, 'd', (lon_name,), compression='zlib'
            ) # lon_var
            lon_var.long_name = lon_long_name
            lon_var.standard_name = lon_long_name
            lon_var.units = lon_units
            lon_var.actual_range = [lon_array.min(), lon_array.max()]
            lon_var[:] = lon_array

            return(lat_name, lon_name, 'crs')
            
        if self.verbose:
            utils.echo_msg('Generateing NetCDF output...')

        # ************  DEM GRP **************
        nc_ds = nc.Dataset(
            os.path.join(
                os.path.dirname(self.fn),
                '{}.nc'.format(os.path.basename(utils.fn_basename2(self.fn)))
            ),
            'w',
            format='NETCDF4'
        )
        #dem_grp = nc_ds.createGroup('dem')    
        lat, lon, crs = generate_dims(self.fn, nc_ds)
        dem_infos = gdalfun.gdal_infos(self.fn, scan=True)
        horz_epsg, vert_epsg = srsfun.epsg_from_input(dem_infos['proj'])
        dem_var = nc_ds.createVariable(
            'z', 'f4', (lat, lon,), compression='zlib', fill_value=dem_infos['ndv']
        ) # dem_h var
        dem_var.units = 'meters'
        dem_var.grid_mapping = 'crs'
        dem_var.long_name = 'z'
        dem_var.standard_name = 'height'
        dem_var.positive = 'up'
        dem_var.vert_crs_epsg = 'EPSG:{}'.format(vert_epsg)
        dem_var.actual_range = dem_infos['zr']
        
        with gdalfun.gdal_datasource(self.fn) as dem_ds:
            dem_band = dem_ds.GetRasterBand(1)
            dem_var[:] = dem_band.ReadAsArray()

        # ************  UNC GRP **************
        if unc_fn is not None:
            if os.path.exists(unc_fn):
                unc_grp = nc_ds.createGroup('uncertainty')
                lat, lon, crs = generate_dims(unc_fn, unc_grp)
                with gdalfun.gdal_datasource(unc_fn) as unc_ds:
                    if unc_ds is not None:
                        for b in range(1, unc_ds.RasterCount+1):
                            unc_infos = gdalfun.gdal_infos(unc_ds, scan=True, band=b)
                            unc_band = unc_ds.GetRasterBand(b)
                            unc_name = unc_band.GetDescription()
                            unc_var = unc_grp.createVariable(
                                unc_name,
                                'f4',
                                (lat, lon,),
                                compression='zlib',
                                fill_value=unc_infos['ndv']
                            )
                            unc_var[:] = unc_band.ReadAsArray()
                            if unc_name != 'proximity':
                                unc_var.units = 'meters'
                            else:
                                unc_var.units = 'cells'

                            unc_var.grid_mapping = "crs"
                            unc_var.actual_range = unc_infos['zr']
            else:
                utils.echo_warning_msg(
                    f'{unc_fn} does not exist'
                )
            
        # ************  MSK GRP ****************
        if mask_fn is not None:
            if os.path.exists(mask_fn):
                mask_grp = nc_ds.createGroup('mask')
                lat, lon, crs = generate_dims(mask_fn, mask_grp)
                with gdalfun.gdal_datasource(mask_fn) as mask_ds:
                    if mask_ds is not None:
                        for b in range(1, mask_ds.RasterCount+1):
                            mask_infos = gdalfun.gdal_infos(mask_ds, scan=True, band=b)
                            mask_band = mask_ds.GetRasterBand(b)
                            mask_var = mask_grp.createVariable(
                                mask_band.GetDescription(),
                                'i4',
                                (lat, lon,),
                                compression='zlib',
                                fill_value=mask_infos['ndv']
                            )
                            mask_var.grid_mapping = "crs"
                            mask_var[:] = mask_band.ReadAsArray()
                            this_band_md = mask_band.GetMetadata()
                            for k in this_band_md.keys():
                                try:
                                    mask_var.__setattr__(
                                        os.path.basename(k), os.path.basename(this_band_md[k])
                                    )
                                except:
                                    mask_var.__setattr__(
                                        '_{}'.format(
                                            os.path.basename(k)
                                        ),
                                        os.path.basename(this_band_md[k])
                                    )
                                    
                            mask_var.actual_range = mask_infos['zr']
            else:
                utils.echo_warning_msg(
                    f'{msk_fn} does not exist'
                )

        # ************  STACK  ***************
        if stack_fn is not None:
            if os.path.exists(stack_fn):
                stack_grp = nc_ds.createGroup('stack')
                lat, lon, crs = generate_dims(stack_fn, stack_grp)
                with gdalfun.gdal_datasource(stack_fn) as stack_ds:
                    for b, n in enumerate(
                            ['stack_h', 'count', 'weight', 'uncertainty',
                             'src_uncertainty', 'x', 'y']
                    ):
                        stack_infos = gdalfun.gdal_infos(stack_ds, scan=True, band=b+1)
                        this_band = stack_ds.GetRasterBand(b+1)
                        this_var = stack_grp.createVariable(
                            n,
                            'f4',
                            (lat, lon,),
                            compression='zlib',
                            fill_value=stack_infos['ndv']
                        )
                        this_var[:] = this_band.ReadAsArray()
                        this_var.short_name = n
                        this_var.actual_range = stack_infos['zr']
                        this_var.grid_mapping = "crs"
            else:
                utils.echo_warning_msg(
                    f'{stack_fn} does not exist'
                )
                                
        nc_ds.close()

        
    def stack2hd5(self, stack_fn=None, mask_fn=None, unc_fn=None, node='pixel'):
        """put all the data into an hdf5 file"""

        if self.verbose:
            utils.echo_msg('Generating H5 output...')
        
        f = h5.File(
            os.path.join(
                os.path.dirname(self.fn),
                '{}.h5'.format(
                    os.path.basename(utils.fn_basename2(self.fn))
                )
            ),
            'w'
        )
        dem_infos = gdalfun.gdal_infos(self.fn)
        nx = dem_infos['nx']
        ny = dem_infos['ny']
        geoT = dem_infos['geoT']

        if node == 'pixel': # pixel-node
            lon_start = geoT[0] + (geoT[1] / 2)
            lat_start = geoT[3] + (geoT[5] / 2)
        else: # grid-node
            lon_start = geoT[0]
            lat_start = geoT[3]

        lon_end = geoT[0] + (geoT[1] * nx)
        lat_end = geoT[3] + (geoT[5] * ny)
        lon_inc = geoT[1]
        lat_inc = geoT[5]
        
        # ************  DEM  GRP *************
        dem_grp = f.create_group('dem')
        
        # ************  STACK  GRP ***********
        if stack_fn is not None:
            stack_grp = f.create_group('stack')
            stack_infos = gdalfun.gdal_infos(stack_fn)

        # ************  MSK  GRP *************
        if mask_fn is not None:
            mask_grp = f.create_group('mask')

        # ************  METADATA  ************
        # f.create_dataset('geo_transform', data=geoT) # geo_transform
        # srs info
        # other metdata

        # ************  Grid Mapping  ********
        crs_dset = f.create_dataset('crs', dtype=h5.string_dtype())
        crs_dset.attrs['GeoTransfrom'] = ' '.join([str(x) for x in geoT])
        crs_dset.attrs['crs_wkt'] = self.ds_config['proj']
        
        # ************  LATITUDE  ************
        lat_array = np.arange(lat_start, lat_end, lat_inc)
        lat_dset = f.create_dataset ('lat', data=lat_array)
        lat_dset.make_scale('latitude')
        lat_dset.attrs["long_name"] = "latitude"
        lat_dset.attrs["units"] = "degrees_north"
        lat_dset.attrs["standard_name"] = "latitude"
        lat_dset.attrs["actual_range"] = [lat_end, lat_start]
        #lat_dset.attrs[''] = ''

        # ************  LONGITUDE  ***********
        lon_array = np.arange(lon_start, lon_end, lon_inc)
        lon_dset = f.create_dataset ('lon', data=lon_array)
        lon_dset.make_scale('longitude')
        lon_dset.attrs["long_name"] = "longitude"
        lon_dset.attrs["units"] = "degrees_east" 
        lon_dset.attrs["standard_name"]= "longitude"
        lon_dset.attrs["actual_range"] = [lon_start, lon_end]

        # ************  DEM  *****************
        with gdalfun.gdal_datasource(self.fn) as dem_ds:
            dem_band = dem_ds.GetRasterBand(1)
            #dem_dset = f.create_dataset ('z', data=dem_band.ReadAsArray())
            dem_dset = dem_grp.create_dataset(
                'dem_h', data=dem_band.ReadAsArray()
            )
            dem_dset.attrs['units'] = 'meters'
            dem_dset.dims[0].attach_scale(lat_dset)
            dem_dset.dims[1].attach_scale(lon_dset)
            dem_dset.attrs['grid_mapping'] = "crs"
            
        # ************  Uncertainty  *********
        if unc_fn is not None:
            with gdalfun.gdal_datasource(unc_fn) as unc_ds:
                unc_band = unc_ds.GetRasterBand(1)
                unc_dset = dem_grp.create_dataset(
                    'uncertainty', data=unc_band.ReadAsArray()
                )
                unc_dset.attrs['units'] = 'meters'
                unc_dset.dims[0].attach_scale(lat_dset)
                unc_dset.dims[1].attach_scale(lon_dset)
                unc_dset.attrs['grid_mapping'] = "crs"
                
        # ************  MASK  ****************
        if mask_fn is not None:
            with gdalfun.gdal_datasource(mask_fn) as mask_ds:
                if mask_ds is not None:
                    for b in range(1, mask_ds.RasterCount+1):
                        mask_band = mask_ds.GetRasterBand(b)
                        mask_dset = mask_grp.create_dataset(
                            mask_band.GetDescription(), data=mask_band.ReadAsArray()
                        )
                        this_band_md = mask_band.GetMetadata()
                        for k in this_band_md.keys():
                            mask_dset.attrs[k] = this_band_md[k]

                        mask_dset.dims[0].attach_scale(lat_dset)
                        mask_dset.dims[1].attach_scale(lon_dset)
                        mask_dset.attrs['grid_mapping'] = "crs"
            
        # ************  STACK  ***************
        if stack_fn is not None:
            with gdalfun.gdal_datasource(stack_fn) as stack_ds:
                for b, n in enumerate(
                        ['stack_h', 'count', 'weight',
                         'uncertainty', 'src_uncertainty']
                ):
                    this_band = stack_ds.GetRasterBand(b+1)
                    this_dset = stack_grp.create_dataset(n, data=this_band.ReadAsArray())
                    this_dset.dims[0].attach_scale(lat_dset)
                    this_dset.dims[1].attach_scale(lon_dset)
                    this_dset.attrs['short_name'] = n
                    this_dset.attrs['grid_mapping'] = "crs"
        f.close()

                
class WaffleFactory(factory.CUDEMFactory):
    """Waffles Factory Settings
    """
    
    _modules = {
        'stacks': {'name': 'stacks', 'stack': True, 'call': WafflesStacks},
        'IDW': {'name': 'IDW', 'stack': True, 'call': WafflesIDW},
        'linear': {'name': 'linear', 'stack': True, 'call': WafflesLinear},
        'cubic': {'name': 'cubic', 'stack': True, 'call': WafflesCubic},
        'nearest': {'name': 'nearest', 'stack': True, 'call': WafflesNearest},
        'gmt-surface': {'name': 'surface', 'stack': True, 'call':GMTSurface},
        'gmt-triangulate': {'name': 'triangulate','stack': True, 'call': GMTTriangulate},
        'gmt-nearneighbor': {'name': 'nearneihbor', 'stack': True, 'call': GMTNearNeighbor},
        'mbgrid': {'name': 'mbgrid', 'stack': True, 'call': WafflesMBGrid},
        'gdal-linear': {'name': 'linear', 'stack': True, 'call': GDALLinear},
        'gdal-nearest': {'name': 'nearest', 'stack': True, 'call': GDALNearest},
        'gdal-average': {'name': 'average', 'stack': True, 'call': GDALMovingAverage},
        'gdal-invdst': {'name': 'invdst', 'stack': True, 'call': GDALInvDst},
        'vdatum': {'name': 'vdatum', 'stack': False, 'call': WafflesVDatum},
        'coastline': {'name': 'coastline', 'stack': False, 'call': WafflesCoastline},
        'lakes': {'name': 'lakes', 'stack': False, 'call': WafflesLakes},
        'cudem': {'name': 'cudem', 'stack': True, 'call': WafflesCUDEM},
        'uncertainty': {'name': 'uncertainty', 'stack': True, 'call': WafflesUncertainty},
        'scratch': {'name': 'scratch', 'stack': True, 'call': WafflesScratch},
        'flatten': {'name': 'flatten', 'stack': True, 'call': WafflesFlatten},
        ## testing
        #'num': {'name': 'num', 'stack': True, 'call': WafflesNum}, # defunct
        #'patch': {'name': 'patch', 'stack': True, 'call': WafflesPatch}, # test
        #'cube': {'name': 'cube', 'stack': True, 'call': WafflesCUBE}, # test
        #'bgrid': {'name': 'bgrid', 'stack': True, 'call': WafflesBGrid}, # test            
    }

    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        
    def _set_modules(self):
        pass

    
def waffle_queue(q):
    """waffle queue for threads

    waffles will run multiple input regions as
    separate threads...

    waffles_args is [waffle_module]
    """
    
    while True:
        waffle_module = q.get()
        if waffle_module is None:
            break
        else:
            try:
                waffle_module[0]()
            except Exception as e:
                utils.echo_error_msg(
                    f'failed to generate {waffle_module}, {e}'
                )
                print(traceback.format_exc())
                pass

            
## ==============================================
## Command-line Interface (CLI)
## $ waffles
##
## waffles cli
## ==============================================
waffles_cli_usage = lambda: """{cmd} ({wf_version}): Generate DEMs and derivatives.

usage: {cmd} [OPTIONS] DATALIST

Options:
  -R, --region\t\t\tSpecifies the desired REGION;
\t\t\t\tWhere a REGION is xmin/xmax/ymin/ymax[/zmin/zmax[/wmin/wmax]]
\t\t\t\tUse '-' to indicate no bounding range; e.g. -R -/-/-/-/-10/10/1/-
\t\t\t\tOR an OGR-compatible vector file with regional polygons. 
\t\t\t\tWhere the REGION is /path/to/vector[:zmin/zmax[/wmin/wmax]].
\t\t\t\tIf a vector file is supplied, will use each region found therein.
\t\t\t\tOptionally, append `:pct_buffer=<value>` to buffer the region(s) by a percentage.
  -E, --increment\t\tGridding INCREMENT and RESAMPLE-INCREMENT in native units.
\t\t\t\tWhere INCREMENT is x-inc[/y-inc][:sample-x-inc/sample-y-inc]
  -M, --module\t\t\tDesired Waffles MODULE and options. (see available Modules below)
\t\t\t\tWhere MODULE is module[:mod_opt=mod_val[:mod_opt1=mod_val1[:...]]]
  -S, --sample_alg\t\tReSAMPLE algorithm to use (from gdalwarp)
\t\t\t\tSet as 'auto' to use 'average' when down-sampling and 'bilinear' when up-sampling
\t\t\t\tThis switch controls resampling of input raster datasets as well as resampling
\t\t\t\tthe final DEM if RESAMPLE-INCREMENT is set in -E
  -X, --extend\t\t\tNumber of cells with which to EXTEND the output DEM REGION and a 
\t\t\t\tpercentage to extend the processing REGION.
\t\t\t\tWhere EXTEND is dem-extend(cell-count)[:processing-extend(percentage)]
\t\t\t\te.g. -X6:10 to extend the DEM REGION by 6 cells and the processing region by 10
\t\t\t\tpercent of the input REGION.
  -T, --filter\t\t\tFILTER the output DEM using one or multiple filters. 
\t\t\t\tWhere FILTER is fltr_name[:opts] (see `grits --modules` for more information)
\t\t\t\tThe -T switch may be set multiple times to perform multiple filters.
\t\t\t\tAppend `:stacks=True` to the filter to perform the filter on the data stack instead 
\t\t\t\tof the final DEM.
\t\t\t\tAvailable FILTERS: {grits_modules}
  -L, --limits\t\t\tLIMIT the output elevation or interpolation values, append 
\t\t\t\t'u<value>' to set the upper elevation limit, 
\t\t\t\t'l<value>' to set the lower elevation limit,
\t\t\t\t'n<value>' to set the count per cell limit (pre-interpolation).
\t\t\t\t'r<value>' to set an uncertainty upper limit (pre-interpolation).
\t\t\t\t'p<value>' to set an interpolation limit by proximity, or 
\t\t\t\t's<value>' to set an interpolation limit by size, or
\t\t\t\t'c<value>' to set an interpolation limit by nodata-size percentile, or
\t\t\t\t'e<value>' to set an interpolation limit by expanded/contracted size.
\t\t\t\te.g. -Lu0 to set all values above 0 to zero, or 
\t\t\t\t-Ls100 to limit interpolation to nodata zones smaller than 100 pixels, or
\t\t\t\t-Ln2 to limit stacked cells to those with at least 2 contributing data records.
  -C, --clip\t\t\tCLIP the output to the clip polygon -C<clip_ply.shp:invert=False>
  -K, --chunk\t\t\tGenerate the DEM in CHUNKs.
  -F, --format\t\t\tOutput grid FORMAT. [GTiff]
  -O, --output-basename\t\tBASENAME for all OUTPUTs.
\t\t\t\tnote: Things will get appended to this OUTPUT_BASENAME, such as output specifiers and format extensions.
  -P, --t_srs\t\t\tProjection of REGION and output DEM.
  -N, --nodata\t\t\tThe NODATA value of output DEM.
  -A, --stack-mode\t\tSet the STACK MODE to 'mean', 'min', 'max', 'mixed' or 'supercede' (higher weighted data supercedes lower weighted data).
  -G, --wg-config\t\tA waffles config JSON file. If supplied, will overwrite all other options.
\t\t\t\tGenerate a waffles_config JSON file using the --config flag.
  -H, --threads\t\t\tSet the number of THREADS (1). Each input region will be run in up to THREADS threads. 
  -D, --cache-dir\t\tCACHE Directory for storing temp data.
\t\t\t\tDefault Cache Directory is ~/.cudem_cache; cache will be cleared after a waffles session.
\t\t\t\tto retain the data, use the --keep-cache flag.
  -CO, --creation-options\tGDAL CREATION OPTIONS to use in raster creation.

  -f, --transform\t\tTransform all data to PROJECTION value set with --t_srs/-P where applicable.
  -p, --prefix\t\t\tAppend OUTPUT_BASENAME (-O) with PREFIX (append <RES>_nYYxYY_wXXxXX_<YEAR>v<VERSION> info to OUTPUT_BASENAME).
\t\t\t\tnote: Set Resolution, Year and Version by setting this as '-p res=X:year=XXXX:version=X', 
\t\t\t\tleave blank for default of <INCREMENT>, <CURRENT_YEAR> and <1>, respectively.
  -r, --grid-node\t\tUse grid-node registration, default is pixel-node.
  -w, --want-weight\t\tUse weights provided in the datalist to weight overlapping data.
  -u, --want-uncertainty\tGenerate/Use uncertainty either calculated or provided in the datalist.
  -m, --want-mask\t\tMask the processed datalist.
  -a, --archive\t\t\tARCHIVE the datalist to the given region.
  -k, --keep-cache\t\tKEEP the cache data intact after run
  -x, --keep-auxiliary\t\tKEEP the auxiliary rasters intact after run (mask, uncertainty, weights, count).
  -s, --spatial-metadata\tGenerate SPATIAL-METADATA.
  -t, --flatten-nodata-values\tFlatten any remaining nodata values.
  -c, --continue\t\tDon't clobber existing files.
  -q, --quiet\t\t\tLower verbosity to a quiet.

  --help\t\t\tPrint the usage text
  --config\t\t\tSave the waffles config JSON and major datalist
  --modules\t\t\tDisplay the module descriptions and usage
  --version\t\t\tPrint the version information

Datalists and data formats:
  A datalist is a file that contains a number of datalist entries, 
  while an entry is a space-delineated line:
  `path [format weight uncertainty [name source date type resolution hdatum vdatum url]]`

Supported datalist formats: 
  {dl_formats}

Modules (see waffles --modules <module-name> for more info):
  {modules}
""".format(
    cmd=os.path.basename(sys.argv[0]),
    dl_formats=factory._cudem_module_name_short_desc(dlim.DatasetFactory._modules),
    modules=factory._cudem_module_short_desc(WaffleFactory._modules),
    wf_version=cudem.__version__,
    grits_modules=factory._cudem_module_short_desc(grits.GritsFactory._modules)
)


def waffles_cli(argv = sys.argv):
    """run waffles from command-line

    See `waffles_cli_usage` for full cli options.
    """

    i = 1
    dls = []
    i_regions = []
    these_regions = []
    module = 'stacks'
    wg_user = None
    want_prefix = False
    prefix_args = {}
    want_config = False
    keep_cache = False
    status = 0
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
    wg['flatten_nodata_values'] = False
    wg['co'] = []

    #waffle_q = queue.Queue()
    processes=[]
    waffle_q = mp.Queue()
    n_threads = 1
        
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
                    
        elif arg == '--output_basename' or arg == '-O':
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
        ## update cache_dir to default to current utils.cache_dir or
        ## if an arg, dont add .cudem_cache!
        elif arg == '--cache-dir' or arg == '-D' or arg == '-cache-dir':
            wg['cache_dir'] = os.path.join(
                utils.str_or(argv[i + 1], os.path.expanduser('~')),
                '.cudem_cache'
            )
            i = i + 1
        elif arg[:2] == '-D': wg['cache_dir'] = os.path.join(
                utils.str_or(arg[2:], os.path.expanduser('~')),
                '.cudem_cache'
        )
        elif arg == '--nodata' or arg == '-N' or arg == '-ndv':
            wg['ndv'] = utils.float_or(argv[i + 1], -9999)
            i = i + 1
        elif arg[:2] == '-N': wg['ndv'] = utils.float_or(arg[2:], -9999)

        elif arg == '--limits' or arg == '-L':
            this_limit = argv[i + 1]
            if this_limit.startswith('u'):
                wg['upper_limit'] = utils.float_or(this_limit[1:])
            elif this_limit.startswith('l'):
                wg['lower_limit'] = utils.float_or(this_limit[1:])
            elif this_limit.startswith('p'):
                wg['proximity_limit'] = utils.int_or(this_limit[1:])
            elif this_limit.startswith('s'):
                wg['size_limit'] = utils.int_or(this_limit[1:])
            elif this_limit.startswith('c'):
                wg['percentile_limit'] = utils.float_or(this_limit[1:])
            elif this_limit.startswith('e'):
                wg['expand_limit'] = utils.float_or(this_limit[1:])
            elif this_limit.startswith('n'):
                wg['count_limit'] = utils.float_or(this_limit[1:])
            elif this_limit.startswith('r'):
                wg['uncertainty_limit'] = utils.float_or(this_limit[1:])
                
            i = i + 1
        elif arg[:2] == '-L':
            this_limit = arg[2:]
            if this_limit.startswith('u'):
                wg['upper_limit'] = utils.float_or(this_limit[1:])
            elif this_limit.startswith('l'):
                wg['lower_limit'] = utils.float_or(this_limit[1:])
            elif this_limit.startswith('p'):
                wg['proximity_limit'] = utils.int_or(this_limit[1:])
            elif this_limit.startswith('s'):
                wg['size_limit'] = utils.int_or(this_limit[1:])
            elif this_limit.startswith('c'):
                wg['percentile_limit'] = utils.float_or(this_limit[1:])
            elif this_limit.startswith('e'):
                wg['expand_limit'] = utils.float_or(this_limit[1:])
            elif this_limit.startswith('n'):
                wg['count_limit'] = utils.float_or(this_limit[1:])
            elif this_limit.startswith('r'):
                wg['uncertainty_limit'] = utils.float_or(this_limit[1:])
                
        elif arg == '-threads' or arg == '--threads' or arg == '-H':
            n_threads = utils.int_or(argv[i + 1], 1)
            i = i + 1

        elif arg[:2] == '-H':
            n_threads = utils.int_or(arg[2:], 1)

        elif arg == '-stack-mode' or arg == '--stack-mode' or arg == '-A':
            wg['stack_mode'] = utils.str_or(argv[i + 1], 'mean')
            i = i + 1

        elif arg[:2] == '-A':
            wg['stack_mode'] = utils.str_or(arg[2:], 'mean')

        elif arg == '--creation-options' or arg == '-CO' or arg == '-co':
            wg['co'].append(argv[i + 1])
            i = i + 1
        elif arg[:3] == '-CO' or arg[:3] == '-co':
            wg['co'].append(arg[3:])
            
        elif arg == '--transform' or arg == '-f' or arg == '-transform':
            wg['srs_transform'] = True
            if wg['dst_srs'] is None:
                wg['dst_srs'] = 'epsg:4326'
                
        elif arg == '-w' or arg == '--want-weight':
            wg['want_weight'] = True
            
        elif arg == '-u' or arg == '--want-uncertainty':
            wg['want_uncertainty'] = True
            
        elif arg == '-p' or arg == '--prefix':
            want_prefix = True
            try:
                prefix_opts = argv[i + 1].split(':')
                prefix_args = factory.args2dict(prefix_opts, prefix_args)
                if len(prefix_args) > 0:
                    i += 1
            except:
                pass

        elif arg == '--want-mask' or arg == '--mask' or arg == '-m':
            #utils.echo_warning_msg('want-mask is depreciated, use --keep-auxiliary instead')
            wg['want_mask'] = True
            
        elif arg == '-k' or arg == '--keep-cache': keep_cache = True
        elif arg == '-x' or arg == '--keep-auxiliary': wg['keep_auxiliary'] = True
        #elif arg == '-t' or arg == '--threads': want_threads = True
        elif arg == '-a' or arg == '--archive': wg['archive'] = True
        elif arg == '-s' or arg == '--spatial-metadata':
            wg['want_mask'] = True
            wg['want_sm'] = True
            
        elif arg == '-c' or arg == '--continue': wg['clobber'] = False
        elif arg == '-r' or arg == '--grid-node': wg['node'] = 'grid'
        elif arg == '-t' or arg == '--flatten-ndata-values': wg['flatten_nodata_values'] = True

        elif arg == '--quiet' or arg == '-q': wg['verbose'] = False
        elif arg == '--config': want_config = True
        elif arg == '--modules':
            factory.echo_modules(
                WaffleFactory._modules,
                None if i+1 >= len(argv) else sys.argv[i+1]
            )
            sys.exit(0)
        elif arg == '--help' or arg == '-h':
            sys.stderr.write(waffles_cli_usage())
            sys.exit(0)
        elif arg == '--version' or arg == '-v':
            sys.stdout.write('{}\n'.format(cudem.__version__))
            sys.exit(0)
        elif arg[0] == '-':
            sys.stderr.write(waffles_cli_usage())
            utils.echo_error_msg(
                f'{arg} is not a valid waffles cli switch'
            )
            sys.exit(0)
        else: dls.append(arg)
        i += 1

    for _ in range(n_threads):
        #t = threading.Thread(target=waffle_queue, args=([waffle_q]))
        t = mp.Process(target=waffle_queue, args=([waffle_q]))
        processes.append(t)
        t.daemon = True
        t.start()

    ## load the user wg json and run waffles with that.
    ## allow input of multiple config files with -G .. -G ..
    if wg_user is not None:
        if os.path.exists(wg_user):
            with open(wg_user, 'r') as wgj:
                wg = json.load(wgj)
                # if wg['kwargs']['src_region'] is not None and \
                #    isinstance(wg['kwargs']['src_region'], list):
                #     wg['kwargs']['src_region'] = regions.Region().from_list(
                #         wg['kwargs']['src_region']
                #     )

            this_waffle = WaffleFactory(
                mod=None,
                mod_name=wg['mod_name'],
                mod_args=wg['mod_args'],
                **wg['kwargs']
            )
            this_waffle_module = this_waffle._acquire_module()
            waffle_q.put([this_waffle_module])
        else:
            utils.echo_error_msg(
                f'specified waffles config file does not exist, {wg_user}'
            )
        waffle_q.put(None)
        [t.join() for t in processes]
        #waffle_q.join()
        sys.exit(0)

    ## Otherwise run from cli options...
    ## set the dem module
    if module.split(':')[0] not in WaffleFactory()._modules.keys():
        utils.echo_error_msg(
            '''{} is not a valid waffles module, available modules are: {}'''.format(
                module.split(':')[0], factory._cudem_module_short_desc(WaffleFactory._modules)
            )
        )
        sys.exit(-1)
        
    if WaffleFactory()._modules[module.split(':')[0]]['stack']:
        if len(dls) == 0:
            sys.stderr.write(waffles_cli_usage())
            utils.echo_error_msg(
                '''must specify a datalist/entry, try `gmrt` or `srtm` for global data.'''
            )
            sys.exit(-1)
    else:
        wg['want_stack'] = True if len(dls) > 0 else False

    ## check the increment
    if 'xinc' in wg.keys():
        if wg['xinc'] is None:
            sys.stderr.write(waffles_cli_usage())
            utils.echo_error_msg('''must specify a gridding increment.''')
            sys.exit(-1)
        else:
            if wg['yinc'] is None:
                wg['yinc'] = wg['xinc']
    else:
        sys.stderr.write(waffles_cli_usage())
        utils.echo_error_msg('''must specify a gridding increment.''')
        sys.exit(-1)      

    ## set the datalists and names
    wg['data'] = dls
    if not i_regions: i_regions = [None]
    these_regions = regions.parse_cli_region(i_regions, wg['verbose'])
    name = wg['name']

    ## parse the regions and add them to the queue, or output as config file
    for i, this_region in enumerate(these_regions):
        ## input region is None, so gather the region from the input data...
        if this_region is None:
            utils.echo_warning_msg(
                'No input region specified, gathering region from input data...'
            )
            this_datalist = dlim.init_data(
                dls,
                region=this_region,
                dst_srs=wg['dst_srs'],
                want_verbose=wg['verbose']
            )
            if this_datalist is not None and this_datalist.valid_p(
                    fmts=dlim.DatasetFactory._modules[this_datalist.data_format]['fmts']
            ):
                this_datalist.initialize()
                this_inf = this_datalist.inf()
                this_region = regions.Region().from_list(this_inf.minmax)
                if wg['dst_srs'] is not None:
                    if this_inf.src_srs is not None:
                        this_region.src_srs = this_inf.src_srs
                        this_region.warp(dst_srs)

            utils.echo_msg('region is {}'.format(this_region))
            if this_region is None: # couldn't gather a region from the data
                break

        wg['src_region'] = this_region

        ## set the output name, appending the region, etc. if wanted/needed
        if want_prefix or len(these_regions) > 1:
            wg['name'] = utils.append_fn(
                name,
                wg['src_region'],
                wg['xsample'] if wg['xsample'] is not None else wg['xinc'],
                **prefix_args
            )
            
        if want_config: # export the waffles module config file
            # this_datalist = dlim.init_data(
            #     dls,
            #     region=this_region,
            #     dst_srs=wg['dst_srs'],
            #     want_verbose=wg['verbose'],
            #     want_weight=True,
            #     xy_inc=[wg['xinc'],wg['yinc']],
                
            # )
            # if this_datalist is not None and this_datalist.valid_p(
            #         fmts=dlim.DatasetFactory._modules[this_datalist.data_format]['fmts']
            # ):
            #     this_datalist.initialize()
            #     wg['data'] = []
            #     for this_entry in this_datalist.parse():
            #         params = this_entry.params

            #         if params['kwargs']['pnt_fltrs'] is not None and \
            #            isinstance(params['kwargs']['pnt_fltrs'], list):
            #             for fltr in params['kwargs']['pnt_fltrs']:
            #                 fltr_dict = factory.fmod2dict(fltr)
            #                 fltr_dict = utils.dict_path2abspath(fltr_dict)
            #                 params['kwargs']['pnt_fltrs'] = factory.dict2fmod(fltr_dict)
                    
            #         params['kwargs'].pop('params')
            #         params['kwargs'].pop('parent')
            #         params['kwargs']['src_region'] = \
            #             params['kwargs']['src_region'].export_as_list()
            #         params = utils.dict_path2abspath(params, except_keys=['mod'])
            #         wg['data'].append(params)
            
            wg['src_region'] = this_region.export_as_list()
            wg['name'] = os.path.abspath(wg['name'])
            this_waffle = WaffleFactory(mod=module, **wg)
            this_waffle.write_parameter_file('{}.json'.format(wg['name']))
        else: # get the waffle module and add it to the queue
            this_waffle = WaffleFactory(mod=module, **wg)
            if this_waffle is not None:
                this_waffle_module = this_waffle._acquire_module()
                if this_waffle_module is not None:
                    waffle_q.put([this_waffle_module])
                    ##this_waffle_module()
                else:
                    if wg['verbose']:
                        utils.echo_error_msg(
                            f'could not acquire waffles module {module}'
                        )

    for _ in range(n_threads):
        waffle_q.put(None)
        
    for t in processes:
        t.join()

    ## remove the cahce dir if not asked to keep
    if not keep_cache:
       utils.remove_glob(wg['cache_dir'])

    return(0)

### End
