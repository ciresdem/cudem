### waffles.py
##
## Copyright (c) 2010 - 2026 Regents of the University of Colorado
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
### Commentary:
##
## Waffles: The CUDEM Gridding and Interpolation Engine.
##
## Waffles is a unified factory interface for generating Digital Elevation Models (DEMs)
## and other geospatial grids from scattered data sources. It abstracts the complexity
## of various interpolation algorithms and data processing steps into a single,
## coherent pipeline.
##
##   * Unified Gridding Interface:
##      - Supports internal algorithms: IDW, Stacks (Weighted Mean), Linear/Cubic/Nearest,
##        Natural Neighbor, Kriging, Machine Learning Regression (RF, KNN, MLP),
##        and Computer Vision Inpainting.
##      - Wraps external tools: GMT Surface/Triangulate, MB-System mbgrid, GDAL Grid, VDatum.
##      - Specialized modules: Coastline masking, Lake filling, and Uncertainty estimation.
##
##   * Robust Data Handling (via DLIM):
##      - Ingests diverse formats: XYZ, LAS/LAZ, GeoTIFF, BAG, OGR, Fetch modules.
##      - "Stacking": Pre-processes overlapping data into a weighted intermediate raster
##        to resolve conflicts and enforce data hierarchy before interpolation.
##
##   * Scalability & Region Management:
##      - Chunking: Automatically splits large regions into tiles for parallel processing
##        and seamless stitching.
##      - Buffering: Handles region buffering to eliminate edge artifacts.
##      - Resolution: Manages independent processing (input) and sampling (output) resolutions.
##
##   * Post-Processing:
##      - Automated clipping to vector polygons.
##      - Generation of auxiliary products: Uncertainty grids, Data Masks, Spatial Metadata.
##      - Filtering: Application of 'grits' filters (smoothing, outlier removal) to the output.
##
## Usage:
##   CLI: waffles -R <region> -E <increment> -M <module> [options] input.datalist
##   API: w = WaffleFactory(mod='surface', data=['input.xyz'], ...).acquire()
##        w.generate()
##
### Code:

import sys
import os
import json
import glob
import traceback
import argparse
import signal
import numpy as np
import multiprocessing as mp
from osgeo import gdal

from cudem.datalists import dlim
from cudem.datalists import gdalfile
from cudem import pointz
from cudem import regions
from cudem import utils
from cudem import gdalfun
from cudem import factory
from cudem.grits import grits
from cudem.globato import spatial_metadata
from cudem import __version__ as __cudem_version__
from . import __version__

## Data cache directory, hold temp data, fetch data, etc here.
waffles_cache = utils.cudem_cache()
gc = utils.config_check()
gdal.UseExceptions()
gdal.SetConfigOption(
    'CPL_LOG', 'NUL' if gc['platform'] == 'win32' else '/dev/null'
)

def _run_chunk_worker(args):
    """Worker function to process a single Waffles chunk in parallel.
    """
    
    from .waffledem import WaffleDEM

    ## Unpack arguments
    (srcwin, dst_gt, xinc, yinc, params, verbose) = args

    ## Calculate Chunk Region Geometry
    ## Pixel -> Geo coordinates
    this_geo_x_origin, this_geo_y_origin = utils._pixel2geo(
        srcwin[0], srcwin[1], dst_gt
    )
    
    this_gt = [
        this_geo_x_origin, float(dst_gt[1]), 0.0,
        this_geo_y_origin, 0.0, float(dst_gt[5])
    ]
    
    this_region = regions.Region().from_geo_transform(
        geo_transform=this_gt, x_count=srcwin[2], y_count=srcwin[3]
    )
    
    ## Buffer the region by 10% (or minimum) to prevent edge artifacts
    this_region.buffer(pct=10, x_inc=xinc, y_inc=yinc)
    
    ## Configure the Sub-Waffle Module
    this_params = params.copy()
    this_params['kwargs']['src_region'] = this_region
    this_params['kwargs']['chunk'] = None
    ## Create a unique cache_dir to hold data, etc.
    this_params['kwargs']['cache_dir'] = os.path.join(this_params['kwargs']['cache_dir'], utils.append_fn('_chunk', this_region, xinc, high_res=True))
    
    ## Generate unique name for chunk to avoid collision
    this_params['kwargs']['name'] = utils.append_fn(
        '_chunk', this_region, xinc, high_res=True
    )
    
    ## Initialize and Generate
    try:
        this_waffle = WaffleFactory().load_parameter_dict(this_params)
        this_waffle_module = this_waffle._acquire_module()
        this_waffle_module.initialize()
        this_waffle_module.generate()
        
        ## Prepare Result Dictionary
        result = {
            'fn': this_waffle_module.fn,
            'stack': getattr(this_waffle_module, 'stack', None),
            'msk_fn': getattr(this_waffle_module, 'msk_fn', None),
            'aux_dems': getattr(this_waffle_module, 'aux_dems', []),
            'cache_dir': this_waffle_module.cache_dir,
            'fmt': this_waffle_module.fmt,
            'name': this_waffle_module.name,
            'valid_dem': False,
            'valid_stack': False
        }

        ## Validate Outputs
        if WaffleDEM(result['fn'], verbose=False).valid_p():
            result['valid_dem'] = True
            
        if result['stack'] and WaffleDEM(result['stack'], verbose=False).valid_p():
            result['valid_stack'] = True
            
        return result

    except Exception as e:
        if verbose:
            utils.echo_warning_msg(f"Chunk failed @ {srcwin}: {e}")
        return None
    

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
        self.fn = f'{self.name}.tif' # output dem filename
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
        return f'<Waffles: {self.name}>'

    
    def __repr__(self):
        return f'<Waffles: {self.name}>'

    
    def initialize(self):
        ## Setup the initialization 'string' to print out in a single blob.
        init_str = ['']
        init_str.append('--------------')
        init_str.append(f'Initializing waffles module < \033[1m{self.params["mod"]}\033[m >')
        
        ## Output dem filename
        self.fn = f'{self.name}.{gdalfun.gdal_fext(self.fmt)}'
        ## Cudem config file holding foriegn programs and versions
        self.gc = utils.config_check()
        ## Initialize regions
        init_str = self._init_regions(init_str=init_str)
        ## Initialize increments
        init_str = self._init_incs(init_str=init_str) 

        if isinstance(self.co, list):
            if len(self.co) == 0:
                self.co = ["COMPRESS=DEFLATE", "TILED=YES"]
                
        else:
            self.co = ["COMPRESS=DEFLATE", "TILED=YES"]

        ## Initialize data, setting set_incs to True will force dlim to process the
        ## data to the set increments
        if self.want_stack:
            self._init_data(set_incs=True)

        ## Set the DataSource Configuration dictionary
        self.xcount, self.ycount, self.dst_gt \
            = self.p_region.geo_transform(
                x_inc=self.xinc, y_inc=self.yinc, node='grid'
            )
        self.ds_config = gdalfun.gdal_set_infos(
            self.xcount, self.ycount, (self.xcount*self.ycount),
            self.dst_gt, gdalfun.osr_wkt(self.dst_srs),
            gdal.GDT_Float32, self.ndv, self.fmt, None, None
        )
        
        init_str.append(f'Output size: {self.ds_config["nx"]}/{self.ds_config["ny"]}')
        init_str.append(f'Output srs: {self.dst_srs}')
        init_str.append(f'Output basename: {self.name}')
        init_str.append(f'Data: {self.data}')
        init_str.append('--------------')

        if self.verbose:
            out_str = '\n'.join(init_str)
            utils.echo_msg(out_str)
            
        ## Use self._init in a sub-module to do custom intializations.
        self.status = self._init()
        return self

    
    def _init(self):
        return 0

    
    def __call__(self):
        self.initialize()
        if self.status == 0:
            return self.generate()
        else:
            utils.echo_warning_msg('Failed to initialize from sub-module')

            
    def _init_regions(self, init_str=[]):
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
                f'Could not parse region: {self.region}'
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

        ## Append some info to the init_str
        init_str.append(f'Input region: {self.region}')
        init_str.append(f'Output region: {self.d_region}')
        init_str.append(f'Processing region: {self.p_region}')
        init_str.append(f'Cache directory is: {self.cache_dir}')
        return init_str

    
    def _init_data(self, set_incs=False, init_str=[]):
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
                elif f.split(':')[0] in pointz.PointFilterFactory._modules.keys():
                    point_filter = pointz.PointFilterFactory(mod=f)._acquire_module()
                    if point_filter is not None:
                        point_fltr.append(f)

        self.data = dlim.init_data(
            self.data,
            src_region=self.p_region,
            src_srs=None,
            dst_srs=self.dst_srs,
            x_inc=self.xinc,
            y_inc=self.yinc,
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
            return None

        
    def _init_incs(self, init_str=[]):
        """Initialize increments
        
        xinc/yinc are the DEM increments, in native units.
        xsample/ysample set the output DEM increments, 
        in native units.
        """
        
        self.xinc = utils.str2inc(self.xinc)
        self.yinc = utils.str2inc(self.yinc)
        self.xsample = utils.str2inc(self.xsample)
        self.ysample = utils.str2inc(self.ysample)

        ## Append some info to the init_str
        init_str.append(f'Gridding increments: {self.xinc}/{self.yinc}')
        init_str.append(f'Output increments: '
                        f'{self.xsample if self.xsample is not None else self.xinc}/'
                        f'{self.ysample if self.ysample is not None else self.yinc}')
        return init_str

    
    def _coast_region(self):
        """Coastline region 
        (extended by percentage self.extend_proc)
        """

        cr = self.d_region.copy()
        return cr.buffer(pct=self.extend_proc, x_inc=self.xinc, y_inc=self.yinc)

    
    def _proc_region(self):
        """Processing region 
        (extended by percentage self.extend_proc)
        """

        pr = self.d_region.copy()
        return pr.buffer(pct=self.extend_proc, x_inc=self.xinc, y_inc=self.yinc)

    
    def _dist_region(self):
        """Distribution region 
        (extended by self.extend).
        """
        
        dr = self.region.copy()
        if self.xsample is None and self.ysample is None:
            return dr.buffer(x_bv=(self.xinc*self.extend), y_bv=(self.yinc*self.extend))
        else:
            return dr.buffer(x_bv=(self.xsample*self.extend), y_bv=(self.ysample*self.extend))

        
    def dump_xyz(self, dst_port=sys.stdout, encode=False):
        """Dump the stacked xyz data to dst_port

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
        """Run and process the WAFFLES module.

        Generate the data 'stack' and use that to run the 
        waffles module and generate the DEM.
        """

        from .waffledem import WaffleDEM
        
        if self.data is None:
            return self

        ## Handle existing files
        if os.path.exists(self.fn):
            if not self.clobber:
                utils.echo_warning_msg(f'DEM {self.fn} already exists, skipping...')
                return self
            else:
                utils.echo_warning_msg(f'DEM {self.fn} exists and will be clobbered.')
                status = gdal.GetDriverByName(self.fmt).Delete(self.fn)
                if status != 0:
                    utils.remove_glob(f'{self.fn}*')
                    ## Remove all named files
                    #utils.remove_glob(f'{self.name}*')
        else:
            if not os.path.exists(os.path.dirname(self.fn)):
                try:
                    os.makedirs(os.path.dirname(self.fn))
                except OSError: pass

        if self.chunk is not None:
            ## Generate in Chunks and Merge
            chunks = []
            stack_chunks = []
            mask_chunks = []
            aux_chunks = []
            
            ## Prepare Task List for multiprocessing
            chunk_tasks = []
            for srcwin in utils.yield_srcwin(
                    (self.ycount, self.xcount),
                    self.chunk,
                    verbose=self.verbose
            ):
                ## Bundle args for the worker
                chunk_args = (
                    srcwin, 
                    self.dst_gt, 
                    self.xinc, 
                    self.yinc, 
                    self.params, 
                    self.verbose
                )
                chunk_tasks.append(chunk_args)

            ## Execution Results
            results = []
            
            ## Check if we can/should run in parallel
            ## We cannot spawn a pool if we are already a daemon (e.g. inside another process)
            is_daemon = mp.current_process().daemon
            num_workers = min(mp.cpu_count(), len(chunk_tasks))
            
            if not is_daemon and num_workers > 1:
                if self.verbose:
                    utils.echo_msg(f"Processing {len(chunk_tasks)} chunks with {num_workers} workers...")
                    
                with mp.Pool(processes=num_workers) as pool:
                    ## Use imap_unordered for potential speed gain on results processing
                    for res in pool.imap_unordered(_run_chunk_worker, chunk_tasks):
                        if res: results.append(res)
            else:
                if self.verbose:
                    utils.echo_msg(f"Processing {len(chunk_tasks)} chunks sequentially...")
                for args in chunk_tasks:
                    res = _run_chunk_worker(args)
                    if res: results.append(res)

            ## Process Results
            for res in results:
                if res['valid_dem']:
                    chunks.append(res['fn'])
                
                if res['valid_stack']:
                    stack_chunks.append(res['stack'])
                
                ## Check Mask presence
                if res['msk_fn'] and os.path.exists(res['msk_fn']):
                    mask_chunks.append(res['msk_fn'])
                else:
                    ## Fallback check manually constructed path from worker
                    mask_name = f"{res['name']}_msk"
                    mask_fn = (f"{os.path.join(res['cache_dir'], mask_name)}."
                               f"{gdalfun.gdal_fext(res['fmt'])}")
                    if os.path.exists(mask_fn):
                        mask_chunks.append(mask_fn)

                if res['aux_dems']:
                    aux_chunks.extend(res['aux_dems'])

            ## Merge DEM Chunks
            if len(chunks) > 0:
                gdalfun.sample_warp(
                    chunks, self.fn,
                    self.xinc, self.yinc,
                    src_srs=self.dst_srs,
                    src_region=self.p_region,
                    sample_alg='near',
                    ndv=-9999,
                    co=["COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES"],
                    ot=gdal.GDT_Float32,
                    verbose=self.verbose
                )
                ## Use VRT for efficient merging
                ## TODO: use gdalfun.create_temp_vrt instead here.
                # vrt_opts = gdal.BuildVRTOptions(resampleAlg='nearest')
                # gdal.BuildVRT(f'{self.fn}.vrt', chunks, options=vrt_opts)
                
                # ## Warp VRT to final TIF                
                # g = gdal.Warp(
                #     self.fn, 
                #     f'{self.fn}.vrt',
                #     format=self.fmt,
                #     resampleAlg='cubicspline',
                #     options=["COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES"]
                # )
                # g = None
                # utils.remove_glob(f'{self.fn}.vrt')

            ## Merge Stack Chunks
            if len(stack_chunks) > 0:
                stack_fn = f'{self.name}_stack.{gdalfun.gdal_fext(self.fmt)}'
                gdalfun.sample_warp(
                    stack_chunks, stack_fn,
                    self.xinc, self.yinc,
                    src_srs=self.dst_srs,
                    src_region=self.p_region,
                    sample_alg='near',
                    ndv=-9999,
                    co=["COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES"],
                    ot=gdal.GDT_Float32,
                    verbose=self.verbose
                )
                # g = gdal.Warp(
                #     stack_fn,
                #     stack_chunks,
                #     format=self.fmt,
                #     resampleAlg='cubicspline',
                #     options=["COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES"]
                # )
                # g = None
                self.stack = stack_fn

            ## Merge Mask Chunks
            if len(mask_chunks) > 0:
                mask_fn = f'{self.name}_msk.{gdalfun.gdal_fext(self.fmt)}'
                gdalfun.sample_warp(
                    mask_chunks, stack_fn,
                    self.xinc, self.yinc,
                    src_srs=self.dst_srs,
                    src_region=self.p_region,
                    sample_alg='near',
                    ndv=-9999,
                    co=["COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES"],
                    ot=gdal.GDT_Float32,
                    verbose=self.verbose
                )
                # g = gdal.Warp(
                #     mask_fn,
                #     mask_chunks,
                #     format=self.fmt,
                #     resampleAlg='cubicspline',
                #     options=["COMPRESS=LZW", "TILED=YES", "NBITS=1"]
                # )
                # g = None
                self.msk_fn = mask_fn
                
            ## Cleanup
            utils.remove_glob(*chunks)
            utils.remove_glob(*stack_chunks)
            utils.remove_glob(*mask_chunks)
            # handle aux cleanup if needed

        else:
            ## Standard Processing (Non-Chunked)
            mask_fn = None
            if self.want_stack:
                stack_name = f'{os.path.basename(self.name)}_stack'
                mask_name = f'{stack_name}_msk'
                mask_fn = (f'{os.path.join(self.cache_dir, mask_name)}'
                           f'.{gdalfun.gdal_fext(self.fmt)}')

                stack_fn = os.path.join(
                    self.cache_dir,
                    f'{stack_name}.{gdalfun.gdal_fext("GTiff")}'
                )
                stack_bn = utils.fn_basename2(stack_fn)
                
                if not self.clobber and os.path.exists(stack_fn):
                    self.stack = stack_fn
                    if not WaffleDEM(
                            self.stack,
                            cache_dir=self.cache_dir,
                            verbose=self.verbose
                    ).initialize().valid_p():
                        utils.echo_warning_msg(f'Existing stack file {stack_fn} is invalid, re-generating')
                        self.stack = self.data.stacks(out_name=stack_bn)
                else:
                    self.stack = self.data.stacks(out_name=stack_bn)
                    
                self.stack_ds = gdalfile.GDALFile(
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

                if self.want_mask:
                    mask_fn = f'{mask_name}.{gdalfun.gdal_fext(self.fmt)}'

                ## Apply limits to stack
                if self.count_limit is not None or self.uncertainty_limit is not None:
                    with gdalfun.gdal_datasource(self.stack, update=True) as stack_ds:
                        stack_infos = gdalfun.gdal_infos(stack_ds)
                        count_mask = None
                        uncertainty_mask = None

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
                            if self.count_limit and count_mask is not None:
                                this_arr[count_mask] = stack_infos['ndv']
                            if self.uncertainty_limit and uncertainty_mask is not None:
                                this_arr[uncertainty_mask] = stack_infos['ndv']
                            this_band.WriteArray(this_arr)
                            
                ## Run Waffles Module
                if WaffleDEM(
                        self.stack,
                        cache_dir=self.cache_dir,
                        verbose=self.verbose
                ).initialize().valid_p():
                    self.run()
                    
            else:
                self.run()    

            ## Post-Processing and validation.
            ## Calculate Uncertainty
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
            
            ## Post-Process DEM
            waffle_dem = WaffleDEM(
                self.fn,
                cache_dir=self.cache_dir,
                verbose=self.verbose,
                co=self.co
            ).initialize()
            
            if waffle_dem.valid_p():
                if mask_fn is not None and os.path.exists(mask_fn):
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
            else:
                utils.echo_warning_msg(f'{self.fn} is invalid')
                
            ## Post-Process Mask & Metadata
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
                            dst_fn=f'{os.path.basename(self.name)}_msk.{gdalfun.gdal_fext(self.fmt)}',
                            dst_dir=os.path.dirname(self.fn)
                        )
                        self.output_files['mask'] = mask_dem.fn
                        
                        if self.want_sm:
                            self.output_files['spatial-metadata'] = []
                            with gdalfun.gdal_datasource(mask_dem.fn) as msk_ds:
                                sm_layer, sm_fmt = spatial_metadata.polygonize_mask_multibands(
                                    msk_ds, verbose=True
                                )
                                sm_files = glob.glob(f'{sm_layer}.*')
                            
                            for f in sm_files:
                                out_sm = f'{self.name}_sm.{f.split(".")[-1]}'
                                if os.path.exists(out_sm):
                                    utils.remove_glob(out_sm)

                                os.rename(f, out_sm)
                                self.output_files['spatial-metadata'].append(out_sm)
                else:
                    utils.echo_warning_msg('Mask DEM is invalid...')
                    
            ## Post-Process Auxiliary Rasters
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

        if self.verbose:
            utils.echo_msg(f'\033[1mOutput files:\033[m {self.output_files}')
            
        return self    

    
    def run(self):
        """run the WAFFLES module (set via sub-module class)."""
        
        raise(NotImplementedError)

    
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
            f'{self.name}.{gdalfun.gdal_fext(self.fmt)}',
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
            utils.echo_msg(f'Stacked data to {self.fn}')            
        return self

    
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

    ## TODO: add parameters for specifying outputs...
    def __init__(self, min_count=None, **kwargs):
        super().__init__(**kwargs)
        self.min_count = min_count

        
    def run(self):
        return self


class WaffleFactory(factory.CUDEMFactory):
    """Waffles Factory Settings"""

    from . import flatten
    from . import idw
    from . import griddata
    from . import gmtgrid
    from . import mbs
    from . import gdalgrid
    from . import vdatum
    from . import coastline
    from . import lakes
    from . import uncertainty
    from . import cudemgrid
    from . import kriging
    from . import ml_interp
    from . import natural_neighbor
    from . import inpaint
    from . import cube
    from . import sdb
    
    _modules = {
        'stacks': {'name': 'stacks', 'stack': True, 'call': WafflesStacks},
        'IDW': {'name': 'IDW', 'stack': True, 'call': idw.WafflesIDW},
        'linear': {'name': 'linear', 'stack': True, 'call': griddata.WafflesLinear},
        'cubic': {'name': 'cubic', 'stack': True, 'call': griddata.WafflesCubic},
        'nearest': {'name': 'nearest', 'stack': True, 'call': griddata.WafflesNearest},
        'gmt-surface': {'name': 'surface', 'stack': True, 'call': gmtgrid.GMTSurface},
        'gmt-triangulate': {'name': 'triangulate','stack': True, 'call': gmtgrid.GMTTriangulate},
        'gmt-nearneighbor': {'name': 'nearneihbor', 'stack': True, 'call': gmtgrid.GMTNearNeighbor},
        'mbgrid': {'name': 'mbgrid', 'stack': True, 'call': mbs.WafflesMBGrid},
        'gdal-linear': {'name': 'linear', 'stack': True, 'call': gdalgrid.GDALLinear},
        'gdal-nearest': {'name': 'nearest', 'stack': True, 'call': gdalgrid.GDALNearest},
        'gdal-average': {'name': 'average', 'stack': True, 'call': gdalgrid.GDALMovingAverage},
        'gdal-invdst': {'name': 'invdst', 'stack': True, 'call': gdalgrid.GDALInvDst},
        'vdatum': {'name': 'vdatum', 'stack': False, 'call': vdatum.WafflesVDatum},
        'coastline': {'name': 'coastline', 'stack': False, 'call': coastline.WafflesCoastline},
        'lakes': {'name': 'lakes', 'stack': False, 'call': lakes.WafflesLakes},
        'cudem': {'name': 'cudem', 'stack': True, 'call': cudemgrid.WafflesCUDEM},
        'uncertainty': {'name': 'uncertainty', 'stack': True, 'call': uncertainty.WafflesSplitSampleUncertainty},
        'uncertainty-ml': {'name': 'uncertainty-ml', 'stack': True, 'call': uncertainty.WafflesMLUncertainty},
        'scratch': {'name': 'scratch', 'stack': True, 'call': WafflesScratch},
        'flatten': {'name': 'flatten', 'stack': True, 'call': flatten.WafflesFlatten},
        ## testing
        'kriging': {'name': 'kriging', 'stack': True, 'call': kriging.WafflesKriging},
        'ml-interp': {'name': 'ml_interp', 'stack': True, 'call': ml_interp.WafflesML},
        'ml-rf': {'name': 'ml-rf', 'stack': True, 'call': ml_interp.WafflesRF},
        'ml-knn': {'name': 'ml-knn', 'stack': True, 'call': ml_interp.WafflesKNN},
        'ml-mlp': {'name': 'ml-mlp', 'stack': True, 'call': ml_interp.WafflesMLP},
        'natural-neighbor': {'name': 'natural-neighbor', 'stack': True, 'call': natural_neighbor.WafflesNaturalNeighbor},
        'inpaint': {'name': 'inpaint', 'stack': True, 'call': inpaint.WafflesInpaint},
        'cube': {'name': 'cube', 'stack': True, 'call': cube.WafflesCUBE},
        'sdb': {'name': 'sdb', 'stack': True, 'call': sdb.WafflesSDB},
        #'num': {'name': 'num', 'stack': True, 'call': WafflesNum}, # defunct
        #'patch': {'name': 'patch', 'stack': True, 'call': WafflesPatch}, # test
        #'bgrid': {'name': 'bgrid', 'stack': True, 'call': WafflesBGrid}, # test            
    }

    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        
    def _set_modules(self):
        pass

    
def waffle_queue(q):
    """Waffle queue for threads from CLI

    waffles will run multiple input regions as
    separate threads...

    waffles_args is [waffle_module]
    """

    signal.signal(signal.SIGINT, signal.SIG_IGN)
    
    while True:
        waffle_module = q.get()
        if waffle_module is None:
            break
        else:
            try:
                waffle_module[0]()
            except Exception as e:
                utils.echo_error_msg(
                    f'Failed to generate {waffle_module}, {e}'
                )
                print(traceback.format_exc())
                pass

            
## ==============================================
## Command-line Interface (CLI)
## $ waffles
##
## waffles cli
## ==============================================
class PrintModulesAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        factory.echo_modules(WaffleFactory._modules, values, md=True if not values else False)
        sys.exit(0)

        
class PrintGritsModulesAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        factory.echo_modules(grits.GritsFactory._modules, values, md=True if not values else False)
        sys.exit(0)

        
class StackModeChoices(list):
    def __contains__(self, item):
        matches = [choice for choice in self if item.split(':')[0] in choice]
        return len(matches) == 1  # Only allow if it's a unique match

    
def waffles_cli():
    """Run waffles from command-line using argparse."""

    parser = argparse.ArgumentParser(
        description=f"Waffles ({__version__}): Generate DEMs and derivatives.",
        epilog=f"Supported Modules:\n{factory.get_module_short_desc(WaffleFactory._modules)}\n\n"\
        "CUDEM home page: <http://cudem.colorado.edu>",
        formatter_class=argparse.RawTextHelpFormatter
    )

    ## --- Data Input ---
    parser.add_argument(
        'datalist', 
        nargs='*', 
        help="Input datalist(s) or data files."
    )
    
    ## --- Core Parameters ---
    ## TODO: -R fails with space in cli: -R -65.0/-63.0/44.85/45.85
    parser.add_argument(
        '-R', '--region', '--aoi', 
        action='append',
        help=regions.region_help_msg()
    )
    parser.add_argument(
        '-E', '--increment',
        required=True,
        help=("Gridding Increment (xinc[/yinc][:xsample/ysample]).\n"
              "Note: Increment units should be the same as the input region.")
    )
    parser.add_argument(
        '-M', '--module', 
        default='stacks', 
        help=("Waffles Module and options (e.g. 'gmt-surface:tension=0.35'). Default: stacks.\n"
              "Where MODULE is module[:mod_opt=mod_val[:mod_opt1=mod_val1[:...]]]")
    )
    
    ## --- Output Control ---
    parser.add_argument(
        '-O', '--output-basename', 
        default='waffles_dem', 
        help="Output Basename. Default: waffles_dem."
    )
    parser.add_argument(
        '-F', '--format', 
        default='GTiff', 
        help="Output Format (GTiff, NetCDF, etc.). Default: GTiff."
    )
    parser.add_argument(
        '-P', '--t_srs', 
        default='epsg:4326', 
        help=("Target Projection (EPSG code or WKT). Default: epsg:4326.\n"
              "Note: The target projection should align with the input region.")
    )
    parser.add_argument(
        '-f', '--transform', 
        action='store_true', 
        help="Transform inputs to target projection."
    )
    parser.add_argument(
        '-N', '--nodata', 
        type=float, 
        default=-9999, 
        help="NoData Value. Default: -9999."
    )
    parser.add_argument(
        '-CO', '--creation-options', 
        action='append', 
        help="GDAL creation options (e.g. COMPRESS=DEFLATE)."
    )
    parser.add_argument(
        '-p', '--prefix',
        action='store_true',
        help="Append prefix info to output name (res=X:year=XXXX:version=X)."
    )

    ## --- Processing Options ---
    parser.add_argument(
        '-X', '--extend', 
        help=("Extend region (cells[:percent]). e.g. '6:10'.\n"
              "Where 'cells' extends the output region and 'percent' extends the processing region.")
    )
    parser.add_argument(
        '-T', '--filter', 
        action='append', 
        help=("Apply Grits Filter (e.g. 'outliers:k=2.5').\n"
              "Where FILTER is fltr_name[:opts] (see `grits --modules` for more information)\n"
              "This may be set multiple times to perform multiple filters.\n"
              "Append `:stacks=True` to the filter to perform the filter on the data stack")
    )
    parser.add_argument(
        '-C', '--clip', 
        help="Clip output to polygon file. e.g. 'clip_ply.shp:invert=False'"
    )
    parser.add_argument(
        '-K', '--chunk', 
        type=int, 
        help="Process in chunks of N x N pixels."
    )
    parser.add_argument(
        '-S', '--sample-alg', 
        default='bilinear', 
        help="Resampling algorithm (bilinear, cubic, nearest, etc.)."
    )
    
    ## --- Stack Control ---
    parser.add_argument(
        '-w', '--want-weight', 
        action='store_true', 
        help="Use weights from datalist."
    )
    parser.add_argument(
        '-u', '--want-uncertainty', 
        action='store_true', 
        help="Generate/Use uncertainty."
    )
    parser.add_argument(
        '-m', '--want-mask', 
        action='store_true', 
        help="Generate data mask."
    )
    parser.add_argument(
        '-A', '--stack-mode', 
        default='mean', 
        choices=StackModeChoices(['mean', 'min', 'max', 'mixed', 'supercede']),
        help="Stacking mode."
    )
    parser.add_argument(
        '-L', '--limits', 
        action='append', 
        help="Set limits (u=upper, l=lower, p=prox, s=size, n=count, r=uncertainty, c=nodata, e=expand). e.g. -Lu0."
    )

    ## --- System / Misc ---
    parser.add_argument(
        '-H', '--threads', 
        type=int, 
        default=1, 
        help="Number of threads for processing regions."
    )
    parser.add_argument(
        '-D', '--cache-dir', 
        help="Cache directory location."
    )
    parser.add_argument(
        '-k', '--keep-cache', 
        action='store_true', 
        help="Keep cache files after run."
    )
    parser.add_argument(
        '-x', '--keep-auxiliary', 
        action='store_true', 
        help="Keep auxiliary rasters (stack, mask, etc.)."
    )
    parser.add_argument(
        '-s', '--spatial-metadata', 
        action='store_true', 
        help="Generate spatial metadata polygons."
    )
    parser.add_argument(
        '-a', '--archive', 
        action='store_true', 
        help="Archive the processed datalist."
    )
    parser.add_argument(
        '-c', '--continue', 
        dest='clobber', 
        action='store_false', 
        help="Don't clobber existing files."
    )
    parser.add_argument(
        '-q', '--quiet', 
        action='store_true', 
        help="Quiet mode."
    )
    parser.add_argument(
        '--config', 
        help="Load from or save to config JSON."
    )
    parser.add_argument(
        '--modules', 
        nargs='?',
        action=PrintModulesAction, 
        help="List available modules."
    )
    parser.add_argument(
        '--filters', 
        nargs='?',
        action=PrintGritsModulesAction, 
        help="List available filter modules. (via grits)"
    )
    parser.add_argument(
        '--version', 
        action='version', 
        version=f'CUDEM {__cudem_version__} :: %(prog)s {__version__}'
    )

    ## Parse
    args = parser.parse_args()

    ## Handle --modules early exit
    if args.modules:
        factory.echo_modules(WaffleFactory._modules)
        sys.exit(0)

    ## Load the user wg json and run waffles with that.
    ## TODO: allow input of multiple config files with -G .. -G ..
    if args.config:
        if os.path.exists(args.config):
            with open(args.config, 'r') as wg_json:
                wg = json.load(wg_json)
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
            if not wg['chunk']:
                waffle_q.put([this_waffle_module])
            else:
                try:
                    this_waffle_module()
                except Exception as e:
                    utils.echo_error_msg(
                        f'Failed to generate {this_waffle_module}, {e}'
                    )
                    print(traceback.format_exc())
                
            sys.exit(0)   
        else:
            utils.echo_error_msg(
                f'specified waffles config file does not exist, {args.config}'
            )
        
    ## Initialize Config Dictionary
    wg = {
        'name': args.output_basename,
        'fmt': args.format,
        'dst_srs': args.t_srs,
        'srs_transform': args.transform,
        'ndv': args.nodata,
        'sample': args.sample_alg,
        'chunk': args.chunk,
        'clip': args.clip,
        'want_weight': args.want_weight,
        'want_uncertainty': args.want_uncertainty,
        'want_mask': args.want_mask or args.spatial_metadata, # SM implies mask
        'want_sm': args.spatial_metadata,
        'archive': args.archive,
        'keep_auxiliary': args.keep_auxiliary,
        'clobber': args.clobber,
        'verbose': not args.quiet,
        'stack_mode': args.stack_mode,
        'fltr': args.filter if args.filter else [],
        'co': args.creation_options if args.creation_options else [],
        'data': args.datalist
    }

    ## --- Post-Process Special Arguments ---
    
    ## Increments (-E)
    if args.increment:
        incs = args.increment.split(':')
        xy_inc = incs[0].split('/')
        wg['xinc'] = utils.str2inc(xy_inc[0])
        wg['yinc'] = utils.str2inc(xy_inc[1]) if len(xy_inc) > 1 else wg['xinc']
        
        if len(incs) > 1:
            xy_samples = incs[1].split('/')
            wg['xsample'] = utils.str2inc(xy_samples[0])
            wg['ysample'] = utils.str2inc(xy_samples[1]) if len(xy_samples) > 1 else wg['xsample']
    else:
        ## Check if config file provided, otherwise error
        if not args.config:
            utils.echo_error_msg("Must specify gridding increment (-E).")
            sys.exit(-1)

    ## Extend (-X)
    if args.extend:
        exts = args.extend.split(':')
        wg['extend'] = utils.int_or(exts[0], 0)
        if len(exts) > 1: wg['extend_proc'] = utils.float_or(exts[1], 10)

    if WaffleFactory._modules[args.module.split(':')[0]]['stack']:
        if len(args.datalist) == 0:
            sys.stderr.write(waffles_cli_usage())
            utils.echo_error_msg(
                '''You must specify a datalist/entry, try `gmrt` or `srtm` for global data.'''
            )
            sys.exit(-1)
    else:
        wg['want_stack'] = True if len(args.datalist) > 0 else False
        
    ## Limits (-L)
    ## Map CLI flags to config keys
    limit_map = {
        'u': 'upper_limit', 'l': 'lower_limit', 
        'p': 'proximity_limit', 's': 'size_limit',
        'c': 'percentile_limit', 'e': 'expand_limit',
        'n': 'count_limit', 'r': 'uncertainty_limit'
    }
    if args.limits:
        for limit_str in args.limits:
            key_char = limit_str[0]
            val = limit_str[1:]
            if key_char in limit_map:
                # Convert to int or float based on expected type (simplified here to float)
                wg[limit_map[key_char]] = utils.float_or(val)

    ## Cache Dir
    if args.cache_dir:
        wg['cache_dir'] = args.cache_dir
    else:
        wg['cache_dir'] = waffles_cache#os.path.join(os.path.expanduser('~'), 'cudem_cache')        
        
    ## --- Execution ---
    
    ## Setup Regions
    ## (Using original parsing logic for consistency)
    i_regions = args.region if args.region else [None]
    these_regions = regions.parse_cli_region(i_regions, wg['verbose'])

    ## Setup Multiprocessing regions
    if not args.chunk:
        waffle_q = mp.Queue()
        processes = []

    try:
        ## For mp regions
        if not args.chunk:
            for _ in range(args.threads):
                t = mp.Process(target=waffle_queue, args=([waffle_q]))
                t.daemon = True
                processes.append(t)
                t.start()

        ## Iterate Regions
        for this_region in these_regions:
            # Auto-detect region from data if missing
            if this_region is None:
                utils.echo_warning_msg('No input region specified, gathering region from input data...')
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

                utils.echo_msg(f'Region is {this_region}')
                if this_region is None: # couldn't gather a region from the data
                    break

            wg['src_region'] = this_region

            ## Prefix (-p)
            if args.prefix or (args.region and len(args.region) > 1):
                prefix_args = {}
                x_sample = wg.get('xsample')
                x_inc = wg.get('xinc')

                wg['name'] = utils.append_fn(
                    args.output_basename,
                    this_region,
                    x_sample if x_sample is not None else x_inc,
                    **prefix_args
                )

            ## Module Init
            mod_name = args.module.split(':')[0]
            if mod_name not in WaffleFactory._modules:
                utils.echo_error_msg(f"Invalid module: {mod_name}")
                continue

            ## Export the waffles module config file
            if args.config: 
                wg['src_region'] = this_region.export_as_list()
                wg['name'] = os.path.abspath(wg['name'])
                ## Create Factory
                ## We pass 'mod' string to factory which parses opts
                utils.echo_msg(wg)
                this_waffle = WaffleFactory(mod=args.module, **wg)
                this_waffle.write_parameter_file(f'{args.config}.json')
                continue

            ## Create Factory
            ## We pass 'mod' string to factory which parses opts
            this_waffle = WaffleFactory(mod=args.module, **wg)
            
            ## Acquire and Queue
            module_instance = this_waffle._acquire_module()
            if module_instance:
                ## for mp regions
                if not args.chunk:
                    waffle_q.put([module_instance])
                else:
                    try:
                        module_instance()
                    except Exception as e:
                        utils.echo_error_msg(
                            f'Failed to generate {module_instance}, {e}'
                        )
                        print(traceback.format_exc())

        ## Cleanup for mp regions
        if not args.chunk:
            for _ in range(args.threads):
                waffle_q.put(None)
            
            for t in processes:
                t.join()

    except KeyboardInterrupt:
        utils.echo_error_msg("Killed by user, Terminating processes...\n")
        ## for mp regions
        if not args.chunk:
            for t in processes:
                if t.is_alive():
                    t.terminate() 
                    t.join()
        sys.exit(1)
        
    if not args.keep_cache:
        utils.remove_glob(wg['cache_dir'])

    return 0

### End
