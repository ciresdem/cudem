### dlim.py - DataLists IMproved
##
## Copyright (c) 2010 - 2026 Regents of the University of Colorado
##
## dlim.py is part of CUDEM
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
## Dlim: The CUDEM Data Ingestion and Processing Engine.
##
## Dlim (Data Lists IMproved) manages the ingestion, processing, and standardization
## of diverse elevation datasets. It utilizes text-based "datalists" to organize
## hierarchical data collections and stream them into unified processing pipelines.
##
##   * Robust Format Handling:
##      - Auto-detects and parses Raster (GeoTIFF, BAG, NetCDF), Point Cloud
##        (XYZ, LAS/LAZ, IceSat2), and Vector (OGR, Shapefile) formats.
##      - Abstracts file access, allowing downstream tools to treat all inputs
##        as a uniform stream of XYZ points or Raster arrays.
##
##   * On-the-Fly Transformation:
##      - Handles vertical and horizontal coordinate transformations using Proj
##        and VDatum grids during data streaming.
##      - Supports spatial windowing (regions) and resolution blocking.
##
##   * Data Management & Archival:
##      - Archives processed datasets into standardized directory structures or
##        tarballs for long-term storage.
##      - Generates spatial metadata (footprints) and JSON/INF catalogs.
##
## Usage:
##   CLI: dlim -R <region> -E <increment> [options] input.datalist
##   API: dl = dlim.DatasetFactory(fn="input.xyz", src_srs="epsg:4326").acquire()
##        for point in dl.yield_xyz(): ...
##
### Code:

import os
import sys
import re
import copy
import json
import glob
import math
import traceback
import warnings
import shlex
import argparse
from typing import Optional, List, Dict, Any, Union

import numpy as np
import pyproj
import h5py as h5
import netCDF4 as nc
from osgeo import gdal, ogr

from cudem import utils
from cudem import regions
from cudem import srsfun
from cudem import gdalfun
from cudem import factory
from cudem import vdatums
from cudem import xyzfun
from cudem import fetches
from cudem import pointz
from cudem.grits import grits
from cudem.globato import globato
from cudem.globato import globato_converter
from cudem.datalists import inf
from cudem.fetches import fetches
from cudem import __version__ as __cudem_version__
from . import __version__
from . import inf

## Config info and setup
gc = utils.config_check()
gdal.UseExceptions()
ogr.UseExceptions()
gdal.SetConfigOption(
    'CPL_LOG', 'NUL' if gc['platform'] == 'win32' else '/dev/null'
) 

## ==============================================
## Datalist Helper Functions
## data_list is a list of dlim supported datasets
## ==============================================
def get_factory_exts():
    fmts = []
    for key in DatasetFactory()._modules.keys():
        if key != '_factory' and int(key) > 0:
            for f in DatasetFactory()._modules[key]['fmts']:
                if f not in fmts:
                    fmts.append(f)                    
    return fmts


def make_datalist(data_list, want_weight, want_uncertainty, region,
                  src_srs, dst_srs, x_inc, y_inc, sample_alg,
                  verbose):
    """Make a datalist object from a list of supported datasets"""

    ## Make the datalist object
    xdl = Datalist(
        fn='scratch',
        data_format=-1,
        weight=None if not want_weight else 1,
        uncertainty=None if not want_uncertainty else 0,
        src_region=region,
        verbose=verbose,
        parent=None,
        src_srs=src_srs,
        dst_srs=dst_srs,
        x_inc=x_inc,
        y_inc=y_inc,
        sample_alg=sample_alg
    )

    ## add the datasets from data_list to the datalist object
    xdl.data_entries = [DatasetFactory(
        fn=" ".join(['-' if x == "" else x for x in dl.split(",")]),
        weight=None if not want_weight else 1,
        uncertainty=None if not want_uncertainty else 0,
        src_region=region,
        verbose=verbose,
        src_srs=src_srs,
        dst_srs=dst_srs,
        x_inc=x_inc,
        y_inc=y_inc,
        sample_alg=sample_alg,
        parent=xdl,
        cache_dir=self.cache_dir
    )._acquire_module() for dl in data_list]
    return xdl


def write_datalist(data_list, outname=None):
    """write a datalist file from a datalist object"""
    
    if outname is None:
        outname = '{}_{}'.format(self.metadata['name'], utils.this_year())
    
    if os.path.exists('{}.datalist'.format(outname)):
        utils.remove_glob('{}.datalist*'.format(outname))
        
    with open('{}.datalist'.format(outname), 'w') as tmp_dl:
        [tmp_dl.write('{}\n'.format(x.format_entry())) for x in data_list]

    return '{}.datalist'.format(outname)


def h5_get_datasets_from_grp(grp, count = 0, out_keys = []):
    if isinstance(grp, h5.Dataset):
        count += 1
        out_keys.append(grp.name)

    elif isinstance(grp, h5.Group):
        for key in grp.keys():
            if isinstance(grp[key], h5.Dataset):
                count += 1
                out_keys.append(grp[key].name)
            elif isinstance(grp[key], h5.Group):
                for sub_key in grp[key].keys():
                    count, out_keys = h5_get_datasets_from_grp(
                        grp[key][sub_key], count, out_keys
                    )

    return count, out_keys


def init_data(data_list, **kwargs):
    """Initialize a datalist object from a list of supported dataset entries.
    
    Args:
        data_list (list): A list containing file paths (str), dicts, or sys.stdin.
        **kwargs: Additional arguments passed to the dataset constructors.
    
    Returns:
        Datalist object or None.
    """

    from . import datalists
    from . import xyzfile
    
    xdls = []
    
    ## Map convenience kwargs to standard internal keys
    kwargs['weight'] = kwargs.get('want_weight', kwargs.get('weight'))
    kwargs['uncertainty'] = kwargs.get('want_uncertainty', kwargs.get('uncertainty'))
    kwargs['verbose'] = kwargs.get('want_verbose', kwargs.get('verbose'))

    try:
        for dl in data_list:
            ## Handle Standard Input (Pipe)
            if dl is sys.stdin:
                xdls.append(xyzfile.XYZFile(fn=dl, data_format=168, **kwargs))

            ## Handle String Entries (e.g. "path/to/data.xyz,168,1")
            elif isinstance(dl, str):
                ## Convert CSV-style input to Space-delimited string for the Factory parser.
                ## Replaces empty CSV fields (,,) with hyphens (-) which Factory recognizes as "inherit/default".
                ## Example: "data.tif,,0.5" -> "data.tif - 0.5"
                entry_str = " ".join(['-' if x.strip() == "" else x for x in dl.split(",")])
                
                ## Pass the formatted string as 'mod' to the factory
                ds = DatasetFactory(mod=entry_str, **kwargs)._acquire_module()
                if ds:
                    xdls.append(ds)

            ## Handle Dictionary Entries
            elif isinstance(dl, dict):
                ## Ensure source region is a Region object if provided
                if dl.get('kwargs', {}).get('src_region') is not None:
                    dl['kwargs']['src_region'] = regions.Region().from_list(
                        dl['kwargs']['src_region']
                    )
                
                ## Merge global kwargs with specific entry kwargs
                ## Entry kwargs take precedence
                combined_kwargs = kwargs.copy()
                combined_kwargs.update(dl.get('kwargs', {}))

                ds = DatasetFactory(
                    mod=dl.get('mod_name'), 
                    **combined_kwargs
                )._acquire_module()
                
                if ds:
                    xdls.append(ds)

        ## Wrap Results
        if len(xdls) > 1:
            ## If multiple datasets, wrap them in a Scratch datalist (in-memory list)
            this_datalist = datalists.Scratch(fn=xdls, data_format=-3, **kwargs)
        elif len(xdls) > 0:
            ## If single dataset, return it directly
            this_datalist = xdls[0]
        else:
            this_datalist = None

        return this_datalist

    except Exception as e:
        utils.echo_error_msg(f'Could not initialize data from {data_list}: {e}')
        return None

    
## ==============================================
## ElevationDataset
## Base Elevation Dataset Class
## ==============================================
class ElevationDataset:
    """Representing an Elevation Dataset.
    
    This is the super class for all datalist (dlim) datasets.
    Each dataset sub-class should define a dataset-specific
    data parser <self.parse> and a <self.generate_inf> function to generate
    inf files.

    Minimally required methods in subclasses:
        - __init__
        - yield_points(self)
    """

    gdal_sample_methods = [
        'near', 'bilinear', 'cubic', 'cubicspline', 'lanczos',
        'average', 'mode', 'max', 'min', 'med', 'Q1', 'Q3', 'sum',
        'auto'
    ]

    stack_modes = [
        'min', 'max', 'mean', 'supercede', 'mixed', 'std', 'var', 'weights'
    ]

    def __init__(self,
                 fn=None,
                 data_format=None,
                 weight=1,
                 uncertainty=0,
                 src_srs=None,
                 mask=None,
                 dst_srs='epsg:4326',
                 src_geoid=None,
                 dst_geoid='g2018',
                 x_inc=None,
                 y_inc=None,
                 want_mask=False,
                 want_sm=False,
                 sample_alg='auto',
                 parent=None,
                 src_region=None,
                 invert_region=False,
                 pnt_fltrs=None,
                 stack_fltrs=None,
                 stack_node=True,
                 stack_mode='mean',
                 cache_dir=None,
                 verbose=False,
                 remote=False,
                 dump_precision=6,
                 upper_limit=None,
                 lower_limit=None,
                 params=None,
                 metadata=None,
                 **kwargs):
        
        self.fn = fn
        self.data_format = data_format
        self.weight = weight
        self.uncertainty = uncertainty
        self.mask = mask
        self.src_srs = src_srs
        self.dst_srs = dst_srs
        self.src_geoid = src_geoid
        self.dst_geoid = dst_geoid
        self.x_inc = utils.str2inc(x_inc)
        self.y_inc = utils.str2inc(y_inc)
        self.stack_fltrs = stack_fltrs if stack_fltrs else []
        self.pnt_fltrs = pnt_fltrs if pnt_fltrs else []
        self.want_mask = want_mask
        self.want_sm = want_sm
        
        if self.want_sm:
            self.want_mask = True
            
        self.sample_alg = sample_alg
        self.metadata = copy.deepcopy(metadata) if metadata else {
            'name': None, 'title': None, 'source': None, 'date': None,
            'data_type': None, 'resolution': None, 'hdatum': None,
            'vdatum': None, 'url': None
        }
        self.parent = parent
        
        self.region = src_region
        if self.region is not None:
            self.region = regions.Region().from_user_input(self.region)
            
        self.invert_region = invert_region
        self.cache_dir = cache_dir
        self.verbose = verbose
        self.remote = remote
        self.dump_precision = dump_precision
        self.stack_node = stack_node
        self.stack_mode = stack_mode
        self._init_stack_mode()

        self.upper_limit = utils.float_or(upper_limit)
        self.lower_limit = utils.float_or(lower_limit)

        self.infos = inf.INF(
            name=self.fn, file_hash='0', numpts=0, fmt=self.data_format
        )

        self.params = params
        if not self.params:
            self.params = {}
            self.params['kwargs'] = self.__dict__.copy()
            self.params['mod'] = self.fn
            self.params['mod_name'] = self.data_format
            self.params['mod_args'] = {}

            
    def __str__(self):
        return f'<Dataset: {self.metadata.get("name")} - {self.fn}>'

    
    def __repr__(self):
        return f'<Dataset: {self.metadata.get("name")} - {self.fn}>'

    
    def __call__(self):
        self.initialize()

        
    def _init_stack_mode(self):
        opts, self.stack_mode_name, self.stack_mode_args = factory.parse_fmod(self.stack_mode)

        if self.stack_mode_name not in self.stack_modes:
            utils.echo_warning_msg(f'Mode: {self.stack_mode_name} is not a valid stack mode')
            self.stack_mode_name = 'mean'
        
        if 'mask_level' not in self.stack_mode_args:
            self.stack_mode_args['mask_level'] = -1

            
    def _init_mask(self, mask):
        opts = factory.fmod2dict(mask, {})
        if 'mask_fn' not in opts:
            if '_module' in opts:
                opts['mask_fn'] = opts['_module']
            else:
                utils.echo_error_msg(f'Could not parse mask {self.mask}')
                return None

        opts.setdefault('invert', False)
        opts.setdefault('verbose', False)
        opts.setdefault('min_z', None)
        opts.setdefault('max_z', None)
        # if 'invert' not in opts:
        #     opts['invert'] = False
        # if 'verbose' not in opts:
        #     opts['verbose'] = False
        # if 'min_z' not in opts:
        #     opts['min_z'] = None
        # if 'max_z' not in opts:
        #     opts['max_z'] = None

        if opts['min_z'] is not None: opts['min_z'] = utils.float_or(opts['min_z'])
        if opts['max_z'] is not None: opts['max_z'] = utils.float_or(opts['max_z'])
            
        # Mask is OGR, rasterize it
        opts['ogr_or_gdal'] = gdalfun.ogr_or_gdal(opts['mask_fn'])
        
        data_mask = opts['mask_fn']
        if opts['ogr_or_gdal'] == 1: 
            if self.region is not None and self.x_inc is not None and self.y_inc is not None:
                data_mask = gdalfun.ogr2gdal_mask(
                    opts['mask_fn'],
                    region=self.region,
                    x_inc=self.x_inc,
                    y_inc=self.y_inc,
                    dst_srs=self.dst_srs,
                    invert=False,
                    verbose=True,
                    temp_dir=self.cache_dir,
                )
                opts['ogr_or_gdal'] = 0 # Now it's a raster
        
        opts['data_mask'] = data_mask
        return opts

    
    def _set_params(self, **kwargs):
        metadata = copy.deepcopy(self.metadata)
        _params = {
            'parent': self,
            'weight': self.weight,
            'uncertainty': self.uncertainty,
            'src_region': self.region,
            'invert_region': self.invert_region,
            'x_inc': self.x_inc,
            'y_inc': self.y_inc,
            'mask': self.mask,            
            'src_srs': self.src_srs,
            'dst_srs': self.dst_srs,
            'stack_fltrs': self.stack_fltrs,
            'pnt_fltrs': self.pnt_fltrs,
            'cache_dir': self.cache_dir,
            'upper_limit': self.upper_limit,
            'lower_limit': self.lower_limit,
            'verbose': self.verbose,
            'metadata': metadata,
            'stack_mode': self.stack_mode,
        }
        for kw in kwargs:
            _params[kw] = kwargs[kw]

        return _params

    
    def _copy_params(self, **kwargs):
        _params = {i: self.params['kwargs'][i] for i in self.params['kwargs'] if i != 'params'}
        for kw in kwargs:
            _params[kw] = kwargs[kw]

        return _params

    
    def _sub_init(self):
        pass

    
    ## initialize transfom
    def parse_srs(self):
        if self.src_srs is not None and self.dst_srs is not None:
            want_vertical = True
            src_geoid = None
            dst_geoid = 'g2018' # dst_geoid is g2018 by default.
            in_vertical_epsg = in_vertical_crs = None
            out_vertical_epsg = out_vertical_crs = None
            
            ## parse out the source geoid, which is not standard in proj,
            ## reset `self.src_srs` without it.
            tmp_src_srs = self.src_srs.split('+geoid:')
            src_srs = tmp_src_srs[0]
            self.src_srs = src_srs
            if len(tmp_src_srs) > 1:
                src_geoid = tmp_src_srs[1]

            ## parse out the destination geoid, which is not standard in proj
            ## reset `self.dst_srs` without it.
            tmp_dst_srs = self.dst_srs.split('+geoid:')
            dst_srs = tmp_dst_srs[0]
            self.dst_srs = dst_srs
            if len(tmp_dst_srs) > 1:
                dst_geoid = tmp_dst_srs[1]

            ## check if this is an ESRI epsg code
            is_esri = False
            in_vertical_epsg_esri = None
            if 'ESRI' in src_srs.upper():
                is_esri = True
                srs_split = src_srs.split('+')
                src_srs = srs_split[0]
                if len(srs_split) > 1:
                    in_vertical_epsg_esri = srs_split[1]

            ## check if the vertical epsg is tidal, fix the epsg code if so
            if utils.int_or(src_srs.split('+')[-1]) in vdatums._tidal_frames.keys():
                src_srs = '{}+{}'.format(
                    src_srs.split('+')[0],
                    vdatums._tidal_frames[utils.int_or(src_srs.split('+')[-1])]['epsg']
                )
                
            ## set the proj crs from the src and dst srs input
            try:
                in_crs = pyproj.CRS.from_user_input(src_srs)
                out_crs = pyproj.CRS.from_user_input(dst_srs)
            except:
                utils.echo_error_msg([src_srs, dst_srs])

            ## if the crs has vertical (compound), parse out the vertical crs
            ## and set the horizontal and vertical crs
            if in_crs.is_compound:
                in_crs_list = in_crs.sub_crs_list
                in_horizontal_crs = in_crs_list[0]
                in_vertical_crs = in_crs_list[1]
                in_vertical_name = in_vertical_crs.name
                in_vertical_epsg = in_vertical_crs.to_epsg()
                if in_vertical_epsg is None:
                    in_vertical_epsg = in_vertical_name.split(' ')[0]
            else:
                in_horizontal_crs = in_crs
                want_vertical=False

            if out_crs.is_compound:            
                out_crs_list = out_crs.sub_crs_list
                out_horizontal_crs = out_crs_list[0]
                out_vertical_crs = out_crs_list[1]
                out_vertical_epsg = out_vertical_crs.to_epsg()
            else:
                out_horizontal_crs = out_crs
                want_vertical=False
                out_vertical_epsg=None

            ## check if esri vertical
            if (in_vertical_epsg_esri is not None and is_esri):
                in_vertical_epsg = in_vertical_epsg_esri
                if out_vertical_epsg is not None:
                    want_vertical = True

            ## make sure the input and output vertical epsg is different
            if want_vertical:
                if (in_vertical_epsg == out_vertical_epsg) and self.src_geoid is None:
                    want_vertical = False

            ## set `self.transform` with the parsed srs info
            self.transform['src_horz_crs'] = in_horizontal_crs
            self.transform['dst_horz_crs'] = out_horizontal_crs
            self.transform['src_vert_crs'] = in_vertical_crs
            self.transform['dst_vert_crs'] = out_vertical_crs
            self.transform['src_vert_epsg'] = in_vertical_epsg
            self.transform['dst_vert_epsg'] = out_vertical_epsg
            self.transform['src_geoid'] = src_geoid
            self.transform['dst_geoid'] = dst_geoid
            self.transform['want_vertical'] = want_vertical        

            
    def set_vertical_transform(self):
        ## set the region for the vdatum transformation grid.
        ## this is either from the input `self.region` or from the input
        ## data's bounding box. Transform it to WGS84.
        if self.region is None:
            vd_region = regions.Region().from_list(self.infos.minmax)
            vd_region.src_srs = self.transform['src_horz_crs'].to_proj4()
        else:
            vd_region = self.region.copy()
            vd_region.src_srs = self.transform['dst_horz_crs'].to_proj4()

        vd_region.zmin = None
        vd_region.zmax = None
        vd_region.warp('epsg:4326')
        vd_region.buffer(pct=10)
        if not vd_region.valid_p():
            utils.echo_warning_msg('failed to generate transformation')
            return
        # else:
        #     utils.echo_msg('generating vertical transformation to region {}'.format(vd_region))

        ## set `self.transform.trans_fn`, which is the transformation grid
        self.transform['trans_fn'] = os.path.join(
            self.cache_dir, '_vdatum_trans_{}_{}_{}.tif'.format(
                self.transform['src_vert_epsg'],
                self.transform['dst_vert_epsg'],
                vd_region.format('fn')
            )
        )
        
        ## if the transformation grid already exists, skip making a new one,
        ## otherwise, make the new one here with `vdatums.VerticalTransform()`
        if not os.path.exists(self.transform['trans_fn']):
            with utils.ccp(
                    desc='generating vertical transformation grid {} from {} to {}'.format(
                        self.transform['trans_fn'], self.transform['src_vert_epsg'],
                        self.transform['dst_vert_epsg']
                    ),
                    leave=self.verbose
            ) as pbar:
                ## set the vertical transformation grid to be 3 arc-seconds. This
                ## is pretty arbitrary, maybe it's too small...
                vd_x_inc = vd_y_inc = utils.str2inc('3s')
                xcount, ycount, dst_gt = vd_region.geo_transform(
                    x_inc=vd_x_inc, y_inc=vd_y_inc, node='grid'
                )

                ## if the input region is so small it creates a tiny grid,
                ## keep increasing the increments until we are at least to
                ## a 10x10 grid.
                while (xcount <=10 or ycount <=10):
                    vd_x_inc /= 2
                    vd_y_inc /= 2
                    xcount, ycount, dst_gt = vd_region.geo_transform(
                        x_inc=vd_x_inc, y_inc=vd_y_inc, node='grid'
                    )

                ## run `vdatums.VerticalTransform()`, grid using `nearest`
                self.transform['trans_fn'], self.transform['trans_fn_unc'] \
                    = vdatums.VerticalTransform(
                        'IDW',
                        vd_region,
                        vd_x_inc,
                        vd_y_inc,
                        self.transform['src_vert_epsg'],
                        self.transform['dst_vert_epsg'],
                        geoid_in=self.transform['src_geoid'],
                        geoid_out=self.transform['dst_geoid'],
                        cache_dir=self.cache_dir,
                        verbose=False
                    ).run(outfile=self.transform['trans_fn'])

        ## set the pyproj.Transformer for both horz+vert and vert only
        ## hack for navd88 datums in feet (6360 is us-feet, 8228 is international-feet
        if utils.str_or(self.transform['src_vert_epsg']) == '6360':
            # or 'us-ft' in utils.str_or(src_vert, ''):
            #out_src_srs = out_src_srs + ' +vto_meter=0.3048006096012192'
            uc = ' +step +proj=unitconvert +z_in=us-ft +z_out=m'
        elif utils.str_or(self.transform['src_vert_epsg']) == '8228':
            uc = ' +step +proj=unitconvert +z_in=ft +z_out=m'
        else:
            uc = ''
            
        if self.transform['trans_fn'] is not None \
           and os.path.exists(self.transform['trans_fn']):
            self.transform['pipeline'] \
                = (f'+proj=pipeline{uc} +step '
                   f'{self.transform["src_horz_crs"].to_proj4()} '
                   '+inv +step +proj=vgridshift '
                   f'+grids="{os.path.abspath(self.transform["trans_fn"])}" '
                   f'+inv +step {self.transform["dst_horz_crs"].to_proj4()}')
            self.transform['vert_transformer'] = pyproj.Transformer.from_pipeline(
                (f'+proj=pipeline{uc} +step +proj=vgridshift '
                 f'+grids="{os.path.abspath(self.transform["trans_fn"])}" +inv')
            )
            utils.echo_debug_msg(self.transform['pipeline'])
        else:
            utils.echo_error_msg(
                ('failed to generate vertical transformation grid between '
                 f'{self.transform["src_vert_epsg"]} and '
                 f'{os.path.abspath(self.transform["dst_vert_epsg"])} '
                 'for this region!')
            )

            
    def set_transform(self):
        """Set the pyproj horizontal and vertical transformations 
        for the dataset
        """

        self.parse_srs()
        if self.transform['src_horz_crs'] is not None \
           and self.transform['dst_horz_crs'] is not None:        
            ## horizontal Transformation
            self.transform['horz_pipeline'] \
                = ('+proj=pipeline +step '
                   f'{self.transform["src_horz_crs"].to_proj4()} '
                   f'+inv +step {self.transform["dst_horz_crs"].to_proj4()}')
            if self.region is not None:
                self.transform['trans_region'] = self.region.copy()
                self.transform['trans_region'].src_srs \
                    = self.transform['dst_horz_crs'].to_proj4()
                self.transform['trans_region'].warp(
                    self.transform['src_horz_crs'].to_proj4()
                )
            else:
                self.transform['trans_region'] \
                    = regions.Region().from_list(self.infos.minmax)
                self.transform['trans_region'].src_srs \
                    = self.infos.src_srs
                self.transform['trans_region'].warp(
                    self.transform['dst_horz_crs'].to_proj4()
                )

            ## vertical Transformation
            if self.transform['want_vertical']:
                self.set_vertical_transform()
            else:
                self.transform['pipeline'] = self.transform['horz_pipeline']

            try:
                self.transform['transformer'] \
                    = pyproj.Transformer.from_pipeline(
                        self.transform['pipeline']
                    )
            except Exception as e:
                utils.echo_warning_msg(
                    ('could not set transformation in: '
                     f'{self.transform["src_horz_crs"].name}, out: '
                     f'{self.transform["dst_horz_crs"].name}, {e}')
                    )

                return

        ## dataset region
        ## mrl: moved out of if block 
        if self.region is not None and self.region.valid_p():
            self.data_region = self.region.copy() \
                if self.transform['trans_region'] is None \
                   else self.transform['trans_region'].copy()
            #inf_region = regions.Region().from_list(self.infos.minmax)
            self.data_region = regions.regions_reduce(
                self.data_region, self.inf_region
            )
            self.data_region.src_srs = self.infos.src_srs

            if not self.data_region.valid_p():
                self.data_region = self.region.copy() \
                    if self.transform['trans_region'] is None \
                       else self.transform['trans_region'].copy()
        else:
            self.data_region = regions.Region().from_list(
                self.infos.minmax
            )
            self.data_region.src_srs = self.infos.src_srs
        #     self.region = self.data_region.copy()
        # self.region.zmax=self.upper_limit
        # self.region.zmin=self.lower_limit
        
    def initialize(self):
        self._sub_init()
        self._fn = None # temp filename holder
        self.data_region = None 
        self.inf_region = None 
        self.archive_datalist = None
        self.data_entries = [] 
        self.data_lists = {} 
        
        self.cache_dir = utils.cudem_cache() if self.cache_dir is None else self.cache_dir
        
        if self.sample_alg not in self.gdal_sample_methods: 
            utils.echo_warning_msg(
                f'Alg: {self.sample_alg} is not a valid gdal warp resample algorithm, falling back to `auto`'
            )
            self.sample_alg = 'auto'

        if utils.fn_url_p(self.fn):
            self.remote = True

        ## Initialize transform
        self.transform = {
            'src_horz_crs': None, 'dst_horz_crs': None,
            'src_vert_crs': None, 'dst_vert_crs': None,
            'src_vert_epsg': None, 'dst_vert_epsg': None,
            'pipeline': None, 'trans_fn': None,
            'trans_fn_unc': None, 'trans_region': None,
            'transformer': None, 'vert_transformer': None,
            'want_vertical': False,
        }

        if self.valid_p():
            self.infos = self.inf(
                check_hash=True if self.data_format == -1 else False, make_block_mean=False
            )
            utils.echo_debug_msg(f'Initialize infos: {self.infos}')
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                try:
                    utils.echo_debug_msg(f'Setting transform: {self}: {self.src_srs}')
                    self.set_transform()
                    # srsfun.set_transform(
                    #     src_srs=self.src_srs,
                    #     dst_srs=self.dst_srs,
                    #     region=self.region,
                    #     infos=self.infos,
                    #     cache_dir=self.cache_dir
                    # )
                except Exception as e:
                    utils.echo_error_msg(f'Could not set transformation on {self.fn}, {e}')

            utils.echo_debug_msg(f'Initialize (region): {self.region}')
            self.set_yield(use_blocks=False)

        if self.pnt_fltrs is not None and isinstance(self.pnt_fltrs, str):
            self.pnt_fltrs = [self.pnt_fltrs]

        if self.mask is not None and isinstance(self.mask, str):
            self.mask = [self.mask]
            
        return self

    
    def fetch(self):
        """Fetch remote data from self.data_entries."""
        
        for entry in self.data_entries:
            if entry.remote:
                if entry._fn is None:
                    entry._fn = os.path.basename(self.fn)
                    
                f = fetches.Fetch(url=entry.fn, verbose=entry.verbose)
                if f.fetch_file(entry._fn) == 0:
                    entry.fn = entry._fn
            else:
                utils.echo_warning_msg('Nothing to fetch')

                
    def valid_p(self, fmts=['scratch']):
        """Check if self appears to be a valid dataset entry."""

        if self.fn is None: return False
        if self.data_format is None: return False

        if self.fn is sys.stdin or self.fn == '-':
            self.fn = sys.stdin
            return True
        
        if self.fn is not None:
            if not isinstance(self.fn, (list, np.ndarray, np.core.records.recarray)):
                if self.fn not in fmts:
                    if not self.fn.startswith(('http', '/vsicurl/', 'BAG')) and ':' not in self.fn:
                        if not utils.fn_url_p(self.fn):
                            if self.data_format > -10:
                                if not os.path.exists(self.fn):
                                    utils.echo_warning_msg(f'Fn: {self.fn} does not exist')
                                    return False

                                if os.path.isfile(self.fn) and os.stat(self.fn).st_size == 0:
                                    utils.echo_warning_msg(f'Fn: {self.fn} is 0 bytes')
                                    return False
                        
        return True

    
    def echo_(self, sep=' ', **kwargs):
        """Print self as a datalist entry string."""
        
        out = []
        for key in self.metadata.keys():
            if key != 'name':
                out.append(str(self.metadata[key]))

        return sep.join([f'"{str(x)}"' for x in out])

    
    def echo(self, sep=' ', **kwargs):
        """Print self.data_entries as a datalist entries."""
        
        out = []
        for entry in self.parse():
            entry_path = os.path.abspath(entry.fn) if not self.remote else entry.fn
            l = [entry_path, entry.data_format]
            if entry.weight is not None:
                l.append(entry.weight)

            entry_string = f'{sep.join([str(x) for x in l])}'
            if entry_string:
                out.append(entry_string)
                utils.print_msg(entry_string)

        return out

    
    def format_data(self, sep=' '):
        out = []        
        for entry in self.parse():
            out.append(entry.format_entry(sep=sep))
        return out

    
    def format_entry(self, sep=' '):
        """Format the dataset information as a `sep` separated string."""
        
        dl_entry = sep.join(
            [str(x) for x in [
                self.fn,
                f'{self.data_format}:{factory.dict2args(self.params.get("mod_args", {}))}',
                self.weight,
                self.uncertainty
            ]]
        )
        metadata = self.echo_()
        return sep.join([dl_entry, metadata])


    def format_metadata(self, **kwargs):
        """Format metadata from self, for use as a datalist entry."""
        
        return self.echo_()

    
    def set_yield(self, use_blocks=False):
        """Set the yield strategy, either default (all points), mask, or stacks."""
        
        self.array_yield = self.mask_and_yield_array()
        self.xyz_yield = self.mask_and_yield_xyz()

        if self.region is None and self.x_inc is not None:
            utils.echo_warning_msg('Must enter a region to output in increment blocks...')

        if self.region and self.x_inc:
            self.x_inc = utils.str2inc(self.x_inc)
            if self.y_inc is None:
                self.y_inc = self.x_inc
            else:
                self.y_inc = utils.str2inc(self.y_inc)

            out_name = os.path.join(self.cache_dir, utils.append_fn('globato', self.region, self.x_inc))
            if not use_blocks:
                self.xyz_yield = self.stacks_yield_xyz(out_name=out_name)
            else:
                self.xyz_yield = self.blocks_yield_xyz(out_name=out_name)
    
                
    def inf(self, check_hash=False, recursive_check=False, write_inf=True, 
            make_grid=True, make_block_mean=False, block_inc=None, **kwargs):
        """Read/Write an INF metadata file.

        If the INF file is not found or is invalid, this will attempt to generate one.
        Crucially, this method attempts to populate self.inf_region immediately upon 
        loading an existing file so that sub-modules can access the region during 
        subsequent processing (e.g., hash checks or grid generation).
        """
        
        inf_path = f'{self.fn}.inf'
        generate_inf = False

        ## Try to load existing INF
        if os.path.exists(inf_path):
            try:
                self.infos.load_inf_file(inf_path)
            except ValueError as e:
                generate_inf = True
                utils.remove_glob(inf_path)
                if self.verbose:
                    utils.echo_warning_msg(f'Failed to parse inf {inf_path}, {e}')
        else:
            generate_inf = True        

        ## Check hash if requested (triggers re-generation if mismatch)
        if check_hash and not generate_inf:
            if self.infos.generate_hash() != self.infos.file_hash:
                generate_inf = True
            
        ## Update SRS info
        if self.src_srs is not None:
            self.infos.src_srs = self.src_srs
        else:                                    
            self.src_srs = self.infos.src_srs
            
        try:
            self.inf_region = regions.Region().from_string(self.infos.wkt)
        except:
            try:
                self.inf_region = regions.Region().from_list(self.infos.minmax)
            except:
                # Fallback to user-specified region if INF data is missing/bad
                self.inf_region = self.region.copy() if self.region else None

        ## Streams/fetches modules cannot be scanned/hashed reliably
        if self.fn is sys.stdin or self.data_format < 0:
            generate_inf = False
               
        ## Generate new INF if needed
        if generate_inf:
            ## generate_inf will handle updating self.inf_region internally between passes
            self.infos = self.generate_inf(
                make_grid=make_grid, 
                make_block_mean=make_block_mean, 
                block_inc=block_inc
            )
            
            # Write to disk
            if self.data_format >= -2 and write_inf:
                self.infos.write_inf_file()

            if recursive_check and self.parent is not None:
                self.parent.inf(check_hash=True)

            try:
                self.inf_region = regions.Region().from_string(self.infos.wkt)
            except:
                self.inf_region = regions.Region().from_list(self.infos.minmax)

        ## Set Transform Region (if projecting)
        if self.transform['transformer'] is not None and self.transform['trans_region'] is None:
            if self.inf_region is not None:
                self.transform['trans_region'] = self.inf_region.copy()
                self.transform['trans_region'].src_srs = self.infos.src_srs
                try:
                    self.transform['trans_region'].warp(self.dst_srs)
                except Exception as e:
                    if self.verbose:
                        utils.echo_warning_msg(f"Failed to warp transform region: {e}")

        return self.infos

    
    def _generate_mini_grid(self, region):
        """Internal helper to generate a 10x10 preview grid.
        
        Args:
            region (Region): The region bounding the data.
            
        Returns:
            list: A 10x10 list of values (or None).
        """
        # Setup Gridder
        pp_grid = pointz.PointPixels(src_region=region, x_size=10, y_size=10, verbose=False)
        grid_sum = np.zeros((10, 10))
        grid_count = np.zeros((10, 10))
        
        try:
            # Scan Data
            for points in self.transform_and_yield_points():
                res, srcwin, _ = pp_grid(points, mode='sums')
                if res['z'] is not None and srcwin is not None:
                    # Accumulate chunk into master grid
                    x_off, y_off, x_s, y_s = srcwin
                    y_slice = slice(y_off, y_off + y_s)
                    x_slice = slice(x_off, x_off + x_s)
                    
                    if (y_off + y_s <= 10) and (x_off + x_s <= 10):
                        grid_sum[y_slice, x_slice] += np.nan_to_num(res['z'])
                        grid_count[y_slice, x_slice] += np.nan_to_num(res['count'])

            # Finalize Mean
            with np.errstate(divide='ignore', invalid='ignore'):
                final_grid = grid_sum / grid_count
            
            # Return as list (None for NaNs)
            return np.where(np.isnan(final_grid), None, final_grid).tolist()

        except Exception as e:
            utils.echo_warning_msg(f"Failed to generate mini-grid for {self.fn}: {e}")
            return None

        
    def _generate_block_mean(self, region, block_inc=None):
        """Internal helper to generate a block-mean XYZ file.
        
        Args:
            region (Region): The region bounding the data.
            block_inc (float): The resolution for the block mean.
            
        Returns:
            str: Path to the generated file (or None).
        """
        # Determine filename
        base = os.path.splitext(self.fn)[0]
        block_out_name = f"{base}_blockmean.xyz"
        
        # Determine Increment (~500px wide if not set)
        if block_inc is None:
            width = region.xmax - region.xmin
            block_inc = width / 500.0 if width > 0 else 0.001
            
        try:
            bx, by, _ = region.geo_transform(x_inc=block_inc, y_inc=block_inc)
            pp_block = pointz.PointPixels(src_region=region, x_size=bx, y_size=by, verbose=False)
            
            # Accumulators
            b_arrs = {
                'z_sum': np.zeros((by, bx)),
                'count': np.zeros((by, bx)),
                'x_sum': np.zeros((by, bx)),
                'y_sum': np.zeros((by, bx)),
            }
            
            # Scan Data
            for points in self.transform_and_yield_points():
                res, srcwin, _ = pp_block(points, mode='sums')
                if res['z'] is not None and srcwin is not None:
                    x_off, y_off, x_s, y_s = srcwin
                    y_slice = slice(y_off, y_off + y_s)
                    x_slice = slice(x_off, x_off + x_s)
                    
                    # Safety check for bounds
                    if (y_off + y_s <= by) and (x_off + x_s <= bx):
                        b_arrs['z_sum'][y_slice, x_slice] += np.nan_to_num(res['z'])
                        b_arrs['x_sum'][y_slice, x_slice] += np.nan_to_num(res['x'])
                        b_arrs['y_sum'][y_slice, x_slice] += np.nan_to_num(res['y'])
                        b_arrs['count'][y_slice, x_slice] += np.nan_to_num(res['count'])

            # Finalize & Write
            with np.errstate(divide='ignore', invalid='ignore'):
                fz = b_arrs['z_sum'] / b_arrs['count']
                fx = b_arrs['x_sum'] / b_arrs['count']
                fy = b_arrs['y_sum'] / b_arrs['count']
            
            valid = (b_arrs['count'] > 0) & np.isfinite(fz)
            
            with open(block_out_name, 'w') as f:
                for x, y, z in zip(fx[valid], fy[valid], fz[valid]):
                    f.write(f"{x:.6f} {y:.6f} {z:.6f}\n")
                    
            if self.verbose:
                utils.echo_msg(f"Generated block-mean: {block_out_name}")
            
            return block_out_name

        except Exception as e:
            utils.echo_warning_msg(f"Failed to generate block-mean for {self.fn}: {e}")
            return None

        
    def _wkt_from_mini_grid(self, region, mini_grid):
        """Generate a WKT polygon footprint from the valid cells of the mini-grid.
        
        Args:
            region (Region): The bounding region of the grid.
            mini_grid (list): The 2D list of grid values (None indicates empty cell).
            
        Returns:
            str: A WKT MultiPolygon or Polygon representing the valid data footprint.
        """
        if not mini_grid or not region:
            return region.export_as_wkt()

        try:
            y_size = len(mini_grid)
            x_size = len(mini_grid[0])
            
            # Calculate cell dimensions
            # region width / cols
            x_inc = (region.xmax - region.xmin) / x_size 
            # region height / rows (absolute value for math)
            y_inc = abs(region.ymax - region.ymin) / y_size 
            
            # Collect valid cell polygons
            multipoly = ogr.Geometry(ogr.wkbMultiPolygon)
            has_valid = False
            
            for y in range(y_size):
                for x in range(x_size):
                    if mini_grid[y][x] is not None:
                        has_valid = True
                        
                        # Calculate cell bounds
                        # Assuming Grid (0,0) is Top-Left (standard raster/PointPixels)
                        cell_ymax = region.ymax - (y * y_inc)
                        cell_ymin = region.ymax - ((y + 1) * y_inc)
                        cell_xmin = region.xmin + (x * x_inc)
                        cell_xmax = region.xmin + ((x + 1) * x_inc)
                        
                        # Create Ring (Closed Loop: TL -> TR -> BR -> BL -> TL)
                        ring = ogr.Geometry(ogr.wkbLinearRing)
                        ring.AddPoint(cell_xmin, cell_ymax)
                        ring.AddPoint(cell_xmax, cell_ymax)
                        ring.AddPoint(cell_xmax, cell_ymin)
                        ring.AddPoint(cell_xmin, cell_ymin)
                        ring.AddPoint(cell_xmin, cell_ymax)
                        
                        poly = ogr.Geometry(ogr.wkbPolygon)
                        poly.AddGeometry(ring)
                        multipoly.AddGeometry(poly)
            
            if not has_valid:
                return region.export_as_wkt()
                
            # Merge adjacent cells into a clean footprint
            # UnionCascaded creates a single geometry from the collection
            union_poly = multipoly.UnionCascaded()
            
            return union_poly.ExportToWkt()

        except Exception as e:
            if self.verbose:
                utils.echo_warning_msg(f"Failed to generate tight WKT from mini-grid: {e}")
            return region.export_as_wkt()

        
    def generate_inf(self, make_grid=True, make_block_mean=False, block_inc=None):
        """Generate metadata (INF) for the data source.
        
        Parses the data to calculate point count, min/max bounds (region), 
        and WKT footprint. Optionally generates a mini-grid preview and 
        a block-mean XYZ file.
        
        If a mini-grid is generated, the WKT footprint will be tightened to 
        wrap only the valid grid cells, providing a more accurate spatial representation.
        """
        
        ## Bounds Scan (Pass 1)
        _region_backup = None
        if self.region is not None:
            _region_backup = self.region.copy()
            self.region = None # Unset to scan full file

        _x_inc_backup = None
        if self.x_inc:
            _x_inc_backup = self.x_inc
            self.x_inc = None
            
        this_region = regions.Region()
        point_count = 0

        for points in self.transform_and_yield_points():
            if point_count == 0:
                this_region.from_list([
                    np.min(points['x']), np.max(points['x']),
                    np.min(points['y']), np.max(points['y']),
                    np.min(points['z']), np.max(points['z'])
                ])
            else:
                if np.min(points['x']) < this_region.xmin: this_region.xmin = np.min(points['x'])
                if np.max(points['x']) > this_region.xmax: this_region.xmax = np.max(points['x'])
                if np.min(points['y']) < this_region.ymin: this_region.ymin = np.min(points['y'])
                if np.max(points['y']) > this_region.ymax: this_region.ymax = np.max(points['y'])
                if np.min(points['z']) < this_region.zmin: this_region.zmin = np.min(points['z'])
                if np.max(points['z']) > this_region.zmax: this_region.zmax = np.max(points['z'])
            point_count += len(points)

        ## Populate INF Basic Stats
        self.infos.numpts = point_count
        if point_count > 0:
            self.infos.minmax = this_region.export_as_list(include_z=True)
            # Default WKT to bounding box (will update below if grid is made)
            self.infos.wkt = this_region.export_as_wkt()
        self.infos.src_srs = self.src_srs

        ## Generate Grids (Pass 2)
        if point_count > 0:
            ## Mini Grid
            if make_grid:
                self.infos.mini_grid = self._generate_mini_grid(this_region)
                
                ## UPDATE WKT: Use the mini-grid to create a tighter footprint
                if self.infos.mini_grid is not None:
                    self.infos.wkt = self._wkt_from_mini_grid(this_region, self.infos.mini_grid)
                
            ## Block Mean
            if make_block_mean:
                self._generate_block_mean(this_region, block_inc=block_inc)

        ## Restore Region and inc Constraint
        if _region_backup is not None:
            self.region = _region_backup.copy()

        if _x_inc_backup:
            self.x_inc = _x_inc_backup
                            
        return self.infos

    
    def _generate_grids_from_children(self, children, region, make_grid, make_block_mean, block_inc):
        """Internal helper to generate grids by iterating over children."""
        
        pp_grid = None
        pp_block = None
        grid_arrays = None
        block_arrays = None
        
        if make_grid:
            pp_grid = pointz.PointPixels(src_region=region, x_size=10, y_size=10, verbose=False)
            grid_arrays = {'sum': np.zeros((10, 10)), 'count': np.zeros((10, 10))}

        if make_block_mean:
            if block_inc is None:
                width = region.xmax - region.xmin
                block_inc = width / 500.0 if width > 0 else 0.001
            try:
                bx, by, _ = region.geo_transform(x_inc=block_inc, y_inc=block_inc)
                pp_block = pointz.PointPixels(src_region=region, x_size=bx, y_size=by, verbose=False)
                block_arrays = {
                    'z_sum': np.zeros((by, bx)), 'count': np.zeros((by, bx)),
                    'x_sum': np.zeros((by, bx)), 'y_sum': np.zeros((by, bx))
                }
            except Exception: pass

        def accumulate(master_dict, chunk_res, chunk_srcwin, keys):
            x_off, y_off, x_s, y_s = chunk_srcwin
            y_slice = slice(y_off, y_off + y_s)
            x_slice = slice(x_off, x_off + x_s)
            for k_m, k_r in keys.items():
                if (y_off + y_s <= master_dict[k_m].shape[0]) and (x_off + x_s <= master_dict[k_m].shape[1]):
                    master_dict[k_m][y_slice, x_slice] += np.nan_to_num(chunk_res[k_r])

        for entry in children:
            for points in entry.yield_points():
                if pp_grid:
                    res, srcwin, _ = pp_grid(points, mode='sums')
                    if res['z'] is not None and srcwin is not None:
                        accumulate(grid_arrays, res, srcwin, {'sum': 'z', 'count': 'count'})
                if pp_block:
                    res, srcwin, _ = pp_block(points, mode='sums')
                    if res['z'] is not None and srcwin is not None:
                        accumulate(block_arrays, res, srcwin, {'z_sum': 'z', 'x_sum': 'x', 'y_sum': 'y', 'count': 'count'})

        if make_grid and grid_arrays is not None:
            with np.errstate(divide='ignore', invalid='ignore'):
                final_grid = grid_arrays['sum'] / grid_arrays['count']
            self.infos.mini_grid = np.where(np.isnan(final_grid), None, final_grid).tolist()
            if self.infos.mini_grid:
                self.infos.wkt = self._wkt_from_mini_grid(region, self.infos.mini_grid)

        if make_block_mean and block_arrays is not None:
            try:
                base = os.path.splitext(self.fn)[0]
                block_out = f"{base}_blockmean.xyz"
                with np.errstate(divide='ignore', invalid='ignore'):
                    fz = block_arrays['z_sum'] / block_arrays['count']
                    fx = block_arrays['x_sum'] / block_arrays['count']
                    fy = block_arrays['y_sum'] / block_arrays['count']
                valid = (block_arrays['count'] > 0) & np.isfinite(fz)
                with open(block_out, 'w') as f:
                    for x, y, z in zip(fx[valid], fy[valid], fz[valid]):
                        f.write(f"{x:.6f} {y:.6f} {z:.6f}\n")
            except Exception: pass

            
    def parse(self):
        """Parse the datasets from the dataset."""
        
        if self.region is not None and self.inf_region is not None:
            check_region = self.region if self.transform['trans_region'] is None else self.transform['trans_region']
            if regions.regions_intersect_p(self.inf_region, check_region):
                self.data_entries.append(self)
                yield self
        else:
            self.data_entries.append(self)
            yield self

            
    def parse_data_lists(self, gather_data=True):
        """Parse the data into a datalist dictionary."""
        
        for e in self.parse():
            if e.parent is not None:
                parent_name = e.parent.metadata.get('name')
                if parent_name in self.data_lists:
                    self.data_lists[parent_name]['data'].append(e)
                else:
                    self.data_lists[parent_name] = {'data': [e], 'parent': e.parent}
            else:
                name = e.metadata.get('name')
                self.data_lists[name] = {'data': [e], 'parent': e}
        return self

    
    def yield_xyz_from_entries(self):
        """Yield XYZ from self.data_entries."""
        
        for this_entry in self.data_entries:
            for xyz in this_entry.xyz_yield:
                yield xyz
                
            if this_entry.remote:
                utils.remove_glob(f'{this_entry.fn}*')

                
    def yield_entries(self):
        """Yield from self.data_entries."""
        
        for this_entry in self.data_entries:
            yield this_entry

            
    def yield_points(self):
        """Yield the numpy xyz rec array points from the dataset."""
        
        for ds in self.parse():
            for points in ds.yield_points():
                yield points

                
    def transform_and_yield_points(self):
        """Yield transformed points."""
        
        #with warnings.catch_warnings():
        #warnings.simplefilter('ignore')
        for points in self.yield_points():
            #utils.echo_debug_msg(f'transformer is: {self.transform}')
            if self.transform['transformer'] is not None or self.transform['vert_transformer'] is not None:
                if self.transform['transformer'] is not None:
                    points['x'], points['y'] = self.transform['transformer'].transform(points['x'], points['y'])

                if self.transform['vert_transformer'] is not None:
                    _, _, points['z'] = self.transform['vert_transformer'].transform(points['x'], points['y'], points['z'])

                points = points[~np.isinf(points['z'])]

            if self.region is not None and self.region.valid_p():
                xyz_region = self.region.copy()
                if self.invert_region:
                    ## Inverted Region
                    points = points[
                        ((points['x'] >= xyz_region.xmax) | (points['x'] <= xyz_region.xmin)) | 
                        ((points['y'] >= xyz_region.ymax) | (points['y'] <= xyz_region.ymin))
                    ]
                    if xyz_region.zmin is not None: points = points[(points['z'] <= xyz_region.zmin)]
                    if xyz_region.zmax is not None: points = points[(points['z'] >= xyz_region.zmax)]
                else:
                    # Standard Region
                    points = points[
                        ((points['x'] <= xyz_region.xmax) & (points['x'] >= xyz_region.xmin)) & 
                        ((points['y'] <= xyz_region.ymax) & (points['y'] >= xyz_region.ymin))
                    ]
                    if xyz_region.zmin is not None: points = points[(points['z'] >= xyz_region.zmin)]
                    if xyz_region.zmax is not None: points = points[(points['z'] <= xyz_region.zmax)]

            if self.upper_limit is not None: points = points[(points['z'] <= self.upper_limit)]
            if self.lower_limit is not None: points = points[(points['z'] >= self.lower_limit)]

            if len(points) > 0:
                if isinstance(self.pnt_fltrs, list):
                    for f in self.pnt_fltrs:
                        ## set verbosity in the mod
                        point_filter = pointz.PointFilterFactory(
                            mod=f,
                            points=points,
                            region=self.region,
                            xyinc=[self.x_inc, self.y_inc],
                            cache_dir=self.cache_dir,
                            verbose=False,
                        )._acquire_module()
                        if point_filter:
                            points = point_filter()

                if len(points) > 0:
                    yield points

        self.transform['transformer'] = self.transform['vert_transformer'] = None

        
    def mask_and_yield_array(self):
        """Mask the incoming array from `self.yield_array` and yield the results.
        Points within the mask will be retained, unless data_mask['invert'] is True.
        """
        
        mask_band = None
        mask_infos = None
        mask_count = 0
        data_masks = []
        if self.mask is not None:
            for mask in self.mask:
                opts = self._init_mask(mask)
                data_masks.append(opts)

        for out_arrays, this_srcwin, this_gt in self.yield_array():
            for data_mask in data_masks:
                if data_mask is None: continue
                
                src_mask = None
                if os.path.exists(data_mask['data_mask']):
                    utils.echo_debug_msg(f'Using mask dataset: {data_mask} to array')                        
                    src_mask = gdal.Open(data_mask['data_mask'])
                    if src_mask is not None:
                        mask_band = src_mask.GetRasterBand(1)
                        mask_infos = gdalfun.gdal_infos(src_mask)
                    else:
                        mask_band = None
                else:
                    utils.echo_warning_msg(f'Could not load mask {data_mask["data_mask"]}')

                ## Z-Mask
                z_mask = np.full(out_arrays['z'].shape, True)
                if data_mask['min_z'] is not None or data_mask['max_z'] is not None:
                    if data_mask['min_z'] is not None and data_mask['max_z'] is not None:
                        z_mask = ((out_arrays['z'] > data_mask['min_z']) & (out_arrays['z'] < data_mask['max_z']))
                    elif data_mask['min_z'] is not None:
                        z_mask = out_arrays['z'] > data_mask['min_z']
                    else:
                        z_mask = out_arrays['z'] < data_mask['max_z']

                if mask_band is not None:
                    mask_data = mask_band.ReadAsArray(*this_srcwin)
                    if mask_data is None or len(mask_data) == 0: continue
                        
                    if not np.isnan(mask_infos['ndv']):
                        mask_data[mask_data == mask_infos['ndv']] = np.nan
                        
                    out_mask = ((~np.isnan(mask_data)) & (z_mask))                    

                    try:
                        for arr in out_arrays.keys():
                            if arr not in ['pixel_x', 'pixel_y']:
                                if out_arrays[arr] is not None:
                                    val = 0 if arr == 'count' else np.nan
                                    if data_mask['invert']:
                                        out_arrays[arr][~out_mask] = val
                                    else:                                    
                                        out_arrays[arr][out_mask] = val
                    except Exception as e:
                        utils.echo_error_msg(f'could not mask array: {arr}, {e}')

                    mask_count += np.count_nonzero(~out_mask) if data_mask['invert'] else np.count_nonzero(out_mask)                        
                    mask_data = None
                    src_mask = mask_band = None

            if mask_count > 0 and self.verbose:
                utils.echo_msg_bold(f'Masked {mask_count} data records from {self.fn}')
            yield out_arrays, this_srcwin, this_gt

            
    def mask_and_yield_xyz(self):
        """Mask the incoming xyz data from `self.yield_xyz` and yield the results."""
        
        for this_entry in self.parse():
            if this_entry.mask is None:
                for this_xyz in this_entry.yield_xyz():
                    yield this_xyz
            else:
                data_masks = []
                mask_count = 0
                for mask in this_entry.mask:
                    opts = self._init_mask(mask)
                    data_masks.append(opts)

                utils.echo_debug_msg(f'Using mask dataset: {data_masks} to xyz')                    
                for this_xyz in this_entry.yield_xyz():
                    masked = False
                    for data_mask in data_masks:                                
                        if data_mask is None: continue

                        z_masked = True
                        if data_mask['min_z'] is not None or data_mask['max_z'] is not None:
                            if data_mask['min_z'] is not None and data_mask['max_z'] is not None:
                                z_masked = (this_xyz.z > data_mask['min_z']) & (this_xyz.z < data_mask['max_z'])
                            elif data_mask['min_z'] is not None:
                                z_masked = this_xyz.z > data_mask['min_z']
                            else:
                                z_masked = this_xyz.z < data_mask['max_z']

                        if data_mask['ogr_or_gdal'] == 0: # Raster Mask
                            src_ds = gdal.Open(data_mask['data_mask'])
                            if src_ds is not None:
                                ds_config = gdalfun.gdal_infos(src_ds)
                                ds_band = src_ds.GetRasterBand(1)
                                
                                xpos, ypos = utils._geo2pixel(
                                    this_xyz.x, this_xyz.y, ds_config['geoT'], node='pixel'
                                )
                                
                                val_masked = False
                                if 0 <= xpos < ds_config['nx'] and 0 <= ypos < ds_config['ny']:
                                    tgrid = ds_band.ReadAsArray(xpos, ypos, 1, 1)
                                    if tgrid is not None:
                                        if not np.isnan(ds_config['ndv']):
                                            tgrid[tgrid == ds_config['ndv']] = np.nan

                                        if np.isnan(tgrid[0][0]):
                                            val_masked = False
                                        else:
                                            val_masked = True
                                
                                if data_mask['invert']:
                                    masked = (not val_masked) & z_masked
                                else:
                                    masked = val_masked & z_masked
                                
                                if masked: mask_count += 1

                        else: # Vector Mask
                            src_ds = ogr.Open(data_mask['data_mask'])
                            if src_ds is not None:
                                layer = src_ds.GetLayer()
                                geomcol = gdalfun.ogr_union_geom(layer, verbose=False)
                                if not geomcol.IsValid():
                                    geomcol = geomcol.Buffer(0)

                                this_point = ogr.CreateGeometryFromWkt(this_xyz.export_as_wkt())
                                is_within = this_point.Within(geomcol)
                                
                                if data_mask['invert']:
                                    masked = (not is_within) & z_masked
                                else:
                                    masked = is_within & z_masked
                                
                                if masked: mask_count += 1
                                src_ds = layer = None

                    if not masked:
                        yield this_xyz

                    if mask_count > 0 and self.verbose:
                        utils.echo_msg_bold(f'Masked {mask_count} data records from {self.fn}')        

                        
    def yield_xyz(self):
        """Yield the data as xyz points."""
        
        count = 0
        for points in self.transform_and_yield_points():
            if 'w' in points.dtype.names:
                points_w = points['w']
            else:
                points_w = np.ones(points['z'].shape).astype(float)

            points_w *= self.weight if self.weight is not None else 1.
            points_w[np.isnan(points_w)] = 1
                
            if 'u' in points.dtype.names:
                points_u = points['u']
            else:
                points_u = np.zeros(points['z'].shape)
                
            points_u = np.sqrt(
                points_u**2 + (self.uncertainty if self.uncertainty is not None else 0)**2
            )
            points_u[np.isnan(points_u)] = 0
            
            dataset = np.vstack((points['x'], points['y'], points['z'], points_w, points_u)).transpose()
            count += len(dataset)
            
            for point in dataset:
                yield xyzfun.XYZPoint(
                    x=point[0], y=point[1], z=point[2], w=point[3], u=point[4]
                )

        if self.verbose:
            utils.echo_msg_bold(f'Parsed {count} data records from {self.fn} @ a weight of {self.weight}')

            
    def yield_array(self, want_sums=True):
        """Yield the data as an array."""
        
        count = 0
        for points in self.transform_and_yield_points():
            count += points.size
            xcount, ycount, _ = self.region.geo_transform(
                x_inc=self.x_inc, y_inc=self.y_inc, node='grid'
            )
            
            point_array = pointz.PointPixels(
                src_region=self.region,
                x_size=xcount,
                y_size=ycount,
                verbose=self.verbose
            )
            yield point_array(
                points,
                weight=self.weight,
                uncertainty=self.uncertainty,
                mode='sums' if want_sums else 'mean'
            )

        if self.verbose:
            utils.echo_msg_bold(f'Parsed {count} data records from {self.fn} @ a weight of {self.weight}')    

            
    def stacks(self, out_name=None, use_blocks=False):
        if not use_blocks:
            gbt = globato.GdalRasterStacker(
                region=self.region,
                x_inc=self.x_inc,
                y_inc=self.y_inc,
                dst_srs=self.dst_srs,
                cache_dir=self.cache_dir
            )
            stacked_fn = gbt.process_stack(self.parse(), out_name=out_name)
        else:
            gbt = globato.GlobatoStacker(
                region=self.region,
                x_inc=self.x_inc,
                y_inc=self.y_inc,
                dst_srs=self.dst_srs,
                cache_dir=self.cache_dir
            )
            blocked_fn = gbt.process_blocks(self.parse(), out_name=out_name)
            stacked_fn = globato_converter.globato_to_gdal(blocked_fn, tif_path=f'{out_name}.tif', verbose=True)
        
        ## Perform any desired grits filters on the stack (making a copy, so we can retain the original).
        if self.stack_fltrs:
            for f in self.stack_fltrs:
                utils.echo_msg(f'Filtering stacks module with {f}')
                grits_filter = grits.GritsFactory(
                    mod=f,
                    src_dem=stacked_fn,
                    uncertainty_mask=4,
                    weight_mask=3,
                    count_mask=2,
                    cache_dir=self.cache_dir,
                    verbose=False
                )._acquire_module()
            
                if grits_filter:
                    grits_filter = grits_filter()
                    if os.path.exists(grits_filter.dst_dem):
                        stacked_fn = grits_filter.dst_dem
                    else:
                        utils.echo_warning_msg('Grits output in invalid: {grits_filter.dst_dem}')

        return stacked_fn

    
    def stacks_yield_xyz(self, out_name=None, ndv=-9999, fmt='GTiff'):
        """Yield the result of `stacks` as an xyz object."""
        
        stacked_fn = self.stacks(out_name=out_name)
        sds = gdal.Open(stacked_fn)
        sds_gt = sds.GetGeoTransform()
        
        ## Bands: 1:Z, 3:W, 4:U, 6:X, 7:Y
        sds_z_band = sds.GetRasterBand(1)
        sds_w_band = sds.GetRasterBand(3)
        sds_u_band = sds.GetRasterBand(4)
        sds_x_band = sds.GetRasterBand(6)
        sds_y_band = sds.GetRasterBand(7)
        
        srcwin = (0, 0, sds.RasterXSize, sds.RasterYSize)
        for y in range(srcwin[1], srcwin[1] + srcwin[3], 1):
            sz = sds_z_band.ReadAsArray(srcwin[0], y, srcwin[2], 1)
            if np.all(sz == ndv): continue
            
            sw = sds_w_band.ReadAsArray(srcwin[0], y, srcwin[2], 1)
            su = sds_u_band.ReadAsArray(srcwin[0], y, srcwin[2], 1)
            sx = sds_x_band.ReadAsArray(srcwin[0], y, srcwin[2], 1)
            sy = sds_y_band.ReadAsArray(srcwin[0], y, srcwin[2], 1)
            
            for x in range(0, sds.RasterXSize):
                z = sz[0, x]
                if z != ndv:
                    if self.stack_node:
                        yield xyzfun.XYZPoint(x=sx[0, x], y=sy[0, x], z=z, w=sw[0, x], u=su[0, x])
                    else: 
                        geo_x, geo_y = utils._pixel2geo(x, y, sds_gt)
                        yield xyzfun.XYZPoint(x=geo_x, y=geo_y, z=z, w=sw[0, x], u=su[0, x])
                    
        sds = None


    def blocks(self, out_name=None):
        gbt = globato.GlobatoStacker(
            region=self.region,
            x_inc=self.x_inc,
            y_inc=self.y_inc,
            dst_srs=self.dst_srs,
            cache_dir=self.cache_dir
        )

        blocked_fn = gbt.process_blocks(self.parse(), out_name=out_name)        
        
        return blocked_fn

    
    def blocks_yield_xyz(self, out_name=None):
        """Yield the result of `blocks` as an xyz object."""

        stacked_fn = blocks(out_name=out_name)
        sds = h5.File(stacked_fn, 'r')
        sds_gt = [float(x) for x in sds['crs'].attrs['GeoTransform'].split()]
        sds_stack = sds['stack'] # Assuming final stack is in 'stack' group
        
        stack_shape = sds_stack['z'].shape
        srcwin = (0, 0, stack_shape[1], stack_shape[0])
        
        for y in range(srcwin[1], srcwin[1] + srcwin[3], 1):
            sz = sds_stack['z'][y:y+1, srcwin[0]:srcwin[0]+srcwin[2]]
            if np.all(np.isnan(sz)): continue

            sw = sds_stack['weights'][y:y+1, srcwin[0]:srcwin[0]+srcwin[2]]
            su = sds_stack['uncertainty'][y:y+1, srcwin[0]:srcwin[0]+srcwin[2]]
            
            ## Use X/Y arrays if available or calc from geo
            sx = sds_stack['x'][y:y+1, srcwin[0]:srcwin[0]+srcwin[2]]
            sy = sds_stack['y'][y:y+1, srcwin[0]:srcwin[0]+srcwin[2]]
            
            for x in range(0, stack_shape[1]):
                z = sz[0, x]
                if not np.isnan(z):
                    if self.stack_node:
                        yield xyzfun.XYZPoint(x=sx[0, x], y=sy[0, x], z=z, w=sw[0, x], u=su[0, x])
                    else:
                        geo_x, geo_y = utils._pixel2geo(x, y, sds_gt)
                        yield xyzfun.XYZPoint(x=geo_x, y=geo_y, z=z, w=sw[0, x], u=su[0, x])
        sds.close()

        
    ## ==============================================
    ## Data Dump/Export/Archive
    ## ==============================================    
    def _xyz_dump(self, this_xyz, dst_port=sys.stdout, encode=False):
        this_xyz.dump(
            include_w=True if self.weight is not None else False,
            include_u=True if self.uncertainty is not None else False,
            dst_port=dst_port,
            encode=encode,
            precision=self.dump_precision
        )

        
    def dump_xyz(self, dst_port=sys.stdout, encode=False):
        """Dump the XYZ data from the dataset."""
        
        for this_xyz in self.xyz_yield:
            self._xyz_dump(this_xyz, dst_port=dst_port, encode=encode)

            
    def dump_xyz_direct(self, dst_port=sys.stdout, encode=False):
        """Dump the XYZ data directly bypassing filters."""
        
        for this_xyz in self.yield_xyz():
            self._xyz_dump(this_xyz, dst_port=dst_port, encode=encode)

            
    def export_xyz_as_list(self, z_only=False):
        """Return the XYZ data as a python list."""
        
        return [xyz.z if z_only else xyz.copy() for xyz in self.xyz_yield]

    
    def export_xyz_as_ogr(self):
        """Make a point vector OGR DataSet Object from src_xyz."""
        
        dst_ogr = self.metadata.get('name', 'dataset')
        driver = gdal.GetDriverByName('Memory')
        ogr_ds = driver.Create('', 0, 0, 0, gdal.GDT_Unknown)
        layer = ogr_ds.CreateLayer(dst_ogr, geom_type=ogr.wkbPoint25D)
        
        for x in ['x', 'y', 'z', 'weight', 'uncertainty']:
            fd = ogr.FieldDefn(x, ogr.OFTReal)
            fd.SetWidth(12)
            fd.SetPrecision(8)
            layer.CreateField(fd)
            
        f_defn = layer.GetLayerDefn()
        for this_xyz in self.xyz_yield:
            f = ogr.Feature(f_defn)
            f.SetField(0, this_xyz.x)
            f.SetField(1, this_xyz.y)
            f.SetField(2, this_xyz.z)
            f.SetField(3, this_xyz.w)
            f.SetField(3, this_xyz.u)
                
            wkt = this_xyz.export_as_wkt(include_z=True)
            g = ogr.CreateGeometryFromWkt(wkt)
            f.SetGeometryDirectly(g)
            layer.CreateFeature(f)
            
        return ogr_ds

    
    def export_xyz_as_pandas(self):
        """Export the point data as a pandas dataframe."""
        
        import pandas as pd
        
        frames = []
        for points in self.transform_and_yield_points():
            points_w = points['w'] if 'w' in points.dtype.names else np.ones(points['z'].shape)
            points_w *= self.weight if self.weight is not None else 1
            points_w[np.isnan(points_w)] = 1
                
            points_u = points['u'] if 'u' in points.dtype.names else np.zeros(points['z'].shape)
            points_u = np.sqrt(points_u**2 + (self.uncertainty if self.uncertainty is not None else 0)**2)
            points_u[np.isnan(points_u)] = 0
            
            frames.append(pd.DataFrame({
                'x': points['x'], 'y': points['y'], 'z': points['z'],
                'weight': points_w, 'uncertainty': points_u
            }))
            
        if frames:
            return pd.concat(frames, ignore_index=True)
        return pd.DataFrame(columns=['x', 'y', 'z', 'weight', 'uncertainty'])

    
    def export_xyz_as_numpy(self):
        """Export the point data as a numpy structured array.
        
        Returns:
            np.ndarray: A structured array with fields:
                        ['x', 'y', 'z', 'weight', 'uncertainty']
        """
        
        chunks = []
        
        ## Define the output data type
        dt = [('x', 'f8'), ('y', 'f8'), ('z', 'f8'), 
              ('weight', 'f8'), ('uncertainty', 'f8')]

        for points in self.transform_and_yield_points():
            ## Handle Weights
            if 'w' in points.dtype.names:
                points_w = points['w'].astype(np.float64)
            else:
                points_w = np.ones(points['z'].shape, dtype=np.float64)

            if self.weight is not None:
                points_w *= self.weight
            
            points_w[np.isnan(points_w)] = 1.0

            ## Handle Uncertainty
            if 'u' in points.dtype.names:
                points_u = points['u'].astype(np.float64)
            else:
                points_u = np.zeros(points['z'].shape, dtype=np.float64)

            global_u = self.uncertainty if self.uncertainty is not None else 0
            points_u = np.sqrt(points_u**2 + global_u**2)
            points_u[np.isnan(points_u)] = 0.0

            ## Create Structured Array Chunk
            count = len(points['z'])
            chunk = np.zeros(count, dtype=dt)
            chunk['x'] = points['x']
            chunk['y'] = points['y']
            chunk['z'] = points['z']
            chunk['weight'] = points_w
            chunk['uncertainty'] = points_u
            
            chunks.append(chunk)

        ## Concatenate
        if chunks:
            return np.concatenate(chunks)
        
        return np.array([], dtype=dt)
    
    
    ## ==============================================
    ## Data archive
    ## ==============================================            
    # def archive_xyz(self, **kwargs):
    #     """Archive data from the dataset to XYZ in the given dataset region.
        
    #     will convert all data to XYZ within the given region and will arrange
    #     the data as a datalist based on inputs.

    #     Data comes from `self.xyz_yield` set in `self.set_yield`. 
    #     So will process data through `stacks` before archival if region, 
    #     x_inc and y_inc are set.
    #     """
        
    #     srs_all = []
    #     a_name = None
    #     if 'dirname' in kwargs.keys() and kwargs['dirname'] is not None:
    #         a_name = kwargs['dirname']
    #     #else:
    #     aa_name = self.metadata['name'].split('/')[0]
    #     if self.region is None:
    #         aa_name = '{}_{}'.format(aa_name, utils.this_year())
    #     else:
    #         aa_name = '{}_{}_{}'.format(
    #             aa_name, self.region.format('fn'), utils.this_year())

    #     if a_name is None:
    #         a_name = aa_name
    #         self.archive_datalist = '{}.datalist'.format(a_name)
    #     else:
    #         self.archive_datalist = '{}.datalist'.format(a_name)
    #         a_name = os.path.join(a_name, aa_name)

    #     utils.echo_msg(self.archive_datalist)
    #     if not os.path.exists(os.path.dirname(self.archive_datalist)):
    #         try:
    #             os.makedirs(os.path.dirname(self.archive_datalist))
    #         except:
    #             pass

    #     archive_keys = []
    #     for this_entry in self.parse():
    #         srs_all.append(this_entry.dst_srs \
    #                        if this_entry.dst_srs is not None \
    #                        else this_entry.src_srs)
    #         datalist_dirname = os.path.join(
    #             a_name, os.path.dirname(this_entry.metadata['name'])
    #         )
    #         this_key = datalist_dirname.split('/')[-1]
    #         if not os.path.exists(datalist_dirname):
    #             os.makedirs(datalist_dirname)

    #         sub_datalist = os.path.join(datalist_dirname, f'{this_key}.datalist')
    #         if utils.fn_basename2(os.path.basename(this_entry.fn)) == '':
    #             sub_xyz_path = '.'.join(
    #                 [utils.fn_basename2(os.path.basename(this_entry.metadata['name'])),
    #                  'xyz']
    #             )
    #         else:
    #             sub_xyz_path = '.'.join(
    #                 [utils.fn_basename2(os.path.basename(this_entry.fn)),
    #                  'xyz']
    #             )

    #         this_xyz_path = os.path.join(datalist_dirname, sub_xyz_path)
    #         if os.path.exists(this_xyz_path):
    #             utils.echo_warning_msg(
    #                 f'{this_xyz_path} already exists, skipping...'
    #             )
    #             continue

    #         with open(sub_datalist, 'a') as sub_dlf:
    #             with open(this_xyz_path, 'w') as xp:
    #                 ## data will be processed independently of each other
    #                 for this_xyz in this_entry.xyz_yield: 
    #                     this_xyz.dump(
    #                         include_w=True if self.weight is not None else False,
    #                         include_u=True if self.uncertainty is not None else False,
    #                         dst_port=xp, encode=False, precision=self.dump_precision
    #                     )

    #             if os.stat(this_xyz_path).st_size > 0:
    #                 sub_dlf.write(f'{sub_xyz_path} 168 1 0\n')
    #             else:
    #                 utils.remove_glob(f'{this_xyz_path}*')

    #         if os.stat(sub_datalist).st_size == 0:
    #             utils.remove_glob(sub_datalist)
    #         else:
    #             with open(self.archive_datalist, 'a') as dlf:
    #                 if not this_key in archive_keys:
    #                     archive_keys.append(this_key)                    
    #                     dlf.write(
    #                         '{name}.datalist -1 {weight} {uncertainty} {metadata}\n'.format(
    #                             name=os.path.join(datalist_dirname, this_key),
    #                             weight=utils.float_or(this_entry.weight, 1),
    #                             uncertainty=utils.float_or(this_entry.uncertainty, 0),
    #                             metadata=this_entry.format_metadata()
    #                         )
    #                     )

    #         if not os.listdir(datalist_dirname):
    #             os.rmdir(datalist_dirname)
                    
    #     ## generate datalist inf/json
    #     srs_set = set(srs_all)
    #     if len(srs_set) == 1:
    #         arch_srs = srs_set.pop()
    #     else:
    #         arch_srs = None

    #     utils.remove_glob(f'{self.archive_datalist}.*')
    #     this_archive = DatasetFactory(
    #         mod=self.archive_datalist,
    #         data_format=-1,
    #         parent=None,
    #         weight=1,
    #         uncertainty=0,
    #         src_srs=arch_srs,
    #         dst_srs=None,
    #         cache_dir=self.cache_dir,
    #     )._acquire_module().initialize()
        
    #     return this_archive.inf()


    def archive_xyz(self, dirname: str = None, **kwargs):
        """Archive data from the dataset to XYZ format within the given region.
        
        This method converts all data to XYZ, organizes it into a directory structure,
        and generates a master datalist pointing to the archived files.
        
        Args:
            dirname/dst_fn (str): The root directory for the archive. Defaults to the dataset name.
                                  If the dirname ends in .tar.gz or .tgz the archive will be tar.gz'd.
            **kwargs: Additional arguments.
            
        Returns:
            INF: The metadata info object for the generated archive datalist.
        """

        import tarfile
        import shutil
                
        ## Determine Archive Root Directory
        archive_root = dirname if dirname is not None else None
        dst_fn = None
        
        ## Default name based on metadata
        dataset_name = self.metadata.get('name', 'dataset').split('/')[0]
        
        if self.region is None:
            archive_name = f"{dataset_name}_{utils.this_year()}"
        else:
            archive_name = f"{dataset_name}_{self.region.format('fn')}_{utils.this_year()}"

        if archive_root is None:            
            archive_root = archive_name
            self.archive_datalist = f"{archive_root}.datalist"
        else:
            ## Check if we are creating a Tarball and set it as the dst_fn
            is_tar = archive_root.endswith('.tar.gz') or archive_root.endswith('.tgz')
            if is_tar:
                dst_fn = archive_root
                archive_root = archive_name

            self.archive_datalist = f"{archive_root}.datalist"
            archive_root = os.path.join(archive_root, archive_name)

        if is_tar:
            utils.echo_msg(f"Archiving to: {dst_fn}")
        else:
            utils.echo_msg(f"Archiving to: {self.archive_datalist}")
            
        if not os.path.exists(os.path.dirname(self.archive_datalist)):
            try:
                os.makedirs(os.path.dirname(self.archive_datalist))
            except OSError:
                pass

        if not os.path.exists(self.archive_datalist):
            utils.touch(self.archive_datalist)

        ## Process Entries
        archive_keys = []
        srs_all = []

        ## Iterate through parsed sub-datasets
        for this_entry in self.parse():
            ## Track SRS for final metadata
            srs = this_entry.dst_srs if this_entry.dst_srs is not None else this_entry.src_srs
            if srs: srs_all.append(srs)

            ## Determine sub-directory for this specific entry
            entry_name = this_entry.metadata.get('name', 'unknown')
            datalist_dirname = os.path.join(archive_root, os.path.dirname(entry_name))
            this_key = os.path.basename(datalist_dirname)
            
            if not os.path.exists(datalist_dirname):
                os.makedirs(datalist_dirname)

            ## Define sub-datalist and xyz paths
            sub_datalist_path = os.path.join(datalist_dirname, f"{this_key}.datalist")
            
            ## Construct XYZ filename
            base_fn = os.path.basename(this_entry.fn)
            clean_name = utils.fn_basename2(base_fn)
            if not clean_name:
                clean_name = utils.fn_basename2(os.path.basename(entry_name))
            
            sub_xyz_filename = f"{clean_name}.xyz"
            this_xyz_path = os.path.join(datalist_dirname, sub_xyz_filename)

            if os.path.exists(this_xyz_path):
                utils.echo_warning_msg(f'Path: {this_xyz_path} already exists, skipping...')
                continue

            ## Write Data
            ## Write XYZ data
            with open(this_xyz_path, 'w') as xp:
                for this_xyz in this_entry.xyz_yield: 
                    this_xyz.dump(
                        include_w=True if self.weight is not None else False,
                        include_u=True if self.uncertainty is not None else False,
                        dst_port=xp, 
                        encode=False, 
                        precision=self.dump_precision
                    )

            ## Update Sub-Datalist
            ## If data was written, append to sub-datalist. Else clean up.
            if os.path.exists(this_xyz_path) and os.stat(this_xyz_path).st_size > 0:
                with open(sub_datalist_path, 'a') as sub_dlf:
                    # Format: path format weight uncertainty
                    sub_dlf.write(f'{sub_xyz_filename} 168 1 0\n')
            else:
                utils.remove_glob(f'{this_xyz_path}*')

            ## Update Master Datalist
            ## If the sub-datalist has content, add it to the master archive list
            if os.path.exists(sub_datalist_path) and os.stat(sub_datalist_path).st_size > 0:
                if this_key not in archive_keys:
                    archive_keys.append(this_key)
                    with open(self.archive_datalist, 'a') as dlf:
                        ## Add entry pointing to the sub-datalist
                        ## Format: path -1 weight uncertainty metadata...
                        rel_path = os.path.join(datalist_dirname, this_key)
                        meta_str = this_entry.format_metadata()
                        w = utils.float_or(this_entry.weight, 1)
                        u = utils.float_or(this_entry.uncertainty, 0)
                        
                        dlf.write(f'{rel_path}.datalist -1 {w} {u} {meta_str}\n')
            else:
                utils.remove_glob(sub_datalist_path)

            ## Cleanup empty directories
            if os.path.exists(datalist_dirname) and not os.listdir(datalist_dirname):
                os.rmdir(datalist_dirname)
                    
        ## Finalize Archive Metadata
        ## Determine common SRS
        srs_set = set(srs_all)
        arch_srs = srs_set.pop() if len(srs_set) == 1 else None

        ## Generate INF for the master datalist
        #utils.remove_glob(f'{self.archive_datalist}.inf')
        
        ## Use Factory to load the new archive and generate its INF
        this_archive = DatasetFactory(
            mod=self.archive_datalist,
            data_format=-1,
            parent=None,
            weight=1,
            uncertainty=0,
            src_srs=arch_srs,
            dst_srs=None,
            cache_dir=self.cache_dir,
        )._acquire_module().initialize()

        this_archive.generate_inf()
        ## Create the Tarball
        if is_tar:
            utils.echo_msg(f"Archiving to {dst_fn}...")
            with tarfile.open(dst_fn, "w:gz") as tar:
                ## Add the Master Datalist to the Root
                if self.archive_datalist and os.path.exists(self.archive_datalist):
                    tar.add(self.archive_datalist, arcname=os.path.basename(self.archive_datalist))
                    
                ## Add the Data directory
                tar.add(data_dir, arcname=archive_root)
                
            utils.echo_msg(f"Archive created: {dst_fn}")
        
        return this_archive.inf()
    
    
    ## ==============================================
    ## Fetching
    ## ==============================================
    def fetch_data(self, fetches_module, check_size=True):
        this_fetches = fetches.FetchesFactory(
            mod=fetches_module,
            src_region=self.region,
            verbose=self.verbose,
            outdir=self.cache_dir,
            callback=fetches.fetches_callback
        )._acquire_module()        
        this_fetches.run()
        fr = fetches.fetch_results(this_fetches, check_size=check_size)
        fr.daemon = True
        fr.start()
        fr.join()

        return fr
    
                
class CUDEMFile(ElevationDataset):
    """CUDEM netcdf raster

    the cudem netcdf/h5 contains uninterpolated elevation data, 
    uncertainty weight and data mask...
    """
    
    def __init__(self, stack=False, **kwargs):
        super().__init__(**kwargs)
        self.stack = stack


    def yield_points_nc(self):

        ## extract z, unc and weight grids
        ## process through gdalfile
        
        nc_data = nc.Dataset(self.fn, 'r')
        stack_lat = nc_data['/stack/latitude'][...,]
        stack_lon = nc_data['/stack/longitude'][...,]
        stack_h = nc_data['/stack/stack_h'][...,]

        stack_u = nc_data['/stack/uncertainty'][...,]
        stack_w = nc_data['/stack/weight'][...,]
        counts = nc_data['/stack/count'][...,]

        nc_data.close()

        dataset = np.column_stack(
            (stack_lon.data, stack_lat.data,
             stack_h.data[0], stack_w.data[0],
             stack_u.data[0])
        )
        points = np.rec.fromrecords(
            dataset, names='x, y, z, w, u'
        )
        #points = points[points.z != -9999]
        
        yield(points)

    def yield_points(self):
        data = h5.File(self.fn)
        stack_lat = data['/stack/y'][...,]
        stack_lon = data['/stack/x'][...,]
        stack_h = data['/stack/z'][...,]

        stack_u = data['/stack/uncertainty'][...,]
        stack_w = data['/stack/weights'][...,]
        counts = data['/stack/count'][...,]

        data.close()

        nan_m = np.isnan(stack_h)
        stack_h = stack_h[~nan_m].flatten()                        
        stack_lat = stack_lat[~nan_m].flatten()                        
        stack_lon = stack_lon[~nan_m].flatten()                        
        stack_u = stack_u[~nan_m].flatten()                        
        stack_w = stack_w[~nan_m].flatten()                        
        
        dataset = np.column_stack(
            (stack_lon, stack_lat,
             stack_h, stack_w,
             stack_u)
        )
        points = np.rec.fromrecords(
            dataset, names='x, y, z, w, u'
        )
        #points = points[points.z != -9999]
        
        yield(points)
                

class FactoryDatalists(ElevationDataset):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        
class CRMDatalist(FactoryDatalists):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)    


class DatasetFactory(factory.CUDEMFactory):
    """Dataset Factory Settings and Generator.
    
    Parses a datalist entry (string) and returns the appropriate 
    dataset object based on the format ID or file extension.
    """

    ## Import dataset modules here to avoid circular imports
    from . import xyzfile
    from . import lasfile
    from . import gdalfile
    from . import mbsfile
    from . import ogrfile
    from . import icesat2file
    from . import swotfile
    from . import ziplistfile
    from . import fetchers
    from . import datalistfile
    from . import globatofile
    
    ## Registry of supported dataset modules
    ## Negative IDs are containers/lists/fetchers
    ## Positive IDs are specific data formats
    _modules = {
        ## --- Utility ---
        0: {'name': 'auto', 'fmts': ['auto'], 'description': 'Automatic Format Detection', 'call': None},
        
        ## --- Lists and Containers ---
        -1: {
            'name': 'datalist',
            'fmts': ['datalist', 'mb-1', 'dl'],
            'description': 'An extended MB-System style datalist containing dlim-compatible datasets',
            'call': datalistfile.Datalist
        },
        -2: {
            'name': 'zip',
            'fmts': ['zip', 'ZIP', 'tar', 'gz', '7z'],
            'description': 'A compressed archive containing dlim-compatible datasets',
            'call': ziplistfile.ZIPlist
        },
        -3: {
            'name': 'scratch',
            'fmts': [],
            'description': 'A scratch dataset, including a python list of dlim-compatible datasets',
            'call': datalistfile.Scratch
        },
        -4: {
            'name': 'points',
            'fmts': [],
            'description': 'A points dataset (numpy rec-array or pandas dataframe)',
            'call': datalistfile.Points
        },
        
        ## --- Fetchers ---
        -100: {'name': 'https', 'fmts': ['https'], 'description': 'URL to dataset', 'call': fetchers.Fetcher},
        -101: {'name': 'gmrt', 'fmts': ['gmrt'], 'description': 'GMRT Fetcher', 'call': fetchers.GMRTFetcher},
        -102: {'name': 'gebco', 'fmts': ['gebco'], 'description': 'GEBCO Fetcher', 'call': fetchers.GEBCOFetcher},
        -103: {'name': 'copernicus', 'fmts': ['copernicus'], 'description': 'Copernicus Fetcher', 'call': fetchers.CopernicusFetcher},
        -104: {'name': 'fabdem', 'fmts': ['fabdem'], 'description': 'FABDEM Fetcher', 'call': fetchers.FABDEMFetcher},
        -105: {'name': 'nasadem', 'fmts': ['nasadem'], 'description': 'NASADEM Fetcher', 'call': fetchers.Fetcher},
        -106: {'name': 'mar_grav', 'fmts': ['mar_grav'], 'description': 'MarGrav Fetcher', 'call': fetchers.MarGravFetcher},
        -107: {'name': 'srtm_plus', 'fmts': ['srtm_plus'], 'description': 'SRTM+ Fetcher', 'call': fetchers.Fetcher},
        -108: {'name': 'synbath', 'fmts': ['synbath'], 'description': 'SynBath Fetcher', 'call': fetchers.Fetcher},
        -109: {'name': 'gedtm30', 'fmts': ['gedtm30'], 'description': 'Global DTM Fetcher', 'call': fetchers.GEDTM30Fetcher},
        -110: {'name': 'swot', 'fmts': ['swot'], 'description': 'SWOT Fetcher', 'call': fetchers.SWOTFetcher},
        -111: {'name': 'icesat2', 'fmts': ['icesat2'], 'description': 'IceSat2 Fetcher', 'call': fetchers.IceSat2Fetcher},
        
        -200: {'name': 'charts', 'fmts': ['charts'], 'description': 'Charts Fetcher', 'call': fetchers.ChartsFetcher},
        -201: {'name': 'multibeam', 'fmts': ['multibeam'], 'description': 'Multibeam Fetcher', 'call': fetchers.MBSFetcher},
        -202: {'name': 'hydronos', 'fmts': ['hydronos'], 'description': 'HydroNOS Fetcher', 'call': fetchers.HydroNOSFetcher},
        -203: {'name': 'ehydro', 'fmts': ['ehydro'], 'description': 'eHydro Fetcher', 'call': fetchers.eHydroFetcher},
        -204: {'name': 'bluetopo', 'fmts': ['bluetopo'], 'description': 'BlueTopo Fetcher', 'call': fetchers.BlueTopoFetcher},
        -205: {'name': 'ngs', 'fmts': ['ngs'], 'description': 'NGS Fetcher', 'call': fetchers.NGSFetcher},
        -206: {'name': 'tides', 'fmts': ['tides'], 'description': 'Tides Fetcher', 'call': fetchers.TidesFetcher},
        -207: {'name': 'digital_coast', 'fmts': ['digital_coast'], 'description': 'Digital Coast Fetcher', 'call': fetchers.Fetcher},
        -208: {'name': 'ncei_thredds', 'fmts': ['ncei_thredds'], 'description': 'NCEI Thredds Fetcher', 'call': fetchers.Fetcher},
        -209: {'name': 'tnm', 'fmts': ['tnm'], 'description': 'The National Map Fetcher', 'call': fetchers.Fetcher},
        
        -210: {'name': 'CUDEM', 'fmts': ['CUDEM'], 'description': 'CUDEM Fetcher', 'call': fetchers.Fetcher},
        -211: {'name': 'CoNED', 'fmts': ['CoNED'], 'description': 'CoNED Fetcher', 'call': fetchers.DAVFetcher_CoNED},
        -212: {'name': 'SLR', 'fmts': ['SLR'], 'description': 'SLR Fetcher', 'call': fetchers.DAVFetcher_SLR},
        -213: {'name': 'waterservices', 'fmts': ['waterservices'], 'description': 'Water Services Fetcher', 'call': fetchers.WaterServicesFetcher},
        -214: {'name': 'ned', 'fmts': ['ned', 'ned1'], 'description': 'NED Fetcher', 'call': fetchers.NEDFetcher},
        -215: {'name': 'csb', 'fmts': ['csb'], 'description': 'CSB Fetcher', 'call': fetchers.Fetcher},
        -216: {'name': 'wa_dnr', 'fmts': ['wa_dnr'], 'description': 'WA DNR Fetcher', 'call': fetchers.DNRFetcher},
        -217: {'name': 'r2r', 'fmts': ['r2r'], 'description': 'R2R Fetcher', 'call': fetchers.R2RFetcher},
        
        -300: {'name': 'emodnet', 'fmts': ['emodnet'], 'description': 'EMODNet Fetcher', 'call': fetchers.EMODNetFetcher},
        -301: {'name': 'chs', 'fmts': ['chs'], 'description': 'CHS Fetcher', 'call': fetchers.Fetcher}, 
        -302: {'name': 'hrdem', 'fmts': ['hrdem'], 'description': 'HRDEM Fetcher', 'call': fetchers.HRDEMFetcher},
        -303: {'name': 'arcticdem', 'fmts': ['arcticdem'], 'description': 'ArcticDEM Fetcher', 'call': fetchers.Fetcher},
        -304: {'name': 'mrdem', 'fmts': ['mrdem'], 'description': 'MRDEM Fetcher', 'call': fetchers.Fetcher},
        
        -500: {'name': 'vdatum', 'fmts': ['vdatum'], 'description': 'VDatum Fetcher', 'call': fetchers.VDatumFetcher},
        -600: {'name': 'hydrolakes', 'fmts': ['hydrolakes'], 'description': 'HydroLAKES Fetcher', 'call': fetchers.HydroLakesFetcher},

        ## --- Data Files ---
        167: {'name': 'yxz', 'fmts': ['yxz'], 'description': 'ASCII y,x,z', 'call': xyzfile.YXZFile},
        168: {'name': 'xyz', 'fmts': ['xyz', 'csv', 'dat', 'ascii', 'txt', 'XYZ'], 'description': 'ASCII x,y,z', 'call': xyzfile.XYZFile},
        200: {'name': 'gdal', 'fmts': ['tif', 'tiff', 'img', 'grd', 'nc', 'vrt'], 'description': 'GDAL Raster', 'call': gdalfile.GDALFile},
        201: {'name': 'bag', 'fmts': ['bag'], 'description': 'BAG Bathymetry', 'call': gdalfile.BAGFile},
        202: {'name': 'swot_pixc', 'fmts': ['h5'], 'description': 'SWOT PIXC HDF5', 'call': swotfile.SWOT_PIXC},
        203: {'name': 'swot_hr_raster', 'fmts': ['nc'], 'description': 'SWOT HR Raster', 'call': swotfile.SWOT_HR_Raster},
        300: {'name': 'las', 'fmts': ['las', 'laz'], 'description': 'LAS/LAZ Lidar', 'call': lasfile.LASFile},
        301: {'name': 'mbs', 'fmts': ['fbt', 'mb'], 'description': 'MB-System', 'call': mbsfile.MBSParser},
        302: {'name': 'ogr', 'fmts': ['000', 'shp', 'geojson', 'gpkg', 'gdb'], 'description': 'OGR Vector', 'call': ogrfile.OGRFile},
        303: {'name': 'icesat2_atl03', 'fmts': ['h5'], 'description': 'IceSat2 ATL03', 'call': icesat2file.IceSat2_ATL03},
        304: {'name': 'icesat2_atl24', 'fmts': ['h5'], 'description': 'IceSat2 ATL24', 'call': icesat2file.IceSat2_ATL24},
        310: {'name': 'globato', 'fmts': ['csg', 'nc', 'h5'], 'description': 'CUDEM Globato Block/Stack file.', 'call': CUDEMFile},
        # 310: {'name': 'cudem', 'fmts': ['csg', 'nc', 'h5'], 'description': 'CUDEM NetCDF/H5', 'call': CUDEMFile},
    }

    _datalist_cols = ['path', 'format', 'weight', 'uncertainty', 'title', 'source',
                      'date', 'type', 'resolution', 'horz', 'vert', 'url']

    _metadata_keys = ['name', 'title', 'source', 'date', 'data_type', 'resolution',
                      'hdatum', 'vdatum', 'url']

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        
    def _parse_mod(self, mod=None):
        """Parse the datalist entry line and configure the factory arguments."""
        
        self.kwargs['fn'] = mod
        if self.kwargs['fn'] is None:
            return self
        
        ## Use shlex with posix=False to PRESERVE internal quotes (e.g. pnt_fltrs="rq:...")
        try:
            if os.path.exists(self.kwargs['fn']):
                ## If the entire string is a valid file, treat it as just a filename
                ## check extension or isfile to distinguish from directory unless strictly file
                if os.path.isfile(self.kwargs['fn']):
                    this_entry = [self.kwargs['fn']]
                else:
                    ## It's a directory or strict path
                    this_entry = [self.kwargs['fn']]
            else:
                ## posix=False preserves quotes, which is critical for nested factory strings
                this_entry = shlex.split(self.kwargs['fn'], posix=False)
        ## Fallback to regex if shlex fails
        except Exception:
            this_entry = re.findall(r"(?:\".*?\"|\S)+", self.kwargs['fn'].rstrip())

        utils.echo_debug_msg(f'Initial split entry: {this_entry}')
        ## Convert tokens to appropriate types
        ## Index mapping: 0:fn, 1:fmt, 2:weight, 3:unc, 4+:metadata
        try:
            entry = []
            for n, x in enumerate(this_entry):
                if n == 0:
                    ## Manually strip quotes from filename if they exist (since posix=False kept them)
                    val = utils.str_or(x)
                    if isinstance(val, str) and val.startswith('"') and val.endswith('"'):
                        val = val[1:-1]
                    entry.append(val)
                    #entry.append(utils.str_or(x)) # Filename
                elif n == 1:
                    entry.append(utils.str_or(x, replace_quote=False)) # Format
                elif n == 2:
                    entry.append(utils.float_or(x)) # Weight
                elif n == 3:
                    entry.append(utils.float_or(x)) # Uncertainty
                else:
                    entry.append(utils.str_or(x)) # Metadata
            
            utils.echo_debug_msg(f'Initial parsed entry: {entry}')
        except Exception as e:
            utils.echo_error_msg(f'Could not parse entry {self.kwargs["fn"]}: {e}')
            return self

        ## Determine Data Format
        if len(entry) < 2:
            if self.kwargs.get('data_format') is not None:
                entry.append(self.kwargs.get('data_format'))
            else:
                entry = self.guess_and_insert_fmt(entry)

            if len(entry) < 2:
                utils.echo_error_msg(f'Could not determine format for entry {self.kwargs["fn"]}')
                return self

        elif entry[1] is None or entry[1] == '-':
            if self.kwargs.get('data_format') is not None:
                entry[1] = self.kwargs['data_format']
            else:
                entry = self.guess_and_insert_fmt(entry)

        ## Parse Format Options (e.g., 168:skip=1)
        ## This handles strings like "auto:skip=1" or "-106:pnt_fltrs=..."
        opts = factory.fmod2dict(str(entry[1]), {})
        utils.echo_debug_msg(f'Initial factory opts: {opts}')

        parsed_fmt_id = None
        if '_module' in opts:
            parsed_fmt_id = utils.int_or(opts['_module'])
            self.mod_args = {i: opts[i] for i in opts if i != '_module'}
            entry[1] = parsed_fmt_id

        ## Handle "Auto" (ID 0)
        ## If the format was "auto" (0), we now have the options in self.mod_args, 
        ## but we still need the *real* format ID (e.g. 168 for xyz) to load the right class.
        if parsed_fmt_id == 0:
            guessed_id = self.guess_data_format(entry[0])
            if guessed_id is not None:
                entry[1] = guessed_id
            else:
                utils.echo_warning_msg(f"Could not auto-detect format for {entry[0]}, defaulting to 168 (XYZ)")
                entry[1] = 168
            
        ## Validate Format ID
        try:
            fmt_id = utils.int_or(entry[1])
            assert isinstance(fmt_id, int)
        except AssertionError:
            utils.echo_error_msg(f'Invalid data format ID in entry: {entry}')
            return self

        ## Special Case: Cloud Optimized GeoTIFF via HTTP
        if str(entry[0]).startswith('http') and fmt_id == 200:
            if not str(entry[0]).startswith('/vsicurl/'):
                entry[0] = f'/vsicurl/{entry[0]}'
        
        self.kwargs['data_format'] = fmt_id
        self.mod_name = fmt_id
        
        ## Resolve File Path (Relative vs Absolute)
        if 'parent' not in self.kwargs:
            self.kwargs['parent'] = None

        if self.kwargs['parent'] is None:
            self.kwargs['fn'] = entry[0]
        else:
            ## If parent exists, dataset is not a fetcher (<-2), 
            ## and path is relative, join with parent dir.
            utils.echo_debug_msg(f'Entry ({entry}) has has a parent, adjusting the filename {self.kwargs["fn"]}')
            parent_fmt = self.kwargs['parent'].data_format
            is_fetcher = self.mod_name < -2
            is_absolute = os.path.isabs(entry[0]) or ':' in entry[0] # ':' check for windows drive or url
            
            if not is_fetcher and not is_absolute and parent_fmt >= -2:
                parent_dir = os.path.dirname(self.kwargs['parent'].fn)
                ## Check if already in the same dir to avoid redundancy
                if parent_dir != os.path.dirname(entry[0]): 
                    self.kwargs['fn'] = os.path.join(parent_dir, entry[0])
                else:
                    self.kwargs['fn'] = os.path.relpath(entry[0]) 
            else:
                self.kwargs['fn'] = entry[0]

            utils.echo_debug_msg(f'Adjusted the filename to {self.kwargs["fn"]}')

        ## Set Weight (Inherit/Multiply from Parent)
        if len(entry) < 3:
            entry.append(self.set_default_weight())
        elif entry[2] is None:
            entry[2] = self.set_default_weight()

        ## Initialize base weight
        if 'weight' not in self.kwargs:
            self.kwargs['weight'] = 1

        if self.kwargs['parent'] is not None:
            ## Multiply parent weight by this entry's weight
            if self.kwargs['weight'] is not None:
                self.kwargs['weight'] *= entry[2]
        else:
            if self.kwargs['weight'] is not None:
                self.kwargs['weight'] = entry[2]

        ## Set Uncertainty (Inherit/RSS from Parent)
        if len(entry) < 4:
            entry.append(self.set_default_uncertainty())
        elif entry[3] is None:
            entry[3] = self.set_default_uncertainty()

        if 'uncertainty' not in self.kwargs:
            self.kwargs['uncertainty'] = 0
            
        if self.kwargs['parent'] is not None:
            ## Root Sum of Squares for uncertainty propagation
            if self.kwargs['uncertainty'] is not None:
                self.kwargs['uncertainty'] = math.sqrt(
                    self.kwargs['uncertainty']**2 + entry[3]**2
                )
        else:
            if self.kwargs['uncertainty'] is not None:
                self.kwargs['uncertainty'] = entry[3]
                    
        ## Set Metadata
        if 'metadata' not in self.kwargs:
            self.kwargs['metadata'] = {}

        ## Initialize keys
        for key in self._metadata_keys:
            if key not in self.kwargs['metadata']:
                self.kwargs['metadata'][key] = None

        ## Inherit metadata from parent
        if self.kwargs['parent'] is not None:
            for key in self._metadata_keys:
                if key in self.kwargs['parent'].metadata:
                    if self.kwargs['metadata'][key] is None or key == 'name':
                        self.kwargs['metadata'][key] = self.kwargs['parent'].metadata[key]

        ## Apply specific metadata from entry line
        for i, key in enumerate(self._metadata_keys):
            if key == 'name':
                ## Name defaults to filename base if not present
                if self.kwargs['metadata'][key] is None:
                    ## Remove potential options like :skip=1 for basename
                    clean_fn = self.kwargs['fn'].split(':')[0]
                    self.kwargs['metadata'][key] = utils.fn_basename2(os.path.basename(clean_fn))
                else:
                    ## If name exists (from parent), append current filename
                    clean_fn = self.kwargs['fn'].split(':')[0]
                    self.set_metadata_entry(
                        utils.fn_basename2(os.path.basename(clean_fn)),
                        key, '/'
                    )
            else:
                ## Parsing extra columns in datalist as metadata
                ## entry index offset: 0=fn, 1=fmt, 2=w, 3=unc, 4=start of metadata
                ## since the above key 'name' is not part of the datalist entry,
                ## we set the idx to i+3 since i will be one ahead.
                entry_idx = i + 3
                if len(entry) < entry_idx + 1:
                    entry.append(self.kwargs['metadata'][key])

                self.set_metadata_entry(entry[entry_idx], key, ', ')

                if key == 'date':
                    ## Normalize date ranges
                    self.kwargs['metadata'][key] = utils.num_strings_to_range(
                        self.kwargs['metadata'][key], entry[entry_idx]
                    )
                
        return self.mod_name, self.mod_args
    
    
    def set_metadata_entry(self, entry, metadata_field, join_string='/'):
        """Safely append or set a metadata field."""
        if entry is not None and str(entry) != '-' and str(entry).lower() != 'none':
            current_val = self.kwargs['metadata'].get(metadata_field)
            
            if current_val is not None:
                ## Don't duplicate if identical
                if str(current_val).lower() != str(entry).lower():
                    self.kwargs['metadata'][metadata_field] = join_string.join(
                        [str(current_val), str(entry)]
                    )
            else:
                self.kwargs['metadata'][metadata_field] = str(entry)

                
    def set_default_weight(self):
        return 1

    
    def set_default_uncertainty(self):
        return 0

    
    def guess_data_format(self, fn):
        """Guess a data format ID based on the file extension."""
        
        if not fn: return None
        
        # Check explicit extension mapping
        ext = fn.split('.')[-1]
        for key, info in self._modules.items():
            if 'fmts' in info and ext in info['fmts']:
                return key
                
        ## Hack for .mbXX files
        if len(fn.split('.')) > 1:
            mb_ext = fn.split('.')[-1][:2]
            for key, info in self._modules.items():
                if 'fmts' in info and mb_ext in info['fmts']:
                    return key
        return None


    def guess_and_insert_fmt(self, entry):
        """Analyze the filename in entry[0] and insert the guessed format ID into entry[1]."""
        
        if len(entry) < 1:
            return entry

        ## Check for stdin
        if entry[0] is None or entry[0] == '-':
            if len(entry) > 1:
                entry[1] = 168 # XYZ
            else:
                entry.append(168)
            return entry

        ## Analyze string
        fname = entry[0]
        ext = None
        
        ## Check if file is remote
        is_remote = fname.startswith(('http', 'ftp', '/vsicurl/'))
        
        if is_remote:
            ## Strip query parameters (e.g., data.tif?token=...) for extension check
            clean_fname = fname.split('?')[0].split('&')[0]
            lower_fname = clean_fname.lower()
            
            ## COPC / LAS (Remote) -> LASFile (300)
            if lower_fname.endswith(('.copc.laz', '.laz', '.las')):
                if len(entry) > 1: entry[1] = 300
                else: entry.append(300)
                return entry
            
            ## GDAL / COG (Remote) -> GDALFile (200)
            ## These will be auto-prefixed with /vsicurl/ in _parse_mod if needed.
            elif lower_fname.endswith(('.tif', '.tiff', '.vrt', '.nc', '.h5', '.bag')):
                if len(entry) > 1: entry[1] = 200
                else: entry.append(200)
                return entry
                
            ## Explicit /vsicurl/ -> GDALFile (200)
            elif fname.startswith('/vsicurl/'):
                if len(entry) > 1: entry[1] = 200
                else: entry.append(200)
                return entry
            
            ## Default Fallback -> Generic Fetcher (-100)
            else:
                ext = 'https'
        else:
            ## Local File Logic
            parts = fname.split('.')
            if len(parts) > 1:
                ext = parts[-1]
            else:
                ext = fname.split(':')[0] # Probably a fetches module

        ## Match extension to module ID in _modules registry
        if ext:
            for key, info in self._modules.items():
                if 'fmts' in info and ext in info['fmts']:
                    if len(entry) > 1:
                        entry[1] = int(key)
                    else:
                        entry.append(int(key))
                    break
        
        return entry
    
    
    def write_parameter_file(self, param_file: str):
        """Dump the factory configuration to a JSON file."""
        
        try:
            with open(param_file, 'w') as outfile:
                ## We can't serialize the 'call' objects (classes), so filter them out or handle distinct dict
                ## For now dumping __dict__ excluding unpicklable parts if necessary
                ## This might need refinement based on what exactly needs to be saved
                json.dump(self.__dict__, outfile, default=lambda o: str(o))
                utils.echo_msg(f'New DatasetFactory file written to {param_file}')
                
        except Exception as e:
            raise ValueError(f'DatasetFactory: Unable to write parameter file to {param_file}: {e}')
        

## ==============================================
## Command-line Interface (CLI)
## $ dlim
##
## datalists cli
## ==============================================
class StackModeChoices(list):
    def __contains__(self, item):
        matches = [choice for choice in self if item.split(':')[0] in choice]
        return len(matches) == 1  # Only allow if it's a unique match
    

def datalists_cli():
    """Run datalists from command-line using argparse."""

    parser = argparse.ArgumentParser(
        description=f"DataLists IMproved %(prog)s v{__version__}: Process and generate datalists.",
        epilog="Examples:\n"
        "  dlim -R-90/-89/30/31/-100/100 *.tif -l -w > tifs_in_region.datalist\n"
        "  dlim tifs_in_region.datalist -R -90/-89/30/31 -E 1s > tifs_1s.xyz\n"
        "  dlim -R my_region.shp my_data.xyz -w -J epsg:4326 -P epsg:3565 > my_data_3565.xyz\n"
        "  dlim -R my_region.shp -w multibeam --archive\n"
        "\nSupported %(prog)s modules (see %(prog)s --modules <module-name> for more info):\n" 
        f"{factory.get_module_name_short_desc(DatasetFactory._modules)}\n\n"
        "CUDEM home page: <http://cudem.colorado.edu>",
        formatter_class=argparse.RawTextHelpFormatter
    )

    ## --- Positional Arguments ---
    parser.add_argument(
        'data', 
        nargs='+', 
        default=[], 
        help="Input datalist(s), file(s), or fetch module name(s)."
    )

    ## --- Processing Options ---
    proc_grp = parser.add_argument_group("Processing Options")
    proc_grp.add_argument(
        '-R', '--region', '--aoi', 
        action='append',
        help=regions.region_help_msg()
    )
    proc_grp.add_argument(
        '-E', '--increment', 
        help="Block data to INCREMENT in native units (x-inc[/y-inc])."
    )
    proc_grp.add_argument(
        '-X', '--extend', 
        default='0',
        help="Extend the output DEM region (cells) and processing region (pct)."
    )
    proc_grp.add_argument(
        '-J', '--s_srs', '--src-srs', 
        help="Set the SOURCE projection (EPSG:XXXX or WKT)."
    )
    proc_grp.add_argument(
        '-P', '--t_srs', '--dst-srs', 
        help="Set the TARGET projection (EPSG:XXXX or WKT)."
    )
    proc_grp.add_argument(
        '-A', '--stack-mode', 
        default='mean',
        choices=StackModeChoices(['mean', 'min', 'max', 'mixed', 'supercede']),
        help="Set the STACK MODE (with -E and -R). Default: mean."
    )
    proc_grp.add_argument(
        '-n', '--stack-node', 
        action='store_true', 
        help="Output stacked x/y data as node (center) rather than pixel."
    )
    proc_grp.add_argument(
        '-Z', '--z-precision', 
        type=int, 
        default=4, 
        help="Set the target precision of dumped z values. Default: 4."
    )

    ## --- Filtering Options ---
    filter_grp = parser.add_argument_group("Filtering Options")
    filter_grp.add_argument(
        '-T', '--stack-filter', 
        action='append', 
        default=[],
        help="FILTER the data stack (e.g. outlier:std=2). See `grits --modules`."
    )
    filter_grp.add_argument(
        '-F', '--point-filter', 
        action='append', 
        default=[],
        help="FILTER the POINT data (e.g. outlierz:std=2). See `--point-filters`."
    )
    filter_grp.add_argument(
        '-m', '--mask', 
        action='store_true', 
        help="MASK the datalist to the given REGION/INCREMENTs."
    )
    filter_grp.add_argument(
        '--invert-region', '-v', 
        action='store_true', 
        help="Invert the given region."
    )

    ## --- Output Options (Mutually Exclusive) ---
    out_grp = parser.add_mutually_exclusive_group()
    out_grp.add_argument(
        '-l', '--list', 
        action='store_true', 
        help="List the associated datasets from the datalist."
    )
    out_grp.add_argument(
        '-i', '--info', 
        action='store_true', 
        help="Generate and return an INFO dictionary of the dataset."
    )
    out_grp.add_argument(
        '-r', '--region-inf', 
        action='store_true', 
        help="Output the calculated region of the data."
    )
    out_grp.add_argument(
        '-g', '--glob', 
        action='store_true', 
        help="GLOB the datasets in the current directory to stdout."
    )
    out_grp.add_argument(
        '-V', '--archive', 
        nargs='?', 
        const='archive', 
        metavar='DIR',
        help="Archive the DATALIST to the given REGION. Optional dirname."
    )
    out_grp.add_argument(
        '-G', '--globato', 
        nargs='?', 
        const='globato', 
        help="Archive the DATALIST to a GLOBATO file. Optional filename."
    )
    out_grp.add_argument(
        '-s', '--separate', 
        action='store_true', 
        help="Process and dump each dataset entry independently."
    )

    ## --- Metadata Flags ---
    meta_grp = parser.add_argument_group("Metadata Flags")
    meta_grp.add_argument(
        '-w', '--weights', 
        action='store_true', 
        help="Output WEIGHT values along with xyz."
    )
    meta_grp.add_argument(
        '-u', '--uncertainties', 
        action='store_true', 
        help="Output UNCERTAINTY values along with xyz."
    )
    meta_grp.add_argument(
        '-sm', '--spatial-metadata', 
        action='store_true', 
        help="Generate SPATIAL METADATA (footprints)."
    )

    ## --- System Options ---
    sys_grp = parser.add_argument_group("System Options")
    sys_grp.add_argument(
        '-D', '--cache-dir', 
        default=utils.cudem_cache(), 
        help="CACHE Directory for storing temp and output data."
    )
    sys_grp.add_argument(
        '-q', '--quiet', 
        action='store_true', 
        help="Lower the verbosity to quiet."
    )
    sys_grp.add_argument(
        '--version', 
        action='version', 
        version=f'CUDEM {__cudem_version__} :: %(prog)s {__version__}'
    )

    ## --- Info Helpers (Exit after printing) ---
    sys_grp.add_argument(
        '--modules', 
        nargs='?', 
        const='all', 
        help="Display dataset descriptions. Optional ID for specific details."
    )
    sys_grp.add_argument(
        '--point-filters', 
        nargs='?', 
        const='all', 
        help="Display point filter descriptions."
    )

    ## Parse Arguments
    args = parser.parse_args()

    ## --- Handle Info Helpers ---
    if args.modules:
        mod_key = None if args.modules == 'all' else utils.int_or(args.modules, str(args.modules))
        factory.echo_modules(DatasetFactory._modules, mod_key)
        sys.exit(0)
    
    if args.point_filters:
        mod_key = None if args.point_filters == 'all' else utils.int_or(args.point_filters, str(args.point_filters))
        factory.echo_modules(pointz.PointFilterFactory._modules, mod_key)
        sys.exit(0)

    if args.glob:
        import glob
        ## Flatten format list
        for key, mod in DatasetFactory._modules.items():
            if key != -1 and key != '_factory':
                for fmt in mod.get('fmts', []):
                    for g in glob.glob(f'*.{fmt}'):
                        print(f"{g} {key} 1 0")
        sys.exit(0)

    ## --- Process Input Data ---
    #dls = args.data if args.data else [sys.stdin]
    dls = args.data if args.data else []
    
    ## Handle Increment
    if args.increment:
        xy_inc = args.increment.split('/')
        if len(xy_inc) < 2: xy_inc.append(xy_inc[0])
    else:
        xy_inc = [None, None]

    ## Handle Extend
    ## Supports "-X6" or "-X6:10"
    extend_cells = 0
    extend_pct = 0
    if args.extend:
        parts = str(args.extend).split(':')
        extend_cells = utils.int_or(parts[0], 0)

    ## Handle Regions
    these_regions = regions.parse_cli_region(args.region, not args.quiet) if args.region else [None]

    ## --- Main Loop over Regions ---
    for this_region in these_regions:
        
        ## Buffer Region if Increment is set
        if xy_inc[0] is not None and this_region is not None and extend_cells != 0:
            this_region.buffer(
                x_bv=(utils.str2inc(xy_inc[0]) * extend_cells),
                y_bv=(utils.str2inc(xy_inc[1]) * extend_cells)
            )

        if not dls:
            utils.echo_error_msg('You must specify some type of data.')
            sys.exit(1)

        ## Initialize Datalist
        try:
            this_datalist = init_data(
                dls,
                src_region=this_region,
                src_srs=args.s_srs,
                dst_srs=args.t_srs,
                x_inc=xy_inc[0],
                y_inc=xy_inc[1],
                sample_alg='auto',
                want_weight=args.weights,
                want_uncertainty=args.uncertainties,
                want_verbose=not args.quiet,
                want_mask=args.mask,
                want_sm=args.spatial_metadata,
                invert_region=args.invert_region,
                cache_dir=args.cache_dir,
                dump_precision=args.z_precision,
                pnt_fltrs=args.point_filter,
                stack_fltrs=args.stack_filter,
                stack_node=args.stack_node,
                stack_mode=args.stack_mode
            )

            ## Validate
            if this_datalist is not None and this_datalist.valid_p():
                this_datalist.initialize()
                
                ## --- Actions ---
                if not args.weights: this_datalist.weight = None
                if not args.uncertainties: this_datalist.uncertainty = None

                if args.info:
                    print(this_datalist.inf())
                
                elif args.list:
                    this_datalist.echo()
                
                elif args.region_inf:
                    inf = this_datalist.inf()
                    r = regions.Region().from_list(inf.minmax)
                    
                    ## Warp if requested
                    if args.t_srs:
                        src = args.s_srs if args.s_srs else inf.src_srs
                        if src:
                            r.src_srs = src
                            r.warp(args.t_srs)
                            
                    print(r.format('gmt'))

                elif args.archive:
                    arc_name = args.archive if args.archive != 'archive' else None
                    this_archive = this_datalist.archive_xyz(dirname=arc_name)
                    if this_archive.numpts == 0:
                        utils.remove_glob(f'{this_archive.name}*')

                elif args.globato:
                    this_globato_fn = this_datalist.blocks(out_name=args.globato)
                    utils.echo_msg(f'Generated GLOBATO file: {this_globato_fn}')
                    
                else:
                    ## Default: Dump Data
                    if args.separate:
                        for entry in this_datalist.parse():
                            entry.dump_xyz()
                    else:
                        this_datalist.dump_xyz()

        except KeyboardInterrupt:
            utils.echo_error_msg('Killed by user')
            sys.exit(1)
        except BrokenPipeError:
            ## Standard unix behavior for SIGPIPE (e.g. piping to head)
            sys.stderr.close()
            sys.exit(0)
        except Exception as e:
            utils.echo_error_msg(f"Error: {e}")
            if not args.quiet:
                print(traceback.format_exc())
            sys.exit(1)

if __name__ == "__main__":
    datalists_cli()


### End                      
