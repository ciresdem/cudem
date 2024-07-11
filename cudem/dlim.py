### dlim.py - DataLists IMproved
##
## Copyright (c) 2010 - 2024 Regents of the University of Colorado
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
## A datalist is similar to an MBSystem datalist; it is a space-delineated file containing the following columns:
## data-path data-format data-weight data-uncertainty data-name data-source data-date data-resolution data-type data-horz data-vert data-url
## Minimally, data-path (column 1) is all that is needed.
##
## an associated inf and geojson file will be gerenated for each datalist
## only an associated inf file will be genereated for individual datasets
##
## Parse various dataset types by region/increments and yield data as xyz or array
## recursive data-structures which point to datasets (datalist, zip, fetches, etc) are negative format numbers, e.g. -1 for datalist
##
## supported datasets include: xyz, gdal, ogr, las/laz (laspy), mbs (MBSystem), fetches (see cudem.fetches)
##
## Initialize a datalist/dataset using init_data(list-of-datalist-entries) where the list of datalist entries can
## be any of the supported dataset formats. init_data will combine all the input datafiles into an internal scratch
## datalist and process that.
##
## If region, x_inc, and y_inc are set, all processing will go through Dataset._stacks() where the data will be combined
## either using the 'supercede', 'weighted-mean', 'min' or 'max' method. Dataset._stacks will output a multi-banded gdal file with the following
## bands: 1: z, 2: count, 3: weight, 4: uncerainty, 5: source-uncertainty
##
## If want_mask is set, _stacks() will also generate a multi-band gdal raster file where each
## mask band contains the mask (0/1) of a specific dataset/datalist, depending on the input data. For a datalist, each band will
## contain a mask for each of the top-level datasets. If want_sm is also set, the multi-band mask
## will be processed to an OGR supported vector, with a feature for each band and the metadata items (cols > 4 in the datalist entry) as fields.
##
## Transform data between horizontal/vertical projections/datums by setting src_srs and dst_srs as 'EPSG:<>'
## if src_srs is not set, but dst_srs is, dlim will attempt to obtain the source srs from the data file itself
## or its respective inf file; otherwise, it will be assumed the source data file is in the same srs as dst_srs
##
## A dataset module should have at least an `__init__` and a `yield_ds` method. `yield_ds` should yield a numpy rec array with at least
## 'x', 'y' and 'z' fields, and optionally 'w' and 'u' fields. Other methods that can be re-written include `parse` which yields
## dlim dataset module(s).
##
### Examples:
##
## my_data_list = ['nos_data.xy', 'my_grid.tif', 'gmrt', 'other_data.datalist'] 
## my_processed_datalist = init_data(my_data_list) # initialize the data
## my_processed_datalist.dump_xyz() # dump all the xyz data
##
## from cudem import regions
## my_region = regions.Region().from_list([-160, -150, 30, 40])
## my_data_list = ['nos_data.xy', 'my_grid.tif', 'gmrt', 'other_data.datalist'] 
## my_processed_datalist = init_data(my_data_list, src_region = my_region, x_inc = 1s, y_inc = 1s)
## my_stack = my_processed_datalist._stacks() # stack the data to the region/increments
## my_processed_datalist.archive_xyz() # archive the data to the region/increments
## my_mask = my_processed_datalist._mask() # mask the data to the region/increments
##
### TODO:
## mask to stacks for supercede
## fetch results class
## speed up xyz parsing, esp in yield_array
## temp files
### Code:

import os
import sys
import re
import copy
import json
import math
from tqdm import tqdm
import warnings
import traceback

# import threading
# import multiprocessing as mp
# mp.set_start_method('spawn')
# try:
#    import Queue as queue
# except: import queue as queue
        
import numpy as np
# from scipy.spatial import ConvexHull
# import lxml.etree

import pyproj
import utm
import laspy as lp
from osgeo import gdal
from osgeo import ogr
import h5py as h5

import cudem
from cudem import utils
from cudem import regions
from cudem import xyzfun
from cudem import gdalfun
from cudem import factory
from cudem import vdatums
from cudem import fetches
from cudem import grits
from cudem import vrbag

# cshelph
import pandas as pd
from cudem import cshelph

## Config info and setup
gc = utils.config_check()
gdal.DontUseExceptions()
ogr.DontUseExceptions()
gdal.SetConfigOption('CPL_LOG', 'NUL' if gc['platform'] == 'win32' else '/dev/null') 

## Datalist convenience functions
## data_list is a list of dlim supported datasets
def make_datalist(data_list, want_weight, want_uncertainty, region,
                  src_srs, dst_srs, x_inc, y_inc, sample_alg, verbose):
    """Make a datalist object from a list of supported datasets"""

    ## Make the datalist object
    xdl = Datalist(fn='scratch', data_format=-1, weight=None if not want_weight else 1,
                   uncertainty=None if not want_uncertainty else 0, src_region=region,
                   verbose=verbose, parent=None, src_srs=src_srs, dst_srs=dst_srs,
                   x_inc=x_inc, y_inc=y_inc, sample_alg=sample_alg)

    ## add the datasets from data_list to the datalist object
    xdl.data_entries = [DatasetFactory(fn=" ".join(['-' if x == "" else x for x in dl.split(",")]),
                                       weight=None if not want_weight else 1, uncertainty=None if not want_uncertainty else 0,
                                       src_region=region, verbose=verbose, src_srs=src_srs, dst_srs=dst_srs,
                                       x_inc=x_inc, y_inc=y_inc, sample_alg=sample_alg, parent=xdl,
                                       cache_dir=self.cache_dir)._acquire_module() for dl in data_list]
    return(xdl)

def write_datalist(data_list, outname=None):
    """write a datalist file from a datalist object"""
    
    if outname is None:
        outname = '{}_{}'.format(self.metadata['name'], utils.this_year())
    
    if os.path.exists('{}.datalist'.format(outname)):
        utils.remove_glob('{}.datalist*'.format(outname))
        
    with open('{}.datalist'.format(outname), 'w') as tmp_dl:
        [tmp_dl.write('{}\n'.format(x.format_entry())) for x in data_list]

    return('{}.datalist'.format(outname))

## initialize a list of datasets into a dataset object
def init_data(data_list, region=None, src_srs=None, dst_srs=None, src_geoid=None, dst_geoid='g2018', xy_inc=(None, None),
              sample_alg='auto', want_weight=False, want_uncertainty=False, want_verbose=True, want_mask=False,
              want_sm=False, invert_region=False, cache_dir=None, dump_precision=4, pnt_fltrs=None, stack_fltrs=None, stack_node=True,
              stack_mode='mean'):
    """initialize a datalist object from a list of supported dataset entries"""

    try:
        xdls = [DatasetFactory(mod=" ".join(['-' if x == "" else x for x in dl.split(",")]), data_format = None,
                               weight=None if not want_weight else 1, uncertainty=None if not want_uncertainty else 0,
                               src_srs=src_srs, dst_srs=dst_srs, src_geoid=src_geoid, dst_geoid=dst_geoid, x_inc=xy_inc[0],
                               y_inc=xy_inc[1], sample_alg=sample_alg, parent=None, src_region=region, invert_region=invert_region,
                               cache_dir=cache_dir, want_mask=want_mask, want_sm=want_sm, verbose=want_verbose,
                               dump_precision=dump_precision, pnt_fltrs=pnt_fltrs, stack_fltrs=stack_fltrs, stack_node=stack_node, stack_mode=stack_mode
                               )._acquire_module() for dl in data_list]

        if len(xdls) > 1:
            this_datalist = Scratch(fn=xdls, data_format=-3, weight=None if not want_weight else 1,
                                    uncertainty=None if not want_uncertainty else 0, src_srs=src_srs, dst_srs=dst_srs,
                                    src_geoid=src_geoid, dst_geoid=dst_geoid, x_inc=xy_inc[0], y_inc=xy_inc[1],
                                    sample_alg=sample_alg, parent=None, src_region=region, invert_region=invert_region,
                                    cache_dir=cache_dir, want_mask=want_mask, want_sm=want_sm, verbose=want_verbose,
                                    dump_precision=dump_precision, pnt_fltrs=pnt_fltrs, stack_fltrs=stack_fltrs, stack_node=stack_node,
                                    stack_mode=stack_mode)
        else:
            this_datalist = xdls[0]

        return(this_datalist)
    
    except Exception as e:
        utils.echo_error_msg('could not initialize data, {}: {}'.format(data_list, e))
        return(None)

class PointFilter:
    def __init__(self, points = None, params = {}, **kwargs):
        self.points = points
        self.params = params
        self.kwargs = kwargs

    def __call__(self):
        utils.echo_msg('filtering points using {}'.format(self))
        return(self.run())

    def run(self):
        raise(NotImplementedError)
    
    def convert_wgs_to_utm(self, lat, lon):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            easting, northing, num, letter = utm.from_latlon(lat, lon)
            if letter >= 'N':
                epsg = 'epsg:326' + str(num)
            elif letter < 'N':
                epsg = 'epsg:327' + str(num)
            else:
                utils.echo_error_msg('Could not find UTM zone')

            return(epsg)

class BinZ(PointFilter):
    def __init__(self, y_res=1, z_res=.5, percentile=50, z_min=None, z_max=None, **kwargs):
        super().__init__(**kwargs)
        self.y_res = utils.float_or(y_res, 1)
        self.z_res = utils.float_or(z_res, .5)
        self.percentile = utils.float_or(percentile, 50)
        self.z_min = utils.float_or(z_min)
        self.z_max = utils.float_or(z_max)
        
    def bin_points(self, points):
        '''Bin data along vertical and horizontal scales for later segmentation'''

        ## Calculate number of bins required both vertically and
        ## horizontally with resolution size
        y_bin_number = round(abs(points['y'].min() - points['y'].max())/self.y_res)
        if y_bin_number == 0:
            y_bin_number = 1
            
        #x_bin_number = round(abs(points['x'].min() - points['x'].max())/self.y_res)
        z_bin_number = round(abs(points['z'].min() - points['z'].max())/self.z_res)
        if z_bin_number == 0:
            z_bin_number = 1

        if (y_bin_number > 0 and z_bin_number > 0):    
            points1 = points
            y_bins = pd.cut(points['y'], y_bin_number, labels = np.array(range(y_bin_number)))
            #x_bins = pd.cut(points['x'], x_bin_number, labels = np.array(range(x_bin_number)))
            points1['y_bins'] = y_bins
            #points1['x_bins'] = x_bins
            #utils.echo_msg('{} {} {}'.format(points['z'].min(), points['z'].max(), z_bin_number))
            z_bins = pd.cut(
                points['z'], z_bin_number, labels = np.round(
                    np.linspace(points['z'].min(), points['z'].max(), num=z_bin_number),
                    decimals = 1
                )
            )
            points1['z_bins'] = z_bins
            points1 = points1.reset_index(drop=True)

            return(points1)

        return(None)

    def get_bin_height(self, binned_data):
        '''Calculate mean sea height for easier calculation of depth and cleaner figures'''

        # Create sea height list
        bin_height = []
        bin_lat = []
        bin_lon = []

        # Group data by latitude
        binned_data_sea = binned_data
        #grouped_data = binned_data_sea.groupby(['y_bins', 'x_bins'], group_keys=True)
        grouped_data = binned_data_sea.groupby(['y_bins'], group_keys=True)
        data_groups = dict(list(grouped_data))

        # Create a percentile threshold of photon counts in each grid, grouped by both x and y axes.
        #count_threshold = np.percentile(binned_data.groupby(['y_bins', 'x_bins', 'z_bins']).size().reset_index().groupby('y_bins')[[0]].max(), percentile)
        count_threshold = np.percentile(binned_data.groupby(['y_bins', 'z_bins']).size().reset_index().groupby('y_bins')[[0]].max(), self.percentile)
        #utils.echo_msg(count_threshold)
        # Loop through groups and return average sea height
        for k,v in data_groups.items():
            # Create new dataframe based on occurance of photons per height bin
            new_df = pd.DataFrame(v.groupby('z_bins').count())

            # Return the bin with the highest count
            largest_h_bin = new_df['z'].argmax()

            # Select the index of the bin with the highest count
            largest_h = new_df.index[largest_h_bin]

            # Set threshold of photon counts per bin
            if new_df.iloc[largest_h_bin]['y'] >= count_threshold:

                [bin_lat.append(x) for x in v.loc[v['z_bins']==largest_h, 'y']]
                [bin_lon.append(x) for x in v.loc[v['z_bins']==largest_h, 'x']]
                [bin_height.append(x) for x in v.loc[v['z_bins']==largest_h, 'z']]
                
                # # Calculate the median value of all values within this bin
                # lat_bin_sea_median = v.loc[v['z_bins']==largest_h, 'z'].median()
                # lat_bin_median = v.loc[v['z_bins']==largest_h, 'y'].median()
                # lon_bin_median = v.loc[v['z_bins']==largest_h, 'x'].median()

                # # Append to sea height list
                # bin_height.append(lat_bin_bin_median)
                # bin_lat.append(lat_bin_median)
                # bin_lon.append(lon_bin_median)
                del new_df
            else:
                del new_df

        # Filter out sea height bin values outside 2 SD of mean.
        if np.all(np.isnan(bin_height)):
            return(None)

        mean = np.nanmean(bin_height, axis=0)
        sd = np.nanstd(bin_height, axis=0)
        bin_height_1 = np.where((bin_height > (mean + 2*sd)) | (bin_height < (mean - 2*sd)), np.nan, bin_height).tolist()
                
        return(bin_lat, bin_lon, bin_height_1)
    
    def run(self):
        try:
            epsg_code = self.convert_wgs_to_utm(self.points['y'][0], self.points['x'][0])
            epsg_num = int(epsg_code.split(':')[-1])
            utm_proj = pyproj.Proj(epsg_code)
            x_utm, y_utm = utm_proj(self.points['x'], self.points['y'])
        except:
            utils.echo_warning_msg('could not transform to utm')
            x_utm, y_utm = self.points['x'], self.points['y']

        points_1 = pd.DataFrame(
            {'y': y_utm,
             'x': x_utm,
             'z': self.points['z']},
            columns=['y', 'x', 'z']
        )

        if utils.float_or(self.z_min) is not None:
            points_1 = points_1[(points_1['z'] > self.z_min)]

        if utils.float_or(self.z_max) is not None:
            points_1 = points_1[(points_1['z'] < self.z_max)]

        #print(points_1)
        #print(len(points_1))
        if len(points_1) > 0:
            binned_points = self.bin_points(points_1)
            points_1 = None
            
            if binned_points is not None:
                ys, xs, zs = self.get_bin_height(binned_points)
                binned_points = None
                
                bin_ds = np.column_stack((xs, ys, zs))
                bin_ds = np.rec.fromrecords(bin_ds, names='x, y, z')
                
                xs = ys = zs = None
                bin_ds = bin_ds[~np.isnan(bin_ds['z'])]
                med_surface_h = np.nanmedian(bin_ds['z'])
                #bin_ds = bin_ds[bin_ds['z'] < med_surface_h + (z_res * 2)]
                #bin_ds = bin_ds[bin_ds['z'] > med_surface_h - (z_res * 2)]

                transformer = pyproj.Transformer.from_crs("EPSG:"+str(epsg_num), "EPSG:4326", always_xy=True)
                lon_wgs84, lat_wgs84 = transformer.transform(bin_ds['x'], bin_ds['y'])
                bin_points = np.column_stack((lon_wgs84, lat_wgs84, bin_ds['z']))
                lon_wgs84 = lat_wgs84 = bin_ds = None
                bin_points = np.rec.fromrecords(bin_points, names='x,y,z')
                #utils.echo_msg('bin_points: {}'.format(bin_points))
                return(bin_points)
            # return(bin_ds)
            
        return(points_1)

class PointFilterFactory(factory.CUDEMFactory):
    _modules = {
        'bin_z': {'name': 'bin_z', 'call': BinZ},
    }
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
    
class INF:
    """INF Files contain information about datasets
    """
    
    def __init__(self, name = None, file_hash = None, numpts = 0, minmax = [], wkt = None,
                 fmt = None, src_srs = None):
        self.name = name
        self.file_hash = file_hash
        self.hash = file_hash
        self.numpts = numpts
        self.minmax = minmax
        self.wkt = wkt
        self.fmt = fmt
        self.format = fmt
        self.src_srs = src_srs

    def __str__(self):
        return('<Dataset Info: {}>'.format(self.__dict__))
    
    def __repr__(self):
        return('<Dataset Info: {}>'.format(self.__dict__))

    def generate_hash(self, fn = None, sha1 = False):
        """generate a hash of the xyz-dataset source file"""

        import hashlib
        fn = self.name if fn is None else fn
        BUF_SIZE = 65536
        if sha1:
            this_hash = hashlib.sha1()
        else:
            this_hash = hashlib.md5()
            
        try:
            with open(fn, 'rb') as f:
                while True:
                    data = f.read(BUF_SIZE)
                    if not data:
                        break

                    this_hash.update(data)

            self.file_hash = this_hash.hexdigest()
        
        except:
            self.file_hash = '0'

        return(self.file_hash)
        
    def generate(self):
        if self.name is None:
            return(self)

        this_ds = DatasetFactory(mod=self.name)
        this_ds_mod = this_ds._acquire_module()        
        self.fmt = this_ds.mod_name            
        
    def load_inf_file(self, inf_fn = None):
        valid_data = False
        data = {}
        if inf_fn is None:
            return(self)
        
        if os.path.exists(inf_fn):
            try:
                with open(inf_fn, 'r') as inf_ob:
                    data = json.load(inf_ob)
            except ValueError:
                try:
                    data = MBSParser(fn=inf_fn).inf_parse().infos.__dict__
                    #self.check_hash = False
                    mb_inf = True
                except Exception as e:
                    raise ValueError('CUDEMFactory: Unable to read data from {} as mb-system inf, {}'.format(inf_fn, e))
            except:
                raise ValueError('CUDEMFactory: Unable to read data from {} as json'.format(inf_fn))

        
        for ky, val in data.items():
            if ky in self.__dict__:
                self.__setattr__(ky, val)
                valid_data = True

        return(self)
            
    def write_inf_file(self, inf_fn = None):
        if inf_fn is None:
            if self.name is not None:
                inf_fn = '{}.inf'.format(self.name)
        try:
            with open(inf_fn, 'w') as outfile:
                json.dump(self.__dict__, outfile)
        except:
            pass
    
class ElevationDataset:
    """representing an Elevation Dataset
    
    This is the super class for all datalist (dlim) datasets .    
    Each dataset sub-class should define a dataset-specific    
    data parser <self.parse> and a <self.generate_inf> function to generate
    inf files. 

    Specifically, each sub-dataset should minimally define the following functions:

    sub_ds.__init__
    sub_ds.yield_ds
    
    Where:
    generate_inf generates a dlim compatible inf file,
    parse yields the dlim dataset module
    yield_ds yields the data as a numpy rec-array with 'x', 'y', 'z' and optionally 'w' and 'u'

    ----
    Parameters:

    fn: dataset filename or fetches module
    data_format: dataset format
    weight: dataset weight
    uncertainty: dataset uncertainty
    mask: mask the dataset
    src_srs: dataset source srs
    dst_srs: dataset target srs
    src_geoid: dataset source geoid (if applicable)
    dst_geoid: dataset target geoid
    x_inc: target dataset x/lat increment
    y_inc: target dataset y/lon increment
    want_mask: generate a mask of the data
    want_sm: generate spatial metadata vector
    sample_alg: the gdal resample algorithm
    parent: dataset parent obj
    region: ROI
    invert_region: invert the region
    stack_node: output x/y data as mean of input x/y instead of cell center
    stack_mode: 'min', 'max', 'supercede', 'mean'
    cache_dir: cache_directory
    verbose: be verbose
    remote: dataset is remote
    dump_precision: specify the float precision of the dumped xyz data
    params: the factory parameters
    metadata: dataset metadata
    """

    gdal_sample_methods = [
        'near', 'bilinear', 'cubic', 'cubicspline', 'lanczos',
        'average', 'mode',  'max', 'min', 'med', 'Q1', 'Q3', 'sum',
        'auto'
    ]

    stack_modes = [
        'min', 'max', 'mean', 'supercede'
    ]

    ## todo: add transformation grid option (stacks += transformation_grid), geoids
    def __init__(self, fn = None, data_format = None, weight = 1, uncertainty = 0, src_srs = None, mask = None,
                 dst_srs = 'epsg:4326', src_geoid = None, dst_geoid = 'g2018', x_inc = None, y_inc = None, want_mask = False,
                 want_sm = False, sample_alg = 'auto', parent = None, src_region = None, invert_region = False,
                 pnt_fltrs=[], stack_fltrs=[], stack_node = True, stack_mode = 'mean', cache_dir = None, verbose = False, remote = False,
                 dump_precision = 6, params = {}, metadata = {'name':None, 'title':None, 'source':None, 'date':None,
                                                              'data_type':None, 'resolution':None, 'hdatum':None,
                                                              'vdatum':None, 'url':None}, **kwargs):
        self.fn = fn # dataset filename or fetches module
        self.data_format = data_format # dataset format
        self.weight = weight # dataset weight
        self.uncertainty = uncertainty # dataset uncertainty
        self.mask = mask #{'mask': mask, 'invert': False}#mask # dataset mask
        self.src_srs = src_srs # dataset source srs
        self.dst_srs = dst_srs # dataset target srs
        self.src_geoid = src_geoid # dataset source geoid
        self.dst_geoid = dst_geoid # dataset target geoid
        self.x_inc = utils.str2inc(x_inc) # target dataset x/lat increment
        self.y_inc = utils.str2inc(y_inc) # target dataset y/lon increment
        self.stack_fltrs = stack_fltrs # pass the stack, if generated, through grits filters
        self.pnt_fltrs = pnt_fltrs # pass the points through filters
        self.want_mask = want_mask # mask the data
        self.want_sm = want_sm # generate spatial metadata vector
        self.sample_alg = sample_alg # the gdal resample algorithm
        self.metadata = copy.deepcopy(metadata) # dataset metadata
        self.parent = parent # dataset parent obj
        self.region = src_region # ROI
        self.invert_region = invert_region # invert the region
        self.cache_dir = cache_dir # cache_directory
        self.verbose = verbose # be verbose
        self.remote = remote # dataset is remote
        self.dump_precision = dump_precision # the precision of the dumped xyz data
        self.stack_node = stack_node # yield avg x/y data instead of center
        self.stack_mode = stack_mode # 'mean', 'min', 'max', 'supercede'
        self.mask_keys = ['mask', 'invert_mask', 'ogr_or_gdal'] # options for input data mask
        if self.mask is not None:
            if isinstance(self.mask, str):
                self.mask = {'mask': self.mask}
                
            for kpam, kval in kwargs.items():
                if kpam not in self.__dict__:
                    self.mask[kpam] = kval

            for kpam, kval in self.mask.items():
                if kpam in kwargs:
                    del kwargs[kpam]

            self.mask['ogr_or_gdal'] = gdalfun.ogr_or_gdal(self.mask['mask'])
            for key in self.mask_keys:
                if key not in self.mask.keys():
                    self.mask[key] = None
            
        self.infos = INF(name=self.fn, file_hash='0', numpts=0, fmt=self.data_format) # infos blob            
        self.params = params # the factory parameters
        if not self.params:
            self.params['kwargs'] = self.__dict__.copy()
            self.params['mod'] = self.fn
            self.params['mod_name'] = self.data_format
            self.params['mod_args'] = {}
            
    def __str__(self):
        return('<Dataset: {} - {}>'.format(self.metadata['name'], self.fn))
    
    def __repr__(self):
        return('<Dataset: {} - {}>'.format(self.metadata['name'], self.fn))

    def __call__(self):
        self.initialize()
        
    def initialize(self):
        self._fn = None # temp filename holder
        self.transformer = None # pyproj Transformer obj
        self.trans_fn = None # vertical datum transformation grid
        self.trans_fn_unc = None # vertical transformation uncertainty
        self.trans_region = None # transformed region
        self.src_proj4 = None # source crs as proj4 string
        self.dst_proj4 = None # destination crs as proj4 string
        self.aux_src_proj4 = None # source horizontal crs as proj4 string
        self.aux_dst_proj4 = None # destination horizontal crs as proj4 string
        self.data_region = None # self.region and inf.region reduced
        self.archive_datalist = None # the datalist of the archived data
        self.data_entries = [] # 
        self.data_lists = {} #
        self.cache_dir = utils.cudem_cache() if self.cache_dir is None else self.cache_dir # cache directory
        if self.sample_alg not in self.gdal_sample_methods: # gdal_warp resmaple algorithm 
            utils.echo_warning_msg(
                '{} is not a valid gdal warp resample algorithm, falling back to auto'.format(
                    self.sample_alg
                )
            )
            self.sample_alg = 'auto'

        if utils.fn_url_p(self.fn):
            self.remote = True

        ## Dataset Mask, if the input mask is a vector, rasterize it if possible.
        if self.mask is not None:
            if self.mask['ogr_or_gdal'] == 1: # mask is ogr, rasterize it
                if self.region is not None and self.x_inc is not None and self.y_inc is not None:
                    self.mask['mask'] = gdalfun.ogr2gdal_mask(
                        self.mask['mask'], region=self.region, x_inc=self.x_inc, y_inc=self.y_inc, dst_srs=self.dst_srs,
                        invert=True, verbose=self.verbose, temp_dir=self.cache_dir
                    )
                    self.mask['ogr_or_gdal'] = 0
                       
        if self.valid_p():
            self.infos = self.inf(check_hash=True if self.data_format == -1 else False)
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                try:
                    self.set_transform()
                except Exception as e:
                    utils.echo_error_msg('could not set transformation, {}'.format(e))
                
            self.set_yield()

        ## initialize filters
        if isinstance(self.stack_fltrs, str):
            self.stack_fltrs = [self.stack_fltrs]

        if isinstance(self.pnt_fltrs, str):
            #self.pnt_fltrs = [self.pnt_fltrs]
            self.pnt_fltrs = [':'.join(self.pnt_fltrs.split('/'))]
            
        # if isinstance(self.fltrs, list):
        #     for f in self.fltrs:
        #         if f.split(':')[0] in grits.GritsFactory()._modules.keys():
        #             self.stack_filters.append(f)
        #         elif f.split(':')[0] in self.dlim_filters:
        #             self.point_filters.append(f)            
        # else:
        #     self.fltrs = []    

        return(self)
    
    def generate_inf(self, callback=lambda: False):
        """generate an inf file for the data source. this is generic and
        will parse through all the data via yield_xyz to get point count
        and region. if the datasource has a better way to do this, such as
        with las/gdal files, re-define this function in its respective
        dataset sub-class.

        todo: hull for wkt
        ## todo: put qhull in a yield and perform while scanning for min/max
        # try:
        #     out_hull = [pts[i] for i in ConvexHull(
        #         pts, qhull_options='Qt'
        #     ).vertices]
        #     out_hull.append(out_hull[0])
        #     self.infos['wkt'] = regions.create_wkt_polygon(out_hull, xpos=0, ypos=1)
        # except:
        #     ...
        """

        _region = None
        if self.region is not None:
            _region = self.region.copy()
            self.region = None

        _x_inc = None
        if self.x_inc is not None:
            _x_inc = self.x_inc
            self.x_inc = None
            
        this_region = regions.Region()
        point_count = 0

        for points in self.yield_points():
            if point_count == 0:
                this_region.from_list(
                    [
                        points['x'].min(), points['x'].max(),
                        points['y'].min(), points['y'].max(),
                        points['z'].min(), points['z'].max()
                    ]
                )
            else:
                if points['x'].min() < this_region.xmin:
                    this_region.xmin = points['x'].min()
                elif points['x'].max() > this_region.xmax:
                    this_region.xmax = points['x'].max()
                    
                if points['y'].min() < this_region.ymin:
                    this_region.ymin = points['y'].min()
                elif points['y'].max() > this_region.ymax:
                    this_region.ymax = points['y'].max()
                    
                if points['z'].min() < this_region.zmin:
                    this_region.zmin = points['z'].min()
                elif points['z'].min() > this_region.zmax:
                    this_region.zmax = points['z'].min()
                
            point_count += len(points)

        self.infos.numpts = point_count
        if point_count > 0:
            self.infos.minmax = this_region.export_as_list(include_z=True)
            self.infos.wkt = this_region.export_as_wkt()
            
        self.infos.src_srs = self.src_srs
        if _region is not None:
            self.region = _region.copy()

        if _x_inc is not None:
            self.x_inc = _x_inc
            
        return(self.infos)

    def yield_ds(self):
        """yield the numpy xyz rec array points from the dataset.

        reset in dataset if needed.
        """

        for ds in self.parse():
            for points in ds.yield_ds():
                yield(points)
        
    def yield_xyz_from_entries(self):
        """yield from self.data_entries, list of datasets"""
        
        for this_entry in self.data_entries:
            for xyz in this_entry.xyz_yield:
                yield(xyz)
                
            if this_entry.remote:
                utils.remove_glob('{}*'.format(this_entry.fn))

    def yield_entries(self):
        """yield from self.data_entries, list of datasets"""
        
        for this_entry in self.data_entries:
            yield(this_entry)

    def set_yield(self):
        """set the yield strategy, either default (all points) or mask or stacks

        This sets both `self.array_yeild` and `self.xyz_yield`.
        
        if region, x_inc and y_inc are set, data will yield through _stacks,
        otherwise, data will yield directly from `self.yield_xyz` and `self.yield_array`
        """

        self.array_yield = self.mask_and_yield_array()
        self.xyz_yield = self.mask_and_yield_xyz()
        #if self.want_archive: # archive only works when yielding xyz data.
        #    self.xyz_yield = self.archive_xyz()
        if self.region is None and self.x_inc is not None:
            utils.echo_warning_msg('must enter a region to output in increment blocks...')

        if self.region is not None and self.x_inc is not None:
            self.x_inc = utils.str2inc(self.x_inc)
            if self.y_inc is None:
                self.y_inc = self.x_inc
            else:
                self.y_inc = utils.str2inc(self.y_inc)

            out_name = utils.make_temp_fn('dlim_stacks', temp_dir=self.cache_dir)
            self.xyz_yield = self.stacks_yield_xyz(out_name=out_name, fmt='GTiff')#, mode=self.stack_mode)
            
    def mask_and_yield_array(self):
        """mask the incoming array from `self.yield_array` and yield the results.
        
        the mask should be either an ogr supported vector or a gdal supported raster.
        """
        
        mask_band = None
        mask_infos = None

        if self.mask is not None:
            if os.path.exists(self.mask['mask']):            
                utils.echo_msg('using mask dataset: {}'.format(self.mask['mask']))
                if self.region is not None and self.x_inc is not None and self.y_inc is not None:
                    #dst_srs = self.dst_srs.split('+geoid')[0]
                    src_mask = gdalfun.sample_warp(
                        self.mask['mask'], None, self.x_inc, self.y_inc,
                        src_region=self.region, sample_alg='nearest', dst_srs=self.dst_srs,
                        ndv=gdalfun.gdal_get_ndv(self.mask['mask']), verbose=self.verbose,
                        co=["COMPRESS=DEFLATE", "TILED=YES"]
                    )[0]
                else:
                    src_mask = gdal.Open(self.mask['mask'])
                    
                mask_band = src_mask.GetRasterBand(1)
                mask_infos = gdalfun.gdal_infos(src_mask)
            else:
                utils.echo_warning_msg('could not load mask {}'.format(self.mask['mask']))

        for out_arrays, this_srcwin, this_gt in self.yield_array():
            if mask_band is not None:
                ycount, xcount = out_arrays['z'].shape
                this_region = regions.Region().from_geo_transform(this_gt, xcount, ycount)
                mask_data = mask_band.ReadAsArray(*this_srcwin)
                if not np.isnan(mask_infos['ndv']):
                    mask_data[mask_data==mask_infos['ndv']] = np.nan
                        
                for arr in out_arrays.keys():
                    if out_arrays[arr] is not None:
                        if self.mask['invert_mask']:
                            out_arrays[arr][~np.isnan(mask_data)] = np.nan
                        else:
                            out_arrays[arr][np.isnan(mask_data)] = np.nan

            yield(out_arrays, this_srcwin, this_gt)

    def mask_and_yield_xyz(self):
        """mask the incoming xyz data from `self.yield_xyz` and yield the results.
        
        the mask should be either an ogr supported vector or a gdal supported raster.
        """

        for this_entry in self.parse():
            if this_entry.mask is None:
                for this_xyz in this_entry.yield_xyz():
                    yield(this_xyz)
            else:
                utils.echo_msg('using mask dataset: {}'.format(this_entry.mask['mask']))
                if this_entry.mask['ogr_or_gdal'] == 0:
                    src_ds = gdal.Open(this_entry.mask['mask'])
                    if src_ds is not None:
                        ds_config = gdalfun.gdal_infos(src_ds)
                        ds_band = src_ds.GetRasterBand(1)
                        ds_gt = ds_config['geoT']
                        ds_nd = ds_config['ndv']
                    
                        for this_xyz in this_entry.yield_xyz():
                            xpos, ypos = utils._geo2pixel(this_xyz.x, this_xyz.y, ds_gt, node='pixel')
                            if xpos < ds_config['nx'] and ypos < ds_config['ny'] and xpos >=0 and ypos >=0:
                                tgrid = ds_band.ReadAsArray(xpos, ypos, 1, 1)
                                if tgrid is not None:
                                    if tgrid[0][0] == ds_nd:
                                        if this_entry.mask['invert_mask']:
                                            yield(this_xyz)
                                    else:
                                        if not this_entry.mask['invert_mask']:
                                            yield(this_xyz)                                            
                else:
                    ## this is very slow! find another way.
                    src_ds = ogr.Open(this_entry.mask['mask'])
                    layer = src_ds.GetLayer()
                    utils.echo_msg(len(layer))

                    geomcol = gdalfun.ogr_union_geom(layer)
                    if not geomcol.IsValid():
                        geomcol = geomcol.Buffer(0)
                    
                    for this_xyz in this_entry.yield_xyz():
                        this_wkt = this_xyz.export_as_wkt()
                        this_point = ogr.CreateGeometryFromWkt(this_wkt)
                        if this_point.Within(geomcol):
                            if not this_entry.mask['invert_mask']:
                                yield(this_xyz)
                        else:
                            if not this_entry.mask['invert_mask']:
                                yield(this_xyz)
            
    def dump_xyz(self, dst_port=sys.stdout, encode=False):
        """dump the XYZ data from the dataset.

        data gets parsed through `self.xyz_yield`. See `set_yield` for more info.
        """

        for this_xyz in self.xyz_yield:
            this_xyz.dump(
                include_w=True if self.weight is not None else False,
                include_u=True if self.uncertainty is not None else False,
                dst_port=dst_port,
                encode=encode,
                precision=self.dump_precision
            )

    def dump_xyz_direct(self, dst_port=sys.stdout, encode=False):
        """dump the XYZ data from the dataset

        data get dumped directly from `self.yield_xyz`, by-passing `self.xyz_yield`.
        """
        
        for this_xyz in self.yield_xyz():
            this_xyz.dump(
                include_w=True if self.weight is not None else False,
                include_u=True if self.uncertainty is not None else False,
                dst_port=dst_port,
                encode=encode,
                precision=self.dump_precision
            )

    def export_xyz_as_list(self, z_only = False):
        """return the XYZ data from the dataset as python list"""

        return([xyz.z if z_only else xyz.copy() for xyz in self.xyz_yield])
            
    def fetch(self):
        """fetch remote data from self.data_entries"""
        
        for entry in self.data_entries:
            if entry.remote:
                if entry._fn is None:
                    entry._fn = os.path.basename(self.fn)
                    
                f = utils.Fetch(
                    url=entry.fn, verbose=entry.verbose
                )
                if f.fetch_file(entry._fn) == 0:
                    entry.fn = entry._fn
            else:
                utils.echo_warning_msg('nothing to fetch')

    def valid_p(self, fmts=['scratch']):
        """check if self appears to be a valid dataset entry"""

        if self.fn is None:
            return(False)

        if self.data_format is None:
            return(False)

        #check_path=True
        #if self.fn.startswith('http') or self.fn.startswith('/vsicurl/') or self.fn.startswith('BAG'):
        #    check_path = False
            
        #if check_path:
        if self.fn is not None:
            if self.fn not in fmts:
                if not isinstance(self.fn, list):
                    if self.fn.startswith('http') or self.fn.startswith('/vsicurl/') or self.fn.startswith('BAG'):
                        if not utils.fn_url_p(self.fn):
                            if self.data_format > -10:
                                if not os.path.exists(self.fn):
                                    return (False)

                                if os.stat(self.fn).st_size == 0:
                                    return(False)
                        
        return(True)
        
    def format_entry(self, sep=' '):
        """format the dataset information as a `sep` separated string."""
        
        dl_entry = sep.join(
            [str(x) for x in [
                self.fn,
                '{}:{}'.format(self.data_format, factory.dict2args(self.params['mod_args'])),
                self.weight,
                self.uncertainty
            ]]
        )
        metadata = self.echo_()
        return(sep.join([dl_entry, metadata]))
        
    def echo_(self, sep=' ', **kwargs):
        """print self as a datalist entry string"""

        out = []
        for key in self.metadata.keys():
            if key != 'name':
                out.append(str(self.metadata[key]))

        return(sep.join(['"{}"'.format(str(x)) for x in out]))
    
    def echo(self, **kwargs):
        """print self.data_entries as a datalist entries."""

        for entry in self.parse():
            entry_path = os.path.abspath(entry.fn) if not self.remote else entry.fn
            l = [entry_path, entry.data_format]
            if entry.weight is not None:
                l.append(entry.weight)
                
            print('{}'.format(" ".join([str(x) for x in l])))

    def format_metadata(self, **kwargs):
        """format metadata from self, for use as a datalist entry."""

        return(self.echo_())
    
    def inf(self, check_hash=False, recursive_check=False, write_inf=True, **kwargs):
        """read/write an inf file

        If the inf file is not found, will attempt to generate one.
        The function `generate_inf` should be defined for each specific
        dataset sub-class.
        """

        inf_path = '{}.inf'.format(self.fn)
        generate_inf = False

        ## try to parse the existing inf file as either a native inf json file
        if os.path.exists(inf_path):
            try:
                self.infos.load_inf_file(inf_path)
                
            except ValueError as e:
                generate_inf = True
                utils.remove_glob(inf_path)
                if self.verbose:
                    utils.echo_warning_msg(
                        'failed to parse inf {}, {}'.format(inf_path, e)
                    )
        else:
            generate_inf = True        

        ## check hash from inf file vs generated hash,
        ## if hashes are different, then generate a new
        ## inf file...only do this if check_hash is set
        ## to True, as this can be time consuming and not
        ## always necessary...
        if check_hash:
            generate_inf = self.infos.generate_hash() != self.infos.file_hash

        ## this being set can break some modules (bags esp)
        # if self.remote:
        #     generate_inf = False

        if generate_inf:
            self.infos = self.generate_inf()
            
            ## update this
            if self.data_format >= -1 and write_inf:
                self.infos.write_inf_file()

            if recursive_check and self.parent is not None:
                self.parent.inf(check_hash=True)

        if self.infos.src_srs is None:
            self.infos.src_srs = self.src_srs

        else:
            if self.src_srs is None:
                self.src_srs = self.infos.src_srs

            if self.transformer is not None and self.trans_region is None:
                self.trans_region = regions.Region().from_list(self.infos.minmax)
                self.trans_region.src_srs = self.infos.src_srs
                self.trans_region.warp(self.dst_srs)

        return(self.infos)

    def set_transform(self):
        """Set the pyproj horizontal and vertical transformations for the dataset"""
        
        want_vertical = True
        src_geoid = None
        dst_geoid = 'g2018'

        if self.src_srs is not None and self.dst_srs is not None:
            tmp_src_srs = self.src_srs.split('+geoid:')
            src_srs = tmp_src_srs[0]
            self.src_srs = src_srs
            if len(tmp_src_srs) > 1:
                src_geoid = tmp_src_srs[1]

            tmp_dst_srs = self.dst_srs.split('+geoid:')
            dst_srs = tmp_dst_srs[0]
            self.dst_srs = dst_srs
            if len(tmp_dst_srs) > 1:
                dst_geoid = tmp_dst_srs[1]
                
            is_esri = False
            in_vertical_epsg_esri = None
            if 'ESRI' in src_srs.upper():
                is_esri = True
                srs_split = src_srs.split('+')
                src_srs = srs_split[0]
                if len(srs_split) > 1:
                    in_vertical_epsg_esri = srs_split[1]
                
            in_crs = pyproj.CRS.from_user_input(src_srs)
            out_crs = pyproj.CRS.from_user_input(dst_srs)

            if in_crs.is_compound:
                in_crs_list = in_crs.sub_crs_list
                in_horizontal_crs = in_crs_list[0]
                in_vertical_crs = in_crs_list[1]
                in_vertical_name = in_vertical_crs.name
                in_vertical_epsg = in_vertical_crs.to_epsg()
                if in_vertical_epsg is None:
                    in_vertical_epsg = in_vertical_name
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

            in_horizontal_epsg = in_horizontal_crs.to_epsg()
            out_horizontal_epsg = out_horizontal_crs.to_epsg()

            if (in_vertical_epsg_esri is not None and is_esri):
                in_vertical_epsg = in_vertical_epsg_esri
                if out_vertical_epsg is not None:
                    want_vertical = True
            
            if want_vertical:
                if (in_vertical_epsg == out_vertical_epsg) and self.src_geoid is None:
                    want_vertical = False
                    
            ## horizontal Transformation
            try:
                self.transformer = pyproj.Transformer.from_crs(in_horizontal_crs, out_horizontal_crs, always_xy=True)
            except Exception as e:
                utils.echo_warning_msg('could not set transformation in: {}, out: {}, {}'.format(
                    in_horizontal_crs.name, out_horizontal_crs.name, e
                ))
                self.transformer = None
                return

            if self.region is not None:
                self.trans_region = self.region.copy()
                self.trans_region.src_srs = out_horizontal_crs.to_proj4() #'epsg:{}'.format(out_horizontal_epsg)
                self.trans_region.warp(in_horizontal_crs.to_proj4())

            self.src_proj4 = in_horizontal_crs.to_proj4()
            self.dst_proj4 = out_horizontal_crs.to_proj4()
                
            ## vertical Transformation
            if want_vertical:
                if self.region is None:
                    vd_region = regions.Region().from_list(self.infos.minmax)
                    vd_region.src_srs = in_horizontal_crs.to_proj4()
                else:
                    vd_region = self.region.copy()
                    vd_region.src_srs = out_horizontal_crs.to_proj4()
                    
                vd_region.warp('epsg:4326')
                if not vd_region.valid_p():
                    utils.echo_warning_msg('failed to generate transformation')
                    return
                
                vd_region.zmin = None
                vd_region.zmax = None
                vd_region.buffer(pct=5)
                
                ## trans_fn is the transformation grid, used in gdalwarp
                self.trans_fn = os.path.join(
                    self.cache_dir, '_vdatum_trans_{}_{}_{}.tif'.format(
                        in_vertical_epsg, out_vertical_epsg, vd_region.format('fn')
                    )
                )

                ## vertical transformation grid is generated in WGS84
                if not os.path.exists(self.trans_fn):
                    with tqdm(
                            desc='generating vertical transformation grid {} from {} to {}'.format(
                                self.trans_fn, in_vertical_epsg, out_vertical_epsg
                            ),
                            leave=self.verbose
                    ) as pbar:
                        vd_x_inc = vd_y_inc = utils.str2inc('3s')
                        xcount, ycount, dst_gt = vd_region.geo_transform(
                            x_inc=vd_x_inc, y_inc=vd_y_inc, node='grid'
                        )

                        while (xcount <=10 or ycount <=10):
                            vd_x_inc /= 2
                            vd_y_inc /= 2
                            xcount, ycount, dst_gt = vd_region.geo_transform(
                                x_inc=vd_x_inc, y_inc=vd_y_inc, node='grid'
                            )

                        self.trans_fn, self.trans_fn_unc = vdatums.VerticalTransform(
                            'IDW', vd_region, vd_x_inc, vd_y_inc, in_vertical_epsg, out_vertical_epsg,
                            geoid_in=src_geoid, geoid_out=dst_geoid, cache_dir=self.cache_dir, verbose=False
                        ).run(outfile=self.trans_fn)                        

                if self.trans_fn is not None and os.path.exists(self.trans_fn):
                    #utils.echo_msg('using vertical tranformation grid {} from {} to {}'.format(self.trans_fn, in_vertical_epsg, out_vertical_epsg))
                    out_src_srs = '{} +geoidgrids={}'.format(in_horizontal_crs.to_proj4(), self.trans_fn)

                    if utils.str_or(in_vertical_epsg) == '6360':# or 'us-ft' in utils.str_or(src_vert, ''):
                        out_src_srs = out_src_srs + ' +vto_meter=0.3048006096012192'
                        self.trans_to_meter = True
                    
                    # if utils.str_or(out_vertical_epsg) == '6360':# or 'us-ft' in utils.str_or(dst_vert, ''):
                    #     out_dst_srs = out_dst_srs + ' +vto_meter=0.3048006096012192'
                    #     self.trans_from_meter = True
                else:
                    utils.echo_error_msg(
                        'failed to generate vertical transformation grid between {} and {} for this region!'.format(
                            in_vertical_epsg, out_vertical_epsg
                        )
                    )
                
                in_vertical_crs = pyproj.CRS.from_user_input(out_src_srs)
                self.src_proj4 = in_vertical_crs.to_proj4()
                self.dst_proj4 = out_horizontal_crs.to_proj4()
                self.aux_src_proj4 = in_horizontal_crs.to_proj4()
                self.aux_dst_proj4 = out_horizontal_crs.to_proj4()
                #utils.echo_msg('{} {}'.format(in_vertical_crs, out_horizontal_crs))
                if self.region is not None:
                    aoi = pyproj.aoi.AreaOfInterest(self.region.xmin, self.region.ymin, self.region.xmax, self.region.ymax)
                else:
                    aoi = None
                    
                self.transformer = pyproj.Transformer.from_crs(in_vertical_crs, out_horizontal_crs, always_xy=True, area_of_interest=aoi)

        ## dataset region
        if self.region is not None and self.region.valid_p():
            self.data_region = self.region.copy() if self.trans_region is None else self.trans_region.copy()
            inf_region = regions.Region().from_list(self.infos.minmax)
            self.data_region = regions.regions_reduce(self.data_region, inf_region)
            self.data_region.src_srs = self.infos.src_srs

            if not self.data_region.valid_p():
                self.data_region = self.region.copy() if self.trans_region is None else self.trans_region.copy()
        else:
            self.data_region = regions.Region().from_list(self.infos.minmax)
            self.data_region.src_srs = self.infos.src_srs
                    
    def parse(self):
        """parse the datasets from the dataset.
        
        Re-define this method when defining a dataset sub-class
        that represents recursive data structures (datalists, zip, etc).

        This will fill self.data_entries

        -------
        Yields:

        dataset object
        """

        if self.region is not None:
            try:
                inf_region = regions.Region().from_string(self.infos.wkt)
            except:
                try:
                    inf_region = regions.Region().from_list(self.infos.minmax)
                except:
                    inf_region = self.region.copy()

            if regions.regions_intersect_p(
                    inf_region,
                    self.region if self.trans_region is None else self.trans_region
            ):
                self.data_entries.append(self)
                yield(self)

        else:
            self.data_entries.append(self)
            yield(self)

    def parse_data_lists(self, gather_data=True):
        """parse the data into a datalist dictionary"""
        
        for e in self.parse():
            if e.parent is not None:
                if e.parent.metadata['name'] in self.data_lists.keys():
                    self.data_lists[e.parent.metadata['name']]['data'].append(e)
                else:
                    self.data_lists[e.parent.metadata['name']] = {'data': [e], 'parent': e.parent}
                    
            else:
                self.data_lists[e.metadata['name']] = {'data': [e], 'parent': e}
                
        return(self)
    
    def archive_xyz(self, **kwargs):
        """Archive data from the dataset to XYZ in the given dataset region.
        
        will convert all data to XYZ within the given region and will arrange
        the data as a datalist based on inputs.

        Data comes from `self.xyz_yield` set in `self.set_yield`. So will process data through
        `_stacks` before archival if region, x_inc and y_inc are set.
        """
        
        def xdl2dir(xdl):
            this_dir = []
            while True:
                if xdl.parent is None:
                    out_dir = xdl.metadata['name']
                    if xdl.remote:
                        out_dir = out_dir.split(':')[0]

                    this_dir.append(out_dir)
                    break
                
                out_dir = xdl.parent.metadata['name']
                if xdl.parent.remote:
                    out_dir = utils.fn_bansename2(out_dir)
                    
                this_dir.append(out_dir)
                xdl = xdl.parent
                
            this_dir.reverse()
            return(this_dir)

        aa_name = self.metadata['name'].split(':')[0]
        if self.region is None:
            a_name = '{}_{}'.format(aa_name, utils.this_year())
        else:
            a_name = '{}_{}_{}'.format(
                aa_name, self.region.format('fn'), utils.this_year())

        self.archive_datalist = '{}.datalist'.format(a_name)
        archive_keys = []
        with tqdm(desc='archiving datasets to {}'.format(self.archive_datalist), leave=self.verbose) as pbar:
            with open('{}.datalist'.format(a_name), 'w') as dlf:
                for this_entry in self.parse():
                    pbar.update()
                    if this_entry.parent is None:
                        this_key = this_entry.metadata['name'].split(':')[0]
                        this_dir = []
                    else:
                        this_key = this_entry.parent.metadata['name'].split(':')[0]
                        this_dir = xdl2dir(this_entry.parent)
                        
                    if self.region is None:
                        a_dir = '{}_{}'.format(this_key, utils.this_year())
                    else:
                        a_dir = '{}_{}_{}'.format(this_key, self.region.format('fn'), utils.this_year())
                        
                    this_dir.append(a_dir)
                    if not this_key in archive_keys:
                        archive_keys.append(this_key)                    
                        dlf.write(
                            '{name}.datalist -1 {weight} {uncertainty} {metadata}\n'.format(
                                name=os.path.join(*(this_dir + [this_dir[-1]])),
                                weight = utils.float_or(this_entry.parent.weight, 1) if this_entry.parent is not None else 1,
                                uncertainty = utils.float_or(this_entry.parent.uncertainty, 0) if this_entry.parent is not None else 0,
                                metadata = this_entry.parent.format_metadata() if this_entry.parent is not None else this_entry.format_metadata()
                            )
                        )
                        
                    this_dir = os.path.join(os.getcwd(), *this_dir)
                    if not os.path.exists(this_dir):
                        os.makedirs(this_dir)

                    ## todo: if metadata['name'] contains a parent directory, make that into its own datalist
                    with open(
                            os.path.join(
                                this_dir, '{}.datalist'.format(os.path.basename(this_dir))
                            ), 'a'
                    ) as sub_dlf:
                        pbar.update()
                        #utils.echo_msg(this_entry)
                        #utils.echo_msg(this_entry.metadata)
                        if this_entry.remote == True:
                            #sub_xyz_path = '{}_{}.xyz'.format(
                            #    this_entry.metadata['name'], self.region.format('fn')
                            #)
                            sub_xyz_path = '{}_{}.xyz'.format(
                                os.path.basename(utils.fn_basename2(this_entry.fn)), self.region.format('fn')
                            )
                                
                        elif len(this_entry.fn.split('.')) > 1:
                            xyz_ext = this_entry.fn.split('.')[-1]
                            sub_xyz_path = '.'.join(
                                [utils.fn_basename(
                                    utils.slugify(
                                        os.path.basename(this_entry.fn)
                                    ),
                                    this_entry.fn.split('.')[-1]),
                                 'xyz']
                            )
                        else:
                            sub_xyz_path = '.'.join([utils.fn_basename2(os.path.basename(this_entry.fn)), 'xyz'])
                            
                        this_xyz_path = os.path.join(this_dir, sub_xyz_path)
                        if os.path.exists(this_xyz_path):
                            utils.echo_warning_msg('{} already exists, skipping...'.format(this_xyz_path))
                            continue
                        
                        if not os.path.exists(os.path.dirname(this_xyz_path)):
                            os.makedirs(os.path.dirname(this_xyz_path))
                                                        
                        sub_dirname = os.path.dirname(sub_xyz_path)
                        if sub_dirname != '.' and sub_dirname != '':
                            sub_sub_dlf_path = os.path.join(this_dir, sub_dirname, '{}.datalist'.format(sub_dirname))
                            if not os.path.exists(os.path.dirname(sub_sub_dlf_path)):
                                os.makedirs(os.path.dirname(sub_sub_dlf_path))

                            sub_sub_dlf = open(sub_sub_dlf_path, 'a')
                            if not sub_dirname in archive_keys:
                                archive_keys.append(sub_dirname)
                                sub_dlf.write('{} -1 1 0\n'.format(os.path.relpath(sub_sub_dlf_path, this_dir)))

                            sub_sub_dlf.write('{} 168 1 0\n'.format(os.path.relpath(this_xyz_path, os.path.dirname(sub_sub_dlf_path))))
                            with open(this_xyz_path, 'w') as xp:
                                for this_xyz in this_entry.xyz_yield: # data will be processed independently of each other
                                    #yield(this_xyz) # don't need to yield data here.
                                    this_xyz.dump(include_w=True if self.weight is not None else False,
                                                  include_u=True if self.uncertainty is not None else False,
                                                  dst_port=xp, encode=False)
                            sub_sub_dlf.close()
                            
                        else:
                            sub_dlf.write('{} 168 1 0\n'.format(sub_xyz_path))
                            with open(this_xyz_path, 'w') as xp:
                                for this_xyz in this_entry.xyz_yield: # data will be processed independently of each other
                                    #yield(this_xyz) # don't need to yield data here.
                                    this_xyz.dump(include_w=True if self.weight is not None else False,
                                                  include_u=True if self.uncertainty is not None else False,
                                                  dst_port=xp, encode=False, precision=self.dump_precision)
                                
        ## generate datalist inf/json
        this_archive = DatasetFactory(mod=self.archive_datalist, data_format=-1, parent=None, weight=1,
                                      uncertainty=0)._acquire_module().initialize()
        this_archive.inf()

    def yield_block_array(self): #*depreciated*
        """yield the xyz data as arrays, for use in `array_yield` or `yield_array`

        Yields:
          tuple: (list of arrays for [weighted]mean, weights mask, uncertainty mask, and count, srcwin, geotransform)
        """

        out_arrays = {'z':None, 'count':None, 'weight':None, 'uncertainty': None, 'mask':None}
        xcount, ycount, dst_gt = self.region.geo_transform(
            x_inc=self.x_inc, y_inc=self.y_inc, node='grid'
        )

        for this_xyz in self.yield_xyz():
            xpos, ypos = utils._geo2pixel(
                this_xyz.x, this_xyz.y, dst_gt, 'pixel'
            )
            if xpos < xcount and ypos < ycount and xpos >= 0 and ypos >= 0:
                w = this_xyz.w if self.weight is not None else 1
                u = this_xyz.u if self.uncertainty is not None else 0
                out_arrays['z'] = np.array([[this_xyz.z]])
                out_arrays['count'] = np.array([[1]])
                out_arrays['weight'] = np.array([[this_xyz.w]])
                out_arrays['uncertainty'] = np.array([[this_xyz.u]])
                this_srcwin = (xpos, ypos, 1, 1)
                yield(out_arrays, this_srcwin, dst_gt)

    ## todo: properly mask supercede mode...
    ## todo: 'separate mode': multi-band z?
    def _stacks(self, out_name = None, ndv = -9999, fmt = 'GTiff', mask_only = False):
        """stack and mask incoming arrays (from `array_yield`) together

        -----------
        Parameters:
        out_name (str): the output stacked raster basename
        ndv (float): the desired no data value
        fmt (str): the output GDAL file format
        mask_only (bool): only generate a mask, don't stack...

        --------
        Returns:
        output-file-name{_msk} of a multi-band raster with a band for each dataset.
        output-file-name of a multi-band raster with the following bands:
          x
          y
          z
          weights
          count
          uncertainty
          src uncertainty
        """

        utils.set_cache(self.cache_dir)
        if self.stack_mode not in ['mean', 'min', 'max', 'supercede']:
            mode = 'mean'
        else:
            mode = self.stack_mode
            
        utils.echo_msg('stacking using {} mode'.format(mode))        
        ## initialize the output rasters
        if out_name is None:
            out_name = os.path.join(self.cache_dir, '{}'.format(
                utils.append_fn('_dlim_stacks', self.region, self.x_inc)
            ))

        out_file = '{}.{}'.format(out_name, gdalfun.gdal_fext(fmt))
        mask_fn = '{}_msk.{}'.format(out_name, gdalfun.gdal_fext(fmt))
        xcount, ycount, dst_gt = self.region.geo_transform(
            x_inc=self.x_inc, y_inc=self.y_inc, node='grid'
        )
        if xcount <= 0 or ycount <=0:
            utils.echo_error_msg(
                'could not create grid of {}x{} cells with {}/{} increments on region: {}'.format(
                    xcount, ycount, self.x_inc, self.y_inc, self.region
                )
            )
            sys.exit(-1)

        gdt = gdal.GDT_Float32
        c_gdt = gdal.GDT_Int32
        driver = gdal.GetDriverByName(fmt)
        if os.path.exists(out_file):
            status = driver.Delete(out_file)
            if status != 0:
                utils.remove_glob('{}*'.format(out_file))
                
        if os.path.exists(mask_fn):
            status = driver.Delete(mask_fn)
            if status != 0:
                utils.remove_glob('{}*'.format(mask_fn))

        dst_ds = driver.Create(out_file, xcount, ycount, 7, gdt,
                               options=['COMPRESS=LZW', 'PREDICTOR=2', 'TILED=YES', 'BIGTIFF=YES'] if fmt != 'MEM' else [])

        if dst_ds is None:
            utils.echo_error_msg('failed to create stack grid {} {} {} {} {}...'.format(out_file, xcount, ycount, gdt, fmt))
            sys.exit(-1)

        dst_ds.SetGeoTransform(dst_gt)
        stacked_bands = {'z': dst_ds.GetRasterBand(1), 'count': dst_ds.GetRasterBand(2),
                         'weights': dst_ds.GetRasterBand(3), 'uncertainty': dst_ds.GetRasterBand(4),
                         'src_uncertainty': dst_ds.GetRasterBand(5), 'x': dst_ds.GetRasterBand(6), 'y': dst_ds.GetRasterBand(7) }
        stacked_data = {'z': None, 'count': None, 'weights': None, 'uncertainty': None, 'src_uncertainty': None, 'x': None, 'y': None}
        
        for key in stacked_bands.keys():
            stacked_bands[key].SetNoDataValue(np.nan)
            stacked_bands[key].SetDescription(key)

        ## incoming arrays arrs['z'], arrs['weight'] arrs['uncertainty'], and arrs['count']
        ## srcwin is the srcwin of the waffle relative to the incoming arrays
        ## gt is the geotransform of the incoming arrays
        ## mask grid
        driver = gdal.GetDriverByName('MEM')
        m_ds = driver.Create(utils.make_temp_fn(out_name), xcount, ycount, 0, gdt)
        m_ds.SetGeoTransform(dst_gt)
        
        ## initialize data mask        
        ## parse each entry and process it
        ## todo: mask here instead of in each dataset module
        for this_entry in self.parse():
            ## MASK
            m_bands = {m_ds.GetRasterBand(i).GetDescription(): i for i in range(1, m_ds.RasterCount + 1)}
            if not this_entry.metadata['name'] in m_bands.keys():
                m_ds.AddBand()
                m_band = m_ds.GetRasterBand(m_ds.RasterCount)
                m_band.SetNoDataValue(0)
                m_band.SetDescription(this_entry.metadata['name'])
                band_md = m_band.GetMetadata()
                for k in this_entry.metadata.keys():
                    band_md[k] = this_entry.metadata[k]

                band_md['weight'] = this_entry.weight
                band_md['uncertainty'] = this_entry.uncertainty
                try:
                    m_band.SetMetadata(band_md)
                except:
                    try:
                        for key in band_md.keys():
                            if band_md[key] is None:
                                del band_md[key]
                                
                        m_band.SetMetadata(band_md)
                    except Exception as e:
                        utils.echo_error_msg('could not set band metadata: {}; {}'.format(band_md, e))
            else:
                m_band = m_ds.GetRasterBand(m_bands[this_entry.metadata['name']])

            ## yield entry arrays for stacks
            #for arrs, srcwin, gt in this_entry.yield_array():
            for arrs, srcwin, gt in this_entry.array_yield:
                ## update the mask
                m_array = m_band.ReadAsArray(srcwin[0], srcwin[1], srcwin[2], srcwin[3])
                m_array[arrs['count'] != 0] = 1
                m_band.WriteArray(m_array, srcwin[0], srcwin[1])
                if mask_only:
                    continue
                
                m_ds.FlushCache()
                ## Read the saved accumulated rasters at the incoming srcwin and set ndv to zero
                for key in stacked_bands.keys():
                    stacked_data[key] = stacked_bands[key].ReadAsArray(srcwin[0], srcwin[1], srcwin[2], srcwin[3])
                    if mode != 'min' and mode != 'max':
                        stacked_data[key][np.isnan(stacked_data[key])] = 0
                    #else:
                    if key == 'count':
                        stacked_data[key][np.isnan(stacked_data[key])] = 0
                    
                ## set incoming np.nans to zero and mask to non-nan count
                arrs['count'][np.isnan(arrs['count'])] = 0
                arrs['weight'][np.isnan(arrs['z'])] = 0
                arrs['uncertainty'][np.isnan(arrs['z'])] = 0
                
                if mode != 'min' and mode != 'max':
                    #arrs['weight'][np.isnan(arrs['z'])] = 0
                    #arrs['uncertainty'][np.isnan(arrs['z'])] = 0
                    arrs['x'][np.isnan(arrs['x'])] = 0
                    arrs['y'][np.isnan(arrs['y'])] = 0
                    arrs['z'][np.isnan(arrs['z'])] = 0
                    for arr_key in arrs:
                        if arrs[arr_key] is not None:
                            arrs[arr_key][np.isnan(arrs[arr_key])] = 0
                #else:


                ## add the count to the accumulated rasters
                stacked_data['count'] += arrs['count']

                ## supercede based on weights, else do weighted mean
                ## todo: do (weighted) mean on cells with same weight
                if mode == 'supercede':
                    ## higher weight supercedes lower weight (first come first served atm)
                    stacked_data['z'][arrs['weight'] > stacked_data['weights']] = arrs['z'][arrs['weight'] > stacked_data['weights']]
                    stacked_data['x'][arrs['weight'] > stacked_data['weights']] = arrs['x'][arrs['weight'] > stacked_data['weights']]
                    stacked_data['y'][arrs['weight'] > stacked_data['weights']] = arrs['y'][arrs['weight'] > stacked_data['weights']]
                    stacked_data['src_uncertainty'][arrs['weight'] > stacked_data['weights']] = arrs['uncertainty'][arrs['weight'] > stacked_data['weights']]
                    stacked_data['weights'][arrs['weight'] > stacked_data['weights']] = arrs['weight'][arrs['weight'] > stacked_data['weights']]
                    #stacked_data['weights'][stacked_data['weights'] == 0] = np.nan
                    ## uncertainty is src_uncertainty, as only one point goes into a cell
                    stacked_data['uncertainty'][:] = stacked_data['src_uncertainty'][:]

                    # ## reset all data where weights are zero to nan
                    # for key in stacked_bands.keys():
                    #     stacked_data[key][np.isnan(stacked_data['weights'])] = np.nan

                elif mode == 'min' or mode == 'max':
                    ## set nodata values in stacked_data to whatever the value is in arrs
                    mask = np.isnan(stacked_data['z'])
                    stacked_data['x'][mask] = arrs['x'][mask]
                    stacked_data['y'][mask] = arrs['y'][mask]
                    stacked_data['src_uncertainty'][mask] = arrs['uncertainty'][mask]
                    stacked_data['uncertainty'][mask] = arrs['uncertainty'][mask]
                    stacked_data['weights'][mask] = arrs['weight'][mask]
                    stacked_data['z'][mask] = arrs['z'][mask]

                    ## mask the min or max and apply it to stacked_data 
                    if mode == 'min':
                        mask = arrs['z'] <= stacked_data['z']
                    else:
                        mask = arrs['z'] >= stacked_data['z']
                        
                    stacked_data['x'][mask] = arrs['x'][mask]
                    stacked_data['y'][mask] = arrs['y'][mask]
                    #stacked_data['src_uncertainty'][mask]  += arrs['uncertainty'][mask]
                    #stacked_data['uncertainty'][mask] = arrs['uncertainty'][mask]
                    #stacked_data['weights'][mask] = arrs['weight'][mask]

                    ## accumulate uncertainty and weight
                    stacked_data['src_uncertainty'][mask] += (arrs['uncertainty'][mask] * arrs['weight'][mask])

                    ## accumulate incoming weights (weight*weight?) and set results to np.nan for calcs
                    stacked_data['weights'][mask] += arrs['weight'][mask]
                    stacked_data['weights'][mask][stacked_data['weights'][mask] == 0] = np.nan
                    
                    ## accumulate variance * weight
                    stacked_data['uncertainty'][mask] += arrs['weight'][mask] * np.power((arrs['z'][mask] - (stacked_data['z'][mask] / stacked_data['weights'][mask])), 2)
                    stacked_data['z'][mask] = arrs['z'][mask]
                    
                # elif mode == 'median':
                #     ## accumulate incoming median(z), etc.
                #     stacked_data['z'] += np.median([arrs['z']])
                #     stacked_data['x'] += np.median(arrs['x'])
                #     stacked_data['y'] += np.median(arrs['y'])
                #     stacked_data['src_uncertainty'] += np.median(arrs['uncertainty'])
                    
                elif mode == 'mean':
                    ## accumulate incoming z*weight and uu*weight
                    stacked_data['z'] += (arrs['z'] * arrs['weight'])
                    stacked_data['x'] += (arrs['x'] * arrs['weight'])
                    stacked_data['y'] += (arrs['y'] * arrs['weight'])
                    stacked_data['src_uncertainty'] += (arrs['uncertainty'] * arrs['weight'])

                    ## accumulate incoming weights (weight*weight?) and set results to np.nan for calcs
                    stacked_data['weights'] += arrs['weight']
                    stacked_data['weights'][stacked_data['weights'] == 0] = np.nan
                    
                    ## accumulate variance * weight
                    stacked_data['uncertainty'] += arrs['weight'] * np.power((arrs['z'] - (stacked_data['z'] / stacked_data['weights'])), 2)

                ## write out results to accumulated rasters
                #stacked_data['count'][stacked_data['count'] == 0] = ndv
                stacked_data['count'][stacked_data['count'] == 0] = np.nan

                for key in stacked_bands.keys():
                    #stacked_data[key][np.isnan(stacked_data[key])] = ndv
                    #stacked_data[key][stacked_data['count'] == ndv] = ndv
                    
                    stacked_data[key][np.isnan(stacked_data['count'])] = np.nan
                    #if mode != 'mean':
                    #    stacked_data[key][np.isnan(stacked_data[key])] = ndv

                    stacked_bands[key].WriteArray(stacked_data[key], srcwin[0], srcwin[1])

        ## Finalize weighted mean rasters and close datasets
        ## incoming arrays have all been processed, if weighted mean the
        ## "z" is the sum of z*weight, "weights" is the sum of weights
        ## "uncertainty" is the sum of variance*weight
        if self.verbose:
            utils.echo_msg('finalizing stacked raster bands...')

        if m_ds.RasterCount > 0:
            #m_ds.FlushCache()
            ## create a new mem ds to hold valid bands
            driver = gdal.GetDriverByName('MEM')
            mm_ds = driver.Create(utils.make_temp_fn(out_name), xcount, ycount, 0, gdt)
            mm_ds.SetGeoTransform(dst_gt)

            for band_num in range(1, m_ds.RasterCount+1):
                band_infos = gdalfun.gdal_infos(m_ds, scan=True, band=band_num)
                if not np.isnan(band_infos['zr'][0]) and not np.isnan(band_infos['zr'][1]):
                    m_band = m_ds.GetRasterBand(band_num)
                    m_band_md = m_band.GetMetadata()
                    mm_ds.AddBand()
                    mm_band = mm_ds.GetRasterBand(mm_ds.RasterCount)
                    mm_band.SetNoDataValue(0)
                    mm_band.SetDescription(m_band.GetDescription())
                    
                    mm_band.SetMetadata(m_band_md)

                    m_array = m_band.ReadAsArray()
                    mm_band.WriteArray(m_array)
                    
                    mm_ds.FlushCache()
                    
            msk_ds = gdal.GetDriverByName(fmt).CreateCopy(mask_fn, mm_ds, 0)
            mm_ds = None
        else:
            if self.verbose:
                utils.echo_msg('no bands found for {}'.format(mask_fn))

        if not mask_only:

            # ## by moving window
            # n_chunk = int(xcount)
            # n_step = n_chunk
            # for srcwin in utils.yield_srcwin(
            #         (dst_ds.RasterYSize, dst_ds.RasterXSize), n_chunk=n_chunk, step=n_step, verbose=True
            # ):
            #     y = srcwin[1]
            #     for key in stacked_bands.keys():
            #         stacked_data[key] = stacked_bands[key].ReadAsArray(*srcwin)
            #         stacked_data[key][stacked_data[key] == ndv] = np.nan

            ## by scan-line
            srcwin = (0, 0, dst_ds.RasterXSize, dst_ds.RasterYSize)
            for y in range(
                    srcwin[1], srcwin[1] + srcwin[3], 1
            ):
                for key in stacked_bands.keys():
                    stacked_data[key] = stacked_bands[key].ReadAsArray(srcwin[0], y, srcwin[2], 1)
                    stacked_data[key][stacked_data[key] == ndv] = np.nan

                if mode == 'mean' or mode == 'min' or mode == 'max':
                    stacked_data['weights'] = stacked_data['weights'] / stacked_data['count']
                    if mode == 'mean':
                        ## average the accumulated arrays for finalization
                        ## x, y, z and u are weighted sums, so divide by weights
                        stacked_data['x'] = (stacked_data['x'] / stacked_data['weights']) / stacked_data['count']
                        stacked_data['y'] = (stacked_data['y'] / stacked_data['weights']) / stacked_data['count']
                        stacked_data['z'] = (stacked_data['z'] / stacked_data['weights']) / stacked_data['count']
                    
                    ## apply the source uncertainty with the sub-cell variance uncertainty
                    ## caclulate the standard error (sqrt( uncertainty / count))
                    stacked_data['src_uncertainty'] = np.sqrt((stacked_data['src_uncertainty'] / stacked_data['weights']) / stacked_data['count'])
                    stacked_data['uncertainty'] = np.sqrt((stacked_data['uncertainty'] / stacked_data['weights']) / stacked_data['count'])
                    stacked_data['uncertainty'] = np.sqrt(np.power(stacked_data['src_uncertainty'], 2) + np.power(stacked_data['uncertainty'], 2))
                    #stacked_data['uncertainty'] = np.sqrt(stacked_data['uncertainty'] / stacked_data['count'])

                ## write out final rasters
                for key in stacked_bands.keys():
                    stacked_data[key][np.isnan(stacked_data[key])] = ndv
                    stacked_bands[key].WriteArray(stacked_data[key], srcwin[0], y)
                
            ## set the final output nodatavalue
            for key in stacked_bands.keys():
                stacked_bands[key].DeleteNoDataValue()
                
            for key in stacked_bands.keys():
                stacked_bands[key].SetNoDataValue(ndv)
                
        ## create a vector of the masks (spatial-metadata)
        if self.want_sm:
            gdalfun.ogr_polygonize_multibands(msk_ds)

        m_ds = msk_ds = dst_ds = None

        ## apply any filters to the stack
        for f in self.stack_fltrs:
            grits_filter = grits.GritsFactory(mod=f, src_dem=out_file, uncertainty_mask=4)._acquire_module()
            if grits_filter is not None:
                #if 'stacks' in grits_filter.kwargs.keys():
                #if grits_filter.kwargs['stacks']:
                grits_filter()
                os.replace(grits_filter.dst_dem, out_file)
            
        # if isinstance(self.fltrs, list):
        #     for f in self.fltrs:
        #         grits_filter = grits.GritsFactory(mod=f, src_dem=out_file, uncertainty_mask=4)._acquire_module()
        #         if grits_filter is not None:
        #             if 'stacks' in grits_filter.kwargs.keys():
        #                 if grits_filter.kwargs['stacks']:
        #                     grits_filter()
        #                     os.replace(grits_filter.dst_dem, out_file)
        
        return(out_file)        
    
    def stacks_yield_xyz(self, out_name = None, ndv = -9999, fmt = 'GTiff'):#, mode = 'mean'):
        """yield the result of `_stacks` as an xyz object"""

        stacked_fn = self._stacks(out_name=out_name, ndv=ndv, fmt=fmt)#, mode=mode)
        sds = gdal.Open(stacked_fn)
        sds_gt = sds.GetGeoTransform()
        sds_z_band = sds.GetRasterBand(1) # the z band from stacks
        sds_w_band = sds.GetRasterBand(3) # the weight band from stacks
        sds_u_band = sds.GetRasterBand(4) # the uncertainty band from stacks
        sds_x_band = sds.GetRasterBand(6) # the x band from stacks
        sds_y_band = sds.GetRasterBand(7) # the y band from stacks
        srcwin = (0, 0, sds.RasterXSize, sds.RasterYSize)
        for y in range(srcwin[1], srcwin[1] + srcwin[3], 1):
            sz = sds_z_band.ReadAsArray(srcwin[0], y, srcwin[2], 1)
            ## skip row if all values are ndv
            if np.all(sz == ndv):
                continue
            
            sw = sds_w_band.ReadAsArray(srcwin[0], y, srcwin[2], 1)
            su = sds_u_band.ReadAsArray(srcwin[0], y, srcwin[2], 1)

            sx = sds_x_band.ReadAsArray(srcwin[0], y, srcwin[2], 1)
            sy = sds_y_band.ReadAsArray(srcwin[0], y, srcwin[2], 1)
            for x in range(0, sds.RasterXSize):
                z = sz[0, x]
                if z != ndv:
                    if self.stack_node: # yield avg x/y data instead of center
                        out_xyz = xyzfun.XYZPoint(
                            x = sx[0, x], y = sy[0, x], z = z, w = sw[0, x], u = su[0, x]
                        )
                    else: # yield center of pixel as x/y
                        geo_x, geo_y = utils._pixel2geo(x, y, sds_gt)
                        out_xyz = xyzfun.XYZPoint(
                            x=geo_x, y=geo_y, z=z, w=sw[0, x], u=su[0, x]
                        )
                        
                    yield(out_xyz)
        sds = None
    
    def vectorize_xyz(self):
        """Make a point vector OGR DataSet Object from src_xyz

        for use in gdal gridding functions
        """

        dst_ogr = '{}'.format(self.metadata['name'])
        ogr_ds = gdal.GetDriverByName('Memory').Create(
            '', 0, 0, 0, gdal.GDT_Unknown
        )
        layer = ogr_ds.CreateLayer(
            dst_ogr,
            geom_type=ogr.wkbPoint25D
        )
        for x in ['long', 'lat', 'elev', 'weight']:
            fd = ogr.FieldDefn(x, ogr.OFTReal)
            fd.SetWidth(12)
            fd.SetPrecision(8)
            layer.CreateField(fd)
            
        f = ogr.Feature(feature_def=layer.GetLayerDefn())        
        for this_xyz in self.xyz_yield:
            f.SetField(0, this_xyz.x)
            f.SetField(1, this_xyz.y)
            f.SetField(2, this_xyz.z)
            f.SetField(3, this_xyz.w)
                
            wkt = this_xyz.export_as_wkt(include_z=True)
            g = ogr.CreateGeometryFromWkt(wkt)
            f.SetGeometryDirectly(g)
            layer.CreateFeature(f)
            
        return(ogr_ds)

    def yield_points(self):
        """points are an array of `points['x']`, `points['y']`, `points['z']`, <`points['w']`, `points['u']`>`

        points will be transformed here, based on `self.transformer`, which is set in `set_transform`
        after the points are transformed, the data will pass through the region, if it exists and finally
        yield the transformed and reduced points.
        """

        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            for points in self.yield_ds(): 
                if self.transformer is not None:
                    points['x'], points['y'], points['z'] = self.transformer.transform(points['x'], points['y'], points['z'])
                    points = points[~np.isinf(points['z'])]
                    
                if self.region is not None and self.region.valid_p():
                    xyz_region = self.region.copy() #if self.trans_region is None else self.trans_region.copy()
                    if self.invert_region:
                        points = points[((points['x'] > xyz_region.xmax) | (points['x'] < xyz_region.xmin)) | \
                                        ((points['y'] > xyz_region.ymax) | (points['y'] < xyz_region.ymin))]
                        if xyz_region.zmin is not None:
                            points =  points[(points['z'] < xyz_region.zmin)]

                        if xyz_region.zmax is not None:
                            points =  points[(points['z'] > xyz_region.zmax)]
                    else:
                        points = points[((points['x'] < xyz_region.xmax) & (points['x'] > xyz_region.xmin)) & \
                                        ((points['y'] < xyz_region.ymax) & (points['y'] > xyz_region.ymin))]
                        if xyz_region.zmin is not None:
                            points =  points[(points['z'] > xyz_region.zmin)]

                        if xyz_region.zmax is not None:
                            points =  points[(points['z'] < xyz_region.zmax)]

                if len(points) > 0:

                    ## apply any dlim filters to the points
                    #if isinstance(self.fltrs, list):
                    #for f in self.fltrs:
                    if self.pnt_fltrs is not None:
                        for f in self.pnt_fltrs:
                            #utils.echo_msg(f)
                            point_filter = PointFilterFactory(mod=f, points=points)._acquire_module()
                            if point_filter is not None:
                                points = point_filter()

                    # for f in self.pnt_fltrs:
                    #     opts = f.split(':')
                    #     if len(opts) > 1:
                    #         fltr_args = utils.args2dict(list(opts[1:]), {})
                    #         fltr_args_1 = {}
                    #         for k in fltr_args.keys():
                    #             if k in ['y_res', 'z_res', 'percentile', 'z_max', 'z_min']:
                    #                 arg_val = utils.float_or(fltr_args[k])
                    #                 if arg_val is not None:
                    #                     fltr_args_1[k] = arg_val
                    #     else:
                    #         fltr_args_1 = {}            

                    #     ## bin-filter the incoming points
                    #     b_points = self.bin_z_points(points, **fltr_args_1)
                    #     if b_points is not None:
                    #         points = b_points

                    if len(points) > 0:
                        yield(points)
        
        self.transformer = None
        
    def yield_xyz(self):
        """Yield the data as xyz points

        incoming data are numpy rec-arrays of x,y,z<w,u> points.

        if 'w' (weight) is not set by the dataset, they will be given the
        default value from `self.weight`, otherwise the incoming weight values will be
        multiplied by `self.weight`

        if 'u' (uncertainty) is not set by the dataset, they will be given the
        default value from `self.uncertainty`, otherwise the incoming uncertainty values
        will be caluculated by `np.sqrt(u**2 + self.uncertainty**2)`
        """
        
        count = 0
        for points in self.yield_points():
            try:
                points_w = points['w']
            except:
                points_w = np.ones(points['z'].shape)

            points_w *= self.weight if self.weight is not None else 1
            points_w[np.isnan(points_w)] = 1
                
            try:
                points_u = points['u']
            except:
                points_u = np.zeros(points['z'].shape)

            points_u = np.sqrt(points_u**2 + (self.uncertainty if self.uncertainty is not None else 0)**2)
            points_u[np.isnan(points_u)] = 0
            dataset = np.vstack((points['x'], points['y'], points['z'], points_w, points_u)).transpose()
            count += len(dataset)
            points = None
            for point in dataset:
                this_xyz = xyzfun.XYZPoint(x=point[0], y=point[1], z=point[2], w=point[3], u=point[4])
                yield(this_xyz)

        if self.verbose:
            utils.echo_msg_bold(
                'parsed {} data records from {}{}'.format(
                    count, self.fn, ' @{}'.format(self.weight) if self.weight is not None else ''
                )
            )
                     
    def yield_array(self):
        """Yield the data as an array which coincides with the desired region, x_inc and y_inc

        incoming data are numpy rec-arrays of x,y,z<w,u> points

        if 'w' (weight) is not set by the dataset, they will be given the
        default value from `self.weight`, otherwise the incoming weight values will be
        multiplied by `self.weight`

        if 'u' (uncertainty) is not set by the dataset, they will be given the
        default value from `self.uncertainty`, otherwise the incoming uncertainty values
        will be caluculated by `np.sqrt(u**2 + self.uncertainty**2)`
        """
        
        out_arrays = {'z':None, 'count':None, 'weight':None, 'uncertainty': None, 'mask':None, 'x': None, 'y': None }
        count = 0
        for points in self.yield_points():
            xcount, ycount, dst_gt = self.region.geo_transform(
                x_inc=self.x_inc, y_inc=self.y_inc, node='grid'
            )
                
            ## convert the points to pixels based on the geotransform
            ## and calculate the local srcwin of the points
            pixel_x = np.floor((points['x'] - dst_gt[0]) / dst_gt[1]).astype(int)
            pixel_y = np.floor((points['y'] - dst_gt[3]) / dst_gt[5]).astype(int)
            # points_x = dst_gt[0] + (pixel_x+0.5) * dst_gt[1] + (pixel_y+0.5) * dst_gt[2]
            # points_y = dst_gt[3] + (pixel_x+0.5) * dst_gt[4] + (pixel_y+0.5) * dst_gt[5]

            points_x = np.array(points['x'])
            points_y = np.array(points['y'])
            
            pixel_z = np.array(points['z'])
            try:
                pixel_w = np.array(points['w'])
            except:
                pixel_w = np.ones(pixel_z.shape)

            try:
                pixel_u = np.array(points['u'])
            except:
                pixel_u = np.zeros(pixel_z.shape)
                
            ## remove pixels that will break the srcwin
            out_idx = np.nonzero((pixel_x >= xcount) | (pixel_x < 0) | (pixel_y >= ycount) | (pixel_y < 0))            
            pixel_x = np.delete(pixel_x, out_idx)
            pixel_y = np.delete(pixel_y, out_idx)
            pixel_z = np.delete(pixel_z, out_idx)
            pixel_w = np.delete(pixel_w, out_idx)
            pixel_u = np.delete(pixel_u, out_idx)
            points_x = np.delete(points_x, out_idx)
            points_y = np.delete(points_y, out_idx)
            
            points = None
            pixel_w[np.isnan(pixel_w)] = 1
            pixel_u[np.isnan(pixel_u)] = 0

            ## set the srcwin of the incoming points
            this_srcwin = (int(min(pixel_x)), int(min(pixel_y)),
                           int(max(pixel_x) - min(pixel_x))+1,
                           int(max(pixel_y) - min(pixel_y))+1)
            count += len(pixel_x)

            ## adjust the pixels to the srcwin and stack together
            pixel_x = pixel_x - this_srcwin[0]
            pixel_y = pixel_y - this_srcwin[1]
            pixel_xy = np.vstack((pixel_y, pixel_x)).T
            
            ## find the non-unique x/y points and mean/min/max their z values together
            ## while calculating the std for uncertainty
            unq, unq_idx, unq_inv, unq_cnt = np.unique(
                pixel_xy, axis=0, return_inverse=True, return_index=True, return_counts=True
            )
            cnt_msk = unq_cnt > 1
            cnt_idx, = np.nonzero(cnt_msk)
            idx_msk = np.in1d(unq_inv, cnt_idx)
            idx_idx, = np.nonzero(idx_msk)
            srt_idx = np.argsort(unq_inv[idx_msk])
            dup_idx = np.split(idx_idx[srt_idx], np.cumsum(unq_cnt[cnt_msk])[:-1])
            zz = pixel_z[unq_idx]
            ww = pixel_w[unq_idx]
            uu = pixel_u[unq_idx]
            xx = points_x[unq_idx]
            yy = points_y[unq_idx]
            #u = np.zeros(zz.shape)
            if np.any([len(dup) for dup in dup_idx]):
                if self.stack_mode == 'min':
                    dup_stack = [np.min(pixel_z[dup]) for dup in dup_idx]
                    dup_stds = np.zeros(dup_stack.shape)
                elif self.stack_mode == 'max':
                    dup_stack = [np.max(pixel_z[dup]) for dup in dup_idx]
                    dup_stds = np.zeros(dup_stack.shape)
                else:
                    dup_stack = [np.mean(pixel_z[dup]) for dup in dup_idx]
                    dup_stack_x = [np.mean(points_x[dup]) for dup in dup_idx]
                    dup_stack_y = [np.mean(points_y[dup]) for dup in dup_idx]
                    dup_stds = [np.std(pixel_z[dup]) for dup in dup_idx]
                    
                zz[cnt_msk] = dup_stack
                uu[cnt_msk] = dup_stds
                
            ## make the output arrays to yield
            out_x = np.zeros((this_srcwin[3], this_srcwin[2]))
            out_x[unq[:,0], unq[:,1]] = xx
            out_x[out_x == 0] = np.nan
            out_arrays['x'] = out_x

            out_y = np.zeros((this_srcwin[3], this_srcwin[2]))
            out_y[unq[:,0], unq[:,1]] = yy
            out_y[out_y == 0] = np.nan
            out_arrays['y'] = out_y
            
            out_z = np.zeros((this_srcwin[3], this_srcwin[2]))
            out_z[unq[:,0], unq[:,1]] = zz
            out_z[out_z == 0] = np.nan
            out_arrays['z'] = out_z
            
            out_arrays['count'] = np.zeros((this_srcwin[3], this_srcwin[2]))
            out_arrays['count'][unq[:,0], unq[:,1]] = unq_cnt
            
            out_arrays['weight'] = np.ones((this_srcwin[3], this_srcwin[2]))
            out_arrays['weight'][:] = self.weight if self.weight is not None else 1
            out_arrays['weight'][unq[:,0], unq[:,1]] *= (ww * unq_cnt)
            #out_arrays['weight'][unq[:,0], unq[:,1]] *= unq_cnt
            
            out_arrays['uncertainty'] = np.zeros((this_srcwin[3], this_srcwin[2]))
            #out_arrays['uncertainty'][:] = self.uncertainty if self.uncertainty is not None else 0
            out_arrays['uncertainty'][unq[:,0], unq[:,1]] = np.sqrt(uu**2 + (self.uncertainty if self.uncertainty is not None else 0)**2)
            
            yield(out_arrays, this_srcwin, dst_gt)

        if self.verbose:
            utils.echo_msg_bold(
                'parsed {} data records from {}{}'.format(
                    count, self.fn, ' @{}'.format(self.weight) if self.weight is not None else ''
                )
            )    
            
class XYZFile(ElevationDataset):
    """representing an ASCII xyz dataset stream.

    Parse data from an xyz file/stdin

    generate_inf - generate an inf file for the xyz data
    yield_xyz - yield the xyz data as xyz
    yield_array - yield the xyz data as an array
                  must set the x_inc/y_inc in the
                  super class
    
    -----------
    Parameters:
    
    delim: the delimiter of the xyz data (str)
    xpos: the position (int) of the x value
    ypos: the position (int) of the y value
    zpos: the position (int) of the z value
    wpos: the position (int) of the w value (weight)
    upos: the position (int) of the u value (uncertainty)
    skip: number of lines to skip
    x_scale: scale the x value
    y_scale: scale the y value
    z_scale: scale the z value
    x_offset: offset the x value
    y_offset: offset the y value    
    """
            
    def __init__(self, delim = None, xpos = 0, ypos = 1, zpos = 2,
                 wpos = None, upos = None, skip = 0, x_scale = 1, y_scale = 1,
                 z_scale = 1, x_offset = 0, y_offset = 0, use_numpy = True,
                 iter_rows = 1000000, **kwargs):
        
        super().__init__(**kwargs)
        self.delim = delim # the file delimiter
        self.xpos = utils.int_or(xpos, 0) # the position of the x/lat value
        self.ypos = utils.int_or(ypos, 1) # the position of the y/lon value
        self.zpos = utils.int_or(zpos, 2) # the position of the z value
        self.wpos = utils.int_or(wpos) # the position of the weight value
        self.upos = utils.int_or(upos) # the position of the uncertainty value
        self.skip = utils.int_or(skip, 0) # number of header lines to skip
        self.x_scale = utils.float_or(x_scale, 1) # multiply x by x_scale
        self.y_scale = utils.float_or(y_scale, 1) # multiply y by y_scale
        self.z_scale = utils.float_or(z_scale, 1) # multiply z by z_scale
        self.x_offset = x_offset # offset x by x_offset
        self.y_offset = utils.int_or(y_offset, 0) # offset y by y_offset
        self.rem = False # x is 360 instead of 180
        self.use_numpy = use_numpy # use numpy.loadtxt to load the xyz points
        self.iter_rows = iter_rows # max rows to process at a time

    def yield_ds(self):
        if self.delim is not None:
            xyzfun._known_delims.insert(0, self.delim)
                
        if self.x_offset == 'REM':
            self.x_offset = 0
            self.rem = True
            
        self.scoff = True if self.x_scale != 1 or self.y_scale != 1 or self.z_scale != 1 \
           or self.x_offset != 0 or self.y_offset != 0 else False

        self.field_names  = [x for x in ['x' if self.xpos is not None else None, 'y' if self.ypos is not None else None,
                                         'z' if self.zpos is not None else None, 'w' if self.wpos is not None else None,
                                         'u' if self.upos is not None else None]
                             if x is not None]
        self.field_formats = [float for x in [self.xpos, self.ypos, self.zpos, self.wpos, self.upos] if x is not None]

        #if self.use_numpy:
        try:
            if self.delim is None:
                self.guess_delim()

            skip_ = self.skip            
            with open(self.fn, 'r') as src_data:
                while True:
                    points = np.loadtxt(
                        src_data, delimiter=self.delim, comments='#', ndmin = 1, skiprows=skip_,
                        usecols=[x for x in [self.xpos, self.ypos, self.zpos, self.wpos, self.upos] if x is not None],
                        dtype={'names': self.field_names, 'formats': self.field_formats}, max_rows=self.iter_rows
                    )
                    skip_ = 0

                    if self.scoff:
                        points['x'] = (points['x'] + self.x_offset) * self.x_scale
                        points['y'] = (points['y'] + self.y_offset) * self.y_scale
                        points['z'] *= self.z_scale

                    if self.rem:
                        points['x'] = np.fmod(points['x'] + 180, 360) - 180 
                    
                    yield(points)
                    
                    if self.iter_rows is None or len(points) < self.iter_rows:
                        break

        ## old processing function used as a fallback for when numpy.loadtxt fails
        except Exception as e:
            utils.echo_warning_msg('could not load xyz data from {}, {}, falling back'.format(self.fn, e))
            if self.fn is not None:
                if os.path.exists(str(self.fn)):
                    self.src_data = open(self.fn, "r")
                else:
                    self.src_data = self.fn
            else:
                self.src_data = sys.stdin

            points_x = []
            points_y = []
            points_z = []
            points_w = []
            points_u = []
            count = 0
            skip = self.skip
            for xyz_line in self.src_data:
                if count >= skip:
                    this_xyz = self.line_delim(xyz_line)
                    if this_xyz is None:
                        continue

                    x = utils.float_or(this_xyz[self.xpos])
                    y = utils.float_or(this_xyz[self.ypos])
                    z = utils.float_or(this_xyz[self.zpos])
                    w = utils.float_or(this_xyz[self.wpos]) if self.wpos is not None else 1
                    u = utils.float_or(this_xyz[self.upos]) if self.upos is not None else 0

                    if x is None or y is None or z is None:
                        continue
                    
                    if self.scoff:
                        x = (x+self.x_offset) * self.x_scale
                        y = (y+self.y_offset) * self.y_scale
                        z *= self.z_scale

                    if self.rem:
                        x = math.fmod(x+180,360)-180 

                    if self.data_region is not None and self.data_region.valid_p():
                        try:
                            this_xyz = xyzfun.XYZPoint(
                                x=this_xyz[self.xpos], y=this_xyz[self.ypos], z=this_xyz[self.zpos]
                            )
                        except Exception as e:
                            utils.echo_error_msg('{} ; {}'.format(e, this_xyz))
                            this_xyz = xyzfun.XYZPoint()
                            
                    points_x.append(x)
                    points_y.append(y)
                    points_z.append(z)
                    points_w.append(w)
                    points_u.append(u)
                    count += 1
                else:
                    skip -= 1

            self.src_data.close()
            dataset = np.column_stack((points_x, points_y, points_z, points_w, points_u))
            points = np.rec.fromrecords(dataset, names='x, y, z, w, u')
            
            yield(points)
        
    def guess_delim(self):
        """guess the xyz delimiter"""

        if self.fn is not None:
            if os.path.exists(str(self.fn)):
                self.src_data = open(self.fn, "r")
            else:
                self.src_data = self.fn
        else:
            self.src_data = sys.stdin

        for xyz_line in self.src_data:
            for delim in xyzfun._known_delims:
                try:
                    this_xyz = xyz_line.split(delim)
                    if len(this_xyz) > 1:
                        self.delim = delim
                        break
                except:
                    pass
            break
            
        self.src_data.close()        


    def line_delim(self, xyz_line):
        """guess a line delimiter and return the split line."""

        for delim in xyzfun._known_delims:
            try:
                this_xyz = xyz_line.split(delim)
                if len(this_xyz) > 1:
                    return(this_xyz)
            except:
                pass

class YXZFile(XYZFile):
    """yxz file shortcut (167 mbdatalist)
    """
    
    def __init__(self, **kwargs):
        super().__init__(xpos = 1, ypos = 0, zpos = 2, **kwargs)
            
class LASFile(ElevationDataset):
    """representing an LAS/LAZ dataset.

    Process LAS/LAZ lidar files using pylas.
    
    get_epsg - attempt to parse the EPSG from the LAS file header
    generate_inf - generate an inf file for the LAS data
    yield_xyz - yield the LAS data as xyz
    yield_array - yield the LAS data as an array
                  must set the x_inc/y_inc in the
                  super class
    
    -----------
    Parameters:

    classes (str): a list of classes to parse, being a string with `/` seperator 
    """

    def __init__(self, classes='2/29/40', **kwargs):
        super().__init__(**kwargs)
        self.classes = [int(x) for x in classes.split('/')] # list of lidar classes to retain
        if self.src_srs is None:
            self.src_srs = self.get_epsg()
            if self.src_srs is None:
                self.src_srs = self.infos.src_srs
                                
    def valid_p(self, fmts = ['scratch']):
        """check if self appears to be a valid dataset entry"""

        if self.fn is None: # and not self.fn.startswith('http'):
            return(False)
        else:
            if os.path.exists(self.fn) :
                if os.stat(self.fn).st_size == 0:
                    return(False)
            else:
                return(False)

            try:
                lp.open(self.fn)
            except:
                utils.echo_warning_msg('{} could not be opened by the lasreader'.format(self.fn))
                return(False)
                        
        return(True)
        
    def get_epsg(self):
        with lp.open(self.fn) as lasf:
            lasf_vlrs = lasf.header.vlrs
            for vlr in lasf_vlrs:
                if vlr.record_id == 2112:
                    src_srs = vlr.string
                    return(src_srs)
            return(None)

    def generate_inf(self, callback=lambda: False):
        """generate an inf file for a lidar dataset."""
        
        with lp.open(self.fn) as lasf:
            self.infos.numpts = lasf.header.point_count
            this_region = regions.Region(xmin=lasf.header.x_min, xmax=lasf.header.x_max,
                                         ymin=lasf.header.y_min, ymax=lasf.header.y_max,
                                         zmin=lasf.header.z_min, zmax=lasf.header.z_max)
            self.infos.minmax = this_region.export_as_list(include_z=True)
            self.infos.wkt = this_region.export_as_wkt()

        self.infos.src_srs = self.src_srs if self.src_srs is not None else self.get_epsg()            
        return(self.infos)

    def yield_ds(self):
        with lp.open(self.fn) as lasf:
            try:
                for points in lasf.chunk_iterator(2_000_000):
                    points = points[(np.isin(points.classification, self.classes))]
                    dataset = np.column_stack((points.x, points.y, points.z))
                    points = np.rec.fromrecords(dataset, names='x, y, z')                    
                    yield(points)
                    
            except Exception as e:
                utils.echo_warning_msg('could not read points from lasfile {}, {}'.format(self.fn, e))
                
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
    
    def __init__(self, weight_mask = None, uncertainty_mask = None,  x_band = None, y_band = None, uncertainty_mask_to_meter = 1,
                 open_options = None, sample = None, check_path = True, super_grid = False, band_no = 1, remove_flat = False,
                 node = None, **kwargs):
        super().__init__(**kwargs)
        self.weight_mask = weight_mask # associated raster file/band holding weight data
        self.uncertainty_mask = uncertainty_mask # associated raster file/band holding uncertainty data
        self.uncertainty_mask_to_meter = uncertainty_mask_to_meter # uncertainty data factor to meters
        self.open_options = open_options # GDAL open-options
        self.sample = sample # GDAL resampling method
        self.check_path = check_path # check for the path to the input data file
        self.super_grid = super_grid # input has super_grids (force their use)
        self.band_no = band_no # band number holding elevation data
        self.tmp_elev_band = None # temporary elevation band
        self.tmp_unc_band = None # temporary uncertainty band
        self.tmp_weight_band = None # temporary weight band
        self.src_ds = None # the GDAL dataset object
        self.remove_flat = remove_flat # remove flattened data from input
        self.x_band = x_band # band holding x values
        self.y_band = y_band # band holding y values
        self.node = node # input is 'pixel' or 'grid' registered (force)

        if self.fn.startswith('http') or self.fn.startswith('/vsicurl/') or self.fn.startswith('BAG'):
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

        try to split the horizontal and vertical and them combine them...
        """
        
        if self.src_srs is None:
            return(gdalfun.gdal_get_srs(src_ds))
        else:
            return(self.src_srs)

    def generate_inf(self, callback=lambda: False):
        with gdalfun.gdal_datasource(self.fn) as src_ds:
            if src_ds is not None:

                if self.src_srs is None:
                    self.infos.src_srs = self.init_srs(src_ds)
                else:
                    self.infos.src_srs = self.src_srs

                ds_infos = gdalfun.gdal_infos(src_ds)
                this_region = regions.Region(src_srs=self.src_srs).from_geo_transform(
                    geo_transform=ds_infos['geoT'],
                    x_count=ds_infos['nx'],
                    y_count=ds_infos['ny']
                )

                #if scan:
                #    zr = src_ds.GetRasterBand(utils.int_or(self.band_no, 1)).ComputeRasterMinMax()
                
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
                if self.transformer is not None:
                    if self.trans_region is not None and self.trans_region.valid_p(
                            check_xy = True
                    ):
                        srcwin = self.trans_region.srcwin(
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
        
    def yield_ds(self):
        """initialize the raster dataset

        if x/y incs are set, will warp raster to that resolution.
        """

        if self.check_path and not os.path.exists(self.fn):
            return(None)

        ## apply any open_options that are specified
        try:
            self.open_options = self.open_options.split('/')
        except AttributeError:
            self.open_options = self.open_options
        except:
            self.open_options = None
            
        ## set up any transformations and other options
        inf_region = regions.Region().from_string(self.infos.wkt)
        self.sample_alg = self.sample if self.sample is not None else self.sample_alg
        self.dem_infos = gdalfun.gdal_infos(self.fn)
        if self.node is None:
            self.node = gdalfun.gdal_get_node(self.fn, 'pixel')
            
        self.resample_and_warp = True        
        if (self.x_inc is None and self.y_inc is None) or self.region is None:
            self.resample_and_warp = False

        if self.node == 'grid':
            self.resample_and_warp = False
            
        ndv = utils.float_or(gdalfun.gdal_get_ndv(self.fn), -9999)
        if self.region is not None:
            self.warp_region = self.region.copy()
        else:
            self.warp_region = inf_region.copy() #regions.Region().from_list(self.infos.minmax)
            if self.transformer is not None:
                self.warp_region.src_srs = self.src_srs
                self.warp_region.warp(self.dst_srs)

        tmp_elev_fn = utils.make_temp_fn('{}'.format(self.fn), temp_dir=self.cache_dir)
        tmp_unc_fn = utils.make_temp_fn('{}'.format(self.fn), temp_dir=self.cache_dir)
        tmp_weight_fn = utils.make_temp_fn('{}'.format(self.fn), temp_dir=self.cache_dir)
        #tmp_elev_fn = tmp_unc_fn = tmp_weight_fn = None

        ## resample/warp src gdal file to specified x/y inc/transformer respectively
        if self.resample_and_warp:
            if self.transformer is not None:
                self.transformer = None

            #utils.echo_msg(self.sample_alg)
            if self.sample_alg == 'auto':
                if self.stack_mode == 'min':
                    self.sample_alg = 'min'
                elif self.stack_mode == 'max':
                    self.sample_alg = 'max'
                elif self.dem_infos['geoT'][1] >= self.x_inc and (self.dem_infos['geoT'][5]*-1) >= self.y_inc:
                    self.sample_alg = 'bilinear'
                else:
                    self.sample_alg = 'average'

            tmp_ds = self.fn
            if self.open_options is not None:
                self.src_ds = gdal.OpenEx(self.fn, open_options=self.open_options)
                if self.src_ds is None:
                    self.src_ds = gdal.Open(self.fn)
            else:
                self.src_ds = gdal.Open(self.fn)

            if self.src_ds is None:
                return(None)

            ## Sample/Warp
            ## resmaple and/or warp dataset based on target srs and x_inc/y_inc
            ## doing this in MEM has a bug, fix if able
            ## extract necessary bands before warping, gadlwarp does not work with multiband rasters!
            tmp_warp = utils.make_temp_fn('{}'.format(tmp_ds), temp_dir=self.cache_dir)
            in_bands = self.src_ds.RasterCount
            src_ds_config = gdalfun.gdal_infos(self.src_ds)
            src_gt = src_ds_config['geoT']
            if in_bands > 1:
                ## the srcwin for to extract data
                if self.trans_region is not None:
                    srcwin_region = self.trans_region.copy()
                elif self.region is not None:
                    srcwin_region = self.region.copy()
                else:
                    srcwin_region = None

                if srcwin_region is not None:
                    srcwin = srcwin_region.srcwin(src_gt, src_ds_config['nx'], src_ds_config['ny'], node='pixel')
                else:
                    srcwin = None

                self.tmp_elev_band = tmp_elev_fn
                if self.verbose:
                    utils.echo_msg('extracting elevation data from {} to {}'.format(self.fn, self.tmp_elev_band))
                    
                self.tmp_elev_band, status = gdalfun.gdal_extract_band(
                    self.src_ds, self.tmp_elev_band, band=self.band_no, exclude=[], srcwin=srcwin, inverse=False
                )
                tmp_ds = self.tmp_elev_band
                if tmp_ds is None:
                    return(None)
                
                if utils.int_or(self.uncertainty_mask) is not None:
                    self.tmp_unc_band = tmp_unc_fn
                    if self.verbose:
                        utils.echo_msg('extracting uncertainty mask from {} to {}'.format(self.fn, self.tmp_unc_band))
                        
                    self.tmp_unc_band, status = gdalfun.gdal_extract_band(
                        self.src_ds, self.tmp_unc_band, band=self.uncertainty_mask, exclude=[], srcwin=srcwin, inverse=False
                    )
                    self.uncertainty_mask = self.tmp_unc_band

                if utils.int_or(self.weight_mask) is not None:                    
                    self.tmp_weight_band = tmp_weight_fn
                    if self.verbose:
                        utils.echo_msg('extracting weight mask from {} to {}'.format(self.fn, self.tmp_weight_band))
                        
                    self.tmp_weight_band, status = gdalfun.gdal_extract_band(
                        self.src_ds, self.tmp_weight_band, band=self.weight_mask, exclude=[], srcwin=srcwin, inverse=False
                    )
                    self.weight_mask = self.tmp_weight_band

            if self.remove_flat:
                #tmp_noflat = utils.make_temp_fn('tmp_flat.tif', temp_dir=self.cache_dir)
                #tmp_ds = gdalfun.gdal_remove_flats(tmp_ds, dst_dem=tmp_noflat, verbose=self.verbose)[0]
                #gdalfun.gdal_flat_to_nan(tmp_ds, verbose=self.verbose)[0]

                # src_dem = os.path.join(self.fetch_module._outdir, result[1])
                grits_filter = grits.GritsFactory(mod='flats', src_dem=tmp_ds, cache_dir=self.cache_dir)._acquire_module()
                grits_filter()
                tmp_ds = grits_filter.dst_dem
                #tmp_ds = gdal.Open(ned_fn)

                
            warp_ = gdalfun.sample_warp(tmp_ds, tmp_warp, self.x_inc, self.y_inc, src_srs=self.src_proj4, dst_srs=self.dst_proj4,
                                        src_region=self.warp_region,#if (self.y_inc is not None and self.x_inc is not None) else None,
                                        sample_alg=self.sample_alg, ndv=ndv, verbose=self.verbose,
                                        co=["COMPRESS=DEFLATE", "TILED=YES"])[0]
            tmp_ds = warp_ = None
            
            ## the following seems to be redundant...
            warp_ds = gdal.Open(tmp_warp)
            if warp_ds is not None:
                ## clip wapred ds to warped srcwin
                warp_ds_config = gdalfun.gdal_infos(warp_ds)
                gt = warp_ds_config['geoT']
                srcwin = self.warp_region.srcwin(gt, warp_ds.RasterXSize, warp_ds.RasterYSize, node='grid')
                dst_gt = (gt[0] + (srcwin[0] * gt[1]), gt[1], 0., gt[3] + (srcwin[1] * gt[5]), 0., gt[5])
                out_ds_config = gdalfun.gdal_set_infos(srcwin[2], srcwin[3], srcwin[2] * srcwin[3], dst_gt, warp_ds_config['proj'],
                                                  warp_ds_config['dt'], warp_ds_config['ndv'], warp_ds_config['fmt'], None, None)

                in_bands = warp_ds.RasterCount
                self.src_ds = gdalfun.gdal_mem_ds(out_ds_config, bands=in_bands)
                if self.src_ds is not None:
                    for band in range(1, in_bands+1):
                        this_band = self.src_ds.GetRasterBand(band)
                        this_band.WriteArray(warp_ds.GetRasterBand(band).ReadAsArray(*srcwin))
                        self.src_ds.FlushCache()

                warp_ds = None
                utils.remove_glob(tmp_warp)
        else:
            if self.open_options:
                self.src_ds = gdal.OpenEx(self.fn, open_options=self.open_options)
                if self.src_ds is None:
                    if self.verbose:
                        utils.echo_warning_msg('could not open file using open options {}'.format(self.open_options))
                        
                    self.src_ds = gdal.Open(self.fn)
            else:
                self.src_ds = gdal.Open(self.fn)

        self.src_dem_infos = gdalfun.gdal_infos(self.src_ds)
        if self.src_ds is None:
            if self.verbose:
                utils.echo_error_msg('could not load raster file {}'.format(self.fn))

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
        src_dem_region = regions.Region().from_geo_transform(self.src_dem_infos['geoT'], self.src_dem_infos['nx'], self.src_dem_infos['ny'])
            
        ## todo: always warp these to src_ds
        ## weight mask, each cell should have the corresponding weight
        ## weight_mask can either be a seperate gdal file or a band number
        ## corresponding to the appropriate band in src_ds
        if self.weight_mask is not None:
            if utils.int_or(self.weight_mask) is not None:
                weight_band = self.src_ds.GetRasterBand(int(self.weight_mask))
            elif os.path.exists(self.weight_mask): # some numbers now return true here (file-descriptors), check for int first!
                src_weight = gdalfun.sample_warp(self.weight_mask, None, src_dem_x_inc, src_dem_y_inc,
                                                 src_srs=self.aux_src_proj4, dst_srs=self.aux_dst_proj4,
                                                 src_region=src_dem_region, sample_alg=self.sample_alg,
                                                 ndv=ndv, verbose=self.verbose, co=["COMPRESS=DEFLATE", "TILED=YES"])[0]
                weight_band = src_weight.GetRasterBand(1)

            else:
                utils.echo_warning_msg('could not load weight mask {}'.format(self.weight_mask))
                weight_band = None

        ## uncertainty mask, each cell should have the corresponding uncertainty
        ## uncertainty_mask can either be a seperate gdal file or a band number
        ## corresponding to the appropriate band in src_ds
        if self.uncertainty_mask is not None:
            if utils.int_or(self.uncertainty_mask):
                uncertainty_band = self.src_ds.GetRasterBand(int(self.uncertainty_mask))
            elif os.path.exists(self.uncertainty_mask):
                src_uncertainty = gdalfun.sample_warp(self.uncertainty_mask, None, src_dem_x_inc, src_dem_y_inc,
                                                      src_srs=self.aux_src_proj4, dst_srs=self.aux_dst_proj4,
                                                      src_region=src_dem_region, sample_alg='bilinear',
                                                      ndv=ndv, verbose=self.verbose, co=["COMPRESS=DEFLATE", "TILED=YES"],)[0]
                uncertainty_band = src_uncertainty.GetRasterBand(1)
            else:
                utils.echo_warning_msg('could not load uncertainty mask {}'.format(self.uncertainty_mask))
                uncertainty_band = None

        ## uncertainty from the vertical transformation
        if self.trans_fn_unc is not None:
            #print(self.trans_fn_unc, src_dem_x_inc, src_dem_y_inc, self.aux_dst_proj4, src_dem_region, self.sample_alg, ndv)
            #utils.echo_msg(self.sample_alg)
            trans_uncertainty = gdalfun.sample_warp(self.trans_fn_unc, None, src_dem_x_inc, src_dem_y_inc,
                                                    src_srs='+proj=longlat +datum=WGS84 +ellps=WGS84', dst_srs=self.aux_dst_proj4,
                                                    src_region=src_dem_region, sample_alg='bilinear',
                                                    ndv=ndv, verbose=self.verbose, co=["COMPRESS=DEFLATE", "TILED=YES"])[0]

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

        ## parse through the data
        srcwin = self.get_srcwin(gt, self.src_ds.RasterXSize, self.src_ds.RasterYSize, node=self.node)
        for y in range(srcwin[1], (srcwin[1] + srcwin[3]), 1):
            band_data = band.ReadAsArray(srcwin[0], y, srcwin[2], 1).astype(float)
            if ndv is not None and not np.isnan(ndv):
                band_data[band_data == ndv] = np.nan

            if np.all(np.isnan(band_data)):
                continue

            ## weights
            if weight_band is not None:
                weight_data = weight_band.ReadAsArray(srcwin[0], y, srcwin[2], 1)
                weight_ndv = float(weight_band.GetNoDataValue())
                if not np.isnan(weight_ndv):
                    weight_data[weight_data==weight_ndv] = np.nan
            else:
                weight_data = np.ones(band_data.shape)

            ## uncertainty
            if uncertainty_band is not None:
                uncertainty_data = uncertainty_band.ReadAsArray(srcwin[0], y, srcwin[2], 1)
                uncertainty_ndv = float(uncertainty_band.GetNoDataValue())
                if not np.isnan(uncertainty_ndv):
                    uncertainty_data[uncertainty_data==uncertainty_ndv] = np.nan

                if self.trans_fn_unc is None:
                    uncertainty_data *= self.uncertainty_mask_to_meter
                    
            else:
                uncertainty_data = np.zeros(band_data.shape)

            ## convert grid array to points
            #utils.echo_msg(band_data[0].shape)
            if self.x_band is None and self.y_band is None:

                x_precision = 12#len(str(gt[0]).split('.')[-1])
                y_precision = 12#len(str(gt[3]).split('.')[-1])
                #utils.echo_msg(gt)
                while True:
                    geo_x_origin, geo_y_origin = utils._pixel2geo(srcwin[0], y, gt, node=self.node, x_precision=x_precision, y_precision=y_precision)
                    geo_x_end, geo_y_end = utils._pixel2geo(srcwin[0] + srcwin[2], y, gt, node='grid', x_precision=x_precision, y_precision=y_precision)
                    lon_array = np.arange(geo_x_origin, geo_x_end, gt[1])

                    if lon_array.shape == band_data[0].shape:
                        break
                    else:
                        if x_precision < 0 or y_precision < 0:
                            break

                        #if self.node == 'grid':
                        x_precision -= 1
                        y_precision -= 1
                        #else:
                        #    x_precision += 1
                        #    y_precision += 1

                #utils.echo_msg('{} {}'.format(x_precision, y_precision))
                #num = round((geo_x_end - geo_x_origin) / gt[1])
                #lon_array = np.linspace(geo_x_origin, geo_x_end, num)
                lat_array = np.zeros((lon_array.shape))
                lat_array[:] = geo_y_origin
                dataset = np.column_stack((lon_array, lat_array, band_data[0], weight_data[0], uncertainty_data[0]))

                try:
                    assert lon_array.shape == lat_array.shape
                    assert lon_array.shape == band_data[0].shape
                    dataset = np.column_stack((lon_array, lat_array, band_data[0], weight_data[0], uncertainty_data[0]))
                except Exception as e:
                    utils.echo_error_msg(e)
                    pass
            else:
                lon_band = self.src_ds.GetRasterBand(self.x_band)
                lon_array = lon_band.ReadAsArray(srcwin[0], y, srcwin[2], 1).astype(float)
                lon_array[np.isnan(band_data)] = np.nan

                lat_band = self.src_ds.GetRasterBand(self.x_band)
                lat_array = lon_band.ReadAsArray(srcwin[0], y, srcwin[2], 1).astype(float)
                lat_array[np.isnan(band_data)] = np.nan
                dataset = np.column_stack((lon_array[0], lat_array[0], band_data[0], weight_data[0], uncertainty_data[0]))

            points = np.rec.fromrecords(dataset, names='x, y, z, w, u')
            points =  points[~np.isnan(points['z'])]
            dataset = band_data = weight_data = uncertainty_data = lat_array = lon_array = None
            utils.remove_glob(tmp_elev_fn, tmp_unc_fn, tmp_weight_fn)            
            yield(points)
            
        src_uncertainty = src_weight = trans_uncertainty = self.src_ds = None
        
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

    def __init__(
            self, explode = False, force_vr = False, vr_resampled_grid = True, vr_strategy = 'MIN', **kwargs
    ):
        super().__init__(**kwargs)
        self.explode = explode
        self.force_vr = force_vr
        self.vr_resampled_grid = vr_resampled_grid
        self.vr_strategy = vr_strategy
        if self.src_srs is None:
            self.src_srs = self.init_srs()

    def init_srs(self):
        if self.src_srs is None:
            src_horz, src_vert = gdalfun.split_srs(gdalfun.gdal_get_srs(self.fn))
            if src_horz is None and src_vert is None:
                return(None)
            
            if src_vert is None:
                src_vert = '5866'

            self.src_srs = gdalfun.combine_epsgs(src_horz, src_vert, name='BAG Combined')
            return(self.src_srs)
        else:
            return(self.src_srs)
                    
    def generate_inf(self, callback=lambda: False):
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

                zr = src_ds.GetRasterBand(1).ComputeRasterMinMax() # bag band 1 is elevation
                
        this_region.zmin, this_region.zmax = zr[0], zr[1]
        self.infos.minmax = this_region.export_as_list(include_z=True)
        self.infos.wkt = this_region.export_as_wkt()
        self.infos.numpts = ds_infos['nb']
        #utils.echo_msg(self.infos)
        return(self.infos)

    def parse(self, resample=True):
        mt = gdal.Info(self.fn, format='json')['metadata']['']
        oo = []

        if self.data_region is not None and self.data_region.valid_p():
            oo.append('MINX={}'.format(self.data_region.xmin))
            oo.append('MAXX={}'.format(self.data_region.xmax))
            oo.append('MINY={}'.format(self.data_region.ymin))
            oo.append('MAXY={}'.format(self.data_region.ymax))
                        
        if ('HAS_SUPERGRIDS' in mt.keys() and mt['HAS_SUPERGRIDS'] == 'TRUE') \
           or self.force_vr \
           or 'MAX_RESOLUTION_X' in mt.keys() \
           or 'MAX_RESOLUTION_Y' in mt.keys():
            if self.explode:
                oo.append("MODE=LIST_SUPERGRIDS")
                src_ds = gdal.OpenEx(self.fn, open_options=oo)
                sub_datasets = src_ds.GetSubDatasets()
                src_ds = None

                with tqdm(
                        total=len(sub_datasets),
                        desc='parsing {} supergrids from BAG file {}'.format(len(sub_datasets), self.fn),
                        leave=self.verbose
                ) as pbar:                
                    for sub_dataset in sub_datasets:
                        pbar.update()
                        sub_ds = GDALFile(fn=sub_dataset[0], data_format=200, band_no=1, src_srs=self.src_srs, dst_srs=self.dst_srs,
                                          weight=self.weight, uncertainty=self.uncertainty, uncertainty_mask_to_meter=0.01, src_region=self.region,
                                          x_inc=self.x_inc, y_inc=self.y_inc, verbose=False, check_path=False, super_grid=True,
                                          uncertainty_mask=2, metadata=copy.deepcopy(self.metadata))
                        self.data_entries.append(sub_ds)
                        sub_ds.initialize()
                        for gdal_ds in sub_ds.parse():
                            yield(gdal_ds)

            elif self.vr_resampled_grid:
                oo.append("MODE=RESAMPLED_GRID")
                oo.append("RES_STRATEGY={}".format(self.vr_strategy))
                sub_ds = GDALFile(fn=self.fn, data_format=200, band_no=1, open_options=oo, src_srs=self.src_srs, dst_srs=self.dst_srs,
                                  weight=self.weight, uncertainty=self.uncertainty, src_region=self.region, x_inc=self.x_inc, y_inc=self.y_inc,
                                  verbose=self.verbose, uncertainty_mask=2, uncertainty_mask_to_meter=0.01, metadata=copy.deepcopy(self.metadata))
                                  #node='pixel')
                self.data_entries.append(sub_ds)
                sub_ds.initialize()
                for gdal_ds in sub_ds.parse():
                    yield(gdal_ds)
            else: # use vrbag.py
                tmp_bag_as_tif = utils.make_temp_fn('{}_tmp.tif'.format(utils.fn_basename2(self.fn)))
                sr_cell_size = None#self.x_inc * 111120 # scale cellsize to meters, todo: check if input is degress/meters/feet
                vrbag.interpolate_vr_bag(self.fn, tmp_bag_as_tif, self.cache_dir, sr_cell_size=sr_cell_size, use_blocks=True, method='linear', nodata=3.4028234663852886e+38)

                sub_ds = GDALFile(fn=tmp_bag_as_tif, data_format=200, band_no=1, src_srs=self.src_srs, dst_srs=self.dst_srs, weight=self.weight,
                                  uncertainty=self.uncertainty, src_region=self.region, x_inc=self.x_inc, y_inc=self.y_inc, verbose=self.verbose,
                                  uncertainty_mask=2, uncertainty_mask_to_meter=0.01, metadata=copy.deepcopy(self.metadata))
                self.data_entries.append(sub_ds)
                sub_ds.initialize()
                for gdal_ds in sub_ds.parse():
                    yield(gdal_ds)
                
        else:
            sub_ds = GDALFile(fn=self.fn, data_format=200, band_no=1, src_srs=self.src_srs, dst_srs=self.dst_srs, weight=self.weight,
                              uncertainty=self.uncertainty, src_region=self.region, x_inc=self.x_inc, y_inc=self.y_inc, verbose=self.verbose,
                              uncertainty_mask=2, uncertainty_mask_to_meter=0.01, metadata=copy.deepcopy(self.metadata))
            self.data_entries.append(sub_ds)
            sub_ds.initialize()
            for gdal_ds in sub_ds.parse():
                yield(gdal_ds)

## NASA SWOT Data class (hdf5)
## uses h5py
class SWOTFile(ElevationDataset):
    """NASA SWOT Data super class

    Uses h5py to parse data. Make a subclass of this to 
    process various types of SWOT data
    """
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        
    def _init_h5File(self, short_name='L2_HR_PIXC'):
        src_h5 = None
        try:
            src_h5 = h5.File(self.fn, 'r')
            if src_h5 is not None:
                if 'short_name' in src_h5.attrs.keys():
                    if src_h5.attrs['short_name'] != short_name.encode('utf-8'):
                        utils.echo_error_msg('{} does not appear to be a SWOT {} file'.format(self.fn, short_name))
                        self._close_h5File(src_h5)
                else:
                    utils.echo_error_msg('{} does not appear to be a SWOT file'.format(self.fn))
                    self._close_h5File(src_h5)
                    
        except Exception as e:
            utils.echo_error_msg(e)

        return(src_h5)
            
    def _close_h5File(self, src_h5):
        if src_h5 is not None:
            src_h5.close()

    def _get_var_arr(self, src_h5, var_path):
        return(src_h5['/{}'.format(var_path)][...,])

# class IceSatFile(SWOTFile):
#     def __init__(self, group = 'pixel_cloud', var = 'height', apply_geoid = True, classes = None, classes_qual = None,
#                  anc_classes = None, remove_class_flags = False, **kwargs):
#         super().__init__(**kwargs)
#         self.group = group
#         self.var = var
#         self.apply_geoid = apply_geoid
#         self.classes = [int(x) for x in classes.split('/')] if classes is not None else []
#         self.classes_qual = [int(x) for x in classes_qual.split('/')] if classes_qual is not None else []
#         self.anc_classes = [int(x) for x in anc_classes.split('/')] if anc_classes is not None else []
#         self.remove_class_flags = remove_class_flags
#         #if self.remove_class_flags:
#         #    self.classes_qual = [1, 2, 4, 8, 16, 2048, 8192, 16384, 32768, 262144, 524288, 134217728, 536870912, 1073741824, 2147483648]

#     def yield_ds(self):
#         src_h5 = self._init_h5File(short_name='L2_HR_PIXC')
#         src_h5_vec = None
        
#         #if self.pixc_vec is not None:
#         #    src_h5_vec = self._init_h5File(short_name='L2_HR_PIXCVec')
        
#         if src_h5 is not None:
#             latitude = self._get_var_arr(src_h5, '{}/latitude'.format(self.group))
#             longitude = self._get_var_arr(src_h5, '{}/longitude'.format(self.group))
#             var_data = self._get_var_arr(src_h5, '{}/{}'.format(self.group, self.var))
#             if self.apply_geoid:
#                 geoid_data = self._get_var_arr(src_h5, '{}/geoid'.format(self.group))
#                 out_data = var_data - geoid_data
#             else:
#                 out_data = var_data
                
#             dataset = np.column_stack((longitude, latitude, out_data))
#             points = np.rec.fromrecords(dataset, names='x, y, z')
#             #points = points[points['z'] != 9.96921e+36]

#             ## Classification Filter
#             if len(self.classes) > 0:
#                 class_data = self._get_var_arr(src_h5, '{}/classification'.format(self.group))
#                 points = points[(np.isin(class_data, self.classes))]

#                 ## Classification Quality Filter
#                 if self.remove_class_flags:
#                     class_qual_data = self._get_var_arr(src_h5, '{}/classification_qual'.format(self.group))
#                     class_qual_data = class_qual_data[(np.isin(class_data, self.classes))]
#                     points = points[class_qual_data == 0]
                                   
#                 elif len(self.classes_qual) > 0:
#                     class_qual_data = self._get_var_arr(src_h5, '{}/classification_qual'.format(self.group))
#                     class_qual_data = class_qual_data[(np.isin(class_data, self.classes))]
#                     points = points[(~np.isin(class_qual_data, self.classes_qual))]
                
#             ## Ancilliary Classification Filter
#             if len(self.anc_classes) > 0:
#                 anc_class_data = self._get_var_arr(src_h5, '{}/ancillary_surface_classification_flag'.format(self.group))
#                 points = points[(np.isin(anc_class_data, self.anc_classes))]
                
#             points = points[points['z'] != -9.969209968386869e+36]
#             self._close_h5File(src_h5)
#             self._close_h5File(src_h5_vec)

#             #tmp_points = points
#             #while len(points) > 0:
#             #tmp_points = points[:1000000]
#             #points = points[1000000:]
#             #utils.echo_msg(len(points))
#             #yield(tmp_points)
#             yield(points)
    
class SWOT_PIXC(SWOTFile):
    """NASA SWOT PIXC data file.

    Extract data from a SWOT PIXC file.

    classes: 1UB, 2UB, 3UB, 4UB, 5UB, 6UB, 7UB
    "land, land_near_water, water_near_land, open_water, dark_water, low_coh_water_near_land, open_low_coh_water"

    classes_qual: 1U, 2U, 4U, 8U, 16U, 2048U, 8192U, 16384U, 32768U, 262144U, 524288U, 134217728U, 536870912U, 1073741824U, 2147483648U
    "no_coherent_gain power_close_to_noise_floor detected_water_but_no_prior_water detected_water_but_bright_land water_false_detection_rate_suspect coherent_power_suspect tvp_suspect sc_event_suspect small_karin_gap in_air_pixel_degraded specular_ringing_degraded coherent_power_bad tvp_bad sc_event_bad large_karin_gap"

    anc_classes: 0UB, 1UB, 2UB, 3UB, 4UB, 5UB, 6UB 
    "open_ocean land continental_water aquatic_vegetation continental_ice_snow floating_ice salted_basin"
  		
    """
    
    def __init__(self, group = 'pixel_cloud', var = 'height', apply_geoid = True, classes = None, classes_qual = None,
                 anc_classes = None, remove_class_flags = False, **kwargs):
        super().__init__(**kwargs)
        self.group = group
        self.var = var
        self.apply_geoid = apply_geoid
        self.classes = [int(x) for x in classes.split('/')] if classes is not None else []
        self.classes_qual = [int(x) for x in classes_qual.split('/')] if classes_qual is not None else []
        self.anc_classes = [int(x) for x in anc_classes.split('/')] if anc_classes is not None else []
        self.remove_class_flags = remove_class_flags
        #if self.remove_class_flags:
        #    self.classes_qual = [1, 2, 4, 8, 16, 2048, 8192, 16384, 32768, 262144, 524288, 134217728, 536870912, 1073741824, 2147483648]

    def yield_ds(self):
        src_h5 = self._init_h5File(short_name='L2_HR_PIXC')
        src_h5_vec = None
        
        #if self.pixc_vec is not None:
        #    src_h5_vec = self._init_h5File(short_name='L2_HR_PIXCVec')
        
        if src_h5 is not None:
            latitude = self._get_var_arr(src_h5, '{}/latitude'.format(self.group))
            longitude = self._get_var_arr(src_h5, '{}/longitude'.format(self.group))
            var_data = self._get_var_arr(src_h5, '{}/{}'.format(self.group, self.var))
            if self.apply_geoid:
                geoid_data = self._get_var_arr(src_h5, '{}/geoid'.format(self.group))
                out_data = var_data - geoid_data
            else:
                out_data = var_data
                
            dataset = np.column_stack((longitude, latitude, out_data))
            points = np.rec.fromrecords(dataset, names='x, y, z')
            #points = points[points['z'] != 9.96921e+36]

            ## Classification Filter
            if len(self.classes) > 0:
                class_data = self._get_var_arr(src_h5, '{}/classification'.format(self.group))
                points = points[(np.isin(class_data, self.classes))]

                ## Classification Quality Filter
                if self.remove_class_flags:
                    class_qual_data = self._get_var_arr(src_h5, '{}/classification_qual'.format(self.group))
                    class_qual_data = class_qual_data[(np.isin(class_data, self.classes))]
                    points = points[class_qual_data == 0]
                                   
                elif len(self.classes_qual) > 0:
                    class_qual_data = self._get_var_arr(src_h5, '{}/classification_qual'.format(self.group))
                    class_qual_data = class_qual_data[(np.isin(class_data, self.classes))]
                    points = points[(~np.isin(class_qual_data, self.classes_qual))]
                
            ## Ancilliary Classification Filter
            if len(self.anc_classes) > 0:
                anc_class_data = self._get_var_arr(src_h5, '{}/ancillary_surface_classification_flag'.format(self.group))
                points = points[(np.isin(anc_class_data, self.anc_classes))]
                
            points = points[points['z'] != -9.969209968386869e+36]
            self._close_h5File(src_h5)
            self._close_h5File(src_h5_vec)

            #tmp_points = points
            #while len(points) > 0:
            #tmp_points = points[:1000000]
            #points = points[1000000:]
            #utils.echo_msg(len(points))
            #yield(tmp_points)
            yield(points)

class SWOT_HR_Raster(ElevationDataset):
    """NASA SWOT HR_Raster data file.

    Extract data from a SWOT HR_Raster file.
    """
        
    def __init__(self, data_set='wse', **kwargs):
        super().__init__(**kwargs)
        self.data_set = data_set

    def parse(self):
        src_ds = gdal.Open(self.fn)
        if src_ds is not None:
            sub_datasets = src_ds.GetSubDatasets()

            idx = 2
            if utils.int_or(self.data_set) is not None:
                idx = utils.int_or(self.data_set)
            else:
                for j, sd in enumerate(sub_datasets):
                    _name = sd[0].split(':')[-1]
                    if self.data_set == _name:
                        idx = j
                        break

            src_ds = None
            src_srs = gdalfun.gdal_get_srs(sub_datasets[idx][0])
            if self.data_set == 'wse':
                src_srs = gdalfun.combine_epsgs(src_srs, '3855', name='SWOT Combined')

            sub_ds = GDALFile(fn=sub_datasets[idx][0], data_format=200, src_srs=src_srs, dst_srs=self.dst_srs,
                              weight=self.weight, uncertainty=self.uncertainty, src_region=self.region,
                              x_inc=self.x_inc, y_inc=self.y_inc, verbose=True, check_path=False,
                              node='grid', metadata=copy.deepcopy(self.metadata))
            
            self.data_entries.append(sub_ds)
            sub_ds.initialize()
            for gdal_ds in sub_ds.parse():
                yield(gdal_ds)
                
class MBSParser(ElevationDataset):
    """providing an mbsystem parser

    Process MB-System supported multibeam data files.
    See MB-System for more information regarding supported
    file formats, etc.
    
    generate_inf - generate an inf file for the MBS data
    yield_xyz - yield the MBS data as xyz
    yield_array - yield the MBS data as an array
 
    -----------
    Parameters:

    mb_fmt=[]
    mb_exclude=[]
    """

    def __init__(
            self, mb_fmt = None, mb_exclude = 'A', want_mbgrid = False,
            want_binned = False, **kwargs
    ):
        super().__init__(**kwargs)
        self.mb_fmt = mb_fmt
        self.mb_exclude = mb_exclude
        self.want_mbgrid = want_mbgrid
        self.want_binned = want_binned
             
    def inf_parse(self):
        self.infos.minmax = [0,0,0,0,0,0]
        this_row = 0
        xinc = 0
        yinc = 0
        dims = []
        inf_fn = '{}.inf'.format(self.fn) if self.fn.split('.')[-1] != 'inf' else self.fn
        with open(inf_fn) as iob:
            for il in iob:
                til = il.split()
                if len(til) > 1:

                    if til[0] == 'Swath':
                        if til[2] == 'File:':
                            self.infos.name = til[3]

                    if ' '.join(til[:-1]) == 'MBIO Data Format ID:':
                        self.mb_fmt = til[-1]
                            
                    if til[0] == 'Number':
                        if til[2] == 'Records:':
                            self.infos.numpts = utils.int_or(til[3])
                            
                    if til[0] == 'Minimum':
                        if til[1] == 'Longitude:':
                            self.infos.minmax[0] = utils.float_or(til[2])
                            self.infos.minmax[1] = utils.float_or(til[5])
                        elif til[1] == 'Latitude:':
                            self.infos.minmax[2] = utils.float_or(til[2])
                            self.infos.minmax[3] = utils.float_or(til[5])
                        elif til[1] == 'Depth:':
                            self.infos.minmax[4] = utils.float_or(til[5]) * -1
                            self.infos.minmax[5] = utils.float_or(til[2]) * -1
                            
                    if til[0] == 'CM':
                        if til[1] == 'dimensions:':
                            dims = [utils.int_or(til[2]), utils.int_or(til[3])]
                            cm_array = np.zeros((dims[0], dims[1]))
                            
                    if til[0] == 'CM:':
                        for j in range(0, dims[0]):
                            cm_array[this_row][j] = utils.int_or(til[j+1])
                        this_row += 1

        mbs_region = regions.Region().from_list(self.infos.minmax)

        if len(dims) > 0:
            xinc = (mbs_region.xmax - mbs_region.xmin) / dims[0]
            yinc = (mbs_region.ymin - mbs_region.ymax) / dims[1]
        
        if abs(xinc) > 0 and abs(yinc) > 0:
            xcount, ycount, dst_gt = mbs_region.geo_transform(
                x_inc=xinc, y_inc=yinc
            )
            ds_config = {'nx': dims[0], 'ny': dims[1], 'nb': dims[1] * dims[0],
                         'geoT': dst_gt, 'proj': gdalfun.osr_wkt(self.src_srs),
                         'dt': gdal.GDT_Float32, 'ndv': 0, 'fmt': 'GTiff'}
            driver = gdal.GetDriverByName('MEM')
            ds = driver.Create(
                'tmp', ds_config['nx'], ds_config['ny'], 1, ds_config['dt']
            )
            ds.SetGeoTransform(ds_config['geoT'])
            if ds_config['proj'] is not None:
                ds.SetProjection(ds_config['proj'])
                
            ds_band = ds.GetRasterBand(1)
            ds_band.SetNoDataValue(ds_config['ndv'])
            ds_band.WriteArray(cm_array)
            tmp_ds = ogr.GetDriverByName('Memory').CreateDataSource('tmp_poly')
            tmp_layer = tmp_ds.CreateLayer('tmp_poly', None, ogr.wkbMultiPolygon)
            tmp_layer.CreateField(ogr.FieldDefn('DN', ogr.OFTInteger))
            gdal.Polygonize(ds_band, ds_band, tmp_layer, 0)
            multi = ogr.Geometry(ogr.wkbMultiPolygon)
            for feat in tmp_layer:
                feat.geometry().CloseRings()
                wkt = feat.geometry().ExportToWkt()
                multi.AddGeometryDirectly(ogr.CreateGeometryFromWkt(wkt))
                
            wkt = multi.ExportToWkt()
            tmp_ds = ds = None
        else:
            wkt = mbs_region.export_as_wkt()

        self.infos.wkt = wkt
        return(self)

    def parse_(self):
        """use mblist to convert data to xyz then set the dataset as xyz and use that to process..."""
        
        mb_fn = os.path.join(self.fn)
        if self.region is not None:
            mb_region = self.region.copy()
            mb_region.buffer(pct=25)
        else:
            mb_region = None

        out_mb = utils.make_temp_fn('{}_.xyz'.format(utils.fn_basename2(mb_fn)), self.cache_dir)
        out, status = utils.run_cmd('mblist -M{}{} -OXYZ -I{} > {}'.format(
            self.mb_exclude, ' {}'.format(
                mb_region.format('gmt') if mb_region is not None else ''
            ), mb_fn, out_mb
        ), verbose=False)

        if status == 0:
            data_set = DatasetFactory(
                mod = out_mb, data_format = 168, weight=self.weight, uncertainty=self.uncertainty, parent=self, src_region=self.region,
                invert_region=self.invert_region, metadata=copy.deepcopy(self.metadata), src_srs=self.src_srs, dst_srs=self.dst_srs, mask=self.mask,
                x_inc=self.x_inc, y_inc=self.y_inc, sample_alg=self.sample_alg, cache_dir=self.cache_dir,
                verbose=self.verbose
            )._acquire_module().initialize()
            
            for ds in data_set.parse(): # fill self.data_entries with each dataset for use outside the yield.
                self.data_entries.append(ds) 
                yield(ds)

        utils.remove_glob('{}*'.format(out_mb))

    def yield_mbgrid_ds(self):
        """process the data through mbgrid and use GDALFile to further process the gridded data"""
        
        mb_fn = os.path.join(self.fn)
        with open('_mb_grid_tmp.datalist', 'w') as tmp_dl:
            tmp_dl.write(
                '{} {} {}\n'.format(
                    self.fn, self.mb_fmt if self.mb_fmt is not None else '', self.weight if self.mb_fmt is not None else ''
                )
            )

        ofn = '_'.join(os.path.basename(mb_fn).split('.')[:-1])
        try:
            utils.run_cmd(
                'mbgrid -I_mb_grid_tmp.datalist {} -E{}/{}/degrees! -O{} -A2 -F1 -C10/1 -S0 -T35'.format(
                    self.data_region.format('gmt'), self.x_inc, self.y_inc, ofn
                ), verbose=True
            )

            gdalfun.gdal2gdal('{}.grd'.format(ofn))
            utils.remove_glob('_mb_grid_tmp.datalist', '{}.cmd'.format(ofn), '{}.mb-1'.format(ofn), '{}.grd*'.format(ofn))
            mbs_ds = GDALFile(fn='{}.tif'.format(ofn), data_format=200, src_srs=self.src_srs, dst_srs=self.dst_srs,
                              weight=self.weight, x_inc=self.x_inc, y_inc=self.y_inc, sample_alg=self.sample_alg,
                              src_region=self.region, verbose=self.verbose, metadata=copy.deepcopy(self.metadata))

            yield(mbs_ds)
            utils.remove_glob('{}.tif*'.format(ofn))
        except:
            pass

    def yield_mblist_ds(self):
        """use mblist to process the multibeam data"""
        
        mb_fn = os.path.join(self.fn)
        xs = []
        ys = []
        zs = []
        ws = []
        us = []
        if self.region is not None:
            mb_region = self.region.copy()
            mb_region.buffer(pct=25)
        else:
            mb_region = None

        for line in utils.yield_cmd(
                'mblist -M{}{} -OXYZDAGgFPpRrS -I{}'.format(
                    self.mb_exclude, ' {}'.format(
                        mb_region.format('gmt') if mb_region is not None else ''
                    ), mb_fn
                ),
                verbose=False,
        ):
            this_line = [float(x) for x in line.strip().split('\t')]
            x = this_line[0]
            y = this_line[1]
            z = this_line[2]
            crosstrack_distance = this_line[3]
            crosstrack_slope = this_line[4]
            flat_bottom_grazing_angle = this_line[5]
            seafloor_grazing_angle = this_line[6]
            beamflag = this_line[7]
            pitch = this_line[8]
            draft = this_line[9]
            roll = this_line[10]
            heave = this_line[11]
            speed = this_line[12]

            #this_xyz = xyzfun.XYZPoint().from_string(line, delim='\t')
            if int(beamflag) == 0:# and abs(this_line[4]) < .15:
                u_depth = ((2+(0.02*(z*-1)))*0.51)
                #u_depth = 0
                #u_depth = math.sqrt(1 + ((.023 * (z * -1))**2))
                u_cd = math.sqrt(1 + ((.023 * abs(crosstrack_distance))**2))
                #u_cd = ((2+(0.02*abs(crosstrack_distance)))*0.51) ## find better alg.
                #u_cd = 0
                u = math.sqrt(u_depth**2 + u_cd**2)
                xs.append(x)
                ys.append(y)
                zs.append(z)
                ws.append(1)
                us.append(u)

        if len(xs) > 0:
            mb_points = np.column_stack((xs, ys, zs, ws, us))
            xs = ys = zs = ws = us = None
            mb_points = np.rec.fromrecords(mb_points, names='x, y, z, w, u')

            if self.want_binned:
                mb_points = self.bin_z_points(mb_points)

            if mb_points is not None:
                yield(mb_points)

    def yield_mblist2_ds(self):
        """use mblist to process the multibeam data"""
        
        mb_fn = os.path.join(self.fn)
        if self.region is None or self.data_region is None:
            self.want_mbgrid = False

        if self.region is not None:
            mb_region = self.region.copy()
            mb_region.buffer(pct=25)
        else:
            mb_region = None

        mb_points = [[float(x) for x in l.strip().split('\t')] for l in utils.yield_cmd(
                'mblist -M{}{} -OXYZ -I{}'.format(
                    self.mb_exclude, ' {}'.format(
                        mb_region.format('gmt') if mb_region is not None else ''
                    ), mb_fn#, '-F{}'.format(self.mb_fmt) if self.mb_fmt is not None else ''
                ),
                verbose=False,
        )]

        if len(mb_points) > 0:
            mb_points = np.rec.fromrecords(mb_points, names='x, y, z')
        else:
            mb_points = None

        if self.want_binned:
            mb_points = self.bin_z_points(mb_points)

        if mb_points is not None:
            yield(mb_points)
                
    def yield_ds(self):        
        mb_fn = os.path.join(self.fn)
        if self.region is None or self.data_region is None:
            self.want_mbgrid = False

        ## update want_mbgrid to output points!
        if self.want_mbgrid and (self.x_inc is not None and self.y_inc is not None):
            for pts in self.yield_mbgrid_ds():
                yield(pts)
        else:
            for pts in self.yield_mblist_ds():
                yield(pts)
                        
class OGRFile(ElevationDataset):
    """providing an OGR 3D point dataset parser.

    Useful for data such as S-57, ENC, E-Hydro, Etc.

    -----------
    Parameters:

    ogr_layer (str/int): the OGR layer containing elevation data
    elev_field (str): the field containing the z values
    weight_field (str): the field containing weight values
    uncertainty_field (str): the field containing uncertainty_values
    """

    _known_layer_names = ['SOUNDG', 'SurveyPoint_HD', 'SurveyPoint']
    _known_elev_fields = ['Elevation', 'elev', 'z', 'height', 'depth', 'topography', 'surveyPointElev', 'Z_depth', 'Z_height']
    
    def __init__(self, ogr_layer = None, elev_field = None, weight_field = None, uncertainty_field = None,
                 z_scale = None, **kwargs):
        super().__init__(**kwargs)
        self.ogr_layer = ogr_layer
        self.elev_field = elev_field
        self.weight_field = weight_field
        self.uncertainty_field = uncertainty_field
        self.z_scale = utils.float_or(z_scale)

    def find_elevation_layer(self, ds_ogr):
        for l in self._known_layer_names:
            test_layer = ds_ogr.GetLayerByName(l)
            if test_layer is not None:
                return((l, test_layer))
        return(None, None)

    def find_elevation_field(self, ogr_feature):
        for f in self._known_elev_fields:
            test_field = ogr_feature.GetField(f)
            if test_field is not None:
                return((l, test_field))
        return(None, None)
    
    def yield_ds(self):
        ds_ogr = ogr.Open(self.fn)
        count = 0
        #utils.echo_msg(self.ogr_layer)
        if ds_ogr is not None:
            layer_name = None
            if self.ogr_layer is None:
                layer_name, layer_s = self.find_elevation_layer(ds_ogr)
            elif utils.str_or(self.ogr_layer) is not None:
                layer_name = self.ogr_layer
                layer_s = ds_ogr.GetLayer(str(self.ogr_layer))
                if layer_s is None:
                    layer_name, layer_s = self.find_elevation_layer(ds_ogr)
            elif utils.int_or(self.ogr_layer) is not None:
                layer_s = ds_ogr.GetLayer(self.ogr_layer)
            else:
                layer_s = ds_ogr.GetLayer()

            if layer_s is not None:
                if self.region is not None:
                    layer_s.SetSpatialFilter(
                        self.region.export_as_geom() if self.transformer is None else self.trans_region.export_as_geom()
                    )

                for f in layer_s:
                    geom = f.GetGeometryRef()
                    g = json.loads(geom.ExportToJson())
                    #utils.echo_msg(g)
                    xyzs = g['coordinates']

                    if not geom.GetGeometryName() == 'MULTIPOINT':
                        xyzs = [xyzs]

                    out_xyzs = []
                    for xyz in xyzs:

                        if not geom.Is3D():
                            if self.elev_field is None:
                                self.elev_field = self.find_elevation_field(f)

                            elev = utils.float_or(f.GetField(self.elev_field))
                            if elev is not None:
                                
                                if isinstance(xyz[0], list):
                                    for x in xyz:
                                        for xx in x:
                                            xx.append(elev)
                                else:
                                    xyz.append(elev)
                            else:
                                continue

                        else:
                            if self.elev_field is not None:
                                elev = utils.float_or(f.GetField(self.elev_field))
                                if isinstance(xyz[0], list):
                                    for x in xyz:
                                        for xx in x:
                                            xx[2] = elev

                    if isinstance(xyzs[0], list):
                        for x in xyzs:
                            if isinstance(x[0], list):
                                for xx in x:
                                    points = np.rec.fromrecords(xx, names='x, y, z')
                            else:
                                #utils.echo_msg(x)
                                points = np.rec.fromrecords([x], names='x, y, z')
                                #utils.echo_msg(points)

                    #points = np.rec.fromrecords(out_xyzs, names='x, y, z')
                    if self.z_scale is not None:
                        points['z'] *= self.z_scale

                    yield(points)
                            
            ds_ogr = layer_s = None

class Scratch(ElevationDataset):
    """Scratch Dataset

    Process a python list of valid dataset entries
    """
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.metadata['name'] = 'scratch'

    def generate_inf(self, callback=lambda: False):
        """generate and return the infos blob of the datalist.

        datasets parsed from the datalist may have variable srs, so
        we want to transform those to whatever self.dst_srs is, if it
        exists, in order to properly fill the datalist/json srs and regions...
        """

        _region = self.region
        self.region = None
        out_regions = []
        out_srs = []
        for entry in self.parse():
            if self.verbose:
                callback()

            entry_minmax = entry.infos.minmax
            ## entry has an srs and dst_srs is set, so lets transform the region to suit
            if entry.src_srs is not None:
                out_srs.append(entry.src_srs)
                if self.dst_srs is not None:
                    entry_region = regions.Region().from_list(entry_minmax)
                    if entry_region.valid_p():
                        entry_region.src_srs = entry.src_srs
                        entry_region.warp(self.dst_srs)
                        entry_minmax = entry_region.export_as_list(include_z=True)
                        
            entry_region = regions.Region().from_list(entry_minmax)
            if entry_region.valid_p():
                out_regions.append(entry_region)
                self.infos.numpts += entry.infos.numpts

        ## merge all the gathered regions
        region_count = 0
        out_region = None
        for this_region in out_regions:
            if this_region.valid_p():
                if region_count == 0:
                    out_region = this_region
                    region_count += 1
                else:
                    out_region = regions.regions_merge(out_region, this_region)
                    
        if out_region is not None:
            self.infos.minmax = out_region.export_as_list(include_z=True)
            self.infos.wkt = out_region.export_as_wkt()

        ## set the epsg for the datalist
        if self.infos.src_srs is None:
            if self.src_srs is not None:
                self.infos.src_srs = self.src_srs
            elif out_srs:
                if all(x == out_srs[0] for x in out_srs):
                    self.infos.src_srs = out_srs[0]

                self.src_srs = self.infos.src_srs
                
        self.region = _region
        return(self.infos)
    
    def parse(self):
        if isinstance(self.fn, list):
            for this_ds in self.fn:
                if this_ds is not None and this_ds.valid_p(
                        fmts=DatasetFactory._modules[this_ds.data_format]['fmts']
                ):
                    this_ds.initialize()
                    for ds in this_ds.parse():
                        self.data_entries.append(ds) # fill self.data_entries with each dataset for use outside the yield.
                        yield(ds)
                        
class Datalist(ElevationDataset):
    """representing a datalist parser
    
    A datalist is an extended MB-System style datalist.
    
    Each datalist consists of datalist-entries, where a datalist entry has the following columns:
    path format weight uncertainty title source date data_type resolution hdatum vdatum url

    the datalist can contain datalist-entries to other datalist files, distributed across a
    file-system.

    see `cudem.dlim.datasets` for superclass ElevationDataset
    """

    _datalist_json_cols = ['path', 'format', 'weight', 'uncertainty', 'name', 'title', 'source',
                           'date', 'data_type', 'resolution', 'hdatum', 'vdatum',
                           'url', 'mod_args']

    def __init__(self, fmt=None, **kwargs):
        super().__init__(**kwargs)
        
    def _init_datalist_vector(self):
        """initialize the datalist geojson vector.

        this vector is used to quickly parse data from distributed files
        across a file-system. Each data file will have it's own entry in
        the geojson datalist vector, containing gathered metadata based on the
        source datalist.
        """

        #self.set_transform()
        self.dst_layer = '{}'.format(self.fn)
        self.dst_vector = self.dst_layer + '.json'

        utils.remove_glob('{}.json'.format(self.dst_layer))
        if self.src_srs is not None:
            gdalfun.osr_prj_file('{}.prj'.format(self.dst_layer), self.src_srs)
            
        self.ds = ogr.GetDriverByName('GeoJSON').CreateDataSource(self.dst_vector)
        if self.ds is not None: 
            self.layer = self.ds.CreateLayer(
                '{}'.format(self.dst_layer), None, ogr.wkbMultiPolygon
            )
            [self.layer.CreateField(
                ogr.FieldDefn('{}'.format(f), ogr.OFTString)
            ) for f in self._datalist_json_cols]
            
            [self.layer.SetFeature(feature) for feature in self.layer]
            return(0)
        else:
            self.layer = None
            return(-1)

    def _create_entry_feature(self, entry, entry_region):
        """create a datalist entry feature and insert it into the datalist-vector geojson

        -----------
        Parameters:

        entry - the datalist entry object
        entry_region - the region of the datalist entry
        """

        entry_path = os.path.abspath(entry.fn) if not entry.remote else entry.fn
        entry_fields = [entry_path, entry.data_format, entry.weight, entry.uncertainty,
                        entry.metadata['name'], entry.metadata['title'], entry.metadata['source'],
                        entry.metadata['date'], entry.metadata['data_type'], entry.metadata['resolution'],
                        entry.metadata['hdatum'], entry.metadata['vdatum'], entry.metadata['url'],
                        utils.dict2args(entry.params['mod_args'])]
        dst_defn = self.layer.GetLayerDefn()
        entry_geom = ogr.CreateGeometryFromWkt(entry_region.export_as_wkt())
        out_feat = ogr.Feature(dst_defn)
        out_feat.SetGeometry(entry_geom)
        for i, f in enumerate(self._datalist_json_cols):
            out_feat.SetField(f, entry_fields[i])
            
        self.layer.CreateFeature(out_feat)

    def generate_inf(self, callback=lambda: False):
        """generate and return the infos blob of the datalist.

        datasets parsed from the datalist may have variable srs, so
        we want to transform those to whatever self.dst_srs is, if it
        exists, in order to properly fill the datalist/json srs and regions...
        """

        _region = self.region
        self.region = None
        out_regions = []
        out_srs = []
        self.infos.file_hash = self.infos.generate_hash()

        ## attempt to generate a datalist-vector geojson and
        ## if successful, fill it wil the datalist entries, using `parse`
        if self._init_datalist_vector() == 0:
            for entry in self.parse():
                #utils.echo_msg_bold(entry.mask)
                if self.verbose:
                    callback()

                entry_minmax = entry.infos.minmax
                #utils.echo_msg(entry.params)
                if entry.mask is not None: ## add all duplicat params???
                    entry.params['mod_args'] = {'mask': entry.mask}
                    
                ## entry has an srs and dst_srs is set, so lets transform the region to suit
                if entry.src_srs is not None:
                    out_srs.append(entry.src_srs)
                    if self.dst_srs is not None:
                        entry_region = regions.Region().from_list(entry_minmax)
                        entry_region.src_srs = entry.src_srs
                        #entry_region.warp(gdalfun.epsg_from_input(self.dst_srs)[0])
                        entry_region.warp(self.dst_srs)
                        entry_minmax = entry_region.export_as_list(include_z=True)

                ## create the feature for the geojson
                entry_region = regions.Region().from_list(entry_minmax)
                if entry_region.valid_p():
                    self._create_entry_feature(entry, entry_region)
                    out_regions.append(entry_region)
                    self.infos.numpts += entry.infos.numpts

            self.ds = self.layer = None # close the geojson ogr dataset

        else:
            utils.echo_warning_msg('could not initialize datalist vector')
            return(self.infos)

        ## merge all the gathered regions
        region_count = 0
        out_region = None
        for this_region in out_regions:
            if this_region.valid_p():
                if region_count == 0:
                    out_region = this_region
                    region_count += 1
                else:
                    out_region = regions.regions_merge(out_region, this_region)
                    
        if out_region is not None:
            self.infos.minmax = out_region.export_as_list(include_z=True)
            self.infos.wkt = out_region.export_as_wkt()

        ## set the epsg for the datalist
        if self.infos.src_srs is None:
            if self.src_srs is not None:
                self.infos.src_srs = self.src_srs
            elif out_srs:
                if all(x == out_srs[0] for x in out_srs):
                    self.infos.src_srs = out_srs[0]

                self.src_srs = self.infos.src_srs

        self.region = _region
        return(self.infos)
        
    def parse(self):
        """parse the datalist using the datalist-vector geojson.

        Quickly find data in a given region using the datalist-vector. The datalist-vector
        must have been previously generated using `parse`. If the datlist-vector is not found
        will fall-back to `parse` and generate a new datalist-vector geojson.

        -------
        Yields:
        dataset object of each dataset found in the datalist-vector
        """

        ## check for the datalist-vector geojson
        status = 0
        count = 0
        
        ## user input to re-gerenate json...?
        if os.path.exists('{}.json'.format(self.fn)):
            driver = ogr.GetDriverByName('GeoJSON')
            dl_ds = driver.Open('{}.json'.format(self.fn))
            if dl_ds is None:
                utils.echo_error_msg('could not open {}.json'.format(self.fn))
                status = -1
        else:
            status = -1

        ## parse the datalist-vector geojson and yield the results
        if status != -1:
            dl_layer = dl_ds.GetLayer()
            ldefn = dl_layer.GetLayerDefn()
            _boundsGeom = None
            if self.region is not None:
                _boundsGeom = self.region.export_as_geom() if self.transformer is None else self.trans_region.export_as_geom()

            dl_layer.SetSpatialFilter(_boundsGeom)
            count = len(dl_layer)
            with tqdm(
                    total=len(dl_layer),
                    desc='parsing {} datasets from datalist json {}.json @ {}'.format(count, self.fn, self.weight),
                    leave=self.verbose
            ) as pbar:
                for l,feat in enumerate(dl_layer):
                    pbar.update()
                    ## filter by input source region extras (weight/uncertainty)
                    if self.region is not None:
                        w_region = self.region.w_region()
                        if w_region[0] is not None:
                            if float(feat.GetField('weight')) < w_region[0]:
                                continue

                        if w_region[1] is not None:
                            if float(feat.GetField('weight')) > w_region[1]:
                                continue

                        u_region = self.region.u_region()
                        if u_region[0] is not None:
                            if float(feat.GetField('uncertainty')) < u_region[0]:
                                continue

                        if u_region[1] is not None:
                            if float(feat.GetField('uncertainty')) > u_region[1]:
                                continue

                    ## extract the module arguments from the datalist-vector
                    try:
                        ds_args = feat.GetField('mod_args')
                        data_set_args = utils.args2dict(list(ds_args.split(':')), {})
                        for kpam, kval in data_set_args.items():
                            if kpam in self.__dict__:
                                self.__dict__[kpam] = kval
                                del data_set_args[kpam]
                    except:
                        data_set_args = {}

                    ## update existing metadata
                    md = copy.deepcopy(self.metadata)
                    for key in self.metadata.keys():
                        md[key] = feat.GetField(key)

                    ## generate the dataset object to yield
                    data_mod = '"{}" {} {} {}'.format(feat.GetField('path'), feat.GetField('format'),
                                                      feat.GetField('weight'), feat.GetField('uncertainty'))
                    #utils.echo_msg_bold(self.src_srs)
                    data_set = DatasetFactory(mod = data_mod, weight=self.weight, uncertainty=self.uncertainty, parent=self, src_region=self.region,
                                              invert_region=self.invert_region, metadata=md, src_srs=self.src_srs, dst_srs=self.dst_srs, mask=self.mask,
                                              x_inc=self.x_inc, y_inc=self.y_inc, sample_alg=self.sample_alg, cache_dir=self.cache_dir,
                                              stack_fltrs=self.stack_fltrs, pnt_fltrs=self.pnt_fltrs, verbose=self.verbose, **data_set_args)._acquire_module()
                    if data_set is not None and data_set.valid_p(
                            fmts=DatasetFactory._modules[data_set.data_format]['fmts']
                    ):
                        data_set.initialize()
                        for ds in data_set.parse(): # fill self.data_entries with each dataset for use outside the yield.
                            self.data_entries.append(ds) 
                            yield(ds)

            dl_ds = dl_layer = None
                
        else:
            ## failed to find/open the datalist-vector geojson, so run `parse` instead and
            ## generate one for future use...
            utils.echo_warning_msg(
                'could not load datalist-vector json {}.json, falling back to parse, generate a json file for the datalist using `dlim -i`'.format(self.fn)
            )
            for ds in self.parse_no_json():
                yield(ds)
                                        
    def parse_no_json(self):
        """parse the datalist file.

        -------
        Yields:
        dataset object of each dataset found in the datalist        
        """

        status = 0
        if os.path.exists(self.fn):
            with open(self.fn, 'r') as f:
                count = sum(1 for _ in f)

            with open(self.fn, 'r') as op:
                with tqdm(desc='parsing datalist {}...'.format(self.fn), leave=self.verbose) as pbar:
                    for l, this_line in enumerate(op):
                        pbar.update()
                        ## parse the datalist entry line
                        if this_line[0] != '#' and this_line[0] != '\n' and this_line[0].rstrip() != '':
                            md = copy.deepcopy(self.metadata)
                            md['name'] = utils.fn_basename2(os.path.basename(self.fn))
                            
                            ## generate the dataset object to yield
                            data_set = DatasetFactory(mod=this_line, weight=self.weight, uncertainty=self.uncertainty, parent=self,
                                                      src_region=self.region, invert_region=self.invert_region, metadata=md, mask=self.mask,
                                                      src_srs=self.src_srs, dst_srs=self.dst_srs, x_inc=self.x_inc, y_inc=self.y_inc,
                                                      stack_fltrs=self.stack_fltrs, pnt_fltrs=self.pnt_fltrs, sample_alg=self.sample_alg,
                                                      cache_dir=self.cache_dir, verbose=self.verbose)._acquire_module()
                            if data_set is not None and data_set.valid_p(
                                    fmts=DatasetFactory._modules[data_set.data_format]['fmts']
                            ):
                                data_set.initialize()
                                ## filter with input source region, if necessary
                                ## check input source region against the dataset region found
                                ## in its inf file.
                                if self.region is not None and self.region.valid_p(check_xy=True):
                                    inf_region = regions.Region().from_list(data_set.infos.minmax)
                                    if inf_region.valid_p():
                                        inf_region.wmin = data_set.weight
                                        inf_region.wmax = data_set.weight
                                        inf_region.umin = data_set.uncertainty
                                        inf_region.umax = data_set.uncertainty

                                        if regions.regions_intersect_p(
                                                inf_region,
                                                self.region if data_set.transformer is None else data_set.trans_region
                                        ):
                                            for ds in data_set.parse(): # fill self.data_entries with each dataset for use outside the yield.
                                                self.data_entries.append(ds) 
                                                yield(ds)

                                    else:
                                        for ds in data_set.parse(): # fill self.data_entries with each dataset for use outside the yield.
                                            self.data_entries.append(ds) 
                                            yield(ds)

                                else:
                                    for ds in data_set.parse(): # fill self.data_entries with each dataset for use outside the yield.
                                        self.data_entries.append(ds)
                                        yield(ds)

        ## self.fn is not a file-name, so check if self.data_entries not empty
        ## and return the dataset objects found there.
        elif len(self.data_entries) > 0:
            for data_set in self.data_entries:
                for ds in data_set.parse():
                    yield(ds)
        else:
            if self.verbose:
                utils.echo_warning_msg(
                    'could not open datalist/entry {}'.format(self.fn)
                )
                
class ZIPlist(ElevationDataset):
    """Zip file parser.

    Parse supported datasets from a zipfile.
    """
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def generate_inf(self, callback=lambda: False):
        """generate and return the infos blob of the datalist.

        datasets parsed from the datalist may have variable srs, so
        we want to transform those to whatever self.dst_srs is, if it
        exists, in order to properly fill the datalist/json srs and regions...
        """

        _region = self.region
        self.region = None
        out_regions = []
        out_srs = []
        self.infos.file_hash = self.infos.generate_hash()
        for entry in self.parse():
            if self.verbose:
                callback()
                
            entry_minmax = entry.infos.minmax
            ## entry has an srs and dst_srs is set, so lets transform the region to suit
            if entry.src_srs is not None:
                out_srs.append(entry.src_srs)
                if self.dst_srs is not None:
                    entry_region = regions.Region().from_list(entry_minmax)
                    entry_region.src_srs = entry.src_srs
                    entry_region.warp(self.dst_srs)
                    entry_minmax = entry_region.export_as_list(include_z=True)

            entry_region = regions.Region().from_list(entry_minmax)
            if entry_region.valid_p():
                out_regions.append(entry_region)
                self.infos.numpts += entry.infos.numpts

        ## merge all the gathered regions
        region_count = 0
        out_region = None
        for this_region in out_regions:
            if this_region.valid_p():
                if region_count == 0:
                    out_region = this_region
                    region_count += 1
                else:
                    out_region = regions.regions_merge(out_region, this_region)
                    
        if out_region is not None:
            self.infos.minmax = out_region.export_as_list(include_z=True)
            self.infos.wkt = out_region.export_as_wkt()

        ## set the epsg for the datalist
        if self.infos.src_srs is None:
            if self.src_srs is not None:
                self.infos.src_srs = self.src_srs
            elif out_srs:
                if all(x == out_srs[0] for x in out_srs):
                    self.infos.src_srs = out_srs[0]

                self.src_srs = self.infos.src_srs
                
        self.region = _region
        return(self.infos)
                
    def parse(self):
        import zipfile
        exts = [DatasetFactory()._modules[x]['fmts'] for x in DatasetFactory()._modules.keys()]
        exts = [x for y in exts for x in y]
        datalist = []
        if self.fn.split('.')[-1].lower() == 'zip':
            with zipfile.ZipFile(self.fn) as z:
                zfs = z.namelist()
                for ext in exts:
                    for zf in zfs:
                        if ext == zf.split('.')[-1]:
                            datalist.append(os.path.basename(zf))
                            
        for this_data in datalist:
            this_line = utils.p_f_unzip(self.fn, [this_data])[0]
            data_set = DatasetFactory(mod=this_line, data_format=None, weight=self.weight, uncertainty=self.uncertainty, parent=self,
                                      src_region=self.region, invert_region=self.invert_region, x_inc=self.x_inc, mask=self.mask,
                                      y_inc=self.y_inc, metadata=copy.deepcopy(self.metadata), src_srs=self.src_srs,
                                      dst_srs=self.dst_srs, cache_dir=self.cache_dir, verbose=self.verbose)._acquire_module()
            if data_set is not None and data_set.valid_p(
                    fmts=DatasetFactory._modules[data_set.data_format]['fmts']
            ):
                data_set.initialize()
                if self.region is not None and self.region.valid_p(check_xy=True):
                    try:
                        inf_region = regions.Region().from_string(
                            data_set.infos.wkt
                        )
                    except:
                        try:
                            inf_region = regions.Region().from_list(
                                data_set.infos.minmax
                            )
                        except:
                            inf_region = self.region.copy()
                            
                    inf_region.wmin = data_set.weight
                    inf_region.wmax = data_set.weight
                    inf_region.umin = data_set.uncertainty
                    inf_region.umax = data_set.uncertainty
                    if inf_region.valid_p(check_xy=True):
                        if regions.regions_intersect_p(inf_region, self.region):
                            for ds in data_set.parse():
                                self.data_entries.append(ds)
                                yield(ds)
                    else:
                        if self.verbose:
                            utils.echo_warning_msg('invalid inf file: {}.inf, skipping'.format(data_set.fn))
                else:
                    for ds in data_set.parse():
                        self.data_entries.append(ds)
                        yield(ds)
            
            utils.remove_glob('{}*'.format(this_data))
        
class Fetcher(ElevationDataset):
    """The default fetches dataset type; dlim Fetcher dataset class

    This is used in waffles/dlim for on-the-fly remote data
    parsing and processing.
    
    See `fetches`, `cudem.fetches` for more information on usage and the
    various fetches modules supported.

    If a fetch module needs special processing define a sub-class
    of Fetcher and redefine the set_ds(self, result) function which yields a
    list of dlim dataset objects, where result is an item from the fetch result list.
    Otherwise, this Fetcher class can be used as default if the fetched data comes
    in a normal sort of way.

    Generally, though not always, if the fetched data is a raster then
    there is no need to redefine set_ds, though if the raster has insufficient
    information, such as with Copernicus, whose nodata value is not
    specified in the geotiff files, it may be best to create a simple
    sub-class for it.
    """

    def __init__(self, keep_fetched_data = True, outdir = None, check_size = True, **kwargs):
        super().__init__(remote=True, **kwargs)        
        self.outdir = outdir if outdir is not None else self.cache_dir # cache directory to store fetched data
        self.fetch_module = fetches.FetchesFactory(
            mod=self.fn, src_region=self.region, verbose=self.verbose, outdir=self.outdir#outdir=self.cache_dir
        )._acquire_module() # the fetches module
        if self.fetch_module is None:
            pass

        self.metadata['name'] = self.fn # set the metadata name from the input fetches module
        self.check_size = check_size # check the size of the fetched data
        self.keep_fetched_data = keep_fetched_data # retain fetched data after processing
        ## breaks when things not set...
        # src_horz, src_vert = gdalfun.epsg_from_input(self.fetch_module.src_srs)
        # self.metadata = {'name':self.fn, 'title':self.fn, 'source':self.fn, 'date':None,
        #                  'data_type':self.data_format, 'resolution':None, 'hdatum':src_horz,
        #                  'vdatum':src_vert, 'url':None}

        
    def generate_inf(self, callback=lambda: False):
        """generate a infos dictionary from the Fetches dataset"""

        tmp_region = self.fetch_module.region if self.region is None else self.region.copy()
        self.infos.minmax = tmp_region.export_as_list()    
        self.infos.wkt = tmp_region.export_as_wkt()
        self.infos.src_srs = self.fetch_module.src_srs
        return(self.infos)

    def parse(self):
        self.fetch_module.run()
        with tqdm(
                total=len(self.fetch_module.results),
                desc='parsing datasets from datalist fetches {} @ {}'.format(self.fetch_module, self.weight),
            leave=self.verbose
        ) as pbar:
            for result in self.fetch_module.results:
                status = self.fetch_module.fetch(result, check_size=self.check_size)
                if status == 0:
                    for this_ds in self.set_ds(result):
                        if this_ds is not None:
                            f_name = os.path.relpath(this_ds.fn, self.fetch_module._outdir)
                            if f_name == '.':
                                f_name = this_ds.fn

                            mod_name = os.path.dirname(utils.fn_basename2(f_name))
                            if mod_name == '':
                                mod_name = self.fetch_module.name
                                
                            this_ds.metadata['name'] = mod_name
                            this_ds.remote = True
                            this_ds.initialize()
                            for ds in this_ds.parse():
                                #ds.metadata['name'] = os.path.basename(ds.fn).split('.')[0]
                                yield(ds)
                                
                            #if not self.keep_fetched_data:
                            #    utils.remove_glob('{}*'.format(this_ds.fn))
                pbar.update()
                
        if not self.keep_fetched_data:
            utils.remove_glob('{}*'.format(self.fn))
            
    def set_ds(self, result):
        ## try to get the SRS info from the result if it's a gdal file
        ## fix this.
        try:
            vdatum = self.fetch_module.vdatum
            src_srs = gdalfun.gdal_get_srs(os.path.join(self.fetch_module._outdir, result[1]))
            horz_epsg, vert_epsg = gdalfun.epsg_from_input(src_srs)
            if vert_epsg is None:
                vert_epsg = vdatum

            #utils.echo_msg('srs: {}+{}'.format(horz_epsg, vert_epsg))
            if vert_epsg is not None:
                self.fetch_module.src_srs = '{}+{}'.format(horz_epsg, vert_epsg)
            else:
                self.fetch_module.src_srs = src_srs
        except:
            pass
        
        ds = DatasetFactory(mod=os.path.join(self.fetch_module._outdir, result[1]), data_format=self.fetch_module.data_format, weight=self.weight,
                            parent=self, src_region=self.region, invert_region=self.invert_region, metadata=copy.deepcopy(self.metadata),
                            mask=self.mask, uncertainty=self.uncertainty, x_inc=self.x_inc, y_inc=self.y_inc, stack_fltrs=self.stack_fltrs,
                            pnt_fltrs=self.pnt_fltrs, src_srs=self.fetch_module.src_srs, dst_srs=self.dst_srs, verbose=self.verbose,
                            cache_dir=self.fetch_module._outdir, remote=True)._acquire_module()
        yield(ds)

class NEDFetcher(Fetcher):
    """National Elevation Dataset from USGS

    This is a wrapper shortcut for fetching NED DEMs from USGS' The National Map (TNM)
    """
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def set_ds(self, result):
        
        src_dem = os.path.join(self.fetch_module._outdir, result[1])
        # grits_filter = grits.GritsFactory(mod='flats', src_dem=src_dem, size_threshold=10, cache_dir=self.fetch_module._outdir)._acquire_module()
        # grits_filter()
        # ned_fn = grits_filter.dst_dem
        
        ds = DatasetFactory(mod=src_dem, data_format=self.fetch_module.data_format, weight=self.weight,
                            parent=self, src_region=self.region, invert_region=self.invert_region, metadata=copy.deepcopy(self.metadata),
                            mask=self.mask, uncertainty=self.uncertainty, x_inc=self.x_inc, y_inc=self.y_inc, 
                            src_srs=self.fetch_module.src_srs, dst_srs=self.dst_srs, verbose=self.verbose, cache_dir=self.fetch_module._outdir,
                            remote=True, remove_flat=True)._acquire_module()
        yield(ds)

        
class DAVFetcher_CoNED(Fetcher):
    """CoNED from the digital coast 

    This is a wrapper shortcut for fetching CoNED DEMs from the Digital Coast,
    mainly so we can pull the vertical datum info from the DAV metadata since the
    CoNED doesn't assign one to their DEMs.
    """

    __doc__ = '''{}
    
    -----------
    Fetches Module: <CoNED> - {}'''.format(__doc__, fetches.CoNED.__doc__)
    
    
    def __init__(self, keep_fetched_data = True, cog = True, **kwargs):
        super().__init__(**kwargs)
        self.keep_fetched_data = keep_fetched_data
        self.cog = cog

    def parse(self):
        self.fetch_module.run()
        for result in self.fetch_module.results:
            if not self.cog:
                status = self.fetch_module.fetch(result, check_size=self.check_size)
                if status != 0:
                    break
                
            for this_ds in self.set_ds(result):
                if this_ds is not None:
                    this_ds.metadata['name'] = utils.fn_basename2(this_ds.fn)
                    #this_ds.remote = True
                    this_ds.initialize()
                    yield(this_ds)
            
    def set_ds(self, result):
        ## try to get the SRS info from the result
        try:
            vdatum = self.fetch_module.vdatum
            if not self.cog:
                src_srs = gdalfun.gdal_get_srs(os.path.join(self.fetch_module._outdir, result[1]))
            else:
                src_srs = gdalfun.gdal_get_srs(result[0])
                
            horz_epsg, vert_epsg = gdalfun.epsg_from_input(src_srs)
            if vert_epsg is None:
                vert_epsg = vdatum

            utils.echo_msg('srs: {}+{}'.format(horz_epsg, vert_epsg))
            if vert_epsg is not None:
                self.fetch_module.src_srs = '{}+{}'.format(horz_epsg, vert_epsg)
            else:
                self.fetch_module.src_srs = src_srs
        except:
            pass

        if not self.cog:
            ds = DatasetFactory(mod=os.path.join(self.fetch_module._outdir, result[1]), data_format=self.fetch_module.data_format, weight=self.weight,
                                parent=self, src_region=self.region, invert_region=self.invert_region, metadata=copy.deepcopy(self.metadata),
                                mask=self.mask, uncertainty=self.uncertainty, x_inc=self.x_inc, y_inc=self.y_inc,
                                src_srs=self.fetch_module.src_srs, dst_srs=self.dst_srs, verbose=self.verbose, cache_dir=self.fetch_module._outdir,
                                remote=True)._acquire_module()
        else:
            ds = DatasetFactory(mod=result[0], data_format='200', weight=self.weight,
                                parent=self, src_region=self.region, invert_region=self.invert_region, metadata=copy.deepcopy(self.metadata),
                                mask=self.mask, uncertainty=self.uncertainty, x_inc=self.x_inc, y_inc=self.y_inc,
                                src_srs=self.fetch_module.src_srs, dst_srs=self.dst_srs, verbose=self.verbose, cache_dir=self.fetch_module._outdir,
                                check_path=False)._acquire_module()
        yield(ds)

class DAVFetcher_SLR(Fetcher):
    """SLR DEM from the digital coast 

    This is a wrapper shortcut for fetching SLR DEMs from the Digital Coast,
    mainly so we can pull the remove the flattened ring around the actual data.
    """

    __doc__ = '''{}
    
    -----------
    Fetches Module: <SLR> - {}'''.format(__doc__, fetches.SLR.__doc__)
    
    def __init__(self, keep_fetched_data = True, **kwargs):
        super().__init__(**kwargs)

    def set_ds(self, result):
        ## this doesn't work in all cases, update to find and remove flattened areas
        #gdalfun.gdal_set_ndv(os.path.join(self.fetch_module._outdir, result[1]), ndv=-99.0000, convert_array=True)
        
        ds = DatasetFactory(mod=os.path.join(self.fetch_module._outdir, result[1]), data_format=self.fetch_module.data_format, weight=self.weight,
                            parent=self, src_region=self.region, invert_region=self.invert_region, metadata=copy.deepcopy(self.metadata),
                            mask=self.mask, uncertainty=self.uncertainty, x_inc=self.x_inc, y_inc=self.y_inc,
                            src_srs=self.fetch_module.src_srs, dst_srs=self.dst_srs, verbose=self.verbose, cache_dir=self.fetch_module._outdir,
                            remote=True, remove_flat=True)._acquire_module()
        yield(ds)
        
class SWOTFetcher(Fetcher):
    """SWOT L2_HR_Raster data from NASA

set `data_set` to one of (for L2_HR_Raster):
    1 - longitude (64-bit floating-point)
    2 - latitude (64-bit floating-point)
    3 - wse (32-bit floating-point)
    4 - status_flag (8-bit unsigned integer)
    5 - status_flag (32-bit unsigned integer)
    6 - wse_uncert (32-bit floating-point)
    7 - water_area (32-bit floating-point)
    8 - status_flag (8-bit unsigned integer)
    9 - status_flag (32-bit unsigned integer)
    10 - water_area_uncert (32-bit floating-point)
    11 - water_frac (32-bit floating-point)
    12 - water_frac_uncert (32-bit floating-point)
    13 - sig0 (32-bit floating-point)
    14 - status_flag (8-bit unsigned integer)
    15 - status_flag (32-bit unsigned integer)
    16 - sig0_uncert (32-bit floating-point)
    17 - inc (32-bit floating-point)
    18 - cross_track (32-bit floating-point)
    19 - time (64-bit floating-point)
    20 - time (64-bit floating-point)
    21 - n_wse_pix (32-bit unsigned integer)
    22 - n_water_area_pix (32-bit unsigned integer)
    23 - n_sig0_pix (32-bit unsigned integer)
    24 - n_other_pix (32-bit unsigned integer)
    25 - dark_frac (32-bit floating-point)
    26 - status_flag (8-bit unsigned integer)
    27 - status_flag (8-bit unsigned integer)
    28 - layover_impact (32-bit floating-point)
    29 - sig0_cor_atmos_model (32-bit floating-point)
    30 - height_cor_xover (32-bit floating-point)
    31 - geoid_height_above_reference_ellipsoid (32-bit floating-point)
    32 - solid_earth_tide (32-bit floating-point)
    33 - load_tide_fes (32-bit floating-point)
    34 - load_tide_got (32-bit floating-point)
    35 - pole_tide (32-bit floating-point)
    36 - model_dry_tropo_cor (32-bit floating-point)
    37 - model_wet_tropo_cor (32-bit floating-point)
    38 - iono_cor_gim_ka (32-bit floating-point)

    currently supported SWOT products and product options:
    - L2_HR_Raster
      apply_geoid
      data_set

    - L2_HR_PIXC
      apply_geoid
      
    """

    __doc__ = '''{}
    
    -----------
    Fetches Module: <swot> - {}'''.format(__doc__, fetches.SWOT.__doc__)
    

    def __init__(self, data_set = 'wse', apply_geoid = True, classes = None, classes_qual = None,
                 anc_classes = None, remove_class_flags = False, **kwargs):
        super().__init__(**kwargs)
        self.data_set = data_set
        self.apply_geoid = apply_geoid
        self.classes = classes
        self.anc_classes = anc_classes
        self.classes_qual = classes_qual
        self.remove_class_flags = remove_class_flags

    # def fetch_pixc_vec(self, swot_fn):
    #     pixc_vec_filter = utils.fn_basename2(swot_fn).split('PIXC_')[1]
    #     this_pixc_vec = fetches.SWOT(
    #         src_region=None, verbose=self.verbose, outdir=self.fetch_module._outdir, product='L2_HR_PIXCVec', filename_filter=pixc_vec_filter
    #     )
    #     this_pixc_vec.run()
    #     if len(this_pixc_vec.results) == 0:
    #         utils.echo_warning_msg('could not locate associated PIXCVec file for {}'.format(pixc_vec_filter))
    #         return(None)
    #     else:
    #         if this_pixc_vec.fetch(this_pixc_vec.results[0], check_size=self.check_size) == 0:
    #             return(os.path.join(this_pixc_vec._outdir, this_pixc_vec.results[0][1]))
        
    def set_ds(self, result):
        #utils.echo_msg(result)
        swot_fn = os.path.join(self.fetch_module._outdir, result[1])
        if 'L2_HR_PIXC_' in result[-1]:
            #pixc_vec_result = self.fetch_pixc_vec(swot_fn)
            #swot_pixc_vec_fn = pixc_vec_result
            
            sub_ds = SWOT_PIXC(fn=swot_fn, data_format=202, apply_geoid=self.apply_geoid, classes=self.classes, anc_classes=self.anc_classes, classes_qual=self.classes_qual, 
                               remove_class_flags=self.remove_class_flags, src_srs='epsg:4326+3855', dst_srs=self.dst_srs, weight=self.weight, uncertainty=self.uncertainty,
                               src_region=self.region, x_inc=self.x_inc, y_inc=self.y_inc, stack_fltrs=self.stack_fltrs, pnt_fltrs=self.pnt_fltrs, verbose=True,
                               metadata=copy.deepcopy(self.metadata))
            
        elif 'L2_HR_Raster' in result[-1]:
            sub_ds = SWOT_HR_Raster(fn=swot_fn, data_format=203, apply_geoid=self.apply_geoid, src_srs=None, dst_srs=self.dst_srs,
                                    weight=self.weight, uncertainty=self.uncertainty, src_region=self.region,
                                    x_inc=self.x_inc, y_inc=self.y_inc, stack_fltrs=self.stack_fltrs, pnt_fltrs=self.pnt_fltrs,
                                    verbose=True, metadata=copy.deepcopy(self.metadata))
        else:
            utils.echo_warning_msg('{} is not a currently supported dlim dataset'.format(result[-1]))
            sub_ds = None
            
        yield(sub_ds)
            
class IceSatFetcher(Fetcher):
    """IceSat data from NASA

    See `fetches --modules icesat` for fetching parameters

    -----------
    Parameters:
    
    want_topo (bool): extract topography
    want_bathy (bool): extract bathymetry
    bathy_thresh (int): bathymetry extraction threshhold.
    topo_thresh (int): topography extraction threshhold.
    water_surface (str): 'mean_tide', 'geoid' or 'ellipsoid' water surface
    """

    __doc__ = '''{}
    
    -----------
    Fetches Module: <icesat> - {}'''.format(__doc__, fetches.IceSat.__doc__)
    
    def __init__(self, want_bathy=True, want_topo=False, water_surface='mean_tide', bathy_thresh = 50, topo_thresh = 30, **kwargs):
        super().__init__(**kwargs)
        self.want_topo = want_topo
        self.want_bathy = want_bathy
        self.water_surface = water_surface
        if self.water_surface not in ['mean_tide', 'geoid', 'ellipsoid']:
            self.water_surface = 'mean_tide'

        self.bathy_thresh = utils.float_or(bathy_thresh, 30)
        self.topo_thresh = utils.float_or(topo_thresh, 30)
        self.atl_fn = None

    def parse(self):
        self.fetch_module.run()
        yield(self)

    def close_atl_h5(self):
        self.atl_03_fn = None
        self.atl_08_fn = None
        self.atl_03_f.close()
        if self.atl_08_f is not None:
            self.atl_08_f.close()
    
    def init_atl_03_h5(self):
        """initialize the atl03 and atl08 h5 files"""
        
        self.atl_03_f = h5.File(self.atl_03_fn, 'r')
        if 'short_name' not in self.atl_03_f.attrs.keys():
            utils.echo_error_msg('this file does not appear to be an ATL file')
            self.atl_03_f.close()
            return(None)
        
    def init_atl_h5(self):
        """initialize the atl03 and atl08 h5 files"""

        self.atl_03_f = h5.File(self.atl_03_fn, 'r')
        if 'short_name' not in self.atl_03_f.attrs.keys():
            utils.echo_error_msg('this file does not appear to be an ATL file')
            self.atl_03_f.close()
        
        if self.atl_08_fn is not None:
            self.atl_08_f = h5.File(self.atl_08_fn, 'r')
        
            if 'short_name' not in self.atl_08_f.attrs.keys():
                utils.echo_error_msg('this file does not appear to be an ATL file')
                self.atl_08_f.close()
        else:
            self.atl_08_f = None

    def fetch_and_yield_results(self, fetch_data = True):
        for result in self.fetch_module.results:
            if fetch_data:
                if self.fetch_module.fetch(result, check_size=self.check_size) == 0:
                    yield(result)
                else:
                    self.fetch_module.results.append(result)
            else:
                yield(result)
            
    def yield_ds(self):
        with tqdm(
                total=len(self.fetch_module.results),
                desc='parsing datasets from datalist fetches {} @ {}'.format(self.fetch_module, self.weight),
                leave=self.verbose
        ) as pbar:
            for result in self.fetch_and_yield_results(fetch_data=True):
                ## load the atl 03 data
                self.atl_03_fn = os.path.join(self.fetch_module._outdir, result[1])
                
                ## load the atl 08 data
                atl_08_result = self.fetch_atl_08()
                self.atl_08_fn = atl_08_result

                try:
                    self.init_atl_h5()
                except Exception as e:
                    utils.echo_error_msg('could not initialize data {}'.format(e))
                    pbar.update()
                    self.close_atl_h5()
                    continue
                                                                                   
                dataset = None

                ## Fetch the associated ATL08 dataset for bare-earth processing
                for i in range(1, 4):
                    ## read the data from the atl03
                    dataset = self.read_atl_data('{}'.format(i))
                    ## Bare-Earth topography (class 1)
                    if self.want_topo:
                        bare_earth_dataset = self.extract_topography(dataset, thresh=self.bathy_thresh)
                        if bare_earth_dataset is not None:
                            bare_earth_dataset = np.column_stack((bare_earth_dataset[1], bare_earth_dataset[0], bare_earth_dataset[2]))
                            bare_earth_points = np.rec.fromrecords(bare_earth_dataset, names='x, y, z')
                            yield(bare_earth_points)

                    ## Bathymetry via C-Shelph
                    if self.want_bathy:
                        bathy_dataset = self.extract_bathymetry(dataset, thresh=self.bathy_thresh)
                        if bathy_dataset is not None:
                            bathy_dataset = np.column_stack((bathy_dataset[1], bathy_dataset[0], bathy_dataset[2]))
                            bathy_points = np.rec.fromrecords(bathy_dataset, names='x, y, z')
                            yield(bathy_points)

                pbar.update()
                self.close_atl_h5()

    def read_atl_data(self, laser_num):
        """Read data from an ATL03 file

        Adapted from 'cshelph' https://github.com/nmt28/C-SHELPh.git 
        and 'iVert' https://github.com/ciresdem/ivert.git

        laser_num is 1, 2 or 3
        surface is 'mean_tide', 'geoid' or 'ellipsoid'
        """

        orientation = self.atl_03_f['/orbit_info/sc_orient'][0]
        # selects the strong beams only [we can include weak beams later on]
        orientDict = {0:'l', 1:'r', 21:'error'}
        laser = 'gt' + laser_num + orientDict[orientation]

        ## Read in the required photon level data
        photon_h = self.atl_03_f['/' + laser + '/heights/h_ph'][...,]
        latitude = self.atl_03_f['/' + laser + '/heights/lat_ph'][...,]
        longitude = self.atl_03_f['/' + laser + '/heights/lon_ph'][...,]
        conf = self.atl_03_f['/' + laser + '/heights/signal_conf_ph/'][...,0]
        dist_ph_along = self.atl_03_f['/' + laser + '/heights/dist_ph_along'][...,]
        this_N = latitude.shape[0]

        ## Read in the geolocation level data
        segment_ph_cnt = self.atl_03_f['/' + laser + '/geolocation/segment_ph_cnt'][...,]
        segment_id = self.atl_03_f['/' + laser + '/geolocation/segment_id'][...,]
        segment_dist_x = self.atl_03_f['/' + laser + '/geolocation/segment_dist_x'][...,]
        seg_ph_count = self.atl_03_f['/' + laser + '/geolocation/segment_ph_cnt'][...,]
        ref_elev = self.atl_03_f['/' + laser + '/geolocation/ref_elev'][...,]
        ref_azimuth = self.atl_03_f['/' + laser + '/geolocation/ref_azimuth'][...,]
        ph_index_beg = self.atl_03_f['/' + laser + '/geolocation/ph_index_beg'][...,]
        altitude_sc = self.atl_03_f['/' + laser + '/geolocation/altitude_sc'][...,]

        ## Read in the geoid data
        photon_geoid = self.atl_03_f['/' + laser + '/geophys_corr/geoid'][...,]
        photon_geoid_f2m = self.atl_03_f['/' + laser + '/geophys_corr/geoid_free2mean'][...,]

        ## Create a dictionary with (segment_id --> index into ATL03 photons)
        ## lookup pairs, for the starting photon of each segment
        segment_indices = np.concatenate(([0], np.cumsum(segment_ph_cnt)[:-1]))
        segment_index_dict = dict(zip(segment_id, segment_indices))

        ## Compute the total along-track distances.
        segment_dist_dict = dict(zip(segment_id, segment_dist_x))

        ## Determine where in the array each segment index needs to look.
        ph_segment_ids = segment_id[np.searchsorted(segment_indices, np.arange(0.5, this_N, 1))-1]
        ph_segment_dist_x = np.array(list(map((lambda pid: segment_dist_dict[pid]), ph_segment_ids)))
        dist_x = ph_segment_dist_x + dist_ph_along

        ## meantide/geoid heights
        h_geoid_dict = dict(zip(segment_id, photon_geoid))
        ph_h_geoid = np.array(list(map((lambda pid: h_geoid_dict[pid]), ph_segment_ids)))

        h_meantide_dict = dict(zip(segment_id, photon_geoid_f2m))
        ph_h_meantide = np.array(list(map((lambda pid: h_meantide_dict[pid]), ph_segment_ids)))

        photon_h_geoid = photon_h - ph_h_geoid
        photon_h_meantide = photon_h - (ph_h_geoid + ph_h_meantide)

        ph_h_classed = np.zeros(photon_h.shape)
        ## Read in the atl08 data
        if self.atl_08_f is not None:
            atl_08_classed_pc_flag  = self.atl_08_f['/' + laser + '/signal_photons/classed_pc_flag'][...,]
            atl_08_ph_segment_id = self.atl_08_f['/' + laser + '/signal_photons/ph_segment_id'][...,]
            atl_08_classed_pc_indx = self.atl_08_f['/' + laser + '/signal_photons/classed_pc_indx'][...,]

            # atl_08_seg_beg = self.atl_08_f['/' + laser + '/land_segments/segment_id_beg'][...,]
            # atl_08_seg_end = self.atl_08_f['/' + laser + '/land_segments/segment_id_end'][...,]
            # atl_08_watermask = self.atl_08_f['/' + laser + '/land_segments/segment_watermask'][...,]
            
            # Type codes:
            # -1 : uncoded
            #  0 : noise
            #  1 : ground
            #  2 : canopy
            #  3 : top of canopy

            dict_success = False
            while not dict_success:
                try:
                    atl_08_ph_segment_indx = np.array(list(map((lambda pid: segment_index_dict[pid]), atl_08_ph_segment_id)))
                except KeyError as e:
                    # One of the atl08_ph_segment_id entries does not exist in the atl03 granule, which
                    # causes problems here. Eliminate it from the list and try again.
                    problematic_id = e.args[0]
                    good_atl08_mask = (atl_08_ph_segment_id != problematic_id)
                    atl_08_classed_pc_flag = atl_08_classed_pc_flag[good_atl08_mask]
                    atl_08_ph_segment_id = atl_08_ph_segment_id[good_atl08_mask]
                    atl_08_classed_pc_indx = atl_08_classed_pc_indx[good_atl08_mask]
                    # Then, try the loop again.
                    continue
                dict_success = True

            atl_08_ph_index = np.array(atl_08_ph_segment_indx + atl_08_classed_pc_indx - 1 , dtype=int)
            ph_h_classed[atl_08_ph_index] = atl_08_classed_pc_flag
            
        return(latitude, longitude, photon_h_meantide if self.water_surface=='mean_tide' else photon_h_geoid if self.water_surface=='geoid' else photon_h,
               conf, ref_elev, ref_azimuth, ph_index_beg, segment_id, altitude_sc, seg_ph_count, ph_h_classed)

    def fetch_atl_08(self):
        atl_08_filter = utils.fn_basename2(self.atl_03_fn).split('ATL03_')[1]
        this_atl08 = fetches.IceSat(
            src_region=None, verbose=self.verbose, outdir=self.fetch_module._outdir, short_name='ATL08', filename_filter=atl_08_filter
        )
        this_atl08.run()
        if len(this_atl08.results) == 0:
            utils.echo_warning_msg('could not locate associated atl08 file for {}'.format(atl_08_filter))
            return(None)
        else:
            if this_atl08.fetch(this_atl08.results[0], check_size=self.check_size) == 0:
                return(os.path.join(this_atl08._outdir, this_atl08.results[0][1]))

    def extract_topography(self, dataset, thresh=None, min_buffer = 0, surface_buffer = 1.1, lat_res = 10 , h_res = .5):
        latitude, longitude, photon_h, conf, ref_elev, ref_azimuth, ph_index_beg, segment_id, alt_sc, seg_ph_count, ph_h_classed = dataset

        epsg_code = cshelph.convert_wgs_to_utm(latitude[0], longitude[0])
        epsg_num = int(epsg_code.split(':')[-1])
        utm_proj = pyproj.Proj(epsg_code)
        lon_utm, lat_utm = utm_proj(longitude, latitude)
        ph_num_per_seg = seg_ph_count[ph_index_beg>0]
        ph_num_per_seg = ph_num_per_seg.astype(np.int64)
        ph_ref_elev = cshelph.ref_linear_interp(ph_num_per_seg, ref_elev[ph_index_beg>0])
        ph_ref_azimuth = cshelph.ref_linear_interp(ph_num_per_seg, ref_azimuth[ph_index_beg>0])
        ph_sat_alt = cshelph.ref_linear_interp(ph_num_per_seg, alt_sc[ph_index_beg>0])

        ## Aggregate data into dataframe
        dataset_sea = pd.DataFrame(
            {'latitude': lat_utm,
             'longitude': lon_utm,
             'photon_height': photon_h,
             'confidence':conf,
             'ref_elevation':ph_ref_elev,
             'ref_azminuth':ph_ref_azimuth,
             'ref_sat_alt':ph_sat_alt,
             'ph_h_classed': ph_h_classed},
            columns=['latitude', 'longitude', 'photon_height', 'confidence', 'ref_elevation', 'ref_azminuth', 'ref_sat_alt', 'ph_h_classed']
        )
        dataset_sea1 = dataset_sea[dataset_sea.confidence == 4]
        dataset_sea1 = dataset_sea1[(dataset_sea1['ph_h_classed'] == 1)]
        dataset_sea1 = dataset_sea1[(dataset_sea1['photon_height'] > min_buffer)]
        if self.region is not None:
            xyz_region = self.region.copy()
            xyz_region.epsg = 'epsg:4326'
            xyz_region.warp('epsg:{}'.format(epsg_num))
            dataset_sea1 = dataset_sea1[(dataset_sea1['latitude'] > xyz_region.ymin) & (dataset_sea1['latitude'] < xyz_region.ymax)]
            dataset_sea1 = dataset_sea1[(dataset_sea1['longitude'] > xyz_region.xmin) & (dataset_sea1['longitude'] < xyz_region.xmax)]

        if len(dataset_sea1) != 0:
            binned_data_sea = cshelph.bin_data(dataset_sea1, lat_res, h_res)
            if binned_data_sea is not None:
                #sea_height = cshelph.get_sea_height(binned_data_sea, surface_buffer)
                lats, lons, sea_height = cshelph.get_bin_height(binned_data_sea, surface_buffer)
                topo_ds = np.column_stack((lons, lats, sea_height))
                topo_ds = np.rec.fromrecords(topo_ds, names='x, y, z')
                topo_ds = topo_ds[~np.isnan(topo_ds['z'])]
                med_water_surface_h = np.nanmedian(topo_ds['z'])
                topo_ds = topo_ds[topo_ds['z'] > med_water_surface_h + (h_res * 2)]

                transformer = pyproj.Transformer.from_crs("EPSG:"+str(epsg_num), "EPSG:4326", always_xy=True)
                lon_wgs84, lat_wgs84 = transformer.transform(topo_ds['x'], topo_ds['y'])

                return(lat_wgs84, lon_wgs84, topo_ds['z'])

        return(None)
            
    ## C-Shelph bathymetric processing
    def get_water_temp(self):
        #import xarray as xr
        #from eosdis_store import EosdisStore

        if self.region is not None:
            this_region = self.region.copy()
        else:
            this_region = regions.Region().from_list(self.infos.minmax)
            
        time_start = self.f.attrs['time_coverage_start'].decode('utf-8')
        time_end = self.f.attrs['time_coverage_end'].decode('utf-8')
        this_sst = fetches.MUR_SST(src_region = this_region, verbose=self.verbose, outdir=self.cache_dir, time_start=time_start, time_end=time_end)
        this_sst.run()

        return(20)
        ## todo do this wityout eosdis-store!
        # datasets = [xr.open_zarr(EosdisStore(x[0]), consolidated=False) for x in this_sst.results]
        # if len(datasets) > 1:
        #     ds = xr.concat(datasets, 'time')
        # else:
        #     ds = datasets[0]

        # lats = slice(this_region.ymin, this_region.ymax)
        # lons = slice(this_region.xmin, this_region.xmax)
        # sea_temp = ds.analysed_sst.sel(lat=lats, lon=lons)
        # sst = round(np.nanmedian(sea_temp.values)-273,2)
        
        # return(sst)
    
    def extract_bathymetry(
            self, dataset, thresh = None, min_buffer = -40, max_buffer = 5,
            start_lat = False, end_lat = False, lat_res = 10 , h_res = .5,
            surface_buffer = -.5, water_temp = None
    ):
        """Extract bathymetry from an ATL03 file. 

        This uses C-Shelph to locate, extract and process bathymetric photons.
        This function is adapted from the C-Shelph CLI
        """
        
        water_temp = utils.float_or(water_temp)
        latitude, longitude, photon_h, conf, ref_elev, ref_azimuth, ph_index_beg, segment_id, alt_sc, seg_ph_count, ph_h_classed = dataset
        epsg_code = cshelph.convert_wgs_to_utm(latitude[0], longitude[0])
        epsg_num = int(epsg_code.split(':')[-1])
        utm_proj = pyproj.Proj(epsg_code)
        lon_utm, lat_utm = utm_proj(longitude, latitude)
        ph_num_per_seg = seg_ph_count[ph_index_beg>0]
        ph_num_per_seg = ph_num_per_seg.astype(np.int64)
        ph_ref_elev = cshelph.ref_linear_interp(ph_num_per_seg, ref_elev[ph_index_beg>0])
        ph_ref_azimuth = cshelph.ref_linear_interp(ph_num_per_seg, ref_azimuth[ph_index_beg>0])
        ph_sat_alt = cshelph.ref_linear_interp(ph_num_per_seg, alt_sc[ph_index_beg>0])

        ## Aggregate data into dataframe
        dataset_sea = pd.DataFrame(
            {'latitude': lat_utm,
             'longitude': lon_utm,
             'photon_height': photon_h,
             'confidence':conf,
             'ref_elevation':ph_ref_elev,
             'ref_azminuth':ph_ref_azimuth,
             'ref_sat_alt':ph_sat_alt,
             'ph_h_classed': ph_h_classed},
            columns=['latitude', 'longitude', 'photon_height', 'confidence', 'ref_elevation', 'ref_azminuth', 'ref_sat_alt', 'ph_h_classed']
        )
        dataset_sea1 = dataset_sea[(dataset_sea.confidence != 0)  & (dataset_sea.confidence != 1)]
        dataset_sea1 = dataset_sea1[(dataset_sea1['photon_height'] > min_buffer) & (dataset_sea1['photon_height'] < max_buffer)]
        if self.region is not None:
            xyz_region = self.region.copy()
            xyz_region.epsg = 'epsg:4326'
            xyz_region.warp('epsg:{}'.format(epsg_num))
            dataset_sea1 = dataset_sea1[(dataset_sea1['latitude'] > xyz_region.ymin) & (dataset_sea1['latitude'] < xyz_region.ymax)]
            dataset_sea1 = dataset_sea1[(dataset_sea1['longitude'] > xyz_region.xmin) & (dataset_sea1['longitude'] < xyz_region.xmax)]

        if len(dataset_sea1) > 0:
            binned_data_sea = cshelph.bin_data(dataset_sea1, lat_res, h_res)
            if binned_data_sea is not None:
                if water_temp is None:
                    try:
                        #water_temp = cshelph.get_water_temp(self.fn, latitude, longitude)
                        ## water_temp via fetches instead of earthaccess
                        water_temp = self.get_water_temp()
                    except Exception as e:
                        #utils.echo_warning_msg('NO SST PROVIDED OR RETRIEVED: 20 degrees C assigned')
                        water_temp = 20

                sea_height = cshelph.get_sea_height(binned_data_sea, surface_buffer)
                if sea_height is not None:
                    med_water_surface_h = np.nanmedian(sea_height)

                    ## Correct for refraction
                    ref_x, ref_y, ref_z, ref_conf, raw_x, raw_y, raw_z, ph_ref_azi, ph_ref_elev = cshelph.refraction_correction(
                        water_temp, med_water_surface_h, 532, dataset_sea1.ref_elevation, dataset_sea1.ref_azminuth, dataset_sea1.photon_height,
                        dataset_sea1.longitude, dataset_sea1.latitude, dataset_sea1.confidence, dataset_sea1.ref_sat_alt
                    )
                    depth = med_water_surface_h - ref_z

                    # Create new dataframe with refraction corrected data
                    dataset_bath = pd.DataFrame({'latitude': raw_y, 'longitude': raw_x, 'cor_latitude':ref_y, 'cor_longitude':ref_x, 'cor_photon_height':ref_z,
                                                 'photon_height': raw_z, 'confidence':ref_conf, 'depth':depth},
                                                columns=['latitude', 'longitude', 'photon_height', 'cor_latitude','cor_longitude', 'cor_photon_height',
                                                         'confidence', 'depth'])

                    # Bin dataset again for bathymetry
                    binned_data = cshelph.bin_data(dataset_bath, lat_res, h_res)

                    if binned_data is not None:
                        # Find bathymetry
                        bath_height, geo_df = cshelph.get_bath_height(binned_data, thresh, med_water_surface_h, h_res)
                        if bath_height is not None:

                            transformer = pyproj.Transformer.from_crs("EPSG:"+str(epsg_num), "EPSG:4326", always_xy=True)
                            lon_wgs84, lat_wgs84 = transformer.transform(geo_df.longitude.values, geo_df.latitude.values)

                            bath_height = [x for x in bath_height if ~np.isnan(x)]
                            return(lat_wgs84, lon_wgs84, geo_df.depth.values*-1)
            
        return(None)
                            
# class DAVFetcher_CUDEM(Fetcher):
#     """CUDEM from the digital coast 

#     This is a wrapper shortcut for fetching CUDEMs from the Digital Coast,
#     mainly so we can pull the remove the flattened ring around the actual data.
#     """
    
#     def __init__(self, keep_fetched_data = True, **kwargs):
#         super().__init__(**kwargs)

#     def set_ds(self, result):

#         # with gdalfun.gdal_datasource(os.path.join(self.fetch_module._outdir, result[1]), update=True) as src_ds:
#         #     if src_ds is not None:
#         #         ds_config = gdalfun.gdal_infos(src_ds)
#         #         curr_nodata = ds_config['ndv']

#         ## this doesn't work in all cases, update to find and remove flattened areas
#         gdalfun.gdal_set_ndv(os.path.join(self.fetch_module._outdir, result[1]), ndv=-99.0000, convert_array=True)
        
#         ds = DatasetFactory(mod=os.path.join(self.fetch_module._outdir, result[1]), data_format=self.fetch_module.data_format, weight=self.weight,
#                             parent=self, src_region=self.region, invert_region=self.invert_region, metadata=copy.deepcopy(self.metadata),
#                             mask=self.mask, uncertainty=self.uncertainty, x_inc=self.x_inc, y_inc=self.y_inc, src_srs=self.fetch_module.src_srs,
#                             dst_srs=self.dst_srs, verbose=self.verbose, cache_dir=self.fetch_module._outdir, remote=True)._acquire_module()
#         yield(ds)        
                
class GMRTFetcher(Fetcher):
    """GMRT Gridded data.

    -----------
    Parameters:
    
    swath_only: onlt return MB swath data
    """
    
    __doc__ = '''{}
    
    -----------
    Fetches Module: <gmrt> - {}'''.format(__doc__, fetches.GMRT.__doc__)
    
    def __init__(self, swath_only = False, **kwargs):
        super().__init__(**kwargs)
        self.swath_only = swath_only

    def set_ds(self, result):
        swath_mask=None
        gmrt_fn = os.path.join(self.fetch_module._outdir, result[1])
        with gdalfun.gdal_datasource(gmrt_fn, update = 1) as src_ds:
            md = src_ds.GetMetadata()
            md['AREA_OR_POINT'] = 'Point'
            src_ds.SetMetadata(md)
            gdalfun.gdal_set_srs(src_ds)
            gdalfun.gdal_set_ndv(src_ds)

        if self.swath_only:
            if fetches.Fetch(
                    self.fetch_module._gmrt_swath_poly_url, verbose=self.verbose
            ).fetch_file(os.path.join(self.fetch_module._outdir, 'gmrt_swath_polygons.zip')) == 0:
                swath_shps = utils.p_unzip(
                    os.path.join(self.fetch_module._outdir, 'gmrt_swath_polygons.zip'),
                    exts=['shp', 'shx', 'prj', 'dbf'],
                    outdir=self.cache_dir
                )
                #swath_shp = None
                #swath_mask = os.path.join(self.fetch_module._outdir, 'gmrt_swath_polygons.tif')
                
                for v in swath_shps:
                    if '.shp' in v:
                        #swath_shp = v
                        swath_mask = {'mask': v, 'invert_mask': True}
                        break
                    
                if not os.path.exists(swath_mask['mask']):
                    utils.echo_error_msg('could not find gmrt swath polygons...')
                    self.swath_only = False
                    swath_mask = None
                # else:
                #     import shutil
                #     tmp_gmrt = '{}_clip.tif'.format(utils.fn_basename2(gmrt_fn))
                #     shutil.copyfile(gmrt_fn, tmp_gmrt)
                    
                #     gi = gdalfun.gdal_infos(gmrt_fn)
                #     gr_cmd = 'gdal_rasterize -burn {} -l {} {} {}'\
                #         .format(gi['ndv'], os.path.basename(swath_shp).split('.')[0], swath_shp, tmp_gmrt)
                #     out, status = utils.run_cmd(gr_cmd, verbose=self.verbose)

                #     with gdalfun.gdal_datasource(tmp_gmrt, update=True) as tmp_ds:
                #         b = tmp_ds.GetRasterBand(1)
                #         a = b.ReadAsArray()
                #         a[a != gi['ndv']] = 0
                #         a[a == gi['ndv']] = 1
                #         b.WriteArray(a)
                #         gdalfun.gdal_set_ndv(tmp_ds, ndv = 0)

        yield(DatasetFactory(mod=gmrt_fn, data_format=200, mask=swath_mask,
                             #data_format='200:mask={}'.format(swath_mask) if self.swath_only else '200', #mask=self.mask, invert_mask=self.invert_mask,
                             src_srs=self.fetch_module.src_srs, dst_srs=self.dst_srs, x_inc=self.x_inc, y_inc=self.y_inc,
                             weight=self.weight, uncertainty=self.uncertainty, src_region=self.region, parent=self,
                             invert_region = self.invert_region, metadata=copy.deepcopy(self.metadata), cache_dir=self.fetch_module._outdir,
                             verbose=self.verbose)._acquire_module())
        
class GEBCOFetcher(Fetcher):
    """GEBCO Gridded data

Note: Fetches entire zip file.
    """
    
    __doc__ = '''{}
    
    -----------
    Fetches Module: <gebco> - {}'''.format(__doc__, fetches.GEBCO.__doc__)
        
    def __init__(self, exclude_tid = None, **kwargs):
        super().__init__(**kwargs)
        self.exclude_tid = []
        if utils.str_or(exclude_tid) is not None:
            for tid_key in exclude_tid.split('/'):
                self.exclude_tid.append(utils.int_or(tid_key))
        
    def set_ds(self, result):
        wanted_gebco_fns = []
        gebco_fns = utils.p_unzip(
            os.path.join(self.fetch_module._outdir, result[1]),
            exts=['tif'],
            outdir=self.fetch_module._outdir
        )
        
        ## fetch the TID zip if needed
        if self.exclude_tid:
            if fetches.Fetch(
                    self.fetch_module._gebco_urls['gebco_tid']['geotiff'], verbose=self.verbose
            ).fetch_file(os.path.join(self.fetch_module._outdir, 'gebco_tid.zip')) == 0:
                tid_fns = utils.p_unzip(
                    os.path.join(self.fetch_module._outdir, 'gebco_tid.zip'),
                    exts=['tif'],
                    outdir=self.cache_dir
                )

                for tid_fn in tid_fns:
                    ds_config = gdalfun.gdal_infos(tid_fn)
                    
                    if self.region is not None and self.region.valid_p(check_xy=True):
                        inf_region = regions.Region().from_geo_transform(ds_config['geoT'], ds_config['nx'], ds_config['ny'])            
                        inf_region.wmin = self.weight
                        inf_region.wmax = self.weight
                        inf_region.umin = self.uncertainty
                        inf_region.umax = self.uncertainty

                        if regions.regions_intersect_p(inf_region, self.region):
                            wanted_gebco_fns.append(tid_fn)
                    else:
                        wanted_gebco_fns.append(tid_fn)

                for tid_fn in wanted_gebco_fns:
                    utils.echo_msg(tid_fn)
                    tmp_tid = utils.make_temp_fn('tmp_tid.tif', temp_dir=self.fetch_module._outdir)
                    with gdalfun.gdal_datasource(tid_fn) as tid_ds:
                        tid_config = gdalfun.gdal_infos(tid_ds)
                        tid_band = tid_ds.GetRasterBand(1)
                        tid_array = tid_band.ReadAsArray().astype(float)

                    tid_config['ndv'] = -9999
                    tid_config['dt'] = gdal.GDT_Float32                        
                    for tid_key in self.exclude_tid:
                        tid_array[tid_array == tid_key] = tid_config['ndv']
                        
                    for tid_key in self.fetch_module.tid_dic.keys():
                        tid_array[tid_array == tid_key] = self.fetch_module.tid_dic[tid_key][1]
                            
                    gdalfun.gdal_write(tid_array, tmp_tid, tid_config)
                    #utils.echo_msg(tmp_tid)
                    #utils.echo_warning_msg('mask: {}'.format(self.mask))
                    if self.mask is not None:
                        new_mask = utils.make_temp_fn('test_tmp_mask')
                        gdalfun.gdal_mask(tmp_tid, self.mask['mask'], new_mask, msk_value = 1, verbose = True)
                        os.replace(new_mask, tmp_tid)
                        
                    yield(DatasetFactory(mod=tid_fn.replace('tid_', ''), data_format='200:mask={tmp_tid}:weight_mask={tmp_tid}'.format(tmp_tid=tmp_tid),
                                         src_srs=self.fetch_module.src_srs, dst_srs=self.dst_srs, x_inc=self.x_inc, y_inc=self.y_inc,
                                         weight=self.weight, uncertainty=self.uncertainty, src_region=self.region, parent=self,
                                         invert_region=self.invert_region, metadata=copy.deepcopy(self.metadata), cache_dir=self.fetch_module._outdir,
                                         verbose=self.verbose)._acquire_module())
                    
                    utils.remove_glob(tmp_tid)
        else:
            for gebco_fn in gebco_fns:
                ds_config = gdalfun.gdal_infos(gebco_fn)
                inf_region = regions.Region().from_geo_transform(ds_config['geoT'], ds_config['nx'], ds_config['ny'])
                if self.region is not None and self.region.valid_p(check_xy=True):                
                    inf_region.wmin = self.weight
                    inf_region.wmax = self.weight
                    inf_region.umin = self.uncertainty
                    inf_region.umax = self.uncertainty

                    if regions.regions_intersect_p(inf_region, self.region):
                        wanted_gebco_fns.append(gebco_fn)
                else:
                    wanted_gebco_fns.append(gebco_fn)    
            
            for gebco_fn in wanted_gebco_fns:
                yield(DatasetFactory(mod=gebco_fn, data_format=200, src_srs=self.fetch_module.src_srs, dst_srs=self.dst_srs,
                                     x_inc=self.x_inc, y_inc=self.y_inc, weight=self.weight, src_region=self.region, mask=self.mask, 
                                     parent=self, invert_region = self.invert_region, metadata=copy.deepcopy(self.metadata),
                                     cache_dir=self.fetch_module._outdir, verbose=self.verbose)._acquire_module())
                
class CopernicusFetcher(Fetcher):
    """Gridded Copernicus sattelite data.
    """
    
    __doc__ = '''{}
    
    -----------
    Fetches Module: <copernicus> - {}'''.format(__doc__, fetches.CopernicusDEM.__doc__)
    
    def __init__(self, datatype=None, **kwargs):
        super().__init__(**kwargs)
        self.check_size=False
        self.datatype=datatype
        
    def set_ds(self, result):
        if self.datatype is None or result[-1] == self.datatype:
            src_cop_dems = utils.p_unzip(
                os.path.join(self.fetch_module._outdir, result[1]),
                exts=['tif'],
                outdir=self.fetch_module._outdir
            )
            for src_cop_dem in src_cop_dems:
                gdalfun.gdal_set_ndv(src_cop_dem, ndv=0, verbose=False)
                yield(DatasetFactory(mod=src_cop_dem, data_format=200, src_srs=self.fetch_module.src_srs, dst_srs=self.dst_srs, mask=self.mask,
                                     x_inc=self.x_inc, y_inc=self.y_inc, weight=self.weight, uncertainty=self.uncertainty, src_region=self.region,
                                     parent=self, invert_region = self.invert_region, metadata=copy.deepcopy(self.metadata),
                                     cache_dir = self.fetch_module._outdir, verbose=self.verbose, node='pixel')._acquire_module())

class FABDEMFetcher(Fetcher):
    """FABDEM Gridded data
    """
    
    __doc__ = '''{}
    
    -----------
    Fetches Module: <fabdem> - {}'''.format(__doc__, fetches.FABDEM.__doc__)

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def set_ds(self, result):
        src_fab_dems = utils.p_unzip(
            os.path.join(self.fetch_module._outdir, result[1]),
            exts=['tif'],
            outdir=self.fetch_module._outdir
        )
        for src_fab_dem in src_fab_dems:
            gdalfun.gdal_set_ndv(src_fab_dem, ndv=0, verbose=False)
            yield(DatasetFactory(mod=src_fab_dem, data_format=200, src_srs=self.fetch_module.src_srs, dst_srs=self.dst_srs, mask=self.mask,
                                 x_inc=self.x_inc, y_inc=self.y_inc, weight=self.weight, uncertainty=self.uncertainty, src_region=self.region,
                                 parent=self, invert_region = self.invert_region, metadata=copy.deepcopy(self.metadata),
                                 cache_dir = self.fetch_module._outdir, verbose=self.verbose)._acquire_module())

class MarGravFetcher(Fetcher):
    """Marine Gravity Bathymetry
    """
    
    __doc__ = '''{}
    
    -----------
    Fetches Module: <mar_grav> - {}'''.format(__doc__, fetches.MarGrav.__doc__)
    
    def __init__(self, rasterize=False, bathy_only=False, upper_limit=None, lower_limit=None, **kwargs):
        super().__init__(**kwargs)
        self.rasterize = rasterize
        self.bathy_only = bathy_only
        self.upper_limit = utils.float_or(upper_limit)
        self.lower_limit = utils.float_or(lower_limit)

        if self.bathy_only:
            self.upper_limit = 0
        
    def set_ds(self, result):
        mg_region = self.region.copy()
        mg_region.zmax = self.upper_limit
        mg_region.zmin = self.lower_limit
            
        if self.rasterize:
            from cudem import waffles
            mar_grav_fn = utils.make_temp_fn('mar_grav')
            _raster = waffles.WaffleFactory(mod='IDW:min_points=16', data=['{},168:x_offset=REM,1'.format(os.path.join(self.fetch_module._outdir, result[1]))],
                                            src_region=mg_region, xinc=utils.str2inc('30s'), yinc=utils.str2inc('30s'), upper_limit = self.upper_limit,
                                            name=mar_grav_fn, node='pixel', verbose=True)._acquire_module()()
            if self.upper_limit is not None or self.lower_limit is not None:
                ds = gdal.Open(_raster.fn)
                ds_band = ds.GetRasterBand(1)
                ds_arr = ds_band.ReadAsArray()
                if self.upper_limit is not None:
                    ds_arr[ds_arr >= self.upper_limit] = ds_band.GetNoDataValue()

                if self.lower_limit is not None:
                    ds_arr[ds_arr <= self.lower_limit] = ds_band.GetNoDataValue()
                ds = None
                
            yield(DatasetFactory(mod=_raster.fn, data_format=200, src_srs=self.fetch_module.src_srs, dst_srs=self.dst_srs, mask=self.mask,
                                 x_inc=self.x_inc, y_inc=self.y_inc, weight=self.weight, uncertainty=self.uncertainty, src_region=mg_region,
                                 parent=self, invert_region = self.invert_region, metadata=copy.deepcopy(self.metadata),
                                 cache_dir = self.fetch_module._outdir, verbose=self.verbose)._acquire_module())
        else:
            yield(DatasetFactory(mod=os.path.join(self.fetch_module._outdir, result[1]), data_format='168:x_offset=REM', src_srs=self.fetch_module.src_srs, dst_srs=self.dst_srs,
                                 x_inc=self.x_inc, y_inc=self.y_inc, weight=self.weight, uncertainty=self.uncertainty, src_region=mg_region, mask=self.mask,
                                 parent=self, invert_region = self.invert_region, metadata=copy.deepcopy(self.metadata),
                                 cache_dir = self.fetch_module._outdir, verbose=self.verbose)._acquire_module())

class ChartsFetcher(Fetcher):
    """NOAA ENC Charts Fetcher

Digital Soundings
    """
    
    __doc__ = '''{}
    
    -----------
    Fetches Module: <charts> - {}'''.format(__doc__, fetches.Charts.__doc__)

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def set_ds(self, result):
        src_000s = utils.p_unzip(
            os.path.join(self.fetch_module._outdir, result[1]),
            exts=['000'],
            outdir=self.fetch_module._outdir
        )
        for src_000 in src_000s:
            usace_ds = DatasetFactory(mod=src_000, data_format="302:ogr_layer=SOUNDG:z_scale=-1", src_srs=self.fetch_module.src_srs, dst_srs=self.dst_srs,
                                      x_inc=self.x_inc, y_inc=self.y_inc, weight=self.weight, uncertainty=self.uncertainty, src_region=self.region,
                                      parent=self, invert_region=self.invert_region, metadata=copy.deepcopy(self.metadata), mask=self.mask, 
                                      cache_dir=self.fetch_module._outdir, verbose=self.verbose)._acquire_module()
            yield(usace_ds)

class MBSFetcher(Fetcher):
    """NOAA Multibeam Fetcher
    """

    __doc__ = '''{}
    
    -----------
    Fetches Module: <multibeam> - {}'''.format(__doc__, fetches.Multibeam.__doc__)

    def __init__(self, mb_exclude = 'A', want_binned = False, **kwargs):
        super().__init__(**kwargs)
        self.mb_exclude = mb_exclude
        self.want_binned = want_binned

    def set_ds(self, result):            
        mb_infos = self.fetch_module.parse_entry_inf(result, keep_inf=True)
        #if mb_infos[2] != '32':
        ds = DatasetFactory(mod=os.path.join(self.fetch_module._outdir, result[1]), data_format=self.fetch_module.data_format, weight=self.weight,
                            parent=self, src_region=self.region, invert_region=self.invert_region, metadata=copy.deepcopy(self.metadata),
                            mask=self.mask, uncertainty=self.uncertainty, x_inc=self.x_inc, y_inc=self.y_inc, want_binned=self.want_binned,
                            src_srs=self.fetch_module.src_srs, dst_srs=self.dst_srs, verbose=self.verbose, cache_dir=self.fetch_module._outdir,
                            remote=True)._acquire_module()

        yield(ds)
                
class HydroNOSFetcher(Fetcher):
    """NOAA HydroNOS Data Fetcher
    """
    
    __doc__ = '''{}
    
    -----------
    Fetches Module: <hydronos> - {}'''.format(__doc__, fetches.HydroNOS.__doc__)

    def __init__(self, explode=False, **kwargs):
        super().__init__(**kwargs)
        self.explode=explode

    def set_ds(self, result):
        if result[2] == 'xyz':
            nos_fns = utils.p_unzip(
                os.path.join(self.fetch_module._outdir, result[1]),
                exts=['xyz', 'dat'],
                outdir=os.path.dirname(os.path.join(self.fetch_module._outdir, result[1]))
            )
            for nos_fn in nos_fns:
                yield(DatasetFactory(mod=nos_fn, data_format='168:skip=1:xpos=2:ypos=1:zpos=3:z_scale=-1:delimiter=,', src_srs='epsg:4326+5866', dst_srs=self.dst_srs,
                                     x_inc=self.x_inc, y_inc=self.y_inc, weight=self.weight, uncertainty=self.uncertainty, src_region=self.region,
                                     parent=self, invert_region=self.invert_region, metadata=copy.deepcopy(self.metadata), mask=self.mask, 
                                     cache_dir=self.fetch_module._outdir, verbose=self.verbose)._acquire_module())

        elif result[2] == 'bag':
            bag_fns = utils.p_unzip(
                os.path.join(self.fetch_module._outdir, result[1]),
                exts=['bag'],
                outdir=os.path.dirname(os.path.join(self.fetch_module._outdir, result[1]))
            )
            for bag_fn in bag_fns:
                if 'ellipsoid' not in bag_fn.lower():
                    yield(DatasetFactory(mod=bag_fn, data_format=201, explode=self.explode, src_srs=None, dst_srs=self.dst_srs,
                                         x_inc=self.x_inc, y_inc=self.y_inc, weight=self.weight, uncertainty=self.uncertainty, src_region=self.region,
                                         parent=self, invert_region=self.invert_region, metadata=copy.deepcopy(self.metadata), mask=self.mask, 
                                         cache_dir=self.fetch_module._outdir, verbose=self.verbose)._acquire_module())

class CSBFetcher(Fetcher):
    """Crowd Sourced Bathymetry data fetcher
    """
    
    __doc__ = '''{}
    
    -----------
    Fetches Module: <hydronos> - {}'''.format(__doc__, fetches.CSB.__doc__)

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def set_ds(self, result):

        csb_fn = os.path.join(self.fetch_module._outdir, result[1])
        yield(DatasetFactory(mod=csb_fn, data_format='168:skip=1:xpos=2:ypos=3:zpos=4:z_scale=-1:delimiter=,', src_srs='epsg:4326+5866', dst_srs=self.dst_srs,
                             x_inc=self.x_inc, y_inc=self.y_inc, weight=self.weight, uncertainty=self.uncertainty, src_region=self.region,
                             parent=self, invert_region=self.invert_region, metadata=copy.deepcopy(self.metadata), mask=self.mask, 
                             cache_dir=self.fetch_module._outdir, verbose=self.verbose)._acquire_module())

class EMODNetFetcher(Fetcher):
    """EMODNet Data Fetcher
    """
    
    __doc__ = '''{}
    
    -----------
    Fetches Module: <emodnet> - {}'''.format(__doc__, fetches.EMODNet.__doc__)

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def set_ds(self, result):
        if result[2] == 'csv':
            yield(DatasetFactory(mod=os.path.join(self.fetch_module._outdir, result[1]), data_format='168:skip=1:xpos=2:ypos=1:zpos=3:delimiter=,',
                                 src_srs='epsg:4326', dst_srs=self.dst_srs, x_inc=self.x_inc, y_inc=self.y_inc, weight=self.weight,
                                 uncertainty=self.uncertainty, src_region=self.region, mask=self.mask, parent=self, invert_region = self.invert_region,
                                 metadata=copy.deepcopy(self.metadata), cache_dir = self.fetch_module._outdir, verbose=self.verbose)._acquire_module())

        elif result[2] == 'nc':
            yield(DatasetFactory(mod=os.path.join(self.fetch_module._outdir, result[1]), data_format=200,
                                 src_srs='epsg:4326', dst_srs=self.dst_srs, x_inc=self.x_inc, y_inc=self.y_inc, weight=self.weight,
                                 uncertainty=self.uncertainty, src_region=self.region, mask=self.mask, parent=self, invert_region = self.invert_region,
                                 metadata=copy.deepcopy(self.metadata), cache_dir = self.fetch_module._outdir, verbose=self.verbose)._acquire_module())
            

class eHydroFetcher(Fetcher):
    """USACE eHydro soundings
    """

    __doc__ = '''{}
    
    -----------
    Fetches Module: <ehydro> - {}'''.format(__doc__, fetches.eHydro.__doc__)
        
    def __init__(self, want_contours = False, **kwargs):
        super().__init__(**kwargs)
        self.want_contours = want_contours
        
    def set_ds(self, result):
        try:
            src_gdb = utils.gdb_unzip(os.path.join(self.fetch_module._outdir, result[1]), outdir=self.fetch_module._outdir, verbose=False)
        except Exception as e:
            utils.echo_error_msg('{}: {}'.format(self.fn, e))
            src_gdb = None

        if src_gdb is not None:
            tmp_gdb = ogr.Open(src_gdb)
            tmp_layer = tmp_gdb.GetLayer('SurveyPoint')
            src_srs = tmp_layer.GetSpatialRef()
            src_epsg = gdalfun.osr_parse_srs(src_srs)
            elev_datum = tmp_layer[1].GetField('elevationDatum')
            tmp_gdb = tmp_layer = None
            v = vdatums.get_vdatum_by_name(elev_datum)
            
            usace_ds = DatasetFactory(mod=src_gdb, data_format="302:ogr_layer=SurveyPoint_HD:elev_field=Z_label:z_scale=-0.3048",
                                      src_srs='{}+{}'.format(src_epsg, v if v is not None else '5866') if src_epsg is not None else None, dst_srs=self.dst_srs,
                                      x_inc=self.x_inc, y_inc=self.y_inc, weight=self.weight, uncertainty=self.uncertainty, src_region=self.region,
                                      parent=self, invert_region = self.invert_region, metadata=copy.deepcopy(self.metadata), mask=self.mask, 
                                      cache_dir = self.fetch_module._outdir, verbose=self.verbose)._acquire_module()
            yield(usace_ds)

            if self.want_contours:
                usace_ds = DatasetFactory(mod=src_gdb, data_format="302:ogr_layer=ElevationContour_ALL:elev_field=contourElevation:z_scale=-0.3048",
                                          src_srs='{}+{}'.format(src_epsg, v if v is not None else '5866') if src_epsg is not None else None, dst_srs=self.dst_srs,
                                          x_inc=self.x_inc, y_inc=self.y_inc, weight=self.weight, uncertainty=self.uncertainty, src_region=self.region,
                                          parent=self, invert_region = self.invert_region, metadata=copy.deepcopy(self.metadata), mask=self.mask, 
                                          cache_dir = self.fetch_module._outdir, verbose=self.verbose)._acquire_module()
                yield(usace_ds)

    def set_ds_XYZ(self, result):
        src_gdb = utils.gdb_unzip(os.path.join(self.fetch_module._outdir, result[1]), outdir=self.fetch_module._outdir, verbose=False)
        if src_gdb is not None:
            tmp_gdb = ogr.Open(src_gdb)
            tmp_layer = tmp_gdb.GetLayer('SurveyPoint')
            src_srs = tmp_layer.GetSpatialRef()
            src_epsg = gdalfun.osr_parse_srs(src_srs)
            tmp_gdb = None

            src_usaces = utils.p_unzip(os.path.join(self.fetch_module._outdir, result[1]), ['XYZ', 'xyz', 'dat'], outdir=self.fetch_module._outdir)
            for src_usace in src_usaces:
                usace_ds = DatasetFactory(mod=src_usace, data_format='168:z_scale=.3048', src_srs='{}+5866'.format(src_epsg) if src_epsg is not None else None,
                                          dst_srs=self.dst_srs, x_inc=self.x_inc, y_inc=self.y_inc, weight=self.weight, uncertainty=self.uncertainty, mask=self.mask, 
                                          src_region=self.region, parent=self, invert_region=self.invert_region,
                                          metadata=copy.deepcopy(self.metadata), cache_dir = self.fetch_module._outdir, verbose=self.verbose)._acquire_module()
                yield(usace_ds)
        
        
class BlueTopoFetcher(Fetcher):
    """BlueTopo Gridded bathymetric data Fetcher

    -----------
    Parameters:
    
    want_interpolation: True/False to include interpolated cells
    unc_weights: use the uncertainty mask as weights
    """

    __doc__ = '''{}
    
    -----------
    Fetches Module: <bluetopo> - {}'''.format(__doc__, fetches.BlueTopo.__doc__)

    def __init__(self, want_interpolation=False, unc_weights=False, **kwargs):
        super().__init__(**kwargs)
        self.want_interpolation = want_interpolation
        self.unc_weights = unc_weights

    def set_ds(self, result):
        sid = None
        if not self.want_interpolation:
            sid = gdalfun.gdal_extract_band(
                os.path.join(self.fetch_module._outdir, result[1]),
                utils.make_temp_fn('tmp_bt_tid.tif', self.fetch_module._outdir),
                band=3,
                exclude=[0]
            )[0]

        if self.mask is not None:
            new_mask = utils.make_temp_fn('test_tmp_mask')
            gdalfun.gdal_mask(sid, self.mask['mask'], new_mask, msk_value = 1, verbose = True)
            os.replace(new_mask, sid)
            
        yield(
            DatasetFactory(
                mod=os.path.join(self.fetch_module._outdir, result[1]),
                data_format='200:band_no=1:mask={}:uncertainty_mask=2{}'.format(sid, ':weight_mask=2' if self.unc_weights else ''),
                src_srs=self.fetch_module.src_srs, dst_srs=self.dst_srs, 
                x_inc=self.x_inc, y_inc=self.y_inc, weight=self.weight, uncertainty=self.uncertainty, src_region=self.region,
                parent=self, invert_region=self.invert_region, metadata=copy.deepcopy(self.metadata),
                cache_dir = self.fetch_module._outdir, verbose=self.verbose
            )._acquire_module()
        ) 

class NGSFetcher(Fetcher):
    """NGS Monument data
    """

    __doc__ = '''{}
    
    -----------
    Fetches Module: <ngs> - {}'''.format(__doc__, fetches.NGS.__doc__)
    
    def __init__(self, datum = 'geoidHt', **kwargs):
        super().__init__(**kwargs)
        self.datum = datum
        if self.datum not in ['orthoHt', 'geoidHt', 'z', 'ellipHeight']:
            utils.echo_warning_msg('could not parse {}, falling back to geoidHt'.format(datum))
            self.datum = 'geoidHt'

    def set_ds(self, result):
        with open(os.path.join(self.fetch_module._outdir, result[1]), 'r') as json_file:
            r = json.load(json_file)
            if len(r) > 0:
                with open(os.path.join(self.fetch_module._outdir, '_tmp_ngs.xyz'), 'w') as tmp_ngs:
                    for row in r:
                        z = utils.float_or(row[self.datum])
                        if z is not None:
                            xyz = xyzfun.XYZPoint(src_srs='epsg:4326').from_list(
                                [float(row['lon']), float(row['lat']), z]
                            )
                            xyz.dump(dst_port=tmp_ngs)

        yield(DatasetFactory(mod=os.path.join(self.fetch_module._outdir, '_tmp_ngs.xyz'), data_format=168, src_srs='epsg:4326', dst_srs=self.dst_srs,
                             x_inc=self.x_inc, y_inc=self.y_inc, weight=self.weight, uncertainty=self.uncertainty, src_region=self.region,
                             parent=self, invert_region=self.invert_region, metadata=copy.deepcopy(self.metadata), mask=self.mask, 
                             cache_dir = self.fetch_module._outdir, verbose=self.verbose)._acquire_module())

class TidesFetcher(Fetcher):
    """NOS Tide Station data
    """

    __doc__ = '''{}
        
    -----------
    Fetches Module: <tides> - {}'''.format(__doc__, fetches.Tides.__doc__)

    def __init__(self, s_datum='mllw', t_datum='msl', units='m', **kwargs):
        super().__init__(**kwargs)
        self.s_datum = s_datum
        self.t_datum = t_datum
        self.units = units
        
    def set_ds(self, result):
        with open(os.path.join(self.fetch_module._outdir, result[1]), 'r') as json_file:
            r = json.load(json_file)
            if len(r) > 0:
                with open(os.path.join(self.fetch_module._outdir, '_tmp_tides.xyz'), 'w') as tmp_ngs:
                    for feature in r['features']:
                        if self.fetch_module.station_id is not None:
                            if self.fetch_module.station_id != feature['attributes']['id']:
                                continue
                        lon = feature['attributes']['longitude']
                        lat = feature['attributes']['latitude']
                        if feature['attributes'][self.s_datum] != -99999.99 and feature['attributes'][self.t_datum] != -99999.99:
                            z = feature['attributes'][self.s_datum] - feature['attributes'][self.t_datum]
                            if self.units == 'm':
                                z = z * 0.3048

                            xyz = xyzfun.XYZPoint(src_srs='epsg:4326').from_list([lon, lat, z])
                            xyz.dump(dst_port=tmp_ngs)

        yield(DatasetFactory(mod=os.path.join(self.fetch_module._outdir, '_tmp_tides.xyz'), data_format=168, src_srs='epsg:4326', dst_srs=self.dst_srs,
                             x_inc=self.x_inc, y_inc=self.y_inc, weight=self.weight, uncertainty=self.uncertainty, src_region=self.region, 
                             parent=self, invert_region=self.invert_region, metadata=copy.deepcopy(self.metadata), mask=self.mask, 
                             cache_dir = self.fetch_module._outdir, verbose=self.verbose)._acquire_module())

class WaterServicesFetcher(Fetcher):
    """USGS Water Services

    -----------
    Parameters:
    
    site_code: the site code to fetch
    units: 'm' for meters

site_codes:
    00065 - Gate Height
    00060 - StreamFlow
    63160 - Stream water level elevation above NAVD 1988
    62611 - Groundwater level above NAVD 1988
    72019 - Depth to water level, units below land surface
    """

    __doc__ = '''{}
    
    -----------
    Fetches Module: <waterservices> - {}'''.format(__doc__, fetches.WaterServices.__doc__)
    
    def __init__(self, site_code='00065', units='m', **kwargs):
        super().__init__(**kwargs)
        self.units = units
        self.site_code = site_code
        
    def set_ds(self, result):
        with open(os.path.join(self.fetch_module._outdir, result[1]), 'r') as json_file:
            r = json.load(json_file)
            if len(r) > 0:
                with open(os.path.join(self.fetch_module._outdir, '_tmp_ws.xyz'), 'w') as tmp_ws:
                    features = r['value']['timeSeries']
                    for feature in features:
                        if feature['variable']['variableCode'][0]['value'] == self.site_code:
                            lon = float(feature['sourceInfo']['geoLocation']['geogLocation']['longitude'])
                            lat = float(feature['sourceInfo']['geoLocation']['geogLocation']['latitude'])
                            z = float(feature['values'][0]['value'][0]['value'])

                            if self.units == 'm':
                                z = z * 0.3048
                            
                            xyz = xyzfun.XYZPoint(src_srs='epsg:4326').from_list([lon, lat, z])
                            xyz.dump(dst_port=tmp_ws)
                        
        yield(DatasetFactory(mod=os.path.join(self.fetch_module._outdir, '_tmp_ws.xyz'), data_format=168, src_srs='epsg:4326', dst_srs=self.dst_srs,
                             x_inc=self.x_inc, y_inc=self.y_inc, weight=self.weight, uncertainty=self.uncertainty, src_region=self.region, 
                             parent=self, invert_region=self.invert_region, metadata=copy.deepcopy(self.metadata), mask=self.mask, 
                             cache_dir = self.fetch_module._outdir, verbose=self.verbose)._acquire_module())

class VDatumFetcher(Fetcher):
    """VDatum transformation grids
    """

    __doc__ = '''{}
    
    -----------
    Fetches Module: <vdatum> - {}'''.format(__doc__, fetches.VDATUM.__doc__)
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def set_ds(self, result):
        src_tif = os.path.join(self.fetch_module._outdir, '{}'.format(utils.fn_basename2(os.path.basename(result[1]))))
        
        if result[1].endswith('.zip'):
            v_gtx = utils.p_f_unzip(os.path.join(self.fetch_module._outdir, result[1]), [result[2]], outdir=self.fetch_module._outdir)[0]
            utils.run_cmd('gdalwarp {} {} --config CENTER_LONG 0'.format(v_gtx, src_tif), verbose=self.verbose)
        
        yield(DatasetFactory(mod=src_tif, data_format=200, node='pixel', src_srs=self.fetch_module.src_srs, dst_srs=self.dst_srs,
                             x_inc=self.x_inc, y_inc=self.y_inc, weight=self.weight, uncertainty=self.uncertainty, src_region=self.region,
                             parent=self, invert_region=self.invert_region, metadata=copy.deepcopy(self.metadata), mask=self.mask, 
                             cache_dir = self.fetch_module._outdir, verbose=self.verbose)._acquire_module())        

## todo: allow lakes bathymetry
## as well as lakes breaklines (shape nodes)
## see: https://www.esri.com/arcgis-blog/products/arcgis-pro/3d-gis/hydro-flattening-of-river-shorelines-in-lidar-based-dem-production/
class HydroLakesFetcher(Fetcher):
    """HydroLakes lake bathymetric data
    """
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def set_ds(self, result):
        pass
    
class DatasetFactory(factory.CUDEMFactory):
    """Dataset Factory Settings and Generator
    
    Parse a datalist entry and return the dataset object
    """
    
    _modules = {
        -1: {'name': 'datalist', 'fmts': ['datalist', 'mb-1', 'dl'], 'call': Datalist},
        -2: {'name': 'zip', 'fmts': ['zip', 'ZIP'], 'call': ZIPlist}, # add other archive formats (gz, tar.gz, 7z, etc.)
        -3: {'name': 'scratch', 'fmts': [''], 'call': Scratch },
        167: {'name': 'yxz', 'fmts': ['yxz'], 'call': YXZFile},
        168: {'name': 'xyz', 'fmts': ['xyz', 'csv', 'dat', 'ascii', 'txt'], 'call': XYZFile},
        200: {'name': 'gdal', 'fmts': ['tif', 'tiff', 'img', 'grd', 'nc', 'vrt'], 'call': GDALFile},
        201: {'name': 'bag', 'fmts': ['bag'], 'call': BAGFile},
        202: {'name': 'swot_pixc', 'fmts': ['h5'], 'call': SWOT_PIXC},
        203: {'name': 'swot_hr_raster', 'fmts': ['nc'], 'call': SWOT_HR_Raster},
        300: {'name': 'las', 'fmts': ['las', 'laz'], 'call': LASFile},
        301: {'name': 'mbs', 'fmts': ['fbt', 'mb'], 'call': MBSParser},
        302: {'name': 'ogr', 'fmts': ['000', 'shp', 'geojson', 'gpkg', 'gdb/'], 'call': OGRFile},
        ## fetches modules
        -100: {'name': 'https', 'fmts': ['https'], 'call': Fetcher},
        -101: {'name': 'gmrt', 'fmts': ['gmrt'], 'call': GMRTFetcher},
        -102: {'name': 'gebco', 'fmts': ['gebco'], 'call': GEBCOFetcher},
        -103: {'name': 'copernicus', 'fmts': ['copernicus'], 'call': CopernicusFetcher},
        -104: {'name': 'fabdem', 'fmts': ['fabdem'], 'call': FABDEMFetcher},
        -105: {'name': 'nasadem', 'fmts': ['nasadem'], 'call': Fetcher},
        -106: {'name': 'mar_grav', 'fmts': ['mar_grav'], 'call': MarGravFetcher},
        -107: {'name': 'srtm_plus', 'fmts': ['srtm_plus'], 'call': Fetcher},
        -200: {'name': 'charts', 'fmts': ['charts'], 'call': ChartsFetcher},
        -201: {'name': 'multibeam', 'fmts': ['multibeam'], 'call': MBSFetcher},
        -202: {'name': 'hydronos', 'fmts': ['hydronos'], 'call': HydroNOSFetcher},
        -203: {'name': 'ehydro', 'fmts': ['ehydro'], 'call': eHydroFetcher},
        -204: {'name': 'bluetopo', 'fmts': ['bluetopo'], 'call': BlueTopoFetcher},
        -205: {'name': 'ngs', 'fmts': ['ngs'], 'call': NGSFetcher},
        -206: {'name': 'tides', 'fmts': ['tides'], 'call': TidesFetcher},
        -207: {'name': 'digital_coast', 'fmts': ['digital_coast'], 'call': Fetcher},
        -208: {'name': 'ncei_thredds', 'fmts': ['ncei_thredds'], 'call': Fetcher},
        -209: {'name': 'tnm', 'fmts': ['tnm'], 'call': Fetcher},
        -210: {'name': "CUDEM", 'fmts': ['CUDEM'], 'call': Fetcher},
        -211: {'name': "CoNED", 'fmts': ['CoNED'], 'call': DAVFetcher_CoNED},
        -212: {'name': "SLR", 'fmts': ['SLR'], 'call': DAVFetcher_SLR},
        -213: {'name': 'waterservies', 'fmts': ['waterservices'], 'call': WaterServicesFetcher},
        -214: {'name': "icesat", 'fmts': ['icesat'], 'call': IceSatFetcher},
        -215: {'name': 'ned', 'fmts': ['ned', 'ned1'], 'call': NEDFetcher},        
        -216: {'name': "swot", 'fmts': ['swot'], 'call': SWOTFetcher},
        -217: {'name': "csb", 'fmts': ['csb'], 'call': CSBFetcher},
        -300: {'name': 'emodnet', 'fmts': ['emodnet'], 'call': EMODNetFetcher},
        -301: {'name': 'chs', 'fmts': ['chs'], 'call': Fetcher}, # chs is broken
        -302: {'name': 'hrdem', 'fmts': ['hrdem'], 'call': Fetcher},
        -303: {'name': 'arcticdem', 'fmts': ['arcticdem'], 'call': Fetcher},
        -500: {'name': 'vdatum', 'fmts': ['vdatum'], 'call': VDatumFetcher},
        -600: {'name': 'hydrolakes', 'fmts': ['hydrolakes'], 'call': HydroLakesFetcher},        
    }
    _datalist_cols = ['path', 'format', 'weight', 'uncertainty', 'title', 'source',
                      'date', 'type', 'resolution', 'horz', 'vert',
                      'url']

    _metadata_keys = ['name', 'title', 'source', 'date', 'data_type', 'resolution',
                      'hdatum', 'vdatum', 'url']
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    ## redefine the factory default _parse_mod function for datasets
    ## the mod in this case is a datalist entry and the format key
    ## becomes the module
    def _parse_mod(self, mod=None):
        """parse the datalist entry line"""

        self.kwargs['fn'] = mod
        if self.kwargs['fn'] is None:
            return(self)

        ## mod exists as a file, no other entry items should occur, so
        ## guess the format and finish there...
        ## the format number becomes the mod_name
        ## check for specified data format as well
        if os.path.exists(self.kwargs['fn']):
            if self.kwargs['data_format'] is None:
                self.mod_name = self.guess_data_format(self.kwargs['fn'])
                self.mod_args = {}
                self.kwargs['data_format'] = self.mod_name
            else:
                opts = str(self.kwargs['data_format']).split(':')
                if len(opts) > 1:
                    self.mod_name = int(opts[0])
                    self.mod_args = utils.args2dict(list(opts[1:]), {})
                else:
                    self.mod_name = int(self.kwargs['data_format'])
                    self.mod_args = {}
                    
                self.kwargs['data_format'] = self.mod_name

            # inherit metadata from parent if available
            self.kwargs['metadata'] = {}
            self.kwargs['metadata']['name'] = utils.fn_basename2(os.path.basename(self.kwargs['fn']))
            
            for key in self._metadata_keys:
                if key not in self.kwargs['metadata'].keys():
                    self.kwargs['metadata'][key] = None

            return(self.mod_name, self.mod_args)

        ## if fn is not a path, parse it as a datalist entry
		## breaks on path with space, e.g. /meda/user/My\ Passport/etc
        #this_entry = re.findall(r'[^"\s]\S*|".+?"', self.kwargs['fn'].rstrip())
        #this_entry = shlex.split(shlex.quote(self.kwargs['fn'].rstrip()))#.replace("'", "\'")))
        #this_entry = [p for p in re.split("( |\\\".*?\\\"|'.*?')", self.kwargs['fn']) if p.strip()]
        #this_entry = [t.strip('"') for t in re.findall(r'[^\s"]+|"[^"]*"', self.kwargs['fn'].rstrip())]
        this_entry = re.findall("(?:\".*?\"|\S)+", self.kwargs['fn'].rstrip())
        try:
            entry = [utils.str_or(x) if n == 0 \
                     else utils.str_or(x) if n < 2 \
                     else utils.float_or(x) if n < 3 \
                     else utils.float_or(x) if n < 4 \
                     else utils.str_or(x) \
                     for n, x in enumerate(this_entry)]

        except Exception as e:
            utils.echo_error_msg('could not parse entry {}'.format(self.kwargs['fn']))
            return(self)

        ## data format - entry[1]
        ## guess format based on fn if not specified otherwise
        ## parse the format for dataset specific opts.
        #utils.echo_msg(entry)
        if len(entry) < 2:
            if self.kwargs['data_format'] is not None:
                entry.append(self.kwargs['data_format'])
            else:
                for key in self._modules.keys():
                    if entry[0].startswith('http') or entry[0].startswith('/vsicurl/'):
                        see = 'https'
                    else:
                        se = entry[0].split('.')
                        see = se[-1] if len(se) > 1 else entry[0].split(":")[0]

                    if 'fmts' in self._modules[key].keys():
                        if see in self._modules[key]['fmts']:
                            entry.append(int(key))
                            break

            if len(entry) < 2:
                utils.echo_error_msg('could not parse entry {}'.format(self.kwargs['fn']))
                return(self)
        #else:
        opts = str(entry[1]).split(':')
        if len(opts) > 1:
            self.mod_args = utils.args2dict(list(opts[1:]), {})
            entry[1] = int(opts[0])
        else:
            self.mod_args = {}
                                                       
        try:
            assert isinstance(utils.int_or(entry[1]), int)
        except:
            utils.echo_error_msg('could not parse datalist entry {}'.format(entry))
            return(self)

        if entry[0].startswith('http') and int(entry[1]) == 200:
            entry[0] = '/vsicurl/{}'.format(entry[0])
        
        self.kwargs['data_format'] = int(entry[1])
        self.mod_name = int(entry[1])
        
        ## file-name (or fetches module name) - entry[0]
        ## set to relative path from parent
        ## don't set relative path if 'fn' is a
        ## fetches module (entry[1] < -3)
        if 'parent' not in self.kwargs:
            self.kwargs['parent'] = None

        if self.kwargs['parent'] is None:
            self.kwargs['fn'] = entry[0]
            self.kwargs['parent'] = None
        else:
            if self.mod_name >= -2:
                self.kwargs['fn'] = os.path.join(
                    os.path.dirname(self.kwargs['parent'].fn), entry[0]
                )
            else:
                self.kwargs['fn'] = entry[0]
            
        ## weight - entry[2]
        ## inherit weight of parent
        if len(entry) < 3:
            entry.append(self.set_default_weight())
        elif entry[2] is None:
            entry[2] = self.set_default_weight()

        if 'weight' not in self.kwargs:
            self.kwargs['weight'] = 1
            
        if self.kwargs['parent'] is not None:
            if self.kwargs['weight'] is not None:
                self.kwargs['weight'] *= entry[2]
        else:
            if self.kwargs['weight'] is not None:
                self.kwargs['weight'] = entry[2]

        ## uncertainty - entry[3]
        ## combine with partent using root sum of squares
        if len(entry) < 4:
            entry.append(self.set_default_uncertainty())
        elif entry[3] is None:
            entry[3] = self.set_default_uncertainty()

        if 'uncertainty' not in self.kwargs:
            self.kwargs['uncertainty'] = 0
            
        if self.kwargs['parent'] is not None:
            if self.kwargs['uncertainty'] is not None:
                self.kwargs['uncertainty'] = math.sqrt(self.kwargs['uncertainty']**2 + entry[3]**2)
        else:
            if self.kwargs['uncertainty'] is not None:
                self.kwargs['uncertainty'] = entry[3]
                    
        ## Optional arguments follow, for metadata generation
        if 'metadata' not in self.kwargs:
            self.kwargs['metadata'] = {}

        for key in self._metadata_keys:
            if key not in self.kwargs['metadata'].keys():
                self.kwargs['metadata'][key] = None

        ## title - entry[4]
        if len(entry) < 5:
            entry.append(self.kwargs['metadata']['title'])
        else:
            self.kwargs['metadata']['title'] = entry[4]

        if self.kwargs['metadata']['name'] is None:
            self.kwargs['metadata']['name'] = utils.fn_basename2(os.path.basename(self.kwargs['fn']))

        ## source - entry[5]
        if len(entry) < 6:
            entry.append(self.kwargs['metadata']['source'])
        else:
            self.kwargs['metadata']['source'] = entry[5]

        ## date - entry[6]
        if len(entry) < 7:
            entry.append(self.kwargs['metadata']['date'])
        else:
            self.kwargs['metadata']['date'] = entry[6]

        ## data type - entry[7]
        if len(entry) < 8:
            entry.append(self.kwargs['metadata']['data_type'])
        else:
            self.kwargs['metadata']['data_type'] = entry[7]

        ## resolution - entry[8]
        if len(entry) < 9:
            entry.append(self.kwargs['metadata']['resolution'])
        else:
            self.kwargs['metadata']['resolution'] = entry[8]

        ## hdatum - entry[9]
        if len(entry) < 10:
            entry.append(self.kwargs['metadata']['hdatum'])
        else:
            self.kwargs['metadata']['hdatum'] = entry[9]

        ## vdatum - entry[10]
        if len(entry) < 11:
            entry.append(self.kwargs['metadata']['vdatum'])
        else:
            self.kwargs['metadata']['vdatum'] = entry[10]

        ## url - entry[11]
        if len(entry) < 12:
            entry.append(self.kwargs['metadata']['url'])
        else:
            self.kwargs['metadata']['url'] = entry[11]

        return(self.mod_name, self.mod_args)

    def set_default_weight(self):
        return(1)

    def set_default_uncertainty(self):
        return(0)
    
    def guess_data_format(self, fn):
        """guess a data format based on the file-name"""
        
        for key in self._modules.keys():
            if fn.split('.')[-1] in self._modules[key]['fmts']:
                return(key)
            ## hack to accept .mb* mb-system files without having to record every one...
            elif fn.split('.')[-1][:2] in self._modules[key]['fmts']:
                return(key)
            
    def write_parameter_file(self, param_file: str):
        try:
            with open(param_file, 'w') as outfile:
                json.dump(self.__dict__, outfile)
                utils.echo_msg('New DatasetFactory file written to {}'.format(param_file))
                
        except:
            raise ValueError('DatasetFactory: Unable to write new parameter file to {}'.format(param_file))            

## ==============================================
## Command-line Interface (CLI)
## $ dlim
##
## datalists cli
## ==============================================
datalists_usage = """{cmd} ({dl_version}): DataLists IMproved; Process and generate datalists

usage: {cmd} [ -acdghijnquwEJPRT [ args ] ] DATALIST,FORMAT,WEIGHT,UNCERTAINTY ...

Options:
  -R, --region\t\tRestrict processing to the desired REGION 
\t\t\tWhere a REGION is xmin/xmax/ymin/ymax[/zmin/zmax[/wmin/wmax/umin/umax]]
\t\t\tUse '-' to indicate no bounding range; e.g. -R -/-/-/-/-10/10/1/-/-/-
\t\t\tOR an OGR-compatible vector file with regional polygons. 
\t\t\tWhere the REGION is /path/to/vector[:zmin/zmax[/wmin/wmax/umin/umax]].
\t\t\tIf a vector file is supplied, will use each region found therein.
  -E, --increment\tBlock data to INCREMENT in native units.
\t\t\tWhere INCREMENT is x-inc[/y-inc]
  -X, --extend\t\tNumber of cells with which to EXTEND the output DEM REGION and a 
\t\t\tpercentage to extend the processing REGION.
\t\t\tWhere EXTEND is dem-extend(cell-count)[:processing-extend(percentage)]
\t\t\te.g. -X6:10 to extend the DEM REGION by 6 cells and the processing region by 10 
\t\t\tpercent of the input REGION.
  -J, --s_srs\t\tSet the SOURCE projection.
  -P, --t_srs\t\tSet the TARGET projection. (REGION should be in target projection) 
  -D, --cache-dir\tCACHE Directory for storing temp and output data.
  -Z, --z-precision\tSet the target precision of dumped z values. (default is 4)
  -A, --stack-mode\tSet the STACK MODE to 'mean', 'min', 'max' or 'supercede' (with -E and -R)
  -T, --stack_filter\tFILTER the data stack using one or multiple filters. 
\t\t\tWhere FILTER is fltr_name[:opts] (see `grits --modules` for more information)
\t\t\tThe -T switch may be set multiple times to perform multiple filters.
\t\t\tAvailable FILTERS: {grits_modules}
  -F, --point_filter\tFILTER the POINT data using one or multiple filters. 
\t\t\tWhere FILTER is fltr_name[:opts] 
\t\t\tThe -F switch may be set multiple times to perform multiple filters.
\t\t\tAvailable FILTERS: {point_filter_modules}

  --mask\t\tMASK the datalist to the given REGION/INCREMENTs
  --spatial-metadata\tGenerate SPATIAL METADATA of the datalist to the given REGION/INCREMENTs
  --archive\t\tARCHIVE the datalist to the given REGION[/INCREMENTs]
  --glob\t\tGLOB the datasets in the current directory to stdout
  --info\t\tGenerate and return an INFO dictionary of the dataset
  --weights\t\tOutput WEIGHT values along with xyz
  --uncertainties\tOutput UNCERTAINTY values along with xyz
  --stack-node\t\tOutput stacked x/y data rather than pixel
  --quiet\t\tLower the verbosity to a quiet

  --modules\t\tDisplay the datatype descriptions and usage
  --help\t\tPrint the usage text
  --version\t\tPrint the version information

Supported datalist formats (see {cmd} --modules <dataset-key> for more info): 
  {dl_formats}

Examples:
  % {cmd} my_data.datalist -R -90/-89/30/31
  % {cmd} -R-90/-89/30/31/-100/100 *.tif -l -w > tifs_in_region.datalist
  % {cmd} -R my_region.shp my_data.xyz -w -s_srs epsg:4326 -t_srs epsg:3565 > my_data_3565.xyz
""".format(cmd=os.path.basename(sys.argv[0]), 
           dl_version=cudem.__version__,
           dl_formats=utils._cudem_module_name_short_desc(DatasetFactory._modules),
           grits_modules=factory._cudem_module_short_desc(grits.GritsFactory._modules),
           point_filter_modules=factory._cudem_module_short_desc(PointFilterFactory._modules))

def datalists_cli(argv=sys.argv):
    """run datalists from command-line

    See `datalists_cli_usage` for full cli options.
    """

    dls = []
    src_srs = None
    dst_srs = None
    i_regions = []
    these_regions = []
    xy_inc = [None, None]
    extend = 0
    want_weights = False
    want_uncertainties = False
    want_mask = False
    want_inf = False
    want_list = False
    want_glob = False
    want_archive = False
    want_verbose = True
    want_region = False
    want_separate = False
    want_sm = False
    invert_region = False
    z_precision = 4
    stack_fltrs = []
    pnt_fltrs = []
    stack_node = False
    stack_mode = 'mean'
    cache_dir = utils.cudem_cache()
    
    ## parse command line arguments.
    i = 1
    while i < len(argv):
        arg = argv[i]
        if arg == '--region' or arg == '-R':
            i_regions.append(str(argv[i + 1]))
            i = i + 1
        elif arg[:2] == '-R':
            i_regions.append(str(arg[2:]))
        elif arg == '--increment' or arg == '-E':
            xy_inc = argv[i + 1].split('/')
            i = i + 1
        elif arg[:2] == '-E':
            xy_inc = arg[2:].split('/')
        elif arg == '--extend' or arg == '-X':
            extend = utils.int_or(argv[i + 1], 0)
            i += 1
        elif arg[:2] == '-X':
            extend = utils.int_or(argv[2:], 0)
        elif arg == '-s_srs' or arg == '--s_srs' or arg == '-J':
            src_srs = argv[i + 1]
            i = i + 1
        elif arg == '-t_srs' or arg == '--t_srs' or arg == '-P':
            dst_srs = argv[i + 1]
            i = i + 1
        elif arg == '-z_precision' or arg == '--z-precision' or arg == '-Z':
            z_precision = utils.int_or(argv[i + 1], 4)
            i = i + 1
        elif arg[:2] == '-Z':
            z_precision = utils.int_or(argv[2:], 4)

        elif arg == '-stack-mode' or arg == '--stack-mode' or arg == '-A':
            stack_mode = utils.str_or(argv[i + 1], 'mean')
            i = i + 1
        elif arg[:2] == '-A':
            stack_mode = utils.str_or(arg[2:], 'mean')
            
        elif arg == '-stack-filter' or arg == '--stack-filter' or arg == '-T':
            stack_fltrs.append(argv[i + 1])
            i = i + 1
        elif arg[:2] == '-T':
            stack_fltrs.append(argv[2:])
            
        elif arg == '-point-filter' or arg == '--point-filter' or arg == '-F':
            pnt_fltrs.append(argv[i + 1])
            i = i + 1
        elif arg[:2] == '-F':
            pnt_fltrs.append(argv[2:])

        elif arg == '--cache-dir' or arg == '-D' or arg == '-cache-dir':
            cache_dir = utils.str_or(argv[i + 1], utils.cudem_cache)
            i = i + 1
        elif arg[:2] == '-D': cache_dir = utils.str_or(argv[i + 1], utils.cudem_cache)
            
        elif arg == '--mask' or arg == '-m':
            want_mask = True
        elif arg == '--invert_region' or arg == '-v':
            invert_region = True
        elif arg == '--archive' or arg == '-a':
            want_archive = True
        elif arg == '--weights' or arg == '-w':
            want_weights = True
        elif arg == '--uncertainties' or arg == '-u':
            want_uncertainties = True
        elif arg == '--info' or arg == '-i':
            want_inf = True
        elif arg == '--region_inf' or arg == '-r':
            want_region = True
        elif arg == '--list' or arg == '-l':
            want_list = True
        elif arg == '--glob' or arg == '-g':
            want_glob = True
        elif arg == '--spatial-metadata' or arg == '-sm':
            want_sm = True
            want_mask = True
        elif arg == '--stack-node' or arg == '-n':
            stack_node = True
            
        elif arg == '--separate' or arg == '-s':
            want_separate = True

        elif arg == '--modules':
            utils.echo_modules(DatasetFactory._modules, None if i+1 >= len(argv) else int(sys.argv[i+1]))
            sys.exit(0)           
        elif arg == '--quiet' or arg == '-q':
            want_verbose = False
        elif arg == '--help' or arg == '-h':
            print(datalists_usage)
            sys.exit(1)
        elif arg == '--version' or arg == '-v':
            print('{}, version {}'.format(
                os.path.basename(sys.argv[0]), cudem.__version__)
                  )
            sys.exit(1)
        elif arg[0] == '-':
            print(datalists_usage)
            sys.exit(0)
        else: dls.append(arg)#'"{}"'.format(arg)) # FIX THIS!!!
        
        i = i + 1

    if len(xy_inc) < 2:
        xy_inc.append(xy_inc[0])
    elif len(xy_inc) == 0:
        xy_inc = [None, None]

    if want_glob:
        import glob
        for key in DatasetFactory()._modules.keys():
            if key != -1 and key != '_factory':
                for f in DatasetFactory()._modules[key]['fmts']:
                    globs = glob.glob('*.{}'.format(f))
                    [sys.stdout.write(
                        '{}\n'.format(
                            ' '.join(
                                [x, str(key), '1', '0']
                            )
                        )
                    ) for x in globs]
                    
        sys.exit(0)

    #stack_fltrs = [f + ':stacks=True' for f in stack_fltrs]
    stack_fltrs = [':'.join(f.split('/')) for f in stack_fltrs]
    pnt_fltrs = [':'.join(f.split('/')) for f in pnt_fltrs]
    if not i_regions: i_regions = [None]
    these_regions = regions.parse_cli_region(i_regions, want_verbose)
    for rn, this_region in enumerate(these_regions):
        ## buffer the region by `extend` if xy_inc is set
        ## this effects the output naming of masks/stacks!
        ## do we want this like in waffles where the output name
        ## does not include the -X extend buffer?
        if xy_inc[0] is not None and xy_inc[1] is not None and this_region is not None:
            this_region.buffer(
                x_bv=(utils.str2inc(xy_inc[0])*extend),
                y_bv=(utils.str2inc(xy_inc[1])*extend)
            )
        if len(dls) == 0:
            sys.stderr.write(datalists_usage)
            utils.echo_error_msg('you must specify some type of data')
        else:
            ## intiialze the input data. Treat data from CLI as a datalist.
            this_datalist = init_data(
                dls, region=this_region, src_srs=src_srs, dst_srs=dst_srs, xy_inc=xy_inc, sample_alg='auto',
                want_weight=want_weights, want_uncertainty=want_uncertainties, want_verbose=want_verbose,
                want_mask=want_mask, want_sm=want_sm, invert_region=invert_region, cache_dir=cache_dir,
                dump_precision=z_precision, pnt_fltrs=pnt_fltrs, stack_fltrs=stack_fltrs, stack_node=stack_node,
                stack_mode=stack_mode
            )

            if this_datalist is not None and this_datalist.valid_p(
                    fmts=DatasetFactory._modules[this_datalist.data_format]['fmts']
            ):
                this_datalist.initialize()                
                if not want_weights:
                    this_datalist.weight = None
                    
                if not want_uncertainties:
                    this_datalist.uncertainty = None
                    
                if want_inf:
                    print(this_datalist.inf()) # output the datalist inf blob
                elif want_list:
                    this_datalist.echo() # output each dataset from the datalist
                elif want_region: # get the region and warp it if necessary
                    this_inf = this_datalist.inf()
                    this_region = regions.Region().from_list(this_inf.minmax)
                    #print(this_region)
                    #print(dst_srs)
                    #print(this_inf.src_srs)
                    if dst_srs is not None:
                        if src_srs is not None:
                            this_region.src_srs = src_srs
                            this_region.warp(dst_srs)
                        elif this_inf.src_srs is not None:
                            this_region.src_srs = this_inf.src_srs
                            this_region.warp(dst_srs, include_z=False)
                    print(this_region.format('gmt'))
                elif want_archive:
                    this_datalist.archive_xyz() # archive the datalist as xyz
                else:
                    #try:
                    if want_separate: # process and dump each dataset independently
                        for this_entry in this_datalist.parse():
                            this_entry.dump_xyz()
                    else: # process and dump the datalist as a whole
                        this_datalist.dump_xyz()
                    # except KeyboardInterrupt:
                    #   utils.echo_error_msg('Killed by user')
                    #   break
                    # except BrokenPipeError:
                    #   utils.echo_error_msg('Pipe Broken')
                    #   break
                    # except Exception as e:
                    #   utils.echo_error_msg(e)
                    #   print(traceback.format_exc())
### End
