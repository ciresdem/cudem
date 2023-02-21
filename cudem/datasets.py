### datasets.py - Datasets
##
## Copyright (c) 2010 - 2023 Regents of the University of Colorado
##
## datasets.py is part of CUDEM
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
## Dataset parsing.
## datasets include: xyz, raster (gdal), las/laz (laspy), mbs (MBSystem), fetches
##
### Code:

import os
import sys
import json
import laspy as lp
import copy
import csv
import math
                        
import numpy as np
from scipy.spatial import ConvexHull

from osgeo import gdal
from osgeo import ogr
from osgeo import osr

from cudem import utils
from cudem import regions
from cudem import xyzfun
from cudem import demfun

## ==============================================
## Elevation Dataset class
## ==============================================
class ElevationDataset():
    """representing an Elevation Dataset
    
    This is the super class for all datalist (dlim) datasets.    
    Each dataset sub-class should define a dataset-specific    
    data parser (to xyz) and a generate_inf function to generate
    inf files. 

    Specifically, each sub-dataset should minimally define the following:
    
    sub_ds.generate_inf()
    sub_ds.yield_xyz()
    
    Where:
    generate_inf() generates a dlim compatible inf file,
    yield_xyz() yields the xyz elevation data (xyzfun.XYZPoint) from the dataset.

    ----
    Parameters:

    """

    gdal_sample_methods = [
        'near', 'bilinear', 'cubic', 'cubicspline', 'lanczos',
        'average', 'mode',  'max', 'min', 'med', 'Q1', 'Q3', 'sum'
    ]
    
    def __init__(
            self,
            fn=None,
            data_format=None,
            weight=1,
            src_srs=None,
            dst_srs='epsg:4326',
            x_inc=None,
            y_inc=None,
            sample_alg='bilinear',
            metadata={
                'name':None,
                'title':None,
                'source':None,
                'date':None,
                'data_type':None,
                'resolution':None,
                'hdatum':None,
                'vdatum':None,
                'url':None
            },
            parent=None,
            src_region=None,
            invert_region=False,
            cache_dir=None,
            verbose=False,
            remote=False
    ):
        self.fn = fn
        self._fn = None
        self.data_format = data_format
        self.weight = weight
        self.src_srs = src_srs
        self.dst_srs = dst_srs
        self.dst_trans = None
        self.trans_region = None
        self.src_trans_srs = None
        self.dst_trans_srs = None
        self.region = src_region
        self.invert_region = invert_region
        self.parent = parent
        self.verbose = verbose
        self.data_entries = []
        self.data_lists = {}        
        self.remote = remote
        self.metadata = copy.deepcopy(metadata)
        self.x_inc = utils.str2inc(x_inc)
        self.y_inc = utils.str2inc(y_inc)
        if sample_alg in self.gdal_sample_methods:
            self.sample_alg = sample_alg
        else:
            utils.echo_warning_msg('{} is not a valid gdal warp resample algorithm, falling back to bilinear'.format(sample_alg))
            self.sample_alg = 'bilinear'

        self.cache_dir = utils.cudem_cache() if cache_dir is None else cache_dir
        if utils.fn_url_p(self.fn):
            self.remote = True
            
        if self.valid_p():
            self.set_yield()
            self.inf(check_hash=True if self.data_format == -1 else False)
            self.set_transform()
            
    def __str__(self):
        return('<Dataset: {} - {}>'.format(self.metadata['name'], self.fn))
    
    def __repr__(self):
        return('<Dataset: {} - {}>'.format(self.metadata['name'], self.fn))

    def generate_inf(self, callback=lambda: False):
        """set in dataset"""
        
        raise(NotImplementedError)

    def yield_xyz(self):
        """set in dataset"""
        
        raise(NotImplementedError)

    def yield_array(self):
        """set in dataset"""
        
        raise(NotImplementedError)
    
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

    def valid_p(self, fmts=['<scratch-datalist>']):
        """check if self appears to be a valid dataset entry"""
        
        if self.fn is None:
            return(False)
        
        if self.data_format is None:
            return(False)
        
        if self.fn is not None:
            if self.fn not in fmts:
                if not utils.fn_url_p(self.fn):
                    if self.data_format != -11:
                        if not os.path.exists(self.fn):
                            return (False)
                        
                        if os.stat(self.fn).st_size == 0:
                            return(False)
                        
        return(True)
        
    def hash(self, sha1=False):
        """generate a hash of the xyz-dataset source file"""

        import hashlib
        BUF_SIZE = 65536
        if sha1:
            this_hash = hashlib.sha1()
        else:
            this_hash = hashlib.md5()
            
        try:
            with open(self.fn, 'rb') as f:
                while True:
                    data = f.read(BUF_SIZE)
                    if not data:
                        break
                    
                    this_hash.update(data)
            return(this_hash.hexdigest())
        except: return('0')

    def format_entry(self, sep=' '):

        dl_entry = sep.join([str(x) for x in [self.fn, self.data_format, self.weight]])
        metadata = self.echo_()
        return(sep.join([dl_entry, metadata]))
        
    def echo_(self, sep=' ', **kwargs):
        """print self as a datalist entry string"""
        
        return(sep.join([ '"{}"'.format(str(self.metadata[x])) for x in self.metadata.keys()]))
    
    def echo(self, **kwargs):
        """print self.data_entries as a datalist entries."""

        for entry in self.parse_json():
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
        The function `generate_inf()` should be defined for each specific
        dataset sub-class.
        """

        inf_path = '{}.inf'.format(self.fn)
        mb_inf = False
        self.infos = {}
        if os.path.exists(inf_path):
            try:
                with open(inf_path) as i_ob:
                    self.infos = json.load(i_ob)
            except ValueError:
                try:
                    self.infos = MBSParser(
                        fn=self.fn, src_srs=self.src_srs).inf_parse().infos
                    self.check_hash = False
                    mb_inf = True
                except:
                    if self.verbose:
                        utils.echo_error_msg(
                            'failed to parse inf {}'.format(inf_path)
                        )
            except:
                if self.verbose:
                    utils.echo_error_msg(
                        'failed to parse inf {}'.format(inf_path)
                    )
        
        if check_hash:
            if 'hash' in self.infos.keys():
                gen_inf = self.hash() != self.infos['hash']
            else:
                gen_inf = True
        elif not mb_inf:
            gen_inf = 'hash' not in self.infos.keys() or 'wkt' not in self.infos.keys()
        else:
            gen_inf = False

        if gen_inf:# and not self.remote:
            if self.verbose:
                _prog = utils.CliProgress('generating inf for {}'.format(self.fn))

            self.infos = self.generate_inf(None if not self.verbose else _prog.update)
            if self.infos is not None:
                if 'minmax' in self.infos:
                    if self.infos['minmax'] is not None:
                        if write_inf:
                            try:
                                with open('{}.inf'.format(self.fn), 'w') as inf:
                                    inf.write(json.dumps(self.infos))
                            except:
                                pass

                if recursive_check and self.parent is not None:
                    self.parent.inf(check_hash=True)

                if self.verbose:
                    _prog.end(0, 'generated inf for {}'.format(self.fn))
            else:
                sys.exit(-1)    
                
        if 'src_srs' not in self.infos.keys() or self.infos['src_srs'] is None:
            self.infos['src_srs'] = self.src_srs
        else:
            if self.src_srs is None:
                self.src_srs = self.infos['src_srs']
                
            if self.dst_trans is not None:
                if self.trans_region is None:
                    self.trans_region = regions.Region().from_list(self.infos['minmax'])
                    self.trans_region.src_srs = self.infos['src_srs']
                    self.trans_region.warp(self.dst_srs)
                    
        if 'format' not in self.infos.keys():
            self.infos['format'] = self.data_format

        return(self.infos)

    def set_yield(self):
        """set the yield strategy, either default block or mask"""

        if self.x_inc is not None and (self.data_format not in [200, 201, -11]):
            self.x_inc = utils.str2inc(self.x_inc)
            if self.y_inc is None:
                self.y_inc = self.x_inc
            else:
                self.y_inc = utils.str2inc(self.y_inc)
                
            self.xyz_yield = self.block_xyz()
        else:
            self.xyz_yield = self.yield_xyz()
        
    def set_transform(self):
        """Set the transformation parameters for the dataset."""
        
        if self.src_srs == '': self.src_srs = None
        if self.dst_srs == '': self.dst_srs = None
        if self.dst_srs is not None \
           and self.src_srs is not None \
           and self.src_srs != self.dst_srs:
            src_horz_epsg, src_vert_epsg = utils.epsg_from_input(self.src_srs)
            dst_horz_epsg, dst_vert_epsg = utils.epsg_from_input(self.dst_srs)
            src_srs = osr.SpatialReference()
            src_srs.SetFromUserInput('epsg:{}'.format(src_horz_epsg))
            dst_srs = osr.SpatialReference()
            dst_srs.SetFromUserInput('epsg:{}'.format(dst_horz_epsg))
            if dst_vert_epsg is not None \
               and src_vert_epsg is not None \
               and dst_vert_epsg != src_vert_epsg:
                vd_region = regions.Region(
                    src_srs=src_srs.ExportToProj4()
                ).from_list(
                    self.infos['minmax']
                ).warp(
                    dst_srs.ExportToProj4()
                ) if self.region is None else self.region.copy()
                
                vd_region.buffer(pct=2)
                trans_fn = os.path.join(
                    self.cache_dir,
                    '_vdatum_trans_{}_{}_{}'.format(
                        src_vert_epsg,
                        dst_vert_epsg,
                        vd_region.format('fn')
                    )
                )
                waffles_cmd = 'waffles -R {} -E 3s -M vdatum:vdatum_in={}:vdatum_out={} -O {} -c -k -D {}'.format(
                    vd_region.format('str'),
                    src_vert_epsg,
                    dst_vert_epsg,
                    trans_fn,
                    self.cache_dir
                )
                if utils.run_cmd(waffles_cmd, verbose=True)[1] == 0:
                    out_src_srs = '{} +geoidgrids={}.tif'.format(src_srs.ExportToProj4(), trans_fn)
                    out_dst_srs = '{}'.format(dst_srs.ExportToProj4())

                    if src_vert_epsg == '6360':
                        out_src_srs = out_src_srs + ' +vto_meter=0.3048006096012192'

                    if dst_vert_epsg == '6360':
                        out_dst_srs = out_dst_srs + ' +vto_meter=0.3048006096012192'
                else:
                    utils.echo_error_msg(
                        'failed to generate vertical transformation grid between {} and {} for this region!'.format(
                            src_vert_epsg, dst_vert_epsg
                        )
                    )
                    out_src_srs = src_srs.ExportToProj4()
                    out_dst_srs = dst_srs.ExportToProj4()
                
            else:
                out_src_srs = src_srs.ExportToProj4()
                out_dst_srs = dst_srs.ExportToProj4()

            src_osr_srs = osr.SpatialReference()
            src_osr_srs.SetFromUserInput(out_src_srs)
            dst_osr_srs = osr.SpatialReference()
            dst_osr_srs.SetFromUserInput(out_dst_srs)
            try:
                src_osr_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
                dst_osr_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
            except:
                pass

            self.src_trans_srs = out_src_srs
            self.dst_trans_srs = out_dst_srs
            self.dst_trans = osr.CoordinateTransformation(src_osr_srs, dst_osr_srs)
            if self.region is not None and self.region.src_srs != self.src_srs:
                self.trans_region = self.region.copy()
                self.trans_region.src_srs = out_dst_srs
                self.trans_region.warp(out_src_srs)
                self.trans_region.src_srs = out_src_srs
            else:
                self.trans_region = regions.Region().from_string(self.infos['wkt'])
                self.trans_region.src_srs = self.src_srs
                self.trans_region.warp(self.dst_srs)

            src_osr_srs = dst_osr_srs = None
            
    def parse_json(self):
        for ds in self.parse():
            yield(ds)
    
    def parse(self):
        """parse the datasets from the dataset.
        
        Re-define this method when defining a dataset sub-class
        that represents recursive data structures (datalists, zip, etc).

        This will fill self.data_entries
        """
        
        if self.region is not None:
            try:
                inf_region = regions.Region().from_string(self.infos['wkt'])
            except:
                try:
                    inf_region = regions.Region().from_list(self.infos['minmax'])
                except:
                    inf_region = self.region.copy()
                
            #if regions.regions_intersect_p(inf_region, self.region):
            if regions.regions_intersect_p(inf_region, self.region if self.dst_trans is None else self.trans_region):
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
        
        will convert all data to XYZ within the given region.
        """
        
        def xdl2dir(xdl):
            this_dir = []
            while True:
                if xdl.parent is None:
                    break
                this_dir.append(xdl.parent.name)
                xdl = xdl.parent
                
            this_dir.reverse()
            return(this_dir)
        
        if self.region is None:
            a_name = '{}_{}'.format(self.metadata['name'], utils.this_year())
        else:
            a_name = '{}_{}_{}'.format(
                self.metadata['name'], self.region.format('fn'), utils.this_year())
            
        self.parse_data_lists()
        with open('{}.datalist'.format(a_name), 'w') as dlf:
            for x in self.data_lists.keys():
                if self.region is None:
                    a_dir = '{}_{}'.format(x, utils.this_year())
                else:
                    a_dir = '{}_{}_{}'.format(x, self.region.format('fn'), utils.this_year())
                    
                this_dir = xdl2dir(self.data_lists[x]['parent'])
                this_dir.append(a_dir)
                tmp_dir = this_dir
                dlf.write(
                    '{}.datalist -1 {} {}\n'.format(
                        os.path.join(*(this_dir + [this_dir[-1]])),
                        1 if self.weight is None else self.weight, self.data_lists[x]['parent'].format_metadata()
                    )
                )
                this_dir = os.path.join(os.getcwd(), *this_dir)
                if not os.path.exists(this_dir):
                    os.makedirs(this_dir)
                    
                with open(
                        os.path.join(
                            this_dir, '{}.datalist'.format(os.path.basename(this_dir))), 'w'
                ) as sub_dlf:
                    for xyz_dataset in self.data_lists[x]['data']:
                        if xyz_dataset.remote == True:
                            ## we will want to split remote datasets into individual files if needed...
                            sub_xyz_path = '{}_{}.xyz'.format(xyz_dataset.fn.split(':')[0], self.region.format('fn'))
                        elif len(xyz_dataset.fn.split('.')) > 1:
                            xyz_ext = xyz_dataset.fn.split('.')[-1]
                            sub_xyz_path = '.'.join(
                                [utils.fn_basename(
                                    os.path.basename(
                                        utils.slugify(xyz_dataset.fn)
                                    ),
                                    xyz_dataset.fn.split('.')[-1]),
                                 'xyz']
                            )
                        else:
                            sub_xyz_path = '.'.join([xyz_dataset.fn, 'xyz'])

                        this_xyz_path = os.path.join(this_dir, sub_xyz_path)
                        sub_dlf.write('{} 168\n'.format(sub_xyz_path))
                        
                        with open(this_xyz_path, 'w') as xp:
                            for this_xyz in xyz_dataset.xyz_yield:
                                yield(this_xyz)
                                this_xyz.dump(
                                    include_w=True if self.weight is not None else False,
                                    dst_port=xp,
                                    encode=False
                                )

    def yield_block_array(self):
        """yield the xyz data as arrays

        Yields:
          list: an array for [weighted]mean, weights mask, and count
        """

        out_arrays = {'z':None, 'count':None, 'weight':None, 'mask':None}
        xcount, ycount, dst_gt = self.region.geo_transform(
            x_inc=self.x_inc, y_inc=self.y_inc, node='grid'
        )

        for this_xyz in self.yield_xyz():
            xpos, ypos = utils._geo2pixel(
                this_xyz.x, this_xyz.y, dst_gt, 'pixel'
            )
            if xpos < xcount and ypos < ycount and xpos >= 0 and ypos >= 0:            
                if self.weight is not None:
                    w = this_xyz.w
                else:
                    w = 1

                out_arrays['z'] = np.array([[this_xyz.z]])
                out_arrays['count'] = np.array([[1]])
                out_arrays['weight'] = np.array([[this_xyz.w]])
                this_srcwin = (xpos, ypos, 1, 1)
                yield(out_arrays, this_srcwin, dst_gt)
                                            
    def block_xyz(self, min_count=None, ndv=-9999, fmt='GTiff'):
        """block the src_xyz data to the mean block value

        Yields:
          list: xyz data for each block with data
        """

        ## initialize output rasters
        out_name = '{}_block'.format(utils.fn_basename2(os.path.basename(self.fn)))
        xcount, ycount, dst_gt = self.region.geo_transform(
            x_inc=self.x_inc, y_inc=self.y_inc, node='grid'
        )
        gdt = gdal.GDT_Float32
        c_gdt = gdal.GDT_Int32
        driver = gdal.GetDriverByName(fmt)

        z_ds = driver.Create('{}_s.tif'.format(out_name), xcount, ycount, 1, gdt,
                             options=['COMPRESS=LZW', 'PREDICTOR=2', 'TILED=YES'])
        z_ds.SetGeoTransform(dst_gt)
        z_band = z_ds.GetRasterBand(1)
        z_band.SetNoDataValue(ndv)
        
        w_ds = driver.Create('{}_w.tif'.format(out_name), xcount, ycount, 1, gdt,
                             options=['COMPRESS=LZW', 'PREDICTOR=2', 'TILED=YES'])
        w_ds.SetGeoTransform(dst_gt)
        w_band = w_ds.GetRasterBand(1)
        w_band.SetNoDataValue(ndv)
                
        c_ds = driver.Create('{}_c.tif'.format(out_name), xcount, ycount, 1, gdt)
        c_ds.SetGeoTransform(dst_gt)
        c_band = c_ds.GetRasterBand(1)
        c_band.SetNoDataValue(ndv)

        for arrs, srcwin, gt  in self.yield_block_array():
            arr = arrs['z']
            w_arr = arrs['weight']
            c_arr = arrs['count']
            c_arr[np.isnan(arr)] = 0
            w_arr[np.isnan(arr)] = 0
            arr[np.isnan(arr)] = 0
            z_array = z_band.ReadAsArray(srcwin[0], srcwin[1], srcwin[2], srcwin[3])
            w_array = w_band.ReadAsArray(srcwin[0], srcwin[1], srcwin[2], srcwin[3])
            c_array = c_band.ReadAsArray(srcwin[0], srcwin[1], srcwin[2], srcwin[3])
            z_array[z_array == ndv] = 0
            w_array[w_array == ndv] = 0
            c_array[c_array == ndv] = 0
            c_array += c_arr
            z_array += (arr * w_arr)
            w_array += w_arr
            c_array[c_array == 0] = ndv
            w_array[w_array == 0] = ndv
            z_array[np.isnan(w_array)] = ndv
            z_array[w_array == ndv] = ndv
            z_band.WriteArray(z_array, srcwin[0], srcwin[1])
            w_band.WriteArray(w_array, srcwin[0], srcwin[1])
            c_band.WriteArray(c_array, srcwin[0], srcwin[1])
            arr = w_arr = c_arr = z_array = w_array = c_array = None

        ## Finalize and close datasets
        srcwin = (0, 0, z_ds.RasterXSize, z_ds.RasterYSize)
        for y in range(
                srcwin[1], srcwin[1] + srcwin[3], 1
        ):
            z_data = z_band.ReadAsArray(srcwin[0], y, srcwin[2], 1)
            c_data = c_band.ReadAsArray(srcwin[0], y, srcwin[2], 1)
            w_data = w_band.ReadAsArray(srcwin[0], y, srcwin[2], 1)
            z_data[z_data == ndv] = np.nan
            w_data[w_data == ndv] = np.nan
            c_data[c_data == ndv] = np.nan
            z_data[np.isnan(w_data)] = np.nan
            w_data = w_data / c_data
            z_data = (z_data/w_data)/c_data
            
            for x in range(0, xcount):
                z = z_data[0, x]
                if not np.isnan(z):
                    geo_x, geo_y = utils._pixel2geo(x, y, gt)
                    out_xyz = xyzfun.XYZPoint(
                        x=geo_x, y=geo_y, z=z, w=w_data[0, x]
                    )
                    yield(out_xyz)
                    
            z_data[np.isnan(z_data)] = ndv
            c_data[np.isnan(c_data)] = ndv
            w_data[np.isnan(w_data)] = ndv
            z_band.WriteArray(z_data, srcwin[0], y)
            c_band.WriteArray(c_data, srcwin[0], y)
            w_band.WriteArray(w_data, srcwin[0], y)
                            
        z_ds = c_ds = w_ds = None
        
    def mask_xyz(self, dst_x_inc, dst_y_inc, dst_format='MEM', ndv=-9999, **kwargs):
        """Create a num grid mask of xyz data. The output grid
        will contain 1 where data exists and 0 where no data exists.
        returns the gdal dataset and config
        """

        xcount, ycount, dst_gt = self.region.geo_transform(x_inc=dst_x_inc, y_inc=dst_y_inc)
        ptArray = np.zeros((ycount, xcount))
        ds_config = demfun.set_infos(
            xcount,
            ycount,
            (xcount*ycount),
            dst_gt,
            utils.sr_wkt(self.dst_srs),
            gdal.GDT_Int32,
            ndv,
            'MEM'
        )
        for this_xyz in self.yield_xyz_from_entries(**kwargs):
            xpos, ypos = utils._geo2pixel(
                this_xyz.x, this_xyz.y, dst_gt, 'pixel'
            )
            try:
                ptArray[ypos, xpos] = 1
            except: pass

        driver = gdal.GetDriverByName(dst_format)
        ds = driver.Create('MEM', ds_config['nx'], ds_config['ny'], 1, ds_config['dt'])
        if ds is not None:
            ds.SetGeoTransform(ds_config['geoT'])
            ds.SetProjection(ds_config['proj'])
            ds.GetRasterBand(1).SetNoDataValue(ds_config['ndv'])
            ds.GetRasterBand(1).WriteArray(ptArray)
                
        return(ds, ds_config)
        
    def mask_and_yield_xyz(self, dst_gdal, dst_inc, dst_format='GTiff', ndv=-9999, **kwargs):
        """Create a num grid mask of xyz data. The output grid
        will contain 1 where data exists and 0 where no data exists.
        yields the xyz data
        """

        xcount, ycount, dst_gt = self.region.geo_transform(x_inc=dst_inc)
        ptArray = np.zeros((ycount, xcount))
        ds_config = demfun.set_infos(
            xcount,
            ycount,
            (xcount*ycount),
            dst_gt,
            utils.sr_wkt(self.src_srs),
            gdal.GDT_Float32,
            ndv,
            'GTiff'
        )

        for this_xyz in self.yield_xyz_from_entries(**kwargs):
            yield(this_xyz)
            xpos, ypos = utils._geo2pixel(
                this_xyz.x, this_xyz.y, dst_gt, 'pixel'
            )
            try:
                ptArray[ypos, xpos] = 1
            except:
                pass

        out, status = utils.gdal_write(ptArray, dst_gdal, ds_config)
                        
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
        fd = ogr.FieldDefn('weight', ogr.OFTReal)
        fd.SetWidth(6)
        fd.SetPrecision(6)
        layer.CreateField(fd)
            
        f = ogr.Feature(feature_def=layer.GetLayerDefn())        
        for this_xyz in self.xyz_yield:
            f.SetField(0, this_xyz.x)
            f.SetField(1, this_xyz.y)
            f.SetField(2, float(this_xyz.z))
            f.SetField(3, this_xyz.w)
                
            wkt = this_xyz.export_as_wkt(include_z=True)
            g = ogr.CreateGeometryFromWkt(wkt)
            f.SetGeometryDirectly(g)
            layer.CreateFeature(f)
            
        return(ogr_ds)
                    
    def dump_xyz(self, dst_port=sys.stdout, encode=False, **kwargs):
        """dump the XYZ data from the dataset"""
        
        for this_xyz in self.xyz_yield:
            this_xyz.dump(
                include_w=True if self.weight is not None else False,
                dst_port=dst_port,
                encode=encode
            )

    def dump_xyz_direct(self, dst_port=sys.stdout, encode=False, **kwargs):
        """dump the XYZ data from the dataset"""
        
        for this_xyz in self.yield_xyz():
            this_xyz.dump(
                include_w=True if self.weight is not None else False,
                dst_port=dst_port,
                encode=encode
            )

    def export_xyz_as_list(self, **kwargs):
        """return the XYZ data from the dataset as python list"""
        
        xyz_l = []
        for this_xyz in self.xyz_yield:
            xyz_l.append(this_xyz.copy())
        return(xyz_l)
    
## ==============================================
## ==============================================
class XYZFile(ElevationDataset):
    """representing an ASCII xyz dataset stream."""

    def __init__(
            self,
            delim=None,
            xpos=0,
            ypos=1,
            zpos=2,
            wpos=None,
            skip=0,
            x_scale=1,
            y_scale=1,
            z_scale=1,
            x_offset=0,
            y_offset=0,
            **kwargs
    ):
        self.ds_args = {
            'delim': delim,
            'xpos': xpos,
            'ypos': ypos,
            'zpos': zpos,
            'wpos': wpos,
            'skip': skip,
            'x_scale': x_scale,
            'y_scale': y_scale,
            'z_scale': z_scale,
            'x_offset': x_offset,
            'y_offset': y_offset,
        }
        
        self.xpos = utils.int_or(xpos, 0)
        self.ypos = utils.int_or(ypos, 1)
        self.zpos = utils.int_or(zpos, 2)
        self.wpos = utils.int_or(wpos)
        self.skip = utils.int_or(skip, 0)
        self.x_scale = utils.float_or(x_scale, 1)
        self.y_scale = utils.float_or(y_scale, 1)
        self.z_scale = utils.float_or(z_scale, 1)
        self.x_offset = utils.int_or(x_offset, 0)
        self.y_offset = utils.int_or(y_offset, 0)
        self._known_delims = [None, ',', '/', ':'] ## space and tab are 'None'
        self.delim = delim
        if delim is not None:
            self._known_delims.insert(0, delim)
            
        self.rem = False
        if x_offset == 'REM':
            x_offset = 0
            self.rem = True
        self.scoff = True if x_scale != 1 or y_scale != 1 or z_scale != 1 or x_offset != 0 or y_offset != 0 else False
            
        super().__init__(**kwargs)

        if self.src_srs is not None:
            self.set_transform()
        
    def generate_inf(self, callback=lambda: False):
        """generate a infos dictionary from the xyz dataset"""
                
        pts = []
        self.infos['name'] = self.fn
        self.infos['hash'] = self.hash()#dl_hash(self.fn)
        self.infos['numpts'] = 0
        self.infos['format'] = self.data_format
        this_region = regions.Region()
        region_ = self.region
        self.region = None

        for i, l in enumerate(self.yield_xyz()):
            if i == 0:
                this_region.from_list([l.x, l.x, l.y, l.y, l.z, l.z])
            else:
                if l.x < this_region.xmin:
                    this_region.xmin = l.x
                elif l.x > this_region.xmax:
                    this_region.xmax = l.x
                if l.y < this_region.ymin:
                    this_region.ymin = l.y
                elif l.y > this_region.ymax:
                    this_region.ymax = l.y
                if l.z < this_region.zmin:
                    this_region.zmin = l.z
                elif l.z > this_region.zmax:
                    this_region.zmax = l.z
            #pts.append(l.export_as_list(include_z = True))
            self.infos['numpts'] = i

        self.infos['minmax'] = this_region.export_as_list(include_z = True)
        if self.infos['numpts'] > 0:
            # try:
            #     out_hull = [pts[i] for i in ConvexHull(
            #         pts, qhull_options='Qt'
            #     ).vertices]
            #     out_hull.append(out_hull[0])
            #     self.infos['wkt'] = regions.create_wkt_polygon(out_hull, xpos=0, ypos=1)
            # except:
            self.infos['wkt'] = this_region.export_as_wkt()
                
        self.region = region_
        self.infos['src_srs'] = self.src_srs
        
        return(self.infos)

    def line_delim(self, xyz_line):
        """guess a line delimiter"""

        for delim in self._known_delims:
            try:
                this_xyz = xyz_line.split(delim)
                if len(this_xyz) > 1:
                    return(this_xyz)
            except:
                pass

    def yield_array(self):        
        for arrs in self.yield_block_array():
            yield(arrs)
            
    def yield_xyz(self):
        """xyz file parsing generator"""
        
        if self.fn is not None:
            if os.path.exists(str(self.fn)):
                self.src_data = open(self.fn, "r")
            else:
                self.src_data = self.fn
        else:
            self.src_data = sys.stdin

        count = 0
        skip = self.skip
        for xyz_line in self.src_data:
            #utils.echo_msg_inline('{} - {}'.format(self.fn, count))
            if count >= skip:
                this_xyz = self.line_delim(xyz_line)
                if this_xyz is None:
                    continue

                if self.wpos is not None:
                    w = float(this_xyz[self.wpos])
                else:
                    w = 1

                try:
                    this_xyz = xyzfun.XYZPoint(
                        x=this_xyz[self.xpos],
                        y=this_xyz[self.ypos],
                        z=this_xyz[self.zpos]
                    )
                except Exception as e:
                    utils.echo_error_msg('{} ; {}'.format(e, this_xyz))
                    this_xyz = xyzfun.XYZPoint()

                if this_xyz.valid_p():
                    if self.scoff:
                        this_xyz.x = (this_xyz.x+self.x_offset) * self.x_scale
                        this_xyz.y = (this_xyz.y+self.y_offset) * self.y_scale
                        this_xyz.z *= self.z_scale

                    if self.rem:
                        this_xyz.x = math.fmod(this_xyz.x+180,360)-180 

                    this_xyz.w = w if self.weight is None else self.weight * w                        
                    if self.dst_trans is not None:
                        this_xyz.transform(self.dst_trans)

                    if self.region is not None and self.region.valid_p():
                        if self.invert_region:
                            if not regions.xyz_in_region_p(this_xyz, self.region):
                                count += 1
                                yield(this_xyz)
                        else:
                            if regions.xyz_in_region_p(this_xyz, self.region):
                                count += 1
                                yield(this_xyz)
                    else:
                        count += 1
                        yield(this_xyz)

            else: skip -= 1

        if self.verbose:
            utils.echo_msg(
                'parsed {} data records from {}{}'.format(
                    count, self.fn, ' @{}'.format(self.weight) if self.weight is not None else ''
                )
            )
            
        self.src_data.close()

## ==============================================
## ==============================================
class LASFile(ElevationDataset):
    """representing an LAS/LAZ dataset."""

    def __init__(self, classes='0/2/29/40', **kwargs):
        self.ds_args = {
            'classes': classes
        }
        self.classes = [int(x) for x in classes.split('/')]
        super().__init__(**kwargs)

        if self.src_srs is None:
            self.src_srs = self.get_epsg()
            if self.src_srs is None:
                self.src_srs = self.infos['src_srs']
            
            self.set_transform()

    def get_epsg(self):
        with lp.open(self.fn) as lasf:
            lasf_vlrs = lasf.header.vlrs
            for vlr in lasf_vlrs:
                #if 'OGC WKT' in vlr.description:
                #utils.echo_msg(vlr.description)
                if vlr.record_id == 2112:
                    src_srs = osr.SpatialReference()
                    src_srs.SetFromUserInput(vlr.string)
                    src_srs.AutoIdentifyEPSG()
                    srs_auth = src_srs.GetAttrValue('AUTHORITY', 1)
                    vert_cs = src_srs.GetAttrValue('vert_cs')
                    proj4 = src_srs.ExportToProj4()
                    vert_auth = src_srs.GetAttrValue('VERT_CS|AUTHORITY',1)
                    src_srs = None
                    
                    if srs_auth is not None:
                        if vert_auth is not None:
                            return('epsg:{}+{}'.format(srs_auth, vert_auth))
                        else:
                            return('epsg:{}'.format(srs_auth))
                    else:
                        return(proj4)
                    break
            return(None)
        
    def generate_inf(self, callback=lambda: False):
        """generate an inf file for a lidar dataset."""
        
        self.infos['name'] = self.fn
        self.infos['hash'] = self.hash()#dl_hash(self.fn)
        self.infos['numpts'] = 0
        self.infos['format'] = self.data_format
        this_region = regions.Region()

        with lp.open(self.fn) as lasf:
            self.infos['numpts'] = lasf.header.point_count
            this_region.from_list(
                [lasf.header.x_min,
                 lasf.header.x_max,
                 lasf.header.y_min,
                 lasf.header.y_max,
                 lasf.header.z_min,
                 lasf.header.z_max]
            )

        self.infos['src_srs'] = self.src_srs if self.src_srs is not None else self.get_epsg()
        self.infos['minmax'] = this_region.export_as_list(
            include_z=True
        )
        self.infos['wkt'] = this_region.export_as_wkt()
        return(self.infos)

    def generate_inf_scan(self, callback=lambda: False):
        """generate an inf file for a lidar dataset.
        ... parse the data to obtain the hull region
        """

        self.infos['name'] = self.fn
        self.infos['hash'] = self.hash()#dl_hash(self.fn)
        self.infos['numpts'] = 0
        this_region = regions.Region()

        pts = []
        region_ = self.region
        self.region = None

        for i, l in enumerate(self.yield_xyz()):
            if i == 0:
                this_region.from_list([l.x, l.x, l.y, l.y, l.z, l.z])
            else:
                if l.x < this_region.xmin:
                    this_region.xmin = l.x
                elif l.x > this_region.xmax:
                    this_region.xmax = l.x
                if l.y < this_region.ymin:
                    this_region.ymin = l.y
                elif l.y > this_region.ymax:
                    this_region.ymax = l.y
                if l.z < this_region.zmin:
                    this_region.zmin = l.z
                elif l.z > this_region.zmax:
                    this_region.zmax = l.z
            pts.append(l.export_as_list(include_z = True))
            self.infos['numpts'] = i

        self.infos['minmax'] = this_region.export_as_list(include_z = True)
        if self.infos['numpts'] > 0:
            try:
                out_hull = [pts[i] for i in ConvexHull(
                    pts, qhull_options='Qt'
                ).vertices]
                out_hull.append(out_hull[0])
                self.infos['wkt'] = regions.create_wkt_polygon(out_hull, xpos=0, ypos=1)
            except:
                self.infos['wkt'] = this_region.export_as_wkt()
                
        self.region = region_
        return(self.infos)
        
    def yield_xyz(self):
        """LAS file parsing generator"""

        count = 0
        with lp.open(self.fn) as lasf:
            for points in lasf.chunk_iterator(2_000_000):
                points = points[(np.isin(points.classification, self.classes))]
                if self.region is not None  and self.region.valid_p():
                    tmp_region = self.region.copy() if self.dst_trans is None else self.trans_region.copy()
                    if self.invert_region:
                        points = points[((points.x > tmp_region.xmax) | (points.x < tmp_region.xmin)) | \
                                        ((points.y > tmp_region.ymax) | (points.y < tmp_region.ymin))]
                        if self.region.zmin is not None:
                            points =  points[(points.z < tmp_region.zmin)]
                            if self.region.zmax is not None:
                                points =  points[(points.z > tmp_region.zmax)]
                    else:
                        points = points[((points.x < tmp_region.xmax) & (points.x > tmp_region.xmin)) & \
                                        ((points.y < tmp_region.ymax) & (points.y > tmp_region.ymin))]
                        if self.region.zmin is not None:
                            points =  points[(points.z > tmp_region.zmin)]
                        if self.region.zmax is not None:
                            points =  points[(points.z < tmp_region.zmax)]
                                              
                dataset = np.vstack((points.x, points.y, points.z)).transpose()
                # if self.region is not None  and self.region.valid_p():
                #     tmp_region = self.region.copy() if self.dst_trans is None else self.trans_region.copy()
                #     if not self.invert_region:
                #         dataset = dataset[dataset[:,0] > tmp_region.xmin,:]
                #         dataset = dataset[dataset[:,0] < tmp_region.xmax,:]
                #         dataset = dataset[dataset[:,1] > tmp_region.ymin,:]
                #         dataset = dataset[dataset[:,1] < tmp_region.ymax,:]
                #         if self.region.zmin is not None:
                #             dataset = dataset[dataset[:,2] > tmp_region.zmin,:]

                #         if self.region.zmax is not None:
                #             dataset = dataset[dataset[:,2] < tmp_region.zmax,:]
                
                count += len(dataset)
                for point in dataset:
                    this_xyz = xyzfun.XYZPoint(
                        x=point[0],
                        y=point[1],
                        z=point[2],
                        w=self.weight
                    )
                    if self.dst_trans is not None:
                        this_xyz.transform(self.dst_trans)
                    
                    yield(this_xyz)
                dataset = None

        if self.verbose:
            utils.echo_msg(
                'parsed {} data records from {}{}'.format(
                    count, self.fn, ' @{}'.format(self.weight) if self.weight is not None else ''
                )
            )

    def yield_array(self):        
        for arrs in self.yield_block_array():
            yield(arrs)
            
    # def yield_array(self):        
    #     for src_arrs, src_gt, src_ds_config in self.yield_block_array():
    #         x_count, y_count, dst_gt = self.region.geo_transform(self.x_inc, self.y_inc)
    #         srcwin_region = regions.Region().from_geo_transform(
    #             geo_transform=src_gt, x_count=src_ds_config['nx'], y_count=src_ds_config['ny']
    #         )
    #         dst_srcwin = srcwin_region.srcwin(dst_gt, x_count, y_count, 'grid')
    #         yield(src_arrs, dst_srcwin, src_gt)
            
## ==============================================
## ==============================================
class RasterFile(ElevationDataset):
    """providing a GDAL raster dataset parser."""
    
    def __init__(
            self,
            mask=None,
            weight_mask=None,
            open_options=None,
            sample=None,
            resample=True,
            check_path=True,
            super_grid=False,
            **kwargs
    ):

        self.ds_args = {
            'mask': mask,
            'weight_mask': weight_mask,
            'open_options': open_options,
            'sample': sample,
            'resample': resample,
            'check_path': check_path,
            'super_grid': super_grid
        }
        
        try:
            self.open_options = open_options.split('/')
        except AttributeError:
            self.open_options = open_options
        except:
            self.open_options = None
        
        super().__init__(**kwargs)

        self.mask = mask
        self.weight_mask = weight_mask

        if self.valid_p() and self.src_srs is None:
            self.src_srs = demfun.get_srs(self.fn)
            
        self.set_transform()
        self.sample_alg = sample if sample is not None else self.sample_alg

        self.resample = resample
        self.check_path = check_path
        self.dem_infos = demfun.infos(self.fn)
        self.super_grid = super_grid

        if self.x_inc is not None and self.y_inc is not None and self.resample:
            self.resample_and_warp = True
        else:
            self.resample_and_warp = False

    def init_ds(self):
        """initialize the raster dataset

        if x/y incs are set, will warp raster to that resolution.
        """

        if self.check_path and not os.path.exists(self.fn):
            return(None)

        ndv = utils.float_or(demfun.get_nodata(self.fn), -9999)

        if self.region is not None:
            self.warp_region = self.region.copy()
        else:
            self.warp_region = regions.Region().from_list(self.infos['minmax'])
            if self.dst_trans is not None:
                self.warp_region.src_srs = self.src_srs
                self.warp_region.warp(self.dst_srs)

        # if self.region.src_srs != self.src_srs:
        #     #if not regions.regions_within_ogr_p(self.warp_region, self.region) or self.invert_region:
        #     #    self.warp_region = self.region.copy()
        #     #else:
        #     if regions.regions_within_ogr_p(self.warp_region, self.region):
        #         self.warp_region.cut(self.region, self.x_inc, self.y_inc)
                
        #     if self.dst_trans is not None:
        #         self.dst_trans = None
                            
        if self.resample_and_warp:
            if self.dst_trans is not None:
                self.dst_trans = None

            #if self.resample:
            dem_inf = demfun.infos(self.fn)

            tmp_ds = self.fn
            if self.open_options is not None:
                src_ds = gdal.OpenEx(self.fn, open_options=self.open_options)
            else:
                src_ds = gdal.Open(self.fn)

            if src_ds is not None:
                remake = False
                mt = src_ds.GetMetadata()            
                ## remake this grid if it's grid-node
                if 'AREA_OR_POINT' in mt.keys():
                    if mt['AREA_OR_POINT'].lower() == 'point':
                        remake = True
                        
                if self.open_options is not None:# or self.super_grid:
                    remake = True

                if remake:
                    ds_config = demfun.gather_infos(src_ds)
                    tmp_ds = demfun.generate_mem_ds(ds_config)
                    band = tmp_ds.GetRasterBand(1)
                    band.WriteArray(src_ds.GetRasterBand(1).ReadAsArray())
                    tmp_ds.SetProjection(utils.sr_wkt(self.src_srs))
                    tmp_ds.FlushCache()
                # else:
                #     tmp_ds = gdal.Open(self.fn)
                    
                src_ds = None

            ## Sample/Warp
            tmp_warp = '_tmp_gdal.tif'
            warp_ = demfun.sample_warp(
                tmp_ds, tmp_warp, self.x_inc, self.y_inc,
                src_srs=self.src_trans_srs, dst_srs=self.dst_trans_srs,
                src_region=self.warp_region, sample_alg=self.sample_alg,
                ndv=ndv, verbose=False
            )[0] 
            tmp_ds = None
            warp_ds = gdal.Open(tmp_warp)
            #utils.remove_glob('_tmp_gdal.tif')

            if warp_ds is not None:
                ## clip
                warp_ds_config = demfun.gather_infos(warp_ds)
                gt = warp_ds_config['geoT']
                srcwin = self.warp_region.srcwin(gt, warp_ds.RasterXSize, warp_ds.RasterYSize, node='grid')
                warp_arr = warp_ds.GetRasterBand(1).ReadAsArray(srcwin[0], srcwin[1], srcwin[2], srcwin[3])
                dst_gt = (gt[0] + (srcwin[0] * gt[1]), gt[1], 0., gt[3] + (srcwin[1] * gt[5]), 0., gt[5])
                out_ds_config = demfun.set_infos(srcwin[2], srcwin[3], srcwin[2] * srcwin[3], dst_gt,
                                                 warp_ds_config['proj'], warp_ds_config['dt'], warp_ds_config['ndv'],
                                                 warp_ds_config['fmt'])

                src_ds = demfun.generate_mem_ds(out_ds_config)
                if src_ds is not None:
                    band = src_ds.GetRasterBand(1)
                    band.WriteArray(warp_arr)
                    src_ds.FlushCache()
                    
                warp_ds = warp_arr = None
                utils.remove_glob(tmp_warp)
            
        else:

            if self.open_options:
                src_ds = gdal.OpenEx(self.fn, open_options=self.open_options)
                if src_ds is None:
                    if self.verbose:
                        utils.echo_warning_msg('could not open file using open options {}'.format(self.open_options))
                        
                    src_ds = gdal.Open(self.fn)
            else:
                src_ds = gdal.Open(self.fn)

        if self.invert_region:
            src_ds_config = demfun.gather_infos(src_ds)
            srcwin = self.warp_region.srcwin(src_ds_config['geoT'], src_ds_config['nx'], src_ds_config['ny'], node='grid')
            
            driver = gdal.GetDriverByName('MEM')
            mem_ds = driver.Create('msk', src_ds_config['nx'], src_ds_config['ny'], 1, src_ds_config['dt'])
            mem_ds.SetGeoTransform(src_ds_config['geoT'])
            mem_ds.SetProjection(src_ds_config['proj'])
            mem_band = mem_ds.GetRasterBand(1)
            mem_band.SetNoDataValue(src_ds_config['ndv'])
            mem_arr = np.ones((src_ds_config['ny'], src_ds_config['nx']))
            mem_arr[srcwin[0]:srcwin[0]+srcwin[2],
                    srcwin[1]:srcwin[1]+srcwin[3]] = src_ds_config['ndv']
            mem_band.WriteArray(mem_arr)

            self.mask = mem_ds
            
        return(src_ds)
        
    def generate_inf(self, check_z=False, callback=lambda: False):
        """generate a infos dictionary from the raster dataset"""
            
        self.infos['name'] = self.fn
        self.infos['hash'] = self.hash()#dl_hash(self.fn)
        self.infos['src_srs'] = self.src_srs if self.src_srs is not None else demfun.get_srs(self.fn)
        self.infos['format'] = self.data_format
        #src_ds = gdal.Open(self.fn)

        if self.open_options is not None:
            src_ds = gdal.OpenEx(self.fn, open_options=self.open_options)
        else:
            src_ds = gdal.Open(self.fn)
        
        if src_ds is not None:
            gt = src_ds.GetGeoTransform()
            this_region = regions.Region(src_srs=self.src_srs).from_geo_transform(
                #geo_transform=self.dem_infos['geoT'],
                #x_count=self.dem_infos['nx'],
                #y_count=self.dem_infos['ny']
                geo_transform=gt,
                x_count=src_ds.RasterXSize,
                y_count=src_ds.RasterYSize
            )
            if check_z:
                try:
                    zr = src_ds.GetRasterBand(1).ComputeRasterMinMax()
                except:
                    zr = [None, None]
            else:
                zr = [None, None]
                
            this_region.zmin, this_region.zmax = zr[0], zr[1]
            self.infos['minmax'] = this_region.export_as_list(include_z=True)
            self.infos['numpts'] = src_ds.RasterXSize * src_ds.RasterYSize
            self.infos['wkt'] = this_region.export_as_wkt()
            src_ds = None
            
        return(self.infos)

    def get_srcwin(self, gt, x_size, y_size, node='grid'):
        if self.region is not None:
            if self.invert_region:
                srcwin = (
                    0, 0, x_size, y_size, node
                )
            else:
                if self.dst_trans is not None:
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
    
    def yield_array(self):
        """parse the data from gdal dataset src_ds (first band only)"""

        src_ds = self.init_ds()
        out_xyz = xyzfun.XYZPoint(w=1)
        if src_ds is None:
            if self.verbose:
                utils.echo_error_msg('could not load raster file {}'.format(self.fn))
        else:
            count = 0
            band = src_ds.GetRasterBand(1)
            gt = src_ds.GetGeoTransform()
            ndv = utils.float_or(band.GetNoDataValue())
            mask_band = None
            weight_band = None
            nodata = ['{:g}'.format(-9999), 'nan', float('nan')]
            if ndv is not None:
                nodata.append('{:g}'.format(ndv))

            band_data = count_data = weight_data = mask_data = None
            out_arrays = {'z':None, 'count':None, 'weight':None, 'mask':None}
            if self.weight_mask is not None:
                if self.x_inc is not None and self.y_inc is not None:
                    src_weight = demfun.sample_warp(
                        self.weight_mask, None, self.x_inc, self.y_inc,
                        src_region=self.warp_region, sample_alg=self.sample_alg,
                         ndv=ndv, verbose=False
                    )[0]
                else:
                    src_weight = gdal.Open(self.weight_mask)
                    
                weight_band = src_weight.GetRasterBand(1)

            if self.mask is not None:
                if self.x_inc is not None and self.y_inc is not None:
                    src_mask = demfun.sample_warp(
                        self.mask, None, self.x_inc, self.y_inc,
                        src_region=self.warp_region, sample_alg=self.sample_alg,
                        ndv=ndv, verbose=False
                    )[0]
                else:
                    if not self.invert_region:
                        src_mask = gdal.Open(self.mask)
                    else:
                        src_mask = self.mask

                mask_band = src_mask.GetRasterBand(1)

            srcwin = self.get_srcwin(gt, src_ds.RasterXSize, src_ds.RasterYSize, node='pixel')
            for y in range(
                    srcwin[1], srcwin[1] + srcwin[3], 1
            ):
                band_data = band.ReadAsArray(
                    srcwin[0], y, srcwin[2], 1
                )
                if ndv is not None and not np.isnan(ndv):
                    band_data[band_data == ndv] = np.nan

                if np.all(np.isnan(band_data)):
                    continue
                
                this_origin = utils._pixel2geo(srcwin[0], y, gt, node='grid')
                this_gt = (this_origin[0], gt[1], 0, this_origin[1], 0, gt[5])
                                    
                if weight_band is not None:
                    weight_data = weight_band.ReadAsArray(
                        srcwin[0], y, srcwin[2], 1
                    )
                    weight_ndv = float(weight_band.GetNoDataValue())
                    if not np.isnan(weight_ndv):
                        weight_data[weight_data==weight_ndv] = np.nan
                else:
                    weight_data = np.ones(band_data.shape)
                    if self.weight:
                        weight_data[:] = self.weight
                                       
                if mask_band is not None:
                    mask_data = mask_band.ReadAsArray(
                        srcwin[0], y, srcwin[2], 1
                    )
                    mask_ndv = float(mask_band.GetNoDataValue())
                    if not np.isnan(mask_ndv):
                        mask_data[mask_data==mask_ndv] = np.nan

                    band_data[np.isnan(mask_data)] = np.nan
                    if weight_band is not None:
                        weight_data[np.isnan(mask_data)] = np.nan
                    
                if self.region is not None and self.region.valid_p():
                    z_region = self.region.z_region()

                    if self.invert_region:
                        if z_region[0] is not None and z_region[1] is not None:
                            band_data[(band_data > z_region[0]) & (band_data < z_region[1])] = np.nan
                        
                        elif z_region[0] is not None:
                            band_data[band_data > z_region[0]] = np.nan
                        
                        elif z_region[1] is not None:
                            band_data[band_data < z_region[1]] = np.nan
                    else:
                        if z_region[0] is not None:
                            band_data[band_data < z_region[0]] = np.nan
                        
                        if z_region[1] is not None:
                            band_data[band_data > z_region[1]] = np.nan

                    if weight_band is not None:
                        w_region = self.region.w_region()
                        if self.invert_region:
                            if w_region[0] is not None and w_region[1] is not None:
                                band_data[(band_data > w_region[0]) & (band_data < w_region[1])] = np.nan

                            elif w_region[0] is not None:
                                band_data[band_data > w_region[0]] = np.nan
                        
                            elif w_region[1] is not None:
                                band_data[band_data < w_region[1]] = np.nan
                        else:
                            if w_region[0] is not None:
                                band_data[weight_data < w_region[0]] = np.nan
                        
                            if w_region[1] is not None:
                                band_data[weight_data > w_region[1]] = np.nan

                count_data = np.zeros(band_data.shape)
                count_data[~np.isnan(band_data)] = 1
                #count += np.count_nonzero(~np.isnan(count_data))
                count += np.count_nonzero(count_data)
                out_arrays = {'z':band_data, 'count':count_data, 'weight':weight_data, 'mask':mask_data}
                ds_config = demfun.set_infos(
                    1,
                    y,
                    y,
                    this_gt,
                    self.dst_srs,
                    gdal.GDT_Float32,
                    np.nan,
                    'GTiff'
                )

                this_srcwin = (srcwin[0], y, srcwin[2], 1)
                yield(out_arrays, this_srcwin, this_gt)
                                            
            band = mask_band = weight_band = src_weight = src_mask = src_ds = None
            if self.verbose:
                utils.echo_msg(
                    'parsed {} data records from {}{}'.format(
                        count, self.fn, ' @{}'.format(self.weight) if self.weight is not None else ''
                    )
                )
                
    def yield_xyz(self):
        for arrs, srcwin, gt in self.yield_array():
            z_array = arrs['z']
            w_array = arrs['weight']
            ycount, xcount = z_array.shape
            for y in range(0, ycount):
                for x in range(0, xcount):
                    z = z_array[y,x]
                    if not np.isnan(z):
                        geo_x, geo_y = utils._pixel2geo(x, y, gt)
                        out_xyz = xyzfun.XYZPoint(
                            x=geo_x, y=geo_y, z=z, w=w_array[y,x]
                        )
                        #print(self.dst_trans)
                        if self.dst_trans is not None and not self.resample_and_warp:
                            out_xyz.transform(self.dst_trans)

                        #print(out_xyz)
                        yield(out_xyz)
            
## ==============================================
## ==============================================
class BAGFile(ElevationDataset):
    """providing a BAG raster dataset parser.

    process supergrids at native resolution if they
    exist, otherwise process as normal grid
    """

    def __init__(self, explode=False, force_vr=False, vr_strategy='MIN', **kwargs):
        self.ds_args = {
            'explode': explode,
            'force_vr': force_vr,
            'vr_strategy': vr_strategy,
        }
        super().__init__(**kwargs)
        self.explode = explode
        self.force_vr = force_vr
        self.vr_strategy = vr_strategy
        
        if self.src_srs is None:
            self.src_srs = demfun.get_srs(self.fn)
            self.set_transform()
        
    def generate_inf(self, callback=lambda: False):
        """generate a infos dictionary from the raster dataset"""
            
        self.infos['name'] = self.fn
        self.infos['hash'] = self.hash()#dl_hash(self.fn)
        self.infos['src_srs'] = self.src_srs if self.src_srs is not None else demfun.get_srs(self.fn)
        self.infos['format'] = self.data_format
        src_ds = gdal.Open(self.fn)
        if src_ds is not None:
            gt = src_ds.GetGeoTransform()
            this_region = regions.Region(src_srs=self.src_srs).from_geo_transform(
                geo_transform=gt,
                x_count=src_ds.RasterXSize,
                y_count=src_ds.RasterYSize
            )
            try:
                zr = src_ds.GetRasterBand(1).ComputeRasterMinMax()
            except:
                zr = [None, None]
                
            this_region.zmin, this_region.zmax = zr[0], zr[1]
            self.infos['minmax'] = this_region.export_as_list(include_z=True)
            self.infos['numpts'] = src_ds.RasterXSize * src_ds.RasterYSize
            self.infos['wkt'] = this_region.export_as_wkt()
            src_ds = None
            
        return(self.infos)

    def parse_(self, resample=True):
        mt = gdal.Info(self.fn, format='json')['metadata']['']
        oo = []

        if self.region is not None and self.region.valid_p():
            bag_region = self.region.copy() if self.dst_trans is None else self.trans_region.copy()
            inf_region = regions.Region().from_list(self.infos['minmax'])
            bag_region = regions.regions_reduce(bag_region, inf_region)
            bag_region.src_srs = self.infos['src_srs']
            #print(bag_region)
            
            oo.append('MINX={}'.format(bag_region.xmin))
            oo.append('MAXX={}'.format(bag_region.xmax))
            oo.append('MINY={}'.format(bag_region.ymin))
            oo.append('MAXY={}'.format(bag_region.ymax))
        else:
            bag_region = regions.Region().from_list(self.infos['minmax'])
            bag_region.src_srs = self.infos['src_srs']

        if self.dst_trans is not None:
            bag_region.warp(self.dst_srs)
            
        if ('HAS_SUPERGRIDS' in mt.keys() and mt['HAS_SUPERGRIDS'] == 'TRUE') \
           or self.force_vr \
           or 'MAX_RESOLUTION_X' in mt.keys() \
           or 'MAX_RESOLUTION_Y' in mt.keys():
            if self.explode:
                oo.append("MODE=LIST_SUPERGRIDS")
                src_ds = gdal.OpenEx(self.fn, open_options=oo)
                sub_datasets = src_ds.GetSubDatasets()
                src_ds = None

                for sub_dataset in sub_datasets:
                    sub_ds = RasterFile(
                        fn=sub_dataset[0],
                        data_format=200,
                        src_srs=self.src_srs,
                        dst_srs=self.dst_srs,
                        weight=self.weight,
                        src_region=bag_region,
                        x_inc=self.x_inc,
                        y_inc=self.y_inc,
                        verbose=self.verbose,
                        resample=resample,
                        check_path=False,
                        super_grid=True
                    )

                    sub_ds.infos = {}
                    sub_ds.generate_inf()
                    yield(sub_ds)

            else:
                oo.append("MODE=RESAMPLED_GRID")
                oo.append("RES_STRATEGY={}".format(self.vr_strategy))

                sub_ds = RasterFile(
                    fn=self.fn,
                    data_format=200,
                    open_options=oo,
                    src_srs=self.src_srs,
                    dst_srs=self.dst_srs,
                    weight=self.weight,
                    src_region=bag_region,
                    x_inc=self.x_inc,
                    y_inc=self.y_inc,
                    verbose=self.verbose,
                    resample=resample
                )

                yield(sub_ds)
        else:
            sub_ds = RasterFile(
                fn=self.fn,
                data_format=200,
                src_srs=self.src_srs,
                dst_srs=self.dst_srs,
                weight=self.weight,
                src_region=self.region,
                x_inc=self.x_inc,
                y_inc=self.y_inc,
                verbose=self.verbose,
                resample=resample
            )
            yield(sub_ds)

    def yield_xyz(self):
        for ds in self.parse_():
            for xyz in ds.yield_xyz():
                yield(xyz)

    def yield_array(self):
        for ds in self.parse_():
            for arr in ds.yield_array():
                yield(arr)
                
## ==============================================
## ==============================================
class MBSParser(ElevationDataset):
    """providing an mbsystem parser"""

    def __init__(self, mb_fmt=None, mb_exclude='A', **kwargs):
        self.ds_args = {
            'mb_fmt': mb_fmt,
            'mb_exclude': mb_exclude,
        }
        super().__init__(**kwargs)
        self.mb_fmt = mb_fmt
        self.mb_exclude = mb_exclude
        self.infos = {}
        
    def generate_inf(self, callback=lambda: False):
        self.infos['name'] = self.fn
        self.infos['hash'] = None
        self.infos['format'] = self.data_format
        try:
            utils.run_cmd('mbdatalist -O -V -I{}'.format(self.fn))
            self.inf_parse()
        except: pass
            
        return(self.infos)
            
    def inf_parse(self):
        self.infos['name'] = self.fn
        self.infos['minmax'] = [0,0,0,0,0,0]
        self.infos['hash'] = None
        this_row = 0
        dims = []
        with open('{}.inf'.format(self.fn)) as iob:
            for il in iob:
                til = il.split()
                if len(til) > 1:
                    if til[0] == 'Swath':
                        if til[2] == 'File:':
                            self.infos['name'] = til[3]
                            
                    if til[0] == 'Number':
                        if til[2] == 'Records:':
                            self.infos['numpts'] = utils.int_or(til[3])
                            
                    if til[0] == 'Minimum':
                        if til[1] == 'Longitude:':
                            self.infos['minmax'][0] = utils.float_or(til[2])
                            self.infos['minmax'][1] = utils.float_or(til[5])
                        elif til[1] == 'Latitude:':
                            self.infos['minmax'][2] = utils.float_or(til[2])
                            self.infos['minmax'][3] = utils.float_or(til[5])
                        elif til[1] == 'Depth:':
                            self.infos['minmax'][4] = utils.float_or(til[5]) * -1
                            self.infos['minmax'][5] = utils.float_or(til[2]) * -1
                            
                    if til[0] == 'CM':
                        if til[1] == 'dimensions:':
                            dims = [utils.int_or(til[2]), utils.int_or(til[3])]
                            cm_array = np.zeros((dims[0], dims[1]))
                            
                    if til[0] == 'CM:':
                        for j in range(0, dims[0]):
                            cm_array[this_row][j] = utils.int_or(til[j+1])
                        this_row += 1

        mbs_region = regions.Region().from_list(self.infos['minmax'])
        xinc = (mbs_region.xmax - mbs_region.xmin) / dims[0]
        yinc = (mbs_region.ymin - mbs_region.ymax) / dims[1]
        if abs(xinc) > 0 and abs(yinc) > 0:
            xcount, ycount, dst_gt = mbs_region.geo_transform(
                x_inc=xinc, y_inc=yinc
            )
            ds_config = {'nx': dims[0], 'ny': dims[1], 'nb': dims[1] * dims[0],
                         'geoT': dst_gt, 'proj': utils.sr_wkt(self.src_srs),
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

        self.infos['wkt'] = wkt
        return(self)

    def parse_(self):
        if self.x_inc is not None and self.y_inc is not None:        
            with open('_mb_grid_tmp.datalist', 'w') as tmp_dl:
                tmp_dl.write('{} {} {}\n'.format(self.fn, self.mb_fmt if self.mb_fmt is not None else '', self.weight if self.mb_fmt is not None else ''))

            ofn = '_'.join(os.path.basename(self.fn).split('.')[:-1])
            mbgrid_region = self.region.copy()
            mbgrid_region = mbgrid_region.buffer(pct=2, x_inc=self.x_inc, y_inc=self.y_inc)
            utils.run_cmd(
                'mbgrid -I_mb_grid_tmp.datalist {} -E{}/{}/degrees! -O{} -A2 -F1 -C10/1 -S0 -T35'.format(
                    mbgrid_region.format('gmt'), self.x_inc, self.y_inc, ofn
                ), verbose=True
            )
            utils.gdal2gdal('{}.grd'.format(ofn))
            utils.remove_glob('_mb_grid_tmp.datalist', '{}.cmd'.format(ofn), '{}.mb-1'.format(ofn), '{}.grd*'.format(ofn))
            mbs_ds = RasterFile(
                fn='{}.tif'.format(ofn),
                data_format=200,
                src_srs=self.src_srs,
                dst_srs=self.dst_srs,
                weight=self.weight,
                x_inc=self.x_inc,
                y_inc=self.y_inc,
                sample_alg=self.sample_alg,
                src_region=self.region,
                verbose=self.verbose
            )

            yield(mbs_ds)
            utils.remove_glob('{}.tif*'.format(ofn))
        else:
            yield(None)
            
    def yield_array(self):
        for ds in self.parse_():
            if ds is not None:
                for arr in ds.yield_array():
                    yield(arr)
            else:
                utils.echo_error_msg('could not parse MBS data {}'.format(self.fn))

    def yield_xyz(self):
        for ds in self.parse_():
            if ds is None:
                for line in utils.yield_cmd(
                        'mblist -M{} -OXYZ -I{}'.format(self.mb_exclude, self.fn),
                        verbose=True,
                ):
                    this_xyz = xyzfun.XYZPoint().from_string(line, delim='\t')
                    this_xyz.weight = self.weight
                    yield(this_xyz)
            else:
                for xyz in ds.yield_xyz():
                    yield(xyz)
                
    def yield_array_DEP(self):
        if self.x_inc is not None and self.y_inc is not None:        
            with open('_mb_grid_tmp.datalist', 'w') as tmp_dl:
                tmp_dl.write('{} {} {}\n'.format(self.fn, self.mb_fmt if self.mb_fmt is not None else '', self.weight if self.mb_fmt is not None else ''))

            ofn = '_'.join(os.path.basename(self.fn).split('.')[:-1])

            mbgrid_region = self.region.copy()
            mbgrid_region = mbgrid_region.buffer(pct=2, x_inc=self.x_inc, y_inc=self.y_inc)

            utils.run_cmd(
                'mbgrid -I_mb_grid_tmp.datalist {} -E{}/{}/degrees! -O{} -A2 -F1 -C10/1 -S0 -T35'.format(
                    mbgrid_region.format('gmt'), self.x_inc, self.y_inc, ofn
                ), verbose=True
            )
            utils.gdal2gdal('{}.grd'.format(ofn))
            utils.remove_glob('_mb_grid_tmp.datalist', '{}.cmd'.format(ofn), '{}.mb-1'.format(ofn), '{}.grd*'.format(ofn))
            #demfun.set_nodata('{}.tif'.format(ofn), nodata=-99999, convert_array=True, verbose=False)
            xyz_ds = RasterFile(
                fn='{}.tif'.format(ofn),
                data_format=200,
                src_srs=self.src_srs,
                dst_srs=self.dst_srs,
                weight=self.weight,
                x_inc=self.x_inc,
                y_inc=self.y_inc,
                sample_alg=self.sample_alg,
                src_region=self.region,
                verbose=self.verbose
            )
            for arr in xyz_ds.yield_array():
                yield(arr)

            utils.remove_glob('{}.tif*'.format(ofn))
            
    def yield_xyz_DEP(self):
        if self.x_inc is not None and self.y_inc is not None:        
            with open('_mb_grid_tmp.datalist', 'w') as tmp_dl:
                tmp_dl.write('{} {} {}\n'.format(self.fn, self.mb_fmt if self.mb_fmt is not None else '', self.weight if self.mb_fmt is not None else ''))

            ofn = '_'.join(os.path.basename(self.fn).split('.')[:-1])

            mbgrid_region = self.region.copy()
            mbgrid_region = mbgrid_region.buffer(pct=2, x_inc=self.x_inc, y_inc=self.y_inc)

            utils.run_cmd(
                'mbgrid -I_mb_grid_tmp.datalist {} -E{}/{}/degrees! -O{} -A2 -F1 -C10/1 -S0 -T35'.format(
                    mbgrid_region.format('gmt'), self.x_inc, self.y_inc, ofn
                ), verbose=True
            )
            utils.gdal2gdal('{}.grd'.format(ofn))
            utils.remove_glob('_mb_grid_tmp.datalist', '{}.cmd'.format(ofn), '{}.mb-1'.format(ofn), '{}.grd*'.format(ofn))
            #demfun.set_nodata('{}.tif'.format(ofn), nodata=-99999, convert_array=True, verbose=False)
            xyz_ds = RasterFile(
                fn='{}.tif'.format(ofn),
                data_format=200,
                src_srs=self.src_srs,
                dst_srs=self.dst_srs,
                weight=self.weight,
                x_inc=self.x_inc,
                y_inc=self.y_inc,
                sample_alg=self.sample_alg,
                src_region=self.region,
                verbose=self.verbose
            )
            for xyz in xyz_ds.yield_xyz():
                yield(xyz)

            utils.remove_glob('{}.tif*'.format(ofn))
        else:
            for line in utils.yield_cmd(
                    'mblist -M{} -OXYZ -I{}'.format(self.mb_exclude, self.fn),
                    verbose=True,
            ):
                this_xyz = xyzfun.XYZPoint().from_string(line, delim='\t')
                this_xyz.weight = self.weight
                yield(this_xyz)
### End
