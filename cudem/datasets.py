### datasets.py - Datasets
##
## Copyright (c) 2010 - 2022 Regents of the University of Colorado
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

    def __init__(
            self,
            fn=None,
            data_format=None,
            weight=1,
            src_srs=None,
            dst_srs=None,
            x_inc=None,
            y_inc=None,
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
        self.region = src_region
        self.parent = parent
        self.verbose = verbose
        self.data_entries = []
        self.data_lists = {}        
        self.remote = remote
        self.metadata = copy.deepcopy(metadata)
        self.x_inc = utils.str2inc(x_inc)
        self.y_inc = utils.str2inc(y_inc)
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

        if gen_inf:
            if self.verbose:
                _prog = utils.CliProgress('generating inf for {}'.format(self.fn))

            self.infos = self.generate_inf(None if not self.verbose else _prog.update)
            if 'minmax' in self.infos:
                if self.infos['minmax'] is not None:
                    if write_inf:
                        try:
                            with open('{}.inf'.format(self.fn), 'w') as inf:
                                inf.write(json.dumps(self.infos))
                        except:
                            pass
                        # if self.region is not None:
                        #     with open('{}_{}.inf'.format(
                        #             'dlim_tmp', self.region.format('fn')), 'w') as inf:
                        #         inf.write(json.dumps(self.infos))
                        # else:
                        #     with open('dlim_tmp.inf', 'w') as inf:
                        #         inf.write(json.dumps(self.infos))

            if recursive_check and self.parent is not None:
                self.parent.inf(check_hash=True)
                
            if self.verbose:
                _prog.end(0, 'generated inf for {}'.format(self.fn))
                
        if 'src_srs' not in self.infos.keys() or self.infos['src_srs'] is None:
            self.infos['src_srs'] = self.src_srs
        else:
            self.src_srs = self.infos['src_srs']

        if 'format' not in self.infos.keys():
            self.infos['format'] = self.data_format

        return(self.infos)

    def set_yield(self):
        """set the yield strategy, either default block or mask"""

        if self.x_inc is not None:
            self.x_inc = utils.str2inc(self.x_inc)
            if self.y_inc is None:
                self.y_inc = self.x_inc
            else:
                self.y_inc = utils.str2inc(self.y_inc)
                
            self.xyz_yield = self.block_xyz()
        else:
            self.xyz_yield = self.yield_xyz()
    
    def set_transform(self):
        """set an srs transform, if needed/wanted."""

        if self.src_srs == '': self.src_srs = None
        if self.dst_srs == '': self.dst_srs = None
        if self.dst_srs is not None and \
           self.src_srs is not None and \
           self.src_srs.split('+')[0] != self.dst_srs.split('+')[0]:
            src_srs = osr.SpatialReference()
            src_srs.SetFromUserInput(self.src_srs)
            dst_srs = osr.SpatialReference()
            dst_srs.SetFromUserInput(self.dst_srs)
            try:
                src_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
                dst_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
            except:
                pass

            self.dst_trans = osr.CoordinateTransformation(src_srs, dst_srs)
            if self.region is not None:
                self.trans_region = self.region.copy()
                self.trans_region.src_srs = self.dst_srs
                self.trans_region.warp(self.src_srs)
    
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
                        if len(xyz_dataset.fn.split('.')) > 1:
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

    def block_array(
            self, want_mean=True, want_count=True, want_weights=True,
            want_mask=True, want_gt=False, want_ds_config=False, min_count=None, out_name=None
    ):
        """block the xyz data to the mean block value

        Yields:
          list: an array for [weighted]mean, weights mask, and count
        """

        out_arrays = {'mean':None, 'count':None, 'weight':None, 'mask':None}
        
        ds_region = regions.Region().from_list(self.infos['minmax'])
        if self.dst_trans is not None:
            ds_region.transform(self.dst_trans)
            
        if self.region is not None and self.region.valid_p():
            block_region = regions.regions_reduce(
                self.region, ds_region
            )
        else:
            block_region = ds_region.copy()

        if block_region.valid_p():
        
            xcount, ycount, dst_gt = block_region.geo_transform(
                x_inc=self.x_inc, y_inc=self.y_inc
            )
            gdt = gdal.GDT_Float32
            count_array = np.zeros((ycount, xcount))
            
            if want_mean:
                z_array = np.zeros((ycount, xcount))
                weight_array = np.zeros((ycount, xcount))

            if self.verbose:
                utils.echo_msg(
                    'blocking data to {}/{} grid'.format(ycount, xcount)
                )
                
            for this_xyz in self.yield_xyz():
                xpos, ypos = utils._geo2pixel(
                    this_xyz.x, this_xyz.y, dst_gt
                )
                try:
                    count_array[ypos, xpos] += 1
                    if want_mean:
                        z_array[ypos, xpos] += this_xyz.z * this_xyz.w
                        
                    if want_weights and self.weight is not None:
                        weight_array[ypos, xpos] += this_xyz.w

                except: pass

            count_array[count_array == 0] = np.nan
            if not self.weight:
                weight_array = np.ones((ycount, xcount))
                
            if want_mean and want_weights:
                weight_array[weight_array == 0] = np.nan
                weight_array = (weight_array/count_array)
                z_array = (z_array/weight_array)/count_array
                
            elif want_mean:
                z_array = (z_array/count_array)

            if min_count is not None and want_mean:
                z_array[count_array < min_count] = np.nan
                weight_array[count_array < min_count] = np.nan

            if want_mean:
                z_array[np.isnan(z_array)] = -9999
                out_arrays['mean'] = z_array
                
            if want_weights:
                weight_array[np.isnan(weight_array)] = -9999
                out_arrays['weight'] = weight_array
                
            if want_mask:
                mask_array = count_array[~np.isnan(count_array)]
                out_arrays['mask'] = mask_array
                
            if want_count:
                count_array[np.isnan(count_array)] = -9999
                out_arrays['count'] = count_array

            ds_config = demfun.set_infos(
                xcount,
                ycount,
                xcount * ycount,
                dst_gt,
                self.dst_srs,
                gdal.GDT_Float32,
                -9999,
                'GTiff'
            )

            return(out_arrays, dst_gt if want_gt else None, ds_config if want_ds_config else None)
        else:
            utils.echo_error_msg('invalid region {}'.format(block_region))
    
    def block_xyz(self, want_gmt=False):
        """block the src_xyz data to the mean block value

        Yields:
          list: xyz data for each block with data
        """
        
        if want_gmt:
            ds_region = regions.Region().from_list(self.infos['minmax'])
            if self.dst_trans is not None:
                ds_region.transform(self.dst_trans)

            if self.region is not None and self.region.valid_p():
                block_region = regions.regions_reduce(
                    self.region, ds_region
                )
            else:
                block_region = ds_region.copy()

            if block_region.valid_p():                
                if utils.config_check()['GMT'] is None:
                    utils.echo_error_msg(
                        'GMT must be installed to use blockmedian; install GMT or set `want_gmt` to False'
                    )
                else:
                    xyz_func = lambda p: self.dump_xyz(dst_port=p, encode=True)
                    for xyz in utils.yield_cmd(
                            'gmt blockmedian -I{:.10f}/{:.10f} {} -r -V'.format(
                                self.x_inc, self.y_inc, block_region.format('gmt')
                            ),
                            verbose=self.verbose,
                            data_fun=xyz_func
                    ):
                        yield(xyzfun.XYZPoint().from_list([float(x) for x in xyz.split()]))
        else:
            block_arrays, dst_gt, ds_config = self.block_array(
                want_mean=True, want_weights=True, want_count=False, want_mask=False,
                want_gt=True, want_ds_config=True
            )
            z_array = block_arrays['mean']
            weight_array = block_arrays['weight']
            for y in range(0, ds_config['ny']):
                for x in range(0, ds_config['nx']):
                    z = z_array[y,x]
                    if not np.isnan(z) and '{:g}'.format(ds_config['ndv']) != '{:g}'.format(z):
                        geo_x, geo_y = utils._pixel2geo(x, y, dst_gt)
                        out_xyz = xyzfun.XYZPoint(
                            x=geo_x, y=geo_y, z=z, w=weight_array[y,x]
                        )
                        yield(out_xyz)

    def mask_xyz(self, dst_x_inc, dst_y_inc, dst_format='MEM', **kwargs):
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
            utils.sr_wkt(self.src_srs),
            gdal.GDT_Int32,
            -9999,
            'MEM'
        )
        for this_xyz in self.yield_xyz_from_entries(**kwargs):
            xpos, ypos = utils._geo2pixel(
                this_xyz.x, this_xyz.y, dst_gt
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
        
    def mask_and_yield_xyz(self, dst_gdal, dst_inc, dst_format='GTiff', **kwargs):
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
            -9999,
            'GTiff'
        )

        for this_xyz in self.yield_xyz_from_entries(**kwargs):
            yield(this_xyz)
            xpos, ypos = utils._geo2pixel(
                this_xyz.x, this_xyz.y, dst_gt
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
        self.xpos = xpos
        self.ypos = ypos
        self.zpos = zpos
        self.wpos = wpos
        self.skip = skip
        self.x_scale = x_scale
        self.y_scale = y_scale
        self.z_scale = z_scale
        self.x_offset = x_offset
        self.y_offset = y_offset
        self._known_delims = [None, ',', '/', ':'] ## space and tab are 'None'
        self.delim = delim
        if delim is not None:
            self._known_delims.insert(0, delim)
        self.scoff = True if x_scale != 1 or y_scale != 1 or z_scale != 1 or x_offset != 0 or y_offset != 0 else False

        super().__init__(**kwargs)
        
    def generate_inf(self, callback=lambda: False):
        """generate a infos dictionary from the xyz dataset"""
                
        pts = []
        self.infos['name'] = self.fn
        self.infos['hash'] = self.hash()#dl_hash(self.fn)
        self.infos['numpts'] = 0
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
        return(self.infos)

    def line_delim(self, xyz_line):
        """guess a line delimiter"""
        
        for delim in self._known_delims:
            this_xyz = xyz_line.split(delim)
            if len(this_xyz) > 1:
                return(this_xyz)

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
            if count >= skip:
                this_xyz = self.line_delim(xyz_line)
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
                    utils.echo_error_msg(e)
                    this_xyz = xyzfun.XYZPoint()

                if this_xyz.valid_p():
                    if self.scoff:
                        this_xyz.x = (this_xyz.x+self.x_offset) * self.x_scale
                        this_xyz.y = (this_xyz.y+self.y_offset) * self.y_scale
                        this_xyz.z *= self.z_scale

                    this_xyz.w = w if self.weight is None else self.weight * w                        
                    if self.dst_trans is not None:
                        this_xyz.transform(self.dst_trans)

                    if self.region is not None and self.region.valid_p():
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

    ## Testing_
    def line_delim_(self, xyz_line):
        """guess a line delimiter"""
        
        for delim in self._known_delims:
            this_xyz = xyz_line.split(delim)
            if len(this_xyz) > 1:
                return(delim)
            
    def yield_xyz_np(self):
        """LAS file parsing generator"""

        count = 0
        self.delim = " "
        # if self.delim is None:
        #     if self.fn is not None:
        #         if os.path.exists(str(self.fn)):
        #             self.src_data = open(self.fn, "r")
            
        #             self.delim = self.line_delim(self.src_data.readline())
        #             self.src_data.close()
                
        #dataset = np.loadtxt(self.fn, delimiter=self.delim, usecols=(self.xpos, self.ypos, self.zpos), unpack=False)
        #dataset = np.vstack((x, y, z)).transpose()

        with open(self.fn, 'r') as f:
            dataset = np.array([l for l in csv.reader(f, delimiter=' ')])
        
        if self.region is not None  and self.region.valid_p():
            if self.dst_trans is not None:
                self.region.src_srs = self.dst_srs
                self.region.warp(self.src_srs)

            dataset = dataset[dataset[:,0] > self.region.xmin,:]
            dataset = dataset[dataset[:,0] < self.region.xmax,:]
            dataset = dataset[dataset[:,1] > self.region.ymin,:]
            dataset = dataset[dataset[:,1] < self.region.ymax,:]
            if self.region.zmin is not None:
                dataset = dataset[dataset[:,2] > self.region.zmin,:]

            if self.region.zmax is not None:
                dataset = dataset[dataset[:,2] < self.region.zmax,:]

        for point in dataset:
            count += 1
            this_xyz = xyzfun.XYZPoint(
                x=point[0], y=point[1], z=point[2], w=self.weight
            )
            if self.dst_trans is not None:
                out_xyz.transform(self.dst_trans)

            yield(this_xyz)
        dataset = None

        if self.verbose:
            utils.echo_msg(
                'parsed {} data records from {}{}'.format(
                    count, self.fn, ' @{}'.format(self.weight) if self.weight is not None else ''
                )
            )            

## ==============================================
## ==============================================
class LASFile(ElevationDataset):
    """representing an LAS/LAZ dataset."""

    def __init__(self, classes='0/2/29/40', **kwargs):
        self.classes = [int(x) for x in classes.split('/')]
        super().__init__(**kwargs)

    def generate_inf(self, callback=lambda: False):
        """generate an inf file for a lidar dataset."""
        
        self.infos['name'] = self.fn
        self.infos['hash'] = self.hash()#dl_hash(self.fn)
        self.infos['numpts'] = 0
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

        # ## gather the projection info if it exists - testing
        # print(lasf.header.vlrs[2].record_id)
        # print(lasf.header.vlrs[2].description)
        # #lasf.header.vlrs[2].parse_record_data(lasf.header.vlrs[2].record_data)
        # #print(lasf.header.vlrs[2].record_data_bytes())
        # print(lasf.header.vlrs[2].geo_keys)
        # [print(x.id) for x in lasf.header.vlrs[2].geo_keys]

        
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
                dataset = np.vstack((points.x, points.y, points.z)).transpose()
                if self.region is not None  and self.region.valid_p():
                    if self.dst_trans is not None:
                        self.region.src_srs = self.dst_srs
                        self.region.warp(self.src_srs)
            
                    dataset = dataset[dataset[:,0] > self.region.xmin,:]
                    dataset = dataset[dataset[:,0] < self.region.xmax,:]
                    dataset = dataset[dataset[:,1] > self.region.ymin,:]
                    dataset = dataset[dataset[:,1] < self.region.ymax,:]
                    if self.region.zmin is not None:
                        dataset = dataset[dataset[:,2] > self.region.zmin,:]
                        
                    if self.region.zmax is not None:
                        dataset = dataset[dataset[:,2] < self.region.zmax,:]
                
                for point in dataset:
                    count += 1
                    this_xyz = xyzfun.XYZPoint(
                        x=point[0], y=point[1], z=point[2], w=self.weight
                    )
                    if self.dst_trans is not None:
                        out_xyz.transform(self.dst_trans)
                    
                    yield(this_xyz)
                dataset = None

        if self.verbose:
            utils.echo_msg(
                'parsed {} data records from {}{}'.format(
                    count, self.fn, ' @{}'.format(self.weight) if self.weight is not None else ''
                )
            )
            
## ==============================================
## ==============================================
class RasterFile(ElevationDataset):
    """providing a GDAL raster dataset parser."""

    def __init__(self, mask=None, weight_mask=None,  open_options=None, **kwargs):
        super().__init__(**kwargs)
        self.open_options = open_options.split('/') if open_options is not None else []
        self.mask = mask
        self.weight_mask = weight_mask
        if self.src_srs is None:
            self.src_srs = demfun.get_srs(self.fn)
            
        self.set_transform()

    def init_ds(self):
        if self.x_inc is not None and self.y_inc is not None:
            
            if self.verbose:
                _warp_prog = utils.CliProgress('warping dataset {} to {}/{}'.format(self.fn, self.x_inc, self.y_inc))
            if self.region is not None:
                dem_region = regions.Region().from_list(self.infos['minmax'])
                warp_region = regions.regions_reduce(self.region, dem_region)
                extent = [warp_region.xmin, warp_region.ymin, warp_region.xmax, warp_region.ymax]
            else:
                extent = [self.infos['minmax'][0], self.infos['minmax'][2], self.infos['minmax'][1], self.infos['minmax'][3]]

            src_ds = gdal.Warp('', self.fn, format='MEM', xRes=self.x_inc, yRes=self.y_inc,
                               dstNodata=-9999, outputBounds=extent, resampleAlg='bilinear', targetAlignedPixels=True,
                               options=["COMPRESS=LZW", "TILED=YES"], callback=lambda x,y,z: _warp_prog.update() if self.verbose else None)#gdal.TermProgress)
            if self.verbose:
                _warp_prog.end(0, 'warped dataset {} to {}/{} @ {}/{}/{}/{}'.format(self.fn, self.x_inc, self.y_inc, *extent))
        else:
            if self.open_options:
                src_ds = gdal.OpenEx(self.fn, open_options=self.open_options)
            else:
                src_ds = gdal.Open(self.fn)

        return(src_ds)
        
    def generate_inf(self, callback=lambda: False):
        """generate a infos dictionary from the raster dataset"""
            
        self.infos['name'] = self.fn
        self.infos['hash'] = self.hash()#dl_hash(self.fn)
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

    def get_srcwin(self, gt, x_size, y_size, node='pixel'):
        if self.region is not None:
            if self.dst_trans is not None:
                if self.trans_region is not None and self.trans_region.valid_p(
                        check_xy = True
                ):
                    srcwin = self.trans_region.srcwin(
                        gt, x_size, y_size, node
                    )
                else:
                    srcwin = (
                        0, 0, x_size, y_size
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
    
    def yield_xyz(self):
        """parse the data from gdal dataset src_ds (first band only)"""

        src_ds = self.init_ds()
        out_xyz = xyzfun.XYZPoint(w=1)
        if src_ds is not None:
            count = 0
            band = src_ds.GetRasterBand(1)
            gt = src_ds.GetGeoTransform()
            ndv = band.GetNoDataValue()
            msk_band = None
            weight_band = None
            if self.mask is not None:
                src_mask = gdal.Open(self.mask)
                msk_band = src_mask.GetRasterBand(1)
            if self.weight_mask is not None:

                if self.x_inc is not None and self.y_inc is not None:
                    src_weight = gdal.Warp('', self.weight_mask, format='MEM', xRes=self.x_inc, yRes=self.y_inc,
                                           dstNodata=-9999, outputBounds=extent, resampleAlg='bilinear', targetAlignedPixels=True,
                                           options=["COMPRESS=LZW", "TILED=YES"])
                else:
                    src_weight = gdal.Open(self.weight_mask)
                    
                weight_band = src_weight.GetRasterBand(1)
            
            nodata = ['{:g}'.format(-9999), 'nan', float('nan')]
            if ndv is not None:
                nodata.append('{:g}'.format(ndv))

            srcwin = self.get_srcwin(gt, src_ds.RasterXSize, src_ds.RasterYSize)
            for y in range(
                    srcwin[1], srcwin[1] + srcwin[3], 1
            ):
                # if self.verbose and self.parent is None:
                #     #if self.verbose:
                #     _prog.update_perc((y, srcwin[1] + srcwin[3]))
                    
                band_data = band.ReadAsArray(
                    srcwin[0], y, srcwin[2], 1
                )
                if weight_band is not None:
                    weight_data = weight_band.ReadAsArray(
                        srcwin[0], y, srcwin[2], 1
                    )
                if self.region is not None and self.region.valid_p():
                    z_region = self.region.z_region()
                    if z_region[0] is not None:
                        band_data[band_data < z_region[0]] = -9999
                        
                    if z_region[1] is not None:
                        band_data[band_data > z_region[1]] = -9999

                    if weight_band is not None:
                        w_region = self.region.w_region()
                        if w_region[0] is not None:
                            band_data[weight_data < w_region[0]] = -9999
                        
                        if w_region[1] is not None:
                            band_data[weight_data > w_region[1]] = -9999
                        
                if msk_band is not None:
                   msk_data = msk_band.ReadAsArray(
                       self.srcwin[0], y, self.srcwin[2], 1
                   )
                   band_data[msk_data==0]=-9999
                   
                band_data = np.reshape(band_data, (srcwin[2], ))
                if weight_band is not None:
                    weight_data = np.reshape(weight_data, (srcwin[2], ))
                for x_i in range(0, srcwin[2], 1):
                    z = band_data[x_i]
                    if '{:g}'.format(z) not in nodata:
                        x = x_i + srcwin[0]
                        out_xyz.x, out_xyz.y = utils._pixel2geo(x, y, gt)
                        out_xyz.z = z
                        if weight_band is not None:
                            out_xyz.w = weight_data[x_i]
                        else:
                            out_xyz.w = utils.float_or(self.weight, 1)
                        count += 1
                        if self.dst_trans is not None:
                            out_xyz.transform(self.dst_trans)

                        yield(out_xyz)
                            
            band = msk_band = weight_band = src_weight = src_mask = None
            if self.verbose:
                utils.echo_msg(
                    'parsed {} data records from {}{}'.format(
                        count, self.fn, ' @{}'.format(self.weight) if self.weight is not None else ''
                    )
                )
                
        src_ds = None
        # if self.verbose and self.parent is None:
        #     _prog.end(0, 'parsed dataset {}{}'.format(self.fn, ' @{}'.format(self.weight) if self.weight is not None else ''))

    ## TODO w region filter
    ## TODO weight band
    ## TODO mask band
    def yield_array(self):
        src_ds = self.init_ds()
        if src_ds is not None:
            band = src_ds.GetRasterBand(1)
            gt = src_ds.GetGeoTransform()
            ndv = float(band.GetNoDataValue())
            dem_infos = demfun.gather_infos(src_ds)
            srcwin = self.get_srcwin(gt, src_ds.RasterXSize, src_ds.RasterYSize)
            src_arr = band.ReadAsArray(srcwin[0],srcwin[1],srcwin[2],srcwin[3])
            x_count, y_count, dst_gt = self.region.geo_transform(self.x_inc, self.y_inc)
            srcwin_region = regions.Region().from_geo_transform(geo_transform=gt, x_count=src_arr.shape[1], y_count=src_arr.shape[0])
            dst_srcwin = srcwin_region.srcwin(dst_gt, x_count, y_count)
            src_arr[src_arr == dem_infos['ndv']] = np.nan

            if self.region is not None and self.region.valid_p():
                z_region = self.region.z_region()
                if z_region[0] is not None:
                    src_arr[src_arr < z_region[0]] = np.nan
                        
                if z_region[1] is not None:
                    src_arr[src_arr > z_region[1]] = np.nan
            
            yield(src_arr, dst_srcwin, dst_gt, self.weight)
        src_ds = None

        if self.verbose:
            utils.echo_msg(
                'parsed {} data records from {}{}'.format(
                    (dst_srcwin[3] * dst_srcwin[2]), self.fn, ' @{}'.format(self.weight) if self.weight is not None else ''
                )
            )
        
    def yield_xyz_np(self):
        """parse the data from gdal dataset src_ds using numpy - testing"""

        if self.x_inc is not None and self.y_inc is not None:
            src_ds = gdal.Open(self.fn)
            gt = src_ds.GetGeoTransform()
            src_ds = None
            if '{:g}'.format(gt[1]) != '{:g}'.format(self.x_inc) and \
               '{:g}'.format(gt[5]*-1) != '{:g}'.format(self.y_inc):
                if self.verbose:
                    _warp_prog = utils.CliProgress('warping dataset {} to {}/{}'.format(self.fn, self.x_inc, self.y_inc))
                if self.region is not None:
                    extent = [self.region.xmin, self.region.ymin, self.region.xmax, self.region.ymax]
                else:
                    extent = [self.infos['minmax'][0], self.infos['minmax'][2], self.infos['minmax'][1], self.infos['minmax'][3]]
                src_ds = gdal.Warp('', self.fn, format='MEM', xRes=self.x_inc, yRes=self.y_inc,
                                   dstNodata=-9999, outputBounds=extent, resampleAlg='bilinear',
                                   options=["COMPRESS=LZW", "TILED=YES"], callback=lambda x,y,z: _warp_prog.update() if self.verbose else None)#gdal.TermProgress)
                if self.verbose:
                    _warp_prog.end(0, 'warped dataset {} to {}/{}'.format(self.fn, self.x_inc, self.y_inc))
        else:
        
            if self.open_options:
                src_ds = gdal.OpenEx(self.fn, open_options=self.open_options)
            else:
                src_ds = gdal.Open(self.fn)

        if src_ds is not None:
            count = 0
            band = src_ds.GetRasterBand(1)
            gt = src_ds.GetGeoTransform()
            ndv = float(band.GetNoDataValue())
            dem_infos = demfun.gather_infos(src_ds)

            if self.weight_mask is not None:
                if self.x_inc is not None and self.y_inc is not None and \
                   '{:g}'.format(gt[1]) != '{:g}'.format(self.x_inc) and \
                       '{:g}'.format(gt[5]*-1) != '{:g}'.format(self.y_inc):
                    if self.verbose:
                        _warp_prog = utils.CliProgress('warping dataset {} to {}/{}'.format(self.weight_mask, self.x_inc, self.y_inc))
                    if self.region is not None:
                        extent = [self.region.xmin, self.region.ymin, self.region.xmax, self.region.ymax]
                    else:
                        extent = [self.infos['minmax'][0], self.infos['minmax'][2], self.infos['minmax'][1], self.infos['minmax'][3]]
                    src_weight = gdal.Warp('', self.weight_mask, format='MEM', xRes=self.x_inc, yRes=self.y_inc,
                                           dstNodata=-9999, outputBounds=extent, resampleAlg='bilinear',
                                           options=["COMPRESS=LZW", "TILED=YES"], callback=lambda x,y,z: _warp.update() if self.verbose else None)#gdal.TermProgress)
                    if self.verbose:
                        _warp.prog.end(0, 'warped dataset {} to {}/{}'.format(self.weight_mask, self.x_inc, self.y_inc))
                else:
                    src_weight = gdal.Open(self.weight_mask)
                    
                weight_band = src_weight.GetRasterBand(1)
            
            for srcwin in utils.yield_srcwin((dem_infos['ny'], dem_infos['nx']), 8000): #dem_infos['nb']/2):
                src_arr = band.ReadAsArray(srcwin[0],srcwin[1],srcwin[2],srcwin[3])
                src_arr = src_arr.flatten()

                this_geo_x_origin, this_geo_y_origin = utils._pixel2geo(srcwin[0], srcwin[1], gt)
                this_geo_x_end, this_geo_y_end = utils._pixel2geo(srcwin[0]+srcwin[2], srcwin[1]+srcwin[3], gt)
                dst_gt = [this_geo_x_origin, float(gt[1]), 0.0, this_geo_y_origin, 0.0, float(gt[5])]

                srcwin_region = regions.Region().from_geo_transform(geo_transform=dst_gt, x_count=srcwin[2], y_count=srcwin[3])
                
                xi = np.linspace(srcwin_region.xmin, srcwin_region.xmax, srcwin[2])
                yi = np.linspace(srcwin_region.ymax, srcwin_region.ymin, srcwin[3])
                xi, yi = np.meshgrid(xi, yi)
                xi, yi = xi.flatten(), yi.flatten()

                if self.weight_mask is not None:
                    weight_arr = weight_band.ReadAsArray(srcwin[0],srcwin[1],srcwin[2],srcwin[3])
                    weight_arr = weight_arr.flatten()
                    dataset = np.vstack((xi, yi, src_arr, weight_arr)).T
                else:
                    dataset = np.vstack((xi, yi, src_arr)).T

                if self.x_inc is not None and self.y_inc is not None:
                    dataset = dataset[dataset[:,2] != -9999.,:]
                else:
                    dataset = dataset[dataset[:,2] != ndv,:]
                    
                if self.region is not None  and self.region.valid_p():
                    if self.dst_trans is not None:
                        self.region.src_srs = self.dst_srs
                        self.region.warp(self.src_srs)

                    dataset = dataset[dataset[:,0] > self.region.xmin,:]
                    dataset = dataset[dataset[:,0] < self.region.xmax,:]
                    dataset = dataset[dataset[:,1] > self.region.ymin,:]
                    dataset = dataset[dataset[:,1] < self.region.ymax,:]
                    if self.region.zmin is not None:
                        dataset = dataset[dataset[:,2] > self.region.zmin,:]

                    if self.region.zmax is not None:
                        dataset = dataset[dataset[:,2] < self.region.zmax,:]

                    if self.weight_mask is not None:
                        
                        if self.region.wmin is not None:
                            dataset = dataset[dataset[:,3] > self.region.wmin,:]

                        if self.region.wmax is not None:
                            dataset = dataset[dataset[:,3] < self.region.wmax,:]

                for point in dataset:
                    count += 1
                    this_xyz = xyzfun.XYZPoint(
                        x=point[0], y=point[1], z=point[2], w=self.weight if self.weight_mask is None else point[3]
                    )
                    if self.dst_trans is not None:
                        out_xyz.transform(self.dst_trans)
                    
                    yield(this_xyz)
                dataset = src_arr = weight_arr = None

            src_ds = src_weight = None
            if self.verbose:
                utils.echo_msg(
                    'parsed {} data records from {}{}'.format(
                        count, self.fn, ' @{}'.format(self.weight) if self.weight is not None else ''
                    )
                )
            
## ==============================================
## ==============================================
class BAGFile(ElevationDataset):
    """providing a BAG raster dataset parser.

    process supergrids at native resolution if they
    exist, otherwise process as normal grid
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        if self.src_srs is None:
            self.src_srs = demfun.get_srs(self.fn)
            
        self.set_transform()
        
    def generate_inf(self, callback=lambda: False):
        """generate a infos dictionary from the raster dataset"""
            
        self.infos['name'] = self.fn
        self.infos['hash'] = self.hash()#dl_hash(self.fn)
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
            #self.infos['src_srs'] = self.src_srs
            
        src_ds = None
        return(self.infos)

    def parse_(self):
        if self.verbose:
            _prog = utils.CliProgress(
                'parsing BAG {}{}'.format(
                    self.fn,
                    ' @{}'.format(self.weight) if self.weight is not None else '')
            )
        
        mt = gdal.Info(self.fn, format='json')['metadata']['']
        if 'HAS_SUPERGRIDS' in mt.keys() and mt['HAS_SUPERGRIDS'] == 'TRUE':
            oo = ["MODE=LIST_SUPERGRIDS"]
            if self.region is not None and self.region.valid_p():
                if self.dst_trans is not None:
                    tmp_region = self.trans_region.copy()
                else:
                    tmp_region = self.region.copy()
                    
                oo.append('MINX={}'.format(tmp_region.xmin))
                oo.append('MAXX={}'.format(tmp_region.xmax))
                oo.append('MINY={}'.format(tmp_region.ymin))
                oo.append('MAXY={}'.format(tmp_region.ymax))
                    
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
                    src_region=self.region,
                    verbose=self.verbose
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
                verbose=self.verbose
            )
            yield(sub_ds)
                    
        if self.verbose:
            _prog.end(0, 'parsed datalist {}{}'.format(
                self.fn, ' @{}'.format(self.weight) if self.weight is not None else ''
            ))

    def yield_xyz(self):
        for ds in self.parse_():
            for xyz in ds.yield_xyz():
                yield(xyz)
            
## ==============================================
## ==============================================
class MBSParser(ElevationDataset):
    """providing an mbsystem parser"""

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.infos = {}
        
    def generate_inf(self, callback=lambda: False):
        self.infos['name'] = self.fn
        self.infos['hash'] = None
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

    def yield_xyz(self):
        for line in utils.yield_cmd(
                'mblist -MA -OXYZ -I{}'.format(self.fn),
                verbose=True,
        ):
            this_xyz = xyzfun.XYZPoint().from_string(line, delim='\t')
            this_xyz.weight = self.weight
            yield(this_xyz)
### End
