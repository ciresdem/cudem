### dlim.py - DataLists IMproved
##
## Copyright (c) 2010 - 2021 CIRES Coastal DEM Team
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
## data-path data-format data-weight data-name data-source data-date data-resolution data-type data-horz data-vert data-url
## Minimally, data-path is all that is needed.
##
### Code:

import os
import sys
import re
import json
import glob
import hashlib
import gdal
import ogr
import osr
import numpy as np
from scipy import spatial
import urllib
import requests
from cudem import utils
from cudem import regions
from cudem import fetches
from cudem import xyzfun
from cudem import demfun

__version__ = '0.1.4'

r_headers = { 'User-Agent': 'DLIM v%s' %(__version__) }

def fetch_req(src_url, params=None, tries=5, timeout=2, read_timeout=10):
    """fetch src_url and return the requests object"""
    
    if tries <= 0:
        utils.echo_error_msg('max-tries exhausted')
        return(None)
    try:
        return(requests.get(src_url, stream=True, params=params, timeout=(timeout,read_timeout), headers=r_headers))
    except:
        return(fetch_req(src_url, params=params, tries=tries - 1, timeout=timeout + 1, read_timeout=read_timeout + 10))

def fetch_file(src_url, dst_fn, params=None, callback=lambda: False, datatype=None, overwrite=False,
               verbose=False, timeout=140, read_timeout=320):
    """fetch src_url and save to dst_fn"""
    
    status = 0
    halt = callback
    req = None

    if verbose:
        progress = utils.CliProgress('fetching remote file: {}...'.format(os.path.basename(src_url)[:20]))
    if not os.path.exists(os.path.dirname(dst_fn)):
        try:
            os.makedirs(os.path.dirname(dst_fn))
        except: pass 
    if not os.path.exists(dst_fn) or overwrite:
        try:
            with requests.get(src_url, stream=True, params=params, headers=r_headers,
                              timeout=(timeout,read_timeout)) as req:
                req_h = req.headers
                if req.status_code == 200:
                    curr_chunk = 0
                    with open(dst_fn, 'wb') as local_file:
                        for chunk in req.iter_content(chunk_size = 8196):
                            if halt():
                                break
                            if verbose:
                                progress.update()
                            if chunk:
                                local_file.write(chunk)
                else:
                    utils.echo_error_msg('server returned: {}'.format(req.status_code))
        except Exception as e:
            utils.echo_error_msg(e)
            status = -1
    if not os.path.exists(dst_fn) or os.stat(dst_fn).st_size ==  0:
        status = -1
    if verbose:
        progress.end(status, 'fetched remote file: {}.'.format(os.path.basename(dst_fn)[:20]))
    return(status)

class XYZDataset():
    """representing an xyz-able parser or a data-list entry."""

    def __init__(self, fn=None, data_format=None, weight=1, epsg=4326, name="<XYZDataset>", title=None,
                 source=None, date=None, data_type=None, resolution=None, vdatum=None, url=None,
                 parent=None, src_region=None, warp=None, verbose=False, remote=False):

        self.name = name
        self.title = title
        self.data_format = data_format
        self.epsg = utils.int_or(epsg)
        self.infos = {}
        self.source = source
        self.date = date
        self.data_type = data_type
        self.resolution = resolution
        self.vdatum = vdatum
        self.url = url
        self.parent = parent
        self.weight = weight
        self.verbose = verbose
        self.fn = fn
        self.data_entries = []
        self.data_lists = {}
        self.region = src_region
        self.warp = utils.int_or(warp)
        self.remote = remote
                
    def valid_p(self):
        """check if self is a valid datalist entry

        Returns:
          bools: True if valid else False
        """
        
        if self.fn is None: return(False)
        if self.data_format is None: return(False)
        if self.fn is not None:
            if self.fn not in DatasetFactory().data_types[self.data_format]['fmts']:
                if not self.fn.startswith('http'):
                    if not os.path.exists(self.fn): return (False)
                    if os.stat(self.fn).st_size == 0: return(False)
        return(True)
        
    def hash(self, sha1=False):
        """generate a hash of the xyz-dataset source file

        Returns:
          str: hexdigest
        """
        
        BUF_SIZE = 65536  # lets read stuff in 64kbchunks!
        if sha1: this_hash = hashlib.sha1()
        else: this_hash = hashlib.md5()

        try:
            with open(self.fn, 'rb') as f:
                while True:
                    data = f.read(BUF_SIZE)
                    if not data:
                        break
                    this_hash.update(data)
                    
            return(this_hash.hexdigest())
        except: return('0')

    def echo_(self, **kwargs):
        """print self as a datalist entry string"""

        return(' '.join(self.fn, self.data_format, self.weight, self.title, self.source,
                        self.date, self.data_type, self.resolution, self.epsg, self.vdatum,
                        self.url))

    def format_metadata(self, **kwargs):
        return('{} {} {} {} {} {} {} {} {}'.format(
            self.weight, self.title, self.source,
            self.date, self.data_type, self.resolution, self.epsg, self.vdatum,
            self.url))
    
    def echo(self, **kwargs):
        """print self as a datalist entry string."""

        for entry in self.data_entries:
            l = [entry.fn, entry.data_format]
            if entry.weight is not None: l.append(entry.weight)
            print('{}'.format(" ".join([str(x) for x in l])))

    def generate_inf(self):
        pass
            
    def inf(self, check_hash = False, **kwargs):
        """read/write an inf file

        Args:
          kwargs (dict): any arguments to pass to the dataset parser
        
        Returns:
          dict: xyz-dataset info dictionary
        """
        
        inf_path = '{}.inf'.format(self.fn)

        if os.path.exists(inf_path):
            try:
                with open(inf_path) as i_ob:
                    self.infos = json.load(i_ob)
            except ValueError:
                try:
                    self.infos = MBSParser(
                        fn = inf_path, epsg = self.epsg).inf_parse().infos
                except:
                    utils.echo_error_msg('failed to parse inf {}'.format(inf_path))
            except:
                utils.echo_error_msg('failed to parse inf {}'.format(inf_path))
        else:
            self.infos = {}
        
        if check_hash:
            if 'hash' in self.infos.keys():
                gen_inf = self.hash() != self.infos['hash']
            else:
                gen_inf = True
        else:
            gen_inf = 'hash' not in self.infos.keys() or 'wkt' not in self.infos.keys()

        if gen_inf:
            utils.echo_msg("generating inf for {}".format(self.fn))
            self.infos = self.generate_inf()
            #print(self.infos)
            self.infos['format'] = self.data_format
            if 'minmax' in self.infos:
                if self.infos['minmax'] is not None:
                    try:
                        with open('{}.inf'.format(self.fn), 'w') as inf:
                            inf.write(json.dumps(self.infos))
                    except:
                        #if self.parent is not None:
                        with open('{}_{}.inf'.format(
                                'dlim_tmp', self.region.format('fn')), 'w') as inf:
                            inf.write(json.dumps(self.infos))
            if self.parent is not None:
                utils.remove_glob('{}.inf'.format(self.parent.fn))
                self.parent.infos = {}
            self.infos['epsg'] = self.epsg
        self.infos['format'] = self.data_format
        return(self.infos)

    def parse(self):
        pass

    def parse_data_lists(self):
        for e in self.data_entries:
            if e.parent is not None:
                if e.parent.name in self.data_lists.keys():
                    self.data_lists[e.parent.name]['data'].append(e)
                else:
                    self.data_lists[e.parent.name] = {'data': [e], 'parent': e.parent}
            else:
                self.data_lists[e.name] = {'data': [e], 'parent': e}
        return(self)
                
    def yield_xyz(self):
        pass

    def archive_xyz(self, **kwargs):

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
            a_name = self.name
        else:
            a_name = '{}_{}_{}'.format(
                self.name, self.region.format('fn'), utils.this_year())

        self.parse_data_lists()
        
        with open('{}.datalist'.format(a_name), 'w') as dlf:
            for x in self.data_lists.keys():
                a_dir = '{}_{}_{}'.format(x, self.region.format('fn'), utils.this_year())
                this_dir = xdl2dir(self.data_lists[x]['parent'])
                this_dir.append(a_dir)
                tmp_dir = this_dir
                #dlf.write('{}.datalist -1 {}\n'.format(os.path.join(*this_dir, this_dir[-1]), self.data_lists[x]['parent'].weight))
                dlf.write('{}.datalist -1 {}\n'.format(os.path.join(*this_dir, this_dir[-1]), self.data_lists[x]['parent'].format_metadata()))
                this_dir = os.path.join(os.getcwd(), *this_dir)
                if not os.path.exists(this_dir):
                    os.makedirs(this_dir)
                with open(os.path.join(this_dir, '{}.datalist'.format(os.path.basename(this_dir))), 'w') as sub_dlf:
                    for xyz_dataset in self.data_lists[x]['data']:
                        sub_xyz_path = '.'.join(
                            [utils.fn_basename(os.path.basename(utils.slugify(xyz_dataset.fn)),
                                               xyz_dataset.fn.split('.')[-1]),
                             'xyz'])
                        
                        this_xyz_path = os.path.join(this_dir, sub_xyz_path)
                        sub_dlf.write('{} 168\n'.format(sub_xyz_path))
                        
                        with open(this_xyz_path, 'w') as xp:
                            for this_xyz in xyz_dataset.yield_xyz(**kwargs):
                                yield(this_xyz)
                                this_xyz.dump(include_w=True if self.weight is not None else False,
                                              dst_port=xp, encode=False)

                                #xyz_dataset.dump_xyz(dst_port=xp, **kwargs)
        Datalist(fn='{}.datalist'.format(a_name)).parse()
            
    def mask_xyz(self, dst_gdal, dst_inc, dst_format='GTiff', **kwargs):
        """Create a num grid mask of xyz data. The output grid
        will contain 1 where data exists and 0 where no data exists.

        yields the xyz data
        """

        xcount, ycount, dst_gt = self.region.geo_transform(x_inc=dst_inc)
        ptArray = np.zeros((ycount, xcount))
        ds_config = demfun.set_infos(
            xcount, ycount, (xcount*ycount), dst_gt, utils.sr_wkt(self.epsg),
            gdal.GDT_Float32, -9999, 'GTiff')
        for this_xyz in self.yield_xyz(**kwargs):
            yield(this_xyz)
            if regions.xyz_in_region_p(this_xyz, self.region):
                xpos, ypos = utils._geo2pixel(this_xyz.x, this_xyz.y, dst_gt)
                try:
                    ptArray[ypos, xpos] = 1
                except:
                    pass
        out, status = utils.gdal_write(ptArray, dst_gdal, ds_config)    
        
    def dump_xyz(self, dst_port=sys.stdout, encode=False, **kwargs):
        for this_xyz in self.yield_xyz(**kwargs):
            this_xyz.dump(include_w=True if self.weight is not None else False,
                          dst_port=dst_port, encode=False)

class MBSParser:
    """providing an mbsystem parser"""

    def __init__(self, fn=None, epsg=None):
        self.fn = fn
        self.epsg = utils.int_or(epsg)
        self.infos = {}
        
    def inf(self):
        pass

    def parse(self):
        pass

    def yield_xyz(self):
        pass
    
    def inf_parse(self):

        self.infos['name'] = self.fn
        self.infos['minmax'] = [0,0,0,0,0,0]
        self.infos['hash'] = None
        dims = []
        this_row = 0

        with open(self.fn) as iob:
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
            xcount, ycount, dst_gt = mbs_region.geo_transform(x_inc=xinc, y_inc=yinc)

            ds_config = {'nx': dims[0], 'ny': dims[1], 'nb': dims[1] * dims[0],
                         'geoT': dst_gt, 'proj': utils.sr_wkt(self.epsg),
                         'dt': gdal.GDT_Float32, 'ndv': 0, 'fmt': 'GTiff'}

            driver = gdal.GetDriverByName('MEM')
            ds = driver.Create('tmp', ds_config['nx'], ds_config['ny'], 1, ds_config['dt'])
            ds.SetGeoTransform(ds_config['geoT'])
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

class Datalist(XYZDataset):
    """representing a datalist parser
    
    A datalist is an MB-System style datalist.
    """

    def __init__(self, fmt=None, **kwargs):
        super().__init__(**kwargs)
        self.name = os.path.basename('.'.join(self.fn.split('.')[:-1]))
        
    def generate_inf(self):
        """return the region of the datalist and generate
        an associated `.inf` file if `inf_file` is True.

        Args:
          dl (str): a datalist pathname
          inf_file (bool): generate an inf file
          epsg (int): EPSG code
          overwrite (bool): overwrite a possibly existing inf_file

        Returns:
          list: the region [xmin, xmax, ymin, ymax]
        """

        _region = self.region
        self.region = None
        self.parse()
        out_regions = []
        out_region = None
        self.infos['name'] = self.fn
        self.infos['numpts'] = 0
        self.infos['hash'] = self.hash()#dl_hash(self.fn)
        for entry in self.data_entries:
            out_regions.append(entry.infos['minmax'])
            self.infos['numpts'] += entry.infos['numpts']

        l = 0

        for this_region in out_regions:
            if l == 0:
                tmp_region = regions.Region().from_list(this_region)
                if tmp_region.valid_p():
                    out_region = regions.Region().from_list(this_region)
                    l += 1
            else:
                tmp_region = regions.Region().from_list(this_region)
                if tmp_region.valid_p():
                    out_region = regions.regions_merge(out_region, tmp_region)
        if out_region is not None:
            self.infos['minmax'] = out_region.export_as_list(include_z=True)
            self.infos['wkt'] = out_region.export_as_wkt()
        else:
            self.infos['minmax'] = None

        self.region = _region
        return(self.infos)
    
    def parse(self):
        """import a datalist entry from a string
    
        Returns:
          datalist_parser: self
        """
        
        _prog = utils.CliProgress('parsing datalist {}'.format(self.fn))
        if os.path.exists(self.fn):
            with open(self.fn, 'r') as op:
                for this_line in op:
                    _prog.update()
                    if this_line[0] != '#' and this_line[0] != '\n' and this_line[0].rstrip() != '':
                        data_set = DatasetFactory(
                            this_line, parent=self, name=self.name, src_region=self.region, source=self.source, date=self.date,
                            data_type=self.data_type, resolution=self.resolution, vdatum=self.vdatum, url=self.url, title=self.title,
                            warp=self.warp, verbose=self.verbose).acquire_dataset()
                        if data_set is not None and data_set.valid_p():
                            data_set.inf()
                            if self.region is not None and self.region.valid_p(check_xy=True):
                                if data_set.infos['minmax'] is not None:
                                    inf_region = regions.Region().from_string(data_set.infos['wkt'])
                                    if regions.regions_intersect_p(inf_region, self.region):
                                        data_set.parse()
                                        for entry in data_set.data_entries:
                                            self.data_entries.append(entry)
                            else:
                                data_set.parse()
                                for entry in data_set.data_entries:
                                    self.data_entries.append(entry)
        else: echo_warning_msg('could not open datalist/entry {}'.format(self.fn))
        self.parse_data_lists()
        _prog.end(0, 'parsed datalist {}'.format(self.fn))
        return(self)
           
    def yield_xyz(self):
        """parse the data from the datalist

        Yields:
          xyz: the parsed xyz data
        """

        for i, this_entry in enumerate(self.data_entries):
            for xyz in this_entry.yield_xyz():
                yield(xyz)
            if this_entry.remote:
                utils.remove_glob('{}*'.format(this_entry.fn))
    
class XYZFile(XYZDataset):
    """representing an xyz dataset stream
    """

    def __init__(self, delim=None, xpos=0, ypos=1, zpos=2, skip=0, x_scale=1, y_scale=1,
                 z_scale=1, x_offset=0, y_offset=0, **kwargs):
        
        super().__init__(**kwargs)
        self.delim = delim
        self.xpos = xpos
        self.ypos = ypos
        self.zpos = zpos
        self.skip = skip
        self.x_scale = x_scale
        self.y_scale = y_scale
        self.z_scale = z_scale
        self.x_offset = x_offset
        self.y_offset = y_offset
        
        self._known_delims = [',', '/', ':']
        self._known_fmts = ['xyz', 'csv', 'dat', 'ascii']

    def generate_inf(self):
        """generate a infos dictionary from the xyz dataset

        Returns:
          dict: a data-entry infos dictionary
        """
                
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
                #try:
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
                #except: pass
            pts.append(l.export_as_list(include_z = True))
            self.infos['numpts'] = i

        self.infos['minmax'] = this_region.export_as_list(include_z = True)
        if self.infos['numpts'] > 0:
            try:
                out_hull = [pts[i] for i in spatial.ConvexHull(pts, qhull_options='Qt').vertices]
                out_hull.append(out_hull[0])
                self.infos['wkt'] = create_wkt_polygon(out_hull, xpos=0, ypos=1)
            except:
                self.infos['wkt'] = this_region.export_as_wkt()
        self.region = region_
        return(self.infos)
        
    def line_delim(self, xyz_line):
        """guess a line delimiter
        Args:
          xyz_line (str): a string representing delimited data.

        Returns:
          str: delimiter (or None)
        """

        for delim in self._known_delims:
            this_xyz = xyz_line.split(delim)
            if len(this_xyz) > 1:
                self.delim = delim

    def parse(self):
        
        if self.region is not None:
            self.inf()
            inf_region = regions.Region().from_string(self.infos['wkt'])
            if regions.regions_intersect_p(inf_region, self.region):
                self.data_entries = [self]
        else: self.data_entries = [self]
        return(self)
                
    def yield_xyz(self):
        """xyz file parsing generator

        Yields:
          xyz: xyz data
        """
        
        if self.fn is not None:
            if os.path.exists(str(self.fn)):
                self.src_data = open(self.fn, "r")
            else:
                self.src_data = sys.stdin
        else:
            self.src_data = sys.stdin
        
        sk = self.skip
        this_xyz = xyzfun.XYZPoint(w = 1)
        if self.region is not None:
            if self.region.epsg != self.epsg:
                if self.warp is not None:
                    if self.region.epsg != self.warp:
                        self.region.warp(warp_epsg=self.epsg)
                else: self.region.warp(warp_epsg=self.epsg)

        warp_epsg = utils.int_or(self.warp)
        if warp_epsg is not  None and self.epsg is not None:
            src_srs = osr.SpatialReference()
            src_srs.ImportFromEPSG(self.epsg)

            dst_srs = osr.SpatialReference()
            dst_srs.ImportFromEPSG(warp_epsg)
            try:
                src_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
                dst_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
            except:
                pass
            dst_trans = osr.CoordinateTransformation(src_srs, dst_srs)
        else: dst_trans = None
        
        ln = 0
        for xyz_line in self.src_data:
            if ln >= sk:
                if self.delim is None: self.line_delim(xyz_line)
                this_xyz.from_string(xyz_line, delim=self.delim,
                                     x_pos=self.xpos, y_pos=self.ypos)
                if this_xyz.valid_p():
                    this_xyz.x = (this_xyz.x+self.x_offset) * self.x_scale
                    this_xyz.y = (this_xyz.y+self.y_offset) * self.y_scale
                    this_xyz.z *= self.z_scale
                    this_xyz.w = self.weight
                    if self.region is not None and self.region.valid_p():
                        if regions.xyz_in_region_p(this_xyz, self.region):
                            if dst_trans is not None:
                                this_xyz.transform(dst_trans)
                            ln += 1
                            yield(this_xyz)
                    else:
                        if dst_trans is not None:
                            this_xyz.transform(dst_trans)
                        ln +=1
                        yield(this_xyz)        
            else: sk -= 1
        if self.verbose: utils.echo_msg('parsed {} data records from {}'.format(ln, self.fn))
        self.src_data.close()

class RasterFile(XYZDataset):
    """providing a raster parser
    """

    def __init__(self, mask=None, step=1, **kwargs):
        super().__init__(**kwargs)
        self.mask = mask
        self.step = 1

        self.src_ds = None
        self.ds_config = None
        self.ds_open_p = True

        if self.remote:
            if self.region is not None:
                outf = '_tmp_dlim_raster_{}.tif'.format(self.region.format('fn'))
            else: outf = '_tmp_dlim_raster.tif'
            if fetch_file(self.fn, outf, verbose = self.verbose) == 0:
                self.fn = outf
            
    def _open_ds(self):
        """open the gdal datasource and gather infos 

        Returns:
          raster_parser: self
        """
        
        if self.fn is not None:
            if os.path.exists(str(self.fn)):
                try:
                    self.src_ds = gdal.Open(self.fn)
                except:
                    self.src_ds = None
            else:
                self.src_ds = None
        else:
            self.src_ds = None

        if self.src_ds is not None:
            self.gather_infos()
            self.ds_open_p = True
        else:
            self.ds_open_p = False
        return(self)

    def _close_ds(self):
        """close the gdal datasource

        Returns:
          raster_parser: self
        """
        
        self.src_ds = None
        self.ds_config = None
        self.ds_open_p = False
        return(self)

    def generate_inf(self):
        """generate a infos dictionary from the raster dataset

        Returns:
          dict: a data-entry infos dictionary
        """
            
        self.infos['name'] = self.fn
        self.infos['hash'] = self.hash()#dl_hash(self.fn)
        self._open_ds()
        gt = self.src_ds.GetGeoTransform()
        #srcwin = (0, 0, self.src_ds.RasterXSize, self.src_ds.RasterYSize)
        if self.ds_open_p:            
            this_region = regions.Region().from_geo_transform(
                geoT=gt, x_count=self.src_ds.RasterXSize, y_count=self.src_ds.RasterYSize)
            try:
                zr = self.src_ds.GetRasterBand(1).ComputeRasterMinMax()
            except:
                zr = [None, None]
            this_region.zmin = zr[0]
            this_region.zmax = zr[1]
            self.infos['minmax'] = this_region.export_as_list(include_z=True)
            self.infos['numpts'] = self.src_ds.RasterXSize * self.src_ds.RasterYSize
            self.infos['wkt'] = this_region.export_as_wkt()
        self._close_ds()
        return(self.infos)

    def set_epsg(self, epsg = 4326):
        if self.ds_open_p:
            self.src_ds.SetProjection(sr_wkt(int(epsg)))

    def cut(self):
        if self.ds_open_p:
            ds_config = demfun.gather_infos()
            gt = ds_config['geoT']
            srcwin = region.srcwin(gt, ds_config['nx'], ds_config['ny'])
            
            ds_arr = self.src_ds.GetRasterBand(1).ReadAsArray(srcwin[0], srcwin[1], srcwin[2], srcwin[3])
            dst_gt = (gt[0] + (srcwin[0] * gt[1]), gt[1], 0., gt[3] + (srcwin[1] * gt[5]), 0., gt[5])
            out_ds_config = self.set_infos(srcwin[2], srcwin[3], srcwin[2] * srcwin[3], dst_gt, ds_config['proj'], ds_config['dt'], ds_config['ndv'], ds_config['fmt'])

            return(utils.gdal_write(ds_arr, dst_fn, out_ds_config))
        else: return(None, -1)
        
    def set_infos(nx, ny, nb, geoT, proj, dt, ndv, fmt):
        """set a datasource config dictionary
            
        returns gdal_config dict."""
        
        return({'nx': nx, 'ny': ny, 'nb': nb, 'geoT': geoT, 'proj': proj, 'dt': dt, 'ndv': ndv, 'fmt': fmt})
        
    def gather_infos(self, scan=False):
        """gather information from `src_ds` GDAL dataset

        Returns:
          raster_parser: self
        """

        if self.ds_open_p:
            self.gt = self.src_ds.GetGeoTransform()
            if self.region is not None and self.region.valid_p(check_xy = True):
                self.srcwin = self.region.srcwin(self.gt, self.src_ds.RasterXSize, self.src_ds.RasterYSize)
            else: self.srcwin = (0, 0, self.src_ds.RasterXSize, self.src_ds.RasterYSize)
            src_band = self.src_ds.GetRasterBand(1)
            #self.gt = (gt[0] + (self.srcwin[0]*gt[1]), gt[1], 0., gt[3] + (self.srcwin[1]*gt[5]), 0., gt[5])
            #self.gt = gt
            
            proj = self.src_ds.GetProjectionRef()
            src_srs = osr.SpatialReference()
            src_srs.ImportFromWkt(proj)
            src_srs.AutoIdentifyEPSG()
            srs_auth = src_srs.GetAuthorityCode(None)

            self.epsg = utils.int_or(srs_auth)
            self.x_count = self.srcwin[2]
            self.y_count = self.srcwin[3]
            self.dt = src_band.DataType
            self.dtn = gdal.GetDataTypeName(src_band.DataType)
            self.ndv = src_band.GetNoDataValue()
            if self.ndv is None: self.ndv = -9999
            self.fmt = self.src_ds.GetDriver().ShortName
            self.zr = None

            if scan:
                src_arr = src_band.ReadAsArray(srcwin[0], self.srcwin[1], self.srcwin[2], self.srcwin[3])
                self.zr = (np.amin(src_arr), np.amax(src_arr))
                src_arr = None
        return(self)

    def parse(self):
        
        if self.region is not None:
            self.inf()
            inf_region = regions.Region().from_list(self.infos['minmax'])
            if regions.regions_intersect_p(inf_region, self.region):
                self.data_entries = [self]
        else: self.data_entries = [self]
        return(self)
    
    def yield_xyz(self):
        """parse the data from gdal dataset src_ds (first band only)

        Yields:
          xyz: the parsed xyz data
        """

        self._open_ds()
        out_xyz = xyzfun.XYZPoint(w = 1)
        if self.src_ds is not None:
            ln = 0
            band = self.src_ds.GetRasterBand(1)
            gt = self.gt
            warp_epsg = self.warp
            
            if warp_epsg is not  None and self.epsg is not None:
                src_srs = osr.SpatialReference()
                src_srs.ImportFromEPSG(self.epsg)

                dst_srs = osr.SpatialReference()
                dst_srs.ImportFromEPSG(warp_epsg)
                try:
                    src_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
                    dst_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
                except:
                    pass
                dst_trans = osr.CoordinateTransformation(src_srs, dst_srs)
            else:
                dst_trans = None
                        
            msk_band = None
            if self.mask is not None:
                src_mask = gdal.Open(self.mask)
                msk_band = src_mask.GetRasterBand(1)

            nodata = ['{:g}'.format(-9999), 'nan', float('nan')]
            if self.ndv is not None:
                nodata.append('{:g}'.format(self.ndv))
            for y in range(self.srcwin[1], self.srcwin[1] + self.srcwin[3], 1):
                band_data = band.ReadAsArray(self.srcwin[0], y, self.srcwin[2], 1)
                if self.region is not None:
                    z_region = self.region.z_region()
                    if z_region[0] is not None:
                        band_data[band_data < z_region[0]] = -9999
                    if z_region[1] is not None:
                        band_data[band_data > z_region[1]] = -9999
                if msk_band is not None:
                   msk_data = msk_band.ReadAsArray(self.srcwin[0], y, self.srcwin[2], 1)
                   band_data[msk_data==0]=-9999
                band_data = np.reshape(band_data, (self.srcwin[2], ))
                for x_i in range(0, self.srcwin[2], 1):
                    x = x_i + self.srcwin[0]
                    z = band_data[x_i]
                    if '{:g}'.format(z) not in nodata:
                        out_xyz.x, out_xyz.y = utils._pixel2geo(x, y, gt)
                        out_xyz.z = z
                        out_xyz.w = self.weight
                        if self.region is not None and self.region.valid_p():
                            if regions.xyz_in_region_p(out_xyz, self.region):
                                ln += 1
                                if dst_trans is not None: out_xyz.transform(dst_trans)
                                yield(out_xyz)
                        else:
                            ln += 1
                            if dst_trans is not None: out_xyz.transform(dst_trans)
                            yield(out_xyz)
            band = None
            src_mask = None
            msk_band = None
            if self.verbose:
                utils.echo_msg('parsed {} data records from {}'.format(ln, self.fn))
        self._close_ds()

class GMRT(XYZDataset):
    """represnting the remote raster data from the GMRT"""
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        
        ## ==============================================
        ## GMRT URLs and directories
        ## ==============================================    
        self._gmrt_grid_url = "https://www.gmrt.org:443/services/GridServer?"
        self._gmrt_grid_urls_url = "https://www.gmrt.org:443/services/GridServer/urls?"
        self._gmrt_grid_metadata_url = "https://www.gmrt.org/services/GridServer/metadata?"
        self._outdir = os.path.join(os.getcwd(), 'gmrt')

        self.name = 'gmrt'
        self.src_gmrt = 'gmrt_{}.tif'.format(self.region.format('fn'))
        self._fn = self.fn
        #self.fn = self.src_gmrt
        
    def generate_inf(self):
        """generate a infos dictionary from the gmrt dataset

        Returns:
          dict: a data-entry infos dictionary
        """
        
        self.infos['name'] = self.fn
        self.infos['hash'] = None
        this_region = regions.Region().from_list([-179,179,-89,89])
        self.infos['minmax'] = this_region.export_as_list()
        self.infos['numpts'] = 0
        self.infos['wkt'] = this_region.export_as_wkt()
        return(self.infos)        

    def parse(self, fmt='geotiff', res='max', layer='topo'):
        """Parse the gmrt module for urls """

        if self.region is not None:
            self.inf()
            if layer != 'topo' and layer != 'topo-mask':
                layer = 'topo'
            inf_region = regions.Region().from_list(self.infos['minmax'])
            if regions.regions_intersect_p(inf_region, self.region):
                _data = {'north':self.region.ymax, 'west':self.region.xmin,
                         'south':self.region.ymin, 'east':self.region.xmax,
                         'mformat':'json', 'resolution':res, 'format':fmt}

                _req = fetch_req(self._gmrt_grid_urls_url, params=_data,
                                 tries=10, timeout=2)
                if _req is not None and _req.status_code == 200:
                    gmrt_urls = _req.json()
                    for url in gmrt_urls:
                        opts = {}
                        url_base = url.split('?')[0]
                        for url_opt in url.split('?')[1].split('&'):
                            opt_kp = url_opt.split('=')
                            opts[opt_kp[0]] = opt_kp[1]
                        opts['layer'] = layer
                        try:
                            url_enc = urllib.urlencode(opts)
                        except:
                            url_enc = urllib.parse.urlencode(opts)
                        this_url = '{}?{}'.format(url_base, url_enc)
                        url_region = regions.Region().from_list(
                            [float(opts['west']), float(opts['east']),
                             float(opts['south']), float(opts['north'])])
                        outf = 'gmrt_{}_{}.tif'.format(opts['layer'], url_region.format('fn'))
                        self.fn = this_url
                        self._fn = outf
                        self.data_format = 201
                        self.data_entries.append(self)
        else:
            utils.echo_error_msg('you must supply a region to use the gmrt datalist')
        return(self)
    
    def yield_xyz(self):
        for i, entry in enumerate(self.data_entries):
            utils.echo_msg('>{}<'.format(self.region.format('fn')))
            gmrt_ds = DatasetFactory(
                fn=self._fn, name=self.name, data_format=200, weight=self.weight, src_region=self.region,
                title=self.title, epsg=self.epsg, source=self.source, date=self.date, data_type=self.data_type,
                resolution=self.resolution, vdatum=self.vdatum, url=self.url, warp=self.warp,
                verbose=self.verbose)
            if fetch_file(entry.fn, entry._fn, verbose = entry.verbose) == 0:
                for xyz in gmrt_ds.acquire_raster_file().parse().yield_xyz():
                    yield(xyz)
            else:
                utils.echo_error_msg('failed to fetch remote file, {}...'.format(self.fn))
            utils.remove_glob(self.fn)
        
class DatasetFactory:

    data_types = {
        -1: {'name': 'datalist',
             'fmts': ['datalist', 'mb-1'],
             'class': lambda k: Datalist(**k),
             },
        167: {'name': 'yxz',
              'fmts': ['yxz'],
              'class': lambda k: XYZFile(**k),
              },
        168: {'name': 'xyz',
              'fmts': ['xyz', 'csv', 'dat', 'ascii', 'txt'],
              'class': lambda k: XYZFile(**k),
              },
        200: {'name': 'raster',
              'fmts': ['tif', 'img', 'grd', 'nc', 'vrt', 'bag'],
              'class': lambda k: RasterFile(**k),
              },
        201: {'name': 'remote-raster',
              'fmts': ['tif', 'img', 'grd', 'nc', 'vrt', 'bag'],
              'class': lambda k: RasterFile(**k),
              },
        -11: {'name': 'gmrt',
              'fmts': ['gmrt'],
              'class': lambda k: GMRT(**k),
              },
    }

    datalist_cols = ['path', 'format', 'weight', 'name', 'source',
                     'date', 'type', 'resolution', 'horz', 'vert',
                     'url']
    
    def __init__(self, fn=None, data_format=None, weight=1, epsg=4326, name="xyz_dataset", title=None,
                 source=None, date=None, data_type=None, resolution=None, vdatum=None, url=None,
                 parent=None, src_region=None, warp=None, verbose=False):
        
        self.name = name
        self.title = title
        self.data_format = data_format
        self.weight = weight
        self.source = source
        self.date = date
        self.data_type = data_type
        self.resolution = resolution
        self.epsg = epsg
        self.vdatum = vdatum
        self.url = url
        self.warp = warp
        self.region = src_region
        self.parent = parent
        self.verbose = verbose
        self.fn = fn
        self.parse_fn()
        if self.data_format is None:
            self.guess_data_format()
        
    def parse_fn(self):

        if self.fn is None: return(self)
        if os.path.exists(self.fn): return(self.fn)

        this_entry = re.findall(r'[^"\s]\S*|".+?"', self.fn.rstrip())
        try:
            entry = [x if n == 0 else utils.int_or(x) if n < 2 else utils.float_or(x) if n < 3 else x \
                     for n, x in enumerate(this_entry)]
        except Exception as e:
            utils.echo_error_msg('could not parse entry {}'.format(self.fn))
            return(self)

        if len(entry) < 2:
            for key in self.data_types.keys():
                se = entry[0].split('.')
                see = se[-1] if len(se) > 1 else entry[0].split(":")[0]
                if see in self.data_types[key]['fmts']:
                    entry.append(int(key))
                    break
            if len(entry) < 2:
                utils.echo_error_msg('could not parse entry {}'.format(self.fn))
                return(self)
            
        if len(entry) < 3:
            entry.append(self.weight)
        elif entry[2] is None:
            entry[2] = self.weight
        if len(entry) < 4:
            entry.append(self.title)
        else: self.title = entry[3]
        if len(entry) < 5:
            entry.append(self.source)
        else: self.source = entry[4]
        if len(entry) < 6:
            entry.append(self.date)
        else: self.date = entry[5]
        if len(entry) < 7:
            entry.append(self.data_type)
        else: self.data_type = entry[6]
        if len(entry) < 8:
            entry.append(self.resolution)
        else: self.resolution = entry[7]
        if len(entry) < 9:
            entry.append(self.epsg)
        else: self.epsg = entry[8]
        if len(entry) < 10:
            entry.append(self.vdatum)
        else: self.vdatum = entry[9]
        if len(entry) < 11:
            entry.append(self.url)
        else: self.url = entry[10]
        
        self.fn = entry[0] if self.parent is None \
            else os.path.join(os.path.dirname(self.parent.fn), entry[0])
        
        self.data_format = entry[1]
        if self.data_format is None:
            self.guess_data_format()
        print(entry)
        if self.weight is not None:
            self.weight *= entry[2]

        return(self)

    def guess_data_format(self):
        if self.fn is not None:
            for key in self.data_types.keys():
                if self.fn.split('.')[-1] in self.data_types[key]['fmts']:
                    self.data_format = key
                    break

    def add_data_type(self, type_def = {}):
        for key in type_def.keys():
            self.data_types[key] = type_def[key]

    def acquire_datalist(self, **kwargs):
        return(Datalist(
            fn=self.fn, data_format=self.data_format, weight=self.weight, src_region=self.region, title=self.title,
            source=self.source, date=self.date, data_type=self.data_type, resolution=self.resolution, vdatum=self.vdatum, url=self.url,
            epsg=self.epsg, warp=self.warp, name=self.name, parent=self.parent, verbose=self.verbose,
            **kwargs))

    def acquire_xyz_file(self, **kwargs):
        return(XYZFile(
            fn=self.fn, data_format=self.data_format, weight=self.weight, src_region=self.region, title=self.title,
            source=self.source, date=self.date, data_type=self.data_type, resolution=self.resolution, vdatum=self.vdatum, url=self.url,
            epsg=self.epsg, warp=self.warp, name=self.name, parent=self.parent, verbose=self.verbose,
            xpos=0, ypos=1, zpos=2, **kwargs))

    def acquire_raster_file(self, **kwargs):
        return(RasterFile(
            fn=self.fn, data_format=self.data_format, weight=self.weight, src_region=self.region, title=self.title,
            source=self.source, date=self.date, data_type=self.data_type, resolution=self.resolution, vdatum=self.vdatum, url=self.url,
            epsg=self.epsg, warp=self.warp, name=self.name, parent=self.parent, verbose=self.verbose,
            **kwargs))

    def acquire_remote_raster_file(self, **kwargs):
        return(RasterFile(
            fn=self.fn, data_format=self.data_format, weight=self.weight, src_region=self.region, title=self.title,
            source=self.source, date=self.date, data_type=self.data_type, resolution=self.resolution, vdatum=self.vdatum, url=self.url,
            epsg=self.epsg, warp=self.warp, name=self.name, parent=self.parent, verbose=self.verbose,
            remote=True, **kwargs))

    def acquire_gmrt(self, **kwargs):
        return(GMRT(
            fn=self.fn, data_format=self.data_format, weight=self.weight, src_region=self.region, title=self.title,
            source=self.source, date=self.date, data_type=self.data_type, resolution=self.resolution, vdatum=self.vdatum, url=self.url,
            epsg=self.epsg, warp=self.warp, name=self.name, parent=self.parent, verbose=self.verbose,
            **kwargs))

    def acquire(self, **kwargs):
        return(self.data_types[self.data_format]['class'](
            fn=self.fn, data_format=self.data_format, weight=self.weight, src_region=self.region, title=self.title,
            source=self.source, date=self.date, data_type=self.data_type, resolution=self.resolution, vdatum=self.vdatum, url=self.url,
            epsg=self.epsg, warp=self.warp, name=self.name, parent=self.parent, verbose=self.verbose,
            **kwargs))
    
    def acquire_dataset(self, **kwargs):
        if self.data_format == -1:
            return(self.acquire_datalist(**kwargs))#.parse())
        
        if self.data_format == 168:
            return(self.acquire_xyz_file(**kwargs))#.parse())
       
        if self.data_format == 200:
            return(self.acquire_raster_file(**kwargs))#.parse())
        
        if self.data_format == 201:
            return(self.acquire_remote_raster_file(**kwargs))#.parse())
        
        if self.data_format == -11:
            return(self.acquire_gmrt(**kwargs))#.parse())

## ==============================================
##
## datalists cli
##
## ==============================================
_datalist_fmts_long_desc = lambda: '\n  '.join(
    ['{}:\t{}'.format(key, DatasetFactory().data_types[key]['name']) for key in DatasetFactory().data_types])
_datalist_fmts_short_desc = lambda: ',  '.join(
    ['{} ({})'.format(DatasetFactory().data_types[key]['name'], key) for key in DatasetFactory().data_types])
datalists_usage = """{cmd} ({dl_version}): DataLists IMproved; Process and generate datalists

usage: {cmd} [ -ghiqwPRW [ args ] ] DATALIST ...

Options:
  -R, --region\t\tRestrict processing to the desired REGION 
\t\t\tWhere a REGION is [ xmin/xmax/ymin/ymax/[ zmin/zmax/[ wmin/wmax ] ] ]
\t\t\tUse '-' to indicate no bounding range; e.g. -R -/-/-/-/-10/10/1/-
\t\t\tOR an OGR-compatible vector file with regional polygons. 
\t\t\tIf a vector file is supplied, will use each region found therein.
\t\t\tAppend :zmin/zmax/[ wmin/wmax ] to the file path to extended REGION.
  -P, --s_epsg\t\tSet the projection EPSG code of the datalist.
  -W, --t_epsg\t\tSet the output warp projection EPSG code.
  -F, --format\t\tOnly process the given data format.

  --archive\t\tARCHIVE the datalist to the given REGION
  --glob\t\tGLOB the datasets in the current directory to stdout
  --info\t\tGenerate and return an INFO dictionary of the dataset
  --weights\t\tOutput WEIGHT values along with xyz
  --quiet\t\tLower the verbosity to a quiet

  --help\t\tPrint the usage text
  --version\t\tPrint the version information

Supported datalist formats: 
  {dl_formats}

Examples:
  % {cmd} my_data.datalist -R -90/-89/30/31
  % {cmd} -R-90/-89/30/31/-100/100 *.tif -l -w > tifs_in_region.datalist
  % {cmd} -R my_region.shp my_data.xyz -w -s_epsg 4326 -t_epsg 3565 > my_data_3565.xyz

CIRES DEM home page: <http://ciresgroups.colorado.edu/coastalDEM>\
""".format(cmd=os.path.basename(sys.argv[0]), 
           dl_version=__version__,
           dl_formats=_datalist_fmts_short_desc())

def datalists_cli(argv = sys.argv):
    """run datalists from command-line

    See `datalists_cli_usage` for full cli options.
    """

    dls = []
    epsg = None
    warp = None
    fmt = None
    i_regions = []
    these_regions = []
    want_weights = False
    want_inf = False
    want_list = False
    want_glob = False
    want_archive = False
    want_verbose = True    
    
    ## ==============================================
    ## parse command line arguments.
    ## ==============================================
    i = 1
    while i < len(argv):
        arg = argv[i]
        if arg == '--region' or arg == '-R':
            i_regions.append(str(argv[i + 1]))
            i = i + 1
        elif arg[:2] == '-R':
            i_regions.append(str(arg[2:]))
        elif arg == '-s_epsg' or arg == '--s_epsg' or arg == '-P':
            epsg = argv[i + 1]
            i = i + 1
        elif arg == '-t_epsg' or arg == '--t_epsg' or arg == '-W':
            warp = argv[i + 1]
            i = i + 1
        elif arg == '--format' or arg == '-F':
            fmt = utils.int_or(argv[i + 1])
            i = i + 1
        elif arg == '--archive' or arg == '-a':
            want_archive = True
        elif arg == '--weights' or arg == '-w':
            want_weights = True
        elif arg == '--info' or arg == '-i':
            want_inf = True
        elif arg == '--list' or arg == '-l':
            want_list = True
        elif arg == '--glob' or arg == '-g':
            want_glob = True
        elif arg == '--quiet' or arg == '-q':
            want_verbose = False
        elif arg == '--help' or arg == '-h':
            print(datalists_usage)
            sys.exit(1)
        elif arg == '--version' or arg == '-v':
            print('{}, version {}'.format(os.path.basename(sys.argv[0]), __version__))
            sys.exit(1)
        elif arg[0] == '-':
            print(datalists_usage)
            sys.exit(0)
        else: dls.append(arg)
        i = i + 1

    if want_glob:
        for key in DatasetFactory().data_types.keys():
            if key != -1:
                for f in DatasetFactory().data_types[key]['fmts']:
                    globs = glob.glob('*.{}'.format(f))
                    [sys.stdout.write('{}\n'.format(' '.join([x, str(key), '1']))) for x in globs]
        sys.exit(0)

    for i_region in i_regions:
        tmp_region = regions.Region().from_string(i_region)
        if tmp_region.valid_p():
            these_regions.append(tmp_region)
        else:
            tmp_region = regions.ogr_wkts(i_region)
            for i in tmp_region:
                if i.valid_p():
                    these_regions.append(i)
                    
    if len(these_regions) == 0:
        these_regions = [None]
    else:
        if want_verbose:
            utils.echo_msg('parsed {} region(s)'.format(len(these_regions)))

    for rn, this_region in enumerate(these_regions):
        if len(dls) == 0:
            print(datalists_usage)
            utils.echo_error_msg('you must specify some type of data')
        xdls = [DatasetFactory(
            fn=" ".join(['-' if x == "" else x for x in dl.split(":")]),
            src_region=this_region, verbose=want_verbose,
            epsg=epsg, warp=warp).acquire_dataset() for dl in dls]
        
        for xdl in xdls:            
            if xdl is not None and xdl.valid_p():
                xdl.parse()
                if not want_weights:
                    xdl.weight = None
                if want_inf:
                    print(xdl.inf())
                elif want_list:
                    xdl.echo()
                elif want_archive:
                    [x for x in xdl.archive_xyz()]
                else: xdl.dump_xyz()
### End
