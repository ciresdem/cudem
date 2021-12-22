### dav.py - DAV fetch
##
## Copyright (c) 2010 - 2021 CIRES Coastal DEM Team
##
## dav.py is part of CUDEM
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
## Fetch data from the Digital Coast
##
### Code:

import os
import json

from osgeo import ogr

from cudem import utils
from cudem import regions
from cudem import datasets
from cudem import xyzfun

import cudem.fetches.utils as f_utils

class DAV(f_utils.FetchModule):
    """Fetch NOAA lidar data from DAV"""
    
    def __init__(self, where='1=1', **kwargs):
        super().__init__(**kwargs)
        self._dav_api_url = 'https://maps.coast.noaa.gov/arcgis/rest/services/DAV/ElevationFootprints/MapServer/0/query?'
        self._outdir = os.path.join(os.getcwd(), 'dav')
        self.name = 'dav'
        self.where = where
        
    def run(self):
        '''Run the DAV fetching module'''
        
        if self.region is None:
            return([])
        
        _data = {
            'where': self.where,
            'outFields': '*',
            #'outFields': 'externalproviderlink',
            'geometry': self.region.format('bbox'),
            'inSR':4326,
            'outSR':4326,
            'f':'pjson',
            'returnGeometry':'False',
        }
        _req = f_utils.Fetch(self._dav_api_url, verbose=self.verbose).fetch_req(params=_data)
        if _req is not None:
            features = _req.json()
            for feature in features['features']:
                links = json.loads(feature['attributes']['ExternalProviderLink'])
                for link in links['links']:
                    if link['serviceID'] == 46:
                        urllist = 'urllist' + str(feature['attributes']['ID']) + '.txt'
                        index_zipfile = 'tileindex.zip'
                        index_zipurl = link['link'] + '/' + index_zipfile
                        if f_utils.Fetch(link['link'] + '/' + urllist, verbose=True).fetch_file(urllist) == 0:
                            with open(urllist, 'r') as ul:
                                for line in ul:
                                    if 'tileindex' in line:
                                        index_zipurl = line.strip()
                                        break
                                    
                            utils.remove_glob(urllist)

                        if f_utils.Fetch(
                                index_zipurl, callback=self.callback, verbose=self.verbose
                        ).fetch_file(index_zipfile) == 0:
                            index_shps = utils.p_unzip('tileindex.zip', ['shp', 'shx', 'dbf', 'prj'])
                            index_shp = None
                            for v in index_shps:
                                if v.split('.')[-1] == 'shp':
                                    index_shp = v
                                    
                            index_ds = ogr.Open(index_shp)
                            index_layer = index_ds.GetLayer(0)
                            for index_feature in index_layer:
                                index_geom = index_feature.GetGeometryRef()
                                if index_geom.Intersects(self.region.export_as_geom()):
                                    tile_url = index_feature.GetField('URL').strip()
                                    self.results.append([tile_url, '{}/{}'.format(
                                        feature['attributes']['ID'], tile_url.split('/')[-1]
                                    ), feature['attributes']['DataType']])
                                    
                            index_ds = index_layer = None
                            utils.remove_glob(index_zipfile, *index_shps)

        return(self)

    def yield_xyz(self, entry):
        src_dc = os.path.basename(entry[1])
        src_ext = src_dc.split('.')[-1].lower()
        if src_ext == 'laz' or src_ext == 'las': dt = 'lidar'
        elif src_ext == 'tif' or src_ext == 'img': dt = 'raster'
        else: dt = None
        if dt == 'lidar':
            if f_utils.Fetch(entry[0], callback=self.callback, verbose=self.verbose).fetch_file(src_dc) == 0:
                _ds = datasets.LASFile(
                    fn=src_dc,
                    data_format=400,
                    dst_srs=self.dst_srs,
                    name=src_dc,
                    src_region=self.region,
                    verbose=self.verbose,
                    remote=True
                )
                if self.inc is not None:
                    b_region = regions.regions_reduce(self.region, regions.Region().from_list(_ds.infos['minmax']))
                    xyz_func = lambda p: _ds.dump_xyz(dst_port=p, encode=True)
                    for xyz in utils.yield_cmd(
                            'gmt blockmedian -I{:.10f} {} -r -V'.format(self.inc, b_region.format('gmt')),
                            verbose=self.verbose,
                            data_fun=xyz_func
                    ):
                        yield(xyzfun.XYZPoint().from_list([float(x) for x in xyz.split()]))
                        
                else:
                    for xyz in _ds.yield_xyz():
                        yield(xyz)                        
        elif dt == 'raster':
            if f_utils.Fetch(entry[0], callback=self.callback, verbose=self.verbose).fetch_file(src_dc) == 0:
                _ds = datasets.RasterFile(
                    fn=src_dc,
                    data_format=200,
                    dst_srs=self.dst_srs,
                    src_srs=None,
                    name=src_dc,
                    src_region=self.region,
                    verbose=self.verbose
                )
                for xyz in _ds.block_xyz(inc=self.inc, want_gmt=True) if self.inc is not None else _ds.yield_xyz():
                    yield(xyz)
                utils.remove_glob(src_dc)

### End
