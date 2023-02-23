### bluetopo.py - NOAA BlueTOPO
##
## Copyright (c) 2022, 2023 CIRES Coastal DEM Team
##
## bluetopo is part of CUDEM
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
## Output 'tiff' files are 3 bands
## 1 - Elevation
## 2 - Uncertainty
## 3 - Data Source Table
##
## yield_xyz outputs elevation (band 1)
## elevation data is in NAVD88
##
## https://nauticalcharts.noaa.gov/data/bluetopo_specs.html
## https://noaa-ocs-nationalbathymetry-pds.s3.amazonaws.com/index.html#
## https://www.nauticalcharts.noaa.gov/data/bluetopo.html
##
### Code:

import os
import json

from osgeo import ogr

from cudem import utils
from cudem import regions
from cudem import datasets

import cudem.fetches.utils as f_utils
import cudem.fetches.FRED as FRED

## boto3 for aws api
import boto3

class BlueTopo(f_utils.FetchModule):
    """BlueTOPO"""
    
    def __init__(self, **kwargs):
        super().__init__(name='bluetopo', **kwargs)
        self._bt_bucket = 'noaa-ocs-nationalbathymetry-pds'
        self._bluetopo_base_url = 'https://noaa-ocs-nationalbathymetry-pds.s3.amazonaws.com/index.html#BlueTopo/'

        s3 = boto3.client('s3', aws_access_key_id='', aws_secret_access_key='')
        s3._request_signer.sign = (lambda *args, **kwargs: None)

        r = s3.list_objects(Bucket = self._bt_bucket, Prefix='BlueTopo/_BlueTopo_Tile_Scheme')
        
        self._bluetopo_index_url = 'https://noaa-ocs-nationalbathymetry-pds.s3.amazonaws.com/{}'.format(r['Contents'][0]['Key'])
        
    def run(self):
        s3 = boto3.client('s3', aws_access_key_id='', aws_secret_access_key='')
        s3._request_signer.sign = (lambda *args, **kwargs: None)
        status = f_utils.Fetch(self._bluetopo_index_url, verbose=self.verbose).fetch_file('bluetopo.gpkg')
        try:
            v_ds = ogr.Open('bluetopo.gpkg')
        except:
            v_ds = None
            status = -1
            
        if v_ds is not None:
            layer = v_ds.GetLayer()
            _boundsGeom = self.region.export_as_geom()
            layer.SetSpatialFilter(_boundsGeom)            
            fcount = layer.GetFeatureCount()
            for feature in layer:
                if feature is None:
                    continue
                
                tile_name = feature.GetField('tile')
                r = s3.list_objects(Bucket = 'noaa-ocs-nationalbathymetry-pds', Prefix='BlueTopo/{}'.format(tile_name))
                if 'Contents' in r:
                    for key in r['Contents']:
                        if key['Key'].split('.')[-1] == 'tiff':
                            data_link = 'https://noaa-ocs-nationalbathymetry-pds.s3.amazonaws.com/{}'.format(key['Key'])
                            self.results.append(
                                [data_link,
                                 os.path.join(self._outdir, data_link.split('/')[-1]),
                                 'raster']
                            )
            v_ds = None
            
        utils.remove_glob('bluetopo.gpkg')
        return(self)

    def yield_ds(self, entry):
        if f_utils.Fetch(entry[0], callback=self.callback, verbose=self.verbose).fetch_file(entry[1]) == 0:
            _ds = datasets.RasterFile(
                fn=entry[1],
                data_format=200,
                dst_srs=self.dst_srs,
                src_region=self.region,
                x_inc=self.x_inc,
                y_inc=self.y_inc,
                verbose=self.verbose
            )
            yield(_ds)

    def yield_xyz(self, entry):
        for _ds in self.yield_ds(entry):
            for xyz in _ds.yield_xyz():
                yield(xyz)

    def yield_xyz(self, entry):
        for _ds in self.yield_ds(entry):
            for arr in _ds.yield_array():
                yield(arr)


class BlueTopo_FRED(f_utils.FetchModule):
    """BlueTOPO"""
    
    def __init__(self, where='1=1', layer=0, **kwargs):
        super().__init__(name='bluetopo', **kwargs)
        self._bluetopo_base_url = 'https://noaa-ocs-nationalbathymetry-pds.s3.amazonaws.com/index.html#BlueTopo/'
        #self._bluetopo_index_url = 'https://noaa-ocs-nationalbathymetry-pds.s3.amazonaws.com/BlueTopo/BlueTopo-Tile-Scheme/BlueTopo_Tile_Scheme_20211214.gpkg'
        self._bluetopo_index_url = 'https://noaa-ocs-nationalbathymetry-pds.s3.amazonaws.com/BlueTopo/BlueTopo-Tile-Scheme/BlueTopo_Tile_Scheme_20230114.gpkg'
        self.where = [where] if len(where) > 0 else []
        
        self.FRED = FRED.FRED(name=self.name, verbose = self.verbose)
        self.update_if_not_in_FRED()
        
    def update_if_not_in_FRED(self):
        self.FRED._open_ds()
        self.FRED._attribute_filter(["DataSource = '{}'".format(self.name)])
        if len(self.FRED.layer) == 0:
            self.FRED._close_ds()
            self.update()
        self.FRED._close_ds()
        
    def update(self):
        self.FRED._open_ds()        
        status = f_utils.Fetch(self._bluetopo_index_url, verbose=self.verbose).fetch_file('bluetopo.gpkg')
        
        shp_regions = regions.gdal_ogr_regions('bluetopo.gpkg')
        shp_region = regions.Region()
        for this_region in shp_regions:
            if shp_region.valid_p(check_xy=True):
                shp_region = regions.regions_merge(shp_region, this_region)
            else: shp_region = this_region
        geom = shp_region.export_as_geom()
        
        self.FRED._attribute_filter(["ID = '{}'".format('BLUETOPO-1')])
        if self.FRED.layer is None or len(self.FRED.layer) == 0:
            self.FRED._add_survey(Name = 'BlueTopo', ID = 'BLUETOPO-1', Agency = 'NOAA', Date = utils.this_year(),
                                  MetadataLink = self._bluetopo_index_url, MetadataDate = utils.this_year(),
                                  DataLink = self._bluetopo_index_url, IndexLink = self._bluetopo_index_url,
                                  DataType = 'raster', DataSource = 'bluetopo', Info = 'Bathy Only', geom = geom)
        utils.remove_glob('bluetopo.gpkg')
        self.FRED._close_ds()

    def run(self):
        s3 = boto3.client('s3', aws_access_key_id='', aws_secret_access_key='')
        s3._request_signer.sign = (lambda *args, **kwargs: None)
        for surv in FRED._filter_FRED(self):
            data_link = None
            status = f_utils.Fetch(self._bluetopo_index_url, verbose=self.verbose).fetch_file('bluetopo.gpkg')
            try:
                v_ds = ogr.Open('bluetopo.gpkg')
            except:
                v_ds = None
                status = -1
            if v_ds is not None:
                layer = v_ds.GetLayer()
                try:
                    self.FRED.layer.SetAttributeFilter("Name = '{}'".format(name))
                except: pass                
                fcount = layer.GetFeatureCount()
                for f in range(0, fcount):
                    feature = layer[f]
                    if feature is None: continue
                    geom = feature.GetGeometryRef()
                    if geom.Intersects(self.region.export_as_geom()):
                        tile_name = feature.GetField('TileName')
                        r = s3.list_objects(Bucket = 'noaa-ocs-nationalbathymetry-pds', Prefix='BlueTopo/{}'.format(tile_name))
                        #print(r)
                        if 'Contents' in r:
                            for key in r['Contents']:
                                print(key)
                                if key['Key'].split('.')[-1] == 'tiff':
                                    data_link = 'https://noaa-ocs-nationalbathymetry-pds.s3.amazonaws.com/{}'.format(key['Key'])
                        if data_link is not None:
                            self.results.append([data_link, os.path.join(self._outdir, data_link.split('/')[-1]), surv['DataType']])
                v_ds = None
            utils.remove_glob('bluetopo.gpkg')
        
        return(self)

    def yield_xyz(self, entry):
        src_dc = os.path.basename(entry[1])
        if f_utils.Fetch(entry[0], callback=self.callback, verbose=self.verbose).fetch_file(src_dc) == 0:
            _ds = datasets.RasterFile(
                fn=src_dc,
                data_format=200,
                dst_srs=self.dst_srs,
                src_srs=None,
                src_region=self.region,
                x_inc=self.x_inc,
                y_inc=self.y_inc,
                verbose=self.verbose
            )
            for xyz in _ds.yield_xyz():
                yield(xyz)

            utils.remove_glob(src_dc)
    
### End
