### FRED.py
##
## Copyright (c) 2010 - 2023 CIRES Coastal DEM Team
##
## FRED.py is part of CUDEM
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
## Fetches Remote Elevation Datalist (FRED)
##
## the fetches reference vector location and related functions
## for generating and parsing FRED
##
### Code:

import os
import json

from osgeo import ogr

from cudem import utils
from cudem import regions

import cudem.fetches.utils as f_utils

this_dir, this_filename = os.path.split(__file__)
fetchdata = os.path.join(this_dir, 'data')

class FRED:
    def __init__(self, name='FRED', verbose=False, local=False):
        self._verbose = verbose
        self.fetchdata = os.path.join(this_dir, 'data')
        self.driver = ogr.GetDriverByName('GeoJSON')
        self.fetch_v = '{}.geojson'.format(name)
        if local:
            self.FREDloc = self.fetch_v
        elif os.path.exists(self.fetch_v):
            self.FREDloc = self.fetch_v
        elif os.path.exists(os.path.join(self.fetchdata, self.fetch_v)):
            self.FREDloc = os.path.join(self.fetchdata, self.fetch_v)
        else: self.FREDloc = self.fetch_v
        if self._verbose: utils.echo_msg('using {}'.format(self.FREDloc))
        self.ds = None
        self.layer = None
        self.open_p = False
        
        self._fields = ['Name', 'ID', 'Date', 'Agency', 'MetadataLink',
                        'MetadataDate', 'DataLink', 'IndexLink', 'Link',
                        'DataType', 'DataSource', 'Resolution', 'HorizontalDatum',
                        'VerticalDatum', 'LastUpdate', 'Etcetra', 'Info']
        
    def _create_ds(self):
        utils.remove_glob(self.FREDloc)
        self.ds = self.driver.CreateDataSource(self.FREDloc)        
        self.layer = self.ds.CreateLayer('FRED', None, ogr.wkbMultiPolygon)
        ldfn = self.layer.GetLayerDefn()
        self.layer.CreateField(ogr.FieldDefn('Name', ogr.OFTString))
        self.layer.CreateField(ogr.FieldDefn('ID', ogr.OFTString))
        self.layer.CreateField(ogr.FieldDefn('Date', ogr.OFTInteger))
        self.layer.CreateField(ogr.FieldDefn('Agency', ogr.OFTString))
        self.layer.CreateField(ogr.FieldDefn('MetadataLink', ogr.OFTString))
        self.layer.CreateField(ogr.FieldDefn('MetadataDate', ogr.OFTString))
        self.layer.CreateField(ogr.FieldDefn('DataLink', ogr.OFTString))
        self.layer.CreateField(ogr.FieldDefn('IndexLink', ogr.OFTString))
        self.layer.CreateField(ogr.FieldDefn('Link', ogr.OFTString))
        self.layer.CreateField(ogr.FieldDefn('DataType', ogr.OFTString))
        self.layer.CreateField(ogr.FieldDefn('DataSource', ogr.OFTString))
        self.layer.CreateField(ogr.FieldDefn('Resolution', ogr.OFTString))
        self.layer.CreateField(ogr.FieldDefn('HorizontalDatum', ogr.OFTString))
        self.layer.CreateField(ogr.FieldDefn('VerticalDatum', ogr.OFTString))
        self.layer.CreateField(ogr.FieldDefn('LastUpdate', ogr.OFTString))
        self.layer.CreateField(ogr.FieldDefn('Etcetra', ogr.OFTString))
        self.layer.CreateField(ogr.FieldDefn('Info', ogr.OFTString))

    def _add_survey(self, Name = 'fetches', ID = None, Date = None, Agency = None, MetadataLink = None,
                    MetadataDate = None, DataLink = None, IndexLink = None, Link = None, DataType = None,
                    DataSource = None, Resolution = None, HorizontalDatum = None, VerticalDatum = None,
                    Lastupdate = None, Etcetra = None, Info = None, geom = None):
        if not self.open_p: return(None)
        if geom is not None and self.layer is not None:
            self._add_feature([
                {'Name': Name, 'ID': ID, 'Agency': Agency,
                 'Date': Date, 'MetadataLink': MetadataLink,
                 'MetadataDate': MetadataDate, 'DataLink': DataLink,
                 'IndexLink': IndexLink, 'Link': Link, 'DataType': DataType,
                 'DataSource': DataSource, 'Resolution': Resolution,
                 'HorizontalDatum': HorizontalDatum, 'VerticalDatum': VerticalDatum,
                 'LastUpdate': utils.this_date(), 'Etcetra': Etcetra, 'Info': Info},
                geom.ExportToJson()
            ])
            [self.layer.SetFeature(ff) for ff in self.layer]
        else: return(None)
        
    def _open_ds(self, mode = 0):
        if not self.open_p:
            try:
                self.ds = self.driver.Open(self.FREDloc, mode)
            except: self.ds = None

            if self.ds is None or len(self.ds.GetLayer()) == 0:
                self.ds = None
                self._create_ds()
                
            self.layer = self.ds.GetLayer()
            self.open_p = True
            return(0)
        else: return(-1)
        
    def _close_ds(self):
        self.layer = self.ds = None
        self.open_p = False
        return(0)
        
    def _get_fields(self):
        if self.open_p:
            schema = []
            ldefn = self.layer.GetLayerDefn()
            for n in range(ldefn.GetFieldCount()):
                fdefn = ldefn.GetFieldDefn(n)
                schema.append(fdefn.name)
            return(schema)
        else: return(-1)
        
    def _add_feature(self, survey):
        '''add a survey to the reference vector layer'''
        if self.open_p:
            layer_defn = self.layer.GetLayerDefn()
            feat = ogr.Feature(layer_defn)
            geom = ogr.CreateGeometryFromJson(survey[1])
            geom_valid = geom.MakeValid()
            feat.SetGeometry(geom_valid)
            for field in self._fields:
                try:
                    feat.SetField(field, survey[0][field])
                except: feat.SetField(field, -1)
            self.layer.CreateFeature(feat)
            feat = None
            return(0)
        else: return(-1)

    def _edit_feature(self, feature, survey):
        if self.open_p:
            geom = ogr.CreateGeometryFromJson(survey[1])
            feature.SetGeometry(geom)
            for field in self._fields:
                try:
                    feature.SetField(field, survey[0][field])
                except: feature.SetField(field, -1)
            self.layer.SetFeature(feature)
            return(0)
        else: return(-1)

    def _add_surveys(self, surveys):
        '''update or create a reference vector using a list of surveys'''
        if self.open_p:
            if self.layer is not None:
                for survey in surveys:
                    self._add_survey(**survey)
            return(0)
        else: return(-1)

    def _get_region(self, where = [], layers = []):
        out_regions = []
        if self._verbose: _prog = _progress('gathering regions from {}...'.format(self.FREDloc))
        for i, layer in enumerate(layers):
            if self._verbose: _prog.update_perc((i, len(layers)))
            this_layer = self.layer
            this_layer.SetAttributeFilter("DataSource = '{}'".format(layer))
            [this_layer.SetAttributeFilter('{}'.format(filt)) for filt in where]
            for feat in this_layer:
                geom = feat.GetGeometryRef()
                wkt = geom.ExportToWkt()
                this_region = ogr.CreateGeometryFromWkt(wkt).GetEnvelope()
                if len(out_regions) > 0:
                    out_regions = regions_merge(out_regions, this_region)
                else: out_regions = this_region
        return(out_regions)

    def _attribute_filter(self, where=[]):
        if self.open_p:
            attr_str = []
            [attr_str.append(f) for f in where]
            wf = ' AND '.join(attr_str)
            self.layer.SetAttributeFilter(wf)
            return(0)
        else: return(-1)
        
    def _filter(self, region=None, where=[], layers=[]):
        """Search for data in the reference vector file"""

        _results = []
        if region is not None:
            _boundsGeom = region.export_as_geom()
        else:
            _boundsGeom = None

        if not self.open_p:
            self._open_ds()
            close_p = True
        else:
            close_p = False

        with utils.CliProgress(
                total=len(layers),
                message='filtering {}'.format(self.FREDloc),
                verbose=self._verbose
        ) as pbar:
            for i, layer in enumerate(layers):
                pbar.update(1)
                #this_layer = self.layer
                where.append("DataSource = '{}'".format(layer))
                if self._verbose:
                    utils.echo_msg('FRED region: {}'.format(region))
                    utils.echo_msg('FRED filter: {}'.format(where))

                self._attribute_filter(where = where)
                for feat in self.layer:
                    if _boundsGeom is not None:
                        geom = feat.GetGeometryRef()
                        if geom is not None:
                            if _boundsGeom.Intersects(geom):
                                _results.append({})
                                f_j = json.loads(feat.ExportToJson())
                                for key in f_j['properties'].keys():
                                    _results[-1][key] = feat.GetField(key)
                    else:
                        _results.append({})
                        f_j = json.loads(feat.ExportToJson())
                        for key in f_j['properties'].keys():
                            _results[-1][key] = feat.GetField(key)
                        
            #this_layer = None
        if close_p:
            self._close_ds()
            
        #clear where
        #where = []
            
        return(_results)

## ==============================================
## lambdas for the FRED using the module object `mod`
## ==============================================
_filter_FRED = lambda mod: mod.FRED._filter(region=mod.wgs_region, where=mod.where, layers=[mod.name])
_update_FRED = lambda mod, s: mod.FRED._add_surveys(s)
_filter_FRED_index = lambda mod: [utils.echo_msg(json.dumps(f, indent = 2)) for f in _filter_FRED(mod)]

### End
