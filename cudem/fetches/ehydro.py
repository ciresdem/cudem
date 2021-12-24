### ehydro.py - eHydro fetch
##
## Copyright (c) 2010 - 2021 CIRES Coastal DEM Team
##
## ehyrdo.py is part of CUDEM
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
## USACE eHydro Fetch
##
## Maintenance responsibility for more than 25,000 miles of navigation channels and 400 ports and 
## harbors throughout the United States requires extensive surveying and mapping services, including 
## boundary, topographic, hydrographic, terrestrial lidar, and multispectral and hyperspectral aerial 
## imagery collection as well as airborne topographic and bathymetric lidar acquisition, project-level 
## GIS implementation, development of file-based geodatabases, and GIS tool development.
##
## Three representative survey and mapping datasets include the National Channel Framework (NCF)—an enterprise 
## geodatabase of information on all 61 USACE-maintained high-tonnage channels —hydrographic surveys, which 
## provide assistance in locating navigable channels, determining dredging requirements, verifying dredging 
## accuracy, and maintaining harbors and rivers —and Inland Electronic Navigational Charts(IENC), accurate 
## navigational charts provided in a highly structured data format for use in navigation systems and to increase 
## overall navigational safety..
##
## Fetch USACE bathymetric surveys via eHydro
##
## Fields:
##
##     objectid (type: esriFieldTypeOID, alias: objectid, SQL Type: sqlTypeOther, length: 0, nullable: false, editable: false)
##     surveyjobidpk (type: esriFieldTypeString, alias: SURVEYJOBIDPK, SQL Type: sqlTypeOther, length: 100, nullable: true, editable: true)
##     sdsid (type: esriFieldTypeString, alias: SDSID, SQL Type: sqlTypeOther, length: 40, nullable: true, editable: true)
##     sdsfeaturename (type: esriFieldTypeString, alias: SDSFEATURENAME, SQL Type: sqlTypeOther, length: 80, nullable: true, editable: true)
##     sdsmetadataid (type: esriFieldTypeString, alias: SDSMETADATAID, SQL Type: sqlTypeOther, length: 80, nullable: true, editable: true)
##     surveytype (type: esriFieldTypeString, alias: SURVEYTYPE, SQL Type: sqlTypeOther, length: 26, nullable: true, editable: true)
##     channelareaidfk (type: esriFieldTypeString, alias: CHANNELAREAIDFK, SQL Type: sqlTypeOther, length: 50, nullable: true, editable: true)
##     dateuploaded (type: esriFieldTypeDate, alias: dateUploaded, SQL Type: sqlTypeOther, length: 8, nullable: true, editable: true)
##     usacedistrictcode (type: esriFieldTypeString, alias: usaceDistrictCode, SQL Type: sqlTypeOther, length: 20, nullable: true, editable: true)
##     surveydatestart (type: esriFieldTypeDate, alias: SURVEYDATESTART, SQL Type: sqlTypeOther, length: 8, nullable: true, editable: true)
##     surveydateend (type: esriFieldTypeDate, alias: SURVEYDATEEND, SQL Type: sqlTypeOther, length: 8, nullable: true, editable: true)
##     sourcedatalocation (type: esriFieldTypeString, alias: SOURCEDATALOCATION, SQL Type: sqlTypeOther, length: 255, nullable: true, editable: true)
##     sourceprojection (type: esriFieldTypeString, alias: SOURCEPROJECTION, SQL Type: sqlTypeOther, length: 75, nullable: true, editable: true)
##     mediaidfk (type: esriFieldTypeString, alias: MEDIAIDFK, SQL Type: sqlTypeOther, length: 100, nullable: true, editable: true)
##     projectedarea (type: esriFieldTypeDouble, alias: PROJECTEDAREA, SQL Type: sqlTypeOther, nullable: true, editable: true)
##     sdsfeaturedescription (type: esriFieldTypeString, alias: SDSFEATUREDESCRIPTION, SQL Type: sqlTypeOther, length: 255, nullable: true, editable: true)
##     dateloadedenterprise (type: esriFieldTypeDate, alias: dateLoadedEnterprise, SQL Type: sqlTypeOther, length: 8, nullable: true, editable: true)
##     datenotified (type: esriFieldTypeDate, alias: dateNotified, SQL Type: sqlTypeOther, length: 8, nullable: true, editable: true)
##     sourcedatacontent (type: esriFieldTypeString, alias: sourceDataContent, SQL Type: sqlTypeOther, length: 1000, nullable: true, editable: true)
##     plotsheetlocation (type: esriFieldTypeString, alias: PLOTSHEETLOCATION, SQL Type: sqlTypeOther, length: 255, nullable: true, editable: true)
##     sourceagency (type: esriFieldTypeString, alias: SOURCEAGENCY, SQL Type: sqlTypeOther, length: 20, nullable: true, editable: true)
##     globalid (type: esriFieldTypeGlobalID, alias: GlobalID, SQL Type: sqlTypeOther, length: 38, nullable: false, editable: false)
##     Shape__Area (type: esriFieldTypeDouble, alias: Shape__Area, SQL Type: sqlTypeDouble, nullable: true, editable: false)
##     Shape__Length (type: esriFieldTypeDouble, alias: Shape__Length, SQL Type: sqlTypeDouble, nullable: true, editable: false)
##
### Code:

import os
import sys
import lxml.etree

from osgeo import ogr

from cudem import utils
from cudem import regions
from cudem import datasets

import cudem.fetches.utils as f_utils
import cudem.fetches.FRED as FRED

class eHydro(f_utils.FetchModule):
    """Fetch USACE bathymetric surveys via eHydro

specify `inc` to blockmedian the data when processing
"""
    
    def __init__(self, where='1=1', inc=None, **kwargs):
        super().__init__(**kwargs)
        self._ehydro_gj_api_url = 'https://opendata.arcgis.com/datasets/80a394bae6b547f1b5788074261e11f1_0.geojson'
        self._ehydro_api_url = 'https://services7.arcgis.com/n1YM8pTrFmm7L4hs/arcgis/rest/services/eHydro_Survey_Data/FeatureServer/0'
        self._ehydro_query_url = '{0}/query?'.format(self._ehydro_api_url)
        self._outdir = os.path.join(os.getcwd(), 'ehydro')
        self.name = 'ehydro'
        self.where = where
        self.inc = utils.str2inc(inc)

    def run(self):
        '''Run the eHydro fetching module'''
        
        if self.region is None:
            return([])
        
        _data = {
            'where': self.where,
            'outFields': '*',
            'geometry': self.region.format('bbox'),
            'inSR':4326,
            'outSR':4326,
            'f':'pjson',
            'returnGeometry':'False',
        }
        _req = f_utils.Fetch(self._ehydro_query_url).fetch_req(params=_data)
        if _req is not None:
            features = _req.json()
            for feature in features['features']:
                fetch_fn = feature['attributes']['sourcedatalocation']
                self.results.append([fetch_fn, fetch_fn.split('/')[-1], 'ehydro'])
                
        return(self)

    def yield_xyz(self, entry):
        src_zip = os.path.basename(entry[1])
        src_epsg = None
        src_region = None
        
        if f_utils.Fetch(
                entry[0], callback=self.callback, verbose=self.verbose
        ).fetch_file(src_zip) == 0:
            src_xmls = utils.p_unzip(src_zip, ['xml', 'XML'])
            for src_xml in src_xmls:
                if src_region is None:
                    this_xml = lxml.etree.parse(src_xml)
                    if this_xml is not None:
                        try:
                            w = this_xml.find('.//westbc').text
                            e = this_xml.find('.//eastbc').text
                            n = this_xml.find('.//northbc').text
                            s = this_xml.find('.//southbc').text
                            src_region = regions.Region().from_list([float(w), float(e), float(s), float(n)])
                        except: pass
                            
                utils.remove_glob(src_xml)

            if src_region is None:
                sys.exit()

            if src_epsg is None:
                this_geom = src_region.export_as_geom()
                sp_fn = os.path.join(FRED.fetchdata, 'stateplane.geojson')
                sp = ogr.Open(sp_fn)
                layer = sp.GetLayer()
                
                for feature in layer:
                    geom = feature.GetGeometryRef()
                    if this_geom.Intersects(geom):
                        src_epsg = feature.GetField('EPSG')
                        break
                    
                sp = None

            src_usaces = utils.p_unzip(src_zip, ['XYZ', 'xyz', 'dat'])
            for src_usace in src_usaces:
                _dl = datasets.XYZFile(
                    fn=src_usace,
                    data_format=168,
                    x_scale=.3048,
                    y_scale=.3048,
                    z_scale=-.3048,
                    src_srs='epsg:{}'.format(src_epsg),
                    dst_srs=self.dst_srs,
                    src_region=src_region,
                    name=src_usace,
                    verbose=self.verbose,
                    remote=True
                )
                for xyz in _dl.yield_xyz():
                    yield(xyz)
                    
                utils.remove_glob(src_usace, src_usace+'.inf')
                
        else:
            utils.echo_error_msg('failed to fetch remote file, {}...'.format(entry[0]))
            
        utils.remove_glob(src_zip)

    def yield_xyz_test(self, entry):
        src_zip = os.path.basename(entry[1])
        src_region = None
        this_proj = None
        
        if f_utils.Fetch(
                entry[0], callback=self.callback, verbose=self.verbose
        ).fetch_file(src_zip) == 0:

            _data = {
                'where': "sourcedatalocation='{}'".format(entry[0]),
                'outFields': '*',
                'geometry': self.region.format('bbox'),
                'inSR':4326,
                'outSR':4326,
                'f':'pjson',
                'returnGeometry':'False',
            }
            _req = f_utils.Fetch(self._ehydro_query_url).fetch_req(params=_data)
            if _req is not None:
                print(_req.text)
                features = _req.json()
                for feature in features['features']:
                    this_proj = feature['attributes']['sourceprojection']
                    print(this_proj)
                    
            src_usaces = utils.p_unzip(src_zip, ['XYZ', 'xyz', 'dat'])
            for src_usace in src_usaces:
                _dl = datasets.XYZFile(
                    fn=src_usace,
                    data_format=168,
                    #x_scale=.3048,
                    #y_scale=.3048,
                    z_scale=-.3048,
                    src_srs=this_proj,
                    dst_srs=self.dst_srs,
                    src_region=src_region,
                    name=src_usace,
                    verbose=self.verbose,
                    remote=True
                )
                for xyz in _dl.yield_xyz():
                    yield(xyz)
                    
                utils.remove_glob(src_usace, src_usace+'.inf')
                
        else:
            utils.echo_error_msg('failed to fetch remote file, {}...'.format(entry[0]))
            
        utils.remove_glob(src_zip)
        
### End
