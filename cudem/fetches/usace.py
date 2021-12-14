### usace.py - USACE fetch
##
## Copyright (c) 2010 - 2021 CIRES Coastal DEM Team
##
## usace.py is part of CUDEM
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
## USACE Fetch
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

class USACE(f_utils.FetchModule):
    """Fetch USACE bathymetric surveys"""
    
    def __init__(self, s_type=None, inc=None, **kwargs):
        super().__init__(**kwargs)
        self._usace_gj_api_url = 'https://opendata.arcgis.com/datasets/80a394bae6b547f1b5788074261e11f1_0.geojson'
        self._usace_gs_api_url = 'https://services7.arcgis.com/n1YM8pTrFmm7L4hs/arcgis/rest/services/eHydro_Survey_Data/FeatureServer/0/query?outFields=*&where=1%3D1'
        self._outdir = os.path.join(os.getcwd(), 'usace')
        self.name = 'usace'
        self.s_type = s_type
        self.inc = utils.str2inc(inc)

    def run(self):
        '''Run the USACE fetching module'''
        
        if self.region is None:
            return([])
        
        _data = {
            'geometry': self.region.format('bbox'),
            'inSR':4326,
            'outSR':4326,
            'f':'pjson',
        }
        _req = f_utils.Fetch(self._usace_gs_api_url).fetch_req(params=_data)
        if _req is not None:
            survey_list = _req.json()
            for feature in survey_list['features']:
                fetch_fn = feature['attributes']['sourcedatalocation']
                if self.s_type is not None:
                    if feature['attributes']['surveytype'].lower() == self.s_type.lower():
                        self.results.append([fetch_fn, fetch_fn.split('/')[-1], 'usace'])
                        
                else:
                    self.results.append([fetch_fn, fetch_fn.split('/')[-1], 'usace'])
                    
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
### End
