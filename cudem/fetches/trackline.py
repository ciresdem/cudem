### trackline.py - NOS fetch
##
## Copyright (c) 2010 - 2021 CIRES Coastal DEM Team
##
## trackline.py is part of CUDEM
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
## TRACKLINE Fetch
##
### Code:

import os

from cudem import utils
from cudem import regions
from cudem import datasets

import cudem.fetches.utils as f_utils

class Trackline(f_utils.FetchModule):
    """trackline"""
    
    def __init__(self, where='1=1', **kwargs):
        super().__init__(**kwargs)
        #self._trackline_url = 'https://gis.ngdc.noaa.gov/arcgis/rest/services/web_mercator/trackline_bathymetry/MapServer/0'
        self._trackline_url = 'https://gis.ngdc.noaa.gov/arcgis/rest/services/web_mercator/trackline_combined_dynamic/MapServer/1'
        self._trackline_query_url = '{0}/query?'.format(self._trackline_url)
        self._outdir = os.path.join(os.getcwd(), 'trackline')
        self.name = 'trackline'
        self.where = where
        
    def run(self):
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
        _req = f_utils.Fetch(self._trackline_query_url, verbose=self.verbose).fetch_req(params=_data)
        if _req is not None:
            features = _req.json()
            ids = []
            for feature in features['features']:
                ids.append(feature['attributes']['SURVEY_ID'])

            print('http://www.ngdc.noaa.gov/trackline/request/?surveyIds={}'.format(','.join(ids)))
### End
