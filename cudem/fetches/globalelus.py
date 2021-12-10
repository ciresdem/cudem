### globalelus.py - USACE fetch
##
## Copyright (c) 2010 - 2021 CIRES Coastal DEM Team
##
## globalelus.py is part of CUDEM
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
## GLOBAL ELUS Fetch (Landsat coast)
##
### Code:

import os

from cudem import utils
from cudem import regions
from cudem import datasets

import cudem.fetches.utils as f_utils

class GlobalELUS(f_utils.FetchModule):
    '''Fetch GLOBAL ELUS data from USGS'''
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        
        self._usace_gs_api_url = 'https://rmgsc.cr.usgs.gov/arcgis/rest/services/globalelus/MapServer/export?'
        
        self._outdir = os.path.join(os.getcwd(), 'globalelus')
        self.name = 'globalelus'

    def run(self):
        '''Run the GLOBALELUS fetching module'''
        
        if self.region is None: return([])
        _data = {
            'bbox': self.region.format('bbox'),
            'bboxSR':4326,
            'imageSR':4326,
            'format':'png',
            'layers':'visible:4',
            'f':'pjson',
        }
        _req = f_utils.Fetch(self._usace_gs_api_url).fetch_req(params=_data)
        print(_req.url)
        if _req is not None:
            survey_list = _req.json()
            print(survey_list)
            fetch_fn = survey_list['href']
            self.results.append([fetch_fn, fetch_fn.split('/')[-1], 'globalelus'])
            
        return(self)

    def yield_xyz(self, entry):
        raise(NotImplementedError)

### End
