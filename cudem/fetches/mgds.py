### mgds.py -- Marine Geo Digital Library
##
## Copyright (c) 2022 Regents of the University of Colorado
##
## mgds.py is part of CUDEM
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
### Code:

import os
import lxml
from cudem import utils
from cudem import regions
from cudem import datasets

import cudem.fetches.utils as f_utils

class MGDS(f_utils.FetchModule):
    '''Fetch marine data from MGDS'''
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs) 

        self._mgds_file_url = "https://www.marine-geo.org/services/FileServer?"
        self._mgds_filedownload_url = "http://www.marine-geo.org/services/FileDownloadServer?"
        self._mgds_filemetadata_url = "http://www.marine-geo.org/services/FileDownloadServer/metadata?"      
        self._outdir = os.path.join(os.getcwd(), 'mgds')
        self.name = 'mgds'
        
    def run(self):
        '''Run the MGDS fetching module'''

        if self.region is None:
            return([])
        
        self.data = {
            'north':self.region.ymax,
            'west':self.region.xmin,
            'south':self.region.ymin,
            'east':self.region.xmax,
            'format':'summary',
            'data_type':'Bathymetry',
        }

        req = f_utils.Fetch(self._mgds_file_url).fetch_req(
            params=self.data, tries=10, timeout=2
        )
        if req is not None:
            req_xml = lxml.etree.fromstring(req.content)
            req_results = req_xml.findall('.//{https://www.marine-geo.org/services/xml/mgdsDataService}file')

            for req_result in req_results:
                name = req_result.attrib['name']
                link = req_result.attrib['download']
                self.results.append([link, name, 'mgds'])            
                
        return(self)

### End
