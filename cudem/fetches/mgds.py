### mgds.py -- Marine Geo Digital Library
##
## Copyright (c) 2022, 2023 Regents of the University of Colorado
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
    """The Marine Geoscience Data System (MGDS)

Fetch marine data from MGDS
    
MGDS is a trusted data repository that provides free public access to a curated collection of marine geophysical 
data products and complementary data related to understanding the formation and evolution 
of the seafloor and sub-seafloor.

https://www.marine-geo.org

data_tpye=[Bathymetry, Bathymetry:Phase, Bathymetry:Swath, Bathymetry:Swath:Ancillary, Bathymetry:Singlebeam, Bathymetry:BPI, Bathymetry:ReferenceSurface, Bathymetry:Paelobathymetry]
            
< mgds:data_type=Bathymetry >"""
    
    def __init__(self, data_type='Bathymetry', **kwargs):
        super().__init__(name='mgds', **kwargs) 
        self._mgds_file_url = "https://www.marine-geo.org/services/FileServer?"
        self._mgds_filedownload_url = "http://www.marine-geo.org/services/FileDownloadServer?"
        self._mgds_filemetadata_url = "http://www.marine-geo.org/services/FileDownloadServer/metadata?"
        self._mgds_archive_url = "http://www.marine-geo.org/services/FileDownloadServer/metadata?"
        self._mgds_search_url = "http://www.marine-geo.org/services/search/datasets??"
        self.data_type = data_type.replace(',', ':')
        
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
            'data_type':'{}'.format(self.data_type),
        }

        # self.data = {
        #     'north':self.region.ymax,
        #     'west':self.region.xmin,
        #     'south':self.region.ymin,
        #     'east':self.region.xmax,
        #     'data_type':'Bathymetry',
        # }
        # req = f_utils.Fetch(self._mgds_search_url).fetch_req(
        #     params=self.data, tries=10, timeout=2
        # )

        # if req is not None:
        #     req_xml = lxml.etree.fromstring(req.content)
                
        #     #req_results = req_xml.findall('.//data_set')
        #     for req_result in req_results:
        #         #rr_results = req_result.findall('.//ds_entry')
        #         for rr in rr_results:
        #             ds_id = rr.attrib['id']
        #             uids = req_result.attrib['uids'].split(',')
        #             fmt = req_result.attrib['file_formats']
        #             if fmt == 'MBSystem': fmt = 'MGD77'
        #             if 'GEOTIFF' in fmt: fmt = 'tif'
        #             print(ds_id, uids, fmt)
        #             #link = req_result.attrib['url']
        #             for uid in uids:
        #                 link = '{}data_uid={}'.format(self._mgds_filedownload_url, uid)
        #                 self.results.append([link, os.path.join(self._outdir, '{}_{}.{}'.format(ds_id, uid, fmt)), 'mgds'])            
        
        req = f_utils.Fetch(self._mgds_file_url).fetch_req(
            params=self.data, tries=10, timeout=2
        )
        if req is not None:
            req_xml = lxml.etree.fromstring(req.content)
            req_results = req_xml.findall('.//{https://www.marine-geo.org/services/xml/mgdsDataService}file')
            for req_result in req_results:
                name = req_result.attrib['name']
                link = req_result.attrib['download']
                self.results.append([link, os.path.join(self._outdir, name), 'mgds'])
                
        return(self)

### End
