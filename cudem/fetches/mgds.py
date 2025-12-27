### mgds.py
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
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
## The above copyright notice and this permission notice shall be included in all
## copies or substantial portions of the Software.
##
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
## INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR 
## PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE 
## FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
## ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
## SOFTWARE.
##
### Commentary:
##
## Fetch marine geophysical data from the Marine Geoscience Data System (MGDS).
##
### Code:

import lxml.etree
from typing import List, Dict, Optional
from cudem.fetches import fetches

## ==============================================
## Constants
## ==============================================
MGDS_FILE_URL = 'https://www.marine-geo.org/services/FileServer?'
MGDS_FILE_DOWNLOAD_URL = 'http://www.marine-geo.org/services/FileDownloadServer?'
MGDS_METADATA_URL = 'http://www.marine-geo.org/services/FileDownloadServer/metadata?'
MGDS_SEARCH_URL = 'http://www.marine-geo.org/services/search/datasets??'

MGDS_NAMESPACE = 'https://www.marine-geo.org/services/xml/mgdsDataService'

## ==============================================
## MGDS Module
## ==============================================
class MGDS(fetches.FetchModule):
    """The Marine Geoscience Data System (MGDS)

    Fetch marine data from MGDS.
    
    MGDS is a trusted data repository that provides free public access 
    to a curated collection of marine geophysical data products.

    https://www.marine-geo.org

    data_type options:
    [Bathymetry, Bathymetry:Phase, Bathymetry:Swath, 
    Bathymetry:Swath:Ancillary, Bathymetry:Singlebeam, Bathymetry:BPI, 
    Bathymetry:ReferenceSurface, Bathymetry:Paelobathymetry]

    Configuration Example:            
    < mgds:data_type=Bathymetry >
    """
    
    def __init__(self, data_type: str = 'Bathymetry', **kwargs):
        super().__init__(name='mgds', **kwargs)
        self.data_type = data_type.replace(',', ':')

        
    def run(self):
        """Run the MGDS fetching module."""

        if self.region is None:
            return []

        search_params = {
            'north': self.region.ymax,
            'west': self.region.xmin,
            'south': self.region.ymin,
            'east': self.region.xmax,
            'format': 'summary',
            'data_type': self.data_type
        }

        ## Fetch Summary XML
        req = fetches.Fetch(MGDS_FILE_URL).fetch_req(
            params=search_params, 
            tries=10, 
            timeout=2
        )

        if req is not None:
            try:
                req_xml = lxml.etree.fromstring(req.content)
                
                ## Parse Results
                ## Using namespace directly in findall
                req_results = req_xml.findall(f'.//{{{MGDS_NAMESPACE}}}file')
                
                for req_result in req_results:
                    name = req_result.attrib.get('name')
                    link = req_result.attrib.get('download')
                    
                    if name and link:
                        self.add_entry_to_results(link, name, 'mgds')
                        
            except lxml.etree.XMLSyntaxError as e:
                pass

        return self

    
### End
