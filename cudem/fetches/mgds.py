### mgds.py
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
##
## fetches.py is part of CUDEM
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
###############################################################################
### Commentary:
##
##
### Code:

import lxml.etree
from cudem.fetches import fetches

## MGDS
class MGDS(fetches.FetchModule):
    """The Marine Geoscience Data System (MGDS)

    Fetch marine data from MGDS
    
    MGDS is a trusted data repository that provides free public access 
    to a curated collection of marine geophysical data products and 
    complementary data related to understanding the formation and evolution 
    of the seafloor and sub-seafloor.

    https://www.marine-geo.org

    data_tpye=[Bathymetry, Bathymetry:Phase, Bathymetry:Swath, 
    Bathymetry:Swath:Ancillary, Bathymetry:Singlebeam, Bathymetry:BPI, 
    Bathymetry:ReferenceSurface, Bathymetry:Paelobathymetry]
            
    < mgds:data_type=Bathymetry >
    """
    
    def __init__(self, data_type = 'Bathymetry', **kwargs):
        super().__init__(name='mgds', **kwargs)
        self.data_type = data_type.replace(',', ':')
        
        ## The various MGDS URLs
        self._mgds_file_url = 'https://www.marine-geo.org/services/FileServer?'
        self._mgds_filedownload_url = 'http://www.marine-geo.org/services/FileDownloadServer?'
        self._mgds_filemetadata_url = ('http://www.marine-geo.org/services/'
                                       'FileDownloadServer/metadata?')
        self._mgds_archive_url = 'http://www.marine-geo.org/services/FileDownloadServer/metadata?'
        self._mgds_search_url = 'http://www.marine-geo.org/services/search/datasets??'

        
    def run(self):
        """Run the MGDS fetching module"""

        if self.region is None:
            return([])

        self.data = {
            'north':self.region.ymax,
            'west':self.region.xmin,
            'south':self.region.ymin,
            'east':self.region.xmax,
            'format':'summary',
            'data_type':'{}'.format(self.data_type)
        }
        req = fetches.Fetch(self._mgds_file_url).fetch_req(
            params=self.data, tries=10, timeout=2
        )
        if req is not None:
            req_xml = lxml.etree.fromstring(req.content)
            req_results = req_xml.findall(
                './/{https://www.marine-geo.org/services/xml/mgdsDataService}file'
            )
            for req_result in req_results:
                name = req_result.attrib['name']
                link = req_result.attrib['download']
                self.add_entry_to_results(link, name, 'mgds')
                
        return(self)

### End
