### chs.py
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
##
## chs.py is part of CUDEM
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
## Fetch bathymetric soundings from the Canadian Hydrographic Service (CHS).
##
### Code:

import lxml.etree
from urllib.parse import urlencode
from typing import List, Dict, Optional, Any

from cudem import regions
from cudem.fetches import fetches

## ==============================================
## Constants
## ==============================================
CHS_WCS_URL = 'https://nonna-geoserver.data.chs-shc.ca/geoserver/wcs?'

## WCS Namespaces
NAMESPACES = {
    'gml': 'http://www.opengis.net/gml/3.2',
    'wcs': 'http://www.opengis.net/wcs/2.0',
    'ows': 'http://www.opengis.net/ows/2.0'
}

## ==============================================
## CHS Module
## ==============================================
class CHS(fetches.FetchModule):
    """Canadian Hydrographic Service Non-Navigational (NONNA) Bathymetric Data.

    Fetch bathymetric soundings from the CHS via WCS.
    
    https://open.canada.ca
    https://open.canada.ca/data/en/dataset/d3881c4c-650d-4070-bf9b-1e00aabf0a1d

    datatypes: 10 or 100 (resolution in meters)

    Configuration Example:
    < chs:datatype=10 >
    """
    
    def __init__(self, datatype: str = '100', **kwargs):
        super().__init__(name='chs', **kwargs)

        self.datatypes = ['10', '100']
        self.datatype = str(datatype) if str(datatype) in self.datatypes else '100'
        self._chs_url = CHS_WCS_URL

        
    def run(self):
        """Run the CHS fetching module."""

        if self.region is None:
            return []

        ## DescribeCoverage to get dataset bounds
        desc_data = {
            'request': 'DescribeCoverage',
            'version': '2.0.1',
            'CoverageID': f'nonna__NONNA {self.datatype} Coverage',
            'service': 'WCS'
        }
        
        req = fetches.Fetch(self._chs_url).fetch_req(params=desc_data)
        if req is None:
            return self

        try:
            ## Parse XML Response
            root = lxml.etree.fromstring(req.text.encode('utf-8'))
            
            ## Find Envelope (Bounding Box)
            envelope = root.find('.//gml:Envelope', namespaces=NAMESPACES)
            
            if envelope is None:
                ## Fallback for some server responses
                envelope = root.find('.//{http://www.opengis.net/gml/3.2}Envelope')

            if envelope is not None:
                lc_node = envelope.find('gml:lowerCorner', namespaces=NAMESPACES)
                uc_node = envelope.find('gml:upperCorner', namespaces=NAMESPACES)
                
                if lc_node is not None and uc_node is not None:
                    lc = [float(x) for x in lc_node.text.split()]
                    uc = [float(x) for x in uc_node.text.split()]
                    
                    ds_region = regions.Region().from_list([lc[1], uc[1], lc[0], uc[0]])
                    
                    ## Check Intersection
                    if regions.regions_intersect_ogr_p(self.region, ds_region):
                        
                        ## Construct GetCoverage URL
                        ## CHS Server Axis labels are 'Lat' and 'Long'.
                        wcs_params = [
                            ('request', 'GetCoverage'),
                            ('version', '2.0.1'),
                            ('CoverageID', f'nonna__NONNA {self.datatype} Coverage'),
                            ('service', 'WCS'),
                            ('subset', [f'Long({self.region.xmin},{self.region.xmax})',
                                        f'Long({self.region.ymin},{self.region.ymax})']),
                            ('subsettingcrs', 'http://www.opengis.net/def/crs/EPSG/0/4326'),
                            ('outputcrs', 'http://www.opengis.net/def/crs/EPSG/0/4326')
                        ]
                        
                        ## Generate URL without making a request
                        query_string = urlencode(wcs_params)
                        full_url = f"{self._chs_url}{query_string}"
                        
                        out_fn = f"chs_nonna{self.datatype}_{self.region.format('fn')}.tif"
                        
                        self.add_entry_to_results(full_url, out_fn, 'chs')

        except Exception as e:
            if self.verbose:
                print(f"Error parsing CHS WCS response: {e}")

        return self

### End
