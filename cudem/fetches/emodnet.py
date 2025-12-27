### emodnet.py
##
## Copyright (c) 2021 - 2025 Regents of the University of Colorado
##
## emodnet.py is part of CUDEM
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
## Fetch elevation data from EMODNet via WCS or ERDDAP.
##
### Code:

from urllib.parse import urlencode
from typing import List, Optional, Any
from cudem.fetches import fetches

## ==============================================
## Constants
## ==============================================
EMODNET_WCS_URL = 'https://ows.emodnet-bathymetry.eu/wcs?'
EMODNET_ERDDAP_BASE = 'https://erddap.emodnet.eu/erddap/griddap/dtm_2020_v2_e0bf_e7e4_5b8f'

## ==============================================
## EMODNet Module
## ==============================================
class EMODNet(fetches.FetchModule):
    """EU elevation data extracts from EMOD DTM.

    Fetch raster data from the EMODNET DTM via WCS (default) or ERDDAP.
    
    https://portal.emodnet-bathymetry.eu/
    https://erddap.emodnet.eu

    Configuration Example:
    < emodnet:want_erddap=False:erddap_format=nc >
    """

    def __init__(self, want_erddap: bool = False, erddap_format: str = 'nc', **kwargs):
        super().__init__(name='emodnet', **kwargs)
        self.want_erddap = want_erddap
        self.erddap_format = erddap_format
        self.src_srs = 'epsg:4326'

        
    def run(self):
        """Run the EMODNET fetching module."""
        
        if self.region is None:
            return []

        if self.want_erddap:
            ## Construct ERDDAP URL
            ## Syntax: .format?variable[(min):1:(max)][(min):1:(max)]
            ## EMODNet DTM variable is 'elevation'            
            query = f"elevation[({self.region.ymin}):1:({self.region.ymax})][({self.region.xmin}):1:({self.region.xmax})]"
            
            erddap_url = f"{EMODNET_ERDDAP_BASE}.{self.erddap_format}?{query}"
            out_fn = f"emodnet_{self.region.format('fn')}.{self.erddap_format}"
            
            self.add_entry_to_results(erddap_url, out_fn, self.erddap_format)
            
        else:
            ## WCS Request
            ## EMODNet DTM is ~1/16 arc minute (~115m).
            
            wcs_params = {
                'service': 'WCS',
                'request': 'GetCoverage',
                'version': '1.0.0',
                'Identifier': 'emodnet:mean',
                'coverage': 'emodnet:mean',
                'format': 'GeoTIFF',
                'bbox': self.region.format('bbox'),
                'crs': 'EPSG:4326',
                'resx': 0.00104166666666667, # Approx 1/16 arc min if needed
                'resy': 0.00104166666666667,
            }

            full_url = f"{EMODNET_WCS_URL}{urlencode(wcs_params)}"
            out_fn = f"emodnet_{self.region.format('fn')}.tif"
            
            self.add_entry_to_results(full_url, out_fn, 'emodnet')
            
        return self

    
### End
