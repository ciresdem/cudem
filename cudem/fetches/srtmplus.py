### srtmplus.py
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
##
## srtmplus.py is part of CUDEM
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
## Fetch SRTM15+ global bathymetry/topography.
##
### Code:

from cudem import utils
from cudem.fetches import fetches

## ==============================================
## Constants
## ==============================================
SRTM_PLUS_CGI_URL = 'https://topex.ucsd.edu/cgi-bin/get_srtm15.cgi'

## ==============================================
## SRTM Plus Module
## ==============================================
class SRTMPlus(fetches.FetchModule):
    """SRTM15+: GLOBAL BATHYMETRY AND TOPOGRAPHY AT 15 ARCSECONDS.

    https://topex.ucsd.edu/WWW_html/srtm15_plus.html
    
    Configuration Example:
    < srtm_plus >
    """
    
    def __init__(self, **kwargs):
        super().__init__(name='srtm_plus', **kwargs)

        ## Metadata
        ## 168 is xyz-file, skip the first line (header)
        self.data_format = '168:skip=1'
        self.src_srs = 'epsg:4326+3855' # WGS84 + EGM2008

        self.title = 'SRTM+'
        self.source = 'Scripps Institution of Oceanography'
        self.date = '2014'
        self.data_type = 'Topographic/Bathymetric Raster'
        self.resolution = '15 arc-seconds' # Updated description to match 15+ name
        self.hdatum = 'WGS84'
        self.vdatum = 'MSL'
        self.url = SRTM_PLUS_CGI_URL

        
    def run(self):
        """Run the SRTM fetching module."""
        
        if self.region is None:
            return []

        ## Prepare CGI Parameters
        data = {
            'north': self.region.ymax,
            'west': self.region.xmin,
            'south': self.region.ymin,
            'east': self.region.xmax
        }
        
        ## Execute Request
        ## Verify=False is often needed for Scripps legacy CGI endpoints due to cert issues
        req = fetches.Fetch(self.url, verify=False).fetch_req(params=data)
        
        if req is not None:
            out_fn = f"srtm_{self.region.format('fn')}.xyz"
            self.add_entry_to_results(
                req.url,
                out_fn,
                'srtm'
            )
            
        return self

### End
