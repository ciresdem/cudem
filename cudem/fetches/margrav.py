### margrav.py
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
##
## margrav.py is part of CUDEM
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
## Fetch Marine Gravity data from Scripps Institution of Oceanography.
##
### Code:

import warnings
from typing import Optional, List
from cudem import utils
from cudem.fetches import fetches

## ==============================================
## Constants
## ==============================================
MARGRAV_CGI_URL = 'https://topex.ucsd.edu/cgi-bin/get_data.cgi'
MARGRAV_GEOTIFF_URL = 'https://topex.ucsd.edu/pub/global_grav_1min/geotiff/'
MARGRAV_IMG_URL = 'https://topex.ucsd.edu/pub/global_topo_1min/topo_27.1.img'

## ==============================================
## MarGrav Module
## ==============================================
class MarGrav(fetches.FetchModule):
    """MARine GRAVity Satellite Altimetry Topography from Scripps.

    Fetch mar_grav satellite altimetry topography.

    https://topex.ucsd.edu/WWW_html/mar_grav.html
    
    Configuration Example:
    < mar_grav:upper_limit=None:lower_limit=None:raster=False:mag=1 >
    """
    
    def __init__(self, where: str = '', mag: float = 1, upper_limit: Optional[float] = None,
                 lower_limit: Optional[float] = None, raster: bool = False,
                 **kwargs):
        super().__init__(name='mar_grav', **kwargs)
        
        self.mag = mag if mag == 1 else 0.1
        self.raster = raster # if True, grid the data and return as a raster in dlim
        self.where = [where] if where else []
        
        ## Region configuration
        if self.region:
            self.grav_region = self.region.copy()
            self.grav_region._wgs_extremes(just_below=True)
            self.grav_region.zmax = utils.float_or(upper_limit)
            self.grav_region.zmin = utils.float_or(lower_limit)
        else:
            self.grav_region = None

        ## DLIM Configuration
        ## 168 is xyz-file, skip the first line,
        ## x_offset=REM handles 0-360 conversion
        self.data_format = '168:x_offset=REM:skip=1'
        self.src_srs = 'epsg:4326+3855'
        
        ## Metadata
        self.title = 'Marine Gravity Model'
        self.source = 'Scripps Institution of Oceanography'
        self.date = '2014'
        self.data_type = 'Bathymetric Raster'
        self.resolution = '60 arc-seconds'
        self.hdatum = 'WGS84'
        self.vdatum = 'MSL'
        self.url = MARGRAV_CGI_URL
        
        self.headers = {
            'User-Agent': ('Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:89.0) '
                           'Gecko/20100101 Firefox/89.0')
        }

    def run(self):
        """Run the mar_grav fetching module."""
        
        if self.region is None:
            return []
        
        ## Prepare CGI Parameters
        data = {
            'north': self.region.ymax,
            'west': self.region.xmin,
            'south': self.region.ymin,
            'east': self.region.xmax,
            'mag': self.mag
        }

        ## Attempt to use the CGI service
        ## Suppress SSL warnings as Topex servers often have cert issues
        req = None
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            try:
                req = fetches.Fetch(MARGRAV_CGI_URL, verify=False, verbose=self.verbose).fetch_req(params=data)
            except Exception:
                req = None

        ## Check if CGI worked
        if req is not None and req.status_code == 200:
            outf = f"mar_grav_{self.region.format('fn_full')}.xyz"
            self.add_entry_to_results(req.url, outf, 'mar_grav')
        else:
            ## Fallback to downloading the full grid IMG file if CGI fails
            if self.verbose:
                utils.echo_warning_msg("MarGrav CGI failed; falling back to full grid download.")
                
            self.add_entry_to_results(
                MARGRAV_IMG_URL,
                'topo_27.1.img',
                'mar_grav_img'
            )

        return self

### End
