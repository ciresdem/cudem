### margrav.py
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

import warnings
from cudem import utils
from cudem.fetches import fetches

## MarGrav - Marine Gravity
class MarGrav(fetches.FetchModule):
    """MARine GRAVity Satellite Altimetry Topography from Scripps.

    Fetch mar_grav sattelite altimetry topography

    https://topex.ucsd.edu/WWW_html/mar_grav.html
    ftp://topex.ucsd.edu/pub/global_grav_1min/
    https://topex.ucsd.edu/marine_grav/explore_grav.html
    https://topex.ucsd.edu/marine_grav/white_paper.pdf
    https://topex.ucsd.edu/pub/global_grav_1min/geotiff/

    < mar_grav:upper_limit=None:lower_limit=None:raster=False:mag=1 >
    """
    
    def __init__(self, where='', mag=1, upper_limit=None,
                 lower_limit=None, raster=False,
                 **kwargs):
        super().__init__(name='mar_grav', **kwargs)
        self.mag = mag if mag == 1 else 0.1
        self.raster = raster # if True, grid the data and return as a raster in dlim
        self.where = [where] if len(where) > 0 else []
        
        ## The mar_grav URl
        self._mar_grav_url = 'https://topex.ucsd.edu/cgi-bin/get_data.cgi'
        self._mar_grav_geotiff_url = 'https://topex.ucsd.edu/pub/global_grav_1min/geotiff/'
        ## `get_data.cgi` is broken atm, use full grid url instrad
        self._mar_grav_27_1_url = 'https://topex.ucsd.edu/pub/global_topo_1min/topo_27.1.img'
        
        ## set up the region, restrict by z-region if desired
        self.grav_region = self.region.copy()
        self.grav_region._wgs_extremes(just_below=True)
        self.grav_region.zmax = utils.float_or(upper_limit)
        self.grav_region.zmin = utils.float_or(lower_limit)

        ## for dlim, data_format of 168 is xyz-file, skip the first line and reset
        ## the x values from 360 to 180
        self.data_format = '168:x_offset=REM:skip=1'
        self.src_srs = 'epsg:4326+3855'
        
        self.title = 'Marine Gravity Model'
        self.source = 'Scripps Institude of Oceanography'
        self.date = '2014'
        self.data_type = 'Bathymetric Raster'
        self.resolution = '60 arc-seconds'
        self.hdatum = 'WGS84'
        self.vdatum = 'MSL'
        self.url = self._mar_grav_url
        
        self.headers = {
            'User-Agent': ('Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:89.0) '
                           'Gecko/20100101 Firefox/89.0')
        }

            
    def run(self):
        """Run the mar_grav fetching module."""
        
        if self.region is None:
            return([])
        
        _data = {
            'north':self.region.ymax,
            'west':self.region.xmin,
            'south':self.region.ymin,
            'east':self.region.xmax,
            'mag':self.mag
        }
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            _req = fetches.Fetch(self._mar_grav_url, verify=False).fetch_req(params=_data)

        ## get_data.cgi is broken, but might not be forever, so we check that first and
        ## if it 404s then we just download the entire raster.
        if _req is not None:
            outf = 'mar_grav_{}.xyz'.format(self.region.format('fn_full'))
            self.add_entry_to_results(_req.url, outf, 'mar_grav')
        else:
            self.add_entry_to_results(
                self._mar_grav_27_1_url,
                'topo_27.1.img',
                'mar_grav_img'
            )

### End
