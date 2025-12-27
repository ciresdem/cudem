### gmrt.py
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
##
## gmrt.py is part of CUDEM
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
## Fetch data from the Global Multi-Resolution Topography (GMRT) synthesis.
##
### Code:

from typing import Optional, Dict, Any, List
from cudem.fetches import fetches

## ==============================================
## Constants
## ==============================================
GMRT_URL = 'https://www.gmrt.org'
GMRT_POINT_URL = "https://www.gmrt.org:443/services/PointServer?"
GMRT_GRID_URL = "https://www.gmrt.org:443/services/GridServer?"
GMRT_GRID_URLS_URL = "https://www.gmrt.org:443/services/GridServer/urls?"
GMRT_METADATA_URL = "https://www.gmrt.org/services/GridServer/metadata?"
GMRT_SWATH_URL = "https://www.gmrt.org/shapefiles/gmrt_swath_polygons.zip"

# GMRT often requires specific user agents (Firefox/Windows)
GMRT_HEADERS = {
    'User-Agent': (
        'Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:89.0) '
        'Gecko/20100101 Firefox/89.0'
    )
}

## ==============================================
## GMRT Functions
## ==============================================
def gmrt_fetch_point(latitude: float, longitude: float) -> Optional[str]:
    """Fetch a single point elevation from GMRT."""
    
    data = {'longitude': longitude, 'latitude': latitude}
    
    req = fetches.Fetch(GMRT_POINT_URL, headers=GMRT_HEADERS).fetch_req(
        params=data, tries=10, timeout=2
    )
    
    if req is not None:
        return req.text
    return None


class GMRT(fetches.FetchModule):
    """The Global Multi-Resolution Topography synthesis.
    
    The Global Multi-Resolution Topography (GMRT) synthesis is a multi-resolutional 
    compilation of edited multibeam sonar data collected by scientists and 
    institutions worldwide, that is reviewed, processed and gridded by the GMRT 
    Team and merged into a single continuously updated compilation of global elevation 
    data.

    Data Formats:
      - GMT v3 Compatible NetCDF (GMT id=cf)
      - COARDS/CF1.6 Compliant NetCDF (GMT id=nd)
      - ESRI ArcASCII
      - GeoTIFF

    Layers: 'topo' or 'topo-mask'
    
    Data is assumed instantaneous MSL.

    Configuration Example:    
    < gmrt:res=max:fmt=geotiff:layer=topo >
    """
    
    def __init__(
            self,
            res: str = 'default',
            fmt: str = 'geotiff',
            layer: str = 'topo',
            want_swath: bool = False,
            **kwargs
    ):
        super().__init__(name='gmrt', **kwargs)
        
        self.res = res 
        self.fmt = fmt 
        self.want_swath = want_swath
        
        ## Validate layer
        self.layer = layer if layer in ['topo', 'topo-mask'] else 'topo'

        ## Buffer the input region and correct to wgs extremes
        ## GMRT specific: 2.33% buffer, 0.0088 increment
        if self.region is not None:
            self.gmrt_region = self.region.copy()
            self.gmrt_region.buffer(pct=2.33, x_inc=0.0088, y_inc=0.0088)
            self.gmrt_region._wgs_extremes(just_below=True)
        else:
            self.gmrt_region = None

        ## Metadata for DLIM/Processing
        self.data_format = 200
        self.src_srs = 'epsg:4326+3855'
        self.title = 'GMRT'
        self.source = 'GMRT'
        self.date = None
        self.data_type = 'Raster'
        self.resolution = None
        self.hdatum = 'wgs84'
        self.vdatum = 'msl'
        self.url = GMRT_URL
        self.headers = GMRT_HEADERS

        
    def run(self):
        """Run the GMRT fetching module."""
        
        if self.region is None or self.gmrt_region is None:
            return []

        self.data = {
            'north': self.gmrt_region.ymax,
            'west': self.gmrt_region.xmin,
            'south': self.gmrt_region.ymin,
            'east': self.gmrt_region.xmax,
            'mformat': 'json',
            'resolution': self.res,
            'format': self.fmt,
            'layer': self.layer
        }

        req = fetches.Fetch(
            GMRT_GRID_URL, 
            headers=self.headers
        ).fetch_req(
            params=self.data, 
            tries=10, 
            timeout=2
        )

        if req is not None:
            ext = 'tif' if self.fmt == 'geotiff' else 'grd'
            outf = f"gmrt_{self.layer}_{self.res}_{self.region.format('fn_full')}.{ext}"
            
            ## Add the fetched URL to results
            self.add_entry_to_results(req.url, outf, 'gmrt')
                
        return self

### End
