### gba.py
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
##
## gba.py is part of CUDEM
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
## Fetch Global Building Atlas (GBA) data via WFS.
##
### Code:

from urllib.parse import urlencode
from cudem.fetches import fetches

## ==============================================
## Constants
## ==============================================
GBA_WFS_URL = 'https://tubvsig-so2sat-vm1.srv.mwn.de/geoserver/ows?'

## ==============================================
## GBA Module
## ==============================================
class GBA(fetches.FetchModule):
    """Global Building Atlas (GBA) fetch module.
    
    Fetches Level of Detail 1 (LOD1) building data via WFS.

    Configuration Example:
    < gba >
    """
    
    def __init__(self, **kwargs):
        super().__init__(name='gba', **kwargs)

        
    def run(self):
        """Run the GBA fetching module."""

        if self.region is None:
            return []

        ## WFS GetFeature Parameters
        wfs_params = {
            'request': 'GetFeature',
            'version': '2.0.0',
            'typenames': 'lod1_global',
            'service': 'WFS',
            'bbox': f"{self.region.format('osm_bbox')},http://www.opengis.net/def/crs/EPSG/0/4326",
            'outputformat': 'json',
            'targetcrs': 'http://www.opengis.net/def/crs/EPSG/0/4326'
        }

        ## Construct URL
        query_string = urlencode(wfs_params)
        full_url = f"{GBA_WFS_URL}{query_string}"
        
        ## Define Output Filename
        out_fn = f"gba_{self.region.format('fn')}.geojson"
        
        # Add to results
        self.add_entry_to_results(full_url, out_fn, 'gba')

        return self

### End
