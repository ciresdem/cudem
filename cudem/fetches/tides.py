### tides.py
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
##
## tides.py is part of CUDEM
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
## Fetch tide station information from NOAA/NOS via ArcGIS REST API.
##
### Code:

from urllib.parse import urlencode
from cudem.fetches import fetches

## ==============================================
## Constants
## ==============================================
TIDES_API_URL = 'https://mapservices.weather.noaa.gov/static/rest/services/NOS_Observations/CO_OPS_Products/FeatureServer/0/query?'

## ==============================================
## Tides Module
## ==============================================
class Tides(fetches.FetchModule):
    """TIDE station information from NOAA/NOS

    Fetch NOS Tide Stations. Fetched file is a GeoJSON (or pjson) 
    containing records for stations within the region.
    
    https://tidesandcurrents.noaa.gov/

    Configuration Example:
    < tides >
    """
    
    def __init__(self, **kwargs):
        super().__init__(name='tides', **kwargs)
        self.src_srs = 'epsg:4326'

    def run(self):
        """Run the TIDES fetching module."""
        
        if self.region is None:
            return []
        
        ## Prepare ArcGIS REST Query
        params = {
            'outFields': '*',
            'units': 'esriSRUnit_Meter',
            'geometry': self.region.format('bbox'),
            'inSR': 4326,
            'outSR': 4326,
            'f': 'pjson',
        }
        
        ## Construct URL
        query_string = urlencode(params)
        full_url = f"{TIDES_API_URL}{query_string}"
        
        ## Define Output Filename
        out_fn = f"tides_results_{self.region.format('fn')}.json"
        
        ## Add to results
        self.add_entry_to_results(
            full_url,
            out_fn,
            'tides'
        )
            
        return self

### End
