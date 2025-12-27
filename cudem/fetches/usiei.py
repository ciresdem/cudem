### usiei.py
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
##
## usiei.py is part of CUDEM
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
## Query the US Interagency Elevation Inventory (USIEI) via ArcGIS REST API.
##
### Code:

import json
from typing import Dict, Any
from cudem import utils
from cudem.fetches import fetches

## ==============================================
## Constants
## ==============================================
USIEI_MAP_SERVER_URL = (
    'https://coast.noaa.gov/arcgis/rest/services/'
    'USInteragencyElevationInventory/USIEIv2/MapServer'
)

## ==============================================
## USIEI Module
## ==============================================
class USIEI(fetches.FetchModule):
    """US Interagency Elevation Inventory (USIEI)

    No data is downloaded with this module. It lists query results from the USIEI.
    Set 'want_geometry' to True to output a GeoJSON formatted vector.

    Layers:
    0 - Lidar-Topobathy
    1 - Lidar-Bathy
    2 - Lidar-Topo
    3 - IfSAR/InSAR
    4 - Other Bathy

    https://coast.noaa.gov/inventory/

    Configuration Example:
    < usiei:where='1=1':layer=0:want_geometry=False >
    """
    
    def __init__(self, where: str = '1=1', want_geometry: bool = False, layer: int = 0, **kwargs):
        super().__init__(name='usiei', **kwargs)
        self.where = where
        self.want_geometry = want_geometry
        self._usiei_query_url = f"{USIEI_MAP_SERVER_URL}/{layer}/query?"

        
    def run(self):
        """Run the USIEI fetches module."""
        
        if self.region is None:
            return []
        
        ## Prepare ArcGIS REST Query
        params = {
            'where': self.where,
            'outFields': '*',
            'geometry': self.region.format('bbox'),
            'inSR': 4326,
            'outSR': 4326,
            'f': 'pjson',
            'returnGeometry': 'True' if self.want_geometry else 'False',
        }

        ## Execute Request
        req = fetches.Fetch(
            self._usiei_query_url,
            verbose=self.verbose
        ).fetch_req(params=params)

        if req is not None:
            if self.want_geometry:
                # Print raw response (likely GeoJSON/PJSON structure)
                print(req.text)
            else:
                try:
                    features = req.json()
                    ## Pretty print the features list
                    if 'features' in features:
                        print(json.dumps(features['features'], indent=4))
                    else:
                        print(json.dumps(features, indent=4))
                except json.JSONDecodeError:
                    utils.echo_error_msg("Failed to parse USIEI response.")
            
        return self

### End
