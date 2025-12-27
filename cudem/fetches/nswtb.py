### nswtb.py
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
##
## nswtb.py is part of CUDEM
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
## Fetch New South Wales (NSW) Topo-Bathy DEM contours via ArcGIS REST API.
##
### Code:

import json
from typing import Optional, Dict, Any
from cudem import utils
from cudem.fetches import fetches

## ==============================================
## Constants
## ==============================================
NSW_MAP_SERVER = (
    'https://mapprod2.environment.nsw.gov.au/arcgis/rest/services/'
    'Coastal_Marine/NSW_Marine_Lidar_Bathymetry_Data_2018/MapServer'
)

## ==============================================
## NSW_TB Module
## ==============================================
class NSW_TB(fetches.FetchModule):
    """New South Wales Topo-Bathy DEM

    Fetches Topo-Bathy data (contours) from NSW environment services.

    < nsw_tb:where='1=1':layer=0:index=False >
    """
    
    def __init__(self, where: str = '1=1', layer: int = 0, index: bool = False, **kwargs):
        super().__init__(name='nsw_tb', **kwargs)
        self.where = where
        self.index = index
        self.src_srs = None

        self._nsw_query_url = f'{NSW_MAP_SERVER}/{layer}/query?'

    def run(self):
        """Run the NSW_TB fetching module."""
        
        if self.region is None:
            return []

        ## Prepare ArcGIS Query
        params = {
            'where': self.where,
            'outFields': '*',
            'geometry': self.region.format('bbox'),
            'inSR': 4326,
            'outSR': 4326,
            'f': 'geojson',
            'returnGeometry': 'True',
            'geometryType': 'esriGeometryEnvelope',
            'spatialRel': 'esriSpatialRelIntersects'
        }
        
        ## Execute Request to verify features exist and get full URL
        req = fetches.Fetch(
            self._nsw_query_url, 
            verbose=self.verbose
        ).fetch_req(params=params)

        if req is not None:
            try:
                features = req.json()
                if features.get('features'):
                    ## Output filename
                    geojson_fn = f"nsw_{self.region.format('fn')}_contours.geojson"
                    
                    ## Add the fully constructed URL to results
                    self.add_entry_to_results(req.url, geojson_fn, 'nsw_contours')
                    
            except json.JSONDecodeError:
                utils.echo_error_msg("Failed to parse NSW_TB response.")

        return self

    
### End
