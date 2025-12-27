### trackline.py
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
##
## trackline.py is part of CUDEM
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
## Query NOAA Trackline bathymetric data.
##
### Code:

import json
from cudem import utils
from cudem.fetches import fetches

## ==============================================
## Constants
## ==============================================
TRACKLINE_BASE_URL = 'https://gis.ngdc.noaa.gov/arcgis/rest/services/web_mercator/trackline_combined_dynamic/MapServer/1'
TRACKLINE_QUERY_URL = f'{TRACKLINE_BASE_URL}/query?'
TRACKLINE_REQUEST_URL = 'http://www.ngdc.noaa.gov/trackline/request/?surveyIds='

## ==============================================
## Trackline Module
## ==============================================
class Trackline(fetches.FetchModule):
    """NOAA TRACKLINE bathymetric data.

    http://www.ngdc.noaa.gov/trackline/

    Note: This module currently does not download files directly. 
    It generates a URL for a "shopping basket" containing the survey IDs 
    found in the region, which must be submitted manually.

    < trackline >
    """
    
    def __init__(self, where: str = '1=1', **kwargs):
        super().__init__(name='trackline', **kwargs)
        self.where = where

        
    def run(self):
        """Run the trackline fetching module."""
        
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
            'returnGeometry': 'False'
        }

        ## Execute Request
        req = fetches.Fetch(
            TRACKLINE_QUERY_URL, 
            verbose=self.verbose
        ).fetch_req(params=params)

        if req is not None:
            try:
                features = req.json()
                
                ## Extract Survey IDs
                ids = []
                if 'features' in features:
                    ids = [
                        f['attributes']['SURVEY_ID'] 
                        for f in features['features'] 
                        if 'attributes' in f and 'SURVEY_ID' in f['attributes']
                    ]

                if ids:
                    basket_link = f"{TRACKLINE_REQUEST_URL}{','.join(ids)}"
                    print(f"Trackline Basket URL: {basket_link}")
                else:
                    if self.verbose:
                        utils.echo_msg("No Trackline surveys found in this region.")

            except json.JSONDecodeError:
                utils.echo_error_msg("Failed to parse Trackline response.")
            except Exception as e:
                utils.echo_error_msg(f"Error processing Trackline data: {e}")

        return self

### End
