### csb.py
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
##
## csb.py is part of CUDEM
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
## Fetch Crowd Sourced Bathymetry (CSB) from NOAA.
##
### Code:

import json
from typing import Optional, Dict, List
from cudem.fetches import fetches

## ==============================================
## Constants
## ==============================================
CSB_DATA_URL = 'https://noaa-dcdb-bathymetry-pds.s3.amazonaws.com'
CSB_MAP_SERVER = 'https://gis.ngdc.noaa.gov/arcgis/rest/services/csb/MapServer'

## ==============================================
## CSB Module
## ==============================================
class CSB(fetches.FetchModule):
    """Crowd Sourced Bathymetry (CSB) from NOAA

    Fetch CSB data, optionally filtering by region and year.
    
    < csb:where='1=1':layer=1:index=False >
    """
    
    def __init__(self, where: str = '1=1', layer: int = 1, index: bool = False, **kwargs):
        super().__init__(name='csb', **kwargs)
        self.where = where
        self.layer = layer
        self.index = index

        self._csb_query_url = f'{CSB_MAP_SERVER}/{self.layer}/query?'
        
        ## Metadata / Processing Defaults
        self.data_format = '168:skip=1:xpos=2:ypos=3:zpos=4:z_scale=-1:delim=,'
        self.src_srs = 'epsg:4326+5866'
        self.title = 'Crowd Sourced Bathymetry'
        self.source = 'NOAA/NCEI'
        self.date = '2024'
        self.data_type = 'Bathymetric Soundings'
        self.resolution = '<10m to several kilometers'
        self.hdatum = 'WGS84'
        self.vdatum = 'MSL'
        self.url = CSB_DATA_URL

    def run(self):
        """Run the CSB fetches module."""
        
        if self.region is None:
            return []

        ## Construct Where Clause with Date Filters
        ## self.min_year/max_year are populated by the parent FetchModule
        ## if passed in kwargs
        current_where = self.where
        
        if self.min_year is not None:
            current_where += f" AND YEAR >= {self.min_year}"
            
        if self.max_year is not None:
            current_where += f" AND YEAR <= {self.max_year}"

        ## Prepare ArcGIS Query
        query_params = {
            'where': current_where,
            'outFields': '*',
            'geometry': self.region.format('bbox'),
            'inSR': 4326,
            'outSR': 4326,
            'f': 'pjson',
            'returnGeometry': 'False'
        }

        ## Fetch Results
        req = fetches.Fetch(
            self._csb_query_url, 
            verbose=self.verbose
        ).fetch_req(params=query_params)

        if req is not None:
            try:
                features_json = req.json()
            except ValueError:
                return self

            if 'features' in features_json:
                for feature in features_json['features']:
                    attributes = feature.get('attributes', {})
                    
                    if self.index:
                        print(json.dumps(attributes, indent=4))
                    else:
                        name = attributes.get('NAME')
                        if not name:
                            continue
                            
                        ## Parse directory structure from NAME (Format: YYYYMMDD...)
                        ## CSB Structure: /csb/csv/YYYY/MM/DD/NAME_pointData.csv
                        try:
                            year = name[:4]
                            dir_a = name[4:6]
                            dir_b = name[6:8]
                            
                            ## Remove extension from name if present (usually .tar.gz)
                            base_name = name[:-7] if len(name) > 7 else name
                            csv_fn = f'{base_name}_pointData.csv'
                            
                            link = f'{CSB_DATA_URL}/csb/csv/{year}/{dir_a}/{dir_b}/{csv_fn}'

                            ## Set date for this specific entry
                            self.date = str(attributes.get('YEAR', self.date))
                            
                            self.add_entry_to_results(link, csv_fn, 'csb')
                        except Exception:
                            continue

        return self

### End
