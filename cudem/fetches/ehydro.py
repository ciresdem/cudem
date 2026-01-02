### ehydro.py
##
## Copyright (c) 2013 - 2026 Regents of the University of Colorado
##
## ehydro.py is part of CUDEM
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
## Fetch USACE eHydro bathymetric data via ArcGIS REST API.
##
### Code:

import json
import datetime
from typing import List, Optional, Any
from cudem import utils
from cudem.fetches import fetches

## ==============================================
## Constants
## ==============================================
EHYDRO_API_URL = 'https://services7.arcgis.com/n1YM8pTrFmm7L4hs/arcgis/rest/services/eHydro_Survey_Data/FeatureServer/0'
EHYDRO_QUERY_URL = f'{EHYDRO_API_URL}/query?'

## ==============================================
## eHydro Module
## ==============================================
class eHydro(fetches.FetchModule):
    """USACE eHydro bathymetric data.
    
    The eHydro dataset supports the USACE navigation mission by providing 
    bathymetric survey data for navigation channels and harbors.

    https://navigation.usace.army.mil/Survey/Hydro

    Configuration Example:
    < ehydro:where='1=1':inc=None:index=False:tables=False >
    """

    def __init__(self, where: str = '1=1', inc: Optional[str] = None, 
                 survey_name: Optional[str] = None, index: bool = False,
                 tables: bool = False, **kwargs):
        super().__init__(name='ehydro', **kwargs)
        self.where = where
        self.survey_name = survey_name
        self.inc = utils.str2inc(inc)
        self.index = index
        self.tables = tables
        self._ehydro_query_url = EHYDRO_QUERY_URL

        
    def _parse_timestamp(self, timestamp: Any) -> Optional[int]:
        """Safely parse ESRI timestamp (milliseconds) to year."""
        
        try:
            ## ESRI timestamps are in milliseconds
            seconds = int(str(timestamp)[:10])
            dt = datetime.datetime.fromtimestamp(seconds)
            return dt.year
        except (ValueError, TypeError):
            return None

    def run(self):
        """Run the eHydro fetching module."""
        
        if self.region is None:
            return []
        
        ## Prepare ArcGIS Query
        params = {
            'where': self.where,
            'outFields': '*',
            'geometry': self.region.format('bbox'),
            'inSR': 4326,
            'outSR': 4326,
            'f': 'pjson',
            'returnGeometry': 'False'
        }
        
        ## Fetch Data
        req = fetches.Fetch(self._ehydro_query_url, verbose=self.verbose).fetch_req(params=params)
        
        if req is None:
            return self

        try:
            response = req.json()
        except json.JSONDecodeError:
            utils.echo_error_msg("Failed to parse eHydro JSON response.")
            return self

        if 'features' not in response:
            utils.echo_error_msg("No features in request.")
            return self

        ## Process Features
        for feature in response['features']:
            attrs = feature.get('attributes', {})
            
            sid = attrs.get('sdsmetadataid')
            fetch_fn = attrs.get('sourcedatalocation')
            
            if not fetch_fn:
                continue

            ## Parse Dates
            start_ts = attrs.get('surveydatestart')
            year = self._parse_timestamp(start_ts)

            ## Filter by Survey Name
            if self.survey_name is not None:
                check_id = sid if sid else fetch_fn
                ## Check if any part of the requested name matches
                matches = [x in check_id for x in self.survey_name.split('/')]
                if not any(matches):
                    continue

            ## Filter by Date (inherited from FetchModule)
            if year:
                if self.min_year is not None and year < self.min_year:
                    continue
                if self.max_year is not None and year > self.max_year:
                    continue

            ## Output Handling
            if self.index:
                ## Dump full metadata
                print(json.dumps(attrs, indent=4))
                
            elif self.tables:
                ## Print CSV-like summary
                line = f"{sid},{year}"
                if sid is not None:
                    print(line)
                    
            else:
                ## Add to download list
                self.add_entry_to_results(
                    fetch_fn,
                    fetch_fn.split('/')[-1],
                    'ehydro'
                )
                
        return self

### End
