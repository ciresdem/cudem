### waterservices.py
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
##
## waterservices.py is part of CUDEM
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
## Fetch WaterServices station information from USGS.
##
### Code:

from urllib.parse import urlencode
from cudem import utils
from cudem.fetches import fetches

## ==============================================
## Constants
## ==============================================
WATER_SERVICES_API_URL = 'https://waterservices.usgs.gov/nwis/iv/?'

## ==============================================
## WaterServices Module
## ==============================================
class WaterServices(fetches.FetchModule):
    """WaterServices station information from USGS

    Fetches active site data in JSON format for the specified region.

    https://waterservices.usgs.gov/

    Configuration Example:
    < waterservices:printout=False >
    """
    
    def __init__(self, printout: bool = False, **kwargs):
        super().__init__(name='waterservices', **kwargs)
        self.printout = printout

        
    def run(self):
        """Run the WATERSERVICES fetching module."""
        
        if self.region is None:
            return []
        
        ## Prepare API Parameters
        params = {
            'bBox': self.region.format('bbox'),
            'siteStatus': 'active',
            'format': 'json',
        }
        
        ## Construct URL
        query_string = urlencode(params)
        full_url = f"{WATER_SERVICES_API_URL}{query_string}"
        
        ## Output filename
        out_fn = f"water_services_results_{self.region.format('fn')}.json"
        
        ## Add to results
        self.add_entry_to_results(full_url, out_fn, 'waterservices')

        ## Optional: Print station information immediately
        if self.printout:
            self._print_station_info(full_url)
            
        return self

    
    def _print_station_info(self, url: str):
        """Fetch and print summary information for stations found."""
        
        try:
            req = fetches.Fetch(url, verbose=self.verbose).fetch_req()
            if req is None:
                return

            data = req.json()
            time_series = data.get('value', {}).get('timeSeries', [])
            
            for item in time_series:
                try:
                    ## Extract Geolocation
                    geo_loc = item['sourceInfo']['geoLocation']['geogLocation']
                    lat = geo_loc['latitude']
                    lon = geo_loc['longitude']

                    ## Extract Variable Info
                    variable = item['variable']
                    var_code = variable['variableCode'][0]['value']
                    var_name = variable['variableName']

                    ## Extract Value
                    values = item['values'][0]['value']
                    value = values[0]['value'] if values else 'N/A'

                    print(f"{var_name} ({var_code}) at {lon}, {lat}: {value}")

                except (KeyError, IndexError) as e:
                    if self.verbose:
                        utils.echo_warning_msg(f"Skipping malformed station record: {e}")
                        
        except Exception as e:
            utils.echo_error_msg(f"Failed to print WaterServices info: {e}")

### End
