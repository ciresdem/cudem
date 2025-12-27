### buoys.py
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
##
## buoys.py is part of CUDEM
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
## Fetch NOAA Buoy Data.
##
### Code:

import lxml.html as lh
from typing import List, Optional
from cudem.fetches import fetches

## ==============================================
## Constants
## ==============================================
NDBC_URL = 'https://www.ndbc.noaa.gov'
BUOY_RADIAL_SEARCH_URL = 'https://www.ndbc.noaa.gov/radial_search.php?'
BUOY_REALTIME_URL = 'https://www.ndbc.noaa.gov/data/realtime2/'

## ==============================================
## Buoys Module
## ==============================================
class BUOYS(fetches.FetchModule):
    """NOAA BUOY data (beta)

    Fetch NOS Buoy Stations data.

    A sustainable and resilient marine observation and monitoring 
    infrastructure which enhances healthy ecosystems, communities, 
    and economies.

    https://www.ndbc.noaa.gov

    Configuration Example:
    < buoys:buoy_id=None >
    """
    
    def __init__(self, buoy_id: Optional[str] = None, **kwargs):
        super().__init__(name='buoys', **kwargs)
        self.buoy_id = buoy_id

        
    def run(self):
        """Run the BUOYS fetching module.
        
        Performs a radial search around the center of the provided region
        to find relevant buoy stations.
        """
        
        if self.region is None:
            return []

        ## Calculate center of the region for radial search
        rc = self.region.center()
        
        ## Search parameters for NDBC radial_search.php
        ## dist is search radius in nm (nautical miles)
        ## API defaults usually imply close proximity. 
        search_params = {
            'lat1': rc[1],
            'lon1': rc[0],
            'uom': 'M',    # Units: Metric
            'ot': 'A',     # Observation Time (A=All)
            'dist': 100,   # Distance
            'time': 0,
        }

        ## Fetch the HTML search results
        req = fetches.Fetch(
            BUOY_RADIAL_SEARCH_URL, 
            verbose=self.verbose
        ).fetch_req(params=search_params)
        
        if req is not None:
            try:
                doc = lh.document_fromstring(req.text)
                spans = doc.xpath('//span')
                current_stations = set()
                
                ## Parse HTML for station links
                for span in spans:
                    links = span.xpath('a')
                    if links:
                        href = links[0].get('href')
                        if href and 'station=' in href:
                            station_id = href.split('=')[-1]
                            current_stations.add(station_id)
                                
                ## Add Realtime Data Links for found stations
                for station_id in current_stations:
                    self.add_entry_to_results(
                        f"{BUOY_REALTIME_URL}{station_id}.txt",
                        f"buoy_results_{station_id}.txt",
                        'buoys'
                    )
            except Exception as e:
                ## Fallback or error logging
                if self.verbose:
                    utils.echo_error_msg(f"couln't parse Buoy results: {e}")

        return self

### End
