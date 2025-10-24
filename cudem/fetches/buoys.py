### buoys.py
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
##
## fetches.py is part of CUDEM
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
###############################################################################
### Commentary:
##
##
### Code:

import lxml.html as lh
from cudem.fetches import fetches

## buoys 
class BUOYS(fetches.FetchModule):
    """NOAA BUOY data (beta)

    Fetch NOS Buoy Stations

    A sustainable and resilient marine observation and monitoring 
    infrastructure which enhances healthy ecosystems, communities, 
    and economies in the face of change and To provide quality observations in 
    the marine environment in a safe and sustainable manner to support 
    the understanding of and predictions  to changes in weather, climate, 
    oceans and coast. 

    https://www.ndbc.noaa.gov

    < buoys:buoy_id=None >
    """
    
    def __init__(self, buoy_id = None, **kwargs):
        super().__init__(name='buoys', **kwargs)
        self.buoy_id = buoy_id
        
        ## various buoy URLs
        self._ndbc_url = 'https://www.ndbc.noaa.gov'
        self._buoy_box_search_url = 'https://www.ndbc.noaa.gov/box_search.php?'
        self._buoy_radial_search_url = 'https://www.ndbc.noaa.gov/radial_search.php?'
        self._buoy_station_url = 'https://www.ndbc.noaa.gov/station_page.php?'
        self._buoy_stations_url = 'https://www.ndbc.noaa.gov/to_station.shtml'
        self._buoy_station_kml = 'https://www.ndbc.noaa.gov/kml/marineobs_by_owner.kml'
        self._buoy_station_realtime = 'https://www.ndbc.noaa.gov/data/realtime2/'

        
    def run(self):
        '''Run the BOUYS fetching module'''
        
        if self.region is None:
            return([])

        _data_box = {
            'lat1': self.region.ymin,
            'lat2': self.region.ymax,
            'lon1': self.region.xmin,
            'lon2': self.region.xmax,
            'uom': 'M',
            'ot': 'A',
            'time': 0,
        }
        rc = self.region.center()
        #print(rc)
        _data = {
            'lon1': rc[0],
            'lat1': rc[1],
            #'lat1': self.region.ymin,
            #'lat2': self.region.ymax,
            'lon1': self.region.xmin,
            #'lon2': self.region.xmax,
            'uom': 'M',
            'ot': 'A',
            'dist': 100,
            'time': 0,
        }

        ## Fetch buoy ids from box search
        _req = fetches.Fetch(
            self._buoy_radial_search_url, verbose=self.verbose
        ).fetch_req(params=_data)
        if _req is not None:
            #print(_req.content)
            #print(_req.url)
            doc = lh.document_fromstring(_req.text)
            sp = doc.xpath('//span')
            current_stations = []
            
            for s in sp:
                #print(s.text_content())
                if len(s.xpath('a')) > 0:
                    station_url = s.xpath('a')[0].get('href')
                    if 'station=' in station_url:
                        station_id = station_url.split('=')[-1]
                        if station_id not in current_stations:
                            current_stations.append(station_id)
                            
            for station_id in current_stations:
                self.add_entry_to_results(
                    self._buoy_station_realtime + station_id + '.txt',
                    'buoy_results_{}.txt'.format(station_id),
                    'buoys'
                )

        return(self)

### End
