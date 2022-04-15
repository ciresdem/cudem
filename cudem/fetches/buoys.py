### buoys.py - STATIONS fetch
##
## Copyright (c) 2010 - 2021 CIRES Coastal DEM Team
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
## The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
##
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
## INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR 
## PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE 
## FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
## ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
##
### Commentary:
##
## Fetch bouy information from NOAA/NOS
##
### Code:

import os
import sys
import json
import lxml.html as lh

from osgeo import ogr

from cudem import utils
from cudem import regions
from cudem import datasets
from cudem import xyzfun

import cudem.fetches.utils as f_utils

## =============================================================================
##
## BUOYS Fetch
##
## Fetch NOAA bouy station information
##
## =============================================================================

class BUOYS(f_utils.FetchModule):
    """Fetch NOS Tide Stations"""
    
    def __init__(self, buoy_id=None, **kwargs):
        super().__init__(**kwargs)
        self._ndbc_url = 'https://www.ndbc.noaa.gov'
        self._buoy_box_search_url = 'https://www.ndbc.noaa.gov/box_search.php?'
        self._buoy_station_url = 'https://www.ndbc.noaa.gov/station_page.php?'
        self._buoy_stations_url = 'https://www.ndbc.noaa.gov/to_station.shtml'
        self._buoy_station_kml = 'https://www.ndbc.noaa.gov/kml/marineobs_by_owner.kml'
        self._buoy_station_realtime = 'https://www.ndbc.noaa.gov/data/realtime2/'
        self._outdir = os.path.join(os.getcwd(), 'buoys')
        self.name = 'buoys'
        self.buoy_id = buoy_id

    def run(self):
        '''Run the BOUYS fetching module'''
        
        if self.region is None:
            return([])

        _data = {
            'lat1': self.region.ymin,
            'lat2': self.region.ymax,
            'lon1': self.region.xmin,
            'lon2': self.region.xmax,
            'uom': 'M',
            'ot': 'A',
            'time': 0,
        }

        ## Fetch buoy ids from box search
        _req = f_utils.Fetch(
            self._buoy_box_search_url, verbose=self.verbose
        ).fetch_req(params=_data)
        if _req is not None:
            #print(_req.content)
            doc = lh.document_fromstring(_req.text)
            sp = doc.xpath('//span')
            current_stations = []
            
            for s in sp:
                #print(s.text_content())
                if len(s.xpath('a')) > 0:
                    station_url = s.xpath('a')[0].get('href')
                    if 'station=' in station_url:
                        station_id = station_url.split('=')[-1]
                        #self.results.append([self._ndbc_url + station_url, 'buoy_results_{}.html'.format(station_id), 'buoys'])
                        if station_id not in current_stations:
                            current_stations.append(station_id)
            for station_id in current_stations:
                self.results.append([self._buoy_station_realtime + station_id + '.txt', os.path.join(self._outdir, 'buoy_results_{}.txt'.format(station_id)), 'buoys'])
            
        return(self)

### End
