### csb.py
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

import json
from cudem.fetches import fetches

class CSB(fetches.FetchModule):
    """crowd sourced bathymetry from NOAA

    < csv:where=None:layer=0:index=False >
    """
    
    def __init__(self, where='1=1', layer=1, index=False, **kwargs):
        super().__init__(name='csb', **kwargs)
        self.where = where
        self.index = index

        ## various CSB URLs
        self._csb_data_url = 'https://noaa-dcdb-bathymetry-pds.s3.amazonaws.com'
        self._csb_map_server = 'https://gis.ngdc.noaa.gov/arcgis/rest/services/csb/MapServer'
        self._csb_query_url = '{0}/{1}/query?'.format(self._csb_map_server, layer)
        
        ## for dlim
        self.data_format = '168:skip=1:xpos=2:ypos=3:zpos=4:z_scale=-1:delim=,'
        self.src_srs = 'epsg:4326+5866'

        self.title = 'Crowd Sourced Bathymetry'
        self.source = 'NOAA/NCEI'
        self.date = '2024'
        self.data_type = 'Bathymetric Soundings'
        self.resolution = '<10m to several kilometers'
        self.hdatum = 'WGS84'
        self.vdatum = 'MSL'
        self.url = self._csb_data_url
        
        ## aws stuff
        #self._bt_bucket = 'noaa-dcdb-bathymetry-pds'
        #self.s3 = boto3.client('s3', aws_access_key_id='', aws_secret_access_key='')
        #self.s3._request_signer.sign = (lambda *args, **kwargs: None)

        
    def run(self):
        """Run the CSB fetches module"""
        
        if self.region is None:
            return([])

        _data = {
            'where': self.where,
            'outFields': '*',
            'geometry': self.region.format('bbox'),
            'inSR':4326,
            'outSR':4326,
            'f':'pjson',
            'returnGeometry':'False'
        }
        _req = fetches.Fetch(
            self._csb_query_url, verbose=self.verbose
        ).fetch_req(params=_data)

        if _req is not None:
            features = _req.json()
            if 'features' in features.keys():
                for feature in features['features']:
                    if self.index:
                        print(json.dumps(feature['attributes'], indent=4))
                    else:
                        _name = feature['attributes']['NAME']
                        _year = _name[:4]
                        _dir_a = _name[4:6]
                        _dir_b = _name[6:8]
                        _csv_fn = '{}_pointData.csv'.format(_name[:-7])
                        link = '{0}/csb/csv/{1}/{2}/{3}/{4}'.format(
                            self._csb_data_url, _year, _dir_a, _dir_b, _csv_fn
                        )
                        if link is None:
                            continue

                        self.date = feature['attributes']['YEAR']
                        self.add_entry_to_results(link, _csv_fn, 'csb')

### End
