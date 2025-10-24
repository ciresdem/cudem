### nswtb.py
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

from cudem.fetches import fetches

class NSW_TB(fetches.FetchModule):
    """New South Wales Topo-Bathy DEM

    **testing**

    < nsw_tb:where=None:layer=0:index=False >
    """
    
    def __init__(self, where = '1=1', layer = 0, index = False, **kwargs):
        super().__init__(name='csb', **kwargs)
        self.where = where
        self.index = index
        self.src_srs = None

        ## The various NSW_TB URLs
        self._nsw_map_server = ('https://mapprod2.environment.nsw.gov.au/arcgis/'
                                'rest/services/Coastal_Marine/'
                                'NSW_Marine_Lidar_Bathymetry_Data_2018/MapServer')
        self._nsw_query_url = '{0}/{1}/query?'.format(self._nsw_map_server, layer)

        
    def run(self):
        """Run the NSW_TB fetching module"""
        
        if self.region is None:
            return([])

        _data = {
            'where': self.where,
            'outFields': '*',
            'geometry': self.region.format('bbox'),
            'inSR':4326,
            'outSR':4326,
            'f':'geojson',
            'returnGeometry':'True',
            'geometryType':'esriGeometryEnvelope',
            'spatialRel':'esriSpatialRelIntersects'
        }
        _req = fetches.Fetch(
            self._nsw_query_url, verbose=self.verbose
        ).fetch_req(params=_data)
        _geojson_fn = 'nsw_{}_contours.geojson'.format(self.region.format('fn'))
        if _req is not None:
            #print(_req.url)
            features = _req.json()
            #print(features)
            if 'features' in features.keys():
                self.add_entry_to_results(_req.url, _geojson_fn, 'nsw_contours')

### End
