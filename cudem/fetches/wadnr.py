### wadnr.py
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
from cudem import regions
from cudem.fetches import fetches

class waDNR(fetches.FetchModule):
    """Washington State Department of Natural Resources lidar data.
    """
    
    def __init__(self, ids=None, **kwargs):
        super().__init__(name='waDNR', **kwargs) 

        ## The various urls to use for GMRT
        self._wa_dnr_url = "https://lidarportal.dnr.wa.gov/download?"
        self._wa_dnr_rest = ('https://lidarportal.dnr.wa.gov/arcgis/rest'
                             '/services/lidar/wadnr_hillshade/MapServer')
        self._wa_dnr_ids = ('https://lidarportal.dnr.wa.gov/arcgis/rest/'
                            'services/lidar/wadnr_hillshade/MapServer/identify')
        self._wa_dnr_layers = ('https://lidarportal.dnr.wa.gov/arcgis/rest/'
                               'services/lidar/wadnr_hillshade/MapServer/layers?f=pjson')
        self._wa_dnr_cap = ('https://lidarportal.dnr.wa.gov/arcgis/services/'
                            'lidar/wadnr_hillshade/MapServer/WmsServer?'
                            'service=WMS&request=GetCapabilities')
        self._wa_dnr_map = ('https://lidarportal.dnr.wa.gov/arcgis/services/'
                            'lidar/wadnr_hillshade/MapServer/WmsServer?'
                            'service=WMS&request=GetMap')
        
        ## for dlim, data format is -2 for a zip file, projections vary
        self.data_format = -2
        self.src_srs = None
        self.title = 'Washington DNR Lidar'
        self.source = 'Washington DNR'
        self.date = None
        self.data_type = 'lidar'
        self.resolution = None
        self.hdatum = None
        self.vdatum = None
        self.url = self._wa_dnr_url
        
        ## Firefox on windows for this one.
        self.headers = {
            'User-Agent': ('Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:89.0) '
                           'Gecko/20100101 Firefox/89.0')
        }

        if ids is not None:
            self.ids = [int(x) for x in ids.split('/')]
        else:
            self.ids = []

            
    def find_ids(self):
        self.data = {
            'geometry': self.region.format('bbox'),
            'format': 'json',
        }

        layers_req = fetches.Fetch(
            self._wa_dnr_layers
        ).fetch_req(
            tries=10, timeout=2
        )

        layer_srs = 'epsg:3857'
        layers_in_region = []
        if layers_req is not None:
            layers_json = layers_req.json()
            for layer in layers_json['layers']:
                layer_name = layer['name']
                layer_id = layer['id']
                layer_sublayers = layer['subLayers']
                    
                if 'parentLayer' in layer.keys():
                    layer_parent_layer = layer['parentLayer']
                else:
                    layer_parent_layer = None
                    
                layer_type = layer['type']
                layer_region = regions.Region(src_srs='epsg:3857').from_list(
                    [layer['extent']['xmin'], layer['extent']['xmax'],
                     layer['extent']['ymin'], layer['extent']['ymax']]
                ).warp('epsg:4326')

                if len(layer_sublayers) > 0:
                    if regions.regions_intersect_ogr_p(layer_region, self.region):
                        layers_in_region.append(
                            [
                                layer_name,
                                min([int(x['name'][:-1])-1 for x in layer_sublayers])
                            ]
                        )
                        
        if len(self.ids) > 0:
            layers_in_region = [x for x in layers_in_region if x[1] in self.ids]
            
        return(layers_in_region)

    
    def run(self):
        '''Run the GMRT fetching module'''

        if self.region is None:
            return([])

        layers_in_region = self.find_ids()        
        data = {}
        for l in layers_in_region:
            if 'ids' in data.keys():
                data['ids'].append(l[1])
            else:
                data['ids'] = [int(l[1])]

            data['geojson'] = json.dumps({
                'type': 'Polygon',
                'coordinates': [self.region.export_as_polygon()]
            })
            
            data_req = fetches.Fetch(
                self._wa_dnr_url
            ).fetch_req(
                params=data, tries=10, timeout=2
            )

            if data_req is not None and data_req.status_code == 200:
                self.add_entry_to_results(
                    data_req.url, '{}_{}.zip'.format(l[0], l[1]), 'wa_dnr'
                )
                
            data = {}
                
        return(self)

### End
