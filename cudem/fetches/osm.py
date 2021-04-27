### osm.py - open street map fetch
##
## Copyright (c) 2010 - 2021 CIRES Coastal DEM Team
##
## osm.py is part of CUDEM
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
### Code:

import os
from cudem import utils
from cudem import regions
from cudem import datasets
import cudem.fetches.utils as f_utils

## =============================================================================
##
## Opens Street Map data (OSM)
##
## Fetch various datasets from OSM/Overpass
##
## =============================================================================
class OpenStreetMap(f_utils.FetchModule):
    """Fetch OSM data"""
    
    def __init__(self, osm_type = 'all', proc=False, **kwargs):
        super().__init__(**kwargs)
        self._osm_api = 'https://lz4.overpass-api.de/api/interpreter'
        self._outdir = os.path.join(os.getcwd(), 'osm')

        self.osm_types = {
            'highway': ['LINESTRING'],
            'waterway': ['LINESTRING'],
            'building': ['POLYGON'],
        }
        
    def run(self):
        '''Run the OSM fetching module.'''
        
        if self.region is None: return([])
        if osm_type == 'all':
            for key in self.osm_types.keys():
                self._fetch_results(osm_type = key, proc = proc)
        else:
            if osm_type in self.osm_types.keys():
                self._fetch_results(osm_type = osm_type, proc = proc)
                
        return(self)
        
    def _fetch_results(self, osm_type = 'highway', proc = False):
        c_bbox = self.region.format('osm_bbox')
        out_fn = 'osm_{}_{}'.format(osm_type, self.region.format('fn'))
        
        osm_q = '''
        [out:json];
        (
        way['{}']({});
        );
        out body;
        '''.format(osm_type, c_bbox)

        osm_data = f_utils.Fetch(self._osm_api).fetch_req(params={'data': osm_q}, timeout=3600)
        osm_json = osm_data.json()

        if self.verbose:
            utils.echo_msg('found {} {} elements'.format(len(osm_json['elements']), osm_type))
        
        if proc:
            if not os.path.exists(self._outdir):
                try:
                    os.makedirs(self._outdir)
                except: pass 

            dst_gmt = open('{}.gmt'.format(os.path.join(self._outdir, out_fn)), 'w')
            dst_gmt.write('# @VGMT1.0 @G{} @Nid\n'.format(self.osm_types[osm_type][0]))
            dst_gmt.write('# @Tdouble\n')
            dst_gmt.write('# @R{}\n'.format(self.region.format('str'))
            dst_gmt.write('# FEATURE_DATA\n')

        for el_n, el in enumerate(osm_json['elements']):
            utils.echo_msg_inline('fetching osm {} data [{}%]'.format(osm_type, float(el_n) / len(osm_json['elements']) * 100))
            if el['type'] == 'way':

                osm_node_q = '''
                [out:json];
                (
                node(id:{});
                );
                out body;
                '''.format(','.join([str(n) for n in el['nodes']]))

                node_data = f_utils.Fetch(self._osm_api).fetch_req(params={'data': osm_node_q})
                if node_data.status_code == 200:
                    self.results.append([node_data.url, '{}.json'.format(out_fn), 'osm'])
                    
                    if proc:
                        node_json = node_data.json()

                        dst_gmt.write('>\n# @D{}\n'.format(el['id']))

                        for node in el['nodes']:
                            for el_node in node_json['elements']:
                                if el_node['id'] == node:
                                    xyzfun.xyz_line([el_node['lon'], el_node['lat']], dst_gmt)
                                    
        if self.verbose: utils.echo_msg_inline('fetching osm {} data [OK]\n'.format(osm_type))
        
        if proc:
            dst_gmt.close()
            utils.run_cmd('ogr2ogr {}.shp {}.gmt'.format(os.path.join(self._outdir, out_fn), os.path.join(self._outdir, out_fn)), verbose=False)
        
### End
