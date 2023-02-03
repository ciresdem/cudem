### osm.py - open street map fetch
##
## Copyright (c) 2010 - 2023 Regents of the University of Colorado
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
##
## Opens Street Map data (OSM)
##
## Fetch various datasets from OSM/Overpass
## https://wiki.openstreetmap.org/wiki/Overpass_API
##
### Code:

import os

from cudem import utils
from cudem import regions
from cudem import datasets
from cudem import xyzfun

import cudem.fetches.utils as f_utils

class OpenStreetMap(f_utils.FetchModule):
    """Fetch OSM data"""
    
    def __init__(self, q=None, fmt='osm', planet=False, chunks=True, min_length=None, **kwargs):
        super().__init__(name='osm', **kwargs)
        self._osm_api = 'https://lz4.overpass-api.de/api/interpreter'
        self._osm_api2 = 'https://overpass.kumi.systems/api/interpreter'
        self._osm_api3 = 'https://overpass.openstreetmap.fr/api/interpreter'
        self._osm_planet_bz2 = 'https://ftpmirror.your.org/pub/openstreetmap/planet/planet-latest.osm.bz2'
        self._osm_planet = 'https://ftpmirror.your.org/pub/openstreetmap/pbf/planet-latest.osm.pbf'
        self._osm_continents = 'https://download.geofabrik.de/'
        self.q = q
        self.fmt = fmt
        self.planet = planet
        self.chunks = chunks

        self.h = ''
        
        if self.q == 'buildings':
            #self.h = '[maxsize:2000000000]'
            self.h = '[timeout:3600]'
            self.q = '''
            (way["building"]{};
            relation["building"]["type"="multipolygon"];
            );
            (._;>;);
            out meta;
            '''.format('(if: length() > {})'.format(min_length) if min_length is not None else '')
        
        self.headers = { 'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:89.0) Gecko/20100101 Firefox/89.0',
                         'referer': 'https://lz4.overpass-api.de/' }
        
    def run(self):
        if self.region is None:
            return([])

        #if self.planet:
        #    self.results.append([self._osm_planet, os.path.join(self._outdir, 'planet-latest.osm.bz2'), 'bz2'])
        
        ## fetch whole planet
        if self.planet:
            self.results.append([self._osm_planet, os.path.join(self._outdir, 'planet-latest.osm.pbf'), 'pbf'])

        ## fetch in chunks
        elif self.chunks:
            x_delta = self.region.xmax - self.region.xmin
            y_delta = self.region.ymax - self.region.ymin
            incs = self.region.increments(1000,1000)

            ## break up the requests into .05 degree chunks for
            ## better usage of the OSM API
            if x_delta > .25 or y_delta > .25:
                xcount, ycount, gt = self.region.geo_transform(x_inc=incs[0], y_inc=incs[1])
                if x_delta >= y_delta:
                    n_chunk = int(xcount*(.25/x_delta))
                elif y_delta > x_delta:
                    n_chunk = int(ycount*(.25/y_delta))
            else:
                n_chunk = None

            these_regions = self.region.chunk(incs[0], n_chunk=n_chunk)
            utils.echo_msg('chunking OSM request into {} regions'.format(len(these_regions)))
            
            for this_region in these_regions:
                c_bbox = this_region.format('osm_bbox')
                out_fn = 'osm_{}'.format(this_region.format('fn_full'))
                osm_q_bbox  = '''
                {1}{2}[bbox:{0}];'''.format(c_bbox, '[out:{}]'.format(self.fmt) if self.fmt != 'osm' else '', self.h)

                osm_q = '''
                (node;
                <;
                >;
                );
                out meta;
                '''
                
                osm_q_ = osm_q_bbox + (osm_q if self.q is None else self.q)

                #utils.echo_msg('using query: {}'.format(osm_q_))
                osm_data = f_utils.urlencode({'data': osm_q_})
                osm_data_url = self._osm_api + '?' + osm_data
                
                self.results.append([osm_data_url, os.path.join(self._outdir, '{}.{}'.format(out_fn, self.fmt)), 'osm'])

        else:
            c_bbox = self.region.format('osm_bbox')
            out_fn = 'osm_{}'.format(self.region.format('fn_full'))
            osm_q_bbox  = '''
            {1}[bbox:{0}];'''.format(c_bbox, '[out:{}]'.format(self.fmt) if self.fmt != 'osm' else '')

            osm_q = '''
            (node;
            <;
            >;
            );
            out meta;
            '''

            osm_q_ = osm_q_bbox + (osm_q if self.q is None else self.q)
            osm_data = f_utils.urlencode({'data': osm_q_})
            osm_data_url = self._osm_api + '?' + osm_data            
            self.results.append([osm_data_url, os.path.join(self._outdir, '{}.{}'.format(out_fn, self.fmt)), 'osm'])
                
                
### End
