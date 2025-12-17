### chs.py
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

import lxml.etree
from cudem import regions
from cudem.fetches import fetches

## CHS - Canada Hydro
class GBA(fetches.FetchModule):
    
    def __init__(self,  **kwargs):
        super().__init__(name='gba', **kwargs)

        ## The various CHS URLs
        self._gba_url = 'https://tubvsig-so2sat-vm1.srv.mwn.de/geoserver/ows?'

        
    def run(self):
        """Run the CHS fetching module"""

        if self.region is None: return([])
        # _data = {
        #     'request': 'DescribeCoverage',
        #     'version': '2.0.1',
        #     'CoverageID': 'nonna__NONNA {} Coverage'.format(self.datatype),
        #     'service': 'WCS'
        # }
        # _req = fetches.Fetch(self._chs_url).fetch_req(params=_data)
        # _results = lxml.etree.fromstring(_req.text.encode('utf-8'))        
        # g_env = _results.findall(
        #     './/{http://www.opengis.net/gml/3.2}GridEnvelope',
        #     namespaces=fetches.namespaces
        # )[0]
        # hl = [float(x) for x in g_env.find(
        #     '{http://www.opengis.net/gml/3.2}high'
        # ).text.split()]
        # g_bbox = _results.findall(
        #     './/{http://www.opengis.net/gml/3.2}Envelope'
        # )[0]
        # lc = [float(x) for x in g_bbox.find(
        #     '{http://www.opengis.net/gml/3.2}lowerCorner'
        # ).text.split()]
        # uc = [float(x) for x in g_bbox.find(
        #     '{http://www.opengis.net/gml/3.2}upperCorner'
        # ).text.split()]
        # ds_region = regions.Region().from_list([lc[1], uc[1], lc[0], uc[0]])
        # resx = (uc[1] - lc[1]) / hl[0]
        # resy = (uc[0] - lc[0]) / hl[1]
        # if regions.regions_intersect_ogr_p(self.region, ds_region):
        # _gba_data = {
        #     'request': 'GetCoverage',
        #     'version': '2.0.1',
        #     #'CoverageID': 'global3D__bf_120m',
        #     #'CoverageID': 'global3D__global_volume_distribution',
        #     'CoverageID': 'global3D__volume_preview',
        #     'service': 'WCS',
        #     'subset': ['Long({},{})'.format(self.region.xmin, self.region.xmax),
        #                'Lat({},{})'.format(self.region.ymin, self.region.ymax)],
        #     'subsettingcrs': 'http://www.opengis.net/def/crs/EPSG/0/4326',
        #     'outputcrs': 'http://www.opengis.net/def/crs/EPSG/0/4326'
        # }
        # _gba_req = fetches.Fetch(self._gba_url).fetch_req(params=_gba_data)
        # outf = 'gba_{}.tif'.format(self.region.format('fn'))
        # self.add_entry_to_results(_gba_req.url, outf, 'gba')

        _gba_data = {
            'request': 'GetFeature',
            'version': '2.0.0',
            #'CoverageID': 'global3D__bf_120m',
            #'CoverageID': 'global3D__global_volume_distribution',
            'typenames': 'lod1_global',
            'service': 'WFS',
            'bbox': '{},http://www.opengis.net/def/crs/EPSG/0/4326'.format(self.region.format('osm_bbox')),
            'outputformat': 'json',
            #'subset': ['Long({},{})'.format(self.region.xmin, self.region.xmax),
            #           'Lat({},{})'.format(self.region.ymin, self.region.ymax)],
            #'subsettingcrs': 'http://www.opengis.net/def/crs/EPSG/0/4326',
            'targetcrs': 'http://www.opengis.net/def/crs/EPSG/0/4326'
        }
        _gba_req = fetches.Fetch(self._gba_url).fetch_req(params=_gba_data)
        #print(_gba_req.url)
        #print(_gba_req.text)
        outf = 'gba_{}.geojson'.format(self.region.format('fn'))
        self.add_entry_to_results(_gba_req.url, outf, 'gba')

        
        return(self)

### End
