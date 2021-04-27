### gmrt.py - GMRT dataset
##
## Copyright (c) 2010 - 2021 CIRES Coastal DEM Team
##
## gmrt.py is part of CUDEM
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
## GMRT Fetch
##
## fetch extracts of the GMRT. - Global Extents
## https://www.gmrt.org/index.php
##
## =============================================================================
class GMRT(f_utils.FetchModule):
    '''Fetch raster data from the GMRT'''
    
    def __init__(self, res='max', fmt='geotiff', **kwargs):
        super().__init__(**kwargs) 
        #if self.verbose: utils.echo_msg('loading GMRT fetch module...')
        self._gmrt_grid_url = "https://www.gmrt.org:443/services/GridServer?"
        self._gmrt_grid_urls_url = "https://www.gmrt.org:443/services/GridServer/urls?"
        self._gmrt_grid_metadata_url = "https://www.gmrt.org/services/GridServer/metadata?"        
        self._outdir = os.path.join(os.getcwd(), 'gmrt')
        self.name = 'gmrt'
        self.res = res
        self.fmt = fmt
        
    def run(self):
        '''Run the GMRT fetching module'''

        if self.region is None: return([])
        self.data = {
            'north':self.region.ymax,
            'west':self.region.xmin,
            'south':self.region.ymin,
            'east':self.region.xmax,
            'mformat':'json',
            'resolution':self.res,
            'format':self.fmt,
        }
        
        req = f_utils.Fetch(self._gmrt_grid_urls_url).fetch_req(params=self.data, tries=10, timeout=2)
        if req is not None:
            gmrt_urls = req.json()
            for url in gmrt_urls:
                opts = {}
                for url_opt in url.split('?')[1].split('&'):
                    opt_kp = url_opt.split('=')
                    opts[opt_kp[0]] = opt_kp[1]
                url_region = regions.Region().from_list([float(opts['west']), float(opts['east']), float(opts['south']), float(opts['north'])])
                outf = 'gmrt_{}_{}.tif'.format(opts['layer'], url_region.format('fn'))
                self.results.append([url, outf, 'gmrt'])
        return(self)

    def yield_xyz(self, entry):
        src_data = 'gmrt_tmp.tif'
        if f_utils.Fetch(entry[0], callback=self.callback, verbose=self.verbose).fetch_file(src_data) == 0:
            gmrt_ds = datasets.RasterFile(fn=src_data, data_format=200, epsg=4326, warp=self.warp,
                                          name=src_data, src_region=self.region, verbose=self.verbose)
            for xyz in gmrt_ds.yield_xyz():
                yield(xyz)
        else: utils.echo_error_msg('failed to fetch remote file, {}...'.format(src_data))
        utils.remove_glob('{}*'.format(src_data))
    
### End
