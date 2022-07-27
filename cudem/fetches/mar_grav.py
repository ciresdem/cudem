### mar_grav.py - mar_grav fetch
##
## Copyright (c) 2010 - 2021 CIRES Coastal DEM Team
##
## mar_grav.py is part of CUDEM
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
## mar_grav - Sattelite Altimetry Topography from Scripps
##
## https://topex.ucsd.edu/WWW_html/mar_grav.html
## ftp://topex.ucsd.edu/pub/global_grav_1min/
## https://topex.ucsd.edu/marine_grav/explore_grav.html
## https://topex.ucsd.edu/marine_grav/white_paper.pdf
##
### Code:

import os

from cudem import utils
from cudem import regions
from cudem import datasets

import cudem.fetches.utils as f_utils

class MarGrav(f_utils.FetchModule):
    '''Fetch mar_grav sattelite altimetry topography'''
    
    def __init__(self, mag=1, **kwargs):
        super().__init__(**kwargs) 
        self._mar_grav_url = 'https://topex.ucsd.edu/cgi-bin/get_data.cgi'
        self._outdir = os.path.join(os.getcwd(), 'mar_grav')
        self.name = 'mar_grav'
        self.mag = mag if mag == 1 else 0.1

    def run(self):
        '''Run the mar_grav fetching module.'''
        
        if self.region is None:
            return([])
        
        _data = {
            'north':self.region.ymax,
            'west':self.region.xmin,
            'south':self.region.ymin,
            'east':self.region.xmax,
            'mag':self.mag,
        }
        _req = f_utils.Fetch(self._mar_grav_url, verify=False).fetch_req(params=_data)
        if _req is not None:
            outf = 'mar_grav_{}.xyz'.format(self.region.format('fn'))
            self.results.append([_req.url, os.path.join(self._outdir, outf), 'mar_grav'])
            
        return(self)
        
    def yield_xyz(self, entry):
        src_data = 'mar_grav_tmp.xyz'
        if f_utils.Fetch(
                entry[0], callback=self.callback, verbose=self.verbose, verify=False
        ).fetch_file(src_data) == 0:
            _ds = datasets.XYZFile(
                fn=src_data,
                data_format=168,
                skip=1,
                x_offset=-360,
                src_srs='epsg:4326+3855',
                dst_srs=self.dst_srs,
                src_region=self.region,
                x_inc=self.x_inc,
                y_inc=self.y_inc,
                verbose=self.verbose,
                remote=True
            )
            for xyz in _ds.yield_xyz():
                yield(xyz)
                
        else:
            utils.echo_error_msg('failed to fetch remote file, {}...'.format(src_data))
            
        utils.remove_glob('{}*'.format(src_data))                

### End
