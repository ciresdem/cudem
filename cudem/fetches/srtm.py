### srtm.py - SRTM fetch
##
## Copyright (c) 2010 - 2021 CIRES Coastal DEM Team
##
## srtm.py is part of CUDEM
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
## SRTM Plus
##
## Fetch srtm+ data
## https://topex.ucsd.edu/WWW_html/srtm15_plus.html
## http://topex.ucsd.edu/sandwell/publications/180_Tozer_SRTM15+.pdf
## https://topex.ucsd.edu/pub/srtm15_plus/
## https://topex.ucsd.edu/pub/srtm15_plus/SRTM15_V2.3.nc
##
### Code:

import os

from cudem import utils
from cudem import regions
from cudem import datasets

import cudem.fetches.utils as f_utils

class SRTMPlus(f_utils.FetchModule):
    '''Fetch SRTM+ data'''
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs) 
        self._srtm_url = 'https://topex.ucsd.edu/cgi-bin/get_srtm15.cgi'
        self._outdir = os.path.join(os.getcwd(), 'srtm_plus')
        self.name = 'srtm_plus'

    def run(self):
        '''Run the SRTM fetching module.'''
        
        if self.region is None: return([])

        self.data = {
            'north':self.region.ymax,
            'west':self.region.xmin,
            'south':self.region.ymin,
            'east':self.region.xmax,
        }

        _req = f_utils.Fetch(self._srtm_url, verify=False).fetch_req(params=self.data)
        if _req is not None:
            outf = 'srtm_{}.xyz'.format(self.region.format('fn'))
            self.results.append([_req.url, outf, 'srtm'])
            
        return(self)

    def yield_xyz(self, entry):
        src_data = 'srtm_plus_tmp.xyz'
        if f_utils.Fetch(
                entry[0], callback=self.callback, verbose=self.verbose, verify=False
        ).fetch_file(src_data) == 0:
            _ds = datasets.XYZFile(
                fn=src_data,
                data_format=168,
                skip=1,
                epsg=4326,
                warp=self.warp,
                name=src_data,
                src_region=self.region,
                verbose=self.verbose,
                remote=True
            )
            for xyz in _ds.yield_xyz():
                yield(xyz)
                
        else:
            utils.echo_error_msg('failed to fetch remote file, {}...'.format(src_data))
            
        utils.remove_glob('{}*'.format(src_data))
### End
