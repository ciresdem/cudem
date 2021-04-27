### ngs.py - SRTM fetch
##
## Copyright (c) 2010 - 2021 CIRES Coastal DEM Team
##
## ngs.py is part of CUDEM
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
import json
from cudem import utils
from cudem import regions
from cudem import datasets
from cudem import xyzfun
import cudem.fetches.utils as f_utils

## =============================================================================
##
## National Geodetic Survey (NGS)
##
## Fetch NGS monuments from NGS - US Only
##
## =============================================================================
class NGS(f_utils.FetchModule):
    """Fetch NGS monuments from NOAA"""
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._ngs_search_url = 'http://geodesy.noaa.gov/api/nde/bounds?'
        self._outdir = os.path.join(os.getcwd(), 'ngs')
        self.name = 'ngs'

    def run(self, csv = False):
        """Run the NGS (monuments) fetching module."""
        
        if self.region is None: return([])
        self.data = {
            'maxlon':self.region.xmax,
            'minlon':self.region.xmin,
            'maxlat':self.region.ymax,
            'minlat':self.region.ymin,
        }

        _req = f_utils.Fetch(self._ngs_search_url).fetch_req(params=self.data)
        if _req is not None:
            self.results.append([_req.url, 'ngs_results_{}.json'.format(self.region.format('fn')), 'ngs'])
        return(self)

    def yield_xyz(self, entry):
        """process ngs monuments"""
        
        src_data = 'ngs_tmp.json'
        src_csv = 'ngs_tmp.csv'
        
        if f_utils.Fetch(entry[0], callback=self.callback, verbose=self.verbose).fetch_file(src_data) == 0:
            with open(src_data, 'r') as json_file:
                r = json.load(json_file)
            
            if len(r) > 0:
                for row in r:
                    yield(xyzfun.XYZPoint().from_list([float(row['lon']), float(row['lat']), float(row['geoidHt'])]))

        else: utils.echo_error_msg('failed to fetch remote file, {}...'.format(src_data))
        utils.remove_glob('{}*'.format(src_data))

### End
