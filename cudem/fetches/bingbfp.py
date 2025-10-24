### bingbfp.py
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

import os
import mercantile
import csv
from cudem import utils
from cudem.fetches import fetches

## BING Building Footprints
class BingBFP(fetches.FetchModule):
    """Bing Building Footprints

    https://github.com/microsoft/GlobalMLBuildingFootprints
    """
    
    def __init__(self, **kwargs):
        super().__init__(name='bingbfp', **kwargs)

        ## The various BING-BFP URLs
        self._bing_bfp_csv = ('https://minedbuildings.z5.web.core.windows.net/'
                              'global-buildings/dataset-links.csv')

        ## Set the user-agent and referer
        self.headers = {
            'User-Agent': ('Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:89.0) '
                           'Gecko/20100101 Firefox/89.0'),
            'referer': 'https://lz4.overpass-api.de/'
        }

        
    def run(self):
        """Run the Bing BFP fetches module"""
        
        if self.region is None:
            return([])

        bbox = self.region.format('bbox')
        quad_keys = set()
        for tile in list(
                mercantile.tiles(
                    self.region.xmin, self.region.ymin,
                    self.region.xmax, self.region.ymax, zooms=9
                )
        ):
            quad_keys.add(int(mercantile.quadkey(tile)))
            
        quad_keys = list(quad_keys)
        #utils.echo_msg('The input area spans {} tiles: {}'.format(len(quad_keys), quad_keys))
        #utils.echo_msg('The input area spans {} tiles.'.format(len(quad_keys)))
        bing_csv = os.path.join(self._outdir, os.path.basename(self._bing_bfp_csv))
        try:
            status = fetches.Fetch(self._bing_bfp_csv, verbose=self.verbose).fetch_file(bing_csv)
        except:
            status = -1
            
        if status == 0 and os.path.exists(bing_csv):
            with open(bing_csv, mode='r') as bc:
                reader = csv.reader(bc)
                next(reader)
                bd = [[row[2], row[1], row[0]] for row in reader if int(row[1]) in quad_keys]

            [self.add_entry_to_results(
                line[0], '{}_{}_{}'.format(line[2], line[1], os.path.basename(line[0])), 'bing'
            ) for line in bd]
        else:
            utils.echo_error_msg('could not fetch BING dataset-links.csv')

### End
