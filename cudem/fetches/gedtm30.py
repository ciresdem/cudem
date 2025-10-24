### gedtm30.py
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
import csv
from io import StringIO
from cudem.fetches import fetches

class GEDTM30(fetches.FetchModule):
    """Global 1-Arc-Second Digital Terrain Model
    """

    def __init__(self, products='Ensemble Digital Terrain Model', **kwargs):
        super().__init__(name='gedtm30', **kwargs)
        if products is not None:
            self.products = products.split('/')
        else:
            self.products = ['Ensemble Digital Terrain Model']
        
        self._gedtm30_url_list = ('https://raw.githubusercontent.com/openlandmap/'
                                  'GEDTM30/refs/heads/main/metadata/cog_list.csv')
        self._gedtm30_url = 'https://github.com/openlandmap/GEDTM30'

        
    def run(self):
        ## fetch COG list
        cog_req = fetches.Fetch(
            self._gedtm30_url_list, verbose=self.verbose
        ).fetch_req()
        cog_url_list = None
        if cog_req is not None:
            cog_url_list = cog_req.text

        if cog_url_list is not None:
            csvfile = StringIO(cog_url_list)
            reader = csv.reader(csvfile)
            header = next(reader)
            #[print(row[0], row[-1]) for row in reader]
            urls = [[row[0], row[-1]] for row in reader if row[0] in self.products]

            for url in urls:
                self.add_entry_to_results(
                    url[1],
                    os.path.basename(url[1]),
                    'gedtm30 - {}'.format(url[0])
                )
                
### End
