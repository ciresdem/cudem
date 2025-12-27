### gedtm30.py
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
##
## gedtm30.py is part of CUDEM
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
### Commentary:
##
## Fetch Global 1-Arc-Second Digital Terrain Model (GEDTM30) data.
##
### Code:

import os
import csv
from io import StringIO
from typing import Optional, List
from cudem import utils
from cudem.fetches import fetches

## ==============================================
## Constants
## ==============================================
GEDTM30_COG_LIST_URL = 'https://raw.githubusercontent.com/openlandmap/GEDTM30/refs/heads/main/metadata/cog_list.csv'
GEDTM30_BASE_URL = 'https://github.com/openlandmap/GEDTM30'

## ==============================================
## GEDTM30 Module
## ==============================================
class GEDTM30(fetches.FetchModule):
    """Global 1-Arc-Second Digital Terrain Model
    
    Fetch data from OpenLandMap's GEDTM30 repository.
    
    https://github.com/openlandmap/GEDTM30

    Configuration Example:
    < gedtm30 >
    """

    def __init__(self, products: str = 'Ensemble Digital Terrain Model', **kwargs):
        super().__init__(name='gedtm30', **kwargs)
        if products is not None:
            self.products = products.split('/')
        else:
            self.products = ['Ensemble Digital Terrain Model']
        
        self._gedtm30_url_list = GEDTM30_COG_LIST_URL

        
    def run(self):
        """Run the GEDTM30 fetching module."""
        
        ## Fetch COG list CSV
        cog_req = fetches.Fetch(
            self._gedtm30_url_list, 
            verbose=self.verbose
        ).fetch_req()
        
        if cog_req is None:
            utils.echo_error_msg("Failed to fetch GEDTM30 COG list.")
            return self

        try:
            ## Parse CSV from response text
            csv_content = StringIO(cog_req.text)
            reader = csv.reader(csv_content)
            
            ## Skip header
            header = next(reader, None)
            if not header:
                return self

            ## Filter rows based on product selection
            ## CSV structure: Product Name (col 0), URL (last col)
            for row in reader:
                if not row:
                    continue
                    
                product_name = row[0]
                url = row[-1]
                
                if product_name in self.products:
                    self.add_entry_to_results(
                        url,
                        os.path.basename(url),
                        f'gedtm30 - {product_name}'
                    )
                    
        except Exception as e:
            utils.echo_error_msg(f"Error parsing GEDTM30 CSV: {e}")

        return self

### End
