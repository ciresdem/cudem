### mrdem.py
##
## Copyright (c) 2023 - 2026 Regents of the University of Colorado
##
## mrdem.py is part of CUDEM
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
## Fetch Multi-Resolution Digital Elevation Model (MRDEM) data (Canada).
##
### Code:

import os
from typing import Optional
from cudem.fetches import fetches

## ==============================================
## Constants
## ==============================================
MRDEM_BASE_URL = 'https://canelevation-dem.s3.ca-central-1.amazonaws.com/mrdem-30'
MRDEM_DTM_URL = f'{MRDEM_BASE_URL}/mrdem-30-dtm.vrt'
MRDEM_DSM_URL = f'{MRDEM_BASE_URL}/mrdem-30-dsm.vrt'

## ==============================================
## MRDEM Module
## ==============================================
class MRDEM(fetches.FetchModule):
    """Multi-Resolution Digital Elevation Model (MRDEM)
    
    Fetches VRT pointers for the Canadian MRDEM-30 dataset.
    
    Args:
        datatype (str): 'dtm' (Digital Terrain Model) or 'dsm' (Digital Surface Model).
                        Defaults to 'dtm'.

    Configuration Example:
    < mrdem:datatype=dtm >
    """
    
    def __init__(self, datatype: str = 'dtm', **kwargs):
        super().__init__(name='mrdem', **kwargs)
        self.datatype = datatype.lower() if datatype else 'dtm'

        
    def run(self):
        """Run the MRDEM fetching module."""
        
        ## Determine which URL to use based on datatype
        if self.datatype == 'dsm':
            url = MRDEM_DSM_URL
        else:
            url = MRDEM_DTM_URL

        self.add_entry_to_results(
            url,
            os.path.basename(url),
            'vrt'
        )
        
        return self

### End
