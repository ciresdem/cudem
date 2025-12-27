### synbath.py
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
##
## synbath.py is part of CUDEM
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
## Fetch UCSD SynBath global bathymetry data.
##
### Code:

from typing import Optional
from cudem import utils
from cudem.fetches import fetches

## ==============================================
## Constants
## ==============================================
SYNBATH_URL = 'https://topex.ucsd.edu/pub/synbath/SYNBATH_V2.0.nc'

## ==============================================
## SynBath Module
## ==============================================
class SynBath(fetches.FetchModule):
    """UCSD SynBath dataset

    Fetch SynBath (Synthetic Bathymetry) global grid.
    Currently only fetches entire grid. Subset in dlim, or elsewhere.
    
    https://topex.ucsd.edu/pub/synbath/
    https://topex.ucsd.edu/pub/synbath/SYNBATH_publication.pdf

    Configuration Example:
    < synbath:upper_limit=None:lower_limit=None >
    """
    
    def __init__(self, upper_limit: Optional[float] = None, lower_limit: Optional[float] = None, **kwargs):
        super().__init__(name='synbath', **kwargs)
        
        ## Set the fetching region and restrict by z-region if desired.
        ## Note: This affects downstream processing, not the fetch itself (which grabs the whole file).
        if self.region is not None:
            self.synbath_region = self.region.copy()
            self.synbath_region.zmax = utils.float_or(upper_limit)
            self.synbath_region.zmin = utils.float_or(lower_limit)
        else:
            self.synbath_region = None

        ## Data format 200 is typically GDAL/NetCDF
        self.data_format = 200
        self.src_srs = 'epsg:4326+3855'

        
    def run(self):
        """Run the SynBath fetching module."""
        
        self.add_entry_to_results(
            SYNBATH_URL,
            'SYNBATH_V2_0.nc',
            'synbath'
        )
        
        return self

    
### End
