### flatten.py
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
##
## flatten.py is part of CUDEM
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
### Code:

from cudem import gdalfun
from cudem.waffles.waffles import Waffle

class WafflesFlatten(Waffle):
    """Stack the data into a DEM and then hydro-flatten all the 
    void areas.

    specify 'size_threshold' to only flatten voids above threshold.

    -----------
    Parameters:
    
    min_count=[val] - only retain data cells if they contain `min_count` 
                      overlapping data
    size_threshold=[val] - the minimum size void to flatten (in cells)
    """
    
    ## todo: add parameters for specifying outputs...
    def __init__(self, min_count=None, size_threshold=1, **kwargs):
        super().__init__(**kwargs)
        self.min_count = min_count
        self.size_threshold = size_threshold

        
    def run(self):
        gdalfun.cudem_flatten_no_data_zones(
            self.stack,
            dst_dem=self.fn,
            band=1,
            size_threshold=self.size_threshold
        )
        
        return(self)

### End
