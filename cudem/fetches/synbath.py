### synbath.py
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

from cudem import utils
from cudem.fetches import fetches

class SynBath(fetches.FetchModule):
    """UCSD SynBath dataset

    Currently only fetches entire grid. Subset in dlim, or elsewhere.
    
    https://topex.ucsd.edu/pub/synbath/
    https://topex.ucsd.edu/pub/synbath/SYNBATH_publication.pdf

    < gebco:upper_limit=None:lower_limit=None >
    """
    
    def __init__(self, upper_limit=None, lower_limit=None, **kwargs):
        super().__init__(name='synbath', **kwargs)
        
        ## various gebco URLs
        self._synbath_url_1_2 = 'https://topex.ucsd.edu/pub/synbath/SYNBATH_V1.2.nc'
        self._synbath_url = 'https://topex.ucsd.edu/pub/synbath/SYNBATH_V2.0.nc'

        ## set the fetching region, restrict by z-region if desired.
        self.synbath_region = self.region.copy()
        self.synbath_region.zmax = utils.float_or(upper_limit)
        self.synbath_region.zmin = utils.float_or(lower_limit)

        ## for dlim
        self.data_format = 200
        self.src_srs = 'epsg:4326+3855'

        
    def run(self):
        """Run the SynBath fetching module"""

        self.add_entry_to_results(
            self._synbath_url,
            'SYNBATH_V2_0.nc',
            'synbath'
        )

### End
