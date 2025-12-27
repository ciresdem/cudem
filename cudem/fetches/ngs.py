### ngs.py
##
## Copyright (c) 2019 - 2025 Regents of the University of Colorado
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
## Fetch NGS Monuments from NOAA.
##
### Code:

from urllib.parse import urlencode
from cudem import utils
from cudem.fetches import fetches

## ==============================================
## Constants
## ==============================================
NGS_SEARCH_URL = 'http://geodesy.noaa.gov/api/nde/bounds?'

## ==============================================
## NGS Module
## ==============================================
class NGS(fetches.FetchModule):
    """NGS Monuments
    
    NGS provides Information about survey marks (including bench marks) 
    in text datasheets or in GIS shapefiles. 

    Note: Some survey markers installed by other organizations may not 
    be available through NGS.

    Fetch NGS monuments from NOAA
    
    http://geodesy.noaa.gov/

    Configuration Example:
    < ngs:datum=geoidHt >
    """

    def __init__(self, datum: str = 'geoidHt', **kwargs):
        super().__init__(name='ngs', **kwargs)
        
        valid_datums = ['orthoHt', 'geoidHt', 'z', 'ellipHeight']
        if datum not in valid_datums:
            utils.echo_warning_msg(
                f'Could not parse datum {datum}, falling back to geoidHt'
            )
            self.datum = 'geoidHt'
        else:
            self.datum = datum

        self.src_srs = 'epsg:4326'

        
    def run(self):
        """Run the NGS (monuments) fetching module."""
        
        if self.region is None:
            return []
        
        ## Construct Query Parameters
        params = {
            'maxlon': self.region.xmax,
            'minlon': self.region.xmin,
            'maxlat': self.region.ymax,
            'minlat': self.region.ymin
        }
        
        ## Generate URL
        query_string = urlencode(params)
        full_url = f"{NGS_SEARCH_URL}{query_string}"
        
        ## Output Filename
        out_fn = f"ngs_results_{self.region.format('fn')}.json"
        
        self.add_entry_to_results(
            full_url,
            out_fn,
            'ngs'
        )
            
        return self

    
### End
