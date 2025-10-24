### ngs.py
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

from cudem.fetches import fetches

## NGS - geodesic monuments
class NGS(fetches.FetchModule):
    """NGS Monuments
    
    NGS provides Information about survey marks (including bench marks) 
    in text datasheets or in GIS shapefiles. 

    Note some survey markers installed by other organizations may not 
    be available through NGS.

    Fetch NGS monuments from NOAA
    
    http://geodesy.noaa.gov/

    < ngs:datum=geoidHt >
    """

    def __init__(self, datum='geoidHt', **kwargs):
        super().__init__(name='ngs', **kwargs)
        if datum not in ['orthoHt', 'geoidHt', 'z', 'ellipHeight']:
            utils.echo_warning_msg(
                'could not parse {}, falling back to geoidHt'.format(datum)
            )
            self.datum = 'geoidHt'
        else:
            self.datum = datum

        ## The various NGS URLs
        self._ngs_search_url = 'http://geodesy.noaa.gov/api/nde/bounds?'

        ## for dlim
        self.src_srs = 'epsg:4326'

        
    def run(self, csv=False):
        """Run the NGS (monuments) fetching module."""
        
        if self.region is None:
            return([])
        
        _data = {
            'maxlon':self.region.xmax,
            'minlon':self.region.xmin,
            'maxlat':self.region.ymax,
            'minlat':self.region.ymin
        }
        _req = fetches.Fetch(self._ngs_search_url).fetch_req(params=_data)
        if _req is not None:
            self.add_entry_to_results(
                _req.url,
                'ngs_results_{}.json'.format(self.region.format('fn')),
                'ngs'
            )
            
        return(self)

### End
