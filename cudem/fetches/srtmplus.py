### srtmplus.py
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

## SRTM Plus
class SRTMPlus(fetches.FetchModule):
    """SRTM15+: GLOBAL BATHYMETRY AND TOPOGRAPHY AT 15 ARCSECONDS.

    https://topex.ucsd.edu/WWW_html/srtm15_plus.html
    http://topex.ucsd.edu/sandwell/publications/180_Tozer_SRTM15+.pdf
    https://topex.ucsd.edu/pub/srtm15_plus/
    https://topex.ucsd.edu/pub/srtm15_plus/SRTM15_V2.3.nc
    
    < srtm_plus >
    """
    
    def __init__(self, **kwargs):
        super().__init__(name='srtm_plus', **kwargs)

        ## The srtm_plus URL
        self._srtm_url = 'https://topex.ucsd.edu/cgi-bin/get_srtm15.cgi'

        ## for dlim, data_format of 168 is xyz-file, skip the first line
        self.data_format = '168:skip=1'
        self.src_srs = 'epsg:4326+3855'

        self.title = 'SRTM+'
        self.source = 'Scripps Institude of Oceanography'
        self.date = '2014'
        self.data_type = 'Topographic/Bathymetric Raster'
        self.resolution = '60 arc-seconds'
        self.hdatum = 'WGS84'
        self.vdatum = 'MSL'
        self.url = self._srtm_url

        
    def run(self):
        '''Run the SRTM fetching module.'''
        
        if self.region is None:
            return([])

        self.data = {
            'north':self.region.ymax,
            'west':self.region.xmin,
            'south':self.region.ymin,
            'east':self.region.xmax
        }
        _req = fetches.Fetch(self._srtm_url, verify=False).fetch_req(params=self.data)
        if _req is not None:
            outf = 'srtm_{}.xyz'.format(self.region.format('fn'))
            self.add_entry_to_results(
                _req.url,
                outf,
                'srtm'
            )
            
### End
