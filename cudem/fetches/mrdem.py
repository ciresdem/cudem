### mrdem.py
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

class MRDEM(fetches.FetchModule):
    def __init__(self, **kwargs):
        super().__init__(name='mrdem', **kwargs)
        # self.mrdem_dtm_url = ('https://datacube-prod-data-public.s3.ca-central-1.'
        #                       'amazonaws.com/store/elevation/mrdem/mrdem-30/'
        #                       'mrdem-30-dtm.vrt')
        # self.mrdem_dsm_url = ('https://datacube-prod-data-public.s3.ca-central-1.'
        #                       'amazonaws.com/store/elevation/mrdem/mrdem-30/'
        #                       'mrdem-30-dsm.vrt')

        self.mrdem_dtm_url = 'https://canelevation-dem.s3.ca-central-1.amazonaws.com/mrdem-30/mrdem-30-dtm.vrt'
        self.mrdem_dsm_url = 'https://canelevation-dem.s3.ca-central-1.amazonaws.com/mrdem-30/mrdem-30-dsm.vrt'

        
    def run(self):
        self.add_entry_to_results(
            self.mrdem_dtm_url,
            self.mrdem_dtm_url.split('/')[-1],
            'vrt'
        )

### End
