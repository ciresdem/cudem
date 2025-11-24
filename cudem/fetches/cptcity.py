### cptcity.py
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

import zipfile
from io import BytesIO
import requests
from cudem.fetches import fetches

class CPTCity(fetches.FetchModule):
    """CPT City

    Fetch various CPT files for DEM hillshades, etc.
    """
    
    def __init__(self, q = None, **kwargs):
        super().__init__(name='cpt_city', **kwargs)
        self.q = q

        ## The various cpt-city URLs
        #self.cpt_pub_url = 'http://soliton.vm.bytemark.co.uk/pub/' # dead url
        self.cpt_pub_url = 'http://seaviewsensing.com/pub/'
        self.cpt_pkg_url = self.cpt_pub_url + 'cpt-city/pkg/'

        
    def run(self):
        """Run the cpt-city fetches module"""
        
        cpt_xml = fetches.iso_xml(self.cpt_pkg_url + "package.xml")
        print(cpt_xml)
        cpt_url_bn = cpt_xml.xml_doc.find('cpt').text
        cpt_zip = requests.get(self.cpt_pkg_url + cpt_url_bn)
        zip_ref = zipfile.ZipFile(BytesIO(cpt_zip.content))
        zip_cpts = zip_ref.namelist()

        if self.q is not None:
            mask = [self.q in x for x in zip_cpts]
            ff = [b for a, b in zip(mask, zip_cpts) if a]
        else:
            ff = zip_cpts
            
        [self.add_entry_to_results(
            self.cpt_pub_url + f, f.split('/')[-1], 'cpt'
        ) for f in ff]

### End
