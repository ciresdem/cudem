### cptcity.py
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
##
## cptcity.py is part of CUDEM
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
## Fetch CPT color palettes from CPT City.
##
### Code:

import zipfile
from io import BytesIO
from typing import Optional
import requests
from cudem import utils
from cudem.fetches import fetches

## ==============================================
## Constants
## ==============================================
CPT_PUB_URL = 'http://seaviewsensing.com/pub/'
CPT_PKG_BASE_URL = 'http://seaviewsensing.com/pub/cpt-city/pkg/'
PACKAGE_XML_URL = f"{CPT_PKG_BASE_URL}package.xml"

## ==============================================
## CPT City Module
## ==============================================
class CPTCity(fetches.FetchModule):
    """CPT City

    Fetch various CPT files for DEM hillshades, etc.
    """
    
    def __init__(self, q: Optional[str] = None, **kwargs):
        super().__init__(name='cpt_city', **kwargs)
        self.q = q

        
    def run(self):
        """Run the cpt-city fetches module."""
        
        ## Fetch Package XML
        try:
            cpt_xml = fetches.iso_xml(PACKAGE_XML_URL)
            if cpt_xml.xml_doc is None:
                utils.echo_error_msg("Failed to parse CPT City package.xml")
                return self
                
            cpt_node = cpt_xml.xml_doc.find('cpt')
            if cpt_node is None or not cpt_node.text:
                 utils.echo_error_msg("Could not find 'cpt' tag in package.xml")
                 return self
                 
            cpt_zip_filename = cpt_node.text
            
            ## Fetch the Main Zip
            zip_url = f"{CPT_PKG_BASE_URL}{cpt_zip_filename}"
            req = requests.get(zip_url)
            req.raise_for_status()
            
            with zipfile.ZipFile(BytesIO(req.content)) as zip_ref:
                zip_cpts = zip_ref.namelist()

            ## Filter results
            if self.q:
                ## Simple substring match
                filtered_files = [x for x in zip_cpts if self.q in x]
            else:
                filtered_files = zip_cpts

            ## Generate download links (pointing to individual files on the server)
            for f in filtered_files:
                f_url = f"{CPT_PUB_URL}{f}"
                f_fn = f.split('/')[-1]
                
                ## Filter for actual files, not directories
                if not f.endswith('/'):
                    self.add_entry_to_results(f_url, f_fn, 'cpt')

        except Exception as e:
            utils.echo_error_msg(f"Error running CPT City fetch: {e}")

        return self

### End
