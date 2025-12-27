### hydrolakes.py
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
##
## hydrolakes.py is part of CUDEM
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
## Fetch HydroLAKES and GLOBathy data.
##
### Code:

import os
from typing import List, Dict, Optional
from cudem.fetches import fetches

## ==============================================
## Constants
## ==============================================
HYDROLAKES_POLY_ZIP_URL = 'https://data.hydrosheds.org/file/hydrolakes/HydroLAKES_polys_v10_shp.zip'
HYDROLAKES_GDB_ZIP_URL = 'https://data.hydrosheds.org/file/hydrolakes/HydroLAKES_polys_v10_gdb.zip'
GLOBATHY_URL = 'https://springernature.figshare.com/ndownloader/files/28919991'

## ==============================================
## HydroLakes Module
## ==============================================
class HydroLakes(fetches.FetchModule):
    """HydroLakes vector and derived elevations.
    
    HydroLAKES provides the shoreline polygons of all global lakes 
    with a surface area of at least 10 ha.
    
    GLOBathy provides global lakes bathymetry.

    < hydrolakes:want_globathy=False >
    """
    
    def __init__(self, where: str = '1=1', want_globathy: bool = False, **kwargs):
        super().__init__(name='hydrolakes', **kwargs)
        self.want_globathy = want_globathy
        self.where = [where] if where else []

        
    def run(self):
        """Run the hydrolakes fetching module."""
        
        ## Add HydroLAKES polygon zip
        self.add_entry_to_results(
            HYDROLAKES_POLY_ZIP_URL,
            os.path.basename(HYDROLAKES_POLY_ZIP_URL),
            'hydrolakes'
        )

        ## Add GLOBathy if requested
        if self.want_globathy:
            self.add_entry_to_results(
                GLOBATHY_URL,
                'globathy_parameters.zip',
                'globathy'
            )
                        
        return self

### End
