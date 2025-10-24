### hydrolakes.py
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

## hydrolakes
class HydroLakes(fetches.FetchModule):
    """HydroLakes vector and derived elevations
    
    HydroLAKES aims to provide the shoreline polygons of all global lakes 
    with a surface area of at least 10 ha. HydroLAKES has been developed 
    using a suite of auxiliary data sources of lake polygons and gridded 
    lake surface areas. All lakes are co-registered to the global river 
    network of the HydroSHEDS database via their lake pour points. The 
    global coverage of HydroLAKES encompasses 1.4 million individual lakes 
    or reservoirs representing a total surface area of 2.67 million km², 
    a total shoreline length of 7.2 million km, and a total storage volume 
    of 181,900 km³.
    
    hydrolakes:
    https://wp.geog.mcgill.ca/hydrolab/data/
    https://www.hydrosheds.org/products/hydrolakes
    https://data.hydrosheds.org/file/hydrolakes/HydroLAKES_polys_v10_shp.zip

    globathy:
    https://springernature.figshare.com/collections/GLOBathy_the_Global_Lakes_Bathymetry_Dataset/5243309

    add globathy output with the 'want_globathy' flag set to True

    < hydrolakes:want_globathy=False >
    """
    
    def __init__(self, where='1=1', want_globathy=False, **kwargs):
        super().__init__(name='hydrolakes', **kwargs)
        self.want_globathy = want_globathy
        self.where = [where] if len(where) > 0 else []
        
        ## The various hydrolakes/globathy URLs
        self._hydrolakes_prods = 'https://www.hydrosheds.org/products/hydrolakes'
        self._hydrolakes_poly_zip = ('https://data.hydrosheds.org/file/hydrolakes/'
                                     'HydroLAKES_polys_v10_shp.zip')
        self._hydrolakes_gdb_zip = ('https://data.hydrosheds.org/file/hydrolakes/'
                                    'HydroLAKES_polys_v10_gdb.zip')
        self._globathy_url = 'https://springernature.figshare.com/ndownloader/files/28919991'

        
    def run(self):
        """Run the hydrolakes URLs"""
        
        self.add_entry_to_results(
            self._hydrolakes_poly_zip,
            self._hydrolakes_poly_zip.split('/')[-1],
            'hydrolakes'
        )

        if self.want_globathy:
            self.add_entry_to_results(
                self._globathy_url,
                'globathy_parameters.zip',
                'globathy'
            )
                        
        return(self)

### End
