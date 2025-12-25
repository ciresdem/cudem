### gmtimage.py 
##
## Copyright (c) 2023 - 2025 Regents of the University of Colorado
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
## Base class for PyGMT-based image generation.
##
### Code:

from cudem import utils
from cudem.perspecto import perspecto

try:
    import pygmt
    HAS_PYGMT = True
except ImportError:
    HAS_PYGMT = False


class GMTImage(perspecto.Perspecto):
    """
    Base class for generating images using PyGMT.
    Ensures input is in NetCDF format and loads the grid.
    """
    
    def __init__(self, **kwargs):
        super().__init__(want_gdal_cpt=False, **kwargs)

        if not HAS_PYGMT:
            utils.echo_error_msg("PyGMT is not installed.")
            return

        self.src_dem_netcdf = None
        # Ensure input is NetCDF for GMT
        if self.dem_infos['fmt'] != 'NetCDF':
            basename = utils.fn_basename2(self.src_dem)
            nc_filename = f"{basename}.nc"
            
            utils.run_cmd(f'gmt grdconvert {self.src_dem} {nc_filename}')
            self.src_dem_netcdf = nc_filename
        
        # Load the grid data
        if self.src_dem_netcdf is not None:
            self.grid = pygmt.load_dataarray(self.src_dem_netcdf)
        
            # Initialize CPT
            self.makecpt(cmap=self.cpt, output=None)
        else:
            return

### End
