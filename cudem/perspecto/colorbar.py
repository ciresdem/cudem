### colorbar.py 
##
## Copyright (c) 2023 - 2026 Regents of the University of Colorado
##
## colorbar.py is part of CUDEM
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
## Generate a standalone colorbar image using PyGMT.
##
### Code:

from cudem import utils
from cudem.perspecto import perspecto

try:
    import pygmt
    HAS_PYGMT = True
except ImportError:
    HAS_PYGMT = False


class Colorbar(perspecto.Perspecto):
    """
    Generate a standalone colorbar image using PyGMT.
    
    Configuration Example:
    < colorbar:colorbar_text='Elevation':width=10:height=2 >
    """
    
    def __init__(self, colorbar_text='Elevation', width=10, height=2, **kwargs):
        # Initialize parent with 'colorbar' module name
        super().__init__(mod='colorbar', want_gdal_cpt=False, **kwargs)
        self.colorbar_text = colorbar_text
        self.width = utils.int_or(width)
        self.height = utils.int_or(height)

        
    def run(self):
        if not HAS_PYGMT:
            utils.echo_error_msg("PyGMT is not installed, cannot generate colorbar.")
            return None

        basename = utils.fn_basename2(self.src_dem)
        outfile = f'{basename}_cbar.png'

        fig = pygmt.Figure()
        
        ## Draw the colorbar
        ## Position: Top Center (jTC), width/height specified in cm (+w), horizontal (+h)
        fig.colorbar(
            cmap=self.cpt,  # Ensure CPT is passed from parent/init
            region=[0, 10, 0, 3],
            projection="X10c/3c",
            position=f'jTC+w{self.width}c/{self.height}c+h',
            frame=[f'x+l{self.colorbar_text}', 'y+1m']
        )
        
        fig.savefig(outfile)
        return outfile

### END
