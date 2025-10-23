### colorbar.py 
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
###############################################################################
### Commentary:
##
##
### Code:

from cudem.perspecto import gmtimage
from cudem import utils
import pygmt

class colorbar(gmtimage.GMTImage):
    """Generate a colorbar

    uses pyGMT

< colorbar:colorbar_text='Elevation':widt=10:height=2 >
    """
    
    def __init__(self, colorbar_text='Elevation', width=10,
                 height=2, **kwargs):
        super().__init__(**kwargs)
        self.colorbar_text=colorbar_text
        self.width = utils.int_or(width)
        self.height = utils.int_or(height)

        
    def run(self):
        fig = pygmt.Figure()
        fig.colorbar(
            region=[0,10,0,3],
            projection="X10c/3c",
            position='jTC+w{}c/{}c+h'.format(self.width, self.height),
            frame=['x+l{}'.format(self.colorbar_text), 'y+1m']
        )
        fig.savefig('{}_cbar.png'.format(utils.fn_basename2(self.src_dem)))
        return('{}_cbar.png'.format(utils.fn_basename2(self.src_dem)))


### END
