### povray.py 
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

from cudem.perspecto import perspecto
from cudem import utils

class POVRay(perspecto.Perspecto):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.cpt_no_slash()
        
        self.rgb_image = '{}_rgb.png'.format(utils.fn_basename2(self.src_dem))
        self.dem_image = '{}_16bit.png'.format(utils.fn_basename2(self.src_dem))

        #if not os.path.exists(self.rgb_image) or not os.path.exists(self.dem_image):
        self.export_as_png()            
        self.output_pov = '{}.pov'.format(utils.fn_basename2(self.src_dem))

        
    def run_povray(self, src_pov_template, pov_width=800, pov_height=800):
        utils.run_cmd(
            'povray {} +W{} +H{} -D'.format(src_pov_template, pov_width, pov_height),
            verbose=True
        )

### End
