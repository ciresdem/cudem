### povray.py 
##
## Copyright (c) 2023 - 2025 Regents of the University of Colorado
##
## povray.py is part of CUDEM
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
## Helper class for POV-Ray rendering.
##
### Code:

import os
from cudem import utils
from . import perspecto

class PovRay(perspecto.Perspecto):
    """Helper class for generating POV-Ray scenes from DEMs."""
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.cpt_no_slash()
        
        basename = utils.fn_basename2(self.src_dem)
        self.rgb_image = f'{basename}_rgb.png'
        self.dem_image = f'{basename}_16bit.png'
        self.output_pov = f'{basename}.pov'

        # Generate intermediate PNGs if they don't exist
        # (Uncomment the check if you want to skip regeneration)
        # if not os.path.exists(self.rgb_image) or not os.path.exists(self.dem_image):
        self.export_as_png()            

        
    def run_povray(self, src_pov_template, pov_width=800, pov_height=800):
        """Executes POV-Ray to render the scene."""
        
        if not os.path.exists(src_pov_template):
            utils.echo_error_msg(f"POV-Ray template not found: {src_pov_template}")
            return

        cmd = f'povray {src_pov_template} +W{pov_width} +H{pov_height} -D'
        utils.run_cmd(cmd, verbose=True)

### End
