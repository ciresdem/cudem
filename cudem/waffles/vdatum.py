### vdatum.py
##
## Copyright (c) 2010 - 2026 Regents of the University of Colorado
##
## vdatum.py is part of CUDEM
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
### Code:

from cudem import utils
from cudem import vdatums
from cudem.waffles.waffles import Waffle

class WafflesVDatum(Waffle):
    """VDATUM transformation grid.
    Generate a Vertical DATUM transformation grid.

    Parameters:
    -----------
    mode=[waffle-module] - the waffles module to use to interpolate/extrapolate 
    vdatum_in=[vdatum] - input vertical datum
    vdatum_out=[vdatum] - output vertical datum    
    """

    def __init__(self, mode='IDW', vdatum_in=None, vdatum_out=None, **kwargs):
        self.valid_modes = [
            'gmt-surface',
            'IDW',
            'linear',
            'cubic',
            'nearest',
            'gmt-triangulate',
            'mbgrid'
        ]
        self.mode = mode
        self.mode_args = {}
        if self.mode not in self.valid_modes:
            self.mode = 'IDW'
            
        tmp_waffles = Waffle()
        tmp_waffles_mode = WaffleFactory(mod=self.mode)._acquire_module()
        for kpam, kval in kwargs.items():
            if kpam not in tmp_waffles.__dict__:
                if kpam in tmp_waffles_mode.__dict__:
                    self.mode_args[kpam] = kval

        for kpam, kval in self.mode_args.items():
            del kwargs[kpam]

        super().__init__(**kwargs)
        self.vdatum_in = vdatum_in
        self.vdatum_out = vdatum_out        

        
    def run(self):
        status = vdatums.VerticalTransform(
            self.mode,
            self.p_region,
            self.xinc,
            self.yinc,
            self.vdatum_in,
            self.vdatum_out,
            node=self.node,
            cache_dir=waffles_cache
        ).run(outfile='{}.tif'.format(self.name))

        return self


### End
