### clip.py - GRId filTerS
##
## Copyright (c) 2024 - 2025 Regents of the University of Colorado
##
## clip.py is part of cudem
##
## Permission is hereby granted, free of charge, to any person obtaining a copy 
## of this software and associated documentation files (the "Software"), to deal 
## in the Software without restriction, including without limitation the rights 
## to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies 
## of the Software, and to permit persons to whom the Software is furnished to do so, 
## subject to the following conditions:
##
## The above copyright notice and this permission notice shall be included in all copies or
## substantial portions of the Software.
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
## Clips a DEM to an OGR-compatible Vector source (Shapefile, etc.).
## Wraps gdalfun.gdal_clip.
##
### Code:

import os
from cudem import utils
from cudem import gdalfun
from cudem import factory
from cudem.grits import grits

class Clip(grits.Grits):
    """Clip the input DEM to a Vector Polygon.
    
    Parameters:
    -----------
    clip_str : str
        The clip configuration string in the format:
        `src_ply:arg1=val1:arg2=val2...`
        
        Where `src_ply` is the path to the vector file.
        Arguments are passed to gdal_clip (e.g., layer, where, invert).
    """
    
    def __init__(self, clip_str=None, **kwargs):
        super().__init__(**kwargs)
        self.clip_str = clip_str

    def run(self):
        """Execute the clip filter."""
        
        if self.clip_str is None:
            return self.src_dem, 0

        ## Parse Clip String
        ## Format: filename.shp:layer=layername:where="id=1"
        clip_args = {}
        parts = self.clip_str.split(':')
        src_ply = parts[0]
        
        ## Parse optional args (layer, where, etc.)
        if len(parts) > 1:
            clip_args = factory.args2dict(parts[1:], clip_args)

        ## Check existence
        if not os.path.exists(src_ply):
            utils.echo_error_msg(f"Clip vector not found: {src_ply}")
            return self.src_dem, -1

        ## Use the destination defined by Grits base
        dst_dem = self.dst_dem

        if self.verbose:
            utils.echo_msg(f"Clipping {self.src_dem} to {src_ply}...")

        ## Perform Clip
        ## gdalfun.gdal_clip handles the gdal.Warp logic with cutlineDSName
        res = gdalfun.gdal_clip(
            self.src_dem, 
            dst_dem, 
            src_ply=src_ply, 
            **clip_args
        )
        
        ## Result tuple is (filename, status)
        if res[1] != 0:
            utils.echo_error_msg(f"Failed to clip data to {src_ply}")
            return self.src_dem, -1
            
        return dst_dem, 0

### End
