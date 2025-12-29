### xyzfile.py 
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
##
## xyzfile.py is part of CUDEM
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
## ASCII XYZ File Parser
##
### Code:

import os
import sys
import math
import numpy as np

from cudem import utils
from cudem import xyzfun
from cudem.datalists.dlim import ElevationDataset

class XYZFile(ElevationDataset):
    """Representing an ASCII xyz dataset stream.
    Parse data from an xyz file/stdin.
    """
            
    def __init__(self, delim=None, xpos=0, ypos=1, zpos=2,
                 wpos=None, upos=None, skip=0, x_scale=1, y_scale=1,
                 z_scale=1, x_offset=0, y_offset=0, use_numpy=True,
                 iter_rows=1000000, **kwargs):
        
        super().__init__(**kwargs)
        
        self.delim = delim
        self.xpos = utils.int_or(xpos, 0)
        self.ypos = utils.int_or(ypos, 1)
        self.zpos = utils.int_or(zpos, 2)
        self.wpos = utils.int_or(wpos)
        self.upos = utils.int_or(upos)
        self.skip = utils.int_or(skip, 0)
        
        self.x_scale = utils.float_or(x_scale, 1)
        self.y_scale = utils.float_or(y_scale, 1)
        self.z_scale = utils.float_or(z_scale, 1)
        
        self.x_offset = x_offset
        self.y_offset = utils.int_or(y_offset, 0)
        
        self.rem = False
        if self.x_offset == 'REM':
            self.x_offset = 0
            self.rem = True
        else:
            self.x_offset = utils.float_or(x_offset, 0)

        self.use_numpy = use_numpy
        self.iter_rows = iter_rows
        
        ## Check if scaling/offsetting is needed
        self.scoff = (self.x_scale != 1 or self.y_scale != 1 or self.z_scale != 1 or
                      self.x_offset != 0 or self.y_offset != 0)

        
    def generate_inf(self, make_grid=True, make_block_mean=False, block_inc=None):
        """Generate metadata (INF) for the XYZ file.
        Ensures delimiter is guessed before scanning.
        """
        if self.delim is None and self.fn is not sys.stdin:
            self.guess_delim()
            
        return super().generate_inf(make_grid, make_block_mean, block_inc)

    
    def guess_delim(self):
        """Guess the xyz delimiter by inspecting the first few lines."""
        
        if self.fn is None: return

        try:
            if os.path.exists(str(self.fn)):
                with open(self.fn, 'r') as f:
                    ## Check first 10 lines
                    for _ in range(10):
                        line = f.readline()
                        if not line: break
                        if line.strip().startswith('#'): continue
                        
                        for delim in xyzfun._known_delims:
                            if len(line.split(delim)) > 1:
                                self.delim = delim
                                return
            elif self.fn == sys.stdin:
                #$ Cannot peek stdin easily without consuming
                pass 
        except Exception:
            pass

        
    def line_delim(self, xyz_line):
        """Split a line using known delimiters."""
        
        for delim in xyzfun._known_delims:
            parts = xyz_line.split(delim)
            if len(parts) > 1:
                return parts
        return None

    
    def yield_points(self):
        """Yield points from the XYZ file.
        Uses numpy.loadtxt for speed if possible, falls back to manual parsing.
        """
        
        ## Register delimiter if provided
        if self.delim is not None and self.delim not in xyzfun._known_delims:
            xyzfun._known_delims.insert(0, self.delim)

        ## Setup Fields for Numpy
        cols = [self.xpos, self.ypos, self.zpos]
        names = ['x', 'y', 'z']
        
        if self.wpos is not None:
            cols.append(self.wpos)
            names.append('w')
            
        if self.upos is not None:
            cols.append(self.upos)
            names.append('u')
            
        ## Cannot use numpy on streams easily
        if self.fn is sys.stdin:
            self.use_numpy = False
            
        ## --- Numpy Approach ---
        if self.use_numpy:
            try:
                if self.delim is None: self.guess_delim()

                current_skip = self.skip            
                with open(self.fn, 'r') as src_data:
                    while True:
                        try:
                            ## Load chunk
                            points = np.loadtxt(
                                src_data,
                                delimiter=self.delim,
                                comments='#',
                                ndmin=1,
                                skiprows=current_skip,
                                usecols=cols,
                                dtype={'names': names, 'formats': [float] * len(names)},
                                max_rows=self.iter_rows
                            )
                            current_skip = 0 # Only skip header once

                            if points.size == 0: break

                            ## Ensure structured array view
                            if points.dtype.names is None:
                                ## Handle case where loadtxt returns unstructured array (single col/row edge cases)
                                pass 
                            
                            ## Apply Transforms
                            if self.scoff:
                                points['x'] = (points['x'] + self.x_offset) * self.x_scale
                                points['y'] = (points['y'] + self.y_offset) * self.y_scale
                                points['z'] *= self.z_scale

                            if self.rem:
                                points['x'] = np.fmod(points['x'] + 180, 360) - 180 

                            yield points.view(np.recarray)

                            if self.iter_rows is None or len(points) < self.iter_rows:
                                break
                                
                        except StopIteration:
                            break
            except Exception as e:
                if self.verbose:
                    utils.echo_warning_msg(f'Numpy load failed for {self.fn}: {e}, falling back to manual parse.')
                self.use_numpy = False

        ## --- Manual Approach ---
        if not self.use_numpy:
            ## Handle source
            src_data = sys.stdin if self.fn == sys.stdin else open(self.fn, 'r')
            
            try:
                chunk_x, chunk_y, chunk_z = [], [], []
                chunk_w, chunk_u = [], []
                count = 0
                
                for i, line in enumerate(src_data):
                    if i < self.skip: continue
                    if line.strip().startswith('#'): continue
                    
                    ## Split
                    parts = None
                    if self.delim:
                        parts = line.split(self.delim)
                    else:
                        parts = self.line_delim(line)
                        
                    if not parts or len(parts) < 3: continue
                    
                    try:
                        x = float(parts[self.xpos])
                        y = float(parts[self.ypos])
                        z = float(parts[self.zpos])
                        
                        w = float(parts[self.wpos]) if self.wpos is not None else 1.0
                        u = float(parts[self.upos]) if self.upos is not None else 0.0
                        
                        chunk_x.append(x)
                        chunk_y.append(y)
                        chunk_z.append(z)
                        chunk_w.append(w)
                        chunk_u.append(u)
                        count += 1
                        
                        ## Yield chunk if size reached
                        if count >= self.iter_rows:
                            points = self._build_recarray(chunk_x, chunk_y, chunk_z, chunk_w, chunk_u)
                            yield points
                            
                            chunk_x, chunk_y, chunk_z = [], [], []
                            chunk_w, chunk_u = [], []
                            count = 0
                            
                    except (ValueError, IndexError):
                        continue

                ## Yield remaining
                if count > 0:
                    yield self._build_recarray(chunk_x, chunk_y, chunk_z, chunk_w, chunk_u)

            finally:
                if src_data is not sys.stdin:
                    src_data.close()

                    
    def _build_recarray(self, x, y, z, w, u):
        """Helper to build recarray from lists and apply transforms."""
        
        x = np.array(x, dtype=float)
        y = np.array(y, dtype=float)
        z = np.array(z, dtype=float)
        w = np.array(w, dtype=float)
        u = np.array(u, dtype=float)
        
        if self.scoff:
            x = (x + self.x_offset) * self.x_scale
            y = (y + self.y_offset) * self.y_scale
            z *= self.z_scale

        if self.rem:
            x = np.fmod(x + 180, 360) - 180
            
        return np.rec.fromarrays([x, y, z, w, u], names='x, y, z, w, u')

    
class YXZFile(XYZFile):
    """yxz file shortcut (167 mbdatalist)"""
    
    def __init__(self, **kwargs):
        super().__init__(xpos=1, ypos=0, zpos=2, **kwargs)

### End
