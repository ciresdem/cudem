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
###############################################################################
### Commentary:
##
### Examples:
##
### TODO:
##
### Code:

import os
import sys
import numpy as np

from cudem import utils
from cudem import xyzfun
from cudem.datalists.dlim import ElevationDataset

class XYZFile(ElevationDataset):
    """representing an ASCII xyz dataset stream.

    Parse data from an xyz file/stdin

    generate_inf - generate an inf file for the xyz data
    yield_xyz - yield the xyz data as xyz
    yield_array - yield the xyz data as an array must set the 
                  x_inc/y_inc in the super class
    
    -----------
    Parameters:
    
    delim: the delimiter of the xyz data (str)
    xpos: the position (int) of the x value
    ypos: the position (int) of the y value
    zpos: the position (int) of the z value
    wpos: the position (int) of the w value (weight)
    upos: the position (int) of the u value (uncertainty)
    skip: number of lines to skip
    x_scale: scale the x value
    y_scale: scale the y value
    z_scale: scale the z value
    x_offset: offset the x value
    y_offset: offset the y value    
    """
            
    def __init__(self, delim=None, xpos=0, ypos=1, zpos=2,
                 wpos=None, upos=None, skip=0, x_scale=1, y_scale=1,
                 z_scale=1, x_offset=0, y_offset=0, use_numpy=True,
                 iter_rows=1000000, **kwargs):
        
        super().__init__(**kwargs)
        self.delim = delim # the file delimiter
        self.xpos = utils.int_or(xpos, 0) # the position of the x/lat value
        self.ypos = utils.int_or(ypos, 1) # the position of the y/lon value
        self.zpos = utils.int_or(zpos, 2) # the position of the z value
        self.wpos = utils.int_or(wpos) # the position of the weight value
        self.upos = utils.int_or(upos) # the position of the uncertainty value
        self.skip = utils.int_or(skip, 0) # number of header lines to skip
        self.x_scale = utils.float_or(x_scale, 1) # multiply x by x_scale
        self.y_scale = utils.float_or(y_scale, 1) # multiply y by y_scale
        self.z_scale = utils.float_or(z_scale, 1) # multiply z by z_scale
        self.x_offset = x_offset # offset x by x_offset
        self.y_offset = utils.int_or(y_offset, 0) # offset y by y_offset
        self.rem = False # x is 360 instead of 180
        self.use_numpy = use_numpy # use numpy.loadtxt to load the xyz points
        self.iter_rows = iter_rows # max rows to process at a time

        
    def yield_points(self):
        if self.delim is not None:
            xyzfun._known_delims.insert(0, self.delim)
                
        if self.x_offset == 'REM':
            self.x_offset = 0
            self.rem = True
            
        self.scoff = True if self.x_scale != 1 \
            or self.y_scale != 1 \
                or self.z_scale != 1 \
                    or self.x_offset != 0 \
                        or self.y_offset != 0 else False


        self.field_names  = [x for x in ['x' if self.xpos is not None else None,
                                         'y' if self.ypos is not None else None,
                                         'z' if self.zpos is not None else None,
                                         'w' if self.wpos is not None else None,
                                         'u' if self.upos is not None else None]
                             if x is not None]
        self.field_formats = [float for x in [self.xpos, self.ypos,
                                              self.zpos, self.wpos,
                                              self.upos] if x is not None]
        if self.fn is sys.stdin:
            self.use_numpy = False
            
        if self.use_numpy:
            try:
                if self.delim is None:
                    self.guess_delim()

                skip_ = self.skip            
                with open(self.fn, 'r') as src_data:
                    while True:
                        points = np.loadtxt(
                            self.fn,
                            delimiter=self.delim,
                            comments='#',
                            ndmin = 1,
                            skiprows=skip_,
                            usecols=[x for x in [self.xpos, self.ypos,
                                                 self.zpos, self.wpos,
                                                 self.upos] if x is not None],
                            dtype={'names': self.field_names,
                                   'formats': self.field_formats},
                            max_rows=self.iter_rows
                        )
                        skip_ = 0

                        if self.scoff:
                            points['x'] = (points['x'] + self.x_offset) * self.x_scale
                            points['y'] = (points['y'] + self.y_offset) * self.y_scale
                            points['z'] *= self.z_scale

                        if self.rem:
                            points['x'] = np.fmod(points['x'] + 180, 360) - 180 

                        points = points.view(np.recarray)                        
                        yield(points)

                        if self.iter_rows is None or len(points) < self.iter_rows:
                            break

            ## old processing function used as a fallback for when numpy.loadtxt fails
            except Exception as e:
                utils.echo_warning_msg(
                    f'could not load xyz data from {self.fn}, {e}, falling back'
                )
                self.use_numpy = False

        if not self.use_numpy:
            if self.fn is not None:
                if os.path.exists(str(self.fn)):
                    self.src_data = open(self.fn, "r")
                else:
                    self.src_data = self.fn
            else:
                self.src_data = sys.stdin

            points_x = []
            points_y = []
            points_z = []
            points_w = []
            points_u = []
            count = 0
            skip = self.skip
            for xyz_line in self.src_data:
                if count >= skip:
                    this_xyz = self.line_delim(xyz_line)
                    if this_xyz is None:
                        continue

                    if len(this_xyz) < 3:
                        continue
                    
                    x = utils.float_or(this_xyz[self.xpos])
                    y = utils.float_or(this_xyz[self.ypos])
                    z = utils.float_or(this_xyz[self.zpos])
                    w = utils.float_or(this_xyz[self.wpos]) \
                        if self.wpos is not None else 1
                    u = utils.float_or(this_xyz[self.upos]) \
                        if self.upos is not None else 0

                    if x is None or y is None or z is None:
                        continue
                    
                    if self.scoff:
                        x = (x+self.x_offset) * self.x_scale
                        y = (y+self.y_offset) * self.y_scale
                        z *= self.z_scale

                    if self.rem:
                        x = math.fmod(x+180,360)-180 

                    if self.data_region is not None \
                       and self.data_region.valid_p():
                        try:
                            this_xyz = xyzfun.XYZPoint(
                                x=this_xyz[self.xpos],
                                y=this_xyz[self.ypos],
                                z=this_xyz[self.zpos]
                            )
                        except Exception as e:
                            utils.echo_error_msg(f'{e} ; {this_xyz}')
                            this_xyz = xyzfun.XYZPoint()
                            
                    points_x.append(x)
                    points_y.append(y)
                    points_z.append(z)
                    points_w.append(w)
                    points_u.append(u)
                    count += 1
                else:
                    skip -= 1

            try:
                self.src_data.close()
            except:
                pass
            
            dataset = np.column_stack(
                (points_x, points_y, points_z, points_w, points_u)
            )
            points = np.rec.fromrecords(
                dataset, names='x, y, z, w, u'
            )

            yield(points)

            
    def guess_delim(self):
        """guess the xyz delimiter"""

        if self.fn is not None:
            if os.path.exists(str(self.fn)):
                self.src_data = open(self.fn, "r")
            else:
                self.src_data = self.fn
        else:
            self.src_data = sys.stdin

        for xyz_line in self.src_data:
            for delim in xyzfun._known_delims:
                try:
                    this_xyz = xyz_line.split(delim)
                    if len(this_xyz) > 1:
                        self.delim = delim
                        break
                except:
                    pass
            break

        try:
            if self.src_data is not sys.stdin:
                self.src_data.close()
        except:
            pass

        
    def line_delim(self, xyz_line):
        """guess a line delimiter and return the split line."""

        for delim in xyzfun._known_delims:
            try:
                this_xyz = xyz_line.split(delim)
                if len(this_xyz) > 1:
                    return(this_xyz)
            except:
                pass

class YXZFile(XYZFile):
    """yxz file shortcut (167 mbdatalist)"""
    
    def __init__(self, **kwargs):
        super().__init__(xpos=1, ypos=0, zpos=2, **kwargs)



### End
