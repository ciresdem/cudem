### xyzfun.py - CUDEM utilities and functions
##
## Copyright (c) 2010 - 2021 CIRES Coastal DEM Team
##
## xyzfun.py is part of CUDEM
##
## Permission is hereby granted, free of charge, to any person obtaining a copy 
## of this software and associated documentation files (the "Software"), to deal 
## in the Software without restriction, including without limitation the rights 
## to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies 
## of the Software, and to permit persons to whom the Software is furnished to do so, 
## subject to the following conditions:
##
## The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
##
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
## INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR 
## PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE 
## FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
## ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
##
### Commentary:
### Code:

import sys

from osgeo import gdal
from osgeo import ogr
from osgeo import osr
        
from cudem import utils

class XYZPoint:
    """represnting an xyz data point"""
    
    def __init__(
            self,
            x=None,
            y=None,
            z=None,
            w=1,
            src_srs='epsg:4326',
            z_units='m',
            z_datum='msl'
    ):
        """
        Args:
          x (float): longitude/x
          y (float): latitude/y
          z (float): elevation/z
          w (float): weight/w
          src_srs (str): the data projection
          z_units (str): units of the z value
          z_datum (str): vertical datum of the z_value
        """
        
        self.x = utils.float_or(x)
        self.y = utils.float_or(y)
        self.z = utils.float_or(z)
        self.w = utils.float_or(w)
        self.src_srs = utils.str_or(src_srs, 'epsg:4326')
        self.z_units = z_units
        self.z_datum = z_datum

    def __repr__(self):
        return('<XYZPoint x: {} y: {} z: {}>'.format(self.x, self.y, self.z))
        
    def __str__(self):
        return('<XYZPoint x: {} y: {} z: {}>'.format(self.x, self.y, self.z))
        
    def copy(self):
        return(
            XYZPoint(
                x=self.x,
                y=self.y,
                z=self.z,
                w=self.w
            )
        )

    def reset(self):
        self.x = self.y = self.z = None
        self.w = 1
        self.src_srs = 'epsg:4326'
        self.z_units = 'm'
        self.z_datum = 'msl'
    
    def valid_p(self):
        
        if self.x is None:
            return(False)
        
        if self.y is None:
            return(False)
        
        if self.z is None:
            return(False)
        
        return(True)
    
    def from_list(
            self,
            xyz_list,
            x_pos=0,
            y_pos=1,
            z_pos=2,
            w_pos=3
    ):
        """load xyz data from a list

        Args:
          xyz_list (list): a list of xyz data [x,y,z,...]

        Returns:
          xyz: self
        """

        if len(xyz_list) > x_pos:
            self.x = utils.float_or(xyz_list[x_pos])
            
        if len(xyz_list) > y_pos:
            self.y = utils.float_or(xyz_list[y_pos])
            
        if len(xyz_list) > z_pos:
            self.z = utils.float_or(xyz_list[z_pos])
            
        if len(xyz_list) > w_pos:
            self.w = utils.float_or(xyz_list[w_pos])
            
        return(self)
    
    def from_string(
            self,
            xyz_str,
            x_pos=0,
            y_pos=1,
            z_pos=2,
            w_pos=3,
            delim=" "
    ):
        """load xyz data from a string

        Args:
          xyz_str (str): a string representing delimited xyz data.

        Returns:
          xyz: self
        """
        
        this_line = xyz_str.strip()
        return(
            self.from_list(
                this_line.split(delim), x_pos, y_pos, z_pos, w_pos
            )
        )
        
    def export_as_list(self, include_z=False, include_w=False):
        """export xyz as a list

        Args:
          include_z: include the z-region in the output
          include_w: include the w-region in the output

        Returns:
          list: the region values in a list
        """
        
        xyz_list = [self.x, self.y]
        if include_z:
            xyz_list.append(self.z)
            
        if include_w:
            xyz_list.append(self.w)
            
        return(xyz_list)

    def export_as_string(self, delim, include_z=False, include_w=False):
        """export xyz data as string

        Args:
          delim (str): the delimiter
        """

        l = self.export_as_list(
            include_z=include_z, include_w=include_w
        )
        
        return('{}\n'.format(delim.join([str(x) for x in l])))

    def export_as_wkt(self, include_z=False):
        if include_z:
            return('POINT ({} {} {})'.format(self.x, self.y, self.z))
        else:
            return('POINT ({} {})'.format(self.x, self.y))
    
    def dump(
            self,
            delim=' ',
            include_z=True,
            include_w=False,
            encode=False,
            dst_port=sys.stdout
    ):
        """dump xyz as a string to dst_port

        Args:
          include_z: include the z-region in the output
          include_w: include the w-region in the output
          dst_port (port): an open destination port
          encode (bool): encode the output
            sys.stdout, etc prefer encoded data, while
            files usually like unencoded...
        """
    
        l = self.export_as_string(
            delim, include_z=include_z, include_w=include_w
        )
        
        if encode:
            l = l.encode('utf-8')
            
        dst_port.write(l)

    def transform(self, dst_trans):
        """transform the x/y using the dst_trans

        Args:
          dst_trans: an srs transformation object

        Returns:
          xyz: self
        """
        
        point = ogr.CreateGeometryFromWkt(
            'POINT ({} {})'.format(self.x, self.y)
        )
        
        point.Transform(dst_trans)
        self.x = point.GetX()
        self.y = point.GetY()

        return(self)
        
    def warp(self, dst_srs='epsg:4326'):
        """transform the x/y using dst_srs"""
        
        dst_srs = utils.str_or(dst_srs, 'epsg:4326')
        if dst_srs is None or self.src_srs is None:
            return(self)

        src_srs = osr.SpatialReference()
        src_srs.SetFromUserInput(self.src_srs)

        dst_srs_ = osr.SpatialReference()
        dst_srs_.SetFromUserInput(dst_srs)
        try:
            src_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
            dst_srs_.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        except: pass
        
        dst_trans = osr.CoordinateTransformation(src_srs, dst_srs_)
        
        if dst_trans is None:
            return(self)
        
        self.transform(dst_trans)
        
        return(self)

def xyz_line(xyz_line, dst_port=sys.stdout, encode=False):
    """write "xyz" `line` to `dst_port`
    `line` should be a list of xyz values [x, y, z, ...].
    
    Args:
      xyz_line (str): a string representing delimited data.
      dst_port (port): an open destination port
      encode (bool): encode the output
        sys.stdout, etc prefer encoded data, while
        files usually like unencoded...
    """
    
    delim = ' '
    
    l = '{}\n'.format(delim.join(['{:.7f}'.format(x) for x in xyz_line]))
    if encode: l = l.encode('utf-8')
    dst_port.write(l)

def xyz_chunks():
    pass

### End
