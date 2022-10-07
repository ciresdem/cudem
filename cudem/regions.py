### regions.py
##
## Copyright (c) 2010 - 2022 CIRES Regents of the University of Colorado
##
## regions.py is part of CUDEM
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

import os
import sys

import numpy as np
from osgeo import ogr
from osgeo import osr

import cudem
from cudem import utils

## ==============================================
##
## regions
## the region class and associated functions.
##
## ==============================================
class Region:
    """Representing a geographic region.
    
    Attributes:
      xmin (float): the minimum x(long) value
      xmax (float): the maximum x(long) value
      ymin (float): the minimum x(lat) value
      ymax (float): the maximum x(lat) value
      zmin (float): the minimum z(elev) value
      zmax (float): the maximum z(elev) value
      wmin (float): the minimum w(weight) value
      wmax (float): the maximum w(weight) value
      wkt (str): the wkt representation of the region
      src_srs (str): the regions projection
    """

    full_region = lambda x: [x.xmin, x.xmax, x.ymin, x.ymax, x.zmin, x.zmax, x.wmin, x.wmax]
    xy_region = lambda x: [x.xmin, x.xmax, x.ymin, x.ymax]
    z_region = lambda x: [x.zmin, x.zmax]
    w_region = lambda x: [x.wmin, x.wmax]
    
    def __init__(self,
                 xmin=None,
                 xmax=None,
                 ymin=None,
                 ymax=None,
                 zmin=None,
                 zmax=None,
                 wmin=None,
                 wmax=None,
                 src_srs='epsg:4326',
                 wkt=None):        
        self.xmin = utils.float_or(xmin)
        self.xmax = utils.float_or(xmax)
        self.ymin = utils.float_or(ymin)
        self.ymax = utils.float_or(ymax)
        self.zmin = utils.float_or(zmin)
        self.zmax = utils.float_or(zmax)
        self.wmin = utils.float_or(wmin)
        self.wmax = utils.float_or(wmax)
        self.src_srs = utils.str_or(src_srs, 'epsg:4326')
        self.wkt = wkt

    def __repr__(self):
        return(str(self.format('fstr')))
        
    def __str__(self):
        return(str(self.format('fstr')))
        
    def valid_p(self, check_xy=False):
        """check the validity of a region
    
        Returns:
          bool: True if region  appears to be valid else False
        """
        
        if check_xy:
            if self.xmin is None: return(False)
            if self.xmax is None: return(False)
            if self.ymin is None: return(False)
            if self.ymax is None: return(False)

        if self.xmin is not None and self.xmax is not None:
            if self.xmin > self.xmax: return(False)
        if self.ymin is not None and self.ymax is not None:
            if self.ymin > self.ymax: return(False)
        if self.zmin is not None and self.zmax is not None:
            if self.zmin > self.zmax: return(False)
        if self.wmin is not None and self.wmax is not None:
            if self.wmin > self.wmax: return(False)

        if not any(self.full_region()): return(False)
        return(True)

    def copy(self):
        """return a copy of the region."""
        
        return(Region(xmin = self.xmin,
                      xmax = self.xmax,
                      ymin = self.ymin,
                      ymax = self.ymax,
                      zmin = self.zmin,
                      zmax = self.zmax,
                      wmin = self.wmin,
                      wmax = self.wmax,
                      src_srs = self.src_srs,
                      wkt = self.wkt))

    def _wgs_extremes(self, just_below=False):
        if self.xmin <= -180:
            self.xmin = -180 if not just_below else -179.99999
        if self.xmax >= 180:
            self.xmax = 180 if not just_below else 179.99999
        if self.ymin <= -90:
            self.ymin = -90 if not just_below else -89.99999
        if self.ymax >= 90:
            self.ymax = 90 if not just_below else 89.99999
            
    
    def from_list(self, region_list):
        """import a region from a region list 

        Args:
          region_list (list): [xmin, xmax, ymin, ymax[, zmin, zmax[, wmin, wmax]]]

        Returns:
          region-object: self
        """
        
        if len(region_list) >= 4:
            self.xmin = utils.float_or(region_list[0])
            self.xmax = utils.float_or(region_list[1])
            self.ymin = utils.float_or(region_list[2])
            self.ymax = utils.float_or(region_list[3])
            if len(region_list) >= 6:
                self.zmin = utils.float_or(region_list[4])
                self.zmax = utils.float_or(region_list[5])
                if len(region_list) >= 8:
                    self.wmin = utils.float_or(region_list[6])
                    self.wmax = utils.float_or(region_list[7])
            if self.wkt is None:
                if self.valid_p(check_xy = True):
                    self.wkt = self.export_as_wkt()
                else: self.wkt = None
        return(self)
        
    def from_string(self, region_str):
        """import a region from a region string 

        Args:
          region_str (str): <-R>xmin/xmax/ymin/ymax/zmin/zmax

        Returns:
          region-object: self
        """
        
        if region_str[:2] == '-R':
            region_str = region_str[2:]
            
        str_list = region_str.strip().split('/')
        if len(str_list) >= 4:
            r_list = [utils.float_or(x) for x in str_list]
            self.from_list(r_list)
        elif region_str.split()[0] == "POLYGON" or region_str.split()[0] == "MULTIPOLYGON":
            self.wkt = region_str
            self.from_list(ogr.CreateGeometryFromWkt(region_str).GetEnvelope())
        return(self)
    
    def from_geo_transform(self, geo_transform=None, x_count=None, y_count=None):
        """import a region from a region string 

        Args:
          geoT (list): a geo transform list
          x_count (int): length of x axis
          y_count (int): length of y axis

        Returns:
          region-object: self
        """
        
        if geo_transform is not None and x_count is not None and y_count is not None:            
            self.xmin = geo_transform[0]
            self.xmax = geo_transform[0] + geo_transform[1] * x_count
            self.ymin = geo_transform[3] + geo_transform[5] * y_count
            self.ymax = geo_transform[3]
        if self.wkt is None:
            self.wkt = self.export_as_wkt()
        return(self)

    def from_region(self, input_region):
        """import a region from another region"""

        return(input_region.copy())
    
    def format(self, t='gmt'):
        """format region to string, defined by `t`

        Args:
          t (str): output mode
            t = 'str': xmin/xmax/ymin/ymax
            t = 'fstr': xmin/xmax/ymin/ymax/zmin/zmax/wmin/wmax
            t = 'sstr': xmin xmax ymin ymax
            t = 'gmt': -Rxmin/xmax/ymin/ymax
            t = 'bbox': xmin,ymin,xmax,ymax
            t = 'osm_bbox': ymin,xmin,ymax,xmax
            t = 'te': xmin ymin xmax ymax
            t = 'ul_lr': xmin ymax xmax ymin
            t = 'fn': ymax_xmin

        Returns:
          str: the formatted region as str or None if region is invalid
        """

        if self.valid_p():
            if t == 'str': return('/'.join([str(self.xmin), str(self.xmax), str(self.ymin), str(self.ymax)]))
            elif t == 'sstr': return(' '.join([str(self.xmin), str(self.xmax), str(self.ymin), str(self.ymax)]))
            elif t == 'fstr': return(' '.join([str(self.xmin), str(self.xmax), str(self.ymin), str(self.ymax), str(self.zmin), str(self.zmax), str(self.wmin), str(self.wmax)]))
            elif t == 'gmt': return('-R' + '/'.join([str(self.xmin), str(self.xmax), str(self.ymin), str(self.ymax)]))
            elif t == 'bbox': return(','.join([str(self.xmin), str(self.ymin), str(self.xmax), str(self.ymax)]))
            elif t == 'osm_bbox': return(','.join([str(self.ymin), str(self.xmin), str(self.ymax), str(self.xmax)]))
            elif t == 'te': return(' '.join([str(self.xmin), str(self.ymin), str(self.xmax), str(self.ymax)]))
            elif t == 'ul_lr': return(' '.join([str(self.xmin), str(self.ymax), str(self.xmax), str(self.ymin)]))
            elif t == 'fn':
                ns = 's' if self.ymax < 0 else 'n'
                ew = 'e' if self.xmin > 0 else 'w'
                return('{}{:02d}x{:02d}_{}{:03d}x{:02d}'.format(ns, abs(int(self.ymax)), abs(int(self.ymax * 100)) % 100, 
                                                                ew, abs(int(self.xmin)), abs(int(self.xmin * 100)) % 100))
            elif t == 'inf': return(' '.join([str(x) for x in self.region]))
            else: return('/'.join([str(x) for x in self.region[:4]]))
        else: return(None)
        
    def geo_transform(self, x_inc=0, y_inc=None, node='grid'):
        """return a count info and a geotransform based on the region and a cellsize

        Args:
          x_inc (float): a x-axis gridding increment
          y_inc (float): a y-axis gridding increment

        Returns:
          list: [xcount, ycount, geot]
        """

        if y_inc is None:
            y_inc = x_inc * -1.
        elif y_inc > 0:
            y_inc = y_inc * -1.

        ## geo_transform is considered in grid-node to properly capture the region
        dst_gt = (self.xmin, x_inc, 0, self.ymax, 0, y_inc)
        this_origin = utils._geo2pixel(self.xmin, self.ymax, dst_gt, node=node)
        this_end = utils._geo2pixel(self.xmax, self.ymin, dst_gt, node=node)
        this_size = (this_end[0] - this_origin[0], this_end[1] - this_origin[1])
        return(this_end[0] - this_origin[0], this_end[1] - this_origin[1], dst_gt)

    def export_as_list(self, include_z=False, include_w=False):
        """export region as a list

        Args:
          include_z: include the z-region in the output
          include_w: include the w-region in the output

        Returns:
          list: the region values in a list
        """
        
        region_list = [self.xmin, self.xmax, self.ymin, self.ymax]
        if include_z:
            region_list.append(self.zmin)
            region_list.append(self.zmax)
        if include_w:
            region_list.append(self.wmin)
            region_list.append(self.wmax)
        return(region_list)

    def export_as_gdal_extent(self, include_z=False, include_w=False):
        """export region as a list

        Args:
          include_z: include the z-region in the output
          include_w: include the w-region in the output

        Returns:
          tuple: the region values in a list
        """

        region_list = [self.xmin, self.ymin, self.xmax, self.ymax]
        return(tuple(region_list))
    
    def export_as_wkt(self):
        """transform a region to wkt

        Returns:
          wkt: a wkt polygon geometry
        """

        eg = [[self.ymin, self.xmin],
              [self.ymin, self.xmax],
              [self.ymax, self.xmax],
              [self.ymax, self.xmin],
              [self.ymin, self.xmin]]

        return(create_wkt_polygon(eg))

    def export_as_geom(self):
        """convert a region to an OGR geometry

        Returns:
          ogr-geom: an ogr polygon geometry
        """

        if self.wkt is None:
            self.wkt = self.export_as_wkt()
            
        if self.wkt is not None:
            return(ogr.CreateGeometryFromWkt(self.wkt))
        else: return(None)

    def export_as_ogr(self, dst_ogr, dst_fmt='ESRI Shapefile', append=False):
        """convert a region string to an OGR vector

        Args:
          dst_ogr (str): destination ogr dataset pathname
          append (bool): Append to existing dataset
        """
        
        wkt = self.export_as_wkt()
        driver = ogr.GetDriverByName(dst_fmt)
        if os.path.exists(dst_ogr):
            driver.DeleteDataSource(dst_ogr)

        dst_ds = driver.CreateDataSource(dst_ogr)
        dst_lyr = dst_ds.CreateLayer(dst_ogr, geom_type = ogr.wkbPolygon)
        dst_lyr.CreateField(ogr.FieldDefn('id', ogr.OFTInteger))
        dst_feat = ogr.Feature(dst_lyr.GetLayerDefn())
        dst_feat.SetGeometryDirectly(ogr.CreateGeometryFromWkt(wkt))
        dst_feat.SetField('id', 1)
        dst_lyr.CreateFeature(dst_feat)
        dst_feat = None
        dst_ds = None

    def increments(self, x_count, y_count):
        """returns x_inc, y_inc"""
        
        x_inc = float((self.xmax - self.xmin) / x_count)
        y_inc = float((self.ymax - self.ymin) / y_count)
        return(x_inc, y_inc)
        
    def srcwin(self, geo_transform, x_count, y_count, node='grid'):
        """output the appropriate gdal srcwin for the region 
        based on the geo_transform and x/y count.

        Returns:
            tuple: the gdal srcwin (xoff, yoff, xsize, ysize)
        """

        ## geo_transform is considered in grid-node to properly capture the region
        this_origin = [0 if x < 0 else x for x in utils._geo2pixel(self.xmin, self.ymax, geo_transform, node=node)]
        this_end = [0 if x < 0 else x for x in utils._geo2pixel(self.xmax, self.ymin, geo_transform, node=node)]
        this_size = [0 if x < 0 else x for x in ((this_end[0] - this_origin[0]), (this_end[1] - this_origin[1]))]
        if this_size[0] > x_count - this_origin[0]:
            this_size[0] = x_count - this_origin[0]
            
        if this_size[1] > y_count - this_origin[1]:
            this_size[1] = y_count - this_origin[1]
            
        if this_size[0] < 0:
            this_size[0] = 0
            
        if this_size[1] < 0:
            this_size[1] = 0
            
        return(this_origin[0], this_origin[1], this_size[0], this_size[1])

    def cut(self, cut_region=None, x_inc=None, y_inc=None):
        """ cut the region based on `cut_region` while accounting for
        x_inc and y_inc so as to align to a potential grid

        Returns:
          region-object: self
        """
        
        if self.valid_p():
            if cut_region is not None and cut_region.valid_p():

                if not regions_intersect_ogr_p(cut_region, self):
                    return(self)
                
                if regions_within_ogr_p(cut_region, self):
                    return(self)
                
                if x_inc is not None and y_inc is not None:
                    x_max_count = int((self.xmax - cut_region.xmax) / x_inc)
                    y_max_count = int((self.ymax - cut_region.ymax) / y_inc)                
                    x_min_count = int((self.xmin - cut_region.xmin) / x_inc)
                    y_min_count = int((self.ymin - cut_region.ymin) / y_inc)
                    
                    self.xmax -= (x_max_count*x_inc)
                    self.ymax -= (y_max_count*y_inc)
                    self.xmin -= (x_min_count*x_inc)
                    self.ymin -= (y_min_count*y_inc)
                else:
                    return(regions_reduce(self, cut_region))
                return(self)
        return(self)
                
    def buffer(self, x_bv=0, y_bv=0, pct=None, x_inc=None, y_inc=None):
        """return the region buffered by buffer-value `bv`

        Args:
          x_bv (float): the x-direction buffer value
          y_bv (float): the y-direction buffer value
          pct (float): attain the buffer-value via percentage

        Returns:
          region-object: self
        """
        
        if self.valid_p():
            if pct is not None:
                ewp = (self.xmax - self.xmin) * (pct * .01)
                nsp = (self.ymax - self.ymin) * (pct * .01)
                x_bv = (ewp + nsp) / 2.
                y_bv = (ewp + nsp) / 2.

                if x_inc is not None:
                    if y_inc is None:
                        y_inc = x_inc
                        
                    tmp_x_bv = int(x_bv/x_inc) * float(x_inc)
                    tmp_y_bv = int(y_bv/y_inc) * float(y_inc)
                    if tmp_x_bv != 0:
                        x_bv = tmp_x_bv
                    if tmp_y_bv != 0:
                        y_bv = tmp_y_bv
                    
            self.xmin -= x_bv
            self.xmax += x_bv
            self.ymin -= y_bv
            self.ymax += y_bv

            self.wkt = self.export_as_wkt()
        return(self)

    def center(self):
        """find the center point of the xy region

        Returns:
          list: the center point [xc, yc]
        """

        if self.valid_p():
            return([self.xmin + (self.xmax-self.xmin/2),
                    self.ymax + (self.ymax-self.ymin/2)])
        else: return(None)        
        
    def chunk(self, inc, n_chunk=10):
        """chunk the xy region [xmin, xmax, ymin, ymax] into 
        n_chunk by n_chunk cell regions, given inc.

        Args:
          inc (float): the chunking increment
          n_chunk (int): number of cells

        Returns:
          list: a list of chunked regions [<regions.region>, ...]
        """

        if n_chunk is None:
            return([self])
        
        i_chunk = 0
        x_i_chunk = 0
        x_chunk = n_chunk
        o_chunks = []
        xcount, ycount, dst_gt = self.geo_transform(x_inc = inc)

        while True:
            y_chunk = n_chunk
            while True:
                this_x_origin = x_chunk - n_chunk
                this_y_origin = y_chunk - n_chunk
                this_x_size = x_chunk - this_x_origin
                this_y_size = y_chunk - this_y_origin

                c_region = Region()
                
                c_region.xmin = self.xmin + this_x_origin * inc
                c_region.xmax = c_region.xmin + this_x_size * inc
                c_region.ymin = self.ymin + this_y_origin * inc
                c_region.ymax = c_region.ymin + this_y_size * inc

                if c_region.ymax > self.ymax: c_region.ymax = self.ymax
                if c_region.ymin < self.ymin: c_region.ymin = self.ymin
                if c_region.xmax > self.xmax: c_region.xmax = self.xmax
                if c_region.xmin < self.xmin: c_region.xmin = self.xmin
                o_chunks.append(c_region)

                if y_chunk < ycount:
                    y_chunk += n_chunk
                    i_chunk += 1
                else: break
            if x_chunk < xcount:
                x_chunk += n_chunk
                x_i_chunk += 1
            else: break

        return(o_chunks)

    def warp(self, dst_srs='epsg:4326'):
        """warp the region from self.epsg to dst_epsg"""

        dst_srs = utils.str_or(dst_srs)
        if dst_srs is None or self.src_srs is None: return(self)

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
            utils.echo_error_msg('could not perform transformation')
            return(self)

        self.src_srs = dst_srs
        self.wkt = None

        src_srs = dst_srs_ = None
        return(self.transform(dst_trans))

    def transform(self, dst_trans=None):
        """transform the region using dst_trans transformation."""

        if dst_trans is None:
            utils.echo_error_msg('could not perform transformation')
            return(self)
            
        pointA = ogr.CreateGeometryFromWkt('POINT ({} {} 0)'.format(self.xmin, self.ymin))
        pointA.Transform(dst_trans)

        pointB = ogr.CreateGeometryFromWkt('POINT ({} {} 0)'.format(self.xmax, self.ymax))
        pointB.Transform(dst_trans)

        pointC = ogr.CreateGeometryFromWkt('POINT ({} {} 0)'.format(self.xmin, self.ymax))
        pointC.Transform(dst_trans)

        pointD = ogr.CreateGeometryFromWkt('POINT ({} {} 0)'.format(self.xmax, self.ymin))
        pointD.Transform(dst_trans)
        
        try:
            if not 'inf' in pointA.ExportToWkt() and not 'inf' in pointB.ExportToWkt():

                xs = [pointA.GetX(), pointB.GetX(), pointC.GetX(), pointD.GetX()]
                ys = [pointA.GetY(), pointB.GetY(), pointC.GetY(), pointD.GetY()]

                self.xmin = min(xs)
                self.xmax = max(xs)
                self.ymin = min(ys)
                self.ymax = max(ys)
                
                #self.xmin = pointA.GetX() if pointA.GetX() < pointB.GetX() else pointB.GetX()
                #self.xmax = pointB.GetX() if pointB.GetX() > pointA.GetX() else pointA.GetX()
                #self.ymin = pointA.GetY() if pointA.GetY() < pointB.GetY() else pointB.GetY()
                #self.ymax = pointB.GetY() if pointB.GetY() > pointA.GetY() else pointA.GetY()
        except Exception as e:
            sys.stderr.write('transform error: {}\n'.format(str(e)))

        dst_trans = pointA = pointB = None
        return(self)
    
## ==============================================
## do things to and with regions 
## ==============================================
def regions_reduce(region_a, region_b):
    """combine two regions and find their minimum combined region.

    if the regions don't overlap, will return an invalid region.
    check the result with valid_p()
    
    Args:
      region_a (region): a region object
      region_b (region): a region object

    Returns:
      region: the minimum region when combining `region_a` and `region_b`
    """
    
    region_c = Region()
    if region_a.valid_p() and region_b.valid_p():
        if region_a.xmin is not None and region_b.xmin is not None:
            region_c.xmin = region_a.xmin if region_a.xmin > region_b.xmin else region_b.xmin
        if region_a.xmax is not None and region_b.xmax is not None:
            region_c.xmax = region_a.xmax if region_a.xmax < region_b.xmax else region_b.xmax
        if region_a.ymin is not None and region_b.ymin is not None:
            region_c.ymin = region_a.ymin if region_a.ymin > region_b.ymin else region_b.ymin
        if region_a.ymax is not None and region_b.ymax is not None:
            region_c.ymax = region_a.ymax if region_a.ymax < region_b.ymax else region_b.ymax
        if region_a.zmin is not None and region_b.zmin is not None:
            region_c.zmin = region_a.zmin if region_a.zmin > region_b.zmin else region_b.zmin
        else:
            if region_a.zmin is not None:
                region_c.zmin = region_a.zmin
            if region_b.zmin is not None:
                region_c.zmin = region_b.zmin
        if region_a.zmax is not None and region_b.zmax is not None:
            region_c.zmax = region_a.zmax if region_a.zmax < region_b.zmax else region_b.zmax
        else:
            if region_a.zmax is not None:
                region_c.zmax = region_a.zmax
            if region_b.zmax is not None:
                region_c.zmax = region_b.zmax
        if region_a.wmin is not None and region_b.wmin is not None:
            region_c.wmin = region_a.wmin if region_a.wmin > region_b.wmin else region_b.wmin
        else:
            if region_a.wmin is not None:
                region_c.wmin = region_a.wmin
            if region_b.wmin is not None:
                region_c.wmin = region_b.wmin
        if region_a.wmax is not None and region_b.wmax is not None:
            region_c.wmax = region_a.wmax if region_a.wmax < region_b.wmax else region_b.wmax
        else:
            if region_a.wmax is not None:
                region_c.wmax = region_a.wmax
            if region_b.wmax is not None:
                region_c.wmax = region_b.wmax
    return(region_c)

def regions_merge(region_a, region_b):
    """combine two regions and find their maximum combined region.
    
    Args:
      region_a (region): a region object
      region_b (region): a region object

    Returns:
      region: the maximum region [xmin, xmax, ymin, ymax] when combining `region_a` `and region_b`
    """
    
    region_c = Region()
    if region_a.valid_p() and region_b.valid_p():
        region_c.xmin = region_a.xmin if region_a.xmin < region_b.xmin else region_b.xmin
        region_c.xmax = region_a.xmax if region_a.xmax > region_b.xmax else region_b.xmax
        region_c.ymin = region_a.ymin if region_a.ymin < region_b.ymin else region_b.ymin
        region_c.ymax = region_a.ymax if region_a.ymax > region_b.ymax else region_b.ymax
        if region_a.zmin is not None and region_b.zmin is not None:
            #if region_a.zmin is None:
            #    region_c.zmin = region_b.zmin
            #elif region_b.zmin is None:
            #    region_c.zmin = region_a.zmin
            #else:
            region_c.zmin = region_a.zmin if region_a.zmin < region_b.zmin else region_b.zmin
        if region_a.zmax is not None and region_b.zmax is not None:
            #if region_a.zmax is None:
            #    region_c.zmax = region_b.zmax
            #elif region_b.zmax is None:
            #    region_c.zmax = region_a.zmax
            #else:
            region_c.zmax = region_a.zmax if region_a.zmax > region_b.zmax else region_b.zmax
    return(region_c)

def regions_intersect_p(region_a, region_b):
    """check if two regions intersect.

    region.valid_p(regions_reduce(region_a, region_b))

    Args:
      region_a (region): a region object
      region_b (region): a region object

    Returns:
      bool: True if `region_a` and `region_b` intersect else False
    """

    if region_a.valid_p() and region_b.valid_p():
        reduced_region = regions_reduce(region_a, region_b)
        if any(reduced_region.full_region()):
            if not reduced_region.valid_p():
                return(False)
            
        if not regions_intersect_ogr_p(region_a, region_b):
            return(False)
        
    return(True)
    
def regions_intersect_ogr_p(region_a, region_b):
    """check if two regions intersect using ogr

    region_a_ogr_geom.Intersects(region_b_ogr_geom)
    
    Args:
      region_a (region): a region 
      region_b (region): a region 

    Returns:
      bool: True if `region_a` and `region_b` intersect else False.
    """

    if region_a.valid_p() and region_b.valid_p():
        if region_a.export_as_geom().Intersects(region_b.export_as_geom()):
            return(True)
    return(False)

def regions_within_ogr_p(region_a, region_b):
    """check if region_b is within region_a using ogr

    Args:
      region_a (region): a region 
      region_b (region): a region 

    Returns:
      bool: True if `region_b` is within `region_a` else False.
    """

    if region_a.valid_p() and region_b.valid_p():
        if region_b.export_as_geom().Within(region_a.export_as_geom()):
            return(True)
    return(False)

def z_region_pass(region, upper_limit=None, lower_limit=None):
    """return True if extended region [xmin, xmax, ymin, ymax, zmin, zmax] is 
    within upper and lower z limits
    
    Args:
      region (list): a long-region [xmin, xmax, ymin, ymax, zmin, zmax]
      upper_limit (float): the z-max
      lower_limit (float): the z-min
    
    Returns:
      bool: True if region z values are within upper and lower limit
    """
    
    if region is not None:
        z_region = region[4:]
        if z_region is not None and len(z_region) >= 2:
            if upper_limit is not None:
                if z_region[0] >= upper_limit:
                    return(False)
            if lower_limit is not None:
                if z_region[1] <= lower_limit:
                    return(False)
    return(True)

def xyz_in_region_p(xyz, this_region):
    """check if xyz point in inside the given region

    Args:
      xyz (xyz): an xyz data object
      this_region (region): a region object

    Returns:
      bool: True if xyz point inside region else False
    """

    pass_d = True
    xyz_wkt = xyz.export_as_wkt()
    p_geom = ogr.CreateGeometryFromWkt(xyz_wkt)
    r_geom = this_region.export_as_geom()

    if r_geom is not None:
        pass_d = p_geom.Within(r_geom)

    if pass_d:
        if this_region.zmin is not None:
            if xyz.z < this_region.zmin:
                pass_d = False
                return(pass_d)
        if this_region.zmax is not None:
            if xyz.z > this_region.zmax:
                pass_d = False
                return(pass_d)
        if this_region.wmin is not None:
            if xyz.w < this_region.wmin:
                pass_d = False
                return(pass_d)
        if this_region.wmax is not None:
            if xyz.w > this_region.wmax:
                pass_d = False
                return(pass_d)

    return(pass_d)

def gdal_ogr_regions(src_ds):
    """return the region(s) of the ogr dataset
    
    Args:
      src_ds (str): source ogr dataset pathname

    Returns:
      list: list of regions
    """
    
    these_regions = []
    if os.path.exists(src_ds):
        poly = ogr.Open(src_ds)
        if poly is not None:
            p_layer = poly.GetLayer(0)
            for pf in p_layer:
                this_region = Region()
                pgeom = pf.GetGeometryRef()
                pwkt = pgeom.ExportToWkt()
                this_region.from_list(ogr.CreateGeometryFromWkt(pwkt).GetEnvelope())
                these_regions.append(this_region)
        poly = None
        
    return(these_regions)

def create_wkt_polygon(coords, xpos=1, ypos=0):
    """convert coords to Wkt

    Args:
      coords (list): x/y geographic coords
      xpos (int): the position of the x value in coords
      ypos (int): the position of the y value in corrds

    Returns:
      wkt: polygon as wkt
    """

    ring = ogr.Geometry(ogr.wkbLinearRing)
    for coord in coords: ring.AddPoint(coord[xpos], coord[ypos])
    poly = ogr.Geometry(ogr.wkbPolygon)
    poly.AddGeometry(ring)
    poly_wkt = poly.ExportToWkt()
    poly = None
    
    return(poly_wkt)

def ogr_wkts(src_ds):
    """return the wkt(s) of the ogr dataset"""
    
    these_regions = []
    src_s = src_ds.split(':')
    if os.path.exists(src_s[0]):
        poly = ogr.Open(src_s[0])
        if poly is not None:
            p_layer = poly.GetLayer(0)
            for pf in p_layer:
                pgeom = pf.GetGeometryRef()
                pwkt = pgeom.ExportToWkt()
                r = Region().from_string(pwkt)
                if len(src_s) > 1:
                    src_r = src_s[1].split('/')
                    if len(src_r) > 0: r.zmin = utils.float_or(src_r[0])
                    if len(src_r) > 1: r.zmax = utils.float_or(src_r[1])
                    if len(src_r) > 2: r.wmin = utils.float_or(src_r[2])
                    if len(src_r) > 3:  r.wmax = utils.float_or(src_r[3])
                these_regions.append(r)

        poly = None
    return(these_regions)

## ==============================================
##
## regions cli
##
## ==============================================
regions_usage = '''{cmd} ({version}): regions; Process and generate regions

usage: {cmd} [ -hqPRW [ args ] ]...

Options:
  -R, --region\t\tThe desired REGION 
\t\t\tWhere a REGION is xmin/xmax/ymin/ymax[/zmin/zmax[/wmin/wmax]]
\t\t\tUse '-' to indicate no bounding range; e.g. -R -/-/-/-/-10/10/1/-
\t\t\tOR an OGR-compatible vector file with regional polygons. 
\t\t\tWhere the REGION is /path/to/vector[:zmin/zmax[/wmin/wmax]].
\t\t\tIf a vector file is supplied, will use each region found therein.
  -P, --s_srs\t\tSet the SOURCE projection.
  -W, --t_srs\t\tSet the TARGET projection.
  -B, --buffer\t\tBUFFER the region with a buffer-value.
  -e, --echo\t\tECHO the <processed> region
  -n, --name\t\tPrint the region as a NAME

  --quiet\t\tLower the verbosity to a quiet

  --help\t\tPrint the usage text
  --version\t\tPrint the version information

Examples:
  % {cmd} -R -90/-89/30/31

CIRES DEM home page: <http://ciresgroups.colorado.edu/coastalDEM>\
'''.format(cmd =  os.path.basename(sys.argv[0]), 
           version = cudem.__version__)

def regions_cli(argv = sys.argv):
    """run regions from command-line

    See `regions_usage` for full cli options.
    """

    src_srs = None
    dst_srs = None
    i_regions = []
    these_regions = []
    want_verbose = True
    echo = False
    echo_fn = False
    bv = None
    te = False
    
    ## ==============================================
    ## parse command line arguments.
    ## ==============================================
    i = 1
    while i < len(argv):
        arg = argv[i]
        if arg == '--region' or arg == '-R':
            i_regions.append(str(argv[i + 1]))
            i = i + 1
        elif arg[:2] == '-R':
            i_regions.append(str(arg[2:]))
        elif arg == '-s_srs' or arg == '--s_srs' or arg == '-P':
            src_srs = argv[i + 1]
            i = i + 1
        elif arg == '-t_srs' or arg == '--t_srs' or arg == '-W':
            dst_srs = argv[i + 1]
            i = i + 1
        elif arg == '-b' or arg == '-B' or arg == '--buffer':
            bv = utils.float_or(argv[i+1])
            i = i + 1
        elif arg == '-e' or arg == '--echo' or arg == '-ee':
            echo = True
        elif arg == '-te':
            echo = True
            te = True
        elif arg == '-n' or arg == '--name':
            echo_fn = True
        elif arg == '--quiet' or arg == '-q':
            want_verbose = False
        elif arg == '--help' or arg == '-h':
            print(datalists_usage)
            sys.exit(1)
        elif arg == '--version' or arg == '-v':
            print('{}, version {}'.format(os.path.basename(sys.argv[0]), cudem.__version__))
            sys.exit(1)
        elif arg[0] == '-':
            print(regions_usage)
            sys.exit(0)
        else: dls.append(arg)
        i = i + 1

    for i_region in i_regions:
        tmp_region = Region().from_string(i_region)
        if tmp_region.valid_p():
            these_regions.append(tmp_region)
        else:
            tmp_region = ogr_wkts(i_region)
            for i in tmp_region:
                if i.valid_p():
                    these_regions.append(i)
                    
    if len(these_regions) == 0:
        print(regions_usage)
        utils.echo_error_msg('you must specify at least one region')
        sys.exit(-1)
    else:
        if want_verbose: utils.echo_msg('parsed {} region(s)'.format(len(these_regions)))
    
    for rn, this_region in enumerate(these_regions):
        if src_srs is not None:
            this_region.src_srs = src_srs
            
        if dst_srs is not None:
            this_region.warp(dst_srs)
        
        if bv is not None:
            this_region.buffer(x_bv=bv, y_bv=bv)
            
        if echo:
            if te:
                print(this_region.format('te'))
            else:
                print(this_region.format('gmt'))
        elif echo_fn:
            print(this_region.format('fn'))
        else:
            this_region.export_as_ogr('region_{}.shp'.format(this_region.format('fn')))
### End
