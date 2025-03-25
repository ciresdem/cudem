### regions.py
##
## Copyright (c) 2010 - 2025 CIRES Regents of the University of Colorado
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
##
## Regions are bounding boxes made up of corner coordinates
##
### Code:

import os
import sys

import warnings
import pyproj
from osgeo import ogr

import cudem
from cudem import utils
from cudem import srsfun

ogr.DontUseExceptions()

## regions
## the region class and associated functions.
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

    full_region = lambda x: [x.xmin, x.xmax, x.ymin, x.ymax, x.zmin, x.zmax, x.wmin, x.wmax, x.umin, x.umax]
    xy_region = lambda x: [x.xmin, x.xmax, x.ymin, x.ymax]
    z_region = lambda x: [x.zmin, x.zmax]
    w_region = lambda x: [x.wmin, x.wmax]
    u_region = lambda x: [x.umin, x.umax]
    
    def __init__(self, xmin = None, xmax = None, ymin = None, ymax = None,
                 zmin = None, zmax = None, wmin = None, wmax = None,
                 umin = None, umax = None, src_srs='epsg:4326', wkt=None):        
        self.xmin = utils.float_or(xmin)
        self.xmax = utils.float_or(xmax)
        self.ymin = utils.float_or(ymin)
        self.ymax = utils.float_or(ymax)
        self.zmin = utils.float_or(zmin)
        self.zmax = utils.float_or(zmax)
        self.wmin = utils.float_or(wmin)
        self.wmax = utils.float_or(wmax)
        self.umin = utils.float_or(umin)
        self.umax = utils.float_or(umax)
        self.src_srs = src_srs
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
            if self.xmin > self.xmax:
                return(False)
            
        if self.ymin is not None and self.ymax is not None:
            if self.ymin > self.ymax:
                return(False)
            
        if self.zmin is not None and self.zmax is not None:
            if self.zmin > self.zmax:
                return(False)
            
        if self.wmin is not None and self.wmax is not None:
            if self.wmin > self.wmax:
                return(False)
            
        if self.umin is not None and self.umax is not None:
            if self.umin > self.umax:
                return(False)

        if not any(self.full_region()):
            return(False)

        return(True)

    def copy(self):
        """return a copy of the region."""
        
        return(Region(xmin=self.xmin, xmax=self.xmax,
                      ymin=self.ymin, ymax=self.ymax,
                      zmin=self.zmin, zmax=self.zmax,
                      wmin=self.wmin, wmax=self.wmax,
                      umin=self.umin, umax=self.umax,
                      src_srs=self.src_srs, wkt=self.wkt))

    def _wgs_extremes(self, just_below = False):
        """adjust the region to make sure it is within WGS extremes..."""
        
        if self.xmin <= -180:
            self.xmin = -180 if not just_below else -179.85
            
        if self.xmax >= 180:
            self.xmax = 180 if not just_below else 179.85
            
        if self.ymin <= -90:
            self.ymin = -90 if not just_below else -89.85
            
        if self.ymax >= 90:
            self.ymax = 90 if not just_below else 89.85            
    
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
                    
                    if len(region_list) >= 10:
                        self.umin = utils.float_or(region_list[8])
                        self.umax = utils.float_or(region_list[9])
                        
            if self.wkt is None:
                if self.valid_p(check_xy = True):
                    self.wkt = self.export_as_wkt()
                    
                else:
                    self.wkt = None
                    
        return(self)

    def from_string(self, region_str):
        """import a region from a region string 

        Args:
          region_str (str): <-R>xmin/xmax/ymin/ymax/zmin/zmax

        Returns:
          region-object: self
        """

        region_str = utils.str_or(region_str)
        if region_str is None:
            return(self)
        
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
            t = 'fstr': xmin/xmax/ymin/ymax/zmin/zmax/wmin/wmax/umin/umax
            t = 'sstr': xmin xmax ymin ymax
            t = 'gmt': -Rxmin/xmax/ymin/ymax
            t = 'waffles': -Rxmin/xmax/ymin/ymax/zmin/zmax/wmin/wmax/umin/umax
            t = 'bbox': xmin,ymin,xmax,ymax
            t = 'osm_bbox': ymin,xmin,ymax,xmax
            t = 'te': xmin ymin xmax ymax
            t = 'ul_lr': xmin ymax xmax ymin
            t = 'fn': ymax_xmin
            t = 'polygon': xmin,ymin,xmin,ymax,xmax,ymax,xmax,ymin,xmin,ymin

        Returns:
          str: the formatted region as str or None if region is invalid
        """

        if self.valid_p():
            if t == 'str': return('/'.join([str(self.xmin), str(self.xmax),
                                            str(self.ymin), str(self.ymax)]))
            elif t == 'sstr': return(' '.join([str(self.xmin), str(self.xmax),
                                               str(self.ymin), str(self.ymax)]))
            elif t == 'fstr': return(' '.join([str(self.xmin), str(self.xmax),
                                               str(self.ymin), str(self.ymax),
                                               str(self.zmin), str(self.zmax),
                                               str(self.wmin), str(self.wmax),
                                               str(self.umin), str(self.umax),]))
            elif t == 'gmt': return('-R{:.16f}/{:.16f}/{:.16f}/{:.16f}'.format(
                    self.xmin, self.xmax,
                    self.ymin, self.ymax))
            elif t == 'cudem': return('-R{:.16f}/{:.16f}/{:.16f}/{:.16f}/{}/{}/{}/{}'.format(
                    self.xmin, self.xmax,
                    self.ymin, self.ymax,
                    utils.float_or(self.wmin, '-'), utils.float_or(self.wmax, '-'),
                    utils.float_or(self.umin, '-'), utils.float_or(self.umax, '-')))
            elif t == 'bbox': return(','.join([str(self.xmin), str(self.ymin),
                                               str(self.xmax), str(self.ymax)]))
            elif t == 'osm_bbox': return(','.join([str(self.ymin), str(self.xmin),
                                                   str(self.ymax), str(self.xmax)]))
            elif t == 'te': return(' '.join([str(self.xmin), str(self.ymin),
                                             str(self.xmax), str(self.ymax)]))
            elif t == 'ul_lr': return(' '.join([str(self.xmin), str(self.ymax),
                                                str(self.xmax), str(self.ymin)]))
            elif t == 'fn':
                ns = 's' if self.ymax < 0 else 'n'
                ew = 'e' if self.xmin > 0 else 'w'
                return('{}{:02d}x{:02d}_{}{:03d}x{:02d}'.format(
                    ns, abs(int(self.ymax)), abs(int(self.ymax * 100)) % 100, 
                    ew, abs(int(self.xmin)), abs(int(self.xmin * 100)) % 100))
            elif t == 'fn_full':
                ns_mx = 's' if self.ymax < 0 else 'n'
                ns_mn = 's' if self.ymin < 0 else 'n'
                ew_mx = 'e' if self.xmax > 0 else 'w'
                ew_mn = 'e' if self.xmin > 0 else 'w'
                return('{}{:02d}x{:6f}_{}{:03d}x{:6f}_{}{:02d}x{:6f}_{}{:03d}x{:6f}'.format(
                    ns_mx, abs(int(self.ymax)), abs(self.ymax * 100) % 100, 
                    ew_mn, abs(int(self.xmin)), abs(self.xmin * 100) % 100,
                    ns_mn, abs(int(self.ymin)), abs(self.ymin * 100) % 100,
                    ew_mx, abs(int(self.xmax)), abs(self.xmax * 100) % 100))
            elif t == 'polygon': return(
                    '{xmin},{ymin},{xmax},{ymin},{xmax},{ymax},{xmin},{ymax},{xmin},{ymin}'.format(
                        xmin=self.xmin, ymin=self.ymin, xmax=self.xmax, ymax=self.ymax
                    )
            )
            elif t == 'inf':
                return(' '.join([str(x) for x in self.region]))
            else:
                return('/'.join([str(x) for x in self.region[:4]]))
            
        else:
            return(None)
        
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

    def export_as_list(self, include_z = False, include_w = False, include_u = False):
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
            
        if include_u:
            region_list.append(self.umin)
            region_list.append(self.umax)
            
        return(region_list)

    def export_as_gdal_extent(self):
        """export region as a list

        Returns:
          tuple: the region values in a list
        """

        region_list = [self.xmin, self.ymin, self.xmax, self.ymax]
        return(tuple(region_list))

    def export_as_polygon(self):
        """convert a region to a polygon list

        REturns:
          list: the region as a 2d polygon
        """
        
        eg = [[self.xmin, self.ymin],
              [self.xmin, self.ymax],
              [self.xmax, self.ymax],
              [self.xmax, self.ymin],
              [self.xmin, self.ymin]]

        return(eg)
    
    def export_as_wkt(self):
        """convert a region to wkt

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

    def export_as_ogr(self, dst_ogr, dst_fmt = 'ESRI Shapefile', append = False):
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
        
    def srcwin(self, geo_transform, x_count, y_count, node = 'grid'):
        """output the appropriate gdal srcwin for the region 
        based on the geo_transform and x/y count.

        Returns:
            tuple: the gdal srcwin (xoff, yoff, xsize, ysize)
        """

        ## geo_transform is considered in grid-node to properly capture the region
        this_origin = [0 if x < 0 else x for x in utils._geo2pixel(
            self.xmin, self.ymax, geo_transform, node=node
        )]
        this_end = [0 if x < 0 else x for x in utils._geo2pixel(
            self.xmax, self.ymin, geo_transform, node=node
        )]
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

    def cut(self, cut_region = None, x_inc = None, y_inc = None):
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
                
    def buffer(self, x_bv = 0, y_bv = 0, pct = None, x_inc = None, y_inc = None):
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

    def round(self, round_val=5):
        self.xmin = round(self.xmin, round_val)
        self.xmax = round(self.xmax, round_val)
        self.ymin = round(self.ymin, round_val)
        self.ymax = round(self.ymax, round_val)
    
    def center(self):
        """find the center point of the xy region

        Returns:
          list: the center point [xc, yc]
        """

        if self.valid_p():
            return([self.xmin + ((self.xmax-self.xmin)/2),
                    self.ymax + ((self.ymax-self.ymin)/2)])
        
        else:
            return(None)        
        
    def chunk(self, inc, n_chunk = 10):
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

                if c_region.ymax > self.ymax:
                    c_region.ymax = self.ymax
                    
                if c_region.ymin < self.ymin:
                    c_region.ymin = self.ymin
                    
                if c_region.xmax > self.xmax:
                    c_region.xmax = self.xmax
                    
                if c_region.xmin < self.xmin:
                    c_region.xmin = self.xmin
                    
                o_chunks.append(c_region)

                if y_chunk < ycount:
                    y_chunk += n_chunk
                    i_chunk += 1
                else:
                    break
                
            if x_chunk < xcount:
                x_chunk += n_chunk
                x_i_chunk += 1
            else:
                break

        return(o_chunks)

    def warp(self, dst_crs = 'epsg:4326', include_z = True):
        ## horizontal only for region
        if utils.str_or(self.src_srs) is None:
            utils.echo_warning_msg('region has no valid associated srs: {}'.format(self.src_srs))
            return(self)

        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            transform = srsfun.parse_srs(src_srs = self.src_srs, dst_srs = dst_crs)

            if transform['src_horz_crs'] is not None and transform['dst_horz_crs'] is not None:        
                ## horizontal Transformation
                transform['horz_pipeline'] = '+proj=pipeline +step {} +inv +step {}'.format(
                    transform['src_horz_crs'].to_proj4(), transform['dst_horz_crs'].to_proj4()
                )
                #transform['trans_region'] = region.copy()
                #transform['trans_region'].src_srs = transform['dst_horz_crs'].to_proj4()
                #transform['trans_region'].warp(transform['src_horz_crs'].to_proj4())
                transformer = pyproj.Transformer.from_pipeline(transform['horz_pipeline'])

                self.src_srs = dst_crs
                self.wkt = None
            
                return(self.transform(transformer, include_z=include_z))
            else:
                return(self)
        
        #transform = srsfun.set_transform(src_srs = self.src_srs, dst_srs = dst_crs, region = self, infos = None)
        
        # if self.src_srs.upper().startswith('EPSG'):
        #     src_srs = self.src_srs.split('+')[0]
        # else:
        #     src_srs = self.src_srs

        # if dst_crs.upper().startswith('EPSG'):
        #     dst_crs = dst_crs.split('+')[0]

        # try:
        #     in_horz_epsg = srsfun.epsg_from_input(src_srs)[0]
        #     in_crs = pyproj.CRS.from_user_input(in_horz_epsg)
        #     out_horz_epsg = srsfun.epsg_from_input(dst_crs)[0]        
        #     out_crs = pyproj.CRS.from_user_input(out_horz_epsg)
        #     transformer = pyproj.Transformer.from_crs(in_crs, out_crs, always_xy=True)
        # except:
        #     transformer = None
        
        # if transformer is None or not self.valid_p():
        #     utils.echo_error_msg('could not perform transformation')
            
        #     return(self)
        
        
    def transform(self, transformer = None, include_z = True):
        if transformer is None or not self.valid_p():
            utils.echo_error_msg('could not perform transformation')
            return(self)

        if include_z and (self.zmin is not None and self.zmax is not None):
            self.xmin, self.ymin, self.zmin = transformer.transform(self.xmin, self.ymin, self.zmin)
            self.xmax, self.ymax, self.zmax = transformer.transform(self.xmax, self.ymax, self.zmax)
        else:
            self.xmin, self.ymin = transformer.transform(self.xmin, self.ymin)
            self.xmax, self.ymax = transformer.transform(self.xmax, self.ymax)

        return(self)
         
## do things to and with regions...
## various region related functions
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

        if region_a.umin is not None and region_b.umin is not None:
            region_c.umin = region_a.umin if region_a.umin > region_b.umin else region_b.umin
        else:
            if region_a.umin is not None:
                region_c.umin = region_a.umin
                
            if region_b.umin is not None:
                region_c.umin = region_b.umin
                
        if region_a.umax is not None and region_b.umax is not None:
            region_c.umax = region_a.umax if region_a.umax < region_b.umax else region_b.umax
        else:
            if region_a.umax is not None:
                region_c.umax = region_a.umax
                
            if region_b.umax is not None:
                region_c.umax = region_b.umax
                
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
            region_c.zmin = region_a.zmin if region_a.zmin < region_b.zmin else region_b.zmin
            
        if region_a.zmax is not None and region_b.zmax is not None:
            region_c.zmax = region_a.zmax if region_a.zmax > region_b.zmax else region_b.zmax
    ## add w and u
            
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

def z_region_pass(region, upper_limit = None, lower_limit = None):
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

    # add w and u
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
            
        if this_region.umin is not None:
            if xyz.u < this_region.umin:
                pass_d = False
                return(pass_d)
            
        if this_region.umax is not None:
            if xyz.u > this_region.umax:
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

def write_shapefile(geom, out_shp):
    driver = ogr.GetDriverByName('Esri Shapefile')
    ds = driver.CreateDataSource(out_shp)
    layer = ds.CreateLayer('', None, ogr.wkbPolygon)
    
    # Add one attribute
    layer.CreateField(ogr.FieldDefn('id', ogr.OFTInteger))
    defn = layer.GetLayerDefn()

    # Create a new feature (attribute and geometry)
    feat = ogr.Feature(defn)
    feat.SetField('id', 123)

    # Make a geometry, from Shapely object
    #geom = ogr.CreateGeometryFromWkb(poly.wkb)
    feat.SetGeometry(geom)

    layer.CreateFeature(feat)
    feat = geom = None  # destroy these

    # Save and close everything
    ds = layer = feat = geom = None

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

def parse_cli_region(region_list, verbose = True, pct_buffer = None):
    """parse a region list into region(s).

    for use in clis

    region_list is a list of regions, where each region is either:
    bbox str ( 'xmin/xmax/ymin/ymax' ),
    tile_set (with options of 'in_region' and 'inc' to define the tiles)
    path to an ogr dataset.

    returns a list of region objects.
    """

    these_regions = []
    for i_region in region_list:
        if i_region is None:
            these_regions.append(None)
            continue
        
        i_region_s = i_region.split(':')
        tmp_region = Region().from_string(i_region_s[0])
        args = utils.args2dict(i_region_s[1:], {})
        if tmp_region.valid_p(check_xy=True):
            if 'pct_buffer' in args.keys():
                tmp_region.buffer(pct=utils.float_or(args['pct_buffer']))
                
            these_regions.append(tmp_region)
        elif str(i_region_s[0]) == 'tile_set':
            these_regions = generate_tile_set(**args)
        else:
            tmp_region = ogr_wkts(i_region_s[0])
            for i in tmp_region:
                if i.valid_p():
                    if len(i_region_s) > 1:
                        region_extender = ''
                        for opts in i_region_s[1:]:
                            if '=' not in opts:
                                region_extender = opts
                                
                        #this_region = Region().from_string('/'.join([i.format('str'), i_region_s[1]]))
                        this_region = Region().from_string('/'.join([i.format('str'), region_extender]))
                        if 'pct_buffer' in args.keys():
                            this_region.buffer(pct=utils.float_or(args['pct_buffer']))
                            
                        these_regions.append(this_region)
                    else:
                        these_regions.append(i)

    if verbose:
        if len(these_regions) > 4:
            utils.echo_msg(
                'parsed {} region(s): {}...{}'.format(len(these_regions), these_regions[:2], these_regions[-2:])
            )
        elif len(these_regions) > 0:
            utils.echo_msg(
                'parsed {} region(s): {}'.format(len(these_regions), these_regions)
            )
        else:
            utils.echo_warning_msg(
                'failed to parse region(s), {}'.format(region_list)
            )
            
    return(these_regions)

def generate_tile_set(in_region = None, inc = .25, pct_buffer = None):
    """Generate a tile-set based on `in_region` and `inc`.

    returns a list of regions
    """
    
    tile_regions = []
    if in_region is not None:
        tmp_region = Region().from_string(in_region)
        if not tmp_region.valid_p(check_xy=True):
            tmp_region = Region(xmin=-180, xmax=180, ymin=-90, ymax=90)
    else:
        tmp_region = Region(xmin=-180, xmax=180, ymin=-90, ymax=90)
        
    inc = utils.float_or(inc, .25)
    this_xmin = tmp_region.xmin
    this_xmax = tmp_region.xmin+inc
    this_ymin = tmp_region.ymin
    this_ymax = tmp_region.ymin+inc
    while this_xmax <= tmp_region.xmax:
        while this_ymax <= tmp_region.ymax: 
            this_region = Region(xmin=this_xmin, xmax=this_xmax, ymin=this_ymin, ymax=this_ymax)
            if pct_buffer is not None:
                this_region.buffer(pct=utils.float_or(pct_buffer))
                
            tile_regions.append(this_region)
            this_ymin = this_ymax
            this_ymax += inc
            
        this_ymin = tmp_region.ymin
        this_ymax = this_ymin + inc
        this_xmin = this_xmax
        this_xmax += inc
            
    return(tile_regions)

def region_list_to_ogr(region_list, dst_ogr, dst_fmt = 'ESRI Shapefile'):
    """convert a region list to an OGR vector
    """

    driver = ogr.GetDriverByName(dst_fmt)
    if os.path.exists(dst_ogr):
        driver.DeleteDataSource(dst_ogr)

    dst_ds = driver.CreateDataSource(dst_ogr)
    dst_lyr = dst_ds.CreateLayer(dst_ogr, geom_type = ogr.wkbPolygon)
    dst_lyr.CreateField(ogr.FieldDefn('id', ogr.OFTInteger))
        
    wkts = []
    for this_id, this_region in enumerate(region_list):
        eg = [[this_region.ymin, this_region.xmin],
              [this_region.ymin, this_region.xmax],
              [this_region.ymax, this_region.xmax],
              [this_region.ymax, this_region.xmin],
              [this_region.ymin, this_region.xmin]]
    
        wkt = create_wkt_polygon(eg)    
        dst_feat = ogr.Feature(dst_lyr.GetLayerDefn())
        dst_feat.SetGeometryDirectly(ogr.CreateGeometryFromWkt(wkt))
        dst_feat.SetField('id', this_id)
        dst_lyr.CreateFeature(dst_feat)
        dst_feat = None
        
    dst_ds = None

## ==============================================
## Command-line Interface (CLI)
## $regions
##
## regions cli
## ==============================================
regions_usage = '''{cmd} ({version}): regions; Process and generate regions

usage: {cmd} [ -hqJPRT [ args ] ]...

Options:
  -R, --region\t\tThe desired REGION 
\t\t\tWhere a REGION is xmin/xmax/ymin/ymax[/zmin/zmax[/wmin/wmax/umin/umax]]
\t\t\tUse '-' to indicate no bounding range; e.g. -R -/-/-/-/-10/10/1/-/-/-
\t\t\tOR an OGR-compatible vector file with regional polygons. 
\t\t\tWhere the REGION is /path/to/vector[:zmin/zmax[/wmin/wmax/umin/umax]].
\t\t\tIf a vector file is supplied, will use each region found therein.
\t\t\tOptionally, append `:pct_buffer=<value>` to buffer the region(s) by a percentage.
  -J, --s_srs\t\tSet the SOURCE projection.
  -P, --t_srs\t\tSet the TARGET projection.
  -T, --tile_set\tGenerate a TILESET from the input region. (set incrememnt here)
  -B, --buffer\t\tBUFFER the region with a buffer-value.
  -e, --echo\t\tECHO the <processed> region
  -n, --name\t\tPrint the region as a NAME
  -m, --merge\t\tMERGE all regions into a single region

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
    tile_set = None
    echo = False
    echo_fn = False
    want_merge = False
    bv = None
    te = False
    
    ## parse command line arguments.
    i = 1
    while i < len(argv):
        arg = argv[i]
        if arg == '--region' or arg == '-R':
            i_regions.append(str(argv[i + 1]))
            i = i + 1
        elif arg[:2] == '-R':
            i_regions.append(str(arg[2:]))
        elif arg == '-s_srs' or arg == '--s_srs' or arg == '-J':
            src_srs = argv[i + 1]
            i = i + 1
        elif arg == '-t_srs' or arg == '--t_srs' or arg == '-P':
            dst_srs = argv[i + 1]
            i = i + 1
        elif arg == '-b' or arg == '-B' or arg == '--buffer':
            bv = utils.float_or(argv[i+1])
            i = i + 1
        elif arg == '-t' or arg == '-T' or arg == '--tile_set':
            tile_set = utils.float_or(argv[i+1])
            i = i + 1
        elif arg == '-e' or arg == '--echo' or arg == '-ee':
            echo = True
        elif arg == '-te':
            echo = True
            te = True
        elif arg == '-n' or arg == '--name':
            echo_fn = True
        elif arg == '-m' or arg == '--merge':
            want_merge = True
        elif arg == '--quiet' or arg == '-q':
            want_verbose = False
        elif arg == '--help' or arg == '-h':
            print(regions_usage)
            sys.exit(1)
        elif arg == '--version' or arg == '-v':
            print('{}, version {}'.format(os.path.basename(sys.argv[0]), cudem.__version__))
            sys.exit(1)
        elif arg[0] == '-':
            print(regions_usage)
            sys.exit(0)
            
        i = i + 1

    ## parse the input region(s) to a list
    these_regions = parse_cli_region(i_regions)
    if not these_regions:
        print(regions_usage)
        utils.echo_error_msg('you must specify at least one region')
        sys.exit(-1)

    ## generate a tile-set from the input region(s)
    if tile_set is not None:
        for rn, this_region in enumerate(these_regions):
            these_tiles = generate_tile_set(this_region.format('gmt'), tile_set)
            region_list_to_ogr(these_tiles, 'regions_tile_set.shp')
            
        these_regions = []

    ## merge the input regions together before further processing,
    ## set the merged region as the default region
    if want_merge:
        region_cnt = len(these_regions)
        merged_region = these_regions[0].copy()
        for r in range(1, region_cnt):
            merged_region = regions_merge(merged_region, these_regions[r])

        these_regions = [merged_region]
        
    for rn, this_region in enumerate(these_regions):
        ## set the source srs if specified
        if src_srs is not None:
            this_region.src_srs = src_srs

        ## warp the region to dst_srs if specified
        if dst_srs is not None:
            this_region.warp(dst_srs)

        ## buffer the region to buffer-value `bv`
        if bv is not None:
            this_region.buffer(x_bv=bv, y_bv=bv)

        ## echo the region or file-name to stdout else export as a shapefile
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
