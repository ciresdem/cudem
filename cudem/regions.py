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
###  Code:

import os
import sys
import warnings
import math
import argparse
from typing import List, Optional, Tuple, Union, Any, Dict

import pyproj
from osgeo import ogr

from cudem import utils
from cudem import factory
from cudem import srsfun
from cudem import __version__ as __cudem_version__

## Suppress OGR exceptions to maintain legacy behavior
ogr.DontUseExceptions()

__version__ = '1.0.0'

class Region:
    """Representing a geographic region.
    
    Attributes:
      xmin (float): the minimum x(long) value
      xmax (float): the maximum x(long) value
      ymin (float): the minimum y(lat) value
      ymax (float): the maximum y(lat) value
      zmin (float): the minimum z(elev) value
      zmax (float): the maximum z(elev) value
      wmin (float): the minimum w(weight) value
      wmax (float): the maximum w(weight) value
      wkt (str): the wkt representation of the region
      src_srs (str): the regions projection
    """

    def __init__(self, xmin=None, xmax=None, ymin=None, ymax=None,
                 zmin=None, zmax=None, wmin=None, wmax=None,
                 umin=None, umax=None, src_srs='epsg:4326', wkt=None):        
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
        return str(self.format('fstr'))

    
    def __str__(self):
        return str(self.format('fstr'))

    
    @property
    def full_region(self):
        return [self.xmin, self.xmax, self.ymin, self.ymax,
                self.zmin, self.zmax, self.wmin, self.wmax,
                self.umin, self.umax]

    
    @property
    def xy_region(self):
        return [self.xmin, self.xmax, self.ymin, self.ymax]

    @property
    def xy_extent(self):
        return [self.xmin, self.ymin, self.xmax, self.ymax]
    
    @property
    def z_region(self):
        return [self.zmin, self.zmax]

    
    @property
    def w_region(self):
        return [self.wmin, self.wmax]

    
    @property
    def u_region(self):
        return [self.umin, self.umax]

    
    def is_valid(self, check_xy: bool = False) -> bool:
        """Check the validity of a region.
    
        Returns:
          bool: True if region appears to be valid, else False.
        """
        
        if check_xy:
            if any(v is None for v in [self.xmin, self.xmax, self.ymin, self.ymax]):
                return False

        if self.xmin is not None and self.xmax is not None and self.xmin >= self.xmax:
            return False
            
        if self.ymin is not None and self.ymax is not None and self.ymin >= self.ymax:
            return False
            
        if self.zmin is not None and self.zmax is not None and self.zmin > self.zmax:
            return False
            
        if self.wmin is not None and self.wmax is not None and self.wmin > self.wmax:
            return False
            
        if self.umin is not None and self.umax is not None and self.umin > self.umax:
            return False

        # Ensure at least one value is set
        if not any(v is not None for v in self.full_region):
            return False

        return True

    ## Alias for backward compatibility
    valid_p = is_valid
    
    def copy(self):
        """Return a deep copy of the region."""
        
        return Region(
            xmin=self.xmin, xmax=self.xmax,
            ymin=self.ymin, ymax=self.ymax,
            zmin=self.zmin, zmax=self.zmax,
            wmin=self.wmin, wmax=self.wmax,
            umin=self.umin, umax=self.umax,
            src_srs=self.src_srs, wkt=self.wkt
        )

    
    def _wgs_extremes(self, just_below: bool = False):
        """Adjust the region to make sure it is within WGS extremes."""
        
        if self.xmin is not None and self.xmin <= -180:
            self.xmin = -180 if not just_below else -179.85
        if self.xmax is not None and self.xmax >= 180:
            self.xmax = 180 if not just_below else 179.85
        if self.ymin is not None and self.ymin <= -90:
            self.ymin = -90 if not just_below else -89.85
        if self.ymax is not None and self.ymax >= 90:
            self.ymax = 90 if not just_below else 89.85            

            
    def from_user_input(self, region_like):
        if isinstance(region_like, list):
            return self.from_list(region_like)
        elif isinstance(region_like, str):
            return self.from_string(region_like)
        elif isinstance(region_like, Region):
            return self.from_region(region_like)
        return self

    
    def from_list(self, region_list: List[float]):
        """Import a region from a list.

        Args:
          region_list (list): [xmin, xmax, ymin, ymax[, zmin, zmax[, wmin, wmax]]]
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
                if self.is_valid(check_xy=True):
                    self.wkt = self.export_as_wkt()
                else:
                    self.wkt = None
                    
        return self

    
    def from_string(self, region_str: str):
        """Import a region from a region string.

        Args:
          region_str (str): <-R>xmin/xmax/ymin/ymax/zmin/zmax
        """
        
        region_str = utils.str_or(region_str)
        if region_str is None:
            return self
        
        if region_str.startswith('-R'):
            region_str = region_str[2:]
            
        str_list = region_str.strip().split('/')
        if len(str_list) >= 4:
            r_list = [utils.float_or(x) for x in str_list]
            self.from_list(r_list)
        elif region_str.startswith("POLYGON") or region_str.startswith("MULTIPOLYGON"):
            self.wkt = region_str
            self.from_list(ogr.CreateGeometryFromWkt(region_str).GetEnvelope())

        return self

    
    def from_geo_transform(self, geo_transform: List[float] = None, 
                           x_count: int = None, y_count: int = None):
        """Import a region from a geotransform list and dimensions."""
        
        if geo_transform is not None and x_count is not None and y_count is not None:            
            self.xmin = geo_transform[0]
            self.xmax = geo_transform[0] + geo_transform[1] * x_count
            self.ymin = geo_transform[3] + geo_transform[5] * y_count
            self.ymax = geo_transform[3]
            
        if self.wkt is None:
            self.wkt = self.export_as_wkt()
            
        return self

    
    def from_region(self, input_region):
        """Import a region from another region."""
        
        return input_region.copy()

    
    def format(self, t: str = 'gmt') -> Optional[str]:
        """Format region to string, defined by `t`."""
        
        if not self.is_valid():
            return None

        if t == 'str': 
            return '/'.join(map(str, [self.xmin, self.xmax, self.ymin, self.ymax]))
        
        elif t == 'sstr': 
            return ' '.join(map(str, [self.xmin, self.xmax, self.ymin, self.ymax]))
        
        elif t == 'fstr': 
            return ' '.join(map(str, [
                self.xmin, self.xmax, self.ymin, self.ymax,
                self.zmin, self.zmax, self.wmin, self.wmax,
                self.umin, self.umax
            ]))
        
        elif t == 'gmt': 
            return f'-R{self.xmin:.16f}/{self.xmax:.16f}/{self.ymin:.16f}/{self.ymax:.16f}'
        
        elif t == 'cudem':
            return (f'-R{self.xmin:.16f}/{self.xmax:.16f}/{self.ymin:.16f}/{self.ymax:.16f}/'
                    f'{utils.float_or(self.wmin, "-")}/{utils.float_or(self.wmax, "-")}/'
                    f'{utils.float_or(self.umin, "-")}/{utils.float_or(self.umax, "-")}')
        
        elif t == 'bbox': 
            return ','.join(map(str, [self.xmin, self.ymin, self.xmax, self.ymax]))
        
        elif t == 'osm_bbox': 
            return ','.join(map(str, [self.ymin, self.xmin, self.ymax, self.xmax]))
        
        elif t == 'te': 
            return ' '.join(map(str, [self.xmin, self.ymin, self.xmax, self.ymax]))
        
        elif t == 'ul_lr': 
            return ' '.join(map(str, [self.xmin, self.ymax, self.xmax, self.ymin]))
        
        elif t == 'fn':
            ns = 's' if self.ymax < 0 else 'n'
            ew = 'e' if self.xmin > 0 else 'w'
            return (f'{ns}{abs(int(self.ymax)):02d}x{abs(int(self.ymax * 100)) % 100:02d}_'
                    f'{ew}{abs(int(self.xmin)):03d}x{abs(int(self.xmin * 100)) % 100:02d}')
        
        elif t == 'fn_full':
            ns_mx = 's' if self.ymax < 0 else 'n'
            ns_mn = 's' if self.ymin < 0 else 'n'
            ew_mx = 'e' if self.xmax > 0 else 'w'
            ew_mn = 'e' if self.xmin > 0 else 'w'
            return (f'{ns_mx}{abs(int(self.ymax)):02d}x{abs(self.ymax * 100) % 100:6f}_'
                    f'{ew_mn}{abs(int(self.xmin)):03d}x{abs(self.xmin * 100) % 100:6f}_'
                    f'{ns_mn}{abs(int(self.ymin)):02d}x{abs(self.ymin * 100) % 100:6f}_'
                    f'{ew_mx}{abs(int(self.xmax)):03d}x{abs(self.xmax * 100) % 100:6f}')
        
        elif t == 'polygon':
            return (f'{self.xmin},{self.ymin},{self.xmin},{self.ymax},'
                    f'{self.xmax},{self.ymax},{self.xmax},{self.ymin},{self.xmin},{self.ymin}')
        
        elif t == 'inf':
            return ' '.join(str(x) for x in self.full_region)
        
        else:
            return '/'.join(str(x) for x in self.xy_region)

        
    def geo_transform(self, x_inc: float = 0, y_inc: float = None, node: str = 'grid'):
        """Return dimensions and a geotransform based on the region and a cellsize.

        Returns:
          list: [xcount, ycount, geot]
        """
        
        if y_inc is None:
            y_inc = -x_inc
        elif y_inc > 0:
            y_inc = -y_inc

        dst_gt = (self.xmin, x_inc, 0, self.ymax, 0, y_inc)
        this_origin = utils._geo2pixel(self.xmin, self.ymax, dst_gt, node=node)
        this_end = utils._geo2pixel(self.xmax, self.ymin, dst_gt, node=node)
        
        return this_end[0] - this_origin[0], this_end[1] - this_origin[1], dst_gt

    def geo_transform_from_count(self, x_count: int = 0, y_count: int = 0):
        x_inc = (self.xmax - self.xmin) / x_count
        y_inc = (self.ymin - self.ymax) / y_count
        dst_gt = (self.xmin, x_inc, 0, self.ymax, 0, y_inc)
        return dst_gt

    
    def export_as_list(self, include_z=False, include_w=False, include_u=False):
        region_list = [self.xmin, self.xmax, self.ymin, self.ymax]
        if include_z:
            region_list.extend([self.zmin, self.zmax])
        if include_w:
            region_list.extend([self.wmin, self.wmax])
        if include_u:
            region_list.extend([self.umin, self.umax])
        return region_list

    
    def export_as_gdal_extent(self):
        """Convert a region to a gdal compatible aoi."""
        
        return (self.xmin, self.ymin, self.xmax, self.ymax)

    
    def export_as_polygon(self):
        """Convert a region to a polygon list."""
        
        return [[self.xmin, self.ymin],
                [self.xmin, self.ymax],
                [self.xmax, self.ymax],
                [self.xmax, self.ymin],
                [self.xmin, self.ymin]]

    
    def export_as_wkt(self, flatten=True):
        """Convert a region to WKT."""
        
        eg = [[self.ymin, self.xmin],
              [self.ymin, self.xmax],
              [self.ymax, self.xmax],
              [self.ymax, self.xmin],
              [self.ymin, self.xmin]]
        return create_wkt_polygon(eg, flatten=flatten)

    
    def export_as_geom(self, flatten=True):
        """Convert a region to an OGR geometry."""
        
        if self.wkt is None:
            self.wkt = self.export_as_wkt(flatten=flatten)
            
        if self.wkt is not None:
            return ogr.CreateGeometryFromWkt(self.wkt)
        return None

    
    def export_as_ogr(self, dst_ogr: str, dst_fmt: str = 'geojson', append: bool = False):
        """Convert a region to an OGR vector file."""
        
        wkt = self.export_as_wkt()
        driver = ogr.GetDriverByName(dst_fmt)
        
        if os.path.exists(dst_ogr):
            driver.DeleteDataSource(dst_ogr)

        srs = None
        if self.src_srs is not None:
            srs = srsfun.osr_srs(srsfun.osr_wkt(self.src_srs))
            
        dst_ds = driver.CreateDataSource(dst_ogr)
        dst_lyr = dst_ds.CreateLayer(dst_ogr, srs, geom_type=ogr.wkbPolygon)
        dst_lyr.CreateField(ogr.FieldDefn('id', ogr.OFTInteger))
        
        dst_feat = ogr.Feature(dst_lyr.GetLayerDefn())
        dst_feat.SetGeometryDirectly(ogr.CreateGeometryFromWkt(wkt))
        dst_feat.SetField('id', 1)
        dst_lyr.CreateFeature(dst_feat)
        
        ## Cleanup
        dst_feat = None
        dst_ds = None

        
    def increments(self, x_count, y_count):
        x_inc = float((self.xmax - self.xmin) / x_count)
        y_inc = float((self.ymax - self.ymin) / y_count)
        return x_inc, y_inc

    
    def srcwin(self, geo_transform, x_count, y_count, node='grid'):
        """Output the appropriate GDAL srcwin (xoff, yoff, xsize, ysize)."""
        
        ## Calculate pixel coordinates
        this_origin = [max(0, x) for x in utils._geo2pixel(
            self.xmin, self.ymax, geo_transform, node=node
        )]
        this_end = [max(0, x) for x in utils._geo2pixel(
            self.xmax, self.ymin, geo_transform, node=node
        )]
        
        this_size = [
            max(0, this_end[0] - this_origin[0]),
            max(0, this_end[1] - this_origin[1])
        ]
        
        ## Clamp to available size
        if this_size[0] > x_count - this_origin[0]:
            this_size[0] = x_count - this_origin[0]
            
        if this_size[1] > y_count - this_origin[1]:
            this_size[1] = y_count - this_origin[1]
            
        return this_origin[0], this_origin[1], this_size[0], this_size[1]

    
    def cut(self, cut_region=None, x_inc=None, y_inc=None):
        """Cut the region based on `cut_region` while accounting for grid alignment."""
        
        if self.is_valid() and cut_region is not None and cut_region.is_valid():
            if not regions_intersect_ogr_p(cut_region, self):
                return self
            
            if regions_within_ogr_p(cut_region, self):
                return self
            
            if x_inc is not None and y_inc is not None:
                x_max_count = int((self.xmax - cut_region.xmax) / x_inc)
                y_max_count = int((self.ymax - cut_region.ymax) / y_inc)                
                x_min_count = int((self.xmin - cut_region.xmin) / x_inc)
                y_min_count = int((self.ymin - cut_region.ymin) / y_inc)
                
                self.xmax -= (x_max_count * x_inc)
                self.ymax -= (y_max_count * y_inc)
                self.xmin -= (x_min_count * x_inc)
                self.ymin -= (y_min_count * y_inc)
            else:
                return regions_reduce(self, cut_region)
            
        return self

    
    def buffer(self, x_bv=0, y_bv=0, pct=None, x_inc=None, y_inc=None, check_if_valid=True):
        """Buffer the region by a value or percentage."""
        
        if check_if_valid and not self.is_valid():
            return self

        if pct is not None:
            ewp = (self.xmax - self.xmin) * (pct * 0.01)
            nsp = (self.ymax - self.ymin) * (pct * 0.01)
            x_bv = (ewp + nsp) / 2.0
            y_bv = (ewp + nsp) / 2.0

        if x_inc is not None:
            if y_inc is None:
                y_inc = x_inc
                    
            tmp_x_bv = int(x_bv / x_inc) * float(x_inc)
            tmp_y_bv = int(y_bv / y_inc) * float(y_inc)
            if tmp_x_bv != 0:
                x_bv = tmp_x_bv
            if tmp_y_bv != 0:
                y_bv = tmp_y_bv

        self.xmin -= x_bv
        self.xmax += x_bv
        self.ymin -= y_bv
        self.ymax += y_bv

        self.wkt = self.export_as_wkt()
        return self

    
    def round(self, round_val=5):
        self.xmin = round(self.xmin, round_val)
        self.xmax = round(self.xmax, round_val)
        self.ymin = round(self.ymin, round_val)
        self.ymax = round(self.ymax, round_val)

        
    def center(self):
        """Find the center point of the xy region."""
        
        if self.is_valid():
            return [(self.xmax + self.xmin) / 2, (self.ymax + self.ymin) / 2]
        return None        

    def chunk(self, inc, n_chunk=10):
        """Chunk the region into n_chunk by n_chunk cell regions."""
        
        if n_chunk is None:
            return [self]
        
        o_chunks = []
        xcount, ycount, _ = self.geo_transform(x_inc=inc)
        
        x_chunk = n_chunk
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

                ## Clamp values
                c_region.ymax = min(c_region.ymax, self.ymax)
                c_region.ymin = max(c_region.ymin, self.ymin)
                c_region.xmax = min(c_region.xmax, self.xmax)
                c_region.xmin = max(c_region.xmin, self.xmin)
                    
                o_chunks.append(c_region)

                if y_chunk < ycount:
                    y_chunk += n_chunk
                else:
                    break
                
            if x_chunk < xcount:
                x_chunk += n_chunk
            else:
                break

        return o_chunks

    
    def warp(self, dst_crs='epsg:4326', include_z=True):
        """Transform region horizontally to a new CRS."""
        
        if utils.str_or(self.src_srs) is None:
            utils.echo_warning_msg(f'Region has no valid associated srs: {self.src_srs}')
            return self

        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            transform = srsfun.parse_srs(src_srs=self.src_srs, dst_srs=dst_crs)

            if transform['src_horz_crs'] and transform['dst_horz_crs']:        
                pipeline_str = '+proj=pipeline +step {} +inv +step {}'.format(
                    transform['src_horz_crs'].to_proj4(), 
                    transform['dst_horz_crs'].to_proj4()
                )
                transformer = pyproj.Transformer.from_pipeline(pipeline_str)

                self.src_srs = dst_crs
                self.wkt = None
                return self.transform(transformer, include_z=include_z)
            else:
                return self

            
    def transform(self, transformer=None, include_z=True):
        if transformer is None or not self.is_valid():
            utils.echo_error_msg(f'Could not perform region transformation; {self}')
            return self

        if include_z and (self.zmin is not None and self.zmax is not None):
            self.xmin, self.ymin, self.zmin = transformer.transform(self.xmin, self.ymin, self.zmin)
            self.xmax, self.ymax, self.zmax = transformer.transform(self.xmax, self.ymax, self.zmax)
        else:
            self.xmin, self.ymin = transformer.transform(self.xmin, self.ymin)
            self.xmax, self.ymax = transformer.transform(self.xmax, self.ymax)

        return self


def regions_reduce(region_a: Region, region_b: Region) -> Region:
    """Combine two regions and find their minimum overlapping region."""
    
    region_c = Region()
    
    if region_a.is_valid() and region_b.is_valid():
        ## Helper to get the inner bounds (max of mins, min of maxes)
        def _get_overlap(val_a, val_b, func):
            if val_a is not None and val_b is not None:
                return func(val_a, val_b)
            return None

        ## Helper to merge optional fields (using first available)
        def _merge_optional(val_a, val_b, func):
            if val_a is not None and val_b is not None:
                return func(val_a, val_b)
            return val_a if val_a is not None else val_b

        region_c.xmin = _get_overlap(region_a.xmin, region_b.xmin, max)
        region_c.xmax = _get_overlap(region_a.xmax, region_b.xmax, min)
        region_c.ymin = _get_overlap(region_a.ymin, region_b.ymin, max)
        region_c.ymax = _get_overlap(region_a.ymax, region_b.ymax, min)

        region_c.zmin = _merge_optional(region_a.zmin, region_b.zmin, max)
        region_c.zmax = _merge_optional(region_a.zmax, region_b.zmax, min)
        
        region_c.wmin = _merge_optional(region_a.wmin, region_b.wmin, max)
        region_c.wmax = _merge_optional(region_a.wmax, region_b.wmax, min)
        
        region_c.umin = _merge_optional(region_a.umin, region_b.umin, max)
        region_c.umax = _merge_optional(region_a.umax, region_b.umax, min)

    return region_c


def regions_merge(region_a: Region, region_b: Region) -> Region:
    """Combine two regions and find their maximum combined extent."""
    
    region_c = Region()
    if region_a.is_valid() and region_b.is_valid():
        region_c.xmin = min(region_a.xmin, region_b.xmin)
        region_c.xmax = max(region_a.xmax, region_b.xmax)
        region_c.ymin = min(region_a.ymin, region_b.ymin)
        region_c.ymax = max(region_a.ymax, region_b.ymax)
        
        if region_a.zmin is not None and region_b.zmin is not None:
            region_c.zmin = min(region_a.zmin, region_b.zmin)
        if region_a.zmax is not None and region_b.zmax is not None:
            region_c.zmax = max(region_a.zmax, region_b.zmax)
            
    return region_c


def regions_intersect_p(region_a: Region, region_b: Region) -> bool:
    """Check if two regions intersect."""
    
    if region_a.is_valid() and region_b.is_valid():
        reduced_region = regions_reduce(region_a, region_b)
        ## Check if the overlap is valid
        if any(reduced_region.full_region) and not reduced_region.is_valid():
            return False
            
        if not regions_intersect_ogr_p(region_a, region_b):
            return False
        
    return True


def regions_intersect_ogr_p(region_a: Region, region_b: Region) -> bool:
    """Check if two regions intersect using OGR geometry."""
    
    if region_a.is_valid() and region_b.is_valid():
        geom_a = region_a.export_as_geom()
        geom_b = region_b.export_as_geom()
        if geom_a and geom_b and geom_a.Intersects(geom_b):
            return True
    return False


def regions_within_ogr_p(region_a: Region, region_b: Region) -> bool:
    """Check if region_b is within region_a using OGR."""
    
    if region_a.is_valid() and region_b.is_valid():
        geom_a = region_a.export_as_geom()
        geom_b = region_b.export_as_geom()
        if geom_a and geom_b and geom_b.Within(geom_a):
            return True
    return False


def create_wkt_polygon(coords, xpos=1, ypos=0, flatten=True):
    """Convert coordinate list to WKT Polygon."""
    
    ring = ogr.Geometry(ogr.wkbLinearRing)
    for coord in coords:
        ring.AddPoint(coord[xpos], coord[ypos])
        
    poly = ogr.Geometry(ogr.wkbPolygon)
    poly.AddGeometry(ring)
    
    if flatten:
        poly.FlattenTo2D()
        
    return poly.ExportToWkt()


def ogr_wkts(src_ds: str) -> List[Region]:
    """Return the regions found in an OGR dataset."""
    
    these_regions = []
    f_no_geom = 0
    src_s = src_ds.split(':')
    
    if os.path.exists(src_s[0]):
        ds = ogr.Open(src_s[0])
        if ds is not None:
            layer = ds.GetLayer(0)
            for feat in layer:
                geom = feat.GetGeometryRef()
                if geom is not None:
                    wkt = geom.ExportToWkt()
                    r = Region().from_string(wkt)
                    
                    ## Apply optional Z/W limits from string if present
                    if len(src_s) > 1:
                        src_r = src_s[1].split('/')
                        if len(src_r) > 0: r.zmin = utils.float_or(src_r[0])
                        if len(src_r) > 1: r.zmax = utils.float_or(src_r[1])
                        if len(src_r) > 2: r.wmin = utils.float_or(src_r[2])
                        if len(src_r) > 3: r.wmax = utils.float_or(src_r[3])

                    these_regions.append(r)
                else:
                    f_no_geom += 1

    if f_no_geom > 0:
        utils.echo_warning_msg(f'{f_no_geom} features had no geometry')
        
    return these_regions


def generate_tile_set(in_region=None, inc=0.25, pct_buffer=None):
    """Generate a tile-set based on `in_region` and `inc`."""
    
    tile_regions = []
    
    if in_region is not None:
        tmp_region = Region().from_string(in_region)
        if not tmp_region.is_valid(check_xy=True):
            tmp_region = Region(xmin=-180, xmax=180, ymin=-90, ymax=90)
    else:
        tmp_region = Region(xmin=-180, xmax=180, ymin=-90, ymax=90)
        
    inc = utils.float_or(inc, 0.25)
    this_xmin = tmp_region.xmin
    this_ymin = tmp_region.ymin
    
    while this_xmin < tmp_region.xmax:
        this_xmax = this_xmin + inc
        current_y = this_ymin
        
        while current_y < tmp_region.ymax: 
            this_ymax = current_y + inc
            this_region = Region(
                xmin=this_xmin,
                xmax=this_xmax,
                ymin=current_y,
                ymax=this_ymax
            )
            if pct_buffer is not None:
                this_region.buffer(pct=utils.float_or(pct_buffer))
                
            tile_regions.append(this_region)
            current_y = this_ymax
            
        this_xmin += inc
            
    return tile_regions


def region_list_to_ogr(region_list: List[Region], dst_ogr: str, dst_fmt: str = 'ESRI Shapefile'):
    """Convert a region list to an OGR vector."""
    
    driver = ogr.GetDriverByName(dst_fmt)
    if os.path.exists(dst_ogr):
        driver.DeleteDataSource(dst_ogr)

    dst_ds = driver.CreateDataSource(dst_ogr)
    dst_lyr = dst_ds.CreateLayer(dst_ogr, geom_type=ogr.wkbPolygon)
    dst_lyr.CreateField(ogr.FieldDefn('id', ogr.OFTInteger))
        
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
    
    ## Cleanup
    dst_ds = None


def parse_cli_region(region_list: List[str], verbose: bool = True):
    """Parse a list of region strings into Region objects."""
    
    these_regions = []
    
    for i_region in region_list:
        if i_region is None:
            continue
        
        i_region_s = i_region.split(':')
        base_region_str = i_region_s[0]
        args = utils.args2dict(i_region_s[1:])
        
        tmp_region = Region().from_string(base_region_str)
        
        if tmp_region.is_valid(check_xy=True):
            if 'pct_buffer' in args:
                tmp_region.buffer(pct=utils.float_or(args['pct_buffer']))
            these_regions.append(tmp_region)
            
        elif base_region_str == 'tile_set':
            these_regions = generate_tile_set(**args)
            
        elif base_region_str == 'coordinates':
            kwargs = factory.args2dict(i_region_s[1:])
            tmp_region = quarter_tile_from_coordinates(**kwargs)
            these_regions.append(tmp_region)
            
        elif base_region_str == 'q':
            if len(i_region_s) > 1:
                coords = fetch_gps_coordinates(i_region_s[1])
                if coords:
                    tmp_region = quarter_tile_from_coordinates(x=coords[0], y=coords[1])
                    these_regions.append(tmp_region)
                    
        else:
            ## Try parsing as OGR source
            tmp_regions_ogr = ogr_wkts(base_region_str)
            for r in tmp_regions_ogr:
                if r.is_valid():
                    if 'pct_buffer' in args:
                        r.buffer(pct=utils.float_or(args['pct_buffer']))
                    these_regions.append(r)

    if verbose:
        if not these_regions:
            utils.echo_warning_msg(f'Failed to parse region(s): {region_list}')
        else:
            msg = f'Parsed {len(these_regions)} region(s)'
            if len(these_regions) > 4:
                msg += f': {these_regions[:2]}...{these_regions[-2:]}'
            else:
                msg += f': {these_regions}'
            utils.echo_msg(msg)
            
    return these_regions


def fetch_gps_coordinates(q):
    from cudem import fetches
    
    fetches.GPSCoordinates(q)
    gpsc_api_url = "http://www.gps-coordinates.net/api/"
    
    if utils.str_or(q) is not None:
        q_url = f'{gpsc_api_url}{q}'
        _req = fetches.Fetch(q_url, verbose=True).fetch_req()
        if _req is not None:
            try:
                results = _req.json()
                if results.get("responseCode") == '200':
                    x = utils.float_or(results["longitude"])
                    y = utils.float_or(results["latitude"])
                    return x, y
            except Exception:
                pass
    return None


def quarter_tile_from_coordinates(x=None, y=None):
    x = utils.float_or(x)
    y = utils.float_or(y)
    
    if x is not None and y is not None:
        x_min = math.floor(x * 4) / 4
        x_max = math.ceil(x * 4) / 4
        y_min = math.floor(y * 4) / 4
        y_max = math.ceil(y * 4) / 4
        return Region(xmin=x_min, xmax=x_max, ymin=y_min, ymax=y_max)
        
    return None


## ==============================================
## Command-line Interface (CLI)
## $ regions
##
## regions cli
## ==============================================
def regions_cli():
    """Run regions from command-line using argparse."""
    
    parser = argparse.ArgumentParser(
        description=f"%(prog)s ({__version__}): Process and generate regions.",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="CUDEM home page: <http://cudem.colorado.edu>"
    )

    parser.add_argument(
        '-R', '--region', '--aoi',
        required=True,
        action='append',
        help="The desired REGION \n"
             "Where a REGION is xmin/xmax/ymin/ymax[/zmin/zmax[/wmin/wmax/umin/umax]]\n"
             "Use '-' to indicate no bounding range; e.g. -R -/-/-/-/-10/10/1/-/-/- \n"
             "OR an OGR-compatible vector file with regional polygons.\n"
             "Where the REGION is /path/to/vector[:zmin/zmax[/wmin/wmax/umin/umax]].\n"
             "If a vector file is supplied, will use each region found therein.\n"
             "Optionally, append `:pct_buffer=<value>` to buffer the region(s) by a percentage."
    )
    
    parser.add_argument(
        '-J', '--s_srs',
        help="Set the SOURCE projection."
    )
    
    parser.add_argument(
        '-P', '--t_srs',
        help="Set the TARGET projection."
    )
    
    parser.add_argument(
        '-T', '--tile_set',
        type=float,
        help="Generate a TILESET from the input region with the specified increment."
    )
    
    parser.add_argument(
        '-B', '--buffer',
        type=float,
        help="BUFFER the region with a buffer-value."
    )
    
    parser.add_argument(
        '-e', '--echo',
        action='store_true',
        help="ECHO the <processed> region."
    )

    parser.add_argument(
        '-t', '--te',
        action='store_true',
        help="ECHO the region in 'xmin ymin xmax ymax' format."
    )
    
    parser.add_argument(
        '-n', '--name',
        action='store_true',
        help="Print the region as a NAME."
    )
    
    parser.add_argument(
        '-m', '--merge',
        action='store_true',
        help="MERGE all regions into a single region."
    )
    
    parser.add_argument(
        '-q', '--quiet',
        action='store_true',
        help="Lower the verbosity to a quiet."
    )
    
    parser.add_argument(
        '--version',
        action='version',
        version=f'CUDEM {__cudem_version__} :: %(prog)s {__version__}'
    )

    args = parser.parse_args()

    ## Handle verbosity
    if not args.quiet:
        pass 

    if not args.region:
        utils.echo_error_msg('You must specify at least one region')
        sys.exit(-1)

    these_regions = parse_cli_region(args.region, verbose=not args.quiet)

    ## Generate tile set if requested
    if args.tile_set is not None:
        for this_region in these_regions:
            these_tiles = generate_tile_set(this_region.format('gmt'), args.tile_set)
            region_list_to_ogr(these_tiles, 'regions_tile_set.shp')
        ## Clear main regions if we just generated tiles
        these_regions = []

    ## Merge regions if requested
    if args.merge and these_regions:
        merged_region = these_regions[0].copy()
        for r in these_regions[1:]:
            merged_region = regions_merge(merged_region, r)
        these_regions = [merged_region]
        
    for this_region in these_regions:
        if args.s_srs is not None:
            this_region.src_srs = args.s_srs

        if args.t_srs is not None:
            this_region.warp(args.t_srs)

        if args.buffer is not None:
            this_region.buffer(x_bv=args.buffer, y_bv=args.buffer)

        if args.te:
            print(this_region.format('te'))
        elif args.echo:
            print(this_region.format('gmt'))
        elif args.name:
            print(this_region.format('fn'))
        else:
            this_region.export_as_ogr(f'region_{this_region.format("fn")}.geojson')


if __name__ == '__main__':
    regions_cli()
