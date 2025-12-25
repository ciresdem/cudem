### xyzfun.py - CUDEM utilities and functions
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
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

import sys
import json
import copy
from typing import List, Dict, Optional, Generator, Union, Any, TextIO

import numpy as np
from osgeo import gdal, ogr, osr

try:
    from scipy import spatial
except ImportError:
    spatial = None

from cudem import utils
from cudem import regions

try:
    from cudem import mbsfun
except ImportError:
    mbsfun = None

from cudem.utils import float_or, str_or


class XYZPoint:
    """Representing an XYZ data point."""
    
    def __init__(self, x=None, y=None, z=None, w=1, u=0,
                 src_srs='epsg:4326', z_units='m', z_datum='msl'):
        self.x = float_or(x)
        self.y = float_or(y)
        self.z = float_or(z)
        self.w = float_or(w, 1)
        self.u = float_or(u, 0)
        self.src_srs = str_or(src_srs, 'epsg:4326')
        self.z_units = z_units
        self.z_datum = z_datum

        
    def __repr__(self):
        return f'<XYZPoint x: {self.x} y: {self.y} z: {self.z}>'

    
    def __str__(self):
        return f'<XYZPoint x: {self.x} y: {self.y} z: {self.z}>'

    
    def copy(self):
        return XYZPoint(
            x=self.x, y=self.y, z=self.z,
            w=self.w, u=self.u,
            src_srs=self.src_srs,
            z_units=self.z_units,
            z_datum=self.z_datum
        )

    
    def reset(self):
        self.x = self.y = self.z = None
        self.w = 1
        self.u = 0
        self.src_srs = 'epsg:4326'
        self.z_units = 'm'
        self.z_datum = 'msl'

        
    def is_valid(self) -> bool:
        return not (self.x is None or self.y is None or self.z is None)

    
    ## Alias for backward compatibility
    valid_p = is_valid

    def from_list(self, xyz_list: List[Any], x_pos=0, y_pos=1,
                  z_pos=2, w_pos=3, u_pos=4):
        """Load XYZ data from a list."""
        
        if len(xyz_list) > x_pos:
            self.x = float_or(xyz_list[x_pos])
        if len(xyz_list) > y_pos:
            self.y = float_or(xyz_list[y_pos])
        if len(xyz_list) > z_pos:
            self.z = float_or(xyz_list[z_pos])
        if len(xyz_list) > w_pos:
            self.w = float_or(xyz_list[w_pos])
        return self

    def from_string(self, xyz_str: str, x_pos=0, y_pos=1, z_pos=2, w_pos=3,
                    u_pos=4, delim=" "):
        """Load XYZ data from a delimited string."""
        
        this_line = xyz_str.strip()
        parts = this_line.split(delim) if delim else this_line.split()
        return self.from_list(parts, x_pos, y_pos, z_pos, w_pos)

    def export_as_list(self, include_z=False, include_w=False, include_u=False) -> List[float]:
        """Export XYZ as a list."""
        
        xyz_list = [self.x, self.y]
        if include_z:
            xyz_list.append(self.z)
        if include_w:
            xyz_list.append(self.w)
        if include_u:
            xyz_list.append(self.u)
        return xyz_list

    def export_as_string(self, delim: str, include_z=False, include_w=False,
                         include_u=False, precision=6) -> str:
        """Export XYZ data as a string."""
        
        data_list = self.export_as_list(include_z, include_w, include_u)
        
        ## Helper to format based on index (coordinates vs values)
        def fmt(idx, val):
            if val is None:
                return "nan"
            ## Assuming first two are X/Y, use higher precision if needed, else user precision
            p = precision if idx > 1 else 8
            return f'{val:.{p}f}'

        return '{}\n'.format(
            delim.join([fmt(j, x) for j, x in enumerate(data_list)])
        )

    
    def export_as_wkt(self, include_z=False) -> str:
        if include_z and self.z is not None:
            return f'POINT ({self.x} {self.y} {self.z})'
        return f'POINT ({self.x} {self.y})'

    
    def dump(self, delim=' ', include_z=True, include_w=False,
             include_u=False, encode=False, dst_port=sys.stdout, precision=6):
        """Dump XYZ as a string to dst_port."""
        
        line = self.export_as_string(
            delim, include_z, include_w, include_u, precision
        )
        dst_port.write(line.encode('utf-8') if encode else line)

        
    def transform(self, dst_trans):
        """Transform the x/y/z using the dst_trans OSR transformation."""
        
        try:
            point = ogr.CreateGeometryFromWkt(self.export_as_wkt(include_z=True))
            point.Transform(dst_trans)
            
            # Check for infinity which indicates projection failure in some cases
            wkt = point.ExportToWkt()
            if 'inf' not in wkt.lower():
                self.x = point.GetX() 
                self.y = point.GetY()
                self.z = point.GetZ()
        except Exception as e:
            sys.stderr.write(f'Transform error: {e}\n')
        return self

    
    def warp(self, dst_srs=None):
        """Transform the x/y using dst_srs string (e.g., 'epsg:3857')."""
        
        dst_srs = str_or(dst_srs)
        if dst_srs is None or self.src_srs is None:
            return self

        src_ref = osr.SpatialReference()
        src_ref.SetFromUserInput(self.src_srs)
        
        dst_ref = osr.SpatialReference()
        dst_ref.SetFromUserInput(dst_srs)

        # Attempt to set traditional axis mapping order (Long/Lat vs Lat/Long)
        try:
            src_ref.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
            dst_ref.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        except AttributeError:
            pass  # OAMS_TRADITIONAL_GIS_ORDER might not be available in older GDAL
        
        dst_trans = osr.CoordinateTransformation(src_ref, dst_ref)
        if dst_trans is not None:
            self.transform(dst_trans)
        
        return self


###############################################################################    
## XYZ Processing Utilities - DEPRECIATED
###############################################################################
_xyz_config = {
    'delim': None, 'xpos': 0, 'ypos': 1, 'zpos': 2, 'wpos': 3,
    'upos': 4, 'skip': 0, 'z-scale': 1, 'x-off': 0, 'y-off': 0,
    'name': '<xyz-data-stream>', 'upper_limit': None,
    'lower_limit': None, 'epsg': 4326, 'warp': None,
    'verbose': False,
}

_known_delims = [None, ',', '/', ':']
_known_xyz_fmts = ['xyz', 'csv', 'dat', 'ascii']

def xyz_line_delim(xyz_line: str) -> Optional[str]:
    """Guess a line delimiter."""
    
    for delim in _known_delims:
        if delim is None:
            ## None implies whitespace splitting in split() default
            if len(xyz_line.split()) > 1:
                return None
        else:
            if len(xyz_line.split(delim)) > 1:
                return delim
    return None


def xyz_warp(xyz: List[float], dst_trans: osr.CoordinateTransformation) -> List[float]:
    """Transform the x/y using the dst_trans."""
    
    if dst_trans is None:
        return xyz
    
    ## Create basic point geometry
    ## Note: OGR/GDAL python bindings usually expect Point(x,y)
    point = ogr.CreateGeometryFromWkt(f'POINT ({xyz[0]} {xyz[1]})')
    point.Transform(dst_trans)
    
    return [point.GetX(), point.GetY(), xyz[2]]


def xyz_parse_line(xyz_line: str, xyz_c: Dict = None) -> Optional[List[float]]:
    """Parse an XYZ line-string using configuration dict."""
    
    if xyz_c is None:
        xyz_c = _xyz_config.copy()

    try:
        this_line = xyz_line.strip()
    except AttributeError as e:
        utils.echo_error_msg(f'Input is not a string: {e}')
        return None
    except Exception as e:
        utils.echo_error_msg(e)
        return None
    
    ## Auto-detect delimiter if not set
    if xyz_c.get('delim') is None:
        xyz_c['delim'] = xyz_line_delim(this_line)
        
    parts = this_line.split(xyz_c['delim'])
    
    try:
        x = float(parts[xyz_c['xpos']]) + float(xyz_c.get('x-off', 0))
        y = float(parts[xyz_c['ypos']]) + float(xyz_c.get('y-off', 0))
        z = float(parts[xyz_c['zpos']]) * float(xyz_c.get('z-scale', 1))
        return [x, y, z]
    except (IndexError, ValueError) as e:
        if xyz_c.get('verbose'):
            utils.echo_error_msg(f'{e}, Line: {parts}')
        return None
    except Exception as e:
        if xyz_c.get('verbose'):
            utils.echo_error_msg(e)
        return None


def xyz_parse(src_xyz, xyz_c: Dict = None, region=None, verbose=False) -> Generator[List[float], None, None]:
    """XYZ file parsing generator.

    Args:
      src_xyz: Iterable (list/generator/file) of XYZ string lines.
      xyz_c (dict): XYZ config dictionary.
      region (list): A region list [xmin, xmax, ymin, ymax].
      verbose (bool): Increase verbosity.

    Yields:
      list: XYZ data [x, y, z].
    """
    
    if xyz_c is None:
        xyz_c = _xyz_config.copy()

    ln = 0
    skip = int(xyz_c.get('skip', 0))
    
    ## Handle Coordinate Transformation Setup
    warp_epsg = xyz_c.get('warp')
    src_epsg = xyz_c.get('epsg')
    
    dst_trans = None
    if warp_epsg and warp_epsg != src_epsg:
        src_srs = osr.SpatialReference()
        src_srs.ImportFromEPSG(int(src_epsg))
        
        dst_srs = osr.SpatialReference()
        dst_srs.ImportFromEPSG(int(warp_epsg))
    
        try:
            src_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
            dst_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        except AttributeError:
            pass
        
        dst_trans = osr.CoordinateTransformation(src_srs, dst_srs)
    else:
        xyz_c['warp'] = None

    for xyz_line_str in src_xyz:
        if skip > 0:
            skip -= 1
            continue

        this_xyz = xyz_parse_line(xyz_line_str, xyz_c)
        
        if this_xyz is None:
            continue
            
        ## Transform
        if dst_trans is not None:
            this_xyz = xyz_warp(this_xyz, dst_trans)
            
        ## Filter: Spatial Region
        if region is not None:
            if region.is_valid():
                if not xyz_in_region_p(this_xyz, region):
                    continue
        
        ## Filter: Z-Range
        upper = xyz_c.get('upper_limit')
        lower = xyz_c.get('lower_limit')
        if upper is not None or lower is not None:
            if not regions.z_region_pass([0, 0, 0, 0, this_xyz[2], this_xyz[2]], 
                                         upper_limit=upper, lower_limit=lower):
                continue

        ln += 1
        yield this_xyz
        
    if verbose:
        utils.echo_msg(f"Parsed {ln} data records from {xyz_c.get('name', 'source')}")


def xyz_dump(src_xyz, xyz_c: Dict = None, region=None,
             verbose=False, dst_port: TextIO = sys.stdout):
    """Dump the XYZ data from the generator to dst_port."""
    
    for xyz in xyz_parse(src_xyz, xyz_c, region, verbose):
        xyz_line(xyz, dst_port, encode=False)


def xyz2py(src_xyz) -> List[List[float]]:
    """Return src_xyz as a Python list."""
    
    return list(xyz_parse(src_xyz))


def xyz_block(src_xyz, region: List[float], inc: float, weights=False, verbose=False):
    """Block the src_xyz data to the mean block value.

    Args:
      src_xyz: Generator of xyz data.
      region: [xmin, xmax, ymin, ymax].
      inc: Blocking increment.
      weights (bool): Block using weights (assumes 4th column is weight).
      verbose (bool): Verbosity.

    Yields:
      list: [geo_x, geo_y, z] for each block.
    """
    
    xcount, ycount, dst_gt = regions.Region().from_list(region).geo_transform(x_inc=inc)
    
    sum_array = np.zeros((ycount, xcount), dtype=np.float64)
    pt_array = np.zeros((ycount, xcount), dtype=np.int32)
    wt_array = np.zeros((ycount, xcount), dtype=np.float64) if weights else None

    if verbose:
        utils.echo_msg(f'Blocking data to {ycount}x{xcount} grid')

    for this_xyz in src_xyz:
        x, y, z = this_xyz[0], this_xyz[1], this_xyz[2]
        
        if region[0] < x < region[1] and region[2] < y < region[3]:
            xpos, ypos = utils._geo2pixel(x, y, dst_gt)
            
            try:
                if 0 <= xpos < xcount and 0 <= ypos < ycount:
                    weight = this_xyz[3] if weights and len(this_xyz) > 3 else 1
                    
                    if weights:
                        sum_array[ypos, xpos] += (z * weight)
                        wt_array[ypos, xpos] += weight
                    else:
                        sum_array[ypos, xpos] += z
                    
                    pt_array[ypos, xpos] += 1
            except IndexError:
                pass

    ## Avoid division by zero
    with np.errstate(divide='ignore', invalid='ignore'):
        if weights:
            ## Weighted Mean = Sum(z * w) / Sum(w)
            out_array = sum_array / wt_array
        else:
            ## Mean = Sum(z) / Count
            out_array = sum_array / pt_array

    ## Cleanup
    del sum_array, pt_array, wt_array

    ## Iterate and yield valid blocks
    for y in range(ycount):
        for x in range(xcount):
            z_val = out_array[y, x]
            if not np.isnan(z_val):
                geo_x, geo_y = utils._pixel2geo(x, y, dst_gt)
                yield [geo_x, geo_y, z_val]


def xyz_block_t(src_xyz, src_region: List[float], inc: float, verbose=False):
    """Block the src_xyz data and return indices of points in blocks.

    Returns:
        blkArray: 2D array of lists, where each list contains indices of points in that block.
        xyzArray: List of all valid points found in region.
    """
    
    xcount, ycount, dst_gt = regions.Region().from_list(src_region).geo_transform(x_inc=inc)
    
    ## Create object array to hold lists
    blk_array = np.empty((ycount, xcount), dtype=object)
    for y in range(ycount):
        for x in range(xcount):
            blk_array[y, x] = []
            
    xyz_array = []
    
    if verbose:
        utils.echo_msg(f'Blocking data to {ycount}x{xcount} grid')
        
    it = 0
    for this_xyz in src_xyz:
        x, y = this_xyz[0], this_xyz[1]
        
        if src_region[0] < x < src_region[1] and src_region[2] < y < src_region[3]:
            xpos, ypos = utils._geo2pixel(x, y, dst_gt)
            
            if 0 <= xpos < xcount and 0 <= ypos < ycount:
                xyz_array.append(this_xyz)
                blk_array[ypos, xpos].append(it)
                it += 1
                    
    return blk_array, xyz_array


def xyz_line(xyz_data: List[Any], dst_port: TextIO = sys.stdout, encode=False, delim=None):
    """Write an XYZ line to dst_port."""
    
    if delim is None:
        delim = _xyz_config['delim'] if _xyz_config['delim'] is not None else ' '
    
    line = '{}\n'.format(delim.join([str(x) for x in xyz_data]))
    if encode:
        line = line.encode('utf-8')
    dst_port.write(line)


def xyz2wkt(xyz: List[float]) -> str:
    """Convert XYZ list to WKT Point string."""
    
    return f'POINT ({xyz[0]} {xyz[1]})'


def xyz_in_region_p(xyz: List[float], region: Union[List[float], str]) -> bool:
    """Check if XYZ point is inside the given region."""
    
    if isinstance(region, list) and len(region) == 4:
        ## Fast Bounding Box check
        x, y = xyz[0], xyz[1]
        if x < region[0] or x > region[1] or y < region[2] or y > region[3]:
            return False
        return True
    
    ## Fallback to OGR geometry check for complex regions (WKT)
    if regions.region_wkt_p(region):
        try:
            p_geom = ogr.CreateGeometryFromWkt(xyz2wkt(xyz))
            r_geom = ogr.CreateGeometryFromWkt(region)
            return p_geom.Within(r_geom)
        except Exception:
            return False
            
    return False


def xyz_inf(src_xyz):
    """Generate metadata/info for an XYZ source.
    
    Calculates min/max bounds and convex hull if scipy is available.
    """
    
    pts = []
    xyzi = {
        'name': getattr(src_xyz, 'name', 'unknown'),
        'numpts': 0,
        'minmax': None
    }
    
    for i, l in enumerate(xyz_parse(src_xyz)):
        if xyzi['minmax'] is None:
            ## Initialize: [xmin, xmax, ymin, ymax, zmin, zmax]
            xyzi['minmax'] = [l[0], l[0], l[1], l[1], l[2], l[2]]
        else:
            ## Update bounds
            xyzi['minmax'][0] = min(xyzi['minmax'][0], l[0])
            xyzi['minmax'][1] = max(xyzi['minmax'][1], l[0])
            xyzi['minmax'][2] = min(xyzi['minmax'][2], l[1])
            xyzi['minmax'][3] = max(xyzi['minmax'][3], l[1])
            xyzi['minmax'][4] = min(xyzi['minmax'][4], l[2])
            xyzi['minmax'][5] = max(xyzi['minmax'][5], l[2])
            
        pts.append(l)
        xyzi['numpts'] = i + 1

    if xyzi['numpts'] > 0 and spatial is not None:
        try:
            ## Calculate 2D Convex Hull
            pts_2d = [[p[0], p[1]] for p in pts]
            hull = spatial.ConvexHull(pts_2d, qhull_options='Qt')
            out_hull = [pts_2d[i] for i in hull.vertices]
            ## Close the polygon
            out_hull.append(out_hull[0])
            
            xyzi['wkt'] = regions.create_wkt_polygon(out_hull, xpos=0, ypos=1)
        except Exception:
             ## Fallback to BBox WKT
            xyzi['wkt'] = regions.Region().from_list(xyzi['minmax']).export_as_wkt()
    elif xyzi['minmax']:
        xyzi['wkt'] = regions.Region().from_list(xyzi['minmax']).export_as_wkt()

    ## Save .inf file
    if xyzi['numpts'] > 0:
        try:
            with open(f"{xyzi['name']}.inf", 'w') as inf:
                inf.write(json.dumps(xyzi))
        except Exception:
            pass

    return xyzi


def xyz_inf_entry(entry):
    """Find the region of the xyz datalist entry."""
    path = entry[0]
    with open(path) as infile:
        if mbsfun:
            try:
                minmax = mbsfun.mb_inf(infile)
                return minmax
            except Exception:
                pass
        
        ## Fallback to manual xyz parsing
        info = xyz_inf(infile)
        return info['minmax']


def xyz_yield_entry(entry: List[Any], region=None, verbose=False, z_region=None, epsg=None):
    """Yield XYZ data from a datalist entry.
    
    Args:
      entry: [path, fmt, weight, ...]
    """
    
    xyzc = copy.deepcopy(_xyz_config)
    xyzc['name'] = entry[0]
    
    if z_region is not None and len(z_region) >= 2:
        xyzc['lower_limit'] = z_region[0]
        xyzc['upper_limit'] = z_region[1]
    
    ## If EPSG is specified for output, or input needs specific parsing
    if epsg:
        xyzc['warp'] = epsg

    weight = entry[2] if len(entry) > 2 else None

    with open(entry[0], 'r') as infile:
        for line_data in xyz_parse(infile, xyz_c=xyzc, region=region, verbose=verbose):
            if weight is not None:
                line_data.append(weight)
            yield line_data


def xyz_dump_entry(entry: List[Any], dst_port=sys.stdout, region=None,
                   verbose=False, z_region=None, epsg=None):
    """Dump XYZ data from a datalist entry to port."""
    
    for xyz in xyz_yield_entry(entry, region, verbose, z_region, epsg):
        xyz_line(xyz, dst_port, encode=True)


def xyz_chunks():
    """Placeholder for future implementation."""
    
    pass


### End
