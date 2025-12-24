### srsfun.py - Projection functions
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
##
## srsfun.py is part of CUDEM
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
## Functions, etc. for common pyproj/osr usage.
##
### Code:

import os
from typing import Optional, Tuple, Any, Dict, Union

import numpy as np
import pyproj
from pyproj import CRS
from osgeo import osr, ogr

from cudem import utils, vdatums, regions 

## Initialize configuration and GDAL settings
gc = utils.config_check()
ogr.DontUseExceptions()
osr.DontUseExceptions()


###############################################################################
## OSR/WKT/proj
###############################################################################

def split_srs(srs: str, as_epsg: bool = False) -> Tuple[Optional[Union[str, int]], Optional[Union[str, int]]]:
    """Split an SRS into horizontal and vertical elements.

    Args:
        srs (str): An srs string for input into SetFromUserInput.
        as_epsg (bool): Try to output EPSG codes if possible.

    Returns:
        tuple: (horz_srs, vert_srs) as WKT or EPSG.
    """
    
    if srs is None:
        return None, None

    ## Strip geoid definition if present
    if '+geoid' in srs:
        srs = '+'.join([x for x in srs.split('+') if 'geoid' not in x])

    # Handle ESRI style strings
    vert_epsg = None
    vert_wkt = None
    if srs.startswith('ESRI'):
        esri_split = srs.split('+')
        if len(esri_split) > 1:
            vert_epsg = esri_split[1]
            vert_srs_obj = osr.SpatialReference()
            vert_srs_obj.SetFromUserInput(f'EPSG:{vert_epsg}')
            vert_wkt = vert_srs_obj.ExportToWkt()
        srs = esri_split[0]

    ## Process standard SRS
    src_srs = osr.SpatialReference()
    if src_srs.SetFromUserInput(srs) != 0:
        return None, None

    try:
        srs_wkt = src_srs.ExportToWkt()
        wkt_crs = CRS.from_wkt(srs_wkt)
    except Exception:
        return None, None

    if wkt_crs.is_compound:
        horz = wkt_crs.sub_crs_list[0]
        vert = wkt_crs.sub_crs_list[1]
        
        horz_epsg = horz.to_epsg()
        horz_wkt = horz.to_wkt()
        
        vert_epsg = vert.to_epsg()
        vert_wkt = vert.to_wkt()
    else:
        horz = wkt_crs
        horz_epsg = horz.to_epsg()
        horz_wkt = horz.to_wkt()
        # Vertical remains None/from ESRI parsing

    if as_epsg:
        return (horz_epsg if horz_epsg is not None else horz_wkt,
                vert_epsg if vert_epsg is not None else vert_wkt)
    else:
        return horz_wkt, vert_epsg


def combine_epsgs(src_horz: str, src_vert: str, name: str = 'Combined') -> Optional[str]:
    """Combine src_horz and src_vert into a CompoundCS WKT."""
    
    if src_horz is None or src_vert is None:
        return None

    try:
        horz_srs = osr.SpatialReference()
        horz_srs.SetFromUserInput(str(src_horz))
        
        vert_srs = osr.SpatialReference()
        vert_srs.SetFromUserInput(f'epsg:{src_vert}')
        
        src_srs = osr.SpatialReference()
        src_srs.SetCompoundCS(f'{name}', horz_srs, vert_srs)
        return src_srs.ExportToWkt()
    except Exception:
        return None


def wkt2geom(wkt: str) -> ogr.Geometry:
    """Transform a WKT string to an OGR geometry."""
    
    return ogr.CreateGeometryFromWkt(wkt)


def osr_srs(src_srs: str) -> Optional[osr.SpatialReference]:
    """Create an osr.SpatialReference from a string."""
    
    try:
        srs = osr.SpatialReference()
        if srs.SetFromUserInput(src_srs) == 0:
            return srs
        return None
    except Exception:
        return None


def osr_wkt(src_srs: str, esri: bool = False) -> Optional[str]:
    """Convert a src_srs to WKT string."""
    
    try:
        sr = osr.SpatialReference()
        if sr.SetFromUserInput(src_srs) != 0:
            return None
        
        if esri:
            sr.MorphToESRI()
        return sr.ExportToWkt()
    except Exception:
        return None


def osr_prj_file(dst_fn: str, src_srs: str) -> int:
    """Generate a .prj file given a src_srs.
    
    Returns:
        0 on success, -1 on failure.
    """
    
    wkt = osr_wkt(src_srs, esri=True)
    if not wkt:
        return -1
        
    try:
        with open(dst_fn, 'w') as out:
            out.write(wkt)
    except IOError:
        return -1

    return 0 if os.path.exists(dst_fn) else -1


def srs_get_cstype(in_srs: str) -> str:
    """Determine if SRS is Geographic or Projected."""
    
    src_srs = osr.SpatialReference()
    src_srs.SetFromUserInput(in_srs)
    
    if src_srs.IsGeographic() == 1:
        return 'GEOGCS'
    return 'PROJCS'


def epsg_from_input(in_srs: str) -> Tuple[Optional[str], Optional[str]]:
    """Get the authority codes (EPSGs) from SRS string.

    Returns:
        list: [horz_code, vert_code]
    """
    
    src_vert = None
    
    ## Clean Geoid
    if '+geoid' in in_srs:
        in_srs = '+'.join([x for x in in_srs.split('+') if 'geoid' not in x])

    ## Handle ESRI
    if in_srs.startswith('ESRI'):
        esri_split = in_srs.split('+')
        if len(esri_split) > 1:
            vert_epsg = esri_split[1]
            # Verify validity
            vert_srs = osr.SpatialReference()
            vert_srs.SetFromUserInput(f'EPSG:{vert_epsg}')
            src_vert = vert_epsg
        in_srs = esri_split[0]

    src_srs = osr.SpatialReference()
    src_srs.SetFromUserInput(in_srs)

    ## Horizontal
    cstype = 'GEOGCS' if src_srs.IsGeographic() else 'PROJCS'
    src_srs.AutoIdentifyEPSG()
    
    an = src_srs.GetAuthorityName(cstype)
    ac = src_srs.GetAuthorityCode(cstype)

    if ac is None:
        src_horz = src_srs.ExportToProj4()
        if not src_horz:
            src_horz = src_srs.ExportToWkt()
    else:
        src_horz = f'{an}:{ac}'

    ## Vertical
    if src_srs.IsVertical() == 1 and src_vert is None:
        csvtype = 'VERT_CS'
        src_vert = src_srs.GetAuthorityCode(csvtype)

    return src_horz, src_vert


def osr_parse_srs(src_srs: osr.SpatialReference, return_vertcs: bool = True) -> Optional[str]:
    """Parse an OSR SRS object and return a proj/auth string."""
    
    if src_srs is None:
        return None

    if src_srs.IsLocal() == 1:
        return src_srs.ExportToWkt()

    cstype = 'GEOGCS' if src_srs.IsGeographic() else 'PROJCS'
    src_srs.AutoIdentifyEPSG()
    
    an = src_srs.GetAuthorityName(cstype)
    ac = src_srs.GetAuthorityCode(cstype)

    vn = vc = None
    if src_srs.IsVertical() == 1:
        csvtype = 'VERT_CS'
        vn = src_srs.GetAuthorityName(csvtype)
        vc = src_srs.GetAuthorityCode(csvtype)

    if an is not None and ac is not None:
        if vn is not None and vc is not None:
            return f'{an}:{ac}+{vc}'
        else:
            return f'{an}:{ac}'
    else:
        dst_srs = src_srs.ExportToProj4()
        return dst_srs if dst_srs else None


def parse_srs(src_srs: str = None, dst_srs: str = None) -> Dict[str, Any]:
    """Parse source and destination SRS to determine transform requirements."""
    
    transform = {}
    
    if src_srs is None or dst_srs is None:
        return transform

    want_vertical = True
    src_geoid = None
    dst_geoid = 'g2018'  # default
    in_vertical_epsg = None
    out_vertical_epsg = None
    in_vertical_epsg_esri = None

    ## Parse Source Geoid
    tmp_src_srs = src_srs.split('+geoid:')
    src_srs = tmp_src_srs[0]
    if len(tmp_src_srs) > 1:
        src_geoid = tmp_src_srs[1]

    ## Parse Destination Geoid
    tmp_dst_srs = dst_srs.split('+geoid:')
    dst_srs = tmp_dst_srs[0]
    if len(tmp_dst_srs) > 1:
        dst_geoid = tmp_dst_srs[1]

    ## Check for ESRI prefix
    is_esri = False
    if 'ESRI' in src_srs.upper():
        is_esri = True
        srs_split = src_srs.split('+')
        src_srs = srs_split[0]
        if len(srs_split) > 1:
            in_vertical_epsg_esri = srs_split[1]

    ## Check for Tidal frames
    last_param = src_srs.split('+')[-1]
    last_int = utils.int_or(last_param)
    if last_int in vdatums._tidal_frames.keys():
        base_srs = src_srs.split('+')[0]
        tidal_epsg = vdatums._tidal_frames[last_int]['epsg']
        src_srs = f'{base_srs}+{tidal_epsg}'

    ## Create PyProj CRS objects
    try:
        in_crs = pyproj.CRS.from_user_input(src_srs)
        out_crs = pyproj.CRS.from_user_input(dst_srs)
    except Exception as e:
        utils.echo_error_msg(f"Error parsing CRS: {e}")
        return transform

    ## Handle Input CRS Component Splitting
    if in_crs.is_compound:
        in_horizontal_crs = in_crs.sub_crs_list[0]
        in_vertical_crs = in_crs.sub_crs_list[1]
        in_vertical_epsg = in_vertical_crs.to_epsg() or in_vertical_crs.name
    else:
        in_horizontal_crs = in_crs
        in_vertical_crs = None
        want_vertical = False

    ## Handle Output CRS Component Splitting
    if out_crs.is_compound:
        out_horizontal_crs = out_crs.sub_crs_list[0]
        out_vertical_crs = out_crs.sub_crs_list[1]
        out_vertical_epsg = out_vertical_crs.to_epsg()
    else:
        out_horizontal_crs = out_crs
        out_vertical_crs = None
        want_vertical = False
        out_vertical_epsg = None

    ## Override with ESRI vertical if present
    if is_esri and in_vertical_epsg_esri is not None:
        in_vertical_epsg = in_vertical_epsg_esri
        if out_vertical_epsg is not None:
            want_vertical = True

    ## Disable vertical transform if EPSGs match and no custom source geoid is set
    if want_vertical:
        if (in_vertical_epsg == out_vertical_epsg) and src_geoid is None:
            want_vertical = False

    ## Populate transform dictionary
    transform.update({
        'src_horz_crs': in_horizontal_crs,
        'dst_horz_crs': out_horizontal_crs,
        'src_vert_crs': in_vertical_crs,
        'dst_vert_crs': out_vertical_crs,
        'src_vert_epsg': in_vertical_epsg,
        'dst_vert_epsg': out_vertical_epsg,
        'src_geoid': src_geoid,
        'dst_geoid': dst_geoid,
        'want_vertical': want_vertical
    })

    return transform


def set_vertical_transform(transform: Dict[str, Any], region=None, infos=None,
                           cache_dir: str = './', verbose: bool = True) -> Dict[str, Any]:
    """Generate the vertical transformation grid and pipeline."""
    
    ## Determine region for transformation
    if region is None and infos is not None:
        vd_region = regions.Region().from_list(infos.minmax)
        vd_region.src_srs = transform['src_horz_crs'].to_proj4()
    elif region is not None:
        vd_region = region.copy()
        vd_region.src_srs = transform['dst_horz_crs'].to_proj4()
    else:
        utils.echo_error_msg("Cannot set vertical transform: No region provided.")
        return transform

    vd_region.zmin = None
    vd_region.zmax = None
    vd_region.warp('epsg:4326')
    vd_region.buffer(pct=10)

    if not vd_region.valid_p():
        utils.echo_warning_msg('Failed to generate transformation: Invalid region')
        return transform

    ## Define transformation grid filename
    trans_fn = os.path.join(
        cache_dir, 
        f"_vdatum_trans_{transform['src_vert_epsg']}_{transform['dst_vert_epsg']}_{vd_region.format('fn')}.tif"
    )
    transform['trans_fn'] = trans_fn

    ## Generate grid if it doesn't exist
    if not os.path.exists(trans_fn):
        with utils.ccp(
            desc=f"Generating v-grid {os.path.basename(trans_fn)}",
            leave=verbose,
            disable=not verbose
        ) as _:
            ## Initial grid resolution (3 arc-seconds)
            vd_x_inc = vd_y_inc = utils.str2inc('3s')
            xcount, ycount, _ = vd_region.geo_transform(
                x_inc=vd_x_inc, y_inc=vd_y_inc, node='grid'
            )

            ## Ensure grid isn't too small (at least 10x10)
            while xcount <= 10 or ycount <= 10:
                vd_x_inc /= 2
                vd_y_inc /= 2
                xcount, ycount, _ = vd_region.geo_transform(
                    x_inc=vd_x_inc, y_inc=vd_y_inc, node='grid'
                )

            try:
                transform['trans_fn'], transform['trans_fn_unc'] = vdatums.VerticalTransform(
                    'IDW',
                    vd_region,
                    vd_x_inc,
                    vd_y_inc,
                    transform['src_vert_epsg'],
                    transform['dst_vert_epsg'],
                    geoid_in=transform['src_geoid'],
                    geoid_out=transform['dst_geoid'],
                    cache_dir=cache_dir,
                    verbose=False
                ).run(outfile=trans_fn)
            except Exception as e:
                utils.echo_error_msg(f"VerticalTransform Run failed: {e}")
                return transform

    ## Setup PyProj pipelines
    if transform['trans_fn'] is not None and os.path.exists(transform['trans_fn']):
        transform['pipeline'] = (
            f"+proj=pipeline +step {transform['src_horz_crs'].to_proj4()} +inv "
            f"+step +proj=vgridshift +grids={os.path.abspath(transform['trans_fn'])} +inv "
            f"+step {transform['dst_horz_crs'].to_proj4()}"
        )
        
        # Standalone vertical transformer
        transform['vert_transformer'] = pyproj.Transformer.from_pipeline(
            f"+proj=pipeline +step +proj=vgridshift +grids={os.path.abspath(transform['trans_fn'])} +inv"
        )
    else:
        utils.echo_error_msg(
            f"Failed to generate vertical transformation grid for region: {vd_region.format('str')}"
        )

    return transform


def set_transform(src_srs: str = None, dst_srs: str = None, 
                  region=None, infos=None, cache_dir: str = './') -> Dict[str, Any]:
    """Set the pyproj horizontal and vertical transformations."""
    
    transform = parse_srs(src_srs=src_srs, dst_srs=dst_srs)
    
    ## Ensure CRSs are valid
    if transform.get('src_horz_crs') is None or transform.get('dst_horz_crs') is None:
        return transform

    ## Horizontal Pipeline
    transform['horz_pipeline'] = (
        f"+proj=pipeline +step {transform['src_horz_crs'].to_proj4()} +inv "
        f"+step {transform['dst_horz_crs'].to_proj4()}"
    )

    ## Determine Transformation Region
    if region is not None:
        transform['trans_region'] = region.copy()
        transform['trans_region'].src_srs = transform['dst_horz_crs'].to_proj4()
        transform['trans_region'].warp(transform['src_horz_crs'].to_proj4())
    elif infos is not None:
        transform['trans_region'] = regions.Region().from_list(infos.minmax)
        transform['trans_region'].src_srs = infos.src_srs
        transform['trans_region'].warp(transform['dst_horz_crs'].to_proj4())
    else:
        transform['trans_region'] = None

    ## Vertical Transformation
    if transform['want_vertical']:
        transform = set_vertical_transform(transform, region=region, infos=infos, cache_dir=cache_dir)
    else:
        transform['pipeline'] = transform['horz_pipeline']

    ## Create Main Transformer
    try:
        transform['transformer'] = pyproj.Transformer.from_crs(
            transform['src_horz_crs'], 
            transform['dst_horz_crs'], 
            always_xy=True
        )
    except Exception as e:
        utils.echo_warning_msg(
            f"Could not set transformation from {transform['src_horz_crs'].name} "
            f"to {transform['dst_horz_crs'].name}: {e}"
        )
        return transform

    ## Define Data Region (Final bounds check logic)
    if region is not None and region.valid_p():
        base_region = transform['trans_region'].copy() if transform['trans_region'] else region.copy()
        
        if infos is not None:
            inf_region = regions.Region().from_list(infos.minmax)
            data_region = regions.regions_reduce(base_region, inf_region)                
            data_region.src_srs = infos.src_srs
        else:
            data_region = base_region

        if not data_region.valid_p():
             # Fallback to the passed region if the reduced region is invalid
             data_region = region.copy()
    elif infos is not None:
        data_region = regions.Region().from_list(infos.minmax)
        data_region.src_srs = infos.src_srs

    return transform


### End
