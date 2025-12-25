### gdalfun.py - OSGEO functions
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
##
## gdalfun.py is part of CUDEM
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

import os
import sys
import shutil
import math
from typing import Optional, List, Tuple, Dict, Union, Any

import numpy as np
from osgeo import gdal, osr, ogr

# Optional dependencies
try:
    import scipy.ndimage
    import scipy.interpolate
except ImportError:
    pass

try:
    from pyproj import CRS
except ImportError:
    CRS = None

from cudem import utils
from cudem import regions
from cudem import xyzfun

# Initialize Configuration
gc = utils.config_check()
gdal.DontUseExceptions()
ogr.DontUseExceptions()
osr.DontUseExceptions()

# Suppress GDAL logging based on platform
log_dev = 'NUL' if gc['platform'] == 'win32' else '/dev/null'
gdal.SetConfigOption('CPL_LOG', log_dev)


###############################################################################
## OSR/WKT/proj
###############################################################################
def split_srs(srs: str, as_epsg: bool = False) -> Tuple[Any, Any]:
    """Split an SRS into horizontal and vertical elements.

    Args:
        srs (str): An SRS string (WKT, EPSG, Proj4).
        as_epsg (bool): Try to output EPSG codes if possible.

    Returns:
        tuple: (horz_srs, vert_srs)
    """
    
    vert_epsg = None
    vert_wkt = None
    
    if srs is None:
        return None, None

    ## Handle basic compound strings manually if possible
    if np.any(['geoid' in x for x in srs.split('+')]):
        srs = '+'.join(srs.split('+')[:-1])
    
    ## Handle ESRI style strings
    if srs.startswith('ESRI'):
        esri_split = srs.split('+')
        if len(esri_split) > 1:
            vert_epsg = srs.split('+')[1]
            vert_srs_obj = osr.SpatialReference()
            vert_srs_obj.SetFromUserInput(f'EPSG:{vert_epsg}')
            vert_wkt = vert_srs_obj.ExportToWkt()
            
        srs = esri_split[0]
        
    ## Process using PyProj CRS if available for robust compound parsing
    src_srs = osr.SpatialReference()
    src_srs.SetFromUserInput(srs)
    srs_wkt = src_srs.ExportToWkt()

    try:
        wkt_crs = CRS.from_wkt(srs_wkt)
    except Exception:
        ## Fallback if pyproj fails or not installed
        return None, None
    
    if wkt_crs.is_compound:
        horz = wkt_crs.sub_crs_list[0]
        horz_epsg = horz.to_epsg()
        horz_wkt = horz.to_wkt()
        
        vert = wkt_crs.sub_crs_list[1]
        vert_epsg = vert.to_epsg()
        vert_wkt = vert.to_wkt()
    else:
        horz = wkt_crs
        horz_epsg = horz.to_epsg()
        horz_wkt = horz.to_wkt()

    if as_epsg:
        return (
            horz_epsg if horz_epsg is not None else horz_wkt,
            vert_epsg if vert_epsg is not None else vert_wkt
        )
    else:
        return horz_wkt, vert_epsg


def combine_epsgs(src_horz, src_vert, name='Combined'):
    """Combine src_horz and src_vert into a CompoundCS."""
    
    if src_horz is None or src_vert is None:
        return None
    
    horz_srs = osr.SpatialReference()
    horz_srs.SetFromUserInput(str(src_horz))
    
    vert_srs = osr.SpatialReference()
    ## Check if src_vert acts like an EPSG int or WKT string
    if utils.int_or(src_vert) is not None:
        vert_srs.SetFromUserInput(f'epsg:{src_vert}')
    else:
        vert_srs.SetFromUserInput(f'{src_vert}')
        
    src_srs = osr.SpatialReference()
    src_srs.SetCompoundCS(f'{name}', horz_srs, vert_srs)
    
    return src_srs.ExportToWkt()


def wkt2geom(wkt: str) -> ogr.Geometry:
    """Transform a WKT string to an OGR geometry."""
    
    return ogr.CreateGeometryFromWkt(wkt)


def osr_wkt(src_srs, esri: bool = False):
    """Convert a src_srs to WKT."""
    
    try:
        sr = osr.SpatialReference()
        sr.SetFromUserInput(str(src_srs))
        if esri:
            sr.MorphToESRI()
        return sr.ExportToWkt()
    except Exception:
        return None


def osr_prj_file(dst_fn: str, src_srs) -> int:
    """Generate a .prj file given a src_srs."""
    
    try:
        with open(dst_fn, 'w') as out:
            out.write(osr_wkt(src_srs, True))
        
        if os.path.exists(dst_fn):
            return 0
    except Exception:
        pass
    return -1


def epsg_from_input(in_srs: str) -> Tuple[Any, Any]:
    """Get the EPSG(s) from SRS suitable as input to SetFromUserInput.

    Returns:
        list: [horz_epsg, vert_epsg]
    """
    
    src_vert = None
    
    ## Pre-cleaning
    if np.any(['geoid' in x for x in in_srs.split('+')]):
        in_srs = '+'.join(in_srs.split('+')[:-1])
        
    if in_srs.startswith('ESRI'):
        esri_split = in_srs.split('+')
        if len(esri_split) > 1:
            vert_epsg = in_srs.split('+')[1]
            src_vert = vert_epsg
        in_srs = esri_split[0]
    
    src_srs = osr.SpatialReference()
    src_srs.SetFromUserInput(in_srs)

    ## Horizontal
    if src_srs.IsGeographic():
        cstype = 'GEOGCS'
    else:
        cstype = 'PROJCS'

    src_srs.AutoIdentifyEPSG()
    an = src_srs.GetAuthorityName(cstype)
    ac = src_srs.GetAuthorityCode(cstype)

    if ac is None:
        ## Fallback to Proj4 or WKT if no EPSG code found
        src_horz = src_srs.ExportToProj4()
        if not src_horz:
            src_horz = src_srs.ExportToWkt()
    else:
        src_horz = f'{an}:{ac}'
        
    ## Vertical
    if src_srs.IsVertical():
        csvtype = 'VERT_CS'
        src_vert = src_srs.GetAuthorityCode(csvtype)

    return src_horz, src_vert


def osr_parse_srs(src_srs, return_vertcs=True):
    """Parse an OSR SRS object and return a proj string."""
    
    if src_srs is None:
        return None

    if src_srs.IsLocal():
        return src_srs.ExportToWkt()
    
    cstype = 'GEOGCS' if src_srs.IsGeographic() else 'PROJCS'
        
    src_srs.AutoIdentifyEPSG()
    an = src_srs.GetAuthorityName(cstype)
    ac = src_srs.GetAuthorityCode(cstype)

    vn, vc = None, None
    if src_srs.IsVertical():
        csvtype = 'VERT_CS'
        vn = src_srs.GetAuthorityName(csvtype)
        vc = src_srs.GetAuthorityCode(csvtype)

    if an is not None and ac is not None:
        if vn is not None and vc is not None:
            return f'{an}:{ac}+{vc}'
        else:
            return f'{an}:{ac}'
    else:
        return src_srs.ExportToProj4() or None


###############################################################################    
## OGR
###############################################################################
def ogr_or_gdal(osgeo_fn: str) -> int:
    """Check if file is OGR (1), GDAL (2), or neither (-1)."""
    
    try:
        ## Try OGR
        src_ds = ogr.Open(osgeo_fn)
        if src_ds is not None:
            src_ds = None
            return 1
        
        ## Try GDAL
        src_ds = gdal.Open(osgeo_fn)
        if src_ds is not None:
            src_ds = None
            return 2
            
        return -1
    except Exception:
        utils.echo_error_msg(f'{osgeo_fn} does not appear to be an osgeo (gdal/ogr) file')
        return -1


def ogr_fext(src_drv_name: str) -> str:
    """Find the common file extension given a OGR driver name."""
    
    try:
        drv = ogr.GetDriverByName(src_drv_name)
        fexts = drv.GetMetadataItem(gdal.DMD_EXTENSIONS)
        if fexts is not None:
            return fexts.split()[0]
    except Exception:
        pass

    ## Fallback map
    mapping = {
        'gtiff': 'tif',
        'hfa': 'img',
        'gmt': 'grd',
        'netcdf': 'nc'
    }
    return mapping.get(src_drv_name.lower(), 'gdal')


def ogr_get_srs(src_ogr: str) -> Optional[str]:
    """Get the SRS (as WKT) from an OGR file."""
    
    src_ds = ogr.Open(src_ogr, 0)
    if src_ds:
        layer = src_ds.GetLayer()
        src_srs = layer.GetSpatialRef()
        src_ds = None
        if src_srs:
            return src_srs.ExportToWkt()
    return None


def ogr_clip(src_ogr_fn, dst_region=None, layer=None, fmt='GPKG',
             overwrite=False, verbose=True):
    """Clip an OGR file to `dst_region` using ogr2ogr CLI."""
    
    if dst_region is None:
        return src_ogr_fn
    
    dst_ogr_bn = os.path.splitext(src_ogr_fn)[0]
    dst_ogr_fn = f"{dst_ogr_bn}_{dst_region.format('fn')}.{ogr_fext(fmt)}"
    
    if not os.path.exists(dst_ogr_fn) or overwrite:
        layer_opt = layer if layer is not None else ''
        cmd = (f"ogr2ogr -nlt PROMOTE_TO_MULTI {dst_ogr_fn} {src_ogr_fn} "
               f"-clipsrc {dst_region.format('te')} {layer_opt}")
        utils.run_cmd(cmd, verbose=verbose)

    return dst_ogr_fn


def ogr_clip2(src_ogr_fn, dst_region=None, layer=None,
              overwrite=False, verbose=True):
    """Clip an OGR file to `dst_region` using Python bindings."""
    
    dst_ogr_bn = os.path.splitext(src_ogr_fn)[0]
    dst_ogr_fn = f"{dst_ogr_bn}_{dst_region.format('fn')}.gpkg"
    
    if not os.path.exists(dst_ogr_fn) or overwrite:
        src_ds = ogr.Open(src_ogr_fn)
        if layer is not None:
            src_layer = src_ds.GetLayer(layer)
        else:
            src_layer = src_ds.GetLayer()
            
        region_ogr = f"region_{dst_region.format('fn')}.shp"
        dst_region.export_as_ogr(region_ogr)
        
        region_ds = ogr.Open(region_ogr)
        region_layer = region_ds.GetLayer()

        driver = ogr.GetDriverByName('GPKG')
        dst_ds = driver.CreateDataSource(dst_ogr_fn)
        dst_layer = dst_ds.CreateLayer(
            layer if layer is not None else 'clipped',
            geom_type=ogr.wkbPolygon
        )

        src_layer.Clip(region_layer, dst_layer)

        ## Cleanup
        src_ds = None
        region_ds = None
        dst_ds = None
        utils.remove_glob(f"{os.path.splitext(region_ogr)[0]}.*")
        
    return dst_ogr_fn


def ogr_polygonize_line_to_region(src_ogr, dst_ogr, region=None, include_landmask=True,
                                  landmask_is_watermask=False, line_buffer=0.0000001):
    """Polygonize an OSM coastline LineString to the given region."""
    
    line_ds = ogr.Open(src_ogr)
    line_layer = line_ds.GetLayer()
    
    ## Determine processing region
    if region is not None and region.is_valid():
        region_geom = region.export_as_geom()
    else:
        region_geom = regions.Region().from_list(line_layer.GetExtent()).export_as_geom()
    
    ## Create the output layer
    driver = ogr.GetDriverByName("ESRI Shapefile")
    output_ds = driver.CreateDataSource(dst_ogr)
    output_layer = output_ds.CreateLayer(
        "split_polygons",
        line_layer.GetSpatialRef(),
        ogr.wkbMultiPolygon
    )
    output_layer.CreateField(ogr.FieldDefn('watermask', ogr.OFTInteger))
    
    for line_feature_layer in line_ds:
        line_type = line_feature_layer.GetGeomType()
        
        ## Process LineStrings
        if line_type == ogr.wkbLineString:
            line_geometries = ogr_union_geom(
                line_feature_layer,
                ogr.wkbMultiLineString
            )
            poly_line = line_geometries.Buffer(line_buffer)
            split_geoms = region_geom.Difference(poly_line)
            
            for split_geom in split_geoms:
                is_water = []
                ## Check intersection logic (ray casting simulation)
                for line_geometry in line_geometries:
                    if split_geom.Intersects(line_geometry.Buffer(line_buffer)):
                        point_count = line_geometry.GetPointCount()
                        poly_center = split_geom.Centroid()
                        pcx = poly_center.GetX()
                        pcy = poly_center.GetY()
                        
                        for point_n in range(point_count - 1):
                            xb = line_geometry.GetX(point_n)
                            yb = line_geometry.GetY(point_n)
                            xe = line_geometry.GetX(point_n + 1)
                            ye = line_geometry.GetY(point_n + 1)

                            s = (xe - xb) * (pcy - yb) > (ye - yb) * (pcx - xb)
                            is_water.append(s)

                ## Determine if polygon is water or land based on voting
                if all(is_water):
                    s_val = True
                elif not any(is_water):
                    s_val = False
                else:
                    s_val = np.count_nonzero(is_water) > len(is_water) / 2

                out_feature = ogr.Feature(output_layer.GetLayerDefn())
                out_feature.SetGeometry(split_geom)

                if landmask_is_watermask:
                    s_val = not s_val
                
                if not s_val:
                    ## Water
                    out_feature.SetField('watermask', 1)
                    output_layer.CreateFeature(out_feature)

                if include_landmask and s_val:
                    ## Land
                    out_feature.SetField('watermask', 0)
                    output_layer.CreateFeature(out_feature)
        
        ## Process Polygons (Already processed areas)
        elif line_type == ogr.wkbMultiPolygon:
            for line_feature in line_feature_layer:
                line_geometry = line_feature.geometry()
                line_geometry = ogr.ForceTo(line_geometry, ogr.wkbLinearRing)
                
                ## Update existing output features
                for feature in output_layer:
                    feature_geom = feature.geometry()
                    if feature_geom.Contains(line_geometry):
                        feature_geoms = feature_geom.Difference(line_geometry)
                        feature.SetGeometry(feature_geoms)
                        output_layer.SetFeature(feature)
                        
                out_feature = ogr.Feature(output_layer.GetLayerDefn())
                out_feature.SetGeometry(line_geometry)
                out_feature.SetField('watermask', 1)
                output_layer.CreateFeature(out_feature)

    line_ds = None
    output_ds = None
    return dst_ogr


def ogr_union_geom(src_layer, geom_type=ogr.wkbMultiPolygon, verbose=True):        
    """Union all geometries in a layer."""
    
    multi = ogr.Geometry(geom_type)
    feats = src_layer.GetFeatureCount()
    
    for f in src_layer:
        multi.AddGeometry(f.geometry())
        
    if verbose:
        utils.echo_msg(f'unioned {feats} features')
    return multi


def ogr_is_empty(src_ogr: str, dn: str = 'ESRI Shapefile') -> bool:
    """Check if the OGR file is empty."""
    
    driver = ogr.GetDriverByName(dn)
    ds = driver.Open(src_ogr, 0)

    if ds is not None:
        layer = ds.GetLayer()
        fc = layer.GetFeatureCount()
        ds = None
        return fc == 0
    return True

## Alias for backward compatibility
ogr_empty_p = ogr_is_empty
    
def ogr_polygonize(src_ds, dst_srs='epsg:4326', ogr_format='ESRI Shapefile',
                   band=1, verbose=True):
    """Polygonize a multi-band raster."""
    
    dst_layer = f"{utils.fn_basename2(src_ds.GetDescription())}_ply"
    dst_vector = f"{dst_layer}.{ogr_fext(ogr_format)}"
    
    utils.remove_glob(f'{dst_layer}.*')
    osr_prj_file(f'{dst_layer}.prj', gdal_infos(src_ds)['proj'])
    
    driver = ogr.GetDriverByName(ogr_format)
    ds = driver.CreateDataSource(dst_vector)
    
    if ds is None:
        return None, ogr_format
        
    layer = ds.CreateLayer(dst_layer, None, ogr.wkbMultiPolygon)
    layer.CreateField(ogr.FieldDefn('DN', ogr.OFTInteger))

    ## Add fields from raster metadata
    this_band = src_ds.GetRasterBand(band)
    this_band_md = this_band.GetMetadata()
    
    ## Capitalize keys and add fields
    formatted_md = {k.title(): v for k, v in this_band_md.items()}
    field_names = [field.name for field in layer.schema]
    
    for k in formatted_md.keys():
        if k[:9] not in field_names:
            layer.CreateField(ogr.FieldDefn(k[:9], ogr.OFTString))

    if verbose:
        utils.echo_msg(f'Polygonizing {this_band.GetDescription()} mask...')

    status = gdal.Polygonize(
        this_band,
        None,
        layer,
        layer.GetLayerDefn().GetFieldIndex('DN'),
        [],
        callback=gdal.TermProgress if verbose else None
    )
    
    ## Clean up features with DN=0 (No Data usually)
    for feat in layer:
        if feat.GetField('DN') == 0:
            layer.DeleteFeature(feat.GetFID())
            
    ds = None
    return dst_layer, ogr_format


def ogr2gdal_mask(mask_fn, region=None, x_inc=None, y_inc=None,
                  dst_srs='epsg:4326', invert=True, verbose=True,
                  temp_dir='./'):
    """Create a raster mask from an OGR vector."""
    
    dst_fn = utils.make_temp_fn(
        f'{mask_fn}.tif', region=region, inc=x_inc, temp_dir=temp_dir
    )

    if os.path.exists(dst_fn):
        return dst_fn
    
    if region is not None and x_inc is not None and y_inc is not None:
        xcount, ycount, dst_gt = region.geo_transform(
            x_inc=x_inc, y_inc=y_inc, node='grid'
        )
        
        if xcount <= 0 or ycount <= 0:
            utils.echo_error_msg(f'Invalid grid size {xcount}x{ycount} for mask generation.')
            sys.exit(-1)

        ds_config = gdal_set_infos(
            xcount, ycount, xcount * ycount, dst_gt,
            dst_srs, gdal.GDT_Float32, 0, 'GTiff', {}, 1
        )
        gdal_nan(ds_config, dst_fn, nodata=0)
        
        clip_layer = ogr_get_layer_name(mask_fn)
        if clip_layer is None:
            clip_layer = os.path.splitext(os.path.basename(mask_fn))[0]
            
        burn_val = 1
        invert_flag = '-i' if invert else ''
        
        ## Use gdal_rasterize
        cmd = f'gdal_rasterize -burn {burn_val} -l {clip_layer} "{mask_fn}" "{dst_fn}" {invert_flag}'
        utils.run_cmd(cmd, verbose=verbose)
        
        return dst_fn
    
    return None


def ogr_geoms2ogr(geoms, out, dst_srs='epsg:4326', ogr_format='ESRI Shapefile'):
    """Save a list of geometries to an OGR file."""
    
    dirname = os.path.dirname(out) or '.'
    dst_layer = os.path.basename(utils.fn_basename2(out))
    dst_vector = os.path.join(dirname, f'{dst_layer}.{ogr_fext(ogr_format)}')
    
    if os.path.exists(out):
        utils.remove_glob(f'{utils.fn_basename2(out)}*')
        
    osr_prj_file(f'{utils.fn_basename2(out)}.prj', dst_srs)
    driver = ogr.GetDriverByName(ogr_format)
    ds = driver.CreateDataSource(dst_vector)
    
    if ds is not None: 
        layer = ds.CreateLayer(dst_layer, None, ogr.wkbMultiPolygon)
        
        for g in geoms:
            out_feature = ogr.Feature(layer.GetLayerDefn())
            out_feature.SetGeometry(g)
            layer.CreateFeature(out_feature)
            
        ds = None
    return out


def ogr_get_layer_name(ogr_fn):
    """Get the name of the first layer in an OGR file."""
    
    ds = ogr.Open(ogr_fn, 0)
    layer_name = None
    if ds is not None:
        layer = ds.GetLayer()
        if layer:
            layer_name = layer.GetName()
        ds = None
    return layer_name


def ogr_geoms(ogr_fn):
    """Extract all geometries from an OGR file as WKT."""
    
    out_geoms = []
    ds = ogr.Open(ogr_fn, 0)
    if ds is not None:
        layer = ds.GetLayer()
        for f in layer:
            geom = f.GetGeometryRef()
            if geom is not None and not geom.IsEmpty():
                out_geoms.append(geom.ExportToWkt())
        ds = None
    return out_geoms

    
###############################################################################        
## GDAL
###############################################################################
class gdal_datasource:
    """Context manager for GDAL datasets.
    
    Usage:
        with gdal_datasource('input.tif') as src_ds:
            # do something
    """
    
    def __init__(self, src_gdal=None, update=False):
        self.src_gdal = src_gdal
        self.update = update
        self.src_ds = None

    def __enter__(self):
        ## If already a dataset object
        if isinstance(self.src_gdal, gdal.Dataset):
            self.src_ds = self.src_gdal
            return self.src_ds

        ## If a path string
        if utils.str_or(self.src_gdal) is not None and (
            os.path.exists(self.src_gdal) or 
            utils.fn_url_p(self.src_gdal) or 
            ':' in self.src_gdal
        ):
            try:
                if self.update:
                    self.src_ds = gdal.OpenEx(
                        self.src_gdal,
                        gdal.OF_RASTER | gdal.OF_UPDATE,
                        open_options=['IGNORE_COG_LAYOUT_BREAK=YES']
                    )
                else:
                    self.src_ds = gdal.Open(self.src_gdal, gdal.GA_ReadOnly)
            except Exception:
                self.src_ds = None

        if self.src_ds is not None and self.update:
             ## Build overviews automatically in update mode if strictly needed
             self.src_ds.BuildOverviews('AVERAGE', [])
            
        return self.src_ds

    def __exit__(self, exc_type, exc_value, exc_traceback):
        ## Only close if we opened it (if string passed)
        if not isinstance(self.src_gdal, gdal.Dataset) and self.src_ds:
            if self.update:
                try:
                    self.src_ds.FlushCache()
                except Exception:
                    pass
            self.src_ds = None


def gdal_get_srs(src_gdal: str) -> Optional[str]:
    """Get the SRS (as WKT) from a GDAL file."""
    
    with gdal_datasource(src_gdal) as src_ds:
        if src_ds is not None:
            src_srs = src_ds.GetSpatialRef()
            if src_srs is not None:
                return src_srs.ExportToWkt()
    return None

        
def gdal_fext(src_drv_name: str) -> str:
    """Find the common file extension given a GDAL driver name."""
    
    try:
        drv = gdal.GetDriverByName(src_drv_name)
        if drv and drv.GetMetadataItem(gdal.DCAP_RASTER):
            fexts = drv.GetMetadataItem(gdal.DMD_EXTENSIONS)
            if fexts:
                return fexts.split()[0]
    except Exception:
        pass

    mapping = {
        'gtiff': 'tif',
        'hfa': 'img',
        'gmt': 'grd',
        'netcdf': 'nc',
        'vrt': 'vrt'
    }
    return mapping.get(src_drv_name.lower(), 'gdal')

    
def gdal_set_infos(nx, ny, nb, geoT, proj, dt, ndv, fmt, md, rc) -> Dict:
    """Create a dictionary of GDAL dataset configuration."""
    
    return {
        'nx': nx, 'ny': ny, 'nb': nb, 'geoT': geoT, 'proj': proj,
        'dt': dt, 'ndv': ndv, 'fmt': fmt, 'metadata': md,
        'raster_count': rc
    }


def gdal_infos(src_gdal, region=None, scan=False, band=1) -> Dict:
    """Gather and return info about a src_gdal file."""
    
    ds_config = {}
    
    with gdal_datasource(src_gdal) as src_ds:
        if src_ds is None:
            utils.echo_warning_msg(f'Could not load raster {src_gdal}')
            return ds_config

        gt = src_ds.GetGeoTransform()
        
        ## Calculate Source Window (srcwin)
        if region is not None:
            srcwin = region.srcwin(gt, src_ds.RasterXSize, src_ds.RasterYSize)
        else:
            srcwin = (0, 0, src_ds.RasterXSize, src_ds.RasterYSize)

        ## Calculate Dest GeoTransform based on srcwin
        dst_gt = (
            gt[0] + (srcwin[0] * gt[1]),
            gt[1],
            0.,
            gt[3] + (srcwin[1] * gt[5]),
            0.,
            gt[5]
        )
        
        src_band = src_ds.GetRasterBand(band)
        if src_band is not None:
            ds_config = {
                'nx': srcwin[2],
                'ny': srcwin[3],
                'nb': srcwin[2] * srcwin[3],
                'geoT': dst_gt,
                'proj': src_ds.GetProjectionRef(),
                'dt': src_band.DataType,
                'dtn': gdal.GetDataTypeName(src_band.DataType),
                'ndv': src_band.GetNoDataValue(),
                'fmt': src_ds.GetDriver().ShortName,
                'metadata': src_ds.GetMetadata(),
                'raster_count': src_ds.RasterCount,
            }
            
            if ds_config['ndv'] is None:
                ds_config['ndv'] = -9999

            if scan:
                src_arr = src_band.ReadAsArray(srcwin[0], srcwin[1], srcwin[2], srcwin[3])
                ## Mask nodata
                if ds_config['ndv'] is not None:
                    # Float check for NaN/Nodata
                    is_nodata = (src_arr == ds_config['ndv'])
                    if np.issubdtype(src_arr.dtype, np.floating):
                        is_nodata |= np.isnan(src_arr)
                    
                    ## Compute min/max ignoring nodata
                    if not np.all(is_nodata):
                        ## Mask array or use gdal compute
                        ds_config['zr'] = src_band.ComputeRasterMinMax()
                    else:
                        ds_config['zr'] = [np.nan, np.nan]
                else:
                    ds_config['zr'] = src_band.ComputeRasterMinMax()
        else:
            utils.echo_warning_msg(f'Invalid band {band} for data source {src_gdal}')

    return ds_config


def gdal_copy_infos(src_config: Dict) -> Dict:
    """Deep copy of src_config dictionary."""
    
    return src_config.copy()


def gdal_set_srs(src_gdal, src_srs='epsg:4326', verbose=True):
    """Set the projection of src_gdal."""
    
    if src_srs and '+geoid' in src_srs:
        src_srs = '+'.join(src_srs.split('+')[:-1])
        
    with gdal_datasource(src_gdal, update=True) as src_ds:    
        if src_ds is not None and src_srs is not None:
            try:
                src_ds.SetProjection(osr_wkt(src_srs))
                return 0
            except Exception:
                if verbose:
                    utils.echo_warning_msg(f'Could not set projection {src_srs}')
    return None


def gdal_get_ndv(src_gdal, band=1):
    """Get the nodata value of src_gdal."""
    
    with gdal_datasource(src_gdal) as src_ds:    
        if src_ds:
            return src_ds.GetRasterBand(band).GetNoDataValue()
    return None


def gdal_set_ndv(src_gdal, ndv=-9999, convert_array=False, verbose=True):
    """Set the nodata value of gdal datasource."""
    
    with gdal_datasource(src_gdal, update=True) as src_ds:
        if src_ds is None:
            return None
            
        ds_config = gdal_infos(src_ds)
        curr_nodata = ds_config['ndv']
        
        if verbose:
            utils.echo_msg(f'Setting nodata value from {curr_nodata} to {ndv}')
                            
        for band_idx in range(1, src_ds.RasterCount + 1):
            this_band = src_ds.GetRasterBand(band_idx)
            ## GDAL doesn't like setting existing NDV directly sometimes, delete first
            this_band.DeleteNoDataValue()
            this_band.SetNoDataValue(ndv)

            if convert_array:
                arr = this_band.ReadAsArray()
                ## Handle NaN if current is nan, else value compare
                if curr_nodata is not None and np.isnan(curr_nodata):
                    arr[np.isnan(arr)] = ndv
                elif curr_nodata is not None:
                    arr[arr == curr_nodata] = ndv
                this_band.WriteArray(arr)            
    return 0


def flatten_no_data_zones(src_arr, src_config, size_threshold=1, verbose=True):
    """Flatten (fill) nodata areas larger than `size_threshold`."""
    
    def expand_for(arr, shiftx=1, shifty=1):
        """Morphological dilation helper"""
        
        ## Note: scipy.ndimage.binary_dilation is more efficient usually
        if 'scipy' in sys.modules:
            struct = scipy.ndimage.generate_binary_structure(2, 2)
            return scipy.ndimage.binary_dilation(arr, structure=struct, iterations=max(shiftx, shifty))
        
        ## Fallback manual implementation
        arr_b = arr.copy().astype(bool)
        rows, cols = arr.shape
        for i in range(rows):
            for j in range(cols):
                if arr[i, j]:
                    i_min, i_max = max(i - shifty, 0), min(i + shifty + 1, rows)
                    j_min, j_max = max(j - shiftx, 0), min(j + shiftx + 1, cols)
                    arr_b[i_min:i_max, j_min:j_max] = True
        return arr_b
    
    ## Normalize nodata to NaN
    src_arr = src_arr.astype(float)
    if src_config['ndv'] is not None:
        src_arr[src_arr == src_config['ndv']] = np.nan
    
    ## Mask array: 1 where NaN (void)
    msk_arr = np.zeros_like(src_arr)
    msk_arr[np.isnan(src_arr)] = 1

    if 'scipy.ndimage' not in sys.modules:
        utils.echo_warning_msg("Scipy not available, skipping void flattening.")
        return src_arr

    ## Label voids
    l, n = scipy.ndimage.label(msk_arr)

    ## Count sizes
    mn = scipy.ndimage.sum_labels(msk_arr, labels=l, index=np.arange(1, n + 1))
    
    with utils.ccp(desc=f'Flattening data voids > {size_threshold}', leave=verbose) as pbar:
        for i in range(n):
            if mn[i] >= size_threshold:
                ## Expand mask to find boundary pixels to interpolate from
                current_void_mask = (l == (i + 1))
                expanded_mask = expand_for(current_void_mask)
                
                ## Calculate fill value (5th percentile of surrounding valid data)
                ## Boundary is expanded minus original void
                boundary_mask = expanded_mask & ~current_void_mask
                
                if np.any(boundary_mask):
                    flat_value = np.nanpercentile(src_arr[boundary_mask], 5)
                    src_arr[current_void_mask] = flat_value
            pbar.update()

    ## Restore Nodata
    src_arr[np.isnan(src_arr)] = src_config['ndv']

    return src_arr


def cudem_flatten_no_data_zones(src_dem, dst_dem=None, band=1,
                                size_threshold=1, verbose=True):
    """Wrapper to flatten nodata zones on a GDAL file."""
    
    update_mode = (dst_dem is None)
    
    with gdal_datasource(src_dem, update=update_mode) as src_ds:
        if src_ds is None:
            return src_dem, -1
            
        src_band = src_ds.GetRasterBand(band)
        src_arr = src_band.ReadAsArray()
        src_config = gdal_infos(src_ds)
        
        src_arr = flatten_no_data_zones(
            src_arr, src_config, size_threshold=size_threshold, verbose=verbose
        )

        if dst_dem is None:
            src_band.WriteArray(src_arr)
        else:
            gdal_write(src_arr, dst_dem, src_config)

    return dst_dem if dst_dem is not None else src_dem, 0


def gdal_mem_ds(ds_config, name='MEM', bands=1, src_srs=None, co=None):
    """Create a temporary GDAL MEM dataset."""
    
    if co is None:
        co = [] # MEM driver options usually empty or different
        
    mem_driver = gdal.GetDriverByName('MEM')
    mem_ds = mem_driver.Create(
        name,
        ds_config['nx'],
        ds_config['ny'],
        bands,
        ds_config['dt'],
        options=co
    )
    
    if mem_ds is not None:
        mem_ds.SetGeoTransform(ds_config['geoT'])
        
        ## Projection
        if src_srs is None:
            if ds_config['proj'] is not None:
                mem_ds.SetProjection(ds_config['proj'])
        else:
            mem_ds.SetProjection(src_srs)

        ## Nodata
        for band in range(1, bands + 1):
            mem_band = mem_ds.GetRasterBand(band)
            mem_band.SetNoDataValue(ds_config['ndv'])
        
    return mem_ds


def generate_mem_ds(ds_config, band_data=None, srcwin=None, return_array=False,
                    interpolation='nearest'):
    """Generate in-memory dataset, optionally interpolating data."""
    
    tmp_band_data = band_data.copy()
    if ds_config['ndv'] is not None:
        tmp_band_data[tmp_band_data == ds_config['ndv']] = np.nan
        
    gt = ds_config['geoT']
    if srcwin is None:
        srcwin = (0, 0, ds_config['nx'], ds_config['ny'])
        
    if np.all(np.isnan(tmp_band_data)):
        return None

    ## Interpolate for neighborhood calculations
    if interpolation in ['nearest', 'linear', 'cubic', 'flats'] and 'scipy' in sys.modules:
        if np.any(np.isnan(band_data)):
            if interpolation == 'flats':
                tmp_band_data = flatten_no_data_zones(
                    tmp_band_data, ds_config, size_threshold=1, verbose=True
                )
            else:
                point_indices = np.nonzero(~np.isnan(tmp_band_data))
                if len(point_indices[0]):
                    point_values = tmp_band_data[point_indices]
                    xi, yi = np.mgrid[0:srcwin[3], 0:srcwin[2]]
                    try:
                        tmp_band_data = scipy.interpolate.griddata(
                            tuple(point_indices), point_values,
                            (xi, yi), method=interpolation
                        )
                    except Exception:
                        pass

    if return_array:
        return tmp_band_data

    ## Generate MEM datasource
    dst_gt = (
        gt[0] + (srcwin[0] * gt[1]),
        gt[1],
        0.,
        gt[3] + (srcwin[1] * gt[5]),
        0.,
        gt[5]
    )
    
    srcwin_config = gdal_set_infos(
        srcwin[2], srcwin[3], srcwin[2]*srcwin[3],
        dst_gt, ds_config['proj'], ds_config['dt'],
        ds_config['ndv'], 'GTiff', {}, 1
    )
    
    srcwin_ds = gdal_mem_ds(srcwin_config, name='MEM', bands=1, src_srs=None)
    srcwin_band = srcwin_ds.GetRasterBand(1)
    srcwin_band.SetNoDataValue(ds_config['ndv'])
    
    ## Restore Nodata for writing
    if ds_config['ndv'] is not None:
        tmp_band_data[np.isnan(tmp_band_data)] = ds_config['ndv']
        
    srcwin_band.WriteArray(tmp_band_data)
    srcwin_ds.FlushCache()
    
    return srcwin_ds


def gdal_nan(ds_config, outfile, nodata=None):
    """Create a nodata grid."""
    
    if nodata is None:
        nodata = np.nan
        
    null_array = np.zeros((ds_config['ny'], ds_config['nx']))
    if nodata != 0:
        null_array[:] = nodata

    ds_config['ndv'] = nodata        
    gdal_write(null_array, outfile, ds_config)


def gdal_extract_band(src_gdal, dst_gdal, band=1, exclude=None,
                      srcwin=None, inverse=False):
    """Extract a band from src_gdal and write to dst_gdal."""
    
    band = utils.int_or(band, 1)
    
    with gdal_datasource(src_gdal) as src_ds:
        if src_ds is None:
            return None, -1
            
        try:
            ds_config = gdal_infos(src_ds)
            ds_band = src_ds.GetRasterBand(band)
            
            if srcwin is not None:
                ds_array = ds_band.ReadAsArray(*srcwin)
                ## Adjust GeoTransform for window
                xy_origin = utils._pixel2geo(
                    srcwin[0], srcwin[1], ds_config['geoT'], node='grid'
                )
                ds_config['geoT'] = (
                    xy_origin[0], ds_config['geoT'][1], ds_config['geoT'][2],
                    xy_origin[1], ds_config['geoT'][4], ds_config['geoT'][5]
                )
                ds_config['nx'] = srcwin[2]
                ds_config['ny'] = srcwin[3]
                ds_config['nb'] = srcwin[2] * srcwin[3]
            else:
                ds_array = ds_band.ReadAsArray()
                
            if ds_array is None:
                raise ValueError("ReadAsArray returned None")
                
        except Exception as e:
            utils.echo_error_msg(f'Could not read datasource {src_gdal}, {e}')
            return None, -1

    if ds_config['ndv'] is None:
        ds_config['ndv'] = -9999
    
    if exclude:
        for key in exclude:
            ds_array[ds_array == key] = ds_config['ndv']

    if inverse:
        ## Avoid division by zero warnings
        with np.errstate(divide='ignore'):
            ds_array = 1 / ds_array
        ds_array[np.isinf(ds_array)] = ds_config['ndv']

    return gdal_write(ds_array, dst_gdal, ds_config)


def gdal_cut(src_gdal, src_region, dst_gdal, node='pixel',
             verbose=True, co=None):
    """Cut src_ds to src_region and output dst_gdal using CreateCopy."""
    
    if co is None:
        co = ["COMPRESS=DEFLATE", "TILED=YES"]
        
    status = -1
    with gdal_datasource(src_gdal) as src_ds:
        if src_ds is not None:
            ds_config = gdal_infos(src_ds)
            gt = ds_config['geoT']
            srcwin = src_region.srcwin(gt, src_ds.RasterXSize, src_ds.RasterYSize, node=node)
            
            ## Destination GeoTransform
            dst_gt = (
                gt[0] + (srcwin[0] * gt[1]),
                gt[1], 0.,
                gt[3] + (srcwin[1] * gt[5]),
                0., gt[5]
            )
            
            out_ds_config = gdal_set_infos(
                srcwin[2], srcwin[3], srcwin[2] * srcwin[3],
                dst_gt, ds_config['proj'], ds_config['dt'],
                ds_config['ndv'], ds_config['fmt'],
                ds_config['metadata'], ds_config['raster_count']
            )

            in_bands = src_ds.RasterCount
            mem_ds = gdal_mem_ds(out_ds_config, bands=in_bands)
            
            if mem_ds is not None:
                for band in range(1, in_bands + 1):
                    this_band = mem_ds.GetRasterBand(band)
                    that_band = src_ds.GetRasterBand(band)
                    this_band.SetDescription(that_band.GetDescription())
                    this_band.SetMetadata(that_band.GetMetadata())
                    this_band.WriteArray(
                        src_ds.GetRasterBand(band).ReadAsArray(*srcwin)
                    )
                    mem_ds.FlushCache()

                driver = gdal.GetDriverByName(ds_config['fmt'])
                driver.CreateCopy(dst_gdal, mem_ds, 0, options=co)
                status = 0
                mem_ds = None

    return dst_gdal, status


def gdal_clip(src_gdal, dst_gdal, src_ply=None, invert=False,
              verbose=True, cache_dir=None):
    """Clip DEM to polygon `src_ply` using rasterize."""
    
    gi = gdal_infos(src_gdal)
    g_region = regions.Region().from_geo_transform(
        geo_transform=gi['geoT'], x_count=gi['nx'], y_count=gi['ny']
    )
    
    tmp_ply = utils.make_temp_fn('tmp_clp_ply.shp', temp_dir=cache_dir)    
    
    ## Clip vector first
    cmd = (f'ogr2ogr "{tmp_ply}" "{src_ply}" -clipsrc {g_region.format("ul_lr")} '
           '-nlt POLYGON -skipfailures -makevalid')
    utils.run_cmd(cmd, verbose=verbose)
    
    if gi is not None and src_ply is not None:
        shutil.copyfile(src_gdal, dst_gdal)
        
        with gdal_datasource(src_gdal) as src_ds:
            band_count = src_ds.RasterCount if src_ds else 1
                
        invert_flag = '-i' if invert else ''
        layer_name = os.path.splitext(os.path.basename(tmp_ply))[0]
        
        ## Burn nodata into destination
        gr_cmd = (f'gdal_rasterize -b {band_count} -burn {gi["ndv"]} '
                  f'-l {layer_name} "{tmp_ply}" "{dst_gdal}" {invert_flag}')
        
        out, status = utils.run_cmd(gr_cmd, verbose=verbose)
        
        ## Cleanup temp vector
        utils.remove_glob(f'{os.path.splitext(tmp_ply)[0]}.*')
        return out, status
    else:
        utils.remove_glob(f'{os.path.splitext(tmp_ply)[0]}.*')
        return None, -1


def crop(src_gdal, dst_gdal):
    """Crop the src_gdal file to dst_gdal by trimming nodata value borders."""
    
    status = -1
    with gdal_datasource(src_gdal) as src_ds: 
        if src_ds is None:
            return None, status
            
        ds_config = gdal_infos(src_ds)
        ds_arr = src_ds.GetRasterBand(1).ReadAsArray()
        
        ## Mask nodata as nan for finding valid area
        if ds_config['ndv'] is not None:
            ds_arr = ds_arr.astype(float)
            ds_arr[ds_arr == ds_config['ndv']] = np.nan
            
        nans = np.isnan(ds_arr)
        
        if np.all(nans):
            utils.echo_warning_msg("Raster is fully nodata/nan.")
            return None, -1

        nancols = np.all(nans, axis=0)
        nanrows = np.all(nans, axis=1)

        firstcol = nancols.argmin()
        firstrow = nanrows.argmin()        
        lastcol = len(nancols) - nancols[::-1].argmin()
        lastrow = len(nanrows) - nanrows[::-1].argmin()

        dst_arr = ds_arr[firstrow:lastrow, firstcol:lastcol]
        
        ## Restore nodata
        if ds_config['ndv'] is not None:
            dst_arr[np.isnan(dst_arr)] = ds_config['ndv']
            
        gt = ds_config['geoT']
        dst_x_origin = gt[0] + (gt[1] * firstcol)
        dst_y_origin = gt[3] + (gt[5] * firstrow)
        
        ds_config['geoT'] = (dst_x_origin, gt[1], 0.0, dst_y_origin, 0.0, gt[5])
        ds_config['nx'] = int(lastcol - firstcol)
        ds_config['ny'] = int(lastrow - firstrow)
        ds_config['nb'] = ds_config['nx'] * ds_config['ny']
        
        _, status = gdal_write(dst_arr, dst_gdal, ds_config)

    return status


def gdal_write(src_arr, dst_gdal, ds_config, dst_fmt='GTiff',
               co=None, max_cache=False, verbose=False):
    """Write src_arr to gdal file."""
    
    if co is None:
        co = ["COMPRESS=DEFLATE", "TILED=YES"]
        
    driver = gdal.GetDriverByName(dst_fmt)
    if os.path.exists(dst_gdal):
        try:
            driver.Delete(dst_gdal)
        except Exception:
            utils.remove_glob(dst_gdal)

    dirname = os.path.dirname(dst_gdal)
    if dirname and not os.path.exists(dirname):
        os.makedirs(dirname)
            
    if max_cache:
        gdal.SetCacheMax(2**30)

    ds = driver.Create(
        dst_gdal,
        ds_config['nx'],
        ds_config['ny'],
        1,
        ds_config['dt'],
        options=co
    )

    if ds is not None:
        ds.SetGeoTransform(ds_config['geoT'])
        try:
            ds.SetProjection(ds_config['proj'])
        except Exception:
            if verbose:
                utils.echo_warning_msg(f"Could not set projection {ds_config['proj']}")
                
        ds.GetRasterBand(1).SetNoDataValue(ds_config['ndv'])
        if src_arr is not None:
            ds.GetRasterBand(1).WriteArray(src_arr)
            
        ds = None        
        return dst_gdal, 0
    else:
        return None, -1


def gdal_yield_srcwin(src_gdal, n_chunk=10, step=5, verbose=False):
    """Yield source windows (srcwin) from dataset."""
    
    ds_config = gdal_infos(src_gdal)
    
    x_chunk = n_chunk
    
    total = math.ceil(ds_config['ny']/step) * math.ceil(ds_config['nx']/step)
    
    with utils.ccp(total=total, desc=f'Chunking {src_gdal}') as pbar:
        while True:
            y_chunk = n_chunk
            while True:
                pbar.update()
                
                this_x_chunk = min(x_chunk, ds_config['nx'])
                this_y_chunk = min(y_chunk, ds_config['ny'])
                
                this_x_origin = x_chunk - n_chunk
                this_y_origin = y_chunk - n_chunk
                
                this_x_size = int(this_x_chunk - this_x_origin)
                this_y_size = int(this_y_chunk - this_y_origin)
                
                if this_x_size <= 0 or this_y_size <= 0:
                    break
                    
                yield (this_x_origin, this_y_origin, this_x_size, this_y_size)

                if y_chunk > ds_config['ny']:
                    break
                y_chunk += step
                    
            if x_chunk > ds_config['nx']:
                break
            x_chunk += step

    
def gdal_chunks(src_gdal, n_chunk, band=1):
    """Split `src_gdal` GDAL file into chunks."""
    
    o_chunks = []
    
    with gdal_datasource(src_gdal) as src_ds:
        if src_ds is None:
            return o_chunks
            
        ds_config = gdal_infos(src_ds)
        src_band = src_ds.GetRasterBand(band)
        gt = ds_config['geoT']
        
        gt_l = list(gt)
        gt_l[0] = gt_l[0] - (gt_l[1] / 2)
        gt_l[3] = gt_l[3] - (gt_l[5] / 2)
        
        c_n = 0
        for srcwin in gdal_yield_srcwin(src_gdal, n_chunk=n_chunk, step=n_chunk):
            this_geo_x_origin, this_geo_y_origin = utils._pixel2geo(
                srcwin[0], srcwin[1], tuple(gt_l)
            )
            
            dst_gt = [
                this_geo_x_origin, float(gt[1]), 0.0,
                this_geo_y_origin, 0.0, float(gt[5])
            ]
            
            band_data = src_band.ReadAsArray(*srcwin)
            
            ## Check if chunk has valid data (not just uniform background)
            if not np.all(band_data == band_data[0, 0]):
                dst_config = gdal_copy_infos(ds_config)
                dst_config['nx'] = srcwin[2]
                dst_config['ny'] = srcwin[3]
                dst_config['geoT'] = dst_gt
                
                base_name = os.path.splitext(os.path.basename(src_gdal))[0]
                o_chunk = f'{base_name}_chnk{c_n}.tif'
                dst_fn = os.path.join(os.path.dirname(src_gdal), o_chunk)
                
                gdal_write(band_data, dst_fn, dst_config)
                o_chunks.append(dst_fn)
                c_n += 1                
    return o_chunks


def gdal_parse(src_ds, dump_nodata=False, srcwin=None, mask=None,
               dst_srs=None, verbose=False, z_region=None, step=1):
    """Parse data from GDAL dataset (first band) and yield XYZ points."""
    
    ln = 0
    band = src_ds.GetRasterBand(1)
    ds_config = gdal_infos(src_ds)
    
    ## Setup Coordinate Transformation
    src_srs_obj = osr.SpatialReference()
    try:
        src_srs_obj.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    except Exception:
        pass
    
    if ds_config['proj']:
        src_srs_obj.ImportFromWkt(ds_config['proj'])
    
    if not src_srs_obj.GetAuthorityCode(None):
        src_srs_obj.ImportFromEPSG(4326)

    dst_trans = None
    if dst_srs is not None:
        dst_srs_obj = osr.SpatialReference()
        dst_srs_obj.SetFromUserInput(dst_srs)
        try:
            dst_srs_obj.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        except Exception:
            pass
        dst_trans = osr.CoordinateTransformation(src_srs_obj, dst_srs_obj)

    ## Masking
    msk_band = None
    src_mask = None
    if mask is not None:
        src_mask = gdal.Open(mask)
        msk_band = src_mask.GetRasterBand(1)
        
    if srcwin is None:
        srcwin = (0, 0, ds_config['nx'], ds_config['ny'])
        
    ndv = band.GetNoDataValue()
    
    gt = ds_config['geoT']

    ## Iterate rows
    for y_i in range(0, srcwin[3], step):
        y = srcwin[1] + y_i
        ## Read one row
        band_data = band.ReadAsArray(srcwin[0], y, srcwin[2], 1)
        if band_data is None:
            continue
            
        band_data = band_data.flatten()
        
        ## Apply Z filters
        if z_region is not None:
            if z_region[0] is not None:
                band_data[band_data < z_region[0]] = ndv if ndv else -9999
            if z_region[1] is not None:
                band_data[band_data > z_region[1]] = ndv if ndv else -9999
                
        if msk_band is not None:
            msk_data = msk_band.ReadAsArray(srcwin[0], y, srcwin[2], 1).flatten()
            band_data[msk_data == 0] = ndv if ndv else -9999
            
        for x_i in range(0, srcwin[2], step):
            z = band_data[x_i]
            
            is_valid = True
            if not dump_nodata:
                if ndv is not None and z == ndv:
                    is_valid = False
                elif np.isnan(z):
                    is_valid = False
            
            if is_valid:
                ln += 1
                x = srcwin[0] + x_i
                geo_x, geo_y = utils._pixel2geo(x, y, gt)
                
                ## Transform point if needed
                if dst_trans:
                    res = dst_trans.TransformPoint(geo_x, geo_y, float(z))
                    yield xyzfun.XYZPoint(x=res[0], y=res[1], z=res[2])
                else:
                    yield xyzfun.XYZPoint(x=geo_x, y=geo_y, z=z)

    if src_mask:
        src_mask = None
        
    if verbose:
        utils.echo_msg(f'Parsed {ln} data records from {src_ds.GetDescription()}')


def sample_warp(
        src_dem, dst_dem, x_sample_inc, y_sample_inc,
        src_srs=None, dst_srs=None, src_region=None, sample_alg='bilinear',
        ndv=-9999, tap=False, size=False, co=None,
        ot=gdal.GDT_Float32, coordinateOperation=None, verbose=False
):
    """Wrapper around gdal.Warp."""
    
    if co is None:
        co = ["COMPRESS=DEFLATE", "TILED=YES"]

    xcount = ycount = 0
    out_region = None
    
    if src_region is not None:
        out_region = [src_region.xmin, src_region.ymin,
                      src_region.xmax, src_region.ymax]        
        
        ## Transform region if SRS differs
        if src_srs is not None and dst_srs is not None:
            trans_region = src_region.copy()
            trans_region.src_srs = dst_srs
            trans_region.warp(src_srs)
        else:
            trans_region = src_region

        if size:
            if x_sample_inc is None and y_sample_inc is None:
                src_infos = gdal_infos(src_dem)
                ## Invert y_inc for geo_transform math
                xcount, ycount, _ = trans_region.geo_transform(
                    x_inc=src_infos['geoT'][1], y_inc=src_infos['geoT'][5] * -1, node='grid'
                )

            if xcount is None or ycount is None:
                xcount, ycount, _ = src_region.geo_transform(
                    x_inc=x_sample_inc, y_inc=y_sample_inc, node='grid'
                )
                x_sample_inc = y_sample_inc = None
        
    if dst_dem is not None:
        dirname = os.path.dirname(dst_dem)
        if dirname and not os.path.exists(dirname):
            os.makedirs(dirname)

    pbar_update = None
    if verbose:
        desc = f'Warping {src_dem}'
        pbar = utils.ccp(desc=desc, total=100, leave=verbose)
        ## Callback wrapper for GDAL
        def pbar_update(pct, msg, data):
            pbar.update((pct * 100) - pbar.n)
            return 1

    try:
        dst_ds = gdal.Warp(
            '' if dst_dem is None else dst_dem,
            src_dem,
            format='MEM' if dst_dem is None else 'GTiff',
            xRes=x_sample_inc,
            yRes=y_sample_inc,
            targetAlignedPixels=tap,
            width=int(xcount),
            height=int(ycount),
            dstNodata=ndv,
            outputBounds=out_region,
            outputBoundsSRS=dst_srs if out_region is not None else None,
            resampleAlg=sample_alg,
            errorThreshold=0,
            creationOptions=co,
            srcSRS=src_srs,
            dstSRS=dst_srs,
            coordinateOperation=coordinateOperation,
            outputType=ot,
            callback=pbar_update
        )
    finally:
        if verbose:
            pbar.close()
        
    if dst_dem is None:
        return dst_ds, 0
    else:
        dst_ds = None
        return dst_dem, 0


### End
