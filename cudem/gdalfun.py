### gdalfun.py - OSGEO functions
##
## Copyright (c) 2010 - 2024 Regents of the University of Colorado
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
## Functions, etc. for common gdal/ogr/osr usage.
##
### Code:

import os
import sys
import shutil
import math
from tqdm import tqdm
from tqdm import trange

from osgeo import gdal
from osgeo import osr
from osgeo import ogr
from pyproj import CRS

import numpy as np

from cudem import utils
from cudem import regions
from cudem import xyzfun

gc = utils.config_check()
gdal.DontUseExceptions()
ogr.DontUseExceptions()
osr.DontUseExceptions()
gdal.SetConfigOption('CPL_LOG', 'NUL' if gc['platform'] == 'win32' else '/dev/null') 

## OSR/WKT/proj
def split_srs(srs, as_epsg = False):
    """split an SRS into the horizontal and vertical elements.

    -----------
    Parameters:
    srs (str): an srs for input into SetFromUserInput
    as_epsg (bool): try to output EPSG codes if possible

    --------
    Returns:
    list: [horz_srs, vert_srs]
    """

    vert_epsg = None
    vert_wkt = None
    if srs is None:
        return(None, None)

    if np.any(['geoid' in x for x in srs.split('+')]):
        srs = '+'.join(srs.split('+')[:-1])
    
    if srs.split(':')[0] == 'ESRI':
        esri_split = srs.split('+')
        if len(esri_split) > 1:
            vert_epsg = srs.split('+')[1]
            vert_srs = osr.SpatialReference()
            vert_srs.SetFromUserInput('EPSG:{}'.format(vert_epsg))
            vert_wkt = vert_srs.ExportToWkt()
            
        srs = esri_split[0]
        
    src_srs = osr.SpatialReference()
    src_srs.SetFromUserInput(srs)
    #try:
    srs_wkt = src_srs.ExportToWkt()
    wkt_CRS = CRS.from_wkt(srs_wkt)
    #except:
    #    return(None, None)
    
    if wkt_CRS.is_compound:
        horz = wkt_CRS.sub_crs_list[0]
        horz_epsg = horz.to_epsg()
        horz_wkt = horz.to_wkt()
        vert = wkt_CRS.sub_crs_list[1]
        vert_epsg = vert.to_epsg()
        vert_wkt = vert.to_wkt()
    else:
        horz = wkt_CRS
        horz_epsg = horz.to_epsg()
        horz_wkt = horz.to_wkt()

    if as_epsg:
        return(horz_epsg if horz_epsg is not None else horz_wkt, vert_epsg if vert_epsg is not None else vert_wkt)
    else:
        return(horz_wkt, vert_epsg)

def combine_epsgs(src_horz, src_vert, name='Combined'):
    """combine src_horz and src_vert into a CompoundCS"""
    
    if src_horz is None or src_vert is None:
        return(None)
    
    horz_srs = osr.SpatialReference()
    horz_srs.SetFromUserInput(src_horz)
    vert_srs = osr.SpatialReference()
    vert_srs.SetFromUserInput('epsg:{}'.format(src_vert))
    src_srs = osr.SpatialReference()
    src_srs.SetCompoundCS('{}'.format(name, src_horz, src_vert), horz_srs, vert_srs)
    return(src_srs.ExportToWkt())
    
def wkt2geom(wkt):
    """transform a wkt to an ogr geometry

    -----------
    Parameters:
    wkt (wkt): a wkt geometry

    --------
    Returns:
    ogr-geom: the ogr geometry
    """
    
    return(ogr.CreateGeometryFromWkt(wkt))

def osr_wkt(src_srs, esri=False):
    """convert a src_srs to wkt"""
    
    try:
        sr = osr.SpatialReference()
        sr.SetFromUserInput(src_srs)
        if esri:
            sr.MorphToESRI()
            
        return(sr.ExportToWkt())
    except:
        return(None)

def osr_prj_file(dst_fn, src_srs):
    """generate a .prj file given a src_srs

    -----------
    Parameters:
    dst_fn (str): the output filename
    src_srs (srs): the output srs
    """
    
    with open(dst_fn, 'w') as out:
        out.write(osr_wkt(src_srs, True))

    if os.path.exists(dst_fn):
        return(0)
    else:
        return(-1)

def epsg_from_input(in_srs):
    """get the epsg(s) from srs suitable as input to SetFromUserInput

    -----------
    Parameters:
    in_srs (str): an srs as a string

    -----------
    Returns:
    list: [horz_epsg, vert_epsg]
    """

    src_vert = None
    if np.any(['geoid' in x for x in in_srs.split('+')]):
        in_srs = '+'.join(in_srs.split('+')[:-1])
        
    if in_srs.split(':')[0] == 'ESRI':
        esri_split = in_srs.split('+')
        if len(esri_split) > 1:
            vert_epsg = in_srs.split('+')[1]
            vert_srs = osr.SpatialReference()
            vert_srs.SetFromUserInput('EPSG:{}'.format(vert_epsg))
            vert_wkt = vert_srs.ExportToWkt()

            src_vert = vert_epsg
            
        in_srs = esri_split[0]
    
    src_srs = osr.SpatialReference()
    src_srs.SetFromUserInput(in_srs)

    ## HORZ
    if src_srs.IsGeographic() == 1:
        cstype = 'GEOGCS'
    else:
        cstype = 'PROJCS'

    src_srs.AutoIdentifyEPSG()
    an = src_srs.GetAuthorityName(cstype)
    src_horz = src_srs.GetAuthorityCode(cstype)

    if src_horz is None:
        src_horz = src_srs.ExportToProj4()
        if src_horz == '':
            src_horz = src_srs.ExportToWkt()
    else:
        src_horz = '{}:{}'.format(an, src_horz)
        
    ## VERT
    if src_srs.IsVertical() == 1:
        csvtype = 'VERT_CS'
        src_vert = src_srs.GetAuthorityCode(csvtype)
    #else:
    #    src_vert = None

    return(src_horz, src_vert)

def osr_parse_srs(src_srs, return_vertcs = True):
    """parse an OSR SRS object and return a proj string"""
    
    if src_srs is not None:
        if src_srs.IsLocal() == 1:
            return(src_srs.ExportToWkt())
        if src_srs.IsGeographic() == 1:
            cstype = 'GEOGCS'
        else:
            cstype = 'PROJCS'
            
        src_srs.AutoIdentifyEPSG()
        an = src_srs.GetAuthorityName(cstype)
        ac = src_srs.GetAuthorityCode(cstype)

        #if return_vertcs:
        if src_srs.IsVertical() == 1:
            csvtype = 'VERT_CS'
            vn = src_srs.GetAuthorityName(csvtype)
            vc = src_srs.GetAuthorityCode(csvtype)
        else:
            csvtype = vc = vn = None

        if an is not None and ac is not None:
            if vn is not None and vc is not None:
                return('{}:{}+{}'.format(an, ac, vc))
            else:
                return('{}:{}'.format(an, ac))
        else:
            dst_srs = src_srs.ExportToProj4()
            if dst_srs:
                    return(dst_srs)
            else:
                return(None)
    else:
        return(None)

## OGR
def ogr_or_gdal(osgeo_fn):
    '''ogr file is 1, gdal file is 2, neither is -1'''
    
    try:
        src_ds = ogr.Open(osgeo_fn)
        if src_ds is not None:
            src_ds = None
            return(1)
        else:
            src_ds = gdal.Open(osgeo_fn)
            if src_ds is not None:
                src_ds = None
                return(2)
            else:
                return(-1)
    except:
        utils.echo_error_msg('{} does not appear to be an osgeo (gdal/ogr) file')
        return(-1)

def ogr_fext(src_drv_name):
    """find the common file extention given a OGR driver name
    older versions of gdal can't do this, so fallback to known standards.

    -----------
    Parameters:
    src_drv_name (str): a source OGR driver name

    -----------
    Returns:
    list: a list of known file extentions or None
    """
    
    fexts = None
    try:
        drv = ogr.GetDriverByName(src_drv_name)
        fexts = drv.GetMetadataItem(gdal.DMD_EXTENSIONS)
        if fexts is not None:
            return(fexts.split()[0])
        else:
            return(None)
    except:
        if src_drv_name.lower() == 'gtiff': fext = 'tif'
        elif src_drv_name == 'HFA': fext = 'img'
        elif src_drv_name == 'GMT': fext = 'grd'
        elif src_drv_name.lower() == 'netcdf': fext = 'nc'
        else: fext = 'gdal'
        
        return(fext)
    
def ogr_get_srs(src_ogr):
    """get the srs (as wkt) from an ogr file"""
    
    src_ds = ogr.Open(src_ogr, 0)
    src_srs = src_ds.GetLayer().GetSpatialRef()
    src_ds = None
    if src_srs is not None:
        return(src_srs.ExportToWkt())
    else:
        return(None)
    
def ogr_clip(src_ogr_fn, dst_region = None, layer = None, overwrite = False, verbose = True):
    """clip an ogr file to `dst_region`"""
    
    dst_ogr_bn = '.'.join(src_ogr_fn.split('.')[:-1])
    dst_ogr_fn = '{}_{}.gpkg'.format(dst_ogr_bn, dst_region.format('fn'))
    
    if not os.path.exists(dst_ogr_fn) or overwrite:
        utils.run_cmd('ogr2ogr -nlt PROMOTE_TO_MULTI {} {} -clipsrc {} {} '.format(
            dst_ogr_fn, src_ogr_fn, dst_region.format('te'), layer if layer is not None else ''
        ), verbose=verbose)
                
    return(dst_ogr_fn)

def ogr_clip2(src_ogr_fn, dst_region=None, layer=None, overwrite=False):
    """clip an ogr file to `dst_region`"""
    
    dst_ogr_bn = '.'.join(src_ogr_fn.split('.')[:-1])
    dst_ogr_fn = '{}_{}.gpkg'.format(dst_ogr_bn, dst_region.format('fn'))
    
    if not os.path.exists(dst_ogr_fn) or overwrite:
    
        src_ds = ogr.Open(src_ogr_fn)
        if layer is not None:
            src_layer = src_ds.GetLayer(layer)
        else:
            src_layer = src_ds.GetLayer()
            
        region_ogr = 'region_{}.shp'.format(dst_region.format('fn'))
        dst_region.export_as_ogr(region_ogr)
        region_ds = ogr.Open(region_ogr)
        region_layer = region_ds.GetLayer()

        driver = ogr.GetDriverByName('GPKG')
        dst_ds = driver.CreateDataSource(dst_ogr_fn)
        dst_layer = dst_ds.CreateLayer(layer if layer is not None else 'clipped', geom_type=ogr.wkbPolygon)

        ogr.Layer.Clip(src_layer, region_layer, dst_layer)

        src_ds = region_ds = dst_ds = None
        utils.remove_glob('{}.*'.format('.'.join(region_ogr.split('.')[:-1])))
        
    return(dst_ogr_fn)

def ogr_clip3(src_ogr, dst_ogr, clip_region=None, dn="ESRI Shapefile"):
    """clip an ogr file to `dst_region`"""
    
    driver = ogr.GetDriverByName(dn)
    ds = driver.Open(src_ogr, 0)
    layer = ds.GetLayer()

    clip_region.export_as_ogr('tmp_clip.{}'.format(ogr_fext(dn)))
    c_ds = driver.Open('tmp_clip.{}'.format(ogr_fext(dn)), 0)
    c_layer = c_ds.GetLayer()
    
    dst_ds = driver.CreateDataSource(dst_ogr)
    dst_layer = dst_ds.CreateLayer(
        dst_ogr.split('.')[0], geom_type=ogr.wkbMultiPolygon
    )

    layer.Clip(c_layer, dst_layer)
    ds = c_ds = dst_ds = None

## in cudem.fetches as `polygonize_osm_coastline`
def ogr_polygonize_line_to_region(
        src_ogr, dst_ogr, region=None, include_landmask=True,
        landmask_is_watermask=False, line_buffer=0.0000001
):
    """Polygonize an OSM coastline LineString to the given region

    if include_landmask is True, polygon(s) will be returned for the land areas 
    with a `watermask` valule of 0, otherwise, only the watermask polygon will
    be returned with a `watermask` value of 1.
    """

    # Open the input line layer
    line_ds = ogr.Open(src_ogr)
    line_layer = line_ds.GetLayer()
    region_geom = regions.Region().from_list(line_layer.GetExtent()).export_as_geom()
    ## todo: check if input region is larger than the line region, if so,
    ##       reduce the region to the size of the line region...
    if region is not None and region.valid_p():
        region_geom = region.export_as_geom()
    
    # Create the output layer
    driver = ogr.GetDriverByName("ESRI Shapefile")
    output_ds = driver.CreateDataSource(dst_ogr)
    output_layer = output_ds.CreateLayer("split_polygons", line_layer.GetSpatialRef(), ogr.wkbMultiPolygon)
    output_layer.CreateField(ogr.FieldDefn('watermask', ogr.OFTInteger))

    has_feature = False
    
    for line_layer in line_ds:
        line_type = line_layer.GetGeomType()
        if line_type == 2:
            has_feature = 1
            line_geometries = ogr_union_geom(line_layer, ogr.wkbMultiLineString if line_type == 2 else ogr.wkbMultiPolygon)##
            poly_line = line_geometries.Buffer(line_buffer)
            split_geoms = region_geom.Difference(poly_line)
            for split_geom in split_geoms:
                ss = []
                for line_geometry in line_geometries:
                    if split_geom.Intersects(line_geometry.Buffer(line_buffer)):
                        point_count = line_geometry.GetPointCount()
                        for point_n in range(0, point_count-1):
                            x_beg = line_geometry.GetX(point_n)
                            y_beg = line_geometry.GetY(point_n)
                            x_end = line_geometry.GetX(point_n+1)
                            y_end = line_geometry.GetY(point_n+1)

                            poly_center = split_geom.Centroid()
                            poly_center_x = poly_center.GetX(0)
                            poly_center_y = poly_center.GetY(0)

                            s = (x_end - x_beg)*(poly_center_y - y_beg) > (y_end - y_beg)*(poly_center_x - x_beg)
                            ss.append(s)

                if all(ss):
                    s = True
                elif not any(ss):
                    s = False
                else:
                    if np.count_nonzero(ss) > len(ss) / 2:
                        s = True
                    else:
                        s = False

                out_feature = ogr.Feature(output_layer.GetLayerDefn())
                out_feature.SetGeometry(split_geom)

                if landmask_is_watermask:
                    s = False if s else True
                
                if s == 0:
                    out_feature.SetField('watermask', 1)
                    output_layer.CreateFeature(out_feature)

                if include_landmask:
                    if s == 1:
                        out_feature.SetField('watermask', 0)
                        output_layer.CreateFeature(out_feature)
                
        if line_type == 6:
            has_feature = 1
            for line_feature in line_layer:
                line_geometry = line_feature.geometry()
                line_geometry = ogr.ForceTo(line_geometry, ogr.wkbLinearRing)
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

    line_ds = output_ds = None
    return(dst_ogr)
    
def ogr_mask_union(src_layer, src_field=None, dst_defn=None):
    """`union` a `src_layer`'s features based on `src_field` where
    `src_field` holds a value of 0 or 1. optionally, specify
    an output layer defn for the unioned feature.

    -----------
    Returns:
    the output feature class
    """
    
    if dst_defn is None:
        dst_defn = src_layer.GetLayerDefn()
        
    multi = ogr.Geometry(ogr.wkbMultiPolygon)
    if src_field is not None:
        src_layer.SetAttributeFilter("{} = 1".format(src_field))
        
    feats = len(src_layer)
    
    if feats > 0:
        with tqdm(total=len(src_layer), desc='unioning {} features...'.format(feats)) as pbar:
            for n, f in enumerate(src_layer):
                pbar.update()
                f_geom = f.geometry()
                f_geom_valid = f_geom.Buffer(0)
                multi.AddGeometry(f_geom_valid)
            
    utils.echo_msg('setting geometry to unioned feature...')
    out_feat = ogr.Feature(dst_defn)
    out_feat.SetGeometry(multi)
    union = multi = None
    
    return(out_feat)

def ogr_union_geom(src_layer, geom_type = ogr.wkbMultiPolygon, verbose = True):
        
    multi = ogr.Geometry(geom_type)
    feats = src_layer.GetFeatureCount()#len(src_layer)
    [multi.AddGeometry(f.geometry()) for f in src_layer]
    if verbose:
        utils.echo_msg('unioned {} features'.format(feats))
        
    # if feats > 0:
    #     with tqdm(total=len(src_layer), desc='unioning {} features...'.format(feats)) as pbar:
    #         for n, f in enumerate(src_layer):
    #             pbar.update()
    #             f_geom = f.geometry()
    #             f_geom_valid = f_geom
    #             multi.AddGeometry(f_geom_valid)

    return(multi)

def ogr_empty_p(src_ogr, dn='ESRI Shapefile'):
    """check if the OGR file is empty or not"""
    
    driver = ogr.GetDriverByName(dn)
    ds = driver.Open(src_ogr, 0)

    if ds is not None:
        layer = ds.GetLayer()
        fc = layer.GetFeatureCount()
        if fc == 0:
            return(True)
        else: return(False)
        
    else: return(True)

def ogr_polygonize(
        src_ds, dst_srs = 'epsg:4326', ogr_format = 'ESRI Shapefile', band = 1, verbose = True
):
    """polygonze a multi-band raster. 
    for use in generating spatial metadata from a waffles mask.

    -----------
    Parameters:
    src_ds (GDALDataset): the source multi-band raster as a gdal dataset object
    dst_srs (str): the output srs
    ogr_format (str): the output OGR format
    verbose (bool): be verbose

    -----------
    Returns:
    tuple: (dst_layer, ogr_format)
    """
    
    dst_layer = '{}_ply'.format(utils.fn_basename2(src_ds.GetDescription()))
    dst_vector = dst_layer + '.{}'.format(ogr_fext(ogr_format))
    utils.remove_glob('{}.*'.format(dst_layer))
    osr_prj_file('{}.prj'.format(dst_layer), gdal_infos(src_ds)['proj'])
    driver = ogr.GetDriverByName(ogr_format)
    ds = driver.CreateDataSource(dst_vector)
    if ds is not None: 
        layer = ds.CreateLayer(
            '{}'.format(dst_layer), None, ogr.wkbMultiPolygon
        )
        [layer.SetFeature(feature) for feature in layer]
    else:
        layer = None

    layer.CreateField(ogr.FieldDefn('DN', ogr.OFTInteger))
    defn = None

    this_band = src_ds.GetRasterBand(band)
    this_band_md = this_band.GetMetadata()
    b_infos = gdal_infos(src_ds, scan=True, band=band)
    field_names = [field.name for field in layer.schema]
    for k in this_band_md.keys():
        this_band_md[k.title()] = this_band_md.pop(k)

    for k in this_band_md.keys():
        if k[:9] not in field_names:
            layer.CreateField(ogr.FieldDefn(k[:9], ogr.OFTString))

    if ds is not None:
        if verbose:
            utils.echo_msg('polygonizing {} mask...'.format(this_band.GetDescription()))

        status = gdal.Polygonize(
            this_band,
            None,
            layer,
            layer.GetLayerDefn().GetFieldIndex('DN'),
            [],
            #callback = gdal.TermProgress if verbose else None
            callback = None
        )
        
        if verbose:
            utils.echo_msg('polygonized {}'.format(this_band.GetDescription()))

        for feat in layer:
            if feat.GetField('DN') == 0:
                layer.DeleteFeature(feat.GetFID())
            
    ds = layer = None
    return(dst_layer, ogr_format)

def ogr2gdal_mask(
        mask_fn, region = None, x_inc = None, y_inc = None, dst_srs = 'epsg:4326',
        invert = True, verbose = True, temp_dir = utils.cudem_cache()
):
    dst_fn = utils.make_temp_fn('{}.tif'.format(mask_fn, temp_dir=temp_dir))
    if os.path.exists(dst_fn):
        return(dst_fn)
    else:
        # if os.path.isdir(self.mask):
        #     dst_layer = os.path.basename('/'.join(self.mask.split('/')[:-1])).split('.')[0]
        #msk_region = self.region if self.region is not None else regions.Region().from_list(self.infos.minmax)

        if region is not None and x_inc is not None and y_inc is not None:
            msk_region = region.copy()
            xcount, ycount, dst_gt = msk_region.geo_transform(
                x_inc=x_inc, y_inc=y_inc, node='grid'
            )
            if xcount <= 0 or ycount <=0:
                utils.echo_error_msg(
                    'could not create grid of {}x{} cells with {}/{} increments on region: {}'.format(
                        xcount, ycount, x_inc, y_inc, region
                    )
                )
                sys.exit(-1)

            ds_config = gdal_set_infos(xcount, ycount, xcount * ycount, dst_gt, dst_srs, gdal.GDT_Float32, 0, 'GTiff', {}, 1)
            gdal_nan(ds_config, dst_fn, nodata=0)
            clip_layer = utils.fn_basename2(os.path.basename(mask_fn))
            mask_fn = ogr_clip(mask_fn, dst_region=region, layer=clip_layer, verbose=verbose)
            
            gr_cmd = 'gdal_rasterize -burn {} -l {} {} {}{}'\
                .format(1, clip_layer, mask_fn, dst_fn, ' -i' if invert else '')
            #utils.echo_msg(gr_cmd)
            out, status = utils.run_cmd(gr_cmd, verbose=verbose)
            return(dst_fn)

        #else:
        #    self.mask = None

## GDAL
class gdal_datasource:
    """
    use as:
    with gdal_datasource('input.tif') as src_ds:
        src_ds.do_things()


    where the input can be either a path-name to a gdal supported
    file or a gdal.Dataset object.

    this is for when you don't like opening and closing gdal datasurces all the time.
    """
    
    def __init__(self, src_gdal = None, update = False):
        self.src_gdal = src_gdal
        self.update = update
        self.src_ds = None

    ## todo: create new datasource if src_gdal is path, but doesn't exists...
    def __enter__(self):
        if isinstance(self.src_gdal, gdal.Dataset):
            try:
                self.src_ds = self.src_gdal
            except:
                self.src_ds = None

        elif utils.str_or(self.src_gdal) is not None \
             and (os.path.exists(self.src_gdal) \
                  or utils.fn_url_p(self.src_gdal) \
                  or len(self.src_gdal.split(':')) > 1):
            try:
                if self.update:
                    self.src_ds = gdal.OpenEx(self.src_gdal, 1, open_options=['IGNORE_COG_LAYOUT_BREAK=YES'])
                else:
                    self.src_ds = gdal.Open(self.src_gdal, 0)
                    
            except:
                self.src_ds = None

        ## remove overviews in update mode
        if self.src_ds is not None:
            if self.update:
                self.src_ds.BuildOverviews('AVERAGE', [])
            
        return(self.src_ds)

    def __exit__(self, exc_type, exc_value, exc_traceback):
        ## don't close src_ds if incoming was already open
        #print(exc_type, exc_value, exc_traceback)
        if not isinstance(self.src_gdal, gdal.Dataset):
            if self.update:
                try:
                    self.src_ds.FlushCache()
                except:
                    pass
                
            self.src_ds = None

def gdal_get_srs(src_gdal):
    """get the srs (as wkt) from a gdal file"""

    with gdal_datasource(src_gdal) as src_ds:
        if src_ds is not None:
            src_srs = src_ds.GetSpatialRef()
            src_ds = None
            if src_srs is not None:
                return(src_srs.ExportToWkt())
            else:
                return(None)
        else:
            return(None)
            
def gdal_multi_mask2single_mask(src_gdal):
    """transformat a mult-banded mask into a single-banded mask"""
    
    gdt = gdal.GDT_Int32
    out_mask = None    
    with gdal_datasource(src_gdal) as src_ds:
        if src_ds is not None:
            ds_config = gdal_infos(src_ds)
            out_mask = '{}_single.{}'.format(utils.fn_basename2(src_ds.GetDescription()), gdal_fext(ds_config['fmt']))
            driver = gdal.GetDriverByName(ds_config['fmt'])
            m_ds = driver.Create(out_mask, ds_config['nx'], ds_config['ny'], 1, gdt)
            m_ds.SetGeoTransform(ds_config['geoT'])
            m_band = m_ds.GetRasterBand(1)
            m_band.SetNoDataValue(ds_config['ndv'])

            for b in range(1, src_ds.RasterCount+1):
                this_band = src_ds.GetRasterBand(b)
                this_arr = this_band.ReadAsArray()
                m_band.WriteArray(this_arr)

            m_ds = None

    return(out_mask)

def gdal_fext(src_drv_name):
    """find the common file extention given a GDAL driver name
    older versions of gdal can't do this, so fallback to known standards.

    -----------
    Parameters:
    src_drv_name (str): a source GDAL driver name

    -----------
    Returns:
    list: a list of known file extentions or None
    """
    
    fexts = None
    try:
        drv = gdal.GetDriverByName(src_drv_name)
        if drv.GetMetadataItem(gdal.DCAP_RASTER):
            fexts = drv.GetMetadataItem(gdal.DMD_EXTENSIONS)
            
        if fexts is not None:
            return(fexts.split()[0])
        else:
            return(None)
    except:
        if src_drv_name.lower() == 'gtiff': fext = 'tif'
        elif src_drv_name == 'HFA': fext = 'img'
        elif src_drv_name == 'GMT': fext = 'grd'
        elif src_drv_name.lower() == 'netcdf': fext = 'nc'
        else: fext = 'gdal'
        
        return(fext)

def gdal_set_infos(nx, ny, nb, geoT, proj, dt, ndv, fmt, md, rc):
    """set a datasource config dictionary

    -----------
    Parameters:
    nx: x size
    ny: y size
    nb: cell count
    geoT: geotransform
    proj: projection
    dt: data type
    ndv: nodata value
    fmt: format
    metadata: metadata
    raster_count: number of bands

    -----------
    Returns:
    gdal_config dict
    """
    
    return({'nx': nx, 'ny': ny, 'nb': nb, 'geoT': geoT, 'proj': proj, 'dt': dt,
            'ndv': ndv, 'fmt': fmt, 'metadata': md, 'raster_count': rc})
    
def gdal_infos(src_gdal, region = None, scan = False, band = 1):
    """gather and return some info about a src_gdal file"""
    
    ds_config = {}
    with gdal_datasource(src_gdal) as src_ds:
        if src_ds is not None:
            gt = src_ds.GetGeoTransform()
            if region is not None:
                srcwin = region.srcwin(gt, src_ds.RasterXSize, src_ds.RasterYSize)
            else:
                srcwin = (0, 0, src_ds.RasterXSize, src_ds.RasterYSize)

            dst_gt = (gt[0] + (srcwin[0] * gt[1]), gt[1], 0., gt[3] + (srcwin[1] * gt[5]), 0., gt[5])
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
                    #'zr': None,
                }
                if ds_config['ndv'] is None:
                    ds_config['ndv'] = -9999

                if scan:
                    src_arr = src_band.ReadAsArray(srcwin[0], srcwin[1], srcwin[2], srcwin[3])
                    try:
                        src_arr[src_arr == ds_config['ndv']] = np.nan
                        if not np.all(np.isnan(src_arr)):
                            ds_config['zr'] = src_band.ComputeRasterMinMax()
                        else:
                            #utils.echo_warning_msg('{} is all nan'.format(src_ds.GetDescription()))
                            ds_config['zr'] = [np.nan, np.nan]
                    except:
                        if not np.all(src_arr == ds_config['ndv']):
                            ds_config['zr'] = src_band.ComputeRasterMinMax()
                        else:
                            #utils.echo_warning_msg('{} is all nan'.format(src_ds.GetDescription()))
                            ds_config['zr'] = [np.nan, np.nan]

                    src_arr = src_band = None
            else:
                utils.echo_warning_msg('invalid band {} for data source {}'.format(band, src_gdal))

    return(ds_config)

def gdal_copy_infos(src_config):
    """copy src_config

    -----------
    Returns:
    copied src_config dict.
    """
    
    dst_config = {}
    for dsc in src_config.keys():
        dst_config[dsc] = src_config[dsc]
        
    return(dst_config)
        
def gdal_set_srs(src_gdal, src_srs = 'epsg:4326', verbose = True):
    """set the src_gdal srs"""
    
    status = None
    if '+geoid' in src_srs:
        src_srs = '+'.join(src_srs.split('+')[:-1])
        
    with gdal_datasource(src_gdal, update=True) as src_ds:    
        if src_ds is not None and src_srs is not None:
            try:
                src_ds.SetProjection(osr_wkt(src_srs))
                status = 0
            except Exception as e:
                if verbose:
                    utils.echo_warning_msg('could not set projection {}'.format(src_srs))
                status = None
        else:
            status = None
            
    return(status)

def gdal_get_ndv(src_gdal, band = 1):
    """get the src_gdal nodata value"""

    with gdal_datasource(src_gdal) as src_ds:    
        if src_ds is not None:
            src_band = src_ds.GetRasterBand(band)
        else:
            return(None)
        
    return(src_band.GetNoDataValue())

def gdal_set_ndv(src_gdal, ndv = -9999, convert_array = False, verbose = True):
    """set the nodata value of gdal datasource"""

    with gdal_datasource(src_gdal, update=True) as src_ds:
        if src_ds is not None:
            ds_config = gdal_infos(src_ds)
            curr_nodata = ds_config['ndv']
            if verbose:
                utils.echo_msg('setting nodata value from {} to {}'.format(curr_nodata, ndv))
                               
            for band in range(1, src_ds.RasterCount+1):
                this_band = src_ds.GetRasterBand(band)
                this_band.DeleteNoDataValue()

            for band in range(1, src_ds.RasterCount+1):
                this_band = src_ds.GetRasterBand(band)
                this_band.SetNoDataValue(ndv)

            if convert_array:
                for band in range(1, src_ds.RasterCount+1):

                    this_band = src_ds.GetRasterBand(band)
                    arr = this_band.ReadAsArray()
                    if np.isnan(curr_nodata):
                        arr[np.isnan(arr)] = ndv
                    else:
                        arr[arr == curr_nodata] = ndv

                    this_band.WriteArray(arr)            
        else:
            return(None)

    return(0)

def gdal_has_ndv(src_gdal, band = 1):
    """check if src_gdal file has a set nodata value"""
    
    with gdal_datasource(src_gdal) as src_ds:        
        if src_ds is not None:
            ds_config = gdal_infos(src_ds)
            ds_arr = ds.GetRasterBand(band).ReadAsArray()
            ds_arr[ds_arr == ds_config['ndv']] = np.nan

            if np.any(np.isnan(ds_arr)):
                return(True)
            
        return(False)
    
def gdal_mem_ds(ds_config, name = 'MEM', bands = 1, src_srs = None, co = ['COMPRESS=DEFLATE', 'TILED=YES']):
    """Create temporary gdal mem dataset"""
        
    mem_driver = gdal.GetDriverByName('MEM')
    mem_ds = mem_driver.Create(name, ds_config['nx'], ds_config['ny'], bands, ds_config['dt'], options=co)
    if mem_ds is not None:
        mem_ds.SetGeoTransform(ds_config['geoT'])
        if src_srs is None:
            if ds_config['proj'] is not None:
                mem_ds.SetProjection(ds_config['proj'])
        else:
            mem_ds.SetProjection(src_srs)

        for band in range(1, bands+1):
            mem_band = mem_ds.GetRasterBand(band)
            mem_band.SetNoDataValue(ds_config['ndv'])
        
    return(mem_ds)

def gdal_nan(ds_config, outfile, nodata = None):
    '''create a nodata grid'''

    if nodata is None:
        nodata = np.nan
        
    nullArray = np.zeros( (ds_config['ny'], ds_config['nx']) )
    if nodata != 0:
        nullArray[nullArray==0]=nodata

    ds_config['ndv'] = nodata        
    gdal_write(nullArray, outfile, ds_config)

def gdal_get_node(src_gdal, node='pixel'):
    with gdal_datasource(src_gdal) as src_ds:
        ds_config = gdal_infos(src_ds)
        if 'metadata' in ds_config.keys():
            mt = ds_config['metadata']
            if 'AREA_OR_POINT' in mt.keys():
                node = mt['AREA_OR_POINT']
            elif 'NC_GLOBAL#node_offset' in mt.keys():
                node = mt['NC_GLOBAL#node_offset']

            elif 'tos#node_offset' in mt.keys():
                node = mt['tos#node_offset']

            elif 'node_offset' in mt.keys():
                node = mt['node_offset']

            if node.upper() == 'AREA' or node == '1':
                node = 'pixel'
            elif node.upper() == 'POINT' or node == '0':
                node = 'grid'

    return(node)
    
def gdal_extract_band(src_gdal, dst_gdal, band = 1, exclude = [], srcwin = None, inverse = False):
    """extract a band from the src_gdal file and put it into dst_gdal"""
    
    band = utils.int_or(band, 1)
    #srcwin = None
    with gdal_datasource(src_gdal) as src_ds:
        try:
            ds_config = gdal_infos(src_ds)
            ds_band = src_ds.GetRasterBand(band)
            if srcwin is not None:
                ds_array = ds_band.ReadAsArray(*srcwin)
                xy_origin = utils._pixel2geo(srcwin[0], srcwin[1], ds_config['geoT'], node='grid')
                ds_config['geoT'] = (xy_origin[0], ds_config['geoT'][1], ds_config['geoT'][2], xy_origin[1], ds_config['geoT'][4], ds_config['geoT'][5])
                ds_config['nx'] = srcwin[2]
                ds_config['ny'] = srcwin[3]
                ds_config['nb'] = srcwin[2] * srcwin[3]
            else:
                ds_array = ds_band.ReadAsArray()
                                
            if ds_array is None:
                utils.echo_error_msg('could not read data array from datasource {}'.format(src_gdal))
                return(None, -1)
            
        except Exception as e:
            utils.echo_error_msg('could not read datasource {}, {}'.format(src_gdal, e))
            return(None, -1)

    if ds_config['ndv'] is None:
        ds_config['ndv'] = -9999
    
    if exclude:
        for key in exclude:
            ds_array[ds_array == key] = ds_config['ndv']

    if inverse:
        ds_array = 1/ds_array

    return(gdal_write(ds_array, dst_gdal, ds_config))

def gdal_get_array(src_gdal, band = 1):
    """get the associated array from the src_gdal file"""
    
    with gdal_datasource(src_gdal) as src_ds:        
        if src_ds is not None:
            infos = gdal_infos(src_ds)
            src_band = src_ds.GetRasterBand(band)
            src_array = src_band.ReadAsArray()
            src_offset = src_band.GetOffset()
            src_scale = src_band.GetScale()
        else:
            utils.echo_error_msg('could not open {}'.format(src_gdal))
            return(None, None)
            
    if src_offset is not None or src_scale is not None:
        src_array = np.ndarray.astype(src_array, dtype=np.float32)

        if src_offset is not None:
            src_array += src_offset

        if src_scale is not None:
            src_array *= src_scale

    return(src_array, infos)

def gdal_cut(src_gdal, src_region, dst_gdal, node='pixel', verbose=True, co=["COMPRESS=DEFLATE", "TILED=YES"]):
    """cut src_ds datasource to src_region and output dst_gdal file

    -----------
    Returns:    
    list: [output-dem, status-code]
    """
    
    status = -1
    with gdal_datasource(src_gdal) as src_ds:
        if src_ds is not None:
            ds_config = gdal_infos(src_ds)
            gt = ds_config['geoT']
            srcwin = src_region.srcwin(gt, src_ds.RasterXSize, src_ds.RasterYSize, node=node)
            dst_gt = (gt[0] + (srcwin[0] * gt[1]), gt[1], 0., gt[3] + (srcwin[1] * gt[5]), 0., gt[5])
            out_ds_config = gdal_set_infos(srcwin[2], srcwin[3], srcwin[2] * srcwin[3], dst_gt,
                                           ds_config['proj'], ds_config['dt'], ds_config['ndv'],
                                           ds_config['fmt'], ds_config['metadata'], ds_config['raster_count'])

            in_bands = src_ds.RasterCount
            mem_ds = gdal_mem_ds(out_ds_config, bands=in_bands)
            if mem_ds is not None:
                for band in range(1, in_bands+1):
                    this_band = mem_ds.GetRasterBand(band)
                    #this_band_md = this_band.GetMetadata()
                    
                    that_band = src_ds.GetRasterBand(band)
                    #that_band_md = that_band.GetMetadata()
                    
                    this_band.SetDescription(that_band.GetDescription())
                    this_band.SetMetadata(that_band.GetMetadata())
                    
                    # for key in that_band_md.keys():
                    #     this_band_md[key] = that_band_md[key]

                    this_band.WriteArray(src_ds.GetRasterBand(band).ReadAsArray(*srcwin))
                    mem_ds.FlushCache()

                dst_ds = gdal.GetDriverByName(ds_config['fmt']).CreateCopy(dst_gdal, mem_ds, 0, options=co)
                status = 0

    return(dst_gdal, status)

## doesn't work with multipolygons and -i
def gdal_clip(src_gdal, dst_gdal, src_ply = None, invert = False, verbose = True, cache_dir = None):
    """clip dem to polygon `src_ply`, optionally invert the clip.

    -----------
    Returns:
    list: [gdal_raserize-output, gdal_rasterize-return-code]
    """

    gi = gdal_infos(src_gdal)
    g_region = regions.Region().from_geo_transform(geo_transform=gi['geoT'], x_count=gi['nx'], y_count=gi['ny'])
    tmp_ply = utils.make_temp_fn('tmp_clp_ply.shp', temp_dir=cache_dir)
    
    out, status = utils.run_cmd('ogr2ogr {} {} -clipsrc {} -nlt MULTIPOLYGON -skipfailures -makevalid'.format(tmp_ply, src_ply, g_region.format('ul_lr')), verbose=verbose)
    if gi is not None and src_ply is not None:
        #if invert:
        #    gr_cmd = 'gdalwarp -cutline {} -cl {} {} {}'.format(tmp_ply, os.path.basename(tmp_ply).split('.')[0], src_dem, dst_dem)
        #    out, status = utils.run_cmd(gr_cmd, verbose=verbose)
        #else:
        shutil.copyfile(src_gdal, dst_gdal)
        band_count = 1
        with gdal_datasource(src_gdal) as src_ds:
            if src_ds is not None:
                band_count = src_ds.RasterCount
                
        gr_cmd = 'gdal_rasterize -b {} -burn {} -l {} {} {}{}'\
            .format(band_count, gi['ndv'], os.path.basename(tmp_ply).split('.')[0], tmp_ply, dst_gdal, ' -i' if invert else '')
        out, status = utils.run_cmd(gr_cmd, verbose=verbose)
        utils.remove_glob(tmp_ply)
        #utils.remove_glob('{}*'.format(utils.fn_basename2(tmp_ply)))#'__tmp_clp_ply.*')
    else:
        utils.remove_glob(tmp_ply)
        return(None, -1)
    
    return(out, status)

def crop(src_gdal, dst_gdal):
    """crop the src_gdal file to dst_gdal by nodata value"""
    
    status = -1
    with gdal_datasource(src_gdal) as src_ds: 
        if src_ds is not None:
            ds_config = gdal_infos(src_ds)
            ds_arr = src_ds.GetRasterBand(1).ReadAsArray()
            ds_arr[ds_arr == ds_config['ndv']] = np.nan
            nans = np.isnan(ds_arr)

            nancols = np.all(nans, axis=0)
            nanrows = np.all(nans, axis=1)

            firstcol = nancols.argmin()
            firstrow = nanrows.argmin()        
            lastcol = len(nancols) - nancols[::-1].argmin()
            lastrow = len(nanrows) - nanrows[::-1].argmin()

            dst_arr = ds_arr[firstrow:lastrow,firstcol:lastcol]
            ds_arr = None

            dst_arr[np.isnan(dst_arr)] = ds_config['ndv']
            GeoT = ds_config['geoT']
            dst_x_origin = GeoT[0] + (GeoT[1] * firstcol)
            dst_y_origin = GeoT[3] + (GeoT[5] * firstrow)
            dst_geoT = [dst_x_origin, GeoT[1], 0.0, dst_y_origin, 0.0, GeoT[5]]
            ds_config['geoT'] = dst_geoT
            ds_config['nx'] = int(lastcol - firstcol)
            ds_config['ny'] = int(lastrow - firstrow)
            ds_config['nb'] = int(ds_config['nx'] * ds_config['ny'])
            status = gdal_write(dst_arr, dst_gdal, ds_config)
        else:
            return(None, status)

    return(status)
        
def gdal_split(src_gdal, split_value = 0, band = 1):
    """split raster file `src_dem`into two files based on z value, 
    or if split_value is a filename, split raster by overlay, where upper is outside and lower is inside.

    -----------
    Returns:
    list: [upper_grid-fn, lower_grid-fn]
    """

    def np_split(src_arr, sv = 0, nd = -9999):
        """split numpy `src_arr` by `sv` (turn u/l into `nd`)
        
        returns [upper_array, lower_array]
        """
        
        sv = utils.int_or(sv, 0)
        u_arr = np.array(src_arr)
        l_arr = np.array(src_arr)
        u_arr[u_arr <= sv] = nd
        l_arr[l_arr >= sv] = nd
        return(u_arr, l_arr)
    
    dst_upper = os.path.join(os.path.dirname(src_gdal), '{}_u.tif'.format(os.path.basename(src_gdal)[:-4]))
    dst_lower = os.path.join(os.path.dirname(src_gdal), '{}_l.tif'.format(os.path.basename(src_gdal)[:-4]))

    with gdal_datasource(src_gdal) as src_ds:        
        if src_ds is not None:
            src_config = gdal_infos(src_ds)
            dst_config = gdal_copy_infos(src_config)
            dst_config['fmt'] = 'GTiff'
            ds_arr = src_ds.GetRasterBand(band).ReadAsArray(0, 0, src_config['nx'], src_config['ny'])
            ua, la = np_split(ds_arr, split_value, src_config['ndv'])
            gdal_write(ua, dst_upper, dst_config)
            gdal_write(la, dst_lower, dst_config)
            ua = la = ds_arr = src_ds = None
            
    return([dst_upper, dst_lower])
        
def gdal_percentile(src_gdal, perc = 95, band = 1):
    """calculate the `perc` percentile of src_gdal file.

    -----------
    Returns:
    the calculated percentile
    """

    p = None
    with gdal_datasource(src_gdal) as src_ds:        
        if src_ds is not None:
            ds_array = src_ds.GetRasterBand(band).ReadAsArray().astype(float)
            ds_array[ds_array == src_ds.GetRasterBand(band).GetNoDataValue()] = np.nan
            x_dim = ds_array.shape[0]
            ds_array_flat = ds_array.flatten()
            ds_array = ds_array_flat[ds_array_flat != 0]
            if len(ds_array) > 0:
                p = np.nanpercentile(ds_array, perc)
                #percentile = 2 if p < 2 else p
            else: p = 2
            
    return(p)

def gdal_slope(src_gdal, dst_gdal, s = 111120):
    """generate a slope grid with GDAL

    -----------
    Returns:
    cmd output and status
    """
    
    gds_cmd = 'gdaldem slope {} {} {} -compute_edges'.format(src_gdal, dst_gdal, '' if s is None else '-s {}'.format(s))
    return(utils.run_cmd(gds_cmd))
    
def gdal_proximity(src_gdal, dst_gdal, band = 1, distunits='pixel'):
    """compute a proximity grid via GDAL

    -----------
    Returns:
    0 if success else None
    """

    distunits = utils.str_or(distunits, 'PIXEL')
    if distunits.upper() not in ['PIXEL', 'GEO']:
        distunits = 'PIXEL'
        
    prog_func = None
    with gdal_datasource(src_gdal) as src_ds:
        if src_ds is not None:
            ds_config = gdal_infos(src_ds)
            mem_ds = gdal_mem_ds(ds_config, name = 'MEM', bands = 1, src_srs = None)
            src_band = src_ds.GetRasterBand(band)
            src_arr = src_band.ReadAsArray()
            src_arr[src_arr != ds_config['ndv']] = 1
            src_arr[src_arr == ds_config['ndv']] = 0
            mem_band = mem_ds.GetRasterBand(1)
            mem_band.WriteArray(src_arr)

    drv = gdal.GetDriverByName('GTiff')
    dst_ds = drv.Create(dst_gdal, ds_config['nx'], ds_config['ny'], 1, gdal.GDT_Int32 if distunits == 'PIXEL' else gdal.GDT_Float32, [])
    dst_ds.SetGeoTransform(ds_config['geoT'])
    dst_ds.SetProjection(ds_config['proj'])
    dst_band = dst_ds.GetRasterBand(1)
    dst_band.SetNoDataValue(ds_config['ndv'])
    gdal.ComputeProximity(mem_band, dst_band, ['VALUES=1', 'DISTUNITS={}'.format(distunits)], callback = prog_func)
    mem_ds = dst_ds = None
    return(dst_gdal)
        
def gdal_polygonize(src_gdal, dst_layer, verbose = False):
    '''run gdal.Polygonize on src_ds and add polygon to dst_layer'''

    status = -1
    with gdal_datasource(src_gdal) as src_ds:
        if src_ds is not None:
            ds_arr = src_ds.GetRasterBand(1)
            if verbose:
                utils.echo_msg('polygonizing {}...'.format(src_gdal))
                
            status = gdal.Polygonize(
                ds_arr, None, dst_layer, 0,
                callback=gdal.TermProgress if verbose else None
            )
            if verbose:
                utils.echo_msg('polygonized {}'.format(src_gdal))

    return(status, status)

def gdal_mask(src_gdal, msk_gdal, out_gdal, msk_value = None, co=["COMPRESS=DEFLATE", "TILED=YES"], verbose = True):
    """mask the src_gdal file with the msk_gdal file to out_gdal"""
    
    with gdal_datasource(src_gdal) as src_ds:
        if src_ds is not None:
            src_config = gdal_infos(src_ds)
            src_band = src_ds.GetRasterBand(1)
            src_array = src_band.ReadAsArray()

            tmp_region = regions.Region().from_geo_transform(src_config['geoT'], src_config['nx'], src_config['ny'])
            #tmp_ds = gdal.Open(msk_dem)
            with gdal_datasource(msk_gdal) as tmp_ds:
                msk_ds = sample_warp(
                    tmp_ds, None, src_config['geoT'][1], src_config['geoT'][5],
                    src_region=tmp_region, sample_alg='bilinear', co=co,
                    verbose=verbose
                )[0] 
                if msk_ds is not None:
                    msk_band = msk_ds.GetRasterBand(1)
                    msk_array = msk_band.ReadAsArray()

                    if msk_value is None:
                        msk_value = msk_band.GetNoDataValue()

                    src_array[msk_array == msk_value] = src_band.GetNoDataValue()

                    gdal_write(src_array, out_gdal, src_config)
                    msk_ds = None
                    
def sample_warp(
        src_dem, dst_dem, x_sample_inc, y_sample_inc,
        src_srs=None, dst_srs=None, src_region=None, sample_alg='bilinear',
        ndv=-9999, tap=False, size=False, co=["COMPRESS=DEFLATE", "TILED=YES"],
        ot=gdal.GDT_Float32, coordinateOperation=None, verbose=False
):
    """sample and/or warp the src_dem"""

    xcount = ycount = 0
    out_region = None
    if src_region is not None:
        out_region = [src_region.xmin, src_region.ymin, src_region.xmax, src_region.ymax]        
        if src_srs is not None and dst_srs is not None:
            trans_region = src_region.copy()
            trans_region.src_srs = dst_srs
            trans_region.warp(src_srs)
        else:
            trans_region = src_region

        if x_sample_inc is None and y_sample_inc is None:
            src_infos = gdal_infos(src_dem)
            xcount, ycount, dst_gt = trans_region.geo_transform(
                x_inc=src_infos['geoT'][1], y_inc=src_infos['geoT'][5]*-1, node='grid'
            )
        
        if size and (xcount is None and ycount is None):
            xcount, ycount, dst_gt = src_region.geo_transform(
                x_inc=x_sample_inc, y_inc=y_sample_inc, node='grid'
            )
            x_sample_inc = y_sample_inc = None
        
    # if verbose:
    #     utils.echo_msg(
    #         'warping DEM: {} :: R:{} E:{}/{}:{}/{} S{} P{} -> T{}'.format(
    #             os.path.basename(str(src_dem)), out_region, x_sample_inc, y_sample_inc, xcount, ycount, sample_alg, src_srs, dst_srs
    #         )
    #     )        
    #utils.echo_msg('gdalwarp {} {} -s_srs {} -t_srs {} -tr {} {} -r {}'.format(src_dem, dst_dem, src_srs, dst_srs, x_sample_inc, y_sample_inc, sample_alg))
    #utils.echo_msg('dst_dem: > {} <'.format(dst_dem))
    
    if dst_dem is not None:
        if not os.path.exists(os.path.dirname(dst_dem)):
            os.makedirs(os.path.dirname(dst_dem))

    if verbose:
        desc = 'warping DEM: {} :: R:{} E:{}/{}:{}/{} S{} P{} -> T{}'.format(
            os.path.basename(str(src_dem)), out_region, x_sample_inc, y_sample_inc, xcount, ycount, sample_alg, src_srs, dst_srs
        )
        pbar = tqdm(desc=desc, total=100, leave=True)
        pbar_update = lambda a,b,c: pbar.update((a*100)-pbar.n)
    else:
        pbar_update = None

    #src_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    dst_ds = gdal.Warp('' if dst_dem is None else dst_dem, src_dem, format='MEM' if dst_dem is None else 'GTiff',
                       xRes=x_sample_inc, yRes=y_sample_inc, targetAlignedPixels=tap, width=xcount, height=ycount,
                       dstNodata=ndv, outputBounds=out_region, outputBoundsSRS=dst_srs if out_region is not None else None,
                       resampleAlg=sample_alg, errorThreshold=0, creationOptions=co, srcSRS=src_srs, dstSRS=dst_srs,
                       coordinateOperation=coordinateOperation, outputType=ot, callback=pbar_update)

    if verbose:
        pbar.close()

    # if dst_srs is not None:
    #     gdal_set_srs(dst_ds, dst_srs)
        
    if dst_dem is None:
        return(dst_ds, 0)
    else:
        dst_ds = None
        return(dst_dem, 0)    
    
def gdal_write(
        src_arr, dst_gdal, ds_config, dst_fmt='GTiff', co=["COMPRESS=DEFLATE", "TILED=YES"],
        max_cache=False, verbose=False
):
    """write src_arr to gdal file dst_gdal using src_config

    -----------
    Returns:
    list: [output-gdal, status-code]
    """
    
    driver = gdal.GetDriverByName(dst_fmt)
    if os.path.exists(dst_gdal):
        try:
            driver.Delete(dst_gdal)
        except Exception as e:
            echo_error_msg(e)
            remove_glob(dst_gdal)

    try:
        if not os.path.exists(os.path.dirname(dst_gdal)):
            os.makedirs(os.path.dirname(dst_gdal))
    except:
        pass
            
    if max_cache:
        gdal.SetCacheMax(2**30)

    #if ds_config['dt'] == 5:
    ds = driver.Create(dst_gdal, ds_config['nx'], ds_config['ny'], 1,
                       ds_config['dt'], options=co)
    #else:
    #ds = driver.Create(dst_gdal, ds_config['nx'], ds_config['ny'], 1,
    #ds_config['dt'], options=['COMPRESS=LZW', 'TILED=YES', 'PREDICTOR=3'])

    if ds is not None:
        ds.SetGeoTransform(ds_config['geoT'])
        try:
            ds.SetProjection(ds_config['proj'])
        except Exception as e:
            if verbose:
                echo_warning_msg('could not set projection {}'.format(ds_config['proj']))
            else: pass
        ds.GetRasterBand(1).SetNoDataValue(ds_config['ndv'])
        if src_arr is not None:
            ds.GetRasterBand(1).WriteArray(src_arr)
            
        ds = src_arr = None        
        return(dst_gdal, 0)
    else:
        return(None, -1)

def gdal2gdal(src_dem, dst_fmt='GTiff', src_srs='epsg:4326', dst_dem=None, co=True):
    """convert the gdal file to gdal using gdal

    -----------
    Returns:
    output-gdal-fn
    """
    
    if os.path.exists(src_dem):
        if dst_dem is None:
            #dst_dem = '{}.{}'.format(os.path.basename(src_dem).split('.')[0], gdal_fext(dst_fmt))
            dst_dem = '{}.{}'.format(utils.fn_basename2(src_dem), gdal_fext(dst_fmt))
            
        if dst_fmt != 'GTiff':
            co = False
            
        if not co: # update co
            gdal2gdal_cmd = ('gdal_translate {} {} -f {}'.format(src_dem, dst_dem, dst_fmt))
        else:
            gdal2gdal_cmd = ('gdal_translate {} {} -f {} -co TILED=YES -co COMPRESS=DEFLATE\
            '.format(src_dem, dst_dem, dst_fmt))
            
        out, status = utils.run_cmd(gdal2gdal_cmd, verbose=False)
        if status == 0:
            return(dst_dem)
        else:
            return(None)
    else:
        return(None)

def gdal_yield_srcwin(src_gdal, n_chunk = 10, step = 5, verbose = False):
    """yield source windows in n_chunks at step"""
    
    ds_config = gdal_infos(src_gdal)
    gt = ds_config['geoT']
    x_chunk = n_chunk
    y_chunk = 0
    i_chunk = 0
    x_i_chunk = 0
    
    with tqdm(total=math.ceil(ds_config['ny']/step) * math.ceil(ds_config['nx']/step), desc='chunking {}'.format(src_gdal)) as pbar:
        while True:
            #pbar.update()                    
            y_chunk = n_chunk
            while True:
                pbar.update()
                this_x_chunk = ds_config['nx'] if x_chunk > ds_config['nx'] else x_chunk
                this_y_chunk = ds_config['ny'] if y_chunk > ds_config['ny'] else y_chunk
                this_x_origin = x_chunk - n_chunk
                this_y_origin = y_chunk - n_chunk
                this_x_size = int(this_x_chunk - this_x_origin)
                this_y_size = int(this_y_chunk - this_y_origin)
                if this_x_size == 0 or this_y_size == 0: break
                srcwin = (this_x_origin, this_y_origin, this_x_size, this_y_size)
                yield(srcwin)

                if y_chunk > ds_config['ny']:
                    break
                else:
                    y_chunk += step
                    i_chunk += 1
                    
                #pbar.update(1)
                #pbar.update()
            if x_chunk > ds_config['nx']:
                break
            else:
                x_chunk += step
                x_i_chunk += 1
    
def gdal_chunks(src_gdal, n_chunk, band = 1):
    """split `src_gdal` GDAL file into chunks with `n_chunk` cells squared.

    -----------
    Returns:
    list of chunked filenames.
    """
    
    o_chunks = []
    with gdal_datasource(src_gdal) as src_ds:
        if src_ds is not None:
            ds_config = gdal_infos(src_ds)
            src_band = src_ds.GetRasterBand(band)
            gt = ds_config['geoT']
            gt = list(gt)
            gt[0] = gt[0] - (gt[1]/2)
            gt[3] = gt[3] - (gt[5]/2)
            gt = tuple(gt)

            c_n = 0
            for srcwin in gdal_yield_srcwin(src_gdal, n_chunk = n_chunk, step = n_chunk):
                this_geo_x_origin, this_geo_y_origin = utils._pixel2geo(srcwin[0], srcwin[1], gt)
                dst_gt = [this_geo_x_origin, float(gt[1]), 0.0, this_geo_y_origin, 0.0, float(gt[5])]
                band_data = src_band.ReadAsArray(srcwin[0], srcwin[1], srcwin[2], srcwin[3])
                if not np.all(band_data == band_data[0,:]):
                    dst_config = gdal_copy_infos(ds_config)
                    dst_config['nx'] = srcwin[2]
                    dst_config['ny'] = srcwin[3]
                    dst_config['geoT'] = dst_gt
                    this_region = regions.Region().from_geo_transform(dst_gt, dst_config['nx'], dst_config['ny'])
                    o_chunk = '{}_chnk{}.tif'.format(os.path.basename(src_gdal).split('.')[0], c_n)
                    dst_fn = os.path.join(os.path.dirname(src_gdal), o_chunk)
                    o_chunks.append(dst_fn)
                    gdal_write(band_data, dst_fn, dst_config)
                    c_n += 1                
        return(o_chunks)

def gdal_parse(src_ds, dump_nodata = False, srcwin = None, mask = None, dst_srs = None,
          verbose = False, z_region = None, step = 1):
    """parse the data from gdal dataset src_ds (first band only)
    optionally mask the output with `mask` or transform the coordinates to `dst_srs`

    yields the parsed xyz data
    """

    #if verbose: sys.stderr.write('waffles: parsing gdal file {}...'.format(src_ds.GetDescription()))
    ln = 0
    band = src_ds.GetRasterBand(1)
    ds_config = gdal_infos(src_ds)
    src_srs = osr.SpatialReference()
    try:
        src_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    except:
        pass
    
    src_srs.ImportFromWkt(ds_config['proj'])
    src_srs.AutoIdentifyEPSG()
    srs_auth = src_srs.GetAuthorityCode(None)

    if srs_auth is None:
        src_srs.ImportFromEPSG(4326)
        src_srs.AutoIdentifyEPSG()
        srs_auth = src_srs.GetAuthorityCode(None)

    #if srs_auth == warp: dst_srs = None

    if dst_srs is not None:
        dst_srs_ = osr.SpatialReference()
        dst_srs_.SetFromUserInput(dst_srs)
        ## GDAL 3+
        try:
            dst_srs_.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        except: pass
        dst_trans = osr.CoordinateTransformation(src_srs, dst_srs_)

    gt = ds_config['geoT']
    msk_band = None
    if mask is not None:
        src_mask = gdal.Open(mask)
        msk_band = src_mask.GetRasterBand(1)
    if srcwin is None:
        srcwin = (0, 0, ds_config['nx'], ds_config['ny'])
        
    nodata = ['{:g}'.format(-9999), 'nan', float('nan')]
    if band.GetNoDataValue() is not None:
        nodata.append('{:g}'.format(band.GetNoDataValue()))
        
    if dump_nodata:
        nodata = []
        
    for y in range(srcwin[1], srcwin[1] + srcwin[3], 1):
        band_data = band.ReadAsArray(srcwin[0], y, srcwin[2], 1)
        if z_region is not None:
            if z_region[0] is not None:
                band_data[band_data < z_region[0]] = -9999
                
            if z_region[1] is not None:
                band_data[band_data > z_region[1]] = -9999
                
        if msk_band is not None:
            msk_data = msk_band.ReadAsArray(srcwin[0], y, srcwin[2], 1)
            band_data[msk_data==0]=-9999
            
        band_data = np.reshape(band_data, (srcwin[2], ))
        for x_i in range(0, srcwin[2], 1):
            x = x_i + srcwin[0]
            z = band_data[x_i]
            if '{:g}'.format(z) not in nodata:
                ln += 1
                geo_x,geo_y = utils._pixel2geo(x, y, gt)
                # if warp is not None:
                #     #point = xyzfun.XYZPoint(x=geo_x, y=geo_y, z=z, epsg=4326)
                #     point = ogr.CreateGeometryFromWkt('POINT ({} {})'.format(geo_x, geo_y))
                #     point.Transform(dst_trans)
                #     pnt = point.GetPoint()
                #     xyz = xyzfun.XYZPoint(x=point.GetX(), y=point.GetY(), z=z)#from_list([pnt[0], pnt[1], z])
                # else:
                xyz = xyzfun.XYZPoint(x=geo_x, y=geo_y, z=z)#line = [geo_x, geo_y, z]
                yield(xyz)
    band = None
    src_mask = None
    msk_band = None
    if verbose: utils.echo_msg('parsed {} data records from {}'.format(ln, src_ds.GetDescription()))

def gdal_point_query(points, src_gdal, x = 'x', y = 'y', z = 'z', band = 1):
    with gdal_datasource(src_gdal) as src_ds:
        if src_ds is not None:
            ds_config = gdal_infos(src_ds)
            ds_band = src_ds.GetRasterBand(band)
            tgrid = ds_band.ReadAsArray()

    t_region = regions.Region().from_geo_transform(ds_config['geoT'], tgrid.shape[0], tgrid.shape[1])
    xcount, ycount, dst_gt = t_region.geo_transform(
        x_inc=ds_config['geoT'][1], y_inc=ds_config['geoT'][5]*-1, node='grid'
    )
    pixel_x = np.floor((points['x'] - dst_gt[0]) / dst_gt[1]).astype(int)
    pixel_y = np.floor((points['y'] - dst_gt[3]) / dst_gt[5]).astype(int)

    out_idx = np.nonzero((pixel_x >= xcount) | (pixel_x < 0) | (pixel_y >= ycount) | (pixel_y < 0))
    out_idx2 = np.nonzero((pixel_x < xcount) | (pixel_x > 0) | (pixel_y < ycount) | (pixel_y > 0))
    #pixel_x = np.delete(pixel_x, out_idx)
    #pixel_y = np.delete(pixel_y, out_idx)
    
    pixel_zz = tgrid[pixel_x[out_idx2], pixel_y[out_idx2]]
    points['z'][out_idx2] = pixel_zz + points['z'][out_idx2]
    return(points)
    
def gdal_yield_query(src_xyz, src_gdal, out_form, band = 1):
    """query a gdal-compatible grid file with xyz data.
    out_form dictates return values

    yields out_form results
    """

    tgrid = None
    with gdal_datasource(src_gdal) as src_ds:
        if src_ds is not None:
            ds_config = gdal_infos(src_ds)
            ds_band = src_ds.GetRasterBand(band)
            ds_gt = ds_config['geoT']
            ds_nd = ds_config['ndv']
            tgrid = ds_band.ReadAsArray()

    if tgrid is not None:
        for xyz in src_xyz:
            x = xyz[0]
            y = xyz[1]
            try: 
                z = xyz[2]
            except:
                z = ds_nd

            if x > ds_gt[0] and y < float(ds_gt[3]):
                xpos, ypos = utils._geo2pixel(x, y, ds_gt, node='pixel')
                try: 
                    g = tgrid[ypos, xpos]
                except: g = ds_nd
                d = c = m = s = ds_nd
                if g != ds_nd:
                    d = z - g
                    m = z + g
                    outs = []
                    for i in out_form:
                        outs.append(vars()[i])

                    #print(g)
                    #print(outs)
                    yield(outs)

def gdal_query(src_xyz, src_gdal, out_form, band = 1):
    """query a gdal-compatible grid file with xyz data.
    out_form dictates return values

    -----------
    Returns:
    array of values
    """
    
    xyzl = []
    for out_q in gdal_yield_query(src_xyz, src_gdal, out_form, band=band):
        xyzl.append(np.array(out_q))
        
    return(np.array(xyzl))
                
### End
