### gdalfun.py - GDAL functions
##
## Copyright (c) 2010 - 2023 Regents of the University of Colorado
##
## demfun.py is part of CUDEM
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
import shutil

from osgeo import gdal
from osgeo import osr
from osgeo import ogr
import numpy as np

from cudem import utils
from cudem import regions
from cudem import xyzfun

## ==============================================
## OSR/WKT
## ==============================================
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
    """generate a .prj file given a src_srs"""
    
    with open(dst_fn, 'w') as out:
        out.write(osr_wkt(src_srs, True))
        
    return(0)

def epsg_from_input(in_srs):
    src_srs = osr.SpatialReference()
    src_srs.SetFromUserInput(in_srs)

    ## HORZ
    if src_srs.IsGeographic() == 1:
        cstype = 'GEOGCS'
    else:
        cstype = 'PROJCS'

    src_srs.AutoIdentifyEPSG()
    an = src_srs.GetAuthorityName(cstype)
    src_horz_epsg = src_srs.GetAuthorityCode(cstype)

    ## VERT
    if src_srs.IsVertical() == 1:
        csvtype = 'VERT_CS'
        src_vert_epsg = src_srs.GetAuthorityCode(csvtype)
    else:
        src_vert_epsg = None

    return(src_horz_epsg, src_vert_epsg)
    
def osr_parse_srs(src_srs):
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

## ==============================================
## OGR
## ==============================================
def ogr_fext(src_drv_name):
    """find the common file extention given a OGR driver name
    older versions of gdal can't do this, so fallback to known standards.

    Args:
      src_drv_name (str): a source OGR driver name

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


def gdal_get_srs(src_gdal):
    src_ds = gdal.Open(src_gdal)
    src_srs = src_ds.GetSpatialRef()
    src_ds = None
    return(osr_parse_srs(src_srs))
    
def ogr_get_srs(src_ogr):
    src_ds = ogr.Open(src_ogr, 0)
    src_srs = src_ds.GetLayer().GetSpatialRef()
    src_ds = None
    return(osr_parse_srs(src_srs))
    
def ogr_clip(src_ogr_fn, dst_region=None, layer=None, overwrite=False):

    dst_ogr_bn = '.'.join(src_ogr_fn.split('.')[:-1])
    dst_ogr_fn = '{}_{}.gpkg'.format(dst_ogr_bn, dst_region.format('fn'))
    
    if not os.path.exists(dst_ogr_fn) or overwrite:
        run_cmd('ogr2ogr -nlt PROMOTE_TO_MULTI {} {} -clipsrc {} {} '.format(dst_ogr_fn, src_ogr_fn, dst_region.format('te'), layer if layer is not None else ''), verbose=True)
                
    return(dst_ogr_fn)

def ogr_clip2(src_ogr_fn, dst_region=None, layer=None, overwrite=False):

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
        remove_glob('{}.*'.format('.'.join(region_ogr.split('.')[:-1])))
        
    return(dst_ogr_fn)

def ogr_clip3(src_ogr, dst_ogr, clip_region=None, dn="ESRI Shapefile"):
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

def ogr_mask_union(src_layer, src_field, dst_defn=None):
    """`union` a `src_layer`'s features based on `src_field` where
    `src_field` holds a value of 0 or 1. optionally, specify
    an output layer defn for the unioned feature.

    returns the output feature class
    """
    
    if dst_defn is None:
        dst_defn = src_layer.GetLayerDefn()
        
    multi = ogr.Geometry(ogr.wkbMultiPolygon)
    src_layer.SetAttributeFilter("{} = 1".format(src_field))
    feats = len(src_layer)
    
    if feats > 0:
        with utils.CliProgress(total=len(src_layer), message='unioning {} features...'.format(feats)) as pbar:
            for n, f in enumerate(src_layer):
                pbar.update()
                f_geom = f.geometry()
                f_geom_valid = f_geom
                multi.AddGeometry(f_geom_valid)
            
    utils.echo_msg('setting geometry to unioned feature...')
    out_feat = ogr.Feature(dst_defn)
    out_feat.SetGeometry(multi)
    union = multi = None
    
    return(out_feat)

def ogr_empty_p(src_ogr, dn='ESRI Shapefile'):
    driver = ogr.GetDriverByName(dn)
    ds = driver.Open(src_ogr, 0)

    if ds is not None:
        layer = ds.GetLayer()
        fc = layer.GetFeatureCount()
        if fc == 0:
            return(True)
        else: return(False)
        
    else: return(True)

def ogr_polygonize_multibands(src_ds, dst_srs = 'epsg:4326', ogr_format='ESRI Shapefile', verbose=True):
    dst_layer = '{}_sm'.format(fn_basename2(src_ds.GetDescription()))
    dst_vector = dst_layer + '.{}'.format(ogr_fext(ogr_format))
    remove_glob('{}.*'.format(dst_layer))
    gdal_prj_file('{}.prj'.format(dst_layer), dst_srs)
    ds = ogr.GetDriverByName(ogr_format).CreateDataSource(dst_vector)
    if ds is not None: 
        layer = ds.CreateLayer(
            '{}'.format(dst_layer), None, ogr.wkbMultiPolygon
        )
        [layer.SetFeature(feature) for feature in layer]
    else:
        layer = None

    layer.CreateField(ogr.FieldDefn('DN', ogr.OFTInteger))
    defn = None

    for b in range(1, src_ds.RasterCount+1):
        this_band = src_ds.GetRasterBand(b)
        this_band_md = this_band.GetMetadata()
        b_infos = gdal_gather_infos(src_ds, scan=True, band=b)
        field_names = [field.name for field in layer.schema]
        for k in this_band_md.keys():
            if k[:9] not in field_names:
                layer.CreateField(ogr.FieldDefn(k[:9], ogr.OFTString))
        
        if b_infos['zr'][1] == 1:
            tmp_ds = ogr.GetDriverByName('Memory').CreateDataSource(
                '{}_poly'.format(this_band.GetDescription())
            )
                    
            if tmp_ds is not None:
                tmp_layer = tmp_ds.CreateLayer(
                    '{}_poly'.format(this_band.GetDescription()), None, ogr.wkbMultiPolygon
                )
                tmp_layer.CreateField(ogr.FieldDefn('DN', ogr.OFTInteger))
                for k in this_band_md.keys():
                    tmp_layer.CreateField(ogr.FieldDefn(k[:9], ogr.OFTString))
                
                if verbose:
                    echo_msg('polygonizing {} mask...'.format(this_band.GetDescription()))
                            
                status = gdal.Polygonize(
                    this_band,
                    None,
                    tmp_layer,
                    tmp_layer.GetLayerDefn().GetFieldIndex('DN'),
                    [],
                    callback = gdal.TermProgress if verbose else None
                )

                if len(tmp_layer) > 0:
                    if defn is None:
                        defn = tmp_layer.GetLayerDefn()

                    out_feat = gdal_ogr_mask_union(tmp_layer, 'DN', defn)
                    echo_msg('creating feature {}...'.format(this_band.GetDescription()))
                    for k in this_band_md.keys():
                        out_feat.SetField(k[:9], this_band_md[k])
                    
                    layer.CreateFeature(out_feat)

            if verbose:
                echo_msg('polygonized {}'.format(this_band.GetDescription()))
            tmp_ds = tmp_layer = None
    ds = None

## ==============================================
## GDAL
## ==============================================
class gdal_datasource:
    """

    use as:
    with gdal_datasource('input.tif') as src_ds:
        src_ds.do_things()

    where the input can be either a path-name to a gdal supported
    file or a gdal.Dataset object.
    """
    
    def __init__(self, src_gdal = None, update = False):
        self.src_gdal = src_gdal
        self.update = update
        self.src_ds = None

    ## todo: create new datasource if src_gdal is path, but doesn't exists...
    def __enter__(self):
        if isinstance(self.src_gdal, gdal.Dataset):
            self.src_ds = self.src_gdal

        elif utils.str_or(self.src_gdal) is not None and os.path.exists(self.src_gdal):
            self.src_ds = gdal.Open(self.src_gdal, 0 if not self.update else 1)

        return(self.src_ds)

    def __exit__(self, exc_type, exc_value, exc_traceback):
        ## don't close src_ds if incoming was already open
        print(exc_type, exc_value, exc_traceback)
        if not isinstance(self.src_gdal, gdal.Dataset):
            self.src_ds = None
            
def gdal_fext(src_drv_name):
    """find the common file extention given a GDAL driver name
    older versions of gdal can't do this, so fallback to known standards.

    Args:
      src_drv_name (str): a source GDAL driver name

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

def gdal_set_infos(nx, ny, nb, geoT, proj, dt, ndv, fmt):
    """set a datasource config dictionary

    returns gdal_config dict.
    """
    
    return({'nx': nx, 'ny': ny, 'nb': nb, 'geoT': geoT, 'proj': proj, 'dt': dt,
            'ndv': ndv, 'fmt': fmt, 'metadata': md, 'raster_count': rc})
    
def gdal_infos(src_gdal, region = None, scan = False, band = 1):
    with gdal_datasource(src_gdal) as src_ds:
        if src_ds is not None:
            if region is not None:
                srcwin = region.srcwin(gt, src_ds.RasterXSize, src_ds.RasterYSize)
            else:
                srcwin = (0, 0, src_ds.RasterXSize, src_ds.RasterYSize)

            gt = src_ds.GetGeoTransform()        
            dst_gt = (gt[0] + (srcwin[0] * gt[1]), gt[1], 0., gt[3] + (srcwin[1] * gt[5]), 0., gt[5])
            src_band = src_ds.GetRasterBand(band)
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
                'zr': None,
            }
            if ds_config['ndv'] is None:
                ds_config['ndv'] = -9999

            if scan:
                src_arr = src_band.ReadAsArray(srcwin[0], srcwin[1], srcwin[2], srcwin[3])
                ds_config['zr'] = src_band.ComputeRasterMinMax()
                src_arr = src_band = None

            return(ds_config)

def gdal_set_srs(src_gdal, src_srs = 'epsg:4326', verbose = True):
    with gdal_datasource(src_gdal) as src_ds:    
        if src_ds is not None and src_srs is not None:
            try:
                src_ds.SetProjection(osr_wkt(src_srs))
                return(0)
            except Exception as e:
                if verbose:
                    utils.echo_warning_msg('could not set projection {}'.format(src_srs))
                return(None)
        else:
            return(None)    

def gdal_get_ndv(src_gdal, band = 1):
    with gdal_datasource(src_gdal) as src_ds:    
        if src_ds is not None:
            src_band = src_ds.GetRasterBand(band)
            return(src_band.GetNoDataValue())

def gdal_set_ndv(src_gdal, ndv = -9999, convert_array = False, verbose  =True):
    """set the nodata value of gdal datasource

    returns 0
    """

    with gdal_datasource(src_gdal) as src_ds:        
        if src_ds is not None:
            ds_config = gdal_gather_infos(src_ds)
            curr_nodata = ds_config['ndv']

            for band in range(1, src_ds.RasterCount+1):
                this_band = src_ds.GetRasterBand(band)
                this_band.DeleteNoDataValue()

            for band in range(1, src_ds.RasterCount+1):
                this_band = src_ds.GetRasterBand(band)
                this_band.SetNoDataValue(nodata)

            if convert_array:
                for band in range(1, src_ds.RasterCount+1):
                    this_band = src_ds.GetRasterBand(band)
                    arr = this_band.ReadAsArray()
                    if np.isnan(curr_nodata):
                        arr[np.isnan(arr)]=nodata
                    else:
                        arr[arr == curr_nodata] = nodata

                    this_band.WriteArray(arr)            
            return(0)
        else:
            return(None)
    
def gdal_generate_mem_ds(ds_config, name = 'MEM', bands = 1):
    """Create temporary gdal mem dataset"""
        
    mem_driver = gdal.GetDriverByName('MEM')
    mem_ds = mem_driver.Create(name, ds_config['nx'], ds_config['ny'], bands, ds_config['dt'])
    if mem_ds is not None:
        mem_ds.SetGeoTransform(ds_config['geoT'])
        if ds_config['proj'] is not None:
            mem_ds.SetProjection(ds_config['proj'])

        for band in range(1, bands+1):
            mem_band = mem_ds.GetRasterBand(band)
            mem_band.SetNoDataValue(ds_config['ndv'])
        
    return(mem_ds)

def gdal_cut(src_gdal, src_region, dst_gdal, node='pixel', verbose=True):
    """cut src_ds datasource to src_region and output dst_gdal file

    returns [output-dem, status-code]
    """

    with gdal_datasource(src_gdal) as src_ds:        
        if src_ds is not None:
            ds_config = gdal_gather_infos(src_ds)
            gt = ds_config['geoT']
            srcwin = src_region.srcwin(gt, src_ds.RasterXSize, src_ds.RasterYSize, node=node)
            dst_gt = (gt[0] + (srcwin[0] * gt[1]), gt[1], 0., gt[3] + (srcwin[1] * gt[5]), 0., gt[5])
            out_ds_config = gdal_set_infos(srcwin[2], srcwin[3], srcwin[2] * srcwin[3], dst_gt,
                                           ds_config['proj'], ds_config['dt'], ds_config['ndv'],
                                           ds_config['fmt'])

            in_bands = src_ds.RasterCount
            mem_ds = gdal_generate_mem_ds(out_ds_config, bands=in_bands)
            for band in range(1, in_bands+1):
                this_band = mem_ds.GetRasterBand(band)
                this_band.WriteArray(src_ds.GetRasterBand(band).ReadAsArray(*srcwin))
                mem_ds.FlushCache()

            dst_ds = gdal.GetDriverByName(ds_config['fmt']).CreateCopy(dst_gdal, mem_ds, 0)                
            dst_ds = None
            return(dst_gdal, 0)
        else:
            return(None, -1)

def gdal_percentile(src_gdal, perc = 95, band = 1):
    """calculate the `perc` percentile of src_fn gdal file.

    return the calculated percentile
    """

    with gdal_datasource(src_gdal) as src_ds:        
        if src_ds is not None:
            ds_array = np.array(src_ds.GetRasterBand(band).ReadAsArray())
            ds_array[ds_array == src_ds.GetRasterBand(band).GetNoDataValue()] = np.nan
            x_dim = ds_array.shape[0]
            ds_array_flat = ds_array.flatten()
            ds_array = ds_array_flat[ds_array_flat != 0]
            if len(ds_array) > 0:
                p = np.nanpercentile(ds_array, perc)
                #percentile = 2 if p < 2 else p
            else: p = 2
            return(p)
        else:
            return(None)
    
def gdal_write(src_arr, dst_gdal, ds_config, dst_fmt='GTiff', max_cache=False, verbose=False):
    """write src_arr to gdal file dst_gdal using src_config

    returns [output-gdal, status-code]
    """
    
    driver = gdal.GetDriverByName(dst_fmt)
    if os.path.exists(dst_gdal):
        try:
            driver.Delete(dst_gdal)
        except Exception as e:
            echo_error_msg(e)
            remove_glob(dst_gdal)

    if max_cache:
        gdal.SetCacheMax(2**30)

    if ds_config['dt'] == 5:
        ds = driver.Create(dst_gdal, ds_config['nx'], ds_config['ny'], 1,
                           ds_config['dt'], options=['COMPRESS=LZW', 'PREDICTOR=2', 'TILED=YES'])
    else:
        ds = driver.Create(dst_gdal, ds_config['nx'], ds_config['ny'], 1,
                           ds_config['dt'], options=['COMPRESS=LZW', 'TILED=YES', 'PREDICTOR=3'])

    if ds is not None:
        ds.SetGeoTransform(ds_config['geoT'])
        try:
            ds.SetProjection(ds_config['proj'])
        except Exception as e:
            if verbose:
                echo_warning_msg('could not set projection {}'.format(ds_config['proj']))
            else: pass
        ds.GetRasterBand(1).SetNoDataValue(ds_config['ndv'])
        ds.GetRasterBand(1).WriteArray(src_arr)
        ds = src_arr = None        
        return(dst_gdal, 0)
    else:
        return(None, -1)

def gdal2gdal(src_dem, dst_fmt='GTiff', src_srs='epsg:4326', dst_dem=None, co=True):
    """convert the gdal file to gdal using gdal

    return output-gdal-fn"""
    
    if os.path.exists(src_dem):
        if dst_dem is None:
            #dst_dem = '{}.{}'.format(os.path.basename(src_dem).split('.')[0], gdal_fext(dst_fmt))
            dst_dem = '{}.{}'.format(fn_basename2(src_dem), gdal_fext(dst_fmt))
            
        if dst_fmt != 'GTiff':
            co = False
            
        if not co:
            gdal2gdal_cmd = ('gdal_translate {} {} -f {}'.format(src_dem, dst_dem, dst_fmt))
        else:
            gdal2gdal_cmd = ('gdal_translate {} {} -f {} -co TILED=YES -co COMPRESS=DEFLATE\
            '.format(src_dem, dst_dem, dst_fmt))
            
        out, status = run_cmd(gdal2gdal_cmd, verbose=False)
        if status == 0:
            return(dst_dem)
        else:
            return(None)
    else:
        return(None)

### End
