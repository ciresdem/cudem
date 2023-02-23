### demfun.py - CUDEM DEM utilities and functions
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
## DEM/GDAL/GMT functions, etc
## ==============================================
def infos(src_dem, region=None, scan=False):
    """scan gdal file src_fn and gather region info.

    returns region dict.
    """
    
    #if os.path.exists(src_dem):
    try:
        ds = gdal.Open(src_dem)
    except: ds = None
    if ds is not None:
        dsc = gather_infos(ds, region = region, scan = scan)#, node='grid' if src_dem.split('.')[-1] == 'nc' else 'pixel')
        ds = None
        return(dsc)
    else: return(None)
    #else: return(None)
    
def gather_infos(src_ds, region=None, scan=False, node='pixel'):
    """gather information from `src_ds` GDAL dataset

    returns gdal_config dict.
    """

    gt = src_ds.GetGeoTransform()
    if region is not None:
        srcwin = region.srcwin(gt, src_ds.RasterXSize, src_ds.RasterYSize)
    else: srcwin = (0, 0, src_ds.RasterXSize, src_ds.RasterYSize)
    src_band = src_ds.GetRasterBand(1)
    dst_gt = (gt[0] + (srcwin[0] * gt[1]), gt[1], 0., gt[3] + (srcwin[1] * gt[5]), 0., gt[5])
    mt = src_ds.GetMetadata()
    
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
        'zr': None,
    }
    if ds_config['ndv'] is None: ds_config['ndv'] = -9999
    if scan:
        src_arr = src_band.ReadAsArray(srcwin[0], srcwin[1], srcwin[2], srcwin[3])
        #ds_config['zr'] = (np.amin(src_arr), np.amax(src_arr))
        ds_config['zr'] = src_band.ComputeRasterMinMax()
        src_arr = None
    return(ds_config)

def set_infos(nx, ny, nb, geoT, proj, dt, ndv, fmt):
    '''set a datasource config dictionary

    returns gdal_config dict.'''
    return({'nx': nx, 'ny': ny, 'nb': nb, 'geoT': geoT, 'proj': proj, 'dt': dt, 'ndv': ndv, 'fmt': fmt})

def copy_infos(src_config):
    """copy src_config

    returns copied src_config dict.
    """
    
    dst_config = {}
    for dsc in src_config.keys():
        dst_config[dsc] = src_config[dsc]
    return(dst_config)

def set_srs(src_dem, src_srs='epsg:4326', verbose=True):
    """set the projection of gdal file src_fn to src_srs"""

    try:
        ds = gdal.Open(src_dem, gdal.GA_Update)
    except: ds = None
    if ds is not None:
        try:
            ds.SetProjection(utils.sr_wkt(src_srs))
        except Exception as e:
            #ds.SetProjection(utils.sr_wkt('epsg:4326')) ## set to default if user input is no good
            if verbose:
                utils.echo_warning_msg('could not set projection {}'.format(src_srs))
        ds = None
        return(0)
    else: return(None)

def get_srs(src_dem):
    src_ds = gdal.Open(src_dem)
    src_srs = src_ds.GetSpatialRef()
    src_ds = None

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

def get_nodata(src_dem):
    """get the nodata valiue of src_dem

    return nodata value
    """
    
    src_ds = gdal.Open(src_dem)
    if src_ds is not None:
        band = src_ds.GetRasterBand(1)
        ndv = band.GetNoDataValue()
        src_ds = None
        return(ndv)
    else:
        return(None)
    
def set_nodata(src_dem, nodata=-9999, convert_array=False, verbose=True):
    """set the nodata value of gdal file src_dem

    returns 0
    """

    if verbose:
        utils.echo_msg('setting nodata value of {} to {}'.format(src_dem, nodata))
        
    try:
        ds = gdal.Open(src_dem, gdal.GA_Update)
    except: ds = None
    if ds is not None:
        ds_config = gather_infos(ds)
        curr_nodata = ds_config['ndv']
        band = ds.GetRasterBand(1)
        if convert_array:
            arr = band.ReadAsArray()
            if np.isnan(curr_nodata):
                arr[np.isnan(arr)]=nodata
            else:
                arr[arr == curr_nodata] = nodata
            band.SetNoDataValue(nodata)
            ds_config['ndv'] = nodata
            band.WriteArray(arr)
        else:
            band.SetNoDataValue(nodata)
        ds = None
        return(0)
    else: return(None)

def set_metadata(src_dem, node='pixel', cudem=False, verbose=True): #, vdatum='NAVD88'):
    """add metadata to the waffled raster

    Args: 
      cudem (bool): add CUDEM metadata
    """

    try:
        ds = gdal.Open(src_dem, gdal.GA_Update)
    except: ds = None
    if ds is not None:
        md = ds.GetMetadata()
        if node == 'pixel':
            md['AREA_OR_POINT'] = 'Area'
            md['NC_GLOBAL#node_offset'] = '1'
            md['tos#node_offset'] = '1'
        else:
            md['AREA_OR_POINT'] = 'Point'
            md['NC_GLOBAL#node_offset'] = '0'
            md['tos#node_offset'] = '0'
        md['TIFFTAG_DATETIME'] = '{}'.format(utils.this_date())
        if cudem:
            md['TIFFTAG_COPYRIGHT'] = 'DOC/NOAA/NESDIS/NCEI > National Centers for Environmental Information, NESDIS, NOAA, U.S. Department of Commerce'
            infos = gather_infos(ds, scan=True)
            if infos['zr'][1] < 0:
                tb = 'Bathymetry'
            elif infos['zr'][0] > 0:
                tb = 'Topography'
            else:
                tb = 'Topography-Bathymetry'

            prj=ds.GetProjection()
            srs=osr.SpatialReference(wkt=prj)
            vdatum=srs.GetAttrValue('vert_cs')
            md['TIFFTAG_IMAGEDESCRIPTION'] = '{}; {}'.format(tb, '' if vdatum is None else vdatum)

        ds.SetMetadata(md)
        ds = None
        return(0)
    else:
        if verbose:
            utils.echo_error_msg('failed to set metadata')
            
        return(None)

def get_array(src_dem):
    ds = gdal.Open(src_dem)
    if ds is not None:
        band = ds.GetRasterBand(1)
        _array = band.ReadAsArray()
        infos = gather_infos(ds)
        ds = None
        return(_array, infos)
    else:
        utils.echo_error_msg('could not open {}'.format(src_dem))
        return(None, None)

def extract_band(data_set, output_data_set, band=1, exclude=[], inverse=False):
    _ds = gdal.Open(data_set)
    _ds_config = gather_infos(_ds)
    _ds_band = _ds.GetRasterBand(band)
    _ds_array = _ds_band.ReadAsArray()
    _ds = None

    if _ds_config['ndv'] is None:
        _ds_config['ndv'] = -9999
    
    if exclude:
        for key in exclude:
            _ds_array[_ds_array == key] = _ds_config['ndv']

    if inverse:
        _ds_array = 1/_ds_array
            
    utils.gdal_write(_ds_array, output_data_set, _ds_config)
    return(output_data_set)
    
def mask_(src_dem, msk_dem, out_dem, msk_value = None):
    src_ds = gdal.Open(src_dem)
    if src_ds is not None:
        src_config = gather_infos(src_ds)
        src_band = src_ds.GetRasterBand(1)
        src_array = src_band.ReadAsArray()

        tmp_region = regions.Region().from_geo_transform(src_config['geoT'], src_config['nx'], src_config['ny'])
        tmp_ds = gdal.Open(msk_dem)
        msk_ds = sample_warp(
            tmp_ds, None, src_config['geoT'][1], src_config['geoT'][5],
            src_region=tmp_region, sample_alg='bilinear',
            verbose=True
        )[0] 
        tmp_ds = None
        
        if msk_ds is not None:
            msk_band = msk_ds.GetRasterBand(1)
            msk_array = msk_band.ReadAsArray()

            if msk_value is None:
                msk_value = msk_band.GetNoDataValue()
            
            src_array[msk_array == msk_value] = src_band.GetNoDataValue()

            utils.gdal_write(src_array, out_dem, src_config)
            msk_ds = None
        src_ds = None
        
def split(src_dem, split_value = 0):
    """split raster file `src_dem`into two files based on z value, 
    or if split_value is a filename, split raster by overlay, where upper is outside and lower is inside.

    returns [upper_gr
    id-fn, lower_grid-fn]
    """

    def np_split(src_arr, sv = 0, nd = -9999):
        """split numpy `src_arr` by `sv` (turn u/l into `nd`)
        
        returns [upper_array, lower_array]
        """
        
        try:
            sv = int(sv)
        except: sv = 0
        u_arr = np.array(src_arr)
        l_arr = np.array(src_arr)
        u_arr[u_arr <= sv] = nd
        l_arr[l_arr >= sv] = nd
        return(u_arr, l_arr)
    
    dst_upper = os.path.join(os.path.dirname(src_dem), '{}_u.tif'.format(os.path.basename(src_dem)[:-4]))
    dst_lower = os.path.join(os.path.dirname(src_dem), '{}_l.tif'.format(os.path.basename(src_dem)[:-4]))
    try:
        src_ds = gdal.Open(src_dem)
    except:
        src_ds = None
    if src_ds is not None:
        src_config = gather_infos(src_ds)
        dst_config = copy_infos(src_config)
        dst_config['fmt'] = 'GTiff'
        ds_arr = src_ds.GetRasterBand(1).ReadAsArray(0, 0, src_config['nx'], src_config['ny'])
        ua, la = np_split(ds_arr, split_value, src_config['ndv'])
        utils.gdal_write(ua, dst_upper, dst_config)
        utils.gdal_write(la, dst_lower, dst_config)
        ua = la = ds_arr = src_ds = None
        src_ds = None
        return([dst_upper, dst_lower])
    else: return(None)

def cut(src_dem, src_region, dst_dem, node='pixel', mode=None, verbose=True):
    """cut src_fn gdal file to srcwin and output dst_fn gdal file

    returns [output-dem, status-code]
    """

    if mode == 'translate':

        src_infos = infos(src_dem)
        this_srcwin = src_region.srcwin(
            geo_transform=src_infos['geoT'],
            x_count=src_infos['nx'],
            y_count=src_infos['ny']
        )
        
        return(utils.run_cmd('gdal_translate {} {} -srcwin {} {} {} {}'.format(
            src_dem, dst_dem, this_srcwin[0], this_srcwin[1], this_srcwin[2], this_srcwin[3]
        ), verbose=verbose))
    else:
        try:
            ds = gdal.Open(src_dem)
        except: ds = None
        if ds is not None:
            ds_config = gather_infos(ds)
            gt = ds_config['geoT']
            srcwin = src_region.srcwin(gt, ds.RasterXSize, ds.RasterYSize, node=node)
            #print(srcwin)
            ds_arr = ds.GetRasterBand(1).ReadAsArray(srcwin[0], srcwin[1], srcwin[2], srcwin[3])
            #print(ds_arr.shape)
            dst_gt = (gt[0] + (srcwin[0] * gt[1]), gt[1], 0., gt[3] + (srcwin[1] * gt[5]), 0., gt[5])
            out_ds_config = set_infos(srcwin[2], srcwin[3], srcwin[2] * srcwin[3], dst_gt,
                                      ds_config['proj'], ds_config['dt'], ds_config['ndv'],
                                      ds_config['fmt'])
            ds = None
            return(utils.gdal_write(ds_arr, dst_dem, out_ds_config, verbose=verbose))
        else: return(None, -1)

def clip(src_dem, dst_dem, src_ply=None, invert=False, verbose=True):
    """clip dem to polygon `src_ply`, optionally invert the clip.

    returns [gdal_raserize-output, gdal_rasterize-return-code]
    """

    gi = infos(src_dem)
    g_region = regions.Region().from_geo_transform(geo_transform=gi['geoT'], x_count=gi['nx'], y_count=gi['ny'])
    tmp_ply = '__tmp_clp_ply.shp'
    
    out, status = utils.run_cmd('ogr2ogr {} {} -clipsrc {} -nlt POLYGON -skipfailures'.format(tmp_ply, src_ply, g_region.format('ul_lr')), verbose=verbose)
    if gi is not None and src_ply is not None:
        #if invert:
        #    gr_cmd = 'gdalwarp -cutline {} -cl {} {} {}'.format(tmp_ply, os.path.basename(tmp_ply).split('.')[0], src_dem, dst_dem)
        #    out, status = utils.run_cmd(gr_cmd, verbose=verbose)
        #else:
        shutil.copyfile(src_dem, dst_dem)
        gr_cmd = 'gdal_rasterize -burn {} -l {} {} {}{}'\
            .format(gi['ndv'], os.path.basename(tmp_ply).split('.')[0], tmp_ply, dst_dem, ' -i' if invert else '')
        out, status = utils.run_cmd(gr_cmd, verbose=verbose)
        utils.remove_glob('__tmp_clp_ply.*')
    else: return(None, -1)
    return(out, status)

def crop(src_dem, dst_dem):

    try:
        ds = gdal.Open(src_dem)
    except: ds = None
    if ds is not None:
        
        ds_config = gather_infos(ds)
        ds_arr = ds.GetRasterBand(1).ReadAsArray()

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
        ds = None
        
        return(utils.gdal_write(dst_arr, dst_dem, ds_config))
    else:
        return(None, -1)

def has_nodata_p(src_gdal):
    try:
        ds = gdal.Open('{}'.format(src_gdal))
    except: ds = None
    if ds is not None:
        ds_config = gather_infos(ds)
        ds_arr = ds.GetRasterBand(1).ReadAsArray()
        ds_arr[ds_arr == ds_config['ndv']] = np.nan

        if np.any(np.isnan(ds_arr)):
            return(True)
    return(False)
    
def polygonize(src_gdal, dst_layer, verbose=False):
    '''run gdal.Polygonize on src_ds and add polygon to dst_layer'''

    try:
        ds = gdal.Open('{}'.format(src_gdal))
    except: ds = None
    if ds is not None:
        ds_arr = ds.GetRasterBand(1)
        if verbose:
            utils.echo_msg('polygonizing {}...'.format(src_gdal))
        status = gdal.Polygonize(
            ds_arr, None, dst_layer, 0,
            callback=gdal.TermProgress if verbose else None
        )
        if verbose:
            utils.echo_msg('polygonized {}'.format(src_gdal))
        ds = ds_arr = None
        return(0, 0)
    else: return(-1, -1)
    
def sample_warp(
        src_dem, dst_dem, x_sample_inc, y_sample_inc,
        src_srs=None, dst_srs=None, src_region=None, sample_alg='bilinear',
        ndv=-9999, tap=False, size=False, verbose=False
):
    
    if size:
        xcount, ycount, dst_gt = src_region.geo_transform(
            x_inc=x_sample_inc, y_inc=y_sample_inc, node='pixel'
        )
        x_sample_inc = y_sample_inc = None
    else:
        xcount = ycount = None

    if src_region is not None:
        out_region = [src_region.xmin, src_region.ymin, src_region.xmax, src_region.ymax]
    else: 
        out_region = None

    if verbose:
        utils.echo_msg(
            'sampling DEM: {} to region: {} and increments: {}/{} using {} from: {} to {}'.format(
                src_dem, out_region, x_sample_inc, y_sample_inc, sample_alg, src_srs, dst_srs
            )
        )

    dst_ds = gdal.Warp('' if dst_dem is None else dst_dem, src_dem, format='MEM' if dst_dem is None else 'GTiff',
                       xRes=x_sample_inc, yRes=y_sample_inc, targetAlignedPixels=tap, width=xcount, height=ycount,
                       dstNodata=ndv, outputBounds=out_region, outputBoundsSRS=dst_srs, resampleAlg=sample_alg, errorThreshold=0,
                       options=["COMPRESS=LZW", "TILED=YES"], srcSRS=src_srs, dstSRS=dst_srs, outputType=gdal.GDT_Float32,
                       callback=gdal.TermProgress if verbose else None)

    #utils.echo_msg('gdalwarp -s_srs {} -t_srs {} -tr {} {} -r bilinear'.format(src_srs, dst_srs, x_sample_inc, y_sample_inc))

    if dst_dem is None:
        return(dst_ds, 0)
    else:
        dst_ds = None
        return(dst_dem, 0)    
    
def chunks(src_dem, n_chunk):
    """split `src_fn` GDAL file into chunks with `n_chunk` cells squared.

    returns a list of chunked filenames.
    """
    
    o_chunks = []
    try:
        src_ds = gdal.Open(src_dem)
    except: src_ds = None
    if src_ds is not None:
        ds_config = gather_infos(src_ds)
        band = src_ds.GetRasterBand(1)
        gt = ds_config['geoT']
        gt = list(gt)
        gt[0] = gt[0] - (gt[1]/2)
        gt[3] = gt[3] - (gt[5]/2)
        gt = tuple(gt)
        
        c_n = 0
        for srcwin in yield_srcwin(src_dem, n_chunk = n_chunk, step = n_chunk):
            this_geo_x_origin, this_geo_y_origin = utils._pixel2geo(srcwin[0], srcwin[1], gt)
            dst_gt = [this_geo_x_origin, float(gt[1]), 0.0, this_geo_y_origin, 0.0, float(gt[5])]
            band_data = band.ReadAsArray(srcwin[0], srcwin[1], srcwin[2], srcwin[3])
            if not np.all(band_data == band_data[0,:]):
                dst_config = copy_infos(ds_config)
                dst_config['nx'] = srcwin[2]
                dst_config['ny'] = srcwin[3]
                dst_config['geoT'] = dst_gt
                this_region = regions.Region().from_geo_transform(dst_gt, dst_config['nx'], dst_config['ny'])
                o_chunk = '{}_chnk{}.tif'.format(os.path.basename(src_dem).split('.')[0], c_n)
                dst_fn = os.path.join(os.path.dirname(src_dem), o_chunk)
                o_chunks.append(dst_fn)
                utils.gdal_write(band_data, dst_fn, dst_config)
                c_n += 1                
    return(o_chunks)

def yield_srcwin(src_dem, n_chunk=10, step=5, verbose=False):
    """yield source windows in n_chunks at step"""
    
    ds_config = infos(src_dem)
    gt = ds_config['geoT']
    x_chunk = n_chunk
    y_chunk = 0
    i_chunk = 0
    x_i_chunk = 0
    
    if verbose:
        _prog = utils.CliProgress(
           'parsing dataset {}'.format(src_dem)
        )
    while True:
        y_chunk = n_chunk
        if verbose:
           _prog.update_perc((i_chunk*step, ds_config['nb']/step))
        
        #utils.echo_msg_inline('[{:.2f}%]'.format(((i_chunk*step)/(ds_config['nb']/step) * 100)))
        while True:
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
        if x_chunk > ds_config['nx']:
            break
        else:
            x_chunk += step
            x_i_chunk += 1
    #utils.echo_msg('[OK]'.format(((i_chunk*step)/(ds_config['nb']/step) * 100)))
    if verbose:
       _prog.end(0, 'parsed dataset {}'.format(src_dem))

def blur(src_dem, dst_dem, sf = 1):
    """gaussian blur on src_dem using a smooth-factor of `sf`
    runs np_gaussian_blur(ds.Array, sf)
    """

    def np_gaussian_blur(in_array, size):
        """blur an array using fftconvolve from scipy.signal
        size is the blurring scale-factor.

        returns the blurred array
        """

        from scipy.signal import fftconvolve
        from scipy.signal import convolve
        padded_array = np.pad(in_array, size, 'symmetric')
        x, y = np.mgrid[-size:size + 1, -size:size + 1]
        g = np.exp(-(x**2 / float(size) + y**2 / float(size)))
        g = (g / g.sum()).astype(in_array.dtype)
        in_array = None
        out_array = fftconvolve(padded_array, g, mode = 'valid')
        return(out_array)
    
    try:
        ds = gdal.Open(src_dem)
    except: ds = None
    
    if ds is not None:
        ds_config = gather_infos(ds)
        ds_array = ds.GetRasterBand(1).ReadAsArray(0, 0, ds_config['nx'], ds_config['ny'])
        ds = None
        msk_array = np.array(ds_array)
        msk_array[msk_array == ds_config['ndv']] = np.nan
        msk_array[~np.isnan(msk_array)] = 1
        #msk_array[np.isnan(msk_array)] = 0
        #msk_array[msk_array != ds_config['ndv']] = 1
        ds_array[np.isnan(msk_array)] = 0
        ds_array[ds_array == ds_config['ndv']] = 0
        smooth_array = np_gaussian_blur(ds_array, int(sf))
        smooth_array = smooth_array * msk_array
        mask_array = ds_array = None
        smooth_array[np.isnan(smooth_array)] = ds_config['ndv']
        return(utils.gdal_write(smooth_array, dst_dem, ds_config))
    else: return([], -1)

def filter_outliers(src_gdal, dst_gdal, threshhold=None, chunk_size=None, chunk_step=None, replace=True):
    """scan a src_gdal file for outliers and remove them"""
    
    try:
        ds = gdal.Open(src_gdal)
    except: ds = None

    if ds is not None:
        tnd = 0
        
        ds_config = gather_infos(ds)
        ds_band = ds.GetRasterBand(1)
        ds_array = ds_band.ReadAsArray(0, 0, ds_config['nx'], ds_config['ny'])
        gt = ds_config['geoT']
        if threshhold is None:
            ds_std = np.std(ds_array)
        else: ds_std = utils.float_or(threshhold)

        driver = gdal.GetDriverByName('MEM')
        mem_ds = driver.Create('tmp', ds_config['nx'], ds_config['ny'], 1, ds_config['dt'])
        mem_ds.SetGeoTransform(gt)
        mem_ds.SetProjection(ds_config['proj'])
        band = mem_ds.GetRasterBand(1)
        band.SetNoDataValue(ds_config['ndv'])
        band.WriteArray(ds_array)

        ds = None
        if chunk_size is None:
            n_chunk = int(ds_config['nx'] * .005)
            n_chunk = 10 if n_chunk < 10 else n_chunk
        else: n_chunk = chunk_size
        if chunk_step is None:
            n_step = int(n_chunk/2)
        else: n_step = chunk_step

        utils.echo_msg('scanning {} for spikes with {}@{} MAX {}...'.format(src_gdal, n_chunk, n_step, ds_std))
        for srcwin in yield_srcwin(src_gdal, n_chunk = n_chunk, step = n_step):
            band_data = band.ReadAsArray(srcwin[0], srcwin[1], srcwin[2], srcwin[3])
            band_data[band_data == ds_config['ndv']] = np.nan
            this_geo_x_origin, this_geo_y_origin = utils._pixel2geo(srcwin[0], srcwin[1], gt)
            dst_gt = [this_geo_x_origin, float(gt[1]), 0.0, this_geo_y_origin, 0.0, float(gt[5])]

            dst_config = copy_infos(ds_config)
            dst_config['nx'] = srcwin[2]
            dst_config['ny'] = srcwin[3]
            dst_config['geoT'] = dst_gt
            
            if not np.all(band_data == band_data[0,:]):
                while True:
                    nd = 0
                    srcwin_std = np.nanstd(band_data)
                    utils.echo_msg(srcwin_std)
                    if srcwin_std < ds_std: break
                    
                    srcwin_perc75 = np.nanpercentile(band_data, 75)
                    srcwin_perc25 = np.nanpercentile(band_data, 25)
                    iqr_p = (srcwin_perc75 - srcwin_perc25) * 1.5
                    upper_limit = srcwin_perc75 + iqr_p
                    lower_limit = srcwin_perc25 - iqr_p

                    for i in range(0, srcwin[2]):
                        for j in range(0, srcwin[3]):
                            bandz = band_data[j][i]
                            if bandz > upper_limit or bandz < lower_limit:
                                ds_array[j+srcwin[1]][i+srcwin[0]] = ds_config['ndv']
                                band.WriteArray(ds_array)
                                band_data[j][i] = np.nan
                                nd += 1
                    tnd += nd
                    if nd == 0: break
                band_data = None
                
        utils.echo_msg('filtering {} spikes...'.format(tnd))
        if tnd > 0 and replace:
            driver = gdal.GetDriverByName('MEM')
            tmp_ds = driver.Create('tmp', ds_config['nx'], ds_config['ny'], 1, ds_config['dt'])
            tmp_ds.SetGeoTransform(ds_config['geoT'])
            tmp_ds.SetProjection(ds_config['proj'])
            ds_band = tmp_ds.GetRasterBand(1)
            ds_band.SetNoDataValue(ds_config['ndv'])
            ds_band.WriteArray(ds_array)
            result = gdal.FillNodata(targetBand=ds_band, maskBand=None, maxSearchDist=100,
                                     smoothingIterations=4, callback=gdal.TermProgress)
        
            ds_array = ds_band.ReadAsArray()
            tmp_ds = None
        else:
            ds_array = ds_band.ReadAsArray()

        out, status = utils.gdal_write(ds_array, dst_gdal, ds_config)
        mem_ds = None
        return(out, status)
    else: return(None)    
    
def filter_outliers_slp(src_dem, dst_dem, chunk_size=None, chunk_step=None, agg_level=1, replace=True):
    """scan a src_dem file for outliers and remove them
    
    aggressiveness depends on the outlier percentiles and the chunk_size/step; 75/25 is default 
    for statistical outlier discovery, 55/45 will be more aggressive, etc. Using a large chunk size 
    will filter more cells and find potentially more or less outliers depending on the data.
    agg_level is 1 to 9
    """

    import scipy

    agg_level = utils.int_or(agg_level)

    if agg_level is not None:
        percentile = 100 - (agg_level*10)
        if percentile < 50: percentle = 50
        if percentile > 100: percentile = 100
        #curv_percentile = 0 + (agg_level*10)

    percentile = agg_level
        
    # if agg_level == 1:
    #     percentile = 95
    #     curv_percentile = 95
    # elif agg_level == 2:
    #     percentile = 85
    #     curv_percentile = 95
    # elif agg_level == 3:
    #     percentile = 75
    #     curv_percentile = 95
    # elif agg_level == 4:
    #     percentile = 65
    #     curv_percentile = 95
    # elif agg_level == 5:
    #     percentile = 55
    #     curv_percentile = 95
    # elif agg_level == 6:
    #     percentile = 55
    #     curv_percentile = 85
    # elif agg_level == 7:
    #     percentile = 55
    #     curv_percentile = 75
    # elif agg_level == 8:
    #     percentile = 55
    #     curv_percentile = 65
    # elif agg_level == 9:
    #     percentile = 55
    #     curv_percentile = 55
    # else:
    #     percentile = 75
    #     curv_percentile = 95

        
    max_percentile = percentile
    min_percentile = 100 - percentile
    #curv_max_percentile = percentile
    #curv_min_percentile = 100 - percentile
    
    try:
        ds = gdal.Open(src_dem)
    except: ds = None

    if ds is not None:
        tnd = 0        
        ds_config = gather_infos(ds)
        ds_band = ds.GetRasterBand(1)
        #ds_array = ds_band.ReadAsArray(0, 0, ds_config['nx'], ds_config['ny'])
        gt = ds_config['geoT']
        gdt = gdal.GDT_Float32
        ndv = ds_band.GetNoDataValue()

        driver = gdal.GetDriverByName('GTiff')
        dst_ds = driver.Create(dst_dem, ds.RasterXSize, ds.RasterYSize, 1, gdt,
                               options=['COMPRESS=LZW', 'PREDICTOR=2', 'TILED=YES'])
        dst_ds.SetGeoTransform(gt)
        dst_band = dst_ds.GetRasterBand(1)
        dst_band.SetNoDataValue(ndv)
        
        #ds = None
        if chunk_size is None:
            n_chunk = int(ds_config['nx'] * .05)
            n_chunk = 10 if n_chunk < 10 else n_chunk
        else:
            n_chunk = chunk_size

        chunk_step = utils.int_or(chunk_step)
        n_step = chunk_step if chunk_step is not None else int(n_chunk)
        #n_step = n_chunk
        n_step = n_chunk/4
        
        utils.echo_msg(
            'scanning {} for outliers with {}@{} using aggression level {} ({}/{})...'.format(
                src_dem, n_chunk, n_step, agg_level, max_percentile, min_percentile
            )
        )

        for srcwin in utils.yield_srcwin(
                (ds.RasterYSize, ds.RasterXSize), n_chunk = n_chunk, step = n_step, verbose=True
        ):
            nd = 0
            band_data = ds_band.ReadAsArray(*srcwin)
            band_data[band_data == ds_config['ndv']] = np.nan

            if np.all(np.isnan(band_data)):
                #dst_band.WriteArray(band_data, srcwin[0], srcwin[1])
                continue
            
            coverage_data = ds_band.ReadAsArray(*srcwin)
            coverage_data[coverage_data == ds_config['ndv']] = np.nan
                        
            this_geo_x_origin, this_geo_y_origin = utils._pixel2geo(srcwin[0], srcwin[1], gt)
            dst_gt = [this_geo_x_origin, float(gt[1]), 0.0, this_geo_y_origin, 0.0, float(gt[5])]
            
            dst_config = copy_infos(ds_config)
            dst_config['nx'] = srcwin[2]
            dst_config['ny'] = srcwin[3]
            dst_config['geoT'] = dst_gt

            if not np.all(band_data == band_data[0,:]):
                
                px, py = np.gradient(band_data, gt[1])
                slp_data = np.sqrt(px ** 2, py ** 2)
                #slp_data = np.degrees(np.arctan(slp_data_))
                
                px, py = np.gradient(slp_data, gt[1])
                curv_data = np.sqrt(px ** 2, py ** 2)
                #curv_data = np.degrees(np.arctan(curv_data_))

                srcwin_perc75 = np.nanpercentile(band_data, max_percentile)
                srcwin_perc25 = np.nanpercentile(band_data, min_percentile)
                iqr_p = (srcwin_perc75 - srcwin_perc25) * 1.5
                upper_limit = srcwin_perc75 + iqr_p
                lower_limit = srcwin_perc25 - iqr_p
                #utils.echo_msg('elev upper limit: {}'.format(upper_limit))
                #utils.echo_msg('elev lower limit: {}'.format(lower_limit))
                
                slp_srcwin_perc75 = np.nanpercentile(slp_data, max_percentile)
                slp_srcwin_perc25 = np.nanpercentile(slp_data, min_percentile)
                slp_iqr_p = (slp_srcwin_perc75 - slp_srcwin_perc25) * 1.5
                slp_upper_limit = slp_srcwin_perc75 + slp_iqr_p
                slp_lower_limit = slp_srcwin_perc25 - slp_iqr_p
                #utils.echo_msg('slp upper limit: {}'.format(slp_upper_limit))
                #utils.echo_msg('slp lower limit: {}'.format(slp_lower_limit))

                curv_srcwin_perc75 = np.nanpercentile(curv_data, max_percentile)
                curv_srcwin_perc25 = np.nanpercentile(curv_data, min_percentile)
                curv_iqr_p = (curv_srcwin_perc75 - curv_srcwin_perc25) * 1.5
                curv_upper_limit = curv_srcwin_perc75 + curv_iqr_p
                curv_lower_limit = curv_srcwin_perc25 - curv_iqr_p
                #utils.echo_msg('curv upper limit: {}'.format(curv_upper_limit))
                #utils.echo_msg('curv lower limit: {}'.format(curv_lower_limit))
                #print(curv_data)
                
                band_data[(curv_data > curv_upper_limit)] = np.nan
                #band_data[((band_data > upper_limit) | (band_data < lower_limit))] = np.nan
                # band_data[((band_data > upper_limit) | (band_data < lower_limit)) \
                #           & ((curv_data > curv_upper_limit) | (curv_data < curv_lower_limit)) \
                #           & ((slp_data > slp_upper_limit) | (slp_data < slp_lower_limit))] = np.nan
                #band_data[((band_data > upper_limit) | (band_data < lower_limit)) \
                #          & (curv_data > curv_upper_limit)] = np.nan
                #(curv_data > curv_upper_limit)] = np.nan
                
                ## fill nodata here...
                if replace:
                    point_indices = np.nonzero(~np.isnan(band_data))
                    if len(point_indices[0]):
                        point_values = band_data[point_indices]
                        xi, yi = np.mgrid[0:srcwin[3], 0:srcwin[2]]

                        try:
                            interp_data = scipy.interpolate.griddata(
                                np.transpose(point_indices), point_values,
                                (xi, yi), method='cubic'
                            )
                            interp_data[np.isnan(coverage_data)] = np.nan
                            dst_band.WriteArray(interp_data, srcwin[0], srcwin[1])
                        except:
                            dst_band.WriteArray(band_data, srcwin[0], srcwin[1])
                else:
                    dst_band.WriteArray(band_data, srcwin[0], srcwin[1])
                    #dst_band.WriteArray(curv_data, srcwin[0], srcwin[1])

        ds = dst_ds = None
        return(dst_dem, 0)
    else:
        return(None, -1)

def filter_outliers_from_array(ds_array, ds_config, chunk_size=None, chunk_step=None, agg_level=1, replace=True):
    """    
    aggressiveness depends on the outlier percentiles and the chunk_size/step; 75/25 is default 
    for statistical outlier discovery, 55/45 will be more aggressive, etc. Using a large chunk size 
    will filter more cells and find potentially more or less outliers depending on the data.
    agg_level is 1 to 5
    """

    def yield_srcwin_(ds_config, n_chunk=10, step=5, verbose=False):
        """yield source windows in n_chunks at step"""

        gt = ds_config['geoT']
        x_chunk = n_chunk
        y_chunk = 0
        i_chunk = 0
        x_i_chunk = 0

        while True:
            y_chunk = n_chunk
            utils.echo_msg_inline('[{:.2f}%]'.format(((i_chunk*step)/(ds_config['nb']/step) * 100)))
            while True:
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
            if x_chunk > ds_config['nx']:
                break
            else:
                x_chunk += step
                x_i_chunk += 1
    
    agg_level = utils.int_or(agg_level)

    if agg_level == 1:
        percentile = 95
    elif agg_level == 2:
        percentile = 85
    elif agg_level == 3:
        percentile = 75
    elif agg_level == 4:
        percentile = 65
    elif agg_level == 5:
        percentile = 55
    else:
        percentile = 75
        
    max_percentile = percentile
    min_percentile = 100 - percentile

    tnd = 0        
    gt = ds_config['geoT']

    if chunk_size is None:
        n_chunk = int(ds_config['nx'] * .05)
        n_chunk = 10 if n_chunk < 10 else n_chunk
    else: n_chunk = chunk_size
    if chunk_step is None:
        n_step = int(n_chunk/4)
    else: n_step = chunk_step

    utils.echo_msg(
        'scanning for outliers with {}@{} using aggression level {} ({}/{})...'.format(
            n_chunk, n_step, agg_level, max_percentile, min_percentile
        )
    )    
    for srcwin in yield_srcwin(ds_config, n_chunk = n_chunk, step = n_step):
        nd = 0
        band_data = ds_array[srcwin[1]:srcwin[1]+srcwin[3],srcwin[0]:srcwin[0]+srcwin[2]]
        band_data[band_data == ds_config['ndv']] = np.nan
        this_geo_x_origin, this_geo_y_origin = utils._pixel2geo(srcwin[0], srcwin[1], gt)
        dst_gt = [this_geo_x_origin, float(gt[1]), 0.0, this_geo_y_origin, 0.0, float(gt[5])]
        dst_config = copy_infos(ds_config)
        dst_config['nx'] = srcwin[2]
        dst_config['ny'] = srcwin[3]
        dst_config['geoT'] = dst_gt

        if np.all(np.isnan(band_data)):
            continue

        if not np.all(band_data == band_data[0,:]):
            px, py = np.gradient(band_data, gt[1])
            slp_data_ = np.sqrt(px ** 2, py ** 2)
            slp_data = np.degrees(np.arctan(slp_data_))

            px, py = np.gradient(slp_data, gt[1])
            curv_data_ = np.sqrt(px ** 2, py ** 2)
            curv_data = np.degrees(np.arctan(curv_data_))

            srcwin_perc75 = np.nanpercentile(band_data, max_percentile)
            srcwin_perc25 = np.nanpercentile(band_data, min_percentile)
            iqr_p = (srcwin_perc75 - srcwin_perc25) * 1.5
            upper_limit = srcwin_perc75 + iqr_p
            lower_limit = srcwin_perc25 - iqr_p

            slp_srcwin_perc75 = np.nanpercentile(slp_data, max_percentile)
            slp_srcwin_perc25 = np.nanpercentile(slp_data, min_percentile)
            slp_iqr_p = (slp_srcwin_perc75 - slp_srcwin_perc25) * 1.5
            slp_upper_limit = slp_srcwin_perc75 + slp_iqr_p
            slp_lower_limit = slp_srcwin_perc25 - slp_iqr_p

            curv_srcwin_perc75 = np.nanpercentile(curv_data, 95)
            curv_srcwin_perc25 = np.nanpercentile(curv_data, 5)
            curv_iqr_p = (curv_srcwin_perc75 - curv_srcwin_perc25) * 1.5
            curv_upper_limit = curv_srcwin_perc75 + curv_iqr_p
            curv_lower_limit = curv_srcwin_perc25 - curv_iqr_p

            band_data[((band_data > upper_limit) | (band_data < lower_limit)) & ((curv_data > curv_upper_limit) | (curv_data < curv_lower_limit))] = ds_config['ndv']
            ds_array[srcwin[1]:srcwin[1]+srcwin[3],srcwin[0]:srcwin[0]+srcwin[2]] = band_data
            tnd += len(band_data[band_data == ds_config['ndv']])

    utils.echo_msg('filtered {} outliers...'.format(tnd))
    return(ds_array, 0)
    
def filter_outliers_slp_(src_dem, dst_dem, chunk_size=None, chunk_step=None, percentile=75):
    """scan a src_dem file for outliers, remove and fill them"""

    try:
        ds = gdal.Open(src_dem)
    except: ds = None

    if ds is not None:
        tnd = 0

        ds_config = gather_infos(ds)
        ds_band = ds.GetRasterBand(1)
        ds_array = ds_band.ReadAsArray(0, 0, ds_config['nx'], ds_config['ny'])
        gt = ds_config['geoT']
        ds_std = np.std(ds_array)
        #slp_std = ds_std
        last_std = ds_std
        driver = gdal.GetDriverByName('MEM')
        mem_ds = driver.Create('tmp', ds_config['nx'], ds_config['ny'], 1, ds_config['dt'])
        mem_ds.SetGeoTransform(gt)
        mem_ds.SetProjection(ds_config['proj'])
        band = mem_ds.GetRasterBand(1)
        band.SetNoDataValue(ds_config['ndv'])
        band.WriteArray(ds_array)

        ds = None
        if chunk_size is None:
            n_chunk = int(ds_config['nx'] * .05)
            n_chunk = 10 if n_chunk < 10 else n_chunk
        else: n_chunk = chunk_size
        if chunk_step is None:
            n_step = int(n_chunk/4)
        else: n_step = chunk_step

        utils.echo_msg('scanning {} for spikes with {}@{}...'.format(src_dem, n_chunk, n_step))
        for srcwin in yield_srcwin(src_dem, n_chunk = n_chunk, step = n_step):
            band_data = band.ReadAsArray(srcwin[0], srcwin[1], srcwin[2], srcwin[3])
            band_data[band_data == ds_config['ndv']] = np.nan

            if np.all(np.isnan(band_data)):
                continue
            this_geo_x_origin, this_geo_y_origin = utils._pixel2geo(srcwin[0], srcwin[1], gt)
            dst_gt = [this_geo_x_origin, float(gt[1]), 0.0, this_geo_y_origin, 0.0, float(gt[5])]

            dst_config = copy_infos(ds_config)
            dst_config['nx'] = srcwin[2]
            dst_config['ny'] = srcwin[3]
            dst_config['geoT'] = dst_gt
            last_std = last_slp_std = None
            
            if not np.all(band_data == band_data[0,:]):
                while True:
                    nd = 0                    
                    srcwin_std = np.nanstd(band_data)                    
                    px, py = np.gradient(band_data, gt[1])
                    slp_data_ = np.sqrt(px ** 2, py ** 2)
                    slp_data = np.degrees(np.arctan(slp_data_))
                    slp_srcwin_std = np.nanstd(slp_data)
                    
                    if (last_std is None or srcwin_std < last_std) and (last_slp_std is None or slp_srcwin_std < last_slp_std):
                        last_std = srcwin_std
                        last_slp_std = slp_srcwin_std
                    else:
                        break

                    srcwin_perc75 = np.nanpercentile(band_data, 75)
                    srcwin_perc25 = np.nanpercentile(band_data, 25)
                    iqr_p = (srcwin_perc75 - srcwin_perc25) * 1.5
                    upper_limit = srcwin_perc75 + iqr_p
                    lower_limit = srcwin_perc25 - iqr_p

                    slp_srcwin_perc75 = np.nanpercentile(slp_data, 75)
                    slp_srcwin_perc25 = np.nanpercentile(slp_data, 25)
                    slp_iqr_p = (slp_srcwin_perc75 - slp_srcwin_perc25) * 1.5
                    slp_upper_limit = slp_srcwin_perc75 + slp_iqr_p
                    slp_lower_limit = slp_srcwin_perc25 - slp_iqr_p

                    for i in range(0, srcwin[2]):
                        for j in range(0, srcwin[3]):
                            bandz = band_data[j][i]
                            slpz = slp_data[j][i]
                            if bandz > upper_limit or bandz < lower_limit:
                                if slpz > slp_upper_limit or slpz < slp_lower_limit:
                                    band_data[j][i] = ds_config['ndv']
                                    nd += 1
                    tnd += nd
                    if nd == 0:
                        break
                    else:
                        driver = gdal.GetDriverByName('MEM')
                        tmp_ds = driver.Create('tmp', srcwin[2], srcwin[3], 1, ds_config['dt'])
                        tmp_ds.SetGeoTransform(dst_gt)
                        tmp_ds.SetProjection(ds_config['proj'])
                        tmp_band = tmp_ds.GetRasterBand(1)
                        tmp_band.SetNoDataValue(ds_config['ndv'])
                        tmp_band.WriteArray(band_data)
                        result = gdal.FillNodata(
                            targetBand=tmp_band,
                            maskBand=None,
                            maxSearchDist=chunk_size,
                            smoothingIterations=4,
                            callback=None
                        )

                        band_data = tmp_band.ReadAsArray()
                        tmp_ds = None
                        band.WriteArray(band_data, xoff=srcwin[0], yoff=srcwin[1])
                band_data = slp_data = None
    else:
        return(None, -1)
    
    utils.echo_msg('filtered {} cells from {} and writing {} to file...'.format(tnd, src_dem, dst_dem))
    ds_array = band.ReadAsArray()
    out, status = utils.gdal_write(ds_array, dst_dem, ds_config)
    mem_ds = None
    return(out, status)

def grdfilter(src_dem, dst_dem, dist='c3s', node='pixel', verbose=False):
    """filter `src_dem` using GMT grdfilter

    Args:
      src_grd (str): pathname to a source grid file
      dst_grd (str): pathname to a destination grid file
      dist (str): a GMT string increment to use in filter
      node (str): `pixel` or `grid`; the grid-node
      verbose (bool): increase verbosity

    Returns:
      list: [cmd-output, cmd-return-code]
    """
    
    #ft_cmd1 = ('gmt grdfilter -V {} -G{} -R{} -Fc{} -D1{}'.format(src_dem, dst_dem, src_dem, dist, ' -rp' if node == 'pixel' else ''))
    ft_cmd1 = ('gmt grdfilter -V {} -G{} -F{} -D1{}'.format(src_dem, dst_dem, dist, ' -rp' if node == 'pixel' else ''))
    return(utils.run_cmd(ft_cmd1, verbose=verbose))

def filter_old(src_dem, dst_dem, fltr=1, fltr_val=None, split_val=None, mask=None, node='pixel'):
    """filter raster using smoothing factor `fltr`; optionally
    only smooth bathymetry (sub-zero) using a split_val of 0.

    Args: 
      fltr (int): the filter to use, 1, 2 or 3
      flt_val (varies): the filter value, varies by filter.
      split_val (float): an elevation value (only filter below this value)

    Returns:
      0 for success or -1 for failure
    """

    utils.echo_msg('filtering DEM {} using {}@{}'.format(src_dem, fltr, fltr_val))

    if os.path.exists(src_dem):
        if split_val is not None:
            dem_u, dem_l = split(src_dem, split_val)
        else: dem_l = src_dem

        if int(fltr) == 1:
            out, status = blur(
                dem_l, 'tmp_fltr.tif', fltr_val if fltr_val is not None else 10)
        elif int(fltr) == 2:
            out, status = grdfilter(
                dem_l, 'tmp_fltr.tif=gd:GTiff', dist = fltr_val if fltr_val is not None else '1s',
                node = node, verbose = True)
        elif int(fltr) == 3:
            out, status = filter_outliers_slp(
                dem_l, 'tmp_fltr.tif', fltr_val if fltr_val is not None else 10)
        else:
            utils.echo_warning_msg('invalid filter {}, defaulting to blur'.format(fltr))
            out, status = blur(dem_l, 'tmp_fltr.tif', fltr_val if utils.int_or(fltr_val) is not None else 10)
        if status != 0: return(status)

        if split_val is not None:
            try:
                ds = gdal.Open(src_dem)
            except:
                ds = None

            if ds is not None:
                ds_config = gather_infos(ds)
                msk_arr = ds.GetRasterBand(1).ReadAsArray()
                msk_arr[msk_arr == ds_config['ndv']] = np.nan
                msk_arr[~np.isnan(msk_arr)] = 1
                
                ds = None
                try:
                    u_ds = gdal.Open(dem_u)
                except:
                    u_ds = None
                    
                if u_ds is not None:
                    try:
                        l_ds = gdal.Open('tmp_fltr.tif')
                    except:
                        l_ds = None
                        
                    if l_ds is not None:
                        u_arr = u_ds.GetRasterBand(1).ReadAsArray()
                        l_arr = l_ds.GetRasterBand(1).ReadAsArray()
                        u_arr[u_arr == ds_config['ndv']] = np.nan
                        u_arr[np.isnan(u_arr)] = 0
                        l_arr[l_arr == ds_config['ndv']] = np.nan
                        l_arr[np.isnan(l_arr)] = 0
                        ds_arr = (u_arr + l_arr) * msk_arr
                        ds_arr[np.isnan(ds_arr)] = ds_config['ndv']
                        utils.gdal_write(ds_arr, '__merged.tif', ds_config)
                        l_ds = None
                        utils.remove_glob(dem_l)
                    u_ds = None
                    utils.remove_glob(dem_u)
                os.rename('__merged.tif', dst_dem)
                utils.remove_glob('tmp_fltr.tif')
        else:
            os.rename('tmp_fltr.tif', dst_dem)
        return(0)
    else: return(-1)

def filter_(src_dem, dst_dem, fltr=1, fltr_val=None, split_val=None, mask=None, node='pixel'):
    """filter raster using smoothing factor `fltr`; optionally
    only smooth bathymetry (sub-zero) using a split_val of 0.

    Args: 
      fltr (int): the filter to use, 1, 2 or 3
      flt_val (varies): the filter value, varies by filter.
      split_val (float): an elevation value (only filter below this value)

    Returns:
      0 for success or -1 for failure
    """

    utils.echo_msg('filtering DEM {} using {}@{}'.format(src_dem, fltr, fltr_val))

    if os.path.exists(src_dem):
        
        if int(fltr) == 1:
            out, status = blur(
                src_dem, 'tmp_fltr.tif', fltr_val if fltr_val is not None else 10)
        elif int(fltr) == 2:
            out, status = grdfilter(
                src_dem, 'tmp_fltr.tif=gd:GTiff', dist = fltr_val if fltr_val is not None else '1s',
                node = node, verbose = True)
        elif int(fltr) == 3:
            out, status = filter_outliers_slp(
                src_dem, 'tmp_fltr.tif', agg_level=fltr_val if fltr_val is not None else 5, replace=True)
        else:
            utils.echo_warning_msg('invalid filter {}, defaulting to blur'.format(fltr))
            out, status = blur(src_dem, 'tmp_fltr.tif', fltr_val if utils.int_or(fltr_val) is not None else 10)
        if status != 0: return(status)

        split_val = utils.float_or(split_val)
        if split_val is not None:
            try:
                ds = gdal.Open(src_dem)
            except:
                ds = None

            if ds is not None:
                ds_config = gather_infos(ds)

                elev_array = ds.GetRasterBand(1).ReadAsArray()
                mask_array = np.zeros((ds_config['ny'],ds_config['nx']))
                
                mask_array[elev_array == ds_config['ndv']] = 0
                mask_array[elev_array < split_val] = 1
                #mask_array[~np.isnan(mask_array)] = 0

                elev_array[elev_array < split_val] = 0

                try:
                    s_ds = gdal.Open('tmp_fltr.tif')
                except:
                    s_ds = None

                if s_ds is not None:
                    s_array = s_ds.GetRasterBand(1).ReadAsArray()
                    s_array = s_array * mask_array
                    smoothed_array = s_array + elev_array
                    elev_array = None

                    utils.gdal_write(smoothed_array, dst_dem, ds_config)
                    s_ds = None
                    
                utils.remove_glob('tmp_fltr.tif')
                ds = None
        else:
            os.rename('tmp_fltr.tif', dst_dem)
        return(0)
    else: return(-1)
    
def slope(src_gdal, dst_gdal, s = 111120):
    """generate a slope grid with GDAL

    return cmd output and status
    """
    
    gds_cmd = 'gdaldem slope {} {} {} -compute_edges'.format(src_gdal, dst_gdal, '' if s is None else '-s {}'.format(s))
    return(utils.run_cmd(gds_cmd))
    
def proximity(src_fn, dst_fn):
    """compute a proximity grid via GDAL

    return 0 if success else None
    """

    prog_func = None
    dst_ds = None
    try:
        src_ds = gdal.Open(src_fn)
    except: src_ds = None
    if src_ds is not None:
        src_band = src_ds.GetRasterBand(1)
        ds_config = gather_infos(src_ds)

        if dst_ds is None:
            drv = gdal.GetDriverByName('GTiff')
            dst_ds = drv.Create(dst_fn, ds_config['nx'], ds_config['ny'], 1, ds_config['dt'], [])
        dst_ds.SetGeoTransform(ds_config['geoT'])
        dst_ds.SetProjection(ds_config['proj'])
        dst_band = dst_ds.GetRasterBand(1)
        dst_band.SetNoDataValue(ds_config['ndv'])
        gdal.ComputeProximity(src_band, dst_band, ['DISTUNITS=PIXEL'], callback = prog_func)
        dst_band = src_band = dst_ds = src_ds = None
        return(0)
    else: return(None)
    
def percentile(src_gdal, perc = 95):
    """calculate the `perc` percentile of src_fn gdal file.

    return the calculated percentile
    """
    
    try:
        ds = gdal.Open(src_gdal)
    except: ds = None
    if ds is not None:
        ds_array = np.array(ds.GetRasterBand(1).ReadAsArray())
        ds_array[ds_array == ds.GetRasterBand(1).GetNoDataValue()] = np.nan
        x_dim = ds_array.shape[0]
        ds_array_flat = ds_array.flatten()
        ds_array = ds_array_flat[ds_array_flat != 0]
        if len(ds_array) > 0:
            p = np.nanpercentile(ds_array, perc)
            #percentile = 2 if p < 2 else p
        else: p = 2
        ds = ds_array = ds_array_flat = None
        return(p)
    else: return(None)

def generate_mem_ds(ds_config, name='MEM'):
    """Create temporary gdal mem dataset"""
        
    mem_driver = gdal.GetDriverByName('MEM')
    mem_ds = mem_driver.Create(name, ds_config['nx'], ds_config['ny'], 1, ds_config['dt'])
    if mem_ds is not None:
        mem_ds.SetGeoTransform(ds_config['geoT'])
        if ds_config['proj'] is not None:
            mem_ds.SetProjection(ds_config['proj'])
            
        mem_band = mem_ds.GetRasterBand(1)
        mem_band.SetNoDataValue(ds_config['ndv'])
        
    return(mem_ds)
    
# def xyz_block_array(src_xyz, region=None, xinc=None, yinc=None, weights=False, min_count=None, out_name=None, verbose=True):
#     """block the src_xyz data to the mean block value
#     src_xyz should be a list or generator of xyzfun.XYZPoint objects.

#     Yields:
#       list: xyz data for each block with data
#     """

#     xcount, ycount, dst_gt = region.geo_transform(
#         x_inc=.xinc, y_inc=.yinc
#     )
#     gdt = gdal.GDT_Float32
#     z_array = np.zeros((ycount, xcount))
#     count_array = np.zeros((ycount, xcount))
#     if weights:
#         weight_array = np.zeros((ycount, xcount))

#     if verbose:
#         utils.echo_msg(
#             'blocking data to {}/{} grid'.format(ycount, xcount)
#         )
#     for this_xyz in src_xyz:
#         #if regions.xyz_in_region_p(this_xyz, self.p_region):
#         #if self.weights:
#         #    this_z = this_xyz.z * this_xyz.w
#         #else:
#         #    this_z = this_xyz.z

#         xpos, ypos = utils._geo2pixel(
#             this_xyz.x, this_xyz.y, dst_gt
#         )
#         try:
#             z_array[ypos, xpos] += this_xyz.z * this_xyz.w
#             count_array[ypos, xpos] += 1
#             if self.weights:
#                 weight_array[ypos, xpos] += this_xyz.w

#         except: pass

#     count_array[count_array == 0] = np.nan
#     if self.weights:
#         weight_array[weight_array == 0] = np.nan
#         weight_array = (weight_array/count_array)
#         z_array = (z_array/weight_array)/count_array
#     else:
#         z_array = (z_array/count_array)
#         weight_array = np.ones((ycount, xcount))

#     if min_count is not None:
#         z_array[count_array < min_count] = np.nan
#         weight_array[count_array < min_count] = np.nan

#     if out_name is not None:
#         ds_config = demfun.set_infos(
#             xcount,
#             ycount,
#             xcount * ycount,
#             dst_gt,
#             self.dst_srs,
#             gdal.GDT_Float32,
#             -9999,
#             'GTiff'
#         )
#         utils.gdal_write(z_array, '{}_n.tif'.format(out_name), ds_config, verbose=True)
#         utils.gdal_write(weight_array, '{}_w.tif'.format(out_name), ds_config, verbose=True)
#         utils.gdal_write(count_array, '{}_c.tif'.format(out_name), ds_config, verbose=True)
#         return('{}_n.tif'.format(out_name), '{}_w.tif'.format(out_name), '{}_c.tif'.format(out_name))
#     else:
#         return(z_array, weight_array, count_array, dst_gt)

# def xyz_block_yield(src_xyz, region=None, xinc=None, yinc=None, weights=False, verbose=True):
#     """block the src_xyz data to the mean block value

#     Yields:
#       list: xyz data for each block with data
#     """

#     z_array, weight_array, count_array, dst_gt = _xyz_block_array(
#         src_xyz, region=region, xinc=xinc, yinc=yinc, weights=weights, verbose=verbose
#     )
#     ycount, xcount = z_array.shape

#     for y in range(0, ycount):
#         for x in range(0, xcount):
#             z = z_array[y,x]
#             if not np.isnan(z):
#                 geo_x, geo_y = utils._pixel2geo(x, y, dst_gt)
#                 out_xyz = xyzfun.XYZPoint(
#                     x=geo_x, y=geo_y, z=z, w=weight_array[y,x]
#                 )
#                 yield(out_xyz)

#     z_array = weight_array = count_array = None
    
def parse(src_ds, dump_nodata=False, srcwin=None, mask=None, dst_srs=None, verbose=False, z_region=None, step=1):
    """parse the data from gdal dataset src_ds (first band only)
    optionally mask the output with `mask` or transform the coordinates to `dst_srs`

    yields the parsed xyz data
    """

    #if verbose: sys.stderr.write('waffles: parsing gdal file {}...'.format(src_ds.GetDescription()))
    ln = 0
    band = src_ds.GetRasterBand(1)
    ds_config = gather_infos(src_ds)
    src_srs = osr.SpatialReference()
    try:
        src_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    except: pass
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
    if srcwin is None: srcwin = (0, 0, ds_config['nx'], ds_config['ny'])
    nodata = ['{:g}'.format(-9999), 'nan', float('nan')]
    if band.GetNoDataValue() is not None: nodata.append('{:g}'.format(band.GetNoDataValue()))
    if dump_nodata: nodata = []
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

def yield_query(src_xyz, src_grd, out_form):
    """query a gdal-compatible grid file with xyz data.
    out_form dictates return values

    yields out_form results
    """
    
    try:
        ds = gdal.Open(src_grd)
    except: ds = None
    if ds is not None:
        ds_config = gather_infos(ds)
        ds_band = ds.GetRasterBand(1)
        ds_gt = ds_config['geoT']
        ds_nd = ds_config['ndv']
        tgrid = ds_band.ReadAsArray()
        dsband = ds = None
        
        for xyz in src_xyz:
            x = xyz[0]
            y = xyz[1]
            try: 
                z = xyz[2]
            except: z = ds_nd

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
                    yield(outs)
        ds = None

def query(src_xyz, src_grd, out_form):
    """query a gdal-compatible grid file with xyz data.
    out_form dictates return values

    returns array of values
    """
    
    xyzl = []
    for out_q in yield_query(src_xyz, src_grd, out_form):
        xyzl.append(np.array(out_q))
    return(np.array(xyzl))
    
### End
