### utils.py - CUDEM utilities and functions
##
## Copyright (c) 2010 - 2021 CIRES Coastal DEM Team
##
## utils.py is part of CUDEM
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
import hashlib
import glob
import time
import datetime
import requests
import subprocess

## ==============================================
## import gdal/numpy
## ==============================================
import numpy as np
import gdal
import ogr
import osr

__version__ = '0.5.0'
## ==============================================
##
## General Utility Functions, definitions, etc.
##
## ==============================================

## heaps of thanks to https://github.com/fitnr/stateplane
FIPS_TO_EPSG = {
    "0101": "26929", "0102": "26930", "0201": "26948", "0202": "26949",
    "0203": "26950", "0301": "26951", "0302": "26952", "0401": "26941",
    "0402": "26942", "0403": "26943", "0404": "26944", "0405": "26945",
    "0406": "26946", "0501": "26953", "0502": "26954", "0503": "26955",
    "0600": "26956", "0700": "26957", "0901": "26958", "0902": "26959",
    "0903": "26960", "1001": "26966", "1002": "26967", "1101": "26968",
    "1102": "26969", "1103": "26970", "1201": "26971", "1202": "26972",
    "1301": "26973", "1302": "26974", "1401": "26975", "1402": "26976",
    "1501": "26977", "1502": "26978", "1600": "03088", "1601": "2205",
    "1602": "26980", "1701": "26981", "1702": "26982", "1703": "32199",
    "1801": "26983", "1802": "26984", "1900": "26985", "2001": "26986",
    "2002": "26987", "2111": "26988", "2112": "26989", "2113": "26990",
    "2201": "26991", "2202": "26992", "2203": "26993", "2301": "26994",
    "2302": "26995", "2401": "26996", "2402": "26997", "2403": "26998",
    "2500": "32100", "2600": "32104", "2701": "32107", "2702": "32108",
    "2703": "32109", "2800": "32110", "2900": "32111", "3001": "32112",
    "3002": "32113", "3003": "32114", "3101": "32115", "3102": "32116",
    "3103": "32117", "3104": "32118", "3200": "32119", "3301": "32120",
    "3302": "32121", "3401": "32122", "3402": "32123", "3501": "32124",
    "3502": "32125", "3601": "32126", "3602": "32127", "3701": "32128",
    "3702": "32129", "3800": "32130", "3900": "32133", "4001": "32134",
    "4002": "32135", "4100": "32136", "4201": "32137", "4202": "32138",
    "4203": "32139", "4204": "32140", "4205": "32141", "4301": "32142",
    "4302": "32143", "4303": "32144", "4400": "32145", "4501": "32146",
    "4502": "32147", "4601": "32148", "4602": "32149", "4701": "32150",
    "4702": "32151", "4801": "32152", "4802": "32153", "4803": "32154",
    "4901": "32155", "4902": "32156", "4903": "32157", "4904": "32158",
    "5001": "26931", "5002": "26932", "5003": "26933", "5004": "26934",
    "5005": "26935", "5006": "26936", "5007": "26937", "5008": "26938",
    "5009": "26939", "5010": "26940", "5101": "26961", "5102": "26962",
    "5103": "26963", "5104": "26964", "5105": "26965", "5200": "32161"
}

def inc2str_inc(inc):
    """convert a WGS84 geographic increment to a str_inc (e.g. 0.0000925 ==> `13`)

    Args:
      inc (float): a gridding increment

    Returns:
      str: a str representation of float(inc)
    """
    
    import fractions
    return(str(fractions.Fraction(str(inc * 3600)).limit_denominator(10)).replace('/', ''))

def this_date():
    """get current data

    Returns:
      str: the current date
    """
    
    return(datetime.datetime.now().strftime('%Y%m%d%H%M%S'))

def this_year():
    """get the current year
    
    Returns
      str: the current year
    """
    
    return(datetime.datetime.now().strftime('%Y'))

def dl_hash(fn, sha1 = False):
    """calculate a hash of a file

    Args:
      fn (str): input filename path
    
    Returns:
      str: hexdigest
    """
    
    BUF_SIZE = 65536  # lets read stuff in 64kbchunks!
    if sha1: this_hash = hashlib.sha1()
    else: this_hash = hashlib.md5()

    with open(fn, 'rb') as f:
        while True:
            data = f.read(BUF_SIZE)
            if not data:
                break
            this_hash.update(data)
            
    return(this_hash.hexdigest())

def args2dict(args, dict_args={}):
    """convert list of arg strings to dict.
    
    Args:
      args (list): a list of ['key=val'] pairs
      dict_args (dict): a dict to append to

    Returns:
      dict: a dictionary of the key/values
    """
    
    for arg in args:
        p_arg = arg.split('=')
        dict_args[p_arg[0]] = False if p_arg[1].lower() == 'false' else True if p_arg[1].lower() == 'true' else None if p_arg[1].lower() == 'none' else p_arg[1]
    return(dict_args)

def remove_glob(*args):
    """glob `glob_str` and os.remove results, pass if error
    
    Args:
      *args (str): any number of pathname or dirname strings

    Returns:
      int: 0
    """
    
    for glob_str in args:
        try:
            globs = glob.glob(glob_str)
            for g in globs:
                if os.path.isdir(g):
                    remove_globs('{}/*'.format(g))
                    os.removedirs(g)
                else: os.remove(g)
        except: pass
    return(0)

def int_or(val, or_val = None):
    """return val if val is integer

    Args:
      val (?): input value to test
      or_val (?): value to return if val is not an int

    Returns:
      ?: val as int otherwise returns or_val
    """
    
    try:
        return(int(val))
    except: return(or_val)

def float_or(val, or_val = None):
    """return val if val is integer

    Args:
      val (?): input value to test
      or_val (?): value to return if val is not an int

    Returns:
      ?: val as int otherwise returns or_val
    """
    
    try:
        return(float(val))
    except: return(or_val)

def _geo2pixel(geo_x, geo_y, geoTransform):
    """convert a geographic x,y value to a pixel location of geoTransform

    Args:
      geo_x (float): geographic x coordinate
      geo_y (float): geographic y coordinate
      geoTransform (list): a geo-transform list describing a raster

    Returns:
      list: a list of the pixel values [pixel-x, pixel-y]
    """
    
    if geoTransform[2] + geoTransform[4] == 0:
        pixel_x = ((geo_x - geoTransform[0]) / geoTransform[1]) + .5
        pixel_y = ((geo_y - geoTransform[3]) / geoTransform[5]) + .5
    else: pixel_x, pixel_y = _apply_gt(geo_x, geo_y, _invert_gt(geoTransform))
    return(int(pixel_x), int(pixel_y))

def _geo2pixel_affine(geo_x, geo_y, geoTransform):
    """convert a geographic x,y value to a pixel location of geoTransform

    note: use _geo2pixel instead

    Args:
      geo_x (float): geographic x coordinate
      geo_y (float): geographic y coordinate
      geoTransform (list): a geo-transform list describing a raster

    Returns:
      list: a list of the pixel values [pixel-x, pixel-y]
    """
    
    import affine
    forward_transform = affine.Affine.from_gdal(*geoTransform)
    reverse_transform = ~forward_transform
    pixel_x, pixel_y = reverse_transform * (geo_x, geo_y)
    pixel_x, pixel_y = int(pixel_x + 0.5), int(pixel_y + 0.5)
    return(pixel_x, pixel_y)

def _pixel2geo(pixel_x, pixel_y, geoTransform):
    """convert a pixel location to geographic coordinates given geoTransform

    Args:
      pixel_x (int): the x pixel value
      pixel_y (int): the y pixel value
      geoTransform (list): a geo-transform list describing a raster
    
    Returns:
      list: [geographic-x, geographic-y]
    """
    
    geo_x, geo_y = _apply_gt(pixel_x, pixel_y, geoTransform)
    return(geo_x, geo_y)

def _apply_gt(in_x, in_y, geoTransform):
    """apply geotransform to in_x,in_y
    
    Args:
      in_x (int): the x pixel value
      in_y (int): the y pixel value
      geoTransform (list): a geo-transform list describing a raster

    Returns:
      list: [geographic-x, geographic-y]
    """
    
    out_x = geoTransform[0] + int(in_x + 0.5) * geoTransform[1] + int(in_y + 0.5) * geoTransform[2]
    out_y = geoTransform[3] + int(in_x + 0.5) * geoTransform[4] + int(in_y + 0.5) * geoTransform[5]

    return(out_x, out_y)

def _invert_gt(geoTransform):
    """invert the geotransform
    
    Args:
      geoTransform (list): a geo-transform list describing a raster

    Returns:
      list: a geo-transform list describing a raster
    """
    
    det = geoTransform[1] * geoTransform[5] - geoTransform[2] * geoTransform[4]
    if abs(det) < 0.000000000000001: return
    invDet = 1.0 / det
    outGeoTransform = [0, 0, 0, 0, 0, 0]
    outGeoTransform[1] = geoTransform[5] * invDet
    outGeoTransform[4] = -geoTransform[4] * invDet
    outGeoTransform[2] = -geoTransform[2] * invDet
    outGeoTransfrom[5] = geoTransform[1] * invDet
    outGeoTransform[0] = (geoTransform[2] * geoTransform[3] - geoTransform[0] * geoTransform[5]) * invDet
    outGeoTransform[3] = (-geoTransform[1] * geoTransform[3] + geoTransform[0] * geoTransform[4]) * invDet
    return(outGeoTransform)

def sr_wkt(epsg, esri = False):
    """convert an epsg code to wkt

    Returns:
      (str): wkt or None
    """
    
    try:
        sr = osr.SpatialReference()
        sr.ImportFromEPSG(int(epsg))
        if esri: sr.MorphToESRI()
        return(sr.ExportToWkt())
    except: return(None)

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
        if drv.GetMetadataItem(gdal.DCAP_RASTER): fexts = drv.GetMetadataItem(gdal.DMD_EXTENSIONS)
        if fexts is not None: return(fexts.split()[0])
        else: return(None)
    except:
        if src_drv_name.lower() == 'gtiff': fext = 'tif'
        elif src_drv_name == 'HFA': fext = 'img'
        elif src_drv_name == 'GMT': fext = 'grd'
        elif src_drv_name.lower() == 'netcdf': fext = 'nc'
        else: fext = 'gdal'
        return(fext)

def gdal_infos(src_fn, region = None, scan = False):
    '''scan gdal file src_fn and gather region info.

    returns region dict.'''
    if os.path.exists(src_fn):
        try:
            ds = gdal.Open(src_fn)
        except: ds = None
        if ds is not None:
            dsc = gdal_gather_infos(ds, region = region, scan = scan)
            ds = None
            return(dsc)
        else: return(None)
    else: return(None)

def gdal_gather_infos(src_ds, region = None, scan = False):
    '''gather information from `src_ds` GDAL dataset

    returns gdal_config dict.'''
    gt = src_ds.GetGeoTransform()
    if region is not None:
        srcwin = region.srcwin(gt, src_ds.RasterXSize, src_ds.RasterYSize)#gdal_srcwin(src_ds, region)
    else: srcwin = (0, 0, src_ds.RasterXSize, src_ds.RasterYSize)
    src_band = src_ds.GetRasterBand(1)
    dst_gt = (gt[0] + (srcwin[0] * gt[1]), gt[1], 0., gt[3] + (srcwin[1] * gt[5]), 0., gt[5])

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
        ds_config['zr'] = (np.amin(src_arr), np.amax(src_arr))
        #ds_config['zr'] = src_band.ComputeRasterMinMax()
        src_arr = None
    return(ds_config)

def gdal_set_epsg(src_fn, epsg = 4326):
    '''set the projection of gdal file src_fn to epsg

    returns status-code (0 == success)'''
    try:
        ds = gdal.Open(src_fn, gdal.GA_Update)
    except: ds = None
    if ds is not None:
        ds.SetProjection(sr_wkt(int(epsg)))
        ds = None
        return(0)
    else: return(None)

def gdal_set_infos(nx, ny, nb, geoT, proj, dt, ndv, fmt):
    '''set a datasource config dictionary

    returns gdal_config dict.'''
    return({'nx': nx, 'ny': ny, 'nb': nb, 'geoT': geoT, 'proj': proj, 'dt': dt, 'ndv': ndv, 'fmt': fmt})

def gdal_copy_infos(src_config):
    """copy src_config

    returns copied src_config dict.
    """
    
    dst_config = {}
    for dsc in src_config.keys():
        dst_config[dsc] = src_config[dsc]
    return(dst_config)

def gdal_split(src_gdal, split_value = 0):
    '''split raster file `src_gdal`into two files based on z value, 
    or if split_value is a filename, split raster by overlay, where upper is outside and lower is inside.

    returns [upper_grid-fn, lower_grid-fn]'''
    
    dst_upper = os.path.join(os.path.dirname(src_gdal), '{}_u.tif'.format(os.path.basename(src_gdal)[:-4]))
    dst_lower = os.path.join(os.path.dirname(src_gdal), '{}_l.tif'.format(os.path.basename(src_gdal)[:-4]))
    src_ds = gdal.Open(src_gdal)
    if src_ds is not None:
        src_config = gdal_gather_infos(src_ds)
        dst_config = gdal_cpy_infos(src_config)
        dst_config['fmt'] = 'GTiff'
        ds_arr = src_ds.GetRasterBand(1).ReadAsArray(0, 0, src_config['nx'], src_config['ny'])
        ua, la = np_split(ds_arr, split_value, src_config['ndv'])
        gdal_write(ua, dst_upper, dst_config)
        gdal_write(la, dst_lower, dst_config)
        ua = la = ds_arr = src_ds = None
        return([dst_upper, dst_lower])
    else: return(None)

def np_gaussian_blur(in_array, size):
    '''blur an array using fftconvolve from scipy.signal
    size is the blurring scale-factor.

    returns the blurred array'''
    
    from scipy.signal import fftconvolve
    from scipy.signal import convolve
    padded_array = np.pad(in_array, size, 'symmetric')
    x, y = np.mgrid[-size:size + 1, -size:size + 1]
    g = np.exp(-(x**2 / float(size) + y**2 / float(size)))
    g = (g / g.sum()).astype(in_array.dtype)
    in_array = None
    out_array = fftconvolve(padded_array, g, mode = 'valid')
    return(out_array)

def gdal_blur(src_gdal, dst_gdal, sf = 1):
    '''gaussian blur on src_gdal using a smooth-factor of `sf`
    runs np_gaussian_blur(ds.Array, sf)'''

    try:
        ds = gdal.Open(src_gdal)
    except: ds = None
    
    if ds is not None:
        ds_config = gdal_gather_infos(ds)
        ds_array = ds.GetRasterBand(1).ReadAsArray(0, 0, ds_config['nx'], ds_config['ny'])
        ds = None
        msk_array = np.array(ds_array)
        msk_array[msk_array != ds_config['ndv']] = 1
        msk_array[msk_array == ds_config['ndv']] = np.nan
        ds_array[ds_array == ds_config['ndv']] = 0
        smooth_array = np_gaussian_blur(ds_array, int(sf))
        smooth_array = smooth_array * msk_array
        mask_array = ds_array = None
        smooth_array[np.isnan(smooth_array)] = ds_config['ndv']
        return(gdal_write(smooth_array, dst_gdal, ds_config))
    else: return([], -1)

def gdal_filter_outliers(src_gdal, dst_gdal, threshhold = None, slp_threshhold = None, chunk_size = None, chunk_step = None, slp = False):
    '''scan a src_gdal file for outliers and remove them'''
    
    try:
        ds = gdal.Open(src_gdal)
    except: ds = None

    if ds is not None:
        tnd = 0
        
        ds_config = gdal_gather_infos(ds)
        ds_band = ds.GetRasterBand(1)
        ds_array = ds_band.ReadAsArray(0, 0, ds_config['nx'], ds_config['ny'])
        gt = ds_config['geoT']
        if threshhold is None:
            ds_std = np.std(ds_array)
        else: ds_std = threshhold
        if slp_threshhold is None:
            slp_std = ds_std
        else: slp_std = slp_threshhold

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
            n_step = int(n_chunk/4)
        else: n_step = chunk_step

        utils.echo_msg('scanning {} for spikes with {}@{} MAX {}/{}...'.format(src_gdal, n_chunk, n_step, ds_std, slp_std))
        for srcwin in gdal_yield_mw_srcwin(src_gdal, n_chunk = n_chunk, step = n_step):
            band_data = band.ReadAsArray(srcwin[0], srcwin[1], srcwin[2], srcwin[3])
            band_data[band_data == ds_config['ndv']] = np.nan
            this_geo_x_origin, this_geo_y_origin = utils._pixel2geo(srcwin[0], srcwin[1], gt)
            dst_gt = [this_geo_x_origin, float(gt[1]), 0.0, this_geo_y_origin, 0.0, float(gt[5])]

            dst_config = gdal_cpy_infos(ds_config)
            dst_config['nx'] = srcwin[2]
            dst_config['ny'] = srcwin[3]
            dst_config['geoT'] = dst_gt
            
            if not np.all(band_data == band_data[0,:]):
                while True:
                    nd = 0                    
                    srcwin_std = np.nanstd(band_data)
                    slp_data = np.gradient(band_data, axis=0)
                    slp_srcwin_std = np.nanstd(slp_data)
                    if srcwin_std < ds_std and slp_srcwin_std < slp_std: break
                    
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
                                    ds_array[j+srcwin[1]][i+srcwin[0]] = ds_config['ndv']
                                    band.WriteArray(ds_array)
                                    band_data[j][i] = np.nan
                                    nd += 1
                    tnd += nd
                    if nd == 0: break
                band_data = slp_data = None
                
        utils.echo_msg('filtering {} spikes...'.format(tnd))
        if tnd > 0:
            driver = gdal.GetDriverByName('MEM')
            tmp_ds = driver.Create('tmp', ds_config['nx'], ds_config['ny'], 1, ds_config['dt'])
            tmp_ds.SetGeoTransform(ds_config['geoT'])
            tmp_ds.SetProjection(ds_config['proj'])
            ds_band = tmp_ds.GetRasterBand(1)
            ds_band.SetNoDataValue(ds_config['ndv'])
            ds_band.WriteArray(ds_array)
            result = gdal.FillNodata(targetBand = ds_band, maskBand = None, maxSearchDist = 100, smoothingIterations = 4, callback = _gdal_progress)
        
            ds_array = ds_band.ReadAsArray()
            tmp_ds = None

        out, status = gdal_write(ds_array, dst_gdal, ds_config)
        mem_ds = None
        return(out, status)
    else: return(None)

## ==============================================
## Write an array to a gdal file
## ==============================================
def gdal_write (src_arr, dst_gdal, ds_config, dst_fmt = 'GTiff'):
    """write src_arr to gdal file dst_gdal using src_config

    returns [output-gdal, status-code]
    """
    
    driver = gdal.GetDriverByName(dst_fmt)
    if os.path.exists(dst_gdal):
        try:
            driver.Delete(dst_gdal)
        except Exception as e:
            utils.echo_error_msg(e)
            utils.remove_glob(dst_gdal)
    ds = driver.Create(dst_gdal, ds_config['nx'], ds_config['ny'], 1, ds_config['dt'])
    if ds is not None:
        ds.SetGeoTransform(ds_config['geoT'])
        ds.SetProjection(ds_config['proj'])
        ds.GetRasterBand(1).SetNoDataValue(ds_config['ndv'])
        ds.GetRasterBand(1).WriteArray(src_arr)
        ds = None
        return(dst_gdal, 0)
    else: return(None, -1)

def mb_inf(src_xyz, src_fmt = 168):
    """generate an info (.inf) file from a src_xyz file using MBSystem.

    Args:
      src_xyz (port): an open xyz port or generator
      src_fmt (int): the datalists format of the source data

    Returns:
      dict: xyz infos dictionary mb_inf_parse(inf_file)
    """

    run_cmd('mbdatalist -O -F{} -I{}'.format(src_fmt, src_xyz.name), verbose = False)
    return(mb_inf_parse('{}.inf'.format(src_xyz.name)))

def gmt_inc2inc(inc_str):
    """convert a GMT-style `inc_str` (6s) to geographic units

    c/s - arc-seconds
    m - arc-minutes

    Args:
      inc_str (str): GMT style increment string

    Returns:
      float: increment value.
    """
    
    if inc_str is None or inc_str.lower() == 'none': return(None)
    units = inc_str[-1]
    if units == 'c': inc = float(inc_str[:-1]) / 3600.
    elif units == 's': inc = float(inc_str[:-1]) / 3600.
    elif units == 'm': inc = float(inc_str[:-1]) / 360.
    else:
        try:
            inc = float(inc_str)
        except ValueError as e:
            echo_error_msg('could not parse increment {}, {}'.format(inc_str, e))
            return(None)
    return(inc)

def gmt_grdfilter(src_grd, dst_grd, dist = '3s', node = 'pixel', verbose = False):
    """filter `src_grd` using GMT grdfilter

    Args:
      src_grd (str): pathname to a source grid file
      dst_grd (str): pathname to a destination grid file
      dist (str): a GMT string increment to use in filter
      node (str): `pixel` or `grid`; the grid-node
      verbose (bool): increase verbosity

    Returns:
      list: [cmd-output, cmd-return-code]
    """
    
    #ft_cmd1 = ('gmt grdfilter -V {} -G{} -R{} -Fc{} -D1{}'.format(src_grd, dst_grd, src_grd, dist, ' -r' if node == 'pixel' else ''))
    ft_cmd1 = ('gmt grdfilter -V {} -G{} -Fc{} -D1{}'.format(src_grd, dst_grd, dist, ' -r' if node == 'pixel' else ''))
    return(utils.run_cmd(ft_cmd1, verbose = verbose))

## ==============================================
##
## system cmd verification and configs.
##
## run_cmd - run a system command as a subprocess
## and return the return code and output
##
## yield_cmd - run a system command as a subprocess
## and yield the output
##
## ==============================================
cmd_exists = lambda x: any(os.access(os.path.join(path, x), os.X_OK) for path in os.environ['PATH'].split(os.pathsep))

def run_cmd(cmd, data_fun=None, verbose=False):
    """Run a system command while optionally passing data.

    `data_fun` should be a function to write to a file-port:
    >> data_fun = lambda p: datalist_dump(wg, dst_port = p, ...)

    Args:
      cmd (str): a system command to run
      data_fun (lambda): a lambda function taking an output port as arg.
      verbose (bool): increase verbosity

    Returns:
      list: [command-output, command-return-code]
    """
    
    if verbose: _prog = CliProgress('running cmd: {}...'.format(cmd.rstrip()[:72]))
    if data_fun is not None:
        pipe_stdin = subprocess.PIPE
    else: pipe_stdin = None
    if verbose:
        p = subprocess.Popen(cmd, shell=True, stdin=pipe_stdin, stdout=subprocess.PIPE, close_fds=True)
    else: p = subprocess.Popen(cmd, shell=True, stdin=pipe_stdin, stdout=subprocess.PIPE, stderr=subprocess.PIPE, close_fds=True)

    if data_fun is not None:
        if verbose: echo_msg('piping data to cmd subprocess...')
        data_fun(p.stdin)
        p.stdin.close()
    
    while p.poll() is None:
        if verbose:
            time.sleep(2)
            _prog.update()

    out = p.stdout.read()
    if not verbose: p.stderr.close()
    p.stdout.close()
    if verbose: _prog.end(p.returncode, 'ran cmd: {}... and returned {}.'.format(cmd.rstrip()[:72], p.returncode))
    return(out, p.returncode)

def yield_cmd(cmd, data_fun=None, verbose=False):
    """Run a system command while optionally passing data.

    `data_fun` should be a function to write to a file-port:
    >> data_fun = lambda p: datalist_dump(wg, dst_port = p, ...)

    Args:
      cmd (str): a system command to run
      data_fun (lambda): a lambda function taking an output port as arg.
      verbose (bool): increase verbosity

    Yields:
      str: each line of output from the cmd
    """
    
    if verbose: echo_msg('running cmd: {}...'.format(cmd.rstrip()))    
    if data_fun is not None:
        pipe_stdin = subprocess.PIPE
    else: pipe_stdin = None
    p = subprocess.Popen(cmd, shell=True, stdin=pipe_stdin, stdout=subprocess.PIPE, close_fds=True)
    
    if data_fun is not None:
        if verbose: echo_msg('piping data to cmd subprocess...')
        data_fun(p.stdin)
        p.stdin.close()

    while p.poll() is None:
        line = p.stdout.readline().decode('utf-8')
        if not line: break
        else: yield(line)
    line = p.stdout.read().decode('utf-8')
    p.stdout.close()
    if verbose: echo_msg('ran cmd: {} and returned {}.'.format(cmd.rstrip(), p.returncode))

def cmd_check(cmd_str, cmd_vers_str):
    """check system for availability of 'cmd_str' 

    Args:
      cmd_str (str): a command string to run
      cmd_vers_str (str): a command that returns a version

    Returns:
      str: the commands version or None
    """
    
    if cmd_exists(cmd_str): 
        cmd_vers, status = run_cmd('{}'.format(cmd_vers_str))
        return(cmd_vers.rstrip())
    else: return("0")

def config_check(chk_vdatum=False, verbose=False):
    """check for needed waffles external software.

    waffles external software: gdal, gmt, mbsystem
    also checks python version and host OS and 
    records waffles version

    Args:
      chk_vdatum (bool): check for vdatum
      verbose (bool): increase verbosity

    Returns:
      dict: a dictionary of gathered results.
    """
    
    _waff_co = {}
    py_vers = str(sys.version_info[0]),
    host_os = sys.platform
    _waff_co['platform'] = host_os
    _waff_co['python'] = py_vers[0]
    ae = '.exe' if host_os == 'win32' else ''

    #if chk_vdatum: _waff_co['VDATUM'] = vdatum(verbose=verbose).vdatum_path
    _waff_co['GDAL'] = cmd_check('gdal_grid{}'.format(ae), 'gdal-config --version').decode()
    _waff_co['GMT'] = cmd_check('gmt{}'.format(ae), 'gmt --version').decode()
    _waff_co['MBGRID'] = cmd_check('mbgrid{}'.format(ae), 'mbgrid -version 2>&1 | grep Version').decode()
    _waff_co['LASTOOLS'] = cmd_check('las2txt{}'.format(ae), 'las2txt -version 2>&1 | awk \'{print $5}\'').decode()
    _waff_co['GEOMODS'] = str(__version__)
    return(_waff_co)
    
## ==============================================
##
## verbosity functions
##
## ==============================================
def echo_warning_msg2(msg, prefix = 'waffles'):
    """echo warning msg to stderr using `prefix`

    >> echo_warning_msg2('message', 'test')
    test: warning, message
    
    Args:
      msg (str): a message
      prefix (str): a prefix for the message
    """
    
    sys.stderr.flush()
    sys.stderr.write('\x1b[2K\r')
    sys.stderr.write('{}: \033[31m\033[1mwarining\033[m, {}\n'.format(prefix, msg))

def echo_error_msg2(msg, prefix = 'waffles'):
    """echo error msg to stderr using `prefix`

    >> echo_error_msg2('message', 'test')
    test: error, message

    Args:
      msg (str): a message
      prefix (str): a prefix for the message
    """
    
    sys.stderr.flush()
    sys.stderr.write('\x1b[2K\r')
    sys.stderr.write('{}: \033[31m\033[1merror\033[m, {}\n'.format(prefix, msg))

def echo_msg2(msg, prefix = 'waffles', nl = True):
    """echo `msg` to stderr using `prefix`

    >> echo_msg2('message', 'test')
    test: message
    
    Args:
      msg (str): a message
      prefix (str): a prefix for the message
      nl (bool): append a newline to the message
    """
    
    sys.stderr.flush()
    sys.stderr.write('\x1b[2K\r')
    sys.stderr.write('{}: {}{}'.format(prefix, msg, '\n' if nl else ''))
    sys.stderr.flush()
    
## ==============================================
## echo message `m` to sys.stderr using
## auto-generated prefix
## lambda runs: echo_msg2(m, prefix = os.path.basename(sys.argv[0]))
## ==============================================
echo_msg = lambda m: echo_msg2(m, prefix = os.path.basename(sys.argv[0]))
echo_msg_inline = lambda m: echo_msg2(m, prefix = os.path.basename(sys.argv[0]), nl = False)

## ==============================================
## echo error message `m` to sys.stderr using
## auto-generated prefix
## ==============================================
echo_error_msg = lambda m: echo_error_msg2(m, prefix = os.path.basename(sys.argv[0]))
echo_warning_msg = lambda m: echo_warning_msg2(m, prefix = os.path.basename(sys.argv[0]))

class CliProgress:
    """geomods minimal progress indicator
    """

    def __init__(self, message = None):
        self.tw = 7
        self.count = 0
        self.pc = self.count % self.tw
        self.opm = message
        self.add_one = lambda x: x + 1
        self.sub_one = lambda x: x - 1
        self.spin_way = self.add_one
        self.spinner = ['*     ', '**    ', '***   ', ' ***  ', '  *** ', '   ***', '    **', '     *']
        
        self.perc = lambda p: ((p[0]/p[1]) * 100.)
        
        if self.opm is not None:
            self._clear_stderr()
            sys.stderr.write('\r {}  {:40}\n'.format(" " * (self.tw - 1), self.opm))
        
    def _switch_way(self):
        self.spin_way = self.sub_one if self.spin_way == self.add_one else self.add_one

    def _clear_stderr(self, slen = 79):
        sys.stderr.write('\x1b[2K\r')
        sys.stderr.flush()

    def update_perc(self, p, msg = None):
        if len(p) == 2:
            self._clear_stderr()
            sys.stderr.write('\r[\033[36m{:^5.2f}%\033[m] {:40}\r'.format(self.perc(p), msg if msg is not None else self.opm))
        else: self.update()
        
    def update(self, msg = None):

        self.pc = (self.count % self.tw)
        self.sc = (self.count % (self.tw + 1))
            
        self._clear_stderr()
        sys.stderr.write('\r[\033[36m{:6}\033[m] {:40}\r'.format(self.spinner[self.sc], msg if msg is not None else self.opm))
        
        if self.count == self.tw: self.spin_way = self.sub_one
        if self.count == 0: self.spin_way = self.add_one
        self.count = self.spin_way(self.count)

    def end(self, status, end_msg = None):
        self._clear_stderr()
        if end_msg is None: end_msg = self.opm
        if status != 0:
            sys.stderr.write('\r[\033[31m\033[1m{:^6}\033[m] {:40}\n'.format('fail', end_msg))
        else: sys.stderr.write('\r[\033[32m\033[1m{:^6}\033[m] {:40}\n'.format('ok', end_msg))
    
### End
