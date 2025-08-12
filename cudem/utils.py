### utils.py - CUDEM utilities and functions
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
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
## General utility functions and classes used in other modules...
## includes, generic gdal functions, fetching classes, progress indicators,
## system command interfaces, etc...
##
### Code:

import os
import sys
import hashlib
import glob
import math
import time
import datetime
import shutil
import subprocess
import zipfile
import tarfile
import gzip
import bz2
import re
import curses
import io
import json
import traceback
import getpass

from tqdm import tqdm
import threading
try:
    import Queue as queue
except: import queue as queue

import numpy as np
import cudem

###############################################################################
##
## General Utility Functions, definitions, etc.
##
###############################################################################

## Cahce directory, for use in dlim/waffles/fetches
this_dir, this_filename = os.path.split(__file__)
cudem_cache = lambda: os.path.abspath('./.cudem_cache')
cudem_data = os.path.join(this_dir, 'data')
def set_cache(cache_dir: str):
    if not os.path.exists(cache_dir):
        os.makedirs(cache_dir)

def remove_cache(cache_dir: str):
    if os.path.exists(cache_dir):
        remove_glob(cache_dir)

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

## General file-name helper functions
def append_fn(fn, src_region, inc, version=None, year=None,
              res=None, high_res=False):
    """append the src_region, inc and version to a string filename"""
    
    return(
        '{}{}_{}_{}v{}'.format(
            fn,
            inc2str(inc) if res is None else res,
            src_region.format('fn' if not high_res else 'fn_full'),
            this_year() if year is None else year,
            1 if version is None else version
        )
    )


def fn_basename(fn, ext):
    """return the basename of fn based on ext"""
    
    if '.' in ext:
        return(fn[:-len(ext)])
    else:
        return(fn[:-(len(ext)+1)])

    
def fn_basename2(fn):
    """return the basename of fn"""
    
    t = fn.split('.')
    if len(t) > 1:
        return('.'.join(fn.split('.')[:-1]))
    else:
        return(fn)

    
def fn_ext(fn):
    """return the extension of fn"""

    ext = None
    t = fn.split('.')
    if len(t) > 1:
        ext = t[-1]

    return(ext)


def make_temp_fn(fn, temp_dir = cudem_cache()):
    """make a temporary unique filename"""
    
    if temp_dir is None:
        temp_dir = cudem_cache()
        
    fn_bn = fn_basename2(os.path.basename(fn))
    fn_et = fn_ext(fn)
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)
                    
    return(os.path.join(
        temp_dir,
        '{}_{}{}'.format(
            fn_bn, datetime.datetime.now().strftime('%Y%m%H%M%S%f'),
            '.{}'.format(fn_et) if fn_et is not None else '')
    ))


def fn_url_p(fn):
    """check if fn is a url"""
    
    url_sw = ['http://', 'https://', 'ftp://', 'ftps://', '/vsicurl']
    for u in url_sw:
        try:
            if fn.startswith(u):
                return(True)
        except:
            return(False)
        
    return(False)


def inc2str(inc):
    """convert a WGS84 geographic increment to a str_inc 
    (e.g. 0.0000925 ==> `13`)

    Args:
      inc (float): a gridding increment

    Returns:
      str: a str representation of float(inc)
    """
    
    import fractions
    return(
        str(
            fractions.Fraction(str(inc * 3600)).limit_denominator(10)
        ).replace('/', '')
    )


def str2inc(inc_str):
    """convert a GMT-style `inc_str` (e.g. 6s) to geographic units

    c/s - arc-seconds
    m - arc-minutes

    Args:
      inc_str (str): GMT style increment string

    Returns:
      float: increment value.
    """

    try:
        inc_str = str(inc_str)
    except:
        return(None)
    
    if inc_str is None or inc_str.lower() == 'none' or len(inc_str) == 0:
        return(None)
    
    units = inc_str[-1]
    if units == 'c':
        inc = float(inc_str[:-1]) / 3600.
    elif units == 's':
        inc = float(inc_str[:-1]) / 3600.
    elif units == 'm':
        inc = float(inc_str[:-1]) / 360.
    else:
        try:
            inc = float_or(inc_str)            
        except ValueError as e:
            echo_error_msg(
                f'could not parse increment {inc_str}, {e}'
            )
            return(None)
        
    return(float(inc))


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


def dl_hash(fn, sha1=False):
    """calculate a hash of a file

    Args:
      fn (str): input filename path
    
    Returns:
      str: hexdigest
    """
    
    BUF_SIZE = 65536
    if sha1:
        this_hash = hashlib.sha1()
    else:
        this_hash = hashlib.md5()

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
        if len(p_arg) > 1:
            dict_args[p_arg[0]] = False if p_arg[1].lower() == 'false' \
                else True if p_arg[1].lower() == 'true' \
                     else None if p_arg[1].lower() == 'none' \
                          else '='.join(p_arg[1:]) if len(p_arg) > 2 \
                               else p_arg[1]
        
    return(dict_args)


def dict2args(in_dict):
    """convert a dictionary to an args string"""
    
    out_args = ''
    for i, key in enumerate(in_dict.keys()):
        out_args += '{}={}{}'.format(
            key, in_dict[key], ':' if i+1 < len(in_dict.keys()) else ''
        )
    return(out_args)


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
                    remove_glob('{}/*'.format(g))
                    remove_glob('{}/.*'.format(g))
                    os.removedirs(g)
                else:
                    os.remove(g)
                    
        except Exception as e:
            echo_error_msg(e)
            return(-1)
       
    return(0)


def slugify(value, allow_unicode=False):
    """Taken from https://github.com/django/django/blob/master/django/utils/text.py
    Convert to ASCII if 'allow_unicode' is False. Convert spaces or repeated
    dashes to single dashes. Remove characters that aren't alphanumerics,
    underscores, or hyphens. Convert to lowercase. Also strip leading and
    trailing whitespace, dashes, and underscores.
    """

    import unicodedata
    import re
    
    value = str(value)
    if allow_unicode:
        value = unicodedata.normalize('NFKC', value)
    else:
        value = unicodedata.normalize('NFKD', value).encode('ascii', 'ignore').decode('ascii')
        
    value = re.sub(r'[^\w\s-]', '', value.lower())
    return(re.sub(r'[-\s]+', '-', value).strip('-_'))


def flatten_recursive(nested_list):
    flat_list = []
    for item in nested_list:
        if isinstance(item, list):
            flat_list.extend(flatten_recursive(item))
        else:
            flat_list.append(item)
            
    return(flat_list)


def dict_path2abspath(d = {}, except_keys = []):
    for key in d.keys():
        if key in except_keys:
            continue
        
        if isinstance(d[key], dict):
            d[key] = dict_path2abspath(d[key])
        elif isinstance(d[key], str):
            # if len(factory.fmod2dict(d[key])) > 1:
            #     dd = factory.fmod2dict(d[key])
            #     dd = dict_path2abspath(dd)
            #     d[key] = factory.dict2fmod(dd)            
            if os.path.exists(d[key]):
                d[key] = os.path.abspath(d[key])
                
        # elif isinstance(d[key], list):
        #     d_list = []
        #     for dd in d[key]:
        #         if isinstance(dd, str):
        #             if len(factory.fmod2dict(dd)) > 1:
        #                 dd = factory.fmod2dict(dd)
        #                 dd = dict_path2abspath(dd)
        #                 dd = factory.dict2fmod(dd)            
        #             elif os.path.exists(d[key]):
        #                 dd = os.path.abspath(dd)
        #         d_list.append(dd)
        #     d[key] = d_list
            
    return(d)
                

def int_or(val, or_val=None):
    """return val if val is integer

    Args:
      val (?): input value to test
      or_val (?): value to return if val is not an int

    Returns:
      ?: val as int otherwise returns or_val
    """
    
    try:
        return(int(float_or(val)))
    except: return(or_val)

    
def float_or(val, or_val=None):
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

    
def str_or(instr, or_val=None, replace_quote=True):
    """return instr if instr is a string, else or_val"""

    if instr is None:
        return(or_val)
    
    try:
        if replace_quote:
            return(str(instr).replace('"', ''))
        else:
            return(str(instr))
    except:
        return(or_val)

    
def convert_size(size_bytes):
   if size_bytes == 0:
       return('0B')
   
   size_name = ('B', 'KB', 'MB', 'GB', 'TB', 'PB', 'EB', 'ZB', 'YB')
   i = int(math.floor(math.log(size_bytes, 1024)))
   p = math.pow(1024, i)
   s = round(size_bytes / p, 2)
   return(f'{s} {size_name[i]}')


def euc_dst(pnt0, pnt1):
    """return the distance between pnt0 and pnt1,
    using the euclidean formula.

    `pnts` are geographic and result is in meters.

    Args:
      pnt0 (list): an xyz data list
      pnt1 (list): an xyz data list

    Returns:
      float: the distance beteween pnt0 and pnt1
    """
    
    rad_m = 637100
    distance = math.sqrt(sum([(a-b) ** 2 for a, b in zip(pnt0, pnt1)]))
    return(rad_m * distance)


def hav_dst(pnt0, pnt1):
    """return the distance between pnt0 and pnt1,
    using the haversine formula.

    `pnts` are geographic and result is in meters.

    Args:
      pnt0 (list): an xyz data list
      pnt1 (list): an xyz data list

    Returns:
      float: the distance beteween pnt0 and pnt1
    """
    
    x0 = float(pnt0[0])
    y0 = float(pnt0[1])
    x1 = float(pnt1[0])
    y1 = float(pnt1[1])
    rad_m = 637100
    dx = math.radians(x1 - x0)
    dy = math.radians(y1 - y0)
    a = math.sin(dx/2)*math.sin(dx/2)+math.cos(math.radians(x0))*math.cos(math.radians(x1))*math.sin(dy/2)*math.sin(dy/2)
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a)) 
    return(rad_m*c)


def wgs_inc2meter(src_inc):
    """return a wgs increment as meters"""
    
    gds_equator = 111321.543
    degree_to_radian = lambda d: math.pi * (d / 180)
    return(math.cos(degree_to_radian(src_inc)) * (gds_equator * src_inc))


def lll(src_lat):
    """return the lon/lat length in meters"""
    
    gds_equator = 111321.543
    gds_pi = 3.14159265358979323846
    degree_to_radian = lambda d: gds_pi * (d / 180)
    lonl_ = math.cos(degree_to_radian(src_lat)) * gds_equator
    latl_ = gds_equator
    return(lonl_, latl_)


def touch(fname, times = None):
    """touch a file to make sure it exists"""
    
    with open(fname, 'a'):
        os.utime(fname, times)
        
    return(fname)


def get_username():
    username = ''

    # For Python 2/3 compatibility:
    try:
        do_input = raw_input  # noqa
    except NameError:
        do_input = input

    while not username:
        username = do_input('username: ')
        
    return username


def get_password():
    password = ''
    while not password:
        password = getpass('password: ')
        
    return password


def get_outliers(
        in_array: any,
        percentile: float = 75,
        k: float = 1.5,
        verbose: bool = False
):
    """get the outliers from in_array based on the percentile

    https://en.wikipedia.org/wiki/Interquartile_range
    """

    if verbose:
        utils.echo_msg(
            f'input percentile: {percentile}'
        )

    if percentile < 0:
        percentile = 0

    if percentile > 100:
        percentile = 100

    max_percentile = percentile
    min_percentile = percentile-50

    if min_percentile < 0:
        min_percentile = 0

    if verbose:
        utils.echo_msg(
            f'percentiles: {min_percentile}>>{max_percentile}'
        )

    if np.all(np.isnan(in_array)):
        upper_limit = np.nan
        lower_limit = np.nan
    else:
        perc_max = np.nanpercentile(in_array, max_percentile)
        perc_min = np.nanpercentile(in_array, min_percentile)
        iqr_p = (perc_max - perc_min) * k
        upper_limit = perc_max + iqr_p
        lower_limit = perc_min - iqr_p

    return(upper_limit, lower_limit)


## ==============================================
##
## Geotransform functions
##
## ==============================================
def _geo2pixel(geo_x, geo_y, geo_transform, node='grid'):
    """convert a geographic x,y value to a pixel location of 
    geoTransform

    Args:
      geo_x (float): geographic x coordinate
      geo_y (float): geographic y coordinate
      geo_transform (list): a geo-transform list describing a raster
      node (str): the registration of the geo_transform

    Returns:
      list: a list of the pixel values [pixel-x, pixel-y]
    """

    if geo_transform[2] + geo_transform[4] == 0:
        ## rounding for edges
        #pixel_x = round((geo_x - geo_transform[0]) / geo_transform[1], 4)
        #pixel_y = round((geo_y - geo_transform[3]) / geo_transform[5], 4)

        ## numpy
        #pixel_x = np.floor((geo_x - geo_transform[0]) / geo_transform[1]).astype(int)
        #pixel_y = np.floor((geo_y - geo_transform[3]) / geo_transform[5]).astype(int)

        ## query
        pixel_x = (geo_x - geo_transform[0]) / geo_transform[1]
        pixel_y = (geo_y - geo_transform[3]) / geo_transform[5]

        ## grid-node geo_transform
        if node == 'grid':
            pixel_x += .5
            pixel_y += .5
    else:
        pixel_x, pixel_y = _apply_gt(geo_x, geo_y, _invert_gt(geo_transform))
        
    return(int(pixel_x), int(pixel_y))


def __geo2pixel(geo_x, geo_y, geo_transform, node='pixel'):
    """convert a geographic x,y value to a pixel location of geoTransform

    note: use _geo2pixel instead

    Args:
      geo_x (float): geographic x coordinate
      geo_y (float): geographic y coordinate
      geo_transform (list): a geo-transform list describing a raster

    Returns:
      list: a list of the pixel values [pixel-x, pixel-y]
    """
    
    import affine
    forward_transform = affine.Affine.from_gdal(*geo_transform)
    reverse_transform = ~forward_transform
    pixel_x, pixel_y = reverse_transform * (geo_x, geo_y)
    pixel_x, pixel_y = int(pixel_x+0.5), int(pixel_y+0.5)
    
    return(pixel_x, pixel_y)


def _pixel2geo(pixel_x, pixel_y, geo_transform, node='pixel',
               x_precision=None, y_precision=None):
    """convert a pixel location to geographic coordinates given geoTransform

    Args:
      pixel_x (int): the x pixel value
      pixel_y (int): the y pixel value
      geo_transform (list): a geo-transform list describing a raster
    
    Returns:
      list: [geographic-x, geographic-y]
    """
    
    geo_x, geo_y = _apply_gt(pixel_x, pixel_y, geo_transform, node)
    x_precision = int_or(x_precision)
    y_precision = int_or(y_precision)
    if x_precision is None:
        if y_precision is None:
            return(geo_x, geo_y)
        else:
            return(geo_x, round(geo_y, y_precision))
    else:
        if y_precision is None:
            return(round(geo_x, x_precision), geo_y)
        else:
            return(round(geo_x, x_precision), round(geo_y, y_precision))

        
def _apply_gt(in_x, in_y, geo_transform, node='pixel'):
    """apply geotransform to in_x,in_y
    
    Args:
      in_x (int): the x pixel value
      in_y (int): the y pixel value
      geoTransform (list): a geo-transform list describing a raster

    Returns:
      list: [geographic-x, geographic-y]
    """

    if node == 'pixel':
        out_x = geo_transform[0] + (in_x+0.5) * geo_transform[1] + (in_y+0.5) * geo_transform[2]
        out_y = geo_transform[3] + (in_x+0.5) * geo_transform[4] + (in_y+0.5) * geo_transform[5]
    else:
        out_x = geo_transform[0] + in_x * geo_transform[1] + in_y * geo_transform[2]
        out_y = geo_transform[3] + in_x * geo_transform[4] + in_y * geo_transform[5]

    return(out_x, out_y)


def _invert_gt(geo_transform):
    """invert the geo_transform
    
    Args:
      geo_transform (list): a geo-transform list describing a raster

    Returns:
      list: a geo-transform list describing a raster
    """
    
    det = (geo_transform[1]*geo_transform[5]) - (geo_transform[2]*geo_transform[4])
    if abs(det) < 0.000000000000001: return
    
    inv_det = 1.0 / det
    out_geo_transform = [0, 0, 0, 0, 0, 0]
    out_geo_transform[1] = geo_transform[5] * inv_det
    out_geo_transform[4] = -geo_transform[4] * inv_det
    out_geo_transform[2] = -geo_transform[2] * inv_det
    out_geo_transfrom[5] = geo_transform[1] * inv_Det
    out_geo_transform[0] = (geo_transform[2] * geo_transform[3] - geo_transform[0] * geo_transform[5]) * inv_det
    out_geo_Transform[3] = (-geo_transform[1] * geo_transform[3] + geo_transform[0] * geo_transform[4]) * inv_det
    return(outGeoTransform)


def x360(x):
    if x == 0:
        return(-180)
    elif x == 360:
        return(180)
    else:
        return(((x + 180) % 360) - 180)


###############################################################################
##
## Archives (zip/gzip/etc.)
##
###############################################################################
def unbz2(bz2_file, outdir='./', overwrite=False):

    newfilepath = '.'.join(bz2_file.split('.')[:-1])

    if not os.path.exists(newfilepath):
        echo_msg('Uncompressing {} to {}...'.format(bz2_file, newfilepath))
        with open(newfilepath, 'wb') as new_file, bz2.BZ2File(bz2_file, 'rb') as file:
            for data in iter(lambda : file.read(100 * 1024), b''):
                new_file.write(data)
    return(newfilepath)


def zip_list(zip_file):
    try:
        zip_ref = zipfile.ZipFile(zip_file)
        zip_files = zip_ref.namelist()
        
        return(zip_files)
    
    except Exception as e:
        echo_error_msg(e)

        return(None)

    
def unzip(zip_file, outdir='./', overwrite=False, verbose=True):
    """unzip (extract) `zip_file`

    Args:
      zip_file (str): a zip file pathname string

    Returns:
      list: a list of extracted file names.
    """

    try:
        zip_ref = zipfile.ZipFile(zip_file)
        zip_files = zip_ref.namelist()
        if not overwrite:
            for fn in zip_files:
                if not os.path.exists(os.path.join(outdir, fn)):
                    if verbose:
                        echo_msg(
                            'Extracting {}'.format(os.path.join(outdir, fn))
                        )
                        
                    zip_ref.extract(fn, outdir)
        else:
            zip_ref.extractall(outdir)
            
        zip_ref.close()
        if outdir != './':
            for i, zf in enumerate(zip_files):
                zip_files[i] = os.path.join(outdir, zf)
                
        return(zip_files)
    
    except Exception as e:
        echo_error_msg(e)
        
        return(None)

    
def gunzip(gz_file, outdir='./'):
    """gunzip `gz_file`

    Args:
      gz_file (str): a gzip file pathname string.

    Returns:
      str: the extracted file name.
    """
    
    if os.path.exists(gz_file):
        guz_file = os.path.join(outdir, os.path.basename(fn_basename2(gz_file)))
        with gzip.open(gz_file, 'rb') as in_gz, \
             open(guz_file, 'wb') as f:
            while True:
                block = in_gz.read(65536)
                if not block:
                    break
                else: f.write(block)
    else:
        echo_error_msg('{} does not exist'.format(gz_file))
        guz_file = None
        
    return(guz_file)


def p_untar(tar_file, exts=None, outdir='./'):
    src_procs = []
    with tarfile.open(tar_file, 'r') as tar:
        tar_fns = tar.getnames()

        for ext in exts:
            for tfn in tar_fns:
                #if ext == tfn.split('.')[-1]:
                if ext == tfn[-len(ext):]:
                    ext_tfn = os.path.join(outdir, os.path.basename(tfn))
                    src_procs.append(ext_tfn)
                    if not os.path.exists(ext_tfn):
                        echo_msg('Extracting {}'.format(ext_tfn))
                        t = tar.extractfile(tfn)
                        with open(ext_tfn, 'wb') as f:
                            f.write(t.read())
    return(src_procs)


def gdb_unzip(src_zip, outdir='./', verbose=True):
    src_gdb = None
    with zipfile.ZipFile(src_zip) as z:
        zfs = z.namelist()
        ext_mask = ['.gdb' in x for x in zfs]

        for i, zf in enumerate(zfs):
            if ext_mask[i]:
                ext_zf = os.path.join(outdir, zf)
                if not os.path.exists(ext_zf) and not ext_zf.endswith('/'):
                    _dirname = os.path.dirname(ext_zf)
                    if not os.path.exists(_dirname):
                        os.makedirs(_dirname)

                    with open(ext_zf, 'wb') as f:
                        f.write(z.read(zf))

                    if verbose:
                        echo_msg('Extracted {}'.format(ext_zf))
                        
                elif ext_zf.endswith('.gdb/'):
                    src_gdb = ext_zf

    return(src_gdb)


def p_unzip(src_file, exts=None, outdir='./', verbose=True):
    """unzip/gunzip src_file based on `exts`
    
    Args:
      src_file (str): a zip/gzip filename string
      exts (list): a list of extensions to extract

    Returns:
      list: a list of the extracted files
    """
    
    src_procs = []
    if src_file.split('.')[-1].lower() == 'zip':
        try:
            with zipfile.ZipFile(src_file) as z:
                zfs = z.namelist()
                for ext in exts:
                    for zf in zfs:
                        if ext == zf.split('.')[-1]:
                            ext_zf = os.path.join(outdir, zf)
                            src_procs.append(ext_zf)
                            if not os.path.exists(ext_zf):
                                #echo_msg('Extracting {}'.format(ext_zf))
                                #pbar.update()
                                _dirname = os.path.dirname(ext_zf)
                                if not os.path.exists(_dirname):
                                    os.makedirs(_dirname)

                                with open(ext_zf, 'wb') as f:
                                    f.write(z.read(zf))

                                if verbose:
                                    echo_msg(f'Extracted {ext_zf}')
        except Exception as e:
            echo_error_msg(
                f'could not process ZIP file {src_file}, {e}'
            )
            return([])
                                
    elif src_file.split('.')[-1] == 'gz':
        try:
            tmp_proc = gunzip(src_file, outdir=outdir)
        except:
            if verbose:
                echo_error_msg(f'unable to gunzip {src_file}')
                
            tmp_proc = None

        if tmp_proc is not None:
            for ext in exts:
                if ext == tmp_proc.split('.')[-1]:
                    src_procs.append(tmp_proc)
                    break
                
                else:
                    remove_glob(tmp_proc)
                    
    else:
        for ext in exts:
            if ext == src_file.split('.')[-1]:
                src_procs.append(src_file)
                break
            
    return(src_procs)


def p_f_unzip(src_file, fns=None, outdir='./', tmp_fn=False, verbose=True):
    """unzip/gunzip src_file based on `fn`
    
    Args:
      src_file (str): a zip/gzip filename string
      exts (list): a list of files to extract

    Returns:
      list: a list of the extracted files
    """
    
    src_procs = []
    if src_file.split('.')[-1].lower() == 'zip':
        with zipfile.ZipFile(src_file) as z:
            zfs = z.namelist()
            for fn in fns:
                for zf in zfs:
                    if fn in os.path.basename(zf):
                        #src_procs.append(os.path.join(outdir, os.path.basename(zf)))
                        out_fn = os.path.join(outdir, os.path.basename(zf))
                        if not zf.endswith('/'):
                            if tmp_fn:
                                out_fn = make_temp_fn(zf, temp_dir = outdir)
                                
                            with open(out_fn, 'wb') as f:
                                f.write(z.read(zf))
                                
                        src_procs.append(out_fn)

    elif src_file.split('.')[-1] == 'gz':
        try:
            tmp_proc = gunzip(src_file)
        except:
            echo_error_msg(f'unable to gunzip {src_file}')
            tmp_proc = None
            
        if tmp_proc is not None:
            for fn in fns:
                if fn == os.path.basename(tmp_proc):
                    src_procs.append(os.path.basename(tmp_proc))
                    break
                else: remove_glob(tmp_proc)
    else:
        for fn in fns:
            if fn == os.path.basename(src_file):
                src_procs.append(src_file)
                break
            
    return(src_procs)


###############################################################################
##
## srcwin functions
##
###############################################################################
def fix_srcwin(srcwin, xcount, ycount):
    ## geo_transform is considered in grid-node to properly capture the region

    out_srcwin = [x for x in srcwin]
    if srcwin[0] + srcwin[2] > xcount:
        out_srcwin[2] = xcount - srcwin[0]

    if srcwin[1] + srcwin[3] > ycount:
        out_srcwin[3] = ycount - srcwin[1]

    return(tuple(out_srcwin))
    

def chunk_srcwin(n_size=(), n_chunk=10, step=None, verbose=True):
    return([s for s in yield_srcwin(n_size, n_chunk, step, verbose)])


def yield_srcwin(
        n_size=(),
        n_chunk=10,
        step=None,
        msg='chunking srcwin',
        end_msg='chunked srcwin',
        start_at_edge=True,
        verbose=True
):
    """yield source windows in n_chunks at step"""
    
    if step is None:
        step = n_chunk

    n_edge = n_chunk if start_at_edge else step
    x_chunk = n_edge
    #x_chunk = step
    y_chunk = 0
    i_chunk = 0
    x_i_chunk = 0

    assert step > 0    
    with tqdm(
            total=(math.ceil((n_size[0] + (n_chunk-n_edge)) / step) \
                   * math.ceil((n_size[1] +  (n_chunk-n_edge)) / step)),
            desc='{}: {} @ chunk:{}/step:{}...'.format(
                _command_name(), msg, int_or(n_chunk), int_or(step)
            ),
            leave=verbose
    ) as pbar:
        while True:
            y_chunk = n_edge
            #y_chunk = step
            while True:
                this_x_chunk = n_size[1] if x_chunk > n_size[1] else x_chunk
                this_y_chunk = n_size[0] if y_chunk > n_size[0] else y_chunk
                this_x_origin = int(x_chunk - n_chunk)
                this_y_origin = int(y_chunk - n_chunk)
                this_x_origin = 0 if this_x_origin < 0 else this_x_origin
                this_y_origin = 0 if this_y_origin < 0 else this_y_origin
                this_x_size = int(this_x_chunk - this_x_origin)
                this_y_size = int(this_y_chunk - this_y_origin)
                if this_x_size <= 0 or this_y_size <= 0:
                    break

                srcwin = (this_x_origin, this_y_origin, this_x_size, this_y_size)
                pbar.update()
                yield(srcwin)
                
                if y_chunk > (n_size[0]*step):
                    break
                else:
                    y_chunk += step
                    i_chunk += 1

            if x_chunk > (n_size[1]*step):
                break
            else:
                x_chunk += step
                x_i_chunk += 1

                
def buffer_srcwin(srcwin=(), n_size=None, buff_size=0, verbose=True):
    """yield source windows in n_chunks at step"""

    if n_size is None:
        n_size = srcwin

    x_origin = srcwin[0] - buff_size
    x_origin = 0 if x_origin < 0 else x_origin
        
    y_origin = srcwin[1] - buff_size
    y_origin = 0 if y_origin < 0 else y_origin

    x_buff_size = buff_size * 2 if x_origin !=0 else buff_size
    y_buff_size = buff_size * 2 if y_origin !=0 else buff_size
    
    x_size = srcwin[3] + x_buff_size#(buff_size*2)
    x_size = (n_size[1] - x_origin) if (x_origin + x_size) > n_size[1] else x_size
    
    y_size = srcwin[2] + y_buff_size#(buff_size*2)
    y_size = (n_size[0] - y_origin) if (y_origin + y_size) > n_size[0] else y_size
    
    return(int(x_origin), int(y_origin), int(x_size), int(y_size))


def expand_for(arr, shiftx=1, shifty=1):
    arr_b = arr.copy().astype(bool)
    for i in range(arr.shape[0]):
        for j in range(arr.shape[1]):
            if(arr[i,j]):
                i_min, i_max = max(i-shifty, 0), min(i+shifty+1, arr.shape[0])
                j_min, j_max = max(j-shiftx, 0), min(j+shiftx+1, arr.shape[1])
                arr_b[i_min:i_max, j_min:j_max] = True
    return arr_b


###############################################################################
##
## MB-System functions
##
###############################################################################
def mb_inf(src_xyz, src_fmt=168):
    """generate an info (.inf) file from a src_xyz file using MBSystem.

    Args:
      src_xyz (port): an open xyz port or generator
      src_fmt (int): the datalists format of the source data

    Returns:
      dict: xyz infos dictionary mb_inf_parse(inf_file)
    """

    run_cmd('mbdatalist -O -F{} -I{}'.format(src_fmt, src_xyz.name), verbose=False)
    return(mb_inf_parse('{}.inf'.format(src_xyz.name)))


###############################################################################
##
## system cmd verification and configs.
##
## run_cmd - run a system command as a subprocess
## and return the return code and output
##
## yield_cmd - run a system command as a subprocess
## and yield the output
##
###############################################################################
cmd_exists = lambda x: any(os.access(os.path.join(path, x), os.X_OK) \
                           for path in os.environ['PATH'].split(os.pathsep))


def run_cmd(cmd, data_fun=None, verbose=False, cwd='.'):
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
    out = None
    width = int(_terminal_width()) - 55
    with tqdm(
            desc='cmd: `{}...`'.format(cmd.rstrip()[:width]),
            leave=verbose
            #desc='cmd: `{}...`'.format(cmd.rstrip())
    ) as pbar:
        
        if data_fun is not None:
            pipe_stdin = subprocess.PIPE
        else:
            pipe_stdin = None

        #encoding='utf-8',
        #, universal_newlines=True, bufsize=1,
        p = subprocess.Popen(
            cmd, shell=True, stdin=pipe_stdin, stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE, close_fds=True, cwd=cwd
        )

        if data_fun is not None:
            if verbose:
                echo_msg('piping data to cmd subprocess...')

            data_fun(p.stdin)
            p.stdin.close()

        io_reader = io.TextIOWrapper(p.stderr, encoding='utf-8')
        while p.poll() is None:
            err_line = io_reader.readline()
            if verbose:
                if err_line:
                    pbar.write(err_line.rstrip())
                    sys.stderr.flush()

            pbar.update()

        out = p.stdout.read()
        p.stderr.close()
        p.stdout.close()
        if verbose:
            echo_msg(
                'ran cmd {} and returned {}'.format(
                    cmd.rstrip(), p.returncode
                )
            )

        ## todo update so no crash if cmd doesn't exist!
        #if p.returncode != 0:
        #    raise Exception(p.returncode)
        
    return(out, p.returncode)


def yield_cmd(cmd, data_fun=None, verbose=False, cwd='.'):
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
    p = subprocess.Popen(
        cmd, shell=True, stdin=pipe_stdin, stdout=subprocess.PIPE,
        close_fds=True, cwd=cwd
    )
    
    if data_fun is not None:
        if verbose: echo_msg('piping data to cmd subprocess...')
        data_fun(p.stdin)
        p.stdin.close()

    while p.poll() is None:
        line = p.stdout.readline().decode('utf-8')
        if not line:
            break
        else:
            yield(line)
        
    line = p.stdout.read().decode('utf-8')
    p.stdout.close()
    if verbose:
        echo_msg(
            'ran cmd: {} and returned {}.'.format(
                cmd.rstrip(), p.returncode)
        )

        
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
    else:
        return("0".encode())

    
def config_check(chk_config_file=True, chk_vdatum=False,
                 generate_config_file=True, verbose=False):
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

    def check_config_file(ccc):
        pass
    
    #cudem_cmd_config = os.path.join(cudem_cache(), '.cudem_cmd_config.json')
    cudem_cmd_config = os.path.join(os.path.expanduser('~'), '.cudem_cmd_config.json')
    
    if chk_config_file and os.path.exists(cudem_cmd_config):
        with open(cudem_cmd_config, 'r') as ccc:
            _waff_co = json.load(ccc)
    else:
        py_vers = str(sys.version_info[0]),
        host_os = sys.platform
        _waff_co = {}
        _waff_co['platform'] = host_os
        _waff_co['python'] = py_vers[0]
        ae = '.exe' if host_os == 'win32' else ''

        #if chk_vdatum: _waff_co['VDATUM'] = vdatum(verbose=verbose).vdatum_path
        #_waff_co['GDAL'] = cmd_check('gdal_grid{}'.format(ae), 'gdal-config --version').decode()
        _waff_co['GDAL'] = cmd_check(
            'gdal_grid{}'.format(ae), 'gdal_translate --version'
        ).decode()
        _waff_co['GMT'] = cmd_check(
            'gmt{}'.format(ae), 'gmt --version'
        ).decode()
        _waff_co['MBGRID'] = cmd_check(
            'mbgrid{}'.format(ae), 'mbgrid -version 2>&1 | grep Version'
        ).decode()
        _waff_co['LASZIP'] = cmd_check(
            'laszip{}'.format(ae), 'laszip -version 2>&1 | awk \'{print $5}\''
        ).decode()
        _waff_co['HTDP'] = cmd_check(
            'htdp{}'.format(ae),
            'echo 0 | htdp 2>&1' \
            if host_os == 'win32' \
            else 'echo 0 | htdp 2>&1 | grep SOFTWARE | awk \'{print $3}\''
        ).decode()
        _waff_co['ImageMagick'] = cmd_check(
            'mogrify{}'.format(ae), 'mogrify --version | grep Version | awk \'{print $3}\''
        ).decode()
        _waff_co['CUDEM'] = str(cudem.__version__)
        _waff_co['conda'] = os.environ.get('CONDA_DEFAULT_ENV', None)

        for key in _waff_co.keys():
            _waff_co[key] = None if _waff_co[key] == '0' else _waff_co[key]

        echo_msg(json.dumps(_waff_co, indent=4, sort_keys=True))
        if generate_config_file:
            with open(cudem_cmd_config, 'w') as ccc:
                ccc.write(json.dumps(_waff_co, indent=4, sort_keys=True))
            
    return(_waff_co)


###############################################################################
##
## verbosity functions
##
## TODO: add threading and verbosity
##
###############################################################################
def _terminal_width():
    #cols = 40
    #try:
    #    cols = shutil.get_terminal_size().columns
    #    print(cols)
    #except:
    # try:
    #     rows, cols = curses.initscr().getmaxyx()
    # finally:
    #     curses.endwin()

    cols = get_terminal_size_stderr()[0]
    return(cols)


def get_terminal_size_stderr(fallback=(80, 24)):
    """
    Unlike shutil.get_terminal_size, which looks at stdout, this looks at stderr.
    """
        
    try:
        size = os.get_terminal_size(sys.__stderr__.fileno())
    except (AttributeError, ValueError, OSError):
        size = os.terminal_size(fallback)
        
    return(size)


def _init_msg2(msg, prefix_len, buff_len=6):
    width = int(_terminal_width()) - (prefix_len+buff_len)
    msg_len = len(msg)
    if msg_len > (width / 2):
        while msg_len > (width / 2):
            msg_beg = msg[:int(msg_len / 3)]
            msg_end = msg[-1*int(msg_len / 3):]
            msg_len = len(msg_beg + msg_end)

        msg = f'{msg_beg}...{msg_end}'
        return(msg)
    else:
        return(f'{msg}')

    
def _init_msg(msg, prefix_len, buff_len=6):
    width = int(_terminal_width()) - (prefix_len+buff_len)
    try:
        if len(msg) > width:
            return('{}...'.format(msg[:width]))
        else:
            return('{}'.format(msg))
    except:
        return('{}'.format(msg))

    
def echo_warning_msg2(msg, prefix='cudem'):
    """echo warning msg to stderr using `prefix`

    >> echo_warning_msg2('message', 'test')
    test: warning, message
    
    Args:
      msg (str): a message
      prefix (str): a prefix for the message
    """

    #msg = _init_msg(msg, len(prefix) + 9)
    #sys.stderr.flush()
    sys.stderr.write('\x1b[2K\r')
    #msg = _init_msg(msg, len(prefix))
    #sys.stderr.write('{}: \033[33m\033[1mwarning\033[m, {}\n'.format(prefix, msg))
    tqdm.write(
        '{}: \033[33m\033[1mwarning\033[m, {}'.format(prefix, msg),
        file=sys.stderr
    )
    sys.stderr.flush()

    
def echo_error_msg2(msg, prefix='cudem'):
    """echo error msg to stderr using `prefix`

    >> echo_error_msg2('message', 'test')
    test: error, message

    Args:
      msg (str): a message
      prefix (str): a prefix for the message
    """

    #msg = _init_msg(msg, len(prefix) + 7)
    #sys.stderr.flush()
    sys.stderr.write('\x1b[2K\r')
    #msg = _init_msg(msg, len(prefix))
    #sys.stderr.write('{}: \033[31m\033[1merror\033[m, {}\n'.format(prefix, msg))
    tqdm.write(
        '{}: \033[31m\033[1merror\033[m, {}'.format(prefix, msg),
        file=sys.stderr
    )
    #tqdm.write(traceback.format_exc())
    sys.stderr.flush()

    
def echo_msg2(msg, prefix='cudem', nl=True, bold=False):
    """echo `msg` to stderr using `prefix`

    >> echo_msg2('message', 'test')
    test: message
    
    Args:
      msg (str): a message
      prefix (str): a prefix for the message
      nl (bool): append a newline to the message
    """
    
    #msg = _init_msg(msg, len(prefix))
    #sys.stderr.flush()
    sys.stderr.write('\x1b[2K\r')
    #msg = _init_msg(msg, len(prefix))
    if bold:
        #sys.stderr.write('{}: \033[1m{}\033[m{}'.format(prefix, msg, '\n' if nl else ''))
        tqdm.write(
            '{}: \033[1m{}\033[m'.format(prefix, msg),
            file=sys.stderr
        )
    else:
        #sys.stderr.write('{}: {}{}'.format(prefix, msg, '\n' if nl else ''))
        tqdm.write(
            '{}: {}'.format(prefix, msg),
            file=sys.stderr
        )
    sys.stderr.flush()


###############################################################################    
##
## echo message `m` to sys.stderr using
## auto-generated prefix
## lambda runs: echo_msg2(m, prefix = os.path.basename(sys.argv[0]))
##
###############################################################################
echo_msg = lambda m: echo_msg2(m, prefix = os.path.basename(sys.argv[0]))
echo_msg_bold = lambda m: echo_msg2(m, prefix = os.path.basename(sys.argv[0]), bold = True)
echo_msg_inline = lambda m: echo_msg2(m, prefix = os.path.basename(sys.argv[0]), nl = False)

###############################################################################
##
## echo error message `m` to sys.stderr using
## auto-generated prefix
##
###############################################################################
echo_error_msg = lambda m: echo_error_msg2(m, prefix = os.path.basename(sys.argv[0]))
echo_warning_msg = lambda m: echo_warning_msg2(m, prefix = os.path.basename(sys.argv[0]))

###############################################################################
##
## echo cudem module options
## modules are a dictionary with the module name
## as the key and at least a 'class' key which
## points to the class/function to call for the module
## uses <class>.__doc__ as description
##
## e.g.
## _cudem_module_long_desc({'module_name': {'class': MyClass}})
##
###############################################################################
_cudem_module_short_desc = lambda m: ', '.join(
    ['{}'.format(key) for key in m])
_cudem_module_name_short_desc = lambda m: ',  '.join(
    ['{} ({})'.format(m[key]['name'], key) for key in m])
_cudem_module_long_desc = lambda m: '{cmd} modules:\n% {cmd} ... <mod>:key=val:key=val...\n\n  '.format(cmd=os.path.basename(sys.argv[0])) + '\n  '.join(
    ['\033[1m{:14}\033[0m{}\n'.format(str(key), m[key]['call'].__doc__) for key in m]) + '\n'

_cudem_module_md_table = lambda m: '\n'.join(['| {:14} | {}'.format(str(key), m[key]['call'].__doc__) for key in m])


def echo_modules(module_dict, key):
    if key is None:
        sys.stderr.write(_cudem_module_long_desc(module_dict))

    else:
        if key in module_dict.keys():
            sys.stderr.write(
                _cudem_module_long_desc(
                    {k: module_dict[k] for k in (key,)}
                )
            )
        else:
            sys.stderr.write(
                'Invalid Module Key: {}\nValid Modules: {}\n'.format(
                    key, _cudem_module_short_desc(module_dict)
                )
            )

    sys.stderr.flush()

_command_name = lambda: os.path.basename(sys.argv[0])

###############################################################################
##
## Progress indicator...
##
## Depreciated, just using tqdm now...
##
###############################################################################
class CliProgress():
    """Cudem Absolute Minimum Progress Indicator

    Minimal progress indication to use with CLI.

    e.g.
    import time
    >>> with CliProgress(message='running thing 1', total=10) as pbar:
    ...     for i in range(10):
    ...         time.sleep(2)
    ...         pbar.update()

             running thing 1
    [  ok  ] running thing 1


    or manually:

    >>> pbar = CliProgress(message='running thing 2', total=10)
             running thing 2
    >>> for i in range(10):
    ...     time.sleep(2)
    ...     pbar.update()
    [   ***] running thing 2
    >>> pbar.end(message='ran thing 2')
    [  ok  ] ran thing 2
    True
    """

    def __init__(self, message=None, end_message=None, total=0,
                 sleep=2, verbose=True): 
        self.thread = threading.Thread(target=self.updates)
        self.thread_is_alive = False
        self.tw = 7
        self.count = 0
        self.pc = self.count % self.tw
        self.p = 0
        self.verbose = verbose
        self.sleep = int_or(sleep)

        self.message = message
        self.end_message = end_message
        self.total = total
        self._init_opm()
        
        self.add_one = lambda x: x + 1
        self.sub_one = lambda x: x - 1
        self.spin_way = self.add_one
        self.spinner = [
            '*     ',
            '**    ',
            '***   ',
            ' ***  ',
            '  *** ',
            '   ***',
            '    **',
            '     *'
        ]
        
        self.perc = lambda p: ((p[0]/p[1]) * 100.)
        
        if self.opm is not None and self.verbose:
            self._clear_stderr()
            sys.stderr.write(
                '\r {}  {}\n'.format(" " * (self.tw - 1), self.opm)
            )

            
    def updates(self):
        while self.thread_is_alive:
            if self.sleep is not None:
                time.sleep(self.sleep)
                self.update(p=None)

                
    def __enter__(self):
        if self.verbose:
            try:
                self.thread_is_alive = True
                self.thread.daemon = True
                self.thread.start()
                #while True:
                #    time.sleep(1)
            except (KeyboardInterrupt, SystemExit) as e:
                self.thread_is_alive = False
                self.__exit__(*sys.exc_info())
                sys.exit(-1)
            
            return(self)
        else:
            return(self)

        
    def __exit__(self, exc_type, exc_value, exc_traceback):
        
        if self.verbose:
            self.thread_is_alive = False
            self.thread.join()
            return(self.end(status=exc_value,
                            exc_type=exc_type,
                            exc_value=exc_value,
                            exc_traceback=exc_traceback))
        else:
            return(True)

    def write(self, msg, nl=True, bold=False):
        sys.stderr.write('\x1b[2K\r')
        if bold:
            sys.stderr.write(
                '\033[1m{}\033[m{}'.format(msg, '\n' if nl else '')
            )
        else:
            sys.stderr.write('{}{}'.format(msg, '\n' if nl else ''))
        sys.stderr.flush()

        
    def _init_opm(self):
        self.width = int(self._terminal_width()) - (self.tw+6)
        if len(self.message) > self.width:
            self.opm = '{}...'.format(self.message[:self.width])
        else:
            self.opm = '{}'.format(self.message)

            
    def _terminal_width(self):
        cols = 40
        try:
            cols = shutil.get_terminal_size().columns
        except:
            try:
                rows, cols = curses.initscr().getmaxyx()
            finally:
                curses.endwin()
                
        return(cols)

    
    def _switch_way(self):
        self.spin_way = self.sub_one if self.spin_way == self.add_one else self.add_one

        
    def _clear_stderr(self, slen = 79):
        sys.stderr.write('\x1b[2K\r')
        sys.stderr.flush()

        
    def update(self, p=1, msg=None):
        if self.verbose:
            self._init_opm()
            self.pc = (self.count % self.tw)
            self.sc = (self.count % (self.tw + 1))

            self._clear_stderr()

            if p is not None:
                self.p += p

            if p is not None and self.p < self.total:
                sys.stderr.write('\r[\033[36m{:^5.2f}%\033[m] {}\r'.format(
                    ((self.p/self.total) * 100.), msg if msg is not None else self.opm
                ))
            else:
                sys.stderr.write('\r[\033[36m{:6}\033[m] {}\r'.format(
                    self.spinner[self.sc], msg if msg is not None else self.opm
                ))

            if self.count == self.tw: self.spin_way = self.sub_one
            if self.count == 0: self.spin_way = self.add_one
            self.count = self.spin_way(self.count)
            sys.stderr.flush()

            
    def end(self, status=0, message=None, exc_type=None,
            exc_value=None, exc_traceback=None):
        self._init_opm()
        self._clear_stderr()
        if message is None:
            if self.end_message is None:
                message = self.message
            else:
                message = self.end_message
            
        if status != 0 and status is not None:
            sys.stderr.write(
                '\r[\033[31m\033[1m{:^6}\033[m] {}\n'.format('fail', message)
            )
            raise(exc_value)
        else:
            sys.stderr.write(
                '\r[\033[32m\033[1m{:^6}\033[m] {}\n'.format('ok', message)
            )
        sys.stderr.flush()

        return(True)

    
def physical_cpu_count():
    """On this machine, get the number of physical cores.

    Not logical cores (when hyperthreading is available), but actual physical cores.
    Things such as multiprocessing.cpu_count often give us the logical cores, which
    means we'll spin off twice as many processes as really helps us when we're
    multiprocessing for performance. We want the physical cores."""

    import multiprocessing as mp
    if sys.platform == "linux" or sys.platform == "linux2":
        # On linux. The "linux2" variant is no longer used but here for backward-compatibility.
        lines = os.popen('lscpu').readlines()
        line_with_sockets = [l for l in lines if l[0:11] == "Socket(s): "][0]
        line_with_cps = [l for l in lines if l[0:20] == "Core(s) per socket: "][0]

        num_sockets = int(line_with_sockets.split()[-1])
        num_cores_per_socket = int(line_with_cps.split()[-1])

        return num_sockets * num_cores_per_socket

    elif sys.platform == "darwin":
        # On a mac
        # TODO: Flesh this out from https://stackoverflow.com/questions/
        # 12902008/python-how-to-find-out-whether-hyperthreading-is-enabled
        return mp.cpu_count()

    elif sys.platform == "win32" or sys.platform == "win64" or sys.platform == "cygwin":
        # On a windows machine.
        # TODO: Flesh this out from https://stackoverflow.com/questions/
        # 12902008/python-how-to-find-out-whether-hyperthreading-is-enabled
        return mp.cpu_count()

    else:
        # If we don't know what platform they're using, just default to the cpu_count()
        # It will only get logical cores, but it's better than nothing.
        return mp.cpu_count()

# A dictionary for converting numpy array dtypes into carray identifiers.
# For integers & floats, does not hangle character/string arrays.
# Reference: https://docs.python.org/3/library/array.html
dtypes_dict = {np.int8:    'b',
               np.uint8:   'B',
               np.int16:   'h',
               np.uint16:  'H',
               np.int32:   'l',
               np.uint32:  'L',
               np.int64:   'q',
               np.uint64:  'Q',
               np.float32: 'f',
               np.float64: 'd',
               # Repeat for these expressions of dtype as well.
               np.dtype('int8'):    'b',
               np.dtype('uint8'):   'B',
               np.dtype('int16'):   'h',
               np.dtype('uint16'):  'H',
               np.dtype('int32'):   'l',
               np.dtype('uint32'):  'L',
               np.dtype('int64'):   'q',
               np.dtype('uint64'):  'Q',
               np.dtype('float32'): 'f',
               np.dtype('float64'): 'd'}


## test
def _err_fit_plot(
        xdata, ydata, out, fitfunc, bins_final, std_final,
        sampling_den, max_int_dist, dst_name='unc',
        xa='distance'
):
    """plot a best fit plot with matplotlib

    Args:
      xdata (list): list of x-axis data
      ydata (list): list of y-axis data

    """

    #try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.offsetbox import AnchoredText

    short_name="All Terrain"

    fig = plt.figure()
    ax = plt.subplot(111)

    plt_data=ax.scatter(
        bins_final,
        std_final,
        zorder=1,
        label='Error St. Dev.',
        marker="o",
        color="black",
        s=30
    )
    #plt_best_fit,=ax.plot(xdata,ydata, zorder=1, linewidth=2.0)
    plt_best_fit,=ax.plot(xdata, fitfunc(out, xdata), '-')

    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])

    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom='on',      # ticks along the bottom edge are off
        top='off',         # ticks along the top edge are off
        labelbottom='on') # labels along the bottom edge are off

    anchored_text = AnchoredText(short_name, loc=2)
    anchored_text.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    ax.add_artist(anchored_text)

    anchored_text2 = AnchoredText(" $y = {%gx}^{%g}$ "%(out[1],out[2]), loc=1)
    #add r2 value using below
    #anchored_text2 = AnchoredText(" $y = {%gx}^{%g}$
    #$r^2=%g$ "%(coeff1,coeff2,rsquared), loc=1)
    anchored_text2.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    ax.add_artist(anchored_text2)

    str_ss_samp_den="Sampling Density = " + str(sampling_den) + " %"
    anchored_text3 = AnchoredText(str_ss_samp_den, loc=4)
    anchored_text3.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    ax.add_artist(anchored_text3)

    plt.legend(
        [plt_data, plt_best_fit],
        ['Interpolation Error St Dev', 'Best-Fit Line'],
        loc='upper center',
        bbox_to_anchor=(0.5, 1.15),
        fancybox=True,
        shadow=True,
        ncol=2,
        fontsize=14)

    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom='on',      # ticks along the bottom edge are off
        top='off',         # ticks along the top edge are off
        labelbottom='on') # labels along the bottom edge are off

    plt.xlabel('Distance from Measurement (cells)', fontsize=14)
    plt.ylabel('Interpolation Error St Dev (m)', fontsize=14)
    plt.xlim(xmin=0)
    plt.xlim(xmax=int(max_int_dist)+1)
    plt.ylim(ymin=0)
    y_max=max(std_final)+(0.25*max(std_final))
    plt.ylim(ymax=y_max)

    # plt.plot(xdata, ydata, 'o')
    # plt.plot(xdata, fitfunc(out, xdata), '-')
    #plt.xlabel(xa)
    #plt.ylabel('Interpolation Error (m)')
    out_png = '{}_bf.png'.format(dst_name)
    plt.savefig(out_png)
    plt.close()

    #except: utils.echo_error_msg('you need to install matplotlib to run uncertainty plots...')

    
def _err_scatter_plot(error_arr, dist_arr, mean, std, max_int_dist,
                      bins_orig, sampling_den, dst_name='unc',
                      xa='distance'):
    """plot a scatter plot with matplotlib

    Args:
      error_arr (array): an array of errors
      dist_arr (array): an array of distances

    """

    #try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.offsetbox import AnchoredText

    short_name="All Terrain"

    fig = plt.figure()
    ax = plt.subplot(111)
    plt_data=ax.scatter(
        dist_arr,
        error_arr,
        zorder=1,
        label="Measurements",
        marker=".",
        color="black",
        s=20
    )
    plt_data_uncert=ax.errorbar(bins_orig, mean, yerr=std, fmt='r-', linewidth=3)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])

    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom='on',      # ticks along the bottom edge are off
        top='off',         # ticks along the top edge are off
        labelbottom='on') # labels along the bottom edge are off

    anchored_text = AnchoredText(short_name, loc=2)
    anchored_text.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    ax.add_artist(anchored_text)

    str_ss_samp_den="Sampling Density = " + str(sampling_den) + " %"
    anchored_text3 = AnchoredText(str_ss_samp_den, loc=4)
    anchored_text3.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    ax.add_artist(anchored_text3)

    plt.legend(
        [plt_data, plt_data_uncert],
        ["Interpolation Error", "Mean +/- St. Deviation"],
        loc="upper center",
        bbox_to_anchor=(0.5, 1.15),
        fancybox=True,
        shadow=True,
        ncol=2,
        fontsize=14
    )
    plt.xlabel("Distance from Measurement (cells)", fontsize=14)
    plt.ylabel("Interpolation Error (m)", fontsize=14)
    plt.xlim(xmin=0)
    plt.xlim(xmax=int(max_int_dist)+1)

    #plt.xlabel(xa)
    #plt.ylabel('Interpolation Error (m)')
    out_png = "{}_scatter.png".format(dst_name)
    plt.savefig(out_png)
    plt.close()

    #xcept: utils.echo_error_msg('you need to install matplotlib to run uncertainty plots...')

    
def _errbin(err_arr, nbins=10):
    """calculate and plot the error coefficient given err_arr which is 
    a 2 col array with `err dist

    Args:
      error_arr (array): an array of errors and distances

    Returns:
      list: [coefficient-list]
    """

    error = err_arr[:,0]
    distance = err_arr[:,1]
    
    n, _ = np.histogram(distance, bins = nbins)
    while 0 in n:
        nbins -= 1
        n, _ = np.histogram(distance, bins=nbins)

    #echo_msg('histogram: {}'.format(n))
    serror, _ = np.histogram(distance, bins=nbins, weights=error)
    serror2, _ = np.histogram(distance, bins=nbins, weights=error**2)

    mean = serror / n
    #echo_msg('mean: {}'.format(mean))
    std = np.sqrt(serror2 / n - mean * mean)
    ydata = np.insert(std, 0, 0)
    bins_orig=(_[1:] + _[:-1]) / 2

    xdata = np.insert(bins_orig, 0, 0)
    #xdata[xdata - 0 < 0.0001] = 0.0001
    #while len(xdata) < 3:
    #    xdata = np.append(xdata, 0)
    #    ydata = np.append(ydata, 0)

    return(xdata, ydata)


def _err2coeff(err_arr, sampling_den, coeff_guess=[0, 0.1, 0.2],
               dst_name='unc', xa='distance', plots=False):
    """calculate and plot the error coefficient given err_arr which is 
    a 2 col array with `err dist

    Args:
      error_arr (array): an array of errors and distances

    Returns:
      list: [coefficient-list]
    """

    from scipy import optimize

    error = err_arr[:,0]
    distance = err_arr[:,1]
    max_err = np.max(error)
    min_err = np.min(error)
    #echo_msg('min/max error: {}/{}'.format(min_err, max_err))
    max_int_dist = int(np.max(distance))
    #echo_msg('max dist: {}'.format(max_int_dist))
    nbins = 10
    n, _ = np.histogram(distance, bins = nbins)
    while 0 in n:
        nbins -= 1
        n, _ = np.histogram(distance, bins=nbins)

    #echo_msg('histogram: {}'.format(n))
    serror, _ = np.histogram(distance, bins=nbins, weights=error)
    serror2, _ = np.histogram(distance, bins=nbins, weights=error**2)

    mean = serror / n
    #echo_msg('mean: {}'.format(mean))
    std = np.sqrt(serror2 / n - mean * mean)
    ydata = np.insert(std, 0, 0)
    bins_orig=(_[1:] + _[:-1]) / 2

    xdata = np.insert(bins_orig, 0, 0)
    xdata[xdata - 0 < 0.0001] = 0.0001
    while len(xdata) < 3:
        xdata = np.append(xdata, 0)
        ydata = np.append(ydata, 0)

    #echo_msg('xdata: {} ydata: {}'.format(xdata, ydata))
    #fitfunc = lambda p, x: (p[0] + p[1]) * (x ** p[2])
    fitfunc = lambda p, x: p[0] + p[1] * (x ** p[2])
    errfunc = lambda p, x, y: y - fitfunc(p, x)
    out, cov, infodict, mesg, ier = optimize.leastsq(
        errfunc, coeff_guess, args=(xdata, ydata), full_output=True
    )

    #echo_msg('{}, {}, {}, {}, {}'.format(out, cov, infodict, mesg, ier))
    
    if plots:
        try:
            utils.echo_msg('plotting error data')
            _err_fit_plot(
                xdata,
                ydata,
                out,
                fitfunc,
                bins_orig,
                std,
                sampling_den,
                max_int_dist,
                dst_name,
                xa
            )
            _err_scatter_plot(
                error,
                distance,
                mean,
                std,
                max_int_dist,
                bins_orig,
                sampling_den,
                dst_name,
                xa
            )
        except Exception as e:
           echo_error_msg(
               f'unable to generate error plots, please check configs. {e}'
           )

    return(out)


### End
