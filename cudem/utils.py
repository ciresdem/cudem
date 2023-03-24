### utils.py - CUDEM utilities and functions
##
## Copyright (c) 2010 - 2023 Regents of the University of Colorado
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
from tqdm import tqdm

import numpy as np
from osgeo import gdal
from osgeo import ogr
from osgeo import osr

import cudem

## ==============================================
##
## General Utility Functions, definitions, etc.
##
## ==============================================

#cudem_cache = os.path.join(os.path.expanduser('~'), '.cudem_cache')
#cudem_cache = os.path.abspath('./.cudem_cache')
cudem_cache = lambda: os.path.abspath('./.cudem_cache')

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

def append_fn(fn, src_region, inc, version=None, year=None, res=None, high_res=False):
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
    """return the basename of fn based on ext"""

    return('.'.join(fn.split('.')[:-1]))

    
def fn_url_p(fn):
    """check if fn is a url"""
    
    url_sw = ['http://', 'https://', 'ftp://', 'ftps://']
    for u in url_sw:
        try:
            if fn.startswith(u):
                return(True)
        except:
            return(False)
        
    return(False)
    
def inc2str(inc):
    """convert a WGS84 geographic increment to a str_inc (e.g. 0.0000925 ==> `13`)

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
    
    if inc_str is None or inc_str.lower() == 'none':
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
            echo_error_msg('could not parse increment {}, {}'.format(inc_str, e))
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
        #this_entry = re.findall(r'[^"\s]\S*|".+?"', arg)
        p_arg = arg.split('=')
        if len(p_arg) > 1:
            dict_args[p_arg[0]] = False if p_arg[1].lower() == 'false' else \
                True if p_arg[1].lower() == 'true' else \
                None if p_arg[1].lower() == 'none' else \
                '='.join(p_arg[1:]) if len(p_arg) > 2 else \
                p_arg[1]
        
    return(dict_args)

def dict2args(in_dict):
    out_args = ''
    for i, key in enumerate(in_dict.keys()):
        out_args += '{}={}{}'.format(key, in_dict[key], ':' if i+1 < len(in_dict.keys()) else '')
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
                    os.removedirs(g)
                else: os.remove(g)
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

def int_or(val, or_val=None):
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

def str_or(instr, or_val=None):
    """return instr if instr is a string, else or_val"""
    
    try:
        return(str(instr).replace('"', ''))
    except:
        return(or_val)

def convert_size(size_bytes):
   if size_bytes == 0:
       return('0B')
   size_name = ('B', 'KB', 'MB', 'GB', 'TB', 'PB', 'EB', 'ZB', 'YB')
   i = int(math.floor(math.log(size_bytes, 1024)))
   p = math.pow(1024, i)
   s = round(size_bytes / p, 2)
   return('{} {}'.format(s, size_name[i]))
    
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

## ==============================================
##
## Geotransform functions
##
## ==============================================
def _geo2pixel(geo_x, geo_y, geo_transform, node='grid'):
    """convert a geographic x,y value to a pixel location of geoTransform

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

def _pixel2geo(pixel_x, pixel_y, geo_transform, node='pixel'):
    """convert a pixel location to geographic coordinates given geoTransform

    Args:
      pixel_x (int): the x pixel value
      pixel_y (int): the y pixel value
      geo_transform (list): a geo-transform list describing a raster
    
    Returns:
      list: [geographic-x, geographic-y]
    """
    
    geo_x, geo_y = _apply_gt(pixel_x, pixel_y, geo_transform, node)
    return(geo_x, geo_y)

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
    
## ==============================================
##
## Archives (zip/gzip/etc.)
##
## ==============================================

def unbz2(bz2_file, outdir='./', overwrite=False):

    newfilepath = '.'.join(bz2_file.split('.')[:-1])

    if not os.path.exists(newfilepath):
        echo_msg('Uncompressing {} to {}...'.format(bz2_file, newfilepath))
        with open(newfilepath, 'wb') as new_file, bz2.BZ2File(bz2_file, 'rb') as file:
            for data in iter(lambda : file.read(100 * 1024), b''):
                new_file.write(data)
    return(newfilepath)

def unzip(zip_file, outdir='./', overwrite=False):
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
                    echo_msg('Extracting {}'.format(os.path.join(outdir, fn)))
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

def gunzip(gz_file):
    """gunzip `gz_file`

    Args:
      gz_file (str): a gzip file pathname string.

    Returns:
      str: the extracted file name.
    """
    
    if os.path.exists(gz_file):
        gz_split = gz_file.split('.')[:-1]
        guz_file = '{}.{}'.format(gz_split[0], gz_split[1])
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
                            
def p_unzip(src_file, exts=None, outdir='./'):
    """unzip/gunzip src_file based on `exts`
    
    Args:
      src_file (str): a zip/gzip filename string
      exts (list): a list of extensions to extract

    Returns:
      list: a list of the extracted files
    """
    
    src_procs = []
    if src_file.split('.')[-1].lower() == 'zip':
        with zipfile.ZipFile(src_file) as z:
            zfs = z.namelist()
            for ext in exts:
                with tqdm(total=len(zfs), desc='{}: unzipping {}...'.format(_command_name(), src_file[:14])) as pbar:
                    for zf in zfs:
                        if ext == zf.split('.')[-1]:
                            ext_zf = os.path.join(outdir, zf)
                            src_procs.append(ext_zf)
                            if not os.path.exists(ext_zf):
                                echo_msg('Extracting {}'.format(ext_zf))
                                pbar.update()
                                _dirname = os.path.dirname(ext_zf)
                                if not os.path.exists(_dirname):
                                    os.makedirs(_dirname)

                                with open(ext_zf, 'wb') as f:
                                    f.write(z.read(zf))
                                    
                            else:
                                pbar.update()

                        else:
                            pbar.update()
                                
    elif src_file.split('.')[-1] == 'gz':
        try:
            tmp_proc = gunzip(src_file)
        except:
            echo_error_msg('unable to gunzip {}'.format(src_file))
            tmp_proc = None
            
        if tmp_proc is not None:
            for ext in exts:
                if ext == tmp_proc.split('.')[-1]:
                    src_procs.append(os.path.basename(tmp_proc))
                    break
                
                else:
                    remove_glob(tmp_proc)
                    
    else:
        for ext in exts:
            if ext == src_file.split('.')[-1]:
                src_procs.append(src_file)
                break
            
    return(src_procs)

def p_f_unzip(src_file, fns=None):
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
                    if fn == os.path.basename(zf):
                        src_procs.append(os.path.basename(zf))
                        with open(os.path.basename(zf), 'wb') as f:
                            f.write(z.read(zf))
    elif src_file.split('.')[-1] == 'gz':
        try:
            tmp_proc = gunzip(src_file)
        except:
            echo_error_msg('unable to gunzip {}'.format(src_file))
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

## ==============================================
##
## srcwin functions
##
## ==============================================

def yield_srcwin(n_size=(), n_chunk=10, step=None, verbose=True):
    """yield source windows in n_chunks at step"""
    
    if step is None:
        step = n_chunk
    x_chunk = n_chunk
    y_chunk = 0
    i_chunk = 0
    x_i_chunk = 0

    with tqdm(total=(n_size[0]*n_size[1])/step, desc='{}: chunking srcwin'.format(_command_name)) as pbar:
        while True:
            y_chunk = n_chunk
            while True:
                this_x_chunk = n_size[1] if x_chunk > n_size[1] else x_chunk
                this_y_chunk = n_size[0] if y_chunk > n_size[0] else y_chunk
                this_x_origin = x_chunk - n_chunk
                this_y_origin = y_chunk - n_chunk
                this_x_size = int(this_x_chunk - this_x_origin)
                this_y_size = int(this_y_chunk - this_y_origin)
                if this_x_size == 0 or this_y_size == 0:
                    break
                
                srcwin = (this_x_origin, this_y_origin, this_x_size, this_y_size)
                yield(srcwin)
                
                if y_chunk > n_size[0]:
                    break
                else:
                    y_chunk += step
                    i_chunk += 1
                    
                pbar.update(step)
                
            if x_chunk > n_size[1]:
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

    x_size = srcwin[3] + (buff_size*2)
    x_size = (n_size[1] - x_origin) if (x_origin + x_size) > n_size[1] else x_size
    
    y_size = srcwin[2] + (buff_size*2)
    y_size = (n_size[0] - y_origin) if (y_origin + y_size) > n_size[0] else y_size
    
    return(x_origin, y_origin, x_size, y_size)

## ==============================================
##
## GDAL/OGR/OSR functions
##
## ==============================================

def wkt2geom(wkt):
    """transform a wkt to an ogr geometry
    
    Args:
      wkt (wkt): a wkt geometry

    Returns:
      ogr-geom: the ogr geometry
    """
    
    return(ogr.CreateGeometryFromWkt(wkt))

def sr_wkt(src_srs, esri=False):
    """convert a src_srs to wkt"""
    
    try:
        sr = osr.SpatialReference()
        sr.SetFromUserInput(src_srs)
        if esri:
            sr.MorphToESRI()
            
        return(sr.ExportToWkt())
    except:
        return(None)

def gdal_prj_file(dst_fn, src_srs):
    """generate a .prj file given a src_srs"""
    
    with open(dst_fn, 'w') as out:
        out.write(sr_wkt(src_srs, True))
        
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
        #src_vert_epsg = src_srs.GetAttrValue('VERT_CS|AUTHORITY', 1)
    else:
        src_vert_epsg = None
                
    #src_srs = osr.SpatialReference()
    #src_srs.SetFromUserInput('epsg:{}'.format(src_horz_epsg))

    return(src_horz_epsg, src_vert_epsg)
    
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

def gdal_write(
        src_arr,
        dst_gdal,
        ds_config,
        dst_fmt='GTiff',
        max_cache=False,
        verbose=False
):
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
        #print(src_arr.size)
        #print(src_arr.shape)
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
            dst_dem = '{}.{}'.format(os.path.basename(src_dem).split('.')[0], gdal_fext(dst_fmt))
            
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

## ==============================================
##
## MB-System functions
##
## ==============================================
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
    
    with tqdm(desc='{}: `{}...`'.format(_command_name(), cmd.rstrip()[:14])) as pbar:
        if data_fun is not None:
            pipe_stdin = subprocess.PIPE
        else:
            pipe_stdin = None

        p = subprocess.Popen(
            cmd, shell=True, stdin=pipe_stdin, stdout=subprocess.PIPE, #encoding='utf-8',
            stderr=subprocess.PIPE, close_fds=True#, universal_newlines=True, bufsize=1,
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
            echo_msg('ran cmd {} and returned {}'.format(cmd.rstrip(), p.returncode))
        
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
    p = subprocess.Popen(
        cmd, shell=True, stdin=pipe_stdin, stdout=subprocess.PIPE, close_fds=True
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
        echo_msg('ran cmd: {} and returned {}.'.format(cmd.rstrip(), p.returncode))

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
    
    py_vers = str(sys.version_info[0]),
    host_os = sys.platform
    _waff_co = {}
    _waff_co['platform'] = host_os
    _waff_co['python'] = py_vers[0]
    ae = '.exe' if host_os == 'win32' else ''

    #if chk_vdatum: _waff_co['VDATUM'] = vdatum(verbose=verbose).vdatum_path
    _waff_co['GDAL'] = cmd_check('gdal_grid{}'.format(ae), 'gdal-config --version').decode()
    _waff_co['GMT'] = cmd_check('gmt{}'.format(ae), 'gmt --version').decode()
    _waff_co['MBGRID'] = cmd_check('mbgrid{}'.format(ae), 'mbgrid -version 2>&1 | grep Version').decode()
    _waff_co['LASZIP'] = cmd_check('laszip{}'.format(ae), 'laszip -version 2>&1 | awk \'{print $5}\'').decode()
    _waff_co['HTDP'] = cmd_check('htdp{}'.format(ae), 'echo 0 | htdp 2>&1 | grep SOFTWARE | awk \'{print $3}\'').decode()
    _waff_co['CUDEM'] = str(cudem.__version__)

    for key in _waff_co.keys():
        _waff_co[key] = None if _waff_co[key] == '0' else _waff_co[key]
            
    return(_waff_co)
    
## ==============================================
##
## verbosity functions
##
## TODO: add threading and verbosity
## ==============================================
def _terminal_width():
    cols = 40
    try:
        cols = shutil.get_terminal_size().columns
    except:
        try:
            rows, cols = curses.initscr().getmaxyx()
        finally:
            curses.endwin()
                
    return(cols)

def _init_msg(msg, prefix_len):
    width = int(_terminal_width()) - (prefix_len+6)
    if len(msg) > width:
        return('{}...'.format(msg[:width]))
    else:
        return('{}'.format(msg))

def echo_warning_msg2(msg, prefix = 'waffles'):
    """echo warning msg to stderr using `prefix`

    >> echo_warning_msg2('message', 'test')
    test: warning, message
    
    Args:
      msg (str): a message
      prefix (str): a prefix for the message
    """

    #msg = _init_msg(msg, len(prefix) + 9)
    sys.stderr.flush()
    sys.stderr.write('\x1b[2K\r')
    sys.stderr.write('{}: \033[33m\033[1mwarning\033[m, {}\n'.format(prefix, msg))
    sys.stderr.flush()

def echo_error_msg2(msg, prefix = 'waffles'):
    """echo error msg to stderr using `prefix`

    >> echo_error_msg2('message', 'test')
    test: error, message

    Args:
      msg (str): a message
      prefix (str): a prefix for the message
    """

    #msg = _init_msg(msg, len(prefix) + 7)
    sys.stderr.flush()
    sys.stderr.write('\x1b[2K\r')
    sys.stderr.write('{}: \033[31m\033[1merror\033[m, {}\n'.format(prefix, msg))
    sys.stderr.flush()
    
def echo_msg2(msg, prefix='waffles', nl=True, bold=False):
    """echo `msg` to stderr using `prefix`

    >> echo_msg2('message', 'test')
    test: message
    
    Args:
      msg (str): a message
      prefix (str): a prefix for the message
      nl (bool): append a newline to the message
    """

    #msg = _init_msg(msg, len(prefix))
    sys.stderr.flush()
    sys.stderr.write('\x1b[2K\r')
    if bold:
        sys.stderr.write('{}: \033[1m{}\033[m{}'.format(prefix, msg, '\n' if nl else ''))
    else:
        sys.stderr.write('{}: {}{}'.format(prefix, msg, '\n' if nl else ''))
    sys.stderr.flush()
    
## ==============================================
## echo message `m` to sys.stderr using
## auto-generated prefix
## lambda runs: echo_msg2(m, prefix = os.path.basename(sys.argv[0]))
## ==============================================
echo_msg = lambda m: echo_msg2(m, prefix = os.path.basename(sys.argv[0]))
echo_msg_bold = lambda m: echo_msg2(m, prefix = os.path.basename(sys.argv[0]), bold = True)
echo_msg_inline = lambda m: echo_msg2(m, prefix = os.path.basename(sys.argv[0]), nl = False)

## ==============================================
## echo error message `m` to sys.stderr using
## auto-generated prefix
## ==============================================
echo_error_msg = lambda m: echo_error_msg2(m, prefix = os.path.basename(sys.argv[0]))
echo_warning_msg = lambda m: echo_warning_msg2(m, prefix = os.path.basename(sys.argv[0]))

## ==============================================
## echo cudem module options
## modules are a dictionary with the module name
## as the key and at least a 'class' key which
## points to the class/function to call for the module
## uses <class>.__doc__ as description
##
## e.g.
## _cudem_module_long_desc({'module_name': {'class': MyClass}})
## ==============================================
_cudem_module_short_desc = lambda m: ', '.join(
    ['{}'.format(key) for key in m])
_cudem_module_name_short_desc = lambda m: ',  '.join(
    ['{} ({})'.format(m[key]['name'], key) for key in m])
_cudem_module_long_desc = lambda m: '{cmd} modules:\n% {cmd} ... <mod>:key=val:key=val...\n\n  '.format(cmd=os.path.basename(sys.argv[0])) + '\n  '.join(
    ['\033[1m{:14}\033[0m{}\n'.format(str(key), m[key]['class'].__doc__) for key in m]) + '\n'

_command_name = lambda: os.path.basename(sys.argv[0])

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
            sys.stderr.write('Invalid Module Key: {}\nValid Modules: {}\n'.format(key, _cudem_module_short_desc(module_dict)))

    sys.stderr.flush()

## ==============================================
## Progress indicator...
##
## just use tqdm!
## ==============================================
class CliProgress:
    '''cudem minimal progress indicator'''

    def __init__(self, message = None):
        self.tw = 7
        self.count = 0
        self.pc = self.count % self.tw

        self.message = message
        self._init_opm()
        
        self.add_one = lambda x: x + 1
        self.sub_one = lambda x: x - 1
        self.spin_way = self.add_one
        self.spinner = ['*     ', '**    ', '***   ', ' ***  ',
                        '  *** ', '   ***', '    **', '     *']
        
        self.perc = lambda p: ((p[0]/p[1]) * 100.)
        
        if self.opm is not None:
            self._clear_stderr()
            sys.stderr.write('\r {}  {}\n'.format(" " * (self.tw - 1), self.opm))

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

    def update_perc(self, p, msg = None):
        
        if len(p) == 2 and p[0] <= p[1]:
            self._init_opm()
            self._clear_stderr()
            sys.stderr.write('\r[\033[36m{:^5.2f}%\033[m] {}\r'.format(
                self.perc(p), msg if msg is not None else self.opm
            ))
        else:
            self.update()
            
        sys.stderr.flush()
        
    def update(self, msg = None):
        self._init_opm()
        self.pc = (self.count % self.tw)
        self.sc = (self.count % (self.tw + 1))
            
        self._clear_stderr()
        sys.stderr.write('\r[\033[36m{:6}\033[m] {}\r'.format(
            self.spinner[self.sc], msg if msg is not None else self.opm
        ))
        
        if self.count == self.tw: self.spin_way = self.sub_one
        if self.count == 0: self.spin_way = self.add_one
        self.count = self.spin_way(self.count)
        sys.stderr.flush()
    
    def end(self, status, end_msg = None):
        self._init_opm()
        self._clear_stderr()
        if end_msg is None:
            end_msg = self.message
            
        if status != 0:
            sys.stderr.write(
                '\r[\033[31m\033[1m{:^6}\033[m] {}\n'.format('fail', end_msg)
            )
        else:
            sys.stderr.write(
                '\r[\033[32m\033[1m{:^6}\033[m] {}\n'.format('ok', end_msg)
            )
        sys.stderr.flush()

import multiprocessing as mp
import numpy

def physical_cpu_count():
    """On this machine, get the number of physical cores.

    Not logical cores (when hyperthreading is available), but actual physical cores.
    Things such as multiprocessing.cpu_count often give us the logical cores, which
    means we'll spin off twice as many processes as really helps us when we're
    multiprocessing for performance. We want the physical cores."""
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
        # TODO: Flesh this out from https://stackoverflow.com/questions/12902008/python-how-to-find-out-whether-hyperthreading-is-enabled
        return mp.cpu_count()

    elif sys.platform == "win32" or sys.platform == "win64" or sys.platform == "cygwin":
        # On a windows machine.
        # TODO: Flesh this out from https://stackoverflow.com/questions/12902008/python-how-to-find-out-whether-hyperthreading-is-enabled
        return mp.cpu_count()

    else:
        # If we don't know what platform they're using, just default to the cpu_count()
        # It will only get logical cores, but it's better than nothing.
        return mp.cpu_count()

# A dictionary for converting numpy array dtypes into carray identifiers.
# For integers & floats, does not hangle character/string arrays.
# Reference: https://docs.python.org/3/library/array.html
dtypes_dict = {numpy.int8:    'b',
               numpy.uint8:   'B',
               numpy.int16:   'h',
               numpy.uint16:  'H',
               numpy.int32:   'l',
               numpy.uint32:  'L',
               numpy.int64:   'q',
               numpy.uint64:  'Q',
               numpy.float32: 'f',
               numpy.float64: 'd',
               # Repeat for these expressions of dtype as well.
               numpy.dtype('int8'):    'b',
               numpy.dtype('uint8'):   'B',
               numpy.dtype('int16'):   'h',
               numpy.dtype('uint16'):  'H',
               numpy.dtype('int32'):   'l',
               numpy.dtype('uint32'):  'L',
               numpy.dtype('int64'):   'q',
               numpy.dtype('uint64'):  'Q',
               numpy.dtype('float32'): 'f',
               numpy.dtype('float64'): 'd'}
            
### End
