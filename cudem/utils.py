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
### Commentary:
##
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
import io
import json
import inspect
import fractions
import threading
import multiprocessing as mp
from typing import List, Dict, Optional, Tuple, Any, Union, Generator

import numpy as np

try:
    from tqdm import tqdm
    USE_TQDM = True
except ImportError:
    USE_TQDM = False

## ==============================================
## General Utility Functions
## ==============================================
THIS_DIR, THIS_FILENAME = os.path.split(__file__)
CUDEM_DATA = os.path.join(THIS_DIR, 'data')

cudem_cache = lambda: os.path.abspath('./cudem_cache')

# def cudem_cache() -> str:
#     return os.path.abspath('./cudem_cache')


def set_cache(cache_dir: str):
    if not os.path.exists(cache_dir):
        os.makedirs(cache_dir)

        
def remove_cache(cache_dir: str):
    if os.path.exists(cache_dir):
        remove_glob(cache_dir)

        
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

def append_fn(fn, src_region, inc, version=None, year=None,
              res=None, high_res=False):
    """Append the src_region, inc and version to a string filename."""
    
    res_str = inc2str(inc) if res is None else res
    reg_str = src_region.format('fn' if not high_res else 'fn_full')
    ver_str = 1 if version is None else version
    year_str = this_year() if year is None else year
    
    return f'{fn}{res_str}_{reg_str}_{year_str}v{ver_str}'


def fn_basename(fn, ext):
    """Return the basename of fn based on ext."""
    
    if '.' in ext:
        return fn[:-len(ext)]
    return fn[:-(len(ext) + 1)]


def fn_basename2(fn):
    """Return the basename of fn."""
    
    t = fn.split('.')
    if len(t) > 1:
        return '.'.join(fn.split('.')[:-1])
    return fn


def fn_ext(fn):
    """Return the extension of fn."""
    
    t = fn.split('.')
    if len(t) > 1:
        return t[-1]
    return None


def make_temp_fn(fn, temp_dir=None, region=None, inc=None):
    """Make a temporary unique filename."""
    
    if temp_dir is None:
        temp_dir = cudem_cache()
        
    fn_bn = fn_basename2(os.path.basename(fn))
    fn_et = fn_ext(fn)
    
    if not os.path.exists(temp_dir):
        try:
            os.makedirs(temp_dir)
        except Exception as e:
            echo_error_msg(f'Could not make dirs "{temp_dir}": {e}')

    ext_str = f'.{fn_et}' if fn_et is not None else ''
    
    if region is None:
        timestamp = datetime.datetime.now().strftime('%Y%m%H%M%S%f')
        return os.path.join(temp_dir, f'{fn_bn}_{timestamp}{ext_str}')
    else:
        inc_str = inc2str(inc) if inc is not None else '0'
        return os.path.join(temp_dir, f'{fn_bn}{inc_str}_{region.format("fn_full")}{ext_str}')

    
def fn_url_p(fn):
    """Check if fn is a URL."""
    
    url_sw = ['http://', 'https://', 'ftp://', 'ftps://', '/vsicurl']
    if str_or(fn):
        try:
            for u in url_sw:
                if fn.startswith(u):
                    return True
        except:
            return False
    return False


def inc2str(inc):
    """Convert a WGS84 geographic increment to a string identifier."""
    
    return str(fractions.Fraction(str(inc * 3600)).limit_denominator(10)).replace('/', '')


def str2inc(inc_str):
    """Convert a GMT-style inc_str (e.g. 6s) to geographic units.

    c/s - arc-seconds
    m - arc-minutes    
    """
    
    if inc_str is None or str(inc_str).lower() == 'none' or len(str(inc_str)) == 0:
        return None
        
    inc_str = str(inc_str)
    units = inc_str[-1]
    
    try:
        if units == 'c' or units == 's':
            return float(inc_str[:-1]) / 3600.
        elif units == 'm':
            return float(inc_str[:-1]) / 360.
        else:
            return float(inc_str)
    except ValueError as e:
        echo_error_msg(f'Could not parse increment {inc_str}: {e}')
        return None

    
def this_date():
    """Get current date."""
    
    return datetime.datetime.now().strftime('%Y%m%d%H%M%S')


def this_year():
    """Get current year."""
    
    return datetime.datetime.now().strftime('%Y')


def dl_hash(fn, sha1=False):
    """Calculate a hash of a file."""
    
    buf_size = 65536
    this_hash = hashlib.sha1() if sha1 else hashlib.md5()

    with open(fn, 'rb') as f:
        while True:
            data = f.read(buf_size)
            if not data:
                break
            this_hash.update(data)
            
    return this_hash.hexdigest()


def args2dict(args, dict_args=None):
    """Convert list of arg strings to dict."""
    
    if dict_args is None:
        dict_args = {}
        
    for arg in args:
        p_arg = arg.split('=')
        if len(p_arg) > 1:
            key = p_arg[0]
            val = p_arg[1]
            
            if val.lower() == 'false':
                dict_args[key] = False
            elif val.lower() == 'true':
                dict_args[key] = True
            elif val.lower() == 'none':
                dict_args[key] = None
            elif len(p_arg) > 2:
                dict_args[key] = '='.join(p_arg[1:])
            else:
                dict_args[key] = val
        
    return dict_args


def dict2args(in_dict):
    """Convert a dictionary to an args string."""
    
    out_args = []
    keys = list(in_dict.keys())
    for i, key in enumerate(keys):
        sep = ':' if i + 1 < len(keys) else ''
        out_args.append(f'{key}={in_dict[key]}{sep}')
    return ''.join(out_args)


def remove_glob(*args):
    """Glob `glob_str` and os.remove results."""
    
    for glob_str in args:
        try:
            globs = glob.glob(glob_str)
            for g in globs:
                if os.path.isdir(g):
                    remove_glob(f'{g}/*')
                    remove_glob(f'{g}/.*')
                    os.removedirs(g)
                else:
                    os.remove(g)
        except Exception as e:
            echo_error_msg(e)
            return -1
    return 0


def slugify(value, allow_unicode=False):
    """Convert to ASCII if 'allow_unicode' is False. Convert spaces/dashes to single dashes."""
    
    import unicodedata
    value = str(value)
    if allow_unicode:
        value = unicodedata.normalize('NFKC', value)
    else:
        value = unicodedata.normalize('NFKD', value).encode('ascii', 'ignore').decode('ascii')
        
    value = re.sub(r'[^\w\s-]', '', value.lower())
    return re.sub(r'[-\s]+', '-', value).strip('-_')


def flatten_recursive(nested_list):
    flat_list = []
    for item in nested_list:
        if isinstance(item, list):
            flat_list.extend(flatten_recursive(item))
        else:
            flat_list.append(item)
    return flat_list


def dict_path2abspath(d=None, except_keys=None):
    if d is None:
        d = {}
    if except_keys is None:
        except_keys = []
        
    for key in d.keys():
        if key in except_keys:
            continue
        
        if isinstance(d[key], dict):
            d[key] = dict_path2abspath(d[key], except_keys)
        elif isinstance(d[key], str):
            if os.path.exists(d[key]):
                d[key] = os.path.abspath(d[key])
    return d


def int_or(val, or_val=None):
    """Return val if val is an integer, else return or_val"""
    
    try:
        return int(float_or(val))
    except:
        return or_val

    
def float_or(val, or_val=None):
    """Return val if val is a float, else return or_val"""
    
    try:
        return float(val)
    except:
        return or_val

    
def str_or(instr, or_val=None, replace_quote=True):
    """Return val if val is a string, else return or_val"""
    
    if instr is None:
        return or_val
    try:
        s = str(instr)
        return s.replace('"', '') if replace_quote else s
    except:
        return or_val

    
def list_str(l: list) -> str:
    """
    Format a list of items into a numbered string representation.
    
    Args:
        l (list): The list of items to format.
        
    Returns:
        str: A string with each item on a new line, prefixed by its index.
             Example: "0: Item A\n1: Item B"
    """
    
    return '\n'.join([f'{i}: {x}' for i, x in enumerate(l)])
    
    
def range_pairs(lst):
    return [(lst[i], lst[i+1]) for i in range(len(lst) - 1)]


def ranges2(lst):
    s = e = None
    r = []
    for i in sorted(lst):
        if s is None:
            s = e = i
        elif i == e or i == e + 1:
            e = i
        else:
            r.append((s, e))
            s = e = i
            
    if s is not None:
        r.append((s, e))
    return r


def convert_size(size_bytes):
   if size_bytes == 0:
       return '0B'
   size_name = ('B', 'KB', 'MB', 'GB', 'TB', 'PB', 'EB', 'ZB', 'YB')
   i = int(math.floor(math.log(size_bytes, 1024)))
   p = math.pow(1024, i)
   s = round(size_bytes / p, 2)
   return f'{s} {size_name[i]}'


def euc_dst(pnt0, pnt1):
    """Euclidean distance in meters (approx).

    Points are geographic and result is in meters.
    """
    
    rad_m = 637100
    distance = math.sqrt(sum([(a - b) ** 2 for a, b in zip(pnt0, pnt1)]))
    return rad_m * distance


def hav_dst(pnt0, pnt1):
    """Haversine distance in meters.

    Points are geographic and result is in meters.
    """
    
    x0, y0 = float(pnt0[0]), float(pnt0[1])
    x1, y1 = float(pnt1[0]), float(pnt1[1])
    rad_m = 637100
    dx = math.radians(x1 - x0)
    dy = math.radians(y1 - y0)
    a = math.sin(dx/2)**2 + math.cos(math.radians(x0)) * math.cos(math.radians(x1)) * math.sin(dy/2)**2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a)) 
    return rad_m * c


def wgs_inc2meter(src_inc):
    """Return a wgs increment as meters"""
    
    gds_equator = 111321.543
    degree_to_radian = lambda d: math.pi * (d / 180)
    return math.cos(degree_to_radian(src_inc)) * (gds_equator * src_inc)


def lll(src_lat):
    """Return the lon/lat length in meters"""
    
    gds_equator = 111321.543
    gds_pi = 3.14159265358979323846
    degree_to_radian = lambda d: gds_pi * (d / 180)
    lonl_ = math.cos(degree_to_radian(src_lat)) * gds_equator
    latl_ = gds_equator
    return lonl_, latl_


def touch(fname, times=None):
    """touch a file to make sure it exists"""
    
    with open(fname, 'a'):
        os.utime(fname, times)
    return fname


def count_data_lines(fn):
    """Count lines in a file, skipping blank lines and comments (starting with #).
    Returns 0 if file does not exist.
    """
    
    if not os.path.exists(fn):
        return 0
        
    try:
        with open(fn, 'r') as f:
            return sum(1 for line in f if line.strip() and not line.strip().startswith('#'))
    except Exception:
        return 0


def get_username():
    username = ''
    while not username:
        username = input('username: ')
    return username


def get_password():
    import getpass
    password = ''
    while not password:
        password = getpass.getpass('password: ')
    return password


def get_outliers(in_array: Any, percentile: float = 75, k: float = 1.5, verbose: bool = False):
    """Get outliers using IQR."""
    
    if verbose:
        echo_msg(f'input percentile: {percentile}')

    percentile = max(0, min(100, percentile))
    max_percentile = percentile
    min_percentile = max(0, percentile - 50)

    if verbose:
        echo_msg(f'percentiles: {min_percentile}>>{max_percentile}')

    if np.all(np.isnan(in_array)):
        return np.nan, np.nan
    else:
        perc_max = np.nanpercentile(in_array, max_percentile)
        perc_min = np.nanpercentile(in_array, min_percentile)
        iqr_p = (perc_max - perc_min) * k
        upper_limit = perc_max + iqr_p
        lower_limit = perc_min - iqr_p

    return upper_limit, lower_limit


def num_strings_to_range(*args):
    """Parse args to a number range string (e.g. '1920-2001')."""
    
    dates = []
    for arg in args:
        if str_or(arg) is not None:
            dd = re.findall(r'[-+]?\d*\.?\d+', str(arg))
            for d in dd:
                dates.append(abs(float(d)))
            
    if len(dates) > 0:
        if min(dates) != max(dates):
            return f'{min(dates)}-{max(dates)}'
        else:
            return str(min(dates))
    return None


## ==============================================
## Geotransform functions
## ==============================================
def _geo2pixel(geo_x, geo_y, geo_transform, node='grid'):
    """Convert a geographic x,y value to a pixel location."""
    
    if geo_transform[2] + geo_transform[4] == 0:
        ## Rounding for edges
        #pixel_x = round((geo_x - geo_transform[0]) / geo_transform[1], 4)
        #pixel_y = round((geo_y - geo_transform[3]) / geo_transform[5], 4)

        ## Numpy
        #pixel_x = np.floor((geo_x - geo_transform[0]) / geo_transform[1]).astype(int)
        #pixel_y = np.floor((geo_y - geo_transform[3]) / geo_transform[5]).astype(int)

        ## Query
        pixel_x = (geo_x - geo_transform[0]) / geo_transform[1]
        pixel_y = (geo_y - geo_transform[3]) / geo_transform[5]

        ## Grid-node geo-transform
        if node == 'grid':
            pixel_x += .5
            pixel_y += .5
    else:
        pixel_x, pixel_y = _apply_gt(geo_x, geo_y, _invert_gt(geo_transform))
        
    return int(pixel_x), int(pixel_y)


def __geo2pixel(geo_x, geo_y, geo_transform, node='pixel'):
    """Convert a geographic x,y value to a pixel location of geoTransform

    Note: use _geo2pixel instead
    """
    
    import affine
    forward_transform = affine.Affine.from_gdal(*geo_transform)
    reverse_transform = ~forward_transform
    pixel_x, pixel_y = reverse_transform * (geo_x, geo_y)
    pixel_x, pixel_y = int(pixel_x+0.5), int(pixel_y+0.5)    
    return pixel_x, pixel_y


def _pixel2geo(pixel_x, pixel_y, geo_transform, node='pixel',
               x_precision=None, y_precision=None):
    """Convert a pixel location to geographic coordinates."""
    
    geo_x, geo_y = _apply_gt(pixel_x, pixel_y, geo_transform, node)
    x_precision = int_or(x_precision)
    y_precision = int_or(y_precision)
    
    gx = round(geo_x, x_precision) if x_precision is not None else geo_x
    gy = round(geo_y, y_precision) if y_precision is not None else geo_y    
    return gx, gy


def _apply_gt(in_x, in_y, geo_transform, node='pixel'):
    """Apply geotransform to in_x, in_y."""
    
    offset = 0.5 if node == 'pixel' else 0
    
    out_x = geo_transform[0] + (in_x + offset) * geo_transform[1] + (in_y + offset) * geo_transform[2]
    out_y = geo_transform[3] + (in_x + offset) * geo_transform[4] + (in_y + offset) * geo_transform[5]
    return out_x, out_y


def _invert_gt(geo_transform):
    """Invert the geo_transform."""
    
    det = (geo_transform[1] * geo_transform[5]) - (geo_transform[2] * geo_transform[4])
    if abs(det) < 1e-15:
        return None
    
    inv_det = 1.0 / det
    out_gt = [0.0] * 6
    out_gt[1] = geo_transform[5] * inv_det
    out_gt[4] = -geo_transform[4] * inv_det
    out_gt[2] = -geo_transform[2] * inv_det
    out_gt[5] = geo_transform[1] * inv_det
    out_gt[0] = (geo_transform[2] * geo_transform[3] - geo_transform[0] * geo_transform[5]) * inv_det
    out_gt[3] = (-geo_transform[1] * geo_transform[3] + geo_transform[0] * geo_transform[4]) * inv_det
    return out_gt


def x360(x):
    if x == 0: return -180
    elif x == 360: return 180
    else: return ((x + 180) % 360) - 180

    
## ==============================================
## Archives (zip/gzip/etc.)
## ==============================================
def unbz2(bz2_file, outdir='./', overwrite=False):
    newfilepath = os.path.splitext(bz2_file)[0]
    if not os.path.exists(newfilepath):
        echo_msg(f'Uncompressing {bz2_file} to {newfilepath}...')
        with open(newfilepath, 'wb') as new_file, bz2.BZ2File(bz2_file, 'rb') as file:
            for data in iter(lambda: file.read(100 * 1024), b''):
                new_file.write(data)
    return newfilepath


def zip_list(zip_file):
    try:
        with zipfile.ZipFile(zip_file) as zip_ref:
            return zip_ref.namelist()
    except Exception as e:
        echo_error_msg(e)
        return None

    
def unzip(zip_file, outdir='./', overwrite=False, verbose=True):
    """Unzip (extract) `zip_file`"""
    
    try:
        with zipfile.ZipFile(zip_file) as zip_ref:
            zip_files = zip_ref.namelist()
            if not overwrite:
                for fn in zip_files:
                    if not os.path.exists(os.path.join(outdir, fn)):
                        if verbose:
                            echo_msg(f'Extracting {os.path.join(outdir, fn)}')
                        zip_ref.extract(fn, outdir)
            else:
                zip_ref.extractall(outdir)
                
            if outdir != './':
                zip_files = [os.path.join(outdir, zf) for zf in zip_files]
                
            return zip_files
    except Exception as e:
        echo_error_msg(f'Could not unzip {zip_file}: {e}')
        return None

    
def gunzip(gz_file, outdir='./'):
    """Gunzip (extract) `gz_file`"""
    
    if os.path.exists(gz_file):
        guz_file = os.path.join(outdir, os.path.basename(fn_basename2(gz_file)))
        with gzip.open(gz_file, 'rb') as in_gz, open(guz_file, 'wb') as f:
            shutil.copyfileobj(in_gz, f)
    else:
        echo_error_msg(f'{gz_file} does not exist')
        guz_file = None
    return guz_file


def p_untar(tar_file, exts=None, outdir='./', verbose=True):
    if exts is None: exts = []
    src_procs = []
    with tarfile.open(tar_file, 'r') as tar:
        tar_fns = tar.getnames()
        for ext in exts:
            for tfn in tar_fns:
                if tfn.endswith(ext):
                    ext_tfn = os.path.join(outdir, os.path.basename(tfn))
                    src_procs.append(ext_tfn)
                    if not os.path.exists(ext_tfn):
                        if verbose:
                            echo_msg(f'Extracting {ext_tfn}')
                        t = tar.extractfile(tfn)
                        with open(ext_tfn, 'wb') as f:
                            f.write(t.read())
    return src_procs


def gdb_unzip(src_zip, outdir='./', verbose=True):
    src_gdb = None
    with zipfile.ZipFile(src_zip) as z:
        zfs = z.namelist()
        for zf in zfs:
            if '.gdb' in zf:
                ext_zf = os.path.join(outdir, zf)
                if not os.path.exists(ext_zf) and not ext_zf.endswith('/'):
                    _dirname = os.path.dirname(ext_zf)
                    if not os.path.exists(_dirname):
                        os.makedirs(_dirname)
                    with open(ext_zf, 'wb') as f:
                        f.write(z.read(zf))
                    if verbose:
                        echo_msg(f'Extracted {ext_zf}')
                elif ext_zf.endswith('.gdb/'):
                    src_gdb = ext_zf
    return src_gdb


def p_unzip(src_file, exts=None, outdir='./', verbose=True):
    """Unzip/gunzip src_file based on `exts`"""
    
    if exts is None: exts = []
    src_procs = []
    
    ext_lower = src_file.split('.')[-1].lower()
    if ext_lower == 'zip':
        try:
            with zipfile.ZipFile(src_file) as z:
                zfs = z.namelist()
                for ext in exts:
                    for zf in zfs:
                        if zf.endswith(ext):
                            ext_zf = os.path.join(outdir, zf)
                            src_procs.append(ext_zf)
                            if not os.path.exists(ext_zf):
                                _dirname = os.path.dirname(ext_zf)
                                if not os.path.exists(_dirname):
                                    os.makedirs(_dirname)
                                with open(ext_zf, 'wb') as f:
                                    f.write(z.read(zf))
                                if verbose:
                                    echo_msg(f'Extracted {ext_zf}')
        except Exception as e:
            echo_error_msg(f'Could not process ZIP file {src_file}: {e}')
            return []
    elif ext_lower == 'gz':
        try:
            tmp_proc = gunzip(src_file, outdir=outdir)
        except:
            if verbose: echo_error_msg(f'Unable to gunzip {src_file}')
            tmp_proc = None

        if tmp_proc is not None:
            for ext in exts:
                if tmp_proc.endswith(ext):
                    src_procs.append(tmp_proc)
                    break
                else:
                    remove_glob(tmp_proc)
    else:
        for ext in exts:
            if src_file.endswith(ext):
                src_procs.append(src_file)
                break
    return src_procs


def p_f_unzip(src_file, fns=None, outdir='./', tmp_fn=False, verbose=True):
    """Unzip/gunzip src_file based on `fn`"""
    
    if fns is None: fns = []
    src_procs = []
    ext_lower = src_file.split('.')[-1].lower()
    
    if ext_lower == 'zip':
        with zipfile.ZipFile(src_file) as z:
            zfs = z.namelist()
            for fn in fns:
                for zf in zfs:
                    if fn in os.path.basename(zf):
                        out_fn = os.path.join(outdir, os.path.basename(zf))
                        if not zf.endswith('/'):
                            if tmp_fn:
                                out_fn = make_temp_fn(zf, temp_dir=outdir)
                            with open(out_fn, 'wb') as f:
                                f.write(z.read(zf))
                        src_procs.append(out_fn)
    elif ext_lower == 'gz':
        try:
            tmp_proc = gunzip(src_file)
        except:
            echo_error_msg(f'Unable to gunzip {src_file}')
            tmp_proc = None
        if tmp_proc is not None:
            for fn in fns:
                if fn == os.path.basename(tmp_proc):
                    src_procs.append(os.path.basename(tmp_proc))
                    break
                else:
                    remove_glob(tmp_proc)
    else:
        for fn in fns:
            if fn == os.path.basename(src_file):
                src_procs.append(src_file)
                break
    return src_procs


## ==============================================
## srcwin functions
## ==============================================
def fix_srcwin(srcwin, xcount, ycount):
    """geo_transform is considered in grid-node to properly capture the region"""
    
    out_srcwin = list(srcwin)
    if srcwin[0] + srcwin[2] > xcount:
        out_srcwin[2] = xcount - srcwin[0]
    if srcwin[1] + srcwin[3] > ycount:
        out_srcwin[3] = ycount - srcwin[1]
    return tuple(out_srcwin)


def chunk_srcwin(n_size=(), n_chunk=10, step=None, verbose=True):
    return list(yield_srcwin(n_size, n_chunk, step, verbose=verbose))


def yield_srcwin(n_size=(), n_chunk=10, step=None, msg='chunking srcwin',
                 end_msg='chunked srcwin', start_at_edge=True, verbose=True):
    """yield source windows in n_chunks at step"""
    
    if step is None:
        step = n_chunk

    n_edge = n_chunk if start_at_edge else step
    x_chunk = n_edge
    
    total_steps = (math.ceil((n_size[0] + (n_chunk - n_edge)) / step) * math.ceil((n_size[1] + (n_chunk - n_edge)) / step))
    
    with ccp(total=total_steps, 
             desc=f'{msg} @ chunk:{int_or(n_chunk)}/step:{int_or(step)}...',
             leave=verbose) as pbar:
        while True:
            y_chunk = n_edge
            while True:
                this_x_chunk = min(x_chunk, n_size[1])
                this_y_chunk = min(y_chunk, n_size[0])
                
                this_x_origin = max(0, int(x_chunk - n_chunk))
                this_y_origin = max(0, int(y_chunk - n_chunk))
                
                this_x_size = int(this_x_chunk - this_x_origin)
                this_y_size = int(this_y_chunk - this_y_origin)
                
                if this_x_size <= 0 or this_y_size <= 0:
                    break

                srcwin = (this_x_origin, this_y_origin, this_x_size, this_y_size)
                pbar.update()
                yield srcwin
                
                if y_chunk > (n_size[0] * step):
                    break
                y_chunk += step

            if x_chunk > (n_size[1] * step):
                break
            x_chunk += step

            
def buffer_srcwin(srcwin=(), n_size=None, buff_size=0, verbose=True):
    """Buffer the srcwin by `buff_size`"""
    
    if n_size is None:
        n_size = srcwin

    x_origin = max(0, srcwin[0] - buff_size)
    y_origin = max(0, srcwin[1] - buff_size)

    x_buff_size = buff_size * 2 if x_origin != 0 else buff_size
    y_buff_size = buff_size * 2 if y_origin != 0 else buff_size
    
    x_size = srcwin[3] + x_buff_size
    if (x_origin + x_size) > n_size[1]:
        x_size = n_size[1] - x_origin
    
    y_size = srcwin[2] + y_buff_size
    if (y_origin + y_size) > n_size[0]:
        y_size = n_size[0] - y_origin
    
    return int(x_origin), int(y_origin), int(x_size), int(y_size)


def expand_for(arr, shiftx=1, shifty=1, revert=False):
    arr_b = arr.copy().astype(bool)
    for i in range(arr.shape[0]):
        for j in range(arr.shape[1]):
            if(arr[i,j]):
                i_min, i_max = max(i-shifty, 0), min(i+shifty+1, arr.shape[0])
                j_min, j_max = max(j-shiftx, 0), min(j+shiftx+1, arr.shape[1])
                arr_b[i_min:i_max, j_min:j_max] = True                
    return arr_b


def fill_for(arr, iterations=3):
    filled_arr = np.copy(arr)
    for _ in range(iterations):
        for i in range(1, arr.shape[0] - 1):
            for j in range(1, arr.shape[1] - 1):
                if not filled_arr[i, j]:
                    if (filled_arr[i-1, j] or filled_arr[i+1, j] or
                        filled_arr[i, j-1] or filled_arr[i, j+1]):
                        filled_arr[i, j] = True                        
    return filled_arr_iter


## ==============================================
## MB-System functions
## ==============================================
def mb_inf(src_xyz, src_fmt=168):
    """Generate an info (.inf) file from a src_xyz file using MBSystem."""
    
    try:
        from cudem.mbsfun import mb_inf_parse
    except ImportError:
        echo_error_msg("mb_inf_parse not available.")
        return None

    run_cmd(f'mbdatalist -O -F{src_fmt} -I{src_xyz.name}', verbose=False)
    return mb_inf_parse(f'{src_xyz.name}.inf')


## ==============================================
## System Command Functions
## ==============================================
cmd_exists = lambda x: any(os.access(os.path.join(path, x), os.X_OK) 
                           for path in os.environ['PATH'].split(os.pathsep))


def run_cmd(cmd, data_fun=None, verbose=False, cwd='.'):
    """Run a system command while optionally passing data.

    `data_fun` should be a function to write to a file-port:
    >> data_fun = lambda p: datalist_dump(wg, dst_port = p, ...)
    """
    
    out = None
    cols, _ = shutil.get_terminal_size()
    width = cols - 55
    
    with ccp(desc=f'`{cmd.rstrip()[:width]}...`', leave=verbose) as pbar:
        pipe_stdin = subprocess.PIPE if data_fun is not None else None

        p = subprocess.Popen(
            cmd, shell=True, stdin=pipe_stdin, stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE, close_fds=True, cwd=cwd
        )

        if data_fun is not None:
            if verbose:
                echo_msg('Piping data to cmd subprocess...')
            data_fun(p.stdin)
            p.stdin.close()

        io_reader = io.TextIOWrapper(p.stderr, encoding='utf-8')
        while p.poll() is None:
            err_line = io_reader.readline()
            if verbose and err_line:
                pbar.write(err_line.rstrip())
                sys.stderr.flush()
            pbar.update()

        out = p.stdout.read()
        p.stderr.close()
        p.stdout.close()
        
        if verbose:
            echo_msg(f'Ran cmd {cmd.rstrip()} and returned {p.returncode}')
        
    return out, p.returncode


def yield_cmd(cmd, data_fun=None, verbose=False, cwd='.'):
    """Yield output from a system command.

    `data_fun` should be a function to write to a file-port:
    >> data_fun = lambda p: datalist_dump(wg, dst_port = p, ...)
    """
    
    if verbose: echo_msg(f'Running cmd {cmd.rstrip()}...')
    
    pipe_stdin = subprocess.PIPE if data_fun is not None else None
    
    p = subprocess.Popen(
        cmd, shell=True, stdin=pipe_stdin, stdout=subprocess.PIPE,
        close_fds=True, cwd=cwd
    )
    
    if data_fun is not None:
        if verbose: echo_msg('Piping data to cmd subprocess...')
        data_fun(p.stdin)
        p.stdin.close()

    while p.poll() is None:
        line = p.stdout.readline().decode('utf-8')
        if not line:
            break
        yield line
        
    p.stdout.close()
    if verbose:
        echo_msg(f'Ran cmd {cmd.rstrip()}, returned {p.returncode}.')

        
def cmd_check(cmd_str, cmd_vers_str):
    """check system for availability of 'cmd_str'"""
    
    if cmd_exists(cmd_str): 
        cmd_vers, status = run_cmd(f'{cmd_vers_str}')
        return cmd_vers.rstrip()
    return b"0"


def config_check(chk_config_file=True, chk_vdatum=False,
                 generate_config_file=True, verbose=False):
    """check for needed waffles external software.
    
    waffles external software: gdal, gmt, mbsystem
    also checks python version and host OS and 
    records waffles version
    """
        
    from cudem import __version__
    
    cudem_cmd_config = os.path.join(os.path.expanduser('~'), '.cudem_cmd_config.json')
    
    if chk_config_file and os.path.exists(cudem_cmd_config):
        with open(cudem_cmd_config, 'r') as ccc:
            _waff_co = json.load(ccc)
    else:
        host_os = sys.platform
        _waff_co = {
            'platform': host_os,
            'python': str(sys.version_info[0])
        }
        ae = '.exe' if host_os == 'win32' else ''

        _waff_co['GDAL'] = cmd_check(
            f'gdal_grid{ae}', 'gdal_translate --version'
        ).decode()
        _waff_co['GMT'] = cmd_check(
            f'gmt{ae}', 'gmt --version'
        ).decode()
        _waff_co['MBGRID'] = cmd_check(
            f'mbgrid{ae}', 'mbgrid -version 2>&1 | grep Version'
        ).decode()
        _waff_co['LASZIP'] = cmd_check(
            f'laszip{ae}', "laszip -version 2>&1 | awk '{print $5}'"
        ).decode()
        
        htdp_cmd = 'echo 0 | htdp 2>&1' if host_os == 'win32' else "echo 0 | htdp 2>&1 | grep SOFTWARE | awk '{print $3}'"
        _waff_co['HTDP'] = cmd_check(f'htdp{ae}', htdp_cmd).decode()
        
        _waff_co['ImageMagick'] = cmd_check(
            f'mogrify{ae}', "mogrify --version | grep Version | awk '{print $3}'"
        ).decode()
        _waff_co['CUDEM'] = str(__version__)
        _waff_co['conda'] = os.environ.get('CONDA_DEFAULT_ENV', None)

        for key in _waff_co.keys():
            if _waff_co[key] == '0':
                _waff_co[key] = None

        if verbose:
            echo_msg(json.dumps(_waff_co, indent=4, sort_keys=True))
            
        if generate_config_file:
            with open(cudem_cmd_config, 'w') as ccc:
                ccc.write(json.dumps(_waff_co, indent=4, sort_keys=True))
            
    return _waff_co


## ==============================================
## Verbosity / Progress
## ==============================================
DST_PORT = sys.stderr
MSG_LEVEL = 1
MSG_LEVELS = {'debug': 0, 'info': 1, 'proc': 2, 'warning': 3, 'error': 4}

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

    return get_terminal_size_stderr()[0]


def get_terminal_size_stderr(fallback=(80, 24)):
    """Unlike shutil.get_terminal_size, which looks at stdout, this looks at stderr."""
        
    try:
        size = os.get_terminal_size(sys.__stderr__.fileno())
    except (AttributeError, ValueError, OSError):
        size = os.terminal_size(fallback)        
    return size


def _init_msg2(msg, prefix_len, buff_len=6):
    width = int(_terminal_width()) - (prefix_len+buff_len)
    msg_len = len(msg)
    if msg_len > (width / 2):
        while msg_len > (width / 2):
            msg_beg = msg[:int(msg_len / 3)]
            msg_end = msg[-1*int(msg_len / 3):]
            msg_len = len(msg_beg + msg_end)

        msg = f'{msg_beg}...{msg_end}'
        return msg
    else:
        return f'{msg}'

    
def _init_msg(msg, prefix_len, buff_len=6):
    width = int(_terminal_width()) - (prefix_len+buff_len)
    try:
        if len(msg) > width:
            return f'{msg[:width]}...'
        else:
            return f'{msg}'
    except:
        return f'{msg}'


def echo_msg2(msg, prefix='cudem', level='info', nl='\n', bold=False, use_tqdm=False, dst_port=sys.stderr):
    """echo `msg` to stderr using `prefix`
    
    >> echo_msg2('message', 'test')
    test: message
    """
        
    if level.lower() in MSG_LEVELS and MSG_LEVELS[level.lower()] >= MSG_LEVEL:
        lvl_color = {
            'warning': '\033[33m\033[1mWARNING\033[m',
            'error': '\033[31m\033[1mERROR\033[m',
            'debug': '\033[35m\033[1mDEBUG\033[m',
        }.get(level.lower(), f'\033[36m\033[1m{level.upper()}\033[m')

        if use_tqdm: nl = ''
        
        dst_port.write('\x1b[2K\r')
        
        message = f'\033[1m{msg}\033[m{nl}' if bold else f'{msg}{nl}'
        
        if use_tqdm:
            tqdm.write(f'[ {lvl_color} ] {prefix}: {message}', file=dst_port)
        else:
            dst_port.write(f'[ {lvl_color} ] {prefix}: {message}')
        
        dst_port.flush()

        
def get_calling_module_name(stack_level=1):
    import inspect
    caller_frame = inspect.stack()[stack_level]
    return inspect.getmodulename(caller_frame.filename)


## echo message `m` to sys.stderr using
## auto-generated prefix
## lambda runs: echo_msg2(m, prefix = os.path.basename(sys.argv[0]))
echo_msg = lambda m: echo_msg2(m, prefix=get_calling_module_name(stack_level=2), level='info', use_tqdm=USE_TQDM, dst_port=DST_PORT)
echo_msg_bold = lambda m: echo_msg2(m, prefix=get_calling_module_name(stack_level=2), level='info', bold=True, use_tqdm=USE_TQDM, dst_port=DST_PORT)
echo_msg_inline = lambda m: echo_msg2(m, prefix=get_calling_module_name(stack_level=2), level='info', nl=False, use_tqdm=USE_TQDM, dst_port=DST_PORT)
echo_debug_msg = lambda m: echo_msg2(m, prefix=get_calling_module_name(stack_level=2), level='debug', use_tqdm=USE_TQDM, dst_port=DST_PORT)
echo_error_msg = lambda m: echo_msg2(m, prefix=get_calling_module_name(stack_level=2), level='error', use_tqdm=USE_TQDM, dst_port=DST_PORT)
echo_warning_msg = lambda m: echo_msg2(m, prefix=get_calling_module_name(stack_level=2), level='warning', use_tqdm=USE_TQDM, dst_port=DST_PORT)

## Module Descriptions
## echo cudem module options
## modules are a dictionary with the module name
## as the key and at least a 'class' key which
## points to the class/function to call for the module
## uses <class>.__doc__ as description
## e.g.
## _cudem_module_long_desc({'module_name': {'class': MyClass}})
_cudem_module_short_desc = lambda m: ', '.join([f'{key}' for key in m])
_cudem_module_name_short_desc = lambda m: ',  '.join([f'{m[key]["name"]} ({key})' for key in m])
_cudem_module_long_desc = lambda m: f'{os.path.basename(sys.argv[0])} modules:\n% {os.path.basename(sys.argv[0])} ... <mod>:key=val:key=val...\n\n  ' + '\n  '.join([f'\033[1m{str(key):14}\033[0m{m[key]["call"].__doc__}\n' for key in m]) + '\n'


def echo_modules(module_dict, key):
    if key is None:
        sys.stderr.write(_cudem_module_long_desc(module_dict))
    else:
        if key in module_dict.keys():
            sys.stderr.write(_cudem_module_long_desc({k: module_dict[k] for k in (key,)}))
        else:
            sys.stderr.write(f'Invalid Module Key: {key}\nValid Modules: {_cudem_module_short_desc(module_dict)}\n')
    sys.stderr.flush()

## ==============================================
## Progress indicator...
## ==============================================
class CudemCommonProgress:
    """Cudem Common Progress

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

    def __init__(self, desc=None, message=None, end_message=None, total=0,
                 sleep=2, verbose=True, leave=True, **kwargs):
        """Minimal progress indicator to use with CLI if tqdm isn't available or preferred."""
        
        self.thread = threading.Thread(target=self.updates)
        self.thread_is_alive = False
        self.tw = 7
        self.count = 0
        self.pc = self.count % self.tw
        self.p = 0
        self.verbose = verbose
        self.sleep = int_or(sleep)
        self.message = desc if desc else message
        if self.message:
            self.message = f'{get_calling_module_name(stack_level=2)}: {self.message}'
            
        self.end_message = end_message

        if self.end_message is not None:
            self.end_message = f'{get_calling_module_name(stack_level=2)}: {self.message}'
            
        self.total = total
        self._init_opm()
        
        self.add_one = lambda x: x + 1
        self.sub_one = lambda x: x - 1
        self.spin_way = self.add_one
        self.spinner = ['*     ', '**    ', '***   ', ' ***  ', '  *** ', '   ***', '    **', '     *']
        
        self.perc = lambda p: ((p[0]/p[1]) * 100.)
        
        # if self.opm is not None and self.verbose:
        #     self._clear_stderr()
        #     sys.stderr.write(
        #         '\r {}  {}\n'.format(" " * (self.tw - 1), self.opm)
        #     )

            
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
            return self
        else:
            return self

        
    def __exit__(self, exc_type, exc_value, exc_traceback):
        
        if self.verbose:
            self.thread_is_alive = False
            self.thread.join()
            return(self.end(status=exc_value,
                            exc_type=exc_type,
                            exc_value=exc_value,
                            exc_traceback=exc_traceback))
        else:
            return True

        
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
        import curses
        
        cols = 40
        try:
            cols = shutil.get_terminal_size().columns
        except:
            try:
                rows, cols = curses.initscr().getmaxyx()
            finally:
                curses.endwin()
                
        return cols

    
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
        return True

    
class SimpleProgress:
    """Minimal progress indicator to use with CLI if tqdm isn't available or preferred."""
    
    def __init__(self, desc=None, message=None, total=0, sleep=None, verbose=True, leave=True): 
        self.tw = 7
        self.count = 0
        self.p = 0
        self.verbose = verbose
        self.message = desc if desc else message
        if self.message:
            self.message = f'{get_calling_module_name(stack_level=2)}: {self.message}'
        self.total = total
        
        self.spinner = ['* ', '** ', '*** ', ' *** ', '  *** ', '   ***', '    **', '     *']
        self.direction = 1
        self.width = shutil.get_terminal_size().columns - (self.tw + 6)

        
    def __enter__(self):
        return self

    
    def __exit__(self, exc_type, exc_value, exc_traceback):
        if self.verbose:
            self.end()

            
    def update(self, p=1, msg=None):
        if self.verbose:
            self.p += p if p else 0
            
            ## Truncate msg
            disp_msg = msg if msg else self.message
            if len(disp_msg) > self.width:
                disp_msg = f'{disp_msg[:self.width]}...'

            sys.stderr.write('\x1b[2K\r')
            if self.total > 0:
                perc = (self.p / self.total) * 100
                sys.stderr.write(f'\r[\033[36m{perc:^5.2f}%\033[m] {disp_msg}\r')
            else:
                sc = self.count % len(self.spinner)
                sys.stderr.write(f'\r[\033[36m{self.spinner[sc]:6}\033[m] {disp_msg}\r')
                
            self.count += self.direction
            if self.count == self.tw or self.count == 0:
                self.direction *= -1
            sys.stderr.flush()

            
    def write(self, msg):
        sys.stderr.write(f'\x1b[2K\r{msg}\n')
        sys.stderr.flush()

        
    def end(self):
        sys.stderr.write('\x1b[2K\r')
        sys.stderr.flush()

        
# Define ccp based on TQDM availability
if USE_TQDM:
    class ccp(tqdm):
        def __init__(self, **kwargs):
            mod_name = get_calling_module_name(stack_level=2)
            orig_desc = kwargs.get("desc", "")
            desc_width = shutil.get_terminal_size().columns - 30
            shortened_desc = orig_desc[:desc_width]
            desc = shortened_desc #if len(orig_desc) > desc_width else orig_desc
            kwargs['desc'] = f'[ \033[32m\033[1mPROC\033[m ] {mod_name}: {desc}'
            #kwargs['dynamic_ncols'] = True
            kwargs['ncols'] = shutil.get_terminal_size().columns
            kwargs['bar_format'] = '{l_bar}{bar}{r_bar}'
            kwargs['file'] = DST_PORT
            super().__init__(**kwargs)
else:
    class ccp(CudemCommonProgress):
        pass

    
def physical_cpu_count():
    """Get number of physical cores."""
    
    if sys.platform == "linux":
        try:
            # Fallback for Linux
            lines = os.popen('lscpu').readlines()
            sockets = int([l for l in lines if 'Socket(s):' in l][0].split()[-1])
            cores = int([l for l in lines if 'Core(s) per socket:' in l][0].split()[-1])
            return sockets * cores
        except:
            return mp.cpu_count()
    elif sys.platform == "darwin" or sys.platform == "win32":
        # psutil would be better here, but avoiding extra dependency
        return mp.cpu_count()
    else:
        return mp.cpu_count()

    
## ==============================================
## Error / Analysis Plotting
## ==============================================
def _err_fit_plot(
        xdata, ydata, out, fitfunc, bins_final, std_final,
        sampling_den, max_int_dist, dst_name='unc',
        xa='distance'
):
    """plot a best fit plot with matplotlib."""
    
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
    #except: echo_error_msg('you need to install matplotlib to run uncertainty plots...')

    
def _err_scatter_plot(error_arr, dist_arr, mean, std, max_int_dist,
                      bins_orig, sampling_den, dst_name='unc',
                      xa='distance'):
    """plot a scatter plot with matplotlib."""
    
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
    #except: echo_error_msg('you need to install matplotlib to run uncertainty plots...')

    
def _errbin(err_arr, nbins=10):
    """calculate and plot the error coefficient given err_arr which is 
    a 2 col array with `err dist
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
    return xdata, ydata

def _err2coeff(err_arr, sampling_den, coeff_guess=None,
               dst_name='unc', xa='distance', plots=False):
    """calculate and plot the error coefficient given err_arr which is 
    a 2 col array with `err dist.
    """
    
    try:
        from scipy import optimize
    except ImportError:
        echo_error_msg("Scipy required for _err2coeff")
        return None

    if coeff_guess is None:
        coeff_guess = [0, 0.1, 0.2]

    error = err_arr[:,0]
    distance = err_arr[:,1]
    
    nbins = 10
    n, _ = np.histogram(distance, bins=nbins)
    while 0 in n and nbins > 1:
        nbins -= 1
        n, _ = np.histogram(distance, bins=nbins)

    serror, _ = np.histogram(distance, bins=nbins, weights=error)
    serror2, _ = np.histogram(distance, bins=nbins, weights=error**2)

    mean = serror / n
    std = np.sqrt(serror2 / n - mean * mean)
    
    ydata = np.insert(std, 0, 0)
    bins_orig = (_[1:] + _[:-1]) / 2
    xdata = np.insert(bins_orig, 0, 0)
    
    # Ensure no zero division in power function
    xdata[xdata < 0.0001] = 0.0001

    fitfunc = lambda p, x: p[0] + p[1] * (x ** p[2])
    errfunc = lambda p, x, y: y - fitfunc(p, x)
    
    out, cov, infodict, mesg, ier = optimize.leastsq(
        errfunc, coeff_guess, args=(xdata, ydata), full_output=True
    )

    if plots:
        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            from matplotlib.offsetbox import AnchoredText

            echo_msg('plotting error data')
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
            
        except ImportError:
            echo_error_msg('Matplotlib required for error plots')
        except Exception as e:
            echo_error_msg(f'Unable to generate error plots: {e}')

    return out


### End
