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
import gzip
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

def append_fn(bn, src_region, inc, version=1):
    return('{}{}_{}_{}v{}'.format(bn, inc2str(inc), src_region.format('fn'), this_year(), version))

def fn_basename(fn, ext):
    if '.' in ext:
        return(fn[:-len(ext)])
    else:
        return(fn[:-(len(ext)+1)])

def fn_url_p(fn):
    url_sw = ['http://', 'https://', 'ftp://', 'ftps://']
    for u in url_sw:
        try:
            if fn.startswith(u):
                return(True)
        except: return(False)
    return(False)
    
def inc2str(inc):
    """convert a WGS84 geographic increment to a str_inc (e.g. 0.0000925 ==> `13`)

    Args:
      inc (float): a gridding increment

    Returns:
      str: a str representation of float(inc)
    """
    
    import fractions
    return(str(fractions.Fraction(str(inc * 3600)).limit_denominator(10)).replace('/', ''))

def str2inc(inc_str):
    """convert a GMT-style `inc_str` (6s) to geographic units

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
        dict_args[p_arg[0]] = False if p_arg[1].lower() == 'false' else True if p_arg[1].lower() == 'true' else None if p_arg[1].lower() == 'none' else '='.join(p_arg[1:]) if len(p_arg) > 2 else p_arg[1]
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
                    remove_glob('{}/*'.format(g))
                    os.removedirs(g)
                else: os.remove(g)
        except: pass
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
    a = math.sin(dx / 2) * math.sin(dx / 2) + math.cos(math.radians(x0)) * math.cos(math.radians(x1)) * math.sin(dy / 2) * math.sin(dy / 2)
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))
    return(rad_m * c)
    
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

    out_x = geoTransform[0] + (int(in_x + 0.5)*geoTransform[1]) + (int(in_y + 0.5)*geoTransform[2])
    out_y = geoTransform[3] + (int(in_x + 0.5)*geoTransform[4]) + (int(in_y + 0.5)*geoTransform[5])

    #out_x = geoTransform[0] + int(in_x + 0.5) * geoTransform[1] + int(in_y + 0.5) * geoTransform[2]
    #out_y = geoTransform[3] + int(in_x + 0.5) * geoTransform[4] + int(in_y + 0.5) * geoTransform[5]

    return(out_x, out_y)

def _invert_gt(geoTransform):
    """invert the geotransform
    
    Args:
      geoTransform (list): a geo-transform list describing a raster

    Returns:
      list: a geo-transform list describing a raster
    """
    
    det = (geoTransform[1]*geoTransform[5]) - (geoTransform[2]*geoTransform[4])
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

def wkt2geom(wkt):
    """transform a wkt to an ogr geometry
    
    Args:
      wkt (wkt): a wkt geometry

    Returns:
      ogr-geom: the ogr geometry
    """
    
    return(ogr.CreateGeometryFromWkt(wkt))

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

def gdal_prj_file(dst_fn, epsg):
    """generate a .prj file given an epsg code

    returns 0
    """
    
    with open(dst_fn, 'w') as out:
        out.write(sr_wkt(int(epsg), True))
    return(0)
    
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

## ==============================================
##
## Archives (zip/gzip/etc.)
##
## ==============================================
def unzip(zip_file):
    """unzip (extract) `zip_file`

    Args:
      zip_file (str): a zip file pathname string

    Returns:
      list: a list of extracted file names.
    """
    
    zip_ref = zipfile.ZipFile(zip_file)
    zip_files = zip_ref.namelist()
    zip_ref.extractall()
    zip_ref.close()
    return(zip_files)

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
    
def p_unzip(src_file, exts=None):
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
                for zf in zfs:
                    if ext == zf.split('.')[-1]:
                        #if ext in zf:
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
            for ext in exts:
                if ext == tmp_proc.split('.')[-1]:
                    src_procs.append(os.path.basename(tmp_proc))
                    break
                else: remove_glob(tmp_proc)
    else:
        for ext in exts:
            if ext == src_file.split('.')[-1]:
                src_procs.append(src_file)
                break
    return(src_procs)
    
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

def gdal2gdal(src_dem, dst_fmt='GTiff', epsg=4326, dst_dem=None, co=True):
    """convert the gdal file to gdal using gdal

    return output-gdal-fn"""
    
    if os.path.exists(src_dem):
        if dst_dem is None:
            dst_dem = '{}.{}'.format(os.path.basename(src_dem).split('.')[0], utils.gdal_fext(dst_fmt))
        if not co:
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
    
    if verbose:
        _prog = CliProgress('running cmd: `{}`'.format(cmd.rstrip()))
        
    if data_fun is not None:
        pipe_stdin = subprocess.PIPE
    else: pipe_stdin = None
    
    if verbose:
        p = subprocess.Popen(
            cmd, shell=True, stdin=pipe_stdin, stdout=subprocess.PIPE, close_fds=True
        )
    else:
        p = subprocess.Popen(
            cmd, shell=True, stdin=pipe_stdin, stdout=subprocess.PIPE, stderr=subprocess.PIPE, close_fds=True
        )

    if data_fun is not None:
        if verbose: echo_msg('piping data to cmd subprocess...')
        data_fun(p.stdin)
        p.stdin.close()
    
    while p.poll() is None:
        if verbose:
            time.sleep(2)
            _prog.update()

    out = p.stdout.read()
    if not verbose:
        p.stderr.close()
        
    p.stdout.close()
    if verbose:
        _prog.end(p.returncode, 'ran cmd: `{}` and returned {}.'.format(cmd.rstrip(), p.returncode))
        
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
    else: return("0".encode())

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
    _waff_co['LASZIP'] = cmd_check('laszip{}'.format(ae), 'laszip -version 2>&1 | awk \'{print $5}\'').decode()
    _waff_co['CUDEM'] = str(cudem.__version__)
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
    sys.stderr.write('{}: \033[31m\033[1mwarning\033[m, {}\n'.format(prefix, msg))

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

        self.message = message
        self._init_opm()
        
        self.add_one = lambda x: x + 1
        self.sub_one = lambda x: x - 1
        self.spin_way = self.add_one
        self.spinner = ['*     ', '**    ', '***   ', ' ***  ', '  *** ', '   ***', '    **', '     *']
        
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
            import curses
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
            sys.stderr.write('\r[\033[36m{:^5.2f}%\033[m] {}\r'.format(self.perc(p), msg if msg is not None else self.opm))
        else: self.update()
        
    def update(self, msg = None):
        self._init_opm()
        self.pc = (self.count % self.tw)
        self.sc = (self.count % (self.tw + 1))
            
        self._clear_stderr()
        sys.stderr.write('\r[\033[36m{:6}\033[m] {}\r'.format(self.spinner[self.sc], msg if msg is not None else self.opm))
        
        if self.count == self.tw: self.spin_way = self.sub_one
        if self.count == 0: self.spin_way = self.add_one
        self.count = self.spin_way(self.count)

    def end(self, status, end_msg = None):
        self._init_opm()
        self._clear_stderr()
        if end_msg is None:
            end_msg = self.message
            
        if status != 0:
            sys.stderr.write('\r[\033[31m\033[1m{:^6}\033[m] {}\n'.format('fail', end_msg))
        else:
            sys.stderr.write('\r[\033[32m\033[1m{:^6}\033[m] {}\n'.format('ok', end_msg))
        
### End
