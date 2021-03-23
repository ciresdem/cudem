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
import requests

## ==============================================
## import gdal/numpy
## ==============================================
import gdal
import ogr
import osr

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
