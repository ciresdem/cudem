#!/usr/bin/env python
### gdal_clip.py
##
## Copyright (c) 2018 - 2021 CIRES Coastal DEM Team
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
## clip a gdal grid using ogr polygon
##
### Code:

import os
import sys
from osgeo import gdal
from cudem import utils
from cudem import demfun

_version = '0.0.6'
_usage = '''gdal_clip.py ({}): clip a gdal grid using vector file.

usage: gdal_clip.py [ src_gdal src_ogr [ OPTIONS ] ]

  src_gdal\t\tThe input raster file-name
  src_ogr\t\tThe input vector file-name

 Options:

  -i, --invert\tInvert the clip
  --help\tPrint the usage text
  --version\tPrint the version information

 Examples:
 % gdal_clip.py input.tif clip.shp

CIRES DEM home page: <http://ciresgroups.colorado.edu/coastalDEM>
'''.format(_version)

if __name__ == '__main__':    
    elev = None
    src_ply = None
    want_invert = False
    i = 1
    while i < len(sys.argv):
        arg = sys.argv[i]
        elif arg == '-i' or arg == '--invert':  want_invert = True
        elif arg == '-help' or arg == '--help' or arg == '-h':
            sys.stderr.write(_usage)
            sys.exit(1)
        elif arg == '-version' or arg == '--version':
            sys.stderr.write('{}\n'.format(_version))
            sys.exit(1)
        elif elev is None: elev = arg
        elif src_ply is None: src_ply = arg
        else:
            sys.stderr.write(_usage)
            sys.exit(1)
        i = i + 1

    if elev is None or not os.path.exists(elev):
        sys.stderr.write(_usage)
        utils.echo_error_msg('you must enter a valid input raster file')
        sys.exit(1)
    else:
        output_name = elev[:-4] + '_cut.tif'

    if src_ply is None or not os.path.exists(src_ply):
        sys.stderr.write(_usage)
        utils.echo_error_msg('you must enter a valid input vector file')
        sys.exit(1)

    out, status = demfun.clip(elev, output_name, src_ply=src_ply, invert=want_invert)
    
### End
