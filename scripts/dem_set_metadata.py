#!/usr/bin/env python
### dem_set_metadata.py
##
## Copyright (c) 2022 CIRES Coastal DEM Team
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
## cut a dem to a region
##
### Code:

import os
import sys
from cudem import demfun
from cudem import utils

_version = '0.0.1'
_usage = '''dem_set_metadata.py ({}): set metadata for a CUDEM DEM

usage: dem_set_meatadata.py [ file ]

 Options:
  file\t\tThe input DEM file-name

  -P, --dst_srs\tSet the projection
  -N, --dst_srs\tSet the nodata value
  -D, --dst_srs\tSet the node
  -F, --dst_srs\tSet the output format

  --help\tPrint the usage text
  --version\tPrint the version information

 Notes:
  by default, will update the given DEM, to output a new file, set the -F switch

 Examples:
 % dem_set_metadata.py input.tif -P epsg:4326 -D grid -N -9999 -F NetCDF

CIRES DEM home page: <http://ciresgroups.colorado.edu/coastalDEM>
'''.format(_version)
        
if __name__ == '__main__':    
    elev = None
    dst_srs = 'epsg:4269'
    nodata = -9999
    node = 'pixel'
    fmt = None
    i = 1
    while i < len(sys.argv):
        arg = sys.argv[i]
        if arg == '--dst_srs' or arg == '-P':
            dst_srs = str(sys.argv[i + 1])
            i = i + 1
        elif arg[:2] == '-P':
            dst_srs = str(arg[2:])
        elif arg == '--nodata' or arg == '-N':
            nodata = float(sys.argv[i + 1])
            i = i + 1
        elif arg[:2] == '-N':
            nodata = float(arg[2:])
        elif arg == '--node' or arg == '-D':
            node = str(sys.argv[i + 1])
            i = i + 1
        elif arg[:2] == '-D':
            node = str(arg[2:])
        elif arg == '--fmt' or arg == '-F':
            fmt = str(sys.argv[i + 1])
            i = i + 1
        elif arg[:2] == '-N':
            fmt = str(arg[2:])
        elif arg == '-help' or arg == '--help' or arg == '-h':
            sys.stderr.write(_usage)
            sys.exit(1)
        elif arg == '-version' or arg == '--version':
            sys.stderr.write('{}\n'.format(_version))
            sys.exit(1)
        elif elev is None: elev = arg
        else:
            sys.stderr.write(_usage)
            sys.exit(1)
        i = i + 1

    if elev is None:
        sys.stderr.write(_usage)
        utils.echo_error_msg('you must enter an input file')
        sys.exit(1)
        
    if os.path.exists(elev):

        if fmt is not  None:
            elev = utils.gdal2gdal(elev, dst_fmt=fmt)
                
        demfun.set_nodata(elev, nodata=nodata, convert_array=False)
        demfun.set_srs(elev, dst_srs)
        demfun.set_metadata(elev, node=node, cudem=True)
        
    else: utils.echo_error_msg('{} is not a valid file'.format(elev))
### End
