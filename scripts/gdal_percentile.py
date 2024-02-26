#!/usr/bin/env python
### gdal_percentile.py
##
## Copyright (c) 2021 - 2024 CIRES Coastal DEM Team
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
## get a percentile from a gdal compatible DEM 
##
### Code:

import sys
from cudem import gdalfun
from cudem import utils

_version = '0.1.0'
_usage = '''gdal_percentile.py ({}): get a percentile from a gdal DEM

usage: gdal_percentile.py [ file ]

 Options:
  file\t\tThe input DEM file-name

  --perc\tPercentile to acquire [95]
  --help\tPrint the usage text
  --version\tPrint the version information

 Examples:
 % gdal_percentile.py input.tif

CIRES DEM home page: <http://ciresgroups.colorado.edu/coastalDEM>
'''.format(_version)

if __name__ == "__main__":

    elev = None
    perc = 95
    i = 1
    while i < len(sys.argv):
        arg = sys.argv[i]
        if arg == '-p' or arg == '-perc' or arg == '--perc':
            perc = float(sys.argv[i+1])
            i += 1
        elif arg == '-help' or arg == '--help' or arg == '-h':
            sys.stderr.write(_usage)
            sys.exit(1)
        elif arg == '-version' or arg == '--version':
            sys.stderr.write('{}\n'.format(_version))
            sys.exit(1)
        elif elev is None: elev = arg
        else:
            print(arg)
            sys.stderr.write(_usage)
            sys.exit(1)
        i = i + 1

    if elev is None:
        sys.stderr.write(_usage)
        utils.echo_error_msg('you must enter an input file')
        sys.exit(1)
  
    print(gdalfun.gdal_percentile(elev, perc=perc))

### End
