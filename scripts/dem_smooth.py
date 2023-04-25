#!/usr/bin/env python
### dem_smooth.py
##
## Copyright (c) 2012 - 2023 CIRES Coastal DEM Team
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
### Code:

import os
import sys
from cudem import demfun
from cudem import utils

_version = '0.1.1'
_usage = '''dem_smooth.py ({}): smooth a DEM
dem_smooth.py: Smooth a DEM with a Gaussian filter.

usage: dem_smooth.py [ si [ args ] ] [ file ]

Options:
  file\t\tThe input DEM file-name
  -s\t\tValue for the smooth factor; default is 10
  -i\t\tA file containing DEM file-names to process

  -help\t\tPrint the usage text
  -version\tPrint the version information

Example:
dem_smooth.py input.tif -s 12
CIRES DEM home page: <http://ciresgroups.colorado.edu/coastalDEM>'''.format(_version)

if __name__ == '__main__':    
    elev = None
    smooth_factor = 10
    in_list = None

    i = 1
    while i < len(sys.argv):
        arg = sys.argv[i]
        if arg == '-s' or arg == '-smooth' or arg == '--smooth':
            smooth_factor = sys.argv[i+1]
            i = i + 1
        elif arg == '-i':
            in_list = sys.argv[i+1]
            i = i + 1
        elif arg == '-help' or arg == '--help' or arg == '-h':
            sys.stderr.write(_usage)
            sys.exit(1)
        elif arg == '-version' or arg == '--version':
            sys.stderr.write('dem_smooth.py v.%s' %(_version))
            sys.stderr.write(_license)
            sys.exit(1)
        elif elev is None:
            elev = arg
        else:
            sys.stderr.write(_usage)
            sys.exit(1)
        i = i + 1

    if elev is None and in_list is None:
        sys.stderr.write(_usage)
        sys.exit(1)

    try: smooth_factor = int(smooth_factor)
    except:
        utils.echo_error_msg('{} is not a valid smooth-factor'.format(smooth_factor))
        sys.stderr.write(_usage)
        sys.exit(1)

    demfun.blur(elev, '{}_smooth{}.tif'.format(utils.fn_basename2(elev), smooth_factor), smooth_factor)
### End
