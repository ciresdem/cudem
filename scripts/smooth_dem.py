#!/usr/bin/env python
### smooth_dem_bathy.py
##
## Copyright (c) 2012 - 2021 CIRES Coastal DEM Team
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
#from geomods import gdalfun
from cudem import demfun

_version = '0.0.7'
_usage = '''smooth_dem_bathy.py ({}): smooth the bathymetry in a DEM
smooth_dem_bathy.py: A script that smooths the bathy areas of DEM (below 0) and merges back with original, unsmoothed topo.

usage: smooth_dem_bathy.py [ si [ args ] ] [ file ]

Options:
  file\t\tThe input DEM file-name
  -s\t\tValue for the smooth factor; default is 10
  -i\t\tA file containing DEM file-names to process

  -help\t\tPrint the usage text
  -version\tPrint the version information

Example:
smooth_dem_bathy.py input.tif -s 12
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
            print(_usage)
            sys.exit(1)
        elif arg == '-version' or arg == '--version':
            print('smooth_dem_bathy.py v.%s' %(_version))
            print(_license)
            sys.exit(1)
        elif elev is None:
            elev = arg
        else:
            print(_usage)
            sys.exit(1)
        i = i + 1

    if elev is None and in_list is None:
        print(_usage)
        sys.exit(1)

    try: smooth_factor = int(smooth_factor)
    except:
        print("Error: %s is not a valid smooth-factor" %(smooth_factor))
        print(_usage)
        sys.exit(1)

    demfun.blur(elev, elev[:-4] + '_smooth{}.tif'.format(smooth_factor), smooth_factor)
