#!/usr/bin/env python
#
# Description: Mask data in a gdal grid
#
#--

import sys
from cudem import gdalfun
from cudem import utils
from cudem import waffles

_version = '0.0.1'

_usage = '''gdal_hydro_flatten.py ({}): Flatten nodata voids

usage: gdal_hydro_flatten.py [ src_dem [dst_dem opts] ]

 Options:
  src_dem\t\tThe input DEM file-name
  dst_dem\t\tThe output DEM file-name

  --band\t\tThe input raster band (1)
  --size_threshold\tThe minimum size (in cells) to flatten
  --help\t\tPrint the usage text
  --version\t\tPrint the version information

 Examples:
 % gdal_hydro_flatten.py input_dem.tif output_dem.tif --size_threshold 10

CIRES DEM home page: <http://ciresgroups.colorado.edu/coastalDEM>
'''.format(_version)
    
# Mainline
if __name__ == "__main__":

    src_dem = None
    dst_dem = None
    band = 1
    size_threshold = 1
    verbose = False
    
    i = 1
    while i < len(sys.argv):
        arg = sys.argv[i]
        if arg == '--size_threshold':
            size_threshold = utils.int_or(sys.argv[i+1], 1)
            i += 1
        elif arg == '--band':
            band = utils.int_or(sys.argv[i+1], 1)
            i += 1
        elif arg == '-verbose':
            verbose = True
        elif src_dem is None:
            src_dem = arg
        elif dst_dem is None:
            dst_dem = arg
        else:
            sys.stderr.write(_usage)
            sys.exit(1)
        i = i + 1

    if src_dem is None:
        sys.stderr.write(_usage)
        sys.exit(1)

    #gdalfun.gdal_hydro_flatten(src_dem, dst_dem=dst_dem, band=band, size_threshold=size_threshold)
    waffles.flatten_no_data_zones(src_dem, dst_dem=dst_dem, band=band, size_threshold=size_threshold)

#--END
