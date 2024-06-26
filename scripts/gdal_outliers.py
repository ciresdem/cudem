#!/usr/bin/env python
### gdal_outliers.py
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
## filter outliers in a gdal compatible DEM
##
### Code:

import sys
import numpy as np
from cudem import gdalfun
from cudem import utils

def get_outliers(in_array, percentile=75):
    """get the outliers from in_array based on the percentile"""

    if percentile <= 50:
        percentile = 51

    if percentile >= 100:
        percentile = 99

    max_percentile = percentile
    min_percentile = 100 - percentile

    perc_max = np.nanpercentile(in_array, max_percentile)
    perc_min = np.nanpercentile(in_array, min_percentile)
    iqr_p = (perc_max - perc_min) * 1.5
    upper_limit = perc_max + iqr_p
    lower_limit = perc_min - iqr_p

    return(upper_limit, lower_limit)

_version = '0.2.1'
_usage = '''gdal_outliers.py ({}): filter z outliers in a DEM

usage: gdal_outliers.py [ file [dst_dem opts] ]

 Options:
  file\t\t\tThe input DEM file-name
  dst_dem\t\tThe output DEM file-name

  --size\t\tThe size in pixels of the moving filter window
  --step\t\tThe step size of the moving filter window
  --percentile\t\tfilter percentile to identify outliers
  --uncertainty\t\tThe associated uncertainty grid

  --elevation_weight\tWeight of elevation outliers
  --curvature_weight\tWeight of curvature outliers
  --roughness_weight\tWeight of roughness outliers
  --tpi_weight\t\tWeight of the tpi outliers
  --tri_weight\t\tWeight of the tri outliers
  --uncertainty_weight\tWeight of uncertainty outliers
  --interpolation\tInterpolation method to fill nodata cells (cubic, linear, nearest)

  --aggressive\t\tBe more aggressive in data removal

  --help\t\tPrint the usage text
  --version\t\tPrint the version information

 Examples:
 % gdal_outliers.py input.tif -size 10 -step 2

CIRES DEM home page: <http://ciresgroups.colorado.edu/coastalDEM>
'''.format(_version)

if __name__ == "__main__":

    elev = None
    dst_gdal = None
    chunk_step = None
    chunk_size = None
    percentile = 75
    replace = False
    unc_mask = None
    i = 1

    elevation_weight = 1
    curvature_weight = 1
    roughness_weight = 1
    tri_weight = 1
    tpi_weight = 1
    uncertainty_weight = 1
    interpolation = 'cubic'
    aggressive = False
    
    argv = sys.argv
    while i < len(sys.argv):
        arg = sys.argv[i]
        if arg == '-size' or arg == '--size' or arg == '-z':
            chunk_size = utils.int_or(argv[i + 1])
            i = i + 1
        elif arg == '-step' or arg == '--step' or arg == '-s':
            chunk_step = utils.int_or(argv[i + 1])
            i = i + 1
        elif arg == '-percentile' or arg == '--percentile' or arg == '-p':
            percentile = utils.int_or(argv[i + 1])
            i = i + 1
        elif arg == '-replace' or arg == '--replace' or arg == '-r':
            replace = True
        elif arg == '-elevation_weight' or arg == '--elevation_weight' or arg == '-ew':
            elevation_weight = utils.float_or(argv[i + 1], 1)
            i = i + 1
        elif arg == '-curvature_weight' or arg == '--curvature_weight' or arg == '-cw':
            curvature_weight = utils.float_or(argv[i + 1], 1)
            i = i + 1
        elif arg == '-uncertainty_weight' or arg == '--uncertainty_weight' or arg == '-uw':
            uncertainty_weight = utils.float_or(argv[i + 1], 1)
            i = i + 1
        elif arg == '-roughness_weight' or arg == '--roughness_weight' or arg == '-rw':
            roughness_weight = utils.float_or(argv[i + 1], 1)
            i = i + 1
        elif arg == '-tri_weight' or arg == '--tri_weight' or arg == '-trw':
            tri_weight = utils.float_or(argv[i + 1], 1)
            i = i + 1
        elif arg == '-tpi_weight' or arg == '--tpi_weight' or arg == '-tw':
            tpi_weight = utils.float_or(argv[i + 1], 1)
            i = i + 1
        elif arg == '-uncertainty' or arg == '--uncertainty' or arg == '-u':
            unc_mask = argv[i + 1]
            i = i + 1
        elif arg == '-interpolation' or arg == '--interpolation' or arg == '-i':
            interpolation = argv[i + 1]
            i = i + 1
        elif arg == '-aggressive' or arg == '--aggressive' or arg == '-a':
            aggressive = True
        elif arg == '-help' or arg == '--help' or arg == '-h':
            sys.stderr.write(_usage)
            sys.exit(1)
        elif arg == '-version' or arg == '--version':
            sys.stderr.write('{}\n'.format(_version))
            sys.exit(1)
        elif elev is None:
            elev = arg
        elif dst_gdal is None:
            dst_gdal = arg
        else:
            utils.echo_warning_msg(arg)
            sys.stderr.write(_usage)
            sys.exit(1)
        i = i + 1

    if elev is None:
        sys.stderr.write(_usage)
        utils.echo_error_msg('you must enter an input file')
        sys.exit(1)

    if dst_gdal is None:
        dst_gdal = elev.split('.')[0] + '_fltr.tif'
        
    # gdalfun.gdal_filter_outliers2(
    #     elev, dst_gdal, unc_mask=unc_mask, chunk_size=chunk_size, chunk_step=chunk_step, percentile=percentile,
    #     elevation_weight=elevation_weight, curvature_weight=curvature_weight, rough_weight=roughness_weight,
    #     unc_weight=uncertainty_weight, tpi_weight=tpi_weight, tri_weight=tri_weight, interpolation=interpolation
    # )
    elev_array = gdalfun.gdal_get_array(elev)
    print(get_outliers(elev_array[0]))

### End
