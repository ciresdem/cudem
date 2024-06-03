#!/usr/bin/env python
### gdal_null.py
##
## Copyright (c) 2018 - 2024 CIRES DEM Team
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
## Create a nodata grid
##
### Code:

import sys
import os
from osgeo import gdal
import numpy as np
from cudem import utils
from cudem import regions
from cudem import gdalfun

_version = '0.1.11'
_usage = '''gdal_null.py ({}): generate a null grid
usage: gdal_null.py [-region xmin xmax ymin ymax] [-cell_size value]
                    [-t_nodata value] [-d_format grid-format] [-overwrite]
                    [-copy grid] [-verbose] output_grid

CIRES DEM home page: <http://ciresgroups.colorado.edu/coastalDEM>'''.format(_version)

def verbosePrint(xcount, ycount, extent, cellsize, outf):
    print('--')
    print('xcols: {} ycols: {}'.format(xcount, ycount))
    print('xmax: {} xmin: {}'.format(extent[1], extent[0]))
    print('ymax: {} ymin: {}'.format(extent[3], extent[2]))
    print('cellsize: {}'.format(cellsize))
    print('output grid format: {}'.format(outf))
    print('--')

def createNullCopy(srcfile, outfile, nodata, outformat, verbose, overwrite):
    '''copy a gdal grid and make a nodata grid'''
    ds = gdal.Open(srcfile)
    ds_config = gdalfun.gdal_infos(ds)
    gt = ds_config['geoT']

    if nodata is None: nodata = ds_config['ndv']
    if nodata is None: nodata = -9999
    ds_config['ndv'] = nodata

    ds = None
    dsArray = np.zeros([ds_config['ny'],ds_config['nx']])
    dsArray[:] = float(nodata)
    gdalfun.gdal_write(dsArray, outfile, ds_config)
    
def createGrid(outfile, extent, cellsize, nodata, outformat, verbose, overwrite):
    '''create a nodata grid'''
    xcount, ycount, gt = extent.geo_transform(x_inc = cellsize)
    ds_config = gdalfun.gdal_set_infos(xcount, ycount, xcount * ycount, gt, gdalfun.osr_wkt(4326), gdal.GDT_Float32, nodata, outformat, {}, 1)
    nullArray = np.zeros( (ycount, xcount) )
    nullArray[nullArray==0]=nodata
    gdalfun.gdal_write(nullArray, outfile, ds_config)

if __name__ == '__main__':

    extent = None
    cellsize = 1
    d_format = "GTiff"
    cpgrd = None
    overwrite = False
    verbose = False
    nodata = -9999
    output = None
    i = 1
    
    while i < len(sys.argv):
        arg = sys.argv[i]
        if arg == '-region' or arg == '-r' or arg == '--region':
            extent = (float(sys.argv[i+1]),float(sys.argv[i+2]),
                      float(sys.argv[i+3]),float(sys.argv[i+4]))
            i = i + 4
        elif arg == '-cell_size' or arg == '-s' or arg == '--cell_size':
            cellsize = float(sys.argv[i+1])
            i = i + 1
        elif arg == '-d_format' or arg == 'd' or arg == '--d_format':
            d_format = str(sys.argv[i+1])
            i = i + 1
        elif arg == '-t_nodata' or arg == '-t' or arg == '--t_nodata':
            nodata = float(sys.argv[i+1])
            i = i + 1
        elif arg == '-copy' or arg == '-c' or arg == '--copy':
            cpgrd = sys.argv[i+1]
            i = i + 1
        elif arg == '-overwrite' or arg == '--overwrite':
            overwrite = True
        elif arg == '-verbose' or arg == '--verbose':
            verbose = True
        elif arg[0] == '-':
            sys.stderr.write(_usage)
            sys.exit(1)
        elif output is None:
            output = arg
        else:
            sys.stderr.write(_usage)
            sys.exit(1)

        i = i + 1

    if output is None:
        sys.stderr.write(_usage)
        utils.echo_error_msg('you must enter an output file name')
        sys.exit(0)
        
    if extent == None:
        extent = '1'
    else:
        this_region = regions.Region().from_list(extent)

    #Run the program given the user input
    if cpgrd is not None:
        createNullCopy(cpgrd, output, nodata, d_format, verbose, overwrite)
    else:
        createGrid(output, this_region, cellsize, nodata, d_format, verbose, overwrite)

### END
