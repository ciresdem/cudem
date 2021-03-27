#!/usr/bin/env python
### gdal_histogram.py
##
## Copyright (c) 2010, 2021 CIRES Coastal DEM Team
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
## Return a histogram of the grid data
##
### Code:

import sys
from osgeo import gdal
from cudem import utils
import numpy as np

_version = 0.1

def Usage():
    print('Usage: gdal_minmax.py grdfile')
    print('')
    print('Will return xmin xmax ymin ymax zmin zmax')
    print('')
    print('gdal_minmax v.%s' %(_version))

def returnMinMax(ingrd):
    # Process the grid file

    try:
        ds = gdal.Open(ingrd)
    except:
        ds = None

    if ds is not None:
        comp_geot = ds.GetGeoTransform()
        ds_band = ds.GetRasterBand(1)
        ndata = ds_band.GetNoDataValue( )
        cellsize = [float(comp_geot[1]), float(comp_geot[5])]
        xmin = float(comp_geot[3])
        ymax = float(comp_geot[0])
        band = ds.GetRasterBand(1)
        tgrid = band.ReadAsArray()
        tgrid2 = tgrid[tgrid != ndata]

        hbin,hist = np.histogram(tgrid2, 1000)

        for i,j in enumerate(hbin):
            tmp_val = str(hist[i]) + " " + str(j)# + "\n"
            print(tmp_val)
    else:
        utils.echo_error_msg('could not access input file {}'.format(ingrd))    

if __name__ == "__main__":

    if len(sys.argv) <= 1:
        Usage()
        sys.exit()
    else:
        infile = sys.argv[1]

    returnMinMax(infile)

### END
