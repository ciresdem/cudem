#!/usr/bin/env python
### gdal_minmax.py
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
## Return the min/max values of a gdal compatible DEM
##
### Code:

import sys
from osgeo import gdal
import numpy as np

grc_version = 0.2

def Usage():
    print('Usage: gdal_minmax.py grdfile')
    print('')
    print('Will return xmin xmax ymin ymax zmin zmax')
    print('')
    print('gdal_minmax v.%s' %(grc_version))

def returnMinMax(ingrd):
    # Process the grid file
    ds = gdal.Open(ingrd)
    comp_geot = ds.GetGeoTransform()

    cellsize = [float(comp_geot[1]), float(comp_geot[5])]
    xmin = float(comp_geot[0])
    ymax = float(comp_geot[3])

    band = ds.GetRasterBand(1)

    nodata = band.GetNoDataValue()
    tgrid = band.ReadAsArray()

    #tgrid[tgrid==nodata]=float('nan')
    #tgrid = band.ReadAsArray()

    #idx=(mask==nodata)
    #tgrid[idx]==
    tgrid=np.ma.MaskedArray(tgrid, mask=(tgrid==nodata))

    #print comp_geot
    xs = ds.RasterXSize
    ys = ds.RasterYSize
    zmin = tgrid.min()
    zmax = tgrid.max()

    xmax = xmin + (cellsize[1] * xs)
    ymin = ymax + (cellsize[0] * ys)
    print(xmin, xmax, ymin, ymax, zmin, zmax, xs, ys)

if __name__ == "__main__":

    if len(sys.argv) <= 1:
        Usage()
        sys.exit()
    else:
        infile = sys.argv[1]

    returnMinMax(infile)

# END
