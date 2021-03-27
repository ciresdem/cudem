#!/usr/bin/env python
### gdal_crop.py
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
## Check if a grid file contains any null values
##
### Code:

import sys
import struct
from osgeo import gdal
from cudem import utils

def Usage():
    print("""

+---------------+  has_nulls

Usage: has_nulls.py grd_datasource_name

This script will check a GDAL compatible grid file for null-values, 
and will print out the location (as y,x cell numbers) of the first
no-data value found, otherwise, will not print out anything.
""")

def isNaN(num):
    return(num != num)

if __name__ == "__main__":

    if len(sys.argv) < 2:
        Usage()
        sys.exit()
    else:
        ingrd = sys.argv[1]

    ds = gdal.Open(ingrd)
    band = ds.GetRasterBand(1)
    ndata = band.GetNoDataValue()

    if isNaN(ndata):
        band.SetNoDataValue(-9999)
        ndata = band.GetNoDataValue()

    for i in range(0, band.YSize):
        scanline = band.ReadRaster(0,i,band.XSize, 1, band.XSize, 1, gdal.GDT_Float32)
        these_values = struct.unpack('f' * band.XSize, scanline)
        tv2 = list(these_values)
        try: 
            tv2.index(ndata)
            utils.echo_msg("There is a NoData Value (%s) at y: %s x: %s." %(ndata, i, tv2.index(ndata)))
            break
        except:
            tv2=None
    ds = None
    
### End
