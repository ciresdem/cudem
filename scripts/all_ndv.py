#!/usr/bin/env python
### all_ndv.py
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
## Check if a grid file contains all nodata values
##
### Code:

import sys
from osgeo import gdal
import numpy as np

if __name__ == "__main__":

    if len(sys.argv) < 2:
        Usage()
        exit(1)
    else:
        ingrd = sys.argv[1]

    ds = gdal.Open(ingrd)
    band = ds.GetRasterBand(1)
    ndata = band.GetNoDataValue()

    band_data = band.ReadAsArray()
    ds = None
    if np.isnan(ndata):
        if np.all(np.isnan(band_data)):
            exit(0)
        else:
            exit(1)
    else:
        if np.all(band_data == ndata):
            exit(0)
        else:
            exit(1)
    exit(1)
   
