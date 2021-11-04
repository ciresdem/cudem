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
import numpy as np
from scipy.signal import fftconvolve
from scipy.signal import convolve
from gdalconst import *
from osgeo import osr
from osgeo import gdal

from geomods import gdalfun

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

def open_file_list(in_list, smooth_factor):
    il = open(in_list, 'r')
    for line in il:
        if line[0] != "#":
            proc_elev(line.strip(), smooth_factor)
    il.close()

def yield_file_list(in_list):
    with open(in_list, 'r') as iob:
        for line in iob:
            if line[0] != '#' and line[0] != '\n':
                yield(line)
    
def gaussian_blur(in_array, size):
    '''blur an array'''
    
    # expand in_array to fit edge of kernel
    padded_array = np.pad(in_array, size, 'symmetric')
    # build kernel
    x, y = np.mgrid[-size:size + 1, -size:size + 1]
    g = np.exp(-(x**2 / float(size) + y**2 / float(size)))
    g = (g / g.sum()).astype(in_array.dtype)
    in_array = None
    # do the Gaussian blur
    try:
        out_array = fftconvolve(padded_array, g, mode='valid')
        # work-around for when fftconvolve returns all nan values for some reason...
        #if np.nan in out_array[0][0]: 1+"A"
    except:
        print('switching to convolve')
        out_array = convolve(padded_array, g, mode='valid')
    return out_array

# Function to read the original file's projection:
def GetGeoInfo(FileName):
    '''get some info from the input gdal file'''
    
    SourceDS = gdal.Open(FileName, GA_ReadOnly)
    NDV = SourceDS.GetRasterBand(1).GetNoDataValue()
    xsize = SourceDS.RasterXSize
    ysize = SourceDS.RasterYSize
    GeoT = SourceDS.GetGeoTransform()
    Projection = osr.SpatialReference()
    Projection.ImportFromWkt(SourceDS.GetProjectionRef())
    DataType = SourceDS.GetRasterBand(1).DataType
    #DataType = gdal.GetDataTypeName(DataType)
    return xsize, ysize, GeoT, Projection, DataType, NDV

# Function to write a new file.
def CreateGeoTiff(Name, Array, driver,
                  xsize, ysize, GeoT, Projection, DataType, NDV):
    '''create a geotiff'''
    
    # if DataType == 'Float32':
    #     DataType = gdal.GDT_Float32
    NewFileName = Name+'.tif'
    # Set nans to the original No Data Value
    #Array[np.isnan(Array)] = NDV
    # Set up the dataset
    #print DataType
    DataSet = driver.Create( NewFileName, xsize, ysize, 1, DataType )
    # the '1' is for band 1.
    DataSet.SetGeoTransform(GeoT)
    
    wkt_proj = Projection.ExportToWkt()
    if wkt_proj.startswith("LOCAL_CS"):
        wkt_proj = wkt_proj[len("LOCAL_CS"):]
        wkt_proj = "PROJCS"+wkt_proj
    DataSet.SetProjection(wkt_proj)
    #DataSet.SetProjection( Projection.ExportToWkt() )
    
    # Write the array
    DataSet.GetRasterBand(1).WriteArray( Array )
    DataSet.GetRasterBand(1).SetNoDataValue(NDV)
    return NewFileName

def proc_elev(elev, smooth_factor):
    '''process the elev array'''
    
    if not os.path.exists(elev):
        print("Error: %s is not a valid file" %(elev))
    else:
        #Create Array
        output_name=elev[:-4]+"_smooth_"+str(smooth_factor)
        xsize, ysize, GeoT, Projection, DataType, NDV = GetGeoInfo(elev)
        
        print("elev is", elev)
        print("size is", xsize, ysize)
        print("nodata is", NDV)
        print("datatype is", DataType)
        print("smooth factor is", smooth_factor)
        print("output_name is", output_name)
        
        elev_g = gdal.Open(elev) #
        elev_array = elev_g.GetRasterBand(1).ReadAsArray(0,0,xsize,ysize) 
        mask_array = elev_array
        elev_array = None
        #Set topo values to zero
        mask_array[mask_array > 0] = 0
        mask_array[mask_array == NDV] = 0
        print("loaded input dem")

        #print mask_array
        #Perform smoothing
        smooth_elev=gaussian_blur(mask_array, smooth_factor)
        #print smooth_elev
        if np.isnan(smooth_elev[0][0]):
            output_name=elev[:-4]+"_smooth_fail"
        #print smooth_elev
        mask_array[mask_array < 0] = 1
        smooth_elev = smooth_elev * mask_array
        mask_array = None
        print("smoothed array")
    
        #Reload original array and merge the topo with the smoothed bathy
        elev_array = elev_g.GetRasterBand(1).ReadAsArray(0,0,xsize,ysize)
        elev_array[elev_array < 0] = 0
        smoothed_array = smooth_elev + elev_array
        elev_g = elev_array = smooth_elev = None
        
        #Export Tif
        driver = gdal.GetDriverByName('GTiff')
        CreateGeoTiff(output_name, smoothed_array, driver, xsize, ysize, GeoT, Projection, DataType, NDV)
        smoothed_array = None
        print("created Smoothed Geotiff")

#def smooth_bathy(src_gdal, smooth_factor):
    # gdal-split by 0
    # smooth lower
    # merge back with upper
        
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

    if in_list:
        for lf in yield_file_list(in_list):
            proc_elev(lf, smooth_factor)
    else:
        #proc_elev(elev, smooth_factor)
        gdalfun.gdal_blur(elev, elev[:-4] + '_smooth.tif', smooth_factor)

### End
