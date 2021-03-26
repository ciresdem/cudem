#!/usr/bin/env python
#
# Description: Convert nan values in a gdal-compatible grid
#
#--

#--
import sys
import os
import struct
import numpy as np
from osgeo import gdal
#--

gnn_version = 0.2

#--
def Usage():
    print('Usage: gdal_nan2null.py ingrd out_nan outgrd')
    print('')
    print('This script will convert nan values in a gdal compatible grid file with')
    print('the user-specified value.')
    print('')
    print('gdal_nan2null v.%s' %(gnn_version))
#--

#--
#
# Mainline
#
#--
if __name__ == "__main__":

    if len(sys.argv) < 2:
        Usage()
        sys.exit()
    else:
        ingrd = sys.argv[1]
        ndata = float(sys.argv[2])
        outgrd = sys.argv[3]

    overwrite = True
    verbose = True
    print >> sys.stderr, "gdal_nan2null: nan--->%s" %(ndata)
    ds = gdal.Open(ingrd)
    band = ds.GetRasterBand(1)
    comp_geot = ds.GetGeoTransform()
    outarray = ds.ReadAsArray()
    outarray[np.isnan(outarray)]=ndata
    
    # Create the output GDAL Raster

    outfile = outgrd
    (ycount,xcount) = outarray.shape[0],outarray.shape[1]
    
    driver = gdal.GetDriverByName("GTiff")

    if os.path.exists(outfile):
        if not overwrite:
            sys.exit("Error - %s already exists!  Use -O or --overwrite to force an \
            overwrite of the output file." %(outfile))
        else:
            driver.Delete(outfile)
            if verbose:
                print >> sys.stderr, "Warning - Overwriting file %s " %(outfile)
                
    dst_ds = driver.Create(outfile, xcount, ycount, 1, gdal.GDT_Float32)

    if dst_ds is None:
        sys.exit("failed to open output file...%s" %(outfile))
            
    gt = comp_geot #(xextent,comp_geot[1],0,yextent,0,comp_geot[5])
    dst_ds.SetGeoTransform(gt)
    dst_band = dst_ds.GetRasterBand(1)
        
    dst_band.SetNoDataValue( ndata )
        
    # Write Numpy Array to the GDAL Raster Band
    dst_band.WriteArray(outarray)
        
    dst_ds = None
    ds = None
#--END
