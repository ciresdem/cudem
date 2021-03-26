#!/usr/bin/env python
#
# Description: Find and replace a value in a gdal-compatible grid file
#
#--

import sys
import os
import struct
import numpy as np
from osgeo import gdal
from cudem import utils
from cudem import demfun

gfr_version = 0.3
def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

def Usage():
    print('Usage: gdal_findreplace.py [-s_value value] [-t_value value] [-nodata]')
    print('                           [-row value] [-column value]')
    print('                           [-overwrite] src_grid dest_grid')
    print('')
    print('gdal_findreplace v.%s' %(gfr_version))

# Mainline
if __name__ == "__main__":

    ingrd = None
    outgrd = None
    fdata = None
    rdata = None
    overwrite = False
    verbose = False
    mk_ndata = False
    ndata = None
    row = None
    column = None
    
    i = 1
    while i < len(sys.argv):
        arg = sys.argv[i]
        if arg == '-s_value':
            fdata = sys.argv[i+1]
            i = i + 1
        elif arg == '-t_value':
            rdata = sys.argv[i+1]
            i = i + 1
        elif arg == '-row':
            row = argv[i+1]
            i = i + 1
        elif arg == '-column':
            column = sys.argv[i+1]
            i = i + 1
        elif arg == '-nodata':
            mk_ndata = True
            ndata = sys.argv[i+1]
            i = i + 1
        elif arg == '-overwrite':
            overwrite = True
        elif arg == '-verbose':
            verbose = True
        elif arg[0] == '-':
            Usage()
        elif ingrd is None:
            ingrd = arg
        elif outgrd is None:
            outgrd = arg
        else: Usage()
        i = i + 1

    if ingrd is None or outgrd is None:
        Usage()
        sys.exit(0)

    ds = gdal.Open(ingrd)
    ds_config = demfun.gather_infos(ds)
    
    band = ds.GetRasterBand(1)
    in_ndata = band.GetNoDataValue()
    if in_ndata is None: in_ndata = -9999
    comp_geot = ds_config['geoT']
    outarray = ds.ReadAsArray()

    if fdata == 'all':
        outarray[outarray!=in_ndata]=rdata
    elif np.isnan(float(fdata)):
        outarray[np.isnan(outarray)]=rdata
    else:
        t = np.isclose(outarray, float(fdata))
        if verbose: sys.stderr.write('{}.'.format(np.any(t)))
    
    # Create the output GDAL Raster
    if np.any(t):
        outarray[t]=float(rdata)
        dst_ds = demfuncopy_infos(ds_config)
        if mk_ndata is True: dst_ds['ndv'] = float(ndata)
        utils.gdal_write(outarray, outgrd, dst_ds)
        
    ds = None
#--END
