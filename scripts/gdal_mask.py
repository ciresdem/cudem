#!/usr/bin/env python
#
# Description: Mask data in a gdal grid
#
#--

import sys
import os
import struct
import numpy as np
from osgeo import gdal
from cudem import utils
from cudem import demfun

gfr_version = 0.1

def Usage():
    print('Usage: gdal_mask.py [-gt value] [-lt value] [-nodata]')
    print('                    [-overwrite] src_grid dest_grid')
    print('')
    print('gdal_mask v.%s' %(gfr_version))

# Mainline
if __name__ == "__main__":

    ingrd = None
    outgrd = None
    gt_data = None
    lt_data = None
    overwrite = False
    verbose = False
    mk_ndata = False
    ndata = None
    
    i = 1
    while i < len(sys.argv):
        arg = sys.argv[i]
        if arg == '-gt':
            gt_data = utils.float_or(sys.argv[i+1])
            i = i + 1
        elif arg == '-lt':
            lt_data = utils.float_or(sys.argv[i+1])
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
    if in_ndata is None:
        in_ndata = -9999
        
    comp_geot = ds_config['geoT']
    outarray = ds.ReadAsArray()

    if gt_data is not None:
        outarray[outarray>gt_data]=ds_config['ndv']
    if lt_data is not None:
        outarray[outarray<lt_data]=ds_config['ndv']

    # Create the output GDAL Raster
    #if np.any(t):
    #outarray[t]=float(rdata)
    dst_ds = demfun.copy_infos(ds_config)
    if mk_ndata is True:
        dst_ds['ndv'] = float(ndata)
        
    utils.gdal_write(outarray, outgrd, dst_ds)
        
    ds = None
#--END
