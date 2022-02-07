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
    print('Usage: gdal_mask.py [-eq value] [-gt value] [-lt value] [-nodata]')
    print('                    [-overwrite] [-binmask] src_grid dest_grid')
    print('')
    print('gdal_mask v.%s' %(gfr_version))

# Mainline
if __name__ == "__main__":

    ingrd = None
    outgrd = None
    gt_data = None
    lt_data = None
    eq_data = None
    overwrite = False
    verbose = False
    mk_ndata = False
    ndata = None
    bin_mask = False
    applygrd = None
    
    i = 1
    while i < len(sys.argv):
        arg = sys.argv[i]
        if arg == '-gt':
            gt_data = utils.float_or(sys.argv[i+1])
            i = i + 1
        elif arg == '-lt':
            lt_data = utils.float_or(sys.argv[i+1])
            i = i + 1
        elif arg == '-eq':
            eq_data = utils.float_or(sys.argv[i+1])
            i = i + 1
        elif arg == '-nodata':
            mk_ndata = True
            ndata = sys.argv[i+1]
            i = i + 1
        elif arg == '-apply':
            bin_mask = True
            applygrd = sys.argv[i+1]
            i = i + 1
        elif arg == '-overwrite':
            overwrite = True
        elif arg == '-binmask':
            bin_mask = True
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

    if bin_mask:
        out_data = 0
    else:
        out_data = ds_config['ndv']
    
    if eq_data is not None:
        outarray[outarray==eq_data]=out_data
    if gt_data is not None:
        outarray[outarray>gt_data]=out_data
    if lt_data is not None:
        outarray[outarray<lt_data]=out_data
        
    if bin_mask:
        outarray[outarray != out_data] = 1
        
    if applygrd is not None:

        a_ds = gdal.Open(applygrd)
        a_ds_config = demfun.gather_infos(a_ds)
    
        a_band = a_ds.GetRasterBand(1)
        msk_array = a_band.ReadAsArray()
        dst_ds = demfun.copy_infos(a_ds_config)
        
        msk_array[outarray == 0] = a_ds_config['ndv']
        outarray = msk_array
        a_ds = None
    else:
        
        dst_ds = demfun.copy_infos(ds_config)
        if mk_ndata is True:
            dst_ds['ndv'] = float(ndata)
        
    utils.gdal_write(outarray, outgrd, dst_ds)
        
    ds = None
#--END
