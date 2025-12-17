#!/usr/bin/env python
#
# Description: polygonize a cudem mask multiband raster for
# spatial metadata
#
#--

import sys
from osgeo import gdal
from cudem import gdalfun
from cudem import dlim
from cudem import regions
from cudem import utils

gfr_version = 0.3

def Usage():
    print('Usage: gdal_polygonize_mask.py src_mask [opts]')
    print('')
    print('gdal_polygonize_mask v.%s' %(gfr_version))

# Mainline
if __name__ == "__main__":

    src_mask = None
    verbose = False
    mask_level = 0
    want_footprint = False
    xy_inc = [None, None]
    region = None
    outfn = None
    
    i = 1
    while i < len(sys.argv):
        arg = sys.argv[i]
        if arg == '--increment' or arg == '-E':
            xy_inc = sys.argv[i + 1].split('/')
            i = i + 1
        elif arg[:2] == '-E':
            xy_inc = arg[2:].split('/')

        elif arg == '--region' or arg == '-R':
            region = str(sys.argv[i + 1])
            i = i + 1
        elif arg[:2] == '-R':
            region = str(arg[2:])

        elif arg == '--outfile' or arg == '-O':
            outfn = str(sys.argv[i + 1])
            i = i + 1
        elif arg[:2] == '-O':
            outfn = str(arg[2:])
            
        elif arg == '--mask_level' or arg == '-mask_level':
            mask_level = utils.int_or(sys.argv[i+1], 0)
            i += 1
        elif arg == '--want_footprint' or arg == '-want_footprint':
            want_footprint = True
        elif arg == '-verbose':
            verbose = True
        elif arg[0] == '-':
            Usage()
        elif src_mask is None:
            src_mask = arg
        else: Usage()
        i = i + 1

    if src_mask is None:
        Usage()
        sys.exit(0)

    if region is None:
        this_region = region
    else:
        this_region = regions.parse_cli_region([region], verbose)[0]
        utils.echo_msg_bold(this_region)
    
    if len(xy_inc) < 2:
        xy_inc.append(xy_inc[0])
    elif len(xy_inc) == 0:
        xy_inc = [None, None]

    if xy_inc[0] is not None and xy_inc[1] is not None:
        desc = src_mask
        src_mask = gdalfun.sample_warp(
            src_mask, None, utils.str2inc(xy_inc[0]), utils.str2inc(xy_inc[1]),
            src_region=this_region,
            sample_alg='average',
            dst_srs=None,
            ndv=gdalfun.gdal_get_ndv(src_mask),
            verbose=True,
            ot=gdal.GDT_Byte,
            co=["COMPRESS=DEFLATE", "TILED=YES"]
        )[0]
        src_mask.SetDescription(desc)
        
    has_gdal_footprint = utils.cmd_exists('gdal_footprint')
    with gdalfun.gdal_datasource(src_mask) as msk_ds:
        if has_gdal_footprint and want_footprint:
            sm_files, sm_fmt = dlim.ogr_mask_footprints(msk_ds, verbose=True, mask_level=mask_level)
        else:
            sm_layer, sm_fmt = dlim.polygonize_mask_multibands(msk_ds, verbose=True, mask_level=mask_level, output=outfn)

#--END
