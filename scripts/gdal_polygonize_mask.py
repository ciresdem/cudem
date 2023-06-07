#!/usr/bin/env python
#
# Description: Mask data in a gdal grid
#
#--

import sys
from cudem import gdalfun

gfr_version = 0.1

def Usage():
    print('Usage: gdal_polygonize_mask.py src_mask')
    print('')
    print('gdal_polygonize_mask v.%s' %(gfr_version))

# Mainline
if __name__ == "__main__":

    ingrd = None
    verbose = False
    
    i = 1
    while i < len(sys.argv):
        arg = sys.argv[i]
        if arg == '-verbose':
            verbose = True
        elif arg[0] == '-':
            Usage()
        elif ingrd is None:
            ingrd = arg
        else: Usage()
        i = i + 1

    if ingrd is None:
        Usage()
        sys.exit(0)

    with gdalfun.gdal_datasource(ingrd) as msk_ds:
        sm_layer, sm_fmt = gdalfun.ogr_polygonize_multibands(msk_ds)

#--END
