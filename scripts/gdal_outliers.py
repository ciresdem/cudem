#!/usr/bin/env python

import sys
from cudem import demfun
from cudem import utils

if __name__ == "__main__":
    src_gdal = sys.argv[1]
    dst_gdal = src_gdal.split('.')[0] + '_fltr.tif'
    chunk_size = 10
    chunk_step = chunk_size
    threshhold = 5
    utils.echo_msg('filtering {} to {}'.format(src_gdal, dst_gdal))
    demfun.filter_outliers(src_gdal, dst_gdal, threshhold=threshhold, chunk_size=chunk_size, chunk_step=chunk_step)

### End
