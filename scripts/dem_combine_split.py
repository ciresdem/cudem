#!/usr/bin/env python
### dem_combine_split.py
##
## Copyright (c) 2018 - 2021 CIRES Coastal DEM Team
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
### Code:

import os
import sys
from osgeo import gdal
from cudem import utils
from cudem import demfun
from cudem import regions

_version = '0.0.1'
_usage = '''dem_combine_split.py ({}):

usage: gdal_combine_split.py [ dems ]

 Options:
  dems\t\tThe input DEM file-name

  --tr\t\tOutput resolution
  --help\tPrint the usage text
  --version\tPrint the version information

 Examples:
 % gdal_combine_split *.tif

CIRES DEM home page: <http://ciresgroups.colorado.edu/coastalDEM>
'''.format(_version)

if __name__ == '__main__':    
    dems = []
    i_region = None
    tr = None
    i = 1
    while i < len(sys.argv):
        arg = sys.argv[i]
        if arg == '--tr' or arg == '-E':
            incs = sys.argv[i + 1].split('/')
            xinc = utils.str2inc(incs[0])
            if len(incs) > 1:
                yinc = utils.str2inc(incs[1])
            else:
                yinc = utils.str2inc(incs[0])
            tr = (xinc, yinc)
            i = i + 1

        elif arg == '-help' or arg == '--help' or arg == '-h':
            sys.stderr.write(_usage)
            sys.exit(1)
        elif arg == '-version' or arg == '--version':
            sys.stderr.write('{}\n'.format(_version))
            sys.exit(1)
        else: dems.append(arg)

        i = i + 1

    if len(dems) < 1:
        sys.stderr.write(_usage)
        utils.echo_error_msg('you must enter at least one input file')
        sys.exit(1)
    else:
        tmp_combined = '__tmp_comb__.vrt'
        # g = gdal.Warp(tmp_combined, dems, format='GTiff',
        #               options=["COMPRESS=LZW", "TILED=YES"], callback=gdal.TermProgress)
        # g = None

        if tr is None:
            utils.run_cmd('gdalbuildvrt -allow_projection_difference -r bilinear {} {}'.format(tmp_combined, ' '.join(dems)), verbose=True)
        else:
            utils.run_cmd('gdalbuildvrt -allow_projection_difference -r bilinear -tr {} {} {} {}'.format(*tr, tmp_combined, ' '.join(dems)), verbose=True)
            
        c_ds = gdal.Open(tmp_combined)
        c_gt = c_ds.GetGeoTransform()
        c_xcount = c_ds.RasterXSize
        c_ycount = c_ds.RasterYSize
        c_region = regions.Region().from_geo_transform(
            geo_transform=c_gt,
            x_count=c_xcount,
            y_count=c_ycount
        )

        c_ds = None

        for dem in dems:

            src_ds = gdal.Open(dem)            
            if src_ds is not None:
                gt = src_ds.GetGeoTransform()

                #print(gt[1], gt[5]*-1)
                #print(tr[0], tr[1])
                #if gt[1] == tr[0] and gt[5]*-1 == tr[1]:

                if '{:g}'.format(gt[1]) == '{:g}'.format(tr[0]) and \
                   '{:g}'.format(gt[5]*-1) == '{:g}'.format(tr[1]):
                    this_region = regions.Region().from_geo_transform(
                        geo_transform=gt,
                        x_count=src_ds.RasterXSize,
                        y_count=src_ds.RasterYSize
                    )

                    output_name = dem[:-4] + '_p.tif'

                    this_srcwin = this_region.srcwin(
                        geo_transform=c_gt,
                        x_count=c_xcount,
                        y_count=c_ycount
                    )

                    src_ds = None

                    utils.run_cmd('gdal_translate {} {} -srcwin {} {} {} {}'.format(
                        tmp_combined, output_name, this_srcwin[0], this_srcwin[1], this_srcwin[2], this_srcwin[3]
                    ), verbose=True)
                    #demfun.cut(tmp_combined, this_region, output_name)
                
        utils.remove_glob(tmp_combined)       
### End
