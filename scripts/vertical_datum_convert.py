#!/usr/bin/env python
### vertical_datum_convert.py
##
## Copyright (c) 2022 Regents of the University of Colorado
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
## tranform a grid between vertical datums
##
## uses htdp and proj cdn
##
### Code:

import os
import sys

from cudem import utils
from cudem import regions
from cudem import demfun
from cudem import vdatums

_version = '0.1.0'
_epsg_desc = lambda t, x: '{}:\n '.format(t) + ' '.join(
    ['\033[1m{}\033[0m\t{}\n'.format(key, x[key]['name']) for key in x])

_usage = """{cmd} ({version}): transform a grid between vertical datums

usage: {cmd} [OPTIONS] input_grid output_grid

  input_grid\t\tThe input raster to transform
  output_grid\t\tThe output transformed raster

 Options:

  -i, --vdatum_in\tthe input vertical datum as EPSG code
  -o, --vdatum_out\tthe output vertical datum as EPSG code

  --list-epsg\t\tList the supported EPSG codes and their names
  --help\t\tPrint the usage text
  --version\t\tPrint the version information

 Examples:
 % {cmd} my_dem_navd88.tif my_dem_wgs84_1674.tif --vdatum_in 5703 --vdatum_out 7662

CIRES DEM home page: <http://ciresgroups.colorado.edu/coastalDEM>
""".format(version=_version, cmd=os.path.basename(sys.argv[0]))
#{htdp_epsg}
#, htdp_epsg=_htdp_epsg_desc(htdpfun.HTDP()._reference_frames))
        
def main():
    src_grid = None
    dst_grid = None
    vdatum_in = 5703
    vdatum_out = 7662
    verbose = False

    i = 1

    argv = sys.argv
    while i < len(argv):
        arg = argv[i]

        if arg == '-i' or arg == '--vdatum_in':
            vdatum_in = argv[i + 1]
            i = i + 1
        elif arg == '-o' or arg == '--vdatum_out':
            vdatum_out = argv[i + 1]
            i = i + 1
        elif arg == '--list-epsg':
            #print(_epsg_desc(htdpfun.HTDP()._reference_frames))
            print(_epsg_desc('htdp epsg', vdatums._htdp_reference_frames))
            print(_epsg_desc('cdn espg', vdatums._cdn_reference_frames))
            print(_epsg_desc('tidal epsg', vdatums._tidal_frames))
            #list_proj_cdn_epsg()
            sys.exit(1)
        elif arg == '--verbose': verbose = True
        elif arg == '-help' or arg == '--help' or arg == '-h':
            print(_usage)
            sys.exit(1)
        elif arg == '-version' or arg == '--version':
            print('vertical_datum_convert.py, version {}'.format(_version))
            sys.exit(1)
        elif src_grid is None:
            src_grid = arg
        elif dst_grid is None:
            dst_grid = arg
        else:
            print(_usage)
            sys.exit(1)
        i = i + 1

    if src_grid is None:
        print(_usage)
        sys.exit(1)

    if dst_grid is None:
        dst_grid = '.'.join(src_grid.split('.')[:-1]) + '_' + str(vdatum_out.replace('(', '_').replace(')', '_')) + '.' + src_grid.split('.')[-1]

    if not os.path.exists(src_grid):
        print('Error: {} is not a valid file'.format(src_grid))
    else:
        src_infos = demfun.infos(src_grid)
        src_region = regions.Region().from_geo_transform(src_infos['geoT'], src_infos['nx'], src_infos['ny'])
        x_inc, y_inc = src_region.increments(src_infos['nx'], src_infos['ny'])
        tmp_x_inc, tmp_y_inc = src_region.increments(src_infos['nx']/10, src_infos['ny']/10)
        vt = vdatums.VerticalTransform(src_region, tmp_x_inc, tmp_y_inc, vdatum_in, vdatum_out)
        _trans_grid = vt.run()
        
        if _trans_grid is not None:
            utils.run_cmd('gdalwarp {} {} -ts {} {}'.format(_trans_grid, '_{}'.format(_trans_grid), src_infos['nx'], src_infos['ny']), verbose=True)
            utils.run_cmd('gdal_calc.py -A {} -B {} --calc "A+B" --outfile {}'.format(src_grid, '_{}'.format(_trans_grid), dst_grid), verbose=True)
            utils.remove_glob(_trans_grid, '_{}'.format(_trans_grid))
        else:
            utils.echo_error_msg('could not parse input/output vertical datums: {} -> {}; check spelling, etc'.format(vdatum_in, vdatum_out))
        
if __name__ == '__main__':
    main()

### End
