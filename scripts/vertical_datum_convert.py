#!/usr/bin/env python
### vertical_datum_convert.py
##
## Copyright (c) 2022, 2023 Regents of the University of Colorado
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

from osgeo import gdal

from cudem import utils
from cudem import regions
from cudem import demfun
from cudem import vdatums

_version = '0.1.1'
_epsg_desc = lambda t, x: '{}:\n '.format(t) + ' '.join(
    ['\033[1m{}\033[0m\t{}\n'.format(key, x[key]['name']) for key in x])

_usage = """{cmd} ({version}): transform a grid between vertical datums

usage: {cmd} [OPTIONS] input_grid output_grid

  input_grid\t\tThe input raster to transform
  output_grid\t\tThe output transformed raster

 Options:

  -i, --vdatum_in\tthe input vertical datum as EPSG code
  -o, --vdatum_out\tthe output vertical datum as EPSG code
  -D, --cache-dir\tCACHE Directory for storing temp data.
\t\tDefault Cache Directory is ~/.cudem_cache; cache will be cleared after a waffles session
\t\tto retain the data, use the --keep-cache flag

  -k, --keep-cache\tKEEP the cache data intact after run
  --list-epsg\t\tList the supported EPSG codes and their names
  --help\t\tPrint the usage text
  --version\t\tPrint the version information

 Examples:
 % {cmd} my_dem_navd88.tif my_dem_wgs84_1674.tif --vdatum_in 5703 --vdatum_out 7662

CIRES DEM home page: <http://ciresgroups.colorado.edu/coastalDEM>
""".format(version=_version, cmd=os.path.basename(sys.argv[0]))
        
def main():
    src_grid = None
    dst_grid = None
    vdatum_in = 5703
    vdatum_out = 7662
    verbose = False
    keep_cache = False
    cache_dir = utils.cudem_cache

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

        elif arg == '--cache-dir' or arg == '-D' or arg == '-cache-dir':
            cache_dir = os.path.join(utils.str_or(argv[i + 1], os.path.expanduser('~')), '.cudem_cache')
            i = i + 1
        elif arg[:2] == '-D': cache_dir = os.path.join(utils.str_or(argv[i + 1], os.path.expanduser('~')), '.cudem_cache')
        elif arg == '--list-epsg':
            #print(_epsg_desc(htdpfun.HTDP()._reference_frames))
            print(_epsg_desc('htdp epsg', vdatums._htdp_reference_frames))
            print(_epsg_desc('cdn espg', vdatums._cdn_reference_frames))
            print(_epsg_desc('tidal epsg', vdatums._tidal_frames))
            #list_proj_cdn_epsg()
            sys.exit(1)
        elif arg == '-k' or arg == '--keep-cache': keep_cache = True
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
            sys.stderr.write(_usage)
            sys.exit(1)
        i = i + 1

    if src_grid is None or not os.path.exists(src_grid):
        sys.stderr.write(_usage)
        sys.exit(1)

    if dst_grid is None:
        dst_grid = '.'.join(src_grid.split('.')[:-1]) + '_' + str(vdatum_out.replace('(', '_').replace(')', '_')) + '.' + src_grid.split('.')[-1]

    if not os.path.exists(src_grid):
        utils.echo_error_msg('Error: {} is not a valid file'.format(src_grid))
    else:
        src_infos = demfun.infos(src_grid)
        
        src_region = regions.Region().from_geo_transform(src_infos['geoT'], src_infos['nx'], src_infos['ny'])
        src_region.src_srs = demfun.get_srs(src_grid)

        trans_region = src_region.copy()
        trans_region.warp()
        trans_region.buffer(pct=2)

        trans_region._wgs_extremes()
        
        #x_inc, y_inc = trans_region.increments(src_infos['nx']/3, src_infos['ny']/3)
        
        x_inc = src_infos['geoT'][1]
        y_inc = -src_infos['geoT'][5]
        tmp_x_inc = 3/3600
        tmp_y_inc = 3/3600
        
        vt = vdatums.VerticalTransform(trans_region, tmp_x_inc, tmp_y_inc, vdatum_in, vdatum_out, cache_dir=cache_dir)
        _trans_grid = vt.run()

        if os.path.exists('_{}'.format(_trans_grid)):
            utils.remove_glob('_{}'.format(_trans_grid))
        
        if _trans_grid is not None:

            out_h, out_v = utils.epsg_from_input(demfun.get_srs(src_grid))
            
            utils.run_cmd('gdalwarp {} {} -te {} -ts {} {} -s_srs epsg:4326 -t_srs epsg:{}'.format(
                _trans_grid, '_{}'.format(_trans_grid),
                src_region.format('te'),
                src_infos['nx'],
                src_infos['ny'],
                out_h), verbose=True)
            
            # utils.run_cmd(
            #     'gdalwarp {} {} -te {} -tr {} {} -s_srs epsg:4326 -t_srs {} -co COMPRESS=LZW -co TILED=YES -co PREDICTOR=3'.format(
            #         _trans_grid,
            #         '_{}'.format(_trans_grid),
            #         src_region.format('te'),
            #         x_inc, y_inc,
            #         demfun.get_srs(src_grid)
            #     ), verbose=True
            # )

            # out, status = utils.run_cmd(
            #     'gdal_calc.py -A {} -B {} --calc "A+B" --outfile {} --co COMPRESS=LZW --co TILED=YES --co PREDICTOR=3 --overwrite'.format(
            #         src_grid.replace(' ', '\ '), '_{}'.format(_trans_grid).replace(' ', '\ '), dst_grid.replace(' ', '\ ')
            #     ),
            #     verbose=True
            # )
            # if status == 0:

            gdc_cmd = 'gdal_calc.py -A {} -B {} --calc "A+B" --outfile {} --co COMPRESS=LZW --co TILED=YES --co PREDICTOR=3 --overwrite'.format(
                src_grid.replace(' ', '\ '), '_{}'.format(_trans_grid).replace(' ', '\ '), dst_grid.replace(' ', '\ '))
            os.system(gdc_cmd)
                
        else:
            utils.echo_error_msg('could not parse input/output vertical datums: {} -> {}; check spelling, etc'.format(vdatum_in, vdatum_out))

        if not keep_cache:
            utils.remove_glob(cache_dir)
            
if __name__ == '__main__':
    main()

### End
