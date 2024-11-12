#!/usr/bin/env python
### fetch_osm_coastline.py
##
## Copyright (c) 2024 CIRES Coastal DEM Team
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
from cudem import utils
from cudem import fetches
from cudem import gdalfun
from cudem import regions
from tqdm import tqdm

def fetch_coastline(this_region):
    if this_region is not None and this_region.valid_p():
        utils.echo_msg('fetching coastline for region {}'.format(this_region))
        cache_dir = utils.cudem_cache()
        this_cst = fetches.OpenStreetMap(src_region=this_region, verbose=True, outdir=cache_dir, q='coastline')#.run()
        fr = fetches.fetch_results(this_cst)
        fr.daemon=True
        fr.start()
        fr.join()
        return(fr)

    return(None)

_version = '0.0.2'
_usage = '''fetch_osm_coastline.py ({}): fetch and process the OSM coastline

usage: fetch_osm_coastline.py [-RBli [output-shape-file]]

 Options:
  -R, --region\t\tRestrict processing to the desired REGION 
\t\t\tWhere a REGION is xmin/xmax/ymin/ymax[/zmin/zmax[/wmin/wmax/umin/umax]]
\t\t\tUse '-' to indicate no bounding range; e.g. -R -/-/-/-/-10/10/1/-/-/-
\t\t\tOR an OGR-compatible vector file with regional polygons. 
\t\t\tWhere the REGION is /path/to/vector[:zmin/zmax[/wmin/wmax/umin/umax]].
\t\t\tIf a vector file is supplied, will use each region found therein.
  -B, --line_buffer\tBuffer the OSM coastline in degrees

  --include_landmask\tinclude the landmask in the output
  --invert_watermask\tinvert the watermask to the landmask

  --help\t\tPrint the usage text
  --version\t\tPrint the version information

 Examples:
 % fetch_osm_coastline.py -R -123.5/-123.25/38.5/38.75

'''.format(_version)

if __name__ == '__main__':    
    region = None
    include_landmask = False
    invert_watermask = False
    line_buffer = 0.0000001
    dst_ogr = None
    i = 1
    
    while i < len(sys.argv):
        arg = sys.argv[i]
        if arg == '--region' or arg == '-R':
            region = str(sys.argv[i + 1])
            i = i + 1
        elif arg[:2] == '-R':
            region = str(arg[2:])
        elif arg == '--line_buffer' or arg == '-B':
            line_buffer = utils.float_or(sys.argv[i + 1], line_buffer)
            i = i + 1
        elif arg[:2] == '-B':
            line_buffer = utils.float_or(arg[2:], line_buffer)
        elif arg == '--include_landmask' or arg == '-l':
            include_landmask = True
        elif arg == '--invert_watermask' or arg == 'i':
            invert_watermask = True
        elif arg == '-help' or arg == '--help' or arg == '-h':
            sys.stderr.write(_usage)
            sys.exit(1)
        elif arg == '-version' or arg == '--version':
            sys.stderr.write('{}\n'.format(_version))
            sys.exit(1)
        elif dst_ogr is None:
            dst_ogr = arg
        else:
            sys.stderr.write(_usage)
            sys.exit(1)
        i = i + 1

    this_region = regions.parse_cli_region([region])[0]
    this_cst = fetch_coastline(this_region)

    if this_region is None or not this_region.valid_p():
        sys.stderr.write(_usage)
        utils.echo_error_msg('invalid region')
        
    if this_cst is not None:
        with tqdm(
                total=len(this_cst.results),
                desc='processing coastline',
                leave=True
        ) as pbar:
            for n, cst_result in enumerate(this_cst.results):                
                if cst_result[-1] == 0:
                    pbar.update()
                    cst_osm = cst_result[1]
                    if dst_ogr is None:
                        dst_ogr = os.path.basename(utils.fn_basename2(cst_osm)) + '_coast.shp'
                        
                    out = fetches.polygonize_osm_coastline(
                        cst_osm, dst_ogr,
                        include_landmask = include_landmask,
                        landmask_is_watermask = invert_watermask,
                        line_buffer = line_buffer,
                    )
### End


