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

_usage = '''fetch_osm_coastline.py ({}): fetch and process the OSM coastline

usage: fetch_osm_coastline.py [-R]

 Options:
  -R, --region\t\tRestrict processing to the desired REGION 
\t\t\tWhere a REGION is xmin/xmax/ymin/ymax[/zmin/zmax[/wmin/wmax/umin/umax]]
\t\t\tUse '-' to indicate no bounding range; e.g. -R -/-/-/-/-10/10/1/-/-/-
\t\t\tOR an OGR-compatible vector file with regional polygons. 
\t\t\tWhere the REGION is /path/to/vector[:zmin/zmax[/wmin/wmax/umin/umax]].
\t\t\tIf a vector file is supplied, will use each region found therein.

  --help\tPrint the usage text
  --version\tPrint the version information

 Examples:
 % fetch_osm_coastline.py -R -123.5/-123.25/38.5/38.75

'''.format(_version)

if __name__ == '__main__':    
    region = None
    remove_value = None
    i = 1
    while i < len(sys.argv):
        arg = sys.argv[i]
        if arg == '-help' or arg == '--help' or arg == '-h':
            sys.stderr.write(_usage)
            sys.exit(1)
        elif arg == '-version' or arg == '--version':
            sys.stderr.write('{}\n'.format(_version))
            sys.exit(1)
        elif arg == '--region' or arg == '-R':
            region = str(argv[i + 1])
            i = i + 1
        elif arg[:2] == '-R':
            region = str(arg[2:])
        else:
            sys.stderr.write(_usage)
            sys.exit(1)
        i = i + 1

    this_region = regions.parse_cli_region([region])[0]
    this_cst = fetch_coastline(this_region)

    if this_region is None or not this_region.valid_p():
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
                    out = gdalfun.ogr_polygonize_line_to_region(
                        cst_osm, os.path.basename(utils.fn_basename2(cst_osm)) + '_coast.shp', include_landmask = True
                    )
### End


