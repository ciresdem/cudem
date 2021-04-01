#!/usr/bin/env python
### spatial_metadata
##
## Copyright (c) 2019 - 2021 CIRES Coastal DEM Team
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
## generate spatial metadata from a datalist
##
### Code:

import os
import sys
from cudem import utils
from cudem import regions
from cudem import waffles

_version = '0.0.1'
_usage = '''spatial_metadata.py ({}): generate spatial metadata from a datalist

usage: spatial_metadata.py [ datalist [ OPTIONS ] ]

  datalist\t\tThe input datalist/entry

 Options:

  -o, --name\toutput name
  --help\tPrint the usage text
  --version\tPrint the version information

 Examples:
  % spatial_metadata.py my_data.datalist

CIRES DEM home page: <http://ciresgroups.colorado.edu/coastalDEM>
'''.format(_version)

if __name__ == '__main__':
    i = 1
    dls = []
    i_regions = []
    these_regions = []
    epsg = 4326
    inc = utils.str2inc('1s')
    node = 'pixel'
    name = 'waffles_spat'
    extend = 0
    want_verbose = True
    want_prefix = False

    argv = sys.argv
    while i < len(argv):
        arg = sys.argv[i]
        if arg == '--region' or arg == '-R':
            i_regions.append(str(argv[i + 1]))
            i = i + 1
        elif arg[:2] == '-R':
            i_regions.append(str(arg[2:]))
        elif arg == '--outname' or arg == '-O':
            name = argv[i + 1]
            i += 1
        elif arg[:2] == '-O': name = arg[2:]
        elif arg == '-s_epsg' or arg == '--s_epsg' or arg == '-P':
            epsg = argv[i + 1]
            i = i + 1
        elif arg == '--increment' or arg == '-E':
            inc = utils.str2inc(argv[i + 1])
            i = i + 1
        elif arg[:2] == '-E':
            inc = utils.str2inc(arg[2:])
        elif arg == '--extend' or arg == '-X':
            extend = utils.int_or(argv[i + 1], 0)
            i = i + 1
        elif arg[:2] == '-X':
            extend = utils.int_or(arg[2:], 0)
        elif arg == '-r' or arg == '--grid-node': node = 'grid'
        elif arg == '-p' or arg == '--prefix': want_prefix = True
        elif arg == '-help' or arg == '--help' or arg == '-h':
            sys.stderr.write(_usage)
            sys.exit(1)
        elif arg == '-version' or arg == '--version':
            sys.stderr.write('{}\n'.format(_version))
            sys.exit(1)
        else: dls.append(arg)
        i = i + 1

    for i_region in i_regions:
        tmp_region = regions.Region().from_string(i_region)
        if tmp_region.valid_p():
            these_regions.append(tmp_region)
        else:
            tmp_region = regions.ogr_wkts(i_region)
            for i in tmp_region:
                if i.valid_p():
                    these_regions.append(i)
                    
    if len(these_regions) == 0:
        these_regions = [None]
    else:
        if want_verbose: utils.echo_msg('parsed {} region(s)'.format(len(these_regions)))

    for rn, this_region in enumerate(these_regions):
        if len(dls) == 0:
            sys.stderr.write(_usage)
            utils.echo_error_msg('you must specify some type of data')
        else:
            if want_prefix or len(these_regions) > 1:
                name = waffles.waffles_append_fn(name, this_region, inc)
            [x for x in waffles.Waffle(data=dls, src_region=this_region, inc=inc, extend=extend, epsg=epsg, node=node, name=name, verbose=want_verbose).spat_meta(yxyz=False)]
        
### End
