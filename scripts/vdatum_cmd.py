#!/usr/bin/env python
### vdatum_cmd.py
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
### Code:

import os
import sys
from cudem import utils
from cudem import vdatumfun

_version = '0.0.4'
_usage = """vdatum_cmd.py ({}): run NOAAs vdatum

usage: vdatum [ file ]

  file\t\tThe input file to translate

 Options:

  -i the input vertical datum (default mllw:m:height) while
  -o specifies the output vertical datum (defafult navd88:m:height).

  -r specifies the input horizontal datum (default nad83:geo:deg) while
  -z specifies the output horizontal datum (default nad83:geo:deg).

  -g specifies the vdatum region (default 3)

  --help\tPrint the usage text
  --version\tPrint the version information

 Examples:
 % vdatum_cmd.py elev_data.xyz -i mllw -o navd88 -g 1
 % vdatum_cmd.py elev_data.xyz -i lmsl:ft:sounding -o navd88:m:height -g 1
 % vdatum_cmd.py elev_data.xyz -i navd88:m:height -o navd88:m:height -r nad83:utm:ft:10 -z nad83:geo:deg

See: vdatum.noaa.gov/docs/userguide_cmd.html for specifics on cli parameters.
CIRES DEM home page: <http://ciresgroups.colorado.edu/coastalDEM>
""".format(_version)

def main():
    src_fn = None
    ihorz = 'NAD83_2011'
    ohorz = 'NAD83_2011'
    ivert = 'navd88:m:height'
    overt = 'mhw:m:height'
    region = '3'
    delim = 'space'
    verbose = False

    i = 1

    argv = sys.argv
    while i < len(argv):
        arg = argv[i]

        if arg == '-i' or arg == '--ivert':
            ivert = argv[i + 1]
            i = i + 1
        elif arg == '-o' or arg == '--overt':
            overt = argv[i + 1]
            i = i + 1
        elif arg == '-r' or arg == '--ihorz':
            ihorz = argv[i + 1]
            i = i + 1
        elif arg == '-z' or arg == '--ohorz':
            ohorz = argv[i + 1]
            i = i + 1
        elif arg == '-e' or arg == '--region':
            region = argv[i + 1]
            i = i + 1
        elif arg == '-d' or arg == '--delimiter':
            delim = argv[i + 1]
            i = i + 1
        elif arg == '-V' or arg == '--verbose': verbose = True
        elif arg == '-help' or arg == '--help' or arg == '-h':
            print(_usage)
            sys.exit(1)
        elif arg == '-version' or arg == '--version':
            print('vdatum_cmd.py, version {}\nVDatum, version {}'.format(_version, utils.vdatum()._version))
            sys.exit(1)
        elif src_fn is None:
            src_fn = arg
        else:
            print(_usage)
            sys.exit(1)
        i = i + 1

    if src_fn is None:
        print(_usage)
        sys.exit(1)

    if not os.path.exists(src_fn):
        print('Error: {} is not a valid file'.format(src_fn))
    else:
        vdatumfun.Vdatum(ivert=ivert, overt=overt, ihorz=ihorz, ohorz=ohorz, delim=delim).run_vdatum(src_fn)
        
if __name__ == '__main__':
    main()

### End
