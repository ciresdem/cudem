#!/usr/bin/env python
### colortable.py
##
## Copyright (c) 2011 - 2021 CIRES Coastal DEM Team
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
### Code:

import sys

_version = '0.1.1'

_usage = '''colortable.py ({}): generate an ETOPO1 colortable

usage: colortable.py [OPTIONS]

Options:
  -n, --min\tmin elevation value of colortable
  -x, --max\tmax elevation value of colortable

  --help\tprint the usage text
  --version\tPrint the version information

CIRES DEM home page: <http://ciresgroups.colorado.edu/coastalDEM>
'''.format(_version)

trs = (-11000, -10500, -10000, -9500, -9000, -8500,
       -8000, -7500, -7000, -6500, -6000, -5500,
       -5000, -4500, -4000, -3500, -3000, -2500,
       -2000, -1500, -1000, -500, -0.001, 0,
       100, 200, 500, 1000, 1500, 2000, 2500,
       3000, 3500, 4000, 4500, 5000, 5500,
       6000, 6500, 7000, 7500, 8000)

colors = ([10,0,121], [26,0,137], [38,0,152], [27,3,166], [16,6,180],
          [5,9,193], [0,14,203], [0,22,210], [0,30,216], [0,39,223],
          [12,68,231], [26,102,240], [19,117,244], [14,133,249], [21,158,252],
          [30,178,255], [43,186,255], [55,193,255], [65,200,255], [79,210,255],
          [94,223,255], [138,227,255], [138,227,255], [51,102,0], [51,204,102],
          [187,228,146], [255,220,185], [243,202,137], [230,184,88], [217,166,39],
          [168,154,31], [164,144,25], [162,134,19], [159,123,13], [156,113,7],
          [153,102,0], [162,89,89], [178,118,118], [183,147,147], [194,176,176],
          [204,204,204], [229,229,229], [138,227,255], [51,102,0])
          
elevs = [-20, -19.09090909, -18.18181818, -17.27272727, -16.36363636,
         -15.45454545, -14.54545455, -13.63636364, -12.72727273, -11.81818182,
         -10.90909091, -10, -9.090909091, -8.181818182, -7.272727273,
         -6.363636364, -5.454545455, -4.545454545, -3.636363636, -2.727272727,
         -1.818181818, -0.909090909, -0.001, 0, 48.06865898,
         96.13731797, 240.3432949, 480.6865898, 721.0298848, 961.3731797,
         1201.716475, 1442.05977, 1682.403064, 1922.746359, 2163.089654,
         2403.432949, 2643.776244, 2884.119539, 3124.462834, 3364.806129,
         3605.149424,  3845.492719]

def scale_el(value, min, max, tr):
    if value > 0 and max > 0:
        return((max * tr) / 8000)
    elif value < 0 and min < 0:
        return((min * tr) / -11000)
    elif value == 0:
        return(0)
    else: return(None)

if __name__ == "__main__":
    gmin = -100
    gmax = 100
    i = 1

    while i < len(sys.argv):
        arg = sys.argv[i]
        if arg == '-n' or arg == '--min':
            gmin = float(sys.argv[i+1])
            i = i + 1
        elif arg == '-x' or arg == '--max':
            gmax = float(sys.argv[i+1])
            i = i + 1
        elif arg == '-help' or arg == '--help' or arg == '-h':
            sys.stderr.write(_usage)
            sys.exit(1)
        elif arg == '-version' or arg == '--version':
            sys.stderr.write('{}\n'.format(_version))
            sys.exit(1)
        else:
            sys.stderr.write(_usage)
            sys.exit(1)
        i = i + 1    

    for i,j in enumerate(elevs):
        if j != None and i + 1 != len(elevs):
            elev = scale_el(j, gmin, gmax, trs[i])
            elev1 = scale_el(elevs[i + 1], gmin, gmax, trs[i + 1])
            print(elev, colors[i][0], colors[i][1], colors[i][2], elev1, colors[i][0], colors[i][1], colors[i][2])

### End
