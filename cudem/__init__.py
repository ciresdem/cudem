## Copyright (c) 2012 - 2025 Regents of the University of Colorado
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

__version__ = "2.4.5"
__author__ = "Matthew Love"
__credits__ = "CIRES"

## Windows support

import os
import sys
from osgeo import gdal
from cudem import utils
from cudem import regions
from cudem import xyzfun
from cudem import gdalfun
from cudem import factory
from cudem import vdatums
from cudem import fetches
from cudem import grits
from cudem import vrbag
from cudem import waffles
from cudem import dlim
from cudem import htdpfun
from cudem import perspecto
from cudem import srsfun

gc = utils.config_check() # cudem config file holding foriegn programs and versions
if gc['platform'] == 'linux':
    #gdal.SetConfigOption('CPL_LOG', '/dev/null') # supress gdal warnings in linux
    pass
else:
    os.system("") # ansi in windows
    gdal.SetConfigOption('CPL_LOG', 'NUL') # supress gdal warnings in windows

__all__ = ['utils', 'regions', 'xyzfun', 'gdalfun', 'factory', 'vdatums',
           'fetches', 'grits', 'vrbag', 'waffles', 'dlim', 'htdpfun',
           'perspecto', 'srsfun']
    
## user module
sys.path.append(os.path.expanduser('~'))
try:
    import cudemrc
    __all__.append('cudemrc')
except:
    pass

#from archook import locate_arcgis, get_arcpy

### End
