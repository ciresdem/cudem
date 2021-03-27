#!/usr/bin/env python
### gdal_outliers.py
##
## Copyright (c) 2021 CIRES Coastal DEM Team
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
## DD = D + (M/60) + (S/3600)
## DMS = int(DD):int((float(DD) - int(DD)) * 60):int((float(DD) - int(DD)) * 60) - float(int((float(DD) - int(DD)) * 60) * 60):
##
### Code:

import os
import sys

ddms_version = '0.1.2'

def dd2dms(dd):
    d = dd.split(".")[0]
    ms = float(dd) - int(d)
    
    mm = ms * 60
    
    m = str(mm).split(".")[0]
    s = (float(mm) - int(m)) * 60
    
    print('%s:%s:%s' %(d,m,s))

def dms2dd(dms):
    dms2 = dms.split(":")
    d = float(dms2[0])
    m = float(dms2[1])
    s = float(dms2[2])
    #print d,m,s
    dd = d + (m/60) + (s/3600)
    print(dd)

def Usage():
    print('Usage: ddms.py [-dd2dms] [-dms2dd]')
    print('               [-help] [-verbose] invalue')
    print('')
    print('ddms.py v.%s' %(ddms_version))

if __name__ == "__main__":
    
    inval=None
    isdd=True
    isdms=False

    i = 1
    while i < len(sys.argv):
        arg = sys.argv[i]

        if arg == '-help' or arg == '--help' or arg == '-h':
            Usage()
            sys.exit(0)

        elif arg == '-dd2dms':
            isdd = True
            isdms = False

        elif arg == '-dms2dd':
            isdms = True
            isdd = False

        elif arg == '-verbose':
            verbose = True

        elif inval is None:
            inval = arg

        elif arg[0] == '-':
            Usage()

        else:
            Usage()

        i = i + 1    

    if inval is None:
        inval = sys.stdin

        for line in inval:
            if isdd is True:
                dd2dms(line)
            if isdms is True:
                dms2dd(line)
    else:
        if isdd:
            dd2dms(inval)
        if isdms:
            dms2dd(inval)

### END
