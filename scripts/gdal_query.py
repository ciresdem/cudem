#!/usr/bin/env python
### gdal_query.py
##
## Copyright (c) 2011 - 2021 CIRES DEM Team
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
## Compare the values found in a source 'xyz' ascii
## file to those of a gdal compatible gridfile, or
## sample a gridfile to xy points in an xy* file.
##
## Depends: GDAL ; NumPy
### Code:

#--
import sys
import os
import math
import re
import struct
import numpy as np
import osgeo.gdal as gdal
#--

gq_version = '1.8.7'

def con_dec(x, dec):
    '''Return a float string with n decimals
    (used for ascii output).'''

    if x is None:
        print >> sys.stderr, "gdal_query: Error, Attempting to convert a 'None' value."
        return
    fstr = "%." + str(dec) + "f"
    return fstr % x

# Convert a geographic x,y value to a pixel location of geoTransform
def geo2pixel( geo_x, geo_y, geoTransform ):
    if geoTransform[2] + geoTransform[4] == 0:
        pixel_x = (geo_x - geoTransform[0]) / geoTransform[1]
        pixel_y = (geo_y - geoTransform[3]) / geoTransform[5]
    else:
        pixel_x, pixel_y = apply_gt( geo_x, geo_y, invert_gt( geoTransform ) )
    return int(pixel_x), int(pixel_y)

# Convert a pixel location to geographic coordinates given geoTransform
def pixel2geo( pixel_x, pixel_y, geoTransform ):
    geo_x, geo_y = apply_gt( pixel_x, pixel_y, geoTransform )
    return geo_x, geo_y

def apply_gt( in_x, in_y, geoTransform ):
    out_x = geoTransform[0] + in_x * geoTransform[1] + in_y * geoTransform[2]
    out_y = geoTransform[3] + in_x * geoTransform[4] + in_y * geoTransform[5]
    return out_x, out_y

def invert_gt(geoTransform):
    det = geoTransform[1] * geoTransform[5] - geoTransform[2] * geoTransform[4]
    if abs(det) < 0.000000000000001:
        return
    invDet = 1.0 / det
    # compute adjoint and divide by determinate
    outGeoTransform = [0,0,0,0,0,0]
    outGeoTransform[1] = geoTransform[5] * invDet
    outGeoTransform[4] = -geoTransform[4] * invDet
    outGeoTransform[2] = -geoTransform[2] * invDet
    outGeoTransfrom[5] = geoTransform[1] * invDet
    outGeoTransform[0] = ( geoTransform[2] * geoTransform[3] - geoTransform[0] * geoTransform[5] ) * invDet
    outGeoTransform[3] = ( -geoTransform[1] * geoTransform[3] + geoTransform[0] * geoTransform[4] ) * invDet
    return outGeoTransform 

def query_gdal(inxyz,ingrd,delim,xloc,yloc,zloc,nodata,out_form,addit,return_all,verbose):
    '''Compare an xyz file to a gdal-compatible grid file.'''

    if inxyz == sys.stdin:
        in_plot = inxyz
    else:
        in_plot = open(inxyz, 'r')

    # Process the grid file
    ds = gdal.Open(ingrd)
    comp_geot = ds.GetGeoTransform()
    band = ds.GetRasterBand(1)

    # Load the grid into a numpy array
    tgrid = band.ReadAsArray()
    #except: xyz_vs_gdal(inxyz,ingrd,delim,xloc,yloc,zloc,nodata,out_form,addit,return_all,verbose)

    # Get some grid info
    if nodata is None: nodata = band.GetNoDataValue()
    if nodata is None: nodata = -9999

    cellsize = [float(comp_geot[1]), float(comp_geot[5])]
    xextent = float(comp_geot[0])
    yextent = float(comp_geot[3])

    # Process the xyz file
    n = 0
    for i in in_plot:
        n+=1
        if verbose:
            sys.stderr.write('gdal_query:' + str(n) + '\r')
        try:
            x = float(i.split(delim)[int(xloc)].strip())
            y = float(i.split(delim)[int(yloc)].strip())
            if zloc == '-':
                z = nodata
            else:
                z = float(i.split(delim)[int(zloc)].strip())
    
        except:
            sys.stderr.write("gdal_query: Failed to read line: {}".format(i))
            x=xextent - 10
            y=yextent + 10
        
        # Continue if values are reasonable.
        if x > xextent and y < yextent:
            xpos,ypos = geo2pixel(x,y,comp_geot)
            #xpos = int(math.fabs(math.ceil((xextent - x)/float(cellsize[0]))))
            #ypos = int(math.fabs(math.ceil((y - yextent)/float(cellsize[1]*-1))))
       
            # Locate the grid cell and get it's value
            try: g = tgrid[ypos,xpos]
            except: g = nodata

            #print g

            d = nodata
            c = nodata
            m = nodata
            s = nodata
            
            if g != nodata:
                d = z - g
                m = z + g
                # FIXME ## as d/g*100 woudl fail if g was zero
                #if g == 0:
                #    g += 0.0000001
                # /FIXME ##
                # c is the percent difference
                c = con_dec(math.fabs(float(d/(g+0.00000001)*100)), 2)
                # s is the 'scaled' difference
                s = con_dec(math.fabs(d / (z + (g+0.00000001))), 4)

            # Print out the results
            d = con_dec(d, 4)

            outs = []
            for i in out_form:
                outs.append(vars()[i])
            print(delim.join(map(str, outs)))
        else:
            if return_all:
                d = nodata
                m = nodata
                c = nodata
                s = nodata
                g = nodata
                outs = []
                for i in out_form:
                    outs.append(vars()[i])
                print(delim.join(map(str, outs)))
    if verbose:
        sys.stderr.write('\n')

def Usage():
    print('Usage: gdal_query.py [-delimiter char] [-s_format "0,1,2"] ')
    print('                     [-d_format "xyzgdps"] [-d_nodata value] [-header]')
    print('                     [-return_all] [-verbose]')
    print('                     grdfile [srcfile]')
    print('')
    print('Input options:')
    print('  grdfile\tSpecifies which gdal-compatible grid file in which to ')
    print('  \t\tquery values.')
    print('  srcfile\tSpecifies which xy* file in which to query grdfile with,')
    print('  \t\tif not given, will read from standard input.')
    print('  -delimiter\tSpecifies the input xy* delimiter')
    print('  -s_format\tSpecifies the input xy* format, a quoted triplet specifying')
    print('  \t\twhere the x,y and z values are located.')
    print('  \t\tIf there is no input z value, specify with a -, i.e. ')
    print('  \t\t-s_format "0,1,-"')
    print('')
    print('Output options:')
    print('  -d_format\tSpecifies the output xy* format, a quoted string specifying ')
    print('  \t\twhich output values to return.')
    print('  \t\tOptions are [x][y][z][g][d][m][c][s]; where z is the input z value, g is the grid ')
    print('  \t\tz value, d is the difference between the input z value and grid')
    print('  \t\tz value, m is the sum of the input z value and the grid z value,')
    print('  \t\tc is the percetage difference and s is the scaled difference.')
    print('  -d_nodata\tSpecifies the output nodata value, will use the nodata value from')
    print('  \t\tthe input grid if not specified.')
    print('  -header\tIndicates that the output will have a header included, default is False.')
    print('  -return_all\tWill prompt the return of all input points, regardless of location')
    print('  \t\twhen compared to the input grid file, with grid values considered as nodata.')
    print('')
    print('General options:')
    print('  -verbose\tWill increase verbosity.')
    print('')
    print('Example:')
    print('To query grd.tif with an ascii xy file (values.xy) with comma delimited fields and ')
    print('no z-value and return the xy values from the input xy file and the z values')
    print('from grd.tif:')
    print('gdal_query.py -delimiter "," -s_format "0,1,-" -d_format "xyg" grd.tif values.xy > values.xyz')
    print('\ngdal_query v.%s' %(gq_version))
    sys.exit( 1 )
#
# Mainline
#
if __name__ == "__main__":

    inxyz=None
    ingrd=None
    delim=" "
    xyzf="0,1,2"
    out_xyzf="xyzg"
    xyzh=False
    addit=False
    verbose=False
    return_all=False
    out_nodata=None

    gdal.AllRegister()
    argv = gdal.GeneralCmdLineProcessor( sys.argv )
    if argv is None:
        sys.exit(0)
        
    # Parse command line arguments.
    i = 1
    while i < len(argv):
        arg = argv[i]

        if arg == '-delimiter':
            delim = argv[i+1]
            i = i + 1

        elif arg == '-s_format':
            xyzf = str(argv[i+1])
            i = i + 1

        elif arg == '-d_format':
            out_xyzf = str(argv[i+1])
            i = i + 1

        elif arg == '-d_nodata':
            out_nodata = float(argv[i+1])
            i = i + 1
            
        elif arg == '-header':
            xyzh = True

        elif arg == '-addition':
            addit = True

        elif arg == '-verbose':
            verbose = True

        elif arg == '-return_all':
            return_all = True

        elif arg[0] == '-':
            Usage()

        elif ingrd is None:
            ingrd = arg

        elif inxyz is None:
            inxyz = arg

        else:
            Usage()

        i = i + 1

    if inxyz == None:
        inxyz = sys.stdin
    if ingrd == None:
        Usage()
        sys.exit(0)

    #-- Parse point locations.
    xloc = xyzf.split(",")[0].strip()
    yloc = xyzf.split(",")[1].strip()
    zloc = xyzf.split(",")[2].strip()

    out_form = tuple(re.findall(r'(\D)',out_xyzf))

    if xyzh == True:
        outs = []
        for i in out_form:
            outs.append(i)
        print(delim.join(outs))

    query_gdal(inxyz,ingrd,delim,xloc,yloc,zloc,out_nodata,out_form,addit,return_all,verbose)
#--END
