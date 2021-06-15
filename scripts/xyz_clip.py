#!/usr/bin/env python
#
#  Description: 
#
#--

#--
import sys
import csv
import math
import struct
import os
import numpy as np
from osgeo import gdal
from cudem import utils
#--

xc_version = '1.3.3'

#--
def xyzClip(inxyz, ingrd, outfile, outfile2, delim, msk_nval, msk_kval, inverse, keep_out):
    '''Compare an xyz infile to a gdal compatible
    ingrd, and output a new xyz outfile.'''

    if inxyz == sys.stdin:
        in_plot = inxyz
    else:
        in_plot = open(inxyz, 'r')

    ## Process the grid file
    ds = gdal.Open(ingrd)
    comp_geot = ds.GetGeoTransform()
    band = ds.GetRasterBand(1)
    nodata = band.GetNoDataValue()
    if verbose:
        print >> sys.stderr, "xyz_clip: " + str(comp_geot)

    cellsize = [float(comp_geot[1]), float(comp_geot[5])]
    xextent = float(comp_geot[0])
    yextent = float(comp_geot[3])
    if outfile is not None:
        writer = csv.writer(open(outfile, 'w'), delimiter=delim, lineterminator='\n')
    if outfile2 is not None:
        writer2 = csv.writer(open(outfile2, 'w'), delimiter=delim, lineterminator='\n')
    pointnum = 0
    
    ## Iterate through the xyz file and compare each point to the input grid file.
    # for record in csv.reader(open(infile, 'rb'), delimiter=delim, skipinitialspace=True):
    for i in in_plot:
        try:
            x = float(i.split(delim)[0].strip())
            y = float(i.split(delim)[1].strip())
            z = float(i.split(delim)[2].strip())
        except:
            print >> sys.stderr, "xyz_clip: Failed to read line"
            x=xextent - 10
            y=yextent + 10
        
        ## Continue if values are reasonable.
        if x > xextent and y < yextent:
            xpos = int(math.fabs(math.ceil((xextent - x)/float(cellsize[0]))))
            ypos = int(math.fabs(math.ceil((y - yextent)/float(cellsize[1]*-1))))
            #xpos, ypos = utils._geo2pixel(x, y, band.)

            scanline = band.ReadRaster(0,ypos,band.XSize, 1, band.XSize, 1, gdal.GDT_Float32)
            these_values = struct.unpack('f' * band.XSize, scanline)
            
            try:
                grdz = these_values[xpos]
                if inverse is True:
                    
                    if grdz == nodata:
                        zd = nodata
                    elif grdz == msk_kval:
                        zd = nodata
                    else:
                        zd = msk_kval
                else:
                    if grdz == nodata:
                        zd = nodata
                    elif grdz == msk_nval:
                        zd = nodata
                    else:
                        zd = msk_kval
            except:
                zd = nodata
                
            if zd == 1:
                # Write the results to the new file.
                if outfile is not None:
                    row = [x,y,z]
                    writer.writerow(row)
                else:
                    #print x,y,z
                    print(delim.join(map(str, [x,y,z])))
            else:
                # Write discared results to file
                if outfile2 is not None:
                    row = [x,y,z]
                    writer2.writerow(row)
                
            if verbose:
                print >> sys.stderr, "xyz_clip: " + str(pointnum) + ' points scanned...\r',
                pointnum += 1
        elif keep_out:
            print >> sys.stderr, "xyz_clip: Out!\r",
            if outfile is not None:
                row = [x,y,z]
                writer.writerow(row)
            else:
                print(delim.join(map(str, [x,y,z])))
                

#--
def xyz_clip(inxyz, ingrd, outfile, outfile2, delim, msk_nval, msk_kval, inverse, keep_out):
    '''Compare an xyz infile to a gdal compatible
    ingrd, and output a new xyz outfile.'''

    if inxyz == sys.stdin:
        in_plot = inxyz
    else:
        in_plot = open(inxyz, 'r')

    ## Process the grid file
    ds = gdal.Open(ingrd)
    comp_geot = ds.GetGeoTransform()
    band = ds.GetRasterBand(1)
    nodata = band.GetNoDataValue()
    if verbose:
        print >> sys.stderr, "xyz_clip: " + str(comp_geot)

    try: tgrid = band.ReadAsArray()
    except: xyzClip(inxyz, ingrd, outfile, outfile2, delim, msk_nval, msk_kval, inverse)

    cellsize = [float(comp_geot[1]), float(comp_geot[5])]
    xextent = float(comp_geot[0])
    yextent = float(comp_geot[3])
    if outfile is not None:
        writer = csv.writer(open(outfile, 'w'), delimiter=delim, lineterminator='\n')
    if outfile2 is not None:
        writer2 = csv.writer(open(outfile2, 'w'), delimiter=delim, lineterminator='\n')
    pointnum = 0
    
    ## Iterate through the xyz file and compare each point to the input grid file.
    # for record in csv.reader(open(infile, 'rb'), delimiter=delim, skipinitialspace=True):
    for i in in_plot:

        try:
            x = float(i.split(delim)[0].strip())
            y = float(i.split(delim)[1].strip())
            z = float(i.split(delim)[2].strip())
        except:
            print >> sys.stderr, "xyz_clip: Failed to read line"
            x=xextent - 10
            y=yextent + 10
        
        ## Continue if values are reasonable.
        if x > xextent and y < yextent:
            xpos = int(math.fabs(math.ceil((xextent - x)/float(cellsize[0]))))
            ypos = int(math.fabs(math.ceil((y - yextent)/float(cellsize[1]*-1))))
            try:
                grdz = tgrid[ypos,xpos]
                if inverse is True:
                    
                    if grdz == nodata:
                        zd = nodata
                    elif grdz == msk_kval:
                        zd = nodata
                    else:
                        zd = msk_kval
                else:
                    if grdz == nodata:
                        zd = nodata
                    elif grdz == msk_nval:
                        zd = nodata
                    else:
                        zd = msk_kval
            except:
                if keep_out:
                    zd = msk_kval
                else:
                    zd = nodata
                
            if zd == msk_kval:
                # Write the results to the new file.
                if outfile is not None:
                    row = [x,y,z]
                    writer.writerow(row)
                else:
                    #print x,y,z
                    print(delim.join(map(str, [x,y,z])))
            else:
                # Write discared results to file
                if outfile2 is not None:
                    row = [x,y,z]
                    writer2.writerow(row)
                
            if verbose:
                print >> sys.stderr, "xyz_clip: " + str(pointnum) + 'points scanned...\r',
                pointnum += 1

        elif keep_out:
            if outfile is not None:
                row = [x,y,z]
                writer.writerow(row)
            else:
                print(delim.join(map(str, [x,y,z])))

#--
def Usage():
    print('Usage: xyz_clip.py [-delimiter char] [-dual outfile] [-msk_keep value]')
    print('                   [-msk_remove value] [-invert] [-quickly] [-verbose]')
    print('                   [-return_all] in_xyz msk_grid [outfile]')
    print('')
    print('xyz_clip v.%s' %(xc_version))

#--
#
# Mainline
#
#--
if __name__ == '__main__':

    infile = None
    outfile = None
    mskgrd = None
    delim = " "
    t_format = ""
    msk_kval = 1
    msk_nval = 0
    invert = False
    dual_out = None
    verbose = False
    quickly = False
    overwrite = False
    rall = False

    gdal.AllRegister()
    argv = gdal.GeneralCmdLineProcessor( sys.argv )
    if argv is None:
        sys.exit(0)
        
    # Parse command line arguments.
    i = 1
    while i < len(argv):
        arg = argv[i]

        if arg == '-delimiter':
            delim = str(argv[i+1])
            i = i + 1

        elif arg == '-dual':
            dual_out = str(argv[i+1])
            i = i + 1

        elif arg == '-msk_keep':
            msk_kval = float(argv[i+1])
            i = i + 1

        elif arg == '-msk_remove':
            msk_nval = float(argv[i+1])
            i = i + 1

        elif arg == '-all':
            rall = True

        elif arg == '-invert':
            invert = True

        elif arg == '-verbose':
            verbose = True

        elif arg == '-quickly':
            quickly = True

        elif arg == '-return_all':
            rall = True

        elif arg[0] == '-':
            Usage()

        elif infile is None:
            infile = arg

        elif mskgrd is None:
            mskgrd = arg

        elif outfile is None:
            outfile = arg

        else:
            Usage()

        i = i + 1
    
    if infile is None:
        infile = sys.stdin

    if mskgrd is None:
        Usage()
        sys.exit(0)

    if quickly:
        xyz_clip(infile, mskgrd, outfile, dual_out, delim, msk_nval, msk_kval, invert, rall)
    else:
        xyzClip(infile, mskgrd, outfile, dual_out, delim, msk_nval, msk_kval, invert, rall)
### End
