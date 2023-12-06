#!/usr/bin/env python
#
# xyz2shp.py
#
# Purpose: Convert xyz text files into ESRI Shapefiles using GDAL/OGR
#
#--

#--IMPORTS
import sys
import os

try:
    import osgeo.ogr as ogr
    import osgeo.gdal as gdal
except ImportError:
    try:
        import ogr
        import gdal
    except ImportError:
        sys.exit('''
xyz2gdal: Sorry: You must have the Python GDAL/OGR bindings for Shapefile support,
Get them here: http://trac.osgeo.org/gdal/wiki/GdalOgrInPython''')
#--

xsversion = '1.3.1'

#--
def checkColumnNumType(inline, delim):	
    colNumType = []
    cols = inline.split(delim)
    for i in range(0, len(cols)):
        try: 
            float(cols[i])
            colgrp = 1
        except:
            colgrp = 0
        colNumType.append(colgrp)
    return colNumType


def proc_xyz(rr, delim, xposition, yposition, zposition, coltypes, f, dim, layer):
    record1 = rr.split(delim)
    #i += 1

    try:
        x = float(record1[xposition].strip())
        y = float(record1[yposition].strip())
        z = float(record1[zposition].strip())
        
        record = map(float,[x,y,z])
        #if verbose:
        #    print('xyz2shp: %s \r' %(i)),
            
        field = 0
        for cell in record:

            if coltypes[field] == 1:
                # try: this_rec = int(cell)
                this_rec = float(cell)
            else:
                this_rec = str(cell)
                
            f.SetField(field, this_rec)
            field += 1

        if dim == 3:
            wkt = 'POINT(%f %f %f)' % (x,y,z)
        else:
            wkt = 'POINT(%f %f)' % (x,y)

        g = ogr.CreateGeometryFromWkt(wkt)
        f.SetGeometryDirectly(g)
        layer.CreateFeature(f)

    except:
        print("xyz2shp: Error: Failed to create shapefile geometry for record: %s please check your source-file and/or command-line options." % (str(rr)))

#--
def xyz2shp(infile,outfile,delim,xposition,yposition,zposition,overwrite,verbose):

    if infile == sys.stdin:
        inxyz = infile
    else:
        inxyz = open(infile, 'r')

    first_line = inxyz.readline()
    coltypes = checkColumnNumType(first_line, delim)
        
    driver = ogr.GetDriverByName('ESRI Shapefile')
    if os.path.exists(outfile):
        if overwrite:
            driver.DeleteDataSource(outfile)
        else:
            print("xyz2shp: Error, the desired output file exists, use -overwrite to overwrite it.")
    ds = driver.CreateDataSource(outfile)

    dim = 3
    
    if dim == 3:
        layer = ds.CreateLayer(outfile, geom_type=ogr.wkbPoint25D)
    else:
        layer = ds.CreateLayer(outfile, geom_type=ogr.wkbPoint)

    for i,j in enumerate(coltypes):

        oft_title = "column_" + str(i)
		
        if i == int(xposition):
            oft_title = "Longitude"
            oft_string = ogr.OFTReal
            fdw = 12
            fdp = 8
        elif i == int(yposition):
            oft_title = "Latitude"
            oft_string = ogr.OFTReal
            fdw = 12
            fdp = 8
        elif i == int(zposition):
            oft_title = "Elevation"
            oft_string = ogr.OFTReal
            fdw = 12
            fdp = 8
        elif j == 1:
            oft_string = ogr.OFTReal
            fdw = 10
            fdp = 4
        else:
            oft_string = ogr.OFTString
            fdw = 25
            fdp = 10

        fd = ogr.FieldDefn(oft_title, oft_string)
        fd.SetWidth(fdw)
        fd.SetPrecision(fdp)
        layer.CreateField(fd)

    f = ogr.Feature(feature_def=layer.GetLayerDefn())
    #i = 0

    #--
    # Process the input xyz file
    if first_line:
        proc_xyz(first_line, delim, xposition, yposition, zposition, coltypes, f, dim, layer)

    for rr in inxyz:
        proc_xyz(rr, delim, xposition, yposition, zposition, coltypes, f, dim, layer)
        
    f.Destroy()
    ds.Destroy()
    if verbose:
        print("\nxyz2shp: Shapefile Created")

#--
def writePrjFile(outfile, epsg_code):
    from osgeo import osr
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(int(epsg_code))
    prj = srs.ExportToWkt()
    otf = open(str(outfile[:-4] + ".prj"), "w")
    otf.write(prj)
    otf.close()
#--

#--
def Usage():
    print('Usage: xyz2shp.py [-d/-delimiter char] [-s/-s_format 0,1,2] [-p/-d_epsg epsg_code]')
    print('                  [-verbose] [-overwrite] outfile [srcfile]')
    print('')
    print('xyz2shp v. %s' %(xsversion))

#--
#
# Mainline
#
#--
if __name__ == '__main__':

    infile = None
    outfile = None
    delim = " "
    xyzf = "0,1,2"
    epsg_code = 4326
    verbose = False
    overwrite = False

    gdal.AllRegister()
    argv = gdal.GeneralCmdLineProcessor( sys.argv )
    if argv is None:
        sys.exit(0)
        
    # Parse command line arguments.
    i = 1
    while i < len(argv):
        arg = argv[i]

        if arg == '-delimiter' or arg == '-d':
            delim = str(argv[i+1])
            i = i + 1

        elif arg == '-s_format' or arg == 's':
            xyzf = argv[i+1]
            i = i + 1

        elif arg == '-d_epsg' or arg == '-p':
            epsg_code = int(argv[i+1])
            i = i + 1

        elif arg == '-verbose':
            verbose = True

        elif arg == '-overwrite':
            overwrite = True

        elif arg[0] == '-':
            Usage()

        elif outfile is None:
            outfile = arg

        elif infile is None:
            infile = arg

        else:
            Usage()

        i = i + 1

    if infile is None:
        infile = sys.stdin

    if outfile is None:
        Usage()
        sys.exit(0)

    xloc = int(xyzf.split(",")[0].strip())
    yloc = int(xyzf.split(",")[1].strip())
    zloc = int(xyzf.split(",")[2].strip())
	
    xyz2shp(infile, outfile, delim, xloc, yloc, zloc, overwrite, verbose)
    writePrjFile(outfile,epsg_code)
#--END
