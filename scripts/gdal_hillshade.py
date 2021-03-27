#!/usr/bin/env python

import sys
import math
import os.path
import optparse
import numpy

try: import osgeo.gdal as gdal
except: sys.exit('Error - You must have the GDAL Python bindings installed')

#--
def Usage():
    print('Usage: xyz_clip.py [-delimiter char] [-dual outfile] [-msk_keep value]')
    print('                   [-msk_remove value] [-invert] [-quickly] [-verbose]')
    print('                   in_xyz msk_grid [outfile]')
    print('')
    print('xyz_clip v.%s' %(xc_version))

NODATA = -1

#def persp(img):

def arr2img(ar):
    '''Convert Numeric.array to PIL.Image.'''

    return PIL.Image.fromstring('L', (ar.shape[1], ar.shape[0]), ar.astype('b').tostring())

def hillshade(cell, xres, yres, nodata=None, azimuth=315.0, altitude=45.0, z=1.0, scale=1.0):
    ''' cell is an array of elevation values. return array will have outermost pixel edge removed.
        logic borrowed from hillshade.cpp, http://www.perrygeo.net/wordpress/?p=7 '''

    if verbose:
        print >> sys.stderr, 'making window...'
    
    window = []
    
    for row in range(3):
        for col in range(3):
            print("cell.shape[0]:" + str(cell.shape))
            print("cell.shape[1]" + str(cell.shape[1]))
            window.append(cell[row:(row + cell.shape[0] - 2), col:(col + cell.shape[1] - 2)])
    
    if verbose:
        print >> sys.stderr, 'calculating slope...'
    
    x = ((z * window[0] + z * window[3] + z * window[3] + z * window[6]) \
       - (z * window[2] + z * window[5] + z * window[5] + z * window[8])) \
      / (8.0 * xres * scale);
    
    y = ((z * window[6] + z * window[7] + z * window[7] + z * window[8]) \
       - (z * window[0] + z * window[1] + z * window[1] + z * window[2])) \
      / (8.0 * yres * scale);

    rad2deg = 180.0 / math.pi
    
    slope = 90.0 - numpy.arctan(numpy.sqrt(x*x + y*y)) * rad2deg
    
    if verbose:
        print >> sys.stderr, 'calculating aspect...'
    
    aspect = numpy.arctan2(x, y)
    
    if verbose:
        print >> sys.stderr, 'calculating shade...'
    
    deg2rad = math.pi / 180.0
    
    shaded = numpy.sin(altitude * deg2rad) * numpy.sin(slope * deg2rad) \
           + numpy.cos(altitude * deg2rad) * numpy.cos(slope * deg2rad) \
           * numpy.cos((azimuth - 90.0) * deg2rad - aspect);

    shaded = shaded * 255
    
    if nodata is not None:
        if verbose:
            print >> sys.stderr, 'removing no-data...'
    
        for pane in window:
            shaded[pane == nodata] = NODATA

    if verbose:
        print >> sys.stderr, 'done.'
    
    return shaded

def getbase(azimuth=315.0, altitude=45.0):
    ''' Shade a flat piece of ground to determine what its color is'''

    cell = numpy.ones((3, 3), dtype=numpy.float32)
    shaded = hillshade(cell, 1, 1, None, azimuth, altitude).astype(numpy.ubyte)
    return shaded[0, 0]

if __name__ == '__main__':

    parser = optparse.OptionParser()
    
    parser.add_option('-i', '--input', dest='input',
                      help='Input file')
    
    parser.add_option('-o', '--output', dest='output',
                      help='Output file')

    parser.add_option('-v', '--verbose', action='store_true', dest='verbose', default=False,
                      help='Increase Verbosity')

    (options, args) = parser.parse_args()
    
    input, output = options.input, options.output
    verbose = options.verbose

    if input == None or output == None:
        parser.print_help()
        sys.exit("Error - You failed to supply an Input and/or Output file")

    if verbose:
        print >> sys.stderr, 'importing dem...'

    dem = gdal.Open(input)
    band = dem.GetRasterBand(1)
    cols, rows = dem.RasterXSize, dem.RasterYSize
    
    if verbose:
        print(cols, 'x', rows,)
        print('floor', getbase())
    
        print >> sys.stderr, 'extracting data...'
    
    data = band.ReadRaster(0, 0, cols, rows, buf_type=gdal.GDT_Float32)
    cell = numpy.fromstring(data, dtype=numpy.float32).reshape(rows, cols)
    
    if verbose:
        print >> sys.stderr, 'shading hills...'
    
    tx = dem.GetGeoTransform()
    assert tx[2] == 0 and tx[4] == 0
    
    xres, yres = tx[1], tx[5]
    
    shaded = hillshade(cell, xres, yres, band.GetNoDataValue())
    
    ## this part is unfortunate, wish it wasn't necessary, seems broken
    #shaded[shaded == NODATA] = getbase()
    
    if verbose:
        print >> sys.stderr, 'compositing image...'
    
        print >> sys.stderr, 'saving jpeg...'
    
    out = numpy.clip(shaded, 0x00, 0xFF).astype(numpy.ubyte)
    #jpg = arr2img(out)
    #jpg.save('out.jpg')
    
    if verbose:
        print >> sys.stderr, 'saving geotiff...'
    
    driver = gdal.GetDriverByName('GTiff')
    
    demtx = list(dem.GetGeoTransform())
    
    # account for a pixel of choke
    demtx[0], demtx[3] = demtx[0] + demtx[1], demtx[3] + demtx[5]
    
    tif = driver.Create(output, shaded.shape[1], shaded.shape[0], 1, gdal.GDT_Int16)
    tif.SetGeoTransform(demtx)
    tif.SetProjection(dem.GetProjection())
    
    tifband = tif.GetRasterBand(1)
    tifband.SetNoDataValue(NODATA)
    
    tifdata = shaded.astype(numpy.int16).tostring()
    tifband.WriteRaster(0, 0, tif.RasterXSize, tif.RasterYSize, tifdata, buf_type=gdal.GDT_Int16)
