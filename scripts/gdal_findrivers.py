#!/usr/bin/env python
#
# Description: Mask data in a gdal grid
#
#--

import os
import sys

import numpy as np
import scipy
from scipy import ndimage
from tqdm import tqdm
from tqdm import trange

from osgeo import ogr

from cudem import gdalfun
from cudem import utils
from cudem import waffles

_version = '0.0.11'

_usage = '''gdal_findrivers.py ({}): Locate possible rivers in a gdal file

usage: gdal_findrivers.py [ src_dem [dst_dem opts] ]

 Options:
  src_dem\t\tThe input DEM file-name
  dst_dem\t\tThe output DEM file-name

  --band\t\tThe input raster band (1)
  --help\t\tPrint the usage text
  --version\t\tPrint the version information

 Examples:
 % gdal_findrivers.py input_dem.tif output_dem.tif

CIRES DEM home page: <http://ciresgroups.colorado.edu/coastalDEM>
'''.format(_version)

#Circularity	Circularity( s Surface ) : Double precision	X	GEOS	computes the Circularity Index from the given Geometry by applying the following formula:
#index = ( 4 * PI * Sum(area) ) / ( Sum(perimeter) * Sum(perimeter) )
#it only applies to Polygons or MultiPolygons with the following interpretation:
#1.0 corresponds to a perfectly circular shape.
#very low values (near zero) correspond to a threadlike shape.
#intermediate values correspond to a more or less flattened shape; lower index values means a stronger flattening effect.
#if the given Geometry does not contains any Polygon but contains at least a Linestring the index will always assume a 0.0 value.
#if the given Geometry only contains one or more Points the index will always assume a NULL value.


def locate_no_data_zones(src_dem, dst_dem = None, band = 1, size_threshold = 100, max_circularity = .25, verbose = True):

    def expand_for(arr, shiftx=1, shifty=1):
        arr_b = arr.copy().astype(bool)
        for i in range(arr.shape[0]):
            for j in range(arr.shape[1]):
                if(arr[i,j]):
                    i_min, i_max = max(i-shifty, 0), min(i+shifty+1, arr.shape[0])
                    j_min, j_max = max(j-shiftx, 0), min(j+shiftx+1, arr.shape[1])
                    arr_b[i_min:i_max, j_min:j_max] = True
        return arr_b
    
    ## load src_dem array
    with gdalfun.gdal_datasource(src_dem, update=True if dst_dem is None else False) as src_ds:
        src_arr = src_ds.GetRasterBand(band).ReadAsArray()
        src_config = gdalfun.gdal_infos(src_ds)        
        src_arr[src_arr == src_config['ndv']] = np.nan

        ## generate the mask array
        msk_arr = np.zeros((src_config['ny'], src_config['nx']))
        msk_arr[np.isnan(src_arr)] = 1
        
        ## group adjacent non-zero cells
        l, n = scipy.ndimage.label(msk_arr)
        
        ## get the total number of cells in each group
        mn = scipy.ndimage.sum_labels(msk_arr, labels=l, index=np.arange(1, n+1))
        msk_arr[:] = 0
        for i in trange(0,
                        n,
                        desc='{}: Scanning data voids greater than {} cells'.format(
                            os.path.basename(sys.argv[0]), size_threshold
                        ),
                        leave=verbose):
            if mn[i] >= size_threshold:
                i += 1                
                ll = expand_for(l==i)
                #flat_value = np.nanpercentile(src_arr[ll], 5)
                #flat_value = np.nanmax(src_arr[ll]) - np.nanmin(src_arr[ll])
                flat_value = np.nanstd(src_arr[ll])
                msk_arr[l==i] = np.ceil(flat_value)#mn[i-1]
                #msk_arr[l==i] = mn[i-1]

        src_arr[np.isnan(src_arr)] = src_config['ndv']

        src_config['ndv'] = 0
        if dst_dem is None:
            src_ds.GetRasterBand(band).WriteArray(src_arr)
        else:
            gdalfun.gdal_write(msk_arr, dst_dem, src_config)
            with gdalfun.gdal_datasource(dst_dem) as nd_ds:
                dst_layer, ogr_format = gdalfun.ogr_polygonize(nd_ds)
                dst_vector = dst_layer + '.{}'.format(gdalfun.ogr_fext(ogr_format))
                
                utils.run_cmd('ogrinfo {vec} -sql "ALTER TABLE {lay} ADD COLUMN Area_ float"'.format(vec=dst_vector,lay=dst_layer))
                utils.run_cmd('ogrinfo {vec} -dialect SQLite -sql "UPDATE {lay} SET Area_ = ST_Area(geometry)"'.format(vec=dst_vector,lay=dst_layer))
                
                utils.run_cmd('ogrinfo {vec} -sql "ALTER TABLE {lay} ADD COLUMN Perim_ float"'.format(vec=dst_vector,lay=dst_layer))
                utils.run_cmd('ogrinfo {vec} -dialect SQLite -sql "UPDATE {lay} SET Perim_ = ST_Perimeter(geometry)"'.format(vec=dst_vector,lay=dst_layer))
                
                #ogrinfo mask1.shp -sql "ALTER TABLE mask1 ADD COLUMN Compact_ float"
                #ogrinfo mask1.shp -dialect SQLite -sql "UPDATE mask1 SET Compact_ = (4*3.14159*Area_)/(Perim_*Perim_)"
                
                utils.run_cmd('ogrinfo {vec} -sql "ALTER TABLE {lay} ADD COLUMN Circul_ float"'.format(vec=dst_vector,lay=dst_layer))
                utils.run_cmd('ogrinfo {vec} -dialect SQLite -sql "UPDATE {lay} SET Circul_ = Circularity(geometry)"'.format(vec=dst_vector,lay=dst_layer))

                d_ds = ogr.Open(dst_vector, 1)
                d_layer = d_ds.GetLayer()
                for feat in d_layer:
                    if feat.GetField('Circul_') > max_circularity:
                        d_layer.DeleteFeature(feat.GetFID())
                        
                d_ds = d_layer = None
            
    return(dst_dem if dst_dem is not None else src_dem, 0)

# Mainline
if __name__ == "__main__":

    src_dem = None
    dst_dem = None
    band = 1
    size_threshold = 1
    max_circularity = .25
    verbose = False
    
    i = 1
    while i < len(sys.argv):
        arg = sys.argv[i]
        if arg == '--size_threshold':
            size_threshold = utils.int_or(sys.argv[i+1], 1)
            i += 1
        elif arg == '--max_circularity':
            max_circularity = utils.float_or(sys.argv[i+1], 1)
            i += 1
        elif arg == '--band':
            band = utils.int_or(sys.argv[i+1], 1)
            i += 1
        elif arg == '-verbose':
            verbose = True
        elif src_dem is None:
            src_dem = arg
        elif dst_dem is None:
            dst_dem = arg
        else:
            sys.stderr.write(_usage)
            sys.exit(1)
        i = i + 1

    if src_dem is None:
        sys.stderr.write(_usage)
        sys.exit(1)

    #gdalfun.gdal_hydro_flatten(src_dem, dst_dem=dst_dem, band=band, size_threshold=size_threshold)
    locate_no_data_zones(src_dem, dst_dem=dst_dem, band=band, size_threshold=size_threshold, max_circularity=max_circularity)

#--END
