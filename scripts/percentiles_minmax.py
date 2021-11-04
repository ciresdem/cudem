#!/usr/bin/python
import numpy as np
from osgeo import gdal
import sys

raster=sys.argv[1]
ds = gdal.Open(raster)
myarray = np.array(ds.GetRasterBand(1).ReadAsArray())
myarray[myarray==0]=np.nan
#print myarray
myarray_flat=myarray.flatten()
#print myarray_pos
min_pct = int(np.nanpercentile(myarray_flat, 0.00001)-0.5)
max_pct = int(np.nanpercentile(myarray_flat, 99.9999)+0.5)

print("0.00001 percentile elevation is", min_pct)
print("99.9999 percentile elevation is", max_pct)
min_pct
max_pct

file = open("thresholds.csv", "w")
file.write(str(raster) + ','+ str(min_pct) + ',' + str(max_pct))
file.write('\n')
file.close() 

