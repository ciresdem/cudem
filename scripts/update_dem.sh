#!/bin/sh
#
# M. Love 12/21
#
#  Update a DEM with some new data
#
#  Adapted from https://topex.ucsd.edu/pub/srtm15_plus/update_grid
#
if [ ${#@} -lt 2 ]; then
   echo " Usage: update_dem.sh topo.xyz dem"
   echo "        topo.xyz   -  file of lon, lat, depth (neg meters); dlim supported"
   echo "        dem        -  existing dem to update"
   echo "Example: update_dem.sh german_iceland.xyz topo_old.tif"

   exit 1
fi
#
#  get the region and increments from the input DEM
#
region=$(dlim -r $2)
xinc=$(gdalinfo $2 | grep Pixel | awk -F= '{print $2}' | awk -F, '{print $1}' | sed 's/ (/''/g')
yinc=$(gdalinfo $2 | grep Pixel | awk -F= '{print $2}' | awk -F, '{print $2}' | sed 's/)/''/g')
yinc=$(echo "$yinc * -1" | bc)
proj='epsg:4269'
#
#  interpolate the new data through the old DEM
#
dlim $1 $region -t_srs $proj | gdal_query.py $2 -d_format "xyd" | gmt blockmedian $region -I$xinc/$yinc > diffs.xyz
#
#  make a difference grid
#
gmt blockmedian diffs.xyz $region -I$xinc/$yinc -rp -V | gmt surface $region -I$xinc/$yinc -Gdiff.tif=gd+n-9999:GTiff -T1. -Z1.2 -C.5 -V -rp -M2c
#
# nan/nodata to 0
#
gdal_findreplace.py -s_value -9999 -t_value 0 diff.tif diff2.tif
#
#  add the two grids
#
gdal_calc.py -A $2 -B diff2.tif --calc "A+B" --outfile $(basename $2 .tif)_update.tif
#
#  clean up the mess
#
rm diffs.xyz* diff*.tif*
