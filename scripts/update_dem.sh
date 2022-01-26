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
radius=$(echo "$xinc * 9" | bc)
proj='epsg:4269'
max_diff=1000
min_diff=-1000
#
#  interpolate the new data through the old DEM
#
dlim $1 $region -t_srs $proj | \
    gdal_query.py $2 -d_format "xyd" | \
    awk '{if ($3 < "$max_diff" && $3 > "$min_diff") {print $1,$2,$3}}' | \
    gmt blockmedian $region -I$xinc/$yinc | \ 
    gmt surface $region -I$xinc/$yinc -G_diff.tif=gd+n-9999:GTiff -T1. -Z1.2 -V -rp -C.5 -Lu${max_diff} -Ll${min_diff} -M${radius}
#
#  make a difference grid
#
#waffles _diff.xyz $region -E $xinc/$yinc -w -M IDW:radius=$radius:min_points=12 -O _diff -P$proj
#gmt blockmedian _diff.xyz $region -I$xinc/$yinc -rp -V | gmt surface $region -I$xinc/$yinc -G_diff.tif=gd+n-9999:GTiff -T1. -Z1.2 -V -rp -C.5 #-M12c
#gmt surface _diff.xyz $region -I$xinc/$yinc -G_diff.tif=gd+n-9999:GTiff -T1. -Z1.2 -V -rp -C.5 -Lu${max_diff} -Ll${min_diff} #-M12c
#
#  add the two grids
#
gdal_findreplace.py -s_value -9999 -t_value 0 _diff.tif _diff2.tif
mv _diff2.tif _diff.tif
gdal_calc.py -A $2 -B _diff.tif --calc "A+B" --outfile $(basename $2 .tif)_update.tif
#
#  clean up the mess
#
rm _diff.*
