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
#radius=$(echo "$xinc * 6" | bc)
radius='3s'
min_weight=.75
#waffles $1 $2,200,${min_weight} $region/-/-/${min_weight}/- -t_srs $proj -E $xinc/$yinc -w -M surface:tension=1 -O $(basename $2 .tif)_update
#waffles $1 $2,200,${min_weight} $region/-/-/${min_weight}/- -t_srs $proj -E $xinc/$yinc -w -M invdst:radius1=$radius:radius2=$radius -O $(basename $2 .tif)_update 
#
#  grid incoming data to DEM res
#
#waffles $1 $region -t_srs $proj -E $xinc/$yinc -w -M IDW:radius=$radius:min_points=24 -O new
#waffles $1 $region -t_srs $proj -E $xinc/$yinc -w -M surface:max_radius=2c:tension=1:convergence=.5 -O new
#waffles $1 $region/-/-/${min_weight}/- -t_srs $proj -E $xinc/$yinc -w -M nearest:radius1=$radius:radius2=$radius:block=True -O new
#waffles $1 $region/-/-/${min_weight}/- -t_srs $proj -E $xinc/$yinc -w -M num:mode=m -O new
#
#  interpolate the new data through the old DEM
#
#gdal_calc.py -A new.tif -B $2 --calc "A-B" --outfile new_diff.tif --NoDataValue -9999
dlim $1 $region/-/-/${min_weight}/- -t_srs $proj | gdal_query.py $2 -d_format "xyd" | gmt blockmedian $region -I$xinc/$yinc > diffs.xyz
#dlim new.tif $region -t_srs $proj | gdal_query.py $2 -d_format "xyd" | gmt blockmedian $region -I$xinc/$yinc > diffs.xyz
#
#  make a difference grid
#
gmt blockmedian diffs.xyz $region -I$xinc/$yinc -rp -V | gmt surface $region -I$xinc/$yinc -Gdiff.tif=gd+n-9999:GTiff -T1. -Z1.2 -C.5 -V -rp -M2c
#
# nan/nodata to 0
#
#mv diff.tif diff2.tif
gdal_findreplace.py -s_value -9999 -t_value 0 diff.tif diff2.tif
#gdal_findreplace.py -s_value -9999 -t_value 0 new_diff.tif new_diff2.tif
#
#  add the two grids
#
gdal_calc.py -A $2 -B diff2.tif --calc "A+B" --outfile $(basename $2 .tif)_update.tif
#gdal_calc.py -A $2 -B new_diff2.tif --calc "A+B" --outfile $(basename $2 .tif)_update.tif
#
#  clean up the mess
#
#rm diffs.xyz* diff*.tif*
