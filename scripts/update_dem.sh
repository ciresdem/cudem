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
echo $region
echo $xinc
echo $yinc
echo $radius
echo $max_diff
echo $min_diff

## grid new data with nearneighbor
#waffles $region/-/-/.75/- $1 -E $xinc/$yinc -O dem_update -M nearneighbor:radius=${radius}:sectors=4+m4 -w -P $proj
waffles $region/-/-/.75/- $1 $2,200,1-e10 -E $xinc/$yinc -O dem_update -M surface -w -P $proj

## crop the grid
#gdal_crop.py dem_update.tif
#gmt grdmath dem_update.tif $2 SUB = _diff.tif=gd:GTiff
##gdal_calc.py -A dem_update.tif -B $2 --calc "B-A" --outfile _diff.tif --NoDataValue -9999.

## grid the differences with surface
#dlim dem_update_crop.tif ${region} -t_srs $proj | gdal_query.py $2 -d_format "xyd" | \
#    gmt blockmedian $region -I$xinc/$yinc | \
#    gmt surface $region -I$xinc/$yinc -G_diff.tif=gd+n-9999:GTiff -T.1 -Z1.2 -V -rp -C.5 -Lu${max_diff} -Ll${min_diff} -M${radius}

##gdal_findreplace.py -s_value -9999 -t_value 0 _diff.tif _diff2.tif
##mv _diff2.tif _diff.tif
##gdal_calc.py -A $2 -B _diff.tif --calc "A+B" --outfile $(basename $2 .tif)_update.tif
#gmt grdmath
# # grid the new data
# waffles $region/-/-/.75/- -P $proj -E $xinc/$yinc -w -O _updates -M nearneighbor:radius=${radius}:sectors=4+m4 $1
# # calculate diffs
# gdal_calc.py -A _updates.tif -B $2 --calc "B-A" --outfile _diff.tif --NoDataValue 0
# #gdal_findreplace.py -s_value -9999 -t_value 0 _diff.tif _diff2.tif
# gdal_calc.py -A $2 -B _diff2.tif --calc "A+B" --outfile $(basename $2 .tif)_update.tif

# dlim $1 ${region}/-/-/.75/- -t_srs $proj | gdal_query.py $2 -d_format "xyd" | \
#     gmt blockmedian $region -I$xinc/$yinc | \
#     gmt surface $region -I$xinc/$yinc -G_diff.tif=gd+n-9999:GTiff -T1. -Z1.2 -V -rp -C.5 -Lu${max_diff} -Ll${min_diff} -M${radius}
# #
# #  make a difference grid
# #
# #waffles _diff.xyz $region -E $xinc/$yinc -w -M IDW:radius=$radius:min_points=12 -O _diff -P$proj
# #gmt blockmedian _diff.xyz $region -I$xinc/$yinc -rp -V | gmt surface $region -I$xinc/$yinc -G_diff.tif=gd+n-9999:GTiff -T1. -Z1.2 -V -rp -C.5 #-M12c
# #gmt surface _diff.xyz $region -I$xinc/$yinc -G_diff.tif=gd+n-9999:GTiff -T1. -Z1.2 -V -rp -C.5 -Lu${max_diff} -Ll${min_diff} #-M12c
# #
# #  add the two grids
# #
# gdal_findreplace.py -s_value -9999 -t_value 0 _diff.tif _diff2.tif
# mv _diff2.tif _diff.tif
# gdal_calc.py -A $2 -B _diff.tif --calc "A+B" --outfile $(basename $2 .tif)_update.tif
# #
# #  clean up the mess
# #
# rm _diff.*

# ## IBCAO REMOVE/RESTORE

# # grid new data by itself using nearneighbor:
# #waffles $region $1 -E $xinc/$yinc -O dem_update -M nearneighbor:radius=${radius} -w -P $proj

# # calculate the difference between the orig DEM and the nn grid with the data:
# #gmt grdmath $2 dem_update.tif SUB = dem_diffs.tif=gd:GTiff

# # generate a boundary of the dem_diffs and apply a buffer
# dlim dem_diffs.tif | bounds --verbose -k $xinc/$(echo ${region} | awk -FR '{print $2}') -g > diffs.gmt
# ogr2ogr -dialect SQLite -sql "select ST_Buffer(geometry, $radius) from diffs" diffs.shp diffs.gmt

# # grid the diffs with surface and clip to buffer
# waffles $region -E$xinc/$yinc dem_diffs.tif -O diffs_buff -M surface -C diffs.shp:invert=True
# #waffles $region -E$xinc/$yinc dem_diffs.tif -O diffs_buff -M surface:max_radius=${radius} -P $proj

# # apply the diffs_buff grid to the original DEM
# gmt grdmath $2 diffs_buff.tif 0 AND SUB = dem_update.tif=gd:GTiff
