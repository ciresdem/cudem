#!/bin/sh
#
# M. Love 12/21
#
#  Update a DEM with some new data
#
#  Adapted from https://topex.ucsd.edu/pub/srtm15_plus/update_grid
#
if [ ${#@} -lt 2 ]; then
   echo " Usage: dem_patch.sh topo.xyz dem"
   echo "        topo.xyz   -  file of lon, lat, depth (neg meters); dlim supported"
   echo "        dem        -  existing dem to update"
   echo "Example: dem_patch.sh german_iceland.xyz topo_old.tif"

   exit 1
fi
#
#  get the region and increments from the input DEM
#
region=$(dlim -r $2)
xinc=$(gdalinfo $2 | grep Pixel | awk -F= '{print $2}' | awk -F, '{print $1}' | sed 's/ (/''/g')
yinc=$(gdalinfo $2 | grep Pixel | awk -F= '{print $2}' | awk -F, '{print $2}' | sed 's/)/''/g')
yinc=$(echo "$yinc * -1" | bc)
radius=$(echo "$xinc * 24" | bc)
#radius2=$(echo "$xinc * 15" | bc)
proj='epsg:4269'
max_diff=100.25
min_weight=1

echo $region
echo $xinc
echo $yinc
echo $radius
echo $max_diff
echo $min_weight
echo PATCHING DEM $2
#
# grid the difference to array using query_dump / num
#
#dlim $1 ${region}/-/-/${min_weight}/- -t_srs $proj | \
dlim $1 ${region} | \
    gdal_query.py $2 -d_format "xyds" | \
    awk -v max_diff="$max_diff" '{if ($4 < max_diff) {print $1,$2,$3}}' | \
    gmt blockmedian $region -I$xinc/$yinc | \
    gmt xyz2grd $region -I$xinc/$yinc -G_diff.tif=gd+n-9999:GTiff -V -rp
#
# polygonize the differences and add small buffer (1% or so)
#
gdal_polygonize.py _diff.tif _diff_ply.shp
ogr2ogr -dialect SQLite -sql "select ST_Buffer(geometry, $radius) from _diff_ply" _diff_ply_buff.shp _diff_ply.shp
#
# make zero array, clipped to buffered polygonized diffs
#
gdal_zeros.py -copy $2 _zeros.tif
gdal_clip.py _zeros.tif _diff_ply_buff.shp
#
# stack clipped zero grid and diff grid
#
waffles $region _diff.tif,200,2 _zeros_cut.tif,200,1 -E $xinc/$yinc -O _update -M stacks:supercede=True -w -P $proj
#
# fill nodata in stacked grid
#
gdal_fillnodata.py _update.tif _update_full.tif -si 2
#waffles $region _update.tif -O _update_full -E $xinc/$yinc -M surface:tension=1
#    
# add full diffs to dem
#
gdal_calc.py -A $2 -B _update_full.tif --calc "A+B" --outfile $(basename $2 .tif)_patch.tif
#
# cleanup
#
rm -rf _diff* _update* _zeros*
#
# End
#
