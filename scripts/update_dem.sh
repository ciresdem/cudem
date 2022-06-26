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
max_diff=100.25
min_weight=1
#mode="diffs"
#mode="ibcao"
mode="patch"
#
#  interpolate the new data through the old DEM
#
echo $region
echo $xinc
echo $yinc
echo $radius
echo $max_diff
echo $min_weight


if [ "$mode" == "diffs" ]; then
    echo UPDATING using DIFFs
    #
    # grid the differences with surface
    #
    dlim $1 ${region}/-/-/${min_weight}/- -t_srs $proj | \

	gdal_query.py $2 -d_format "xyds" | \
	awk -v max_diff="$max_diff" '{if ($4 < max_diff) {print $1,$2,$3}}' | \
	gmt blockmedian $region -I$xinc/$yinc | \
	gmt surface $region -I$xinc/$yinc -G_diff.tif=gd+n-9999:GTiff -T.1 -Z1.2 -V -rp -C.5 -Lud -Lld -M${radius}    
    #
    #  add the two grids
    #
    gdal_findreplace.py -s_value -9999 -t_value 0 _diff.tif _diff2.tif
    mv _diff2.tif _diff.tif
    dem_smooth.py _diff.tif -s 5
    gdal_calc.py -A $2 -B _diff_smooth5.tif --calc "A+B" --outfile $(basename $2 .tif)_update.tif
    #
    #
    #  clean up the mess
    #
    rm _diff*

elif [ "$mode" == "patch" ]; then
    echo UPDATING using PATCH
    ## grid the difference to array using query_dump / num
    dlim $1 ${region}/-/-/${min_weight}/- -t_srs $proj | \
	gdal_query.py $2 -d_format "xyds" | \
	awk -v max_diff="$max_diff" '{if ($4 < max_diff) {print $1,$2,$3}}' | \
	gmt blockmedian $region -I$xinc/$yinc | \
	gmt xyz2grd $region -I$xinc/$yinc -G_diff.tif=gd+n-9999:GTiff -V -rp
    
    ## polygonize the differences and add small buffer (1% or so)
    gdal_polygonize.py _diff.tif _diff_ply.shp
    ogr2ogr -dialect SQLite -sql "select ST_Buffer(geometry, $radius) from _diff_ply" _diff_ply_buff.shp _diff_ply.shp
    
    ## make zero array, inverse clipped to buffered polygonized diffs
    gdal_null.py -copy $2 _zeros.tif
    gdal_clip.py _zeros.tif _diff_ply_buff.shp
    #gdal_clip.py $2 _diff_ply_buff.shp-i

    waffles $region _diff.tif,200,2 _zeros_cut.tif,200,1 -E $xinc/$yinc -O _update -M stacks:supercede=True -w -P $proj

    ## surface is nicer output, but slower than gdal_fillnodata
    ## surface zero array and diffs...
    gdal_fillnodata.py _update.tif _update_full.tif
    #waffles $region _update.tif -E $xinc/$yinc -O _update_full -M surface:tension=1 -w -P $proj 
    #waffles $region _diff.tif,200,2 _zeros_cut.tif,200,1 -E $xinc/$yinc -O _update -M surface:tension=1 -w -P $proj 

    # dlim $1 ${region} | \
    # 	gmt blockmedian $region -I$xinc/$yinc | \
    # 	gmt surface $region -I$xinc/$yinc -G_update.tif=gd+n-9999:GTiff -T.1 -Z1.2 -V -rp -C.5 -Lud -Lld
    
    ## add surfaced diffs to self.dem
    gdal_calc.py -A $2 -B _update_full.tif --calc "A+B" --outfile $(basename $2 .tif)_update.tif
    
else
    echo IBCAO REMOVE/RESTORE

    # grid new data by itself using mean num/nearneighbor:
    #waffles $region $1 -E $xinc/$yinc -O dem_update -M num:mode=m:min_count=4 -w -P $proj
    #waffles $region/-/-/${min_weight}/- $1 -E $xinc/$yinc -O dem_update -M nearneighbor:sectors=12+m10:radius=$radius -w -P $proj -T 1:5:split_level=0
    waffles $region/-/-/${min_weight}/- $1,-1,2 $2,200,1 -E $xinc/$yinc -O dem_update -M surface:tension=1 -w -P $proj -T 1:5:split_level=0

    # generate a boundary of the dem_diffs and apply a buffer
    #dlim dem_update.tif | bounds --verbose -k $xinc/$(echo ${region} | awk -FR '{print $2}') -g > diffs.gmt
    #dlim $region/-/-/${min_weight}/- $1 | bounds --verbose -k $xinc/$(echo ${region} | awk -FR '{print $2}') -g > diffs.gmt
    #ogr2ogr -dialect SQLite -sql "select ST_Buffer(geometry, $radius) from diffs" diffs.shp diffs.gmt
    #gdal_clip.py dem_update.tif diffs.shp

    gdal_calc.py -A dem_update_cut.tif -B $2 --calc "B-A" --outfile _diff.tif --NoDataValue -9999.
    
    # dlim dem_update.tif ${region}/-/-/${min_weight}/- -t_srs $proj --weights | \
    # 	gdal_query.py $2 -d_format "xyds" | \
    # 	awk -v max_diff="$max_diff" '{if ($4 < max_diff) {print $1,$2,$3}}' | \
    # 	gmt blockmedian $region -I$xinc/$yinc -W | \
    # 	gmt surface $region -I$xinc/$yinc -G_diff.tif=gd+n-9999:GTiff -T.1 -Z1.2 -V -rp -C.5 -Lud -Lld -M${radius}
    # #gmt nearneighbor $region -I$xinc/$yinc -Gdem_update.tif=gd+n-9999:GTiff -W -V -rp -N10+m7 -S${radius}
    
    # calculate the difference between the orig DEM and the nn grid with the data:
    #gmt grdmath -N $2 dem_update.tif SUB = dem_diffs.tif=gd:GTiff
    
    # grid the diffs with surface and clip to buffer
    #waffles $region -E$xinc/$yinc dem_diffs.tif -O diffs_buff -M surface -C diffs.shp:invert=True
    #waffles $region -E$xinc/$yinc dem_diffs.tif -O diffs_buff -M surface:max_radius=${radius} -P $proj
    
    # apply the diffs_buff grid to the original DEM
    #gmt grdmath $2 diffs_buff.tif 0 AND SUB = dem_update.tif=gd:GTiff

    #
    #  add the two grids
    #
    gdal_findreplace.py -s_value -9999 -t_value 0 _diff.tif _diff2.tif
    mv _diff2.tif _diff.tif
    gdal_calc.py -A $2 -B _diff.tif --calc "A+B" --outfile $(basename $2 .tif)_update.tif
    #smooth_dem_bathy.py $(basename $2 .tif)_update.tif -s 5
    #
    #  clean up the mess
    #
    rm _diff.*
    
fi
