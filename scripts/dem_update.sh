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
radius=$(echo "$xinc * 3" | bc)
proj='epsg:4269'
max_diff=.25
min_weight=1
smoothing=5
mode='soest'
#
#  interpolate the new data through the old DEM
#
echo REGION is $region
echo INCREMENT is $xinc $yinc
echo RADIUS is $radius
echo MAXIMUM SCALED DIFF is$max_diff
echo MINIMUM WEIGHT is $min_weight
echo SMOOTHING FACTOR is $smoothing
echo UPDATING $2 with $1 using DIFFs

if [ "$mode" == "soest" ]; then
    #
    # grid the differences with surface
    #
    dlim $1 ${region}/-/-/${min_weight}/- -t_srs $proj | \
	gdal_query.py $2 -d_format "xyds" | \
	awk -v max_diff="$max_diff" '{if ($4 < max_diff) {print $1,$2,$3}}' | \
	gmt blockmedian $region -I$xinc/$yinc | \
	gmt surface $region -I$xinc/$yinc -G_diff.tif=gd+n-9999:GTiff -T.1 -Z1.2 -V -rp -C.5 -Lud -Lld -M${radius}

    dem_smooth.py _diff.tif -s $smoothing
    
    # # generate a boundary of the dem_diffs and apply a buffer
    # dlim $region/-/-/${min_weight}/- $1 | bounds --verbose -k $xinc/$(echo ${region} | awk -FR '{print $2}') -g > diffs.gmt
    # ogr2ogr -dialect SQLite -sql "select ST_Buffer(geometry, $radius) from diffs" diffs.shp diffs.gmt
    # gdal_clip.py _diff_smooth${smoothing}.tif diffs.shp
    # mv _diff_smooth${smoothing}_cut.tif _diff.tif
    
    #
    #  add the two grids
    #
    gdal_findreplace.py -s_value -9999 -t_value 0 _diff.tif _diff2.tif
    mv _diff2.tif _diff.tif

    #gdal_calc.py -A $2 -B _diff_smooth${smoothing}.tif --calc "A+B" --outfile $(basename $2 .tif)_update.tif
    gdal_calc.py -A $2 -B _diff.tif --calc "A+B" --outfile $(basename $2 .tif)_update.tif
    #
    #
    #  clean up the mess
    #
    rm _diff*
elif [ "$mode" == "waffles"]; then
    echo WAFFLES

    waffles -M num:mode=m $region/-/-/${min_weight}/- $1 -O _tmp -E ${xinc} -w
    waffles -M IDW:radius=${radius} $region/-/-/${min_weight}/- $1 -O _tmp -E ${xinc} -w -K 1000
    waffles -M $region surface:tension=1 -w -P $proj -O $(basename $2 .tif)_u -E ${xinc}/${yinc} _tmp.tif,200:weight_mask=_tmp_w.tif,1 $2,200,$min_weight

elif [ "$mode" == "ibcao" ]; then
    echo IBCAO REMOVE/RESTORE

    # grid new data by itself using mean num/nearneighbor:
    waffles $region/-/-/${min_weight}/- $1 -E $xinc/$yinc -O dem_update -M nearneighbor:sectors=12+m10:radius=$radius -w -P $proj -T 1:5:split_level=0 -K 2000

    # generate a boundary of the dem_diffs and apply a buffer
    # not needed since using nearneighbor
    #dlim dem_update.tif | bounds --verbose -k $xinc/$(echo ${region} | awk -FR '{print $2}') -g > diffs.gmt
    #dlim $region/-/-/${min_weight}/- $1 | bounds --verbose -k $xinc/$(echo ${region} | awk -FR '{print $2}') -g > diffs.gmt
    #ogr2ogr -dialect SQLite -sql "select ST_Buffer(geometry, $radius) from diffs" diffs.shp diffs.gmt
    #gdal_clip.py dem_update.tif diffs.shp

    # calculate the difference between 
    gdal_calc.py -A dem_update_cut.tif -B $2 --calc "B-A" --outfile _diff.tif --NoDataValue -9999.
        
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
else
    echo please choose an update method - soest, ibcao or waffles
fi
### End
