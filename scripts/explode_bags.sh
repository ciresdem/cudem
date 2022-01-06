#!/bin/sh
function help () {
echo "explode_bags.sh - Explode all supergrids from a VR BAG"
echo "Usage: $0 bag"
}

if [ ${#@} == 1 ]; 
then
    bag=$1
    echo "Processing" $bag
    echo -- Checking to see if bag is VR and contains supergrids
    bag_check=`gdalinfo $bag | grep -e "HAS_SUPERGRIDS" | cut -d "=" -f2`
    echo "bag_check is" $bag_check
    if [[ -z "$bag_check" ]];
    then
	echo "BAG is NOT VALID VR"
	echo "Reprojecting to NAD83 and converting to tif"
	gdalwarp $bag -dstnodata -999999 -r bilinear -t_srs EPSG:4269 $(basename $bag .bag)".tif" -overwrite
    else
	echo "BAG is VALID VR"
	all_vr_bags=`gdalinfo $bag -oo MODE=LIST_SUPERGRIDS | grep -e "NAME" | cut -d "=" -f2`
	count=0
	for this_bag in $all_vr_bags;
	do
	    #gdal_translate $this_bag $(basename $bag .bag)"_"$(count)".tif"
	    echo "Reprojecting to NAD83 and converting to tif"
	    gdalwarp $this_bag -dstnodata -999999 -r bilinear -t_srs EPSG:4269 $(basename $bag .bag)"_"${count}".tif" -overwrite
	    count=$((count + 1))
	done
    fi
else
    help
fi
