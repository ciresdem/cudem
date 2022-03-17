#!/bin/sh

## testing xyz data cleaning by comparing it to a cleaned DEM

in_xyz=$1
in_dem=$2
thresh=$3
block=$4

if [ ! $thresh ]; then
    thresh=.25
fi

if [ ! $block ]; then
    block=3s
fi

if [ "$in_dem" == "gmrt" ]; then
    ## fetch the DEM (currently from GMRT)
    in_dem="gmrt:layer=topo-mask"
    fetches $(regions3 $(dlim -r $in_xyz) -B .25 -e) ${in_dem}
    gmrt_name_region=$(regions3 $(dlim -r $in_xyz) -B .25 --name)
    q_dem=$(ls gmrt/gmrt_topo-mask_${gmrt_name_region}.tif)
elif [ "$in_dem" == "block" ]; then
    ## grid the data with gmt xyz2grd -Am @ $block cell-size
    waffles $(dlim -r $in_xyz) -E $block -O _tmp_block -M num:mode=Am $1
    q_dem=_tmp_block.tif
else
    q_dem=$in_dem
fi

## query the DEM with the input points and return
## in format of xyz(diff)(percentage diff)(scaled diff)
dlim $in_xyz | gdal_query.py $q_dem -d_format xyzgdcs > ${in_xyz}.xyd

## remove points where the 'scaled difference' is above .2
awk -v thresh="$thresh" '{if ($4 == "nan" || $7 < thresh) {print $1,$2,$3}}' ${in_xyz}.xyd > $(basename ${in_xyz} .xyz)_clean.xyz

rm -rf ${in_xyz}.xyd
if [ "$in_dem" == "gmrt" ]; then
    rm -rf gmrt
elif [ "$in_dem" == "block" ]; then
    rm -rf $q_dem
fi

### End
