#!/bin/sh
## ------------------------------------------------------------
### coastline2xyz.sh  
## Copyright (c) 2019 - 2020 Matthew Love <matthew.love@colorado.edu>
## This file is liscensed under the GPL v.2 or later and 
## is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details. 
## <http://www.gnu.org/licenses/> 
##--------------------------------------------------------------
##
### Code:

c2x_version="0.0.1"
usage="usage: castline2xyz.sh [options]
\n\n
Options:\n
-I\tinput coastline shapefile\n
-Z\toutput z value (default is 0)\n
\n
coastline2xyz.sh v.${c2x_version}\n"

## Getopts
while getopts ":I:Z:W:E:S:N:O:" options; do
    case $options in
	I ) in_coast=$OPTARG;;
	Z ) out_zed=$OPTARG;;
	W ) w=$OPTARG;;
	E ) e=$OPTARG;;
	S ) s=$OPTARG;;
	N ) n=$OPTARG;;
	O ) out_coast=$OPTARG;;
	h ) echo -e $usage;;
	\? ) echo -e $usage
	exit 1;;
	* ) echo -e $usage
	    exit 1;;
    esac
done

if [ ! $in_coast ] ; then
    echo -e $usage
    exit 1
fi

if [ ! $out_coast ] ; then
    out_coast=$(basename $in_coast .shp).xyz
fi

#mkdir .c2x_tmp
cp $(dirname $in_coast)/$(basename $in_coast .shp).shp ${out_coast}.shp
cp $(dirname $in_coast)/$(basename $in_coast .shp).shx ${out_coast}.shx
#cd .c2x_tmp
tmp_coast=${out_coast}.csv
if [ ! $w ] ; then
    ogr2ogr $tmp_coast $(basename $in_coast .shp).shp -f CSV -lco GEOMETRY=AS_WKT -overwrite
else
    ogr2ogr $tmp_coast $in_coast -f CSV -lco GEOMETRY=AS_WKT -overwrite -clipsrc $w $n $e $s
    #ogr2ogr $tmp_coast $in_coast -f CSV -lco GEOMETRY=AS_WKT -overwrite -clipdst $w $s $e $n -clipsrc $w $s $e $n -spat $w $s $e $n 
fi
if [ ! $out_zed ] ; then
    cat $tmp_coast | sed -e '1,1d' | tr ',' '\n' | sed 's/[A-Za-z"()]*//g' | tr ' ' ',' | sed 's/^,//' | awk -F, '{if (NR != 1) {print $1,$2,$3}}' > $out_coast
else
    if [ -f "$out_zed" ]; then
	cat $tmp_coast | sed -e '1,1d' | tr ',' '\n' | sed 's/[A-Za-z"()]*//g' | tr ' ' ',' | sed 's/^,//' | awk -F, '{if (NR != 1 && $1 != 0) {print $1,$2}}' | gdal_query.py $out_zed -d_format 'xyg' -s_format '0,1,-' | grep -v "\-9999" | awk '{print $1,$2,$3-.5}' >> $out_coast
    else
	cat $tmp_coast | sed -e '1,1d' | tr ',' '\n' | sed 's/[A-Za-z"()]*//g' | tr ' ' ',' | sed 's/^,//' | awk -v zed="$out_zed" -F, '{if (NR != 1 && $1 != 0) {print $1,$2,zed}}' > $out_coast
    fi
fi
rm -rf $tmp_coast
cp $(dirname $in_coast)/$(basename $in_coast .shp).shx ${out_coast}.shp
cp $(dirname $in_coast)/$(basename $in_coast .shp).shx ${out_coast}.shx
#cd ..
#rm -rf .c2x_tmp
### End
