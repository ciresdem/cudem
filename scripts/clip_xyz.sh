#!/bin/sh
### clip_xyz.sh
##
## Copyright (c) 2013 - 2021 CIRES Coastal DEM Team
##
## Permission is hereby granted, free of charge, to any person obtaining a copy 
## of this software and associated documentation files (the "Software"), to deal 
## in the Software without restriction, including without limitation the rights 
## to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies 
## of the Software, and to permit persons to whom the Software is furnished to do so, 
## subject to the following conditions:
##
## The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
##
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
## INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR 
## PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE 
## FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
## ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
##
### Commentary: 
## clip xyz data to ogr
##
### Code:

usage="clip_xyz.sh - clip xyz data to ogr\n\n\
usage: clip_xyz.sh [ -DIMOSiVQ [ args ] ]
\n\n\
Clip an xyz data file according to the given OGR compatible polygon file.\n\
The default will return xyz data points which are inside of the given polygon, \n\
to return the data points which are outside the given polygon, use the -i switch to invert the clipping.\n\n\
 -I input xyz file\n\
 -M input polygon mask\n\
\n\
 -O output xyz file\n\
 -S cellsize mask\n\
 -D xyz delimiter \n\
 -i invert the mask\n\
 -V increase verbosity\n\
 -Q run quickly\n\
\n\
 Examples:\n\
 % clip_xyz.sh -I input.xyz -O input_clip.xyz -M coast_ply.shp -i\n\
\n\
Report bugs to <matthew.love@colorado.edu>\n\
CIRES DEM home page: <http://ciresgroups.colorado.edu/coastalDEM>
"
while getopts ":I:M:O:iS:D:VQ" options; do
    case $options in
	I ) inxyz=$OPTARG;;
	O ) outxyz=$OPTARG;;
	M ) maskpoly=$OPTARG;;
	D ) delim=$OPTARG;;
	S ) inc=$OPTARG;;
	i ) invert=True;;
	V ) verbose=True;;
	Q ) quickly=True;;
	h ) echo -e $usage;;
	\? ) echo -e $usage
	exit 1;;
	* ) echo -e $usage
	    exit 1;;
    esac
done

# Check for variables and set them to default 
# values if they aren't set by the user.
if [ ! $inxyz ] || [ ! $maskpoly ] ; then
    echo -e $usage
    exit 1
fi

if [ ! $outxyz ] ; then
    outxyz=${inxyz}_clip.xyz
fi

if [ ! $delim ] ; then
    delim=""
else
    delim="-delimiter ${delim}"
fi

if [ ! $inc ] ; then
    inc=0.00009259259
fi

msk=$(basename $inxyz .xyz)_msk.tif

if [ $verbose ] ; then
    echo $msk
fi

# Generate the null grid

if [ $verbose ] ; then
    gdal_null.py -region $(gmt gmtinfo -C $inxyz | awk '{print $1,$2,$3,$4}') -cell_size $inc -verbose -overwrite $msk
else
    gdal_null.py -region $(gmt gmtinfo -C $inxyz | awk '{print $1,$2,$3,$4}') -cell_size $inc -overwrite $msk
fi

# Burn the polygon onto the null grid - use -i to invert the polygon

if [ $invert ] ; then
    gdal_rasterize -i -burn 1 -l $(basename $maskpoly .shp) $maskpoly $msk
else
    gdal_rasterize -burn 1 -l $(basename $maskpoly .shp) $maskpoly $msk
fi

if [ $quickly ] ; then
    xyz_clip.py $delim -quickly $inxyz $msk $outxyz
else
    xyz_clip.py $delim $inxyz $msk $outxyz
fi

if [ $verbose ] ; then
    echo "Completed"
fi
## END
