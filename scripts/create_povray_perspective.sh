#!/bin/sh
### create_povray_perspective.sh
##
## Copyright (c) 2011 - 2023 CIRES Coastal DEM Team
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
### Code:  

_version="0.0.1"
_usage="Usage: create_povray_perspective.sh [options]

\n\n\
Options:\n\
-I\tInput gdal file\n\
-O\tOutput Hillshade image file\n\
-C\tThe cpt color palette file to use, if not given\n\
\twill attempt to generate one.\n\
\n\
create_povray_perspective v.${_version}
"

## ==============================================
## Getopts
## ==============================================
while getopts ":I:O:Z:P:C:A:T:i" options; do
    case $options in
	I ) ingrd=$OPTARG;;
	O ) outgrd=$OPTARG;;
	Z ) zex=$OPTARG;;
	P ) proj=$OPTARG;;
	C ) cpt=$OPTARG;;
	A ) aziumth=$OPTARG;;
	T ) altitude=$OPTARG;;
	i ) invert="True";;
	h ) echo -e $_usage;;
	\? ) echo -e $_usage
	exit 1;;
	* ) echo -e $_usage
	    exit 1;;
    esac
done

if [ ! $ingrd ] ; then
    echo -e $_usage
    exit 1
fi

if [ ! $outgrd ] ; then
    outgrd=${ingrd}_persp.tif
fi

if [ ! $cpt ] ; then
    gmin=$(gdalinfo $ingrd -stats | grep MINIMUM | cut -d= -f2)
    gmax=$(gdalinfo $ingrd -stats | grep MAXIMUM | cut -d= -f2)
    colortable.py --min $gmin --max $gmax > colors.txt
    cpt=colors.txt
fi

gdal_translate -ot UInt16 -of PNG -scale $(gdalinfo $ingrd -stats | grep Minimum | sed 's/,/=/g' | cut -d= -f2) $(gdalinfo $ingrd -stats | grep Minimum | sed 's/,/=/g' | cut -d= -f4) 0 65535 $ingrd temp.png
row=$(gdalinfo $ingrd | grep "Size is" | awk '{print $3}' | awk -F, '{print $1}')
col=$(gdalinfo $ingrd | grep "Size is" | awk '{print $4}')
gdal_translate -srcwin 1 1 $row $col -of PNG temp.png dem_16bit.png
gdaldem color-relief $ingrd $cpt temp2.tif
gdal_translate -srcwin 1 1 $row $col -of PNG temp2.tif rgb.png
rm temp.png temp2.tif
convert -size 10000 dem_16bit.png dem_16bit_10000.png
convert -size 10000 rgb.png rgb_10000.png

## pov-ray

#make the povray script:
create_povray_template.sh $ingrd > ${outgrd}.pov
povray ${outgrd}.pov +W1000 +H800
