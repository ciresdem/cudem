#!/bin/sh
### hillshade.sh
##
## Copyright (c) 2011 - 2021 CIRES Coastal DEM Team
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

_version="0.1.7"
_usage="Usage: hillshade.sh [options]

\n\n\
Options:\n\
-I\tInput gdal file\n\
-O\tOutput Hillshade image file\n\
-Z\tThe Z Exageration [1]\n\
-A\tThe azimuth to use [315]\n\
-T\tThe altitude to use [45]\n\
-P\tThe output tif projection ESPG code [4326]\n\
-C\tThe cpt color palette file to use, if not given\n\
\twill attempt to generate one.\n\
\n\
hillshade v.${_version}
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
    outgrd=${ingrd}_hs.tif
fi

if [ ! $zex ] ; then
    zex=1
fi

if [ ! $proj ] ; then
    proj=4326
fi

if [ ! $cpt ] ; then
    gmin=$(gdalinfo $ingrd -stats | grep MINIMUM | cut -d= -f2)
    gmax=$(gdalinfo $ingrd -stats | grep MAXIMUM | cut -d= -f2)
    colortable.py --min $gmin --max $gmax > colors.txt
    cpt=colors.txt
fi

if [ ! $azimuth ] ; then
    azimuth=315
fi

if [ ! $azltitude ] ; then
    altitude=45
fi

## ==============================================
# Generate the Hillshade
## ==============================================
gdaldem hillshade -s 111120 -z $zex -az $azimuth -alt $altitude $ingrd hillshade.tif
# Generate the color-releif
gdaldem color-relief $ingrd $cpt colors.tif
# Composite the hillshade and the color-releif
composite -compose multiply -depth 8 colors.tif hillshade.tif output.tif
mogrify -modulate 115 -depth 8 output.tif
# Generate the combined georeferenced tif
gdal_translate -co "TFW=YES" $ingrd temp.tif
mv temp.tfw output.tfw
gdal_translate -a_srs epsg:$proj output.tif temp2.tif
# Cleanup
rm output.tif temp.tif hillshade.tif colors.tif output.tfw
#subtract 2 cells from rows and columns 
#gdal_translate -srcwin 1 1 $(gdal_rowcol.py $ingrd t) temp2.tif $outgrd
gdal_translate temp2.tif ${outgrd}.tif
rm temp2.tif

### End
