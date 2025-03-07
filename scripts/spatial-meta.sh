#!/bin/sh
### spatial-meta.sh
##
## Copyright (c) 2019 - 2021 CIRES Coastal DEM Team
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
##
## Generate a spatial metadata shapefile for each region found in `itiles`
## which is a GMT fortmatted vector file containing regional tile(s).
## `idatalist` is a datalist that holds all the data for which the spatial
## metadata will be generated. 
##
## Uses `regions, `bounds, `gmt, `gdal
##
### Code:

usage="spatial-meta.sh [ -I all_data.datalist -R all_tiles.gmt [ -EFOPX ]  ]\n\
Generate spatial metadata for a DEM tileset.\n\n\
where 'all_data.datalist' is a MBSystem/dlim style datalist containing elevation data\n\
and 'all_tiles.gmt' is a GMT formatted vector containing regional tile(s).\n\
  note: the tile vector can be any ogr compatible format (such as shp)
\n\
Output will be a GMT formatted vector and a ESRI shapefile for each tile.\n\
\n\
Options:
  -I - input datalist
  -R - input tileset vector
  -O - output basename (meta)
  -X - extend region by number of cells (6)
  -E - cellsize of gridded boundary (0.000277777777)
  -F - output vector format (GeoJSON)
  -P - EPSG code of output vector (4326)
  -v - make output shapefile valid
 Examples:
 %% spatial-meta.sh -I master.datalist -R master_regions.gmt -O meta -X 6 -P 4326 \n\
\n\
CIRES DEM home page: <http://ciresgroups.colorado.edu/coastalDEM>
"

while getopts ":I:R:X:O:E:F:P:" options; do
    case $options in
	I ) idatalist=$OPTARG;;
	R ) itiles=$OPTARG;;
	X ) extend=$OPTARG;;
	O ) oname=$OPTARG;;
	E ) iinc=$OPTARG;;
	P ) epsg=$OPTARG;;
	F ) frmt=$OPTARG;;
	v ) makevalid=True;;
	h ) echo -e $usage;;
	\? ) echo -e $usage
	exit 1;;
	* ) echo -e $usage
	    exit 1;;
    esac
done

if [ -z "$idatalist" -o -z "$itiles" ]; then
    printf "$usage"
    exit
fi

## Set the Increment
## 1/9 0.000030864194
## 1/3 - .0000925925925
## 1 - 0.000277777777
if [ -z "$iinc" ]; then
    iinc=0.000277777777
fi

#if [ $iinc -lt .0000925 ]; then
#    iinc=.0000925925925
#fi

## Set the cell extend
if [ -z "$iextend" ]; then
    iextend=6
fi

## Set the basename
if [ -z "$oname" ]; then
    oname="meta"
fi

## Set the output EPSG
if [ -z "$epsg" ]; then
    epsg=4326
fi

## Set the output Format
if [ -z "$frmt" ]; then
    frmt="GeoJSON"
fi

tcount=$(gmt gmtinfo $itiles -F | awk '{print $2}')
printf "\nspatial-meta: Generating Spatial Metadata for %d regions using %s\n" $tcount $idatalist

## Loop through each region found in the $itiles gmt vector
for region in $(gmt gmtinfo $itiles -I- -As); do 

    # The proc region is the region extended by 6 cells at 1/9 (form 1/9), 
    # 6 cells at 1/3 (for 1/3) and 2 cells at 1/3 (for 1/9), based on the increment value
    proc_region=$(regions -b $(echo $iextend*$iinc | bc) $region -e)
    gmt_region=$(regions $proc_region -ee)
    out_name="${oname}_"$(regions $region -n)_$(date +%Y)v1_sm
    printf "\nspatial-meta: %s < %s >\n" $out_name $proc_region

    if [ "$frmt" = "GMT" ]; then
	# Initialize the output gmt vector
	cat /dev/null | bounds -g > ${out_name}.gmt
	sed -i 's/\@NName/\@NName|Title|Agency|Date|Type|Resolution|HDatum|VDatum|URL/g' ${out_name}.gmt
	sed -i 's/\@Tstring/\@Tstring|string|string|string|string|string|string|string|string/g' ${out_name}.gmt
	sed -i "s,\@GMULTIPOLYGON,\@GMULTIPOLYGON\ \@R${gmt_region},g" ${out_name}.gmt
	dlim ${idatalist} -c ${proc_region} > tmp.csv
    elif [ "$frmt" = "GeoJSON" ]; then
	printf "{ \"type\": \"FeatureCollection\",\n\"features\": [\n" > ${out_name}.json
    fi
    
    # Loop through each datalist entry
    cnt=0
    for datalist in $(dlim ${idatalist} -d ${proc_region} | awk '{print $1}') ; do
	dlbn=$(basename $datalist .datalist)
	
	# Generate the boundary of the data found in this datalist and add it to the output gmt vector
	printf "spatial-meta: using datalist: %s < %s >\n" $dlbn $datalist

	if [ "$frmt" = "GeoJSON" ]; then
	    meta=$(grep "${dlbn}" $(basename $idatalist .datalist).json)
	    
	    if [ $cnt -gt 0 ]; then
		printf "," >> ${out_name}.json
	    fi
	    printf "{ \"type\": \"Feature\", \"properties\": ${meta}," >> ${out_name}.json
	    printf " \"geometry\": { \"type\": \"MultiPolygon\",\n \"coordinates\": [[[" >> ${out_name}.json
	    
	    dlim $datalist $proc_region | bounds -k ${iinc}/$(regions $proc_region -ee) -jjj --verbose >> ${out_name}.json
	    printf "]]]}}" >> ${out_name}.json
	elif [ "$frmt" = "GMT" ]; then
	    meta=$(grep "${dlbn}" tmp.csv)
	    
	    dlim $datalist $proc_region | bounds -k ${iinc}/$(regions $proc_region -ee) -n $dlbn -gg --verbose | \
		sed "0,/\@D${dlbn}/{s^\@D${dlbn}^\@D${meta}^;}" >> ${out_name}.gmt;
	fi
	cnt=$(echo $cnt+1 | bc)
    done

    if [ "$frmt" = "GeoJSON" ]; then
	printf "]}\n" >> ${out_name}.json
    fi
    
    # # Convert gmt vector to shapefile
    # printf "spatial-meta: converting boundary to shapefile: %s.shp\n" $out_name
    # if [ $makevalid ] ; then
    # 	ogr2ogr ${out_name}.shp ${out_name}.gmt -overwrite -a_srs EPSG:${epsg} -makevalid
    # else
    # 	ogr2ogr ${out_name}.shp ${out_name}.gmt -overwrite -a_srs EPSG:${epsg}
    # fi
done

### End
