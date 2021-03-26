#!/bin/sh
##
### Commentary:
##
## Generate a spatial metadata shapefile for each region found in `itiles`
## which is a GMT fortmatted vector file containing regional tile(s).
## `idatalist` is a datalist that holds all the data for which the spatial
## metadata will be generated. 
##
## It is assumed that a file `idatalist`.csv exists
## and that file will be used to fill the attribute data for each datalist.
##
## The $datalist.csv file is formatted like this:
##
##  Title,Agency,Date,Type,Resolution,HDatum,VDatum,URL,Unc
##  "NOAA CRM",NOAA,2009,DEM,0.00027777,WGS84,NAVD88,coast.noaa.gov,NA
## 
## Uses `regions, `bounds, `mbsystem, `gmt, `gdal
##
### Code:

usage="spatial-meta.sh [ all_data.datalist all_tiles.gmt ]\n\
Generate spatial metadata for a DEM tileset.\n\n\
where 'all_data.datalist' is a MBSystem style datalist containing xyz data\n\
and 'all_tiles.gmt' is a GMT formatted vector containing regional tile(s).\n\
\n\
Output will be a GMT formatted vector and a ESRI shapefile for each tile.\n\
\n\
 Examples:
 %% spatial-meta.sh master.datalist master_regions.gmt\n\
\n\
Report bugs to <matthew.love@colorado.edu>\n\
CIRES DEM home page: <http://ciresgroups.colorado.edu/coastalDEM>
"

idatalist=$1
itiles=$2

if [ -z "$idatalist" -o -z "$itiles" ]; then
    printf "$usage"
    exit
fi

## Set the Increment
iinc=.0000925925925
#iinc=0.000277777777
#iinc=0.000030864194

tcount=$(gmt gmtinfo $itiles -F | awk '{print $2}')
printf "\n: Generating Spatial Metadata for %d regions using %s\n" $tcount $idatalist

## Loop through each region found in the $itiles gmt vector
for region in $(gmt gmtinfo $itiles -I- -As); do 

    # The proc region is the region extended by 6 cells at 1/9 (form 1/9), 
    # 6 cells at 1/3 (for 1/3) and 2 cells at 1/3 (for 1/9), based on the increment value
    proc_region=$(regions -b $(echo 2*$iinc | bc) $region -e)
    out_name="meta_"$(regions $region -n)_$(date +%Y)
    printf "\n: %s < %s >\n" $out_name $proc_region

    # Initialize the output gmt vector
    cat /dev/null | bounds -g > ${out_name}.gmt
    sed -i 's/\@NName/\@NName|Title|Agency|Date|Type|Resolution|HDatum|VDatum|URL|UNC/g' ${out_name}.gmt
    sed -i 's/\@Tstring/\@Tstring|string|string|string|string|string|string|string|string|string/g' ${out_name}.gmt

    # Loop through each datalist entry
    for datalist in $(cat ${idatalist} | awk '{print $1}' | grep -v "\#") ; do 

	# Gather the data from the datalist entry that is within our region
	dlbn=$(basename $datalist .datalist)
	data=$(mbdatalist -I${datalist} $proc_region -F-1 | awk '{print $1}')
	dinfo="NA|NA|NA|NA|NA|NA|NA|NA|NA"

	# If there is data in the region, continue
	if [ -n "$data" ]; then

	    # Generate the boundary of the data found in this datalist and add it to the output gmt vector
	    printf "  datalist: %s < %s >\n" $dlbn $datalist
	    cat $data | bounds -k ${iinc}/$(regions $proc_region -ee) -n $dlbn -gg --verbose >> ${out_name}.gmt;

	    # Add extra attribute data (from $datalist.csv).
	    if [ -f ${datalist}.csv ]; then
		dinfo=$(sed -e '1d' ${datalist}.csv | sed -e 's/,/|/g')
		#dinfo=$(awk -F, 'BEGIN {OFS="|"} {if (NR!=1) {print $1,$2,$3,$4,$5,$6,$7,$8,$9}}' ${datalist}.csv)
	    fi
	    sed -i "s/\@D${dlbn}/\@D${dlbn}|${dinfo}/g" ${out_name}.gmt
	fi
	unset data
	unset dinfo
    done

    # Convert gmt vector to shapefile
    printf ": boundary: %s.shp\n" $out_name
    ogr2ogr ${out_name}.shp ${out_name}.gmt -overwrite;
done

### End
