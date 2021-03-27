#!/bin/sh
function help () {
echo "tif2chunks2xyz.sh - A script that converts all tif files in a directory to smaller chunks and then to xyz."
	echo "Usage: $0 chunk_dims resamp cellsize"
	echo "* chunk_dims: <number of rows/columns per chunk in resamp resolution>"
	echo "* resamp: <resample the tif, yes or no>"
	echo "* cellsize: <resampled cell size in arc-seconds>
	0.0000102880663 = 1/27 arc-second
	0.000030864199 = 1/9th arc-second 
	0.000092592596 = 1/3rd arc-second
	1 = 1 arc-second"

}

if [ ${#@} == 3 ]; 
then


mkdir -p tif
mkdir -p xyz

#Example, create chunks 150 rows by 150 cols
chunk_dim_x_int=$1
chunk_dim_y_int=$1
resamp_input=$2
resamp_res=$3

#get count of tif files
total_files=$(ls -1 | grep '\.tif$' | wc -l)
echo "Total number of tif files to process:" $total_files

file_num=1
for i in *.tif;
do
	echo "Processing File" $file_num "out of" $total_files
	echo "Processing" $i

	if [ "$resamp_input" == "yes" ];
	then
		echo -- Resampling to target resolution
		#exit 1
		gdalwarp $i -dstnodata -999999 -r cubicspline -tr $resamp_res $resamp_res -t_srs EPSG:4269 $(basename $i .tif)"_resamp.tif" -overwrite
	else
		echo -- Keeping orig resolution
		#exit 1
		cp $i $(basename $i .tif)"_resamp.tif"
	fi

	#get input grid dimensions
	x_dim=`gdalinfo $(basename $i .tif)"_resamp.tif" | grep -e "Size is" | awk '{print $3}' | sed 's/.$//'`
	y_dim=`gdalinfo $(basename $i .tif)"_resamp.tif" | grep -e "Size is" | awk '{print $4}'`
	
	echo chunk x_dim is $chunk_dim_x_int
	echo chunk y_dim is $chunk_dim_y_int

	echo
	echo -- Starting Chunk Analysis
	echo

	#initiate chunk names with numbers, starting with 1
	chunk_name="1"
	#remove file extension to get basename from input file
	input_name=${i::-4}
	#starting point for tiling
	xoff=0
	yoff=0

	while [ "$(bc <<< "$xoff < $x_dim")" == "1"  ]; do
	    yoff=0
	    while [ "$(bc <<< "$yoff < $y_dim")" == "1"  ]; do
	    chunk_name_full=$input_name"_chunk_"$chunk_name".tif"
	    echo creating chunk $chunk_name_full
	    echo xoff is $xoff
	    echo yoff is $yoff
	    echo chunk_dim_x_int is $chunk_dim_x_int
	    echo chunk_dim_y_int $chunk_dim_y_int
	    gdal_translate -of GTiff -srcwin $xoff $yoff $chunk_dim_x_int $chunk_dim_y_int $(basename $i .tif)"_resamp.tif" $chunk_name_full -stats
	    
	    valid_check=
	    valid_check=`gdalinfo $chunk_name_full | grep -e "STATISTICS_MAXIMUM"`
	    echo "Valid check is" $valid_check
	    #exit 1

	    if [[ -z "$valid_check" ]];
		then
		      echo "chunk has no data, deleting..."
		      rm $chunk_name_full
		else
		      echo "chunk has data, keeping..."
		      echo -- Converting to xyz
		      gdal_translate -of XYZ $chunk_name_full $chunk_name_full"_tmp.xyz"
		      awk '{if ($3 > -99999 && $3 < 99999 ) {printf "%.8f %.8f %.2f\n", $1,$2,$3}}' $chunk_name_full"_tmp.xyz" > $chunk_name_full".xyz"
		      echo -- Converted to xyz
		      rm $chunk_name_full"_tmp.xyz" 
		      mv $chunk_name_full tif/$chunk_name_full
		      mv $chunk_name_full".xyz" xyz/$chunk_name_full".xyz"
		fi
		yoff=$(echo "$yoff+$chunk_dim_y_int" | bc)
	    chunk_name=$((chunk_name+1))
	    done
	  xoff=$(echo "$xoff+$chunk_dim_x_int" | bc)
	done
	file_num=$((file_num + 1))
	rm $(basename $i .tif)"_resamp.tif"
	echo
	echo
done

else
	help
fi

