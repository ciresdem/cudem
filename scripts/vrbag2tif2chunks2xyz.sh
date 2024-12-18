#!/bin/bash
function help () {
echo "bag2tif2chunks2xyz.sh - A script that unzips any .gz and converts all .bag files in a directory to .tif files and then chunks and then converts to xyz."
echo "Usage: $0 chunk_dims resamp cellsize"
echo "* chunk_dims: <number of rows/columns per chunk in resamp resolution>"
echo "* resamp: <resample the tif, yes or no>"
echo "* cellsize: <resampled cell size in arc-seconds>
     	0.0000102880663 = 1/27 arc-second
	0.00003086420 = 1/9th arc-second 
	0.00009259259 = 1/3rd arc-second
	0.00027777777 = 1 arc-second"

}

if [ ${#@} == 3 ]; 
then

    mkdir -p tif
    mkdir -p xyz
    mkdir -p ellip
    
    chunk_dim_x_int=$1
    chunk_dim_y_int=$1
    resamp_input=$2
    resamp_res=$3

    echo "Moving all ellipsoid files to separate folder"
    mv *Ellipsoid* ellip/

    echo unzipping files
    for i in *.gz;
    do
	echo "Unzipping" $i
	gunzip $i
    done

    mkdir -p tif
    mkdir -p xyz
 	
    #get count of bag files
    total_files=$(ls -1 | grep '\.bag$' | wc -l)
    echo "Total number of bag files to process:" $total_files

    file_num=1
    for i in *.bag;
    do
	echo "Processing File" $file_num "out of" $total_files
	echo "Processing" $i

	if [ "$resamp_input" == "yes" ];
	then
	    echo -- Resampling to target resolution
	    echo -- Checking to see if bag is VR and contains supergrids
	    bag_check=`gdalinfo $i -oo MODE=LIST_SUPERGRIDS | grep -e "SUBDATASET_1_NAME" | cut -d "=" -f2`
	    echo "bag_check is" $bag_check
   	    if [[ -z "$bag_check" ]];
	    then
		echo "BAG is NOT VALID VR"
		echo "Reprojecting to NAD83 with provided resamp res"
		gdalwarp $i -dstnodata -999999 -r bilinear -tr $resamp_res $resamp_res -t_srs EPSG:4269 $(basename $i .tif)"_resamp.tif" -overwrite
	    else
		echo "BAG is VALID VR"
		all_vr_bags=`gdalinfo $i -oo MODE=LIST_SUPERGRIDS | grep -e "NAME" | cut -d "=" -f2`
		for this_bag in $all_vr_bags; do
		    gdal_translate $this_bag $i"_resamp_vr".tif
		    echo "Reprojecting to NAD83 with provided resamp res"
		    gdalwarp $i"_resamp_vr".tif -dstnodata -999999 -r bilinear -tr $resamp_res $resamp_res -t_srs EPSG:4269 $(basename $i .tif)"_resamp.tif" -overwrite
		    rm $i"_resamp_vr".tif		    
		done		
		x_res=`gdalinfo $i -oo MODE=LIST_SUPERGRIDS | grep -e "MIN_RESOLUTION_X" | cut -d "=" -f2`
		echo "min x_res is" $x_res
		y_res=`gdalinfo $i -oo MODE=LIST_SUPERGRIDS | grep -e "MIN_RESOLUTION_Y" | cut -d "=" -f2`
		echo "min y_res is" $y_res
		
		if (( $(echo "$x_res < 2" | bc -l) )); 
		then
		    echo "x res is less than 2, setting to 2"
		    x_res=2
		else
		    echo "using default res"
		fi
		
		if (( $(echo "$y_res < 2" | bc -l) )); 
		then
		    echo "y res is less than 2, setting to 2"
		    y_res=2
		else
		    echo "using default res"
		fi
		
		echo "min x_res is" $x_res
		echo "min y_res is" $y_res 
		#exit 1
		gdal_translate $i -oo MODE=RESAMPLED_GRID -oo RESX=$x_res -oo RESY=$y_res  $i"_resamp_vr".tif
		echo "Reprojecting to NAD83 with provided resamp res"
		gdalwarp $i"_resamp_vr".tif -dstnodata -999999 -r bilinear -tr $resamp_res $resamp_res -t_srs EPSG:4269 $(basename $i .tif)"_resamp.tif" -overwrite
		rm $i"_resamp_vr".tif
	    fi
	    
	    input_file=$(basename $i .tif)"_resamp.tif"
	    echo input_file is $input_file
	    
	else
	    echo -- Keeping orig resolution
	    input_file=${i}
	    echo input_file is $input_file
	fi

	x_dim=`gdalinfo $input_file | grep -e "Size is" | awk '{print $3}' | sed 's/.$//'`
	y_dim=`gdalinfo $input_file | grep -e "Size is" | awk '{print $4}'`
	echo x_dim is $x_dim
	echo y_dim is $y_dim

	echo chunk x_dim is $chunk_dim_x_int
	echo chunk y_dim is $chunk_dim_y_int

	#exit 1

	echo
	echo -- Starting Chunk Analysis
	echo

	#initiate chunk names with numbers, starting with 1
	chunk_name="1"
	#remove file extension to get basename from input file
	input_name=${i::-4}
	#starting point for chunks
	xoff=0
	yoff=0

	while [ "$(bc <<< "$xoff < $x_dim")" == "1"  ]; do
	    yoff=0
	    while [ "$(bc <<< "$yoff < $y_dim")" == "1"  ]; do
	    chunk_name_full_tmp=$input_name"_chunk_"$chunk_name"_tmp.tif"
	    chunk_name_full_tmp2=$input_name"_chunk_"$chunk_name"_tmp2.tif"
	    chunk_name_full=$input_name"_chunk_"$chunk_name".tif"
	    echo creating chunk $chunk_name_full
	    echo xoff is $xoff
	    echo yoff is $yoff
	    echo chunk_dim_x_int is $chunk_dim_x_int
	    echo chunk_dim_y_int is $chunk_dim_y_int
	    if [ "$resamp_input" == "yes" ];
	    then
		echo -- Already converted to target coord system, chunking...
		gdal_translate -of GTiff -srcwin $xoff $yoff $chunk_dim_x_int $chunk_dim_y_int $input_file $chunk_name_full -stats
		
	    else
		echo -- Chunking, and then converting to target coord system 
		gdal_translate -b 1 -of GTiff -srcwin $xoff $yoff $chunk_dim_x_int $chunk_dim_y_int $input_file $chunk_name_full_tmp
		gdalwarp $chunk_name_full_tmp -dstnodata -999999 -t_srs EPSG:4269 $chunk_name_full_tmp2
		echo -- Re-Calculating Stats
		gdal_translate -stats $chunk_name_full_tmp2 $chunk_name_full
		
		rm $chunk_name_full_tmp
		rm $chunk_name_full_tmp2
	    fi

	    valid_check=`gdalinfo $chunk_name_full | grep -e "STATISTICS_MAXIMUM"`
	    echo "Valid check is" $valid_check
	    #exit 1

	    if [[ -z "$valid_check" ]];
	    then
		echo "chunk has no data, deleting..."
		rm $chunk_name_full
		#rm $chunk_name_full_tmp
		#rm $chunk_name_full_tmp2
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
	#delete resampled file if it exists
	[ -e $(basename $i .tif)"_resamp.tif" ] && rm $(basename $i .tif)"_resamp.tif"
	echo
	echo
    done
    
else
    help
fi

