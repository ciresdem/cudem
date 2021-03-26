#!/bin/bash

function help () {
echo "usace_interp - A simple script that converts all .xyz files in a directory to .grd files for a provided cell size, fills in no data to a specificed amounts, resamples to finer resolution, and converts to xyz"
	echo "Usage: $0 lastools_path initial_cellsize nodata_fill final_cellsize"
	echo "* lastools_path: <path to lastools bin with lasboundary64.exe>
	e.g., /media/sf_C_win_lx/software/LAStools/bin"
	echo "* initial_cellsize: <cell size in arc-seconds>
	0.00003086420 = 1/9th arc-second 
	0.00009259259 = 1/3rd arc-second
	0.00027777777 = 1 arc-second"
	echo "* nodata_fill: <number of nodata cells to fill>"
	echo "* final_cellsize: <cell size in arc-seconds>
	0.00003086420 = 1/9th arc-second 
	0.00009259259 = 1/3rd arc-second
	0.00027777777 = 1 arc-second"
}

#see if 4 parameters were provided
#show help if not
if [ ${#@} == 4 ]; 
then
	#User inputs    	
	lastools_path=$1
	initial_cellsize=$2
	nodata_fill=$3
	final_cellsize=$4

	echo "Copying LAStools license"
	cp $lastools_path/lastoolslicense.txt lastoolslicense.txt 

	echo "Copying lasboundary64.exe to current directory"
	cp $lastools_path/lasboundary64.exe lasboundary64.exe

	mkdir -p interp

	for i in *.xyz;
	do
		#Create tmp text file of minmax for each xyz file
		#minmax $i > minmax_tmp.txt
	
		#Get minx, maxx, miny, maxy from temporary file
		minx="$(gmt gmtinfo  $i -C | awk '{print $1}')"
		maxx="$(gmt gmtinfo  $i -C | awk '{print $2}')"
		miny="$(gmt gmtinfo  $i -C | awk '{print $3}')"
		maxy="$(gmt gmtinfo  $i -C | awk '{print $4}')"
		#short_name = 

		echo "minx is $minx"
		echo "maxx is $maxx"
		echo "miny is $miny"
		echo "maxy is $maxy"

		echo "Creating Polygon for " $i
		wine ./lasboundary64.exe -i $i -o $(basename $i .xyz).shp -concavity 0.0015

		echo "coverting $i to grd"
		gmt xyz2grd $i -R${minx}/${maxx}/${miny}/${maxy} -G$(basename $i .xyz).grd -I$initial_cellsize/$initial_cellsize
		
		echo "Converting $i to tif"
		gmt grdconvert $(basename $i .xyz).grd $(basename $i .xyz).tif=gd:GTiff
		
		echo "Filling in nodata values up to $nodata_fill cells from data"
		gdal_fillnodata.py -md $nodata_fill $(basename $i .xyz).tif $(basename $i .xyz)_fill.tif
		
		echo "Resampling to target resolution"
		gdalwarp $(basename $i .xyz)_fill.tif -r cubicspline -tr $final_cellsize $final_cellsize‬ $(basename $i .xyz)_fill_resamp.tif

		echo "Clipping to Polygon"
		gdal_rasterize -i -burn nan -l $(basename $i .xyz) $(basename $i .xyz).shp $(basename $i .xyz)_fill_resamp.tif

		echo "Converting to XYZ"
		gdal_translate -of xyz $(basename $i .xyz)_fill_resamp.tif $(basename $i .xyz)_fill_resamp.xyz
		
		echo "Removing Nodata values"
		grep -v nan $(basename $i .xyz)_fill_resamp.xyz | awk '{printf "%.8f %.8f %.2f\n", $1,$2,$3}' > interp/$(basename $i .xyz)_interp.xyz

		rm $(basename $i .xyz).grd
		rm $(basename $i .xyz).tif
		rm $(basename $i .xyz).tif.aux.xml
		rm $(basename $i .xyz)_fill.tif
		rm $(basename $i .xyz)_fill_resamp.tif
		rm $(basename $i .xyz)_fill_resamp.xyz
	
	done

echo "Creating datalist"
cd interp
./create_datalist.sh usace_dredge_interp

else
	help

fi































