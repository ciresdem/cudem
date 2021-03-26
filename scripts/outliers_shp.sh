#!/bin/bash

function help () {
echo "outliers_shp - A script that reclassifies DEM based on user-defined thresholds to identify outlier cells with an option to create shapefile of output"
	echo "Usage: $0 DEM max_threshold min_threshold shp"
	echo "* DEM: <DEM to reclassify>"
	echo "* min_threshold: <find values less than vertical threshold in DEM units>"
	echo "* max_threshold: <find values greater than vertical threshold in DEM units>"
	echo "* shp: <yes or no; create shapefile of outliers>"
	echo
}

DEM=$1
min_threshold=$2
max_threshold=$3
shp=$4
#replace decimal with _ and negative sign with neg
min_threshold_str_tmp=${min_threshold/./_}
max_threshold_str_tmp=${max_threshold/./_}
replace="neg"
min_threshold_str=${min_threshold_str_tmp/-/$replace}
max_threshold_str=${max_threshold_str_tmp/-/$replace}
#echo $min_threshold_str 
threshold_str=$min_threshold_str"_"$max_threshold_str
#echo $threshold_str

#see if 4 parameters were provided
#show help if not
mkdir -p outliers

if [ ${#@} == 4 ]; 
then
	#reclassify
	echo "Reclassifying DEM to 1s and 0s based on max threshold of" $max_threshold
	gdal_calc.py -A $DEM --outfile=$(basename $DEM .tif)"_rc_"$max_threshold_str".tif" --calc="1*(A>$max_threshold)+0*(A<=$max_threshold)"
	echo "Reclassifying DEM to 1s and 0s based on min threshold of" $min_threshold
	gdal_calc.py -A $DEM --outfile=$(basename $DEM .tif)"_rc_"$min_threshold_str".tif" --calc="1*(A<$min_threshold)+0*(A>=$min_threshold)"
	echo "Adding Reclassifed DEMs"
	gdal_calc.py -A $(basename $DEM .tif)"_rc_"$max_threshold_str".tif" $DEM -B $(basename $DEM .tif)"_rc_"$min_threshold_str".tif" --outfile=$(basename $DEM .tif)"_rc_"$threshold_str".tif" --calc="A + B"
	
	#shp creation
	if [ "$shp" == "yes" ];
		then
		echo "Creating shp from sum raster"
		gdal_polygonize.py $(basename $DEM .tif)"_rc_"$threshold_str".tif" -8 -f "ESRI Shapefile" $(basename $DEM .tif)"_"$threshold_str"_outliers_all.shp"
		rm $(basename $DEM .tif)"_rc_"$threshold_str".tif"
		rm $(basename $DEM .tif)"_rc_"$max_threshold_str".tif"
		rm $(basename $DEM .tif)"_rc_"$min_threshold_str".tif"
		
		echo "Creating polygons of only outlier cells"
		ogr2ogr -dialect SQLITE -sql "SELECT * FROM $(basename $DEM .tif)"_"$threshold_str"_outliers_all" WHERE DN='1'" outliers/$(basename $DEM .tif)"_"$threshold_str"_outliers.shp" $(basename $DEM .tif)"_"$threshold_str"_outliers_all.shp"
		rm $(basename $DEM .tif)"_"$threshold_str"_outliers_all.shp"
		rm $(basename $DEM .tif)"_"$threshold_str"_outliers_all.shx"
		rm $(basename $DEM .tif)"_"$threshold_str"_outliers_all.dbf"
		rm $(basename $DEM .tif)"_"$threshold_str"_outliers_all.prj"
	else
		echo "Not creating shp, so keeping reclassified rasters"
		mv $(basename $DEM .tif)"_rc_"$threshold_str".tif" outliers/$(basename $DEM .tif)"_rc_"$threshold_str".tif"
		mv $(basename $DEM .tif)"_rc_"$max_threshold_str".tif" outliers/$(basename $DEM .tif)"_rc_"$max_threshold_str".tif"
		mv $(basename $DEM .tif)"_rc_"$min_threshold_str".tif" outliers/$(basename $DEM .tif)"_rc_"$min_threshold_str".tif"
	fi
else
	help
fi
