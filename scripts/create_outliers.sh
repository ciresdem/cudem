for i in *.tif;
do
echo "Identifying Potential Outliers"
echo "DEM is" $i
rm -f thresholds.csv
echo "Calculating min and max thresholds from percentiles"
percentiles_minmax.py $i
min_threshold=`awk -F, '{print $2}' thresholds.csv`
max_threshold=`awk -F, '{print $3}' thresholds.csv`
echo "Creating outliers shapefiles"
outliers_shp.sh $i $min_threshold $max_threshold yes
echo
done
