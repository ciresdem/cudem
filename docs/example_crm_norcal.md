## Tiled CRMs in Northern California

Generate a .25 degree tileset for an area in northern California

```
regions -R-123.5/-122.75/37.75/38.75 -T .25
```

Output:

regions_tile_set.shp

Generate DEMs of all the tiles

```
waffles -R regions_tile_set.shp -E 1s -M cudem:pre_count=1:filter_outliers=85:min_weight=.6 -O nocal -P epsg:4326+5703 -p -w -X0:5 mar_grav,-106:bathy_only=True,.001 charts,-200,.01 hydronos:datatype=xyz,-202,.1 copernicus,-103,.6 ehydro,-203,.6
```

Output:

nocal1_n38x00_w123x00_2024v1.tif

nocal1_n38x00_w123x50_2024v1.tif

nocal1_n38x25_w123x25_2024v1.tif

nocal1_n38x50_w123x00_2024v1.tif

nocal1_n38x50_w123x50_2024v1.tif

nocal1_n38x75_w123x25_2024v1.tif

nocal1_n38x00_w123x25_2024v1.tif

nocal1_n38x25_w123x00_2024v1.tif

nocal1_n38x25_w123x50_2024v1.tif

nocal1_n38x50_w123x25_2024v1.tif

nocal1_n38x75_w123x00_2024v1.tif

nocal1_n38x75_w123x50_2024v1.tif

Generate hillshades of all the DEMs, setting min/max z of the CPT so the hillshade edges match:

```
for i in nocal1*v1.tif; do perspecto -M hillshade $i --min_z -360 --max_z 800; done
```

![](/media/nocal_hs_test.png)
**Figure 1.** DEM hillshade generated from CUDEM code example.