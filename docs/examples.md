# Examples

## Coastal Releif Model (CRM) at 1 arc-second DEM off the coast of Malibu

Generate a CRM of the Point Mugu area of Malibu using remote datalists provided from fetches.
Datasets included in this example are Marine Gravity (mar_grav), Nautical Charts (charts),
NOS Soundings (hydronos:datatype=xyz), Copernicus (copernicus), USACE surveys (ehydro) and
NOS BAG surveys (hydronos:datatype=bag)

```
waffles -R -119.25/-119/34/34.25 -E 1s -M cudem:pre_count=1:filter_outliers=85:min_weight=.6 -O mugu_crm_test -P epsg:4326+5703 -w -X0:5 mar_grav,-106:bathy_only=True,.001 charts,-200,.01 hydronos:datatype=xyz,-202,.1 copernicus,-103,.6 ehydro,-203,.6 hydronos:datatype=bag,-202,.75
```

Output:

mugu_crm_test.tif

Generate a hillshade of the generated DEM:

```
perspecto -M hillshade mugu_crm_test.tif
```

Output:

mugu_crm_test_hs.tif

![](media/mugu_crm_test_hs.png)
**Figure 1.** DEM hillshade generated from CUDEM code example.

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

![](media/nocal_hs_test.png)
**Figure 2.** DEM hillshade generated from CUDEM code example.

## DEM and Uncertainty for the Bahamas

We provide an example of code to download and process depth and
elevation data from the "gebco", "copernicus", and "mar_grav" data
sources in Table 1, and then generate a DEM and accompanying
interpolation uncertainty grid. First, we download all the data, and
then process the data by masking out GEBCO data where TID is equal to 0
(Land) or 40 (Predicted based on satellite-derived gravity data - depth
value is an interpolated value guided by satellite-derived gravity
data). See the GEBCO Cookbook chapter on GEBCO TID generation for more
information.

We then stack the raster data sources with higher weighted datasets
masking out lower weighted datasets and apply spline interpolation to
the bathymetry by limiting the interpolation to an upper limit value of
zero and clipping the resulting grid to an automatically generated
coastline from the Copernicus data.

Lastly, we combine this bathymetric surface with the other data sources
that have a weighting above the min_weight specification (in this
example, min_weight=.15), and apply spline interpolation with a weighted
averages of these data sources to generate the final integrated
bathymetric-topographic DEM.

**To generate the configuration file that is then used to generate the
15 arc-second resolution DEM and the interpolation uncertainty grid,
execute the following command:**

```
waffles -R -79.5/-76.5/25/28 -E 15s -O gebcop -p -w -m -P epsg:4326 -M cudem:min_weight=.15:landmask=True:pre_count=1 gebco,-102:exclude_tid=0/40,.015 copernicus,-103,10 mar_grav,-106:bathy_only=True,.001 gmrt,-101:swath_only=True,1.5 --config
```

The contents of the generated config file configuration file
(gebcop15_n28x00_w079x50_2023v1.json) are as follows:

```
{
    "mod": "cudem:min_weight=.15:landmask=True:pre_count=1",
    "mod_name": "cudem",
    "mod_args": {
        "min_weight": ".15",
        "landmask": true,
        "pre_count": "1"
    },
    "kwargs": {
        "verbose": true,
        "sample": "bilinear",
        "xsample": null,
        "ysample": null,
        "dst_srs": "epsg:4326",
        "srs_transform": false,
        "fltr": [],
        "name": "gebcop15_n28x00_w079x50_2023v1",
        "cache_dir": "/home/ncei/Projects/crm/software/.cudem_cache",
        "ndv": -9999,
        "xinc": 0.004166666666666667,
        "yinc": 0.004166666666666667,
        "want_weight": true,
        "want_mask": true,
        "data": [
            "gebco,-102:exclude_tid=0/40,.015,0",
            "copernicus,-103,10,0",
            "mar_grav,-106:bathy_only=True,.001,0",
            "gmrt,-101:swath_only=True,1.5,0"
        ],
        "src_region": [
            -79.5,
            -76.5,
            25.0,
            28.0
        ]
    }
}
```

**To generate the DEM, execute the following command:**

```
waffles -G gebcop15_n28x00_w079x50_2023v1.json
```

Output:

gebcop3_n28x00_w079x50_2023v1.json - json config file

gebcop3_n28x00_w079x50_2023v1_msk.tif - data mask

gebcop3_n28x00_w079x50_2023v1_stack.tif - stacks data \'stacked\'
output

gebcop3_n28x00_w079x50_2023v1_u.tif - uncertainty grid

gebcop3_n28x00_w079x50_2023v1.tif - final DEM

![](media/gebcop3_n28x00_w079x50_2023v1_figure1.png)

**Figure 3.** Final DEM generated from CUDEM code example.

![](media/gebcop3_n28x00_w079x50_2023v1_u_figure1.png)

**Figure 4.** Interpolation uncertainty grid generated from the best-fit
interpolation uncertainty equation applied to the distance to the
nearest measurement raster.

## Automatic gridding and filtering of Multibeam data

First, we use the 'multibeam' fetches module as a datalist entry in a waffles command to fetch and grid the
multibeam data in our region, the using that result as a datalist entry in a waffles command to grid the data
using the 'IDW' waffles gridding module.

```
waffles -R-123.25/-123/48.25/48.5 -E .11111111s -O mb -p -P epsg:4326 -m -u -M stacks multibeam
```

```
waffles -R-123.25/-123/48.25/48.5 -E .11111111s -O mb_idw -p -P epsg:4326 -m -u -M IDW:radius=10 mb19_n48x50_w123x25_2024v1.tif
```

![](media/mb_idw19_n48x50_w123x25_2024v1_nofilter.png)

**Figure 5.** Auto-gridded raw multibeam data

The above results show some major artifacts from the raw multibeam data, so we then filter that output using
the grits 'outlier' filter and re-generate the IDW grid with waffles:

```
grits mb19_n48x50_w123x25_2024v1_filtered.tif -M outliers
```

```
waffles $(dlim -r mb19_n48x50_w123x25_2024v1.tif) -M IDW:radius=5 -O mb_idw -E .11111111s -p -P epsg:4326 mb_data19_n_w_filtered.tif
```

![](media/mb_idw19_n48x50_w123x25_2024v1_out75.png)

**Figure 6.** Auto-grided filtered multibeam data

We can combine all the above steps using a single waffles command to fetch the multibeam, filter it using the
grits 'outlier' module with the waffles -T switch and generate an IDW DEM with the waffles 'IDW' module.

```
waffles -R-123.25/-123/48.25/48.5 -E .11111111s -O mb -p -P epsg:4326 -m -u -M IDW:radius=10 -T outliers:stacks=True multibeam
```
