# DEM and Uncertainty for the Bahamas

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

![](/media/gebcop3_n28x00_w079x50_2023v1_figure1.png)

**Figure 1.** Final DEM generated from CUDEM code example.

![](/media/gebcop3_n28x00_w079x50_2023v1_u_figure1.png)

**Figure 2.** Interpolation uncertainty grid generated from the best-fit
interpolation uncertainty equation applied to the distance to the
nearest measurement raster.