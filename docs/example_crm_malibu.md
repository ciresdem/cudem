# Coastal Releif Model (CRM) at 1 arc-second DEM off the coast of Malibu

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

![](/media/mugu_crm_test_hs.png)
**Figure 1.** DEM hillshade generated from CUDEM code example.