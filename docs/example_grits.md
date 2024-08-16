# Automatic gridding and filtering of Multibeam data

First, we use the 'multibeam' fetches module as a datalist entry in a waffles command to fetch and grid the
multibeam data in our region, then using that result as a datalist entry in a waffles command to grid the data
using the 'IDW' waffles gridding module.

```
waffles -R-123.25/-123/48.25/48.5 -E .11111111s -O mb -p -P epsg:4326 -m -u -M stacks multibeam
```

```
waffles -R-123.25/-123/48.25/48.5 -E .11111111s -O mb_idw -p -P epsg:4326 -m -u -M IDW:radius=10 mb19_n48x50_w123x25_2024v1.tif
```

![](/media/mb_idw19_n48x50_w123x25_2024v1_nofilter.png)

**Figure 1.** Auto-gridded raw multibeam data

The above results show some major artifacts from the raw multibeam data, so we then filter that output using
the grits 'outlier' filter and re-generate the IDW grid with waffles:

```
grits mb19_n48x50_w123x25_2024v1_filtered.tif -M outliers
```

```
waffles $(dlim -r mb19_n48x50_w123x25_2024v1.tif) -M IDW:radius=5 -O mb_idw -E .11111111s -p -P epsg:4326 mb_data19_n_w_filtered.tif
```

![](/media/mb_idw19_n48x50_w123x25_2024v1_out75.png)

**Figure 2.** Auto-grided filtered multibeam data

We can combine all the above steps using a single waffles command to fetch the multibeam, filter it using the
grits 'outlier' module with the waffles -T switch and generate an IDW DEM with the waffles 'IDW' module.

```
waffles -R-123.25/-123/48.25/48.5 -E .11111111s -O mb -p -P epsg:4326 -m -u -M IDW:radius=10 -T outliers:stacks=True multibeam
```
