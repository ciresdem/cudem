# Overview

DEM generation example of Northern Washington.

## Specifications

| Region | Tile-Size | Cell-size | Horz Projection | Vert Projection |
|---|---|---|---|---|
| -R-125/-124/47/48.5 | .25 degrees | 1/9 Arc-Second (~3m) | NAD83 | NAVD88 |


## Generate the region vectors

### Generate a full-region buffered vector

This is for data fetching, etc.

Add a slight buffer (.1 degree) to ensure coverage of fetched data.

```bash
regions -R -125/-124/47/48.5 -B .1
```

will output [region_n48x60_w125x10.shp](region_n48x60_w125x10.geojson)

### Generate the .25 degree tiles

these will be the extents of each DEM generated.

```regions -R -121/-114/31.5/35 -T .25```

will output a shapefile named "regions_tile_set.shp" with 392 .25 degree tiles.

## Generate a Coastline vector

Check the waffles coastline module options to determine what we need:

```
$ waffles --modules coastline
waffles modules:
% waffles ... <mod>:key=val:key=val...

  coastline     COASTLINE (land/etc-mask) generation
    
    Generate a coastline (land/etc-mask) using a variety of sources. 
    User data can be provided to provide source for further land masking. 
    Output raster will mask land-areas as 1 and oceans/(lakes/buildings) as 0.
    Output vector will polygonize land-areas.
    
    -----------
    Parameters:
    
    want_gmrt=[True/False] - use GMRT to fill background (will use Copernicus otherwise)
    want_copernicus=[True/False] - use COPERNICUS to fill background
    want_nhd=[True/False] - use high-resolution NHD to fill US coastal zones
    want_lakes=[True/False] - mask LAKES using HYDROLAKES
    invert_lakes=[True/False] - invert the lake mask (invert=True to remove lakes from the waterbodies)
    want_buildings=[True/False] - mask BUILDINGS using OSM
    osm_tries=[val] - OSM max server attempts
    min_building_length=[val] - only use buildings larger than val
    want_wsf=[True/False] - mask BUILDINGS using WSF
    invert=[True/False] - invert the output results
    polygonize=[True/False] - polygonize the output
    min_weight=[val] - weight applied to fetched coastal data

    < coastline:want_gmrt=False:want_nhd=True:want_lakes=False:want_buildings=False:invert=False:polygonize=True >
```

```waffles -R region_n35x10_w121x10.shp -E 1s -M coastline -O socal_coastline -P epsg:4326```

generates socal_coastline.shp and socal_coastline.tif for masking, etc.

## Make a datalist

- see fetches --help for available datasets, or gather and use your own data.


## Generate a test tile

pick a tile and generate an on-the-fly DEM to see what it looks like

either, use the region dimensions of the desired tile or select and export the tile to a new vector using a GIS.

