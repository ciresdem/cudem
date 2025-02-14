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

Outputs: [region_n48x60_w125x10.shp](region_n48x60_w125x10.geojson)

![](wa_region.png)

### Generate the .25 degree tiles

these will be the extents of each DEM generated.

```regions -R -121/-114/31.5/35 -T .25```

Outputs: [regions_tile_set.shp](regions_tile_set.geojson) with 24 .25 degree tiles.

![](wa_region_tiles.png)

### Edit the tile set to only include desired tiles (using a GIS)

[tiles_1_9.shp](tiles_1_9.geojson)

![](wa_tiles_1_9.png)

## Generate a Coastline vector (optional)

```bash
fetch_osm_coastline.py -R region_n48x60_w125x10.geojson wa_coast.shp
```
Outputs: [wa_coast.shp](wa_coast.geojson)

![](wa_coast.png)

## Fetch common datasets

Use the [fetches](/docs/fetches.md) command to download common datasets. Use the `-H` switch to fetch data in multiple threads. The fetched data will be located in the current working directory in a directory named after the fetch module. Using the tiled vector file generated above, `tiles_1_9.shp`, we will fetch data for each of the tiles. The data won't be duplicated, so a file fetched for one tile will be skipped if it's present in other tile.

### Bathymetry
#### HydroNOS

```bash
fetches -R tiles_1_9.shp hydronos -H3
```

#### Nautical Charts
```bash
fetches -R tiles_1_9.shp charts -H3
```

#### Multibeam
```bash
fetches -R tiles_1_9.shp multibeam -H3
```

#### EHydro
```bash
fetches -R tiles_1_9.shp ehydro -H3
```

#### Crowd-Sourced Bathymetry
```bash
fetches -R tiles_1_9.shp csb -H3
```

### Topography / Near-shore Bathymetry

#### Digital Coast Lidar

There is a lot of data on the Digital Coast Access Viewer, we don't necessarily need it all. Running the following command will fetch all the available lidar in the AOI and unless you have a lot of disk space, it may fill up and fail at some point. You can search and discover the survey ID for desired lidar surveys using the [DAV](https://coast.noaa.gov/dataviewer/#/lidar/search/-13929884.03469052,5827559.036461838,-13604568.042308811,6204240.711851186)

```bash
fetches -R tiles_1_9.shp digital_coast:datatype=lidar -H3
```

Otherwise, we can determine which surveys we want specifically and fetch them using the ':where="ID=<survey-ID>"' option in the digital_coast fetches module. We can either run a command for each survey ID we want, or string multiple together with 'OR': ':where="ID=<ID1> OR ID=<ID2> ..."'

```bash
fetches -R tiles_1_9.shp digital_coast:where="ID=9703 OR ID=" -H3
```

#### USGS Lidar

#### USGS DEMs
```bash
fetches -R tiles_1_9.shp ned1
```

```bash
fetches -R tiles_1_9.shp CoNED
```

#### CUDEMs
```bash
fetches -R tiles_1_9.shp CUDEM -H3
```

## Make a datalist

- see fetches --help for available datasets, or gather and use your own data.

## Generate a test tile

pick a tile and generate an on-the-fly DEM to see what it looks like

either, use the region dimensions of the desired tile or select and export the tile to a new vector using a GIS.

