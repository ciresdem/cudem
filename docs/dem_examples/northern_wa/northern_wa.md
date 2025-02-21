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

## Fetch common datasets and create datalists

Use the [fetches](/docs/fetches.md) command to download common datasets. Use the `-H` switch to fetch data in multiple threads. The fetched data will be located in the current working directory in a directory named after the fetch module. Using the tiled vector file generated above, `tiles_1_9.shp`, we will fetch data for each of the tiles. The data won't be duplicated, so a file fetched for one tile will be skipped if it's present in other tile.

Fetched data can then be either processed to XYZ or used as-is in a datalist. Data originating in raster format will be used as-is, while some other datasets will be processed to XYZ and common datums for use in the datalists.

### Bathymetry
#### HydroNOS

```bash
fetches -R tiles_1_9.shp hydronos -H3
```

#### Nautical Charts
```bash
fetches -R tiles_1_9.shp charts -H3
```

Process the fetched data to XYZ and NAD83/NAVD88
```bash
dlim -R tiles_1_9.shp --archive -P epsg:4269+5703 charts
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

There is a lot of data on the Digital Coast Access Viewer, we don't necessarily need it all. Running the following command will fetch all the available lidar in the AOI and unless you have a lot of disk space, it may fill up and fail at some point. You can search and discover the survey ID for desired lidar surveys using the [DAV](https://coast.noaa.gov/dataviewer/#/lidar/search/-13932941.515821926,5948023.793039274,-13788628.406419514,6198737.245814653).

```bash
fetches -R tiles_1_9.shp digital_coast:datatype=lidar -H3
```

Otherwise, we can determine which surveys we want specifically and fetch them using the ':where="ID=<survey-ID>"' option in the digital_coast fetches module. We can either run a command for each survey ID we want, or string multiple together with 'OR': ':where="ID=<ID1> OR ID=<ID2> ..."'

```bash
fetches -R tiles_1_9.shp digital_coast:where="ID=9703 OR ID=" -H3
```

We can also gather lidar surveys by date, if we want to just use lidar since 2018, for example, we can do:

```bash
fetches -R tiles_1_9.shp digital_coast:datatype=lidar:where="YEAR>2017" -H3
```

See `fetches --modules digital_coast` to see all the possible query fields.

For this example, we'll fetch specific lidar surveys, with the following command:

```bash
fetches -R tiles_1_9.shp digital_coast:where="ID=9703 OR ID=10116 OR ID=9072 OR ID=9512 OR ID=4989 OR ID=6263 OR ID=5008 OR ID=2492 OR ID=2508 OR ID=2603 OR ID=8607 OR ID=2584 OR ID=2482" -H3
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

