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

will output [regions_tile_set.shp](regions_tile_set.geojson) with 24 .25 degree tiles.

### Edit the tile set to only include desired tiles (using a GIS)

[tiles_1_9.shp](tiles_1_9.geojson)

## Generate a Coastline vector (optional)

```bash
fetch_osm_coastline.py -R region_n48x60_w125x10.geojson wa_coast.shp
```

will output [wa_coast.shp](wa_coast.geojson)

## Fetch common datasets

### Bathymetry
#### HydroNOS

```bash
fetches -R tiles_1_9.shp hydronos
```

#### Nautical Charts
```bash
fetches -R tiles_1_9.shp charts
```

#### Multibeam
```bash
fetches -R tiles_1_9.shp multibeam
```

#### Crowd-Sourced Bathymetry
```bash
fetches -R tiles_1_9.shp csb
```

### Topography / Near-shore Bathymetry

#### Digital Coast Lidar
```bash
fetches -R tiles_1_9.shp digital_coast:datatype=lidar
```

#### USGS Lidar

#### USGS DEMs

#### CUDEMs
```bash
fetches -R tiles_1_9.shp CUDEM
```

```
digital_coast/
├── 8483
│   ├── ncei19_n47x00_w0124x25_2018v1.tif
│   ├── ncei19_n47x25_w0124x25_2018v1.tif
│   ├── ncei19_n47x50_w0124x25_2018v1.tif
│   ├── ncei19_n47x50_w0124x50_2018v1.tif
│   ├── ncei19_n47x75_w0124x50_2018v1.tif
│   ├── ncei19_n48x00_w0124x50_2018v1.tif
│   ├── ncei19_n48x00_w0124x75_2018v1.tif
│   ├── ncei19_n48x25_w0124x75_2018v1.tif
│   ├── ncei19_n48x25_w124x25_2021v1.tif
│   ├── ncei19_n48x25_w124x50_2021v1.tif
│   ├── ncei19_n48x50_w0124x75_2018v1.tif
│   └── ncei19_n48x50_w124x50_2021v1.tif
└── 8580
    ├── ncei13_n47x25_w0124x50_2018v1.tif
    ├── ncei13_n47x75_w0124x75_2018v1.tif
    ├── ncei13_n48x25_w0125x00_2018v1.tif
    └── ncei13_n48x50_w0125x00_2018v1.tif
```


## Make a datalist

- see fetches --help for available datasets, or gather and use your own data.

## Generate a test tile

pick a tile and generate an on-the-fly DEM to see what it looks like

either, use the region dimensions of the desired tile or select and export the tile to a new vector using a GIS.

