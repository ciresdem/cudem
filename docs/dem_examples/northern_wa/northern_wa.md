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

will output [regions_tile_set.shp](region_tile_set.geojson) with 24 .25 degree tiles.

### Edit the tile set to only include desired tiles (using a GIS)

[tiles_1_9.shp](tiles_1_9.geojson)

## Generate a Coastline vector (optional)

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

```waffles -R region_n48x60_w125x10.shp -E .11111111s -M coastline:want_nhd=True:polygonize=True -O wa_coastline -P epsg:4326```

generates wa_coastline.shp and wa_coastline.tif for masking, etc.

## Fetch common data

### Hydronos

```bash
fetches -R tiles_1_9.shp hydronos
```

```
hydronos/
├── bag
│   ├── D00165_MB_VR_MSL_1of1.bag
│   ├── H12219_MB_1m_MLLW_1of2.bag
│   ├── H12219_MB_2m_MLLW_2of2.bag
│   ├── H12220_MB_2m_MLLW_1of3.bag
│   ├── H12220_MB_4m_MLLW_2of3.bag
│   ├── H12220_MB_8m_MLLW_3of3.bag
│   ├── H12221_MB_1m_MLLW_1of3.bag
│   ├── H12221_MB_2m_MLLW_2of3.bag
│   ├── H12221_MB_4m_MLLW_3of3.bag
│   ├── H12222_MB_2m_MLLW_1of3.bag
│   ├── H12222_MB_4m_MLLW_2of3.bag
│   ├── H12222_MB_8m_MLLW_3of3.bag
│   ├── H12223_MB_VR_MLLW_1of1.bag
│   ├── H13412_MB_VR_Ellipsoid.bag
│   ├── H13412_MB_VR_MLLW.bag
│   ├── W00262_MB_VR_MLLW_1of1.bag
│   ├── W00442_MB_128m_MLLW_5of5.bag
│   ├── W00442_MB_128m_MLLW_Combined.bag
│   ├── W00442_MB_16m_MLLW_2of5.bag
│   ├── W00442_MB_32m_MLLW_3of5.bag
│   ├── W00442_MB_64m_MLLW_4of5.bag
│   ├── W00442_MB_8m_MLLW_1of5.bag
│   ├── W00445_MB_128m_MLLW_4of4.bag
│   ├── W00445_MB_128m_MLLW_Combined.bag
│   ├── W00445_MB_16m_MLLW_1of4.bag
│   ├── W00445_MB_32m_MLLW_2of4.bag
│   └── W00445_MB_64m_MLLW_3of4.bag
└── geodas
    ├── H02096.xyz.gz
    ├── H02170.xyz.gz
    ├── H02869.xyz.gz
    ├── H04715.xyz.gz
    ├── H04716.xyz.gz
    ├── H04728.xyz.gz
    ├── H04729.xyz.gz
    ├── H04735.xyz.gz
    ├── H05068.xyz.gz
    ├── H05069.xyz.gz
    ├── H05070.xyz.gz
    ├── H05107.xyz.gz
    ├── H05108.xyz.gz
    ├── H05109.xyz.gz
    ├── H05110.xyz.gz
    ├── H05111.xyz.gz
    ├── H05114.xyz.gz
    ├── H05146.xyz.gz
    ├── H05147.xyz.gz
    ├── H05148.xyz.gz
    ├── H05155.xyz.gz
    ├── H05157.xyz.gz
    ├── H07036.xyz.gz
    ├── H07037.xyz.gz
    ├── H08241.xyz.gz
    ├── H08242.xyz.gz
    ├── H09413.xyz.gz
    ├── H09415.xyz.gz
    ├── H09416.xyz.gz
    ├── H09418.xyz.gz
    ├── H11083.xyz.gz
    ├── H11086.xyz.gz
    ├── H12219.xyz.gz
    ├── H12220.xyz.gz
    ├── H12221.xyz.gz
    ├── H12222.xyz.gz
    └── H12223.xyz.gz
```

### Multibeam
### Nautical Charts
### Digital Coast Lidar
### 

## Make a datalist

- see fetches --help for available datasets, or gather and use your own data.

## Generate a test tile

pick a tile and generate an on-the-fly DEM to see what it looks like

either, use the region dimensions of the desired tile or select and export the tile to a new vector using a GIS.

