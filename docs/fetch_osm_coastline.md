# fetch_osm_coastline.py

Fetch and process an OpenStreetMap coastline to polygon(s)

## Synopsis

```
fetch_osm_coastline.py [-RBli [output-shape-file]]
```

## Description

Fetch a coastline from OpenStreetMap and process it to (a) polygon(s).

## Options

`-R, --region`
Restrict processing to the desired REGION 
Where a REGION is xmin/xmax/ymin/ymax[/zmin/zmax[/wmin/wmax/umin/umax]]
Use '-' to indicate no bounding range; e.g. -R -/-/-/-/-10/10/1/-/-/-
OR an OGR-compatible vector file with regional polygons. 
Where the REGION is /path/to/vector[:zmin/zmax[/wmin/wmax/umin/umax]].
If a vector file is supplied, will use each region found therein.

`-B, --line_buffer`
Buffer the OSM coastline in degrees

`--include_landmask`
include the landmask in the output

`--invert_watermask`
invert the watermask to the landmask

`--help`
Print the usage text

`--version`
Print the version informatio

## Examples

```bash
fetch_osm_coastline.py -R -123.02/-122.73/37.98/38.27 test_msk
```

![](/media/osm_coastpoly.png)

