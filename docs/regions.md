# regions

Process and generate regions

## Synopsis

```
regions [ -hqJPRT [ args ] ]...
```

## Description

Process regions, which are bounding-boxes made up of corner coordiates in the format `x-min/x-max/y-min/y-max`.

## Options

`-R, --region`
The desired REGION 
Where a REGION is xmin/xmax/ymin/ymax[/zmin/zmax[/wmin/wmax/umin/umax]]
Use '-' to indicate no bounding range; e.g. -R -/-/-/-/-10/10/1/-/-/-
OR an OGR-compatible vector file with regional polygons. 
Where the REGION is /path/to/vector[:zmin/zmax[/wmin/wmax/umin/umax]].
If a vector file is supplied, will use each region found therein.

`-J, --s_srs`
Set the SOURCE projection.

`-P, --t_srs`
Set the TARGET projection.

`-T, --tile_set`
Generate a TILESET from the input region. (set incrememnt here)

`-B, --buffer`
BUFFER the region with a buffer-value.

`-e, --echo`
ECHO the <processed> region

`-n, --name`
Print the region as a NAME

`-m, --merge`
MERGE all regions into a single region

`--quiet`
Lower the verbosity to a quiet

`--help`
Print the usage text

`--version`
Print the version information
  
## Python API

## Examples