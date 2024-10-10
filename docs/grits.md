# grits

Filter Digital Elevation Models.

## Synopsis

```
grits [ -hvCMNUWX [ args ] ] DEM ...
```

## Description

Filter DEMs using a variety of metods (Table 1).

## Options
`-M, --module`
Desired grits MODULE and options. (see available Modules below)
Where MODULE is module[:mod_opt=mod_val[:mod_opt1=mod_val1[:...]]]
This can be set multiple times to perform multiple filters.

`-N, --min_z`
Minimum z value (filter data above this value)

`-X, --max_z`
Maximum z value (filter data below this value)

`-U, --uncertainty_mask`
An associated uncertainty raster or band number

`-W, --weight_mask`
An associated weight raster or band number

`-C, --count_mask`
An associated count raster or band number

`--help`
Print the usage text

`--modules`
Display the module descriptions and usage

`--version`
Print the version information

**Table 1.** Filtering odules available in the CUDEM software tool "grits"

|  ***Name***  |  ***Description*** |
|----------------------|----------------------------------|
| blur | apply a gaussian blur to a DEM |
| grdfilter | filter a DEM using the GMT tool `grdfilter` |
| outliers | discover and remove outliers from a DEM |
| flats | remove flattened areas from a DEM |

## Python API

## Examples

### Filter outliers with a weight below 1.2 from a 1/9 arc-second `stacks` DEM (generated with `waffles`)

```bash
$ grits -M outliers --max_weight 1.2 -U 4 -W 3 -C 2 ncei19_n38x25_w123x50_2024v1_stack.tif
```

### Filter 'flattened' areas from a USGS NED DEM

```bash
$ grits -M flats USGS_1_n45w125_20130911.tif
```