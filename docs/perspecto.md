# perspecto

Generate images of DEMs

## Synopsis

```
perspecto [ -hvCM [ args ] ] DEM ...
```

## Description

Generate images of DEMs, including perspectives, hillshades, etc. (Table 1)

## Options
`-C, --cpt`

> Color Pallette file (if not specified will auto-generate ETOPO CPT)

`-M, --module`

> Desired perspecto MODULE and options. (see available Modules below)\
> Where MODULE is module[:mod_opt=mod_val[:mod_opt1=mod_val1[:...]]]

`--min_z`

> Minimum z value to use in CPT

`--max_z`

> Maximum z value to use in CPT

`--help`

> Print the usage text

`--modules`

> Display the module descriptions and usage

`--version`

> Print the version information

## Modules

**Table 1.** Modules available in the CUDEM software tool "perspecto"

|  ***Name***  |  ***Description*** |
|----------------------|----------------------------------|
| hillshade | generate a DEM hillshade (req. gdal/imagemagick) |
| perspective | generate a DEM perspective (req. POVRay) |
| sphere | generate a DEM on a sphere |
| figure1 | generate a DEM figure (req. GMT) |
| colorbar | generate a colorbar image based on the input DEM/CPT |

### hillshade
### perspective
### sphere
### figure1
### colorbar

## Python API

```python
from cudem import perspecto

dem_path = '/my_dems/dem.tif'
p = perspecto.PerspectoFactory(mod='hillshade', src_dem=dem_path, min_z=-1000, max_z=100)._acquire_module()
p.run()
```

## Examples

```bash
perspecto my_dem.tif -M hillshade -C GMT_wysiwyg
```

## Gallery

### ETOPO 2022 *sphere*
![](/media/etopo22_northAmerica.png)