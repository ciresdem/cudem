# perspecto

Generate images of DEMs

## Synopsis

```
perspecto [ -hvCM [ args ] ] DEM ...

perspecto ... <mod>:key=val:key=val...
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

|  ***Name***  |  ***Description*** | ***Module Options*** |
|----------------------|----------------------------------|----------------------------------|
| hillshade | generate a DEM hillshade (req. gdal/imagemagick) | :vertical_exaggeration=1 :projection=4326 :azimuth=315 :altitude=45 |
| perspective | generate a DEM perspective (req. POVRay) | :cam_azimuth=-130 :cam_elevation=30 :cam_distance=265 :cam_view_angle=40 :light_elevation=20 :light_distance=10000 :vertical_exaggeration=2 |
| sphere | generate a DEM on a sphere | :cam_azimuth=310 :cam_elevation=27 :cam_distance=8 :cam_view_angle=33 :center_lat=None :center_long=None |
| figure1 | generate a DEM figure (req. GMT) | :perspective=False :vertical_exaggeration=1.5 :interval=100 :azimuth=-130 :elevation=30 | 
| colorbar | generate a colorbar image based on the input DEM/CPT | :colorbar_text='Elevation' :width=10 :height=2 |

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

### Hillshade image of CRM volume 7

- Download the CRM Volume 7 DEM

```bash
wget https://www.ngdc.noaa.gov/thredds/fileServer/crm/cudem/crm_vol7_2024.nc
```

- Generate the hillshade image

```bash
perspecto -C crmVol7.cpt -M hillshade --min_z -4500 --max_z 2260 crm_vol7_2024_3as.tif
```

![](/media/crm_vol7_2024_hs.png)

### Perspective image of CRM volume 7

- Download the CRM Volume 7 DEM

```bash
wget https://www.ngdc.noaa.gov/thredds/fileServer/crm/cudem/crm_vol7_2024.nc
```

- Generate the hillshade image

```bash
perspecto -C crmVol7.cpt -M perspective --min_z -4500 --max_z 2260 crm_vol7_2024_3as.tif
```

![](/media/crm_vol7_2024_perspective.png)

### ETOPO 2022 *sphere*

- Download the ETOPO 2022 DEM

```bash
wget https://www.ngdc.noaa.gov/thredds/fileServer/global/ETOPO2022/30s/30s_surface_elev_netcdf/ETOPO_2022_v1_30s_N90W180_surface.nc
```

- Warp the DEM to 60 arc-seconds for memory management in sphere generation

```bash
gdalwarp ETOPO_2022_v1_30s_N90W180_surface.nc ETOPO_2022_v1_60s_N90W180_surface.tif -tr 0.016666666666666666 0.016666666666666666
```

- Generate the Sphere image

```bash
perspecto -M sphere ETOPO_2022_v1_60s_N90W180_surface.tif
```

![](/media/etopo22_northAmerica.png)