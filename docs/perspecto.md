# Perspecto (Visualization)

**Perspecto** is the visualization module within the CUDEM software suite, designed to generate high-quality images and 2D/3D representations from Digital Elevation Models (DEMs). It serves as a unified interface for various rendering backends, including GDAL, POV-Ray, and GMT/PyGMT.

Designed to work seamlessly with the `waffles` module, Perspecto can accept either an existing raster file or a Waffles JSON configuration file as input. If provided a JSON config, it can automatically generate the requested DEM before visualizing it.

## Summary

### Visualization Modules

`perspecto` uses a modular factory system to support different rendering types:

* **Hillshading:**
  * `hillshade`: Generates standard hillshade images to visualize terrain relief.
  * `hillshade2`: An alternative command-line driven hillshade implementation.


* **3D Rendering:**
  * `perspective`: Creates 3D perspective views of the terrain, utilizing POV-Ray ray-tracing capabilities.
   * `sphere`: Generates spherical visualizations, useful for global or planetary datasets.


* **PyGMT Integration:**
  * `figure1`: Generates GMT-based figures (requires PyGMT).
  * `colorbar`: Generates color bars for maps (requires PyGMT).



### Advanced Color Palette (CPT) Management

Perspecto includes robust tools for managing color tables, essential for creating color relief maps:

* **Auto-Generation:** Can automatically generate an "ETOPO" style color palette based on the min/max elevation of the input DEM.
* **Fetching:** Can fetch named CPTs from **cpt-city**.
* **Processing:** Supports re-scaling CPTs to specific Z-ranges and "splitting" CPTs (e.g., creating distinct color ramps above and below zero for land/sea distinctions).


## Usage

The module is accessed via the command line, requiring a module selection and an input DEM (or Waffles config).

```bash
perspecto input_dem.tif -M hillshade

```

**Common Options:**

* `-M, --module`: Select the visualization type (e.g., `hillshade`, `perspective`).
* `-C, --cpt`: Specify a custom Color Palette Table file.
* `-Z, --split-cpt`: Split the color palette at a specific value (e.g., 0 for coastlines).
* `--min_z / --max_z`: Force specific elevation ranges for color scaling.

If `input_dem` is a JSON file (Waffles config), Perspecto will attempt to parse it and generate the DEM using the `waffles` module if the file does not already exist.

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