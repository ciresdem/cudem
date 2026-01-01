# Dlim (Data Lists IMproved)

**Dlim** is the data management and processing module within the CUDEM software suite. It is designed to handle the complex logistics of organizing, cleaning, and standardizing diverse geospatial data sources before they enter the gridding pipeline. Dlim acts as a universal adapter, converting various input formats—from simple XYZ text files to complex lidar point clouds or remote datasets—into a unified stream of data ready for interpolation.

The module operates on the concept of **Data Lists** (`.datalist`), which are text-based configuration files that define a collection of datasets, their locations (local paths or remote URLs), and specific processing instructions for each entry.

## Summary

### Unified Data Abstraction

* **Format Agnostic:** Dlim can ingest and normalize a wide variety of geospatial formats, including:
  * **Raster:** GeoTIFF, NetCDF, BAG, HDF5, VRT.
  * **Point Cloud:** LAS, LAZ, COPC (Cloud Optimized Point Cloud).
  * **Vector:** Shapefile, GeoJSON, OGR-compatible formats.
  * **Text:** XYZ, CSV, DAT.
  * **Remote:** URLs pointing to data fetched via the `fetches` module.


* **Standardized Output:** regardless of the input format, Dlim processes and yields data in a consistent structure (typically XYZ or array chunks), abstracting away file-specific parsing logic from downstream tools like `waffles`.

### Advanced Data Processing

Dlim applies a suite of operations to data on-the-fly as it is read:

* **Spatial Filtering:**
* **Region Clipping:** Automatically filters data to a specified bounding box (`-R`) or vector polygon, ensuring only relevant data is processed.
* **Masking:** Can apply complex polygon masks (e.g., land/water masks) to exclude unwanted areas.


* **Vertical Transformation:**
  * **Datum Conversion:** Integrates with `vdatums` to perform vertical datum transformations (e.g., NAVD88 to MLLW) during ingestion.
  * **Z-Scaling/Shifting:** Allows applying scalar weights, offsets, or unit conversions (e.g., meters to feet) to elevation values.


* **Outlier Detection:** Can pre-filter noise or outliers from point clouds before gridding, improving the quality of the final DEM.

### Data List Management (Recursive Processing)

* **Hierarchical DataLists:** `dlim` supports nested datalists, where a principal list can reference other lists. This allows for organized, multi-scale data management (e.g., a "Global" list pointing to "Regional" lists).
* **Weighted Ranking:** Data entries can be assigned weights. When datasets overlap, `dlim` (in conjunction with `waffles`) uses these weights to determine which data source takes precedence, allowing high-quality modern surveys to supercede older, lower-resolution data.
* **Archive Generation:** The module can package referenced data into a portable archive, facilitating data sharing and reproducibility.

### Integration

`dlim` serves as the bridge between **`fetches`** and **`waffles`**. It can execute `fetches` modules to retrieve remote data dynamically and then pipe that data directly into `waffles` for gridding. This enables "dataless" workflows where local storage of massive source datasets is not required; data is streamed, processed, and gridded in a single pipeline.

## Usage Example

### A typical use case involves creating a datalist file (`my_project.datalist`) that references a local lidar file and a remote NOAA survey:

```text
# my_project.datalist
/data/local/lidar_2020.laz -1 1
https://fetch_noaa_survey_123.laz 10 1

```

The `dlim` CLI can then process this list to inspect, filter, or export the data:

#### Archive Data: Gather all referenced data in a datalist into a portable archive:

```bash
dlim my_project.datalist -R -90/-89/29/30 --archive my_project_data.tar.gz

```

### Other Examples

#### Process fetches dataset `hydronos` to input region and output to an xyz file, including weights and uncertainty

```bash
dlim -R-119.25/-119/34/34.25 hydronos -w -u > my_hydronos.xyz
```

#### Process all local `*.laz` lidar files, transform them from UTM zone 10N to WGS84 and block them to the input region at 1 arc-second increments

```bash
dlim -R-123.458/-123.448/38.705/38.711 -E1s -J epsg:26910 -P epsg:4326 *.laz > my_lidar_1as.xyz
```

#### Generate a data-mask (raster) and spatial-metadata (vector) of a datalist:

```bash
dlim -Rmy_region.shp -E1s -Pepsg:4269 --mask --spatial-metadata my_data.datalist > /dev/null
```