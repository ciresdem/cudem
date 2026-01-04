# Waffles (Gridding & Interpolation)

**Waffles** is the DEM generation and interpolation engine within the CUDEM software suite. It is designed to ingest scattered elevation data (XYZ, LAS/LAZ, Raster) and convert them into integrated, regularly spaced Digital Elevation Models (DEMs) using a variety of algorithms.

The module operates on a "stacking" principle, where multiple input datasets are first combined into a weighted mean surface (the "stack") before being passed to an interpolation module to fill gaps and generate the final surface.

## Summary

### Flexible Interpolation Modules

`waffles` utilizes a factory system to support numerous gridding algorithms, allowing users to choose the best method for their specific terrain and data density, and is extensible via a factory pattern:

* **`stacks`**: Simple weighted averaging of data (no interpolation). Useful for compositing existing grids.
* **GMT Wrappers**: `gmt-surface` (continuous curvature splines), `gmt-triangulate` (Delaunay), `gmt-nearneighbor` (nearest neighbor averaging).
* **Scipy Wrappers**: `linear` (TIN), `cubic` (Spline), `nearest` (Nearest Neighbor).
* **Advanced Methods** `cudem` (Multi-resolution integration), _testing_: `kriging` (Geostats), `natural_neighbor` (Sibson), `ml-interp` (Machine Learning), `inpaint` (Void Filling).
* **Others**: `IDW` (Inverse Distance Weighting), `mbgrid` (MB-System wrapper), `gdal-linear`, `gdal-nearest`, `gdal-average`, `gdal-invdst` (GDAL grid wrappers).

### Smart Data Handling (via `dlim`)

* **Accepts diverse inputs**: ASCII XYZ, LAS/LAZ, GeoTIFF, BAG, OGR vectors, MB-System datalists, and fetch modules.
* **Stacking**: Before interpolation, data is "stacked" into a weighted intermediate raster. This handles overlapping datasets by calculating weighted means or allowing high-quality data to supersede lower-quality data.

### Region & Resolution Management

* **Chunking**: Capable of processing massive datasets by splitting the region into smaller chunks (-K), processing them in parallel, and stitching the results back together.
* **Buffering/Extension**: Automatically buffers regions during processing (-X) to prevent edge artifacts, then crops the result to the desired extent.
* **Increments**: Supports independent input processing resolution (-E xinc) and output sampling resolution (xsample).

### Advanced Post-Processing

Once a raw DEM is generated, `waffles` (via the `WaffleDEM` class) performs extensive post-processing to ensure quality:

* **Clipping & Cutting**: Clips the output to vector polygons (e.g., coastlines) or precise bounding boxes.
* **Filtering**: Applies `grits` filters (e.g., smoothing, outlier removal) to the generated grid.
* **Limiting**: Constrains interpolation based on proximity to valid data or the size of data gaps, preventing artifacts in sparse areas.
* **Uncertainty**: Can generate a corresponding uncertainty grid, estimating the error of the interpolation.

### Output Management

* **Formats**: Supports exporting to GeoTIFF, NetCDF, and HDF5.
* **Metadata**: Automatically injects spatial metadata (ISO 19115 tags) into the output files.
* **Masking**: Generates validity masks to distinguish interpolated areas from observed data.

## Usage Example

The `waffles` CLI is used to generate a DEM. For example, to generate a DEM using GMT's surface spline algorithm, clipped to a coastline:

```bash
waffles my_data.datalist -R -90/-89/28/29 -E .111111111s -M gmt-surface:tension=0.35 -C coast.shp -O output_dem

```

**Common Options:**

* `-M, --module`: The gridding algorithm to use (e.g., `surface`, `linear`, `IDW`).
* `-R, --region`: The output bounding box.
* `-E, --increment`: The grid resolution.
* `-w, --want-weight`: Use weights from the input datalist.
* `-u, --want-uncertainty`: Generate an uncertainty layer.