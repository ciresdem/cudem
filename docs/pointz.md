# Pointz (Point Cloud Filtering)

**Pointz** is a specialized module within the CUDEM software suite for filtering and manipulating raw point cloud data (XYZ, LAS/LAZ, etc.) before it is gridded into a DEM. While `grits` operates on raster grids, `pointz` operates directly on vector point data, allowing for precise cleaning of source data at the individual sounding level.

The module provides tools for statistical outlier removal, reference-based quality checking, spatial masking using vector and raster data, and density-based thinning. It serves as a pre-processing step to ensure that only valid, high-quality measurements enter the interpolation pipeline.

## Core Capabilities

### 1. Statistical Outlier Detection (`outlierz`)

* **Block-Based Statistics:** The `PointZOutlier` class divides the point cloud into local blocks (spatial bins) and calculates statistics (mean, standard deviation) for each block.
* **Residual Analysis:** It computes the residual of each point relative to the local block mean.
* **Thresholding:** Points with residuals exceeding a user-defined percentile (e.g., 95th or 99th percentile) are flagged as outliers. This is effective for removing noise spikes or gross errors that deviate significantly from the local trend.
* **Multi-Pass Filtering:** Supports iterative filtering at multiple resolutions (`multipass`), gradually refining the outlier detection from coarse to fine scales to catch different types of errors.

### 2. Reference-Based Quality Control (`rq`)

* **External Reference Validation:** The `RQOutlierZ` class compares input points against a trusted reference surface (e.g., a low-resolution regional grid like GMRT, ETOPO, or a previous CUDEM version).
* **Difference Thresholding:** Points that deviate from the reference surface by more than a specified threshold (absolute difference or percentage) are flagged. This is useful for identifying gross vertical errors or datum shifts in new datasets.
* **Dynamic Fetching:** Can automatically fetch and use standard reference datasets (GMRT, CUDEM) if a local reference raster is not provided.

### 3. Spatial Masking

* **Vector Masking (`vector_mask`):** The `PointZVectorMask` class filters points based on their inclusion within vector polygons (e.g., Shapefiles). It supports both "keep inside" and "remove inside" operations via the `invert` parameter.
* **Raster Masking (`raster_mask`):** The `PointZRasterMask` class filters points based on the values of an underlying raster file (e.g., Land/Water mask). It keeps points on valid (non-zero) pixels by default.

### 4. Thinning and Density Control

* **Grid-Based Decimation (`block_minmax`):** A highly optimized filter that divides the domain into grid cells and retains only the minimum (shoalest) or maximum (deepest) point per cell. This is critical for hydrographic applications requiring shoal-biased thinning.
* **Density Control (`density`):** Thins the point cloud to a specific resolution using various selection modes (`random`, `median`, `mean`, `center`).

### 5. Advanced Geometric Filtering

* **Range Filtering (`range`):** Strictly clips points to a specific vertical Z-range (e.g., removing all points above sea level).
* **Plane Fitting (`coplanar`):** Identifies points that deviate from a locally fitted plane, useful for cleaning noise from generally flat surfaces like roads or the seafloor.
* **Difference Filtering (`diff`):** Calculates the signed difference between points and a reference surface, useful for bias detection or removing points that are specifically "above" or "below" a reference.

### 6. Point Gridding Utilities

* **`PointPixels` Class:** A helper utility that bins scattered points into a regular grid structure, calculating aggregated statistics (mean, min, max, sum, uncertainty) for each cell. This is the foundational logic used by the outlier filters to establish local statistics.

## Usage and Examples

The `pointz` module is accessible via a CLI tool that takes an input point file, applies a chain of filters, and outputs the cleaned data.

### Remove outliers and mask out water

```bash
pointz input.xyz output.xyz -M outlierz:percentile=98:res=50 -M vector_mask:mask_fn=land_poly.shp:invert=True

```

### Basic Cleaning: Remove statistical outliers from a noisy dataset.

```bash

pointz raw.xyz clean.xyz -M outlierz:percentile=95:res=50

```

### Hydrographic Thinning: Shoal-bias thin a dataset to 10m resolution.

```bash

pointz dense.laz thinned.xyz -M block_minmax:res=10:mode=min

```

### Complex Workflow: Clip to a range, mask out land, and remove outliers.

```bash

pointz input.xyz output.xyz \
  -M range:max_z=0 \
  -M vector_mask:mask_fn=land.shp:invert=True \
  -M outlierz:percentile=98

```

**Common Options:**

* `-M, --module`: Select the filter module and parameters (e.g., `outlierz:percentile=95`).
* `outlierz` parameters: `percentile` (threshold), `res` (block resolution), `multipass` (number of iterations).
* `rq` parameters: `threshold` (difference limit), `raster` (path to reference grid).
* `vector_mask` parameters: `mask_fn` (path to vector file), `invert` (bool).
* `block_minmax` parameters: `res` (cell size), `mode` (min/max).