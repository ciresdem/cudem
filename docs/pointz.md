# Pointz (Point Cloud Filtering)

**Pointz** is a specialized module within the CUDEM software suite for filtering and manipulating raw point cloud data (XYZ, LAS/LAZ, etc.) before it is gridded into a DEM. While `grits` operates on raster grids, `pointz` operates directly on vector point data, allowing for precise cleaning of source data at the individual sounding level.

The module provides tools for statistical outlier removal, reference-based quality checking, and spatial masking using vector polygons. It serves as a pre-processing step to ensure that only valid, high-quality measurements enter the interpolation pipeline.

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

### 3. Spatial Masking (`vector_mask`)

* **Polygon Filtering:** The `PointZVectorMask` class filters points based on their inclusion within vector polygons (e.g., Shapefiles).
* **Inversion Support:** Supports both "keep inside" (e.g., bounding box) and "remove inside" (e.g., land masking for bathymetry) operations via the `invert` parameter.

### 4. Point Gridding Utilities

* **`PointPixels` Class:** A helper utility that bins scattered points into a regular grid structure, calculating aggregated statistics (mean, min, max, sum, uncertainty) for each cell. This is the foundational logic used by the outlier filters to establish local statistics.

## Usage

The `pointz` module is accessible via a CLI tool that takes an input point file, applies a chain of filters, and outputs the cleaned data.

```bash
pointz input.xyz output.xyz -M outlierz:percentile=98:res=50 -M vector_mask:mask_fn=land_poly.shp:invert=True

```

**Common Options:**

* `-M, --module`: Select the filter module and parameters (e.g., `outlierz:percentile=95`).
* `outlierz` parameters: `percentile` (threshold), `res` (block resolution), `multipass` (number of iterations).
* `rq` parameters: `threshold` (difference limit), `raster` (path to reference grid).
* `vector_mask` parameters: `mask_fn` (path to vector file), `invert` (bool).