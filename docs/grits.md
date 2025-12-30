# Grits (GRId filTerS)

**Grits** is the elevation grid processing and filtering engine within the CUDEM software suite. It provides a standardized framework for manipulating raster datasets, allowing users to clean artifacts, smooth noise, blend overlapping datasets, and enforce hydrological or morphological constraints.

Designed to handle massive datasets, Grits operates on a chunk-by-chunk basis (tiling), ensuring that operations on multi-gigabyte DEMs do not exceed system memory.

## Core Capabilities

### 1. Framework Features

* **Chunked Processing:** Automatically handles large files by processing them in windowed chunks (with edge buffering) to prevent memory overflows.
* **Auxiliary Data Loading:** Can load and align arbitrary auxiliary rasters (e.g., weight masks, uncertainty grids, count grids) to guide the filtering process.
* **Standardized Masking:** All filters support "Split" operations, allowing users to apply filters only to pixels within specific Z, Weight, or Uncertainty ranges.
* **Chaining:** Filters can be chained together in the CLI (e.g., `-M outliers -M fill -M blur`) to create complex processing pipelines.

### 2. Smoothing & Denoising

* **`blur`:** Applies a **Gaussian Blur** to smooth the entire DEM. Useful for general noise reduction.
* **`denoise`:**
* **Median Filter:** Removes "salt-and-pepper" noise (single-pixel spikes or pits) while preserving sharp edges like cliffs.
* **Bilateral Filter:** A non-linear filter that smooths flat areas but preserves strong edges (requires `scikit-image`).


* **`gmtfilter`:** A wrapper for the Generic Mapping Tools (GMT) `grdfilter` module, providing access to robust geodetic filtering (if GMT/PyGMT is installed).

### 3. Artifact & Outlier Removal

* **`outliers` (LSPOutliers):** A sophisticated statistical filter. It calculates derived Land Surface Parameters (Slope, TPI, TRI, Roughness) and uses Tukey's Fences (IQR) to identify and remove statistical outliers in terrain complexity. Supports multi-pass scanning at different scales.
* **`zscore`:** Detects local anomalies by calculating the Z-Score of a pixel relative to its immediate neighborhood. Excellent for finding subtle spikes or pits that don't violate global min/max thresholds.
* **`flats`:** Identifies and removes artificial "terracing" or flat plateaus often caused by integer-precision data or interpolation artifacts.

### 4. Blending & Integration

* **`blend`:** Smooths the transition between a source DEM and auxiliary data (or high-weight foreground data). It creates a buffer zone where values are interpolated to seamlessly stitch datasets together.
* **`weights`:** Buffers around "high-quality" (high weight) data pixels and removes lower-quality pixels found within that buffer. This cleans up the "halos" of bad data often found at the edges of swath sonar coverage.
* **`diff`:** Calculates the difference between the source DEM and a reference grid. Can output the difference grid or mask pixels where the change exceeds a threshold.

### 5. Morphology & Hydrology

* **`morphology`:** Applies grayscale morphological operations:
* **Erosion:** Widens valleys/channels and removes small peaks (useful for removing vegetation/buildings).
* **Dilation:** Widens peaks/ridges and fills small pits.
* **Opening/Closing:** Combinations used to remove specific noise features while preserving overall shape.


* **`hydro`:** Provides basic hydrological enforcement, such as **Sink Filling**, to remove local depressions and ensure continuous flow across the surface.

### 6. Masking & Geometry

* **`cut`:** Masks the DEM to a specific bounding box (Region). Supports inversion (hole punching).
* **`clip`:** Clips the DEM to an OGR-compatible vector polygon (Shapefile, GeoJSON).
* **`slope` (SlopeFilter):** Masks or reverts pixels based on derived terrain metrics (e.g., "Remove all data where Slope > 45 degrees").

### 7. Void Filling

* **`fill`:** Fills NoData voids (inpainting) using:
* **IDW (Inverse Distance Weighting):** Fast filling for small gaps.
* **Spline/Linear Interpolation:** Fits a smooth surface across larger voids to maintain terrain trends.



## Usage Example

The `grits` CLI allows you to chain these modules. For example, to clean a noisy DEM by removing outliers, filling the resulting holes, and then clipping it to a coastline vector:

```bash
grits input_dem.tif \
  -M outliers:k=2.5:aggressive=True \
  -M fill:method=spline:max_dist=50 \
  -M clip:src_ply=coastline.shp \
  -O cleaned_output.tif

```