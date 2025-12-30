# Regions (Spatial Management)

**Regions** is the spatial management module within the CUDEM software suite. It provides a robust `Region` class and associated utilities to handle geographic bounding boxes, coordinate systems, and tiling logic. This module is fundamental to CUDEM's chunked processing architecture, defining the spatial extent of operations for modules like `fetches`, `waffles`, and `grits`.

## Core Capabilities

### 1. Unified Region Representation

The `Region` class standardizes spatial bounds across the suite. It stores more than just X/Y coordinates; it tracks data range (Z), weights (W), and uncertainty (U) limits:

* **Spatial Bounds:** `xmin`, `xmax`, `ymin`, `ymax`.
* **Data Limits:** `zmin`/`zmax` (elevation), `wmin`/`wmax` (weight), `umin`/`umax` (uncertainty).
* **Projection Awareness:** Tracks the source SRS (`src_srs`) and supports warping regions between coordinate systems (e.g., transforming a bounding box from NAD83 to WGS84) via `pyproj`.

### 2. Flexible Input Parsing

The module can parse region definitions from a wide variety of formats:

* **Strings:** GMT-style strings (`-R-90/-89/28/29`), space-separated lists (`xmin xmax ymin ymax`), or even coordinate strings.
* **Vector Files:** Can read OGR-compatible vector files (Shapefile, GeoJSON) and extract the bounding box of every feature within the file.
* **Geotransforms:** Can derive regions from GDAL geotransform arrays and pixel dimensions.

### 3. Manipulation & Tiling

Regions includes powerful tools for manipulating spatial extents:

* **Buffering:** Expand regions by a fixed unit value or a percentage.
* **Chunking/Tiling:** Decompose a large region into smaller, manageable tiles (e.g., generate a set of 0.25-degree tiles covering a larger area).
* **Intersection/Union:** Calculate the intersection (overlap) or union (merge) of multiple regions.
* **Cutting:** "Cut" one region out of another, handling grid alignment to ensuring pixels remain snapped to the grid.

### 4. Output Formatting

Regions can be exported to numerous standard formats for interoperability with other tools:

* **Strings:** GMT format (`-R...`), PROJ (`-te...`), or filename-safe strings (e.g., `n40x00_w105x00`).
* **Vectors:** Export region extents directly to OGR vector files (GeoJSON, Shapefile) for visualization in GIS software.
* **WKT:** Well-Known Text polygon representation.

## Usage Example

The `regions` CLI tool allows for quick manipulation of bounding boxes. For example, to buffer a region by 5% and output the result in GMT format:

```bash
regions -R -90/-89/28/29:pct_buffer=5 -e

```

Or to generate a 0.25-degree tile set covering a larger area:

```bash
regions -R -90/-88/28/30 -T 0.25

```