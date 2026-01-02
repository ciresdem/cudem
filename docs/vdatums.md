# Vdatums (Vertical Datums)

**Vdatums** is the vertical datum transformation module within the CUDEM software suite. It provides a comprehensive framework for transforming elevation data between different vertical reference frames, a critical step when integrating datasets from diverse sources (e.g., bathymetry in Mean Lower Low Water vs. topography in NAVD88).

Vdatums generates vertical transformation grids that represent the offset between two specified vertical datums. These grids are then added to or subtracted from the source DEM to shift it into the desired target datum.

## Summary

### Extensive Datum Support

Vdatums supports transformations between three major categories of vertical reference frames:

* **Tidal Datums:** Supports transformations involving MLLW, MHHW, MHW, MSL, MTL, and others (EPSG codes 1089, 5866, 1091, 5869, 5868, 5714, 5713). It leverages NOAA's VDatum software to generate the necessary separation grids.
* **Ellipsoidal Datums:** Supports transformations between various geometric reference frames used by GPS/GNSS, such as NAD83 variants (2011, PA11, MA11), WGS84 realizations (G730, G873, G1150, G1674, G1762, G2139), and ITRF realizations (ITRF88-2020). It utilizes the **HTDP** (Horizontal Time-Dependent Positioning) utility for these geodetic shifts.
* **Orthometric Datums:** Supports transformations involving gravity-based datums like NAVD88, NGVD29, EGM96, EGM2008, CGVD2013, and PRVD02. It accesses the **CDN** (Cumulative Distribution Network) or internal geoid grids (e.g., GEOID12B, GEOID18) to perform these conversions.

### Grid-Based Transformation Generation

Rather than transforming points one by one, Vdatums generates a spatially varying transformation grid (raster) that covers the extent of the input data. This ensures that local variations in the separation between datums (e.g., the geoid slope or tidal variance) are accurately captured across the entire dataset.

### Integrated Workflow

* **Seamless Chaining:** The module intelligently chains transformations. If a direct conversion doesn't exist (e.g., NAVD88 to MLLW), it routes the transformation through intermediate steps (e.g., NAVD88 -> NAD83 -> MLLW) automatically.
* **Uncertainty Estimation:** Vdatums can estimate the uncertainty introduced by the transformation process, providing an error grid alongside the transformation grid.

### External Tool Integration

Vdatums acts as a wrapper and orchestrator for industry-standard geodesy tools:

* **NOAA VDatum:** Used for U.S. tidal transformations.
* **HTDP:** Used for rigorous crustal motion and reference frame transformations.
* **PROJ/GDAL:** Used for grid handling and certain projection-based shifts.

## Usage Example

The `vdatums` CLI can be used to transform a DEM from NAVD88 (EPSG:5703) to WGS84 (G1674) Ellipsoid Height (EPSG:7662):

```bash
vdatums input_navd88.tif output_wgs84.tif -i 5703 -o 7662

```

This command will:

* Identify the transformation path (NAVD88 -> NAD83 -> WGS84).
* Fetch or generate the necessary separation grids (Geoid18, etc.).
* Compute the composite offset grid.
* Apply the offset to `input_navd88.tif` to create `output_wgs84.tif`.