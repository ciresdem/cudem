# Fetches (Data Acquisition)

**Fetches** is the data acquisition and download module within the CUDEM software suite. It automates the process of identifying, locating, and retrieving topographic and bathymetric source data from over forty different publicly available datasets.

Designed to streamline the initial phase of DEM generation, Fetches allows users to target specific geographic regions and automatically pull relevant data without manual searching or downloading.

## Summary

### Automated Data Retrieval

* **Region-Based fetching:** Users specify a geographic extent (bounding box or vector polygon), and the module identifies all available data tiles or granules from supported services that intersect that region.
* **Smart Buffering:** It typically downloads data in an area slightly larger (~5%) than the requested extent to prevent edge effects during interpolation.
* **Protocol Handling:** Manages various remote access protocols (HTTP, FTP, S3, APIs) to retrieve files seamlessly.

### On-the-Fly Processing & Parsing

The module does not just download files; it acts as an interface to standardize diverse data formats for the CUDEM pipeline. Through the `Fetcher` classes defined in `datalists.fetchers`:

* **Format Conversion:** It can parse complex formats (e.g., BAG, HDF5, NetCDF) and yield them through `dlim` as standardized datasets (XYZ or Raster) for processing.
* **Metadata Extraction:** It extracts critical metadata such as horizontal/vertical datums, resolution, and collection dates from the source files.
* **Masking:** Can automatically apply coastline or water masks to global grids (e.g., masking land in bathymetry grids).
* **Filtering:** Can pre-filter point clouds (via `dlim`, `pointz` and `grits`)(e.g., removing specific classifications from ICESat-2 data) before they enter the gridding pipeline.

### Modular & Extensible

Fetches uses a factory system, allowing specific modules to be written for different data providers. If a dataset requires special API calls or post-download processing (like unzipping or converting datums), a dedicated datalists `Fetcher` subclass handles it.

## Supported Data Sources

The module supports a wide array of global and regional datasets, including but not limited to:

* **NOAA (National Oceanic and Atmospheric Administration):**
  * **NOS Hydrographic Surveys:** Bathymetric sounding data (`HydroNOS`).
  * **Digital Coast:** CoNED Topobathy and Sea Level Rise (SLR) DEMs.
  * **Multibeam:** Raw and processed swath sonar data (`MBS`).
  * **BlueTopo:** High-resolution target detection bathymetry.
  * **Electronic Navigational Charts (ENC):** Digital soundings and contours (`Charts`).
  * **Geodesy:** NGS Monuments (`NGS`) and VDatum grids.


* **USGS (United States Geological Survey):**
  * **The National Map (TNM):** National Elevation Dataset (NED/3DEP).
  * **Water Services:** River and stream gauge data.


* **NASA (National Aeronautics and Space Administration):**
  * **ICESat-2:** Satellite laser altimetry (ATL03/ATL24).
  * **SWOT:** Surface Water and Ocean Topography data.


* **Global & Regional Grids:**
  * **GEBCO:** General Bathymetric Chart of the Oceans.
  * **GMRT:** Global Multi-Resolution Topography.
  * **Copernicus:** European global DEM.
  * **FABDEM:** Forest And Buildings removed Copernicus DEM.
  * **EMODnet:** European Marine Observation and Data Network.


* **Other Sources:**
  * **Crowd Sourced Bathymetry (CSB):** Citizen science depth data.
  * **USACE:** eHydro hydrographic surveys.
  * **MarGrav:** Satellite-derived marine gravity bathymetry.



## Integration

While `fetches` can be run as a standalone command-line tool to download files, it is tightly integrated with the **`dlim`** module. `dlim` can directly call `fetches` modules to stream remote data into the **`waffles`** gridding engine without creating intermediate local files, enabling efficient end-to-end processing workflows.


## Python API

## Examples

### Fetch IceSat2 data with the fetches Python API

```python
from cudem import regions
from cudem import fetches
r = regions.Region().from_list([-121.25, -121.0, 37.5, 37.75])
f = fetches.IceSat2(subset=False, time_start='2020-05-04', time_end='2021-05-04', src_region=r).run()
fr = fetches.fetch_results(f)
fr.daemon=True
fr.start()
fr.join()
```

### Fetch IceSat2 data as a subset with the fetches CLI

```bash
$ fetches -R -121.25/-121/37.5/37.75 icesat2:subset=True:time_start=2020-05-04:time_end=2021-05-04
```

### List the available HydroNOS data files in a given region

```bash
$ fetches -R-119.25/-119/34/34.25 hydronos -l
fetches: parsed 1 region(s): [-119.25 -119.0 34.0 34.25 None None None None None None]
fetches: running fetch module hydronos on region -119.25/-119.0/34.0/34.25...
fetches: found 34 data files.
https://data.ngdc.noaa.gov/platforms/ocean/nos/coast/H04001-H06000/H05392/GEODAS/H05392.xyz.gz
https://data.ngdc.noaa.gov/platforms/ocean/nos/coast/H04001-H06000/H05420/GEODAS/H05420.xyz.gz
https://data.ngdc.noaa.gov/platforms/ocean/nos/coast/H04001-H06000/H05425/GEODAS/H05425.xyz.gz
https://data.ngdc.noaa.gov/platforms/ocean/nos/coast/H04001-H06000/H05426/GEODAS/H05426.xyz.gz
https://data.ngdc.noaa.gov/platforms/ocean/nos/coast/H04001-H06000/H05507/GEODAS/H05507.xyz.gz
https://data.ngdc.noaa.gov/platforms/ocean/nos/coast/H04001-H06000/H05446/GEODAS/H05446.xyz.gz
https://data.ngdc.noaa.gov/platforms/ocean/nos/coast/H04001-H06000/H05851/GEODAS/H05851.xyz.gz
https://data.ngdc.noaa.gov/platforms/ocean/nos/coast/H08001-H10000/H09600/GEODAS/H09600.xyz.gz
https://data.ngdc.noaa.gov/platforms/ocean/nos/coast/H08001-H10000/H09666/GEODAS/H09666.xyz.gz
https://data.ngdc.noaa.gov/platforms/ocean/nos/coast/H08001-H10000/H09667/GEODAS/H09667.xyz.gz
https://data.ngdc.noaa.gov/platforms/ocean/nos/coast/H08001-H10000/H09725/GEODAS/H09725.xyz.gz
https://data.ngdc.noaa.gov/platforms/ocean/nos/coast/H10001-H12000/H11501/BAG/H11501_MB_1m_MLLW_2of4.bag
https://data.ngdc.noaa.gov/platforms/ocean/nos/coast/H10001-H12000/H11501/BAG/H11501_MB_4m_MLLW_3of4.bag
https://data.ngdc.noaa.gov/platforms/ocean/nos/coast/H10001-H12000/H11501/BAG/H11501_MB_50cm_MLLW_1of4.bag
https://data.ngdc.noaa.gov/platforms/ocean/nos/coast/H10001-H12000/H11501/BAG/H11501_MB_8m_MLLW_4of4.bag
https://data.ngdc.noaa.gov/platforms/ocean/nos/coast/H10001-H12000/H11501/GEODAS/H11501.xyz.gz
https://data.ngdc.noaa.gov/platforms/ocean/nos/coast/H10001-H12000/H11891/BAG/H11891_MB_16m_MLLW_5of6.bag
https://data.ngdc.noaa.gov/platforms/ocean/nos/coast/H10001-H12000/H11891/BAG/H11891_MB_1m_MLLW_1of6.bag
https://data.ngdc.noaa.gov/platforms/ocean/nos/coast/H10001-H12000/H11891/BAG/H11891_MB_2m_MLLW_2of6.bag
https://data.ngdc.noaa.gov/platforms/ocean/nos/coast/H10001-H12000/H11891/BAG/H11891_MB_32m_MLLW_6of6.bag
https://data.ngdc.noaa.gov/platforms/ocean/nos/coast/H10001-H12000/H11891/BAG/H11891_MB_4m_MLLW_3of6.bag
https://data.ngdc.noaa.gov/platforms/ocean/nos/coast/H10001-H12000/H11891/BAG/H11891_MB_8m_MLLW_4of6.bag
https://data.ngdc.noaa.gov/platforms/ocean/nos/coast/W00001-W02000/W00431/BAG/W00431_MB_128m_MLLW_5of5.bag
https://data.ngdc.noaa.gov/platforms/ocean/nos/coast/W00001-W02000/W00431/BAG/W00431_MB_128m_MLLW_Combined.bag
https://data.ngdc.noaa.gov/platforms/ocean/nos/coast/W00001-W02000/W00431/BAG/W00431_MB_16m_MLLW_2of5.bag
https://data.ngdc.noaa.gov/platforms/ocean/nos/coast/W00001-W02000/W00431/BAG/W00431_MB_32m_MLLW_3of5.bag
https://data.ngdc.noaa.gov/platforms/ocean/nos/coast/W00001-W02000/W00431/BAG/W00431_MB_64m_MLLW_4of5.bag
https://data.ngdc.noaa.gov/platforms/ocean/nos/coast/W00001-W02000/W00431/BAG/W00431_MB_8m_MLLW_1of5.bag
https://data.ngdc.noaa.gov/platforms/ocean/nos/coast/W00001-W02000/W00432/BAG/W00432_MB_128m_MLLW_5of5.bag
https://data.ngdc.noaa.gov/platforms/ocean/nos/coast/W00001-W02000/W00432/BAG/W00432_MB_128m_MLLW_Combined.bag
https://data.ngdc.noaa.gov/platforms/ocean/nos/coast/W00001-W02000/W00432/BAG/W00432_MB_16m_MLLW_2of5.bag
https://data.ngdc.noaa.gov/platforms/ocean/nos/coast/W00001-W02000/W00432/BAG/W00432_MB_32m_MLLW_3of5.bag
https://data.ngdc.noaa.gov/platforms/ocean/nos/coast/W00001-W02000/W00432/BAG/W00432_MB_64m_MLLW_4of5.bag
https://data.ngdc.noaa.gov/platforms/ocean/nos/coast/W00001-W02000/W00432/BAG/W00432_MB_8m_MLLW_1of5.bag

```

### Generate a Datalist: Instead of downloading immediately, generate a datalist of URLs for later processing:

```bash

fetches -R -90/-89/28/29 -l nos:datatype=bag > my_survey_data.datalist

```

### Fetch Multiple Sources:
```bash

fetches -R -90/-89/28/29 'noaa_lidar' 'usace' '3dep'

```