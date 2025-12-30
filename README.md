Matthew Love[^1][^2], Christopher Amante[^1][^2], Kelly Carignan[^1][^2], Elliot Lim[^1][^2], Michael MacFerrin[^1][^2]

[^1]: Cooperative Institute for Research in Environmental Sciences (CIRES)
at the University of Colorado Boulder

[^2]: National Oceanic and Atmospheric Administration (NOAA) National
Centers for Environmental Information (NCEI)

# Overview

The National Oceanic and Atmospheric Administration (NOAA) National
Centers for Environmental Information (NCEI), through its collaboration
with the Cooperative Institute for Research in Environmental Sciences
(CIRES) at the University of Colorado Boulder, develops digital
elevation models (DEMs) that range from the local to global scale.
Collectively, these elevation models are essential to determining the
timing and extent of coastal inundation and improving community
preparedness, event forecasting, and warning systems. We initiated a
comprehensive framework at NCEI, the Continuously-Updated DEM (CUDEM)
Program, to systematically generate DEMs from the local coastal
community to the global scale.

We generate the CUDEMs through a standardized process using free and
open-source software (FOSS) and provide open-access to our code
repository
([https://github.com/ciresdem](https://github.com/ciresdem)) for
consistency, transparency, and to promote accessibility. The CUDEM
framework consists of systematic tiled geographic extents, spatial
resolutions, and horizontal and vertical datums to facilitate rapid
updates of targeted areas with new data collections, especially
post-storm and tsunami events. The CUDEM Program is also enabling the
rapid incorporation of high-resolution data collections ingested into
local-scale DEMs into NOAA NCEI's suite of regional and global DEMs. The
CUDEMs are a shift from project-based DEM specifications, to a
comprehensive program that systematically and continuously develops and
updates DEMs across all spatial scales.

## CUDEM Framework

The CUDEM framework provides a set of powerful open-source command-line tools as well as a Python-3-based application programming interface (API) for:
1.  **Fetching**: Identifying and downloading topographic and bathymetric source data from more than forty different publicly available datasets (e.g., NOAA, USGS, USACE).
2.  **Processing**: Filtering, cleaning, and masking source data (point clouds, rasters, vectors) to remove artifacts and enforce quality standards.
3.  **Gridding**: Ranking and merging multiple datasets using a variety of interpolation algorithms to generate seamless Digital Elevation Models (DEMs).
4.  **Analysis**: Generating auxiliary data products including uncertainty grids, data masks, and spatial metadata.
5.  **Transformation**: converting data between various horizontal and vertical datums.

Individual components are modular and may be run standalone; for example, the `fetches` module can download airborne lidar data or nautical charts covering specific geographic regions; the `waffles` module can generate DEMs from scattered data; and the `grits` module can filter and post-process those DEMs. The CUDEM tools may also be conjoined in end-to-end workflows ensuring permanent traceability and reproducibility in DEM generation.

# Installation and Setup

## Requirements
* **Python 3.8+**
* **GDAL** (and Python bindings) is required.
* **Numpy, Scipy**

Other useful external programs needed for full functionality include:
* **GMT** / **PyGMT**: For robust gridding and filtering options.
* **MB-System**: For processing raw bathymetric sonar formats.
* **VDatum**: For vertical datum transformations.
* **PDAL**: For processing point cloud formats (LAS/LAZ/COPC).

## Install via Pip
Download and install git (If you have not already): [git installation](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git)

```bash
pip install git+[https://github.com/ciresdem/cudem.git#egg=cudem](https://github.com/ciresdem/cudem.git#egg=cudem)
