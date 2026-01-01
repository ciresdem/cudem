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

```

## Install via Conda

It is recommended to use a conda environment to manage GDAL and other geospatial dependencies.

```bash
conda create -n cudem -c conda-forge gdal gmt pygmt numpy scipy pandas pyproj utm requests lxml matplotlib laspy h5py boto3 tqdm mercantile git netCDF4 h5netcdf libgdal-hdf5 libgdal-netcdf pyhdf pip scikit-image
conda activate cudem
pip install laspy[laszip]
pip install --upgrade --no-deps git+[https://github.com/ciresdem/cudem.git](https://github.com/ciresdem/cudem.git)

```

## Extras

* **Windows**: Install `windows-curses` via pip: `pip install windows-curses`
* **MB-System**: [Installation Instructions](https://www.mbari.org/technology/mb-system/installation/)
* **HTDP / VDatum**: See [NOAA Geodesy Tools](https://geodesy.noaa.gov/TOOLS/Htdp/Htdp.shtml).

# Modules

The CUDEM suite is composed of several key modules, each accessible via command-line tools or Python API:

### CUDEM Modules

| Module | Description | CLI Command |
| --- | --- | --- |
| **[fetches](/docs/fetches.md)** | **Data Acquisition**: Fetch elevation data from a wide variety of public online sources (NOAA, USGS, USACE, etc.). | `fetches` |
| **[pointz](/docs/pointz.md)** | **Point-Cloud Filtering**: A powerful filtering engine to filter Point Clouds. | `pointz` |
| **[grits](/docs/grits.md)** | **Grid Filtering & Post-Processing**: A powerful filtering engine to smooth, clean, blend, and analyze raster DEMs. | `grits` |
| **[dlim](/docs/dlim.md)** | **Data Lists & Processing**: Process, clean, and standardize diverse data types (XYZ, Raster, LAS, COPC, etc.) into unified formats. | `dlim` |
| **[waffles](/docs/waffles.md)** | **Gridding & interpolation**: Generate Digital Elevation Models (DEMs) from scattered data using various algorithms (IDW, Spline, GMT Surface, etc.). | `waffles` |
| **[regions](/docs/regions.md)** | **Spatial Management**: Process and manipulate bounding box regions and vector polygons. | `regions` |
| **[vdatums](docs/vdatums.md)** | **Vertical Transformation**: Generate vertical transformation grids and manage datums. | `vdatums` |
| **[perspecto](/docs/perspecto.md)** | **Visualization**: Generate perspective images, hillshades, and color-relief maps of DEMs. | `perspecto` |

---

## Module Synopsis:

### Fetches (Data Acquisition)

The **[fetches](/docs/fetches.md)** module automates the discovery and retrieval of geospatial data. It interfaces with over forty public APIs and repositories to download topographic and bathymetric data.

**Key Capabilities:**

* **Sources**: Access data from federal agencies including NOAA (NOS, NCEI, OCM), USGS (3DEP, CoNED), USACE, and more.
* **Spatial Querying**: Download data specific to a bounding box or vector polygon.
* **Integration**: Seamlessly pipes downloaded data into processing workflows.

### Dlim (Data Lists & Processing)

The **[dlim](/docs/dlim.md)** module serves as the primary data ingestion engine. It utilizes "datalists"—text-based configuration files—to manage, process, and standardize diverse datasets into a unified stream for gridding.

**Key Capabilities:**

* **Abstraction**: Handles mixed formats (Raster, XYZ, LAS/LAZ, COPC) seamlessly, abstracting file-specific parsing.
* **Manipulation**: Supports on-the-fly vertical datum transformation, spatial clipping, and masking.
* **Hierarchy**: Supports nested datalists for organized, multi-scale data management.

### Pointz (Point Cloud Filtering)

The **[pointz](/docs/pointz.md)** module provides specialized tools for filtering and manipulating raw point cloud data (XYZ, LAS/LAZ) prior to gridding.

**Key Capabilities:**

* **Outlier Removal**: Statistical filters (`outlierz`) to remove noise spikes and gross errors.
* **Quality Control**: Reference-based checks (`rq`) to validate data against known baselines.
* **Masking**: Filter points using vector polygons (`vector_mask`) or raster masks.
* **Thinning**: Reduce data density (`block_minmax`, `density`) for efficient processing.

### Grits (Grid Filters)

The **[grits](/docs/grits.md)** module provides a standardized framework for raster manipulation. It supports chunked processing for large files and can chain multiple filters together.

**Available Filters:**

* **Smoothing**: `blur` (Gaussian), `denoise` (Median/Bilateral), `gmtfilter` (GMT wrapper).
* **Cleaning**: `outliers` (Statistical outlier removal), `zscore` (Local anomaly detection), `flats` (Remove artifacts).
* **Restoration**: `fill` (Void filling/inpainting via IDW or Spline), `morphology` (Erosion/Dilation).
* **Integration**: `blend` (Seamless blending of datasets), `weights` (Quality-based buffering), `diff` (Change detection).
* **Hydrology**: `hydro` (Sink filling).
* **Geometry**: `cut` (Mask to region), `clip` (Clip to vector), `slope` (Slope-based filtering).

### Waffles (Gridding & Interpolation)

The **[waffles](/docs/waffles.md)** module serves as the primary DEM generation engine within the framework. It provides a unified factory interface to transform scattered elevation data into seamless raster surfaces using a wide range of gridding algorithms. It manages the entire production pipeline, from data ingestion and spatial partitioning (chunking) to interpolation and post-processing.

**Key Capabilities:**

* **Unified Interface**: Users can switch between different interpolation methods simply by changing the module flag (e.g., `-M gmt-surface` for spline, `-M IDW` for inverse distance weighting).
* **Algorithms**:
  * **Internal**: `IDW` (Inverse Distance Weighting), `stacks` (Weighted Mean/Supercedure), `linear`/`cubic`/`nearest` (Scipy-based), `natural_neighbor` (Sibson), `kriging` (Geostatistical), and `ml_interp` (Machine Learning Regression).
  * **Wrappers**: `surface` / `triangulate` (GMT), `mbgrid` (MB-System), and `gdal` (GDAL Grid).
  * **Specialized**: `coastline` (Land/Water masking), `vdatum` (Vertical datum grids), `inpaint` (Void filling), and `uncertainty` (Interpolation error estimation).


* **Scalability**: Supports "chunking" (`-K`) to process massive datasets by splitting the region into tiles, processing them in parallel, and merging the results.
* **Smart Stacking**: Pre-processes overlapping datasets into a weighted intermediate raster ("stack") to resolve conflicts before interpolation.
* **Post-Processing**: Automatically handles buffering (`-X`) to prevent edge artifacts, clipping (`-C`) to vector polygons, and generating auxiliary products like uncertainty grids (`-u`) and data masks (`-m`).

# Usage Examples

### General Wokflows

* [Generate a CRM of Point Mugu in California](/docs/example_crm_malibu.md)
* [Generate a set of Tiled CRMs of Northern California](/docs/example_crm_norcal.md)
* [Generate a set of Tiled CUDEMs of North-West Washington State](/docs/dem_examples/northern_wa/northern_wa.md)

### Specific Tasks

* **Fetching Data**: `fetches -R -90/-89/29/30 CUDEM`
* **Processing Data**: `dlim input.datalist -R -90/-89/29/30 -E .3s`
* **Generating DEM**: `waffles -M surface:tension=.35 -R -90/-89/29/30 -E .1s -O my_dem.tif input.xyz`
* **Filtering DEM**: `grits my_dem.tif -M outliers:k=2.5 -M fill:method=spline -M clip:src_ply=coast.shp -O clean_dem.tif`
* **Estimating Uncertainty**: [Docs](/docs/uncertainty.md)

# Additional Information

The CUDEM code repository is frequently updated, and code syntax is
subject to change. Please see the code help function (`--help`) for the latest code
syntax and examples. See Eakins and Grothe (2014) for more information
on the challenges of building integrated DEMs and Eakins et al. (2015)
for the initial specifications of the comprehensive DEM development
framework. See Hare et al. (2011), Amante and Eakins (2016), and Amante
(2018) for additional information on the DEM uncertainty.

For additional questions, please contact:

[Matthew.Love@colorado.edu](mailto:matthew.love@colorado.edu)

[Christopher.Amante@colorado.edu](mailto:christopher.amante@colorado.edu)

## License

---

```
     MIT License

```

---

```
     Copyright (c) 2010 - 2025 Regents of the University of Colorado
     
     Permission is hereby granted, free of charge, to any person
     obtaining a copy
     of this software and associated documentation files (the
     \"Software\"), to deal
     in the Software without restriction, including without
     limitation the rights
     to use, copy, modify, merge, publish, distribute, sublicense,
     and/or sell
     copies of the Software, and to permit persons to whom the
     Software is
     furnished to do so, subject to the following conditions:
  
     The above copyright notice and this permission notice shall be
     included in all
     copies or substantial portions of the Software.
     
     THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY
     KIND, EXPRESS OR
     IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
     MERCHANTABILITY,
     FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO
     EVENT SHALL THE
     AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
     OTHER
     LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
     ARISING FROM,
     OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
     DEALINGS IN THE
     SOFTWARE.

```

---

## Code Citation

Love, M., Amante, C., Carignan, K., MacFerrin, M., & Lim, E. (2025).
CUDEM (Version 2.0.0) 

.
[[https://github.com/ciresdem/cudem]{.ul}](https://github.com/ciresdem/cudem)

# References

Amante CJ, Love M, Carignan K, Sutherland MG, MacFerrin M, Lim E. Continuously Updated Digital Elevation Models (CUDEMs) to Support Coastal Inundation Modeling. Remote Sensing. 2023; 15(6):1702. https://doi.org/10.3390/rs15061702

Amante, C. J. (2018). Estimating coastal digital elevation model
uncertainty. Journal of Coastal Research, 34(6), 1382-1397. https://doi.org/10.2112/JCOASTRES-D-17-00211.1

Amante, C. J., & Eakins, B. W. (2016). Accuracy of interpolated
bathymetry in digital elevation models. Journal of Coastal Research, (76
(10076)), 123-133. https://doi.org/10.2112/SI76-011

Eakins, B. W., & Grothe, P. R. (2014). Challenges in building coastal
digital elevation models. Journal of Coastal Research, 30(5), 942-953.

Eakins, B. W., Danielson, J. J., Sutherland, M. G., & Mclean, S. J.
(2015). A framework for a seamless depiction of merged bathymetry and
topography along US coasts. In Proc. US Hydro. Conf (pp. 16-19).

Hare, R., Eakins, B., & Amante, C. J. (2011). Modelling bathymetric
uncertainty. The International Hydrographic Review.

### [Project Board](https://github.com/ciresdem/cudem/discussions).
