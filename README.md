Matthew Love[^1][^2], Christopher Amante[^1][^2]

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

The CUDEM framework provides a set of powerful open-source command-line tools as well as a Python-3-based application programming interface (API) for (1) identifying and downloading topographic and bathymetric source data from more than forty different publicly available datasets, (2) filtering data, (3) ranking and merging multiple datasets using a variety of algorithms, (4) generating DEMs from those source datasets in a variety of data formats and projections/datums, and (5) generating metadata documentation of those DEMs. Individual components are modular and may be run standalone; for example, the CUDEM “fetches” module can download airborne lidar data, multi-beam sonar data, or nautical charts covering specific geographic regions; or the “stacks” module can efficiently combine and gap-fill raster datasets of varying qualities using customizable data lists. The CUDEM tools may also be conjoined in end-to-end workflows ensuring permanent traceability and reproducibility in DEM generation. The Coastal DEM Team uses the CUDEM framework to accurately inform and map risk assessments of coastal regions (e.g. Amante et al.)

# Installation and Setup

## pip
GDAL and GDAL-Python are required for use.

Other useful external programs needed for full functionality include:
GMT, MB-System, HTDP and VDatum.

Installation of the above programs are system dependent. Only GDAL and GDAL-Python are required.

Download and install git (If you have not already): [git installation](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git)

```
pip install git+https://github.com/ciresdem/cudem.git#egg=cudem
```

## conda
- Setup a conda environment and install the dependencies:

```
conda create -n cudem -c conda-forge gdal gmt numpy scipy pandas pyproj utm requests lxml matplotlib laspy h5py boto3 tqdm mercantile git
conda activate cudem
(cudem) pip install laspy[laszip]
(cudem) pip install --no-deps git+https://github.com/ciresdem/cudem.git
```

## Extras
- If on Windows, install windows-curses with pip

```
git config --global http.sslVerify false
pip install windows-curses
```

- install [HTDP](https://geodesy.noaa.gov/TOOLS/Htdp/Htdp.shtml) and [MB-System](https://www.mbari.org/technology/mb-system/installation/)

Installation of HTDP and MB-System are system dependent, see their respective installation instructions for your system.



# Modules

The open-access code includes command-line tools and a Python
application programming interface (API) for automated data download,
processing, DEM gridding, and interpolation uncertainty grid generation
with three main software tools: "fetches", "waffles", and "dlim".
"Fetches" is the data download tool for obtaining publicly available
elevation data froma variety of sources and can optionally list, download or
process thefetched data for use in DEM generation. We download a variety of
data types, e.g., topographic-bathymetry lidar, multibeam swath sonar
bathymetry, hydrographic soundings, compiled grids, etc., from a variety
of sources, e.g., NOAA Office for Coastal Management (OCM) Digital
Coast, NOAA NCEI NOS Hydro Surveys, NOAA NCEI Multibeam, USGS The
National Map, and U.S. Army Corps of Engineers (USACE) Navigation
Condition Surveys. Other data sources include digitized bathymetric
charts or topographic maps, shorelines, satellite-derived elevations,
and precisely surveyed geodetic monuments (Table 1). We typically
download data in an area slightly larger (\~5%) than the DEM extents.
This data "buffer" ensures that interpolative gridding occurs across
rather than along the DEM boundaries to prevent edge effects, which is
especially important with sparse bathymetric data with large
interpolation distances. Data buffers also minimize artificial offsets
between adjacent DEM tiles.

**Main Console Programs and Python APIs provided with CUDEM:**

| Module | Description |
|---|---|
| [dlim](/docs/dlim.md) | process data from a variety of data types |
| [waffles](/docs/waffles.md) | generate Digital Elevation Models from scattered data using a variety of methods|
| [fetches](/docs/fetches.md) | fetch elevation data from a variety of public online sources |
| [regions](/docs/regions.md) | process REGIONS |
| [vdatums](/docs/vdatums.md) | generate vertical transformation grids |
| [perspecto](/docs/perspecto.md) | generate images of DEMs |
| [grits](/docs/grits.md) | filter DEMs |
| cudem   | run CUDEM cli programs and scripts |

# Usage Examples

- [Generate a CRM of Point Mugu in California](/docs/example_crm_malibu.md)
- [Generate a set of Tile CRMs of Northern California](/docs/example_crm_norcal.md)
- [Generate a DEM and Uncertainty grid of The Bahamas](/docs/example_uncertainty.md)
- [Filter Multibeam data](/docs/example_grits.md)

- [Canada-US elevation model collaboration to improve tsunami inundation mapping](https://www.youtube.com/watch?v=fc4ibBim_k0)
- [CSSP DEM Workshop](/docs/CSSP_DEM_Workshop_Report.pdf)


# Additional Information

The CUDEM code repository is frequently updated, and code syntax is
subject to change. Please see the code help function for the latest code
syntax and examples. See Eakins and Grothe (2014) for more information
on the challenges of building integrated DEMs and Eakins et al. (2015)
for the initial specifications of the comprehensive DEM development
framework. See Hare et al. (2011), Amante and Eakins (2016), and Amante
(2018) for additional information on the DEM uncertainty.

For additional questions, please contact:

[Matthew.Love@noaa.gov](mailto:Matthew.Love@noaa.gov)

[Christopher.Amante@noaa.gov](mailto:Christopher.Amante@noaa.gov)

## License

  -----------------------------------------------------------------------
         MIT License
  -----------------------------------------------------------------------
         
         Copyright (c) 2010 - 2024 Regents of the University of Colorado
         
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
  -----------------------------------------------------------------------

## Code Citation

Love, M., Amante, C., Carignan, K., MacFerrin, M., & Lim, E. (2023).
CUDEM (Version 1.10.5) \[Computer software\].
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