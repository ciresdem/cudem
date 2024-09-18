# fetches

Access remote elevation datasets

## Synopsis

```
fetches [ -hlqzAHR [ args ] ] MODULE ...
```

## Description

Fetches is a data download tool for obtaining publicly available elevation data from a variety of sources and can optionally list or download the fetched data for use in DEM generation. Common elevation datasets are available for a variety of data types, e.g., topographic-bathymetry lidar, multibeam swath sonar bathymetry, hydrographic soundings, compiled grids, etc., from a variety of sources, e.g., NOAA Office for Coastal Management (OCM) Digital Coast, NOAA NCEI NOS Hydro Surveys, NOAA NCEI Multibeam, USGS The National Map, and U.S. Army Corps of Engineers (USACE) Navigation Condition Surveys. Other data sources include digitized bathymetric charts or topographic maps, shorelines, satellite-derived elevations, and precisely surveyed geodetic monuments (Table 1).

## Options

`-R, --region`
Restrict processing to the desired REGION 
Where a REGION is xmin/xmax/ymin/ymax[/zmin/zmax[/wmin/wmax]]
Use '-' to indicate no bounding range; e.g. -R -/-/-/-/-10/10/1/-
OR an OGR-compatible vector file with regional polygons. 
Where the REGION is /path/to/vector[:zmin/zmax[/wmin/wmax]].
If a vector file is supplied, will use each region found therein.

`-H, --threads`
Set the number of threads (1)

`-A, --attempts`
Set the number of fetching attempts (5)

`-l, --list`
Return a list of fetch URLs in the given region.

`-z, --no_check_size`
Don't check the size of remote data if local data exists.

`-q, --quiet`
Lower the verbosity to a quiet

`--modules`
Display the module descriptions and usage

`--help`
Print the usage text

`--version`
Print the version information

**Table 1.** Data source modules available in the CUDEM software tool
"fetches*.*"

| ***Name***            | ***Description***     | ***URL***             |
|---|---|---|
| arcticdem | Arctic DEM | https://www.pgc.umn.edu/data/arcticdem/ |
| bluetopo | A curated collection of high resolution seafloor models from NOAA. | https://www.nauticalcharts.noaa.gov/data/bluetopo.html |
| buoys | Buoy information from NOAA | https://www.ndbc.noaa.gov |
| charts | NOS Nautical Charts, including electronic Nautical Charts and Raster Nautical Charts | https://www.charts.noaa.gov/ |
| chs | Canadian Hydrographic Surveys | https://open.canada.ca |
| copernicus | Copernicus elevation data | https://ec.europa.eu/eurostat/web/gisco/geodata/reference-data/elevation/copernicus-dem/elevation |
| cpt_city | CPT files from CPT City | http://soliton.vm.bytemark.co.uk/pub/ |
| csb | NOAA Crowd Sourced Bathymetry | https://noaa-dcdb-bathymetry-pds.s3.amazonaws.com |
| digital_coast | Lidar and Raster data from NOAA’s Digital Coast | https://coast.noaa.gov |
| earthdata | NASA Earthdata | https://cmr.earthdata.nasa.gov |
| ehydro | USACE hydrographic surveys | https://navigation.usace.army.mil/Survey/Hydro |
| emodnet | EmodNET European Bathymetric/Topographic DEM | https://portal.emodnet-bathymetry.eu/ |
| etopo | The ETOPO Global Relief Model integrates topography, bathymetry, and shoreline data from regional and global datasets to enable comprehensive, high resolution renderings of geophysical characteristics of the earth’s surface. | https://www.ncei.noaa.gov/products/etopo-global-relief-model |
| fabdem | FABDEM  (Forest And Buildings removed Copernicus DEM) | https://data.bris.ac.uk/data/dataset/s5hqmjcdj8yo2ibzi9b4ew3sn |
| gebco | A global continuous terrain model for ocean and land with a spatial resolution of 15 arc seconds. | https://www.gebco.net/data_and_products/gridded_bathymetry_data/ |
| gmrt | The Global MultiResolution Topography synthesis | https://www.gmrt.org |
| hrdem | High-Resolution DEMs from Canada | https://open.canada.ca |
| hydrolakes | HydroLakes vector and derived elevations | https://www.hydrosheds.org/products/hydrolakes |
| icesat2 | IceSat2 granules from NASA (requires NASA credentials) | https://cmr.earthdata.nasa.gov |
| mar_grav | Marine Gravity Satellite Altimetry Topography from Scripps. | https://topex.ucsd.edu/WWW_html/mar_grav.html |
| mgds | Marine Geoscience Data System | https://www.marine-geo.org |
| multibeam | NOAA Multibeam bathymetric data | https://data.ngdc.noaa.gov/platforms/ |
| mur_sst | Sea Surface Tempuratures from NASA | https://cmr.earthdata.nasa.gov |
| nasadem | NASA Digital Elevation Model | https://www.earthdata.nasa.gov/esds/competitive-programs/measures/nasadem |
| ncei_thredds | NCEI DEM THREDDS Catalog | https://www.ngdc.noaa.gov/thredds/demCatalog.xml |
| ngs | NGS monuments | http://geodesy.noaa.gov/ |
| nos | NOS Hydrographic DataBase (NOSHDB) | https://www.ngdc.noaa.gov/mgg/bathymetry/hydro.html |
| osm | Open Street Map | https://wiki.openstreetmap.org/ |
| srtm_plus | SRTM15+: Global bathymetry and topography at 15 arc-seconds. | https://topex.ucsd.edu/WWW_html/srtm15_plus.html |
| swot | SWOT data from NASA | https://podaac.jpl.nasa.gov/SWOT?tab=datasets-information |
| tides | Tide station information from NOAA | https://tidesandcurrents.noaa.gov/ |
| tnm | USGS National Map | http://tnmaccess.nationalmap.gov/ |
| trackline | NOAA trackline bathymetry data | http://www.ngdc.noaa.gov/trackline/ |
| usiei | US Interagency Elevation Inventory | https://coast.noaa.gov/inventory/ |
| waterservices | WaterServices station information from USGS | https://waterservices.usgs.gov/ |
| wsf | World Settlement Footprint from DLR (German Aerospace Center) | https://download.geoservice.dlr.de/WSF2019/ |
| vdatum | Vertical Datum transformation grids | https://vdatum.noaa.gov https://cdn.proj.org/ |

## Python API

## Examples

### Fetch IceSat2 data with the fetches Python API

```python
from cudem import regions
from cudem import fetches
r = regions.Region().from_list([-121.25, -121.0, 37.5, 37.75])
f = fetches.IceSat2(subset=False, time_start='2020-05-04', time_end='2021-05-04', src_region=r)
f.run()
f.results
[['https://data.nsidc.earthdatacloud.nasa.gov/nsidc-cumulus-prod-protected/ATLAS/ATL03/006/2020/05/19/ATL03_20200519193744_08290702_006_01.h5', 'ATL03_20200519193744_08290702_006_01.h5', 'ATL03'], ['https://data.nsidc.earthdatacloud.nasa.gov/nsidc-cumulus-prod-protected/ATLAS/ATL03/006/2020/06/16/ATL03_20200616062933_12480706_006_01.h5', 'ATL03_20200616062933_12480706_006_01.h5', 'ATL03'], ['https://data.nsidc.earthdatacloud.nasa.gov/nsidc-cumulus-prod-protected/ATLAS/ATL03/006/2020/07/15/ATL03_20200715050536_03030806_006_01.h5', 'ATL03_20200715050536_03030806_006_01.h5', 'ATL03'], ['https://data.nsidc.earthdatacloud.nasa.gov/nsidc-cumulus-prod-protected/ATLAS/ATL03/006/2020/07/20/ATL03_20200720164127_03870802_006_01.h5', 'ATL03_20200720164127_03870802_006_01.h5', 'ATL03'], ['https://data.nsidc.earthdatacloud.nasa.gov/nsidc-cumulus-prod-protected/ATLAS/ATL03/006/2020/09/15/ATL03_20200915020921_12480806_006_02.h5', 'ATL03_20200915020921_12480806_006_02.h5', 'ATL03'], ['https://data.nsidc.earthdatacloud.nasa.gov/nsidc-cumulus-prod-protected/ATLAS/ATL03/006/2020/10/14/ATL03_20201014004522_03030906_006_01.h5', 'ATL03_20201014004522_03030906_006_01.h5', 'ATL03'], ['https://data.nsidc.earthdatacloud.nasa.gov/nsidc-cumulus-prod-protected/ATLAS/ATL03/006/2020/10/19/ATL03_20201019122113_03870902_006_01.h5', 'ATL03_20201019122113_03870902_006_01.h5', 'ATL03'], ['https://data.nsidc.earthdatacloud.nasa.gov/nsidc-cumulus-prod-protected/ATLAS/ATL03/006/2020/11/17/ATL03_20201117105718_08290902_006_01.h5', 'ATL03_20201117105718_08290902_006_01.h5', 'ATL03'], ['https://data.nsidc.earthdatacloud.nasa.gov/nsidc-cumulus-prod-protected/ATLAS/ATL03/006/2020/12/14/ATL03_20201214214913_12480906_006_01.h5', 'ATL03_20201214214913_12480906_006_01.h5', 'ATL03'], ['https://data.nsidc.earthdatacloud.nasa.gov/nsidc-cumulus-prod-protected/ATLAS/ATL03/006/2021/01/18/ATL03_20210118080103_03871002_006_01.h5', 'ATL03_20210118080103_03871002_006_01.h5', 'ATL03'], ['https://data.nsidc.earthdatacloud.nasa.gov/nsidc-cumulus-prod-protected/ATLAS/ATL03/006/2021/03/15/ATL03_20210315172905_12481006_006_01.h5', 'ATL03_20210315172905_12481006_006_01.h5', 'ATL03'], ['https://data.nsidc.earthdatacloud.nasa.gov/nsidc-cumulus-prod-protected/ATLAS/ATL03/006/2021/04/19/ATL03_20210419034059_03871102_006_02.h5', 'ATL03_20210419034059_03871102_006_02.h5', 'ATL03']]

fr = fetches.fetch_results(f)
fr.daemon=True
fr.start()
fr.join()
```

### Fetch IceSat2 data as a subset with the fetches CLI

```bash
fetches -R -121.25/-121/37.5/37.75 icesat2:subset=True:time_start=2020-05-04:time_end=2021-05-04
```