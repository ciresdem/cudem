# dlim

process various elevation datasets.

## Synopsis

```
dlim [ -acdghijnquwAEJPRT [ args ] ] DATALIST,FORMAT,WEIGHT,UNCERTAINTY ...
```

## Description

`dlim` is the elevation data processing tool using various dataset modules (Table 1). `dlim`'s native dataset format is a "datalist". A datalist is similar to an MBSystem datalist; it is a space-delineated file containing the following columns:

```data-path data-format data-weight data-uncertainty data-name data-source data-date data-resolution data-type data-horz data-vert data-url```

Minimally, `data-path` (column 1) is all that is needed.

An associated `inf` and geojson file will be gerenated for each datalist while only an associated `inf` file will be genereated for individual datasets

Parse various dataset types by region/increments and yield data as xyz or array. Recursive data-structures, which point to 
datasets (datalist, zip, fetches, etc), are negative format numbers, e.g. -1 for datalist. Fetches modules are <= -100.

### Options

`-R, --region`

> Restrict processing to the desired REGION \
> Where a REGION is xmin/xmax/ymin/ymax[/zmin/zmax[/wmin/wmax/umin/umax]]\
> Use '-' to indicate no bounding range; e.g. -R -/-/-/-/-10/10/1/-/-/-\
> OR an OGR-compatible vector file with regional polygons. \
> Where the REGION is /path/to/vector[:zmin/zmax[/wmin/wmax/umin/umax]].\
> If a vector file is supplied, will use each region found therein.
     
`-E, --increment`

> Block data to INCREMENT in native units.\
> Where INCREMENT is x-inc[/y-inc]

`-X, --extend`

> Number of cells with which to EXTEND the output DEM REGION and a \
> percentage to extend the processing REGION.\
> Where EXTEND is dem-extend(cell-count)[:processing-extend(percentage)]\
> e.g. -X6:10 to extend the DEM REGION by 6 cells and the processing region by 10 \
> percent of the input REGION.
     
`-J, --s_srs`

> Set the SOURCE projection.

`-P, --t_srs`

> Set the TARGET projection. (REGION should be in target projection)
     
`-D, --cache-dir`

> CACHE Directory for storing temp and output data.
     
`-Z, --z-precision`

> Set the target precision of dumped z values. (default is 4)
     
`-A, --stack-mode`

> Set the STACK MODE to 'mean', 'min', 'max', 'mixed' or 'supercede' (with -E and -R)
     
`-T, --stack_filter`

> FILTER the data stack using one or multiple filters. \
> Where FILTER is fltr_name[:opts] (see `grits --modules` for more information)\
> The -T switch may be set multiple times to perform multiple filters.\
> Available FILTERS: blur, grdfilter, outliers, flats\
     
`-F, --point_filter`

> FILTER the POINT data using one or multiple filters. \
> Where FILTER is fltr_name[:opts] (See dlim --point-filters for more information)\
> The -F switch may be set multiple times to perform multiple filters.\
> Available FILTERS: outlierz, rq

`-V, --archive`

> ARCHIVE the datalist to the given REGION[/INCREMENTs]

`-m, --mask`

> MASK the datalist to the given REGION/INCREMENTs
	
`-s, --spatial-metadata`

> Generate SPATIAL METADATA of the datalist to the given REGION/INCREMENTs
	
`-g, --glob`

> GLOB the datasets in the current directory to stdout
	
`-i, --info`

> Generate and return an INFO dictionary of the dataset
	
`-w, --weights`

> Output WEIGHT values along with xyz
	
`-u --uncertainties`

> Output UNCERTAINTY values along with xyz
	
`-n, --stack-node`

> Output stacked x/y data rather than pixel
	
`-q, --quiet`

> Lower the verbosity to a quiet

`--point-filters`

> Display the POINT FILTER descriptions and usage

`--modules`

> Display the datatype descriptions and usage
	
`--help`

> Print the usage text
	
`--version`

> Print the version information


**Table 1.** Dataset modules available in the CUDEM software tool "dlim"

|  ***Name***  |  ***dlim dataset module code*** | ***Description*** |
|----------------------|----------------------------------|----------------------------------|
| datalist | -1 | An extended MB-System style datalist containting dlim-compatible datasets |
| zip | -2 | A zipfile containing dlim-compatible datasets |
| scratch | -3 | A scratch dataset, including a python list of dlim-compatible datasets |
| yxz | 167 | An ascii DSV datafile formatted as y,x,z |
| xyz | 168 | An ascii DSV datafile formatted as x,y,z |
| gdal | 200 | A gdal-compatible raster dataset |
| bag | 201 | A BAG bathymetry dataset |
| swot_pixc | 202 | An HDF5 SWOT PIXC datafile |
| swot_hr_raster | 203 | an HDF5 SWOT HR Raster datafile |
| las | 300 | An las or laz lidar datafile |
| mbs | 301 | An MB-System-compatible multibeam datafile |
| ogr | 302 | An ogr-compatible vector datafile |
| icesat2_atl | 303 | An HDF5 IceSat2 ATL03 datafile |
| https | -100 | a URL pointing to a dlim-compatible datafile |
| gmrt | -101 | The GMRT fetches module |
| gebco | -102 | The GEBCO fetches module |
| copernicus | -103 | The Copernicus fetches module |
| fabdem | -104 | The FABDEM fetches module |
| nasadem | -105 | The NASADEM fetches module |
| mar_grav | -106 | The mar_grav fetches module |
| srtm_plus | -107 | The srtm_plus fetches module |
| synbath | -108 | The synbath fetches module |
| gedtm30 | -109 | The Global DTM fetches module |
| swot | -110 | The SWOT fetches module |
| icesat2 | -111 | The IceSat2 fetches module |
| charts | -200 | The charts fetches module |
| multibeam | -201 | The multibeam fetches module | 
| hydronos | -202 | The hydronos fetches module |
| ehydro | -203 | The ehydro fetches module |
| bluetopo | -204 | The bluetopo fetches module |
| ngs | -205 | The ngs fetches module |
| tides | -206 | The tids fetches module |
| digital_coast | -207 | The digital_coast fetches module |
| ncei_thredds | -208 | The ncei_thredds fetches module |
| tnm | -209 | The TNM fetches module |
| CUDEM | -210 | The CUDEM fetches module |
| CoNED | -211 | The CoNED fetches module |
| SLR | -212 | The SLR fetches module |
| waterservies | -213 | The waterservices fetches module |
| ned | -214 | The NED fetches module |
| csb | -215 | the CSB fetches module |
| wa_dnr | -216 | The Washington State DNR fetches module |
| r2r | -217 | Rolling deck to repository fetches module |
| emodnet | -300 | The emodnet fetches module |
| chs | -301 | The chs fetches module |
| hrdem | -302 | The hrdem fetches module | 
| arcticdem | -303 | The arcticdem fetches module |
| vdatum | -500 | The vdatum fetches module | 
| hydrolakes | -600 | The hydrolakes fetches module | 

## Python API

## Examples

### Load and process data to xyz using the Python API

```python
from cudem import regions, dlim

# set the region
r = regions.Region().from_string('-R-119.25/-119/34/34.25')

# acquire and initialize the dlim dataset module for xyz file `my_hydronos.xyz`
d = dlim.DatasetFactory(mod='my_hydronos.xyz', src_region=r)._acquire_module().initialize()

# dump the xyz data from the dataset to stdout
d.dump_xyz()

# dump the xyz data to a file
with open('test_dump.xyz', 'w') as dump_file:
     d.dump_xyz(dump_file)

# make and initialize a datalist object from multiple datasets with the given region at 1 arc-second output blocks
# this fetches NED data from the USGS using the `ned` fetches module as a dataset
# set the output projection to WGS84 (epsg:4326)
d = dlim.init_data(['my_hydronos.xyz', 'ned'], region=r, xyinc=('1s','1s'), dst_srs='epsg:4326').initialize()

# dump the xyz data, include the uncertainty
# the dumped xyz data will be in wgs84 at 1 arc-second blocks
d.dump_xyz(include_u=True)
```

### Process fetches dataset `hydronos` to input region and output to an xyz file, including weights and uncertainty

```bash
dlim -R-119.25/-119/34/34.25 hydronos -w -u > my_hydronos.xyz
```

### Process all local `*.laz` lidar files, transform them from UTM zone 10N to WGS84 and block them to the input region at 1 arc-second increments

```bash
dlim -R-123.458/-123.448/38.705/38.711 -E1s -J epsg:26910 -P epsg:4326 *.laz > my_lidar_1as.xyz
```

### Generate a data-mask (raster) and spatial-metadata (vector) of a datalist:

```bash
dlim -Rmy_region.shp -E1s -Pepsg:4269 --mask --spatial-metadata my_data.datalist > /dev/null
```