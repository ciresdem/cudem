# Overview

"dlim" is the elevation data processing tool using various dataset modules (Table 1). dlim's native dataset format is a "datalist". A datalist is similar to an MBSystem datalist; it is a space-delineated file containing the following columns:
data-path data-format data-weight data-uncertainty data-name data-source data-date data-resolution data-type data-horz data-vert data-url
Minimally, data-path (column 1) is all that is needed.

an associated inf and geojson file will be gerenated for each datalist
only an associated inf file will be genereated for individual datasets

Parse various dataset types by region/increments and yield data as xyz or array
recursive data-structures which point to datasets (datalist, zip, fetches, etc) are negative format numbers, e.g. -1 for datalist

supported datasets include: xyz, gdal, ogr, las/laz (laspy), mbs (MBSystem), fetches (see cudem.fetches)

Initialize a datalist/dataset using init_data(list-of-datalist-entries) where the list of datalist entries can
be any of the supported dataset formats. init_data will combine all the input datafiles into an internal scratch
datalist and process that.

If region, x_inc, and y_inc are set, all processing will go through Dataset._stacks() where the data will be combined
either using the 'supercede', 'weighted-mean', 'min' or 'max' method. Dataset._stacks will output a multi-banded gdal file with the following
bands: 1: z, 2: count, 3: weight, 4: uncerainty, 5: source-uncertainty

If want_mask is set, _stacks() will also generate a multi-band gdal raster file where each
mask band contains the mask (0/1) of a specific dataset/datalist, depending on the input data. For a datalist, each band will
contain a mask for each of the top-level datasets. If want_sm is also set, the multi-band mask
will be processed to an OGR supported vector, with a feature for each band and the metadata items (cols > 4 in the datalist entry) as fields.

Transform data between horizontal/vertical projections/datums by setting src_srs and dst_srs as 'EPSG:<>'
if src_srs is not set, but dst_srs is, dlim will attempt to obtain the source srs from the data file itself
or its respective inf file; otherwise, it will be assumed the source data file is in the same srs as dst_srs

A dataset module should have at least an `__init__` and a `yield_ds` method. `yield_ds` should yield a numpy rec array with at least
'x', 'y' and 'z' fields, and optionally 'w' and 'u' fields. Other methods that can be re-written include `parse` which yields
dlim dataset module(s).

### Code Syntax:
```
dlim (2.3.6): DataLists IMproved; Process and generate datalists

usage: dlim [ -acdghijnquwEJPRT [ args ] ] DATALIST,FORMAT,WEIGHT,UNCERTAINTY ...

Options:
  -R, --region		Restrict processing to the desired REGION 
			Where a REGION is xmin/xmax/ymin/ymax[/zmin/zmax[/wmin/wmax/umin/umax]]
			Use '-' to indicate no bounding range; e.g. -R -/-/-/-/-10/10/1/-/-/-
			OR an OGR-compatible vector file with regional polygons. 
			Where the REGION is /path/to/vector[:zmin/zmax[/wmin/wmax/umin/umax]].
			If a vector file is supplied, will use each region found therein.
  -E, --increment	Block data to INCREMENT in native units.
			Where INCREMENT is x-inc[/y-inc]
  -X, --extend		Number of cells with which to EXTEND the output DEM REGION and a 
			percentage to extend the processing REGION.
			Where EXTEND is dem-extend(cell-count)[:processing-extend(percentage)]
			e.g. -X6:10 to extend the DEM REGION by 6 cells and the processing region by 10 
			percent of the input REGION.
  -J, --s_srs		Set the SOURCE projection.
  -P, --t_srs		Set the TARGET projection. (REGION should be in target projection) 
  -D, --cache-dir	CACHE Directory for storing temp and output data.
  -Z, --z-precision	Set the target precision of dumped z values. (default is 4)
  -A, --stack-mode	Set the STACK MODE to 'mean', 'min', 'max' or 'supercede' (with -E and -R)
  -T, --stack_filter	FILTER the data stack using one or multiple filters. 
			Where FILTER is fltr_name[:opts] (see `grits --modules` for more information)
			The -T switch may be set multiple times to perform multiple filters.
			Available FILTERS: blur, grdfilter, outliers, flats
  -F, --point_filter	FILTER the POINT data using one or multiple filters. 
			Where FILTER is fltr_name[:opts] 
			The -F switch may be set multiple times to perform multiple filters.
			Available FILTERS: bin_z

  --mask		MASK the datalist to the given REGION/INCREMENTs
  --spatial-metadata	Generate SPATIAL METADATA of the datalist to the given REGION/INCREMENTs
  --archive		ARCHIVE the datalist to the given REGION[/INCREMENTs]
  --glob		GLOB the datasets in the current directory to stdout
  --info		Generate and return an INFO dictionary of the dataset
  --weights		Output WEIGHT values along with xyz
  --uncertainties	Output UNCERTAINTY values along with xyz
  --stack-node		Output stacked x/y data rather than pixel
  --quiet		Lower the verbosity to a quiet

  --modules		Display the datatype descriptions and usage
  --help		Print the usage text
  --version		Print the version information

Supported datalist formats (see dlim --modules <dataset-key> for more info): 
  datalist (-1),  zip (-2),  scratch (-3),  yxz (167),  xyz (168),  gdal (200),  bag (201),
  swot_pixc (202),  swot_hr_raster (203),  las (300),  mbs (301),  ogr (302),  icesat_atl (303),
  https (-100),  gmrt (-101),  gebco (-102),  copernicus (-103),  fabdem (-104),  nasadem (-105),
  mar_grav (-106),  srtm_plus (-107),  charts (-200),  multibeam (-201),  hydronos (-202),  ehydro (-203),
  bluetopo (-204),  ngs (-205),  tides (-206),  digital_coast (-207),  ncei_thredds (-208),  tnm (-209),
  CUDEM (-210),  CoNED (-211),  SLR (-212),  waterservies (-213),  icesat (-214),  ned (-215),  swot (-216),
  csb (-217),  emodnet (-300),  chs (-301),  hrdem (-302),  arcticdem (-303),  vdatum (-500),  hydrolakes (-600)
```