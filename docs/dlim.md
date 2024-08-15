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