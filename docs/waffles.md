### Code Syntax:
```
waffles (2.3.6): Generate DEMs and derivatives.

usage: waffles [OPTIONS] DATALIST

Options:
  -R, --region			Specifies the desired REGION;
				Where a REGION is xmin/xmax/ymin/ymax[/zmin/zmax[/wmin/wmax]]
				Use '-' to indicate no bounding range; e.g. -R -/-/-/-/-10/10/1/-
				OR an OGR-compatible vector file with regional polygons. 
				Where the REGION is /path/to/vector[:zmin/zmax[/wmin/wmax]].
				If a vector file is supplied, will use each region found therein.
  -E, --increment		Gridding INCREMENT and RESAMPLE-INCREMENT in native units.
				Where INCREMENT is x-inc[/y-inc][:sample-x-inc/sample-y-inc]
  -M, --module			Desired Waffles MODULE and options. (see available Modules below)
				Where MODULE is module[:mod_opt=mod_val[:mod_opt1=mod_val1[:...]]]
  -S, --sample_alg		ReSAMPLE algorithm to use (from gdalwarp)
				Set as 'auto' to use 'average' when down-sampling and 'bilinear' when up-sampling
				This switch controls resampling of input raster datasets as well as resampling
				the final DEM if RESAMPLE-INCREMENT is set in -E
  -X, --extend			Number of cells with which to EXTEND the output DEM REGION and a 
				percentage to extend the processing REGION.
				Where EXTEND is dem-extend(cell-count)[:processing-extend(percentage)]
				e.g. -X6:10 to extend the DEM REGION by 6 cells and the processing region by 10
				percent of the input REGION.
  -T, --filter			FILTER the output DEM using one or multiple filters. 
				Where FILTER is fltr_name[:opts] (see `grits --modules` for more information)
				The -T switch may be set multiple times to perform multiple filters.
				Append `:stacks=True` to the filter to perform the filter on the data stack instead 
				of the final DEM.
				Available FILTERS: blur, grdfilter, outliers, flats
  -L, --limits			LIMIT the output elevation or interpolation values, append 
				'u<value>' to set the upper elevation limit, 
				'l<value>' to set the lower elevation limit,
				'n<value>' to set the count per cell limit.
				'p<value>' to set an interpolation limit by proximity, or 
				's<value>' to set an interpolation limit by size, or
				'c<value>' to set an interpolation limit by nodata-size percentile.
				e.g. -Lu0 to set all values above 0 to zero, or 
				-Ls100 to limit interpolation to nodata zones smaller than 100 pixels, or
				-Ln2 to limit stacked cells to those with at least 2 contributing data records.
  -C, --clip			CLIP the output to the clip polygon -C<clip_ply.shp:invert=False>
  -K, --chunk			Generate the DEM in CHUNKs.
  -F, --format			Output grid FORMAT. [GTiff]
  -O, --output-name		BASENAME for all outputs.
  -P, --t_srs			Projection of REGION and output DEM.
  -N, --nodata			The NODATA value of output DEM.
  -A, --stack-mode		Set the STACK MODE to 'mean', 'min', 'max' or 'supercede' (higher weighted data supercedes lower weighted data).
  -G, --wg-config		A waffles config JSON file. If supplied, will overwrite all other options.
				Generate a waffles_config JSON file using the --config flag.
  -H, --threads			Set the number of THREADS (1). Each input region will be run in up to THREADS threads. 
  -D, --cache-dir		CACHE Directory for storing temp data.
				Default Cache Directory is ~/.cudem_cache; cache will be cleared after a waffles session.
				to retain the data, use the --keep-cache flag.
  -CO, --creation-options	GDAL CREATION OPTIONS to use in raster creation.

  -f, --transform		Transform all data to PROJECTION value set with --t_srs/-P where applicable.
  -p, --prefix			Set BASENAME (-O) to PREFIX (append <RES>_nYYxYY_wXXxXX_<YEAR>v<VERSION> info to output BASENAME).
				note: Set Resolution, Year and Version by setting this to 'res=X:year=XXXX:version=X', 
				leave blank for default of <INCREMENT>, <CURRENT_YEAR> and <1>, respectively.
  -r, --grid-node		Use grid-node registration, default is pixel-node.
  -w, --want-weight		Use weights provided in the datalist to weight overlapping data.
  -u, --want-uncertainty	Generate/Use uncertainty either calculated or provided in the datalist.
  -m, --want-mask		Mask the processed datalist.
  -a, --archive			ARCHIVE the datalist to the given region.
  -k, --keep-cache		KEEP the cache data intact after run
  -x, --keep-auxiliary		KEEP the auxiliary rastesr intact after run (mask, uncertainty, weights, count).
  -s, --spatial-metadata	Generate SPATIAL-METADATA.
  -c, --continue		Don't clobber existing files.
  -q, --quiet			Lower verbosity to a quiet.

  --help			Print the usage text
  --config			Save the waffles config JSON and major datalist
  --modules			Display the module descriptions and usage
  --version			Print the version information

Datalists and data formats:
  A datalist is a file that contains a number of datalist entries, 
  while an entry is a space-delineated line:
  `path [format weight uncertainty [name source date type resolution hdatum vdatum url]]`

Supported datalist formats: 
  datalist (-1),  zip (-2),  scratch (-3),  yxz (167),  xyz (168),  gdal (200),  bag (201),
  swot_pixc (202),  swot_hr_raster (203),  las (300),  mbs (301),  ogr (302),  icesat_atl (303),
  https (-100),  gmrt (-101),  gebco (-102),  copernicus (-103),  fabdem (-104),  nasadem (-105),
  mar_grav (-106),  srtm_plus (-107),  charts (-200),  multibeam (-201),  hydronos (-202),
  ehydro (-203),  bluetopo (-204),  ngs (-205),  tides (-206),  digital_coast (-207),  ncei_thredds (-208),
  tnm (-209),  CUDEM (-210),  CoNED (-211),  SLR (-212),  waterservies (-213),  icesat (-214),  ned (-215),
  swot (-216),  csb (-217),  emodnet (-300),  chs (-301),  hrdem (-302),  arcticdem (-303),  vdatum (-500),
  hydrolakes (-600)

Modules (see waffles --modules <module-name> for more info):
  stacks, IDW, linear, cubic, nearest, gmt-surface, gmt-triangulate, gmt-nearneighbor, mbgrid, gdal-linear,
  gdal-nearest, gdal-average, gdal-invdst, vdatum, coastline, lakes, cudem, uncertainty, scratch, flatten
```