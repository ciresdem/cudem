# Overview

"Waffles" is the data gridding tool for DEM generation using various
gridding modules (Table 1). We leverage existing open-source software
packages from a variety of sources including GMT, GDAL, and MB-System.
From these sources, there are also multiple gridding algorithms, e.g.,
spline, inverse distance weighting (IDW), triangulate, average,
near-neighbor, etc. Previous research at NOAA NCEI indicates spline
interpolation is the most accurate gridding algorithm for generating
DEMs (Amante and Eakins, 2016). However, the measurement density and
terrain characteristics (e.g. terrain slope and curvature) may influence
the accuracy of the various gridding algorithms and multiple algorithms
should be evaluated. We then generate DEMs by a combination of data
masking, weighting, and interpolation.

**Table 1.** Gridding modules available in the CUDEM software tool
"waffles."

|  ***Name***          |                ***Description*** |
|----------------------|----------------------------------|
|  gdal-average        |                Generate an average DEM using GDAL\'s *gdal_grid* command. |
|  coastline           |                Generate a coastline (land/sea mask) using a variety of data sources. |
|  cudem               |                CUDEM integrated DEM generation. Generate an topographic-bathymetric integrated DEM using a variety of data sources. |
|  IDW                 |                Generate a DEM using an Inverse Distance Weighted algorithm. If weights are used, will generate a UIDW DEM, using weight values as inverse uncertainty, as described [here](https://ir.library.oregonstate.edu/concern/graduate_projects/79407x932), and [here](https://stackoverflow.com/questions/3104781/inverse-distance-weighted-idw-interpolation-with-python) |
|  gdal-invdst         |                Generate an inverse distance DEM using GDAL\'s *gdal_grid* command. |
|  gdal-linear         |                Generate a linear DEM using GDAL\'s *gdal_grid* command. |
|  mbgrid              |                Generate a DEM using MB-System\'s *mbgrid* command. |
|  gdal-nearest        |                Generate a nearest DEM using GDAL\'s *gdal_grid* command. |
|  gmt-nearneighbor    |                Generate a DEM using GMT\'s *nearneighbor* command |
|  num                 |                Generate an uninterpolated DEM using various gridding modes, including options from GMT's *xyz2grd* command. |
|  linear               |               Generate a DEM using Scipy's linear gridding algorithm) |
|  cubic               |                Generate a DEM using Scipy's cubic gridding algorithm |
|  nearest              |               Generate a DEM using Scipy's nearest gridding algorithm |
|  stacks              |                Generate a DEM using a raster stacking method. By default, calculate the \[weighted\]-mean where overlapping cells occur. Set supercede to True to overwrite overlapping cells with higher weighted data. | 
|  gmt-surface         |               Generate a DEM using GMT\'s *surface* command | 
|  gmt-triangulate     |               Generate a DEM using GMT\'s *triangulate* command | 
|  vdatum              |                Generate a vertical datum conversion grid. |
|  flatten              |               Generate a DEM by flattening nodata areas. |
|  lakes              |               Estimate/interpolate Lake bathymetry. |
|  uncertainty         |               Estimate interpolation uncertainty. |

# Uncertainty Grid Generation

We generate uncertainty grids that represent an estimate of potential
interpolation errors based on a split-sample approach and distance to
the nearest measurement. Using a split-sample approach, a percentage of
the data is omitted, an interpolation method is applied, and the
differences between the interpolated elevations and the original omitted
elevations are calculated. We omit a percentage of the measurements,
apply an interpolation method, and calculate the differences between the
interpolated values and the omitted elevations. We repeat this process
and aggregate the differences between the original measurements and the
interpolated elevations and then we derive a best-fit equation of
interpolation uncertainty as a function of distance to the nearest
measurement.

# Spatial Metadata

We generate spatial metadata for our DEMs that document the data sources
that went into the DEM generation process.

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