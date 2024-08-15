## Code Syntax
```
fetches (2.3.6): Fetches; Fetch and process remote elevation data

usage: fetches [ -hlqzAHR [ args ] ] MODULE ...

Options:
  -R, --region		Restrict processing to the desired REGION 
			Where a REGION is xmin/xmax/ymin/ymax[/zmin/zmax[/wmin/wmax]]
			Use '-' to indicate no bounding range; e.g. -R -/-/-/-/-10/10/1/-
			OR an OGR-compatible vector file with regional polygons. 
			Where the REGION is /path/to/vector[:zmin/zmax[/wmin/wmax]].
			If a vector file is supplied, will use each region found therein.
  -H, --threads		Set the number of threads (1)
  -A, --attempts	Set the number of fetching attempts (5)
  -l, --list		Return a list of fetch URLs in the given region.
  -z, --no_check_size	Don't check the size of remote data if local data exists.
  -q, --quiet		Lower the verbosity to a quiet

  --modules		Display the module descriptions and usage
  --help		Print the usage text
  --version		Print the version information

Supported FETCHES modules (see fetches --modules <module-name> for more info): 
  gmrt, mar_grav, srtm_plus, charts, digital_coast, SLR, CoNED, CUDEM, multibeam, gebco,
  mgds, trackline, ehydro, ngs, hydronos, ncei_thredds, etopo, tnm, ned, ned1, emodnet,
  chs, hrdem, arcticdem, bluetopo, osm, copernicus, fabdem, nasadem, tides, vdatum, buoys,
  earthdata, icesat, mur_sst, swot, usiei, wsf, hydrolakes, https, bing_bfp, waterservices,
  csb, cpt_city  
```