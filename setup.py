### setup.py
##
## Copyright (c) 2020 - 2025 Regents of the University of Colorado
##
## Permission is hereby granted, free of charge, to any person obtaining a copy 
## of this software and associated documentation files (the "Software"), to deal 
## in the Software without restriction, including without limitation the rights 
## to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies 
## of the Software, and to permit persons to whom the Software is furnished to do so, 
## subject to the following conditions:
##
## The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
##
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
## INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR 
## PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE 
## FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
## ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
##
### Code:
import setuptools
import site
import sys

site.ENABLE_USER_SITE = "--user" in sys.argv[1:]

with open('README.md', 'r', encoding='utf-8') as fh:
    long_description = fh.read()

setuptools.setup(
    name = 'cudem',
    version = '2.4.2',
    description = 'Modules and scripts for utilizing geographic data Digital Elevation Models',
    long_description = long_description,
    long_description_content_type = 'text/markdown',
    license = 'MIT',
    author = 'CIRES Coastal DEM Team',
    author_email = 'matthew.love@colorado.edu',
    url = 'http://ciresgroups.colorado.edu/coastalDEM',
    packages = ['cudem'],
    package_data = {'cudem': ['data/*.geojson', 'data/*_errs.dat.gz']},
    classifiers = [
        'Programming Language :: Python :: 3',
        'License :: OSI APPROVED :: MIT License',
        'Operating System :: OS Independent',
    ],
    install_requires = [
        'numpy', # all
        'scipy', # waffles/convex hulls
        'requests[security]', # fetches/nsidc
        'lxml', # fetches xml parsing
        'matplotlib', # uncertainty plots
        'laspy[lazrs,laszip]', # dlim/lasfile support
        'pyproj', # all
        'h5py', # nsidc
        'boto3', # for amazon/bluetopo
        'tqdm', # progress bar
        'mercantile', # for MS building footprints
        'utm', # for cshelph
        'pandas', # for cshelph
        'netCDF4', # for netcdf outputs
        # 'windows-curses', # for windows
        # 'shapely',
        #'xarray', # for cshelph
        #'geopandas', # for cshelph
        #'cshelph', # icesat bathymetry extraction        
        #'dask', # for cshelph sst
        #'earthaccess', # for earthdata access sst
        # 'h5netcdf', # for cshelph sst
        #'zarr',
        #'zarr-eosdis-store',
        
        # 'pydap',
    ], 
    entry_points = {
        'console_scripts': [
            'regions = cudem.regions:regions_cli',
            'dlim = cudem.dlim:datalists_cli',
            'waffles = cudem.waffles:waffles_cli',
            'fetches = cudem.fetches:fetches_cli',
            'perspecto = cudem.perspecto:perspecto_cli',
            'vdatums = cudem.vdatums:vdatums_cli',
            'grits = cudem.grits:grits_cli',
            'vrbag = cudem.vrbag:vrbag_cli',
            'cudem = cudem.cudem_cli:cudem_cli',
            #'vdatumfun = cudem.vdatums2:vdatums_cli',
        ],
    },
    scripts = [
        'scripts/gdal_null.py',
        'scripts/gdal_zeros.py',
        'scripts/gdal_outliers.py',
        'scripts/gdal_nan2null.py',
        'scripts/gdal_findreplace.py',
        'scripts/gdal_query.py',
        'scripts/gdal_chunk.py',
        'scripts/gdal_crop.py',
        'scripts/gdal_cut.py',
        'scripts/gdal_clip.py',
        'scripts/gdal_split.py',
        'scripts/gdal_percentile.py',
        'scripts/gdal_histogram.py',
        'scripts/gdal_hillshade.py',
        'scripts/gdal_minmax.py',
        'scripts/gdal_mask.py',
        'scripts/gdal_polygonize_mask.py',
        'scripts/cudem_polygonize_mask.py',
        'scripts/gdal_hydro_flatten.py',
        'scripts/gdal_remove_flats.py',
        'scripts/gdal_findrivers.py',
        'scripts/ddms.py',
        'scripts/has_nulls.py',
        'scripts/all_ndv.py',
        'scripts/xyz_clip.py',
        'scripts/clip_xyz.sh',
        'scripts/percentiles_minmax.py',
        'scripts/outliers_shp.sh',
        'scripts/usace_interp.sh',
        'scripts/spatial-meta.sh',
        'scripts/hillshade.sh',
        'scripts/colortable.py',
        'scripts/coastline2xyz.sh',
        'scripts/coastline_mask.py',
        'scripts/create_datalist.sh',
        'scripts/create_outliers.sh',
        'scripts/create_povray_template.sh',
        'scripts/create_povray_perspective.sh',
        'scripts/lll.py',
        'scripts/tif2chunks2xyz.sh',
        'scripts/bag2tif2chunks2xyz.sh',
        'scripts/ogr_edit_field.py',
        'scripts/create_coastline.py',
        'scripts/error_distance_plots.py',
        'scripts/grd2mesh.py',
        'scripts/nsidc_download.py',
        'scripts/rename_shp.py',
        'scripts/smooth_dem_bathy.py',
        'scripts/vdatum_cmd.py',
        'scripts/x360.py',
        'scripts/xyztindex.py',
        'scripts/vrbag2tif2chunks2xyz.sh',
        'scripts/VR_to_image.py',
        'scripts/explode_bags.sh',
        'scripts/dem_smooth.py',
        'scripts/dem_combine_split.py',
        'scripts/dem_update.sh',
        'scripts/xyz_vs_dem.sh',
        'scripts/dem_set_metadata.py',
        'scripts/dem_patch.sh',
        'scripts/xyz2shp.py',
        'scripts/fetch_osm_coastline.py',
        'scripts/icesat2dump.py',
    ],
    python_requires = '>=3.0',
    project_urls = {
        'Source': 'https://github.com/ciresdem/cudem',
    },
)

### End
