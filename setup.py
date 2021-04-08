### setup.py
##
## Copyright (c) 2020 CIRES Coastal DEM Team
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

with open('README', 'r') as fh:
    long_description = fh.read()

setuptools.setup(
    name = 'cudem',
    version = '1.0.0',
    description = 'Modules and scripts for utilizing geographic data Digital Elevation Models',
    long_description = long_description,
    long_description_content_type = 'text/markdown',
    license = 'MIT',
    author = 'CIRES Coastal DEM Team',
    author_email = 'matthew.love@colorado.edu',
    url = 'http://ciresgroups.colorado.edu/coastalDEM',
    #packages = setuptools.find_packages(),#['cudem'],  #same as name
    packages = ['cudem'],
    package_data = {'cudem': ['data/*.geojson']},
    classifiers = [
        'Programming Language :: Python :: 3',
        'License :: OSI APPROVED :: MIT License',
        'Operating System :: OS Independent',
    ],
    install_requires = [
        'GDAL',
        'numpy',
        'scipy',
        'requests',
        'lxml',
        'matplotlib',
    ], 
    entry_points = {
        'console_scripts': [
            'dlim = cudem.dlim:datalists_cli',
            'regions = cudem.regions:regions_cli',
            'waffles = cudem.waffles:waffles_cli',
            'fetches = cudem.fetches:fetches_cli',
            'spatial_metadata = cudem.metadata:spat_meta_cli',
        ],
    },
    scripts = [
        'scripts/gdal_null.py',
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
        'scripts/ddms.py',
        'scripts/has_nulls.py',
        'scripts/xyz_clip.py',
        'scripts/clip_xyz.sh',
        'scripts/percentiles_minmax.py',
        'scripts/outliers_shp.sh',
        'scripts/usace_interp.sh',
        'scripts/spatial-meta.sh',
        'scripts/hillshade.sh',
        'scripts/colortable.py',
        'scripts/coastline2xyz.sh',
        'scripts/create_datalist.sh',
        'scripts/create_outliers.sh',
        'scripts/create_povray_template.sh',
        'scripts/tif2chunks2xyz.sh',
        'scripts/bag2tif2chunks2xyz.sh',
    ],
    python_requires = '>=2.7',
    project_urls = {
        'Source': 'https://github.com/ciresdem/cudem',
    },
)
