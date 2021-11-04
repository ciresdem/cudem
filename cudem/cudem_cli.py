### cudem_cli.py - DataLists IMproved
##
## Copyright (c) 2010 - 2021 CIRES Coastal DEM Team
##
## cudem_cli.py is part of CUDEM
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
### Commentary:
##
## Run cudem programs and scripts from the command-line
##
### Code:

import sys
import os
import cudem
from cudem import utils
from cudem import waffles
from cudem import dlim
from cudem import regions
from cudem import metadata
from cudem import uncertainties
import cudem.fetches.fetches as fetches

class CUDEMFactory:

    mods = {
        'dlim': {
            'description': 'process data'
        },
        'waffles': {
            'description': 'generate Digital Elevation Models'
        },
        'regions': {
            'description': 'process REGIONS'
        },
        'fetches': {
            'description': 'fetch elevation data'
        },
        'spatial_metadata': {
            'description': 'generate spatial metadata'
        },
        'uncertainties': {
            'description': 'calculate DEM uncertainties'
        },
        'bag2tif2chunks2xyz.sh': {
            'description': 'convert a BAG to chunked XYZ'
        },
        'clip_xyz.sh': {
            'description': 'clip an xyz file based on a vector'
        },
        'coastline2xyz.sh': {
            'description': 'convert a coastline vector to XYZ'
        },
        'colortable.py': {
            'description': 'generate a colortable'
        },
        'create_datalist.sh': {
            'description': 'create a datalist from data in the current directory'
        },
        'create_outliers.sh': {
            'description': 'identify outliers in a DEM'
        },
        'create_povray_template.sh': {
            'description': 'generate a POVray template from a DEM'
        },
        'create_coastline.py': {
            'description': 'generate a coastline'
        },
        'ddms.py': {
            'description': 'convert between dd and dms'
        },
        'error_distance_plots.py': {
            'description': 'generate an error/distance plots'
        },
        'gdal_null.py': {
            'description': 'generate a null grid'
        },
        'gdal_outliers.py': {
            'description': 'filter vertical outliers from a grid'
        },
        'gdal_nan2null.py': {
            'description': 'convert NaN values from a grid'
        },
        'gdal_findreplace.py': {
            'description': 'find/replace values in a grid'
        },
        'gdal_query.py': {
            'description': 'query values from a grid'
        },
        'gdal_chunk.py': {
            'description': 'parse a grid into chunks'
        },
        'gdal_crop.py': {
            'description': 'crop a grid by its nodata value'
        },
        'gdal_cut.py': {
            'description': 'cut a grid to a given region'
        },
        'gdal_clip.py': {
            'description': 'clip a grid to a vector'
        },
        'gdal_split.py': {
            'description': 'split a grid by z value'
        },
        'gdal_percentile.py': {
            'description': 'get a percentile from a grid'
        },
        'gdal_histogram.py': {
            'description': 'generate an historgram from a grid'
        },
        'gdal_hillshade.py': {
            'description': 'generate a hillshade image from a grid'
        },
        'gdal_minmax.py': {
            'description': 'get min/max values from a grid'
        },
        'grd2mesh.py': {
            'description': 'generate an unstructured grid'
        },
        'has_nulls.py': {
            'description': 'check if a grid has nodata values'
        },
        'hillshade.sh': {
            'description': 'generate a hillshade from a DEM'
        },
        'nsidc_download.py': {
            'description': 'downlaod nsidc data'
        },
        'ogr_edit_field.py': {
            'description': 'edit OGR field values'
        },
        'outliers_shp.sh': {
            'description': 'identify outliers in a DEM'
        },
        'percentiles_minmax.py': {
            'description': 'get percentiles from a grid'
        },
        'rename_shp.py': {
            'description': 'rename a shapefile'
        },
        'smooth_dem_bathy.py': {
            'description': 'smooth a DEM < 0 with a Gaussian filter'
        },
        'spatial-meta.sh': {
            'description': 'generate spatial metadata using BOUNDS'
        },
        'tif2chunks2xyz.sh': {
            'description': 'chunk a DEM and output as chunked XYZ'
        },
        'usace_interp.sh': {
            'description': 'interpolate usace cross surveys'
        },
        'vdatum_cmd.py': {
            'description': 'run NOAAs vdatum from command-line'
        },
        'x360.py': {
            'description': 'flip a DEM'
        },
        'xyz_clip.py': {
            'description': 'clip an xyz file based on a raster mask'
        },
        'xyztindex.py': {
            'description': 'generate a tile index of xyz files.'
        }
    }
    
    def __init__(self, mod=None, mod_args=None):
        self.mod = mod
        self.mod_args = mod_args

    def acquire(self):
        if self.mod == 'dlim':
            return(dlim.datalists_cli(self.mod_args))

        if self.mod == 'waffles':
            return(waffles.waffles_cli(self.mod_args))

        if self.mod == 'regions':
            return(regions.regions_cli(self.mod_args))
        
        if self.mod == 'spatial_metadata':
            return(metadata.spat_meta_cli(self.mod_args))

        if self.mod == 'uncertainties':
            return(uncertainties.uncertainties_cli(self.mod_args))
        
        if self.mod == 'fetches':
            return(fetches.fetches_cli(self.mod_args))

        if self.mod in self.mods.keys():
            return(utils.run_cmd(' '.join(self.mod_args), verbose=True))

_cudem_module_short_desc = lambda: ', '.join(
    ['{}'.format(key) for key in CUDEMFactory().mods])

_cudem_module_long_desc = lambda x: 'cudem commands:\n ' + ' '.join(
    ['\033[1m{:26}\033[0m{}\n'.format(key, x[key]['description']) for key in x])

## ==============================================
## Command-line Interface (CLI)
## $ cudem
##
## cudem cli
## ==============================================

cudem_cli_usage = """{cmd} ({c_version}): CUDEM; run CUDEM programs and scripts

  --help\t\tPrint the usage text
  --version\t\tPrint the version information

{c_formats}

CIRES DEM home page: <http://ciresgroups.colorado.edu/coastalDEM>
""".format(cmd=os.path.basename(sys.argv[0]), 
           c_version=cudem.__version__,
           c_formats=_cudem_module_long_desc(CUDEMFactory.mods))

def cudem_cli(argv=sys.argv):
    """run cudem programs and scripts

    See `cudem_cli_usage` for full cli options.
    """
    
    mod = None
    
    ## ==============================================
    ## Parse command-line arguments
    ## ==============================================

    i = 1
    while i < len(argv):
        arg = argv[i]
        if arg == '--help' or arg == '-h':
            sys.stderr.write(cudem_cli_usage)
            sys.exit(0)
        elif arg == '--version' or arg == '-v':
            sys.stdout.write('{}\n'.format(cudem.__version__))
            sys.exit(0)
        elif arg[0] == '-':
            sys.stderr.write(cudem_cli_usage)
            utils.echo_error_msg('{} is not a valid cudem cli argument'.format(arg))
            sys.exit(0)
        else:
            mod = arg
            break
            
        i += 1

    if mod is not None:
        CUDEMFactory(mod=mod, mod_args=argv[i:]).acquire()
    else:
        sys.stderr.write(cudem_cli_usage)

### End
