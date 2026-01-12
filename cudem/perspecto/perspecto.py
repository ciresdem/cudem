### perspecto.py
##
## Copyright (c) 2023 - 2026 Regents of the University of Colorado
##
## perspecto.py is part of cudem
##
## Permission is hereby granted, free of charge, to any person obtaining a copy 
## of this software and associated documentation files (the "Software"), to deal 
## in the Software without restriction, including without limitation the rights 
## to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies 
## of the Software, and to permit persons to whom the Software is furnished to do so, 
## subject to the following conditions:
##
## The above copyright notice and this permission notice shall be included in all
## copies or substantial portions of the Software.
##
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
## INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR 
## PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE 
## FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
## ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
## SOFTWARE.
##
### Commentary:
##
## Perspecto: The CUDEM Visualization Engine.
##
## Perspecto automates the generation of high-quality 2D and 3D visualizations
## from DEMs. It wraps standard tools like GDAL, POV-Ray, and GMT to create
## aesthetic maps, hillshades, and data analysis plots.
##
##   * Hillshading & Relief:
##      - Generates standard, multidirectional, and color-relief hillshades.
##      - Manages Color Palette Tables (CPTs) automatically based on data range.
##
##   * 3D Rendering (Using POVRay):
##      - 'Perspective' module uses POV-Ray to create ray-traced 3D views.
##      - 'Sphere' module generates orthographic global views.
##
##   * Data Analysis:
##      - 'Histogram' generates hypsometric curves and CDF plots.
##
## Usage:
##   CLI: perspecto input.tif -M hillshade -C globe.cpt -O output.tif
##   API: p = PerspectoFactory(mod='perspective', src_dem='dem.tif').acquire()
##        p.run()
##
### Code:

import os
import sys
import json
import argparse

from cudem import utils
from cudem import gdalfun
from cudem import regions
from cudem import factory
from cudem import __version__ as __cudem_version__
from . import cpt
from . import __version__

try:
    import pygmt
    HAS_PYGMT = True
except ImportError:
    HAS_PYGMT = False

## ==============================================
## Main Perspecto Class
## ==============================================
class Perspecto:
    def __init__(
        self,
        mod=None,
        src_dem=None,
        cpt_file=None,
        min_z=None,
        max_z=None,
        callback=lambda: False,
        outfile=None,
        outdir=None,
        verbose=True,
        split_cpt=False,
        want_gdal_cpt=True,
        params=None
    ):
        self.mod = mod
        self.mod_args = {}
        self.src_dem = src_dem
        self.outfile = outfile
        self.cpt = cpt_file
        self.min_z = utils.float_or(min_z)
        self.max_z = utils.float_or(max_z)
        self.split_cpt = split_cpt
        self.outdir = outdir
        self.callback = callback
        self.verbose = verbose
        self.params = params if params is not None else {}
        
        self.dem_infos = gdalfun.gdal_infos(self.src_dem, scan=True)
        self.dem_region = regions.Region().from_geo_transform(
            self.dem_infos['geoT'], self.dem_infos['nx'], self.dem_infos['ny']
        )
        
        if self.outfile is None:
            self.outfile = f"{utils.fn_basename2(self.src_dem)}_{self.mod if self.mod else 'pp'}.tif"

        self.init_cpt(want_gdal=want_gdal_cpt)

        
    def __call__(self):
        return self.run()

    
    def init_cpt(self, want_gdal=False):
        """Initialize and process the CPT file."""
        
        min_z = self.min_z if self.min_z is not None else self.dem_infos['zr'][0]
        max_z = self.max_z if self.max_z is not None else self.dem_infos['zr'][1]


        if self.cpt is None:
            self.cpt = cpt.generate_etopo_cpt(min_z, max_z)
            
        elif os.path.exists(self.cpt):
            if self.verbose:
                utils.echo_msg(
                    f"processing cpt {self.cpt}, want_gdal is {want_gdal}, split_cpt: {self.split_cpt}"
                )
            self.cpt = cpt.process_cpt(
                self.cpt, min_z, max_z, gdal=want_gdal, split_cpt=self.split_cpt
            )
        else:
            ## Attempt to fetch city CPT
            self.cpt = cpt.process_cpt(
                cpt.fetch_cpt_city(q=self.cpt),
                min_z,
                max_z,
                gdal=want_gdal,
                split_cpt=self.split_cpt
            )

        ## If cpt is still None, generate the default ETOPO cpt
        if self.cpt is None:
            self.cpt = cpt.generate_etopo_cpt(min_z, max_z)                
            
    def makecpt(self, cmap='etopo1', color_model='r', output=None):
        min_z = self.min_z if self.min_z is not None else self.dem_infos['zr'][0]
        max_z = self.max_z if self.max_z is not None else self.dem_infos['zr'][1]
        self.cpt = cpt.make_cpt(
            cmap=cmap, 
            color_model=color_model, 
            output=output, 
            min_z=min_z, 
            max_z=max_z
        )

        
    def cpt_no_slash(self):
        """Removes backslashes from the CPT file."""
        
        try:
            with open(self.cpt, 'r') as file:
                filedata = file.read()
            
            filedata = filedata.replace('\\', ' ')
                
            with open(self.cpt, 'w') as file:
                file.write(filedata)
        except OSError as e:
            utils.echo_error_msg(f"Error processing CPT file: {e}")

            
    def export_as_png(self, rgb=True, dem=True):
        basename = utils.fn_basename2(self.src_dem)
        
        if dem:
            # Scale DEM to 16-bit PNG
            utils.run_cmd(
                f"gdal_translate -ot UInt16 -of PNG -scale "
                f"{self.dem_infos['zr'][0]} {self.dem_infos['zr'][1]} 0 65535 "
                f"{self.src_dem} _dem_temp.png",
                verbose=True
            )
            # Crop 1 pixel border
            utils.run_cmd(
                f"gdal_translate -srcwin 1 1 {self.dem_infos['nx']-1} {self.dem_infos['ny']-1} "
                f"-of PNG _dem_temp.png {basename}_16bit.png",
                verbose=True
            )
            utils.remove_glob('_dem_temp*')

        if rgb:
            self.init_cpt(want_gdal=True)
            # Generate Color Relief
            utils.run_cmd(
                f"gdaldem color-relief {self.src_dem} {self.cpt} _rgb_temp.tif -alpha",
                verbose=True
            )
            # Crop 1 pixel border
            utils.run_cmd(
                f"gdal_translate -srcwin 1 1 {self.dem_infos['nx']-1} {self.dem_infos['ny']-1} "
                f"-of PNG _rgb_temp.tif {basename}_rgb.png",
                verbose=True
            )
            utils.remove_glob('_rgb_temp*')
            r
            
    def run(self):
        raise NotImplementedError("Perspecto is an abstract base class; run() must be implemented by subclasses.")


## ==============================================
## Perspecto Factory
## ==============================================
class PerspectoFactory(factory.CUDEMFactory):
    from . import hillshade
    from . import perspective
    from . import sphere
    from . import figure1
    from . import colorbar
    from . import histogram
    #from . import joyplot

    _modules = {
        'hillshade': {'call': hillshade.Hillshade},
        'hillshade2': {'call': hillshade.Hillshade_cmd},
        'perspective': {'call': perspective.Perspective},
        'sphere': {'call': sphere.Sphere},
        #'joyplot': {'call': joyplot.Joyplot},
        'histogram': {'call': histogram.Histogram},
    }

    if HAS_PYGMT:
        _modules['figure1'] = {'call': figure1.Figure1}
        _modules['colorbar'] = {'call': colorbar.Colorbar}

        
    def __init__(self, **kwargs):
        super().__init__(**kwargs)


## ==============================================
## Command-line Interface (CLI)
## $ perspecto
##
## perspecto cli
## ==============================================
class PrintModulesAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        factory.echo_modules(PerspectoFactory._modules, values, md=True if not values else False)
        sys.exit(0)

        
def perspecto_cli():
    parser = argparse.ArgumentParser(
        description=f"%(prog)s ({__version__}): Generate iMAGEs from a DEM",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=f"""
Supported %(prog)s modules (see %(prog)s --modules <module-name> for more info): 
{factory.get_module_short_desc(PerspectoFactory._modules)}
        
CUDEM home page: <http://cudem.colorado.edu>
        """
    )
    
    parser.add_argument(
        'dem', 
        help="Input DEM file or Waffles config JSON"
    )
    
    parser.add_argument(
        '-M', '--module',
        default='hillshade',
        help="Desired perspecto MODULE and options.\n"
             "Format: module[:opt=val[:opt=val...]]"
    )
    
    parser.add_argument(
        '-C', '--cpt',
        help="Color Palette file (if not specified will auto-generate ETOPO CPT)"
    )
    
    parser.add_argument(
        '-Z', '--split-cpt',
        type=float,
        metavar='VALUE',
        help="Split the CPT values at specified value (usually zero)"
    )
    
    parser.add_argument(
        '--min_z', 
        type=float, 
        help="Minimum z value to use in CPT"
    )
    
    parser.add_argument(
        '--max_z', 
        type=float, 
        help="Maximum z value to use in CPT"
    )
    
    parser.add_argument(
        '-m', '--modules',
        nargs='?',
        default=None,
        action=PrintModulesAction,
        help="Display the module descriptions and usage"
    )
    
    parser.add_argument(
        '--version',
        action='version',
        version=f'CUDEM {__cudem_version__} :: %(prog)s {__version__}'
    )

    args = parser.parse_args()

    # --- Validation ---
    
    # Check for Module
    if not args.module:
        parser.print_usage()
        utils.echo_error_msg("must specify a perspecto -M module.")
        sys.exit(1)

    module_name = args.module.split(':')[0]
    if module_name not in PerspectoFactory._modules:
        utils.echo_error_msg(
            f"{module_name} is not a valid perspecto module, available modules are: "
            f"{factory._cudem_module_short_desc(PerspectoFactory._modules)}"
        )
        sys.exit(1)

    # Check for DEM input
    if not args.dem:
        parser.print_usage()
        utils.echo_error_msg(
            "you must supply a waffles config file or an existing DEM; "
            "see waffles --help for more information."
        )
        sys.exit(1)

    # --- Processing ---

    src_dem = None
    wg_user = args.dem

    if os.path.exists(wg_user):
        try:
            # Attempt to load as JSON (waffles config)
            with open(wg_user, 'r') as wgj:
                wg = json.load(wgj)
                
            if wg.get('src_region') is not None:
                wg['src_region'] = regions.Region().from_list(wg['src_region'])

            # Assuming 'waffles' is available in current scope or factory
            import cudem.waffles.waffles as waffles 
            
            this_waffle = waffles.WaffleFactory(**wg).acquire()
            this_waffle.mask = True
            this_waffle.clobber = False

            if not this_waffle.valid_p():
                this_waffle.generate()

            src_dem = this_waffle.fn
            
        except (json.JSONDecodeError, UnicodeDecodeError):
            # If JSON fails, assume it is a direct DEM path
            src_dem = wg_user
        except Exception as e:
            utils.echo_error_msg(f"Error processing input file: {e}")
            sys.exit(1)
    else:
        # If file doesn't exist on disk
        utils.echo_error_msg(
            f"specified waffles config file/DEM does not exist: {wg_user}"
        )
        sys.exit(1)

    # --- Execution ---

    this_perspecto = PerspectoFactory(
        mod=args.module,
        src_dem=src_dem,
        cpt_file=args.cpt,
        min_z=args.min_z,
        max_z=args.max_z,
        split_cpt=args.split_cpt
    )

    if this_perspecto is not None:
        this_perspecto_module = this_perspecto._acquire_module()
        if this_perspecto_module is not None:
            out_fig = this_perspecto_module()
            utils.echo_msg(f'generated {out_fig}')
        else:
            utils.echo_error_msg(f'could not acquire perspecto module {args.module}')

### End
