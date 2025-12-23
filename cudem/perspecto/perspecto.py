### perspecto.py 
##
## Copyright (c) 2023 - 2025 Regents of the University of Colorado
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
###############################################################################
### Commentary:
##
## Generate iMAGEs from a DEm
##
## uses:
##   povray
##   gdal
##   ImageMagick
##   GMT
##
### Code:

import os, sys
import math

from cudem import utils
from cudem import gdalfun
from cudem import regions
from cudem import factory
from . import cpt
from . import __version__

try:
    import pygmt
    has_pygmt = True
except ImportError as e:
    has_pygmt = False
    

class Perspecto:
    def __init__(
            self,
            mod=None,
            src_dem=None,
            cpt=None,
            min_z=None,
            max_z=None,
            callback=lambda: False,
            outfile=None,
            outdir=None,
            verbose=True,
            split_cpt=False,
            params={}
    ):
        self.mod = mod
        self.mod_args = {}
        self.src_dem = src_dem
        self.outfile = outfile
        self.cpt = cpt
        self.min_z = utils.float_or(min_z)
        self.max_z = utils.float_or(max_z)
        self.split_cpt = split_cpt
        self.outdir = outdir
        self.callback = callback
        self.verbose = verbose
        self.params = params    
        self.dem_infos = gdalfun.gdal_infos(self.src_dem, scan=True)
        self.dem_region = regions.Region().from_geo_transform(
            self.dem_infos['geoT'], self.dem_infos['nx'], self.dem_infos['ny']
        )
        if self.outfile is None:
            self.outfile = '{}_{}.tif'.format(
                utils.fn_basename2(self.src_dem), self.mod
            )
        
        self.init_cpt(want_gdal=True)

        
    def __call__(self):
        return(self.run())

    
    def init_cpt(self, want_gdal=False):
        min_z = self.min_z if self.min_z is not None else self.dem_infos['zr'][0]
        max_z = self.max_z if self.max_z is not None else self.dem_infos['zr'][1]
        if self.cpt is None:
            # self.makecpt(
            #     'etopo1',
            #     output='{}_etopo1.cpt'.format(utils.fn_basename2(self.src_dem))
            # )
            self.cpt = cpt.generate_etopo_cpt(min_z, max_z)
        #else:
        #self.cpt = process_cpt(self.cpt, min_z, max_z)
        elif os.path.exists(self.cpt):
            #if has_pygmt:
            #    self.makecpt(cmap=self.cpt, color_model='r', output='{}.cpt'.format(utils.fn_basename2(self.src_dem)))
            #else:
            utils.echo_msg(
                'processing cpt {}, want_gdal is {}, split_cpt: {}'.format(
                    self.cpt, want_gdal, self.split_cpt
                )
            )
            self.cpt = cpt.process_cpt(
                self.cpt, min_z, max_z, gdal=want_gdal, split_cpt=self.split_cpt
            )
        else:
            self.cpt = cpt.process_cpt(
                cpt.fetch_cpt_city(q=self.cpt),
                min_z,
                max_z,
                gdal=want_gdal,
                split_cpt=self.split_cpt
            )
        
    
    def makecpt(self, cmap='etopo1', color_model='r', output=None):
        min_z = self.min_z if self.min_z is not None else self.dem_infos['zr'][0]
        max_z = self.max_z if self.max_z is not None else self.dem_infos['zr'][1]
        self.cpt = cpt.makecpt(cmap = cmap, color_model = color_model, output = output, min_z = min_z, max_z = max_z)

        
    def cpt_no_slash(self):
        # Read in the file
        with open(self.cpt, 'r') as file:
            filedata = file.read()
            
        # Replace the target string
        filedata = filedata.replace('\\', ' ')
            
        # Write the file out again
        with open(self.cpt, 'w') as file:
            file.write(filedata)
            
    
    def export_as_png(self, rgb=True, dem=True):
        if dem:
            utils.run_cmd(
                'gdal_translate -ot UInt16 -of PNG -scale {} {} 0 65535 {} _dem_temp.png'.format(
                    self.dem_infos['zr'][0],
                    self.dem_infos['zr'][1],
                    self.src_dem
                ),
                verbose=True
            )
            utils.run_cmd(
                'gdal_translate -srcwin 1 1 {} {} -of PNG _dem_temp.png {}_16bit.png'.format(
                    self.dem_infos['nx']-1,
                    self.dem_infos['ny']-1,
                    utils.fn_basename2(self.src_dem)
                ),
                verbose=True
            )
            utils.remove_glob('_dem_temp*')
            
        if rgb:
            self.init_cpt(want_gdal=True)
            utils.run_cmd(
                'gdaldem color-relief {} {} _rgb_temp.tif -alpha'.format(
                    self.src_dem, self.cpt
                ),
                verbose=True
            )
            utils.run_cmd(
                'gdal_translate -srcwin 1 1 {} {} -of PNG _rgb_temp.tif {}_rgb.png'.format(
                    self.dem_infos['nx']-1,
                    self.dem_infos['ny']-1,
                    utils.fn_basename2(self.src_dem)
                ),
                verbose=True
            )
            utils.remove_glob('_rgb_temp*')

            
    def run(self):
        raise(NotImplementedError)


class PerspectoFactory(factory.CUDEMFactory):
    from . import hillshade
    from . import perspective
    from . import sphere
    from . import figure1
    from . import colorbar

    _modules = {
        'hillshade': {'call': hillshade.Hillshade},
        'hillshade2': {'call': hillshade.Hillshade_cmd},
        'perspective': {'call': perspective.perspective},
        'sphere': {'call': sphere.sphere},
    }

    if has_pygmt:
        _modules['figure1'] = {'call': figure1.figure1}
        _modules['colorbar'] = {'call': colorbar.colorbar}
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        

## ==============================================
## Command-line Interface (CLI)
## $ perspecto
##
## perspecto cli
## ==============================================
perspecto_cli_usage = lambda: """{cmd}

usage: {cmd} [ -hvCMZ [ args ] ] DEM ...

Options:

  -C, --cpt\t\tColor Pallette file (if not specified will auto-generate ETOPO CPT)
  -M, --module\t\tDesired perspecto MODULE and options. (see available Modules below)
\t\t\tWhere MODULE is module[:mod_opt=mod_val[:mod_opt1=mod_val1[:...]]]
  -Z, --split-cpt\tSplit the CPT values at zero

  --min_z\t\tMinimum z value to use in CPT
  --max_z\t\tMaximum z value to use in CPT

  --help\t\tPrint the usage text
  --modules\t\tDisplay the module descriptions and usage
  --version\t\tPrint the version information

Supported PERSPECTO modules (see perspecto --modules <module-name> for more info): 
  {d_formats}
""".format(cmd=os.path.basename(sys.argv[0]),
           d_formats=factory.get_module_short_desc(PerspectoFactory._modules))
        
#if __name__ == '__main__':
def perspecto_cli(argv = sys.argv):
    i = 1
    src_dem = None
    wg_user = None
    module = None
    src_cpt = None
    min_z = None
    max_z = None
    split_cpt = None
    
    while i < len(argv):
        arg = argv[i]
        if arg == '--module' or arg == '-M':
            module = str(argv[i + 1])
            i += 1
        elif arg[:2] == '-M':
            module = str(arg[2:])

        elif arg == '--cpt' or arg == '-C':
            src_cpt = str(argv[i + 1])
            i += 1
        elif arg[:2] == '-C':
            src_cpt = str(arg[2:])

        elif arg == '--min_z':
            min_z = utils.float_or(argv[i + 1])
            i += 1
            
        elif arg == '--max_z':
            max_z = utils.float_or(argv[i + 1])
            i += 1
        elif arg == '--split-cpt' or arg == '-Z':
            split_cpt = utils.float_or(argv[i + 1])
            i += 1
        
        elif arg == '--modules' or arg == '-m':
            factory.echo_modules(
                PerspectoFactory._modules, None if i+1 >= len(argv) else sys.argv[i+1]
            )
            sys.exit(0)            
        elif arg == '--help' or arg == '-h':
            sys.stderr.write(perspecto_cli_usage())
            sys.exit(0)
        elif arg == '--version' or arg == '-v':
            sys.stdout.write('{}\n'.format(__version__))
            sys.exit(0)
        elif arg[0] == '-':
            sys.stdout.write(perspecto_cli_usage())
            utils.echo_error_msg('{} is not a valid perspecto cli switch'.format(arg))
            sys.exit(0)
        else:
            wg_user = arg
        i += 1

    if module is None:
        sys.stderr.write(perspecto_cli_usage())
        utils.echo_error_msg(
            '''must specify a perspecto -M module.'''
        )
        sys.exit(-1)

    if module.split(':')[0] not in PerspectoFactory()._modules.keys():
        utils.echo_error_msg(
            '''{} is not a valid perspecto module, available modules are: {}'''.format(
                module.split(':')[0], factory._cudem_module_short_desc(PerspectoFactory._modules)
            )
        )
        sys.exit(-1)
        
    ## load the user wg json and run perspecto with that.
    if wg_user is not None:
        if os.path.exists(wg_user):
            try:
                with open(wg_user, 'r') as wgj:
                    wg = json.load(wgj)
                    if wg['src_region'] is not None:
                        wg['src_region'] = regions.Region().from_list(
                            wg['src_region']
                        )

                    this_waffle = waffles.WaffleFactory(**wg).acquire()
                    this_waffle.mask = True
                    this_waffle.clobber = False

                    if not this_waffle.valid_p():
                        this_waffle.generate()

                    src_dem = this_waffle.fn
            except:
                src_dem = wg_user
        else:
            sys.stderr.write(perspecto_cli_usage())
            utils.echo_error_msg(
                'specified waffles config file/DEM does not exist, {}'.format(wg_user)
            )
            sys.exit(-1)
    else:
        sys.stderr.write(perspecto_cli_usage())
        utils.echo_error_msg(
            ('you must supply a waffles config file or an existing DEM; '
             'see waffles --help for more information.')
        )
        sys.exit(-1)

    this_perspecto = PerspectoFactory(
        mod=module,
        src_dem=src_dem,
        cpt=src_cpt,
        min_z=min_z,
        max_z=max_z,
        split_cpt=split_cpt
    )
    if this_perspecto is not None:
        this_perspecto_module = this_perspecto._acquire_module()
        if this_perspecto_module is not None:
            out_fig = this_perspecto_module()
            utils.echo_msg(f'generated {out_fig}')
            #this_perspecto_module.run()
        else:
            utils.echo_error_msg(
                f'could not acquire perspecto module {module}'
            )

            
### End
