### perspecto_cli.py 
##
## Copyright (c) 2023 - 2025 Regents of the University of Colorado
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
##
### Code:

import os, sys
from cudem.perspecto import *
from cudem import utils
from cudem import factory

try:
    import pygmt
    has_pygmt = True
except ImportError as e:
    has_pygmt = False

class PerspectoFactory(factory.CUDEMFactory):
    if has_pygmt:
        _modules = {
            'hillshade': {'call': hillshade.Hillshade},
            'hillshade2': {'call': hillshade.Hillshade_cmd},
            'perspective': {'call': perspective.perspective},
            'sphere': {'call': sphere.sphere},
            'figure1': {'call': figure1.figure1},
            'colorbar': {'call': colorbar.colorbar},
        }
    else:
        _modules = {
            'hillshade': {'call': hillshade.Hillshade},
            'hillshade2': {'call': hillshade.Hillshade_cmd},
            'perspective': {'call': perspective.perspective},
            'sphere': {'call': sphere.sphere},
        }    
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
           d_formats=factory._cudem_module_short_desc(PerspectoFactory._modules))
        
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
            split_cpt = 0
        
        elif arg == '--modules' or arg == '-m':
            factory.echo_modules(
                PerspectoFactory._modules, None if i+1 >= len(argv) else sys.argv[i+1]
            )
            sys.exit(0)            
        elif arg == '--help' or arg == '-h':
            sys.stderr.write(perspecto_cli_usage())
            sys.exit(0)
        elif arg == '--version' or arg == '-v':
            sys.stdout.write('{}\n'.format(cudem.__version__))
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
