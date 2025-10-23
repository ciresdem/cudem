### grits_factory.py - GRId filTerS
##
## Copyright (c) 2024 - 2025 Regents of the University of Colorado
##
## Permission is hereby granted, free of charge, to any person obtaining a copy 
## of this software and associated documentation files (the "Software"), to deal 
## in the Software without restriction, including without limitation the rights 
## to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies 
## of the Software, and to permit persons to whom the Software is furnished to do so, 
## subject to the following conditions:
##
## The above copyright notice and this permission notice shall be included in all copies or
## substantial portions of the Software.
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
import traceback
#import cudem
from cudem import factory
from cudem import utils
from cudem import grits
from cudem.grits import blur
from cudem.grits import gmtfilter
from cudem.grits import lspoutliers
from cudem.grits import weights
from cudem.grits import flats

class GritsFactory(factory.CUDEMFactory):
    """Grits Factory Settings and Generator
    
    Parse a grits module and return the filtering object
    """
    
    _modules = {
        'blur': {
            'name': 'blur',
            'description': 'Filter with a Gaussian Blur',
            'call': blur.Blur
        },
        'grdfilter': {
            'name': 'grdfilter',
            'description': 'Filter using GMTs `grdfilter` command',
            'call': gmtfilter.GMTgrdfilter
        },
        'outliers': {
            'name': 'outliers',
            'description': 'Remove outliers from the DEM',
            'call': lspoutliers.LSPOutliers
        },
        'lsp': {
            'name': 'lsp-outliers',
            'description': 'Remove outliers from the DEM',
            'call': lspoutliers.LSPOutliers
        },
        'flats': {
            'name': 'flats',
            'description': 'Remove flat areas from the DEM',
            'call': flats.Flats
        },
        'weights': {
            'name': 'weights',
            'description': 'Make a NDV buffer around the weight threshold',
            'call': weights.Weights
        },
        'weight_zones': {
            'name': 'weight_zones',
            'description': 'Make a NDV buffer around the weight threshold',
            'call': weights.WeightZones
        },
    }

    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)


###############################################################################        
## Command-line Interface (CLI)
## $ grits
##
## grits cli
###############################################################################
grits_cli_usage = """{cmd} ({version}): grits; GRId filTerS

usage: {cmd} [ -hvCMNUWX [ args ] ] DEM ...

Options:

  -M, --module\t\t\tDesired grits MODULE and options. (see available Modules below)
\t\t\t\tWhere MODULE is module[:mod_opt=mod_val[:mod_opt1=mod_val1[:...]]]
\t\t\t\tThis can be set multiple times to perform multiple filters.

  -N, --min_z\t\t\tMinimum z value (filter data above this value)
  -X, --max_z\t\t\tMaximum z value (filter data below this value)
  -Wn, --min_weight\t\tMinimum weight value (filter data above this value)
  -Wx, --max_weight\t\tMaximum weight value (filter data below this value)
  -U, --uncertainty_mask\tAn associated uncertainty raster or band number
  -W, --weight_mask\t\tAn associated weight raster or band number
  -C, --count_mask\t\tAn associated count raster or band number

  --help\t\t\tPrint the usage text
  --modules\t\t\tDisplay the module descriptions and usage
  --version\t\t\tPrint the version information

Supported GRITS modules (see grits --modules <module-name> for more info): 
  {d_formats}

Examples:
  % {cmd} input_dem.tif -M blur
  % {cmd} input_dem.tif --uncertainty_mask input_dem_u.tif --max_z 0 -M outliers:percentile=65
""".format(cmd=os.path.basename(sys.argv[0]),
           version=grits.__version__,
           d_formats=factory._cudem_module_short_desc(GritsFactory._modules))
        
#if __name__ == '__main__':
def grits_cli(argv = sys.argv):
    i = 1
    src_dem = None
    dst_dem = None
    wg_user = None
    min_z = None
    max_z = None
    min_weight = None
    max_weight = None
    uncertainty_mask = None
    weight_mask = None
    count_mask = None
    filters = []
    
    while i < len(argv):
        arg = argv[i]
        if arg == '--module' or arg == '-M':
            module = str(argv[i + 1])
            if module.split(':')[0] not in GritsFactory()._modules.keys():
                utils.echo_warning_msg(
                    '''{} is not a valid grits module, available modules are: {}'''.format(
                        module.split(':')[0],
                        factory._cudem_module_short_desc(GritsFactory._modules)
                    )
                )
            else:
                filters.append(module)
                
            i += 1
        elif arg[:2] == '-M':
            module = str(arg[2:])
            if module.split(':')[0] not in GritsFactory()._modules.keys():
                utils.echo_warning_msg(
                    '''{} is not a valid grits module, available modules are: {}'''.format(
                        module.split(':')[0],
                        factory._cudem_module_short_desc(GritsFactory._modules)
                    )
                )
            else:
                filters.append(module)
                
        elif arg == '--min_z' or arg == '-N':
            min_z = utils.float_or(argv[i + 1])
            i += 1
        elif arg[:2] == '-N':
            min_z = utils.float_or(arg[2:])
        elif arg == '--max_z' or arg == '-X':
            max_z = utils.float_or(argv[i + 1])
            i += 1
        elif arg == '--min_weight' or arg == '-Wn':
            min_weight = utils.float_or(argv[i + 1])
            i += 1
        elif arg[:3] == '-Wn':
            min_weight = utils.float_or(arg[3:])
        elif arg == '--max_weight' or arg == '-Wx':
            max_weight = utils.float_or(argv[i + 1])
            i += 1
        elif arg[:3] == 'Wx':
            max_weight = utils.float_or(arg[3:])            
        elif arg == '--uncertainty_mask' or arg == '-U':
            uncertainty_mask = argv[i + 1]
            i += 1
        elif arg[:2] == '-U':
            uncertainty_mask = arg[2:]
        elif arg == '--weight_mask' or arg == '-W':
            weight_mask = argv[i + 1]
            i += 1
        elif arg[:2] == '-W':
            weight_mask = arg[2:]
        elif arg == '--count_mask' or arg == '-C':
            count_mask = argv[i + 1]
            i += 1
        elif arg[:2] == '-C':
            count_mask = arg[2:]
        elif arg == '--outfile' or arg == '-O':
            dst_dem = argv[i + 1]
            i += 1
        elif arg[:2] == '-O':
            dst_dem = arg[2:]
        
        elif arg == '--modules' or arg == '-m':
            factory.echo_modules(
                GritsFactory._modules, None if i+1 >= len(argv) else sys.argv[i+1]
            )
            sys.exit(0)            
        elif arg == '--help' or arg == '-h':
            sys.stderr.write(grits_cli_usage)
            sys.exit(0)
        elif arg == '--version' or arg == '-v':
            sys.stdout.write('{}\n'.format(cudem.__version__))
            sys.exit(0)
        elif arg[0] == '-':
            sys.stdout.write(grits_cli_usage)
            utils.echo_error_msg(
                f'{arg} is not a valid grits cli switch'
            )
            sys.exit(0)
        else:
            wg_user = arg
        i += 1

    #if module is None:
    if len(filters) == 0:
        sys.stderr.write(grits_cli_usage)
        utils.echo_error_msg(
            '''must specify a grits -M module.'''
        )
        sys.exit(-1)
        
    ## load the user wg json and run grits with that.
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
            sys.stderr.write(grits_cli_usage)
            utils.echo_error_msg(
                f'specified waffles config file/DEM does not exist, {wg_user}'
            )
            sys.exit(-1)
    else:
        sys.stderr.write(grits_cli_usage)
        utils.echo_error_msg(
            ('you must supply a waffles config file or an existing DEM; '
             'see waffles --help for more information.')
        )
        sys.exit(-1)

    if dst_dem is None:
        dst_dem = utils.make_temp_fn('{}_filtered.{}'.format(
            utils.fn_basename2(src_dem), utils.fn_ext(src_dem)
        ))
        
    # src_dem = dst_dem        
    for module in filters:
        if module.split(':')[0] not in GritsFactory()._modules.keys():
            utils.echo_error_msg(
                '''{} is not a valid grits module, available modules are: {}'''.format(
                    module.split(':')[0], factory._cudem_module_short_desc(GritsFactory._modules)
                )
            )
            continue
        
        this_grits = GritsFactory(
            mod=module,
            src_dem=src_dem,
            min_z=min_z,
            max_z=max_z,
            min_weight=min_weight,
            max_weight=max_weight,
            uncertainty_mask=uncertainty_mask,
            weight_mask=weight_mask,
            count_mask=count_mask,
            dst_dem=dst_dem
        )
        try:
            this_grits_module = this_grits._acquire_module()
            if this_grits_module is not None:
                out_dem = this_grits_module()
                utils.echo_msg(f'filtered DEM to {out_dem.dst_dem}')
                #os.replace(out_dem.dst_dem, src_dem)
            else:
                utils.echo_error_msg(
                    f'could not acquire grits module {module}'
                )
                
        except KeyboardInterrupt:
            utils.echo_error_msg('Killed by user')
            
        except Exception as e:
            utils.echo_error_msg(e)
            print(traceback.format_exc())

### End
