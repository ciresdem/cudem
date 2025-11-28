### grits.py - GRId filTerS
##
## Copyright (c) 2024 - 2025 Regents of the University of Colorado
##
## grits.py is part of cudem
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
## filter grids using various methods
##
## Grits modules are sub-classes of the Grits class.
## Define a new sub-class to create a new DEM filter.
##
### Code:

import os, sys
import traceback
import numpy as np
from osgeo import gdal

from cudem import utils
from cudem import gdalfun
from cudem import factory
from . import __version__

class Grits:
    """DEM Filtering.

    Filter DEMs using various filtering methods. 

    Define a sub-class to make a new grits filter.
    """
    
    def __init__(
            self,
            src_dem: str = None,
            dst_dem: str = None,
            band: int = 1,
            min_z: float = None,
            max_z: float = None,
            min_weight: float = None,
            max_weight: float = None,
            count_mask: any = None,
            weight_mask: any = None,
            uncertainty_mask: any = None,
            cache_dir: str = './',
            verbose: bool = True,
            params: dict = {},
            **kwargs: any
    ):
        self.src_dem = src_dem
        self.dst_dem = dst_dem
        self.band = utils.int_or(band, 1)
        self.min_z = utils.float_or(min_z)
        self.max_z = utils.float_or(max_z)
        self.min_weight = utils.float_or(min_weight)
        self.max_weight = utils.float_or(max_weight)
        self.count_mask = count_mask
        self.weight_mask = weight_mask
        self.uncertainty_mask = uncertainty_mask
        self.cache_dir = cache_dir
        self.verbose = verbose
        self.params = params
        self.kwargs = kwargs

        if self.dst_dem is None:
            if self.src_dem is not None:
                self.dst_dem = utils.make_temp_fn('{}_filtered.{}'.format(
                    utils.fn_basename2(self.src_dem),
                    utils.fn_ext(self.src_dem)
                ), temp_dir=self.cache_dir)
            else:
                self.dst_dem = 'grits_filtered.tif'

                
    def __call__(self):
        return(self.generate())

    
    def init_ds(self, src_ds: any = None):
        self.ds_config = gdalfun.gdal_infos(src_ds)
        self.ds_band = src_ds.GetRasterBand(self.band)
        self.gt = self.ds_config['geoT']

        if self.ds_band.GetNoDataValue() is None:
            self.ds_band.SetNoDataValue(self.ds_config['ndv'])
        
        ## setup the associted mask data (uncertainty, weights, counts)
        self.weight_is_fn = False
        self.weight_is_band = False
        self.uncertainty_is_fn = False
        self.uncertainty_is_band = False
        self.count_is_fn = False
        self.count_is_band = False
        if self.weight_mask is not None:
            self.weight_is_band = False
            self.weight_is_fn = False
            if utils.int_or(self.weight_mask) is not None:
                self.weight_is_band = True
                self.weight_mask = utils.int_or(self.weight_mask)
            elif os.path.exists(self.weight_mask):
                self.weight_is_fn = True
            else:
                self.weight_mask = None        

        if self.uncertainty_mask is not None:
            self.unc_is_band = False
            self.unc_is_fn = False
            if utils.int_or(self.uncertainty_mask) is not None:
                self.unc_is_band = True
                self.uncertainty_mask = utils.int_or(self.uncertainty_mask)
            elif os.path.exists(self.uncertainty_mask):
                self.unc_is_fn = True
            else:
                self.uncertainty_mask = None
                
        if self.count_mask is not None:
            self.cnt_is_band = False
            self.cnt_is_fn = False
            if utils.int_or(self.count_mask) is not None:
                self.cnt_is_band = True
                self.count_mask = utils.int_or(self.count_mask)
            elif os.path.exists(self.count_mask):
                self.cnt_is_fn = True
            else:
                self.count_mask = None

                
    def generate(self):
        if self.verbose:
            utils.echo_msg(
                f'filtering {self.src_dem} using {self}'
            )
            
        self.run()
        self.split_by_z()
        self.split_by_weight()
        return(self)        

    
    def run(self):
        raise(NotImplementedError)

    
    def copy_src_dem(self):
        with gdalfun.gdal_datasource(
                self.src_dem, update=False
        ) as src_ds:
            if src_ds is not None:
                src_infos = gdalfun.gdal_infos(src_ds)
                driver = gdal.GetDriverByName(src_infos['fmt'])
                copy_ds = driver.CreateCopy(
                    self.dst_dem, src_ds, 1, options=['COMPRESS=DEFLATE']
                )
            else:
                copy_ds = None
                
        return(copy_ds)

    
    def _density(self, src_arr):
        nonzero = np.count_nonzero(~np.isnan(src_arr))
        dd = nonzero / src_arr.size
        
        return(dd)

    
    def split_by_z(self):
        """Split the filtered DEM by z-value"""

        if self.max_z is not None or self.min_z is not None:
            utils.echo_msg(
                f'split by z:{self.min_z} {self.max_z}'
            )
            with gdalfun.gdal_datasource(self.src_dem) as src_ds:
                if src_ds is not None:
                    self.init_ds(src_ds)
                    elev_array = self.ds_band.ReadAsArray()
                    mask_array = np.zeros(
                        (self.ds_config['ny'], self.ds_config['nx'])
                    )
                    mask_array[elev_array == self.ds_config['ndv']] = np.nan
                    if self.min_z is not None:
                        mask_array[elev_array > self.min_z] = 1
                        if self.max_z is not None:
                            mask_array[elev_array > self.max_z] = 0
                        
                    elif self.max_z is not None:
                        mask_array[elev_array < self.max_z] = 1
                        if self.min_z is not None:
                            mask_array[elev_array < self.min_z] = 0
                        
                    elev_array[mask_array == 1] = 0

                    ## todo: all bands
                    with gdalfun.gdal_datasource(self.dst_dem, update=True) as s_ds:
                        if s_ds is not None:
                            #for b in range(1, s_ds.RasterCount+1):
                            s_band = s_ds.GetRasterBand(1)
                            s_array = s_band.ReadAsArray()
                            s_array = s_array * mask_array
                            smoothed_array = s_array + elev_array
                            elev_array = None
                            s_band.WriteArray(smoothed_array)
                            
        return(self)

    
    def split_by_weight(self):
        """Split the filtered DEM by z-value"""

        if self.max_weight is not None \
           or self.min_weight is not None:
            if self.weight_mask is not None:
                utils.echo_msg(
                    f'split by weight: {self.min_weight} {self.max_weight}'
                )
                with gdalfun.gdal_datasource(self.src_dem) as src_ds:
                    if src_ds is not None:
                        self.init_ds(src_ds)
                        elev_array = self.ds_band.ReadAsArray()

                        # uncertainty ds
                        weight_band = None
                        if self.weight_is_fn:
                            weight_ds = gdal.Open(self.weight_mask)
                            weight_band = weight_ds.GetRasterBand(1)
                        elif self.weight_is_band:
                            weight_band = src_ds.GetRasterBand(self.weight_mask)

                        if weight_band is not None:
                            weight_array = weight_band.ReadAsArray()
                            weight_array[(weight_array == self.ds_config['ndv'])] = 0

                        mask_array = np.zeros(
                            (self.ds_config['ny'], self.ds_config['nx'])
                        )
                        mask_array[elev_array == self.ds_config['ndv']] = np.nan
                        mask_array[weight_array == self.ds_config['ndv']] = np.nan

                        if self.min_weight is not None:
                            mask_array[weight_array > self.min_weight] = 1
                            if self.max_weight is not None:
                                mask_array[weight_array > self.max_weight] = 0

                        elif self.max_weight is not None:
                            mask_array[weight_array < self.max_weight] = 1
                            if self.min_weight is not None:
                                mask_array[weight_array < self.min_weight] = 0

                        elev_array[mask_array == 1] = 0

                        ## todo: all bands
                        with gdalfun.gdal_datasource(
                                self.dst_dem, update=True
                        ) as s_ds:
                            if s_ds is not None:
                                s_band = s_ds.GetRasterBand(1)
                                s_array = s_band.ReadAsArray()
                                s_array = s_array * mask_array
                                smoothed_array = s_array + elev_array
                                elev_array = None
                                s_band.WriteArray(smoothed_array)
        return(self)

    
    def get_outliers(self, in_array: any, percentile: float = 75,
                     k: float = 1.5, verbose: bool = False):
        """get the outliers from in_array based on the percentile

        https://en.wikipedia.org/wiki/Interquartile_range
        """

        if verbose:
            utils.echo_msg(
                f'input percentile: {percentile}'
            )

        # if np.isnan(percentile):
        #     percentile = 75
            
        if percentile < 0:
            percentile = 0
            
        if percentile > 100:
            percentile = 100

        max_percentile = percentile
        #min_percentile = 100 - percentile
        min_percentile = percentile-50

        if min_percentile < 0:
            min_percentile = 0

        if verbose:
            utils.echo_msg(
                f'percentiles: {min_percentile}>>{max_percentile}'
            )

        if np.all(np.isnan(in_array)):
            upper_limit = np.nan
            lower_limit = np.nan
        else:
            perc_max = np.nanpercentile(in_array, max_percentile)
            perc_min = np.nanpercentile(in_array, min_percentile)
            iqr_p = (perc_max - perc_min) * k
            upper_limit = perc_max + iqr_p
            lower_limit = perc_min - iqr_p

        return(upper_limit, lower_limit)


class GritsFactory(factory.CUDEMFactory):
    """Grits Factory Settings and Generator
    
    Parse a grits module and return the filtering object
    """

    from . import blur
    from . import gmtfilter
    from . import lspoutliers
    from . import weights
    from . import flats
    
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
           version=__version__,
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
            sys.stdout.write('{}\n'.format(__version__))
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
