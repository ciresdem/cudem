### globato.py 
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
##
## globato.py is part of CUDEM
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
## GLOBATO BLOCKS
##
### Examples:
##
### TODO:
## add temporal
##
### Code:
###############################################################################

from cudem import utils
from cudem import factory

class GlobatoFile(ElevationDataset):

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

class Globato:

    stack_modes = [
        'min', 'max', 'mean', 'supercede', 'mixed'
    ]

    def __init__(
            self, src_region=None, x_inc=None, y_inc=None, dst_srs=None,
            verbose=True, cache_dir=None, stack_mode='mean', dst_fn=None,
    ):
        self.region = src_region
        self.x_inc = x_inc
        self.y_inc = y_inc
        self.dst_srs = dst_srs
        self.verbose = verbose
        self.cache_dir = cache_dir
        self.stack_mode = stack_mode
        
        self.block_ds = None
        self.datasets_data = {
            'z': None,
            'count': None,
            'weight': None,
            'uncertainty': None,
            'x': None,
            'y': None,
            #'datetime': None
        }
        self.blocked_data = self.datasets_data.copy()
        self.blocked_data['src_uncertainty'] = None
        self.sums_data = self.blocked_data.copy()

        self.stack_keys = list(self.blocked_data.keys())
        
        self.status = 0

        self._init_stack_mode()
        if 'mask_level' in self.stack_mode_args.keys():
            self.mask_level = utils.int_or(
                self.stack_mode_args['mask_level'], 0
            )
        else:
            self.mask_level = 0

        utils.set_cache(self.cache_dir)
        self.fn = os.path.join(self.cache_dir, '{}.h5'.format(
            utils.append_fn('globato', self.region, self.x_inc)
        ))

        utils.set_cache(self.cache_dir)

        self.dst_fn = dst_fn
        if self.dst_fn is None:
            self.dst_fn = os.path.join(
                self.cache_dir, utils.append_fn(
                    'globato', self.region, self.x_inc
                )
            )

        if self.verbose:
            utils.echo_msg(
                f'globato using {self.stack_mode_name} with {self.stack_mode_args} to {self.dst_fn)'
            )

            
    def _init_stack_mode(self):
        """from `self.stack_mode` will create:
        `self.stack_mode_name` and `self.stack_mode_args`
        """
        
        opts, self.stack_mode_name, self.stack_mode_args \
            = factory.parse_fmod(self.stack_mode)

        if self.stack_mode_name not in self.stack_modes:
            utils.echo_warning_msg(
                f'{self.stack_mode_name} is not a valid stack mode'
            )
            self.stack_mode_name = 'mean'
        
        if 'mask_level' not in self.stack_mode_args.keys():
            self.stack_mode_args['mask_level'] = -1


    def _load(self):
        raise(NotImplementedError)
            
    def average(self, weight_above, stacked_data, arrs):
        ## average of incoming data with existing data above weight_threshold
        stacked_data['count'][weight_above] += arrs['count'][weight_above]
        stacked_data['z'][weight_above] \
            += arrs['z'][weight_above]# * arrs['weight'][weight_above])
        stacked_data['x'][weight_above] \
            += arrs['x'][weight_above]# * arrs['weight'][weight_above])
        stacked_data['y'][weight_above] \
            += arrs['y'][weight_above]# * arrs['weight'][weight_above])
        stacked_data['src_uncertainty'][weight_above] \
            = np.sqrt(np.power(stacked_data['src_uncertainty'][weight_above], 2) \
                      + np.power(arrs['uncertainty'][weight_above], 2))
        stacked_data['weights'][weight_above] \
            += arrs['weight'][weight_above]
        
        ## accumulate variance * weight
        stacked_data['uncertainty'][weight_above] \
            += arrs['weight'][weight_above] \
            * np.power(
                (((arrs['z'][weight_above] / arrs['weight'][weight_above]) \
                  / arrs['count'][weight_above])\
                 - ((stacked_data['z'][weight_above] / stacked_data['weights'][weight_above]) \
                    / stacked_data['count'][weight_above])), 2)
        
        return(stacked_data)

    
    def supercede(self, weight_above_sup, stacked_data, arrs):
        stacked_data['count'][weight_above_sup] = arrs['count'][weight_above_sup]
        stacked_data['z'][weight_above_sup] \
            = arrs['z'][weight_above_sup]
        stacked_data['x'][weight_above_sup] \
            = arrs['x'][weight_above_sup]
        stacked_data['y'][weight_above_sup] \
            = arrs['y'][weight_above_sup]
        stacked_data['src_uncertainty'][weight_above_sup] \
            = arrs['uncertainty'][weight_above_sup]
        stacked_data['weights'][weight_above_sup] \
            = arrs['weight'][weight_above_sup]
        stacked_data['uncertainty'][weight_above_sup] \
            = np.array(stacked_data['src_uncertainty'][weight_above_sup])

                    
            return(stacked_data)

        
            
    def accumulate_sums_from_array(self, arrs, srcwin, gt):            
        raise(NotImplementedError)


    def finalize_arrs(self, stacked_data, srcwin):
        #######################################################################
        ## finalize the stacked_data weighted sums
        #######################################################################
        for y in range(
                srcwin[1], srcwin[1] + srcwin[3], 1
        ):
            for key in self.stack_keys:
                stacked_data[key] = stacked_bands[key].ReadAsArray(
                    srcwin[0], y, srcwin[2], 1
                )
                stacked_data[key][stacked_data[key] == ndv] = np.nan

            stacked_data['weights'][stacked_data['weights'] == 0] = 1
            #if mode != 'supercede':
            stacked_data['weights'] = stacked_data['weights'] / stacked_data['count']
            #if mode == 'mean' or mode == 'mixed':
            ## average the accumulated arrays for finalization
            ## x, y, z and u are weighted sums, so divide by weights
            stacked_data['x'] = (stacked_data['x'] / stacked_data['weights']) \
                / stacked_data['count']
            stacked_data['y'] = (stacked_data['y'] / stacked_data['weights']) \
                / stacked_data['count']
            stacked_data['z'] = (stacked_data['z'] / stacked_data['weights']) \
                / stacked_data['count']

            ## apply the source uncertainty with the sub-cell variance uncertainty
            ## caclulate the standard error (sqrt( uncertainty / count))
            stacked_data['uncertainty'] = np.sqrt(
                (stacked_data['uncertainty'] / stacked_data['weights']) \
                / stacked_data['count']
            )
            stacked_data['uncertainty'] = np.sqrt(
                np.power(stacked_data['src_uncertainty'], 2) \
                + np.power(stacked_data['uncertainty'], 2)
            )

    
class Stacks(Globato):
    """ multi-banded gdal raster
    """
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)


class Blocks(Globato):
    """h5 output
    """
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
