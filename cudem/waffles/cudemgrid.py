### cudemgrid.py
##
## Copyright (c) 2015 - 2025 Regents of the University of Colorado
##
## cudemgrid.py is part of CUDEM
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
### Code:

import os

from cudem import utils
from cudem import gdalfun
from cudem import factory
from cudem.waffles.waffles import Waffle

class WafflesCUDEM(Waffle):
    """CUDEM integrated DEM generation.
    
    Generate an topo/bathy integrated DEM using a variety of data sources.
    Will iterate <pre_count> pre-surfaces at lower-resolutions, 
    specified with <inc_levels>, e.g. .3s/1s

    Each pre-surface will be clipped to <landmask>, if it exists, and smoothed 
    with <pre_smoothing> factor

    Each pre-surface is used in subsequent pre-surface(s)/final DEM at each 
    iterative weight specified in <weight_levels>

    generate a DEM with `pre_surface`s which are generated at lower resolution(s) 
    and with various weight threshholds

    To generate a typical CUDEM tile, generate 1 pre-surface ('bathy_surface'), 
    clipped to a coastline

    Use a <weight_levels>, e.g. 1/.5 that excludes low-resolution bathymetry data 
    from being used as input in the final DEM generation

    e.g.
    cudem:pre_mode=gmt-surface:tension=.65:geographic=True:landmask=True:polygonize=2:weight_levels=1/.5:inc_levels=.3333333s/1s:pre_count=2:pre_upper_limit=-.1:pre_smoothing=2:flatten=95

    Parameters:
    -----------    
    pre_mode (str): the waffles module to perform the initial pre-surface
    pre_count (int): number of pre-surface iterations to perform
    weight_levels (): the minimum weights for each of the pre-surface iterations
    inc_levels ():  the increments for each of the pre-surface iterations
    pre_upper_limit (float): the upper elevation limit of the pre-surfaces 
                              (used with landmask)
    pre_smoothing (float): the smoothing (blur) factor to apply to the pre-surfaces
    pre_verbose (bool): increase the verbosity of pre-surface generation
    landmask (bool): path to coastline vector mask or set as `coastline` to auto-generate
    invert_landmask (bool): invert the coastline mask, invert assumes mask over ocean, while
                            non invert assumes mask over land.
    want_supercede (bool): supercede subsquent pre-surfaces
    flatten (float): the nodata-size percentile above which to flatten
    exclude_lakes (bool): exclude lakes from the landmask
    filter_outliers (float): percentile at which to filter outliers
    <mode-opts>: options for the waffles module specified in 'mode'
    <coastline-opts>: options for the coastline module when generating the landmask
    """

    def __init__(
            self,
            pre_mode='gmt-surface',
            pre_count=2,
            pre_upper_limit=-0.1,
            pre_smoothing=None,
            weight_levels=None,
            inc_levels=None,
            landmask=True,
            invert_landmask=True,
            filter_outliers=None,
            want_supercede=False,
            flatten=None,
            exclude_lakes=False,
            mode=None,
            min_weight='1',
            pre_verbose=False,
            final_mode='IDW',
            keep_pre_surfaces=False,
            **kwargs
    ):
        from .coastline import WafflesCoastline
        from .waffles import WaffleFactory
        
        self.valid_modes = [
            'gmt-surface', 'IDW', 'linear', 'cubic',
            'nearest', 'gmt-triangulate', 'mbgrid'
        ]
        self.coastline_args = {}
        tmp_waffles = Waffle()
        tmp_coastline = WafflesCoastline()
        for kpam, kval in kwargs.items():
            if kpam not in tmp_waffles.__dict__:
                if kpam in tmp_coastline.__dict__:
                    self.coastline_args[kpam] = kval

        for kpam, kval in self.coastline_args.items():
            del kwargs[kpam]

        if exclude_lakes:
            self.coastline_args['want_lakes'] = True
            self.coastline_args['invert_lakes'] = True

        if mode is not None:
            pre_mode = mode
            
        self.pre_mode = pre_mode
        self.pre_mode_args = {}
        if self.pre_mode not in self.valid_modes:
            utils.echo_warning_msg(
                (f'{self.pre_mode} is not a valid waffles module! '
                 'falling back to `gmt-surface`')
            )
            self.pre_mode = 'gmt-surface'

        tmp_waffles_mode = WaffleFactory(mod=self.pre_mode)._acquire_module()
        for kpam, kval in kwargs.items():
            if kpam not in tmp_waffles.__dict__:
                if kpam in tmp_waffles_mode.__dict__:
                    self.pre_mode_args[kpam] = kval

        for kpam, kval in self.pre_mode_args.items():
            del kwargs[kpam]

        self.final_mode = final_mode
        self.final_mode_args = {}
        if self.final_mode not in self.valid_modes:
            self.final_mode = 'IDW'
            
        super().__init__(**kwargs)
        self.pre_count = utils.int_or(pre_count, 1)
        self.min_weight = min_weight        
        if min_weight is not None:
            weight_levels = min_weight
        
        if weight_levels is not None:
            self.weight_levels = [utils.float_or(w) for w in weight_levels.split('/')]
            self.weight_levels = self.weight_levels[:self.pre_count]
        else:
            self.weight_levels = []
            
        if inc_levels is not None:
            self.inc_levels = [utils.str2inc(i) for i in inc_levels.split('/')]
            self.inc_levels = self.inc_levels[:self.pre_count]
        else:
            self.inc_levels = []
            
        self.landmask = landmask
        self.invert_landmask = invert_landmask
        self.pre_upper_limit = utils.float_or(pre_upper_limit, -0.1) \
            if landmask \
               else None
        
        self.filter_outliers = filter_outliers
        self.want_supercede = want_supercede
        self.pre_smoothing = utils.float_or(pre_smoothing)
        if self.pre_smoothing is not None:
            self.pre_smoothing = ['blur:blur_factor={}'.format(self.pre_smoothing)]
            
        self.flatten = utils.float_or(flatten)
        self.want_weight = True
        self.pre_verbose = pre_verbose
        self.keep_pre_surfaces = keep_pre_surfaces

        
    ## todo: remove coastline after processing...
    def generate_coastline(self, pre_data=None):
        from .waffles import WaffleFactory
        
        cst_region = self.p_region.copy()
        cst_region.wmin = self.weight_levels[0]
        cst_fn = '{}_cst'.format(
            os.path.join(self.cache_dir, os.path.basename(self.name))
        )
        this_coastline = 'coastline:{}'.format(
            factory.dict2args(self.coastline_args)
        )

        ## update to call wafflesCoastline directly
        coastline = WaffleFactory(
            mod=this_coastline,
            data=pre_data,
            src_region=cst_region,
            want_weight=True,
            min_weight=self.weight_levels[0],
            xinc=self.xinc,
            yinc=self.yinc,
            name=cst_fn,
            node=self.node,
            dst_srs=self.dst_srs,
            srs_transform=self.srs_transform,
            clobber=True,
            cache_dir=self.cache_dir,
            verbose=self.pre_verbose)._acquire_module()
        coastline.initialize()
        coastline.generate()

        if coastline is not None:
            return('{}.shp:invert={}'.format(coastline.name, self.invert_landmask))
        else:
            return(None)

        
    def run(self):
        from .waffles import WaffleFactory
        
        ## set the weights if not already set correctly
        self.weight_levels.sort(reverse=True)
        if len(self.weight_levels) < self.pre_count+1:
            tmp_pre = self.pre_count - len(self.weight_levels)
            while tmp_pre >= 0:
                if len(self.weight_levels) == 0:
                    # self.stack isn't set yet!
                    tmp_weight = gdalfun.gdal_percentile(
                        self.stack, perc=99, band=3
                    ) 
                    tmp_weight = utils.float_or(tmp_weight, 1)
                else:
                    tmp_weight = self.weight_levels[-1]/(tmp_pre + 1)
                    if tmp_weight == 0:
                        tmp_weight = 1-e20

                if tmp_pre == 0:
                    tmp_weight = 0
                    
                self.weight_levels.append(tmp_weight)                
                tmp_pre -= 1

        ## set the increments if not already set correctly
        self.inc_levels.sort()
        if len(self.inc_levels) < self.pre_count+1:
            tmp_pre = self.pre_count - len(self.inc_levels)
            if len(self.inc_levels) == 0:
                tmp_xinc = self.xinc
                self.inc_levels.append(tmp_xinc)
                
            for t in range(1, tmp_pre+1):
                tmp_xinc = float(self.inc_levels[-1] * 3)
                tmp_yinc = float(self.inc_levels[-1] * 3)
                self.inc_levels.append(tmp_xinc)

        if self.inc_levels[0] != self.xinc or len(self.inc_levels) < self.pre_count+1:
            self.inc_levels.insert(0, self.xinc)
            self.inc_levels = self.inc_levels[:self.pre_count+1]
            self.inc_levels.sort()
        
        if self.verbose:
            init_str = ['\n==============================================\n']
            init_str.append(f'cudem generating {self.pre_count} pre-surface(s):\n')
            for i in range(len(self.inc_levels)-1, -1, -1):
                init_str.append((f'* {"pre_surface" if i !=0 else "final_surface"} {i} '
                                 f'<{self.pre_mode if i == self.pre_count else "stacks" if i != 0 else self.final_mode}> '
                                 f'@ {self.inc_levels[1]} using data with a min weight of {self.weight_levels[i]}'))

            init_str.append(f'\ncudem using stack file: {self.stack}')
            if self.landmask:
                if isinstance(self.landmask, str):
                    if os.path.exists(self.landmask.split(':')[0]):
                        init_str.append(f'cudem using coastline: {self.landmask}')
                    else:
                        self.landmask = True
                        
                if isinstance(self.landmask, bool):
                    init_str.append('cudem using coastline: Waffles COASTLINE module')

            init_str.append(f'cudem flattening: {self.flatten}')
            init_str.append(f'cudem output DEM: {self.name}')
            init_str.append('\n==============================================')
            
            out_str = '\n'.join(init_str)
            utils.echo_msg(out_str)
            
        orig_stack = self.stack
        pre = self.pre_count
        pre_weight = 0 # initial run will use all data weights
        pre_region = self.p_region.copy()
        pre_region.wmin = None

        # _pre_name_minus = os.path.join(
        #     self.cache_dir, utils.append_fn('_pre_surface', pre_region, pre)
        # )
        # ## initial data to pass through surface (stack)
        stack_data_entry = (f'"{self.stack}",200:band_no=1:weight_mask=3:'
                            'uncertainty_mask=4:sample=average,1')
        pre_data = [stack_data_entry]
         
        ## generate coastline
        pre_clip = None
        if self.landmask:            
            if isinstance(self.landmask, str):
                if os.path.exists(self.landmask.split(':')[0]):
                    # todo: update to make 'invert' an option
                    pre_clip = '{}:invert={}'.format(self.landmask, self.invert_landmask) 

            if pre_clip is None:
                coast_data = [
                    (f'"{self.stack}",200:band_no=1:weight_mask=3:'
                     'uncertainty_mask=4:sample=cubicspline,1')]
                coastline = self.generate_coastline(pre_data=coast_data)
                pre_clip = coastline

        ## Grid/Stack the data `pre` times concluding in full
        ## resolution with data > min_weight
        pre_surfaces = []
        with utils.ccp(
                total=self.pre_count+1,
                desc='generating CUDEM surfaces',
                leave=self.verbose
        ) as pbar:
            for pre in range(self.pre_count, -1, -1): # pre_count - 1
                pre_xinc = self.inc_levels[pre]
                pre_yinc = self.inc_levels[pre]
                xsample = self.inc_levels[pre-1] if pre != 0 else self.xinc
                ysample = self.inc_levels[pre-1] if pre != 0 else self.yinc

                ## initial data to pass through surface (stack)

                
                ## if not final or initial output, setup the configuration
                ## for the pre-surface
                if pre != self.pre_count:
                    pre_weight = self.weight_levels[pre]
                    _pre_name_plus = os.path.join(
                        self.cache_dir, utils.append_fn('_pre_surface', pre_region, pre+1)
                    )
                    _pre_unc_name = f'{_pre_name_plus}_u.tif' if self.want_uncertainty else None
                    pre_data_entry = (f'"{_pre_name_plus}.tif",200'
                                      f':uncertainty_mask={"_pre_unc_name" if _pre_unc_name else None}'
                                      f':sample=cubicspline:check_path=True'
                                      f',{pre_weight-.1}')

                    pre_data = [stack_data_entry, pre_data_entry]
                    pre_region.wmin = None#pre_weight

                ## reset pre_region for final grid
                if pre == 0:
                    pre_region = self.p_region.copy()
                    pre_region.wmin = None#self.weight_levels[pre]

                _pre_name = os.path.join(
                    self.cache_dir, utils.append_fn('_pre_surface', pre_region, pre)
                )

                last_fltr = [
                    (f'weights:stacks=True:weight_threshold={pre_weight}'
                     ':buffer_cells=1:verbose=False')
                     ]
                waffles_mod = '{}:{}'.format(
                    self.pre_mode,
                    factory.dict2args(self.pre_mode_args)
                ) if pre==self.pre_count \
                else 'stacks'\
                     if pre != 0 \
                        else 'IDW'

                utils.echo_msg(
                    'cudem gridding surface {} @ {} {}/{} to {} using {}...'.format(
                        pre, pre_region.format('str'), pre_xinc, pre_yinc, _pre_name, waffles_mod
                    )
                )
                pre_surface = WaffleFactory(
                    mod=waffles_mod,
                    data=pre_data,
                    src_region=pre_region,
                    xinc=pre_xinc,
                    yinc=pre_yinc,
                    xsample=xsample,
                    ysample=ysample,
                    name=_pre_name,
                    node='pixel',
                    want_weight=True,
                    want_uncertainty=self.want_uncertainty,
                    dst_srs=self.dst_srs,
                    srs_transform=self.srs_transform,
                    clobber=True,
                    verbose=self.pre_verbose,
                    clip=pre_clip if pre !=0 else None,
                    stack_mode='mixed:weight_threshold={}'.format(pre_weight),
                    upper_limit=self.pre_upper_limit if pre != 0 else None,
                    keep_auxiliary=False,
                    fltr=self.pre_smoothing if pre != 0 else None,
                    percentile_limit=self.flatten if pre == 0 else None,
                    cache_dir=self.cache_dir,
                )._acquire_module()
                pre_surface.initialize()
                pre_surface.generate()
                pbar.update()
                pre_surfaces.append(pre_surface.name)

        ## todo add option to flatten here...or move flatten up
        gdalfun.cudem_flatten_no_data_zones(
            pre_surface.fn, dst_dem=self.fn, band=1, size_threshold=1
        )

        ## reset the stack for uncertainty
        self.stack = orig_stack
        if not self.keep_pre_surfaces:
            utils.remove_glob(*[f"{x}*" for x in pre_surfaces])
        
        return self


### End
