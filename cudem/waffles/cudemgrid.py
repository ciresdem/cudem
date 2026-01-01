### cudemgrid.py
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
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
## CUDEM integrated DEM generation.
##
### Code:

import os
from cudem import utils
from cudem import gdalfun
from cudem import factory
from cudem.waffles.waffles import Waffle
from cudem.waffles.coastline import WafflesCoastline

class WafflesCUDEM(Waffle):
    """CUDEM integrated DEM generation.
    
    Generate a topo/bathy integrated DEM using a variety of data sources.
    Iterates <pre_count> pre-surfaces at lower-resolutions, specified with <inc_levels>.
    
    Each pre-surface is used in subsequent pre-surface(s)/final DEM at each 
    iterative weight specified in <weight_levels>.

    Parameters:
    -----------
    pre_mode (str) : The waffles module to perform the initial pre-surface
    pre_count (int) : Number of pre-surface iterations to perform
    weight_levels (str) : Slash-delimited list of min weights for each pre-surface (e.g. "1/0.5")
    inc_levels (str) : Slash-delimited list of increments for each pre-surface (e.g. "0.3s/1s")
    pre_upper_limit (float) : The upper elevation limit of pre-surfaces (used with landmask)
    pre_smoothing (float) : The smoothing (blur) factor to apply to the pre-surfaces
    pre_verbose (bool) : Increase the verbosity of pre-surface generation
    landmask (bool/str) : Path to coastline vector mask or True to auto-generate
    invert_landmask (bool) : Invert the coastline mask (True = mask ocean)
    want_supercede (bool) : Supercede subsequent pre-surfaces
    flatten (float) : The nodata-size percentile above which to flatten
    exclude_lakes (bool) : Exclude lakes from the landmask
    filter_outliers (float) : Percentile at which to filter outliers
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
        from cudem.waffles.waffles import WaffleFactory
        
        self.valid_modes = [
            'gmt-surface', 'IDW', 'linear', 'cubic',
            'nearest', 'gmt-triangulate', 'mbgrid'
        ]
        
        ## Extract Coastline Args
        self.coastline_args = {}
        tmp_waffles = Waffle()
        tmp_coastline = WafflesCoastline()
        
        ## Move relevant kwargs to coastline_args
        for kpam, kval in list(kwargs.items()):
            if kpam not in tmp_waffles.__dict__ and kpam in tmp_coastline.__dict__:
                self.coastline_args[kpam] = kval
                del kwargs[kpam]

        if exclude_lakes:
            self.coastline_args['want_lakes'] = True
            self.coastline_args['invert_lakes'] = True

        ## Extract Pre-Mode Args
        if mode is not None:
            pre_mode = mode
            
        self.pre_mode = pre_mode
        if self.pre_mode not in self.valid_modes:
            utils.echo_warning_msg(
                f'pre_mode: {self.pre_mode} is not a valid waffles module! falling back to `gmt-surface`'
            )
            self.pre_mode = 'gmt-surface'

        self.pre_mode_args = {}
        tmp_waffles_mode = WaffleFactory(mod=self.pre_mode)._acquire_module()
        
        ## Move relevant kwargs to pre_mode_args
        for kpam, kval in list(kwargs.items()):
            if kpam not in tmp_waffles.__dict__ and kpam in tmp_waffles_mode.__dict__:
                self.pre_mode_args[kpam] = kval
                del kwargs[kpam]

        self.final_mode = final_mode
        if self.final_mode not in self.valid_modes:
            self.final_mode = 'IDW'
            
        super().__init__(**kwargs)
        
        ## Setup CUDEM specific attributes
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
        self.pre_upper_limit = utils.float_or(pre_upper_limit, -0.1) if landmask else None
        
        self.filter_outliers = filter_outliers
        self.want_supercede = want_supercede
        self.pre_smoothing = utils.float_or(pre_smoothing)
        if self.pre_smoothing is not None:
            self.pre_smoothing = [f'blur:blur_factor={self.pre_smoothing}']
            
        self.flatten = utils.float_or(flatten)
        self.want_weight = True
        self.pre_verbose = pre_verbose
        self.keep_pre_surfaces = keep_pre_surfaces

        
    def generate_coastline(self, pre_data=None):
        """Generate a coastline mask for the region."""

        from cudem.waffles.waffles import WaffleFactory
        
        cst_region = self.p_region.copy()
        if self.weight_levels:
            cst_region.wmin = self.weight_levels[0]
            
        utils.echo_msg(f'Coast region is: {cst_region}')
        cst_fn = f"{os.path.join(self.cache_dir, os.path.basename(self.name))}_cst"
        
        this_coastline = f"coastline:{factory.dict2args(self.coastline_args)}"
        coastline = WaffleFactory(
            mod=this_coastline,
            data=pre_data,
            src_region=cst_region,
            want_weight=True,
            min_weight=self.weight_levels[0] if self.weight_levels else 0,
            xinc=self.xinc,
            yinc=self.yinc,
            name=cst_fn,
            node=self.node,
            dst_srs=self.dst_srs,
            srs_transform=self.srs_transform,
            clobber=True,
            cache_dir=self.cache_dir,
            verbose=True
        )._acquire_module()
        
        coastline.initialize()
        coastline.generate()

        if coastline is not None:
            return f"{coastline.name}.shp:invert={self.invert_landmask}"
        else:
            return None

        
    def run(self):
        """Execute the CUDEM generation process."""

        from cudem.waffles.waffles import WaffleFactory
        
        ## Configure Weight Levels
        self.weight_levels.sort(reverse=True)
        if len(self.weight_levels) < self.pre_count + 1:
            tmp_pre = self.pre_count - len(self.weight_levels)
            while tmp_pre >= 0:
                if len(self.weight_levels) == 0:
                    ## Default if stack isn't set yet
                    tmp_weight = gdalfun.gdal_percentile(self.stack, perc=99, band=3) 
                    tmp_weight = utils.float_or(tmp_weight, 1)
                else:
                    tmp_weight = self.weight_levels[-1] / (tmp_pre + 1)
                    if tmp_weight == 0:
                        tmp_weight = 1e-20

                if tmp_pre == 0:
                    tmp_weight = 0
                    
                self.weight_levels.append(tmp_weight)                
                tmp_pre -= 1

        ## Configure Increment Levels
        self.inc_levels.sort()
        if len(self.inc_levels) < self.pre_count + 1:
            tmp_pre = self.pre_count - len(self.inc_levels)
            if len(self.inc_levels) == 0:
                self.inc_levels.append(self.xinc)
                
            for t in range(1, tmp_pre + 1):
                tmp_xinc = float(self.inc_levels[-1] * 3)
                self.inc_levels.append(tmp_xinc)

        if self.inc_levels[0] != self.xinc or len(self.inc_levels) < self.pre_count + 1:
            self.inc_levels.insert(0, self.xinc)
            self.inc_levels = self.inc_levels[:self.pre_count + 1]
            self.inc_levels.sort()
        
        if self.verbose:
            utils.echo_msg_bold('==============================================')
            utils.echo_msg_bold(f'CUDEM generating {self.pre_count} pre-surface(s):')
            for i in range(len(self.inc_levels) - 1, -1, -1):
                mode_str = self.pre_mode if i == self.pre_count \
                           else 'stacks' if i != 0 \
                           else self.final_mode
                
                utils.echo_msg_bold(
                    f"* {'pre_surface' if i != 0 else 'final surface'} {i} "
                    f"<{mode_str}> @ {self.inc_levels[i]} "
                    f"using data with a minimum weight of {self.weight_levels[i]}"
                )
                
            utils.echo_msg_bold(f'CUDEM using stack file: {self.stack}')
            if self.landmask:
                if isinstance(self.landmask, str):
                    if os.path.exists(self.landmask.split(':')[0]):
                        utils.echo_msg_bold(
                            f'CUDEM using coastline: {self.landmask}'
                        )
                    else:
                        self.landmask = True
                        
                if isinstance(self.landmask, bool):
                    utils.echo_msg_bold(
                        'CUDEM using coastline: Waffles COASTLINE module'
                    )
                
            utils.echo_msg_bold(f'CUDEM flattening: {self.flatten}')
            utils.echo_msg_bold(f'CUDEM output DEM: {self.name}')
            utils.echo_msg_bold('==============================================')

        orig_stack = self.stack
        pre_region = self.p_region.copy()
        pre_region.wmin = None
        
        ## Initial data
        stack_data_entry = (f'"{self.stack}",200:band_no=1:weight_mask=3:'
                            'uncertainty_mask=4:sample=average,1')
        pre_data = [stack_data_entry]
         
        ## Generate Coastline / Clip
        pre_clip = None
        if self.landmask:            
            if isinstance(self.landmask, str) and os.path.exists(self.landmask.split(':')[0]):
                pre_clip = f'{self.landmask}:invert={self.invert_landmask}'

            if pre_clip is None:
                coast_data = [(f'"{self.stack}",200:band_no=1:weight_mask=3:'
                               'uncertainty_mask=4:sample=cubicspline,1')]
                coastline = self.generate_coastline(pre_data=coast_data)
                pre_clip = coastline

        ## Grid/Stack the data `pre` times
        pre_surfaces = []
        pre_surface = None
        
        with utils.ccp(total=self.pre_count + 1, desc='generating CUDEM surfaces', leave=self.verbose) as pbar:
            for pre in range(self.pre_count, -1, -1): 
                pre_xinc = self.inc_levels[pre]
                pre_yinc = self.inc_levels[pre]
                xsample = self.inc_levels[pre - 1] if pre != 0 else self.xinc
                ysample = self.inc_levels[pre - 1] if pre != 0 else self.yinc

                ## Configure inputs for intermediate pre-surfaces
                if pre != self.pre_count:
                    pre_weight = self.weight_levels[pre]
                    _pre_name_plus = os.path.join(
                        self.cache_dir, utils.append_fn('_pre_surface', pre_region, pre + 1)
                    )
                    _pre_unc_name = f'{_pre_name_plus}_u.tif' if self.want_uncertainty else None
                    
                    pre_data_entry = (f'"{_pre_name_plus}.tif",200'
                                      f':uncertainty_mask="{_pre_unc_name}"'
                                      f':sample=cubicspline:check_path=True'
                                      f',{pre_weight - .1}')
                    
                    pre_data = [stack_data_entry, pre_data_entry]
                    pre_region.wmin = None

                ## Reset region for final grid
                if pre == 0:
                    pre_region = self.p_region.copy()
                    pre_region.wmin = None

                _pre_name = os.path.join(
                    self.cache_dir, utils.append_fn('_pre_surface', pre_region, pre)
                )

                ## Determine Waffles Module
                if pre == self.pre_count:
                    waffles_mod = f"{self.pre_mode}:{factory.dict2args(self.pre_mode_args)}"
                elif pre != 0:
                    waffles_mod = 'stacks'
                else:
                    waffles_mod = 'IDW'

                utils.echo_msg(f"CUDEM gridding surface {pre} to {_pre_name}...")
                
                ## Stack Mode
                stack_mode = f'mixed:weight_threshold={pre_weight}' if pre != self.pre_count else None

                ## Generate the pre-surface
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
                    clip=pre_clip if pre != 0 else None,
                    stack_mode=stack_mode,
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

        ## Flatten NoData Zones in final output
        if pre_surface and os.path.exists(pre_surface.fn):
             gdalfun.cudem_flatten_no_data_zones(
                pre_surface.fn, dst_dem=self.fn, band=1, size_threshold=1
            )

        ## Cleanup
        self.stack = orig_stack
        if not self.keep_pre_surfaces:
            utils.remove_glob(*[f"{x}*" for x in pre_surfaces])
        
        return self


### End
