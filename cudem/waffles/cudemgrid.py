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
## CUDEM Integrated Gridding Module.
##
## This module implements a multi-resolution "step-down" gridding strategy.
## It generates a DEM by iterating through lower resolutions (pre-surfaces) 
## to fill gaps, progressively increasing resolution and weight thresholds 
## until the final target resolution is reached.
##
## 1. Generates a low-res background surface using all data (low weights allowed).
## 2. Generates intermediate surfaces at increasing resolutions, filtering out low-weight data.
## 3. Uses the previous surface as a background fill for the next.
## 4. Generates the final surface at target resolution using only high-quality data.
##
### Code:

import os
import numpy as np
from osgeo import ogr

from cudem import utils
from cudem import gdalfun
from cudem import factory
from cudem.grits import grits
from cudem.waffles.waffles import Waffle
from cudem.waffles.coastline import WafflesCoastline

class WafflesCUDEM(Waffle):
    """CUDEM Integrated DEM Generation.
    
    Parameters:
    -----------    
    pre_mode (str): Waffles module for initial pre-surface (default: 'gmt-surface').
    pre_count (int): Number of iterations/pre-surfaces to generate (default: 2).
    weight_levels (str): Slash-sep list of min weights for iterations (e.g. "1/0.5").
    inc_levels (str): Slash-sep list of increments for iterations (e.g. "1s/3s").
    landmask (bool/str): Use/Path to coastline mask.
    invert_landmask (bool): Invert coastline mask (Default True = Mask Ocean).
    pre_upper_limit (float): Max Z for pre-surfaces (used with landmask).
    pre_smoothing (float): Smoothing factor (buffer cells) for pre-surfaces. 
                           Uses 'grits:blend' to blend high/low weight data.
    flatten (float): Percentile threshold to flatten NoData zones in final output.
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

        ## Separate arguments for sub-modules
        ## We pass Waffle explicitly to _extract so it knows what to protect
        self.coastline_args = self._extract_submodule_args(kwargs, WafflesCoastline)
        
        if exclude_lakes:
            self.coastline_args['want_lakes'] = True
            self.coastline_args['invert_lakes'] = True

        ## Validate Modes
        self.valid_modes = [
            'gmt-surface', 'IDW', 'linear', 'cubic',
            'nearest', 'gmt-triangulate', 'mbgrid'
        ]
        
        self.pre_mode = mode if mode else pre_mode
        if self.pre_mode not in self.valid_modes:
            utils.echo_warning_msg(f'{self.pre_mode} invalid, falling back to `gmt-surface`')
            self.pre_mode = 'gmt-surface'
            
        self.final_mode = final_mode if final_mode in self.valid_modes else 'IDW'

        ## Extract args for pre-mode
        from .waffles import WaffleFactory
        tmp_pre = WaffleFactory(mod=self.pre_mode)._acquire_module()
        
        ## We pass type(tmp_pre) to extract specific args for that module
        self.pre_mode_args = self._extract_submodule_args(kwargs, type(tmp_pre))
        
        ## 4. Initialize Base Waffle (Now kwargs is safe to pass)
        super().__init__(**kwargs)

        ## 5. Configuration (Rest of init...)
        self.pre_count = utils.int_or(pre_count, 1)
        self.weight_config = min_weight if min_weight else weight_levels
        self.inc_config = inc_levels
        self.weight_levels = []
        self.inc_levels = []
        self.landmask = landmask
        self.invert_landmask = invert_landmask
        self.pre_upper_limit = utils.float_or(pre_upper_limit, -0.1) if landmask else None
        self.filter_outliers = filter_outliers
        self.want_supercede = want_supercede
        self.pre_smoothing = utils.float_or(pre_smoothing)
        self.flatten = utils.float_or(flatten)
        self.want_weight = True
        self.pre_verbose = pre_verbose
        self.keep_pre_surfaces = keep_pre_surfaces

        
    def _extract_submodule_args(self, kwargs, module_class):
        """Extract arguments relevant to a specific sub-module class,
        protecting base Waffle arguments from being stolen."""
        
        extracted = {}
        
        ## Instantiate dummies to check keys
        target_keys = set(module_class().__dict__.keys())
        
        ## We must protect keys that belong to the base Waffle class
        ## to ensure WafflesCUDEM (self) initializes correctly.
        base_keys = set(Waffle().__dict__.keys())
        ## Explicitly protect 'params' and 'mod' which are critical
        base_keys.update(['params', 'mod'])

        ## Find intersection in kwargs
        for k in list(kwargs.keys()):
            ## Take it if it belongs to target AND is NOT a base key
            if k in target_keys and k not in base_keys:
                extracted[k] = kwargs.pop(k)
                
        return extracted
    
    
    def _setup_levels(self):
        """Calculate the resolution and weight steps for iterations."""
        
        ## Weight Levels
        if self.weight_config:
            self.weight_levels = [utils.float_or(w) for w in self.weight_config.split('/')]
            self.weight_levels = self.weight_levels[:self.pre_count]
        
        ## Auto-fill missing weights
        ## We assume descending weights (Highest quality last)
        self.weight_levels.sort(reverse=True)
        
        while len(self.weight_levels) <= self.pre_count:
            if not self.weight_levels:
                ## If no weights provided, try to guess max from stack (99th percentile)
                ## or default to 1.
                try:
                    tmp_weight = gdalfun.gdal_percentile(self.stack, perc=99, band=3)
                except:
                    tmp_weight = 1
                self.weight_levels.append(utils.float_or(tmp_weight, 1))
            else:
                ## Halve the previous weight for the next lower level
                next_weight = self.weight_levels[-1] / 2
                if next_weight == 0: next_weight = 1e-20 # Avoid absolute zero
                self.weight_levels.append(next_weight)
                
        ## Ensure final level (0) is 0 to include everything in first pass
        if len(self.weight_levels) > self.pre_count:
             self.weight_levels[-1] = 0

        ## Increment Levels
        if self.inc_config:
            self.inc_levels = [utils.str2inc(i) for i in self.inc_config.split('/')]
            self.inc_levels = self.inc_levels[:self.pre_count]
            
        ## Auto-fill missing increments
        ## We assume ascending increments (Low res -> High res)
        self.inc_levels.sort()
        
        ## Make sure target resolution is at index 0
        if not self.inc_levels or self.inc_levels[0] != self.xinc:
            self.inc_levels.insert(0, self.xinc)
            
        while len(self.inc_levels) <= self.pre_count:
            ## Triple the resolution for each step down (1s -> 3s -> 9s)
            self.inc_levels.append(self.inc_levels[-1] * 3)
            
        ## Sort so index 0 is target (smallest) and index -1 is coarsest
        self.inc_levels.sort()

        
    def generate_coastline(self, pre_data=None):
        """Generate and robustly dissolve the coastline mask."""
        
        from .waffles import WaffleFactory
        
        cst_region = self.p_region.copy()
        cst_region.wmin = self.weight_levels[-1] # Use lowest weight for coast
        
        cst_name = os.path.join(self.cache_dir, f"{os.path.basename(self.name)}_cst")
        
        ## Generate Raw Coastline (Chunks)
        coast_conf = f'coastline:{factory.dict2args(self.coastline_args)}'

        coast_waffle = WaffleFactory(
            mod=coast_conf,
            data=pre_data,
            src_region=cst_region,
            want_weight=True,
            min_weight=self.weight_levels[-1],
            xinc=self.xinc,
            yinc=self.yinc,
            name=cst_name,
            node=self.node,
            dst_srs=self.dst_srs,
            srs_transform=self.srs_transform,
            clobber=True,
            cache_dir=self.cache_dir,
            verbose=self.pre_verbose
        )._acquire_module()
        
        coast_waffle.initialize()
        coast_waffle.generate()

        if coast_waffle is None:
            return None

        ## Dissolve and Buffer to fix seams
        cst_shp = f'{coast_waffle.name}.shp'
        if os.path.exists(cst_shp):
            try:
                self._dissolve_coastline(cst_shp)
            except Exception as e:
                utils.echo_warning_msg(f"Failed to dissolve coastline: {e}")

        return f'{cst_shp}:invert={self.invert_landmask}'

    
    def _dissolve_coastline(self, shp_path):
        """Dissolve chunked coastline polygons and buffer to remove artifacts."""
        
        ds = ogr.Open(shp_path, 0) # Read-only
        if not ds: return

        layer = ds.GetLayer()
        srs = layer.GetSpatialRef()
        
        ## Collect geometries
        geom_col = ogr.Geometry(ogr.wkbGeometryCollection)
        for feature in layer:
            ref = feature.GetGeometryRef()
            if ref: geom_col.AddGeometry(ref.Clone())
        
        ds = None # Close file

        if geom_col.IsEmpty(): return

        ## Union (Dissolve)
        dissolved = geom_col.UnionCascaded()
        
        ## Buffer (Fix micro-gaps)
        ## Use 10% of pixel size or default small value
        buff_dist = abs(float(self.xinc)) * 0.1 if self.xinc else 1e-5
        clean_geom = dissolved.Buffer(buff_dist)

        if not clean_geom or clean_geom.IsEmpty(): return

        ## Rewrite Shapefile
        driver = ogr.GetDriverByName('ESRI Shapefile')
        if os.path.exists(shp_path):
            driver.DeleteDataSource(shp_path)
            
        out_ds = driver.CreateDataSource(shp_path)
        out_layer = out_ds.CreateLayer('coastline', srs=srs, geom_type=clean_geom.GetGeometryType())
        
        feat = ogr.Feature(out_layer.GetLayerDefn())
        feat.SetGeometry(clean_geom)
        out_layer.CreateFeature(feat)
        
        feat = None
        out_ds = None
        
        if self.pre_verbose:
            utils.echo_msg(f"Dissolved and buffered coastline: {shp_path}")

            
    def run(self):
        from .waffles import WaffleFactory
        
        ## Setup Resolution/Weight Ladders
        self._setup_levels()

        if self.verbose:
            self._print_plan()

        orig_stack = self.stack

        if self.pre_smoothing:
            for w in self.weight_levels:
                grits_filter = grits.GritsFactory(
                    mod=f'blend:weight_threshold{w}:buffer_cells={int(self.pre_smoothing)*4}:sub_buffer_cells={int(self.pre_smoothing)}:random_scale=.5',
                    src_dem=orig_stack,
                    uncertainty_mask=4,
                    weight_mask=3,
                    count_mask=2,
                    cache_dir=self.cache_dir,
                    verbose=True
                )._acquire_module()

                if grits_filter:
                    grits_filter = grits_filter()
                    if os.path.exists(grits_filter.dst_dem):
                        orig_stack = grits_filter.dst_dem
                    else:
                        utils.echo_warning_msg('Grits output in invalid: {grits_filter.dst_dem}')

        
        ## Initial Data Setup
        ## Define the base data entry (the stack)
        stack_entry = (f'"{self.stack}",200:band_no=1:weight_mask=3:'
                       'uncertainty_mask=4:sample=average,1')
        
        ##  Generate Coastline (if needed)
        pre_clip = None
        if self.landmask:
            if isinstance(self.landmask, str) and os.path.exists(self.landmask.split(':')[0]):
                ## User provided file
                pre_clip = f'{self.landmask}:invert={self.invert_landmask}'
            else:
                # Auto-generate
                pre_clip = self.generate_coastline(pre_data=[stack_entry])

        ## Iterative Gridding Loop
        ## We iterate backwards from pre_count down to 0
        ## pre_count (Highest index) = Coarsest Resolution / Lowest Weight
        ## 0 (Lowest index) = Final Target Resolution / Highest Weight        
        pre_surfaces_generated = []
        
        with utils.ccp(total=self.pre_count+1, desc='Generating CUDEM', leave=self.verbose) as pbar:
            for i in range(self.pre_count, -1, -1):
                
                ## Configuration for this level
                current_inc = self.inc_levels[i]
                current_weight = self.weight_levels[i]
                
                ## Determine Resampling Resolution (Sample at next level up, or target if final)
                ## i-1 is the next finer resolution. If i=0, we are at target.
                sample_idx = max(0, i-1)
                xsample = self.inc_levels[sample_idx]
                ysample = self.inc_levels[sample_idx]

                current_region = self.p_region.copy()

                if i > 0:
                    # Buffer by ~10 cells of the current resolution to be safe
                    # (e.g. 10 * 3s = 30s buffer)
                    x_buff = float(current_inc) * 10
                    y_buff = float(current_inc) * 10
                    current_region.buffer(x_bv=x_buff, y_bv=y_buff)
                
                ## Input Data Setup
                ## Start with the base stack
                current_data = [stack_entry]
                
                ## If not the very first (coarsest) run, add the previous surface as background
                if i < self.pre_count:
                    prev_surface = pre_surfaces_generated[-1]
                    
                    ## Add previous surface to input list with a slightly lower weight
                    ## so real data supercedes it.
                    ## We subtract 0.1 from current threshold to ensure it acts as fill.
                    bg_weight = max(0, current_weight - 0.1)
                    
                    bg_entry = (f'"{prev_surface}",200:sample=cubicspline:check_path=True,'
                                f'{bg_weight}')
                    current_data.append(bg_entry)

                ## Output Naming
                ## i=0 is final pre-surface before post-process
                suffix = f'_pre_{i}' if i > 0 else '_final_raw'
                out_name = os.path.join(
                    self.cache_dir, utils.append_fn(suffix, current_region, i)
                )

                ## Select Module
                ## Iterations use 'pre_mode' (e.g. gmt-surface) to fill gaps.
                ## Final uses 'IDW' or 'stacks' to enforce hard data, unless override.
                if i == self.pre_count:
                    ## Initial coarse fill
                    mod_str = f'{self.pre_mode}:{factory.dict2args(self.pre_mode_args)}'
                elif i == 0:
                    ## Final run
                    mod_str = self.final_mode
                else:
                    ## Intermediate runs (usually just stacking/blending)
                    mod_str = 'stacks'

                ## Use 'blend' filter for intermediate steps (excluding first and last)
                ## blend blends the current high-weight data into the previous low-weight background
                this_fltr = None
                if self.pre_smoothing and i < self.pre_count:# and i > 0:
                    this_fltr = [
                        f'blend:stacks=True:weight_threshold={current_weight}:sub_buffer_cells={int(self.pre_smoothing)}:buffer_cells={int(self.pre_smoothing)*4}:verbose=True'
                    ]
                    
                if self.verbose:
                    utils.echo_msg(f'Step {i}: {mod_str} @ {current_inc} (Min Weight: {current_weight})')
                    
                ## Run Waffles
                waffle_step = WaffleFactory(
                    mod=mod_str,
                    data=current_data,
                    src_region=current_region,
                    xinc=current_inc,
                    yinc=current_inc,
                    xsample=xsample,
                    ysample=ysample,
                    name=out_name,
                    node='pixel',
                    want_weight=True,
                    want_uncertainty=self.want_uncertainty,
                    dst_srs=self.dst_srs,
                    srs_transform=self.srs_transform,
                    clobber=True,
                    verbose=self.pre_verbose,
                    ## Only clip intermediate surfaces, usually leave final unclipped to fill edges
                    clip=pre_clip if i > 0 else None, 
                    ## Set stacking mode to exclude low-weight data
                    stack_mode=f'mixed:weight_threshold={current_weight}',
                    ## Upper limit only applies to initial coarse fills (e.g. remove land noise)
                    upper_limit=self.pre_upper_limit if i > 0 else None,
                    keep_auxiliary=False,
                    ## Apply smoothing to intermediate steps
                    #fltr=self.pre_smoothing if i > 0 else None,
                    fltr=this_fltr,
                    ## Apply flattening only on final step
                    percentile_limit=self.flatten if i == 0 else None,
                    cache_dir=self.cache_dir,
                )._acquire_module()
                
                waffle_step.initialize()
                waffle_step.generate()
                
                pre_surfaces_generated.append(waffle_step.fn)
                pbar.update()

        ## Finalize
        final_raw = pre_surfaces_generated[-1]
        
        if self.flatten:
            ## If flattening was requested as a post-process
            gdalfun.cudem_flatten_no_data_zones(
                final_raw, dst_dem=self.fn, band=1, size_threshold=1
            )
        else:
            ## Just move the last result to the requested output filename
            if os.path.exists(self.fn):
                try: os.remove(self.fn)
                except: pass
            os.rename(final_raw, self.fn)

        ## Cleanup
        self.stack = orig_stack # Restore stack pointer
        if not self.keep_pre_surfaces:
            ## Remove all intermediate files except the one we just renamed
            to_remove = [f"{x}*" for x in pre_surfaces_generated if x != final_raw]
            utils.remove_glob(*to_remove)
        
        return self

    
    def _print_plan(self):
        """Debug print of the generation plan."""
        
        utils.echo_msg('==============================================')
        utils.echo_msg(f'CUDEM Generation Plan ({self.pre_count + 1} steps)')
        utils.echo_msg(f'Output: {self.name}')
        utils.echo_msg(f'Landmask: {self.landmask}')
        
        for i in range(self.pre_count, -1, -1):
            w = self.weight_levels[i]
            inc = self.inc_levels[i]
            label = "Final" if i == 0 else "Pre-Surface"
            utils.echo_msg(f'  * Step {i} [{label}]: Res={inc}, MinWeight={w}')
        
        utils.echo_msg('==============================================')


### End
