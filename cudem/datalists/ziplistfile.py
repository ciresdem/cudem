### ziplistfile.py 
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
##
## ziplistfile.py is part of CUDEM
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
## Zip file dataset parser.
##
### Code:

import os
import zipfile
import numpy as np
        
from cudem import utils
from cudem import regions
from cudem import pointz
from cudem.datalists.dlim import ElevationDataset

class ZIPlist(ElevationDataset):
    """
    Zip file parser.
    Parse supported datasets from a zipfile.
    """
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)


    def generate_inf(self, make_grid=True, make_block_mean=False, block_inc=None):
        """Generate metadata (INF) for the zip container.
        
        Iterates through contained datasets to calculate aggregated bounds and stats.
        Optionally generates a merged mini-grid and block-mean from the children.
        """
        
        ## Setup Accumulators
        _region_backup = self.region
        self.region = None # Unset region to scan everything
        
        self.infos.file_hash = self.infos.generate_hash()
        
        point_count = 0
        out_regions = []
        out_srs = []
        
        ## Iterate Children (Pass 1: Gather Extents & Stats)
        valid_children = []
        for entry in self.parse():
            ## Trigger child INF generation if missing (recursive scan of child)
            ## We disable grid generation for children to save time during this pass
            entry.inf(make_grid=False, make_block_mean=False) 
            
            if entry.infos.minmax:
                valid_children.append(entry)
                entry_region = regions.Region().from_list(entry.infos.minmax)
                
                ## Handle SRS transforms for bounds aggregation
                if entry.src_srs is not None:
                    out_srs.append(entry.src_srs)
                    if self.dst_srs is not None:
                        entry_region.src_srs = entry.src_srs
                        entry_region.warp(self.dst_srs)

                if entry_region.valid_p():
                    out_regions.append(entry_region)
                    point_count += entry.infos.numpts

        ## Aggregate Regions
        master_region = None
        for r in out_regions:
            if master_region is None: master_region = r
            else: master_region = regions.regions_merge(master_region, r)
            
        if master_region is not None:
            self.infos.minmax = master_region.export_as_list(include_z=True)
            self.infos.wkt = master_region.export_as_wkt()
            self.infos.numpts = point_count
            
            ## Determine Master SRS
            if self.infos.src_srs is None:
                if self.src_srs is not None:
                    self.infos.src_srs = self.src_srs
                elif out_srs:
                    ## If all children have same SRS, use it
                    if all(x == out_srs[0] for x in out_srs):
                        self.infos.src_srs = out_srs[0]
                    self.src_srs = self.infos.src_srs

            ## Generate Grids (Pass 2 - Aggregation)
            if (make_grid or make_block_mean) and point_count > 0:
                self._generate_grids_from_children(
                    valid_children, master_region, 
                    make_grid, make_block_mean, block_inc
                )

        self.region = _region_backup
        return self.infos

    
    def _generate_grids_from_children(self, children, region, make_grid, make_block_mean, block_inc):
        """Internal helper to generate grids by iterating over children.
        """
        
        ## Setup Gridders
        pp_grid = None
        pp_block = None
        grid_arrays = None
        block_arrays = None
        
        if make_grid:
            pp_grid = pointz.PointPixels(src_region=region, x_size=10, y_size=10, verbose=False)
            grid_arrays = {'sum': np.zeros((10, 10)), 'count': np.zeros((10, 10))}

        if make_block_mean:
            if block_inc is None:
                width = region.xmax - region.xmin
                block_inc = width / 500.0 if width > 0 else 0.001
            
            try:
                bx, by, _ = region.geo_transform(x_inc=block_inc, y_inc=block_inc)
                pp_block = pointz.PointPixels(src_region=region, x_size=bx, y_size=by, verbose=False)
                block_arrays = {
                    'z_sum': np.zeros((by, bx)), 'count': np.zeros((by, bx)),
                    'x_sum': np.zeros((by, bx)), 'y_sum': np.zeros((by, bx))
                }
            except Exception:
                pp_block = None

        ## Helper to accumulate chunk results into master grid
        def accumulate(master_dict, chunk_res, chunk_srcwin, keys):
            x_off, y_off, x_s, y_s = chunk_srcwin
            y_slice = slice(y_off, y_off + y_s)
            x_slice = slice(x_off, x_off + x_s)
            
            for k_m, k_r in keys.items():
                if (y_off + y_s <= master_dict[k_m].shape[0]) and \
                   (x_off + x_s <= master_dict[k_m].shape[1]):
                    master_dict[k_m][y_slice, x_slice] += np.nan_to_num(chunk_res[k_r])

        ## Scan Children
        for entry in children:
            ## We must yield points from the child to grid them
            ## This triggers re-extraction of the zip contents for this entry
            for points in entry.yield_points():
                if pp_grid:
                    res, srcwin, _ = pp_grid(points, mode='sums')
                    if res['z'] is not None and srcwin is not None:
                        accumulate(grid_arrays, res, srcwin, {'sum': 'z', 'count': 'count'})
                
                if pp_block:
                    res, srcwin, _ = pp_block(points, mode='sums')
                    if res['z'] is not None and srcwin is not None:
                        accumulate(block_arrays, res, srcwin, {
                            'z_sum': 'z', 'x_sum': 'x', 'y_sum': 'y', 'count': 'count'
                        })

        ## Finalize Mini-Grid
        if make_grid and grid_arrays is not None:
            with np.errstate(divide='ignore', invalid='ignore'):
                final_grid = grid_arrays['sum'] / grid_arrays['count']
            self.infos.mini_grid = np.where(np.isnan(final_grid), None, final_grid).tolist()
            
            ## Update WKT with tighter grid footprint
            if self.infos.mini_grid:
                self.infos.wkt = self._wkt_from_mini_grid(region, self.infos.mini_grid)

        ## Finalize Block-Mean
        if make_block_mean and block_arrays is not None:
            base = os.path.splitext(self.fn)[0]
            block_out = f"{base}_blockmean.xyz"
            try:
                with np.errstate(divide='ignore', invalid='ignore'):
                    fz = block_arrays['z_sum'] / block_arrays['count']
                    fx = block_arrays['x_sum'] / block_arrays['count']
                    fy = block_arrays['y_sum'] / block_arrays['count']
                
                valid = (block_arrays['count'] > 0) & np.isfinite(fz)
                with open(block_out, 'w') as f:
                    for x, y, z in zip(fx[valid], fy[valid], fz[valid]):
                        f.write(f"{x:.6f} {y:.6f} {z:.6f}\n")
            except Exception:
                pass

            
    def parse(self):
        """Parse the zipfile and yield supported datasets.
        """

        from cudem.datalists.dlim import DatasetFactory
        
        ## Determine valid extensions from Factory
        valid_exts = set()
        for mod in DatasetFactory._modules.values():
            if 'fmts' in mod:
                valid_exts.update(mod['fmts'])

        ## Identify files to extract
        files_to_process = []
        if self.fn.lower().endswith('.zip'):
            try:
                with zipfile.ZipFile(self.fn) as z:
                    for zf in z.namelist():
                        ext = zf.split('.')[-1]
                        if ext in valid_exts:
                            files_to_process.append(os.path.basename(zf))
            except Exception as e:
                utils.echo_error_msg(f'Could not unzip {self.fn}: {e}')
                return

        ## Extract and Process
        ## We extract files individually to temporary locations to handle large zip archives
        for filename in files_to_process:
            try:
                extracted_path = utils.p_f_unzip(
                    self.fn,
                    fns=[filename],
                    outdir=os.path.normpath(os.path.dirname(self.fn)),
                    tmp_fn=True
                )[0]
                
                ## Instantiate Dataset via Factory
                ds = DatasetFactory(
                    **self._set_params(
                        mod=extracted_path,
                        data_format=None,
                        src_srs=self.src_srs,
                        parent=self
                    )
                )._acquire_module()
                
                # ds = DatasetFactory(
                #     mod=extracted_path,
                #     data_format=None, # Auto-detect
                #     src_srs=self.src_srs,
                #     parent=self,
                #     weight=self.weight,
                #     uncertainty=self.uncertainty,
                #     verbose=self.verbose
                # )._acquire_module()

                if ds is not None and ds.valid_p():
                    ds.initialize()
                    
                    ## Region Intersection Check
                    if self.region is not None and self.region.valid_p(check_xy=True):
                        ## Load INF to check bounds before full parsing
                        ds.inf(check_hash=False) 
                        
                        check_region = ds.transform['trans_region'] if ds.transform['trans_region'] else ds.region
                        if check_region is None:
                             try:
                                 check_region = regions.Region().from_list(ds.infos.minmax)
                             except:
                                 check_region = None

                        if check_region and regions.regions_intersect_p(check_region, self.region):
                            for sub_ds in ds.parse():
                                self.data_entries.append(sub_ds)
                                yield sub_ds
                    else:
                        for sub_ds in ds.parse():
                            self.data_entries.append(sub_ds)
                            yield sub_ds
                            
                ## Cleanup extracted file
                utils.remove_glob(f'{extracted_path}*')
                
            except Exception as e:
                if self.verbose:
                    utils.echo_warning_msg(f"Failed to process {filename} in {self.fn}: {e}")
                continue

### End
