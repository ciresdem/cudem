### globatofile.py - DataLists IMproved
##
## Copyright (c) 2010 - 2026 Regents of the University of Colorado
##
## globatofile.py is part of CUDEM
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
## CUDEM Globato file parsing...
##
### Code:

import h5py
import numpy as np
from cudem import utils
from cudem import regions
from cudem.datalists.dlim import ElevationDataset

class GlobatoFile(ElevationDataset):
    """Parser for Globato HDF5 Stack files.
    Reads the 'stack' group (finalized data) or 'datasets' group.
    """
    
    def __init__(self, layer='z', group='stack', use_group_xy=False, **kwargs):
        super().__init__(**kwargs)
        self.layer = layer # 'z', 'uncertainty', 'weight', etc.
        self.group = group # 'stack' or 'datasets'
        self.use_group_xy = use_group_xy # If True, prefer 'x'/'y' datasets in the group over global lat/lon
        self.data_format = 320

        
    def generate_inf(self, make_grid=True, make_block_mean=False, block_inc=None):
        """Generate metadata (INF) for the Globato file.
        
        This overrides the base implementation to optimize the initial bounds scan
        using Globato HDF5 attributes (if available), but still supports 
        mini-grid and block-mean generation if requested.
        """

        ## ==============================================
        ## Quick Metadata Parse from HDF5 Attributes
        ## Globato files store GeoTransform in the 'crs' group attributes.
        ## This is much faster than scanning points.
        ##
        ## If use_group_xy is True, we CANNOT rely on global attributes
        ## because the group's specific X/Y arrays might differ (e.g. node vs pixel).
        ## We must skip this optimization and scan the actual points.
        ## ==============================================
        found_attrs = False
        
        if not self.use_group_xy:
            try:
                with h5py.File(self.fn, 'r') as f:
                    if 'crs' in f and 'GeoTransform' in f['crs'].attrs:
                        gt = [float(x) for x in f['crs'].attrs['GeoTransform'].split()]
                        ## GT: [minx, x_inc, 0, maxy, 0, y_inc]
                        
                        target_shape = None
                        # Determine shape from stack or first dataset
                        if 'stack' in f and 'z' in f['stack']:
                            target_shape = f['stack']['z'].shape # (y, x)
                        elif 'datasets' in f:
                            keys = list(f['datasets'].keys())
                            ## If the first key is a data group, grab shape. 
                            ## If nested, this simple check might fail, but full scan fallback handles it.
                            if keys and 'z' in f['datasets'][keys[0]]:
                                target_shape = f['datasets'][keys[0]]['z'].shape
                        
                        if target_shape:
                            shape = target_shape
                            minx = gt[0]
                            maxy = gt[3]
                            maxx = minx + (gt[1] * shape[1])
                            miny = maxy + (gt[5] * shape[0])
                            
                            ## Set basic INF properties
                            self.infos.minmax = [minx, maxx, miny, maxy, 0, 0] # Z is unknown without scan, set 0
                            self.infos.numpts = shape[0] * shape[1]
                            
                            ## SRS
                            if 'crs_wkt' in f['crs'].attrs:
                                self.infos.src_srs = f['crs'].attrs['crs_wkt'].decode('utf-8')
                                
                            found_attrs = True

            except Exception as e:
                if self.verbose:
                    utils.echo_warning_msg(f"Could not read Globato attributes: {e}")

        ## ==============================================
        ## Full Scan / Grid Generation
        ## If attributes were missing OR if grids are requested OR if use_group_xy is True,
        ## we fall back to the base class logic which scans the points.
        ## ==============================================
        if not found_attrs or make_grid or make_block_mean:
            return super().generate_inf(
                make_grid=make_grid, 
                make_block_mean=make_block_mean, 
                block_inc=block_inc
            )
            
        return self.infos


    def _yield_group_points(self, grp, lats, lons):
        """Internal generator to yield points from a specific h5 group/dataset recursively.
        
        This traverses the HDF5 group structure until it finds the group containing 
        the target data layer (self.layer).
        """
        
        ## The current group contains the target data layer
        if self.layer in grp:
            try:
                data = grp[self.layer][...] # 2D array
                
                ## Handle variable naming differences between stack and datasets
                ## stack usually uses 'weights', datasets uses 'weight'
                weights = None
                if 'weights' in grp: 
                    weights = grp['weights'][...]
                elif 'weight' in grp: 
                    weights = grp['weight'][...]

                uncertainty = None
                if 'uncertainty' in grp: 
                    uncertainty = grp['uncertainty'][...]

                counts = None
                if 'count' in grp: 
                    counts = grp['count'][...]

                ## Check for local XY override (Node positions from stacker)
                local_x = None
                local_y = None
                if self.use_group_xy and 'x' in grp and 'y' in grp:
                    local_x = grp['x'][...]
                    local_y = grp['y'][...]

                #height, width = data.shape

                local_z = None
                if 'z' in grp:
                    local_z = grp['z'][...]

                valid_mask = counts > 0
                weights = weights / counts
                local_x = local_x / weights / counts
                local_y = local_y / weights / counts
                local_z = local_z / weights / counts
                
                ds = np.rec.fromarrays(
                    [local_x[valid_mask], local_y[valid_mask], local_z[valid_mask],
                     weights[valid_mask], uncertainty[valid_mask]],
                    names=['x', 'y', 'z', 'w', 'u']
                )
                yield ds
                
                # ## Iterate scanlines (chunks) to yield points
                # for i in range(height):
                #     row_data = data[i, :]
                #     row_valid = ~np.isnan(row_data)
                    
                #     if not np.any(row_valid): continue
                    
                #     ## --- Coordinate Processing ---
                #     if local_x is not None and local_y is not None:
                #         ## Use group-specific coordinates
                #         ## If 2D (swath/stack-node), slice the row
                #         if local_x.ndim == 2:
                #             x_vals = local_x[i, :][row_valid]
                #             y_vals = local_y[i, :][row_valid]
                #         ## If 1D (grid axes), broadcast Y
                #         else:
                #             y_val = local_y[i]
                #             x_vals = local_x[row_valid]
                #             y_vals = np.full(x_vals.shape, y_val)
                #     else:
                #         ## Use global lat/lon vectors
                #         y_val = lats[i]
                #         x_vals = lons[row_valid]
                #         y_vals = np.full(x_vals.shape, y_val)
                    
                #     z_vals = row_data[row_valid]
                    
                #     ## --- Weights & Uncertainty ---
                #     if weights is not None:
                #         w_vals = weights[i, :][row_valid]
                #     else:
                #         w_vals = np.ones(x_vals.shape)
                    
                #     if uncertainty is not None:
                #         u_vals = uncertainty[i, :][row_valid]
                #     else:
                #         u_vals = np.zeros(x_vals.shape)

                #     if counts is not None:
                #         c_vals = counts[i, :][row_valid]
                #     else:
                #         c_vals = np.zeros(x_vals.shape)

                #     w_vals = w_vals / c_vals
                #     x_vals = x_vals / w_vals / c_vals
                #     y_vals = y_vals / w_vals / c_vals
                #     z_vals = z_vals / w_vals / c_vals
                        
                #     ds = np.rec.fromarrays(
                #         [x_vals, y_vals, z_vals, w_vals, u_vals],
                #         names=['x', 'y', 'z', 'w', 'u']
                #     )
                #     yield ds

            except Exception as e:
                utils.echo_error_msg(f"Error parsing group {grp.name}: {e}")
                utils.echo_msg(grp.keys())

        ## Recursive Step: Search subgroups
        ## If the layer isn't here, check if this is a group and recurse into children
        elif isinstance(grp, h5py.Group):
            for key in grp.keys():
                yield from self._yield_group_points(grp[key], lats, lons)

                
    def yield_points(self):
        """Yields non-nan points from the Globato stack."""
        
        with h5py.File(self.fn, 'r') as f:
            try:
                ## Globato stores 1D lat/lon arrays at root
                lats = f['lat'][...]
                lons = f['lon'][...]
                
                if self.group == 'datasets':
                    if 'datasets' in f:
                        for ds_name in f['datasets']:
                            ## Yield points from each dataset sub-group
                            yield from self._yield_group_points(f['datasets'][ds_name], lats, lons)
                else:
                    ## Default to stack
                    if 'stack' in f:
                        yield from self._yield_group_points(f['stack'], lats, lons)

            except Exception as e:
                utils.echo_error_msg(f"Globato read error: {e}")

### End
