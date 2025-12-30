### fill.py - GRId filTerS
##
## Copyright (c) 2024 - 2025 Regents of the University of Colorado
##
## fill.py is part of cudem
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
### Commentary:
##
## Module for filling NoData voids (inpainting) in DEMs.
## Supports GDAL's FillNodata (IDW) and Scipy's GridData (Spline/Linear).
##
### Code:

import numpy as np
from scipy import interpolate
from osgeo import gdal

from cudem import utils
from cudem import gdalfun
from cudem.grits import grits

class Fill(grits.Grits):
    """Fill NoData voids in the input DEM.

    Parameters:
    -----------
    method : str
        The filling method to use:
        - 'idw': Inverse Distance Weighting (via gdal.FillNodata). Best for small holes.
        - 'spline': Cubic spline interpolation. Good for smooth surfaces.
        - 'linear': Linear tessellation.
        - 'nearest': Nearest neighbor (blocky).
    max_dist : float
        The maximum distance (in pixels) to search for valid values to fill a hole.
        Default is 100.
    smoothing_iterations : int
        Number of smoothing iterations to apply after filling (IDW only).
    """
    
    def __init__(self, method='idw', max_dist=100, smoothing_iterations=0, **kwargs):
        super().__init__(**kwargs)
        self.method = method.lower() if method else 'idw'
        self.max_dist = utils.float_or(max_dist, 100)
        self.smoothing_iterations = utils.int_or(smoothing_iterations, 0)

        
    def _fill_idw(self, dst_ds):
        """Use GDAL's FillNodata algorithm.
        This operates directly on the GDAL Band.
        """
        
        if self.verbose:
            utils.echo_msg(f"Filling voids using GDAL IDW (Max Dist: {self.max_dist})...")
            
        band = dst_ds.GetRasterBand(self.band)
        
        ## gdal.FillNodata(targetBand, maskBand, maxSearchDist, smoothingIterations)
        ## If maskBand is None, it treats all non-NoData values as valid.
        gdal.FillNodata(
            targetBand=band,
            maskBand=None,
            maxSearchDist=self.max_dist,
            smoothingIterations=self.smoothing_iterations
        )
        band.FlushCache()

        
    def _fill_scipy(self, dst_ds):
        """Use Scipy griddata to interpolate voids.
        Processes in chunks to manage memory.
        """
        
        ## Mapping common names to scipy methods
        scipy_method = 'cubic' if self.method == 'spline' else self.method
        if scipy_method not in ['linear', 'cubic', 'nearest']:
            scipy_method = 'linear'

        if self.verbose:
            utils.echo_msg(f"Filling voids using Scipy {scipy_method} interpolation...")

        ## We need a buffer to ensure interpolation has context at chunk edges
        ## Buffer should be at least as large as the max hole radius we expect to fill
        buffer = int(self.max_dist)
        
        ## Iterate over chunks
        for srcwin in utils.yield_srcwin(
            (dst_ds.RasterYSize, dst_ds.RasterXSize), 
            n_chunk=2048, 
            verbose=self.verbose
        ):
            ## Buffer the window
            srcwin_buff = utils.buffer_srcwin(
                srcwin, (dst_ds.RasterYSize, dst_ds.RasterXSize), buffer
            )
            
            band = dst_ds.GetRasterBand(self.band)
            data = band.ReadAsArray(
                srcwin_buff[0], srcwin_buff[1], srcwin_buff[2], srcwin_buff[3]
            )
            
            if data is None: continue
            
            ndv = self.ds_config['ndv']
            
            ## Identify Voids and Valid Data
            is_valid = data != ndv
            is_void = data == ndv
            
            ## Optimization: If no voids in this chunk (including buffer), skip
            if not np.any(is_void):
                continue
                
            ## Optimization: If chunk is completely empty, we can't fill it
            if not np.any(is_valid):
                continue

            ## Check if there are voids specifically in the unbuffered center
            ## (If voids are only in the buffer, we don't need to process this chunk)
            y_off = srcwin[1] - srcwin_buff[1]
            x_off = srcwin[0] - srcwin_buff[0]
            center_slice = data[y_off:y_off+srcwin[3], x_off:x_off+srcwin[2]]
            
            if not np.any(center_slice == ndv):
                continue

            ## Prepare Interpolation Data
            ## Get coordinates of valid pixels
            valid_y, valid_x = np.nonzero(is_valid)
            valid_values = data[valid_y, valid_x]
            
            ## Get coordinates of void pixels we want to fill
            ## (We construct a grid for the whole chunk or just the voids?)
            ## griddata constructs the whole grid usually.
            grid_y, grid_x = np.mgrid[0:srcwin_buff[3], 0:srcwin_buff[2]]
            
            try:
                ## Interpolate
                ## Note: griddata can be slow for very large numbers of points.
                ## Linear/Cubic require qhull.
                filled_data = interpolate.griddata(
                    (valid_y, valid_x), 
                    valid_values, 
                    (grid_y, grid_x), 
                    method=scipy_method,
                    fill_value=ndv
                )
                
                ## We only want to update the VOIDS, preserving original valid data exactly
                ## (Interpolation might slightly alter valid points due to precision)
                ## Also, we only write back the unbuffered center region.
                
                out_arr = center_slice.copy()
                
                ## Extract fill result for center
                fill_center = filled_data[y_off:y_off+srcwin[3], x_off:x_off+srcwin[2]]
                
                ## Identify voids in the original center
                void_mask = (out_arr == ndv)
                
                ## Apply fill where original was void AND fill succeeded (not nan/ndv)
                ## Note: griddata returns nan for fill_value usually unless specified
                valid_fill = (fill_center != ndv) & (~np.isnan(fill_center))
                
                update_mask = void_mask & valid_fill
                out_arr[update_mask] = fill_center[update_mask]
                
                ## Write to disk
                band.WriteArray(out_arr, srcwin[0], srcwin[1])
                
            except Exception as e:
                if self.verbose:
                    utils.echo_warning_msg(f"Interpolation failed for chunk {srcwin}: {e}")

                    
    def run(self):
        """Execute the fill filter."""
        
        ## Create Output (Copy of Source)
        dst_ds = self.copy_src_dem()
        if dst_ds is None: return self.src_dem, -1

        ## Initialize config for NDV checks
        with gdalfun.gdal_datasource(self.src_dem) as src_ds:
            if src_ds: self.init_ds(src_ds)

        ## Select Method
        if self.method == 'idw':
            ## GDAL FillNodata handles buffering/chunking internally (mostly)
            ## but operates on the whole dataset object.
            self._fill_idw(dst_ds)
        else:
            ## Scipy based interpolation
            self._fill_scipy(dst_ds)

        dst_ds = None
        return self.dst_dem, 0

### End
