### natural_neighbor.py
##
## Copyright (c) 2010 - 2026 Regents of the University of Colorado
##
## natural_neighbor.py is part of CUDEM
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
## Natural Neighbor (Sibson) Interpolation module.
## Uses the 'naturalneighbor' python package (Discrete Sibson).
##
### Code:

import numpy as np
from osgeo import gdal

from cudem import utils
from cudem import gdalfun
from cudem.waffles.waffles import Waffle

## Optional Dependency
try:
    import naturalneighbor
    HAS_NATURALNEIGHBOR = True
except ImportError:
    HAS_NATURALNEIGHBOR = False

class WafflesNaturalNeighbor(Waffle):
    """NATURAL NEIGHBOR DEM via Discrete Sibson Interpolation.
    
    Generates a DEM using Natural Neighbor interpolation. 
    This method is "smoother" than linear/nearest interpolation and does 
    not require tuning parameters (like power/radius in IDW).
    
    It handles irregular data clusters well but stays strictly within the 
    convex hull of the input data (no extrapolation).

    Parameters:
    -----------
    max_points (int) : Maximum number of input points to use. Default: 50000.
    depth (float) : The 'depth' of the discrete grid for the library (rarely needs tuning). Default: 1.
    chunk_size (int) : Processing chunk size in pixels. Default: 512.

    < natural_neighbor:max_points=50000 >
    """
    
    def __init__(
            self,
            max_points=50000,
            depth=1,
            chunk_size=None,
            **kwargs
    ):
        super().__init__(**kwargs)
        self.max_points = utils.int_or(max_points, 50000)
        self.depth = utils.float_or(depth, 1)
        self.chunk_size = chunk_size
        self.chunk_step = None

        
    def run(self):
        if not HAS_NATURALNEIGHBOR:
            utils.echo_error_msg("The 'naturalneighbor' package must be installed to use this module.")
            utils.echo_error_msg("Try: pip install naturalneighbor")
            return self

        if self.verbose:
            utils.echo_msg(
                f'Generating Natural Neighbor Grid @ {self.ycount}/{self.xcount} '
                f'using max {self.max_points} points...'
            )

        if self.chunk_size is None:
            n_chunk = 512
        else:
            n_chunk = self.chunk_size

        ## Load Data from Stack
        with gdalfun.gdal_datasource(self.stack) as stack_ds:
            points_band = stack_ds.GetRasterBand(1)
            points_no_data = points_band.GetNoDataValue()
            
            ## Read all points into memory
            points_array = points_band.ReadAsArray()
            
            ## Mask valid data
            valid_mask = points_array != points_no_data
            if np.count_nonzero(valid_mask) == 0:
                utils.echo_warning_msg("No valid data points found in stack.")
                return self

            ## Convert to coordinates
            gt = stack_ds.GetGeoTransform()
            py, px = np.nonzero(valid_mask)
            
            z_train = points_array[py, px]
            x_train, y_train = utils._pixel2geo(px, py, gt, node='pixel')
            
            ## Decimate if necessary
            if len(z_train) > self.max_points:
                if self.verbose:
                    utils.echo_msg(f"Decimating input from {len(z_train)} to {self.max_points} points.")
                indices = np.random.choice(len(z_train), self.max_points, replace=False)
                x_train = x_train[indices]
                y_train = y_train[indices]
                z_train = z_train[indices]

            ## Prepare Input Points (Nx3 array) for naturalneighbor
            ## [[x, y, z], ...]
            input_points = np.column_stack((x_train, y_train, z_train))

            ## Prepare Output Dataset
            driver = gdal.GetDriverByName(self.fmt)
            dst_ds = driver.Create(
                self.fn,
                self.xcount,
                self.ycount,
                1,
                gdal.GDT_Float32,
                options=self.co
            )
            
            dst_ds.SetGeoTransform(self.dst_gt)
            dst_ds.SetProjection(gdalfun.osr_wkt(self.dst_srs))
            
            elev_band = dst_ds.GetRasterBand(1)
            elev_band.SetNoDataValue(self.ndv)
            elev_band.SetDescription("Natural Neighbor Interpolation")

            ## Interpolate (Chunked)
            ## The 'naturalneighbor' library supports grid interpolation defined by ranges.
            
            for srcwin in utils.yield_srcwin(
                    (self.ycount, self.xcount),
                    n_chunk=n_chunk,
                    msg='Natural Neighbor Interpolation',
                    verbose=self.verbose
            ):
                x_off, y_off, x_size, y_size = srcwin
                
                ## Define the grid range for this chunk
                ## [[start, stop, steps], [start, stop, steps], [start, stop, steps]]
                ## Note: Z range (3rd dim) is handled internally by the library for 3D, 
                ## but for 2D interp we pass grid_ranges for X/Y.
                
                ## Get Extents
                x_start, y_start = utils._pixel2geo(x_off, y_off, self.dst_gt, node='pixel')
                ## Determine end based on pixel count and resolution
                ## Ensure directionality matches GT
                x_end = x_start + (x_size * self.dst_gt[1])
                y_end = y_start + (y_size * self.dst_gt[5])
                
                ## Construct Grid Ranges
                ## naturalneighbor expects: [[min, max, step], [min, max, step]]
                ## Note: It handles ordering, so we ensure min < max
                grid_ranges = [
                    [min(x_start, x_end), max(x_start, x_end), abs(self.dst_gt[1])],
                    [min(y_start, y_end), max(y_start, y_end), abs(self.dst_gt[5])],
                ]

                try:
                    ## Execute Interpolation
                    ## The library returns a 3D grid if we give it 3D points, 
                    ## but we are doing 2D surface interpolation (Z=f(X,Y)).
                    ## Actually, naturalneighbor.griddata takes known points and target ranges.
                    
                    z_chunk = naturalneighbor.griddata(
                        input_points, 
                        grid_ranges
                    )
                    
                    ## The output shape might need squeezing or transposition depending on library version
                    ## Usually returns (Y, X) or (X, Y)
                    
                    ## If the library returns a result larger than requested (due to step alignment),
                    ## crop it to match x_size/y_size
                    if z_chunk.shape[0] != y_size or z_chunk.shape[1] != x_size:
                        ## Simple resize/crop logic if floating point errors cause mismatch
                        z_chunk = z_chunk[:y_size, :x_size]

                    ## Fill NaNs with NDV
                    z_chunk = np.nan_to_num(z_chunk, nan=self.ndv)
                    
                    ## Write
                    elev_band.WriteArray(z_chunk, x_off, y_off)
                    
                except Exception as e:
                    ## If chunk falls outside convex hull, it might return all NaNs or fail
                    if self.verbose:
                        utils.echo_warning_msg(f"Chunk failed or empty @ {srcwin}: {e}")

            dst_ds = None
            
        return self


### End
