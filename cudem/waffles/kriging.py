### kriging.py
##
## Copyright (c) 2025 - 2026 Regents of the University of Colorado
##
## kriging.py is part of CUDEM
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
## Geostatistical Gridding via PyKrige.
##
### Code:

import numpy as np
from osgeo import gdal

from cudem import utils
from cudem import gdalfun
from cudem.waffles.waffles import Waffle

## Optional Dependency
try:
    from pykrige.ok import OrdinaryKriging
    HAS_PYKRIGE = True
except ImportError:
    HAS_PYKRIGE = False

    
class WafflesKriging(Waffle):
    """GEOSTATISTICAL DEM via Kriging (PyKrige)
    
    Generate a DEM using Ordinary Kriging.
    This method provides both an interpolated surface and an estimation 
    variance (uncertainty) surface.

    Note: Kriging is computationally expensive ($O(N^3)$). 
    For large datasets, use 'max_points' to decimate the input training data 
    or ensure your 'chunk_size' is small.

    Parameters:
    -----------
    model (str) : Variogram model ('linear', 'power', 'gaussian', 'spherical', 
                  'exponential', 'hole-effect'). Default: 'linear'.
    nlags (int) : Number of averaging bins for the variogram. Default: 6.
    max_points (int) : Maximum number of input points to use for training 
                       (randomly decimated if exceeded). Default: 10000.
    chunk_size (int) : Processing chunk size in pixels. Default: 512.

    < kriging:model=linear:nlags=6:max_points=10000 >
    """
    
    def __init__(
            self,
            model='linear',
            nlags=6,
            max_points=10000,
            chunk_size=None,
            **kwargs
    ):
        super().__init__(**kwargs)
        self.model = model
        self.nlags = utils.int_or(nlags, 6)
        self.max_points = utils.int_or(max_points, 10000)
        self.chunk_size = chunk_size
        self.chunk_step = None

        
    def run(self):
        if not HAS_PYKRIGE:
            utils.echo_error_msg("PyKrige must be installed to use the KRIGING module.")
            return self

        if self.verbose:
            utils.echo_msg(
                f'Generating Kriging Grid ({self.model}) @ {self.ycount}/{self.xcount} '
                f'using max {self.max_points} training points...'
            )

        if self.chunk_size is None:
            n_chunk = 512
        else:
            n_chunk = self.chunk_size

        ## Load Data from Stack
        ## We read the stack to get X, Y, Z coordinates for training the Krige model.
        with gdalfun.gdal_datasource(self.stack) as stack_ds:
            points_band = stack_ds.GetRasterBand(1)
            points_no_data = points_band.GetNoDataValue()
            
            ## Read Array
            points_array = points_band.ReadAsArray()
            
            ## Extract valid points
            valid_mask = points_array != points_no_data
            if np.count_nonzero(valid_mask) == 0:
                utils.echo_warning_msg("No valid data points found in stack.")
                return self

            ## Create Training Coordinate Arrays
            ## Convert grid indices to georeferenced coordinates
            gt = stack_ds.GetGeoTransform()
            py, px = np.nonzero(valid_mask)
            
            z_train = points_array[py, px]
            #x_train, y_train = utils._pixel2geo(px, py, gt, node='pixel')
            x_train, y_train = px, py
            
            ## Decimate if too many points to prevent memory crash
            if len(z_train) > self.max_points:
                if self.verbose:
                    utils.echo_msg(f"Decimating input from {len(z_train)} to {self.max_points} points.")
                indices = np.random.choice(len(z_train), self.max_points, replace=False)
                x_train = x_train[indices]
                y_train = y_train[indices]
                z_train = z_train[indices]

            ## Prepare Output Dataset
            driver = gdal.GetDriverByName(self.fmt)
            dst_ds = driver.Create(
                self.fn,
                self.xcount,
                self.ycount,
                2, # Band 1: Elevation, Band 2: Variance
                gdal.GDT_Float32,
                options=self.co
            )
            
            dst_ds.SetGeoTransform(self.dst_gt)
            dst_ds.SetProjection(gdalfun.osr_wkt(self.dst_srs))
            
            elev_band = dst_ds.GetRasterBand(1)
            elev_band.SetNoDataValue(self.ndv)
            elev_band.SetDescription("Elevation")
            
            var_band = dst_ds.GetRasterBand(2)
            var_band.SetNoDataValue(self.ndv)
            var_band.SetDescription("Kriging Variance")

            ## Train Kriging Model
            if self.verbose:
                utils.echo_msg("Training Variogram model...")

            try:
                OK = OrdinaryKriging(
                    x_train, 
                    y_train, 
                    z_train, 
                    variogram_model=self.model,
                    nlags=self.nlags,
                    verbose=True,
                    enable_plotting=False
                )
            except Exception as e:
                utils.echo_error_msg(f"Kriging training failed: {e}")
                return self

            # ## Execute Interpolation (Chunked)
            # ## PyKrige execute('grid') takes unique 1D arrays for x and y axes.
            # ## We iterate through srcwin chunks to manage memory.            
            # for srcwin in utils.yield_srcwin(
            #         (self.ycount, self.xcount),
            #         n_chunk=n_chunk,
            #         msg='Kriging Interpolation',
            #         verbose=self.verbose
            # ):
            #     ## Calculate coordinate vectors for this chunk
            #     x_off, y_off, x_size, y_size = srcwin
                
            #     ## Get the Geo coordinates for the chunk axes
            #     ## Pixel centers
            #     #x_start, y_start = utils._pixel2geo(x_off, y_off, self.dst_gt, node='pixel')
            #     #x_end, y_end = utils._pixel2geo(x_off + x_size, y_off + y_size, self.dst_gt, node='pixel')
                
            #     ## Generate 1D coordinate arrays for PyKrige
            #     ## Note: linspace endpoint=False is slightly safer for raster align
            #     ## But PyKrige expects exact coords. Using arange/linspace based on GT.
                
            #     ## X coordinates (Columns)
            #     #grid_x = np.linspace(x_start, x_start + (x_size * self.dst_gt[1]), x_size, endpoint=False)
            #     #grid_x = np.linspace(x_off, x_off + x_size, x_size)
            grid_x = np.linspace(0, self.xcount, self.xcount)
                
            #     ## Y coordinates (Rows)
            #     ## Note: dst_gt[5] is usually negative.
            #     #grid_y = np.linspace(y_start, y_start + (y_size * self.dst_gt[5]), y_size, endpoint=False)
            #     #grid_y = np.linspace(y_off, y_off + y_size, y_size)
            grid_y = np.linspace(0, self.ycount, self.ycount)

            #     # ## srcwin: (x_off, y_off, x_size, y_size)
            #     # xi, yi = np.mgrid[srcwin[0]:srcwin[0] + srcwin[2],
            #     #                   srcwin[1]:srcwin[1] + srcwin[3]]

            #     #utils.echo_msg(f'{grid_x}, {grid_y}')
            #     ## Execute Kriging
            #     ## returns: z_values, sigma_squared (variance)
            try:
                z_chunk, ss_chunk = OK.execute('grid', grid_x, grid_y, backend='loop', n_closest_points=10)
                
                ## PyKrige returns masked arrays if points fall outside domain, 
                ## fill with NDV.
                # if np.ma.is_masked(z_chunk):
                #     z_chunk = z_chunk.filled(self.ndv)
                # if np.ma.is_masked(ss_chunk):
                #     ss_chunk = ss_chunk.filled(self.ndv)
                
                utils.echo_msg(z_chunk)
                    
                ## Write to disk
                elev_band.WriteArray(z_chunk)#, x_off, y_off)
                var_band.WriteArray(ss_chunk)#, x_off, y_off)
                    
            except Exception as e:
                utils.echo_warning_msg(f"Krige failed; {e}")
                    
            dst_ds = None
            
        return self

### End
