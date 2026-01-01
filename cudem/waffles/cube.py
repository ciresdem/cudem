### cube.py
##
## Copyright (c) 2025 - 2026 Regents of the University of Colorado
##
## cube.py is part of CUDEM
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
## CUBE (Combined Uncertainty and Bathymetry Estimator) module.
## Wraps the 'bathycube' library.
##
### Code:

import numpy as np
from osgeo import gdal

from cudem import utils
from cudem import gdalfun
from cudem.waffles.waffles import Waffle

class WafflesCUBE(Waffle):
    """
    BathyCUBE interpolation module.
    https://github.com/noaa-ocs-hydrography/bathycube
    
    Parameters:
    -----------
    method (str) : CUBE method ['local', 'propagate']. Default: 'local'.
    order (str) : CUBE order ['order1a', ...]. Default: 'order1a'.
    """
    
    def __init__(
            self,
            method='local',
            order='order1a',
            chunk_size=None,
            **kwargs):
        """Generate a `CUBE` DEM."""
        
        super().__init__(**kwargs)
        self.method = method
        self.order = order
        self.chunk_size = utils.int_or(chunk_size)

        
    def run(self):
        ## Check for library availability
        try:
            import bathycube.cube as cube
        except ImportError:
            utils.echo_error_msg(
                'Could not import bathycube. Please install it from: '
                'https://github.com/noaa-ocs-hydrography/bathycube'
            )
            return self

        if self.verbose:
            utils.echo_msg(f'Running CUBE ({self.method}/{self.order})...')

        ## Open the Stack (Input Data)
        ds = gdal.Open(self.stack)
        if ds is None:
            utils.echo_error_msg(f'Could not open stack file: {self.stack}')
            return self

        ## Setup Output Dataset (Elevation)
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
        elev_band.SetDescription("CUBE Depth")
        
        ## Prepare Data Arrays
        z_band = ds.GetRasterBand(1) # Z
        u_band = ds.GetRasterBand(4) # Uncertainty (TVU)
        
        z_arr = z_band.ReadAsArray()
        u_arr = u_band.ReadAsArray()
        
        ## Mask NoData
        ndv = z_band.GetNoDataValue()
        mask = z_arr != ndv
        
        if np.count_nonzero(mask) == 0:
            utils.echo_warning_msg("No valid points in stack.")
            return self

        ## Extract valid points
        z_vals = z_arr[mask]
        tvu_vals = u_arr[mask]
        
        ## Convert Pixel Indices to Georeferenced Coordinates
        py, px = np.nonzero(mask)
        x_vals, y_vals = utils._pixel2geo(px, py, self.dst_gt, node='pixel')

        ## Generate generic THU if not available
        thu_vals = np.full(z_vals.shape, self.xinc * 0.5)

        ## Run CUBE
        try:
            depth_grid, unc_grid, ratio_grid, numhyp_grid = cube.run_cube_gridding(
                z_vals, 
                thu_vals, 
                tvu_vals, 
                x_vals, 
                y_vals, 
                self.ycount,
                self.xcount, 
                self.dst_gt[0], # min_x
                self.dst_gt[3], # max_y
                self.method, 
                self.order, 
                abs(self.xinc), 
                abs(self.yinc)
            )
            
            ## Write Depth
            depth_grid = np.nan_to_num(depth_grid, nan=self.ndv)
            elev_band.WriteArray(depth_grid)
            
            ## Write Uncertainty (Optional)
            if self.want_uncertainty:
                unc_fn = f'{self.name}_u.tif'
                if self.verbose:
                    utils.echo_msg(f'Writing CUBE uncertainty to {unc_fn}...')
                    
                unc_ds = driver.Create(
                    unc_fn,
                    self.xcount,
                    self.ycount,
                    1,
                    gdal.GDT_Float32,
                    options=self.co
                )
                unc_ds.SetGeoTransform(self.dst_gt)
                unc_ds.SetProjection(gdalfun.osr_wkt(self.dst_srs))
                
                unc_band = unc_ds.GetRasterBand(1)
                unc_band.SetNoDataValue(self.ndv)
                unc_band.SetDescription("CUBE Uncertainty")
                
                unc_grid = np.nan_to_num(unc_grid, nan=self.ndv)
                unc_band.WriteArray(unc_grid)
                
                unc_ds = None # Close/Flush
            
        except Exception as e:
            utils.echo_error_msg(f"CUBE execution failed: {e}")

        ## Cleanup
        ds = None
        dst_ds = None
        
        return self

    
### End
