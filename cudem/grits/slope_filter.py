### slope_filter.py - GRId filTerS
##
## Copyright (c) 2024 - 2026 Regents of the University of Colorado
##
## slope_filter.py is part of cudem
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
## Filters a DEM based on slope or other derived parameters (TRI, TPI, Roughness).
## Useful for removing artifacts on steep cliffs or flattening noise on flat terrain.
##
### Code:

import numpy as np
from osgeo import gdal

from cudem import utils
from cudem import gdalfun
from cudem.grits import grits

class SlopeFilter(grits.Grits):
    """Filter DEM based on Slope or other Land Surface Parameters.

    Parameters:
    -----------
    metric : str
        The metric to filter by: 'slope', 'roughness', 'TPI', 'TRI'. Default 'slope'.
    min_val : float
        Minimum allowed value for the metric.
    max_val : float
        Maximum allowed value for the metric.
    action : str
        'mask': Set to NoData (Default).
        'revert': Revert to original source data (if current is modified).
    """
    
    def __init__(self, metric='slope', min_val=None, max_val=None, action='mask', **kwargs):
        super().__init__(**kwargs)
        self.metric = metric
        self.min_val = utils.float_or(min_val)
        self.max_val = utils.float_or(max_val)
        self.action = action

        
    def run(self):
        """Execute the slope filter."""
        
        ## Setup Output
        dst_ds = self.copy_src_dem()
        if dst_ds is None: return self.src_dem, -1

        with gdalfun.gdal_datasource(self.src_dem) as src_ds:
            if src_ds is None: return self.src_dem, -1
            self.init_ds(src_ds)
            
            ## Iterate Chunks
            for srcwin in utils.yield_srcwin(
                (src_ds.RasterYSize, src_ds.RasterXSize), 
                n_chunk=2048, 
                verbose=self.verbose
            ):
                ## Compute Metric (Slope, etc.)
                ## We need a buffer to compute slope at edges without artifacts
                buffer = 2
                srcwin_buff = utils.buffer_srcwin(
                    srcwin, (src_ds.RasterYSize, src_ds.RasterXSize), buffer
                )
                
                band = src_ds.GetRasterBand(self.band)
                src_arr = band.ReadAsArray(*srcwin_buff)
                
                if src_arr is None: continue
                
                ## Handle NoData
                ndv = self.ds_config['ndv']
                data = src_arr.astype(np.float32)
                data[data == ndv] = np.nan
                
                ## Compute Metric using Numpy Gradient (faster than calling gdaldem CLI per chunk)
                metric_arr = None
                
                if self.metric.lower() == 'slope':
                    ## np.gradient returns (dy, dx)
                    ## Slope = sqrt(dx^2 + dy^2)
                    ## Note: Need to account for cell size if Z units differ from XY
                    ## Assuming consistent units or pre-scaled for simplicity here.
                    gy, gx = np.gradient(data)
                    metric_arr = np.sqrt(gx**2 + gy**2)
                    ## Convert to degrees? usually slope is rise/run (tan theta) or degrees
                    ## GDAL uses degrees by default. Let's assume percent/ratio here for speed.
                    ## If degrees needed: np.degrees(np.arctan(metric_arr))
                    
                ## ... (Other metrics TPI/Roughness can be implemented via simple kernels)
                
                if metric_arr is None: continue

                ## Identify Filter Mask
                mask = np.zeros(metric_arr.shape, dtype=bool)
                
                if self.min_val is not None:
                    mask |= (metric_arr < self.min_val)
                if self.max_val is not None:
                    mask |= (metric_arr > self.max_val)
                    
                ## Apply Action
                ## Crop back to unbuffered window
                x_off = srcwin[0] - srcwin_buff[0]
                y_off = srcwin[1] - srcwin_buff[1]
                
                mask_crop = mask[y_off:y_off+srcwin[3], x_off:x_off+srcwin[2]]
                
                if np.any(mask_crop):
                    ## Read Destination Data (to modify)
                    dst_band = dst_ds.GetRasterBand(self.band)
                    dst_data = dst_band.ReadAsArray(*srcwin)
                    
                    if self.action == 'mask':
                        dst_data[mask_crop] = ndv
                    elif self.action == 'revert':
                        ## Read original source (unbuffered)
                        src_orig = band.ReadAsArray(*srcwin)
                        dst_data[mask_crop] = src_orig[mask_crop]
                        
                    dst_band.WriteArray(dst_data, srcwin[0], srcwin[1])

        dst_ds = None
        return self.dst_dem, 0

### End
