### diff.py - GRId filTerS
##
## Copyright (c) 2025 - 2026 Regents of the University of Colorado
##
## diff.py is part of cudem
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
## Calculates difference between Source and Auxiliary DEM.
## Usage:
## * Filter: Remove Source pixels where abs(Source - Aux) > Threshold.
## * Output: Generate a difference grid (Output = Source - Aux).
##
### Code:

import numpy as np
from osgeo import gdal

from cudem import utils
from cudem import gdalfun
from cudem.grits import grits

class Diff(grits.Grits):
    """Difference Filter.
    
    Parameters:
    -----------
    aux_dem : str
        Path to the comparison/reference DEM.
    diff_threshold : float
        Maximum allowed difference. Pixels exceeding this are masked.
        If None, the module outputs the difference grid instead of filtering.
    mode : str
        'absolute': Filter based on abs(src - aux).
        'relative': Filter based on (src - aux).
    """
    
    def __init__(self, aux_dem=None, diff_threshold=None, mode='absolute', **kwargs):
        super().__init__(**kwargs)
        self.aux_dem = aux_dem
        self.diff_threshold = utils.float_or(diff_threshold)
        self.mode = mode

        # If no aux_dem provided in init, check if passed via aux_data list
        if self.aux_dem is None and len(self.aux_data) > 0:
            self.aux_dem = self.aux_data[0]

            
    def run(self):
        """Execute the difference filter."""
        
        if not self.aux_dem:
            utils.echo_error_msg("Diff filter requires an auxiliary DEM.")
            return self.src_dem, -1

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
                ## Read Source
                src_band = src_ds.GetRasterBand(self.band)
                src_arr = src_band.ReadAsArray(*srcwin).astype(np.float32)
                
                ## Load Aux Data (Aligned)
                ## Note: We need a Region object for load_aux_data.
                ## srcwin is (xoff, yoff, xsize, ysize)
                ## Calculate Region for this chunk
                gt = self.ds_config['geoT']
                
                x_min = gt[0] + (srcwin[0] * gt[1])
                x_max = x_min + (srcwin[2] * gt[1])
                y_max = gt[3] + (srcwin[1] * gt[5])
                y_min = y_max + (srcwin[3] * gt[5])
                
                ## Depending on gt[5] sign, ymin/ymax might be flipped logic
                ## Standard GDAL GT: top-left origin, dy is negative.
                chunk_region = regions.Region().from_list([x_min, x_max, min(y_min,y_max), max(y_min,y_max)])
                
                aux_arr = self.load_aux_data(
                    aux_source_list=[self.aux_dem],
                    src_region=chunk_region,
                    src_gt=gt, # Use source transform to ensure pixel alignment
                    x_count=srcwin[2],
                    y_count=srcwin[3]
                )
                
                if aux_arr is None: continue
                
                ## Handle NoData
                ndv = self.ds_config['ndv']
                src_arr[src_arr == ndv] = np.nan
                aux_arr[aux_arr == ndv] = np.nan # Assuming same NDV or handled by load_aux
                
                ## Calculate Difference
                diff_arr = src_arr - aux_arr
                
                ## Filtering (Masking)
                if self.diff_threshold is not None:
                    mask = np.zeros(diff_arr.shape, dtype=bool)
                    
                    if self.mode == 'absolute':
                        mask = np.abs(diff_arr) > self.diff_threshold
                    else:
                        mask = diff_arr > self.diff_threshold
                        
                    # Apply Mask to Destination
                    if np.any(mask):
                        dst_band = dst_ds.GetRasterBand(self.band)
                        dst_data = dst_band.ReadAsArray(*srcwin)
                        dst_data[mask] = ndv
                        dst_band.WriteArray(dst_data, srcwin[0], srcwin[1])
                        
                ## Output Difference Grid
                else:
                    diff_arr[np.isnan(diff_arr)] = ndv
                    dst_band = dst_ds.GetRasterBand(self.band)
                    dst_band.WriteArray(diff_arr, srcwin[0], srcwin[1])

        dst_ds = None
        return self.dst_dem, 0

### End
