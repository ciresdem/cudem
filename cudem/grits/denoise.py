### denoise.py - GRId filTerS
##
## Copyright (c) 2025 - 2026 Regents of the University of Colorado
##
## denoise.py is part of cudem
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
## Denoising filters for DEMs.
## - Median: Removes salt-and-pepper noise (spikes/pits) while preserving edges.
## - Bilateral: Smooths surfaces while preserving sharp edges (requires scikit-image).
##
### Code:

import numpy as np
import scipy.ndimage
from osgeo import gdal

from cudem import utils
from cudem import gdalfun
from cudem.grits import grits

try:
    from skimage.restoration import denoise_bilateral
    HAS_SKIMAGE = True
except ImportError:
    HAS_SKIMAGE = False

    
class Denoise(grits.Grits):
    """Apply denoising filters to the DEM.

    Parameters:
    -----------
    method : str
        'median' or 'bilateral'. Default is 'median'.
    kernel_size : int
        Size of the filter kernel (pixels). For median filter. Default 3.
    sigma_spatial : float
        Standard deviation for range distance (Bilateral only). Default 1.
    """
    
    def __init__(self, method='median', kernel_size=3, sigma_spatial=1, **kwargs):
        super().__init__(**kwargs)
        self.method = method.lower() if method else 'median'
        self.kernel_size = utils.int_or(kernel_size, 3)
        self.sigma_spatial = utils.float_or(sigma_spatial, 1)

        
    def run(self):
        """Execute the denoise filter."""
        
        ## Setup Output
        dst_ds = self.copy_src_dem()
        if dst_ds is None: return self.src_dem, -1

        with gdalfun.gdal_datasource(self.src_dem) as src_ds:
            if src_ds is None: return self.src_dem, -1
            
            self.init_ds(src_ds)
            
            ## Process chunks to handle large files
            ## (Median/Bilateral need context, so we process with overlap if tiling)
            ## For simplicity in this implementation, we process scanlines or full chunks
            ## utilizing the base class extract/write logic pattern if applicable, 
            ## or just iterate chunks here.
            
            ## Utilizing simple chunking iterator
            ## Note: For Median filters, edge artifacts occur at chunk boundaries unless buffered.
            ## Using a simple overlap buffer logic would be ideal. 
            ## Here we assume file fits in memory or accept slight seam artifacts for simplicity,
            ## but ideally, one uses `utils.yield_srcwin` with a buffer.
            
            buffer = self.kernel_size * 2
            
            for srcwin in utils.yield_srcwin(
                (src_ds.RasterYSize, src_ds.RasterXSize), 
                n_chunk=2048, # Reasonable chunk size
                verbose=self.verbose,
                start_at_edge=False 
            ):
                ## Buffer the source window
                srcwin_buff = utils.buffer_srcwin(
                    srcwin, (src_ds.RasterYSize, src_ds.RasterXSize), buffer
                )
                
                ## Read Data
                band = src_ds.GetRasterBand(self.band)
                src_arr = band.ReadAsArray(
                    srcwin_buff[0], srcwin_buff[1], srcwin_buff[2], srcwin_buff[3]
                )
                
                if src_arr is None: continue
                
                ## Handle NoData (Convert to NaN for processing)
                ndv = self.ds_config['ndv']
                data = src_arr.astype(np.float32)
                data[data == ndv] = np.nan
                
                ## Filter Logic
                result = None
                
                if self.method == 'median':
                    ## Scipy median_filter doesn't handle NaNs well (they propagate).
                    ## We usually fill NaNs before filtering or use a generic_filter (slow).
                    ## Strategy: Fill NaNs with nearest valid value, filter, then restore NaNs.
                    
                    ## Identify valid data
                    valid_mask = ~np.isnan(data)
                    
                    if np.any(valid_mask):
                        ## Simple hole filling for the filter step (nearest neighbor via indices)
                        ## or just use 0 if assumption holds. 
                        ## Better: scipy.ndimage.grey_dilation can fill small holes?
                        ## Here we use a quick fill of mean/0 just to stabilize the filter.
                        filled_data = data.copy()
                        filled_data[np.isnan(filled_data)] = np.nanmean(data)
                        
                        result = scipy.ndimage.median_filter(
                            filled_data, size=self.kernel_size
                        )
                        
                        ## Restore original NoData (don't invent data in large voids)
                        result[~valid_mask] = ndv
                    else:
                        result = data

                elif self.method == 'bilateral':
                    if not HAS_SKIMAGE:
                        utils.echo_warning_msg("scikit-image not installed. Skipping bilateral filter.")
                        result = data
                    else:
                        ## Bilateral needs valid data range.
                        ## Fill NaNs
                        valid_mask = ~np.isnan(data)
                        if np.any(valid_mask):
                            filled_data = data.copy()
                            ## Replace NaN with local mean or 0
                            fill_val = np.nanmean(data)
                            filled_data[~valid_mask] = fill_val
                            
                            ## Apply Filter
                            ## sigma_color: Standard deviation for gray value distance.
                            ## sigma_spatial: Standard deviation for range distance.
                            ## We estimate sigma_color based on data range or roughness? 
                            ## Using default or user param would be best. 
                            ## Adding sigma_color param to class might be needed.
                            
                            try:
                                result = denoise_bilateral(
                                    filled_data, 
                                    sigma_color=None, # Default
                                    sigma_spatial=self.sigma_spatial,
                                    mode='edge'
                                )
                                result[~valid_mask] = ndv
                            except Exception as e:
                                if self.verbose:
                                    utils.echo_warning_msg(f"Bilateral failed: {e}")
                                result = data
                        else:
                            result = data
                
                ## Handling Result
                if result is not None:
                    ## Crop buffer back to original window
                    ## Offset of unbuffered window relative to buffered window
                    x_off = srcwin[0] - srcwin_buff[0]
                    y_off = srcwin[1] - srcwin_buff[1]
                    
                    out_arr = result[y_off:y_off+srcwin[3], x_off:x_off+srcwin[2]]
                    
                    ## Convert NaNs back to NDV
                    out_arr[np.isnan(out_arr)] = ndv
                    
                    ## Write
                    dst_band = dst_ds.GetRasterBand(self.band)
                    dst_band.WriteArray(out_arr, srcwin[0], srcwin[1])

        dst_ds = None
        return self.dst_dem, 0

### End
