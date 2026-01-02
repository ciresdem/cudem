### blur.py - GRId filTerS
##
## Copyright (c) 2024 - 2026 Regents of the University of Colorado
##
## blur.py is part of cudem
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
## Gaussian Blur Filter
##
### Code:

import numpy as np
from scipy.signal import fftconvolve
from osgeo import gdal

from cudem import utils
from cudem import gdalfun
from cudem.grits import grits

class Blur(grits.Grits):
    """Blur DEM values using a Gaussian Blur.

    Parameters:
    -----------
    blur_factor : float
        The sigma (standard deviation) for the Gaussian kernel.
        Controls the amount of blurring. Default is 1.
    """
    
    def __init__(self, blur_factor: float = 1, **kwargs: any):
        super().__init__(**kwargs)
        self.blur_factor = utils.int_or(blur_factor, 1)

        
    def np_gaussian_blur(self, in_array: np.ndarray, size: int) -> np.ndarray:
        """Blur an array using fftconvolve from scipy.signal.
        
        Args:
            in_array (np.ndarray): Input 2D array.
            size (int): The blurring scale-factor (kernel radius/sigma).
            
        Returns:
            np.ndarray: The blurred array.
        """
        
        ## Pad array to handle edges
        padded_array = np.pad(in_array, size, 'symmetric')
        
        ## Create Gaussian Kernel
        x, y = np.mgrid[-size:size + 1, -size:size + 1]
        g = np.exp(-(x**2 / float(size) + y**2 / float(size)))
        g = (g / g.sum()).astype(in_array.dtype)
        
        ## Convolve
        out_array = fftconvolve(padded_array, g, mode='valid')
        
        return out_array

    
    def run(self):
        """Run the Gaussian Blur filter on the source DEM.
        """
        
        ## Create Output Dataset (Copy of Source)
        dst_ds = self.copy_src_dem()
        if dst_ds is None: 
            return

        with gdalfun.gdal_datasource(self.src_dem) as src_ds:
            if src_ds is None: return
            
            self.init_ds(src_ds)
            
            ## Read Source Array
            src_band = src_ds.GetRasterBand(self.band)
            ds_array = src_band.ReadAsArray()
            
            ## Handle NoData
            ndv = self.ds_config['ndv']
            
            ## Create a validity mask (1=Valid, NaN=Invalid)
            ## We convert NDV to NaN for processing logic
            float_array = ds_array.astype(np.float32)
            float_array[float_array == ndv] = np.nan
            
            valid_mask = ~np.isnan(float_array)
            
            ## Fill NaNs with 0 (or nearest neighbor?) for FFT to avoid artifact spread
            ## Simple 0-fill can pull down edges. 
            ## 'symmetric' padding handles image edges, but internal holes (lakes) are tricky.
            ## A common strategy is to fill holes, blur, then re-mask.
            
            ## Here we just zero out NaNs for the convolution input
            fft_input = float_array.copy()
            fft_input[~valid_mask] = 0
            
            ## Perform Blur
            blurred_array = self.np_gaussian_blur(fft_input, self.blur_factor)
            
            ## Restore Original Mask (don't fill in holes with blur data)
            ## The blur spreads value into NoData areas; we usually want to clip that back.
            blurred_array[~valid_mask] = ndv
            
            ## Write Result
            dst_band = dst_ds.GetRasterBand(self.band)
            dst_band.WriteArray(blurred_array)
            
        dst_ds = None

        
### End
