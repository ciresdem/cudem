### morphology.py - GRId filTerS
##
## Copyright (c) 2024 - 2025 Regents of the University of Colorado
##
## morphology.py is part of cudem
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
## Morphological filters for DEMs (Erosion, Dilation, Opening, Closing).
## Useful for removing small features (trees, noise) or filling sinks.
##
### Code:

import numpy as np
import scipy.ndimage
from osgeo import gdal

from cudem import utils
from cudem import gdalfun
from cudem.grits import grits

class Morphology(grits.Grits):
    """Apply morphological operations to the DEM.

    Operations:
    - Erosion: Replaces pixel with minimum in neighborhood (widens valleys, removes peaks).
    - Dilation: Replaces pixel with maximum in neighborhood (widens peaks, fills pits).
    - Opening: Erosion followed by Dilation (removes small bright spots/peaks).
    - Closing: Dilation followed by Erosion (fills small dark spots/pits).

    Parameters:
    -----------
    operation : str
        'erosion', 'dilation', 'opening', 'closing'. Default 'erosion'.
    kernel_size : int
        Size of the structuring element (pixels). Default 3.
    """
    
    def __init__(self, operation='erosion', kernel_size=3, **kwargs):
        super().__init__(**kwargs)
        self.operation = operation.lower()
        self.kernel_size = utils.int_or(kernel_size, 3)

        
    def run(self):
        """Execute the morphological filter."""
        
        dst_ds = self.copy_src_dem()
        if dst_ds is None: return self.src_dem, -1

        with gdalfun.gdal_datasource(self.src_dem) as src_ds:
            if src_ds is None: return self.src_dem, -1
            
            self.init_ds(src_ds)
            
            ## Setup structuring element (square/rect kernel)
            ## Could add option for circular/diamond later if needed
            structure = np.ones((self.kernel_size, self.kernel_size))
            
            ## Buffer size needs to cover the kernel radius
            buffer = self.kernel_size
            
            for srcwin in utils.yield_srcwin(
                (src_ds.RasterYSize, src_ds.RasterXSize), 
                n_chunk=2048,
                verbose=self.verbose
            ):
                ## Buffer window
                srcwin_buff = utils.buffer_srcwin(
                    srcwin, (src_ds.RasterYSize, src_ds.RasterXSize), buffer
                )
                
                band = src_ds.GetRasterBand(self.band)
                src_arr = band.ReadAsArray(
                    srcwin_buff[0], srcwin_buff[1], srcwin_buff[2], srcwin_buff[3]
                )
                
                if src_arr is None: continue
                
                ndv = self.ds_config['ndv']
                data = src_arr.astype(np.float32)
                
                ## Identify NoData
                ## Morphology operations propagate min/max values.
                ## If we have NDV, we must handle it carefully.
                ## For Erosion (Min): NDV (often very low or high) can dominate.
                ## Strategy: Mask NDV, process, re-mask.
                
                valid_mask = data != ndv
                
                ## Fill invalid data with neutral values so they don't affect min/max
                ## Erosion (Min): Fill with Max possible
                ## Dilation (Max): Fill with Min possible
                ## For floats, use inf/-inf or max/min of data
                
                data_min = np.nanmin(data[valid_mask]) if np.any(valid_mask) else 0
                data_max = np.nanmax(data[valid_mask]) if np.any(valid_mask) else 0
                
                fill_val = data_max if self.operation in ['erosion', 'opening'] else data_min
                
                data[~valid_mask] = fill_val
                
                ## Execute Operation
                result = None
                
                try:
                    if self.operation == 'erosion':
                        result = scipy.ndimage.grey_erosion(data, structure=structure)
                    elif self.operation == 'dilation':
                        result = scipy.ndimage.grey_dilation(data, structure=structure)
                    elif self.operation == 'opening':
                        result = scipy.ndimage.grey_opening(data, structure=structure)
                    elif self.operation == 'closing':
                        result = scipy.ndimage.grey_closing(data, structure=structure)
                    else:
                        utils.echo_warning_msg(f"Unknown operation: {self.operation}")
                        result = data
                except Exception as e:
                    if self.verbose:
                        utils.echo_error_msg(f"Morphology failed: {e}")
                    result = data

                ## Restore NoData
                ## Strict: If original was NDV, keep it NDV
                ## (Optional: Morphology can technically "fill" holes, 
                ## but standard behavior preserves valid data extents unless specified otherwise)
                result[~valid_mask] = ndv
                
                # Write Result (Crop buffer)
                x_off = srcwin[0] - srcwin_buff[0]
                y_off = srcwin[1] - srcwin_buff[1]
                
                out_arr = result[y_off:y_off+srcwin[3], x_off:x_off+srcwin[2]]
                
                dst_band = dst_ds.GetRasterBand(self.band)
                dst_band.WriteArray(out_arr, srcwin[0], srcwin[1])

        dst_ds = None
        return self.dst_dem, 0

### End
