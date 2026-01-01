### inpaint.py
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
##
## inpaint.py is part of CUDEM
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
## Void filling module using OpenCV Inpainting algorithms.
## Best used for filling small gaps or repairing artifacts.
##
### Code:

import numpy as np
from osgeo import gdal

from cudem import utils
from cudem import gdalfun
from cudem.waffles.waffles import Waffle

## Optional Dependency
try:
    import cv2
    HAS_OPENCV = True
except ImportError:
    HAS_OPENCV = False

class WafflesInpaint(Waffle):
    """VOID FILLING DEM via OpenCV Inpainting.
    
    Fills NoData gaps in a raster using neighboring pixel information.
    This is effectively a computer-vision based "repair" tool.

    Algorithms:
      - 'ns': Navier-Stokes based inpainting (Fluid dynamics).
      - 'telea': Alexandru Telea's Fast Marching Method.

    Parameters:
    -----------
    method (str) : Inpainting algorithm ['ns', 'telea']. Default: 'ns'.
    radius (float) : Radius of a circular neighborhood of each point inpainted 
                     that is considered by the algorithm. Default: 3.0.
    chunk_size (int) : Processing chunk size in pixels.

    < inpaint:method=ns:radius=3 >
    """
    
    def __init__(
            self,
            method='ns',
            radius=3.0,
            chunk_size=None,
            **kwargs
    ):
        super().__init__(**kwargs)
        self.method = method.lower()
        self.radius = utils.float_or(radius, 3.0)
        self.chunk_size = chunk_size
        self.chunk_step = None

        
    def run(self):
        if not HAS_OPENCV:
            utils.echo_error_msg("OpenCV (cv2) must be installed to use the INPAINT module.")
            utils.echo_error_msg("Try: pip install opencv-python-headless")
            return self

        if self.verbose:
            utils.echo_msg(
                f'Inpainting Grid ({self.method}) @ {self.ycount}/{self.xcount} '
                f'with radius {self.radius}...'
            )

        if self.chunk_size is None:
            n_chunk = 1024 # Inpainting is fast, can handle larger chunks
        else:
            n_chunk = self.chunk_size

        ## Determine OpenCV flag
        if self.method == 'telea':
            inpaint_mode = cv2.INPAINT_TELEA
        else:
            inpaint_mode = cv2.INPAINT_NS

        ## Initialize Output Dataset
        ## We copy the source structure since we are essentially modifying it
        with gdalfun.gdal_datasource(self.stack) as stack_ds:
            stack_gt = stack_ds.GetGeoTransform()
            stack_proj = stack_ds.GetProjection()
            
            ## Prepare Output
            driver = gdal.GetDriverByName(self.fmt)
            dst_ds = driver.Create(
                self.fn,
                self.xcount,
                self.ycount,
                1,
                gdal.GDT_Float32,
                options=self.co
            )
            
            dst_ds.SetGeoTransform(self.dst_gt) # Use requested GT, not stack GT
            dst_ds.SetProjection(gdalfun.osr_wkt(self.dst_srs))
            
            elev_band = dst_ds.GetRasterBand(1)
            elev_band.SetNoDataValue(self.ndv)
            elev_band.SetDescription(f"Inpainted ({self.method})")

            ## Process in Chunks
            ## Note: Inpainting at chunk boundaries can create artifacts.
            ## Ideally, one would use a buffer, but for simplicity here we process standard chunks.
            ## Users should avoid extremely small chunks.
            
            band = stack_ds.GetRasterBand(1)
            nodata = band.GetNoDataValue()
            if nodata is None: nodata = -9999

            for srcwin in utils.yield_srcwin(
                    (self.ycount, self.xcount),
                    n_chunk=n_chunk,
                    msg='Inpainting',
                    verbose=self.verbose
            ):
                x_off, y_off, x_size, y_size = srcwin
                
                ## Read Source Data
                ## Note: We are reading from the 'stack', which acts as the source "image"
                ## If the stack resolution != output resolution, this naive read might mismatch.
                ## However, waffles generally ensures stack is built at target resolution (xinc/yinc).
                try:
                    data_chunk = band.ReadAsArray(x_off, y_off, x_size, y_size).astype(np.float32)
                except Exception as e:
                    ## Handle edge cases where stack might be slightly different size or alignment
                    ## Fallback: fill with NDV
                    utils.echo_warning_msg(f"Read error @ {srcwin}: {e}")
                    continue

                ## Identify Holes (NoData)
                ## Create mask: 0 = Valid Data, 1 = Hole to Fill
                ## OpenCV expects uint8 mask where non-zero pixels indicate the area to be inpainted.
                if np.isnan(nodata):
                    mask = np.isnan(data_chunk).astype(np.uint8)
                else:
                    mask = (data_chunk == nodata).astype(np.uint8)

                ## Check if there is anything to inpaint
                if np.count_nonzero(mask) == 0:
                    ## No holes, just write original data
                    elev_band.WriteArray(data_chunk, x_off, y_off)
                    continue
                
                if np.count_nonzero(mask) == mask.size:
                    ## All holes (empty chunk), write NDV
                    elev_band.WriteArray(np.full(data_chunk.shape, self.ndv), x_off, y_off)
                    continue

                ## Prepare Data for OpenCV
                ## Replace NDV with 0 (or mean) temporarily for the algo input, 
                ## though the mask tells it what to ignore/fill.
                ## data_chunk[mask == 1] = 0 

                try:
                    ## Execute Inpainting
                    ## src: Input image (8-bit, 16-bit unsigned or 32-bit float 1-channel)
                    ## inpaintMask: 8-bit 1-channel mask. Non-zero pixels indicate the area that needs to be inpainted.
                    inpainted_chunk = cv2.inpaint(
                        data_chunk, 
                        mask, 
                        self.radius, 
                        inpaint_mode
                    )
                    
                    ## Write result
                    elev_band.WriteArray(inpainted_chunk, x_off, y_off)
                    
                except Exception as e:
                    utils.echo_warning_msg(f"Inpaint failed @ {srcwin}: {e}")
                    ## Fallback to original
                    elev_band.WriteArray(data_chunk, x_off, y_off)

            dst_ds = None
            
        return self

### End
