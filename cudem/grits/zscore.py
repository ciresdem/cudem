### zscore.py - GRId filTerS
##
## Copyright (c) 2024 - 2025 Regents of the University of Colorado
##
## zscore.py is part of cudem
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
## Filters data based on local Z-Score (Standard Score).
## Z = (Value - LocalMean) / LocalStdDev
##
### Code:

import numpy as np
import scipy.ndimage
from osgeo import gdal

from cudem import utils
from cudem import gdalfun
from cudem.grits import grits

class ZScore(grits.Grits):
    """Local Z-Score Filter.

    Parameters:
    -----------
    threshold : float
        The Z-Score threshold. Pixels with |Z| > threshold are masked.
        Typical values: 2.0 (95%), 3.0 (99.7%).
    kernel_size : int
        Size of the neighborhood window (pixels).
    """

    def __init__(self, threshold=3.0, kernel_size=5, **kwargs):
        super().__init__(**kwargs)
        self.threshold = utils.float_or(threshold, 3.0)
        self.kernel_size = utils.int_or(kernel_size, 5)


    def run(self):
        """Execute the Z-Score filter."""

        dst_ds = self.copy_src_dem()
        if dst_ds is None: return self.src_dem, -1

        with gdalfun.gdal_datasource(self.src_dem) as src_ds:
            if src_ds is None: return self.src_dem, -1
            self.init_ds(src_ds)

            ## Use chunks with buffer for window stats
            buffer = self.kernel_size

            for srcwin in utils.yield_srcwin(
                (src_ds.RasterYSize, src_ds.RasterXSize),
                n_chunk=2048,
                verbose=self.verbose
            ):
                srcwin_buff = utils.buffer_srcwin(
                    srcwin, (src_ds.RasterYSize, src_ds.RasterXSize), buffer
                )

                band = src_ds.GetRasterBand(self.band)
                src_arr = band.ReadAsArray(*srcwin_buff)

                if src_arr is None: continue

                ndv = self.ds_config['ndv']
                data = src_arr.astype(np.float32)
                data[data == ndv] = np.nan

                ## Calculate Local Mean and StdDev
                ## Uniform filter gives the mean
                ## StdDev = sqrt(E[x^2] - (E[x])^2)

                ## Fill NaNs for convolution stability
                ## (Simple 0-fill, or mean-fill would be better but slower)
                filled_data = data.copy()
                filled_data[np.isnan(data)] = np.nanmean(data)

                local_mean = scipy.ndimage.uniform_filter(
                    filled_data, size=self.kernel_size, mode='reflect'
                )

                local_sq_mean = scipy.ndimage.uniform_filter(
                    filled_data**2, size=self.kernel_size, mode='reflect'
                )

                local_var = local_sq_mean - local_mean**2
                local_std = np.sqrt(np.maximum(0, local_var)) # Ensure non-negative

                ## Avoid div by zero
                local_std[local_std == 0] = 1e-6

                ## Calculate Z
                z_score = np.abs((filled_data - local_mean) / local_std)

                ## Crop Buffer
                x_off = srcwin[0] - srcwin_buff[0]
                y_off = srcwin[1] - srcwin_buff[1]

                z_crop = z_score[y_off:y_off+srcwin[3], x_off:x_off+srcwin[2]]
                data_crop = data[y_off:y_off+srcwin[3], x_off:x_off+srcwin[2]]

                ## Identify Outliers
                ## Mask pixels where Z > Threshold AND pixel was originally valid
                mask = (z_crop > self.threshold) & (~np.isnan(data_crop))

                if np.any(mask):
                    dst_band = dst_ds.GetRasterBand(self.band)
                    dst_data = dst_band.ReadAsArray(*srcwin)
                    dst_data[mask] = ndv
                    dst_band.WriteArray(dst_data, srcwin[0], srcwin[1])

        dst_ds = None
        return self.dst_dem, 0

### End
