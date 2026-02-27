### flats.py - GRId filTerS
##
## Copyright (c) 2024 - 2026 Regents of the University of Colorado
##
## flats.py is part of cudem
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
## Removes flat areas from a DEM by identifying contiguous areas of identical values
## that exceed a specified size threshold.
##
### Code:

import numpy as np
from osgeo import gdal

from cudem import utils
from cudem import gdalfun
from cudem.grits import grits

class Flats(grits.Grits):
    """Remove flat areas from the input DEM.
    Identified flat areas are set to NoData.

    Parameters:
    -----------
    size_threshold : int, optional
        The minimum number of pixels required to define a "flat" area to be removed.
        If None, it is auto-calculated using outlier detection on value counts.
    n_chunk : int, optional
        The processing chunk size in pixels. Defaults to full file if None.
    """

    def __init__(self, size_threshold=None, n_chunk=None, **kwargs):
        super().__init__(**kwargs)
        self.size_threshold = utils.int_or(size_threshold)
        self.n_chunk = utils.int_or(n_chunk)


    def run(self):
        """Execute the flat removal filter."""

        ## Setup Output
        dst_ds = self.copy_src_dem()
        if dst_ds is None: return

        count_removed = 0

        with gdalfun.gdal_datasource(self.src_dem) as src_ds:
            if src_ds is None:
                utils.echo_error_msg('Could not load {self.src_dem}')
                return

            self.init_ds(src_ds)

            ## Default to processing the whole file at once if no chunk size given
            ## (Flats detection works best on larger areas to see full extent of flat)
            if self.n_chunk is None:
                #self.n_chunk = self.ds_config['ny'] # scanline processing often safest
                #self.n_chunk = self.ds_config['nb'] # total pixels

                try:
                    block_size = src_ds.GetRasterBand(1).GetBlockSize()
                    bx, by = block_size[0], block_size[1]

                    if bx > 1 and by > 1:
                        target_size = 4096
                        self.n_chunk = (target_size // bx) * bx
                        if self.n_chunk == 0: self.n_chunk = bx
                    else:
                        self.n_chunk = 4096
                except:
                    ## Fallback
                    self.n_chunk = 4096


            ## Iterate chunks
            for srcwin in gdalfun.gdal_yield_srcwin(
                    src_ds, n_chunk=self.n_chunk, step=self.n_chunk, verbose=self.verbose
            ):
                ## Read Data
                src_arr = self.ds_band.ReadAsArray(srcwin[0], srcwin[1], srcwin[2], srcwin[3]).astype(float)

                ## Identify Unique Values and their Counts
                uv, uv_counts = np.unique(src_arr, return_counts=True)

                ## Determine Threshold
                threshold = self.size_threshold
                if threshold is None:
                    ## Auto-detect: Assume flats are statistical outliers in terms of frequency
                    ## (e.g. background noise is random, flats are highly frequent identical values)
                    outlier_val, _ = self.get_outliers(uv_counts, percentile=99)
                    threshold = outlier_val if not np.isnan(outlier_val) else 100

                ## Identify values that appear too frequently (Flats)
                ## Filter out NoData from being "removed" (it's already removed)
                ndv = self.ds_config['ndv']

                flat_values = uv[uv_counts > threshold]

                ## Create Mask
                ## (Exclude NDV from mask calculation to avoid re-masking it)
                if ndv is not None:
                    flat_values = flat_values[flat_values != ndv]

                if flat_values.size > 0:
                    mask = np.isin(src_arr, flat_values)

                    ## Count and Apply
                    n_removed = np.count_nonzero(mask)
                    if n_removed > 0:
                        count_removed += n_removed
                        src_arr[mask] = ndv

                        ## Write back to Destination
                        dst_band = dst_ds.GetRasterBand(self.band)
                        dst_band.WriteArray(src_arr, srcwin[0], srcwin[1])
                        dst_ds.FlushCache()
        dst_ds = None
        return self.dst_dem, 0

### End
