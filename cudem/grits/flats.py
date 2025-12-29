### flats.py - GRId filTerS
##
## Copyright (c) 2024 - 2025 Regents of the University of Colorado
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
###############################################################################
### Commentary:
##
##
### Code:

import numpy as np

from cudem import utils
from cudem import gdalfun
from cudem.grits import grits

class Flats(grits.Grits):
    """Remove flat areas from the input DEM

    Parameters:

    size_threshold(int) - the minimum flat area in pixels to remove
    n_chunk(int) - the moving window size in pixels
    """
    
    def __init__(self, size_threshold=None, n_chunk=None, **kwargs):
        super().__init__(**kwargs)
        self.size_threshold = utils.int_or(size_threshold)
        self.n_chunk = utils.int_or(n_chunk)

        
    def run(self):
        """Discover and remove flat zones"""
        
        count = 0
        dst_ds = self.copy_src_dem()
        with gdalfun.gdal_datasource(self.src_dem) as src_ds:
            if src_ds is not None:
                self.init_ds(src_ds)
                if self.n_chunk is None:
                    self.n_chunk = self.ds_config['nb']

                for srcwin in gdalfun.gdal_yield_srcwin(
                        src_ds, n_chunk=self.n_chunk,
                        step=self.n_chunk, verbose=True
                ):
                    src_arr = self.ds_band.ReadAsArray(*srcwin).astype(float)
                    uv, uv_counts = np.unique(src_arr, return_counts=True)
                    if self.size_threshold is None:
                        _size_threshold = self.get_outliers(uv_counts, 99)[0]
                    else:
                        _size_threshold = self.size_threshold

                    uv_ = uv[uv_counts > _size_threshold]
                    mask = np.isin(src_arr, uv_)
                    count += np.count_nonzero(mask)
                    src_arr[mask] = self.ds_config['ndv']
                    dst_band = dst_ds.GetRasterBand(self.band)
                    dst_band.WriteArray(src_arr, srcwin[0], srcwin[1])
                
        dst_ds = None
        if self.verbose:
            utils.echo_msg(f'removed {count} flats.')            
        return self.dst_dem, 0

### End
