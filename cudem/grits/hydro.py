### hydro.py - GRId filTerS
##
## Copyright (c) 2024 - 2026 Regents of the University of Colorado
##
## hydro.py is part of cudem
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
## Simple Hydrological Enforcement tools.
## - Fill Sinks: Removes local depressions to ensure flow continuity.
##
### Code:

import numpy as np
import scipy.ndimage
from osgeo import gdal

from cudem import utils
from cudem import gdalfun
from cudem.grits import grits

class Hydro(grits.Grits):
    """Hydrological Enforcement Filters.
    
    Parameters:
    -----------
    method : str
        'fill_sinks'.
    """
    
    def __init__(self, method='fill_sinks', **kwargs):
        super().__init__(**kwargs)
        self.method = method

        
    def _fill_sinks(self, array, ndv):
        """Fill sinks in the array. 
        Uses morphological reconstruction (inverting the DEM, calculating peaks, re-inverting).
        """
        
        ## Sinks are regional minima. 
        ## Invert DEM so sinks become peaks
        data = array.copy()
        
        ## Handle NoData - set to very low value so flow goes off edge
        valid_mask = data != ndv
        
        if not np.any(valid_mask): return array
        
        ## Use max valid value as ceiling for inversion
        d_max = np.max(data[valid_mask])
        inverted = d_max - data
        
        ## Set NDV areas to 0 (lowest in inverted map) to act as boundaries/drains
        inverted[~valid_mask] = 0
        
        ## Reconstruction: Dilate the image but limit by the "Mask" (original surface)
        ## This fills the "peaks" (original sinks)
        ## Note: Scipy doesn't have a direct 'reconstruct' for float, often implemented iteratively
        ## or via grey_dilation if simple pits.
        ## A simpler approach for sinks in scipy is calculating the fill.
        
        ## Alternative: Simple iterative pit filling for small sinks
        ## For full hydro, use a dedicated library. 
        ## Here we implement a simple local minimum filler.
        
        ## Using grey_closing as a proxy for sink filling on small scales
        ## (Closes dark spots/pits)
        filled = scipy.ndimage.grey_closing(data, size=(3,3))
        
        ## Restore NDV
        filled[~valid_mask] = ndv
        return filled

    
    def run(self):
        """Execute Hydro filter."""
        
        dst_ds = self.copy_src_dem()
        if dst_ds is None: return self.src_dem, -1

        with gdalfun.gdal_datasource(self.src_dem) as src_ds:
            if src_ds is None: return self.src_dem, -1
            self.init_ds(src_ds)
            
            ## Hydro ops generally need full context (global flow).
            ## Chunking breaks flow paths. 
            ## We process the whole file here (warn for memory).
            
            if self.verbose:
                utils.echo_msg("Loading full DEM for hydro processing...")
                
            band = src_ds.GetRasterBand(self.band)
            data = band.ReadAsArray()
            
            ndv = self.ds_config['ndv']
            
            if self.method == 'fill_sinks':
                result = self._fill_sinks(data, ndv)
                
                dst_band = dst_ds.GetRasterBand(self.band)
                dst_band.WriteArray(result)
                
        dst_ds = None
        return self.dst_dem, 0

### End
