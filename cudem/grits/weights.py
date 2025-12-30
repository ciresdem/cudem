### weights.py - GRId filTerS
##
## Copyright (c) 2024 - 2025 Regents of the University of Colorado
##
## weights.py is part of cudem
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
## Filter that buffers around areas of high weight (quality) and removes lower-quality data 
## that falls within that buffer zone. Useful for cleaning up "halos" around good data.
##
### Code:

import numpy as np
import scipy.ndimage
from osgeo import gdal

from cudem import utils
from cudem import gdalfun
from cudem.grits import grits

class Weights(grits.Grits):
    """Create a nodata buffer around data at and above `weight_threshold`.
    Removes data from the source DEM that falls within this buffer but has a lower weight.

    Parameters:
    -----------
    weight_thresholds : str or list
        The weight threshold(s); separate levels with slash '/'.
    buffer_sizes : str or list
        The number of cells to buffer for each level; separate levels with slash '/'.
    gap_fill_sizes : str or list
        The number of cells to gap-fill (close holes) before buffering; separate levels with slash '/'.
    binary_dilation : bool
        Use binary dilation (faster) instead of distance transform expansion.
    fill_holes : bool
        Fill holes in the weight level mask before processing.
    """
    
    def __init__(self, weight_threshold=None, buffer_cells=1, gap_fill_cells=None,
                 weight_thresholds=None, buffer_sizes=None, gap_fill_sizes=None,
                 binary_dilation=True, binary_pulse=False,
                 fill_holes=False, **kwargs):
        super().__init__(**kwargs)
        
        ## Base defaults
        self.base_weight_threshold = utils.float_or(weight_threshold, 1)
        self.base_buffer_cells = utils.int_or(buffer_cells, 1)
        self.base_gap_fill_cells = utils.int_or(gap_fill_cells, 0)
        
        self.weight_thresholds_str = weight_thresholds
        self.buffer_sizes_str = buffer_sizes
        self.gap_fill_sizes_str = gap_fill_sizes
        
        self.binary_dilation = binary_dilation
        self.binary_pulse = binary_pulse
        self.fill_holes = fill_holes
        
        self._init_thresholds()

        
    def _init_thresholds(self):
        """Parse and align threshold/buffer lists."""
        
        if self.weight_thresholds_str is not None:
            self.weight_thresholds = [float(x) for x in str(self.weight_thresholds_str).split('/')]

            if self.buffer_sizes_str is not None:
                self.buffer_sizes = [int(x) for x in str(self.buffer_sizes_str).split('/')]
            else:
                self.buffer_sizes = [self.base_buffer_cells] * len(self.weight_thresholds)

            if self.gap_fill_sizes_str is not None:
                self.gap_fill_sizes = [int(x) for x in str(self.gap_fill_sizes_str).split('/')]
            else:
                self.gap_fill_sizes = [self.base_gap_fill_cells] * len(self.weight_thresholds)
                
        else:
            self.weight_thresholds = [self.base_weight_threshold]
            self.buffer_sizes = [self.base_buffer_cells]
            self.gap_fill_sizes = [self.base_gap_fill_cells]

        ## Ensure lists match length of thresholds
        len_diff = len(self.weight_thresholds) - len(self.buffer_sizes)
        if len_diff > 0:
            self.buffer_sizes.extend([1] * len_diff)
            
        len_diff = len(self.weight_thresholds) - len(self.gap_fill_sizes)
        if len_diff > 0:
            self.gap_fill_sizes.extend([0] * len_diff)

        if self.verbose:
            utils.echo_msg(f"Weights Config: Thresh={self.weight_thresholds}, Buff={self.buffer_sizes}, Gap={self.gap_fill_sizes}")

            
    def binary_closed_dilation(self, arr, iterations=1, closing_iterations=1):
        """Perform binary closing followed by dilation."""
        
        ## Note: Original implementation did dilation THEN erosion for 'closing', 
        ## which is technically a 'Close' operation (dilate->erode).
        
        ## Structure for connectivity (8-connected usually preferred for rasters)
        struct = scipy.ndimage.generate_binary_structure(2, 2)
        
        ## Closing (fill gaps)
        temp = arr.copy()
        if closing_iterations > 0:
            temp = scipy.ndimage.binary_dilation(temp, iterations=closing_iterations, structure=struct)
            temp = scipy.ndimage.binary_erosion(temp, iterations=closing_iterations, border_value=1, structure=struct)

        ## Dilation (expand buffer)
        if iterations > 0:
            return scipy.ndimage.binary_dilation(temp, iterations=iterations, structure=struct)
        return temp

    
    def binary_reversion(self, arr, iterations, closing_iterations):
        """Alternative expansion logic."""
        
        reversion_arr = arr.copy()
        total_dilation = closing_iterations + iterations
        
        struct = scipy.ndimage.generate_binary_structure(2, 2)
        
        ## Dilate out
        if total_dilation > 0:
            reversion_arr = scipy.ndimage.binary_dilation(reversion_arr, iterations=total_dilation, structure=struct)
            
        ## Erode back (Revert)
        erosion_iterations = max(0, closing_iterations - iterations) # Logic from original seems to be (closing - buffer)?
        ## Actually logic was: max(closing_iterations-iterations, 0). 
        ## If closing=10, buffer=4 -> erode 6. Net effect: +4 buffer, but holes < 20 filled?
        
        if erosion_iterations > 0:
            reversion_arr = scipy.ndimage.binary_erosion(reversion_arr, iterations=erosion_iterations, border_value=1, structure=struct)

        return reversion_arr

    
    def run(self):
        """Execute the weights filter."""
        
        ## Check for Weights
        if self.weight_mask is None:
            utils.echo_error_msg("No weight mask provided for Weights filter.")
            return self.src_dem, -1
        
        ## Setup Output
        dst_ds = self.copy_src_dem()
        if dst_ds is None: return self.src_dem, -1
        
        ## Load Data
        with gdalfun.gdal_datasource(self.src_dem) as src_ds:
            if src_ds is None: return self.src_dem, -1
            
            self.init_ds(src_ds)
            
            ## Load Weight Array using base class helper
            w_arr = self._get_mask_array('weight_mask', 'weight_is_fn', 'weight_is_band')
            
            if w_arr is None:
                utils.echo_error_msg("Failed to load weight array.")
                return self.src_dem, -1
                
            ## Handle NDV in weights
            ## Assuming NDV in weights means "No Data", effectively 0 weight
            w_arr[np.isnan(w_arr)] = -9999
            
            ## Processing Loop
            with utils.ccp(
                    total=len(self.weight_thresholds),
                    desc=f'Buffering around weights {self.weight_thresholds}',
                    leave=self.verbose
            ) as w_pbar:
                
                for n, weight_threshold in enumerate(self.weight_thresholds):
                    ## Identify High Quality Pixels
                    ## (Pixels >= threshold)
                    high_quality_mask = w_arr >= weight_threshold
                    
                    ## Expand this mask to create the buffer zone
                    if self.binary_dilation:                            
                        if self.binary_pulse:
                            expanded_mask = self.binary_closed_dilation(
                                high_quality_mask, 
                                iterations=self.buffer_sizes[n],
                                closing_iterations=self.gap_fill_sizes[n]
                            )                            
                        else:
                            expanded_mask = self.binary_reversion(
                                high_quality_mask, 
                                self.buffer_sizes[n],
                                self.gap_fill_sizes[n]
                            )
                    else:
                        ## Euclidean/legacy expansion fallback
                        ## (Simplified here to use scipy distance transform for robustness if binary_dilation is False)
                        dt = scipy.ndimage.distance_transform_edt(~high_quality_mask)
                        expanded_mask = dt <= self.buffer_sizes[n]

                    if self.fill_holes:
                        expanded_mask = scipy.ndimage.binary_fill_holes(expanded_mask)
                        
                    ## Identify pixels to REMOVE
                    ## Condition: Inside Expanded Buffer AND Weight < Threshold
                    ## i.e., "Bad data near Good data"
                    remove_mask = expanded_mask & (w_arr < weight_threshold)
                    
                    if np.count_nonzero(remove_mask) > 0:
                        ## Apply to Output
                        for b in range(1, dst_ds.RasterCount + 1):
                            band = dst_ds.GetRasterBand(b)
                            data = band.ReadAsArray()
                            
                            ## Set to NDV
                            ## Note: self.ds_config['ndv'] might be None/float
                            ndv = self.ds_config['ndv']
                            data[remove_mask] = ndv
                            
                            band.WriteArray(data)
                            band.FlushCache()
                            
                        #w_arr[remove_mask] = np.nan
                        w_arr[remove_mask] = -9999 # Set to low/nodata

                    w_pbar.update()

        dst_ds = None        
        return self.dst_dem, 0

    
class WeightZones(Weights):
    """Experimental: Filter based on contiguous zone size and density.
    Removes disconnected islands of data that are too small or sparse.
    """
    
    def __init__(self, size_threshold=None, **kwargs):
        super().__init__(**kwargs)
        self.size_threshold = utils.float_or(size_threshold)

        
    def run(self):
        if self.weight_mask is None: return self.src_dem, -1

        dst_ds = self.copy_src_dem()
        if dst_ds is None: return self.src_dem, -1

        with gdalfun.gdal_datasource(self.src_dem) as src_ds:
            if src_ds is None: return self.src_dem, -1
            
            self.init_ds(src_ds)
            w_arr = self._get_mask_array('weight_mask', 'weight_is_fn', 'weight_is_band')
            
            if w_arr is None: return self.src_dem, -1
            
            # Treat NDV
            # w_arr[np.isnan(w_arr)] = 0
            
            # Calculate Density
            n_den = self._density(w_arr) # Base class method
            
            # Determine Size Threshold
            # Heuristic: 25000 units area scaled by density?
            cell_size_x = self.ds_config['geoT'][1]
            # Assuming meters for heuristic default
            # If degrees, convert approx:
            if cell_size_x < 0.1: cell_size_x *= 111120 

            if self.size_threshold is None:
                m_size = 25000
                size_threshold = (m_size * (1/n_den if n_den > 0 else 1)) / cell_size_x
            else:
                size_threshold = self.size_threshold
            
            if self.verbose:
                utils.echo_msg(f"WeightZones: Density={n_den:.4f}, SizeThreshold={size_threshold}")

            # Binarize based on weight threshold
            binary_mask = (w_arr >= self.base_weight_threshold).astype(int)
            
            # Label contiguous regions
            labeled_arr, num_features = scipy.ndimage.label(binary_mask)
            
            # Sum size of each label
            # (bincount is faster than sum_labels for simple count)
            region_sizes = np.bincount(labeled_arr.ravel())
            
            # Identify small regions (Label 0 is background, ignore it usually)
            # Mask: True where region is small
            # Note: region_sizes indices correspond to label IDs
            small_regions_mask = region_sizes < size_threshold
            small_regions_mask[0] = False # Don't mask the background itself (value 0)
            
            # Create a removal mask
            # map the boolean mask back to the array shape
            remove_mask = small_regions_mask[labeled_arr]
            
            count_removed = np.count_nonzero(remove_mask)
            
            if count_removed > 0:
                if self.verbose:
                    utils.echo_msg(f"Removing {count_removed} pixels in small zones.")
                    
                for b in range(1, dst_ds.RasterCount + 1):
                    band = dst_ds.GetRasterBand(b)
                    data = band.ReadAsArray()
                    data[remove_mask] = self.ds_config['ndv']
                    band.WriteArray(data)
                    band.FlushCache()

        dst_ds = None        
        return self.dst_dem, 0

### End
