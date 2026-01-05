### rivers.py
##
## Copyright (c) 2023 - 2026 Regents of the University of Colorado
##
## rivers.py is part of CUDEM
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
##
### Commentary:
##
## Grits module to locate potential rivers/water bodies by analyzing 
## NoData voids, calculating their circularity, and flatness of banks.
##
### Code:

import os
import sys
import numpy as np
import scipy.ndimage
from osgeo import ogr
from cudem import utils
from cudem import gdalfun
from cudem.grits import grits

class Rivers(grits.Grits):
    """
    Locate possible rivers in a DEM based on void geometry and border statistics.
    
    Analysis Steps:
    1. Identifies NoData voids in the source DEM.
    2. Filters small voids based on `size_threshold`.
    3. Calculates the Standard Deviation of the pixels immediately bordering the void.
       (Rivers typically have flat/low-std-dev banks).
    4. Polygonizes the voids.
    5. Calculates Area, Perimeter, and Circularity Index.
    6. Filters polygons that are too circular (keeping thread-like features).
    
    Args:
        size_threshold (int): Minimum size (in pixels) of a void to be considered.
        max_circularity (float): Maximum circularity index (0-1). 
                                 1.0 is a circle, ~0.0 is a line.
                                 Rivers should be low.
    """

    def __init__(self, size_threshold=100, max_circularity=0.25, **kwargs):
        super().__init__(**kwargs)
        self.size_threshold = utils.int_or(size_threshold, 100)
        self.max_circularity = utils.float_or(max_circularity, 0.25)

        
    def run(self):
        """Run the river finder algorithm."""
        
        if self.src_dem is None:
            return

        ## Setup Output
        if self.dst_dem is None:
            self.dst_dem = utils.append_fn(self.src_dem, '_rivers')
            
        utils.copy_fn(self.src_dem, self.dst_dem)

        ## Load Data
        with gdalfun.gdal_datasource(self.dst_dem, update=True) as src_ds:
            src_band = src_ds.GetRasterBand(1)
            src_arr = src_band.ReadAsArray()
            src_config = gdalfun.gdal_infos(src_ds)
            
            ## Mask NoData
            ndv = src_config['ndv']
            if ndv is not None:
                src_arr[src_arr == ndv] = np.nan
            
            ## Generate binary mask of voids (NaNs)
            msk_arr = np.zeros(src_arr.shape)
            msk_arr[np.isnan(src_arr)] = 1
            
            ## Label connected components (voids)
            labeled_arr, num_features = scipy.ndimage.label(msk_arr)
            
            ## Get sizes of features
            feature_sizes = scipy.ndimage.sum(msk_arr, labeled_arr, index=range(1, num_features + 1))
            
            ## Reset mask to store output values
            msk_arr[:] = 0
            
            ## Structuring element for dilation (3x3 square, equivalent to expand_for shift=1)
            dilation_structure = np.ones((3,3))

            ## Iterate and Analyze Voids
            ## Using utils.ccp for progress bar style consistency
            with utils.ccp(total=num_features, desc='Scanning voids', leave=self.verbose) as pbar:
                for i in range(num_features):
                    pbar.update(1)
                    
                    ## Skip small voids
                    if feature_sizes[i] < self.size_threshold:
                        continue
                        
                    ## Create boolean mask for this specific void
                    label_id = i + 1
                    void_mask = (labeled_arr == label_id)
                    
                    ## Dilation: Expand mask to include boundary pixels
                    expanded_mask = scipy.ndimage.binary_dilation(void_mask, structure=dilation_structure)
                    
                    ## Calculate StdDev of the BANK pixels.
                    ## src_arr has NaNs in the void center.
                    ## src_arr[expanded_mask] includes the void (NaNs) and the immediate valid neighbors.
                    flat_value = np.nanstd(src_arr[expanded_mask])
                    
                    ## Assign this roughness value to the void mask
                    msk_arr[void_mask] = np.ceil(flat_value)

            ## Restore original NaNs where appropriate.
            src_arr[np.isnan(src_arr)] = src_config['ndv']
            
            ## Write the mask array (containing std_dev values of river banks)
            src_band.SetNoDataValue(0)
            src_band.WriteArray(msk_arr)
            src_band.FlushCache()
            
        with gdalfun.gdal_datasource(self.dst_dem) as nd_ds:
            dst_layer, ogr_format = gdalfun.ogr_polygonize(nd_ds)
            
            ## Define vector filename
            dst_vector = f"{dst_layer}.{gdalfun.ogr_fext(ogr_format)}"
            
            if self.verbose:
                utils.echo_msg(f"Vectorizing to {dst_vector} for geometric analysis...")

            ## OGR SQL Geometry Calculations
            ## Calculate Area
            utils.run_cmd(f'ogrinfo {dst_vector} -sql "ALTER TABLE {dst_layer} ADD COLUMN Area_ float"')
            utils.run_cmd(f'ogrinfo {dst_vector} -dialect SQLite -sql "UPDATE {dst_layer} SET Area_ = ST_Area(geometry)"')
            
            ## Calculate Perimeter
            utils.run_cmd(f'ogrinfo {dst_vector} -sql "ALTER TABLE {dst_layer} ADD COLUMN Perim_ float"')
            utils.run_cmd(f'ogrinfo {dst_vector} -dialect SQLite -sql "UPDATE {dst_layer} SET Perim_ = ST_Perimeter(geometry)"')
            
            ## Calculate Circularity
            ## Formula: (4 * PI * Area) / (Perimeter^2)
            utils.run_cmd(f'ogrinfo {dst_vector} -sql "ALTER TABLE {dst_layer} ADD COLUMN Circul_ float"')
            utils.run_cmd(f'ogrinfo {dst_vector} -dialect SQLite -sql "UPDATE {dst_layer} SET Circul_ = Circularity(geometry)"')

            ## Filter Vectors
            d_ds = ogr.Open(dst_vector, 1) # 1 = Update
            if d_ds:
                d_layer = d_ds.GetLayer()
                delete_ids = []
                
                for feat in d_layer:
                    ## Filter out features that are too circular (lakes/ponds)
                    ## We want rivers (low circularity)
                    c_val = feat.GetField('Circul_')
                    if c_val is not None and c_val > self.max_circularity:
                        delete_ids.append(feat.GetFID())
                
                ## Batch delete
                if len(delete_ids) > 0:
                    if self.verbose:
                        utils.echo_msg(f"Removing {len(delete_ids)} non-river features...")
                    
                    for fid in delete_ids:
                        d_layer.DeleteFeature(fid)
                        
                d_ds = None

        return self.dst_dem

### End
