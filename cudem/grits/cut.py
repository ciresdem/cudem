### cut.py - GRId filTerS
##
## Copyright (c) 2024 - 2026 Regents of the University of Colorado
##
## cut.py is part of cudem
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
## Cuts (masks) a DEM to a specific region. 
## Cells outside the region are set to NoData.
## Use invert=True to mask cells INSIDE the region instead.
##
### Code:

import numpy as np
from osgeo import gdal

from cudem import utils
from cudem import regions
from cudem import gdalfun
from cudem.grits import grits

class Cut(grits.Grits):
    """Masks the input DEM to a specific region.
    Pixels outside this region are set to NoData.

    Parameters:
    -----------
    cut_region : str or list or Region
        The region to cut to (xmin/xmax/ymin/ymax).
    invert : bool
        If True, masks pixels INSIDE the region instead (hole punch).
    """
    
    def __init__(self, cut_region=None, invert=False, **kwargs):
        super().__init__(**kwargs)
        self.cut_region = cut_region
        self.invert = invert

        
    def run(self):
        """Execute the cut filter."""
        
        if self.cut_region is None:
            ## If no region provided, do nothing or warn
            utils.echo_warning_msg("No cut region provided.")
            return self.src_dem, 0

        ## Setup Output
        dst_ds = self.copy_src_dem()
        if dst_ds is None: return self.src_dem, -1

        ## Parse Cut Region
        region_obj = None
        if isinstance(self.cut_region, regions.Region):
            region_obj = self.cut_region
        elif isinstance(self.cut_region, (list, tuple)):
            region_obj = regions.Region().from_list(self.cut_region)
        elif isinstance(self.cut_region, str):
            region_obj = regions.Region().from_string(self.cut_region)
            
        if region_obj is None or not region_obj.valid_p():
            utils.echo_error_msg(f"Invalid cut region: {self.cut_region}")
            return self.src_dem, -1

        ## Apply Mask
        with gdalfun.gdal_datasource(self.dst_dem) as ds:
            if ds is not None:
                self.init_ds(ds)
                
                ## Get window of the cut region within the raster
                ## (x_off, y_off, x_size, y_size)
                cut_win = region_obj.srcwin(
                    self.gt, self.ds_config['nx'], self.ds_config['ny']
                )
                
                ## Create coordinate grids for the raster
                ## For large rasters, this might be memory intensive; 
                ## strictly speaking we only need to mask the array indices.
                
                band = ds.GetRasterBand(self.band)
                arr = band.ReadAsArray()
                
                ## Initialize a full "Keep" mask (True = Keep, False = Mask)
                ## If invert=False (Standard Cut): Keep INSIDE region (True inside, False outside)
                ## If invert=True (Hole Punch): Keep OUTSIDE region (False inside, True outside)
                
                ## Create boolean mask for the region window
                ## Initialize with opposite of desired default
                mask = np.zeros(arr.shape, dtype=bool) 
                
                ## Set the window area to True
                ## Clip window to array dims to be safe
                y_start = max(0, cut_win[1])
                y_end = min(arr.shape[0], cut_win[1] + cut_win[3])
                x_start = max(0, cut_win[0])
                x_end = min(arr.shape[1], cut_win[0] + cut_win[2])
                
                if x_start < x_end and y_start < y_end:
                    mask[y_start:y_end, x_start:x_end] = True
                
                ## Determine what to set to NoData
                ## We want to set data to NDV where we are NOT keeping it.
                
                if not self.invert:
                    ## Cut Mode: Keep Inside (Mask=True). Set NDV where Mask=False.
                    arr[~mask] = self.ds_config['ndv']
                else:
                    ## Invert Mode: Remove Inside (Mask=True). Set NDV where Mask=True.
                    arr[mask] = self.ds_config['ndv']

                band.WriteArray(arr)
                
        dst_ds = None
        return self.dst_dem, 0

### End
