### globatofile.py - DataLists IMproved
##
## Copyright (c) 2010 - 2026 Regents of the University of Colorado
##
## globatofile.py is part of CUDEM
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
## CUDEM Globato file parsing...
##
### Code:

import h5py
import numpy as np
from cudem import utils
from cudem import regions
from cudem.datalists.dlim import ElevationDataset

class GlobatoFile(ElevationDataset):
    """Parser for Globato HDF5 Stack files.
    Reads the 'stack' group (finalized data).
    """
    
    def __init__(self, layer='z', **kwargs):
        super().__init__(**kwargs)
        self.layer = layer # 'z', 'uncertainty', 'weight', etc.
        self.data_format = 320

    def generate_inf(self):
        """Globato files store bounds in attributes, fast parse."""
        try:
            with h5py.File(self.fn, 'r') as f:
                if 'crs' in f and 'GeoTransform' in f['crs'].attrs:
                    gt = [float(x) for x in f['crs'].attrs['GeoTransform'].split()]
                    ## GT: [minx, x_inc, 0, maxy, 0, y_inc]
                    ## Dimensions needed
                    if 'stack' in f and 'z' in f['stack']:
                        shape = f['stack']['z'].shape # (y, x)
                        
                        minx = gt[0]
                        maxy = gt[3]
                        maxx = minx + (gt[1] * shape[1])
                        miny = maxy + (gt[5] * shape[0])
                        
                        ## Z Bounds (scan or attribute? Globato doesn't seem to store zmin/max in attrs)
                        ## We'll stick to XY bounds for INF speed
                        ## TODO: Add zmin/zmax in globato attrs
                        self.infos.minmax = [minx, maxx, miny, maxy, 0, 0]
                        self.infos.numpts = shape[0] * shape[1]
                        
                        ## SRS
                        if 'crs_wkt' in f['crs'].attrs:
                            self.infos.src_srs = f['crs'].attrs['crs_wkt'].decode('utf-8')

        except Exception as e:
            utils.echo_warning_msg(f"Error reading Globato INF: {e}")
            
        return self.infos

    def yield_points(self):
        """Yields non-nan points from the Globato stack.
        """
        
        with h5py.File(self.fn, 'r') as f:
            try:
                ## Globato stores 1D lat/lon arrays
                lats = f['lat'][...]
                lons = f['lon'][...]
                data = f['stack'][self.layer][...] # 2D array
                
                ## Check for mask if available
                ## Or just filter NaNs
                valid_mask = ~np.isnan(data)
                
                ## Create coordinate meshgrids for valid pixels only? 
                ## Or iterate chunks to save RAM.
                
                ## Iterate scanlines (chunks) to yield points
                ## HDF5 is chunked, let's assume row-based iteration is safe
                
                height, width = data.shape
                
                ## Globato 'lat' array corresponds to rows
                ## Globato 'lon' array corresponds to cols                
                for i in range(height):
                    row_data = data[i, :]
                    row_valid = ~np.isnan(row_data)
                    
                    if not np.any(row_valid): continue
                    
                    y_val = lats[i]
                    x_vals = lons[row_valid]
                    z_vals = row_data[row_valid]
                    
                    ## Create records
                    ## Expand y_val to match size
                    y_vals = np.full(x_vals.shape, y_val)
                    
                    ds = np.rec.fromarrays(
                        [x_vals, y_vals, z_vals],
                        names=['x', 'y', 'z']
                    )
                    yield ds

            except Exception as e:
                utils.echo_error_msg(f"Globato read error: {e}")


### End
