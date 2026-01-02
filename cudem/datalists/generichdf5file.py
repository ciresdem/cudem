### generichdf5file.py - OSGEO functions
##
## Copyright (c) 2010 - 2026 Regents of the University of Colorado
##
## generichdf5file.py is part of CUDEM
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
## Generic hdf5 file parsing.
##
### Code:

import h5py
import numpy as np
from cudem import utils
from cudem.datalists.dlim import ElevationDataset

class GenericHDF5File(ElevationDataset):
    """
    Generic HDF5/NetCDF4 Parser.
    User must specify internal paths to variables.
    
    Parameters:
    -----------
    x_path : str
        Path to x/longitude variable (e.g. '/gt1l/geolocation/lon')
    y_path : str
        Path to y/latitude variable
    z_path : str
        Path to z/elevation variable
    """
    
    def __init__(self, x_path=None, y_path=None, z_path=None, **kwargs):
        super().__init__(**kwargs)
        self.x_path = x_path
        self.y_path = y_path
        self.z_path = z_path
        self.data_format = 310

    def generate_inf(self):
        """Scan header for simple stats (if variables are 1D and simple)."""
        if not all([self.x_path, self.y_path, self.z_path]):
            return super().generate_inf() # Fallback to scan

        try:
            with h5py.File(self.fn, 'r') as f:
                # Quick bounds check if datasets exist
                if self.x_path in f and self.y_path in f:
                    x = f[self.x_path]
                    y = f[self.y_path]
                    
                    # Heuristic: If 1D, easy to get min/max
                    # If 2D (grids), min/max of whole array
                    self.infos.minmax = [
                        np.min(x), np.max(x),
                        np.min(y), np.max(y),
                        0, 0 # Z is harder without reading
                    ]
                    self.infos.numpts = x.size
        except:
            pass
            
        return self.infos

    def yield_points(self):
        if not all([self.x_path, self.y_path, self.z_path]):
            utils.echo_error_msg("Missing variable paths for GenericHDF5File")
            return

        with h5py.File(self.fn, 'r') as f:
            try:
                # Read Data
                # Note: This loads fully into RAM. For massive files, 
                # we should implement slicing/chunking here.
                x = f[self.x_path][...]
                y = f[self.y_path][...]
                z = f[self.z_path][...]
                
                # Handle shapes (flatten if grid)
                if x.ndim > 1: x = x.flatten()
                if y.ndim > 1: y = y.flatten()
                if z.ndim > 1: z = z.flatten()
                
                # Handle Coordinates if 1D axes for a 2D grid
                # (e.g. NetCDF convention: lat(ny), lon(nx), z(ny, nx))
                if x.size != z.size and z.ndim == 2:
                    # Meshgrid
                    x, y = np.meshgrid(x, y)
                    x = x.flatten()
                    y = y.flatten()
                    z = z.flatten()

                # Yield
                ds = np.rec.fromarrays([x, y, z], names=['x', 'y', 'z'])
                yield ds

            except KeyError as e:
                utils.echo_error_msg(f"Variable not found in {self.fn}: {e}")

### End
