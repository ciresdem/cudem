### ncfile.py - DataLists IMproved
##
## Copyright (c) 2010 - 2026 Regents of the University of Colorado
##
## ncfile.py is part of CUDEM
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
## NetCDF Raster Functions
##
### Code:

import numpy as np
import os
try:
    import netCDF4 as nc
except ImportError:
    nc = None

def netcdf_write(src_arr, dst_fn, infos):
    """Write a numpy array to a CF-Compliant NetCDF file."""
    
    if nc is None:
        raise ImportError("netCDF4 module required for grid output.")

    ## Remove existing
    if os.path.exists(dst_fn):
        os.remove(dst_fn)

    ## Dimensions
    ny, nx = src_arr.shape
    
    ## Calculate Coordinate Arrays (Node or Pixel)
    ## infos['geoT'] = (x_min, x_inc, 0, y_max, 0, y_inc)
    gt = infos['geoT']
    x_range = np.linspace(gt[0], gt[0] + (gt[1] * nx), nx, endpoint=False)
    ## Adjust for half-pixel shift if pixel-registered, essentially:
    x_range += (gt[1] / 2.0) 
    
    ## Y is typically top-down in GT, NetCDF usually prefers bottom-up, 
    ## but we can match GT order (Top-Down) and set attributes accordingly.
    y_range = np.linspace(gt[3], gt[3] + (gt[5] * ny), ny, endpoint=False)
    y_range += (gt[5] / 2.0)

    with nc.Dataset(dst_fn, 'w', format='NETCDF4') as ds:
        ## Dimensions
        ds.createDimension('lat', ny)
        ds.createDimension('lon', nx)
        
        ## Variables
        lats = ds.createVariable('lat', 'f8', ('lat',))
        lons = ds.createVariable('lon', 'f8', ('lon',))
        z = ds.createVariable('z', 'f4', ('lat', 'lon'), zlib=True, complevel=4, fill_value=infos['ndv'])
        
        ## Attributes (CF Conventions)
        ds.Conventions = "CF-1.7"
        ds.title = "CUDEM Waffles Output"
        
        lats.units = "degrees_north"
        lats.standard_name = "latitude"
        lats[:] = y_range
        
        lons.units = "degrees_east"
        lons.standard_name = "longitude"
        lons[:] = x_range
        
        z.units = "meters"
        z.standard_name = "height"
        
        ## Write Data
        z[:] = src_arr
        
        ## Write CRS (WKT)
        crs = ds.createVariable('crs', 'i4')
        crs.spatial_ref = infos['wkt']
        z.grid_mapping = "crs"


### End
