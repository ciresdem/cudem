### uncertainties.py
##
## Copyright (c) 2010 - 2021 CIRES Coastal DEM Team
##
## Permission is hereby granted, free of charge, to any person obtaining a copy 
## of this software and associated documentation files (the "Software"), to deal 
## in the Software without restriction, including without limitation the rights 
## to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies 
## of the Software, and to permit persons to whom the Software is furnished to do so, 
## subject to the following conditions:
##
## The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
##
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
## INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR 
## PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE 
## FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
## ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
##
### Commentary:
### Code:

from cudem import utils
from cudem import regions
from cudem import dlim
from cudem import waffles

## ==============================================
## Waffles Interpolation Uncertainty module
## ==============================================
class InterpolationUncertainty: #(waffles.Waffle):

    def __init__(self, dem=None, percentile=95, sims=None, chnk_lvl=None):
        """calculate cell-level interpolation uncertainty

        Args:
          dem (WaffledRaster): a waffles generated DEM (or constructed WaffledRaster object)
          percentile (int): max percentile
          sims (int): number of split-sample simulations
          chnk_lvl (int): the 'chunk-level'
        """
        
        self.dem = dem
        self.percentile = percentile
        self.sims = sims
        self.chnk_lvl = chnk_lvl
        
        self._zones = ['low-dens','mid-dens','high-dens','low-slp','mid-slp','high-slp']

## ==============================================
## testing
## ==============================================
region = regions.Region().from_list([480000.01, 482999.99, 4395000.0, 4397999.99])
epsg = 26913

module = 'triangulate'
module_args = ()

inc = 1
name='CO_SoPlatteRiver_MergeTriangulate1'

datasets = [
    'USGS_LPC_CO_SoPlatteRiver_Lot5_2013_13SDD480395_LAS_2015.xyz',
    'USGS_LPC_CO_SoPlatteRiver_Lot5_2013_13SDD481395_LAS_2015.xyz',
    'USGS_LPC_CO_SoPlatteRiver_Lot5_2013_13SDD480396_LAS_2015.xyz',
    'USGS_LPC_CO_SoPlatteRiver_Lot5_2013_13SDD481396_LAS_2015.xyz',
]

waffles_dem = waffles.WaffleFactory(
    data = datasets,
    src_region = region,
    inc = inc,
    name = name,
    epsg = epsg,
).acquire_module_by_name(module, *module_args)

waffles_dem.fn = 'CO_SoPlatteRiver_MergeTriangulate1.tif'
waffles_dem.waffled = True

print(waffles_dem.valid_p())

i = InterpolationUncertainty(dem=waffles_dem)

print(i.dem)
