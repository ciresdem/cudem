### gebco.py
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
##
## gebco.py is part of CUDEM
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
## Fetch global bathymetry from GEBCO.
##
### Code:

from typing import Optional, Union, List
from cudem import utils
from cudem.fetches import fetches

## ==============================================
## Constants
## ==============================================
GEBCO_URLS = {
    'gebco_ice': {
        'netcdf': 'https://www.bodc.ac.uk/data/open_download/gebco/gebco_2022/zip/',
        'geotiff': 'https://www.bodc.ac.uk/data/open_download/gebco/gebco_2022/geotiff/',
        'ascii': 'https://www.bodc.ac.uk/data/open_download/gebco/gebco_2022/esri_ascii/'
    },
    'gebco_sub_ice': {
        'netcdf': 'https://www.bodc.ac.uk/data/open_download/gebco/gebco_2022_sub_ice_topo/zip/',
        'geotiff': 'https://www.bodc.ac.uk/data/open_download/gebco/gebco_2022_sub_ice_topo/geotiff/',
        'ascii': 'https://www.bodc.ac.uk/data/open_download/gebco/gebco_2022_sub_ice_topo/esri_ascii/'
    },
    'gebco_tid': {
        'netcdf': 'https://www.bodc.ac.uk/data/open_download/gebco/gebco_2022_tid/zip/',
        'geotiff': 'https://www.bodc.ac.uk/data/open_download/gebco/gebco_2022_tid/geotiff/',
        'ascii': 'https://www.bodc.ac.uk/data/open_download/gebco/gebco_2022_tid/esri_ascii/'
    }
}

## Type Identifier (TID) Dictionary: {TID: [Description, Weight]}
TID_INFO = {
    0: ['Land', 0.21],
    10: ['Singlebeam - depth value collected by a single beam echo-sounder', 0.65],
    11: ['Multibeam - depth value collected by a multibeam echo-sounder', 0.95],
    12: ['Seismic - depth value collected by seismic methods', 0.3],
    13: ['Isolated sounding - depth value that is not part of a regular survey or trackline', 0.4],
    14: ['ENC sounding - depth value extracted from an Electronic Navigation Chart (ENC)', 0.5],
    15: ['Lidar - depth derived from a bathymetric lidar sensor', 1.0],
    16: ['Depth measured by optical light sensor', 0.3],
    17: ['Combination of direct measurement methods', 0.2],
    40: ['Predicted based on satellite-derived gravity data', 0.19],
    41: ['Interpolated based on a computer algorithm', 0.18],
    42: ['Digital bathymetric contours from charts', 0.17],
    43: ['Digital bathymetric contours from ENCs', 0.16],
    44: ['Bathymetric sounding - constrained by bathymetric sounding(s)', 0.15],
    45: ['Predicted based on helicopter/flight-derived gravity data', 0.14],
    46: ['Depth estimated by calculating the draft of a grounded iceberg', 0.13],
    70: ['Pre-generated grid - mixed source data types', 0.12],
    71: ['Unknown source', 0.11],
    72: ['Steering points', 0.1]
}

## ==============================================
## GEBCO Module
## ==============================================
class GEBCO(fetches.FetchModule):
    """GEneral Bathymetric Chart of the Oceans (GEBCO)
    
    GEBCO’s current gridded bathymetric data set, the GEBCO_2022 Grid, 
    is a global terrain model for ocean and land, providing elevation data 
    on a 15 arc-second interval grid.

    Currently only fetches entire grid (zip files).

    Configuration Example:
    < gebco:want_ice=geotiff:want_sub_ice=False:want_tid=False:exclude_tid=None:upper_limit=None:lower_limit=None >
    """
    
    def __init__(
            self,
            want_ice: Union[str, bool] = 'geotiff',
            want_sub_ice: Union[str, bool] = False,
            want_tid: Union[str, bool] = False,
            exclude_tid: Optional[str] = None,
            upper_limit: Optional[float] = None,
            lower_limit: Optional[float] = None,
            **kwargs
    ):
        super().__init__(name='gebco', **kwargs)
        
        ## Configuration for data types (defaulting to 'geotiff' if True is passed)
        self.want_ice = utils.str_or(want_ice, 'geotiff') if want_ice else False
        self.want_sub_ice = utils.str_or(want_sub_ice, 'geotiff') if want_sub_ice else False
        self.want_tid = utils.str_or(want_tid, 'geotiff') if want_tid else False 

        ## Process excluded TIDs
        self.exclude_tid = []
        if exclude_tid:
            self.exclude_tid = [utils.int_or(tid) for tid in str(exclude_tid).split('/')]
            self.exclude_tid = [tid for tid in self.exclude_tid if tid is not None]

        ## Region configuration
        if self.region:
            self.gebco_region = self.region.copy()
            self.gebco_region.zmax = utils.float_or(upper_limit)
            self.gebco_region.zmin = utils.float_or(lower_limit)
        else:
            self.gebco_region = None

        ## Data format (-2 is zipfile) and SRS
        self.data_format = -2
        self.src_srs = 'epsg:4326+3855'
        
        ## Reference constant for use in processing
        self.tid_dic = TID_INFO

        
    def run(self):
        """Run the GEBCO fetching module."""
        
        if self.want_ice:
            self.add_entry_to_results(
                GEBCO_URLS['gebco_ice'][self.want_ice],
                'gebco_ice.zip',
                'gebco'
            )
            
        if self.want_sub_ice:
            self.add_entry_to_results(
                GEBCO_URLS['gebco_sub_ice'][self.want_sub_ice],
                'gebco_sub_ice.zip',
                'gebco'
            )
            
        if self.want_tid:
            self.add_entry_to_results(
                GEBCO_URLS['gebco_tid'][self.want_tid],
                'gebco_tid.zip',
                'gebco'
            )
                
        return self

### End
