### gebco.py
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
##
## fetches.py is part of CUDEM
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
###############################################################################
### Commentary:
##
##
### Code:

from cudem import utils
from cudem.fetches import fetches

## GEBCO
class GEBCO(fetches.FetchModule):
    """GEneral Bathymetric Chart of the Oceans (GEBCO)
    
    GEBCO’s current gridded bathymetric data set, the GEBCO_2022 Grid, 
    is a global terrain model for ocean and land, providing elevation data, 
    in meters, on a 15 arc-second interval grid. It is accompanied by a 
    Type Identifier (TID) Grid that gives information on the types of 
    source data that the GEBCO_2022 Grid is based. 

    https://www.gebco.net

    Currently only fetches entire grid. Subset in dlim, or elsewhere.

    < gebco:want_ice=geotiff:want_sub_ice=False:want_tid=False:exclude_tid=None:upper_limit=None:lower_limit=None >
    """
    
    def __init__(
            self,
            want_ice='geotiff',
            want_sub_ice=False,
            want_tid=False,
            exclude_tid=None,
            upper_limit=None,
            lower_limit=None,
            **kwargs
    ):
        super().__init__(name='gebco', **kwargs)
        
        ## various gebco URLs
        self._gebco_urls = {
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

        # ice surface
        self.want_ice = utils.str_or(want_ice, 'geotiff') if want_ice else False
        # sub-ice surface
        self.want_sub_ice = utils.str_or(want_sub_ice, 'geotiff') if want_sub_ice else False
        # source id grid
        self.want_tid = utils.str_or(want_tid, 'geotiff') if want_tid else False 

        ## see tid_dic for a list of the tid values/descriptions.
        exclude_tid = utils.str_or(exclude_tid)
        self.exclude_tid = []
        if exclude_tid is not None:
            for tid_key in exclude_tid.split('/'):
                self.exclude_tid.append(utils.int_or(tid_key))
                
        try:
            self.exclude_tid.remove(None)
        except: pass

        ## this dictionary holds the TID values as the keys, and the values are a list of [description, weight]
        ## the weights are used in dlim/waffles for processing different TID values appropriately.
        self.tid_dic = {
            0: ['Land', .21],
            10: ['Singlebeam - depth value collected by a single beam echo-sounder', .65],
            11:	['Multibeam - depth value collected by a multibeam echo-sounder', .95],
            12:	['Seismic - depth value collected by seismic methods', .3],
            13:	[('Isolated sounding - depth value that is not part of a regular '
                  'survey or trackline'), .4],
            14:	[('ENC sounding - depth value extracted from an '
                  'Electronic Navigation Chart (ENC)'), .5],
            15:	['Lidar - depth derived from a bathymetric lidar sensor', 1],
            16:	['Depth measured by optical light sensor', .3],
            17:	['Combination of direct measurement methods', .2],
            40:	[('Predicted based on satellite-derived gravity data - '
                  'depth value is an interpolated value guided by '
                  'satellite-derived gravity data'), .19],
            41:	[('Interpolated based on a computer algorithm - '
                  'depth value is an interpolated value based on a computer '
                  'algorithm (e.g. Generic Mapping Tools)'), .18],
            42:	[('Digital bathymetric contours from charts - '
                  'depth value taken from a bathymetric contour data set'), .17],
            43:	[('Digital bathymetric contours from ENCs - '
                  'depth value taken from bathymetric contours from an '
                  'Electronic Navigation Chart (ENC)'), .16],
            44:	[('Bathymetric sounding - depth value at this location is '
                  'constrained by bathymetric sounding(s) within a gridded '
                  'data set where interpolation between sounding points is '
                  'guided by satellite-derived gravity data'), .15],
            45:	['Predicted based on helicopter/flight-derived gravity data', .14],
            46:	[('Depth estimated by calculating the draft of a grounded '
                  'iceberg using satellite-derived freeboard measurement.'), .13],
            70:	[('Pre-generated grid - depth value is taken from a pre-generated '
                  'grid that is based on mixed source data types, e.g. single beam, '
                  'multibeam, interpolation etc.'), .12],
            71:	['Unknown source - depth value from an unknown source', .11],
            72:	[('Steering points - depth value used to constrain the grid '
                  'in areas of poor data coverage'), .1]
        }

        ## set the fetching region, restrict by z-region if desired.
        self.gebco_region = self.region.copy()
        self.gebco_region.zmax = utils.float_or(upper_limit)
        self.gebco_region.zmin = utils.float_or(lower_limit)

        ## for dlim, data format is -2 for a zip file.
        self.data_format = -2
        self.src_srs = 'epsg:4326+3855'
        
    def run(self):
        """Run the GEBCO fetching module"""

        outf = 'gebco.zip'
        if self.want_ice:
            self.add_entry_to_results(
                self._gebco_urls['gebco_ice'][self.want_ice],
                'gebco_ice.zip',
                'gebco'
            )
            
        if self.want_sub_ice:
            self.add_entry_to_results(
                self._gebco_urls['gebco_sub_ice'][self.want_sub_ice],
                'gebco_sub_ice.zip',
                'gebco'
            )
            
        if self.want_tid:
            self.add_entry_to_results(
                self._gebco_urls['gebco_tid'][self.want_tid],
                'gebco_tid.zip',
                'gebco'
            )
                
        return(self)

### End
