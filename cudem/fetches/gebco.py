### gebco.py - GEBCO dataset
##
## Copyright (c) 2010 - 2022 Regents of the University of Colorado
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
## The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
##
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
## INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR 
## PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE 
## FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
## ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
##
### Commentary:
##
## GEBCO Fetch
##
### Code:

import os

from cudem import utils
from cudem import regions
from cudem import datasets

import cudem.fetches.utils as f_utils

class GEBCO(f_utils.FetchModule):
    '''Fetch raster data from GEBCO'''
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs) 

        self._gebco_url = 'https://www.bodc.ac.uk/data/open_download/gebco/gebco_2022/'
        self._gebco_subice_url = 'https://www.bodc.ac.uk/data/open_download/gebco/gebco_2022_sub_ice_topo/'
        self._nc = 'zip/'
        self._asc = 'esri_ascii/'
        self._tiff = 'geotiff/'
        self._outdir = os.path.join(os.getcwd(), 'gebco')
        self.name = 'gebco'
        
    def run(self):
        '''Run the GEBCO fetching module'''

        outf = 'gebco.zip'
        self.results.append([self._gebco_url + self._tiff, os.path.join(self._outdir, outf), 'gebco'])
                
        return(self)

    # def yield_xyz(self, entry):
    #     src_data = 'gebco_tmp.tif'
    #     if f_utils.Fetch(
    #             entry[0], callback=self.callback, verbose=self.verbose
    #     ).fetch_file(src_data) == 0:
    #         gebco_ds = datasets.RasterFile(
    #             fn=src_data,
    #             data_format=200,
    #             src_srs='epsg:4326',
    #             dst_srs=self.dst_srs,
    #             x_inc=self.x_inc,
    #             y_inc=self.y_inc,
    #             weight=self.weight,
    #             src_region=self.region,
    #             verbose=self.verbose
    #         )
    #         if self.bathy_only:
    #             for xyz in gmrt_ds.yield_xyz():
    #                 if xyz.z < 0:
    #                     yield(xyz)
    #         else:
    #             for xyz in gmrt_ds.yield_xyz():
    #                 yield(xyz)
                    
    #     else:
    #         utils.echo_error_msg('failed to fetch remote file, {}...'.format(src_data))
            
    #     utils.remove_glob('{}*'.format(src_data))

    # def yield_array(self, entry):
    #     src_data = 'gmrt_tmp.tif'
    #     if f_utils.Fetch(
    #             entry[0], callback=self.callback, verbose=self.verbose
    #     ).fetch_file(src_data) == 0:
    #         if self.bathy_only:
    #             ds = gdal.Open(src_data)
    #             ds_config = demfun.gather_infos(ds)
    #             band = ds.GetRasterBand(1)
    #             comp_geot = ds_config['geoT']
    #             outarray = ds.ReadAsArray()
    #             outarray[outarray > 0] = band.GetNoDataValue()
    #             band.WriteArray(outarray)
    #             ds = None

    #         gmrt_ds = datasets.RasterFile(
    #             fn=src_data,
    #             data_format=200,
    #             src_srs='epsg:4326',
    #             dst_srs=self.dst_srs,
    #             x_inc=self.x_inc,
    #             y_inc=self.y_inc,
    #             weight=self.weight,
    #             src_region=self.region,
    #             verbose=self.verbose
    #         )
                
    #         for arr in gmrt_ds.yield_array():
    #             yield(arr)
                    
    #     else:
    #         utils.echo_error_msg('failed to fetch remote file, {}...'.format(src_data))
            
    #     utils.remove_glob('{}*'.format(src_data))
        
### End
