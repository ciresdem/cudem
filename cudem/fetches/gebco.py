### gebco.py - GEBCO dataset
##
## Copyright (c) 2010 - 2023 Regents of the University of Colorado
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

from osgeo import gdal

from cudem import utils
from cudem import regions
from cudem import datasets
from cudem import dlim
from cudem import demfun

import cudem.fetches.utils as f_utils

class GEBCO(f_utils.FetchModule):
    '''Fetch raster data from GEBCO'''
    
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

        self.want_ice = utils.str_or(want_ice, 'geotiff') if want_ice else False
        self.want_sub_ice = utils.str_or(want_sub_ice, 'geotiff') if want_sub_ice else False
        self.want_tid = utils.str_or(want_tid, 'geotiff') if want_tid else False

        exclude_tid = utils.str_or(exclude_tid)
        self.exclude_tid = []
        if exclude_tid is not None:
            for tid_key in exclude_tid.split('/'):
                self.exclude_tid.append(utils.int_or(tid_key))
                
        try:
            self.exclude_tid.remove(None)
        except:
            pass

        self.tid_dic = {
            0: ['Land', .21],
            10: ['Singlebeam - depth value collected by a single beam echo-sounder', .65],
            11:	['Multibeam - depth value collected by a multibeam echo-sounder', .95],
            12:	['Seismic - depth value collected by seismic methods', .3],
            13:	['Isolated sounding - depth value that is not part of a regular survey or trackline', .4],
            14:	['ENC sounding - depth value extracted from an Electronic Navigation Chart (ENC)', .5],
            15:	['Lidar - depth derived from a bathymetric lidar sensor', 1],
            16:	['Depth measured by optical light sensor', .3],
            17:	['Combination of direct measurement methods', .2],
            40:	['Predicted based on satellite-derived gravity data - depth value is an interpolated value guided by satellite-derived gravity data', .19],
            41:	['Interpolated based on a computer algorithm - depth value is an interpolated value based on a computer algorithm (e.g. Generic Mapping Tools)', .18],
            42:	['Digital bathymetric contours from charts - depth value taken from a bathymetric contour data set', .17],
            43:	['Digital bathymetric contours from ENCs - depth value taken from bathymetric contours from an Electronic Navigation Chart (ENC)', .16],
            44:	['Bathymetric sounding - depth value at this location is constrained by bathymetric sounding(s) within a gridded data set where interpolation between sounding points is guided by satellite-derived gravity data', .15],
            45:	['Predicted based on helicopter/flight-derived gravity data', .14],
            46:	['Depth estimated by calculating the draft of a grounded iceberg using satellite-derived freeboard measurement.', .13],
            70:	['Pre-generated grid - depth value is taken from a pre-generated grid that is based on mixed source data types, e.g. single beam, multibeam, interpolation etc.', .12],
            71:	['Unknown source - depth value from an unknown source', .11],
            72:	['Steering points - depth value used to constrain the grid in areas of poor data coverage', .1]
        }

        self.gebco_region = self.region.copy()
        self.gebco_region.zmax = utils.float_or(upper_limit)
        self.gebco_region.zmin = utils.float_or(lower_limit)
        
    def run(self):
        '''Run the GEBCO fetching module'''

        outf = 'gebco.zip'
        if self.want_ice:
            self.results.append([self._gebco_urls['gebco_ice'][self.want_ice], os.path.join(self._outdir, 'gebco_ice.zip'), 'gebco'])
        if self.want_sub_ice:
            self.results.append([self._gebco_urls['gebco_sub_ice'][self.want_sub_ice], os.path.join(self._outdir, 'gebco_sub_ice.zip'), 'gebco'])
        if self.want_tid:
            self.results.append([self._gebco_urls['gebco_tid'][self.want_tid], os.path.join(self._outdir, 'gebco_tid.zip'), 'gebco'])
                
        return(self)

    def yield_ds(self, entry):
        src_data = 'tmp_gebco.tif'
        if f_utils.Fetch(
                entry[0], callback=self.callback, verbose=self.verbose
        ).fetch_file(entry[1]) == 0:

            gebco_fns = utils.p_unzip(entry[1], ['tif'], self._outdir)
            ## fetch the TID zip if needed
            if self.exclude_tid:
                if f_utils.Fetch(
                        self._gebco_urls['gebco_tid']['geotiff'], callback=self.callback, verbose=self.verbose
                ).fetch_file(os.path.join(self._outdir, 'gebco_tid.zip')) == 0:

                    ## only extract the file(s) needed for the region...
                    tid_fns = utils.p_unzip(os.path.join(self._outdir, 'gebco_tid.zip'), ['tif'], self._outdir)
                    for tid_fn in tid_fns:
                        tmp_tid = os.path.join(self._outdir, 'tmp_tid.tif')
                        tid_ds = gdal.Open(tid_fn)
                        tid_config = demfun.gather_infos(tid_ds)
                        tid_band = tid_ds.GetRasterBand(1)
                        tid_array = tid_band.ReadAsArray().astype(float)
                        tid_ds = None
                        tid_config['ndv'] = -9999
                        tid_config['dt'] = gdal.GDT_Float32                        
                        for tid_key in self.exclude_tid:
                            tid_array[tid_array == tid_key] = tid_config['ndv']
                        
                        for tid_key in self.tid_dic.keys():
                            tid_array[tid_array == tid_key] = self.tid_dic[tid_key][1]
                            
                        utils.gdal_write(tid_array, tmp_tid, tid_config)
                        gebco_ds = datasets.RasterFile(
                            fn=tid_fn.replace('tid_', ''),
                            data_format=200,
                            src_srs='epsg:4326+3855',
                            dst_srs=self.dst_srs,
                            x_inc=self.x_inc,
                            y_inc=self.y_inc,
                            weight=self.weight,
                            src_region=self.gebco_region,
                            verbose=self.verbose,
                            mask=tmp_tid,
                            weight_mask=tmp_tid
                        )
                        yield(gebco_ds)
            else:
                for gebco_fn in gebco_fns:                    
                    gebco_ds = datasets.RasterFile(
                        fn=gebco_fn,
                        data_format=200,
                        src_srs='epsg:4326+3855',
                        dst_srs=self.dst_srs,
                        x_inc=self.x_inc,
                        y_inc=self.y_inc,
                        weight=self.weight,
                        src_region=self.gebco_region,
                        verbose=self.verbose,
                        #weight_mask=tmp_tid
                    )
                    yield(gebco_ds)
        
    def yield_array(self, entry):
        for ds in self.yield_ds(entry):
            for arr in ds.yield_array():
                yield(arr)

    def yield_xyz(self, entry):
        for ds in self.yield_ds(entry):
            for xyz in ds.yield_xyz():
                yield(xyz)
            
    def yield_zip_xyz(self, entry):
        if f_utils.Fetch(
                entry[0], callback=self.callback, verbose=self.verbose
        ).fetch_file(entry[1]) == 0:
            gebco_ds = dlim.ZIPlist(
                fn=entry[1],
                data_format=-2,
                src_srs='epsg:4326',
                dst_srs=self.dst_srs,
                x_inc=self.x_inc,
                y_inc=self.y_inc,
                weight=self.weight,
                src_region=self.region,
                verbose=self.verbose
            )
            for xyz in gebco_ds.yield_xyz():
                yield(xyz)
                    
        else:
            utils.echo_error_msg('failed to fetch remote file, {}...'.format(entry[0]))

    def yield_zip_array(self, entry):
        if f_utils.Fetch(
                entry[0], callback=self.callback, verbose=self.verbose
        ).fetch_file(entry[1]) == 0:
            gebco_ds = dlim.ZIPlist(
                fn=entry[1],
                data_format=-2,
                src_srs='epsg:4326',
                dst_srs=self.dst_srs,
                x_inc=self.x_inc,
                y_inc=self.y_inc,
                weight=self.weight,
                src_region=self.region,
                verbose=self.verbose
            )
            for arr in gebco_ds.yield_array():
                yield(arr)
                    
        else:
            utils.echo_error_msg('failed to fetch remote file, {}...'.format(entry[0]))
                              
### End
