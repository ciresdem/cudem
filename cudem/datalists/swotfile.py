### swotfile.py - DataLists IMproved
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
##
## swotfile.py is part of CUDEM
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
### Examples:
##
### TODO:
##
### Code:

import os
import numpy as np
import h5py as h5

from osgeo import gdal

from cudem import utils
from cudem import regions
from cudem import gdalfun
from cudem.datalists.dlim import ElevationDataset

## NASA SWOT Data class (hdf5)
## uses h5py
class SWOTFile(ElevationDataset):
    """NASA SWOT Data super class

    Uses h5py to parse data. Make a subclass of this to 
    process various types of SWOT data
    """
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        
    def _init_h5File(self, short_name='L2_HR_PIXC'):
        src_h5 = None
        try:
            src_h5 = h5.File(self.fn, 'r')
            if src_h5 is not None:
                if 'short_name' in src_h5.attrs.keys():
                    if src_h5.attrs['short_name'] != short_name.encode('utf-8'):
                        utils.echo_error_msg(
                            (f'{self.fn} does not appear to be a '
                             f'SWOT {short_name} file')
                        )
                        self._close_h5File(src_h5)
                else:
                    utils.echo_error_msg(
                        f'{self.fn} does not appear to be a SWOT file'
                    )
                    self._close_h5File(src_h5)
                    
        except Exception as e:
            utils.echo_error_msg(e)

        return(src_h5)

    
    def _close_h5File(self, src_h5):
        if src_h5 is not None:
            src_h5.close()

            
    def _get_var_arr(self, src_h5, var_path):
        return(src_h5['/{}'.format(var_path)][...,])

    
class SWOT_PIXC(SWOTFile):
    """NASA SWOT PIXC data file.

    Extract data from a SWOT PIXC file.

    classes: 1UB, 2UB, 3UB, 4UB, 5UB, 6UB, 7UB
    "land, land_near_water, water_near_land, open_water, 
    dark_water, low_coh_water_near_land, open_low_coh_water"

    classes_qual: 1U, 2U, 4U, 8U, 16U, 2048U, 8192U, 16384U, 
    32768U, 262144U, 524288U, 134217728U, 536870912U, 
    1073741824U, 2147483648U

    "no_coherent_gain power_close_to_noise_floor 
    detected_water_but_no_prior_water detected_water_but_bright_land 
    water_false_detection_rate_suspect coherent_power_suspect 
    tvp_suspect sc_event_suspect small_karin_gap in_air_pixel_degraded 
    specular_ringing_degraded coherent_power_bad tvp_bad sc_event_bad 
    large_karin_gap"

    anc_classes: 0UB, 1UB, 2UB, 3UB, 4UB, 5UB, 6UB 
    "open_ocean land continental_water aquatic_vegetation 
    continental_ice_snow floating_ice salted_basin"  		
    """
    
    def __init__(self,
                 group='pixel_cloud',
                 var='height',
                 apply_geoid=True,
                 classes=None,
                 classes_qual=None,
                 anc_classes=None,
                 remove_class_flags=False,
                 **kwargs):
        super().__init__(**kwargs)
        self.group = group
        self.var = var
        self.apply_geoid = apply_geoid
        self.classes = [int(x) for x in classes.split('/')] \
            if classes is not None \
               else []
        self.classes_qual = [int(x) for x in classes_qual.split('/')] \
            if classes_qual is not None \
               else []
        self.anc_classes = [int(x) for x in anc_classes.split('/')] \
            if anc_classes is not None \
               else []
        self.remove_class_flags = remove_class_flags
        # if self.remove_class_flags:
        #     self.classes_qual = [1, 2, 4, 8, 16, 2048, 8192, 16384, 32768,
        #                          262144, 524288, 134217728, 536870912,
        #                          1073741824, 2147483648]

        
    def yield_points(self):
        src_h5 = self._init_h5File(short_name='L2_HR_PIXC')
        src_h5_vec = None
        
        #if self.pixc_vec is not None:
        #    src_h5_vec = self._init_h5File(short_name='L2_HR_PIXCVec')
        
        if src_h5 is not None:
            latitude = self._get_var_arr(
                src_h5, '{}/latitude'.format(self.group)
            )
            longitude = self._get_var_arr(
                src_h5, '{}/longitude'.format(self.group)
            )
            var_data = self._get_var_arr(
                src_h5, '{}/{}'.format(self.group, self.var)
            )
            if self.apply_geoid:
                geoid_data = self._get_var_arr(
                    src_h5, '{}/geoid'.format(self.group)
                )
                out_data = var_data - geoid_data
            else:
                out_data = var_data
                
            dataset = np.column_stack(
                (longitude, latitude, out_data)
            )
            points = np.rec.fromrecords(
                dataset, names='x, y, z'
            )
            #points = points[points['z'] != 9.96921e+36]

            ## Classification Filter
            if len(self.classes) > 0:
                class_data = self._get_var_arr(
                    src_h5, f'{self.group}/classification'
                )
                points = points[(np.isin(class_data, self.classes))]

                ## Classification Quality Filter
                if self.remove_class_flags:
                    class_qual_data = self._get_var_arr(
                        src_h5, f'{self.group}/classification_qual'
                    )
                    class_qual_data = class_qual_data[
                        (np.isin(class_data, self.classes))
                    ]
                    points = points[class_qual_data == 0]
                                   
                elif len(self.classes_qual) > 0:
                    class_qual_data = self._get_var_arr(
                        src_h5, f'{self.group}/classification_qual'
                    )
                    class_qual_data = class_qual_data[
                        (np.isin(class_data, self.classes))
                    ]
                    points = points[
                        (~np.isin(class_qual_data, self.classes_qual))
                    ]
                
            ## Ancilliary Classification Filter
            if len(self.anc_classes) > 0:
                anc_class_data = self._get_var_arr(
                    src_h5, f'{self.group}/ancillary_surface_classification_flag'
                )
                points = points[(np.isin(anc_class_data, self.anc_classes))]
                
            points = points[points['z'] != -9.969209968386869e+36]
            self._close_h5File(src_h5)
            self._close_h5File(src_h5_vec)

            yield(points)

            
## todo: update to h5
class SWOT_HR_Raster(ElevationDataset):
    """NASA SWOT HR_Raster data file.

    Extract data from a SWOT HR_Raster file.
    """
        
    def __init__(self, data_set='wse', **kwargs):
        super().__init__(**kwargs)
        self.data_set = data_set

        
    def parse(self):
        from .dlim import DatasetFactory
        
        src_ds = gdal.Open(self.fn)
        if src_ds is not None:
            sub_datasets = src_ds.GetSubDatasets()
            idx = 2
            if utils.int_or(self.data_set) is not None:
                idx = utils.int_or(self.data_set)
            else:
                for j, sd in enumerate(sub_datasets):
                    _name = sd[0].split(':')[-1]
                    if self.data_set == _name:
                        idx = j
                        break

            src_ds = None
            src_srs = gdalfun.gdal_get_srs(sub_datasets[idx][0])
            if self.data_set == 'wse':
                src_srs = gdalfun.combine_epsgs(
                    src_srs, '3855', name='SWOT Combined'
                )
            sub_ds = DatasetFactory(
                **self._set_params(
                    mod=sub_datasets[idx][0],
                    data_format=200,
                    node='grid',
                    check_path=False
                )
            )._acquire_module()
            self.data_entries.append(sub_ds)
            sub_ds.initialize()
            for gdal_ds in sub_ds.parse():
                yield(gdal_ds)                                  


### End
