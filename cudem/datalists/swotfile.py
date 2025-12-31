### swotfile.py 
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
### Commentary:
##
## NASA SWOT Data Parser (HDF5 / NetCDF)
##
### Code:

import os
import numpy as np
import h5py as h5
from osgeo import gdal

from cudem import utils
from cudem import regions
from cudem import gdalfun
from cudem import srsfun
from cudem.datalists.dlim import ElevationDataset

class SWOTFile(ElevationDataset):
    """NASA SWOT Data super class (HDF5).
    Uses h5py to parse data. Subclass to process specific SWOT products.
    """
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        
    def _init_h5File(self, short_name='L2_HR_PIXC'):
        """Open HDF5 file and validate short_name."""
        
        src_h5 = None
        try:
            src_h5 = h5.File(self.fn, 'r')
            if 'short_name' in src_h5.attrs:
                ## Handle bytes vs string for attribute
                attr_name = src_h5.attrs['short_name']
                if isinstance(attr_name, bytes):
                    attr_name = attr_name.decode('utf-8')
                    
                if attr_name != short_name:
                    utils.echo_warning_msg(
                        f'{self.fn} short_name ({attr_name}) does not match expected ({short_name})'
                    )
            else:
                utils.echo_warning_msg(f'{self.fn} missing "short_name" attribute')
                
        except Exception as e:
            utils.echo_error_msg(f"Failed to open SWOT file {self.fn}: {e}")
            if src_h5: 
                src_h5.close()
                src_h5 = None

        return src_h5

    
    def _close_h5File(self, src_h5):
        if src_h5 is not None:
            try:
                src_h5.close()
            except Exception:
                pass

            
    def _get_var_arr(self, src_h5, var_path):
        """Safely retrieve variable array."""
        
        try:
            return src_h5[f'/{var_path}'][...]
        except KeyError:
            return None

        
class SWOT_PIXC(SWOTFile):
    """NASA SWOT Pixel Cloud (L2_HR_PIXC) data parser.
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
        self.remove_class_flags = remove_class_flags
        
        self.classes = self._parse_list(classes)
        self.classes_qual = self._parse_list(classes_qual)
        self.anc_classes = self._parse_list(anc_classes)

        
    def _parse_list(self, input_str):
        if input_str:
            return [int(x) for x in input_str.split('/')]
        return []

    
    def yield_points(self):
        """Yield filtered points from PIXC HDF5."""
        
        src_h5 = self._init_h5File(short_name='L2_HR_PIXC')
        if src_h5 is None: return

        try:
            ## Load Coordinates & Data
            lat = self._get_var_arr(src_h5, f'{self.group}/latitude')
            lon = self._get_var_arr(src_h5, f'{self.group}/longitude')
            val = self._get_var_arr(src_h5, f'{self.group}/{self.var}')
            
            if lat is None or lon is None or val is None:
                raise ValueError("Missing core variables (lat/lon/height)")

            ## Apply Geoid Correction
            if self.apply_geoid:
                geoid = self._get_var_arr(src_h5, f'{self.group}/geoid')
                if geoid is not None:
                    val = val - geoid

            ## Create Recarray
            ## Filter NaNs / Fill Values immediately
            valid_mask = val > -1e30 
            
            dataset = np.column_stack((lon[valid_mask], lat[valid_mask], val[valid_mask]))
            points = np.rec.fromrecords(dataset, names='x, y, z')
            
            ## --- Classification Filtering ---
            ## Indices must be aligned, so we apply the initial valid_mask to aux arrays too
            
            ## Classification
            if self.classes:
                cls_data = self._get_var_arr(src_h5, f'{self.group}/classification')
                if cls_data is not None:
                    cls_data = cls_data[valid_mask]
                    points = points[np.isin(cls_data, self.classes)]
                    
                    # Update masks for subsequent filters to stay aligned with 'points'
                    # Actually, better to accumulate a boolean mask on the subset?
                    # For simplicity in this stream logic, we filter 'points' directly
                    # but this assumes subsequent filter arrays are re-sliced. 
                    # Optimization: slice qual arrays *after* initial filtering if possible, 
                    # but here we might lose alignment.
                    # Strategy: Re-slice aux arrays based on the same mask logic.
                    
                    # NOTE: Re-slicing aux arrays is expensive. 
                    # Let's build a master boolean mask first.
            
            ## Optimization: Build composite mask first
            final_mask = np.ones(len(points), dtype=bool) # All valid from valid_mask
            
            ## Load Aux Data if needed (already masked by valid_mask)
            if self.classes or self.remove_class_flags or self.classes_qual:
                cls_data = self._get_var_arr(src_h5, f'{self.group}/classification')
                if cls_data is not None: cls_data = cls_data[valid_mask]
                
                qual_data = self._get_var_arr(src_h5, f'{self.group}/classification_qual')
                if qual_data is not None: qual_data = qual_data[valid_mask]

                if self.classes and cls_data is not None:
                    final_mask &= np.isin(cls_data, self.classes)
                
                if self.remove_class_flags and qual_data is not None:
                    final_mask &= (qual_data == 0)
                elif self.classes_qual and qual_data is not None:
                    final_mask &= (~np.isin(qual_data, self.classes_qual))

            if self.anc_classes:
                anc_data = self._get_var_arr(src_h5, f'{self.group}/ancillary_surface_classification_flag')
                if anc_data is not None:
                    anc_data = anc_data[valid_mask]
                    final_mask &= np.isin(anc_data, self.anc_classes)

            ## Apply Final Mask
            yield points[final_mask]

        except Exception as e:
            utils.echo_error_msg(f"Error processing SWOT PIXC: {e}")
        finally:
            self._close_h5File(src_h5)

            
class SWOT_HR_Raster(ElevationDataset):
    """NASA SWOT HR_Raster data parser (NetCDF/GeoTIFF via GDAL).
    """
        
    def __init__(self, data_set='wse', **kwargs):
        super().__init__(**kwargs)
        self.data_set = data_set

        
    def parse(self):
        """Parse subdatasets from the SWOT Raster."""

        from cudem.datasets import DatasetFactory
        
        src_ds = gdal.Open(self.fn)
        if src_ds is None: return

        sub_datasets = src_ds.GetSubDatasets()
        target_sub = None
        
        ## Find requested subdataset
        if utils.int_or(self.data_set) is not None:
            idx = int(self.data_set)
            if 0 <= idx < len(sub_datasets):
                target_sub = sub_datasets[idx]
        else:
            for sd in sub_datasets:
                ## Name format often "NETCDF:filename:varname"
                if sd[0].endswith(f":{self.data_set}") or self.data_set in sd[0]:
                    target_sub = sd
                    break
        
        ## Default fallback
        if target_sub is None and len(sub_datasets) > 2:
            target_sub = sub_datasets[2] # Often WSE

        if target_sub:
            ## Determine SRS override for WSE (Water Surface Elevation)
            ## SWOT WSE is usually relative to ellipsoid, often we want it processed
            src_srs = gdalfun.gdal_get_srs(target_sub[0])
            
            if self.data_set == 'wse':
                src_srs = srsfun.combine_epsgs(
                    src_srs, '3855', name='SWOT Combined'
                )
                
            sub_ds = DatasetFactory(
                **self._set_params(
                    mod=target_sub[0],
                    data_format=200,
                    node='grid',
                    check_path=False
                )
            )._acquire_module()
            
            # sub_ds = DatasetFactory(
            #     mod=target_sub[0],
            #     data_format=200, # Treat as GDAL Raster
            #     node='grid',
            #     check_path=False,
            #     src_srs=src_srs
            # )._acquire_module()
            
            self.data_entries.append(sub_ds)
            sub_ds.initialize()
            
            ## Delegate to GDALFile parser
            for gdal_ds in sub_ds.parse():
                yield gdal_ds

### End
