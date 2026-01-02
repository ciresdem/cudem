### grits.py - GRId filTerS
##
## Copyright (c) 2024 - 2025 Regents of the University of Colorado
##
## grits.py is part of cudem
##
## Permission is hereby granted, free of charge, to any person obtaining a copy 
## of this software and associated documentation files (the "Software"), to deal 
## in the Software without restriction, including without limitation the rights 
## to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies 
## of the Software, and to permit persons to whom the Software is furnished to do so, 
## subject to the following conditions:
##
## The above copyright notice and this permission notice shall be included in all copies or
## substantial portions of the Software.
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
## Grits: The CUDEM Grid Filtering Engine.
##
## Grits (GRId filTerS) provides a standardized framework for manipulating and
## filtering raster DEMs. It is designed to clean artifacts, smooth noise, and
## perform morphological operations on massive datasets via chunked processing.
##
##   * Smoothing & Denoising:
##      - Gaussian Blur, Median/Bilateral denoising, and GMT filter wrappers
##        to remove high-frequency noise.
##
##   * Cleaning & Restoration:
##      - Statistical outlier removal (IQR-based).
##      - Identification and removal of "flats" (artifacts).
##      - Gap filling (inpainting) using interpolation or spline methods.
##
##   * Advanced Morphology & Hydrology:
##      - Erosion, Dilation, Opening, and Closing for structural manipulation.
##      - Hydro-enforcement (breaching/filling) and river discovery.
##
##   * Integration:
##      - Seamless blending of overlapping datasets.
##      - Quality-based masking and cutting/clipping to vector polygons.
##
## Usage:
##   CLI: grits input.tif -M outliers:k=2.5 -M fill:method=spline -O clean.tif
##   API: g = GritsFactory(mod='blur:sigma=2', src_dem='dem.tif').acquire()
##        g.run()
##
### Code:

import os
import sys
import traceback
import argparse
import numpy as np
from osgeo import gdal

from cudem import utils
from cudem import regions
from cudem import gdalfun
from cudem import factory
from cudem.datalists.dlim import DatasetFactory
from cudem import __version__ as __cudem_version__
from . import __version__

class Grits:
    """DEM Filtering Base Class.
    Define a sub-class to create a new filter.
    """
    
    def __init__(
            self,
            src_dem: str = None,
            dst_dem: str = None,
            src_region: any = None,
            band: int = 1,
            min_z: float = None,
            max_z: float = None,
            min_weight: float = None,
            max_weight: float = None,
            min_uncertainty: float = None,
            max_uncertainty: float = None,
            count_mask: any = None,
            weight_mask: any = None,
            uncertainty_mask: any = None,
            aux_data: list = None,
            cache_dir: str = './',
            verbose: bool = True,
            params: dict = None,
            **kwargs: any
    ):
        self.src_dem = src_dem
        self.dst_dem = dst_dem
        self.src_region = src_region
        self.band = utils.int_or(band, 1)
        self.ds_config = None
        
        self.min_z = utils.float_or(min_z)
        self.max_z = utils.float_or(max_z)
        self.min_weight = utils.float_or(min_weight)
        self.max_weight = utils.float_or(max_weight)
        self.min_uncertainty = utils.float_or(min_uncertainty)
        self.max_uncertainty = utils.float_or(max_uncertainty)
        
        self.count_mask = count_mask
        self.weight_mask = weight_mask
        self.uncertainty_mask = uncertainty_mask
        self.aux_data = aux_data if aux_data else []
        
        self.cache_dir = cache_dir
        self.verbose = verbose
        self.params = params if params else {}
        self.kwargs = kwargs

        ## Output Filename
        if self.dst_dem is None:
            if self.src_dem:
                base = utils.fn_basename2(self.src_dem)
                ext = utils.fn_ext(self.src_dem)
                self.dst_dem = utils.make_temp_fn(f'{base}_filtered.{ext}', temp_dir=self.cache_dir)
            else:
                self.dst_dem = 'grits_filtered.tif'

        ## Auto-detect Globato HDF5 features
        if self.src_dem and self.src_dem.endswith('.h5'):
            self._init_globato()
            
        ## Initialize flags
        self.weight_is_fn = False
        self.weight_is_band = False
        self.uncertainty_is_fn = False
        self.uncertainty_is_band = False
        self.count_is_fn = False
        self.count_is_band = False
        
        ## Analyze masks
        self._analyze_mask('weight_mask', 'weight_is_fn', 'weight_is_band')
        self._analyze_mask('uncertainty_mask', 'uncertainty_is_fn', 'uncertainty_is_band')
        self._analyze_mask('count_mask', 'count_is_fn', 'count_is_band')

        ## initialize the souce dataset info
        #self.init_ds(self.src_dem)

        
    def _analyze_mask(self, attr_name, fn_flag, band_flag):
        """Helper to determine if a mask is a file or a band index."""
        
        val = getattr(self, attr_name)
        if val is not None:
            if utils.int_or(val) is not None:
                setattr(self, band_flag, True)
                setattr(self, attr_name, int(val))
            elif os.path.exists(str(val)):
                setattr(self, fn_flag, True)

                
    def _init_globato(self):
        """Attempt to auto-load weight/uncertainty from Globato HDF5."""
        pass 

    
    def __call__(self):
        return self.generate()

    
    def init_region(self, src_ds=None):
        """Initialize region from dataset."""
        
        if src_ds is None and self.src_dem:
            src_ds = gdal.Open(self.src_dem)
            
        if src_ds:
            self.ds_config = gdalfun.gdal_infos(src_ds)
            self.src_region = regions.Region().from_geo_transform(
                self.ds_config['geoT'], self.ds_config['nx'], self.ds_config['ny']
            )
            self.gt = self.ds_config['geoT']
            return self.src_region, self.ds_config
        return None, None

    
    def init_ds(self, src_ds=None):
        """Initialize GDAL dataset properties."""
        
        if src_ds is None: return

        self.ds_config = gdalfun.gdal_infos(src_ds)
        self.ds_band = src_ds.GetRasterBand(self.band)
        self.gt = self.ds_config['geoT']
        
        if self.ds_band.GetNoDataValue() is None:
            self.ds_band.SetNoDataValue(self.ds_config['ndv'])

            
    def load_aux_data(self, aux_source_list=None, src_region=None, src_gt=None, x_count=None, y_count=None):
        """Load auxiliary datasets into aligned array."""
        
        if aux_source_list is None: aux_source_list = self.aux_data
        if not aux_source_list: return None

        if src_region is None: src_region = self.src_region
        if src_gt is None: src_gt = self.gt
        
        if x_count is None: x_count = self.ds_config['nx']
        if y_count is None: y_count = self.ds_config['ny']
        
        x_inc = src_gt[1]
        y_inc = src_gt[5]
        dst_srs = self.ds_config.get('proj')

        combined_arr = np.full((y_count, x_count), np.nan, dtype=np.float32)
        
        for aux_fn in aux_source_list:
            try:
                aux_ds = DatasetFactory(
                    mod=aux_fn,
                    src_region=src_region,
                    x_inc=x_inc,
                    y_inc=y_inc,
                    dst_srs=dst_srs,
                    verbose=self.verbose
                )._acquire_module()
                
                if aux_ds is None: continue
                aux_ds.initialize()

                for arrs, srcwin, gt in aux_ds.yield_array(want_sums=False):
                    if arrs['z'] is None: continue
                    
                    y_start, y_end = srcwin[1], srcwin[1] + srcwin[3]
                    x_start, x_end = srcwin[0], srcwin[0] + srcwin[2]
                    
                    if y_end > combined_arr.shape[0] or x_end > combined_arr.shape[1]:
                        pass 

                    chunk_data = arrs['z']
                    valid_mask = ~np.isnan(chunk_data)
                    
                    target_slice = combined_arr[y_start:y_end, x_start:x_end]
                    
                    if target_slice.shape == chunk_data.shape:
                        target_slice[valid_mask] = chunk_data[valid_mask]

            except Exception as e:
                if self.verbose:
                    utils.echo_warning_msg(f"Failed to load aux data {aux_fn}: {e}")
                    
        return combined_arr

    
    def generate(self):
        """Execute the filter."""
        
        if self.verbose:
            utils.echo_msg(f'Filtering {self.src_dem} using {self.__class__.__name__}')
            
        self.run()
        self.split_by_z()
        self.split_by_weight()
        self.split_by_uncertainty()

        ## Calculate final statistics
        self.calculate_difference()

        if self.verbose:
            utils.echo_msg(f'Filtered {self.src_dem} to {self.dst_dem}')
            
        self.run()
        self.split_by_z()
        self.split_by_weight()
        self.split_by_uncertainty()

        ## Calculate final statistics

        return self        

    
    def run(self):
        raise NotImplementedError("Subclasses must implement run()")


    def calculate_difference(self):
        """Calculate statistics on modification:

        - Removed Cells: Valid -> NoData
        - Changed Cells: Valid -> Valid (New Value)
        """
        
        if not self.src_dem or not os.path.exists(self.src_dem): return {}
        if not self.dst_dem or not os.path.exists(self.dst_dem): return {}

        stats = {}
        try:
            src_ds = gdal.Open(self.src_dem)
            dst_ds = gdal.Open(self.dst_dem)
            
            if not src_ds or not dst_ds: return {}
            
            src_band = src_ds.GetRasterBand(self.band)
            dst_band = dst_ds.GetRasterBand(1)
            
            ## Read Arrays (Warning: Large files read fully into RAM)
            src_arr = src_band.ReadAsArray().astype(np.float32)
            dst_arr = dst_band.ReadAsArray().astype(np.float32)
            
            src_ndv = src_band.GetNoDataValue()
            dst_ndv = dst_band.GetNoDataValue()
            
            ## Normalize to NaN
            if src_ndv is not None: src_arr[src_arr == src_ndv] = np.nan
            if dst_ndv is not None: dst_arr[dst_arr == dst_ndv] = np.nan
            
            ## Identify Valid Sets
            src_valid = ~np.isnan(src_arr)
            dst_valid = ~np.isnan(dst_arr)
            
            total_src_valid = np.count_nonzero(src_valid)
            
            ## Removed Cells (Source Valid -> Dest Invalid)
            removed_mask = src_valid & ~dst_valid
            removed_count = np.count_nonzero(removed_mask)
            
            ## Changed Cells (Source Valid -> Dest Valid AND Value Changed)
            ## Intersection of valid pixels
            intersection = src_valid & dst_valid
            
            changed_count = 0
            if np.count_nonzero(intersection) > 0:
                diff = np.abs(dst_arr[intersection] - src_arr[intersection])
                
                ## Check for significant difference
                changed_mask = diff > 1e-5
                changed_count = np.count_nonzero(changed_mask)
                
                ## Calc stats on CHANGED values only
                if changed_count > 0:
                    real_diffs = dst_arr[intersection][changed_mask] - src_arr[intersection][changed_mask]
                    stats['mean_diff'] = np.mean(real_diffs)
                    stats['min_diff'] = np.min(real_diffs)
                    stats['max_diff'] = np.max(real_diffs)
                    stats['std_diff'] = np.std(real_diffs)

            ## Populate Stats
            total_modified = removed_count + changed_count
            stats['total_valid_source'] = total_src_valid
            stats['removed_cells'] = removed_count
            stats['changed_cells'] = changed_count
            stats['modified_cells'] = total_modified
            stats['percent_modified'] = (total_modified / total_src_valid) * 100 if total_src_valid > 0 else 0
            
            if self.verbose:
                utils.echo_msg("Filter Statistics:")
                utils.echo_msg(f"  Total Valid Cells: {total_src_valid}")
                utils.echo_msg(f"  Removed Cells: {removed_count} ({(removed_count/total_src_valid)*100:.2f}%)" if total_src_valid else "  Removed Cells: 0")
                utils.echo_msg(f"  Changed Cells: {changed_count} ({(changed_count/total_src_valid)*100:.2f}%)" if total_src_valid else "  Changed Cells: 0")
                if changed_count > 0:
                    utils.echo_msg(f"  Difference (Mean +/- Std): {stats.get('mean_diff',0):.4f} +/- {stats.get('std_diff', 0):.4f}")

        except Exception as e:
            if self.verbose:
                utils.echo_warning_msg(f"Failed to calculate stats: {e}")
            
        return stats
    
    
    def copy_src_dem(self):
        with gdalfun.gdal_datasource(self.src_dem, update=False) as src_ds:
            if src_ds:
                info = gdalfun.gdal_infos(src_ds)
                driver = gdal.GetDriverByName(info['fmt'])
                return driver.CreateCopy(self.dst_dem, src_ds, 0, options=['COMPRESS=DEFLATE'])
        return None


    def extract_src_array(self, band=1, buffer_cells=0):
        """load the src dem as an array"""

        if not self.src_region or not self.ds_config:
            self.init_region()

        x_inc = self.ds_config['geoT'][1]
        y_inc = self.ds_config['geoT'][5] * -1

        self.src_region.buffer(
            x_bv=(x_inc * self.buffer_cells),
            y_bv=(y_inc * self.buffer_cells)
        )

        xcount, ycount, dst_gt = self.src_region.geo_transform(
            x_inc=self.ds_config['geoT'][1],
            y_inc=self.ds_config['geoT'][5],
            node='grid'
        )
        
        src_arr = np.full((ycount, xcount), np.nan)

        with gdalfun.gdal_datasource(self.src_dem) as src_ds:
            if src_ds is not None:
                #self.init_ds(src_ds)
                ds_band = src_ds.GetRasterBand(band)
                if ds_band.GetNoDataValue() is None:
                    ds_band.SetNoDataValue(self.ds_config['ndv'])
                    
                srcwin = (self.buffer_cells, self.buffer_cells,
                          src_ds.RasterXSize, src_ds.RasterYSize)

                src_arr[
                    srcwin[1]:srcwin[1] + srcwin[3],
                    srcwin[0]:srcwin[0] + srcwin[2]
                ] = ds_band.ReadAsArray()

                src_arr[src_arr == self.ds_config['ndv']] = np.nan
                return src_arr
        return None
        
    
    def _get_mask_array(self, mask_attr, is_fn_attr, is_band_attr):
        """Load a mask array based on config."""
        
        mask_val = getattr(self, mask_attr)
        if mask_val is None: return None
        
        if getattr(self, is_fn_attr):
            ## Load from external file using aux loader
            return self.load_aux_data(aux_source_list=[mask_val])
        elif getattr(self, is_band_attr):
            ## Load from internal band
            with gdalfun.gdal_datasource(self.src_dem) as ds:
                if ds:
                    return ds.GetRasterBand(mask_val).ReadAsArray()
        return None


    def _apply_split_mask(self, mask_array, min_val, max_val, label="Value"):
        """Apply a masking operation to the destination DEM.
        
        If pixels fall OUTSIDE the valid range (min_val, max_val),
        they are reverted to the ORIGINAL source data values 
        """
        
        if mask_array is None: return

        ## Open Destination (Filtered Result)
        ds = gdal.Open(self.dst_dem, gdal.GA_Update)
        if not ds: return

        band = ds.GetRasterBand(1)
        arr = band.ReadAsArray()
        
        ## Treat NDV in mask as 0
        mask_array[np.isnan(mask_array)] = 0 
        
        ## Determine Pixels to Revert
        filter_mask = np.zeros(arr.shape, dtype=bool)
        
        if min_val is not None:
            filter_mask |= (mask_array < min_val)
        if max_val is not None:
            filter_mask |= (mask_array > max_val)
            
        if np.any(filter_mask):
            if self.verbose:
                count = np.count_nonzero(filter_mask)
                utils.echo_msg(f'Reverting {count} pixels based on {label} thresholds.')
                
            ## Load Original Source Data
            with gdalfun.gdal_datasource(self.src_dem) as src_ds:
                if src_ds:
                    src_band = src_ds.GetRasterBand(self.band)
                    src_arr = src_band.ReadAsArray()
                    
                    ## Apply Reversion
                    ## Only where the filter mask is True, replace filtered data with source data
                    arr[filter_mask] = src_arr[filter_mask]
            
            band.WriteArray(arr)
            
        ds = None
        
        
    def split_by_z(self):
        """Filter dst_dem by Z thresholds."""
        
        if self.max_z is None and self.min_z is None: return
        if self.verbose: utils.echo_msg(f'Splitting by Z: {self.min_z} to {self.max_z}')
        
        ds = gdal.Open(self.dst_dem, gdal.GA_Update)
        if not ds: return

        band = ds.GetRasterBand(1)
        arr = band.ReadAsArray()
        ndv = band.GetNoDataValue()
        
        ## Self-masking (mask array IS the data array)
        self._apply_split_mask(arr, self.min_z, self.max_z, label="Z")

        
    def split_by_weight(self):
        """Filter dst_dem by Weight thresholds."""
        
        if self.max_weight is None and self.min_weight is None: return
        if self.verbose: utils.echo_msg(f'Splitting by Weight: {self.min_weight} to {self.max_weight}')
        
        weight_arr = self._get_mask_array('weight_mask', 'weight_is_fn', 'weight_is_band')
        self._apply_split_mask(weight_arr, self.min_weight, self.max_weight, label="Weight")

        
    def split_by_uncertainty(self):
        """Filter dst_dem by Uncertainty thresholds."""
        
        if self.max_uncertainty is None and self.min_uncertainty is None: return
        if self.verbose: utils.echo_msg(f'Splitting by Uncertainty: {self.min_uncertainty} to {self.max_uncertainty}')
        
        unc_arr = self._get_mask_array('uncertainty_mask', 'uncertainty_is_fn', 'uncertainty_is_band')
        self._apply_split_mask(unc_arr, self.min_uncertainty, self.max_uncertainty, label="Uncertainty")
        

    def get_outliers(self, in_array, percentile=75, k=1.5, verbose=False):
        if np.all(np.isnan(in_array)): return np.nan, np.nan
        p_max = np.nanpercentile(in_array, percentile)
        p_min = np.nanpercentile(in_array, 100 - percentile)
        iqr = (p_max - p_min) * k
        return p_max + iqr, p_min - iqr


    def _density(self, src_arr):
        """Calculate the density of valid data cells.
        
        If a 'count_mask' is configured, density is calculated as:
        (Cells with Count > 0) / (Total Cells).
        
        Otherwise, it falls back to:
        (Non-NaN Cells in src_arr) / (Total Cells).
        
        Args:
            src_arr (np.array): The source data array (elevation) for fallback dimensions.
            
        Returns:
            float: Density value (0.0 to 1.0).
        """
        
        count_arr = None

        ## Try to load Count Mask
        if self.count_mask is not None:
            
            ## External File
            if hasattr(self, 'count_is_fn') and self.count_is_fn:
                ## Use load_aux_data to get aligned array
                count_arr = self.load_aux_data(
                    aux_source_list=[self.count_mask],
                    x_count=src_arr.shape[1],
                    y_count=src_arr.shape[0]
                )
                
            ## Band in Source DEM
            elif hasattr(self, 'count_is_band') and self.count_is_band:
                with gdalfun.gdal_datasource(self.src_dem) as ds:
                    if ds:
                        if not self.src_region: self.init_region(ds)
                        x_count = src_arr.shape[1]
                        y_count = src_arr.shape[0]
                        srcwin = self.src_region.srcwin(self.gt, x_count, y_count, node='grid')
                        
                        if srcwin[2] > 0 and srcwin[3] > 0:
                            band = ds.GetRasterBand(self.count_mask)
                            count_arr = band.ReadAsArray(srcwin[0], srcwin[1], srcwin[2], srcwin[3])

        ## Calculate Density
        if count_arr is not None:
            if count_arr.shape == src_arr.shape:
                # Density of populated cells (count > 0)
                populated = np.count_nonzero(count_arr > 0)
                total = count_arr.size
                return populated / total if total > 0 else 0.0

        ## Fallback: Density of valid pixels in source array
        if src_arr is not None and src_arr.size > 0:
            populated = np.count_nonzero(~np.isnan(src_arr))
            return populated / src_arr.size
            
        return 0.0

    
class GritsFactory(factory.CUDEMFactory):
    from . import blur
    from . import gmtfilter
    from . import lspoutliers
    from . import weights
    from . import flats
    from . import blend
    from . import morphology
    from . import denoise
    from . import fill
    from . import slope_filter
    from . import diff
    from . import zscore
    from . import hydro
    from . import cut
    from . import clip
    from . import rivers
    
    _modules = {
        'blur': {'name': 'blur', 'desc': 'Gaussian Blur', 'call': blur.Blur},
        'grdfilter': {'name': 'grdfilter', 'desc': 'GMT grdfilter', 'call': gmtfilter.GMTgrdfilter},
        'outliers': {'name': 'outliers', 'desc': 'Remove outliers (IQR)', 'call': lspoutliers.LSPOutliers},
        'flats': {'name': 'flats', 'desc': 'Remove flat areas', 'call': flats.Flats},
        'weights': {'name': 'weights', 'desc': 'Weight buffering', 'call': weights.Weights},
        'blend': {'name': 'blend', 'desc': 'Blend aux data', 'call': blend.Blend},
        'morphology': {'name': 'morphology', 'desc': 'Morphological filtering', 'call': morphology.Morphology},
        'denoise': {'name': 'denoise', 'desc': 'Denoising filtering', 'call': denoise.Denoise},
        'fill': {'name': 'fill', 'desc': 'Fill nodata values', 'call': fill.Fill},
        'slope': {'name': 'slope', 'desc': 'Filter by Slope/LSP', 'call': slope_filter.SlopeFilter},
        'diff': {'name': 'diff', 'desc': 'Difference Filter', 'call': diff.Diff},
        'zscore': {'name': 'zscore', 'desc': 'Z-Score Anomaly Filter', 'call': zscore.ZScore},
        'hydro': {'name': 'hydro', 'desc': 'Hydro Enforcement', 'call': hydro.Hydro},
        'cut': {'name': 'cut', 'desc': 'Cut/Mask to region', 'call': cut.Cut},
        'clip': {'name': 'cut', 'desc': 'Clip to vector', 'call': clip.Clip},
        'rivers': {'name': 'rivers', 'desc': 'Discover rivers from ndv', 'call': rivers.Rivers}, 
    }

    def __init__(self, **kwargs):
        super().__init__(**kwargs)


## ==============================================
## Command-line Interface (CLI)
## $ grits
##
## grits cli
## ==============================================        
class PrintModulesAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        factory.echo_modules(GritsFactory._modules, values)
        sys.exit(0)

        
def grits_cli():
    parser = argparse.ArgumentParser(
        description=f'Grits ({__version__}): GRId filTerS - DEM Filtering Tool',
        epilog="Supported %(prog)s modules (see %(prog)s --modules <module-name> for more info):\n" 
        f"{factory.get_module_short_desc(GritsFactory._modules)}\n\n"
        "CUDEM home page: <http://cudem.colorado.edu>",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    ## --- Processing Options ---
    proc_grp = parser.add_argument_group("Processing Options")
    proc_grp.add_argument('dems', nargs='+', help='Input DEM files or Waffles JSON configs')
    proc_grp.add_argument('-M', '--module', action='append', required=True, help='Grits module')
    proc_grp.add_argument('-O', '--outfile', help='Output filename')
    proc_grp.add_argument('-R', '--region', action='append', help='Restrict processing region')
    
    # --- Threshold Options ---
    thresh_grp = parser.add_argument_group("Threshold Options")
    thresh_grp.add_argument('-N', '--min-z', type=float, help='Minimum Z value')
    thresh_grp.add_argument('-X', '--max-z', type=float, help='Maximum Z value')
    thresh_grp.add_argument('-Wn', '--min-weight', type=float, help='Minimum Weight')
    thresh_grp.add_argument('-Wx', '--max-weight', type=float, help='Maximum Weight')
    thresh_grp.add_argument('-Un', '--min-uncertainty', type=float, help='Minimum Uncertainty')
    thresh_grp.add_argument('-Ux', '--max-uncertainty', type=float, help='Maximum Uncertainty')
    
    # --- Mask Otions ---
    mask_grp = parser.add_argument_group("Mask Options")
    mask_grp.add_argument('-U', '--uncertainty-mask', help='Uncertainty raster/band')
    mask_grp.add_argument('-W', '--weight-mask', help='Weight raster/band')
    mask_grp.add_argument('-C', '--count-mask', help='Count raster/band')
    mask_grp.add_argument('-A', '--aux-data', action='append', help='Auxiliary data files')

    ## --- System Options ---
    sys_grp = parser.add_argument_group("System Options")
    sys_grp.add_argument('--version', action='version', version=f'CUDEM {__cudem_version__} :: %(prog)s {__version__}')
    sys_grp.add_argument('-q', '--quiet', action='store_true', help='Quiet mode')

    ## --- Info Helpers (Exit after printing) ---
    sys_grp.add_argument('-m', '--modules', nargs='?', action=PrintModulesAction, help='List available modules')
    
    args = parser.parse_args()

    for src_dem in args.dems:
        current_dem = src_dem
        for mod_str in args.module:
            if args.outfile:
                dst_dem = args.outfile
            else:
                base = os.path.basename(current_dem)
                dst_dem = f"{os.path.splitext(base)[0]}_{mod_str.split(':')[0]}.tif"

            fac = GritsFactory(
                mod=mod_str,
                src_dem=current_dem,
                dst_dem=dst_dem,
                min_z=args.min_z,
                max_z=args.max_z,
                min_weight=args.min_weight,
                max_weight=args.max_weight,
                min_uncertainty=args.min_uncertainty,
                max_uncertainty=args.max_uncertainty,
                uncertainty_mask=args.uncertainty_mask,
                weight_mask=args.weight_mask,
                count_mask=args.count_mask,
                aux_data=args.aux_data,
                verbose=not args.quiet
            )
            
            try:
                filter_obj = fac._acquire_module()
                if filter_obj:
                    res = filter_obj()
                    if res: current_dem = res.dst_dem
                else:
                    utils.echo_error_msg(f"Could not load module: {mod_str}")
            except Exception as e:
                utils.echo_error_msg(f"Error executing {mod_str}: {e}")
                if not args.quiet: traceback.print_exc()

                
if __name__ == '__main__':
    grits_cli()

### End
