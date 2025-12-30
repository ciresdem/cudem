### lspoutliers.py - GRId filTerS
##
## Copyright (c) 2024 - 2025 Regents of the University of Colorado
##
## lspoutliers.py is part of cudem
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
## Filter grids using Land Surface Parameters (LSP) to detect and remove outliers.
## Uses Tukey's fences method on various derived metrics (Slope, TRI, Roughness, etc.)
##
### Code:

import os
import math
import numpy as np
import scipy.ndimage
import scipy.interpolate
from osgeo import gdal

from cudem import utils
from cudem import gdalfun
from cudem.grits import grits

class LSPOutliers(grits.Grits):
    """Discover and remove outliers from the input DEM using LSP metrics.

    Uses derived parameters (Slope, TRI, Roughness, TPI) alongside Elevation
    to calculate outlier scores. Cells exceeding a combined score threshold are masked.

    Parameters:
    -----------
    percentile : float
        The base percentile to define Q3 (e.g. 75). Auto-calculated if None.
    k : float
        The outlier factor (Q3 + k*IQR). Default 1.5.
    multipass : int
        Number of passes at varying scales.
    accumulate : bool
        If True, accumulates scores across passes before masking.
        If False, masks outliers immediately after each pass.
    """
    
    def __init__(self, percentile=None, max_percentile=None, 
                 chunk_size=None, chunk_step=None, max_chunk=None, max_step=None,
                 k=None, max_k=None, return_mask=False,
                 elevation_weight=1, curvature_weight=1, slope_weight=1,
                 tpi_weight=1, unc_weight=1, rough_weight=1, tri_weight=1,
                 count_weight=0, multipass=1, accumulate=False, interpolation='nearest',
                 aggressive=True, units_are_degrees=True, size_is_step=True,
                 mode='scaled', fill_removed_data=False, **kwargs):
        
        super().__init__(**kwargs)
        
        ## Outlier Parameters
        self.percentile = utils.float_or(percentile)
        self.max_percentile = utils.float_or(max_percentile)
        self.k = utils.float_or(k)
        self.max_k = utils.float_or(max_k)
        
        ## Windowing Parameters
        self.chunk_size = chunk_size
        self.chunk_step = chunk_step
        self.max_chunk = max_chunk
        self.max_step = max_step
        self.size_is_step = size_is_step
        
        ## Method Options
        self.multipass = utils.int_or(multipass, 1)
        self.accumulate = accumulate
        self.interpolation = interpolation
        self.return_mask = return_mask
        self.mode = mode
        self.fill_removed_data = fill_removed_data
        
        ## Weights
        self.elevation_weight = utils.float_or(elevation_weight, 0)
        self.slope_weight = utils.float_or(slope_weight, .25)
        self.tpi_weight = utils.float_or(tpi_weight, .5)
        self.unc_weight = utils.float_or(unc_weight, 1)
        self.rough_weight = utils.float_or(rough_weight, .25)
        self.tri_weight = utils.float_or(tri_weight, .5)
        self.curvature_weight = utils.float_or(curvature_weight, 1)
        self.count_weight = utils.float_or(count_weight, 1)
        
        self.units_are_degrees = units_are_degrees
        self.aggressive = aggressive

        
    def init_percentiles(self, src_ds=None):
        """Initialize outlier detection thresholds."""
        
        ## If not set, try to auto-calculate based on roughness
        if self.percentile is None or self.k is None:
            p, k, pp = self.rough_q(src_ds)
            if self.percentile is None:
                self.percentile = p if p is not None else 75
            if self.k is None:
                self.k = k if k is not None else 1.5
                
        ## Defaults
        if self.percentile is None or np.isnan(self.percentile):
            self.percentile = 75
        if self.k is None:
            self.k = 1.5

        ## Max values for multipass
        self.max_percentile = utils.float_or(
            self.max_percentile, ((100 - self.percentile) / 2) + self.percentile
        )
        self.max_k = utils.float_or(self.max_k, self.k**2)

        if self.verbose:
            utils.echo_msg(f"Outlier Config: Perc={self.percentile}-{self.max_percentile}, K={self.k}-{self.max_k}")

            
    def init_chunks(self, src_ds=None):
        """Initialize processing chunk sizes."""
        
        if src_ds is None: return
        
        ## Get pixel size
        gt = src_ds.GetGeoTransform()
        cell_size = gt[1]
        
        ## Scale to meters if degrees (approx)
        if self.units_are_degrees and cell_size < 1:
            cell_size *= 111120 

        ## Default sizes (meters)
        m_size = 1000
        mm_size = 10000
        
        ## Parse inputs (handle 'm' suffix)
        if self.chunk_size and str(self.chunk_size).endswith('m'):
            m_size = utils.int_or(self.chunk_size[:-1])
        elif self.chunk_size:
            m_size = float(self.chunk_size) * cell_size

        if self.max_chunk and str(self.max_chunk).endswith('m'):
            mm_size = utils.int_or(self.max_chunk[:-1])
        elif self.max_chunk:
            mm_size = float(self.max_chunk) * cell_size

        ## Convert back to pixels
        ## Simple heuristic based on density (default 1)
        n_den = 1.0 
        self.n_chunk = int((m_size * (1/n_den)) / cell_size)
        self.max_chunk = int((mm_size * (1/n_den)) / cell_size)
        
        ## Safety limits
        self.n_chunk = max(15, self.n_chunk)
        self.max_chunk = max(self.n_chunk * 2, self.max_chunk)
        
        if self.size_is_step:
            self.n_step = self.n_chunk
            self.max_step = self.max_chunk
        else:
            self.n_step = int(self.n_chunk / 2)
            self.max_step = int(self.max_chunk / 2)

            
    def _generate_mask_ds(self, src_ds):
        """Create temporary mask dataset on disk."""
        
        self.mask_mask_fn = utils.make_temp_fn(f'{utils.fn_basename2(self.src_dem)}_outliers.tif', temp_dir=self.cache_dir)
        
        driver = gdal.GetDriverByName('GTiff')
        if os.path.exists(self.mask_mask_fn):
            driver.Delete(self.mask_mask_fn)
            
        self.mask_mask_ds = driver.Create(
            self.mask_mask_fn,
            src_ds.RasterXSize, src_ds.RasterYSize,
            2, gdal.GDT_Float32,
            options=['COMPRESS=DEFLATE', 'TILED=YES', 'BIGTIFF=YES']
        )
        self.mask_mask_ds.SetGeoTransform(src_ds.GetGeoTransform())
        self.mask_mask_ds.SetProjection(src_ds.GetProjection())
        
        self.mask_mask_band = self.mask_mask_ds.GetRasterBand(1)
        self.mask_count_band = self.mask_mask_ds.GetRasterBand(2)
        
        ## Init with 0
        self.mask_mask_band.Fill(0)
        self.mask_count_band.Fill(0)

        
    def gdal_dem(self, input_ds=None, var=None):
        """Helper to run gdaldem via API."""
        
        if var == 'curvature':
            return self.gdal_dem_curvature(input_ds)
        elif var == 'easterliness':
            return self.gdal_dem_trig(input_ds, np.sin)
        elif var == 'northerliness':
            return self.gdal_dem_trig(input_ds, np.cos)
        
        tmp_fn = utils.make_temp_fn(f'gdaldem_{var}.tif', self.cache_dir)
        tmp_ds = gdal.DEMProcessing(
            tmp_fn, input_ds, var, computeEdges=True, scale=111120 if self.units_are_degrees else 1.0
        )
        return tmp_ds, tmp_fn

    
    def gdal_dem_curvature(self, input_ds):
        slp_ds, slp_fn = self.gdal_dem(input_ds, 'slope')
        curv_ds, curv_fn = self.gdal_dem(slp_ds, 'slope')
        utils.remove_glob(slp_fn)
        return curv_ds, curv_fn

    
    def gdal_dem_trig(self, input_ds, func):
        asp_ds, asp_fn = self.gdal_dem(input_ds, 'aspect')
        with gdalfun.gdal_datasource(asp_ds, update=True) as ds:
            b = ds.GetRasterBand(1)
            a = b.ReadAsArray()
            ## Aspect is degrees, convert to rad
            res = func(np.radians(a))
            b.WriteArray(res)
        return asp_ds, asp_fn

    
    def get_pk(self, src_ds, var='roughness'):
        """Calculate percentile (p) and outlier factor (k) stats from a derived surface."""
        
        if var == 'elevation':
            ds_ds, ds_fn = src_ds, None
        else:
            ds_ds, ds_fn = self.gdal_dem(src_ds, var)
            
        ## Use TPI of the derived surface as a roughness proxy
        pk_ds, pk_fn = self.gdal_dem(ds_ds, 'TPI')
        pk_arr = pk_ds.GetRasterBand(1).ReadAsArray().astype(float)
        
        ## Clean
        pk_arr[pk_arr == -9999] = np.nan
        pk_arr[np.isinf(pk_arr)] = np.nan
        
        if np.all(np.isnan(pk_arr)):
            utils.remove_glob(pk_fn, ds_fn)
            return None, None, None

        std_pk = np.nanstd(pk_arr)
        pk_range = np.nanmax(pk_arr) - np.nanmin(pk_arr)
        
        if pk_range == 0:
            utils.remove_glob(pk_fn, ds_fn)
            return 75, 1.5, 50

        ## Heuristic for k and p
        pkrm = (std_pk - np.nanmin(pk_arr)) / pk_range
        k = (1 - pkrm) * 4.5
        p = (pkrm * (.5 - 1) + 1) * 100
        pp = pkrm * 100
        
        utils.remove_glob(pk_fn, ds_fn)
        return pp, k, p

    
    def rough_q(self, src_ds):
        """Aggregate roughness stats from multiple metrics."""
        
        vars = ['roughness', 'slope', 'TPI', 'elevation'] # curvature omitted for speed
        ks = []
        ps = []
        pps = []
        
        for v in vars:
            pp, k, p = self.get_pk(src_ds, v)
            if k is not None:
                ks.append(k)
                ps.append(p)
                pps.append(pp)
                
        if not ks: return 75, 1.5, 50
        
        return np.mean(ps), np.mean(ks), np.mean(pps)

    
    def mask_outliers(self, src_data, mask_data, count_data, percentile=75, 
                      upper_only=False, src_weight=1, k=1.5):
        """Identify outliers in src_data and accumulate scores into mask_data/count_data."""
        
        ## Clean data
        src_data[np.isinf(src_data)] = np.nan
        
        ## Calculate limits
        upper_lim, lower_lim = self.get_outliers(src_data, percentile, k)
        
        ## Process Upper
        upper_mask = src_data > upper_lim
        if np.any(upper_mask):
            score = 0
            if self.mode == 'average':
                ## Linear score distance from limit
                m = np.nanmax(src_data[upper_mask])
                if m != upper_lim:
                    score = np.abs((src_data[upper_mask] - upper_lim) / (m - upper_lim))
            elif self.mode == 'scaled':
                ## Root sum square accumulation
                pass # Handled below
                
            ## Generic Score Accumulation (Scaled Mode)
            m = np.nanmax(src_data[upper_mask])
            diff = m - upper_lim
            val = np.zeros_like(src_data[upper_mask])
            if diff != 0:
                val = np.abs((src_data[upper_mask] - upper_lim) / diff)
                
            ## Accumulate
            # mask_data[upper_mask] = sqrt(current^2 + (weight*score)^2)
            current = mask_data[upper_mask]
            mask_data[upper_mask] = np.sqrt(current**2 + (src_weight * val)**2)
            count_data[upper_mask] += src_weight

        ## Process Lower
        if not upper_only:
            lower_mask = src_data < lower_lim
            if np.any(lower_mask):
                m = np.nanmin(src_data[lower_mask])
                diff = m - lower_lim
                val = np.zeros_like(src_data[lower_mask])
                if diff != 0:
                    val = np.abs((src_data[lower_mask] - lower_lim) / diff)
                    
                current = mask_data[lower_mask]
                mask_data[lower_mask] = np.sqrt(current**2 + (src_weight * val)**2)
                count_data[lower_mask] += src_weight

                
    def apply_mask(self, perc=75):
        """Apply the accumulated outlier mask to the destination DEM."""
        
        ## Read accumulators
        mask_arr = self.mask_mask_band.ReadAsArray()
        count_arr = self.mask_count_band.ReadAsArray()
        
        ## Filter NaNs
        mask_arr[np.isnan(mask_arr)] = 0
        count_arr[np.isnan(count_arr)] = 0
        
        if np.all(mask_arr == 0):
            if self.verbose: utils.echo_msg("No outliers found.")
            return 0

        ## Calculate Thresholds on the Score Arrays themselves
        ## (Outliers of the Outlier Scores)
        count_upper, _ = self.get_outliers(count_arr, perc, self.k)
        mask_upper, _ = self.get_outliers(mask_arr, perc, self.k)
        
        ## Final Decision Mask
        final_mask = (mask_arr > mask_upper) & (count_arr > count_upper)
        
        removed_count = np.count_nonzero(final_mask)
        
        ## Apply to Destination
        ## Open dest (copy of source)
        ds = gdal.Open(self.dst_dem, gdal.GA_Update)
        band = ds.GetRasterBand(self.band)
        data = band.ReadAsArray()
        ndv = self.ds_config['ndv']
        
        if self.fill_removed_data:
            ## Interpolate holes
            data[final_mask] = np.nan
            ## Simple hole fill logic (e.g. gdal_fillnodata or scipy)
            ## For now, just set to NDV as filling is complex without context
            data[final_mask] = ndv 
        else:
            data[final_mask] = ndv
            
        band.WriteArray(data)
        ds = None
        
        if self.verbose:
            utils.echo_msg(f"Removed {removed_count} outliers.")
            
        return removed_count

    
    def run(self):
        """Execute the outlier filter."""
        
        ## Setup
        dst_ds = self.copy_src_dem() # self.dst_dem
        if not dst_ds: return self.src_dem, -1
        
        with gdalfun.gdal_datasource(self.src_dem) as src_ds:
            if not src_ds: return self.src_dem, -1
            
            self.init_ds(src_ds)
            self.init_chunks(src_ds)
            self.init_percentiles(src_ds)
            
            ## Setup Accumulation Mask
            self._generate_mask_ds(src_ds)
            
            ## Iterate Passes (Multipass)
            ## Generate stepping arrays
            chunks = np.linspace(self.n_chunk, self.max_chunk, self.multipass).astype(int)
            percs = np.linspace(self.percentile, self.max_percentile, self.multipass)
            ks = np.linspace(self.k, self.max_k, self.multipass)
            
            for i in range(self.multipass):
                chunk = chunks[i]
                perc = percs[i]
                k = ks[i]
                
                ## Weight attenuation (later passes weigh less)
                w_scale = 1.0 / (i + 1)
                
                if self.verbose:
                    utils.echo_msg(f"Pass {i+1}/{self.multipass}: Chunk={chunk}, Perc={perc:.1f}, K={k:.2f}")

                ## If not accumulating, reset mask each pass
                if not self.accumulate:
                    self.mask_mask_band.Fill(0)
                    self.mask_count_band.Fill(0)

                ## Process Chunks
                for srcwin in utils.yield_srcwin(
                    (src_ds.RasterYSize, src_ds.RasterXSize), 
                    n_chunk=chunk, step=chunk, 
                    verbose=self.verbose, 
                    msg=f'Scanning outliers'
                ):
                    ## Read Data
                    band_data = self.ds_band.ReadAsArray(*srcwin)
                    band_data[band_data == self.ds_config['ndv']] = np.nan
                    
                    if np.all(np.isnan(band_data)): continue
                    
                    ## Read Mask Accumulators
                    m_mask = self.mask_mask_band.ReadAsArray(*srcwin)
                    m_count = self.mask_count_band.ReadAsArray(*srcwin)
                    
                    ## Create Memory DS for GDAL processing
                    mem_ds = gdalfun.generate_mem_ds(self.ds_config, band_data, srcwin)
                    
                    ## --- Compute LSPs and Mask ---
                    
                    ## Elevation
                    self.mask_outliers(band_data, m_mask, m_count, perc, False, self.elevation_weight * w_scale, k)
                    
                    ## Slope
                    ds, fn = self.gdal_dem(mem_ds, 'slope')
                    d = ds.GetRasterBand(1).ReadAsArray()
                    self.mask_outliers(d, m_mask, m_count, perc, False, self.slope_weight * w_scale, k)
                    utils.remove_glob(fn)
                    
                    ## Roughness
                    ds, fn = self.gdal_dem(mem_ds, 'roughness')
                    d = ds.GetRasterBand(1).ReadAsArray()
                    self.mask_outliers(d, m_mask, m_count, perc, False, self.rough_weight * w_scale, k)
                    utils.remove_glob(fn)
                    
                    ## TRI
                    ds, fn = self.gdal_dem(mem_ds, 'TRI')
                    d = ds.GetRasterBand(1).ReadAsArray()
                    self.mask_outliers(d, m_mask, m_count, perc, False, self.tri_weight * w_scale, k)
                    utils.remove_glob(fn)
                    
                    ## TPI (on Slope)
                    ## Often useful to detect ridges/valleys in slope map
                    slp_ds, slp_fn = self.gdal_dem(mem_ds, 'slope')
                    ds, fn = self.gdal_dem(slp_ds, 'TPI')
                    d = ds.GetRasterBand(1).ReadAsArray()
                    self.mask_outliers(d, m_mask, m_count, perc, False, self.tpi_weight * w_scale, k)
                    utils.remove_glob(fn, slp_fn)

                    ## Write Accumulators Back
                    self.mask_mask_band.WriteArray(m_mask, srcwin[0], srcwin[1])
                    self.mask_count_band.WriteArray(m_count, srcwin[0], srcwin[1])

                ## Apply Mask (If not accumulating for final pass only)
                if not self.accumulate:
                    self.apply_mask(perc)

            ## Final Apply (If accumulating)
            if self.accumulate:
                self.apply_mask(self.percentile)

        ## Cleanup
        self.mask_mask_ds = None
        if not self.return_mask:
            utils.remove_glob(self.mask_mask_fn)
            
        return self.dst_dem, 0

### End

# ### lspoutliers.py - GRId filTerS
# ##
# ## Copyright (c) 2024 - 2025 Regents of the University of Colorado
# ##
# ## Permission is hereby granted, free of charge, to any person obtaining a copy 
# ## of this software and associated documentation files (the "Software"), to deal 
# ## in the Software without restriction, including without limitation the rights 
# ## to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies 
# ## of the Software, and to permit persons to whom the Software is furnished to do so, 
# ## subject to the following conditions:
# ##
# ## The above copyright notice and this permission notice shall be included in all copies or
# ## substantial portions of the Software.
# ##
# ## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
# ## INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR 
# ## PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE 
# ## FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
# ## ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# ## SOFTWARE.
# ##
# ### Commentary:
# ##
# ##
# ### Code:

# import os
# import math
# import numpy as np
# import scipy
# from osgeo import gdal

# from cudem import utils
# from cudem import gdalfun
# from cudem.grits import grits

# class LSPOutliers(grits.Grits):
#     """Discover and remove outliers from the input DEM.

#     Uses various LSPs as well as elevation and uncertainty to discover 
#     possible outliers and remove data cells where there are outliers 
#     of calculated outliers.

#     Outliers are calculated using the `Tukey's fences` method, 
#     where outliers are outside the range:
#     [Q1 - k(Q3 - Q1), Q3 + k(Q3 - Q1)]

#     In general, k=1.5 for general outliers and k=3 for 'far out' outliers. 
#     Since this is fairly arbitrary, k can be set to any positive number, 
#     when auto-generated it will fall between 0 and 4.5

#     Q1 and Q3 above are generally defined as the 25th and 75th percentile, 
#     respectively.

#     If `percentile` is not set, we will calculate a reasonable value to 
#     use in both the gathering of the sub-regional LSP outliers as well as 
#     the final removal of the outliers of the outliers.

#     if `chunk_size` is not set, we calculate the inital moving window size 
#     to be about 5000 meters squared, plus or minus depending on the density 
#     of the input raster.

#     Using multipass, we can scan the data at multiple scales to get a better 
#     representation of possible features in the data. multipass will increase 
#     the moving window size and outlier percentile on each pass from `chunk_size` 
#     to `max_chunk` and `percentile` to `max_percentile`, respectively. By default, 
#     multipass will pass through the data at varying scales and  continually 
#     accumulate outlier scores and calculate the final outliers to remove based 
#     on the accumulated outlier mask. 

#     If `accumulate` is set to `False`, we will instead calculate and remove the 
#     final outliers at each pass.
    
#     Most datasets should be fine with using only the default parameters, however, 
#     closer refinement of the parameters may yield more acceptable results.

#     Parameters:

#     percentile(float) - the percentile to use to define Q3
#     k(float) - the outlier factor ( Q1 - k(iqr), Q3 + k(iqr) )
#     chunk_size(int) - the moving window size in pixels (append `m` to set the size in meters)
#     chunk_step(int) - the moving window step in pixels (append `m` to set the size in meters)
#     interpolation(str) - interpolation method to use for neighborhood calculations 
#                          (None, `linear`, `cubic` or `nearest`)
#     multipass(int) - pass through the data at multiple `chunk_size/step`s, `percentiles` and `k`s, 
#                      set the number of passes here
#     max_percentile(float) - the maximum percentile to use to define Q1 and Q3 (multipass)
#     max_k(float) - the maximum outlier factor ( Q1 - k(iqr), Q3 + k(iqr) ) (multipass)
#     max_chunk(int) - the maximum moving window size in pixels (multipass) 
#                      (append `m` to set the size in meters)
#     max_step(int) - the maximum moving window step in pixels (multipass) 
#                     (append `m` to set the size in meters)
#     acuumulate(bool) - accumulate outliers with multipass, set to `False` to remove 
#                        outliers after each pass in multipass
#     return_mask(bool) - save the generated outlier mask
#     size_is_step(bool) - the chunk_step and max_step will be set to the chunk_size 
#                          and chunk_step, respectively
#     mode(str) - the mode to use when calculating the outliers (average, scaled, None)
#     """
    
#     def __init__(self, percentile=None, max_percentile=None, outlier_percenitle=None,
#                  chunk_size=None, chunk_step=None, max_chunk=None, max_step=None,
#                  k=None, max_k=None, outlier_k=None, return_mask=False,
#                  elevation_weight=1, curvature_weight=1, slope_weight=1,
#                  tpi_weight=1, unc_weight=1, rough_weight=1, tri_weight=1,
#                  count_weight=0, multipass=1, accumulate=False, interpolation='nearest',
#                  aggressive=True, units_are_degrees=True, size_is_step=True,
#                  mode='scaled', fill_removed_data=False, **kwargs):
        
#         super().__init__(**kwargs)
#         self.percentile = utils.float_or(percentile)
#         self.max_percentile = utils.float_or(max_percentile)
#         self.chunk_size = chunk_size
#         self.chunk_step = chunk_step
#         self.max_chunk = max_chunk
#         self.max_step = max_step
#         self.return_mask = return_mask
#         self.interpolation = interpolation
#         self.multipass = utils.int_or(multipass, 1)
#         self.accumulate = accumulate
#         self.elevation_weight = utils.float_or(elevation_weight, 0)
#         self.curvature_weight = utils.float_or(curvature_weight, 1)
#         self.slope_weight = utils.float_or(slope_weight, .25)
#         self.tpi_weight = utils.float_or(tpi_weight, .5)
#         self.unc_weight = utils.float_or(unc_weight, 1)
#         self.count_weight = utils.float_or(count_weight, 1)
#         self.rough_weight = utils.float_or(rough_weight, .25)
#         self.tri_weight = utils.float_or(tri_weight, .5)
#         self.k = utils.float_or(k)
#         self.max_k = utils.float_or(max_k)
#         self.aggressive = aggressive
#         self.units_are_degrees = units_are_degrees
#         self.size_is_step = size_is_step
#         self.mode = mode
#         self.fill_removed_data = fill_removed_data

        
#     def init_percentiles(self, src_ds=None):
#         if self.percentile is None or self.k is None:
#             p, k, pp = self.rough_q(src_ds)
#             if self.percentile is None:
#                 self.percentile = p

#             utils.echo_msg_bold(self.percentile)
#             if self.k is None:
#                 if np.isnan(k):
#                     k = 1.5
#                 else:
#                     self.k = k
#             utils.echo_msg_bold(self.k)
            
#         if self.percentile is None or np.isnan(self.percentile):
#             self.percentile = 75
        
#         self.max_percentile = utils.float_or(
#             self.max_percentile, ((100 - self.percentile) / 2) + self.percentile
#         )
#         if self.k is None:
#             self.k = 1.5

#         self.max_k = utils.float_or(self.max_k, self.k**2)
#         if self.verbose:
#             utils.echo_msg(
#                 'outlier percentiles: {} {} < {} {}'.format(
#                     self.percentile, self.k, self.max_percentile, self.max_k
#                 )
#             )
#             utils.echo_msg(
#                 '[Q1 - k(iqr), Q3 + k(iqr)]; Q1:{q1} Q3:{q3} k:{k}'.format(
#                     q1=100-self.percentile, q3=self.percentile, k=self.k
#                 )
#             )

            
#     def init_chunks(self, src_ds=None):
#         src_arr, src_den = self.gdal_density(src_ds)
#         self._chunks(src_arr)
#         src_arr = src_den = None

        
#     def _chunks(self, src_arr):
#         n_den = self._density(src_arr)
#         cell_size = self.ds_config['geoT'][1]
#         if self.units_are_degrees:
#             # scale cellsize to meters,
#             # todo: check if input is degress/meters/feet
#             cell_size *= 111120 

#         m_size = 1000
#         mm_size = 10000
#         #m_size = (src_arr.shape[0] / n_den) / 24#800#500 # 1000
#         #mm_size = (src_arr.shape[1] / n_den) / 12#8000 # 10000
#         if self.chunk_size is not None:
#             if self.chunk_size[-1] == 'm':
#                 m_size = utils.int_or(self.chunk_size[:-1])
#                 c_size = m_size / cell_size

#         if self.chunk_step is not None:
#             if self.chunk_step[-1] == 'm':
#                 m_step = utils.int_or(self.chunk_step[:-1])
#                 c_step = m_step / cell_size
                
#         if self.max_chunk is not None:
#             if self.max_chunk[-1] == 'm':
#                 mm_size = utils.int_or(self.max_chunk[:-1])
#                 cm_size = mm_size / cell_size

#         if self.max_step is not None:
#             if self.max_step[-1] == 'm':
#                 mm_step = utils.int_or(self.max_step[:-1])
#                 cm_step = mm_step / cell_size

#         self.n_chunk = utils.int_or(
#             self.chunk_size, ((m_size * (1/n_den)) / cell_size)
#         )
#         self.max_chunk = utils.int_or(
#             self.max_chunk, ((mm_size * (1/n_den)) / cell_size)
#         )
#         if self.n_chunk < 15:
#             self.n_chunk = 15
#         if self.max_chunk < 15:
#             self.max_chunk = self.n_chunk * 2
            
#         if self.size_is_step:
#             self.n_step = utils.int_or(self.chunk_step, self.n_chunk)
#             self.max_step = utils.int_or(self.max_step, self.max_chunk)
#         else:
#             self.n_step = utils.int_or(self.chunk_step, math.ceil(self.n_chunk / 2))
#             self.max_step = utils.int_or(self.max_step, math.ceil(self.max_chunk / 2))
#             if self.max_step > self.max_chunk:
#                 self.max_step = self.max_chunk

#         if self.verbose:
#             utils.echo_msg(
#                 'outlier chunks ({} {}): {} {} < {} {}'.format(
#                     n_den, src_arr.shape, self.n_chunk, self.n_step,
#                     self.max_chunk, self.max_step
#                 )
#             )

            
#     def _density(self, src_arr):
#         nonzero = np.count_nonzero(~np.isnan(src_arr))
#         dd = nonzero / src_arr.size
        
#         return(dd)

    
#     def gdal_density(self, src_ds=None):
#         src_arr, src_config = gdalfun.gdal_get_array(src_ds)
#         src_arr[src_arr == src_config['ndv']] = np.nan

#         return(src_arr, self._density(src_arr))

    
#     def _generate_mask_ds(self, src_ds=None):
#         ## to hold the mask data
#         self.mask_mask_fn = '{}{}'.format(utils.fn_basename2(self.src_dem), '_outliers.tif')
#         mask_mask = np.zeros((src_ds.RasterYSize, src_ds.RasterXSize))
#         driver = gdal.GetDriverByName('GTiff')

#         if os.path.exists(self.mask_mask_fn):
#             status = driver.Delete(self.mask_mask_fn)
#             if status != 0:
#                 utils.remove_glob('{}*'.format(self.mask_mask_fn))
        
#         self.mask_mask_ds = driver.Create(
#             self.mask_mask_fn,
#             self.ds_config['nx'],
#             self.ds_config['ny'],
#             2,
#             gdal.GDT_Float32,
#             options=['COMPRESS=DEFLATE',
#                      'PREDICTOR=1',
#                      'TILED=YES',
#                      'BIGTIFF=YES']
#         )
#         self.mask_mask_ds.SetGeoTransform(self.ds_config['geoT'])
#         self.mask_mask_band = self.mask_mask_ds.GetRasterBand(1)
#         self.mask_count_band = self.mask_mask_ds.GetRasterBand(2)
        
#         self.mask_mask_band.SetNoDataValue(0)        
#         self.mask_mask_band.WriteArray(mask_mask)
#         self.mask_count_band.WriteArray(mask_mask)
#         mask_mask = None

        
#     def generate_mask_ds(self, src_ds=None):
#         ## to hold the mask data
#         self.mask_mask_fn = '{}{}'.format(
#             utils.fn_basename2(self.src_dem), '_outliers.tif'
#         )
#         mask_mask = np.zeros((src_ds.RasterYSize, src_ds.RasterXSize))
#         mask_count = np.zeros((src_ds.RasterYSize, src_ds.RasterXSize))
#         driver = gdal.GetDriverByName('GTiff')

#         if os.path.exists(self.mask_mask_fn):
#             status = driver.Delete(self.mask_mask_fn)
#             if status != 0:
#                 utils.remove_glob('{}*'.format(self.mask_mask_fn))
        
#         self.mask_mask_ds = driver.Create(
#             self.mask_mask_fn,
#             self.ds_config['nx'],
#             self.ds_config['ny'],
#             2,
#             gdal.GDT_Float32,
#             options=['COMPRESS=DEFLATE',
#                      'PREDICTOR=1',
#                      'TILED=YES',
#                      'BIGTIFF=YES']
#         )
#         self.mask_mask_ds.SetGeoTransform(self.ds_config['geoT'])
#         self.mask_mask_band = self.mask_mask_ds.GetRasterBand(1)
#         self.mask_count_band = self.mask_mask_ds.GetRasterBand(2)

#         # uncertainty ds
#         unc_band = None
#         if self.uncertainty_mask is not None:
#             if self.unc_is_fn:
#                 unc_ds = gdal.Open(self.uncertainty_mask)
#                 unc_band = unc_ds.GetRasterBand(1)
#             elif self.unc_is_band:
#                 unc_band = src_ds.GetRasterBand(self.uncertainty_mask)

#             if unc_band is not None:
#                 unc_data = unc_band.ReadAsArray()
#                 unc_data[(unc_data == self.ds_config['ndv'])] = 0
#                 self.mask_outliers(
#                     src_data=unc_data,
#                     mask_data=mask_mask,
#                     count_data=mask_count,
#                     percentile=self.percentile,
#                     upper_only=True,
#                     k=self.k,
#                     src_weight=self.unc_weight
#                 )
#                 unc_data = None

#         # # uncertainty ds
#         # cnt_band = None
#         # if self.count_mask is not None:
#         #     if self.cnt_is_fn:
#         #         cnt_ds = gdal.Open(self.count_mask)
#         #         cnt_band = cnt_ds.GetRasterBand(1)
#         #     elif self.cnt_is_band:
#         #         cnt_band = src_ds.GetRasterBand(self.count_mask)

#         #     if cnt_band is not None:
#         #         cnt_data = cnt_band.ReadAsArray()
#         #         #cnt_data[(cnt_data == self.ds_config['ndv'])] = np.nan
#         #         tmp_ds = self.generate_mem_ds(band_data=cnt_data, srcwin=None, return_array=False)
#         #         self.mask_gdal_dem_outliers(srcwin_ds=tmp_ds, band_data=cnt_data, mask_mask_data=mask_mask,
#         #                                     mask_count_data=mask_count, percentile=self.percentile,
#         #                                     upper_only=True, src_weight=self.count_weight, var='TPI', k=self.k)
#         #         # self.mask_outliers(
#         #         #     src_data=cnt_data,
#         #         #     mask_data=mask_mask,
#         #         #     count_data=mask_count,
#         #         #     percentile=self.percentile,
#         #         #     upper_only=True,
#         #         #     k=self.k,
#         #         #     src_weight=self.count_weight
#         #         # )
#         #         cnt_data = tmp_ds = None
                
#         self.mask_mask_band.SetNoDataValue(0)        
#         self.mask_mask_band.WriteArray(mask_mask)
#         self.mask_count_band.WriteArray(mask_count)
#         mask_mask = mask_count = None

        
#     def generate_mem_ds(self, band_data=None, srcwin=None, return_array=False):
#         tmp_band_data = band_data
#         if np.all(np.isnan(tmp_band_data)):
#             return(None)

#         ## interpolate the srcwin for neighborhood calculations
#         if self.interpolation is not None \
#            and self.interpolation in ['nearest', 'linear', 'cubic']:
#             if np.any(np.isnan(band_data)):                        
#                 point_indices = np.nonzero(~np.isnan(tmp_band_data))
#                 if len(point_indices[0]):
#                     point_values = tmp_band_data[point_indices]
#                     xi, yi = np.mgrid[0:srcwin[3], 0:srcwin[2]]

#                     try:
#                         tmp_band_data = scipy.interpolate.griddata(
#                             np.transpose(point_indices), point_values,
#                             (xi, yi), method=self.interpolation
#                         )
#                     except:
#                         pass

#                     point_values = xi = yi = None
#                 point_indices = None

#         if return_array:
#             return(tmp_band_data)
#         else:                
#             ## generate a mem datasource to feed into gdal.DEMProcessing
#             dst_gt = (self.gt[0] + (srcwin[0] * self.gt[1]),
#                       self.gt[1],
#                       0.,
#                       self.gt[3] + (srcwin[1] * self.gt[5]),
#                       0.,
#                       self.gt[5])
#             srcwin_config = gdalfun.gdal_set_infos(
#                 srcwin[2],
#                 srcwin[3],
#                 srcwin[2]*srcwin[3],
#                 dst_gt,
#                 self.ds_config['proj'],
#                 self.ds_band.DataType,
#                 self.ds_config['ndv'],
#                 'GTiff',
#                 {},
#                 1
#             )
#             srcwin_ds = gdalfun.gdal_mem_ds(
#                 srcwin_config, name='MEM', bands=1, src_srs=None
#             )
#             srcwin_band = srcwin_ds.GetRasterBand(1)
#             srcwin_band.SetNoDataValue(self.ds_config['ndv'])
#             tmp_band_data[np.isnan(tmp_band_data)] = self.ds_config['ndv']
#             srcwin_band.WriteArray(tmp_band_data)
#             srcwin_ds.FlushCache()
#             tmp_band_data = None

#             return(srcwin_ds)

        
#     def gdal_dem(self, input_ds=None, var=None):
#         """use gdal to generate various LSPs"""
        
#         if var == 'curvature':
#             return(self.gdal_dem_curvature(input_ds=input_ds))
#         elif var == 'easterliness':
#             return(self.gdal_dem_easterliness(input_ds=input_ds))
#         elif var == 'northerliness':
#             return(self.gdal_dem_northerliness(input_ds=input_ds))
        
#         tmp_ = utils.make_temp_fn(
#             f'gdaldem_{var}.tif', self.cache_dir
#         )
#         tmp_ds = gdal.DEMProcessing(
#             tmp_, input_ds, var, computeEdges=True, scale=111120
#         )

#         return(tmp_ds, tmp_)

    
#     def gdal_dem_curvature(self, input_ds=None):
#         slp_ds, slp_fn = self.gdal_dem(input_ds=input_ds, var='slope')
#         curv_ds, curv_fn = self.gdal_dem(input_ds=slp_ds, var='slope')
#         slp_ds = None
#         utils.remove_glob(slp_fn)
        
#         return(curv_ds, curv_fn)

    
#     def gdal_dem_easterliness(self, input_ds=None):
#         ee_ds, ee_fn = self.gdal_dem(input_ds=input_ds, var='aspect')
#         with gdalfun.gdal_datasource(ee_ds, update=True) as src_ds:
#             b = src_ds.GetRasterBand(1)
#             a = b.ReadAsArray()
#             aa = a * (np.pi/180)
#             ee = np.sin(aa)
#             b.WriteArray(ee)
#             b.FlushCache()
        
#         return(ee_ds, ee_fn)

    
#     def gdal_dem_northerliness(self, input_ds=None):
#         nn_ds, nn_fn = self.gdal_dem(input_ds=input_ds, var='aspect')
#         with gdalfun.gdal_datasource(nn_ds, update=True) as src_ds:
#             b = src_ds.GetRasterBand(1)
#             a = b.ReadAsArray()
#             aa = a * (np.pi/180)
#             nn = np.cos(aa)
#             b.WriteArray(nn)
#             b.FlushCache()
        
#         return(nn_ds, nn_fn)

    
#     def mask_outliers(
#             self, src_data=None, mask_data=None, count_data=None,
#             percentile=75, upper_only=False, src_weight=1, k=1.5,
#             verbose=False, other_data=None, mask_clumps=False
#     ):
#         """mask outliers and assign an outlier 'score' to each affected cell.

#         outlier scores are normalized to UL*2 and LL*2 and combined with 
#         existing outlier scores using the root sum squared.

#         outlier counts are adjusted to the respective `src_weight`.

#         `percentile` and `k` can be adjusted to modify the outlier calculations.

#         Set `upper_only` to True to only mask data which falls above the UL.
#         """
#         #mask_clumps=True
#         if src_data is not None \
#            and mask_data is not None \
#            and count_data is not None:
#             src_data[((np.isnan(src_data)) \
#                       | (np.isinf(src_data)) \
#                       | (src_data == self.ds_config['ndv']))] = np.nan
#             upper_limit, lower_limit = self.get_outliers(
#                 src_data, percentile, k
#             )
#             # if np.isnan(upper_limit) or np.isnan(lower_limit):
#             #upper_limit = np.nanpercentile(src_data, percentile)
#             #lower_limit = np.nanpercentile(src_data, 100-percentile)
                
#             src_size = src_data.size
                
#             if verbose:
#                 utils.echo_msg(f'{upper_limit} {lower_limit}')

#             src_upper = src_data[src_data > upper_limit]
#             src_upper_mask = (src_data > upper_limit)
#             if src_upper.size > 0:
#                 if (upper_limit - lower_limit) != 0:
#                     # if mask_clumps:
#                     l, n = scipy.ndimage.label(src_upper_mask)
#                     uv, uv_counts = np.unique(l[l!=0], return_counts=True)
#                     _size_threshold = src_size * .001
#                     uv_ = uv[uv_counts > _size_threshold]
#                     mask = np.isin(l, uv_)
#                     src_upper_mask[mask] = False
                        
#                     try:
#                         #     src_max = np.nanmax(src_data[src_upper_mask])
#                         if self.mode == 'average':
#                             src_max = np.nanmax(src_data[src_upper_mask])
#                             mask_data[src_upper_mask] += (src_weight * np.abs((src_data[src_upper_mask] - upper_limit) / (src_max - upper_limit)))
#                             #mask_data[src_upper_mask] += (src_weight * np.abs((src_data[src_upper_mask] - lower_limit) / (upper_limit - lower_limit)))
#                             count_data[src_upper_mask] += 1#src_weight
#                         elif self.mode == 'scaled':
#                             src_max = np.nanmax(src_data[src_upper_mask])
#                             mask_data[src_upper_mask] = np.sqrt(
#                                 (np.power(mask_data[src_upper_mask], 2) +
#                                  np.power(src_weight * np.abs((src_data[src_upper_mask] - upper_limit) / (src_max - upper_limit)), 2))
#                             )
#                             count_data[src_upper_mask] += src_weight
#                         else:                            
#                             mask_data[src_upper_mask] = np.sqrt(
#                                 (np.power(mask_data[src_upper_mask], 2) +
#                                  np.power(src_weight * np.abs((src_data[src_upper_mask] - lower_limit) / (upper_limit - lower_limit)), 2))
#                             )
#                             count_data[src_upper_mask] += src_weight

#                         # mask_data[(src_data > upper_limit)] = np.sqrt(
#                         #     (np.power(mask_data[(src_data > upper_limit)], 2) +
#                         #      np.power(src_weight * np.abs((src_upper - lower_limit) / (upper_limit - lower_limit)), 2))
#                         # )
#                         # count_data[(src_data > upper_limit)] += src_weight
#                     except ValueError as e:
#                         pass
#                     except Exception as e:
#                         raise(e)

#             if not upper_only:
#                 src_lower_mask = src_data < lower_limit
#                 src_lower = src_data[src_data < lower_limit]
#                 if src_lower.size > 0:
#                     if (lower_limit - upper_limit) != 0:
#                         # if mask_clumps:
#                         l, n = scipy.ndimage.label(src_lower_mask)
#                         uv, uv_counts = np.unique(l[l!=0], return_counts=True)
#                         _size_threshold = src_size * .001
#                         uv_ = uv[uv_counts > _size_threshold]              
#                         mask = np.isin(l, uv_)
#                         src_lower_mask[mask] = False

#                         try:
#                             #     src_min = np.nanmin(src_data[src_lower_mask])
#                             if self.mode == 'average':
#                                 src_min = np.nanmin(src_data[src_lower_mask])
#                                 mask_data[src_lower_mask] += (src_weight * np.abs((src_data[src_lower_mask] - lower_limit) / (src_min - lower_limit)))
#                                 #mask_data[src_lower_mask] += (src_weight * np.abs((src_data[src_lower_mask] - upper_limit) / (lower_limit - upper_limit)))
#                                 count_data[src_lower_mask] += 1#src_weight
#                             elif self.mode == 'scaled':
#                                 src_min = np.nanmin(src_data[src_lower_mask])
#                                 mask_data[src_lower_mask] = np.sqrt(
#                                     (np.power(mask_data[src_lower_mask], 2) +
#                                      np.power(src_weight * np.abs((src_data[src_lower_mask] - lower_limit) / (src_min - lower_limit)), 2))
#                                 )
#                                 count_data[src_lower_mask] += src_weight
#                             else:
#                                 mask_data[src_lower_mask] = np.sqrt(
#                                     (np.power(mask_data[src_lower_mask], 2) +
#                                      np.power(src_weight * np.abs((src_data[src_lower_mask] - upper_limit) / (lower_limit - upper_limit)), 2))
#                                 )
#                                 count_data[src_lower_mask] += src_weight

#                             # mask_data[(src_data < lower_limit)] = np.sqrt(
#                             #     (np.power(mask_data[(src_data < lower_limit)], 2) +
#                             #      np.power(src_weight * np.abs((src_lower - upper_limit) / (lower_limit - upper_limit)), 2))
#                             # )
#                             # count_data[(src_data < lower_limit)] += src_weight 
#                         except ValueError as e:
#                             pass
#                         except Exception as e:
#                             raise(e)

                        
#     def mask_gdal_dem_outliers(
#             self, srcwin_ds=None, band_data=None, mask_mask_data=None,
#             mask_count_data=None, var=None, percentile=75, upper_only=False,
#             src_weight=None, k=1.5
#     ):
#         """Generate an LSP with gdal DEMProcessing and send the result to `mask_outliers()`"""

#         tmp_ds, tmp_fn = self.gdal_dem(input_ds=srcwin_ds, var=var)
#         tmp_data = tmp_ds.GetRasterBand(1).ReadAsArray()
#         tmp_data[((np.isnan(band_data)) \
#                   | (np.isinf(band_data)) \
#                   | (tmp_data == self.ds_config['ndv']))] = np.nan
#         self.mask_outliers(
#             src_data=tmp_data, mask_data=mask_mask_data, count_data=mask_count_data,
#             percentile=percentile, upper_only=upper_only, src_weight=src_weight, k = k
#         )
#         tmp_ds = tmp_data = None
#         utils.remove_glob(tmp_fn)
#         return(0)

    
#     def get_pk(self, src_ds, var='roughness', invert=True):
#         if var == 'elevation':
#             ds_ds = src_ds
#             ds_fn = ''
#         else:
#             ds_ds, ds_fn = self.gdal_dem(input_ds=src_ds, var=var)
            
#         pk_ds, pk_fn = self.gdal_dem(input_ds=ds_ds, var='TPI')
#         pk_arr = pk_ds.GetRasterBand(1).ReadAsArray().astype(float)
#         pk_arr[(pk_arr == self.ds_config['ndv']) | (pk_arr == -9999) ] = np.nan
#         pk_arr[np.isinf(pk_arr)] = np.nan
        
#         if np.all(np.isnan(pk_arr)):
#             ds_ds = pk_ds = pk_arr = None
#             utils.remove_glob(pk_fn, ds_fn)
#             return(None, None, None)

#         med_pk = np.nanmedian(pk_arr)
#         m_pk = np.nanmean(pk_arr)
#         std_pk = np.nanstd(pk_arr)

#         # if std_pk == 0 or med_pk == 0:
#         #     pk_arr = pk_ds = ds_ds = None
#         #     utils.remove_glob(pk_fn, ds_fn)
#         #     return(75,3)

#         pkr = (np.nanmax(pk_arr) - np.nanmin(pk_arr))
#         #pkrm = pkr / std_pk #* .01
#         pkrm = (std_pk - np.nanmin(pk_arr)) / (np.nanmax(pk_arr) - np.nanmin(pk_arr))
#         kk = pkrm
#         k = (1-kk) * 4.5
#         p = kk * (.5 - 1) + 1
#         p = p * 100
#         pp = pkrm * 100
#         pk_arr = pk_ds = ds_ds = None
#         #print(med_pk, m_pk, std_pk, pkr, pkrm, kk, k, p, pp)
#         utils.remove_glob(pk_fn, ds_fn)

#         return(pp, k, p)

    
#     def rough_q(self, src_ds):
#         rp,rk,rpp = self.get_pk(src_ds, var='roughness')
#         #print('rough', rp,rk,rpp)
#         sp,sk,spp = self.get_pk(src_ds, var='slope')
#         #print('slope', sp,sk,spp)
#         tp,tk,tpp = self.get_pk(src_ds, var='TPI')
#         #print('tpi', tp,tk,tpp)
#         ep,ek,epp = self.get_pk(src_ds, var='elevation')
#         #print('elevation', ep,ek,epp)
#         cp,ck,cpp = self.get_pk(src_ds, var='curvature')
#         #print('curvature', cp,ck,cpp)

#         try:
#             k = (rk + sk + tk + ek + ck) / 5
#             p = (rp + sp + tp + ep + cp) / 5
#             pp = (rpp + spp + tpp + epp + cpp) / 5
#             #print(p, k, pp)
#             return(p, k, pp)
#         except:
#             return(np.nan, np.nan, np.nan)

        
#     def apply_mask(self, perc=75, src_ds=None):
#         """apply the generated outlier mask to the source DEM data.

#         This will calculate the outliers in the outlier mask and remove the corresponding data
#         from the input source DEM.

#         set perc to adjust the outlier calculations. k in this function is automatically 
#         generated based on the roughness of the outlier mask.
#         """

#         src_data = self.ds_band.ReadAsArray()
#         mask_mask_data = self.mask_mask_band.ReadAsArray()        
#         mask_count_data = self.mask_count_band.ReadAsArray()
#         mask_mask_data[np.isinf(mask_mask_data)] = 0
#         mask_mask_data[mask_mask_data == 0] = np.nan
#         mask_count_data[np.isinf(mask_count_data)] = 0
#         mask_count_data[mask_count_data == 0] = np.nan

#         if np.all(np.isnan(mask_mask_data)):
#             utils.echo_warning_msg('could not find any outliers')
#             return(0)
        
#         if self.mode == 'average':
#             mask_mask_data = mask_mask_data / mask_count_data
#             mask_mask_data[np.isnan(mask_mask_data)] = 0
#             self.mask_mask_band.WriteArray(mask_mask_data)
#             self.mask_mask_band.FlushCache()
#             mask_mask_data[mask_mask_data == 0] = np.nan

#         k = self.k
#         perc,self.k,perc1 = self.get_pk(self.mask_mask_ds, var='elevation')
#         if perc is None or np.isnan(perc):
#             #perc = 75
#             perc = self.percentile

#         if perc1 is None or np.isnan(perc1) or perc1 > 100 or perc1 < 0:
#             perc1 = 65

#         if self.k is None or np.isnan(self.k):
#             self.k = k
            
#         # if self.aggressive:
#         #     perc = perc1
            
#         #if self.mode == 'average':
#         #    count_upper_limit = np.nanpercentile(mask_count_data, perc)
#         #    mask_upper_limit = np.nanpercentile(mask_mask_data, perc)
#         #    #outlier_mask = (mask_mask_data > mask_upper_limit)
#         #else:
#         count_upper_limit, count_lower_limit = self.get_outliers(
#             mask_count_data, perc, k=self.k, verbose=False
#         )
#         mask_upper_limit, mask_lower_limit = self.get_outliers(
#             mask_mask_data, perc, k=self.k, verbose=False
#         )
#         outlier_mask = ((mask_mask_data > mask_upper_limit) \
#                         & (mask_count_data > count_upper_limit))

#         ## expand the outliers by a cell in each direction
#         #outlier_mask = utils.expand_for(outlier_mask, shiftx=1, shifty=1)
        
#         if self.fill_removed_data:
#             utils.echo_msg('filling filtered data from stack')
#             src_data[outlier_mask] = np.nan
#             src_data[src_data == self.ds_config['ndv']] = np.nan

#             for srcwin in utils.yield_srcwin(
#                     (self.ds_config['ny'], self.ds_config['nx']),
#                     n_chunk=self.ds_config['nx']/4,
#                     step=None, verbose=self.verbose, start_at_edge=True,
#                     msg=f'filling removed outliers with {self.interpolation}'
#             ):
#                 srcwin_src_data = src_data[srcwin[1]:srcwin[1]+srcwin[3],
#                                            srcwin[0]:srcwin[0]+srcwin[2]]
#                 srcwin_outliers = outlier_mask[srcwin[1]:srcwin[1]+srcwin[3],
#                                                srcwin[0]:srcwin[0]+srcwin[2]]

#                 #interp_data = self.generate_mem_ds(band_data=srcwin_src_data, srcwin=srcwin, return_array=True)
#                 interp_data = gdalfun.generate_mem_ds(
#                     self.ds_config,
#                     band_data=srcwin_src_data,
#                     srcwin=srcwin,
#                     return_array=True
#                 )
#                 srcwin_src_data[srcwin_outliers] = interp_data[srcwin_outliers]

#                 src_data[srcwin[1]:srcwin[1]+srcwin[3],
#                          srcwin[0]:srcwin[0]+srcwin[2]] = srcwin_src_data
                         
#             #(0,0,src_data.shape[1],src_data.shape[0]), return_array=True)
#             #interp_data = self.generate_mem_ds(band_data=src_data, srcwin=(0,0,self.ds_config['ny'],self.ds_config['nx']), return_array=True)
#             #src_data[outlier_mask] = interp_data[outlier_mask]
#             src_data[np.isnan(src_data)] = self.ds_config['ndv']
#         else:
#             src_data[outlier_mask] = self.ds_config['ndv']

#             # for b in range(1, src_ds.RasterCount+1):
#             #     this_band = src_ds.GetRasterBand(b)
#             #     this_arr = this_band.ReadAsArray()
#             #     m_band.WriteArray(this_arr)

            
#         self.ds_band.WriteArray(src_data)
#         self.ds_band.FlushCache()
        
#         # for b in range(1, src_ds.RasterCount+1):
#         #     this_band = src_ds.GetRasterBand(b)
#         #     this_arr = this_band.ReadAsArray()
#         #     this_arr[outlier_mask] = self.ds_config['ndv']
#         #     this_band.WriteArray(this_arr)
        
#         if self.verbose:
#             utils.echo_msg_bold('removed {} outliers @ <{}:{}>{}:{}.'.format(
#                 np.count_nonzero(outlier_mask), perc, self.k,
#                 mask_upper_limit, count_upper_limit
#             ))
#         self.ds_band.SetNoDataValue(self.ds_config['ndv'])
#         src_data = mask_mask_data = mask_count_data = None

#         return(np.count_nonzero(outlier_mask))

    
#     def run(self):
#         """Run the outlier module and scan a source DEM file for 
#         outliers and remove them.
#         """

#         src_ds = self.copy_src_dem()
#         if src_ds is not None:
#             self.init_ds(src_ds=src_ds)
#             self.init_chunks(src_ds=src_ds)
#             self.init_percentiles(src_ds=src_ds)
#             src_config = gdalfun.gdal_infos(src_ds)
#             # uncertainty ds
#             unc_band = None
#             if self.uncertainty_mask is not None:
#                 if self.unc_is_fn:
#                     unc_ds = gdal.Open(self.uncertainty_mask)
#                     unc_band = unc_ds.GetRasterBand(1)
#                 elif self.unc_is_band:
#                     unc_band = src_ds.GetRasterBand(self.uncertainty_mask)

#             chunks_it = np.ceil(np.linspace(self.n_chunk, self.max_chunk, self.multipass))
#             steps_it = np.ceil(np.linspace(self.n_step, self.max_step, self.multipass))
#             percs_it = np.linspace(self.percentile, self.max_percentile, self.multipass)
#             ks_it = np.linspace(self.k, self.max_k, self.multipass)
#             weights_it = np.linspace(1, 1/self.multipass, self.multipass)
#             if self.accumulate:
#                 self.generate_mask_ds(src_ds=src_ds)
                
#             for n, chunk in enumerate(chunks_it):
#                 step = steps_it[n]
#                 perc = percs_it[n]
#                 k = ks_it[n]
#                 elevation_weight = self.elevation_weight * weights_it[n]#(1/n)
#                 curvature_weight = self.curvature_weight * weights_it[n]#(1/n)
#                 slope_weight = self.slope_weight * weights_it[n]#(1/n)
#                 rough_weight = self.rough_weight * weights_it[n]#(1/n)
#                 tri_weight = self.tri_weight# * weights_it[n]#(1/n)
#                 tpi_weight = self.tpi_weight# * weights_it[n]#(1/n)
#                 n+=1
#                 if not self.accumulate:
#                     self.generate_mask_ds(src_ds=src_ds)
               
#                 for srcwin in utils.yield_srcwin(
#                         (src_ds.RasterYSize, src_ds.RasterXSize), n_chunk=chunk,
#                         step=step, verbose=self.verbose, start_at_edge=True,
#                         msg='scanning for outliers ({}:{})'.format(perc, k),
#                 ):
#                     band_data = self.ds_band.ReadAsArray(*srcwin)
#                     band_data[band_data == self.ds_config['ndv']] = np.nan
#                     if np.all(np.isnan(band_data)):
#                         band_data = None
#                         continue

#                     ## read in the mask data for the srcwin
#                     # read in the mask id data
#                     mask_mask_data = self.mask_mask_band.ReadAsArray(*srcwin)
#                     # read in the count data
#                     mask_count_data = self.mask_count_band.ReadAsArray(*srcwin) 
                        
#                     ## generate a mem datasource to feed into gdal.DEMProcessing
#                     #srcwin_ds = self.generate_mem_ds(band_data=band_data, srcwin=srcwin) # possibly interpolated
#                     srcwin_ds = gdalfun.generate_mem_ds(
#                         self.ds_config,
#                         band_data=band_data,
#                         srcwin=srcwin,
#                         return_array=False
#                     )
#                     if srcwin_ds is None:
#                         band_data = None
#                         continue
                    
#                     slp_ds, slp_fn = self.gdal_dem(input_ds=srcwin_ds, var='slope')
#                     #curv_ds, curv_fn = self.gdal_dem(input_ds=slp_ds, var='slope')
#                     rough_ds, rough_fn = self.gdal_dem(input_ds=srcwin_ds, var='roughness')
#                     #p, k = self.rough_q(srcwin_ds)
#                     # if k is None:
#                     #     srcwin_ds = slp_ds = rough_ds = None
#                     #     utils.remove_glob(slp_fn, rough_fn)
#                     #     continue

#                     ## apply elevation outliers
#                     #perc,k,p = self.get_pk(srcwin_ds, var='elevation')
#                     self.mask_outliers(
#                         src_data=band_data,
#                         mask_data=mask_mask_data,
#                         count_data=mask_count_data,
#                         percentile=perc, 
#                         src_weight=elevation_weight,
#                         k=k
#                     )                    
#                     ## apply slope outliers
#                     #perc,k,p = self.get_pk(srcwin_ds, var='slope')                    
#                     self.mask_gdal_dem_outliers(
#                         srcwin_ds=srcwin_ds,
#                         band_data=band_data,
#                         mask_mask_data=mask_mask_data,
#                         mask_count_data=mask_count_data,
#                         percentile=perc, 
#                         upper_only=False,
#                         src_weight=slope_weight,
#                         var='slope',
#                         k=k
#                     )
#                     ## apply tri outliers
#                     #perc,k,p = self.get_pk(srcwin_ds, var='tri')                    
#                     self.mask_gdal_dem_outliers(
#                         srcwin_ds=srcwin_ds,
#                         band_data=band_data,
#                         mask_mask_data=mask_mask_data,
#                         mask_count_data=mask_count_data,
#                         percentile=perc,
#                         upper_only=False,
#                         src_weight=tri_weight,
#                         var='TRI',
#                         k=k
#                     )
#                     ## apply curvature outliers ## tmp
#                     #perc,k,p = self.get_pk(slp_ds, var='slope')
#                     # self.mask_gdal_dem_outliers(srcwin_ds=srcwin_ds, band_data=band_data, mask_mask_data=mask_mask_data,
#                     #                             mask_count_data=mask_count_data, percentile=perc,
#                     #                             upper_only=True, src_weight=curvature_weight, var='curvature', k=k)
#                     ## apply roughness outliers
#                     #perc,k,p = self.get_pk(srcwin_ds, var='roughness')                    
#                     self.mask_gdal_dem_outliers(
#                         srcwin_ds=srcwin_ds,
#                         band_data=band_data,
#                         mask_mask_data=mask_mask_data,
#                         mask_count_data=mask_count_data,
#                         percentile=perc,
#                         upper_only=False,
#                         src_weight=rough_weight,
#                         var='roughness',
#                         k=k
#                     )
#                     ## apply TPI outliers ## tmp
#                     # self.mask_gdal_dem_outliers(srcwin_ds=curv_ds, band_data=band_data, mask_mask_data=mask_mask_data,
#                     #                             mask_count_data=mask_count_data, percentile=perc,
#                     #                             upper_only=False, src_weight=curvature_weight, var='TPI', k=k)
#                     # self.mask_gdal_dem_outliers(srcwin_ds=rough_ds, band_data=band_data, mask_mask_data=mask_mask_data,
#                     #                             mask_count_data=mask_count_data, percentile=perc,
#                     #                             upper_only=False, src_weight=rough_weight, var='TPI', k=k)
#                     #perc,k,p = self.get_pk(slp_ds, var='elevation')                    
#                     self.mask_gdal_dem_outliers(
#                         srcwin_ds=slp_ds,
#                         band_data=band_data,
#                         mask_mask_data=mask_mask_data,
#                         mask_count_data=mask_count_data,
#                         percentile=perc,
#                         upper_only=False,
#                         src_weight=slope_weight,
#                         var='TPI',
#                         k=k
#                     )
#                     #perc,k,p = self.get_pk(srcwin_ds, var='TPI')                    
#                     self.mask_gdal_dem_outliers(
#                         srcwin_ds=srcwin_ds,
#                         band_data=band_data,
#                         mask_mask_data=mask_mask_data,
#                         mask_count_data=mask_count_data,
#                         percentile=perc,
#                         upper_only=False,
#                         src_weight=tpi_weight,
#                         var='TPI',
#                         k=k
#                     )

#                     srcwin_ds = slp_ds = rough_ds = None
#                     utils.remove_glob(slp_fn, rough_fn)

#                     ## write the mask data to file
#                     self.mask_mask_band.WriteArray(
#                         mask_mask_data, srcwin[0], srcwin[1]
#                     )
#                     self.mask_count_band.WriteArray(
#                         mask_count_data, srcwin[0], srcwin[1]
#                     )
#                     band_data = mask_mask_data = mask_count_data = None

#                 if not self.accumulate:
#                     outliers = self.apply_mask(self.percentile)
#                     #if outliers == 0:
#                     #    break
                    
#                     self.mask_mask_ds = None
                    
#             if self.accumulate:
#                 outliers = self.apply_mask(self.percentile, src_ds)
#                 self.mask_mask_ds = None
                
#             unc_ds = src_ds = None
#             if not self.return_mask and not self.accumulate:
#                 utils.remove_glob(self.mask_mask_fn)

#             return(self.dst_dem, 0)
#         else:
#             return(None, -1)

# ### End
