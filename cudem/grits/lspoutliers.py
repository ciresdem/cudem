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
