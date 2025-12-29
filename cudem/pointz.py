### pointz.py 
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
##
## pointz.py is part of CUDEM
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
## Point Cloud Filtering and Manipulation Tools.
##
### Code:

import os
import sys
import argparse
import numpy as np
from osgeo import ogr, gdal

from cudem import utils
from cudem import gdalfun
from cudem import regions
from cudem import factory
from cudem import xyzfun
from cudem.fetches import fetches
from cudem import __version__ as __cudem_version__

__version__ = '0.1.0'

# =============================================================================
# Point Gridding Helper
# =============================================================================
class PointPixels:
    """Bins point cloud data into a grid coinciding with a desired region.
    Returns aggregated values (Z, Weights, Uncertainty) for each grid cell.
    
    Incoming data are numpy structured arrays (rec-arrays) of x, y, z, <w, u>.
    """

    def __init__(self, src_region=None, x_size=None, y_size=None, 
                 verbose=True, ppm=False, **kwargs):
        self.src_region = src_region
        self.x_size = utils.int_or(x_size, 10)
        self.y_size = utils.int_or(y_size, 10)
        self.verbose = verbose
        self.ppm = ppm  # Pixel-Per-Metric flag? (Unused logic in original kept for legacy)
        self.dst_gt = None

        
    def init_region_from_points(self, points):
        """Initialize the source region based on point extents."""
        
        if self.src_region is None:
            self.src_region = regions.Region().from_list([
                np.min(points['x']), np.max(points['x']),
                np.min(points['y']), np.max(points['y'])
            ])
        ## If the points create an invalid region, lets buffer it
        ## until it becomes valid...
        if not self.src_region.is_valid():
            self.src_region.buffer(pct=2, check_if_valid=False)
        self.init_gt()
        
    def init_gt(self):
        """Initialize the GeoTransform based on region and size."""
        
        if self.src_region is not None:
            self.dst_gt = self.src_region.geo_transform_from_count(
                x_count=self.x_size, y_count=self.y_size
            )

            
    def __call__(self, points, weight=1, uncertainty=0, mode='mean'):
        """Process points into a gridded array.
        
        Args:
            points (np.recarray): Input data containing 'x', 'y', 'z'.
            weight (float): Global weight multiplier.
            uncertainty (float): Global uncertainty value.
            mode (str): Aggregation mode ('mean', 'min', 'max', 'sums').
        """
        
        out_arrays = {
            'z': None, 'count': None, 'weight': None, 'uncertainty': None,
            'mask': None, 'x': None, 'y': None, 'pixel_x': None, 'pixel_y': None
        }

        if points is None or len(points) == 0:
            return out_arrays, None, None

        ## Ensure region and geotransform are set
        if self.src_region is None:
            self.init_region_from_points(points)
        elif self.dst_gt is None:
            self.init_gt()

        ## Sanitize inputs
        weight = utils.float_or(weight, 1)
        uncertainty = utils.float_or(uncertainty, 0)

        ## Extract data columns
        points_x = np.array(points['x'])
        points_y = np.array(points['y'])
        pixel_z = np.array(points['z'])

        ## Handle optional fields
        pixel_w = np.array(points['w']) if 'w' in points.dtype.names else np.ones_like(pixel_z)
        pixel_u = np.array(points['u']) if 'u' in points.dtype.names else np.zeros_like(pixel_z)

        pixel_w[np.isnan(pixel_w)] = 1
        pixel_u[np.isnan(pixel_u)] = 0

        ## Convert to pixel coordinates
        ## dst_gt: [origin_x, pixel_width, 0, origin_y, 0, pixel_height]
        pixel_x = np.floor((points_x - self.dst_gt[0]) / self.dst_gt[1]).astype(int)
        pixel_y = np.floor((points_y - self.dst_gt[3]) / self.dst_gt[5]).astype(int)

        ## Filter pixels outside window
        valid_mask = (
            (pixel_x >= 0) & (pixel_x < self.x_size) &
            (pixel_y >= 0) & (pixel_y < self.y_size)
        )
        
        if not np.any(valid_mask):
            return out_arrays, None, None

        # Apply mask
        pixel_x = pixel_x[valid_mask]
        pixel_y = pixel_y[valid_mask]
        pixel_z = pixel_z[valid_mask]
        pixel_w = pixel_w[valid_mask]
        pixel_u = pixel_u[valid_mask]
        points_x = points_x[valid_mask]
        points_y = points_y[valid_mask]

        if len(pixel_x) == 0:
            return out_arrays, None, None

        ## Local Source Window Calculation
        min_px, max_px = int(np.min(pixel_x)), int(np.max(pixel_x))
        min_py, max_py = int(np.min(pixel_y)), int(np.max(pixel_y))
        
        this_srcwin = (min_px, min_py, max_px - min_px + 1, max_py - min_py + 1)

        ## Shift to local coordinates
        local_px = pixel_x - min_px
        local_py = pixel_y - min_py
        
        ## Unique pixel identification (row-major: y, x)
        pixel_xy = np.vstack((local_py, local_px)).T

        unq, unq_idx, unq_inv, unq_cnt = np.unique(
            pixel_xy, axis=0, return_inverse=True,
            return_index=True, return_counts=True
        )

        ## Initial population with first occurrences
        if mode == 'sums':
            ww = pixel_w[unq_idx] * weight
            zz = pixel_z[unq_idx] * ww
            xx = points_x[unq_idx] * ww
            yy = points_y[unq_idx] * ww
        else:
            zz = pixel_z[unq_idx]
            ww = pixel_w[unq_idx]
            xx = points_x[unq_idx]
            yy = points_y[unq_idx]
            
        uu = pixel_u[unq_idx]
        
        ## --- Handle Duplicates (Aggregation) ---
        cnt_msk = unq_cnt > 1
        if np.any(cnt_msk):
            # Sort indices to group by pixel
            srt_idx = np.argsort(unq_inv)
            split_indices = np.cumsum(unq_cnt)[:-1]
            grouped_indices = np.split(srt_idx, split_indices)
            
            # Filter groups with duplicates
            dup_indices = [grouped_indices[i] for i in np.flatnonzero(cnt_msk)]
            dup_stds = []

            if mode == 'min':
                zz[cnt_msk] = [np.min(pixel_z[idx]) for idx in dup_indices]
                xx[cnt_msk] = [np.min(pixel_x[idx]) for idx in dup_indices] # Use min X/Y? Or corresponding X/Y? Using min for consistency with logic
                yy[cnt_msk] = [np.min(pixel_y[idx]) for idx in dup_indices]
                dup_stds = np.zeros(len(dup_indices))

            elif mode == 'max':
                zz[cnt_msk] = [np.max(pixel_z[idx]) for idx in dup_indices]
                xx[cnt_msk] = [np.max(pixel_x[idx]) for idx in dup_indices]
                yy[cnt_msk] = [np.max(pixel_y[idx]) for idx in dup_indices]
                dup_stds = np.zeros(len(dup_indices))

            elif mode == 'mean':
                zz[cnt_msk] = [np.mean(pixel_z[idx]) for idx in dup_indices]
                xx[cnt_msk] = [np.mean(points_x[idx]) for idx in dup_indices]
                yy[cnt_msk] = [np.mean(points_y[idx]) for idx in dup_indices]
                dup_stds = [np.std(pixel_z[idx]) for idx in dup_indices]

            elif mode == 'sums':
                zz[cnt_msk] = [np.sum(pixel_z[idx] * pixel_w[idx] * weight) for idx in dup_indices]
                xx[cnt_msk] = [np.sum(points_x[idx] * pixel_w[idx] * weight) for idx in dup_indices]
                yy[cnt_msk] = [np.sum(points_y[idx] * pixel_w[idx] * weight) for idx in dup_indices]
                ww[cnt_msk] = [np.sum(pixel_w[idx] * weight) for idx in dup_indices]
                dup_stds = [np.std(pixel_z[idx]) for idx in dup_indices]

            # Aggregate uncertainty
            uu[cnt_msk] = np.sqrt(np.power(uu[cnt_msk], 2) + np.power(dup_stds, 2))

        ## --- Fill Output Grids ---
        grid_shape = (this_srcwin[3], this_srcwin[2]) # rows, cols

        def fill_grid(values, fill_val=np.nan):
            grid = np.full(grid_shape, fill_val)
            grid[unq[:, 0], unq[:, 1]] = values
            return grid

        out_arrays['z'] = fill_grid(zz)
        out_arrays['x'] = fill_grid(xx)
        out_arrays['y'] = fill_grid(yy)
        out_arrays['count'] = fill_grid(unq_cnt, fill_val=0)
        
        # Uncertainty
        out_arrays['uncertainty'] = fill_grid(
            np.sqrt(uu**2 + (uncertainty)**2), fill_val=0
        )

        # Weights
        out_arrays['weight'] = np.ones(grid_shape)
        if mode == 'sums':
            out_arrays['weight'][unq[:, 0], unq[:, 1]] = ww
        else:
            out_arrays['weight'][:] = weight
            out_arrays['weight'][unq[:, 0], unq[:, 1]] *= (ww * unq_cnt)

        # Helper coords for calling class to map back
        out_arrays['pixel_x'] = local_px
        out_arrays['pixel_y'] = local_py

        return out_arrays, this_srcwin, self.dst_gt

    
## =============================================================================
## PointZ Base Class
## =============================================================================
class PointZ:
    """Base Class for Point Data Filters.

    Wraps an array of xyz data and provides methods to manipulate 
    and filter that array.
    """

    def __init__(self, points=None, region=None, verbose=False,
                 xyinc=None, cache_dir='.', **kwargs):
        
        self.points = points
        self.region = region
        
        # Auto-calculate region if not provided and points exist
        if self.points is not None and len(self.points) > 0 and self.region is None:
            self.region = self.init_region()

        if xyinc is not None:
            self.xyinc = [utils.str2inc(x) for x in xyinc]
        else:
            self.xyinc = None
            
        self.verbose = verbose
        self.cache_dir = cache_dir
        self.kwargs = kwargs

        
    def __call__(self):
        """Execute the filter."""
        if self.verbose:
            utils.echo_msg(f'Filtering points using {self.__class__.__name__}...')

        if self.points is None or len(self.points) == 0:
            return self.points

        outliers_mask = self.run()
        
        ## If run() returns a boolean mask (True=Outlier), apply it. 
        if isinstance(outliers_mask, np.ndarray) and outliers_mask.dtype == bool:
            if self.verbose:
                utils.echo_msg(f'Filtered {np.count_nonzero(outliers_mask)} records.')
            return self.points[~outliers_mask]
        
        ## Backward compatibility if run returns points directly
        return outliers_mask

    
    def run(self):
        """Override in subclasses to perform filtering."""
        
        return np.zeros(len(self.points), dtype=bool)

    
    def init_region(self):
        """Initialize the data-region AOI from point extents."""
        
        if self.points is None:
            return None
        return regions.Region().from_list([
            np.min(self.points['x']), np.max(self.points['x']),
            np.min(self.points['y']), np.max(self.points['y'])
        ])

    
    def fetch_data(self, fetches_module, check_size=True):
        """Fetch data from a fetches module for the data-region."""
        
        this_fetches = fetches.FetchesFactory(
            mod=fetches_module,
            src_region=self.region,
            verbose=self.verbose,
            callback=fetches.fetches_callback
        )._acquire_module()
        
        if this_fetches is None:
            return None

        this_fetches.run()
        fr = fetches.fetch_results(this_fetches, check_size=check_size)
        fr.daemon = True
        fr.start()
        fr.join()
        return fr

    
    def point_pixels(self, points, x_size=50, y_size=50):
        """Bin the points to a grid of x_size/y_size and return the 
        associated pixel-z-data at the x/y locations of the points.
        Used for residual analysis.
        """
        
        try:
            pa = PointPixels(x_size=x_size, y_size=y_size, ppm=True)
            point_arrays, _, _ = pa(points, mode='mean')

            ## Map pixels back to point indices
            ## The PointPixels __call__ returns arrays aligned to the grid.
            ## We need to map point locations to grid values.
            ## Re-calculating pixel indices here for mapping since PA returns aggregated grids.
            
            ## Use the dst_gt calculated by PA
            dst_gt = pa.dst_gt
            if dst_gt is None: return None

            pixel_x = np.floor((points['x'] - dst_gt[0]) / dst_gt[1]).astype(int)
            pixel_y = np.floor((points['y'] - dst_gt[3]) / dst_gt[5]).astype(int)
            
            ## Clip indices to grid size to prevent index errors
            grid_h, grid_w = point_arrays['z'].shape
            
            ## Mask valid indices
            valid_mask = (pixel_x >= 0) & (pixel_x < grid_w) & (pixel_y >= 0) & (pixel_y < grid_h)
            
            result = np.full(len(points), np.nan)
            
            ## Extract Z from grid at point locations
            result[valid_mask] = point_arrays['z'][pixel_y[valid_mask], pixel_x[valid_mask]]
            
            return result

        except Exception as e:
            if self.verbose:
                utils.echo_error_msg(f"Error calculating point pixels: {e}")
            return None

        
    def vectorize_points(self):
        """Convert internal points to an OGR Memory Layer."""
        
        dst_ogr = 'points_dataset'
        driver = gdal.GetDriverByName('Memory')
        ogr_ds = driver.Create('', 0, 0, 0, gdal.GDT_Unknown)
        layer = ogr_ds.CreateLayer(dst_ogr, geom_type=ogr.wkbPoint)
        
        fd = ogr.FieldDefn('index', ogr.OFTInteger)
        layer.CreateField(fd)
        f_defn = layer.GetLayerDefn()
        
        with utils.ccp(total=len(self.points), desc='Vectorizing points', leave=False) as pbar:
            for index, this_row in enumerate(self.points):
                feat = ogr.Feature(f_defn)
                feat.SetField(0, index)
                geom = ogr.CreateGeometryFromWkt(f'POINT ({this_row["x"]} {this_row["y"]})')
                feat.SetGeometryDirectly(geom)
                layer.CreateFeature(feat)
                pbar.update()
        return ogr_ds

    
    def mask_to_raster(self, mask_fn):
        """Rasterize a vector mask for GDAL querying."""
        
        return gdalfun.ogr2gdal_mask(
            mask_fn,
            region=self.region,
            x_inc=self.xyinc[0] if self.xyinc else 0.0000925,
            y_inc=self.xyinc[1] if self.xyinc else 0.0000925,
            verbose=self.verbose,
            temp_dir=self.cache_dir
        )

    
## =============================================================================
## Vector Mask Filter
## =============================================================================
class PointZVectorMask(PointZ):
    """Filter data using a vector mask (Shapefile, GeoJSON, etc).

    Config: <vector_mask:mask_fn=path:invert=False>
    """
    
    def __init__(self, mask_fn=None, invert=False, **kwargs):
        super().__init__(**kwargs)
        self.mask_fn = mask_fn
        self.invert = invert
        if self.verbose:
            utils.echo_msg(f'Masking with {mask_fn}')

            
    def mask_points(self, points, invert=False):
        """Mask points by rasterizing the input vector mask."""
        
        if self.verbose:
            utils.echo_msg(f'Using mask dataset: {self.mask_fn}')
            
        if self.mask_fn is None or not os.path.exists(self.mask_fn):
            utils.echo_error_msg(f"Mask file not found: {self.mask_fn}")
            return np.zeros(len(points), dtype=bool)

        ogr_or_gdal = gdalfun.ogr_or_gdal(self.mask_fn)
        mask_raster = self.mask_to_raster(self.mask_fn) if ogr_or_gdal == 1 else self.mask_fn
            
        ## Sample the mask at point locations
        sampled_values = gdalfun.gdal_query(points, mask_raster, 'g').flatten()

        ## mask raster: 1=Inside, 0=Outside
        is_outside = sampled_values == 0
        return is_outside

    
    def run(self):
        outliers = self.mask_points(self.points, invert=self.invert)
        ## If invert is True, we want to remove Inside points.
        ## mask_points returns True for Outside.
        ## If invert, return ~outliers (True for Inside)
        return ~outliers if self.invert else outliers

    
## =============================================================================
## Outlier Filters
## =============================================================================
class PointZOutlier(PointZ):
    """XYZ outlier filter based on local block statistics.

    Config: <outlierz:percentile=98:multipass=4:invert=False:res=50>
    """
    
    def __init__(self, percentile=98, max_percentile=99.9,
                 multipass=4, percentage=False, invert=False,
                 res=50, max_res=5000, **kwargs):
        super().__init__(**kwargs)
        self.percentile = utils.float_or(percentile, 98)
        self.max_percentile = utils.float_or(max_percentile, 99.9)
        self.multipass = utils.int_or(multipass, 1)
        self.percentage = percentage
        self.invert = invert
        self.res = utils.int_or(res, 50)
        self.max_res = utils.int_or(max_res, 5000)

        
    def point_residuals(self, points, percentage=False, res=50):
        """Calculate residuals relative to local block mean."""
        
        point_pixels = self.point_pixels(points, x_size=res, y_size=res)
        if point_pixels is None:
            return None

        diff = points['z'] - point_pixels
        
        if percentage:
            with np.errstate(divide='ignore', invalid='ignore'):
                residuals = np.abs(diff / point_pixels) * 100
                residuals[~np.isfinite(residuals)] = 0
        else:
            residuals = np.abs(diff)
        return residuals

    
    def find_outliers(self, residuals, percentile=98, percentile_is_threshold=False):
        """Identify outliers based on percentile threshold."""
        
        if percentile_is_threshold:
            outlier_threshold = percentile
        else:
            outlier_threshold = np.nanpercentile(residuals, percentile)

        outliers = residuals > outlier_threshold
        if self.verbose:
            utils.echo_msg(f'Found {np.count_nonzero(outliers)} outliers > {outlier_threshold:.8f}')
        return outliers

    
    def filter_points(self, points, percentile=92, res=50, percentage=False):
        residuals = self.point_residuals(points, percentage=percentage, res=res)
        if residuals is not None:
            return self.find_outliers(residuals, percentile=percentile, percentile_is_threshold=percentage)
        return np.zeros(len(points), dtype=bool)
    

    def run(self):
        """Execute multipass filtering with master mask tracking."""
        
        percs_it = np.linspace(self.percentile, self.max_percentile, self.multipass)
        res_it = np.linspace(self.max_res, self.res, self.multipass)
        
        ## mask: True = Remove
        total_outliers_mask = np.zeros(len(self.points), dtype=bool)
        
        for mpass in range(self.multipass):
            if self.verbose:
                utils.echo_msg(f"Pass {mpass+1}/{self.multipass}: Res={res_it[mpass]:.1f}, Thresh={percs_it[mpass]:.2f}")
            
            valid_indices = np.nonzero(~total_outliers_mask)[0]
            if len(valid_indices) == 0: break

            current_subset = self.points[valid_indices]
            
            subset_outliers = self.filter_points(
                current_subset, 
                percentile=percs_it[mpass],
                res=res_it[mpass], 
                percentage=self.percentage
            )
            
            if np.count_nonzero(subset_outliers) == 0:
                if self.verbose: utils.echo_msg("No new outliers found.")
                break
                
            indices_to_remove = valid_indices[subset_outliers]
            total_outliers_mask[indices_to_remove] = True
            
        return ~total_outliers_mask if self.invert else total_outliers_mask


class RQOutlierZ(PointZOutlier):
    """XYZ outlier filter using a Reference Raster (RQ).

    Config: <rq:threshold=5:raster=None>
    """
    
    def __init__(self, threshold=10, raster=None, scaled_percentile=False,
                 resample_raster=True, **kwargs):
        ## Remove args intended for parent class manual override
        kwargs.pop('percentile', None)
        kwargs.pop('percentage', None)
        kwargs.pop('multipass', None)
            
        super().__init__(
            percentile=threshold, percentage=True, multipass=1,
            **kwargs
        )
        self.threshold = threshold
        self.resample_raster = resample_raster
        self.fetches_modules = ['gmrt', 'CUDEM', 'etopo:datatype=surface']
        self.raster = self.init_raster(raster)
        self.scaled_percentile = scaled_percentile

        
    def mask_gmrt(self, raster_path):
        out_fn = f'{raster_path}_swath.tif'
        if os.path.exists(out_fn): return out_fn

        from cudem.fetches.gmrt import GMRT_SWATH_URL
        fetcher = fetches.FetchesFactory(mod='gmrt', src_region=self.region)._acquire_module()
        swath_zip = os.path.join(fetcher._outdir, 'gmrt_swath_polygons.zip')
        
        if fetches.Fetch(GMRT_SWATH_URL, verbose=self.verbose).fetch_file(swath_zip) == 0:
            extracted = utils.p_unzip(swath_zip, exts=['shp'], outdir=fetcher._outdir)
            swath_shp = next((f for f in extracted if f.endswith('.shp')), None)
            if swath_shp:
                gdalfun.gdal_clip(
                    raster_path, out_fn, src_ply=swath_shp, invert=True,
                    verbose=self.verbose
                )
                return out_fn
        return None

    
    def init_raster(self, raster):
        """Initialize the reference raster list."""
        
        if raster and os.path.exists(raster) and os.path.isfile(raster):
            return [raster]
        
        ## Fetch Defaults
        ## GMRT/ETOPO/CUDEM
        rasters = []
        gmrt_fetch = self.fetch_data('gmrt')
        if gmrt_fetch:
            for r in gmrt_fetch.results:
                fn = r[1]
                if fn.endswith('.grd'): fn = gdalfun.gmt_grd2gdal(fn, verbose=False)
                #masked_fn = self.mask_gmrt(fn)
                masked_fn = None
                rasters.append(masked_fn if masked_fn else fn)

        if not rasters:
            cudem_fetch = self.fetch_data('CUDEM:datatype=13')
            if cudem_fetch:
                rasters.extend([r[1] for r in cudem_fetch.results])

        return rasters

    
    def point_residuals(self, points, percentage=True, res=50):
        """Calculate residuals against the loaded raster(s)."""
        
        if not self.raster:
            return None

        ## Handle Multiple Rasters via Warping/Merging
        ref_raster = self.raster[0]
        temp_merged = None
        
        if len(self.raster) > 1 or self.resample_raster:
            tmp_raster_fn = os.path.join(self.cache_dir, f'rq_merged_{utils.this_year()}.vrt')
            if self.region is None: self.region = self.init_region()
            
            x_res = self.xyinc[0] if self.xyinc else 0.0000925
            y_res = self.xyinc[1] if self.xyinc else 0.0000925

            try:
                warped = gdalfun.sample_warp(
                    self.raster, tmp_raster_fn, x_res, y_res, 
                    sample_alg='bilinear', src_region=self.region, 
                    verbose=self.verbose
                )
                if warped:
                    ref_raster = warped[0]
                    temp_merged = ref_raster
            except Exception as e:
                utils.echo_error_msg(f"Failed to merge rasters: {e}")
                return None

        ## Interpolate
        sampled_z = gdalfun.gdal_query(points, ref_raster, 'g').flatten()
        
        if temp_merged and os.path.exists(temp_merged):
            if temp_merged.endswith('.vrt'):
                try: os.remove(temp_merged)
                except OSError: pass

        diff = points['z'] - sampled_z
        
        if percentage:
            with np.errstate(divide='ignore', invalid='ignore'):
                if self.scaled_percentile:
                    vals = np.abs(diff / (points['z'] + sampled_z)) * 100
                else:
                    vals = np.abs(diff / sampled_z) * 100
                vals[~np.isfinite(vals)] = 0
            return vals
        else:
            return np.abs(diff)

## =============================================================================
## Factory & CLI
## =============================================================================
class PointFilterFactory(factory.CUDEMFactory):
    _modules = {
        'outlierz': {'name': 'outlierz', 'call': PointZOutlier},
        'rq': {'name': 'rq', 'call': RQOutlierZ},
        'vector_mask': {'name': 'vector_mask', 'call': PointZVectorMask},
    }

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        
class PrintModulesAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        factory.echo_modules(PointFilterFactory._modules, values)
        sys.exit(0)

        
def pointz_cli():
    parser = argparse.ArgumentParser(
        description=f"%(prog)s ({__version__}): Filter Point Cloud Data (XYZ)",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=f"""
Supported %(prog)s modules (see %(prog)s --modules <module-name> for more info): 
{factory.get_module_short_desc(PointFilterFactory._modules)}
        
CUDEM home page: <http://cudem.colorado.edu>
        """

    )
    parser.add_argument("input", help="Input XYZ file")
    parser.add_argument("output", help="Output filtered XYZ file")
    parser.add_argument("-M", "--module", action="append", help="Filter module string (e.g. 'outlierz:percentile=95')")
    parser.add_argument("-q", "--quiet", action="store_true", help="Suppress output")
    parser.add_argument(
        '-m', '--modules',
        nargs='?',
        default=None,
        action=PrintModulesAction,
        help="Display the module descriptions and usage"
    )
    parser.add_argument(
        '--version',
        action='version',
        version=f'CUDEM {__cudem_version__} :: %(prog)s {__version__}'
    )

    
    args = parser.parse_args()
    verbose = not args.quiet

    from cudem.datalists.dlim import DatasetFactory
    from numpy.lib.recfunctions import stack_arrays
    
    if verbose: utils.echo_msg(f"Loading {args.input}...")
    try:
        xyz_file = DatasetFactory(
            mod=args.input,
            verbose=True,
        )._acquire_module().initialize()
        ## Load all points into memory
        points = None
        for p in xyz_file.yield_points():
            if points is None:
                points = p
            else:
                points = stack_arrays((points, p), asrecarray=True, usemask=False)

        points = points[np.isfinite(points['z'])]
        #points = xyz_file.export_xyz_as_pandas()
        #points = points.to_records(index=False)
    except Exception as e:
        utils.echo_error_msg(f"Failed to load input: {e}")
        sys.exit(1)

    if args.module:
        for f_str in args.module:
            if verbose: utils.echo_msg(f"Applying filter: {f_str}")
            f_parts = f_str.split(':')
            f_name = f_parts[0]
            f_opts = factory.fmod2dict(':'.join(f_parts[1:]), {}) if len(f_parts) > 1 else {}
            
            try:
                f_opts['verbose'] = verbose
                f_opts['points'] = points
                filt_obj = PointFilterFactory(mod=f_name, **f_opts)._acquire_module()
                if filt_obj:
                    points = filt_obj()
                else:
                    utils.echo_warning_msg(f"Filter '{f_name}' not found.")
            except Exception as e:
                utils.echo_error_msg(f"Filter execution failed: {e}")

    if verbose: utils.echo_msg(f"Saving to {args.output}...")
    try:
        with open(args.output, 'w') as f:
            for p in points:
                f.write(f"{p['x']} {p['y']} {p['z']}\n")
    except Exception as e:
        utils.echo_error_msg(f"Failed to write output: {e}")
        sys.exit(1)

        
if __name__ == "__main__":
    pointz_cli()

    
### End
