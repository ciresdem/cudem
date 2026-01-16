### pointz.py 
##
## Copyright (c) 2010 - 2026 Regents of the University of Colorado
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
## Pointz: The CUDEM Point Cloud Filtering Engine.
##
## Pointz provides specialized tools for filtering, thinning, and manipulating
## raw point cloud data (XYZ, LAS/LAZ) prior to gridding. It ensures data quality
## by removing noise and reducing redundancy.
##
##   * Outlier Removal:
##      - Statistical filters (Local Z-Score/IQR) to remove spikes.
##      - 'Coplanar' filters to detect points deviating from local planes.
##      - 'RQ' (Reference Quotient) filters to validate data against a known basemap.
##
##   * Data Thinning:
##      - Grid-based thinning (block_thin) to normalize point density.
##      - Shoal-biased thinning (block_minmax) for hydrographic safety.
##      - Random, Median, or Center-based decimation.
##
##   * Masking & Selection:
##      - Filter points using vector polygons or raster masks.
##      - 'RangeZ' and 'DiffZ' filters to select points based on absolute elevation
##        or difference from a reference surface.
##
## Usage:
##   CLI: pointz input.xyz output.xyz -M outlierz:percentile=95 -M block_thin:res=10
##   API: p = PointFilterFactory(mod='outlierz', points=my_points).acquire()
##        clean_points = p.run()
##
### Code:

import os
import sys
import argparse
import warnings
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
        self.ppm = ppm 
        self.dst_gt = None

        
    # def init_region_from_points(self, points):
    #     """Initialize the source region based on point extents."""
        
    #     if self.src_region is None:
    #         self.src_region = regions.Region().from_list([
    #             np.min(points['x']), np.max(points['x']),
    #             np.min(points['y']), np.max(points['y'])
    #         ])
    #     if not self.src_region.is_valid():
    #         self.src_region.buffer(pct=2, check_if_valid=False)
    #     self.init_gt()

    def init_region_from_points(self, points):
        """Initialize the source region based on point extents."""
        
        if self.src_region is None:
            self.src_region = regions.Region().from_list([
                np.min(points['x']), np.max(points['x']),
                np.min(points['y']), np.max(points['y'])
            ])
            
        ## If percentage buffering results in 0 (because width is 0),
        ## fallback to a small fixed epsilon.
        if not self.src_region.is_valid():
            self.src_region.buffer(pct=2, check_if_valid=False)
            
            # If still invalid (width/height is 0), force a tiny buffer
            if not self.src_region.is_valid():
                epsilon = 0.00001
                self.src_region.buffer(x_bv=epsilon, y_bv=epsilon, check_if_valid=False)
                
        self.init_gt()

        
    def init_gt(self):
        """Initialize the GeoTransform based on region and size."""
        
        if self.src_region is not None:
            self.dst_gt = self.src_region.geo_transform_from_count(
                x_count=self.x_size, y_count=self.y_size
            )

            
    def __call__(self, points, weight=1.0, uncertainty=0.0, mode='mean'):
        """Process points into a gridded array.
        
        Args:
            points (np.recarray): Input data containing 'x', 'y', 'z'.
            weight (float): Global weight multiplier.
            uncertainty (float): Global uncertainty value.
            mode (str): Aggregation mode.
                        Options: 'mean', 'min', 'max', 'median', 'std', 'var', 'sums'.
        """

        ## mrl: removed 'mask': None
        out_arrays = {
            'z': None, 'count': None, 'weight': None, 'uncertainty': None,
            'x': None, 'y': None, 'pixel_x': None, 'pixel_y': None
        }

        if points is None or len(points) == 0:
            return out_arrays, None, None

        ## If input points are pandas dataframe, tranform it to recarray
        if hasattr(points, 'to_records'):
            points = points.to_records(index=False)
            
        ## Ensure region and geotransform are set
        if self.src_region is None:
            self.init_region_from_points(points)
        elif self.dst_gt is None:
            self.init_gt()

        ## Sanitize inputs
        weight = utils.float_or(weight, 1)
        uncertainty = utils.float_or(uncertainty, 0.0)
        mode = mode.lower()

        ## Extract data columns
        points_x = np.array(points['x'])
        points_y = np.array(points['y'])
        pixel_z = np.array(points['z'])

        ## Handle optional fields
        ## This still gives a warning sometimes:
        ##   RuntimeWarning: invalid value encountered in divide
        ##   pixel_x = np.floor((points_x - self.dst_gt[0]) / self.dst_gt[1]).astype(int)
        ##   RuntimeWarning: invalid value encountered in cast
        ## TODO: Figure this out and fix.
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

        ## Apply mask
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
            ## Sort indices to group by pixel
            srt_idx = np.argsort(unq_inv)
            split_indices = np.cumsum(unq_cnt)[:-1]
            grouped_indices = np.split(srt_idx, split_indices)
            
            ## Filter groups with duplicates
            dup_indices = [grouped_indices[i] for i in np.flatnonzero(cnt_msk)]
            dup_stds = []

            if mode == 'min':
                zz[cnt_msk] = [np.min(pixel_z[idx]) for idx in dup_indices]
                xx[cnt_msk] = [np.min(points_x[idx]) for idx in dup_indices]
                yy[cnt_msk] = [np.min(points_y[idx]) for idx in dup_indices]
                dup_stds = np.zeros(len(dup_indices))

            elif mode == 'max':
                zz[cnt_msk] = [np.max(pixel_z[idx]) for idx in dup_indices]
                xx[cnt_msk] = [np.max(points_x[idx]) for idx in dup_indices]
                yy[cnt_msk] = [np.max(points_y[idx]) for idx in dup_indices]
                dup_stds = np.zeros(len(dup_indices))

            elif mode == 'mean':
                zz[cnt_msk] = [np.mean(pixel_z[idx]) for idx in dup_indices]
                xx[cnt_msk] = [np.mean(points_x[idx]) for idx in dup_indices]
                yy[cnt_msk] = [np.mean(points_y[idx]) for idx in dup_indices]
                dup_stds = [np.std(pixel_z[idx]) for idx in dup_indices]
                
            elif mode == 'median':
                zz[cnt_msk] = [np.median(pixel_z[idx]) for idx in dup_indices]
                xx[cnt_msk] = [np.mean(points_x[idx]) for idx in dup_indices]
                yy[cnt_msk] = [np.mean(points_y[idx]) for idx in dup_indices]
                dup_stds = [np.std(pixel_z[idx]) for idx in dup_indices]

            elif mode == 'std':
                zz[cnt_msk] = [np.std(pixel_z[idx]) for idx in dup_indices]
                xx[cnt_msk] = [np.mean(points_x[idx]) for idx in dup_indices]
                yy[cnt_msk] = [np.mean(points_y[idx]) for idx in dup_indices]
                dup_stds = np.zeros(len(dup_indices))

            elif mode == 'var':
                zz[cnt_msk] = [np.var(pixel_z[idx]) for idx in dup_indices]
                xx[cnt_msk] = [np.mean(points_x[idx]) for idx in dup_indices]
                yy[cnt_msk] = [np.mean(points_y[idx]) for idx in dup_indices]
                dup_stds = np.zeros(len(dup_indices))

            elif mode == 'sums':
                zz[cnt_msk] = [np.sum(pixel_z[idx] * pixel_w[idx] * weight) for idx in dup_indices]
                xx[cnt_msk] = [np.sum(points_x[idx] * pixel_w[idx] * weight) for idx in dup_indices]
                yy[cnt_msk] = [np.sum(points_y[idx] * pixel_w[idx] * weight) for idx in dup_indices]
                ww[cnt_msk] = [np.sum(pixel_w[idx] * weight) for idx in dup_indices]
                dup_stds = [np.std(pixel_z[idx]) for idx in dup_indices]

            ## Aggregate uncertainty
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
        
        ## Uncertainty
        out_arrays['uncertainty'] = fill_grid(
            np.sqrt(uu**2 + (uncertainty)**2), fill_val=0.0
        )

        ## Weights
        out_arrays['weight'] = np.ones(grid_shape)
        if mode == 'sums':
            out_arrays['weight'][unq[:, 0], unq[:, 1]] = ww
        else:
            out_arrays['weight'][:] = weight
            out_arrays['weight'][unq[:, 0], unq[:, 1]] *= (ww * unq_cnt)

        ## Helper coords for calling class to map back
        out_arrays['pixel_x'] = local_px
        out_arrays['pixel_y'] = local_py

        return out_arrays, this_srcwin, self.dst_gt


class PointPixels_:
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
        self.ppm = ppm 
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

            
    def __call__(self, points, weight=1, uncertainty=0.0, mode='mean'):
        """Process points into a gridded array.
        
        Args:
            points (np.recarray): Input data containing 'x', 'y', 'z'.
            weight (float): Global weight multiplier.
            uncertainty (float): Global uncertainty value.
            mode (str): Aggregation mode ('mean', 'min', 'max', 'sums').
        """

        ## mrl: removed 'mask': None
        out_arrays = {
            'z': None, 'count': None, 'weight': None, 'uncertainty': None,
            'x': None, 'y': None, 'pixel_x': None, 'pixel_y': None
        }

        if points is None or len(points) == 0:
            return out_arrays, None, None

        ## Ensure region and geotransform are set
        if self.src_region is None:
            self.init_region_from_points(points)
        elif self.dst_gt is None:
            self.init_gt()

        ## Sanitize inputs
        weight = utils.float_or(weight, 1.)
        uncertainty = utils.float_or(uncertainty, 0.)

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
            np.sqrt(uu**2 + (uncertainty)**2), fill_val=0.0
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
                utils.echo_msg(f'Passed {len(self.points[~outliers_mask])} records.')
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
            pa = PointPixels(x_size=x_size, y_size=y_size, ppm=False)
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
class VectorMask(PointZ):
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
class OutlierZ(PointZ):
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
            if np.any(np.isfinite(residuals)):
                outlier_threshold = np.nanpercentile(residuals, percentile)
            else:
                outlier_threshold = np.nan

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


class RQOutlierZ(OutlierZ):
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

        if self.verbose:
            utils.echo_msg('Using Raster: {self.raster}')

            
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
        if raster and isinstance(raster, str):
            if os.path.exists(raster) and os.path.isfile(raster):
                return [raster]

            elif any(raster in item for item in self.fetches_modules):
                _raster = [item for item in self.fetches_modules if raster in item][0]
                this_fetch = self.fetch_data(_raster, self.region.copy().buffer(pct=1))
                raster = [x[1] for x in this_fetch.results]
                
                if self.xyinc is not None and self.resample_raster:
                    with warnings.catch_warnings():
                        warnings.simplefilter('ignore')
                        raster = [gdalfun.sample_warp(
                            raster, None, self.xyinc[0], self.xyinc[1],
                            sample_alg='bilinear', src_region=self.region,
                            verbose=False,
                            co=["COMPRESS=DEFLATE", "TILED=YES"]
                        )[0]]
                return raster
            
        elif raster is None:
            if (self.region is not None or self.xyinc is not None):# and self.resample_raster:
                _raster = utils.append_fn(
                    'rq_merged', self.region, 
                    self.xyinc[0], res=1 if not all(self.xyinc) else None
                )
                _raster = os.path.join(self.cache_dir, f'{_raster}.tif')
            else:
                _raster = utils.make_temp_fn('rq_merged.tif', self.cache_dir)


            # If the merged raster from a previous run exists, return it immediately.
            if os.path.exists(_raster) and os.path.isfile(_raster):
                return [_raster]
                
            raster = []
            ## Try gmrt all
            this_fetch = self.fetch_data(
                'gmrt', self.region.copy().buffer(pct=1)
            )
            raster_ = [x[1] for x in this_fetch.results]
            raster.extend([gdalfun.gmt_grd2gdal(x, verbose=False) \
                           if x.split('.')[-1] == 'grd' else x for x in raster_])
            
            ## Try etopo
            this_fetch = self.fetch_data(
                'etopo:datatype=surface', self.region.copy().buffer(pct=1)
            )
            raster.extend([x[1] for x in this_fetch.results])

            ## Try gmrt swath
            this_fetch = self.fetch_data('gmrt', self.region.copy().buffer(pct=1))
            raster_ = [x[1] for x in this_fetch.results]
            raster_ = [gdalfun.gmt_grd2gdal(x, verbose=False) \
                       if x.split('.')[-1] == 'grd' else x for x in raster_]
            gmrt_swath = self.mask_gmrt(raster_[0])
            if gmrt_swath is not None:
                raster.extend([gmrt_swath])

            ## Try cudem 1/3
            this_fetch = self.fetch_data(
                'CUDEM:datatype=13:keep_footprints=True', self.region.copy().buffer(pct=1)
            )
            raster.extend([x[1] for x in this_fetch.results])        

            ## Try cudem 1/9
            this_fetch = self.fetch_data(
                'CUDEM:datatype=19:keep_footprints=True', self.region.copy().buffer(pct=1)
            )
            raster.extend([x[1] for x in this_fetch.results])                     

        else:
            utils.echo_warning_msg(f'could not parse rq raster {raster}')
        return raster

    
    def init_raster_fallback(self, raster):
        """Initialize the reference raster list."""
        
        if raster and os.path.exists(raster) and os.path.isfile(raster):
            return [raster]
        
        ## Fetch Default GMRT
        rasters = []
        gmrt_fetch = self.fetch_data('gmrt')
        if gmrt_fetch:
            for r in gmrt_fetch.results:
                fn = r[1]
                if fn.endswith('.grd'): fn = gdalfun.gmt_grd2gdal(fn, verbose=False)
                #masked_fn = self.mask_gmrt(fn)
                masked_fn = None
                rasters.append(masked_fn if masked_fn else fn)

        ## Fallback
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
            if self.region is None:
                self.region = self.init_region()
            
            ## Ensure we check the cache_dir for the merged file
            tmp_raster_fn = utils.append_fn(
                'rq_merged', self.region, 
                self.xyinc[0] if self.xyinc[0] else 1, 
                res=1 if not all(self.xyinc) else None
            )
            tmp_raster_fn = os.path.join(self.cache_dir, f'{tmp_raster_fn}.tif')

            ## Create a temporary region buffered by 5%
            warp_region = self.region.copy()
            warp_region.buffer(pct=5) 

            x_res = self.xyinc[0] if self.xyinc else 0.0000925
            y_res = self.xyinc[1] if self.xyinc else 0.0000925

            try:
                ## Check if the file exists in the cache_dir before attempting to warp
                if os.path.exists(tmp_raster_fn):
                    ref_raster = tmp_raster_fn
                    temp_merged = ref_raster
                else:
                    warped = gdalfun.sample_warp(
                        self.raster, tmp_raster_fn, x_res, y_res,
                        sample_alg='bilinear', src_region=warp_region, 
                        verbose=self.verbose
                    )
                    if warped:
                        ref_raster = warped[0]
                        temp_merged = ref_raster
            except Exception as e:
                utils.echo_error_msg(f"Failed to merge rasters: {e}")
                return None

        ## Interpolate / Sample Directly
        ## We use direct array indexing to ensure perfect 1:1 alignment with input points.
        ds = gdal.Open(ref_raster)
        if ds is not None:
            gt = ds.GetGeoTransform()
            band = ds.GetRasterBand(1)
            ndv = band.GetNoDataValue()
            
            ## Load raster data (assuming chunk size fits in memory)
            raster_data = band.ReadAsArray()
            
            ## Convert Geo Coordinates to Pixel Indices
            ## px = (geo_x - origin_x) / pixel_width
            ## py = (geo_y - origin_y) / pixel_height
            px = ((points['x'] - gt[0]) / gt[1]).astype(int)
            py = ((points['y'] - gt[3]) / gt[5]).astype(int)
            
            ## Identify valid indices within raster bounds
            rows, cols = raster_data.shape
            valid_mask = (px >= 0) & (px < cols) & (py >= 0) & (py < rows)
            
            ## Initialize output array with NaNs
            sampled_z = np.full(len(points), np.nan)
            
            if np.any(valid_mask):
                ## Extract values using numpy advanced indexing
                z_vals = raster_data[py[valid_mask], px[valid_mask]]
                
                ## Filter NoData values
                if ndv is not None:
                    z_vals[z_vals == ndv] = np.nan
                    
                sampled_z[valid_mask] = z_vals
            
            ds = None
        else:
            ## Fallback if raster open fails
            sampled_z = np.full(len(points), np.nan)
        
        ## Cleanup: Only delete if it's a VRT (temp). 
        ## TIFs are kept for reuse.
        if temp_merged and os.path.exists(temp_merged):
            if temp_merged.endswith('.vrt'):
                try: os.remove(temp_merged)
                except OSError: pass

        ## Calculate differences
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
        
    
    # def point_residuals(self, points, percentage=True, res=50):
    #     """Calculate residuals against the loaded raster(s)."""
        
    #     if not self.raster:
    #         return None

    #     ## Handle Multiple Rasters via Warping/Merging
    #     ref_raster = self.raster[0]
    #     temp_merged = None
        
    #     if len(self.raster) > 1 or self.resample_raster:
    #         if self.region is None: self.region = self.init_region()
            
    #         ## Ensure we check the cache_dir for the merged file
    #         tmp_raster_fn = utils.append_fn(
    #             'rq_merged', self.region, 
    #             self.xyinc[0] if self.xyinc[0] else 1, 
    #             res=1 if not all(self.xyinc) else None
    #         )
    #         tmp_raster_fn = os.path.join(self.cache_dir, f'{tmp_raster_fn}.tif')

    #         ## Create a temporary region buffered by 5%
    #         warp_region = self.region.copy()
    #         warp_region.buffer(pct=5) 

    #         x_res = self.xyinc[0] if self.xyinc else 0.0000925
    #         y_res = self.xyinc[1] if self.xyinc else 0.0000925

    #         try:
    #             ## Check if the file exists in the cache_dir before attempting to warp
    #             if os.path.exists(tmp_raster_fn):
    #                 ref_raster = tmp_raster_fn
    #                 temp_merged = ref_raster
    #             else:
    #                 warped = gdalfun.sample_warp(
    #                     self.raster, tmp_raster_fn, x_res, y_res,
    #                     sample_alg='bilinear', src_region=warp_region, 
    #                     verbose=self.verbose
    #                 )
    #                 if warped:
    #                     ref_raster = warped[0]
    #                     temp_merged = ref_raster
    #         except Exception as e:
    #             utils.echo_error_msg(f"Failed to merge rasters: {e}")
    #             return None

    #     ## Interpolate
    #     sampled_z = gdalfun.gdal_query(points, ref_raster, 'g').flatten()
        
    #     ## Cleanup: Only delete if it's a VRT (temp). 
    #     ## TIFs are kept for reuse.
    #     if temp_merged and os.path.exists(temp_merged):
    #         if temp_merged.endswith('.vrt'):
    #             try: os.remove(temp_merged)
    #             except OSError: pass

    #     ## Validation check to catch alignment errors early
    #     if len(sampled_z) != len(points):
    #         if self.verbose:
    #             utils.echo_warning_msg(
    #                 f"Shape mismatch in RQ filter: Points {len(points)} vs Sampled {len(sampled_z)}. "
    #                 "Padding with NaNs."
    #             )
    #         new_sampled = np.full(len(points), np.nan)
    #         limit = min(len(points), len(sampled_z))
    #         new_sampled[:limit] = sampled_z[:limit]
    #         sampled_z = new_sampled
                
    #     diff = points['z'] - sampled_z
        
    #     if percentage:
    #         with np.errstate(divide='ignore', invalid='ignore'):
    #             if self.scaled_percentile:
    #                 vals = np.abs(diff / (points['z'] + sampled_z)) * 100
    #             else:
    #                 vals = np.abs(diff / sampled_z) * 100
    #             vals[~np.isfinite(vals)] = 0
    #         return vals
    #     else:
    #         return np.abs(diff)
        

class BlockThin(PointZ):
    """Thin point cloud by keeping one representative point per grid cell.
    Useful for reducing data volume or shoal-biasing (hydrography).
    
    Config: <blockthin:res=10:mode=min>
    """
    
    def __init__(self, res=10, mode='min', **kwargs):
        super().__init__(**kwargs)
        self.res = utils.float_or(res, 10)
        self.mode = mode if mode in ['min', 'max', 'mean', 'median'] else 'min'

        
    def run(self):
        if self.points is None: return None
        ## Calculate grid indices
        x_idx = np.floor((self.points['x'] - np.min(self.points['x'])) / self.res).astype(np.int64)
        y_idx = np.floor((self.points['y'] - np.min(self.points['y'])) / self.res).astype(np.int64)
        
        ## Hash/Encode X,Y into a single unique ID for grouping
        width = np.max(x_idx) + 1
        grid_ids = y_idx * width + x_idx
        
        ## Sort points by Grid ID
        sort_indices = np.argsort(grid_ids)
        sorted_ids = grid_ids[sort_indices]
        
        ## Find unique block boundaries
        unique_ids, unique_indices = np.unique(sorted_ids, return_index=True)
        
        ## Indices of the points we want to KEEP
        keep_indices = []
        
        ## Iterate over blocks
        ## unique_indices points to the start of each block in the sorted array
        for i in range(len(unique_ids)):
            start = unique_indices[i]
            end = unique_indices[i+1] if i+1 < len(unique_ids) else len(sorted_ids)
            
            ## Get indices for points in this block
            block_indices = sort_indices[start:end]
            block_points = self.points[block_indices]
            
            ## Select winner based on mode
            if self.mode == 'min':
                ## Argmin of Z
                local_idx = np.argmin(block_points['z'])
                keep_indices.append(block_indices[local_idx])
                
            elif self.mode == 'max':
                ## Argmax of Z
                local_idx = np.argmax(block_points['z'])
                keep_indices.append(block_indices[local_idx])
                
            elif self.mode == 'median':
                med_z = np.median(block_points['z'])
                local_idx = np.argmin(np.abs(block_points['z'] - med_z))
                keep_indices.append(block_indices[local_idx])

        mask = np.ones(len(self.points), dtype=bool) # All True (Remove all)
        mask[keep_indices] = False # Set Keepers to False
        
        return mask # Returns True for points to DELETE
    

## Dlim does this automatically with the -R region set with z-values,
## but this may be useful if using pointz outside of dlim/waffles
class RangeZ(PointZ):
    """Filter points based on vertical Z-range.
    
    By default, keeps points WITHIN the specified range [min_z, max_z].
    Use invert=True to remove points within the range.
    
    Config: <rangez:min_z=-100:max_z=0:invert=False>
    """
    
    def __init__(self, min_z=None, max_z=None, invert=False, **kwargs):
        super().__init__(**kwargs)
        self.min_z = utils.float_or(min_z)
        self.max_z = utils.float_or(max_z)
        self.invert = invert

        
    def run(self):
        """Generate a mask of outliers."""
        if self.points is None or len(self.points) == 0:
            return None
            
        z_vals = self.points['z']
        
        ## Start with all False (Keep everything)
        outliers = np.zeros(len(z_vals), dtype=bool)
        
        ## "Keep Inside" (Standard)
        ## Outliers are those OUTSIDE the range
        if not self.invert:
            if self.min_z is not None:
                outliers |= (z_vals < self.min_z)
            if self.max_z is not None:
                outliers |= (z_vals > self.max_z)
                
        ## "Remove Inside" (Invert)
        ## Outliers are those INSIDE the range
        else:
            inside_mask = np.ones(len(z_vals), dtype=bool)
            
            if self.min_z is not None:
                inside_mask &= (z_vals >= self.min_z)
            if self.max_z is not None:
                inside_mask &= (z_vals <= self.max_z)
                
            outliers = inside_mask

        return outliers


class CoplanarZ(PointZ):
    """Filter points that deviate from a locally fitted plane.
    Useful for removing noise from generally flat features (roads, water, plains).

    Config: <coplanar:radius=10:threshold=0.5:min_neighbors=3:invert=False>
    """
    
    def __init__(self, radius=10, threshold=0.5, min_neighbors=3, invert=False, **kwargs):
        super().__init__(**kwargs)
        self.radius = utils.float_or(radius, 10)
        self.threshold = utils.float_or(threshold, 0.5)
        self.min_neighbors = utils.int_or(min_neighbors, 3)
        self.invert = invert

        
    def run(self):
        """Generate outlier mask based on plane fitting."""
        
        if self.points is None or len(self.points) == 0:
            return None
            
        try:
            from scipy.spatial import cKDTree
        except ImportError:
            utils.echo_error_msg("scipy.spatial.cKDTree required for coplanar filter.")
            return None

        ## Build KDTree for efficient neighbor search
        coords = np.column_stack((self.points['x'], self.points['y']))
        tree = cKDTree(coords)
        
        ## Query neighbors within radius
        if self.verbose:
            utils.echo_msg(f"Querying neighbors (radius={self.radius})...")
        indices_list = tree.query_ball_point(coords, self.radius)
        
        outliers = np.zeros(len(self.points), dtype=bool)
        z_vals = self.points['z']
        
        ## Iterate points and fit planes
        with utils.ccp(total=len(self.points), desc='Plane Fitting', leave=False) as pbar:
            for i, neighbors in enumerate(indices_list):
                pbar.update()
                
                ## Check neighbor count (including self)
                if len(neighbors) < self.min_neighbors + 1:
                    # Treat isolated points as outliers (noise)
                    outliers[i] = True 
                    continue
                
                ## Get neighbor coordinates
                nb_coords = coords[neighbors]
                nb_z = z_vals[neighbors]
                
                center_x, center_y = coords[i]
                
                ## Setup Least Squares: Z = a*X + b*Y + c
                ## A matrix columns: [x_rel, y_rel, 1]
                A = np.column_stack((
                    nb_coords[:, 0] - center_x,
                    nb_coords[:, 1] - center_y,
                    np.ones(len(neighbors))
                ))
                
                try:
                    ## Fit plane
                    ## c (coeffs[2]) is the fitted Z at (0,0) relative coordinates (the query point)
                    coeffs, residuals, rank, s = np.linalg.lstsq(A, nb_z, rcond=None)
                    
                    fitted_z = coeffs[2] # Intercept at relative (0,0)
                    
                    ## Calculate deviation of the point from the fitted plane
                    deviation = abs(z_vals[i] - fitted_z)
                    
                    if deviation > self.threshold:
                        outliers[i] = True
                        
                except np.linalg.LinAlgError:
                    ## Collinear points or singular matrix -> Outlier
                    outliers[i] = True

        return ~outliers if self.invert else outliers
    

class BlockMinMax(PointZ):
    """Thin point cloud by keeping only the Min or Max Z point per grid block.
    Commonly used in hydrography for "shoal-biased" thinning.
    
    Config: <block_minmax:res=10:mode=min:invert=False>
    """
    
    def __init__(self, res=10, mode='min', invert=False, **kwargs):
        super().__init__(**kwargs)
        self.res = utils.float_or(res, 10)
        self.mode = mode.lower() if mode else 'min'
        self.invert = invert

        
    def run(self):
        """Generate mask of points to remove (non-min/max points)."""
        
        if self.points is None or len(self.points) == 0:
            return None
            
        ## Calculate Block Indices
        x_idx = np.floor(self.points['x'] / self.res).astype(np.int64)
        y_idx = np.floor(self.points['y'] / self.res).astype(np.int64)
        z_vals = self.points['z']
        
        ## Sort Data
        if self.mode == 'max':
            ## For Max, we want Z descending, so we sort by -Z
            sort_order = np.lexsort((-z_vals, y_idx, x_idx))
        else:
            # For Min, we want Z ascending
            sort_order = np.lexsort((z_vals, y_idx, x_idx))
            
        ## Apply sort to indices
        sorted_x = x_idx[sort_order]
        sorted_y = y_idx[sort_order]
        
        ## Find Unique Blocks
        ## Identify where the block index changes (flag=True at start of new block)
        ## Prepend True for the first element
        change_mask = np.concatenate(
            ([True], (sorted_x[1:] != sorted_x[:-1]) | (sorted_y[1:] != sorted_y[:-1]))
        )
        
        ## Indices in the SORTED array that are the "Keepers" (Min/Max points)
        keeper_sorted_indices = np.nonzero(change_mask)[0]
        
        ## Map back to ORIGINAL array indices
        keeper_original_indices = sort_order[keeper_sorted_indices]
        
        ## Create Outlier Mask
        ## Initialize to True (Remove All)
        outliers = np.ones(len(self.points), dtype=bool)
        
        ## Set Keepers to False (Don't Remove)
        outliers[keeper_original_indices] = False
        
        ## If invert=True, we return the Keepers as "Points to Remove" (dropping the min/max).
        return ~outliers if self.invert else outliers


class RasterMask(PointZ):
    """Filter points using a raster mask.
    
    By default, KEEPS points where the raster value is non-zero (Inside).
    Use invert=True to REMOVE points where the raster value is non-zero.
    
    Config: <raster_mask:mask_fn=path/to/mask.tif:invert=False>
    """
    
    def __init__(self, mask_fn=None, invert=False, **kwargs):
        super().__init__(**kwargs)
        self.mask_fn = mask_fn
        self.invert = invert
        if self.verbose:
            utils.echo_msg(f'Raster masking with {mask_fn}')

            
    def run(self):
        """Generate mask of outliers."""
        if self.points is None or len(self.points) == 0:
            return None
            
        if self.mask_fn is None or not os.path.exists(self.mask_fn):
            utils.echo_error_msg(f"Mask file not found: {self.mask_fn}")
            return None

        ## Query the raster at point locations
        sampled_values = gdalfun.gdal_query(self.points, self.mask_fn, 'g').flatten()

        ## Determine "Inside" vs "Outside"
        ## We assume Non-Zero = Inside (Valid Mask Area)
        ## We assume Zero = Outside (Background/NoData)
        is_inside = (sampled_values != 0)
        
        ## Determine Outliers (Points to Remove)
        ## Default: Remove Outside (Keep Inside) -> outliers = ~is_inside
        ## Invert: Remove Inside (Keep Outside) -> outliers = is_inside        
        if self.invert:
            ## Remove "Inside" points (e.g. masking OUT land)
            outliers = is_inside
        else:
            ## Remove "Outside" points (e.g. cropping TO a region)
            outliers = ~is_inside
            
        if self.verbose:
            action = "Removed" if self.invert else "Kept"
            count = np.count_nonzero(is_inside)
            utils.echo_msg(f"Mask analysis: {count} points inside mask area. ({action} inside)")

        return outliers
    

class DensityZ(PointZ):
    """Thin point cloud to a specific resolution/density.
    
    Modes:
    - random: Keep the first point found in the cell (fastest).
    - median: Keep the point with the median Z value.
    - mean: Keep the point with Z closest to the cell mean.
    - center: Keep the point closest to the 2D center of the cell.
    
    Config: <density:res=10:mode=random>
    """
    
    def __init__(self, res=10, mode='random', **kwargs):
        super().__init__(**kwargs)
        self.res = utils.float_or(res, 10)
        self.mode = mode.lower() if mode else 'random'

        
    def run(self):
        if self.points is None or len(self.points) == 0:
            return None
            
        ## Bin points into grid cells
        x_idx = np.floor(self.points['x'] / self.res).astype(np.int64)
        y_idx = np.floor(self.points['y'] / self.res).astype(np.int64)
        
        ## Sort based on Mode
        if self.mode == 'median':
            ## Sort by Block, then Z. The median index is roughly the middle.
            sort_keys = (self.points['z'], y_idx, x_idx)
        elif self.mode == 'mean':
            ## We need to calculate means first, which is expensive. 
            ## Approximate 'mean' by sorting by Z and picking middle (median) is usually sufficient.
            sort_keys = (self.points['z'], y_idx, x_idx)
        elif self.mode == 'center':
            ## Calculate distance to cell center:
            #center_x = (x_idx * res) + (res/2)
            dist = np.sqrt(
                (self.points['x'] - ((x_idx * self.res) + (self.res/2)))**2 + 
                (self.points['y'] - ((y_idx * self.res) + (self.res/2)))**2
            )
            sort_keys = (dist, y_idx, x_idx)
        else: # Random (or default)
            sort_keys = (y_idx, x_idx)

        ## Apply Sort
        sort_order = np.lexsort(sort_keys)
        sorted_x = x_idx[sort_order]
        sorted_y = y_idx[sort_order]
        
        ## Identify Blocks
        ## Find indices where the block ID changes
        change_mask = np.concatenate(
            ([True], (sorted_x[1:] != sorted_x[:-1]) | (sorted_y[1:] != sorted_y[:-1]))
        )
        
        ## Select Representatives
        block_start_indices = np.nonzero(change_mask)[0]
        
        ## If Random/Center, we just pick the first one (since we sorted by distance for center)
        if self.mode in ['random', 'center']:
            keeper_sorted_indices = block_start_indices
            
        ## If Median/Mean, we pick the middle index of the block
        elif self.mode in ['median', 'mean']:
            ## Calculate block sizes
            ## Append total length to calc size of last block
            block_ends = np.concatenate((block_start_indices[1:], [len(self.points)]))
            block_sizes = block_ends - block_start_indices
            
            ## The "middle" index relative to the start
            offsets = block_sizes // 2
            keeper_sorted_indices = block_start_indices + offsets

        ## Map back to original indices
        keeper_original_indices = sort_order[keeper_sorted_indices]
        
        ## Create Outlier Mask (True = Remove)
        outliers = np.ones(len(self.points), dtype=bool)
        outliers[keeper_original_indices] = False
        
        return outliers
    

class DiffZ(PointZ):
    """Filter points based on the signed difference from a reference raster.
    Diff = Point_Z - Raster_Z
    
    Useful for:
    - Bias filtering (e.g. remove points > 1m above reference)
    - Change detection (keep points within a specific change band)
    
    Config: <diff:raster=path.tif:min_diff=-5:max_diff=5:invert=False>
    """
    
    def __init__(self, raster=None, min_diff=None, max_diff=None, invert=False, **kwargs):
        super().__init__(**kwargs)
        self.raster = raster
        self.min_diff = utils.float_or(min_diff)
        self.max_diff = utils.float_or(max_diff)
        self.invert = invert

        
    def run(self):
        if self.points is None or len(self.points) == 0:
            return None
        
        if not self.raster or not os.path.exists(self.raster):
            utils.echo_error_msg(f"Reference raster not found: {self.raster}")
            return None

        ## Sample Raster
        sampled_z = gdalfun.gdal_query(self.points, self.raster, 'g').flatten()
        
        ## Calculate Signed Difference
        ## Positive diff = Point is ABOVE raster
        ## Negative diff = Point is BELOW raster
        diffs = self.points['z'] - sampled_z
        
        ## Determine "Inside Range" Mask
        ## We define "valid" as being BETWEEN min and max.
        keep_mask = np.ones(len(diffs), dtype=bool)
        
        if self.min_diff is not None:
            keep_mask &= (diffs >= self.min_diff)
            
        if self.max_diff is not None:
            keep_mask &= (diffs <= self.max_diff)
            
        ## Handle Outliers
        ## If invert=False (default): We want to KEEP points inside the range.
        ## So outliers are those NOT in the keep_mask.
        if not self.invert:
            outliers = ~keep_mask
        else:
            ## If invert=True: We want to REMOVE points inside the range.
            outliers = keep_mask
            
        if self.verbose:
            utils.echo_msg(f"Diff filter: {np.count_nonzero(outliers)} points flagged.")
            
        return outliers

    
## ==============================================
## Command-line Interface (CLI)
## $ pointz
##
## pointz cli
## ==============================================
class PointFilterFactory(factory.CUDEMFactory):
    _modules = {
        'outlierz': {'name': 'outlierz', 'call': OutlierZ},
        'rq': {'name': 'rq', 'call': RQOutlierZ},
        'vector_mask': {'name': 'vector_mask', 'call': VectorMask},
        'raster_mask': {'name': 'raster_mask', 'call': RasterMask},
        'block_thin': {'name': 'block_thin', 'call': BlockThin},
        'rangez': {'name': 'rangez', 'call': RangeZ},
        'coplanar': {'name': 'coplanar', 'call': CoplanarZ},
        'block_minmax': {'name': 'block_minmax', 'call': BlockMinMax},
        'density': {'name': 'density', 'call': DensityZ},
        'diff': {'name': 'diff', 'call': DiffZ},
    }

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        
class PrintModulesAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        factory.echo_modules(PointFilterFactory._modules, values, md=True if not values else False)
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
