### blend.py - GRId filTerS
##
## Copyright (c) 2024 - 2025 Regents of the University of Colorado
##
## blend.py is part of cudem
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
## This grits module will blend the input data with the input raster. Using a weight_threshold
## this will take all cells over the weight threshold (either from a weight band in the input
## or a seperate raster) and buffer around those cells. Within that buffer it will interpolate
## from the higher weighted cells to the lower weighted cells and apply those interpolated
## values to the source raster tapering from the higher weighted cells to the edge of the buffer.
##
### Code:

import numpy as np
import scipy
from scipy import interpolate
from cudem import utils
from cudem import gdalfun
from cudem.grits import grits


def interpolate_array(
        points_arr, ndv=-9999, chunk_size=None, chunk_step=None,
        chunk_buffer=0, method='cubic', verbose=True
):

    ycount, xcount = points_arr.shape
    if chunk_size is None:
        n_chunk = int(xcount * .05)
        n_chunk = 10 if n_chunk < 10 else n_chunk
    else:
        n_chunk = chunk_size

    if chunk_step is None:
        #n_step = int(n_chunk/2)
        n_step = n_chunk
    else:
        n_step = chunk_step

    if verbose and chunk_buffer > 0:
        utils.echo_msg(f'buffering srcwin by {chunk_buffer} pixels')

    interp_array = np.full(points_arr.shape, np.nan)
    for srcwin in utils.yield_srcwin(
            (ycount, xcount),
            n_chunk=n_chunk,
            msg='generating {} interpolated array'.format(method),
            verbose=verbose,
            step=n_step
    ):
        srcwin_buff = utils.buffer_srcwin(
            srcwin, (ycount, xcount), chunk_buffer
        )

        points_array = points_arr[
            srcwin_buff[1]:srcwin_buff[1] + srcwin_buff[3],
            srcwin_buff[0]:srcwin_buff[0] + srcwin_buff[2]
        ]
        
        points_array[points_array == ndv] = np.nan
        point_indices = np.nonzero(~np.isnan(points_array))
        
        y_origin = srcwin[1] - srcwin_buff[1]
        x_origin = srcwin[0] - srcwin_buff[0]
        y_size = y_origin + srcwin[3]
        x_size = x_origin + srcwin[2]
        
        if np.all(np.isnan(points_array)):
            continue
        
        elif np.count_nonzero(np.isnan(points_array)) == 0:
            points_array = points_array[y_origin:y_size,x_origin:x_size]
            interp_array[
                srcwin[1]:srcwin[1]+points_array.shape[0],
                srcwin[0]:srcwin[0]+points_array.shape[1]
            ] = points_array

        elif len(point_indices[0]) <= 3:
            # elif len(point_indices[0]):
            continue
        
        else:
            point_values = points_array[point_indices]
            xi, yi = np.mgrid[
                0:srcwin_buff[2],
                0:srcwin_buff[3]
            ]

            interp_data = interpolate.griddata(
                np.transpose(point_indices),
                point_values,
                (xi, yi),
                method=method
            )
            interp_data = interp_data[y_origin:y_size, x_origin:x_size]
            interp_array[
                srcwin[1]:srcwin[1] + interp_data.shape[0],
                srcwin[0]:srcwin[0] + interp_data.shape[1]
            ] = interp_data

    point_values = None    
    return interp_array


def generate_slope_array(
        points_arr, ndv=-9999, chunk_size=None, chunk_step=None,
        chunk_buffer=0, normalize=True, buffer_mask=None,
        verbose=True
):

    ycount, xcount = points_arr.shape
    if chunk_size is None:
        n_chunk = int(xcount * .05)
        n_chunk = 10 if n_chunk < 10 else n_chunk
    else:
        n_chunk = chunk_size

    if chunk_step is None:
        #n_step = int(n_chunk/2)
        n_step = n_chunk
    else:
        n_step = chunk_step

    if verbose and chunk_buffer > 0:
        utils.echo_msg(f'buffering srcwin by {chunk_buffer} pixels')

    slope_array = np.full(points_arr.shape, np.nan)
    for srcwin in utils.yield_srcwin(
            (ycount, xcount),
            n_chunk=n_chunk,
            msg='generating slope array',
            verbose=verbose,
            step=n_step
    ):
        srcwin_buff = utils.buffer_srcwin(
            srcwin, (ycount, xcount), chunk_buffer
        )

        points_array = points_arr[
            srcwin_buff[1]:srcwin_buff[1] + srcwin_buff[3],
            srcwin_buff[0]:srcwin_buff[0] + srcwin_buff[2]
        ]
        
        points_array[points_array == ndv] = np.nan
        #point_indices = np.nonzero(~np.isnan(points_array))
        
        y_origin = srcwin[1] - srcwin_buff[1]
        x_origin = srcwin[0] - srcwin_buff[0]
        y_size = y_origin + srcwin[3]
        x_size = x_origin + srcwin[2]
        
        if np.all(np.isnan(points_array)):
            continue

        gy, gx = np.gradient(points_array)
        slope_data = np.sqrt(gx * gx + gy * gy)
        
        if slope_data.size == 0:
            continue
        
        slope_data = slope_data[y_origin:y_size, x_origin:x_size]
        slope_array[
            srcwin[1]:srcwin[1] + slope_data.shape[0],
            srcwin[0]:srcwin[0] + slope_data.shape[1]
        ] = slope_data


    vals = slope_array[buffer_mask]
    vals = vals[np.isfinite(vals)]

    if vals.size == 0:
        return None

    m = np.nanmax(np.abs(vals))
    if m == 0 or not np.isfinite(m):
        return None

    slope_norm = np.abs(slope_array[buffer_mask]) / m
    slope_norm[np.isnan(slope_norm)] = 0.0
    return slope_norm


class Blend(grits.Grits):
    """Blend data together

        slope_scale:
            0.0 -> disable slope gating (classic behavior)
            0–1 -> slope threshold in normalized slope units; higher = only
                   steeper regions allowed to keep random src points.

        random_scale:
            0.0 -> disable random point insertion
            0-1 -> higher = inlcude more random points
    """
    
    def __init__(
            self,
            weight_threshold=None,
            buffer_cells=1,
            binary_dilation=True,
            binary_pulse=False,
            sub_buffer_cells=0,
            slope_scale=.5,
            random_scale=.025,
            interp_method='linear',
            aux_dems=None,
            **kwargs
    ):        
        super().__init__(**kwargs)
        self.weight_threshold = utils.float_or(weight_threshold, 1)        
        self.buffer_cells = utils.int_or(buffer_cells, 1)
        
        self.binary_dilation = binary_dilation
        self.binary_pulse = binary_pulse
        
        self.sub_buffer_cells = utils.int_or(sub_buffer_cells, 0)
        self.slope_scale = np.clip(utils.float_or(slope_scale, 0.5), 0.0, 1.0)
        self.random_scale = np.clip(utils.float_or(random_scale, .025), 0.0, 1.0)
        if utils.str_or(aux_dems) is not None:
            self.aux_dems = aux_dems.split(',')
        else:
            self.aux_dems = []

        self.interp_method = 'linear'
        if interp_method in ['linear', 'cubic', 'nearest']:
            self.interp_method = interp_method

        
    def binary_closed_dilation(self, arr, iterations=1, closing_iterations=1):
        closed_and_dilated_arr = arr.copy()
        struct_ = scipy.ndimage.generate_binary_structure(2, 2)
        for i in range(closing_iterations):
            closed_and_dilated_arr = scipy.ndimage.binary_dilation(
                closed_and_dilated_arr, iterations=1, structure=struct_
            )
            closed_and_dilated_arr = scipy.ndimage.binary_erosion(
                closed_and_dilated_arr, iterations=1, border_value=1, structure=struct_
            )

        closed_and_dilated_arr = scipy.ndimage.binary_dilation(
            closed_and_dilated_arr, iterations=iterations, structure=struct_
        )
            
        return closed_and_dilated_arr


    def binary_reversion(self, arr, iterations, closing_iterations):
        reversion_arr = arr.copy()
        closing_iterations += iterations
            
        struct_ = scipy.ndimage.generate_binary_structure(2, 2)
        reversion_arr = scipy.ndimage.binary_dilation(
            reversion_arr, iterations=closing_iterations, structure=struct_
        )
        erosion_iterations = max(closing_iterations-iterations, 0)
        if erosion_iterations > 0:
            reversion_arr = scipy.ndimage.binary_erosion(
                reversion_arr, iterations=erosion_iterations, border_value=1, structure=struct_
            )

        return reversion_arr

    
    def init_weight(self, src_ds=None):
        weight_band = None
        if self.weight_is_fn:
            if self.weight_threshold is None:
                self.weight_threshold = gdalfun.gdal_percentile(
                    self.weight_mask, perc=75
                )

            weight_ds = gdal.Open(self.weight_mask)
            weight_band = weight_ds.GetRasterBand(1)

        elif self.weight_is_band and src_ds is not None:
            if self.weight_threshold is None:
                self.weight_threshold = gdalfun.gdal_percentile(
                    self.src_dem, perc=50, band=self.weight_mask
                )

            weight_band = src_ds.GetRasterBand(self.weight_mask)

        return weight_band

    
    def init_data(self):
        from cudem.datalists.dlim import DatasetFactory
        
        aux_arrs = []
        src_arr = None

        src_region, ds_config = self.init_region(self.src_dem)
        x_inc = ds_config['geoT'][1]
        y_inc = ds_config['geoT'][5] * -1

        src_region.buffer(
            x_bv=(x_inc * self.buffer_cells),
            y_bv=(y_inc * self.buffer_cells)
        )

        xcount, ycount, dst_gt = src_region.geo_transform(
            x_inc=ds_config['geoT'][1],
            y_inc=ds_config['geoT'][5],
            node='grid'
        )

        ## Load the src DEM.
        src_arr = np.full((ycount, xcount), np.nan)
        w_arr = np.full_like(src_arr, np.nan)
        #combined_arr = np.full_like(src_arr, np.nan)

        with gdalfun.gdal_datasource(self.src_dem) as src_ds:
            if src_ds is not None:
                self.init_ds(src_ds)
                srcwin = (self.buffer_cells, self.buffer_cells,
                          src_ds.RasterXSize, src_ds.RasterYSize)

                weight_band = self.init_weight(src_ds)
                if weight_band is not None:
                    w_arr[
                        srcwin[1]:srcwin[1] + srcwin[3],
                        srcwin[0]:srcwin[0] + srcwin[2]
                    ] = weight_band.ReadAsArray()
                    w_arr[w_arr == self.ds_config['ndv']] = np.nan
                    weights = np.unique(w_arr)[::-1] 
                
                src_arr[
                    srcwin[1]:srcwin[1] + srcwin[3],
                    srcwin[0]:srcwin[0] + srcwin[2]
                ] = self.ds_band.ReadAsArray()

                src_arr[src_arr == ds_config['ndv']] = np.nan

        #combined_arr = src_arr[w_arr >= self.weight_thresholds[0]]
        src_mask = np.isfinite(src_arr)

        ## Load aux data
        combined_arr = self.load_aux_data(aux_source_list=self.aux_data, src_region=src_region, src_gt=self.gt, x_count=xcount, y_count=ycount)
        if combined_arr is None:
            combined_arr = np.full_like(src_arr, np.nan)
            
        if np.any(w_arr):
            combined_arr = np.where(w_arr >= self.weight_threshold, src_arr, combined_arr)

        combined_mask = ~np.isnan(combined_arr)
                
        if self.sub_buffer_cells is not None:
            combined_sub_mask = self.binary_closed_dilation(
                combined_mask, iterations=self.sub_buffer_cells
            )
        else:
            combined_sub_mask = np.zeros_like(combined_mask, dtype=bool)
            
        if self.buffer_cells is not None:
            combined_mask = self.binary_closed_dilation(
                combined_mask, iterations=self.buffer_cells
            )
            
        return combined_arr, src_arr, combined_mask, src_mask, combined_sub_mask

    
    def _compute_slope(self, src_arr, buffer_mask=None):
        """Compute normalized slope over the whole array."""
        
        tmp = src_arr.copy()
        tmp[~np.isfinite(tmp)] = np.nan
        if np.all(np.isnan(tmp)):
            return(None)
        
        gy, gx = np.gradient(tmp)
        slope = np.sqrt(gx * gx + gy * gy)
        tmp = None
        
        vals = slope[buffer_mask]
        vals = vals[np.isfinite(vals)]

        if vals.size == 0:
            return(None)
        
        m = np.nanmax(np.abs(vals))
        if m == 0 or not np.isfinite(m):
            return(None)

        slope_norm = np.abs(slope[buffer_mask]) / m
        slope_norm[np.isnan(slope_norm)] = 0.0
        return slope_norm

    
    def blend_data(self):
        ## combined_arr is the aux_data, src_arr is the src dem data
        ## combined_mask is the combined data + buffer
        ## src_mask is the src data
        combined_arr, src_arr, combined_mask, src_mask, combined_sub_mask = self.init_data()
        buffer_mask = combined_mask & src_mask
        sub_buffer_mask = combined_sub_mask & src_mask   # core seam: no randomization
        slope_norm = None
        
        ## combined_arr now has a nan buffer between aux and src data               
        combined_arr[~combined_mask] = src_arr[~combined_mask]

        ## initial distance transform
        dt = scipy.ndimage.distance_transform_cdt(
            combined_mask, metric='taxicab'
        )

        random_arr = np.random.rand(*combined_arr.shape)
        random_mask = random_arr < self.random_scale  # base density of random picks

        ## Never randomize in core seam
        random_mask[sub_buffer_mask] = False

        ## if slope gating is enabled, only allow flips in steeper areas
        ## slope-aware random src mask in the OUTER buffer only
        if self.slope_scale > 0:
            # slope_norm = self._compute_slope(
            #     src_arr, buffer_mask=(buffer_mask & (~sub_buffer_mask))
            # )
            outer_mask = buffer_mask & (~sub_buffer_mask)
            slope_norm = generate_slope_array(
                src_arr,
                self.ds_config['ndv'],
                buffer_mask=outer_mask,
                chunk_size=self.buffer_cells * 5,
                chunk_step=self.buffer_cells * 5,
            )

            if slope_norm is not None:
                low_slope = slope_norm < self.slope_scale                
                #random_mask[low_slope & outer_mask] = False
                random_mask[outer_mask][low_slope] = False
                slope_norm = outer_mask = low_slope = None
                
        ## Only apply randomization where both src + combined_mask are valid (buffer region)
        valid_random = random_mask & buffer_mask

        ## Inject source values in the outer buffer; this helps preserve shape.
        combined_arr[valid_random] = src_arr[valid_random]
        combined_mask[valid_random] = True

        ## Recompute distance transform after randomization.
        dt = scipy.ndimage.distance_transform_cdt(
            combined_mask, metric='taxicab'
        )
        
        ## extract and normalize dt in the seam/buffer region
        dt_vals = dt[buffer_mask].astype(np.float32)
        if dt_vals.size > 0:
            dt_min = np.min(dt_vals)
            dt_max = np.max(dt_vals)
            if dt_max > dt_min:
                dt_norm = (dt_vals - dt_min) / (dt_max - dt_min)
            else:
                dt_norm = np.zeros_like(dt_vals, dtype=np.float32)
        else:
            dt_norm = np.zeros(0, dtype=np.float32)

        ## Interpolate the buffer and extract just the buffer area and calculate the
        ## difference between the src data and the interp data.
        interp_arr = interpolate_array(
            combined_arr,
            ndv=self.ds_config['ndv'],
            method=self.interp_method,
            chunk_size=self.buffer_cells * 5,
            chunk_step=self.buffer_cells * 2.5,
        )
        interp_data_in_buffer = interp_arr[buffer_mask]
        src_data_in_buffer = src_arr[buffer_mask]
        interp_arr = None

        buffer_diffs = interp_data_in_buffer - src_data_in_buffer
        buffer_diffs[np.isnan(buffer_diffs)] = 0

        ## Apply the normalized results to the differences and set and write out the results
        combined_arr[buffer_mask] = src_arr[buffer_mask] + (buffer_diffs * dt_norm)
        combined_arr[~src_mask] = np.nan
        combined_arr[np.isnan(combined_arr)] = self.ds_config['ndv']

        srcwin = (
            self.buffer_cells,
            self.buffer_cells,
            self.ds_config['nx'],
            self.ds_config['ny']
        )
        
        return combined_arr, srcwin

        
    def run(self):
        dst_ds = self.copy_src_dem()
        if dst_ds is None:
            return(self.src_dem, -1)

        combined_arr, srcwin = self.blend_data()        
        this_band = dst_ds.GetRasterBand(1)
        this_band.WriteArray(
            combined_arr[srcwin[1]:srcwin[1] + srcwin[3],
                         srcwin[0]:srcwin[0] + srcwin[2]]
        )

        dst_ds = None
        return self.dst_dem, 0


### End
