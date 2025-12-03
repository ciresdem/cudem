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
###############################################################################
### Commentary:
##
## This grits module will blend the input data with the input raster.
##
### Code:

import numpy as np
import scipy
from scipy import interpolate
from tqdm import trange
from tqdm import tqdm
from cudem import utils
from cudem import gdalfun
from cudem.grits import grits

def interpolate_array(
        points_arr, ndv=-9999, chunk_size=None, chunk_step=None, chunk_buffer=0,
        method='cubic', verbose=True
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

    if verbose:
        utils.echo_msg(
            f'buffering srcwin by {chunk_buffer} pixels'
        )

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

        points_array = points_arr[srcwin_buff[1]:srcwin_buff[1]+srcwin_buff[3],
                                  srcwin_buff[0]:srcwin_buff[0]+srcwin_buff[2]]
        
        points_array[points_array == ndv] = np.nan
        point_indices = np.nonzero(~np.isnan(points_array))
        
        y_origin = srcwin[1]-srcwin_buff[1]
        x_origin = srcwin[0]-srcwin_buff[0]
        y_size = y_origin + srcwin[3]
        x_size = x_origin + srcwin[2]
        
        if np.all(np.isnan(points_array)):
            continue
        elif np.count_nonzero(np.isnan(points_array)) == 0:

            points_array = points_array[y_origin:y_size,x_origin:x_size]
            interp_array[srcwin[1]:srcwin[1]+points_array.shape[0],
                         srcwin[0]:srcwin[0]+points_array.shape[1]] = points_array

        elif len(point_indices[0]):
            point_values = points_array[point_indices]
            xi, yi = np.mgrid[0:srcwin_buff[2],
                              0:srcwin_buff[3]]

            interp_data = interpolate.griddata(
                np.transpose(point_indices), point_values,
                (xi, yi), method=method
            )
            interp_data = interp_data[y_origin:y_size,x_origin:x_size]
            interp_array[srcwin[1]:srcwin[1]+interp_data.shape[0],
                         srcwin[0]:srcwin[0]+interp_data.shape[1]] = interp_data

    point_values = None

    return(interp_array)


class Blend(grits.Grits):
    def __init__(self,  weight_threshold=None, buffer_cells=1, gap_fill_cells=None,
                 weight_thresholds=None, buffer_sizes=None, gap_fill_sizes=None,
                 binary_dilation=True, binary_pulse=False,
                 fill_holes=False, aux_dems=None, random_buffer=False, **kwargs):
        super().__init__(**kwargs)
        self.weight_threshold = utils.float_or(weight_threshold, 1)
        self.buffer_cells = utils.int_or(buffer_cells, 1)
        self.gap_fill_cells = utils.int_or(gap_fill_cells, 0)
        
        self.weight_thresholds = weight_thresholds
        self.buffer_sizes = buffer_sizes
        self.gap_fill_sizes = gap_fill_sizes
        
        self.binary_dilation = binary_dilation
        self.binary_pulse = binary_pulse
        
        self.fill_holes = fill_holes
        self.random_buffer = random_buffer

        self.aux_dems = aux_dems.split(',')
        

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
            
        return(closed_and_dilated_arr)

        
    def init_data(self):
        from cudem.datalists.dlim import DatasetFactory
        
        aux_arrs = []
        src_arr = None
        combined_arr = None

        src_region, ds_config = self.init_region(self.src_dem)
        x_inc=ds_config['geoT'][1]
        y_inc=ds_config['geoT'][5]*-1
        src_region.buffer(
            x_bv=(x_inc*self.buffer_cells),
            y_bv=(y_inc*self.buffer_cells)
        )

        xcount, ycount, dst_gt = src_region.geo_transform(
            x_inc=ds_config['geoT'][1], y_inc=ds_config['geoT'][5], node='grid'
        )
        combined_arr = np.full((ycount, xcount), np.nan)

        for aux_fn in self.aux_dems:
            aux_ds = DatasetFactory(
                mod=aux_fn, src_region=src_region, x_inc=x_inc, y_inc=y_inc
            )._acquire_module().initialize()
            for arrs, srcwin, gt in aux_ds.array_yield:
                ## supercede
                combined_arr[srcwin[1]:srcwin[1]+srcwin[3],
                             srcwin[0]:srcwin[0]+srcwin[2]] = arrs['z']

        combined_mask = ~np.isnan(combined_arr)
        if self.buffer_cells is not None:
            combined_mask = self.binary_closed_dilation(combined_mask, iterations=self.buffer_cells)
            
        ## load the src dem
        src_arr = np.full((ycount, xcount), np.nan)        
        with gdalfun.gdal_datasource(self.src_dem) as src_ds:
            if src_ds is not None:
                self.init_ds(src_ds)
                srcwin = (self.buffer_cells, self.buffer_cells, src_ds.RasterXSize, src_ds.RasterYSize)
                src_arr[srcwin[1]:srcwin[1]+srcwin[3],
                        srcwin[0]:srcwin[0]+srcwin[2]] = self.ds_band.ReadAsArray()

                src_arr[src_arr == ds_config['ndv']] = np.nan
                src_mask = ~np.isnan(src_arr)
                
        return(combined_arr, src_arr, combined_mask, src_mask)
        

    def blend_data(self):
        ## combined_arr is the aux_data, src_arr is the src dem data
        ## combined mask is the combined data + buffer
        ## src_mask is the src data
        combined_arr, src_arr, combined_mask, src_mask = self.init_data()
            
        ## combined_arr now has a nan buffer between aux and src data               
        combined_arr[~combined_mask] = src_arr[~combined_mask]

        ## generate a distance transform and normalize the results
        ## to values between 0 (near combined) and 1 (near src)
        dt = scipy.ndimage.distance_transform_cdt(combined_mask, metric='taxicab')
        dt = dt[(src_mask) & (combined_mask)]
        dt = (dt - np.min(dt)) / (np.max(dt) - np.min(dt))
        
        ## random src mask
        if self.random_buffer:
            utils.echo_msg('rand!')
            random_mask = np.random.rand(combined_arr.shape[0], combined_arr.shape[1]) > 0.985
            combined_arr[(combined_mask) & (src_mask) & (random_mask)] = src_arr[(combined_mask) & (src_mask) & (random_mask)]
        
        ## interpolate the buffer and extract just the buffer area and calculate the
        ## difference between the src data and the interp data
        interp_arr = interpolate_array(
            combined_arr, ndv=self.ds_config['ndv'], method='linear',
            chunk_size=self.buffer_cells*5, chunk_step=self.buffer_cells*5,
            #chunk_buffer=self.buffer_cells
        )
        interp_data_in_buffer = interp_arr[(combined_mask) & (src_mask)]
        src_data_in_buffer = src_arr[(combined_mask) & (src_mask)]
        interp_arr = None
        buffer_diffs = interp_data_in_buffer - src_data_in_buffer
        buffer_diffs[np.isnan(buffer_diffs)] = 0

        ## testings
        # med_buffer_diffs = np.median(buffer_diffs)
        # diff_threshold = .5
        # buffer_diffs[np.abs(buffer_diffs) <= diff_threshold] = 0
        # buffer_diffs[np.abs(buffer_diffs) > diff_threshold] = med_buffer_diffs
        
        ## apply the normalize results to the differences and set and write out the results
        combined_arr[(combined_mask) & (src_mask)] = src_arr[(combined_mask) & (src_mask)] + (buffer_diffs*dt)
        combined_arr[~src_mask] = np.nan
        combined_arr[np.isnan(combined_arr)] = self.ds_config['ndv']
        srcwin = (self.buffer_cells, self.buffer_cells, self.ds_config['nx'], self.ds_config['ny'])
        
        return(combined_arr, srcwin)

        
    def run(self):
        dst_ds = self.copy_src_dem()
        if dst_ds is None:
            return(self.src_dem, -1)

        combined_arr, srcwin = self.blend_data()        
        this_band = dst_ds.GetRasterBand(1)
        this_band.WriteArray(combined_arr[srcwin[1]:srcwin[1]+srcwin[3],
                                          srcwin[0]:srcwin[0]+srcwin[2]])
        
        dst_ds = None
        return(self.dst_dem, 0)


### End
