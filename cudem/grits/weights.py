### weights.py - GRId filTerS
##
## Copyright (c) 2024 - 2025 Regents of the University of Colorado
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
##
### Code:

import numpy as np
import scipy
from tqdm import trange
from tqdm import tqdm
from cudem import utils
from cudem import gdalfun
from cudem.grits import grits

class Weights(grits.Grits):
    """Create a nodata buffer around data at and above `weight_threshold`

    You must supply a weight-mask raster along with the DEM to filter.

    Parameters:

    weight_thresholds(float) - the weight threshold; separate levels with slash '/'
    buffer_sizes(int) - the number of cells to buffer; separate levels with slash '/'
    gap_fill_sizes(int) - the number of cells to gap-fill; separate levels with slash '/'
    fill_holes(bool) - fill holes in the weight level mask

    <weights:weight_thresholds=2/20:buffer_sizes=4/2:gap_fill_sizes=10/8>
    """
    
    def __init__(self,  weight_threshold=None, buffer_cells=1, gap_fill_cells=None,
                 weight_thresholds=None, buffer_sizes=None, gap_fill_sizes=None,
                 binary_dilation=True, binary_pulse=False,
                 fill_holes=False, **kwargs):
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
        self.init_thresholds()


    def init_thresholds(self):
        if self.weight_thresholds is not None:
            self.weight_thresholds = [float(x) for x in self.weight_thresholds.split('/')]

            if self.buffer_sizes is not None:
                self.buffer_sizes = [int(x) for x in self.buffer_sizes.split('/')]
            else:
                self.buffer_sizes = [self.buffer_cells]*len(self.weight_thresholds)

            if self.gap_fill_sizes is not None:
                self.gap_fill_sizes = [int(x) for x in self.gap_fill_sizes.split('/')]
            else:
                self.gap_fill_sizes = [self.gap_fill_cells]*len(self.weight_thresholds)
                
        else:
            self.weight_thresholds = [self.weight_threshold]
            self.buffer_sizes = [self.buffer_cells]
            self.gap_fill_sizes = [self.gap_fill_cells]

        if len(self.buffer_sizes) < len(self.weight_thresholds):
            len_diff = len(self.weight_thresholds) - len(self.buffer_sizes)
            self.buffer_sizes.extend([1] * len_diff)
            
        if len(self.gap_fill_sizes) < len(self.weight_thresholds):
            len_diff = len(self.weight_thresholds) - len(self.gap_fill_sizes)
            self.gap_fill_sizes.extend([0] * len_diff)

        utils.echo_msg([self.weight_thresholds, self.buffer_sizes, self.gap_fill_sizes])

            
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

        return(weight_band)


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

        return(reversion_arr)

        
    def run(self):
        if self.weight_mask is None:
            return(self.src_dem, -1)
        
        dst_ds = self.copy_src_dem()
        if dst_ds is None:
            return(self.src_dem, -1)
        
        with gdalfun.gdal_datasource(self.src_dem) as src_ds:
            if src_ds is not None:
                self.init_ds(src_ds)                
                weight_band = self.init_weight(src_ds)
                # srcwin = (0, 0, dst_ds.RasterXSize, dst_ds.RasterYSize)
                # with tqdm(
                #         total=srcwin[1] + srcwin[3],
                #         desc=f'buffering around weights over {self.weight_threshold}',
                #         leave=self.verbose
                # ) as pbar:
                #     for y in range(
                #             srcwin[1], srcwin[1] + srcwin[3], 1
                #     ):
                #         pbar.update()
                #n_chunk=self.ds_config['nx']/9, step=self.ds_config['nx']/27
                #n_chunk = max(self.buffer_cells**2, self.ds_config['ny']/9)
                #n_step = max(self.buffer_cells, self.ds_config['ny']/27)
                # n_chunk = self.ds_config['ny']
                # for srcwin in utils.yield_srcwin(
                #         (self.ds_config['ny'], self.ds_config['nx']),
                #         n_chunk=n_chunk, #step=n_step,
                #         verbose=self.verbose, start_at_edge=True,
                #         msg=f'buffering around weights over {self.weight_thresholds}'
                # ):
                #w_arr = weight_band.ReadAsArray(srcwin[0], y, srcwin[2], 1)
                w_arr = weight_band.ReadAsArray()#*srcwin)
                w_arr[w_arr == self.ds_config['ndv']] = np.nan
                weights = np.unique(w_arr)[::-1]
                this_w_arr = w_arr.copy()

                with tqdm(
                        total=len(self.weight_thresholds),
                        desc=f'buffering around weights over {self.weight_thresholds}',
                        leave=self.verbose
                ) as w_pbar:
                    for n, weight_threshold in enumerate(self.weight_thresholds):
                        this_w_arr[this_w_arr < weight_threshold] = np.nan
                        if self.binary_dilation:                            
                            if self.binary_pulse:
                                weight_mask = this_w_arr >= weight_threshold
                                expanded_w_arr = self.binary_closed_dilation(
                                    weight_mask, iterations=self.buffer_sizes[n],
                                    closing_iterations=self.gap_fill_sizes[n]
                                )                            
                            else:
                                weight_mask = this_w_arr >= weight_threshold
                                expanded_w_arr = self.binary_reversion(
                                    weight_mask, self.buffer_sizes[n],
                                    self.gap_fill_sizes[n]
                                )
                        else:
                            ## mrl use ndimage.binary_dilation rather than the slower expand_for below
                            shiftx = self.buffer_sizes[n] + self.gap_fill_sizes[n]#self.buffer_sizes[n] * self.revert_multiplier
                            shifty = self.buffer_sizes[n] + self.gap_fill_sizes[n]#self.buffer_sizes[n] * self.revert_multiplier
                            utils.echo_msg([shiftx, shifty])
                            expanded_w_arr = utils.expand_for(
                                this_w_arr >= weight_threshold,
                                shiftx=int(shiftx), shifty=int(shifty)
                            )

                            shiftx = max(0, shiftx - self.buffer_sizes[n])
                            shifty = max(0, shifty - self.buffer_sizes[n])
                            if shiftx > 0 and shifty > 0:
                                utils.echo_msg([shiftx, shifty])
                                contracted_w_arr = np.invert(utils.expand_for(
                                    np.invert(expanded_w_arr),
                                    shiftx=int(shiftx), shifty=int(shifty),
                                ))
                                expanded_w_arr = contracted_w_arr.copy()

                        if self.fill_holes:
                            expanded_w_arr = scipy.ndimage.binary_fill_holes(expanded_w_arr)
                            
                        mask = (w_arr < weight_threshold) & expanded_w_arr
                        this_w_arr = w_arr.copy()
                        this_w_arr[mask] = np.nan
                        ## mask out any values in other bands of the input/output raster
                        for b in range(1, dst_ds.RasterCount+1):
                            this_band = dst_ds.GetRasterBand(b)
                            this_arr = this_band.ReadAsArray()#*srcwin)
                            this_arr[mask] = self.ds_config['ndv']
                            this_band.WriteArray(this_arr)#, srcwin[0], srcwin[1])
                            this_band.FlushCache()
                            this_band = this_arr = None

                        w_pbar.update()

                if self.weight_is_fn:
                    weight_ds = None

            else:
                utils.echo_msg('failed')
                dst_ds = None
                return(self.src_dem, -1)

        dst_ds = None        
        return(self.dst_dem, 0)


## uses way too much memory on large dems
class WeightZones(Weights):
    def __init__(self, size_threshold=None, **kwargs):
        super().__init__(**kwargs)
        self.size_threshold = utils.float_or(size_threshold)

        
    def run(self):
        if self.weight_mask is None:
            return(self.src_dem, -1)

        dst_ds = self.copy_src_dem()
        with gdalfun.gdal_datasource(self.src_dem) as src_ds:
            if src_ds is not None:
                self.init_ds(src_ds)
                weight_band = self.init_weight(src_ds)                    
                this_w_arr = weight_band.ReadAsArray()
                #_mask = np.zeros(this_w_arr.shape)
                #_mask[(this_w_arr < self.weight_threshold)] = 1
                this_w_arr[this_w_arr == self.ds_config['ndv']] = np.nan
                n_den = self._density(this_w_arr)
                cell_size = self.ds_config['geoT'][1]
                #if self.units_are_degrees:
                # scale cellsize to meters,
                # todo: check if input is degress/meters/feet
                cell_size *= 111120 

                m_size = 25000
                utils.echo_msg([m_size, n_den, cell_size])
                if self.size_threshold is None:
                    size_threshold = (m_size * (1/n_den)) / cell_size
                else:
                    size_threshold = self.size_threshold
                
                size_mask = (this_w_arr < self.weight_threshold)
                this_w_arr[size_mask] = 1
                this_w_arr[~size_mask] = 0

                ## group adjacent non-zero cells
                l, n = scipy.ndimage.label(this_w_arr)
                
                ## get the total number of cells in each group
                mn = scipy.ndimage.sum_labels(
                    this_w_arr, labels=l, index=np.arange(1, n+1)
                )
                utils.echo_msg(mn)
                #utils.echo_msg(self.ds_config['nb'])
                #size_threshold = np.nanpercentile(mn, 99)
                
                #size_threshold = self.ds_config['nb']  * .1
                #size_threshold = 10
                utils.echo_msg(size_threshold)
                mn_size_indices = np.where(mn < size_threshold)
                for i in mn_size_indices:
                    i+=1
                    
                this_w_arr[np.isin(l, mn_size_indices)] = 2
                utils.echo_msg(np.count_nonzero(this_w_arr==2))
                for b in range(1, dst_ds.RasterCount+1):
                    this_band = dst_ds.GetRasterBand(b)
                    this_arr = this_band.ReadAsArray()
                    this_arr[this_w_arr==2] = self.ds_config['ndv']
                    this_band.WriteArray(this_arr)
                    this_band.FlushCache()
                    
                if self.weight_is_fn:
                    weight_ds = None
                    
            else:
                utils.echo_msg('failed')
                dst_ds = None
                return(self.src_dem, -1)

        dst_ds = None        
        return(self.dst_dem, 0)

### End
