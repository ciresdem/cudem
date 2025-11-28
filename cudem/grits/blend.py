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
#from cudem.datalists import GDALFile

class Blend(grits.Grits):
    def __init__(self,  weight_threshold=None, buffer_cells=1, gap_fill_cells=None,
                 weight_thresholds=None, buffer_sizes=None, gap_fill_sizes=None,
                 binary_dilation=True, binary_pulse=False,
                 fill_holes=False, aux_dems=None, **kwargs):
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
        #self.init_thresholds()

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

        # aux_params = self._set_params(
        #     data_format=200,
        #     #src_srs=self.fetch_module.src_srs,
        #     cache_dir=self.cache_dir,
        # )
        for aux_fn in self.aux_dems:
            #DatasetFactory(**aux_params)._acquire_module())
            aux_ds = DatasetFactory(mod=aux_fn, src_region=src_region, x_inc=x_inc, y_inc=y_inc)._acquire_module().initialize()
            for arrs, srcwin, gt in aux_ds.array_yield:
                ## supercede
                combined_arr[srcwin[1]:srcwin[1]+srcwin[3],
                             srcwin[0]:srcwin[0]+srcwin[2]] = arrs['z']



        #print(combined_arr)
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

                #src_mask = ~np.isnan(combined_arr[~combined_mask])
                src_mask = ~np.isnan(src_arr)
                
        return(combined_arr, src_arr, combined_mask, src_mask)
        

    def run(self):
        dst_ds = self.copy_src_dem()
        if dst_ds is None:
            return(self.src_dem, -1)

        ## combined_arr is the aux_data, src_arr is the src dem data
        ## combined mask is the combined data + buffer
        ## src_mask is the src data
        combined_arr, src_arr, combined_mask, src_mask = self.init_data()

        ## combined_arr now has a nan buffer between aux and src data
        src_arr[combined_mask] = np.nan
        combined_arr[~combined_mask] = src_arr[~combined_mask]

        ## interpolate combined_arr
        interp_arr = np.full((self.ds_config['ny'], self.ds_config['nx']), np.nan)
        self.method = 'cubic'
        self.chunk_size = None
        self.chunk_step = None
        self.chunk_buffer = 20
        if self.chunk_size is None:
            n_chunk = int(self.ds_config['nx'] * .01)
            n_chunk = 10 if n_chunk < 10 else n_chunk
        else:
            n_chunk = self.chunk_size
            
        if self.chunk_step is None:
            n_step = int(n_chunk/2)
            #n_step = n_chunk
        else:
            n_step = self.chunk_step

        for srcwin in utils.yield_srcwin(
                (combined_arr.shape[0], combined_arr.shape[1]),
                n_chunk=n_chunk,
                msg='generating {} grid'.format(self.method),
                verbose=self.verbose,
                step=n_step
        ):
            chunk_buffer = self.chunk_buffer
            srcwin_buff = utils.buffer_srcwin(
                srcwin, (combined_arr.shape[0], combined_arr.shape[1]), chunk_buffer
            )

            tmp_arr = combined_arr[srcwin_buff[1]:srcwin_buff[1]+srcwin_buff[3],
                                   srcwin_buff[0]:srcwin_buff[0]+srcwin_buff[2]]
            
            tmp_indices = np.nonzero(~np.isnan(tmp_arr))
            if np.count_nonzero(np.isnan(tmp_arr)) == 0:
                y_origin = srcwin[1]-srcwin_buff[1]
                x_origin = srcwin[0]-srcwin_buff[0]
                y_size = y_origin + srcwin[3]
                x_size = x_origin + srcwin[2]
                tmp_arr = tmp_arr[y_origin:y_size,x_origin:x_size]
                #interp_arr[y_origin:y_size,x_origin:x_size] = tmp_arr#[y_origin:y_size,x_origin:x_size]
                interp_arr[srcwin[1]:y_size,srcwin[0]:x_size] = tmp_arr#[y_origin:y_size,x_origin:x_size]
                #continue

            elif len(tmp_indices[0]):
                tmp_values = tmp_arr[tmp_indices]
                xi, yi = np.mgrid[0:srcwin_buff[2],
                                  0:srcwin_buff[3]]

                try:
                    interp_data = interpolate.griddata(
                        np.transpose(tmp_indices), tmp_values,
                        (xi, yi), method=self.method
                    )
                    y_origin = srcwin[1]-srcwin_buff[1]
                    x_origin = srcwin[0]-srcwin_buff[0]
                    y_size = y_origin + srcwin[3]
                    x_size = x_origin + srcwin[2]
                    #combined_arr = combined_arr[y_origin:y_size,x_origin:x_size]
                    interp_data = interp_data[y_origin:y_size,x_origin:x_size]
                    #interp_arr[y_origin:y_size,x_origin:x_size] = interp_data#[y_origin:y_size,x_origin:x_size]
                    interp_arr[srcwin[1]:ysize, srcwin[0]:xsize] = interp_data
                    
                except Exception as e:
                    continue
        
        
        srcwin = (self.buffer_cells, self.buffer_cells, self.ds_config['nx'], self.ds_config['ny'])        
        this_band = dst_ds.GetRasterBand(1)
        #this_arr = this_band.ReadAsArray()
        #this_arr[mask] = self.ds_config['ndv']
        # this_band.WriteArray(combined_arr[srcwin[1]:srcwin[1]+srcwin[3],
        #                                   srcwin[0]:srcwin[0]+srcwin[2]])
        this_band.WriteArray(interp_arr)


        dst_ds = None        
        return(self.dst_dem, 0)

        
    # def init_aux(self):
    #     aux_arr = None
    #     src_region, ds_config = self.init_region(self.src_dem)
    #     src_region.buffer(
    #         x_bv=(ds_config['GeoTransform'][1]*self.buffer_cells),
    #         y_bv=(ds_config['GeoTransform'][3])*self.buffer_cells)
    #     )
        
    #     with gdalfun.gdal_datasource(self.aux_fn) as aux_ds
    #         if aux_ds is not None:
    #             aux_ds_config = gdal_infos(aux_ds)
    #             gt = aux_ds_config['geoT']
    #             srcwin = src_region.srcwin(
    #                 gt, aux_ds.RasterXSize, aux_ds.RasterYSize, node='grid'
    #             )
                
    #             aux_band = aux_ds.GetRasterBand(1)
    #             aux_arr = aux_band.ReadAsArray(*srcwin)

    #     return(aux_arr)


    # def init_src(self):
    #     src_arr = None
    #     aux_region, ds_config = self.init_region(self.aux_fn)
    #     aux_region.buffer(
    #         x_bv=(ds_config['GeoTransform'][1]*self.buffer_cells),
    #         y_bv=(ds_config['GeoTransform'][3])*self.buffer_cells)
    #     )
        
    #     with gdalfun.gdal_datasource(self.src_dem) as src_ds
    #         if src_ds is not None:
    #             src_ds_config = gdal_infos(src_ds)
    #             gt = src_ds_config['geoT']
    #             srcwin = aux_region.srcwin(
    #                 gt, src_ds.RasterXSize, src_ds.RasterYSize, node='grid'
    #             )
                
    #             src_band = src_ds.GetRasterBand(1)
    #             src_arr = src_band.ReadAsArray(*srcwin)

    #     return(src_arr)


    # def init_thresholds(self):
    #     if self.weight_thresholds is not None:
    #         self.weight_thresholds = [float(x) for x in self.weight_thresholds.split('/')]

    #         if self.buffer_sizes is not None:
    #             self.buffer_sizes = [int(x) for x in self.buffer_sizes.split('/')]
    #         else:
    #             self.buffer_sizes = [self.buffer_cells]*len(self.weight_thresholds)

    #         if self.gap_fill_sizes is not None:
    #             self.gap_fill_sizes = [int(x) for x in self.gap_fill_sizes.split('/')]
    #         else:
    #             self.gap_fill_sizes = [self.gap_fill_cells]*len(self.weight_thresholds)
                
    #     else:
    #         self.weight_thresholds = [self.weight_threshold]
    #         self.buffer_sizes = [self.buffer_cells]
    #         self.gap_fill_sizes = [self.gap_fill_cells]

    #     if len(self.buffer_sizes) < len(self.weight_thresholds):
    #         len_diff = len(self.weight_thresholds) - len(self.buffer_sizes)
    #         self.buffer_sizes.extend([1] * len_diff)
            
    #     if len(self.gap_fill_sizes) < len(self.weight_thresholds):
    #         len_diff = len(self.weight_thresholds) - len(self.gap_fill_sizes)
    #         self.gap_fill_sizes.extend([0] * len_diff)

    #     utils.echo_msg([self.weight_thresholds, self.buffer_sizes, self.gap_fill_sizes])

            
    # def init_weight(self, src_ds=None):
    #     weight_band = None
    #     if self.weight_is_fn:
    #         if self.weight_threshold is None:
    #             self.weight_threshold = gdalfun.gdal_percentile(
    #                 self.weight_mask, perc=75
    #             )

    #         weight_ds = gdal.Open(self.weight_mask)
    #         weight_band = weight_ds.GetRasterBand(1)

    #     elif self.weight_is_band and src_ds is not None:
    #         if self.weight_threshold is None:
    #             self.weight_threshold = gdalfun.gdal_percentile(
    #                 self.src_dem, perc=50, band=self.weight_mask
    #             )

    #         weight_band = src_ds.GetRasterBand(self.weight_mask)

    #     return(weight_band)


    # def binary_closed_dilation(self, arr, iterations=1, closing_iterations=1):
    #     closed_and_dilated_arr = arr.copy()
    #     struct_ = scipy.ndimage.generate_binary_structure(2, 2)
    #     for i in range(closing_iterations):
    #         closed_and_dilated_arr = scipy.ndimage.binary_dilation(
    #             closed_and_dilated_arr, iterations=1, structure=struct_
    #         )
    #         closed_and_dilated_arr = scipy.ndimage.binary_erosion(
    #             closed_and_dilated_arr, iterations=1, border_value=1, structure=struct_
    #         )

    #     closed_and_dilated_arr = scipy.ndimage.binary_dilation(
    #         closed_and_dilated_arr, iterations=iterations, structure=struct_
    #     )
            
    #     return(closed_and_dilated_arr)


    # def binary_reversion(self, arr, iterations, closing_iterations):
    #     reversion_arr = arr.copy()
    #     closing_iterations += iterations
            
    #     struct_ = scipy.ndimage.generate_binary_structure(2, 2)
    #     reversion_arr = scipy.ndimage.binary_dilation(
    #         reversion_arr, iterations=closing_iterations, structure=struct_
    #     )
    #     erosion_iterations = max(closing_iterations-iterations, 0)
    #     if erosion_iterations > 0:
    #         reversion_arr = scipy.ndimage.binary_erosion(
    #             reversion_arr, iterations=erosion_iterations, border_value=1, structure=struct_
    #         )

    #     return(reversion_arr)


    # def buffer_to_nd(self):
    #     with gdalfun.gdal_datasource(self.src_dem) as src_ds:
    #         if src_ds is not None:
    #             self.init_ds(src_ds)                
    #             weight_band = self.init_weight(src_ds)
    #             w_arr = weight_band.ReadAsArray()#*srcwin)
    #             w_arr[w_arr == self.ds_config['ndv']] = np.nan
    #             weights = np.unique(w_arr)[::-1]
    #             this_w_arr = w_arr.copy()

    #             with tqdm(
    #                     total=len(self.weight_thresholds),
    #                     desc=f'buffering around weights over {self.weight_thresholds}',
    #                     leave=self.verbose
    #             ) as w_pbar:
    #                 for n, weight_threshold in enumerate(self.weight_thresholds):
    #                     this_w_arr[this_w_arr < weight_threshold] = np.nan
    #                     if self.binary_dilation:
    #                         if self.binary_pulse:
    #                             weight_mask = this_w_arr >= weight_threshold
    #                             expanded_w_arr = self.binary_closed_dilation(
    #                                 weight_mask, iterations=self.buffer_sizes[n],
    #                                 closing_iterations=self.gap_fill_sizes[n]
    #                             )

    #                         else:
    #                             weight_mask = this_w_arr >= weight_threshold
    #                             expanded_w_arr = self.binary_reversion(
    #                                 weight_mask, self.buffer_sizes[n],
    #                                 self.gap_fill_sizes[n]
    #                             )
    #                     else:
    #                         ## mrl use ndimage.binary_dilation rather than the slower expand_for below
    #                         shiftx = self.buffer_sizes[n] + self.gap_fill_sizes[n]#self.buffer_sizes[n] * self.revert_multiplier
    #                         shifty = self.buffer_sizes[n] + self.gap_fill_sizes[n]#self.buffer_sizes[n] * self.revert_multiplier
    #                         utils.echo_msg([shiftx, shifty])
    #                         expanded_w_arr = utils.expand_for(
    #                             this_w_arr >= weight_threshold,
    #                             shiftx=int(shiftx), shifty=int(shifty)
    #                         )

    #                         shiftx = max(0, shiftx - self.buffer_sizes[n])
    #                         shifty = max(0, shifty - self.buffer_sizes[n])
    #                         if shiftx > 0 and shifty > 0:
    #                             utils.echo_msg([shiftx, shifty])
    #                             contracted_w_arr = np.invert(utils.expand_for(
    #                                 np.invert(expanded_w_arr),
    #                                 shiftx=int(shiftx), shifty=int(shifty),
    #                             ))
    #                             expanded_w_arr = contracted_w_arr.copy()

    #                     if self.fill_holes:
    #                         expanded_w_arr = scipy.ndimage.binary_fill_holes(expanded_w_arr)

    #                     mask = (w_arr < weight_threshold) & expanded_w_arr
    #                     this_w_arr = w_arr.copy()
    #                     this_w_arr[mask] = np.nan

    #                     w_pbar.update()

    #                 return(this_w_arr)
                
    #             if self.weight_is_fn:
    #                 weight_ds = None

    #         else:
    #             utils.echo_msg('could not load {self.src_dem}')
    #             dst_ds = None
    #             return(self.src_dem, -1)

        
    # def run(self):
    #     if self.weight_mask is None:
    #         return(self.src_dem, -1)
        
    #     dst_ds = self.copy_src_dem()
    #     if dst_ds is None:
    #         return(self.src_dem, -1)

    #     this_masked_w_arr = self.buffer_to_nd(self)
    #     mask = np.isnan(this_masked_w_arr)

    #     ## interpolate masked data
    #     ## weighted average of interpolation and original data
        
    #     ## mask out any values in other bands of the input/output raster
    #     for b in range(1, dst_ds.RasterCount+1):
    #         this_band = dst_ds.GetRasterBand(b)
    #         this_arr = this_band.ReadAsArray()#*srcwin)
    #         this_arr[mask] = self.ds_config['ndv']
    #         this_band.WriteArray(this_arr)#, srcwin[0], srcwin[1])
    #         this_band.FlushCache()
    #         this_band = this_arr = None
        
    #     dst_ds = None        
    #     return(self.dst_dem, 0)


### End
