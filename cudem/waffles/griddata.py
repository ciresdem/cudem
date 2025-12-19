### griddata.py
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
##
## griddata.py is part of CUDEM
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
###############################################################################
### Commentary:
##
### Code:


import threading

from osgeo import gdal

import numpy as np
from scipy import interpolate

from cudem import utils
from cudem.waffles.waffles import Waffle

def interpolate_array(points_arr,
                      ndv=-9999,
                      chunk_size=None,
                      chunk_step=None,
                      chunk_buffer=0,
                      method='cubic',
                      verbose=True):

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
    
    return(interp_array)


def write_array_queue(wq, q, m):
    while True:
        wq_args = wq.get()
        m.interp_band.WriteArray(wq_args[0], wq_args[1][0], wq_args[1][1])
        wq.task_done()

        
def scipy_queue(q, wq, m, p):
    while True:
        this_srcwin = q.get()
        p.update()
        try:
            interp_array = m.grid_srcwin(this_srcwin)
        except Exception as e:
            utils.echo_msg(e)
            utils.echo_warning_msg(
                f'failed to grid srcwin {this_srcwin}, placing back into queue'
            )
            q.put(this_srcwin)
            q.task_done()
            continue

        wq.put([interp_array, this_srcwin])
        q.task_done()

        
class grid_scipy(threading.Thread):
    """Scipy interpolate.griddata gridding 
    (linear, cubic, nearest)

    https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.griddata.html
    
    set 'num_threads' to int over 1 to generate in multiple threads...
    """
    
    def __init__(self, mod, n_threads=3):
        threading.Thread.__init__(self)
        self.mod = mod
        self.scipy_q = queue.Queue()
        self.grid_q = queue.Queue()
        self.n_threads = n_threads

        
    def run(self):
        for this_srcwin in utils.yield_srcwin(
                n_size=(self.mod.ycount, self.mod.xcount),
                n_chunk=self.mod.chunk_size,
                step=self.mod.chunk_size,
                verbose=True
        ):
            self.scipy_q.put(this_srcwin)

        self.pbar = utils.ccp(
            desc=f'gridding data to {self.mod}',
            total=self.scipy_q.qsize()
        )
        for _ in range(1):
            tg = threading.Thread(
                target=write_array_queue,
                args=(self.grid_q, self.scipy_q, self.mod)
            )
            tg.daemon = True
            tg.start()

        for _ in range(self.n_threads):
            t = threading.Thread(
                target=scipy_queue,
                args=(self.scipy_q, self.grid_q, self.mod, self.pbar)
            )
            t.daemon = True
            t.start()

        self.grid_q.join()
        self.scipy_q.join()
        self.pbar.close()

        
class WafflesSciPy(Waffle):
    """Generate DEM using Scipy gridding interpolation
    
    Generate a DEM using Scipy's gridding interpolation
    Optional gridding methods are 'linear', 'cubic' and 'nearest'
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.griddata.html

    -----------
    Parameters:
    
    method=[linear/cubic/nearest] - interpolation method to use
    chunk_size=[val] - size of chunks in pixels
    chunk_buffer=[val] - size of the chunk buffer in pixels
    chunk_step=[val] - size of the chunk step in pixels
    num_threads=[val] - number of threads to use in interpolation, 
                        use None to process as a single thread

    < scipy:method=<method>:chunk_size=None:chunk_buffer=40 >
    """
    
    def __init__(
            self,
            method='linear',
            chunk_size=None,
            chunk_buffer=20,
            chunk_step=None,
            num_threads=None,
            **kwargs
    ):
        """generate a `scipy` dem"""
        
        super().__init__(**kwargs)
        self.methods = ['linear', 'cubic', 'nearest']
        self.method = method
        self.chunk_size = utils.int_or(chunk_size)
        self.chunk_step = utils.int_or(chunk_step)
        self.chunk_buffer = utils.int_or(chunk_buffer)
        self.num_threads = utils.int_or(num_threads)

        
    ## this 'run' command runs the scipy module in a multiple threads
    def _run(self):
        if self.num_threads is None:
            return(self._run())
            
        self.open()
        try:
            gs = grid_scipy(self, n_threads=self.num_threads)
            gs.daemon = True
    
            gs.start()
            gs.join()
        except (KeyboardInterrupt, SystemExit):
            utils.echo_error_msg(
                'user breakage...please wait while fetches exits.'
            )
            stop_threads = True
            while not gs.scipy_q.empty():
                try:
                    gs.scipy_q.get(False)
                except Empty:
                    continue
                        
                gs.scipy_q.task_done()
            
        self.close()
        
        return(self)

    
    def open(self):
        if self.method not in self.methods:
            utils.echo_error_msg(
                (f'{self.method} is not a valid interpolation method, '
                 f'options are {self.methods}')
            )
            return(self)
        
        if self.chunk_size is None:
            self.chunk_size = int(self.ds_config['nx'] * .1)
            self.chunk_size = 10 \
                if self.chunk_size < 10 \
                   else self.chunk_size
            
        if self.chunk_step is None:
            self.chunk_step = int(self.chunk_size/2)

        self.stack_ds = gdal.Open(self.stack)
        self.points_band = self.stack_ds.GetRasterBand(1)
        self.points_no_data = self.points_band.GetNoDataValue()        
        self.interp_ds = self.stack_ds.GetDriver().Create(
            self.fn,
            self.stack_ds.RasterXSize,
            self.stack_ds.RasterYSize,
            bands=1,
            eType=self.points_band.DataType,
            options=["BLOCKXSIZE=256",
                     "BLOCKYSIZE=256",
                     "TILED=YES",
                     "COMPRESS=LZW",
                     "BIGTIFF=YES"]
        )
        if self.interp_ds is not None:
            self.interp_ds.SetProjection(self.stack_ds.GetProjection())
            self.interp_ds.SetGeoTransform(self.stack_ds.GetGeoTransform())
            self.interp_band = self.interp_ds.GetRasterBand(1)
            self.interp_band.SetNoDataValue(np.nan)
        else:
            utils.echo_error_msg(
                f'could not create {self.fn}...'
            )
            return(self)
        
        if self.verbose:
            utils.echo_msg(
                f'buffering srcwin by {self.chunk_buffer} pixels'
            )

        self.points_array = self.points_band.ReadAsArray()
        self.stack_ds = None

        
    def close(self):
        self.interp_ds = self.stack_ds = None            

        
    def grid_srcwin(self, srcwin):
        srcwin_buff = utils.buffer_srcwin(
            srcwin, (self.ycount, self.xcount), self.chunk_buffer
        )
        points_array = self.points_array[
            srcwin_buff[1]:srcwin_buff[1]+srcwin_buff[3],
            srcwin_buff[0]:srcwin_buff[0]+srcwin_buff[2]
        ]
        points_array[points_array == self.points_no_data] = np.nan
        point_indices = np.nonzero(~np.isnan(points_array))
        if np.count_nonzero(np.isnan(points_array)) == 0:
            y_origin = srcwin[1]-srcwin_buff[1]
            x_origin = srcwin[0]-srcwin_buff[0]
            y_size = y_origin + srcwin[3]
            x_size = x_origin + srcwin[2]
            points_array = points_array[y_origin:y_size,x_origin:x_size]
            return(points_array)

        elif len(point_indices[0]):
            point_values = points_array[point_indices]
            xi, yi = np.mgrid[0:srcwin_buff[2],
                              0:srcwin_buff[3]]
            #try:
            interp_data = interpolate.griddata(
                np.transpose(point_indices), point_values,
                (xi, yi), method=self.method
            )
            # while np.any(interp_data[np.isnan(interp_data)]):
            #     utils.echo_msg('nodata in {}'.format(srcwin))
            #     point_indices = np.nonzero(~np.isnan(interp_data))
            #     point_values = interp_data[point_indices]
            #     interp_data = interpolate.griddata(
            #         np.transpose(point_indices), point_values,
            #         (xi, yi), method=self.method
            #     )

            y_origin = srcwin[1]-srcwin_buff[1]
            x_origin = srcwin[0]-srcwin_buff[0]
            y_size = y_origin + srcwin[3]
            x_size = x_origin + srcwin[2]
            interp_data = interp_data[y_origin:y_size,x_origin:x_size]
            return(interp_data)
            #except Exception as e:
            #    return(points_array)
                
        return(None)

    
    ## this 'run' command runs the scipy module in a single thread
    def run(self):
        if self.method not in self.methods:
            utils.echo_error_msg(
                (f'{self.method} is not a valid interpolation method, '
                 f'options are {self.methods}')
            )
            return(self)
        
        if self.chunk_size is None:
            n_chunk = int(self.ds_config['nx'] * .1)
            n_chunk = 10 if n_chunk < 10 else n_chunk
        else:
            n_chunk = self.chunk_size
            
        if self.chunk_step is None:
            n_step = int(n_chunk/2)
            #n_step = n_chunk
        else:
            n_step = self.chunk_step

        stack_ds = gdal.Open(self.stack)
        points_band = stack_ds.GetRasterBand(1)
        points_no_data = points_band.GetNoDataValue()
        interp_ds = stack_ds.GetDriver().Create(
            self.fn,
            stack_ds.RasterXSize,
            stack_ds.RasterYSize,
            bands=1,
            eType=points_band.DataType,
            options=["BLOCKXSIZE=256",
                     "BLOCKYSIZE=256",
                     "TILED=YES",
                     "COMPRESS=LZW",
                     "BIGTIFF=YES"]
        )
        if interp_ds is not None:
            interp_ds.SetProjection(stack_ds.GetProjection())
            interp_ds.SetGeoTransform(stack_ds.GetGeoTransform())
            interp_band = interp_ds.GetRasterBand(1)
            interp_band.SetNoDataValue(np.nan)
        else:
            utils.echo_error_msg(
                f'could not create {self.fn}...'
            )
            return(self)
        
        if self.verbose:
            utils.echo_msg(
                f'buffering srcwin by {self.chunk_buffer} pixels'
            )
       
        for srcwin in utils.yield_srcwin(
                (self.ycount, self.xcount),
                n_chunk=n_chunk,
                msg='generating {} grid'.format(self.method),
                verbose=self.verbose,
                step=n_step
        ):
            chunk_buffer = self.chunk_buffer
            srcwin_buff = utils.buffer_srcwin(
                srcwin, (self.ycount, self.xcount), chunk_buffer
            )
            points_array = points_band.ReadAsArray(*srcwin_buff)
            points_array[points_array == points_no_data] = np.nan
            point_indices = np.nonzero(~np.isnan(points_array))
            if np.count_nonzero(np.isnan(points_array)) == 0:
                y_origin = srcwin[1]-srcwin_buff[1]
                x_origin = srcwin[0]-srcwin_buff[0]
                y_size = y_origin + srcwin[3]
                x_size = x_origin + srcwin[2]
                points_array = points_array[y_origin:y_size,x_origin:x_size]
                interp_band.WriteArray(points_array, srcwin[0], srcwin[1])

            elif len(point_indices[0]):
                point_values = points_array[point_indices]
                xi, yi = np.mgrid[0:srcwin_buff[2],
                                  0:srcwin_buff[3]]

                try:
                    interp_data = interpolate.griddata(
                        np.transpose(point_indices), point_values,
                        (xi, yi), method=self.method
                    )
                    y_origin = srcwin[1]-srcwin_buff[1]
                    x_origin = srcwin[0]-srcwin_buff[0]
                    y_size = y_origin + srcwin[3]
                    x_size = x_origin + srcwin[2]
                    interp_data = interp_data[y_origin:y_size,x_origin:x_size]
                    interp_band.WriteArray(interp_data, srcwin[0], srcwin[1])
                except Exception as e:
                    continue
                
        interp_ds = stack_ds = point_values = weight_values = None
        
        return(self)

    
class WafflesLinear(WafflesSciPy):
    """LINEAR (triangulated) DEM
    
    -----------
    Parameters:
    
    chunk_size (int): size of chunks in pixels
    chunk_step (iint):  size of chunks in pixels
    chunk_buffer (int):  size of the chunk buffer in pixels
    """
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.method = 'linear'

        
class WafflesCubic(WafflesSciPy):
    """CUBIC (triangulated) DEM
    
    -----------
    Parameters:
    
    chunk_size=[val] - size of chunks in pixels
    chunk_buffer=[val] - size of the chunk buffer in pixels
    """
        
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.method = 'cubic'

        
class WafflesNearest(WafflesSciPy):
    """NEAREST neighbor DEM
    
    -----------
    Parameters:
    
    chunk_size=[val] - size of chunks in pixels
    chunk_buffer=[val] - size of the chunk buffer in pixels
    """
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.method = 'nearest'


### End
