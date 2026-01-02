### idw.py
##
## Copyright (c) 2022 - 2026 Regents of the University of Colorado
##
## idw.py is part of CUDEM
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
##
### Commentary:
##
## Inverse Distance Weighted (IDW) interpolation module.
##
### Code:

import numpy as np
from scipy import spatial

from cudem import utils
from cudem import gdalfun
from cudem.waffles.waffles import Waffle

class Invdisttree:
    """Inverse-distance-weighted interpolation using KDTree.
    
    Adapted from: https://stackoverflow.com/questions/3104781/inverse-distance-weighted-idw-interpolation-with-python
    """

    def __init__(self, X, z, leafsize=10, stat=0):
        assert len(X) == len(z), f"len(X) {len(X)} != len(z) {len(z)}"
        self.tree = spatial.cKDTree(X, leafsize=leafsize)
        self.z = z
        self.stat = stat
        self.wn = 0
        self.wsum = None

        
    def __call__(
            self, q,
            nnear=6,
            eps=0,
            p=1,
            dub=np.inf,
            weights=None,
            uncertainties=None
    ):
        # q may be one point, or a batch of points.
        q = np.asarray(q)
        qdim = q.ndim
        if qdim == 1:
            q = np.array([q])
            
        if self.wsum is None:
            self.wsum = np.zeros(nnear)

        eps = utils.float_or(eps, .1)
        self.distances, self.ix = self.tree.query(
            q, k=nnear, eps=eps, distance_upper_bound=dub
        )
        
        interpol = np.zeros((len(self.distances),) + np.shape(self.z[0]))
        jinterpol = 0
        
        for dist, ix in zip(self.distances, self.ix):
            if np.any(np.isinf(dist)):
                wz = np.nan
            elif nnear == 1:
                wz = self.z[ix]
            elif dist[0] < 1e-10:
                wz = self.z[ix[0]]
            else: 
                w = 1 / dist**p
                if weights is not None:
                    w *= weights[ix]

                # Inverse Uncertainty Weighting support (future)
                # if uncertainties is not None:
                #     w *= (1/uncertainties[ix])
                    
                w /= np.sum(w)
                wz = np.dot(w, self.z[ix])
                
                if self.stat:
                    self.wn += 1
                    self.wsum += w
                    
            interpol[jinterpol] = wz
            jinterpol += 1
            
        return interpol if qdim > 1 else interpol[0]

    
class WafflesIDW(Waffle):
    """INVERSE DISTANCE WEIGHTED DEM
    
    Generate a DEM using an Inverse Distance Weighted algorithm.
    Interpolation occurs in pixel-space (indices), so anisotropic 
    behavior may occur if pixels are not square.

    Parameters:
    -----------
    power (float): weight**power (default: 1)
    min_points (int): minimum neighbor points for IDW (default: 8)
    radius (float): search radius in cells (default: inf)
    chunk_size (int): processing chunk size in pixels
    """
    
    def __init__(
            self,
            power=1,
            min_points=8,
            radius=None,
            chunk_size=None,
            **kwargs
    ):
        super().__init__(**kwargs)
        self.power = utils.float_or(power)
        self.min_points = utils.int_or(min_points)
        self.radius = np.inf if radius is None else utils.str2inc(radius) 
        self.chunk_size = chunk_size
        self.chunk_step = None

        
    def run(self):
        if self.verbose:
            msg = f'generating IDW grid @ {self.ycount}/{self.xcount}'
            if self.min_points:
                msg += f' looking for at least {self.min_points} neighbors within {self.radius} pixels'
            utils.echo_msg(msg)

        if self.chunk_size is None:
            n_chunk = int(self.ds_config['nx'] * .1)
            n_chunk = 10 if n_chunk < 10 else n_chunk
        else:
            n_chunk = self.chunk_size
            
        ## Initialize Output Dataset
        with gdalfun.gdal_datasource(self.stack) as stack_ds:
            points_band = stack_ds.GetRasterBand(1)
            points_no_data = points_band.GetNoDataValue()
            weights_band = stack_ds.GetRasterBand(3)
            uncertainty_band = stack_ds.GetRasterBand(4)

            interp_ds = stack_ds.GetDriver().Create(
                self.fn,
                stack_ds.RasterXSize,
                stack_ds.RasterYSize,
                bands=1,
                eType=points_band.DataType,
                options=["BLOCKXSIZE=256", "BLOCKYSIZE=256", "TILED=YES", 
                         "COMPRESS=LZW", "BIGTIFF=YES"]
            )
            interp_ds.SetProjection(stack_ds.GetProjection())
            interp_ds.SetGeoTransform(stack_ds.GetGeoTransform())
            interp_band = interp_ds.GetRasterBand(1)
            interp_band.SetNoDataValue(np.nan)

            ## NOTE: This reads the entire raster into memory.
            ## Large stacks may cause MemoryError.
            points_array = points_band.ReadAsArray()
            point_indices = np.nonzero(points_array != points_no_data)
            point_values = points_array[point_indices]
            points_array = None

            if self.want_weight:
                weights_array = weights_band.ReadAsArray()
                weight_values = weights_array[point_indices]
                weights_array = None
            else:
                weight_values = None

            if self.want_uncertainty:
                uncertainty_array = uncertainty_band.ReadAsArray()
                uncertainty_values = uncertainty_array[point_indices]
                uncertainty_array = None
            else:
                uncertainty_values = None

        ## Build KDTree
        invdisttree = Invdisttree(
            np.transpose(point_indices),
            point_values,
            leafsize=10,
            stat=1
        )
        
        ## Process Chunks
        for srcwin in utils.yield_srcwin(
                (self.ycount, self.xcount),
                n_chunk=n_chunk,
                msg='Generating IDW DEM',
                end_msg='Generated IDW DEM',
                verbose=self.verbose
        ):
            if len(point_indices[0]):
                ## Create grid coordinates for the current chunk
                ## srcwin: (x_off, y_off, x_size, y_size)
                ## indices: [rows (y), cols (x)]
                xi, yi = np.mgrid[srcwin[0]:srcwin[0]+srcwin[2],
                                  srcwin[1]:srcwin[1]+srcwin[3]]
                
                # Interpolate
                interp_data = invdisttree(
                    np.vstack((yi.flatten(), xi.flatten())).T,
                    nnear=self.min_points,
                    eps=.1,
                    p=self.power,
                    dub=self.radius,
                    weights=weight_values,
                    uncertainties=uncertainty_values
                )
                
                ## Reshape and write (Transpose needed because mgrid order vs GDAL write order)
                interp_data = np.reshape(interp_data, (srcwin[2], srcwin[3]))
                interp_band.WriteArray(interp_data.T, srcwin[0], srcwin[1])
                
        interp_ds = point_values = weight_values = None
        
        return self

### End
