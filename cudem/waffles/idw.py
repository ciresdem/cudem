### idw.py
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
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
###############################################################################
### Commentary:
##
### Code:


import numpy as np
from scipy import spatial

from cudem import utils
from cudem import gdalfun
from cudem.waffles.waffles import Waffle

class Invdisttree():
    """ inverse-distance-weighted interpolation using KDTree:
    @Denis via https://stackoverflow.com/questions/3104781/inverse-distance-weighted-idw-interpolation-with-python
    https://creativecommons.org/licenses/by-nc-sa/3.0/

invdisttree = Invdisttree( X, z )  -- data points, values
interpol = invdisttree( q, nnear=3, eps=0, p=1, weights=None, stat=0 )
    interpolates z from the 3 points nearest each query point q;
    For example, interpol[ a query point q ]
    finds the 3 data points nearest q, at distances d1 d2 d3
    and returns the IDW average of the values z1 z2 z3
        (z1/d1 + z2/d2 + z3/d3)
        / (1/d1 + 1/d2 + 1/d3)
        = .55 z1 + .27 z2 + .18 z3  for distances 1 2 3

    q may be one point, or a batch of points.
    eps: approximate nearest, dist <= (1 + eps) * true nearest
    p: use 1 / distance**p
    weights: optional multipliers for 1 / distance**p, of the same shape as q
    stat: accumulate wsum, wn for average weights

How many nearest neighbors should one take ?
a) start with 8 11 14 .. 28 in 2d 3d 4d .. 10d; see Wendel's formula
b) make 3 runs with nnear= e.g. 6 8 10, and look at the results --
    |interpol 6 - interpol 8| etc., or |f - interpol*| if you have f(q).
    I find that runtimes don't increase much at all with nnear -- ymmv.

p=1, p=2 ?
    p=2 weights nearer points more, farther points less.
    In 2d, the circles around query points have areas ~ distance**2,
    so p=2 is inverse-area weighting. For example,
        (z1/area1 + z2/area2 + z3/area3)
        / (1/area1 + 1/area2 + 1/area3)
        = .74 z1 + .18 z2 + .08 z3  for distances 1 2 3
    Similarly, in 3d, p=3 is inverse-volume weighting.

Scaling:
    if different X coordinates measure different things, Euclidean distance
    can be way off.  For example, if X0 is in the range 0 to 1
    but X1 0 to 1000, the X1 distances will swamp X0;
    rescale the data, i.e. make X0.std() ~= X1.std() .

A nice property of IDW is that it's scale-free around query points:
if I have values z1 z2 z3 from 3 points at distances d1 d2 d3,
the IDW average
    (z1/d1 + z2/d2 + z3/d3)
    / (1/d1 + 1/d2 + 1/d3)
is the same for distances 1 2 3, or 10 20 30 -- only the ratios matter.
In contrast, the commonly-used Gaussian kernel exp( - (distance/h)**2 )
is exceedingly sensitive to distance and to h.

quite heavy on memory when large grid-size...

    """
# anykernel( dj / av dj ) is also scale-free
# error analysis, |f(x) - idw(x)| ? todo: regular grid, nnear ndim+1, 2*ndim

    def __init__( self, X, z, leafsize=10, stat=0 ):
        assert len(X) == len(z), "len(X) %d != len(z) %d" % (len(X), len(z))
        self.tree = spatial.cKDTree( X, leafsize=leafsize )  # build the tree
        self.z = z
        self.stat = stat
        self.wn = 0
        self.wsum = None;

    def __call__(
            self, q,
            nnear=6,
            eps=0,
            p=1,
            dub=np.inf,
            weights=None,
            uncertainties=None
    ):
        # nnear nearest neighbours of each query point --
        q = np.asarray(q)
        qdim = q.ndim
        if qdim == 1:
            q = np.array([q])
            
        if self.wsum is None:
            self.wsum = np.zeros(nnear)

        eps = utils.float_or(eps, .1)
        self.distances, self.ix \
            = self.tree.query(q, k=nnear, eps=eps, distance_upper_bound=dub)
        interpol = np.zeros((len(self.distances),) + np.shape(self.z[0]))
        jinterpol = 0
        for dist, ix in zip( self.distances, self.ix ):
            if np.any(np.isinf(dist)):
                wz = np.nan
            elif nnear == 1:
                wz = self.z[ix]
            elif dist[0] < 1e-10:
                wz = self.z[ix[0]]
            else:  # weight z s by 1/dist --
                w = 1 / dist**p
                if weights is not None:
                    w *= weights[ix]  # >= 0

                # if uncertainties is not None:
                #     w *= (1/uncertainties[ix])
                    
                w /= np.sum(w)
                wz = np.dot( w, self.z[ix] )
                if self.stat:
                    self.wn += 1
                    self.wsum += w
                    
            interpol[jinterpol] = wz
            jinterpol += 1
        return(interpol if qdim > 1  else interpol[0])

    
class WafflesIDW(Waffle):
    """INVERSE DISTANCE WEIGHTED DEM
    
    Generate a DEM using an Inverse Distance Weighted algorithm.
    If weights are used, will generate a UIDW DEM, using weight 
    values as inverse uncertainty, as described here: 
    https://ir.library.oregonstate.edu/concern/graduate_projects/79407x932
    and here: 
    https://stackoverflow.com/questions/3104781/inverse-distance-weighted-idw-interpolation-with-python

    -----------
    Parameters:
    
    power=[val] - weight**power
    min_points=[val] - minimum neighbor points for IDW
    radius=[val] - search radius (in cells), only fill data cells 
                   within radius from data
    chunk_size=[val] - size of chunks in pixels

    < IDW:min_points=8:radius=inf:power=1:chunk_size=None >
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
            if self.min_points:
                utils.echo_msg(
                    (f'generating IDW grid @ {self.ycount}/{self.xcount} '
                     f'looking for at least {self.min_points} neighbors '
                     f'within {self.radius} pixels')
                )
            else:
                utils.echo_msg(
                    f'generating IDW grid @ {self.ycount}/{self.xcount}'
                )
            i=0

        if self.chunk_size is None:
            n_chunk = int(self.ds_config['nx'] * .1)
            n_chunk = 10 if n_chunk < 10 else n_chunk
        else:
            n_chunk = self.chunk_size
            
        if self.chunk_step is None:
            n_step = int(n_chunk/2)
        else:
            n_step = self.chunk_step

        with gdalfun.gdal_datasource(self.stack) as stack_ds:
            #stack_ds = gdal.Open(self.stack)
            points_band = stack_ds.GetRasterBand(1)
            points_no_data = points_band.GetNoDataValue()
            weights_band = stack_ds.GetRasterBand(3)
            weights_no_data = weights_band.GetNoDataValue()
            uncertainty_band = stack_ds.GetRasterBand(4)
            uncertainty_no_data = uncertainty_band.GetNoDataValue()

            #try:
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
            interp_ds.SetProjection(stack_ds.GetProjection())
            interp_ds.SetGeoTransform(stack_ds.GetGeoTransform())
            interp_band = interp_ds.GetRasterBand(1)
            interp_band.SetNoDataValue(np.nan)
            #except:
            #    return(self)

            points_array = points_band.ReadAsArray()
            #points_array[points_array == points_no_data] = np.nan
            #point_indices = np.nonzero(~np.isnan(points_array))
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

            #stack_ds = None
            
        invdisttree = Invdisttree(
            np.transpose(point_indices),
            point_values,
            leafsize=10,
            stat=1
        )
        for srcwin in utils.yield_srcwin(
                (self.ycount, self.xcount),
                n_chunk=n_chunk,
                msg='Generating IDW DEM',
                end_msg='Generated IDW DEM',
                verbose=self.verbose
        ):
            # if np.count_nonzero(np.isnan(points_array)) == 0:
            #     #if not np.all(~np.nan(points_array)):
            #     # y_origin = srcwin[1]-srcwin_buff[1]
            #     # x_origin = srcwin[0]-srcwin_buff[0]
            #     # y_size = y_origin + srcwin[3]
            #     # x_size = x_origin + srcwin[2]
            #     # points_array = points_array[y_origin:y_size,x_origin:x_size]
            #     interp_band.WriteArray(points_array, srcwin[0], srcwin[1])

            if len(point_indices[0]):
                xi, yi = np.mgrid[srcwin[0]:srcwin[0]+srcwin[2],
                                  srcwin[1]:srcwin[1]+srcwin[3]]
                interp_data = invdisttree(
                    np.vstack((yi.flatten(), xi.flatten())).T,
                    nnear=self.min_points,
                    eps=.1,
                    p=self.power,
                    dub=self.radius,
                    weights=weight_values,
                    uncertainties=uncertainty_values
                )
                interp_data = np.reshape(interp_data, (srcwin[2], srcwin[3]))
                interp_band.WriteArray(interp_data.T, srcwin[0], srcwin[1])
                
        interp_ds = point_values = weight_values = None
        
        return(self)    

    
### End
