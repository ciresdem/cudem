### gdalgrid.py
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
##
## gdalgrid.py is part of CUDEM
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

from osgeo import gdal
from osgeo import ogr
from cudem import utils
from cudem import gdalfun
from cudem.waffles.waffles import Waffle

class WafflesGDALGrid(Waffle):
    """Waffles GDAL_GRID module.

    see gdal_grid for more info and gridding algorithms
    """
    
    def __init__(self, **kwargs):
        """run gdal grid using alg_str

        parse the data through xyzfun.xyz_block to get weighted 
        mean before building the GDAL dataset to pass into gdal_grid
        
        Args: 
          alg_str (str): the gdal_grid algorithm string
        """
        
        super().__init__(**kwargs)
        self.alg_str = 'linear:radius=-1'
        self.mod = self.alg_str.split(':')[0]

        
    def _vectorize_stack(self):
        """Make a point vector OGR DataSet Object from src_xyz

        for use in gdal gridding functions.
        """

        dst_ogr = '{}'.format(self.name)
        ogr_ds = gdal.GetDriverByName('Memory').Create(
            '', 0, 0, 0, gdal.GDT_Unknown
        )
        layer = ogr_ds.CreateLayer(
            dst_ogr, geom_type=ogr.wkbPoint25D
        )
        for x in ['long', 'lat', 'elev', 'weight']:
            fd = ogr.FieldDefn(x, ogr.OFTReal)
            fd.SetWidth(12)
            fd.SetPrecision(8)
            layer.CreateField(fd)
            
        f = ogr.Feature(feature_def=layer.GetLayerDefn())        
        with utils.ccp(desc='vectorizing stack', leave=self.verbose) as pbar:
            for this_xyz in self.stack_ds.yield_xyz():
                pbar.update()
                f.SetField(0, this_xyz.x)
                f.SetField(1, this_xyz.y)
                f.SetField(2, float(this_xyz.z))
                if self.want_weight:
                    f.SetField(3, this_xyz.w)

                wkt = this_xyz.export_as_wkt(include_z=True)
                g = ogr.CreateGeometryFromWkt(wkt)
                f.SetGeometryDirectly(g)
                layer.CreateFeature(f)
            
        return(ogr_ds)

    
    def run(self):
        if self.verbose:
            utils.echo_msg(
                'running GDAL GRID {} algorithm @ {} and {}/{}...'.format(
                    self.alg_str.split(':')[0],
                    self.p_region.format('fn'),
                    self.xcount,
                    self.ycount
                )
            )
        _prog = utils.ccp(
            desc=f'running GDAL GRID {self.alg_str} algorithm',
            leave=self.verbose
        )
        _prog_update = lambda x, y, z: _prog.update()
        ogr_ds = self._vectorize_stack()
        if ogr_ds.GetLayer().GetFeatureCount() == 0:
            utils.echo_error_msg('no input data')
            
        gd_opts = gdal.GridOptions(
            outputType = gdal.GDT_Float32,
            noData = self.ndv,
            format = 'GTiff',
            width = self.xcount,
            height = self.ycount,
            algorithm = self.alg_str,
            callback = _prog_update if self.verbose else None,
            outputBounds = [
                self.p_region.xmin, self.p_region.ymax,
                self.p_region.xmax, self.p_region.ymin
            ]
        )
        gdal.Grid(
            '{}.tif'.format(self.name),
            ogr_ds,
            options=gd_opts
        )
        gdalfun.gdal_set_ndv(
            '{}.tif'.format(
                self.name, ndv=self.ndv, convert_array=False
            )
        )
        ogr_ds = None
        return(self)

    
class GDALLinear(WafflesGDALGrid):
    """LINEAR DEM via gdal_grid
    Generate a DEM using GDAL's gdal_grid command.
    see gdal_grid --help for more info

    -----------
    Parameters:
    
    radius=[val] - search radius
    
    < gdal-linear:radius=-1 >
    """
    
    def __init__(self, radius=None, nodata=-9999, **kwargs):
        super().__init__(**kwargs)        
        radius = self.xinc * 4 if radius is None else utils.str2inc(radius)
        self.alg_str = 'linear:radius={}:nodata={}'.format(radius, nodata)

        
class GDALInvDst(WafflesGDALGrid):
    """INVERSE DISTANCE DEM via gdal_grid
    
    Generate a DEM using GDAL's gdal_grid command.
    see gdal_grid --help for more info

    -----------
    Parameters:
    
    radius1=[val] - search radius 1
    radius2=[val] - search radius 2
    power=[val] - weight**power
    min_points=[val] - minimum points per IDW bucket (use to fill entire DEM)

    < gdal-invdst:power=2.0:smoothing=0.0:radius1=0:radius2=0:angle=0:max_points=0:min_points=0:nodata=0 >
    """
    
    def __init__(
            self,
            power=2.0,
            smoothing=0.0,
            radius1=None,
            radius2=None,
            angle=0.0,
            max_points=0,
            min_points=0,
            nodata=-9999,
            **kwargs
    ):
        super().__init__(**kwargs)
        radius1 = self.xinc * 2 if radius1 is None else utils.str2inc(radius1)
        radius2 = self.yinc * 2 if radius2 is None else utils.str2inc(radius2)
        self.alg_str = (f'invdist:power={power}'
                        f':smoothing={smoothing}'
                        f':radius1={radius1}'
                        f':radius2={radius2}'
                        f':angle={angle}'
                        f':max_points={max_points}'
                        f':min_points={min_points}'
                        f':nodata={nodata}')

        
class GDALMovingAverage(WafflesGDALGrid):
    """Moving AVERAGE DEM via gdal_grid
    
    Generate a DEM using GDAL's gdal_grid command.
    see gdal_grid --help for more info
    
    -----------
    Parameters:
    
    radius1=[val] - search radius 1
    radius2=[val] - search radius 2
    min_points=[val] - minimum points per bucket 
                       (use to fill entire DEM)

    < gdal-average:radius1=0:radius2=0:angle=0:min_points=0:nodata=0 >
    """
    
    def __init__(
            self,
            radius1=None,
            radius2=None,
            angle=0.0,
            min_points=0,
            nodata=-9999,
            **kwargs
    ):
        super().__init__(**kwargs)
        radius1 = self.xinc * 2 if radius1 is None else utils.str2inc(radius1)
        radius2 = self.yinc * 2 if radius2 is None else utils.str2inc(radius2)
        self.alg_str = (f'average:radius1={radius1}'
                        f':radius2={radius2}'
                        f':angle={angle}'
                        f':min_points={min_points}'
                        f':nodata={nodata}')

        
class GDALNearest(WafflesGDALGrid):
    """NEAREST DEM via gdal_grid
    
    Generate a DEM using GDAL's gdal_grid command.
    see gdal_grid --help for more info

    -----------
    Parameters:
    
    radius1=[val] - search radius 1
    radius2=[val] - search radius 2
    angle=[val] - angle
    nodata=[val] - nodata

    < gdal-nearest:radius1=0:radius2=0:angle=0:nodata=0 >
    """
    
    def __init__(
            self,
            radius1=None,
            radius2=None,
            angle=0.0,
            nodata=-9999,
            **kwargs
    ):
        super().__init__(**kwargs)
        radius1 = self.xinc * 2 if radius1 is None else utils.str2inc(radius1)
        radius2 = self.yinc * 2 if radius2 is None else utils.str2inc(radius2)
        self.alg_str = (f'nearest:radius1={radius1}'
                        f':radius2={radius2}'
                        f':angle={angle}'
                        f':nodata={nodata}')


### End
