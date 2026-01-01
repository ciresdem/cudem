### gdalgrid.py
##
## Copyright (c) 2010 - 2026 Regents of the University of Colorado
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
### Commentary:
##
## Wrapper for GDAL Grid (gdal_grid) functionality.
##
### Code:

from osgeo import gdal
from osgeo import ogr
from cudem import utils
from cudem import gdalfun
from cudem.waffles.waffles import Waffle

class WafflesGDALGrid(Waffle):
    """Waffles GDAL_GRID module.

    Wraps the gdal_grid functionality for interpolation.
    """
    
    def __init__(self, **kwargs):
        """Run gdal grid using alg_str.

        Args: 
          alg_str (str): the gdal_grid algorithm string (e.g. 'linear:radius=-1')
        """
        
        super().__init__(**kwargs)
        if not hasattr(self, 'alg_str'):
            self.alg_str = 'linear:radius=-1'
        self.mod = self.alg_str.split(':')[0]

        
    def _vectorize_stack(self):
        """Make a point vector OGR DataSet Object from src_xyz.
        
        Converts the processing stack into an OGR Memory layer 
        required by gdal.Grid.
        """

        dst_ogr = '{}'.format(self.name)
        drv = gdal.GetDriverByName('Memory')
        ogr_ds = drv.Create('', 0, 0, 0, gdal.GDT_Unknown)
        
        layer = ogr_ds.CreateLayer(dst_ogr, geom_type=ogr.wkbPoint25D)
        
        # Define fields
        for x in ['long', 'lat', 'elev', 'weight']:
            fd = ogr.FieldDefn(x, ogr.OFTReal)
            fd.SetWidth(12)
            fd.SetPrecision(8)
            layer.CreateField(fd)
            
        feat_def = layer.GetLayerDefn()
        
        with utils.ccp(desc='vectorizing stack', leave=self.verbose) as pbar:
            for this_xyz in self.stack_ds.yield_xyz():
                pbar.update()
                
                f = ogr.Feature(feat_def)
                f.SetField(0, this_xyz.x)
                f.SetField(1, this_xyz.y)
                f.SetField(2, float(this_xyz.z))
                if self.want_weight:
                    f.SetField(3, this_xyz.w)

                geom = ogr.Geometry(ogr.wkbPoint25D)
                geom.SetPoint(0, this_xyz.x, this_xyz.y, float(this_xyz.z))
                f.SetGeometryDirectly(geom)
                
                layer.CreateFeature(f)
                f = None
            
        return ogr_ds

    
    def run(self):
        """Execute the GDAL Grid algorithm."""
        
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
        # Callback wrapper for GDAL
        def _prog_update(pct, msg, user_data):
            _prog.update_pct(pct * 100)
            return 1

        ogr_ds = self._vectorize_stack()
        if ogr_ds.GetLayer().GetFeatureCount() == 0:
            utils.echo_error_msg('no input data')
            return self
            
        gd_opts = gdal.GridOptions(
            outputType=gdal.GDT_Float32,
            noData=self.ndv,
            format='GTiff',
            width=self.xcount,
            height=self.ycount,
            algorithm=self.alg_str,
            callback=_prog_update if self.verbose else None,
            outputBounds=[
                self.p_region.xmin, self.p_region.ymax,
                self.p_region.xmax, self.p_region.ymin
            ]
        )
        
        out_fn = '{}.tif'.format(self.name)
        gdal.Grid(out_fn, ogr_ds, options=gd_opts)
        
        gdalfun.gdal_set_ndv(out_fn, ndv=self.ndv, convert_array=False)
        
        ogr_ds = None
        return self

    
class GDALLinear(WafflesGDALGrid):
    """LINEAR DEM via gdal_grid.
    
    Parameters:
    -----------
    radius (float) : search radius (default: 4 * xinc)
    nodata (float) : nodata value
    """
    
    def __init__(self, radius=None, nodata=-9999, **kwargs):
        super().__init__(**kwargs)        
        radius = self.xinc * 4 if radius is None else utils.str2inc(radius)
        self.alg_str = 'linear:radius={}:nodata={}'.format(radius, nodata)

        
class GDALInvDst(WafflesGDALGrid):
    """INVERSE DISTANCE DEM via gdal_grid.
    
    Parameters:
    -----------
    radius1, radius2 (float) : search radii
    power (float) : weight**power
    smoothing (float) : smoothing parameter
    angle (float) : search ellipse angle
    max_points, min_points (int) : point counts per bucket
    nodata (float) : nodata value
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
    """Moving AVERAGE DEM via gdal_grid.
    
    Parameters:
    -----------
    radius1, radius2 (float) : search radii
    angle (float) : search ellipse angle
    min_points (int) : minimum points per bucket
    nodata (float) : nodata value
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
    """NEAREST DEM via gdal_grid.
    
    Parameters:
    -----------
    radius1, radius2 (float) : search radii
    angle (float) : search ellipse angle
    nodata (float) : nodata value
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
