### pointz.py 
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
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
###############################################################################
### Commentary:
##
### Examples:
##
### TODO:
##
### Code:

import os
#from tqdm import tqdm
import numpy as np
from osgeo import ogr
from osgeo import gdal
from cudem import utils
from cudem import gdalfun
from cudem import regions
from cudem import fetches
from cudem import factory


###############################################################################
## PointZ filters
##
## filter points and return the result
###############################################################################
class PointZ:
    """Point Data.

    points is an array of xyz data.
    in this class we manipulate such arrays and return
    the manipulated array
    """

    def __init__(
            self, points=None, region=None, verbose=False,
            xyinc=None, cache_dir='.', **kwargs
    ):
        # if isinstance(points, np.ndarray):
        #     self.points = np.rec.fromrecords(points, names='x, y, z')
        # elif isinstance(points, np.core.records.recarray):
        #     self.points = points
        # elif isinstance(points, pd.DataFrame):
        #     self.points = points
        
        self.points = points
        if self.points is not None and len(self.points) > 0:
            self.region = self.init_region(region)
        else:
            self.region = region

        if xyinc is not None:
            self.xyinc = [utils.str2inc(x) for x in xyinc]
            
        self.verbose = verbose
        self.cache_dir = cache_dir
        self.kwargs = kwargs

        
    def __call__(self):
        if self.verbose:
            utils.echo_msg(f'filtering points using {self}')

        if len(self.points) == 0 or self.points is None:
            return(self.points)

        outliers = self.run()
        # if self.verbose:
        #     utils.echo_msg(f'filtered {len(outliers)} records')
            
        return(self.points)


    def run(self):
        pass

    
    def init_region(self, region):
        """Initialize the data-region AOI
        """
        
        if region is None:
            region = regions.Region().from_list(
                [np.min(self.points['x']), np.max(self.points['x']),
                 np.min(self.points['y']), np.max(self.points['y'])]
            )
            
        return(region)

    
    def fetch_data(self, fetches_module, check_size=True):
        """Fetch data from a fetches module for the data-region
        """
        
        this_fetches = fetches.fetches.FetchesFactory(
            mod=fetches_module,
            src_region=self.region,
            verbose=self.verbose,
            #outdir='./',
            callback=fetches.fetches.fetches_callback
        )._acquire_module()        
        this_fetches.run()
        fr = fetches.fetches.fetch_results(this_fetches, check_size=check_size)
        fr.daemon = True
        fr.start()
        fr.join()
        
        return(fr)

    
    def point_pixels(self, points, x_size = 50, y_size = 50):
        """bin the points to a grid of x_size/y_size and return the 
        associated pixel-z-data at the x/y locations of the points
        """
        
        pa = PointPixels(x_size=x_size, y_size=y_size, ppm=True)
        point_arrays, _, _ = pa(points, mode='mean')
        point_pixels = point_arrays['z'][point_arrays['pixel_y'],
                                         point_arrays['pixel_x']]
        
        return(point_pixels)

    
    def vectorize_points(self):
        """Make a point vector OGR DataSet Object from points
        """

        dst_ogr = 'points_dataset'
        ogr_ds = gdal.GetDriverByName('Memory').Create(
            '', 0, 0, 0, gdal.GDT_Unknown
        )
        layer = ogr_ds.CreateLayer(
            dst_ogr, geom_type=ogr.wkbPoint
        )
        fd = ogr.FieldDefn('index', ogr.OFTInteger)
        layer.CreateField(fd)
        f = ogr.Feature(feature_def=layer.GetLayerDefn())        
        with utils.ccp(
                total=len(self.points),
                desc='vectorizing points dataset', leave=False
        ) as pbar:
            for index, this_row in enumerate(self.points):
                pbar.update()
                f.SetField(0, index)
                g = ogr.CreateGeometryFromWkt(
                    f'POINT ({this_row["x"]} {this_row["y"]})'
                )
                f.SetGeometryDirectly(g)
                layer.CreateFeature(f)
            
        return(ogr_ds)


    def mask_to_raster(self):
        data_mask = gdalfun.ogr2gdal_mask(
            self.mask_fn,
            region=self.region,
            x_inc=self.xyinc[0],
            y_inc=self.xyinc[1],
            #dst_srs=self.dst_srs,
            #invert=True,
            verbose=self.verbose,
            temp_dir=self.cache_dir
        )
        
        return(data_mask)

    
## check if regions overlap before vectorizing points
class PointZVectorMask(PointZ):
    """Filter data using a vector mask

    <mask:mask_fn=path:invert=False>
    """
    
    def __init__(self, mask_fn=None, invert=False, **kwargs):
        super().__init__(**kwargs)
        self.mask_fn = mask_fn
        self.invert = invert
        #self.vectorize_points = vectorize_points
        #if self.xyinc is None:
        #    self.vectorize_points = True
        
        if self.verbose:
            utils.echo_msg(f'masking with {mask_fn}')

            
    def mask_points(self, points, invert=False):
        """mask points by rasterizing the input vector mask
        and querying that with the point data
        """

        if self.verbose:
            utils.echo_msg(
                f'using mask dataset: {self.mask_fn} to xyz'
            )
        ogr_or_gdal = gdalfun.ogr_or_gdal(self.mask_fn)
        if ogr_or_gdal == 1:
            mask_raster = self.mask_to_raster()
        else:
            mask_raster = self.mask_fn
            
        smoothed_depth = gdalfun.gdal_query(
            points, mask_raster, 'g'
        ).flatten()

        outliers = smoothed_depth == 0

        if invert:
            return(points[outliers], outliers)
        else:
            return(points[~outliers], outliers)

        
    def filter_mask(self, points, invert = False):
        os.environ["OGR_OSM_OPTIONS"] = "INTERLEAVED_READING=YES"
        os.environ["OGR_OSM_OPTIONS"] = "OGR_INTERLEAVED_READING=YES"
        outliers = np.zeros(points.shape)

        points_region = regions.Region().from_list(
            [points.x.min(), points.y.min(),
             points.x.max(), points.y.max()]
        )
        
        if self.mask_fn is not None:
            ## vectorize the points
            mask_ds = ogr.Open(self.mask_fn, 0)
            if mask_ds is not None:
                mask_layer = mask_ds.GetLayer()
                mask_geom = gdalfun.ogr_union_geom(
                    mask_layer, verbose=False
                )
                mask_region = regions.Region().from_list(
                    mask_geom.GetEnvelope()
                )

                if regions.regions_intersect_p(points_region, mask_region):                
                    ogr_ds = self.vectorize_points()
                    if ogr_ds is not None:
                        for f in mask_layer:
                            mask_geom = f.geometry()                        
                            points_layer = ogr_ds.GetLayer()
                            points_layer.SetSpatialFilter(mask_geom)
                            for f in points_layer:
                                idx = f.GetField('index')
                                outliers[idx] = 1

                            points_layer.SetSpatialFilter(None)
                        ogr_ds = points_layer = None
                        
                    else:
                        utils.echo_warning_msg(
                            f'could not vectorize {len(self.points)} points for masking'
                        )
                        
                mask_ds = mask_layer = None
                
            else:
                utils.echo_warning_msg(
                    f'could not load mask {self.mask_fn}'
                )

        else:
            utils.echo_warning_msg(
                f'no vector mask was specified {self.mask_fn}'
            )
                
        outliers = outliers == 1
        if self.verbose:
            utils.echo_msg_bold(
                f'found {np.count_nonzero(outliers)} outliers @ {self.mask_fn}'
            )
            
        if invert:
            return(points[outliers], outliers)
        else:
            return(points[~outliers], outliers)

        
    def run(self):
        #if self.vectorize_points:
        self.points, outliers = self.filter_mask(
            self.points, invert=self.invert
        )

        #else:
        #    self.points, outliers = self.mask_points(
        #        self.points, invert=self.invert
        #    )

        #return(self.points)
        return(outliers)

    
class PointZOutlier(PointZ):
    """XYZ outlier filter.

    Find and remove outliers from the points dataset based on
    their residual percentile. 

    <outlierz:percentile=98:multipass=4:invert=False:res=50>
    """
    
    def __init__(
            self, percentile=98, max_percentile=99.9,
            multipass=4, percentage=False, invert=False,
            res=50, max_res=5000, **kwargs
    ):
        super().__init__(**kwargs)
        self.percentile = utils.float_or(percentile, 98)
        self.max_percentile = utils.float_or(max_percentile, 99)
        self.multipass = utils.int_or(multipass, 1)
        self.percentage = percentage
        self.invert = invert
        self.res = utils.int_or(res, 50)
        self.max_res = utils.int_or(max_res, 5000)

        
    def point_residuals(self, points, percentage=False, res=50):
        point_pixels = self.point_pixels(
            points, x_size=res, y_size=res
        )
        if percentage:
            residuals =  np.abs(
                (points['z'] - point_pixels) / point_pixels
            ) * 100
        else:
            residuals = np.abs(points['z'] - point_pixels)

        return(residuals)

    
    def find_outliers(
            self, residuals, percentile=98, percentile_is_threshold=False
    ):
        if percentile_is_threshold:
            outlier_threshold = percentile
        else:
            outlier_threshold = np.percentile(residuals, percentile)

        outliers = residuals > outlier_threshold
        if self.verbose:
            utils.echo_msg_bold(
                f'found {np.count_nonzero(outliers)} outliers @ {percentile}'
            )
        
        return(outliers)

    
    def filter_points(
            self, points, percentile=92, res=50, percentage=False, invert=False
    ):
        residuals = self.point_residuals(
            points, percentage=percentage, res=res
        )
        if residuals is not None:
            outliers = self.find_outliers(
                residuals, percentile=percentile,
                percentile_is_threshold=percentage
            )

            if invert:
                return(points[outliers], outliers)
            else:
                return(points[~outliers], outliers)
        else:
            return(points, None)
        
        
    def run(self):
        percs_it = np.linspace(
            self.percentile, self.max_percentile, self.multipass
        )
        res_it = np.linspace(self.max_res, self.res, self.multipass)
        for mpass in range(0, self.multipass):
            self.points, outliers = self.filter_points(
                self.points, percentile=percs_it[mpass],
                res=res_it[mpass], percentage=self.percentage,
                invert=self.invert
            )
            # if np.count_nonzero(outliers) == 0:
            #     break
            
        #return(self.points)
        return(outliers)

    
## todo: remove the gmrt or other fetched rasters after processing...
class RQOutlierZ(PointZOutlier):
    """xyz outlier filter, using a reference raster

    This will use a reference raster dataset, GMRT by default, to determine
    residual percentages of the input points dataset and remove points which
    have a residual percentage above the given threshold.

    <rq:threshold=5:raster=None>
    """
    
    def __init__(self, threshold=10, raster=None, scaled_percentile=False,
                 resample_raster=True, **kwargs):
        if 'percentile' in kwargs.keys():
            del kwargs['percentile']
        if 'percentage' in kwargs.keys():
            del kwargs['percentage']
        if 'multipass' in kwargs.keys():
            del kwargs['multipass']
            
        super().__init__(
            percentile=threshold, percentage=True, multipass=1,
            **kwargs
        )
        self.threshold = threshold
        self.resample_raster = resample_raster
        self.fetches_modules = ['gmrt', 'CUDEM', 'etopo:datatype=surface']
        self.raster = self.init_raster(raster)            
        self.scaled_percentile = scaled_percentile

        
    def __str__(self):
        return(f'< rq >: {self.raster}')

    
    def mask_gmrt(self, raster):
        #this_fetch = self.fetch_data('gmrt', self.region.copy().buffer(pct=1))
        if os.path.exists(f'{raster}_swath.tif'):
            return(f'{raster}_swath.tif')

        this_fetch = fetches.fetches.FetchesFactory(
            mod='gmrt',
            src_region=self.region,
            verbose=self.verbose,
            callback=fetches.fetches.fetches_callback
        )._acquire_module()        

        swath_masks = []
        if fetches.fetches.Fetch(
                this_fetch._gmrt_swath_poly_url,
                verbose=self.verbose
        ).fetch_file(
                os.path.join(
                    this_fetch._outdir, 'gmrt_swath_polygons.zip'
                )
            ) == 0:
                swath_shps = utils.p_unzip(
                    os.path.join(
                        this_fetch._outdir, 'gmrt_swath_polygons.zip'
                    ),
                    exts=['shp', 'shx', 'prj', 'dbf'],
                    outdir=this_fetch._outdir,
                    verbose=self.verbose
                )

                for v in swath_shps:
                    if '.shp' in v:
                        swath_masks.append(v)
                        break
                    
        if len(swath_masks) > 0:
            for swath_ply in swath_masks:
                gdalfun.gdal_clip(
                    raster, f'{raster}_swath.tif', src_ply=swath_ply, invert=True,
                    verbose=True, cache_dir=this_fetch._outdir
                )

            return(f'{raster}_swath.tif')

        return(None)

    
    def set_raster_fn(self, raster):
        ## raster is not a local file, create a unique name to use
        if (self.region is not None or self.xyinc is not None) and self.resample_raster:
            _raster = utils.append_fn(
                f'rq_raster_{raster}', self.region,
                self.xyinc[0], res=1 if not all(self.xyinc) else None
            )
            _raster = os.path.join(self.cache_dir, f'{_raster}.tif')
            if not os.path.exists(os.path.dirname(_raster)):
                os.makedirs(os.path.dirname(_raster))
            
            if os.path.exists(_raster) and os.path.isfile(_raster):
                return([_raster])

        return(raster)                

        
    def init_raster(self, raster):

        if raster is not None and isinstance(raster, str):
            if os.path.exists(raster) and os.path.isfile(raster):
                return([raster])
            
        elif raster is None:
            if (self.region is not None or self.xyinc is not None) and self.resample_raster:
                _raster = utils.append_fn(
                    f'rq_raster_{raster}', self.region,
                    self.xyinc[0], res=1 if not all(self.xyinc) else None
                )
                _raster = os.path.join(self.cache_dir, f'{_raster}.tif')
                if not os.path.exists(os.path.dirname(_raster)):
                    os.makedirs(os.path.dirname(_raster))

                if os.path.exists(_raster) and os.path.isfile(_raster):
                    return([_raster])
                
            raster = []
            # try gmrt all
            this_fetch = self.fetch_data(
                'gmrt', self.region.copy().buffer(pct=1)
            )
            raster_ = [x[1] for x in this_fetch.results]
            raster.extend([gdalfun.gmt_grd2gdal(x, verbose=False) \
                           if x.split('.')[-1] == 'grd' else x for x in raster_])
            
            # try etopo
            this_fetch = self.fetch_data(
                'etopo:datatype=surface', self.region.copy().buffer(pct=1)
            )
            raster.extend([x[1] for x in this_fetch.results])
            # raster.extend([gdalfun.gmt_grd2gdal(x, verbose=False) \
            #                if x.split('.')[-1] == 'grd' else x for x in raster_])

            # try gmrt swath
            this_fetch = self.fetch_data('gmrt', self.region.copy().buffer(pct=1))
            raster_ = [x[1] for x in this_fetch.results]
            raster_ = [gdalfun.gmt_grd2gdal(x, verbose=False) \
                       if x.split('.')[-1] == 'grd' else x for x in raster_]
            gmrt_swath = self.mask_gmrt(raster_[0])
            if gmrt_swath is not None:
                raster.extend([gmrt_swath])

            # try cudem 1/3
            this_fetch = self.fetch_data(
                'CUDEM:datatype=13:keep_footprints=True', self.region.copy().buffer(pct=1)
            )
            raster.extend([x[1] for x in this_fetch.results])        

            # try cudem 1/9
            this_fetch = self.fetch_data(
                'CUDEM:datatype=19:keep_footprints=True', self.region.copy().buffer(pct=1)
            )
            raster.extend([x[1] for x in this_fetch.results])        
            
            #utils.echo_msg_bold(raster)
            if (self.region is not None or self.xyinc is not None) and self.resample_raster:
                # vrt_options = gdal.BuildVRTOptions(resampleAlg='cubic') # Example option
                # vrt_ds = gdal.BuildVRT('tmp.vrt', raster, options=vrt_options)

                ## todo: use cell count instead of xyinc here
                ##       to make sure the output aligns with the src_dem
                try:
                    raster = [gdalfun.sample_warp(
                        raster, _raster, self.xyinc[0], self.xyinc[1],
                        sample_alg='cubic', src_region=self.region,
                        verbose=self.verbose,
                        co=["COMPRESS=DEFLATE", "TILED=YES"]
                    )[0]]
                except Exception as e:
                    utils.echo_warning_msg(
                        f'failed to process stacked rasters, falling back to GMRT, {e}'
                    )
                    this_fetch = self.fetch_data('gmrt', self.region.copy().buffer(pct=1))
                    raster = [x[1] for x in this_fetch.results]
                    raster = [gdalfun.gmt_grd2gdal(x, verbose=False) \
                              if x.split('.')[-1] == 'grd' else x for x in raster]
                    if self.xyinc is not None and self.resample_raster:
                        raster = [gdalfun.sample_warp(
                            raster[0], _raster, self.xyinc[0], self.xyinc[1],
                            sample_alg='bilinear',
                            verbose=self.verbose,
                            co=["COMPRESS=DEFLATE", "TILED=YES"]
                        )[0]]

        elif any(raster in item for item in self.fetches_modules):
            _raster = [item for item in self.fetches_modules if raster in item][0]
            #elif raster.split(':')[0] in self.fetches_modules:
            this_fetch = self.fetch_data(_raster, self.region)
            raster = [x[1] for x in this_fetch.results]
            raster = [gdalfun.gmt_grd2gdal(x, verbose=False) \
                      if x.split('.')[-1] == 'grd' else x for x in raster]
            if self.xyinc is not None and self.resample_raster:
                raster = [gdalfun.sample_warp(
                    raster, _raster, self.xyinc[0], self.xyinc[1],
                    sample_alg='bilinear', src_region=self.region,
                    verbose=self.verbose,
                    co=["COMPRESS=DEFLATE", "TILED=YES"]
                )[0]]

        else:
            utils.echo_warning_msg(f'could not parse rq raster {raster}')

        #utils.echo_msg_bold(raster)
        return(raster)


    ## todo: allow for multiple rasters
    def point_residuals(self, points, percentage=True, res=50):
        if len(self.raster) == 0:
            return(None)

        #smoothed_depth = []
        #utils.echo_msg(self.raster)
        #for r in self.raster:
        smoothed_depth = gdalfun.gdal_query(
            points, self.raster[0], 'g'
        ).flatten()
        #utils.echo_msg(smoothed_depth)
        #smoothed_depth += smoothed_depth.flatten()

        for x in self.raster:
            x = None
            
        if len(smoothed_depth) == 0:
            return(None)
        
        if percentage:
            if self.scaled_percentile:
                residuals =  np.abs(
                    (points['z'] - smoothed_depth) / (points['z'] + smoothed_depth)
                ) * 100
            else:
                residuals =  np.abs(
                    (points['z'] - smoothed_depth) / smoothed_depth
                ) * 100
        else:
            residuals = np.abs(points['z'] - smoothed_depth)

        #utils.echo_msg(residuals)
        return(residuals)

    
class PointFilterFactory(factory.CUDEMFactory):
    _modules = {
        'outlierz': {
            'name': 'outlierz', 'call': PointZOutlier
        },
        'rq': {
            'name': 'rq', 'call': RQOutlierZ
        },
        'vector_mask': {
            'name': 'vector_mask', 'call': PointZVectorMask
        },
    }

    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        
class PointPixels():
    """Return the point cloud data as an array which coincides with the 
    desired region, x_size and y_size

    incoming data are numpy rec-arrays of x,y,z<w,u> points
    """

    def __init__(
            self,
            src_region=None,
            x_size=None,
            y_size=None,
            verbose=True,
            ppm=False,
            **kwargs
    ):

        self.src_region = src_region
        self.x_size = utils.int_or(x_size, 10)
        self.y_size = utils.int_or(y_size, 10)
        self.verbose = verbose
        self.ppm = ppm

        
    def init_region_from_points(self, points):
        if self.src_region is None:
            self.src_region = regions.Region().from_list(
                [np.min(points['x']), np.max(points['x']),
                 np.min(points['y']), np.max(points['y'])]
            )        

        self.init_gt()

        
    def init_gt(self):
        self.dst_gt = self.src_region.geo_transform_from_count(
            x_count=self.x_size, y_count=self.y_size
        )


    ### todo: allow input of 'out_arrays' for appending
    ### and concluding the stack here instead of _stacks() to file.
    def __call__(
            self, points, weight=1, uncertainty=0,
            mode='mean'
    ):
        out_arrays = {
            'z':None,
            'count':None,
            'weight':None,
            'uncertainty': None,
            'mask':None,
            'x': None,
            'y': None,
            'pixel_x': None,
            'pixel_y': None
        }
        count = 0
        if len(points) == 0:
            return(out_arrays, None, None)
        
        if self.src_region is None:
            self.init_region_from_points(points)
        else:
            self.init_gt()

        weight = utils.float_or(weight, 1)
        uncertainty = utils.float_or(weight, 0)
            
        #######################################################################
        ## convert the points to pixels based on the geotransform
        ## and calculate the local srcwin of the points
        ## pixel_x and pixel_y are the location of the points_* of the srcwin
        #######################################################################            
        pixel_x = np.floor(
            (points['x'] - self.dst_gt[0]) / self.dst_gt[1]
        ).astype(int)
        pixel_y = np.floor(
            (points['y'] - self.dst_gt[3]) / self.dst_gt[5]
        ).astype(int)
        points_x = np.array(points['x'])
        points_y = np.array(points['y'])
        pixel_z = np.array(points['z'])
        try:
            pixel_w = np.array(points['w'])
        except:
            pixel_w = np.ones(pixel_z.shape)

        try:
            pixel_u = np.array(points['u'])
        except:
            pixel_u = np.zeros(pixel_z.shape)

        #######################################################################
        ## remove pixels that will break the srcwin
        #######################################################################
        out_idx = np.nonzero((pixel_x >= self.x_size) \
                             | (pixel_x < 0) \
                             | (pixel_y >= self.y_size) \
                             | (pixel_y < 0))

        if not self.ppm:
            ##commented out for pmm! maybe put back...
            pixel_x = np.delete(pixel_x, out_idx)
            pixel_y = np.delete(pixel_y, out_idx)
            pixel_z = np.delete(pixel_z, out_idx)
            pixel_w = np.delete(pixel_w, out_idx)
            pixel_u = np.delete(pixel_u, out_idx)
            points_x = np.delete(points_x, out_idx)
            points_y = np.delete(points_y, out_idx)
            #if len(pixel_x) == 0 or len(pixel_y) == 0:
            #    continue

        points = None
        pixel_w[np.isnan(pixel_w)] = 1
        pixel_u[np.isnan(pixel_u)] = 0

        #######################################################################
        ## set the srcwin of the incoming points
        #######################################################################
        this_srcwin = (int(min(pixel_x)), int(min(pixel_y)),
                       int(max(pixel_x) - min(pixel_x))+1,
                       int(max(pixel_y) - min(pixel_y))+1)
        count += len(pixel_x)

        #######################################################################
        ## adjust the pixels to the srcwin and stack together
        #######################################################################
        pixel_x = pixel_x - this_srcwin[0]
        pixel_y = pixel_y - this_srcwin[1]
        pixel_xy = np.vstack((pixel_y, pixel_x)).T

        #utils.echo_msg_bold(pixel_x.shape)
        #out_arrays['pixel_x'] = pixel_x
        #out_arrays['pixel_y'] = pixel_y

        #######################################################################
        ## find the non-unique x/y points and mean/min/max
        ## their z values together while calculating the std
        ## for uncertainty
        #######################################################################
        unq, unq_idx, unq_inv, unq_cnt = np.unique(
            pixel_xy, axis=0, return_inverse=True,
            return_index=True, return_counts=True
        )
        cnt_msk = unq_cnt > 1
        cnt_idx, = np.nonzero(cnt_msk)
        idx_msk = np.in1d(unq_inv, cnt_idx)
        idx_idx, = np.nonzero(idx_msk)
        srt_idx = np.argsort(unq_inv[idx_msk])
        dup_idx = np.split(
            idx_idx[srt_idx],
            np.cumsum(unq_cnt[cnt_msk])[:-1]
        )

        #######################################################################
        ## set the output arrays; set to the target grid
        ## arrays will hold weighted-sums
        #######################################################################
        if mode == 'sums':
            ww = pixel_w[unq_idx] * weight
            zz = pixel_z[unq_idx] * ww 
            uu = pixel_u[unq_idx]
            xx = points_x[unq_idx] * ww
            yy = points_y[unq_idx] * ww
        else:
            zz = pixel_z[unq_idx]
            ww = pixel_w[unq_idx]
            uu = pixel_u[unq_idx]
            xx = points_x[unq_idx]
            yy = points_y[unq_idx]
            #u = np.zeros(zz.shape)

        px = pixel_x[unq_idx]
        py = pixel_y[unq_idx]
            
        #######################################################################
        ## min/max/sum the duplicates (values in the same cell)
        ## min and max get the count reset to 1, the count is accumulated
        ## in mean, mixed
        #######################################################################
        if np.any([len(dup) for dup in dup_idx]):
            if mode == 'min':
                dup_stack = [np.min(pixel_z[dup]) for dup in dup_idx]
                dup_stds = np.zeros(dup_stack.shape)
                dup_stack_x = [np.min(pixel_x[dup]) for dup in dup_idx]
                dup_stack_y = [np.min(pixel_y[dup]) for dup in dup_idx]
                dup_counts = [1 for dup in dup_idx]
            elif mode == 'max':
                dup_stack = [np.max(pixel_z[dup]) for dup in dup_idx]
                dup_stds = np.zeros(dup_stack.shape)
                dup_stack_x = [np.max(pixel_x[dup]) for dup in dup_idx]
                dup_stack_y = [np.max(pixel_y[dup]) for dup in dup_idx]
                dup_counts = [1 for dup in dup_idx]
            elif mode == 'mean':
                dup_stack = [np.mean(pixel_z[dup]) for dup in dup_idx]
                dup_stack_x = [np.mean(points_x[dup]) for dup in dup_idx]
                dup_stack_y = [np.mean(points_y[dup]) for dup in dup_idx]
                dup_stds = [np.std(pixel_z[dup]) for dup in dup_idx]
            elif mode == 'sums':
                dup_stack = [np.sum((pixel_z[dup] * (pixel_w[dup] * weight))) \
                             for dup in dup_idx]
                dup_stack_x = [np.sum((points_x[dup] * (pixel_w[dup] * weight))) \
                               for dup in dup_idx]
                dup_stack_y = [np.sum((points_y[dup] * (pixel_w[dup] * weight))) \
                               for dup in dup_idx]
                dup_stds = [np.std(pixel_z[dup]) for dup in dup_idx]
                dup_w = [np.sum((pixel_w[dup] * weight)) for dup in dup_idx]
                ww[cnt_msk] = dup_w
                
            zz[cnt_msk] = dup_stack
            yy[cnt_msk] = dup_stack_y
            xx[cnt_msk] = dup_stack_x
            uu[cnt_msk] = np.sqrt(
                np.power(uu[cnt_msk], 2) + np.power(dup_stds, 2)
            )
                
        ## make the output arrays to yield
        out_x = np.zeros((this_srcwin[3], this_srcwin[2]))
        out_x[unq[:,0], unq[:,1]] = xx
        out_x[out_x == 0] = np.nan
        out_arrays['x'] = out_x

        out_y = np.zeros((this_srcwin[3], this_srcwin[2]))
        out_y[unq[:,0], unq[:,1]] = yy
        out_y[out_y == 0] = np.nan
        out_arrays['y'] = out_y

        out_z = np.zeros((this_srcwin[3], this_srcwin[2]))
        out_z[unq[:,0], unq[:,1]] = zz
        out_z[out_z == 0] = np.nan
        out_arrays['z'] = out_z

        out_arrays['count'] = np.zeros((this_srcwin[3], this_srcwin[2]))
        out_arrays['count'][unq[:,0], unq[:,1]] = unq_cnt

        out_arrays['pixel_x'] = np.zeros((this_srcwin[3], this_srcwin[2]))
        out_arrays['pixel_x'][unq[:,0], unq[:,1]] = px
        
        out_arrays['pixel_y'] = np.zeros((this_srcwin[3], this_srcwin[2]))
        out_arrays['pixel_y'][unq[:,0], unq[:,1]] = py
        
        if mode == 'sums':
            out_arrays['weight'] = np.ones((this_srcwin[3], this_srcwin[2]))
            out_arrays['weight'][unq[:,0], unq[:,1]] = ww #* weight)

        else:
            out_arrays['weight'] = np.ones((this_srcwin[3], this_srcwin[2]))
            out_arrays['weight'][:] = weight if weight is not None else 1
            out_arrays['weight'][unq[:,0], unq[:,1]] *= (ww * unq_cnt)
            #out_arrays['weight'][unq[:,0], unq[:,1]] *= unq_cnt

        out_arrays['uncertainty'] = np.zeros((this_srcwin[3], this_srcwin[2]))
        #out_arrays['uncertainty'][:] = self.uncertainty if self.uncertainty is not None else 0
        out_arrays['uncertainty'][unq[:,0], unq[:,1]] \
            = np.sqrt(uu**2 + (uncertainty if uncertainty is not None else 0)**2)                

        return(out_arrays, this_srcwin, self.dst_gt)


### End
