### hydrolakes.py - hydrolakes
##
## Copyright (c) 2022 CIRES Coastal DEM Team
##
## hydrolakes is part of CUDEM
##
## Permission is hereby granted, free of charge, to any person obtaining a copy 
## of this software and associated documentation files (the "Software"), to deal 
## in the Software without restriction, including without limitation the rights 
## to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies 
## of the Software, and to permit persons to whom the Software is furnished to do so, 
## subject to the following conditions:
##
## The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
##
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
## INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR 
## PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE 
## FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
## ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
##
### Commentary:
##
## hydrolakes:
## https://wp.geog.mcgill.ca/hydrolab/data/
## https://www.hydrosheds.org/products/hydrolakes
## https://data.hydrosheds.org/file/hydrolakes/HydroLAKES_polys_v10_shp.zip
##
## globathy:
## https://springernature.figshare.com/collections/GLOBathy_the_Global_Lakes_Bathymetry_Dataset/5243309
## 
##
### Code:

import os
import json

from osgeo import ogr

from cudem import utils
from cudem import regions
from cudem import datasets

import cudem.fetches.utils as f_utils
import cudem.fetches.FRED as FRED

class HydroLakes(f_utils.FetchModule):
    """hydrolakes"""
    
    def __init__(self, where='1=1', layer=0, want_globathy=False, **kwargs):
        super().__init__(**kwargs)
        self.want_globathy = want_globathy
        self._hydrolakes_prods = 'https://www.hydrosheds.org/products/hydrolakes'
        self._hydrolakes_poly_zip = 'https://data.hydrosheds.org/file/hydrolakes/HydroLAKES_polys_v10_shp.zip'
        self._hydrolakes_gdb_zip = 'https://data.hydrosheds.org/file/hydrolakes/HydroLAKES_polys_v10_gdb.zip'
        self._globathy_url = 'https://springernature.figshare.com/ndownloader/files/28919991'
        self._outdir = os.path.join(os.getcwd(), 'hydrolakes')
        self.name = 'hydrolakes'
        self.where = [where] if len(where) > 0 else []
        
        #self.FRED = FRED.FRED(name=self.name, verbose = self.verbose)
        #self.update_if_not_in_FRED()
        
    def update_if_not_in_FRED(self):
        self.FRED._open_ds()
        self.FRED._attribute_filter(["DataSource = '{}'".format(self.name)])
        if len(self.FRED.layer) == 0:
            self.FRED._close_ds()
            self.update()
        self.FRED._close_ds()
        
    def update(self):
        self.FRED._open_ds()
        v_shp = None
        v_zip = os.path.basename(self._hydrolakes_poly_zip)
        status = f_utils.Fetch(self._hydrolakes_poly_zip, verbose=self.verbose).fetch_file(v_zip)
        v_shps = utils.p_unzip(v_zip, ['shp', 'shx', 'dbf', 'prj'])
        for v in v_shps:
            if '.shp' in v:
                v_shp = v
            
        shp_regions = regions.gdal_ogr_regions(v_shp)
        shp_region = regions.Region()
        for this_region in shp_regions:
            if shp_region.valid_p(check_xy=True):
                shp_region = regions.regions_merge(shp_region, this_region)
            else:
                shp_region = this_region
                
        geom = shp_region.export_as_geom()
        self.FRED._attribute_filter(["ID = '{}'".format('HYDROLAKES')])
        if self.FRED.layer is None or len(self.FRED.layer) == 0:
            self.FRED._add_survey(Name='HydorLakes', ID='HYDROLAKES', Agency='HydroSheds', Date=utils.this_year(),
                                  MetadataLink=self._hydrolakes_prods, MetadataDate=utils.this_year(),
                                  DataLink=self._hydrolakes_poly_zip, IndexLink=self._hydrolakes_poly_zip,
                                  DataType='vector', DataSource='hydrolakes', Info='World-wide lakes', geom=geom)
            
        utils.remove_glob(v_zip, *v_shps)
        self.FRED._close_ds()

    def extract_region(self, in_ogr, out_ogr=None):
        if out_ogr is None:
            out_ogr = '{}_{}.shp'.format('.'.join(in_ogr.split('.')[:-1]), self.region.format('fn'))
            
        out, status = utils.run_cmd(
            'ogr2ogr {} {} -clipsrc {} -nlt POLYGON -skipfailures'.format(out_ogr, in_ogr, self.region.format('ul_lr')),
            verbose=True
        )

        return(out_ogr)

    def extract_shp(self, in_zip, out_ogr=None):
        v_shp = None
        v_shps = utils.p_unzip(in_zip, ['shp', 'shx', 'dbf', 'prj'])
        for v in v_shps:
            if v.split('.')[-1] == 'shp':
                v_shp = v
                break
        if v_shp is not None:
            r_shp = self.extract_region(v_shp, out_ogr)
            return(r_shp)

    ## From GLOBathy
    def apply_calculation(self, shore_distance_arr, max_depth=0, shore_elev=0):
        """
        Apply distance calculation, which is each pixel's distance to shore, multiplied
        by the maximum depth, all divided by the maximum distance to shore. This provides
        a smooth slope from shore to lake max depth.

        shore_distance - Input numpy array containing array of distances to shoreline.
            Must contain positive values for distance away from shore and 0 elsewhere.
        max_depth - Input value with maximum lake depth.
        NoDataVal - value to assign to all non-lake pixels (zero distance to shore).
        """
        
        # Process the distance to shore and max depth into bathymetry
        bathy_arr = (shore_distance_arr * max_depth) / float(shore_distance_arr.max())
        bathy_arr[bathy_arr == 0] = np.nan
        if self.apply_elevations:
            bathy_arr = shore_elev - bathy_arr
        bathy_arr[np.isnan(bathy_arr)] = 0
        return(bathy_arr)

    # def process_bathy(self, polys, depths):
        

    #     globathy_zip = self._fetch_globathy()
    #     glob_csvs = utils.unzip(globathy_zip, self.cache_dir)
    #     #print(glob_csvs)
    #     globathy_csv = os.path.join(self.cache_dir, 'GLOBathy_basic_parameters/GLOBathy_basic_parameters(ALL_LAKES).csv')


    #     ## generate proximity array
    #     prox_srs = osr.SpatialReference()
    #     prox_srs.SetFromUserInput(self.dst_srs)
    #     prox_driver = gdal.GetDriverByName('MEM')
    #     prox_ds = prox_driver.Create('MEM', self.ds_config['nx'], self.ds_config['ny'], 1, self.ds_config['dt'])
    #     prox_ds.SetGeoTransform(self.ds_config['geoT'])
    #     prox_ds.SetProjection(self.ds_config['proj'])
    #     prox_band = prox_ds.GetRasterBand(1)
    #     prox_band.SetNoDataValue(self.ds_config['ndv'])

    #     ## generate mask array
    #     dst_srs = osr.SpatialReference()
    #     dst_srs.SetFromUserInput(self.dst_srs)
    #     driver = gdal.GetDriverByName('MEM')
    #     out_ds = driver.Create('MEM', self.ds_config['nx'], self.ds_config['ny'], 1, self.ds_config['dt'])
    #     out_ds.SetGeoTransform(self.ds_config['geoT'])
    #     out_ds.SetProjection(self.ds_config['proj'])
    #     msk_band = out_ds.GetRasterBand(1)
    #     msk_band.SetNoDataValue(self.ds_config['ndv'])
        
    #     ## process polygons
    #     ll = ogr.Open(lakes_shp)
    #     lll = ll.GetLayer()

    #     lll.SetSpatialFilter(self.p_region.export_as_geom())
    #     lake_features = lll.GetFeatureCount()
    #     _prog = utils.CliProgress('processing {} lakes'.format(lake_features))
    #     bathy_arr = np.zeros( (self.ds_config['ny'], self.ds_config['nx']) )
        
    #     for i, feat in enumerate(lll):
    #         _prog.update_perc((i, lake_features))
    #         lk_id = feat.GetField('Hylak_id')
    #         lk_name = feat.GetField('Lake_name')
    #         lk_elevation = feat.GetField('Elevation')
    #         lk_depth = feat.GetField('Depth_avg')
    #         lk_depth_glb = globd[lk_id]
            
    #         utils.echo_msg('lake: {}'.format(lk_name))
    #         utils.echo_msg('lake elevation: {}'.format(lk_elevation))
    #         utils.echo_msg('lake_depth: {}'.format(lk_depth))
    #         utils.echo_msg('lake_depth (globathy): {}'.format(lk_depth_glb))
    #         #feat_geom = feat.GetGeometryRef()

    #         # Get feature geometry and extent envelope
    #         geom = feat.GetGeometryRef()
            
    #         # Create temporary polygon vector layer from feature geometry so that we can rasterize it (Rasterize needs a layer)
    #         ds = ogr.GetDriverByName('MEMORY').CreateDataSource('')
    #         Layer = ds.CreateLayer('', geom_type=ogr.wkbPolygon, srs=lll.GetSpatialRef())  # Use projection from input vector layer
    #         outfeature = ogr.Feature(Layer.GetLayerDefn())      # Create the output feature
    #         outfeature.SetGeometry(geom)                        # Set geometry of output feature
    #         Layer.SetFeature(outfeature)
            
    #         gdal.RasterizeLayer(out_ds, [1], Layer, burn_values=[1])
    #         msk_band = out_ds.GetRasterBand(1)
    #         msk_band.SetNoDataValue(self.ds_config['ndv'])
    #         #lakes_ds_arr = msk_band.ReadAsArray()
            
    #         proximity_options = ["VALUES=0", "DISTUNITS=PIXEL"]

    #         #stats = msk_band.GetStatistics(0, 1)                                 # Calculate statistics
    #         #print(stats)
    #         gdal.ComputeProximity(msk_band, prox_band, options=proximity_options)
    #         #stats = prox_band.GetStatistics(0, 1)                                 # Calculate statistics
    #         #print(stats)
    #         prox_arr = prox_band.ReadAsArray()
    #         msk_arr = msk_band.ReadAsArray()
            
    #         bathy_arr += self.apply_calculation(prox_band.ReadAsArray(), max_depth=lk_depth_glb, shore_elev=lk_elevation)
            
    #         prox_arr[msk_arr == 1] = 0
    #         prox_ds.GetRasterBand(1).WriteArray(prox_arr)
    #         msk_arr[msk_arr == 1] = 0
    #         out_ds.GetRasterBand(1).WriteArray(msk_arr)
            
    #     _prog.end(0, 'processed {} lakes'.format(lake_features))
    #     bathy_arr[bathy_arr == 0] = self.ndv        
    #     utils.gdal_write(
    #         bathy_arr, '{}.tif'.format(self.name), self.ds_config,
    #     )            
    #     out_ds = lakes_ds_arr = None
    #     return(self)
    
    def run(self):
        self.results.append(
            [self._hydrolakes_poly_zip,
             os.path.join(
                 self._outdir,
                 self._hydrolakes_poly_zip.split('/')[-1]
             ),
             'hydrolakes']
        )

        if self.want_globathy:
            self.results.append(
                [self._globathy_url,
                 os.path.join(
                     self._outdir,
                     'globathy_parameters.zip'
                 ),
                 'globathy']
            )
        return(self)

    def yield_xyz(self, entry):

        ## use globathy.py to get lake depths and return those
        #pass
        
        print(entry)
        if f_utils.Fetch(entry[0], callback=self.callback, verbose=self.verbose, headers=self.headers).fetch_file(entry[1]) == 0:
            if entry[2] == 'hydrolakes1':
                v_shp = None
                v_zip = entry[1]
                v_shps = utils.p_unzip(v_zip, ['shp', 'shx', 'dbf', 'prj'])
                for v in v_shps:
                    if v.split('.')[-1] == 'shp':
                        v_shp = v
                        break

                if v_shp is not None:
                    r_shp = self.extract_region(v_shp)
                                
                utils.remove_glob(v_zip, *v_shps)
                
            elif entry[2] == 'globathy':
                glob_csvs = utils.unzip(entry[1], self._outdir)
                globathy_all = 'GLOBathy_basic_parameters/GLOBathy_basic_parameters(ALL_LAKES).csv'
                globathy_csv = os.path.join(self._outdir, globathy_all)
                
                _ds = datasets.XYZFile(
                    fn=globathy_csv,
                    data_format=168,
                    delim=',',
                    skip=1,
                    xpos=3,
                    ypos=4,
                    zpos=-1,
                    src_srs='epsg:4326',
                    dst_srs=self.dst_srs,
                    src_region=self.region,
                    verbose=self.verbose,
                    remote=True
                )
                
                for xyz in _ds.yield_xyz():
                    yield(xyz)
            else:
                pass
                
### End
