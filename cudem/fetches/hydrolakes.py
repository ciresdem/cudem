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
import csv

from osgeo import ogr
from osgeo import osr
from osgeo import gdal

import numpy as np

from cudem import utils
from cudem import regions
from cudem import datasets
from cudem import demfun

import cudem.fetches.utils as f_utils
import cudem.fetches.FRED as FRED
import cudem.fetches.copernicus
        
class HydroLakes(f_utils.FetchModule):
    """hydrolakes"""
    
    def __init__(self, where='1=1', want_globathy=False, **kwargs):
        super().__init__(**kwargs)
        self.want_globathy = want_globathy
        self._hydrolakes_prods = 'https://www.hydrosheds.org/products/hydrolakes'
        self._hydrolakes_poly_zip = 'https://data.hydrosheds.org/file/hydrolakes/HydroLAKES_polys_v10_shp.zip'
        self._hydrolakes_gdb_zip = 'https://data.hydrosheds.org/file/hydrolakes/HydroLAKES_polys_v10_gdb.zip'
        self._globathy_url = 'https://springernature.figshare.com/ndownloader/files/28919991'
        self._outdir = os.path.join(os.getcwd(), 'hydrolakes')
        self.name = 'hydrolakes'
        self.where = [where] if len(where) > 0 else []
        
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

    def _load_bathy(self):
        """create a nodata grid"""
        
        xcount, ycount, gt = self.region.geo_transform(x_inc=self.x_inc, y_inc=self.y_inc)
        self.ds_config = demfun.set_infos(
            xcount,
            ycount,
            xcount * ycount,
            gt,
            utils.sr_wkt(self.dst_srs),
            gdal.GDT_Float32,
            -9999,
            'GTiff'
        )        
        self.bathy_arr = np.zeros( (ycount, xcount) )

    def _fetch_lakes(self):
        """fetch hydrolakes polygons"""

        this_lakes = cudem.fetches.hydrolakes.HydroLakes(
            src_region=self.region, verbose=self.verbose
        )
        this_lakes._outdir = self._outdir
        this_lakes.run()

        fr = cudem.fetches.utils.fetch_results(this_lakes, want_proc=False)
        fr.daemon = True
        fr.start()
        fr.join()

        lakes_shp = None
        lakes_zip = this_lakes.results[0][1]
        lakes_shps = utils.unzip(lakes_zip, self._outdir)
        for i in lakes_shps:
            if i.split('.')[-1] == 'shp':
                lakes_shp = i

        return(lakes_shp)
        
    def _fetch_globathy(self):
        """fetch globathy csv data and process into dict"""
        
        _globathy_url = 'https://springernature.figshare.com/ndownloader/files/28919991'
        globathy_zip = os.path.join(self._outdir, 'globathy_parameters.zip')
        f_utils.Fetch(_globathy_url, verbose=self.verbose).fetch_file(globathy_zip)
        globathy_csvs = utils.unzip(globathy_zip, self._outdir)        
        globathy_csv = os.path.join(self._outdir, 'GLOBathy_basic_parameters/GLOBathy_basic_parameters(ALL_LAKES).csv')
        with open(globathy_csv, mode='r') as globc:
            reader = csv.reader(globc)
            next(reader)
            globd = {int(row[0]):float(row[-1]) for row in reader}
        
        return(globd)

    def _fetch_copernicus(self):
        """copernicus"""

        this_cop = cudem.fetches.copernicus.CopernicusDEM(
            src_region=self.region, verbose=self.verbose, datatype='1'
        )
        this_cop._outdir = self._outdir
        this_cop.run()

        fr = f_utils.fetch_results(this_cop, want_proc=False)
        fr.daemon = True
        fr.start()
        fr.join()

        if self.dst_srs is not None:
            dst_srs = osr.SpatialReference()
            dst_srs.SetFromUserInput(self.dst_srs)
        else:
            dst_srs = None
            
        cop_ds = demfun.generate_mem_ds(self.ds_config, name='copernicus')        
        for i, cop_tif in enumerate(this_cop.results):
            gdal.Warp(cop_ds, cop_tif[1], dstSRS=dst_srs, resampleAlg='bilinear')
            
        return(cop_ds)
    
    def generate_mem_ogr(self, geom, srs):
        """Create temporary polygon vector layer from feature geometry 
        so that we can rasterize it (Rasterize needs a layer)
        """
        
        ds = ogr.GetDriverByName('MEMORY').CreateDataSource('')
        Layer = ds.CreateLayer('', geom_type=ogr.wkbPolygon, srs=srs)
        outfeature = ogr.Feature(Layer.GetLayerDefn())
        outfeature.SetGeometry(geom)
        Layer.SetFeature(outfeature)

        return(ds)

    ##
    ## From GLOBathy
    def apply_calculation(self, shore_distance_arr, max_depth=0, shore_arr=None, shore_elev=0):
        """
        Apply distance calculation, which is each pixel's distance to shore, multiplied
        by the maximum depth, all divided by the maximum distance to shore. This provides
        a smooth slope from shore to lake max depth.

        shore_distance - Input numpy array containing array of distances to shoreline.
            Must contain positive values for distance away from shore and 0 elsewhere.
        max_depth - Input value with maximum lake depth.
        NoDataVal - value to assign to all non-lake pixels (zero distance to shore).
        """

        max_dist = float(shore_distance_arr.max())
        max_dist = .1 if max_dist <= 0 else max_dist
        bathy_arr = (shore_distance_arr * max_depth) / max_dist
        bathy_arr[bathy_arr == 0] = np.nan
        if shore_arr is None or np.all(np.isnan(bathy_arr)) or shore_arr[~np.isnan(bathy_arr)].max() == 0:
            bathy_arr = shore_elev - bathy_arr
        else:
            bathy_arr = shore_arr - bathy_arr
            
        bathy_arr[np.isnan(bathy_arr)] = 0
        return(bathy_arr)
        
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

    def generate_elevations(self, entry):
        if f_utils.Fetch(entry[0], callback=self.callback, verbose=self.verbose, headers=self.headers).fetch_file(entry[1]) == 0:
            if entry[2] == 'hydrolakes':
                self._load_bathy()
                lakes_shp = None
                lakes_zip = entry[1]
                lakes_shps = utils.unzip(lakes_zip, self._outdir)
                for i in lakes_shps:
                    if i.split('.')[-1] == 'shp':
                        lakes_shp = i
                        break

                globd = self._fetch_globathy()
                cop_ds = self._fetch_copernicus()
                cop_band = cop_ds.GetRasterBand(1)
                prox_ds = demfun.generate_mem_ds(self.ds_config, name='prox')
                msk_ds = demfun.generate_mem_ds(self.ds_config, name='msk')

                lk_ds = ogr.Open(lakes_shp)
                lk_layer = lk_ds.GetLayer()
                lk_layer.SetSpatialFilter(self.region.export_as_geom())
                lk_features = lk_layer.GetFeatureCount()
                lk_prog = utils.CliProgress('processing {} lakes'.format(lk_features))

                for i, feat in enumerate(lk_layer):
                    lk_id = feat.GetField('Hylak_id')
                    lk_prog.update_perc((i, lk_features), msg='processing {} lakes: {}'.format(lk_features, lk_id))
                    lk_elevation = feat.GetField('Elevation')
                    lk_depth = feat.GetField('Depth_avg')
                    lk_depth_glb = globd[lk_id]

                    tmp_ds = self.generate_mem_ogr(feat.GetGeometryRef(), lk_layer.GetSpatialRef())
                    tmp_layer = tmp_ds.GetLayer()            
                    gdal.RasterizeLayer(msk_ds, [1], tmp_layer, burn_values=[1])            
                    msk_band = msk_ds.GetRasterBand(1)
                    msk_band.SetNoDataValue(self.ds_config['ndv'])

                    prox_band = prox_ds.GetRasterBand(1)
                    proximity_options = ["VALUES=0", "DISTUNITS=PIXEL"]
                    gdal.ComputeProximity(msk_band, prox_band, options=proximity_options)

                    prox_arr = prox_band.ReadAsArray()
                    msk_arr = msk_band.ReadAsArray()
                    
                    self.bathy_arr += self.apply_calculation(
                        prox_band.ReadAsArray(),
                        max_depth=lk_depth_glb,
                        shore_arr=cop_band.ReadAsArray(),
                        shore_elev=lk_elevation
                    )

                    prox_arr[msk_arr == 1] = 0
                    prox_ds.GetRasterBand(1).WriteArray(prox_arr)
                    msk_arr[msk_arr == 1] = 0
                    msk_ds.GetRasterBand(1).WriteArray(msk_arr)
                    tmp_ds = None

                lk_prog.end(0, 'processed {} lakes'.format(lk_features))
                self.bathy_arr[self.bathy_arr == 0] = self.ds_config['ndv']
                utils.gdal_write(
                    self.bathy_arr, '{}.tif'.format(self.name), self.ds_config,
                )            
                prox_ds = msk_ds = None

                return('{}.tif'.format(self.name))

    def yield_xyz(self, entry):
        lk_elev = self.generate_elevations(entry)

        _ds = datasets.RasterFile(
            fn=lk_elev,
            data_format=200,
            dst_srs=self.dst_srs,
            src_srs='epsg:4326',
            src_region=self.region,
            x_inc=self.x_inc,
            y_inc=self.y_inc,
            weight=self.weight,
            verbose=self.verbose
        )
        for xyz in _ds.yield_xyz():
            yield(xyz)

    def yield_array(self, entry):
        lk_elev = self.generate_elevations(entry)

        _ds = datasets.RasterFile(
            fn=lk_elev,
            data_format=200,
            dst_srs=self.dst_srs,
            src_srs='epsg:4326',
            src_region=self.region,
            x_inc=self.x_inc,
            y_inc=self.y_inc,
            weight=self.weight,
            verbose=self.verbose
        )
        for arr in _ds.yield_array():
            yield(arr)
            
                    
### End
