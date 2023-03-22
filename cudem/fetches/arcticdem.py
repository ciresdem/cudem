### arcticdem.py - Arctic DEM
##
## Copyright (c) 2022, 2023 CIRES Coastal DEM Team
##
## arctic.py is part of CUDEM
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
##  objectid (Integer64)
##  name (String)
##  tile (String)
##  nd_value (Real)
##  resolution (Real)
##  creationda (Date)
##  raster (String)
##  fileurl (String)
##  spec_type (String)
##  qual (Real)
##  reg_src (String)
##  num_gcps (Integer)
##  meanresz (Real)
##  active (Integer)
##  qc (Integer)
##  rel_ver (String)
##  num_comp (Integer)
##
##
## TODO: fix projection issues...
##
### Code:

import os
import json

from osgeo import ogr

from cudem import utils
from cudem import regions
from cudem import datasets
from cudem import demfun

import cudem.fetches.utils as f_utils
import cudem.fetches.FRED as FRED

class ArcticDEM(f_utils.FetchModule):
    """Arctic DEM
ArcticDEM is an NGA-NSF public-private initiative to automatically produce a high-resolution, 
high quality, digital surface model (DSM) of the Arctic using optical stereo imagery, 
high-performance computing, and open source photogrammetry software.

https://www.pgc.umn.edu/data/arcticdem/

< arcticdem >"""
    
    def __init__(self, where='1=1', layer=0, **kwargs):
        super().__init__(name='arcticdem', **kwargs)
        self._arctic_dem_index_url = 'https://data.pgc.umn.edu/elev/dem/setsm/ArcticDEM/indexes/ArcticDEM_Tile_Index_Rel7.zip'
        self.where = [where] if len(where) > 0 else []
        self.arctic_region = self.region.copy()
        self.arctic_region.warp('epsg:3413')
        
        self.FRED = FRED.FRED(name=self.name, verbose = self.verbose)
        self.update_if_not_in_FRED()

    def update_if_not_in_FRED(self):
        self.FRED._open_ds()
        self.FRED._attribute_filter(["DataSource = '{}'".format(self.name)])
        if len(self.FRED.layer) == 0:
            self.FRED._close_ds()
            self.update()
        self.FRED._close_ds()
        
    def update(self):
        self.FRED._open_ds()        
        v_zip = os.path.basename(self._arctic_dem_index_url)
        status = f_utils.Fetch(self._arctic_dem_index_url, verbose=self.verbose).fetch_file(v_zip)
        v_shps = utils.p_unzip(v_zip, ['shp', 'shx', 'dbf', 'prj'])

        v_shp = None
        for v in v_shps:
            if '.shp' in v:
                v_shp = v
                break

        utils.run_cmd('ogr2ogr arctic_tmp.shp {} -t_srs epsg:4326'.format(v_shp), verbose=True)
        utils.remove_glob(v_zip, *v_shps)
        v_shp = 'arctic_tmp.shp'
        v_shps = ['arctic_tmp.shp','arctic_tmp.dbf','arctic_tmp.shx','arctic_tmp.prj']
        shp_regions = regions.gdal_ogr_regions(v_shp)
        shp_region = regions.Region()
        for this_region in shp_regions:
            #this_region.src_srs = 'epsg:3413'
            #this_region.warp('epsg:4326')
            if shp_region.valid_p(check_xy=True):
                shp_region = regions.regions_merge(shp_region, this_region)
            else: shp_region = this_region
        geom = shp_region.export_as_geom()
        
        self.FRED._attribute_filter(["ID = '{}'".format('ARCTICDEM-1')])
        if self.FRED.layer is None or len(self.FRED.layer) == 0:
            self.FRED._add_survey(Name = 'ArcticDEM', ID = 'ARCTICDEM-1', Agency = 'UMN', Date = utils.this_year(),
                                  MetadataLink = 'https://data.pgc.umn.edu/elev/dem/setsm/ArcticDEM/', MetadataDate = utils.this_year(),
                                  DataLink = self._arctic_dem_index_url, IndexLink = self._arctic_dem_index_url,
                                  DataType = 'raster', DataSource = 'arcticdem', Info = 'Arctic Only', geom = geom)
        utils.remove_glob(*v_shps)
        self.FRED._close_ds()

    def run(self):
        #for surv in FRED._filter_FRED(self):

        v_zip = os.path.join(self._outdir, os.path.basename(self._arctic_dem_index_url))
        status = f_utils.Fetch(self._arctic_dem_index_url, verbose=self.verbose).fetch_file(v_zip)
        v_shps = utils.p_unzip(v_zip, ['shp', 'shx', 'dbf', 'prj'], outdir=self._outdir)
        v_shp = None
        for v in v_shps:
            if v.split('.')[-1] == 'shp':
                v_shp = v
                break

        #v_out_shp = os.path.join(self._outdir, 'arctic_tmp.shp')
        #utils.run_cmd('ogr2ogr {} {} -t_srs epsg:4326'.format(v_out_shp, v_shp), verbose=True)
        #utils.remove_glob(v_zip, *v_shps)
        #v_shps = ['arctic_tmp.shp','arctic_tmp.dbf','arctic_tmp.shx','arctic_tmp.prj']
        try:
            v_ds = ogr.Open(v_shp)
        except:
            v_ds = None
            status = -1

        if v_ds is not None:
            layer = v_ds.GetLayer()
            _boundsGeom = self.arctic_region.export_as_geom()
            layer.SetSpatialFilter(_boundsGeom)
            fcount = layer.GetFeatureCount()
            utils.echo_msg('filtered {} arcticdem features'.format(fcount))
            for f in range(0, fcount):
                feature = layer[f]
                data_link = feature.GetField('fileurl')
                self.results.append([data_link, os.path.join(self._outdir, data_link.split('/')[-1]), 'raster'])

            v_ds = None
            
        return(self)

    def yield_xyz(self, entry):
        """yield the xyz data from the arcticdem fetch module"""

        status = f_utils.Fetch(entry[0], callback=self.callback, verbose=self.verbose).fetch_file(entry[1])
        if status == 0:
            src_dems = utils.p_untar(entry[1], ['reg_dem.tif'])
            for src_dem in src_dems:
                #src_srs = demfun.get_srs(src_dem)
                _ds = datasets.RasterFile(
                    fn=src_dem,
                    data_format=200,
                    src_srs='epsg:3413',
                    dst_srs=self.dst_srs,
                    src_region=self.region,
                    x_inc=self.x_inc,
                    y_inc=self.y_inc,
                    verbose=self.verbose,
                    #remote=True
                )
                
            for xyz in _ds.yield_xyz():
                yield(xyz)
        else:
            utils.echo_error_msg(status)

    def yield_array(self, entry):
        """yield the array data from the arcticdem fetch module"""

        status = f_utils.Fetch(entry[0], callback=self.callback, verbose=self.verbose).fetch_file(entry[1])
        if status == 0:
            src_dems = utils.p_untar(entry[1], ['reg_dem.tif'])
            for src_dem in src_dems:
                #src_srs = demfun.get_srs(src_dem)
                _ds = datasets.RasterFile(
                    fn=src_dem,
                    data_format=200,
                    src_srs='epsg:3413',
                    dst_srs=self.dst_srs,
                    src_region=self.region,
                    x_inc=self.x_inc,
                    y_inc=self.y_inc,
                    verbose=self.verbose,
                    remote=True
                )
                
            for arr in _ds.yield_array():
                yield(arr)
        else:
            utils.echo_error_msg(status)
            
### End
