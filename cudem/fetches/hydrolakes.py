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
    
    def __init__(self, where='1=1', layer=0, **kwargs):
        super().__init__(**kwargs)
        self._hydrolakes_prods = 'https://www.hydrosheds.org/products/hydrolakes'
        self._hydrolakes_poly_zip = 'https://data.hydrosheds.org/file/hydrolakes/HydroLAKES_polys_v10_shp.zip'
        self._hydrolakes_gdb_zip = 'https://data.hydrosheds.org/file/hydrolakes/HydroLAKES_polys_v10_gdb.zip'
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
            
    def run(self):
        self.results.append(
            [self._hydrolakes_poly_zip,
             os.path.join(
                 self._outdir,
                 self._hydrolakes_poly_zip.split('/')[-1]
             ),
             'hydrolakes']
        )
        return(self)

    def yield_xyz(self, entry):

        ## use globathy.py to get lake depths and return those
        #pass
        
        if f_utils.Fetch(entry[0], callback=self.callback, verbose=self.verbose, headers=self.headers).fetch_file(entry[1]) == 0:
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

    
### End
