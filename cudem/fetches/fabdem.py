### fabdem.py
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
##
## fetches.py is part of CUDEM
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
##
### Code:

import os
from osgeo import ogr
from cudem import utils
from cudem import regions
from cudem.fetches import fetches
from cudem.fetches import FRED

## FABDEM
class FABDEM(fetches.FetchModule):
    """FABDEM elevation data
    
    FABDEM (Forest And Buildings removed Copernicus DEM) is a global elevation 
    map that removes building and tree height biases from the Copernicus GLO 30 
    Digital Elevation Model (DEM). The data is available at 1 arc second grid 
    spacing (approximately 30m at the equator) for the globe.
    
    https://data.bris.ac.uk/data/dataset/s5hqmjcdj8yo2ibzi9b4ew3sn
    
    < fabdem >
    """
    
    def __init__(self, **kwargs):
        super().__init__(name='fabdem', **kwargs)

        ## The various FABDEM URLs
        self._fabdem_footprints_url = ('https://data.bris.ac.uk/datasets/'
                                       's5hqmjcdj8yo2ibzi9b4ew3sn/FABDEM_v1-2_tiles.geojson')
        self._fabdem_info_url = 'https://data.bris.ac.uk/data/dataset/s5hqmjcdj8yo2ibzi9b4ew3sn'
        self._fabdem_data_url = 'https://data.bris.ac.uk/datasets/s5hqmjcdj8yo2ibzi9b4ew3sn'

        ## for dlim, data_format of -2 is zipfile
        self.data_format = -2
        self.src_srs = 'epsg:4326+3855'

        self.headers = {
            'User-Agent': ('Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:89.0) '
                           'Gecko/20100101 Firefox/89.0'),
            'referer': self._fabdem_info_url
        }

        
    def run(self):
        """run the FABDEM fetches module"""
        ## use the fabdem footprints vector to find the files to fetch
        v_json = os.path.join(self._outdir, os.path.basename(self._fabdem_footprints_url))
        #try:
        status = fetches.Fetch(
            self._fabdem_footprints_url, verbose=self.verbose
        ).fetch_file(v_json)
        v_ds = ogr.Open(v_json)
        #except Exception as e:
        #utils.echo_warning_msg(f'could not fetch fabdem footprints, {e}')
        #v_ds = None
        #status = -1
            
        if v_ds is not None:
            layer = v_ds.GetLayer()
            fcount = layer.GetFeatureCount()
            for f in range(0, fcount):
                feature = layer[f]
                geom = feature.GetGeometryRef()
                if geom.Intersects(self.region.export_as_geom()):
                    zipfile_name = feature.GetField('zipfile_name')
                    zipfile_url = '/'.join([self._fabdem_data_url, zipfile_name])
                    if zipfile_url not in [x['url'] for x in self.results]:
                        self.add_entry_to_results(zipfile_url, zipfile_name, 'raster')
                        
            v_ds = None
                        
        utils.remove_glob(v_json)

        
class FABDEM_FRED(fetches.FetchModule):
    """Fetch FABDEM data via FRED. 
    
    This module is depreciated, now using the remote footprints 
    instead of FRED.
    """
    
    def __init__(self, where='', **kwargs):
        super().__init__(name='fabdem', **kwargs)
        self.where = [where] if len(where) > 0 else []

        ## The various FABDEM urls
        self._fabdem_footprints_url = ('https://data.bris.ac.uk/datasets/'
                                       's5hqmjcdj8yo2ibzi9b4ew3sn/FABDEM_v1-2_tiles.geojson')
        self._fabdem_info_url = 'https://data.bris.ac.uk/data/dataset/s5hqmjcdj8yo2ibzi9b4ew3sn'        

        ## for dlim, data_format of -2 is zipfile
        self.data_format = -2
        self.src_srs = 'epsg:4326+3855'

        ## FABDEM is in FRED, so set that up here.
        self.FRED = FRED.FRED(name=self.name, verbose=self.verbose)
        self.update_if_not_in_FRED()

        
    def update_if_not_in_FRED(self):
        """update the fetches module in FRED if it's not 
        already in there.
        """
        
        self.FRED._open_ds()
        self.FRED._attribute_filter(
            [f"DataSource = '{self.name}'"]
        )
        if len(self.FRED.layer) == 0:
            self.FRED._close_ds()
            self.update()
            
        self.FRED._close_ds()

            
    def update(self):
        """Crawl the FABDEM database and update/generate the 
        FABDEM reference vector in FRED.
        """
        
        self.FRED._open_ds()
        v_json = os.path.basename(self._fabdem_footprints_url)
        try:
            status = fetches.Fetch(
                self._fabdem_footprints_url,
                verbose=self.verbose
            ).fetch_file(v_json)
        except:
            status = -1
            
        shp_regions = regions.gdal_ogr_regions(v_json)
        shp_region = regions.Region()
        for this_region in shp_regions:
            if shp_region.valid_p(check_xy=True):
                shp_region = regions.regions_merge(shp_region, this_region)
            else:
                shp_region = this_region
                
        geom = shp_region.export_as_geom()        
        self.FRED._attribute_filter(["ID = '{}'".format('FABDEM-1')])
        if self.FRED.layer is None or len(self.FRED.layer) == 0:
            self.FRED._add_survey(
                Name = 'FABDEM',
                ID = 'FABDEM-1',
                Agency = 'Univesity of Bristol',
                Date = utils.this_year(),
                MetadataLink = self._fabdem_info_url,
                MetadataDate = utils.this_year(),
                DataLink = self._fabdem_footprints_url,
                IndexLink = self._fabdem_footprints_url,
                DataType = 'raster',
                DataSource = 'fabdem',
                Info = 'Bare Earth Copernicus',
                geom = geom
            )
        utils.remove_glob(v_json)
        self.FRED._close_ds()

            
    def run(self):
        """Run the FABDEM fetches module"""
        
        for surv in FRED._filter_FRED(self):
            v_json = os.path.basename(self._fabdem_footprints_url)
            status = fetches.Fetch(surv['IndexLink']).fetch_file(
                v_json, verbose=self.verbose
            )
            try:
                v_ds = ogr.Open(v_json)
            except Exception as e:
                utils.echo_warning_msg(f'could not open fabdem footprints, {e}')
                v_ds = None
                status = -1
                
            if v_ds is not None:
                layer = v_ds.GetLayer()
                fcount = layer.GetFeatureCount()
                for f in range(0, fcount):
                    feature = layer[f]
                    geom = feature.GetGeometryRef()
                    if geom.Intersects(self.region.export_as_geom()):
                        zipfile_name = feature.GetField('zipfile_name')
                        self.add_entry_to_results(
                            '/'.join([self._fabdem_data_url, zipfile_name]),
                            zipfile_name,
                            'raster'
                        )
                        
            utils.remove_glob(v_zip)

### End
