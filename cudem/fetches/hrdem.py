### hrdem.py - NOAA Digital Coast fetch
##
## Copyright (c) 2010 - 2023 CIRES Coastal DEM Team
##
## hrdem.py is part of CUDEM
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
## HRDEM Fetch - Canada High Resolution DEM dataset
##
## Fetch Canadian HRDEM data.
## https://open.canada.ca/data/en/dataset/957782bf-847c-4644-a757-e383c0057995#wb-auto-6
##
### Code:

import os

from osgeo import ogr

from cudem import utils
from cudem import regions
from cudem import datasets

import cudem.fetches.utils as f_utils
import cudem.fetches.FRED as FRED

class HRDEM(f_utils.FetchModule):
    """Fetch HRDEM data from Canada (NRCAN)"""
    
    def __init__(self, where='', **kwargs):
        super().__init__(name='hrdem', **kwargs)
        self._hrdem_footprints_url = 'ftp://ftp.maps.canada.ca/pub/elevation/dem_mne/highresolution_hauteresolution/Datasets_Footprints.zip'
        self._hrdem_info_url = 'https://open.canada.ca/data/en/dataset/957782bf-847c-4644-a757-e383c0057995#wb-auto-6'
        self.where = [where] if len(where) > 0 else []
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
        v_zip = os.path.basename(self._hrdem_footprints_url)
        status = f_utils.Fetch(self._hrdem_footprints_url, verbose=self.verbose).fetch_ftp_file(v_zip)
        v_shps = utils.p_unzip(v_zip, ['shp', 'shx', 'dbf', 'prj'])
        v_shp = None
        for v in v_shps:
            if '.shp' in v: v_shp = v
        shp_regions = regions.gdal_ogr_regions(v_shp)
        shp_region = regions.Region()
        for this_region in shp_regions:
            if shp_region.valid_p(check_xy=True):
                shp_region = regions.regions_merge(shp_region, this_region)
            else: shp_region = this_region
        geom = shp_region.export_as_geom()
        
        self.FRED._attribute_filter(["ID = '{}'".format('HRDEM-1')])
        if self.FRED.layer is None or len(self.FRED.layer) == 0:
            self.FRED._add_survey(Name = 'High-Resolution DEM (Canada)', ID = 'HRDEM-1', Agency = 'NRCAN', Date = utils.this_year(),
                                  MetadataLink = self._hrdem_info_url, MetadataDate = utils.this_year(),
                                  DataLink = self._hrdem_footprints_url, IndexLink = self._hrdem_footprints_url,
                                  DataType = 'raster', DataSource = 'hrdem', Info = 'Canada Only', geom = geom)
        utils.remove_glob(v_zip, *v_shps)
        self.FRED._close_ds()

    def run(self):
        for surv in FRED._filter_FRED(self):
            status = f_utils.Fetch(surv['IndexLink']).fetch_ftp_file(v_zip, verbose=self.verbose)
            v_shps = utils.p_unzip(v_zip, ['shp', 'shx', 'dbf', 'prj'])
            v_shp = None
            for v in v_shps:
                if v.split('.')[-1] == 'shp':
                    v_shp = v
                    break
            try:
                v_ds = ogr.Open(v_shp)
            except:
                v_ds = None
                status = -1
            if v_ds is not None:
                layer = v_ds.GetLayer()
                try:
                    self.FRED.layer.SetAttributeFilter("Name = '{}'".format(name))
                except: pass                
                fcount = layer.GetFeatureCount()
                for f in range(0, fcount):
                    feature = layer[f]
                    if data_link is not None:
                        geom = feature.GetGeometryRef()
                        if geom.Intersects(self.region.export_as_geom()):
                            data_link = feature.GetField('Ftp_dtm').replace('http', 'ftp')
                            self.results.append([data_link, os.path.join(self._outdir, data_link.split('/')[-1]), surv['DataType']])
            utils.remove_glob(v_zip, *v_shps)
            
    ## ==============================================
    ## _update_all() and _parse_results_all() will update FRED with all the data
    ## from the hrdem footprints, which can take a long time and use a lot of FRED
    ## space, which is mostly unncessary, as the footprints are fairly quick to download
    ## and process on the spot for the most part.
    ## use thse *_all functions to revert back to having all the data in the FRED...
    ## ==============================================
    def _update_all(self):
        self.FRED._open_ds(1)
        v_zip = os.path.basename(self._hrdem_footprints_url)
        status = f_utils.Fetch(self._hrdem_footprints_url, verbose=self.verbose).fetch_ftp_file(v_zip)
        v_shps = utils.p_unzip(v_zip, ['shp', 'shx', 'dbf', 'prj'])
        v_shp = None
        for v in v_shps:
            if '.shp' in v: v_shp = v
        try:
            v_ds = ogr.Open(v_shp)
        except:
            v_ds = None
            status = -1
        if v_ds is not None:
            layer = v_ds.GetLayer()
            fcount = layer.GetFeatureCount()
            if self.verbose:
                _prog = utils.CliProgress('scanning {} datasets...'.format(fcount))
            for f in range(0, fcount):
                feature = layer[f]
                name = feature.GetField('Tile_name')
                if self.verbose:
                    _prog.update_perc((f, fcount))
                try:
                    self.FRED.layer.SetAttributeFilter("Name = '{}'".format(name))
                except: pass
                if self.FRED.layer is None or len(self.FRED.layer) == 0:
                    data_link = feature.GetField('Ftp_dtm')
                    if data_link is not None:
                        geom = feature.GetGeometryRef()
                        self.FRED._add_survey(Name = name, ID = feature.GetField('Project'), Agency = 'NRCAN', Date = utils.this_year(),
                                              MetadataLink = feature.GetField('Meta_dtm'), MetadataDate = utils.this_year(),
                                              DataLink = data_link.replace('http', 'ftp'), IndexLink = self._hrdem_footprints_url,
                                              DataType = 'raster', DataSource = 'hrdem', HorizontalDatum = feature.GetField('Coord_Sys').split(':')[-1],
                                              Info = feature.GetField('Provider'), geom = geom)

            if self.verbose:
                _prog.end('scanned {} datasets.'.format(fcount))
        utils.remove_glob(v_zip, *v_shps)
        self.FRED._close_ds()

    def _parse_results_all(self):
        """Parse the results of a filtered FRED"""
        
        for surv in _filter_FRED(self):
            for d in surv['DataLink'].split(','):
                if d != '':
                    self.results.append([d, d.split('/')[-1], surv['DataType']])
                    
    def yield_xyz(self, entry):
        src_dc = os.path.basename(entry[1])
        src_ext = src_dc.split('.')[-1]
        try:
            src_ds = gdal.Open(entry[0])
            src_dc = entry[0]
        except Exception as e:
            f_utils.Fetch(entry[0], callback=self.callback, verbose=self.verbose).fetch_file(src_dc)
            try:
                src_ds = gdal.Open(src_dc)
            except Exception as e:
                echo_error_msg('could not read dc raster file: {}, {}'.format(entry[0], e))
                src_ds = None
        except Exception as e:
            echo_error_msg('could not read dc raster file: {}, {}'.format(entry[0], e))
            src_ds = None

        if src_ds is not None:

            _ds = datasets.RasterFile(
                fn=src_dc,
                data_format=200,
                src_srs='epsg:4326',
                dst_srs=self.dst_srs,
                src_region=self.region,
                x_inc=self.x_inc,
                y_inc=self.y_inc,
                verbose=self.verbose
            )
            _ds.src_ds = src_ds
            _ds.ds_open_p = True
            for xyz in _ds.yield_xyz():
                yield(xyz)
                
        src_ds = None
        utils.remove_glob(src_dc)
                    
### End
