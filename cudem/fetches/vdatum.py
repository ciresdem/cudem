### vdatum.py - vdatum conversion grids
##
## Copyright (c) 2010 - 2022 Regents of the University of Colorado
##
## vdatum.py is part of CUDEM
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
## VDATUM Fetch - NOAA VDatum conversion grids
##
## Fetch vertical datum conversion grids from NOAA's VDATUM project and PROJ
##
## also here (has more)
## https://www.agisoft.com/downloads/geoids/
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

def proc_vdatum_inf(vdatum_inf, name='vdatum'):
    
    _inf = open(vdatum_inf, 'r')
    _inf_areas = {}

    for line in _inf:
        line_list = line.split('.')
        if len(line_list) > 1:
            if line_list[0] not in _inf_areas.keys():
                _inf_areas[line_list[0]] = {}
            if len(line_list) > 1:
                line_val = '.'.join(line_list[1:]).strip()
                utils.args2dict([line_val], _inf_areas[line_list[0]])
    _inf.close()

    _inf_areas_fmt = {}

    for key in _inf_areas.keys():
        if name is not None:
            _out_key = '{}_{}'.format(name, key)
        else:
            _out_key = key
            
        if _out_key not in _inf_areas_fmt.keys():
            _inf_areas_fmt[_out_key] = {}
            
        xmin = utils.x360(float(_inf_areas[key]['minlon']))
        xmax = utils.x360(float(_inf_areas[key]['maxlon']))
        ymin = float(_inf_areas[key]['minlat'])
        ymax = float(_inf_areas[key]['maxlat'])
        _inf_areas_fmt[_out_key]['region'] = [xmin, xmax, ymin, ymax]
        _inf_areas_fmt[_out_key]['grid'] = _inf_areas[key]['source'].split('\\')[-1]

        if 'AK' in _out_key: print(_inf_areas_fmt[_out_key])
        
    return(_inf_areas_fmt)

def search_proj_cdn(region, epsg=None, crs_name=None, name=None, verbose=True, cache_dir='./'):
    """Search PROJ CDN for transformation grids:
    the PROJ CDN holds transformation grids from around the
    world, including global transformations such as EGM
    """
    
    _proj_vdatum_index = 'https://cdn.proj.org/files.geojson'
    cdn_index = os.path.join(cache_dir, 'proj_cdn_files.geojson')
    #cdn_headers = { 'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:89.0) Gecko/20100101 Firefox/89.0' }
    cdn_headers = {'Host': 'cdn.proj.org',
                   'User-Agent': 'Mozilla/5.0 (X11; Linux x86_64; rv:45.0) Gecko/20100101 Firefox/45.0',
                   'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8',
                   'Accept-Language': 'pt-BR,pt;q=0.8,en-US;q=0.5,en;q=0.3',
                   'Accept-Encoding': 'gzip, deflate',
                   'Connection': 'keep-alive',
                   'Pragma': 'no-cache',
                   'Cache-Control': 'no-cache'}
    
    if f_utils.Fetch(_proj_vdatum_index, headers=cdn_headers, verbose=verbose).fetch_file(cdn_index, timeout=5, read_timeout=5) == 0:
        cdn_driver = ogr.GetDriverByName('GeoJSON')
        cdn_ds = cdn_driver.Open(cdn_index, 0)
        cdn_layer = cdn_ds.GetLayer()
        _boundsGeom = region.export_as_geom()
        _results = []

        if crs_name is not None:
            cdn_layer.SetAttributeFilter(
                "type != 'HORIZONTAL_OFFSET' AND (target_crs_name LIKE '%{}%' OR source_crs_name LIKE '%{}%')".format(
                    name.upper(), name.upper()
                )
            )
        elif epsg is not None:
            cdn_layer.SetAttributeFilter(
                "type != 'HORIZONTAL_OFFSET' AND (target_crs_code LIKE '%{}%' OR source_crs_code LIKE '%{}%')".format(
                    epsg, epsg
                )
            )
        elif name is not None:
            cdn_layer.SetAttributeFilter("type != 'HORIZONTAL_OFFSET' AND name LIKE '%{}%'".format(name))
        else:
            cdn_layer.SetAttributeFilter("type != 'HORIZONTAL_OFFSET'")

        for feat in cdn_layer:
            if _boundsGeom is not None:
                geom = feat.GetGeometryRef()
                if geom is not None:
                    if _boundsGeom.Intersects(geom):
                        _results.append({})
                        f_j = json.loads(feat.ExportToJson())
                        for key in f_j['properties'].keys():
                            _results[-1][key] = feat.GetField(key)
            else:
                _results.append({})
                f_j = json.loads(feat.ExportToJson())
                for key in f_j['properties'].keys():
                    _results[-1][key] = feat.GetField(key)

        cdn_ds = None
        #utils.remove_glob(cdn_index)
        return(_results)

class VDATUM(f_utils.FetchModule):
    """Fetch vertical datum conversion grids from NOAA, etc."""

    _tidal_references = {
        1089: {'name': 'mllw',
               'description': 'Mean Lower Low Water Height',
               'grid': 'mllw.gtx'},
        1091: {'name': 'mlw',
               'description': 'Mean Low Water Height',
               'grid': 'mlw.gtx'},
        5868: {'name': 'mhw',
               'description': 'Mean High Water',
               'grid': 'mhw.gtx'},
        5869: {'name': 'mhhw',
               'description': 'Mean Higher High Water',
               'grid': 'mhhw.gtx'},
        5703: {'name': 'tss',
               'description': 'NAVD88 tss geoid',
               'grid': 'tss.gtx'},
        6641: {'name': 'tss',
               'description': 'PRVD02 tss geoid',
               'grid': 'tss.gtx'},
        6642: {'name': 'tss',
               'description': 'VIVD09 tss geoid',
               'grid': 'tss.gtx'},
    }
    
    def __init__(self, where=[], datatype=None, gtx=False, epsg=None, **kwargs):
        super().__init__(**kwargs)
        
        self._vdatum_data_url = 'https://vdatum.noaa.gov/download/data/'
        self._proj_vdatum_index = 'https://cdn.proj.org/files.geojson'
        
        self._outdir = os.path.join(os.getcwd(), 'vdatum')
        
        ## add others IGLD85
        #self._vdatums = ['VERTCON', 'EGM1984', 'EGM1996', 'EGM2008', 'GEOID03', 'GEOID06', 'GEOID09', 'GEOID12A', 'GEOID12B', 'GEOID96', 'GEOID99', 'TIDAL']
        self._vdatums = ['TIDAL']
        self._tidal_datums = ['mhw', 'mhhw', 'mlw', 'mllw', 'tss', 'mtl', 'msl']
        self.where = where
        self.datatype = datatype
        self.epsg = utils.int_or(epsg)
        self.gtx = gtx
        self.name = 'vdatum'
        self.v_datum = 'varies'
        
        self.FRED = FRED.FRED(name=self.name, verbose=self.verbose)
        self.update_if_not_in_FRED()
        
    def update_if_not_in_FRED(self):
        self.FRED._open_ds()
        self.FRED._attribute_filter(["DataSource = '{}'".format(self.name)])
        if len(self.FRED.layer) == 0:
            self.FRED._close_ds()
            self.update()
            
        self.FRED._close_ds()

    def update(self):
        """Update or create the reference vector file"""

        self.FRED._open_ds(1)
        for vd in self._vdatums:
            surveys = []

            if vd == 'TIDAL' or vd == 'IGLD85':
                ## ==============================================
                ## All tidal inf data is in each one, so we only
                ## have to download one of the tidal zips to process
                ## them all; lets use the smallest one
                ## Keep this link up-to-date!
                ## ==============================================
                if vd == 'TIDAL':
                    vd_ = 'DEVAemb12_8301'
                else:
                    vd_ = vd
                    
                vd_zip_url = '{}{}.zip'.format(self._vdatum_data_url, vd_)
                v_inf = 'tidal_area.inf'
            elif vd == 'VERTCON':
                vd_zip_url = '{}vdatum_{}.zip'.format(self._vdatum_data_url, vd)
                v_inf = 'vcn.inf'
            else:
                vd_zip_url = '{}vdatum_{}.zip'.format(self._vdatum_data_url, vd)
                v_inf = '{}.inf'.format(vd.lower())
                
            if f_utils.Fetch(vd_zip_url, verbose=True).fetch_file('{}.zip'.format(vd)) == 0:
                v_infs = utils.p_unzip('{}.zip'.format(vd), ['inf'])
                v_dict = proc_vdatum_inf(v_inf, name=vd if vd != 'TIDAL' else None)#, loff=-360 if vd =='TIDAL' else -180)
                v_dict = proc_vdatum_inf(v_inf, name=vd if vd != 'TIDAL' else None)#, loff=-360)

                for key in v_dict.keys():
                    v_dict[key]['vdatum'] = vd
                    v_dict[key]['remote'] = vd_zip_url
                
                ## tidal datums:
                if vd == 'TIDAL':
                    v_dict_ = {}
                    for tidal_key in v_dict.keys():                        
                        for t in self._tidal_datums:
                            key_ = '{}_{}'.format(t, tidal_key)
                            v_dict_[key_] = {}
                            v_dict_[key_]['region'] = v_dict[tidal_key]['region']
                            v_dict_[key_]['vdatum'] = t
                            v_dict_[key_]['grid'] = '{}.gtx'.format(t)
                            v_dict_[key_]['remote'] = '{}{}.zip'.format(self._vdatum_data_url, tidal_key)
                            
                    v_dict = v_dict_
                    
                print(v_dict)

                for key in v_dict.keys():
                    self.FRED._attribute_filter(["ID = '{}'".format(key)])
                    if self.FRED.layer is None or len(self.FRED.layer) == 0:
                        geom = regions.Region().from_list(v_dict[key]['region']).export_as_geom()
                        if geom is not None:
                            surveys.append({'Name': v_dict[key]['grid'], 'ID': key, 'Agency': 'NOAA', 'Date': utils.this_date(),
                                            'MetadataLink': "", 'MetadataDate': utils.this_date(),
                                            'DataLink': v_dict[key]['remote'], 'Link': self._vdatum_data_url, 'DataType': v_dict[key]['vdatum'],
                                            'DataSource': 'vdatum', 'HorizontalDatum': 4326, 'VerticalDatum': v_dict[key]['vdatum'],
                                            'Info': "", 'geom': geom})
            self.FRED._add_surveys(surveys)
            #utils.remove_glob(*v_infs, '{}.zip'.format(vd))
            utils.remove_glob(*v_infs)
            
        self.FRED._close_ds()

    def run(self):
        """Search for data in the reference vector file.
        If self.gtx is true, will download the zip and extract the 
        appropriate gtx file.
        """

        w = []
        if self.datatype is not None:
            w.append("DataType = '{}'".format(self.datatype))
        elif self.epsg is not None:
            w.append("DataType = '{}'".format(self._tidal_references[self.epsg]['name']))

        ## ==============================================
        ## Search FRED for VDATUM TIDAL TRANSFORMATION GRIDS
        ## FRED holds the VDATUM tidal grids,
        ## mllw, mlw, mhhw, mhw, tss, mtl
        ## where all convert to MSL and tss represents the
        ## current geoid
        ## ==============================================
        for surv in self.FRED._filter(self.region, w, [self.name]):
            if self.gtx:
                dst_zip = '{}.zip'.format(surv['ID'])
                if f_utils.Fetch(surv['DataLink'], callback=self.callback, verbose=self.verbose).fetch_file(dst_zip) == 0:
                    v_gtxs = utils.p_f_unzip(dst_zip, [surv['Name']])
                    for v_gtx in v_gtxs:
                        os.rename(v_gtx, '{}.gtx'.format(surv['ID']))
                    #utils.remove_glob(dst_zip)
            else:
                self.results.append([surv['DataLink'], os.path.join(self._outdir, '{}.zip'.format(surv['ID'])), surv['Name']])

        ## ==============================================
        ## Search PROJ CDN for all other transformation grids:
        ## the PROJ CDN holds transformation grids from around the
        ## world, including global transformations such as EGM
        ## ==============================================
        cdn_index = 'proj_cdn_files.geojson'
        if f_utils.Fetch(self._proj_vdatum_index, callback=self.callback, verbose=self.verbose).fetch_file(cdn_index) == 0:
            cdn_driver = ogr.GetDriverByName('GeoJSON')
            cdn_ds = cdn_driver.Open(cdn_index, 0)
            cdn_layer = cdn_ds.GetLayer()
            _boundsGeom = self.region.export_as_geom()
            _results = []

            if self.datatype is not None:
                cdn_layer.SetAttributeFilter(
                    "type != 'HORIZONTAL_OFFSET' AND (target_crs_name LIKE '%{}%' OR source_crs_name LIKE '%{}%')".format(
                        self.datatype.upper(), self.datatype.upper()
                    )
                )
            elif self.epsg is not None:
                cdn_layer.SetAttributeFilter(
                    "type != 'HORIZONTAL_OFFSET' AND (target_crs_code LIKE '%{}%' OR source_crs_code LIKE '%{}%')".format(
                        self.epsg, self.epsg
                    )
                )
            else:
                cdn_layer.SetAttributeFilter("type != 'HORIZONTAL_OFFSET'")
             
            for feat in cdn_layer:
                if _boundsGeom is not None:
                    geom = feat.GetGeometryRef()
                    if geom is not None:
                        if _boundsGeom.Intersects(geom):
                            _results.append({})
                            f_j = json.loads(feat.ExportToJson())
                            for key in f_j['properties'].keys():
                                _results[-1][key] = feat.GetField(key)
                else:
                    _results.append({})
                    f_j = json.loads(feat.ExportToJson())
                    for key in f_j['properties'].keys():
                        _results[-1][key] = feat.GetField(key)
                        
            for _result in _results:
                self.results.append([_result['url'], os.path.join(self._outdir, _result['name']), _result['source_id']])
                
            cdn_ds = None
            utils.remove_glob(cdn_index)
            
    def yield_xyz(self, entry):
        src_zip = entry[1]
        src_gtx = entry[2]
        src_tif = '{}.tif'.format(utils.fn_basename(src_zip, 'zip'))
        if f_utils.Fetch(entry[0], callback=self.callback, verbose=self.verbose).fetch_file(src_zip) == 0:
            v_gtxs = utils.p_f_unzip(src_zip, [src_gtx])
            utils.run_cmd('gdalwarp {} {} --config CENTER_LONG 0'.format(src_gtx, src_tif), verbose=True)
            _ds = datasets.RasterFile(
                fn=src_tif,
                data_format=200,
                dst_srs=self.dst_srs,
                src_srs=None,
                src_region=self.region,
                x_inc=self.x_inc,
                y_inc=self.y_inc,
                verbose=self.verbose
            )

            for xyz in _ds.yield_xyz():
                yield(xyz)
            
        utils.remove_glob(*v_gtxs, src_tif, src_zip)
            
### End

