### vdatum.py - vdatum conversion grids
##
## Copyright (c) 2010 - 2021 CIRES Coastal DEM Team
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
## Fetch VDATUM transformation grids from NOAA
##
### Code:

import os
import json

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
            
        xmin = float(_inf_areas[key]['minlon'])
        xmax = float(_inf_areas[key]['maxlon'])
        ymin = float(_inf_areas[key]['minlat'])
        ymax = float(_inf_areas[key]['maxlat'])

        if xmin == 0:
            xmin = -180
        else:
            xmin = ((xmin + 180) % 360) - 180

        if xmax == 360:
            xmax = 180
        else:
            xmax = ((xmax + 180) % 360) - 180
        
        _inf_areas_fmt[_out_key]['region'] = [xmin, xmax, ymin, ymax]
        _inf_areas_fmt[_out_key]['grid'] = _inf_areas[key]['source'].split('\\')[-1]

            
    return(_inf_areas_fmt)

## =============================================================================
##
## VDATUM Fetch - NOAA VDatum conversion grids
##
## Fetch vertical datum conversion grids from NOAA's VDATUM project
##
## =============================================================================
class VDATUM(f_utils.FetchModule):
    """Fetch vertical datum conversion grids from NOAA"""

    def __init__(self, where=[], datatype=None, gtx=False, **kwargs):
        super().__init__(**kwargs)
        
        self._vdatum_data_url = 'https://vdatum.noaa.gov/download/data/'
        self._outdir = os.path.join(os.getcwd(), 'vdatum')

        ## add others IGLD85
        self._vdatums = ['VERTCON', 'EGM1984', 'EGM1996', 'EGM2008', 'GEOID03', 'GEOID06', 'GEOID09', 'GEOID12A', 'GEOID12B', 'GEOID96', 'GEOID99', 'TIDAL']
        #self._vdatums = ['GEOID03', 'TIDAL']
        self._tidal_datums = ['mhw', 'mhhw', 'mlw', 'mllw', 'tss', 'mtl']
        self.where = where
        self.datatype = datatype
        self.gtx = gtx
        self.name = 'vdatum'
        self._info = '''Vertical datum transformation grids'''
        self._title = '''NOAA VDatum'''
        self._usage = '''< vdatum >'''
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

        #vertical_datums = {}
        
        self.FRED._open_ds(1)
        for vd in self._vdatums:
            surveys = []

            if vd == 'TIDAL' or vd == 'IGLD85':
                #continue
                ## All tidal inf data is in each one, so we only
                ## have to download one of the tidal zips to process
                ## them all; lets use the smallest one
                ## Keep this link up-tod-date!
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
            utils.remove_glob(*v_infs, '{}.zip'.format(vd))
            
        self.FRED._close_ds()

    def run(self):
        """Search for data in the reference vector file.
        If self.gtx is true, will download the zip and extract the 
        appropriate gtx file.
        """

        if self.datatype is not None:
            self.where.append("DataType = '{}'".format(self.datatype))

        for surv in FRED._filter_FRED(self):
            if self.gtx:
                dst_zip = '{}.zip'.format(surv['ID'])
                if f_utils.Fetch(surv['DataLink'], callback=self.callback, verbose=self.verbose).fetch_file(dst_zip) == 0:
                    v_gtxs = utils.p_f_unzip(dst_zip, [surv['Name']])
                    for v_gtx in v_gtxs:
                        os.rename(v_gtx, '{}.gtx'.format(surv['ID']))
                    utils.remove_glob(dst_zip)
            
            self.results.append([surv['DataLink'], '{}.zip'.format(surv['ID']), surv['Name']])
            
    def yield_xyz(self, entry):
        src_zip = entry[1]
        src_gtx = entry[2]
        src_tif = '{}.tif'.format(utils.fn_basename(src_zip, 'zip'))
        if f_utils.Fetch(entry[0], callback=self.callback, verbose=self.verbose).fetch_file(src_zip) == 0:
            v_gtxs = utils.p_f_unzip(src_zip, [src_gtx])
            utils.run_cmd('gdalwarp {} {} --config CENTER_LONG 0'.format(src_gtx, src_tif))
            _ds = datasets.RasterFile(
                fn=src_tif,
                data_format=200,
                warp=self.warp,
                epsg=None,
                name=src_tif,
                src_region=self.region,
                verbose=self.verbose
            )

            for xyz in _ds.yield_xyz():
                yield(xyz)
            
        utils.remove_glob(*v_gtxs, src_tif, src_zip)
            
### End

