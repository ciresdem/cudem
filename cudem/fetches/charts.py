### charts.py - NOAA Nautical CHARTS fetch
##
## Copyright (c) 2010 - 2021 CIRES Coastal DEM Team
##
## charts.py is part of CUDEM
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
## Chart Fetch - ENC & RNC
##
## Fetch digital charts from NOAA, including ENC and RNC
##
### Code:

import os
import json

from osgeo import ogr
from cudem import utils
from cudem import regions
from cudem import datasets
from cudem import xyzfun
#from cudem import vdatums

import cudem.fetches.utils as f_utils
import cudem.fetches.FRED as FRED

class NauticalCharts(f_utils.FetchModule):
    """Fetch digital chart data from NOAA"""

    def __init__(self, where=[], datatype=None, **kwargs):
        super().__init__(**kwargs)
        
        self._charts_url = 'https://www.charts.noaa.gov/'
        self._enc_data_catalog = 'https://charts.noaa.gov/ENCs/ENCProdCat_19115.xml'
        self._rnc_data_catalog = 'https://charts.noaa.gov/RNCs/RNCProdCat_19115.xml'
        self._urls = [self._enc_data_catalog, self._rnc_data_catalog]
        self._outdir = os.path.join(os.getcwd(), 'charts')
        self._dt_xml = { 'ENC':self._enc_data_catalog,
                         'RNC':self._rnc_data_catalog }

        self.where = where
        self.datatype = datatype
        self.name = 'charts'
        self.v_datum = 'mhw'
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
        for dt in self._dt_xml.keys():
            surveys = []
            this_xml = f_utils.iso_xml(self._dt_xml[dt], timeout=1000, read_timeout=2000)
            charts = this_xml.xml_doc.findall('.//{*}has', namespaces = this_xml.namespaces)
            if self.verbose:
                _prog = utils.CliProgress('scanning {} surveys in {}.'.format(len(charts), dt))
                
            for i, chart in enumerate(charts):
                this_xml.xml_doc = chart
                title = this_xml.title()
                if self.verbose:
                    _prog.update_perc((i, len(charts)))
                    
                self.FRED._attribute_filter(["ID = '{}'".format(title)])
                if self.FRED.layer is None or len(self.FRED.layer) == 0:
                    h_epsg, v_epsg = this_xml.reference_system()
                    this_data = this_xml.linkages()
                    geom = this_xml.polygon(geom=True)
                    if geom is not None:
                        surveys.append({'Name': title, 'ID': title, 'Agency': 'NOAA', 'Date': this_xml.date(),
                                        'MetadataLink': this_xml.url, 'MetadataDate': this_xml.xml_date(),
                                        'DataLink': this_data, 'Link': self._charts_url, 'DataType': dt,
                                        'DataSource': 'charts', 'HorizontalDatum': h_epsg, 'VerticalDatum': v_epsg,
                                        'Info': this_xml.abstract, 'geom': geom})
                        
            self.FRED._add_surveys(surveys)
            if self.verbose:
                _prog.end(0, 'scanned {} surveys in {}'.format(len(charts), dt))
                utils.echo_msg('added {} surveys from {}'.format(len(surveys), dt))
                
        self.FRED._close_ds()

    def generate_tidal_vdatum(self, src_vdatum, dst_vdatum):
        self.vdatum_grid = vdatums._tidal_transform(self.region, src_vdatum, dst_vdatum)
        
    def run(self):
        """Search for data in the reference vector file"""

        if self.datatype is not None:
            self.where.append("DataType = '{}'".format(self.datatype))
        
        for surv in FRED._filter_FRED(self):
            for i in surv['DataLink'].split(','):
                self.results.append([i, i.split('/')[-1], surv['DataType']])

        #self.generate_tidal_vdatum('mhw', 'tss')

    def yield_xyz(self, entry):
        """ENC data comes as a .000 file in a zip.

The data is referenced to MHW and is represente as a depth.
In U.S. waters, MHW can be transformed to MSL or the local GEOID using
VDatum and/or it's associated grids (mhw.gtx or tss.gtx)"""

        ## create the tidal transformation grid from mhw to geoid
        src_zip = os.path.basename(entry[1])
        if f_utils.Fetch(entry[0], callback=self.callback, verbose=self.verbose).fetch_file(src_zip) == 0:
            if entry[2].lower() == 'enc':
                src_encs = utils.p_unzip(src_zip, ['000'])
                for src_ch in src_encs:
                    dst_xyz = src_ch.split('.')[0] + '.xyz'
                    try:
                        ds_ogr = ogr.Open(src_ch)
                        layer_s = ds_ogr.GetLayerByName('SOUNDG')
                        if layer_s is not None:
                            with open(dst_xyz, 'w') as o_xyz:
                                for f in layer_s:
                                    g = json.loads(f.GetGeometryRef().ExportToJson())
                                    for xyz in g['coordinates']:
                                        xyzfun.XYZPoint().from_list([float(x) for x in xyz]).dump(dst_port=o_xyz, encode=False)
                        ds_ogr = layer_s = None
                    except:
                        utils.echo_warning_msg('could not parse {}'.format(src_ch))

                    _ds = datasets.XYZFile(
                        fn=dst_xyz,
                        data_format=168,
                        z_scale=-1,
                        src_srs='epsg:4326',
                        #src_srs='+proj=longlat +datum=WGS84 +geoidgrids=./{}'.format(vdatum_grid),
                        dst_srs=self.dst_srs,
                        name=dst_xyz,
                        src_region=self.region,
                        verbose=self.verbose,
                        remote=True
                    )
                    for xyz in _ds.yield_xyz():
                        yield(xyz)

                utils.remove_glob(dst_xyz, o_xyz, *src_encs)
        utils.remove_glob(src_zip)
    
### End
