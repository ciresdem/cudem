### charts.py
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

from tqdm import tqdm
from cudem.fetches import fetches
from cudem.fetches import FRED

## Charts - ENC/RNC
## the arcgis rest server doesn't filter by bbox for some reason, always returns all data
class Charts(fetches.FetchModule):
    """NOAA Nautical CHARTS

    Here we can fetch either ENC or RNC
    
    Use the NauticalCharts module instead of this one. The arcgis rest 
    server doesn't always work as expected.
    
    https://www.charts.noaa.gov/
    
    < charts:want_rnc=False >
    """
    
    def __init__(self, where='1=1', want_rnc=False, **kwargs):
        super().__init__(name='charts', **kwargs)
        self.where = where
        self.want_rnc = want_rnc

        ## charts URLs
        self._charts_url = ('https://gis.charttools.noaa.gov/arcgis/rest/'
                            'services/MCS/ENCOnline/MapServer/exts/'
                            'MaritimeChartService/MapServer')
        self._charts_query_url = '{0}/queryDatasets?'.format(self._charts_url)

        
    def run(self):
        if self.region is None:
            return([])

        _data = {
            'where': self.where,
            'outFields': '*',
            'geometry': self.region.format('bbox'),
            'inSR':4326,
            'outSR':4326,
            'f':'pjson',
            'returnGeometry':'True',
        }
        _req = fetches.Fetch(
            self._charts_query_url, verbose=self.verbose
        ).fetch_req(params=_data)
        if _req is not None:
            _req_json = _req.json()
            print(len(_req_json['charts']))
            print(_req_json['charts'][0])
            print(_req.text)

            
class NauticalCharts(fetches.FetchModule):
    """NOAA Nautical CHARTS

    Fetch digital chart data from NOAA
    
    set the 'want_rnc' flag to True to fetch RNC along with ENC data
    
    https://www.charts.noaa.gov/
    
    < charts:want_rnc=False >
    """

    def __init__(self, where='', want_rnc=False, tables=False, **kwargs):
        super().__init__(name='charts', **kwargs)
        self.where = [where] if len(where) > 0 else []
        self.want_rnc = want_rnc
        self.tables = tables
        
        ## various charts URLs
        self._charts_url = 'https://www.charts.noaa.gov/'
        self._enc_data_catalog = 'https://charts.noaa.gov/ENCs/ENCProdCat_19115.xml'
        self._rnc_data_catalog = 'https://charts.noaa.gov/RNCs/RNCProdCat_19115.xml'
        self._ienc_charts_data_catalog = ('https://ienccloud.us/ienc/products/'
                                          'catalog/IENCU37ProductsCatalog.xml')
        self._ienc_buoys_data_catalog = ('https://ienccloud.us/ienc/products/'
                                         'catalog/IENCBuoyProductsCatalog.xml')
        self._urls = [self._enc_data_catalog, self._rnc_data_catalog]
        self._dt_xml = {'ENC':self._enc_data_catalog, 'RNC':self._rnc_data_catalog}
        
        ## for dlim, ENC data comes as .000 files, parse with OGR
        self.src_srs='epsg:4326+5866'
        self.title = 'NOAA OCS electronic navigational chart (ENC) extracted soundings'
        self.source = 'NOAA/NOS'
        self.date = '1966 - 2024'
        self.data_type = 'Digitized Bathymetric Charts'
        self.resolution = '<10m to several kilometers'
        self.hdatum = 'WGS84'
        self.vdatum = 'MLLW'
        self.url = self._charts_url

        ## Charts is in FRED, set that up here.
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
        """Update or create the reference vector file in FRED"""
        
        self.FRED._open_ds(1)
        for dt in self._dt_xml.keys():
            surveys = []
            this_xml = fetches.iso_xml(
                self._dt_xml[dt], timeout=1000, read_timeout=2000
            )
            charts = this_xml.xml_doc.findall(
                './/{*}has',
                namespaces=this_xml.namespaces
            )
            with tqdm(
                    total=len(charts),
                    desc=f'scanning for CHARTS ({dt}) datasets',
                    leave=self.verbose
            ) as pbar:
                for i, chart in enumerate(charts):
                    pbar.update(1)
                    this_xml.xml_doc = chart
                    title = this_xml.title()
                    self.FRED._attribute_filter(["ID = '{}'".format(title)])
                    if self.FRED.layer is None or len(self.FRED.layer) == 0:
                        h_epsg, v_epsg = this_xml.reference_system()
                        this_data = this_xml.linkages()
                        geom = this_xml.polygon(geom=True)
                        if geom is not None:
                            surveys.append(
                                {
                                    'Name': title,
                                    'ID': title,
                                    'Agency': 'NOAA',
                                    'Date': this_xml.date(),
                                    'MetadataLink': this_xml.url,
                                    'MetadataDate': this_xml.xml_date(),
                                    'DataLink': this_data,
                                    'Link': self._charts_url,
                                    'DataType': dt,
                                    'DataSource': 'charts',
                                    'HorizontalDatum': h_epsg,
                                    'VerticalDatum': v_epsg,
                                    'Info': this_xml.abstract,
                                    'geom': geom
                                }
                            )
                        
            self.FRED._add_surveys(surveys)
                
        self.FRED._close_ds()

        
    def generate_tidal_vdatum(self, src_vdatum, dst_vdatum):
        self.vdatum_grid = vdatums._tidal_transform(
            self.region, src_vdatum, dst_vdatum
        )

        
    def run(self):
        """Run the NauticalCharts fetches module. 

        Search for data in the reference vector file (FRED).
        """

        if self.want_rnc:
            self.where.append("DataType = 'RNC'")
        else:
            self.where.append("DataType = 'ENC'")

        _results = FRED._filter_FRED(self)
        with tqdm(
                total=len(_results),
                desc='scanning CHARTS datasets',
                leave=self.verbose
        ) as pbar:
            for surv in _results:
                pbar.update(1)
                for i in surv['DataLink'].split(','):
                    if self.tables:
                        pass
                    else:
                        self.add_entry_to_results(
                            i, i.split('/')[-1], surv['DataType']
                        )

### End
