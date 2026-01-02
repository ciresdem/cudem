### charts.py
##
## Copyright (c) 2010 - 2026 Regents of the University of Colorado
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
### Commentary:
##
## Fetch NOAA Nautical Charts (ENC/RNC).
##
### Code:

from typing import List, Optional, Dict, Any
from cudem import utils
from cudem.fetches import fetches
from cudem.fetches import FRED

## Optional import for vdatums if available in the package structure
try:
    from cudem import vdatums
except ImportError:
    vdatums = None

## ==============================================
## Constants
## ==============================================
ARCGIS_CHARTS_URL = (
    'https://gis.charttools.noaa.gov/arcgis/rest/services/MCS/'
    'ENCOnline/MapServer/exts/MaritimeChartService/MapServer'
)

NOAA_CHARTS_URL = 'https://www.charts.noaa.gov/'
ENC_CATALOG_URL = 'https://charts.noaa.gov/ENCs/ENCProdCat_19115.xml'
RNC_CATALOG_URL = 'https://charts.noaa.gov/RNCs/RNCProdCat_19115.xml'

## ==============================================
## Charts Classes
## ==============================================
class Charts(fetches.FetchModule):
    """NOAA Nautical CHARTS (ArcGIS REST)

    Here we can fetch either ENC or RNC.
    
    Note: The ArcGIS REST server does not always filter by bbox correctly.
    It is recommended to use the NauticalCharts module instead.
   
    Configuration Example:
    < charts:want_rnc=False >
    """
    
    def __init__(self, where: str = '1=1', want_rnc: bool = False, **kwargs):
        super().__init__(name='charts', **kwargs)
        self.where = where
        self.want_rnc = want_rnc
        self._charts_query_url = f'{ARCGIS_CHARTS_URL}/queryDatasets?'

        
    def run(self):
        """Run the ArcGIS REST Charts fetcher."""
        
        if self.region is None:
            return []

        data = {
            'where': self.where,
            'outFields': '*',
            'geometry': self.region.format('bbox'),
            'inSR': 4326,
            'outSR': 4326,
            'f': 'pjson',
            'returnGeometry': 'True',
        }
        
        req = fetches.Fetch(
            self._charts_query_url, 
            verbose=self.verbose
        ).fetch_req(params=data)
        
        if req is not None:
            req_json = req.json()
            if self.verbose:
                utils.echo_msg(f"Found {len(req_json.get('charts', []))} charts.")
            
        return self

            
class NauticalCharts(fetches.FetchModule):
    """NOAA Nautical CHARTS (Catalog XML)

    Fetch digital chart data from NOAA using XML catalogs and FRED.
    
    Set the 'want_rnc' flag to True to fetch RNC along with ENC data.
    
    https://www.charts.noaa.gov/

    Configuration Example:    
    < charts:want_rnc=False >
    """

    def __init__(self, where: str = '', want_rnc: bool = False, tables: bool = False, **kwargs):
        super().__init__(name='charts', **kwargs)
        self.where = [where] if where else []
        self.want_rnc = want_rnc
        self.tables = tables
        
        ## ENC/RNC Catalogs
        self._dt_xml = {
            'ENC': ENC_CATALOG_URL, 
            'RNC': RNC_CATALOG_URL
        }
        
        ## Metadata
        self.src_srs = 'epsg:4326+5866'
        self.title = 'NOAA OCS electronic navigational chart (ENC) extracted soundings'
        self.source = 'NOAA/NOS'
        self.date = '1966 - 2024'
        self.data_type = 'Digitized Bathymetric Charts'
        self.resolution = '<10m to several kilometers'
        self.hdatum = 'WGS84'
        self.vdatum = 'MLLW'
        self.url = NOAA_CHARTS_URL

        ## Initialize FRED
        self.FRED = FRED.FRED(name=self.name, verbose=self.verbose)
        self.update_if_not_in_FRED()

        
    def update_if_not_in_FRED(self):
        """Update the fetches module in FRED if it's not already there."""
        
        self.FRED._open_ds()
        self.FRED._attribute_filter([f"DataSource = '{self.name}'"])
        
        if len(self.FRED.layer) == 0:
            self.FRED._close_ds()
            self.update()
        else:
            self.FRED._close_ds()

            
    def update(self):
        """Update or create the reference vector file in FRED."""
        
        self.FRED._open_ds(1)
        
        try:
            for dt, catalog_url in self._dt_xml.items():
                surveys = []
                
                ## Fetch XML Catalog
                this_xml = fetches.iso_xml(
                    catalog_url, timeout=1000, read_timeout=2000
                )
                
                if this_xml.xml_doc is None:
                    continue
                    
                charts = this_xml.xml_doc.findall(
                    './/{*}has',
                    namespaces=this_xml.namespaces
                )
                
                with utils.ccp(
                        total=len(charts),
                        desc=f'Scanning for CHARTS ({dt}) datasets',
                        leave=self.verbose
                ) as pbar:
                    for chart in charts:
                        pbar.update(1)
                        
                        ## Update XML context to current chart node
                        this_xml.xml_doc = chart
                        title = this_xml.title()
                        
                        ## check if exists
                        self.FRED._attribute_filter([f"ID = '{title}'"])
                        if self.FRED.layer is None or len(self.FRED.layer) == 0:
                            h_epsg, v_epsg = this_xml.reference_system()
                            this_data = this_xml.linkages()
                            geom = this_xml.polygon(geom=True)
                            
                            if geom is not None:
                                surveys.append({
                                    'Name': title,
                                    'ID': title,
                                    'Agency': 'NOAA',
                                    'Date': this_xml.date(),
                                    'MetadataLink': this_xml.url,
                                    'MetadataDate': this_xml.xml_date(),
                                    'DataLink': this_data,
                                    'Link': NOAA_CHARTS_URL,
                                    'DataType': dt,
                                    'DataSource': 'charts',
                                    'HorizontalDatum': h_epsg,
                                    'VerticalDatum': v_epsg,
                                    'Info': this_xml.abstract(),
                                    'geom': geom
                                })
                            
                self.FRED._add_surveys(surveys)
        except:
            utils.echo_error_msg('Could not update FRED')                
        finally:
            self.FRED._close_ds()

            
    def generate_tidal_vdatum(self, src_vdatum: str, dst_vdatum: str):
        """Generate a vertical datum grid (Requires vdatums module)."""
        
        if vdatums:
            self.vdatum_grid = vdatums._tidal_transform(
                self.region, src_vdatum, dst_vdatum
            )
        else:
            utils.echo_warning_msg("vdatums module not available.")

            
    def run(self):
        """Run the NauticalCharts fetches module."""
        
        ## Add Data Type filter
        if self.want_rnc:
            self.where.append("DataType = 'RNC'")
        else:
            self.where.append("DataType = 'ENC'")

        ## Filter FRED
        _results = FRED._filter_FRED(self)
        
        with utils.ccp(
                total=len(_results),
                desc='Scanning CHARTS datasets',
                leave=self.verbose
        ) as pbar:
            for surv in _results:
                pbar.update(1)
                
                if not self.tables and surv.get('DataLink'):
                    ## Handle multiple links if comma-separated
                    links = surv['DataLink'].split(',')
                    for link in links:
                        link = link.strip()
                        if link:
                            self.add_entry_to_results(
                                link, 
                                link.split('/')[-1], 
                                surv['DataType']
                            )
        
        return self

### End
