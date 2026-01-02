### hydronos.py
##
## Copyright (c) 2010 - 2026 Regents of the University of Colorado
##
## hydronos.py is part of CUDEM
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
## Fetch NOS Soundings (BAG/Hydro) from NOAA.
##
### Code:

import os
import json
from typing import Optional, Dict, Any, List
from cudem import utils
from cudem.fetches import fetches

## ==============================================
## Constants
## ==============================================
NOS_DYNAMIC_URL = 'https://gis.ngdc.noaa.gov/arcgis/rest/services/web_mercator/nos_hydro_dynamic/MapServer'
NOS_DATA_URL = 'https://data.ngdc.noaa.gov/platforms/ocean/nos/coast/'

## ==============================================
## HydroNOS Module
## ==============================================
class HydroNOS(fetches.FetchModule):
    """NOS Soundings (bag/hydro)
    
    NCEI maintains the National Ocean Service Hydrographic Data Base 
    (NOSHDB) and Hydrographic Survey Meta Data Base (HSMDB).

    Layers:
    0: Surveys with BAGs available.
    1: Surveys with digital sounding data available (including BAGs).

    https://www.ngdc.noaa.gov/mgg/bathymetry/hydro.html
    
    Configuration Example:
    < nos:where='1=1':layer=1:datatype=None:index=False:tables=False >
    """
    
    def __init__(self, where: str = '1=1', layer: int = 1, datatype: Optional[str] = None, 
                 index: bool = False, tables: bool = False, survey_id: Optional[str] = None, 
                 exclude_survey_id: Optional[str] = None, **kwargs):
        super().__init__(name='hydronos', **kwargs)
        self.where = where
        self.datatype = datatype
        self.index = index
        self.tables = tables
        self.survey_id = survey_id
        self.exclude_survey_id = exclude_survey_id

        self._nos_query_url = f'{NOS_DYNAMIC_URL}/{layer}/query?'

        ## Metadata defaults
        self.data_format = None 
        self.src_srs = None
        self.title = 'NOAA Hydrographic Surveys'
        self.source = 'NOAA/NOS'
        self.date = '1933 - 2024'
        self.data_type = 'varies'
        self.resolution = '<1m to several km'
        self.hdatum = 'WGS84'
        self.vdatum = 'MLLW'
        self.url = 'https://gis_ngdc.noaa.gov/'

        
    def run(self):
        """Run the hydronos fetches module."""
        
        if self.region is None:
            return []

        ## Prepare ArcGIS Query
        params = {
            'where': self.where,
            'outFields': '*',
            'geometry': self.region.format('bbox'),
            'inSR': 4326,
            'outSR': 4326,
            'f': 'pjson',
            'returnGeometry': 'False'
        }

        req = fetches.Fetch(
            self._nos_query_url,
            verbose=self.verbose
        ).fetch_req(params=params)

        if req is None:
            return self

        try:
            response = req.json()
        except json.JSONDecodeError:
            utils.echo_error_msg("Failed to parse HydroNOS response.")
            return self

        features = response.get('features', [])
        for feature in features:
            attrs = feature.get('attributes', {})
            
            ## Filter by Year
            year_val = attrs.get('SURVEY_YEAR')
            year = int(year_val) if year_val else 0
            
            if self.min_year is not None and year < self.min_year:
                continue
            if self.max_year is not None and year > self.max_year:
                continue

            ## Handle Index/Table Output
            if self.index:
                print(json.dumps(attrs, indent=4))
                continue
            
            if self.tables:
                self._print_table_row(attrs, year)
                continue

            ## Process Download Links
            self._process_download(attrs, year)

        return self

    
    def _print_table_row(self, attrs: Dict, year: int):
        """Print table formatted output."""
        
        sid = attrs.get('SURVEY_ID')
        has_bag = attrs.get('BAGS_EXIST')
        dtype = 'bag' if has_bag else 'hydro'
        
        line = f"{sid},{year}"
        
        if self.datatype is None or 'bag' in self.datatype.lower():
            if dtype == 'bag':
                print(line)
        else:
            if dtype != 'bag':
                print(line)

                
    def _process_download(self, attrs: Dict, year: int):
        """Process download URL."""
        
        survey_id = attrs.get('SURVEY_ID')
        download_url = attrs.get('DOWNLOAD_URL')
        
        if not download_url:
            return

        ## Filter by Survey ID
        if self.survey_id:
            if survey_id not in self.survey_id.split('/'):
                return
        
        if self.exclude_survey_id:
            if survey_id in self.exclude_survey_id.split('/'):
                return
        
        ## Construct Base Data Link
        ## Example Link: .../platforms/ocean/nos/coast/H12001-H14000/H12345/
        ## Extract directory from original URL
        try:
            nos_dir = download_url.split('/')[-2]
            data_link = f'{NOS_DATA_URL}{nos_dir}/{survey_id}/'
        except IndexError:
            return

        ## Fetch BAGs
        if self.datatype is None or 'bag' in self.datatype.lower():
            bags_exist = str(attrs.get('BAGS_EXIST', '')).upper()
            if bags_exist in ['TRUE', 'Y']:
                self.title = 'NOAA BAG Surveys'
                self.date = str(year)
                
                bag_page = fetches.Fetch(f"{data_link}BAG").fetch_html()
                if bag_page is not None:
                    bags = bag_page.xpath('//a[contains(@href, ".bag")]/@href')
                    for bag in bags:
                        self.add_entry_to_results(
                            f'{data_link}BAG/{bag}',
                            os.path.join('bag', bag),
                            'bag'
                        )

        ## Fetch XYZ (GEODAS)
        if self.datatype is None or 'xyz' in self.datatype.lower():
            self.title = 'NOAA NOS Hydrographic Surveys'
            self.date = str(year)
            
            xyz_page = fetches.Fetch(data_link).fetch_html()
            if xyz_page is not None:
                geodas = xyz_page.xpath('//a[contains(@href, "GEODAS")]/@href')
                if geodas:
                    xyz_filename = f'{survey_id}.xyz.gz'
                    xyz_link = f'{data_link}GEODAS/{xyz_filename}'
                    self.add_entry_to_results(
                        xyz_link,
                        os.path.join('geodas', xyz_filename),
                        'xyz'
                    )

### End
