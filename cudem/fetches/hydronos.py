### hydronos.py
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
from cudem import utils
from cudem.fetches import fetches

## NOAA NOS
class HydroNOS(fetches.FetchModule):
    """NOS Soundings (bag/hydro)
    
    NCEI maintains the National Ocean Service Hydrographic Data Base 
    (NOSHDB) and Hydrographic Survey Meta Data Base (HSMDB). Both 
    are populated by the Office of Coast Survey and National Geodetic 
    Service, and provide coverage of coastal waters and the U.S. 
    exclusive economic zone and its territories. 

    Fields:

    SURVEY_ID ( type: esriFieldTypeString, alias: Survey ID, length: 10 )
    DATE_SURVEY_BEGIN ( type: esriFieldTypeDate, alias: Begin Date, length: 8 )
    DATE_SURVEY_END ( type: esriFieldTypeDate, alias: End Date, length: 8 )
    DATE_MODIFY_DATA ( type: esriFieldTypeDate, alias: Modify Data Date, length: 8 )
    DATE_SURVEY_APPROVAL ( type: esriFieldTypeDate, alias: Survey Approval Date, length: 8 )
    DATE_ADDED ( type: esriFieldTypeDate, alias: Date Added, length: 8 )
    SURVEY_YEAR ( type: esriFieldTypeDouble, alias: Survey Year )
    DIGITAL_DATA ( type: esriFieldTypeString, alias: Digital Data?, length: 15 )
    LOCALITY ( type: esriFieldTypeString, alias: Locality, length: 150 )
    SUBLOCALITY ( type: esriFieldTypeString, alias: Sublocality, length: 150 )
    PLATFORM ( type: esriFieldTypeString, alias: Platform Name, length: 150 )
    PRODUCT_ID ( type: esriFieldTypeString, alias: Product ID, length: 24 )
    BAGS_EXIST ( type: esriFieldTypeString, alias: BAGS_EXIST, length: 4 )
    DOWNLOAD_URL ( type: esriFieldTypeString, alias: Download URL, length: 256 )
    DECADE ( type: esriFieldTypeDouble, alias: Decade )
    PUBLISH ( type: esriFieldTypeString, alias: PUBLISH, length: 1 )
    OBJECTID ( type: esriFieldTypeOID, alias: OBJECTID )
    SHAPE ( type: esriFieldTypeGeometry, alias: SHAPE )
    
    Layer 0: Surveys with BAGs available (Bathymetric Attributed Grids).
    Layer 1: Surveys with digital sounding data available for 
             download (including those with BAGs).

    https://www.ngdc.noaa.gov/mgg/bathymetry/hydro.html
    
    < nos:where=None:layer=0:datatype=None:index=False:tables=False >
    """
    
    def __init__(self, where='1=1', layer=1, datatype=None, index=False,
                 tables=False, survey_id=None, exclude_survey_id=None,
                 **kwargs):
        super().__init__(name='hydronos', **kwargs)
        self.where = where
        self.datatype = datatype
        self.index = index
        self.tables = tables
        self.survey_id = survey_id
        self.exclude_survey_id = exclude_survey_id

        ## various NOS URLs
        self._nos_dynamic_url = ('https://gis.ngdc.noaa.gov/arcgis/rest/services/'
                                 'web_mercator/nos_hydro_dynamic/MapServer')
        self._nos_url = ('https://gis.ngdc.noaa.gov/arcgis/rest/services/'
                         'web_mercator/nos_hydro/MapServer')
        self._nos_data_url = 'https://data.ngdc.noaa.gov/platforms/ocean/nos/coast/'
        self._nos_query_url = '{0}/{1}/query?'.format(self._nos_dynamic_url, layer)

        ## for dlim
        self.data_format = None # bag/xyz data are different, reset later
        self.src_srs = None # bag/xyz data are different, reset later

        self.title = 'NOAA Hydrographic Surveys'
        self.source = 'NOAA/NOS'
        self.date = '1933 - 2024'
        self.data_type = 'varies'
        self.resolution = '<1m to several km'
        self.hdatum = 'WGS84'
        self.vdatum = 'MLLW'
        self.url = 'https://gis_ngdc.noaa.gov/'

        
    def run(self):
        """Run the hydronos fetches module"""
        
        if self.region is None:
            return([])

        _data = {
            'where': self.where,
            'outFields': '*',
            'geometry': self.region.format('bbox'),
            'inSR':4326,
            'outSR':4326,
            'f':'pjson',
            'returnGeometry':'False'
        }
        _req = fetches.Fetch(
            self._nos_query_url,
            verbose=self.verbose
        ).fetch_req(params=_data)
        if _req is not None:
            features = _req.json()
            if 'features' in features.keys():
                for feature in features['features']:
                    if self.index:
                        print(json.dumps(feature['attributes'], indent=4))
                    elif self.tables:
                        sid = feature['attributes']['SURVEY_ID']
                        year = utils.int_or(feature['attributes']['SURVEY_YEAR'])
                        url = feature['attributes']['DOWNLOAD_URL']
                        has_bag = feature['attributes']['BAGS_EXIST']
                        dtype = 'bag' if has_bag is not None else 'hydro'
                        line = '{},{}'.format(sid,year)
                        print(line)
                    else:
                        ID = feature['attributes']['SURVEY_ID']
                        link = feature['attributes']['DOWNLOAD_URL']
                        year = utils.int_or(feature['attributes']['SURVEY_YEAR'])
                        if link is None:
                            continue

                        if self.survey_id is not None:
                            if ID not in self.survey_id.split('/'):
                                continue

                        if self.exclude_survey_id is not None:
                            if ID in self.exclude_survey_id.split('/'):
                                continue
                        
                        nos_dir = link.split('/')[-2]
                        data_link = '{}{}/{}/'.format(self._nos_data_url, nos_dir, ID)
                        
                        if self.datatype is None or 'bag' in self.datatype.lower():
                            self.title = 'NOAA BAG Surveys'
                            self.date = year
                            if feature['attributes']['BAGS_EXIST'] == 'TRUE' \
                               or feature['attributes']['BAGS_EXIST'] == 'Y':
                                page = fetches.Fetch(data_link + 'BAG').fetch_html()
                                bags = page.xpath('//a[contains(@href, ".bag")]/@href')
                                [self.add_entry_to_results(
                                    '{0}BAG/{1}'.format(data_link, bag),
                                    os.path.join('bag', bag),
                                    'bag'
                                ) for bag in bags]

                        if self.datatype is None or 'xyz' in self.datatype.lower():
                            self.title = 'NOAA NOS Hydrographic Surveys'
                            self.date = year
                            page = fetches.Fetch(data_link).fetch_html()
                            if page is not None:
                                geodas = page.xpath('//a[contains(@href, "GEODAS")]/@href')
                                if geodas:
                                    xyz_link = data_link + 'GEODAS/{0}.xyz.gz'.format(ID)
                                    self.add_entry_to_results(
                                        xyz_link,
                                        os.path.join('geodas', xyz_link.split('/')[-1]),
                                        'xyz'
                                    )

### End
