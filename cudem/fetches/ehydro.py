### ehydro.py
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

import time
import json
from cudem import utils
from cudem.fetches import fetches

## eHydro (USACE)
class eHydro(fetches.FetchModule):
    """USACE eHydro bathymetric data.
    
    Maintenance responsibility for more than 25,000 miles of navigation 
    channels and 400 ports and harbors throughout the United States requires 
    extensive surveying and mapping services, including boundary, topographic, 
    hydrographic, terrestrial lidar, and multispectral and hyperspectral aerial 
    imagery collection as well as airborne topographic and bathymetric lidar 
    acquisition, project-level GIS implementation, development of file-based 
    geodatabases, and GIS tool development.

    Three representative survey and mapping datasets include the National Channel 
    Framework (NCF)—an enterprise geodatabase of information on all 61 
    USACE-maintained high-tonnage channels —hydrographic surveys, which provide 
    assistance in locating navigable channels, determining dredging requirements, 
    verifying dredging accuracy, and maintaining harbors and rivers —and Inland 
    Electronic Navigational Charts(IENC), accurate navigational charts provided 
    in a highly structured data format for use in navigation systems and to increase 
    overall navigational safety.. 

    https://navigation.usace.army.mil/Survey/Hydro

    Fields:

    objectid (type: esriFieldTypeOID, alias: objectid, SQL Type: sqlTypeOther, length: 0, nullable: false, editable: false)
    surveyjobidpk (type: esriFieldTypeString, alias: SURVEYJOBIDPK, SQL Type: sqlTypeOther, length: 100, nullable: true, editable: true)
    sdsid (type: esriFieldTypeString, alias: SDSID, SQL Type: sqlTypeOther, length: 40, nullable: true, editable: true)
    sdsfeaturename (type: esriFieldTypeString, alias: SDSFEATURENAME, SQL Type: sqlTypeOther, length: 80, nullable: true, editable: true)
    sdsmetadataid (type: esriFieldTypeString, alias: SDSMETADATAID, SQL Type: sqlTypeOther, length: 80, nullable: true, editable: true)
    surveytype (type: esriFieldTypeString, alias: SURVEYTYPE, SQL Type: sqlTypeOther, length: 26, nullable: true, editable: true)
    channelareaidfk (type: esriFieldTypeString, alias: CHANNELAREAIDFK, SQL Type: sqlTypeOther, length: 50, nullable: true, editable: true)
    dateuploaded (type: esriFieldTypeDate, alias: dateUploaded, SQL Type: sqlTypeOther, length: 8, nullable: true, editable: true)
    usacedistrictcode (type: esriFieldTypeString, alias: usaceDistrictCode, SQL Type: sqlTypeOther, length: 20, nullable: true, editable: true)
    surveydatestart (type: esriFieldTypeDate, alias: SURVEYDATESTART, SQL Type: sqlTypeOther, length: 8, nullable: true, editable: true)
    surveydateend (type: esriFieldTypeDate, alias: SURVEYDATEEND, SQL Type: sqlTypeOther, length: 8, nullable: true, editable: true)
    sourcedatalocation (type: esriFieldTypeString, alias: SOURCEDATALOCATION, SQL Type: sqlTypeOther, length: 255, nullable: true, editable: true)
    sourceprojection (type: esriFieldTypeString, alias: SOURCEPROJECTION, SQL Type: sqlTypeOther, length: 75, nullable: true, editable: true)
    mediaidfk (type: esriFieldTypeString, alias: MEDIAIDFK, SQL Type: sqlTypeOther, length: 100, nullable: true, editable: true)
    projectedarea (type: esriFieldTypeDouble, alias: PROJECTEDAREA, SQL Type: sqlTypeOther, nullable: true, editable: true)
    sdsfeaturedescription (type: esriFieldTypeString, alias: SDSFEATUREDESCRIPTION, SQL Type: sqlTypeOther, length: 255, nullable: true, editable: true)
    dateloadedenterprise (type: esriFieldTypeDate, alias: dateLoadedEnterprise, SQL Type: sqlTypeOther, length: 8, nullable: true, editable: true)
    datenotified (type: esriFieldTypeDate, alias: dateNotified, SQL Type: sqlTypeOther, length: 8, nullable: true, editable: true)
    sourcedatacontent (type: esriFieldTypeString, alias: sourceDataContent, SQL Type: sqlTypeOther, length: 1000, nullable: true, editable: true)
    plotsheetlocation (type: esriFieldTypeString, alias: PLOTSHEETLOCATION, SQL Type: sqlTypeOther, length: 255, nullable: true, editable: true)
    sourceagency (type: esriFieldTypeString, alias: SOURCEAGENCY, SQL Type: sqlTypeOther, length: 20, nullable: true, editable: true)
    globalid (type: esriFieldTypeGlobalID, alias: GlobalID, SQL Type: sqlTypeOther, length: 38, nullable: false, editable: false)
    Shape__Area (type: esriFieldTypeDouble, alias: Shape__Area, SQL Type: sqlTypeDouble, nullable: true, editable: false)
    Shape__Length (type: esriFieldTypeDouble, alias: Shape__Length, SQL Type: sqlTypeDouble, nullable: true, editable: false)
        
    < ehydro:where=None:inc=None:index=False:tables=False >
    """

    def __init__(self, where='1=1', inc=None, survey_name=None, index=False,
                 tables=False, min_year=None, max_year=None, **kwargs):
        super().__init__(name='ehydro', **kwargs)
        self.where = where
        self.survey_name = survey_name
        self.inc = utils.str2inc(inc)
        self.index = index
        self.tables = tables
        self.min_year = utils.int_or(min_year)
        self.max_year = utils.int_or(max_year)
        
        ## Various EHydro URLs
        self._ehydro_gj_api_url = ('https://opendata.arcgis.com/datasets/'
                                   '80a394bae6b547f1b5788074261e11f1_0.geojson')
        self._ehydro_api_url = ('https://services7.arcgis.com/n1YM8pTrFmm7L4hs/'
                                'arcgis/rest/services/eHydro_Survey_Data/FeatureServer/0')
        self._ehydro_query_url = '{0}/query?'.format(self._ehydro_api_url)

        
    def run(self):
        """Run the eHydro fetching module"""
        
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
        _req = fetches.Fetch(self._ehydro_query_url).fetch_req(params=_data)
        if _req is not None:
            features = _req.json()
            if 'features' in features.keys():
                for feature in features['features']:
                    sid = feature['attributes']['sdsmetadataid']
                    fetch_fn = feature['attributes']['sourcedatalocation']
                    year = time.gmtime(
                        int(str(feature['attributes']['surveydatestart'])[:10])
                    ).tm_year
                    if self.survey_name is not None:
                        if sid is None:
                            sid = fetch_fn

                        s = [x in sid for x in self.survey_name.split('/')]
                        if not any(s):
                            continue

                    if self.min_year is not None and int(year) < self.min_year:
                        continue

                    if self.max_year is not None and int(year) > self.max_year:
                        continue
                    
                    if self.index:
                        print(json.dumps(feature['attributes'], indent=4))
                    elif self.tables:
                        sid = feature['attributes']['sdsmetadataid']
                        year = time.gmtime(
                            int(str(feature['attributes']['surveydatestart'])[:10])
                        ).tm_year
                        url = feature['attributes']['sourcedatalocation']
                        dtype = feature['attributes']['surveytype']
                        line = '{},{}'.format(sid,year)
                        if sid is not None:
                            print(line)                            
                    else:
                        self.add_entry_to_results(
                            fetch_fn,
                            fetch_fn.split('/')[-1],
                            'ehydro'
                        )
                
        return(self)

### End
