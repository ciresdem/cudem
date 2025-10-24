### usiei.py
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

import json
from cudem.fetches import fetches

## USIEI
class USIEI(fetches.FetchModule):
    """US Interagency Elevation Inventory

    No data is fetched with this module. Will list out query results from the USIEI.
    Set 'want_geometry' to True to output a geojson formatted vector.

    Fields:

    OBJECTID ( type: esriFieldTypeOID, alias: OBJECTID )
    Shape ( type: esriFieldTypeGeometry, alias: Shape )
    ID ( type: esriFieldTypeInteger, alias: ID )
    Title ( type: esriFieldTypeString, alias: Title, length: 250 )
    DataType ( type: esriFieldTypeString, alias: DataType, length: 20 )
    Status ( type: esriFieldTypeString, alias: Status, length: 20 )
    Links ( type: esriFieldTypeString, alias: Links, length: 8000 )
    pointspacing ( type: esriFieldTypeString, alias: pointspacing, length: 50 )
    verticalaccuracy ( type: esriFieldTypeString, alias: verticalaccuracy, length: 8000 )
    horizontalaccuracy ( type: esriFieldTypeString, alias: horizontalaccuracy, length: 300 )
    collectiondate ( type: esriFieldTypeString, alias: collectiondate, length: 200 )
    collectionyear ( type: esriFieldTypeSmallInteger, alias: collectionyear )
    InfoContact ( type: esriFieldTypeString, alias: InfoContact, length: 500 )
    qualitylevel ( type: esriFieldTypeSmallInteger, alias: qualitylevel )
    meets3dep ( type: esriFieldTypeString, alias: meets3dep, length: 50 )
    reasons3dep ( type: esriFieldTypeString, alias: reasons3dep, length: 150 )
    meets3dep_lpc ( type: esriFieldTypeString, alias: meets3dep_lpc, length: 255 )
    productsavailable ( type: esriFieldTypeString, alias: productsavailable, length: 300 )
    pointclasses ( type: esriFieldTypeString, alias: pointclasses, length: 1000 )
    verticaldatum ( type: esriFieldTypeString, alias: verticaldatum, length: 300 )
    horizontaldatum ( type: esriFieldTypeString, alias: horizontaldatum, length: 300 )
    restrictions ( type: esriFieldTypeString, alias: restrictions, length: 20 )
    leafOnOff ( type: esriFieldTypeString, alias: leafOnOff, length: 255 )
    AlternateTitle ( type: esriFieldTypeString, alias: AlternateTitle, length: 500 )
    notes ( type: esriFieldTypeString, alias: notes, length: 500 )
    RecordOwner ( type: esriFieldTypeString, alias: RecordOwner, length: 100 )
    ContractSpec ( type: esriFieldTypeString, alias: ContractSpec, length: 255 )
    Shape_Length ( type: esriFieldTypeDouble, alias: Shape_Length )
    Shape_Area ( type: esriFieldTypeDouble, alias: Shape_Area )
    
    layers:

    0 - Lidar-Topobathy
    1 - Lidar-Bathy
    2 - Lidar-Topo
    3 - IfSAR/InSAR
    4 - Other Bathy

    https://coast.noaa.gov/inventory/

    < usiei:where=None:layer=0:want_geometry=False >
    """
    
    def __init__(self, where = '1=1', want_geometry = False, layer = 0, **kwargs):
        super().__init__(name='usiei', **kwargs)
        self.where = where
        self.want_geometry = want_geometry

        ## The various USIEI URLs
        self._usiei_api_url = ('https://coast.noaa.gov/arcgis/rest/services/'
                               'USInteragencyElevationInventory/USIEIv2/MapServer')
        self._usiei_query_url = '{0}/{1}/query?'.format(self._usiei_api_url, layer)

        
    def run(self):
        """Run the USIEI fetches module"""
        
        if self.region is None:
            return([])
        
        _data = {
            'where': self.where,
            'outFields': '*',
            'geometry': self.region.format('bbox'),
            'inSR':4326,
            'outSR':4326,
            'f':'pjson',
            'returnGeometry':'True' if self.want_geometry else 'False',
        }
        _req = fetches.Fetch(
            self._usiei_query_url,
            verbose=self.verbose
        ).fetch_req(params=_data)
        if _req is not None:
            if self.want_geometry:
                print(_req.text)
            else:
                features = _req.json()
                print(json.dumps(features['features'], indent=True))
            
        return(self)

### End
