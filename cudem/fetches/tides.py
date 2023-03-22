### tides.py - TIDE stations fetch
##
## Copyright (c) 2018 - 2023 Regents of the University of Colorado
##
## tides.py is part of CUDEM
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
## Fetch tidal datum station information from NOAA/NOS
##
## Fields:
##     objectid ( type: esriFieldTypeOID , alias: objectid , editable: false , nullable: false )
##     id ( type: esriFieldTypeString , alias: id , editable: true , nullable: true , length: 50 )
##     name ( type: esriFieldTypeString , alias: name , editable: true , nullable: true , length: 50 )
##     affil ( type: esriFieldTypeString , alias: affil , editable: true , nullable: true , length: 50 )
##     latitude ( type: esriFieldTypeDouble , alias: latitude , editable: true , nullable: true )
##     longitude ( type: esriFieldTypeDouble , alias: longitude , editable: true , nullable: true )
##     data ( type: esriFieldTypeString , alias: data , editable: true , nullable: true , length: 200 )
##     dataapi ( type: esriFieldTypeString , alias: dataapi , editable: true , nullable: true , length: 200 )
##     accepted ( type: esriFieldTypeString , alias: accepted , editable: true , nullable: true , length: 50 )
##     epoch ( type: esriFieldTypeString , alias: epoch , editable: true , nullable: true , length: 50 )
##     units ( type: esriFieldTypeString , alias: units , editable: true , nullable: true , length: 50 )
##     orthodatum ( type: esriFieldTypeString , alias: orthodatum , editable: true , nullable: true , length: 50 )
##     mhhw ( type: esriFieldTypeDouble , alias: mhhw , editable: true , nullable: true )
##     mhw ( type: esriFieldTypeDouble , alias: mhw , editable: true , nullable: true )
##     mtl ( type: esriFieldTypeDouble , alias: mtl , editable: true , nullable: true )
##     msl ( type: esriFieldTypeDouble , alias: msl , editable: true , nullable: true )
##     dtl ( type: esriFieldTypeDouble , alias: dtl , editable: true , nullable: true )
##     mlw ( type: esriFieldTypeDouble , alias: mlw , editable: true , nullable: true )
##     mllw ( type: esriFieldTypeDouble , alias: mllw , editable: true , nullable: true )
##     stnd ( type: esriFieldTypeDouble , alias: stnd , editable: true , nullable: true )
##     mn ( type: esriFieldTypeDouble , alias: mn , editable: true , nullable: true )
##     dhq ( type: esriFieldTypeDouble , alias: dhq , editable: true , nullable: true )
##     dlq ( type: esriFieldTypeDouble , alias: dlq , editable: true , nullable: true )
##     hwi ( type: esriFieldTypeDouble , alias: hwi , editable: true , nullable: true )
##     lwi ( type: esriFieldTypeDouble , alias: lwi , editable: true , nullable: true )
##     gt ( type: esriFieldTypeDouble , alias: gt , editable: true , nullable: true )
##     navd88 ( type: esriFieldTypeDouble , alias: navd88 , editable: true , nullable: true )
##     wl_max ( type: esriFieldTypeDouble , alias: wl_max , editable: true , nullable: true )
##     max_date ( type: esriFieldTypeString , alias: max_date , editable: true , nullable: true , length: 50 )
##     wl_min ( type: esriFieldTypeDouble , alias: wl_min , editable: true , nullable: true )
##     min_date ( type: esriFieldTypeString , alias: min_date , editable: true , nullable: true , length: 50 )
#3
### Code:

import os
import json

from cudem import utils
from cudem import regions
from cudem import datasets
from cudem import xyzfun

import cudem.fetches.utils as f_utils

class Tides(f_utils.FetchModule):
    """TIDE station information from NOAA/NOS

Fetch NOS Tide Stations
    
https://tidesandcurrents.noaa.gov/

< tides:station_id=None:s_datum=mllw:t_datum=msl:units=m >"""

    
    def __init__(self, s_datum='mllw', t_datum='msl', units='m', station_id=None, **kwargs):
        super().__init__(name='tides', **kwargs)
        self._stations_api_url = 'https://idpgis.ncep.noaa.gov/arcgis/rest/services/NOS_Observations/CO_OPS_Products/FeatureServer/0/query?'
        self.s_datum = s_datum
        self.t_datum = t_datum
        self.units = units
        self.station_id = station_id
        
    def run(self):
        '''Run the TIDES fetching module'''
        
        if self.region is None:
            return([])
        
        _data = {
            'outFields': '*',
            'units': 'esriSRUnit_Meter',
            'geometry': self.region.format('bbox'),
            'inSR':4326,
            'outSR':4326,
            'f':'pjson',
        }
        _req = f_utils.Fetch(self._stations_api_url, verbose=self.verbose).fetch_req(params=_data)
        if _req is not None:
            self.results.append([_req.url, os.path.join(self._outdir, 'tides_results_{}.json'.format(self.region.format('fn'))), 'tides'])
            
        return(self)

    def yield_xyz(self, entry):
        """process stations"""
        
        src_data = 'tides_tmp.json'
        src_csv = 'tides_tmp.csv'
        ln = 0
        
        if f_utils.Fetch(entry[0], callback=self.callback, verbose=self.verbose).fetch_file(src_data) == 0:
            with open(src_data, 'r') as json_file:
                r = json.load(json_file)

            if len(r) > 0:
                for feature in r['features']:
                    if self.station_id is not None:
                        if self.station_id != feature['attributes']['id']:
                            continue
                    lon = feature['attributes']['longitude']
                    lat = feature['attributes']['latitude']
                    if feature['attributes'][self.s_datum] != -99999.99 and feature['attributes'][self.t_datum] != -99999.99:
                        z = feature['attributes'][self.s_datum] - feature['attributes'][self.t_datum]
                        if self.units == 'm':
                            z = z * 0.3048
                        
                        xyz = xyzfun.XYZPoint(src_srs='epsg:4326').from_list([lon, lat, z])
                        if self.dst_srs is not None:
                            xyz.warp(dst_srs=self.dst_srs)

                        ln += 1
                        yield(xyz)
                    else:
                        utils.echo_error_msg('could not parse difference between {} and {} for {}'.format(self.s_datum, self.t_datum, feature['attributes']['name']))
                    
        else:
            utils.echo_error_msg('failed to fetch remote file, {}...'.format(src_data))

        if self.verbose:
            utils.echo_msg(
                'parsed {} data records from {}'.format(ln, src_data)
            )
            
        utils.remove_glob('{}*'.format(src_data))
### End
