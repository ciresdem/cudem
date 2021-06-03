### copernicus.py - COPERNICUS fetch
##
## Copyright (c) 2010 - 2021 CIRES Coastal DEM Team
##
## copernicus.py is part of CUDEM
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
### Code:

import os
import sys
import lxml.etree
from cudem import utils
from cudem import regions
from cudem import datasets
import cudem.fetches.utils as f_utils
import cudem.fetches.FRED as FRED

## =============================================================================
##
## COPERNICUS Fetch ()
##
## =============================================================================
class CopernicusDEM(f_utils.FetchModule):
    '''Fetch COPERNICUS data'''
    
    def __init__(self, where=[], **kwargs):
        super().__init__(**kwargs)

        self.cop30_url = 'https://opentopography.s3.sdsc.edu/minio/raster/COP30/COP30_hh/'
        self.cop90_url = ''
        self.cop_10_url = 'https://gisco-services.ec.europa.eu/dem/copernicus/outD/'

        self.where = where
        
        self._outdir = os.path.join(os.getcwd(), 'copernicus')
        self.name = 'copernicus'
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
        """Crawl the COP30 database and update/generate the COPERNICUS reference vector."""
        
        self.FRED._open_ds(1)
        surveys = []
        page = f_utils.Fetch(self.cop_10_url).fetch_html()
        rows = page.xpath('//a[contains(@href, ".zip")]/@href')

        for row in rows:
            spat = row.split('.')[0].split('_')[-1]
            x = int(spat.split('x')[-1])
            y = int(spat.split('x')[0].split('y')[-1])
            this_region = regions.Region().from_list([x, x + 10, y, y + 10])

            #geom = this_xml.bounds(geom=True)
            geom = this_region.export_as_geom()
            if geom is not None:
                surveys.append({'Name': row.split('.')[0], 'ID': row.split('.')[0], 'Agency': 'EU', 'Date': utils.this_date(),
                                'MetadataLink': '', 'MetadataDate': utils.this_date(), 'DataLink': self.cop_10_url + row,
                                'DataType': 'dem', 'DataSource': 'copernicus', 'HorizontalDatum': 4326,
                                'VerticalDatum': 'msl', 'Info': '', 'geom': geom})
        self.FRED._add_surveys(surveys)
        self.FRED._close_ds()

    def run(self):
        '''Run the COPERNICUS DEM fetching module'''
        
        for surv in FRED._filter_FRED(self):
            for i in surv['DataLink'].split(','):
                if i != '':
                    self.results.append([i, i.split('/')[-1], surv['DataType']])
        return(self)

    def yield_xyz(self, entry):
        """yield the xyz data from the copernicus fetch module"""
        
        if f_utils.Fetch(entry[0], callback=self.callback, verbose=self.verbose).fetch_file(entry[1]) == 0:
            src_cop_dems = p_unzip(entry[1], ['tif'])
            for src_cop_dem in src_cop_dems:
                _ds = datasets.RasterFile(fn=src_cop_dem, data_format=200, epsg=4326, warp=self.warp,
                                          name=src_cop_dem, src_region=self.region, verbose=self.verbose)
                for xyz in _ds.yield_xyz():
                    if xyz.z != 0:
                        yield(xyz)
                utils.remove_glob(src_cop_dem)
        utils.remove_glob(entry[1])

if __name__ == '__main__':
    cop_dem = CopernicusDEM()
    
### End
