### wsf.py - WSF from DLR (German Aerospace Center)
##
## Copyright (c) 2022, 2023 Regents of the University of Colorado
##
## wsf.py is part of CUDEM
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
## WSF Fetch
##
## https://www.dlr.de/EN/Home/home_node.html
## https://geoservice.dlr.de/web/services
##
### Code:

import os
import sys

from osgeo import gdal

from cudem import utils
from cudem import regions
from cudem import datasets
from cudem import demfun

import cudem.fetches.utils as f_utils
import cudem.fetches.FRED as FRED

class WSF(f_utils.FetchModule):
    """WSF from German Aerospace Service (DLR)

World Settlement Footprint (WSF) 2019

< usiei:where=None:datatype=None >"""

    def __init__(self, where='', datatype=None, **kwargs):
        super().__init__(name='wsf', **kwargs)

        self._wsf_url = 'https://download.geoservice.dlr.de/WSF2019/files/'
        self.where = [where] if len(where) > 0 else []
        self.datatype = datatype
        self.FRED = FRED.FRED(name=self.name, verbose=self.verbose)

        self.headers = { 'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:89.0) Gecko/20100101 Firefox/89.0' }
        self.update_if_not_in_FRED()
        
    def update_if_not_in_FRED(self):
        self.FRED._open_ds()
        self.FRED._attribute_filter(["DataSource = '{}'".format(self.name)])
        if len(self.FRED.layer) == 0:
            self.FRED._close_ds()
            self.update()
            
        self.FRED._close_ds()

    def update(self):
        """Crawl the SWF database and update/generate the COPERNICUS reference vector."""
        
        self.FRED._open_ds(1)
        surveys = []
        page = f_utils.Fetch(self._wsf_url, verbose=True).fetch_html()
        rows = page.xpath('//a[contains(@href, ".tif")]/@href')
        if self.verbose:
            _prog = utils.CliProgress('scanning {} tiles in {}...'.format(len(rows), self._wsf_url))
        
        for i, row in enumerate(rows):
            sid = row.split('.')[0]
            if sid == 'WSF2019_cog':
                continue
            
            if self.verbose:
                _prog.update_perc((i, len(rows)))
                
            self.FRED._attribute_filter(["ID = '{}'".format(sid)])
            if self.FRED.layer is None or len(self.FRED.layer) == 0:
                spat = row.split('.')[0].split('_')
                x = int(spat[-2])
                y = int(spat[-1])
                this_region = regions.Region().from_list(
                    [x, x + 2, y, y + 2]
                )
                geom = this_region.export_as_geom()
                if geom is not None:
                    surveys.append({'Name': row.split('.')[0], 'ID': sid, 'Agency': 'DLR', 'Date': utils.this_date(),
                                    'MetadataLink': row.split('.')[0] + '_stac.json', 'MetadataDate': utils.this_date(), 'DataLink': self._wsf_url + row,
                                    'DataType': 'WSF', 'DataSource': 'WSF', 'HorizontalDatum': 'epsg:4326',
                                    'VerticalDatum': 'None', 'Info': '', 'geom': geom})

        if self.verbose:
            _prog.end(0, 'scanned {} tiles in {}.'.format(len(rows), self._wsf_url))                    
            utils.echo_msg('added {} WSF tiles'.format(len(surveys)))
            
        self.FRED._add_surveys(surveys)
        self.FRED._close_ds()

    def run(self):
        '''Run the WSF fetching module'''

        if self.datatype is not None:
            self.where.append("DataType = '{}'".format(self.datatype))

        for surv in FRED._filter_FRED(self):
            for i in surv['DataLink'].split(','):
                self.results.append([i, os.path.join(self._outdir, i.split('/')[-1].split('?')[0]), surv['DataType']])
                
        return(self)
        
if __name__ == '__main__':
    wsf_fetch = WSF()
    
### End
        
