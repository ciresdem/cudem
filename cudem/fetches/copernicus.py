### copernicus.py - COPERNICUS fetch
##
## Copyright (c) 2021 CIRES Coastal DEM Team
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
##
## COPERNICUS Fetch
##
## The Copernicus DEM is a Digital Surface Model (DSM) which represents the surface of the Earth including buildings, 
## infrastructure and vegetation.
##
## Fetches from opentopography
##
### Code:

import os
import sys

from cudem import utils
from cudem import regions
from cudem import datasets

import cudem.fetches.utils as f_utils
import cudem.fetches.FRED as FRED

class CopernicusDEM(f_utils.FetchModule):
    '''Fetch COPERNICUS data'''
    
    def __init__(self, where='', datatype=None, **kwargs):
        super().__init__(**kwargs)
        self.cop30_rurl = 'https://opentopography.s3.sdsc.edu/minio/raster/COP30/COP30_hh/'
        self.cop30_url = 'https://opentopography.s3.sdsc.edu/minio/download/raster/COP30/COP30_hh/'
        self.cop30_vrt_url = 'https://opentopography.s3.sdsc.edu/minio/download/raster/COP30/COP30_hh.vrt?token='
        self.cop_10_url = 'https://gisco-services.ec.europa.eu/dem/copernicus/outD/'
        self.cop_10_aux_url = 'https://gisco-services.ec.europa.eu/dem/copernicus/outA/'
        self.cop_10_web = 'https://ec.europa.eu/eurostat/web/gisco/geodata/reference-data/elevation/copernicus-dem/elevation'
        self.where = [where] if len(where) > 0 else []
        self.datatype = datatype        
        self._outdir = os.path.join(os.getcwd(), 'copernicus')
        self.name = 'copernicus'
        self.FRED = FRED.FRED(name=self.name, verbose=self.verbose)

        self.headers = { 'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:89.0) Gecko/20100101 Firefox/89.0',
                         'referer': 'https://opentopography.s3.sdsc.edu/minio/raster/COP30/COP30_hh/' }
        
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
        page = f_utils.Fetch(self.cop_10_url, verbose=True).fetch_html()
        rows = page.xpath('//a[contains(@href, ".zip")]/@href')
        if self.verbose:
            _prog = utils.CliProgress('scanning {} tiles in {}...'.format(len(rows), self.cop_10_url))
        
        for i, row in enumerate(rows):
            sid = row.split('.')[0]
            if self.verbose:
                _prog.update_perc((i, len(rows)))
                
            self.FRED._attribute_filter(["ID = '{}'".format(sid)])
            if self.FRED.layer is None or len(self.FRED.layer) == 0:
                spat = row.split('.')[0].split('_')[-1]
                x = int(spat.split('x')[-1])
                y = int(spat.split('x')[0].split('y')[-1])
                this_region = regions.Region().from_list(
                    [x, x + 10, y, y + 10]
                )
                geom = this_region.export_as_geom()
                if geom is not None:
                    surveys.append({'Name': row.split('.')[0], 'ID': sid, 'Agency': 'EU', 'Date': utils.this_date(),
                                    'MetadataLink': self.cop_10_aux_url, 'MetadataDate': utils.this_date(), 'DataLink': self.cop_10_url + row,
                                    'DataType': '3', 'DataSource': 'copernicus', 'HorizontalDatum': 'epsg:4326',
                                    'VerticalDatum': 'msl', 'Info': '', 'geom': geom})

        if self.verbose:
            _prog.end(0, 'scanned {} tiles in {}.'.format(len(rows), self.cop_10_url))
                    
        f = f_utils.Fetch(self.cop30_vrt_url, headers=self.headers, verbose=True)
        page = f.fetch_xml()
        fns = page.findall('.//SourceFilename')
        if self.verbose:
            _prog = utils.CliProgress('scanning {} tiles in {}...'.format(len(fns), self.cop30_url))
        
        for i, fn in enumerate(fns):
            sid = fn.text.split('/')[-1].split('.')[0]
            if self.verbose:
                _prog.update_perc((i, len(fns)))
                
            self.FRED._attribute_filter(["ID = '{}'".format(sid)])
            if self.FRED.layer is None or len(self.FRED.layer) == 0:            
                spat = fn.text.split('_10_')[-1].split('_DEM')[0]
                xsplit = '_E' if 'E' in spat else '_W'
                ysplit = 'S' if 'S' in spat else 'N'
                x = int(spat.split(xsplit)[-1].split('_')[0])
                y = int(spat.split(xsplit)[0].split(ysplit)[-1].split('_')[0])

                if xsplit == '_W':
                    x = x * -1
                if ysplit == 'S':
                    y = y * -1

                this_region = regions.Region().from_list([x, x + 1, y, y + 1])
                geom = this_region.export_as_geom()
                if geom is not None:
                    surveys.append({'Name': fn.text.split('.')[0].split('/')[-1], 'ID': sid, 'Agency': 'EU', 'Date': utils.this_date(),
                                    'MetadataLink': '', 'MetadataDate': utils.this_date(), 'DataLink': self.cop30_url + fn.text.split('/')[-1] + '?token=',
                                    'DataType': '1', 'DataSource': 'copernicus', 'HorizontalDatum': 'epsg:4326', 'Etcetra': self.cop30_rurl,
                                    'VerticalDatum': 'msl', 'Info': '', 'geom': geom})

        if self.verbose:
            _prog.end(0, 'scanned {} tiles in {}.'.format(len(fns), self.cop30_url))
            utils.echo_msg('added {} COPERNICUS DEM tiles'.format(len(surveys)))
            
        self.FRED._add_surveys(surveys)
        self.FRED._close_ds()

    def run(self):
        '''Run the COPERNICUS DEM fetching module'''

        if self.datatype is not None:
            self.where.append("DataType = '{}'".format(self.datatype))

        for surv in FRED._filter_FRED(self):
            for i in surv['DataLink'].split(','):
                self.results.append([i, i.split('/')[-1].split('?')[0], surv['DataType']])
                
        return(self)

    def yield_xyz(self, entry):
        """yield the xyz data from the copernicus fetch module"""
        
        if f_utils.Fetch(entry[0], callback=self.callback, verbose=self.verbose, headers=self.headers).fetch_file(entry[1]) == 0:
            src_cop_dems = utils.p_unzip(entry[1], ['tif'])
            for src_cop_dem in src_cop_dems:
                _ds = datasets.RasterFile(
                    fn=src_cop_dem,
                    data_format=200,
                    src_srs='epsg:4326',
                    dst_srs=self.dst_srs,
                    #name=src_cop_dem,
                    src_region=self.region,
                    verbose=self.verbose
                )
                for xyz in _ds.yield_xyz():
                    if xyz.z != 0:
                        yield(xyz)
                        
                utils.remove_glob(src_cop_dem, src_cop_dem + '.inf')
        utils.remove_glob(entry[1])

if __name__ == '__main__':
    cop_dem = CopernicusDEM()
    #cop_dem.update()
    
### End
