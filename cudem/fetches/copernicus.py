### copernicus.py
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
from cudem import utils
from cudem import regions
from cudem.fetches import fetches
from cudem.fetches import FRED

## Copernicus
class CopernicusDEM(fetches.FetchModule):
    """COPERNICUS sattelite elevation data
    
    The Copernicus DEM is a Digital Surface Model (DSM) which 
    represents the surface of the Earth including buildings, 
    infrastructure and vegetation.

    datatype of 1 is 10 m and datatype of 3 is 30 m
    
    https://ec.europa.eu/eurostat/web/gisco/geodata/reference-data/elevation/copernicus-dem/elevation

    < copernicus:datatype=None >
    """
    
    def __init__(self, where='', datatype=None, **kwargs):
        super().__init__(name='copernicus', **kwargs)
        self.where = [where] if len(where) > 0 else []
        self.datatype = datatype        

        ## various Copernicus URLs
        self.cop30_rurl = 'https://opentopography.s3.sdsc.edu/minio/raster/COP30/COP30_hh/'
        self.cop30_url = 'https://opentopography.s3.sdsc.edu/minio/download/raster/COP30/COP30_hh/'
        self.cop30_vrt_url = ('https://opentopography.s3.sdsc.edu/minio/download/'
                              'raster/COP30/COP30_hh.vrt?token=')
        self.cop_10_url = 'https://gisco-services.ec.europa.eu/dem/copernicus/outD/'
        self.cop_10_aux_url = 'https://gisco-services.ec.europa.eu/dem/copernicus/outA/'
        self.cop_10_web = ('https://ec.europa.eu/eurostat/web/gisco/geodata/'
                           'reference-data/elevation/copernicus-dem/elevation')

        ## for dlim, data_format of -2 is zipfile
        self.data_format = -2
        self.src_srs = 'epsg:4326+3855'

        ## Firefox on Windows here, referer from opentopography
        self.headers = {
            'User-Agent': ('Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:89.0) '
                           'Gecko/20100101 Firefox/89.0'),
            'referer': 'https://opentopography.s3.sdsc.edu/minio/raster/COP30/COP30_hh/'
        }
        ## Copernicus is in FRED, so set that up here.
        self.FRED = FRED.FRED(
            name=self.name, verbose=self.verbose
        )
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
        """Crawl the COP30 database and update/generate the 
        COPERNICUS reference vector."""
        
        self.FRED._open_ds(1)
        surveys = []
        page = fetches.Fetch(self.cop_10_url, verbose=True).fetch_html()
        rows = page.xpath('//a[contains(@href, ".zip")]/@href')
        with tqdm(
                desc='scanning for COPERNICUS COP-10 datasets',
                leave=self.verbose
        ) as pbar:            
            for i, row in enumerate(rows):
                pbar.update()
                sid = row.split('.')[0]
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
                        surveys.append(
                            {
                                'Name': row.split('.')[0],
                                'ID': sid,
                                'Agency': 'EU',
                                'Date': utils.this_date(),
                                'MetadataLink': self.cop_10_aux_url,
                                'MetadataDate': utils.this_date(),
                                'DataLink': self.cop_10_url + row,
                                'DataType': '3',
                                'DataSource': 'copernicus',
                                'HorizontalDatum': 'epsg:4326',
                                'VerticalDatum': 'msl',
                                'Info': '',
                                'geom': geom
                            }
                        )

        f = fetches.Fetch(self.cop30_vrt_url, headers=self.headers, verbose=True)
        page = f.fetch_xml()
        fns = page.findall('.//SourceFilename')
        with tqdm(
                total=len(fns),
                desc='scanning for COPERNICUS COP-30 datasets',
                leave=self.verbose
        ) as pbar:
            for i, fn in enumerate(fns):
                pbar.update()
                sid = fn.text.split('/')[-1].split('.')[0]
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
                        surveys.append(
                            {
                                'Name': fn.text.split('.')[0].split('/')[-1],
                                'ID': sid,
                                'Agency': 'EU',
                                'Date': utils.this_date(),
                                'MetadataLink': '',
                                'MetadataDate': utils.this_date(),
                                'DataLink': self.cop30_url + fn.text.split('/')[-1] + '?token=',
                                'DataType': '1',
                                'DataSource': 'copernicus',
                                'HorizontalDatum': 'epsg:4326',
                                'Etcetra': self.cop30_rurl,
                                'VerticalDatum': 'msl',
                                'Info': '',
                                'geom': geom
                            }
                        )

        self.FRED._add_surveys(surveys)
        self.FRED._close_ds()

        
    def run(self):
        """Run the COPERNICUS DEM fetching module"""
        
        if self.datatype is not None:
            self.where.append(
                f"DataType = '{self.datatype}'"
            )

        _results = FRED._filter_FRED(self)
        with tqdm(
                total=len(_results),
                desc='scanning COPERNICUS datasets',
                leave=self.verbose
        ) as pbar:
            for surv in _results:
                pbar.update()
                for i in surv['DataLink'].split(','):
                    self.add_entry_to_results(
                        i,
                        i.split('/')[-1].split('?')[0],
                        surv['DataType']
                    )

        return(self)

### End
