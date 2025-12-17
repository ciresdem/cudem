### nasadem.py
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

## NASADEM
class NASADEM(fetches.FetchModule):
    """NASA Digital Elevation Model
    
    Our objective is to provide the scientific and civil communities with 
    a state-of-the-art global digital elevation model (DEM) derived from 
    a combination of Shuttle Radar Topography Mission (SRTM) processing 
    improvements, elevation control, void-filling and merging with data 
    unavailable at the time of the original SRTM production.

    https://www.earthdata.nasa.gov/esds/competitive-programs/measures/nasadem
    
    This module fetches NASADEM via OpenTopography. 
    You can also use the EarthData module to fetch this data.

    < nasadem:datatype=None >
    """
    
    def __init__(self, where='', datatype=None, **kwargs):
        super().__init__(name='nasadem', **kwargs)
        self.where = [where] if len(where) > 0 else []
        self.datatype = datatype
        
        ## various NASADEM urls in opentopography
        self.nasadem_rurl = 'https://opentopography.s3.sdsc.edu/minio/raster/NASADEM/NASADEM_be/'
        self.nasadem_url = ('https://opentopography.s3.sdsc.edu/minio/download/'
                            'raster/NASADEM/NASADEM_be/')
        self.nasadem_vrt_url = ('https://opentopography.s3.sdsc.edu/minio/download/'
                                'raster/NASADEM/NASADEM_be.vrt?token=')

        ## for dlim, data_format of 200 is a GDAL file.
        self.data_format = 200
        self.src_srs = 'epsg:4326+3855'

        ## NASADEM is in FRED, so set that up here.
        self.FRED = FRED.FRED(
            name=self.name, verbose=self.verbose
        )
        self.update_if_not_in_FRED()

        ## We need an opentopography referer here.
        self.headers = {
            'User-Agent': ('Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:89.0) '
                           'Gecko/20100101 Firefox/89.0'),
            'referer': 'https://opentopography.s3.sdsc.edu/minio/raster/NASADEM/NASADEM_be/'
        }

            
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
        NASADEM reference vector.
        """
        
        self.FRED._open_ds(1)
        surveys = []                    
        f = fetches.Fetch(self.nasadem_vrt_url, headers=self.headers, verbose=True)
        page = f.fetch_xml()
        fns = page.findall('.//SourceFilename')
        with tqdm(
                total=len(fns),
                desc='scanning NASADEM datasets',
                leave=self.verbose
        ) as pbar:        
            for i, fn in enumerate(fns):
                sid = fn.text.split('/')[-1].split('.')[0]
                pbar.update()
                self.FRED._attribute_filter(["ID = '{}'".format(sid)])
                if self.FRED.layer is None or len(self.FRED.layer) == 0:
                    spat = fn.text.split('_HGT_')[-1].split('.')[0]
                    xsplit = 'e' if 'e' in spat else 'w'
                    ysplit = 's' if 's' in spat else 'n'
                    x = int(spat.split(xsplit)[-1])
                    y = int(spat.split(xsplit)[0].split(ysplit)[-1])
                    if xsplit == 'w':
                        x = x * -1
                        
                    if ysplit == 's':
                        y = y * -1

                    this_region = regions.Region().from_list([x, x + 1, y, y + 1])
                    geom = this_region.export_as_geom()
                    if geom is not None:
                        surveys.append(
                            {
                                'Name': fn.text.split('.')[0].split('/')[-1],
                                'ID': sid,
                                'Agency': 'NASA',
                                'Date': utils.this_date(),
                                'MetadataLink': '',
                                'MetadataDate': utils.this_date(),
                                'DataLink': self.nasadem_url + fn.text.split('/')[-1] + '?token=',
                                'DataType': '1',
                                'DataSource': 'nasadem',
                                'HorizontalDatum': 4326,
                                'Etcetra': self.nasadem_rurl,
                                'VerticalDatum': 'msl',
                                'Info': '',
                                'geom': geom
                            }
                        )

        self.FRED._add_surveys(surveys)
        self.FRED._close_ds()

            
    def run(self):
        """Run the NASADEM DEM fetching module"""

        for surv in FRED._filter_FRED(self):
            for i in surv['DataLink'].split(','):
                self.add_entry_to_results(
                    i,
                    i.split('/')[-1].split('?')[0],
                    surv['DataType']
                )
                
        return(self)

### End
