### etopo.py
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

## ETOPO
class ETOPO(fetches.FetchModule):
    """Fetch ETOPO 2022 data. 

    The ETOPO Global Relief Model integrates topography, bathymetry, and 
    shoreline data from regional and global datasets to enable comprehensive, 
    high resolution renderings of geophysical characteristics of the earth’s 
    surface. 

    The model is designed to support tsunami forecasting, modeling, and warning, 
    as well as ocean circulation modeling and Earth visualization.  

    The current version, ETOPO 2022, is available in Ice Surface and Bedrock 
    versions that portray either the top layer of the ice sheets covering 
    Greenland and Antarctica, or the bedrock below. 

    For more information, email dem.info@noaa.gov

    We have bedrock or surface in both geotiff and netcdf. 
    Use `datatype` to specify which to fetch.

    datatype options are:
    'bed', 'bed_sid', 'surface', 'surface_sid', 'bed_netcdf', 
    'bed_sid_netcdf', 'surface_netcdf', 'surface_sid_netcdf'

    e.g.  datatype=surface_netcdf

    https://www.ncei.noaa.gov/products/etopo-global-relief-model
    
    < etopo:datatype=None >
    """
    
    def __init__(self, where='', datatype=None, **kwargs):
        super().__init__(name='etopo', **kwargs)        
        self.where = [where] if len(where) > 0 else []
        self.datatype = datatype
        
        ## The various etopo URLs.
        self.etopo_urls = {
            'netcdf': {
                'bed': 'https://www.ngdc.noaa.gov/thredds/fileServer/global/ETOPO2022/15s/15s_bed_elev_netcdf/',
                'bed_sid': 'https://www.ngdc.noaa.gov/thredds/fileServer/global/ETOPO2022/15s/15s_bed_sid_netcdf/',
                'surface': 'https://www.ngdc.noaa.gov/thredds/fileServer/global/ETOPO2022/15s/15s_surface_elev_netcdf/',
                'surface_sid': 'https://www.ngdc.noaa.gov/thredds/fileServer/global/ETOPO2022/15s/15s_surface_sid_netcdf/',
                },
            15: {
                'bed': 'https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO2022/data/15s/15s_bed_elev_gtif/',
                'bed_sid': 'https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO2022/data/15s/15s_bed_sid_gtif/',
                'surface': 'https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO2022/data/15s/15s_surface_elev_gtif/',
                'surface_sid': 'https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO2022/data/15s/15s_surface_sid_gtif/',
            },
            30: {
                'bed': 'https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO2022/data/30s/30s_bed_elev_gtif/ETOPO_2022_v1_30s_N90W180_bed.tif',
                'ice_elevation': 'https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO2022/data/30s/30s_surface_elev_gtif/ETOPO_2022_v1_30s_N90W180_surface.tif',
            },
            60: {
                'bed_elevation': 'https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO2022/data/60s/60s_bed_elev_gtif/ETOPO_2022_v1_60s_N90W180_bed.tif',
                'ice_elevation': 'https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO2022/data/60s/60s_surface_elev_gtif/ETOPO_2022_v1_60s_N90W180_surface.tif',
            },
        }
        self.etopo_aux_url = 'https://data.noaa.gov/metaview/page?xml=NOAA/NESDIS/NGDC/MGG/DEM//iso/xml/etopo_2022.xml&view=getDataView&header=none'

        ## for dlim, data_format is -2 for zip files.
        self.data_format = -2
        self.src_srs = 'epsg:4326+3855'

        ## Firefox on Windows here, 
        self.headers = {
            'User-Agent': ('Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:89.0) '
                           'Gecko/20100101 Firefox/89.0')
        }

        ## etopo is in FRED, so we set that up here.
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
            ["DataSource = '{}'".format(self.name)]
        )
        if len(self.FRED.layer) == 0:
            self.FRED._close_ds()
            self.update()
            
        self.FRED._close_ds()

        
    def update(self):
        """Crawl the ETOPO database and update/generate the 
        ETOPO reference vector in FRED.
        """
        
        self.FRED._open_ds(1)
        surveys = []
        res = 15
        for dtype in self.etopo_urls[res].keys():
            this_url = self.etopo_urls[res][dtype]
            netcdf_url = self.etopo_urls['netcdf'][dtype]
            page = fetches.Fetch(this_url, verbose=True).fetch_html()
            rows = page.xpath('//a[contains(@href, ".tif")]/@href')
            with tqdm(
                    desc=f'scanning for ETOPO {dtype} datasets',
                    leave=self.verbose
            ) as pbar:            
                for i, row in enumerate(rows):
                    pbar.update()
                    sid = row.split('.')[0]
                    self.FRED._attribute_filter(["ID = '{}'".format(sid)])
                    if self.FRED.layer is None or len(self.FRED.layer) == 0:
                        spat = row.split('.')[0].split('_{}'.format(dtype))[0].split('_')[-1]
                        xsplit = 'E' if 'E' in spat else 'W'
                        ysplit = 'S' if 'S' in spat else 'N'
                        x = int(spat.split(xsplit)[-1])
                        y = int(spat.split(xsplit)[0].split(ysplit)[-1])
                        if xsplit == 'W':
                            x = x * -1
                            
                        if ysplit == 'S':
                            y = y * -1

                        this_region = regions.Region().from_list(
                            [x, x + 15, y, y - 15]
                        )
                        geom = this_region.export_as_geom()
                        if geom is not None:
                            surveys.append(
                                {
                                    'Name': row.split('.')[0],
                                    'ID': sid,
                                    'Agency': 'NOAA',
                                    'Date': utils.this_date(),
                                    'MetadataLink': self.etopo_aux_url,
                                    'MetadataDate': utils.this_date(),
                                    'DataLink': this_url + row,
                                    'DataType': dtype,
                                    'DataSource': 'etopo',
                                    'HorizontalDatum': 'epsg:4326',
                                    'VerticalDatum': 'EGM2008',
                                    'Info': '',
                                    'geom': geom
                                }
                            )
                            # netcdf (thredds)
                            surveys.append(
                                {
                                    'Name': row.split('.')[0],
                                    'ID': sid,
                                    'Agency': 'NOAA',
                                    'Date': utils.this_date(),
                                    'MetadataLink': self.etopo_aux_url,
                                    'MetadataDate': utils.this_date(),
                                    'DataLink': netcdf_url + row.split('.')[0] + '.nc',
                                    'DataType': '{}_netcdf'.format(dtype),
                                    'DataSource': 'etopo',
                                    'HorizontalDatum': 'epsg:4326',
                                    'VerticalDatum': 'EGM2008',
                                    'Info': '',
                                    'geom': geom
                                }
                            )

        self.FRED._add_surveys(surveys)
        self.FRED._close_ds()

        
    def run(self):
        '''Run the ETOPO DEM fetching module'''

        if self.datatype is not None:
            self.where.append(
                f"DataType = '{self.datatype}'"
            )

        _results = FRED._filter_FRED(self)
        with tqdm(
                total=len(_results),
                desc='scanning ETOPO datasets',
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
