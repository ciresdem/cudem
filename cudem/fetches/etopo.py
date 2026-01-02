### etopo.py
##
## Copyright (c) 2022 - 2026 Regents of the University of Colorado
##
## etopo.py is part of CUDEM
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
### Commentary:
##
## Fetch ETOPO 2022 global relief data.
##
### Code:

from typing import List, Dict, Optional, Any
from cudem import utils
from cudem import regions
from cudem.fetches import fetches
from cudem.fetches import FRED

## ==============================================
## Constants
## ==============================================
ETOPO_BASE_URL_15S_GTIF = 'https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO2022/data/15s/'
ETOPO_BASE_URL_15S_NC = 'https://www.ngdc.noaa.gov/thredds/fileServer/global/ETOPO2022/15s/'
ETOPO_BASE_URL_30S_GTIF = 'https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO2022/data/30s/'
ETOPO_BASE_URL_60S_GTIF = 'https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO2022/data/60s/'
ETOPO_METADATA_URL = 'https://data.noaa.gov/metaview/page?xml=NOAA/NESDIS/NGDC/MGG/DEM//iso/xml/etopo_2022.xml&view=getDataView&header=none'

ETOPO_URLS = {
    'netcdf': {
        'bed': f'{ETOPO_BASE_URL_15S_NC}15s_bed_elev_netcdf/',
        'bed_sid': f'{ETOPO_BASE_URL_15S_NC}15s_bed_sid_netcdf/',
        'surface': f'{ETOPO_BASE_URL_15S_NC}15s_surface_elev_netcdf/',
        'surface_sid': f'{ETOPO_BASE_URL_15S_NC}15s_surface_sid_netcdf/',
    },
    15: {
        'bed': f'{ETOPO_BASE_URL_15S_GTIF}15s_bed_elev_gtif/',
        'bed_sid': f'{ETOPO_BASE_URL_15S_GTIF}15s_bed_sid_gtif/',
        'surface': f'{ETOPO_BASE_URL_15S_GTIF}15s_surface_elev_gtif/',
        'surface_sid': f'{ETOPO_BASE_URL_15S_GTIF}15s_surface_sid_gtif/',
    },
    30: {
        'bed': f'{ETOPO_BASE_URL_30S_GTIF}30s_bed_elev_gtif/ETOPO_2022_v1_30s_N90W180_bed.tif',
        'surface': '{ETOPO_BASE_URL_30S_GTIF}30s_surface_elev_gtif/ETOPO_2022_v1_30s_N90W180_surface.tif',
    },
    60: {
        'bed': '{ETOPO_BASE_URL_60S_GTIF}60s_bed_elev_gtif/ETOPO_2022_v1_60s_N90W180_bed.tif',
        'surface': '{ETOPO_BASE_URL_60S_GTIF}60s_surface_elev_gtif/ETOPO_2022_v1_60s_N90W180_surface.tif',
    }
}

## ==============================================
## ETOPO Module
## ==============================================
class ETOPO(fetches.FetchModule):
    """Fetch ETOPO 2022 data. 

    The ETOPO Global Relief Model integrates topography, bathymetry, and 
    shoreline data.

    Datatype options:
    'bed', 'bed_sid', 'surface', 'surface_sid', 'bed_netcdf', 
    'bed_sid_netcdf', 'surface_netcdf', 'surface_sid_netcdf'

    Configuration Example:
    < etopo:datatype=None >
    """
    
    def __init__(self, where: str = '', datatype: Optional[str] = None, **kwargs):
        super().__init__(name='etopo', **kwargs)        
        self.where = [where] if where else []
        self.datatype = datatype
        
        self.data_format = -2 # zip files
        self.src_srs = 'epsg:4326+3855' # WGS84 + EGM2008

        self.headers = {
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:89.0) Gecko/20100101 Firefox/89.0'
        }

        ## Initialize FRED
        self.FRED = FRED.FRED(name=self.name, verbose=self.verbose)
        self.update_if_not_in_FRED()

        
    def update_if_not_in_FRED(self):
        """Update the fetches module in FRED if it's not already in there."""
        
        self.FRED._open_ds()
        self.FRED._attribute_filter([f"DataSource = '{self.name}'"])
        
        if len(self.FRED.layer) == 0:
            self.FRED._close_ds()
            self.update()
        else:
            self.FRED._close_ds()

            
    def update(self):
        """Crawl the ETOPO database and update/generate the ETOPO reference vector in FRED."""
        
        self.FRED._open_ds(1)
        
        try:
            surveys = []
            res = 15 # Currently only crawling 15s tiles
            
            ## Iterate through the GeoTIFF directories
            for dtype, this_url in ETOPO_URLS[res].items():
                netcdf_url = ETOPO_URLS['netcdf'][dtype]
                
                ## Fetch directory listing
                page = fetches.Fetch(this_url, verbose=self.verbose).fetch_html()
                if page is None:
                    continue
                    
                rows = page.xpath('//a[contains(@href, ".tif")]/@href')
                
                with utils.ccp(total=len(rows), desc=f'Scanning ETOPO {dtype} datasets', leave=self.verbose) as pbar:            
                    for row in rows:
                        pbar.update()
                        filename = row.split('/')[-1]
                        sid = filename.split('.')[0]
                        
                        ## Check existence in FRED
                        self.FRED._attribute_filter([f"ID = '{sid}'"])
                        if self.FRED.layer is not None and len(self.FRED.layer) > 0:
                            continue

                        ## Parse Geometry from Filename
                        ## Format example: ETOPO_2022_v1_15s_N90W180_bed.tif
                        try:
                            ## Split out the coordinate string. 
                            parts = sid.split('_')
                            ## Look for the part matching N/S and E/W structure
                            ## e.g., N90W180
                            spat = next((p for p in parts if ('N' in p or 'S' in p) and ('E' in p or 'W' in p)), None)
                            
                            if not spat:
                                continue

                            xsplit = 'E' if 'E' in spat else 'W'
                            ysplit = 'S' if 'S' in spat else 'N'
                            
                            parts_geo = spat.split(xsplit)
                            y_str = parts_geo[0].split(ysplit)[-1]
                            x_str = parts_geo[-1]
                            
                            x = int(x_str)
                            y = int(y_str)

                            if xsplit == 'W':
                                x = -x
                            if ysplit == 'S':
                                y = -y

                            ## 15s tiles are 15x15 degrees. 
                            ## N90W180 is Top Left corner.
                            ## Region takes [xmin, xmax, ymin, ymax]
                            this_region = regions.Region().from_list(
                                [x, x + 15, y - 15, y]
                            )
                            
                            geom = this_region.export_as_geom()
                            
                            if geom:
                                ## Add GeoTIFF Entry
                                surveys.append({
                                    'Name': sid,
                                    'ID': sid,
                                    'Agency': 'NOAA',
                                    'Date': utils.this_date(),
                                    'MetadataLink': ETOPO_METADATA_URL,
                                    'MetadataDate': utils.this_date(),
                                    'DataLink': f"{this_url}{row}",
                                    'DataType': dtype,
                                    'DataSource': 'etopo',
                                    'HorizontalDatum': 'epsg:4326',
                                    'VerticalDatum': 'EGM2008',
                                    'Info': 'ETOPO 2022 GeoTIFF',
                                    'geom': geom
                                })
                                
                                ## Add NetCDF Entry
                                nc_filename = f"{sid}.nc"
                                surveys.append({
                                    'Name': sid,
                                    'ID': sid,
                                    'Agency': 'NOAA',
                                    'Date': utils.this_date(),
                                    'MetadataLink': ETOPO_METADATA_URL,
                                    'MetadataDate': utils.this_date(),
                                    'DataLink': f"{netcdf_url}{nc_filename}",
                                    'DataType': f"{dtype}_netcdf",
                                    'DataSource': 'etopo',
                                    'HorizontalDatum': 'epsg:4326',
                                    'VerticalDatum': 'EGM2008',
                                    'Info': 'ETOPO 2022 NetCDF',
                                    'geom': geom
                                })
                        except Exception as e:
                            if self.verbose:
                                utils.echo_warning_msg(f"Error parsing ETOPO tile {row}: {e}")

            self.FRED._add_surveys(surveys)
            
        except Exception as e:
            utils.echo_error_msg(f"Error updating ETOPO FRED: {e}")
        finally:
            self.FRED._close_ds()

            
    def run(self):
        """Run the ETOPO DEM fetching module."""

        if self.datatype is not None:
            self.where.append(f"DataType = '{self.datatype}'")

        _results = FRED._filter_FRED(self)
        
        with utils.ccp(total=len(_results), desc='Scanning ETOPO datasets', leave=self.verbose) as pbar:
            for surv in _results:
                pbar.update()
                ## Handle comma-separated links if present
                for link in surv['DataLink'].split(','):
                    self.add_entry_to_results(
                        link,
                        link.split('/')[-1].split('?')[0],
                        surv['DataType']
                    )
                
        return self

### End
