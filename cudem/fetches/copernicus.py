### copernicus.py
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
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
## Fetch data from the Copernicus Digital Elevation Model (DEM).
##
### Code:

from typing import List, Dict, Optional, Any
from cudem import utils
from cudem import regions
from cudem.fetches import fetches
from cudem.fetches import FRED

## ==============================================
## Constants
## COP-30 (Global 30m) - Hosted by OpenTopography
## (requires API Key logic usually, or public bucket access)
##
## COP-10 (European 10m) - Hosted by Eurostat/GISCO
## ==============================================

COP30_BUCKET_URL = 'https://opentopography.s3.sdsc.edu/minio/raster/COP30/COP30_hh/'
COP30_DOWNLOAD_URL = 'https://opentopography.s3.sdsc.edu/minio/download/raster/COP30/COP30_hh/'
COP30_VRT_URL = 'https://opentopography.s3.sdsc.edu/minio/download/raster/COP30/COP30_hh.vrt?token='

COP10_URL = 'https://gisco-services.ec.europa.eu/dem/copernicus/outD/'
COP10_AUX_URL = 'https://gisco-services.ec.europa.eu/dem/copernicus/outA/'
COP10_INFO_URL = 'https://ec.europa.eu/eurostat/web/gisco/geodata/reference-data/elevation/copernicus-dem/elevation'

## Headers for scraping OpenTopography
HEADERS = {
    'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:89.0) Gecko/20100101 Firefox/89.0',
    'referer': COP30_BUCKET_URL
}

## ==============================================
## Copernicus Module
## ==============================================
class CopernicusDEM(fetches.FetchModule):
    """COPERNICUS satellite elevation data
    
    The Copernicus DEM is a Digital Surface Model (DSM) which 
    represents the surface of the Earth including buildings, 
    infrastructure and vegetation.

    Datatypes:
      '1' = COP-30 (Global ~30m)
      '3' = COP-10 (Europe ~10m)
    
    https://ec.europa.eu/eurostat/web/gisco/geodata/reference-data/elevation/copernicus-dem/elevation

    Configuration Example:
    < copernicus:datatype=None >
    """
    
    def __init__(self, where: str = '', datatype: Optional[str] = None, **kwargs):
        super().__init__(name='copernicus', **kwargs)
        self.where = [where] if where else []
        self.datatype = datatype        

        ## DLIM/Processing constants
        self.data_format = -2 # zipfile
        self.src_srs = 'epsg:4326+3855'
        self.headers = HEADERS

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

            
    def _update_cop10(self) -> List[Dict]:
        """Scrape and parse COP-10 (European 10m) datasets."""
        
        surveys = []
        page = fetches.Fetch(COP10_URL, verbose=self.verbose).fetch_html()
        if page is None:
            return surveys

        rows = page.xpath('//a[contains(@href, ".zip")]/@href')
        
        with utils.ccp(total=len(rows), desc='Scanning COP-10 datasets', leave=self.verbose) as pbar:
            for row in rows:
                pbar.update()
                sid = row.split('.')[0]
                
                ## Check if already exists in FRED to avoid re-parsing geometry
                self.FRED._attribute_filter([f"ID = '{sid}'"])
                if self.FRED.layer is not None and len(self.FRED.layer) > 0:
                    continue

                ## Parse Geometry from filename (e.g., ..._x30y40.zip)
                try:
                    spat = sid.split('_')[-1]
                    ## Expected format parts usually denote bottom-left corner
                    x_str = spat.split('x')[-1]
                    y_str = spat.split('x')[0].split('y')[-1]
                    x = int(x_str)
                    y = int(y_str)
                    
                    ## COP-10 tiles are typically 10x10 degrees.
                    this_region = regions.Region().from_list([x, x + 10, y, y + 10])
                    geom = this_region.export_as_geom()

                    if geom:
                        surveys.append({
                            'Name': sid,
                            'ID': sid,
                            'Agency': 'EU',
                            'Date': utils.this_date(),
                            'MetadataLink': COP10_AUX_URL,
                            'MetadataDate': utils.this_date(),
                            'DataLink': f"{COP10_URL}{row}",
                            'DataType': '3', # Mapped to '3' for COP-10
                            'DataSource': 'copernicus',
                            'HorizontalDatum': 'epsg:4326',
                            'VerticalDatum': 'msl',
                            'Info': 'COP-10',
                            'geom': geom
                        })
                except Exception as e:
                    if self.verbose:
                        utils.echo_warning_msg(f"Failed to parse COP-10 file {row}: {e}")
        return surveys

    
    def _update_cop30(self) -> List[Dict]:
        """Parse COP-30 (Global 30m) datasets from VRT."""
        
        surveys = []
        f = fetches.Fetch(COP30_VRT_URL, headers=self.headers, verbose=self.verbose)
        page = f.fetch_xml()
        
        if page is None:
            return surveys

        fns = page.findall('.//SourceFilename')
        
        with utils.ccp(total=len(fns), desc='Scanning COP-30 datasets', leave=self.verbose) as pbar:
            for fn in fns:
                pbar.update()
                
                ## Filename example: COP30_hh_10_N30_00_W120_00_DEM.tif
                raw_fn = fn.text
                sid = raw_fn.split('/')[-1].split('.')[0]

                # Check existence
                self.FRED._attribute_filter([f"ID = '{sid}'"])
                if self.FRED.layer is not None and len(self.FRED.layer) > 0:
                    continue

                try:
                    ## Extract spatial info
                    ## Format: ..._10_Nxx_00_Wxxx_00_DEM
                    spat = raw_fn.split('_10_')[-1].split('_DEM')[0]
                    
                    xsplit = '_E' if 'E' in spat else '_W'
                    ysplit = 'S' if 'S' in spat else 'N'
                    
                    parts = spat.split(xsplit)
                    y_part = parts[0].split(ysplit)[-1].split('_')[0]
                    x_part = parts[-1].split('_')[0]
                    
                    y = int(y_part)
                    x = int(x_part)

                    if xsplit == '_W':
                        x = x * -1
                    if ysplit == 'S':
                        y = y * -1

                    ## COP-30 tiles are 1x1 degree
                    this_region = regions.Region().from_list([x, x + 1, y, y + 1])
                    geom = this_region.export_as_geom()

                    if geom:
                        surveys.append({
                            'Name': sid,
                            'ID': sid,
                            'Agency': 'EU',
                            'Date': utils.this_date(),
                            'MetadataLink': '',
                            'MetadataDate': utils.this_date(),
                            'DataLink': f"{COP30_DOWNLOAD_URL}{raw_fn.split('/')[-1]}?token=",
                            'DataType': '1', # Mapped to '1' for COP-30
                            'DataSource': 'copernicus',
                            'HorizontalDatum': 'epsg:4326',
                            'VerticalDatum': 'msl',
                            'Etcetra': COP30_BUCKET_URL,
                            'Info': 'COP-30',
                            'geom': geom
                        })
                except Exception as e:
                     if self.verbose:
                        utils.echo_warning_msg(f"Failed to parse COP-30 file {raw_fn}: {e}")

        return surveys

    
    def update(self):
        """Crawl data sources and update the COPERNICUS reference vector."""
        
        self.FRED._open_ds(1)
        
        try:
            ## Update COP-10
            surveys_10 = self._update_cop10()
            if surveys_10:
                self.FRED._add_surveys(surveys_10)
            
            ## Update COP-30
            surveys_30 = self._update_cop30()
            if surveys_30:
                self.FRED._add_surveys(surveys_30)

        except Exception as e:
            utils.echo_error_msg(f"Error updating Copernicus FRED: {e}")
        finally:
            self.FRED._close_ds()

            
    def run(self):
        """Run the COPERNICUS DEM fetching module."""
        
        if self.datatype is not None:
            self.where.append(f"DataType = '{self.datatype}'")

        _results = FRED._filter_FRED(self)
        
        with utils.ccp(total=len(_results), desc='Scanning Copernicus datasets', leave=self.verbose) as pbar:
            for surv in _results:
                pbar.update()
                
                ## Handle comma-separated links if present
                for link in surv['DataLink'].split(','):
                    ## Clean URL (remove potential query params if needed or just filename)
                    clean_name = link.split('/')[-1].split('?')[0]
                    self.add_entry_to_results(
                        link,
                        clean_name,
                        surv['DataType']
                    )

        return self

    
### End
