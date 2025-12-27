### nasadem.py
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
##
## nasadem.py is part of CUDEM
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
## Fetch NASA Digital Elevation Model (NASADEM) data via OpenTopography.
##
### Code:

from typing import Optional
from cudem import utils
from cudem import regions
from cudem.fetches import fetches
from cudem.fetches import FRED

## ==============================================
## Constants
## ==============================================
NASADEM_BASE_URL = 'https://opentopography.s3.sdsc.edu/minio/raster/NASADEM/NASADEM_be/'
NASADEM_DOWNLOAD_URL = 'https://opentopography.s3.sdsc.edu/minio/download/raster/NASADEM/NASADEM_be/'
NASADEM_VRT_URL = 'https://opentopography.s3.sdsc.edu/minio/download/raster/NASADEM/NASADEM_be.vrt?token='

## ==============================================
## NASADEM Module
## ==============================================
class NASADEM(fetches.FetchModule):
    """NASA Digital Elevation Model (NASADEM)
    
    NASADEM is derived from SRTM processing improvements, elevation control, 
    and void-filling. This module fetches NASADEM via OpenTopography.
    
    https://www.earthdata.nasa.gov/esds/competitive-programs/measures/nasadem

    Configuration Example:
    < nasadem:datatype=None >
    """
    
    def __init__(self, where: str = '', datatype: Optional[str] = None, **kwargs):
        super().__init__(name='nasadem', **kwargs)
        self.where = [where] if where else []
        self.datatype = datatype
        
        ## Metadata defaults
        self.data_format = 200 # GDAL
        self.src_srs = 'epsg:4326+5773' # WGS84 + EGM96 

        ## Headers for OpenTopography
        self.headers = {
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:89.0) Gecko/20100101 Firefox/89.0',
            'referer': NASADEM_BASE_URL
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
        """Crawl the NASADEM VRT and update/generate the reference vector."""
        
        self.FRED._open_ds(1)
        
        try:
            surveys = []                    
            
            ## Fetch VRT which acts as the catalog
            f = fetches.Fetch(NASADEM_VRT_URL, headers=self.headers, verbose=self.verbose)
            page = f.fetch_xml()
            
            if page is None:
                utils.echo_error_msg("Failed to fetch NASADEM VRT catalog.")
                return

            fns = page.findall('.//SourceFilename')
            
            with utils.ccp(total=len(fns), desc='Scanning NASADEM datasets', leave=self.verbose) as pbar:        
                for fn in fns:
                    pbar.update()
                    filename = fn.text
                    sid = filename.split('/')[-1].split('.')[0]
                    
                    ## Check if already exists
                    self.FRED._attribute_filter([f"ID = '{sid}'"])
                    if self.FRED.layer is not None and len(self.FRED.layer) > 0:
                        continue

                    ## Parse spatial info from filename
                    ## Example format expected: NASADEM_HGT_n00e006.hgt
                    try:
                        spat = filename.split('_HGT_')[-1].split('.')[0]
                        
                        xsplit = 'e' if 'e' in spat else 'w'
                        ysplit = 's' if 's' in spat else 'n'
                        
                        ## Split string to get lat/lon integers
                        parts = spat.split(xsplit)
                        y_part = parts[0].split(ysplit)[-1]
                        x_part = parts[-1]
                        
                        x = int(x_part)
                        y = int(y_part)
                        
                        if xsplit == 'w':
                            x = -x
                        if ysplit == 's':
                            y = -y

                        ## Create 1x1 degree tile region
                        this_region = regions.Region().from_list([x, x + 1, y, y + 1])
                        geom = this_region.export_as_geom()
                        
                        if geom:
                            surveys.append({
                                'Name': sid,
                                'ID': sid,
                                'Agency': 'NASA',
                                'Date': utils.this_date(),
                                'MetadataLink': '',
                                'MetadataDate': utils.this_date(),
                                'DataLink': f"{NASADEM_DOWNLOAD_URL}{filename.split('/')[-1]}?token=",
                                'DataType': '1',
                                'DataSource': 'nasadem',
                                'HorizontalDatum': 4326,
                                'Etcetra': NASADEM_BASE_URL,
                                'VerticalDatum': 'msl',
                                'Info': 'NASADEM_be',
                                'geom': geom
                            })
                    except Exception as e:
                        if self.verbose:
                            utils.echo_warning_msg(f"Failed to parse NASADEM tile {filename}: {e}")

            if surveys:
                self.FRED._add_surveys(surveys)
        except Exception as e:
            utils.echo_error_msg(f"Error updating NASADEM FRED: {e}")                
        finally:
            self.FRED._close_ds()

            
    def run(self):
        """Run the NASADEM DEM fetching module."""
        
        _results = FRED._filter_FRED(self)
        
        for surv in _results:
            ## Handle comma-separated links if present
            for link in surv['DataLink'].split(','):
                ## Clean URL (remove query params for filename)
                filename = link.split('/')[-1].split('?')[0]
                self.add_entry_to_results(
                    link,
                    filename,
                    surv['DataType']
                )
                
        return self

    
### End
