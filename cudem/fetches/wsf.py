### wsf.py
##
## Copyright (c) 2010 - 2026 Regents of the University of Colorado
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
## Fetch World Settlement Footprint (WSF) 2019 data from DLR.
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
WSF_BASE_URL = 'https://download.geoservice.dlr.de/WSF2019/files/'

## ==============================================
## WSF Module
## ==============================================
class WSF(fetches.FetchModule):
    """WSF from German Aerospace Service (DLR)

    World Settlement Footprint (WSF) 2019.
    Data is organized in 2x2 degree tiles.

    https://www.dlr.de/EN/Home/home_node.html
    https://geoservice.dlr.de/web/services

    Configuration Example:    
    < wsf >
    """

    def __init__(self, where: str = '', datatype: Optional[str] = None, **kwargs):
        super().__init__(name='wsf', **kwargs)
        self.where = [where] if where else []
        self.datatype = datatype

        self.headers = {
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:89.0) Gecko/20100101 Firefox/89.0'
        }

        ## Initialize FRED
        self.FRED = FRED.FRED(name=self.name, verbose=self.verbose)
        self.update_if_not_in_FRED()

        
    def update_if_not_in_FRED(self):
        """Update the fetches module in FRED if it's not already in there."""
        
        self.FRED._open_ds()
        try:
            self.FRED._attribute_filter([f"DataSource = '{self.name}'"])
            if len(self.FRED.layer) == 0:
                self.FRED._close_ds()
                self.update()
        finally:
            if self.FRED.ds is not None:
                self.FRED._close_ds()

                
    def update(self):
        """Crawl the WSF directory and update/generate the WSF reference vector."""
        
        self.FRED._open_ds(1)
        
        try:
            surveys = []
            page = fetches.Fetch(WSF_BASE_URL, verbose=self.verbose).fetch_html()
            
            if page is None:
                return

            rows = page.xpath('//a[contains(@href, ".tif")]/@href')
            
            with utils.ccp(total=len(rows), desc='Scanning WSF datasets', leave=self.verbose) as pbar:
                for row in rows:
                    pbar.update(1)
                    sid = row.split('.')[0]
                    
                    ## Skip COG overview files if present
                    if sid == 'WSF2019_cog':
                        continue

                    ## Check existence in FRED
                    self.FRED._attribute_filter([f"ID = '{sid}'"])
                    if self.FRED.layer is not None and len(self.FRED.layer) > 0:
                        continue
                    
                    ## Parse spatial info from filename
                    ## Format expected: WSF2019_v1_{minx}_{miny}.tif
                    ## e.g., WSF2019_v1_-74_40.tif
                    try:
                        spat = sid.split('_')
                        x = int(spat[-2])
                        y = int(spat[-1])
                        
                        ## WSF 2019 tiles are 2x2 degrees
                        this_region = regions.Region().from_list(
                            [x, x + 2, y, y + 2]
                        )
                        geom = this_region.export_as_geom()
                        
                        if geom is not None:
                            surveys.append({
                                'Name': sid,
                                'ID': sid,
                                'Agency': 'DLR',
                                'Date': utils.this_date(),
                                'MetadataLink': f"{sid}_stac.json",
                                'MetadataDate': utils.this_date(),
                                'DataLink': f"{WSF_BASE_URL}{row}",
                                'DataType': 'WSF',
                                'DataSource': 'WSF',
                                'HorizontalDatum': 'epsg:4326',
                                'VerticalDatum': 'None',
                                'Info': 'World Settlement Footprint 2019',
                                'geom': geom
                            })
                    except (IndexError, ValueError):
                        if self.verbose:
                            utils.echo_warning_msg(f"Could not parse spatial info for {row}")

            self.FRED._add_surveys(surveys)
        except Exception as e:
            utils.echo_error_msg(f"Error updating WSF FRED: {e}")            
        finally:
            self.FRED._close_ds()

            
    def run(self):
        """Run the WSF fetching module."""

        if self.datatype is not None:
            self.where.append(f"DataType = '{self.datatype}'")

        _results = FRED._filter_FRED(self)
        
        for surv in _results:
            ## Handle comma-separated links if present
            for link in surv['DataLink'].split(','):
                filename = link.split('/')[-1].split('?')[0]
                self.add_entry_to_results(
                    link, 
                    filename, 
                    surv['DataType']
                )
                
        return self
    
### End
