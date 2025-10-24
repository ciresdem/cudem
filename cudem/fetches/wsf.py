### wsf.py
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

## wsf
class WSF(fetches.FetchModule):
    """WSF from German Aerospace Service (DLR)

    World Settlement Footprint (WSF) 2019

    https://www.dlr.de/EN/Home/home_node.html
    https://geoservice.dlr.de/web/services
    
    < wsf >
    """

    def __init__(self, where = '', datatype = None, **kwargs):
        super().__init__(name='wsf', **kwargs)
        self.where = [where] if len(where) > 0 else []
        self.datatype = datatype

        ## The various WSF URLs
        self._wsf_url = 'https://download.geoservice.dlr.de/WSF2019/files/'

        ## WSF is in FRED, set that up here
        self.FRED = FRED.FRED(name=self.name, verbose=self.verbose)
        self.update_if_not_in_FRED()

        ## set user agent
        self.headers = {
            'User-Agent': ('Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:89.0) '
                           'Gecko/20100101 Firefox/89.0')
        }

        
    def update_if_not_in_FRED(self):
        self.FRED._open_ds()
        self.FRED._attribute_filter(["DataSource = '{}'".format(self.name)])
        if len(self.FRED.layer) == 0:
            self.FRED._close_ds()
            self.update()
            
        self.FRED._close_ds()

        
    def update(self):
        """Crawl the SWF database and update/generate the WSF 
        reference vector.
        """
        
        self.FRED._open_ds(1)
        surveys = []
        page = fetches.Fetch(self._wsf_url, verbose=True).fetch_html()
        rows = page.xpath('//a[contains(@href, ".tif")]/@href')
        with tqdm(
                total=len(rows),
                desc='scanning WSF datasets',
                leave=self.verbose
        ) as pbar:
            for i, row in enumerate(rows):
                pbar.update(1)
                sid = row.split('.')[0]
                if sid == 'WSF2019_cog':
                    continue

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
                        surveys.append(
                            {'Name': row.split('.')[0],
                             'ID': sid,
                             'Agency': 'DLR',
                             'Date': utils.this_date(),
                             'MetadataLink': row.split('.')[0] + '_stac.json',
                             'MetadataDate': utils.this_date(),
                             'DataLink': self._wsf_url + row,
                             'DataType': 'WSF',
                             'DataSource': 'WSF',
                             'HorizontalDatum': 'epsg:4326',
                             'VerticalDatum': 'None',
                             'Info': '',
                             'geom': geom,}
                        )

        self.FRED._add_surveys(surveys)
        self.FRED._close_ds()

        
    def run(self):
        """Run the WSF fetching module"""

        if self.datatype is not None:
            self.where.append("DataType = '{}'".format(self.datatype))

        for surv in FRED._filter_FRED(self):
            for i in surv['DataLink'].split(','):
                self.add_entry_to_results(
                    i, i.split('/')[-1].split('?')[0], surv['DataType']
                )
                
        return(self)

### End
