### nos.py - SRTM fetch
##
## Copyright (c) 2010 - 2021 CIRES Coastal DEM Team
##
## nos.py is part of CUDEM
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
### Code:

import os
from cudem import utils
from cudem import regions
from cudem import datasets
import cudem.fetches.utils as f_utils
import cudem.fetches.FRED as FRED

## =============================================================================
##
## NOS Fetch
##
## fetch NOS BAG and XYZ sounding data from NOAA
## BAG data is in projected units and MLLW (height)
## XYZ data is CSV in MLLW (Sounding)
##
## https://www.ngdc.noaa.gov/mgg/bathymetry/hydro.html
##
## =============================================================================
class NOS(f_utils.FetchModule):
    """Fetch NOS BAG and XYZ sounding data from NOAA"""

    def __init__(self, where=[], datatype=None, **kwargs):
        super().__init__(**kwargs)
        self._nos_url = 'https://www.ngdc.noaa.gov/mgg/bathymetry/hydro.html'
        self._nos_xml_url = lambda nd: 'https://data.noaa.gov/waf/NOAA/NESDIS/NGDC/MGG/NOS/%siso_u/xml/' %(nd)
        self._nos_directories = [
            "B00001-B02000/", "D00001-D02000/", "F00001-F02000/", \
            "H00001-H02000/", "H02001-H04000/", "H04001-H06000/", \
            "H06001-H08000/", "H08001-H10000/", "H10001-H12000/", \
            "H12001-H14000/", "L00001-L02000/", "L02001-L04000/", \
            "T00001-T02000/", "W00001-W02000/" \
        ]

        self._outdir = os.path.join(os.getcwd(), 'nos')
        self._nos_fmts = ['.xyz.gz', '.bag.gz', '.bag']

        self.where = where
        self.datatype = datatype

        if self.datatype is not None:
            self.where.append("DataType LIKE '%{}%'".format(self.datatype.upper()))
        
        self.name = 'nos'
        self._info = '''Bathymetry surveys and data (xyz & BAG)'''
        self._title = '''NOAA NOS Bathymetric Data'''
        self._usage = '''< nos >'''
        self._urls = [self._nos_url]
        self.FRED = FRED.FRED(name=self.name,verbose=self.verbose)
        self.update_if_not_in_FRED()
        
    def update_if_not_in_FRED(self):
        self.FRED._open_ds()
        self.FRED._attribute_filter(["DataSource = '{}'".format(self.name)])
        if len(self.FRED.layer) == 0:
            self.FRED._close_ds()
            self.update()
        self.FRED._close_ds()
        
    def update(self):
        """Crawl the NOS database and update/generate the NOS reference vector."""
        
        self.FRED._open_ds(1)
        for nosdir in self._nos_directories:
            if self.callback(): break
            surveys = []
            xml_catalog = self._nos_xml_url(nosdir)
            page = f_utils.Fetch(xml_catalog).fetch_html()
            rows = page.xpath('//a[contains(@href, ".xml")]/@href')
            if self.verbose: _prog = utils.CliProgress('scanning {} surveys in {}...'.format(len(rows), nosdir))

            for i, survey in enumerate(rows):
                if self.callback(): break
                sid = survey[:-4]
                if self.verbose:
                    _prog.update_perc((i, len(rows)))
                self.FRED._attribute_filter(["ID = '{}'".format(sid)])
                if self.FRED.layer is None or len(self.FRED.layer) == 0:
                    this_xml = FRED.iso_xml(xml_catalog + survey)
                    h_epsg, v_epsg = this_xml.reference_system()
                    this_data = this_xml.data_links()
                    d_links = []
                    d_types = []

                    for key in this_data.keys():
                        if key in ['GEODAS_XYZ', 'BAG', 'GRID_BAG']:
                            d_links.append(this_data[key])
                            d_types.append(key)

                    geom = this_xml.bounds(geom=True)
                    if geom is not None:
                        surveys.append({'Name': this_xml.title(), 'ID': sid, 'Agency': 'NOAA/NOS', 'Date': this_xml.date(),
                                        'MetadataLink': this_xml.url, 'MetadataDate': this_xml.xml_date(), 'DataLink': ','.join([','.join(x) for x in d_links]),
                                        'DataType': ','.join(list(set(d_types))), 'DataSource': 'nos', 'HorizontalDatum': h_epsg,
                                        'VerticalDatum': v_epsg, 'Info': this_xml.abstract(), 'geom': geom})
            if self.verbose:
                _prog.end(0, 'scanned {} surveys in {}.'.format(len(rows), nosdir))
                utils.echo_msg('added {} surveys from {}'.format(len(surveys), nosdir))
            self.FRED._add_surveys(surveys)
        self.FRED._close_ds()

    def _data_type(self, src_nos):
        src_ext = os.path.basename(src_nos).split('.')
        if len(src_ext) > 2:
            if src_ext[-2] == 'bag': dt = 'grid_bag'
            elif src_ext[-2] == 'xyz': dt = 'geodas_xyz'
            else: dt = None
        elif len(src_ext) == 2:
            if src_ext[-1] == 'bag': dt = 'grid_bag'
            elif src_ext[-1] == 'xyz': dt = 'geodas_xyz'
            else: dt = None
        else:
            dt = None
            #print(src_nos)
            #print(dt)
        return(dt)
        
    def run(self):
        """Search the NOS reference vector and append the results
        to the results list."""
        
        for surv in FRED._filter_FRED(self):
            for i in surv['DataLink'].split(','):
                if i != '':
                    dt = self._data_type(i)
                    if self.datatype is not None:
                        if self.datatype.lower() in dt:
                            self.results.append([i, i.split('/')[-1], surv['DataType']])
                    else:
                        self.results.append([i, i.split('/')[-1], surv['DataType']])

    def yield_xyz(self, entry):
        src_nos = os.path.basename(entry[1])
        dt = None
        if f_utils.Fetch(entry[0], callback=self.callback, verbose=self.verbose).fetch_file(src_nos) == 0:
            dt = self._data_type(src_nos)
            
            if dt == 'geodas_xyz':
                nos_fns = utils.p_unzip(src_nos, ['xyz', 'dat'])
                for nos_f_r in nos_fns:
                    _ds = datasets.XYZFile(fn=nos_f_r, data_format=168, skip=1, xpos=2, ypos=1, zpos=3, z_scale=-1, epsg=4326, warp=self.warp,
                                           name=nos_f_r, src_region=self.region, verbose=self.verbose, remote=True)
                    for xyz in _ds.yield_xyz():
                        yield(xyz)
                utils.remove_glob(*nos_fns)

            elif dt == 'grid_bag':
                src_bags = utils.p_unzip(src_nos, exts=['bag'])
                for src_bag in src_bags:

                    _ds = datasets.RasterFile(fn=src_bag, data_format=200, warp=self.warp,
                                              name=src_bag, src_region=self.region, verbose=self.verbose)
                    for xyz in _ds.yield_xyz():
                        yield(xyz)
                utils.remove_glob(*src_bags)
        utils.remove_glob(src_nos)
    
### End
