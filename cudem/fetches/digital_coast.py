### digital_coast.py - NOAA Digital Coast fetch
##
## Copyright (c) 2010 - 2021 CIRES Coastal DEM Team
##
## digital_coast.py is part of CUDEM
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
from osgeo import ogr
from osgeo import gdal
from cudem import utils
from cudem import regions
from cudem import datasets
import cudem.fetches.utils as f_utils
import cudem.fetches.FRED as FRED

## =============================================================================
##
## Digital Coast ('dc')
## fetches-dc
## Raster and Lidar data from NOAAs Digital Coast
##
## FRED holds the dataset level, fetch the index shapefile to parse the individual
## data files...
##
## - check bounds in xml for slr dems
## =============================================================================
class DigitalCoast(f_utils.FetchModule):

    def __init__(self, where=[], datatype=None, inc=None, **kwargs):
        super().__init__(**kwargs)
        self._dc_url = 'https://coast.noaa.gov'
        self._dc_htdata_url = 'https://coast.noaa.gov/htdata/'
        self._dc_dirs = ['lidar1_z', 'lidar2_z', 'lidar3_z', 'lidar4_z', 'raster1', 'raster2', 'raster5']
        self._outdir = os.path.join(os.getcwd(), 'digital_coast')
        self.where = where
        self.datatype = datatype
        self.inc = utils.str2inc(inc)
        self.name = 'dc'
        self._info = '''Lidar and Raster data from NOAA's Digital Coast'''
        self._title = '''NOAA Digital Coast'''
        self._usage = '''< digital_coast >'''
        self._urls = [self._dc_url, self._dc_htdata_url]
        self.FRED = FRED.FRED(name=self.name, verbose = self.verbose)
        self._update_if_not_in_FRED()
        
    def _update_if_not_in_FRED(self):
        self.FRED._open_ds()
        self.FRED._attribute_filter(["DataSource = '{}'".format(self.name)])
        if len(self.FRED.layer) == 0:
            self.FRED._close_ds()
            self._update()
        self.FRED._close_ds()
        
    def _update(self):
        """Update the FRED reference vector after scanning
        the relevant metadata from Digital Coast.
        """
        
        #self.FRED = FRED(verbose=self.verbose, local=True)
        self.FRED._open_ds(1)
        for ld in self._dc_dirs:
            cols = []
            surveys = []
            page = f_utils.Fetch(self._dc_htdata_url + ld).fetch_html()
            if page is None: continue
            tr = page.xpath('//table')[0].xpath('.//tr')
            if len(tr) <= 0: continue
            [cols.append(i.text_content()) for i in tr[0]]
            
            if self.verbose:
                _prog = utils.CliProgress('scanning {} datasets in {}...'.format(len(tr), ld))
                
            for i in range(1, len(tr)):
                if self.callback(): break
                if self.verbose:
                    _prog.update_perc((i, len(tr))) #dc['ID #']))
                    
                cells = tr[i].getchildren()
                dc = {}
                for j, cell in enumerate(cells):
                    cl = cell.xpath('a')
                    if len(cl) > 0:
                        if cols[j] == 'Dataset Name':
                            dc[cols[j]] = cell.text_content()
                            dc['Metadata'] = cl[0].get('href')
                        else: dc[cols[j]] = cl[0].get('href')
                    else: dc[cols[j]] = cell.text_content()
                self.FRED._attribute_filter(["ID = '{}'".format(dc['ID #'])])
                if self.FRED.layer is None or len(self.FRED.layer) == 0:
                    if 'Metadata' in dc.keys():
                        this_xml = FRED.iso_xml(dc['Metadata'])
                        h_epsg, v_epsg = this_xml.reference_system()
                        geom = this_xml.bounds(geom=True)
                        if geom is not None:
                            if self.verbose:
                                _prog.update_perc((i, len(tr)), msg = '{} ** adding: {} **'.format(_prog.opm, dc['ID #']))
                                
                            surveys.append({'Name': dc['Dataset Name'], 'ID': dc['ID #'], 'Date': this_xml.date(),
                                            'MetadataLink': dc['Metadata'], 'MetadataDate': this_xml.xml_date(),
                                            'DataLink': dc['https'], 'IndexLink': dc['Tile Index'], 'Link': self._dc_url,
                                            'DataType': ld.split('_')[0], 'DataSource': 'dc', 'HorizontalDatum': h_epsg,
                                            'VerticalDatum': v_epsg, 'Info': this_xml.abstract(), 'geom': geom})
            self.FRED._add_surveys(surveys)
            if self.verbose:
                _prog.end(0, 'scanned {} datasets in {}.'.format(len(tr), ld))
                utils.echo_msg('added {} surveys from {}'.format(len(surveys), ld))
        self.FRED._close_ds()

    def run(self):
        
        if self.datatype is not None:
            self.where.append("DataType LIKE '%{}%'".format(self.datatype))
        
        for surv in FRED._filter_FRED(self):
            if self.callback(): break
            surv_shp_zip = os.path.basename(surv['IndexLink'])
            if f_utils.Fetch(surv['IndexLink'], callback=self.callback, verbose=self.verbose).fetch_file(surv_shp_zip) == 0:
                v_shps = utils.p_unzip(surv_shp_zip, ['shp', 'shx', 'dbf', 'prj'])
                v_shp = None
                for v in v_shps:
                    if v.split('.')[-1] == 'shp':
                        v_shp = v
                #try:
                v_ds = ogr.Open(v_shp)
                slay1 = v_ds.GetLayer(0)
                for sf1 in slay1:
                    geom = sf1.GetGeometryRef()
                    if geom.Intersects(self.region.export_as_geom()):
                        tile_url = sf1.GetField('URL').strip()
                        self.results.append([tile_url, '{}/{}'.format(surv['ID'], tile_url.split('/')[-1]), surv['DataType']])
                v_ds = slay1 = None
                #except: pass
                utils.remove_glob(surv_shp_zip, *v_shps)

    def yield_xyz(self, entry):
        src_dc = os.path.basename(entry[1])
        src_ext = src_dc.split('.')[-1].lower()
        if src_ext == 'laz' or src_ext == 'las': dt = 'lidar'
        elif src_ext == 'tif' or src_ext == 'img': dt = 'raster'
        else: dt = None
        if dt == 'lidar':
            if f_utils.Fetch(entry[0], callback=self.callback, verbose=self.verbose).fetch_file(src_dc) == 0:
                xyz_dat = utils.yield_cmd('las2txt -stdout -parse xyz -keep_xy {} -keep_class {} -i {}\
                '.format(region_format(self.region, 'te'), '2 29', src_dc), verbose = False)
                _ds = datasets.XYZFile(fn=xyz_dat, data_format=168, warp=self.warp,
                                       name=xyz_dat, src_region=self.region, verbose=self.verbose, remote=True)

                # if self.inc is not None:
                #     xyz_func = lambda p: _ds.dump_xyz(dst_port=p, encode=True)
                #     for xyz in utils.yield_cmd('gmt blockmedian -I{:.10f} {} -r -V'.format(self.inc, self.region.format('gmt')), verbose=self.verbose, data_fun=xyz_func):
                #         yield(xyzfun.XYZPoint().from_list([float(x) for x in xyz.split()]))
                # else:
                #     for xyz in _ds.yield_xyz():
                #         yield(xyz)

                y = _ds.block_xyz if self.inc is not None else _ds.yield_xyz

                for xyz in y():
                    yield(xyz)
                        
        elif dt == 'raster':
            try:
                src_ds = gdal.Open(entry[0])
                src_dc = entry[0]
            except Exception as e:
                f_utils.Fetch(entry[0], callback=self.callback, verbose=self.verbose).fetch_file(src_dc)
                try:
                    src_ds = gdal.Open(src_dc)
                except Exception as e:
                    utils.echo_error_msg('could not read dc raster file: {}, {}'.format(entry[0], e))
                    src_ds = None
            except Exception as e:
                utils.echo_error_msg('could not read dc raster file: {}, {}'.format(entry[0], e))
                src_ds = None
            
            if src_ds is not None:

                _ds = datasets.RasterFile(fn=src_dc, data_format=200, warp=self.warp,
                                          name=src_dc, src_region=self.region, verbose=self.verbose)
                _ds.src_ds = src_ds
                _ds.ds_open_p = True

                for xyz in _ds.block_xyz(self.inc) if self.inc is not None else _ds.yield_xyz():
                    yield(xyz)
            src_ds = None
        utils.remove_glob(src_dc)    
### End
