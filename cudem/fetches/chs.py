### chs.py - chs fetch - Canada
##
## Copyright (c) 2010 - 2021 CIRES Coastal DEM Team
##
## chs.py is part of CUDEM
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
##
## CHS Fetch
##
## fetch bathymetric soundings from the Canadian Hydrographic Service (CHS) - Canada Only
## https://open.canada.ca/data/en/dataset/d3881c4c-650d-4070-bf9b-1e00aabf0a1d
##
## NONNA 10 and NONNA 100
##
### Code:

import os
import lxml.etree

from cudem import utils
from cudem import regions
from cudem import datasets

import cudem.fetches.utils as f_utils
import cudem.fetches.FRED as FRED

class CHS(f_utils.FetchModule):
    """Fetch bathymetric soundings from the CHS"""
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._chs_api_url = "https://geoportal.gc.ca/arcgis/rest/services/FGP/CHS_NONNA_100/MapServer/0/query?"
        self._chs_url = 'https://data.chs-shc.ca/geoserver/wcs?'
        self._outdir = os.path.join(os.getcwd(), 'chs')
        self.name = 'chs'
        
    def run(self):
        """Run the CHS fetching module"""
        
        if self.region is None: return([])
        _data = {
            'request': 'DescribeCoverage',
            'version': '2.0.1',
            'CoverageID': 'caris:NONNA 100',
            'service': 'WCS',
            }
        _req = f_utils.Fetch(self._chs_url).fetch_req(params=_data)
        _results = lxml.etree.fromstring(_req.text.encode('utf-8'))
        
        g_env = _results.findall('.//{http://www.opengis.net/gml/3.2}GridEnvelope', namespaces=f_utils.namespaces)[0]
        hl = [float(x) for x in g_env.find('{http://www.opengis.net/gml/3.2}high').text.split()]

        g_bbox = _results.findall('.//{http://www.opengis.net/gml/3.2}Envelope')[0]
        lc = [float(x) for x in g_bbox.find('{http://www.opengis.net/gml/3.2}lowerCorner').text.split()]
        uc = [float(x) for x in g_bbox.find('{http://www.opengis.net/gml/3.2}upperCorner').text.split()]

        ds_region = regions.Region().from_list(
            [lc[1], uc[1], lc[0], uc[0]]
        )
        resx = (uc[1] - lc[1]) / hl[0]
        resy = (uc[0] - lc[0]) / hl[1]

        if regions.regions_intersect_ogr_p(self.region, ds_region):
            chs_wcs = '{}service=WCS&request=GetCoverage&version=1.0.0&Identifier=caris:NONNA+100&coverage=caris:NONNA+100&format=GeoTIFF&bbox={}&resx={}&resy={}&crs=EPSG:4326'\
                                  .format(self._chs_url, self.region.format('bbox'), resx, resy)
            outf = 'chs_{}.tif'.format(self.region.format('fn'))
            self.results.append([chs_wcs, outf, 'chs'])
        return(self)

    def yield_xyz(self, entry):
        src_chs = 'chs_tmp.tif'
        if f_utils.Fetch(entry[0], callback=self.callback, verbose=self.verbose).fetch_file(src_chs) == 0:
            _ds = datasets.RasterFile(
                fn=src_chs,
                data_format=200,
                src_srs='epsg:4326',
                dst_srs=self.dst_srs,
                name=src_chs,
                src_region=self.region,
                verbose=self.verbose
            )
            for xyz in _ds.yield_xyz():
                yield(xyz)
                
        else:
            utils.echo_error_msg('failed to fetch remote file, {}...'.format(src_chs))
            
        utils.remove_glob(src_chs)

## ==============================================
## class CHS_FRED attempts to add chs data to FRED.
## Use class CHS instead
## ==============================================
class CHS_FRED(f_utils.FetchModule):
    """Fetch raster data from CHS"""

    def __init__(self, where = [], **kwargs):
        super().__init__(**kwargs)
        self._chs_url = 'https://data.chs-shc.ca/geoserver/wcs?'
        self._chs_grid_cap = 'https://data.chs-shc.ca/geoserver/wcs?request=GetCapabilities&service=WMS'
        self._chs_info_url = 'https://open.canada.ca/data/en/dataset/d3881c4c-650d-4070-bf9b-1e00aabf0a1d'
        self._chs_api_url = "https://geoportal.gc.ca/arcgis/rest/services/FGP/CHS_NONNA_100/MapServer/0/query?"
        self._outdir = os.path.join(os.getcwd(), 'chs')
        self.where = where
        self.FRED = FRED.FRED(verbose = self.verbose)

        self.name = 'chs'
        self._info = '''CHS NONNA 10m and 100m Bathymetric Survey Grids; Non-Navigational gridded bathymetric data based on charts and soundings.'''
        self._title = '''Bathymetric data from CHS'''
        self._usage = '''< chs >'''
        self._urls = [self._chs_info_url, self._chs_url, self._chs_grid_cap]
        self._pdate_if_not_in_FRED()
        
    def update_if_not_in_FRED(self):
        self.FRED._open_ds()
        self.FRED._attribute_filter(["DataSource = '{}'".format(self.name)])
        if len(self.FRED.layer) == 0:
            self.FRED._close_ds()
            self.update()
        self.FRED._close_ds()
        
    def update(self):
        self.FRED._open_ds(1)
        chs_wcs = f_utils.WCS(self._chs_url)
        contents = chs_wcs._contents()
        if self.verbose:
            _prog = utils.CliProgress('Scanning {} WCS coverages from {}...'.format(len(contents), self._chs_url))
            
        for i, layer in enumerate(contents):
            if self.verbose:
                _prog.update_perc((i, len(contents)))
            if 'Tiles' not in layer['CoverageId'][0]:
                self.FRED._attribute_filter(["ID = '{}'".format(layer['CoverageId'][0])])
            if self.FRED.layer is None or len(self.FRED.layer) == 0:                
                d = chs_wcs._describe_coverage(layer['CoverageId'][0])
                if d is not None:
                    ds_region = chs_wcs._get_coverage_region(d)
                    geom = ds_region.export_as_geom()
                    url = chs_wcs._get_coverage_url(layer['CoverageId'][0], region=ds_region)
                    try:
                        name = d['name'][0]
                    except: name = d['CoverageId'][0]
                    try:
                        meta = layer['Metadata']
                    except: meta = None
                    try:
                        info = layer['Abstract']
                    except: info = None
                    self.FRED._add_survey(Name = name, ID = layer['CoverageId'][0], Date = this_year(), MetadataLink = meta,
                                          MetadataDate = this_year(), DataLink = url, DataType = 'raster',
                                          DataSource = 'chs', HorizontalDatum = 4326, VerticalDatum = 1092,
                                          Info = info, geom = geom)
        if self.verbose:
            _prog.end(0, 'Scanned {} WCS coverages from {}'.format(len(contents), self._chs_url))
        self.FRED._close_ds()
        
    def run(self):        
        chs_wcs = f_utils.WCS(self._chs_url)
        for surv in FRED._filter_FRED(self):
            d = chs_wcs._describe_coverage(surv['ID'])
            if d is not None:
                ds_region = chs_wcs._get_coverage_region(d)
                if regions_intersect_ogr_p(self.region, ds_region):
                    chs_url = chs_wcs._get_coverage_url(chs_wcs.fix_coverage_id(surv['ID']), region=self.region)
                    outf = '{}_{}.tif'.format(surv['ID'].replace(' ', '_').replace('caris__', 'chs_'), self.region.format('fn'))
                    self.results.append([chs_url, outf, surv['DataType']])

    def _yield_xyz(self, entry):
        src_chs = 'chs_tmp.tif'
        if f_utils.Fetch(entry[0], callback=self.callback, verbose=self.verbose).fetch_file(src_chs) == 0:
            _ds = datasets.RasterFile(fn=src_chs, data_format=200, src_srs='epsg:4326', dst_srs=self.dst_srs,
                                      name=src_chs, src_region=self.region, verbose=self.verbose)
            for xyz in _ds.yield_xyz():
                yield(xyz)
        else: utils.echo_error_msg('failed to fetch remote file, {}...'.format(src_chs))
        utils.remove_glob(src_chs)
        
### End
