### emodnet.py - emodnet fetch - Europe
##
## Copyright (c) 2010 - 2021 CIRES Coastal DEM Team
##
## emodnet.py is part of CUDEM
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
## EMODNET Fetch
##
## fetch extracts of the EMOD DTM - Mostly European extents
## https://portal.emodnet-bathymetry.eu/
##
### Code:

import os
import lxml.etree

from cudem import utils
from cudem import regions
from cudem import datasets

import cudem.fetches.utils as f_utils
import cudem.fetches.FRED as FRED

class EMODNet(f_utils.FetchModule):
    """Fetch raster data from the EMODNET DTM"""

    def __init__(self, **kwargs):
        super().__init__(**kwargs) 
        self._emodnet_grid_url = 'https://ows.emodnet-bathymetry.eu/wcs?'
        self._outdir = os.path.join(os.getcwd(), 'emodnet')
        self.name = 'emodnet'

    def run(self):
        """Run the EMODNET fetching module"""
        
        if self.region is None: return([])

        desc_data = {
            'request': 'DescribeCoverage',
            'version': '2.0.1',
            'CoverageID': 'emodnet:mean',
            'service': 'WCS',
            }
        desc_req = f_utils.Fetch(self._emodnet_grid_url).fetch_req(params=desc_data)
        desc_results = lxml.etree.fromstring(desc_req.text.encode('utf-8'))
        g_env = desc_results.findall('.//{http://www.opengis.net/gml/3.2}GridEnvelope', namespaces = f_utils.namespaces)[0]
        hl = [float(x) for x in g_env.find('{http://www.opengis.net/gml/3.2}high').text.split()]

        g_bbox = desc_results.findall('.//{http://www.opengis.net/gml/3.2}Envelope')[0]
        lc = [float(x) for x in  g_bbox.find('{http://www.opengis.net/gml/3.2}lowerCorner').text.split()]
        uc = [float(x) for x in g_bbox.find('{http://www.opengis.net/gml/3.2}upperCorner').text.split()]
        
        ds_region = regions.Region().from_list([lc[1], uc[1], lc[0], uc[0]])
        resx = (uc[1] - lc[1]) / hl[0]
        resy = (uc[0] - lc[0]) / hl[1]

        if regions.regions_intersect_ogr_p(self.region, ds_region):
            emodnet_wcs = '{}service=WCS&request=GetCoverage&version=1.0.0&Identifier=emodnet:mean&coverage=emodnet:mean&format=GeoTIFF&bbox={}&resx={}&resy={}&crs=EPSG:4326'\
                                      .format(self._emodnet_grid_url, self.region.format('bbox'), resx, resy)
            outf = 'emodnet_{}.tif'.format(self.region.format('fn'))
            self.results.append([emodnet_wcs, outf, 'emodnet'])
        return(self)

    def yield_xyz(self, entry):
        src_emodnet = 'emodnet_tmp.tif'
        if f_utils.Fetch(entry[0], callback=self.callback, verbose=self.verbose).fetch_file(src_emodnet) == 0:
            _ds = datasets.RasterFile(fn=src_emodnet, data_format=200, epsg=4326, warp=self.warp,
                                      name=src_emodnet, src_region=self.region, verbose=self.verbose)
            for xyz in _ds.yield_xyz():
                yield(xyz)
                
        else:
            utils.echo_error_msg('failed to fetch remote file, {}...'.format(src_emodnet))
            
        utils.remove_glob(src_emodnet)

## ==============================================
## the EMODNetFRED class attempts to store EMDNet data in FRED
## rather than use the EMODNet API as above.
## ==============================================
class EMODNetFRED(f_utils.FetchModule):
    """Fetch raster data from the EMODNET DTM"""
    
    def __init__(self, where = [], **kwargs):
        super().__init__(**kwargs)
        self._emodnet_grid_url = 'https://ows.emodnet-bathymetry.eu/wcs?'
        self._emodnet_grid_cap = 'https://ows.emodnet-bathymetry.eu/wms?request=GetCapabilities&service=WMS'
        self._emodnet_help_url = 'https://portal.emodnet-bathymetry.eu/help/help.html'
        self._outdir = os.path.join(os.getcwd(), 'emodnet')
        self.where = where
        self.FRED = FRED.FRED(verbose = self.verbose)

        self.name = 'emodnet'
        self._info = 'European Bathymetry/Topographic data from EMODNET'
        self._title = '''EMODNET Elevation Data.'''
        self._usage = '''< emodnet >'''
        self._urls = [self._emodnet_help_url, self._emodnet_grid_url, self._emodnet_grid_cap]
        self.update_if_not_in_FRED()
        
    def update_if_not_in_FRED(self):
        self.FRED._open_ds()
        self.FRED._attribute_filter(["DataSource = '{}'".format(self._name)])
        if len(self.FRED.layer) == 0:
            self.FRED._close_ds()
            self.update()
        self.FRED._close_ds()

    def update(self):
        self.FRED._open_ds(1)
        emod_wcs = FRED.WCS(self._emodnet_grid_url)
        contents = emod_wcs._contents()
        if self.verbose:
            _prog = utils.CliProgress('Scanning {} WCS coverages from {}...'.format(len(contents), self._emodnet_grid_url))
        for i, layer in enumerate(contents):
            if self.verbose:
                _prog.update_perc((i, len(contents)))
            self.FRED._attribute_filter(["ID = '{}'".format(layer['CoverageId'][0])])
            if self.FRED.layer is None or len(self.FRED.layer) == 0:
                d = emod_wcs._describe_coverage(layer['CoverageId'][0])
                if d is not None:
                    ds_region = emod_wcs._get_coverage_region(d)
                    geom = ds_region.export_as_geom()
                    url = emod_wcs._get_coverage_url(layer['CoverageId'][0], region = ds_region)
                    self.FRED._add_survey(Name = d['name'][0], ID = layer['CoverageId'][0], Date = this_year(), MetadataLink = layer['Metadata'],
                                          MetadataDate = this_year(), DataLink = url, DataType = 'raster',
                                          DataSource = 'emodnet', HorizontalDatum = 4326, VerticalDatum = 1092,
                                          Info = layer['Abstract'], geom = geom)
        if self.verbose:
            _prog.end(0, 'Scanned {} WCS coverages from {}'.format(len(contents), self._emodnet_grid_url))
        self.FRED._close_ds()
        
    def run(self):        
        emod_wcs = FRED.WCS(self._emodnet_grid_url)
        for surv in FRED._filter_FRED(self):
            d = emod_wcs._describe_coverage(surv['ID'])
            if d is not None:
                ds_region = emod_wcs._get_coverage_region(d)
                if regions_intersect_ogr_p(self.region, ds_region):
                    emod_url = emod_wcs._get_coverage_url(surv['ID'], region=self.region)
                    outf = 'emodnet_{}.tif'.format(self.region.format('fn'))
                    self.results.append([emod_url, outf, surv['DataType']])

    def yield_xyz(self, entry):
        src_emodnet = 'emodnet_tmp.tif'
        if f_utils.Fetch(entry[0], callback=self.callback, verbose=self.verbose).fetch_file(src_emodnet) == 0:
            _ds = datasets.RasterFile(fn=src_emodnet, data_format=200, epsg=4326, warp=self.warp,
                                      name=src_emodnet, src_region=self.region, verbose=self.verbose)
            for xyz in _ds.yield_xyz():
                yield(xyz)
        else: utils.echo_error_msg('failed to fetch remote file, {}...'.format(src_emodnet))
        utils.remove_glob(src_emodnet)
        
### End
