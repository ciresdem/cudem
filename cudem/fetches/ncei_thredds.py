### ncei_thredds.py - NOAA NCEI THREDDS DEM fetch
##
## Copyright (c) 2010 - 2021 CIRES Coastal DEM Team
##
## ncei_thredds.py is part of CUDEM
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
## NCEI THREDDS Catalog (CUDEM)
##
## =============================================================================
class NCEIThreddsCatalog(f_utils.FetchModule):
    """Fetch DEMs from NCEI THREDDS Catalog"""

    def __init__(self, where=[], **kwargs):
        super().__init__(**kwargs)
        self._nt_catalog = "https://www.ngdc.noaa.gov/thredds/demCatalog.xml"
        self._ngdc_url = "https://www.ngdc.noaa.gov"
        self._outdir = os.path.join(os.getcwd(), 'ncei_thredds')        
        self.where = where

        self.name = 'ncei_thredds'
        self._info = '''NOAA NCEI THREDDS DEM Catalog Access.
Digital Elevation Models around the world at various resolutions and extents.
NCEI builds and distributes high-resolution, coastal digital elevation models (DEMs) that integrate ocean 
bathymetry and land topography supporting NOAA's mission to understand and predict changes in Earth's environment, 
and conserve and manage coastal and marine resources to meet our Nation's economic, social, and environmental needs.
\nDEMs are used for coastal process modeling (tsunami inundation, storm surge, sea-level rise, contaminant dispersal, 
etc.), ecosystems management and habitat research, coastal and marine spatial planning, and hazard mitigation and 
community preparedness.'''
        self._title = '''NCEI THREDDS'''
        self._usage = '''< ncei_thredds >'''
        self._urls = [self._nt_catalog, self._ngdc_url]
        self.FRED = FRED.FRED(name=self.name, verbose = self.verbose)
        self.update_if_not_in_FRED()
        
    def update_if_not_in_FRED(self):
        self.FRED._open_ds()
        self.FRED._attribute_filter(["DataSource = '{}'".format(self.name)])
        if len(self.FRED.layer) == 0:
            self.FRED._close_ds()
            self.update()
        self.FRED._close_ds()
        
    def _parse_catalog(self, catalog_url):
        ntCatalog = FRED.iso_xml(catalog_url)
        ntCatRefs = ntCatalog.xml_doc.findall('.//th:catalogRef', namespaces = ntCatalog.namespaces)
        for ntCatRef in ntCatRefs:
            ntCatHref =ntCatRef.attrib['{http://www.w3.org/1999/xlink}href']
            if ntCatHref[0] == "/":
                ntCatUrl = '{}{}'.format(self._ngdc_url, ntCatHref)
            else: ntCatUrl = '{}/{}'.format(os.path.dirname(catalog_url), ntCatHref)
            self._parse_dataset(ntCatUrl)
            
    def _parse_dataset(self, catalog_url):
        ntCatXml = FRED.iso_xml(catalog_url)
        this_ds = ntCatXml.xml_doc.findall('.//th:dataset', namespaces = ntCatXml.namespaces)
        this_ds_services = ntCatXml.xml_doc.findall('.//th:service', namespaces = ntCatXml.namespaces)
        if self.verbose:
            _prog = utils.CliProgress('scanning {} datasets in {}...'.format(len(this_ds), this_ds[0].attrib['name']))
        surveys = []
        
        for i, node in enumerate(this_ds):
            this_title = node.attrib['name']
            this_id = node.attrib['ID']
            if self.verbose:
                _prog.update_perc((i, len(this_ds)))
                
            self.FRED._attribute_filter(["ID = '{}'".format(this_id)])
            if self.FRED.layer is None or len(self.FRED.layer) == 0:
                subCatRefs = node.findall('.//th:catalogRef', namespaces = ntCatXml.namespaces)
                if len(subCatRefs) > 0:
                    self._parse_catalog(catalog_url)
                    break
                try:
                    ds_path = node.attrib['urlPath']
                except: continue

                iso_url = False
                wcs_url = False
                http_url = False
                for service in this_ds_services:
                    service_name = service.attrib['name']
                    if service_name == 'iso': iso_url = '{}{}{}'.format(self._ngdc_url, service.attrib['base'], ds_path)
                    if service_name == 'wcs': wcs_url = '{}{}{}'.format(self._ngdc_url, service.attrib['base'], ds_path)
                    if service_name == 'http': http_url = '{}{}{}'.format(self._ngdc_url, service.attrib['base'], ds_path)

                this_xml = FRED.iso_xml(iso_url)
                title = this_xml.title()
                h_epsg, v_epsg = this_xml.reference_system()

                zv = this_xml.xml_doc.findall('.//gmd:dimension/gmd:MD_Band/gmd:sequenceIdentifier/gco:MemberName/gco:aName/gco:CharacterString', namespaces = this_xml.namespaces)
                if zv is not None:
                    for zvs in zv:
                        print(zvs.text)
                        if zvs.text == 'bathy' or zvs.text == 'Band1' or zvs.text == 'z':
                            zvar = zvs.text
                            break
                        else:
                            zvar = 'z'
                geom = this_xml.bounds(geom=True)
                if geom is not None:
                    surveys.append({'Name': title, 'ID': this_id, 'Agency': 'NOAA', 'Date': this_xml.date(),
                                    'MetadataLink': this_xml.url, 'MetadataDate': this_xml.xml_date(),
                                    'DataLink': http_url, 'IndexLink': wcs_url, 'Link': self._nt_catalog,
                                    'DataType': 'raster', 'DataSource': 'ncei_thredds', 'HorizontalDatum': h_epsg,
                                    'VerticalDatum': v_epsg, 'Etcetra': zvar, 'Info': this_xml.abstract(), 'geom': geom})
        self.FRED._add_surveys(surveys) 
        if self.verbose:
            _prog.end(0, 'scanned {} datasets in {}.'.format(len(this_ds), this_ds[0].attrib['name']))
            utils.echo_msg('added {} surveys from {}'.format(len(surveys), this_ds[0].attrib['name']))
        
    def update(self):
        self.FRED._open_ds(1)
        self._parse_catalog(self._nt_catalog)
        self.FRED._close_ds()
    
    def run(self):
        """Search for data in the reference vector file"""
        
        for surv in FRED._filter_FRED(self):
            #wcs_url = "{}?request=GetCoverage&version=1.0.0&service=WCS&coverage={}&bbox={}&format=NetCDF3"\
            #    .format(surv['IndexLink'], surv['Etcetra'], self.region.format('bbox'))
            #print(wcs_url)
            #self.results.append([wcs_url, surv['DataLink'].split(',')[0].split('/')[-1], surv['DataType']])
            for d in surv['DataLink'].split(','):
                if d != '':
                    self.results.append([d, d.split('/')[-1], surv['DataType']])

    def yield_xyz(self, entry):
        src_ncei = os.path.basename(entry[1])
        #try:
        #    src_ds = gdal.Open(entry[0])
        #    src_dc = entry[0]
        #except Exception as e:
        f_utils.Fetch(entry[0], callback=self.callback, verbose=self.verbose).fetch_file(src_ncei)
        try:
            src_ds = gdal.Open(src_ncei)
        except Exception as e:
            utils.echo_error_msg('could not read ncei raster file: {}, {}'.format(entry[0], e))
            src_ds = None
        #except Exception as e:
        #    utils.echo_error_msg('could not read ncei raster file: {}, {}'.format(entry[0], e))
        #    src_ds = None
            
        if src_ds is not None:

            _ds = datasets.RasterFile(fn=src_ncei, data_format=200, epsg=4326, warp=self.warp,
                                      name=src_ncei, src_region=self.region, verbose=self.verbose)
            _ds.src_ds = src_ds
            _ds.ds_open_p = True
            for xyz in _ds.yield_xyz():
                yield(xyz)
        src_ds = None
        utils.remove_glob(src_ncei)    
                    
### End
