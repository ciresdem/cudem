### nos.py - NOS fetch
##
## Copyright (c) 2010 - 2022 Regents of the University of Colorado
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
##
## NOS Fetch
##
## fetch NOS BAG and XYZ sounding data from NOAA
## BAG data is in projected units and MLLW (height)
## XYZ data is CSV in MLLW (Sounding)
##
## NCEI maintains the National Ocean Service Hydrographic Data Base (NOSHDB) and Hydrographic 
## Survey Meta Data Base (HSMDB). Both are populated by the Office of Coast Survey and National 
## Geodetic Service, and provide coverage of coastal waters and the U.S. exclusive economic zone 
## and its territories.
##
## https://www.ngdc.noaa.gov/mgg/bathymetry/hydro.html
##
## Update for use of NOS Mapserver
##
## A map service showing the location and coverage of National Ocean Service (NOS) Hydrographic Surveys.
## The NOS Hydrographic Database (NOSHDB) and Hydrographic Survey Metadata Database (HSMDB), both maintained
## by NOS and NOAA's National Centers for Environmental Information (NCEI), provide extensive survey 
## coverage and ISO metadata of the coastal waters and Exclusive Economic Zone (EEZ) of the United States 
## and its territories. The NOSHDB contains digitized data from smooth sheets of hydrographic surveys 
## completed between 1837 and 1965, and from survey data acquired digitally on NOS survey vessels since 1965. 
## Data products from NOS surveys, including Bathymetric Attributed Grid (BAG) files, Descriptive Reports, 
## smooth sheet images, survey data images, textual gridded data, and geo-referenced sidescan sonar mosaics, 
## ISO metadata, and survey statistics are available for download from NCEI.
##
##Fields:
##

## SURVEY_ID ( type: esriFieldTypeString, alias: Survey ID, length: 10 )
## DATE_SURVEY_BEGIN ( type: esriFieldTypeDate, alias: Begin Date, length: 8 )
## DATE_SURVEY_END ( type: esriFieldTypeDate, alias: End Date, length: 8 )
## DATE_MODIFY_DATA ( type: esriFieldTypeDate, alias: Modify Data Date, length: 8 )
## DATE_SURVEY_APPROVAL ( type: esriFieldTypeDate, alias: Survey Approval Date, length: 8 )
## DATE_ADDED ( type: esriFieldTypeDate, alias: Date Added, length: 8 )
## SURVEY_YEAR ( type: esriFieldTypeDouble, alias: Survey Year )
## DIGITAL_DATA ( type: esriFieldTypeString, alias: Digital Data?, length: 15 )
## LOCALITY ( type: esriFieldTypeString, alias: Locality, length: 150 )
## SUBLOCALITY ( type: esriFieldTypeString, alias: Sublocality, length: 150 )
## PLATFORM ( type: esriFieldTypeString, alias: Platform Name, length: 150 )
## PRODUCT_ID ( type: esriFieldTypeString, alias: Product ID, length: 24 )
## BAGS_EXIST ( type: esriFieldTypeString, alias: BAGS_EXIST, length: 4 )
## DOWNLOAD_URL ( type: esriFieldTypeString, alias: Download URL, length: 256 )
## DECADE ( type: esriFieldTypeDouble, alias: Decade )
## PUBLISH ( type: esriFieldTypeString, alias: PUBLISH, length: 1 )
## OBJECTID ( type: esriFieldTypeOID, alias: OBJECTID )
## SHAPE ( type: esriFieldTypeGeometry, alias: SHAPE )
##
## Layer 0: Surveys with BAGs available (Bathymetric Attributed Grids).
## Layer 1: Surveys with digital sounding data available for download (including those with BAGs).
##
### Code:

import os
import json

from osgeo import osr
from osgeo import gdal

from cudem import utils
from cudem import regions
from cudem import datasets
from cudem import demfun

import cudem.fetches.utils as f_utils
import cudem.fetches.FRED as FRED

class HydroNOS(f_utils.FetchModule):
    """NOSHydro"""
    
    def __init__(self, where='1=1', layer=1, datatype=None, index=False, **kwargs):
        super().__init__(**kwargs)
        self._nos_dynamic_url = 'https://gis.ngdc.noaa.gov/arcgis/rest/services/web_mercator/nos_hydro_dynamic/MapServer'
        self._nos_url = 'https://gis.ngdc.noaa.gov/arcgis/rest/services/web_mercator/nos_hydro/MapServer'
        self._nos_data_url = 'https://data.ngdc.noaa.gov/platforms/ocean/nos/coast/'
        self._nos_query_url = '{0}/{1}/query?'.format(self._nos_dynamic_url, layer)
        self._outdir = os.path.join(os.getcwd(), 'hydronos')
        self.name = 'hydronos'
        self.where = where
        self.datatype = datatype
        self.index = index
        
    def run(self):
        if self.region is None:
            return([])

        _data = {
            'where': self.where,
            'outFields': '*',
            'geometry': self.region.format('bbox'),
            'inSR':4326,
            'outSR':4326,
            'f':'pjson',
            'returnGeometry':'False',
        }
        _req = f_utils.Fetch(self._nos_query_url, verbose=self.verbose).fetch_req(params=_data)
        if _req is not None:
            features = _req.json()
            for feature in features['features']:
                if self.index:
                    print(json.dumps(feature['attributes'], indent=4))
                else:
                    ID = feature['attributes']['SURVEY_ID']
                    link = feature['attributes']['DOWNLOAD_URL']
                    nos_dir = link.split('/')[-2]
                    data_link = '{}{}/{}/'.format(self._nos_data_url, nos_dir, ID)

                    if self.datatype is None or 'bag' in self.datatype.lower():
                        if feature['attributes']['BAGS_EXIST'] == 'TRUE':
                            page = f_utils.Fetch(data_link + 'BAG').fetch_html()
                            bags = page.xpath('//a[contains(@href, ".bag")]/@href')
                            [self.results.append(['{0}BAG/{1}'.format(data_link, bag), os.path.join(self._outdir, bag), 'bag']) for bag in bags]

                    if self.datatype is None or 'xyz' in self.datatype.lower():
                        page = f_utils.Fetch(data_link).fetch_html()
                        if page is not None:
                            geodas = page.xpath('//a[contains(@href, "GEODAS")]/@href')
                            if geodas:
                                xyz_link = data_link + 'GEODAS/{0}.xyz.gz'.format(ID)
                                self.results.append([xyz_link, os.path.join(self._outdir, xyz_link.split('/')[-1]), 'xyz'])                

    def yield_xyz(self, entry):
        src_nos = os.path.basename(entry[1])
        if 'ellipsoid' not in src_nos.lower():
            if f_utils.Fetch(entry[0], callback=self.callback, verbose=self.verbose).fetch_file(src_nos) == 0:
                if entry[2] == 'xyz':
                    nos_fns = utils.p_unzip(src_nos, ['xyz', 'dat'])
                    for src_xyz in nos_fns:
                        _ds = datasets.XYZFile(
                            fn=src_xyz,
                            data_format=168,
                            skip=1,
                            xpos=2,
                            ypos=1,
                            zpos=3,
                            z_scale=-1,
                            src_srs='epsg:4326+5866',
                            dst_srs=self.dst_srs,
                            src_region=self.region,
                            x_inc=self.x_inc,
                            y_inc=self.y_inc,
                            verbose=self.verbose,
                            remote=True
                        )
                        for xyz in _ds.yield_xyz():
                            yield(xyz)

                    utils.remove_glob(*nos_fns, *[x+'.inf' for x in nos_fns])

                elif entry[2] == 'bag':
                    src_bags = utils.p_unzip(src_nos, exts=['bag'])
                    for src_bag in src_bags:
                        if 'ellipsoid' not in src_bag.lower() and 'vb' not in src_bag.lower():
                            bag_ds = gdal.Open(src_bag)
                            bag_srs = bag_ds.GetSpatialRef()
                            bag_srs.AutoIdentifyEPSG()
                            bag_srs.SetAuthority('VERT_CS', 'EPSG', 5866)
                            bag_ds = None

                            _ds = datasets.BAGFile(
                                fn=src_bag,
                                data_format=201,
                                src_srs=bag_srs.ExportToWkt(),
                                dst_srs=self.dst_srs,
                                src_region=self.region,
                                x_inc=self.x_inc,
                                y_inc=self.y_inc,
                                verbose=self.verbose
                            )
                            for xyz in _ds.yield_xyz():
                                yield(xyz)

                    utils.remove_glob(*src_bags)
            utils.remove_glob(src_nos)
                        
## ==============================================
## the NOS class is the old NOS fetches module.
## This module scrapes the data from NOAA and generates
## a reference vector to discover dataset footprints. This is prone
## to possible error and can miss newer datasets if the reference
## vector is not up-to-date. Use HydroNOS class instead, which uses
## the NOAA MapServer to discover dataset footprints.
## ==============================================
class NOS(f_utils.FetchModule):
    """Fetch NOS BAG and XYZ sounding data from NOAA"""

    def __init__(self, where='', datatype=None, update=False, **kwargs):
        super().__init__(**kwargs)
        self._nos_url = 'https://www.ngdc.noaa.gov/mgg/bathymetry/hydro.html'
        self._nos_xml_url = lambda nd: 'https://data.noaa.gov/waf/NOAA/NESDIS/NGDC/MGG/NOS/%siso_u/xml/' %(nd)
        self._nos_iso_xml_url = lambda nd: 'https://data.noaa.gov/waf/NOAA/NESDIS/NGDC/MGG/NOS/%siso/xml/' %(nd)
        self._nos_directories = [
            "B00001-B02000/", "D00001-D02000/", "F00001-F02000/", \
            "H00001-H02000/", "H02001-H04000/", "H04001-H06000/", \
            "H06001-H08000/", "H08001-H10000/", "H10001-H12000/", \
            "H12001-H14000/", "L00001-L02000/", "L02001-L04000/", \
            "T00001-T02000/", "W00001-W02000/" \
        ]

        self._outdir = os.path.join(os.getcwd(), 'nos')
        self._nos_fmts = ['.xyz.gz', '.bag.gz', '.bag']
        self.where = [where] if len(where) > 0 else []
        self.datatype = datatype

        if self.datatype is not None:
            self.where.append("DataType LIKE '%{}%'".format(self.datatype.upper()))
        
        self.name = 'nos'
        self._urls = [self._nos_url]
        self.FRED = FRED.FRED(name=self.name, verbose=self.verbose)
        self.want_update = update
        
        if self.want_update:
            self.update()
        else:
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
            if self.callback():
                break
            
            surveys = []
            xml_catalog = self._nos_xml_url(nosdir)
            page = f_utils.Fetch(xml_catalog).fetch_html()
            if page is None:
                xml_catalog = self._nos_iso_xml_url(nosdir)
                page = f_utils.Fetch(xml_catalog).fetch_html()
                
            if page is None:
                utils.echo_error_msg('failed to retrieve {}'.format(nosdir))
                break
            
            rows = page.xpath('//a[contains(@href, ".xml")]/@href')
            if self.verbose:
                _prog = utils.CliProgress('scanning {} surveys in {}...'.format(len(rows), nosdir))

            for i, survey in enumerate(rows):
                if self.callback():
                    break
                
                sid = survey[:-4]
                if self.verbose:
                    _prog.update_perc((i, len(rows)))
                    
                self.FRED._attribute_filter(["ID = '{}'".format(sid)])
                if self.FRED.layer is None or len(self.FRED.layer) == 0:
                    this_xml = f_utils.iso_xml(xml_catalog + survey)
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
                            self.results.append([i, os.path.join(self._outdir, i.split('/')[-1]), surv['DataType']])
                    else:
                        self.results.append([i, os.path.join(self._outdir, i.split('/')[-1]), surv['DataType']])

    def yield_xyz(self, entry):
        src_nos = os.path.basename(entry[1])
        dt = None
        if f_utils.Fetch(entry[0], callback=self.callback, verbose=self.verbose).fetch_file(src_nos) == 0:
            dt = self._data_type(src_nos)
            
            if dt == 'geodas_xyz':
                nos_fns = utils.p_unzip(src_nos, ['xyz', 'dat'])
                for nos_f_r in nos_fns:
                    _ds = datasets.XYZFile(
                        fn=nos_f_r,
                        data_format=168,
                        skip=1,
                        xpos=2,
                        ypos=1,
                        zpos=3,
                        z_scale=-1,
                        #src_srs='epsg:4326+1089',
                        src_srs='epsg:4326',
                        dst_srs=self.dst_srs,
                        src_region=self.region,
                        verbose=self.verbose,
                        remote=True
                    )
                    for xyz in _ds.yield_xyz():
                        yield(xyz)
                        
                utils.remove_glob(*nos_fns, *[x+'.inf' for x in nos_fns])

            elif dt == 'grid_bag':
                src_bags = utils.p_unzip(src_nos, exts=['bag'])
                for src_bag in src_bags:
                    ## get bag proj from bag itself
                    _ds = datasets.RasterFile(
                        fn=src_bag,
                        data_format=200,
                        dst_srs=self.dst_srs,
                        #name=src_bag,
                        src_region=self.region,
                        verbose=self.verbose
                    )
                    for xyz in _ds.yield_xyz():
                        yield(xyz)
                        
                utils.remove_glob(*src_bags)
        utils.remove_glob(src_nos)
    
### End
