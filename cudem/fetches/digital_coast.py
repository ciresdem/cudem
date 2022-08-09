### digital_coast.py - NOAA Digital Coast fetch
##
## Copyright (c) 2010 - 2022 Regents of the University of Colorado
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
##
## Digital Coast ('dc')
## fetches-dc
## Raster and Lidar data from NOAAs Digital Coast
##
## FRED holds the dataset level, fetch the index shapefile to parse the individual
## data files...
##
## - check bounds in xml for slr dems
##
## update to DAV to remove dependence on FRED.
##
## This map service presents spatial information about Elevation Data Access Viewer services across the United States
## and Territories in the Web Mercator projection. The service was developed by the National Oceanic and Atmospheric
## Administration (NOAA), but may contain data and information from a variety of data sources, including non-NOAA data.
## NOAA provides the information “as-is” and shall incur no responsibility or liability as to the completeness or accuracy
## of this information. NOAA assumes no responsibility arising from the use of this information. The NOAA Office for Coastal
## Management will make every effort to provide continual access to this service but it may need to be taken down during
## routine IT maintenance or in case of an emergency. If you plan to ingest this service into your own application and would
## like to be informed about planned and unplanned service outages or changes to existing services, please register for our
## Data Services Newsletter (http://coast.noaa.gov/digitalcoast/publications/subscribe). For additional information, please
## contact the NOAA Office for Coastal Management (coastal.info@noaa.gov).
##
##
## Fields:
##
##     OBJECTID ( type: esriFieldTypeOID, alias: OBJECTID )
##     Shape ( type: esriFieldTypeGeometry, alias: Shape )
##     OBJECTID_1 ( type: esriFieldTypeInteger, alias: OBJECTID_1 )
##     ID ( type: esriFieldTypeInteger, alias: ID )
##     FileSize ( type: esriFieldTypeDouble, alias: FileSize )
##     pixbytes ( type: esriFieldTypeInteger, alias: pixbytes )
##     DataTypeID ( type: esriFieldTypeInteger, alias: DataTypeID )
##     provider_results ( type: esriFieldTypeString, alias: provider_results, length: 1000 )
##     provider_details ( type: esriFieldTypeString, alias: provider_details, length: 1000 )
##     licStatus ( type: esriFieldTypeInteger, alias: licStatus )
##     Name ( type: esriFieldTypeString, alias: Name, length: 200 )
##     provider_results_name ( type: esriFieldTypeString, alias: provider_results_name, length: 2147483647 )
##     provider_details_name ( type: esriFieldTypeString, alias: provider_details_name, length: 2147483647 )
##     DataType ( type: esriFieldTypeString, alias: DataType, length: 7 )
##     DataBin ( type: esriFieldTypeString, alias: DataBin, length: 9 )
##     Year ( type: esriFieldTypeInteger, alias: Year )
##     ProjectID ( type: esriFieldTypeInteger, alias: ProjectID )
##     Project ( type: esriFieldTypeString, alias: Project, length: 150 )
##     Project_Description ( type: esriFieldTypeString, alias: Project_Description, length: 8000 )
##     dclink ( type: esriFieldTypeString, alias: dclink, length: 200 )
##     Metalink ( type: esriFieldTypeString, alias: Metalink, length: 4000 )
##     licLink ( type: esriFieldTypeString, alias: licLink, length: 256 )
##     imgname ( type: esriFieldTypeString, alias: imgname, length: 250 )
##     InfoLink ( type: esriFieldTypeString, alias: InfoLink, length: 200 )
##     SpecialNote ( type: esriFieldTypeString, alias: SpecialNote, length: 8000 )
##     ProvisioningDetails ( type: esriFieldTypeString, alias: ProvisioningDetails, length: 8000 )
##     ExternalProviderLink ( type: esriFieldTypeString, alias: ExternalProviderLink, length: 2147483647 )
##     ExternalProviderLinkLabel ( type: esriFieldTypeString, alias: ExternalProviderLinkLabel, length: 14 )
##     ExternalParameters ( type: esriFieldTypeString, alias: ExternalParameters, length: 100 )
##     ExternalParametersAlias ( type: esriFieldTypeString, alias: ExternalParametersAlias, length: 100 )
##     Vertical_Accuracy ( type: esriFieldTypeString, alias: Vertical_Accuracy, length: 313 )
##     Horizontal_Accuracy ( type: esriFieldTypeString, alias: Horizontal_Accuracy, length: 313 )
##     Nominal_Ground_Spacing ( type: esriFieldTypeDouble, alias: Nominal_Ground_Spacing )
##     Data_Classes_Available ( type: esriFieldTypeString, alias: Data_Classes_Available, length: 2147483647 )
##     TideControlled ( type: esriFieldTypeString, alias: TideControlled, length: 3 )
##     NativeVdatum ( type: esriFieldTypeString, alias: NativeVdatum, length: 20 )
##     Classified ( type: esriFieldTypeString, alias: Classified, length: 2147483647 )
##     ReturnsOption ( type: esriFieldTypeInteger, alias: ReturnsOption )
##     AncillaryData ( type: esriFieldTypeString, alias: AncillaryData, length: 100 )
##     AncillaryOpt ( type: esriFieldTypeInteger, alias: AncillaryOpt )
##     AllowLAS ( type: esriFieldTypeInteger, alias: AllowLAS )
##     CellSizeFt ( type: esriFieldTypeDouble, alias: CellSizeFt )
##     CellSizeM ( type: esriFieldTypeDouble, alias: CellSizeM )
##     MinCellSizeFt ( type: esriFieldTypeDouble, alias: MinCellSizeFt )
##     MinCellSizeMeters ( type: esriFieldTypeDouble, alias: MinCellSizeMeters )
##     MinContourIntervalFt ( type: esriFieldTypeString, alias: MinContourIntervalFt, length: 100 )
##     ImageService_Server ( type: esriFieldTypeString, alias: ImageService_Server, length: 4000 )
##     ImageService_Service ( type: esriFieldTypeString, alias: ImageService_Service, length: 200 )
##     ImageService_Key ( type: esriFieldTypeString, alias: ImageService_Key, length: 50 )
##     ImageService_Value ( type: esriFieldTypeString, alias: ImageService_Value, length: 50 )
##     ImageService_FullURL ( type: esriFieldTypeString, alias: ImageService_FullURL, length: 4000 )
##     Shape_Length ( type: esriFieldTypeDouble, alias: Shape_Length )
##     Shape_Area ( type: esriFieldTypeDouble, alias: Shape_Area )
##
### Code:

import os
import json

from osgeo import ogr
from osgeo import gdal

from cudem import utils
from cudem import regions
from cudem import datasets
from cudem import demfun
from cudem import vdatums

import cudem.fetches.utils as f_utils
import cudem.fetches.FRED as FRED

class DAV(f_utils.FetchModule):
    """Fetch NOAA lidar data from DAV

Uses Digital Coasts Data Access Viewer Mapserver to discover
dataset footprints.

Use where=SQL_QUERY to query the MapServer to filter datasets

If inc is set, upon processing the data, will blockmedian it at `inc`
increment to save space.
"""
    
    def __init__(self, where='1=1', index=False, datatype=None, **kwargs):
        super().__init__(**kwargs)
        self._dav_api_url = 'https://maps.coast.noaa.gov/arcgis/rest/services/DAV/ElevationFootprints/MapServer/0/query?'
        self._outdir = os.path.join(os.getcwd(), 'digital_coast')
        self.name = 'digital_coast'
        self.where = where
        self.index = index
        self.datatype = datatype
        
    def run(self):
        '''Run the DAV fetching module'''
        
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
        _req = f_utils.Fetch(self._dav_api_url, verbose=self.verbose).fetch_req(params=_data)
        print(_req.url)
        if _req is not None:
            features = _req.json()
            if not 'features' in features.keys():
                utils.echo_error_msg('DAV failed to execute query...try again later')
            else:
                for feature in features['features']:

                    if self.datatype is not None:
                        if self.datatype.lower() != feature['attributes']['DataType'].lower():
                            continue

                    links = json.loads(feature['attributes']['ExternalProviderLink'])

                    ept_infos = None
                    ## get ept link to gather datum infos...for lidar only apparently...
                    for link in links['links']:
                        if link['serviceID'] == 167:
                            if link['label'] == 'EPT NOAA':
                                ept_req = f_utils.Fetch(link['link'], verbose=True).fetch_req()
                                ept_infos = ept_req.json()

                    ## get metadata for datum infos...for raster data

                    if self.index:
                        feature['attributes']['ExternalProviderLink'] = links
                        print(json.dumps(feature['attributes'], indent=4))
                    else:
                        for link in links['links']:
                            if link['serviceID'] == 46:
                                urllist = 'urllist' + str(feature['attributes']['ID']) + '.txt'
                                surv_name = '_'.join(link['link'].split('/')[-1].split('_')[:-1])
                                #index_zipfile = 'tileindex.zip'
                                index_zipfile = 'tileindex_{}.zip'.format(surv_name)
                                index_zipurl = link['link'] + '/' + index_zipfile
                                #urllist_url = '/'.join(link['link'].split('/')[:-1]) + '/' + urllist
                                urllist_url = link['link'] + '/' + urllist
                                while True:
                                    if f_utils.Fetch(urllist_url, verbose=True).fetch_file(urllist) != 0:
                                        if urllist_url == '/'.join(link['link'].split('/')[:-1]) + '/' + urllist:
                                            break
                                        urllist_url = '/'.join(link['link'].split('/')[:-1]) + '/' + urllist
                                    else:
                                        break

                                with open(urllist, 'r') as ul:
                                    for line in ul:
                                        if 'tileindex' in line:
                                            index_zipurl = line.strip()
                                            break

                                utils.remove_glob(urllist)                            
                                if f_utils.Fetch(
                                        index_zipurl, callback=self.callback, verbose=self.verbose
                                ).fetch_file(index_zipfile) == 0:
                                    index_shps = utils.p_unzip(index_zipfile, ['shp', 'shx', 'dbf', 'prj'])
                                    index_shp = None
                                    for v in index_shps:
                                        if v.split('.')[-1] == 'shp':
                                            index_shp = v

                                    index_ds = ogr.Open(index_shp)
                                    index_layer = index_ds.GetLayer(0)
                                    for index_feature in index_layer:
                                        index_geom = index_feature.GetGeometryRef()
                                        if index_geom.Intersects(self.region.export_as_geom()):
                                            tile_url = index_feature.GetField('URL').strip()
                                            ## add vertical datum to output;
                                            ## field is NativeVdatum
                                            ## must get from metadata
                                            if ept_infos is None:
                                                this_epsg = vdatums.get_vdatum_by_name(feature['attributes']['NativeVdatum'])
                                            else:
                                                ## horizontal datum is wrong in ept, most seem to be nad83
                                                #this_epsg = 'epsg:{}+{}'.format(ept_infos['srs']['horizontal'], ept_infos['srs']['vertical'])
                                                this_epsg = 'epsg:4269+{}'.format(ept_infos['srs']['vertical'])

                                            self.results.append(
                                                [tile_url,
                                                 os.path.join(
                                                     self._outdir, '{}/{}'.format(
                                                         feature['attributes']['ID'], tile_url.split('/')[-1]
                                                     )
                                                 ),
                                                 this_epsg,
                                                 feature['attributes']['DataType']])

                                    index_ds = index_layer = None
                                    utils.remove_glob(index_zipfile, *index_shps)

        return(self)

    def yield_xyz(self, entry):
        src_dc = os.path.basename(entry[1])
        src_ext = src_dc.split('.')[-1].lower()
        if src_ext == 'laz' or src_ext == 'las':
            dt = 'lidar'
        elif src_ext == 'tif' or src_ext == 'img':
            dt = 'raster'
        else:
            dt = None

        if dt == 'lidar':
            if f_utils.Fetch(entry[0], callback=self.callback, verbose=self.verbose).fetch_file(src_dc) == 0:
                ds_epsg = str(entry[2])
                ## no ept epsg's...assume nad83, NAVD88
                if len(ds_epsg.split(':')) == 1:
                    ds_epsg = 'epsg:4269+{}'.format(ds_epsg)

                _ds = datasets.LASFile(
                    fn=src_dc,
                    data_format=300,
                    src_srs=ds_epsg,
                    dst_srs=self.dst_srs,
                    weight=self.weight,
                    src_region=self.region,
                    x_inc=self.x_inc,
                    y_inc=self.y_inc,
                    verbose=self.verbose,
                    remote=True
                )
                # if self.inc is not None:
                #     b_region = regions.regions_reduce(self.region, regions.Region().from_list(_ds.infos['minmax']))
                #     xyz_func = lambda p: _ds.dump_xyz(dst_port=p, encode=True)
                #     for xyz in utils.yield_cmd(
                #             'gmt blockmedian -I{:.10f} {} -r -V'.format(self.inc, b_region.format('gmt')),
                #             verbose=self.verbose,
                #             data_fun=xyz_func
                #     ):
                #         yield(xyzfun.XYZPoint().from_list([float(x) for x in xyz.split()]))
                        
                # else:
                #for xyz in _ds.yield_xyz():
                for xyz in _ds.yield_xyz():
                    yield(xyz)
                        
                utils.remove_glob('{}*'.format(src_dc))
        elif dt == 'raster':
            if f_utils.Fetch(entry[0], callback=self.callback, verbose=self.verbose).fetch_file(src_dc) == 0:
                ds_epsg = str(entry[2])
                if len(ds_epsg.split(':')) == 1:
                    ds_epsg = '{}+{}'.format(demfun.get_srs(src_dc), ds_epsg)
                _ds = datasets.RasterFile(
                    fn=src_dc,
                    data_format=200,
                    src_srs=ds_epsg,
                    dst_srs=self.dst_srs,
                    weight=self.weight,
                    src_region=self.region,
                    x_inc=self.x_inc,
                    y_inc=self.y_inc,
                    verbose=self.verbose
                )
                for xyz in _ds.yield_xyz():
                    yield(xyz)
                    
                utils.remove_glob('{}.*'.format(src_dc))

## ==============================================
## the DigitalCoast class is the old digital coast fetches module.
## This module scrapes the data from digital coast and generates
## a reference vector to discover dataset footprints. This is prone
## to possible error and can miss newer datasets if the reference
## vector is not up-to-date. Use DAV class instead, which uses
## the Digital Coast MapServer to discover dataset footprints.
## ==============================================
class DigitalCoast(f_utils.FetchModule):

    def __init__(self, where='', datatype=None, inc=None, **kwargs):
        super().__init__(**kwargs)
        self._dc_url = 'https://coast.noaa.gov'
        self._dc_htdata_url = 'https://coast.noaa.gov/htdata/'
        #self._dc_dirs = ['lidar1_z', 'lidar2_z', 'lidar3_z', 'lidar4_z', 'raster1', 'raster2', 'raster5']
        self._dc_dirs = ['lidar1_z', 'raster1']
        self._outdir = os.path.join(os.getcwd(), 'digital_coast')
        self.where = [where] if len(where) > 0 else []
        self.datatype = datatype
        self.inc = utils.str2inc(inc)
        self.name = 'dc'
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
                        this_xml = f_utils.iso_xml(dc['Metadata'])
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
                        self.results.append([tile_url, os.path.join(self._outdir, '{}/{}'.format(surv['ID'], tile_url.split('/')[-1])), surv['DataType']])
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
                # xyz_dat = utils.yield_cmd('las2txt -stdout -parse xyz -keep_xy {} -keep_class {} -i {}\
                # '.format(self.region.format('te'), '2 29', src_dc), verbose = False)
                # _ds = datasets.XYZFile(fn=xyz_dat, data_format=168, dst_srs=self.dst_Srs,
                #                        name=xyz_dat, src_region=self.region, verbose=self.verbose, remote=True)
                #xyz_dat = utils.yield_cmd('las2txt -stdout -parse xyz -keep_xy {} -keep_class {} -i {}\
                #'.format(self.region.format('te'), '2 29', src_dc), verbose = False)
                _ds = datasets.LASFile(
                    fn=src_dc,
                    data_format=400,
                    dst_srs=self.dst_srs,
                    name=src_dc,
                    src_region=self.region,
                    verbose=self.verbose,
                    remote=True
                )
                if self.inc is not None:
                    b_region = regions.regions_reduce(self.region, regions.Region().from_list(_ds.infos['minmax']))
                    xyz_func = lambda p: _ds.dump_xyz(dst_port=p, encode=True)
                    for xyz in utils.yield_cmd(
                            'gmt blockmedian -I{:.10f} {} -r -V'.format(self.inc, b_region.format('gmt')),
                            verbose=self.verbose,
                            data_fun=xyz_func
                    ):
                        yield(xyzfun.XYZPoint().from_list([float(x) for x in xyz.split()]))
                        
                else:
                    for xyz in _ds.yield_xyz():
                        yield(xyz)

                #for xyz in _ds.block_xyz(inc=self.inc, want_gmt=True) if self.inc is not None else _ds.yield_xyz():
                #    yield(xyz)
                #y = _ds.block_xyz if self.inc is not None else _ds.yield_xyz

                #for xyz in y():
                #    yield(xyz)
                        
        elif dt == 'raster':
            #try:
            #    src_ds = gdal.Open(entry[0])
            #    src_dc = entry[0]
            #except Exception as e:
            if f_utils.Fetch(entry[0], callback=self.callback, verbose=self.verbose).fetch_file(src_dc) == 0:
                #try:
                #    src_ds = gdal.Open(src_dc)
                #except Exception as e:
                #    utils.echo_error_msg('could not read dc raster file: {}, {}'.format(entry[0], e))
                #    src_ds = None
                #except Exception as e:
                #    utils.echo_error_msg('could not read dc raster file: {}, {}'.format(entry[0], e))
                #    src_ds = None
                
                #if src_ds is not None:
                #src_ds = None
                _ds = datasets.RasterFile(
                    fn=src_dc,
                    data_format=200,
                    dst_srs=self.dst_srs,
                    src_srs=None,
                    name=src_dc,
                    src_region=self.region,
                    verbose=self.verbose
                )
                #_ds.src_ds = src_ds
                #_ds.ds_open_p = True
                for xyz in _ds.block_xyz(inc=self.inc, want_gmt=True) if self.inc is not None else _ds.yield_xyz():
                    yield(xyz)
                #src_ds = None
                utils.remove_glob(src_dc)    
### End
