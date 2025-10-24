### dav.py
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

import os
import json
from osgeo import ogr

from cudem import utils
from cudem import gdalfun
from cudem import vdatums
from cudem.fetches import fetches

## Digital Coast - Data Access Viewer
class DAV_old(fetches.FetchModule):
    """Fetch NOAA lidar data from DAV

    Uses Digital Coasts Data Access Viewer Mapserver to discover
    dataset footprints.

    This map service presents spatial information about Elevation Data Access 
    Viewer services across the United States and Territories in the Web 
    Mercator projection. The service was developed by the National Oceanic 
    and Atmospheric Administration (NOAA), but may contain data and information 
    from a variety of data sources, including non-NOAA data. 

    NOAA provides the information “as-is” and shall incur no responsibility or 
    liability as to the completeness or accuracy of this information. NOAA 
    assumes no responsibility arising from the use of this information. 
    The NOAA Office for Coastal Management will make every effort to provide 
    continual access to this service but it may need to be taken down during
    routine IT maintenance or in case of an emergency. If you plan to ingest 
    this service into your own application and would like to be informed about 
    planned and unplanned service outages or changes to existing services, 
    please register for our Data Services Newsletter 
    (http://coast.noaa.gov/digitalcoast/publications/subscribe). 

    For additional information, please contact the NOAA Office for Coastal 
    Management (coastal.info@noaa.gov).

    Fields:

    OBJECTID ( type: esriFieldTypeOID, alias: OBJECTID )
    Shape ( type: esriFieldTypeGeometry, alias: Shape )
    OBJECTID_1 ( type: esriFieldTypeInteger, alias: OBJECTID_1 )
    ID ( type: esriFieldTypeInteger, alias: ID )
    FileSize ( type: esriFieldTypeDouble, alias: FileSize )
    pixbytes ( type: esriFieldTypeInteger, alias: pixbytes )
    DataTypeID ( type: esriFieldTypeInteger, alias: DataTypeID )
    provider_results ( type: esriFieldTypeString, alias: provider_results, length: 1000 )
    provider_details ( type: esriFieldTypeString, alias: provider_details, length: 1000 )
    licStatus ( type: esriFieldTypeInteger, alias: licStatus )
    Name ( type: esriFieldTypeString, alias: Name, length: 200 )
    provider_results_name ( type: esriFieldTypeString, alias: provider_results_name, length: 2147483647 )
    provider_details_name ( type: esriFieldTypeString, alias: provider_details_name, length: 2147483647 )
    DataType ( type: esriFieldTypeString, alias: DataType, length: 7 )
    DataBin ( type: esriFieldTypeString, alias: DataBin, length: 9 )
    Year ( type: esriFieldTypeInteger, alias: Year )
    ProjectID ( type: esriFieldTypeInteger, alias: ProjectID )
    Project ( type: esriFieldTypeString, alias: Project, length: 150 )
    Project_Description ( type: esriFieldTypeString, alias: Project_Description, length: 8000 )
    dclink ( type: esriFieldTypeString, alias: dclink, length: 200 )
    Metalink ( type: esriFieldTypeString, alias: Metalink, length: 4000 )
    licLink ( type: esriFieldTypeString, alias: licLink, length: 256 )
    imgname ( type: esriFieldTypeString, alias: imgname, length: 250 )
    InfoLink ( type: esriFieldTypeString, alias: InfoLink, length: 200 )
    SpecialNote ( type: esriFieldTypeString, alias: SpecialNote, length: 8000 )
    ProvisioningDetails ( type: esriFieldTypeString, alias: ProvisioningDetails, length: 8000 )
    ExternalProviderLink ( type: esriFieldTypeString, alias: ExternalProviderLink, length: 2147483647 )
    ExternalProviderLinkLabel ( type: esriFieldTypeString, alias: ExternalProviderLinkLabel, length: 14 )
    ExternalParameters ( type: esriFieldTypeString, alias: ExternalParameters, length: 100 )
    ExternalParametersAlias ( type: esriFieldTypeString, alias: ExternalParametersAlias, length: 100 )
    Vertical_Accuracy ( type: esriFieldTypeString, alias: Vertical_Accuracy, length: 313 )
    Horizontal_Accuracy ( type: esriFieldTypeString, alias: Horizontal_Accuracy, length: 313 )
    Nominal_Ground_Spacing ( type: esriFieldTypeDouble, alias: Nominal_Ground_Spacing )
    Data_Classes_Available ( type: esriFieldTypeString, alias: Data_Classes_Available, length: 2147483647 )
    TideControlled ( type: esriFieldTypeString, alias: TideControlled, length: 3 )
    NativeVdatum ( type: esriFieldTypeString, alias: NativeVdatum, length: 20 )
    Classified ( type: esriFieldTypeString, alias: Classified, length: 2147483647 )
    ReturnsOption ( type: esriFieldTypeInteger, alias: ReturnsOption )
    AncillaryData ( type: esriFieldTypeString, alias: AncillaryData, length: 100 )
    AncillaryOpt ( type: esriFieldTypeInteger, alias: AncillaryOpt )
    AllowLAS ( type: esriFieldTypeInteger, alias: AllowLAS )
    CellSizeFt ( type: esriFieldTypeDouble, alias: CellSizeFt )
    CellSizeM ( type: esriFieldTypeDouble, alias: CellSizeM )
    MinCellSizeFt ( type: esriFieldTypeDouble, alias: MinCellSizeFt )
    MinCellSizeMeters ( type: esriFieldTypeDouble, alias: MinCellSizeMeters )
    MinContourIntervalFt ( type: esriFieldTypeString, alias: MinContourIntervalFt, length: 100 )
    ImageService_Server ( type: esriFieldTypeString, alias: ImageService_Server, length: 4000 )
    ImageService_Service ( type: esriFieldTypeString, alias: ImageService_Service, length: 200 )
    ImageService_Key ( type: esriFieldTypeString, alias: ImageService_Key, length: 50 )
    ImageService_Value ( type: esriFieldTypeString, alias: ImageService_Value, length: 50 )
    ImageService_FullURL ( type: esriFieldTypeString, alias: ImageService_FullURL, length: 4000 )
    Shape_Length ( type: esriFieldTypeDouble, alias: Shape_Length )
    Shape_Area ( type: esriFieldTypeDouble, alias: Shape_Area )
    
    Layers:
    Elevation > 1:20M (0)
    Elevation 1:12.5M - 1:20M (1)
    Elevation 1:1.25M - 1:12.5M (2)
    Elevation < 1:1.25M scale (3)

    https://coast.noaa.gov
    
    Use where=SQL_QUERY to query the MapServer to filter datasets

    * For CUDEM tiles, use where="ID=8483" (1/9) or where="ID=8580" (1/3) or where="Name LIKE '%CUDEM%'" for all available.
    * For OCM SLR DEMs, use where="ID=6230" or where="Name LIKE '%Sea Level Rise%'"
    * For USGS CoNED DEMs, use where="ID=9181" or where="Name LIKE '%CoNED%'"
    * To only return lidar data, use datatype=lidar, for only raster, use datatype=dem
    * datatype is either 'lidar', 'dem' or 'sm'

    < digital_coast:where=None:datatype=None:footprints_only=False >
    """
    
    def __init__(
            self,
            where='1=1',
            index=False,
            datatype=None,
            layer=0,
            name='digital_coast',
            want_footprints=False,
            footprints_only=False,
            **kwargs
    ):
        super().__init__(name=name, **kwargs)
        self.where = where
        self.index = index
        self.datatype = datatype
        self.layer = utils.int_or(layer)
    
        ## The various DAV URLs
        # self._dav_api_url = ('https://maps.coast.noaa.gov/arcgis/rest/services/'
        #                      f'DAV/ElevationFootprints/MapServer/{layer}/query?')
        self._dav_api_url = ('https://maps.coast.noaa.gov/arcgis/rest/services/'
                             'DAV/ElevationFootprints/MapServer/')

        ## data formats vary
        self.data_format = None
        self.want_footprints = want_footprints
        self.footprints_only = footprints_only
        if self.footprints_only:
            self.want_footprints = True

            
    def run(self):
        '''Run the DAV fetching module'''
        
        if self.region is None:
            return([])

        _data = {
            'where': self.where,
            'outFields': '*',
            'geometry': self.region.format('bbox'),
            'inSR':4326,
            #'outSR':4326,
            'f':'pjson',
            'returnGeometry':'False',
        }
        _req = fetches.Fetch(
            self._dav_api_url + str(self.layer) + '/query?', verbose=self.verbose
        ).fetch_req(params=_data)
        if _req is not None:
            features = _req.json()
            if not 'features' in features.keys():
                utils.echo_error_msg(
                    'DAV failed to execute query...try again later'
                )
                utils.echo_msg(_req.text)
                utils.echo_msg(_req.url)
            else:
                for feature in features['features']:
                    if self.datatype is not None:
                        if self.datatype.lower() != feature['attributes']['DataType'].lower():
                            if self.datatype.lower() != 'sm':
                                continue

                    self.vdatum = vdatums.get_vdatum_by_name(
                        feature['attributes']['NativeVdatum']
                    ) 
                    links = json.loads(
                        feature['attributes']['ExternalProviderLink']
                    )
                    ept_infos = None
                    ## get ept link to gather datum infos...
                    ## for lidar only apparently...
                    for link in links['links']:
                        if link['serviceID'] == 167:
                            if link['label'] == 'EPT NOAA':
                                ept_req = fetches.Fetch(link['link'], verbose=True).fetch_req()
                                ept_infos = ept_req.json()

                    ## get metadata for datum infos...for raster data
                    if self.index:
                        feature['attributes']['ExternalProviderLink'] = links
                        utils.echo_msg(json.dumps(feature['attributes'], indent=4))
                    else:
                        for link in links['links']:
                            if link['serviceID'] == 46 \
                               and (self.datatype == 'lidar' \
                                    or self.datatype == 'dem' \
                                    or self.datatype is None):
                                surv_name = '_'.join(link['link'].split('/')[-2].split('_'))
                                index_zipfile = os.path.join(
                                    self._outdir, 'tileindex_{}.zip'.format(surv_name)
                                )
                                page = fetches.Fetch(link['link'], verbose=False).fetch_html()
                                rows = page.xpath('//a[contains(@href, ".txt")]/@href')
                                for l in rows:
                                    if 'urllist' in l:
                                        urllist = l
                                        break

                                if 'http' in urllist:
                                    urllist_url = urllist
                                else:
                                    urllist_url = os.path.dirname(link['link']) + '/' + urllist

                                urllist = os.path.join(self._outdir, os.path.basename(urllist))
                                status = fetches.Fetch(urllist_url, verbose=True).fetch_file(urllist)
                                if not os.path.exists(urllist):
                                    continue
                                
                                with open(urllist, 'r') as ul:
                                    for line in ul:
                                        if 'tileindex' in line and 'zip' in line:
                                            index_zipurl = line.strip()
                                            break
                                        
                                utils.remove_glob(urllist)
                                if self.want_footprints:
                                    self.add_entry_to_results(
                                        index_zipurl,
                                        os.path.join(
                                            '{}/{}'.format(
                                                feature['attributes']['ID'],
                                                index_zipurl.split('/')[-1]
                                            )
                                        ),
                                        'footprint'
                                    )
                                    if self.footprints_only:
                                        continue
                                    
                                try:
                                    status = fetches.Fetch(
                                        index_zipurl,
                                        callback=self.callback,
                                        verbose=self.verbose
                                    ).fetch_file(index_zipfile)
                                except:
                                    status = -1
                                    
                                #utils.echo_msg(index_zipfile)
                                if status == 0:
                                    index_shps = utils.p_unzip(
                                        index_zipfile, ['shp', 'shx', 'dbf', 'prj'],
                                        outdir=self._outdir,
                                        verbose=True
                                    )
                                    index_shp = None
                                    for v in index_shps:
                                        if v.split('.')[-1] == 'shp':
                                            index_shp = v

                                    warp_region = self.region.copy()
                                    index_prj = None
                                    for v in index_shps:
                                        if v.split('.')[-1] == 'prj':
                                            index_prj = v
                                            
                                    if index_prj is not None:
                                        prj_ds = open(index_prj)
                                        prj_wkt = prj_ds.read()
                                        warp_region.warp(dst_crs = prj_wkt)
                                        prj_ds.close()
                                        
                                    index_ds = ogr.Open(index_shp)
                                    index_layer = index_ds.GetLayer(0)
                                    for index_feature in index_layer:
                                        index_geom = index_feature.GetGeometryRef()
                                        known_name_fields = ['Name', 'location', 'filename']
                                        known_url_fields = ['url', 'URL']
                                        
                                        if index_geom.Intersects(warp_region.export_as_geom()):
                                            tile_name = None
                                            tile_url = None
                                            field_names = [field.name for field in index_layer.schema]
                                            
                                            for f in known_name_fields:
                                                if f in field_names:
                                                    tile_name = index_feature.GetField(f).strip()
                                                    
                                                if tile_name is not None:
                                                    break

                                            for f in known_url_fields:
                                                if f in field_names:
                                                    tile_url = index_feature.GetField(f).strip()
                                                    
                                                if tile_url is not None:
                                                    break

                                            if tile_name is None or tile_url is None:
                                                utils.echo_warning_msg('could not parse index fields')
                                                continue
                                            
                                            tile_url = '/'.join(tile_url.split('/')[:-1]) + '/' + tile_name.split('/')[-1]
                                            ## add vertical datum to output;
                                            ## field is NativeVdatum
                                            ## must get from metadata
                                            if ept_infos is None:
                                                this_epsg = vdatums.get_vdatum_by_name(
                                                    feature['attributes']['NativeVdatum']
                                                )
                                            else:
                                                #print(ept_infos['srs'])
                                                ## horizontal datum is wrong in ept, most seem to be nad83
                                                #this_epsg = 'epsg:{}+{}'.format(ept_infos['srs']['horizontal'], ept_infos['srs']['vertical'])
                                                if 'vertical' in ept_infos['srs'].keys():
                                                    vertical_epsg = ept_infos['srs']['vertical']
                                                    horizontal_epsg = ept_infos['srs']['horizontal']
                                                    this_epsg = 'epsg:{}+{}'.format(
                                                        horizontal_epsg, vertical_epsg
                                                    )
                                                    #this_epsg = 'epsg:4269+{}'.format(ept_infos['srs']['vertical'])
                                                else:
                                                    # try to extract the vertical datum from the wkt
                                                    #horizontal_epsg = ept_infos['srs']['horizontal']
                                                    this_wkt = ept_infos['srs']['wkt']
                                                    dst_horz, dst_vert = gdalfun.epsg_from_input(
                                                        this_wkt
                                                    )
                                                    this_epsg = '{}+{}'.format(
                                                        dst_horz, dst_vert
                                                    )

                                            self.add_entry_to_results(
                                                tile_url,
                                                os.path.join(
                                                    '{}/{}'.format(
                                                        feature['attributes']['ID'],
                                                        tile_url.split('/')[-1]
                                                    )
                                                ),
                                                feature['attributes']['DataType'],
                                                this_epsg=this_epsg,
                                            )

                                    index_ds = index_layer = None
                                    #if not self.want_footprints:
                                    utils.remove_glob(index_zipfile, *index_shps)
                                    #utils.remove_glob(*index_shps)

                            # spatial_metadata
                            elif link['serviceID'] == 166 and self.datatype == 'sm': 
                                self.add_entry_to_results(
                                    link['link'],
                                    os.path.join(
                                        '{}/{}'.format(
                                            feature['attributes']['ID'],
                                            link['link'].split('/')[-1]
                                        )
                                    ),
                                    link['label'],
                                    this_epsg=None
                                )

        #self.results = [x for x in np.unique(self.results, axis=0)]
        return(self)


class DAV(fetches.FetchModule):
    """Fetch NOAA lidar data from DAV

    Uses Digital Coasts Data Access Viewer Mapserver to discover
    dataset footprints.

    This map service presents spatial information about Elevation Data Access 
    Viewer services across the United States and Territories in the Web 
    Mercator projection. The service was developed by the National Oceanic 
    and Atmospheric Administration (NOAA), but may contain data and information 
    from a variety of data sources, including non-NOAA data. 

    NOAA provides the information “as-is” and shall incur no responsibility or 
    liability as to the completeness or accuracy of this information. NOAA 
    assumes no responsibility arising from the use of this information. 
    The NOAA Office for Coastal Management will make every effort to provide 
    continual access to this service but it may need to be taken down during
    routine IT maintenance or in case of an emergency. If you plan to ingest 
    this service into your own application and would like to be informed about 
    planned and unplanned service outages or changes to existing services, 
    please register for our Data Services Newsletter 
    (http://coast.noaa.gov/digitalcoast/publications/subscribe). 

    For additional information, please contact the NOAA Office for Coastal 
    Management (coastal.info@noaa.gov).

    https://coast.noaa.gov
    
    Use where=SQL_QUERY to query the MapServer to filter datasets

    * For CUDEM tiles, use where="ID=8483" (1/9) or where="ID=8580" (1/3) or where="Name LIKE '%CUDEM%'" for all available.
    * For OCM SLR DEMs, use where="ID=6230" or where="Name LIKE '%Sea Level Rise%'"
    * For USGS CoNED DEMs, use where="ID=9181" or where="Name LIKE '%CoNED%'"
    * To only return lidar data, use datatype=lidar, for only raster, use datatype=dem
    * datatype is either 'lidar', 'dem' or 'sm'

    < digital_coast:where=None:datatype=None:footprints_only=False >
    """
    
    def __init__(
            self,
            where='1=1',
            index=False,
            datatype=None,
            layer=0,
            name='digital_coast',
            want_footprints=False,
            keep_footprints=False,
            #want_urllist=False,
            footprints_only=False,
            **kwargs
    ):
        super().__init__(name=name, **kwargs)
        self.where = where
        self.index = index
        self.datatype = datatype
        self.layer = utils.int_or(layer)
    
        ## The various DAV URLs
        self._dav_api_url = ('https://maps.coast.noaa.gov/arcgis/rest/services/'
                             'DAV/DAV_footprints/MapServer/')

        ## data formats vary
        self.data_format = None
        #self.want_urllist = want_urllist
        self.want_footprints = want_footprints
        self.keep_footprints = keep_footprints
        self.footprints_only = footprints_only
        if self.footprints_only:
            self.want_footprints = True

            
    def run(self):
        '''Run the DAV fetching module'''
        
        if self.region is None:
            return([])

        _data = {
            'where': self.where,
            'outFields': '*',
            'geometry': self.region.format('bbox'),
            'inSR':4326,
            #'outSR':4326,
            'f':'pjson',
            'returnGeometry':'False',
        }
        _req = fetches.Fetch(
            self._dav_api_url + str(self.layer) + '/query?', verbose=self.verbose
        ).fetch_req(params=_data)
        if _req is not None:
            features = _req.json()
            if not 'features' in features.keys():
                utils.echo_error_msg(
                    'DAV failed to execute query...try again later'
                )
                utils.echo_msg(_req.text)
                utils.echo_msg(_req.url)
            else:
                for feature in features['features']:
                    if self.datatype is not None:
                        if self.datatype.lower() != feature['attributes']['data_type'].lower():
                            if self.datatype.lower() != 'sm':
                                continue

                    fid = feature['attributes']['id']
                    metadata_link = feature['attributes']['metadata']
                    page = fetches.Fetch(metadata_link, verbose=False).fetch_html()
                    if page is None:
                        utils.echo_msg(metadata_link)
                        continue
                    
                    index_urls = page.xpath('//a[contains(@href, "index.html")]/@href')
                    if len(index_urls) == 0:
                        #utils.echo_msg(metadata_link)
                        #page = fetches.Fetch(metadata_link, verbose=False).fetch_xml()
                        #index_urls = page.findall('.//distributorTransferOptions')
                        index_urls = page.xpath(f'//a[contains(@href, "_{fid}/")]/@href')

                    if len(index_urls) == 0:
                        continue
                    
                    # for link in links['links']:
                    #     if link['serviceID'] == 46 \
                    #        and (self.datatype == 'lidar' \
                    #             or self.datatype == 'dem' \
                    #             or self.datatype is None):
                    surv_name = '_'.join(index_urls[0].split('/')[-2].split('_'))
                    index_zipfile = os.path.join(
                        self._outdir, 'tileindex_{}.zip'.format(surv_name)
                    )
                    page = fetches.Fetch(index_urls[0], verbose=False).fetch_html()
                    if page is None:
                        utils.echo_warning_msg(f'could not fetch {index_urls}')
                        continue
                    
                    rows = page.xpath('//a[contains(@href, ".txt")]/@href')
                    #rows = page.xpath('//a[contains(@href, "u")]/@href')
                    #utils.echo_msg_bold(page)
                    #utils.echo_msg_bold(rows)
                    urllist = None
                    for l in rows:
                        #utils.echo_msg_bold(l)
                        if 'urllist' in l:
                            urllist = l
                            break

                    if urllist is None:
                        utils.echo_warning_msg(f'could not find urllist from {feature}')
                        utils.echo_msg(index_urls)
                        continue
                     
                    #utils.echo_msg_bold(urllist)
                    if 'http' in urllist:
                        urllist_url = urllist
                    else:
                        urllist_url = os.path.dirname(index_urls[0]) + '/' + urllist

                    urllist = os.path.join(self._outdir, os.path.basename(urllist))
                    status = fetches.Fetch(urllist_url, verbose=True).fetch_file(urllist)
                    if not os.path.exists(urllist):
                        continue

                    with open(urllist, 'r') as ul:
                        for line in ul:
                            if 'tileindex' in line and 'zip' in line:
                                index_zipurl = line.strip()
                                break

                    #if not self.want_urllist:
                    utils.remove_glob(urllist)                        
                    if self.want_footprints:
                        self.add_entry_to_results(
                            index_zipurl,
                            os.path.join(
                                '{}/{}'.format(
                                    feature['attributes']['id'],
                                    index_zipurl.split('/')[-1]
                                )
                            ),
                            'footprint'
                        )
                        if self.footprints_only:
                            continue

                    try:
                        status = fetches.Fetch(
                            index_zipurl,
                            callback=self.callback,
                            verbose=self.verbose
                        ).fetch_file(index_zipfile)
                    except:
                        status = -1

                    #utils.echo_msg(index_zipurl)
                    #utils.echo_msg(status)
                    if status == 0:
                        index_shps = utils.p_unzip(
                            index_zipfile, ['shp', 'shx', 'dbf', 'prj'],
                            outdir=self._outdir,
                            verbose=True
                        )
                        index_shp = None
                        for v in index_shps:
                            if v.split('.')[-1] == 'shp':
                                index_shp = v

                        warp_region = self.region.copy()
                        index_prj = None
                        for v in index_shps:
                            if v.split('.')[-1] == 'prj':
                                index_prj = v

                        if index_prj is not None:
                            prj_ds = open(index_prj)
                            prj_wkt = prj_ds.read()
                            warp_region.warp(dst_crs = prj_wkt)
                            prj_ds.close()

                        index_ds = ogr.Open(index_shp)
                        index_layer = index_ds.GetLayer(0)
                        for index_feature in index_layer:
                            index_geom = index_feature.GetGeometryRef()
                            known_name_fields = ['Name', 'location', 'filename']
                            known_url_fields = ['url', 'URL']

                            if index_geom.Intersects(warp_region.export_as_geom()):
                                tile_name = None
                                tile_url = None
                                field_names = [field.name for field in index_layer.schema]

                                for f in known_name_fields:
                                    if f in field_names:
                                        tile_name = index_feature.GetField(f).strip()

                                    if tile_name is not None:
                                        break

                                for f in known_url_fields:
                                    if f in field_names:
                                        tile_url = index_feature.GetField(f).strip()

                                    if tile_url is not None:
                                        break

                                if tile_name is None or tile_url is None:
                                    utils.echo_warning_msg('could not parse index fields')
                                    continue

                                tile_url = '/'.join(tile_url.split('/')[:-1]) + '/' + tile_name.split('/')[-1]
                                ## add vertical datum to output;
                                ## field is NativeVdatum
                                ## must get from metadata
                                # if ept_infos is None:
                                #     this_epsg = vdatums.get_vdatum_by_name(
                                #         feature['attributes']['NativeVdatum']
                                #     )
                                # else:
                                #     #print(ept_infos['srs'])
                                #     ## horizontal datum is wrong in ept, most seem to be nad83
                                #     #this_epsg = 'epsg:{}+{}'.format(ept_infos['srs']['horizontal'], ept_infos['srs']['vertical'])
                                #     if 'vertical' in ept_infos['srs'].keys():
                                #         vertical_epsg = ept_infos['srs']['vertical']
                                #         horizontal_epsg = ept_infos['srs']['horizontal']
                                #         this_epsg = 'epsg:{}+{}'.format(
                                #             horizontal_epsg, vertical_epsg
                                #         )
                                #         #this_epsg = 'epsg:4269+{}'.format(ept_infos['srs']['vertical'])
                                #     else:
                                #         # try to extract the vertical datum from the wkt
                                #         #horizontal_epsg = ept_infos['srs']['horizontal']
                                #         this_wkt = ept_infos['srs']['wkt']
                                #         dst_horz, dst_vert = gdalfun.epsg_from_input(
                                #             this_wkt
                                #         )
                                #         this_epsg = '{}+{}'.format(
                                #             dst_horz, dst_vert
                                #         )

                                self.add_entry_to_results(
                                    tile_url,
                                    os.path.join(
                                        '{}/{}'.format(
                                            feature['attributes']['id'],
                                            tile_url.split('/')[-1]
                                        )
                                    ),
                                    feature['attributes']['data_type'],
                                    #this_epsg=this_epsg,
                                )

                        index_ds = index_layer = None
                        #if not self.keep_footprints:
                        utils.remove_glob(index_zipfile, *index_shps)
                        #utils.remove_glob(*index_shps)

                # # spatial_metadata
                # elif link['serviceID'] == 166 and self.datatype == 'sm': 
                #     self.add_entry_to_results(
                #         link['link'],
                #         os.path.join(
                #             '{}/{}'.format(
                #                 feature['attributes']['ID'],
                #                 link['link'].split('/')[-1]
                #             )
                #         ),
                #         link['label'],
                #         this_epsg=None
                #     )
 
        #self.results = [x for x in np.unique(self.results, axis=0)]
        return(self)
    
## Digital Coast - Data Access Viewer - SLR shortcut
class SLR(DAV):
    """Sea Level Rise DEMs via Digital Coast.

    < SLR >
    """
    
    def __init__(self, **kwargs):
        super().__init__(name='SLR', where='ID=6230', **kwargs)

        
## Digital Coast - Data Access Viewer CoNED shortcut
class CoNED(DAV):
    """Coastal NED (CoNED) DEMs via Digital Coast

    < CoNED >
    """
    
    def __init__(self, **kwargs):
        super().__init__(
            name='CoNED', where="title LIKE '%CoNED%'", **kwargs
        )

        
## Digital Coast - Data Access Viewer - CUDEM shortcut
class CUDEM(DAV):
    """CUDEM Tiled DEMs via Digital Coast

    < CUDEM:datatype=None >
    """
    
    def __init__(self, datatype = None, **kwargs):
        datatype = utils.str_or(datatype, 'all')

        if datatype == '19' or datatype.lower() == 'ninth':
            where="title LIKE '%CUDEM%Ninth%'"
        elif datatype == '13' or datatype.lower() == 'third':
            where="title LIKE '%CUDEM%Third%'"
        else:
            where="title LIKE '%CUDEM%'"
            
        super().__init__(name='CUDEM', where=where, **kwargs)

### End
