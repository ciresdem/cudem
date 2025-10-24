### tnm.py
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
from cudem import utils
from cudem.fetches import fetches
from cudem.fetches import FRED

class TheNationalMap(fetches.FetchModule):
    """The National Map:

    Fetch elevation data from The National Map
        
    Various datasets from USGS's National Map. The National Map is a 
    collaborative effort among the USGS and other Federal, State, and 
    local partners to improve and deliver topographic information for 
    the Nation.
    
    http://tnmaccess.nationalmap.gov/

    < tnm:datasets=None:formats=None:extents=None:q=None:date_type=None:date_start=None:date_end=None >
    """
   
    dataset_codes = [
        "National Boundary Dataset (NBD)",
        "National Elevation Dataset (NED) 1 arc-second",
        "Digital Elevation Model (DEM) 1 meter",
        "National Elevation Dataset (NED) 1/3 arc-second",
        "National Elevation Dataset (NED) 1/9 arc-second",
        "National Elevation Dataset (NED) Alaska 2 arc-second",
        "Alaska IFSAR 5 meter DEM",
        "National Elevation Dataset (NED) 1/3 arc-second - Contours",
        "Original Product Resolution (OPR) Digital Elevation Model (DEM)",
        "Ifsar Digital Surface Model (DSM)",
        "Ifsar Orthorectified Radar Image (ORI)",
        "Lidar Point Cloud (LPC)",
        "Historical Topographic Maps",
        "National Hydrography Dataset Plus High Resolution (NHDPlus HR)",
        "National Hydrography Dataset (NHD) Best Resolution",
        "National Watershed Boundary Dataset (WBD)",
        "Map Indices",
        "National Geographic Names Information System (GNIS)",
        "Small-scale Datasets - Boundaries",
        "Small-scale Datasets - Contours",
        "Small-scale Datasets - Hydrography",
        "Small-scale Datasets - Transportation",
        "National Structures Dataset (NSD)",
        "Combined Vector",
        "National Transportation Dataset (NTD)",
        "US Topo Current",
        "US Topo Historical",
        "Land Cover - Woodland",
        "3D Hydrography Program (3DHP)",
    ]

    format_keywords = [
        "ArcExport", "ArcGrid", "BIL", "FileGDB", "FileGDB 10.1", "FileGDB 10.2",
        "GeoPDF", "GeoTIFF", "GridFlow", "IMG", "JPEG2000", "LAS,LAZ", "NLAPS",
        "PDF", "SDE Export", "Shapefile", "Text", "TIFF", "TXT (pipes)",
    ]

    dataset_keywords = [
        "7.5", "Airports", "Bend", "Bridge", "Building", "Canal", "Cape", "Cave",
        "Cemetery", "Census", "Channel", "Chamber", "Church", "Civil", "Cliff",
        "Coast", "Coastline", "Conduit", "Contour", "Crossing", "Dam", "Dams",
        "Ditch", "Elevation", "Falls", "Flat", "Flume", "Forest", "Gaging", "Gaging Station",
        "Gap", "Gate", "Glacier", "Gut", "HUC", "Harbor", "Hospital", "Hydrography",
        "Hydrologic", "Images", "Intake", "Island", "Isthmus", "Lake", "Lakes", "Lava",
        "Levee", "Lidar", "Lock", "Map", "Maps", "Military", "Mine", "Oilfield", "Outflow",
        "Park", "Pillar", "Pipeline", "Plain", "Populated Place", "Post Office", "Range",
        "Rapids", "Reef", "Reserve", "Reservoir", "Ridge", "Rise", "River", "Rock",
        "School", "Sea", "Seep", "Shore", "Sink", "Slope", "Spring", "Stream", "Summit",
        "Swamp", "Tower", "Trail", "Tunnel", "Valley", "Water", "Waterfall", "Watershed",
        "Weir", "Well", "Woods",
    ]

    dataset_extents = [
        "10000 x 10000 meter", "1500 x 1500 meter", "15 x 15 minute", "1 x 1 degree",
        "1 x 2 degree", "1 x 3 degree", "1 x 4 degree", "2 x 1 degree", "30 x 30 minute",
        "30 x 60 minute", "3.75 minute x 3.75 minute", "3 x 3 degree", "7.5 x 15 minute",
        "7.5 x 7.5 minute", "Contiguous US", "HU-2 Region", "HU-4 Subregion", "HU-8 Subbasin",
        "National", "North America", "State", "Varies",
    ]

    date_types = ["dateCreated", "lastUpdated", "Publication"]

    __doc__ = '''{}\n dataset codes (datasets): \n{}\
    \n\n format keywords (formats): \n{}\
    \n\n dataset keywords (q): \n{}\
    \n\n dataset extents (extents): \n{}\
    \n\n date types (date_types): \n{}'''.format(
        __doc__,
        '\n'.join(['{}: {}'.format(i, x) for i, x in enumerate(dataset_codes)]),
        format_keywords, dataset_keywords, dataset_extents, date_types
    )

    
    def __init__(
            self, datasets=None, formats=None, extents=None, q=None,
            date_type=None, date_start=None, date_end=None, **kwargs
    ):
        super().__init__(name='tnm', **kwargs)
        self.q = q
        self.f = formats
        self.e = extents
        self.extents = extents
        self.formats = formats
        self.datasets = datasets
        self.date_type = date_type
        self.date_start = date_start
        self.date_end = date_end
        
        self._tnm_api_url = 'http://tnmaccess.nationalmap.gov/api/v1'
        self._tnm_api_products_url = 'http://tnmaccess.nationalmap.gov/api/v1/products?'
        #self.headers['Host'] = 	'tnmaccess.nationalmap.gov'

        
    def run(self):
        offset = 0
        total = 0
        while True:
            _data = {
                'bbox': self.region.format('bbox'),
                'max': 100,
                'offset': offset
            }
            if self.datasets is not None:
                datasets = self.datasets.split('/')
                try:
                    datasets = [self.dataset_codes[i] for i in [int(x) for x in datasets]]
                    _data['datasets'] =  ','.join(datasets)
                except Exception as e:
                    utils.echo_warning_msg(
                        f'could not parse datasets: {datasets}, {e}'
                    )
                    _data['datasets'] = "National Elevation Dataset (NED) 1 arc-second"

            if self.q is not None: _data['q'] = str(self.q)
            if self.f is not None: _data['prodFormats'] = ','.join(self.f.split('/'))
            if self.e is not None: _data['prodExtents'] = ','.join(self.e.split('/'))
            if self.date_start is not None:
                _data['start'] = self.date_start

                if self.date_end is not None:
                    _data['end'] = self.date_end
                else:
                    _data['end'] = datetime.datetime.now().strftime('%Y-%m-%d')

                if self.date_type is not None:
                    _data['dateType'] = self.date_type
                else:
                    _data['dateType'] = 'dateCreated'

            _req = fetches.Fetch(
                self._tnm_api_products_url, verbose=self.verbose
            ).fetch_req(params=_data, timeout=60, read_timeout=60)
            if _req is not None and _req.status_code == 200:
                utils.echo_msg(_req.url)
                response_text = _req.text
                if response_text.startswith("{errorMessage"):
                    continue
                
                features = _req.json()                
                if 'total' in features.keys():
                    total = features['total']
                    for feature in features['items']:
                        self.add_entry_to_results(
                            feature['downloadURL'],
                            feature['downloadURL'].split('/')[-1],
                            feature['format']
                        )

            offset += 100
            if offset >= total:
                break                

            
## The National Map - NED (1 & 1/3) shortcut
class NED(TheNationalMap):
    """National Elevation Dataset (NED) via The National Map (TNM)

    < NED >
    """
    
    def __init__(self, **kwargs):
        super().__init__(datasets='1/3/4/5', **kwargs)
        self.data_format = 200

        
## The National Map - NED (1m) shortcut
class NED1(TheNationalMap):
    """National Elevation Dataset (NED) (1 meter) via The National Map (TNM)

    < NED1 >
    """
    
    def __init__(self, **kwargs):
        super().__init__(datasets='2', **kwargs)
        self.data_format = 200

        
class TNM_LAZ(TheNationalMap):
    def __init__(self, **kwargs):
        super().__init__(formats="LAZ", **kwargs)
        self.data_format = 300

        
## The National Map
## update is broken! fix this.
class TheNationalMapOLD(fetches.FetchModule):
    """USGS' The National Map

    Fetch elevation data from The National Map
        
    Various datasets from USGS's National Map. The National Map is a 
    collaborative effort among the USGS and other Federal, State, and 
    local partners to improve and deliver topographic information for 
    the Nation.
    
    http://tnmaccess.nationalmap.gov/

    < tnm:formats=None:extents=None:q=None >
    """

    def __init__(self, where=[], formats=None, extents=None, q=None, **kwargs):
        super().__init__(name='tnm', **kwargs)
        self.where = [where] if len(where) > 0 else []        
        self.formats = formats
        self.extents = extents
        self.q = q

        ## various TNM URLs
        self._tnm_api_url = 'http://tnmaccess.nationalmap.gov/api/v1'
        self._tnm_dataset_url = 'https://tnmaccess.nationalmap.gov/api/v1/datasets?'
        self._tnm_product_url = 'https://tnmaccess.nationalmap.gov/api/v1/products?'
        self._tnm_meta_base = 'https://www.sciencebase.gov/catalog/item/'
        self._urls = [self._tnm_api_url]
        
        ## The relevant TNM datasets
        self._elev_ds = ['National Elevation Dataset (NED) 1 arc-second',
                         'Digital Elevation Model (DEM) 1 meter',
                         'National Elevation Dataset (NED) 1/3 arc-second',
                         'National Elevation Dataset (NED) 1/9 arc-second',
                         'National Elevation Dataset (NED) Alaska 2 arc-second',
                         'Alaska IFSAR 5 meter DEM',
                         'Original Product Resolution (OPR) Digital Elevation Model (DEM)',
                         'Ifsar Digital Surface Model (DSM)',
                         'Ifsar Orthorectified Radar Image (ORI)',
                         'Lidar Point Cloud (LPC)',
                         'National Hydrography Dataset Plus High Resolution (NHDPlus HR)',
                         'National Hydrography Dataset (NHD) Best Resolution',
                         'National Watershed Boundary Dataset (WBD)',
                         'USDA National Agriculture Imagery Program (NAIP)',
                         'Topobathymetric Lidar DEM',
                         'Topobathymetric Lidar Point Cloud']

        ## TNM is in FRED, set that up here
        self.FRED = FRED.FRED(name=self.name, verbose = self.verbose)
        self.update_if_not_in_FRED()

        
    def update_if_not_in_FRED(self):
        self.FRED._open_ds()
        self.FRED._attribute_filter(["DataSource = '{}'".format(self.name)])
        if len(self.FRED.layer) == 0:
            self.FRED._close_ds()
            self.update()
        self.FRED._close_ds()

        
    def _fred_datasets(self):
        ids = {}
        for r in self.FRED._filter(self.region, self.where, ['tnm']):
            ids[r['ID']] = r['DataType']
        return(ids)

    
    def _datasets(self, dataset = None):
        _req = fetches.Fetch(self._tnm_dataset_url).fetch_req()
        if _req is not None and _req.status_code == 200:
            try:
                _datasets = _req.json()
                if dataset is not None:
                    for ds in _datasets:
                        tags = ds['tags']
                        if len(tags) > 0:
                            for t in tags:
                                if dataset == t['sbDatasetTag']:
                                    _datasets = t
                                    break
                        else:
                            if dataset == ds['sbDatasetTag']:
                                _datasets = ds
                                break
            except Exception as e:
                utils.echo_error_msg('try again, {}'.format(e))
        else: _datasets = None
        return(_datasets)

    
    def _dataset_tags():
        tnm_ds = self._datasets()
        dsTags = []
        for ds in tnm_ds:
            tags = ds['tags']
            if len(tags) > 0:
                for t in tags:
                    dsTags.append(t['sbDatasetTag'])
            else: dsTags.append(ds['sbDatasetTag'])
        return(dsTags)

    
    def _update_dataset(self, ds, fmt, geom, h_epsg, v_epsg):
        self.FRED._attribute_filter(["ID = '{}'".format(ds['id'])])
        if self.FRED.layer is None or len(self.FRED.layer) == 0:
            if 'IMG' in fmt or 'TIFF' in fmt:
                datatype = 'raster'
            elif 'LAS' in fmt or 'LAZ' in fmt:
                datatype = 'lidar'
            else:
                datatype = 'tnm'
                
            url_enc = urlencode({'datasets': ds['sbDatasetTag']})
            try:
                pubDate = ds['lastPublishedDate']
            except:
                pubDate = utils.this_year()
                
            try:
                metadataDate = ds['lastUpdatedDate']
            except:
                metadataDate = utils.this_year()
                
            if geom is not None:
                self.FRED._add_survey(
                    Name = ds['sbDatasetTag'],
                    ID = ds['id'],
                    Agency = 'USGS',
                    Date = pubDate[-4:],
                    MetadataLink = ds['infoUrl'],
                    MetadataDate = metadataDate,
                    DataLink = '{}{}'.format(self._tnm_product_url, url_enc),
                    Link = ds['dataGovUrl'],
                    Resolution = ','.join(ds['extents']),
                    DataType = datatype,
                    DataSource = 'tnm',
                    HorizontalDatum = None,#h_epsg,
                    VerticalDatum = None,#v_epsg,
                    Etcetra = fmt,
                    Info = ds['refreshCycle'],
                    geom = geom
                )

                
    ## update the FRED geojson with each TNM dataset
    ## each dataset will have any number of products, which get parsed for the data-link
    ## in _parse_results().
    def update(self):
        """Update FRED with each dataset in TNM"""
        
        datasets = self._datasets()
        self.FRED._open_ds(1)
        with tqdm(
                total=len(datasets),
                desc='scanning TNM datasets',
                leave=self.verbose
        ) as pbar:
            for i, ds in enumerate(datasets):
                pbar.update(1)
                for fmt in ds['formats']:
                    if 'isDefault' in fmt.keys():
                        fmt = fmt['value']
                        break
                    #print(ds)
                #print(len(ds))
                #this_xml = FRED.iso_xml('{}{}?format=iso'.format(self._tnm_meta_base, ds['id']))

                tags = ds['tags']
                if len(tags) > 0:
                    for tag in tags:
                        #utils.echo_msg(tag)
                        #this_xml = iso_xml('{}?format=iso'.format(tag['infoUrl']))
                        info_url = tag['infoUrl']
                        if info_url == '':
                            continue
                        #    info_url = tag['dataGovUrl']
                        
                        #utils.echo_msg(info_url)
                        _req = fetches.Fetch(info_url+'?format=json').fetch_req()
                        if _req is not None and _req.status_code == 200:
                            try:
                                _results = _req.json()
                                utils.echo_msg(_results.keys())
                                bbox = _results['spatial']['boundingBox']
                                geom =  regions.Region().from_list(
                                    [float(bbox['minX']), float(bbox['maxX']),
                                     float(bbox['minY']), float(bbox['maxY'])]
                                ).export_as_geom()
                                h_epsg = None
                                v_epsg = None
                                #print(_results)
                                #this_xml = iso_xml('{}?format=atom'.format(info_url))
                                #geom = this_xml.bounds(geom=True)
                                #h_epsg, v_epsg = this_xml.reference_system()                                
                                self._update_dataset(tag, fmt, geom, h_epsg, v_epsg)
                            except:
                                utils.echo_warning_msg(tag)
                else:
                    info_url = ds['infoUrl']
                    if info_url == '':
                        continue
                    #if info_url == '':
                    #    info_url = ds['dataGovUrl']

                    _req = fetches.Fetch(info_url+'?format=json').fetch_req()
                    if _req is not None and _req.status_code == 200:
                        try:
                            _results = _req.json()
                            bbox = _results['spatial']['boundingBox']
                            geom =  regions.Region().from_list(
                                [float(bbox['minX']), float(bbox['maxX']),
                                 float(bbox['minY']), float(bbox['maxY'])]
                            ).export_as_geom()
                            h_epsg = None
                            v_epsg = None
                            #this_xml = iso_xml('{}?format=iso'.format(ds['infoUrl']))
                            #this_xml = iso_xml('{}?format=atom'.format(info_url))
                            #geom = this_xml.bounds(geom = True)
                            #h_epsg, v_epsg = this_xml.reference_system()
                            self._update_dataset(ds, fmt, geom, h_epsg, v_epsg)
                        except:
                            utils.echo_warning_msg(tag)
        self.FRED._close_ds()

        
    def run(self):#, f = None, e = None, q = None):
        """parse the tnm results from FRED"""
        
        e = self.extents.split(',') if self.extents is not None else None
        f = self.formats.split(',') if self.formats is not None else None
        q = self.q
        _results = FRED._filter_FRED(self)
        with tqdm(
                total=len(_results),
                desc='scanning for TNM datasets',
                leave=self.verbose
        ) as pbar:
            for surv in _results:
                offset = 0
                total = 0
                pbar.update(1)
                while True:
                    _dataset_results = []
                    _data = {
                        'bbox': self.region.format('bbox'),
                        'max': 100,
                        'offset': offset
                    }
                    if q is not None: _data['q'] = str(q)
                    if f is None:
                        _data['prodFormats'] = surv['Etcetra']
                    else:
                        _data['prodFormats'] = ','.join(f)

                    if e is None:
                        e = []

                    _req = fetches.Fetch(surv['DataLink']).fetch_req(params=_data)
                    if _req is not None and _req.status_code == 200:
                        try:
                            _dataset_results = _req.json()
                            total = _dataset_results['total']
                        except ValueError:
                            utils.echo_error_msg(
                                'tnm server error resulting in {}, try again'.format(e)
                            )
                        except Exception as e:
                            utils.echo_error_msg('error, {}'.format(e))

                    if len(_dataset_results) > 0:
                        #utils.echo_msg(_dataset_results)
                        for item in _dataset_results['items']:
                            p_dir = '_'.join(item['title'].split(' '))
                            if _data['prodFormats'] is None:
                                fmts = []
                            else:
                                fmts = _data['prodFormats'].split(',')

                            f_url = None
                            if len(e) > 0:
                                for extent in e:
                                    if item['extent'] == extent:
                                        for fmt in fmts:
                                            if fmt in item['urls'].keys():
                                                f_url = item['urls'][fmt]
                                                break

                                        if f_url is None:
                                            f_url = item['downloadURL']

                                        self.add_entry_to_results(
                                            f_url, f_url.split('/')[-1], surv['DataType']
                                        )
                            else:
                                for fmt in fmts:
                                    if fmt in item['urls'].keys():
                                        f_url = item['urls'][fmt]
                                        break

                                if f_url is None:
                                    f_url = item['downloadURL']

                                self.add_entry_to_results(
                                    f_url, f_url.split('/')[-1], surv['DataType']
                                )

                    offset += 100
                    if offset >= total:
                        break
                
        return(self)

    
    ## _update_prods() and _parse_prods_results() will update FRED with every
    ## product as a feature, rather than the default of each feature being a
    ## TNM dataset. _update_prods() takes much longer time to gather the
    ## products for each dataset and recording them in FRED, though the parsing
    ## of results is much faster.
    ## For our purposes, we wont be using most of what's available on TNM, so
    ## it is a bit of a waste to store all their datasets, which are already
    ## stored online, in FRED. This means user-time for fetches TNM is a
    ## bit slower, however storage costs are minimal and fewer updates may be
    ## necesary...
    def _update_prods(self):
        """updated FRED with each product file available from TNM"""
        
        for dsTag in self._elev_ds:
            offset = 0
            utils.echo_msg('processing TNM dataset {}...'.format(dsTag))
            _req = fetches.Fetch(
                self._tnm_product_url
            ).fetch_req(params={'max': 1, 'datasets': dsTag})
            try:
                _dsTag_results = _req.json()
            except ValueError:
                utils.echo_error_msg('tnm server error, try again')
            except Exception as e:
                utils.echo_error_msg('error, {}'.format(e))
                
            total = _dsTag_results['total']
            ds = self._datasets(dataset = dsTag)
            #this_xml = iso_xml('{}{}?format=iso'.format(self._tnm_meta_base, ds['id']))
            this_xml = iso_xml('{}?format=iso'.format(ds['infoUrl']))
            h_epsg, v_epsg = this_xml.reference_system()
            ##offset for prog
            while True:
                _data = {'max': 100, 'datasets': dsTag, 'offset': offset}
                _req = fetches.Fetch(self._tnm_product_url).fetch_req(params = _data)
                try:
                    _dsTag_results = _req.json()
                except ValueError:
                    utils.echo_error_msg('tnm server error, try again')
                except Exception as e:
                    utils.echo_error_msg('error, {}'.format(e))
                
                for i, item in enumerate(_dsTag_results['items']):
                    if self.verbose:
                        _prog.update_perc(
                            (i+offset,total),
                            msg='gathering {} products from {}...'.format(
                                total, dsTag
                            )
                        )
                    try:
                        self.FRED.layer.SetAttributeFilter(
                            "ID = '{}'".format(item['sourceId'])
                        )
                    except: pass
                    if self.FRED.layer is None or len(self.FRED.layer) == 0:
                        bbox = item['boundingBox']
                        geom = regions.Region().from_list(
                            [bbox['minX'], bbox['maxX'],
                             bbox['minY'], bbox['maxY']]
                        ).export_as_geom()

                        if item['format'] == 'IMG' or item['format'] == 'GeoTIFF':
                            tnm_ds = 'raster'
                        elif item['format'] == 'LAZ' or item['format'] == 'LAS':
                            tnm_ds = 'lidar'
                        else: tnm_ds = 'tnm'

                        if geom is not None:
                            self.FRED._add_survey(
                                Name = item['title'],
                                ID = item['sourceId'],
                                Agency = 'USGS',
                                Date = item['publicationDate'],
                                MetadataLink = item['metaUrl'],
                                MetadataDate = item['dateCreated'],
                                DataLink = item['downloadURL'],
                                Link = item['sourceOriginId'],
                                Resolution = item['extent'],
                                DataType = tnm_ds,
                                DataSource = 'tnm',
                                HorizontalDatum = h_epsg,
                                VerticalDatum = v_epsg,
                                Etcetra = dsTag,
                                Info = item['moreInfo'],
                                geom = geom
                            )
                offset += 100
                if total - offset <= 0: break

                
    def _parse_prods_results(self, r, f=None, e=None, q=None):
        for surv in FRED._filter_FRED(self):
            for d in surv['DataLink'].split(','):
                if d != '':
                    self.add_entry_to_results(
                        d,
                        os.path.join(self._outdir, d.split('/')[-1]),
                        surv['DataType']
                    )

                    
## The National Map - NED (1 & 1/3) shortcut
class NEDOLD(TheNationalMap):
    """National Elevation Dataset (NED) via The National Map (TNM)

    < NED >
    """
    
    def __init__(self, **kwargs):
        super().__init__(where="NAME LIKE '%NED%'", **kwargs)
        self.data_format = 200

        
## The National Map - NED (1m) shortcut
class NED1OLD(TheNationalMap):
    """National Elevation Dataset (NED) (1 meter) via The National Map (TNM)

    < NED1 >
    """
    
    def __init__(self, **kwargs):
        super().__init__(where="NAME LIKE '%DEM%'", **kwargs)
        self.data_format = 200

        
class TNM_LAZOLD(TheNationalMap):
    def __init__(self, **kwargs):
        super().__init__(formats="LAZ", **kwargs)
        self.data_format = 300


### End
