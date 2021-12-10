### tnm.py - USGS The National Map fetch
##
## Copyright (c) 2010 - 2021 CIRES Coastal DEM Team
##
## tnm.py is part of CUDEM
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
## The National Map (TNM) - USGS
##
## Fetch elevation data from The National Map
## NED, 3DEP, NHD, Etc.
##
## Various datasets from USGS's National Map. The National Map is a 
## collaborative effort among the USGS and other Federal, State, and local partners to improve
## and deliver topographic information for the Nation.
##
### Code:

import os
import sys

from cudem import utils
from cudem import regions
from cudem import datasets

import cudem.fetches.utils as f_utils
import cudem.fetches.FRED as FRED

class TheNationalMap(f_utils.FetchModule):
    """Fetch elevation data from The National Map"""

    def __init__(self, where=[], formats=None, extents=None, q=None, **kwargs):
        super().__init__(**kwargs)
        self._tnm_api_url = 'http://tnmaccess.nationalmap.gov/api/v1'
        self._tnm_dataset_url = 'https://tnmaccess.nationalmap.gov/api/v1/datasets?'
        self._tnm_product_url = 'https://tnmaccess.nationalmap.gov/api/v1/products?'
        self._tnm_meta_base = 'https://www.sciencebase.gov/catalog/item/'
        self._elev_ds = ['National Elevation Dataset (NED) 1 arc-second', 'Digital Elevation Model (DEM) 1 meter',
                         'National Elevation Dataset (NED) 1/3 arc-second', 'National Elevation Dataset (NED) 1/9 arc-second',
                         'National Elevation Dataset (NED) Alaska 2 arc-second', 'Alaska IFSAR 5 meter DEM',
                         'Original Product Resolution (OPR) Digital Elevation Model (DEM)', 'Ifsar Digital Surface Model (DSM)',
                         'Ifsar Orthorectified Radar Image (ORI)', 'Lidar Point Cloud (LPC)',
                         'National Hydrography Dataset Plus High Resolution (NHDPlus HR)', 'National Hydrography Dataset (NHD) Best Resolution',
                         'National Watershed Boundary Dataset (WBD)', 'USDA National Agriculture Imagery Program (NAIP)',
                         'Topobathymetric Lidar DEM', 'Topobathymetric Lidar Point Cloud']
        self._outdir = os.path.join(os.getcwd(), 'tnm')
        self.where = where
        
        self.name = 'tnm'
        self._urls = [self._tnm_api_url]
        
        self.FRED = FRED.FRED(name=self.name, verbose = self.verbose)
        self.update_if_not_in_FRED()

        self.formats = formats
        self.extents = extents
        self.q = q
        
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
        _req = f_utils.Fetch(self._tnm_dataset_url).fetch_req()
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
            else: datatype = 'tnm'
            url_enc = f_utils.urlencode({'datasets': ds['sbDatasetTag']})
            try:
                pubDate = ds['lastPublishedDate']
            except: pubDate = utils.this_year()
            try:
                metadataDate = ds['lastUpdatedDate']
            except: metadataDate = utils.this_year()            
            if geom is not None:
                self.FRED._add_survey(Name = ds['sbDatasetTag'], ID = ds['id'], Agency = 'USGS', Date = pubDate[-4:],
                                      MetadataLink = ds['infoUrl'], MetadataDate = metadataDate,
                                      DataLink = '{}{}'.format(self._tnm_product_url, url_enc), Link = ds['dataGovUrl'], Resolution = ','.join(ds['extents']),
                                      DataType = datatype, DataSource = 'tnm', HorizontalDatum = h_epsg,
                                      VerticalDatum = v_epsg, Etcetra = fmt, Info = ds['refreshCycle'], geom = geom)

    ## ==============================================
    ## update the FRED geojson with each TNM dataset
    ## each dataset will have any number of products, which get parsed for the data-link
    ## in _parse_results().
    ## ==============================================
    def update(self):
        """Update FRED with each dataset in TNM"""
        
        datasets = self._datasets()
        self.FRED._open_ds(1)
        if self.verbose:
            _prog = utils.CliProgress('scanning {} datasets from TNM...'.format(len(datasets)))
        for i, ds in enumerate(datasets):
            if self.verbose:
                _prog.update_perc((i, len(datasets)))
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
                    print(tag)
                    this_xml = FRED.iso_xml('{}?format=iso'.format(tag['infoUrl']))
                    geom = this_xml.bounds(geom=True)
                    h_epsg, v_epsg = this_xml.reference_system()
                    self._update_dataset(tag, fmt, geom, h_epsg, v_epsg)
            else:
                this_xml = FRED.iso_xml('{}?format=iso'.format(ds['infoUrl']))
                geom = this_xml.bounds(geom = True)
                h_epsg, v_epsg = this_xml.reference_system()
                self._update_dataset(ds, fmt, geom, h_epsg, v_epsg)
        if self.verbose:
            _prog.end(0, 'scanned {} datasets from TNM'.format(len(datasets)))
        self.FRED._close_ds()

    def run(self):#, f = None, e = None, q = None):
        """parse the tnm results from FRED"""
        
        e = self.extents.split(',') if self.extents is not None else None
        f = self.formats.split(',') if self.formats is not None else None
        q = self.q
        for surv in FRED._filter_FRED(self):
            #print(surv)
            offset = 0
            total = 0
            while True:
                _dataset_results = []
                _data = {'bbox': self.region.format('bbox'), 'max': 100, 'offset': offset}
                if q is not None: _data['q'] = str(q)
                if f is None:
                    _data['prodFormats'] = surv['Etcetra']
                else: _data['prodFormats'] = ','.join(f)
                if e is None: e = []

                _req = f_utils.Fetch(surv['DataLink']).fetch_req(params=_data)
                if _req is not None and _req.status_code == 200:
                    try:
                        _dataset_results = _req.json()
                        total = _dataset_results['total']
                    except ValueError:
                        utils.echo_error_msg('tnm server error, try again')
                    except Exception as e:
                        utils.echo_error_msg('error, {}'.format(e))
                
                if len(_dataset_results) > 0:
                    for item in _dataset_results['items']:
                        if _data['prodFormats'] is None:
                            fmts = []
                        else: fmts = _data['prodFormats'].split(',')
                        f_url = None
                        if len(e) > 0:
                            for extent in e:
                                if item['extent'] == extent:
                                    for fmt in fmts:
                                        if fmt in item['urls'].keys():
                                            f_url = item['urls'][fmt]
                                            break
                                    if f_url is None: f_url = item['downloadURL']
                                    self.results.append([f_url, f_url.split('/')[-1], surv['DataType']])
                        else:
                            for fmt in fmts:
                                if fmt in item['urls'].keys():
                                    f_url = item['urls'][fmt]
                                    break
                            if f_url is None:  f_url = item['downloadURL']
                            self.results.append([f_url, f_url.split('/')[-1], surv['DataType']])
                offset += 100
                if offset >= total: break
        return(self)
    
    ## ==============================================
    ## _update_prods() and _parse_prods_results() will update FRED with every product as a feature, rather than
    ## the default of each feature being a TNM dataset. _update_prods() takes much longer time to gather the
    ## products for each dataset and recording them in FRED, though the parsing of results is much faster.
    ## For our purposes, we wont be using most of what's available on TNM, so it is a bit of a waste to store
    ## all their datasets, which are already stored online, in FRED. This means user-time for fetches TNM is a
    ## bit slower, however storage costs are minimal and fewer updates may be necesary...
    ## ==============================================                
    def _update_prods(self):
        """updated FRED with each product file available from TNM"""
        
        for dsTag in self._elev_ds:
            offset = 0
            utils.echo_msg('processing TNM dataset {}...'.format(dsTag))
            _req = f_utils.Fetch(self._tnm_product_url).fetch_req(params={'max': 1, 'datasets': dsTag})
            try:
                _dsTag_results = _req.json()
            except ValueError:
                utils.echo_error_msg('tnm server error, try again')
            except Exception as e:
                utils.echo_error_msg('error, {}'.format(e))
                
            total = _dsTag_results['total']
            if self.verbose:
                _prog = utils.CliProgress('gathering {} products from {}...'.format(total, dsTag))
            
            ds = self._datasets(dataset = dsTag)
            #this_xml = FRED.iso_xml('{}{}?format=iso'.format(self._tnm_meta_base, ds['id']))
            this_xml = FRED.iso_xml('{}?format=iso'.format(ds['infoUrl']))
            h_epsg, v_epsg = this_xml.reference_system()
            
            while True:
                _data = {'max': 100, 'datasets': dsTag, 'offset': offset}
                _req = f_utils.Fetch(self._tnm_product_url).fetch_req(params = _data)
                try:
                    _dsTag_results = _req.json()
                except ValueError:
                    utils.echo_error_msg('tnm server error, try again')
                except Exception as e:
                    utils.echo_error_msg('error, {}'.format(e))
                if self.verbose:
                    _prog.update_perc((offset,total), msg = 'gathering {} products from {}...'.format(total, dsTag))
                
                for i, item in enumerate(_dsTag_results['items']):
                    if self.verbose:
                        _prog.update_perc((i+offset,total), msg = 'gathering {} products from {}...'.format(total, dsTag))
                    try:
                        self.FRED.layer.SetAttributeFilter("ID = '{}'".format(item['sourceId']))
                    except: pass
                    if self.FRED.layer is None or len(self.FRED.layer) == 0:
                        bbox = item['boundingBox']
                        geom = regions.Region().from_list([bbox['minX'], bbox['maxX'], bbox['minY'], bbox['maxY']]).export_as_geom()

                        if item['format'] == 'IMG' or item['format'] == 'GeoTIFF':
                            tnm_ds = 'raster'
                        elif item['format'] == 'LAZ' or item['format'] == 'LAS':
                            tnm_ds = 'lidar'
                        else: tnm_ds = 'tnm'

                        if geom is not None:
                            self.FRED._add_survey(Name = item['title'], ID = item['sourceId'], Agency = 'USGS', Date = item['publicationDate'],
                                                  MetadataLink = item['metaUrl'], MetadataDate = item['dateCreated'],
                                                  DataLink = item['downloadURL'], Link = item['sourceOriginId'], Resolution = item['extent'],
                                                  DataType = tnm_ds, DataSource = 'tnm', HorizontalDatum = h_epsg,
                                                  VerticalDatum = v_epsg, Etcetra = dsTag, Info = item['moreInfo'], geom = geom)
                offset += 100
                if total - offset <= 0: break
            if self.verbose:
                _prog.end(0, 'gathered {} products from {}'.format(total, dsTag))
                           
    def _parse_prods_results(self, r, f = None, e = None, q = None):
        for surv in FRED._filter_FRED(self):
            for d in surv['DataLink'].split(','):
                if d != '':
                    self.results.append([d, d.split('/')[-1], surv['DataType']])
        
    def yield_xyz(self, entry):
        """yield the xyz data from the tnm fetch module"""
        
        if f_utils.Fetch(entry[0], callback=self.callback, verbose=self.verbose).fetch_file(entry[1]) == 0:
            datatype = entry[-1]
            if datatype == 'raster':
                src_tnms = utils.p_unzip(entry[1], ['tif', 'img', 'gdal', 'asc', 'bag'])
                for src_tnm in src_tnms:
                    _ds = datasets.RasterFile(fn=src_tnm, data_format=200, epsg=4326, warp=self.warp,
                                              name=src_tnm, src_region=self.region, verbose=self.verbose)
                    for xyz in _ds.yield_xyz():
                        if xyz.z != 0:
                            yield(xyz)
                    utils.remove_glob(src_tnm)
        utils.remove_glob(entry[1])

## ==============================================
## class TNM is the old tnm api, doesn't currently work.
## held for now for historical purposes
## ==============================================
class TNM(f_utils.FetchModule):
    """Fetch elevation data from The National Map"""

    def __init__(self, index=False, ds=1, sub_ds=None, formats=None, extent=None, q=None, **kwargs):
        super().__init__(**kwargs)
        self._tnm_api_url = 'http://tnmaccess.nationalmap.gov/api/v1'
        self._tnm_dataset_url = 'https://tnmaccess.nationalmap.gov/api/v1/datasets?'
        self._tnm_product_url = 'https://tnmaccess.nationalmap.gov/api/v1/products?'
        self._outdir = os.path.join(os.getcwd(), 'tnm')        
        self._dataset_results = {}
        self._tnm_ds = []
        self._tnm_df = []
        self.name = 'tnm'

        self.index = index
        self.ds = ds
        self.sub_ds = sub_ds
        self.formats = formats
        self.extent = extent
        self.q = q

    def run(self):
        '''Run the TNM (National Map) fetching module.'''
        
        if self.region is None: return([])
        _req = f_utils.Fetch(self._tnm_dataset_url).fetch_req()
        if _req is not None:
            try:
                self._datasets = _req.json()
            except Exception as e:
                utils.echo_error_msg('try again, {}'.format(e))
                self.status = -1
        else: self.status = -1

        if self.status == 0:
            if self.index:
                self.print_dataset_index()
                sys.exit()

            this_ds = [int(self.ds)]
            if self.sub_ds is not None: this_ds.append(int(self.sub_ds))
            self._tnm_ds = [this_ds]
            self._tnm_df = [] if self.formats is None else [x.strip() for x in self.formats.split(',')]
            self._extents = [] if self.extent is None else [x.strip() for x in self.extent.split(',')]
            self.filter_datasets(self.q)
        return(self)

    def filter_datasets(self, q = None):
        """Filter TNM datasets."""
        
        utils.echo_msg('filtering TNM dataset results...')
        sbDTags = []        
        for ds in self._tnm_ds:
            dtags = self._datasets[ds[0]]['tags']
            if len(ds) > 1:
                if len(dtags) == 0:
                    sbDTags.append( self._datasets[ds[0]]['sbDatasetTag'])
                else:
                    dtag = list(dtags.keys())[ds[1]]
                    sbDTag = self._datasets[ds[0]]['tags'][dtag]['sbDatasetTag']
                    if len(self._tnm_df) == 0:
                        formats = self._datasets[ds[0]]['tags'][dtag]['formats']
                        self._tnm_df = formats
                    sbDTags.append(sbDTag)
            else:
                all_formats = True if len(self._tnm_df) == 0 else False
                if len(dtags) == 0:
                    sbDTags.append( self._datasets[ds[0]]['sbDatasetTag'])
                else:
                    print(dtags)
                    for dtag in list(dtags.keys()):
                        sbDTag = self._datasets[ds[0]]['tags'][dtag]['sbDatasetTag']
                        if all_formats:
                            formats = self._datasets[ds[0]]['tags'][dtag]['formats']
                            for ff in formats:
                                if ff not in self._tnm_df:
                                    if ff == 'FileGDB':
                                        self._tnm_df.append('FileGDB 10.1')
                                    self._tnm_df.append(ff)                            
                        sbDTags.append(sbDTag)
        self.data = {
            'bbox': regions.region_format(self.region, 'bbox'),
            'max': 10000,
        }
        if q is not None: self.data['q'] = str(q)
        self.data['datasets'] = ','.join(sbDTags)
        if len(self._tnm_df) > 0: self.data['prodFormats'] = ','.join(self._tnm_df)
        req = f_utils.Fetch(self._tnm_product_url).fetch_req(params=self.data)

        if req is not None:
            try:
                self._dataset_results = req.json()
            except ValueError:
                utils.echo_error_msg('tnm server error, try again')
            except Exception as e:
                utils.echo_error_msg('error, {}'.format(e))                
        else: self._status = -1

        if len(self._dataset_results) > 0:
            for item in self._dataset_results['items']:
                if len(self._extents) > 0:
                    for extent in self._extents:
                        if item['extent'] == extent:
                            #try:
                            f_url = item['downloadURL']
                            if item['format'] == 'IMG' or item['format'] == 'GeoTIFF':
                                tnm_ds = 'ned'
                            elif item['format'] == 'LAZ' or item['format'] == 'LAS':
                                tnm_ds = 'lidar'
                            else: tnm_ds = 'tnm'
                            self.results.append([f_url, f_url.split('/')[-1], tnm_ds])
                            #except: pass
                else:
                    #try:
                    f_url = item['downloadURL']
                    if item['format'] == 'IMG' or item['format'] == 'GeoTIFF':
                        tnm_ds = 'ned'
                    elif item['format'] == 'LAZ' or item['format'] == 'LAS':
                        tnm_ds = 'lidar'
                    else: tnm_ds = 'tnm'
                    self.results.append([f_url, f_url.split('/')[-1], tnm_ds])
                    #except: pass
        if self.verbose:
            utils.echo_msg('filtered \033[1m{}\033[m data files from TNM dataset results.'.format(len(self.results)))

    def print_dataset_index(self):
        for i,j in enumerate(self._datasets):
            try: fmts = ', '.join(j['formats'])
            except: fmts = ''
            try: exts = ', '.join(j['extents'])
            except: exts = ''
            print('%s: %s [ %s ]' %(i, j['title'], fmts))
            for m, n in enumerate(j['tags']):
                print('\t{}: {} [ {} ] [ {} ]'.format(m, n['title'], j['formats'][0]['value'], ', '.join(n['extents'])))

    def yield_xyz(self, entry):
        """yield the xyz data from the tnm fetch module"""
        
        if f_utils.Fetch(entry[0], callback=self.callback, verbose=self.verbose).fetch_file(entry[1]) == 0:
            datatype = entry[-1]
            if datatype == 'ned':
                src_tnm, src_zips = utils.procs_unzip(entry[1], ['tif', 'img', 'gdal', 'asc', 'bag'])

                _ds = datasets.RasterFile(
                    fn=src_tnm,
                    data_format=200,
                    epsg=4326,
                    warp=self.warp,
                    name=src_tnm,
                    src_region=self.region,
                    verbose=self.verbose
                )
                for xyz in _ds.yield_xyz():
                    if xyz.z != 0:
                        yield(xyz)
                        
                utils.remove_glob(src_tnm, *src_zips)
        utils.remove_glob(entry[1])
        
### End
