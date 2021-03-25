### fetches.py
##
## Copyright (c) 2012 - 2021 CIRES Coastal DEM Team
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
## current fetch modules: dc, nos, mb, charts, usace, srtm_cgiar, srtm_plus, tnm, gmrt, emodnet, mar_grav, osm
##
### Code:

import os
import sys
import time
import copy

## ==============================================
## networking/etc
## ==============================================
import requests
import ftplib
import urllib
import lxml.html as lh
import lxml.etree
import json

## ==============================================
## queues and threading
## ==============================================
import threading
try:
    import Queue as queue
except: import queue as queue

## ==============================================
## import gdal/numpy
## ==============================================
import ogr
import gdal

## ==============================================
## import CUDEM
## ==============================================
from cudem import utils
from cudem import regions
from cudem import dlim

__version__ = '0.7.0'

## =============================================================================
##
## Fetching Functions
##
## Generic fetching and processing functions, etc.
##
## The `FetchResults` class will fetch a list of fetch results [[url, file-name, data-type]...]
## in a queue `fetch_queue` using 3 threads; set `p` and `s` as a fetch module object to processes
## and dump XYZ data from the fetched results, respectively. Use `fetch_file` to fetch single files.
##
## =============================================================================
r_headers = { 'User-Agent': 'GeoMods: Fetches v%s' %(__version__) }

def fetch_queue(q, fo, fg):
    """fetch queue `q` of fetch results"""
    
    while True:
        fetch_args = q.get()
        this_region = fetch_args[2]
        this_dt = fetch_args[4].lower()
        fetch_args[2] = None
        if not fetch_args[3]():
            if not fg['proc_p'] and not fg['dump_p']:
                if fetch_args[0].split(':')[0] == 'ftp':
                    fetch_ftp_file(*tuple(fetch_args))
                else: fetch_file(*tuple(fetch_args))
            else:
                if not fg['dump_p']:
                    if not os.path.exists(os.path.dirname(fetch_args[1])):
                        try:
                            os.makedirs(os.path.dirname(fetch_args[1]))
                        except: pass
                    #o_x_fn = '.'.join(fetch_args[1].split('.')[:-1]) + '.xyz'
                    o_x_fn = fetch_args[1] + '.xyz'
                    utils.echo_msg('processing local file: {}'.format(o_x_fn))
                    if not os.path.exists(o_x_fn):
                        with open(o_x_fn, 'w') as out_xyz:
                            fetch_dump_xyz(fetch_args, module=fo, epsg=4326,
                                           z_region=fg['z_region'], inc=fg['inc'],
                                           dst_port=out_xyz)
                        if os.stat(o_x_fn).st_size == 0: utils.remove_glob(o_x_fn)
                else: fetch_dump_xyz(fetch_args, module=fo, epsg=4326)
        q.task_done()

def fetch_ftp_file(src_url, dst_fn, params=None, callback=None,
                   datatype=None, overwrite=False, verbose=False):
    """fetch an ftp file via urllib"""
    
    status = 0
    f = None
    halt = callback
    
    if verbose: utils.echo_msg('fetching remote ftp file: {}...'.format(src_url[:20]))
    if not os.path.exists(os.path.dirname(dst_fn)):
        try:
            os.makedirs(os.path.dirname(dst_fn))
        except: pass 
    try:
        f = urllib.request.urlopen(src_url)
    except: status - 1
    
    if f is not None:
        with open(dst_fn, 'wb') as local_file:
             local_file.write(f.read())
    if verbose: utils.echo_msg('fetched remote ftp file: {}.'.format(os.path.basename(src_url)))
    return(status)

def fetch_file(src_url, dst_fn, params=None, callback=lambda: False, datatype=None,
               overwrite=False, verbose=False, timeout=140, read_timeout=320):
    """fetch src_url and save to dst_fn"""
    
    status = 0
    req = None
    halt = callback

    if verbose: progress = utils.CliProgress('fetching remote file: {}...'.format(os.path.basename(src_url)[:20]))
    if not os.path.exists(os.path.dirname(dst_fn)):
        try:
            os.makedirs(os.path.dirname(dst_fn))
        except: pass 
    if not os.path.exists(dst_fn) or overwrite:
        try:
            with requests.get(src_url, stream=True, params=params, headers=r_headers, timeout=(timeout, read_timeout)) as req:
                req_h = req.headers
                if req.status_code == 200:
                    curr_chunk = 0
                    with open(dst_fn, 'wb') as local_file:
                        for chunk in req.iter_content(chunk_size=8196):
                            if halt(): break
                            if verbose: progress.update()
                            if chunk: local_file.write(chunk)
                else: utils.echo_error_msg('server returned: {}'.format(req.status_code))
        except Exception as e:
            utils.echo_error_msg(e)
            status = -1
    if not os.path.exists(dst_fn) or os.stat(dst_fn).st_size ==  0: status = -1
    if verbose: progress.end(status, 'fetched remote file: {}.'.format(os.path.basename(dst_fn)[:20]))
    return(status)

def fetch_req(src_url, params=None, tries=5, timeout=2, read_timeout=10):
    """fetch src_url and return the requests object"""
    
    if tries <= 0:
        utils.echo_error_msg('max-tries exhausted')
        return(None)
    try:
        return(requests.get(src_url, stream=True, params=params, timeout=(timeout, read_timeout), headers=r_headers))
    except: return(fetch_req(src_url, params=params, tries=(tries-1), timeout=(timeout+1), read_timeout=(read_timeout+10)))

def fetch_nos_xml(src_url, timeout=2, read_timeout=10, verbose=False):
    """fetch src_url and return it as an XML object"""
    
    results = lxml.etree.fromstring('<?xml version="1.0"?><!DOCTYPE _[<!ELEMENT _ EMPTY>]><_/>'.encode('utf-8'))
    try:
        req = fetch_req(src_url, timeout=timeout, read_timeout=read_timeout)
        results = lxml.etree.fromstring(req.text.encode('utf-8'))
    except:
        if verbose: utils.echo_error_msg('could not access {}'.format(src_url))
    return(results)
        
def fetch_html(src_url):
    """fetch src_url and return it as an HTML object"""
    
    req = fetch_req(src_url, timeout=2)
    if req:
        return(lh.document_fromstring(req.text))
    else: return(None)

class FetchResults(threading.Thread):
    """fetch results gathered from a fetch module.

    results is a list of URLs with data type
    """
    
    def __init__(self, results, region, out_dir, fetch_obj=None, fg=None, callback=lambda: False):
        threading.Thread.__init__(self)
        self.fetch_q = queue.Queue()
        self.results = results
        self.region = region
        self._outdir = out_dir
        self.stop_threads = callback
        self.fetch_obj = fetch_obj
        self.fg = fg
        
    def run(self):
        for _ in range(3):
            t = threading.Thread(target=fetch_queue, args=(self.fetch_q, self.fetch_obj, self.fg))
            t.daemon = True
            t.start()
        for row in self.results:
            if self.fg['list_p']: print(row[0])
            if self.fg['fetch_p']: self.fetch_q.put([row[0], os.path.join(self._outdir, row[1]), self.region, self.stop_threads, row[2], False, True])
        self.fetch_q.join()
        
## =============================================================================
##
## XML Metadata parsing
##
## =============================================================================
def xml2py(node):
    """parse an xml file into a python dictionary"""
    
    texts = {}
    if node is None: return(None)
    for child in list(node):
        child_key = lxml.etree.QName(child).localname
        if 'name' in child.attrib.keys(): child_key = child.attrib['name']
        if '{http://www.w3.org/1999/xlink}href' in child.attrib.keys():
            href = child.attrib['{http://www.w3.org/1999/xlink}href']
        else: href = None
        if child.text is None or child.text.strip() == '':
            if href is not None:
                if child_key in texts.keys():
                    texts[child_key].append(href)
                else: texts[child_key] = [href]
            else:
                if child_key in texts.keys():
                    ck = xml2py(child)
                    texts[child_key][list(ck.keys())[0]].update(ck[list(ck.keys())[0]])
                else: texts[child_key] = xml2py(child)
        else:
            if child_key in texts.keys():
                texts[child_key].append(child.text)
            else: texts[child_key] = [child.text]
    return(texts)

class IsoXml:
    def __init__(self, xml_url, timeout=2, read_timeout=10):
        self.url = xml_url
        self.xml_doc = self.fetch(timeout=timeout, read_timeout=read_timeout)
        self.namespaces = {
            'gmd': 'http://www.isotc211.org/2005/gmd', 
            'gmi': 'http://www.isotc211.org/2005/gmi', 
            'gco': 'http://www.isotc211.org/2005/gco',
            'gml': 'http://www.isotc211.org/2005/gml',
            'th': 'http://www.unidata.ucar.edu/namespaces/thredds/InvCatalog/v1.0',
            'wms': 'http://www.opengis.net/wms',
        }
        
    def fetch(self, timeout=2, read_timeout=10):
        
        return(fetch_nos_xml(self.url, timeout=timeout, read_timeout=read_timeout))

    def title(self):
        
        t = self.xml_doc.find('.//gmd:MD_DataIdentification/gmd:citation/gmd:CI_Citation/gmd:title/gco:CharacterString',
                              namespaces=self.namespaces)
        return(t.text if t is not None else 'Unknown')
        
    def bounds(self, geom = True):
        
        wl = self.xml_doc.find('.//gmd:westBoundLongitude/gco:Decimal', namespaces = self.namespaces)
        el = self.xml_doc.find('.//gmd:eastBoundLongitude/gco:Decimal', namespaces = self.namespaces)
        sl = self.xml_doc.find('.//gmd:southBoundLatitude/gco:Decimal', namespaces = self.namespaces)                            
        nl = self.xml_doc.find('.//gmd:northBoundLatitude/gco:Decimal', namespaces = self.namespaces)                           
        if wl is not None and el is not None and sl is not None and nl is not None:
            region = [float(wl.text), float(el.text), float(sl.text), float(nl.text)]
            if geom: return(regions.region2geom([float(wl.text), float(el.text), float(sl.text), float(nl.text)]))
            else: return(region)
        else: return(None)

    def polygon(self, geom=True):
        
        opoly = []
        polygon = self.xml_doc.find('.//{*}Polygon', namespaces=self.namespaces)
        if polygon is not None:
            nodes = polygon.findall('.//{*}pos', namespaces = self.namespaces)
            [opoly.append([float(x) for x in node.text.split()]) for node in nodes]
            if geom: return(regions.wkt2geom(utils.create_wkt_polygon(opoly)))
            else: return(opoly)
        else: return(None)
        
    def date(self):
        
        dt = self.xml_doc.find('.//gmd:date/gco:Date', namespaces=self.namespaces)
        if dt is None:
            dt = self.xml_doc.find('.//gmd:MD_DataIdentification/gmd:citation/gmd:CI_Citation/gmd:date/gmd:CI_Date/gmd:date/gco:Date', namespaces=self.namespaces)
        return(dt.text[:4] if dt is not None else '0000')

    def xml_date(self):
        
        mddate = self.xml_doc.find('.//gmd:dateStamp/gco:DateTime', namespaces=self.namespaces)
        return(utils.this_date() if mddate is None else mddate.text)
        
    def reference_system(self):
        
        ref_s = self.xml_doc.findall('.//gmd:MD_ReferenceSystem', namespaces=self.namespaces)
        if ref_s is None or len(ref_s) == 0: return(None, None)
        h_epsg = ref_s[0].find('.//gmd:code/gco:CharacterString', namespaces=self.namespaces)
        if h_epsg is not None: h_epsg = h_epsg.text.split(':')[-1]
        if len(ref_s) > 1:
            v_epsg = ref_s[1].find('.//gmd:code/gco:CharacterString', namespaces=self.namespaces)
            if v_epsg is not None: v_epsg = v_epsg.text.split(':')[-1]
        else: v_epsg = None
            
        return(h_epsg, v_epsg)

    def abstract(self):
        
        try:
            abstract = self.xml_doc.find('.//gmd:abstract/gco:CharacterString', namespaces=self.namespaces)
            abstract = '' if abstract is None else abstract.text
        except: abstract = ''
        return(abstract)

    def linkages(self):
        
        linkage = self.xml_doc.find('.//{*}linkage/{*}URL', namespaces=self.namespaces)
        if linkage is not None: linkage = linkage.text
        return(linkage)
    
    def data_links(self):
        
        dd = {}        
        dfs = self.xml_doc.findall('.//gmd:MD_Format/gmd:name/gco:CharacterString', namespaces=self.namespaces)
        dus = self.xml_doc.findall('.//gmd:onLine/gmd:CI_OnlineResource/gmd:linkage/gmd:URL', namespaces=self.namespaces)
        if dfs is not None:
            for i,j in enumerate(dfs):
                dd[j.text] = dus[i].text
        return(dd)

class WCS:
    
    def __init__(self, url):
        
        self.url = url
        self.namespaces = {
            'wms': 'http://www.opengis.net/wms', 'wcs': 'http://www.opengis.net/wcs/2.0',
            'ows': 'http://www.opengis.net/ows/2.0', 'gml': 'http://www.opengis.net/gml/3.2',
            'gmlcov': 'http://www.opengis.net/gmlcov/1.0'}
        self._get_capabilities()
        self._s_version = self._si()['ServiceTypeVersion'][0]

    def _get_capabilities(self):
        
        _data = {'request': 'GetCapabilities', 'service': 'WCS'}
        c = fetch_req(self.url, params = _data)
        cx = lxml.etree.fromstring(c.text.encode('utf-8'))
        self.service_provider = cx.find('.//ows:ServiceProvider', namespaces=self.namespaces)
        self.service_identification = cx.find('.//ows:ServiceIdentification', namespaces=self.namespaces)
        self.operations_metadata = cx.find('.//ows:OperationsMetadata', namespaces=self.namespaces)
        self.service_metadata = cx.find('.//wcs:ServiceMetadata', namespaces=self.namespaces)
        self.contents = cx.find('.//wcs:Contents', namespaces=self.namespaces)

    def _contents(self):
        
        c = []
        for coverage in self.contents.xpath('//wcs:CoverageSummary', namespaces=self.namespaces):
            c.append(xml2py(coverage))
        return(c)

    def _om(self):
        
        return(xml2py(self.operations_metadata))

    def _sp(self):
        
        return(xml2py(self.service_provider))
    
    def _si(self):
        
        return(xml2py(self.service_identification))
    
    def fix_coverage_id(self, coverage):
        
        return(':'.join(coverage.split('__')))

    def unfix_coverage_id(self, coverage):
        
        return('__'.join(coverage.split(':')))

    def _describe_coverage(self, coverage):
        
        c_d = {}
        valid = False
        c = self._contents()
        for cc in c:
            if coverage == cc['CoverageId']:
                valid = True
                c_d = cc
                break

        om = self._om()
        url = om['DescribeCoverage']['DCP']['HTTP']['Get'][0]
        _data = {'request': 'DescribeCoverage', 'service': 'WCS',
            'version': self._s_version, 'CoverageID': self.unfix_coverage_id(coverage)}
        d = fetch_req(url, params = _data)
        d_r = lxml.etree.fromstring(d.text.encode('utf-8'))
        cd = d_r.find('.//wcs:CoverageDescription', namespaces = self.namespaces)
        return(xml2py(d_r.find('.//wcs:CoverageDescription', namespaces = self.namespaces)))

    def _get_coverage_region(self, cov_desc):
        
        uc = [float(x) for x in cov_desc["boundedBy"]["Envelope"]["upperCorner"][0].split()]
        lc = [float(x) for x in cov_desc["boundedBy"]["Envelope"]["lowerCorner"][0].split()]
        return([lc[1], uc[1], lc[0], uc[0]])
    
    def _get_coverage_url(self, coverage, region = None):
        
        dl_coverage = self.fix_coverage_id(coverage)
        cov_desc = self._describe_coverage(coverage)
        fmt = cov_desc["ServiceParameters"]["nativeFormat"][0]        
        hl = [float(x) for x in cov_desc["domainSet"]["RectifiedGrid"]["limits"]["GridEnvelope"]['high'][0].split()]
        uc = [float(x) for x in cov_desc["boundedBy"]["Envelope"]["upperCorner"][0].split()]
        lc = [float(x) for x in cov_desc["boundedBy"]["Envelope"]["lowerCorner"][0].split()]
        ds_region = [lc[1], uc[1], lc[0], uc[0]]
        resx = (uc[1] - lc[1]) / hl[0]
        resy = (uc[0] - lc[0]) / hl[1]
        data = {'request': 'GetCoverage', 'version': '1.0.0', 'service': 'WCS',
                'resx': resx, 'resy': resy, 'crs': 'EPSG:4326', 'format': fmt,
                'coverage': coverage, 'Identifier': coverage}
        if region is not None: data['bbox'] = regions.region_format(region, 'bbox')        
        try:
            enc_data = urllib.urlencode(data)
        except: enc_data = urllib.parse.urlencode(data)
        return('{}{}'.format(self.url, enc_data))
    
    def fetch_coverage(coverage, region = None):
        
        c_url = self._get_coverage_url(coverage, region)
        return(fetch_file(c_url, '{}_{}.tif'.format(coverage, regions.region_format(region, 'fn')), params=data, verbose=True))

# class fetches_dataset(dlim.xyz_dataset):
#     """providing a fetches datalist
#     """
    
#     def __init__(self, **kwargs):
#         super().__init__(**kwargs)
            
# ## =============================================================================
# ##
# ## GMRT Fetch
# ##
# ## fetch extracts of the GMRT. - Global Extents
# ## https://www.gmrt.org/index.php
# ##
# ## =============================================================================
# class gmrt(fetches_dataset):
#     """represnting the remote raster data from the GMRT
#     """
    
#     def __init__(self, **kwargs):
#         super().__init__(**kwargs)
        
#         ## ==============================================
#         ## GMRT URLs and directories
#         ## ==============================================    
#         self._gmrt_grid_url = "https://www.gmrt.org:443/services/GridServer?"
#         self._gmrt_grid_urls_url = "https://www.gmrt.org:443/services/GridServer/urls?"
#         self._gmrt_grid_metadata_url = "https://www.gmrt.org/services/GridServer/metadata?"
#         self._outdir = os.path.join(os.getcwd(), 'gmrt')

#     def generate_inf(self):
#         """generate a infos dictionary from the gmrt dataset

#         Returns:
#           dict: a data-entry infos dictionary
#         """
        
#         self.infos['name'] = self.fn
#         self.infos['hash'] = None
#         this_region = regions.region().from_list([-179,179,-89,89])
#         self.infos['minmax'] = this_region.export_as_list()
#         self.infos['numpts'] = 0
#         self.infos['wkt'] = this_region.export_as_wkt()
#         return(self.infos)        

#     def parse(self, fmt = 'geotiff', res = 'max', layer = 'topo'):
#         """Parse the gmrt module for urls
#         """

#         if self.region is not None:
#             self.inf()
#             if layer != 'topo' and layer != 'topo-mask': layer = 'topo'
#             inf_region = regions.region().from_list(self.infos['minmax'])
#             if regions.regions_intersect_p(inf_region, self.region):
#                 _data = {'north':self.region.ymax, 'west':self.region.xmin,
#                          'south':self.region.ymin, 'east':self.region.xmax,
#                          'mformat':'json', 'resolution':res, 'format':fmt}

#                 _req = fetch_req(self._gmrt_grid_urls_url, params = _data, tries = 10, timeout = 2)
#                 if _req is not None and _req.status_code == 200:
#                     gmrt_urls = _req.json()
#                     for url in gmrt_urls:
#                         opts = {}
#                         url_base = url.split('?')[0]
#                         for url_opt in url.split('?')[1].split('&'):
#                             opt_kp = url_opt.split('=')
#                             opts[opt_kp[0]] = opt_kp[1]
#                         opts['layer'] = layer
#                         try:
#                             url_enc = urllib.urlencode(opts)
#                         except: url_enc = urllib.parse.urlencode(opts)
#                         this_url = '{}?{}'.format(url_base, url_enc)
#                         url_region = regions.region().from_list([float(opts['west']), float(opts['east']), float(opts['south']), float(opts['north'])])
#                         outf = 'gmrt_{}_{}.{}'.format(opts['layer'], url_region.format('fn'), utils.gdal_fext(opts['format']))                      
#                         self.fn = this_url
#                         self.data_entries.append(self)
#         return(self)
    
#     def yield_xyz(self):
#         for i, entry in enumerate(self.data_entries):
#             src_gmrt = 'gmrt_{}_{}.tif'.format(i, self.region.format('fn'))
#             utils.echo_msg('>{}<'.format(self.region.format('fn')))
#             gmrt_ds = dlim.dataset_factory(fn = src_gmrt, data_format = 200, weight = self.weight, src_region = self.region, epsg = self.epsg, warp = self.warp, verbose = self.verbose)
#             if fetch_file(entry.fn, src_gmrt, verbose = self.verbose) == 0:
#                 for xyz in gmrt_ds.acquire_raster_file().parse().yield_xyz():
#                     yield(xyz)
#             else: utils.echo_error_msg('failed to fetch remote file, {}...'.format(src_gmrt))
#             utils.remove_glob(src_gmrt)

# class fetches_factory(dlim.dataset_factory):
    
#     _fetch_modules = {
#         -11: {'name': 'gmrt',
#               'fmts': ['gmrt'],
#               'class': gmrt,
#               },
#     }

#     def __init__(self, fn = None, data_format = None, weight = 1, epsg = 4326, name = "<xyz-dataset>", parent = None, src_region = None, warp = None, verbose = False):
#         self.name = name
#         self.data_format = data_format
#         self.epsg = epsg
#         self.warp = warp
#         self.region = src_region
#         self.parent = parent
#         self.weight = weight
#         self.verbose = verbose
#         self.fn = fn
#         if self.data_format is None:
#             self.guess_data_format()

#         for key in self._fetch_modules.keys():
#             self.data_types[key] = self._fetch_modules[key]

        
#     def acquire_gmrt(self, **kwargs):
#         return(gmrt(fn = self.fn, data_format = self.data_format, weight = self.weight, src_region = self.region,
#                     epsg = self.epsg, warp = self.warp, name = self.name, parent = self.parent, verbose = self.verbose,
#                     **kwargs))
            
# ## ==============================================
# ## Run the Fetch Module(s)
# ## ==============================================
# def fetch(fg):
#     out_results = []
#     stop_threads = False
#     if fg is None:
#         utils.echo_error_msg('invalid configuration, {}'.format(fg))
#         sys.exit(-1)
        
#     for fetch_mod in fg['mods'].keys():
#         if stop_threads: break
#         status = 0
#         args = tuple(fg['mods'][fetch_mod])
        
#         if fg['region'] is None or fg['region'][0] is None:
#             this_region = None
#         else: this_region = regions.region_buffer(fg['region'], 5, pct = True)
#         fl = _fetch_modules[fetch_mod](this_region, fg['where'], lambda: stop_threads, fg['verbose'])
#         args_d = utils.args2dict(args)
#         if fg['verbose']: _prog = utils._progress('running FETCHES module < {} > with {} [{}]...'.format(fetch_mod, args, args_d))

#         ## ==============================================
#         ## Run update in a thread to cleanly close FRED
#         ## ============================================== 
#         if fg['update_p']:
#             t = threading.Thread(target = fl._update, args = ())
#             t.daemon = True
#             try:
#                 t.start()
#                 while True:
#                     time.sleep(2)
#                     if fg['verbose']: _prog.update(msg='updating FETCHES module < {} >...'.format(fetch_mod))
#                     if not t.is_alive():
#                         break
#             except (KeyboardInterrupt, SystemExit):
#                 utils.echo_error_msg('user breakage...please wait while fetches exits.')
#                 stop_threads = True
#                 status = -1
#             except Exception as e:
#                 utils.echo_error_msg(e)
#                 stop_threads = True
#                 status = -1
#             t.join()
#             if fg['verbose']: _prog.end(status, 'updated FETCHES module < {} >.'.format(fetch_mod))
            
#         if this_region is not None:
#             if fg['index_p']:
#                 _filter_FRED_index(fl)
#                 continue
#             ## ==============================================    
#             ## Fetch the data
#             ## fetching will be done in a queue with 3 threads fetching the data at a time.
#             ##
#             ## fetch_p must be true to fetch the data to _outdir
#             ## list_p will list the urls
#             ## dump_p will dump the xyz data from the fetch module to stdout
#             ## proc_p will output the xyz data to file in _outdir
#             ## ==============================================
#             fr = FetchResults(fl._parse_results(**args_d), this_region, fl._outdir, fl, fg, lambda: stop_threads)
#             fr.daemon = True
#             try:
#                 fr.start()
#                 while True:
#                     time.sleep(2)
#                     if fg['verbose']: _prog.update()
#                     if not fr.is_alive():
#                         break
#             except (KeyboardInterrupt, SystemExit):
#                 utils.echo_error_msg('user breakage...please wait for while fetches exits.')
#                 stop_threads = True
#                 status = -1
#                 while not fr.fetch_q.empty():
#                     try:
#                         fr.fetch_q.get(False)
#                     except Empty: continue
#                     fr.fetch_q.task_done()
#             except Exception as e:
#                 stop_threads = True
#                 status = -1
#                 utils.echo_error_msg(e)
#             fr.join()
#         if fg['verbose']: _prog.end(status, 'ran FETCHES module < {} > with {} [{}]...'.format(fetch_mod, args, fg))

# src_region = regions.region().from_list([-90,-89.75,28.25,28.5])

# #df = dlim.dataset_factory()
# #for key in _fetch_modules.keys():
# #    df.data_types[key] = _fetch_modules[key]
# #df.fn = 'gmrt'


# #d = dlim.dataset_factory(fn = 'gmrt', src_region = src_region, verbose = True)
# #d.add_data_type(fetches_factory()._fetch_modules, gmrt)
# #d.acquire_gmrt = 

# #print(d.data_types)
# #print(df.infos)
# f = fetches_factory(fn = 'gmrt', src_region = src_region, verbose = True)
# f.guess_data_format()
# print(f.data_format)
# g = f.acquire()
# print(g)
# #acquire_dataset()
# #g = f.acquire_gmrt()
# #g = gmrt(fn = 'gmrt', src_region = src_region)
# g.parse()
# #g.inf()
# g.echo()
# #for xyz in g.yield_xyz():
# #    xyz.dump()

# #fdl = fetches_dataset().from_string('gmrt')

# #fdl.dump_xyz(src_region = src_region)
        
### End
