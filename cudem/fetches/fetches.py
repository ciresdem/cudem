### fetches.py
##
## Copyright (c) 2010 - 2021 CIRES Coastal DEM Team
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
### Code:

import os
import sys
import re
import time
import requests
import urllib
import lxml.etree
import threading
try:
    import Queue as queue
except: import queue as queue
from osgeo import ogr
from cudem import utils
from cudem import regions
from cudem import datasets
from cudem import xyzfun

__version__ = '0.1.0'
this_dir, this_filename = os.path.split(__file__)
fetchdata = os.path.join(this_dir, 'data')

## =============================================================================
##
## Fetching Functions
##
## Generic fetching and processing functions, etc.
##
## The `fetch_results` class will fetch a list of fetch results [[url, file-name, data-type]...]
## in a queue `fetch_queue` using 3 threads; set `p` and `s` as a fetch module object to processes
## and dump XYZ data from the fetched results, respectively. Use `fetch_file` to fetch single files.
##
## =============================================================================
r_headers = { 'User-Agent': 'Fetches v%s' %(__version__) }

def urlencode(opts):
    try:
        url_enc = urllib.urlencode(opts)
    except:
        url_enc = urllib.parse.urlencode(opts)
    return(url_enc)

class Fetch:

    def __init__(self, url=None, callback=lambda: False, verbose=None):
        self.url = url
        self.callback = callback
        self.verbose = verbose

    def fetch_req(self, params=None, tries=5, timeout=2, read_timeout=10):
        """fetch src_url and return the requests object"""
        
        if tries <= 0:
            utils.echo_error_msg('max-tries exhausted')
            return(None)
        try:
            return(requests.get(self.url, stream=True, params=params, timeout=(timeout,read_timeout), headers=r_headers))
        except:
            return(self.fetch_req(params=params, tries=tries - 1, timeout=timeout + 1, read_timeout=read_timeout + 10))
            
    def fetch_as_html(self, timeout=2):
        """fetch src_url and return it as an HTML object"""
    
        req = self.fetch_req(timeout=timeout)
        if req:
            return(lh.document_fromstring(req.text))
        else:
            return(None)

    def fetch_xml(self, timeout=2, read_timeout=10):
        """fetch src_url and return it as an XML object"""
    
        results = lxml.etree.fromstring('<?xml version="1.0"?><!DOCTYPE _[<!ELEMENT _ EMPTY>]><_/>'.encode('utf-8'))
        try:
            req = fetch_req(src_url, timeout=timeout, read_timeout=read_timeout)
            results = lxml.etree.fromstring(req.text.encode('utf-8'))
        except:
            utils.echo_error_msg('could not access {}'.format(self.url))
        return(results)

    def fetch_file(self, dst_fn, params=None, datatype=None, overwrite=False, timeout=140, read_timeout=320):
        """fetch src_url and save to dst_fn"""
    
        status = 0
        req = None

        if self.verbose:
            progress = utils.CliProgress('fetching remote file: {}...'.format(os.path.basename(self.url)[:20]))
        if not os.path.exists(os.path.dirname(dst_fn)):
            try:
                os.makedirs(os.path.dirname(dst_fn))
            except: pass 
        if not os.path.exists(dst_fn) or overwrite:
            try:
                with requests.get(self.url, stream=True, params=params, headers=r_headers,
                                  timeout=(timeout,read_timeout)) as req:
                    req_h = req.headers
                    if req.status_code == 200:
                        curr_chunk = 0
                        with open(dst_fn, 'wb') as local_file:
                            for chunk in req.iter_content(chunk_size = 8196):
                                if self.callback():
                                    break
                                if self.verbose:
                                    progress.update()
                                if chunk:
                                    local_file.write(chunk)
                    else:
                        utils.echo_error_msg('server returned: {}'.format(req.status_code))
            except Exception as e:
                utils.echo_error_msg(e)
                status = -1
        if not os.path.exists(dst_fn) or os.stat(dst_fn).st_size ==  0:
            status = -1
        if self.verbose:
            progress.end(status, 'fetched remote file: {}.'.format(os.path.basename(dst_fn)[:20]))
        return(status)

    def fetch_ftp_file(self, dst_fn, params=None, datatype=None, overwrite=False):
        """fetch an ftp file via urllib"""

        status = 0
        f = None

        if self.verbose:
            utils.echo_msg('fetching remote ftp file: {}...'.format(self.url[:20]))
        if not os.path.exists(os.path.dirname(dst_fn)):
            try:
                os.makedirs(os.path.dirname(dst_fn))
            except:
                pass 
        try:
            f = urllib.request.urlopen(self.url)
        except:
            f = None
            status - 1

        if f is not None:
            with open(dst_fn, 'wb') as local_file:
                 local_file.write(f.read())
            if self.verbose:
                utils.echo_msg('fetched remote ftp file: {}.'.format(os.path.basename(self.url)))
        return(status)

def fetch_queue(q, m, p=False):
    """fetch queue `q` of fetch results\
    each fetch queue should be a list of the following:
    [remote_data_url, local_data_path, regions.region, lambda: stop_p, data-type]
    if region is defined, will prompt the queue to process the data to the given region.
    """
    
    while True:
        fetch_args = q.get()
        #this_region = fetch_args[2]
        if not m.callback():
            if not os.path.exists(os.path.dirname(fetch_args[1])):
                try:
                    os.makedirs(os.path.dirname(fetch_args[1]))
                except: pass
            
            #if this_region is None:
            if not p:
                if fetch_args[0].split(':')[0] == 'ftp':
                    Fetch(url=fetch_args[0], callback=m.callback, verbose=m.verbose).fetch_ftp_file(fetch_args[1])
                else:
                    Fetch(url=fetch_args[0], callback=m.callback, verbose=m.verbose).fetch_file(fetch_args[1])
            else:
                o_x_fn = fetch_args[1] + '.xyz'
                utils.echo_msg('processing local file: {}'.format(o_x_fn))
                if not os.path.exists(o_x_fn):
                    with open(o_x_fn, 'w') as out_xyz:
                        m.dump_xyz(fetch_args, dst_port=out_xyz)
                    if os.stat(o_x_fn).st_size == 0:
                        utils.remove_glob(o_x_fn)
        q.task_done()

class fetch_results(threading.Thread):
    """fetch results gathered from a fetch module.
    results is a list of URLs with data type
    e.g. results = [[http://data/url.xyz.gz, /home/user/data/url.xyz.gz, data-type], ...]
    """
    
    #def __init__(self, results, out_dir, region=None, fetch_module=None, callback=lambda: False):
    def __init__(self, mod, want_proc=False):
        threading.Thread.__init__(self)
        self.fetch_q = queue.Queue()
        self.mod = mod
        self.outdir_ = self.mod._outdir
        self.want_proc = want_proc
        if len(self.mod.results) == 0:
            self.mod.run()
        
    def run(self):
        for _ in range(3):
            t = threading.Thread(target=fetch_queue, args=(self.fetch_q, self.mod, self.want_proc))
            t.daemon = True
            t.start()
        for row in self.mod.results:
            #self.fetch_q.put([row[0], os.path.join(self.outdir_, row[1]), self.stop_threads, row[2], self.fetch_module])
            self.fetch_q.put([row[0], os.path.join(self.outdir_, row[1]), row[2], self.mod])
        self.fetch_q.join()

class FetchModule:

    def __init__(self, src_region=None, callback=lambda: False, weight=None, verbose=True, warp=None):
        self.region = src_region
        self.callback = callback
        self.weight = weight
        self.verbose = verbose
        self.status = 0
        self.results = []
        self.warp = warp
        self.name = None

    def run(self):
        raise(NotImplementedError)

    def yield_xyz(self, **kwargs):
        raise(NotImplementedError)

    def dump_xyz(self, entry, dst_port=sys.stdout, **kwargs):
        for xyz in self.yield_xyz(entry, **kwargs):
            xyz.dump(include_w=True if self.weight is not None else False,
                     dst_port=dst_port, encode=False)
            
    def yield_results_to_xyz(self, **kwargs):
        if len(self.results == 0):
            self.run()
        for entry in self.results:
            for xyz in self.yield_xyz(**kwargs):
                yield(xyz)
                
    def dump_results_to_xyz(self, dst_port=sys.stdout, **kwargs):
        for xyz in self.yield_results_to_xyz(**kwargs):
            xyz.dump(include_w=True if self.weight is not None else False,
                     dst_port=dst_port, encode=False)

## =============================================================================
##
## MB Fetch
##
## Fetch Multibeam bathymetric surveys from NOAA
## MBSystem is required to process the resulting data
##
## =============================================================================
class MB(FetchModule):

    def __init__(self, processed=False, inc=None, **kwargs):
        self._mb_data_url = "https://data.ngdc.noaa.gov/platforms/"
        self._mb_metadata_url = "https://data.noaa.gov/waf/NOAA/NESDIS/NGDC/MGG/Multibeam/iso/"
        self._mb_search_url = "https://maps.ngdc.noaa.gov/mapviewer-support/multibeam/files.groovy?"
        self._mb_autogrid = "https://www.ngdc.noaa.gov/maps/autogrid/"
        self._mb_html = "https://www.ngdc.noaa.gov/"
        self._outdir = os.path.join(os.getcwd(), 'mb')
        self._urls = [self._mb_data_url, self._mb_metadata_url, self._mb_autogrid]
        self.name = 'mb'
        self.processed_p = processed
        self.inc = utils.str2inc(inc)
        super().__init__(**kwargs)

    def mb_inf_data_format(self, src_inf):
        """extract the data format from the mbsystem inf file.

        Args:
          src_inf (str): the source mbsystem .inf file pathname

        Returns:
          str: the mbsystem datalist format number
        """

        with open(src_inf) as iob:
            for il in iob:
                til = il.split()
                if len(til) > 1:
                    if til[0] == 'MBIO':
                        return(til[4])
        
    def run(self):
        these_surveys = {}
        these_versions = {}
        if self.region is None: return([])
        _req = Fetch(self._mb_search_url).fetch_req(params={'geometry': self.region.format('bbox')}, timeout=20)
        if _req is not None and _req.status_code == 200:
            survey_list = _req.text.split('\n')[:-1]
            for r in survey_list:
                dst_pfn = r.split(' ')[0]
                dst_fn = dst_pfn.split('/')[-1:][0]
                survey = dst_pfn.split('/')[6]
                dn = r.split(' ')[0].split('/')[:-1]
                version = dst_pfn.split('/')[9][-1]
                data_url = self._mb_data_url + '/'.join(r.split('/')[3:])
                if survey in these_surveys.keys():
                    if version in these_surveys[survey].keys():
                        these_surveys[survey][version].append([data_url.split(' ')[0], '/'.join([survey, dst_fn]), 'mb'])
                    else: these_surveys[survey][version] = [[data_url.split(' ')[0], '/'.join([survey, dst_fn]), 'mb']]
                else: these_surveys[survey] = {version: [[data_url.split(' ')[0], '/'.join([survey, dst_fn]), 'mb']]}
        else: utils.echo_error_msg('{}'.format(_req.reason))
                    
        for key in these_surveys.keys():
            if self.processed_p:
                if '2' in these_surveys[key].keys():
                    for v2 in these_surveys[key]['2']:
                        self.results.append(v2)
                else:
                    for v1 in these_surveys[key]['1']:
                        self.results.append(v1)
            else:
                for keys in these_surveys[key].keys():
                    for survs in these_surveys[key][keys]:
                        self.results.append(survs)

    def yield_xyz(self, entry):
        src_data = os.path.basename(entry[1])
        if Fetch(entry[0], callback=self.callback, verbose=self.verbose).fetch_file(src_data) == 0:
            src_xyz = os.path.basename(src_data) + '.xyz'
            out, status = utils.run_cmd('mblist -MA -OXYZ -I{}  > {}'.format(src_data, src_xyz), verbose=False)
            if status != 0:
                if Fetch('{}.inf'.format(entry[0]), callback=self.callback, verbose=self.verbose).fetch_file('{}.inf'.format(src_data)) == 0:
                    mb_fmt = self.mb_inf_data_format('{}.inf'.format(src_mb))
                    remove_glob('{}.inf'.format(entry[0]))
                    out, status = run_cmd('mblist -F{} -MA -OXYZ -I{}  > {}'.format(mb_fmt, src_data, src_xyz), verbose=False)
            if status == 0:
                _ds = datasets.XYZFile(fn=src_xyz, delim='\t', data_format=168, epsg=4326, warp=self.warp,
                                       name=os.path.basename(entry[1]), src_region=self.region, verbose=self.verbose, remote=True)

                if self.inc is not None:
                    xyz_func = lambda p: _ds.dump_xyz(dst_port=p, encode=True)
                    for xyz in utils.yield_cmd('gmt blockmedian -I{:.10f} {} -r -V'.format(self.inc, self.region.format('gmt')), verbose=self.verbose, data_fun=xyz_func):
                        yield(xyzfun.XYZPoint().from_list([float(x) for x in xyz.split()]))
                    
                for xyz in _ds.yield_xyz():
                    yield(xyz)

                utils.remove_glob(src_data, src_xyz)
            else:
                echo_error_msg('failed to process local file, {} [{}]...'.format(src_data, entry[0]))
                with open('{}'.format(os.path.join(self._outdir, 'fetch_{}_{}.err'.format(self._name, region_format(self.region, 'fn')))), 'a') as mb_err:
                    mb_err.write('{}\n'.format(','.join([src_mb, entry[0]])))
                os.rename(src_data, os.path.join(self._outdir, src_data))
                utils.remove_glob(src_xyz)
        else:
            utils.echo_error_msg('failed to fetch remote file, {}...'.format(src_data))

## =============================================================================
##
## USACE Fetch
##
## Fetch USACE bathymetric surveys via eHydro
##
## =============================================================================
class USACE(FetchModule):
    '''Fetch USACE bathymetric surveys'''
    
    def __init__(self, s_type=None, inc=None, **kwargs):
        self._usace_gj_api_url = 'https://opendata.arcgis.com/datasets/80a394bae6b547f1b5788074261e11f1_0.geojson'
        self._usace_gs_api_url = 'https://services7.arcgis.com/n1YM8pTrFmm7L4hs/arcgis/rest/services/eHydro_Survey_Data/FeatureServer/0/query?outFields=*&where=1%3D1'
        self._outdir = os.path.join(os.getcwd(), 'usace')
        self.name = 'usace'
        self.s_type = s_type
        self.inc = utils.str2inc(inc)
        super().__init__(**kwargs)

    def run(self):
        '''Run the USACE fetching module'''
        
        if self.region is None: return([])
        _data = {
            'geometry': self.region.format('bbox'),
            'inSR':4326,
            'outSR':4326,
            'f':'pjson',
        }
        _req = Fetch(self._usace_gs_api_url).fetch_req(params=_data)
        if _req is not None:
            survey_list = _req.json()
            for feature in survey_list['features']:
                fetch_fn = feature['attributes']['SOURCEDATALOCATION']
                if self.s_type is not None:
                    if feature['attributes']['SURVEYTYPE'].lower() == self.s_type.lower():
                        self.results.append([fetch_fn, fetch_fn.split('/')[-1], 'usace'])
                else: self.results.append([fetch_fn, fetch_fn.split('/')[-1], 'usace'])
        return(self)

    def yield_xyz(self, entry):
        src_zip = os.path.basename(entry[1])
        src_epsg = None
        src_region = None
        #xyzc['warp'] = epsg
        
        if Fetch(entry[0], callback=self.callback, verbose=self.verbose).fetch_file(src_zip) == 0:
            src_xmls = utils.p_unzip(src_zip, ['xml', 'XML'])
            for src_xml in src_xmls:
                if src_region is None:
                    this_xml = lxml.etree.parse(src_xml)
                    if this_xml is not None:
                        try:
                            w = this_xml.find('.//westbc').text
                            e = this_xml.find('.//eastbc').text
                            n = this_xml.find('.//northbc').text
                            s = this_xml.find('.//southbc').text
                            src_region = regions.Region().from_list([float(w), float(e), float(s), float(n)])
                        except:
                            utils.echo_warning_msg('could not determine survey bb from {}'.format(src_xml))
                if src_epsg is None:
                    try:
                        prj = this_xml.find('.//gridsysn').text
                        szone = this_xml.find('.//spcszone').text
                        utils.echo_msg('zone: {}'.format(szone))                        
                        src_epsg = int(utils.FIPS_TO_EPSG[szone])
                    except:
                        utils.echo_warning_msg('could not determine state plane zone from {}'.format(src_xml))
                utils.remove_glob(src_xml)
                
            if src_epsg is None:
                this_geom = src_region.export_as_geom()
                sp_fn = os.path.join(fetchdata, 'stateplane.geojson')
                try:
                    sp = ogr.Open(sp_fn)
                    layer = sp.GetLayer()
                
                    for feature in layer:
                        geom = feature.GetGeometryRef()
                        if this_geom.Intersects(geom):
                            src_epsg = feature.GetField('EPSG')
                    sp = None
                except: pass

            src_usaces = utils.p_unzip(src_zip, ['XYZ', 'xyz', 'dat'])
            for src_usace in src_usaces:
                _dl = datasets.XYZFile(fn=src_usace, epsg=src_epsg, warp=self.warp, src_region=src_region, name=src_usace, verbose=self.verbose, remote=True)
                for xyz in _dl.yield_xyz():
                    yield(xyz)
                utils.remove_glob(src_usace)
        else: utils.echo_error_msg('failed to fetch remote file, {}...'.format(entry[0]))
        utils.remove_glob(src_zip)
            
## =============================================================================
##
## GMRT Fetch
##
## fetch extracts of the GMRT. - Global Extents
## https://www.gmrt.org/index.php
##
## =============================================================================
class GMRT(FetchModule):
    '''Fetch raster data from the GMRT'''
    
    def __init__(self, res='max', fmt='geotiff', **kwargs):
        super().__init__(**kwargs) 
        #if self.verbose: utils.echo_msg('loading GMRT fetch module...')
        self._gmrt_grid_url = "https://www.gmrt.org:443/services/GridServer?"
        self._gmrt_grid_urls_url = "https://www.gmrt.org:443/services/GridServer/urls?"
        self._gmrt_grid_metadata_url = "https://www.gmrt.org/services/GridServer/metadata?"        
        self._outdir = os.path.join(os.getcwd(), 'gmrt')
        self.name = 'gmrt'
        self.res = res
        self.fmt = fmt
        
    def run(self):
        '''Run the GMRT fetching module'''

        if self.region is None: return([])
        self.data = {
            'north':self.region.ymax,
            'west':self.region.xmin,
            'south':self.region.ymin,
            'east':self.region.xmax,
            'mformat':'json',
            'resolution':self.res,
            'format':self.fmt,
        }
        
        req = Fetch(self._gmrt_grid_urls_url).fetch_req(params=self.data, tries=10, timeout=2)
        if req is not None:
            gmrt_urls = req.json()
            for url in gmrt_urls:
                opts = {}
                for url_opt in url.split('?')[1].split('&'):
                    opt_kp = url_opt.split('=')
                    opts[opt_kp[0]] = opt_kp[1]
                url_region = regions.Region().from_list([float(opts['west']), float(opts['east']), float(opts['south']), float(opts['north'])])
                outf = 'gmrt_{}_{}.tif'.format(opts['layer'], url_region.format('fn'))
                self.results.append([url, outf, 'gmrt'])
        return(self)

    def yield_xyz(self, entry):
        src_data = 'gmrt_tmp.tif'
        if Fetch(entry[0], callback=self.callback, verbose=self.verbose).fetch_file(src_data) == 0:
            gmrt_ds = datasets.RasterFile(fn=src_data, data_format=200, epsg=4326, warp=self.warp,
                                          name=src_data, src_region=self.region, verbose=self.verbose)
            for xyz in gmrt_ds.yield_xyz():
                yield(xyz)
        else: utils.echo_error_msg('failed to fetch remote file, {}...'.format(src_data))
        utils.remove_glob('{}*'.format(src_data))

## =============================================================================
##
## mar_grav - Sattelite Altimetry Topography from Scripps
##
## https://topex.ucsd.edu/WWW_html/mar_grav.html
## ftp://topex.ucsd.edu/pub/global_grav_1min/
## https://topex.ucsd.edu/marine_grav/explore_grav.html
## https://topex.ucsd.edu/marine_grav/white_paper.pdf
##
## =============================================================================
class MarGrav(FetchModule):
    '''Fetch mar_grav sattelite altimetry topography'''
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs) 
        self._mar_grav_url = 'https://topex.ucsd.edu/cgi-bin/get_data.cgi'
        self._outdir = os.path.join(os.getcwd(), 'mar_grav')
        self.name = 'mar_grav'

    def run(self):
        '''Run the mar_grav fetching module.'''
        
        if self.region is None: return([])
        self.data = {
            'north':self.region.ymax,
            'west':self.region.xmin,
            'south':self.region.ymin,
            'east':self.region.xmax,
            'mag':1,
        }
        _req = Fetch(self._mar_grav_url).fetch_req(params=self.data, tries=10, timeout=2)
        if _req is not None:
            url = _req.url
            outf = 'mar_grav_{}.xyz'.format(self.region.format('fn'))
            self.results.append([url, outf, 'mar_grav'])
        return(self)
        
    def yield_xyz(self, entry):
        src_data = 'mar_grav_tmp.xyz'
        if Fetch(entry[0], callback=self.callback, verbose=self.verbose).fetch_file(src_data) == 0:
            _ds = datasets.XYZFile(fn=src_data, data_format=168, skip=1, x_offset=-360, epsg=4326, warp=self.warp,
                                   name=src_data, src_region=self.region, verbose=self.verbose, remote=True)
            for xyz in _ds.yield_xyz():
                yield(xyz)
        else: utils.echo_error_msg('failed to fetch remote file, {}...'.format(src_data))
        utils.remove_glob('{}*'.format(src_data))                
                
## =============================================================================
##
## SRTM Plus
##
## Fetch srtm+ data
## https://topex.ucsd.edu/WWW_html/srtm15_plus.html
## http://topex.ucsd.edu/sandwell/publications/180_Tozer_SRTM15+.pdf
##
## =============================================================================
class SRTMPlus(FetchModule):
    '''Fetch SRTM+ data'''
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs) 
        self._srtm_url = 'https://topex.ucsd.edu/cgi-bin/get_srtm15.cgi'
        self._outdir = os.path.join(os.getcwd(), 'srtm_plus')
        self.name = 'srtm_plus'

    def run(self):
        '''Run the SRTM fetching module.'''
        
        if self.region is None: return([])

        self.data = {
            'north':self.region.ymax,
            'west':self.region.xmin,
            'south':self.region.ymin,
            'east':self.region.xmax,
        }

        _req = Fetch(self._srtm_url).fetch_req(params=self.data, tries=10, timeout=2)
        if _req is not None:
            url = _req.url
            outf = 'srtm_{}.xyz'.format(self.region.format('fn'))
            self.results.append([url, outf, 'srtm'])
        return(self)

    def yield_xyz(self, entry):
        src_data = 'srtm_plus_tmp.xyz'
        if Fetch(entry[0], callback=self.callback, verbose=self.verbose).fetch_file(src_data) == 0:
            _ds = datasets.XYZFile(fn=src_data, data_format=168, skip=1, epsg=4326, warp=self.warp,
                                   name=src_data, src_region=self.region, verbose=self.verbose, remote=True)

            for xyz in _ds.yield_xyz():
                yield(xyz)
        else: utils.echo_error_msg('failed to fetch remote file, {}...'.format(src_data))
        utils.remove_glob('{}*'.format(src_data))
        
## ==============================================
## Fetches Module Parser
## ==============================================
class FetchesFactory:

    mods = {
        'gmrt': {'description': """The Global Multi-Resolution Topography (GMRT) synthesis is a multi-resolutional 
compilation of edited multibeam sonar data collected by scientists and institutions worldwide, that is 
reviewed, processed and gridded by the GMRT Team and merged into a single continuously updated compilation 
of global elevation data. The synthesis began in 1992 as the Ridge Multibeam Synthesis (RMBS), was expanded 
to include multibeam bathymetry data from the Southern Ocean, and now includes bathymetry from throughout 
the global and coastal oceans."""},
        'srtm_plus': {'description': """"""},
        'mar_grav': {'description': """"""},
        'mb': {'description': """"""},
        'usace': {'description': """"""},
    }
    
    def __init__(self, mod=None, src_region=None, callback=lambda: False, weight=None, verbose=True):
        self.mod=mod
        self.mod_args=()
        self.region = src_region
        self.callback = callback
        self.weight = weight
        self.verbose = verbose
        self.status = 0
        self.results = []

        if self.mod is not None:
            self.parse_mod()
        
    def parse_mod(self):
        opts = self.mod.split(':')
        if opts[0] in FetchesFactory.mods.keys():
            self.mod_args = utils.args2dict(list(opts[1:]))
            self.mod = opts[0]
        else:
            utils.echo_error_msg('invalid module name `{}`'.format(opts[0]))
        return(self)

    def add_module(self, type_def={}):
        for key in type_def.keys():
            self.mods[key] = type_def[key]

    def acquire_mb(self, **kwargs):
        return(MB(
            src_region=self.region, callback=self.callback, weight=self.weight, verbose=self.verbose, **kwargs, **self.mod_args))

    def acquire_usace(self, **kwargs):
        return(USACE(
            src_region=self.region, callback=self.callback, weight=self.weight, verbose=self.verbose, **kwargs, **self.mod_args))

    def acquire_gmrt(self, **kwargs):
        return(GMRT(
            src_region=self.region, callback=self.callback, weight=self.weight, verbose=self.verbose, **kwargs, **self.mod_args))

    def acquire_srtm_plus(self, **kwargs):
        return(SRTMPlus(
            src_region=self.region, callback=self.callback, weight=self.weight, verbose=self.verbose, **kwargs, **self.mod_args))

    def acquire_mar_grav(self, **kwargs):
        return(MarGrav(
            src_region=self.region, callback=self.callback, weight=self.weight, verbose=self.verbose, **kwargs, **self.mod_args))

    def acquire(self, **kwargs):
        if self.mod == 'mb':
            return(self.acquire_mb(**kwargs))

        if self.mod == 'usace':
            return(self.acquire_usace(**kwargs))
        
        if self.mod == 'gmrt':
            return(self.acquire_gmrt(**kwargs))

        if self.mod == 'srtm_plus':
            return(self.acquire_srtm_plus(**kwargs))
        
        if self.mod == 'mar_grav':
            return(self.acquire_mar_grav(**kwargs))

_fetches_module_short_desc = lambda: ', '.join(
    ['{}'.format(key) for key in FetchesFactory().mods])
_fetches_module_long_desc = lambda x: 'fetches modules:\n% fetches ... <mod>:key=val:key=val...\n\n  ' + '\n  '.join(
    ['\033[1m{:14}\033[0m{}\n'.format(key, x[key]['description']) for key in x]) + '\n'

## ==============================================
## Command-line Interface (CLI)
## $ fetches
##
## fetches cli
## ==============================================
fetches_usage = """{cmd} ({f_version}): Fetches; Fetch and process remote elevation data

usage: {cmd} [ -hiqwPRW [ args ] ] MODULE ...

Options:
  -R, --region\t\tRestrict processing to the desired REGION 
\t\t\tWhere a REGION is [ xmin/xmax/ymin/ymax/[ zmin/zmax/[ wmin/wmax ] ] ]
\t\t\tUse '-' to indicate no bounding range; e.g. -R -/-/-/-/-10/10/1/-
\t\t\tOR an OGR-compatible vector file with regional polygons. 
\t\t\tIf a vector file is supplied, will use each region found therein.
\t\t\tAppend :zmin/zmax/[ wmin/wmax ] to the file path to extended REGION.
  -l, --list\t\tReturn a list of fetch URLs in the given region.
  -p, --process\t\tProcess fetched elevation data to ASCII XYZ format in WGS84. <beta>

  --weights\t\tOutput WEIGHT values along with xyz
  --quiet\t\tLower the verbosity to a quiet

  --help\t\tPrint the usage text
  --version\t\tPrint the version information

Supported FETCHES modules: 
  {f_formats}

CIRES DEM home page: <http://ciresgroups.colorado.edu/coastalDEM>
""".format(cmd=os.path.basename(sys.argv[0]), 
           f_version=__version__,
           f_formats=_fetches_module_short_desc())

def fetches_cli(argv = sys.argv):
    """run fetches from command-line

    See `fetches_cli_usage` for full cli options.
    """

    i_regions = []
    these_regions = []
    mods = []
    mod_opts = {}
    want_weights = False
    want_list = False
    want_proc = False
    want_verbose = True
    stop_threads = False
    
    ## ==============================================
    ## parse command line arguments.
    ## ==============================================
    i = 1
    while i < len(argv):
        arg = argv[i]
        if arg == '--region' or arg == '-R':
            i_regions.append(str(argv[i + 1]))
            i = i + 1
        elif arg[:2] == '-R':
            i_regions.append(str(arg[2:]))
        elif arg == '--weights' or arg == '-w':
            want_weights = True
        elif arg == '--list' or arg == '-l':
            want_list = True
        elif arg == '--process' or arg == '-p':
            want_proc = True
        elif arg == '--quiet' or arg == '-q':
            want_verbose = False
        elif arg == '--help' or arg == '-h':
            sys.stderr.write(fetches_usage)
            sys.exit(1)
        elif arg == '--version' or arg == '-v':
            print('{}, version {}'.format(os.path.basename(sys.argv[0]), __version__))
            sys.exit(1)
        elif arg == '--modules' or arg == '-m':
            try:
                if argv[i + 1] in FetchesFactory.mods.keys():
                    sys.stderr.write(_fetches_module_long_desc({k: FetchesFactory.mods[k] for k in (argv[i + 1],)}))
                else: sys.stderr.write(_fetches_module_long_desc(FetchesFactory.mods))
            except: sys.stderr.write(_fetches_module_long_desc(FetchesFactory.mods))
            sys.exit(0)
        elif arg[0] == '-':
            sys.stderr.write(fetches_usage)
            sys.exit(0)
        else:
            mods.append(arg)
            
        i = i + 1

    if len(mods) == 0:
        sys.stderr.write(fetches_usage)
        utils.echo_error_msg('you must select at least one fetch module')
        sys.exit(-1)
        
    for i_region in i_regions:
        tmp_region = regions.Region().from_string(i_region)
        if tmp_region.valid_p():
            these_regions.append(tmp_region)
        else:
            tmp_region = regions.ogr_wkts(i_region)
            for i in tmp_region:
                if i.valid_p():
                    these_regions.append(i)
                    
    if len(these_regions) == 0:
        these_regions = [None]
    else:
        if want_verbose:
            utils.echo_msg('parsed {} region(s)'.format(len(these_regions)))

    for rn, this_region in enumerate(these_regions):
        if stop_threads: return
        x_fs = [FetchesFactory(mod=mod, src_region=this_region, verbose=want_verbose).acquire(warp=4326) for mod in mods]
        for x_f in x_fs:
            utils.echo_msg('running fetch module {} on region {}...\
            '.format(x_f.name, this_region.format('str')))
            
            #r = x_f.run().results
            x_f.run()
            utils.echo_msg('found {} data files.'.format(len(x_f.results)))
            if want_list:
                for result in x_f.results:
                    print(result[0])
            else:
                fr = fetch_results(x_f, want_proc=want_proc)
                fr.daemon = True
                _p = utils.CliProgress('fetching {} remote data files'.format(len(x_f.results)))
                try:
                    fr.start()
                    while True:
                        time.sleep(2)
                        sys.stderr.write('\x1b[2K\r')
                        perc = float((len(x_f.results) - fr.fetch_q.qsize())) / len(x_f.results) * 100 if len(x_f.results) > 0 else 1
                        #if want_verbose: sys.stderr.write('fetches: fetching remote data files [{}%]'.format(perc))
                        if want_verbose:
                            _p.update_perc((len(x_f.results) - fr.fetch_q.qsize(), len(x_f.results)))
                        sys.stderr.flush()
                        if not fr.is_alive():
                            break
                except (KeyboardInterrupt, SystemExit):
                    utils.echo_error_msg('user breakage...please wait for while fetches exits.')
                    x_f.status = -1
                    stop_threads = True
                    while not fr.fetch_q.empty():
                        try:
                            fr.fetch_q.get(False)
                        except Empty: continue
                        fr.fetch_q.task_done()
                fr.join()
                _p.end(x_f.status, 'fetched {} remote data files'.format(len(x_f.results)))
            if want_verbose:
                utils.echo_msg('ran fetch module {} on region {}...\
            '.format(x_f.name, this_region.format('str')))
### End
