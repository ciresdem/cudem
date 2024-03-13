### fetches.py
##
## Copyright (c) 2010 - 2024 Regents of the University of Colorado
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
## Fetch elevation and related data from a variety of sources.
## Use CLI command 'fetches'
## or use FetchesFactory() to acquire and use a fetch module.
##
### Code:

import os
import sys
import time
import json
import requests
import urllib
import lxml.etree
import lxml.html as lh
from tqdm import tqdm
import mercantile
import csv

import threading
try:
    import Queue as queue
except: import queue as queue

import boto3 # boto3 for aws api
from osgeo import ogr

import cudem
from cudem import utils
from cudem import regions
from cudem import factory
from cudem import FRED
from cudem import vdatums
from cudem import gdalfun

# for get_credentials
import base64
import netrc
from getpass import getpass
try:
    from urllib.parse import urlparse
    from urllib.request import urlopen, Request, build_opener, HTTPCookieProcessor
    from urllib.error import HTTPError, URLError
except ImportError:
    from urlparse import urlparse
    from urllib2 import urlopen, Request, HTTPError, URLError, build_opener, HTTPCookieProcessor

## Some servers don't like custom user agents...
#r_headers = { 'User-Agent': 'Fetches v%s' %(fetches.__version__) }
r_headers = { 'User-Agent': 'Mozilla/5.0 (X11; Linux x86_64; rv:45.0) Gecko/20100101 Firefox/45.0' }
namespaces = {
    'gmd': 'http://www.isotc211.org/2005/gmd', 
    'gmi': 'http://www.isotc211.org/2005/gmi', 
    'gco': 'http://www.isotc211.org/2005/gco',
    'gml': 'http://www.isotc211.org/2005/gml',
}

thredds_namespaces = {
    'th': 'http://www.unidata.ucar.edu/namespaces/thredds/InvCatalog/v1.0',
}

def urlencode(opts):
    try:
        url_enc = urllib.urlencode(opts)
    except:
        url_enc = urllib.parse.urlencode(opts)
    return(url_enc)

def xml2py(node):
    """parse an xml file into a python dictionary"""
    
    texts = {}
    if node is None:
        return(None)

    for child in list(node):
        child_key = lxml.etree.QName(child).localname
        if 'name' in child.attrib.keys():
            child_key = child.attrib['name']
        
        if '{http://www.w3.org/1999/xlink}href' in child.attrib.keys():
            href = child.attrib['{http://www.w3.org/1999/xlink}href']
        else: href = None
        
        if child.text is None or child.text.strip() == '':
            if href is not None:
                if child_key in texts.keys():
                    texts[child_key].append(href)
                else:
                    texts[child_key] = [href]
                    
            else:
                if child_key in texts.keys():
                    ck = xml2py(child)
                    texts[child_key][list(ck.keys())[0]].update(ck[list(ck.keys())[0]])
                else:
                    texts[child_key] = xml2py(child)
                    
        else:
            if child_key in texts.keys():
                texts[child_key].append(child.text)
            else:
                texts[child_key] = [child.text]
                
    return(texts)

def get_credentials(url, authenticator_url = 'https://urs.earthdata.nasa.gov'):
    """Get user credentials from .netrc or prompt for input."""
    
    credentials = None
    errprefix = ''
    try:
        info = netrc.netrc()
        username, account, password = info.authenticators(urlparse(authenticator_url).hostname)
        errprefix = 'netrc error: '
    except Exception as e:
        if (not ('No such file' in str(e))):
            print('netrc error: {0}'.format(str(e)))
        username = None
        password = None

    while not credentials:
        if not username:
            username = utils.get_username()
            password = utils.get_password()
            
        credentials = '{0}:{1}'.format(username, password)
        credentials = base64.b64encode(credentials.encode('ascii')).decode('ascii')

        if url:
            try:
                req = Request(url)
                req.add_header('Authorization', 'Basic {0}'.format(credentials))
                opener = build_opener(HTTPCookieProcessor())
                opener.open(req)
            except HTTPError:
                print(errprefix + 'Incorrect username or password')
                errprefix = ''
                credentials = None
                username = None
                password = None

    return(credentials)

## ==============================================
## a few convenience functions to parse common iso xml files
## ==============================================
class iso_xml:
    def __init__(self, xml_url, timeout=2, read_timeout=10):
        self.url = xml_url
        self.xml_doc = self._fetch(timeout=timeout, read_timeout=read_timeout)
        self.namespaces = {
            'gmd': 'http://www.isotc211.org/2005/gmd', 
            'gmi': 'http://www.isotc211.org/2005/gmi', 
            'gco': 'http://www.isotc211.org/2005/gco',
            'gml': 'http://www.isotc211.org/2005/gml',
            'th': 'http://www.unidata.ucar.edu/namespaces/thredds/InvCatalog/v1.0',
            'wms': 'http://www.opengis.net/wms',
        }
        
    def _fetch(self, timeout = 2, read_timeout = 10):
        return(Fetch(self.url).fetch_xml(timeout=timeout, read_timeout=read_timeout))

    def title(self):
        t = self.xml_doc.find('.//gmd:MD_DataIdentification/gmd:citation/gmd:CI_Citation/gmd:title/gco:CharacterString', namespaces = self.namespaces)
        return(t.text if t is not None else 'Unknown')
        
    def bounds(self, geom = True):
        wl = self.xml_doc.find('.//gmd:westBoundLongitude/gco:Decimal', namespaces = self.namespaces)
        el = self.xml_doc.find('.//gmd:eastBoundLongitude/gco:Decimal', namespaces = self.namespaces)
        sl = self.xml_doc.find('.//gmd:southBoundLatitude/gco:Decimal', namespaces = self.namespaces)                            
        nl = self.xml_doc.find('.//gmd:northBoundLatitude/gco:Decimal', namespaces = self.namespaces)                           
        if wl is not None and el is not None and sl is not None and nl is not None:
            region = [float(wl.text), float(el.text), float(sl.text), float(nl.text)]
            if geom: return(regions.Region().from_list([float(wl.text), float(el.text), float(sl.text), float(nl.text)]).export_as_geom())
            else: return(region)
        else: return(None)

    def polygon(self, geom = True):
        opoly = []
        polygon = self.xml_doc.find('.//{*}Polygon', namespaces = self.namespaces)
        if polygon is not None:
            nodes = polygon.findall('.//{*}pos', namespaces = self.namespaces)
            [opoly.append([float(x) for x in node.text.split()]) for node in nodes]
            if opoly[0][0] != opoly[-1][0] or opoly[0][1] != opoly[-1][1]:
                opoly.append(opoly[0])
            if geom:
                return(gdalfun.wkt2geom(regions.create_wkt_polygon(opoly)))
            else:
                return(opoly)
        else:
            return(None)
        
    def date(self):
        dt = self.xml_doc.find('.//gmd:date/gco:Date', namespaces = self.namespaces)
        if dt is None:
            dt = self.xml_doc.find('.//gmd:MD_DataIdentification/gmd:citation/gmd:CI_Citation/gmd:date/gmd:CI_Date/gmd:date/gco:Date', namespaces = self.namespaces)
            
        return(dt.text[:4] if dt is not None else '0000')

    def xml_date(self):
        mddate = self.xml_doc.find('.//gmd:dateStamp/gco:DateTime', namespaces = self.namespaces)
        
        return(utils.this_date() if mddate is None else mddate.text)
        
    def reference_system(self):
        ref_s = self.xml_doc.findall('.//gmd:MD_ReferenceSystem', namespaces = self.namespaces)
        if ref_s is None or len(ref_s) == 0:
            return(None, None)
        
        h_epsg = ref_s[0].find('.//gmd:code/gco:CharacterString', namespaces = self.namespaces)
        if h_epsg is not None:
            h_epsg = h_epsg.text.split(':')[-1]
        
        if len(ref_s) > 1:
            v_epsg = ref_s[1].find('.//gmd:code/gco:CharacterString', namespaces = self.namespaces)
            if v_epsg is not None:
                v_epsg = v_epsg.text.split(':')[-1]
                
        else:
            v_epsg = None
            
        return(h_epsg, v_epsg)

    def abstract(self):
        try:
            abstract = self.xml_doc.find('.//gmd:abstract/gco:CharacterString', namespaces = self.namespaces)
            abstract = '' if abstract is None else abstract.text
        except:
            abstract = ''
            
        return(abstract)

    def linkages(self):
        linkage = self.xml_doc.find('.//{*}linkage/{*}URL', namespaces = self.namespaces)
        if linkage is not None:
            linkage = linkage.text
        
        return(linkage)
    
    def data_links(self):
        dd = {}        
        dfs = self.xml_doc.findall('.//gmd:MD_Format/gmd:name/gco:CharacterString', namespaces = self.namespaces)
        dus = self.xml_doc.findall('.//gmd:onLine/gmd:CI_OnlineResource/gmd:linkage/gmd:URL', namespaces =  self.namespaces)

        if dfs is not None:
            for i,j in enumerate(dfs):
                if j.text in dd.keys():
                    dd[j.text].append(dus[i].text)
                else:
                    dd[j.text] = [dus[i].text]
                    
        return(dd)

## ==============================================
## ==============================================
class Fetch:
    def __init__(self, url=None, callback=lambda: False, verbose=None, headers=r_headers, verify=True):
        self.url = url
        self.callback = callback
        self.verbose = verbose
        self.headers = headers
        self.verify = verify

    def fetch_req(self, params=None, tries=5, timeout=2, read_timeout=10):
        """fetch src_url and return the requests object"""
        
        if tries <= 0:
            utils.echo_error_msg('max-tries exhausted {}'.format(self.url))
            return(None)
            #raise ConnectionError('Maximum attempts at connecting have failed.')
        
        try:
            req = requests.get(
                self.url, stream=True, params=params, timeout=(timeout,read_timeout), headers=self.headers, verify=self.verify
            )
        except:
            req = self.fetch_req(params=params, tries=tries - 1, timeout=timeout + 1, read_timeout=read_timeout + 10)

        if req.status_code == 416:
            if 'Range' in self.headers.keys():
                del self.headers['Range']

                req = self.fetch_req(params=params, tries=tries - 1, timeout=timeout + 1, read_timeout=read_timeout + 10)
                    
        if req is not None and req.status_code != 200:
            utils.echo_error_msg('request from {} returned {}'.format(req.url, req.status_code))
            req = None
            
        return(req)
            
    def fetch_html(self, timeout=2):
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
            req = self.fetch_req(timeout=timeout, read_timeout=read_timeout)
            results = lxml.etree.fromstring(req.text.encode('utf-8'))
        except:
            utils.echo_error_msg('could not access {}'.format(self.url))
        return(results)

    def fetch_file(
            self, dst_fn, params=None, datatype=None, overwrite=False,
            timeout=10, read_timeout=300, tries=5, check_size=True
    ):
        """fetch src_url and save to dst_fn"""

        def retry(resume=False):
            if 'Range' in self.headers:
                del self.headers['Range']

            self.fetch_file(
                dst_fn,
                params=params,
                datatype=datatype,
                overwrite=overwrite,
                timeout=timeout+5,
                read_timeout=read_timeout+50,
                tries=tries-1
            )
        
        status = 0
        dst_fn_size = 0
        if 'Range' in self.headers:
            del self.headers['Range']

        req = None
        if not os.path.exists(os.path.dirname(dst_fn)):
            try:
                os.makedirs(os.path.dirname(dst_fn))
            except: pass
            
        try:
            try:
                if not overwrite and os.path.exists(dst_fn):
                    if not check_size:
                        raise UnboundLocalError('{} exists, '.format(dst_fn))
                    else:
                        dst_fn_size = os.stat(dst_fn).st_size
                        resume_byte_pos = dst_fn_size
                        self.headers['Range'] = 'bytes={}-'.format(resume_byte_pos)
            except OSError:
                pass

            with requests.get(self.url, stream=True, params=params, headers=self.headers,
                              timeout=(timeout,read_timeout), verify=self.verify) as req:
                req_h = req.headers
                if 'Content-Length' in req_h:
                    req_s = int(req_h['Content-Length'])
                else:
                    req_s = -1

                try:
                    if not overwrite and check_size and req_s == os.path.getsize(dst_fn):
                        raise UnboundLocalError('{} exists, '.format(dst_fn))
                    elif req_s == -1 or req_s == 0 or req_s == 49:
                        raise UnboundLocalError('{} exists, '.format(dst_fn))
                except OSError:
                    pass

                if req.status_code == 300:
                    pass

                ## ==============================================
                ## hack for earthdata credential redirect...
                ## recursion here may never end with incorrect user/pass
                ## ==============================================
                timed_out = False
                if req.status_code == 401:
                    ## ==============================================
                    ## we're hoping for a redirect url here.
                    ## ==============================================
                    if self.url == req.url:
                        raise UnboundLocalError('Incorrect Authentication')

                    Fetch(url=req.url, headers=self.headers, verbose=self.verbose).fetch_file(
                        dst_fn,
                        params=params,
                        datatype=datatype,
                        overwrite=overwrite,
                        timeout=timeout,
                        read_timeout=read_timeout
                    )

                elif req.status_code == 200 or req.status_code == 206:
                    curr_chunk = 0
                    total_size = int(req.headers.get('content-length', 0))
                    with open(dst_fn, 'ab' if req.status_code == 206 else 'wb') as local_file:
                        with tqdm(
                                desc='fetching: {}'.format(utils._init_msg(self.url, len('fetching: '), 40)),
                                total=total_size,
                                unit='iB',
                                unit_scale=True
                        ) as pbar:
                            try:
                                for chunk in req.iter_content(chunk_size = 8196):
                                    if self.callback():
                                        break

                                    pbar.update(len(chunk))
                                    if not chunk:
                                        break

                                    local_file.write(chunk)
                                    local_file.flush()
                            except Exception as e:
                                #utils.echo_warning_msg(e)
                                if self.verbose:
                                    utils.echo_warning_msg(
                                        'server returned: {}, and an exception occured: {}, taking a nap and trying again (attempts left: {})...'.format(req.status_code, e, tries)
                                    )

                                    time.sleep(2)
                                    
                                Fetch(url=self.url, headers=self.headers, verbose=self.verbose).fetch_file(
                                    dst_fn,
                                    params=params,
                                    datatype=datatype,
                                    overwrite=overwrite,
                                    timeout=timeout+5,
                                    read_timeout=read_timeout+50,
                                    tries=tries-1
                                )
                                self.verbose=False

                            # if chunk:
                            #     local_file.write(chunk)
                            # except TimeoutError as e:
                            #     timed_out = True

                    if check_size and total_size != os.stat(dst_fn).st_size:
                        raise UnboundLocalError('sizes do not match!')
                        timed_out = True
                        
                elif req.status_code == 429 or req.status_code == 504 or timed_out:

                    #elif tries > 0:
                    # if self.verbose:
                    #     utils.echo_warning_msg('max tries exhausted...')
                    #     status = -1
                    #else:
                    ## ==============================================
                    ## pause a bit and retry...
                    ## ==============================================
                    if self.verbose:
                        utils.echo_warning_msg(
                            'server returned: {}, taking a nap and trying again (attempts left: {})...'.format(req.status_code, tries)
                        )

                    time.sleep(10)
                    Fetch(url=self.url, headers=self.headers, verbose=self.verbose).fetch_file(
                        dst_fn,
                        params=params,
                        datatype=datatype,
                        overwrite=overwrite,
                        timeout=timeout+5,
                        read_timeout=read_timeout+50,
                        tries=tries-1
                    )
                    self.verbose=False
                else:
                    if self.verbose:
                        utils.echo_error_msg('server returned: {} ({})'.format(req.status_code, req.url))
                        raise UnboundLocalError(req.status_code)
                    status = -1

        # except requests.exceptions.Timeout:
        #     utils.echo_error_msg('Connection Timed Out!')
            
        except UnboundLocalError as e:
            #utils.echo_error_msg(e)
            pass
            
        except Exception as e:
            utils.echo_error_msg(e)
            status = -1
                
        if not os.path.exists(dst_fn) or os.stat(dst_fn).st_size ==  0:
            #utils.echo_error_msg('data not fetched...')
            status = -1
            if check_size:
                raise UnboundLocalError('data not fetched')
            #status = -1
                        
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

## ==============================================
## fetch queue for threads
## ==============================================
def fetch_queue(q, c=True):
    """fetch queue `q` of fetch results

    each fetch queue `q` should be a list of the following:
    [remote_data_url, local_data_path, data-type, fetches-module]    

    set c to False to skip size-checking
    """

    while True:
        fetch_args = q.get()
        if not fetch_args[3].callback():
            if not os.path.exists(os.path.dirname(fetch_args[1])):
                try:
                    os.makedirs(os.path.dirname(fetch_args[1]))
                except: pass
            
            if fetch_args[0].split(':')[0] == 'ftp':
                Fetch(
                    url=fetch_args[0],
                    callback=fetch_args[3].callback,
                    verbose=fetch_args[3].verbose,
                    headers=fetch_args[3].headers
                ).fetch_ftp_file(fetch_args[1])
            else:
                try:
                    Fetch(
                        url=fetch_args[0],
                        callback=fetch_args[3].callback,
                        verbose=fetch_args[3].verbose,
                        headers=fetch_args[3].headers,
                        verify=False if fetch_args[2] == 'srtm' or fetch_args[2] == 'mar_grav' else True
                    ).fetch_file(fetch_args[1], check_size=c)
                except:
                    utils.echo_warning_msg('fetch of {} failed...'.format(fetch_args[0]))
                        
        q.task_done()
        
class fetch_results(threading.Thread):
    """fetch results gathered from a fetch module.

    results is a list of URLs with data type

    when a fetch module is run with {module}.run() it will fill {module}.results with a list of urls, e.g.
    {module}.results = [[http://data/url.xyz.gz, /home/user/data/url.xyz.gz, data-type], ...]
    where each result in is [data_url, data_fn, data_type]

    run this on an initialized fetches module:
    >>> fetch_result(fetches_module, n_threads=3).run()
    and this will fill a queue for data fetching, using 'n_threads' threads.
    """
    
    def __init__(self, mod, check_size=True, n_threads=3):
        threading.Thread.__init__(self)
        self.fetch_q = queue.Queue()
        self.mod = mod
        self.check_size = check_size
        self.n_threads = n_threads
        if len(self.mod.results) == 0:
            self.mod.run()
        
    def run(self):
        for _ in range(self.n_threads):
            t = threading.Thread(target=fetch_queue, args=(self.fetch_q, self.check_size))
            t.daemon = True
            t.start()
        for row in self.mod.results:
            self.fetch_q.put([row[0], os.path.join(self.mod._outdir, row[1]), row[2], self.mod])
            
        self.fetch_q.join()

## ==============================================
## Fetch Modules
## ==============================================
class FetchModule:
    def __init__(self, src_region = None, callback = lambda: False, verbose = True,
                 outdir = None, name = 'fetches', params = {}):
        self.region = src_region
        self.callback = callback
        self.verbose = verbose
        self.outdir = outdir
        self.params = params
        self.status = 0
        self.results = []
        self.name = name
        #self.headers = { 'User-Agent': 'Fetches v%s' %(fetches.__version__) }
        self.headers = { 'User-Agent': 'Mozilla/5.0 (X11; Linux x86_64; rv:45.0) Gecko/20100101 Firefox/45.0' }

        if self.outdir is None:
            self._outdir = os.path.join(os.getcwd(), self.name)
        else:
            self._outdir = self.outdir

        ## for dlim support
        self.data_format = None
        self.src_srs = None

        if self.region is None or not self.region.valid_p():
            self.region = regions.Region().from_list([-180,180,-90,90])
        
    def run(self):
        raise(NotImplementedError)

    def fetch(self, entry, check_size = True):
        ## make sure this fetches correctly and retry if it doesn't!
        try:
            status = Fetch(
                entry[0],
                callback=self.callback,
                verbose=self.verbose,
                headers=self.headers,
            ).fetch_file(os.path.join(self._outdir, entry[1]), timeout=5, read_timeout=5, check_size=check_size)
        except:
            status = -1

        return(status)

    def fetch_results(self):
        for entry in self.results:
            self.fetch(entry)

## ==============================================
## GMRT
## ==============================================
class GMRT(FetchModule):
    """The Global Multi-Resolution Topography synthesis.
    
The Global Multi-Resolution Topography (GMRT) synthesis is a multi-resolutional 
compilation of edited multibeam sonar data collected by scientists and institutions worldwide, that is 
reviewed, processed and gridded by the GMRT Team and merged into a single continuously updated compilation 
of global elevation data. The synthesis began in 1992 as the Ridge Multibeam Synthesis (RMBS), was expanded 
to include multibeam bathymetry data from the Southern Ocean, and now includes bathymetry from throughout 
the global and coastal oceans.

Data Formats
    GMT v3 Compatible NetCDF (GMT id=cf)
    COARDS/CF1.6 Compliant NetCDF (GMT id=nd)
    ESRI ArcASCII
    GeoTIFF
Metadata Formats
    XML (metadata)
    JSON (metadata)
    Plain text (metadata)
    
layers: 'topo' or 'topo-mask'
fmt: 'geotiff', 'netcdf'
    
Data is assumed instantaneous MSL (5773?)

https://www.gmrt.org

< gmrt:res=max:fmt=geotiff:layer=topo >
    """
    
    def __init__(self, res = 'default', fmt = 'geotiff', layer = 'topo', want_swath = False, **kwargs):
        super().__init__(name='gmrt', **kwargs) 
        self.res = res
        self.fmt = fmt
        self.want_swath = want_swath
        if layer != 'topo' and layer != 'topo-mask':
            self.layer = 'topo'
        else:
            self.layer = layer
            
        self._gmrt_grid_url = "https://www.gmrt.org:443/services/GridServer?"
        self._gmrt_grid_urls_url = "https://www.gmrt.org:443/services/GridServer/urls?"
        self._gmrt_grid_metadata_url = "https://www.gmrt.org/services/GridServer/metadata?"        
        self._gmrt_swath_poly_url = "https://www.gmrt.org/shapefiles/gmrt_swath_polygons.zip"
        
        self.gmrt_region = self.region.copy()
        self.gmrt_region.buffer(pct=2.33,x_inc=.0088,y_inc=.0088)
        self.gmrt_region._wgs_extremes(just_below=True)

        self.data_format = 200
        self.src_srs = 'epsg:4326+3855'

        self.headers = { 'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:89.0) Gecko/20100101 Firefox/89.0'}
        
    def run(self):
        '''Run the GMRT fetching module'''

        if self.region is None:
            return([])

        self.data = {'north':self.gmrt_region.ymax, 'west':self.gmrt_region.xmin,
                     'south':self.gmrt_region.ymin, 'east':self.gmrt_region.xmax,
                     'mformat':'json', 'resolution':self.res, 'format':self.fmt,
                     'layer':self.layer}
        req = Fetch(self._gmrt_grid_url).fetch_req(params=self.data, tries=10, timeout=2)
        if req is not None:
            outf = 'gmrt_{}_{}_{}.{}'.format(self.layer, self.res, self.region.format('fn_full'), 'tif' if self.fmt == 'geotiff' else 'grd')
            self.results.append([req.url, outf, 'gmrt']) 
        else:            
            gmrt_urls = req.json()
            for url in gmrt_urls:
                if self.layer == 'topo-mask':
                    url = url.replace('topo', 'topo-mask')

                opts = {}
                for url_opt in url.split('?')[1].split('&'):
                    opt_kp = url_opt.split('=')
                    opts[opt_kp[0]] = opt_kp[1]

                url_region = regions.Region().from_list([
                    float(opts['west']),
                    float(opts['east']),
                    float(opts['south']),
                    float(opts['north'])
                ])
                outf = 'gmrt_{}_{}.{}'.format(opts['layer'], url_region.format('fn'), 'tif' if self.fmt == 'geotiff' else 'grd')
                self.results.append([url, outf, 'gmrt'])
                if self.want_swath:
                    self.results.append([self._gmrt_swath_poly_url, 'gmrt_swath_polygons.zip', 'gmrt'])

                
        return(self)

## ==============================================
## GEBCO
## ==============================================
class GEBCO(FetchModule):
    """GEBCO
    
GEBCO’s current gridded bathymetric data set, the GEBCO_2022 Grid, is a global terrain model for ocean and land, 
providing elevation data, in meters, on a 15 arc-second interval grid. It is accompanied by a Type Identifier 
(TID) Grid that gives information on the types of source data that the GEBCO_2022 Grid is based. 

https://www.gebco.net

< gebco >
"""
    
    def __init__(self, want_ice = 'geotiff', want_sub_ice = False, want_tid = False,
                 exclude_tid = None, upper_limit = None, lower_limit = None, **kwargs):
        super().__init__(name='gebco', **kwargs) 
        self._gebco_urls = {
            'gebco_ice': {
                'netcdf': 'https://www.bodc.ac.uk/data/open_download/gebco/gebco_2022/zip/',
                'geotiff': 'https://www.bodc.ac.uk/data/open_download/gebco/gebco_2022/geotiff/',
                'ascii': 'https://www.bodc.ac.uk/data/open_download/gebco/gebco_2022/esri_ascii/'
            },
            'gebco_sub_ice': {
                'netcdf': 'https://www.bodc.ac.uk/data/open_download/gebco/gebco_2022_sub_ice_topo/zip/',
                'geotiff': 'https://www.bodc.ac.uk/data/open_download/gebco/gebco_2022_sub_ice_topo/geotiff/',
                'ascii': 'https://www.bodc.ac.uk/data/open_download/gebco/gebco_2022_sub_ice_topo/esri_ascii/'
            },
            'gebco_tid': {
                'netcdf': 'https://www.bodc.ac.uk/data/open_download/gebco/gebco_2022_tid/zip/',
                'geotiff': 'https://www.bodc.ac.uk/data/open_download/gebco/gebco_2022_tid/geotiff/',
                'ascii': 'https://www.bodc.ac.uk/data/open_download/gebco/gebco_2022_tid/esri_ascii/'
            }
        }

        self.want_ice = utils.str_or(want_ice, 'geotiff') if want_ice else False
        self.want_sub_ice = utils.str_or(want_sub_ice, 'geotiff') if want_sub_ice else False
        self.want_tid = utils.str_or(want_tid, 'geotiff') if want_tid else False

        exclude_tid = utils.str_or(exclude_tid)
        self.exclude_tid = []
        if exclude_tid is not None:
            for tid_key in exclude_tid.split('/'):
                self.exclude_tid.append(utils.int_or(tid_key))
                
        try:
            self.exclude_tid.remove(None)
        except:
            pass

        self.tid_dic = {
            0: ['Land', .21],
            10: ['Singlebeam - depth value collected by a single beam echo-sounder', .65],
            11:	['Multibeam - depth value collected by a multibeam echo-sounder', .95],
            12:	['Seismic - depth value collected by seismic methods', .3],
            13:	['Isolated sounding - depth value that is not part of a regular survey or trackline', .4],
            14:	['ENC sounding - depth value extracted from an Electronic Navigation Chart (ENC)', .5],
            15:	['Lidar - depth derived from a bathymetric lidar sensor', 1],
            16:	['Depth measured by optical light sensor', .3],
            17:	['Combination of direct measurement methods', .2],
            40:	['Predicted based on satellite-derived gravity data - depth value is an interpolated value guided by satellite-derived gravity data', .19],
            41:	['Interpolated based on a computer algorithm - depth value is an interpolated value based on a computer algorithm (e.g. Generic Mapping Tools)', .18],
            42:	['Digital bathymetric contours from charts - depth value taken from a bathymetric contour data set', .17],
            43:	['Digital bathymetric contours from ENCs - depth value taken from bathymetric contours from an Electronic Navigation Chart (ENC)', .16],
            44:	['Bathymetric sounding - depth value at this location is constrained by bathymetric sounding(s) within a gridded data set where interpolation between sounding points is guided by satellite-derived gravity data', .15],
            45:	['Predicted based on helicopter/flight-derived gravity data', .14],
            46:	['Depth estimated by calculating the draft of a grounded iceberg using satellite-derived freeboard measurement.', .13],
            70:	['Pre-generated grid - depth value is taken from a pre-generated grid that is based on mixed source data types, e.g. single beam, multibeam, interpolation etc.', .12],
            71:	['Unknown source - depth value from an unknown source', .11],
            72:	['Steering points - depth value used to constrain the grid in areas of poor data coverage', .1]
        }

        self.gebco_region = self.region.copy()
        self.gebco_region.zmax = utils.float_or(upper_limit)
        self.gebco_region.zmin = utils.float_or(lower_limit)

        self.data_format = -2
        self.src_srs = 'epsg:4326+3855'
        
    def run(self):
        '''Run the GEBCO fetching module'''

        outf = 'gebco.zip'
        if self.want_ice:
            self.results.append([self._gebco_urls['gebco_ice'][self.want_ice], 'gebco_ice.zip', 'gebco'])
        if self.want_sub_ice:
            self.results.append([self._gebco_urls['gebco_sub_ice'][self.want_sub_ice], 'gebco_sub_ice.zip', 'gebco'])
        if self.want_tid:
            self.results.append([self._gebco_urls['gebco_tid'][self.want_tid], 'gebco_tid.zip', 'gebco'])
                
        return(self)

## ==============================================
## Copernicus
## ==============================================
class CopernicusDEM(FetchModule):
    """COPERNICUS sattelite elevation data
    
The Copernicus DEM is a Digital Surface Model (DSM) which represents the surface of the Earth including buildings, 
infrastructure and vegetation.

datatype of 1 is 10 m and datatype of 3 is 30 m
    
https://ec.europa.eu/eurostat/web/gisco/geodata/reference-data/elevation/copernicus-dem/elevation

< copernicus >"""
    
    def __init__(self, where='', datatype=None, **kwargs):
        super().__init__(name='copernicus', **kwargs)
        self.cop30_rurl = 'https://opentopography.s3.sdsc.edu/minio/raster/COP30/COP30_hh/'
        self.cop30_url = 'https://opentopography.s3.sdsc.edu/minio/download/raster/COP30/COP30_hh/'
        self.cop30_vrt_url = 'https://opentopography.s3.sdsc.edu/minio/download/raster/COP30/COP30_hh.vrt?token='
        self.cop_10_url = 'https://gisco-services.ec.europa.eu/dem/copernicus/outD/'
        self.cop_10_aux_url = 'https://gisco-services.ec.europa.eu/dem/copernicus/outA/'
        self.cop_10_web = 'https://ec.europa.eu/eurostat/web/gisco/geodata/reference-data/elevation/copernicus-dem/elevation'
        self.where = [where] if len(where) > 0 else []
        self.datatype = datatype        
        self.FRED = FRED.FRED(name=self.name, verbose=self.verbose)

        self.data_format = -2
        self.src_srs = 'epsg:4326+3855'
        
        self.headers = { 'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:89.0) Gecko/20100101 Firefox/89.0',
                         'referer': 'https://opentopography.s3.sdsc.edu/minio/raster/COP30/COP30_hh/' }
        
        self.update_if_not_in_FRED()

    def update_if_not_in_FRED(self):
        self.FRED._open_ds()
        self.FRED._attribute_filter(["DataSource = '{}'".format(self.name)])
        if len(self.FRED.layer) == 0:
            self.FRED._close_ds()
            self.update()
            
        self.FRED._close_ds()

    def update(self):
        """Crawl the COP30 database and update/generate the COPERNICUS reference vector."""
        
        self.FRED._open_ds(1)
        surveys = []
        page = Fetch(self.cop_10_url, verbose=True).fetch_html()
        rows = page.xpath('//a[contains(@href, ".zip")]/@href')
        with tqdm(
                desc='scanning for COPERNICUS COP-10 datasets',
                leave=self.verbose
        ) as pbar:            
            for i, row in enumerate(rows):
                pbar.update()
                sid = row.split('.')[0]
                self.FRED._attribute_filter(["ID = '{}'".format(sid)])
                if self.FRED.layer is None or len(self.FRED.layer) == 0:
                    spat = row.split('.')[0].split('_')[-1]
                    x = int(spat.split('x')[-1])
                    y = int(spat.split('x')[0].split('y')[-1])
                    this_region = regions.Region().from_list(
                        [x, x + 10, y, y + 10]
                    )
                    geom = this_region.export_as_geom()
                    if geom is not None:
                        surveys.append(
                            {
                                'Name': row.split('.')[0],
                                'ID': sid,
                                'Agency': 'EU',
                                'Date': utils.this_date(),
                                'MetadataLink': self.cop_10_aux_url,
                                'MetadataDate': utils.this_date(),
                                'DataLink': self.cop_10_url + row,
                                'DataType': '3',
                                'DataSource': 'copernicus',
                                'HorizontalDatum': 'epsg:4326',
                                'VerticalDatum': 'msl',
                                'Info': '',
                                'geom': geom
                            }
                        )

        f = Fetch(self.cop30_vrt_url, headers=self.headers, verbose=True)
        page = f.fetch_xml()
        fns = page.findall('.//SourceFilename')
        with tqdm(
                total=len(fns),
                desc='scanning for COPERNICUS COP-30 datasets',
                leave=self.verbose
        ) as pbar:
            for i, fn in enumerate(fns):
                pbar.update()
                sid = fn.text.split('/')[-1].split('.')[0]
                self.FRED._attribute_filter(["ID = '{}'".format(sid)])
                if self.FRED.layer is None or len(self.FRED.layer) == 0:            
                    spat = fn.text.split('_10_')[-1].split('_DEM')[0]
                    xsplit = '_E' if 'E' in spat else '_W'
                    ysplit = 'S' if 'S' in spat else 'N'
                    x = int(spat.split(xsplit)[-1].split('_')[0])
                    y = int(spat.split(xsplit)[0].split(ysplit)[-1].split('_')[0])

                    if xsplit == '_W':
                        x = x * -1
                    if ysplit == 'S':
                        y = y * -1

                    this_region = regions.Region().from_list([x, x + 1, y, y + 1])
                    geom = this_region.export_as_geom()
                    if geom is not None:
                        surveys.append(
                            {
                                'Name': fn.text.split('.')[0].split('/')[-1],
                                'ID': sid,
                                'Agency': 'EU',
                                'Date': utils.this_date(),
                                'MetadataLink': '',
                                'MetadataDate': utils.this_date(),
                                'DataLink': self.cop30_url + fn.text.split('/')[-1] + '?token=',
                                'DataType': '1',
                                'DataSource': 'copernicus',
                                'HorizontalDatum': 'epsg:4326',
                                'Etcetra': self.cop30_rurl,
                                'VerticalDatum': 'msl',
                                'Info': '',
                                'geom': geom
                            }
                        )

        self.FRED._add_surveys(surveys)
        self.FRED._close_ds()

    def run(self):
        '''Run the COPERNICUS DEM fetching module'''

        if self.datatype is not None:
            self.where.append("DataType = '{}'".format(self.datatype))

        _results = FRED._filter_FRED(self)
        with tqdm(
                total=len(_results),
                desc='scanning COPERNICUS datasets',
                leave=self.verbose
        ) as pbar:
            for surv in _results:
                pbar.update()
                for i in surv['DataLink'].split(','):
                    self.results.append([i, i.split('/')[-1].split('?')[0], surv['DataType']])
                    #self.results.append([i, i.split('/')[-1].split('?')[0], surv['DataType']])
                
        return(self)

## ==============================================
## FABDEM
## ==============================================
class FABDEM(FetchModule):
    """FABDEM elevation data
    
FABDEM (Forest And Buildings removed Copernicus DEM) is a global elevation map that removes building and tree height
biases from the Copernicus GLO 30 Digital Elevation Model (DEM). The data is available at 1 arc second
grid spacing (approximately 30m at the equator) for the globe.

https://data.bris.ac.uk/data/dataset/s5hqmjcdj8yo2ibzi9b4ew3sn

< fabdem >"""
    
    def __init__(self, **kwargs):
        super().__init__(name='fabdem', **kwargs)
        self._fabdem_footprints_url = 'https://data.bris.ac.uk/datasets/s5hqmjcdj8yo2ibzi9b4ew3sn/FABDEM_v1-2_tiles.geojson'
        self._fabdem_info_url = 'https://data.bris.ac.uk/data/dataset/s5hqmjcdj8yo2ibzi9b4ew3sn'
        self._fabdem_data_url = 'https://data.bris.ac.uk/datasets/s5hqmjcdj8yo2ibzi9b4ew3sn'
        self.headers = { 'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:89.0) Gecko/20100101 Firefox/89.0',
                         'referer': self._fabdem_info_url }

        self.data_format = -2
        self.src_srs = 'epsg:4326+3855'
        
    def run(self):
        v_json = os.path.basename(self._fabdem_footprints_url)
        try:
            status = Fetch(self._fabdem_footprints_url, verbose=self.verbose).fetch_file(v_json)
            v_ds = ogr.Open(v_json)
        except:
            v_ds = None
            status = -1
        if v_ds is not None:
            layer = v_ds.GetLayer()
            fcount = layer.GetFeatureCount()
            for f in range(0, fcount):
                feature = layer[f]
                geom = feature.GetGeometryRef()
                if geom.Intersects(self.region.export_as_geom()):
                    zipfile_name = feature.GetField('zipfile_name')
                    zipfile_url = '/'.join([self._fabdem_data_url, zipfile_name])
                    if zipfile_url not in [x[0] for x in self.results]:
                        self.results.append([zipfile_url, zipfile_name, 'raster']
                    )
            v_ds = None
                        
        utils.remove_glob(v_json)
        
class FABDEM_FRED(FetchModule):
    """Fetch FABDEM data"""
    
    def __init__(self, where='', **kwargs):
        super().__init__(name='fabdem', **kwargs)
        self._fabdem_footprints_url = 'https://data.bris.ac.uk/datasets/s5hqmjcdj8yo2ibzi9b4ew3sn/FABDEM_v1-2_tiles.geojson'
        self._fabdem_info_url = 'https://data.bris.ac.uk/data/dataset/s5hqmjcdj8yo2ibzi9b4ew3sn'
        self.where = [where] if len(where) > 0 else []
        self.FRED = FRED.FRED(name=self.name, verbose=self.verbose)
        self.update_if_not_in_FRED()
        self.data_format = -2
        self.src_srs = 'epsg:4326+3855'
        
    def update_if_not_in_FRED(self):
        self.FRED._open_ds()
        self.FRED._attribute_filter(["DataSource = '{}'".format(self.name)])
        if len(self.FRED.layer) == 0:
            self.FRED._close_ds()
            self.update()
        self.FRED._close_ds()
        
    def update(self):
        self.FRED._open_ds()
        v_json = os.path.basename(self._fabdem_footprints_url)
        try:
            status = Fetch(self._fabdem_footprints_url, verbose=self.verbose).fetch_file(v_json)
        except:
            status = -1
            
        shp_regions = regions.gdal_ogr_regions(v_json)
        shp_region = regions.Region()
        for this_region in shp_regions:
            if shp_region.valid_p(check_xy=True):
                shp_region = regions.regions_merge(shp_region, this_region)
            else: shp_region = this_region
        geom = shp_region.export_as_geom()
        
        self.FRED._attribute_filter(["ID = '{}'".format('FABDEM-1')])
        if self.FRED.layer is None or len(self.FRED.layer) == 0:
            self.FRED._add_survey(Name = 'FABDEM', ID = 'FABDEM-1', Agency = 'Univesity of Bristol', Date = utils.this_year(),
                                  MetadataLink = self._fabdem_info_url, MetadataDate = utils.this_year(),
                                  DataLink = self._fabdem_footprints_url, IndexLink = self._fabdem_footprints_url,
                                  DataType = 'raster', DataSource = 'fabdem', Info = 'Bare Earth Copernicus', geom = geom)
        utils.remove_glob(v_json)
        self.FRED._close_ds()

    def run(self):
        for surv in FRED._filter_FRED(self):
            v_json = os.path.basename(self._fabdem_footprints_url)
            status = Fetch(surv['IndexLink']).fetch_file(v_json, verbose=self.verbose)
            try:
                v_ds = ogr.Open(v_json)
            except:
                v_ds = None
                status = -1
                
            if v_ds is not None:
                layer = v_ds.GetLayer()
                fcount = layer.GetFeatureCount()
                for f in range(0, fcount):
                    feature = layer[f]
                    geom = feature.GetGeometryRef()
                    if geom.Intersects(self.region.export_as_geom()):
                        zipfile_name = feature.GetField('zipfile_name')
                        self.results.append(['/'.join([self._fabdem_data_url, zipfile_name]), zipfile_name, 'raster'])
            utils.remove_glob(v_zip)

## ==============================================
## NASADEM
## ==============================================
class NASADEM(FetchModule):
    """NASA Digital Elevation Model
    
Our objective is to provide the scientific and civil communities with a state-of-the-art global 
digital elevation model (DEM) derived from a combination of Shuttle Radar Topography Mission (SRTM) 
processing improvements, elevation control, void-filling and merging with data unavailable at the 
time of the original SRTM production.

https://www.earthdata.nasa.gov/esds/competitive-programs/measures/nasadem

< nasadem >"""
    
    def __init__(self, where='', datatype=None, **kwargs):
        super().__init__(name='nasadem', **kwargs)
        self.nasadem_rurl = 'https://opentopography.s3.sdsc.edu/minio/raster/NASADEM/NASADEM_be/'
        self.nasadem_url = 'https://opentopography.s3.sdsc.edu/minio/download/raster/NASADEM/NASADEM_be/'
        self.nasadem_vrt_url = 'https://opentopography.s3.sdsc.edu/minio/download/raster/NASADEM/NASADEM_be.vrt?token='
        self.where = [where] if len(where) > 0 else []
        self.datatype = datatype
        self.FRED = FRED.FRED(name=self.name, verbose=self.verbose)

        self.headers = { 'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:89.0) Gecko/20100101 Firefox/89.0',
                         'referer': 'https://opentopography.s3.sdsc.edu/minio/raster/NASADEM/NASADEM_be/' }
        
        self.update_if_not_in_FRED()
        self.data_format = 200
        self.src_srs = 'epsg:4326+3855'

    def update_if_not_in_FRED(self):
        self.FRED._open_ds()
        self.FRED._attribute_filter(["DataSource = '{}'".format(self.name)])
        if len(self.FRED.layer) == 0:
            self.FRED._close_ds()
            self.update()
        self.FRED._close_ds()

    def update(self):
        """Crawl the COP30 database and update/generate the NASADEM reference vector."""
        
        self.FRED._open_ds(1)
        surveys = []                    
        f = Fetch(self.nasadem_vrt_url, headers=self.headers, verbose=True)
        page = f.fetch_xml()
        fns = page.findall('.//SourceFilename')

        with tqdm(total=len(fns), desc='scanning NASADEM datasets', leave=self.verbose) as pbar:        
            for i, fn in enumerate(fns):
                sid = fn.text.split('/')[-1].split('.')[0]
                pbar.update()
                self.FRED._attribute_filter(["ID = '{}'".format(sid)])
                if self.FRED.layer is None or len(self.FRED.layer) == 0:
                    spat = fn.text.split('_HGT_')[-1].split('.')[0]
                    xsplit = 'e' if 'e' in spat else 'w'
                    ysplit = 's' if 's' in spat else 'n'
                    x = int(spat.split(xsplit)[-1])
                    y = int(spat.split(xsplit)[0].split(ysplit)[-1])

                    if xsplit == 'w':
                        x = x * -1
                    if ysplit == 's':
                        y = y * -1

                    this_region = regions.Region().from_list([x, x + 1, y, y + 1])
                    geom = this_region.export_as_geom()
                    if geom is not None:
                        surveys.append(
                            {
                                'Name': fn.text.split('.')[0].split('/')[-1],
                                'ID': sid,
                                'Agency': 'NASA',
                                'Date': utils.this_date(),
                                'MetadataLink': '',
                                'MetadataDate': utils.this_date(),
                                'DataLink': self.nasadem_url + fn.text.split('/')[-1] + '?token=',
                                'DataType': '1',
                                'DataSource': 'nasadem',
                                'HorizontalDatum': 4326,
                                'Etcetra': self.nasadem_rurl,
                                'VerticalDatum': 'msl',
                                'Info': '',
                                'geom': geom
                            }
                        )

        self.FRED._add_surveys(surveys)
        self.FRED._close_ds()

    def run(self):
        '''Run the NASADEM DEM fetching module'''

        for surv in FRED._filter_FRED(self):
            for i in surv['DataLink'].split(','):
                self.results.append([i, i.split('/')[-1].split('?')[0], surv['DataType']])
        return(self)
            
## ==============================================
## MarGrav - Marine Gravity
## ==============================================
class MarGrav(FetchModule):
    """MARine GRAVity Satellite Altimetry Topography from Scripps.

Fetch mar_grav sattelite altimetry topography'''
https://topex.ucsd.edu/WWW_html/mar_grav.html
ftp://topex.ucsd.edu/pub/global_grav_1min/
https://topex.ucsd.edu/marine_grav/explore_grav.html
https://topex.ucsd.edu/marine_grav/white_paper.pdf
    
< mar_grav:upper_limit=None:lower_limit=None:raster=False:mag=1 >"""
    
    def __init__(self, mag=1, upper_limit=None, lower_limit=None, raster=False, **kwargs):
        super().__init__(name='mar_grav', **kwargs) 
        self._mar_grav_url = 'https://topex.ucsd.edu/cgi-bin/get_data.cgi'
        self.mag = mag if mag == 1 else 0.1
        self.grav_region = self.region.copy()
        self.grav_region._wgs_extremes(just_below=True)
        self.grav_region.zmax = utils.float_or(upper_limit)
        self.grav_region.zmin = utils.float_or(lower_limit)
        self.raster = raster

        self.data_format = '168:x_offset=REM:skip=1'
        self.src_srs = 'epsg:4326+3855'
        self.headers = { 'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:89.0) Gecko/20100101 Firefox/89.0'}
        
    def run(self):
        '''Run the mar_grav fetching module.'''
        
        if self.region is None:
            return([])
        
        _data = {'north':self.region.ymax, 'west':self.region.xmin,
                 'south':self.region.ymin, 'east':self.region.xmax,
                 'mag':self.mag}
        _req = Fetch(self._mar_grav_url, verify=False).fetch_req(params=_data)
        if _req is not None:
            outf = 'mar_grav_{}.xyz'.format(self.region.format('fn_full'))
            self.results.append([_req.url, outf, 'mar_grav'])
            
## ==============================================
## SRTM Plus
## ==============================================
class SRTMPlus(FetchModule):
    """SRTM15+: GLOBAL BATHYMETRY AND TOPOGRAPHY AT 15 ARCSECONDS.

https://topex.ucsd.edu/WWW_html/srtm15_plus.html
http://topex.ucsd.edu/sandwell/publications/180_Tozer_SRTM15+.pdf
https://topex.ucsd.edu/pub/srtm15_plus/
https://topex.ucsd.edu/pub/srtm15_plus/SRTM15_V2.3.nc
    
< srtm_plus >
"""
    
    def __init__(self, **kwargs):
        super().__init__(name='srtm_plus', **kwargs) 
        self._srtm_url = 'https://topex.ucsd.edu/cgi-bin/get_srtm15.cgi'

        self.data_format = '168:skip=1'
        self.src_srs = 'epsg:4326+3855'

    def run(self):
        '''Run the SRTM fetching module.'''
        
        if self.region is None: return([])

        self.data = {'north':self.region.ymax, 'west':self.region.xmin,
                     'south':self.region.ymin, 'east':self.region.xmax}

        _req = Fetch(self._srtm_url, verify=False).fetch_req(params=self.data)
        if _req is not None:
            outf = 'srtm_{}.xyz'.format(self.region.format('fn'))
            self.results.append([_req.url, outf, 'srtm'])
            
## ==============================================
## Charts - ENC/RNC
## ==============================================
## the arcgis rest server doesn't filter by bbox for some reason, always returns all data
class Charts(FetchModule):
    """NOAA Nautical CHARTS

https://www.charts.noaa.gov/

< charts >"""
    
    def __init__(self, where = '1=1', want_rnc = False, **kwargs):
        super().__init__(name='charts', **kwargs)
        self._charts_url = 'https://gis.charttools.noaa.gov/arcgis/rest/services/MCS/ENCOnline/MapServer/exts/MaritimeChartService/MapServer'
        self._charts_query_url = '{0}/queryDatasets?'.format(self._charts_url)
        self.where = where
        self.want_rnc = want_rnc
        
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
            'returnGeometry':'True',
        }
        print(_data)
        print(self.region.format('bbox'))
        _req = Fetch(self._charts_query_url, verbose=self.verbose).fetch_req(params=_data)
        if _req is not None:
            _req_json = _req.json()
            print(len(_req_json['charts']))
            print(_req_json['charts'][0])
            #print(_req.text)
            
class NauticalCharts(FetchModule):
    """NOAA Nautical CHARTS

Fetch digital chart data from NOAA

set the 'want_rnc' flag to True to fetch RNC along with ENC data

https://www.charts.noaa.gov/

< charts:want_rnc=False >"""

    def __init__(self, where='', want_rnc = False, **kwargs):
        super().__init__(name='charts', **kwargs)
        self._charts_url = 'https://www.charts.noaa.gov/'
        self._enc_data_catalog = 'https://charts.noaa.gov/ENCs/ENCProdCat_19115.xml'
        self._rnc_data_catalog = 'https://charts.noaa.gov/RNCs/RNCProdCat_19115.xml'
        self._ienc_charts_data_catalog = 'https://ienccloud.us/ienc/products/catalog/IENCU37ProductsCatalog.xml'
        self._ienc_buoys_data_catalog = 'https://ienccloud.us/ienc/products/catalog/IENCBuoyProductsCatalog.xml'
        self._urls = [self._enc_data_catalog, self._rnc_data_catalog]
        self._outdir = os.path.join(os.getcwd(), 'charts')
        self._dt_xml = {'ENC':self._enc_data_catalog, 'RNC':self._rnc_data_catalog}
        self.where = [where] if len(where) > 0 else []
        #self.datatype = datatype # 'enc' or 'rnc'
        self.want_rnc = want_rnc
        self.name = 'charts'
        self.v_datum = 'mllw'
        #self.data_format = 302
        self.src_srs='epsg:4326+5866'
        
        self.FRED = FRED.FRED(name=self.name, verbose=self.verbose)
        self.update_if_not_in_FRED()
        
    def update_if_not_in_FRED(self):
        self.FRED._open_ds()
        self.FRED._attribute_filter(["DataSource = '{}'".format(self.name)])
        if len(self.FRED.layer) == 0:
            self.FRED._close_ds()
            self.update()
            
        self.FRED._close_ds()

    def update(self):
        """Update or create the reference vector file"""
        
        self.FRED._open_ds(1)
        for dt in self._dt_xml.keys():
            surveys = []
            this_xml = iso_xml(self._dt_xml[dt], timeout=1000, read_timeout=2000)
            charts = this_xml.xml_doc.findall('.//{*}has', namespaces = this_xml.namespaces)
            with tqdm(
                    total=len(charts),
                    desc='scanning for CHARTS ({}) datasets'.format(dt),
                    leave=self.verbose
            ) as pbar:
                for i, chart in enumerate(charts):
                    pbar.update(1)
                    this_xml.xml_doc = chart
                    title = this_xml.title()
                    self.FRED._attribute_filter(["ID = '{}'".format(title)])
                    if self.FRED.layer is None or len(self.FRED.layer) == 0:
                        h_epsg, v_epsg = this_xml.reference_system()
                        this_data = this_xml.linkages()
                        geom = this_xml.polygon(geom=True)
                        if geom is not None:
                            surveys.append(
                                {
                                    'Name': title,
                                    'ID': title,
                                    'Agency': 'NOAA',
                                    'Date': this_xml.date(),
                                    'MetadataLink': this_xml.url,
                                    'MetadataDate': this_xml.xml_date(),
                                    'DataLink': this_data,
                                    'Link': self._charts_url,
                                    'DataType': dt,
                                    'DataSource': 'charts',
                                    'HorizontalDatum': h_epsg,
                                    'VerticalDatum': v_epsg,
                                    'Info': this_xml.abstract,
                                    'geom': geom
                                }
                            )
                        
            self.FRED._add_surveys(surveys)
                
        self.FRED._close_ds()

    def generate_tidal_vdatum(self, src_vdatum, dst_vdatum):
        self.vdatum_grid = vdatums._tidal_transform(self.region, src_vdatum, dst_vdatum)
        
    def run(self):
        """Search for data in the reference vector file"""

        #if self.datatype is not None:
        if self.want_rnc:
            self.where.append("DataType = 'RNC'")
        else:
            self.where.append("DataType = 'ENC'")

        _results = FRED._filter_FRED(self)
        with tqdm(
                total=len(_results),
                desc='scanning CHARTS datasets',
                leave=self.verbose
        ) as pbar:
            for surv in _results:
                pbar.update(1)
                for i in surv['DataLink'].split(','):
                    self.results.append([i, i.split('/')[-1], surv['DataType']])
    
## ==============================================
## NCEI Multibeam
## ==============================================
## MapServer testing
class MBDB(FetchModule):
    """NOSHydro"""
    
    def __init__(self, where='1=1', **kwargs):
        super().__init__(name='multibeam', **kwargs)
        #self._mb_dynamic_url = 'https://gis.ngdc.noaa.gov/arcgis/rest/services/web_mercator/multibeam_dynamic/MapServer/0'
        #self._mb_url = 'https://gis.ngdc.noaa.gov/arcgis/rest/services/web_mercator/multibeam/MapServer/0'
        #self._nos_data_url = 'https://data.ngdc.noaa.gov/platforms/ocean/nos/coast/'
        self._mb_dynamic_url = 'https://gis.ngdc.noaa.gov/arcgis/rest/services/multibeam_footprints/MapServer/0'
        self._mb_query_url = '{0}/query?'.format(self._mb_dynamic_url)
        self.where = where
        
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
        _req = Fetch(self._mb_query_url, verbose=self.verbose).fetch_req(params=_data)
        if _req is not None:
            #print(_req.text)
            #print(_req.url)
            features = _req.json()
            for feature in features['features']:
                #pass
                print(feature)

class Multibeam(FetchModule):
    """NOAA MULTIBEAM bathymetric data.

Fetch multibeam data from NOAA NCEI
        
NCEI is the U.S. national archive for multibeam bathymetric data and holds more than 9 million 
nautical miles of ship trackline data recorded from over 2400 cruises and received from sources 
worldwide.

https://data.ngdc.noaa.gov/platforms/

< multibeam:process=False:min_year=None:survey_id=None:exclude=None >"""

    
    def __init__(self, processed=True, process=False, min_year=None, survey_id=None, exclude=None, make_datalist=False, **kwargs):
        super().__init__(name='multibeam', **kwargs)
        self._mb_data_url = "https://data.ngdc.noaa.gov/platforms/"
        self._mb_metadata_url = "https://data.noaa.gov/waf/NOAA/NESDIS/NGDC/MGG/Multibeam/iso/"
        self._mb_search_url = "https://maps.ngdc.noaa.gov/mapviewer-support/multibeam/files.groovy?"
        self._mb_autogrid = "https://www.ngdc.noaa.gov/maps/autogrid/"
        self._mb_html = "https://www.ngdc.noaa.gov/"
        #self._outdir = os.path.join(os.getcwd(), 'mb')
        self._urls = [self._mb_data_url, self._mb_metadata_url, self._mb_autogrid]
        self.name = 'multibeam'
        self.processed_p = processed
        self.process = process
        self.min_year = utils.int_or(min_year)
        self.survey_id = survey_id
        self.exclude = exclude
        self.make_datalist = make_datalist
        self.data_format = 301
        self.src_srs = 'epsg:4326+3855'

    def mb_inf_data_format(self, src_inf):
        """extract the data format from the mbsystem inf file."""

        with open(src_inf, errors='ignore') as iob:
            for il in iob:
                til = il.split()
                if len(til) > 1:
                    if til[0] == 'MBIO':
                        return('mb'.format(til[4]))

    def mb_inf_data_date(self, src_inf):
        """extract the data format from the mbsystem inf file."""

        with open(src_inf, errors='ignore') as iob:
            for il in iob:
                til = il.split()
                if len(til) > 1:
                    if til[0] == 'Time:':
                        return(til[3])

    def mb_inf_perc_good(self, src_inf):
        """extract the data format from the mbsystem inf file."""

        with open(src_inf, errors='ignore') as iob:
            for il in iob:
                til = il.split(':')
                if len(til) > 1:
                    if til[0].strip() == 'Number of Good Beams':
                        return(til[1].split()[-1].split('%')[0])
                    
    def run(self):
        these_surveys = {}
        these_versions = {}
        if self.region is None:
            return([])
        else:
            fetch_region = self.region.copy()
            fetch_region.buffer(pct=25)
            
        _req = Fetch(self._mb_search_url).fetch_req(params={'geometry': fetch_region.format('bbox')}, timeout=20)
        if _req is not None and _req.status_code == 200:
            survey_list = _req.text.split('\n')[:-1]
            for r in survey_list:
                dst_pfn = r.split(' ')[0]
                dst_fn = dst_pfn.split('/')[-1:][0]
                survey = dst_pfn.split('/')[6]
                dn = r.split(' ')[0].split('/')[:-1]
                version = dst_pfn.split('/')[9][-1]
                data_url = self._mb_data_url + '/'.join(r.split('/')[3:])
                if self.survey_id is not None:
                    if survey != self.survey_id:
                        continue
                    
                if survey in these_surveys.keys():
                    if version in these_surveys[survey].keys():
                        these_surveys[survey][version].append([data_url.split(' ')[0], os.path.join(self._outdir, '/'.join([survey, dst_fn])), 'mb'])
                    else:
                        these_surveys[survey][version] = [[data_url.split(' ')[0], '/'.join([survey, dst_fn]), 'mb']]
                        
                else:
                    these_surveys[survey] = {version: [[data_url.split(' ')[0], '/'.join([survey, dst_fn]), 'mb']]}
                    
        else:
            utils.echo_error_msg('failed to fetch multibeam request')
                    
        for key in these_surveys.keys():
            if self.processed_p:
                if '2' in these_surveys[key].keys():
                    for v2 in these_surveys[key]['2']:
                        if self.min_year is not None:
                            try:
                                s,d,f,p,t = self.parse_entry_inf(v2)
                                if int(t) >= self.min_year:
                                    self.results.append(v2)
                            except: pass
                        else:
                            self.results.append(v2)
                else:
                    for v1 in these_surveys[key]['1']:
                        if self.min_year is not None:
                            try:
                                s,d,f,p,t = self.parse_entry_inf(v1)
                                if int(t) >= self.min_year:
                                    self.results.append(v1)
                            except:
                                pass
                        else:
                            self.results.append(v1)
                        
            else:
                for keys in these_surveys[key].keys():
                    for survs in these_surveys[key][keys]:
                        if self.min_year is not None:
                            try:
                                s,d,f,p,t = self.parse_entry_inf(survs)
                                if int(t) >= self.min_year:
                                    self.results.append(survs)
                            except: pass
                        else:
                            self.results.append(survs)

        if self.make_datalist:
            s_got = []
            with open('mb_inf.txt', 'w') as mb_inf_txt:
                for entry in self.results:
                    try:
                        survey, src_data, mb_fmt, mb_perc, mb_date = self.parse_entry_inf(entry)
                        if survey in s_got:
                            continue
                        
                        this_year = int(utils.this_year()) if self.min_year is None else self.min_year
                        this_weight = float(mb_perc) * ((int(mb_date)-2000)/(this_year-2000))/100.
                        mb_inf_txt.write('{} -1 {}\n'.format(survey, this_weight))
                        s_got.append(survey)
                        #mb_inf_txt.write('\n')
                        #self.echo_inf(entry)
                    except:
                        pass

    def echo_inf(self, entry):
        print(self.parse_entry_inf(entry))
        
    def parse_entry_inf(self, entry, keep_inf=False, out_dir=None):
        src_data = os.path.basename(entry[1])
        if src_data[-3:] == 'fbt':
            src_mb = utils.fn_basename2(src_data)
            inf_url = utils.fn_basename2(entry[0])
        else:
            inf_url = entry[0]
            src_mb = src_data
            
        survey = entry[0].split('/')[7]
        #if out_dir is not None:
            
        src_inf = os.path.join(self._outdir, '{}.inf'.format(entry[1]))
        #utils.echo_msg(src_inf)
        try:
            status = Fetch('{}.inf'.format(inf_url), callback=self.callback, verbose=True).fetch_file(src_inf)
        except:
            utils.echo_warning_msg('failed to fetch inf file: {}.inf'.format(inf_url))
            status = -1
            
        if status == 0:
            mb_fmt = self.mb_inf_data_format(src_inf)
            mb_date = self.mb_inf_data_date(src_inf)
            mb_perc = self.mb_inf_perc_good(src_inf)
            if not keep_inf:
                utils.remove_glob(src_inf)
                
            return(survey, src_data, mb_fmt, mb_perc, mb_date)

## ==============================================
## NOAA NOS
## ==============================================
class HydroNOS(FetchModule):
    """NOS Soundings (bag/hydro)
    
NCEI maintains the National Ocean Service Hydrographic Data Base (NOSHDB) and Hydrographic 
Survey Meta Data Base (HSMDB). Both are populated by the Office of Coast Survey and National 
Geodetic Service, and provide coverage of coastal waters and the U.S. exclusive economic zone 
and its territories. 

Fields:

SURVEY_ID ( type: esriFieldTypeString, alias: Survey ID, length: 10 )
DATE_SURVEY_BEGIN ( type: esriFieldTypeDate, alias: Begin Date, length: 8 )
DATE_SURVEY_END ( type: esriFieldTypeDate, alias: End Date, length: 8 )
DATE_MODIFY_DATA ( type: esriFieldTypeDate, alias: Modify Data Date, length: 8 )
DATE_SURVEY_APPROVAL ( type: esriFieldTypeDate, alias: Survey Approval Date, length: 8 )
DATE_ADDED ( type: esriFieldTypeDate, alias: Date Added, length: 8 )
SURVEY_YEAR ( type: esriFieldTypeDouble, alias: Survey Year )
DIGITAL_DATA ( type: esriFieldTypeString, alias: Digital Data?, length: 15 )
LOCALITY ( type: esriFieldTypeString, alias: Locality, length: 150 )
SUBLOCALITY ( type: esriFieldTypeString, alias: Sublocality, length: 150 )
PLATFORM ( type: esriFieldTypeString, alias: Platform Name, length: 150 )
PRODUCT_ID ( type: esriFieldTypeString, alias: Product ID, length: 24 )
BAGS_EXIST ( type: esriFieldTypeString, alias: BAGS_EXIST, length: 4 )
DOWNLOAD_URL ( type: esriFieldTypeString, alias: Download URL, length: 256 )
DECADE ( type: esriFieldTypeDouble, alias: Decade )
PUBLISH ( type: esriFieldTypeString, alias: PUBLISH, length: 1 )
OBJECTID ( type: esriFieldTypeOID, alias: OBJECTID )
SHAPE ( type: esriFieldTypeGeometry, alias: SHAPE )
    
Layer 0: Surveys with BAGs available (Bathymetric Attributed Grids).
Layer 1: Surveys with digital sounding data available for download (including those with BAGs).

https://www.ngdc.noaa.gov/mgg/bathymetry/hydro.html

< nos:where=None:layer=0:datatype=None:index=False >"""
    
    def __init__(self, where='1=1', layer=1, datatype=None, index=False, **kwargs):
        super().__init__(name='hydronos', **kwargs)
        self._nos_dynamic_url = 'https://gis.ngdc.noaa.gov/arcgis/rest/services/web_mercator/nos_hydro_dynamic/MapServer'
        self._nos_url = 'https://gis.ngdc.noaa.gov/arcgis/rest/services/web_mercator/nos_hydro/MapServer'
        self._nos_data_url = 'https://data.ngdc.noaa.gov/platforms/ocean/nos/coast/'
        self._nos_query_url = '{0}/{1}/query?'.format(self._nos_dynamic_url, layer)
        self.where = where
        self.datatype = datatype
        self.index = index
        self.data_format = None # bag/xyz data are different, reset later
        self.src_srs = None # bag/xyz data are different, reset later
        
    def run(self):
        if self.region is None:
            return([])

        _data = {'where': self.where, 'outFields': '*', 'geometry': self.region.format('bbox'),
                 'inSR':4326, 'outSR':4326, 'f':'pjson', 'returnGeometry':'False'}
        _req = Fetch(self._nos_query_url, verbose=self.verbose).fetch_req(params=_data)

        if _req is not None:
            features = _req.json()
            if 'features' in features.keys():
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
                                page = Fetch(data_link + 'BAG').fetch_html()
                                bags = page.xpath('//a[contains(@href, ".bag")]/@href')
                                #[self.results.append(['{0}BAG/{1}'.format(data_link, bag), os.path.join(self._outdir, 'bag', bag), 'bag']) for bag in bags]
                                [self.results.append(['{0}BAG/{1}'.format(data_link, bag), os.path.join('bag', bag), 'bag']) for bag in bags]

                        if self.datatype is None or 'xyz' in self.datatype.lower():
                            page = Fetch(data_link).fetch_html()
                            if page is not None:
                                geodas = page.xpath('//a[contains(@href, "GEODAS")]/@href')
                                if geodas:
                                    xyz_link = data_link + 'GEODAS/{0}.xyz.gz'.format(ID)
                                    self.results.append([xyz_link, os.path.join('geodas', xyz_link.split('/')[-1]), 'xyz'])                

## ==============================================
## NOAA DEMs
## doesn't really work well, use ncei_thredds or digital_coast instead...
## ==============================================
class DEMMosaic(FetchModule):
    """
Fields:

    OBJECTID ( type: esriFieldTypeOID, alias: OBJECTID )
    Shape ( type: esriFieldTypeGeometry, alias: Shape )
    Name ( type: esriFieldTypeString, alias: Name, length: 200 )
    MinPS ( type: esriFieldTypeDouble, alias: MinPS )
    MaxPS ( type: esriFieldTypeDouble, alias: MaxPS )
    LowPS ( type: esriFieldTypeDouble, alias: LowPS )
    HighPS ( type: esriFieldTypeDouble, alias: HighPS )
    Category ( type: esriFieldTypeInteger, alias: Category , Coded Values: [0: Unknown] , [1: Primary] , [2: Overview] , ...6 more... )
    Tag ( type: esriFieldTypeString, alias: Tag, length: 100 )
    GroupName ( type: esriFieldTypeString, alias: GroupName, length: 100 )
    ProductName ( type: esriFieldTypeString, alias: ProductName, length: 100 )
    CenterX ( type: esriFieldTypeDouble, alias: CenterX )
    CenterY ( type: esriFieldTypeDouble, alias: CenterY )
    ZOrder ( type: esriFieldTypeInteger, alias: ZOrder )
    Shape_Length ( type: esriFieldTypeDouble, alias: Shape_Length )
    Shape_Area ( type: esriFieldTypeDouble, alias: Shape_Area )
    ZOrder_1 ( type: esriFieldTypeInteger, alias: ZOrder_1 )
    DEM_ID ( type: esriFieldTypeSmallInteger, alias: DEM_ID )
    DateCompleted ( type: esriFieldTypeDate, alias: DateCompleted, length: 8 )
    CellsizeArcseconds ( type: esriFieldTypeSingle, alias: CellsizeArcseconds )
    DemName ( type: esriFieldTypeString, alias: DemName, length: 100 )
    MetadataURL ( type: esriFieldTypeString, alias: MetadataURL, length: 250 )
    VerticalDatum ( type: esriFieldTypeString, alias: VerticalDatum, length: 50 )

    https://gis.ngdc.noaa.gov/arcgis/rest/services/DEM_mosaics/DEM_global_mosaic/ImageServer

    """
    
    def __init__(self, where='1=1', layer=1, index=False, **kwargs):
        super().__init__(name='hydronos', **kwargs)
        self._dem_mosaic_url = 'https://gis.ngdc.noaa.gov/arcgis/rest/services/DEM_mosaics/DEM_all/ImageServer'
        #self._dem_mosaic_query_url = '{0}/{1}/query?'.format(self._dem_mosaic_url, layer)
        self._dem_mosaic_query_url = '{0}/query?'.format(self._dem_mosaic_url)
        self.where = where
        self.index = index
        self.data_format = None # bag/xyz data are different, reset later
        self.src_srs = None # dems vary, set later
        
    def run(self):
        if self.region is None:
            return([])

        _data = {'where': self.where, 'outFields': '*', 'geometry': self.region.format('bbox'),
                 'inSR':4326, 'outSR':4326, 'f':'pjson', 'returnGeometry':'False',
                 'geometryType':'esriGeometryEnvelope', 'spatialRel':'esriSpatialRelIntersects'}
        _req = Fetch(self._dem_mosaic_query_url, verbose=self.verbose).fetch_req(params=_data)

        if _req is not None:
            #utils.echo_msg(_req.url)
            features = _req.json()
            for feature in features['features']:
                if self.index:
                    print(json.dumps(feature['attributes'], indent=4))
                else:
                    #utils.echo_msg(feature['attributes']['MetadataURL'])
                    print(feature)
                    Name = feature['attributes']['Name']
                    ID = feature['attributes']['DEM_ID']
                    link = feature['attributes']['MetadataURL']

                    if link is not None:
                        utils.echo_msg(Name)
                        utils.echo_msg(ID)
                        utils.echo_msg(link)

                        page = Fetch(link).fetch_xml()
                        print(page)
                        sys.exit()
                        
                    # nos_dir = link.split('/')[-2]
                    # data_link = '{}{}/{}/'.format(self._nos_data_url, nos_dir, ID)

                    # if self.datatype is None or 'bag' in self.datatype.lower():
                    #     if feature['attributes']['BAGS_EXIST'] == 'TRUE':
                    #         page = Fetch(data_link + 'BAG').fetch_html()
                    #         bags = page.xpath('//a[contains(@href, ".bag")]/@href')
                    #         #[self.results.append(['{0}BAG/{1}'.format(data_link, bag), os.path.join(self._outdir, 'bag', bag), 'bag']) for bag in bags]
                    #         [self.results.append(['{0}BAG/{1}'.format(data_link, bag), os.path.join('bag', bag), 'bag']) for bag in bags]

                    # if self.datatype is None or 'xyz' in self.datatype.lower():
                    #     page = Fetch(data_link).fetch_html()
                    #     if page is not None:
                    #         geodas = page.xpath('//a[contains(@href, "GEODAS")]/@href')
                    #         if geodas:
                    #             xyz_link = data_link + 'GEODAS/{0}.xyz.gz'.format(ID)
                    #             self.results.append([xyz_link, os.path.join('geodas', xyz_link.split('/')[-1]), 'xyz'])                
                                
## ==============================================
## NOAA Trackline
## ==============================================
class Trackline(FetchModule):
    """NOAA TRACKLINE bathymetric data.

http://www.ngdc.noaa.gov/trackline/

< trackline >"""
    
    def __init__(self, where='1=1', **kwargs):
        super().__init__(name='trackline', **kwargs)
        #self._trackline_url = 'https://gis.ngdc.noaa.gov/arcgis/rest/services/web_mercator/trackline_bathymetry/MapServer/0'
        self._trackline_url = 'https://gis.ngdc.noaa.gov/arcgis/rest/services/web_mercator/trackline_combined_dynamic/MapServer/1'
        #self._trackline_url = 'https://services2.arcgis.com/C8EMgrsFcRFL6LrL/arcgis/rest/services/trackline_geophysics/FeatureServer/8'
        self._trackline_query_url = '{0}/query?'.format(self._trackline_url)
        self.where = where
        
    def run(self):
        if self.region is None:
            return([])

        _data = {'where': self.where, 'outFields': '*', 'geometry': self.region.format('bbox'),
                 'inSR':4326, 'outSR':4326, 'f':'pjson', 'returnGeometry':'False'}
                 #'geometryType':'esriGeometryEnvelope', 'spatialRel':'esriSpatialRelIntersects'}
        _req = Fetch(self._trackline_query_url, verbose=self.verbose).fetch_req(params=_data)
        if _req is not None:
            features = _req.json()
            ids = []
            for feature in features['features']:
                ids.append(feature['attributes']['SURVEY_ID'])
                #data_link = feature['attributes']['DOWNLOAD_URL']
                #utils.echo_msg(data_link)

            print('http://www.ngdc.noaa.gov/trackline/request/?surveyIds={}'.format(','.join(ids)))
            #xyz_link = 'http://www.ngdc.noaa.gov/trackline/request/?surveyIds={}'.format(','.join(ids))
            #self.results.append([xyz_link, os.path.join(self._outdir, xyz_link.split('/')[-1]), 'trackline'])
            #self.results.append([xyz_link, os.path.join(self._outdir, 'tmp_trackline.xyz'), 'trackline'])

## ==============================================
## eHydro (USACE)
## ==============================================
class eHydro(FetchModule):
    """USACE eHydro bathymetric data.
    
Maintenance responsibility for more than 25,000 miles of navigation channels and 400 ports and 
harbors throughout the United States requires extensive surveying and mapping services, including 
boundary, topographic, hydrographic, terrestrial lidar, and multispectral and hyperspectral aerial 
imagery collection as well as airborne topographic and bathymetric lidar acquisition, project-level 
GIS implementation, development of file-based geodatabases, and GIS tool development.

Three representative survey and mapping datasets include the National Channel Framework (NCF)—an enterprise 
geodatabase of information on all 61 USACE-maintained high-tonnage channels —hydrographic surveys, which 
provide assistance in locating navigable channels, determining dredging requirements, verifying dredging 
accuracy, and maintaining harbors and rivers —and Inland Electronic Navigational Charts(IENC), accurate 
navigational charts provided in a highly structured data format for use in navigation systems and to increase 
overall navigational safety.. 

https://navigation.usace.army.mil/Survey/Hydro

Fields:

    objectid (type: esriFieldTypeOID, alias: objectid, SQL Type: sqlTypeOther, length: 0, nullable: false, editable: false)
    surveyjobidpk (type: esriFieldTypeString, alias: SURVEYJOBIDPK, SQL Type: sqlTypeOther, length: 100, nullable: true, editable: true)
    sdsid (type: esriFieldTypeString, alias: SDSID, SQL Type: sqlTypeOther, length: 40, nullable: true, editable: true)
    sdsfeaturename (type: esriFieldTypeString, alias: SDSFEATURENAME, SQL Type: sqlTypeOther, length: 80, nullable: true, editable: true)
    sdsmetadataid (type: esriFieldTypeString, alias: SDSMETADATAID, SQL Type: sqlTypeOther, length: 80, nullable: true, editable: true)
    surveytype (type: esriFieldTypeString, alias: SURVEYTYPE, SQL Type: sqlTypeOther, length: 26, nullable: true, editable: true)
    channelareaidfk (type: esriFieldTypeString, alias: CHANNELAREAIDFK, SQL Type: sqlTypeOther, length: 50, nullable: true, editable: true)
    dateuploaded (type: esriFieldTypeDate, alias: dateUploaded, SQL Type: sqlTypeOther, length: 8, nullable: true, editable: true)
    usacedistrictcode (type: esriFieldTypeString, alias: usaceDistrictCode, SQL Type: sqlTypeOther, length: 20, nullable: true, editable: true)
    surveydatestart (type: esriFieldTypeDate, alias: SURVEYDATESTART, SQL Type: sqlTypeOther, length: 8, nullable: true, editable: true)
    surveydateend (type: esriFieldTypeDate, alias: SURVEYDATEEND, SQL Type: sqlTypeOther, length: 8, nullable: true, editable: true)
    sourcedatalocation (type: esriFieldTypeString, alias: SOURCEDATALOCATION, SQL Type: sqlTypeOther, length: 255, nullable: true, editable: true)
    sourceprojection (type: esriFieldTypeString, alias: SOURCEPROJECTION, SQL Type: sqlTypeOther, length: 75, nullable: true, editable: true)
    mediaidfk (type: esriFieldTypeString, alias: MEDIAIDFK, SQL Type: sqlTypeOther, length: 100, nullable: true, editable: true)
    projectedarea (type: esriFieldTypeDouble, alias: PROJECTEDAREA, SQL Type: sqlTypeOther, nullable: true, editable: true)
    sdsfeaturedescription (type: esriFieldTypeString, alias: SDSFEATUREDESCRIPTION, SQL Type: sqlTypeOther, length: 255, nullable: true, editable: true)
    dateloadedenterprise (type: esriFieldTypeDate, alias: dateLoadedEnterprise, SQL Type: sqlTypeOther, length: 8, nullable: true, editable: true)
    datenotified (type: esriFieldTypeDate, alias: dateNotified, SQL Type: sqlTypeOther, length: 8, nullable: true, editable: true)
    sourcedatacontent (type: esriFieldTypeString, alias: sourceDataContent, SQL Type: sqlTypeOther, length: 1000, nullable: true, editable: true)
    plotsheetlocation (type: esriFieldTypeString, alias: PLOTSHEETLOCATION, SQL Type: sqlTypeOther, length: 255, nullable: true, editable: true)
    sourceagency (type: esriFieldTypeString, alias: SOURCEAGENCY, SQL Type: sqlTypeOther, length: 20, nullable: true, editable: true)
    globalid (type: esriFieldTypeGlobalID, alias: GlobalID, SQL Type: sqlTypeOther, length: 38, nullable: false, editable: false)
    Shape__Area (type: esriFieldTypeDouble, alias: Shape__Area, SQL Type: sqlTypeDouble, nullable: true, editable: false)
    Shape__Length (type: esriFieldTypeDouble, alias: Shape__Length, SQL Type: sqlTypeDouble, nullable: true, editable: false)
        
< ehydro:where=None:inc=None >"""

    def __init__(self, where='1=1', inc=None, survey_name=None, **kwargs):
        super().__init__(name='ehydro', **kwargs)
        self._ehydro_gj_api_url = 'https://opendata.arcgis.com/datasets/80a394bae6b547f1b5788074261e11f1_0.geojson'
        self._ehydro_api_url = 'https://services7.arcgis.com/n1YM8pTrFmm7L4hs/arcgis/rest/services/eHydro_Survey_Data/FeatureServer/0'
        self._ehydro_query_url = '{0}/query?'.format(self._ehydro_api_url)
        self.where = where
        self.survey_name = survey_name
        self.inc = utils.str2inc(inc)

    def run(self):
        '''Run the eHydro fetching module'''
        
        if self.region is None:
            return([])
        
        _data = {'where': self.where, 'outFields': '*', 'geometry': self.region.format('bbox'),
                 'inSR':4326, 'outSR':4326, 'f':'pjson', 'returnGeometry':'False'}
        _req = Fetch(self._ehydro_query_url).fetch_req(params=_data)
        if _req is not None:
            features = _req.json()
            if 'features' in features.keys():
                for feature in features['features']:
                    fetch_fn = feature['attributes']['sourcedatalocation']
                    if self.survey_name is not None:
                        if self.survey_name in fetch_fn:
                            self.results.append([fetch_fn, fetch_fn.split('/')[-1], 'ehydro'])
                    else:
                        self.results.append([fetch_fn, fetch_fn.split('/')[-1], 'ehydro'])
                
        return(self)

## ==============================================
## BlueTopo
## ==============================================
class BlueTopo(FetchModule):
    """BlueTOPO DEM
    
BlueTopo is a compilation of the nation's best available bathymetric data. 
In the same way that topographic map details the height of land, BlueTopo details the depth of 
lake beds and seafloor beneath navigationally significant U.S. waters. Created as part of the 
Office of Coast Survey nautical charting mission and its National Bathymetric Source project, 
BlueTopo is curated bathymetric source data to provide a definitive nationwide model of the seafloor 
and the Great Lakes.

Output 'tiff' files are 3 bands
1 - Elevation
2 - Uncertainty
3 - Data Source Table

yield_xyz outputs elevation (band 1)
elevation data is in NAVD88


https://nauticalcharts.noaa.gov/data/bluetopo_specs.html
https://noaa-ocs-nationalbathymetry-pds.s3.amazonaws.com/index.html#
https://www.nauticalcharts.noaa.gov/data/bluetopo.html
https://noaa-ocs-nationalbathymetry-pds.s3.amazonaws.com/index.html#BlueTopo/
    
https://www.nauticalcharts.noaa.gov/data/bluetopo.html

< bluetopo:want_interpolation=False:unc_weights=False:keep_index=False >"""
    
    def __init__(self, want_interpolation=False, unc_weights=False, keep_index=False, **kwargs):
        super().__init__(name='bluetopo', **kwargs)
        self.unc_weights = unc_weights
        self.want_interpolation = want_interpolation
        self.keep_index = keep_index
        self._bt_bucket = 'noaa-ocs-nationalbathymetry-pds'
        s3 = boto3.client('s3', aws_access_key_id='', aws_secret_access_key='')
        s3._request_signer.sign = (lambda *args, **kwargs: None)
        r = s3.list_objects(Bucket = self._bt_bucket, Prefix='BlueTopo/_BlueTopo_Tile_Scheme')
        self._bluetopo_index_url = 'https://{}.s3.amazonaws.com/{}'.format(self._bt_bucket, r['Contents'][0]['Key'])
        self._bluetopo_index = self._bluetopo_index_url.split('/')[-1]
        
    def run(self):
        s3 = boto3.client('s3', aws_access_key_id='', aws_secret_access_key='')
        s3._request_signer.sign = (lambda *args, **kwargs: None)
        try:
            status = Fetch(self._bluetopo_index_url, verbose=self.verbose).fetch_file(self._bluetopo_index)
            v_ds = ogr.Open(self._bluetopo_index)
        except:
            v_ds = None
            status = -1
            
        if v_ds is not None:
            layer = v_ds.GetLayer()
            _boundsGeom = self.region.export_as_geom()
            layer.SetSpatialFilter(_boundsGeom)            
            fcount = layer.GetFeatureCount()
            for feature in layer:
                if feature is None:
                    continue
                
                tile_name = feature.GetField('tile')
                r = s3.list_objects(Bucket = 'noaa-ocs-nationalbathymetry-pds', Prefix='BlueTopo/{}'.format(tile_name))
                if 'Contents' in r:
                    for key in r['Contents']:
                        if key['Key'].split('.')[-1] == 'tiff':
                            data_link = 'https://noaa-ocs-nationalbathymetry-pds.s3.amazonaws.com/{}'.format(key['Key'])
                            self.results.append([data_link, data_link.split('/')[-1], 'raster'])
            v_ds = None

        if not self.keep_index:
            utils.remove_glob(self._bluetopo_index)
            
        return(self)
    
## ==============================================
## MGDS
## ==============================================
class MGDS(FetchModule):
    """The Marine Geoscience Data System (MGDS)

Fetch marine data from MGDS
    
MGDS is a trusted data repository that provides free public access to a curated collection of marine geophysical 
data products and complementary data related to understanding the formation and evolution 
of the seafloor and sub-seafloor.

https://www.marine-geo.org

data_tpye=[Bathymetry, Bathymetry:Phase, Bathymetry:Swath, Bathymetry:Swath:Ancillary, Bathymetry:Singlebeam, Bathymetry:BPI, Bathymetry:ReferenceSurface, Bathymetry:Paelobathymetry]
            
< mgds:data_type=Bathymetry >"""
    
    def __init__(self, data_type='Bathymetry', **kwargs):
        super().__init__(name='mgds', **kwargs) 
        self._mgds_file_url = "https://www.marine-geo.org/services/FileServer?"
        self._mgds_filedownload_url = "http://www.marine-geo.org/services/FileDownloadServer?"
        self._mgds_filemetadata_url = "http://www.marine-geo.org/services/FileDownloadServer/metadata?"
        self._mgds_archive_url = "http://www.marine-geo.org/services/FileDownloadServer/metadata?"
        self._mgds_search_url = "http://www.marine-geo.org/services/search/datasets??"
        self.data_type = data_type.replace(',', ':')
        
    def run(self):
        '''Run the MGDS fetching module'''

        if self.region is None:
            return([])

        self.data = {'north':self.region.ymax, 'west':self.region.xmin,
                     'south':self.region.ymin, 'east':self.region.xmax,
                     'format':'summary', 'data_type':'{}'.format(self.data_type)}        
        req = Fetch(self._mgds_file_url).fetch_req(
            params=self.data, tries=10, timeout=2
        )
        if req is not None:
            req_xml = lxml.etree.fromstring(req.content)
            req_results = req_xml.findall('.//{https://www.marine-geo.org/services/xml/mgdsDataService}file')
            for req_result in req_results:
                name = req_result.attrib['name']
                link = req_result.attrib['download']
                self.results.append([link, name, 'mgds'])
                
        return(self)

## ==============================================
## NGS - geodesic monuments
## ==============================================
class NGS(FetchModule):
    """NGS Monuments
    
NGS provides Information about survey marks (including bench marks) in text datasheets or in GIS shapefiles. 
Note some survey markers installed by other organizations may not be available through NGS.

Fetch NGS monuments from NOAA
    
http://geodesy.noaa.gov/

< ngs:datum=geoidHt >"""

    
    def __init__(self, datum='geoidHt', **kwargs):
        super().__init__(name='ngs', **kwargs)
        self._ngs_search_url = 'http://geodesy.noaa.gov/api/nde/bounds?'
        if datum not in ['orthoHt', 'geoidHt', 'z', 'ellipHeight']:
            utils.echo_warning_msg('could not parse {}, falling back to geoidHt'.format(datum))
            self.datum = 'geoidHt'
        else:
            self.datum = datum

    def run(self, csv=False):
        """Run the NGS (monuments) fetching module."""
        
        if self.region is None:
            return([])
        
        _data = {'maxlon':self.region.xmax, 'minlon':self.region.xmin,
                 'maxlat':self.region.ymax, 'minlat':self.region.ymin}

        _req = Fetch(self._ngs_search_url).fetch_req(params=_data)
        if _req is not None:
            self.results.append([_req.url, 'ngs_results_{}.json'.format(self.region.format('fn')), 'ngs'])
            
        return(self)

## ==============================================
## TIDES
## ==============================================
class Tides(FetchModule):
    """TIDE station information from NOAA/NOS

Fetch NOS Tide Stations

Fields:
    objectid ( type: esriFieldTypeOID , alias: objectid , editable: false , nullable: false )
    id ( type: esriFieldTypeString , alias: id , editable: true , nullable: true , length: 50 )
    name ( type: esriFieldTypeString , alias: name , editable: true , nullable: true , length: 50 )
    affil ( type: esriFieldTypeString , alias: affil , editable: true , nullable: true , length: 50 )
    latitude ( type: esriFieldTypeDouble , alias: latitude , editable: true , nullable: true )
    longitude ( type: esriFieldTypeDouble , alias: longitude , editable: true , nullable: true )
    data ( type: esriFieldTypeString , alias: data , editable: true , nullable: true , length: 200 )
    dataapi ( type: esriFieldTypeString , alias: dataapi , editable: true , nullable: true , length: 200 )
    accepted ( type: esriFieldTypeString , alias: accepted , editable: true , nullable: true , length: 50 )
    epoch ( type: esriFieldTypeString , alias: epoch , editable: true , nullable: true , length: 50 )
    units ( type: esriFieldTypeString , alias: units , editable: true , nullable: true , length: 50 )
    orthodatum ( type: esriFieldTypeString , alias: orthodatum , editable: true , nullable: true , length: 50 )
    mhhw ( type: esriFieldTypeDouble , alias: mhhw , editable: true , nullable: true )
    mhw ( type: esriFieldTypeDouble , alias: mhw , editable: true , nullable: true )
    mtl ( type: esriFieldTypeDouble , alias: mtl , editable: true , nullable: true )
    msl ( type: esriFieldTypeDouble , alias: msl , editable: true , nullable: true )
    dtl ( type: esriFieldTypeDouble , alias: dtl , editable: true , nullable: true )
    mlw ( type: esriFieldTypeDouble , alias: mlw , editable: true , nullable: true )
    mllw ( type: esriFieldTypeDouble , alias: mllw , editable: true , nullable: true )
    stnd ( type: esriFieldTypeDouble , alias: stnd , editable: true , nullable: true )
    mn ( type: esriFieldTypeDouble , alias: mn , editable: true , nullable: true )
    dhq ( type: esriFieldTypeDouble , alias: dhq , editable: true , nullable: true )
    dlq ( type: esriFieldTypeDouble , alias: dlq , editable: true , nullable: true )
    hwi ( type: esriFieldTypeDouble , alias: hwi , editable: true , nullable: true )
    lwi ( type: esriFieldTypeDouble , alias: lwi , editable: true , nullable: true )
    gt ( type: esriFieldTypeDouble , alias: gt , editable: true , nullable: true )
    navd88 ( type: esriFieldTypeDouble , alias: navd88 , editable: true , nullable: true )
    wl_max ( type: esriFieldTypeDouble , alias: wl_max , editable: true , nullable: true )
    max_date ( type: esriFieldTypeString , alias: max_date , editable: true , nullable: true , length: 50 )
    wl_min ( type: esriFieldTypeDouble , alias: wl_min , editable: true , nullable: true )
    min_date ( type: esriFieldTypeString , alias: min_date , editable: true , nullable: true , length: 50 )
    
https://tidesandcurrents.noaa.gov/

< tides:station_id=None:s_datum=mllw:t_datum=msl:units=m >"""

    
    def __init__(self, s_datum='mllw', t_datum='msl', units='m', station_id=None, **kwargs):
        super().__init__(name='tides', **kwargs)
        #self._stations_api_url = 'https://idpgis.ncep.noaa.gov/arcgis/rest/services/NOS_Observations/CO_OPS_Products/FeatureServer/0/query?'
        self._stations_api_url = 'https://mapservices.weather.noaa.gov/static/rest/services/NOS_Observations/CO_OPS_Products/FeatureServer/0/query?'
        #self._stations_api_url = 'https://api.tidesandcurrents.noaa.gov/api/prod/datagetter?'
        self.s_datum = s_datum
        self.t_datum = t_datum
        self.units = units
        self.station_id = station_id
        
    def run(self):
        '''Run the TIDES fetching module'''
        
        if self.region is None:
            return([])
        
        _data = {
            'outFields': '*',
            'units': 'esriSRUnit_Meter',
            'geometry': self.region.format('bbox'),
            'inSR':4326,
            'outSR':4326,
            'f':'pjson',
        }
        _req = Fetch(self._stations_api_url, verbose=self.verbose).fetch_req(params=_data)
        if _req is not None:
            self.results.append([_req.url, 'tides_results_{}.json'.format(self.region.format('fn')), 'tides'])
            
        return(self)

## ==============================================
## WaterServices (USGS)
## ==============================================
class WaterServices(FetchModule):
    """WaterServices station information from USGS

https://waterservices.usgs.gov/

< waterservices >"""

    
    def __init__(self, **kwargs):
        super().__init__(name='waterservices', **kwargs)
        self._water_services_api_url = 'https://waterservices.usgs.gov/nwis/iv/?'
        
    def run(self):
        '''Run the WATERSERVICES fetching module'''
        
        if self.region is None:
            return([])
        
        _data = {
            'bBox': self.region.format('bbox'),
            'siteStatus': 'active',
            'format':'json',
        }
        _req = Fetch(self._water_services_api_url, verbose=self.verbose).fetch_req(params=_data)
        if _req is not None:
            self.results.append([_req.url, 'water_services_results_{}.json'.format(self.region.format('fn')), 'waterservices'])
            
        return(self)
    
## ==============================================
## buoys 
## ==============================================
class BUOYS(FetchModule):
    """NOAA BUOY data (beta)

Fetch NOS Tide Stations

A sustainable and resilient marine observation and monitoring infrastructure which enhances healthy 
ecosystems, communities, and economies in the face of change and To provide quality observations in 
the marine environment in a safe and sustainable manner to support the understanding of and predictions 
to changes in weather, climate, oceans and coast. 

https://www.ndbc.noaa.gov

< buoys:buoy_id=None >"""
    
    def __init__(self, buoy_id=None, **kwargs):
        super().__init__(name='buoys', **kwargs)
        self._ndbc_url = 'https://www.ndbc.noaa.gov'
        self._buoy_box_search_url = 'https://www.ndbc.noaa.gov/box_search.php?'
        self._buoy_station_url = 'https://www.ndbc.noaa.gov/station_page.php?'
        self._buoy_stations_url = 'https://www.ndbc.noaa.gov/to_station.shtml'
        self._buoy_station_kml = 'https://www.ndbc.noaa.gov/kml/marineobs_by_owner.kml'
        self._buoy_station_realtime = 'https://www.ndbc.noaa.gov/data/realtime2/'
        self.buoy_id = buoy_id

    def run(self):
        '''Run the BOUYS fetching module'''
        
        if self.region is None:
            return([])

        _data = {
            'lat1': self.region.ymin,
            'lat2': self.region.ymax,
            'lon1': self.region.xmin,
            'lon2': self.region.xmax,
            'uom': 'M',
            'ot': 'A',
            'time': 0,
        }

        ## Fetch buoy ids from box search
        _req = Fetch(
            self._buoy_box_search_url, verbose=self.verbose
        ).fetch_req(params=_data)
        if _req is not None:
            #print(_req.content)
            doc = lh.document_fromstring(_req.text)
            sp = doc.xpath('//span')
            current_stations = []
            
            for s in sp:
                #print(s.text_content())
                if len(s.xpath('a')) > 0:
                    station_url = s.xpath('a')[0].get('href')
                    if 'station=' in station_url:
                        station_id = station_url.split('=')[-1]
                        #self.results.append([self._ndbc_url + station_url, 'buoy_results_{}.html'.format(station_id), 'buoys'])
                        if station_id not in current_stations:
                            current_stations.append(station_id)
            for station_id in current_stations:
                self.results.append([self._buoy_station_realtime + station_id + '.txt', 'buoy_results_{}.txt'.format(station_id), 'buoys'])
            
        return(self)

## ==============================================
## Digital Coast - Data Access Viewer
## ==============================================
class DAV(FetchModule):
    """Fetch NOAA lidar data from DAV

Uses Digital Coasts Data Access Viewer Mapserver to discover
dataset footprints.

This map service presents spatial information about Elevation Data Access Viewer services across the United States
and Territories in the Web Mercator projection. The service was developed by the National Oceanic and Atmospheric
Administration (NOAA), but may contain data and information from a variety of data sources, including non-NOAA data.
NOAA provides the information “as-is” and shall incur no responsibility or liability as to the completeness or accuracy
of this information. NOAA assumes no responsibility arising from the use of this information. The NOAA Office for Coastal
Management will make every effort to provide continual access to this service but it may need to be taken down during
routine IT maintenance or in case of an emergency. If you plan to ingest this service into your own application and would
like to be informed about planned and unplanned service outages or changes to existing services, please register for our
Data Services Newsletter (http://coast.noaa.gov/digitalcoast/publications/subscribe). For additional information, please
contact the NOAA Office for Coastal Management (coastal.info@noaa.gov).


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
    
https://coast.noaa.gov
    
Use where=SQL_QUERY to query the MapServer to filter datasets

* For CUDEM tiles, use where="ID=8483" (1/9) or where="ID=8580" (1/3) or where="Name LIKE '%CUDEM%'" for all available.
* For OCM SLR DEMs, use where="ID=6230" or where="Name LIKE '%Sea Level Rise%'"
* For USGS CoNED DEMs, use where="ID=9181" or where="Name LIKE '%CoNED%'"
* To only return lidar data, use datatype=lidar, for only raster, use datatype=dem
* datatype is either 'lidar', 'dem' or 'sm'

< digital_coast:where=None:datatype=None >"""
    
    def __init__(self, where='1=1', index=False, datatype=None, **kwargs):
        super().__init__(name='digital_coast', **kwargs)
        self._dav_api_url = 'https://maps.coast.noaa.gov/arcgis/rest/services/DAV/ElevationFootprints/MapServer/0/query?'
        self.where = where
        self.index = index
        self.datatype = datatype

        ## data formats vary
        self.data_format = None
        
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
        _req = Fetch(self._dav_api_url, verbose=self.verbose).fetch_req(params=_data)
        utils.echo_msg(_req.url)
        if _req is not None:
            features = _req.json()
            if not 'features' in features.keys():
                utils.echo_error_msg('DAV failed to execute query...try again later')
            else:
                for feature in features['features']:
                    if self.datatype is not None:
                        if self.datatype.lower() != feature['attributes']['DataType'].lower():
                            if self.datatype.lower() != 'sm':
                                continue

                    self.vdatum = vdatums.get_vdatum_by_name(feature['attributes']['NativeVdatum'])                        
                    links = json.loads(feature['attributes']['ExternalProviderLink'])
                    ept_infos = None
                    ## get ept link to gather datum infos...for lidar only apparently...
                    for link in links['links']:
                        if link['serviceID'] == 167:
                            if link['label'] == 'EPT NOAA':
                                ept_req = Fetch(link['link'], verbose=True).fetch_req()
                                ept_infos = ept_req.json()

                    ## get metadata for datum infos...for raster data
                    if self.index:
                        feature['attributes']['ExternalProviderLink'] = links
                        utils.echo_msg(json.dumps(feature['attributes'], indent=4))
                    else:
                        for link in links['links']:
                            if link['serviceID'] == 46 and (self.datatype == 'lidar' or self.datatype == 'dem' or self.datatype is None):
                                urllist = 'urllist' + str(feature['attributes']['ID']) + '.txt'
                                surv_name = '_'.join(link['link'].split('/')[-1].split('_')[:-1])
                                #index_zipfile = 'tileindex.zip'
                                index_zipfile = 'tileindex_{}.zip'.format(surv_name)
                                index_zipurl = link['link'] + '/' + index_zipfile
                                #urllist_url = '/'.join(link['link'].split('/')[:-1]) + '/' + urllist
                                urllist_url = link['link'] + '/' + urllist
                                while True:
                                    try:
                                        status = Fetch(urllist_url, verbose=True).fetch_file(urllist)
                                    except:
                                        status = -1
                                        
                                    if status != 0:
                                        #if Fetch(urllist_url, verbose=True).fetch_file(urllist) != 0:
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
                                try:
                                    status = Fetch(
                                        index_zipurl, callback=self.callback, verbose=self.verbose
                                    ).fetch_file(index_zipfile)
                                except:
                                    status = -1
                                    
                                if status == 0:
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
                                            tile_name = None
                                            try:
                                                tile_name = index_feature.GetField('Name').strip()
                                            except:
                                                tile_name = index_feature.GetField('location').strip()
                                                
                                            tile_url = index_feature.GetField('URL').strip()
                                            tile_url = '/'.join(tile_url.split('/')[:-1]) + '/' + tile_name.split('/')[-1]
                                            ## add vertical datum to output;
                                            ## field is NativeVdatum
                                            ## must get from metadata
                                            if ept_infos is None:
                                                this_epsg = vdatums.get_vdatum_by_name(feature['attributes']['NativeVdatum'])
                                            else:
                                                #print(ept_infos['srs'])
                                                ## horizontal datum is wrong in ept, most seem to be nad83
                                                #this_epsg = 'epsg:{}+{}'.format(ept_infos['srs']['horizontal'], ept_infos['srs']['vertical'])
                                                if 'vertical' in ept_infos['srs'].keys():
                                                    vertical_epsg = ept_infos['srs']['vertical']
                                                    horizontal_epsg = ept_infos['srs']['horizontal']
                                                    this_epsg = 'epsg:{}+{}'.format(horizontal_epsg, vertical_epsg)
                                                    #this_epsg = 'epsg:4269+{}'.format(ept_infos['srs']['vertical'])
                                                else:
                                                    # try to extract the vertical datum from the wkt
                                                    #horizontal_epsg = ept_infos['srs']['horizontal']
                                                    this_wkt = ept_infos['srs']['wkt']
                                                    dst_horz, dst_vert = gdalfun.epsg_from_input(this_wkt)
                                                    this_epsg = '{}+{}'.format(dst_horz, dst_vert)

                                            self.results.append(
                                                [tile_url,
                                                 os.path.join('{}/{}'.format(feature['attributes']['ID'], tile_url.split('/')[-1])),
                                                 this_epsg,
                                                 feature['attributes']['DataType']]
                                            )

                                    index_ds = index_layer = None
                                    utils.remove_glob(index_zipfile, *index_shps)
                            elif link['serviceID'] == 166 and self.datatype == 'sm': # spatial_metadata
                                self.results.append(
                                    [link['link'],
                                     os.path.join('{}/{}'.format(feature['attributes']['ID'], link['link'].split('/')[-1])),
                                     None,
                                     link['label']
                                    ]
                                )

        return(self)

## ==============================================
## Digital Coast - Data Access Viewer - SLR shortcut
## ==============================================
class SLR(DAV):
    def __init__(self, **kwargs):
        super().__init__(where='ID=6230', **kwargs)

## ==============================================
## Digital Coast - Data Access Viewer CoNED shortcut
## ==============================================
class CoNED(DAV):
    def __init__(self, **kwargs):
        super().__init__(where="NAME LIKE '%CoNED%'", **kwargs)

## ==============================================
## Digital Coast - Data Access Viewer - CUDEM shortcut
## ==============================================
class CUDEM(DAV):
    def __init__(self, **kwargs):
        super().__init__(where="NAME LIKE '%CUDEM%'", **kwargs)
    
## ==============================================
## NCEI THREDDS Catalog
## ==============================================
class NCEIThreddsCatalog(FetchModule):
    """NOAA NCEI DEMs via THREDDS

Fetch DEMs from NCEI THREDDS Catalog
    
Digital Elevation Models around the world at various resolutions and extents.
NCEI builds and distributes high-resolution, coastal digital elevation models (DEMs) that integrate ocean 
bathymetry and land topography supporting NOAA's mission to understand and predict changes in Earth's environment, 
and conserve and manage coastal and marine resources to meet our Nation's economic, social, and environmental needs.

DEMs are used for coastal process modeling (tsunami inundation, storm surge, sea-level rise, contaminant dispersal, 
etc.), ecosystems management and habitat research, coastal and marine spatial planning, and hazard mitigation and 
community preparedness.

https://www.ngdc.noaa.gov/thredds/demCatalog.xml

< ncei_thredds:where=None:want_wcs=False >"""

    def __init__(self, where=[], want_wcs=False, datatype=None, **kwargs):
        super().__init__(name='ncei_thredds', **kwargs)
        self._nt_catalog = "https://www.ngdc.noaa.gov/thredds/catalog/demCatalog.xml"
        self._ngdc_url = "https://www.ngdc.noaa.gov"
        self.where = [where] if len(where) > 0 else []
        self.want_wcs = want_wcs
        self.datatype = datatype
        if self.datatype is not None:
            self.where.append("ID LIKE '%{}%'".format(datatype))
            
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
        ntCatalog = iso_xml(catalog_url)
        ntCatRefs = ntCatalog.xml_doc.findall('.//th:catalogRef', namespaces = ntCatalog.namespaces)
        for ntCatRef in ntCatRefs:
            ntCatHref =ntCatRef.attrib['{http://www.w3.org/1999/xlink}href']
            if ntCatHref[0] == "/":
                ntCatUrl = '{}{}'.format(self._ngdc_url, ntCatHref)
            else:
                ntCatUrl = '{}/{}'.format(os.path.dirname(catalog_url), ntCatHref)
                
            self._parse_dataset(ntCatUrl)
            
    def _parse_dataset(self, catalog_url):
        ntCatXml = iso_xml(catalog_url)
        this_ds = ntCatXml.xml_doc.findall('.//th:dataset', namespaces = ntCatXml.namespaces)
        this_ds_services = ntCatXml.xml_doc.findall('.//th:service', namespaces = ntCatXml.namespaces)
        surveys = []
        with tqdm(
                total=len(this_ds),
                desc='scanning NCEI THREDDS datasets in {}'.format(this_ds[0].attrib['name']),
                leave=self.verbose
        ) as pbar:
            for i, node in enumerate(this_ds):
                this_title = node.attrib['name']
                try:
                    this_id = node.attrib['ID']
                except:
                    this_id = None

                pbar.update(1)
                self.FRED._attribute_filter(["ID = '{}'".format(this_id)])
                if self.FRED.layer is None or len(self.FRED.layer) == 0:
                    subCatRefs = node.findall('.//th:catalogRef', namespaces=ntCatXml.namespaces)
                    if len(subCatRefs) > 0:
                        self._parse_catalog(catalog_url)
                        break

                    try:
                        ds_path = node.attrib['urlPath']
                    except:
                        continue

                    iso_url = False
                    wcs_url = False
                    http_url = False
                    for service in this_ds_services:
                        service_name = service.attrib['name']
                        if service_name == 'iso': iso_url = '{}{}{}'.format(self._ngdc_url, service.attrib['base'], ds_path)
                        if service_name == 'wcs': wcs_url = '{}{}{}'.format(self._ngdc_url, service.attrib['base'], ds_path)
                        if service_name == 'http': http_url = '{}{}{}'.format(self._ngdc_url, service.attrib['base'], ds_path)

                    this_xml = iso_xml(iso_url)
                    title = this_xml.title()
                    h_epsg, v_epsg = this_xml.reference_system()
                    zv = this_xml.xml_doc.findall(
                        './/gmd:dimension/gmd:MD_Band/gmd:sequenceIdentifier/gco:MemberName/gco:aName/gco:CharacterString',
                        namespaces=this_xml.namespaces
                    )
                    if zv is not None:
                        for zvs in zv:
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
        
    def update(self):
        self.FRED._open_ds(1)
        self._parse_catalog(self._nt_catalog)
        self.FRED._close_ds()
    
    def run(self):
        """Search for data in the reference vector file"""
        
        for surv in FRED._filter_FRED(self):
            wcs_url = "{}?request=GetCoverage&version=1.0.0&service=WCS&coverage={}&bbox={}&format=geotiff_float"\
                .format(surv['IndexLink'], surv['Etcetra'], self.region.format('bbox'))
            if self.want_wcs:
                self.results.append(
                    [wcs_url, surv['DataLink'].split(',')[0].split('/')[-1].replace('.nc', '.tif'), surv['DataType']]
                )
            else:
                for d in surv['DataLink'].split(','):
                    if d != '':
                        self.results.append(
                            [d, d.split('/')[-1], surv['DataType']]
                        )

## ==============================================
## The National Map
## update is broken! fix this.
## ==============================================
class TheNationalMap(FetchModule):
    """USGS' The National Map

Fetch elevation data from The National Map
        
Various datasets from USGS's National Map. The National Map is a 
collaborative effort among the USGS and other Federal, State, and local partners to improve
and deliver topographic information for the Nation.

http://tnmaccess.nationalmap.gov/

< tnm:formats=None:extents=None:q=None >"""

    def __init__(self, where=[], formats=None, extents=None, q=None, **kwargs):
        super().__init__(name='tnm', **kwargs)
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

        self.where = [where] if len(where) > 0 else []        
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
        _req = Fetch(self._tnm_dataset_url).fetch_req()
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
                    Name = ds['sbDatasetTag'], ID = ds['id'], Agency = 'USGS', Date = pubDate[-4:],
                    MetadataLink = ds['infoUrl'], MetadataDate = metadataDate,
                    DataLink = '{}{}'.format(self._tnm_product_url, url_enc),
                    Link = ds['dataGovUrl'], Resolution = ','.join(ds['extents']),
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
                        #print(tag)
                        this_xml = iso_xml('{}?format=iso'.format(tag['infoUrl']))
                        geom = this_xml.bounds(geom=True)
                        h_epsg, v_epsg = this_xml.reference_system()
                        self._update_dataset(tag, fmt, geom, h_epsg, v_epsg)
                else:
                    this_xml = iso_xml('{}?format=iso'.format(ds['infoUrl']))
                    geom = this_xml.bounds(geom = True)
                    h_epsg, v_epsg = this_xml.reference_system()
                    self._update_dataset(ds, fmt, geom, h_epsg, v_epsg)
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
                    _data = {'bbox': self.region.format('bbox'), 'max': 100, 'offset': offset}
                    if q is not None: _data['q'] = str(q)
                    if f is None:
                        _data['prodFormats'] = surv['Etcetra']
                    else:
                        _data['prodFormats'] = ','.join(f)

                    if e is None: e = []

                    _req = Fetch(surv['DataLink']).fetch_req(params=_data)
                    if _req is not None and _req.status_code == 200:
                        try:
                            _dataset_results = _req.json()
                            total = _dataset_results['total']
                        except ValueError:
                            utils.echo_error_msg('tnm server error resulting in {}, try again'.format(e))
                        except Exception as e:
                            utils.echo_error_msg('error, {}'.format(e))

                    if len(_dataset_results) > 0:
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

                                        #self.results.append([f_url, os.path.join(self._outdir, os.path.join(*f_url.split('/')[:-1][3:]), f_url.split('/')[-1]), surv['DataType']])
                                        #self.results.append([f_url, os.path.join(self._outdir, surv['ID'].replace('-', '_'), f_url.split('/')[-1]), surv['DataType']])
                                        self.results.append([f_url, f_url.split('/')[-1], surv['DataType']])
                            else:
                                for fmt in fmts:
                                    if fmt in item['urls'].keys():
                                        f_url = item['urls'][fmt]
                                        break

                                if f_url is None:
                                    f_url = item['downloadURL']

                                #self.results.append([f_url, os.path.join(self._outdir, os.path.join(*f_url.split('/')[:-1][3:]), f_url.split('/')[-1]), surv['DataType']])
                                #self.results.append([f_url, os.path.join(self._outdir, surv['ID'].replace('-', '_'), f_url.split('/')[-1]), surv['DataType']])
                                self.results.append([f_url, f_url.split('/')[-1], surv['DataType']])

                    offset += 100
                    if offset >= total:
                        break
                
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
            _req = Fetch(self._tnm_product_url).fetch_req(params={'max': 1, 'datasets': dsTag})
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
                _req = Fetch(self._tnm_product_url).fetch_req(params = _data)
                try:
                    _dsTag_results = _req.json()
                except ValueError:
                    utils.echo_error_msg('tnm server error, try again')
                except Exception as e:
                    utils.echo_error_msg('error, {}'.format(e))
                
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
                           
    def _parse_prods_results(self, r, f = None, e = None, q = None):
        for surv in FRED._filter_FRED(self):
            for d in surv['DataLink'].split(','):
                if d != '':
                    self.results.append([d, os.path.join(self._outdir, d.split('/')[-1]), surv['DataType']])

## ==============================================
## EMODNet - EU data
## ==============================================
class EMODNet(FetchModule):
    """EU elevation data extracts from EMOD DTM.

Fetch raster data from the EMODNET DTM
    
https://portal.emodnet-bathymetry.eu/
https://erddap.emodnet.eu

erddap formats (default is csv):
https://erddap.emodnet.eu/erddap/griddap/documentation.html#fileType

Data
fileTypes	Description
.asc	View OPeNDAP-style ISO-8859-1 comma-separated text.
.csv	Download a ISO-8859-1 comma-separated text table (line 1: names; line 2: units; ISO 8601 times).
.csvp	Download a ISO-8859-1 .csv file with line 1: name (units). Times are ISO 8601 strings.
.csv0	Download a ISO-8859-1 .csv file without column names or units. Times are ISO 8601 strings.
.das	View the dataset's metadata via an ISO-8859-1 OPeNDAP Dataset Attribute Structure (DAS).
.dds	View the dataset's structure via an ISO-8859-1 OPeNDAP Dataset Descriptor Structure (DDS).
.dods	OPeNDAP clients use this to download the data in the DODS binary format.
.esriAscii	Download an ISO-8859-1 ESRI ASCII file (latitude longitude data only; longitude must be all below or all above 180).
.fgdc	View the dataset's UTF-8 FGDC .xml metadata.
.graph	View a Make A Graph web page.
.help	View a web page with a description of griddap.
.html	View an OPeNDAP-style HTML Data Access Form.
.htmlTable	View a UTF-8 .html web page with the data in a table. Times are ISO 8601 strings.
.iso19115	View the dataset's ISO 19115-2/19139 UTF-8 .xml metadata.
.itx	Download an ISO-8859-1 Igor Text File. Each axis variable and each data variable becomes a wave.
.json	View a table-like UTF-8 JSON file (missing value = 'null'; times are ISO 8601 strings).
.jsonlCSV1	View a UTF-8 JSON Lines CSV file with column names on line 1 (mv = 'null'; times are ISO 8601 strings).
.jsonlCSV	View a UTF-8 JSON Lines CSV file without column names (mv = 'null'; times are ISO 8601 strings).
.jsonlKVP	View a UTF-8 JSON Lines file with Key:Value pairs (missing value = 'null'; times are ISO 8601 strings).
.mat	Download a MATLAB binary file.
.nc	Download a NetCDF-3 binary file with COARDS/CF/ACDD metadata.
.ncHeader	View the UTF-8 header (the metadata) for the NetCDF-3 .nc file.
.ncml	View the dataset's structure and metadata as a UTF-8 NCML .xml file.
.nccsv	Download a NetCDF-3-like 7-bit ASCII NCCSV .csv file with COARDS/CF/ACDD metadata.
.nccsvMetadata	View the dataset's metadata as the top half of a 7-bit ASCII NCCSV .csv file.
.ncoJson	Download a UTF-8 NCO lvl=2 JSON file with COARDS/CF/ACDD metadata.
.odvTxt	Download time,latitude,longitude,otherVariables as an ODV Generic Spreadsheet File (.txt).
.timeGaps	View a UTF-8 list of gaps in the time values which are larger than the median gap.
.tsv	Download a ISO-8859-1 tab-separated text table (line 1: names; line 2: units; ISO 8601 times).
.tsvp	Download a ISO-8859-1 .tsv file with line 1: name (units). Times are ISO 8601 strings.
.tsv0	Download a ISO-8859-1 .tsv file without column names or units. Times are ISO 8601 strings.
.wav	Download a .wav audio file. All columns must be numeric and of the same type.
.xhtml	View a UTF-8 XHTML (XML) file with the data in a table. Times are ISO 8601 strings.

< emodnet >"""

    def __init__(self, want_erddap=True, erddap_format='csv', **kwargs):
        super().__init__(name='emodnet', **kwargs)
        self.want_erddap = want_erddap
        self.erddap_format = erddap_format
        self._emodnet_grid_url = 'https://ows.emodnet-bathymetry.eu/wcs?'
        self._emodnet_grid_url_erddap = 'https://erddap.emodnet.eu/erddap/griddap/dtm_2020_v2_e0bf_e7e4_5b8f.{}?'.format(erddap_format)
        
    def run(self):
        """Run the EMODNET fetching module"""
        
        if self.region is None: return([])

        if self.want_erddap:
            suff = 'elevation%5B({}):1:({})%5D%5B({}):1:({})%5D'.format(self.region.ymin, self.region.ymax, self.region.xmin, self.region.xmax)
            erddap_url = self._emodnet_grid_url_erddap + suff
            outf = 'emodnet_{}.{}'.format(self.region.format('fn'), self.erddap_format)
            self.results.append([erddap_url, outf, self.erddap_format])
            
        else:
            _data = {'request': 'DescribeCoverage', 'version': '2.0.1',
                     'CoverageID': 'emodnet:mean', 'service': 'WCS'}
            _req = Fetch(self._emodnet_grid_url).fetch_req(params=_data)
            _results = lxml.etree.fromstring(_req.text.encode('utf-8'))
            g_env = _results.findall('.//{http://www.opengis.net/gml/3.2}GridEnvelope', namespaces=namespaces)[0]
            hl = [float(x) for x in g_env.find('{http://www.opengis.net/gml/3.2}high').text.split()]

            g_bbox = _results.findall('.//{http://www.opengis.net/gml/3.2}Envelope')[0]
            lc = [float(x) for x in  g_bbox.find('{http://www.opengis.net/gml/3.2}lowerCorner').text.split()]
            uc = [float(x) for x in g_bbox.find('{http://www.opengis.net/gml/3.2}upperCorner').text.split()]

            ds_region = regions.Region().from_list(
                [lc[1], uc[1], lc[0], uc[0]]
            )
            resx = (uc[1] - lc[1]) / hl[0]
            resy = (uc[0] - lc[0]) / hl[1]
            if regions.regions_intersect_ogr_p(self.region, ds_region):
                emodnet_wcs = '{}service=WCS&request=GetCoverage&version=1.0.0&Identifier=emodnet:mean&coverage=emodnet:mean&format=GeoTIFF&bbox={}&resx={}&resy={}&crs=EPSG:4326'\
                                          .format(self._emodnet_grid_url, self.region.format('bbox'), resx, resy)
                outf = 'emodnet_{}.tif'.format(self.region.format('fn'))
                self.results.append([emodnet_wcs, outf, 'emodnet'])
            
        return(self)

## ==============================================
## CHS - Canada Hydro
## ==============================================
class CHS(FetchModule):
    """High-Resolution Digital Elevation Model data for Canada

Fetch bathymetric soundings from the CHS
    
https://open.canada.ca

< chs >"""
    
    def __init__(self, **kwargs):
        super().__init__(name='chs', **kwargs)
        self._chs_api_url = "https://geoportal.gc.ca/arcgis/rest/services/FGP/CHS_NONNA_100/MapServer/0/query?"
        #self._chs_url = 'https://data.chs-shc.ca/geoserver/wcs?'
        self._chs_url = 'https://nonna-geoserver.data.chs-shc.ca/geoserver/wcs?' # new
        
    def run(self):
        """Run the CHS fetching module"""
        
        if self.region is None: return([])
        _data = {'request': 'DescribeCoverage', 'version': '2.0.1', 'CoverageID': 'nonna__NONNA 10 Coverage',
                 'service': 'WCS'}
        _req = Fetch(self._chs_url).fetch_req(params=_data)
        _results = lxml.etree.fromstring(_req.text.encode('utf-8'))
        print(lxml.etree.tostring(_results))
        
        g_env = _results.findall('.//{http://www.opengis.net/gml/3.2}GridEnvelope', namespaces=namespaces)[0]
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
            _wcs_data = {'request': 'GetCoverage', 'version': '2.0.1', 'CoverageID': 'nonna__NONNA 10 Coverage',
                         'service': 'WCS', 'crs': '4326', 'bbox': self.region.format('bbox'),
                         'resx': resx, 'resy': resy}

            utils.echo_msg(_wcs_data)
            _wcs_req = Fetch(self._chs_url).fetch_req(params=_wcs_data)
            #chs_wcs = '{}service=WCS&request=GetCoverage&version=2.0.1&CoverageID=nonna__NONNA+10+Coverage&format=BAG&bbox={}&resx={}&resy={}&crs=EPSG:4326'\
            #chs_wcs = '{}service=WCS&request=GetCoverage&version=2.0.1&CoverageID=nonna__NONNA+10+Coverage&bbox={}&crs=EPSG:4326'\
            #                 .format(self._chs_url, self.region.format('bbox'), resx, resy)
            outf = 'chs_{}.tif'.format(self.region.format('fn'))
            self.results.append([_wcs_req.url, outf, 'chs'])
        return(self)

## ==============================================
## HRDEM - Canada Topo
## ==============================================    
class HRDEM(FetchModule):
    """High-Resolution Digital Elevation Model data for Canada

Fetch HRDEM data from Canada (NRCAN)
    
https://open.canada.ca

< hrdem >"""
    
    def __init__(self, **kwargs):
        super().__init__(name='hrdem', **kwargs)
        self._hrdem_footprints_url = 'ftp://ftp.maps.canada.ca/pub/elevation/dem_mne/highresolution_hauteresolution/Datasets_Footprints.zip'
        self._hrdem_info_url = 'https://open.canada.ca/data/en/dataset/957782bf-847c-4644-a757-e383c0057995#wb-auto-6'

    def run(self):
        v_zip = os.path.join(self._outdir, 'Datasets_Footprints.zip')
        status = Fetch(self._hrdem_footprints_url, verbose=self.verbose).fetch_ftp_file(v_zip)
        v_shps = utils.p_unzip(v_zip, ['shp', 'shx', 'dbf', 'prj'])
        v_shp = None
        for v in v_shps:
            if v.split('.')[-1] == 'shp':
                v_shp = v
                break
        try:
            v_ds = ogr.Open(v_shp)
        except:
            v_ds = None
            status = -1
                
        if v_ds is not None:
            layer = v_ds.GetLayer()
            fcount = layer.GetFeatureCount()
            for f in range(0, fcount):
                feature = layer[f]
                geom = feature.GetGeometryRef()
                if geom.Intersects(self.region.export_as_geom()):
                    data_link = feature.GetField('Ftp_dtm')
                    self.results.append([data_link, data_link.split('/')[-1], 'raster'])
            v_ds = None

        utils.remove_glob(v_zip, *v_shps)

## ==============================================
## ArcticDEM
## ==============================================
class ArcticDEM(FetchModule):
    """Arctic DEM
ArcticDEM is an NGA-NSF public-private initiative to automatically produce a high-resolution, 
high quality, digital surface model (DSM) of the Arctic using optical stereo imagery, 
high-performance computing, and open source photogrammetry software.

 objectid (Integer64)
 name (String)
 tile (String)
 nd_value (Real)
 resolution (Real)
 creationda (Date)
 raster (String)
 fileurl (String)
 spec_type (String)
 qual (Real)
 reg_src (String)
 num_gcps (Integer)
 meanresz (Real)
 active (Integer)
 qc (Integer)
 rel_ver (String)
 num_comp (Integer)
    
https://www.pgc.umn.edu/data/arcticdem/

< arcticdem >"""
    
    def __init__(self, where='1=1', layer=0, **kwargs):
        super().__init__(name='arcticdem', **kwargs)
        self._arctic_dem_index_url = 'https://data.pgc.umn.edu/elev/dem/setsm/ArcticDEM/indexes/ArcticDEM_Tile_Index_Rel7.zip'
        self.where = [where] if len(where) > 0 else []
        self.arctic_region = self.region.copy()
        self.arctic_region.warp('epsg:3413')
        
        self.FRED = FRED.FRED(name=self.name, verbose = self.verbose)
        self.update_if_not_in_FRED()

    def update_if_not_in_FRED(self):
        self.FRED._open_ds()
        self.FRED._attribute_filter(["DataSource = '{}'".format(self.name)])
        if len(self.FRED.layer) == 0:
            self.FRED._close_ds()
            self.update()
        self.FRED._close_ds()
        
    def update(self):
        self.FRED._open_ds()        
        v_zip = os.path.basename(self._arctic_dem_index_url)
        try:
            status = Fetch(self._arctic_dem_index_url, verbose=self.verbose).fetch_file(v_zip)
        except:
            status = -1
            
        v_shps = utils.p_unzip(v_zip, ['shp', 'shx', 'dbf', 'prj'])

        v_shp = None
        for v in v_shps:
            if '.shp' in v:
                v_shp = v
                break

        utils.run_cmd('ogr2ogr arctic_tmp.shp {} -t_srs epsg:4326'.format(v_shp), verbose=self.verbose)
        utils.remove_glob(v_zip, *v_shps)
        v_shp = 'arctic_tmp.shp'
        v_shps = ['arctic_tmp.shp','arctic_tmp.dbf','arctic_tmp.shx','arctic_tmp.prj']
        shp_regions = regions.gdal_ogr_regions(v_shp)
        shp_region = regions.Region()
        for this_region in shp_regions:
            #this_region.src_srs = 'epsg:3413'
            #this_region.warp('epsg:4326')
            if shp_region.valid_p(check_xy=True):
                shp_region = regions.regions_merge(shp_region, this_region)
            else: shp_region = this_region
        geom = shp_region.export_as_geom()
        
        self.FRED._attribute_filter(["ID = '{}'".format('ARCTICDEM-1')])
        if self.FRED.layer is None or len(self.FRED.layer) == 0:
            self.FRED._add_survey(Name = 'ArcticDEM', ID = 'ARCTICDEM-1', Agency = 'UMN', Date = utils.this_year(),
                                  MetadataLink = 'https://data.pgc.umn.edu/elev/dem/setsm/ArcticDEM/', MetadataDate = utils.this_year(),
                                  DataLink = self._arctic_dem_index_url, IndexLink = self._arctic_dem_index_url,
                                  DataType = 'raster', DataSource = 'arcticdem', Info = 'Arctic Only', geom = geom)
        utils.remove_glob(*v_shps)
        self.FRED._close_ds()

    def run(self):
        #for surv in FRED._filter_FRED(self):

        v_zip = os.path.join(self._outdir, os.path.basename(self._arctic_dem_index_url))
        try:
            status = Fetch(self._arctic_dem_index_url, verbose=self.verbose).fetch_file(v_zip)
        except:
            status = -1
            
        v_shps = utils.p_unzip(v_zip, ['shp', 'shx', 'dbf', 'prj'], outdir=self._outdir)
        v_shp = None
        for v in v_shps:
            if v.split('.')[-1] == 'shp':
                v_shp = v
                break

        #v_out_shp = os.path.join(self._outdir, 'arctic_tmp.shp')
        #utils.run_cmd('ogr2ogr {} {} -t_srs epsg:4326'.format(v_out_shp, v_shp), verbose=True)
        #utils.remove_glob(v_zip, *v_shps)
        #v_shps = ['arctic_tmp.shp','arctic_tmp.dbf','arctic_tmp.shx','arctic_tmp.prj']
        try:
            v_ds = ogr.Open(v_shp)
        except:
            v_ds = None
            status = -1

        if v_ds is not None:
            layer = v_ds.GetLayer()
            _boundsGeom = self.arctic_region.export_as_geom()
            layer.SetSpatialFilter(_boundsGeom)
            fcount = layer.GetFeatureCount()
            utils.echo_msg('filtered {} arcticdem features'.format(fcount))
            for f in range(0, fcount):
                feature = layer[f]
                data_link = feature.GetField('fileurl')
                self.results.append([data_link, data_link.split('/')[-1], 'raster'])

            v_ds = None
            
        return(self)

## ==============================================
## OSM - Open Street Map
## ==============================================
class OpenStreetMap(FetchModule):
    """OpenStreetMap data.
    
OpenStreetMap is a free, editable map of the whole world that is 
being built by volunteers largely from scratch and released with an 
open-content license.

https://wiki.openstreetmap.org/

< osm:q=None:fmt=osm:planet=False:chunks=True:min_length=None >"""
    
    def __init__(self, q=None, fmt='osm', planet=False, chunks=True, min_length=None, **kwargs):
        super().__init__(name='osm', **kwargs)
        self._osm_api = 'https://lz4.overpass-api.de/api/interpreter'
        self._osm_api2 = 'https://overpass.kumi.systems/api/interpreter'
        self._osm_api3 = 'https://overpass.openstreetmap.fr/api/interpreter'
        self._osm_planet_bz2 = 'https://ftpmirror.your.org/pub/openstreetmap/planet/planet-latest.osm.bz2'
        self._osm_planet = 'https://ftpmirror.your.org/pub/openstreetmap/pbf/planet-latest.osm.pbf'
        self._osm_continents = 'https://download.geofabrik.de/'
        self.q = q
        self.fmt = fmt
        self.planet = planet
        self.chunks = chunks

        self.h = ''
        
        if self.q == 'buildings':
            #self.h = '[maxsize:2000000000]'
            self.h = '[timeout:3600]'
            self.q = '''
            (way["building"]{};
            relation["building"]["type"="multipolygon"];
            );
            (._;>;);
            out meta;
            '''.format('(if: length() > {})'.format(min_length) if min_length is not None else '')
        
        self.headers = { 'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:89.0) Gecko/20100101 Firefox/89.0',
                         'referer': 'https://lz4.overpass-api.de/' }
        
    def run(self):
        if self.region is None:
            return([])

        #if self.planet:
        #    self.results.append([self._osm_planet, os.path.join(self._outdir, 'planet-latest.osm.bz2'), 'bz2'])
        
        ## fetch whole planet
        if self.planet:
            self.results.append([self._osm_planet, os.path.join(self._outdir, 'planet-latest.osm.pbf'), 'pbf'])

        ## fetch in chunks
        elif self.chunks:
            x_delta = self.region.xmax - self.region.xmin
            y_delta = self.region.ymax - self.region.ymin
            incs = self.region.increments(1000,1000)

            ## break up the requests into .05 degree chunks for
            ## better usage of the OSM API
            if x_delta > .25 or y_delta > .25:
                xcount, ycount, gt = self.region.geo_transform(x_inc=incs[0], y_inc=incs[1])
                if x_delta >= y_delta:
                    n_chunk = int(xcount*(.25/x_delta))
                elif y_delta > x_delta:
                    n_chunk = int(ycount*(.25/y_delta))
            else:
                n_chunk = None

            these_regions = self.region.chunk(incs[0], n_chunk=n_chunk)
            utils.echo_msg('chunking OSM request into {} regions'.format(len(these_regions)))
            
            for this_region in these_regions:
                c_bbox = this_region.format('osm_bbox')
                out_fn = 'osm_{}'.format(this_region.format('fn_full'))
                osm_q_bbox  = '''
                {1}{2}[bbox:{0}];'''.format(c_bbox, '[out:{}]'.format(self.fmt) if self.fmt != 'osm' else '', self.h)

                osm_q = '''
                (node;
                <;
                >;
                );
                out meta;
                '''
                
                osm_q_ = osm_q_bbox + (osm_q if self.q is None else self.q)

                #utils.echo_msg('using query: {}'.format(osm_q_))
                osm_data = urlencode({'data': osm_q_})
                osm_data_url = self._osm_api + '?' + osm_data
                
                self.results.append([osm_data_url, '{}.{}'.format(out_fn, self.fmt), 'osm'])

        else:
            c_bbox = self.region.format('osm_bbox')
            out_fn = 'osm_{}'.format(self.region.format('fn_full'))
            osm_q_bbox  = '''
            {1}[bbox:{0}];'''.format(c_bbox, '[out:{}]'.format(self.fmt) if self.fmt != 'osm' else '')

            osm_q = '''
            (node;
            <;
            >;
            );
            out meta;
            '''

            osm_q_ = osm_q_bbox + (osm_q if self.q is None else self.q)
            osm_data = urlencode({'data': osm_q_})
            osm_data_url = self._osm_api + '?' + osm_data            
            self.results.append([osm_data_url, '{}.{}'.format(out_fn, self.fmt), 'osm'])

## ==============================================
## BING Building Footprints
## ==============================================
class BingBFP(FetchModule):
    """Bing Building Footprints

https://github.com/microsoft/GlobalMLBuildingFootprints
"""
    
    def __init__(self, **kwargs):
        super().__init__(name='bingbfp', **kwargs)
        self._bing_bfp_csv = 'https://minedbuildings.blob.core.windows.net/global-buildings/dataset-links.csv'
        
        self.headers = { 'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:89.0) Gecko/20100101 Firefox/89.0',
                         'referer': 'https://lz4.overpass-api.de/' }
        
    def run(self):
        if self.region is None:
            return([])

        bbox = self.region.format('bbox')

        quad_keys = set()
        for tile in list(mercantile.tiles(self.region.xmin, self.region.ymin, self.region.xmax, self.region.ymax, zooms=9)):
            quad_keys.add(int(mercantile.quadkey(tile)))
            
        quad_keys = list(quad_keys)
        print(f"The input area spans {len(quad_keys)} tiles: {quad_keys}")

        bing_csv = os.path.basename(self._bing_bfp_csv)
        try:
            status = Fetch(self._bing_bfp_csv, verbose=self.verbose).fetch_file(bing_csv)
        except:
            status = -1

        with open(bing_csv, mode='r') as bc:
            reader = csv.reader(bc)
            next(reader)
            bd = [row[2] for row in reader if int(row[1]) in quad_keys]

        self.results = [[url, os.path.basename(url), 'bing'] for url in bd]
        
## ==============================================
## VDATUM
## ==============================================
def proc_vdatum_inf(vdatum_inf, name='vdatum'):
    
    _inf = open(vdatum_inf, 'r')
    _inf_areas = {}

    for line in _inf:
        line_list = line.split('.')
        if len(line_list) > 1:
            if line_list[0] not in _inf_areas.keys():
                _inf_areas[line_list[0]] = {}
            if len(line_list) > 1:
                line_val = '.'.join(line_list[1:]).strip()
                utils.args2dict([line_val], _inf_areas[line_list[0]])
    _inf.close()

    _inf_areas_fmt = {}

    for key in _inf_areas.keys():
        if name is not None:
            _out_key = '{}_{}'.format(name, key)
        else:
            _out_key = key
            
        if _out_key not in _inf_areas_fmt.keys():
            _inf_areas_fmt[_out_key] = {}
            
        xmin = utils.x360(float(_inf_areas[key]['minlon']))
        xmax = utils.x360(float(_inf_areas[key]['maxlon']))
        ymin = float(_inf_areas[key]['minlat'])
        ymax = float(_inf_areas[key]['maxlat'])
        _inf_areas_fmt[_out_key]['region'] = [xmin, xmax, ymin, ymax]
        _inf_areas_fmt[_out_key]['grid'] = _inf_areas[key]['source'].split('\\')[-1]

        if 'AK' in _out_key: print(_inf_areas_fmt[_out_key])
        
    return(_inf_areas_fmt)

def search_proj_cdn(region, epsg=None, crs_name=None, name=None, verbose=True, cache_dir='./'):
    """Search PROJ CDN for transformation grids:
    the PROJ CDN holds transformation grids from around the
    world, including global transformations such as EGM
    """
    
    _proj_vdatum_index = 'https://cdn.proj.org/files.geojson'
    cdn_index = os.path.join(cache_dir, 'proj_cdn_files.geojson')
    #cdn_headers = { 'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:89.0) Gecko/20100101 Firefox/89.0' }
    cdn_headers = {'Host': 'cdn.proj.org',
                   'User-Agent': 'Mozilla/5.0 (X11; Linux x86_64; rv:45.0) Gecko/20100101 Firefox/45.0',
                   'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8',
                   'Accept-Language': 'pt-BR,pt;q=0.8,en-US;q=0.5,en;q=0.3',
                   'Accept-Encoding': 'gzip, deflate',
                   'Connection': 'keep-alive',
                   'Pragma': 'no-cache',
                   'Cache-Control': 'no-cache'}

    try:
        status = Fetch(_proj_vdatum_index, headers=cdn_headers, verbose=verbose).fetch_file(cdn_index, timeout=5, read_timeout=5, check_size=False)
    except:
        status = -1

    if status == 0:
        cdn_driver = ogr.GetDriverByName('GeoJSON')
        cdn_ds = cdn_driver.Open(cdn_index, 0)
        cdn_layer = cdn_ds.GetLayer()
        _boundsGeom = region.export_as_geom()
        _results = []

        if crs_name is not None:
            cdn_layer.SetAttributeFilter(
                "type != 'HORIZONTAL_OFFSET' AND (target_crs_name LIKE '%{}%' OR source_crs_name LIKE '%{}%')".format(
                    name.upper(), name.upper()
                )
            )
        elif epsg is not None:
            cdn_layer.SetAttributeFilter(
                "type != 'HORIZONTAL_OFFSET' AND (target_crs_code LIKE '%{}%' OR source_crs_code LIKE '%{}%')".format(
                    epsg, epsg
                )
            )
        elif name is not None:
            cdn_layer.SetAttributeFilter("type != 'HORIZONTAL_OFFSET' AND name LIKE '%{}%'".format(name))
        else:
            cdn_layer.SetAttributeFilter("type != 'HORIZONTAL_OFFSET'")

        for feat in cdn_layer:
            if _boundsGeom is not None:
                geom = feat.GetGeometryRef()
                if geom is not None:
                    if _boundsGeom.Intersects(geom):
                        _results.append({})
                        f_j = json.loads(feat.ExportToJson())
                        for key in f_j['properties'].keys():
                            _results[-1][key] = feat.GetField(key)
            else:
                _results.append({})
                f_j = json.loads(feat.ExportToJson())
                for key in f_j['properties'].keys():
                    _results[-1][key] = feat.GetField(key)

        cdn_ds = None
        #utils.remove_glob(cdn_index)
        return(_results)

class VDATUM(FetchModule):
    """NOAA's VDATUM transformation grids

Fetch vertical datum conversion grids from NOAA, etc.
    
VDatum is a free software tool being developed jointly by NOAA's National Geodetic Survey (NGS), 
Office of Coast Survey (OCS), and Center for Operational Oceanographic Products and Services (CO-OPS). 

VDatum is designed to vertically transform geospatial data among a variety of tidal, orthometric and 
ellipsoidal vertical datums - allowing users to convert their data from different horizontal/vertical 
references into a common system and enabling the fusion of diverse geospatial data in desired reference 
levels.

https://vdatum.noaa.gov
https://cdn.proj.org

< vdatum:datatype=None:gtx=False >"""
    

    _tidal_references = {
        1089: {'name': 'mllw',
               'description': 'Mean Lower Low Water Height',
               'grid': 'mllw.gtx'},
        1091: {'name': 'mlw',
               'description': 'Mean Low Water Height',
               'grid': 'mlw.gtx'},
        5868: {'name': 'mhw',
               'description': 'Mean High Water',
               'grid': 'mhw.gtx'},
        5869: {'name': 'mhhw',
               'description': 'Mean Higher High Water',
               'grid': 'mhhw.gtx'},
        5703: {'name': 'tss',
               'description': 'NAVD88 tss geoid',
               'grid': 'tss.gtx'},
        6641: {'name': 'tss',
               'description': 'PRVD02 tss geoid',
               'grid': 'tss.gtx'},
        6642: {'name': 'tss',
               'description': 'VIVD09 tss geoid',
               'grid': 'tss.gtx'},
        5714: {'name': 'tss',
               'description': 'to MSL tss geoid',
               'grid': 'tss.gtx'},
    }
    
    def __init__(self, where=[], datatype=None, gtx=False, epsg=None, **kwargs):
        super().__init__(name='vdatum', **kwargs)
        
        self._vdatum_data_url = 'https://vdatum.noaa.gov/download/data/'
        self._proj_vdatum_index = 'https://cdn.proj.org/files.geojson'
        
        ## add others IGLD85
        #self._vdatums = ['VERTCON', 'EGM1984', 'EGM1996', 'EGM2008', 'GEOID03', 'GEOID06', 'GEOID09', 'GEOID12A', 'GEOID12B', 'GEOID96', 'GEOID99', 'TIDAL']
        self._vdatums = ['TIDAL']
        self._tidal_datums = ['mhw', 'mhhw', 'mlw', 'mllw', 'tss', 'mtl']
        self.where = where
        self.datatype = datatype
        self.epsg = utils.int_or(epsg)
        self.gtx = gtx
        self.v_datum = 'varies'
        self.FRED = FRED.FRED(name=self.name, verbose=self.verbose)
        self.update_if_not_in_FRED()
        self.src_srs = 'epsg:4326'
        
    def update_if_not_in_FRED(self):
        self.FRED._open_ds()
        self.FRED._attribute_filter(["DataSource = '{}'".format(self.name)])
        if len(self.FRED.layer) == 0:
            self.FRED._close_ds()
            self.update()
            
        self.FRED._close_ds()

    def update(self):
        """Update or create the reference vector file"""

        self.FRED._open_ds(1)
        for vd in self._vdatums:
            surveys = []

            if vd == 'TIDAL' or vd == 'IGLD85':
                ## ==============================================
                ## All tidal inf data is in each one, so we only
                ## have to download one of the tidal zips to process
                ## them all; lets use the smallest one
                ## Keep this link up-to-date!
                ## ==============================================
                if vd == 'TIDAL':
                    vd_ = 'DEVAemb12_8301'
                else:
                    vd_ = vd
                    
                vd_zip_url = '{}{}.zip'.format(self._vdatum_data_url, vd_)
                v_inf = 'tidal_area.inf'
            elif vd == 'VERTCON':
                vd_zip_url = '{}vdatum_{}.zip'.format(self._vdatum_data_url, vd)
                v_inf = 'vcn.inf'
            else:
                vd_zip_url = '{}vdatum_{}.zip'.format(self._vdatum_data_url, vd)
                v_inf = '{}.inf'.format(vd.lower())

            try:
                status = Fetch(vd_zip_url, verbose=True).fetch_file('{}.zip'.format(vd))
            except:
                status = -1
                
            if status == 0:
                v_infs = utils.p_unzip('{}.zip'.format(vd), ['inf'])
                utils.echo_msg(v_infs)
                v_dict = proc_vdatum_inf(v_inf, name=vd if vd != 'TIDAL' else None)#, loff=-360 if vd =='TIDAL' else -180)
                v_dict = proc_vdatum_inf(v_inf, name=vd if vd != 'TIDAL' else None)#, loff=-360)

                for key in v_dict.keys():
                    v_dict[key]['vdatum'] = vd
                    v_dict[key]['remote'] = vd_zip_url
                
                ## tidal datums:
                if vd == 'TIDAL':
                    v_dict_ = {}
                    for tidal_key in v_dict.keys():                        
                        for t in self._tidal_datums:
                            key_ = '{}_{}'.format(t, tidal_key)
                            v_dict_[key_] = {}
                            v_dict_[key_]['region'] = v_dict[tidal_key]['region']
                            v_dict_[key_]['vdatum'] = t
                            v_dict_[key_]['grid'] = '{}.gtx'.format(t)
                            v_dict_[key_]['remote'] = '{}{}.zip'.format(self._vdatum_data_url, tidal_key)
                            
                    v_dict = v_dict_
                    
                print(v_dict)

                for key in v_dict.keys():
                    self.FRED._attribute_filter(["ID = '{}'".format(key)])
                    if self.FRED.layer is None or len(self.FRED.layer) == 0:
                        geom = regions.Region().from_list(v_dict[key]['region']).export_as_geom()
                        if geom is not None:
                            surveys.append({'Name': v_dict[key]['grid'], 'ID': key, 'Agency': 'NOAA', 'Date': utils.this_date(),
                                            'MetadataLink': "", 'MetadataDate': utils.this_date(),
                                            'DataLink': v_dict[key]['remote'], 'Link': self._vdatum_data_url, 'DataType': v_dict[key]['vdatum'],
                                            'DataSource': 'vdatum', 'HorizontalDatum': 4326, 'VerticalDatum': v_dict[key]['vdatum'],
                                            'Info': "", 'geom': geom})
            self.FRED._add_surveys(surveys)
            #utils.remove_glob(*v_infs, '{}.zip'.format(vd))
            utils.remove_glob(*v_infs)
            
        self.FRED._close_ds()

    def run(self):
        """Search for data in the reference vector file.
        If self.gtx is true, will download the zip and extract the 
        appropriate gtx file.
        """

        w = []
        if self.datatype is not None:
            w.append("DataType = '{}'".format(self.datatype))
        elif self.epsg is not None:
            w.append("DataType = '{}'".format(self._tidal_references[self.epsg]['name']))
            
        ## ==============================================
        ## Search FRED for VDATUM TIDAL TRANSFORMATION GRIDS
        ## FRED holds the VDATUM tidal grids,
        ## mllw, mlw, mhhw, mhw, tss, mtl
        ## where all convert to MSL and tss represents the
        ## current geoid
        ## ==============================================
        for surv in self.FRED._filter(self.region, w, [self.name]):
            if self.gtx:
                dst_zip = '{}.zip'.format(surv['ID'])
                try:
                    status = Fetch(surv['DataLink'], callback=self.callback, verbose=self.verbose).fetch_file(dst_zip)
                except:
                    status = -1
                    
                if status == 0:
                    v_gtxs = utils.p_f_unzip(dst_zip, [surv['Name']])
                    for v_gtx in v_gtxs:
                        os.replace(v_gtx, '{}.gtx'.format(surv['ID']))
                    #utils.remove_glob(dst_zip)
            else:
                self.results.append([surv['DataLink'], '{}.zip'.format(surv['ID']), surv['Name']])

        ## ==============================================
        ## Search PROJ CDN for all other transformation grids:
        ## the PROJ CDN holds transformation grids from around the
        ## world, including global transformations such as EGM
        ## ==============================================
        cdn_index = 'proj_cdn_files.geojson'
        try:
            status = Fetch(self._proj_vdatum_index, callback=self.callback, verbose=self.verbose).fetch_file(cdn_index)
        except:
            status = -1
            
        if status == 0:
            cdn_driver = ogr.GetDriverByName('GeoJSON')
            cdn_ds = cdn_driver.Open(cdn_index, 0)
            cdn_layer = cdn_ds.GetLayer()
            _boundsGeom = self.region.export_as_geom()
            _results = []

            if self.datatype is not None:
                cdn_layer.SetAttributeFilter(
                    "type != 'HORIZONTAL_OFFSET' AND (target_crs_name LIKE '%{}%' OR source_crs_name LIKE '%{}%')".format(
                        self.datatype.upper(), self.datatype.upper()
                    )
                )
            elif self.epsg is not None:
                cdn_layer.SetAttributeFilter(
                    "type != 'HORIZONTAL_OFFSET' AND (target_crs_code LIKE '%{}%' OR source_crs_code LIKE '%{}%')".format(
                        self.epsg, self.epsg
                    )
                )
            else:
                cdn_layer.SetAttributeFilter("type != 'HORIZONTAL_OFFSET'")
             
            for feat in cdn_layer:
                if _boundsGeom is not None:
                    geom = feat.GetGeometryRef()
                    if geom is not None:
                        if _boundsGeom.Intersects(geom):
                            _results.append({})
                            f_j = json.loads(feat.ExportToJson())
                            for key in f_j['properties'].keys():
                                _results[-1][key] = feat.GetField(key)
                else:
                    _results.append({})
                    f_j = json.loads(feat.ExportToJson())
                    for key in f_j['properties'].keys():
                        _results[-1][key] = feat.GetField(key)
                        
            for _result in _results:
                self.results.append([_result['url'], _result['name'], _result['source_id']])
                
            cdn_ds = None
            utils.remove_glob(cdn_index)

## ==============================================
## EarthData - NASA (requires login credentials)
## ==============================================
class EarthData(FetchModule):
    """ACCESS NASA EARTH SCIENCE DATA
    
NASA promotes the full and open sharing of all its data to research and applications communities, 
private industry, academia, and the general public. In order to meet the needs of these different 
communities, NASA’s Earth Observing System Data and Information System (EOSDIS) has provided various 
ways to discover, access, and use the data.

If version is omitted, will fetch all versions
Use wildcards in 'short_name' to return granules for all matching short_name entries.

Commentary from nsidc_download.py
Tested in Python 2.7 and Python 3.4, 3.6, 3.7

To run the script at a Linux, macOS, or Cygwin command-line terminal:
  $ python nsidc-data-download.py

On Windows, open Start menu -> Run and type cmd. Then type:
    python nsidc-data-download.py

The script will first search Earthdata for all matching files.
You will then be prompted for your Earthdata username/password
and the script will download the matching files.

If you wish, you may store your Earthdata username/password in a .netrc
file in your $HOME directory and the script will automatically attempt to
read this file. The .netrc file should have the following format:
   machine urs.earthdata.nasa.gov login myusername password mypassword
where 'myusername' and 'mypassword' are your Earthdata credentials.

you might need to `chmod 0600 ~/.netrc`

NASA promotes the full and open sharing of all its data to research and applications communities, 
private industry, academia, and the general public. In order to meet the needs of these different 
communities, NASA’s Earth Observing System Data and Information System (EOSDIS) has provided various 
ways to discover, access, and use the data.

nsidc_download.py updated for fetches integration 12/21
Updated from nsidc.py to earthdata.py to support all earthdata datasets

some notable datasets:
    ATL03
    ATL06
    ATL07
    ATL08
    GLAH06
    GLAH14
    GEDI01_B
    GEDI02_A
    GEDI02_B
    ASTGTM
    ASTGTM_NC
    ASTL1A
    ASTL1T
    SRTMGL1
    SRTMGL1_NC
    SRTMGL3
    SRTMGL3S
    SRTMGL30
    NASADEM_HGT
    NASADEM_NC
    NASADEM_SHHP
    NASADEM_NUMNC
    NASADEM_SIM
    
https://cmr.earthdata.nasa.gov

< earthdata:short_name=ATL08:version=004:time_start='':time_end='':filename_filter='' >"""

    def __init__(self, short_name='ATL03', provider='', time_start='', time_end='', version='', filename_filter=None, **kwargs):
        super().__init__(name='cmr', **kwargs)
        self._cmr_url = 'https://cmr.earthdata.nasa.gov/search/granules.json?'
        self.short_name = short_name
        self.provider = provider
        self.time_start = time_start
        self.time_end = time_end
        self.version = version
        self.filename_filter = filename_filter

        credentials = get_credentials(None)
        self.headers = {'Authorization': 'Basic {0}'.format(credentials)}

    def add_wildcards_to_str(self, in_str):
        if not in_str.startswith('*'):
            in_str = '*' + in_str
            
        if not in_str.endswith('*'):
            in_str = in_str + '*'
            
        return(in_str)
        
    def run(self):
        if self.region is None:
            return([])

        _data = {
            'provider': self.provider,
            'short_name': self.short_name,
            'bounding_box': self.region.format('bbox'),
            'temporal': f'{self.time_start}, {self.time_end}',
            'page_size': 2000,
        }

        if '*' in self.short_name:
            _data['options[short_name][pattern]'] = 'true'
        
        if self.filename_filter is not None:
            _data['options[producer_granule_id][pattern]'] = 'true'
            filename_filters = self.filename_filter.split(',')
            for filename_filter in filename_filters:
                _data['producer_granule_id'] = self.add_wildcards_to_str(filename_filter)
                
        _req = Fetch(self._cmr_url).fetch_req(params=_data)
        if _req is not None:
            features = _req.json()['feed']['entry']
            for feature in features:
                if 'polygons' in feature.keys():
                    poly = feature['polygons'][0][0]
                    cc = [float(x) for x in poly.split()]
                    gg = [x for x in zip(cc[::2], cc[1::2])]
                    ogr_geom = ogr.CreateGeometryFromWkt(regions.create_wkt_polygon(gg))
                    ## ==============================================
                    ## uncomment below to output shapefiles of the feature polygons
                    ## ==============================================
                    #regions.write_shapefile(ogr_geom, '{}.shp'.format(feature['title']))
                else:
                    ogr_geom = self.region.export_as_geom()

                if self.region.export_as_geom().Intersects(ogr_geom):
                    links = feature['links']
                    for link in links:
                        if link['rel'].endswith('/data#') and 'inherited' not in link.keys():
                            if not any([link['href'].split('/')[-1] in res for res in self.results]):
                                self.results.append([link['href'], link['href'].split('/')[-1], self.short_name])

## ==============================================
## IceSat2 from EarthData shortcut - NASA (requires login credentials)
##
## This module allows us to use icesat2 data in dlim/waffles
##
## todo: dmrpp
## ==============================================
class IceSat(EarthData):
    """Access IceSat2 data.

By default access ATL03 data, specify 'short_name' to fetch specific ATL data.

If you wish, you may store your Earthdata username/password in a .netrc
file in your $HOME directory and the script will automatically attempt to
read this file. The .netrc file should have the following format:
   machine urs.earthdata.nasa.gov login myusername password mypassword
where 'myusername' and 'mypassword' are your Earthdata credentials.

you might need to `chmod 0600 ~/.netrc`
   
< icesat:short_name=ATL03:time_start='':time_end='':filename_filter='' >"""
    
    def __init__(self, short_name='ATL03', **kwargs):
        if short_name is not None:
            short_name = short_name.upper()
            if not short_name.startswith('ATL'):
                utils.echo_warning_msg('{} is not a valid icesat short_name, using ATL03'.format(short_name))
                short_name = 'ATL03'
                
        super().__init__(short_name=short_name, **kwargs)   
        self.data_format = 202
        self.src_srs = 'epsg:4326+3855'

## ==============================================
## SST from EarthData
##
## This module allows us to use SST data in dlim/waffles
##
## ==============================================
class MUR_SST(EarthData):
    """Access SST data.

< mur_sst:time_start='':time_end='':filename_filter='' >"""
    
    def __init__(self, **kwargs):
        super().__init__(short_name='MUR-JPL-L4-GLOB-v4.1', **kwargs)   
        
## ==============================================
## USIEI
## ==============================================
class USIEI(FetchModule):
    """US Interagency Elevation Inventory

No data is fetched with this module. Will list out query results from the USIEI.
Set 'want_geometry' to True to output a geojson formatted vector.

Fields:
OBJECTID ( type: esriFieldTypeOID, alias: OBJECTID )
Shape ( type: esriFieldTypeGeometry, alias: Shape )
ID ( type: esriFieldTypeInteger, alias: ID )
Title ( type: esriFieldTypeString, alias: Title, length: 250 )
DataType ( type: esriFieldTypeString, alias: DataType, length: 20 )
Status ( type: esriFieldTypeString, alias: Status, length: 20 )
Links ( type: esriFieldTypeString, alias: Links, length: 8000 )
pointspacing ( type: esriFieldTypeString, alias: pointspacing, length: 50 )
verticalaccuracy ( type: esriFieldTypeString, alias: verticalaccuracy, length: 8000 )
horizontalaccuracy ( type: esriFieldTypeString, alias: horizontalaccuracy, length: 300 )
collectiondate ( type: esriFieldTypeString, alias: collectiondate, length: 200 )
collectionyear ( type: esriFieldTypeSmallInteger, alias: collectionyear )
InfoContact ( type: esriFieldTypeString, alias: InfoContact, length: 500 )
qualitylevel ( type: esriFieldTypeSmallInteger, alias: qualitylevel )
meets3dep ( type: esriFieldTypeString, alias: meets3dep, length: 50 )
reasons3dep ( type: esriFieldTypeString, alias: reasons3dep, length: 150 )
meets3dep_lpc ( type: esriFieldTypeString, alias: meets3dep_lpc, length: 255 )
productsavailable ( type: esriFieldTypeString, alias: productsavailable, length: 300 )
pointclasses ( type: esriFieldTypeString, alias: pointclasses, length: 1000 )
verticaldatum ( type: esriFieldTypeString, alias: verticaldatum, length: 300 )
horizontaldatum ( type: esriFieldTypeString, alias: horizontaldatum, length: 300 )
restrictions ( type: esriFieldTypeString, alias: restrictions, length: 20 )
leafOnOff ( type: esriFieldTypeString, alias: leafOnOff, length: 255 )
AlternateTitle ( type: esriFieldTypeString, alias: AlternateTitle, length: 500 )
notes ( type: esriFieldTypeString, alias: notes, length: 500 )
RecordOwner ( type: esriFieldTypeString, alias: RecordOwner, length: 100 )
ContractSpec ( type: esriFieldTypeString, alias: ContractSpec, length: 255 )
Shape_Length ( type: esriFieldTypeDouble, alias: Shape_Length )
Shape_Area ( type: esriFieldTypeDouble, alias: Shape_Area )
    
layers:
  0 - Lidar-Topobathy
  1 - Lidar-Bathy
  2 - Lidar-Topo
  3 - IfSAR/InSAR
  4 - Other Bathy

https://coast.noaa.gov/inventory/

< usiei:where=None:layer=0:want_geometry=False >"""
    
    def __init__(self, where='1=1', want_geometry=False, layer=0, **kwargs):
        super().__init__(name='usiei', **kwargs)
        self._usiei_api_url = 'https://coast.noaa.gov/arcgis/rest/services/USInteragencyElevationInventory/USIEIv2/MapServer'
        self._usiei_query_url = '{0}/{1}/query?'.format(self._usiei_api_url, layer)
        self.where = where
        self.want_geometry = want_geometry
        
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
            'returnGeometry':'True' if self.want_geometry else 'False',
        }
        _req = Fetch(self._usiei_query_url, verbose=self.verbose).fetch_req(params=_data)
        if _req is not None:
            if self.want_geometry:
                print(_req.text)
            else:
                features = _req.json()
                print(json.dumps(features['features'], indent=True))
            
        return(self)
    
## ==============================================
## wsf
## ==============================================
class WSF(FetchModule):
    """WSF from German Aerospace Service (DLR)

World Settlement Footprint (WSF) 2019

https://www.dlr.de/EN/Home/home_node.html
https://geoservice.dlr.de/web/services
    
< wsf >"""

    def __init__(self, where='', datatype=None, **kwargs):
        super().__init__(name='wsf', **kwargs)

        self._wsf_url = 'https://download.geoservice.dlr.de/WSF2019/files/'
        self.where = [where] if len(where) > 0 else []
        self.datatype = datatype
        self.FRED = FRED.FRED(name=self.name, verbose=self.verbose)

        self.headers = { 'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:89.0) Gecko/20100101 Firefox/89.0' }
        self.update_if_not_in_FRED()
        
    def update_if_not_in_FRED(self):
        self.FRED._open_ds()
        self.FRED._attribute_filter(["DataSource = '{}'".format(self.name)])
        if len(self.FRED.layer) == 0:
            self.FRED._close_ds()
            self.update()
            
        self.FRED._close_ds()

    def update(self):
        """Crawl the SWF database and update/generate the WSF reference vector."""
        
        self.FRED._open_ds(1)
        surveys = []
        page = Fetch(self._wsf_url, verbose=True).fetch_html()
        rows = page.xpath('//a[contains(@href, ".tif")]/@href')
        with tqdm(
                total=len(rows),
                desc='scanning WSF datasets',
                leave=self.verbose
        ) as pbar:
            for i, row in enumerate(rows):
                pbar.update(1)
                sid = row.split('.')[0]
                if sid == 'WSF2019_cog':
                    continue

                self.FRED._attribute_filter(["ID = '{}'".format(sid)])
                if self.FRED.layer is None or len(self.FRED.layer) == 0:
                    spat = row.split('.')[0].split('_')
                    x = int(spat[-2])
                    y = int(spat[-1])
                    this_region = regions.Region().from_list(
                        [x, x + 2, y, y + 2]
                    )
                    geom = this_region.export_as_geom()
                    if geom is not None:
                        surveys.append(
                            {
                                'Name': row.split('.')[0],
                                'ID': sid,
                                'Agency': 'DLR',
                                'Date': utils.this_date(),
                                'MetadataLink': row.split('.')[0] + '_stac.json',
                                'MetadataDate': utils.this_date(),
                                'DataLink': self._wsf_url + row,
                                'DataType': 'WSF',
                                'DataSource': 'WSF',
                                'HorizontalDatum': 'epsg:4326',
                                'VerticalDatum': 'None',
                                'Info': '',
                                'geom': geom,
                            }
                        )

        self.FRED._add_surveys(surveys)
        self.FRED._close_ds()

    def run(self):
        '''Run the WSF fetching module'''

        if self.datatype is not None:
            self.where.append("DataType = '{}'".format(self.datatype))

        for surv in FRED._filter_FRED(self):
            for i in surv['DataLink'].split(','):
                self.results.append(
                    [i, i.split('/')[-1].split('?')[0], surv['DataType']]
                )
                
        return(self)

## ==============================================
## hydrolakes
## ==============================================
class HydroLakes(FetchModule):
    """HydroLakes vector and derived elevations
    
HydroLAKES aims to provide the shoreline polygons of all global lakes with a surface area 
of at least 10 ha. HydroLAKES has been developed using a suite of auxiliary data sources of 
lake polygons and gridded lake surface areas. All lakes are co-registered to the global 
river network of the HydroSHEDS database via their lake pour points. The global coverage of 
HydroLAKES encompasses 1.4 million individual lakes or reservoirs representing a total 
surface area of 2.67 million km², a total shoreline length of 7.2 million km, and a total 
storage volume of 181,900 km³.

hydrolakes:
https://wp.geog.mcgill.ca/hydrolab/data/
https://www.hydrosheds.org/products/hydrolakes
https://data.hydrosheds.org/file/hydrolakes/HydroLAKES_polys_v10_shp.zip

globathy:
https://springernature.figshare.com/collections/GLOBathy_the_Global_Lakes_Bathymetry_Dataset/5243309

add globathy output with the 'want_globathy' flag set to True

< hydrolakes:want_globathy=False >"""
    
    def __init__(self, where='1=1', want_globathy=False, **kwargs):
        super().__init__(name='hydrolakes', **kwargs)
        self.want_globathy = want_globathy
        self._hydrolakes_prods = 'https://www.hydrosheds.org/products/hydrolakes'
        self._hydrolakes_poly_zip = 'https://data.hydrosheds.org/file/hydrolakes/HydroLAKES_polys_v10_shp.zip'
        self._hydrolakes_gdb_zip = 'https://data.hydrosheds.org/file/hydrolakes/HydroLAKES_polys_v10_gdb.zip'
        self._globathy_url = 'https://springernature.figshare.com/ndownloader/files/28919991'
        self.where = [where] if len(where) > 0 else []

    def run(self):
        self.results.append(
            [self._hydrolakes_poly_zip, self._hydrolakes_poly_zip.split('/')[-1], 'hydrolakes']
        )

        if self.want_globathy:
            self.results.append(
                [self._globathy_url, 'globathy_parameters.zip', 'globathy']
            )
        return(self)
        
class ShallowBathyEverywhere(FetchModule):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.sbe_dem_url - 'https://shallowbathymetryeverywhere.com/data/dem/'

class HttpDataset(FetchModule):
    """fetch an http file"""
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def run(self):
        self.results.append([self.params['mod'], os.path.basename(self.params['mod']), 'https'])
        return(self)
        
## ==============================================
## Fetches Module Parser
## ==============================================
class FetchesFactory(factory.CUDEMFactory):
    """Acquire a fetches module."""
    
    _modules = {
        'gmrt': {'call': GMRT},
        'mar_grav': {'call': MarGrav},
        'srtm_plus': {'call': SRTMPlus},
        'charts': {'call': NauticalCharts}, # 'Charts' isn't working! :(
	    'digital_coast': {'call': DAV},
        'SLR': {'call': SLR},
        'CoNED': {'call': CoNED},
        'CUDEM': {'call': CUDEM},
        'multibeam': {'call': Multibeam}, # MBDB isn't working!
        'gebco': {'call': GEBCO},
        'mgds': {'call': MGDS},
        'trackline': {'call': Trackline},
        'ehydro': {'call': eHydro},
        'ngs': {'call': NGS},
        'hydronos': {'call': HydroNOS},
        'ncei_thredds': {'call': NCEIThreddsCatalog},
        'tnm': {'call': TheNationalMap},
        'emodnet': {'call': EMODNet},
        'chs': {'call': CHS}, # chs isn't working! awaiting IT assistance from CA
        'hrdem': {'call': HRDEM},
        'arcticdem': {'call': ArcticDEM},
        'bluetopo': {'call': BlueTopo},
        'osm': {'call': OpenStreetMap},
        'copernicus': {'call': CopernicusDEM},
        'fabdem': {'call': FABDEM},
        'nasadem': {'call': NASADEM},
        'tides': {'call': Tides},
        'vdatum': {'call': VDATUM},
        'buoys': {'call': BUOYS},
        'earthdata': {'call': EarthData},
        'icesat': {'call': IceSat},
        'mur_sst': {'call': MUR_SST},
        'usiei': {'call': USIEI},
        'wsf': {'call': WSF},
        'hydrolakes': {'call': HydroLakes},
        'https': {'call': HttpDataset},
        'bing_bfp': {'call': BingBFP},
        'waterservices': {'call': WaterServices},
    }

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
    
## ==============================================
## Command-line Interface (CLI)
## $ fetches
##
## fetches cli
## ==============================================
fetches_usage = """{cmd} ({version}): Fetches; Fetch and process remote elevation data

usage: {cmd} [ -hlqzHR [ args ] ] MODULE ...

Options:
  -R, --region\t\tRestrict processing to the desired REGION 
\t\t\tWhere a REGION is xmin/xmax/ymin/ymax[/zmin/zmax[/wmin/wmax]]
\t\t\tUse '-' to indicate no bounding range; e.g. -R -/-/-/-/-10/10/1/-
\t\t\tOR an OGR-compatible vector file with regional polygons. 
\t\t\tWhere the REGION is /path/to/vector[:zmin/zmax[/wmin/wmax]].
\t\t\tIf a vector file is supplied, will use each region found therein.
  -H, --threads\t\tSet the number of threads (1)
  -l, --list\t\tReturn a list of fetch URLs in the given region.
  -z, --no_check_size\tDon't check the size of remote data if local data exists.
  -q, --quiet\t\tLower the verbosity to a quiet

  --modules\t\tDisplay the module descriptions and usage
  --help\t\tPrint the usage text
  --version\t\tPrint the version information

Supported FETCHES modules (see fetches --modules <module-name> for more info): 
  {f_formats}
""".format(cmd=os.path.basename(sys.argv[0]), 
           version=cudem.__version__,
           f_formats=utils._cudem_module_short_desc(FetchesFactory._modules))

def fetches_cli(argv = sys.argv):
    """run fetches from command-line

See `fetches_cli_usage` for full cli options.
    """

    i_regions = []
    these_regions = []
    mods = []
    mod_opts = {}
    want_list = False
    want_verbose = True
    stop_threads = False
    check_size = True
    num_threads = 1
    
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
        elif arg == '-threads' or arg == '--threads' or arg == '-H':
            num_threads = utils.int_or(argv[i + 1], 1)
            i = i + 1
        elif arg == '--list' or arg == '-l':
            want_list = True
        elif arg == '--no_check_size' or arg == '-z':
            check_size = False
        elif arg == '--quiet' or arg == '-q':
            want_verbose = False
        elif arg == '--help' or arg == '-h':
            sys.stderr.write(fetches_usage)
            sys.exit(1)
        elif arg == '--version' or arg == '-v':
            print('{}, version {}'.format(
                os.path.basename(sys.argv[0]), cudem.__version__)
                  )
            sys.exit(1)
        elif arg == '--modules' or arg == '-m':
            utils.echo_modules(FetchesFactory._modules, None if i+1 >= len(argv) else sys.argv[i+1])
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

    these_regions = regions.parse_cli_region(i_regions, want_verbose)
    if not these_regions:
        these_regions = [regions.Region().from_string('-R-180/180/-90/90')]
        
    for rn, this_region in enumerate(these_regions):
        if stop_threads:
            return
        
        x_fs = [FetchesFactory(mod=mod, src_region=this_region, verbose=want_verbose)._acquire_module() for mod in mods]
        for x_f in x_fs:
            if x_f is None:
                continue
            
            if want_verbose:
                utils.echo_msg(
                    'running fetch module {} on region {}...'.format(
                        x_f.name, this_region.format('str')
                    )
                )

            try:
                x_f.run()
            except (KeyboardInterrupt, SystemExit):
                utils.echo_error_msg('user breakage...please wait while fetches exits.')
                sys.exit(-1)
                
            if want_verbose:
                utils.echo_msg(
                    'found {} data files.'.format(len(x_f.results))
                )
                
            if len(x_f.results) == 0:
                break
            
            if want_list:
                for result in x_f.results:
                    print(result[0])
            else:
                try:
                    fr = fetch_results(x_f, n_threads=num_threads, check_size=check_size)
                    fr.daemon = True
                
                    fr.start()
                    fr.join()                
                except (KeyboardInterrupt, SystemExit):
                    utils.echo_error_msg('user breakage...please wait while fetches exits.')
                    x_f.status = -1
                    stop_threads = True
                    while not fr.fetch_q.empty():
                        try:
                            fr.fetch_q.get(False)
                        except Empty:
                            continue
                        
                        fr.fetch_q.task_done()                        
### End
