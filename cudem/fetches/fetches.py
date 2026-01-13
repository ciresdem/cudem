### fetches.py
##
## Copyright (c) 2010 - 2026 Regents of the University of Colorado
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
### Commentary:
##
## Fetches: The CUDEM Data Acquisition Engine.
##
## Fetches provides a unified interface for discovering and downloading geospatial
## data from a wide variety of public repositories and APIs. It automates the 
## complex logic of querying remote services, parsing metadata, and managing 
## concurrent downloads.
##
##   * Extensive Data Support:
##      - Integrates with 40+ data sources including NOAA (NOS, NCEI, OCM), 
##        USGS (TNM, 3DEP), USACE (eHydro), NASA (EarthData), and OSM.
##      - Supports diverse formats: Bathymetry (BAG, XYZ), Topography (DEM, LiDAR),
##        Vectors (SHP, PBF), and Oceanographic data (Tides, Buoys).
##
##   * Intelligent Querying:
##      - Spatially filters datasets using bounding boxes or vector polygons.
##      - Parses ISO XML and HTML metadata to extract valid download links 
##        and spatial extents.
##      - Handles authentication (Netrc, EarthData tokens) transparently.
##
##   * Robust Downloading:
##      - Multi-threaded download engine for high-throughput data retrieval.
##      - Built-in resiliency with automatic retries, timeout handling, and 
##        resume capability (via Range headers).
##      - "List mode" to generate URL lists without downloading.
##
##   * Modular Architecture:
##      - Uses a Factory pattern to easily extend support for new datasets.
##      - Each module (e.g., 'srtm_plus', 'multibeam') encapsulates the specific
##        logic for querying its respective API.
##
## Usage:
##   CLI: fetches -R <region> <module> [options] (e.g., fetches -R -90/-89/29/30 srtm_plus)
##   API: f = FetchesFactory(mod='gmrt', src_region=region).acquire()
##        f.run()
##
### Code:

import os
import sys
import time
import base64
import threading
import queue
import netrc
import argparse
import urllib.request
import urllib.parse
from urllib.error import HTTPError, URLError
from urllib.request import Request, build_opener, HTTPCookieProcessor
from typing import List, Dict, Optional, Union, Any, Tuple

import requests
import lxml.etree
import lxml.html as lh

from cudem import utils
from cudem import regions
from cudem import gdalfun
from cudem import factory
from cudem.fetches import __version__ as __cudem_version__
from . import __version__

## ==============================================
## Constants & Configuration
## ==============================================
CUDEM_USER_AGENT = f'Fetches/{__version__} CUDEM/{__cudem_version__}'
DEFAULT_USER_AGENT = 'Mozilla/5.0 (X11; Linux x86_64; rv:146.0) Gecko/20100101 Firefox/146.0'
R_HEADERS = {'User-Agent': DEFAULT_USER_AGENT}

NAMESPACES = {
    'gmd': 'http://www.isotc211.org/2005/gmd', 
    'gmi': 'http://www.isotc211.org/2005/gmi', 
    'gco': 'http://www.isotc211.org/2005/gco',
    'gml': 'http://www.isotc211.org/2005/gml',
    'th': 'http://www.unidata.ucar.edu/namespaces/thredds/InvCatalog/v1.0',
    'wms': 'http://www.opengis.net/wms',
}

## ==============================================
## Helper Functions
## ==============================================
def fetches_callback(r: List[Any]):
    """Default callback for fetches processes.
    r: [url, local-fn, data-type, fetch-status-or-error-code]
    """
    
    pass


def urlencode_(opts: Dict) -> str:
    """Encode `opts` for use in a URL."""
    
    return urllib.parse.urlencode(opts)


def urlencode(opts: Dict, doseq: bool = True) -> str:
    """Encode `opts` for use in a URL.
    
    Args:
        opts: Dictionary of query parameters.
        doseq: If True, lists in values are encoded as separate parameters 
               (e.g., {'a': [1, 2]} -> 'a=1&a=2').
    """
    return urllib.parse.urlencode(opts, doseq=doseq)


def xml2py(node) -> Optional[Dict]:
    """Parse an xml file into a python dictionary."""
    
    texts = {}
    if node is None:
        return None

    for child in list(node):
        child_key = lxml.etree.QName(child).localname
        if 'name' in child.attrib:
            child_key = child.attrib['name']
        
        href = child.attrib.get('{http://www.w3.org/1999/xlink}href')
        
        if child.text is None or child.text.strip() == '':
            if href is not None:
                if child_key in texts:
                    texts[child_key].append(href)
                else:
                    texts[child_key] = [href]
            else:
                if child_key in texts:
                    ck = xml2py(child)
                    if ck:
                        first_key = list(ck.keys())[0]
                        texts[child_key][first_key].update(ck[first_key])
                else:
                    texts[child_key] = xml2py(child)
        else:
            if child_key in texts:
                texts[child_key].append(child.text)
            else:
                texts[child_key] = [child.text]
                
    return texts


def get_userpass(authenticator_url: str) -> Tuple[Optional[str], Optional[str]]:
    """Retrieve username and password from netrc for a given URL."""
    
    try:
        info = netrc.netrc()
        username, _, password = info.authenticators(urllib.parse.urlparse(authenticator_url).hostname)
    except Exception as e:
        if 'No such file' not in str(e):
            utils.echo_error_msg(f'Failed to parse netrc: {e}')
        username = None
        password = None

    return username, password


def get_credentials(url: str, authenticator_url: str = 'https://urs.earthdata.nasa.gov') -> Optional[str]:
    """Get user credentials from .netrc or prompt for input. 
    Used for EarthData.
    """
    
    credentials = None
    errprefix = ''
    
    username, password = get_userpass(authenticator_url)

    while not credentials:
        if not username:
            username = utils.get_username()
            password = utils.get_password()
            
        cred_str = f'{username}:{password}'
        credentials = base64.b64encode(cred_str.encode('ascii')).decode('ascii')

        if url:
            try:
                req = Request(url)
                req.add_header('Authorization', f'Basic {credentials}')
                opener = build_opener(HTTPCookieProcessor())
                opener.open(req)
            except HTTPError:
                print(f'{errprefix}Incorrect username or password')
                errprefix = ''
                credentials = None
                username = None
                password = None

    return credentials


## ==============================================
## ISO XML Parser
## ==============================================
class iso_xml:
    def __init__(self, xml_url: str, timeout: int = 2, read_timeout: int = 10):
        self.url = xml_url
        self.timeout = timeout
        self.read_timeout = read_timeout
        self.namespaces = NAMESPACES
        self.xml_doc = self._fetch()

        
    def _fetch(self):
        return Fetch(self.url).fetch_xml(timeout=self.timeout, read_timeout=self.read_timeout)

    
    def title(self) -> str:
        t = self.xml_doc.find(
            './/gmd:MD_DataIdentification/gmd:citation/gmd:CI_Citation/gmd:title/gco:CharacterString',
            namespaces=self.namespaces
        )
        return t.text if t is not None else 'Unknown'

    
    def bounds(self, geom: bool = True):
        wl = self.xml_doc.find('.//gmd:westBoundLongitude/gco:Decimal', namespaces=self.namespaces)
        el = self.xml_doc.find('.//gmd:eastBoundLongitude/gco:Decimal', namespaces=self.namespaces)
        sl = self.xml_doc.find('.//gmd:southBoundLatitude/gco:Decimal', namespaces=self.namespaces)
        nl = self.xml_doc.find('.//gmd:northBoundLatitude/gco:Decimal', namespaces=self.namespaces)
        
        if all(x is not None for x in [wl, el, sl, nl]):
            coords = [float(wl.text), float(el.text), float(sl.text), float(nl.text)]
            if geom:
                return regions.Region().from_list(coords).export_as_geom()
            else:
                return coords          
        return None

    
    def polygon(self, geom: bool = True):
        opoly = []
        polygon = self.xml_doc.find('.//{*}Polygon', namespaces=self.namespaces)
        if polygon is not None:
            nodes = polygon.findall('.//{*}pos', namespaces=self.namespaces)
            for node in nodes:
                opoly.append([float(x) for x in node.text.split()])
                
            ## Close polygon if open
            if opoly and (opoly[0][0] != opoly[-1][0] or opoly[0][1] != opoly[-1][1]):
                opoly.append(opoly[0])
            
            if geom:
                return gdalfun.wkt2geom(regions.create_wkt_polygon(opoly))
            return opoly
        return None

    
    def date(self) -> str:
        dt = self.xml_doc.find('.//gmd:date/gco:Date', namespaces=self.namespaces)
        if dt is None:
            dt = self.xml_doc.find(
                './/gmd:MD_DataIdentification/gmd:citation/gmd:CI_Citation/gmd:date/gmd:CI_Date/gmd:date/gco:Date',
                namespaces=self.namespaces
            )
        return dt.text[:4] if dt is not None else '0000'

    
    def xml_date(self) -> str:
        mddate = self.xml_doc.find('.//gmd:dateStamp/gco:DateTime', namespaces=self.namespaces)
        return utils.this_date() if mddate is None else mddate.text

    
    def reference_system(self) -> Tuple[Optional[str], Optional[str]]:
        ref_s = self.xml_doc.findall('.//gmd:MD_ReferenceSystem', namespaces=self.namespaces)
        if not ref_s:
            return None, None
        
        h_epsg = ref_s[0].find('.//gmd:code/gco:CharacterString', namespaces=self.namespaces)
        h_epsg = h_epsg.text.split(':')[-1] if h_epsg is not None else None
        
        v_epsg = None
        if len(ref_s) > 1:
            v_val = ref_s[1].find('.//gmd:code/gco:CharacterString', namespaces=self.namespaces)
            if v_val is not None:
                v_epsg = v_val.text.split(':')[-1]
                
        return h_epsg, v_epsg

    
    def abstract(self) -> str:
        try:
            abstract = self.xml_doc.find('.//gmd:abstract/gco:CharacterString', namespaces=self.namespaces)
            return '' if abstract is None else abstract.text
        except:
            return ''

        
    def linkages(self) -> Optional[str]:
        linkage = self.xml_doc.find('.//{*}linkage/{*}URL', namespaces=self.namespaces)
        return linkage.text if linkage is not None else None

    
    def data_links(self) -> Dict[str, List[str]]:
        dd = {}        
        dfs = self.xml_doc.findall('.//gmd:MD_Format/gmd:name/gco:CharacterString', namespaces=self.namespaces)
        dus = self.xml_doc.findall('.//gmd:onLine/gmd:CI_OnlineResource/gmd:linkage/gmd:URL', namespaces=self.namespaces)

        if dfs:
            for i, j in enumerate(dfs):
                if i < len(dus):
                    key = j.text
                    val = dus[i].text
                    if key in dd:
                        dd[key].append(val)
                    else:
                        dd[key] = [val]
        return dd

    
## ==============================================
## Fetch 
## ==============================================
class Fetch:
    """Fetch class to fetch ftp/http data files"""
    
    def __init__(
            self,
            url: str = None,
            callback = fetches_callback,
            verbose: bool = None,
            headers: Dict = R_HEADERS,
            verify: bool = True,
            allow_redirects: bool = True
    ):
        self.url = url
        self.callback = callback
        self.verbose = verbose
        self.headers = headers
        self.verify = verify
        self.allow_redirects = allow_redirects


    def fetch_req(self, params=None, json=None, tries=5, timeout=None, read_timeout=None) -> Optional[requests.Response]:
        """Fetch src_url and return the requests object (iterative retry)."""
        
        req = None
        current_timeout = timeout
        current_read_timeout = read_timeout

        for attempt in range(tries):
            try:
                ## Calculate timeouts for this attempt
                to_tuple = (
                    current_timeout if current_timeout else None,
                    current_read_timeout if current_read_timeout else None
                )

                req = requests.get(
                    self.url,
                    stream=True,
                    params=params,
                    json=json,
                    timeout=to_tuple,
                    verify=self.verify,
                    allow_redirects=self.allow_redirects,
                    headers=self.headers
                )
                
                ## Check status codes
                if req.status_code == 504: # Gateway Timeout
                    time.sleep(2)
                    ## Increase timeouts next loop
                    if current_timeout: current_timeout += 1
                    if current_read_timeout: current_read_timeout += 10
                    continue

                elif req.status_code == 416: # Range Not Satisfiable
                    ## If range fails, try fetching whole file
                    if 'Range' in self.headers:
                        del self.headers['Range']
                        continue
                
                elif 200 <= req.status_code <= 299:
                    return req
                
                else:
                    if self.verbose:
                        utils.echo_error_msg(f'Request from {req.url} returned {req.status_code}')
                    return req

            except Exception as e:
                utils.echo_warning_msg(f"Attempt {attempt + 1}/{tries} failed: {e}")
                if current_timeout: current_timeout *= 2
                if current_read_timeout: current_read_timeout *= 2
                time.sleep(1)

        utils.echo_error_msg(f'Max-tries exhausted {self.url}')
        raise ConnectionError('Maximum attempts at connecting have failed.')
    
    
    def fetch_html(self, timeout=2):
        """Fetch src_url and return it as an HTML object."""
        
        req = self.fetch_req(timeout=timeout)
        if req:
            return lh.document_fromstring(req.text)
        return None

    
    def fetch_xml(self, timeout=2, read_timeout=10):
        """Fetch src_url and return it as an XML object."""
        
        try:
            req = self.fetch_req(timeout=timeout, read_timeout=read_timeout)
            results = lxml.etree.fromstring(req.text.encode('utf-8'))
        except Exception as e:
            ## Fallback empty XML
            results = lxml.etree.fromstring(
                '<?xml version="1.0"?><!DOCTYPE _[<!ELEMENT _ EMPTY>]><_/>'.encode('utf-8')
            )
        return results

    
    def fetch_file(
            self,
            dst_fn: str,
            params=None,
            datatype=None,
            overwrite=False,
            timeout=30,
            read_timeout=None,
            tries=5,
            check_size=True
    ) -> int:
        """Fetch src_url and save to dst_fn."""
        
        status = 0
        
        ## Ensure directory exists
        dst_dir = os.path.abspath(os.path.dirname(dst_fn))
        if not os.path.exists(dst_dir):
            try:
                os.makedirs(dst_dir)
            except Exception as e:
                utils.echo_error_msg(e)

        ## Handle Ranges for resuming
        if check_size and 'Range' in self.headers:
            del self.headers['Range']
        
        if not overwrite and os.path.exists(dst_fn):
            if not check_size:
                ## File exists and we don't care about size -> Success
                return 0
            else:
                dst_fn_size = os.stat(dst_fn).st_size
                self.headers['Range'] = f'bytes={dst_fn_size}-'

        try:
            for attempt in range(tries):
                try:
                    with requests.get(
                            self.url, stream=True, params=params, headers=self.headers,
                            timeout=(timeout, read_timeout), verify=self.verify
                    ) as req:

                        ## Determine content length
                        req_h = req.headers
                        content_length = 0
                        if 'Content-Range' in req_h:
                            content_length = int(req_h['Content-Range'].split('/')[-1])
                        elif 'Content-Length' in req_h:
                            content_length = int(req_h['Content-Length'])
                        else:
                            content_length = int(req_h.get('content-length', 0))

                        ## Check if file is already complete
                        if not overwrite and check_size and os.path.exists(dst_fn):
                            if content_length == os.path.getsize(dst_fn) or content_length == 0:
                                return 0 # Success, already done

                        # Handle Status Codes
                        if req.status_code == 416: # Range error
                            overwrite = True # Force overwrite next try
                            raise FileExistsError(f'{dst_fn} exists, invalid Range requested.')
                        
                        elif req.status_code == 401: # Auth error
                            ## Earthdata hack: they redirect 401s sometimes
                            if self.url == req.url:
                                raise UnboundLocalError('Incorrect Authentication')
                            ## Recursion for redirect url
                            return Fetch(
                                url=req.url, headers=self.headers, verbose=self.verbose
                            ).fetch_file(dst_fn, params, datatype, overwrite, timeout, read_timeout)

                        elif req.status_code in [429, 504]: # Rate limit / Gateway
                            if attempt < tries - 1:
                                time.sleep(10)
                                continue
                            else:
                                raise ConnectionError(f'{req.url}: {req.status_code}')

                        elif req.status_code in [200, 206]: # OK or Partial Content
                            mode = 'ab' if req.status_code == 206 else 'wb'
                            with open(dst_fn, mode) as local_file:
                                with utils.ccp(
                                        desc=self.url,
                                        total=content_length,
                                        leave=self.verbose,
                                        unit='iB',
                                        unit_scale=True,
                                ) as pbar:
                                    for chunk in req.iter_content(chunk_size=8196):
                                        if chunk:
                                            pbar.update(len(chunk))
                                            local_file.write(chunk)
                                            local_file.flush()
                            
                            ## Verify Size
                            if check_size and (content_length != 0):
                                if content_length != os.stat(dst_fn).st_size:
                                    raise UnboundLocalError(f'Size mismatch! Expected {content_length}, got {os.stat(dst_fn).st_size}')
                            
                            return 0 # Success

                        else:
                            ## Unknown error code
                            if self.verbose:
                                utils.echo_error_msg(f'Server returned: {req.status_code} ({req.url})')
                            return -1

                except (requests.exceptions.RequestException, UnboundLocalError) as e:
                    if attempt < tries - 1:
                        if self.verbose:
                            utils.echo_warning_msg(f'Exception: {e}. Retrying ({tries - attempt - 1} left)...')
                        time.sleep(2)
                        continue
                    else:
                        raise e

        except FileExistsError:
            return 0
        except Exception as e:
            if self.verbose:
                utils.echo_error_msg(str(e))
            status = -1

        ## Final check
        if check_size and (not os.path.exists(dst_fn) or os.stat(dst_fn).st_size == 0):
            status = -1
            
        return status

    def fetch_ftp_file(self, dst_fn, params=None, datatype=None, overwrite=False):
        """Fetch an ftp file via urllib."""
        
        status = 0
        
        if self.verbose:
            utils.echo_msg(f'Fetching remote ftp file: {self.url[:20]}...')
            
        if not os.path.exists(os.path.dirname(dst_fn)):
            try:
                os.makedirs(os.path.dirname(dst_fn))
            except:
                pass
            
        try:
            with urllib.request.urlopen(self.url) as f:
                with open(dst_fn, 'wb') as local_file:
                     local_file.write(f.read())
            
            if self.verbose:
                utils.echo_msg(f'Fetched remote ftp file: {os.path.basename(self.url)}.')
        except Exception as e:
            utils.echo_error_msg(e)
            status = -1

        return status


## ==============================================
## Threading & Queues
## ==============================================
def fetch_queue(q: queue.Queue, stop_event: threading.Event, c: bool = True):
    """Worker for the fetch queue.
    q items: [remote_data_url, local_data_path, data-type, fetches-module, attempts, results-list]    
    """
    
    ## Modules that bypass SSL verification
    no_verify = ['mar_grav', 'srtm_plus']

    while not stop_event.is_set():
        fetch_args = q.get()
        url = fetch_args[0]
        local_path = fetch_args[1]
        data_type = fetch_args[2]
        module = fetch_args[3]
        retries = fetch_args[4]
        results_list = fetch_args[5]

        if stop_event.is_set():
            q.task_done()
            continue
        
        ## Ensure dir exists
        if not os.path.exists(os.path.dirname(local_path)):
            try:
                os.makedirs(os.path.dirname(local_path))
            except: pass

        parsed_url = urllib.parse.urlparse(url)
        
        try:
            if parsed_url.scheme == 'ftp':
                status = Fetch(
                    url=url,
                    callback=module.callback,
                    verbose=module.verbose,
                    headers=module.headers
                ).fetch_ftp_file(local_path)
            else:
                verify_ssl = False if module.name in no_verify else True
                status = Fetch(
                    url=url,
                    callback=module.callback,
                    verbose=module.verbose,
                    headers=module.headers,
                    verify=verify_ssl
                ).fetch_file(local_path, check_size=c)

            ## Record result
            fetch_results_entry = [url, local_path, data_type, status]
            results_list.append(fetch_results_entry)

            if callable(module.callback):
                module.callback(fetch_results_entry)
        
        except Exception as e:
            if retries > 0:
                fetch_args[4] -= 1
                utils.remove_glob(local_path)
                q.put(fetch_args) # Put back in queue
            else:
                utils.echo_error_msg(f'Fetch of {url} failed...')
                module.status = -1
                fetch_results_entry = [url, local_path, data_type, str(e)]
                results_list.append(fetch_results_entry)
                
                if callable(module.callback):
                    module.callback(fetch_results_entry)

        q.task_done()

        
class fetch_results(threading.Thread):
    """Threaded fetch runner."""
    
    def __init__(self, mod, check_size=True, n_threads=3, attempts=5):
        threading.Thread.__init__(self)
        self.fetch_q = queue.Queue()
        self.stop_event = threading.Event()
        self.mod = mod
        self.check_size = check_size
        self.n_threads = n_threads
        self.attempts = attempts
        self.results = []
        
        ## Ensure module has run to generate URLs
        if len(self.mod.results) == 0:
            self.mod.run()

            
    def run(self):
        ## Start workers
        for _ in range(self.n_threads):
            t = threading.Thread(
                target=fetch_queue,
                args=(self.fetch_q, self.stop_event, self.check_size)
            )
            t.daemon = True
            t.start()

        ## Populate queue
        ## Queue item: [url, path, type, module, retries, results_ptr]
        for row in self.mod.results:
            if row['dst_fn']:
                self.fetch_q.put([
                    row['url'],
                    os.path.join(self.mod._outdir, row['dst_fn']),
                    row['data_type'],
                    self.mod,
                    self.attempts,
                    self.results
                ])

        ## Wait
        while not self.fetch_q.empty() and not self.stop_event.is_set():
             time.sleep(0.1)
        
        if not self.stop_event.is_set():
            self.fetch_q.join()
            
        #self.fetch_q.join()

    def stop(self):
        """Stop all threads"""
        
        self.stop_event.set()
        
## ==============================================
## Fetch Modules (Base & Implementations)
## ==============================================
class FetchModule:
    """Base class for all fetch modules."""
    
    def __init__(
            self,
            src_region=None,
            callback=fetches_callback,
            verbose=True,
            outdir=None,
            name='fetches',
            min_year=None,
            max_year=None,
            params={}
    ):
        self.region = src_region
        self.callback = callback
        self.verbose = verbose
        self.outdir = outdir
        self.params = params
        self.status = 0
        self.results = []
        self.name = name
        self.min_year = utils.int_or(min_year)
        self.max_year = utils.int_or(max_year)
        self.headers = R_HEADERS.copy()

        if self.outdir is None:
            self._outdir = os.path.join(os.getcwd(), self.name)
        else:
            self._outdir = os.path.join(self.outdir, self.name)

        ## For dlim support, we can check these variables for
        ## to do the proper processing. Set these to their correct
        ## values in the sub-class
        self.data_format = None
        self.src_srs = None
        self.title = None
        self.source = None
        self.date = None
        self.data_type = None
        self.resolution = None
        self.hdatum = None
        self.vdatum = None
        self.url = None
            
        ## Default to whole world if region is invalid/missing
        ## Set a generic region of the entire world in WGS84 if no region
        ## was specified or if its an invalid region...this will result in quite
        ## a lot of downloading on global datasets, so be careful with this.
        if self.region is None or not self.region.valid_p():
            self.region = regions.Region().from_list([-180, 180, -90, 90])

            
    def run(self):
        raise NotImplementedError

    
    def fetch_entry(self, entry, check_size=True, retries=5):
        try:
            status = Fetch(
                url=entry['url'],
                verbose=self.verbose,
                headers=self.headers,
            ).fetch_file(
                os.path.join(self._outdir, entry['dst_fn']),
                check_size=check_size,
                tries=retries
            )
        except:
            status = -1
        return status

    
    def fetch_results(self):
        """fetch the gathered `results` from the sub-class"""
        
        for entry in self.results:
            status = self.fetch(entry)

            
    def fill_results(self, entry):
        """fill self.results with the fetch module entry"""
        
        self.results.append(
            {'url': entry[0], 'dst_fn': entry[1], 'data_type': entry[2]}
        )

        
    def add_entry_to_results(self, url, dst_fn, data_type, **kwargs):
        entry = {'url': url, 'dst_fn': dst_fn, 'data_type': data_type}
        entry.update(kwargs)
        self.results.append(entry)             

        
class HttpDataset(FetchModule):
    """Fetch an http file directly."""
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        
    def run(self):
        if 'mod' in self.params:
            self.add_entry_to_results(
                self.params['mod'],
                os.path.basename(self.params['mod']),
                'https'
            )


class Nominatim(FetchModule):
    """Fetch coordinates from OpenStreetMap's Nominatim service.
    
    Nominatim requires a valid User-Agent identifying the application.
    """
    
    def __init__(self, q='boulder', **kwargs):
        super().__init__(name='nominatim', **kwargs)
        self.q = q
        self._nom_url = 'https://nominatim.openstreetmap.org/search'
        
        ## Nominatim usage policy requires a custom User-Agent and Referer.
        self.headers = {
            'User-Agent': 'CUDEM/Fetches 1.0 (cudem.colorado.edu)',
            'Referer': 'https://cudem.colorado.edu'
        }

        
    def run(self):
        if utils.str_or(self.q) is not None:
            ## Construct parameters using the module's urlencode helper
            params = {
                'q': self.q,
                'format': 'jsonv2',
                'limit': 1,
                'addressdetails': 1
            }
            query_str = urlencode(params)
            q_url = f'{self._nom_url}?{query_str}'

            _req = Fetch(q_url, headers=self.headers, verbose=self.verbose).fetch_req()
            
            if _req is not None and _req.status_code == 200:
                try:
                    results = _req.json()
                    if results and isinstance(results, list) and len(results) > 0:
                        ## Parse coordinates
                        x = utils.float_or(results[0].get("lon"))
                        y = utils.float_or(results[0].get("lat"))
                        
                        if self.verbose:
                            ## Print the display name found (helpful for debugging vague queries)
                            disp_name = results[0].get("display_name", "Unknown Location")
                            utils.echo_msg(f"Resolved '{self.q}' to: {disp_name}")
                            
                        ## Standard output for CLI piping: "lon, lat"
                        #print(f'{x}, {y}')
                        
                        ## Populate results list in case this is called programmatically
                        self.results.append({
                            'url': q_url,
                            'dst_fn': None,
                            'data_type': 'coords',
                            'metadata': results[0],
                            'x': x,
                            'y': y
                        })
                    else:
                        utils.echo_warning_msg(f"Nominatim: No results found for query '{self.q}'")
                except Exception as e:
                    utils.echo_error_msg(f"Nominatim parse error: {e}")
            else:
                status = _req.status_code if _req else "Connection Failed"
                utils.echo_error_msg(f"Nominatim request failed: {status}")                

                    
class GPSCoordinates(FetchModule):
    """Fetch various coordinates for places via gps-coordinates.net."""
    
    def __init__(self, q='boulder', **kwargs):
        super().__init__(name='gps_coordinates', **kwargs)
        self.q = q
        self.gpsc_api_url = "http://www.gps-coordinates.net/api/"

        
    def run(self):
        if utils.str_or(self.q) is not None:
            q_url = f'{self.gpsc_api_url}{self.q}'
            _req = Fetch(q_url, verbose=self.verbose).fetch_req()
            if _req is not None:
                results = _req.json()
                if results.get("responseCode") == '200':
                    x = utils.float_or(results["longitude"])
                    y = utils.float_or(results["latitude"])
                    print(f'{x}, {y}')
                else:
                    print(results)

                    
## ==============================================
## Fetches Factory
## ==============================================
class FetchesFactory(factory.CUDEMFactory):
    """Factory to acquire and initialize specific fetch modules."""
    
    from . import gmrt, margrav, srtmplus, synbath, charts, dav, multibeam, gebco, gedtm30, \
        mgds, trackline, ehydro, ngs, hydronos, nceithredds, etopo, tnm, emodnet, chs, \
        hrdem, mrdem, arcticdem, bluetopo, osm, copernicus, fabdem, nasadem, tides, vdatum, \
        buoys, earthdata, usiei, wsf, hydrolakes, bingbfp, waterservices, csb, cptcity, \
        wadnr, nswtb, cdse, gba
    
    _modules = {
        'https': {'call': HttpDataset},
        'gmrt': {'call': gmrt.GMRT},
        'mar_grav': {'call': margrav.MarGrav},
        'srtm_plus': {'call': srtmplus.SRTMPlus},
        'synbath': {'call': synbath.SynBath},
        'charts': {'call': charts.NauticalCharts},
        'digital_coast': {'call': dav.DAV},
        'SLR': {'call': dav.SLR},
        'CoNED': {'call': dav.CoNED},
        'CUDEM': {'call': dav.CUDEM},
        'multibeam': {'call': multibeam.Multibeam}, 
        'mbdb': {'call': multibeam.MBDB},
        'r2r': {'call': multibeam.R2R},
        'gebco': {'call': gebco.GEBCO},
        'gedtm30': {'call': gedtm30.GEDTM30},
        'mgds': {'call': mgds.MGDS},
        'trackline': {'call': trackline.Trackline},
        'ehydro': {'call': ehydro.eHydro},
        'ngs': {'call': ngs.NGS},
        'hydronos': {'call': hydronos.HydroNOS},
        'ncei_thredds': {'call': nceithredds.NCEIThreddsCatalog},
        'etopo': {'call': etopo.ETOPO},
        'tnm': {'call': tnm.TheNationalMap},
        'ned': {'call': tnm.NED},
        'ned1': {'call': tnm.NED1},
        'tnm_laz': {'call': tnm.TNM_LAZ},
        'emodnet': {'call': emodnet.EMODNet},
        'chs': {'call': chs.CHS},
        'hrdem': {'call': hrdem.HRDEM},
        'mrdem': {'call': mrdem.MRDEM},
        'arcticdem': {'call': arcticdem.ArcticDEM},
        'bluetopo': {'call': bluetopo.BlueTopo},
        'osm': {'call': osm.OpenStreetMap},
        'copernicus': {'call': copernicus.CopernicusDEM},
        'fabdem': {'call': fabdem.FABDEM},
        'nasadem': {'call': nasadem.NASADEM},
        'tides': {'call': tides.Tides},
        'vdatum': {'call': vdatum.VDATUM},
        'buoys': {'call': buoys.BUOYS},
        'earthdata': {'call': earthdata.EarthData},
        'icesat2': {'call': earthdata.IceSat2},
        'mur_sst': {'call': earthdata.MUR_SST},
        'swot': {'call': earthdata.SWOT},
        'usiei': {'call': usiei.USIEI},
        'wsf': {'call': wsf.WSF},
        'hydrolakes': {'call': hydrolakes.HydroLakes},
        'bingbfp': {'call': bingbfp.BingBFP},
        'waterservices': {'call': waterservices.WaterServices},
        'csb': {'call': csb.CSB},
        'cpt_city': {'call': cptcity.CPTCity},
        'gps_coordinates': {'call': GPSCoordinates},
        'nominatim': {'call': Nominatim},
        'wa_dnr': {'call': wadnr.WADNR},
        'nsw_tb': {'call': nswtb.NSW_TB},
        'sentinel2': {'call': cdse.Sentinel2},
        'gba': {'call': gba.GBA},
    }
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)


## ==============================================
## Command-line Interface (CLI)
## $ fetches
##
## fetches cli
## ==============================================
class PrintModulesAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        factory.echo_modules(FetchesFactory._modules, values, md=True if not values else False)
        sys.exit(0)

        
def fetches_cli():
    """Run fetches from command-line using argparse."""
    
    parser = argparse.ArgumentParser(
        description=f"%(prog)s ({__version__}): Fetch and process remote elevation data",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=f"""
Supported %(prog)s modules (see %(prog)s --modules <module-name> for more info): 
{factory.get_module_short_desc(FetchesFactory._modules)}
        
CUDEM home page: <http://cudem.colorado.edu>
        """

    )

    parser.add_argument(
        '-R', '--region', '--aoi',
        action='append',
        help=regions.region_help_msg()
    )
    parser.add_argument(
        '-H', '--threads',
        type=int,
        default=1,
        help="Set the number of threads (default: 1)"
    )
    parser.add_argument(
        '-A', '--attempts',
        type=int,
        default=5,
        help="Set the number of fetching attempts (default: 5)"
    )
    parser.add_argument(
        '-l', '--list',
        action='store_true',
        help="Return a list of fetch URLs in the given region."
    )
    parser.add_argument(
        '-z', '--no_check_size',
        action='store_true',
        help="Don't check the size of remote data if local data exists."
    )
    parser.add_argument(
        '-q', '--quiet',
        action='store_true',
        help="Lower the verbosity to a quiet"
    )
    parser.add_argument(
        '-m', '--modules',
        nargs='?',
        default=None,
        action=PrintModulesAction,
        help="Display the module descriptions and usage"
    )
    parser.add_argument(
        '-v', '--version',
        action='version',
        version=f'CUDEM {__cudem_version__} :: %(prog)s {__version__}'
    )
    parser.add_argument(
        'modules_to_run',
        nargs='+',
        help="The modules to run (e.g., srtm_plus, gmrt, etc.)"
    )

    ## Parse arguments
    args = parser.parse_args()

    ## Validate Positional Arguments
    if not args.modules_to_run:
        parser.print_help()
        utils.echo_error_msg('You must select at least one fetch module')
        sys.exit(-1)

    ## Process Flags
    want_verbose = not args.quiet
    check_size = not args.no_check_size

    ## Parse Regions
    ## If no region provided, default to world.
    if not args.region:
        these_regions = [regions.Region().from_string('-R-180/180/-90/90')]
    else:
        these_regions = regions.parse_cli_region(args.region, want_verbose)
        
    ## Execution Loop by region
    for this_region in these_regions:
        x_fs = [
            FetchesFactory(
                mod=mod,
                src_region=this_region,
                verbose=want_verbose
            )._acquire_module() for mod in args.modules_to_run
        ]
        
        for x_f in x_fs:
            if x_f is None:
                continue
            
            if want_verbose:
                utils.echo_msg(f'Running fetch module {x_f.name} on region {this_region.format("str")}...')

            try:
                x_f.run()
            except (KeyboardInterrupt, SystemExit):
                utils.echo_error_msg('User breakage... exiting.')
                sys.exit(-1)
                
            if want_verbose:
                utils.echo_msg(f'Found {len(x_f.results)} data files.')
                
            if not x_f.results:
                continue
            
            if args.list:
                for result in x_f.results:
                    print(result['url'])
            else:
                try:
                    fr = fetch_results(
                        x_f,
                        n_threads=args.threads,
                        check_size=check_size,
                        attempts=args.attempts
                    )
                    fr.daemon = True                
                    fr.start()
                    fr.join()         
                except (KeyboardInterrupt, SystemExit):
                    utils.echo_error_msg('User breakage... please wait while fetches exits.')
                    x_f.status = -1
                    ## Drain queue
                    while not fr.fetch_q.empty():
                        try:
                            fr.fetch_q.get(False)
                            fr.fetch_q.task_done()
                        except queue.Empty:
                            break
### End
