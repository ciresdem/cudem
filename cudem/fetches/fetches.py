### fetches.py
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
## Fetch elevation and related data from a variety of sources.
##
## Use CLI command 'fetches'
## or use FetchesFactory() to acquire and use a fetch module in python.
##
### TODO:
##
### Code:

import os
import time
import requests
import urllib
import lxml.etree
import lxml.html as lh
from tqdm import tqdm
import threading
try:
    import Queue as queue
except: import queue as queue

from cudem import utils
from cudem import regions
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
r_headers = {
    'User-Agent': 'Mozilla/5.0 (X11; Linux x86_64; rv:45.0) Gecko/20100101 Firefox/45.0'
}

## Namespaces for vaious XML files
namespaces = {
    'gmd': 'http://www.isotc211.org/2005/gmd', 
    'gmi': 'http://www.isotc211.org/2005/gmi', 
    'gco': 'http://www.isotc211.org/2005/gco',
    'gml': 'http://www.isotc211.org/2005/gml',
    }
thredds_namespaces = {
    'th': 'http://www.unidata.ucar.edu/namespaces/thredds/InvCatalog/v1.0',
}

# http = urllib3.PoolManager(
#     cert_reqs="CERT_REQUIRED",
#     ca_certs=certifi.where()
# )

## callback for use in fetches. currently, only the coastline and hydrolakes modules use
## fetches processes...Change this function to better handle failed fetches. `r` is the
## fetches results as a list: [url, local-fn, data-type, fetch-status-or-error-code]
## this code is run in fetches.fetch_results() after a successful or failed download.
def fetches_callback(r):
    pass


def urlencode(opts):
    """encode `opts` for use in a URL"""
    
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


def get_userpass(authenticator_url):
    try:
        info = netrc.netrc()
        username, account, password \
            = info.authenticators(urlparse(authenticator_url).hostname)
        errprefix = 'netrc error: '
    except Exception as e:
        if (not ('No such file' in str(e))):
            print('netrc error: {0}'.format(str(e)))
        username = None
        password = None

    return(username, password)


def get_credentials(url, authenticator_url='https://urs.earthdata.nasa.gov'):
    """Get user credentials from .netrc or prompt for input. 
    Used for EarthData.
    """
    
    credentials = None
    errprefix = ''
    try:
        info = netrc.netrc()
        username, account, password \
            = info.authenticators(urlparse(authenticator_url).hostname)
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
                req.add_header(
                    'Authorization', 'Basic {0}'.format(
                        credentials
                    )
                )
                opener = build_opener(HTTPCookieProcessor())
                opener.open(req)
            except HTTPError:
                print(errprefix + 'Incorrect username or password')
                errprefix = ''
                credentials = None
                username = None
                password = None

    return(credentials)


## a few convenience functions to parse common iso xml files
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
        t = self.xml_doc.find(
            ('.//gmd:MD_DataIdentification/gmd:citation/'
             'gmd:CI_Citation/gmd:title/gco:CharacterString'),
            namespaces=self.namespaces
        )
        return(t.text if t is not None else 'Unknown')
        
    def bounds(self, geom = True):
        wl = self.xml_doc.find(
            './/gmd:westBoundLongitude/gco:Decimal',
            namespaces=self.namespaces
        )
        el = self.xml_doc.find(
            './/gmd:eastBoundLongitude/gco:Decimal',
            namespaces=self.namespaces
        )
        sl = self.xml_doc.find(
            './/gmd:southBoundLatitude/gco:Decimal',
            namespaces=self.namespaces
        )
        nl = self.xml_doc.find(
            './/gmd:northBoundLatitude/gco:Decimal',
            namespaces=self.namespaces
        )
        if wl is not None \
           and el is not None \
           and sl is not None \
           and nl is not None:
            region = [float(wl.text), float(el.text),
                      float(sl.text), float(nl.text)]
            if geom:
                return(
                    regions.Region().from_list(
                        [float(wl.text), float(el.text),
                         float(sl.text), float(nl.text)]
                    ).export_as_geom()
                )
            else:
                return(region)            
        else:
            return(None)

        
    def polygon(self, geom=True):
        opoly = []
        polygon = self.xml_doc.find(
            './/{*}Polygon', namespaces=self.namespaces
        )
        if polygon is not None:
            nodes = polygon.findall('.//{*}pos', namespaces=self.namespaces)
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
        dt = self.xml_doc.find(
            './/gmd:date/gco:Date', namespaces=self.namespaces
        )
        if dt is None:
            dt = self.xml_doc.find(
                ('.//gmd:MD_DataIdentification/gmd:citation/'
                 'gmd:CI_Citation/gmd:date/gmd:CI_Date/gmd:date/gco:Date'),
                namespaces=self.namespaces
            )
            
        return(dt.text[:4] if dt is not None else '0000')

    
    def xml_date(self):
        mddate = self.xml_doc.find(
            './/gmd:dateStamp/gco:DateTime',
            namespaces=self.namespaces
        )
        
        return(utils.this_date() if mddate is None else mddate.text)

    
    def reference_system(self):
        ref_s = self.xml_doc.findall(
            './/gmd:MD_ReferenceSystem',
            namespaces=self.namespaces
        )
        if ref_s is None or len(ref_s) == 0:
            return(None, None)
        
        h_epsg = ref_s[0].find(
            './/gmd:code/gco:CharacterString',
            namespaces=self.namespaces
        )
        if h_epsg is not None:
            h_epsg = h_epsg.text.split(':')[-1]
        
        if len(ref_s) > 1:
            v_epsg = ref_s[1].find(
                './/gmd:code/gco:CharacterString',
                namespaces=self.namespaces
            )
            if v_epsg is not None:
                v_epsg = v_epsg.text.split(':')[-1]
                
        else:
            v_epsg = None
            
        return(h_epsg, v_epsg)

    
    def abstract(self):
        try:
            abstract = self.xml_doc.find(
                './/gmd:abstract/gco:CharacterString',
                namespaces=self.namespaces
            )
            abstract = '' if abstract is None else abstract.text
        except:
            abstract = ''
            
        return(abstract)

    
    def linkages(self):
        linkage = self.xml_doc.find(
            './/{*}linkage/{*}URL',
            namespaces=self.namespaces
        )
        if linkage is not None:
            linkage = linkage.text
        
        return(linkage)

    
    def data_links(self):
        dd = {}        
        dfs = self.xml_doc.findall(
            './/gmd:MD_Format/gmd:name/gco:CharacterString',
            namespaces=self.namespaces
        )
        dus = self.xml_doc.findall(
            './/gmd:onLine/gmd:CI_OnlineResource/gmd:linkage/gmd:URL',
            namespaces=self.namespaces
        )

        if dfs is not None:
            for i,j in enumerate(dfs):
                if j.text in dd.keys():
                    dd[j.text].append(dus[i].text)
                else:
                    dd[j.text] = [dus[i].text]
                    
        return(dd)

    
class Fetch:
    """Fetch class to fetch ftp/http data files"""
    
    def __init__(
            self,
            url=None,
            callback=fetches_callback,
            verbose=None,
            headers=r_headers,
            verify=True,
            allow_redirects=True
    ):
        self.url = url
        self.callback = callback
        self.verbose = verbose
        self.headers = headers
        self.verify = verify
        self.allow_redirects = allow_redirects

        
    def fetch_req(
            self,
            params=None,
            json=None,
            tries=5,
            timeout=None,
            read_timeout=None
    ):
        """fetch src_url and return the requests object"""
        
        if tries <= 0:
            utils.echo_error_msg(
                f'max-tries exhausted {self.url}'
            )
            raise ConnectionError(
                'Maximum attempts at connecting have failed.'
            )
        
        try:
            req = requests.get(
                self.url,
                stream=True,
                params=params,
                json=json,
                timeout=(timeout,read_timeout),
                #headers=self.headers,
                verify=self.verify,
                allow_redirects=self.allow_redirects
            )
        except Exception as e:
            ## there was an exception and we'll try again until tries is less than 1
            utils.echo_warning_msg(e)
            req = self.fetch_req(
                params=params,
                json=json,
                tries=tries - 1,
                #headers=self.headers,
                timeout=timeout * 2 if timeout is not None else None,
                read_timeout=read_timeout * 2 if read_timeout is not None else None
            )

        if req is not None:
            #utils.echo_msg(req.status_code)
            if req.status_code == 504:
                time.sleep(2)
                req = self.fetch_req(
                    params=params,
                    json=json,
                    tries=tries - 1,
                    #headers=self.headers,
                    timeout=timeout + 1 if timeout is not None else None,
                    read_timeout=read_timeout + 10 if read_timeout is not None else None
                )

            ## server says we have a band Range in the header, so we will
            ## remove the Range from the header and just try again.
            elif req.status_code == 416:
                if 'Range' in self.headers.keys():
                    del self.headers['Range']                    
                    req = self.fetch_req(
                        params=params,
                        json=json,
                        tries=tries - 1,
                        #headers=self.headers,
                        timeout=timeout + 1 if timeout is not None else None,
                        read_timeout=read_timeout + 10 if read_timeout is not None else None
                    )

            ## some unaccounted for return status code, report and exit.
            #req.status_code != 200 and req.status_code != 201:
            elif req.status_code < 200 or req.status_code > 299: 
                if self.verbose:
                    utils.echo_error_msg(
                        'request from {} returned {}'.format(
                            req.url, req.status_code
                        )
                    )
                #req = None
            
            return(req)

        
    def fetch_html(self, timeout=2):
        """fetch src_url and return it as an HTML object"""
    
        req = self.fetch_req(timeout=timeout)
        if req:
            return(lh.document_fromstring(req.text))
        else:
            return(None)

        
    def fetch_xml(self, timeout = 2, read_timeout = 10):
        """fetch src_url and return it as an XML object"""
    
        results = lxml.etree.fromstring(
            '<?xml version="1.0"?><!DOCTYPE _[<!ELEMENT _ EMPTY>]><_/>'.encode('utf-8')
        )
        try:
            req = self.fetch_req(
                timeout=timeout, read_timeout=read_timeout, headers=self.headers
            )
            results = lxml.etree.fromstring(req.text.encode('utf-8'))
        except:
            utils.echo_error_msg(
                f'could not access {self.url}'
            )
        return(results)

    
    def fetch_file(
            self,
            dst_fn,
            params=None,
            datatype=None,
            overwrite=False,
            timeout=30,
            read_timeout=None,
            tries=5,
            check_size=True
    ):
        """fetch src_url and save to dst_fn"""

        def retry(resume = False):
            if 'Range' in self.headers:
                del self.headers['Range']

            status = self.fetch_file(
                dst_fn,
                params=params,
                datatype=datatype,
                overwrite=overwrite,
                timeout=timeout+5 if timeout is not None else None,
                read_timeout=read_timeout+50 if read_timeout is not None else None,
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
            ## Set the Range header to the size of the existing file,
            ## unless check_size is False or overwrite is True.
            ## Otherwise, exit, FileExistsError will set status to 0,
            ## as if it completed a fetch.
            if not overwrite and os.path.exists(dst_fn):
                if not check_size:
                    raise FileExistsError(
                        f'{dst_fn} exists, '
                    )
                else:
                    dst_fn_size = os.stat(dst_fn).st_size
                    resume_byte_pos = dst_fn_size
                    #utils.echo_msg_bold('{} {}'.format(dst_fn_size, resume_byte_pos))
                    self.headers['Range'] = 'bytes={}-'.format(resume_byte_pos)

            with requests.get(
                    self.url, stream=True, params=params, headers=self.headers,
                    timeout=(timeout,read_timeout), verify=self.verify
            ) as req:

                ## requested range is not satisfiable, most likely the
                ## requestedrange is the complete size of the file,
                ## we'll skip here and assume that is the case.
                ## Set overwrite to True to overwrite instead
                req_h = req.headers
                if 'Content-Range' in req_h:
                    # this is wrong/hack
                    content_length = int(req_h['Content-Range'].split('/')[-1]) 
                    content_range = int(req_h['Content-Range'].split('/')[-1])
                    
                elif 'Content-Length' in req_h:
                    content_length = int(req_h['Content-Length'])                    
                else:
                    content_length = int(req_h.get('content-length', 0))

                req_s = content_length
                # if req_s != 0:
                #     utils.echo_msg(req_h)
                ## raise FileExistsError here if the file exists and the
                ## header Range value is the same as the requested
                ## content-length, unless overwrite is True or
                ## check_size is False.
                if not overwrite and check_size:
                    if os.path.exists(dst_fn):
                        if req_s == os.path.getsize(dst_fn) or req_s == 0:
                            raise FileExistsError(
                                f'{dst_fn} exists, '
                            )

                ## server returned bad content-length
                elif req_s == -1 or req_s == 0 or req_s == 49:
                    req_s = 0
                    
                if req.status_code == 416:
                    overwrite = True
                    raise FileExistsError(
                        '{} exists, and requested Range is invalid, {}'.format(
                            dst_fn, self.headers['Range']
                        )
                    )
                    
                ## redirect response. pass
                if req.status_code == 300:
                    pass

                ## hack for earthdata credential redirect...
                ## recursion here may never end with incorrect user/pass
                if req.status_code == 401:
                    ## we're hoping for a redirect url here.
                    if self.url == req.url:
                        raise UnboundLocalError('Incorrect Authentication')

                    ## re-run the Fetch with the new URL
                    status = Fetch(
                        url=req.url, headers=self.headers, verbose=self.verbose
                    ).fetch_file(
                        dst_fn,
                        params=params,
                        datatype=datatype,
                        overwrite=overwrite,
                        timeout=timeout,
                        read_timeout=read_timeout
                    )

                ## got a good response, so we'll attempt to fetch the file now.
                # or (req.status_code :
                elif (req.status_code == 200) or (req.status_code == 206):
                    curr_chunk = 0
                    total_size = int(req.headers.get('content-length', 0))                    
                    with open(
                            dst_fn, 'ab' if req.status_code == 206 else 'wb'
                    ) as local_file:
                        with tqdm(
                                desc='fetching: {}'.format(
                                    utils._init_msg(
                                        self.url, len('fetching: '), 40
                                    )
                                ),
                                total=req_s,
                                unit='iB',
                                unit_scale=True,
                                leave=self.verbose
                        ) as pbar:
                            try:
                                for chunk in req.iter_content(chunk_size = 8196):
                                    pbar.update(len(chunk))
                                    if not chunk:
                                        break

                                    local_file.write(chunk)
                                    local_file.flush()
                                    
                            except Exception as e:
                                ## reset the Fetch here if there was an exception
                                ## in fetching.
                                ## We'll attempt this `tries` times.
                                if tries != 0:
                                    if self.verbose:
                                        utils.echo_warning_msg(
                                            (f'server returned: {req.status_code}, '
                                             f'and an exception occured: {e}, '
                                             f'(attempts left: {tries})...')
                                        )
                                        
                                    time.sleep(2)
                                    status = Fetch(
                                        url=self.url,
                                        headers=self.headers,
                                        verbose=self.verbose
                                    ).fetch_file(
                                        dst_fn,
                                        params=params,
                                        datatype=datatype,
                                        overwrite=overwrite,
                                        timeout=timeout+5 \
                                        if timeout is not None \
                                        else None,
                                        read_timeout=read_timeout+50 \
                                        if read_timeout is not None \
                                        else None,
                                        tries=tries-1
                                    )
                                    self.verbose = False
                                else:
                                    raise e

                    ## something went wrong here and the size of the
                    ## fetched file does not match the requested
                    ## content-length
                    if check_size \
                       and (total_size != 0) \
                       and (total_size != os.stat(dst_fn).st_size):
                        raise UnboundLocalError(f'sizes do not match! {total_size} {os.stat(dst_fn).st_size}')

                ## 429: "Too Many Requests!"
                ## 416: "Bad header Range!"
                ## 504: "Gateway Timeout!"
                ## lets try again if we haven't already, these might resolve.
                elif (req.status_code == 429) \
                     or (req.status_code == 416) \
                     or (req.status_code == 504):
                    if tries != 0:
                        if self.verbose:
                            utils.echo_warning_msg(
                                'server returned: {}, (attempts left: {})...'.format(
                                    req.status_code, tries
                                )
                            )
                            
                        time.sleep(10)
                        status = Fetch(
                            url=self.url,
                            headers=self.headers,
                            verbose=self.verbose
                        ).fetch_file(
                            dst_fn,
                            params=params,
                            datatype=datatype,
                            overwrite=overwrite,
                            timeout=timeout+5 \
                            if timeout is not None \
                            else None,
                            read_timeout=read_timeout+50 \
                            if read_timeout is not None \
                            else None,
                            tries=tries-1
                        )
                else:
                    ## server returned some non-accounted-for status,
                    ## report and exit...
                    if self.verbose:
                        utils.echo_error_msg(
                            'server returned: {} ({})'.format(
                                req.status_code, req.url
                            )
                        )
                        
                    status = -1
                    raise ConnectionError(req.status_code)

        ## file exists, so we return status of 0,
        ## as if we were successful!
        except FileExistsError as e:
            status = 0

        ## other exceptions will return a status of -1,
        ## failure.
        except ConnectionError as e:
            status = -1
            raise e
        
        except requests.exceptions.ConnectionError as e:
            status = -1
            raise e#ConnectionError('Connection Aborted!')
        
        except UnboundLocalError as e:
            status = -1
            #utils.echo_msg(e)
            raise e#UnboundLocalError(e)
        
        except Exception as e:
            status = -1
            raise e#Exception(e)

        ## if the file exists now after all the above, make sure
        ## the size of that file is not zero, if `check_size`
        ## is True.
        if not os.path.exists(dst_fn) or os.stat(dst_fn).st_size ==  0:
            if check_size:
                status = -1
                raise UnboundLocalError('data not fetched')

        return(status)

    
    def fetch_ftp_file(self, dst_fn, params=None, datatype=None, overwrite=False):
        """fetch an ftp file via urllib"""

        status = 0
        f = None
        if self.verbose:
            utils.echo_msg(
                'fetching remote ftp file: {}...'.format(
                    self.url[:20]
                )
            )
            
        if not os.path.exists(os.path.dirname(dst_fn)):
            try:
                os.makedirs(os.path.dirname(dst_fn))
            except:
                pass
            
        try:
            f = urllib.request.urlopen(self.url)
        except Exception as e:
            utils.echo_error_msg(e)
            f = None
            status - 1

        if f is not None:
            with open(dst_fn, 'wb') as local_file:
                 local_file.write(f.read())
                 
            if self.verbose:
                utils.echo_msg(
                    'fetched remote ftp file: {}.'.format(
                        os.path.basename(self.url)
                    )
                )
                
        return(status)

    
## fetch queue for threads
def fetch_queue(q, c = True):
    """fetch queue `q` of fetch results

    each fetch queue `q` should be a list of the following:
    [remote_data_url, local_data_path, data-type, fetches-module, 
    number-of-attempts, results-list]    

    set c to False to skip size-checking
    """

    ## temporary bypass of ssl for certain modules...
    #no_verify = ['tnm', 'mar_grav', 'srtm_plus']
    no_verify = ['mar_grav', 'srtm_plus']
    while True:
        fetch_args = q.get()
        if not os.path.exists(os.path.dirname(fetch_args[1])):
            try:
                os.makedirs(os.path.dirname(fetch_args[1]))
            except: pass

        ## fetch either FTP or HTTP
        if fetch_args[0].split(':')[0] == 'ftp':
            Fetch(
                url=fetch_args[0],
                callback=fetch_args[3].callback,
                verbose=fetch_args[3].verbose,
                headers=fetch_args[3].headers
            ).fetch_ftp_file(fetch_args[1])
            #fetch_args[5].append([fetch_args[0], fetch_args[1], fetch_args[2]])
            fetch_results_ = [fetch_args[0], fetch_args[1], fetch_args[2], 0]
            fetch_args[5].append(fetch_results_)
            
            ## call the fetches callback function, does nothing
            ## unless reset by user, must be defined with a single
            ## argument, which is the fetch_results just populated
            if callable(fetch_args[3].callback):
                fetch_args[3].callback(fetch_results_)
        else:
            try:
                status = Fetch(
                    url=fetch_args[0],
                    callback=fetch_args[3].callback,
                    verbose=fetch_args[3].verbose,
                    headers=fetch_args[3].headers,
                    verify=False if fetch_args[3].name in no_verify else True
                ).fetch_file(fetch_args[1], check_size=c)
                fetch_results_ = [fetch_args[0], fetch_args[1], fetch_args[2], status]
                fetch_args[5].append(fetch_results_)

                ## call the fetches callback function, does nothing
                ## unless reset by user, must be defined with a single
                ## argument, which is the fetch_results just populated
                if callable(fetch_args[3].callback):
                    fetch_args[3].callback(fetch_results_)
                    
            except Exception as e:
                ## There was an exception in fetch_file, we'll put the request back into
                ## the queue to attempt to try again, fetch_args[4] is the number of times
                ## we will try to do this, once exhausted, we will give up.
                # and (utils.int_or(str(e), 0) < 400 or utils.int_or(str(e), 0) >= 500):
                if fetch_args[4] > 0:
                    # utils.echo_warning_msg(
                    #     'fetch of {} failed...putting back in the queue: {}'.format(
                    #         fetch_args[0], e
                    #     )
                    # )
                    fetch_args[4] -= 1
                    utils.remove_glob(fetch_args[1])
                    q.put(fetch_args)
                else:
                    utils.echo_error_msg(
                        'fetch of {} failed...'.format(fetch_args[0])
                    )
                    fetch_args[3].status = -1
                    fetch_results_ = [fetch_args[0], fetch_args[1], fetch_args[2], e]
                    fetch_args[5].append(fetch_results_)

                    ## call the fetches callback function, does nothing
                    ## unless reset by user, must be defined with a single
                    ## argument, which is the fetch_results just populated
                    if callable(fetch_args[3].callback):
                        fetch_args[3].callback(fetch_results_)

        q.task_done()

        
class fetch_results(threading.Thread):
    """fetch results gathered from a fetch module.

    results is a list of URLs with data type

    when a fetch module is run with {module}.run() it will 
    fill {module}.results with a list of urls, e.g.
    {module}.results = [[http://data/url.xyz.gz, /home/user/data/url.xyz.gz, data-type], ...]
    where each result in is [data_url, data_fn, data_type]

    run this on an initialized fetches module:
    >>> fetch_result(fetches_module, n_threads=3).run()
    and this will fill a queue for data fetching, using 
    'n_threads' threads.
    """
    
    def __init__(self, mod, check_size=True, n_threads=3, attempts=5):
        threading.Thread.__init__(self)
        self.fetch_q = queue.Queue()
        self.mod = mod
        self.check_size = check_size
        self.n_threads = n_threads
        self.attempts = attempts
        ## results holds the info from mod.results as a list,
        ## with the addition of the fetching status at the end.
        self.results = []
        if len(self.mod.results) == 0:
            self.mod.run()
            
    def run(self):
        for _ in range(self.n_threads):
            t = threading.Thread(
                target=fetch_queue,
                args=(self.fetch_q, self.check_size)
            )
            t.daemon = True
            t.start()

        ## fetch_q data is [fetch_results, fetch_path,
        ## fetch_dt, fetch_module, retries, results]
        attempts = self.attempts
        while True:
            for row in self.mod.results:
                self.fetch_q.put(
                    [row['url'],
                     os.path.join(
                         self.mod._outdir, row['dst_fn']
                     ),
                     row['data_type'],
                     self.mod,
                     self.attempts,
                     self.results]
                )

            self.fetch_q.join()
            status = [x[3]==0 for x in self.results]
            all_ok = len(status) == len(self.mod.results)                
            if (all(status) and all_ok) or attempts < 0:
                break
            else:
                attempts-=1

                
## Fetch Modules
class FetchModule:
    """The FetchModule super class to hold all the fetch modules.

    Make a sub-class from this to add a new fetch module, and 
    add it to the FetchesFactory to include it in the factory 
    for CLI or API.

    Each Fetch Module (sub-class) should define a `run` function 
    that will gather a list of `results`. 

    The `results` should be a list of 
    [remote-url, destination-file, data-type]
    """
    
    def __init__(
            self,
            src_region=None,
            callback=fetches_callback,
            verbose=True,
            outdir=None,
            name='fetches',
            params={}
    ):
        self.region = src_region # fetching region
        self.callback = callback # callback, run after a fetch attempt
        self.verbose = verbose # verbosity
        self.outdir = outdir # the directoy to place the fetched data
        self.params = params # FetchesFactory parameters
        self.status = 0 # fetching status
        self.results = [] # fetching results, a list of dicts
        self.name = name # the name of the fetch module
        
        ## some servers don't like us, or any 'bot' at all, so let's
        ## pretend we're just a Mozilla user on Linux. 
        #self.headers = { 'User-Agent': 'Fetches v%s' %(fetches.__version__) }
        self.headers = {
            'User-Agent': 'Mozilla/5.0 (X11; Linux x86_64; rv:45.0) Gecko/20100101 Firefox/45.0'
        }

        if self.outdir is None:
            self._outdir = os.path.join(os.getcwd(), self.name)
        else:
            self._outdir = os.path.join(self.outdir, self.name)

        ## for dlim support, we can check these variables for
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

        ## set a generic region of the entire world in WGS84 if no region
        ## was specified or if its an invalid region...this will result in quite
        ## a lot of downloading on global datasets, so be careful with this.
        if self.region is None or not self.region.valid_p():
            self.region = regions.Region().from_list([-180, 180, -90, 90])

            
    def run(self):
        """define the `run` function in the sub-class"""
        
        raise(NotImplementedError)

    
    def fetch_entry(self, entry, check_size=True, retries=5):
        try:
            status = Fetch(
                url=entry['url'],
                verbose=self.verbose,
                headers=self.headers,
            ).fetch_file(
                os.path.join(self._outdir, entry['dst_fn']),
                check_size=check_size
            )
        except:
            status = -1
        
        return(status)

    
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
        """add a fetches results entry to the results list"""
        
        entry = {'url': url, 'dst_fn': dst_fn, 'data_type': data_type}
        for key in kwargs.keys():
            entry[key] = kwargs[key]

        self.results.append(entry)             
        
### End
