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
import requests
import urllib
import lxml.etree
import lxml.html as lh
import threading
try:
    import Queue as queue
except: import queue as queue
from cudem import utils
from cudem import fetches

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
r_headers = { 'User-Agent': 'Fetches v%s' %(fetches.__version__) }
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
        if len(self.results) == 0:
            self.run()
        for entry in self.results:
            for xyz in self.yield_xyz(entry, **kwargs):
                yield(xyz)
                
    def dump_results_to_xyz(self, dst_port=sys.stdout, **kwargs):
        for xyz in self.yield_results_to_xyz(**kwargs):
            xyz.dump(include_w=True if self.weight is not None else False,
                     dst_port=dst_port, encode=False)

### End
