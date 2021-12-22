### utils.py
##
## Copyright (c) 2010 - 2021 CIRES Coastal DEM Team
##
## utils.py is part of CUDEM
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
## Fetching Utility Functions
##
## Generic fetching and processing functions, etc.
##
## The `fetch_results` class will fetch a list of fetch results [[url, file-name, data-type]...]
## in a queue `fetch_queue` using 3 threads; set `p` and `s` as a fetch module object to processes
## and dump XYZ data from the fetched results, respectively. Use `fetch_file` to fetch single files.
##
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
            if geom: return(utils.wkt2geom(regions.create_wkt_polygon(opoly)))
            else: return(opoly)
        else: return(None)
        
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
        c = Fetch(self.url).fetch_req(params=_data)
        cx = lxml.etree.fromstring(c.text.encode('utf-8'))
        self.service_provider = cx.find('.//ows:ServiceProvider', namespaces = self.namespaces)
        self.service_identification = cx.find('.//ows:ServiceIdentification', namespaces = self.namespaces)
        self.operations_metadata = cx.find('.//ows:OperationsMetadata', namespaces = self.namespaces)
        self.service_metadata = cx.find('.//wcs:ServiceMetadata', namespaces = self.namespaces)
        self.contents = cx.find('.//wcs:Contents', namespaces = self.namespaces)

    def _contents(self):
        c = []
        for coverage in self.contents.xpath('//wcs:CoverageSummary', namespaces = self.namespaces):
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
        d = Fetch(url).fetch_req(params=_data)
        d_r = lxml.etree.fromstring(d.text.encode('utf-8'))
        cd = d_r.find('.//wcs:CoverageDescription', namespaces = self.namespaces)
        return(xml2py(d_r.find('.//wcs:CoverageDescription', namespaces = self.namespaces)))

    def _get_coverage_region(self, cov_desc):
        uc = [float(x) for x in cov_desc["boundedBy"]["Envelope"]["upperCorner"][0].split()]
        lc = [float(x) for x in cov_desc["boundedBy"]["Envelope"]["lowerCorner"][0].split()]
        return(regions.Region().from_list([lc[1], uc[1], lc[0], uc[0]]))
    
    def _get_coverage_url(self, coverage, region = None):
        dl_coverage = self.fix_coverage_id(coverage)
        cov_desc = self._describe_coverage(coverage)
        fmt = cov_desc["ServiceParameters"]["nativeFormat"][0]        
        hl = [float(x) for x in cov_desc["domainSet"]["RectifiedGrid"]["limits"]["GridEnvelope"]['high'][0].split()]
        uc = [float(x) for x in cov_desc["boundedBy"]["Envelope"]["upperCorner"][0].split()]
        lc = [float(x) for x in cov_desc["boundedBy"]["Envelope"]["lowerCorner"][0].split()]
        ds_region = regions.Region().from_list([lc[1], uc[1], lc[0], uc[0]])
        resx = (uc[1] - lc[1]) / hl[0]
        resy = (uc[0] - lc[0]) / hl[1]
        data = {'request': 'GetCoverage', 'version': '1.0.0', 'service': 'WCS',
                'resx': resx, 'resy': resy, 'crs': 'EPSG:4326', 'format': fmt,
                'coverage': coverage, 'Identifier': coverage}
        if region is not None: data['bbox'] = region.format('bbox')
        enc_data = urlencode(data)
        #try:
        #    enc_data = urllib.urlencode(data)
        #except: enc_data = urllib.parse.urlencode(data)
        return('{}{}'.format(self.url, enc_data))
    
    def fetch_coverage(coverage, region = None):
        c_url = self._get_coverage_url(coverage, region)
        return(Fetch(c_url, verbose=True).fetch_file('{}_{}.tif'.format(coverage, region.format('fn')), params=data))

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
            utils.echo_error_msg('max-tries exhausted')
            return(None)
        try:
            return(requests.get(self.url, stream=True, params=params, timeout=(timeout,read_timeout), headers=self.headers, verify=self.verify))
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
            progress = utils.CliProgress('fetching remote file: {}...'.format(self.url))
            
        if not os.path.exists(os.path.dirname(dst_fn)):
            try:
                os.makedirs(os.path.dirname(dst_fn))
            except: pass
            
        if not os.path.exists(dst_fn) or overwrite:
            try:
                with requests.get(self.url, stream=True, params=params, headers=self.headers,
                                  timeout=(timeout,read_timeout), verify=self.verify) as req:

                    req_h = req.headers
                    if 'Content-length' in req_h:
                        req_s = int(req_h['Content-length'])
                    else: req_s = -1

                    if req.status_code == 300:
                        #if req_h['Location']
                        pass
                    
                    ## ==============================================
                    ## hack for earthdata credential redirect...
                    ## recursion here may never end with incorrect user/pass
                    ## ==============================================
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
                        
                    elif req.status_code == 200:
                        curr_chunk = 0
                        with open(dst_fn, 'wb') as local_file:
                            for chunk in req.iter_content(chunk_size = 8196):
                                if self.callback():
                                    break
                                if self.verbose:
                                    progress.update_perc((curr_chunk, req_s))
                                    curr_chunk += 8196
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
            progress.end(status, 'fetched remote file as: {}.'.format(dst_fn))
            
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
                    Fetch(
                        url=fetch_args[0],
                        callback=m.callback,
                        verbose=m.verbose,
                        headers=m.headers
                    ).fetch_ftp_file(fetch_args[1])
                else:
                    Fetch(
                        url=fetch_args[0],
                        callback=m.callback,
                        verbose=m.verbose,
                        headers=m.headers,
                        verify=False if fetch_args[2] == 'srtm' or fetch_args[2] == 'mar_grav' else True
                    ).fetch_file(fetch_args[1])
            else:
                if m.region is not None:
                    o_x_fn = fetch_args[1] + m.region.format('fn') + '.xyz'
                else: o_x_fn = fetch_args[1] + '.xyz'
                
                utils.echo_msg('processing local file: {}'.format(o_x_fn))
                if not os.path.exists(o_x_fn):
                    with open(o_x_fn, 'w') as out_xyz:
                        m.dump_xyz(fetch_args, dst_port=out_xyz)
                        
                    try:
                        if os.path.exists(o_x_fn):
                            if os.stat(o_x_fn).st_size == 0:
                                utils.remove_glob(o_x_fn)
                                
                    except: pass
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

    def __init__(self, src_region=None, callback=lambda: False, weight=None, verbose=True, dst_srs=None):
        self.region = src_region
        self.callback = callback
        self.weight = weight
        self.verbose = verbose
        self.status = 0
        self.results = []
        self.dst_srs = dst_srs
        self.name = None
        self.headers = { 'User-Agent': 'Fetches v%s' %(fetches.__version__) }

    def run(self):
        raise(NotImplementedError)

    def yield_xyz(self, **kwargs):
        raise(NotImplementedError)
    
    def fetch(self, entry):
        Fetch(
            entry[0],
            callback=self.callback,
            verbose=self.verbose,
            headers=self.headers
        ).fetch_file(entry[1])

    def fetch_results(self):
        for entry in self.results:
            self.fetch(entry)
            
    def dump_xyz(self, entry, dst_port=sys.stdout, **kwargs):
        for xyz in self.yield_xyz(entry, **kwargs):
            xyz.dump(
                include_w=True if self.weight is not None else False,
                dst_port=dst_port,
                encode=False
            )
            
    def yield_results_to_xyz(self, **kwargs):
        if len(self.results) == 0:
            self.run()
            
        for entry in self.results:
            for xyz in self.yield_xyz(entry, **kwargs):
                yield(xyz)
                
    def dump_results_to_xyz(self, dst_port=sys.stdout, **kwargs):
        for xyz in self.yield_results_to_xyz(**kwargs):
            xyz.dump(
                include_w=True if self.weight is not None else False,
                dst_port=dst_port,
                encode=False
            )

### End
