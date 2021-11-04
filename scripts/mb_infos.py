#!/usr/bin/env python

import os
import requests
__version__ = '0.0.1'
r_headers = { 'User-Agent': 'Fetches v%s' %(__version__) }

def fetch_req(src_url, params = None, tries = 5, timeout = 2, read_timeout = 10):
    '''fetch src_url and return the requests object'''
    
    if tries <= 0:
        sys.stderr.write('error: max-tries exhausted')
        return(None)
    try:
        return(requests.get(src_url, stream = True, params = params, timeout = (timeout,read_timeout), headers = r_headers))
    except: return(fetch_req(src_url, params = params, tries = tries - 1, timeout = timeout + 1, read_timeout = read_timeout + 10))


def fetch_file(src_url, dst_fn, params = None, callback = lambda: False, datatype = None, overwrite = False, verbose = False):
    '''fetch src_url and save to dst_fn'''
    
    status = 0
    req = None
    halt = callback

    if not os.path.exists(os.path.dirname(dst_fn)):
        try:
            os.makedirs(os.path.dirname(dst_fn))
        except: pass 

    if not os.path.exists(dst_fn) or overwrite:
        try:
            with requests.get(src_url, stream = True, params = params, headers = r_headers, timeout=(45,320)) as req:
                req_h = req.headers
                if req.status_code == 200:
                    curr_chunk = 0
                    with open(dst_fn, 'wb') as local_file:
                        for chunk in req.iter_content(chunk_size = 8196):
                            if chunk:
                                if halt(): 
                                    status = -1
                                    break
                                local_file.write(chunk)
                else: utils.echo_error_msg('server returned: {}'.format(req.status_code))
                
        except Exception as e:
            utils.echo_error_msg(e)
            status = -1
    else:
        if os.path.exists(dst_fn): return(status)
        status = -1
    if not os.path.exists(dst_fn) or os.stat(dst_fn).st_size ==  0: status = -1
    if verbose: utils.echo_msg('fetched remote file: {}.'.format(os.path.basename(dst_fn)))
    return(status)


def get_mb_info_from_survey(survey_id):
    _req = fetch_req('ngdc.noaa.gov/')
    _req = f_utils.Fetch(self._mb_search_url).fetch_req(params={'geometry': self.region.format('bbox')}, timeout=20)

