### EarthData Download Script
##
## Copyright (c) 2021 - 2022 Regents of the University of Colorado
##
## Permission is hereby granted, free of charge, to any person obtaining
## a copy of this software and associated documentation files (the "Software"),
## to deal in the Software without restriction, including without limitation
## the rights to use, copy, modify, merge, publish, distribute, sublicense,
## and/or sell copies of the Software, and to permit persons to whom the
## Software is furnished to do so, subject to the following conditions:
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
## https://earthdata.nasa.gov/
##
## Commentary from nsidc_download.py
## Tested in Python 2.7 and Python 3.4, 3.6, 3.7
##
## To run the script at a Linux, macOS, or Cygwin command-line terminal:
##   $ python nsidc-data-download.py
##
## On Windows, open Start menu -> Run and type cmd. Then type:
##     python nsidc-data-download.py
##
## The script will first search Earthdata for all matching files.
## You will then be prompted for your Earthdata username/password
## and the script will download the matching files.
##
## If you wish, you may store your Earthdata username/password in a .netrc
## file in your $HOME directory and the script will automatically attempt to
## read this file. The .netrc file should have the following format:
##    machine urs.earthdata.nasa.gov login myusername password mypassword
## where 'myusername' and 'mypassword' are your Earthdata credentials.
##
## you might need to `chmod 0600 ~/.netrc`
##
## NASA promotes the full and open sharing of all its data to research and applications communities, 
## private industry, academia, and the general public. In order to meet the needs of these different 
## communities, NASA’s Earth Observing System Data and Information System (EOSDIS) has provided various 
## ways to discover, access, and use the data.
##
## nsidc_download.py updated for fetches integration 12/21
## Updated from nsidc.py to earthdata.py to support all earthdata datasets
##
## some notable datasets:
##     ATL03
##     ATL06
##     ATL07
##     ATL08
##     GLAH06
##     GLAH14
##     GEDI01_B
##     GEDI02_A
##     GEDI02_B
##     ASTGTM
##     ASTGTM_NC
##     ASTL1A
##     ASTL1T
##     SRTMGL1
##     SRTMGL1_NC
##     SRTMGL3
##     SRTMGL3S
##     SRTMGL30
##     NASADEM_HGT
##     NASADEM_NC
##     NASADEM_SHHP
##     NASADEM_NUMNC
##     NASADEM_SIM
##
### Code:

from __future__ import print_function

import base64
import getopt
import itertools
import json
import math
import netrc
import os.path
import ssl
import sys
import time
from getpass import getpass

try:
    from urllib.parse import urlparse
    from urllib.request import urlopen, Request, build_opener, HTTPCookieProcessor
    from urllib.error import HTTPError, URLError
except ImportError:
    from urlparse import urlparse
    from urllib2 import urlopen, Request, HTTPError, URLError, build_opener, HTTPCookieProcessor

import numpy as np

from cudem import utils
from cudem import regions
from cudem import datasets
from cudem import xyzfun

import cudem.fetches.utils as f_utils

class EarthData(f_utils.FetchModule):

    def __init__(
            self,
            short_name='ATL08',
            provider=None,
            version=None,
            time_start='',
            time_end='',
            filename_filter = '',
            **kwargs
    ):
        super().__init__(**kwargs)
        self.CMR_URL = 'https://cmr.earthdata.nasa.gov'
        self.URS_URL = 'https://urs.earthdata.nasa.gov'
        self.CMR_PAGE_SIZE = 2000
        self.CMR_FILE_URL = ('{0}/search/granules.json?'
                             'sort_key[]=start_date&sort_key[]=producer_granule_id'
                             '&scroll=true&page_size={1}'.format(CMR_URL, CMR_PAGE_SIZE))
        self.name = 'earthdata'
        self._outdir = os.path.join(os.getcwd(), 'earthdata')
        
        self.short_name = short_name
        self.provider = provider
        self.version = version
        self.time_start = time_start
        self.time_end = time_end
        self.bounding_box = ''
        self.polygon = ''
        self.filename_filter = ''
        self.url_list = []
        self.dest_dir = None
        self.quiet = True

        credentials = get_credentials(None)

        #'User-Agent': 'Mozilla/5.0 (X11; Linux ppc64le; rv:75.0) Gecko/20100101 Firefox/75.0',
        self.headers = {'Authorization': 'Basic {0}'.format(credentials)}
        
    def run(self):
        url_list = cmr_search(
            self.short_name,
            self.provider,
            self.version,
            self.time_start,
            self.time_end,
            bounding_box=self.region.format('bbox'),
            polygon=self.polygon,
            filename_filter=self.filename_filter,
            quiet=self.quiet
        )

        for url in url_list:
            self.results.append([url, url.split('/')[-1], 'earthdata'])

    def yield_xyz(self, entry):

        ## IceSat Data
        import h5py
        ln = 0
        
        if entry[1].split('.')[-1] == 'h5' and 'ATL' in self.short_name:
            if f_utils.Fetch(
                    entry[0],
                    callback=self.callback,
                    verbose=self.verbose,
                    headers=self.headers
            ).fetch_file(entry[1]) == 0:
                h5_file = entry[1]
                h5 = h5py.File(h5_file, 'r')

                for g in h5['/']:
                    if 'gt' in g:
                        this_xyz = xyzfun.XYZPoint(w=1, src_srs='epsg:4326')
                        try:
                            h_ph = h5['{}/land_segments/dem_h'.format(g)]
                            lon_ph = h5['{}/land_segments/longitude'.format(g)]
                            lat_ph = h5['{}/land_segments/latitude'.format(g)]
                        except:
                            continue
                        
                        ofn = h5_file[:-3]
                        dataset = np.vstack((lon_ph, lat_ph, h_ph)).transpose()
                        if self.region is not None  and self.region.valid_p():        
                            dataset = dataset[dataset[:,0] > self.region.xmin,:]
                            dataset = dataset[dataset[:,0] < self.region.xmax,:]
                            dataset = dataset[dataset[:,1] > self.region.ymin,:]
                            dataset = dataset[dataset[:,1] < self.region.ymax,:]
                            if self.region.zmin is not None:
                                dataset = dataset[dataset[:,2] > self.region.zmin,:]
                                
                            if self.region.zmax is not None:
                                dataset = dataset[dataset[:,2] < self.region.zmax,:]
                
                        for point in dataset:
                            this_xyz.x = point[0]
                            this_xyz.y = point[1]
                            this_xyz.z = point[2]
                            this_xyz.w = self.weight

                            if self.dst_srs is not None:
                                this_xyz.warp(dst_srs=self.dst_srs)
                            
                            ln += 1
                            yield(this_xyz)

                        dataset = None

            if self.verbose:
                utils.echo_msg(
                    'parsed {} data records from {}'.format(ln, entry[1])
                )
            
            utils.remove_glob('{}*'.format(entry[1]))

## ==============================================
## nsidc_download.py from Mike McFerrin
## run nsidc.py directly to use non-fetches CLI
## ==============================================
CMR_URL = 'https://cmr.earthdata.nasa.gov'
URS_URL = 'https://urs.earthdata.nasa.gov'
CMR_PAGE_SIZE = 2000
# leave provider out of url, so we can search other datasets:
#CMR_FILE_URL = ('{0}/search/granules.json?provider=NSIDC_ECS'
#                '&sort_key[]=start_date&sort_key[]=producer_granule_id'
#                '&scroll=true&page_size={1}'.format(CMR_URL, CMR_PAGE_SIZE))
CMR_FILE_URL = ('{0}/search/granules.json?'
                'sort_key[]=start_date&sort_key[]=producer_granule_id'
                '&scroll=true&page_size={1}'.format(CMR_URL, CMR_PAGE_SIZE))

# Comments from Mike MacFerrin:
# Right now, the get_username() and get_password() just prompts the user for credentials, every time, which is fine (for now)
# In the future, we could (a) prompt the user once, then save an encrypted version of those credentials on the local file system.
# Perhaps use the computer name and the username as "salt" entries in the encryption.
#     It isn't a wholly secure system (if a hacker got ahold of those encrypted files
#     locally on your machine, they could read the code and figure out how to decrypt
#      the files to get your credentials), some thought should be put into how securely
#       do you want to make the system locally.
# Perhaps create a NOAA-DEM-group set of credentials and (interally) distribute those to use among us.
# There are other approaches that could be used, with varying levels of security.
def get_username():
    username = ''

    # For Python 2/3 compatibility:
    try:
        do_input = raw_input  # noqa
    except NameError:
        do_input = input

    while not username:
        username = do_input('Earthdata username: ')
    return username

def get_password():
    password = ''
    while not password:
        password = getpass('password: ')
    return password

def get_credentials(url):
    """Get user credentials from .netrc or prompt for input."""
    credentials = None
    errprefix = ''
    try:
        info = netrc.netrc()
        username, account, password = info.authenticators(urlparse(URS_URL).hostname)
        errprefix = 'netrc error: '
    except Exception as e:
        if (not ('No such file' in str(e))):
            print('netrc error: {0}'.format(str(e)))
        username = None
        password = None

    while not credentials:
        if not username:
            username = get_username()
            password = get_password()
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

    return credentials


def build_version_query_params(version):    
    desired_pad_length = 3
    if len(version) > desired_pad_length:
        print('Version string too long: "{0}"'.format(version))
        quit()
        
    version = str(int(version))  # Strip off any leading zeros
    query_params = ''
    while len(version) <= desired_pad_length:
        padded_version = version.zfill(desired_pad_length)
        query_params += '&version={0}'.format(padded_version)
        desired_pad_length -= 1
        
    return query_params


def filter_add_wildcards(filter):
    if not filter.startswith('*'):
        filter = '*' + filter
    if not filter.endswith('*'):
        filter = filter + '*'
    return filter


def build_filename_filter(filename_filter):
    filters = filename_filter.split(',')
    result = '&options[producer_granule_id][pattern]=true'
    for filter in filters:
        result += '&producer_granule_id[]=' + filter_add_wildcards(filter)
    return result


def build_cmr_query_url(short_name, version, provider, time_start, time_end,
                        bounding_box=None, polygon=None,
                        filename_filter=None):
    if '*' in short_name:
        params = '&short_name={0}&options[short_name][pattern]=true'.format(short_name)
    else:
        params = '&short_name={0}'.format(short_name)
        
    if provider is not None:
        params += '&provider={0}'.format(provider)
        
    if version is not None:
        params += build_version_query_params(version)
        
    params += '&temporal[]={0},{1}'.format(time_start, time_end)
    if polygon:
        params += '&polygon={0}'.format(polygon)
    elif bounding_box:
        params += '&bounding_box={0}'.format(bounding_box)
        
    if filename_filter:
        params += build_filename_filter(filename_filter)
        
    return CMR_FILE_URL + params


def get_speed(time_elapsed, chunk_size):
    if time_elapsed <= 0:
        return ''
    speed = chunk_size / time_elapsed
    if speed <= 0:
        speed = 1
    size_name = ('', 'k', 'M', 'G', 'T', 'P', 'E', 'Z', 'Y')
    i = int(math.floor(math.log(speed, 1000)))
    p = math.pow(1000, i)
    return '{0:.1f}{1}B/s'.format(speed / p, size_name[i])


def output_progress(count, total, status='', bar_len=60):
    if total <= 0:
        return
    fraction = min(max(count / float(total), 0), 1)
    filled_len = int(round(bar_len * fraction))
    percents = int(round(100.0 * fraction))
    bar = '=' * filled_len + ' ' * (bar_len - filled_len)
    fmt = '  [{0}] {1:3d}%  {2}   '.format(bar, percents, status)
    print('\b' * (len(fmt) + 4), end='')  # clears the line
    sys.stdout.write(fmt)
    sys.stdout.flush()


def cmr_read_in_chunks(file_object, chunk_size=1024 * 1024):
    """Read a file in chunks using a generator. Default chunk size: 1Mb."""
    while True:
        data = file_object.read(chunk_size)
        if not data:
            break
        yield data

# Comments from Mike MacFerrin:
# Put a local directory, either in a configuration file or hard-coded here,
# where you would like to store the data.
# Right now, just put it in a "data" directory in the parent directory of this script.
# Not perfect, but you/we can decide the best way to handle this.
def get_default_dest_dir(product_short_name="ATL08"):
    if product_short_name.upper() == "ATL06":
        dest_dir = os.path.join("..", "data", "ATL06")
    elif product_short_name.upper() == "ATL08":
        dest_dir = os.path.join("..", "data", "ATL08")
    else:
        dest_dir = os.path.join("..", "data", product_short_name)

    # Create the data directory if it doesn't yet exist. This assumes you have local permissions.
    if not os.path.exists(dest_dir):
        os.makedirs(dest_dir)

    return dest_dir

def cmr_download(urls, force=False, quiet=False, short_name="ATL08", dest_dir=None):
    """Download files from list of urls."""
    if not urls:
        return
    urls.sort()

    url_count = len(urls)
    if not quiet:
        print('Downloading {0} files...'.format(url_count))
    credentials = None

    if dest_dir is None:
        output_dir = get_default_dest_dir(short_name)

    for index, url in enumerate(urls, start=1):
        if not credentials and urlparse(url).scheme == 'https':
            credentials = get_credentials(url)

        filename = url.split('/')[-1]
        if not quiet:
            print('{0}/{1}: {2}'.format(str(index).zfill(len(str(url_count))),
                                        url_count, filename))

        try:
            if output_dir != None:
                filename_out = os.path.join(output_dir, os.path.split(filename)[1])

            req = Request(url)
            if credentials:
                req.add_header('Authorization', 'Basic {0}'.format(credentials))
            opener = build_opener(HTTPCookieProcessor())
            response = opener.open(req)
            length = int(response.headers['content-length'])
            try:
                if not force and length == os.path.getsize(filename_out):
                    if not quiet:
                        print('  File exists, skipping')
                    continue
            except OSError:
                pass
            count = 0
            chunk_size = min(max(length, 1), 1024 * 1024)
            max_chunks = int(math.ceil(length / chunk_size))
            time_initial = time.time()

            with open(filename_out, 'wb') as out_file:
                for data in cmr_read_in_chunks(response, chunk_size=chunk_size):
                    out_file.write(data)
                    if not quiet:
                        count = count + 1
                        time_elapsed = time.time() - time_initial
                        download_speed = get_speed(time_elapsed, count * chunk_size)
                        output_progress(count, max_chunks, status=download_speed)
            if not quiet:
                print()
        except HTTPError as e:
            print('HTTP error {0}, {1}'.format(e.code, e.reason))
        except URLError as e:
            print('URL error: {0}'.format(e.reason))
        except IOError:
            raise


def cmr_filter_urls(search_results):
    """Select only the desired data files from CMR response."""
    if 'feed' not in search_results or 'entry' not in search_results['feed']:
        return []

    entries = [e['links']
               for e in search_results['feed']['entry']
               if 'links' in e]
    # Flatten "entries" to a simple list of links
    links = list(itertools.chain(*entries))

    urls = []
    unique_filenames = set()
    for link in links:
        if 'href' not in link:
            # Exclude links with nothing to download
            continue
        if 'inherited' in link and link['inherited'] is True:
            # Why are we excluding these links?
            continue
        if 'rel' in link and 'data#' not in link['rel']:
            # Exclude links which are not classified by CMR as "data" or "metadata"
            continue

        if 'title' in link and 'opendap' in link['title'].lower():
            # Exclude OPeNDAP links--they are responsible for many duplicates
            # This is a hack; when the metadata is updated to properly identify
            # non-datapool links, we should be able to do this in a non-hack way
            continue

        filename = link['href'].split('/')[-1]
        if filename in unique_filenames:
            # Exclude links with duplicate filenames (they would overwrite)
            continue
        unique_filenames.add(filename)

        urls.append(link['href'])

    return urls


def cmr_search(short_name, provider, version, time_start, time_end,
               bounding_box='', polygon='', filename_filter='', quiet=False):
    """Perform a scrolling CMR query for files matching input criteria."""
    cmr_query_url = build_cmr_query_url(short_name=short_name, version=version,
                                        provider=provider, time_start=time_start, time_end=time_end,
                                        bounding_box=bounding_box,
                                        polygon=polygon, filename_filter=filename_filter)
    if not quiet:
        print('Querying for data:\n\t{0}\n'.format(cmr_query_url))

    cmr_scroll_id = None
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE

    urls = []
    hits = 0
    while True:
        req = Request(cmr_query_url)
        if cmr_scroll_id:
            req.add_header('cmr-scroll-id', cmr_scroll_id)
        response = urlopen(req, context=ctx)
        if not cmr_scroll_id:
            # Python 2 and 3 have different case for the http headers
            headers = {k.lower(): v for k, v in dict(response.info()).items()}
            cmr_scroll_id = headers['cmr-scroll-id']
            hits = int(headers['cmr-hits'])
            if not quiet:
                if hits > 0:
                    print('Found {0} matches.'.format(hits))
                else:
                    print('Found no matches.')
        search_page = response.read()
        search_page = json.loads(search_page.decode('utf-8'))
        url_scroll_results = cmr_filter_urls(search_page)
        if not url_scroll_results:
            break
        if not quiet and hits > CMR_PAGE_SIZE:
            print('.', end='')
            sys.stdout.flush()
        urls += url_scroll_results

    if not quiet and hits > CMR_PAGE_SIZE:
        print()
    return urls


def download_granules(short_name = 'ATL08',
                      provider = 'NSIDC_ECS',
                      version = '004',
                      time_start = '',
                      time_end = '',
                      bounding_box = '',
                      polygon = '',
                      filename_filter = '',
                      url_list = [],
                      dest_dir = None,
                      argv=None):
    """Download NSIDC data granules.

    The function *should* work with nearly any NSIDC dataset, but it hasn't been
    thoroughly tested on them all. It comes with a variety of optional parameters to use.

    Parameters
    ----------
    short_name: the NSIDC short name of the datset to use. Documentation of various
                datasets can be found at https://nsidc.org/data/
                For ICESat-2 data, derived data products from the ATLAS instrument
                all have the format "ATLXX", with XX being the 2-digist product data
                number, e.g. "ATL03". ATLAS product documentation can be found here:
                https://icesat-2.gsfc.nasa.gov/science/data-products

    provider: dataset Provider. e.g. "NSIDC_ECS"

    version: Version number of the dataset, as a 3-digit integer string. e.g. "004"

    time_start: A Zulu-time date string. e.g. '2020-05-04T00:00:00Z'

    time_end:   A Zulu-time date string. e.g. '2020-06-20T00:00:00Z'
                Leaving either time_start or time_end as a blank string ('') will
                default to searching from the start and/or end of the entire
                dataset collection, respectively.

    bounding_box: A string of lat/lon bounding box coordinates, using WGS84 decimal degrees.
                Defined as 'lon_min, lat_min, lon_max, lat_max'.
                e.g. '-114.12,42.91,-101.1,48.74'
                A blank string searches the whole globe, unless the polygon parameter is used.

    polygon:    A polygon of coordinates, using WGS84 decimal degrees.
                Defined as 'p0_lon, p0_lat, p1_lon, p1_lat, ... , pN_lon, pN_lat, p0_lon, p0_lat'
                Coordinates should be in a loop, with the final coordinate equal to the first.
                Coordinates can be clockwise or counter-clockwise.
                Holes are not supported.
                e.g. '-109.89,41.63,-117.02,36.30,-107.41,35.21,-108.85,39.41,-109.89,41.63'
                A blank string searches the whole globe, unless the bounding_box parameter is used.

                bounding_box and polygon parameters should be used exclusively, not together.

    filename_filter: A string filter for the filenames.
                open-ended * flags can be used.
                e.g. '*20200101202129_00930611*'

    url_list:   A specific list of urls to acquire. Defaults to a blank list.
                This parameter should be used exclusively with any other search
                parameters listed above.

    dest_dir:   A destination directory in which to put the downloaded files.
                If a directory is given, the directory should already exist.
                Defaults to None, which will choose a default directory
                in the location "../data/[short_name]". If this directory doesn't
                exist, it will be created.

    argv:       A string of argv command-line options. If None, argv will be
                collected from the actual command-line.
                Available options:
                [--help, -h] [--force, -f] [--quiet, -q]'

    Returns
    -------
    None




    """

    if argv is None:
        argv = sys.argv[1:]

    force = False
    quiet = False
    usage = 'usage: nsidc-download.py [--help, -h] [--force, -f] [--quiet, -q]'

    try:
        opts, args = getopt.getopt(argv, 'hfq', ['help', 'force', 'quiet'])
        for opt, _arg in opts:
            if opt in ('-f', '--force'):
                force = True
            elif opt in ('-q', '--quiet'):
                quiet = True
            elif opt in ('-h', '--help'):
                print(usage)
                sys.exit(0)
    except getopt.GetoptError as e:
        print(e.args[0])
        print(usage)
        sys.exit(1)

    if dest_dir is None:
        dest_dir = get_default_dest_dir(short_name)

    try:
        if not url_list:
            url_list = cmr_search(short_name, provider, version, time_start, time_end,
                                  bounding_box=bounding_box, polygon=polygon,
                                  filename_filter=filename_filter, quiet=quiet)

        cmr_download(url_list, force=force, quiet=quiet, output_dir=dest_dir)
    except KeyboardInterrupt:
        quit()


if __name__ == '__main__':
    download_granules()

    # Get the ATL08, v4 granules for 2020. (Land/Veg elevation)
    # ATL08 is the global Land and Canopy elevation dataset.
    # download_granules(short_name='ATL08',
    #                   version='004',
    #                   time_start='2020-05-04T00:00:00Z',
    #                   time_end='2021-04-20T00:00:00Z')

    # Other parameters to use, examples:
    # filename_filter if you have an idea of specific files you want to use. Can put * strings in there (often you want both the .h5 and the .iso.xml files).
    #                   # filename_filter="*20200101202129_00930611*")
    # bounding_box is defined as 'lon_min, lat_min, lon_max, lat_max', as a string.
                        # bounding_box = '-114.12,42.91,-101.1,48.74',
    # polygon is defined as 'p0_lon, p0_lat, p1_lon, p1_lat, ... , pN_lon, pN_lat, p0_lon, p0_lat', as a string.
                        # polygon = '-109.8957,41.6345,-117.0209,36.3045,-107.4134,35.2118,-110.0575,37.7867,-99.5624,37.2008,-108.8538,39.4172,-109.8957,41.6345',


    # Also get the ATL06 granules for 2020. (Land Ice measurements)
    # ATL06 is the Global Land Ice elevation dataset.
    # download_granules(short_name = 'ATL06',
    #                   version='004',
    #                   time_start='2020-11-23T00:00:00Z',
    #                   time_end='2021-04-20T00:00:00Z')

### End
