### earthdata.py
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
##
### Code:

import os
import time
import datetime
from tqdm import tqdm
from osgeo import ogr
from cudem import utils
from cudem import regions
from cudem.fetches import fetches

## EarthData - NASA (requires login credentials)
class EarthData(fetches.FetchModule):
    """ACCESS NASA EARTH SCIENCE DATA
    
    NASA promotes the full and open sharing of all its data to research 
    and applications communities, private industry, academia, and the 
    general public. In order to meet the needs of these different communities, 
    NASA’s Earth Observing System Data and Information System (EOSDIS) has 
    provided various  ways to discover, access, and use the data.

    If version is omitted, will fetch all versions
    Use wildcards in 'short_name' to return granules for all matching 
    short_name entries.

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

    NASA promotes the full and open sharing of all its data to research and 
    applications communities, private industry, academia, and the general public. 
    In order to meet the needs of these different communities, NASA’s Earth 
    Observing System Data and Information System (EOSDIS) has provided various 
    ways to discover, access, and use the data.

    nsidc_download.py updated for fetches integration 12/21
    Updated from nsidc.py to earthdata.py to support all earthdata datasets

    time_start: A Zulu-time date string. e.g. '2020-05-04T00:00:00Z'

    time_end:   A Zulu-time date string. e.g. '2020-06-20T00:00:00Z'
                Leaving either time_start or time_end as a blank string ('') will
                default to searching from the start and/or end of the entire
                dataset collection, respectively.

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

    < earthdata:short_name=ATL08:version=004:time_start='':time_end='':filename_filter='' >
    """

    def __init__(self, short_name='ATL03', provider='', time_start='', time_end='',
                 version='', filename_filter=None, subset=False, **kwargs):
        super().__init__(name='cmr', **kwargs)
        self.short_name = short_name
        self.provider = provider
        self.time_start = time_start
        self.time_end = time_end
        self.version = version
        self.filename_filter = filename_filter
        self.subset = subset

        ## The various EarthData URLs
        self._cmr_url = 'https://cmr.earthdata.nasa.gov/search/granules.json?'
        self._harmony_url = f'https://harmony.earthdata.nasa.gov/ogc-api-edr/1.1.0/collections/{short_name}/cube?'

        ## Set up the earthdata credentials, and add it to our headers
        credentials = fetches.get_credentials(None)
        self.headers = {
            'Authorization': 'Basic {0}'.format(credentials),
            'User-Agent': ('Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:89.0) '
                           'Gecko/20100101 Firefox/89.0')
        }

        
    def add_wildcards_to_str(self, in_str):
        if not in_str.startswith('*'):
            in_str = '*' + in_str
            
        if not in_str.endswith('*'):
            in_str = in_str + '*'
            
        return(in_str)

    
    def run(self):
        """Run the earthdata fetches module"""
        
        if self.region is None:
            return([])

        _data = {
            'provider': self.provider,
            'short_name': self.short_name,
            'bounding_box': self.region.format('bbox'),
            'temporal': f'{self.time_start}, {self.time_end}',
            'page_size': 2000,
        }

        # if self.version != '':
        #     _data['version'] = f'self.version'
            
        if '*' in self.short_name:
            _data['options[short_name][pattern]'] = 'true'
        
        if self.filename_filter is not None:
            _data['options[producer_granule_id][pattern]'] = 'true'
            filename_filters = self.filename_filter.split(',')
            for filename_filter in filename_filters:
                _data['producer_granule_id'] = self.add_wildcards_to_str(filename_filter)

        if self.subset:
            _harmony_data = {
                'bbox': self.region.format('bbox'),
            }
            ## add time if specified
            if self.time_start != '' or self.time_end != '':
                start_time = datetime.datetime.fromisoformat(self.time_start).isoformat() + 'Z' if self.time_start != '' else '..'
                end_time = datetime.datetime.fromisoformat(self.time_end).isoformat() + 'Z' if self.time_end != '' else '..'
                _harmony_data['datetime'] = f'{start_time}/{end_time}'
                #'time': f'{self.time_start}, {self.time_end}',

            utils.echo_msg('requesting data subsets, please wait...')
            #utils.echo_msg(self._harmony_url)
            #utils.echo_msg(_harmony_data)
            status_url = None
            _req = fetches.Fetch(
                self._harmony_url
            ).fetch_req(
                params=_harmony_data, timeout=None, read_timeout=None
            )
            if _req is not None and _req.status_code == 200:
                status_json = _req.json()
                #utils.echo_msg(status_json['message'])
                #utils.echo_msg(status_json['Error'])
                utils.echo_msg(status_json)
                for link in status_json['links']:
                    if link['title'] == 'Job Status' or link['title'] == 'The current page':
                        status_url = link['href']

                if status_url is None:
                    if 'request' in status_json.keys():
                        status_url = status_json['request']

                if status_url is not None:                
                    with tqdm(
                            total=100,
                            desc='processing IceSat2 data',
                            leave=self.verbose
                    ) as pbar:
                    
                        while True:
                            _req = fetches.Fetch(status_url).fetch_req(timeout=None, read_timeout=None)
                            #utils.echo_msg(_req.status_code)
                            if _req is not None and _req.status_code == 200:
                                status = _req.json()
                                #utils.echo_msg(status)
                                pbar.n = status['progress']
                                pbar.refresh()
                                if status['status'] == 'successful':
                                    for link in status['links']:
                                        if link['href'].endswith('.h5'):
                                            self.add_entry_to_results(
                                                link['href'],
                                                os.path.basename(link['href']),
                                                f'{self.short_name} subset'
                                            )

                                    break
                                
                                elif status['status'] == 'running':
                                    time.sleep(10)
                                else:
                                    time.sleep(10)
                                    
                            else:
                                break
            else:
                utils.echo_warning_msg(f'failed to make subset request: {_req.status_code}')
            
        else:
            _req = fetches.Fetch(self._cmr_url).fetch_req(params=_data)
            if _req is not None:
                features = _req.json()['feed']['entry']
                for feature in features:
                    if 'polygons' in feature.keys():
                        poly = feature['polygons'][0][0]
                        cc = [float(x) for x in poly.split()]
                        gg = [x for x in zip(cc[::2], cc[1::2])]
                        ogr_geom = ogr.CreateGeometryFromWkt(regions.create_wkt_polygon(gg))
                        ## uncomment below to output shapefiles of the feature polygons
                        #regions.write_shapefile(ogr_geom, '{}.shp'.format(feature['title']))
                    else:
                        ogr_geom = self.region.export_as_geom()

                    if self.region.export_as_geom().Intersects(ogr_geom):
                        links = feature['links']
                        for link in links:
                            if link['rel'].endswith('/data#') and 'inherited' not in link.keys():
                                if not any([link['href'].split('/')[-1] in res for res in self.results]):
                                    self.add_entry_to_results(
                                        link['href'], link['href'].split('/')[-1], self.short_name
                                    )

                                    
## IceSat2 from EarthData shortcut - NASA (requires login credentials)
##
## This module allows us to use icesat2 data in dlim/waffles
##
## todo: dmrpp
class IceSat2(EarthData):
    """Access IceSat2 data.

    By default access ATL03 data, specify 'short_name' to fetch specific ATL data.

    If you wish, you may store your Earthdata username/password in a .netrc
    file in your $HOME directory and the script will automatically attempt to
    read this file. The .netrc file should have the following format:
    machine urs.earthdata.nasa.gov login myusername password mypassword
    where 'myusername' and 'mypassword' are your Earthdata credentials.

    you might need to `chmod 0600 ~/.netrc`
    
    < icesat2:short_name=ATL03:time_start='':time_end='':filename_filter='' >
    """
    
    def __init__(self, short_name='ATL03', subset=False, version='007', **kwargs):
        if short_name is not None:
            short_name = short_name.upper()
            if not short_name.startswith('ATL'):
                utils.echo_warning_msg(
                    '{} is not a valid icesat2 short_name, using ATL03'.format(short_name)
                )
                short_name = 'ATL03'

        if subset:
            atl08_v06_id = 'C2613553260-NSIDC_CPRD'
            atl03_v06_id = 'C2596864127-NSIDC_CPRD'
            atl03_v07_id = 'C3326974349-NSIDC_CPRD'
            if version == '007':
                short_name = atl03_v07_id
            elif version == '006':
                short_name = atl03_v06_id
                
        super().__init__(short_name=short_name, subset=subset, **kwargs)

        ## for dlim
        self.data_format = 303
        self.src_srs = 'epsg:4326+3855'
        #self.subset = subset

        
## SWOT from EarthData
## This module allows us to use SWOT data in dlim/waffles
class SWOT(EarthData):
    """Access SWOT data.

    https://podaac.jpl.nasa.gov/SWOT?tab=datasets-information

    set `product` to one of:

    Land-based KaRIn (HR mode)
    L2_HR_PIXC	netCDF
    L2_HR_PIXC_2	netCDF
    L2_HR_PIXC_1	netCDF
    L2_HR_PIXCVec	netCDF
    L2_HR_Raster	netCDF
    L2_HR_Raster_2	netCDF
    L2_HR_Raster_1	netCDF
    L2_HR_Raster_100m	netCDF
    L2_HR_Raster_250m	netCDF
    L2_HR_RiverSP	shapefile
    L2_HR_RiverSP_node	shapefile
    L2_HR_RiverSP_reach	shapefile
    L2_HR_LakeSP	shapefile
    L2_HR_LakeSP_obs	shapefile
    L2_HR_LakeSP_prior	shapefile
    L2_HR_LakeSP_unassigned	shapefile
    L2_HR_RiverAvg	shapefile
    L2_HR_LakeAvg	shapefile
    
    Ocean-based KaRIn (LR mode)
    L2_LR_SSH (2 km grid)	netCDF
    L2_LR_SSH_BASIC	netCDF
    L2_LR_SSH_WINDWAVE	netCDF
    L2_LR_SSH_EXPERT	netCDF
    L2_LR_SSH_UNSMOOTHED (250m)	netCDF

    Continuation of nadir altimetry
    L2_NALT_OGDR	netCDF
    L2_NALT_OGDR_SSHA	netCDF
    L2_NALT_OGDR_GDR	netCDF
    L2_NALT_IGDR	netCDF
    L2_NALT_IGDR_SSHA	netCDF
    L2_NALT_IGDR_GDR	netCDF
    L2_NALT_IGDR_SGDR	netCDF
    L2_NALT_GDR	netCDF
    L2_NALT_GDR_SSHA	netCDF
    L2_NALT_GDR_GDR	netCDF
    L2_NALT_GDR_SGDR	netCDF

    < swot:time_start='':time_end='':filename_filter='':product='' >
    """
    
    def __init__(self, product='L2_HR_Raster_2', **kwargs):
        super().__init__(short_name='SWOT_{}*'.format(product), **kwargs)
        self.src_srs = 'epsg:4326+3855'

        
## SST from EarthData
## This module allows us to use SST data in dlim/waffles
class MUR_SST(EarthData):
    """Access SST data.

    < mur_sst:time_start='':time_end='':filename_filter='' >
    """
    
    def __init__(self, **kwargs):
        super().__init__(short_name='MUR-JPL-L4-GLOB-v4.1', **kwargs)   


### End
