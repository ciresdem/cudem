### earthdata.py
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
##
## earthdata.py is part of CUDEM
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
## Fetch data from NASA's EarthData CMR and Harmony API.
##
### Code:

import os
import time
import datetime
from typing import List, Dict, Optional, Any, Union
from osgeo import ogr
from cudem import utils
from cudem import regions
from cudem.fetches import fetches

## ==============================================
## Constants
## ==============================================
CMR_SEARCH_URL = 'https://cmr.earthdata.nasa.gov/search/granules.json?'
HARMONY_BASE_URL = 'https://harmony.earthdata.nasa.gov'

## ==============================================
## EarthData Module
## ==============================================
class EarthData(fetches.FetchModule):
    """ACCESS NASA EARTH SCIENCE DATA
    
    NASA promotes the full and open sharing of all its data.
    Requires ~/.netrc credentials.

    time_start: A Zulu-time date string. e.g. '2020-05-04T00:00:00Z'
    time_end:   A Zulu-time date string. e.g. '2020-06-20T00:00:00Z'

    < earthdata:short_name=ATL03:version='':time_start='':time_end='':filename_filter='' >
    """

    def __init__(self, 
                 short_name: str = 'ATL03', 
                 provider: str = '', 
                 time_start: str = '', 
                 time_end: str = '',
                 version: str = '', 
                 filename_filter: Optional[str] = None, 
                 subset: bool = False, 
                 subset_job_id: Optional[str] = None,
                 harmony_ping: Optional[str] = None, 
                 **kwargs):
        super().__init__(name='cmr', **kwargs)
        self.short_name = short_name
        self.provider = provider
        self.time_start = time_start
        self.time_end = time_end
        self.version = version
        self.filename_filter = filename_filter
        self.subset = subset
        self.subset_job_id = subset_job_id
        self.harmony_ping = harmony_ping  # 'status', 'pause', 'resume' or 'cancel'
        
        ## URLs
        self._cmr_url = CMR_SEARCH_URL
        self._harmony_url = f'{HARMONY_BASE_URL}/ogc-api-edr/1.1.0/collections/{short_name}/cube?'

        ## Authentication
        credentials = fetches.get_credentials(None)
        if credentials:
            self.headers = {
                'Authorization': f'Basic {credentials}',
                'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:89.0) Gecko/20100101 Firefox/89.0'
            }
        else:
            self.headers = {}
            utils.echo_warning_msg("Could not retrieve EarthData credentials.")

            
    def add_wildcards_to_str(self, in_str: str) -> str:
        """Ensure wildcards exist at start/end of string."""
        
        if not in_str.startswith('*'):
            in_str = f'*{in_str}'
        if not in_str.endswith('*'):
            in_str = f'{in_str}*'
        return in_str


    def _format_date(self, date_str: str) -> str:
        """Formats an ISO date string for filtering."""
        
        if not date_str:
            return '..'
        try:
            dt = datetime.datetime.fromisoformat(date_str)
            return dt.isoformat(timespec='milliseconds') + 'Z'
        except ValueError:
            return '..'

    
    def harmony_ping_for_status(self, job_id: str, ping_request: str = 'status') -> Optional[Dict]:
        """Check status of a Harmony Job."""
        
        valid_requests = ['status', 'pause', 'resume', 'cancel']
        base_url = f'{HARMONY_BASE_URL}/jobs/{job_id}'

        if ping_request in valid_requests[1:]:
            status_url = f'{base_url}/{ping_request}'
        else:
            status_url = base_url

        req = fetches.Fetch(status_url, headers=self.headers).fetch_req(timeout=None, read_timeout=None)
        if req and req.status_code == 200:
            return req.json()
        return None

    
    def harmony_make_request(self) -> Optional[Dict]:
        """Initiate a Harmony Subset Request."""
        
        if not self.region:
            return None

        harmony_data = {
            'bbox': self.region.format('bbox'),
        }

        ## Format Temporal Params for Harmony
        start_t = self._format_date(self.time_start)
        end_t = self._format_date(self.time_end)
        
        if start_t != '..' or end_t != '..':
            harmony_data['datetime'] = f'{start_t}/{end_t}'

        req = fetches.Fetch(
            self._harmony_url, headers=self.headers
        ).fetch_req(
            params=harmony_data, timeout=None, read_timeout=None
        )

        if req and req.status_code in [200, 201, 202]:
            return req.json()
        
        utils.echo_error_msg(f"Harmony request failed: {req.status_code if req else 'No Response'}")
        if req: utils.echo_msg(req.text)
        return None

    
    def earthdata_set_config(self) -> Dict:
        """Configure CMR Search Parameters."""
        
        data = {
            'provider': self.provider,
            'short_name': self.short_name,
            'bounding_box': self.region.format('bbox'),
            'temporal': f'{self.time_start},{self.time_end}',
            'page_size': 2000,
        }

        if self.version:
            data['version'] = self.version
            
        if '*' in self.short_name:
            data['options[short_name][pattern]'] = 'true'
        
        if self.filename_filter:
            data['options[producer_granule_id][pattern]'] = 'true'
            filters = self.filename_filter.split(',')
            for f in filters:
                data['producer_granule_id'] = self.add_wildcards_to_str(f)

        return data

    
    def _run_cmr_search(self):
        """Execute standard CMR Granule Search."""
        
        params = self.earthdata_set_config()
        req = fetches.Fetch(self._cmr_url).fetch_req(params=params)
        
        if not req:
            return

        try:
            feed = req.json().get('feed', {})
            entries = feed.get('entry', [])
        except Exception as e:
            utils.echo_error_msg(f"Error parsing CMR response: {e}")
            return

        for entry in entries:
            ## Spatial Filtering
            ## Refine BBox search with exact polygon intersection if available.
            geom_valid = True
            if 'polygons' in entry:
                try:
                    poly_str = entry['polygons'][0][0]
                    coords = [float(x) for x in poly_str.split()]
                    points = list(zip(coords[::2], coords[1::2])) # Lon, Lat
                    ogr_geom = ogr.CreateGeometryFromWkt(regions.create_wkt_polygon(points))
                    ## Ucomment line below to output shapefiles of the polygons
                    #regions.write_shapefile(ogr_geom, '{}.shp'.format(feature['title']))
                    if not self.region.export_as_geom().Intersects(ogr_geom):
                        geom_valid = False
                except Exception:
                    pass
            
            if geom_valid:
                for link in entry.get('links', []):
                    ## Filter for data links
                    #utils.echo_msg(link)
                    if link.get('rel', '').endswith('/data#') and 'inherited' not in link:
                        href = link.get('href')
                        if href:
                            fname = href.split('/')[-1]
                            self.add_entry_to_results(href, fname, self.short_name)

                            
    def _run_harmony_subset(self):
        """Execute Harmony Subset Job."""
        
        if not self.subset_job_id:
            status = self.harmony_make_request()
            if status and 'jobID' in status:
                self.subset_job_id = status['jobID']
                utils.echo_msg(f"Harmony Job Initiated: {self.subset_job_id}")
            else:
                return

        if self.subset_job_id:
            with utils.ccp(total=100, desc=f'Harmony Job ({self.subset_job_id})', leave=self.verbose) as pbar: 
                while True:
                    try:
                        status = self.harmony_ping_for_status(self.subset_job_id)
                        if not status:
                            time.sleep(10)
                            continue

                        progress = status.get('progress', 0)
                        state = status.get('status', 'unknown')
                        
                        pbar.n = int(progress)
                        pbar.set_description(f'Harmony Job {self.subset_job_id} -- ({state})')
                        pbar.refresh()

                        if state == 'successful':
                            for link in status.get('links', []):
                                href = link.get('href', '')
                                if href.endswith('.h5') or href.endswith('.nc'):
                                    base_name = utils.fn_basename2(os.path.basename(href))
                                    out_fn = f'{utils.append_fn(base_name, self.region, 1)}.h5'
                                    
                                    self.add_entry_to_results(
                                        href,
                                        out_fn,
                                        f'{self.short_name} subset'
                                    )
                            break
                        
                        elif state in ['failed', 'canceled']:
                            utils.echo_error_msg(f"Harmony Job {state}: {status.get('message', '')}")
                            break
                        
                        else:
                            time.sleep(15)
                                
                    except Exception as e:
                        utils.echo_error_msg(f'Harmony polling failed: {e}')
                        time.sleep(15)

                        
    def run(self):
        """Run the EarthData fetch module."""
        
        ## Handle manual ping if requested
        if self.harmony_ping:
            if self.subset_job_id:
                status = self.harmony_ping_for_status(self.subset_job_id, self.harmony_ping)
                if status:
                    utils.echo_msg(status)
                else:
                    utils.echo_warning_msg(f'Bad Harmony ping: {self.harmony_ping}')
            return []
        
        if self.region is None:
            return []

        if not self.subset:
            self._run_cmr_search()
        else:
            self._run_harmony_subset()
            
        return self

    
## ==============================================
## Subclasses / Shortcuts
## ==============================================
class IceSat2(EarthData):
    """Access IceSat2 data (Shortcuts for ATL03/ATL08)."""
    
    def __init__(self, short_name: str = 'ATL03', subset: bool = False, version: str = '007', **kwargs):
        
        ## Normalize Short Name
        if short_name:
            short_name = short_name.upper()
            if not short_name.startswith('ATL'):
                utils.echo_warning_msg(f'{short_name} is invalid, defaulting to ATL03')
                short_name = 'ATL03'

        if short_name == 'ATL24' and version == '':
            version = '001'
                
        ## Handle Subset Collection IDs for Harmony
        if subset:
            ## Collection IDs map short names to specific provider/version IDs required by Harmony
            collection_map = {
                '007': {'ATL03': 'C3326974349-NSIDC_CPRD', 'ATL08': 'C3565574177-NSIDC_CPRD'},
                '006': {'ATL03': 'C2596864127-NSIDC_CPRD', 'ATL08': 'C2613553260-NSIDC_CPRD'}
            }
            
            if version in collection_map and short_name in collection_map[version]:
                short_name = collection_map[version][short_name]
                
        super().__init__(short_name=short_name, subset=subset, version=version, **kwargs)
        
        ## DLIM Defaults
        self.data_format = 303
        self.src_srs = 'epsg:4326+3855'

        
class SWOT(EarthData):
    """Access SWOT data."""
    
    def __init__(self, product: str = 'L2_HR_Raster_2', **kwargs):
        super().__init__(short_name=f'SWOT_{product}*', **kwargs)
        self.src_srs = 'epsg:4326+3855'

        
class MUR_SST(EarthData):
    """Access SST data."""
    
    def __init__(self, **kwargs):
        super().__init__(short_name='MUR-JPL-L4-GLOB-v4.1', **kwargs)   

        
### End
