### multibeam.py
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
##
## multibeam.py is part of CUDEM
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
## Fetch Multibeam bathymetry from NOAA NCEI, MBDB, and R2R.
##
### Code:

import os
import re
from io import StringIO
from typing import Optional, List, Dict, Any, Tuple
from cudem import utils
from cudem import regions
from cudem.fetches import fetches

## ==============================================
## Constants
## ==============================================
## NOAA NCEI
NCEI_DATA_URL = "https://data.ngdc.noaa.gov/platforms/"
NCEI_METADATA_URL = "https://data.noaa.gov/waf/NOAA/NESDIS/NGDC/MGG/Multibeam/iso/"
NCEI_SEARCH_URL = "https://gis.ngdc.noaa.gov/mapviewer-support/multibeam/files.groovy?"
NCEI_AUTOGRID_URL = "https://www.ngdc.noaa.gov/maps/autogrid/"
NCEI_HTML_URL = "https://www.ngdc.noaa.gov/"

## MBDB (ArcGIS)
MBDB_DYNAMIC_URL = 'https://gis.ngdc.noaa.gov/arcgis/rest/services/multibeam_footprints/MapServer'
MBDB_FEATURES_URL = 'https://gis.ngdc.noaa.gov/arcgis/rest/services/multibeam_datasets/FeatureServer'

## R2R
R2R_API_URL = 'https://service.rvdata.us/api/fileset/keyword/multibeam?'
R2R_PRODUCT_URL = 'https://service.rvdata.us/api/product/?'

## ==============================================
## Helper Functions
## ==============================================
def _parse_mbsystem_inf_bounds(inf_text: StringIO) -> 'regions.Region':
    """Parse spatial bounds from an MBSystem .inf file content."""
    
    minmax = [0, 0, 0, 0, 0, 0]
    for line in inf_text:
        parts = line.split()
        if len(parts) > 1 and parts[0] == 'Minimum':
            try:
                if parts[1] == 'Longitude:':
                    minmax[0] = utils.float_or(parts[2]) # xmin
                    minmax[1] = utils.float_or(parts[5]) # xmax
                elif parts[1] == 'Latitude:':
                    minmax[2] = utils.float_or(parts[2]) # ymin
                    minmax[3] = utils.float_or(parts[5]) # ymax
                elif parts[1] == 'Depth:':
                    minmax[4] = utils.float_or(parts[5]) * -1 # zmin
                    minmax[5] = utils.float_or(parts[2]) * -1 # zmax
            except (IndexError, ValueError):
                continue
    
    return regions.Region().from_list(minmax)


## ==============================================
## Multibeam Module (NCEI)
## ==============================================
class Multibeam(fetches.FetchModule):
    """NOAA MULTIBEAM bathymetric data.

    Fetch multibeam data from NOAA NCEI.
    
    Configuration Example:
    < multibeam:processed=True:min_year=None:max_year=None:survey_id=None:ship_id=None >
    """
    
    def __init__(
            self, 
            processed: bool = True, 
            survey_id: Optional[str] = None, 
            exclude_survey_id: Optional[str] = None, 
            ship_id: Optional[str] = None,
            exclude_ship_id: Optional[str] = None, 
            exclude: Optional[str] = None, 
            make_datalist: bool = False, 
            want_inf: bool = True,
            want_vdatum: bool = False, 
            inf_only: bool = False, 
            tables: bool = False, 
            **kwargs
    ):
        super().__init__(name='multibeam', **kwargs)
        self.processed_p = processed
        self.survey_id = survey_id
        self.exclude_survey_id = exclude_survey_id
        self.ship_id = ship_id
        self.exclude_ship_id = exclude_ship_id
        self.exclude = exclude
        self.make_datalist = make_datalist
        
        self.inf_only = inf_only
        self.want_inf = True if inf_only else want_inf
            
        self.want_vdatum = want_vdatum
        self.tables = tables
        
        ## Metadata
        self.data_format = 301 # MBSystem
        self.src_srs = 'epsg:4326+3855' if self.want_vdatum else 'epsg:4326'
        self.title = 'NOAA NCEI Multibeam bathymetric surveys'
        self.source = 'NOAA/NCEI'
        self.date = '1966 - 2022'
        self.data_type = 'Bathymetric Soundings'
        self.resolution = '~1m to ~30m'
        self.hdatum = 'WGS84'
        self.vdatum = 'LMSL'
        self.url = 'https://www.ngdc.noaa.gov/mgg/bathymetry/multibeam.html'

        
    def _extract_inf_metadata(self, src_inf: str, key: str) -> Optional[str]:
        """Generic extractor for local .inf files."""
        
        if not os.path.exists(src_inf):
            return None
            
        try:
            with open(src_inf, 'r', errors='ignore') as f:
                for line in f:
                    parts = line.split() if key != 'Number of Good Beams' else line.split(':')
                    if len(parts) > 1:
                        if key == 'MBIO' and parts[0] == 'MBIO':
                            return str(parts[4])
                        elif key == 'Time:' and parts[0] == 'Time:':
                            return parts[3]
                        elif key == 'Number of Good Beams' and parts[0].strip() == 'Number of Good Beams':
                            return parts[1].split()[-1].split('%')[0]
        except Exception:
            return None
        return None

    
    def mb_inf_data_format(self, src_inf: str) -> Optional[str]:
        return self._extract_inf_metadata(src_inf, 'MBIO')

    
    def mb_inf_data_date(self, src_inf: str) -> Optional[str]:
        return self._extract_inf_metadata(src_inf, 'Time:')

    
    def mb_inf_perc_good(self, src_inf: str) -> Optional[str]:
        return self._extract_inf_metadata(src_inf, 'Number of Good Beams')

    
    def check_for_generated_data(self, base_url: str) -> bool:
        """Check if a 'generated' directory exists for processed data."""
        
        ## Quick check on the specific file URL to see if it exists directly or needs 'generated' path        
        req = fetches.Fetch(base_url, verbose=False).fetch_req()
        if req is None or req.status_code == 404:
            parts = base_url.split('/')
            parts.insert(-1, 'generated')
            gen_url = '/'.join(parts)
            
            req_gen = fetches.Fetch(gen_url, verbose=False).fetch_req()
            if req_gen is not None and req_gen.status_code != 404:
                return True
        return False

    
    def run(self):
        """Run the multibeam fetches module."""
        
        if self.region is None:
            return []
        
        ## Search NCEI groovy script
        params = {'geometry': self.region.format('bbox')}
        req = fetches.Fetch(NCEI_SEARCH_URL).fetch_req(params=params, timeout=20)
        
        if req is None or req.status_code != 200:
            utils.echo_error_msg(f'Failed to fetch multibeam request: {req.status_code if req else "None"}')
            return []

        if self.verbose:
            utils.echo_msg(f"Query URL: {req.url}")
        
        ## Parse Results
        ## Structure: Survey -> Version -> List of files
        surveys_found = {} # {survey_name: {'date': date, 'versions': {ver: [files]}}}

        lines = req.text.split('\n')
        for line in lines:
            if not line.strip(): continue
            
            ## Line format is path-like, e.g.: data/../survey/ship/..
            parts = line.split(' ')[0].split('/') # Split path from potential trailing data
            if len(parts) < 10: continue

            ## Extract Metadata from path
            ## Typical Path: .../platforms/ocean/mgg/multibeam/data/version/SHIP/SURVEY/...
            survey = parts[6]
            ship = parts[5]
            version = parts[9][-1] # '1' or '2' usually
            filename = parts[-1]
            
            ## Construct Data URL
            data_url = f"{NCEI_DATA_URL}{'/'.join(line.split('/')[3:]).split(' ')[0]}"
            
            ## Date Extraction
            date_match = re.search(r"([0-9]{8})", filename)
            date_str = date_match.group(0) if date_match else None
            year = int(date_str[:4]) if date_str else None

            ## Filters
            if self.survey_id and survey not in self.survey_id.split('/'): continue
            if self.exclude_survey_id and survey in self.exclude_survey_id.split('/'): continue
            if self.ship_id and ship.lower() not in [x.lower() for x in self.ship_id.split('/')]: continue
            if self.exclude_ship_id and ship.lower() in [x.lower() for x in self.exclude_ship_id.split('/')]: continue
            
            if self.min_year and year and year < self.min_year: continue
            if self.max_year and year and year > self.max_year: continue

            ## Store
            if survey not in surveys_found:
                surveys_found[survey] = {'date': date_str, 'versions': {}}
            
            if version not in surveys_found[survey]['versions']:
                surveys_found[survey]['versions'][version] = []

            local_path = os.path.join(self._outdir, survey, filename)
            surveys_found[survey]['versions'][version].append([data_url, local_path, 'mb'])

        ## Process Survey List
        with utils.ccp(total=len(surveys_found), desc='Scanning NCEI Multibeam datasets', leave=self.verbose) as pbar:
            for survey, data in surveys_found.items():
                pbar.update()
                
                versions = data['versions']
                ## Prefer version '2' (Processed) if processed flag is True
                target_version = '2' if self.processed_p and '2' in versions else '1'
                
                if target_version not in versions:
                    ## Fallback if preferred missing.
                    if not self.processed_p:
                        ## Add all versions
                        for v in versions:
                            self._add_version_files(versions[v])
                        continue
                    else:
                        continue # Processed requested but v2 missing

                ## Process specific version files
                file_list = versions[target_version]
                if not file_list: continue

                ## Check for 'generated' directory (often holds the actual processed grids/data)
                use_generated = self.check_for_generated_data(file_list[0][0])
                
                final_files = []
                for f_entry in file_list:
                    url, dst, fmt = f_entry
                    if use_generated:
                        u_parts = url.split('/')
                        u_parts.insert(-1, 'generated')
                        url = '/'.join(u_parts)
                    
                    final_files.append([url, dst, fmt])

                self._add_version_files(final_files)

        ## Tables / Datalist Output
        if self.tables:
            for entry in self.results:
                inf_data = self.parse_entry_inf(entry, keep_inf=True) # Don't delete for table view?
                if inf_data:
                    print(f"{inf_data[0]} {inf_data[4]}")
                                
        if self.make_datalist:
            self._generate_datalist()
            
        return self

    
    def _add_version_files(self, file_list: List[List[str]]):
        """Helper to add files to results."""
        
        for entry in file_list:
            url, dst, fmt = entry
            
            ## Determine INF url
            ## .fbt files usually have .inf files with specific naming conventions
            inf_url = utils.fn_basename2(url) + '.inf' if url.endswith('fbt') else url + '.inf'
            
            if not self.inf_only:
                self.add_entry_to_results(url, dst, fmt)
            
            if self.want_inf:
                ## Add metadata file
                self.add_entry_to_results(inf_url, dst + '.inf', 'mb_inf')

                
    def parse_entry_inf(self, entry: Dict, keep_inf: bool = False) -> Optional[Tuple]:
        """Download and parse a local .inf file for an entry."""
        
        dst_fn = entry['dst_fn']
        ## If entry is the MB data, we need to find the INF. If it IS the inf, we read it.
        
        inf_dst = dst_fn + '.inf'
        
        if entry['url'].endswith('.inf'):
            inf_url = entry['url']
        else:
            base_url = utils.fn_basename2(entry['url']) if entry['url'].endswith('fbt') else entry['url']
            inf_url = base_url + '.inf'

        ## Fetch
        if fetches.Fetch(inf_url, callback=self.callback, verbose=False).fetch_file(inf_dst) == 0:
            survey = entry['url'].split('/')[7] 
            mb_fmt = self.mb_inf_data_format(inf_dst)
            mb_date = self.mb_inf_data_date(inf_dst)
            mb_perc = self.mb_inf_perc_good(inf_dst)
            
            if not keep_inf:
                utils.remove_glob(inf_dst)
            
            return (survey, os.path.basename(dst_fn), mb_fmt, mb_perc, mb_date)
        return None

    
    def _generate_datalist(self):
        """Generate a weighted datalist based on metadata."""
        
        seen = set()
        with open('mb_inf.txt', 'w') as f:
            for entry in self.results:
                if entry['format'] == 'mb_inf': continue # Skip INF entries themselves
                
                res = self.parse_entry_inf(entry)
                if not res: continue
                
                survey, _, _, mb_perc, mb_date = res
                if survey in seen: continue
                
                try:
                    yr = int(mb_date) if mb_date else 0
                    perc = float(mb_perc) if mb_perc else 0
                    
                    ## Weighting
                    ## (Year - 2000) factor
                    current_year = self.min_year if self.min_year else int(utils.this_year())
                    weight = perc * ((yr - 2000) / (current_year - 2000)) / 100.0
                    
                    f.write(f'{survey} -1 {weight}\n')
                    seen.add(survey)
                except (ValueError, TypeError):
                    pass

                
## ==============================================
## MBDB Module (ArcGIS)
## ==============================================
class MBDB(fetches.FetchModule):
    """MBDB fetching (ArcGIS REST Services).

    Configuration Example:
    < mbdb >
    """
    
    def __init__(self, where: str = '1=1', layer: int = 1, list_surveys: bool = False, want_inf: bool = True, **kwargs):
        super().__init__(name='mbdb', **kwargs)
        self.where = where        
        self.list_surveys = list_surveys
        self.want_inf = want_inf
        self._mb_features_query_url = f'{MBDB_FEATURES_URL}/{layer}/query?'

        
    def check_inf_region(self, mb_url: str) -> Tuple[str, Optional['regions.Region']]:
        """Fetch remote .inf file and parse its region."""
        
        ## Try finding the inf file
        src_mb = mb_url
        inf_url = f"{utils.fn_basename2(src_mb)}.inf"
        
        req = fetches.Fetch(inf_url, callback=self.callback, verbose=False).fetch_req()
        if req is None:
            ## Fallback
            inf_url = f"{src_mb}.inf"
            req = fetches.Fetch(inf_url, callback=self.callback, verbose=False).fetch_req()
        
        inf_region = None
        if req is not None and req.status_code == 200:
            with StringIO(req.text) as f:
                inf_region = _parse_mbsystem_inf_bounds(f)
                
        return inf_url, inf_region

    
    def run(self):
        """Run the MBDB fetching module."""
        
        if self.region is None:
            return []

        params = {
            'where': self.where,
            'outFields': '*',
            'geometry': self.region.format('bbox'),
            'inSR': 4326,
            'outSR': 4326,
            'f': 'pjson',
            'returnGeometry': 'false',
        }
        
        req = fetches.Fetch(self._mb_features_query_url, verbose=self.verbose).fetch_req(params=params)
        if req is None:
            return []
            
        features = req.json().get('features', [])
        
        with utils.ccp(total=len(features), desc='Parsing MBDB surveys', leave=self.verbose) as pbar:
            for feature in features:
                pbar.update()
                
                attrs = feature.get('attributes', {})
                download_url = attrs.get('DOWNLOAD_URL')
                
                if not download_url: continue
                
                if self.list_surveys:
                    print(download_url)
                    continue

                ## Scrape the download page for MB files
                page = fetches.Fetch(download_url, verbose=False).fetch_html()
                if page is None: continue
                
                ## Look for /MB/ links
                mb_links = page.xpath('//a[contains(@href, "/MB/")]/@href')
                
                for mb in mb_links:
                    ## Resolve full URL if relative. Usually absolute in these indexes, 
                    if 'gz' in mb:
                        mb_check_url = f"{utils.fn_basename2(mb)}.fbt"
                    else:
                        mb_check_url = mb

                    ## Check spatial bounds via INF
                    inf_url, inf_region = self.check_inf_region(mb_check_url)
                    
                    if inf_region and regions.regions_intersect_ogr_p(inf_region, self.region):
                        self.add_entry_to_results(
                            mb, 
                            os.path.join(self._outdir, os.path.basename(mb)),
                            'mbs'
                        )
                        
                        if self.want_inf:
                            self.add_entry_to_results(
                                inf_url,
                                os.path.join(self._outdir, os.path.basename(inf_url)),
                                'mb_inf'
                            )
        return self

    
## ==============================================
## R2R Module
## ==============================================
class R2R(fetches.FetchModule):
    """R2R Multibeam Fetching.

    Configuration Example:
    < r2r >
    """
    
    def __init__(self, check_inf: bool = False, **kwargs):
        super().__init__(name='R2R', **kwargs)
        self.check_inf = check_inf

        
    def run(self):
        """Run the R2R fetching module."""
        
        if self.region is None:
            return []

        params = {'spatial_bounds': self.region.export_as_wkt()}
        
        req = fetches.Fetch(R2R_API_URL, verbose=self.verbose).fetch_req(params=params)
        if req is None:
            return []
            
        data = req.json().get('data', [])
        
        with utils.ccp(total=len(data), desc='Parsing R2R datasets', leave=self.verbose) as pbar:
            for item in data:
                pbar.update()
                
                cruise_id = item.get('cruise_id')
                if not cruise_id: continue
                
                ## Fetch Products for Cruise
                prod_url = f"{R2R_PRODUCT_URL}cruise_id={cruise_id}"
                prod_req = fetches.Fetch(prod_url, verbose=False).fetch_req()
                
                if prod_req and prod_req.status_code == 200:
                    products = prod_req.json().get('data', [])
                    for prod in products:
                        if prod.get('datatype_name') == 'Bathymetry':
                            actual_url = prod.get('actual_url')
                            if actual_url:
                                self.add_entry_to_results(
                                    actual_url,
                                    os.path.basename(actual_url),
                                    'r2rBathymetry'
                                )
        return self

### End
