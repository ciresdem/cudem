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
##
### Code:

import os
import re
from tqdm import tqdm
from io import StringIO
from cudem import utils
from cudem import regions
from cudem.fetches import fetches

class Multibeam(fetches.FetchModule):
    """NOAA MULTIBEAM bathymetric data.

    Fetch multibeam data from NOAA NCEI
        
    NCEI is the U.S. national archive for multibeam bathymetric data 
    and holds more than 9 million nautical miles of ship trackline 
    data recorded from over 2400 cruises and received from sources 
    worldwide.

    https://data.ngdc.noaa.gov/platforms/

    <exclude_>survey_id and <exclude_>ship_id can be lists of surveys or ships, 
    repsectively, using a '/' as a seperator.

    < multibeam:processed=True:min_year=None:max_year=None:survey_id=None:ship_id=None:exclude_survey_id=None:exclude_ship_id=None >
    """
    
    def __init__(
            self, processed=True, survey_id=None, exclude_survey_id=None, ship_id=None,
            exclude_ship_id=None, min_year=None, max_year=None, exclude=None,
            make_datalist=False, want_inf=True, want_vdatum=False, inf_only=False, **kwargs
    ):
        super().__init__(name='multibeam', **kwargs)
        self.processed_p = processed
        self.min_year = utils.int_or(min_year)
        self.max_year = utils.int_or(max_year)
        self.survey_id = survey_id
        self.exclude_survey_id = exclude_survey_id
        self.ship_id = ship_id
        self.exclude_ship_id = exclude_ship_id
        self.exclude = exclude
        self.make_datalist = make_datalist
        self.want_inf = want_inf
        self.inf_only = inf_only
        if self.inf_only:
            self.want_inf = True
            
        self.want_vdatum = want_vdatum

        ## various multibeam URLs
        self._mb_data_url = "https://data.ngdc.noaa.gov/platforms/"
        self._mb_metadata_url = "https://data.noaa.gov/waf/NOAA/NESDIS/NGDC/MGG/Multibeam/iso/"
        self._mb_search_url = "https://gis.ngdc.noaa.gov/mapviewer-support/multibeam/files.groovy?"
        self._mb_autogrid = "https://www.ngdc.noaa.gov/maps/autogrid/"
        self._mb_html = "https://www.ngdc.noaa.gov/"
        self._mb_json = "https://data.ngdc.noaa.gov/platforms/ocean/mgg/mb-archive.json"
        self._urls = [self._mb_data_url, self._mb_metadata_url, self._mb_autogrid]        

        ## for dlim, data_format of 301 is multibeam data parsed with MBSystem.
        self.data_format = 301
        if self.want_vdatum:
            self.src_srs = 'epsg:4326+3855'
        else:
            self.src_srs = 'epsg:4326'

        self.title = 'NOAA NCEI Multibeam bathymetric surveys'
        self.source = 'NOAA/NCEI'
        self.date = '1966 - 2022'
        self.data_type = 'Bathymetric Soundings'
        self.resolution = '~1m to ~30m'
        self.hdatum = 'WGS84'
        self.vdatum = 'LMSL'
        self.url = 'https://www.ngdc.noaa.gov/mgg/bathymetry/multibeam.html'

        
    def mb_inf_data_format(self, src_inf):
        """extract the data format from the mbsystem inf file."""

        with open(src_inf, errors='ignore') as iob:
            for il in iob:
                til = il.split()
                if len(til) > 1:
                    if til[0] == 'MBIO':
                        return('{}'.format(til[4]))

                    
    def mb_inf_data_date(self, src_inf):
        """extract the date from the mbsystem inf file."""

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


    def want_generated(self, tmp_url):
        #tmp_url = these_surveys[key]['1'][0][0]
        _req = fetches.Fetch(
            tmp_url,
            verbose=False
        ).fetch_req()
        if _req is None or _req.status_code == 404:
            tmp_url = tmp_url.split('/')
            tmp_url.insert(-1, 'generated')
            tmp_url = '/'.join(tmp_url)
            
            _req = fetches.Fetch(
                tmp_url,
                verbose=False
            ).fetch_req()
            if _req is not None and _req.status_code != 404:
                want_generated = True
        else:
            want_generated = False    

        return(want_generated)
                    
    def run(self):
        """Run the multibeam fetches module"""
        
        these_surveys = {}
        these_versions = {}
        if self.region is None:
            return([])
        else:
            fetch_region = self.region.copy()

        _req = fetches.Fetch(
            self._mb_search_url
        ).fetch_req(
            params={'geometry': fetch_region.format('bbox')},
            timeout=20
        )
        if _req is not None and _req.status_code == 200:
            utils.echo_msg(_req.url)
            survey_list = _req.text.split('\n')[:-1]
            #utils.echo_msg(survey_list)
            for r in survey_list:
                dst_pfn = r.split(' ')[0]
                dst_p = dst_pfn.split('/')
                dst_fn = dst_p[-1:][0]
                survey = dst_p[6]
                dn = dst_p[:-1]
                version = dst_p[9][-1]
                ship = dst_p[5]
                data_url = self._mb_data_url + '/'.join(r.split('/')[3:])
                date = re.search("([0-9]{8})", dst_pfn)
                if self.survey_id is not None:
                    if survey not in self.survey_id.split('/'):
                        continue

                if self.exclude_survey_id is not None:
                    if survey in self.exclude_survey_id.split('/'):
                        continue

                if self.ship_id is not None:
                    if ship.lower() not in [x.lower() for x in self.ship_id.split('/')]:
                        continue

                if self.exclude_ship_id is not None:
                    if ship.lower() in [x.lower() for x in self.exclude_ship_id.split('/')]:
                        continue

                if date is not None:
                    date = date[0]
                    if self.min_year is not None and int(date[:4]) < self.min_year:
                        continue
                
                    if self.max_year is not None and int(date[:4]) > self.max_year:
                        continue

                self.date = date                    
                if survey in these_surveys.keys():
                    mod_url = data_url.split(' ')[0]
                    if version in these_surveys[survey].keys():
                        these_surveys[survey][version].append(
                            [data_url.split(' ')[0],
                             os.path.join(self._outdir, '/'.join([survey, dst_fn])),
                             'mb']
                        )
                    else:
                        these_surveys[survey][version] = [
                            [data_url.split(' ')[0], '/'.join([survey, dst_fn]), 'mb']
                        ]
                        
                else:
                    these_surveys[survey] = {version: [
                        [data_url.split(' ')[0], '/'.join([survey, dst_fn]), 'mb']
                    ]}

        else:
            utils.echo_error_msg(
                'failed to fetch multibeam request'
            )

        with tqdm(
                total=len(these_surveys.keys()),
                desc='scanning NCEI Multibeam datasets',
                leave=self.verbose
        ) as pbar:

            for key in these_surveys.keys():
                pbar.update()
                want_generated = False
                if self.processed_p:
                    if '2' in these_surveys[key].keys():
                        want_generated = self.want_generated(these_surveys[key]['2'][0][0])
                        for v2 in these_surveys[key]['2']:
                            if want_generated:
                                v2_url = v2[0].split('/')
                                v2_url.insert(-1, 'generated')
                                v2_url = '/'.join(v2_url)
                                v2_gen = [v2_url, v2[1], v2[2]]
                                if not self.inf_only:
                                    self.add_entry_to_results(*v2_gen)
                                    
                                inf_url = self.inf_url(v2_gen)
                            else:
                                if not self.inf_only:
                                    self.add_entry_to_results(*v2)
                                    
                                inf_url = self.inf_url(v2)

                            if self.want_inf:
                                self.add_entry_to_results(
                                    '{}.inf'.format(inf_url),
                                    '{}.inf'.format(v2[1]),
                                    'mb_inf'
                                )
                    else:
                        want_generated = self.want_generated(these_surveys[key]['1'][0][0])
                        for v1 in these_surveys[key]['1']:
                            if want_generated:
                                v1_url = v1[0].split('/')
                                v1_url.insert(-1, 'generated')
                                v1_url = '/'.join(v1_url)
                                v1_gen = [v1_url, v1[1], v1[2]]
                                if not self.inf_only:
                                    self.add_entry_to_results(*v1_gen)
                                    
                                inf_url = self.inf_url(v1_gen)
                            else:
                                if not self.inf_only:
                                    self.add_entry_to_results(*v1)
                                    
                                inf_url = self.inf_url(v1)

                            if self.want_inf:
                                #inf_url = self.inf_url(v1)
                                self.add_entry_to_results(
                                    '{}.inf'.format(inf_url),
                                    '{}.inf'.format(v1[1]),
                                    'mb_inf'
                                )
                else:
                    for keys in these_surveys[key].keys():
                        want_generated = self.want_generated(these_surveys[key][keys][0][0])
                        for survs in these_surveys[key][keys]:
                            if want_generated:
                                survs_url = survs[0].split('/')
                                survs_url.insert(-1, 'generated')
                                survs_url = '/'.join(survs_url)
                                survs_gen = [survs_url, survs[1], survs[2]]
                                if not self.inf_only:
                                    self.add_entry_to_results(*survs_gen)
                                    
                                inf_url = self.inf_url(survs_gen)
                            else:
                                if not self.inf_only:
                                    self.add_entry_to_results(*survs)
                                    
                                inf_url = self.inf_url(survs)

                            if self.want_inf:
                                #inf_url = self.inf_url(survs)
                                self.add_entry_to_results(
                                    '{}.inf'.format(inf_url),
                                    '{}.inf'.format(survs[1]),
                                    'mb_inf'
                                )

        if self.make_datalist:
            s_got = []
            with open('mb_inf.txt', 'w') as mb_inf_txt:
                for entry in self.results:
                    try:
                        survey, src_data, mb_fmt, mb_perc, mb_date \
                            = self.parse_entry_inf(entry)
                        if survey in s_got:
                            continue
                        
                        this_year = int(utils.this_year()) \
                            if self.min_year is None \
                               else self.min_year
                        this_weight = float(mb_perc) \
                            * ((int(mb_date)-2000)/(this_year-2000))/100.
                        mb_inf_txt.write(
                            '{} -1 {}\n'.format(survey, this_weight)
                        )
                        s_got.append(survey)
                    except:
                        pass

                    
    def inf_url(self, entry):
        if entry[0][-3:] == 'fbt':
            inf_url = utils.fn_basename2(entry[0])
        else:
            inf_url = entry[0]
        return(inf_url)            

    
    def echo_inf(self, entry):
        print(self.parse_entry_inf(entry))

        
    def parse_entry_inf(self, entry, keep_inf=False, out_dir=None):
        src_data = os.path.basename(entry['dst_fn'])
        if src_data[-3:] == 'fbt':
            src_mb = utils.fn_basename2(src_data)
            inf_url = utils.fn_basename2(entry['url'])
        else:
            inf_url = entry['url']
            src_mb = src_data
            
        survey = entry['url'].split('/')[7]
        src_inf = os.path.join(self._outdir, '{}.inf'.format(entry['dst_fn']))
        try:
            status = fetches.Fetch(
                '{}.inf'.format(inf_url), callback=self.callback, verbose=True
            ).fetch_file(src_inf)
        except:
            utils.echo_warning_msg(
                f'failed to fetch inf file: {inf_url}.inf'
            )
            status = -1
            
        if status == 0:
            mb_fmt = self.mb_inf_data_format(src_inf)
            mb_date = self.mb_inf_data_date(src_inf)
            mb_perc = self.mb_inf_perc_good(src_inf)
            if not keep_inf:
                utils.remove_glob(src_inf)
                
            return(survey, src_data, mb_fmt, mb_perc, mb_date)


class MBDB(fetches.FetchModule):
    """MBDB fetching. This is a test module and does not work."""
    
    def __init__(self, where='1=1', layer=1, list_surveys=False, want_inf=True, **kwargs):
        super().__init__(name='mbdb', **kwargs)
        self.where = where        
        self._mb_dynamic_url = ('https://gis.ngdc.noaa.gov/arcgis/rest/'
                                'services/multibeam_footprints/MapServer')
        self._mb_features_url = ('https://gis.ngdc.noaa.gov/arcgis/rest/'
                                 'services/multibeam_datasets/FeatureServer')
        self._mb_features_products_url = f'{self._mb_features_url}/0'
        self._mb_features_processed_url = f'{self._mb_features_url}/1'
        self._mb_features_raw_url = f'{self._mb_features_url}/2'

        self._mb_dynamic_processed_url = f'{self._mb_dynamic_url}1'

        self._mb_dynamic_query_url = '{0}/{1}/query?'.format(self._mb_dynamic_url, layer)
        self._mb_features_query_url = '{0}/{1}/query?'.format(self._mb_features_url, layer)

        self.list_surveys = list_surveys
        self.want_inf = want_inf

    def inf_parse(self, inf_text):
        this_row = 0
        minmax = [0,0,0,0,0,0]
        for il in inf_text:
            til = il.split()
            if len(til) > 1:
                if til[0] == 'Minimum':
                    if til[1] == 'Longitude:':
                        minmax[0] = utils.float_or(til[2])
                        minmax[1] = utils.float_or(til[5])
                    elif til[1] == 'Latitude:':
                        minmax[2] = utils.float_or(til[2])
                        minmax[3] = utils.float_or(til[5])
                    elif til[1] == 'Depth:':
                        minmax[4] = utils.float_or(til[5]) * -1
                        minmax[5] = utils.float_or(til[2]) * -1

        mbs_region = regions.Region().from_list(minmax)
        return(mbs_region)

        
    def check_inf_region(self, mb_url, keep_inf=False):
        src_mb = mb_url
        inf_url = '{}.inf'.format(utils.fn_basename2(src_mb))
        inf_region = None
        #try:
        #utils.echo_msg(inf_url)
        _req = fetches.Fetch(inf_url, callback=self.callback, verbose=False).fetch_req()
        if _req is None:
            inf_url = '{}.inf'.format(src_mb)
            _req = fetches.Fetch(inf_url, callback=self.callback, verbose=False).fetch_req()
        #except:
        #    utils.echo_warning_msg('failed to fetch inf file: {}'.format(inf_url))
        #    status = -1
            
        #if status == 0:
        # if not keep_inf:
        #     utils.remove_glob(src_inf)
        if _req is not None:
            src_inf = _req.text
            inffile = StringIO(src_inf)
            inffile.read()
            inffile.seek(0)
            inf_region = self.inf_parse(inffile)
            inffile.close()
            
        return(inf_url, inf_region)
        
        
    def run(self):
        """Run the MBDB fetching module"""
        
        if self.region is None:
            return([])

        _data = {
            'where': self.where,
            'outFields': '*',
            'geometry': self.region.format('bbox'),
            'inSR':4326,
            'outSR':4326,
            'f':'pjson',
            'returnGeometry':'false',
        }
        _req = fetches.Fetch(
            self._mb_features_query_url, verbose=self.verbose
        ).fetch_req(params=_data)
        if _req is not None:
            #print(_req.url)
            #print(_req.text)
            features = _req.json()
            with tqdm(
                    total=len(features['features']),
                    desc='parsing mb surveys..',
                    leave=self.verbose
            ) as pbar:
                #for feature in features['data']:
                for feature in features['features']:
                    feature_set_url = feature['attributes']['DOWNLOAD_URL']
                    if feature_set_url is not None:
                        feature_set_id = feature['attributes']['SURVEY_ID']
                        dataset_id = feature['attributes']['DATASET_ID']
                        survey_name = feature['attributes']['SURVEY_NAME']

                        if self.list_surveys:
                            print(feature_set_url)
                        else:
                            page = fetches.Fetch(feature_set_url, verbose=True).fetch_html()
                            page_mb_links = page.xpath(f'//a[contains(@href, "/MB/")]/@href')
                            with tqdm(
                                    total=len(page_mb_links),
                                    desc='parsing mb survey datasets..',
                                    leave=self.verbose
                            ) as pbar_inf:

                                for mb in page_mb_links:
                                    if 'gz' in mb:
                                        mb = '{}.fbt'.format(utils.fn_basename2(mb))                                        
                                        #mb_inf = f'{mb[:-3]}.inf'
                                        
                                    mb_inf, mb_inf_region = self.check_inf_region(mb)
                                    #utils.echo_msg(mb_inf)
                                    #utils.echo_msg(mb_inf_region)
                                    if mb_inf_region is not None:
                                        if regions.regions_intersect_ogr_p(mb_inf_region, self.region):
                                            self.add_entry_to_results(mb, os.path.join(self._outdir, os.path.basename(mb)),'mbs')
                                            
                                            if self.want_inf:
                                                self.add_entry_to_results(
                                                    '{}.inf'.format(mb_inf),
                                                    '{}.inf'.format(os.path.join(self._outdir, os.path.basename(mb_inf))),
                                                    'mb_inf'
                                                )                                            

                                    pbar_inf.update()
                            # mbs = [f'{x[:-3]}.fbt' for x in page.xpath(f'//a[contains(@href, "/MB/")]/@href')]
                            # [self.check_inf_region(x) for x in mbs]
                            # [self.add_entry_to_results(mb, os.path.join(self._outdir, os.path.basename(mb)),'mbs') for mb in mbs]
                            #print(mbs)
                        pbar.update()

            
            # for feature in features['features']:
            #     print(feature)



class R2R(fetches.FetchModule):
    def __init__(self, check_inf=False, **kwargs):
        super().__init__(name='R2R', **kwargs)
        self.check_inf = check_inf
        self.r2r_api_url = 'https://service.rvdata.us/api/fileset/keyword/multibeam?'
        self.r2r_api_manifest_url = 'https://service.rvdata.us/api/file_manifest/?'
        self.r2r_api_product_url = 'https://service.rvdata.us/api/product/?'


    def inf_parse(self, inf_text):
        this_row = 0
        minmax = [0,0,0,0,0,0]
        for il in inf_text:
            til = il.split()
            if len(til) > 1:
                if til[0] == 'Minimum':
                    if til[1] == 'Longitude:':
                        minmax[0] = utils.float_or(til[2])
                        minmax[1] = utils.float_or(til[5])
                    elif til[1] == 'Latitude:':
                        minmax[2] = utils.float_or(til[2])
                        minmax[3] = utils.float_or(til[5])
                    elif til[1] == 'Depth:':
                        minmax[4] = utils.float_or(til[5]) * -1
                        minmax[5] = utils.float_or(til[2]) * -1

        mbs_region = regions.Region().from_list(minmax)
        return(mbs_region)


    def check_inf_region(self, mb_url, keep_inf=False):
        src_mb = mb_url
        inf_url = '{}.inf'.format(utils.fn_basename2(src_mb))
        inf_region = None
        #try:
        #utils.echo_msg(inf_url)
        _req = fetches.Fetch(inf_url, callback=self.callback, verbose=False).fetch_req()
        if _req is None:
            inf_url = '{}.inf'.format(src_mb)
            _req = fetches.Fetch(inf_url, callback=self.callback, verbose=False).fetch_req()
        #except:
        #    utils.echo_warning_msg('failed to fetch inf file: {}'.format(inf_url))
        #    status = -1
            
        #if status == 0:
        # if not keep_inf:
        #     utils.remove_glob(src_inf)
        if _req is not None:
            src_inf = _req.text
            inffile = StringIO(src_inf)
            inffile.read()
            inffile.seek(0)
            inf_region = self.inf_parse(inffile)
            inffile.close()
            
        return(inf_url, inf_region)

    
    def run(self):

        _data = {
            'spatial_bounds': self.region.export_as_wkt(),
        }
        _req = fetches.Fetch(
            self.r2r_api_url, verbose=self.verbose
        ).fetch_req(params=_data)
        if _req is not None:
            features = _req.json()
            #utils.echo_msg(features)
            #utils.echo_msg('found {} cruises'.format(len(features['data'])))
            with tqdm(
                    total=len(features['data']),
                    desc='parsing cruise datasets..',
                    leave=self.verbose
            ) as pbar:
                for feature in features['data']:
                    pbar.update()
                    feature_set_url = feature['download_url']
                    if feature_set_url is not None:
                        feature_set_id = feature['fileset_id']
                        cruise_id = feature['cruise_id']
                        
                        #product_url = '{}cruise_id={}?datatype=Bathymetry'.format(
                        product_url = '{}cruise_id={}'.format(
                            self.r2r_api_product_url, cruise_id
                        )
                        #utils.echo_msg(product_url)
                        page = fetches.Fetch(product_url, verbose=True).fetch_req()
                        page_json = page.json()
                        if page_json is not None:
                            page_data = page_json['data']
                            if page_data is not None:
                                for data in page_data:
                                    #print(data['datatype_name'])
                                    if data['datatype_name'] == 'Bathymetry':
                                        #utils.echo_msg(data)
                                        #utils.echo_msg(page_json)
                                        actual_url = data['actual_url']
                                        if actual_url is not None:
                                            self.add_entry_to_results(
                                                actual_url,
                                                os.path.basename(actual_url),
                                                'r2rBathymetry'
                                            )

### End
