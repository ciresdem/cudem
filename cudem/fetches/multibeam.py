### multibeam.py - NCEI Multibeam
##
## Copyright (c) 2010 - 2022 Regents of the University of Colorado
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
## MB Fetch
##
## Fetch Multibeam bathymetric surveys from NOAA
## MBSystem is required to process the resulting data
##
## NCEI is the U.S. national archive for multibeam bathymetric data and holds more than 9 million 
## nautical miles of ship trackline data recorded from over 2400 cruises and received from sources 
## worldwide.
##
## Uses NCEI multibeam groovy script to discover multibeam surveys.
##
### Code:

import os

from cudem import utils
from cudem import regions
from cudem import datasets
from cudem import xyzfun

import cudem.fetches.utils as f_utils

## ==============================================
## MapServer testing
## ==============================================
class MBDB(f_utils.FetchModule):
    """NOSHydro"""
    
    def __init__(self, where='1=1', **kwargs):
        super().__init__(**kwargs)
        self._mb_dynamic_url = 'https://gis.ngdc.noaa.gov/arcgis/rest/services/web_mercator/multibeam_dynamic/MapServer/0'
        self._mb_url = 'https://gis.ngdc.noaa.gov/arcgis/rest/services/web_mercator/multibeam/MapServer/0'
        #self._nos_data_url = 'https://data.ngdc.noaa.gov/platforms/ocean/nos/coast/'
        self._mb_query_url = '{0}/query?'.format(self._mb_dynamic_url)
        self._outdir = os.path.join(os.getcwd(), 'multibeam')
        self.name = 'multibeam'
        self.where = where
        
    def run(self):
        if self.region is None:
            return([])

        _data = {
            'where': self.where,
            'outFields': '*',
            'geometry': self.region.format('bbox'),
            'inSR':4326,
            'outSR':4326,
            'f':'pjson',
            'returnGeometry':'False',
        }
        _req = f_utils.Fetch(self._mb_query_url, verbose=self.verbose).fetch_req(params=_data)
        if _req is not None:
            print(_req.text)
            features = _req.json()
            for feature in features['features']:
                pass

class Multibeam(f_utils.FetchModule):
    """Fetch multibeam data from NOAA NCEI"""
    
    def __init__(self, processed=False, inc=None, process=False, **kwargs):
        super().__init__(**kwargs)
        self._mb_data_url = "https://data.ngdc.noaa.gov/platforms/"
        self._mb_metadata_url = "https://data.noaa.gov/waf/NOAA/NESDIS/NGDC/MGG/Multibeam/iso/"
        self._mb_search_url = "https://maps.ngdc.noaa.gov/mapviewer-support/multibeam/files.groovy?"
        self._mb_autogrid = "https://www.ngdc.noaa.gov/maps/autogrid/"
        self._mb_html = "https://www.ngdc.noaa.gov/"
        self._outdir = os.path.join(os.getcwd(), 'mb')
        self._urls = [self._mb_data_url, self._mb_metadata_url, self._mb_autogrid]
        self.name = 'multibeam'
        self.processed_p = processed
        self.process = process
        self.inc = utils.str2inc(inc)

    def mb_inf_data_format(self, src_inf):
        """extract the data format from the mbsystem inf file."""

        with open(src_inf, errors='ignore') as iob:
            for il in iob:
                til = il.split()
                if len(til) > 1:
                    if til[0] == 'MBIO':
                        return(til[4])

    def mb_inf_data_date(self, src_inf):
        """extract the data format from the mbsystem inf file."""

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
                    
    def run(self):
        these_surveys = {}
        these_versions = {}
        if self.region is None: return([])
        _req = f_utils.Fetch(self._mb_search_url).fetch_req(params={'geometry': self.region.format('bbox')}, timeout=20)
        if _req is not None and _req.status_code == 200:
            survey_list = _req.text.split('\n')[:-1]
            for r in survey_list:
                dst_pfn = r.split(' ')[0]
                dst_fn = dst_pfn.split('/')[-1:][0]
                survey = dst_pfn.split('/')[6]
                dn = r.split(' ')[0].split('/')[:-1]
                version = dst_pfn.split('/')[9][-1]
                data_url = self._mb_data_url + '/'.join(r.split('/')[3:])
                if survey in these_surveys.keys():
                    if version in these_surveys[survey].keys():
                        these_surveys[survey][version].append([data_url.split(' ')[0], '/'.join([survey, dst_fn]), 'mb'])
                    else:
                        these_surveys[survey][version] = [[data_url.split(' ')[0], '/'.join([survey, dst_fn]), 'mb']]
                        
                else:
                    these_surveys[survey] = {version: [[data_url.split(' ')[0], '/'.join([survey, dst_fn]), 'mb']]}
                    
        else: utils.echo_error_msg('{}'.format(_req.reason))
                    
        for key in these_surveys.keys():
            if self.processed_p:
                if '2' in these_surveys[key].keys():
                    for v2 in these_surveys[key]['2']:
                        self.results.append(v2)
                else:
                    for v1 in these_surveys[key]['1']:
                        self.results.append(v1)
                        
            else:
                for keys in these_surveys[key].keys():
                    for survs in these_surveys[key][keys]:
                        self.results.append(survs)

        with open('mb_inf.txt', 'w') as mb_inf_txt:
            for entry in self.results:
                mb_inf_txt.write(self.parse_entry_inf(entry))
                mb_inf_txt.write('\n')
                #self.echo_inf(entry)

    def echo_inf(self, entry):
        print(self.parse_entry_inf(entry))
        
    def parse_entry_inf(self, entry, keep_inf=False):
        src_data = os.path.basename(entry[1])
        src_mb = src_data[:-4]
        survey = entry[0].split('/')[7]
        if f_utils.Fetch('{}.inf'.format(entry[0][:-4]), callback=self.callback, verbose=self.verbose).fetch_file('{}.inf'.format(src_mb)) == 0:
            mb_fmt = self.mb_inf_data_format('{}.inf'.format(src_mb))
            mb_date = self.mb_inf_data_date('{}.inf'.format(src_mb))
            mb_perc = self.mb_inf_perc_good('{}.inf'.format(src_mb))
            if not keep_inf:
                utils.remove_glob('{}.inf'.format(src_mb))
            return(survey, src_data, mb_fmt, mb_perc, mb_date)
            
    def yield_xyz(self, entry):
        src_data = os.path.basename(entry[1])
        src_mb = src_data[:-4]
        try:
            survey, src_data, mb_fmt, mb_perc, mb_date = self.parse_entry_inf(entry)
        except TypeError:
            return
        this_inf = self.parse_entry_inf(entry)

        if f_utils.Fetch(entry[0], callback=self.callback, verbose=self.verbose).fetch_file(src_data) == 0:
            src_xyz = os.path.basename(src_data) + '.xyz'
            if not self.process:
                this_weight = self.weight
                out, status = utils.run_cmd('mblist -OXYZ -I{} -Ma > {}'.format(src_data, src_xyz), verbose=True)
            else:
                this_weight = (float(mb_perc) * (1 + (2*((int(mb_date)-2015)/100))))/100.
                out, status = utils.run_cmd('mblist -OXYZ -I{} -MX{} > {}'.format(src_data, str(100-float(mb_perc)), src_xyz), verbose=True)
            
                if status != 0:
                    if f_utils.Fetch('{}.inf'.format(entry[0]), callback=self.callback, verbose=self.verbose).fetch_file('{}.inf'.format(src_mb)) == 0:
                        mb_fmt = self.mb_inf_data_format('{}.inf'.format(src_mb))
                        mb_date = self.mb_inf_data_date('{}.inf'.format(src_mb))
                        out, status = utils.run_cmd('mblist -F{} -OXYZ -I{} -MX{}  > {}'.format(mb_fmt, src_data, str(100-float(mb_perc)), src_xyz), verbose=True)
                    
            if status == 0:
                _ds = datasets.XYZFile(
                    fn=src_xyz,
                    delim='\t',
                    data_format=168,
                    src_srs='epsg:4326',
                    dst_srs=self.dst_srs,
                    #name=os.path.basename(entry[1]),
                    src_region=self.region,
                    verbose=self.verbose,
                    weight=this_weight,
                    remote=True
                )

                if self.inc is not None:
                    xyz_func = lambda p: _ds.dump_xyz(dst_port=p, encode=True)
                    for xyz in utils.yield_cmd(
                            'gmt blockmedian -I{:.10f} {} -r -V'.format(
                                self.inc, self.region.format('gmt')
                            ),
                            verbose=self.verbose,
                            data_fun=xyz_func
                    ):
                        yield(xyzfun.XYZPoint().from_list([float(x) for x in xyz.split()]))
                else:
                    for xyz in _ds.yield_xyz():
                        yield(xyz)

                utils.remove_glob(src_data, '{}*'.format(src_xyz), '{}*.inf'.format(src_mb))
            else:
                utils.echo_error_msg('failed to process local file, {} [{}]...'.format(src_data, entry[0]))
                with open(
                        '{}'.format(os.path.join(self._outdir, 'fetch_{}_{}.err'.format(self.name, self.region.format('fn')))),
                        'a'
                ) as mb_err:
                    mb_err.write('{}\n'.format(','.join([src_mb, entry[0]])))
                    
                os.rename(src_data, os.path.join(self._outdir, src_data))
                utils.remove_glob(src_xyz)
        else:
            utils.echo_error_msg(
                'failed to fetch remote file, {}...'.format(src_data)
            )

### End
