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
import time
from cudem import utils
from cudem import regions
from cudem import datasets

import cudem.fetches.utils as f_utils
import cudem.fetches.multibeam as mb
import cudem.fetches.usace as usace
import cudem.fetches.gmrt as gmrt
import cudem.fetches.srtm as srtm
import cudem.fetches.mar_grav as mar_grav
import cudem.fetches.ngs as ngs
import cudem.fetches.nos as nos
        
## ==============================================
## Fetches Module Parser
## ==============================================
class FetchesFactory:

    mods = {
        'gmrt': {'description': """The Global Multi-Resolution Topography (GMRT) synthesis is a multi-resolutional 
compilation of edited multibeam sonar data collected by scientists and institutions worldwide, that is 
reviewed, processed and gridded by the GMRT Team and merged into a single continuously updated compilation 
of global elevation data. The synthesis began in 1992 as the Ridge Multibeam Synthesis (RMBS), was expanded 
to include multibeam bathymetry data from the Southern Ocean, and now includes bathymetry from throughout 
the global and coastal oceans."""},
        'srtm_plus': {'description': """"""},
        'mar_grav': {'description': """"""},
        'multibeam': {'description': """"""},
        'usace': {'description': """"""},
        'ngs': {'description': """"""},
        'nos': {'description': """"""},
    }
    
    def __init__(self, mod=None, src_region=None, callback=lambda: False, weight=None, verbose=True):
        self.mod = mod
        self.mod_args = {}
        self.region = src_region
        self.callback = callback
        self.weight = weight
        self.verbose = verbose
        self.status = 0
        self.results = []
        if self.mod is not None:
            self.parse_mod()
            
    def parse_mod(self):
        opts = self.mod.split(':')
        if opts[0] in FetchesFactory.mods.keys():
            if len(opts) > 1:
                self.mod_args = utils.args2dict(list(opts[1:]))
            self.mod = opts[0]
        else:
            utils.echo_error_msg('invalid module name `{}`'.format(opts[0]))
        return(self)

    def add_module(self, type_def={}):
        for key in type_def.keys():
            self.mods[key] = type_def[key]

    def acquire_mb(self, **kwargs):
        return(mb.Multibeam(
            src_region=self.region, callback=self.callback, weight=self.weight, verbose=self.verbose, **kwargs, **self.mod_args))

    def acquire_usace(self, **kwargs):
        return(usace.USACE(
            src_region=self.region, callback=self.callback, weight=self.weight, verbose=self.verbose, **kwargs, **self.mod_args))

    def acquire_gmrt(self, **kwargs):
        return(gmrt.GMRT(
            src_region=self.region, callback=self.callback, weight=self.weight, verbose=self.verbose, **kwargs, **self.mod_args))

    def acquire_srtm_plus(self, **kwargs):
        return(srtm.SRTMPlus(
            src_region=self.region, callback=self.callback, weight=self.weight, verbose=self.verbose, **kwargs, **self.mod_args))

    def acquire_mar_grav(self, **kwargs):
        return(mar_grav.MarGrav(
            src_region=self.region, callback=self.callback, weight=self.weight, verbose=self.verbose, **kwargs, **self.mod_args))

    def acquire_ngs(self, **kwargs):
        return(ngs.NGS(
            src_region=self.region, callback=self.callback, weight=self.weight, verbose=self.verbose, **kwargs, **self.mod_args))

    def acquire_nos(self, **kwargs):
        return(nos.NOS(
            src_region=self.region, callback=self.callback, weight=self.weight, verbose=self.verbose, **kwargs, **self.mod_args))

    def acquire(self, **kwargs):
        if self.mod == 'mb':
            return(self.acquire_mb(**kwargs))

        if self.mod == 'usace':
            return(self.acquire_usace(**kwargs))
        
        if self.mod == 'gmrt':
            return(self.acquire_gmrt(**kwargs))

        if self.mod == 'srtm_plus':
            return(self.acquire_srtm_plus(**kwargs))
       
        if self.mod == 'mar_grav':
            return(self.acquire_mar_grav(**kwargs))
        
        if self.mod == 'ngs':
            return(self.acquire_ngs(**kwargs))

        if self.mod == 'nos':
            return(self.acquire_nos(**kwargs))

class Fetcher(datasets.XYZDataset):
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.fetch_module = FetchesFactory(mod=self.fn, src_region=self.region, verbose=self.verbose).acquire(warp=4326)
        
    def generate_inf(self):
        """generate a infos dictionary from the multibeam dataset

        Returns:
          dict: a data-entry infos dictionary
        """

        self.infos['name'] = self.fn
        self.infos['hash'] = None
        self.infos['minmax'] = self.region.export_as_list()
        self.infos['numpts'] = 0
        self.infos['wkt'] = self.region.export_as_wkt()
        return(self.infos)

    def yield_xyz(self):
        for xyz in self.fetch_module.yield_results_to_xyz():
            yield(xyz)
        
_fetches_module_short_desc = lambda: ', '.join(
    ['{}'.format(key) for key in FetchesFactory().mods])
_fetches_module_long_desc = lambda x: 'fetches modules:\n% fetches ... <mod>:key=val:key=val...\n\n  ' + '\n  '.join(
    ['\033[1m{:14}\033[0m{}\n'.format(key, x[key]['description']) for key in x]) + '\n'

## ==============================================
## Command-line Interface (CLI)
## $ fetches
##
## fetches cli
## ==============================================
fetches_usage = """{cmd} ({f_version}): Fetches; Fetch and process remote elevation data

usage: {cmd} [ -hiqwPRW [ args ] ] MODULE ...

Options:
  -R, --region\t\tRestrict processing to the desired REGION 
\t\t\tWhere a REGION is [ xmin/xmax/ymin/ymax/[ zmin/zmax/[ wmin/wmax ] ] ]
\t\t\tUse '-' to indicate no bounding range; e.g. -R -/-/-/-/-10/10/1/-
\t\t\tOR an OGR-compatible vector file with regional polygons. 
\t\t\tIf a vector file is supplied, will use each region found therein.
\t\t\tAppend :zmin/zmax/[ wmin/wmax ] to the file path to extended REGION.
  -l, --list\t\tReturn a list of fetch URLs in the given region.
  -p, --process\t\tProcess fetched elevation data to ASCII XYZ format in WGS84. <beta>

  --weights\t\tOutput WEIGHT values along with xyz
  --quiet\t\tLower the verbosity to a quiet

  --help\t\tPrint the usage text
  --version\t\tPrint the version information

Supported FETCHES modules: 
  {f_formats}

CIRES DEM home page: <http://ciresgroups.colorado.edu/coastalDEM>
""".format(cmd=os.path.basename(sys.argv[0]), 
           f_version=f_utils.__version__,
           f_formats=_fetches_module_short_desc())

def fetches_cli(argv = sys.argv):
    """run fetches from command-line

    See `fetches_cli_usage` for full cli options.
    """

    i_regions = []
    these_regions = []
    mods = []
    mod_opts = {}
    want_weights = False
    want_list = False
    want_proc = False
    want_verbose = True
    stop_threads = False
    
    ## ==============================================
    ## parse command line arguments.
    ## ==============================================
    i = 1
    while i < len(argv):
        arg = argv[i]
        if arg == '--region' or arg == '-R':
            i_regions.append(str(argv[i + 1]))
            i = i + 1
        elif arg[:2] == '-R':
            i_regions.append(str(arg[2:]))
        elif arg == '--weights' or arg == '-w':
            want_weights = True
        elif arg == '--list' or arg == '-l':
            want_list = True
        elif arg == '--process' or arg == '-p':
            want_proc = True
        elif arg == '--quiet' or arg == '-q':
            want_verbose = False
        elif arg == '--help' or arg == '-h':
            sys.stderr.write(fetches_usage)
            sys.exit(1)
        elif arg == '--version' or arg == '-v':
            print('{}, version {}'.format(os.path.basename(sys.argv[0]), __version__))
            sys.exit(1)
        elif arg == '--modules' or arg == '-m':
            try:
                if argv[i + 1] in FetchesFactory.mods.keys():
                    sys.stderr.write(_fetches_module_long_desc({k: FetchesFactory.mods[k] for k in (argv[i + 1],)}))
                else: sys.stderr.write(_fetches_module_long_desc(FetchesFactory.mods))
            except: sys.stderr.write(_fetches_module_long_desc(FetchesFactory.mods))
            sys.exit(0)
        elif arg[0] == '-':
            sys.stderr.write(fetches_usage)
            sys.exit(0)
        else:
            mods.append(arg)
            
        i = i + 1

    if len(mods) == 0:
        sys.stderr.write(fetches_usage)
        utils.echo_error_msg('you must select at least one fetch module')
        sys.exit(-1)
        
    for i_region in i_regions:
        tmp_region = regions.Region().from_string(i_region)
        if tmp_region.valid_p():
            these_regions.append(tmp_region)
        else:
            tmp_region = regions.ogr_wkts(i_region)
            for i in tmp_region:
                if i.valid_p():
                    these_regions.append(i)
                    
    if len(these_regions) == 0:
        these_regions = [None]
    else:
        if want_verbose:
            utils.echo_msg('parsed {} region(s)'.format(len(these_regions)))

    for rn, this_region in enumerate(these_regions):
        if stop_threads: return
        x_fs = [FetchesFactory(mod=mod, src_region=this_region, verbose=want_verbose).acquire(warp=4326) for mod in mods]
        for x_f in x_fs:
            utils.echo_msg('running fetch module {} on region {}...\
            '.format(x_f.name, this_region.format('str')))
            
            #r = x_f.run().results
            x_f.run()
            utils.echo_msg('found {} data files.'.format(len(x_f.results)))
            if want_list:
                for result in x_f.results:
                    print(result[0])
            else:
                fr = f_utils.fetch_results(x_f, want_proc=want_proc)
                fr.daemon = True
                _p = utils.CliProgress('fetching {} remote data files'.format(len(x_f.results)))
                try:
                    fr.start()
                    while True:
                        time.sleep(2)
                        sys.stderr.write('\x1b[2K\r')
                        perc = float((len(x_f.results) - fr.fetch_q.qsize())) / len(x_f.results) * 100 if len(x_f.results) > 0 else 1
                        #if want_verbose: sys.stderr.write('fetches: fetching remote data files [{}%]'.format(perc))
                        if want_verbose:
                            _p.update_perc((len(x_f.results) - fr.fetch_q.qsize(), len(x_f.results)))
                        sys.stderr.flush()
                        if not fr.is_alive():
                            break
                except (KeyboardInterrupt, SystemExit):
                    utils.echo_error_msg('user breakage...please wait for while fetches exits.')
                    x_f.status = -1
                    stop_threads = True
                    while not fr.fetch_q.empty():
                        try:
                            fr.fetch_q.get(False)
                        except Empty: continue
                        fr.fetch_q.task_done()
                fr.join()
                _p.end(x_f.status, 'fetched {} remote data files'.format(len(x_f.results)))
            if want_verbose:
                utils.echo_msg('ran fetch module {} on region {}...\
            '.format(x_f.name, this_region.format('str')))
### End