### fetches.py
##
## Copyright (c) 2010 - 2023 Regents of the University of Colorado
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
##
## Fetch elevation and related data from a variety of sources.
## Use CLI command 'fetches'
## or use FetchesFactory() to acquire and use a fetch module.
##
### Code:

import os
import sys
import time

from cudem import utils
from cudem import regions
#from cudem import datasets
from cudem import fetches

import cudem.fetches.utils as f_utils
import cudem.fetches.multibeam as mb
import cudem.fetches.ehydro as ehydro
import cudem.fetches.gmrt as gmrt
import cudem.fetches.srtm as srtm
import cudem.fetches.mar_grav as mar_grav
import cudem.fetches.ngs as ngs
import cudem.fetches.nos as nos
import cudem.fetches.charts as charts
import cudem.fetches.digital_coast as digital_coast
import cudem.fetches.ncei_thredds as ncei_thredds
import cudem.fetches.tnm as tnm
import cudem.fetches.emodnet as emodnet
import cudem.fetches.chs as chs
import cudem.fetches.hrdem as hrdem
import cudem.fetches.osm as osm
import cudem.fetches.globalelus as globalelus
import cudem.fetches.copernicus as copernicus
import cudem.fetches.nasadem as nasadem
import cudem.fetches.tides as tides
import cudem.fetches.vdatum as vdatum
import cudem.fetches.buoys as buoys
import cudem.fetches.earthdata as earthdata
import cudem.fetches.usiei as usiei
import cudem.fetches.trackline as trackline
import cudem.fetches.mgds as mgds
import cudem.fetches.arcticdem as arcticdem
import cudem.fetches.bluetopo as bluetopo
import cudem.fetches.hydrolakes as hydrolakes
import cudem.fetches.gebco as gebco
import cudem.fetches.wsf as wsf
import cudem.fetches.fabdem as fabdem

## ==============================================
## Fetches Module Parser
## ==============================================
class FetchesFactory:
    """Acquire a fetches module."""
    
    mods = {
        'gmrt': {
            'class': gmrt.GMRT,
        },
        'gebco': {
            'class': gebco.GEBCO,
        },
        'srtm_plus': {
            'class': srtm.SRTMPlus,
        },
        'mar_grav': {
            'class': mar_grav.MarGrav,
        },
        'multibeam': {
            'class': mb.Multibeam,
        },
        'mgds': {
            'class': mgds.MGDS,
        },
        'trackline': {
            'class': trackline.Trackline,
        },
        'ehydro': {
            'class': ehydro.eHydro,
        },
        'ngs': {
            'class': ngs.NGS,
        },
        'nos': {
            'class': nos.HydroNOS,
        },
        'charts': {
            'class': charts.NauticalCharts,
        },
        'digital_coast': {
            'class': digital_coast.DAV,
        },
        'ncei_thredds': {
            'class': ncei_thredds.NCEIThreddsCatalog,
        },
        'tnm': {
            'class': tnm.TheNationalMap,
        },
        'emodnet': {
            'class': emodnet.EMODNet,
        },
        'chs': {
            'class': chs.CHS,
        },
        'hrdem': {
            'class': hrdem.HRDEM,
        },
        'arcticdem': {
            'class': arcticdem.ArcticDEM,
        },

        'bluetopo': {
            'class': bluetopo.BlueTopo,
        },
        'osm': {
            'class': osm.OpenStreetMap,
        },
        'copernicus': {
            'class': copernicus.CopernicusDEM,
        },
        'fabdem': {
            'class': fabdem.FABDEM,
        },
        'nasadem': {
            'class': nasadem.NASADEM,
        },
        'tides': {
            'class': tides.Tides,
        },
        'vdatum': {
            'class': vdatum.VDATUM,
        },
        'buoys': {
            'class': buoys.BUOYS,
        },
        'earthdata': {
            'class': earthdata.EarthData,
        },
        'usiei': {
            'class': usiei.USIEI,
        },
        'wsf': {
            'class': wsf.WSF,
        },
        'hydrolakes': {
            'class': hydrolakes.HydroLakes,
        },
    }
    
    def __init__(
            self,
            mod=None,
            src_region=None,
            callback=lambda: False,
            weight=None,
            dst_srs=None,
            x_inc=None,
            y_inc=None,
            outdir=None,
            verbose=True
    ):
        self.mod = mod
        self.mod_args = {}
        self.region = src_region
        self.callback = callback
        self.weight = weight
        self.verbose = verbose
        self.status = 0
        self.dst_srs = dst_srs
        self.x_inc = utils.str2inc(x_inc)
        self.y_inc = utils.str2inc(y_inc)
        self.outdir = outdir
        self.results = []
        if self.mod is not None:
            this_mod = self.parse_mod()
        
    def parse_mod(self):
        opts = self.mod.split(':')
        if opts[0] in FetchesFactory.mods.keys():
            if len(opts) > 1:
                self.mod_args = utils.args2dict(list(opts[1:]), {})
            self.mod = opts[0]
        else:
            utils.echo_error_msg('could not parse fetch module `{}`'.format(opts[0]))
            return(None)
        
        return(self)

    def add_module(self, type_def={}):
        for key in type_def.keys():
            self.mods[key] = type_def[key]

    def acquire(self, **kwargs):
        """Acquire the fetches module self.mod"""

        if self.mod in self.mods.keys():
            return(
                self.mods[self.mod]['class'](
                    src_region=self.region,
                    callback=self.callback,
                    weight=self.weight,
                    dst_srs=self.dst_srs,
                    x_inc=self.x_inc,
                    y_inc=self.y_inc,
                    verbose=self.verbose,
                    outdir=self.outdir,
                    **kwargs,
                    **self.mod_args
                )
            )
        else:
            return(None)

## ==============================================
## Command-line Interface (CLI)
## $ fetches
##
## fetches cli
## ==============================================
fetches_usage = """{cmd} ({f_version}): Fetches; Fetch and process remote elevation data

usage: {cmd} [ -hlmpqERW [ args ] ] MODULE ...

Options:
  -R, --region\t\tRestrict processing to the desired REGION 
\t\t\tWhere a REGION is xmin/xmax/ymin/ymax[/zmin/zmax[/wmin/wmax]]
\t\t\tUse '-' to indicate no bounding range; e.g. -R -/-/-/-/-10/10/1/-
\t\t\tOR an OGR-compatible vector file with regional polygons. 
\t\t\tWhere the REGION is /path/to/vector[:zmin/zmax[/wmin/wmax]].
\t\t\tIf a vector file is supplied, will use each region found therein.
  -W, --t_srs\t\tSet the TARGET projection (for use with `-p`).
  -E, --increment\tBlock data to INCREMENT in native units.
\t\t\tWhere INCREMENT is x-inc[/y-inc] (for use with `-p`).
  -H, --threads\t\tSet the number of threads (1)
  -l, --list\t\tReturn a list of fetch URLs in the given region.
  -p, --process\t\tProcess fetched elevation data to ASCII XYZ format. <beta>
  -z, --no_check_size\tDon't check the size of remote data if local data exists.
  -q, --quiet\t\tLower the verbosity to a quiet

  --modules\t\tDisplay the module descriptions and usage
  --help\t\tPrint the usage text
  --version\t\tPrint the version information

Supported FETCHES modules (see fetches --modules <module-name> for more info): 
  {f_formats}

CIRES DEM home page: <http://ciresgroups.colorado.edu/coastalDEM>
""".format(cmd=os.path.basename(sys.argv[0]), 
           f_version=fetches.__version__,
           f_formats=utils._cudem_module_short_desc(FetchesFactory.mods))

def fetches_cli(argv = sys.argv):
    """run fetches from command-line

See `fetches_cli_usage` for full cli options.
    """

    i_regions = []
    these_regions = []
    mods = []
    mod_opts = {}
    want_list = False
    want_proc = False
    want_verbose = True
    stop_threads = False
    check_size = True
    num_threads = 1
    dst_srs = None
    xy_inc = [None, None]
    
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
        elif arg == '--increment' or arg == '-E':
            xy_inc = argv[i + 1].split('/')
            i = i + 1
        elif arg[:2] == '-E':
            xy_inc = arg[2:].split('/')
        elif arg == '-t_srs' or arg == '--t_srs' or arg == '-W':
            dst_srs = argv[i + 1]
            i = i + 1
        elif arg == '-threads' or arg == '--threads' or arg == '-H':
            num_threads = utils.int_or(argv[i + 1], 1)
            i = i + 1
        elif arg == '--list' or arg == '-l':
            want_list = True
        elif arg == '--process' or arg == '-p':
            want_proc = True
        elif arg == '--no_check_size' or arg == '-z':
            check_size = False
        elif arg == '--quiet' or arg == '-q':
            want_verbose = False
        elif arg == '--help' or arg == '-h':
            sys.stderr.write(fetches_usage)
            sys.exit(1)
        elif arg == '--version' or arg == '-v':
            print('{}, version {}'.format(
                os.path.basename(sys.argv[0]), fetches.__version__)
                  )
            sys.exit(1)
        elif arg == '--modules' or arg == '-m':
            utils.echo_modules(FetchesFactory.mods, None if i+1 >= len(argv) else sys.argv[i+1])
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

    if len(xy_inc) < 2:
        xy_inc.append(xy_inc[0])
        
    elif len(xy_inc) == 0:
        xy_inc = [None, None]
        
    these_regions = regions.parse_cli_region(i_regions, want_verbose)
    if not these_regions:
        these_regions = [regions.Region().from_string('-R-180/180/-90/90')]

    if want_proc:
        if None in xy_inc:
            utils.echo_error_msg('You must supply an increment -E in order to process this dataset')
            sys.exit(-1)
        else:
            if want_verbose:
                utils.echo_msg('Process fetched data to region(s) {} @ {}x{}'.format(these_regions, xy_inc[0], xy_inc[1]))
        
    for rn, this_region in enumerate(these_regions):
        if stop_threads:
            return
        
        x_fs = [FetchesFactory(
            mod=mod,
            src_region=this_region,
            dst_srs=dst_srs,
            x_inc=xy_inc[0],
            y_inc=xy_inc[1],
            verbose=want_verbose
        ).acquire() for mod in mods]
        for x_f in x_fs:
            if x_f is None:
                continue
            
            if want_verbose:
                utils.echo_msg(
                    'running fetch module {} on region {}...'.format(
                        x_f.name, this_region.format('str')
                    )
                )
            
            x_f.run()
            if want_verbose:
                utils.echo_msg(
                    'found {} data files.'.format(len(x_f.results))
                )
                
            if len(x_f.results) == 0:
                break
            
            if want_list:
                for result in x_f.results:
                    print(result[0])
            else:
                fr = f_utils.fetch_results(x_f, want_proc=want_proc, n_threads=num_threads, check_size=check_size)
                fr.daemon = True
                try:
                    fr.start()
                except (KeyboardInterrupt, SystemExit):
                    utils.echo_error_msg('user breakage...please wait for while fetches exits.')
                    x_f.status = -1
                    stop_threads = True
                    while not fr.fetch_q.empty():
                        try:
                            fr.fetch_q.get(False)
                        except Empty:
                            continue
                        
                        fr.fetch_q.task_done()
                        
                fr.join()
                
### End
