### fetches_factory.py
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

import os, sys
from cudem import utils
from cudem import regions
from cudem import factory
from cudem.fetches import fetches
from cudem.fetches import __version__

## Fetches Module Parser
class FetchesFactory(factory.CUDEMFactory):
    """Acquire a fetches module. Add a new fetches module here to 
    expose it in the  CLI or API via FetchesFactory.
    
    Use the Factory in python by calling: FetchesFactory()"""

    from . import gmrt
    from . import margrav
    from . import srtmplus
    from . import synbath
    from . import charts
    from . import dav
    from . import multibeam
    from . import gebco
    from . import gedtm30
    from . import mgds
    from . import trackline
    from . import ehydro
    from . import ngs
    from . import hydronos
    from . import nceithredds
    from . import etopo
    from . import tnm
    from . import emodnet
    from . import chs
    from . import hrdem
    from . import mrdem
    from . import arcticdem
    from . import bluetopo
    from . import osm
    from . import copernicus
    from . import fabdem
    from . import nasadem
    from . import tides
    from . import vdatum
    from . import buoys
    from . import earthdata
    from . import usiei
    from . import wsf
    from . import hydrolakes
    from . import bingbfp
    from . import waterservices
    from . import csb
    from . import cptcity
    from . import wadnr
    from . import nswtb
    from . import cdse
    
    _modules = {
        'https': {'call': fetches.HttpDataset},
        'gmrt': {'call': gmrt.GMRT},
        'mar_grav': {'call': margrav.MarGrav},
        'srtm_plus': {'call': srtmplus.SRTMPlus},
        'synbath': {'call': synbath.SynBath},
        'charts': {'call': charts.NauticalCharts},
	    'digital_coast': {'call': dav.DAV},
        'SLR': {'call': dav.SLR},
        'CoNED': {'call': dav.CoNED},
        'CUDEM': {'call': dav.CUDEM},
        'multibeam': {'call': multibeam.Multibeam}, 
        'mbdb': {'call': multibeam.MBDB},
        'r2r': {'call': multibeam.R2R},
        'gebco': {'call': gebco.GEBCO},
        'gedtm30': {'call': gedtm30.GEDTM30},
        'mgds': {'call': mgds.MGDS},
        'trackline': {'call': trackline.Trackline},
        'ehydro': {'call': ehydro.eHydro},
        'ngs': {'call': ngs.NGS},
        'hydronos': {'call': hydronos.HydroNOS},
        'ncei_thredds': {'call': nceithredds.NCEIThreddsCatalog},
        'etopo': {'call': etopo.ETOPO},
        'tnm': {'call': tnm.TheNationalMap},
        'ned': {'call': tnm.NED},
        'ned1': {'call': tnm.NED1},
        'tnm_laz': {'call': tnm.TNM_LAZ},
        'emodnet': {'call': emodnet.EMODNet},
        'chs': {'call': chs.CHS},
        'hrdem': {'call': hrdem.HRDEM},
        'mrdem': {'call': mrdem.MRDEM},
        'arcticdem': {'call': arcticdem.ArcticDEM},
        'bluetopo': {'call': bluetopo.BlueTopo},
        'osm': {'call': osm.OpenStreetMap},
        'copernicus': {'call': copernicus.CopernicusDEM},
        'fabdem': {'call': fabdem.FABDEM},
        'nasadem': {'call': nasadem.NASADEM},
        'tides': {'call': tides.Tides},
        'vdatum': {'call': vdatum.VDATUM},
        'buoys': {'call': buoys.BUOYS},
        'earthdata': {'call': earthdata.EarthData},
        'icesat2': {'call': earthdata.IceSat2},
        'mur_sst': {'call': earthdata.MUR_SST},
        'swot': {'call': earthdata.SWOT},
        'usiei': {'call': usiei.USIEI},
        'wsf': {'call': wsf.WSF},
        'hydrolakes': {'call': hydrolakes.HydroLakes},
        'bing_bfp': {'call': bingbfp.BingBFP},
        'waterservices': {'call': waterservices.WaterServices},
        'csb': {'call': csb.CSB},
        'cpt_city': {'call': cptcity.CPTCity},
        'gps_coordinates': {'call': fetches.GPSCoordinates},
        'wa_dnr': {'call': wadnr.waDNR},
        'nsw_tb': {'call': nswtb.NSW_TB},
        'sentinel2': {'call': cdse.Sentinel2},
    }

    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        
## ==============================================
## Command-line Interface (CLI)
## $ fetches
##
## fetches cli
## ==============================================
fetches_usage = lambda: """{cmd} ({version}): Fetches; Fetch and process remote elevation data

usage: {cmd} [ -hlqzAHR [ args ] ] MODULE ...

Options:
  -R, --region\t\tRestrict processing to the desired REGION 
\t\t\tWhere a REGION is xmin/xmax/ymin/ymax[/zmin/zmax[/wmin/wmax]]
\t\t\tUse '-' to indicate no bounding range; e.g. -R -/-/-/-/-10/10/1/-
\t\t\tOR an OGR-compatible vector file with regional polygons. 
\t\t\tWhere the REGION is /path/to/vector[:zmin/zmax[/wmin/wmax]].
\t\t\tIf a vector file is supplied, will use each region found therein.
\t\t\tOptionally, append `:pct_buffer=<value>` to buffer the region(s) by a percentage.
  -H, --threads\t\tSet the number of threads (1)
  -A, --attempts\tSet the number of fetching attempts (5)
  -l, --list\t\tReturn a list of fetch URLs in the given region.
  -z, --no_check_size\tDon't check the size of remote data if local data exists.
  -q, --quiet\t\tLower the verbosity to a quiet

  --modules\t\tDisplay the module descriptions and usage
  --help\t\tPrint the usage text
  --version\t\tPrint the version information

Supported FETCHES modules (see fetches --modules <module-name> for more info): 
  {f_formats}
""".format(cmd=os.path.basename(sys.argv[0]), 
           version=__version__,
           f_formats=utils._cudem_module_short_desc(FetchesFactory._modules))

def fetches_cli(argv = sys.argv):
    """run fetches from command-line

See `fetches_cli_usage` for full cli options.
    """

    i_regions = []
    these_regions = []
    mods = []
    mod_opts = {}
    want_list = False
    want_verbose = True
    stop_threads = False
    check_size = True
    num_threads = 1
    fetch_attempts = 5
    
    ## parse command line arguments.
    i = 1
    while i < len(argv):
        arg = argv[i]
        if arg == '--region' or arg == '-R':
            i_regions.append(str(argv[i + 1]))
            i = i + 1
        elif arg[:2] == '-R':
            i_regions.append(str(arg[2:]))
        elif arg == '-threads' or arg == '--threads' or arg == '-H':
            num_threads = utils.int_or(argv[i + 1], 1)
            i = i + 1
        elif arg[:2] == '-H':
            num_threads = utils.int_or(argv[2:], 1)
        elif arg == '-attempts' or arg == '--attempts' or arg == '-A':
            fetch_attempts = utils.int_or(argv[i + 1], 1)
            i = i + 1
        elif arg[:2] == '-A':
            fetch_attempts = utils.int_or(argv[2:], 1)
        elif arg == '--list' or arg == '-l':
            want_list = True
        elif arg == '--no_check_size' or arg == '-z':
            check_size = False
        elif arg == '--quiet' or arg == '-q':
            want_verbose = False
        elif arg == '--help' or arg == '-h':
            sys.stderr.write(fetches_usage())
            sys.exit(1)
        elif arg == '--version' or arg == '-v':
            print('{}, version {}'.format(
                os.path.basename(sys.argv[0]), __version__)
                  )
            sys.exit(1)
        elif arg == '--modules' or arg == '-m':
            utils.echo_modules(
                FetchesFactory._modules,
                None if i+1 >= len(argv) else sys.argv[i+1]
            )
            sys.exit(0)
        elif arg[0] == '-':
            sys.stderr.write(fetches_usage())
            sys.exit(0)
        else:
            mods.append(arg)
            
        i = i + 1

    if len(mods) == 0:
        sys.stderr.write(fetches_usage())
        utils.echo_error_msg('you must select at least one fetch module')
        sys.exit(-1)

    these_regions = regions.parse_cli_region(i_regions, want_verbose)
    if not these_regions:
        these_regions = [regions.Region().from_string('-R-180/180/-90/90')]
        
    for rn, this_region in enumerate(these_regions):
        if stop_threads:
            return
        
        x_fs = [
            FetchesFactory(
                mod=mod,
                src_region=this_region,
                verbose=want_verbose
            )._acquire_module() for mod in mods
        ]
        for x_f in x_fs:
            if x_f is None:
                continue
            
            if want_verbose:
                utils.echo_msg(
                    'running fetch module {} on region {}...'.format(
                        x_f.name, this_region.format('str')
                    )
                )

            try:
                x_f.run()
            except (KeyboardInterrupt, SystemExit):
                utils.echo_error_msg(
                    'user breakage...please wait while fetches exits.'
                )
                sys.exit(-1)
                
            if want_verbose:
                utils.echo_msg(
                    'found {} data files.'.format(len(x_f.results))
                )
                
            if len(x_f.results) == 0:
                break
            
            if want_list:
                for result in x_f.results:
                    print(result['url'])
            else:
                try:
                    fr = fetches.fetch_results(
                        x_f,
                        n_threads=num_threads,
                        check_size=check_size,
                        attempts=fetch_attempts
                    )
                    fr.daemon = True                
                    fr.start()
                    fr.join()         
                except (KeyboardInterrupt, SystemExit):
                    utils.echo_error_msg(
                        'user breakage...please wait while fetches exits.'
                    )
                    x_f.status = -1
                    stop_threads = True
                    while not fr.fetch_q.empty():
                        try:
                            fr.fetch_q.get(False)
                        except Empty:
                            continue
                        
                        fr.fetch_q.task_done()                        

### End
