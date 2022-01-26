### fetches.py
##
## Copyright (c) 2010 - 2022 Regents of the University of Colorado
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

## ==============================================
## Fetches Module Parser
## ==============================================
class FetchesFactory:
    """Acquire a fetches module."""
    
    mods = {
        'gmrt': {
            'class': gmrt.GMRT,
            'description': """The Global Multi-Resolution Topography synthesis.
The Global Multi-Resolution Topography (GMRT) synthesis is a multi-resolutional 
compilation of edited multibeam sonar data collected by scientists and institutions worldwide, that is 
reviewed, processed and gridded by the GMRT Team and merged into a single continuously updated compilation 
of global elevation data. The synthesis began in 1992 as the Ridge Multibeam Synthesis (RMBS), was expanded 
to include multibeam bathymetry data from the Southern Ocean, and now includes bathymetry from throughout 
the global and coastal oceans.

layers: 'topo' or 'topo-mask'

< gmrt:res=max:fmt=geotiff:bathy_only=False:layer=topo >""",
        },
        'srtm_plus': {
            'class': srtm.SRTMPlus,
            'description': """SRTM15+: GLOBAL BATHYMETRY AND TOPOGRAPHY AT 15 ARCSECONDS.

< srtm_plus >""",
        },
        'mar_grav': {
            'class': mar_grav.MarGrav,
            'description': """MARine GRAVity Sattelite Altimetry Topography from Scripps.

< mar_grav >""",
        },
        'multibeam': {
            'class': mb.Multibeam,
            'description': """NOAA MULTIBEAM bathymetric data.
NCEI is the U.S. national archive for multibeam bathymetric data and holds more than 9 million 
nautical miles of ship trackline data recorded from over 2400 cruises and received from sources 
worldwide.

< multibeam:inc=None:process=False:min_year=None:survey_id=None:exclude=None >""",
        },
        'trackline': {
            'class': trackline.Trackline,
            'description': """NOAA TRACKLINE bathymetric data.

< trackline >""",
        },
        'ehydro': {
            'class': ehydro.eHydro,
            'description': """USACE eHydro bathymetric data.
Maintenance responsibility for more than 25,000 miles of navigation channels and 400 ports and 
harbors throughout the United States requires extensive surveying and mapping services, including 
boundary, topographic, hydrographic, terrestrial lidar, and multispectral and hyperspectral aerial 
imagery collection as well as airborne topographic and bathymetric lidar acquisition, project-level 
GIS implementation, development of file-based geodatabases, and GIS tool development.

Three representative survey and mapping datasets include the National Channel Framework (NCF)—an enterprise 
geodatabase of information on all 61 USACE-maintained high-tonnage channels —hydrographic surveys, which 
provide assistance in locating navigable channels, determining dredging requirements, verifying dredging 
accuracy, and maintaining harbors and rivers —and Inland Electronic Navigational Charts(IENC), accurate 
navigational charts provided in a highly structured data format for use in navigation systems and to increase 
overall navigational safety.. 

< ehydro:where=None >""",
        },
        'ngs': {
            'class': ngs.NGS,
            'description': """NGS Monuments
NGS provides Information about survey marks (including bench marks) in text datasheets or in GIS shapefiles. 
Note some survey markers installed by other organizations may not be available through NGS.

< ngs:datum=geoidHt >""",
        },
        'nos': {
            'class': nos.HydroNOS,
            'description': """NOS Soundings (bag/hydro)
NCEI maintains the National Ocean Service Hydrographic Data Base (NOSHDB) and Hydrographic 
Survey Meta Data Base (HSMDB). Both are populated by the Office of Coast Survey and National 
Geodetic Service, and provide coverage of coastal waters and the U.S. exclusive economic zone 
and its territories. 

Layer 0: Surveys with BAGs available (Bathymetric Attributed Grids).
Layer 1: Surveys with digital sounding data available for download (including those with BAGs).

< nos:where=None:layer=0 >""",
        },
        'charts': {
            'class': charts.NauticalCharts,
            'description': """NOAA Nautical CHARTS

< charts >""",
        },
        'digital_coast': {
            'class': digital_coast.DAV,
            'description': """NOAA DIGITAL COAST elevation data

< digital_coast:where=None:inc=None >""",
        },
        'ncei_thredds': {
            'class': ncei_thredds.NCEIThreddsCatalog,
            'description': """NOAA NCEI DEMs via THREDDS
Digital Elevation Models around the world at various resolutions and extents.
NCEI builds and distributes high-resolution, coastal digital elevation models (DEMs) that integrate ocean 
bathymetry and land topography supporting NOAA's mission to understand and predict changes in Earth's environment, 
and conserve and manage coastal and marine resources to meet our Nation's economic, social, and environmental needs.

DEMs are used for coastal process modeling (tsunami inundation, storm surge, sea-level rise, contaminant dispersal, 
etc.), ecosystems management and habitat research, coastal and marine spatial planning, and hazard mitigation and 
community preparedness.

< ncei_thredds:where=None >""",
        },
        'tnm': {
            'class': tnm.TheNationalMap,
            'description': """USGS' The National Map
Various datasets from USGS's National Map. The National Map is a 
collaborative effort among the USGS and other Federal, State, and local partners to improve
and deliver topographic information for the Nation.

< tnm:formats=None:extents=None:q=None >""",
        },
        'emodnet': {
            'class': emodnet.EMODNet,
            'description': """EU elevation data extracts from EMOD DTM.

< emodnet >""",
        },
        'chs': {
            'class': chs.CHS,
            'description': """Canadian Hydrographic Surveys
CHS NONNA 10m and 100m Bathymetric Survey Grids; Non-Navigational gridded bathymetric data based on charts and soundings.

< chs >""",
        },
        'hrdem': {
            'class': hrdem.HRDEM,
            'description': """High-Resolution Digital Elevation Model data for Canada

< hrdem >""",
        },
        'osm': {
            'class': osm.OpenStreetMap,
            'description': """OpenStreetMap data. 
OpenStreetMap is a free, editable map of the whole world that is 
being built by volunteers largely from scratch and released with an 
open-content license.

< osm:q=None:fmt=osm >""",
        },
        'copernicus': {
            'class': copernicus.CopernicusDEM,
            'description': """COPERNICUS sattelite elevation data
The Copernicus DEM is a Digital Surface Model (DSM) which represents the surface of the Earth including buildings, 
infrastructure and vegetation.

< copernicus >""",
        },
        'nasadem': {
            'class': nasadem.NASADEM,
            'description': """NASA Digital Elevation Model

< nasadem >""",
        },
        'tides': {
            'class': tides.Tides,
            'description': """TIDE station information from NOAA/NOS

< tides:station_id=None:s_datum=mllw:t_datum=msl:units=m >""",
        },
        'vdatum': {
            'class': vdatum.VDATUM,
            'description': """NOAA's VDATUM transformation grids
VDatum is a free software tool being developed jointly by NOAA's National Geodetic Survey (NGS), 
Office of Coast Survey (OCS), and Center for Operational Oceanographic Products and Services (CO-OPS). 

VDatum is designed to vertically transform geospatial data among a variety of tidal, orthometric and 
ellipsoidal vertical datums - allowing users to convert their data from different horizontal/vertical 
references into a common system and enabling the fusion of diverse geospatial data in desired reference 
levels.

< vdatum:datatype=None:gtx=False >""",
        },
        'buoys': {
            'class': buoys.BUOYS,
            'description': """NOAA BUOY data (beta)
A sustainable and resilient marine observation and monitoring infrastructure which enhances healthy 
ecosystems, communities, and economies in the face of change and To provide quality observations in 
the marine environment in a safe and sustainable manner to support the understanding of and predictions 
to changes in weather, climate, oceans and coast. 

< buoys >""",
        },
        'earthdata': {
            'class': earthdata.EarthData,
            'description': """ACCESS NASA EARTH SCIENCE DATA
NASA promotes the full and open sharing of all its data to research and applications communities, 
private industry, academia, and the general public. In order to meet the needs of these different 
communities, NASA’s Earth Observing System Data and Information System (EOSDIS) has provided various 
ways to discover, access, and use the data.

If version is omitted, will fetch all versions
Use wildcards in 'short_name' to return granules for all matching short_name entries.

< earthdata:short_name=ATL08:version=004:time_start='':time_end='':filename_filter='' >""",
        },
        'usiei': {
            'class': usiei.USIEI,
            'description': """US Interagency Elevation Inventory

No data is fetched with this module. Will list out query results from the USIEI.

layers:
  0 - Lidar-Topobathy
  1 - Lidar-Bathy
  2 - Lidar-Topo
  3 - IfSAR/InSAR
  4 - Other Bathy

< usiei:where=None:layer=0 >""",
        },
    }
    
    def __init__(
            self,
            mod=None,
            src_region=None,
            callback=lambda: False,
            weight=None,
            verbose=True
    ):
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
                self.mod_args = utils.args2dict(list(opts[1:]), {})
            self.mod = opts[0]
        else:
            utils.echo_error_msg('could not parse module `{}`'.format(opts[0]))
            return(None)
        
        return(self)

    def add_module(self, type_def={}):
        for key in type_def.keys():
            self.mods[key] = type_def[key]

    def acquire(self, **kwargs):
        """Acquire the fetches module self.mod"""
        
        return(
            self.mods[self.mod]['class'](
                src_region=self.region,
                callback=self.callback,
                weight=self.weight,
                verbose=self.verbose,
                **kwargs,
                **self.mod_args
            )
        )
        
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

usage: {cmd} [ -hlmpqR [ args ] ] MODULE ...

Options:
  -R, --region\t\tRestrict processing to the desired REGION 
\t\t\tWhere a REGION is xmin/xmax/ymin/ymax[/zmin/zmax[/wmin/wmax]]
\t\t\tUse '-' to indicate no bounding range; e.g. -R -/-/-/-/-10/10/1/-
\t\t\tOR an OGR-compatible vector file with regional polygons. 
\t\t\tWhere the REGION is /path/to/vector[:zmin/zmax[/wmin/wmax]].
\t\t\tIf a vector file is supplied, will use each region found therein.
  -l, --list\t\tReturn a list of fetch URLs in the given region.
  -p, --process\t\tProcess fetched elevation data to ASCII XYZ format in WGS84. <beta>
  -q, --quiet\t\tLower the verbosity to a quiet

  --modules\t\tDisply the module descriptions and usage
  --help\t\tPrint the usage text
  --version\t\tPrint the version information

Supported FETCHES modules (see fetches --modules <module-name> for more info): 
  {f_formats}

CIRES DEM home page: <http://ciresgroups.colorado.edu/coastalDEM>
""".format(cmd=os.path.basename(sys.argv[0]), 
           f_version=fetches.__version__,
           f_formats=_fetches_module_short_desc())

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
            print('{}, version {}'.format(
                os.path.basename(sys.argv[0]), fetches.__version__)
                  )
            sys.exit(1)
        elif arg == '--modules' or arg == '-m':
            try:
                if argv[i + 1] in FetchesFactory.mods.keys():
                    sys.stderr.write(
                        _fetches_module_long_desc(
                            {k: FetchesFactory.mods[k] for k in (argv[i + 1],)}
                        )
                    )
                else:
                    sys.stderr.write(_fetches_module_long_desc(FetchesFactory.mods))
                    
            except:
                sys.stderr.write(_fetches_module_long_desc(FetchesFactory.mods))
                
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
        if tmp_region.valid_p(check_xy=True):
            these_regions.append(tmp_region)
        else:
            i_region_s = i_region.split(':')
            tmp_region = regions.ogr_wkts(i_region_s[0])
            for i in tmp_region:
                if i.valid_p():
                    if len(i_region_s) > 1:
                        these_regions.append(
                            regions.Region().from_string(
                                '/'.join([i.format('str'), i_region_s[1]])
                            )
                        )
                    else:
                        these_regions.append(i)

    if not these_regions:
        these_regions = [regions.Region().from_string('-R-180/180/-90/90')]
    if want_verbose:
        utils.echo_msg('parsed {} region(s)'.format(len(these_regions)))

    for rn, this_region in enumerate(these_regions):
        if stop_threads:
            return
        
        x_fs = [FetchesFactory(
            mod=mod,
            src_region=this_region,
            verbose=want_verbose
        ).acquire(dst_srs='epsg:4326') for mod in mods]
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
                fr = f_utils.fetch_results(x_f, want_proc=want_proc)
                fr.daemon = True
                _p = utils.CliProgress('fetching {} remote data files'.format(len(x_f.results)))
                try:
                    fr.start()
                    while True:
                        time.sleep(2)
                        sys.stderr.write('\x1b[2K\r')
                        perc = float((len(x_f.results)-fr.fetch_q.qsize())) / len(x_f.results)*100 if len(x_f.results) > 0 else 1
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
                        except Empty:
                            continue
                        
                        fr.fetch_q.task_done()
                fr.join()
                _p.end(x_f.status, 'fetched {} remote data files'.format(len(x_f.results)))
            if want_verbose:
                utils.echo_msg('ran fetch module {} on region {}...\
            '.format(x_f.name, this_region.format('str')))
                
### End
