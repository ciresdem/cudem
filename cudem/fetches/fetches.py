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
            'description': """The Global Multi-Resolution Topography synthesis.
The Global Multi-Resolution Topography (GMRT) synthesis is a multi-resolutional 
compilation of edited multibeam sonar data collected by scientists and institutions worldwide, that is 
reviewed, processed and gridded by the GMRT Team and merged into a single continuously updated compilation 
of global elevation data. The synthesis began in 1992 as the Ridge Multibeam Synthesis (RMBS), was expanded 
to include multibeam bathymetry data from the Southern Ocean, and now includes bathymetry from throughout 
the global and coastal oceans.

layers: 'topo' or 'topo-mask'

https://www.gmrt.org

< gmrt:res=max:fmt=geotiff:bathy_only=False:layer=topo >""",
        },
        'gebco': {
            'class': gebco.GEBCO,
            'description': """GEBCO
GEBCO’s current gridded bathymetric data set, the GEBCO_2022 Grid, is a global terrain model for ocean and land, 
providing elevation data, in meters, on a 15 arc-second interval grid. It is accompanied by a Type Identifier 
(TID) Grid that gives information on the types of source data that the GEBCO_2022 Grid is based. 

https://www.gebco.net

< gebco >""",
        },
        'srtm_plus': {
            'class': srtm.SRTMPlus,
            'description': """SRTM15+: GLOBAL BATHYMETRY AND TOPOGRAPHY AT 15 ARCSECONDS.

https://topex.ucsd.edu/WWW_html/srtm15_plus.html

< srtm_plus >""",
        },
        'mar_grav': {
            'class': mar_grav.MarGrav,
            'description': """MARine GRAVity Satellite Altimetry Topography from Scripps.

https://topex.ucsd.edu/WWW_html/mar_grav.html

< mar_grav >""",
        },
        'multibeam': {
            'class': mb.Multibeam,
            'description': """NOAA MULTIBEAM bathymetric data.
NCEI is the U.S. national archive for multibeam bathymetric data and holds more than 9 million 
nautical miles of ship trackline data recorded from over 2400 cruises and received from sources 
worldwide.

https://data.ngdc.noaa.gov/platforms/

< multibeam:process=False:min_year=None:survey_id=None:exclude=None >""",
        },
        'mgds': {
            'class': mgds.MGDS,
            'description': """The Marine Geoscience Data System (MGDS) is a trusted data 
repository that provides free public access to a curated collection of marine geophysical 
data products and complementary data related to understanding the formation and evolution 
of the seafloor and sub-seafloor.

https://www.marine-geo.org

data_tpye=[Bathymetry, Bathymetry:Phase, Bathymetry:Swath, Bathymetry:Swath:Ancillary, Bathymetry:Singlebeam, Bathymetry:BPI, Bathymetry:ReferenceSurface, Bathymetry:Paelobathymetry]
            
< mgds:data_type=Bathymetry >""",
        },
        'trackline': {
            'class': trackline.Trackline,
            'description': """NOAA TRACKLINE bathymetric data.

http://www.ngdc.noaa.gov/trackline/

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

https://navigation.usace.army.mil/Survey/Hydro

< ehydro:where=None >""",
        },
        'ngs': {
            'class': ngs.NGS,
            'description': """NGS Monuments
NGS provides Information about survey marks (including bench marks) in text datasheets or in GIS shapefiles. 
Note some survey markers installed by other organizations may not be available through NGS.

http://geodesy.noaa.gov/

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

https://www.ngdc.noaa.gov/mgg/bathymetry/hydro.html

< nos:where=None:layer=0 >""",
        },
        'charts': {
            'class': charts.NauticalCharts,
            'description': """NOAA Nautical CHARTS

https://www.charts.noaa.gov/

< charts >""",
        },
        'digital_coast': {
            'class': digital_coast.DAV,
            'description': """NOAA DIGITAL COAST elevation data

https://coast.noaa.gov

< digital_coast:where=None >""",
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

https://www.ngdc.noaa.gov/thredds/demCatalog.xml

< ncei_thredds:where=None >""",
        },
        'tnm': {
            'class': tnm.TheNationalMap,
            'description': """USGS' The National Map
Various datasets from USGS's National Map. The National Map is a 
collaborative effort among the USGS and other Federal, State, and local partners to improve
and deliver topographic information for the Nation.

http://tnmaccess.nationalmap.gov/

< tnm:formats=None:extents=None:q=None >""",
        },
        'emodnet': {
            'class': emodnet.EMODNet,
            'description': """EU elevation data extracts from EMOD DTM.

https://portal.emodnet-bathymetry.eu/

< emodnet >""",
        },
        'chs': {
            'class': chs.CHS,
            'description': """Canadian Hydrographic Surveys
CHS NONNA 10m and 100m Bathymetric Survey Grids; Non-Navigational gridded bathymetric data based on charts and soundings.

https://open.canada.ca

< chs >""",
        },
        'hrdem': {
            'class': hrdem.HRDEM,
            'description': """High-Resolution Digital Elevation Model data for Canada

https://open.canada.ca

< hrdem >""",
        },
        'arcticdem': {
            'class': arcticdem.ArcticDEM,
            'description': """Arctic DEM
ArcticDEM is an NGA-NSF public-private initiative to automatically produce a high-resolution, 
high quality, digital surface model (DSM) of the Arctic using optical stereo imagery, 
high-performance computing, and open source photogrammetry software.

https://www.pgc.umn.edu/data/arcticdem/

< arcticdem >""",
        },

        'bluetopo': {
            'class': bluetopo.BlueTopo,
            'description': """BlueTOPO DEM
BlueTopo is a compilation of the nation's best available bathymetric data. 
In the same way that topographic map details the height of land, BlueTopo details the depth of 
lake beds and seafloor beneath navigationally significant U.S. waters. Created as part of the 
Office of Coast Survey nautical charting mission and its National Bathymetric Source project, 
BlueTopo is curated bathymetric source data to provide a definitive nationwide model of the seafloor 
and the Great Lakes.

https://www.nauticalcharts.noaa.gov/data/bluetopo.html

< bluetopo:want_interpolation=False:unc_weights=False:keep_index=False >""",
        },
'osm': {
            'class': osm.OpenStreetMap,
            'description': """OpenStreetMap data. 
OpenStreetMap is a free, editable map of the whole world that is 
being built by volunteers largely from scratch and released with an 
open-content license.

https://wiki.openstreetmap.org/

< osm:q=None:fmt=osm >""",
        },
        'copernicus': {
            'class': copernicus.CopernicusDEM,
            'description': """COPERNICUS sattelite elevation data
The Copernicus DEM is a Digital Surface Model (DSM) which represents the surface of the Earth including buildings, 
infrastructure and vegetation.

https://ec.europa.eu/eurostat/web/gisco/geodata/reference-data/elevation/copernicus-dem/elevation

< copernicus >""",
        },
        'fabdem': {
            'class': fabdem.FABDEM,
            'description': """FABDEM elevation data
FABDEM (Forest And Buildings removed Copernicus DEM) is a global elevation map that removes building and tree height
biases from the Copernicus GLO 30 Digital Elevation Model (DEM). The data is available at 1 arc second
grid spacing (approximately 30m at the equator) for the globe.

https://data.bris.ac.uk/data/dataset/s5hqmjcdj8yo2ibzi9b4ew3sn

< fabdem >""",
        },
        'nasadem': {
            'class': nasadem.NASADEM,
            'description': """NASA Digital Elevation Model
Our objective is to provide the scientific and civil communities with a state-of-the-art global 
digital elevation model (DEM) derived from a combination of Shuttle Radar Topography Mission (SRTM) 
processing improvements, elevation control, void-filling and merging with data unavailable at the 
time of the original SRTM production.

https://www.earthdata.nasa.gov/esds/competitive-programs/measures/nasadem

< nasadem >""",
        },
        'tides': {
            'class': tides.Tides,
            'description': """TIDE station information from NOAA/NOS

https://tidesandcurrents.noaa.gov/

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

https://vdatum.noaa.gov
https://cdn.proj.org

< vdatum:datatype=None:gtx=False >""",
        },
        'buoys': {
            'class': buoys.BUOYS,
            'description': """NOAA BUOY data (beta)
A sustainable and resilient marine observation and monitoring infrastructure which enhances healthy 
ecosystems, communities, and economies in the face of change and To provide quality observations in 
the marine environment in a safe and sustainable manner to support the understanding of and predictions 
to changes in weather, climate, oceans and coast. 

https://www.ndbc.noaa.gov

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

https://cmr.earthdata.nasa.gov

< earthdata:short_name=ATL08:version=004:time_start='':time_end='':filename_filter='' >""",
        },
        'usiei': {
            'class': usiei.USIEI,
            'description': """US Interagency Elevation Inventory

No data is fetched with this module. Will list out query results from the USIEI.
Set 'want_geometry' to True to output a geojson formatted vector.

layers:
  0 - Lidar-Topobathy
  1 - Lidar-Bathy
  2 - Lidar-Topo
  3 - IfSAR/InSAR
  4 - Other Bathy

https://coast.noaa.gov/inventory/

< usiei:where=None:layer=0:want_geometry=False >""",
        },
        'wsf': {
            'class': wsf.WSF,
            'description': """WSF from German Aerospace Service (DLR)

World Settlement Footprint (WSF) 2019

< usiei:where=None:datatype=None >""",
        },
        'hydrolakes': {
            'class': hydrolakes.HydroLakes,
            'description': """HydroLakes vector and derived elevations
HydroLAKES aims to provide the shoreline polygons of all global lakes with a surface area 
of at least 10 ha. HydroLAKES has been developed using a suite of auxiliary data sources of 
lake polygons and gridded lake surface areas. All lakes are co-registered to the global 
river network of the HydroSHEDS database via their lake pour points. The global coverage of 
HydroLAKES encompasses 1.4 million individual lakes or reservoirs representing a total 
surface area of 2.67 million km², a total shoreline length of 7.2 million km, and a total 
storage volume of 181,900 km³.

https://www.hydrosheds.org/products/hydrolakes

< hydrolakes >""",
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
            #if this_mod is None:
            #    self.mod = None
        
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

    if len(xy_inc) < 2:
        xy_inc.append(xy_inc[0])
        
    elif len(xy_inc) == 0:
        xy_inc = [None, None]
        
    these_regions = regions.parse_cli_region(i_regions, want_verbose)
    if not these_regions:
        these_regions = [regions.Region().from_string('-R-180/180/-90/90')]
        
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
