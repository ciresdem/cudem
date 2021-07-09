### dlim.py - DataLists IMproved
##
## Copyright (c) 2010 - 2021 CIRES Coastal DEM Team
##
## dlim.py is part of CUDEM
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
## A datalist is similar to an MBSystem datalist; it is a space-delineated file containing the following columns:
## data-path data-format data-weight data-name data-source data-date data-resolution data-type data-horz data-vert data-url
## Minimally, data-path is all that is needed.
##
### Code:

import os
import sys
import re
import glob
import cudem
from cudem import utils
from cudem import regions
from cudem import datasets
import cudem.fetches.fetches as fetches

## ==============================================
## Datalist Class - Recursive data structure
## ==============================================
class Datalist(datasets.XYZDataset):
    """representing a datalist parser
    
    A datalist is an extended MB-System style datalist.
    """

    def __init__(self, fmt=None, **kwargs):
        super().__init__(**kwargs)
        self.name = os.path.basename('.'.join(self.fn.split('.')[:-1]))
        
    def generate_inf(self, callback=lambda: False):
        """return the region of the datalist and generate
        an associated `.inf` file if `inf_file` is True.

        Args:
          dl (str): a datalist pathname
          inf_file (bool): generate an inf file
          epsg (int): EPSG code
          overwrite (bool): overwrite a possibly existing inf_file

        Returns:
          list: the region [xmin, xmax, ymin, ymax]
        """
        
        _region = self.region
        self.region = None
        out_regions = []
        out_region = None
        self.infos['name'] = self.fn
        self.infos['numpts'] = 0
        self.infos['hash'] = self.hash()#dl_hash(self.fn)
        for entry in self.parse():
            if self.verbose:
                callback()
                
            out_regions.append(entry.infos['minmax'])
            if 'numpts' in self.infos.keys():
                self.infos['numpts'] += entry.infos['numpts']

        l = 0

        for this_region in out_regions:
            if l == 0:
                tmp_region = regions.Region().from_list(this_region)
                if tmp_region.valid_p():
                    out_region = regions.Region().from_list(this_region)
                    l += 1
            else:
                tmp_region = regions.Region().from_list(this_region)
                if tmp_region.valid_p():
                    out_region = regions.regions_merge(out_region, tmp_region)
                    
        if out_region is not None:
            self.infos['minmax'] = out_region.export_as_list(include_z=True)
            self.infos['wkt'] = out_region.export_as_wkt()
        else:
            self.infos['minmax'] = None

        self.region = _region
        return(self.infos)
    
    def parse(self):
        """import a datalist entry from a string
    
        Returns:
          datalist_parser: self
        """
        
        if self.verbose:
            _prog = utils.CliProgress('parsing datalist {}{}'.format(self.fn, ' @{}'.format(self.weight) if self.weight is not None else ''))
        if os.path.exists(self.fn):
            with open(self.fn, 'r') as f:
                count = sum(1 for _ in f)
                
            with open(self.fn, 'r') as op:
                for l, this_line in enumerate(op):
                    if self.verbose:
                        _prog.update_perc((l, count))
                    if this_line[0] != '#' and this_line[0] != '\n' and this_line[0].rstrip() != '':
                        data_set = DatasetFactory(
                            this_line,
                            parent=self,
                            name=self.name,
                            src_region=self.region,
                            source=self.source,
                            date=self.date,
                            data_type=self.data_type,
                            resolution=self.resolution,
                            hdatum=self.hdatum,
                            vdatum=self.vdatum,
                            url=self.url,
                            title=self.title,
                            warp=self.warp,
                            weight=self.weight,
                            verbose=self.verbose
                        ).acquire_dataset()
                        if data_set is not None and data_set.valid_p(
                                fmts=DatasetFactory.data_types[data_set.data_format]['fmts']
                        ):
                            if self.region is not None and self.region.valid_p(check_xy=True):
                                if data_set.infos['minmax'] is not None:
                                    # inf_region = regions.Region().from_string(
                                    #     data_set.infos['wkt']
                                    # )
                                    inf_region = regions.Region().from_list(
                                        data_set.infos['minmax']
                                    )
                                    inf_region.wmin = data_set.weight
                                    inf_region.wmax = data_set.weight
                                    #print(data_set.weight)
                                    #inf.region.z_min = data_set.infos['zmin'
                                    #print(inf_region.format('fstr'))
                                    #print(self.region.format('fstr'))
                                    if regions.regions_intersect_p(inf_region, self.region):
                                        for ds in data_set.parse():
                                            self.data_entries.append(ds)
                                            #self.parse_data_lists()
                                            yield(ds)
                            else:
                                for ds in data_set.parse():
                                    self.data_entries.append(ds)
                                    #self.parse_data_lists()
                                    yield(ds)
            #self.parse_data_lists()
        else:
            utils.echo_warning_msg(
                'could not open datalist/entry {}'.format(self.fn)
            )
        if self.verbose:
            _prog.end(0, 'parsed datalist {}{}'.format(self.fn, ' @{}'.format(self.weight) if self.weight is not None else ''))
           
    def yield_xyz(self):
        """parse the data from the datalist

        Yields:
          xyz: the parsed xyz data
        """

        for this_entry in self.parse():
            for xyz in this_entry.yield_xyz():
                yield(xyz)
                
            if this_entry.remote:
                utils.remove_glob('{}*'.format(this_entry.fn))

    def yield_xyz_from_entries(self):
        """parse the data from the datalist

        Yields:
          xyz: the parsed xyz data
        """

        for this_entry in self.data_entries:
            #for this_entry in self.parse():
            for xyz in this_entry.yield_xyz():
                yield(xyz)
                
            if this_entry.remote:
                utils.remove_glob('{}*'.format(this_entry.fn))
                
## ==============================================
## Dataset generator
## ==============================================
class DatasetFactory:

    data_types = {
        -1: {'name': 'datalist',
             'fmts': ['datalist', 'mb-1'],
             'class': lambda k: Datalist(**k),
             },
        167: {'name': 'yxz',
              'fmts': ['yxz'],
              'class': lambda k: datasets.XYZFile(**k),
              },
        168: {'name': 'xyz',
              'fmts': ['xyz', 'csv', 'dat', 'ascii', 'txt'],
              'class': lambda k: datasets.XYZFile(**k),
              },
        200: {'name': 'raster',
              'fmts': ['tif', 'img', 'grd', 'nc', 'vrt', 'bag'],
              'class': lambda k: datasets.RasterFile(**k),
              },
        300: {'name': 'las',
              'fmts': ['las', 'laz'],
              'class': lambda k: datasets.LASFile(**k),
              },
        -11: {'name': 'fetches',
              'fmts': [
                  'gmrt',
                  'multibeam',
                  'usace',
                  'mar_grav',
                  'srtm_plus',
                  'ngs',
                  'nos',
                  'charts',
                  'digital_coast',
                  'ncei_thredds',
                  'tnm',
                  'emodnet',
                  'chs',
                  'hrdem',
                  'copernicus',
                  'nasadem'
              ],
              'class': lambda k: fetches.Fetcher(remote=True, **k),
              },
    }

    datalist_cols = ['path', 'format', 'weight', 'name', 'source',
                     'date', 'type', 'resolution', 'horz', 'vert',
                     'url']
    
    def __init__(
            self,
            fn=None,
            data_format=None,
            weight=1,
            epsg=4326,
            name="xyz_dataset",
            title=None,
            source=None,
            date=None,
            data_type=None,
            resolution=None,
            hdatum=None,
            vdatum=None,
            url=None,
            parent=None,
            src_region=None,
            warp=None,
            verbose=False
    ):
        
        self.name = name
        self.title = title
        self.data_format = data_format
        self.weight = weight
        self.source = source
        self.date = date
        self.data_type = data_type
        self.resolution = resolution
        self.epsg = epsg
        self.hdatum = hdatum
        self.vdatum = vdatum
        self.url = url
        self.warp = warp
        self.region = src_region
        self.parent = parent
        self.verbose = verbose
        self.fn = fn
        self.parse_fn()

        #if self.weight is not None:
        #    self.weight = self.weight**4
        
        if self.data_format is None:
            self.guess_data_format()
        
    def parse_fn(self):

        if self.fn is None: return(self)
        if os.path.exists(self.fn): return(self.fn)

        this_entry = re.findall(r'[^"\s]\S*|".+?"', self.fn.rstrip())

        try:
            entry = [x if n == 0 else utils.int_or(x) if n < 2 else utils.float_or(x) if n < 3 else x \
                     for n, x in enumerate(this_entry)]
        except Exception as e:
            utils.echo_error_msg('could not parse entry {}'.format(self.fn))
            return(self)

        if len(entry) < 2:
            for key in self.data_types.keys():
                se = entry[0].split('.')
                see = se[-1] if len(se) > 1 else entry[0].split(":")[0]
                if see in self.data_types[key]['fmts']:
                    entry.append(int(key))
                    break
                
            if len(entry) < 2:
                utils.echo_error_msg('could not parse entry {}'.format(self.fn))
                return(self)

        if len(entry) < 3:
            entry.append(1)
        elif entry[2] is None:
            #entry[2] = self.weight
            entry[2] = 1
            #else:
        #entry[2] *= self.weight
        if self.parent is not None:
            #entry[2] *= self.parent.weight
            #self.weight *= self.parent.weight
            if self.weight is not None:
                self.weight *= entry[2]
            else:
                self.weight = entry[2]

        #print(self.weight)
        
        if len(entry) < 4:
            entry.append(self.title)
        else:
            self.title = entry[3]
            
        if len(entry) < 5:
            entry.append(self.source)
        else:
            self.source = entry[4]
            
        if len(entry) < 6:
            entry.append(self.date)
        else:
            self.date = entry[5]
            
        if len(entry) < 7:
            entry.append(self.data_type)
        else:
            self.data_type = entry[6]
            
        if len(entry) < 8:
            entry.append(self.resolution)
        else:
            self.resolution = entry[7]
            
        if len(entry) < 9:
            entry.append(self.hdatum)
        else:
            self.hdatum = entry[8]
            
        if len(entry) < 10:
            entry.append(self.vdatum)
        else:
            self.vdatum = entry[9]
            
        if len(entry) < 11:
            entry.append(self.url)
        else:
            self.url = entry[10]

        if self.parent is None or entry[1] == -11:
            self.fn = entry[0]
        else:
            self.fn = os.path.join(
                os.path.dirname(self.parent.fn), entry[0]
            )
        
        self.data_format = entry[1]
        if self.data_format is None:
            self.guess_data_format()

        return(self)

    def guess_data_format(self):
        if self.fn is not None:
            for key in self.data_types.keys():
                if self.fn.split('.')[-1] in self.data_types[key]['fmts']:
                    self.data_format = key
                    break

    def add_data_type(self, type_def = {}):
        for key in type_def.keys():
            self.data_types[key] = type_def[key]

    def acquire_datalist(self, **kwargs):
        return(
            Datalist(
                fn=self.fn,
                data_format=self.data_format,
                weight=self.weight,
                src_region=self.region,
                title=self.title,
                source=self.source,
                date=self.date,
                data_type=self.data_type,
                resolution=self.resolution,
                hdatum=self.hdatum,
                vdatum=self.vdatum,
                url=self.url,
                epsg=self.epsg,
                warp=self.warp,
                name=self.name,
                parent=self.parent,
                verbose=self.verbose,
                **kwargs
            )
        )

    def acquire_xyz_file(self, **kwargs):
        return(
            datasets.XYZFile(
                fn=self.fn,
                data_format=self.data_format,
                weight=self.weight,
                src_region=self.region,
                title=self.title,
                source=self.source,
                date=self.date,
                data_type=self.data_type,
                resolution=self.resolution,
                hdatum=self.hdatum,
                vdatum=self.vdatum,
                url=self.url,
                epsg=self.epsg,
                warp=self.warp,
                name=self.name,
                parent=self.parent,
                verbose=self.verbose,
                xpos=0,
                ypos=1,
                zpos=2,
                **kwargs
            )
        )

    def acquire_las_file(self, **kwargs):
        return(
            datasets.LASFile(
                fn=self.fn,
                data_format=self.data_format,
                weight=self.weight,
                src_region=self.region,
                title=self.title,
                source=self.source,
                date=self.date,
                data_type=self.data_type,
                resolution=self.resolution,
                hdatum=self.hdatum,
                vdatum=self.vdatum,
                url=self.url,
                epsg=self.epsg,
                warp=self.warp,
                name=self.name,
                parent=self.parent,
                verbose=self.verbose,
                classes=[2,29],
                **kwargs
            )
        )    
    
    def acquire_raster_file(self, **kwargs):
        return(
            datasets.RasterFile(
                fn=self.fn,
                data_format=self.data_format,
                weight=self.weight,
                src_region=self.region,
                title=self.title,
                source=self.source,
                date=self.date,
                data_type=self.data_type,
                resolution=self.resolution,
                hdatum=self.hdatum,
                vdatum=self.vdatum,
                url=self.url,
                epsg=self.epsg,
                warp=self.warp,
                name=self.name,
                parent=self.parent,
                verbose=self.verbose,
                **kwargs
            )
        )

    def acquire_fetcher(self, **kwargs):
        return(
            fetches.Fetcher(
                fn=self.fn,
                data_format=self.data_format,
                weight=self.weight,
                src_region=self.region,
                title=self.title,
                source=self.source,
                date=self.date,
                data_type=self.data_type,
                resolution=self.resolution,
                hdatum=self.hdatum,
                vdatum=self.vdatum,
                url=self.url,
                epsg=self.epsg,
                warp=self.warp,
                name=self.name,
                parent=self.parent,
                verbose=self.verbose,
                remote=True,
                **kwargs
            )
        )

    def acquire(self, **kwargs):
        return(
            self.data_types[self.data_format]['class'](
                fn=self.fn,
                data_format=self.data_format,
                weight=self.weight,
                src_region=self.region,
                title=self.title,
                source=self.source,
                date=self.date,
                data_type=self.data_type,
                resolution=self.resolution,
                hdatum=self.hdatum,
                vdatum=self.vdatum,
                url=self.url,
                epsg=self.epsg,
                warp=self.warp,
                name=self.name,
                parent=self.parent,
                verbose=self.verbose,
                **kwargs
            )
        )
    
    def acquire_dataset(self, **kwargs):
        if self.data_format == -1:
            return(self.acquire_datalist(**kwargs))#.parse())
        
        if self.data_format == 168:
            return(self.acquire_xyz_file(**kwargs))#.parse())
        
        if self.data_format == 200:
            return(self.acquire_raster_file(**kwargs))#.parse())
        
        if self.data_format == 300:
            return(self.acquire_las_file(**kwargs))#.parse())
        
        if self.data_format == -11:
            return(self.acquire_fetcher(**kwargs))#.parse())

_datalist_fmts_long_desc = lambda: '\n  '.join(
    ['{}:\t{}'.format(key, DatasetFactory().data_types[key]['name']) for key in DatasetFactory().data_types])
_datalist_fmts_short_desc = lambda: ',  '.join(
    ['{} ({})'.format(DatasetFactory().data_types[key]['name'], key) for key in DatasetFactory().data_types])

## ==============================================
## Command-line Interface (CLI)
## $ dlim
##
## datalists cli
## ==============================================
datalists_usage = """{cmd} ({dl_version}): DataLists IMproved; Process and generate datalists

usage: {cmd} [ -ghiqwPRW [ args ] ] DATALIST,FORMAT,WEIGHT ...

Options:
  -R, --region\t\tRestrict processing to the desired REGION 
\t\t\tWhere a REGION is [ xmin/xmax/ymin/ymax/[ zmin/zmax/[ wmin/wmax ] ] ]
\t\t\tUse '-' to indicate no bounding range; e.g. -R -/-/-/-/-10/10/1/-
\t\t\tOR an OGR-compatible vector file with regional polygons. 
\t\t\tIf a vector file is supplied, will use each region found therein.
\t\t\tAppend :zmin/zmax/[ wmin/wmax ] to the file path to extended REGION.
  -P, --s_epsg\t\tSet the projection EPSG code of the datalist.
  -W, --t_epsg\t\tSet the output warp projection EPSG code.
  -F, --format\t\tOnly process the given data format.

  --archive\t\tARCHIVE the datalist to the given REGION
  --glob\t\tGLOB the datasets in the current directory to stdout
  --info\t\tGenerate and return an INFO dictionary of the dataset
  --weights\t\tOutput WEIGHT values along with xyz
  --quiet\t\tLower the verbosity to a quiet

  --help\t\tPrint the usage text
  --version\t\tPrint the version information

Supported datalist formats: 
  {dl_formats}

Examples:
  % {cmd} my_data.datalist -R -90/-89/30/31
  % {cmd} -R-90/-89/30/31/-100/100 *.tif -l -w > tifs_in_region.datalist
  % {cmd} -R my_region.shp my_data.xyz -w -s_epsg 4326 -t_epsg 3565 > my_data_3565.xyz

CIRES DEM home page: <http://ciresgroups.colorado.edu/coastalDEM>\
""".format(cmd=os.path.basename(sys.argv[0]), 
           dl_version=cudem.__version__,
           dl_formats=_datalist_fmts_short_desc())

def datalists_cli(argv = sys.argv):
    """run datalists from command-line

    See `datalists_cli_usage` for full cli options.
    """

    dls = []
    epsg = None
    warp = None
    fmt = None
    i_regions = []
    these_regions = []
    want_weights = False
    want_inf = False
    want_list = False
    want_glob = False
    want_archive = False
    want_verbose = True
    want_region = False
    
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
        elif arg == '-s_epsg' or arg == '--s_epsg' or arg == '-P':
            epsg = argv[i + 1]
            i = i + 1
        elif arg == '-t_epsg' or arg == '--t_epsg' or arg == '-W':
            warp = argv[i + 1]
            i = i + 1
        elif arg == '--format' or arg == '-F':
            fmt = utils.int_or(argv[i + 1])
            i = i + 1
        elif arg == '--archive' or arg == '-a':
            want_archive = True
        elif arg == '--weights' or arg == '-w':
            want_weights = True
        elif arg == '--info' or arg == '-i':
            want_inf = True
        elif arg == '--region_inf' or arg == '-r':
            want_region = True
        elif arg == '--list' or arg == '-l':
            want_list = True
        elif arg == '--glob' or arg == '-g':
            want_glob = True
        elif arg == '--quiet' or arg == '-q':
            want_verbose = False
        elif arg == '--help' or arg == '-h':
            print(datalists_usage)
            sys.exit(1)
        elif arg == '--version' or arg == '-v':
            print('{}, version {}'.format(
                os.path.basename(sys.argv[0]), cudem.__version__)
                  )
            sys.exit(1)
        elif arg[0] == '-':
            print(datalists_usage)
            sys.exit(0)
        else: dls.append(arg)
        
        i = i + 1

    if want_glob:
        for key in DatasetFactory().data_types.keys():
            if key != -1:
                for f in DatasetFactory().data_types[key]['fmts']:
                    globs = glob.glob('*.{}'.format(f))
                    [sys.stdout.write(
                        '{}\n'.format(
                            ' '.join(
                                [x, str(key), '1']
                            )
                        )
                    ) for x in globs]
                    
        sys.exit(0)

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

    if len(these_regions) == 0:
        these_regions = [None]
    else:
        if want_verbose:
            utils.echo_msg('parsed {} region(s)'.format(len(these_regions)))

    for rn, this_region in enumerate(these_regions):
        if len(dls) == 0:
            print(datalists_usage)
            utils.echo_error_msg('you must specify some type of data')
        xdls = [DatasetFactory(
            fn=" ".join(['-' if x == "" else x for x in dl.split(",")]),
            src_region=this_region, verbose=want_verbose,
            epsg=epsg, warp=warp).acquire_dataset() for dl in dls]

        for xdl in xdls:
            if xdl is not None and xdl.valid_p(
                    fmts=DatasetFactory.data_types[xdl.data_format]['fmts']
            ):
                #xdl.parse()
                if not want_weights:
                    xdl.weight = None
                    
                if want_inf:
                    print(xdl.inf())
                elif want_list:
                    xdl.echo()
                elif want_archive:
                    [x for x in xdl.archive_xyz()]
                elif want_region:
                    print(regions.Region().from_list(xdl.inf()['minmax']).format('gmt'))
                else: xdl.dump_xyz()
### End
