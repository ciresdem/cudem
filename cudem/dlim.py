### dlim.py - DataLists IMproved
##
## Copyright (c) 2010 - 2022 Regents of the University of Colorado
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
## an associated inf and geojson file will be gerenated for each datalist
## only an associated inf file will be genereated for individual datasets
##
## recursive data-structures which point to datasets (see cudem.datasets) (datalist, zip, fetches, etc)
##
### Code:

import os
import sys
import re
import copy
from osgeo import ogr

import cudem
from cudem import utils
from cudem import regions
from cudem import datasets
from cudem.fetches import fetches

def make_datalist(data_list, weight, region, src_srs, dst_srs, x_inc, y_inc, sample_alg, verbose):
    xdl = Datalist(
        fn='<scratch-datalist>',
        data_format=-1,
        weight=None if not weight else 1,
        src_region=region,
        verbose=verbose,
        parent=None,
        src_srs=src_srs,
        dst_srs=dst_srs,
        x_inc=x_inc,
        y_inc=y_inc,
        sample_alg=sample_alg
    )
    xdl.data_entries = [DatasetFactory(
        fn=" ".join(['-' if x == "" else x for x in dl.split(",")]),
        src_region=region,
        verbose=verbose,
        src_srs=src_srs,
        dst_srs=dst_srs,
        x_inc=x_inc,
        y_inc=y_inc,
        sample_alg=sample_alg,
        parent=xdl
    ).acquire() for dl in data_list]
    return(xdl)

def write_datalist(data_list, outname=None):
    if outname is None:
        outname = '{}_{}'.format(self.metadata['name'], utils.this_year())
    
    if os.path.exists('{}.datalist'.format(outname)):
        utils.remove_glob('{}.datalist*'.format(outname))
        
    with open('{}.datalist'.format(outname), 'w') as tmp_dl:
        [tmp_dl.write('{}\n'.format(x.format_entry())) for x in data_list]

    return('{}.datalist'.format(outname))

def init_data(data_list, this_region, src_srs, dst_srs, xy_inc, sample_alg, want_verbose):

    xdls = [DatasetFactory(
        fn=" ".join(['-' if x == "" else x for x in dl.split(",")]),
        src_region=this_region,
        verbose=want_verbose,
        src_srs=src_srs,
        dst_srs=dst_srs,
        x_inc=xy_inc[0],
        y_inc=xy_inc[1],
        sample_alg=sample_alg
    ).acquire() for dl in data_list]

    if len(xdls) > 1:
        dl_fn = write_datalist(xdls, 'dlim')
        this_datalist = DatasetFactory(
            fn=dl_fn,
            src_region=this_region,
            verbose=want_verbose,
            src_srs=src_srs,
            dst_srs=dst_srs,
            x_inc=xy_inc[0],
            y_inc=xy_inc[1],
            sample_alg=sample_alg
        ).acquire()
    else:
        this_datalist = xdls[0]

    return(this_datalist)
    
## ==============================================
## Datalist Class - Recursive data structure
## ==============================================
class Datalist(datasets.ElevationDataset):
    """representing a datalist parser
    
    A datalist is an extended MB-System style datalist.
    """

    def __init__(self, fmt=None, **kwargs):
        self.v_fields = [
            'Path',
            'Format',
            'Weight',
        ]
        self.t_fields = [
            ogr.OFTString,
            ogr.OFTString,
            ogr.OFTString,
        ]
        super().__init__(**kwargs)
        self.metadata['name'] = os.path.basename('.'.join(self.fn.split('.')[:-1]))

    ##TODO:
    ## add more dataset parameters into json fields
    ## including dataset specific options (ds_args)
    def _init_datalist_vector(self):
        self.dst_layer = '{}'.format(self.fn)
        self.dst_vector = self.dst_layer + '.json'

        utils.remove_glob('{}.json'.format(self.dst_layer))
        if self.src_srs is not None:
            utils.gdal_prj_file('{}.prj'.format(self.dst_layer), self.src_srs)
            
        self.ds = ogr.GetDriverByName('GeoJSON').CreateDataSource(self.dst_vector)
        if self.ds is not None: 
            self.layer = self.ds.CreateLayer(
                '{}'.format(self.dst_layer), None, ogr.wkbMultiPolygon
            )
            [self.layer.CreateField(
                ogr.FieldDefn('{}'.format(f), self.t_fields[i])
            ) for i, f in enumerate(self.v_fields)]
            
            [self.layer.SetFeature(feature) for feature in self.layer]
        else:
            self.layer = None

    def _create_entry_feature(self, entry, entry_region):
        entry_path = os.path.abspath(entry.fn) if not entry.remote else entry.fn
        o_v_fields = [
            entry_path,
            entry.data_format,
            entry.weight,
        ]
        dst_defn = self.layer.GetLayerDefn()
        entry_geom = ogr.CreateGeometryFromWkt(entry_region.export_as_wkt())
        out_feat = ogr.Feature(dst_defn)
        out_feat.SetGeometry(entry_geom)
        for i, f in enumerate(self.v_fields):
            out_feat.SetField(f, o_v_fields[i])
        self.layer.CreateFeature(out_feat)
        
    def generate_inf(self, callback=lambda: False):
        """return the region of the datalist and generate
        an associated `.inf` file if `inf_file` is True.
        """
        
        _region = self.region
        out_region = None
        out_regions = []
        out_srs = []
        self.region = None
        self.infos['name'] = self.fn
        self.infos['numpts'] = 0
        self.infos['hash'] = self.hash()#dl_hash(self.fn)
        self._init_datalist_vector()
        
        dst_srs_ = self.dst_srs
        if self.dst_srs is not None:
            self.dst_srs = self.dst_srs.split('+')[0]
            
        for entry in self.parse():
            if self.verbose:
                callback()

            if entry.src_srs is not None:
                out_srs.append(entry.src_srs)

                if self.dst_srs is not None:
                    #self.infos['src_srs'] = self.dst_srs
                    e_region = regions.Region().from_list(entry.infos['minmax'])
                    e_region.src_srs = entry.src_srs
                    e_region.warp(self.dst_srs)
                    entry_region = e_region.export_as_list(include_z=True)
                else:
                    entry_region = entry.infos['minmax']

            else:
                out_srs.append(None)
                entry_region = entry.infos['minmax']

            if regions.Region().from_list(entry_region).valid_p():
                self._create_entry_feature(entry, regions.Region().from_list(entry_region))
                out_regions.append(entry_region)
                if 'numpts' in self.infos.keys():
                    self.infos['numpts'] += entry.infos['numpts']

        self.ds = self.layer = None
        count = 0
        for this_region in out_regions:
            tmp_region = regions.Region().from_list(this_region)
            if tmp_region.valid_p():
                if count == 0:
                    out_region = tmp_region
                    count += 1
                else:
                    out_region = regions.regions_merge(out_region, tmp_region)
                    
        if out_region is not None:
            self.infos['minmax'] = out_region.export_as_list(include_z=True)
            self.infos['wkt'] = out_region.export_as_wkt()
        else:
            self.infos['minmax'] = None
            
        self.region = _region
        self.dst_srs = dst_srs_
            
        if 'src_srs' not in self.infos.keys() or self.infos['src_srs'] is None:
            if self.src_srs is not None:
                self.infos['src_srs'] = self.src_srs
            else:
                if all(x == out_srs[0] for x in out_srs):
                    self.infos['src_srs'] = out_srs[0]
        else:
            self.src_srs = self.infos['src_srs']
        
        return(self.infos)

    def parse_json(self):
        status = 0
        count = 0
        if self.verbose:
            _prog = utils.CliProgress(
                'parsing datalist json {}{}'.format(
                    self.fn,
                    ' @{}'.format(self.weight) if self.weight is not None else '')
            )

        if os.path.exists('{}.json'.format(self.fn)):
            driver = ogr.GetDriverByName('GeoJSON')
            dl_ds = driver.Open('{}.json'.format(self.fn))
            dl_layer = dl_ds.GetLayer()
            ldefn = dl_layer.GetLayerDefn()
            if self.region is not None:
                _boundsGeom = self.region.export_as_geom() if self.dst_trans is None else self.trans_region.export_as_geom()
            else:
                _boundsGeom = None

            dl_layer.SetSpatialFilter(_boundsGeom)
            count = len(dl_layer)
            
            for l,feat in enumerate(dl_layer):
                _prog.update_perc((l, count))
                if self.region is not None:
                    w_region = self.region.w_region()
                    if w_region[0] is not None:
                        if float(feat.GetField('Weight')) < w_region[0]:
                            continue
                    if w_region[1] is not None:
                        if float(feat.GetField('Weight')) > w_region[1]:
                            continue
                        
                data_set = DatasetFactory(
                    '{} {} {}'.format(feat.GetField('Path'),feat.GetField('Format'),feat.GetField('Weight')),
                    weight=self.weight,
                    parent=self,
                    src_region=self.region,
                    metadata=copy.deepcopy(self.metadata),
                    src_srs=self.src_srs,
                    dst_srs=self.dst_srs,
                    x_inc=self.x_inc,
                    y_inc=self.y_inc,
                    sample_alg=self.sample_alg,
                    verbose=self.verbose
                ).acquire()
                if data_set is not None and data_set.valid_p(
                        fmts=DatasetFactory.data_types[data_set.data_format]['fmts']
                ):
                    for ds in data_set.parse():
                        self.data_entries.append(ds)
                        yield(ds)
            dl_ds = dl_layer = None
        else:
            status = -1
            for ds in self.parse():
                yield(ds)
            
        if self.verbose:
            _prog.end(
                status, 'parsed {} datasets from datalist {}{}'.format(
                    count, self.fn, ' @{}'.format(self.weight) if self.weight is not None else ''
                )
            )
                            
    def parse(self):
        """import a datalist entry from a string"""

        status = 0
        if self.verbose:
            _prog = utils.CliProgress(
                'parsing datalist {}{}'.format(
                    self.fn,
                    ' @{}'.format(self.weight) if self.weight is not None else '')
            )
                    
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
                            weight=self.weight,
                            parent=self,
                            src_region=self.region,
                            metadata=copy.deepcopy(self.metadata),
                            src_srs=self.src_srs,
                            dst_srs=self.dst_srs,
                            x_inc=self.x_inc,
                            y_inc=self.y_inc,
                            sample_alg=self.sample_alg,
                            verbose=self.verbose
                        ).acquire()
                        if data_set is not None and data_set.valid_p(
                                fmts=DatasetFactory.data_types[data_set.data_format]['fmts']
                        ):
                            if self.region is not None and self.region.valid_p(check_xy=True):
                                try:
                                    inf_region = regions.Region().from_list(
                                        data_set.infos['minmax']
                                    )
                                except:
                                    inf_region = self.region.copy()

                                inf_region.wmin = data_set.weight
                                inf_region.wmax = data_set.weight
                                if regions.regions_intersect_p(
                                        inf_region,
                                        self.region if data_set.dst_trans is None else data_set.trans_region
                                ):
                                    for ds in data_set.parse():
                                        self.data_entries.append(ds)
                                        yield(ds)
                            else:
                                for ds in data_set.parse():
                                    self.data_entries.append(ds)
                                    yield(ds)
        elif len(self.data_entries) > 0:
            for data_set in self.data_entries:
                for ds in data_set.parse():
                    yield(ds)
        else:
            if self.verbose:
                utils.echo_warning_msg(
                    'could not open datalist/entry {}'.format(self.fn)
                )
            status = -1
            
        if self.verbose:
            _prog.end(status, 'parsed datalist {}{}'.format(
                self.fn, ' @{}'.format(self.weight) if self.weight is not None else ''
            ))
           
    def yield_xyz(self):
        """parse the data from the datalist and yield as xyz"""
        
        for this_entry in self.parse_json():
            #for this_entry in self.parse():
            for xyz in this_entry.yield_xyz():
                yield(xyz)
                
            if this_entry.remote:
                utils.remove_glob('{}*'.format(this_entry.fn))

    def yield_array(self):
        """parse the data from the datalist and yield as array"""
        
        for this_entry in self.parse_json():
            #for this_entry in self.parse():
            for arr in this_entry.yield_array():
                yield(arr)

            if this_entry.remote:
                utils.remove_glob('{}*'.format(this_entry.fn))
                
## ==============================================
## ZIPlist Class - Recursive data structure - testing
## ==============================================
class ZIPlist(datasets.ElevationDataset):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def generate_inf(self, callback=lambda: False):
        """return the region of the datalist and generate
        an associated `.inf` file if `inf_file` is True.
        """
        
        _region = self.region
        out_region = None
        out_regions = []
        self.region = None
        self.infos['name'] = self.fn
        self.infos['numpts'] = 0
        self.infos['hash'] = self.hash()#dl_hash(self.fn)
        self.infos['format'] = self.data_format
        
        for entry in self.parse_():
            if self.verbose:
                callback()

            out_regions.append(entry.infos['minmax'])
            if 'numpts' in self.infos.keys():
                self.infos['numpts'] += entry.infos['numpts']

        count = 0
        for this_region in out_regions:
            tmp_region = regions.Region().from_list(this_region)
            if tmp_region.valid_p():
                if count == 0:
                    out_region = tmp_region
                    count += 1
                else:
                    out_region = regions.regions_merge(out_region, tmp_region)
                    
        if out_region is not None:
            self.infos['minmax'] = out_region.export_as_list(include_z=True)
            self.infos['wkt'] = out_region.export_as_wkt()
        else:
            self.infos['minmax'] = None

        self.region = _region
        self.infos['src_srs'] = self.src_srs
        
        return(self.infos)
        
    def parse_(self):
        import zipfile
        exts = [DatasetFactory().data_types[x]['fmts'] for x in DatasetFactory().data_types.keys()]
        exts = [x for y in exts for x in y]
        datalist = []
        if self.fn.split('.')[-1].lower() == 'zip':
            with zipfile.ZipFile(self.fn) as z:
                zfs = z.namelist()
                for ext in exts:
                    for zf in zfs:
                        if ext == zf.split('.')[-1]:
                            datalist.append(os.path.basename(zf))
                            
        for this_data in datalist:
            this_line = utils.p_f_unzip(self.fn, [this_data])[0]
            data_set = DatasetFactory(
                this_line,
                weight=self.weight,
                parent=self,
                src_region=self.region,
                metadata=copy.deepcopy(self.metadata),
                src_srs=self.src_srs,
                dst_srs=self.dst_srs,
                verbose=self.verbose
            ).acquire()
            if data_set is not None and data_set.valid_p(
                    fmts=DatasetFactory.data_types[data_set.data_format]['fmts']
            ):                            
                if self.region is not None and self.region.valid_p(check_xy=True):
                    try:
                        inf_region = regions.Region().from_string(
                            data_set.infos['wkt']
                        )
                    except:
                        try:
                            inf_region = regions.Region().from_list(
                                data_set.infos['minmax']
                            )
                        except:
                            inf_region = self.region.copy()
                            
                    inf_region.wmin = data_set.weight
                    inf_region.wmax = data_set.weight
                    if regions.regions_intersect_p(inf_region, self.region):
                        for ds in data_set.parse():
                            self.data_entries.append(ds)
                            yield(ds)
                else:
                    for ds in data_set.parse():
                        self.data_entries.append(ds)
                        yield(ds)
            
            utils.remove_glob('{}*'.format(this_data))

    def yield_xyz(self):
        for ds in self.parse_():
            for xyz in ds.yield_xyz():
                yield(xyz)

## ==============================================
## dlim Fetcher dataset class
## ==============================================
class Fetcher(datasets.ElevationDataset):
    """The fetches dataset type.

This is used in waffles/dlim for on-the-fly remote data
parsing and processing.
"""
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.remote=True
        self.metadata['name'] = self.fn
        self.fetch_module = fetches.FetchesFactory(
            mod=self.fn,
            src_region=self.region,
            dst_srs=self.dst_srs,
            verbose=self.verbose,
            weight=self.weight,
            x_inc=self.x_inc,
            y_inc=self.y_inc
        ).acquire()

        if self.fetch_module is None:
            pass
        
    def generate_inf(self, callback=lambda: False):
        """generate a infos dictionary from the Fetches dataset"""

        self.infos['name'] = self.fn
        self.infos['hash'] = None
        self.infos['numpts'] = 0
        if self.region is None:
            #self.region = self.fetch_module.region
            self.region = regions.Region().from_list([-180,180,-90,90])
            
        self.infos['minmax'] = self.region.export_as_list()
        self.infos['wkt'] = self.region.export_as_wkt()
        self.infos['format'] = self.data_format
        return(self.infos)

    def parse_(self):
        self.fetch_module.run()
        for result in self.fetch_module.results:
            data_set = DatasetFactory(
                result[0],
                weight=self.weight,
                parent=self,
                src_region=self.region,
                metadata=copy.deepcopy(self.metadata),
                src_srs=self.src_srs,
                dst_srs=self.dst_srs,
                verbose=self.verbose,
            ).acquire(remote=True)
            
            yield(data_set)
    
    def yield_xyz(self):
        for xyz in self.fetch_module.yield_results_to_xyz():
            yield(xyz)

    def yield_array(self):
        """yield the data array
        set x_inc and y_inc to warp data to target resolution
        """
        
        for arr in self.fetch_module.yield_results_to_array():
            yield(arr)
            
## ==============================================
## Dataset generator
## ==============================================
class DatasetFactory:

    data_types = {
        -1: {'name': 'datalist',
             'fmts': ['datalist', 'mb-1', 'dl'],
             'opts': '< -1 >',
             'class': Datalist,
             },
        -2: {'name': 'zip',
             'fmts': ['zip'],
             'opts': '< -2 >',
             'class': ZIPlist,
             },
        168: {'name': 'xyz',
              'fmts': ['xyz', 'csv', 'dat', 'ascii', 'txt'],
              'opts': '< 168:delim=None:xpos=0:ypos=1:zpos=2:skip=0:x_scale=1:y_scale=1:z_scale=1:x_offset=0:y_offset=0 >',
              'class': datasets.XYZFile,
              },
        200: {'name': 'raster',
              'fmts': ['tif', 'tiff', 'img', 'grd', 'nc', 'vrt'],
              'opts': '< 200:open_options="OPT1=KEY"/"OPT2=KEY2" >',
              'class': datasets.RasterFile,
              },
        201: {'name': 'bag',
              'fmts': ['bag'],
              'opts': '< 201:open_options="OPT1=KEY"/"OPT2=KEY2" >',
              'class': datasets.BAGFile,
              },
        300: {'name': 'las',
              'fmts': ['las', 'laz'],
              'opts': '< 300:classes=0/2/29/40 >',
              'class': datasets.LASFile,
              },
        301: {'name': 'mbs',
              'fmts': ['fbt'],
              'opts': '< 301 >',
              'class': datasets.MBSParser,
              },
        -11: {'name': 'fetches',
              'fmts': [
                  'gmrt',
                  'multibeam',
                  'ehydro',
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
                  'nasadem',
                  'tides',
                  'vdatum',
                  'earthdata',
                  'bluetopo',
                  'hydrolakes'
              ],
              'opts': '< -11 >',
              'class': Fetcher,
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
            src_srs=None,
            dst_srs=None,
            x_inc=None,
            y_inc=None,
            sample_alg='bilinear',
            metadata={
                'name':'xyz dataset',
                'title':None,
                'source':None,
                'date':None,
                'data_type':None,
                'resolution':None,
                'hdatum':None,
                'vdatum':None,
                'url':None
            },
            parent=None,
            src_region=None,
            cache_dir=None,
            verbose=False,
            remote=False
    ):
        self.fn = fn
        self.data_format = data_format
        self.weight = weight
        self.src_srs = src_srs
        self.dst_srs = dst_srs
        self.metadata = copy.deepcopy(metadata)
        self.region = src_region
        self.parent = parent
        self.verbose = verbose
        self.ds_args = {}
        self.x_inc = x_inc
        self.y_inc = y_inc
        self.sample_alg = sample_alg
        self.remote = remote
        self.cache_dir = utils.cudem_cache() if cache_dir is None else cache_dir
        
        self.parse_fn()
            
    def parse_fn(self):
        """parse the datalist entry line"""
        
        if self.fn is None: return(self)
        if os.path.exists(self.fn):
            self.guess_data_format()
            self.metadata['name'] = self.fn.split('.')[0]
            return(self.fn)
        
        this_entry = re.findall(r'[^"\s]\S*|".+?"', self.fn.rstrip())

        try:
            entry = [utils.str_or(x) if n == 0 else utils.int_or(x) if n < 2 else utils.float_or(x) if n < 3 else utils.str_or(x) \
                     for n, x in enumerate(this_entry)]
        except Exception as e:
            utils.echo_error_msg('could not parse entry {}'.format(self.fn))
            return(self)

        ## ==============================================
        ## data format
        ## guess format based on fn if not specified and
        ## parse format for dataset specific opts.
        ## ==============================================
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

        else:
            opts = this_entry[1].split(':')
            if len(opts) > 1:
                self.ds_args = utils.args2dict(list(opts[1:]), {})
                this_entry[1] = opts[0]
            else:
                self.ds_args = {}

        ## ==============================================
        ## weight
        ## inherit weight of parent
        ## ==============================================
        if len(entry) < 3:
            entry.append(1)
        elif entry[2] is None:
            entry[2] = 1

        if self.parent is not None:
            if self.weight is not None:
                self.weight *= entry[2]
        else:
            if self.weight is not None:
                self.weight = entry[2]

        ## ==============================================
        ## title
        ## ==============================================
        if len(entry) < 4:
            entry.append(self.metadata['title'])
        else:
            self.metadata['title'] = entry[3]

        ## ==============================================
        ## source
        ## ==============================================
        if len(entry) < 5:
            entry.append(self.metadata['source'])
        else:
            self.metadata['source'] = entry[4]

        ## ==============================================
        ## date
        ## ==============================================
        if len(entry) < 6:
            entry.append(self.metadata['date'])
        else:
            self.metadata['date'] = entry[5]

        ## ==============================================
        ## data type
        ## ==============================================
        if len(entry) < 7:
            entry.append(self.metadata['data_type'])
        else:
            self.metadata['data_type'] = entry[6]

        ## ==============================================
        ## resolution
        ## ==============================================
        if len(entry) < 8:
            entry.append(self.metadata['resolution'])
        else:
            self.metadata['resolution'] = entry[7]

        ## ==============================================
        ## hdatum
        ## ==============================================
        if len(entry) < 9:
            entry.append(self.metadata['hdatum'])
        else:
            self.metadata['hdatum'] = entry[8]

        ## ==============================================
        ## vdatum
        ## ==============================================
        if len(entry) < 10:
            entry.append(self.metadata['vdatum'])
        else:
            self.metadata['vdatum'] = entry[9]

        ## ==============================================
        ## url
        ## ==============================================
        if len(entry) < 11:
            entry.append(self.metadata['url'])
        else:
            self.metadata['url'] = entry[10]

        ## ==============================================
        ## file-name
        ## ==============================================
        if self.parent is None or entry[1] == -11:
            self.fn = entry[0]
        else:
            self.fn = os.path.join(
                os.path.dirname(self.parent.fn), entry[0]
            )
        
        self.data_format = entry[1]
        self.guess_data_format()
        return(self)

    def guess_data_format(self):
        """guess a data format based on the file-name"""

        if self.data_format is not None:
            return
        
        if self.fn is not None:
            for key in self.data_types.keys():
                if self.fn.split('.')[-1] in self.data_types[key]['fmts']:
                    self.data_format = key
                    break

    def add_data_type(self, type_def={}):
        for key in type_def.keys():
            self.data_types[key] = type_def[key]

    def acquire(self, **kwargs):
        return(
            self.data_types[self.data_format]['class'](
                fn=self.fn,
                data_format=self.data_format,
                weight=self.weight,
                src_region=self.region,
                metadata=copy.deepcopy(self.metadata),
                src_srs=self.src_srs,
                dst_srs=self.dst_srs,
                parent=self.parent,
                verbose=self.verbose,
                x_inc=self.x_inc,
                y_inc=self.y_inc,
                sample_alg=self.sample_alg,
                **self.ds_args,
                **kwargs
            )
        )

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

usage: {cmd} [ -acdghijqwEPRW [ args ] ] DATALIST,FORMAT,WEIGHT ...

Options:
  -R, --region\t\tRestrict processing to the desired REGION 
\t\t\tWhere a REGION is xmin/xmax/ymin/ymax[/zmin/zmax[/wmin/wmax]]
\t\t\tUse '-' to indicate no bounding range; e.g. -R -/-/-/-/-10/10/1/-
\t\t\tOR an OGR-compatible vector file with regional polygons. 
\t\t\tWhere the REGION is /path/to/vector[:zmin/zmax[/wmin/wmax]].
\t\t\tIf a vector file is supplied, will use each region found therein.
  -E, --increment\tBlock data to INCREMENT in native units.
\t\t\tWhere INCREMENT is x-inc[/y-inc]
  -P, --s_srs\t\tSet the SOURCE projection.
  -W, --t_srs\t\tSet the TARGET projection.

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
  % {cmd} -R my_region.shp my_data.xyz -w -s_srs epsg:4326 -t_srs epsg:3565 > my_data_3565.xyz

CIRES DEM home page: <http://ciresgroups.colorado.edu/coastalDEM>\
""".format(cmd=os.path.basename(sys.argv[0]), 
           dl_version=cudem.__version__,
           dl_formats=_datalist_fmts_short_desc())

def datalists_cli(argv=sys.argv):
    """run datalists from command-line

    See `datalists_cli_usage` for full cli options.
    """

    dls = []
    src_srs = None
    dst_srs = None
    i_regions = []
    these_regions = []
    xy_inc = [None, None]
    want_weights = False
    want_inf = False
    want_list = False
    want_glob = False
    want_archive = False
    want_verbose = True
    want_region = False
    want_csv = False
    want_json = False
    want_datalists=False
    want_separate=False
    
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
        elif arg == '-s_srs' or arg == '--s_srs' or arg == '-P':
            src_srs = argv[i + 1]
            i = i + 1
        elif arg == '-t_srs' or arg == '--t_srs' or arg == '-W':
            dst_srs = argv[i + 1]
            i = i + 1
        elif arg == '--separate' or arg == '-s':
            want_separate = True
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
        elif arg == '--datalists' or arg == '-d':
            want_datalists = True
        elif arg == '--csv' or arg == '-c':
            want_csv = True
        elif arg == '--json' or arg == '-j':
            want_json = True
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

    if len(xy_inc) < 2:
        xy_inc.append(xy_inc[0])
        
    elif len(xy_inc) == 0:
        xy_inc = [None, None]

    if want_glob:
        import glob
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
            sys.stderr.write(datalists_usage)
            utils.echo_error_msg('you must specify some type of data')
        else:
            this_datalist = init_data(dls, this_region, src_srs, dst_srs, xy_inc, 'bilinear', want_verbose)
            if this_datalist is not None and this_datalist.valid_p(
                    fmts=DatasetFactory.data_types[this_datalist.data_format]['fmts']
            ):
                if not want_weights:
                    this_datalist.weight = None

                if want_inf:
                    print(this_datalist.inf())
                elif want_list:
                    this_datalist.echo()            
                elif want_archive:
                    [x for x in this_datalist.archive_xyz()]
                elif want_region:
                    this_inf = this_datalist.inf()
                    this_region = regions.Region().from_list(this_inf['minmax'])
                    if dst_srs is not None:
                        if src_srs is not None:
                            this_region.src_srs = src_srs
                            this_region.warp(dst_srs)
                        elif 'src_srs' in this_inf and this_inf['src_srs'] is not None:
                            this_region.src_srs = this_inf['src_srs']
                            this_region.warp(dst_srs)

                    print(this_region.format('gmt'))
                    #print(regions.Region().from_list(this_datalist.inf()['minmax']).format('gmt'))
                elif want_csv:
                    this_datalist.parse_data_lists()
                    for x in this_datalist.data_lists.keys():
                        this_datalist.data_entries = this_datalist.data_lists[x]['data']
                        p = this_datalist.data_lists[x]['parent']
                        print(
                            '|'.join(
                                [
                                    '"{}"'.format(str(y)) for y in [
                                        x,
                                        p.metadata['title'] if p.metadata['title'] is not None else x,
                                        p.metadata['source'],
                                        p.metadata['date'],
                                        p.metadata['data_type'],
                                        p.metadata['resolution'],
                                        p.metadata['hdatum'],
                                        p.metadata['vdatum'],
                                        p.metadata['url']
                                    ]
                                ]
                            )
                        )
                elif want_datalists:
                    #import json
                    #j = open('{}.json'.format(this_datalist.metadata['name']), 'w')
                    this_datalist.parse_data_lists()
                    # for x in this_datalist.data_lists.keys():
                    #     p = this_datalist.data_lists[x]['parent']

                    #     out_json = {
                    #         "Name": x,
                    #         "Title": p.title if p.title is not None else x,
                    #         "Source": p.source,
                    #         "Date": p.date,
                    #         "DataType": p.data_type,
                    #         "Resolution": p.resolution,
                    #         "HDatum": p.hdatum,
                    #         "VDatum": p.vdatum,
                    #         "URL": p.url
                    #     }
                    #     j.write(json.dumps(out_json))
                    #     j.write('\n')
                    # j.close()
                    
                    for x in this_datalist.data_lists.keys():
                        p = this_datalist.data_lists[x]['parent']
                        print('{} ({})|{}'.format(p.metadata['title'], p.metadata['name'], p.weight))
                        #print(xdl.data_lists[x]['parent'].echo_())
                        #print('{}'.format(this_datalist.data_lists[x]['parent'].fn))
                        #print(this_datalist.data_lists[x])

                        #print('{} ({})|{}'.format(this_datalist.data_lists[x]['parent'].metadata['title'], this_datalist.data_lists[x]['parent'].metadata['name'], this_datalist.data_lists[x]['parent'].weight))
                else:
                    if want_separate:
                        try:
                            for this_entry in this_datalist.parse():
                                this_entry.dump_xyz()
                        except KeyboardInterrupt:
                            utils.echo_error_msg('Killed by user')
                            break
                        except BrokenPipeError:
                            utils.echo_error_msg('Pipe Broken')
                            break
                    else:
                        try:
                            this_datalist.dump_xyz()
                        except KeyboardInterrupt:
                            utils.echo_error_msg('Killed by user')
                            break
                        except BrokenPipeError:
                            utils.echo_error_msg('Pipe Broken')
                            break
### End
