### datalistfile.py - DataLists IMproved
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
##
## datalistfile.py is part of CUDEM
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
##
### Code:

import os
import copy
import numpy as np
import pandas as pd
from osgeo import ogr

from cudem import utils
from cudem import regions
from cudem import gdalfun
from cudem import srsfun
from cudem.datalists.dlim import ElevationDataset

class Points(ElevationDataset):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.metadata['name'] = 'points'
        self.data_format = -4

    def yield_points(self):
        points = self.fn
        if isinstance(points, np.ndarray):
            points = np.rec.fromrecords(points, names='x, y, z')
        elif isinstance(points, np.core.records.recarray):
            pass # already recarray
        elif isinstance(points, pd.DataFrame):
            points = points.to_records(index=False)
        
        yield points

        
class Scratch(ElevationDataset):
    """Scratch Dataset
    Process a python list of valid dataset entries
    """
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.metadata['name'] = 'scratch'
        self.data_format = -3

        
    def generate_inf(self, make_grid=True, make_block_mean=False, block_inc=None):
        """Generate INF for the scratch list."""
        
        _region_backup = self.region
        self.region = None # Unset to scan all children
        
        point_count = 0
        out_regions = []
        out_srs = []
        valid_children = []

        ## Pass 1: Gather Stats from Children
        for entry in self.parse():
            ## Trigger child INF generation
            entry.inf(make_grid=False, make_block_mean=False)
            
            if entry.infos.minmax:
                valid_children.append(entry)
                entry_region = regions.Region().from_list(entry.infos.minmax)
                
                if entry.src_srs is not None:
                    out_srs.append(entry.src_srs)
                    if self.dst_srs is not None:
                        entry_region.src_srs = entry.src_srs
                        entry_region.warp(self.dst_srs)

                if entry_region.valid_p():
                    out_regions.append(entry_region)
                    point_count += entry.infos.numpts

        ## Aggregate Region
        master_region = None
        for r in out_regions:
            if master_region is None: master_region = r
            else: master_region = regions.regions_merge(master_region, r)
            
        if master_region is not None:
            self.infos.minmax = master_region.export_as_list(include_z=True)
            self.infos.wkt = master_region.export_as_wkt()
            self.infos.numpts = point_count

            ## Determine SRS
            if self.infos.src_srs is None:
                if self.src_srs is not None:
                    self.infos.src_srs = self.src_srs
                elif out_srs:
                    if all(x == out_srs[0] for x in out_srs):
                        self.infos.src_srs = out_srs[0]
                    self.src_srs = self.infos.src_srs

            ## Pass 2: Grids from Children
            if (make_grid or make_block_mean) and point_count > 0:
                self._generate_grids_from_children(
                    valid_children, master_region, 
                    make_grid, make_block_mean, block_inc
                )

        self.region = _region_backup
        return self.infos

    
    def parse(self):
        if isinstance(self.fn, list):
            for this_ds in self.fn:
                if this_ds is not None and this_ds.valid_p():
                    this_ds.initialize()
                    for ds in this_ds.parse():
                        self.data_entries.append(ds)
                        yield ds

                        
class Datalist(ElevationDataset):
    """Representing a datalist parser (Extended MB-System style).
    Can parse recursive datalists and sidecar JSON indices.
    """

    _datalist_json_cols = [
        'path', 'format', 'weight', 'uncertainty', 'name', 'title', 'source',
        'date', 'data_type', 'resolution', 'hdatum', 'vdatum',
        'url', 'mod_args'
    ]
    
    def __init__(self, fmt=None, **kwargs):
        super().__init__(**kwargs)
        self.kwargs = kwargs

        
    def _init_datalist_vector(self):
        """Initialize the datalist geojson vector index."""
        
        self.dst_layer = f'{self.fn}'
        self.dst_vector = f'{self.dst_layer}.json'

        utils.remove_glob(self.dst_vector)
        
        self.ds = ogr.GetDriverByName('GeoJSON').CreateDataSource(self.dst_vector)
        srs = srsfun.osr_srs(self.src_srs)
        
        if self.ds is not None: 
            self.layer = self.ds.CreateLayer(self.dst_layer, srs, ogr.wkbMultiPolygon)
            for f in self._datalist_json_cols:
                self.layer.CreateField(ogr.FieldDefn(f, ogr.OFTString))
            return 0
        return -1

    
    def _create_entry_feature(self, entry, entry_region):
        """Create a datalist entry feature in the geojson."""
        
        entry_path = os.path.abspath(entry.fn) if not entry.remote else entry.fn
        entry_fields = [
            entry_path,
            entry.data_format,
            entry.weight,
            entry.uncertainty,
            entry.metadata.get('name'),
            entry.metadata.get('title'),
            entry.metadata.get('source'),
            entry.metadata.get('date'),
            entry.metadata.get('data_type'),
            entry.metadata.get('resolution'),
            entry.metadata.get('hdatum'),
            entry.metadata.get('vdatum'),
            entry.metadata.get('url'),
            utils.dict2args(entry.params.get('mod_args', {}))
        ]
        
        out_feat = ogr.Feature(self.layer.GetLayerDefn())
        out_feat.SetGeometry(ogr.CreateGeometryFromWkt(entry_region.export_as_wkt()))
        
        for i, f in enumerate(self._datalist_json_cols):
            out_feat.SetField(f, str(entry_fields[i]))
            
        self.layer.CreateFeature(out_feat)

        
    def generate_inf(self, make_grid=True, make_block_mean=False, block_inc=None):
        """Generate metadata (INF) for the datalist.
        
        Iterates through datalist entries, aggregates bounds, and generates 
        a JSON index. Optionally generates aggregated grids from children.
        """
        
        _region_backup = self.region
        self.region = None # Unset region to scan all entries
        
        self.infos.file_hash = self.infos.generate_hash()
        
        point_count = 0
        out_regions = []
        out_srs = []
        valid_children = []

        ## Initialize JSON Vector Index
        json_init = self._init_datalist_vector()

        ## Parse Entries, Gather Stats, Populate JSON
        for entry in self.parse():
            ## Generate Child INF (recursively if needed)
            entry.inf(make_grid=False, make_block_mean=False)
            
            if entry.infos.minmax:
                valid_children.append(entry)
                entry_region = regions.Region().from_list(entry.infos.minmax)
                
                ## Transform bounds if needed
                if entry.src_srs is not None:
                    out_srs.append(entry.src_srs)
                    if self.dst_srs is not None:
                        entry_region.src_srs = entry.src_srs
                        entry_region.warp(self.dst_srs)

                if entry_region.valid_p():
                    out_regions.append(entry_region)
                    point_count += entry.infos.numpts
                    
                    ## Add to JSON Index
                    if json_init == 0:
                        self._create_entry_feature(entry, entry_region)

        ## Close JSON DataSource
        self.ds = self.layer = None

        ## Aggregate Regions
        master_region = None
        for r in out_regions:
            if master_region is None: master_region = r
            else: master_region = regions.regions_merge(master_region, r)
            
        if master_region is not None:
            self.infos.minmax = master_region.export_as_list(include_z=True)
            self.infos.wkt = master_region.export_as_wkt()
            self.infos.numpts = point_count

            ## Determine Master SRS
            if self.infos.src_srs is None:
                if self.src_srs is not None:
                    self.infos.src_srs = self.src_srs
                elif out_srs:
                    if all(x == out_srs[0] for x in out_srs):
                        self.infos.src_srs = out_srs[0]
                    self.src_srs = self.infos.src_srs

            ## Grids from Children
            if (make_grid or make_block_mean) and point_count > 0:
                self._generate_grids_from_children(
                    valid_children, master_region, 
                    make_grid, make_block_mean, block_inc
                )

        self.region = _region_backup
        return self.infos

    
    def parse_json(self):
        """Parse datalist using the JSON sidecar index."""

        from cudem.datalists.dlim import DatasetFactory
        
        if not os.path.exists(f'{self.fn}.json'):
            ## Fallback to standard parse if JSON missing
            for ds in self.parse():
                yield ds
            return

        ds = ogr.Open(f'{self.fn}.json')
        if ds is None: return

        layer = ds.GetLayer()
        ## Spatial Filter
        if self.region is not None:
            check_region = self.transform['trans_region'] if self.transform['trans_region'] else self.region
            layer.SetSpatialFilter(ogr.CreateGeometryFromWkt(check_region.export_as_wkt()))

        for feat in layer:
            ## Metadata Filter (Weight/Uncertainty)
            if self.region:
                wr = self.region.w_region()
                ur = self.region.u_region()
                w, u = float(feat.GetField('weight')), float(feat.GetField('uncertainty'))
                
                if (wr[0] is not None and w < wr[0]) or (wr[1] is not None and w > wr[1]): continue
                if (ur[0] is not None and u < ur[0]) or (ur[1] is not None and u > ur[1]): continue

            ## Reconstruct Entry
            ds_args_str = feat.GetField('mod_args')
            try:
                ds_args = feat.GetField('mod_args')
                data_set_args = utils.args2dict(list(ds_args.split(':')), {})
                # for kpam in list(data_set_args.keys()):
                #     if kpam in self.__dict__:
                #         kval = data_set_args[kpam]
                #         self.__dict__[kpam] = kval
                #         del data_set_args[kpam]
                
                ds_args = utils.dict2args(data_set_args)
                # for kpam, kval in data_set_args.items():
                #     if kpam in self.__dict__:
                #         utils.echo_msg(kpam)
                #         self.__dict__[kpam] = kval
                #         #del data_set_args[kpam]
            except:
                ds_args = None
                data_set_args = {}
            
            md = copy.deepcopy(self.metadata)
            for key in self.metadata.keys():
                val = feat.GetField(key)
                if val: md[key] = val

            data_mod = f'"{feat.GetField("path")}" {feat.GetField("format")} {ds_args_str or ""} {feat.GetField("weight")} {feat.GetField("uncertainty")}'
            
            data_set = DatasetFactory(
                **self._set_params(mod=data_mod, metadata=md)#, **data_set_args)
            )._acquire_module()
            
            if data_set and data_set.valid_p():
                data_set.initialize()
                for ds in data_set.parse():
                    self.data_entries.append(ds)
                    yield ds

    def parse(self):
        """Parse the datalist text file."""

        from cudem.datalists.dlim import DatasetFactory
        
        if os.path.exists(self.fn):
            ## Get the number of lines in the datalist
            with open(self.fn, 'r') as f:
                count = sum(1 for _ in f)
                
            with open(self.fn, 'r') as f:
                with utils.ccp(
                        total=count,
                        desc=f'Parsing datalist {self.fn}...',
                        leave=False
                ) as pbar:
                    for line in f:
                        pbar.update()
                        line = line.strip()
                        if line and not line.startswith('#'):
                            md = copy.deepcopy(self.metadata)
                            md['name'] = utils.fn_basename2(os.path.basename(self.fn))

                            ds_params = self._set_params(
                                mod=line,
                                metadata=md,
                                src_srs=self.src_srs,
                                parent=self,
                                fn=None,
                            )
                            data_set = DatasetFactory(**ds_params)._acquire_module()
                            utils.echo_debug_msg(f'Parsed entry from line ({line}): {data_set}, {data_set.fn}')
                            # data_set = DatasetFactory(
                            #     mod=line,
                            #     metadata=md,
                            #     src_srs=self.src_srs,
                            #     parent=self
                            # )._acquire_module()

                            if data_set and data_set.valid_p(fmts=DatasetFactory._modules[data_set.data_format]['fmts']):
                                data_set.initialize()

                                ## Spatial Filter
                                if self.region is not None and self.region.valid_p(check_xy=True):
                                    ## Load INF to check bounds
                                    data_set.inf(check_hash=False)

                                    check_region = data_set.transform['trans_region']
                                    if check_region is None and data_set.infos.minmax:
                                        check_region = regions.Region().from_list(data_set.infos.minmax)

                                    if check_region and not regions.regions_intersect_p(check_region, self.region):
                                        continue

                                for ds in data_set.parse():
                                    self.data_entries.append(ds)
                                    yield ds

        elif self.data_entries:
            for data_set in self.data_entries:
                for ds in data_set.parse():
                    yield ds
        else:
            if self.verbose:
                utils.echo_warning_msg(f'could not open datalist/entry {self.fn}')

### End
