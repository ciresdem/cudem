### datalistfile.py - DataLists IMproved
##
## Copyright (c) 2012 - 2026 Regents of the University of Colorado
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
## DatalistFile: The CUDEM Datalist Parsing and Indexing Engine.
##
## This module defines the core classes for managing CUDEM datalists—recursive
## lists of geospatial data entries. It handles the parsing, validation, and
## metadata generation for diverse data sources.
##
##   * Recursive Parsing: Efficiently traverses nested datalists to build
##     flat lists of processable data entries.
##   * Spatial Indexing: Generates sidecar spatial indices to optimize
##     data discovery and region-based filtering:
##       - GeoJSON (.json): Standard vector footprints for OGR/GIS compatibility.
##       - H3 Hexagons (.h3.json): Hierarchical hexagonal indices for fast O(1)
##         spatial lookups and variable-resolution coverage masks.
##
### Code:

import os
import sys
import copy
import json
import argparse
import numpy as np
import pandas as pd
from osgeo import ogr

from cudem import utils
from cudem import regions
from cudem import gdalfun
from cudem import srsfun
from cudem.datalists.dlim import ElevationDataset

## Optional H3 Dependency
try:
    import h3
    HAS_H3 = True
except ImportError:
    HAS_H3 = False

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

        ## Gather Stats from Children
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

            ## Grids from Children
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
    Can parse recursive datalists and sidecar JSON/H3 indices.
    """

    _datalist_json_cols = [
        'path', 'format', 'weight', 'uncertainty', 'name', 'title', 'source',
        'date', 'data_type', 'resolution', 'hdatum', 'vdatum',
        'url', 'mod_args'
    ]
    
    def __init__(self, fmt=None, h3_res=6, want_h3=True, want_json=True, **kwargs):
        super().__init__(**kwargs)
        self.kwargs = kwargs
        self.h3_res = h3_res # Default H3 Resolution (Res 6 ~= 36 km2)
        self.want_h3 = want_h3
        self.want_json = want_json
        
        ## H3 Index Storage
        ## Structure: { "entries": [{meta}], "index": { "h3_cell": [entry_idx, ...] } }
        self.h3_data = {"resolution": self.h3_res, "entries": [], "index": {}}
        self.data_format = -1
        
        
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

    
    def _update_h3_index(self, entry, entry_region, entry_idx):
        """Update the H3 Inverted Index with a new entry."""
        if not HAS_H3: return

        ## Store Entry Metadata
        entry_meta = {
            'path': os.path.abspath(entry.fn) if not entry.remote else entry.fn,
            'format': entry.data_format,
            'weight': entry.weight,
            'uncertainty': entry.uncertainty,
            'mod_args': utils.dict2args(entry.params.get('mod_args', {}))
            # Add other metadata fields as needed
        }
        self.h3_data['entries'].append(entry_meta)

        ## Convert Region to GeoJSON Polygon
        ## H3 expects (lat, lon) ? No, h3-py v4+ uses (lat, lon) tuples or GeoJSON (lon, lat)
        ## We will use h3.polygon_to_cells which expects GeoJSON-like dictionary
        if not entry_region.valid_p(): return

        geo_json_poly = {
            "type": "Polygon",
            "coordinates": [[
                [entry_region.xmin, entry_region.ymin],
                [entry_region.xmin, entry_region.ymax],
                [entry_region.xmax, entry_region.ymax],
                [entry_region.xmax, entry_region.ymin],
                [entry_region.xmin, entry_region.ymin]
            ]]
        }

        try:
            ## Polyfill region to H3 cells
            cells = h3.polygon_to_cells(geo_json_poly, self.h3_res)
            
            ## Update Inverted Index: Cell -> List of Entry IDs
            for cell in cells:
                if cell not in self.h3_data['index']:
                    self.h3_data['index'][cell] = []
                self.h3_data['index'][cell].append(entry_idx)
                
        except Exception as e:
            utils.echo_warning_msg(f"H3 Indexing Error for {entry.fn}: {e}")

            
    def _save_h3_index(self):
        """Save the H3 index to a sidecar JSON file."""
        
        if not HAS_H3 or not self.h3_data['entries']: return
        
        h3_fn = f'{self.fn}.h3.json'
        try:
            with open(h3_fn, 'w') as f:
                json.dump(self.h3_data, f)
        except Exception as e:
            utils.echo_warning_msg(f"Failed to save H3 index: {e}")

            
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
        a JSON/H3 index. Optionally generates aggregated grids from children.
        """

        _region_backup = self.region
        self.region = None # Unset region to scan all entries
        
        self.infos.file_hash = self.infos.generate_hash()
        
        point_count = 0
        out_regions = []
        out_srs = []
        valid_children = []

        json_init = -1
        if self.want_json:
            ## Initialize JSON Vector Index
            json_init = self._init_datalist_vector()
        
        ## Reset H3 Index
        self.h3_data = {"resolution": self.h3_res, "entries": [], "index": {}}
        entry_idx = 0

        #self.data_format = -1
        #self.metadata['name'] = self.fn
        #self.initialize()
        #utils.echo_msg(self.metadata)
        ## Parse Entries, Gather Stats, Populate JSON/H3
        for entry in self.parse():
            utils.echo_msg_bold(entry)
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
                    
                    ## Add to OGR JSON Index
                    if json_init == 0:
                        self._create_entry_feature(entry, entry_region)
                        
                    ## Add to H3 Index
                    if HAS_H3 and self.want_h3:
                        self._update_h3_index(entry, entry_region, entry_idx)
                        entry_idx += 1

        ## Close JSON DataSource
        self.ds = self.layer = None
        
        ## Save H3 Index
        if HAS_H3 and self.want_h3:
            self._save_h3_index()

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

    
    def parse_h3(self):
        """Parse datalist using the H3 sidecar index.
        
        Performs O(1) spatial filtering using H3 cells.
        """
        
        from cudem.datalists.dlim import DatasetFactory
        
        h3_fn = f'{self.fn}.h3.json'
        if not HAS_H3 or not os.path.exists(h3_fn):
            ## Fallback to standard parse
            for ds in self.parse():
                yield ds
            return

        ## Load Index
        try:
            with open(h3_fn, 'r') as f:
                h3_index = json.load(f)
        except:
            return

        ## Determine Matching Entries
        matched_indices = set()
        
        if self.region is None:
            ## If no spatial filter, yield all
            matched_indices = set(range(len(h3_index['entries'])))
        else:
            ## Convert Query Region to H3 Cells
            ## We use the resolution specified in the index file
            res = h3_index.get('resolution', 6)
            
            ## Expand region slightly to catch edge hexes
            ## Or assume the region is a bbox
            geo_json_poly = {
                "type": "Polygon",
                "coordinates": [[
                    [self.region.xmin, self.region.ymin],
                    [self.region.xmin, self.region.ymax],
                    [self.region.xmax, self.region.ymax],
                    [self.region.xmax, self.region.ymin],
                    [self.region.xmin, self.region.ymin]
                ]]
            }
            
            try:
                ## Get cells for query region
                query_cells = h3.polygon_to_cells(geo_json_poly, res)
                
                ## Check intersection with index
                ## (Look up each query cell in the index)
                idx_map = h3_index['index']
                for cell in query_cells:
                    if cell in idx_map:
                        ## Add all entry IDs found in this cell
                        matched_indices.update(idx_map[cell])
                        
            except Exception as e:
                utils.echo_warning_msg(f"H3 Query Error: {e}")
                ## Fallback to yielding all might be safer? 
                ## Or yielding none? Let's just return to avoid bad data.
                return

        ## Yield Matched Entries
        for i in sorted(list(matched_indices)):
            entry_meta = h3_index['entries'][i]
            
            ## Reconstruct Entry Object
            ## We construct a mock mod string: "path format mod_args weight uncertainty"
            data_mod = (f'"{entry_meta["path"]}" {entry_meta["format"]} '
                        f'{entry_meta.get("mod_args", "")} '
                        f'{entry_meta["weight"]} {entry_meta["uncertainty"]}')
            
            ## Reconstruct Metadata (Deep copy from parent, apply entry overrides if saved)
            md = copy.deepcopy(self.metadata)
            
            data_set = DatasetFactory(
                **self._set_params(mod=data_mod, metadata=md)
            )._acquire_module()
            
            if data_set and data_set.valid_p():
                data_set.initialize()
                for ds in data_set.parse():
                    self.data_entries.append(ds)
                    yield ds

                    
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
                ds_args = utils.dict2args(data_set_args)
            except:
                ds_args = None
            
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

        # # [Priority Check] Check for H3 Index First
        # if HAS_H3 and os.path.exists(f'{self.fn}.h3.json'):
        #      for ds in self.parse_h3():
        #          yield ds
        #      return

        # # [Secondary Check] Check for OGR JSON Index
        # if os.path.exists(f'{self.fn}.json'):
        #     for ds in self.parse_json():
        #         yield ds
        #     return

        # Standard Text Parse
        if os.path.exists(self.fn):
            ## Get the number of lines in the datalist
            count = utils.count_data_lines(self.fn)
            
            with open(self.fn, 'r') as f:
                with utils.ccp(
                        total=count,
                        desc=f'Parsing datalist {self.fn}...',
                        leave=self.verbose
                ) as pbar:
                    for line in f:
                        line = line.strip()
                        if line and not line.startswith('#'):
                            pbar.update()
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


def datalist2json_cli():
    """Command-line interface for datalist2json."""
    
    parser = argparse.ArgumentParser(
        description="datalist2json: Create spatial indices (JSON/H3) for datalists.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        'datalist',
        help="Input datalist file path."
    )
    
    ## Indexing Options
    parser.add_argument(
        '--json',
        dest='want_json',
        action='store_true',
        default=True,
        help="Generate standard OGR/GeoJSON spatial index (.json)."
    )
    parser.add_argument(
        '--no-json',
        dest='want_json',
        action='store_false',
        help="Do not generate standard OGR/GeoJSON spatial index."
    )
    
    parser.add_argument(
        '--h3',
        dest='want_h3',
        action='store_true',
        default=True,
        help="Generate H3 hexagonal spatial index (.h3.json)."
    )
    parser.add_argument(
        '--no-h3',
        dest='want_h3',
        action='store_false',
        help="Do not generate H3 hexagonal spatial index."
    )
    parser.add_argument(
        '--h3-res',
        type=int,
        default=6,
        help="H3 resolution level (0-15). 6 is ~36km2, 9 is ~0.1km2."
    )

    ## General Options
    parser.add_argument(
        '-q', '--quiet',
        action='store_true',
        help="Suppress output."
    )

    args = parser.parse_args()

    ## Validate Input
    if not os.path.exists(args.datalist):
        utils.echo_error_msg(f"Datalist not found: {args.datalist}")
        sys.exit(1)

    ## Initialize Datalist
    dl = Datalist(
        fn=args.datalist,
        want_json=args.want_json,
        want_h3=args.want_h3,
        h3_res=args.h3_res,
        verbose=not args.quiet
    ).initialize()

    if not args.quiet:
        utils.echo_msg(f"Processing datalist: {args.datalist}")
        if args.want_json: utils.echo_msg(f" - Generating GeoJSON Index")
        if args.want_h3:   utils.echo_msg(f" - Generating H3 Index (Res {args.h3_res})")

    ## Generate Metadata & Indices
    try:
        ## calling generate_inf() triggers the index creation side-effects
        dl.generate_inf()
    
        if not args.quiet:
            utils.echo_msg("Done.")
            
    except Exception as e:
        utils.echo_error_msg(f"Failed to process datalist: {e}")
        sys.exit(1)

        
if __name__ == "__main__":
    datalist2json_cli()
    
                
### End
