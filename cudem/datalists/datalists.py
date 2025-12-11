### datalists.py - DataLists IMproved
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
##
## datalists.py is part of CUDEM
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
### Examples:
##
### TODO:
##
### Code:

import os
import copy
import numpy as np
import pandas as pd
#from tqdm import tqdm

from osgeo import ogr

from cudem import utils
from cudem import regions
from cudem import gdalfun
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
            points = points
        elif isinstance(points, pd.DataFrame):
            points = points
        
        yield(points)

        
class Scratch(ElevationDataset):
    """Scratch Dataset

    Process a python list of valid dataset entries
    """
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.metadata['name'] = 'scratch'
        self.data_format = -3

    def generate_inf(self):
        """generate and return the infos blob of the datalist.

        datasets parsed from the datalist may have variable srs, so
        we want to transform those to whatever self.dst_srs is, if it
        exists, in order to properly fill the datalist/json srs 
        and regions...
        """

        _region = self.region
        self.region = None
        out_regions = []
        out_srs = []
        for entry in self.parse():            
            entry_minmax = entry.infos.minmax
            ## entry has an srs and dst_srs is set,
            ## so lets transform the region to suit
            if entry.src_srs is not None:
                out_srs.append(entry.src_srs)
                if self.dst_srs is not None:
                    entry_region = regions.Region().from_list(entry_minmax)
                    if entry_region.valid_p():
                        entry_region.src_srs = entry.src_srs
                        entry_region.warp(self.dst_srs)
                        entry_minmax = entry_region.export_as_list(
                            include_z=True
                        )
                        
            entry_region = regions.Region().from_list(entry_minmax)
            if entry_region.valid_p():
                out_regions.append(entry_region)
                self.infos.numpts += entry.infos.numpts

        ## merge all the gathered regions
        region_count = 0
        out_region = None
        for this_region in out_regions:
            if this_region.valid_p():
                if region_count == 0:
                    out_region = this_region
                    region_count += 1
                else:
                    out_region = regions.regions_merge(
                        out_region, this_region
                    )
                    
        if out_region is not None:
            self.infos.minmax = out_region.export_as_list(include_z=True)
            self.infos.wkt = out_region.export_as_wkt()

        ## set the epsg for the datalist
        if self.infos.src_srs is None:
            if self.src_srs is not None:
                self.infos.src_srs = self.src_srs
            elif out_srs:
                if all(x == out_srs[0] for x in out_srs):
                    self.infos.src_srs = out_srs[0]

                self.src_srs = self.infos.src_srs
                
        self.region = _region
        return(self.infos)

    
    def parse(self):
        from .dlim import DatasetFactory
        
        if isinstance(self.fn, list):
            for this_ds in self.fn:
                if this_ds is not None and this_ds.valid_p(
                        fmts=DatasetFactory._modules[this_ds.data_format]['fmts']
                ):
                    this_ds.initialize()
                    #self.data_entries.append(this_ds)
                    #yield(this_ds)
                    for ds in this_ds.parse():
                    # fill self.data_entries with each dataset for
                    # use outside the yield.
                        self.data_entries.append(ds)
                        yield(ds)

                        
class Datalist(ElevationDataset):
    """representing a datalist parser
    
    A datalist is an extended MB-System style datalist.
    
    Each datalist consists of datalist-entries, where a datalist 
    entry has the following columns:
    `path format weight uncertainty title source date data_type 
    resolution hdatum vdatum url`

    the datalist can contain datalist-entries to other datalist files, 
    distributed across a file-system.

    see `cudem.dlim.datasets` for superclass ElevationDataset
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
        """initialize the datalist geojson vector.
        
        this vector is used to quickly parse data from distributed files
        across a file-system. Each data file will have it's own entry in
        the geojson datalist vector, containing gathered metadata based on the
        source datalist.
        """

        #self.set_transform()
        self.dst_layer = '{}'.format(self.fn)
        self.dst_vector = self.dst_layer + '.json'

        utils.remove_glob('{}.json'.format(self.dst_layer))
        #if self.src_srs is not None:
        #    gdalfun.osr_prj_file('{}.prj'.format(self.dst_layer), self.src_srs)
            
        self.ds = ogr.GetDriverByName('GeoJSON').CreateDataSource(
            self.dst_vector
        )
        srs = srsfun.osr_srs(self.src_srs)
        
        if self.ds is not None: 
            self.layer = self.ds.CreateLayer(
                '{}'.format(self.dst_layer), srs, ogr.wkbMultiPolygon
            )
            [self.layer.CreateField(
                ogr.FieldDefn('{}'.format(f), ogr.OFTString)
            ) for f in self._datalist_json_cols]
            
            [self.layer.SetFeature(feature) for feature in self.layer]
            return(0)
        else:
            self.layer = None
            return(-1)

        
    def _create_entry_feature(self, entry, entry_region):
        """create a datalist entry feature and insert it into the
        datalist-vector geojson

        -----------
        Parameters:
        
        entry - the datalist entry object
        entry_region - the region of the datalist entry
        """

        #utils.echo_msg(entry)
        #utils.echo_msg(entry.params['kwargs'])
        entry_path = os.path.abspath(entry.fn) if not entry.remote else entry.fn
        entry_fields = [entry_path,
                        entry.data_format,
                        entry.weight,
                        entry.uncertainty,
                        entry.metadata['name'],
                        entry.metadata['title'],
                        entry.metadata['source'],
                        entry.metadata['date'],
                        entry.metadata['data_type'],
                        entry.metadata['resolution'],
                        entry.metadata['hdatum'],
                        entry.metadata['vdatum'],
                        entry.metadata['url'],
                        utils.dict2args(entry.params['mod_args'])]
        dst_defn = self.layer.GetLayerDefn()
        entry_geom = ogr.CreateGeometryFromWkt(entry_region.export_as_wkt())
        out_feat = ogr.Feature(dst_defn)
        out_feat.SetGeometry(entry_geom)
        for i, f in enumerate(self._datalist_json_cols):
            out_feat.SetField(f, entry_fields[i])
            
        self.layer.CreateFeature(out_feat)

        
    def generate_inf(self):
        """generate and return the infos blob of the datalist.

        datasets parsed from the datalist may have variable srs, so
        we want to transform those to whatever self.dst_srs is, if it
        exists, in order to properly fill the datalist/json srs and 
        regions...
        """

        self.region = None
        self.infos.file_hash = self.infos.generate_hash()
        _region = self.region
        out_regions = []
        out_srs = []

        ## attempt to generate a datalist-vector geojson and
        ## if successful, fill it wil the datalist entries, using `parse`
        if self._init_datalist_vector() == 0:
            for entry in self.parse():
                entry_minmax = entry.infos.minmax
                if entry.mask is not None: ## add all duplicate params???
                    entry.params['mod_args']['mask'] = entry.mask['mask']
                    for key in entry.mask.keys():
                        entry.params['mod_args'][key] = entry.mask[key]
                        #entry.params['mod_args'] = {'mask': entry.mask}

                ## entry has an srs and dst_srs is set, so lets transform
                ## the region to suit
                if entry.src_srs is not None:
                    out_srs.append(entry.src_srs)
                    if self.dst_srs is not None:
                        entry_region = regions.Region(
                            src_srs=entry.src_srs
                        ).from_list(entry_minmax).warp(self.dst_srs)
                        entry_minmax = entry_region.export_as_list(include_z=True)

                ## create the feature for the geojson
                entry_region = regions.Region().from_list(entry_minmax)
                if entry_region.valid_p():
                    out_regions.append(entry_region)
                    self._create_entry_feature(entry, entry_region)
                    self.infos.numpts += entry.infos.numpts

            self.ds = self.layer = None # close the geojson ogr dataset

        else:
            utils.echo_warning_msg(
                'could not initialize datalist vector'
            )
            return(self.infos)

        ## merge all the gathered regions
        region_count = 0
        out_region = None
        for this_region in out_regions:
            if this_region.valid_p():
                if region_count == 0:
                    out_region = this_region
                    region_count += 1
                else:
                    out_region = regions.regions_merge(out_region, this_region)
                    
        if out_region is not None:
            self.infos.minmax = out_region.export_as_list(include_z=True)
            self.infos.wkt = out_region.export_as_wkt()

        ## set the epsg for the datalist
        if self.infos.src_srs is None:
            if self.src_srs is not None:
                self.infos.src_srs = self.src_srs
            elif out_srs:
                if all(x == out_srs[0] for x in out_srs):
                    self.infos.src_srs = out_srs[0]

                self.src_srs = self.infos.src_srs

        self.region = _region
        return(self.infos)

    
    def parse_json(self):
        """parse the datalist using the datalist-vector geojson.

        Quickly find data in a given region using the datalist-vector. 
        The datalist-vector must have been previously generated using `parse`. 
        If the datlist-vector is not found will fall-back to `parse` and 
        generate a new datalist-vector geojson.

        -------
        Yields:
        dataset object of each dataset found in the datalist-vector
        """

        from .dlim import DatasetFactory
        
        ## check for the datalist-vector geojson
        status = 0
        count = 0
        
        ## user input to re-gerenate json...?
        if os.path.exists('{}.json'.format(self.fn)):
            driver = ogr.GetDriverByName('GeoJSON')
            dl_ds = driver.Open('{}.json'.format(self.fn))
            if dl_ds is None:
                # utils.echo_warning_msg(
                #     f'could not open {self.fn}.json'
                # )
                status = -1
        else:
            status = -1

        ## parse the datalist-vector geojson and yield the results
        if status != -1:
            dl_layer = dl_ds.GetLayer()
            ldefn = dl_layer.GetLayerDefn()
            _boundsGeom = None
            if self.region is not None:
                _boundsGeom = self.region.export_as_geom() \
                    if self.transform['transformer'] is None \
                       else self.transform['trans_region'].export_as_geom()

            dl_layer.SetSpatialFilter(_boundsGeom)
            count = len(dl_layer)
            with utils.ccp(
                    total=len(dl_layer),
                    desc=f'parsing {count} datasets from datalist json {self.fn}.json @ {self.weight}',
                    leave=self.verbose
            ) as pbar:
                for l,feat in enumerate(dl_layer):
                    pbar.update()
                    ## filter by input source region extras (weight/uncertainty)
                    if self.region is not None:
                        w_region = self.region.w_region()
                        if w_region[0] is not None:
                            if float(feat.GetField('weight')) < w_region[0]:
                                continue

                        if w_region[1] is not None:
                            if float(feat.GetField('weight')) > w_region[1]:
                                continue

                        u_region = self.region.u_region()
                        if u_region[0] is not None:
                            if float(feat.GetField('uncertainty')) < u_region[0]:
                                continue

                        if u_region[1] is not None:
                            if float(feat.GetField('uncertainty')) > u_region[1]:
                                continue

                    ## extract the module arguments from the datalist-vector
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

                    ## update existing metadata
                    md = copy.deepcopy(self.metadata)
                    for key in self.metadata.keys():
                        md[key] = feat.GetField(key)

                    ## generate the dataset object to yield
                    data_mod = '"{}" {}{} {} {}'.format(
                        feat.GetField('path'),
                        feat.GetField('format'),
                        ':{}'.format(ds_args) if ds_args is not None else '',
                        feat.GetField('weight'),
                        feat.GetField('uncertainty')
                    )
                    data_set = DatasetFactory(
                        **self._set_params(mod=data_mod, metadata=md)#, **data_set_args)
                    )._acquire_module()
                    if data_set is not None and data_set.valid_p(
                            fmts=DatasetFactory._modules[data_set.data_format]['fmts']
                    ):
                        data_set.initialize()
                        #utils.echo_msg(data_set.params)
                        ## fill self.data_entries with each dataset for use outside the yield.
                        for ds in data_set.parse(): 
                            self.data_entries.append(ds)
                            yield(ds)

            dl_ds = dl_layer = None
                
        else:
            ## failed to find/open the datalist-vector geojson, so run `parse` instead and
            ## generate one for future use...
            # utils.echo_warning_msg(
            #     'could not load datalist-vector json {}.json,
            # falling back to parse, generate a json file for the datalist using
            # `dlim -i`'.format(self.fn)
            # )
            for ds in self.parse_no_json():
                yield(ds)

                
    def parse(self):
        """parse the datalist file.

        -------
        Yields:
        dataset object of each dataset found in the datalist        
        """

        from .dlim import DatasetFactory
        
        status = 0
        if os.path.exists(self.fn):
            with open(self.fn, 'r') as f:
                count = sum(1 for _ in f)

            with open(self.fn, 'r') as op:
                with utils.ccp(
                        total=count,
                        desc=f'parsing datalist {self.fn}...',
                        leave=False
                ) as pbar:
                    for l, this_line in enumerate(op):
                        pbar.update()
                        ## parse the datalist entry line
                        if this_line[0] != '#' \
                           and this_line[0] != '\n' \
                               and this_line[0].rstrip() != '':
                            md = copy.deepcopy(self.metadata)
                            md['name'] = utils.fn_basename2(
                                os.path.basename(self.fn)
                            )

                            ## generate the dataset object to yield
                            # ds_kwargs = self.kwargs.copy()
                            # ds_kwargs['mod'] = this_line
                            # ds_kwargs['metadata'] = md
                            # ds_kwargs['src_srs'] = self.src_srs
                            # ds_kwargs['parent'] = self
                            # ds_kwargs['fn'] = None
                            ds_params = self._set_params(
                                #**ds_kwargs
                                mod=this_line,
                                metadata=md,
                                src_srs=self.src_srs,
                                parent=self,
                                fn=None,
                            )
                            #utils.echo_msg_bold(ds_params)
                            data_set = DatasetFactory(**ds_params)._acquire_module()
                            if data_set is not None and data_set.valid_p(
                                    fmts=DatasetFactory._modules[
                                        data_set.data_format
                                    ]['fmts']
                            ):
                                data_set.initialize()
                                ## filter with input source region, if necessary
                                ## check input source region against the dataset
                                ## region found in its inf file.
                                if self.region is not None \
                                   and self.region.valid_p(check_xy=True):
                                    # inf_region = regions.Region().from_list(
                                    #     data_set.infos.minmax
                                    # )
                                    # if inf_region.valid_p():
                                    #     inf_region.wmin = data_set.weight
                                    #     inf_region.wmax = data_set.weight
                                    #     inf_region.umin = data_set.uncertainty
                                    #     inf_region.umax = data_set.uncertainty

                                    if not regions.regions_intersect_p(
                                            self.inf_region,
                                            self.region \
                                            if data_set.transform['trans_region'] is None \
                                            else data_set.transform['trans_region']
                                    ):
                                        continue
                                        
                                ## fill self.data_entries with each dataset for use
                                ## outside the yield and yield the dataset object.
                                #yield(data_set)
                                for ds in data_set.parse(): 
                                    self.data_entries.append(ds)
                                    yield(ds)

        ## self.fn is not a file-name, so check if self.data_entries not empty
        ## and return the dataset objects found there.
        elif len(self.data_entries) > 0:
            for data_set in self.data_entries:
                for ds in data_set.parse():
                    yield(ds)
        else:
            if self.verbose:
                utils.echo_warning_msg(
                    f'could not open datalist/entry {self.fn}'
                )


### End
