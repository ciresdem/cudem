### ziplists.py 
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
##
## ziplists.py is part of CUDEM
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
import zipfile
        
from cudem import utils
from cudem import regions
from cudem.datalists import dlim

class ZIPlist(dlim.ElevationDataset):
    """Zip file parser.

    Parse supported datasets from a zipfile.
    """
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        
    def generate_inf(self):
        """generate and return the infos blob of the datalist.

        datasets parsed from the datalist may have variable srs, so
        we want to transform those to whatever self.dst_srs is, if it
        exists, in order to properly fill the datalist/json srs and regions...
        """

        _region = self.region
        self.region = None
        out_regions = []
        out_srs = []
        self.infos.file_hash = self.infos.generate_hash()
        for entry in self.parse():            
            entry_minmax = entry.infos.minmax

            ## entry has an srs and dst_srs is set, so lets transform
            ## the region to suit
            if entry.src_srs is not None:
                out_srs.append(entry.src_srs)
                if self.dst_srs is not None:
                    entry_region = regions.Region().from_list(entry_minmax)
                    entry_region.src_srs = entry.src_srs
                    entry_region.warp(self.dst_srs)
                    entry_minmax = entry_region.export_as_list(include_z=True)

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

    
    def parse(self):
        from cudem.datalists.dlim_factory import DatasetFactory
        
        exts = [DatasetFactory()._modules[x]['fmts'] \
                for x in DatasetFactory()._modules.keys()]
        exts = [x for y in exts for x in y]
        datalist = []
        if self.fn.split('.')[-1].lower() == 'zip':
            try:
                with zipfile.ZipFile(self.fn) as z:
                    zfs = z.namelist()
                    for ext in exts:
                        for zf in zfs:
                            if ext == zf.split('.')[-1]:
                                datalist.append(os.path.basename(zf))
            except Exception as e:
                utils.echo_error_msg(
                    f'could not unzip {self.fn}, {e}'
                )
                            
        for this_data in datalist:
            this_line = utils.p_f_unzip(
                self.fn,
                fns=[this_data],
                outdir=os.path.normpath(os.path.dirname(self.fn)),
                tmp_fn=True
            )[0]
            data_set = DatasetFactory(
                **self._set_params(
                    mod=os.path.basename(this_line),
                    data_format=None,
                    src_srs=self.src_srs,
                    parent=self
                )
            )._acquire_module()
            if data_set is not None and data_set.valid_p(
                    fmts=DatasetFactory._modules[data_set.data_format]['fmts']
            ):
                data_set.initialize()
                if self.region is not None \
                   and self.region.valid_p(check_xy=True):
                    inf_region = self.inf_region.copy()
                    # try:
                    #     inf_region = regions.Region().from_string(
                    #         data_set.infos.wkt
                    #     )
                    # except:
                    #     try:
                    #         inf_region = regions.Region().from_list(
                    #             data_set.infos.minmax
                    #         )
                    #     except:
                    #         inf_region = self.region.copy()
                            
                    inf_region.wmin = data_set.weight
                    inf_region.wmax = data_set.weight
                    inf_region.umin = data_set.uncertainty
                    inf_region.umax = data_set.uncertainty
                    if inf_region.valid_p(check_xy=True):
                        if regions.regions_intersect_p(
                                inf_region,
                                self.region \
                                if data_set.transform['transformer'] is None \
                                else data_set.transform['trans_region']
                        ):
                            for ds in data_set.parse():
                                self.data_entries.append(ds)
                                yield(ds)
                    else:
                        if self.verbose:
                            utils.echo_warning_msg(
                                f'invalid inf file: {data_set.fn}.inf, skipping'
                            )
                            utils.remove_glob(f'{data_set.fn}*')
                            
                else:
                    for ds in data_set.parse():
                        self.data_entries.append(ds)
                        yield(ds)
                        
            utils.remove_glob(f'{data_set.fn}*')


### End
