### mbsfile.py 
##
## Copyright (c) 2010 - 2026 Regents of the University of Colorado
##
## mbsfile.py is part of CUDEM
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
## MB-System Multibeam Parser
##
### Code:

import os
import math
import numpy as np
import pandas as pd

from osgeo import gdal, ogr

from cudem import utils
from cudem import regions
from cudem import gdalfun
from cudem.grits import grits
from cudem.datalists.dlim import ElevationDataset

class MBSParser(ElevationDataset):
    """Providing an mbsystem parser.
    Process MB-System supported multibeam data files.
    """

    def __init__(self,
                 mb_fmt=None,
                 mb_exclude='A',
                 want_mbgrid=False,
                 want_binned=False,
                 min_year=None,
                 auto_weight=True,
                 auto_uncertainty=True,
                 **kwargs):
        
        super().__init__(**kwargs)
        
        self.mb_fmt = mb_fmt
        self.mb_exclude = mb_exclude
        self.want_mbgrid = want_mbgrid
        self.want_binned = want_binned
        self.min_year = min_year
        self.auto_weight = auto_weight
        self.auto_uncertainty = auto_uncertainty
        
        if self.src_srs is None:
            self.src_srs = 'epsg:4326'

            
    def generate_inf(self, make_grid=True, make_block_mean=False, block_inc=None):
        """Generate metadata (INF) for the MBS data.
        
        Attempts to parse existing MB-System .inf files first for speed.
        If grids are requested or parsing fails, performs a full scan via parent class.
        """
        
        ## Try to parse existing MB-System .inf file
        inf_fn = f"{self.fn}.inf" #if not self.fn.endswith('.inf') else self.fn
        parsed_ok = False
        
        if os.path.exists(inf_fn) and not (make_grid or make_block_mean):
            try:
                self._parse_mbs_inf_file(inf_fn)
                parsed_ok = True
            except Exception:
                parsed_ok = False

        ## Fallback to full scan if needed
        if not parsed_ok or make_grid or make_block_mean:
            return super().generate_inf(
                make_grid=make_grid, 
                make_block_mean=make_block_mean, 
                block_inc=block_inc
            )
            
        return self.infos

    
    def _parse_mbs_inf_file(self, inf_fn):
        """Internal helper to parse native MB-System .inf files."""
        
        self.infos.minmax = [0, 0, 0, 0, 0, 0]
        dims = []
        cm_array = None
        this_row = 0
        xinc, yinc = 0, 0

        with open(inf_fn) as iob:
            for line in iob:
                parts = line.split()
                if not parts: continue

                if parts[0] == 'Swath' and parts[2] == 'File:':
                    self.infos.name = parts[3]

                if ' '.join(parts[:-1]) == 'MBIO Data Format ID:':
                    self.mb_fmt = parts[-1]
                        
                if parts[0] == 'Number' and parts[2] == 'Records:':
                    self.infos.numpts = utils.int_or(parts[3])
                        
                if parts[0] == 'Minimum':
                    val = utils.float_or(parts[2])
                    val2 = utils.float_or(parts[5])
                    
                    if parts[1] == 'Longitude:':
                        self.infos.minmax[0] = val
                        self.infos.minmax[1] = val2
                    elif parts[1] == 'Latitude:':
                        self.infos.minmax[2] = val
                        self.infos.minmax[3] = val2
                    elif parts[1] == 'Depth:':
                        # Convert depth to elevation (negate)
                        self.infos.minmax[4] = val2 * -1
                        self.infos.minmax[5] = val * -1
                        
                if parts[0] == 'CM' and parts[1] == 'dimensions:':
                    dims = [utils.int_or(parts[2]), utils.int_or(parts[3])]
                    cm_array = np.zeros((dims[0], dims[1]))
                        
                if parts[0] == 'CM:' and cm_array is not None:
                    for j in range(dims[0]):
                        if j+1 < len(parts):
                            cm_array[this_row][j] = utils.int_or(parts[j+1])
                    this_row += 1

        ## Generate Region & WKT
        mbs_region = regions.Region().from_list(self.infos.minmax)
        
        if dims:
            xinc = (mbs_region.xmax - mbs_region.xmin) / dims[0]
            yinc = (mbs_region.ymin - mbs_region.ymax) / dims[1]
        
        if abs(xinc) > 0 and abs(yinc) > 0 and cm_array is not None:
            ## Generate Polygon from CM array
            _, _, dst_gt = mbs_region.geo_transform(x_inc=xinc, y_inc=yinc)
            
            ## Create in-memory raster for polygonization
            mem_driver = gdal.GetDriverByName('MEM')
            ds = mem_driver.Create('', dims[0], dims[1], 1, gdal.GDT_Float32)
            ds.SetGeoTransform(dst_gt)
            ds.SetProjection(gdalfun.osr_wkt(self.src_srs))
            ds.GetRasterBand(1).WriteArray(cm_array)
            
            ## Polygonize
            tmp_ds = ogr.GetDriverByName('Memory').CreateDataSource('tmp_poly')
            tmp_layer = tmp_ds.CreateLayer('tmp_poly', None, ogr.wkbMultiPolygon)
            tmp_layer.CreateField(ogr.FieldDefn('DN', ogr.OFTInteger))
            
            gdal.Polygonize(ds.GetRasterBand(1), None, tmp_layer, 0)
            
            multi = ogr.Geometry(ogr.wkbMultiPolygon)
            for feat in tmp_layer:
                ## Assuming DN != 0 is valid data
                if feat.GetField('DN') != 0:
                    feat.geometry().CloseRings()
                    multi.AddGeometry(feat.geometry())
            
            ## Union/Merge if needed, currently just dumping all valid parts
            self.infos.wkt = multi.ExportToWkt()
            ds = None
        else:
            self.infos.wkt = mbs_region.export_as_wkt()

        if self.infos.src_srs is None:
            self.infos.src_srs = self.src_srs

            
    def _get_mbs_meta(self, src_inf):
        """Extract metadata from mbsystem inf file."""
        
        meta = {'format': None, 'date': None, 'perc_good': None}
        
        if not os.path.exists(src_inf): return meta
        
        try:
            with open(src_inf, errors='ignore') as f:
                for line in f:
                    parts = line.split()
                    if not parts: continue
                    
                    if parts[0].strip() == 'MBIO':
                        meta['format'] = parts[4]
                    elif parts[0] == 'Time:':
                        meta['date'] = parts[3]
                    
                    if ':' in line:
                        p = line.split(':')
                        if p[0].strip() == 'Number of Good Beams':
                            meta['perc_good'] = p[1].split()[-1].split('%')[0]
        except Exception:
            pass
            
        return meta

    
    def yield_mbgrid_ds(self):
        """Process data through mbgrid -> GDAL -> Points."""

        from cudem.datalists.dlim import DatasetFactory
        
        mb_fn = self.fn
        tmp_dl = '_mb_grid_tmp.datalist'
        
        with open(tmp_dl, 'w') as f:
            fmt_str = f" {self.mb_fmt}" if self.mb_fmt else ""
            f.write(f"{self.fn}{fmt_str} {self.weight}\n")

        base_name = os.path.splitext(os.path.basename(mb_fn))[0]
        out_name = f"{base_name}_mbgrid"
        
        try:
            ## Run mbgrid
            cmd = (f'mbgrid -I{tmp_dl} {self.data_region.format("gmt")} '
                   f'-E{self.x_inc}/{self.y_inc}/degrees! -O{out_name} '
                   f'-A2 -F1 -C2/1 -S0 -T35 -W.1')
            utils.run_cmd(cmd, verbose=self.verbose)

            ## Convert to GTiff
            if os.path.exists(f"{out_name}.grd"):
                gdalfun.gdal2gdal(f"{out_name}.grd", dst_fmt='GTiff')
            
            ## Filter Outliers (Optional)
            tif_fn = f"{out_name}.tif"
            if os.path.exists(tif_fn):
                grits_filter = grits.GritsFactory(
                    mod='outliers:multipass=2',
                    src_dem=tif_fn,
                    cache_dir=self.cache_dir,
                )._acquire_module()
                
                if grits_filter:
                    grits_filter()
                    if os.path.exists(grits_filter.dst_dem):
                        os.replace(grits_filter.dst_dem, tif_fn)

                ## Use Factory to yield points from the resulting grid
                #mbs_ds = DatasetFactory(mod=tif_fn, data_format=200)._acquire_module()
                mbs_ds = DatasetFactory(
                    **self._set_params(
                        mod=tif_fn,
                        data_format=200
                    )
                )._acquire_module()
                mbs_ds.initialize()
                for gdal_ds in mbs_ds.parse():
                    for pts in gdal_ds.transform_and_yield_points():
                        yield pts

        except Exception as e:
            if self.verbose:
                utils.echo_error_msg(f"mbgrid failed: {e}")
        finally:
            utils.remove_glob(tmp_dl, f"{out_name}.cmd", f"{out_name}.mb-1", 
                              f"{out_name}.grd*", f"{out_name}.tif*")


    #def read_mblist_ds(self):
        
    def yield_mblist_ds(self):
        """Use mblist to stream points."""
        
        ## Determine format/metadata
        src_inf = f"{self.fn}.inf"
        meta = self._get_mbs_meta(src_inf)
        
        mb_format = meta['format']
        ## .fbt files override format detection
        if self.fn.endswith('.fbt'): mb_format = None
        
        ## Calculate Weight
        this_weight = 1.0
        if self.auto_weight:
            if meta['date']:
                ## Decaying weight based on age
                this_weight = min(0.99, max(0.01, 1 - ((2024 - int(meta['date'])) / (2024 - 1980))))
            
            if self.weight is not None:
                self.weight *= this_weight

        ## Build Command
        mb_region = None
        if self.region is not None:
            mb_region = self.region.copy()
            mb_region.buffer(pct=5)
            
        #mb_region = self.region.copy().buffer(pct=25) if self.region else None
        region_arg = f" {mb_region.format('gmt')}" if mb_region else ""
        fmt_arg = f" -F{mb_format}" if mb_format else ""
        
        cmd = f'mblist -M{self.mb_exclude} -OXYZDS{region_arg} -I{self.fn}{fmt_arg}'
        
        mb_points = []
        BATCH_SIZE = 1000000

        for line in utils.yield_cmd(cmd, verbose=False):
            try:
                ## X Y Z CrossTrack Speed
                vals = [float(x) for x in line.strip().split('\t')]
                if len(vals) < 5: continue
                
                x, y, z, xtrack, speed = vals
                if speed == 0: continue

                w = 1.0
                u = 0.0

                if self.auto_weight:
                    ## Uncertainty Calculation
                    ## U_depth = 0.51 * (0.25 + 0.02 * depth)
                    u_depth = (0.25 + (0.02 * abs(z))) * 0.51
                    u_xtrack = 0.005 * abs(xtrack)
                    
                    ## Speed penalty
                    tmp_speed = min(14, abs(speed - 14))
                    u_speed = tmp_speed * 0.51
                    
                    u = math.sqrt(u_depth**2 + u_xtrack**2 + u_speed**2)
                    w = 1.0 / u if u > 0 else 1.0
                
                mb_points.append([x, y, z, w, u])

                if len(mb_points) >= BATCH_SIZE:
                    yield np.rec.fromrecords(mb_points, names='x, y, z, w, u')
                    mb_points = []

            except (ValueError, IndexError):
                continue

        if mb_points:
            yield np.rec.fromrecords(mb_points, names='x, y, z, w, u')

            
    def yield_points(self):        
        # If region/incs not set, force mblist
        if self.region is None or self.x_inc is None or self.y_inc is None:
            self.want_mbgrid = False

        if self.want_mbgrid:
            for pts in self.yield_mbgrid_ds():
                yield pts
        else:
            for pts in self.yield_mblist_ds():
                yield pts

### End
