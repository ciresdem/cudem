### waffles.py
##
## Copyright (c) 2010 - 2021 CIRES Coastal DEM Team
##
## waffles.py is part of CUDEM
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
### Code:

import sys
import os
import math
import json
import numpy as np
from osgeo import gdal
from osgeo import ogr
from osgeo import osr
import cudem
from cudem import dlim
from cudem import datasets
from cudem import regions
from cudem import utils
from cudem import xyzfun
from cudem import demfun
from cudem import metadata
from cudem import vdatumfun

class Waffle:
    """Representing a WAFFLES DEM/MODULE.
    Specific Gridding modules are sub-classes of this class.
    See WaffleFactory for module specifications and generation.

    Procedures:
      yield_xyz() - yield the xyz data from self.data
      dump_xyz() - dump the xyz data from self.data to port
      run() - run the WAFFLES module (function set via module sub-class)
      generate() - run and process the WAFFLES module

    """
    
    def __init__(
            self,
            data=[],
            src_region=None,
            inc=None,
            xinc=None,
            yinc=None,
            xsize=None,
            ysize=None,
            name='waffles_dem',
            node='pixel',
            fmt='GTiff',
            extend=0,
            extend_proc=20,
            weights=None,
            fltr=[],
            sample=None,
            xsample=None,
            ysample=None,
            clip=None,
            chunk=None,
            epsg=4326,
            verbose=False,
            archive=False,
            mask=False,
            spat=False,
            clobber=True
    ):

        self.data = data
        self.region = src_region
        self.inc = inc
        self.xinc = xinc
        self.yinc = yinc
        self.sample = sample
        self.xsample = xsample
        self.ysample = ysample
        #self.size = size
        self.name = name
        self.node = node
        self.fmt = fmt
        self.extend = extend
        self.extend_proc = extend_proc
        self.weights = weights
        self.fltr = fltr
        self.clip = clip
        self.chunk = chunk
        self.epsg = utils.int_or(epsg)
        self.mod = None
        self.mod_args = {}
        self.archive = archive
        self.mask = mask
        #self.data_mask = None
        self.clobber = clobber
        self.verbose = verbose
        self.gc = utils.config_check()
        self.spat = spat
        self.block_t = None
        self.ogr_ds = None
        
        if self.node == 'grid':
            self.region = self.region.buffer(x_bv=self.xinc*.5, y_bv=self.yinc*.5)
            #self._init_regions()
        self.p_region = self._proc_region()
        self.d_region = self._dist_region()
        self.c_region = self._coast_region()
        
        self.data_ = data
        self._init_data()
        # self.data = [dlim.DatasetFactory(
        #     fn=" ".join(['-' if x == "" else x for x in dl.split(",")]),
        #     src_region=self.p_region,
        #     verbose=self.verbose,
        #     weight=self.weights,
        #     warp=self.epsg
        # ).acquire_dataset() for dl in self.data]

        # self.data = [d for d in self.data if d is not None]
        
        self.fn = '{}.tif'.format(self.name)
        self.mask_fn = '{}_m.tif'.format(self.name)
        self.waffled = False

    def _init_regions(self):    
        self.p_region = self._proc_region()
        self.d_region = self._dist_region()
        self.c_region = self._coast_region()
        
    def _init_data(self):

        self.data = [dlim.DatasetFactory(
            fn=" ".join(['-' if x == "" else x for x in dl.split(",")]),
            src_region=self.p_region, verbose=self.verbose,
            epsg=self.epsg, weight=self.weights).acquire_dataset() for dl in self.data_]

        self.data = [d for d in self.data if d is not None]
        
    def _set_config(self):
        """export the waffles config info as a dictionary"""
        
        self._config = {
            self.mod: {
                'data': self.data_,
                'src_region': self.region.export_as_list(
                    include_z=True, include_w=True
                ) if self.region is not None else None,
                'inc': self.inc,
                'xinc': self.xinc,
                'yinc': self.yinc,
                'xsample': self.xsample,
                'ysample': self.ysample,
                'name': self.name,
                'node': self.node,
                'fmt': self.fmt,
                'extend': self.extend,
                'extend_proc': self.extend_proc,
                'weights': self.weights,
                'fltr': self.fltr,
                'sample': self.sample,
                'xsample': self.xsample,
                'ysample': self.ysample,
                'clip': self.clip,
                'chunk': self.chunk,
                'epsg': self.epsg,
                'verbose': self.verbose,
                'archive': self.archive,
                'mask': self.mask,
                'spat': self.spat,
                'clobber': self.clobber,
            }
        }

        for key in self.mod_args.keys():
            if key not in self._config[self.mod].keys():
                self._config[self.mod][key] = self.mod_args[key]

        return(self._config)

    def _coast_region(self):
        """processing region (extended by self.extend and self.extend_proc."""
        
        cr = regions.Region().from_region(self.region)
        #return(cr.buffer((self.inc*self.extend_proc) + (self.inc*self.extend) + (self.inc*10)))
        return(cr.buffer(x_bv=(self.xinc*self.extend)+ (self.xinc*self.extend) + (self.xinc*10), y_bv=(self.yinc*self.extend) + (self.yinc*self.extend) + (self.yinc*10)))
    
    def _proc_region(self):
        """processing region (extended by self.extend and self.extend_proc."""
        
        pr = regions.Region().from_region(self.region)
        #return(pr.buffer((self.inc*self.extend_proc) + (self.inc*self.extend)))
        return(pr.buffer(x_bv=(self.xinc*self.extend), y_bv=(self.yinc*self.extend)))
    
    def _dist_region(self):
        """distribution region (extended by self.extend)."""
        
        dr = regions.Region().from_region(self.region)
        return(dr.buffer(x_bv=(self.xinc*self.extend), y_bv=(self.yinc*self.extend)))

    def _xyz_block_t(self, src_xyz):
        """block the src_xyz data for fast lookup"""

        xcount, ycount, dst_gt = self.p_region.geo_transform(x_inc=self.xinc, y_inc=self.yinc)
        self.block_t = np.empty((ycount, xcount), dtype=object)
        for y in range(0, ycount):
            for x in range(0, xcount):
                self.block_t[y,x] = []

        if self.verbose:
            utils.echo_msg(
                'blocking data into {} buckets'.format(ycount*xcount)
            )
            
        for this_xyz in src_xyz:
            if regions.xyz_in_region_p(this_xyz, self.p_region):
                xpos, ypos = utils._geo2pixel(
                    this_xyz.x, this_xyz.y, dst_gt
                )
                if xpos < xcount and ypos < ycount:
                    self.block_t[ypos,xpos].append(this_xyz.copy())

    def _xyz_ds(self, src_xyz):
        """Make a point vector OGR DataSet Object from src_xyz"""

        dst_ogr = '{}'.format(self.name)
        self.ogr_ds = gdal.GetDriverByName('Memory').Create(
            '', 0, 0, 0, gdal.GDT_Unknown
        )
        layer = self.ogr_ds.CreateLayer(dst_ogr, geom_type = ogr.wkbPoint25D)
        fd = ogr.FieldDefn('long', ogr.OFTReal)
        fd.SetWidth(10)
        fd.SetPrecision(8)
        layer.CreateField(fd)
        fd = ogr.FieldDefn('lat', ogr.OFTReal)
        fd.SetWidth(10)
        fd.SetPrecision(8)
        layer.CreateField(fd)
        fd = ogr.FieldDefn('elev', ogr.OFTReal)
        fd.SetWidth(12)
        fd.SetPrecision(12)
        layer.CreateField(fd)
        
        if self.weights:
            fd = ogr.FieldDefn('weight', ogr.OFTReal)
            fd.SetWidth(6)
            fd.SetPrecision(6)
            layer.CreateField(fd)
            
        f = ogr.Feature(feature_def = layer.GetLayerDefn())
        
        for this_xyz in self.yield_xyz():#src_xyz:
            #print(this_xyz.x, this_xyz.y, this_xyz.z)
            f.SetField(0, this_xyz.x)
            f.SetField(1, this_xyz.y)
            f.SetField(2, float(this_xyz.z))
            if self.weights:
                f.SetField(3, this_xyz.w)
                
            wkt = this_xyz.export_as_wkt(include_z=True)
            g = ogr.CreateGeometryFromWkt(wkt)
            f.SetGeometryDirectly(g)
            layer.CreateFeature(f)
            
        return(self.ogr_ds)

    def _xyz_block(self, src_xyz):
        """block the src_xyz data to the mean block value

        Yields:
          list: xyz data for each block with data
        """

        xcount, ycount, dst_gt = self.p_region.geo_transform(
            x_inc=self.xinc, y_inc=self.yinc
        )
        gdt = gdal.GDT_Float32
        sum_array = np.zeros((ycount, xcount))
        count_array = np.zeros((ycount, xcount))
        if self.weights:
            weight_array = np.zeros((ycount, xcount))
            
        if self.verbose:
            utils.echo_msg(
                'blocking data to {}/{} grid'.format(ycount, xcount)
            )
            
        for this_xyz in src_xyz:
            if regions.xyz_in_region_p(this_xyz, self.p_region):
                if self.weights:
                    this_z = this_xyz.z * this_xyz.w
                else:
                    this_z = this_xyz.z
                
                xpos, ypos = utils._geo2pixel(
                    this_xyz.x, this_xyz.y, dst_gt
                )
                try:
                    sum_array[ypos, xpos] += this_z
                    count_array[ypos, xpos] += 1
                    if self.weights:
                        weight_array[ypos, xpos] += this_xyz.w
                        
                except: pass

        count_array[count_array == 0] = np.nan
        if self.weights:
            weight_array[weight_array == 0] = np.nan
            out_weight_array = (weight_array/count_array)
            out_array = (sum_array/out_weight_array)/count_array
        else:
            out_array = (sum_array/count_array)
            out_weight_array = np.ones((ycount, xcount))
            
        sum_array = count_array = None
        if self.weights:
            weight_array = None
            
        for y in range(0, ycount):
            for x in range(0, xcount):
                z = out_array[y,x]
                if not np.isnan(z):
                    geo_x, geo_y = utils._pixel2geo(x, y, dst_gt)
                    out_xyz = xyzfun.XYZPoint(
                        x=geo_x, y=geo_y, z=z, w=out_weight_array[y,x]
                    )
                    yield(out_xyz)

        out_array = out_weight_array = None
                        
    def _xyz_mask(self, src_xyz, dst_gdal, dst_format='GTiff'):
        """Create a num grid mask of xyz data. The output grid
        will contain 1 where data exists and 0 where no data exists.

        yields the xyz data
        """

        xcount, ycount, dst_gt = self.region.geo_transform(x_inc=self.xinc, y_inc=self.yinc)
        ptArray = np.zeros((ycount, xcount))
        ds_config = demfun.set_infos(
            xcount, ycount, (xcount*ycount), dst_gt, utils.sr_wkt(self.epsg),
            gdal.GDT_Float32, -9999, 'GTiff')
        for this_xyz in src_xyz:
            yield(this_xyz)
            
            if regions.xyz_in_region_p(this_xyz, self.region):
                xpos, ypos = utils._geo2pixel(
                    this_xyz.x, this_xyz.y, dst_gt
                )
                try:
                    ptArray[ypos, xpos] = 1
                except:
                    pass
                
        out, status = utils.gdal_write(ptArray, dst_gdal, ds_config)    

    ## TODO: allow spat-meta and archive at same time...
    def yield_xyz(self, block=False, **kwargs):
        """yields the xyz data"""
        
        for xdl in self.data:
            if self.spat:
                xyz_yield = metadata.SpatialMetadata(
                    data=[xdl.fn],
                    src_region=self.p_region,
                    inc=self.xinc,
                    extend=self.extend,
                    epsg=self.epsg,
                    node=self.node,
                    name=self.name,
                    verbose=self.verbose,
                    make_valid=True
                ).yield_xyz()
                
            elif self.archive:
                xyz_yield = xdl.archive_xyz()
            else:
                xyz_yield = xdl.yield_xyz()
                
            if self.mask:
                xyz_yield = self._xyz_mask(xyz_yield, self.mask_fn)
                
            if block:
                xyz_yield = self._xyz_block(xyz_yield)
                
            for xyz in xyz_yield:
                yield(xyz)

        if self.spat:
            dst_layer = '{}_sm'.format(self.name)
            dst_vector = dst_layer + '.shp'
            utils.run_cmd(
                'ogr2ogr -clipsrc {} __tmp_clip.shp {} -overwrite -nlt POLYGON -skipfailures'.format(
                    self.d_region.format('ul_lr'), dst_vector
                ),
                verbose=True
            )
            utils.run_cmd(
                'ogr2ogr {} __tmp_clip.shp -overwrite'.format(dst_vector),
                verbose=True
            )
            utils.remove_glob('__tmp_clip.*')

    def dump_xyz(self, dst_port=sys.stdout, encode=False, **kwargs):
        """dump the xyz data to dst_port"""
        
        for xyz in self.yield_xyz(**kwargs):
            xyz.dump(
                include_w = True if self.weights is not None else False,
                dst_port=dst_port,
                encode=encode,
                **kwargs
            )

    def _process(self, fn=None, filter_=False):
        """process the outpt WAFFLES DEM (filter, clip, cut, set)."""
        
        if fn is None:
            fn = self.fn
        
        demfun.set_nodata(fn, nodata=-9999, convert_array=True)
        
        if filter_:
            if len(self.fltr) > 0:
                for f in self.fltr:
                    fltr_val = None
                    split_val = None
                    fltr_opts = f.split(':')
                    fltr = fltr_opts[0]
                    if len(fltr_opts) > 1:
                        fltr_val = fltr_opts[1]
                        
                    if len(fltr_opts) > 2:
                        split_val= fltr_opts[2]
                    
                    if demfun.filter_(
                            fn, '__tmp_fltr.tif', fltr=fltr, fltr_val=fltr_val, split_val=split_val,
                    ) == 0:
                        os.rename('__tmp_fltr.tif', fn)
            
        if self.xsample is not None or self.ysample is not None:
            if demfun.sample(fn, '__tmp_sample.tif', self.xsample, self.ysample, self.p_region)[1] == 0:
                os.rename('__tmp_sample.tif', fn)
            
        if self.clip is not None:

            clip_args = {}
            cp = self.clip.split(':')
            clip_args['src_ply'] = cp[0]
            clip_args = utils.args2dict(cp[1:], clip_args)
            #print(clip_args)
            
            if clip_args['src_ply'] == 'coastline':
                self.coast = WaffleFactory(
                    mod='coastline:invert=True',
                    data=[],
                    src_region=self.c_region,
                    xinc=self.xinc,
                    yinc=self.yinc,
                    name='tmp_coast',
                    node=self.node,
                    extend=self.extend+12,
                    weights=self.weights,
                    epsg=self.epsg,
                    clobber=True,
                    verbose=self.verbose,
                ).acquire().generate()
                clip_args['src_ply'] = 'tmp_coast.shp'

            if demfun.clip(fn, '__tmp_clip__.tif', **clip_args)[1] == 0:
                os.rename('__tmp_clip__.tif', '{}'.format(fn))
                
            # gdalfun.gdal_clip(this_dem, **clip_args)
            # if this_wg['mask']:
            #     if this_wg['verbose']: utils.echo_msg('clipping {}...'.format(this_dem_msk))
            #     gdalfun.gdal_clip(this_dem_msk, **clip_args)
                
        if demfun.cut(fn, self.d_region, '__tmp_cut__.tif')[1] == 0:
            os.rename('__tmp_cut__.tif', '{}'.format(fn))
            
        demfun.set_metadata(fn, node=self.node)
        demfun.set_epsg(fn, self.epsg)
        
        if self.fmt != 'GTiff':
            out_dem = utils.gdal2gdal(fn, dst_fmt=self.fmt)
            if out_dem is not None:
                utils.remove_glob(fn)
                
        return(self)
    
    def generate(self):
        """run and process the WAFFLES module"""
        
        if os.path.exists(self.fn):
            if not self.clobber:
                utils.echo_msg(
                    'DEM {} already exists, skipping...'.format(self.fn)
                )
                if self.mask:
                    if not os.path.exists(self.mask_fn):
                        [x for x in self.yield_xyz()]
                        self._process(fn=self.mask_fn, filter_=False)
                return(self)
            
        self.run()
        if self.mask:
            if os.path.exists(self.mask_fn):
                self._process(fn=self.mask_fn, filter_=False)
                
        self.waffled = True
        if self.valid_p():
            return(self._process(filter_=True))
        else:
            return(self)
        
    def valid_p(self):
        """check if the output WAFFLES DEM is valid"""

        if not os.path.exists(self.fn):
            return(False)
        
        if self.mask:
            if not os.path.exists(self.mask_fn):
                return(False)
            
        gdi = demfun.infos(self.fn, scan=True)
        
        if gdi is not None:
            if np.isnan(gdi['zr'][0]):
                return(False)
        else:
            return(False)
        
        return(True)            

    def run(self):
        """run the WAFFLES module (set via sub-module class)."""
        
        raise(NotImplementedError)
    
class GMTSurface(Waffle):
    """Waffles GMT surface module.
    Grid the data using GMT 'surface'.
    Data passes through GMT 'blockmean' using weighted mean value if self.weights is True
    """
    
    def __init__(self, tension = .35, relaxation = 1.2,
                 lower_limit = 'd', upper_limit = 'd',
                 breakline=None, **kwargs):
        """generate a DEM with GMT surface"""

        super().__init__(**kwargs)
        
        if self.gc['GMT'] is None:
            utils.echo_error_msg(
                'GMT must be installed to use the SURFACE module'
            )
            return(None, -1)
        
        self.tension = tension
        self.relaxation = relaxation
        self.lower_limit = lower_limit
        self.upper_limit = upper_limit
        self.breakline = breakline

        out, status = utils.run_cmd(
            'gmt gmtset IO_COL_SEPARATOR = SPACE',
            verbose = False
        )
        
        self.mod = 'surface'
        self.mod_args = {
            'tension': self.tension,
            'relaxation': self.relaxation,
            'lower_limit': self.lower_limit,
            'upper_limit': self.upper_limit,
            'breakline': self.breakline,
        }
        self._set_config()
        
    def run(self):        
        dem_surf_cmd = (
            'gmt blockmean {} -I{:.10f}/{:.10f}{} -V -r | gmt surface -V {} -I{:.10f}/{:.10f} -G{}.tif=gd+n-9999:GTiff -T{} -Z{} -Ll{} -Lu{} -r{}'.format(
                self.p_region.format('gmt'),
                self.xinc,
                self.yinc,
                ' -W' if self.weights else '',
                self.p_region.format('gmt'),
                self.xinc,
                self.yinc,
                self.name,
                self.tension,
                self.relaxation,
                self.lower_limit,
                self.upper_limit,
                ' -D{}'.format(self.breakline) if self.breakline is not None else '',
            )
        )
        
        out, status = utils.run_cmd(
            dem_surf_cmd,
            verbose=self.verbose,
            data_fun=lambda p: self.dump_xyz(
                dst_port=p, encode=True
            )
        )
        
        return(self)

class GMTTriangulate(Waffle):
    """Waffles GMT triangulate module
    Grid the data using GMT 'triangulate'.
    """
    
    def __init__(self, **kwargs):
        """generate a DEM with GMT triangulate"""

        super().__init__(**kwargs)
        
        if self.gc['GMT'] is None:
            utils.echo_error_msg(
                'GMT must be installed to use the TRIANGULATE module'
            )
            return(None, -1)

        out, status = utils.run_cmd(
            'gmt gmtset IO_COL_SEPARATOR = SPACE',
            verbose = False
        )
        
        self.mod = 'triangulate'
        self.mod_args = {}
        self._set_config()
        
    def run(self):
        dem_tri_cmd = 'gmt triangulate -V {} -I{:.10f}/{:.10f} -G{}.tif=gd:GTiff -r'.format(
            self.p_region.format('gmt'), self.xinc, self.yinc, self.name
        )
        out, status = utils.run_cmd(
            dem_tri_cmd,
            verbose = self.verbose,
            data_fun = lambda p: self.dump_xyz(
                dst_port=p, encode=True
            )
        )
        
        return(self)

class WafflesMBGrid(Waffle):
    """Waffles MBGrid module
    Does not use internal datalists processing,
    input datalist must be compatible with MB-SYSTEM.
    """
    
    def __init__(self, dist='10/3', tension=35, use_datalists=False, nc=False, **kwargs):
        super().__init__(**kwargs)

        if self.gc['MBGRID'] is None:
            utils.echo_error_msg(
                'MB-System must be installed to use the MBGRID module'
            )
            return(None, -1)

        if self.gc['GMT'] is None:
            utils.echo_error_msg(
                'GMT must be installed to use the MBGRID module'
            )
            return(None, -1)

        self.nc = nc
        self.dist = dist
        self.tension = tension
        self.use_datalists = use_datalists
        self.mod = 'mbgird'
        self.mod_args = {
            'dist': self.dist,
            'tension': self.tension,
            'use_datalists': self.use_datalists
        }
        self._set_config()
                
    def _gmt_num_msk(self, num_grd, dst_msk):
        """generate a num-msk from a NUM grid using GMT grdmath

        Args:
          num_grd (str): pathname to a source `num` grid file
          dst_msk (str): pathname to a destination `msk` grid file

        Returns:
          list: [cmd-output, cmd-return-code]
        """

        num_msk_cmd = 'gmt grdmath -V {} 0 MUL 1 ADD 0 AND = {}'.format(
            num_grd, dst_msk
        )
        return(
            utils.run_cmd(
                num_msk_cmd, verbose=self.verbose
            )
        )

    def _gmt_grd2gdal(self, src_grd, dst_fmt='GTiff'):
        """convert the grd file to tif using GMT

        Args:
          src_grd (str): a pathname to a grid file
          dst_fmt (str): the output GDAL format string

        Returns:
          str: the gdal file name or None
        """

        dst_gdal = '{}.{}'.format(
            os.path.basename(src_grd).split('.')[0], utils.gdal_fext(dst_fmt)
        )
        
        grd2gdal_cmd = 'gmt grdconvert {} {}=gd+n-9999:{} -V'.format(
            src_grd, dst_gdal, dst_fmt
        )
        
        out, status = utils.run_cmd(
            grd2gdal_cmd, verbose=self.verbose
        )
        
        if status == 0:
            return(dst_gdal)
        
        else: return(None)
    
    def run(self):
        # if use_datalists:
        #     #datalist_archive(wg, arch_dir = '.mb_tmp_datalist', verbose = True)
        #     archive = wg['archive']
        #     wg['archive'] = True
        #     for xyz in waffles_yield_datalist(wg): pass
        #     wg['datalist'] = datalists.datalist_major(['archive/{}.datalist'.format(wg['name'])])
        #     wg['archive'] = archive

        mb_region = self.p_region
        if self.node == 'pixel':
            #mb_region = mb_region.buffer(self.xinc * -.5)
            mb_region = mb_region.buffer(x_bv=(self.xinc*-.5), y_bv=(self.yinc*-.5))
        xsize, ysize, gt = mb_region.geo_transform(x_inc=self.xinc)

        ## -G100 breaks mbgrid >= 5.7.8
        mbgrid_cmd = 'mbgrid -I{} {} -D{}/{} -O{} -A2 -F1 -N -C{} -S0 -X0.1 -T{} {}'.format(
            self.data[0].fn,
            mb_region.format('gmt'),
            xsize,
            ysize,
            self.name,
            self.dist,
            self.tension,
            '-M' if self.mask else ''
        )
        for out in utils.yield_cmd(mbgrid_cmd, verbose=self.verbose):
            sys.stderr.write('{}'.format(out))
        #out, status = utils.run_cmd(mbgrid_cmd, verbose=self.verbose)

        if not self.node == 'pixel':
            self._gmt_grd2gdal('{}.grd'.format(self.name))
            utils.remove_glob('*.cmd', '*.mb-1', '{}.grd'.format(self.name))
        else:
            utils.remove_glob('*.cmd', '*.mb-1')
        if self.use_datalists and not self.archive:
            utils.remove_glob('archive')

        if self.mask:
            num_grd = '{}_num.grd'.format(self.name)
            dst_msk = '{}_m.tif=gd+n-9999:GTiff'.format(self.name)
            self.mask_fn = dst_msk
            out, status = self._gmt_num_msk(
                num_grd, dst_msk, verbose=self.verbose
            )
            
            utils.remove_glob(num_grd, '*_sd.grd')
            if not self.use_datalists:
                if self.spat or self.archive:
                    [x for x in self.yield_xyz()]
                    
        return(self)
    
class WafflesNum(Waffle):
    """Generate an Uninterpolated grid from XYZ data;
    num methods include: mask, mean, num, landmask and any gmt grd2xyz -A option.
    """
    
    def __init__(self, mode='n', **kwargs):
        """generate an uninterpolated Grid
        `mode` of `n` generates a num grid
        `mode` of `m` generates a mean grid
        `mode` of `k` generates a mask grid
        `mode` of `w` generates a wet/dry mask grid
        `mode` of `A` generates via GMT 'xyz2grd'
        """

        super().__init__(**kwargs)        
        self.mode = mode
        self.mod = 'num'
        self.mod_args = {'mode': self.mode}
        self._set_config()
        
    def _xyz_num(self, src_xyz):
        """Create a GDAL supported grid from xyz data """

        xcount, ycount, dst_gt = self.p_region.geo_transform(x_inc=self.xinc, y_inc=self.yinc)
        if self.verbose:
            progress = utils.CliProgress(
                'generating uninterpolated NUM grid `{}` @ {}/{}'.format(
                    self.mode, ycount, xcount
                )
            )
            
        if self.mode == 'm' or self.mode == 'w':
            sum_array = np.zeros((ycount, xcount))
            
        count_array = np.zeros((ycount, xcount))
        gdt = gdal.GDT_Float32
        ds_config = demfun.set_infos(
            xcount,
            ycount,
            xcount * ycount,
            dst_gt,
            utils.sr_wkt(self.epsg),
            gdt,
            -9999,
            'GTiff'
        )
        
        for this_xyz in src_xyz:
            #this_xyz = t.copy()
            if regions.xyz_in_region_p(this_xyz, self.p_region):
                xpos, ypos = utils._geo2pixel(
                    this_xyz.x, this_xyz.y, dst_gt
                )
                
                try:
                    if self.mode == 'm' or self.mode == 'w':
                        sum_array[ypos, xpos] += this_xyz.z
                    if self.mode == 'n' or self.mode == 'm' or self.mode == 'w':
                        count_array[ypos, xpos] += 1
                    else:
                        count_array[ypos, xpos] = 1
                        
                except Exception as e:
                    pass
                
        if self.mode == 'm' or self.mode == 'w':
            count_array[count_array == 0] = np.nan
            out_array = (sum_array/count_array)
            if self.mode == 'w':
                out_array[out_array >= 0] = 1
                out_array[out_array < 0] = 0
        else:
            out_array = count_array

        out_array[np.isnan(out_array)] = -9999
        if self. verbose:
            progress.end(
                0,
                'generated uninterpolated num grid `{}` @ {}/{}'.format(
                    self.mode, ycount, xcount
                )
            )
            
        utils.gdal_write(out_array, self.fn, ds_config)
        
    def run(self):
        if self.mode.startswith('A'):
            if self.gc['GMT'] is None:
                utils.echo_error_msg(
                    'GMT must be installed to use the Mode `A` with the NUM module'
                )
                return(None, -1)

            out, status = utils.run_cmd('gmt gmtset IO_COL_SEPARATOR = SPACE', verbose = False)
            
            dem_xyz2grd_cmd = 'gmt xyz2grd -{} -V {} -I{:.10f}/{:.10f} -G{}.tif=gd:GTiff -r{}'.format(
                self.mode,
                self.p_region.format('gmt'),
                self.xinc,
                self.yinc,
                self.name,
                ' -Wi' if self.weights else ''
            )
            
            out, status = utils.run_cmd(
                dem_xyz2grd_cmd,
                verbose=self.verbose,
                data_fun=lambda p: self.dump_xyz(
                    dst_port=p, encode=True
                )
            )
        else:
            #self._xyz_num(self.yield_xyz(block=True))
            self._xyz_num(self.yield_xyz(block=False))
            
        return(self)

class WafflesIDW(Waffle):
    """Generate an IDW DEM.
    If self.weights is True, use weights as uncertainty values
    to generate a UIDW DEM as described here: 
    https://ir.library.oregonstate.edu/concern/graduate_projects/79407x932
    """
    
    def __init__(
            self, radius='1', power=2, block=False, min_points=None, **kwargs
    ):
        super().__init__(**kwargs)
        self.radius = utils.str2inc(radius)
        self.power = utils.float_or(power)
        self.block_p = block
        self.min_points = utils.int_or(min_points)
        self.mod = 'IDW'
        self.mod_args = {
            'radius':self.radius,
            'power':self.power,
            'block':self.block_p,
            'min_points':self.min_points,
        }
        self._set_config()
        
    def _distance(self, pnt0, pnt1):
        return(math.sqrt(sum([(a-b)**2 for a, b in zip(pnt0, pnt1)])))
    
    def run(self):
        xcount, ycount, dst_gt = self.p_region.geo_transform(x_inc=self.xinc, yinc=self.yinc)
        ds_config = demfun.set_infos(
            xcount,
            ycount,
            xcount * ycount,
            dst_gt,
            utils.sr_wkt(self.epsg),
            gdal.GDT_Float32,
            -9999,
            self.fmt
        )
        
        outArray = np.empty((ycount, xcount))
        outArray[:] = np.nan

        if self.verbose:
            if self.min_points:
                progress = utils.CliProgress(
                    'generating IDW grid @ {} and {}/{} looking for at least {} volunteers'.format(
                        self.radius, ycount, xcount, self.min_points
                    )
                )
            else:
                progress = utils.CliProgress(
                    'generating IDW grid @ {} and {}/{}'.format(
                        self.radius, ycount, xcount
                    )
                )
            i=0

        self._xyz_block_t(self.yield_xyz(block=self.block_p))
        for y in range(0, ycount):
            if self.verbose:
                i+=1
                progress.update_perc((i, ycount))
                
            for x in range(0, xcount):
                xyz_bucket = []
                z_list = []
                dw_list = []
                if self.weights:
                    ww_list = []

                if self.min_points is not None:
                    l_radius = 0
                    while len(xyz_bucket) < self.min_points:
                        l_radius += self.radius
                        xg, yg = utils._pixel2geo(x, y, dst_gt)
                        block_region = regions.Region(
                            xmin=xg-l_radius,
                            xmax=xg+l_radius,
                            ymin=yg-l_radius,
                            ymax=yg+l_radius
                        )
                        
                        srcwin = block_region.srcwin(dst_gt, xcount, ycount)
                        for y_i in range(srcwin[1], srcwin[1] + srcwin[3], 1):
                            for x_i in range(srcwin[0], srcwin[0] + srcwin[2], 1):
                                for b in self.block_t[y_i, x_i]:
                                    xyz_bucket.append(b)
                else:
                    xg, yg = utils._pixel2geo(x, y, dst_gt)
                    block_region = regions.Region(
                        xmin=xg-self.radius,
                        xmax=xg+self.radius,
                        ymin=yg-self.radius,
                        ymax=yg+self.radius
                    )
                    
                    srcwin = block_region.srcwin(dst_gt, xcount, ycount)
                    for y_i in range(srcwin[1], srcwin[1] + srcwin[3], 1):
                        for x_i in range(srcwin[0], srcwin[0] + srcwin[2], 1):
                            for b in self.block_t[y_i, x_i]:
                                xyz_bucket.append(b)
                                
                for this_xyz in xyz_bucket:
                    d = self._distance([this_xyz.x, this_xyz.y], [xg, yg])
                    z_list.append(this_xyz.z)
                    if d > 0:
                        dw_list.append(1./(d**self.power))
                    else: dw_list.append(0)
                    
                    if self.weights:
                        w = this_xyz.w
                        if w > 0:
                            ww_list.append(w**self.power)
                        else: ww_list.append(0)

                # if 0 not in dw_list:
                #     if len(dw_list) <= 3:
                #         dw_list = []
                        
                if len(dw_list) > 0:
                    dwt = np.transpose(dw_list)
                    if self.weights:
                        wwt = np.transpose(ww_list)
                        den = sum(np.array(dw_list)*np.array(ww_list))
                        if den > 0:
                            outArray[y,x] = np.dot(
                                z_list, (np.array(dwt)*np.array(wwt))
                            )/den
                        else: outArray[y,x] = sum(z_list)/len(z_list)
                        
                    else:
                        dw_list_sum = sum(dw_list)
                        if dw_list_sum > 0:
                            outArray[y,x] = np.dot(z_list, dwt)/dw_list_sum
                        else: outArray[y,x] = sum(z_list)/len(z_list)
        if self.verbose:
            progress.end(
                0,
                'generated IDW grid {}/{}'.format(
                    ycount, xcount
                )
            )
            
        ds = None            
        outArray[np.isnan(outArray)] = -9999
        out, status = utils.gdal_write(
            outArray, '{}.tif'.format(self.name), ds_config
        )
        
        return(self)

class WafflesIDW2(Waffle):
    """Generate an IDW DEM.
    If self.weights is True, use weights as uncertainty values
    to generate a UIDW DEM as described here: 
    https://ir.library.oregonstate.edu/concern/graduate_projects/79407x932
    """
    
    def __init__(
            self, radius='1', power=2, block=True, min_points=None, **kwargs
    ):
        super().__init__(**kwargs)
        self.radius = utils.str2inc(radius)
        self.power = utils.float_or(power)
        self.block_p = block
        self.min_points = utils.int_or(min_points)
        self.mod = 'IDW'
        self.mod_args = {
            'radius':self.radius,
            'power':self.power,
            'block':self.block_p,
            'min_points':self.min_points,
        }
        self._set_config()
        
    def _distance(self, pnt0, pnt1):
        return(math.sqrt(sum([(a-b)**2 for a, b in zip(pnt0, pnt1)])))
    
    def run(self):
        xcount, ycount, dst_gt = self.p_region.geo_transform(x_inc=self.xinc, y_iinc=self.yinc)
        ds_config = demfun.set_infos(
            xcount,
            ycount,
            xcount * ycount,
            dst_gt,
            utils.sr_wkt(self.epsg),
            gdal.GDT_Float32,
            -9999,
            self.fmt
        )
        
        outArray = np.empty((ycount, xcount))
        outArray[:] = np.nan

        if self.verbose:
            if self.min_points:
                progress = utils.CliProgress(
                    'generating IDW grid @ {} and {}/{} looking for at least {} volunteers'.format(
                        self.radius, ycount, xcount, self.min_points
                    )
                )
            else:
                progress = utils.CliProgress(
                    'generating IDW grid @ {} and {}/{}'.format(
                        self.radius, ycount, xcount
                    )
                )
            i=0

        self._xyz_block_t(self.yield_xyz(block=self.block_p))
        for y in range(0, ycount):
            if self.verbose:
                i+=1
                progress.update_perc((i, ycount))
                
            for x in range(0, xcount):
                xyz_bucket = []
                z_list = []
                dw_list = []
                if self.weights:
                    ww_list = []

                if self.min_points is not None:
                    l_radius = 0
                    while len(xyz_bucket) < self.min_points:
                        l_radius += self.radius
                        xg, yg = utils._pixel2geo(x, y, dst_gt)
                        block_region = regions.Region(
                            xmin=xg-l_radius,
                            xmax=xg+l_radius,
                            ymin=yg-l_radius,
                            ymax=yg+l_radius
                        )
                        
                        srcwin = block_region.srcwin(dst_gt, xcount, ycount)
                        for y_i in range(srcwin[1], srcwin[1] + srcwin[3], 1):
                            for x_i in range(srcwin[0], srcwin[0] + srcwin[2], 1):
                                for b in self.block_t[y_i, x_i]:
                                    xyz_bucket.append(b)
                else:
                    xg, yg = utils._pixel2geo(x, y, dst_gt)
                    block_region = regions.Region(
                        xmin=xg-self.radius,
                        xmax=xg+self.radius,
                        ymin=yg-self.radius,
                        ymax=yg+self.radius
                    )
                    
                    srcwin = block_region.srcwin(dst_gt, xcount, ycount)
                    for y_i in range(srcwin[1], srcwin[1] + srcwin[3], 1):
                        for x_i in range(srcwin[0], srcwin[0] + srcwin[2], 1):
                            for b in self.block_t[y_i, x_i]:
                                xyz_bucket.append(b)
                                
                for this_xyz in xyz_bucket:
                    d = self._distance([this_xyz.x, this_xyz.y], [xg, yg])
                    z_list.append(this_xyz.z)
                    if d > 0:
                        dw_list.append(1./(d**self.power))
                    else: dw_list.append(0)
                    
                    if self.weights:
                        w = this_xyz.w
                        if w > 0:
                            ww_list.append(w**self.power)
                        else: ww_list.append(0)

                if len(dw_list) > 0:
                    dwt = np.transpose(dw_list)
                    if self.weights:
                        wwt = np.transpose(ww_list)
                        den = sum(np.array(dw_list)*np.array(ww_list))
                        if den > 0:
                            outArray[y,x] = np.dot(
                                z_list, (np.array(dwt)*np.array(wwt))
                            )/den
                        else: outArray[y,x] = sum(z_list)/len(z_list)
                        
                    else:
                        dw_list_sum = sum(dw_list)
                        if dw_list_sum > 0:
                            outArray[y,x] = np.dot(z_list, dwt)/dw_list_sum
                        else: outArray[y,x] = sum(z_list)/len(z_list)
        if self.verbose:
            progress.end(
                0,
                'generated IDW grid {}/{}'.format(
                    ycount, xcount
                )
            )
            
        ds = None            
        outArray[np.isnan(outArray)] = -9999
        out, status = utils.gdal_write(
            outArray, '{}.tif'.format(self.name), ds_config
        )
        
        return(self)
    
class WafflesVdatum(Waffle):
    """vertical datum transformation grid via NOAA's VDATUM.
    U.S. and territories only.
    """
    
    def __init__(self, ivert='navd88', overt='mhw', region='4', jar=None, **kwargs):
        """generate a 'conversion-grid' with vdatum.
        
        output will be the differences (surfaced) between 
        `ivert` and `overt` for the region
        
        Args: 
          ivert (str): input vertical datum string
          overt (str): output vertical datum string
          region (str): vdatum grid region
          jar (path): path to vdatum .jar file
        
        Returns:
          list: [{'dem': ['dem-fn', 'raster']}, status]
        """

        super().__init__(**kwargs)
        
        #if self.gc['VDATUM'] is None:
        #    utils.echo_error_msg('VDatum must be installed to use the VDATUM module')
        #    return(None, -1)
        
        self.mod = 'vdatum'
        self.mod_args = {
            'ivert': ivert,
            'overt': overt,
            'region': region,
            'jar': jar
        }
        self.vdc = vdatumfun.Vdatum(
            ivert=ivert,
            overt=overt,
            region=region,
            jar=jar,
            verbose=True
        )
        self._set_config()
        
    def _create_null(
            self, outfile, extent, x_inc, y_inc, nodata, outformat, verbose, overwrite
    ):
        """create a nodata grid"""
        
        xcount, ycount, gt = extent.geo_transform(x_inc=x_inc, y_inc=y_inc)
        ds_config = demfun.set_infos(
            xcount,
            ycount,
            xcount * ycount,
            gt,
            utils.sr_wkt(4326),
            gdal.GDT_Float32,
            nodata,
            outformat
        )
        
        nullArray = np.zeros( (ycount, xcount) )
        nullArray[nullArray==0]=nodata
        utils.gdal_write(nullArray, outfile, ds_config)
        
    def run(self):
        self._create_null(
            'empty.tif', self.p_region, 0.00083333, 0.00083333, 0, 'GTiff', self.verbose, True
        )
        demfun.set_nodata('empty.tif')
        empty = datasets.RasterFile(
            fn='empty.tif',
            data_format=200,
            src_region=self.region,
            epsg=self.epsg,
            name=self.name,
            verbose=self.verbose
        )

        with open('empty.xyz', 'w') as mt_xyz:
            empty.dump_xyz(mt_xyz)
            
        self.vdc.run_vdatum('empty.xyz')
        if os.path.exists('result/empty.xyz') and os.stat('result/empty.xyz').st_size != 0:
            empty_xyz = datasets.XYZFile(
                fn='result/empty.xyz',
                data_format=168,
                src_region=self.region,
                epsg=self.epsg,
                name=self.name,
                verbose=self.verbose
            )
            
            empty_infos = empty_xyz.inf()
            ll = 'd' if empty_infos['minmax'][4] < 0 else '0'
            lu = 'd' if empty_infos['minmax'][5] > 0 else '0'

            GMTSurface(
                data=['result/empty.xyz'],
                name=self.name,
                src_region=self.region,
                xinc=self.xinc,
                yinc=self.yinc,
                epsg=self.epsg,
                verbose=self.verbose,
                tension=0,
                upper_limit=lu,
                lower_limit=ll
            ).generate()
        else:
            utils.echo_error_msg(
                'failed to generate VDatum grid, check settings'
            )
            vd_out = {}
            status = -1

        utils.remove_glob('empty.*', 'result/*', 'result')
        return(self)

class WafflesGDALGrid(Waffle):
    """Waffles GDAL_GRID module.
    see gdal_grid for more info and gridding algorithms
    """
    
    def __init__(self, block=False, **kwargs):
        """run gdal grid using alg_str

        parse the data through xyzfun.xyz_block to get weighted mean before
        building the GDAL dataset to pass into gdal_grid
        
        Args: 
          alg_str (str): the gdal_grid algorithm string
        """
        
        super().__init__(**kwargs)
        self.block_p = block
        self.alg_str = 'linear:radius=-1'
        self.mod = self.alg_str.split(':')[0]
        self.mod_args = {}
        self._set_config()
                
    def run(self):
        
        _prog = utils.CliProgress(
            'running GDAL GRID {} algorithm @ {}...'.format(
                self.alg_str.split(':')[0], self.p_region.format('fn')
            )
        )
        
        _prog_update = lambda x, y, z: _prog.update()
        ds = self._xyz_ds(self.yield_xyz(block=self.block_p))
        
        if ds.GetLayer().GetFeatureCount() == 0:
            utils.echo_error_msg('no input data')
            
        xcount, ycount, dst_gt = self.p_region.geo_transform(x_inc=self.xinc, y_inc=self.yinc)
        
        gd_opts = gdal.GridOptions(
            outputType = gdal.GDT_Float32,
            noData = -9999,
            format = 'GTiff',
            width = xcount,
            height = ycount,
            algorithm = self.alg_str,
            callback = _prog_update if self.verbose else None,
            outputBounds = [
                self.p_region.xmin,
                self.p_region.ymax,
                self.p_region.xmax,
                self.p_region.ymin
            ]
        )
        
        gdal.Grid(
            '{}.tif'.format(self.name), ds, options = gd_opts
        )
        
        demfun.set_nodata(
            '{}.tif'.format(self.name, nodata=-9999, convert_array=False)
        )
        
        _prog.end(
            0,
            'ran GDAL GRID {} algorithm @ {}.'.format(
                self.alg_str.split(':')[0], self.p_region.format('fn')
            )
        )
        
        ds = None
        return(self)

class WafflesLinear(WafflesGDALGrid):

    def __init__(self, radius=None, nodata=-9999, **kwargs):
        super().__init__(**kwargs)
        
        radius = self.xinc * 4 if radius is None else utils.str2inc(radius)
        self.alg_str = 'linear:radius={}:nodata={}'.format(radius, nodata)
        self.mod_args = {'radius':radius, 'nodata':nodata}
        self._set_config()
        
class WafflesInvDst(WafflesGDALGrid):

    def __init__(self, power = 2.0, smoothing = 0.0, radius1 = None, radius2 = None, angle = 0.0,
                   max_points = 0, min_points = 0, nodata = -9999, **kwargs):
        super().__init__(**kwargs)

        radius1 = self.xinc * 2 if radius1 is None else utils.str2inc(radius1)
        radius2 = self.yinc * 2 if radius2 is None else utils.str2inc(radius2)
        self.alg_str = 'invdist:power={}:smoothing={}:radius1={}:radius2={}:angle={}:max_points={}:min_points={}:nodata={}'\
            .format(power, smoothing, radius1, radius2, angle, max_points, min_points, nodata)

        self.mod_args = {
            'power':power,
            'smoothing':smoothing,
            'radius1':radius1,
            'radius2': radius2,
            'angle': angle,
            'max_points': max_points,
            'min_points': min_points,
            'nodata':nodata
        }
        self._set_config()
                
class WafflesMovingAverage(WafflesGDALGrid):

    def __init__(self, radius1=None, radius2=None, angle=0.0, min_points=0, nodata=-9999, **kwargs):
        super().__init__(**kwargs)

        radius1 = self.xinc * 2 if radius1 is None else utils.str2inc(radius1)
        radius2 = self.yinc * 2 if radius2 is None else utils.str2inc(radius2)
        self.alg_str = 'average:radius1={}:radius2={}:angle={}:min_points={}:nodata={}'\
            .format(radius1, radius2, angle, min_points, nodata)
        self.mod_args = {
            'radius1':radius1,
            'radius2': radius2,
            'angle': angle,
            'min_points': min_points,
            'nodata':nodata
        }
        self._set_config()
                
class WafflesNearest(WafflesGDALGrid):

    def __init__(self, radius1=None, radius2=None, angle=0.0, nodata=-9999, **kwargs):
        super().__init__(**kwargs)

        radius1 = self.xinc * 2 if radius1 is None else utils.str2inc(radius1)
        radius2 = self.yinc * 2 if radius2 is None else utils.str2inc(radius2)
        self.alg_str = 'nearest:radius1={}:radius2={}:angle={}:nodata={}'\
            .format(radius1, radius2, angle, nodata)
        self.mod_args = {
            'radius1':radius1,
            'radius2': radius2,
            'angle': angle,
            'nodata':nodata
        }
        self._set_config()

class WafflesCUDEM(Waffle):
    def __init__(
            self,
            want_landmask=False,
            coastline=None,
            min_weight=1,
            smoothing=10,
            radius=None,
            bathy_xinc=None,
            bathy_yinc=None,
            idw=False,
            mask_z=None,
            **kwargs
    ):

        try:
            super().__init__(**kwargs)
        except Exception as e:
            utils.echo_error_msg(e)
            sys.exit()
        
        self.coast_xyz = '{}_coast.xyz'.format(self.name)
        self.want_landmask = want_landmask
        self.coastline = coastline
        self.mask_p = True if self.coastline is not None or self.want_landmask else False
        self.mask_z = utils.float_or(mask_z)
        
        self.min_weight = utils.float_or(min_weight)
        self.smoothing = utils.int_or(smoothing)
        #self.radius = self.inc*9 if radius is None else radius
        self.radius = self.xinc*9 if radius is None else radius
        self.bathy_xinc = self.xinc*3 if bathy_xinc is None else bathy_xinc
        self.bathy_yinc = self.yinc*3 if bathy_yinc is None else bathy_yinc
        self.idw = idw
        self.bathy_data = self.data_
        self.mod = 'cudem'
        self.mod_args = {}
        self._set_config()
        
    def run(self):

        surface_region = self.p_region.copy()
        surface_region.wmin = self.min_weight

        bathy_region = self.p_region.copy()
        if self.mask_z is not None:
            bathy_region.zmax = self.mask_z + 10
        bathy_clip = None
        #bathy_data = self.data_

        if self.idw:
            idw_region = self.p_region.copy()
            if self.mask_p:
                idw_region.zmax = self.mask_z

            self.idw = WaffleFactory(
                mod='IDW:radius={}'.format(self.radius),
                data=self.data_,
                src_region=idw_region,
                xinc=utils.str2inc(self.bathy_xinc),
                yinc=utils.str2inc(self.bathy_yinc),
                name='tmp_idw',
                node=self.node,
                extend=self.extend+6,
                extend_proc=self.extend_proc+6,
                #fltr=['1:{}'.format(self.smoothing)] if self.smoothing is not None and not self.mask_p else [],
                weights=1,
                epsg=self.epsg,
                clobber=True,
                verbose=self.verbose,
            ).acquire().generate()

            #self.bathy_data = self.data_ + ['{},200,{}'.format(self.idw.fn, self.min_weight)],
            self.bathy_data = ['{},200,{}'.format(self.idw.fn, self.min_weight)]

        if self.mask_p:

            #c_cmd = 'coastline2xyz.sh -I {} -O {} -Z 0 -W {} -E {} -S {} -N {}'.format(
            #    self.coastline,
            #    self.coast_xyz,
            #    self.c_region.xmin,
            #    self.c_region.xmax,
            #    self.c_region.ymin,
            #    self.c_region.ymax)
            #out, status = utils.run_cmd(c_cmd, verbose=True)
            #bathy_data = self.data_ + ['{},168,0.1'.format(self.coast_xyz)]

            #c_cmd = 'coastline2xyz.sh -I {} -O {} -Z 0 -W {} -E {} -S {} -N {}'.format(
            #    'tmp_coast.shp',
            #    self.coast_xyz,
            #    self.c_region.xmin,
            #    self.c_region.xmax,
            #    self.c_region.ymin,
            #    self.c_region.ymax)
            #out, status = utils.run_cmd(c_cmd, verbose=True)
            #bathy_data = self.data_ + ['{},168,0.1'.format(self.coast_xyz)]
            
            if self.coastline is not None:            
                bathy_clip = '{}:invert=True'.format(self.coastline)
            elif self.want_landmask:
                self.coast = WaffleFactory(
                    mod='coastline',
                    data=self.data_,
                    src_region=self.c_region,
                    xinc=self.xinc,
                    yinc=self.yinc,
                    name='tmp_coast',
                    node=self.node,
                    extend=self.extend+12,
                    weights=self.weights,
                    epsg=self.epsg,
                    clobber=True,
                    verbose=self.verbose,
                ).acquire().generate()
                bathy_clip = '{}:invert=True'.format(self.coast.name + '.shp')
        
            self.bathy = WaffleFactory(
                mod='surface:tension=1:upper_limit={}'.format(self.mask_z),
                #data=bathy_data + ['{},200,{}'.format(self.idw.fn, self.min_weight)],
                #data=['{},200,{}'.format(self.idw.fn, self.min_weight)],
                data=self.bathy_data,
                src_region=bathy_region,
                xinc=utils.str2inc(self.bathy_xinc),
                yinc=utils.str2inc(self.bathy_yinc),
                name='tmp_bathy',
                node=self.node,
                extend=self.extend+2,
                extend_proc=self.extend_proc+2,
                fltr=['1:{}'.format(self.smoothing)] if self.smoothing is not None else [],
                weights=1,
                epsg=self.epsg,
                clobber=True,
                verbose=self.verbose,
                xsample=utils.str2inc(self.xinc),
                ysample=utils.str2inc(self.yinc),
                clip=bathy_clip,
            ).acquire().generate()
        else:
            self.bathy = WaffleFactory(
                mod='surface:tension=1',
                #data=bathy_data + ['{},200,{}'.format(self.idw.fn, self.min_weight)],
                #data=['{},200,{}'.format(self.idw.fn, self.min_weight)],
                data=self.bathy_data,
                src_region=bathy_region,
                xinc=utils.str2inc(self.bathy_xinc),
                yinc=utils.str2inc(self.bathy_yinc),
                name='tmp_bathy',
                node=self.node,
                extend=self.extend+2,
                extend_proc=self.extend_proc+2,
                fltr=['1:{}'.format(self.smoothing)] if self.smoothing is not None else [],
                weights=1,
                epsg=self.epsg,
                clobber=True,
                verbose=self.verbose,
                #sample=utils.str2inc(self.inc),
                clip=bathy_clip,
            ).acquire().generate()
            
            #self.bathy = self.idw
        
        self.surface = WaffleFactory(
            mod='surface:tension=1',
            data=self.data_ + ['{},200,{}'.format(self.bathy.fn, self.min_weight)],
            src_region=surface_region,
            xinc=self.xinc,
            yinc=self.yinc,
            name=self.name,
            node=self.node,
            extend=self.extend,
            extend_proc=self.extend_proc,
            weights=self.weights,
            epsg=self.epsg,
            clobber=True,
            verbose=self.verbose,
        ).acquire().generate()

        utils.remove_glob('tmp_bathy*', 'tmp_coast*', 'tmp_idw*', '{}*'.format(self.coast_xyz))
        
        return(self)
        
class WafflesCoastline(Waffle):
    def __init__(self, want_nhd=True, want_gmrt=False, invert=False, polygonize=True, **kwargs):
        """Generate a landmask from various sources."""

        super().__init__(**kwargs)
        self.want_nhd = want_nhd
        self.want_gmrt = want_gmrt
        self.invert = invert
        self.polygonize = polygonize

        self.coast_array = None
        self.ds_config = None

        import cudem.fetches.copernicus
        import cudem.fetches.tnm
        import cudem.fetches.gmrt
        import cudem.fetches.utils

        self.f_region = self.p_region.copy()
        self.f_region.buffer(x_bv=(self.xinc*10), y_bv=(self.yinc*10))
        self.f_region.epsg = self.epsg
        self.wgs_region = self.f_region.copy()
        self.wgs_region.warp(4326)
        
        self.mod = 'coastline'
        self.mod_args = {
            'want_nhd': want_nhd,
            'want_gmrt': want_gmrt,
            'polygonize': polygonize,
            'invert': invert,
        }

        self._set_config()

    def run(self):
        self._load_coast_mask()

        if self.want_gmrt:
            self._load_gmrt()
        else:
            self._load_copernicus()

        if self.want_nhd:
            self._load_nhd()
        
        if len(self.data) > 0:
            self._load_data()

        utils.echo_msg(
            'finanlizing array for region {} at {} {}...'.format(
                self.p_region.format('gmt'), self.ds_config['nx'], self.ds_config['ny']
            )
        )
        self._finalize_array()
        utils.echo_msg('writing array to {}.tif...'.format(self.name))
        
        self._write_coast_array()
        
        if self.polygonize:
            self._write_coast_poly()
            
        return(self)

    def _finalize_array(self):
        #self.coast_array[self.coast_array > 0] = 0 if self.invert else 1
        #self.coast_array[self.coast_array <= 0] = 1 if self.invert else 0

        self.coast_array[self.coast_array > 0] = 1
        self.coast_array[self.coast_array <= 0] = 0

        if self.invert:
            self.coast_array[self.coast_array == 0] = 2
            self.coast_array[self.coast_array == 1] = 0
            self.coast_array[self.coast_array == 2] = 1
            
    def _load_coast_mask(self):
        """create a nodata grid"""
        
        xcount, ycount, gt = self.p_region.geo_transform(x_inc=self.xinc, y_inc=self.yinc)
        self.ds_config = demfun.set_infos(
            xcount,
            ycount,
            xcount * ycount,
            gt,
            utils.sr_wkt(self.epsg),
            gdal.GDT_Int32,
            -9999,
            'GTiff'
        )
        
        self.coast_array = np.zeros( (ycount, xcount) )

    def _load_gmrt(self):
        """GMRT - Global low-res
        Used to fill un-set cells.
        """
        
        this_gmrt = cudem.fetches.gmrt.GMRT(
            src_region=self.f_region, weight=self.weights, verbose=self.verbose, layer='topo-mask'
        )
        this_gmrt._outdir = './'
        this_gmrt.run()

        fr = cudem.fetches.utils.fetch_results(this_gmrt, want_proc=False)
        fr.daemon = True
        fr.start()
        fr.join()
        
        #this_gmrt.fetch_results()
        gmrt_tif = this_gmrt.results[0]

        dst_srs = osr.SpatialReference()
        dst_srs.ImportFromEPSG(int(self.epsg))
        
        driver = gdal.GetDriverByName('MEM')
        out_ds = driver.Create('MEM', self.ds_config['nx'], self.ds_config['ny'], 1, self.ds_config['dt'])
            
        out_ds.SetGeoTransform(self.ds_config['geoT'])
        out_ds.SetProjection(self.ds_config['proj'])
        out_ds.GetRasterBand(1).SetNoDataValue(self.ds_config['ndv'])
            
        gdal.Warp(out_ds, gmrt_tif[1], dstSRS = dst_srs, resampleAlg = gdal.GRA_CubicSpline)
            
        gmrt_ds_arr = out_ds.GetRasterBand(1).ReadAsArray()

        gmrt_ds_arr[gmrt_ds_arr > 0] = 1
        gmrt_ds_arr[gmrt_ds_arr <= 0] = 0
        self.coast_array += gmrt_ds_arr

        out_ds = gmrt_ds_arr = None

        utils.remove_glob(gmrt_tif[1])
        
    def _load_copernicus(self):
        """copernicus"""

        this_cop = cudem.fetches.copernicus.CopernicusDEM(
            src_region=self.wgs_region, weight=self.weights, verbose=self.verbose, datatype='1'
        )

        this_cop._outdir = './'
        this_cop.run()

        fr = cudem.fetches.utils.fetch_results(this_cop, want_proc=False)
        fr.daemon = True
        fr.start()
        fr.join()

        #this_cop.fetch_results()
        
        dst_srs = osr.SpatialReference()
        dst_srs.ImportFromEPSG(int(self.epsg))
        
        for i, cop_tif in enumerate(this_cop.results):

            driver = gdal.GetDriverByName('MEM')
            out_ds = driver.Create('MEM', self.ds_config['nx'], self.ds_config['ny'], 1, self.ds_config['dt'])
            
            out_ds.SetGeoTransform(self.ds_config['geoT'])
            out_ds.SetProjection(self.ds_config['proj'])
            out_ds.GetRasterBand(1).SetNoDataValue(self.ds_config['ndv'])
            
            gdal.Warp(out_ds, cop_tif[1], dstSRS = dst_srs, resampleAlg = gdal.GRA_CubicSpline)
            
            c_ds_arr = out_ds.GetRasterBand(1).ReadAsArray()

            c_ds_arr[c_ds_arr != 0] = 1
            self.coast_array += c_ds_arr

            out_ds = c_ds_arr = None
            
            utils.remove_glob(cop_tif[1]) #, '{}*'.format(out))            
    
    def _load_nhd(self):
        """USGS NHD (HIGH-RES U.S. Only)
        Fetch NHD (NHD High/Plus) data from TNM to fill in near-shore areas. 
        High resoultion data varies by location...
        """

        this_tnm = cudem.fetches.tnm.TheNationalMap(
            src_region=self.wgs_region,
            weight=self.weights,
            verbose=self.verbose,
            where=["Name LIKE '%Hydro%'"],
            extents='HU-8 Subbasin,HU-4 Subregion'
        )
        this_tnm._outdir = './'
        this_tnm.run()

        fr = cudem.fetches.utils.fetch_results(this_tnm, want_proc=False)
        fr.daemon = True
        fr.start()
        fr.join()
        #this_tnm.fetch_results()
            
        if len(this_tnm.results) > 0:

            for i, tnm_zip in enumerate(this_tnm.results):
                tnm_zips = utils.unzip(tnm_zip[1])
                gdb = tnm_zip[1].split('.')[:-1][0] + '.gdb'
                utils.run_cmd(
                    'ogr2ogr -update -append nhdArea_merge.shp {} NHDArea -where "FType=312 OR FType=336 OR FType=445 OR FType=460 OR FType=537" -clipdst {} 2>/dev/null'.format(
                        gdb, self.p_region.format('ul_lr')
                    ),
                    verbose=True
                )

                utils.run_cmd(
                    'ogr2ogr -update -append nhdArea_merge.shp {} NHDPlusBurnWaterBody -clipdst {} 2>/dev/null'.format(
                        gdb, self.p_region.format('ul_lr')
                    ),
                    verbose=True
                )

                utils.run_cmd(
                    'ogr2ogr -update -append nhdArea_merge.shp {} NHDWaterBody -where "FType=493 OR FType=466" -clipdst {} 2>/dev/null'.format(
                        gdb, self.p_region.format('ul_lr')
                    ),
                    verbose=True
                )

                utils.remove_glob(tnm_zip[1], *tnm_zips)

            utils.run_cmd(
                'gdal_rasterize -burn 1 nhdArea_merge.shp nhdArea_merge.tif -te {} -ts {} {} -ot Int32'.format(
                    self.p_region.format('te'),
                    self.ds_config['nx'],
                    self.ds_config['ny'],
                ),
                verbose=True
            )

            tnm_ds = gdal.Open('nhdArea_merge.tif')
            if tnm_ds is not None:
                tnm_ds_arr = tnm_ds.GetRasterBand(1).ReadAsArray()

                tnm_ds_arr[tnm_ds_arr != 1] = -1

                self.coast_array -= tnm_ds_arr

                tnm_ds = tnm_ds_arr = None
            utils.remove_glob('nhdArea_merge.*')
        
    def _load_data(self):
        """load data from user datalist and amend coast_array"""
        
        for this_data in self.data:
            for this_xyz in this_data.yield_xyz():
                xpos, ypos = utils._geo2pixel(
                    this_xyz.x, this_xyz.y, self.ds_config['geoT']
                )
                try:
                    if this_xyz.z >= 0:
                        self.coast_array[ypos, xpos] += 1
                    else:
                        self.coast_array[ypos, xpos] -= 1
                except: pass
        
    def _write_coast_array(self):
        """write coast_array to file
        1 = land
        0 = water
        """
        
        utils.gdal_write(
            self.coast_array, '{}.tif'.format(self.name), self.ds_config,
        )

    def _write_coast_poly(self):
        """convert to coast_array vector"""

        tmp_ds = ogr.GetDriverByName('ESRI Shapefile').CreateDataSource(
            'tmp_c_{}.shp'.format(self.name)
        )
        
        if tmp_ds is not None:
            tmp_layer = tmp_ds.CreateLayer('tmp_c_{}'.format(self.name), None, ogr.wkbMultiPolygon)
            tmp_layer.CreateField(ogr.FieldDefn('DN', ogr.OFTInteger))
            demfun.polygonize('{}.tif'.format(self.name), tmp_layer, verbose=self.verbose)
            tmp_ds = None
        utils.run_cmd(
            'ogr2ogr -dialect SQLITE -sql "SELECT * FROM tmp_c_{} WHERE DN=0 order by ST_AREA(geometry) desc limit 4" {}.shp tmp_c_{}.shp'.format(
                self.name, self.name, self.name),
            verbose=True
        )
        
        utils.remove_glob('tmp_c_{}.*'.format(self.name))
        utils.run_cmd(
            'ogrinfo -dialect SQLITE -sql "UPDATE {} SET geometry = ST_MakeValid(geometry)" {}.shp'.format(
                self.name, self.name)
        )

        utils.gdal_prj_file(self.name + '.prj', self.epsg)
                
class WaffleFactory():
    """Find and generate a WAFFLE object for DEM generation."""
    
    _modules = {
        'surface': {
            'name': 'surface',
            'datalist-p': True,
            'description': '''SPLINE DEM via GMT surface\n
Generate a DEM using GMT's surface command
        
< surface:tension=.35:relaxation=1.2:lower_limit=d:upper_limit=d >
 :tension=[0-1] - Spline tension.
 :relaxation=[val] - Spline relaxation factor.
 :lower_limit=[val] - Constrain interpolation to lower limit.
 :upper_limit=[val] - Constrain interpolation to upper limit.''',
            'class': lambda k: GMTSurface(**k),
        },
        'triangulate': {
            'name': 'triangulate',
            'datalist-p': True,
            'description': '''TRIANGULATION DEM via GMT triangulate\n
Generate a DEM using GMT's triangulate command.
        
< triangulate >''',
            'class': lambda k: GMTTriangulate(**k),
        },
        'num': {
            'name': 'num',
            'datalist-p': True,
            'description': """Uninterpolated DEM populated by <mode>.\n
Generate an uninterpolated DEM using <mode> option.
Using mode of 'A<mode>' uses GMT's xyz2grd command, 
see gmt xyz2grd --help for more info.

< num:mode=n >
 :mode=[key] - specify mode of grid population: 
  keys: k (mask), m (mean), n (num), w (wet), A<mode> (gmt xyz2grd)""",
            'class': lambda k: WafflesNum(**k),
        },
        'linear': {
            'name': 'linear',
            'datalist-p': True,
            'description': """LINEAR DEM via gdal_grid\n
Generate a DEM using GDAL's gdal_grid command.
see gdal_grid --help for more info

< linear:radius=-1 >
 :radius=[val] - search radius""",
            'class': lambda k: WafflesLinear(**k),
        },
        'nearest': {
            'name': 'nearest',
            'datalist-p': True,
            'description': """NEAREST DEM via gdal_grid\n
Generate a DEM using GDAL's gdal_grid command.
see gdal_grid --help for more info

< nearest:radius1=0:radius2=0:angle=0:nodata=0 >
 :radius1=[val] - search radius
 :radius2=[val] - search radius""",
            'class': lambda k: WafflesNearest(**k),
        },
        'average': {
            'name': 'average',
            'datalist-p': True,
            'description': """Moving AVERAGE DEM via gdal_grid\n
Generate a DEM using GDAL's gdal_grid command.
see gdal_grid --help for more info

< nearest:radius1=0:radius2=0:angle=0:min_points=0:nodata=0 >
 :radius1=[val] - search radius
 :radius2=[val] - search radius
 :min_points=[val] - minimum points per IDW bucket (use to fill entire DEM)""",
            'class': lambda k: WafflesMovingAverage(**k),
        },
        'invdst': {
            'name': 'invdst',
            'datalist-p': True,
            'description': """INVERSE DISTANCE DEM via gdal_grid\n
Generate a DEM using GDAL's gdal_grid command.
see gdal_grid --help for more info

< invdst:power=2.0:smoothing=0.0:radius1=0:radius2=0:angle=0:max_points=0:min_points=0:nodata=0 >
 :radius1=[val] - search radius
 :radius2=[val] - search radius
 :power=[val] - weight**power
 :min_points=[val] - minimum points per IDW bucket (use to fill entire DEM)""",
            'class': lambda k: WafflesInvDst(**k),
        },
        'IDW': {
            'name': 'IDW',
            'datalist-p': True,
            'description': """INVERSE DISTANCE WEIGHTED DEM\n
Generate a DEM using an Inverse Distance Weighted algorithm.
If weights are used, will generate a UIDW DEM, using weight values as inverse uncertainty,
as described here: https://ir.library.oregonstate.edu/concern/graduate_projects/79407x932

< IDW:radius=0:min_points=0:power=2:block=True >
 :radius=[val] - search radius
 :power=[val] - weight**power
 :min_points=[val] - minimum points per IDW bucket (use to fill entire DEM)
 :block=[True/False] - block the data before performing the IDW routine""",
            'class': lambda k: WafflesIDW(**k),
        },
        'vdatum': {
            'name': 'vdatum',
            'datalist-p': False,
            'description': """VDATUM conversion grid.
            """,
            'class': lambda k: WafflesVDatum(**k),
        },
        'mbgrid': {
            'name': 'mbgrid',
            'datalist-p': True,
            'description': """SPLINE DEM via MB-System's mbgrid.\n
Generate a DEM using MB-System's mbgrid command.
see mbgrid --help for more info

< mbgrid:dist='10/3':tension=35:use_datalists=False >
 :dist=[val] - the dist variable to use in mbgrid
 :tension=[val] - the spline tension value (0-inf)
 :use_datalist=[True/False] - use built-in datalists rather than mbdatalist""",
            'class': lambda k: WafflesMBGrid(**k),
        },
        'coastline': {
            'name': 'coastline',
            'datalist-p': False,
            'description': """COASTLINE generation.
Generate a coastline (land/sea mask) using a variety of sources.

< coastline:wet=None:dry=None:want_gmrt=False:invert=False >
 :want_gmrt=[True/False] - use GMRT to fill background
 :want_nhd=[True/False] - use NHD
 :invert=[True/False] - invert the output results
 :polygonize=[True/False] - polygonize the output""",
            'class': lambda k: WafflesCoastline(**k),
        },
        'cudem': {
            'name': 'cudem',
            'datalist-p': True,
            'description': """CUDEM integrated DEM generation.
Generate a DEM using a variety of sources.

< cudem >""",
            'class': lambda k: WafflesCUDEM(**k),
        },
    }
    
    def __init__(
            self,
            mod=None,
            data=[],
            src_region=None,
            inc=None,
            xinc=None,
            yinc=None,
            name='waffles_dem',
            node='pixel',
            fmt='GTiff',
            extend=0,
            extend_proc=20,
            weights=None,
            fltr=[],
            sample=None,
            xsample=None,
            ysample=None,
            clip=None,
            chunk=None,
            epsg=4326,
            verbose=False,
            archive=False,
            mask=False,
            spat=False,
            clobber=True
    ):
        self.mod = mod
        self.mod_args = {}
        self.data = data
        self.region = src_region
        self.inc = inc
        self.xinc = xinc
        self.yinc = yinc
        self.sample = sample
        self.xsample = xsample
        self.ysample = ysample
        self.name = name
        self.node = node
        self.fmt = fmt
        self.extend = extend
        self.extend_proc = extend_proc
        self.weights = weights
        self.fltr = fltr
        self.clip = clip
        self.chunk = chunk
        self.epsg = utils.int_or(epsg)
        self.archive = archive
        self.mask = mask
        self.spat = spat
        self.clobber = clobber
        self.verbose = verbose

        if self.mod is not None:
            self._parse_mod(self.mod)

        self._set_config()
        
    def _parse_mod(self, mod):
        opts = mod.split(':')
        if opts[0] in self._modules.keys():
            if len(opts) > 1:
                self.mod_args = utils.args2dict(list(opts[1:]), {})
            self.mod_name = opts[0]
        else:
            utils.echo_error_msg(
                'invalid module name `{}`'.format(opts[0])
            )
        return(self.mod_name, self.mod_args)
        
    def _set_config(self):
        """export the waffles config info as a dictionary"""
        
        self._config = {
            'mod': self.mod,
            'data': self.data,
            'src_region': self.region.export_as_list(
                include_z=True, include_w=True
            ) if self.region is not None else None,
            'inc': self.inc,
            'xinc': self.xinc,
            'yinc': self.yinc,
            'sample': self.sample,
            'xsample': self.xsample,
            'ysample': self.ysample,
            'name': self.name,
            'node': self.node,
            'fmt': self.fmt,
            'extend': self.extend,
            'extend_proc': self.extend_proc,
            'weights': self.weights,
            'fltr': self.fltr,
            'clip': self.clip,
            'chunk': self.chunk,
            'epsg': self.epsg,
            'verbose': self.verbose,
            'archive': self.archive,
            'mask': self.mask,
            'spat': self.spat,
            'clobber': self.clobber,
        }

        for key in self.mod_args.keys():
            if key not in self._config.keys():
                self._config[key] = self.mod_args[key]

        return(self._config)
            
    def acquire_surface(self, **kwargs):
        return(
            GMTSurface(
                data=self.data,
                src_region=self.region,
                inc=self.inc,
                xinc=self.xinc,
                yinc=self.yinc,
                sample=self.sample,
                xsample=self.xsample,
                ysample=self.ysample,
                name=self.name,
                node=self.node,
                fmt=self.fmt,
                extend=self.extend,
                extend_proc=self.extend_proc,
                weights=self.weights,
                fltr=self.fltr,
                clip=self.clip,
                chunk=self.chunk,
                epsg=self.epsg,
                archive=self.archive,
                mask=self.mask,
                spat=self.spat,
                clobber=self.clobber,
                verbose=self.verbose,
                **kwargs
            )
        )

    def acquire_triangulate(self, **kwargs):
        return(
            GMTTriangulate(
                data=self.data,
                src_region=self.region,
                inc=self.inc,
                xinc=self.xinc,
                yinc=self.yinc,
                xsample=self.xsample,
                ysample=self.ysample,
                name=self.name,
                node=self.node,
                fmt=self.fmt,
                extend=self.extend,
                extend_proc=self.extend_proc,
                weights=self.weights,
                fltr=self.fltr,
                sample=self.sample,
                clip=self.clip,
                chunk=self.chunk,
                epsg=self.epsg,
                archive=self.archive,
                mask=self.mask,
                spat=self.spat,
                clobber=self.clobber,
                verbose=self.verbose,
                **kwargs
            )
        )

    def acquire_num(self, **kwargs):
        return(
            WafflesNum(
                data=self.data,
                src_region=self.region,
                inc=self.inc,
                xinc=self.xinc,
                yinc=self.yinc,
                xsample=self.xsample,
                ysample=self.ysample,
                name=self.name,
                node=self.node,
                fmt=self.fmt,
                extend=self.extend,
                extend_proc=self.extend_proc,
                weights=self.weights,
                fltr=self.fltr,
                sample=self.sample,
                clip=self.clip,
                chunk=self.chunk,
                epsg=self.epsg,
                archive=self.archive,
                mask=self.mask,
                spat=self.spat,
                clobber=self.clobber,
                verbose=self.verbose,
                **kwargs
            )
        )

    def acquire_linear(self, **kwargs):
        return(
            WafflesLinear(
                data=self.data,
                src_region=self.region,
                inc=self.inc,
                xinc=self.xinc,
                yinc=self.yinc,
                xsample=self.xsample,
                ysample=self.ysample,
                name=self.name,
                node=self.node,
                fmt=self.fmt,
                extend=self.extend,
                extend_proc=self.extend_proc,
                weights=self.weights,
                fltr=self.fltr,
                sample=self.sample,
                clip=self.clip,
                chunk=self.chunk,
                epsg=self.epsg,
                archive=self.archive,
                mask=self.mask,
                spat=self.spat,
                clobber=self.clobber,
                verbose=self.verbose,
                **kwargs
            )
        )

    def acquire_average(self, **kwargs):
        return(
            WafflesMovingAverage(
                data=self.data,
                src_region=self.region,
                inc=self.inc,
                xinc=self.xinc,
                yinc=self.yinc,
                xsample=self.xsample,
                ysample=self.ysample,
                name=self.name,
                node=self.node,
                fmt=self.fmt,
                extend=self.extend,
                extend_proc=self.extend_proc,
                weights=self.weights,
                fltr=self.fltr,
                sample=self.sample,
                clip=self.clip,
                chunk=self.chunk,
                epsg=self.epsg,
                archive=self.archive,
                mask=self.mask,
                spat=self.spat,
                clobber=self.clobber,
                verbose=self.verbose,
                **kwargs
            )
        )

    def acquire_invdst(self, **kwargs):
        return(
            WafflesInvDst(
                data=self.data,
                src_region=self.region,
                inc=self.inc,
                xinc=self.xinc,
                yinc=self.yinc,
                xsample=self.xsample,
                ysample=self.ysample,
                name=self.name,
                node=self.node,
                fmt=self.fmt,
                extend=self.extend,
                extend_proc=self.extend_proc,
                weights=self.weights,
                fltr=self.fltr,
                sample=self.sample,
                clip=self.clip,
                chunk=self.chunk,
                epsg=self.epsg,
                archive=self.archive,
                mask=self.mask,
                spat=self.spat,
                clobber=self.clobber,
                verbose=self.verbose,
                **kwargs
            )
        )
    
    def acquire_nearest(self, **kwargs):
        return(
            WafflesNearest(
                data=self.data,
                src_region=self.region,
                inc=self.inc,
                xinc=self.xinc,
                yinc=self.yinc,
                xsample=self.xsample,
                ysample=self.ysample,
                name=self.name,
                node=self.node,
                fmt=self.fmt,
                extend=self.extend,
                extend_proc=self.extend_proc,
                weights=self.weights,
                fltr=self.fltr,
                sample=self.sample,
                clip=self.clip,
                chunk=self.chunk,
                epsg=self.epsg,
                archive=self.archive,
                mask=self.mask,
                spat=self.spat,
                clobber=self.clobber,
                verbose=self.verbose,
                **kwargs
            )
        )

    def acquire_IDW(self, **kwargs):
        return(
            WafflesIDW(
                data=self.data,
                src_region=self.region,
                inc=self.inc,
                xinc=self.xinc,
                yinc=self.yinc,
                xsample=self.xsample,
                ysample=self.ysample,
                name=self.name,
                node=self.node,
                fmt=self.fmt,
                extend=self.extend,
                extend_proc=self.extend_proc,
                weights=self.weights,
                fltr=self.fltr,
                sample=self.sample,
                clip=self.clip,
                chunk=self.chunk,
                epsg=self.epsg,
                archive=self.archive,
                mask=self.mask,
                spat=self.spat,
                clobber=self.clobber,
                verbose=self.verbose,
                **kwargs
            )
        )
    
    def acquire_vdatum(self, **kwargs):
        return(
            WafflesVdatum(
                data=self.data,
                src_region=self.region,
                inc=self.inc,
                xinc=self.xinc,
                yinc=self.yinc,
                xsample=self.xsample,
                ysample=self.ysample,
                name=self.name,
                node=self.node,
                fmt=self.fmt,
                extend=self.extend,
                extend_proc=self.extend_proc,
                weights=self.weights,
                fltr=self.fltr,
                sample=self.sample,
                clip=self.clip,
                chunk=self.chunk,
                epsg=self.epsg,
                archive=self.archive,
                mask=self.mask,
                spat=self.spat,
                clobber=self.clobber,
                verbose=self.verbose,
                **kwargs
            )
        )

    def acquire_mbgrid(self, **kwargs):
        return(
            WafflesMBGrid(
                data=self.data,
                src_region=self.region,
                inc=self.inc,
                xinc=self.xinc,
                yinc=self.yinc,
                xsample=self.xsample,
                ysample=self.ysample,
                name=self.name,
                node=self.node,
                fmt=self.fmt,
                extend=self.extend,
                extend_proc=self.extend_proc,
                weights=self.weights,
                fltr=self.fltr,
                sample=self.sample,
                clip=self.clip,
                chunk=self.chunk,
                epsg=self.epsg,
                archive=self.archive,
                mask=self.mask,
                spat=self.spat,
                clobber=self.clobber,
                verbose=self.verbose,
                **kwargs
            )
        )

    def acquire_coastline(self, **kwargs):
        return(
            WafflesCoastline(
                data=self.data,
                src_region=self.region,
                inc=self.inc,
                xinc=self.xinc,
                yinc=self.yinc,
                xsample=self.xsample,
                ysample=self.ysample,
                name=self.name,
                node=self.node,
                fmt=self.fmt,
                extend=self.extend,
                extend_proc=self.extend_proc,
                weights=self.weights,
                fltr=self.fltr,
                sample=self.sample,
                clip=self.clip,
                chunk=self.chunk,
                epsg=self.epsg,
                archive=self.archive,
                mask=self.mask,
                spat=self.spat,
                clobber=self.clobber,
                verbose=self.verbose,
                **kwargs
            )
        )
    
    def acquire_cudem(self, **kwargs):
        return(
            WafflesCUDEM(
                data=self.data,
                src_region=self.region,
                inc=self.inc,
                xinc=self.xinc,
                yinc=self.yinc,
                xsample=self.xsample,
                ysample=self.ysample,
                name=self.name,
                node=self.node,
                fmt=self.fmt,
                extend=self.extend,
                extend_proc=self.extend_proc,
                weights=self.weights,
                fltr=self.fltr,
                sample=self.sample,
                clip=self.clip,
                chunk=self.chunk,
                epsg=self.epsg,
                archive=self.archive,
                mask=self.mask,
                spat=self.spat,
                clobber=self.clobber,
                verbose=self.verbose,
                **kwargs
            )
        )

    def acquire(self):
        if self.mod_name == 'surface':
            return(self.acquire_surface(**self.mod_args))
        
        if self.mod_name == 'triangulate':
            return(self.acquire_triangulate(**self.mod_args))
        
        if self.mod_name == 'num':
            return(self.acquire_num(**self.mod_args))
        
        if self.mod_name == 'linear':
            return(self.acquire_linear(**self.mod_args))
        
        if self.mod_name == 'average':
            return(self.acquire_average(**self.mod_args))
        
        if self.mod_name == 'nearest':
            return(self.acquire_nearest(**self.mod_args))
        
        if self.mod_name == 'invdst':
            return(self.acquire_invdst(**self.mod_args))
        
        if self.mod_name == 'IDW':
            return(self.acquire_IDW(**self.mod_args))
        
        if self.mod_name == 'vdatum':
            return(self.acquire_vdatum(**self.mod_args))
        
        if self.mod_name == 'mbgrid':
            return(self.acquire_mbgrid(**self.mod_args))
        
        if self.mod_name == 'coastline':
            return(self.acquire_coastline(**self.mod_args))

        if self.mod_name == 'cudem':
            return(self.acquire_cudem(**self.mod_args))
    
    def acquire_module_by_name(self, mod_name, **mod_args):

        #mod_name, mod_args = self._parse_mod(mod)
        
        if mod_name == 'surface':
            return(self.acquire_surface(**mod_args))
        
        if mod_name == 'triangulate':
            return(self.acquire_triangulate(**mod_args))
        
        if mod_name == 'num':
            return(self.acquire_num(**mod_args))
        
        if mod_name == 'linear':
            return(self.acquire_linear(**mod_args))
        
        if mod_name == 'average':
            return(self.acquire_average(**mod_args))
        
        if mod_name == 'nearest':
            return(self.acquire_nearest(**mod_args))
        
        if mod_name == 'invdst':
            return(self.acquire_invdst(**mod_args))
        
        if mod_name == 'IDW':
            return(self.acquire_IDW(**mod_args))
        
        if mod_name == 'vdatum':
            return(self.acquire_vdatum(**mod_args))
        
        if mod_name == 'mbgrid':
            return(self.acquire_mbgrid(**mod_args))
        
        if mod_name == 'coastline':
            return(self.acquire_coastline(**mod_args))

        if mod_name == 'cudem':
            return(self.acquire_cudem(**mod_args))
            
    def acquire_from_config(self, config):
        def waffles_config(wc):
            print(wc)
            if wc['src_region'] is not None:
                wc['src_region'] = regions.Region().from_list(
                    wc['src_region']
                )
            return(wc)
        
        mod_name = list(config.keys())[0]
        args = waffles_config(config[mod_name])
        return(self._modules[mod_name]['class'](args))
    
## ==============================================
## Command-line Interface (CLI)
## $ waffles
##
## waffles cli
## ==============================================
_waffles_module_short_desc = lambda: ', '.join(
    ['{}'.format(key) for key in WaffleFactory()._modules])
_waffles_module_long_desc = lambda x: 'waffles modules:\n% waffles ... -M <mod>:key=val:key=val...\n\n  ' + '\n  '.join(
    ['\033[1m{:14}\033[0m{}\n'.format(key, x[key]['description']) for key in x]) + '\n'

waffles_cli_usage = """{cmd} ({wf_version}): Generate DEMs and derivatives.

usage: {cmd} [OPTIONS] DATALIST

Options:
  -R, --region\t\tSpecifies the desired REGION;
\t\t\tThis can either be a GMT-style region ( -R xmin/xmax/ymin/ymax )
\t\t\tor an OGR-compatible vector file with regional polygons. 
\t\t\tIf a vector file is supplied it will search each region found therein.
\t\t\tIf omitted, use the region gathered from the data in DATALIST.
  -E, --increment\tGridding CELL-SIZE in native units or GMT-style increments.
\t\t\tappend :<inc> to resample the output to the given <inc>: -E.3333333s/.3333333s:.1111111s/.1111111s
  -F, --format\t\tOutput grid FORMAT. [GTiff]
  -M, --module\t\tDesired Waffles MODULE and options. (see available Modules below)
\t\t\tsyntax is -M module:mod_opt=mod_val:mod_opt1=mod_val1:...
  -O, --output-name\tBASENAME for all outputs.
  -P, --epsg\t\tHorizontal projection of data as EPSG code [4326]
  -X, --extend\t\tNumber of cells with which to EXTEND the REGION.
\t\t\tappend :<num> to extend the processing region: -X6:12
  -T, --filter\t\tFILTER the output DEM using one or multiple filters. <fltr:fltr_val:split_value=z>
\t\t\tAvailable filters:
\t\t\t1: perform a Gaussian filter at -T1:<factor>.
\t\t\t2: use a Cosine Arch Filter at -T2:<dist(km)> search distance.
\t\t\t3: Spike Filter at -T3:<stand-dev. threshhold>.
\t\t\tThe -T switch may be set multiple times to perform multiple filters.
\t\t\tAppend :split_value=<num> to only filter values below z-value <num>.
\t\t\te.g. -T1:10:split_value=0 to smooth bathymetry (z<0) using Gaussian filter
  -C, --clip\t\tCLIP the output to the clip polygon -C<clip_ply.shp:invert=False>
  -G, --wg-config\tA waffles config JSON file. If supplied, will overwrite all other options.
\t\t\tgenerate a waffles_config JSON file using the --config flag.

  -p, --prefix\t\tSet BASENAME (-O) to PREFIX (append inc_nYYxYY_wXXxXX_YEAR info to output BASENAME).
  -r, --grid-node\tUse grid-node registration, default is pixel-node
  -w, --weights\t\tUse weights provided in the datalist to weight overlapping data.

  -a, --archive\t\tARCHIVE the datalist to the given region.
  -m, --mask\t\tGenerate a data MASK raster.
  -s, --spat\t\tGenerate SPATIAL-METADATA.
  -c, --continue\tDon't clobber existing files.
  -q, --quiet\t\tLower verbosity to a quiet.

  --help\t\tPrint the usage text
  --config\t\tSave the waffles config JSON and major datalist
  --modules\t\tDisply the module descriptions and usage
  --version\t\tPrint the version information

Datalists and data formats:
  A datalist is a file that contains a number of datalist entries, 
  while an entry is a space-delineated line:
  `path [format weight [name source date type resolution hdatum vdatum url]]`

Supported datalist formats: 
  {dl_formats}

Modules (see waffles --modules <module-name> for more info):
  {modules}

CIRES DEM home page: <http://ciresgroups.colorado.edu/coastalDEM>
""".format(cmd=os.path.basename(sys.argv[0]),
           dl_formats=dlim._datalist_fmts_short_desc(),
           modules=_waffles_module_short_desc(),
           wf_version=cudem.__version__)

def waffles_cli(argv = sys.argv):
    """run waffles from command-line

    See `waffles_cli_usage` for full cli options.
    """
    
    dls = []
    i_regions = []
    these_regions = []
    module = None
    wg_user = None
    want_prefix = False
    want_config = False
    status = 0
    i = 1
    wg = {}
    wg['verbose'] = True
    wg['xsample'] = None
    wg['ysample'] = None
    wg['fltr'] = []
    wg['name'] = 'waffles'
    
    while i < len(argv):
        arg = argv[i]
        if arg == '--region' or arg == '-R':
            i_regions.append(str(argv[i + 1]))
            i += 1
        elif arg[:2] == '-R': i_regions.append(str(arg[2:]))
        elif arg == '--module' or arg == '-M':
            module = str(argv[i + 1])
            i += 1
        elif arg[:2] == '-M': module = str(arg[2:])
        elif arg == '--increment' or arg == '-E':
            incs = argv[i + 1].split(':')
            xy_inc = incs[0].split('/')
            #wg['inc'] = utils.str2inc(incs[0])
            wg['xinc'] = utils.str2inc(xy_inc[0])
            if len(xy_inc) > 1:
                wg['yinc'] = utils.str2inc(xy_inc[1])
            else:
                wg['yinc'] = utils.str2inc(xy_inc[0])
            if len(incs) > 1:
                xy_samples = incs[1].split('/')
                wg['xsample'] = utils.str2inc(xy_samples[0])
                if len(xy_samples) > 1:
                    wg['ysample'] = utils.str2inc(xy_samples[1])
                else:
                    wg['ysample'] = utils.str2inc(xy_samples[0])
            i = i + 1
        elif arg[:2] == '-E':
            incs = arg[2:].split(':')
            xy_inc = incs[0].split('/')
            wg['xinc'] = utils.str2inc(xy_inc[0])
            if len(xy_inc) > 1:
                wg['yinc'] = utils.str2inc(xy_inc[1])
            else:
                wg['yinc'] = utils.str2inc(xy_inc[0])
            if len(incs) > 1:
                xy_samples = incs[1].split('/')
                wg['xsample'] = utils.str2inc(xy_samples[0])
                if len(xy_samples) > 1:
                    wg['ysample'] = utils.str2inc(xy_samples[1])
                else:
                    wg['ysample'] = utils.str2inc(xy_samples[0])

            #incs = arg[2:].split(':')
            #wg['inc'] = utils.str2inc(arg[2:].split(':')[0])
            #if len(incs) > 1: wg['sample'] = utils.str2inc(incs[1])
        elif arg == '--outname' or arg == '-O':
            wg['name'] = argv[i + 1]
            i += 1
        elif arg[:2] == '-O': wg['name'] = arg[2:]
        elif arg == '--format' or arg == '-F':
            wg['fmt'] = argv[i + 1]
            i += 1
        elif arg[:2] == '-F': wg['fmt'] = arg[2:]
        elif arg == '--filter' or arg == '-T':
            wg['fltr'].append(argv[i + 1])
            i += 1
        elif arg[:2] == '-T': wg['fltr'].append(arg[2:])
        elif arg == '--extend' or arg == '-X':
            exts = argv[i + 1].split(':')
            wg['extend'] = utils.int_or(exts[0], 0)
            if len(exts) > 1: wg['extend_proc'] = utils.int_or(exts[1], 10)
            i += 1
        elif arg[:2] == '-X':
            exts = arg[2:].split(':')
            wg['extend'] = utils.int_or(exts[0], 0)
            if len(exts) > 1: wg['extend_proc'] = utils.int_or(exts[1], 10)
        elif arg == '--wg-config' or arg == '-G':
            wg_user = argv[i + 1]
            i += 1
        elif arg[:2] == '-G': wg_user = arg[2:]
        elif arg == '--clip' or arg == '-C':
            wg['clip'] = argv[i + 1]
            i = i + 1
        elif arg[:2] == '-C': wg['clip'] = arg[2:]
        elif arg == '--chunk' or arg == '-K':
            wg['chunk'] = utils.int_or(argv[i + 1], None)
            i = i + 1
        elif arg[:2] == '-K': wg['chunk'] = utils.int_or(arg[2:], None)
        elif arg == '--epsg' or arg == '-P':
            wg['epsg'] = utils.int_or(argv[i + 1], 4326)
            i = i + 1
        elif arg[:2] == '-P': wg['epsg'] = utils.int_or(arg[2:], 4326)
        
        elif arg == '-w' or arg == '--weights':
            if 'weights' not in wg.keys():
                wg['weights'] = 1
            else:
                wg['weights'] += 1
        elif arg == '-t' or arg == '--threads': want_threads = True
        elif arg == '-p' or arg == '--prefix': want_prefix = True
        elif arg == '-a' or arg == '--archive': wg['archive'] = True
        elif arg == '-m' or arg == '--mask': wg['mask'] = True
        elif arg == '-s' or arg == '--spat': wg['spat'] = True
        elif arg == '-c' or arg == '--continue': wg['clobber'] = False
        elif arg == '-r' or arg == '--grid-node': wg['node'] = 'grid'

        elif arg == '--quiet' or arg == '-q': wg['verbose'] = False
        elif arg == '--config': want_config = True
        elif arg == '--modules' or arg == '-m':
            try:
                if argv[i + 1] in WaffleFactory()._modules.keys():
                    sys.stderr.write(
                        _waffles_module_long_desc({k: WaffleFactory()._modules[k] for k in (argv[i + 1],)})
                    )
                else: sys.stderr.write(
                        _waffles_module_long_desc(WaffleFactory()._modules)
                )
            except: sys.stderr.write(
                    _waffles_module_long_desc(WaffleFactory()._modules)
            )
            sys.exit(0)
        elif arg == '--help' or arg == '-h':
            sys.stderr.write(waffles_cli_usage)
            sys.exit(0)
        elif arg == '--version' or arg == '-v':
            sys.stdout.write('{}\n'.format(cudem.__version__))
            sys.exit(0)
        elif arg[0] == '-':
            sys.stderr.write(waffles_cli_usage)
            utils.echo_error_msg('{} is not a valid waffles cli switch'.format(arg))
            sys.exit(0)
        else: dls.append(arg)
        i += 1

    ## ==============================================
    ## load the user wg json and run waffles with that.
    ## ==============================================
    if wg_user is not None:
        if os.path.exists(wg_user):
            #try:
            with open(wg_user, 'r') as wgj:
                wg = json.load(wgj)

                for key in wg.keys():
                    #this_waffle = WaffleFactory()._modules[key]['class'](wg[key])
                    this_waffle = WaffleFactory().acquire_from_config(wg)
                    this_waffle.generate()

                sys.exit(0)
            # except Exception as e:
            #     utils.echo_error_msg(e)
            #     sys.exit(-1)
        else:
            utils.echo_error_msg(
                'specified waffles config file does not exist, {}'.format(wg_user)
            )
            sys.stderr.write(waffles_cli_usage)
            sys.exit(-1)

    ## ==============================================
    ## Otherwise run from cli options...
    ## set the dem module
    ## ==============================================
        
    if module is None:
        sys.stderr.write(waffles_cli_usage)
        utils.echo_error_msg(
            '''must specify a waffles -M module.'''
        )
        sys.exit(-1)
                
    # if WaffleFactory()._modules[mod]['datalist-p']:
    #     if len(dls) == 0:
    #         sys.stderr.write(waffles_cli_usage)
    #         utils.echo_error_msg('''must specify a datalist/entry, try gmrt or srtm for global data.''')
    #         sys.exit(-1)

    ## ==============================================
    ## check the increment
    ## ==============================================
    if 'xinc' in wg.keys():
        if wg['xinc'] is None:
            sys.stderr.write(waffles_cli_usage)
            utils.echo_error_msg('''must specify a gridding increment.''')
            sys.exit(-1)
        else:
            if wg['yinc'] is None:
                wg['yinc'] = wg['xinc']
    else:
        sys.stderr.write(waffles_cli_usage)
        utils.echo_error_msg('''must specify a gridding increment.''')
        sys.exit(-1)      
    
    ## ==============================================
    ## set the datalists and names
    ## ==============================================
    wg['data'] = dls

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
        these_regions = []
        utils.echo_error_msg('failed to parse region(s), {}'.format(i_regions))
    else:
        if wg['verbose']: utils.echo_msg('parsed {} region(s)'.format(len(these_regions)))

    name = wg['name']
    for i, this_region in enumerate(these_regions):
        wg['src_region'] = this_region
        if want_prefix or len(these_regions) > 1:
            wg['name'] = utils.append_fn(
                name, wg['src_region'], wg['xsample'] if wg['xsample'] is not None else wg['xinc']
            )
        this_waffle = WaffleFactory(mod=module, **wg).acquire()
        if want_config:
            this_wg = this_waffle._config
            utils.echo_msg(json.dumps(this_wg, indent = 4, sort_keys = True))
            with open('{}.json'.format(this_waffle.name), 'w') as wg_json:
                utils.echo_msg('generating waffles config file: {}.json'.format(this_waffle.name))
                wg_json.write(json.dumps(this_wg, indent = 4, sort_keys = True))
        else:
            this_waffle.generate()
            #utils.echo_msg('generated DEM: {} @ {}/{}'.format(wf.fn, wf.))
### End
