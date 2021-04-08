### waffles.py
##
## Copyright (c) 2010 - 2021 CIRES Coastal DEM Team
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
import numpy as np
import math
import json
import gdal
import ogr
import osr
from cudem import dlim
from cudem import regions
from cudem import utils
from cudem import xyzfun
from cudem import demfun
from cudem import metadata
from cudem import vdatumfun

__version__ = '0.10.0'

## ==============================================
## Waffles XYZ blocking
## ==============================================        
def xyz_block(src_xyz, src_region, inc, weights=False, verbose=False):
    """block the src_xyz data to the mean block value

    Args:
      src_xyz (generataor): list/generator of xyz data
      region (list): a `region` list [xmin, xmax, ymin, ymax]
      inc (float): blocking increment, in native units
      weights (bool): block using weights
      verbose (bool): increase verbosity    

    Yields:
      list: xyz data for each block with data
    """

    xcount, ycount, dst_gt = src_region.geo_transform(x_inc=inc)
    sumArray = np.zeros((ycount, xcount))
    gdt = gdal.GDT_Float32
    ptArray = np.zeros((ycount, xcount))
    if weights: wtArray = np.zeros((ycount, xcount))
    if verbose: utils.echo_msg('blocking data to {}/{} grid'.format(ycount, xcount))
    for this_xyz in src_xyz:
        if weights: this_xyz.z = this_xyz.z * this_xyz.w
        if regions.xyz_in_region_p(this_xyz, src_region):
            xpos, ypos = utils._geo2pixel(this_xyz.x, this_xyz.y, dst_gt)
            try:
                sumArray[ypos, xpos] += this_xyz.z
                ptArray[ypos, xpos] += 1
                if weights: wtArray[ypos, xpos] += this_xyz.w
            except: pass
            
    ptArray[ptArray == 0] = np.nan
    if weights:
        wtArray[wtArray == 0] = 1
        outarray = (sumArray / wtArray) / ptArray
    else: outarray = sumArray / ptArray

    sumArray = ptArray = None
    if weights: wtArray = None

    outarray[np.isnan(outarray)] = -9999
    
    for y in range(0, ycount):
        for x in range(0, xcount):
            geo_x, geo_y = utils._pixel2geo(x, y, dst_gt)
            out_xyz = xyzfun.XYZPoint(x=geo_x, y=geo_y, z=outarray[y,x])
            if out_xyz.z != -9999:
                yield(out_xyz)

def xyz2gdal(src_xyz, dst_gdal, src_region, inc, dst_format = 'GTiff', mode = 'n', epsg = 4326, verbose = False):
    '''Create a GDAL supported grid from xyz data 
    `mode` of `n` generates a num grid
    `mode` of `m` generates a mean grid
    `mode` of `k` generates a mask grid
    `mode` of 'w' generates a wet/dry mask grid

    returns output, status'''

    xcount, ycount, dst_gt = src_region.geo_transform(x_inc=inc)
    if verbose:
        progress = utils.CliProgress('generating uninterpolated num grid {} @ {}/{}'.format(mode, ycount, xcount))
    if mode == 'm' or mode == 'w':
        sumArray = np.zeros((ycount, xcount))
    gdt = gdal.GDT_Float32
    ptArray = np.zeros((ycount, xcount))
    ds_config = demfun.set_infos(xcount, ycount, xcount * ycount, dst_gt, utils.sr_wkt(epsg), gdt, -9999, dst_format)
    for this_xyz in src_xyz:
        if regions.xyz_in_region_p(this_xyz, src_region):
            xpos, ypos = utils._geo2pixel(this_xyz.x, this_xyz.y, dst_gt)
            try:
                if mode == 'm' or mode == 'w':
                    sumArray[ypos, xpos] += this_xyz.z
                if mode == 'n' or mode == 'm':
                    ptArray[ypos, xpos] += 1
                else: ptArray[ypos, xpos] = 1
            except Exception as e:
                pass
    if mode == 'm' or mode == 'w':
        ptArray[ptArray == 0] = np.nan
        outarray = sumArray / ptArray
        if mode == 'w':
            outarray[outarray >= 0] = 1
            outarray[outarray < 0] = 0
    elif mode == 'n': outarray = ptArray
    else: outarray = ptArray
    outarray[np.isnan(outarray)] = -9999
    if verbose: progress.end(0, 'generated uninterpolated num grid {} @ {}/{}'.format(mode, ycount, xcount))
    return(utils.gdal_write(outarray, dst_gdal, ds_config))

def gdal_ogr_mask_union(src_layer, src_field, dst_defn = None):
    '''`union` a `src_layer`'s features based on `src_field` where
    `src_field` holds a value of 0 or 1. optionally, specify
    an output layer defn for the unioned feature.

    returns the output feature class'''
    
    if dst_defn is None: dst_defn = src_layer.GetLayerDefn()
    multi = ogr.Geometry(ogr.wkbMultiPolygon)
    feats = len(src_layer)
    utils.echo_msg('unioning {} features'.format(feats))
    for n, f in enumerate(src_layer):
        gdal.TermProgress_nocb((n+1 / feats) * 100)
        if f.GetField(src_field) == 0:
            src_layer.DeleteFeature(f.GetFID())
        elif f.GetField(src_field) == 1:
            f.geometry().CloseRings()
            wkt = f.geometry().ExportToWkt()
            multi.AddGeometryDirectly(ogr.CreateGeometryFromWkt(wkt))
            src_layer.DeleteFeature(f.GetFID())
    #union = multi.UnionCascaded() ## slow on large multi...
    out_feat = ogr.Feature(dst_defn)
    out_feat.SetGeometryDirectly(multi)
    #union = multi = None
    return(out_feat)

def ogr_clip(src_ogr, dst_ogr, clip_region = None, dn = "ESRI Shapefile"):
    driver = ogr.GetDriverByName(dn)
    ds = driver.Open(src_ogr, 0)
    layer = ds.GetLayer()

    clip_region.export_as_ogr('tmp_clip.shp')
    c_ds = driver.Open('tmp_clip.shp', 0)
    c_layer = c_ds.GetLayer()
    
    dst_ds = driver.CreateDataSource(dst_ogr)
    dst_layer = dst_ds.CreateLayer(dst_ogr.split('.')[0], geom_type=ogr.wkbMultiPolygon)

    layer.Clip(c_layer, dst_layer)

    ds = c_ds = dst_ds = None

def ogr_empty_p(src_ogr):
    driver = ogr.GetDriverByName('ESRI Shapefile')
    ds = driver.Open(src_ogr, 0)

    if ds is not None:
        layer = ds.GetLayer()
        fc = layer.GetFeatureCount()
        if fc == 0:
            return(True)
        else: return(False)
    else: return(True)

class Waffled:
    """Providing a waffled object."""
    
    def __init__(self, fn=None, src_region=None, epsg=None, warp=None, fmt=None, clip_poly=None, cudem=False):
        self.fn = fn
        self.region = src_region
        self.epsg = epsg
        self.warp = warp
        self.fmt = fmt
        self.cudem = cudem
        self.clip_poly = clip_poly
        
    def process(self):
        pass

class WaffledRaster(Waffled):
    """Providing a waffled RASTER file processor."""
    
    def __init__(self, sample_inc=None, filters=[], node='pixel', **kwargs):
        super().__init__(**kwargs)
        self.src_ds = None
        self.ds_open_p = False
        self.node = node
        self.sample_inc = sample_inc
        self.filters = filters
            
    def valid_p(self):
        if not os.path.exists(self.fn): return(False)
        gdi = demfun.infos(this_dem, scan=True)
        
        if gdi is not None:
            if np.isnan(gdi['zr'][0]):
                return(False)
        else:
            return(False)
        return(True)
            
    def process(self):

        demfun.set_nodata(self.fn, nodata=-9999, convert_array=True)
        
        if len(self.filters) > 0:
            for f in self.filters:
                fltr_val = None
                split_val = None
                fltr_opts = f.split(':')
                fltr = fltr_opts[0]
                if len(fltr_opts) > 1:
                    fltr_val = fltr_opts[1]
                if len(fltr_opts) > 2:
                    split_val= fltr_opts[2]

                if demfun.filter_(self.fn, '__tmp_fltr.tif', fltr=fltr, fltr_val=fltr_val) == 0:
                    os.rename('__tmp_fltr.tif', self.fn)
            
        if self.sample_inc is not None:
            if demfun.sample(self.fn, '__tmp_sample.tif', self.sample_inc, self.region)[1] == 0:
                os.rename('__tmp_sample.tif', self.fn)
            
        if self.clip_poly is not None:
            if demfun.clip(self.fn, src_ply=self.clip_poly)[1] == 0:
                os.rename('__tmp_cut__.tif', '{}'.format(self.fn))
                
        if demfun.cut(self.fn, self.region, '__tmp_cut__.tif')[1] == 0:
            os.rename('__tmp_cut__.tif', '{}'.format(self.fn))
            
        demfun.set_metadata(self.fn, node=self.node, cudem=self.cudem)
        demfun.set_epsg(self.fn, self.epsg)
        
        if self.fmt != 'GTiff':
            out_dem = utils.gdal2gdal(self.fn, dst_fmt=self.fmt)
            if out_dem is not None:
                utils.remove_glob(self.fn)
            
class WaffledVector(Waffled):
    """Providing a waffled VECTOR file processor."""
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        pass

    def set_metadata(self):
        pass
    
    def process(self):
        pass

class Waffle:
    
    def __init__(self, data=[], src_region=None, inc=None, name='waffles_dem',
                 node='pixel', fmt='GTiff', extend=0, extend_proc=20, weights=None,
                 fltr=[], sample=None, clip=None, chunk=None, epsg=4326,
                 verbose=False, archive=False, mask=False, spat=False, clobber=True):

        self.data = data
        self.region = src_region
        self.inc = inc
        self.name = name
        self.node = node
        self.fmt = fmt
        self.extend = extend
        self.extend_proc = extend_proc
        self.weights = weights
        self.fltr = fltr
        self.sample = sample
        self.clip = clip
        self.chunk = chunk
        self.epsg = utils.int_or(epsg)
        self.mod = None
        self.archive = archive
        self.mask = mask
        self.clobber = clobber
        self.verbose = verbose
        self.gc = utils.config_check()
        self.spat = spat
        
        if self.node == 'grid':
            self.region = self.region.buffer(self.inc*.5)
        self.p_region = self.proc_region()
        self.d_region = self.dist_region()
        
        self.data = [dlim.DatasetFactory(
            fn=" ".join(['-' if x == "" else x for x in dl.split(":")]),
            src_region=self.p_region, verbose=self.verbose, weight=self.weights,
            epsg=self.epsg).acquire_dataset() for dl in self.data]

        self.data = [d for d in self.data if d is not None]
        
        for d in self.data:
            d.parse()

        self.waffle_out = {}
        self.waffled = False
        out, status = utils.run_cmd('gmt gmtset IO_COL_SEPARATOR = SPACE', verbose = False)
        
    def proc_region(self):
        
        pr = regions.Region().from_region(self.region)
        return(pr.buffer((self.inc*self.extend_proc) + (self.inc*self.extend)))

    def dist_region(self):
            
        dr = regions.Region().from_region(self.region)
        return(dr.buffer((self.inc*self.extend)))

    def run(self):
        pass

    def generate(self):

        self.run()
        if self.mask:
            WaffledRaster(fn='{}_m.tif'.format(self.name), epsg=self.epsg, fmt=self.fmt, sample_inc=self.sample, src_region = self.region).process()
        self.waffled = True

    def xyz_archive(self, src_xyz):
        ## update to new methods...!!
        if self.region is None:
            a_name = self.name
        else: a_name = '{}_{}_{}'.format(self.name, self.region.format('fn'), utils.this_year())
        i_dir = os.path.dirname(entry[0])
        i_xyz = os.path.basename(entry[0]).split('.')[0]
        i_xyz = ''.join(x for x in i_xyz if x.isalnum())
        a_dir = os.path.join(dirname, a_name, 'data', entry[-1])
        a_xyz_dir = os.path.join(a_dir, 'xyz')
        a_xyz = os.path.join(a_xyz_dir, i_xyz + '.xyz')
        a_dl = os.path.join(a_xyz_dir, '{}.datalist'.format(entry[-1]))
        if not os.path.exists(a_dir):
            os.makedirs(a_dir)
        if not os.path.exists(a_xyz_dir):
            os.makedirs(a_xyz_dir)
        with open(a_xyz, 'w') as fob:
            for xyz in src_xyz:
                xyz.dump(fob)
                yield(xyz)
        inf = utils.mb_inf(a_xyz)
        datalist_append_entry([i_xyz + '.xyz', 168, entry[2] if entry[2] is not None else 1], a_dl)
    
    def xyz_mask(self, src_xyz, dst_gdal, dst_format='GTiff'):
        """Create a num grid mask of xyz data. The output grid
        will contain 1 where data exists and 0 where no data exists.

        yields the xyz data
        """

        xcount, ycount, dst_gt = self.region.geo_transform(x_inc=self.inc)
        ptArray = np.zeros((ycount, xcount))
        ds_config = demfun.set_infos(
            xcount, ycount, (xcount*ycount), dst_gt, utils.sr_wkt(self.epsg),
            gdal.GDT_Float32, -9999, 'GTiff')
        for this_xyz in src_xyz:
            yield(this_xyz)
            if regions.xyz_in_region_p(this_xyz, self.region):
                xpos, ypos = utils._geo2pixel(this_xyz.x, this_xyz.y, dst_gt)
                try:
                    ptArray[ypos, xpos] = 1
                except:
                    pass
        out, status = utils.gdal_write(ptArray, dst_gdal, ds_config)    

    def spat_meta(self, yxyz=False, **kwargs):
        self.dst_layer = '{}_sm'.format(self.name)
        self.dst_vector = self.dst_layer + '.shp'
        self.v_fields = ['Name', 'Agency', 'Date', 'Type', 'Resolution', 'HDatum', 'VDatum', 'URL']
        self.t_fields = [ogr.OFTString, ogr.OFTString, ogr.OFTString, ogr.OFTString,
                         ogr.OFTString, ogr.OFTString, ogr.OFTString, ogr.OFTString]
        utils.remove_glob('{}.*'.format(self.dst_layer))
        utils.gdal_prj_file('{}.prj'.format(self.dst_layer), self.epsg)
        sm_inc = self.inc if self.inc >= utils.str2inc('.3333333s') else utils.str2inc('.3333333s')
        
        self.ds = ogr.GetDriverByName('ESRI Shapefile').CreateDataSource(self.dst_vector)
        if self.ds is not None: 
            self.layer = self.ds.CreateLayer('{}'.format(self.dst_layer), None, ogr.wkbMultiPolygon)
            [self.layer.CreateField(ogr.FieldDefn('{}'.format(f), self.t_fields[i])) for i, f in enumerate(self.v_fields)]
            [self.layer.SetFeature(feature) for feature in self.layer]
        else: self.layer = None

        for xdl in self.data:
            #xdl.parse()
            for x in xdl.data_lists.keys():
                xdl.data_entries = xdl.data_lists[x]
                dl_name = x
                o_v_fields = [dl_name, 'Unknown', '0', 'xyz_elevation', 'Unknown', 'WGS84', 'NAVD88', 'URL']
                defn = None if self.layer is None else self.layer.GetLayerDefn()
                for xyz in xdl.mask_xyz('{}.tif'.format(dl_name), sm_inc):
                    if yxyz:
                        yield(xyz)
                    else:
                        pass

                if demfun.infos('{}.tif'.format(dl_name), scan=True)['zr'][1] == 1:
                    tmp_ds = ogr.GetDriverByName('ESRI Shapefile').CreateDataSource('{}_poly.shp'.format(dl_name))
                    if tmp_ds is not None:
                        tmp_layer = tmp_ds.CreateLayer('{}_poly'.format(dl_name), None, ogr.wkbMultiPolygon)
                        tmp_layer.CreateField(ogr.FieldDefn('DN', ogr.OFTInteger))
                        demfun.polygonize('{}.tif'.format(dl_name), tmp_layer, verbose=self.verbose)

                        if len(tmp_layer) > 1:
                            if defn is None: defn = tmp_layer.GetLayerDefn()
                            out_feat = gdal_ogr_mask_union(tmp_layer, 'DN', defn)
                            [out_feat.SetField(f, o_v_fields[i]) for i, f in enumerate(self.v_fields)]
                            self.layer.CreateFeature(out_feat)
                    tmp_ds = None
                    utils.remove_glob('{}_poly.*'.format(dl_name), '{}.tif'.format(dl_name))
        self.ds = None
        
        dr = regions.Region().from_region(self.region)
        dr.buffer((sm_inc*self.extend))
        #print(dr.format('gmt'))
        #ogr_clip('{}'.format(self.dst_vector), '__tmp_clip.shp', clip_region = dr)
        utils.run_cmd('ogrinfo -spat {} -dialect SQLITE -sql "UPDATE {} SET geometry = ST_MakeValid(geometry)" {}\
        '.format(dr.format('ul_lr'), self.dst_layer, self.dst_vector))
        utils.run_cmd('ogr2ogr -clipsrc {} __tmp_clip.shp {} -overwrite -nlt POLYGON -skipfailures'.format(dr.format('ul_lr'), self.dst_vector), verbose=True)
        utils.run_cmd('ogr2ogr {} __tmp_clip.shp -overwrite'.format(self.dst_vector), verbose=True)
        utils.remove_glob('__tmp_clip.*')
        utils.run_cmd('ogrinfo -spat {} -dialect SQLITE -sql "UPDATE {} SET geometry = ST_MakeValid(geometry)" {}\
        '.format(dr.format('ul_lr'), self.dst_layer, self.dst_vector))
                 
    def yield_xyz(self, **kwargs):
        """yields the xyz data"""
        
        for xdl in self.data:
            if self.spat:
                xyz_yield = metadata.SpatialMetadata(
                    data=self.data, src_regtion=self.p_region, inc=self.inc, extend=self.extend, epsg=self.epsg,
                    node=self.node, name=self.name, verbose=self.verbose).yield_xyz()
            else:
                xyz_yield = xdl.yield_xyz()
            if self.mask:
                xyz_yield = self.xyz_mask(xyz_yield, '{}_m.tif'.format(self.name))
            if self.archive:
                xyz_yield = self.xyz_archive(xyz_yield)
        
            for xyz in xyz_yield:
                yield(xyz)

        if self.spat:
            dst_layer = '{}_sm'.format(self.name)
            dst_vector = dst_layer + '.shp'
            utils.run_cmd('ogr2ogr -clipsrc {} __tmp_clip.shp {} -overwrite -nlt POLYGON -skipfailures'.format(dr.format('ul_lr'), dst_vector), verbose=True)
            utils.run_cmd('ogr2ogr {} __tmp_clip.shp -overwrite'.format(dst_vector), verbose=True)
            utils.remove_glob('__tmp_clip.*')
                
        if self.archive:
            a_dl = os.path.join('archive', '{}.datalist'.format(self.name))
            
            for dir_, _, files in os.walk('archive'):
                for f in files:
                    if '.datalist' in f:
                        rel_dir = os.path.relpath(dir_, 'archive')
                        rel_file = os.path.join(rel_dir, f)
                        datalists.datalist_append_entry([rel_file, -1, 1], a_dl)
                    
    def dump_xyz(self, dst_port=sys.stdout, encode=False, **kwargs):
        for xyz in self.yield_xyz(**kwargs):
            xyz.dump(include_w = True if self.weights is not None else False, dst_port=dst_port, encode=encode, **kwargs)
            
## ==============================================
## Waffles GMT surface module
## ==============================================
class GMTSurface(Waffle):
    
    def __init__(self, tension = .35, relaxation = 1.2,
                 lower_limit = 'd', upper_limit = 'd', **kwargs):
        """generate a DEM with GMT surface"""

        super().__init__(**kwargs)
        
        if self.gc['GMT'] is None:
            utils.echo_error_msg('GMT must be installed to use the SURFACE module')
            return(None, -1)
        
        self.tension = tension
        self.relaxation = relaxation
        self.lower_limit = lower_limit
        self.upper_limit = upper_limit

        self.mod = 'surface'
        
    def run(self):
        
        dem_surf_cmd = ('gmt blockmean {} -I{:.10f}{} -V -r | gmt surface -V {} -I{:.10f} -G{}.tif=gd+n-9999:GTiff -T{} -Z{} -Ll{} -Lu{} -r\
        '.format(self.p_region.format('gmt'), self.inc, ' -Wi' if self.weights else '', self.p_region.format('gmt'),
                 self.inc, self.name, self.tension, self.relaxation, self.lower_limit, self.upper_limit))
        out, status = utils.run_cmd(dem_surf_cmd, verbose = self.verbose, data_fun = lambda p: self.dump_xyz(dst_port=p, encode=True))
        
        WaffledRaster(fn='{}.tif'.format(self.name), epsg=self.epsg, fmt=self.fmt, sample_inc=self.sample, src_region=self.d_region, filters=self.fltr).process()
        
## ==============================================
## Waffles GMT triangulate module
## ==============================================
class GMTTriangulate(Waffle):
    
    def __init__(self, **kwargs):
        """generate a DEM with GMT triangulate"""

        super().__init__(**kwargs)
        
        if self.gc['GMT'] is None:
            utils.echo_error_msg('GMT must be installed to use the TRIANGULATE module')
            return(None, -1)
        
        self.mod = 'triangulate'
        
    def run(self):
        
        #dem_tri_cmd = ('gmt blockmean {} -I{:.10f}{} -V -r | gmt triangulate -V {} -I{:.10f} -G{}.tif=gd:GTiff -r\
        #'.format(self.p_region.format('gmt'), self.inc, ' -Wi' if self.weights else '', self.p_region.format('gmt'), self.inc, self.name))
        dem_tri_cmd = ('gmt triangulate -V {} -I{:.10f} -G{}.tif=gd:GTiff -r'.format(self.p_region.format('gmt'), self.inc, self.name))
        out, status = utils.run_cmd(dem_tri_cmd, verbose = self.verbose, data_fun = lambda p: self.dump_xyz(dst_port=p, encode=True))

        WaffledRaster(fn='{}.tif'.format(self.name), epsg=self.epsg, fmt=self.fmt, sample_inc=self.sample, src_region = self.d_region, filters=self.fltr).process()

## ==============================================
## Waffles 'NUM grid' module
## Uninterpolated grdi from data;
## num methods include: mask, mean, num, landmask and any gmt grd2xyz -A option.
## ==============================================
class WafflesNum(Waffle):
    
    def __init__(self, mode='n', **kwargs):
        """generate a DEM with GMT triangulate"""

        super().__init__(**kwargs)
        
        if self.gc['GMT'] is None:
            utils.echo_error_msg('GMT must be installed to use the TRIANGULATE module')
            return(None, -1)
        
        self.mod = 'num'
        self.mode = mode
        
    def run(self):

        if self.mode.startswith('A'):
            if self.gc['GMT'] is None:
                utils.echo_error_msg('GMT must be installed to use the Mode `A` with the NUM module')
                return(None, -1)

            dem_xyz2grd_cmd = ('gmt xyz2grd -{} -V {} -I{:.10f} -G{}.tif=gd:GTiff -r{}'.format(
                self.mode, self.p_region.format('gmt'), self.inc, self.name, ' -Wi' if self.weights else '',))
            out, status = utils.run_cmd(dem_xyz2grd_cmd, verbose=self.verbose, data_fun=lambda p: self.dump_xyz(dst_port=p, encode=True))
        else:
            dly = self.yield_xyz()
            if self.weights: dly = xyz_block(dly, self.p_region, self.inc, weights=True)
            out, status = xyz2gdal(dly, '{}.tif'.format(self.name), self.p_region, self.inc, dst_format=self.fmt, mode=self.mode, verbose=self.verbose)

        WaffledRaster(fn='{}.tif'.format(self.name), epsg=self.epsg, fmt=self.fmt, sample_inc=self.sample, src_region=self.d_region).process()        

## ==============================================
## Waffles 'Vdatum conversion grid' module
## vertical datum transformation grid;
## note: U.S. Only
## ==============================================
class WafflesVdatum(Waffle):
    
    def __init__(self, ivert='navd88', overt='mhw', region='3', jar=None, **kwargs):
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
        self.vdc = vdatumfun.Vdatum(ivert=ivert, overt=overt, region=region, jar=jar)

    def create_null(self, outfile, extent, cellsize, nodata, outformat, verbose, overwrite):
        """create a nodata grid"""
        
        xcount, ycount, gt = extent.geo_transform(x_inc = cellsize)
        ds_config = demfun.set_infos(xcount, ycount, xcount * ycount, gt, utils.sr_wkt(4326), gdal.GDT_Float32, nodata, outformat)
        nullArray = np.zeros( (ycount, xcount) )
        nullArray[nullArray==0]=nodata
        utils.gdal_write(nullArray, outfile, ds_config)
        
    def run(self):

        self.create_null('empty.tif', self.p_region, 0.00083333, 0, 'GTiff', self.verbose, True)
        demfun.set_nodata('empty.tif')
        empty = dlim.RasterFile(fn='empty.tif', data_format=200, src_region=self.region, epsg=self.epsg,
                                name=self.name, verbose=self.verbose)

        empty.parse()
        with open('empty.xyz', 'w') as mt_xyz:
            empty.dump_xyz(mt_xyz)
        self.vdc.run_vdatum('empty.xyz')

        if os.path.exists('result/empty.xyz') and os.stat('result/empty.xyz').st_size != 0:
            empty_xyz = dlim.XYZFile(fn='result/empty.xyz', data_format=168, src_region=self.region, epsg=self.epsg,
                                     name=self.name, verbose=self.verbose)
            empty_infos = empty_xyz.inf()
            
            ll = 'd' if empty_infos['minmax'][4] < 0 else '0'
            lu = 'd' if empty_infos['minmax'][5] > 0 else '0'

            GMTSurface(data=['result/empty.xyz'], name=self.name, src_region=self.region, inc=self.inc, epsg=self.epsg, verbose=self.verbose, tension=0, upper_limit=lu, lower_limit=ll).generate()
        else:
            utils.echo_error_msg('failed to generate VDatum grid, check settings')
            vd_out = {}
            status = -1

        utils.remove_glob('empty.*', 'result/*', 'result')        
        WaffledRaster(fn='{}.tif'.format(self.name), epsg=self.epsg, fmt=self.fmt, sample_inc=self.sample, src_region=self.d_region).process()        

## ==============================================
## Waffles GDAL_GRID module
## see gdal_grid for more info and gridding algorithms
## ==============================================
class WafflesGDALGrid(Waffle):
    
    def __init__(self, **kwargs):
        """run gdal grid using alg_str

        parse the data through xyzfun.xyz_block to get weighted mean before
        building the GDAL dataset to pass into gdal_grid
        
        Args: 
          alg_str (str): the gdal_grid algorithm string
        """
        
        super().__init__(**kwargs)
        
        self.alg_str = 'linear:radius=-1'
        self.mod = self.alg_str.split(':')[0]
        
    def run(self):
        
        _prog = utils.CliProgress('running GDAL GRID {} algorithm @ {}...\
        '.format(self.alg_str.split(':')[0], self.p_region.format('fn')))
        _prog_update = lambda x, y, z: _prog.update()
        dly = xyz_block(self.yield_xyz(), self.p_region, self.inc, weights = True if self.weights is not None else False, verbose=self.verbose)
        ds = xyzfun.xyz2gdal_ds(dly, '{}'.format(self.name))
        
        if ds.GetLayer().GetFeatureCount() == 0:
            utils.echo_error_msg('no input data')
            
        xcount, ycount, dst_gt = self.p_region.geo_transform(x_inc=self.inc)
        
        gd_opts = gdal.GridOptions(outputType = gdal.GDT_Float32, noData = -9999, format = 'GTiff', width = xcount,
                                   height = ycount, algorithm = self.alg_str, callback = _prog_update if self.verbose else None,
                                   outputBounds = [self.p_region.xmin, self.p_region.ymax, self.p_region.xmax, self.p_region.ymin])
        
        gdal.Grid('{}.tif'.format(self.name), ds, options = gd_opts)
        ds = None
        demfun.set_nodata('{}.tif'.format(self.name, nodata=-9999, convert_array=False))
        _prog.end(0, 'ran GDAL GRID {} algorithm @ {}.'.format(self.alg_str.split(':')[0], self.p_region.format('fn')))
        WaffledRaster(fn='{}.tif'.format(self.name), epsg=self.epsg, fmt=self.fmt, sample_inc=self.sample, src_region=self.d_region, filters=self.fltr).process()

class WafflesLinear(WafflesGDALGrid):

    def __init__(self, radius=None, nodata=-9999, **kwargs):
        super().__init__(**kwargs)
        
        radius = self.inc * 4 if radius is None else utils.str2inc(radius)
        self.alg_str = 'linear:radius={}:nodata={}'.format(radius, nodata)

class WafflesInvDst(WafflesGDALGrid):

    def __init__(self, power = 2.0, smoothing = 0.0, radius1 = None, radius2 = None, angle = 0.0,
                   max_points = 0, min_points = 0, nodata = -9999, **kwargs):
        super().__init__(**kwargs)

        radius1 = self.inc * 2 if radius1 is None else utils.str2inc(radius1)
        radius2 = self.inc * 2 if radius2 is None else utils.str2inc(radius2)
        self.alg_str = 'invdist:power={}:smoothing={}:radius1={}:radius2={}:angle={}:max_points={}:min_points={}:nodata={}'\
            .format(power, smoothing, radius1, radius2, angle, max_points, min_points, nodata)

class WafflesMovingAverage(WafflesGDALGrid):

    def __init__(self, radius1=None, radius2=None, angle=0.0, min_points=0, nodata=-9999, **kwargs):
        super().__init__(**kwargs)

        radius1 = self.inc * 2 if radius1 is None else utils.str2inc(radius1)
        radius2 = self.inc * 2 if radius2 is None else utils.str2inc(radius2)
        self.alg_str = 'average:radius1={}:radius2={}:angle={}:min_points={}:nodata={}'\
            .format(radius1, radius2, angle, min_points, nodata)
        
class WafflesNearest(WafflesGDALGrid):

    def __init__(self, radius1=None, radius2=None, angle=0.0, nodata=-9999, **kwargs):
        super().__init__(**kwargs)

        radius1 = self.inc * 2 if radius1 is None else utils.str2inc(radius1)
        radius2 = self.inc * 2 if radius2 is None else utils.str2inc(radius2)
        self.alg_str = 'nearest:radius1={}:radius2={}:angle={}:nodata={}'\
            .format(radius1, radius2, angle, min_points, nodata)
        
        
class WaffleFactory:

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
        },
        'triangulate': {
            'name': 'triangulate',
            'datalist-p': True,
            'description': '''TRIANGULATION DEM via GMT triangulate\n
Generate a DEM using GMT's triangulate command.
        
< triangulate >''',
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
        },
        'linear': {
            'name': 'linear',
            'datalist-p': True,
            'description': """LINEAR DEM via gdal_grid\n
Generate a DEM using GDAL's gdal_grid command.
see gdal_grid --help for more info

< linear:radius=-1 >
 :radius=[val] - search radius""",
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
        },
        'average': {
            'name': 'average',
            'datalist-p': True,
            'description': """Moving AVERAGE DEM via gdal_grid\n
Generate a DEM using GDAL's gdal_grid command.
see gdal_grid --help for more info

< nearest:radius1=0:radius2=0:angle=0:min_points=0:nodata=0 >
 :radius1=[val] - search radius
 :radius2=[val] - search radius""",
        },
        'invdst': {
            'name': 'invdst',
            'datalist-p': True,
            'description': """INVERSE DISTANCE DEM via gdal_grid\n
Generate a DEM using GDAL's gdal_grid command.
see gdal_grid --help for more info

< nearest:power=2.0:smoothing=0.0:radius1=0:radius2=0:angle=0:max_points=0:min_points=0:nodata=0 >
 :radius1=[val] - search radius
 :radius2=[val] - search radius""",
        },
        'vdatum': {
            'name': 'vdatum',
            'datalist-p': False,
            'description': """VDATUM conversion grid.
            """
        },
    }
    
    def __init__(self, data=[], src_region=None, inc=None, name='waffles_dem',
                 node='pixel', fmt='GTiff', extend=0, extend_proc=20, weights=None,
                 fltr=[], sample = None, clip=None, chunk=None, epsg=4326,
                 verbose = False, archive=False, mask=False, spat=False, clobber=True):

        self.data = data
        self.region = src_region
        self.inc = inc
        self.name = name
        self.node = node
        self.fmt = fmt
        self.extend = extend
        self.extend_proc = extend_proc
        self.weights = weights
        self.fltr = fltr
        self.sample = sample
        self.clip = clip
        self.chunk = chunk
        self.epsg = utils.int_or(epsg)
        self.archive = archive
        self.mask = mask
        self.spat = spat
        self.clobber = clobber
        self.verbose = verbose
        
    def acquire_surface(self, **kwargs):
        return(GMTSurface(data=self.data, src_region=self.region, inc=self.inc, name=self.name, node=self.node,
                          fmt=self.fmt, extend=self.extend, extend_proc=self.extend_proc, weights=self.weights,
                          fltr=self.fltr, sample=self.sample, clip=self.clip, chunk=self.chunk, epsg=self.epsg,
                          archive=self.archive, mask=self.mask, spat=self.spat, clobber=self.clobber, verbose=self.verbose, **kwargs))

    def acquire_triangulate(self, **kwargs):
        return(GMTTriangulate(data=self.data, src_region=self.region, inc=self.inc, name=self.name, node=self.node,
                              fmt=self.fmt, extend=self.extend, extend_proc=self.extend_proc, weights=self.weights,
                              fltr=self.fltr, sample=self.sample, clip=self.clip, chunk=self.chunk, epsg=self.epsg,
                              archive=self.archive, mask=self.mask, spat=self.spat, clobber=self.clobber, verbose=self.verbose, **kwargs))

    def acquire_num(self, **kwargs):
        return(WafflesNum(data=self.data, src_region=self.region, inc=self.inc, name=self.name, node=self.node,
                          fmt=self.fmt, extend=self.extend, extend_proc=self.extend_proc, weights=self.weights,
                          fltr=self.fltr, sample=self.sample, clip=self.clip, chunk=self.chunk, epsg=self.epsg,
                          archive=self.archive, mask=self.mask, spat=self.spat, clobber=self.clobber, verbose=self.verbose, **kwargs))

    def acquire_linear(self, **kwargs):
        return(WafflesLinear(data=self.data, src_region=self.region, inc=self.inc, name=self.name, node=self.node,
                             fmt=self.fmt, extend=self.extend, extend_proc=self.extend_proc, weights=self.weights,
                             fltr=self.fltr, sample=self.sample, clip=self.clip, chunk=self.chunk, epsg=self.epsg,
                             archive=self.archive, mask=self.mask, spat=self.spat, clobber=self.clobber, verbose=self.verbose, **kwargs))

    def acquire_average(self, **kwargs):
        return(WafflesMovingAverage(data=self.data, src_region=self.region, inc=self.inc, name=self.name, node=self.node,
                                    fmt=self.fmt, extend=self.extend, extend_proc=self.extend_proc, weights=self.weights,
                                    fltr=self.fltr, sample=self.sample, clip=self.clip, chunk=self.chunk, epsg=self.epsg,
                                    archive=self.archive, mask=self.mask, spat=self.spat, clobber=self.clobber, verbose=self.verbose, **kwargs))

    def acquire_invdst(self, **kwargs):
        return(WafflesInvDst(data=self.data, src_region=self.region, inc=self.inc, name=self.name, node=self.node,
                             fmt=self.fmt, extend=self.extend, extend_proc=self.extend_proc, weights=self.weights,
                             fltr=self.fltr, sample=self.sample, clip=self.clip, chunk=self.chunk, epsg=self.epsg,
                             archive=self.archive, mask=self.mask, spat=self.spat, clobber=self.clobber, verbose=self.verbose, **kwargs))
    
    def acquire_nearest(self, **kwargs):
        return(WafflesNearest(data=self.data, src_region=self.region, inc=self.inc, name=self.name, node=self.node,
                              fmt=self.fmt, extend=self.extend, extend_proc=self.extend_proc, weights=self.weights,
                              fltr=self.fltr, sample=self.sample, clip=self.clip, chunk=self.chunk, epsg=self.epsg,
                              archive=self.archive, mask=self.mask, spat=self.spat, clobber=self.clobber, verbose=self.verbose, **kwargs))

    def acquire_vdatum(self, **kwargs):
        return(WafflesVdatum(data=self.data, src_region=self.region, inc=self.inc, name=self.name, node=self.node,
                             fmt=self.fmt, extend=self.extend, extend_proc=self.extend_proc, weights=self.weights,
                             fltr=self.fltr, sample=self.sample, clip=self.clip, chunk=self.chunk, epsg=self.epsg,
                             archive=self.archive, mask=self.mask, spat=self.spat, clobber=self.clobber, verbose=self.verbose, **kwargs))
    
    def acquire_module_by_name(self, mod_name, **kwargs):
        if mod_name == 'surface':
            return(self.acquire_surface(**kwargs))
        if mod_name == 'triangulate':
            return(self.acquire_triangulate(**kwargs))
        if mod_name == 'num':
            return(self.acquire_num(**kwargs))
        if mod_name == 'linear':
            return(self.acquire_linear(**kwargs))
        if mod_name == 'average':
            return(self.acquire_average(**kwargs))
        if mod_name == 'nearest':
            return(self.acquire_nearest(**kwargs))
        if mod_name == 'invdst':
            return(self.acquire_invdst(**kwargs))
        if mod_name == 'vdatum':
            return(self.acquire_vdatum(**kwargs))
        
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
\t\t\tappend :<inc> to resample the output to the given <inc>: -E.3333333s:.1111111s
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

  -p, --prefix\t\tSet BASENAME (-O) to PREFIX (append inc/region/year info to output BASENAME).
  -r, --grid-node\tUse grid-node registration, default is pixel-node
  -w, --weights\t\tUse weights provided in the datalist to weight overlapping data.

  -a, --archive\t\tArchive the datalist to the given region.
  -m, --mask\t\tGenerate a data mask raster.
  -s, --spat\t\tGenerate spatial metadata.
  -c, --continue\tDon't clobber existing files.
  -q, --quiet\t\tLower verbosity to a quiet. (overrides --verbose)

  --help\t\tPrint the usage text
  --config\t\tSave the waffles config JSON and major datalist
  --modules\t\tDisply the module descriptions and usage
  --version\t\tPrint the version information

Datalists and data formats:
  A datalist is a file that contains a number of datalist entries, while an entry is a space-delineated line:
  `/path/to/data format weight data,meta,data`

Supported datalist formats: 
  {dl_formats}

Modules (see waffles --modules <module-name> for more info):
  {modules}

CIRES DEM home page: <http://ciresgroups.colorado.edu/coastalDEM>
""".format(cmd=os.path.basename(sys.argv[0]),
           dl_formats=dlim._datalist_fmts_short_desc(),
           modules=_waffles_module_short_desc(),
           wf_version=__version__)

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
    status = 0
    i = 1
    wg = {}
    wg['verbose'] = True
    wg['sample'] = None
    wg['fltr'] = []
    
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
            wg['inc'] = utils.str2inc(incs[0])
            if len(incs) > 1: wg['sample'] = utils.str2inc(incs[1])
            i = i + 1
        elif arg[:2] == '-E':
            incs = arg[2:].split(':')
            wg['inc'] = utils.str2inc(arg[2:].split(':')[0])
            if len(incs) > 1: wg['sample'] = utils.str2inc(incs[1])
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
        
        elif arg == '-w' or arg == '--weights': wg['weights'] = 1
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
                    sys.stderr.write(_waffles_module_long_desc({k: WaffleFactory()._modules[k] for k in (argv[i + 1],)}))
                else: sys.stderr.write(_waffles_module_long_desc(WaffleFactory()._modules))
            except: sys.stderr.write(_waffles_module_long_desc(WaffleFactory()._modules))
            sys.exit(0)
        elif arg == '--help' or arg == '-h':
            sys.stderr.write(waffles_cli_usage)
            sys.exit(0)
        elif arg == '--version' or arg == '-v':
            sys.stdout.write('{}\n'.format(__version__))
            sys.exit(0)
        elif arg[0] == '-':
            print(waffles_cli_usage)
            utils.echo_error_msg('{} is not a valid waffles cli switch'.format(arg))
            sys.exit(0)
        else: dls.append(arg)
        i += 1

    ## ==============================================
    ## Otherwise run from cli options...
    ## set the dem module
    ## ==============================================
    if module is not None:
        mod_opts = {}
        opts = module.split(':')
        if opts[0] in WaffleFactory()._modules.keys():
            mod_opts[opts[0]] = list(opts[1:])
        else:
            utils.echo_error_msg('invalid module name `{}`'.format(opts[0]))
            sys.exit(-1)

        for key in mod_opts.keys():
            mod_opts[key] = [None if x == '' else x for x in mod_opts[key]]
        mod = opts[0]
        try:
            mod_args = utils.args2dict(tuple(mod_opts[mod]))
        except:
            utils.echo_error_msg('invalid module string: "{}"'.format(module))
            sys.exit(-1)
    else: 
        sys.stderr.write(waffles_cli_usage)
        utils.echo_error_msg('''must specify a waffles -M module.''')
        sys.exit(-1)
        
    if WaffleFactory()._modules[mod]['datalist-p']:
        if len(dls) == 0:
            sys.stderr.write(waffles_cli_usage)
            utils.echo_error_msg('''must specify a datalist/entry, try gmrt or srtm for global data.''')
            sys.exit(-1)

    ## ==============================================
    ## check the increment
    ## ==============================================
    if wg['inc'] is None:
        sys.stderr.write(waffles_cli_usage)
        utils.echo_error_msg('''must specify a gridding increment.''')
        sys.exit(-1)
    
    ## ==============================================
    ## set the datalists and names
    ## ==============================================
    wg['data'] = dls

    for i_region in i_regions:
        tmp_region = regions.Region().from_string(i_region)
        if tmp_region.valid_p():
            these_regions.append(tmp_region)
        else:
            tmp_region = regions.ogr_wkts(i_region)
            for i in tmp_region:
                if i.valid_p():
                    these_regions.append(i)
                    
    if len(these_regions) == 0:
        these_regions = []
        utils.echo_error_msg('failed to parse region(s), {}'.format(i_regions))
    else:
        if wg['verbose']: utils.echo_msg('parsed {} region(s)'.format(len(these_regions)))

    for i, this_region in enumerate(these_regions):
        wg['src_region'] = this_region
        if want_prefix or len(these_regions) > 1:
            wg['name'] = utils.append_fn(wg['name'], wg['src_region'], wg['sample'] if wg['sample'] is not None else wg['inc'])

        wf = WaffleFactory(**wg).acquire_module_by_name(mod, **mod_args).generate()
### End
