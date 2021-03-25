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
from cudem import xyzs

__version__ = '0.10.0'

def waffles_append_fn(bn, src_region, inc):
    return('{}{}_{}_{}v1'.format(bn, utils.inc2str_inc(inc), src_region.format('fn'), utils.this_year()))

def xyz_block(src_xyz, src_region, inc, weights = False, verbose = False):
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
                sumArray[ypos, xpos] += z
                ptArray[ypos, xpos] += 1
                if weights: wtArray[ypos, xpos] += this_xyz[3]
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
            out_xyz = xyzs.XYZPoint(x=geo_x, y=geo_y, z=outarray[y,x])
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
    ds_config = utils.gdal_set_infos(xcount, ycount, xcount * ycount, dst_gt, utils.sr_wkt(epsg), gdt, -9999, dst_format)
    for this_xyz in src_xyz:
        if regions.xyz_in_region_p(this_xyz, src_region):
            xpos, ypos = utils._geo2pixel(this_xyz.x, this_xyz.y, dst_gt)
            try:
                if mode == 'm' or mode == 'w':
                    sumArray[ypos, xpos] += z
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
        
    def open(self, update=False):
        try:
            if update:
                self.src_ds = gdal.Open(self.fn, gdal.GA_Update)
            else:
                self.src_ds = gdal.Open(self.fn)

            gt = self.src_ds.GetGeoTransform()

            if self.region is not None and self.region.valid_p(check_xy = True):
                self.srcwin = self.region.srcwin(gt, self.src_ds.RasterXSize, self.src_ds.RasterYSize)
            else: self.srcwin = (0, 0, self.src_ds.RasterXSize, self.src_ds.RasterYSize)

            self.gt = (gt[0] + (self.srcwin[0]*gt[1]), gt[1], 0., gt[3] + (self.srcwin[1]*gt[5]), 0., gt[5])

            self.ds_open_p = True
        except:
            utils.echo_error_msg('could not open raster file {}'.format(self.fn))
            self.src_ds = None
            self.ds_open_p = False
        return(self)

    def close(self):
        self.src_ds = None
        self.ds_open_p = False

    def gather_infos(self, scan=False):
        """gather information from `src_ds` GDAL dataset

        Returns:
          raster_parser: self
        """

        if self.ds_open_p:
            src_band = self.src_ds.GetRasterBand(1)
            self.x_count = self.srcwin[2]
            self.y_count = self.srcwin[3]
            self.dt = src_band.DataType
            self.dtn = gdal.GetDataTypeName(src_band.DataType)
            self.ndv = src_band.GetNoDataValue()
            if self.ndv is None: self.ndv = -9999
            self.fmt = self.src_ds.GetDriver().ShortName
            self.zr = None

            if scan:
                src_arr = src_band.ReadAsArray(srcwin[0], self.srcwin[1], self.srcwin[2], self.srcwin[3])
                self.zr = (np.amin(src_arr), np.amax(src_arr))
                src_arr = None
        return(self)

    def set_nodata(self, nodata=-9999):
        """set the nodata value of gdal file src_fn
        
        returns 0
        """
        
        if self.ds_open_p:
            band = self.src_ds.GetRasterBand(1)
            band.SetNoDataValue(nodata)
            self.ndv = nodata
    
    def set_epsg(self):
        if self.ds_open_p:
            proj = self.src_ds.GetProjectionRef()
            src_srs = osr.SpatialReference()
            src_srs.ImportFromWkt(proj)
            src_srs.AutoIdentifyEPSG()
            srs_auth = src_srs.GetAuthorityCode(None)
            epsg = utils.int_or(srs_auth)
            if epsg is None:
                self.src_ds.SetProjection(utils.sr_wkt(int(self.epsg)))
            else: self.epsg = epsg

    def set_infos(self):
        """set a datasource config dictionary
            
        returns gdal_config dict."""
        if self.ds_open_p:
            src_band = self.src_ds.GetRasterBand(1)
            return({'nx': self.srcwin[2], 'ny': self.srcwin[3], 'nb': self.srcwin[2] * self.srcwin[3], 'geoT': self.gt,
                    'proj': self.src_ds.GetProjectionRef(), 'dt': src_band.DataType, 'ndv': src_band.GetNoDataValue(),
                    'fmt': self.src_ds.GetDriver().ShortName})
        else: return(None)
    
    def set_metadata(self, cudem=False):
        """add metadata to the waffled raster
    
        Args: 
          cudem (bool): add CUDEM metadata
        """

        if self.ds_open_p:
            md = self.src_ds.GetMetadata()
            if self.node == 'pixel':
                md['AREA_OR_POINT'] = 'Area'
            else: md['AREA_OR_POINT'] = 'Point'
            md['TIFFTAG_DATETIME'] = '{}'.format(utils.this_date())
            if cudem:
                md['TIFFTAG_COPYRIGHT'] = 'DOC/NOAA/NESDIS/NCEI > National Centers for Environmental Information, NESDIS, NOAA, U.S. Department of Commerce'
                md['TIFFTAG_IMAGEDESCRIPTION'] = 'Topography-Bathymetry; NAVD88'
            self.src_ds.SetMetadata(md)
        else: utils.echo_error_msg('failed to set metadata')

    def valid_p(self):
        if not os.path.exists(self.fn): return(False)
        self.open()
        gdi = self.gather_infos(this_dem, scan = True)
        self.close()
        
        if gdi is not None:
            if np.isnan(gdi['zr'][0]):
                #utils.remove_glob(self.fn)
                return(False)
        else:
            return(False)
        return(True)

    def filter_(self, fltr=1, fltr_val=None, split_val=None, mask=None):
        """filter raster using smoothing factor `fltr`; optionally
        only smooth bathymetry (sub-zero) using a split_val of 0.

        Args: 
          fltr (int): the filter to use, 1, 2 or 3
          flt_val (varies): the filter value, varies by filter.
          split_val (float): an elevation value (only filter below this value)

        Returns:
          0 for success or -1 for failure
        """

        was_open = False
        if self.ds_open_p:
            self.close()
            was_open = True

        if os.path.exists(self.fn):
            if split_val is not None:
                dem_u, dem_l = gdalfun.gdal_split(self.fn, split_val)
            else: dem_l = self.fn

            if int(fltr) == 1: out, status = utils.gdal_blur(dem_l, 'tmp_fltr.tif', fltr_val if fltr_val is not None else 10)
            elif int(fltr) == 2: out, status = utils.gmt_grdfilter(dem_l, 'tmp_fltr.tif=gd+n-9999:GTiff',
                                                                    dist = fltr_val if fltr_val is not None else '1s',
                                                                    node = node, verbose = False)
            elif int(fltr) == 3: out, status = utils.gdal_filter_outliers(dem_l, 'tmp_fltr.tif',
                                                                            fltr_val if fltr_val is not None else 10)
            else: out, status = utils.gdal_blur(dem_l, 'tmp_fltr.tif', fltr_val if fltr_val is not None else 10)
            if status != 0: return(status)

            if split_val is not None:
                ds = gdal.Open(self.fn)
                ds_config = utils.gdal_gather_infos(ds)
                msk_arr = ds.GetRasterBand(1).ReadAsArray()
                msk_arr[msk_arr != ds_config['ndv']] = 1
                msk_arr[msk_arr == ds_config['ndv']] = np.nan
                ds = None
                u_ds = gdal.Open(dem_u)
                if u_ds is not None:
                    l_ds = gdal.Open('tmp_fltr.tif')
                    if l_ds is not None:
                        u_arr = u_ds.GetRasterBand(1).ReadAsArray()
                        l_arr = l_ds.GetRasterBand(1).ReadAsArray()
                        u_arr[u_arr == ds_config['ndv']] = 0
                        l_arr[l_arr == ds_config['ndv']] = 0
                        ds_arr = (u_arr + l_arr) * msk_arr
                        ds_arr[np.isnan(ds_arr)] = ds_config['ndv']
                        utils.gdal_write(ds_arr, 'merged.tif', ds_config)
                        l_ds = None
                        utils.remove_glob(dem_l)
                    u_ds = None
                    utils.remove_glob(dem_u)
                os.rename('merged.tif', 'tmp_fltr.tif')
            os.rename('tmp_fltr.tif', self.fn)
            if was_open:
                self.open()
            return(0)
        else: return(-1)

    
    def sample(self):
        was_open = False
        if self.ds_open_p:
            self.close()
            was_open = True

        out, status = utils.run_cmd('gdalwarp -tr {:.10f} {:.10f} {} -r bilinear -te {} __tmp__.tif\
        '.format(self.sample_inc, self.sample_inc, self.fn, self.region.format('te')))

        if status == 0:
            os.rename('__tmp__.tif', '{}'.format(self.fn))
                                    
        if was_open:
            self.open()
            
    def clip(self):
        pass

    def cut(self, dst_fn):
        if self.ds_open_p:
            ds_arr = self.src_ds.GetRasterBand(1).ReadAsArray(self.srcwin[0], self.srcwin[1], self.srcwin[2], self.srcwin[3])
            out_ds_config = self.set_infos()
            return(utils.gdal_write(ds_arr, dst_fn, out_ds_config))
        else: return(None, -1)

    def translate(self):
        #if self.fmt is not None:
        #    pass
        pass

    def process(self):

        self.open()
        
        if len(self.filters) > 0:
            utils.echo_msg('filtering')
            for f in self.filters:
                fltr_val = None
                split_val = None
                fltr_opts = f.split(':')
                fltr = fltr_opts[0]
                if len(fltr_opts) > 1:
                    fltr_val = fltr_opts[1]
                if len(fltr_opts) > 2:
                    split_val= fltr_opts[2]
                    
                self.filter_(fltr=fltr, fltr_val=fltr_val)
            
        if self.sample_inc is not None:
            self.sample()
            
        if self.clip_poly is not None:
            self.clip()
            
        if self.srcwin[0] != 0 or self.srcwin[1] != 0 or self.srcwin[2] != self.src_ds.RasterXSize or self.srcwin[3] != self.src_ds.RasterYSize:
            self.cut('__tmp_cut__.tif')
            self.close()
            os.rename('__tmp_cut__.tif', '{}'.format(self.fn))
        else:
            self.close()
            
        self.open(update=True)

        self.set_nodata()
        self.set_epsg()
        self.set_metadata(cudem=self.cudem)
        self.close()
        
        if self.fmt != 'GTiff':
            self.translate()
            
class WaffledVector(Waffled):
    """Providing a waffled VECTOR file processor."""
    
    def __init__(self, **kwargs):
        pass

    def set_metadata(self):
        pass
    
    def process(self):
        pass

class Waffle:
    
    def __init__(self, data=[], src_region=None, inc=None, name='waffles_dem',
                 node='pixel', fmt='GTiff', extend=0, extend_proc=20, weights=None,
                 fltr=[], sample=None, clip=None, chunk=None, epsg=4326,
                 verbose=False, archive=False, mask=False, clobber=True):

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

        if self.node == 'grid':
            self.region = self.region.buffer(self.inc*.5)
        self.p_region = self.proc_region()
        self.d_region = self.dist_region()

        self.data = [dlim.DatasetFactory(
            fn=" ".join(['-' if x == "" else x for x in dl.split(":")]),
            src_region=self.p_region, verbose=self.verbose, weight=self.weights,
            epsg=self.epsg).acquire_dataset() for dl in self.data]

        for d in self.data:
            d.parse()

        self.waffle_out = {}
        self.waffled = False

        out, status = utils.run_cmd('gmt gmtset IO_COL_SEPARATOR = SPACE', verbose = False)
        
    def proc_region(self):
        
        pr = regions.Region().from_region(self.region)
        return(pr.buffer((self.inc*self.extend_proc) + (self.inc*self.extend)))

    def dist_region(self):
            
        pr = regions.Region().from_region(self.region)
        return(pr.buffer((self.inc*self.extend)))

    def run(self):
        pass

    def generate(self):

        self.run()
        if self.mask:
            WaffledRaster(fn='{}_m.tif'.format(self.name), epsg=self.epsg, fmt=self.fmt, sample_inc=self.sample, src_region = self.region).process()
        self.waffled = True

        # for out in self.waffle_out:
        #     out.fmt = self.fmt
        #     out.sample_inc = self.sample
        #     out.filters = self.fltr
        #     out.region = self.region
        #     utils.echo_msg('processing {}'.format(out.fn))
        #     out.process()
                    
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
        ds_config = utils.gdal_set_infos(
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
        
    def yield_xyz(self, **kwargs):

        for xdl in self.data:
            xyz_yield = xdl.yield_xyz()
            if self.mask:
                xyz_yield = self.xyz_mask(xyz_yield, '{}_m.tif'.format(self.name))
            if self.archive:
                xyz_yield = self.xyz_archive(xyz_yield)
        
            for xyz in xyz_yield:
                yield(xyz)
                
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
        
        WaffledRaster(fn='{}.tif'.format(self.name), epsg=self.epsg, fmt=self.fmt, sample_inc=self.sample, src_region = self.region).process()
        
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
        
        dem_tri_cmd = ('gmt blockmean {} -I{:.10f}{} -V -r | gmt triangulate -V {} -I{:.10f} -G{}.tif=gd:GTiff -r\
        '.format(self.p_region.format('gmt'), self.inc, ' -Wi' if self.weights else '', self.p_region.format('gmt'), self.inc, self.name))
        out, status = utils.run_cmd(dem_tri_cmd, verbose = self.verbose, data_fun = lambda p: self.dump_xyz(dst_port=p, encode=True))

        WaffledRaster(fn='{}.tif'.format(self.name), epsg=self.epsg, fmt=self.fmt, sample_inc=self.sample, src_region = self.region).process()

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
        
        self.mod = 'triangulate'
        self.mode = mode
        
    def run(self):

        if self.mode.startswith('A'):
            if self.gc['GMT'] is None:
                utils.echo_error_msg('GMT must be installed to use the Mode `A` with the NUM module')
                return(None, -1)

            dem_xyz2grd_cmd = ('gmt xyz2grd -{} -V {} -I{:.10f} -G{}.tif=gd:GTiff -r\
            '.format(self.p_region.format('gmt'), self.inc, ' -Wi' if self.weights else '',
                     self.mode, self.p_region.format('gmt'), self.inc, self.name))
            out, status = utils.run_cmd(dem_xyz2grd_cmd, verbose=self.verbose, data_fun=lambda p: self.dump_xyz(dst_port=p, encode=True))
        else:
            dly = self.yield_xyz()
            
            if self.weights: dly = xyz_block(dly, self.p_region, self.inc, weights = True)
            out, status = xyz2gdal(dly, '{}.tif'.format(self.name), self.p_region, self.inc, dst_format=self.fmt, mode=self.mode, verbose=self.verbose)

        WaffledRaster(fn='{}.tif'.format(self.name), epsg=self.epsg, fmt=self.fmt, sample_inc=self.sample, src_region = self.region).process()        
        
class WaffleFactory:

    _modules = {
        'surface': {
            'name': 'surface',
            'datalist-p': True,
        },
        'triangulate': {
            'name': 'triangulate',
            'datalist-p': True
        },
        'num': {
            'name': 'num',
            'datalist-p': True
        },
    }
    
    def __init__(self, data=[], src_region=None, inc=None, name='waffles_dem',
                 node='pixel', fmt='GTiff', extend=0, extend_proc=20, weights=None,
                 fltr=[], sample = None, clip=None, chunk=None, epsg=4326,
                 verbose = False, archive=False, mask=False, clobber=True):

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
        self.clobber = clobber
        self.verbose = verbose
        
    def acquire_surface(self, **kwargs):
        return(GMTSurface(data=self.data, src_region=self.region, inc=self.inc, name=self.name, node=self.node,
                          fmt=self.fmt, extend=self.extend, extend_proc=self.extend_proc, weights=self.weights,
                          fltr=self.fltr, sample=self.sample, clip=self.clip, chunk=self.chunk, epsg=self.epsg,
                          archive=self.archive, mask=self.mask, clobber=self.clobber, verbose=self.verbose, **kwargs))

    def acquire_triangulate(self, **kwargs):
        return(GMTTriangulate(data=self.data, src_region=self.region, inc=self.inc, name=self.name, node=self.node,
                              fmt=self.fmt, extend=self.extend, extend_proc=self.extend_proc, weights=self.weights,
                              fltr=self.fltr, sample=self.sample, clip=self.clip, chunk=self.chunk, epsg=self.epsg,
                              archive=self.archive, mask=self.mask, clobber=self.clobber, verbose=self.verbose, **kwargs))

    def acquire_num(self, **kwargs):
        return(WafflesNum(data=self.data, src_region=self.region, inc=self.inc, name=self.name, node=self.node,
                          fmt=self.fmt, extend=self.extend, extend_proc=self.extend_proc, weights=self.weights,
                          fltr=self.fltr, sample=self.sample, clip=self.clip, chunk=self.chunk, epsg=self.epsg,
                          archive=self.archive, mask=self.mask, clobber=self.clobber, verbose=self.verbose, **kwargs))

    
    def acquire_module_by_name(self, mod_name, **kwargs):
        if mod_name == 'surface':
            return(self.acquire_surface(**kwargs))
        if mod_name == 'triangulate':
            return(self.acquire_triangulate(**kwargs))
        if mod_name == 'num':
            return(self.acquire_num(**kwargs))

## ==============================================
## Command-line Interface (CLI)
## $ waffles
##
## waffles cli
## ==============================================
_waffles_module_short_desc = lambda: ', '.join(
    ['{}'.format(key) for key in WaffleFactory()._modules])
#_waffles_module_long_desc = lambda x: 'waffles modules:\n% waffles ... -M <mod>:key=val:key=val...\n\n  ' + '\n  '.join(['\033[1m{:14}\033[0m{}\n'.format(key, x[key]['description']) for key in x]) + '\n'
#_waffles_module_short_desc = lambda x: ', '.join(['{}'.format(key) for key in x])

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
  -M, --module\t\tDesired DEM MODULE and options. (see available Modules below)
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
\t\t\te.g. -T1:10:split_value=0 to smooth bathymetry using Gaussian filter
  -Z, --z-region\t\tRestrict data processing to records that fall within the z-region
\t\t\tUse '-' to indicate no bounding range; e.g. -Z-/0 will restrict processing to data
\t\t\trecords whose z value is below zero.
  -C, --clip\t\tCLIP the output to the clip polygon -C<clip_ply.shp:invert=False>
  -K, --chunk\t\tProcess the region in CHUNKs. -K<chunk-level>
  -W, --w-region\tRestrict data processing to records that fall within the w-region (weight).
\t\t\tUse '-' to indicate no bounding range; e.g. -W1/- will restrict processing to data
\t\t\trecords whose weight value is at least 1.
  -G, --wg-config\tA waffles config JSON file. If supplied, will overwrite all other options.
\t\t\tgenerate a waffles_config JSON file using the --config flag.

  -p, --prefix\t\tSet BASENAME to PREFIX (append inc/region/year info to output BASENAME).
  -r, --grid-node\tUse grid-node registration, default is pixel-node
  -w, --weights\t\tUse weights provided in the datalist to weight overlapping data.

  -a, --archive\t\tArchive the datalist to the given region.
  -m, --mask\t\tGenerate a data mask raster.
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
    status = 0
    i = 1
    wg = {}
    wg['verbose'] = True
    wg['sample'] = None
    
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
            wg['inc'] = utils.gmt_inc2inc(incs[0])
            if len(incs) > 1: wg['sample'] = utils.gmt_inc2inc(incs[1])
            i = i + 1
        elif arg[:2] == '-E':
            incs = arg[2:].split(':')
            wg['inc'] = utils.gmt_inc2inc(arg[2:].split(':')[0])
            if len(incs) > 1: wg['sample'] = utils.gmt_inc2inc(incs[1])
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
        elif arg == '-c' or arg == '--continue': wg['clobber'] = False
        elif arg == '-r' or arg == '--grid-node': wg['node'] = 'grid'

        elif arg == '--quiet' or arg == '-q': wg['verbose'] = False
        elif arg == '--config': want_config = True
        elif arg == '--modules' or arg == '-m':
            try:
                if argv[i + 1] in _waffles_modules.keys():
                    sys.stderr.write(_waffles_module_long_desc({k: _waffles_modules[k] for k in (argv[i + 1],)}))
                else: sys.stderr.write(_waffles_module_long_desc(_waffles_modules))
            except: sys.stderr.write(_waffles_module_long_desc(_waffles_modules))
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
        mod_args = utils.args2dict(tuple(mod_opts[mod]))
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
            wg['name'] = waffles_append_fn(wg['name'], wg['src_region'], wg['sample'] if wg['sample'] is not None else wg['inc'])

        wf = WaffleFactory(**wg).acquire_module_by_name(mod, **mod_args).generate()
### End
