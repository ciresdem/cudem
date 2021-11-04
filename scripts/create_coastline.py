#!/usr/bin/env python
### create_coastline.py
##
## Copyright (c) 2020 - 2021 CIRES Coastal DEM Team
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
## ==============================================
## Waffles Coastline module
## generate a coastline (wet/dry mask)
##
## Using various sources (local datalists/USGS NHD/GMRT/ETc.), generate
## a detailed coastline at the given resolution.
## using a local datalist can take a long time if there is a lot
## of lidar to process. Speed things up by limiting the data
## by weight or z (w-range, z-range), and only including the best
## coastal data.
##
## GMT or GMRT will fill in the gaps where local-data and NHD can't fill;
## hopefully this is just to fill areas off-shore and far-inland.
## ==============================================
### Code:

import os
from osgeo import gdal
from osgeo import ogr
from cudem import utils
from cudem import waffles

class WafflesCoastline(waffles.Waffle):
    def __init__(self, want_nhd=True, want_gmrt=False, **kwargs):
        """Generate a coastline polygon from various sources."""

        super().__init__(**kwargs)
        self.want_nhd = want_nhd
        self.want_gmrt = want_gmrt

        self.w_name = '{}_w'.format(self.name)
        self.w_mask = '{}.tif'.format(self.w_name)

        self.u_name = '{}_u'.format(self.name)
        self.u_mask = '{}.tif'.format(self.u_name)

        self.g_name = '{}_g'.format(self.name)
        self.g_mask = '{}.tif'.format(self.g_name)
        
        self.coast_array = None

        self.ds_config = None

    def run(self):
        self._burn_region()
        self._load_coast_mask()
        self._load_nhd()
        self._load_background()
        self._write_coast_array()
        return(self)
        
    def _burn_region(self):
        """wet/dry datalist mask or burn region."""

        if len(self.data) < 0:
            c_region = self.d_region
            c_region.zmin = -1
            c_region.zmax = 1

            waffles.WafflesNum(data=self.data, src_region=c_region, inc=self.inc, name=self.w_name,
                               extend=self.extend, extend_proc=self.extend_proc, weights=self.weights,
                               sample=self.sample, clip=self.clip, epsg=self.epsg, verbose=self.verbose,
                               **kwargs, mode='w').run()
        else:
            c_region.export_as_ogr('region_buff.shp')
            xsize, ysize, gt = c_region.geo_transform(x_inc=self.inc)
            
            utils.run_cmd('gdal_rasterize -ts {} {} -te {} -burn -9999 -a_nodata -9999 \
            -ot Int32 -co COMPRESS=DEFLATE -a_srs EPSG:{} region_buff.shp {}\
            '.format(xsize, ysize, c_region.format('te'), self.epsg, self.w_mask), verbose=self.verbose)

    def _load_coast_mask(self):
        """load the wet/dry mask array"""

        ds = gdal.Open(self.w_mask)
        if ds is not None:
            self.ds_config = demfun.gather_infos(ds)
            this_region = regions.Region().from_geo_transform(self.ds_config['geoT'], self.ds_config['nx'], self.ds_config['ny'])
            self.coast_array = ds.GetRasterBand(1).ReadAsArray(0, 0, self.ds_config['nx'], self.ds_config['ny'])
            ds = None
        else:
            utils.echo_error_msg('could not open {}'.format(self.w_mask))
            sys.exit()
        utils.remove_glob('{}*'.format(self.w_mask))

    def _load_coast_shape(self):
        """Input coastline shapefile `coastpoly`"""

        raise(NotImplementedError)

    def _load_nhd(self):
        """USGS NHD (HIGH-RES U.S. Only)
        Fetch NHD (NHD High/Plus) data from TNM to fill in near-shore areas. 
        High resoultion data varies by location...
        """

        self.p_region.export_as_ogr('region_buff.shp')
        xsize, ysize, gt = self.p_region.geo_transform(x_inc=self.inc)

        utils.run_cmd('gdal_rasterize -ts {} {} -te {} -burn -9999 -a_nodata -9999 \
        -ot Int32 -co COMPRESS=DEFLATE -a_srs EPSG:{} region_buff.shp {}\
        '.format(xsize, ysize, self.p_region.format('te'), self.epsg, self.u_mask), verbose=self.verbose)
        utils.remove_glob('region_buff.*')

        this_tnm = fetches.tnm.TheNationalMap(src_region=self.p_region, weight=self.weight, verbose=self.verbose, where="Name LIKE '%Hydro%'", extents='HU-4 Subregion,HU-8 Subbasin').run()

        #fl = fetches._fetch_modules['tnm'](waffles_proc_region(wg), ["Name LIKE '%Hydro%'"], None, True)
        r_shp = []
        for result in this_tnm.results:
            #fl._parse_results(e = 'HU-2 Region,HU-4 Subregion,HU-8 Subbasin'):
            if f_utils.Fetch(result[0], verbose=self.verbose).fetch_file(os.path.join(result[2], result[1])) == 0:
                gdb_zip = os.path.join(result[2], result[1])
                gdb_files = utils.unzip(gdb_zip)
                gdb_bn = os.path.basename('.'.join(gdb_zip.split('.')[:-1]))
                gdb = gdb_bn + '.gdb'

                utils.run_cmd('ogr2ogr {}_NHDArea.shp {} NHDArea -clipdst {} -overwrite 2>&1\
                '.format(gdb_bn, gdb, regions.region_format(wg['region'], 'ul_lr')), verbose=False)
                if os.path.exists('{}_NHDArea.shp'.format(gdb_bn)):
                    r_shp.append('{}_NHDArea.shp'.format(gdb_bn))
                utils.run_cmd('ogr2ogr {}_NHDPlusBurnWaterBody.shp {} NHDPlusBurnWaterBody -clipdst {} -overwrite 2>&1\
                '.format(gdb_bn, gdb, self.p_region.format('ul_lr')), verbose=False)
                if os.path.exists('{}_NHDPlusBurnWaterBody.shp'.format(gdb_bn)):
                    r_shp.append('{}_NHDPlusBurnWaterBody.shp'.format(gdb_bn))
                utils.run_cmd('ogr2ogr {}_NHDWaterBody.shp {} NHDWaterBody -where "FType = 390" -clipdst {} -overwrite 2>&1\
                '.format(gdb_bn, gdb, self.p_region.format('ul_lr')), verbose=False)
                if os.path.exists('{}_NHDWaterBody.shp'.format(gdb_bn)):
                    r_shp.append('{}_NHDWaterBody.shp'.format(gdb_bn))
                utils.remove_glob(gbd)
            else: utils.echo_error_msg('unable to fetch {}'.format(result))

            [utils.run_cmd('ogr2ogr -skipfailures -update -append nhdArea_merge.shp {} 2>&1\
            '.format(shp), verbose=False) for shp in r_shp]
            utils.run_cmd('gdal_rasterize -burn 1 nhdArea_merge.shp {}'.format(self.u_mask), verbose=True)
            utils.remove_glob('nhdArea_merge.*', 'NHD_*', *r_shp)

        ## ==============================================
        ## update wet/dry mask with nhd data
        ## ==============================================
        utils.echo_msg('filling the coast mask with NHD data...')
        c_ds = gdal.Open(self.u_mask)
        c_ds_arr = c_ds.GetRasterBand(1).ReadAsArray()
        #c_ds = gdal.Open(u_mask)
        for this_xyz in demfun.parse(c_ds):
            xpos, ypos = utils._geo2pixel(this_xyz.x, this_xyz.y, self.dst_gt)
            try:
                if self.coast_array[ypos, xpos] == self.ds_config['ndv']:
                    if this_xyz.z == 1: self.coast_array[ypos, xpos] = 0
            except: pass
        c_ds = None            
        utils.remove_glob('{}*'.format(self.u_mask))

    def _load_background(self):
        """GSHHG/GMRT - Global low-res
        Used to fill un-set cells.
        """
        
        if wg['gc']['GMT'] is not None and not self.want_gmrt:
            utils.run_cmd('gmt grdlandmask {} -I{} -r -Df -G{}=gd:GTiff -V -N1/0/1/0/1\
            '.format(self.p_region.format('gmt'), self.inc, self.g_mask), verbose=self.verbose)
        else:
            this_gmrt = gmrt.GMRT(src_region=self.p_region, weight=self.weight, verbose=self.verbose, layer='topo-mask').run()
            #gmrt_tif = this_gmrt.results[0]
            this_gmrt.fetch_results()
            
            utils.run_cmd('gdalwarp {} {} -tr {} {} -overwrite'.format(gmrt_tif, g_mask, wg['inc'], wg['inc']), verbose = True)
            #utils.remove_glob(gmrt_tif)

        ## ==============================================
        ## update wet/dry mask with gsshg/gmrt data
        ## speed up!
        ## ==============================================
        utils.echo_msg('filling the coast mask with gsshg/gmrt data...')
        c_ds = gdal.Open(self.g_mask)
        c_ds_arr = c_ds.GetRasterBand(1).ReadAsArray()
        #c_ds = gdal.Open(self.g_mask)
        for this_xyz in gdalfun.gdal_parse(c_ds):
            xpos, ypos = utils._geo2pixel(this_xyz.x, this_xyz.y, self.dst_gt)
            try:
                if self.coast_array[ypos, xpos] == self.ds_config['ndv']:
                    if this_xyz.z == 1:
                        self.coast_array[ypos, xpos] = 0
                    elif this_xyz.z == 0:
                        self.coast_array[ypos, xpos] = 1
            except: pass
        c_ds = None
        utils.remove_glob('{}*'.format(self.g_mask))

    def _write_coast_array(self):
        """write coast_array to file
        1 = land
        0 = water
        """
        
        utils.gdal_write(coast_array, '{}.tif'.format(self.name), self.ds_config)

    def _write_coast_poly(self):
        """convert to coast_array vector"""

        tmp_ds = ogr.GetDriverByName('ESRI Shapefile').CreateDataSource('tmp_c_{}.shp'.format(self.name))
        if tmp_ds is not None:
            tmp_layer = tmp_ds.CreateLayer('tmp_c_{}'.format(self.name), None, ogr.wkbMultiPolygon)
            tmp_layer.CreateField(ogr.FieldDefn('DN', ogr.OFTInteger))
            demfun.polygonize('{}.tif'.format(self.name), tmp_layer, verbose=self.verbose)
            tmp_ds = None
        utils.run_cmd('ogr2ogr -dialect SQLITE -sql "SELECT * FROM tmp_c_{} WHERE DN=0 order by ST_AREA(geometry) desc limit 8"\
        {}.shp tmp_c_{}.shp'.format(self.name, self.name, self.name), verbose = True)
        utils.remove_glob('tmp_c_{}.*'.format(self.name))
        utils.run_cmd('ogrinfo -dialect SQLITE -sql "UPDATE {} SET geometry = ST_MakeValid(geometry)" {}.shp\
        '.format(self.name, self.name))


#if __name__ == '__main__':
    
### End
