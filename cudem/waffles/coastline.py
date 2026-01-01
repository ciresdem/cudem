### coastline.py
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
##
## coastline.py is part of CUDEM
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
## Generate a coastline (land/etc-mask) using a variety of sources.
##
### Code:

import os
import numpy as np
from osgeo import gdal
from osgeo import ogr

from cudem import utils
from cudem import gdalfun
from cudem.fetches import fetches
from cudem.fetches import osm
from cudem.waffles.waffles import Waffle

class WafflesCoastline(Waffle):
    """COASTLINE (land/etc-mask) generation
    
    Generate a coastline (land/etc-mask) using a variety of sources. 
    Output raster will mask land-areas as 1 and oceans/(lakes/buildings) as 0.
    Output vector will polygonize land-areas.
    
    Parameters:
    -----------
    want_osm_coast (bool): use OSM to fill background
    want_enc_coast (bool): use ENC charts to fill background
    want_gmrt (bool): use GMRT to fill background
    want_copernicus (bool): use COPERNICUS to fill background
    want_cudem (bool): use CUDEM to fill background
    want_nhd (bool): use high-resolution NHD to fill US coastal zones
    want_lakes (bool): mask LAKES using HYDROLAKES
    invert_lakes (bool): invert the lake mask
    want_osm_buildings (bool): mask BUILDINGS using OSM
    want_bing_buildings (bool): mask BUILDINGS using BING
    want_wsf_buildings (bool): mask BUILDINGS using WSF
    osm_tries (int): OSM max server attempts
    min_building_length (float): only use buildings larger than val
    invert (bool): invert the output results
    polygonize (bool/int): polygonize the output (int limits poly count)
    min_weight (float): weight applied to fetched coastal data
    """
    
    def __init__(
            self,
            want_osm_coast=True,
            want_enc_coast=False,
            want_nhd=True,
            want_nhd_plus=False,
            want_copernicus=False,
            want_cudem=False,
            want_gmrt=False,
            want_lakes=False,
            invert_lakes=False,
            want_buildings=False,
            min_building_length=None,
            want_osm_planet=False,
            want_bing_buildings=False,
            want_osm_buildings=False,
            want_wsf_buildings=False,
            invert=False,
            polygonize=True,
            osm_tries=5,
            want_wsf=False,
            min_weight=1,
            **kwargs
    ):
        super().__init__(**kwargs)
        self.want_osm_coast = want_osm_coast
        self.want_enc_coast = want_enc_coast
        self.want_nhd = want_nhd
        self.want_nhd_plus = want_nhd_plus
        self.want_gmrt = want_gmrt
        self.want_copernicus = want_copernicus
        self.want_cudem = want_cudem
        self.want_lakes = want_lakes
        self.invert_lakes = invert_lakes
        self.want_buildings = want_buildings
        self.want_wsf = want_wsf
        self.want_bing_buildings = want_bing_buildings
        self.want_osm_buildings = want_osm_buildings
        self.want_wsf_buildings = want_wsf_buildings
        self.min_building_length = min_building_length
        self.want_osm_planet = want_osm_planet
        self.invert = invert
        self.polygonize = polygonize
        self.osm_tries = utils.int_or(osm_tries, 5)
        self.coast_array = None
        self.ds_config = None
        self.min_weight = utils.float_or(min_weight, 1)

        
    def fetch_data(self, fetches_module, check_size=True):
        """Fetch remote data using the FetchesFactory."""
        
        this_fetches = fetches.FetchesFactory(
            mod=fetches_module,
            src_region=self.wgs_region,
            verbose=self.verbose,
            outdir=self.cache_dir,
            callback=fetches.fetches_callback
        )._acquire_module()        
        
        this_fetches.run()
        
        fr = fetches.fetch_results(this_fetches, check_size=check_size)
        fr.daemon = True
        fr.start()
        fr.join()

        return fr

    
    def run(self):
        """Execute the coastline generation."""
        
        ## Setup Regions
        self.f_region = self.p_region.copy()
        self.f_region.buffer(pct=5, x_inc=self.xinc, y_inc=self.yinc)
        self.f_region.src_srs = self.dst_srs
        self.wgs_region = self.f_region.copy()
        self.wgs_srs = 'epsg:4326'
        
        if self.dst_srs is not None:
            self.wgs_region.warp(self.wgs_srs)
        else:
            self.dst_srs = self.wgs_srs
            
        horz_epsg, _ = gdalfun.epsg_from_input(self.dst_srs)
        self.cst_srs = horz_epsg

        self._load_coast_mask()
        self.mem_config = self.ds_config.copy()
        self.mem_config['dt'] = gdal.GDT_Byte
        
        self.want_nhd = True if self.want_nhd_plus else self.want_nhd

        ## Load Sources
        if self.want_osm_coast: self._load_osm_coast()
        if self.want_enc_coast: self._load_enc_coast()
        if self.want_gmrt: self._load_gmrt()
        if self.want_copernicus: self._load_copernicus()
        if self.want_cudem: self._load_cudem()
        if self.want_nhd: self._load_nhd()
        if self.want_lakes: self._load_lakes()

        if self.want_buildings or self.want_bing_buildings:
            self._load_bing_bldgs()

        if self.want_osm_buildings:
            self._load_osm_bldgs()

        if self.want_wsf or self.want_wsf_buildings:
            self._load_wsf_bldgs()
            
        if self.want_stack:
            self._load_data()
            
        if self.verbose:
            utils.echo_msg(
                'finalizing array for region {} at {} {}...'.format(
                    self.p_region.format('gmt'),
                    self.ds_config['nx'],
                    self.ds_config['ny']
                )
            )
            
        self._finalize_array()
        self._write_coast_array()
        
        if self.polygonize:
            if utils.int_or(self.polygonize) is not None:
                self._write_coast_poly(poly_count=self.polygonize)
            else:
                self._write_coast_poly()

        if np.all(self.coast_array == 0) or np.all(self.coast_array == 1):
            return None
        else:
            return self

        
    def _finalize_array(self):
        """Finalize binary mask values."""
        
        self.coast_array[self.coast_array > 0] = 1
        self.coast_array[self.coast_array <= 0] = 0
        
        if self.invert:
            ## Swap 0 and 1 using a temp value of 2
            self.coast_array[self.coast_array == 0] = 2
            self.coast_array[self.coast_array == 1] = 0
            self.coast_array[self.coast_array == 2] = 1

            
    def _load_coast_mask(self):
        """Create initial zero-grid."""
        
        xcount, ycount, gt = self.p_region.geo_transform(
            x_inc=self.xinc, y_inc=self.yinc
        )
        self.ds_config = gdalfun.gdal_set_infos(
            xcount, ycount, xcount * ycount, gt, self.cst_srs,
            gdal.GDT_Float32, self.ndv, 'GTiff', None, None
        )        
        self.coast_array = np.zeros((ycount, xcount))

            
    def _load_osm_coast(self):
        """Load OpenStreetMap Coastline."""
        
        this_cst = self.fetch_data('osm:q=coastline', check_size=False)
        if this_cst is not None:
            with utils.ccp(
                    total=len(this_cst.results),
                    desc='Processing osm coastline',
                    leave=True
            ) as pbar:
                for n, cst_result in enumerate(this_cst.results):                
                    if cst_result[-1] == 0:
                        pbar.update()
                        cst_osm = cst_result[1]
                        out = osm.polygonize_osm_coastline(
                            cst_osm, utils.make_temp_fn(
                                utils.fn_basename2(cst_osm) + '_coast.shp',
                                temp_dir=self.cache_dir
                            ),
                            region=self.wgs_region,
                            include_landmask=False,
                            landmask_is_watermask=True,
                            line_buffer=0.0000001
                        )

                        cst_ds = gdalfun.gdal_mem_ds(
                            self.mem_config, name='osm_coast', src_srs=self.wgs_srs, co=self.co
                        )
                        cst_ds_ = ogr.Open(out)
                        if cst_ds_ is not None:
                            cst_layer = cst_ds_.GetLayer()
                            gdal.RasterizeLayer(cst_ds, [1], cst_layer, burn_values=[1])

                            for srcwin in utils.yield_srcwin(
                                    (self.ds_config['ny'], self.ds_config['nx']),
                                    1000, verbose=self.verbose
                            ):
                                cst_ds_arr = cst_ds.GetRasterBand(1).ReadAsArray(*srcwin)
                                self.coast_array[
                                    srcwin[1]:srcwin[1]+srcwin[3],
                                    srcwin[0]:srcwin[0]+srcwin[2]
                                ] += cst_ds_arr
                            cst_ds = cst_ds_ = None
                        else:
                            utils.echo_error_msg(f'could not open {out}')

                            
    def _load_enc_coast(self):
        """Load ENC Charts Coastline."""
        
        this_cst = self.fetch_data('charts')
        if this_cst is not None:
            with utils.ccp(
                    total=len(this_cst.results),
                    desc='Processing charts coastline',
                    leave=True
            ) as pbar:
                for n, cst_result in enumerate(this_cst.results):
                    if cst_result[-1] == 0:
                        pbar.update()
                        cst_enc = cst_result[1]
                        src_000s = utils.p_unzip(
                            os.path.join(cst_enc), exts=['000'],
                            outdir=self.cache_dir, verbose=self.verbose            
                        )
                        for src_000 in src_000s:
                            enc_level = utils.int_or(os.path.basename(src_000)[2], 0)
                            if enc_level <= 3:
                                continue

                            enc_ds = gdalfun.gdal_mem_ds(
                                self.mem_config, name='enc', src_srs=self.wgs_srs, co=self.co
                            )
                            cst_ds_ = ogr.Open(src_000)
                            if cst_ds_ is not None:
                                cst_layer = cst_ds_.GetLayer('SEAARE')                            
                                _ = gdal.RasterizeLayer(
                                    enc_ds, [1], cst_layer, burn_values=[1]
                                )
                                cst_ds_ = cst_layer = None
                                
                            if enc_ds is not None:
                                for srcwin in utils.yield_srcwin(
                                        (self.ds_config['ny'], self.ds_config['nx']),
                                        1000, verbose=self.verbose
                                ):
                                    enc_ds_arr = enc_ds.GetRasterBand(1).ReadAsArray(*srcwin)
                                    self.coast_array[
                                        srcwin[1]:srcwin[1]+srcwin[3],
                                        srcwin[0]:srcwin[0]+srcwin[2]
                                    ] -= (enc_ds_arr * self.min_weight)
                                    
                                enc_ds = enc_ds_arr = None

                                
    def _load_gmrt(self):
        """GMRT - Global low-res."""
        
        this_gmrt = self.fetch_data('gmrt:layer=topo')        
        gmrt_result = this_gmrt.results[0]
        if gmrt_result[-1] == 0:
            gmrt_tif = gmrt_result[1]
            gmrt_ds = gdalfun.gdal_mem_ds(
                self.ds_config, name='gmrt', src_srs=self.wgs_srs, co=self.co
            )
            gdal.Warp(
                gmrt_ds, gmrt_tif, dstSRS=self.cst_srs, resampleAlg=self.sample
            )
            gmrt_ds_arr = gmrt_ds.GetRasterBand(1).ReadAsArray()
            gmrt_ds_arr[gmrt_ds_arr > 0] = 1
            gmrt_ds_arr[gmrt_ds_arr < 0] = 0
            self.coast_array += (gmrt_ds_arr * self.min_weight)
            gmrt_ds = gmrt_ds_arr = None

            
    def _load_copernicus(self):
        """Copernicus Global DEM."""

        this_cop = self.fetch_data('copernicus:datatype=1', check_size=False)
        for i, cop_result in enumerate(this_cop.results):
            if cop_result[-1] == 0:
                cop_tif = cop_result[1]
                gdalfun.gdal_set_ndv(cop_tif, 0, verbose=False)
                cop_ds = gdalfun.gdal_mem_ds(
                    self.ds_config, name='copernicus', src_srs=self.wgs_srs, co=self.co
                )
                gdal.Warp(
                    cop_ds, cop_tif, dstSRS=self.cst_srs,
                    resampleAlg=self.sample, callback=False, srcNodata=0
                )
                for srcwin in utils.yield_srcwin(
                        (self.ds_config['ny'], self.ds_config['nx']),
                        1000, verbose=self.verbose
                ):
                    cop_ds_arr = cop_ds.GetRasterBand(1).ReadAsArray(*srcwin)
                    cop_ds_arr[cop_ds_arr != 0] = 1
                    self.coast_array[
                        srcwin[1]:srcwin[1]+srcwin[3],
                        srcwin[0]:srcwin[0]+srcwin[2]
                    ] += (cop_ds_arr * self.min_weight)
                    
                cop_ds = cop_ds_arr = None

                
    def _load_cudem(self):
        """CUDEM (Standard)."""

        this_cudem = self.fetch_data('CUDEM:datatype=ninth', check_size=False)
        for i, cudem_result in enumerate(this_cudem.results):
            if cudem_result[-1] == 0:
                cudem_tif = cudem_result[1]
                cudem_ds = gdalfun.gdal_mem_ds(
                    self.mem_config, name='cudem', src_srs=self.wgs_srs, co=self.co
                )
                gdal.Warp(
                    cudem_ds, cudem_tif, dstSRS=self.cst_srs,
                    resampleAlg=self.sample, callback=False, srcNodata=0
                )
                for srcwin in utils.yield_srcwin(
                        (self.ds_config['ny'], self.ds_config['nx']),
                        1000, verbose=self.verbose
                ):
                    cudem_ds_arr = cudem_ds.GetRasterBand(1).ReadAsArray(*srcwin)
                    cudem_ds_arr[cudem_ds_arr != 0] = 1
                    self.coast_array[
                        srcwin[1]:srcwin[1]+srcwin[3],
                        srcwin[0]:srcwin[0]+srcwin[2]
                    ] += (cudem_ds_arr * self.min_weight)
                    
                cudem_ds = cudem_ds_arr = None

                
    def _load_nhd(self):
        """USGS NHD (High-Res US Only)."""

        this_tnm = self.fetch_data(
            "tnm:datasets=13/14:extents='HU-8 Subbasin,HU-4 Subregion'"
        )
        if len(this_tnm.results) > 0:
            tnm_ds = gdalfun.gdal_mem_ds(
                self.ds_config, name='nhd', src_srs=self.wgs_srs, co=self.co
            )
            for i, tnm_result in enumerate(this_tnm.results):
                if tnm_result[-1] == 0:
                    tnm_zip = tnm_result[1]
                    if not os.path.exists(tnm_zip):
                        break

                    # Unzip logic (assumes gdb structure)
                    utils.unzip(tnm_zip, outdir=os.path.dirname(tnm_zip), verbose=False)
                    gdb = '.'.join(tnm_zip.split('.')[:-1]) + '.gdb'
                    if 'GDB' not in gdb:
                        continue
                    
                    # Merge relevant layers into temp shapefile
                    utils.run_cmd(
                        f'ogr2ogr -update -append nhdArea_merge.shp "{gdb}" NHDArea -where "FType=312 OR FType=336 OR FType=445 OR FType=460 OR FType=537" -clipdst {self.p_region.format("ul_lr")} 2>/dev/null',
                        verbose=self.verbose
                    )
                    utils.run_cmd(
                        f'ogr2ogr -update -append nhdArea_merge.shp "{gdb}" NHDWaterbody -where "FType=493 OR FType=466" -clipdst {self.p_region.format("ul_lr")} 2>/dev/null',
                        verbose=self.verbose
                    )

                    if self.want_nhd_plus:
                        utils.run_cmd(
                            f'ogr2ogr -update -append nhdArea_merge.shp "{gdb}" NHDPlusBurnWaterbody -clipdst {self.p_region.format("ul_lr")} 2>/dev/null',
                            verbose=self.verbose
                        )

                try:
                    utils.run_cmd(
                        f'gdal_rasterize -burn 1 nhdArea_merge.shp nhdArea_merge.tif -te {self.p_region.format("te")} -ts {self.ds_config["nx"]} {self.ds_config["ny"]} -ot Int32',
                        verbose=self.verbose
                    )

                    gdal.Warp(
                        tnm_ds, 'nhdArea_merge.tif',
                        dstSRS=self.cst_srs, resampleAlg=self.sample
                    )
                except Exception:
                    tnm_ds = None

                if tnm_ds is not None:
                    tnm_ds_arr = tnm_ds.GetRasterBand(1).ReadAsArray()
                    tnm_ds_arr[tnm_ds_arr < 1] = 0
                    self.coast_array -= (tnm_ds_arr * self.min_weight)
                    tnm_ds = tnm_ds_arr = None

                utils.remove_glob('nhdArea_merge.*')

                
    def _load_lakes(self):
        """HydroLakes -- Global Lakes."""
        
        this_lakes = self.fetch_data('hydrolakes')
        lakes_shp = None
        if this_lakes.results and this_lakes.results[0][-1] == 0:
            lakes_zip = this_lakes.results[0][1]
            lakes_shps = utils.unzip(lakes_zip, self.cache_dir)
            for i in lakes_shps:
                if i.split('.')[-1] == 'shp':
                    lakes_shp = i

            lakes_ds = gdalfun.gdal_mem_ds(
                self.ds_config, name='lakes', src_srs=self.wgs_srs, co=self.co
            )
            lakes_warp_ds = gdalfun.gdal_mem_ds(
                self.ds_config, name='lakes_warp', src_srs=self.wgs_srs, co=self.co
            )
            
            lk_ds = ogr.Open(lakes_shp)
            if lk_ds is not None:
                lk_layer = lk_ds.GetLayer()
                lk_layer.SetSpatialFilter(self.f_region.export_as_geom())
                gdal.RasterizeLayer(lakes_ds, [1], lk_layer, burn_values=[-1])
                gdal.Warp(
                    lakes_warp_ds, lakes_ds,
                    dstSRS=self.cst_srs, resampleAlg=self.sample
                )
                lakes_ds_arr = lakes_warp_ds.GetRasterBand(1).ReadAsArray()
                
                self.coast_array[lakes_ds_arr == -1] = 0 if not self.invert_lakes else 1 
                lakes_ds = lk_ds = lakes_warp_ds = None
            else:
                utils.echo_error_msg(f'could not open {lakes_shp}')

                
    def _load_wsf_bldgs(self):
        """WSF Buildings (World Settlement Footprint)."""

        this_wsf = self.fetch_data('wsf', check_size=False)        
        for i, wsf_result in enumerate(this_wsf.results):
            if wsf_result[-1] == 0:
                wsf_tif = wsf_result[1]
                gdalfun.gdal_set_ndv(wsf_tif, 0, verbose=False)
                wsf_ds = gdalfun.gdal_mem_ds(
                    self.ds_config, name='wsf', src_srs=self.wgs_srs, co=self.co
                )
                gdal.Warp(
                    wsf_ds, wsf_tif, dstSRS=self.cst_srs,
                    resampleAlg='cubicspline', callback=gdal.TermProgress,
                    srcNodata=0, outputType=gdal.GDT_Float32
                )
                wsf_ds_arr = wsf_ds.GetRasterBand(1).ReadAsArray()
                wsf_ds_arr[wsf_ds_arr != 0 ] = -1
                self.coast_array += (wsf_ds_arr * self.min_weight)
                wsf_ds = wsf_ds_arr = None

                
    def _load_osm_bldgs(self):
        """Load buildings from OSM."""

        this_osm = self.fetch_data(
            "osm:q=buildings:fmt=osm:chunks=True:min_length={}:planet={}".format(
                self.min_building_length, self.want_osm_planet
            )
        )
        os.environ["OGR_OSM_OPTIONS"] = "OGR_INTERLEAVED_READING=YES"
        
        with utils.ccp(
                total=len(this_osm.results),
                desc='Processing OSM buildings',
                leave=self.verbose
        ) as pbar:
            for n, osm_result in enumerate(this_osm.results):
                if osm_result[-1] == 0:
                    pbar.update()
                    osm_z = osm_result[1]
                    
                    if osm_result[-1] == 'bz2':
                        osm_planet = utils.unbz2(osm_z, self.cache_dir)
                        osm_file = utils.ogr_clip(osm_planet, self.wgs_region)
                    elif osm_result[-1] == 'pbf':
                        osm_file = utils.ogr_clip(osm_z, self.wgs_region, 'multipolygons')
                    else:
                        osm_file = osm_z

                    if os.path.getsize(osm_file) == 366:
                        continue

                    out, status = utils.run_cmd(
                        (f'gdal_rasterize -burn -1 -l multipolygons "{osm_file}" '
                         f'bldg_osm.tif -te {self.p_region.format("te")} '
                         f'-ts {self.ds_config["nx"]} '
                         f'{self.ds_config["ny"]} -ot Int32 -q'),
                        verbose=False
                    )

                    if status == 0:
                        bldg_ds = gdal.Open('bldg_osm.tif')
                        if bldg_ds is not None:
                            bldg_ds_arr = bldg_ds.GetRasterBand(1).ReadAsArray()
                            self.coast_array[bldg_ds_arr == -1] = 0 
                            bldg_ds = bldg_ds_arr = None
                        bldg_ds = None
                    else:
                        utils.echo_warning_msg('could not parse buildings!')

                    utils.remove_glob('bldg_osm.tif*')

        bldg_ds = None

        
    def _load_bing_bldgs(self):
        """Load buildings from BING."""

        this_bing = self.fetch_data("bing_bfp", check_size=True)
        os.environ["OGR_OSM_OPTIONS"] = "OGR_INTERLEAVED_READING=YES"
        
        with utils.ccp(
                total=len(this_bing.results),
                desc='Processing BING buildings',
                leave=self.verbose
        ) as pbar:
            for n, bing_result in enumerate(this_bing.results):                
                if bing_result[-1] == 0:
                    pbar.update()
                    bing_gz = bing_result[1]
                    try:
                        bing_gj = utils.gunzip(bing_gz, self.cache_dir)
                        os.rename(bing_gj, bing_gj + '.geojson')
                        bing_gj = bing_gj + '.geojson'

                        out, status = utils.run_cmd(
                            ('gdal_rasterize -burn -1 -l '
                             f'{os.path.basename(utils.fn_basename2(bing_gj))} '
                             f'{bing_gj} bldg_bing.tif '
                             f'-te {self.p_region.format("te")} '
                             f'-ts {self.ds_config["nx"]} '
                             f'{self.ds_config["ny"],} -ot Int32'),
                            verbose=True
                        )
                    except Exception as e:
                        utils.remove_glob(bing_gz)
                        utils.echo_error_msg(f'Could not process bing bfp, {e}')
                        status = -1
                    
                    if status == 0:
                        bldg_ds = gdal.Open('bldg_bing.tif')
                        if bldg_ds is not None:
                            bldg_ds_arr = bldg_ds.GetRasterBand(1).ReadAsArray()
                            self.coast_array[bldg_ds_arr == -1] = 0 
                            bldg_ds = bldg_ds_arr = None
                        bldg_ds = None
                    else:
                        utils.echo_warning_msg('could not parse buildings!')

                    utils.remove_glob('bldg_bing.tif*')

        bldg_ds = None

        
    def _load_data(self):
        """Load data from user datalist and amend coast_array."""
        
        for this_arr in self.stack_ds.yield_array():
            data_arr = this_arr[0]['z']
            weight_arr = this_arr[0]['weight']
            weight_arr[np.isnan(weight_arr)] = 0
            srcwin = this_arr[1]
            data_arr[np.isnan(data_arr)] = 0
            data_arr[data_arr > 0] = 1
            data_arr[data_arr < 0] = -1
            self.coast_array[
                srcwin[1]:srcwin[1]+srcwin[3],
                srcwin[0]:srcwin[0]+srcwin[2]
            ] += (data_arr * weight_arr)

            
    def _write_coast_array(self):
        """Write coast_array to file."""

        if self.verbose:
            utils.echo_msg(f'writing array to {self.name}.tif...')
            
        gdalfun.gdal_write(
            self.coast_array, f'{self.name}.tif', self.ds_config
        )

        
    def _write_coast_poly(self, poly_count=None):
        """Convert to coast_array vector."""

        if self.verbose:
            utils.echo_msg(
                f'polygonizing {poly_count if poly_count else "all"} features from array to {self.name}.shp...'
            )
        
        poly_count = utils.int_or(poly_count)
        tmp_ds = ogr.GetDriverByName('ESRI Shapefile').CreateDataSource(
            '{}_tmp_c.shp'.format(self.name)
        )
        
        if tmp_ds is not None:
            tmp_layer = tmp_ds.CreateLayer(
                '{}_tmp_c'.format(os.path.basename(self.name)),
                None,
                ogr.wkbMultiPolygon
            )
            tmp_layer.CreateField(ogr.FieldDefn('DN', ogr.OFTInteger))
            gdalfun.gdal_polygonize(
                '{}.tif'.format(self.name), tmp_layer, verbose=self.verbose
            )
            tmp_ds = None
            
        # SQL filtering for largest polygons
        sql_limit = f" order by ST_AREA(geometry) desc limit {poly_count}" if poly_count is not None else ""
        
        utils.run_cmd(
            f'ogr2ogr -dialect SQLITE -sql "SELECT * FROM {os.path.basename(self.name)}_tmp_c WHERE DN=0{sql_limit}" "{self.name}.shp" "{self.name}_tmp_c.shp"',
            verbose=self.verbose
        )        
        utils.remove_glob(f'{self.name}_tmp_c.*')
        utils.run_cmd(
            f'ogrinfo -dialect SQLITE -sql "UPDATE {os.path.basename(self.name)} SET geometry = ST_MakeValid(geometry)" "{self.name}.shp"',
            verbose=self.verbose
        )        
        gdalfun.osr_prj_file(self.name + '.prj', self.dst_srs)


### End
