### fetchers.py 
##
## Copyright (c) 2010 - 2026 Regents of the University of Colorado
##
## fetchers.py is part of CUDEM
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
## Data Fetcher Modules (Remote Data Sources)
##
### Code:

import os
import json
import copy
from osgeo import ogr, gdal

from cudem import utils
from cudem import regions
from cudem import gdalfun
from cudem import srsfun
from cudem import xyzfun
from cudem import vdatums
from cudem import fetches
from cudem.datalists.dlim import ElevationDataset

## ==============================================
## dlim Fetcher classes
## ==============================================
class Fetcher(ElevationDataset):
    """The default fetches dataset type.
    Used for on-the-fly remote data parsing and processing.
    """
    
    def __init__(self,
                 keep_fetched_data=True,
                 mask_coast=False,
                 invert_coast=True,
                 coast_buffer=0.00001,
                 outdir=None,
                 check_size=True,
                 want_single_metadata_name=False,
                 callback=fetches.fetches.fetches_callback,
                 **kwargs):
        
        super().__init__(**kwargs)
        
        ## Setup Region for WGS84 (Remote APIs usually expect Lat/Lon)
        #if self.region:
        self.wgs_region = self.region.copy()
        #else:
        #    self.wgs_region = regions.Region(xmin=-180, xmax=180, ymin=-90, ymax=90)
            
        self.wgs_region.src_srs = self.dst_srs
        self.wgs_srs = 'epsg:4326'
        if self.dst_srs is not None:
            self.wgs_region.warp(self.wgs_srs)

        if outdir is None:
            outdir = self.cache_dir
            
        ## Initialize the underlying Fetch Module (from cudem.fetches)
        self.fetch_module = fetches.fetches.FetchesFactory(
            mod=self.fn, src_region=self.wgs_region,
            callback=callback, verbose=False,
            outdir=outdir,
        )._acquire_module()

        if self.fetch_module is None:
            utils.echo_warning_msg(f'Fetch module {self.fn} returned None')

        self.mask_coast = mask_coast
        self.invert_coast = invert_coast
        self.coast_buffer = utils.float_or(coast_buffer, 0.00001)
        self.want_single_metadata_name = want_single_metadata_name
        self.check_size = check_size
        self.keep_fetched_data = keep_fetched_data 
        
        if self.fetch_module:
            self.outdir = os.path.abspath(self.fetch_module._outdir)
            if not os.path.exists(self.outdir):
                os.makedirs(self.outdir)
        else:
            self.outdir = outdir

        self._reset_params()

        
    def _reset_params(self):
        """Prepare parameters for child datasets."""
        
        if not self.fetch_module: return

        md = copy.deepcopy(self.metadata)
        md.update({
            'title': self.fetch_module.title,
            'source': self.fetch_module.source,
            'date': self.fetch_module.date,
            'resolution': self.fetch_module.resolution,
            'hdatum': self.fetch_module.hdatum,
            'vdatum': self.fetch_module.vdatum,
            'url': self.fetch_module.url
        })
        
        self.fetches_params = self._set_params(
            data_format=self.fetch_module.data_format,
            src_srs=self.fetch_module.src_srs,
            cache_dir=self.fetch_module._outdir,
            remote=True,
            metadata=md,
        )

        ## Coastline Masking
        if self.mask is None and self.mask_coast:
            from cudem.fetches import osm
            
            ## Fetch OSM Coastline
            this_osm = osm.osmCoastline(
                region=self.region, chunks=True, verbose=False, attempts=5,
                cache_dir=self.cache_dir, landmask_is_watermask=True
            )
            coast_mask = this_osm(return_geom=False, overwrite=True)

            # ## Fetch OSM Water
            # this_osm_water = osm.osmCoastline(
            #     region=self.region, chunks=True, verbose=False, attempts=5, cache_dir=self.cache_dir, q='water'
            # )
            # water_mask = this_osm_water(return_geom=False, overwrite=False)
            
            masks = []
            if coast_mask:
                masks.append(f'mask_fn={coast_mask}:invert={self.invert_coast}')
            # if water_mask:
            #     masks.append(f'mask_fn={water_mask}:invert={self.invert_coast}')
                
            if masks:
                self.fetches_params['mask'] = masks

                
    def generate_inf(self, make_grid=True, make_block_mean=False, block_inc=None):
        """Generate metadata (INF) from the Fetcher results."""
        
        # Note: We usually generate INF based on the Search Region because
        # actually fetching all data to determine exact bounds is expensive/slow.
        
        if self.fetch_module:
            tmp_region = self.fetch_module.region if self.region is None else self.region.copy()
            self.infos.minmax = tmp_region.export_as_list()    
            self.infos.wkt = tmp_region.export_as_wkt()
            self.infos.src_srs = self.fetch_module.src_srs
        utils.echo_msg(self.infos)
        return self.infos

    
    def parse(self):
        """Run the fetcher, iterate results, and yield Datasets."""
        
        if not self.fetch_module: return

        ## Run the Fetcher
        if len(self.fetch_module.results) == 0:
            try:
                self.fetch_module.run()
            except Exception as e:
                utils.echo_warning_msg(f'Fetch module {self.fn} failed: {e}')
                return

        ## Iterate Results
        with utils.ccp(total=len(self.fetch_module.results), 
                       desc=f'Parsing {self.fetch_module.name}', leave=False) as pbar:
            
            for result in self.fetch_module.results:
                ## Download File
                status = self.fetch_module.fetch_entry(result, check_size=self.check_size)
                
                if status == 0:
                    ## Update Mod Param to local path
                    self._reset_params()
                    self.fetches_params['mod'] = os.path.join(self.fetch_module._outdir, result['dst_fn'])

                    ## Yield Child Datasets
                    for this_ds in self.yield_ds(result):
                        if this_ds is not None:
                            utils.echo_debug_msg(f'{this_ds}: {self.fetches_params}')
                            ## Update Name/Metadata
                            f_name = os.path.basename(this_ds.fn)
                            mod_name = self.fetch_module.name
                            
                            this_ds.metadata['name'] = f"{this_ds.metadata.get('name', '')}/{f_name}"
                            this_ds.initialize()
                            
                            ## Recursively Parse Child
                            for ds in this_ds.parse():
                                self.data_entries.append(ds)
                                yield ds
                        else:
                            utils.echo_warning_msg(f'Could not configure dataset for {result}')
                else:
                    utils.echo_warning_msg(f'Failed to fetch {result}: Status {status}')
                
                pbar.update()
                
        ## Cleanup
        if not self.keep_fetched_data:
            utils.remove_glob(f'{self.fn}*')

            
    def yield_ds(self, result):
        """Factory method to generate ElevationDataset from a fetch result.
        Subclasses should override or extend this.
        """
        
        from cudem.datalists.dlim import DatasetFactory
        
        ## Try to detect SRS from file if GDAL supported
        try:
            vdatum = self.fetch_module.vdatum
            local_fn = os.path.join(self.fetch_module._outdir, result['dst_fn'])
            src_srs = gdalfun.gdal_get_srs(local_fn)
            
            if src_srs:
                horz_epsg, vert_epsg = gdalfun.epsg_from_input(src_srs)
                if vert_epsg is None: vert_epsg = vdatum
                
                if vert_epsg is not None:
                    self.fetch_module.src_srs = f'{horz_epsg}+{vert_epsg}'
                else:
                    self.fetch_module.src_srs = src_srs
        except Exception:
            pass

        yield DatasetFactory(**self.fetches_params)._acquire_module()


class NEDFetcher(Fetcher):
    """National Elevation Dataset (USGS) Fetcher."""

    __doc__ = f"{__doc__}\nFetches Module: <ned> - {fetches.tnm.NED.__doc__}"
    
    def __init__(self, coast_buffer=0.00001, remove_flat=True, **kwargs):
        super().__init__(**kwargs)
        self.remove_flat = remove_flat
        self.coast_buffer = utils.float_or(coast_buffer, 0.00001)

        
    def yield_ds(self, result):
        from cudem.datalists.dlim import DatasetFactory
        
        src_dem = os.path.join(self.fetch_module._outdir, result['dst_fn'])
        self.fetches_params['mod'] = src_dem
        self.fetches_params['mask'] = self.mask
        self.fetches_params['remove_flat'] = self.remove_flat
        yield DatasetFactory(**self.fetches_params)._acquire_module()


class DNRFetcher(Fetcher):
    """WA DNR Fetcher."""
    
    __doc__ = f"{__doc__}\nFetches Module: <wadnr> - {fetches.wadnr.WADNR.__doc__}"
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.data_format = -2

    def yield_ds(self, result):
        from cudem.datalists.dlim import DatasetFactory
        
        src_dnr_dems = utils.p_unzip(
            os.path.join(self.fetch_module._outdir, result['dst_fn']),
            exts=['tif'], outdir=self.fetch_module._outdir, verbose=self.verbose
        )
        for src_dnr_dem in src_dnr_dems:
            self.fetches_params['mod'] = src_dnr_dem
            self.fetches_params['remove_flat'] = True
            yield DatasetFactory(**self.fetches_params)._acquire_module()


class DAVFetcher_CoNED(Fetcher):
    """Digital Coast CoNED Fetcher."""
    
    __doc__ = f"{__doc__}\nFetches Module: <CoNED> - {fetches.dav.CoNED.__doc__}"
    
    def __init__(self, keep_fetched_data=True, cog=True, **kwargs):
        super().__init__(**kwargs)
        self.keep_fetched_data = keep_fetched_data
        self.cog = cog

    def parse(self):
        ## Override parse to support COG direct access (no fetch needed)
        for result in self.fetch_module.results:
            if not self.cog:
                status = self.fetch_module.fetch_entry(result, check_size=self.check_size)
                if status != 0: continue
                
            for this_ds in self.yield_ds(result):
                if this_ds is not None:
                    this_ds.initialize()
                    yield this_ds

                    
    def yield_ds(self, result):
        from cudem.datalists.dlim import DatasetFactory
        
        ## SRS Detection
        try:
            vdatum = self.fetch_module.vdatum
            if not self.cog:
                src_srs = gdalfun.gdal_get_srs(os.path.join(self.fetch_module._outdir, result['dst_fn']))
            else:
                src_srs = gdalfun.gdal_get_srs(result['url'])
                
            horz_epsg, vert_epsg = gdalfun.epsg_from_input(src_srs)
            if vert_epsg is None: vert_epsg = vdatum

            if vert_epsg is not None:
                self.fetch_module.src_srs = f'{horz_epsg}+{vert_epsg}'
            else:
                self.fetch_module.src_srs = src_srs
        except: pass

        self.fetches_params['mod'] = os.path.join(self.fetch_module._outdir, result['dst_fn']) if not self.cog else result['url']
        self.fetches_params['check_path'] = not self.cog
        self.fetches_params['src_srs'] = self.fetch_module.src_srs
        self.fetches_params['data_format'] = 200
        
        yield DatasetFactory(**self.fetches_params)._acquire_module()


class DAVFetcher_SLR(Fetcher):
    """Digital Coast SLR Fetcher."""
    
    __doc__ = f"{__doc__}Fetches Module: <SLR> - {fetches.dav.SLR.__doc__}"
    
    def yield_ds(self, result):
        from cudem.datalists.dlim import DatasetFactory
        self.fetches_params['remove_flat'] = True
        yield DatasetFactory(**self.fetches_params)._acquire_module()    


class SWOTFetcher(Fetcher):
    """NASA SWOT Data Fetcher."""
    
    __doc__ = f"{__doc__}\nFetches Module: <swot> - {fetches.earthdata.SWOT.__doc__}"

    def __init__(self,
                 data_set='wse',
                 apply_geoid=True,
                 classes=None,
                 classes_qual=None,
                 anc_classes=None,
                 remove_class_flags=False,
                 **kwargs):
        super().__init__(**kwargs)
        self.fetches_params.update({
            'var': data_set,
            'apply_geoid': apply_geoid,
            'classes': classes,
            'anc_classes': anc_classes,
            'classes_qual': classes_qual,
            'remove_class_flags': remove_class_flags
        })

        
    def yield_ds(self, result):
        from cudem.datalists.dlim import DatasetFactory
        
        if 'L2_HR_PIXC_' in result['data_type']:
            self.fetches_params['data_format'] = 202
            yield DatasetFactory(**self.fetches_params)._acquire_module()
        elif 'L2_HR_Raster' in result['data_type']:
            self.fetches_params['data_format'] = 203
            yield DatasetFactory(**self.fetches_params)._acquire_module()
        else:
            utils.echo_warning_msg(f'{result["data_type"]} is not a supported SWOT dataset type')


class IceSat2Fetcher(Fetcher):
    """NASA ICESat-2 Data Fetcher.
    Supports ATL03 (Photons) and ATL24 (Bathymetry).
    """
    
    from cudem.datalists.icesat2file import IceSat2_ATL03
    __doc__ = f"{__doc__}\nDLIM Module: <IceSat2File> - {IceSat2_ATL03.__doc__}"
    __doc__ += f"\nFetches Module: <icesat2> - {fetches.earthdata.IceSat2.__doc__}"

    def __init__(self,
                 water_surface=None,
                 classes=None,
                 confidence_levels=None,
                 columns={},
                 classify_bathymetry=False,
                 classify_buildings=True,
                 classify_water=True,
                 classify_inland_water=True,
                 reject_failed_qa=True,
                 min_bathy_confidence=None,
                 dataset='atl03',
                 **kwargs):
        
        super().__init__(**kwargs)
        self.water_surface = water_surface
        self.classes = classes
        self.confidence_levels = confidence_levels
        self.columns = columns
        
        self.classify_bathymetry = classify_bathymetry
        self.classify_buildings = classify_buildings
        self.classify_water = classify_water
        self.classify_inland_water = classify_inland_water
        self.reject_failed_qa = reject_failed_qa
        self.min_bathy_confidence = min_bathy_confidence
        self.data_format = -111

        from cudem.fetches import bingbfp
        from cudem.fetches import osm
        
        ## Prepare Masks (fetch once here so they can be passed to all granules)
        if self.classify_buildings:
            bfp = bingbfp.bingBuildings(region=self.region, verbose=self.verbose, cache_dir=self.cache_dir)
            self.fetches_params['classify_buildings'] = bfp(return_geom=True)

        if self.classify_water:
            osmc = osm.osmCoastline(region=self.region, chunks=True, verbose=self.verbose, attempts=5, cache_dir=self.cache_dir)
            coast_mask = osmc(return_geom=True)
            self.fetches_params['classify_water'] = coast_mask
            
        if self.classify_inland_water:
            osml = osm.osmCoastline(region=self.region, chunks=True, verbose=self.verbose, attempts=5, cache_dir=self.cache_dir, q='water')
            lakes_mask = osml(return_geom=True)
            self.fetches_params['classify_inland_water'] = lakes_mask

            
    def yield_ds(self, result):
        from cudem.datalists.dlim import DatasetFactory
        from cudem.datalists.icesat2file import IceSat2_ATL03
        
        icesat2_fn = os.path.join(self.fetch_module._outdir, result['dst_fn'])

        if result['data_type'].lower() == 'atl03' or 'NSIDC_CPRD' in result['data_type'].upper():
            ## Setup ATL03 Parameters
            self.fetches_params.update({
                'water_surface': self.water_surface or 'geoid',
                'classes': self.classes,
                'confidence_levels': self.confidence_levels,
                'columns': self.columns,
                'reject_failed_qa': self.reject_failed_qa,
                'min_bathy_confidence': self.min_bathy_confidence,
                'classify_bathymetry': self.classify_bathymetry
                # classify_buildings/water already in fetches_params from init
            })

            if 'processed_zip' in result['data_type']:
                ## Handle ZIPs from NSIDC
                icesat2_h5s = utils.p_unzip(icesat2_fn, exts=['h5'], outdir=self.cache_dir, verbose=self.verbose)
                for h5_fn in icesat2_h5s:
                    self.fetches_params['fn'] = h5_fn
                    yield IceSat2_ATL03(**self.fetches_params)._acquire_module()
            else:
                self.fetches_params['fn'] = icesat2_fn
                self.fetches_params['data_format'] = 303
                yield DatasetFactory(**self.fetches_params)._acquire_module()
                
        elif result['data_type'].lower() == 'atl24':
            ## Setup ATL24 Parameters
            self.fetches_params.update({
                'water_surface': self.water_surface or 'ortho',
                'classes': self.classes,
                'min_confidence': self.min_bathy_confidence,
                'fn': icesat2_fn,
                'data_format': 304
            })
            yield DatasetFactory(**self.fetches_params)._acquire_module()
            
        else:
            utils.echo_warning_msg(f'{icesat2_fn} cannot be processed (Unknown Type)')


class GMRTFetcher(Fetcher):
    """GMRT Gridded data Fetcher."""
    
    __doc__ = f"{__doc__}\nFetches Module: <gmrt> - {fetches.gmrt.GMRT.__doc__}"
    
    def __init__(self, swath_only=False, **kwargs):
        super().__init__(**kwargs)
        self.swath_only = swath_only

    def yield_ds(self, result):
        from cudem.datalists.dlim import DatasetFactory
        
        gmrt_fn = os.path.join(self.fetch_module._outdir, result['dst_fn'])
        
        ## Prepare GMRT file
        with gdalfun.gdal_datasource(gmrt_fn, update=1) as src_ds:
            md = src_ds.GetMetadata()
            md['AREA_OR_POINT'] = 'Point'
            src_ds.SetMetadata(md)
            gdalfun.gdal_set_srs(src_ds)
            gdalfun.gdal_set_ndv(src_ds, verbose=False)

        ## Swath Masking
        if self.swath_only:
            # Fetch mask if needed
            if fetches.fetches.Fetch(self.fetch_module._gmrt_swath_poly_url, verbose=self.verbose).fetch_file(
                os.path.join(self.fetch_module._outdir, 'gmrt_swath_polygons.zip')
            ) == 0:
                swath_shps = utils.p_unzip(
                    os.path.join(self.fetch_module._outdir, 'gmrt_swath_polygons.zip'),
                    exts=['shp', 'shx', 'prj', 'dbf'], outdir=self.cache_dir, verbose=self.verbose
                )
                
                swath_mask = next((f for f in swath_shps if f.endswith('.shp')), None)
                if swath_mask and os.path.exists(swath_mask):
                    self.fetches_params['mask'] = {'mask': swath_mask, 'invert_mask': True}
                else:
                    self.swath_only = False

        yield DatasetFactory(**self.fetches_params)._acquire_module()


class GEBCOFetcher(Fetcher):
    """GEBCO Gridded data Fetcher."""
    
    __doc__ = f"{__doc__}\nFetches Module: <gebco> - {fetches.gebco.GEBCO.__doc__}"

    def __init__(self, exclude_tid=None, **kwargs):
        super().__init__(**kwargs)
        self.exclude_tid = []
        if exclude_tid:
            self.exclude_tid = [utils.int_or(x) for x in str(exclude_tid).split('/')]

            
    def yield_ds(self, result):
        from cudem.datalists.dlim import DatasetFactory
        
        wanted_gebco_fns = []
        gebco_fns = utils.p_unzip(
            os.path.join(self.fetch_module._outdir, result['dst_fn']),
            exts=['tif'], outdir=self.fetch_module._outdir, verbose=self.verbose
        )
        
        ## TID (Type Identifier) Masking
        if self.exclude_tid:
            ## Fetch TID Grid
            if fetches.fetches.Fetch(self.fetch_module._gebco_urls['gebco_tid']['geotiff'], verbose=self.verbose).fetch_file(
                os.path.join(self.fetch_module._outdir, 'gebco_tid.zip')
            ) == 0:
                tid_fns = utils.p_unzip(
                    os.path.join(self.fetch_module._outdir, 'gebco_tid.zip'),
                    exts=['tif'], outdir=self.cache_dir, verbose=self.verbose
                )
                
                ## Filter TID files
                for tid_fn in tid_fns:
                    ## Spatial Check
                    ds_config = gdalfun.gdal_infos(tid_fn)
                    inf_region = regions.Region().from_geo_transform(ds_config['geoT'], ds_config['nx'], ds_config['ny'])
                    
                    if self.region is None or self.region.valid_p(check_xy=True) and regions.regions_intersect_p(inf_region, self.region):
                        wanted_gebco_fns.append(tid_fn)

                ## Process Masking
                for tid_fn in wanted_gebco_fns:
                    tmp_tid = utils.make_temp_fn('tmp_tid.tif', temp_dir=self.fetch_module._outdir)
                    
                    ## Create Mask
                    with gdalfun.gdal_datasource(tid_fn) as tid_ds:
                        tid_band = tid_ds.GetRasterBand(1)
                        tid_array = tid_band.ReadAsArray().astype(float)
                        
                        ## Apply Exclusions
                        for tid_key in self.exclude_tid:
                            tid_array[tid_array == tid_key] = -9999
                        
                        ## Write Temp Mask
                        gdalfun.gdal_write(tid_array, tmp_tid, {'geoT': tid_ds.GetGeoTransform(), 'proj': tid_ds.GetProjection()})
                        
                    self.fetches_params['mod'] = tid_fn.replace('tid_', '')
                    self.fetches_params['data_format'] = 200
                    self.fetches_params['mask'] = tmp_tid
                    yield DatasetFactory(**self.fetches_params)._acquire_module()
                    utils.remove_glob(tmp_tid)
        else:
            ## Standard GEBCO
            for gebco_fn in gebco_fns:
                self.fetches_params['mod'] = gebco_fn
                self.fetches_params['data_format'] = 200
                yield DatasetFactory(**self.fetches_params)._acquire_module()


class CopernicusFetcher(Fetcher):
    """Copernicus Fetcher."""
    
    __doc__ = f"{__doc__}\nFetches Module: <copernicus> - {fetches.copernicus.CopernicusDEM.__doc__}"
    
    def __init__(self, datatype=None, **kwargs):
        super().__init__(**kwargs)
        self.check_size = False
        self.datatype = datatype

        
    def yield_ds(self, result):
        from cudem.datalists.dlim import DatasetFactory
        
        if self.datatype is None or result['data_type'] == self.datatype:
            src_cop_dems = utils.p_unzip(
                os.path.join(self.fetch_module._outdir, result['dst_fn']),
                exts=['tif'], outdir=self.fetch_module._outdir, verbose=self.verbose
            )
            for src in src_cop_dems:
                gdalfun.gdal_set_ndv(src, ndv=0, verbose=False)
                self.fetches_params['mod'] = src
                self.fetches_params['data_format'] = 200
                self.fetches_params['node'] = 'pixel'
                yield DatasetFactory(**self.fetches_params)._acquire_module()


class FABDEMFetcher(Fetcher):
    """FABDEM Fetcher."""
    
    __doc__ = f"{__doc__}\nFetches Module: <fabdem> - {fetches.fabdem.FABDEM.__doc__}"
    
    def yield_ds(self, result):
        from cudem.datalists.dlim import DatasetFactory
        
        src_fab_dems = utils.p_unzip(
            os.path.join(self.fetch_module._outdir, result['dst_fn']),
            exts=['tif'], outdir=self.fetch_module._outdir, verbose=self.verbose
        )
        for src in src_fab_dems:
            gdalfun.gdal_set_ndv(src, ndv=0, verbose=False)
            self.fetches_params['mod'] = src
            self.fetches_params['data_format'] = 200
            yield DatasetFactory(**self.fetches_params)._acquire_module()


class MarGravFetcher(Fetcher):
    """Marine Gravity Fetcher."""
    
    __doc__ = f"{__doc__}\nFetches Module: <mar_grav> - {fetches.margrav.MarGrav.__doc__}"

    def __init__(self, rasterize=False, bathy_only=False, **kwargs):
        super().__init__(**kwargs)
        self.rasterize = rasterize
        self.bathy_only = bathy_only

    def yield_ds(self, result):
        from cudem.datalists.dlim import DatasetFactory
        
        if result['data_type'] == 'mar_grav_img':
            ## Convert IMG to GRD/NC
            nc_fn = utils.make_temp_fn(f"{utils.fn_basename2(result['dst_fn'])}.nc", temp_dir=self.fetch_module._outdir)
            
            cmd = f"gmt img2grd {os.path.join(self.fetch_module._outdir, result['dst_fn'])} {self.region.format('gmt')} -G{nc_fn} -D -T1 -I1m -E"
            utils.run_cmd(cmd, verbose=self.verbose)
            utils.run_cmd(f'gmt grdedit {nc_fn} -T')
            
            self.fetches_params['mod'] = nc_fn
            self.fetches_params['data_format'] = 200
            self.fetches_params['resample_and_warp'] = False
            self.fetches_params['node'] = 'grid'
            if self.bathy_only: self.fetches_params['upper_limit'] = 0
            
        elif self.rasterize:
            ## Waffles Rasterization
            from cudem import waffles
            mg_region = self.region.copy()
            if self.bathy_only: mg_region.zmax = 0
            
            mar_grav_fn = utils.make_temp_fn('mar_grav')
            _raster = waffles.WaffleFactory(
                mod='IDW:min_points=16',
                data=[f"{os.path.join(self.fetch_module._outdir, result['dst_fn'])},168:x_offset=REM,1"],
                src_region=mg_region,
                xinc=utils.str2inc('30s'), yinc=utils.str2inc('30s'),
                upper_limit=self.upper_limit,
                name=mar_grav_fn, node='pixel', verbose=self.verbose
            )._acquire_module()()
            
            self.fetches_params['mod'] = _raster.fn
            self.fetches_params['data_format'] = 200
            
        yield DatasetFactory(**self.fetches_params)._acquire_module()


class ChartsFetcher(Fetcher):
    """NOAA ENC Charts Fetcher."""
    
    __doc__ = f"{__doc__}\nFetches Module: <charts> - {fetches.charts.Charts.__doc__}"

    def __init__(self, want_soundings=True, want_contours=False, **kwargs):
        super().__init__(**kwargs)
        self.want_soundings = want_soundings
        self.want_contours = want_contours

    def yield_ds(self, result):
        from cudem.datalists.dlim import DatasetFactory
        
        src_000s = utils.p_unzip(
            os.path.join(self.fetch_module._outdir, result['dst_fn']),
            exts=['000'], outdir=self.fetch_module._outdir, verbose=self.verbose
        )
        
        for src_000 in src_000s:
            if self.want_soundings:
                self.fetches_params['mod'] = src_000
                self.fetches_params['data_format'] = 302
                self.fetches_params['ogr_layer'] = 'SOUNDG'
                self.fetches_params['z_scale'] = -1
                yield DatasetFactory(**self.fetches_params)._acquire_module()

            if self.want_contours:
                ## Basic Level Check
                try:
                    if int(os.path.basename(src_000)[2]) < 5: continue
                except: pass
                
                self.fetches_params['mod'] = src_000
                self.fetches_params['data_format'] = '302:ogr_layer=DEPCNT:elev_field=VALDCO:z_scale=-1'
                yield DatasetFactory(**self.fetches_params)._acquire_module()


class MBSFetcher(Fetcher):
    """NOAA Multibeam Fetcher."""
    
    __doc__ = f"{__doc__}\nFetches Module: <multibeam> - {fetches.multibeam.Multibeam.__doc__}"

    def __init__(self, mb_exclude='A', want_binned=False, want_mbgrid=False, auto_weight=True, **kwargs):
        super().__init__(**kwargs)
        self.fetches_params.update({
            'mb_exclude': mb_exclude,
            'want_binned': want_binned,
            'want_mbgrid': want_mbgrid,
            'auto_weight': auto_weight
        })

    def yield_ds(self, result):
        from cudem.datalists.dlim import DatasetFactory
        
        if not result['url'].endswith('.inf'):
            ## Ensure INFO file is parsed/fetched
            self.fetch_module.parse_entry_inf(result, keep_inf=True)
            yield DatasetFactory(**self.fetches_params)._acquire_module()


class HydroNOSFetcher(Fetcher):
    """NOAA HydroNOS Fetcher."""
    
    __doc__ = f"{__doc__}\nFetches Module: <hydronos> - {fetches.hydronos.HydroNOS.__doc__}"

    def __init__(self, explode=False, min_weight=0, **kwargs):
        super().__init__(**kwargs)
        self.explode = explode
        self.min_weight = min_weight

    def yield_ds(self, result):
        from cudem.datalists.dlim import DatasetFactory
        
        base_dir = os.path.dirname(os.path.join(self.fetch_module._outdir, result['dst_fn']))
        
        if result['data_type'] == 'xyz':
            nos_fns = utils.p_unzip(
                os.path.join(self.fetch_module._outdir, result['dst_fn']),
                exts=['xyz', 'dat'], outdir=base_dir, verbose=self.verbose
            )
            for nos_fn in nos_fns:
                self.fetches_params['mod'] = nos_fn
                self.fetches_params['data_format'] = '168:skip=1:xpos=2:ypos=1:zpos=3:z_scale=-1:delim=,'
                self.fetches_params['src_srs'] = 'epsg:4326+5866'
                yield DatasetFactory(**self.fetches_params)._acquire_module()

        elif result['data_type'] == 'bag':
            bag_fns = utils.p_unzip(
                os.path.join(self.fetch_module._outdir, result['dst_fn']),
                exts=['bag'], outdir=base_dir, verbose=self.verbose
            )
            for bag_fn in bag_fns:
                if 'ellipsoid' not in bag_fn.lower():
                    self.fetches_params['mod'] = bag_fn
                    self.fetches_params['data_format'] = 201
                    self.fetches_params['explode'] = self.explode
                    self.fetches_params['min_weight'] = self.min_weight
                    yield DatasetFactory(**self.fetches_params)._acquire_module()


class CSBFetcher(Fetcher):
    """CSB Fetcher."""
    
    __doc__ = f"{__doc__}\nFetches Module: <csb> - {fetches.csb.CSB.__doc__}"
    
    def yield_ds(self, result):
        from cudem.datalists.dlim import DatasetFactory
        yield DatasetFactory(**self.fetches_params)._acquire_module()


class R2RFetcher(Fetcher):
    """R2R Fetcher."""
    
    __doc__ = f"{__doc__}\nFetches Module: <r2r> - {fetches.multibeam.R2R.__doc__}"
    
    def yield_ds(self, result):
        from cudem.datalists.dlim import DatasetFactory
        
        r2r_fns = utils.p_untar(
            os.path.join(self.fetch_module._outdir, result['dst_fn']),
            exts=['geoCSV'],
            outdir=os.path.dirname(os.path.join(self.fetch_module._outdir, result['dst_fn'])),
            verbose=self.verbose
        )
        for fn in r2r_fns:
            self.fetches_params['mod'] = fn
            self.fetches_params['data_format'] = '168:skip=16:xpos=1:ypos=2:zpos=3:z_scale=-1:delim=,'
            self.fetches_params['src_srs'] = 'epsg:4326'
            yield DatasetFactory(**self.fetches_params)._acquire_module()


class EMODNetFetcher(Fetcher):
    """EMODNet Fetcher."""
    
    __doc__ = f"{__doc__}\nFetches Module: <emodnet> - {fetches.emodnet.EMODNet.__doc__}"
    
    def yield_ds(self, result):
        from cudem.datalists.dlim import DatasetFactory
        
        if result['data_type'] == 'csv':
            self.fetches_params['data_format'] = '168:skip=1:xpos=2:ypos=1:zpos=3:delim=,'
        elif result['data_type'] == 'nc':
            self.fetches_params['data_format'] = 200
            
        yield DatasetFactory(**self.fetches_params)._acquire_module()


class GEDTM30Fetcher(Fetcher):
    """GEDTM30 Data Fetcher"""
    
    __doc__ = '''{}\nFetches Module: <gedtm30> - {}'''.format(
        __doc__, fetches.gedtm30.GEDTM30.__doc__
    )

                                        
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

                                        
    def parse(self):
        for result in self.fetch_module.results:
            for this_ds in self.yield_ds(result):
                if this_ds is not None:
                    this_ds.remote = True
                    this_ds.initialize()
                    for ds in this_ds.parse():
                        yield ds

                                        
    def yield_ds(self, result):
        from cudem.datalists.dlim import DatasetFactory
        
        self.fetches_params['mod'] = '/vsicurl/{}'.format(result['url'])
        self.fetches_params['data_format'] = 200            
        yield DatasetFactory(**self.fetches_params)._acquire_module()


class HRDEMFetcher(Fetcher):
    """HRDEM Data Fetcher"""
    
    __doc__ = '''{}\nFetches Module: <hrdem> - {}'''.format(
        __doc__, fetches.hrdem.HRDEM.__doc__
    )

                                        
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        
    def yield_ds(self, result):
        from cudem.datalists.dlim import DatasetFactory
        
        yield DatasetFactory(**self.fetches_params)._acquire_module()

        
class eHydroFetcher(Fetcher):
    """USACE eHydro Fetcher."""
    
    __doc__ = f"{__doc__}\nFetches Module: <ehydro> - {fetches.ehydro.eHydro.__doc__}"

    def __init__(self, want_soundings=True, want_contours=True, **kwargs):
        super().__init__(**kwargs)
        self.want_soundings = want_soundings
        self.want_contours = want_contours

        
    def _get_projinfo_from_gdb(self, gdb_path):
        try:
            ds = ogr.Open(gdb_path)
            lyr = ds.GetLayer('SurveyPoint')
            src_srs = lyr.GetSpatialRef()
            datum = lyr[1].GetField('elevationDatum')
        except Exception as e:
            utils.echo_error_msg(f'Could not obtain srs info for {gdb_path}')
            #if not src_srs:
            src_srs = None
            #if not datum:
            datum = None
        finally:
            ds = lyr = None
            v_datum = vdatums.get_vdatum_by_name(datum)
            return src_srs, v_datum
        
        
    def yield_ds(self, result):
        from cudem.datalists.dlim import DatasetFactory
        
        try:
            src_gdb = utils.gdb_unzip(os.path.join(self.fetch_module._outdir, result['dst_fn']), 
                                      outdir=self.fetch_module._outdir, verbose=False)
        except Exception as e:
            utils.echo_error_msg('Could not extract {result}')
            src_gdb = None

        # if src_gdb.endswith('/'):
        #     src_gdb = src_gdb[:-1]

        if src_gdb:
            h, v = self._get_projinfo_from_gdb(src_gdb)

            #utils.debug_echo_msg(f'{src_gdb} projinfo: {h}, {v}')
            # gdb_ds = ogr.Open(src_gdb)
            # survey_point_layer = gdb_ds.GetLayer('SurveyPoint')
            # src_srs = survey_point_layer.GetSpatialRef()
            # elev_datum = survey_point_layer[1].GetField('elevationDatum')
            # v = vdatums.get_vdatum_by_name(elev_datum)
            # gdb_ds = survey_point_layer = None

            src_epsg = gdalfun.osr_parse_srs(h)
            utils.echo_debug_msg(f'{src_gdb} projinfo:\nhorz: {h}\nvert: {v}\nepsg: {src_epsg}')
            self.fetches_params['mod'] = src_gdb
            self.fetches_params['src_srs'] = f'{src_epsg}+{v or "5866"}' if src_epsg else None
            self.src_srs = f'{src_epsg}+{v or "5866"}' if src_epsg else None
            
            if self.want_soundings:
                self.metadata['name'] = src_gdb
                self.fetches_params['data_format'] = '302:ogr_layer=SurveyPoint_HD:elev_field=Z_label:z_scale=-0.3048006096012192'
                yield(DatasetFactory(**self.fetches_params)._acquire_module())            

            if self.want_contours:
                self.metadata['name'] = f'{utils.fn_basename2(src_gdb)}_contours'
                self.fetches_params['data_format']  = '302:ogr_layer=ElevationContour_ALL:elev_field=contourElevation:z_scale=-0.3048006096012192'
                yield(DatasetFactory(**self.fetches_params)._acquire_module())   
            
        # if src_gdb:
        #     ## Determine VDatum from GDB
        #     v = self._get_vdatum_from_gdb(src_gdb)
            
        #     ## Soundings
        #     ## Pass data to datalists.dlim.ogrfile
        #     if self.want_soundings:
        #         self.fetches_params['mod'] = src_gdb
        #         self.fetches_params['src_srs'] = f'epsg:4326+{v or "5866"}'
        #         self.fetches_params['data_format'] = '302:ogr_layer=SurveyPoint_HD:elev_field=Z_label:z_scale=-0.3048006096012192'
        #         self.metadata['name'] = self.fn
        #         yield DatasetFactory(**self.fetches_params)._acquire_module()            

        #     ## Contours
        #     ## Pass data to datalists.dlim.ogrfile
        #     if self.want_contours:
        #         self.metadata['name'] = f'{utils.fn_basename2(self.fn)}_contours'
        #         self.fetches_params['data_format'] = '302:ogr_layer=ElevationContour_ALL:elev_field=contourElevation:z_scale=-0.3048006096012192'
        #         yield DatasetFactory(**self.fetches_params)._acquire_module() 

        
    def yield_ds_XYZ(self, result):
        from cudem.datalists.dlim import DatasetFactory
        
        src_gdb = utils.gdb_unzip(
            os.path.join(self.fetch_module._outdir, result['dst_fn']),
            outdir=self.fetch_module._outdir,
            verbose=False
        )
        if src_gdb is not None:
            tmp_gdb = ogr.Open(src_gdb)
            tmp_layer = tmp_gdb.GetLayer('SurveyPoint')
            src_srs = tmp_layer.GetSpatialRef()
            src_epsg = gdalfun.osr_parse_srs(src_srs)
            tmp_gdb = None
            src_usaces = utils.p_unzip(
                os.path.join(self.fetch_module._outdir, result['dst_fn']),
                ['XYZ', 'xyz', 'dat'],
                outdir=self.fetch_module._outdir,
                verbose=self.verbose
            )
            for src_usace in src_usaces:
                self.fetches_params['mod'] = src_usace
                self.fetches_params['data_format'] = '168:z_scale=.3048'
                self.fetches_params['src_srs'] \
                    = '{}+{}'.format(src_epsg, v if v is not None else '5866') \
                    if src_epsg is not None \
                       else None
                yield(DatasetFactory(**self.fetches_params)._acquire_module())   
        

class BlueTopoFetcher(Fetcher):
    """BlueTopo Fetcher."""
    
    __doc__ = f"{__doc__}\nFetches Module: <bluetopo> - {fetches.bluetopo.BlueTopo.__doc__}"

    def __init__(self, want_interpolation=False, unc_weights=False, **kwargs):
        super().__init__(**kwargs)
        self.want_interpolation = want_interpolation
        self.unc_weights = unc_weights

    def yield_ds(self, result):
        from cudem.datalists.dlim import DatasetFactory
        
        fn = os.path.join(self.fetch_module._outdir, result['dst_fn'])
        sid = None
        
        if not self.want_interpolation:
            # Extract non-interpolated band (Band 3 usually?)
            sid = gdalfun.gdal_extract_band(
                fn, utils.make_temp_fn('tmp_bt_tid.tif', self.fetch_module._outdir),
                band=3, exclude=[0]
            )[0]

        if self.mask:
            new_mask = utils.make_temp_fn('test_tmp_mask')
            gdalfun.gdal_mask(
                sid, self.mask['mask'], new_mask, msk_value=1, verbose=True
            )
            os.replace(new_mask, sid)


        self.fetches_params['data_format'] = f'200:band_no=1:mask={sid}:uncertainty_mask=2' + (':weight_mask=2' if self.unc_weights else '')
        yield DatasetFactory(**self.fetches_params)._acquire_module()        


class NGSFetcher(Fetcher):
    """NGS Monument Fetcher."""
    
    __doc__ = f"{__doc__}\nFetches Module: <ngs> - {fetches.ngs.NGS.__doc__}"
    
    def __init__(self, datum='geoidHt', **kwargs):
        super().__init__(**kwargs)
        self.datum = datum if datum in ['orthoHt', 'geoidHt', 'z', 'ellipHeight'] else 'geoidHt'

    def yield_ds(self, result):
        from cudem.datalists.dlim import DatasetFactory
        
        # Parse JSON to XYZ
        xyz_fn = os.path.join(self.fetch_module._outdir, '_tmp_ngs.xyz')
        with open(os.path.join(self.fetch_module._outdir, result['dst_fn']), 'r') as f:
            data = json.load(f)
            with open(xyz_fn, 'w') as out:
                for row in data:
                    z = utils.float_or(row.get(self.datum))
                    if z is not None:
                        out.write(f"{row['lon']} {row['lat']} {z}\n")

        self.fetches_params['mod'] = xyz_fn
        self.fetches_params['data_format'] = 168
        yield DatasetFactory(**self.fetches_params)._acquire_module()
        utils.remove_glob(xyz_fn)


class TidesFetcher(Fetcher):
    """NOS Tides Fetcher."""
    
    __doc__ = f"{__doc__}\nFetches Module: <tides> - {fetches.tides.Tides.__doc__}"
    
    def __init__(self, s_datum='mllw', t_datum='msl', units='m', **kwargs):
        super().__init__(**kwargs)
        self.s_datum = s_datum
        self.t_datum = t_datum
        self.units = units

    def yield_ds(self, result):
        from cudem.datalists.dlim import DatasetFactory
        
        xyz_fn = os.path.join(self.fetch_module._outdir, '_tmp_tides.xyz')
        with open(os.path.join(self.fetch_module._outdir, result['dst_fn']), 'r') as f:
            data = json.load(f)
            with open(xyz_fn, 'w') as out:
                for feat in data.get('features', []):
                    attrs = feat['attributes']
                    v1 = attrs.get(self.s_datum, -99999.99)
                    v2 = attrs.get(self.t_datum, -99999.99)
                    
                    if v1 != -99999.99 and v2 != -99999.99:
                        z = v1 - v2
                        if self.units == 'm': z *= 0.3048
                        out.write(f"{attrs['longitude']} {attrs['latitude']} {z}\n")

        self.fetches_params['mod'] = xyz_fn
        self.fetches_params['data_format'] = 168
        yield DatasetFactory(**self.fetches_params)._acquire_module()
        utils.remove_glob(xyz_fn)


class WaterServicesFetcher(Fetcher):
    """USGS Water Services Fetcher."""
    __doc__ = f"{__doc__}\nFetches Module: <waterservices> - {fetches.waterservices.WaterServices.__doc__}"
    
    def __init__(self, site_code='00065', units='m', **kwargs):
        super().__init__(**kwargs)
        self.units = units
        self.site_code = site_code

    def yield_ds(self, result):
        from cudem.datalists.dlim import DatasetFactory
        
        xyz_fn = os.path.join(self.fetch_module._outdir, '_tmp_ws.xyz')
        with open(os.path.join(self.fetch_module._outdir, result['dst_fn']), 'r') as f:
            data = json.load(f)
            with open(xyz_fn, 'w') as out:
                for feat in data.get('value', {}).get('timeSeries', []):
                    if feat['variable']['variableCode'][0]['value'] == self.site_code:
                        loc = feat['sourceInfo']['geoLocation']['geogLocation']
                        val = float(feat['values'][0]['value'][0]['value'])
                        if self.units == 'm': val *= 0.3048
                        out.write(f"{loc['longitude']} {loc['latitude']} {val}\n")

        self.fetches_params['mod'] = xyz_fn
        self.fetches_params['data_format'] = 168
        yield DatasetFactory(**self.fetches_params)._acquire_module()
        utils.remove_glob(xyz_fn)

        
class VDatumFetcher(Fetcher):
    """VDatum transformation grids.
    """

    __doc__ = '''{}\nFetches Module: <vdatum> - {}'''.format(__doc__, fetches.vdatum.VDATUM.__doc__)

    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        
    def yield_ds(self, result):
        from cudem.datalists.dlim import DatasetFactory
        
        src_tif = os.path.join(
            self.fetch_module._outdir, '{}.tif'.format(
                utils.fn_basename2(os.path.basename(result['dst_fn']))
            )
        )
        if result['dst_fn'].endswith('.zip'):
            v_gtx = utils.p_f_unzip(
                os.path.join(
                    self.fetch_module._outdir, result['dst_fn']
                ), [result['data_type']], outdir=self.fetch_module._outdir, verbose=self.verbose
            )[0]
            status = utils.run_cmd(
                f'gdalwarp "{v_gtx}" "{src_tif}" -t_srs epsg:4269 --config CENTER_LONG 0',
                verbose=self.verbose
            )

        self.fetches_params['mod'] = src_tif
        self.fetches_params['data_format'] = 200
        self.fetches_params['node'] = 'pixel'
        yield DatasetFactory(**self.fetches_params)._acquire_module()

        
## todo: allow lakes bathymetry
## as well as lakes breaklines (shape nodes)
## see: https://www.esri.com/arcgis-blog/products/arcgis-pro/3d-gis/hydro-flattening-of-river-shorelines-in-lidar-based-dem-production/
class HydroLakesFetcher(Fetcher):
    """HydroLakes lake bathymetric data
    """
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        
    def yield_ds(self, result):
        pass    
        
### End
