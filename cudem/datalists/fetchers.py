### fetchers.py 
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
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
###############################################################################
### Commentary:
##
### Examples:
##
### TODO:
##
### Code:

import os
import json
import copy
#from tqdm import tqdm
from osgeo import ogr

from cudem import utils
from cudem import regions
from cudem import gdalfun
from cudem import vdatums
from cudem import fetches
from cudem.datalists.dlim import ElevationDataset

class Fetcher(ElevationDataset):
    """The default fetches dataset type; dlim Fetcher dataset class

    This is used in waffles/dlim for on-the-fly remote data
    parsing and processing.
    
    See `fetches`, `cudem.fetches` for more information on usage and the
    various fetches modules supported.

    If a fetch module needs special processing define a sub-class
    of Fetcher and redefine the yield_ds(self, result) function which 
    yields a list of dlim dataset objects, where result is an item from 
    the fetch result list. Otherwise, this Fetcher class can be used as 
    default if the fetched data comes in a normal sort of way.

    Generally, though not always, if the fetched data is a raster then
    there is no need to redefine yield_ds, though if the raster has 
    insufficient information, such as with Copernicus, whose nodata value 
    is not specified in the geotiff files, it may be best to create a 
    simple sub-class for it.
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
                 **kwargs
    ):
        super().__init__(**kwargs)
        ## cache directory to store fetched data
        #self.outdir = outdir if outdir is not None else self.cache_dir        
        self.wgs_region = self.region.copy()
        self.wgs_region.src_srs = self.dst_srs
        self.wgs_srs = 'epsg:4326'
        if self.dst_srs is not None:
            self.wgs_region.warp(self.wgs_srs)

        if outdir is None:
            outdir = self.cache_dir
            
        self.fetch_module = fetches.fetches.FetchesFactory(
            mod=self.fn, src_region=self.wgs_region,
            callback=callback, verbose=True,
            outdir=outdir,
        )._acquire_module() # the fetches module
        if self.fetch_module is None:
            utils.echo_warning_msg(
                f'fetches modules {self.fn} returned None'
            )
            pass

        self.mask_coast = mask_coast
        self.invert_coast = invert_coast
        self.coast_buffer = utils.float_or(coast_buffer, 0.00001)
        self.want_single_metadata_name = want_single_metadata_name
        # check the size of the fetched data
        self.check_size = check_size
        # retain fetched data after processing
        self.keep_fetched_data = keep_fetched_data 
        self.outdir = self.fetch_module._outdir
        if not os.path.exists(os.path.dirname(self.outdir)):
            os.makedirs(os.path.dirname(self.outdir))

        self.outdir = os.path.abspath(self.outdir)        
        ##self.cache_dir = self.outdir                
        # try:
        #     self.fetch_module.run()
        # except Exception as e:
        #     utils.echo_warning_msg(
        #         f'fetch module {self.fn} returned zero results, {e}'
        #     )
        #     self.fetch_module.results = []

        ## breaks when things not set...
        # src_horz, src_vert = gdalfun.epsg_from_input(self.fetch_module.src_srs)

        self._reset_params()


    def _sub_init(self):
        if len(self.fetch_module.results) == 0:
            try:
                self.fetch_module.run()
            except Exception as e:
                utils.echo_warning_msg(
                    f'fetch module {self.fn} returned zero results, {e}'
                )
                self.fetch_module.results = []
        
    def init_fetch_module(self):
        self.fetch_module = fetches.fetches.FetchesFactory(
            mod=self.fn,
            src_region=self.region,
            callback=fetches.fetches.fetches_callback,
            verbose=self.verbose,
            outdir=outdir
        )._acquire_module() # the fetches module
        if self.fetch_module is None:
            utils.echo_warning_msg(
                f'fetches modules {self.fn} returned None'
            )
            pass

        try:
            self.fetch_module.run()
        except Exception as e:
            utils.echo_warning_msg(
                f'fetch module {self.fn} returned zero results, {e}'
            )
            self.fetch_module.results = []

            
    def _reset_params(self):
        ## set the metadata from the fetches module
        md = copy.deepcopy(self.metadata)
        md['name'] = self.metadata['name']
        md['title'] = self.fetch_module.title
        md['source'] = self.fetch_module.source
        md['date'] = self.fetch_module.date
        md['data_type'] = self.data_format
        md['resolution'] = self.fetch_module.resolution
        md['hdatum'] = self.fetch_module.hdatum
        md['vdatum'] = self.fetch_module.vdatum
        md['url'] = self.fetch_module.url
        
        #self.fetches_params = self._set_params(**self.fetches_params)
        self.fetches_params = self._set_params(
            data_format=self.fetch_module.data_format,
            src_srs=self.fetch_module.src_srs,
            cache_dir=self.fetch_module._outdir,
            remote=True,
            metadata=md,
            #parent=self,
        )

        if self.mask is None and self.mask_coast:
            from cudem.fetches import osm
            this_osm = osm.osmCoastline(
                region=self.region, chunks=True, verbose=self.verbose, attempts=5, cache_dir=self.cache_dir
            )
            coast_mask = this_osm(return_geom=False, overwrite=False)

            this_osm_water = osm.osmCoastline(
                region=self.region, chunks=True, verbose=self.verbose, attempts=5, cache_dir=self.cache_dir, q='water'
            )
            water_mask = this_osm_water(return_geom=False, overwrite=False)
            if coast_mask is not None:
                self.fetches_params['mask'] = [f'mask_fn={coast_mask}:invert={self.invert_coast}']

            if water_mask is not None:
                if self.fetches_params['mask'] is not None:                    
                    self.fetches_params['mask'].append(f'mask_fn={water_mask}:invert={self.invert_coast}')
                else:
                    self.fetches_params['mask'] = [f'mask_fn={water_mask}:invert={self.invert_coast}']
        
        
    def generate_inf(self):
        """generate a infos dictionary from the Fetches dataset"""

        tmp_region = self.fetch_module.region \
            if self.region is None \
               else self.region.copy()
        self.infos.minmax = tmp_region.export_as_list()    
        self.infos.wkt = tmp_region.export_as_wkt()
        self.infos.src_srs = self.fetch_module.src_srs
        return(self.infos)

    
    def parse(self):
        #self.init_fetch_module()
        with utils.ccp(
                total=len(self.fetch_module.results),
                desc='parsing datasets from datalist fetches {} @ {}'.format(
                    self.fetch_module, self.weight
                ),
            leave=False
        ) as pbar:
            for result in self.fetch_module.results:
                ## mrl commented out to set params in sub-modules
                #self._reset_params()
                status = self.fetch_module.fetch_entry(
                    result, check_size=self.check_size
                )
                if status == 0:
                    self.fetches_params['mod'] = os.path.join(
                        self.fetch_module._outdir, result['dst_fn']
                    )

                    #self.fetches_params['mod'] = result['dst_fn']
                    #utils.echo_msg(self.fetches_params)
                    for this_ds in self.yield_ds(result):
                        if this_ds is not None:
                            #this_ds.initialize()
                            f_name = os.path.relpath(
                                this_ds.fn.split(':')[0], self.fetch_module._outdir
                            )
                            if f_name == '.':
                                f_name = this_ds.fn.split(':')[0]

                            mod_name = os.path.dirname(utils.fn_basename2(f_name))
                            if mod_name == '':
                                mod_name = self.fetch_module.name

                            #if self.want_single_metadata_name:
                            #    this_ds.metadata['name'] = mod_name
                            #else:
                            this_ds.metadata['name'] \
                                = '/'.join(
                                    ['/'.join(
                                        this_ds.metadata['name'].split('/')[:-1]
                                    ),
                                     f_name]
                                )
                            #this_ds.remote = True
                            this_ds.initialize()
                            for ds in this_ds.parse():
                                #ds.initialize()
                                self.data_entries.append(ds)
                                yield(ds)
                        else:
                            utils.echo_warning_msg(
                                f'could not set fetches datasource {result}'
                            )
                else:
                    utils.echo_warning_msg(
                        f'data not fetched {status}:{result}'
                    )
                        
                pbar.update()
                
        if not self.keep_fetched_data:
            utils.remove_glob(f'{self.fn}*')

            
    def yield_ds(self, result):
        from cudem.datalists.dlim import DatasetFactory
        
        ## try to get the SRS info from the result if it's a gdal file
        ## fix this.
        try:
            vdatum = self.fetch_module.vdatum
            src_srs = gdalfun.gdal_get_srs(
                os.path.join(
                    self.fetch_module._outdir, result['dst_fn']
                )
            )
            horz_epsg, vert_epsg = gdalfun.epsg_from_input(src_srs)
            if vert_epsg is None:
                vert_epsg = vdatum

            if vert_epsg is not None:
                self.fetch_module.src_srs = f'{horz_epsg}+{vert_epsg}'
            else:
                self.fetch_module.src_srs = src_srs
                
        except:
            pass

        yield(DatasetFactory(**self.fetches_params)._acquire_module())

        
class NEDFetcher(Fetcher):
    """National Elevation Dataset from USGS

    This is a wrapper shortcut for fetching NED DEMs from USGS' 
    The National Map (TNM)
    """

    __doc__ = '''{}\nFetches Module: <ned> - {}'''.format(
        __doc__, fetches.tnm.NED.__doc__
    )

    
    def __init__(self, coast_buffer=0.00001, remove_flat=True, **kwargs):
        super().__init__(**kwargs)
        self.remove_flat = remove_flat
        self.coast_buffer = utils.float_or(coast_buffer, 0.00001)

        
    def yield_ds(self, result):
        from cudem.datalists.dlim import DatasetFactory
        
        ## todo: merge the coast mask with user input self.mask
        coast_mask = None
        ned_mask = self.mask
        src_dem = os.path.join(self.fetch_module._outdir, result['dst_fn'])
            
        self.fetches_params['mod'] = src_dem
        self.fetches_params['mask'] = ned_mask
        self.fetches_params['remove_flat'] = self.remove_flat
        yield(DatasetFactory(**self.fetches_params)._acquire_module())

        if coast_mask is not None:
            utils.remove_glob(
                '{}.*'.format(utils.fn_basename2(coast_mask))
            )

            
class DNRFetcher(Fetcher):
    """
    """

    __doc__ = '''{}\nFetches Module: <wadnr> - {}'''.format(
        __doc__, fetches.wadnr.waDNR.__doc__
    )
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        self.data_format = -2

        
    def yield_ds(self, result):
        from cudem.datalists.dlim import DatasetFactory
        
        src_dnr_dems = utils.p_unzip(
            os.path.join(self.fetch_module._outdir, result['dst_fn']),
            exts=['tif'],
            outdir=self.fetch_module._outdir,
            verbose=self.verbose
        )
        for src_dnr_dem in src_dnr_dems:
            self.fetches_params['mod'] = src_dnr_dem
            #self.fetches_params['data_format'] = 200
            #self.fetches_params['node'] = 'pixel'
            self.fetches_params['remove_flat'] = True
            yield(DatasetFactory(**self.fetches_params)._acquire_module())

            
class DAVFetcher_CoNED(Fetcher):
    """CoNED from the digital coast 

    This is a wrapper shortcut for fetching CoNED DEMs from the Digital Coast,
    mainly so we can pull the vertical datum info from the DAV metadata 
    since the CoNED doesn't assign one to their DEMs.
    """

    __doc__ = '''{}\nFetches Module: <CoNED> - {}'''.format(
        __doc__, fetches.dav.CoNED.__doc__
    )
    
    
    def __init__(self, keep_fetched_data=True, cog=True, **kwargs):
        super().__init__(**kwargs)
        self.keep_fetched_data = keep_fetched_data
        self.cog = cog

        
    def parse(self):
        #self.fetch_module.run()
        for result in self.fetch_module.results:
            if not self.cog:
                status = self.fetch_module.fetch_entry(
                    result, check_size=self.check_size
                )
                if status != 0:
                    break
                
            for this_ds in self.yield_ds(result):
                if this_ds is not None:
                    this_ds.initialize()
                    yield(this_ds)

                    
    def yield_ds(self, result):
        from cudem.datalists.dlim import DatasetFactory
        
        ## try to get the SRS info from the result
        try:
            vdatum = self.fetch_module.vdatum
            if not self.cog:
                src_srs = gdalfun.gdal_get_srs(
                    os.path.join(self.fetch_module._outdir, result['dst_fn'])
                )
            else:
                src_srs = gdalfun.gdal_get_srs(result['url'])
                
            horz_epsg, vert_epsg = gdalfun.epsg_from_input(src_srs)
            if vert_epsg is None:
                vert_epsg = vdatum

            #utils.echo_msg('srs: {}+{}'.format(horz_epsg, vert_epsg))
            if vert_epsg is not None:
                self.fetch_module.src_srs = f'{horz_epsg}+{vert_epsg}'
            else:
                self.fetch_module.src_srs = src_srs
                
        except:
            pass

        self.fetches_params['mod'] = os.path.join(
            self.fetch_module._outdir, result['dst_fn']
        ) if not self.cog else result['url']
        
        self.fetches_params['check_path'] = False if self.cog else True
        self.fetches_params['src_srs'] = self.fetch_module.src_srs
        self.fetches_params['data_format'] = 200
        
        yield(DatasetFactory(**self.fetches_params)._acquire_module())

        
class DAVFetcher_SLR(Fetcher):
    """SLR DEM from the digital coast 

    This is a wrapper shortcut for fetching SLR DEMs from the Digital Coast,
    mainly so we can pull and remove the flattened ring around the actual data.
    """

    __doc__ = '''{}Fetches Module: <SLR> - {}'''.format(
        __doc__, fetches.dav.SLR.__doc__
    )

    
    def __init__(self, keep_fetched_data = True, **kwargs):
        super().__init__(**kwargs)

        
    def yield_ds(self, result):
        from cudem.datalists.dlim import DatasetFactory
        
        ## this doesn't work in all cases,
        ## update to find and remove flattened areas
        # gdalfun.gdal_set_ndv(
        #     os.path.join(
        #         self.fetch_module._outdir, result[1]
        #     ), ndv=-99.0000, convert_array=True
        # )
        self.fetches_params['remove_flat'] = True
        yield(DatasetFactory(**self.fetches_params)._acquire_module())    

        
class SWOTFetcher(Fetcher):
    """SWOT L2_HR_Raster data from NASA

    set `data_set` to one of (for L2_HR_Raster):
    1 - longitude (64-bit floating-point)
    2 - latitude (64-bit floating-point)
    3 - wse (32-bit floating-point)
    4 - status_flag (8-bit unsigned integer)
    5 - status_flag (32-bit unsigned integer)
    6 - wse_uncert (32-bit floating-point)
    7 - water_area (32-bit floating-point)
    8 - status_flag (8-bit unsigned integer)
    9 - status_flag (32-bit unsigned integer)
    10 - water_area_uncert (32-bit floating-point)
    11 - water_frac (32-bit floating-point)
    12 - water_frac_uncert (32-bit floating-point)
    13 - sig0 (32-bit floating-point)
    14 - status_flag (8-bit unsigned integer)
    15 - status_flag (32-bit unsigned integer)
    16 - sig0_uncert (32-bit floating-point)
    17 - inc (32-bit floating-point)
    18 - cross_track (32-bit floating-point)
    19 - time (64-bit floating-point)
    20 - time (64-bit floating-point)
    21 - n_wse_pix (32-bit unsigned integer)
    22 - n_water_area_pix (32-bit unsigned integer)
    23 - n_sig0_pix (32-bit unsigned integer)
    24 - n_other_pix (32-bit unsigned integer)
    25 - dark_frac (32-bit floating-point)
    26 - status_flag (8-bit unsigned integer)
    27 - status_flag (8-bit unsigned integer)
    28 - layover_impact (32-bit floating-point)
    29 - sig0_cor_atmos_model (32-bit floating-point)
    30 - height_cor_xover (32-bit floating-point)
    31 - geoid_height_above_reference_ellipsoid (32-bit floating-point)
    32 - solid_earth_tide (32-bit floating-point)
    33 - load_tide_fes (32-bit floating-point)
    34 - load_tide_got (32-bit floating-point)
    35 - pole_tide (32-bit floating-point)
    36 - model_dry_tropo_cor (32-bit floating-point)
    37 - model_wet_tropo_cor (32-bit floating-point)
    38 - iono_cor_gim_ka (32-bit floating-point)

    currently supported SWOT products and product options:
    - L2_HR_Raster
      apply_geoid
      data_set

    - L2_HR_PIXC
      apply_geoid
      
    """

    __doc__ = '''{}\nFetches Module: <swot> - {}'''.format(
        __doc__, fetches.earthdata.SWOT.__doc__
    )

    
    def __init__(self,
                 data_set='wse',
                 apply_geoid=True,
                 classes=None,
                 classes_qual=None,
                 anc_classes=None,
                 remove_class_flags=False,
                 **kwargs):
        super().__init__(**kwargs)
        self.fetches_params['var'] = data_set
        self.fetches_params['apply_geoid'] = apply_geoid
        self.fetches_params['classes'] = classes
        self.fetches_params['anc_classes'] = anc_classes
        self.fetches_params['classes_qual'] = classes_qual
        self.fetches_params['remove_class_flags'] = remove_class_flags

        
    # def fetch_pixc_vec(self, swot_fn):
    #     pixc_vec_filter = utils.fn_basename2(swot_fn).split('PIXC_')[1]
    #     this_pixc_vec = fetches.SWOT(
    #         src_region=None,
    #         verbose=self.verbose,
    #         outdir=self.fetch_module._outdir,
    #         product='L2_HR_PIXCVec',
    #         filename_filter=pixc_vec_filter
    #     )
    #     this_pixc_vec.run()
    #     if len(this_pixc_vec.results) == 0:
    #         utils.echo_warning_msg(
    #             'could not locate associated PIXCVec file for {}'.format(pixc_vec_filter)
    #         )
    #         return(None)
    #     else:
    #         if this_pixc_vec.fetch(this_pixc_vec.results[0], check_size=self.check_size) == 0:
    #             return(os.path.join(this_pixc_vec._outdir, this_pixc_vec.results[0][1]))

    
    def yield_ds(self, result):
        from cudem.datalists.dlim import DatasetFactory
        
        #swot_fn = os.path.join(self.fetch_module._outdir, result[1])
        if 'L2_HR_PIXC_' in result['data_type']:
            #pixc_vec_result = self.fetch_pixc_vec(swot_fn)
            #swot_pixc_vec_fn = pixc_vec_result
            self.fetches_params['data_format'] = 202
            yield(DatasetFactory(**self.fetches_params)._acquire_module())
        elif 'L2_HR_Raster' in result['data_type']:
            self.fetches_params['data_format'] = 203
            yield(DatasetFactory(**self.fetches_params)._acquire_module())
        else:
            utils.echo_warning_msg(
                (f'{result["data_type"]} is not a currently supported '
                 'dlim dataset')
            )

            
class IceSat2Fetcher(Fetcher):
    """IceSat2 data from NASA
    
    """

    from cudem.datalists.icesat2file import IceSat2_ATL03
    from cudem.fetches import osm
            
    __doc__ = '''{}\nDLIM Module: <IceSat2File> - {}'''.format(
        __doc__, IceSat2_ATL03.__doc__
    )
    
    __doc__ = '''{}\nFetches Module: <icesat2> - {}'''.format(
        __doc__, fetches.earthdata.IceSat2.__doc__
    )

    
    def __init__(self,
                 water_surface=None,
                 classes=None,
                 confidence_levels=None,
                 columns={},
                 classify_bathymetry=False,
                 classify_buildings=True,
                 classify_water=True,
                 reject_failed_qa=True,
                 min_bathy_confidence=None,
                 **kwargs):
        super().__init__(**kwargs)
        self.water_surface = water_surface
        self.classes = classes
        self.confidence_levels = confidence_levels
        self.columns = columns
        self.classify_bathymetry = classify_bathymetry
        self.classify_buildings = classify_buildings
        self.classify_water = classify_water
        self.reject_failed_qa = reject_failed_qa
        self.min_bathy_confidence = min_bathy_confidence        
        self.data_format = -111
        self.fetches_params['classify_bathymetry'] = self.classify_bathymetry
        self.fetches_params['classify_buildings'] = self.classify_buildings
        self.fetches_params['classify_water'] = self.classify_water

        if self.fetches_params['classify_buildings']:
            this_bfp = bingbfp.bingBuildings(
                region=self.region, verbose=self.verbose, cache_dir=self.cache_dir
            )
            bldg_mask = this_bfp(return_geom=True)
            self.fetches_params['classify_buildings'] = bldg_mask

        if self.fetches_params['classify_water']:
            this_osm = osm.osmCoastline(
                region=self.region, chunks=True, verbose=self.verbose, attempts=5, cache_dir=self.cache_dir
            )
            coast_mask = this_osm(return_geom=True)

            this_osm_water = osm.osmCoastline(
                region=self.region, chunks=True, verbose=self.verbose, attempts=5, cache_dir=self.cache_dir, q='water'
            )
            coast_mask.extend(this_osm_water(return_geom=True))
            self.fetches_params['classify_water'] = coast_mask

        
    def yield_ds(self, result):
        from cudem.datalists.dlim import DatasetFactory
        
        icesat2_fn= os.path.join(
            self.fetch_module._outdir, result['dst_fn']
        )

        if result['data_type'].lower() == 'atl03' or \
           'NSIDC_CPRD' in result['data_type'].upper():
            self.fetches_params['water_suface'] = self.water_surface if self.water_surface is not None else 'geoid'
            self.fetches_params['classes'] = self.classes
            self.fetches_params['confidence_levels'] = self.confidence_levels
            self.fetches_params['columns'] = self.columns
            self.fetches_params['reject_failed_qa'] = self.reject_failed_qa
            self.fetches_params['min_bathy_confidence'] = self.min_bathy_confidence

            if 'processed_zip' in result['data_type']:
                icesat2_h5s = utils.p_unzip(
                    icesat2_fn,
                    exts=['h5'],
                    outdir=self.cache_dir,
                    verbose=self.verbose
                )
                for icesat2_h5 in icesat2_h5s:
                    self.fetches_params['fn'] = icesat2_h5
                    yield(datalists.IceSat2File(**self.fetches_params)._acquire_module())
                    yield(sub_ds)

            else:
                self.fetches_params['fn'] = icesat2_fn
                self.fetches_params['data_format'] = 303
                #yield(IceSat2File(**self.fetches_params)._acquire_module())
                yield(DatasetFactory(**self.fetches_params)._acquire_module())
                
        elif result['data_type'].lower() == 'atl24':
            self.fetches_params['water_suface'] = self.water_surface if self.water_surface is not None else 'ortho'
            self.fetches_params['classes'] = self.classes
            self.fetches_params['min_confidence'] = self.min_bathy_confidence
            self.fetches_params['fn'] = icesat2_fn
            self.fetches_params['data_format'] = 304
            yield(DatasetFactory(**self.fetches_params)._acquire_module())
        else:
            utils.echo_warning_msg(f'{icesat2_fn}({result["dst_fn"]} - {result["data_type"]} cannot be processed')

            
class GMRTFetcher(Fetcher):
    """GMRT Gridded data.

    -----------
    Parameters:
    
    swath_only: only return MB swath data
    """
    
    __doc__ = '''{}\nFetches Module: <gmrt> - {}'''.format(
        __doc__, fetches.gmrt.GMRT.__doc__
    )

    
    def __init__(self, swath_only=False, **kwargs):
        super().__init__(**kwargs)
        self.swath_only = swath_only

        
    def yield_ds(self, result):
        from cudem.datalists.dlim import DatasetFactory
        
        swath_mask=None
        gmrt_fn = os.path.join(self.fetch_module._outdir, result['dst_fn'])
        with gdalfun.gdal_datasource(gmrt_fn, update = 1) as src_ds:
            md = src_ds.GetMetadata()
            md['AREA_OR_POINT'] = 'Point'
            src_ds.SetMetadata(md)
            gdalfun.gdal_set_srs(src_ds)
            gdalfun.gdal_set_ndv(src_ds, verbose=False)

        if self.swath_only:
            if fetches.fetches.Fetch(
                    self.fetch_module._gmrt_swath_poly_url,
                    verbose=self.verbose
            ).fetch_file(
                os.path.join(
                    self.fetch_module._outdir, 'gmrt_swath_polygons.zip'
                )
            ) == 0:
                swath_shps = utils.p_unzip(
                    os.path.join(
                        self.fetch_module._outdir, 'gmrt_swath_polygons.zip'
                    ),
                    exts=['shp', 'shx', 'prj', 'dbf'],
                    outdir=self.cache_dir,
                    verbose=self.verbose
                )
                for v in swath_shps:
                    if '.shp' in v:
                        #swath_shp = v
                        ## upate to new masking
                        swath_mask = {'mask': v, 'invert_mask': True}
                        break
                    
                if not os.path.exists(swath_mask['mask']):
                    utils.echo_error_msg(
                        'could not find gmrt swath polygons...'
                    )
                    self.swath_only = False
                    swath_mask = None
                else:
                    self.fetches_params['mask'] = swath_mask

        yield(DatasetFactory(**self.fetches_params)._acquire_module())

        
class GEBCOFetcher(Fetcher):
    """GEBCO Gridded data

    Note: Fetches entire zip file.
    """
    
    __doc__ = '''{}\nFetches Module: <gebco> - {}'''.format(
        __doc__, fetches.gebco.GEBCO.__doc__
    )

    
    def __init__(self, exclude_tid=None, **kwargs):
        super().__init__(**kwargs)
        self.exclude_tid = []
        if utils.str_or(exclude_tid) is not None:
            for tid_key in exclude_tid.split('/'):
                self.exclude_tid.append(utils.int_or(tid_key))

                
    def yield_ds(self, result):
        from cudem.datalists.dlim import DatasetFactory
        
        wanted_gebco_fns = []
        gebco_fns = utils.p_unzip(
            os.path.join(self.fetch_module._outdir, result['dst_fn']),
            exts=['tif'],
            outdir=self.fetch_module._outdir,
            verbose=self.verbose
        )
        
        ## fetch the TID zip if needed
        if self.exclude_tid:
            if fetches.fetches.Fetch(
                    self.fetch_module._gebco_urls['gebco_tid']['geotiff'],
                    verbose=self.verbose
            ).fetch_file(
                os.path.join(self.fetch_module._outdir, 'gebco_tid.zip')
            ) == 0:
                tid_fns = utils.p_unzip(
                    os.path.join(self.fetch_module._outdir, 'gebco_tid.zip'),
                    exts=['tif'],
                    outdir=self.cache_dir,
                    verbose=self.verbose
                )
                for tid_fn in tid_fns:
                    ds_config = gdalfun.gdal_infos(tid_fn)                    
                    if self.region is not None \
                       and self.region.valid_p(check_xy=True):
                        inf_region = regions.Region().from_geo_transform(
                            ds_config['geoT'], ds_config['nx'], ds_config['ny']
                        )
                        inf_region.wmin = self.weight
                        inf_region.wmax = self.weight
                        inf_region.umin = self.uncertainty
                        inf_region.umax = self.uncertainty
                        if regions.regions_intersect_p(inf_region, self.region):
                            wanted_gebco_fns.append(tid_fn)
                    else:
                        wanted_gebco_fns.append(tid_fn)

                for tid_fn in wanted_gebco_fns:
                    tmp_tid = utils.make_temp_fn(
                        'tmp_tid.tif', temp_dir=self.fetch_module._outdir
                    )
                    with gdalfun.gdal_datasource(tid_fn) as tid_ds:
                        tid_config = gdalfun.gdal_infos(tid_ds)
                        tid_band = tid_ds.GetRasterBand(1)
                        tid_array = tid_band.ReadAsArray().astype(float)

                    tid_config['ndv'] = -9999
                    tid_config['dt'] = gdal.GDT_Float32                        
                    for tid_key in self.exclude_tid:
                        tid_array[tid_array == tid_key] = tid_config['ndv']
                        
                    for tid_key in self.fetch_module.tid_dic.keys():
                        tid_array[tid_array == tid_key] \
                            = self.fetch_module.tid_dic[tid_key][1]
                            
                    gdalfun.gdal_write(tid_array, tmp_tid, tid_config)
                    if self.mask is not None:
                        new_mask = utils.make_temp_fn('test_tmp_mask')
                        gdalfun.gdal_mask(
                            tmp_tid, self.mask['mask'], new_mask,
                            msk_value=1, verbose=True
                        )
                        os.replace(new_mask, tmp_tid)

                    self.fetches_params['mod'] = tid_fn.replace('tid_', '')
                    self.fetches_params['data_format'] = 200
                    self.fetches_params['mask'] = tmp_tid
                    self.fetches_params['weight_mask'] = tmp_tid
                    yield(DatasetFactory(**self.fetches_params)._acquire_module())
                    utils.remove_glob(tmp_tid)
        else:
            for gebco_fn in gebco_fns:
                ds_config = gdalfun.gdal_infos(gebco_fn)
                inf_region = regions.Region().from_geo_transform(
                    ds_config['geoT'], ds_config['nx'], ds_config['ny']
                )
                if self.region is not None \
                   and self.region.valid_p(check_xy=True):                
                    inf_region.wmin = self.weight
                    inf_region.wmax = self.weight
                    inf_region.umin = self.uncertainty
                    inf_region.umax = self.uncertainty

                    if regions.regions_intersect_p(inf_region, self.region):
                        wanted_gebco_fns.append(gebco_fn)
                else:
                    wanted_gebco_fns.append(gebco_fn)    
            
            for gebco_fn in wanted_gebco_fns:
                self.fetches_params['mod'] = gebco_fn
                self.fetches_params['data_format'] = 200
                yield(DatasetFactory(
                    **self.fetches_params
                )._acquire_module())

                
class CopernicusFetcher(Fetcher):
    """Gridded Copernicus sattelite data.
    """
    
    __doc__ = '''{}\nFetches Module: <copernicus> - {}'''.format(
        __doc__, fetches.copernicus.CopernicusDEM.__doc__
    )

    
    def __init__(self, datatype=None, **kwargs):
        super().__init__(**kwargs)
        self.check_size=False
        self.datatype=datatype

        
    def yield_ds(self, result):
        from cudem.datalists.dlim import DatasetFactory
        
        if self.datatype is None \
           or result['data_type'] == self.datatype:
            src_cop_dems = utils.p_unzip(
                os.path.join(
                    self.fetch_module._outdir, result['dst_fn']
                ),
                exts=['tif'],
                outdir=self.fetch_module._outdir,
                verbose=self.verbose
            )
            for src_cop_dem in src_cop_dems:
                gdalfun.gdal_set_ndv(src_cop_dem, ndv=0, verbose=False)
                self.fetches_params['mod'] = src_cop_dem
                self.fetches_params['data_format'] = 200
                self.fetches_params['node'] = 'pixel'
                yield(DatasetFactory(
                    **self.fetches_params
                )._acquire_module())

                
class FABDEMFetcher(Fetcher):
    """FABDEM Gridded data
    """
    
    __doc__ = '''{}
    Fetches Module: <fabdem> - {}'''.format(
        __doc__, fetches.fabdem.FABDEM.__doc__
    )

    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        
    def yield_ds(self, result):
        from cudem.datalists.dlim import DatasetFactory
        
        src_fab_dems = utils.p_unzip(
            os.path.join(self.fetch_module._outdir, result['dst_fn']),
            exts=['tif'],
            outdir=self.fetch_module._outdir,
            verbose=self.verbose
        )
        for src_fab_dem in src_fab_dems:
            gdalfun.gdal_set_ndv(src_fab_dem, ndv=0, verbose=False)
            self.fetches_params['mod'] = src_fab_dem
            self.fetches_params['data_format'] = 200
            yield(DatasetFactory(**self.fetches_params)._acquire_module())

            
class MarGravFetcher(Fetcher):
    """Marine Gravity Bathymetry
    """
    
    __doc__ = '''{}\nFetches Module: <mar_grav> - {}'''.format(
        __doc__, fetches.margrav.MarGrav.__doc__
    )

    
    def __init__(self,
                 rasterize=False,
                 bathy_only=False,
                 #upper_limit=None,
                 #lower_limit=None,
                 **kwargs):
        super().__init__(**kwargs)
        self.rasterize = rasterize
        self.bathy_only = bathy_only
        #self.upper_limit = utils.float_or(upper_limit)
        #self.lower_limit = utils.float_or(lower_limit)
        #if self.bathy_only:
        #    self.upper_limit = 0

        #self.region.zmax=self.upper_limit
        #self.region.zmin=self.lower_limit

        
    def yield_ds(self, result):
        from cudem.datalists.dlim import DatasetFactory
        
        if result['data_type'] == 'mar_grav_img':
            nc_fn = utils.make_temp_fn(
                '{}.nc'.format(utils.fn_basename2(result['dst_fn'])),
                temp_dir=self.fetch_module._outdir)
            img2grd_cmd = 'gmt img2grd {} {} -G{} -D -T1 -I1m -E'.format(
                os.path.join(
                    self.fetch_module._outdir,
                    result['dst_fn']
                ), self.region.format('gmt'), nc_fn
            )
            out, status = utils.run_cmd(img2grd_cmd, verbose=True)
            out, status = utils.run_cmd(f'gmt grdedit {nc_fn} -T')
            self.fetches_params['mod'] = nc_fn
            self.fetches_params['data_format'] = 200
            self.fetches_params['resample_and_warp'] = False
            self.fetches_params['node'] = 'grid'
            if self.bathy_only:
                self.fetches_params['upper_limit'] = 0
            
        elif self.rasterize:
            from cudem import waffles
            mg_region = self.region.copy()
            if self.bathy_only:
                mg_region.zmax = 0
                #mg_region.zmin = self.lower_limit        
            mar_grav_fn = utils.make_temp_fn('mar_grav')
            _raster = waffles.WaffleFactory(
                mod='IDW:min_points=16',
                data=['{},168:x_offset=REM,1'.format(
                    os.path.join(self.fetch_module._outdir, result['dst_fn'])
                )],
                src_region=mg_region,
                xinc=utils.str2inc('30s'),
                yinc=utils.str2inc('30s'),
                upper_limit = self.upper_limit,
                name=mar_grav_fn,
                node='pixel',
                verbose=True
            )._acquire_module()()            
            #if self.upper_limit is not None or self.lower_limit is not None:
            ds = gdal.Open(_raster.fn)
            ds_band = ds.GetRasterBand(1)
            ds_arr = ds_band.ReadAsArray()
            # if self.upper_limit is not None:
            #     ds_arr[ds_arr >= self.upper_limit] = ds_band.GetNoDataValue()

            # if self.lower_limit is not None:
            #     ds_arr[ds_arr <= self.lower_limit] = ds_band.GetNoDataValue()
                    
            ds = None

            self.fetches_params['mod'] = _raster.fn
            self.fetches_params['data_format'] = 200
            
        yield(DatasetFactory(**self.fetches_params)._acquire_module())

                                        
class ChartsFetcher(Fetcher):
    """NOAA ENC Charts Fetcher

    Digital Soundings
    """
    
    __doc__ = '''{}\nFetches Module: <charts> - {}'''.format(
        __doc__, fetches.charts.Charts.__doc__
    )

                                        
    def __init__(self, want_soundings=True, want_contours=False, **kwargs):
        super().__init__(**kwargs)
        self.want_soundings = want_soundings
        self.want_contours = want_contours

                                        
    def yield_ds(self, result):
        from cudem.datalists.dlim import DatasetFactory
        
        src_000s = utils.p_unzip(
            os.path.join(self.fetch_module._outdir, result['dst_fn']),
            exts=['000'],
            outdir=self.fetch_module._outdir,
            verbose=self.verbose            
        )
        for src_000 in src_000s:
            if self.want_soundings:
                self.fetches_params['mod'] = src_000
                self.fetches_params['data_format'] = 302
                self.fetches_params['ogr_layer'] = 'SOUNDG'
                self.fetches_params['z_scale'] = -1
                yield(DatasetFactory(**self.fetches_params)._acquire_module())

            if self.want_contours:
                enc_level = utils.int_or(
                    os.path.basename(src_000)[2], 0
                )
                if enc_level < 5:
                    continue
                
                self.fetches_params['mod'] = src_000
                self.metadata['name'] = '{}_contours'.format(
                    utils.fn_basename2(self.fn)
                )
                self.fetches_params['DEPCNT'] = 'SOUNDG'
                self.fetches_params['data_format'] \
                    = '302:ogr_layer=DEPCNT:elev_field=VALDCO:z_scale=-1'
                yield(DatasetFactory(**self.fetches_params)._acquire_module())            

                                        
class MBSFetcher(Fetcher):
    """NOAA Multibeam Fetcher
    """

    __doc__ = '''{}\nFetches Module: <multibeam> - {}'''.format(
        __doc__, fetches.multibeam.Multibeam.__doc__
    )

                                        
    def __init__(self, mb_exclude='A', want_binned=False,
                 want_mbgrid=False, auto_weight=True, **kwargs):
        super().__init__(**kwargs)
        self.fetches_params['mb_exclude'] = mb_exclude
        self.fetches_params['want_binned'] = want_binned
        self.fetches_params['want_mbgrid'] = want_mbgrid
        self.fetches_params['auto_weight'] = auto_weight

        
    def yield_ds(self, result):
        from cudem.datalists.dlim import DatasetFactory
        
        if not result['url'].endswith('.inf'):
            mb_infos = self.fetch_module.parse_entry_inf(
                result, keep_inf=True
            )
            yield(DatasetFactory(**self.fetches_params)._acquire_module())

                                        
class HydroNOSParser(ElevationDataset):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

                                        
class HydroNOSFetcher(Fetcher):
    """NOAA HydroNOS Data Fetcher
    """
    
    __doc__ = '''{}\nFetches Module: <hydronos> - {}'''.format(
        __doc__, fetches.hydronos.HydroNOS.__doc__
    )

                                        
    def __init__(self, explode=False, min_weight=0, **kwargs):
        super().__init__(**kwargs)
        self.explode=explode
        self.min_weight = min_weight

                                        
    def yield_ds(self, result):
        from cudem.datalists.dlim import DatasetFactory
        
        if result['data_type'] == 'xyz':
            nos_fns = utils.p_unzip(
                os.path.join(self.fetch_module._outdir, result['dst_fn']),
                exts=['xyz', 'dat'],
                outdir=os.path.dirname(
                    os.path.join(self.fetch_module._outdir, result['dst_fn'])
                ),
                verbose=self.verbose
            )
            for nos_fn in nos_fns:
                self.fetches_params['mod'] = nos_fn
                self.fetches_params['data_format'] \
                    = ('168:skip=1:xpos=2:ypos=1'
                       ':zpos=3:z_scale=-1:delim=,')
                self.fetches_params['src_srs'] = 'epsg:4326+5866'
                yield(DatasetFactory(**self.fetches_params)._acquire_module())

        elif result['data_type'] == 'bag':
            bag_fns = utils.p_unzip(
                os.path.join(self.fetch_module._outdir, result['dst_fn']),
                exts=['bag'],
                outdir=os.path.dirname(
                    os.path.join(self.fetch_module._outdir, result['dst_fn'])
                ),
                verbose=self.verbose
            )
            for bag_fn in bag_fns:
                # bag_fn = os.path.join(self.fetch_module._outdir, result[1])
                if 'ellipsoid' not in bag_fn.lower():
                    src_srs = gdalfun.gdal_get_srs(
                        os.path.join(
                            self.fetch_module._outdir, result['dst_fn']
                        )
                    )
                    #utils.echo_msg(src_srs)
                    #self.src_srs = src_srs
                    self.fetches_params['mod'] = bag_fn
                    self.fetches_params['data_format'] = 201
                    #self.fetches_params['src_srs'] = src_srs
                    self.fetches_params['explode'] = self.explode
                    self.fetches_params['min_weight'] = self.min_weight
                    yield(DatasetFactory(
                        **self.fetches_params
                    )._acquire_module())

                                        
class CSBFetcher(Fetcher):
    """Crowd Sourced Bathymetry data fetcher
    """
    
    __doc__ = '''{}\nFetches Module: <csb> - {}'''.format(
        __doc__, fetches.csb.CSB.__doc__
    )

                                        
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

                                        
    def yield_ds(self, result):
        from cudem.datalists.dlim import DatasetFactory
        yield(DatasetFactory(**self.fetches_params)._acquire_module())


class R2RFetcher(Fetcher):
    """R2R Bathymetry data fetcher
    """
    
    __doc__ = '''{}\nFetches Module: <r2r> - {}'''.format(__doc__, fetches.multibeam.R2R.__doc__)

                                        
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

                                        
    def yield_ds(self, result):
        #utils.echo_msg_bold(result)
        from cudem.datalists.dlim import DatasetFactory
        
        r2r_fns = utils.p_untar(
            os.path.join(self.fetch_module._outdir, result['dst_fn']),
            exts=['geoCSV'],
            outdir=os.path.dirname(
                os.path.join(self.fetch_module._outdir, result['dst_fn'])
            ),
            verbose=self.verbose
        )
        #utils.echo_msg_bold(r2r_fns)
        for r2r_fn in r2r_fns:
            self.fetches_params['mod'] = r2r_fn
            self.fetches_params['data_format'] \
                = ('168:skip=16:xpos=1:ypos=2'
                   ':zpos=3:z_scale=-1:delim=,')
            self.fetches_params['src_srs'] = 'epsg:4326'
            yield(DatasetFactory(**self.fetches_params)._acquire_module())

                                        
class EMODNetFetcher(Fetcher):
    """EMODNet Data Fetcher
    """
    
    __doc__ = '''{}
    Fetches Module: <emodnet> - {}'''.format(
        __doc__, fetches.emodnet.EMODNet.__doc__
    )

                                        
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

                                        
    def yield_ds(self, result):
        from cudem.datalists.dlim import DatasetFactory
        
        if result['data_type'] == 'csv':
            self.fetches_params['data_format'] \
                = '168:skip=1:xpos=2:ypos=1:zpos=3:delim=,'
        elif result['data_type'] == 'nc':
            self.fetches_params['data_format'] = 200
            
        yield(DatasetFactory(**self.fetches_params)._acquire_module())

                                        
class GEDTM30Fetcher(Fetcher):
    """GEDTM30 Data Fetcher
    """
    
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
                        yield(ds)

                                        
    def yield_ds(self, result):
        from cudem.datalists.dlim import DatasetFactory
        
        self.fetches_params['mod'] = '/vsicurl/{}'.format(result['url'])
        self.fetches_params['data_format'] = 200            
        yield(DatasetFactory(**self.fetches_params)._acquire_module())


class HRDEMFetcher(Fetcher):
    """HRDEM Data Fetcher
    """
    
    __doc__ = '''{}\nFetches Module: <hrdem> - {}'''.format(
        __doc__, fetches.hrdem.HRDEM.__doc__
    )

                                        
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        
    def yield_ds(self, result):
        from cudem.datalists.dlim import DatasetFactory
        
        yield(DatasetFactory(**self.fetches_params)._acquire_module())

                                        
class eHydroFetcher(Fetcher):
    """USACE eHydro soundings
    """

    __doc__ = '''{}\nFetches Module: <ehydro> - {}'''.format(__doc__, fetches.ehydro.eHydro.__doc__)

                                        
    def __init__(self, want_contours=True, **kwargs):
        super().__init__(**kwargs)
        self.want_contours = want_contours

                                        
    def yield_ds(self, result):
        from cudem.datalists.dlim import DatasetFactory
        
        try:
            src_gdb = utils.gdb_unzip(
                os.path.join(self.fetch_module._outdir, result['dst_fn']),
                outdir=self.fetch_module._outdir,
                verbose=False
            )
        except Exception as e:
            utils.echo_error_msg(f'{self.fn}: {e}')
            src_gdb = None

        if src_gdb is not None:
            tmp_gdb = ogr.Open(src_gdb)
            tmp_layer = tmp_gdb.GetLayer('SurveyPoint')
            src_srs = tmp_layer.GetSpatialRef()
            src_epsg = gdalfun.osr_parse_srs(src_srs)
            elev_datum = tmp_layer[1].GetField('elevationDatum')
            tmp_gdb = tmp_layer = None
            v = vdatums.get_vdatum_by_name(elev_datum)
            self.fetches_params['mod'] = src_gdb
            self.fetches_params['src_srs'] = f'{src_epsg}+{v if v is not None else "5866"}' \
                if src_epsg is not None else None
            self.fetches_params['data_format'] \
                = ('302:ogr_layer=SurveyPoint_HD'
                   ':elev_field=Z_label'
                   ':z_scale=-0.3048006096012192')
            self.metadata['name'] = self.fn
            yield(DatasetFactory(**self.fetches_params)._acquire_module())            

            if self.want_contours:
                self.metadata['name'] = f'{utils.fn_basename2(self.fn)}_contours'
                self.fetches_params['data_format'] \
                    = ('302:ogr_layer=ElevationContour_ALL'
                       ':elev_field=contourElevation'
                       ':z_scale=-0.3048006096012192')
                yield(DatasetFactory(**self.fetches_params)._acquire_module())            

                
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
    """BlueTopo Gridded bathymetric data Fetcher

    -----------
    Parameters:
    
    want_interpolation: True/False to include interpolated cells
    unc_weights: use the uncertainty mask as weights
    """

    __doc__ = '''{}\nFetches Module: <bluetopo> - {}'''.format(
        __doc__, fetches.bluetopo.BlueTopo.__doc__
    )

    
    def __init__(self, want_interpolation=False, unc_weights=False,
                 **kwargs):
        super().__init__(**kwargs)
        self.want_interpolation = want_interpolation
        self.unc_weights = unc_weights

        
    def yield_ds(self, result):
        from cudem.datalists.dlim import DatasetFactory
        
        sid = None
        if not self.want_interpolation:
            sid = gdalfun.gdal_extract_band(
                os.path.join(self.fetch_module._outdir, result['dst_fn']),
                utils.make_temp_fn(
                    'tmp_bt_tid.tif', self.fetch_module._outdir
                ),
                band=3,
                exclude=[0]
            )[0]

        if self.mask is not None:
            new_mask = utils.make_temp_fn('test_tmp_mask')
            gdalfun.gdal_mask(
                sid, self.mask['mask'], new_mask, msk_value=1, verbose=True
            )
            os.replace(new_mask, sid)

        self.fetches_params['data_format'] \
            = '200:band_no=1:mask={}:uncertainty_mask=2{}'.format(
                sid, ':weight_mask=2' if self.unc_weights else ''
            )
        
        yield(DatasetFactory(**self.fetches_params)._acquire_module())        

        
class NGSFetcher(Fetcher):
    """NGS Monument data
    """

    __doc__ = '''{}\nFetches Module: <ngs> - {}'''.format(__doc__, fetches.ngs.NGS.__doc__)

    
    def __init__(self, datum = 'geoidHt', **kwargs):
        super().__init__(**kwargs)
        self.datum = datum
        if self.datum not in ['orthoHt', 'geoidHt', 'z', 'ellipHeight']:
            utils.echo_warning_msg(
                f'could not parse {datum}, falling back to geoidHt'
                )
            self.datum = 'geoidHt'

            
    def yield_ds(self, result):
        from cudem.datalists.dlim import DatasetFactory
        
        with open(
                os.path.join(
                    self.fetch_module._outdir, result['dst_fn']
                ), 'r'
        ) as json_file:
            r = json.load(json_file)
            if len(r) > 0:
                with open(
                        os.path.join(
                            self.fetch_module._outdir, '_tmp_ngs.xyz'
                        ), 'w'
                ) as tmp_ngs:
                    for row in r:
                        z = utils.float_or(row[self.datum])
                        if z is not None:
                            xyz = xyzfun.XYZPoint(src_srs='epsg:4326').from_list(
                                [float(row['lon']), float(row['lat']), z]
                            )
                            xyz.dump(dst_port=tmp_ngs)

        self.fetches_params['mod'] = os.path.join(
            self.fetch_module._outdir, '_tmp_ngs.xyz'
        )
        self.fetches_params['data_format'] = 168
        yield(DatasetFactory(**self.fetches_params)._acquire_module())

        utils.remove_glob(
            os.path.join(
                self.fetch_module._outdir, '_tmp_ngs.xyz'
            )
        )

        
class TidesFetcher(Fetcher):
    """NOS Tide Station data
    """

    __doc__ = '''{}\nFetches Module: <tides> - {}'''.format(__doc__, fetches.tides.Tides.__doc__)

    
    def __init__(self, s_datum='mllw', t_datum='msl', units='m', **kwargs):
        super().__init__(**kwargs)
        self.s_datum = s_datum
        self.t_datum = t_datum
        self.units = units

        
    def yield_ds(self, result):
        from cudem.datalists.dlim import DatasetFactory
        
        with open(
                os.path.join(self.fetch_module._outdir, result['dst_fn']), 'r'
        ) as json_file:
            r = json.load(json_file)
            if len(r) > 0:
                with open(
                        os.path.join(
                            self.fetch_module._outdir, '_tmp_tides.xyz'
                        ), 'w'
                ) as tmp_ngs:
                    for feature in r['features']:
                        # if self.fetch_module.station_id is not None:
                        #     if self.fetch_module.station_id != feature['attributes']['id']:
                        #         continue
                            
                        lon = feature['attributes']['longitude']
                        lat = feature['attributes']['latitude']
                        if feature['attributes'][self.s_datum] != -99999.99 \
                           and feature['attributes'][self.t_datum] != -99999.99:
                            z = feature['attributes'][self.s_datum] \
                                - feature['attributes'][self.t_datum]
                            if self.units == 'm':
                                z = z * 0.3048

                            xyz = xyzfun.XYZPoint(
                                src_srs='epsg:4326'
                            ).from_list([lon, lat, z])
                            xyz.dump(dst_port=tmp_ngs)

        self.fetches_params['mod'] \
            = os.path.join(self.fetch_module._outdir, '_tmp_tides.xyz')
        self.fetches_params['data_format'] = 168
        yield(DatasetFactory(**self.fetches_params)._acquire_module())

        utils.remove_glob(
            os.path.join(self.fetch_module._outdir, '_tmp_tides.xyz')
        )

        
class WaterServicesFetcher(Fetcher):
    """USGS Water Services

    -----------
    Parameters:
    
    site_code: the site code to fetch
    units: 'm' for meters

    site_codes:
    00065 - Gate Height
    00060 - StreamFlow
    63160 - Stream water level elevation above NAVD 1988
    62611 - Groundwater level above NAVD 1988
    72019 - Depth to water level, units below land surface
    """

    __doc__ = '''{}\nFetches Module: <waterservices> - {}'''.format(
        __doc__, fetches.waterservices.WaterServices.__doc__
    )

    
    def __init__(self, site_code='00065', units='m', **kwargs):
        super().__init__(**kwargs)
        self.units = units
        self.site_code = site_code

        
    def yield_ds(self, result):
        from cudem.datalists.dlim import DatasetFactory
        
        with open(
                os.path.join(self.fetch_module._outdir, result['dst_fn']), 'r'
        ) as json_file:
            r = json.load(json_file)
            if len(r) > 0:
                with open(
                        os.path.join(self.fetch_module._outdir, '_tmp_ws.xyz'), 'w'
                ) as tmp_ws:
                    features = r['value']['timeSeries']
                    for feature in features:
                        if feature['variable']['variableCode'][0]['value'] \
                           == self.site_code:
                            lon = float(
                                feature['sourceInfo']['geoLocation']['geogLocation']['longitude']
                            )
                            lat = float(
                                feature['sourceInfo']['geoLocation']['geogLocation']['latitude']
                            )
                            z = float(
                                feature['values'][0]['value'][0]['value']
                            )

                            if self.units == 'm':
                                z = z * 0.3048
                            
                            xyz = xyzfun.XYZPoint(
                                src_srs='epsg:4326'
                            ).from_list([lon, lat, z])
                            xyz.dump(dst_port=tmp_ws)

        self.fetches_params['mod'] = os.path.join(
            self.fetch_module._outdir, '_tmp_ws.xyz'
        )
        self.fetches_params['data_format'] = 168
        yield(DatasetFactory(**self.fetches_params)._acquire_module())

        utils.remove_glob(
            os.path.join(self.fetch_module._outdir, '_tmp_ws.xyz')
        )

        
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
        yield(DatasetFactory(**self.fetches_params)._acquire_module())

        
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
