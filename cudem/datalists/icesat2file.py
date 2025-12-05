### icesat2file.py - DataLists IMproved
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
##
## icesat2file.py is part of CUDEM
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
import glob
import numpy as np
import h5py as h5
import pandas as pd
from tqdm import tqdm

from osgeo import gdal
from osgeo import ogr

from cudem import utils
from cudem import regions
from cudem import gdalfun
from cudem.fetches import fetches
from cudem.fetches import earthdata
from cudem.datalists.dlim import ElevationDataset

class IceSat2_ATL24(ElevationDataset):
    def __init__(self, min_confidence=None, classes='40', water_surface='ortho', **kwargs):
        super().__init__(**kwargs)
        self.orientDict = {0:'l', 1:'r', 21:'error'}
        self.lasers = ['gt1l', 'gt2l', 'gt3l', 'gt1r', 'gt2r', 'gt3r']

        self.data_format = 304
        self.water_surface = water_surface
        if self.water_surface not in ['surface', 'ortho', 'ellipse']:
            self.water_surface = 'ortho'
        
        self.water_surface
        self.min_confidence = utils.float_or(min_confidence)
        self.classes = [int(x) for x in classes.split('/')] \
            if classes is not None \
               else []        
        
    def generate_inf(self):
        if self.src_srs is None:
            self.infos.src_srs = 'epsg:4326+3855'
        else:
            self.infos.src_srs = self.src_srs

        f = h5.File(self.fn, 'r')
        these_regions = []
        numpts = 0
        for b in range(1, 4):
            for p in ['l', 'r']:
                if f'gt{b}{p}' not in list(f.keys()):
                    continue
                
                y = f[f'gt{b}{p}/lat_ph'][...,]
                x = f[f'gt{b}{p}/lon_ph'][...,]
                z = f[f'gt{b}{p}/ortho_h'][...,]
                this_region = regions.Region(src_srs=self.infos.src_srs).from_list(
                    [np.min(x), np.max(x), np.min(y), np.max(y), np.min(z), np.max(z)]
                )
                these_regions.append(this_region)
                numpts += len(x)

        region_cnt = len(these_regions)
        if region_cnt > 0:
            merged_region = these_regions[0].copy()
            for r in range(1, region_cnt):
                merged_region = regions.regions_merge(merged_region, these_regions[r])

            self.infos.minmax = merged_region.export_as_list(include_z=True)
            self.infos.wkt = merged_region.export_as_wkt()
            
        self.infos.numpts = numpts 
        
        f.close()
        return(self.infos)
    
        
    def yield_points(self):
        dataset = None

        f = h5.File(self.fn)
        for b in range(1, 4):
            for p in ['l', 'r']:
        
                if f'gt{b}{p}' not in list(f.keys()):
                    continue

                lat_ph = f[f'gt{b}{p}/lat_ph'][...,]
                lon_ph = f[f'gt{b}{p}/lon_ph'][...,]
                ortho_h = f[f'gt{b}{p}/ortho_h'][...,]
                surface_h = f[f'gt{b}{p}/surface_h'][...,]
                ellipse_h = f[f'gt{b}{p}/ellipse_h'][...,]
                class_ph = f[f'gt{b}{p}/class_ph'][...,]
                index_seg = f[f'gt{b}{p}/index_seg'][...,]
                index_ph = f[f'gt{b}{p}/index_ph'][...,]
                conf_ph = f[f'gt{b}{p}/confidence'][...,]
                low_confidence_ph = f[f'gt{b}{p}/low_confidence_flag'][...,]

                ## append the laser to each record
                laser_arr = np.empty(ortho_h.shape, dtype='object')
                laser_arr[:] = f'gt{b}{p}'

                ## append the filename to each record
                fn_arr = np.empty(ortho_h.shape, dtype='object')
                fn_arr[:] = self.fn

                if self.water_surface == 'surface':
                    ph_height = surface_h
                elif self.water_surface == 'ortho':
                    ph_height = ortho_h
                else:
                    ph_height = ellipse_h
                
                dataset = pd.DataFrame(
                    {'latitude': lat_ph,
                     'longitude': lon_ph,
                     'photon_height': ph_height,
                     'laser': laser_arr,
                     'fn': fn_arr,
                     'confidence': conf_ph,
                     'low_confidence_ph': low_confidence_ph,
                     'ph_h_classed': class_ph},
                    columns=['latitude', 'longitude', 'photon_height', 'laser', 'fn',
                             'confidence', 'low_confidence_ph', 'ph_h_classed']
                )

                ## keep only photons with a classification mentioned in `self.classes`
                if len(self.classes) > 0:
                    dataset = dataset[
                        (np.isin(dataset['ph_h_classed'], self.classes))
                    ]

                if self.min_confidence is not None:
                    dataset = dataset[dataset['confidence'] >= self.min_confidence]
                    
                if dataset is None or len(dataset) == 0:
                    continue
                
                ## rename the x,y,z columns for `transform_and_yield_points`
                dataset.rename(
                    columns={
                        'longitude': 'x', 'latitude': 'y', 'photon_height': 'z'
                    },
                    inplace=True
                )
                
                yield(dataset)                
                
        f.close()

        
class IceSat2_ATL03_old(ElevationDataset):
    """IceSat2 data from NASA

    Parameters:
    
    water_surface: 'geoid' # this is the vertical datum, can be 'geoid', 
    'ellipsoid' or 'mean_tide'

    classes: None # return only data with the specified classes, e.g. '2/3/4'
    # return only data with the specified confidence levels, e.g. '2/3/4'
    confidence_levels: None 
    columns: {} # the additional columns to export in transform_and_yield_points 
             {'/h5/atl/path', column_name}
    classify_bathymetry: True # extract bathymetry with CShelph
    classify_buildings: True # classify buildings using BING BFP
    classify_water: True # classify water using OSM
    reject_failed_qa: True # skip granules that failed QA
    min_bathy_confidence: <0-1> the minimum passable ATL24 confidence level

    Classes:
    -1 - no classification (ATL08)
    0 - noise / atmosphere (ATL08)
    1 - ground surface (ATL08)
    2 - canopy (ATL08)
    3 - canopy top (ATL08)
    40 - bathymetry floor surface (CShelph, ATL24)
    41 - bathymetry water surface (OSM coastline, ATL24)
    6 - ice surface (ATL06) 
        (unused for now, just planning ahead for possible future 
        ATL06 integration)
    7 - built structure (OSM or Bing)
    8 - "urban" (WSF, if used)
    9 - inland water surface

    Confidence Levels:
    0, 1, 2, 3, 4
    """
    
    def __init__(self,
                 water_surface='geoid',
                 classes=None,
                 confidence_levels='2/3/4',
                 columns={},
                 classify_bathymetry=False,
                 classify_buildings=True,
                 classify_inland_water=True,
                 reject_failed_qa=True,
                 classify_water=True,
                 append_atl24=False,
                 min_bathy_confidence=None,
                 **kwargs):
        super().__init__(**kwargs)
        self.data_format = 303
        self.water_surface = water_surface
        if self.water_surface not in ['mean_tide', 'geoid', 'ellipsoid']:
            self.water_surface = 'mean_tide'

        self.classes = [int(x) for x in classes.split('/')] \
            if classes is not None \
               else []
        self.confidence_levels = [int(x) for x in confidence_levels.split('/')] \
            if confidence_levels is not None \
               else []
        self.columns = columns
        #self.columns = {'/gtx/heights/delta_time': 'delta_time_b'}
        self.atl_fn = None
        self.want_bathymetry = classify_bathymetry
        self.want_buildings = classify_buildings
        self.want_watermask = classify_water
        self.want_inland_water = classify_inland_water
        self.reject_failed_qa = reject_failed_qa
        self.append_atl24 = append_atl24
        self.min_bathy_confidence = utils.float_or(min_bathy_confidence)

        self.orientDict = {0:'l', 1:'r', 21:'error'}
        self.lasers = ['gt1l', 'gt2l', 'gt3l', 'gt1r', 'gt2r', 'gt3r']

        
    def generate_inf(self):
        if self.src_srs is None:
            self.infos.src_srs = 'epsg:4326+3855'
        else:
            self.infos.src_srs = self.src_srs

        f = h5.File(self.fn, 'r')
        these_regions = []
        numpts = 0
        for b in range(1, 4):
            for p in ['l', 'r']:
                if f'gt{b}{p}' not in list(f.keys()):
                    continue

                if 'heights' not in list(f[f'gt{b}{p}'].keys()):
                    continue
                
                y = f[f'gt{b}{p}/heights/lat_ph'][...,]
                x = f[f'gt{b}{p}/heights/lon_ph'][...,]
                z = f[f'gt{b}{p}/heights/h_ph'][...,]
                this_region = regions.Region(src_srs=self.infos.src_srs).from_list(
                    [np.min(x), np.max(x), np.min(y), np.max(y), np.min(z), np.max(z)]
                )
                these_regions.append(this_region)
                numpts += len(x)

        region_cnt = len(these_regions)
        if region_cnt > 0:
            merged_region = these_regions[0].copy()
            for r in range(1, region_cnt):
                merged_region = regions.regions_merge(merged_region, these_regions[r])

            self.infos.minmax = merged_region.export_as_list(include_z=True)
            self.infos.wkt = merged_region.export_as_wkt()
            
        self.infos.numpts = numpts 
        
        f.close()
        return(self.infos)

        
    def init_atl_h5(self):
        """initialize all the relevant ATL h5 files"""

        self.atl03_f = None
        self.atl08_f = None
        self.atl24_f = None
        #self.atl12_f = None
        #self.atl13_f = None

        if self.atl03_fn is not None and os.path.exists(self.atl03_fn):
            self.atl03_f = h5.File(self.atl03_fn, 'r')            
            if 'short_name' not in self.atl03_f.attrs.keys():
                raise UnboundLocalError(
                    f'{self.atl03_fn} does not appear to be an ATL03 file'
                )
        
            if self.reject_failed_qa:
                if self.atl03_f[
                        '/quality_assessment/qa_granule_pass_fail'
                ][...,][0] != 0:
                    raise UnboundLocalError(
                        'this granule has failed qa {}'.format(
                            self.atl03_f[
                                '/quality_assessment/qa_granule_pass_fail'
                            ][...,][0]
                        )
                    )

            if not np.any(np.isin(list(self.atl03_f.keys()), self.lasers)):
                raise UnboundLocalError(
                    f'this granule ({self.atl03_fn}) has no lasers'
                )
            
        if self.atl08_fn is not None and os.path.exists(self.atl08_fn):
            #utils.echo_msg_bold(self.atl08_fn)
            self.atl08_f = h5.File(self.atl08_fn, 'r')
            if self.verbose:
                utils.echo_msg(f'Using ATL08 file: {self.atl08_fn}')
                
            if 'short_name' not in self.atl08_f.attrs.keys():
                utils.echo_warning_msg(
                    f'{self.atl08_fn} does not appear to be an ATL file'
                )
                self.atl08_f.close()

        ## atl24 doesn't have the short_name attr
        if self.atl24_fn is not None and os.path.exists(self.atl24_fn):
            self.atl24_f = h5.File(self.atl24_fn, 'r')
            if self.verbose:
                utils.echo_msg(f'Using ATL24 file: {self.atl24_fn}')
                
            # if 'short_name' not in self.atl24_f.attrs.keys():
            #     utils.echo_warning_msg(
            #         f'{self.atl24_fn} does not appear to be an ATL file'
            #     )
            #     self.atl24_f.close()

        else:
            #utils.echo_warning_msg('falling back to CSphelph to classify bathymetry')
            utils.echo_warning_msg('falling back to classify water points only')
            #self.want_bathymetry = True

            #if self.want_watermask:
            if isinstance(self.want_watermask, bool):
                self.want_watermask = True
                this_wm = self.process_coastline(
                    self.fetch_coastline(chunks=False, verbose=self.verbose),
                    return_geom=True,
                    verbose=self.verbose
                )
            elif isinstance(self.want_watermask, list):
                this_wm = self.want_watermask
                

        # if self.atl12_fn is not None and os.path.exists(self.atl12_fn):
        #     self.atl12_f = h5.File(self.atl12_fn, 'r')
        #     if self.verbose:
        #         utils.echo_msg(f'Using ATL12 file: {self.atl12_fn}')
                
        #     if 'short_name' not in self.atl13_f.attrs.keys():
        #         utils.echo_warning_msg(
        #             f'{self.atl12_fn} does not appear to be an ATL file'
        #         )
        #         self.atl12_f.close()                

        
        # if self.atl13_fn is not None and os.path.exists(self.atl13_fn):
        #     self.atl13_f = h5.File(self.atl13_fn, 'r')
        #     if self.verbose:
        #         utils.echo_msg(f'Using ATL13 file: {self.atl13_fn}')
                
        #     if 'short_name' not in self.atl13_f.attrs.keys():
        #         utils.echo_warning_msg(
        #             f'{self.atl13_fn} does not appear to be an ATL file'
        #         )
        #         self.atl13_f.close()                

                
    def close_atl_h5(self):
        """close all open atl files"""
        
        self.atl03_fn = self.atl08_fn = self.atl24_fn = None
        if self.atl03_f is not None:
            self.atl03_f.close()
            
        if self.atl08_f is not None:
            self.atl08_f.close()

        if self.atl24_f is not None:
            self.atl24_f.close()

        # if self.atl12_f is not None:
        #     self.atl12_f.close()
            
        # if self.atl13_f is not None:
        #     self.atl13_f.close()

            
    def yield_points(self):
        """yield the points from the dataset.

        In this case, it will yield a pandas dataframe
        """

        dataset = None
        self.atl03_fn = self.fn
        self.atl08_fn = None
        self.atl24_fn = None
        # self.atl12_fn = None
        # self.atl13_fn = None

        ## only fetch atl08/atl24 if classes are desired
        if len(self.classes) > 0:
            atl08_result = self.fetch_atlxx(short_name='ATL08')
            self.atl08_fn = atl08_result
            atl24_result = self.fetch_atlxx(short_name='ATL24')
            self.atl24_fn = atl24_result

            # atl12_result = self.fetch_atlxx(short_name='ATL12')
            # self.atl12_fn = atl12_result
            # atl13_result = self.fetch_atlxx(short_name='ATL13')
            # self.atl13_fn = atl13_result

        #try:
        self.init_atl_h5()
        # except Exception as e:
        #     utils.echo_error_msg(
        #         f'could not initialize data {e}'
        #     )
        #     self.close_atl_h5()
        #     return

        ## fetch and process buildings, if wanted
        #this_bing = None
        if self.want_buildings:
            if isinstance(self.want_buildings, bool):
                this_bing = self.process_buildings(
                    self.fetch_buildings(verbose=self.verbose),
                    verbose=self.verbose
                )
            elif isinstance(self.want_buildings, list):
                this_bing = self.want_buildings         

        ## fetch and process watermask, if wanted
        #this_wm = None
        if self.want_watermask:
            if isinstance(self.want_watermask, bool):
                this_wm = self.process_coastline(
                    self.fetch_coastline(chunks=False, verbose=self.verbose),
                    return_geom=True,
                    verbose=self.verbose
                )
            elif isinstance(self.want_watermask, list):
                this_wm = self.want_watermask

        ## fetch and process watermask, if wanted
        #this_wm = None
        # if self.want_watermask:
        #     this_wm = self._load_gmrt()
        #     # if isinstance(self.want_watermask, bool):
        #     #     this_wm = self.process_coastline(
        #     #         self.fetch_coastline(chunks=False, verbose=self.verbose),
        #     #         return_geom=True,
        #     #         verbose=self.verbose
        #     #     )
        #     # elif isinstance(self.want_watermask, list):
        #     #     this_wm = self.want_watermask                
                

        ## parse through the icesat2 file by laser number
        for i in range(1, 4):
            for orient in range(2):
                #try:
                ## selects the strong beams only [we can include weak beams later on]
                if f'gt{i}{self.orientDict[orient]}' not in list(self.atl03_f.keys()):
                    #utils.echo_warning_msg(f'could not find laser gt{i}{self.orientDict[orient]} in {self.atl03_fn}')
                    continue

                if 'heights' not in list(self.atl03_f[f'gt{i}{self.orientDict[orient]}'].keys()):
                    #utils.echo_warning_msg(f'laser gt{i}{self.orientDict[orient]} in {self.atl03_fn} has no heights')
                    continue
                
                dataset = self.read_atl_data('{}'.format(i), orientation=orient)
                if dataset is None or len(dataset) == 0:
                    continue

                ## keep only photons with confidence levels mentioned
                ## in `self.confidence_levels`
                if len(self.confidence_levels) > 0:
                    dataset = dataset[
                        (np.isin(dataset['confidence'], self.confidence_levels))
                    ]

                ## reduce the dataset to the input region for faster masking
                #dataset = dataset[(dataset['longitude'] >= self.region.xmin) & (dataset['longitude'] <= self.region.xmax)]
                #dataset = dataset[(dataset['latitude'] >= self.region.ymin) & (dataset['latitude'] <= self.region.ymax)]
                    
                if dataset is None or len(dataset) == 0:
                    continue
                
                ## re-classify photons based on buildings/watermask/bathymetry
                if self.want_buildings and this_bing is not None:
                    #dataset = self.classify_buildings(dataset, this_bing)
                    dataset = self.classify_by_mask_geoms(dataset, mask_geoms=this_bing, classification=7)

                #utils.echo_msg(f'this_wm: {this_wm}')
                if self.want_watermask and this_wm is not None:
                    #dataset = self.classify_water(dataset, this_wm)
                    dataset = self.classify_by_mask_geoms(dataset, mask_geoms=this_wm, classification=41, except_classes=[40])

                if dataset is None or len(dataset) == 0:
                    continue
                
                # ## bathymetry is classified in `read_atl_data` using ATL24 now...
                # if self.want_bathymetry:
                #     dataset = self.classify_bathymetry(dataset)
                    
                # if dataset is None or len(dataset) == 0:
                #     continue

                ## keep only photons with a classification mentioned in `self.classes`
                if len(self.classes) > 0:
                    dataset = dataset[
                        (np.isin(dataset['ph_h_classed'], self.classes))
                    ]

                if dataset is None or len(dataset) == 0:
                    continue

                ## rename the x,y,z columns for `transform_and_yield_points`
                dataset.rename(
                    columns={
                        'longitude': 'x', 'latitude': 'y', 'photon_height': 'z'
                    },
                    inplace=True
                )
                yield(dataset)

        self.close_atl_h5()

        
    def fetch_data(self, fetches_module, check_size=True):
        this_fetches = fetches.FetchesFactory(
            mod=fetches_module,
            src_region=self.region,
            verbose=self.verbose,
            outdir=self.cache_dir,
            callback=fetches.fetches_callback
        )._acquire_module()        
        this_fetches.run()
        fr = fetches.fetch_results(this_fetches, check_size=check_size)
        fr.daemon = True
        fr.start()
        fr.join()

        return(fr)


    def _load_gmrt(self):
        """GMRT - Global low-res.
        """
        
        this_gmrt = self.fetch_data('gmrt:layer=topo')        
        gmrt_result = this_gmrt.results[0]
        if gmrt_result[-1] == 0:
            gmrt_tif = gmrt_result[1]

            return(gmrt_tif)

        
    def fetch_atlxx(self, short_name='ATL08'):
        """fetch an associated ATLxx file"""

        atlxx_filter = short_name + '_' + '_'.join(utils.fn_basename2(self.atl03_fn).split('_')[1:4])
        atlxx_filter_no_version = short_name + '_' + '_'.join(utils.fn_basename2(self.atl03_fn).split('_')[1:3])
        #utils.echo_msg(f'{os.path.join(self.cache_dir, atlxx_filter)}*')
        atlxx_paths = glob.glob(f'{os.path.join(self.cache_dir, atlxx_filter)}*.h5')
        if len(atlxx_paths) > 0:
            return(atlxx_paths[0])
        else:
            atlxx_paths = glob.glob(f'{os.path.join(self.cache_dir, atlxx_filter_no_version)}*.h5')
            if len(atlxx_paths) > 0:
                return(atlxx_paths[0])

        atlxx_filter = '_'.join(utils.fn_basename2(self.atl03_fn).split('_')[1:4])
        atlxx_filter_no_version = '_'.join(utils.fn_basename2(self.atl03_fn).split('_')[1:3])            
        this_atlxx = earthdata.IceSat2(
            src_region=None,
            verbose=self.verbose,
            outdir=self.cache_dir,
            short_name=short_name,
            filename_filter=atlxx_filter,
            version='',
            #subset=True if 'subset' in self.fn else False,
        )
        this_atlxx.run()
        if len(this_atlxx.results) == 0:
            #atlxx_filter = '_'.join(utils.fn_basename2(self.atl03_fn).split('_')[1:3])
            this_atlxx = earthdata.IceSat2(
                src_region=None,
                verbose=self.verbose,
                outdir=self.cache_dir,
                short_name=short_name,
                filename_filter=atlxx_filter_no_version,
                version='',
            )
            this_atlxx.run()
            if len(this_atlxx.results) == 0:
                utils.echo_warning_msg(
                    (f'could not locate associated {short_name} '
                     f'file for {atlxx_filter}')
                )
                return(None)
        else:
            if this_atlxx.fetch_entry(
                    this_atlxx.results[0], check_size=True
            ) == 0:
                return(os.path.join(
                    this_atlxx._outdir, this_atlxx.results[0]['dst_fn']
                ))

            
    def apply_atl08_classifications(
            self, ph_h_classed, laser, segment_id, segment_index_dict, ph_segment_ids
    ):
        if self.atl08_f is not None:
            if laser in self.atl08_f.keys():
                ## classed flag (signal_photons)
                atl08_classed_pc_flag  = self.atl08_f[
                    f'/{laser}/signal_photons/classed_pc_flag'
                ][...,]
                atl08_ph_segment_id = self.atl08_f[
                    f'/{laser}/signal_photons/ph_segment_id'
                ][...,] # photon src 20 m seg id
                atl08_classed_pc_indx = self.atl08_f[
                    f'/{laser}/signal_photons/classed_pc_indx'
                ][...,]

                ## set the classifications from ATL08
                atl08_segment_id_msk = [
                    True if x in segment_id else False for x in atl08_ph_segment_id
                ]
                atl08_ph_segment_indx = np.array(
                    list(
                        map((lambda pid: segment_index_dict[pid]),
                            atl08_ph_segment_id[atl08_segment_id_msk])
                    )
                )
                atl08_ph_index = np.array(
                    atl08_ph_segment_indx + (atl08_classed_pc_indx[atl08_segment_id_msk] - 1),
                    dtype=int
                )
                class_mask = atl08_ph_index < len(ph_segment_ids)
                ph_h_classed[atl08_ph_index[class_mask]] \
                    = atl08_classed_pc_flag[atl08_segment_id_msk][class_mask]

                #utils.echo_msg_bold(np.count_nonzero(ph_h_classed==1))
            
        return(ph_h_classed)


    def apply_atl24_classifications(
            self, ph_h_classed, ph_h_bathy_conf, latitude, longitude, photon_h, photon_h_meantide,
            photon_h_geoid, laser, geoseg_beg, geoseg_end, ph_segment_ids
    ):
        if self.atl24_f is not None:
            ## some atl24 files don't have all the lasers, so make sure
            ## it exists before proceding
            if laser in self.atl24_f.keys():            
                lat_ph = self.atl24_f[f'{laser}/lat_ph'][...,]
                lon_ph = self.atl24_f[f'{laser}/lon_ph'][...,]
                ortho_h = self.atl24_f[f'{laser}/ortho_h'][...,]
                ellipse_h = self.atl24_f[f'{laser}/ellipse_h'][...,]
                surface_h = self.atl24_f[f'{laser}/surface_h'][...,]
                class_ph = self.atl24_f[f'{laser}/class_ph'][...,]
                index_seg = self.atl24_f[f'{laser}/index_seg'][...,]
                index_ph = self.atl24_f[f'{laser}/index_ph'][...,]
                conf_ph = self.atl24_f[f'{laser}/confidence'][...,]

                #if 'subset' in self.atl03_fn:
                ## determine the original `geosed_id` to determine the subset indices
                orig_segment_id = np.arange(geoseg_beg, geoseg_end+1, 1)
                orig_segment_id_msk = np.isin(orig_segment_id, ph_segment_ids)
                index_seg_orig = orig_segment_id[index_seg]
                segment_id_msk = np.isin(index_seg_orig, ph_segment_ids)
                if len(index_ph[segment_id_msk]) > 0:
                    ## check for previous index
                    index_ph = index_ph - np.min(index_ph[segment_id_msk])
                    index_msk = (index_ph[segment_id_msk] > 0) & (index_ph[segment_id_msk] < np.max(index_ph[segment_id_msk]))
                    class_msk = class_ph[segment_id_msk][index_msk] >= 40
                    ph_msk = (class_msk)
                    if self.min_bathy_confidence is not None:
                        confidence_msk = conf_ph[segment_id_msk][index_msk] >= self.min_bathy_confidence
                        ph_msk = (class_msk) & (confidence_msk)

                    ph_h_classed[index_ph[segment_id_msk][index_msk][ph_msk]] = class_ph[segment_id_msk][index_msk][ph_msk]
                    ph_h_bathy_conf[index_ph[segment_id_msk][index_msk][ph_msk]] = conf_ph[segment_id_msk][index_msk][ph_msk]
                    # we also need to change the lon/lat/height values to the
                    # updated/refracted bathymetry values (we'll just do it to class 40)
                    class_msk = class_ph[segment_id_msk][index_msk] == 40
                    ph_class = (class_msk)
                    if self.min_bathy_confidence is not None:
                        ph_msk = (class_msk) & (confidence_msk)

                    longitude[index_ph[segment_id_msk][index_msk][ph_msk]] = lon_ph[segment_id_msk][index_msk][ph_msk]
                    latitude[index_ph[segment_id_msk][index_msk][ph_msk]] = lat_ph[segment_id_msk][index_msk][ph_msk]
                    photon_h[index_ph[segment_id_msk][index_msk][ph_msk]] = ellipse_h[segment_id_msk][index_msk][ph_msk]
                    photon_h_meantide[index_ph[segment_id_msk][index_msk][ph_msk]] = surface_h[segment_id_msk][index_msk][ph_msk]
                    photon_h_geoid[index_ph[segment_id_msk][index_msk][ph_msk]] = ortho_h[segment_id_msk][index_msk][ph_msk]
                    #conf[index_ph[segment_id_msk][index_msk][ph_msk]] = conf_ph[segment_id_msk][index_msk][ph_msk]

        return(ph_h_classed, ph_h_bathy_conf, latitude, longitude, photon_h, photon_h_meantide, photon_h_geoid)
    
    
    def read_atl_data(self, laser_num, orientation=None):
        """Read data from an ATL03 file

        Adapted from 'cshelph' https://github.com/nmt28/C-SHELPh.git 
        and 'iVert' https://github.com/ciresdem/ivert.git

        laser_num is 1, 2 or 3
        surface is 'mean_tide', 'geoid' or 'ellipsoid'
        """

        if orientation is None:
            orientation = self.atl03_f['/orbit_info/sc_orient'][0]
            
        ## selects the strong beams only [we can include weak beams later on]
        #orientDict = {0:'l', 1:'r', 21:'error'}
        laser = 'gt' + laser_num + self.orientDict[orientation]

        ## for 'subsets', where heights don't come through
        if 'heights' not in self.atl03_f['/{}'.format(laser)].keys():
            return(None)
        
        ## Read in the required photon level data
        photon_h = self.atl03_f['/' + laser + '/heights/h_ph'][...,]
        latitude = self.atl03_f['/' + laser + '/heights/lat_ph'][...,]
        longitude = self.atl03_f['/' + laser + '/heights/lon_ph'][...,]
        ph_count = self.atl03_f['/' + laser + '/heights/ph_id_count'][...,]
        conf = self.atl03_f['/' + laser + '/heights/signal_conf_ph/'][...,0]
        qual = self.atl03_f['/' + laser + '/heights/quality_ph/'][...,0]
        dist_ph_along = self.atl03_f['/' + laser + '/heights/dist_ph_along'][...,]
        delta_time = self.atl03_f['/' + laser + '/heights/delta_time'][...,]
        this_N = latitude.shape[0]

        ## Read in the geolocation level data
        segment_ph_cnt = self.atl03_f['/' + laser + '/geolocation/segment_ph_cnt'][...,]
        segment_id = self.atl03_f['/' + laser + '/geolocation/segment_id'][...,]
        segment_dist_x = self.atl03_f['/' + laser + '/geolocation/segment_dist_x'][...,]
        ref_elev = self.atl03_f['/' + laser + '/geolocation/ref_elev'][...,]
        ref_azimuth = self.atl03_f['/' + laser + '/geolocation/ref_azimuth'][...,]
        ref_photon_index = self.atl03_f['/' + laser + '/geolocation/reference_photon_index'][...,]
        ph_index_beg = self.atl03_f['/' + laser + '/geolocation/ph_index_beg'][...,]
        altitude_sc = self.atl03_f['/' + laser + '/geolocation/altitude_sc'][...,]
        surf_type = self.atl03_f['/' + laser + '/geolocation/surf_type'][...,]
        
        ## Read in the geophysical correction
        photon_geoid = self.atl03_f['/' + laser + '/geophys_corr/geoid'][...,]
        photon_geoid_f2m = self.atl03_f['/' + laser + '/geophys_corr/geoid_free2mean'][...,]
        photon_dem_h = self.atl03_f['/' + laser + '/geophys_corr/dem_h'][...,]

        ## Read in the geosegment info for subsets
        geoseg_beg = self.atl03_f['ancillary_data/start_geoseg'][...,][0]
        geoseg_end = self.atl03_f['ancillary_data/end_geoseg'][...,][0]
                    
        ## Create a dictionary with (segment_id --> index into ATL03 photons)
        ## lookup pairs, for the starting photon of each segment
        #segment_indices = ph_index_beg+(ref_photon_index-1)
        #segment_index_dict = dict(zip(segment_id, ph_index_beg+ref_photon_index-1))
        segment_indices = np.concatenate(([0], np.cumsum(segment_ph_cnt)[:-1]))
        segment_index_dict = dict(zip(segment_id, segment_indices))
        ph_segment_ids = segment_id[
            np.searchsorted(segment_indices, np.arange(0.5, len(photon_h), 1))-1
        ]
        
        ## Compute the total along-track distances.
        segment_dist_dict = dict(zip(segment_id, segment_dist_x))

        ## Determine where in the array each segment index needs to look.
        ph_segment_dist_x = np.array(
            list(map((lambda pid: segment_dist_dict[pid]), ph_segment_ids))
        )
        dist_x = ph_segment_dist_x + dist_ph_along

        ## meantide/geoid/dem heights
        h_geoid_dict = dict(zip(segment_id, photon_geoid))
        ph_h_geoid = np.array(list(map((lambda pid: h_geoid_dict[pid]), ph_segment_ids)))        
        h_meantide_dict = dict(zip(segment_id, photon_geoid_f2m))
        ph_h_meantide = np.array(list(map((lambda pid: h_meantide_dict[pid]), ph_segment_ids)))
        photon_h_geoid = photon_h - ph_h_geoid
        photon_h_meantide = photon_h - (ph_h_geoid + ph_h_meantide)

        h_dem_dict = dict(zip(segment_id, photon_dem_h))
        ph_h_dem = np.array(list(map((lambda pid: h_dem_dict[pid]), ph_segment_ids)))
        photon_h_dem = ph_h_dem - (ph_h_geoid + ph_h_meantide)
        
        ## setup classifications asn bathy confidence
        ph_h_classed = np.zeros(photon_h.shape)
        ph_h_classed[:] = -1

        ph_h_bathy_conf = np.zeros(photon_h.shape)
        ph_h_bathy_conf[:] = -1

        ## append the laser to each record
        laser_arr = np.empty(photon_h.shape, dtype='object')
        laser_arr[:] = laser

        ## append the filename to each record
        fn_arr = np.empty(photon_h.shape, dtype='object')
        fn_arr[:] = self.fn
        
        ## ref values
        h_ref_elev_dict = dict(zip(segment_id, ref_elev))
        ph_ref_elev = np.array(
            list(map((lambda pid: h_ref_elev_dict[pid]), ph_segment_ids))
        )#.astype(float)        
        h_ref_azimuth_dict = dict(zip(segment_id, ref_azimuth))
        ph_ref_azimuth = np.array(
            list(map((lambda pid: h_ref_azimuth_dict[pid]), ph_segment_ids))
        )#.astype(float)
        h_altitude_sc_dict = dict(zip(segment_id, altitude_sc))
        ph_altitude_sc = np.array(
            list(map((lambda pid: h_altitude_sc_dict[pid]), ph_segment_ids))
        )#.astype(float)
        h_surf_type_dict = dict(zip(segment_id, surf_type))
        ph_surf_type = np.array(
            list(map((lambda pid: h_surf_type_dict[pid]), ph_segment_ids))
        )#.astype(float)
        
        ## classify the photons
        ## Read in the ATL08 data and classify points based on the ATL08 classifications
        ph_h_classed = self.apply_atl08_classifications(
            ph_h_classed, laser, segment_id, segment_index_dict, ph_segment_ids
        )
        
        ## classify phtons as water surface if dem_h is below zero
        #ph_h_classed[photon_h_dem < 0] = 41
        
        ## classify bathymetry and water-surface points based on ATL24
        ## classifications...
        ph_h_classed, ph_h_bathy_conf, latitude, longitude, photon_h, photon_h_meantide, photon_h_geoid = self.apply_atl24_classifications(
            ph_h_classed, ph_h_bathy_conf, latitude, longitude, photon_h, photon_h_meantide,
            photon_h_geoid, laser, geoseg_beg, geoseg_end, ph_segment_ids
        )

        ## classify inland water based on ATL13 classifications
        # if self.atl13_f is not None:
        #     atl13_refid = self.atl13_f['/' + laser + '/segment_id_beg'][...,]
        #     ph_h_classed[atl13_refid] = 44

        ## classify water based on ATL12 classifications
        # if self.atl13_f is not None:
        #     atl13_refid = self.atl13_f['/' + laser + '/segment_id_beg'][...,]
        #     ph_h_classed[atl13_refid] = 44
        
        ## classify water points based on the 'surf_type' parameter to
        ## remove ATL08 bare-earth clasification over water
        ph_surf_type = [x[1] for x in ph_surf_type]
            
        ## set the photon height, either 'mean_tide' or 'geoid', else ellipsoid
        if self.water_surface == 'mean_tide':
            ph_height = photon_h_meantide
        elif self.water_surface == 'geoid':
            ph_height = photon_h_geoid
        else:
            ph_height = photon_h
            
        ## create the pandas dataframe            
        dataset = pd.DataFrame(
            {'latitude': latitude,
             'longitude': longitude,
             'photon_height': ph_height,
             'laser': laser_arr,
             'fn': fn_arr,
             'confidence': conf,
             'ref_elevation':ph_ref_elev,
             'ref_azimuth':ph_ref_azimuth,
             'ref_sat_alt':ph_altitude_sc,
             'delta_time':delta_time,
             'bathy_confidence':ph_h_bathy_conf,
             'ph_surf_type':ph_surf_type,
             'photon_h_dem':photon_h_dem,
             'ph_h_classed': ph_h_classed},
            columns=['latitude', 'longitude', 'photon_height', 'laser', 'fn',
                     'confidence', 'ref_elevation', 'ref_azimuth', 'ref_sat_alt',
                     'delta_time', 'bathy_confidence', 'ph_surf_type', 'photon_h_dem', 'ph_h_classed']
        )

        ## Process extra columns specified in `self.columns`
        for col in self.columns.keys():
            try:
                if 'gtx' in col:
                    _col = col.replace('/gtx', '/' + laser)
                    col_arr = self.atl03_f[_col][...,]                
                    if 'heights' not in _col: 
                        col_dict = dict(zip(segment_id, col_arr))
                        col_arr = np.array(
                            list(map((lambda pid: col_dict[pid]), ph_segment_ids))
                        )
                else:
                    col_arr = np.empty(photon_h.shape, dtype='object')
                    col_arr[:] = self.atl03_f[col][...,]

                extra_dataset = pd.DataFrame(
                    {self.columns[col]: col_arr},
                    columns = [self.columns[col]]
                )
                dataset = dataset.join(extra_dataset)
            except:
                utils.echo_warning_msg(
                    f'could not find and/or process {col}'
                )

        return(dataset)        

    
    def append_atl24_to_dataset(dataset):
        if self.atl24_f is not None:
            orientDict = {0:'l', 1:'r', 21:'error'}
            laser = 'gt' + f'{i}' + orientDict[orient]
            atl24_classed_pc_flag  = self.atl24_f[
                f'/{laser}/class_ph'
            ][...,]
            atl24_lat_ph = self.atl24_f[
                f'/{laser}/lat_ph'
            ][...,]
            atl24_lon_ph = self.atl24_f[
                f'/{laser}/lon_ph'
            ][...,]
            atl24_suface_h = self.atl24_f[
                f'/{laser}/surface_h'
            ][...,]
            atl24_ortho_h = self.atl24_f[
                f'/{laser}/ortho_h'
            ][...,]
            atl24_ellipse_h = self.atl24_f[
                f'/{laser}/ellipse_h'
            ][...,]

            atl24_ph_segment_id = self.atl24_f[
                f'/{laser}/index_seg'
            ][...,] # photon src 20 m seg id
            atl24_classed_pc_indx = self.atl24_f[
                f'/{laser}/index_ph'
            ][...,]
            atl24_conf = self.atl24_f[
                f'/{laser}/confidence'
            ][...,]

            if self.water_surface == 'mean_tide':
                atl24_ph_height = atl24_surface_h
            elif self.water_surface == 'geoid':
                atl24_ph_height = atl24_ortho_h
            else:
                atl24_ph_height = atl24_ellipse_h

            ## append the laser to each record

            laser_arr = np.empty(atl24_ortho_h.shape, dtype='object')
            laser_arr[:] = laser

            ## append the filename to each record
            fn_arr = np.empty(atl24_ortho_h.shape, dtype='object')
            fn_arr[:] = self.fn

            class_mask = (atl24_classed_pc_flag == 40) | (atl24_classed_pc_flag == 41)
            if self.region is not None:
                lat_mask = (atl24_lat_ph >= self.region.ymin) & (atl24_lat_ph <= self.region.ymax)
                lon_mask = (atl24_lon_ph >= self.region.xmin) & (atl24_lon_ph <= self.region.xmax)
                dataset_mask = class_mask & lat_mask & lon_mask
            else:
                dataset_mask = class_mask

            atl24_dataset = pd.DataFrame(
                {'latitude': atl24_lat_ph[dataset_mask],
                 'longitude': atl24_lon_ph[dataset_mask],
                 'photon_height': atl24_ph_height[dataset_mask],
                 'laser': laser_arr[dataset_mask],
                 'fn': fn_arr[dataset_mask],
                 'confidence': atl24_conf[dataset_mask],
                 'ph_h_classed': atl24_classed_pc_flag[dataset_mask]},
                columns=['latitude', 'longitude', 'photon_height', 'laser', 'fn',
                         'confidence', 'ph_h_classed']
            )
            for column in dataset.columns:
                if column not in atl24_dataset.columns:
                    tmp_data = np.zeros(atl24_ortho_h.shape)
                    tmp_data[:] = -1

                    atl24_dataset[column] = tmp_data

            dataset = pd.concat([dataset, atl24_dataset], ignore_index=True)
    
    
    def classify_bathymetry(
            self, dataset, thresh=95, min_buffer=-40, max_buffer=5,
            start_lat=False, end_lat=False, lat_res=10, h_res=.5,
            surface_buffer=-.5, water_temp=None
    ):
        """Classify bathymetry in an ATL03 file. 

        This uses C-Shelph to locate, extract and process bathymetric photons.
        This function is adapted from the C-Shelph CLI
        Depreciated once ATL24 are online
        """
        
        water_temp = utils.float_or(water_temp)
        epsg_code = cshelph.convert_wgs_to_utm(
            dataset.latitude.iloc[0], dataset.longitude.iloc[0]
        )
        epsg_num = int(epsg_code.split(':')[-1])
        utm_proj = pyproj.Proj(epsg_code)
        lon_utm, lat_utm = utm_proj(dataset.longitude, dataset.latitude)
                       
        ## Aggregate data into dataframe
        dataset_sea = pd.DataFrame(
            {'latitude': lat_utm,
             'longitude': lon_utm,
             'photon_height': dataset.photon_height,
             'laser': dataset.laser,
             'fn': dataset.fn,
             'confidence': dataset.confidence,
             'ref_elevation': dataset.ref_elevation,
             'ref_azimuth': dataset.ref_azimuth,
             'ref_sat_alt': dataset.ref_sat_alt,
             'ph_h_classed': dataset.ph_h_classed},
            columns=['latitude', 'longitude', 'photon_height',
                     'laser', 'fn', 'confidence', 'ref_elevation',
                     'ref_azimuth', 'ref_sat_alt', 'ph_h_classed']
        )
        dataset_sea1 = dataset_sea[
            (dataset_sea['photon_height'] > min_buffer) \
            & (dataset_sea['photon_height'] < max_buffer)
        ]
        dataset_sea1 = dataset_sea1[(dataset_sea1['ph_h_classed'] == 5)]
        
        if self.region is not None:
            xyz_region = self.region.copy()
            xyz_region.epsg = 'epsg:4326'
            xyz_region.warp('epsg:{}'.format(epsg_num))
            dataset_sea1 = dataset_sea1[
                (dataset_sea1['latitude'] > xyz_region.ymin) \
                & (dataset_sea1['latitude'] < xyz_region.ymax)
            ]
            dataset_sea1 = dataset_sea1[
                (dataset_sea1['longitude'] > xyz_region.xmin) \
                & (dataset_sea1['longitude'] < xyz_region.xmax)
            ]

        if len(dataset_sea1) > 0:
            binned_data_sea = cshelph.bin_data(dataset_sea1, lat_res, h_res)
            if binned_data_sea is not None:
                if water_temp is None:
                    try:
                        # water_temp = cshelph.get_water_temp(
                        #     self.fn, latitude, longitude
                        # )
                        ## water_temp via fetches instead of earthaccess
                        water_temp = self.get_water_temp()
                    except Exception as e:
                        # utils.echo_warning_msg(
                        #     'NO SST PROVIDED OR RETRIEVED: 20 degrees C assigned'
                        # )
                        water_temp = 20

                #utils.echo_msg('water temp is {}'.format(water_temp))
                sea_height = cshelph.get_sea_height(binned_data_sea, surface_buffer)
                sea_height1 = cshelph.get_bin_height(binned_data_sea, 60, surface_buffer)
                
                if sea_height is not None:
                    med_water_surface_h = np.nanmedian(sea_height) #* -1
                    
                    ## Correct for refraction
                    ref_x, ref_y, ref_z, ref_conf, raw_x, raw_y, raw_z, ph_ref_azi, ph_ref_elev = cshelph.refraction_correction(
                        water_temp, med_water_surface_h, 532, dataset_sea1.ref_elevation,
                        dataset_sea1.ref_azimuth, dataset_sea1.photon_height,
                        dataset_sea1.longitude, dataset_sea1.latitude, dataset_sea1.confidence,
                        dataset_sea1.ref_sat_alt
                    )
                    depth = med_water_surface_h - ref_z

                    # Create new dataframe with refraction corrected data
                    dataset_bath = pd.DataFrame(
                        {'latitude': raw_y, 'longitude': raw_x,
                         'cor_latitude':ref_y, 'cor_longitude':ref_x,
                         'cor_photon_height':ref_z, 'photon_height': raw_z,
                         'confidence':ref_conf, 'depth':depth},
                        columns=['latitude', 'longitude', 'photon_height',
                                 'cor_latitude','cor_longitude',
                                 'cor_photon_height', 'confidence',
                                 'depth']
                    )

                    # Bin dataset again for bathymetry
                    if len(dataset_bath) > 0:
                        binned_data = cshelph.bin_data(dataset_bath, lat_res, h_res)
                        if binned_data is not None:
                            # Find bathymetry
                            bath_height, geo_df = cshelph.get_bath_height(
                                binned_data, thresh, med_water_surface_h, h_res
                            )
                            if bath_height is not None:
                                transformer = pyproj.Transformer.from_crs(
                                    "EPSG:"+str(epsg_num),
                                    "EPSG:4326",
                                    always_xy=True
                                )
                                lon_wgs84, lat_wgs84 = transformer.transform(
                                    geo_df.longitude.values, geo_df.latitude.values
                                )
                                bath_height = [x for x in bath_height if ~np.isnan(x)]
                                
                                for n, id in enumerate(geo_df.ids.values):
                                    dataset.at[id, 'ph_h_classed'] = 40
                                    dataset.at[id, 'latitude'] = lat_wgs84[n]
                                    dataset.at[id, 'longitude'] = lon_wgs84[n]
                                    dataset.at[id, 'photon_height'] = geo_df.depth.values[n] * -1
                            
                                return(dataset)
            
        return(dataset)

    
    def _vectorize_df(self, dataset):
        """Make a point vector OGR DataSet Object from a pandas dataframe

        This is to allow spatial filtering for watermask and buildings.
        """

        dst_ogr = '{}'.format('icesat_dataframe')
        ogr_ds = gdal.GetDriverByName('Memory').Create(
            '', 0, 0, 0, gdal.GDT_Unknown
        )
        layer = ogr_ds.CreateLayer(
            dst_ogr, geom_type=ogr.wkbPoint
        )
        fd = ogr.FieldDefn('index', ogr.OFTInteger)
        layer.CreateField(fd)
        # for x in dataset.columns.to_list():
        #     fd = ogr.FieldDefn(x, ogr.OFTReal)
        #     fd.SetWidth(12)
        #     fd.SetPrecision(8)
        #     layer.CreateField(fd)
            
        f = ogr.Feature(feature_def=layer.GetLayerDefn())        
        with tqdm(desc='vectorizing dataframe', leave=False) as pbar:
            for index, this_row in dataset.iterrows():
                pbar.update()
                f.SetField(0, index)
                #[f.SetField(n+1, this_row[x]) for n,x in enumerate(dataset.columns.to_list())]
                g = ogr.CreateGeometryFromWkt(
                    f'POINT ({this_row["longitude"]} {this_row["latitude"]})'
                )
                f.SetGeometryDirectly(g)
                layer.CreateFeature(f)
            
        return(ogr_ds)

    
    def classify_buildings(self, dataset, these_bldgs):
        """classify building photons using BING building footprints 
        """

        ## vectorize the icesate2 photons
        ogr_df = self._vectorize_df(dataset)
        
        os.environ["OGR_OSM_OPTIONS"] = "INTERLEAVED_READING=YES"
        os.environ["OGR_OSM_OPTIONS"] = "OGR_INTERLEAVED_READING=YES"
        with tqdm(
                total=len(these_bldgs),
                desc='classifying building photons',
                leave=False
        ) as pbar:
            for n, bldg_geom in enumerate(these_bldgs):                
                pbar.update()
                icesat_layer = ogr_df.GetLayer()
                icesat_layer.SetSpatialFilter(bldg_geom)
                #[dataset.at[f.GetField('index'), 'ph_h_classed'] = 7 for f in icesat_layer]
                for f in icesat_layer:
                    idx = f.GetField('index')
                    dataset.at[idx, 'ph_h_classed'] = 7

                icesat_layer.SetSpatialFilter(None)

        ogr_df = None
        return(dataset)

    
    def classify_water(self, dataset, these_wms):
        """classify water photons using OSM coastline 
        """

        ## vectorize the photons
        ogr_df = self._vectorize_df(dataset)
        
        os.environ["OGR_OSM_OPTIONS"] = "INTERLEAVED_READING=YES"
        os.environ["OGR_OSM_OPTIONS"] = "OGR_INTERLEAVED_READING=YES"
        with tqdm(
                total=len(these_wms),
                desc='classifying water photons',
                leave=False
        ) as pbar:
            for n, wm_geom in enumerate(these_wms):                
                pbar.update()
                icesat_layer = ogr_df.GetLayer()
                icesat_layer.SetSpatialFilter(wm_geom)
                #[dataset.at[f.GetField('index'), 'ph_h_classed'] = 7 for f in icesat_layer]
                for f in icesat_layer:
                    idx = f.GetField('index')
                    #bathy_mask = dataset.at[idx, 'ph_h_classed'] == 40
                    #utils.echo_msg_bold(idx)
                    #utils.echo_msg_bold(dataset.at[idx, 'ph_h_classed'])
                    if dataset.at[idx, 'ph_h_classed'] != 40:
                        dataset.at[idx, 'ph_h_classed'] = 41

                icesat_layer.SetSpatialFilter(None)

        ogr_df = None
        return(dataset)

    def classify_by_mask_geoms(self, dataset, mask_geoms=[], classification=-1, except_classes=[]):
        """classify water photons using OSM coastline 
        """

        ## vectorize the photons
        ogr_df = self._vectorize_df(dataset)
        
        os.environ["OGR_OSM_OPTIONS"] = "INTERLEAVED_READING=YES"
        os.environ["OGR_OSM_OPTIONS"] = "OGR_INTERLEAVED_READING=YES"
        with tqdm(
                total=len(mask_geoms),
                desc='classifying water photons',
                leave=False
        ) as pbar:
            for n, wm_geom in enumerate(mask_geoms):                
                pbar.update()
                icesat_layer = ogr_df.GetLayer()
                icesat_layer.SetSpatialFilter(wm_geom)
                #utils.echo_msg(len(icesat_layer))
                for f in icesat_layer:
                    idx = f.GetField('index')
                    #except_mask = dataset.at[idx, 'ph_h_classed']
                    if dataset.at[idx, 'ph_h_classed'] not in except_classes:
                        dataset.at[idx, 'ph_h_classed'] = classification

                icesat_layer.SetSpatialFilter(None)

        ogr_df = None
        return(dataset)


class IceSat2_ATL03(ElevationDataset):
    """IceSat2 data from NASA

    Parameters:
    
    water_surface: 'geoid' # this is the vertical datum, can be 'geoid', 
    'ellipsoid' or 'mean_tide'

    classes: None # return only data with the specified classes, e.g. '2/3/4'
    # return only data with the specified confidence levels, e.g. '2/3/4'
    confidence_levels: None 
    columns: {} # the additional columns to export in transform_and_yield_points 
             {'/h5/atl/path', column_name}
    classify_buildings: True # classify buildings using BING BFP
    classify_water: True # classify water using OSM
    reject_failed_qa: True # skip granules that failed QA
    min_bathy_confidence: <0-1> the minimum passable ATL24 confidence level

    Classes:
    -1 - no classification (ATL08)
    0 - noise / atmosphere (ATL08)
    1 - ground surface (ATL08)
    2 - canopy (ATL08)
    3 - canopy top (ATL08)
    40 - bathymetry floor surface (CShelph, ATL24)
    41 - bathymetry water surface (OSM coastline, ATL24)
    6 - ice surface (ATL06) 
        (unused for now, just planning ahead for possible future 
        ATL06 integration)
    7 - built structure (OSM or Bing)
    8 - "urban" (WSF, if used)
    9 - inland water surface

    Confidence Levels:
    0, 1, 2, 3, 4
    """

    def __init__(self,
                 water_surface='geoid',
                 classes=None,
                 confidence_levels='2/3/4',
                 columns={},
                 reject_failed_qa=True,
                 classify_buildings=True,
                 classify_water=True,
                 append_atl24=False,
                 min_bathy_confidence=None,
                 **kwargs):
        super().__init__(**kwargs)
        self.data_format = 303
        self.water_surface = water_surface
        if self.water_surface not in ['mean_tide', 'geoid', 'ellipsoid']:
            self.water_surface = 'mean_tide'

        self.classes = [int(x) for x in classes.split('/')] \
            if classes is not None \
               else []
        self.confidence_levels = [int(x) for x in confidence_levels.split('/')] \
            if confidence_levels is not None \
               else []
        self.columns = columns
        #self.columns = {'/gtx/heights/delta_time': 'delta_time_b'}
        self.atl_fn = None
        self.classify_buildings = classify_buildings
        self.classify_water = classify_water
        self.reject_failed_qa = reject_failed_qa
        self.append_atl24 = append_atl24
        self.min_bathy_confidence = utils.float_or(min_bathy_confidence)

        self.orientDict = {0:'l', 1:'r', 21:'error'}
        self.lasers = ['gt1l', 'gt2l', 'gt3l', 'gt1r', 'gt2r', 'gt3r']
        

    def generate_inf(self):
        if self.src_srs is None:
            self.infos.src_srs = 'epsg:4326+3855'
        else:
            self.infos.src_srs = self.src_srs

        these_regions = []
        numpts = 0
        with h5.File(self.fn, 'r') as f:
            for b in range(1, 4):
                for p in ['l', 'r']:
                    if f'gt{b}{p}' not in list(f.keys()):
                        continue

                    if 'heights' not in list(f[f'gt{b}{p}'].keys()):
                        continue

                    y = f[f'gt{b}{p}/heights/lat_ph'][...,]
                    x = f[f'gt{b}{p}/heights/lon_ph'][...,]
                    z = f[f'gt{b}{p}/heights/h_ph'][...,]
                    this_region = regions.Region(src_srs=self.infos.src_srs).from_list(
                        [np.min(x), np.max(x), np.min(y), np.max(y), np.min(z), np.max(z)]
                    )
                    these_regions.append(this_region)
                    numpts += len(x)

        region_cnt = len(these_regions)
        if region_cnt > 0:
            merged_region = these_regions[0].copy()
            for r in range(1, region_cnt):
                merged_region = regions.regions_merge(merged_region, these_regions[r])

            self.infos.minmax = merged_region.export_as_list(include_z=True)
            self.infos.wkt = merged_region.export_as_wkt()
            
        self.infos.numpts = numpts 
        return(self.infos)        


    def _load_gmrt(self):
        """GMRT - Global low-res.
        """
        
        this_gmrt = self.fetch_data('gmrt:layer=topo')        
        gmrt_result = this_gmrt.results[0]
        if gmrt_result[-1] == 0:
            gmrt_tif = gmrt_result[1]

            return(gmrt_tif)

        
    def fetch_atlxx(self, atl03_fn, short_name='ATL08'):
        """fetch an associated ATLxx file"""

        ## check for local file in atl03 dir and cache
        atlxx_filter = '_'.join(utils.fn_basename2(atl03_fn).split('_')[1:4])
        atlxx_filter_no_version = '_'.join(utils.fn_basename2(atl03_fn).split('_')[1:3])

        for this_dir in [os.path.dirname(atl03_fn), self.cache_dir]:
            tmp_bn = os.path.join(this_dir, f'{short_name}_{atlxx_filter}')
            atlxx_paths = glob.glob(f'{tmp_bn}*.h5')
            if len(atlxx_paths) > 0:
                return(atlxx_paths[0])
            else:
                tmp_bn = os.path.join(this_dir, f'{short_name}_{atlxx_filter_no_version}')
                atlxx_paths = glob.glob(f'{tmp_bn}*.h5')
                if len(atlxx_paths) > 0:
                    return(atlxx_paths[0])

        ## check with earthdata for aux files
        this_atlxx = earthdata.IceSat2(
            src_region=None,
            verbose=self.verbose,
            outdir=self.cache_dir,
            short_name=short_name,
            filename_filter=atlxx_filter,
            version='',
        )
        this_atlxx.run()
        if len(this_atlxx.results) == 0:
            this_atlxx = earthdata.IceSat2(
                src_region=None,
                verbose=self.verbose,
                outdir=self.cache_dir,
                short_name=short_name,
                filename_filter=atlxx_filter_no_version,
                version='',
            )
            this_atlxx.run()
            if len(this_atlxx.results) == 0:
                utils.echo_warning_msg(
                    (f'could not locate associated {short_name} '
                     f'file for {atlxx_filter}')
                )
                return(None)
        else:
            if this_atlxx.fetch_entry(
                    this_atlxx.results[0], check_size=True
            ) == 0:
                return(os.path.join(
                    this_atlxx._outdir, this_atlxx.results[0]['dst_fn']
                ))

        return(None)

    
    def apply_atl08_classifications(
            self, atl08_fn, ph_h_classed, laser, segment_id, segment_index_dict, ph_segment_ids
    ):
        with h5.File(atl08_fn, 'r') as atl08_f:
            if laser in atl08_f.keys():
                ## classed flag (signal_photons)
                atl08_classed_pc_flag  = atl08_f[
                    f'/{laser}/signal_photons/classed_pc_flag'
                ][...,]
                atl08_ph_segment_id = atl08_f[
                    f'/{laser}/signal_photons/ph_segment_id'
                ][...,] # photon src 20 m seg id
                atl08_classed_pc_indx = atl08_f[
                    f'/{laser}/signal_photons/classed_pc_indx'
                ][...,]

                ## set the classifications from ATL08
                atl08_segment_id_msk = [
                    True if x in segment_id else False for x in atl08_ph_segment_id
                ]
                atl08_ph_segment_indx = np.array(
                    list(
                        map((lambda pid: segment_index_dict[pid]),
                            atl08_ph_segment_id[atl08_segment_id_msk])
                    )
                )
                atl08_ph_index = np.array(
                    atl08_ph_segment_indx + (atl08_classed_pc_indx[atl08_segment_id_msk] - 1),
                    dtype=int
                )
                class_mask = atl08_ph_index < len(ph_segment_ids)
                ph_h_classed[atl08_ph_index[class_mask]] \
                    = atl08_classed_pc_flag[atl08_segment_id_msk][class_mask]
            
        return(ph_h_classed)


    def apply_atl24_classifications(
            self, atl24_fn, ph_h_classed, ph_h_bathy_conf, latitude, longitude, photon_h, photon_h_meantide,
            photon_h_geoid, laser, geoseg_beg, geoseg_end, ph_segment_ids
    ):
        with h5.File(atl24_fn, 'r') as atl24_f:
            ## some atl24 files don't have all the lasers, so make sure
            ## it exists before proceding
            if laser in atl24_f.keys():            
                lat_ph = atl24_f[f'{laser}/lat_ph'][...,]
                lon_ph = atl24_f[f'{laser}/lon_ph'][...,]
                ortho_h = atl24_f[f'{laser}/ortho_h'][...,]
                ellipse_h = atl24_f[f'{laser}/ellipse_h'][...,]
                surface_h = atl24_f[f'{laser}/surface_h'][...,]
                class_ph = atl24_f[f'{laser}/class_ph'][...,]
                index_seg = atl24_f[f'{laser}/index_seg'][...,]
                index_ph = atl24_f[f'{laser}/index_ph'][...,]
                conf_ph = atl24_f[f'{laser}/confidence'][...,]

                #if 'subset' in self.atl03_fn:
                ## determine the original `geosed_id` to determine the subset indices
                orig_segment_id = np.arange(geoseg_beg, geoseg_end+1, 1)
                orig_segment_id_msk = np.isin(orig_segment_id, ph_segment_ids)
                index_seg_orig = orig_segment_id[index_seg]
                segment_id_msk = np.isin(index_seg_orig, ph_segment_ids)
                if len(index_ph[segment_id_msk]) > 0:
                    ## check for previous index
                    index_ph = index_ph - np.min(index_ph[segment_id_msk])
                    index_msk = (index_ph[segment_id_msk] > 0) & (index_ph[segment_id_msk] < np.max(index_ph[segment_id_msk]))
                    class_msk = class_ph[segment_id_msk][index_msk] >= 40
                    ph_msk = (class_msk)
                    if self.min_bathy_confidence is not None:
                        confidence_msk = conf_ph[segment_id_msk][index_msk] >= self.min_bathy_confidence
                        ph_msk = (class_msk) & (confidence_msk)

                    ph_h_classed[index_ph[segment_id_msk][index_msk][ph_msk]] = class_ph[segment_id_msk][index_msk][ph_msk]
                    ph_h_bathy_conf[index_ph[segment_id_msk][index_msk][ph_msk]] = conf_ph[segment_id_msk][index_msk][ph_msk]
                    # we also need to change the lon/lat/height values to the
                    # updated/refracted bathymetry values (we'll just do it to class 40)
                    class_msk = class_ph[segment_id_msk][index_msk] == 40
                    ph_class = (class_msk)
                    if self.min_bathy_confidence is not None:
                        ph_msk = (class_msk) & (confidence_msk)

                    longitude[index_ph[segment_id_msk][index_msk][ph_msk]] = lon_ph[segment_id_msk][index_msk][ph_msk]
                    latitude[index_ph[segment_id_msk][index_msk][ph_msk]] = lat_ph[segment_id_msk][index_msk][ph_msk]
                    photon_h[index_ph[segment_id_msk][index_msk][ph_msk]] = ellipse_h[segment_id_msk][index_msk][ph_msk]
                    photon_h_meantide[index_ph[segment_id_msk][index_msk][ph_msk]] = surface_h[segment_id_msk][index_msk][ph_msk]
                    photon_h_geoid[index_ph[segment_id_msk][index_msk][ph_msk]] = ortho_h[segment_id_msk][index_msk][ph_msk]
                    #conf[index_ph[segment_id_msk][index_msk][ph_msk]] = conf_ph[segment_id_msk][index_msk][ph_msk]

        return(ph_h_classed, ph_h_bathy_conf, latitude, longitude, photon_h, photon_h_meantide, photon_h_geoid)


    def read_atl03(self, atl03_f, laser_num, orientation=None, atl08_fn=None, atl24_fn=None):
        """Read data from an ATL03 file

        Adapted from 'cshelph' https://github.com/nmt28/C-SHELPh.git 
        and 'iVert' https://github.com/ciresdem/ivert.git

        laser_num is 1, 2 or 3
        surface is 'mean_tide', 'geoid' or 'ellipsoid'
        """

        if orientation is None:
            orientation = atl03_f['/orbit_info/sc_orient'][0]
            
        ## selects the strong beams only [we can include weak beams later on]
        #orientDict = {0:'l', 1:'r', 21:'error'}
        laser = 'gt' + laser_num + self.orientDict[orientation]

        ## for 'subsets', where heights don't come through
        if 'heights' not in atl03_f['/{}'.format(laser)].keys():
            return(None)
        
        ## Read in the required photon level data
        photon_h = atl03_f['/' + laser + '/heights/h_ph'][...,]
        latitude = atl03_f['/' + laser + '/heights/lat_ph'][...,]
        longitude = atl03_f['/' + laser + '/heights/lon_ph'][...,]
        ph_count = atl03_f['/' + laser + '/heights/ph_id_count'][...,]
        conf = atl03_f['/' + laser + '/heights/signal_conf_ph/'][...,0]
        qual = atl03_f['/' + laser + '/heights/quality_ph/'][...,0]
        dist_ph_along = atl03_f['/' + laser + '/heights/dist_ph_along'][...,]
        delta_time = atl03_f['/' + laser + '/heights/delta_time'][...,]
        this_N = latitude.shape[0]

        ## Read in the geolocation level data
        segment_ph_cnt = atl03_f['/' + laser + '/geolocation/segment_ph_cnt'][...,]
        segment_id = atl03_f['/' + laser + '/geolocation/segment_id'][...,]
        segment_dist_x = atl03_f['/' + laser + '/geolocation/segment_dist_x'][...,]
        ref_elev = atl03_f['/' + laser + '/geolocation/ref_elev'][...,]
        ref_azimuth = atl03_f['/' + laser + '/geolocation/ref_azimuth'][...,]
        ref_photon_index = atl03_f['/' + laser + '/geolocation/reference_photon_index'][...,]
        ph_index_beg = atl03_f['/' + laser + '/geolocation/ph_index_beg'][...,]
        altitude_sc = atl03_f['/' + laser + '/geolocation/altitude_sc'][...,]
        surf_type = atl03_f['/' + laser + '/geolocation/surf_type'][...,]
        
        ## Read in the geophysical correction
        photon_geoid = atl03_f['/' + laser + '/geophys_corr/geoid'][...,]
        photon_geoid_f2m = atl03_f['/' + laser + '/geophys_corr/geoid_free2mean'][...,]
        photon_dem_h = atl03_f['/' + laser + '/geophys_corr/dem_h'][...,]

        ## Read in the geosegment info for subsets
        geoseg_beg = atl03_f['ancillary_data/start_geoseg'][...,][0]
        geoseg_end = atl03_f['ancillary_data/end_geoseg'][...,][0]
                    
        ## Create a dictionary with (segment_id --> index into ATL03 photons)
        ## lookup pairs, for the starting photon of each segment
        #segment_indices = ph_index_beg+(ref_photon_index-1)
        #segment_index_dict = dict(zip(segment_id, ph_index_beg+ref_photon_index-1))
        segment_indices = np.concatenate(([0], np.cumsum(segment_ph_cnt)[:-1]))
        segment_index_dict = dict(zip(segment_id, segment_indices))
        ph_segment_ids = segment_id[
            np.searchsorted(segment_indices, np.arange(0.5, len(photon_h), 1))-1
        ]
        
        ## Compute the total along-track distances.
        segment_dist_dict = dict(zip(segment_id, segment_dist_x))

        ## Determine where in the array each segment index needs to look.
        ph_segment_dist_x = np.array(
            list(map((lambda pid: segment_dist_dict[pid]), ph_segment_ids))
        )
        dist_x = ph_segment_dist_x + dist_ph_along

        ## meantide/geoid/dem heights
        h_geoid_dict = dict(zip(segment_id, photon_geoid))
        ph_h_geoid = np.array(list(map((lambda pid: h_geoid_dict[pid]), ph_segment_ids)))        
        h_meantide_dict = dict(zip(segment_id, photon_geoid_f2m))
        ph_h_meantide = np.array(list(map((lambda pid: h_meantide_dict[pid]), ph_segment_ids)))
        photon_h_geoid = photon_h - ph_h_geoid
        photon_h_meantide = photon_h - (ph_h_geoid + ph_h_meantide)

        h_dem_dict = dict(zip(segment_id, photon_dem_h))
        ph_h_dem = np.array(list(map((lambda pid: h_dem_dict[pid]), ph_segment_ids)))
        photon_h_dem = ph_h_dem - (ph_h_geoid + ph_h_meantide)
        
        ## setup classifications asn bathy confidence
        ph_h_classed = np.zeros(photon_h.shape)
        ph_h_classed[:] = -1

        ph_h_bathy_conf = np.zeros(photon_h.shape)
        ph_h_bathy_conf[:] = -1

        ## append the laser to each record
        laser_arr = np.empty(photon_h.shape, dtype='object')
        laser_arr[:] = laser

        ## append the filename to each record
        fn_arr = np.empty(photon_h.shape, dtype='object')
        fn_arr[:] = self.fn
        
        ## ref values
        h_ref_elev_dict = dict(zip(segment_id, ref_elev))
        ph_ref_elev = np.array(
            list(map((lambda pid: h_ref_elev_dict[pid]), ph_segment_ids))
        )#.astype(float)        
        h_ref_azimuth_dict = dict(zip(segment_id, ref_azimuth))
        ph_ref_azimuth = np.array(
            list(map((lambda pid: h_ref_azimuth_dict[pid]), ph_segment_ids))
        )#.astype(float)
        h_altitude_sc_dict = dict(zip(segment_id, altitude_sc))
        ph_altitude_sc = np.array(
            list(map((lambda pid: h_altitude_sc_dict[pid]), ph_segment_ids))
        )#.astype(float)
        h_surf_type_dict = dict(zip(segment_id, surf_type))
        ph_surf_type = np.array(
            list(map((lambda pid: h_surf_type_dict[pid]), ph_segment_ids))
        )#.astype(float)
        
        ## classify the photons
        ## Read in the ATL08 data and classify points based on the ATL08 classifications
        if atl08_fn is not None:
            ph_h_classed = self.apply_atl08_classifications(
                atl08_fn, ph_h_classed, laser, segment_id, segment_index_dict, ph_segment_ids
            )
        
        ## classify phtons as water surface if dem_h is below zero
        #ph_h_classed[photon_h_dem < 0] = 41
        
        ## classify bathymetry and water-surface points based on ATL24
        ## classifications...
        if atl24_fn is not None:
            ph_h_classed, ph_h_bathy_conf, latitude, longitude, photon_h, photon_h_meantide, photon_h_geoid = self.apply_atl24_classifications(
                atl24_fn, ph_h_classed, ph_h_bathy_conf, latitude, longitude, photon_h, photon_h_meantide,
                photon_h_geoid, laser, geoseg_beg, geoseg_end, ph_segment_ids
            )

        ## classify inland water based on ATL13 classifications
        # if self.atl13_f is not None:
        #     atl13_refid = self.atl13_f['/' + laser + '/segment_id_beg'][...,]
        #     ph_h_classed[atl13_refid] = 44

        ## classify water based on ATL12 classifications
        # if self.atl13_f is not None:
        #     atl13_refid = self.atl13_f['/' + laser + '/segment_id_beg'][...,]
        #     ph_h_classed[atl13_refid] = 44
        
        ## classify water points based on the 'surf_type' parameter to
        ## remove ATL08 bare-earth clasification over water
        ph_surf_type = [x[1] for x in ph_surf_type]
            
        ## set the photon height, either 'mean_tide' or 'geoid', else ellipsoid
        if self.water_surface == 'mean_tide':
            ph_height = photon_h_meantide
        elif self.water_surface == 'geoid':
            ph_height = photon_h_geoid
        else:
            ph_height = photon_h
            
        ## create the pandas dataframe            
        dataset = pd.DataFrame(
            {'latitude': latitude,
             'longitude': longitude,
             'photon_height': ph_height,
             'laser': laser_arr,
             'fn': fn_arr,
             'confidence': conf,
             'ref_elevation':ph_ref_elev,
             'ref_azimuth':ph_ref_azimuth,
             'ref_sat_alt':ph_altitude_sc,
             'delta_time':delta_time,
             'bathy_confidence':ph_h_bathy_conf,
             'ph_surf_type':ph_surf_type,
             'photon_h_dem':photon_h_dem,
             'ph_h_classed': ph_h_classed},
            columns=['latitude', 'longitude', 'photon_height', 'laser', 'fn',
                     'confidence', 'ref_elevation', 'ref_azimuth', 'ref_sat_alt',
                     'delta_time', 'bathy_confidence', 'ph_surf_type', 'photon_h_dem', 'ph_h_classed']
        )

        ## Process extra columns specified in `self.columns`
        for col in self.columns.keys():
            try:
                if 'gtx' in col:
                    _col = col.replace('/gtx', '/' + laser)
                    col_arr = self.atl03_f[_col][...,]                
                    if 'heights' not in _col: 
                        col_dict = dict(zip(segment_id, col_arr))
                        col_arr = np.array(
                            list(map((lambda pid: col_dict[pid]), ph_segment_ids))
                        )
                else:
                    col_arr = np.empty(photon_h.shape, dtype='object')
                    col_arr[:] = self.atl03_f[col][...,]

                extra_dataset = pd.DataFrame(
                    {self.columns[col]: col_arr},
                    columns = [self.columns[col]]
                )
                dataset = dataset.join(extra_dataset)
            except:
                utils.echo_warning_msg(
                    f'could not find and/or process {col}'
                )

        return(dataset)        


    def yield_points(self):
        dataset = None

        ## only fetch atl08/atl24 if classes are desired
        atl08_fn = None
        atl24_fn = None
        if len(self.classes) > 0:
            atl08_fn = self.fetch_atlxx(self.fn, short_name='ATL08')
            atl24_fn = self.fetch_atlxx(self.fn, short_name='ATL24')

        ## fetch and process buildings, if wanted
        if self.classify_buildings:
            if isinstance(self.classify_buildings, bool):
                this_bing = self.process_buildings(
                    self.fetch_buildings(verbose=self.verbose),
                    verbose=self.verbose
                )
            elif isinstance(self.classify_buildings, list):
                this_bing = self.classify_buildings
            elif isinstance(self.classify_buildings, str):
                if os.path.exists(self.classify_buildings):
                    this_bing = gdalfun.ogr_geoms(self.classify_buildings)

        ## fetch and process watermask, if wanted
        if self.classify_water:
            if isinstance(self.classify_water, bool):
                this_wm = self.process_coastline(
                    self.fetch_coastline(chunks=False, verbose=self.verbose),
                    return_geom=True,
                    verbose=self.verbose
                )
            elif isinstance(self.classify_water, list):
                this_wm = self.classify_water
            elif isinstance(self.classify_water, str):
                if os.path.exists(self.classify_water):
                    this_wm = gdalfun.ogr_geoms(self.classify_water)

        ## parse through the icesat2 file by laser number
        with h5.File(self.fn, 'r') as atl03_f:
            if 'short_name' in atl03_f.attrs.keys():
                if self.reject_failed_qa:
                    if atl03_f[
                            '/quality_assessment/qa_granule_pass_fail'
                    ][...,][0] == 0:
            
                        for i in range(1, 4):
                            for orient in range(2):
                                ## selects the strong beams only [we can include weak beams later on]
                                if f'gt{i}{self.orientDict[orient]}' not in list(atl03_f.keys()):
                                    continue

                                if 'heights' not in list(atl03_f[f'gt{i}{self.orientDict[orient]}'].keys()):
                                    continue

                                dataset = self.read_atl03(atl03_f, '{}'.format(i), orientation=orient, atl08_fn=atl08_fn, atl24_fn=atl24_fn)
                                if dataset is None or len(dataset) == 0:
                                    continue

                                ## keep only photons with confidence levels mentioned
                                ## in `self.confidence_levels`
                                if len(self.confidence_levels) > 0:
                                    dataset = dataset[
                                        (np.isin(dataset['confidence'], self.confidence_levels))
                                    ]

                                ## reduce the dataset to the input region for faster masking
                                dataset = dataset[(dataset['longitude'] >= self.region.xmin) & (dataset['longitude'] <= self.region.xmax)]
                                dataset = dataset[(dataset['latitude'] >= self.region.ymin) & (dataset['latitude'] <= self.region.ymax)]

                                if dataset is None or len(dataset) == 0:
                                    continue

                                ## re-classify photons based on buildings/watermask/bathymetry
                                if self.classify_buildings and this_bing is not None:
                                    #dataset = self.classify_buildings(dataset, this_bing)
                                    dataset = self.classify_by_mask_geoms(dataset, mask_geoms=this_bing, classification=7)

                                #utils.echo_msg(f'this_wm: {this_wm}')
                                if self.classify_water and this_wm is not None:
                                    #dataset = self.classify_water(dataset, this_wm)
                                    dataset = self.classify_by_mask_geoms(dataset, mask_geoms=this_wm, classification=41, except_classes=[40])

                                if dataset is None or len(dataset) == 0:
                                    continue

                                # ## bathymetry is classified in `read_atl_data` using ATL24 now...
                                # if self.want_bathymetry:
                                #     dataset = self.classify_bathymetry(dataset)

                                # if dataset is None or len(dataset) == 0:
                                #     continue

                                ## keep only photons with a classification mentioned in `self.classes`
                                if len(self.classes) > 0:
                                    dataset = dataset[
                                        (np.isin(dataset['ph_h_classed'], self.classes))
                                    ]

                                if dataset is None or len(dataset) == 0:
                                    continue

                                ## rename the x,y,z columns for `transform_and_yield_points`
                                dataset.rename(
                                    columns={
                                        'longitude': 'x', 'latitude': 'y', 'photon_height': 'z'
                                    },
                                    inplace=True
                                )
                                yield(dataset)

        
    def _vectorize_df(self, dataset):
        """Make a point vector OGR DataSet Object from a pandas dataframe

        This is to allow spatial filtering for watermask and buildings.
        """

        dst_ogr = '{}'.format('icesat_dataframe')
        ogr_ds = gdal.GetDriverByName('Memory').Create(
            '', 0, 0, 0, gdal.GDT_Unknown
        )
        layer = ogr_ds.CreateLayer(
            dst_ogr, geom_type=ogr.wkbPoint
        )
        fd = ogr.FieldDefn('index', ogr.OFTInteger)
        layer.CreateField(fd)
        f = ogr.Feature(feature_def=layer.GetLayerDefn())        
        with tqdm(desc='vectorizing dataframe', leave=False) as pbar:
            for index, this_row in dataset.iterrows():
                pbar.update()
                f.SetField(0, index)
                g = ogr.CreateGeometryFromWkt(
                    f'POINT ({this_row["longitude"]} {this_row["latitude"]})'
                )
                f.SetGeometryDirectly(g)
                layer.CreateFeature(f)
            
        return(ogr_ds)


    def classify_by_mask_geoms(self, dataset, mask_geoms=[], classification=-1, except_classes=[]):
        """classify water photons using OSM coastline 
        """

        ## vectorize the photons
        ogr_df = self._vectorize_df(dataset)
        
        os.environ["OGR_OSM_OPTIONS"] = "INTERLEAVED_READING=YES"
        os.environ["OGR_OSM_OPTIONS"] = "OGR_INTERLEAVED_READING=YES"
        with tqdm(
                total=len(mask_geoms),
                desc='classifying water photons',
                leave=False
        ) as pbar:
            for n, wm_geom in enumerate(mask_geoms):                
                pbar.update()
                icesat_layer = ogr_df.GetLayer()
                if isinstance(wm_geom,  str):
                    wm_geom = ogr.CreateGeometryFromWkt(wm_geom)
                    
                icesat_layer.SetSpatialFilter(wm_geom)
                for f in icesat_layer:
                    idx = f.GetField('index')
                    #except_mask = dataset.at[idx, 'ph_h_classed']
                    if dataset.at[idx, 'ph_h_classed'] not in except_classes:
                        dataset.at[idx, 'ph_h_classed'] = classification

                icesat_layer.SetSpatialFilter(None)

        ogr_df = None
        return(dataset)

    
### End
