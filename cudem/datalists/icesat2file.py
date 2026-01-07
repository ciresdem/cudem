### icesat2file.py - DataLists IMproved
##
## Copyright (c) 2023 - 2026 Regents of the University of Colorado
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
### Commentary:
##
## ICESat-2 Data Parser (ATL03, ATL24)
##
## Product,Type,Useful For,Class
## ATL08,Land/Canopy,"Ground (1), Canopy (2), Top Canopy (3)","1, 2, 3"
## ATL24,Bathymetry,"Near-shore Bathymetry (40), Coastline (41)","40, 41"
## ATL06,Land Ice,"Glaciers, Ice Sheets",6
## ATL13,Inland Water,"Lakes, Rivers, Reservoirs",42
## ATL12,Ocean,Open Ocean Surface,44
## ATL09,Atmosphere,Cloud Layers (Filtering Noise),0 (Noise) (todo)
##
### Code:

import os
import glob
import numpy as np
import h5py as h5
import pandas as pd
from tqdm import tqdm

from osgeo import gdal, ogr

from cudem import utils
from cudem import regions
from cudem import gdalfun
from cudem.fetches import earthdata
from cudem.fetches import osm
from cudem.fetches import bingbfp
from cudem.datalists.dlim import ElevationDataset


## ==============================================
## ATL-24 Dataset (Bathymetry)
## ==============================================
class IceSat2_ATL24(ElevationDataset):
    """ICESat-2 ATL24 (Bathymetry) Data Parser.
    """
    
    def __init__(self, min_confidence=None, classes='40', water_surface='ortho', **kwargs):
        super().__init__(**kwargs)
        self.orientDict = {0:'l', 1:'r', 21:'error'}
        self.lasers = ['gt1l', 'gt2l', 'gt3l', 'gt1r', 'gt2r', 'gt3r']

        self.data_format = 304
        self.water_surface = water_surface if water_surface in ['surface', 'ortho', 'ellipse'] else 'ortho'
        self.min_confidence = utils.float_or(min_confidence)
        
        self.classes = []
        if classes is not None:
            self.classes = [int(x) for x in str(classes).split('/')]

            
    def generate_inf(self, make_grid=False, make_block_mean=False, block_inc=None):
        """Generate metadata (INF) for the ATL24 data.
        """
        
        ## If grids are requested, delegate to parent to perform full scan
        if make_grid or make_block_mean:
            return super().generate_inf(make_grid, make_block_mean, block_inc)

        if self.src_srs is None:
            self.infos.src_srs = 'epsg:4326+3855'
        else:
            self.infos.src_srs = self.src_srs

        these_regions = []
        numpts = 0
        
        try:
            with h5.File(self.fn, 'r') as f:
                for b in range(1, 4):
                    for p in ['l', 'r']:
                        beam = f'gt{b}{p}'
                        if beam not in f: continue
                        
                        try:
                            ## Read bounds from HDF5 datasets directly for speed
                            y = f[f'{beam}/lat_ph'][...]
                            x = f[f'{beam}/lon_ph'][...]
                            z = f[f'{beam}/ortho_h'][...]
                            
                            if len(x) > 0:
                                this_region = regions.Region(src_srs=self.infos.src_srs).from_list(
                                    [np.min(x), np.max(x), np.min(y), np.max(y), np.min(z), np.max(z)]
                                )
                                these_regions.append(this_region)
                                numpts += len(x)
                        except KeyError:
                            continue

            ## Merge regions
            if these_regions:
                merged_region = these_regions[0].copy()
                for r in these_regions[1:]:
                    merged_region = regions.regions_merge(merged_region, r)

                self.infos.minmax = merged_region.export_as_list(include_z=True)
                self.infos.wkt = merged_region.export_as_wkt()
                
            self.infos.numpts = numpts 
            
        except Exception as e:
            utils.echo_warning_msg(f"Error generating INF for {self.fn}: {e}")
            ## Fallback to scan if HDF5 read fails
            return super().generate_inf(make_grid, make_block_mean, block_inc)

        return self.infos

    
    def yield_points(self):
        """Yield points from ATL24 HDF5 file."""
        
        with h5.File(self.fn, 'r') as f:
            for b in range(1, 4):
                for p in ['l', 'r']:
                    beam = f'gt{b}{p}'
                    if beam not in f: continue

                    try:
                        # Read Arrays
                        lat_ph = f[f'{beam}/lat_ph'][...]
                        lon_ph = f[f'{beam}/lon_ph'][...]
                        class_ph = f[f'{beam}/class_ph'][...]
                        conf_ph = f[f'{beam}/confidence'][...]
                        
                        ## Select Height Type
                        if self.water_surface == 'surface':
                            ph_height = f[f'{beam}/surface_h'][...]
                        elif self.water_surface == 'ellipse':
                            ph_height = f[f'{beam}/ellipse_h'][...]
                        else:
                            ph_height = f[f'{beam}/ortho_h'][...]

                        ## Metadata columns
                        laser_arr = np.full(ph_height.shape, beam, dtype='object')
                        fn_arr = np.full(ph_height.shape, self.fn, dtype='object')
                        
                        dataset = pd.DataFrame({
                            'latitude': lat_ph,
                            'longitude': lon_ph,
                            'photon_height': ph_height,
                            'laser': laser_arr,
                            'fn': fn_arr,
                            'confidence': conf_ph,
                            'ph_h_classed': class_ph
                        })

                        ## Filter by Class
                        if self.classes:
                            dataset = dataset[dataset['ph_h_classed'].isin(self.classes)]

                        ## Filter by Confidence
                        if self.min_confidence is not None:
                            dataset = dataset[dataset['confidence'] >= self.min_confidence]
                            
                        if dataset.empty: continue
                        
                        ## Normalize columns for CUDEM
                        dataset.rename(columns={
                            'longitude': 'x', 'latitude': 'y', 'photon_height': 'z'
                        }, inplace=True)
                        
                        yield dataset

                    except KeyError as e:
                        if self.verbose:
                            utils.echo_warning_msg(f"Missing dataset in {beam}: {e}")
                        continue


## ==============================================
## ATL03 Dataset (full) - classified
## ==============================================
def apply_atl06_classifications(df, atl06_fn, laser, verbose=True):
    """Map ATL06 Land Ice segments to ATL03 photons (Class 6)."""

    try:
        with h5.File(atl06_fn, 'r') as f:
            if laser not in f: return df

            li = f[f'/{laser}/land_ice_segments']
            if 'segment_id' not in li: return df

            atl06_seg = li['segment_id'][...]
            atl06_h = li['h_li'][...]

            valid_mask = atl06_h < 1e30
            valid_segs = atl06_seg[valid_mask]

            ## Identify photons belonging to these valid segments
            is_ice = df['ph_segment_id'].isin(valid_segs)

            ## Class 6 (Ice), override unclassified/ground/canopy/noise (<7)
            mask = is_ice & (df['ph_h_classed'] < 7)
            if verbose:
                utils.echo_msg(f'Classified {np.count_nonzero(mask)} photons from {laser} based on ATL06')

            df.loc[mask, 'ph_h_classed'] = 6

    except Exception as e:
        utils.echo_warning_msg(f"Failed to apply ATL06 classifications from {atl06_fn}: {e}")

    return df

    
def apply_atl08_classifications(df, atl08_fn, laser, segment_index_dict, verbose=True):
    """Map ATL08 Land/Canopy classes to ATL03 photons."""

    try:
        with h5.File(atl08_fn, 'r') as f:
            if laser not in f: return df

            sig = f[f'/{laser}/signal_photons']
            atl08_flag = sig['classed_pc_flag'][...]
            atl08_seg = sig['ph_segment_id'][...]
            atl08_idx = sig['classed_pc_indx'][...] # 1-based index within segment

            ## Filter to segments present in our current dataframe
            ## We use the unique segments from the DF to avoid processing unrelated ATL08 data
            relevant_segments = df['ph_segment_id'].unique()
            mask = np.isin(atl08_seg, relevant_segments)

            if not np.any(mask): return df

            ## Calculate Absolute Index in ATL03 Arrays
            ## Global_Index = Segment_Start_Index + (Photon_Index_In_Segment - 1)

            ## Get start indices for the relevant segments
            ## Map segment IDs to their start index in the DF (using the passed dict)
            seg_starts = np.array([segment_index_dict.get(s, -1) for s in atl08_seg[mask]])

            ## Calculate absolute indices
            valid_seg_mask = seg_starts != -1
            atl03_indices = seg_starts[valid_seg_mask] + (atl08_idx[mask][valid_seg_mask] - 1)

            ## Filter indices that might be out of bounds (truncated granule)
            valid_idx_mask = (atl03_indices >= 0) & (atl03_indices < len(df))
            final_indices = atl03_indices[valid_idx_mask]

            ## Apply Classification
            values_to_assign = atl08_flag[mask][valid_seg_mask][valid_idx_mask]

            if verbose:
                utils.echo_msg(f'Classified {len(final_indices)} photons from {laser} based on ATL08')
            ## Using numpy indexing on the underlying column for speed, or .iloc
            #df.iloc[final_indices, df.columns.get_loc('ph_h_classed')] = values_to_assign
            ## A safer pandas way:
            df.loc[df.index[final_indices], 'ph_h_classed'] = values_to_assign

    except Exception as e:
        utils.echo_warning_msg(f"Failed to apply ATL08 classifications from {atl08_fn}: {e}")

    return df


def apply_atl12_classifications(df, atl12_fn, laser, verbose=True):
    """Map ATL12 Ocean Surface segments to ATL03 photons (Class 44)."""

    try:
        with h5.File(atl12_fn, 'r') as f:
            if laser not in f: return df

            path = f'/{laser}/ssh_segments/stats'
            if path not in f or 'segment_id_beg' not in f[path]: return df

            atl12_seg = f[f'{path}/segment_id_beg'][...]

            is_ocean = df['ph_segment_id'].isin(atl12_seg)

            ## Prioritize Bathy (40) and Inland (42) over Ocean (44)
            mask = is_ocean & (df['ph_h_classed'] < 40)
            if verbose:
                utils.echo_msg(f'Classified {np.count_nonzero(final_indices)} photons from {laser} based on ATL12')

            df.loc[mask, 'ph_h_classed'] = 44

    except Exception as e:
        utils.echo_warning_msg(f"Failed to apply ATL12 classifications from {atl12_fn}: {e}")

    return df


def apply_atl13_classifications(df, atl13_fn, laser, verbose=True):
    """Map ATL13 Inland Water segments to ATL03 photons (Class 42)."""

    try:
        with h5.File(atl13_fn, 'r') as f:
            if laser not in f: return df

            if f'/{laser}/segment_id_beg' not in f: return df

            atl13_seg = f[f'/{laser}/segment_id_beg'][...]

            is_water = df['ph_segment_id'].isin(atl13_seg)

            ## Prioritize Bathy (40) over inland water
            mask = is_water & (df['ph_h_classed'] < 40)
            utils.echo_msg(f'Classed {np.count_nonzero(mask)} photons from {laser} based on ATL13')
            df.loc[mask, 'ph_h_classed'] = 42

    except Exception as e:
        utils.echo_warning_msg(f"Failed to apply ATL13 classifications from {atl13_fn}: {e}")

    return df

    
class IceSat2_ATL03(ElevationDataset):
    """ICESat-2 ATL03 (Global Geolocated Photon Data) Parser.
    
    Can integrate classifications from ATL08 (Land/Canopy), ATL24 (Bathymetry),
    and ATL06 (Land Ice), and ATL13 (Inland Water). 
    Can also classify using external vectors (OSM/Bing) or Algorithms (Known Bathy/DBSCAN).

    Ground (1), Canopy (2), Top Canopy (3), Bathymetry (40), Water Surface (41), Land Ice (6), Inland Water (42)
    Open Ocean Surface (44), Atmosphere (0), Buildings (7), Unclassified (-1)

    """

    def __init__(self,
                 water_surface='geoid', # height
                 classes=None, # classify the data
                 confidence_levels='2/3/4', # filter data by confidence level
                 columns=None,
                 reject_failed_qa=True, # reject bad granules
                 classify_buildings=True, # classify buildings
                 classify_water=True, # classify water surface 
                 classify_inland_water=True, # classify inland water surface
                 append_atl24=False, # append atl24 data to yield_points
                 min_bathy_confidence=None, # minimum desired bathy confidence (from atl24)
                 known_bathymetry=None,  # New: Path to raster or 'gmrt'
                 known_bathy_threshold=5.0, # New: Vertical window (+/- meters)
                 use_dbscan=False,       # New: Enable DBSCAN clustering
                 dbscan_eps=1.5,         # New: DBSCAN epsilon (neighbor distance)
                 dbscan_min_samples=10,  # New: DBSCAN min points for core
                 **kwargs):
        
        super().__init__(**kwargs)
        
        self.data_format = 303
        self.water_surface = water_surface if water_surface in ['mean_tide', 'geoid', 'ellipsoid'] else 'mean_tide'

        self.classes = [int(x) for x in str(classes).split('/')] if classes is not None else []
        self.confidence_levels = [int(x) for x in str(confidence_levels).split('/')] if confidence_levels is not None else []
        self.columns = columns if columns is not None else {}
        
        self.classify_buildings = classify_buildings
        self.classify_water = classify_water
        self.classify_inland_water = classify_inland_water
        self.reject_failed_qa = reject_failed_qa
        self.append_atl24 = append_atl24
        self.min_bathy_confidence = utils.float_or(min_bathy_confidence)

        ## --- SciKit Algo Classification Options ---
        self.known_bathymetry = known_bathymetry
        self.known_bathy_threshold = utils.float_or(known_bathy_threshold, 5.0)
        self.use_dbscan = use_dbscan
        self.dbscan_eps = utils.float_or(dbscan_eps, 1.5)
        self.dbscan_min_samples = utils.int_or(dbscan_min_samples, 10)
        
        self.orientDict = {0:'l', 1:'r', 21:'error'}


    def _init_region_from_df(self, df):
        """Initialize a region object from the dataframe extents."""
        
        if df is None or df.empty:
            return None
            
        return regions.Region().from_list([
            df['longitude'].min(), df['longitude'].max(),
            df['latitude'].min(), df['latitude'].max()
        ])
    

    def _get_processing_region(self):
        """Determine the processing region for fetches."""
        
        if self.region is not None:
            return self.region.copy()

        #utils.echo_msg(self.infos.minmax)
        #if isinstance(self.infos, dict) and 'minmax' in self.infos and self.infos.minmax is not None:
        if self.infos.minmax is not None:
            # minmax format: [xmin, xmax, ymin, ymax, zmin, zmax]
            return regions.Region().from_list(self.infos.minmax)
            
        return None
    
    
    def generate_inf(self, make_grid=False, make_block_mean=False, block_inc=None):
        """Generate metadata (INF) for the ATL03 data."""

        ## Delegate to parent for grids
        make_grid = False
        make_block_mean = False
        # if make_grid or make_block_mean:
        #     return super().generate_inf(make_grid, make_block_mean, block_inc)

        if self.src_srs is None:
            self.infos.src_srs = 'epsg:4326+3855'
        else:
            self.infos.src_srs = self.src_srs

        these_regions = []
        numpts = 0
        try:
            #utils.echo_msg(self.fn)
            with h5.File(self.fn, 'r') as f:
                for b in range(1, 4):
                    for p in ['l', 'r']:
                        beam = f'gt{b}{p}'
                        if beam not in f or 'heights' not in f[beam]: continue

                        y = f[f'{beam}/heights/lat_ph'][...]
                        x = f[f'{beam}/heights/lon_ph'][...]
                        z = f[f'{beam}/heights/h_ph'][...]
                        if len(x) > 0:
                            this_region = regions.Region(src_srs=self.infos.src_srs).from_list(
                                [np.min(x), np.max(x), np.min(y), np.max(y), np.min(z), np.max(z)]
                            )
                            these_regions.append(this_region)
                            numpts += len(x)

            if these_regions:
                merged = these_regions[0].copy()
                for r in these_regions[1:]:
                    merged = regions.regions_merge(merged, r)

                self.infos.minmax = merged.export_as_list(include_z=True)
                self.infos.wkt = merged.export_as_wkt()

            self.infos.numpts = numpts 
            
        except Exception as e:
            ## Fallback to scan
            utils.echo_warning_msg(f'Failed to generate inf from h5, falling back to parsing all points, {e}')
            return super().generate_inf(make_grid, make_block_mean, block_inc)

        return self.infos

    
    def fetch_atlxx(self, atl03_fn, short_name='ATL08'):
        """Fetch associated ATLxx file (ATL08, ATL24, ATL06, etc)."""
        
        ## Build filename patterns
        bn = utils.fn_basename2(atl03_fn)
        parts = bn.split('_')
        if len(parts) < 4: return None
        
        atlxx_filter = '_'.join(parts[1:4])
        atlxx_filter_no_ver = '_'.join(parts[1:3])

        ## Check Local/Cache
        for d in [os.path.dirname(atl03_fn), self.cache_dir]:
            for filt in [atlxx_filter, atlxx_filter_no_ver]:
                matches = glob.glob(os.path.join(d, f'{short_name}_{filt}*.h5'))
                if matches: return matches[0]

        ## Fetch from Earthdata
        for filt in [atlxx_filter, atlxx_filter_no_ver]:
            fetcher = earthdata.IceSat2(
                src_region=None, verbose=self.verbose, outdir=self.cache_dir,
                short_name=short_name, filename_filter=filt, version=''
            )
            fetcher.run()
            if fetcher.results:
                if fetcher.fetch_entry(fetcher.results[0], check_size=True) == 0:
                    return os.path.join(fetcher._outdir, fetcher.results[0]['dst_fn'])
        
        return None

    
    def apply_atl08_classifications(self, df, atl08_fn, laser, segment_index_dict):
        """Map ATL08 Land/Canopy classes to ATL03 photons."""
        
        try:
            with h5.File(atl08_fn, 'r') as f:
                if laser not in f: return df
                
                sig = f[f'/{laser}/signal_photons']
                atl08_flag = sig['classed_pc_flag'][...]
                atl08_seg = sig['ph_segment_id'][...]
                atl08_idx = sig['classed_pc_indx'][...] # 1-based index within segment

                ## Filter to segments present in our current dataframe
                ## We use the unique segments from the DF to avoid processing unrelated ATL08 data
                relevant_segments = df['ph_segment_id'].unique()
                mask = np.isin(atl08_seg, relevant_segments)
                
                if not np.any(mask): return df

                ## Calculate Absolute Index in ATL03 Arrays
                ## Global_Index = Segment_Start_Index + (Photon_Index_In_Segment - 1)
                
                ## Get start indices for the relevant segments
                ## Map segment IDs to their start index in the DF (using the passed dict)
                seg_starts = np.array([segment_index_dict.get(s, -1) for s in atl08_seg[mask]])
                
                ## Calculate absolute indices
                valid_seg_mask = seg_starts != -1
                atl03_indices = seg_starts[valid_seg_mask] + (atl08_idx[mask][valid_seg_mask] - 1)
                
                ## Filter indices that might be out of bounds (truncated granule)
                valid_idx_mask = (atl03_indices >= 0) & (atl03_indices < len(df))
                final_indices = atl03_indices[valid_idx_mask]
                
                ## Apply Classification
                values_to_assign = atl08_flag[mask][valid_seg_mask][valid_idx_mask]

                if self.verbose:
                    utils.echo_msg(f'Classified {len(final_indices)} photons from {laser} based on ATL08')
                ## Using numpy indexing on the underlying column for speed, or .iloc
                #df.iloc[final_indices, df.columns.get_loc('ph_h_classed')] = values_to_assign
                ## A safer pandas way:
                df.loc[df.index[final_indices], 'ph_h_classed'] = values_to_assign

        except Exception as e:
            utils.echo_warning_msg(f"Failed to apply ATL08 classifications from {atl08_fn}: {e}")
            
        return df

    
    def apply_atl06_classifications(self, df, atl06_fn, laser):
        """Map ATL06 Land Ice segments to ATL03 photons (Class 6)."""
        
        try:
            with h5.File(atl06_fn, 'r') as f:
                if laser not in f: return df
                
                li = f[f'/{laser}/land_ice_segments']
                if 'segment_id' not in li: return df

                atl06_seg = li['segment_id'][...]
                atl06_h = li['h_li'][...]
                
                valid_mask = atl06_h < 1e30
                valid_segs = atl06_seg[valid_mask]

                ## Identify photons belonging to these valid segments
                is_ice = df['ph_segment_id'].isin(valid_segs)
                
                ## Class 6 (Ice), override unclassified/ground/canopy/noise (<7)
                mask = is_ice & (df['ph_h_classed'] < 7)
                if self.verbose:
                    utils.echo_msg(f'Classified {np.count_nonzero(mask)} photons from {laser} based on ATL06')
                    
                df.loc[mask, 'ph_h_classed'] = 6

        except Exception as e:
            utils.echo_warning_msg(f"Failed to apply ATL06 classifications from {atl06_fn}: {e}")
            
        return df

    
    def apply_atl12_classifications(self, df, atl12_fn, laser):
        """Map ATL12 Ocean Surface segments to ATL03 photons (Class 44)."""
        
        try:
            with h5.File(atl12_fn, 'r') as f:
                if laser not in f: return df
                
                path = f'/{laser}/ssh_segments/stats'
                if path not in f or 'segment_id_beg' not in f[path]: return df
                
                atl12_seg = f[f'{path}/segment_id_beg'][...]
                
                is_ocean = df['ph_segment_id'].isin(atl12_seg)
                
                ## Prioritize Bathy (40) and Inland (42) over Ocean (44)
                mask = is_ocean & (df['ph_h_classed'] < 40)
                if self.verbose:
                    utils.echo_msg(f'Classified {np.count_nonzero(final_indices)} photons from {laser} based on ATL12')
                    
                df.loc[mask, 'ph_h_classed'] = 44

        except Exception as e:
            utils.echo_warning_msg(f"Failed to apply ATL12 classifications from {atl12_fn}: {e}")
            
        return df

    
    def apply_atl13_classifications(self, df, atl13_fn, laser):
        """Map ATL13 Inland Water segments to ATL03 photons (Class 42)."""
        
        try:
            with h5.File(atl13_fn, 'r') as f:
                if laser not in f: return df
                
                if f'/{laser}/segment_id_beg' not in f: return df
                
                atl13_seg = f[f'/{laser}/segment_id_beg'][...]
                
                is_water = df['ph_segment_id'].isin(atl13_seg)

                ## Prioritize Bathy (40) over inland water
                mask = is_water & (df['ph_h_classed'] < 40)
                utils.echo_msg(f'Classed {np.count_nonzero(mask)} photons from {laser} based on ATL13')
                df.loc[mask, 'ph_h_classed'] = 42

        except Exception as e:
            utils.echo_warning_msg(f"Failed to apply ATL13 classifications from {atl13_fn}: {e}")
            
        return df
    
    
    def apply_atl24_classifications_full(self, df, atl24_fn, laser, geoseg_beg, geoseg_end, segment_index_dict):
        """Map ATL24 Bathymetry classes using index mapping."""
        
        try:
            with h5.File(atl24_fn, 'r') as f:
                if laser not in f: return df

                grp = f[laser]
                atl24_lat = grp['lat_ph'][...]
                atl24_lon = grp['lon_ph'][...]
                atl24_ortho = grp['ortho_h'][...]
                atl24_ellipse = grp['ellipse_h'][...]
                atl24_surface = grp['surface_h'][...]
                atl24_class = grp['class_ph'][...]
                atl24_seg = grp['index_seg'][...] 
                atl24_idx = grp['index_ph'][...] # 1-based index
                atl24_conf = grp['confidence'][...]

                ## Map internal ATL24 segment index to Real Segment ID
                orig_segs = np.arange(geoseg_beg, geoseg_end + 1)
                atl24_real_seg_ids = orig_segs[atl24_seg]
                
                ## Filter to relevant data
                mask = np.isin(atl24_real_seg_ids, df['ph_segment_id'].unique())
                if not np.any(mask): return df
                
                ## Calculate Global Indices in DF
                ## Get start index for each photon's segment
                seg_starts = np.array([segment_index_dict.get(s, -1) for s in atl24_real_seg_ids[mask]])
                
                valid_seg_lookup = seg_starts != -1
                
                ## Global Index = SegmentStart + (PhotonIndex - 1)
                ## Filter down to only those where lookup succeeded
                subset_idx = atl24_idx[mask][valid_seg_lookup]
                subset_starts = seg_starts[valid_seg_lookup]
                
                target_indices = subset_starts + (subset_idx - 1)
                
                ## Boundary check
                valid_bounds = (target_indices >= 0) & (target_indices < len(df))
                final_indices = target_indices[valid_bounds]

                ## ==============================================
                ## Filter for Bathy Criteria (Class >= 40)
                ## We need to filter the source arrays by [mask][valid_seg_lookup][valid_bounds]
                ## AND check the class/confidence
                ## ==============================================
                
                ## Extract source values corresponding to final_indices
                src_class = atl24_class[mask][valid_seg_lookup][valid_bounds]
                src_conf = atl24_conf[mask][valid_seg_lookup][valid_bounds]
                
                is_bathy = src_class >= 40
                if self.min_bathy_confidence is not None:
                    is_bathy &= (src_conf >= self.min_bathy_confidence)
                
                ## Apply Bathy Filter to indices and values
                bathy_indices = final_indices[is_bathy]
                utils.echo_msg(f'Classed {len(bathy_indices)} photons from {laser} based on ATL24')
                if len(bathy_indices) == 0: return df
                
                ## Update DataFrame
                ## Coordinates (Refracted)
                df.loc[df.index[bathy_indices], 'latitude'] = atl24_lat[mask][valid_seg_lookup][valid_bounds][is_bathy]
                df.loc[df.index[bathy_indices], 'longitude'] = atl24_lon[mask][valid_seg_lookup][valid_bounds][is_bathy]
                df.loc[df.index[bathy_indices], 'photon_height'] = atl24_ortho[mask][valid_seg_lookup][valid_bounds][is_bathy] # Assuming Ortho default
                
                ## Extra heights, etc.
                ## df.loc[bathy_indices, 'h_ellipse'] = atl24_ellipse[...]
                
                ## Classification & Confidence
                df.loc[df.index[bathy_indices], 'ph_h_classed'] = src_class[is_bathy]
                df.loc[df.index[bathy_indices], 'bathy_confidence'] = src_conf[is_bathy]

        except Exception as e:
            utils.echo_warning_msg(f"Failed to apply ATL24 data from {atl24_fn}: {e}")

        return df


    def apply_atl24_classifications(self, df, atl24_fn, laser, geoseg_beg, geoseg_end):
        """Map ATL24 Bathymetry classes to ATL03 using Pandas Merge."""
        
        try:
            with h5.File(atl24_fn, 'r') as f:
                if laser not in f: return df

                grp = f[laser]
                
                ## Read ATL24 Arrays
                try:
                    atl24_class = grp['class_ph'][...]
                    atl24_seg = grp['index_seg'][...] 
                    atl24_idx = grp['index_ph'][...]
                    atl24_conf = grp['confidence'][...]
                    
                    # Refracted Coordinates
                    atl24_lat = grp['lat_ph'][...]
                    atl24_lon = grp['lon_ph'][...]
                    atl24_z = grp['ortho_h'][...] # Default to Ortho
                except KeyError:
                    return df

                ## Reconstruct Real Segment IDs
                ## ATL24 'index_seg' is relative to the granule's geoseg_beg
                orig_segs = np.arange(geoseg_beg, geoseg_end + 1)
                
                ## Ensure index_seg doesn't exceed orig_segs bounds (subsetting issue)
                ## If harmony subsets, index_seg might still be 0-based relative to the subset.
                ## Usually it's relative to the original granule start.
                ## If mismatch persists, we rely on the fact that we only keep what matches 'df'.
                try:
                    atl24_real_seg_ids = orig_segs[atl24_seg]
                except IndexError:
                    ## Fallback
                    ## If orig_segs is too small, the subset metadata might be wrong.
                    ## We skip this file if we can't align segments.
                    utils.echo_warning_msg(f"ATL24 Segment Index mismatch in {atl24_fn}")
                    return df

                ## Build ATL24 DataFrame
                atl24_df = pd.DataFrame({
                    'ph_segment_id': atl24_real_seg_ids,
                    'ph_index_within_seg': atl24_idx,
                    'atl24_class': atl24_class,
                    'atl24_conf': atl24_conf,
                    'atl24_lat': atl24_lat,
                    'atl24_lon': atl24_lon,
                    'atl24_z': atl24_z
                })

                ## Filter ATL24 DF to only relevant segments
                atl24_df = atl24_df[atl24_df['ph_segment_id'].isin(df['ph_segment_id'].unique())]

                if self.verbose:
                    utils.echo_msg(f'Classified {len(atl24_df)} photons from {laser} based on ATL24')
                
                if atl24_df.empty: return df

                ## Filter for Bathy Criteria (Class >= 40)
                ## We do this BEFORE merge to keep the merge small
                is_bathy = atl24_df['atl24_class'] >= 40
                if self.min_bathy_confidence is not None:
                    is_bathy &= (atl24_df['atl24_conf'] >= self.min_bathy_confidence)
                
                atl24_df = atl24_df[is_bathy]
                
                ## --- MERGE ---
                ## We merge 'left' to keep all ATL03 points, adding info where ATL24 matches
                ## Keys: SegmentID + PhotonIndex
                merged = df.merge(
                    atl24_df, 
                    on=['ph_segment_id', 'ph_index_within_seg'], 
                    how='left'
                )
                
                ## Update Main DF where match found (atl24_class is not NaN)
                mask = merged['atl24_class'].notna()
                
                if np.any(mask):
                    # Direct assignment using the mask ensures alignment
                    df.loc[mask, 'ph_h_classed'] = merged.loc[mask, 'atl24_class'].astype(int)
                    df.loc[mask, 'bathy_confidence'] = merged.loc[mask, 'atl24_conf'].astype(int)
                    
                    # Update Coordinates (Refraction Correction)
                    df.loc[mask, 'latitude'] = merged.loc[mask, 'atl24_lat']
                    df.loc[mask, 'longitude'] = merged.loc[mask, 'atl24_lon']
                    df.loc[mask, 'photon_height'] = merged.loc[mask, 'atl24_z']

        except Exception as e:
            utils.echo_warning_msg(f"Failed to apply ATL24 data from {atl24_fn}: {e}")

        return df

    
    def classify_bathymetry_algo(self, df):
        """Apply algorithmic bathymetry classification (Known Raster or DBSCAN)."""
        
        ## Known Bathymetry
        if self.known_bathymetry:
            try:
                if self.verbose:
                    utils.echo_msg(f"Classifying using known bathymetry: {self.known_bathymetry}")
                
                ## Create a temporary point object for querying
                ## We use 'x' and 'y' which were renamed from longitude/latitude in yield_points
                ## But inside read_atl03/classify pipeline, columns are still 'longitude'/'latitude'
                pts = np.rec.fromarrays(
                    [df['longitude'].values, df['latitude'].values], 
                    names=['x', 'y']
                )
                
                ## Sample the raster
                ref_z = gdalfun.gdal_query(pts, self.known_bathymetry, 'g').flatten()
                
                ## Calculate Difference (Z_photon - Z_known)
                ## Ensure Z datum matches (Usually requires 'geoid' or 'mean_tide' alignment)
                diff = np.abs(df['photon_height'].values - ref_z)
                
                ## Keep points within threshold
                ## Only update points that aren't already classified as something specific (like land/canopy)
                ## Or, forcibly overwrite if they look like bathy.
                ## Here we target unclassified (-1) or noise (0/1)
                is_bathy = diff <= self.known_bathy_threshold
                
                ## Assign Class 40 (Bathymetry)
                #update_mask = is_bathy & (ph_h_classed < 40)
                #ph_h_classed[update_mask] = 40

                mask = is_bathy & (df['ph_h_classed'] >= 40)
                df.loc[mask, 'ph_h_classed'] = 40
                
                if self.verbose:
                    utils.echo_msg(f"  Identified {np.count_nonzero(mask)} bathy points via known raster.")
                    
            except Exception as e:
                utils.echo_warning_msg(f"Known bathymetry classification failed: {e}")

        ## DBSCAN Clustering
        if self.use_dbscan:
            try:
                from sklearn.cluster import DBSCAN
                from sklearn.preprocessing import StandardScaler
            except ImportError:
                utils.echo_warning_msg("scikit-learn required for DBSCAN. Skipping.")
                return df

            if self.verbose:
                utils.echo_msg("Classifying using DBSCAN...")

            ## Filter to potential subsurface points to speed up clustering
            ## e.g., only look at points with valid confidence or below water surface
            ## For now, we cluster everything that isn't already Land (1) or Canopy (2/3)
            ## Maybe, if we do this after setting to 44, only use those, since they are
            ## photons over water.
            mask_candidates = df['ph_h_classed'] <= 1
            if np.count_nonzero(mask_candidates) < self.dbscan_min_samples:
                return df
                
            subset = df[mask_candidates].copy()
            
            ## Along-track distance (approx via index or lat) and Height
            ## We must normalize because Dist is in km/m and Height is in m
            X = np.column_stack((
                np.arange(len(subset)), # Proxy for along-track distance
                subset['photon_height'].values
            ))
            
            ## Standardize features
            X_scaled = StandardScaler().fit_transform(X)
            
            ## Run DBSCAN
            ## eps is in "standard deviations" now due to scaling. 
            ## Adjust eps logic if physical meters are required (requires custom metric).
            ## For simplicity with StandardScaler, eps=0.3 to 0.5 usually captures dense lines.
            db = DBSCAN(eps=0.3, min_samples=self.dbscan_min_samples).fit(X_scaled)
            labels = db.labels_
            
            ## Identify the "Bathymetry" Cluster
            ## Bathymetry is the largest cluster (mode) found below sea level, maybe...
            ## Or maybe all dense clusters are signal.
            unique_labels = set(labels)
            if -1 in unique_labels: unique_labels.remove(-1) # Remove noise label
            
            for k in unique_labels:
                class_member_mask = (labels == k)
                
                ## Check if the cluster is "bathy-like"?
                ## Check mean height (should be < 0 for bathy usually)
                cluster_z_mean = np.mean(subset.loc[class_member_mask, 'photon_height'])
                
                ## If it looks like bathy (e.g. below 0), assign it
                ## Here we aggressively classify ALL dense clusters as signal (Class 40)
                ## Users can filter later.
                if cluster_z_mean < 0: 
                    original_indices = subset.index[class_member_mask]
                    df.loc[original_indices, 'ph_h_classed'] = 40
                
            if self.verbose:
                utils.echo_msg(f"  DBSCAN found {len(unique_labels)} dense clusters.")

        return df

    
    def classify_buildings_algo(self, df, min_height=2.5, max_roughness=0.6, ground_window=50):
        """Classify buildings by identifying 'Flat & Elevated' features.
        
        Args:
            df (pd.DataFrame): Photon dataframe.
            min_height (float): Minimum height above ground (m) to consider.
            max_roughness (float): Max StdDev (m) for a feature to be a 'Roof'.
                                   Trees are usually > 1.0m. Roofs < 0.5m.
            ground_window (int): Number of photons for rolling ground estimation.
        """
        
        try:
            if self.verbose:
                utils.echo_msg("Running Dynamic Building Classification...")

            ## Estimate Ground Surface (Rolling 5th Percentile)
            ## This is a rough "Digital Terrain Model" (DTM) derived from the data itself.
            ## We assume the lowest points in a window are ground.
            ## Using a centered window handles slopes better.
            ground_proxy = df['photon_height'].rolling(window=ground_window, center=True, min_periods=5).quantile(0.05)
            
            ## Fill edges where rolling is NaN
            ground_proxy = ground_proxy.bfill().ffill()
            
            ## Calculate Height Above Ground (HAG)
            hag = df['photon_height'] - ground_proxy
            
            ## Identify Elevated Candidates (Exclude Ground & Low Noise)
            ## Must be high confidence to avoid atmospheric noise
            is_elevated = (hag >= min_height) & (df['confidence'] >= 3)
            
            if not np.any(is_elevated): return df

            ## Calculate Roughness (StdDev) on Elevated Points ONLY
            ## We create a subset to calculate roughness, otherwise the gaps between 
            ## buildings and ground would artificially inflate the std dev.
            elevated_df = df[is_elevated].copy()
            
            ## Use a small window (e.g., 10 photons) to check local surface consistency
            elevated_df['local_roughness'] = elevated_df['photon_height'].rolling(window=10, center=True).std()
            
            ## Classify Buildings (Class 7)
            ## levated AND Smooth (Low Roughness)
            ## We filter out NaN roughness (edges)
            is_building = (elevated_df['local_roughness'] <= max_roughness)
            
            building_indices = elevated_df.index[is_building]
            
            if len(building_indices) > 0:
                ## Apply to Main DataFrame
                ## Only overwrite if not already Bathy(40), Coast(41), Water(42)
                ## We overwrite Land(1), Canopy(2/3), or Unclassified
                mask = df.index.isin(building_indices) & (df['ph_h_classed'] < 40)

                ## 7 is standard LAS class for Low Point/Noise, but standard LAS for Building is actually 6. 
                ## But ATL06 uses 6 for Land Ice. 
                df.loc[mask, 'ph_h_classed'] = 7
                
                if self.verbose:
                    utils.echo_msg(f"  Classified {len(building_indices)} photons as Buildings (Algo).")

        except Exception as e:
            utils.echo_warning_msg(f"Building classification failed: {e}")
            
        return df
    

    def classify_nearshore_roughness(self, df, height_window=5.0, max_roughness=1.5):
        """Classify nearshore water by analyzing segment roughness (StdDev).
        
        Segments near the Geoid (0m) with low roughness are likely water.
        
        Args:
            df (pd.DataFrame): The photon dataframe.
            height_window (float): Max vertical distance from 0m (Geoid) to consider.
                                   Use ~5m to account for tides/surge/setup.
            max_roughness (float): Max Standard Deviation (m) to be considered water.
                                   Calm ocean is < 0.3m. Surf zone can be 1.0-1.5m.
        """
        
        try:
            if self.verbose:
                utils.echo_msg("Running Dynamic Nearshore Roughness Classification...")

            ## Select "Signal" Photons for Statistics
            ## We must exclude noise (conf < 2) so it doesn't inflate the std dev
            signal_mask = df['confidence'] >= 3  # Medium or High confidence
            if not np.any(signal_mask): return df
            
            signal_df = df[signal_mask]

            ## Calculate Segment Statistics (20m resolution)
            ## Group by 'ph_segment_id' to analyze each 20m strip of the track
            grouped = signal_df.groupby('ph_segment_id')['photon_height']
            #grouped = signal_df.groupby('ph_segment_id')['photon_meantide']
            
            ## We need Median (robust center) and Std (roughness)
            seg_stats = grouped.agg(['median', 'std'])
            
            ## Identify Water Segments
            ## Height is within window of 0m (Geoid)
            ## Roughness is low (Flat surface)
            is_water_seg = (
                (seg_stats['median'].abs() <= height_window) & 
                (seg_stats['std'] <= max_roughness)
            )
            
            valid_water_segs = seg_stats.index[is_water_seg]
            
            if len(valid_water_segs) == 0: return df

            ## Apply Classification to Original DataFrame
            ## We map the identified segments back to all photons (including lower confidence ones)
            is_water_photon = df['ph_segment_id'].isin(valid_water_segs)
            
            ## Update Class (Use 41 for Nearshore/Coastline)
            ## Only update if not already classified as Bathy (40) or Inland Water (42)
            ## We overwrite Unclassified (-1), Noise (0), Land (1), etc.
            mask = is_water_photon & (df['ph_h_classed'] < 40)
            
            ## Assign Class 41 (Coastline/Nearshore)
            df.loc[mask, 'ph_h_classed'] = 41
            
            if self.verbose:
                utils.echo_msg(f"  Classified {np.count_nonzero(mask)} photons as Nearshore Water (41).")

        except Exception as e:
            utils.echo_warning_msg(f"Nearshore roughness classification failed: {e}")
            
        return df
    
    
    def read_atl03(self, f, laser_num, orientation=None, atl08_fn=None, atl24_fn=None, atl06_fn=None, atl12_fn=None, atl13_fn=None):
        """Read and classify data from an ATL03 file/laser."""
        
        if orientation is None:
            orientation = f['/orbit_info/sc_orient'][0]
            
        laser = 'gt' + laser_num + self.orientDict[orientation]
        if laser not in f or 'heights' not in f[laser]: return None

        try:
            h_grp = f[f'/{laser}/heights']
            geo_grp = f[f'/{laser}/geolocation']
            geophys_grp = f[f'/{laser}/geophys_corr']
            anc = f['ancillary_data']
            
            ## Read Base Arrays
            lat = h_grp['lat_ph'][...]
            lon = h_grp['lon_ph'][...]
            h_ph = h_grp['h_ph'][...]
            conf = h_grp['signal_conf_ph'][..., 0] 
            dt = h_grp['delta_time'][...]
            
            ## Geolocation / Indexing
            seg_ph_cnt = geo_grp['segment_ph_cnt'][...]
            seg_id = geo_grp['segment_id'][...]
            geoseg_beg = anc['start_geoseg'][0]
            geoseg_end = anc['end_geoseg'][0]

            ## shape is (n_segments, 5). Columns: 0=Land, 1=Ocean, 2=SeaIce, 3=LandIce, 4=InlandWater
            surf_type = geo_grp['surf_type'][...]
            
            ## Corrections
            geoid = geophys_grp['geoid'][...]
            geoid_f2m = geophys_grp['geoid_free2mean'][...]
            dem_h = geophys_grp['dem_h'][...]

        except KeyError: return None

        ## --- Indexing ---
        ## mrl: move these down for merging.
        #seg_starts = np.concatenate(([0], np.cumsum(seg_ph_cnt)[:-1]))
        #seg_idx_dict = dict(zip(seg_id, seg_starts))

        ## Map Photon -> Segment ID
        ph_seg_ids = np.repeat(seg_id, seg_ph_cnt)

        ## Extract Ocean Column (Index 1) and expand to photon length
        ## 1 means the mask indicates Ocean
        seg_is_ocean = surf_type[:, 1]
        ph_is_ocean = np.repeat(seg_is_ocean, seg_ph_cnt)
        
        ## Truncate
        min_len = min(len(ph_seg_ids), len(h_ph))
        ph_seg_ids = ph_seg_ids[:min_len]
        lat = lat[:min_len]; lon = lon[:min_len]; h_ph = h_ph[:min_len]
        conf = conf[:min_len]; dt = dt[:min_len]

        ## --- Generate "Index Within Segment" for Merging ---
        ## This is the magic key that allows us to merge with ATL24 reliably.
        ## Since ph_seg_ids is sorted (repeated), we can generate counters.
        ## * Find where segment changes
        ## * Reset counter
        ##
        ## Calculate indices
        ## We want: 1, 2, 3... for seg A; 1, 2... for seg B
        ## sequence = np.arange(len(ph_seg_ids))
        ## Let's use the segment counts we already have: seg_ph_cnt
        ## We need to repeat [1..count] for each segment.
        
        ## Adjust seg_ph_cnt to match the truncation we did above
        ## If we truncated ph_seg_ids, we need to fix the count of the last segment
        if len(ph_seg_ids) < np.sum(seg_ph_cnt):
            ## Re-calculate counts based on the truncated array
            unique, counts = np.unique(ph_seg_ids, return_counts=True)
            ## This ensures 'counts' matches the actual data we kept
            ph_index_counters = np.concatenate([np.arange(1, c + 1) for c in counts])
        else:
            ph_index_counters = np.concatenate([np.arange(1, c + 1) for c in seg_ph_cnt])
            ## Safety clip if original counts exceeded data
            ph_index_counters = ph_index_counters[:min_len]
        
        ## --- Calculate Heights ---
        h_geoid_map = dict(zip(seg_id, geoid))
        h_f2m_map = dict(zip(seg_id, geoid_f2m))
        h_dem_map = dict(zip(seg_id, dem_h))

        p_geoid = np.array([h_geoid_map.get(s, 0) for s in ph_seg_ids])
        p_f2m = np.array([h_f2m_map.get(s, 0) for s in ph_seg_ids])
        p_dem = np.array([h_dem_map.get(s, 0) for s in ph_seg_ids])

        h_ortho = h_ph - p_geoid
        h_meantide = h_ph - (p_geoid + p_f2m)
        h_dem = p_dem - (p_geoid + p_f2m)

        ## --- Initialize Classifications ---
        ## Start with -1 (Unclassified)
        ph_class = np.full(len(h_ph), -1, dtype=int)
        
        ## --- Apply Fallback Ocean Classification ---
        ## Onboard Mask says Ocean (ph_is_ocean == 1)
        ## Height is near Sea Level (abs(h_ortho) < threshold). 
        ## Open ocean should be near 0. We use 30m to be very safe against tides/surge/geoid error.
        is_open_ocean = (ph_is_ocean == 1) & (np.abs(h_ortho) < 30)
        might_be_open_ocean = (np.abs(h_ortho) < 5)                    
        
        ## --- Create DataFrame ---
        ## Select Height
        if self.water_surface == 'mean_tide': z_out = h_meantide
        elif self.water_surface == 'geoid': z_out = h_ortho
        else: z_out = h_ph

        df = pd.DataFrame({
            'latitude': lat, 'longitude': lon, 'photon_height': z_out,
            'laser': laser, 'fn': self.fn, 'confidence': conf, 'delta_time': dt,
            'photon_h_dem': h_dem, 'photon_meantide': h_meantide,
            'ph_h_classed': -1, # Init Unclassified
            'bathy_confidence': -1, # Init
            'ph_segment_id': ph_seg_ids, # Critical for linking
            'ph_index_within_seg': ph_index_counters # join key
        })

        ## --- Apply Classifications based on ATLXX data ---
        ## ATL08,Land/Canopy, Ground (1), Canopy (2), Top Canopy (3)
        seg_starts = np.concatenate(([0], np.cumsum(seg_ph_cnt)[:-1]))
        seg_idx_dict = dict(zip(seg_id, seg_starts))
        
        if atl08_fn:
            df = self.apply_atl08_classifications(df, atl08_fn, laser, seg_idx_dict)

        ## Assign Class 44 (Ocean)
        if self.verbose:
            utils.echo_msg(f'Classified {np.count_nonzero(is_open_ocean)} photons to open ocean')            
        df.loc[is_open_ocean, 'ph_h_classed'] = 44
            
        ## ATL24,Bathymetry, Near-shore Bathymetry (40), Coastline (41)
        if atl24_fn:
            df = self.apply_atl24_classifications(
                df, atl24_fn, laser, geoseg_beg, geoseg_end
            )
            
        ## ATL13,Inland Water,"Lakes, Rivers, Reservoirs",42
        if atl13_fn:
            df = self.apply_atl13_classifications(df, atl13_fn, laser)

        ## ATL12,Ocean,Open Ocean Surface,44
        if atl12_fn:
            df = self.apply_atl12_classifications(df, atl12_fn, laser)

        # ## ATL06,Land Ice,"Glaciers, Ice Sheets",6
        # if atl06_fn:
        #     df = self.apply_atl06_classifications(df, atl06_fn, laser)
            
        ## --- Apply Dynamic Nearshore Classification ---
        ## This fills gaps between ATL12 (Ocean) and ATL08 (Land)
        ## It is robust in the surf zone where ATL12 often drops out.
        if self.classify_water:
            ## Use wider window (5m) and higher roughness (1.5m) for surf zones
            df = self.classify_nearshore_roughness(df, height_window=5.0, max_roughness=1.5)

        if self.classify_buildings:
            df = self.classify_buildings_algo(df)
            
        ## --- Apply Algorithmic Bathymetry Classification (Known Bathy / DBSCAN) ---
        if self.known_bathymetry or self.use_dbscan:
            df = self.classify_bathymetry_algo(df)

        return df

    
    def yield_points(self):
        """Pipeline to yield classified points."""

        ## Initialize external mask geometries
        bing_geom = None
        osm_geom = None
        osm_lakes = None
        
        ## Setup Region for Aux Fetches
        ## We try to use the user-supplied region or the INF metadata bounds
        # proc_region = self._get_processing_region()
        # if proc_region is not None:
        #     proc_region.buffer(pct=1) # Slight buffer for safety        

        #     ## External Masks (Buildings/Water)
        #     if self.classify_buildings:
        #         bfp = bingbfp.BingBuildings(region=proc_region, verbose=self.verbose, cache_dir=self.cache_dir)
        #         bing_geom = bfp(return_geom=True)

        #     if self.classify_water:
        #         osmc = osm.osmCoastline(region=proc_region, q='coastline', verbose=self.verbose, cache_dir=self.cache_dir)
        #         osm_geom = osmc(return_geom=True)

        #     if self.classify_inland_water:
        #         osml = osm.osmCoastline(region=proc_region, q='water', verbose=self.verbose, cache_dir=self.cache_dir)
        #         osm_lakes = osml(return_geom=True)

        ## Fetch Aux ATLXX Data
        atl08_fn = self.fetch_atlxx(self.fn, 'ATL08') if self.classes else None
        atl24_fn = self.fetch_atlxx(self.fn, 'ATL24') if self.classes else None
        # atl06_fn = self.fetch_atlxx(self.fn, 'ATL06') if self.classes else None
        atl12_fn = self.fetch_atlxx(self.fn, 'ATL12') if self.classes else None
        atl13_fn = self.fetch_atlxx(self.fn, 'ATL13') if self.classes else None
        atl06_fn = None
        #atl12_fn = None
        #atl13_fn = None
        if self.verbose:
            if atl08_fn:
                utils.echo_msg(f'Using {atl08_fn} for classifications')
            if atl24_fn:
                utils.echo_msg(f'Using {atl24_fn} for classifications')
            if atl06_fn:
                utils.echo_msg(f'Using {atl06_fn} for classifications')
            if atl12_fn:
                utils.echo_msg(f'Using {atl12_fn} for classifications')
            if atl13_fn:
                utils.echo_msg(f'Using {atl13_fn} for classifications')
            
        ## Process Granule
        with h5.File(self.fn, 'r') as f:
            ## QA Check
            if self.reject_failed_qa and 'quality_assessment' in f:
                if f['/quality_assessment/qa_granule_pass_fail'][0] != 0:
                    utils.echo_warning_msg(f"Skipping failed granule {self.fn}")
                    return

            for i in range(1, 4):
                for orient in range(2):
                    ## Read_atl03 reads all atl03 data and classifies the
                    ## photons (in 'ph_h_classed') based on aux ATL data.
                    ## This creates the pandas dataframe.
                    #utils.echo_debug_msg('Loading ATL data into data-frame...')
                    dataset = self.read_atl03(
                        f, str(i), orientation=orient,
                        atl08_fn=atl08_fn, atl24_fn=atl24_fn, atl06_fn=atl06_fn
                    )
                    
                    if dataset is None or dataset.empty: continue

                    ## Filter by Confidence
                    if self.confidence_levels:
                        dataset = dataset[dataset['confidence'].isin(self.confidence_levels)]

                    ## Spatial Subset (Coarse)
                    ## Only apply if user supplied a region (self.region).
                    if self.region:
                        dataset = dataset[
                            (dataset['longitude'] >= self.region.xmin) & 
                            (dataset['longitude'] <= self.region.xmax) & 
                            (dataset['latitude'] >= self.region.ymin) & 
                            (dataset['latitude'] <= self.region.ymax)
                        ]

                    if dataset.empty: continue
                    
                    ## We classify buildings (7), ocean (41) and inland water (42) - maybe 44 for ocean?
                    ## using Bing BFP and OpenStreetMap
                    if self.classes:
                        if bing_geom:
                            dataset = self.classify_by_mask_geoms(dataset, bing_geom, 7)
                        if osm_geom:
                            dataset = self.classify_by_mask_geoms(dataset, osm_geom, 41, except_classes=[40])
                        if osm_lakes:
                            dataset = self.classify_by_mask_geoms(dataset, osm_lakes, 42, except_classes=[40])

                        ## Final Class Filter
                        dataset = dataset[dataset['ph_h_classed'].isin(self.classes)]

                    if dataset.empty: continue

                    dataset.rename(columns={'longitude': 'x', 'latitude': 'y', 'photon_height': 'z'}, inplace=True)
                    yield dataset

                    
    def _vectorize_df(self, dataset):
        """Convert dataframe to OGR Memory Layer."""
        
        ds = gdal.GetDriverByName('Memory').Create('', 0, 0, 0, gdal.GDT_Unknown)
        lyr = ds.CreateLayer('pts', geom_type=ogr.wkbPoint)
        lyr.CreateField(ogr.FieldDefn('index', ogr.OFTInteger))
        
        ## Use transaction for speed
        lyr.StartTransaction()
        feat = ogr.Feature(lyr.GetLayerDefn())
        
        for idx, row in dataset.iterrows():
            feat.SetField(0, int(idx))
            feat.SetGeometryDirectly(ogr.CreateGeometryFromWkt(f"POINT ({row['longitude']} {row['latitude']})"))
            lyr.CreateFeature(feat)
            
        lyr.CommitTransaction()
        return ds

    
    def classify_by_mask_geoms(self, dataset, mask_geoms, classification, except_classes=[]):
        """Spatial classification using OGR masks."""
        
        if not mask_geoms: return dataset
        
        ogr_ds = self._vectorize_df(dataset)
        layer = ogr_ds.GetLayer()
        
        ## Iterate masks, set filter, update dataframe indices found
        for geom in mask_geoms:
            if isinstance(geom, str):
                geom = ogr.CreateGeometryFromWkt(geom)

            layer.SetSpatialFilter(geom)
            
            indices_to_update = []
            for feat in layer:
                indices_to_update.append(feat.GetField('index'))
            
            if indices_to_update:
                # Update only if not in protected classes
                mask = dataset.index.isin(indices_to_update) & (~dataset['ph_h_classed'].isin(except_classes))
                utils.echo_msg(f'Classed {np.count_nonzero(mask)} photons from based to {classification}')
                dataset.loc[mask, 'ph_h_classed'] = classification
                
            layer.SetSpatialFilter(None)
            
        return dataset

### End
