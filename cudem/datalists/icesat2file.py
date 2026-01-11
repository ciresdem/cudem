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
## This module provides a robust pipeline for parsing, classifying, and gridding 
## ICESat-2 photon data. It handles the full hierarchy of ICESat-2 products, 
## merging raw photon geolocation (ATL03) with higher-level classifications 
## (ATL08, ATL12, ATL24) and augmenting them with dynamic, physics-based 
## algorithms for feature detection.
##
## Supported Data Products:
## - ATL03: Global Geolocated Photon Data (Base Input)
## - ATL08: Land and Vegetation Height (Classifies Ground vs Canopy)
## - ATL09: Atmospheric Characteristics (Provides Apparent Surface Reflectance)
## - ATL12: Ocean Surface Height (Open Ocean classification)
## - ATL13: Inland Water Body Height (Lakes, Rivers, Reservoirs)
## - ATL24: Shallow Water Bathymetry (Refracted Subsurface Photon Data)
##
## Dynamic Algorithmic Classification:
## When standard products are unavailable or insufficient, this module applies 
## self-contained algorithms to detect features based on photon geometry and radiometry:
##
## 1. Nearshore Water (Surf Zone):
##    Identifies water surfaces in complex surf zones where standard algorithms fail.
##    Uses a "Roughness + Radiometry" approach: segments near the Geoid that are 
##    either "Flat" (calm) or "Rough & Bright" (whitewater/surf) are classified as water.
##
## 2. Inland Water (Flat & Dark):
##    Detects lakes and rivers without external masks by identifying segments that 
##    are geometrically flat (low StdDev) and radiometrically dark (low Reflectance),
##    distinguishing them from bright flat features like roads.
##
## 3. Buildings (Flat, Elevated & Bright):
##    Classifies structures without external footprints. Identifies photons that are:
##    - Elevated above a calculated ground surface (HAG > 2.5m).
##    - Geometrically "Flat" (Roof) vs "Rough" (Vegetation).
##    - Radiometrically "Bright" (Concrete/Metal) to rescue complex roofs.
##
## 4. Bathymetry (Subsurface):
##    - Supervised: Validates potential bathymetry against a known raster (e.g., GMRT/GEBCO).
##    - Unsupervised: Uses DBSCAN clustering to identify dense subsurface signal 
##      clusters within noisy environments.
##
## 5. Pseudo-Reflectance:
##    If ATL09 is missing, estimates a "Pseudo-Reflectance" proxy based on 
##    signal photon density per shot, normalized by beam strength (Strong/Weak),
##    enabling radiometric classification in standalone ATL03 files.
##
## Classes:
## 0: Noise (if enabled)
## 1: Ground (ATL08)
## 2: Canopy (ATL08)
## 3: Top Canopy (ATL08)
## 6: Land Ice (ATL06)
## 7: Buildings (Dynamic Algo / Bing Mask)
## 40: Bathymetry (ATL24 / Dynamic Algo)
## 41: Coastline / Nearshore Water (ATL24 / Dynamic Algo)
## 42: Inland Water (ATL13 / Dynamic Algo)
## 44: Open Ocean (ATL12 / Geoid Fallback)
## -1: Unclassified
##
### Code:

import os
import glob
import numpy as np
import h5py as h5
import pandas as pd

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
class IceSat2_ATL03(ElevationDataset):
    """ICESat-2 ATL03 (Global Geolocated Photon Data) Parser.
    
    This module ingests ATL03 photon data and fuses it with higher-level products 
    (ATL08, ATL12, ATL24) or external references (known bathymetry, building masks). 
    It also applies dynamic, physics-based algorithms to classify features (surf, 
    inland water, buildings) when standard products are unavailable.

    Ground (1), Canopy (2), Top Canopy (3), Bathymetry (40), Water Surface (41), Land Ice (6), Inland Water (42)
    Open Ocean Surface (44), Atmosphere (0), Buildings (7), Unclassified (-1)

    Parameters:
    -----------
        water_surface (str): Vertical datum for photon heights. 
                             Options: 'geoid' (EGM2008), 'mean_tide', 'ellipsoid'. Default: 'mean_tide'.
        classes (list/str): Filter output to specific classification codes. 
                            Can be a list of ints or a '/' separated string (e.g., '1/40/41').
                            Classes: 1=Ground, 2=Canopy, 6=Ice, 7=Building, 40=Bathy, 41=Nearshore, 
                                     42=Inland Water, 44=Ocean.
        confidence_levels (list/str): Filter photons by ATL03 signal confidence.
                                      Options: 0=Noise, 1=Buffer, 2=Low, 3=Medium, 4=High. Default: '2/3/4'.
        reject_failed_qa (bool): If True, skip granules where 'qa_granule_pass_fail' is 0. Default: True.
        append_atl24 (bool): If True, append ATL24 bathymetry data if found. Default: False.
        min_bathy_confidence (int/float): Minimum confidence threshold for ATL24 photons.
        known_bathymetry (str): Path to a GDAL-readable raster (e.g., 'gebco.tif') or 'gmrt' 
                                to use for Supervised Bathymetry classification.
        known_bathy_threshold (float): Vertical window (+/- meters) around the known bathymetry 
                                       surface to classify photons as Class 40 (Bathymetry). Default: 5.0.
        use_dbscan (bool): If True, use unsupervised DBSCAN clustering to detect subsurface 
                           features (bathymetry) in noisy data. Default: False.
        dbscan_eps (float): The maximum distance between two samples for one to be considered 
                            as in the neighborhood of the other (DBSCAN parameter). Default: 1.5.
        dbscan_min_samples (int): The number of samples (or total weight) in a neighborhood 
                                  for a point to be considered as a core point (DBSCAN parameter). Default: 10.
        **kwargs: Additional keyword arguments passed to the parent ElevationDataset.
    """

    def __init__(self,
                 water_surface='geoid', # height
                 classes=None, # classify the data
                 confidence_levels='2/3/4', # filter data by confidence level
                 columns=None,
                 reject_failed_qa=True, # reject bad granules
                 append_atl24=False, # append atl24 data to yield_points
                 min_bathy_confidence=None, # minimum desired bathy confidence (from atl24)
                 use_external_masks=False, # use external masks, e.g. bingbfp, osm
                 known_bathymetry=None, # Path to raster or 'gmrt'
                 known_bathy_threshold=5.0, # Vertical window (+/- meters)
                 use_dbscan=False, #  Enable DBSCAN clustering
                 dbscan_eps=1.5, #  DBSCAN epsilon (neighbor distance)
                 dbscan_min_samples=10, #  DBSCAN min points for core
                 **kwargs):
        
        super().__init__(**kwargs)
        
        self.data_format = 303
        self.water_surface = water_surface if water_surface in ['mean_tide', 'geoid', 'ellipsoid'] else 'mean_tide'

        self.classes = [int(x) for x in str(classes).split('/')] if classes is not None else []
        self.confidence_levels = [int(x) for x in str(confidence_levels).split('/')] if confidence_levels is not None else []
        self.columns = columns if columns is not None else {}
        
        self.reject_failed_qa = reject_failed_qa
        self.append_atl24 = append_atl24
        self.min_bathy_confidence = utils.float_or(min_bathy_confidence)

        self.use_external_masks = use_external_masks
        
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
    
    
    ## ==============================================
    ## Fetch AUX ATL* Data
    ## ==============================================    
    def fetch_atlxx(self, atl03_fn, short_name='ATL08'):
        """Fetch associated ATLxx file (ATL08, ATL24, ATL06, etc)."""
        
        ## Build filename patterns
        bn = utils.fn_basename2(atl03_fn)
        parts = bn.split('_')
        if len(parts) < 4: return None
        
        atlxx_filter = '_'.join(parts[1:4])
        atlxx_filter_no_ver = '_'.join(parts[1:3])

        if self.verbose:
            utils.echo_msg(f'Fetching {short_name}: {atlxx_filter}') 
        
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
            else:
                if self.verbose:
                    utils.echo_warning_msg(f'  Could not locate {short_name}: {atlxx_filter}')
        
        return None


    ## ==============================================
    ## Add 'reflectance' to ATL03, either using ATL09 if it's
    ## available, or calculated 'pseudo-reflectance' based
    ## on the density of returns...
    ## ==============================================    
    def apply_atl09_data(self, df, atl09_fn, laser):
        """Map ATL09 Apparent Surface Reflectance to ATL03 photons.
        
        ATL09 is organized by 'Profile' (1, 2, 3), typically representing 
        the 3 Strong Beams or Ground Tracks. We must map the current 
        ATL03 laser (e.g., 'gt1l') to the correct ATL09 Profile.
        
        Note: If ATL09 only contains Strong beams, Weak beam photons 
        will receive the reflectance of their corresponding Strong beam 
        (which is usually a valid approximation for atmosphere/reflectance).
        """
        
        try:
            with h5.File(atl09_fn, 'r') as f:
                ## Identify which Profile matches our Laser
                ## We iterate profiles and check their metadata
                target_profile = None
                
                ## Extract GT number (gt1, gt2, gt3) from laser string 'gt1l'
                current_gt = laser[:3] 
                
                for p_num in range(1, 4):
                    profile = f'profile_{p_num}'
                    if profile not in f: continue
                    
                    ## Check Lineage/Metadata to find Ground Track
                    ## Usually /profile_x/metadata/orbit_info/ground_track
                    ## Or we check attributes
                    try:
                        ## Some versions store it in /orbit_info/beam_mapping or similar
                        ## But a robust proxy is often just matching the index if standard
                        ## Let's try to find a specific attribute if possible.
                        ## For now, we assume standard mapping:
                        ## Profile 1 = GT1, Profile 2 = GT2, Profile 3 = GT3
                        ## This holds for most nominal operations.
                        if f'{current_gt}' in f'{profile}': # Heuristic check not reliable
                            pass 
                            
                        ## Use simple index mapping: profile_1 -> gt1
                        if f'profile_{current_gt[-1]}' == profile:
                            target_profile = profile
                            break
                    except: pass
                
                ## Map gt1->profile_1, gt2->profile_2, gt3->profile_3
                if target_profile is None:
                    target_profile = f'profile_{laser[2]}'
                
                if target_profile not in f: return df

                ## Read High Rate Data (25Hz)
                ## /profile_x/high_rate/
                grp = f[f'{target_profile}/high_rate']
                
                if 'apparent_surf_reflec' not in grp: return df
                
                reflec = grp['apparent_surf_reflec'][...]
                seg_beg = grp['ds_segment_id_beg'][...] # Start Segment of this 25Hz shot
                seg_end = grp['ds_segment_id_end'][...] # End Segment
                
                ## Handle NoData (usually > 3.4e38)
                reflec[reflec > 1e30] = np.nan
                
                ## Create Lookup Table (Segment ID -> Reflectance)
                ## We essentially "fill" the gaps between beg and end with the value
                
                ## Filter to relevant segments in our dataframe
                min_seg = df['ph_segment_id'].min()
                max_seg = df['ph_segment_id'].max()
                
                ## Filter ATL09 data to overlapping range
                ## Use seg_end >= min and seg_beg <= max
                overlap_mask = (seg_end >= min_seg) & (seg_beg <= max_seg)
                if not np.any(overlap_mask): return df
                
                r_sub = reflec[overlap_mask]
                b_sub = seg_beg[overlap_mask]
                e_sub = seg_end[overlap_mask]
                
                ## Create a dense lookup array for the range of segments in this chunk
                ## Offset by min_seg to keep array size manageable
                lookup_len = max_seg - min_seg + 1
                lookup_arr = np.full(lookup_len, np.nan, dtype=np.float32)
                
                ## Fill lookup
                ## Since multiple 25Hz shots might cover overlapping segments (rare but possible),
                ## or gaps exist, we iterate. Vectorized fill is tricky with variable windows.
                ## Fast loop approach:
                for r_val, s_start, s_stop in zip(r_sub, b_sub, e_sub):
                    # Clip to current window
                    start_idx = max(0, s_start - min_seg)
                    stop_idx = min(lookup_len, s_stop - min_seg + 1)
                    
                    if start_idx < stop_idx:
                        lookup_arr[start_idx:stop_idx] = r_val
                        
                ## --- Map to Photons --- 
                ## Calculate offsets for DF segment IDs
                ph_offsets = df['ph_segment_id'].values - min_seg
                
                # Safety clip (shouldn't be needed due to min/max logic, but safety first)
                valid_offsets = (ph_offsets >= 0) & (ph_offsets < lookup_len)
                
                ## Create column if not exists
                if 'reflectance' not in df.columns:
                    df['reflectance'] = np.nan
                    
                ## We use numpy indexing for speed
                mapped_values = lookup_arr[ph_offsets[valid_offsets]]
                
                ## df.loc[valid_offsets, 'reflectance'] = mapped_values # Logic error in indexing
                df_indices = df.index[valid_offsets]
                df.loc[df_indices, 'reflectance'] = mapped_values
                
        except Exception as e:
            utils.echo_warning_msg(f"Failed to apply ATL09 data: {e}")
            
        return df

    
    def is_strong_beam(self, laser_name, orientation):
        """Check if laser is Strong based on spacecraft orientation.
        Orientation 0 (Backward): Left is Strong.
        Orientation 1 (Forward): Right is Strong.
        """
        
        ## laser_name is like 'gt1l', 'gt2r'
        side = laser_name[-1] # 'l' or 'r'

        if orientation == 0: # Backward
            return side == 'l'
        elif orientation == 1: # Forward
            return side == 'r'
        return False


    def calculate_pseudo_reflectance(self, df, is_strong=True):
        """Estimate 'Pseudo-Reflectance' from Signal Photon density.
        
        - Uses STRICT confidence to prevent solar noise from corrupting 'Dark' water detection.
        - Uses RELAXED clipping (4.0) to capture 'Specular Saturation' (Mirror water).
        """
        
        try:
            ## Determine Thresholds & Divisors
            if is_strong:
                ## Strong Beams: Use Medium+ Confidence (3, 4)
                signal_mask = df['confidence'] >= 3
                divisor = 29.0
                min_photons = 0 
            else:
                ## Weak Beams: Use High Confidence ONLY (4)
                ## Weak beams are noisy; we need strict filtering to ensure 
                ## any detected signal is actually surface, not background noise.
                signal_mask = df['confidence'] >= 4
                divisor = 7.25
                min_photons = 5 # Require cluster to prove existence

            if not np.any(signal_mask): 
                df['reflectance'] = np.nan
                return df
            
            ## Count Signal Photons
            segment_counts = df.loc[signal_mask].groupby('ph_segment_id').size()
            
            ## Apply Floor (Weak Beam Protection)
            ## If a weak beam has < 5 high-conf photons, it's indistinguishable from noise.
            if min_photons > 0:
                segment_counts = segment_counts.where(segment_counts >= min_photons, np.nan)
            
            ## Calculate Density
            pseudo_reflec = segment_counts / divisor
            
            ## Clip Specular Returns
            ## We clip at 4.0 (instead of 1.2) to capture "Saturation Events".
            ## Normal Grass = ~0.6. Specular Water > 1.5. Saturation > 3.0.
            pseudo_reflec = pseudo_reflec.clip(upper=4.0) 
            
            ## Map back to DataFrame
            df['reflectance'] = df['ph_segment_id'].map(pseudo_reflec)
            
            if self.verbose:
                utils.echo_msg(f"Calculated Pseudo-Reflectance (Strong={is_strong}).")

        except Exception as e:
            utils.echo_warning_msg(f"Pseudo-reflectance calculation failed: {e}")
            
        return df
    
        
    ## ==============================================
    ## Apply ATL* Classifications
    ## ==============================================    
    def apply_atl08_classifications(self, df, atl08_fn, laser, segment_index_dict):
        """Map ATL08 Land/Canopy classes to ATL03 photons. 
        This sets the initial classifications."""
        
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

                ## ==============================================
                ## Calculate Absolute Index in ATL03 Arrays
                ## Global_Index = Segment_Start_Index + (Photon_Index_In_Segment - 1)
                ## ==============================================
                
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

                ## Only consider ground, water surface and open ocean
                _classes = [1, 41, 44]
                
                #mask = is_water & (df['ph_h_classed'] < 40)
                mask = is_water & (df['ph_h_classed'].isin(_classes))
                utils.echo_msg(f'Classed {np.count_nonzero(mask)} photons from {laser} based on ATL13')
                df.loc[mask, 'ph_h_classed'] = 42

        except Exception as e:
            utils.echo_warning_msg(f"Failed to apply ATL13 classifications from {atl13_fn}: {e}")
            
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
                
                if atl24_df.empty: return df

                ## Filter for Bathy Criteria (Class >= 40)
                ## We do this BEFORE merge to keep the merge small
                is_bathy = atl24_df['atl24_class'] == 40
                if self.min_bathy_confidence is not None:
                    is_bathy &= (atl24_df['atl24_conf'] >= self.min_bathy_confidence)
                
                atl24_df = atl24_df[is_bathy]

                if self.verbose:
                    utils.echo_msg(f'Classified {np.count_nonzero(is_bathy)} photons from ATL24')
                
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


    ## ==============================================
    ## Apply Algo Classifications
    ## ==============================================
    def classify_outliers_algo(self, df, multiplier=3.0):
        """Detect and classify outliers using Segment-based IQR.
        
        Outliers (noise/birds/clouds) are reclassified to 0 (Noise).
        This cleans the data before other algorithms run.
        
        Args:
            multiplier (float): IQR threshold. 
                                1.5 is standard. 
                                3.0 is conservative (safest for preserving terrain).
        """
        try:
            if self.verbose:
                utils.echo_msg("Running Dynamic Outlier Classification (IQR)...")

            ## Select Candidates for Check
            valid_classes = [-1,1,2,3]
            candidate_mask = (
                (df['confidence'] >= 3) & 
                (df['ph_h_classed'] != 0)#.isin(valid_classes))
            )
            
            if not np.any(candidate_mask): return df
            
            subset = df[candidate_mask]

            ## Calculate Statistics per Segment (Fast Aggregation)
            grouped = subset.groupby('ph_segment_id')['photon_height']
            
            ## Calculate Q1 (25%) and Q3 (75%)
            ## This returns a Series indexed by ph_segment_id
            q1 = grouped.quantile(0.25)
            q3 = grouped.quantile(0.75)
            
            ## Map Stats back to Photon Level
            ## Map the segment-level Q1/Q3 values to each photon based on its Segment ID
            mapped_q1 = subset['ph_segment_id'].map(q1)
            mapped_q3 = subset['ph_segment_id'].map(q3)
            
            ## Identify Outliers
            iqr = mapped_q3 - mapped_q1
            lower_bound = mapped_q1 - (multiplier * iqr)
            upper_bound = mapped_q3 + (multiplier * iqr)
            
            is_outlier = (subset['photon_height'] < lower_bound) | \
                         (subset['photon_height'] > upper_bound)
            
            outlier_indices = subset.index[is_outlier]
            
            if len(outlier_indices) > 0:
                ## Update Main DataFrame
                ## Set class to 0 (Noise)
                df.loc[outlier_indices, 'ph_h_classed'] = 0
                
                if self.verbose:
                    utils.echo_msg(f"  Classified {len(outlier_indices)} outliers as Noise (0).")

        except Exception as e:
            utils.echo_warning_msg(f"Outlier classification failed: {e}")
            
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

            ## ==============================================
            ## Filter to potential subsurface points to speed up clustering
            ## e.g., only look at points with valid confidence or below water surface
            ## For now, we cluster everything that isn't already Land (1) or Canopy (2/3)
            ## Maybe, if we do this after setting to 44, only use those, since they are
            ## photons over water.
            ## ==============================================
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

    
    def classify_buildings_algo(self, df, min_height=4, max_roughness=0.25, max_range=2.0,
                                max_thickness=0.5, roughness_window=35, ground_window=60,
                                min_reflectance=0.6, dark_veto_threshold=0.25,
                                max_building_length=150):
        """Classify buildings using Geometry (Flat & Elevated) AND Radiometry (Reflectance).
        
        Args:
            df (pd.DataFrame): Photon dataframe.
            min_height (float): Min height above ground (m).
            max_roughness (float): Max StdDev (m) for standard flat roofs.
            ground_window (int): Window for ground estimation.
            roughness_window (int): Window for roughness estimation.
            min_reflectance (float): Reflectance threshold for "Bright" objects (Concrete/Metal).
                                     Vegetation is typically < 0.3. Bright roofs > 0.4.
            max_building_length (float): Max along-track length (m) for a continuous
                                         building segment. Features longer than this
                                         are re-classified as Ground (False Positives).
            max_thickness (float): Max vertical distance between 90th and 10th percentile.
                                   Buildings are thin/solid (<0.5m). Trees are thick (>1.0m).
        """
        
        try:
            if self.verbose:
                utils.echo_msg("Running Dynamic Building Classification (Geometry + Radiometry)...")

            signal_mask = (
                (df['confidence'] >= 3) & 
                (df['ph_h_classed'] != 0) 
            )
            if not np.any(signal_mask): return df
            
            signal_df = df[signal_mask]
                
            ## Estimate Ground & Calculate HAG # min_periods=5
            ground_proxy = df['photon_height'].rolling(
                window=ground_window, center=True, min_periods=ground_window//3
            ).quantile(0.05)
            ground_proxy = ground_proxy.bfill().ffill()
            hag = df['photon_height'] - ground_proxy
            
            ## Identify Elevated Candidates
            is_elevated = (hag >= min_height) & (df['confidence'] >= 3) & (df['ph_h_classed'] != 0)
            
            if not np.any(is_elevated): return df

            ## Calculate Local Roughness (StdDev)
            elevated_df = df[is_elevated].copy()

            roller = elevated_df['photon_height'].rolling(window=roughness_window, center=True, min_periods=5)
            
            elevated_df['local_roughness'] = roller.std()
            elevated_df['local_range'] = roller.max() - roller.min()

            ## Thickness / Penetration Check
            ## We check the spread between the top (90%) and bottom (10%) of the signal
            ## Buildings are "Skin Deep". Forests have "Depth".
            elevated_df['local_thickness'] = roller.quantile(0.90) - roller.quantile(0.10)
            
            ## ==============================================
            ## --- CLASSIFICATION ---            
            ## Geometric (Standard Flat Roofs)
            ## Elevated + Low Roughness (e.g. < 0.6m) + SOLID
            ## This captures most commercial flat roofs (concrete/asphalt)
            ## ==============================================
            mask_geo = (
                (elevated_df['local_roughness'] <= max_roughness) &
                (elevated_df['local_range'] <= max_range) &
                (elevated_df['local_thickness'] <= max_thickness)
            )
            
            ## Filter out Flat Trees
            if 'reflectance' in df.columns:
                ## If it's flat but very dark, it's likely a dense hedge/tree
                is_too_dark = elevated_df['reflectance'] < dark_veto_threshold
                mask_geo = mask_geo & (~is_too_dark)
                 
            ## ==============================================
            ## Radiometric (Bright/Complex Roofs)
            ## Elevated + Modest Roughness (allow up to 1.5m for pitched roofs) + High Reflectance
            ## This captures bright metal/concrete roofs that are too rough
            ## ==============================================
            mask_rad = np.zeros(len(elevated_df), dtype=bool)            
            if 'reflectance' in df.columns:
                ## Check if we actually have valid reflectance data (not all NaNs)
                if elevated_df['reflectance'].notna().any():
                    mask_rad = (
                        (elevated_df['local_roughness'] <= 1.5) &
                        (elevated_df['reflectance'] >= min_reflectance) &
                        (elevated_df['local_range'] <= 3.0) &
                        (elevated_df['local_thickness'] <= 1.5)
                    )
            
            is_building = mask_geo | mask_rad
            building_candidates = elevated_df[is_building].copy()
            
            if len(building_candidates) == 0: return df
            
            ## =========================================================
            ## GROUPING & FILTERING (Length + Wall Check)
            ## =========================================================            
            idx_series = building_candidates.index.to_series()

            ## Identify breaks between potential buildings
            gap_check = idx_series.diff() > 20 # 10 photons ~ 7m gap triggers new group
            group_ids = gap_check.cumsum()
            
            ## Length Check
            max_photon_span = max_building_length / 0.7            
            groups = idx_series.groupby(group_ids)
            group_spans = groups.max() - groups.min()

            ## Wall Check
            ## We calculate the 'Diff' of the full dataset to find jumps
            ## Then we check if any point in the group (or its boundary) is a 'Jump Point'
            full_diffs = df['photon_height'].diff().abs()
            is_wall_jump = full_diffs > min_height  # Threshold: 2 meter vertical jump between neighbors
            
            ## We map this 'Wall Property' to our candidates
            candidate_has_wall = is_wall_jump.loc[building_candidates.index]
            
            ## Check if any point in the group is a wall (or adjacent to one)
            group_has_wall = candidate_has_wall.groupby(group_ids).any()
            
            # Valid = Not Too Long AND Has A Wall
            valid_group_ids = group_spans.index[
                (group_spans <= max_photon_span) & 
                (group_has_wall)
            ]
            
            ## Filter candidates
            final_mask = group_ids.isin(valid_group_ids)
            final_indices = building_candidates.index[final_mask]

            if len(final_indices) > 0:                  
                ## ==============================================
                ## Apply to Main DataFrame
                ## Overwrite Land(1), Canopy(2/3), Unclassified(-1)
                ## Do NOT overwrite Water(42), Bathy(40), Coast(41) or Ocean(44)
                ## Only consider the ground class here (from atl08) and unclassed (-1)
                ## ==============================================
                protected_classes = [40, 41, 42, 44]                
                mask = df.index.isin(final_indices) & (~df['ph_h_classed'].isin(protected_classes))
                df.loc[mask, 'ph_h_classed'] = 7
                
                if self.verbose:
                    count = np.count_nonzero(mask)
                    utils.echo_msg(f"  Classified {count} photons as Buildings.")
                    rejected_len = len(building_candidates) - len(final_indices)
                    if rejected_len > 0:
                        utils.echo_msg(f"    (Rejected {rejected_len} candidates due to Length or Missing Walls)")

        except Exception as e:
            utils.echo_warning_msg(f"Building classification failed: {e}")
            
        return df
        

    def classify_nearshore_roughness(self, df, height_window=2.5, max_roughness=1.5, 
                                     use_reflectance=True, surf_reflectance=0.4):
        """Classify nearshore water, including rough surf zones, using Geometry + Radiometry.
        
        Calm Water: Near Geoid (+-2m) + Flat (StdDev < 0.3).
        Surf Zone: Near Geoid (+-2m) + Rough (StdDev < 1.5) + Bright (Reflectance > 0.4).
        
        Args:
            df (pd.DataFrame): Photon dataframe.
            height_window (float): Max vertical distance from 0m (Geoid).
            max_roughness (float): Max Roughness for Surf.
            use_reflectance (bool): Enable radiometric check for surf.
            surf_reflectance (float): Min reflectance to consider a rough segment as 'Whitewater'.
        """
        
        try:
            if self.verbose:
                utils.echo_msg("Running Dynamic Nearshore Classification (Surf-Aware)...")

            ## Select Signal Photons
            signal_mask = df['confidence'] >= 3
            signal_mask = (
                (df['confidence'] >= 3) & 
                (df['ph_h_classed'] != 0) 
            )
            if not np.any(signal_mask): return df
            
            signal_df = df[signal_mask]

            ## Calculate Segment Statistics
            ## We need Median (Height), Std (Roughness), and Median (Reflectance)
            aggs = {'photon_height': ['median', 'std']}
            if use_reflectance and 'reflectance' in df.columns:
                aggs['reflectance'] = 'median'
                
            grouped = signal_df.groupby('ph_segment_id')
            seg_stats = grouped.agg(aggs)
            
            ## Flatten MultiIndex columns for easier access
            seg_stats.columns = ['_'.join(col).strip() for col in seg_stats.columns.values]
            ## columns: photon_height_median, photon_height_std, reflectance_median (optional)

            ## Identify Water Segments
            
            ## Base Criteria: Elevation near Geoid
            is_near_geoid = seg_stats['photon_height_median'].abs() <= height_window
            
            ## Calm Water (Very Flat, Roughness < 0.3)
            ## Reflectance doesn't matter (can be dark or specular bright)
            is_calm = is_near_geoid & (seg_stats['photon_height_std'] <= 0.3)
            
            ## Rough Surf (Roughness 0.3 - 1.5)
            ## MUST be Bright to distinguish from land/vegetation
            is_surf = np.zeros(len(seg_stats), dtype=bool)
            
            if use_reflectance and 'reflectance_median' in seg_stats.columns:
                is_surf = (
                    is_near_geoid &
                    (seg_stats['photon_height_std'] > 0.3) &
                    (seg_stats['photon_height_std'] <= max_roughness) &
                    (seg_stats['reflectance_median'] >= surf_reflectance)
                )
            elif not use_reflectance:
                ## If no reflectance, we just accept the roughness (riskier)
                is_surf = is_near_geoid & (seg_stats['photon_height_std'] <= max_roughness)

            ## Determine water and Combine
            valid_water_segs = seg_stats.index[is_calm | is_surf]
            
            if len(valid_water_segs) == 0: return df

            is_water_segment = df['ph_segment_id'].isin(valid_water_segs)

            ## Photon ITSELF is within the vertical window
            ## This prevents a 10m tall tree in a water segment from being classed as water
            is_within_window = df['photon_height'].abs() <= height_window
            
            ## Update Class 41 (Coastline/Nearshore)
            #mask = is_water_photon & (df['ph_h_classed'] < 40)
            mask = is_water_segment & is_within_window & (df['ph_h_classed'] < 40)
            df.loc[mask, 'ph_h_classed'] = 41
            
            if self.verbose:
                count = np.count_nonzero(mask)
                utils.echo_msg(f"  Classified {count} photons as Nearshore Water.")
                if np.any(is_surf):
                    surf_count = np.count_nonzero(is_surf)
                    utils.echo_msg(f"    (Captured {surf_count} rough surf segments via reflectance)")

        except Exception as e:
            utils.echo_warning_msg(f"Nearshore classification failed: {e}")
            
        return df

    
    def classify_inland_water_algo(self, df,
                                   max_roughness=0.3,
                                   max_reflectance=0.25,
                                   max_range=0.75,
                                   fill_gaps=True, 
                                   gap_window=15,
                                   fill_threshold=0.3,
                                   debug=False):        
        """Classify inland water by detecting 'Flat & Dark' segments.
        
        Rough & Dark: Vegetation (Excluded)        
        Rough & Bright: Surf / Ice (Included)
        Flat & Dark: Calm Water (Included)
        Flat & Bright: Specular Water / Concrete (Included, but distinguished by height)

        Args:
            df (pd.DataFrame): Photon dataframe.
            max_roughness (float): Max StdDev (m). Inland water is very flat (<0.3m).
            max_reflectance (float): Max Apparent Reflectance. Water is dark (<0.25).
                                     Roads/Salt Flats are flat but bright (>0.4).
            max_range (float): Max Peak-to-Peak height (m) within a segment.
                               Water must be level. Tilted flat surfaces (roads/hills)
                               will have low roughness but high range.
            fill_gaps (bool): If True, coerce 'Land' segments to 'Water' if they 
                              are surrounded by water segments (Majority Filter).
            gap_window (int): Size of the rolling window (segments). 
                              15 segments ~= 300m.
            fill_threshold (float): Percent of water neighbors required to trigger fill.
                                    0.3 = Aggressive (bridges wide gaps).
                                    0.5 = Conservative.
        """
        
        try:
            if self.verbose:
                utils.echo_msg("Running Dynamic Inland Water Classification (Flat & Dark)...")

            if 'reflectance' not in df.columns:
                utils.echo_warning_msg("  Skipping: 'reflectance' data required for inland water algo.")
                return df

            ## Select Signal Photons
            ## We want high confidence and no noise (class 0)
            signal_mask = (
                (df['confidence'] >= 3) & 
                (df['ph_h_classed'] != 0) 
            )
            if not np.any(signal_mask): return df
            
            signal_df = df[signal_mask].copy()

            ## Calculate Segment Statistics
            ## We aggregate by segment to get a robust 'surface' property
            grouped = signal_df.groupby('ph_segment_id')
            
            ## We need Roughness (StdDev of Height) and Brightness (Median Reflectance)
            seg_stats = grouped.agg({
                'photon_height': ['std', 'count', 'max', 'min'],
                'reflectance': 'median'
            })
            seg_stats.columns = ['_'.join(col).strip() for col in seg_stats.columns.values]
            
            ## Calculate Peak-to-Peak (Approximate Slope)
            seg_stats['height_ptp'] = seg_stats['photon_height_max'] - seg_stats['photon_height_min']

            if debug:
                ## --- DEBUG BLOCK ---
                target_id = 249036  # <--- REPLACE THIS with the segment_id you want to track

                if target_id in seg_stats.index:
                    #utils.echo_msg(signal_df[signal_df['ph_segment_id'] == target_id])
                    d_row = seg_stats.loc[target_id]
                    utils.echo_msg(f"\n>>> DEBUG SEGMENT: {target_id} <<<")

                    ## Print Physical Stats
                    ## Note: Adjust column names if using the rolling/algo version (e.g. 'photon_height_std')
                    utils.echo_msg(f"  Roughness (StdDev): {d_row['photon_height_std']:.3f}m  (Limit: {max_roughness})")
                    utils.echo_msg(f"  Reflectance:        {d_row['reflectance_median']:.3f}   (Limit: {max_reflectance})")
                    utils.echo_msg(f"  Slope (PTP):        {d_row['height_ptp']:.3f}m      (Limit: {max_range})")
                    utils.echo_msg(f"  Photon Count:       {d_row['photon_height_count']}")
                    utils.echo_msg(f"  Min Photon Count:   {d_row['photon_height_min']}")
                    utils.echo_msg(f"  Max Photon Count:   {d_row['photon_height_max']}")

                    ## Print Classification Logic Results
                    ## You can manually evaluate the booleans to see exactly which check failed
                    check_dark = (d_row['photon_height_std'] <= max_roughness) and (d_row['reflectance_median'] <= max_reflectance) and (d_row['photon_height_count'] > 3)
                    check_ptp = (d_row['height_ptp'] <= max_range)
                    check_specular = (d_row['reflectance_median'] > 1.5) and (d_row['photon_height_std'] <= max_roughness)
                    check_sparse = (d_row['photon_height_count'] >= 3) and (d_row['photon_height_count'] <= 10) and (d_row['photon_height_std'] <= 0.1)

                    utils.echo_msg(f"  [Passes Dark Water Check?]: {check_dark}")
                    utils.echo_msg(f"  [Passes Slope Check?]:      {check_ptp}")
                    utils.echo_msg(f"  [Passes Specular Check?]:   {check_specular}")
                    utils.echo_msg(f"  [Passes Sparse Check?]:     {check_sparse}")
                    utils.echo_msg("---------------------------------------\n")
            
            ## Dark Water (Standard)
            #is_dark_water = (seg_stats['photon_height_std'] < max_roughness) & (seg_stats['reflectance_median'] < max_reflectance)
            is_dark_water = (
                (seg_stats['photon_height_std'] <= max_roughness) & 
                (seg_stats['reflectance_median'] <= max_reflectance) &
                (seg_stats['photon_height_count'] > 3)
                #(seg_stats['height_ptp'] <= max_range)  # Must be Level
            )
            
            ## Specular Water (Mirror (water/ice/surf)
            ## Very Flat + Very Bright (Saturation)
            is_specular_water = (
                (seg_stats['reflectance_median'] > 1.5) &  # Very Bright threshold (high count implied)
                (seg_stats['photon_height_std'] <= max_roughness) #&   # STRICT Roughness
                #(seg_stats['height_ptp'] <= .3)  # Must be Level
            )

            ## Sparse Water
            ## Extremely Flat, Low Density
            is_sparse_water = (
                (seg_stats['photon_height_count'] >= 3) &   # At least 3 points to form a line
                (seg_stats['photon_height_count'] <= 10) &  # Low density (Weak beam char)
                (seg_stats['photon_height_std'] <= 0.1) #&    # Extremely coherent/flat
                #(seg_stats['reflectance_median'] > 1.5)  # Very Bright threshold
                #(seg_stats['height_ptp'] <= 0.3)        # Must be Level
            )

            seg_stats['is_water'] = (is_dark_water | is_specular_water | is_sparse_water).astype(int)

            ## Gap Filling (Majority Filter)
            if fill_gaps:
                ## Sort by index (Segment ID) to ensure spatial continuity
                seg_stats.sort_index(inplace=True)

                ## Calculate percent of water neighbors
                neighbor_water_rate = seg_stats['is_water'].rolling(
                    window=gap_window, center=True, min_periods=1
                ).mean()
                
                ## Check for "Islands" (Valid Land)
                ## If a segment is a cliff (PTP > 2.0) or super rough (Std > 1.0), it is REAL land.
                ## Don't overwrite it just because it's in a lake.
                #is_valid_land = (seg_stats['height_ptp'] > 2.0) | (seg_stats['photon_height_std'] > 1.0)

                ## This lets us fill over choppy waves or reeds that looked like land.
                ## Real islands (cliffs/trees) usually exceed 5m+ slope/roughness.
                is_safe_to_overwrite = (
                    (seg_stats['height_ptp'] <= 4.0) &       # Increased from 2.0
                    (seg_stats['photon_height_std'] <= 2.5)  # Increased from 1.0
                )
                
                ## Identify Gaps to Fill
                is_gap_fill = (
                    (neighbor_water_rate > fill_threshold) &
                    (is_safe_to_overwrite) &
                    (seg_stats['is_water'] == 0)
                )
                
                ## Apply Fill
                seg_stats.loc[is_gap_fill, 'is_water'] = 1
                
                if self.verbose:
                    filled_count = np.count_nonzero(is_gap_fill)
                    if filled_count > 0:
                        utils.echo_msg(f"    Filled {filled_count} segment gaps in water bodies.")

            valid_water_segs = seg_stats.index[seg_stats['is_water'] == 1]
            if len(valid_water_segs) == 0: return df

            ## Apply to Main DataFrame
            is_water_photon = df['ph_segment_id'].isin(valid_water_segs)
            
            ## ==============================================
            ## Update Class 42 (Inland Water)
            ## Do not overwrite Bathy(40) or Ocean(44) if already set
            ## We overwrite Land(1), Unclassified(-1), etc.
            ## ==============================================
            mask = is_water_photon & (df['ph_h_classed'] < 40)
            df.loc[mask, 'ph_h_classed'] = 42
            
            if self.verbose:
                utils.echo_msg(f"  Classified {np.count_nonzero(mask)} photons as Inland Water (42).")
                if np.any(is_specular_water):
                    utils.echo_msg(f"    (Rescued {np.count_nonzero(is_specular_water)} bright/flat segments)")
                if np.any(is_sparse_water):
                     utils.echo_msg(f"    (Rescued {np.count_nonzero(is_sparse_water)} sparse segments via geometric coherence)")


        except Exception as e:
            utils.echo_warning_msg(f"Inland water classification failed: {e}")
            
        return df    
        
    
    ## ==============================================
    ## Classify by Masks
    ## ==============================================
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

    
    ## ==============================================
    ## ATL03 - Read and Classify
    ## ==============================================    
    def read_atl03(self, f, laser_num, orientation=None, atl08_fn=None, atl09_fn=None, atl24_fn=None, atl06_fn=None, atl12_fn=None, atl13_fn=None):
        """Read and classify data from an ATL03 file/laser."""
        
        if orientation is None:
            orientation = f['/orbit_info/sc_orient'][0]
            
        laser = 'gt' + laser_num + self.orientDict[orientation]
        if laser not in f or 'heights' not in f[laser]: return None

        try:
            # sc_orient: 0 = Backward, 1 = Forward
            true_sc_orient = f['/orbit_info/sc_orient'][0]
        except KeyError:
            true_sc_orient = 0 # Default

        side = laser[-1]
        #is_strong = self.is_strong_beam(laser, orientation)
        # Backward (0) -> Left is Strong
        # Forward (1)  -> Right is Strong
        is_strong = (true_sc_orient == 0 and side == 'l') or \
                    (true_sc_orient == 1 and side == 'r')
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

        ## ==============================================
        ## --- Indexing ---
        ## mrl: move these down for merging.
        #seg_starts = np.concatenate(([0], np.cumsum(seg_ph_cnt)[:-1]))
        #seg_idx_dict = dict(zip(seg_id, seg_starts))
        ## ==============================================
        
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

        ## ==============================================
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
        ## ==============================================
        
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

        ## ==============================================
        ## --- Calculate Heights ---
        ## ==============================================
        h_geoid_map = dict(zip(seg_id, geoid))
        h_f2m_map = dict(zip(seg_id, geoid_f2m))
        h_dem_map = dict(zip(seg_id, dem_h))

        p_geoid = np.array([h_geoid_map.get(s, 0) for s in ph_seg_ids])
        p_f2m = np.array([h_f2m_map.get(s, 0) for s in ph_seg_ids])
        p_dem = np.array([h_dem_map.get(s, 0) for s in ph_seg_ids])

        h_ortho = h_ph - p_geoid
        h_meantide = h_ph - (p_geoid + p_f2m)
        h_dem = p_dem - (p_geoid + p_f2m)

        ## ==============================================
        ## --- Initialize Classifications ---
        ## Start with -1 (Unclassified)
        ## ==============================================
        ph_class = np.full(len(h_ph), -1, dtype=int)        

        ## ==============================================
        ## --- Create DataFrame ---
        ## ==============================================
        
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

        # ## Filter by Confidence
        # if self.confidence_levels:
        #     df = df[df['confidence'].isin(self.confidence_levels)]
        #     if df.empty: return None
        
        seg_starts = np.concatenate(([0], np.cumsum(seg_ph_cnt)[:-1]))
        seg_idx_dict = dict(zip(seg_id, seg_starts))

        ## ==============================================
        ## --- Apply Classifications based on ATLXX data ---
        ## ==============================================
        
        ## ATL08,Land/Canopy, Ground (1), Canopy (2), Top Canopy (3)
        if atl08_fn:
            df = self.apply_atl08_classifications(df, atl08_fn, laser, seg_idx_dict)

        ## ATL09 to get reflectance for classifications
        if atl09_fn:
            df = self.apply_atl09_data(df, atl09_fn, laser)

        ## Fallback to 'pseudo-reflectance' if ATL09 isn't available.
        #is_strong = self.is_strong_beam(laser, orientation)
        if 'reflectance' not in df.columns or df['reflectance'].isna().all():
            df = self.calculate_pseudo_reflectance(df, is_strong=is_strong)

        ## ==============================================
        ## --- Apply Fallback Ocean Classification --- (44 - Ocean)
        ## Onboard Mask says Ocean (ph_is_ocean == 1)
        ## Height is near Sea Level (abs(h_ortho) < threshold). 
        ## Open ocean should be near 0. We can use 30m to be very safe against tides/surge/geoid error.
        ## We use 2m to be less safe, to avoid classifying too much near-shore land.
        ## ==============================================
        is_open_ocean = (ph_is_ocean == 1) & (np.abs(h_ortho) < 2)
        if self.verbose:
            utils.echo_msg(f'Classified {np.count_nonzero(is_open_ocean)} photons to open ocean')            
        df.loc[is_open_ocean, 'ph_h_classed'] = 44

        ## This cleans spikes so they don't ruin downstream stats (Roughness/Slope)
        ## We use a safe 3.0 multiplier to avoid killing steep terrain features.
        df = self.classify_outliers_algo(df, multiplier=3.0)
        #df = self.classify_outliers_algo(df, multiplier=1.5)
        
        ## ATL24,Bathymetry, Near-shore Bathymetry (40), Coastline (41)
        ## This will over-ride atl08/open-ocean classes.
        if atl24_fn:
            df = self.apply_atl24_classifications(
                df, atl24_fn, laser, geoseg_beg, geoseg_end
            )
            
        ## ATL13,Inland Water, (42)
        if atl13_fn:
            df = self.apply_atl13_classifications(df, atl13_fn, laser)

        ## ATL12,Ocean,Open Ocean Surface, (44)
        if atl12_fn:
            df = self.apply_atl12_classifications(df, atl12_fn, laser)

        # ## ATL06,Land Ice, Glaciers, Ice Sheets, (6)
        # if atl06_fn:
        #     df = self.apply_atl06_classifications(df, atl06_fn, laser)

        ## ==============================================
        ## --- Apply Dynamic Nearshore Classification ---
        ## This fills gaps between ATL12 (Ocean) and ATL08 (Land)
        ## It is robust in the surf zone where ATL12 often drops out.
        ## Use wider window (3.5m) and higher roughness (1.5m) for surf zones
        ## ==============================================
        df = self.classify_nearshore_roughness(df)#, height_window=3.5, max_roughness=1.5)

        ## ==============================================
        ## --- Apply Dynamic Inland Water Classification ---
        ## This checks for roughness (want flat areas) and also checks for
        ## reflectance (inland water is dark compared to other flat features).
        ## This requires ATL09 to get reflectance.
        ## ==============================================
        #df = self.classify_inland_water_rolling(df, window=40, max_roughness=0.4, max_range=0.75)#, max_roughness=.45, max_reflectance=.2)
        df = self.classify_inland_water_algo(df, max_roughness=0.45, max_reflectance=0.2, max_range=1, fill_gaps=True)

        ## ==============================================
        ## --- Apply Dynamic Building Classification ---
        ## This classifies buildings based on the height above ground, the roughness of those heights
        ## and the reflectance (from ATL09) of the surface.
        ## ==============================================
        df = self.classify_buildings_algo(df)#, min_height=5, max_roughness=.6, min_reflectance=.8, roughness_window=60)

        ## ==============================================
        ## --- Apply Algorithmic Bathymetry Classification (Known Bathy / DBSCAN) ---
        ## ==============================================
        if self.known_bathymetry or self.use_dbscan:
            df = self.classify_bathymetry_algo(df)

        return df


    ## ==============================================
    ## Yield (classified) points
    ## This goes through the entire fetchings/classifying process
    ## ==============================================    
    def yield_points(self):
        """Pipeline to yield classified points."""

        ## Initialize external mask geometries
        bing_geom = None
        osm_geom = None
        osm_lakes = None

        ## ==============================================    
        ## --- External Masks ---
        ## BingBFP for buildings
        ## OSM for land/water mask (including inland water bodies)
        ##
        ## Fallback for missing/misclassified ATL* data.
        ## This ensures we get 'bare-earth' when needed...
        ## ==============================================    
        if self.use_external_masks:        
            ## Setup Region for Aux Mask Fetches
            ## We try to use the user-supplied region or the INF metadata bounds
            proc_region = self._get_processing_region()
            if proc_region is not None:
                proc_region.buffer(pct=1) # Slight buffer for safety        

                if self.classes:
                    bfp = bingbfp.BingBuildings(region=proc_region, verbose=self.verbose, cache_dir=self.cache_dir)
                    bing_geom = bfp(return_geom=True)

                if self.classes:
                    osmc = osm.osmCoastline(region=proc_region, q='coastline', verbose=self.verbose, cache_dir=self.cache_dir)
                    osm_geom = osmc(return_geom=True)

                if self.classes:
                    osml = osm.osmCoastline(region=proc_region, q='water', verbose=self.verbose, cache_dir=self.cache_dir)
                    osm_lakes = osml(return_geom=True)

        ## Fetch Aux ATLXX Data
        atl08_fn = self.fetch_atlxx(self.fn, 'ATL08') if self.classes else None
        atl24_fn = self.fetch_atlxx(self.fn, 'ATL24') if self.classes else None
        # atl06_fn = self.fetch_atlxx(self.fn, 'ATL06') if self.classes else None
        atl12_fn = self.fetch_atlxx(self.fn, 'ATL12') if self.classes else None
        atl13_fn = self.fetch_atlxx(self.fn, 'ATL13') if self.classes else None
        atl09_fn = self.fetch_atlxx(self.fn, 'ATL09') if self.classes else None
        atl06_fn = None
        
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
                
            utils.echo_msg('Loading ATL data into data-frame...')
            for i in range(1, 4):
                for orient in range(2):
                    ## ==============================================    
                    ## Read_atl03 reads all atl03 data and classifies the
                    ## photons (in 'ph_h_classed') based on aux ATL data.
                    ## This creates the pandas dataframe with classifications
                    ## from aux ATL* files and algos.
                    ## ==============================================    
                    utils.echo_msg(f'  Loading laser: {i}{self.orientDict[orient]}')
                    dataset = self.read_atl03(
                        f, str(i), orientation=orient,
                        atl08_fn=atl08_fn, atl09_fn=atl09_fn, atl24_fn=atl24_fn, atl06_fn=atl06_fn
                    )

                    if dataset is None or dataset.empty: continue

                    ## Filter by confidence
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

                    ## ==============================================    
                    ## We classify buildings (7), ocean (41) and inland water (42) - maybe 44 for ocean?
                    ## using Bing BFP and OpenStreetMap. These will over-ride (for the most part)
                    ## previously classed data.
                    ##
                    ## These get skipped unless `use_external_masks` is set to True.
                    ## ==============================================    
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
                

### End
