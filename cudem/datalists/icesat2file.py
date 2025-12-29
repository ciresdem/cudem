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
## ATL09,Atmosphere,Cloud Layers (Filtering Noise),0 (Noise)
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

            
    def generate_inf(self, make_grid=True, make_block_mean=False, block_inc=None):
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

                    
class IceSat2_ATL03(ElevationDataset):
    """ICESat-2 ATL03 (Global Geolocated Photon Data) Parser.
    
    Can integrate classifications from ATL08 (Land/Canopy), ATL24 (Bathymetry),
    and ATL06 (Land Ice), and ATL13 (Inland Water). 
    Can also classify using external vectors (OSM/Bing).
    """

    def __init__(self,
                 water_surface='geoid',
                 classes=None,
                 confidence_levels='2/3/4',
                 columns=None,
                 reject_failed_qa=True,
                 classify_buildings=True,
                 classify_water=True,
                 classify_inland_water=True,
                 append_atl24=False,
                 min_bathy_confidence=None,
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

        self.orientDict = {0:'l', 1:'r', 21:'error'}

        
    def generate_inf(self, make_grid=True, make_block_mean=False, block_inc=None):
        """Generate metadata (INF) for the ATL03 data."""
        
        ## Delegate to parent for grids
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
            
        except Exception:
            ## Fallback to scan
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

    
    def apply_atl08_classifications(self, atl08_fn, ph_h_classed, laser, segment_id, segment_index_dict, ph_segment_ids):
        """Map ATL08 Land/Canopy classes to ATL03 photons."""
        
        try:
            with h5.File(atl08_fn, 'r') as f:
                if laser not in f: return ph_h_classed
                
                sig = f[f'/{laser}/signal_photons']
                atl08_flag = sig['classed_pc_flag'][...]
                atl08_seg = sig['ph_segment_id'][...]
                atl08_idx = sig['classed_pc_indx'][...]

                ## Mask for segments present in current chunk
                mask = np.isin(atl08_seg, ph_segment_ids)
                
                ## Map to ATL03 indices
                ## ATL03 Global Index = Segment_Start_Index + (ATL08_Offset - 1)
                
                ## Get start index for each ATL08 segment in the chunk
                atl08_seg_starts = np.array([segment_index_dict[sid] for sid in atl08_seg[mask]])
                
                ## Calculate absolute photon indices
                atl03_indices = atl08_seg_starts + (atl08_idx[mask] - 1)
                
                ## Boundary check
                valid_idx = atl03_indices < len(ph_h_classed)
                
                ## Apply Classification
                ph_h_classed[atl03_indices[valid_idx]] = atl08_flag[mask][valid_idx]
                
        except Exception as e:
            utils.echo_warning_msg(f"Failed to apply ATL08 classifications from {atl08_fn}: {e}")
            
        return ph_h_classed

    
    def apply_atl06_classifications(self, atl06_fn, ph_h_classed, laser, segment_id, segment_index_dict, ph_segment_ids):
        """Map ATL06 Land Ice segments to ATL03 photons (Class 6)."""
        
        try:
            with h5.File(atl06_fn, 'r') as f:
                if laser not in f: return ph_h_classed
                
                li = f[f'/{laser}/land_ice_segments']
                if 'segment_id' not in li: return ph_h_classed

                atl06_seg = li['segment_id'][...]
                atl06_h = li['h_li'][...]
                
                ## Filter valid ATL06 segments (valid height)
                valid_mask = atl06_h < 1e30
                valid_segs = atl06_seg[valid_mask]

                ## Identify photons belonging to these valid segments
                is_ice = np.isin(ph_segment_ids, valid_segs)
                
                ## Update Classification (Class 6 = Ice)
                ## Override unclassified/noise, preserve Bathy/Urban if already set
                update_mask = is_ice & (ph_h_classed < 7)
                ph_h_classed[update_mask] = 6

        except Exception as e:
            utils.echo_warning_msg(f"Failed to apply ATL06 classifications from {atl06_fn}: {e}")
            
        return ph_h_classed


    def apply_atl12_classifications(self, atl12_fn, ph_h_classed, laser, segment_id, segment_index_dict, ph_segment_ids):
        """Map ATL12 Ocean Surface segments to ATL03 photons (Class 44).
        
        ATL12 provides ocean surface heights.
        Class 44 is assigned to photons belonging to valid ATL12 ocean segments.
        """
        
        try:
            with h5.File(atl12_fn, 'r') as f:
                if laser not in f: return ph_h_classed
                
                ## Check for segment arrays
                ## ATL12 structure: /gtx/ssh_segments/stats/segment_id_beg
                ## (Note: ATL12 groups by ssh_segments, similar to ATL13)
                
                path = f'/{laser}/ssh_segments/stats'
                if path not in f or 'segment_id_beg' not in f[path]: return ph_h_classed
                
                atl12_seg = f[f'{path}/segment_id_beg'][...]
                
                ## Filter by valid height if needed (h_ssh)
                ## h_ssh = f[f'/{laser}/ssh_segments/heights/h_ssh'][...]
                
                ## Assume presence in ATL12 means it is ocean
                is_ocean = np.isin(ph_segment_ids, atl12_seg)
                
                ## Update Mask: Is Ocean Segment AND Not already special class (like Bathy 40/41)
                ## We prioritize Bathy (40) and Inland Water (42) over Ocean (44) if there's overlap.
                ## Usually: Bathy > Inland > Ocean > Unclassified
                update_mask = is_ocean & (ph_h_classed < 40)
                ph_h_classed[update_mask] = 44

        except Exception as e:
            utils.echo_warning_msg(f"Failed to apply ATL12 classifications from {atl12_fn}: {e}")
            
        return ph_h_classed
    
    
    def apply_atl13_classifications(self, atl13_fn, ph_h_classed, laser, segment_id, segment_index_dict, ph_segment_ids):
        """Map ATL13 Inland Water segments to ATL03 photons (Class 42).
        
        ATL13 provides water surface heights at the segment level (short or long segments).
        We map these segments back to the individual photons.
        Class 42 is assigned to photons belonging to valid ATL13 water segments.
        """
        
        try:
            with h5.File(atl13_fn, 'r') as f:
                if laser not in f: return ph_h_classed
                
                if f'/{laser}/segment_id_beg' not in f: return ph_h_classed
                
                atl13_seg = f[f'/{laser}/segment_id_beg'][...]
                #atl13_h = f[f'/{laser}/ht_water_surf'][...] # could filter by valid height
                
                ## Assume presence in ATL13 means it is water
                is_water = np.isin(ph_segment_ids, atl13_seg)
                
                ## Update Mask: Is Water Segment AND Not already special class (like Bathy 40/41)
                update_mask = is_water & (ph_h_classed < 40)
                ph_h_classed[update_mask] = 42

        except Exception as e:
            utils.echo_warning_msg(f"Failed to apply ATL13 classifications from {atl13_fn}: {e}")
            
        return ph_h_classed

    
    def apply_atl24_classifications(self, atl24_fn, ph_h_classed, ph_h_bathy_conf, 
                                    lat, lon, h, h_mean, h_geoid, 
                                    laser, geoseg_beg, geoseg_end, ph_segment_ids):
        """Map ATL24 Bathymetry classes and Refracted Heights to ATL03 photons."""
        
        try:
            with h5.File(atl24_fn, 'r') as f:
                if laser not in f: return ph_h_classed, ph_h_bathy_conf, lat, lon, h, h_mean, h_geoid

                grp = f[laser]
                atl24_lat = grp['lat_ph'][...]
                atl24_lon = grp['lon_ph'][...]
                atl24_ortho = grp['ortho_h'][...]
                atl24_ellipse = grp['ellipse_h'][...]
                atl24_surface = grp['surface_h'][...]
                atl24_class = grp['class_ph'][...]
                atl24_seg = grp['index_seg'][...] # relative segment index? ATL24 uses index_seg differently
                atl24_idx = grp['index_ph'][...]
                atl24_conf = grp['confidence'][...]

                ## Reconstruct original segment IDs for mapping
                ## ATL24 index_seg refers to the range [geoseg_beg, geoseg_end]
                orig_segs = np.arange(geoseg_beg, geoseg_end + 1)
                
                ## Filter ATL24 points relevant to current chunk's segments
                ## Map ATL24 internal segment index -> Real Segment ID
                atl24_real_seg_ids = orig_segs[atl24_seg]
                
                mask = np.isin(atl24_real_seg_ids, ph_segment_ids)
                if not np.any(mask): return ph_h_classed, ph_h_bathy_conf, lat, lon, h, h_mean, h_geoid
                
                chunk_min_idx = np.min(atl24_idx[mask]) # Simplified alignment assumption
                rel_idx = atl24_idx[mask] - chunk_min_idx
                
                ## Safety Clip
                valid_rel = (rel_idx >= 0) & (rel_idx < len(ph_h_classed))
                
                ## Filter for Bathy (Class >= 40)
                is_bathy = atl24_class[mask][valid_rel] >= 40
                if self.min_bathy_confidence is not None:
                    is_bathy &= (atl24_conf[mask][valid_rel] >= self.min_bathy_confidence)

                target_indices = rel_idx[valid_rel][is_bathy]
                
                ## Update Arrays
                ph_h_classed[target_indices] = atl24_class[mask][valid_rel][is_bathy]
                ph_h_bathy_conf[target_indices] = atl24_conf[mask][valid_rel][is_bathy]
                
                ## Update Coordinates (Refracted)
                lat[target_indices] = atl24_lat[mask][valid_rel][is_bathy]
                lon[target_indices] = atl24_lon[mask][valid_rel][is_bathy]
                h[target_indices] = atl24_ellipse[mask][valid_rel][is_bathy]
                h_mean[target_indices] = atl24_surface[mask][valid_rel][is_bathy]
                h_geoid[target_indices] = atl24_ortho[mask][valid_rel][is_bathy]

        except Exception as e:
            utils.echo_warning_msg(f"Failed to apply ATL24 data from {atl24_fn}: {e}")

        return ph_h_classed, ph_h_bathy_conf, lat, lon, h, h_mean, h_geoid

    def read_atl03(self, f, laser_num, orientation=None, atl08_fn=None, atl24_fn=None, atl06_fn=None, atl12_fn=None, atl13_fn=None):
        """Read and classify data from an ATL03 file/laser."""
        
        if orientation is None:
            orientation = f['/orbit_info/sc_orient'][0]
            
        laser = 'gt' + laser_num + self.orientDict[orientation]
        if laser not in f or 'heights' not in f[laser]: return None

        ## --- Read Data ---
        try:
            h_grp = f[f'/{laser}/heights']
            geo_grp = f[f'/{laser}/geolocation']
            geophys_grp = f[f'/{laser}/geophys_corr']
            
            # Photon Data
            lat = h_grp['lat_ph'][...]
            lon = h_grp['lon_ph'][...]
            h_ph = h_grp['h_ph'][...]
            conf = h_grp['signal_conf_ph'][..., 0] # Column 0 = Land
            dt = h_grp['delta_time'][...]
            dist_ph = h_grp['dist_ph_along'][...]

            ## Geolocation Data
            seg_ph_cnt = geo_grp['segment_ph_cnt'][...]
            seg_id = geo_grp['segment_id'][...]
            seg_dist_x = geo_grp['segment_dist_x'][...]
            
            ## Corrections
            geoid = geophys_grp['geoid'][...]
            geoid_f2m = geophys_grp['geoid_free2mean'][...]
            dem_h = geophys_grp['dem_h'][...]

            ## Ancillary
            anc = f['ancillary_data']
            geoseg_beg = anc['start_geoseg'][0]
            geoseg_end = anc['end_geoseg'][0]

        except KeyError:
            return None

        ## --- Indexing ---
        ## Map Segment ID -> Start Index in Photon Array
        ## Cumulative sum of counts gives start indices
        seg_starts = np.concatenate(([0], np.cumsum(seg_ph_cnt)[:-1]))
        seg_idx_dict = dict(zip(seg_id, seg_starts))
        
        ## Map Photon -> Segment ID
        ## np.searchsorted to find which segment bin each photon belongs to
        ## Creating array of segment IDs for every photon:
        ## Repeat segment_id[i] count[i] times
        ph_seg_ids = np.repeat(seg_id, seg_ph_cnt)
        
        ## Truncate if mismatch (sometimes happens with last partial segment)
        min_len = min(len(ph_seg_ids), len(h_ph))
        ph_seg_ids = ph_seg_ids[:min_len]
        lat = lat[:min_len]; lon = lon[:min_len]; h_ph = h_ph[:min_len]
        conf = conf[:min_len]; dt = dt[:min_len]

        ## --- Calculate Heights ---
        ## Map segment corrections to photons
        h_geoid_map = dict(zip(seg_id, geoid))
        h_f2m_map = dict(zip(seg_id, geoid_f2m))
        h_dem_map = dict(zip(seg_id, dem_h))

        ## Vectorized map using search/indexing is faster than list(map lambda)
        ## But dictionary lookups are robust for sparse segments.
        p_geoid = np.array([h_geoid_map.get(s, 0) for s in ph_seg_ids])
        p_f2m = np.array([h_f2m_map.get(s, 0) for s in ph_seg_ids])
        p_dem = np.array([h_dem_map.get(s, 0) for s in ph_seg_ids])

        h_ortho = h_ph - p_geoid
        h_meantide = h_ph - (p_geoid + p_f2m)
        h_dem = p_dem - (p_geoid + p_f2m)

        ## Initialize Classification (-1 = Unclassified)
        ph_class = np.full(h_ph.shape, -1, dtype=int)
        ph_bathy_conf = np.full(h_ph.shape, -1, dtype=int)

        ## --- Apply Classifications ---
        if atl08_fn:
            ph_class = self.apply_atl08_classifications(
                atl08_fn, ph_class, laser, seg_id, seg_idx_dict, ph_seg_ids
            )
            
        if atl06_fn:
            ph_class = self.apply_atl06_classifications(
                atl06_fn, ph_class, laser, seg_id, seg_idx_dict, ph_seg_ids
            )

        if atl12_fn:
            ph_class = self.apply_atl12_classifications(
                atl12_fn, ph_class, laser, seg_id, seg_idx_dict, ph_seg_ids
            )
            
        if atl13_fn:
            ph_class = self.apply_atl13_classifications(
                atl13_fn, ph_class, laser, seg_id, seg_idx_dict, ph_seg_ids
            )
            
        if atl24_fn:
            ph_class, ph_bathy_conf, lat, lon, h_ph, h_meantide, h_ortho = self.apply_atl24_classifications(
                atl24_fn, ph_class, ph_bathy_conf, lat, lon, h_ph, h_meantide, h_ortho,
                laser, geoseg_beg, geoseg_end, ph_seg_ids
            )

        ## Select Output Height
        if self.water_surface == 'mean_tide': z_out = h_meantide
        elif self.water_surface == 'geoid': z_out = h_ortho
        else: z_out = h_ph

        ## Build DataFrame
        df = pd.DataFrame({
            'latitude': lat, 'longitude': lon, 'photon_height': z_out,
            'laser': laser, 'fn': self.fn, 'confidence': conf, 'delta_time': dt,
            'bathy_confidence': ph_bathy_conf, 'photon_h_dem': h_dem, 'ph_h_classed': ph_class
        })
        
        ## Extra columns
        for col_path, col_name in self.columns.items():
            try:
                ## Handle simple path
                if col_path in f:
                    val = f[col_path][...]
                    pass 
            except: pass

        return df

    
    def yield_points(self):
        """Pipeline to yield classified points."""
        
        ## Fetch Aux Data
        atl08_fn = self.fetch_atlxx(self.fn, 'ATL08') if self.classes else None
        atl24_fn = self.fetch_atlxx(self.fn, 'ATL24') if self.classes else None
        atl06_fn = self.fetch_atlxx(self.fn, 'ATL06') if self.classes else None
        atl12_fn = self.fetch_atlxx(self.fn, 'ATL12') if self.classes else None
        atl13_fn = self.fetch_atlxx(self.fn, 'ATL13') if self.classes else None

        ## External Masks (Buildings/Water)
        bing_geom = None
        if self.classify_buildings:
            bfp = bingbfp.BingBuildings(region=self.region, verbose=self.verbose, cache_dir=self.cache_dir)
            bing_geom = bfp(return_geom=True)

        osm_geom = None
        if self.classify_water:
            osmc = osm.osmCoastline(region=self.region, verbose=self.verbose, cache_dir=self.cache_dir)
            osm_geom = osmc(return_geom=True)
            
        osm_lakes = None
        if self.classify_inland_water:
            osml = osm.osmCoastline(region=self.region, q='water', verbose=self.verbose, cache_dir=self.cache_dir)
            osm_lakes = osml(return_geom=True)

        ## Process Granule
        with h5.File(self.fn, 'r') as f:
            ## QA Check
            if self.reject_failed_qa and 'quality_assessment' in f:
                if f['/quality_assessment/qa_granule_pass_fail'][0] != 0:
                    utils.echo_warning_msg(f"Skipping failed granule {self.fn}")
                    return

            for i in range(1, 4):
                for orient in range(2):
                    dataset = self.read_atl03(
                        f, str(i), orientation=orient,
                        atl08_fn=atl08_fn, atl24_fn=atl24_fn, atl06_fn=atl06_fn
                    )
                    
                    if dataset is None or dataset.empty: continue

                    ## Filter by Confidence
                    if self.confidence_levels:
                        dataset = dataset[dataset['confidence'].isin(self.confidence_levels)]

                    ## Spatial Subset (Coarse)
                    if self.region:
                        dataset = dataset[
                            (dataset['longitude'] >= self.region.xmin) & 
                            (dataset['longitude'] <= self.region.xmax) & 
                            (dataset['latitude'] >= self.region.ymin) & 
                            (dataset['latitude'] <= self.region.ymax)
                        ]

                    if dataset.empty: continue

                    ## Spatial Classification (Masks)
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
        
        ## Optimizing: Iterate masks, set filter, update dataframe indices found
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
                dataset.loc[mask, 'ph_h_classed'] = classification
                
            layer.SetSpatialFilter(None)
            
        return dataset

### End
