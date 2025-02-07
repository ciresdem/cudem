#!/usr/bin/env python
### icesat2dump.py
##
## Copyright (c) 2024, 2025 Matthew Love <matthew.love@colorado.edu>
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
### Code:

import os
import sys
import h5py as h5
import numpy as np
import pandas as pd
from cudem import utils
from cudem import fetches
from cudem import xyzfun

def read_atl_data(f, laser_num, orientation, f08 = None, f24 = None, water_surface = 'geoid'):
    """Read data from an ATL03 file

    Adapted from 'cshelph' https://github.com/nmt28/C-SHELPh.git 
    and 'iVert' https://github.com/ciresdem/ivert.git

    laser_num is 1, 2 or 3
    surface is 'mean_tide', 'geoid' or 'ellipsoid'
    """

    #orientation = f['/orbit_info/sc_orient'][0]

    ## selects the strong beams only [we can include weak beams later on]
    orientDict = {0:'l', 1:'r', 21:'error'}
    laser = 'gt' + laser_num + orientDict[orientation]

    ## for 'subsets', where heights don't come through
    if 'heights' not in f['/{}'.format(laser)].keys():
        return(None)

    ## Read in the required photon level data
    photon_h = f['/' + laser + '/heights/h_ph'][...,]
    latitude = f['/' + laser + '/heights/lat_ph'][...,]
    longitude = f['/' + laser + '/heights/lon_ph'][...,]
    ph_count = f['/' + laser + '/heights/ph_id_count'][...,]
    conf = f['/' + laser + '/heights/signal_conf_ph/'][...,0]
    qual = f['/' + laser + '/heights/quality_ph/'][...,0]
    dist_ph_along = f['/' + laser + '/heights/dist_ph_along'][...,]
    this_N = latitude.shape[0]

    ## Read in the geolocation level data
    segment_ph_cnt = f['/' + laser + '/geolocation/segment_ph_cnt'][...,]
    segment_id = f['/' + laser + '/geolocation/segment_id'][...,]
    segment_dist_x = f['/' + laser + '/geolocation/segment_dist_x'][...,]
    ref_elev = f['/' + laser + '/geolocation/ref_elev'][...,]
    ref_azimuth = f['/' + laser + '/geolocation/ref_azimuth'][...,]
    ref_photon_index = f['/' + laser + '/geolocation/reference_photon_index'][...,]
    ph_index_beg = f['/' + laser + '/geolocation/ph_index_beg'][...,]
    altitude_sc = f['/' + laser + '/geolocation/altitude_sc'][...,]

    ## Read in the geoid data
    photon_geoid = f['/' + laser + '/geophys_corr/geoid'][...,]
    photon_geoid_f2m = f['/' + laser + '/geophys_corr/geoid_free2mean'][...,]

    ## Create a dictionary with (segment_id --> index into ATL03 photons)
    ## lookup pairs, for the starting photon of each segment
    segment_indices = np.concatenate(([0], np.cumsum(segment_ph_cnt)[:-1]))
    segment_index_dict = dict(zip(segment_id, segment_indices))
    ph_segment_ids = segment_id[np.searchsorted(segment_indices, np.arange(0.5, len(photon_h), 1))-1]

    ## meantide/geoid heights
    h_geoid_dict = dict(zip(segment_id, photon_geoid))
    ph_h_geoid = np.array(list(map((lambda pid: h_geoid_dict[pid]), ph_segment_ids)))        
    h_meantide_dict = dict(zip(segment_id, photon_geoid_f2m))
    ph_h_meantide = np.array(list(map((lambda pid: h_meantide_dict[pid]), ph_segment_ids)))
    photon_h_geoid = photon_h - ph_h_geoid
    photon_h_meantide = photon_h - (ph_h_geoid + ph_h_meantide)

    ## setup classifications
    ph_h_classed = np.zeros(photon_h.shape)
    ph_h_classed[:] = -1

    ## append the laser to each record
    laser_arr = np.empty(photon_h.shape, dtype='object')
    laser_arr[:] = laser

    ## ref values
    h_ref_elev_dict = dict(zip(segment_id, ref_elev))
    ph_ref_elev = np.array(list(map((lambda pid: h_ref_elev_dict[pid]), ph_segment_ids)))#.astype(float)        
    h_ref_azimuth_dict = dict(zip(segment_id, ref_azimuth))
    ph_ref_azimuth = np.array(list(map((lambda pid: h_ref_azimuth_dict[pid]), ph_segment_ids)))#.astype(float)
    h_altitude_sc_dict = dict(zip(segment_id, altitude_sc))
    ph_altitude_sc = np.array(list(map((lambda pid: h_altitude_sc_dict[pid]), ph_segment_ids)))#.astype(float)

    ## Read in the atl08 data
    if f08 is not None:
        ## classed flag (signal_photons)
        atl08_classed_pc_flag  = f08['/' + laser + '/signal_photons/classed_pc_flag'][...,]
        atl08_ph_segment_id = f08['/' + laser + '/signal_photons/ph_segment_id'][...,] # photon src 20 m seg id
        atl08_classed_pc_indx = f08['/' + laser + '/signal_photons/classed_pc_indx'][...,]

        ## set the classifications from atl08
        atl08_segment_id_msk = [True if x in segment_id else False for x in atl08_ph_segment_id]
        atl08_ph_segment_indx = np.array(list(map((lambda pid: segment_index_dict[pid]), atl08_ph_segment_id[atl08_segment_id_msk])))
        atl08_ph_index = np.array(atl08_ph_segment_indx + (atl08_classed_pc_indx[atl08_segment_id_msk] - 1), dtype=int)
        class_mask = atl08_ph_index < len(ph_segment_ids)
        ph_h_classed[atl08_ph_index[class_mask]] = atl08_classed_pc_flag[atl08_segment_id_msk][class_mask]

    
    ## Read in the atl24 data
    if f24 is not None:        
        atl24_classed_pc_flag  = f24['/' + laser + '/class_ph'][...,]
        atl24_classed_pc_indx = f24['/' + laser + '/index_ph'][...,]
        atl24_longitude = f24['/' + laser + '/lon_ph'][...,]
        atl24_latitude = f24['/' + laser + '/lat_ph'][...,]
        atl24_surface_h = f24['/' + laser + '/surface_h'][...,]
        atl24_ellipse_h = f24['/' + laser + '/ellipse_h'][...,]
        atl24_ortho_h = f24['/' + laser + '/ortho_h'][...,]

        class_40_mask = atl24_classed_pc_flag == 40
        ph_h_classed[atl24_classed_pc_indx] = atl24_classed_pc_flag

        # we also need to change the lon/lat/height values to the updated bathymetry values (we'll just do it to class 40)
        longitude[atl24_classed_pc_indx[class_40_mask]] = atl24_longitude[class_40_mask]
        latitude[atl24_classed_pc_indx[class_40_mask]] = atl24_latitude[class_40_mask]
        photon_h[atl24_classed_pc_indx[class_40_mask]] = atl24_ellipse_h[class_40_mask]
        photon_h_geoid[atl24_classed_pc_indx[class_40_mask]] = atl24_ortho_h[class_40_mask]
        
    
    ## set the photon height, either 'mean_tide' or 'geoid', else ellipsoid
    if water_surface == 'mean_tide':
        ph_height = photon_h_meantide
    elif water_surface == 'geoid':
        ph_height = photon_h_geoid
    else:
        ph_height = photon_h

    ## create the pandas dataframe            
    dataset = pd.DataFrame(
        {'latitude': latitude,
         'longitude': longitude,
         'photon_height': ph_height,
         'laser': laser_arr,
         'confidence': conf,
         'ref_elevation':ph_ref_elev,
         'ref_azimuth':ph_ref_azimuth,
         'ref_sat_alt':ph_altitude_sc,
         'ph_h_classed': ph_h_classed},
        columns=['latitude', 'longitude', 'photon_height', 'laser',
                 'confidence', 'ref_elevation', 'ref_azimuth', 'ref_sat_alt',
                 'ph_h_classed']
    )

    return(dataset)        

def fetch_atl08(atl03_fn):
    """fetch the associated atl08 file"""

    atl08_filter = utils.fn_basename2(atl03_fn).split('ATL03_')[1]
    this_atl08 = fetches.IceSat2(
        src_region=None, verbose=True, short_name='ATL08',
        filename_filter=atl08_filter
    )
    this_atl08.run()
    if len(this_atl08.results) == 0:
        utils.echo_warning_msg('could not locate associated atl08 file for {}'.format(atl08_filter))
        return(None)
    else:
        if this_atl08.fetch(this_atl08.results[0], check_size=True) == 0:
            return(os.path.join(this_atl08._outdir, this_atl08.results[0][1]))

def fetch_atl24(atl03_fn):
    """fetch the associated atl08 file"""

    atl24_filter = utils.fn_basename2(atl03_fn).split('ATL03_')[1]
    this_atl24 = fetches.IceSat2(
        src_region=None, verbose=True, short_name='ATL24',
        filename_filter=atl24_filter
    )
    this_atl24.run()
    if len(this_atl24.results) == 0:
        utils.echo_warning_msg('could not locate associated atl24 file for {}'.format(atl24_filter))
        return(None)
    else:
        if this_atl24.fetch(this_atl24.results[0], check_size=True) == 0:
            return(os.path.join(this_atl24._outdir, this_atl24.results[0][1]))        

# atl03
atl03_fn = sys.argv[1]
atl03 = h5.File(atl03_fn)

# atl08
atl08_fn = fetch_atl08(atl03_fn)
if os.path.exists(atl08_fn):
    utils.echo_msg('using ATL08: {}'.format(atl08_fn))
    atl08 = h5.File(atl08_fn)
else:
    atl08 = None

# atl24
#atl24_fn = sys.argv[2]
#atl08_fn = fetch_atl08(atl03_fn) # not available on earthdata yet, using test data
atl24_fn = 'ATL24_' + utils.fn_basename2(atl03_fn).split('ATL03_')[1].split('.')[0] + '_001_01.h5'
if os.path.exists(atl24_fn):
    utils.echo_msg('using ATL24: {}'.format(atl24_fn))
    atl24 = h5.File(atl24_fn)
else:
    utils.echo_msg('could not locate ATL24: {}'.format(atl24_fn))
    atl24 = None

# we'll return ground (1) and bathymetry (40)
classes = [1, 40]

for laser_num in range(1, 3):
    for orientation_num in range(0, 1):
        dataset = read_atl_data(atl03, str(laser_num), orientation_num, f08=atl08, f24=atl24)

        ## return only classes mentioned in `classes`
        dataset = dataset[(np.isin(dataset['ph_h_classed'], classes))]

        ## dump the data to xyz
        dataset = np.vstack((dataset['longitude'], dataset['latitude'], dataset['photon_height'])).transpose()
        for point in dataset:
            this_xyz = xyzfun.XYZPoint(x=point[0], y=point[1], z=point[2])
            this_xyz.dump(dst_port=sys.stdout, encode=False, include_w=False, include_u=False)

# close h5 objs
atl03.close()

if atl08 is not None:
    atl08.close()
    
if atl24 is not None:
    atl24.close()

### End
