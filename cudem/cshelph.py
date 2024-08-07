#!/usr/bin/env python
# coding: utf-8

'''
THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXpRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.'''

import numpy as np
import h5py as h5
from pyproj import Transformer
import pandas as pd
import utm

def read_atl03(h5_file, laser_num):
    # Read File
    f = h5.File(h5_file,'r')
    
    # Select a laser
    orientation = f['/orbit_info/sc_orient'][0]
    
    # selects the strong beams only [we can include weak beams later on]
    orientDict = {0:'l', 1:'r', 21:'error'}
    laser = 'gt' + laser_num + orientDict[orientation]
    
    # Read in the required photon level data
    photon_h = f['/' + laser + '/heights/h_ph'][...,]
    latitude = f['/' + laser + '/heights/lat_ph'][...,]
    longitude = f['/' + laser + '/heights/lon_ph'][...,]
    conf = f['/' + laser + '/heights/signal_conf_ph/'][...,0]

    # params needed for refraction correction
    
    ref_elev = f['/' + laser + '/geolocation/ref_elev'][...,]
    ref_azimuth = f['/' + laser + '/geolocation/ref_azimuth'][...,]
    ph_index_beg = f['/' + laser + '/geolocation/ph_index_beg'][...,]
    segment_id = f['/' + laser + '/geolocation/segment_id'][...,]
    altitude_sc = f['/' + laser + '/geolocation/altitude_sc'][...,]
    seg_ph_count = f['/' + laser + '/geolocation/segment_ph_cnt'][...,]
    
    return latitude, longitude, photon_h, conf, ref_elev, ref_azimuth, ph_index_beg, segment_id, altitude_sc, seg_ph_count
    
def convert_wgs_to_utm(lat, lon):
	easting, northing, num, letter = utm.from_latlon(lat, lon)
	if letter >= 'N':
		epsg = 'epsg:326' + str(num)
	elif letter < 'N':
		epsg = 'epsg:327' + str(num)
	else:
		print('Error Finding UTM')
		
	return epsg
    
def ref_linear_interp(photon_count, ref_elev):

    arr = []
    for i in range(len(ref_elev)):
        try:
            min = ref_elev[i-1]
            max = ref_elev[i]
        except:
            min = ref_elev[i]
            max = ref_elev[i]
            
        try:
            min = ref_elev[i]
            max = ref_elev[i+1]
        except:
            min = ref_elev[i]
            max = ref_elev[i]
        
        if min==max:
            sub = np.full((photon_count[i]), min)
            arr.append(sub)
        else:
            sub_tmp = np.linspace(min, max, photon_count[i]+1)
            if len(sub_tmp)>1:
                sub = np.linspace(sub_tmp[1], sub_tmp[-1], photon_count[i])
                arr.append(sub)
            else:
                arr.append(sub_tmp)

    return np.concatenate(arr, axis=None).ravel()

def bin_data(dataset, lat_res, height_res):
    '''Bin data along vertical and horizontal scales for later segmentation'''
    
    # Calculate number of bins required both vertically and horizontally with resolution size
    lat_bin_number = round(abs(dataset['latitude'].min() - dataset['latitude'].max())/lat_res)
    height_bin_number = round(abs(dataset['photon_height'].min() - dataset['photon_height'].max())/height_res)

    if (lat_bin_number > 0 and height_bin_number > 0):    
         # Duplicate dataframe
        dataset1 = dataset

        # Cut lat bins
        lat_bins = pd.cut(dataset['latitude'], lat_bin_number, labels = np.array(range(lat_bin_number)))

        # Add bins to dataframe
        dataset1['lat_bins'] = lat_bins

        # Cut height bins
        height_bins = pd.cut(
            dataset['photon_height'], height_bin_number, labels = np.round(
                np.linspace(dataset['photon_height'].min(), dataset['photon_height'].max(), num=height_bin_number),
                decimals = 1
            )
        )

        # Add height bins to dataframe
        dataset1['height_bins'] = height_bins
        dataset1 = dataset1.reset_index(drop=True)

        return dataset1

    return None

def get_sea_height(binned_data, surface_buffer=-0.5):
    '''Calculate mean sea height for easier calculation of depth and cleaner figures'''
    
    # Create sea height list
    sea_height = []
    
    # Group data by latitude
    binned_data_sea = binned_data[(binned_data['photon_height'] > surface_buffer)] # Filter out subsurface data
    grouped_data = binned_data_sea.groupby(['lat_bins'], group_keys=True)
    data_groups = dict(list(grouped_data))
    
    # Loop through groups and return average sea height
    for k,v in data_groups.items():
        # Create new dataframe based on occurance of photons per height bin
        new_df = pd.DataFrame(v.groupby('height_bins').count())
        
        # Return the bin with the highest count
        largest_h_bin = new_df['latitude'].argmax()
        
        # Select the index of the bin with the highest count
        largest_h = new_df.index[largest_h_bin]
        
        # Calculate the median value of all values within this bin
        lat_bin_sea_median = v.loc[v['height_bins']==largest_h, 'photon_height'].median()
        
        # Append to sea height list
        sea_height.append(lat_bin_sea_median)
        del new_df
        
    # Filter out sea height bin values outside 2 SD of mean.
    if np.all(np.isnan(sea_height)):
        return(None)
    
    mean = np.nanmean(sea_height, axis=0)
    sd = np.nanstd(sea_height, axis=0)
    sea_height_1 = np.where((sea_height > (mean + 2*sd)) | (sea_height < (mean - 2*sd)), np.nan, sea_height).tolist()
    
    return sea_height_1

def refraction_correction(water_temp, water_surface, wavelength, photon_ref_elev, ph_ref_azimuth, photon_z, photon_x, photon_y, ph_conf, satellite_altitude):
    
    '''
    WTemp; there is python library that pulls water temp data
    water_surface is the value surface height
    Wavelength is fixed
    '''
    
    # Only process photons below water surface model
    photon_x = photon_x[photon_z<=water_surface]
    photon_y = photon_y[photon_z<=water_surface]
    photon_ref_elev = photon_ref_elev[photon_z<=water_surface]
    satellite_altitude = satellite_altitude[photon_z<=water_surface]
    ph_ref_azimuth = ph_ref_azimuth[photon_z<=water_surface]#.astype(float)
    ph_conf = ph_conf[photon_z<=water_surface]
    photon_z = photon_z[photon_z<=water_surface]
    
    # Refraction coefficient #
    a = -0.000001501562500
    b = 0.000000107084865
    c = -0.000042759374989
    d = -0.000160475520686
    e = 1.398067112092424
    wl = wavelength
    
    # refractive index of air
    n1 = 1.00029
    
    # refractive index of water
    n2 = (a*water_temp**2) + (b*wl**2) + (c*water_temp) + (d*wl) + e
    
    # assumption is 0.25416
    # This example is refractionCoef = 0.25449
    # 1.00029 is refraction of air constant
    #correction_coef = (1-(n1/n2))
    #########################
    
    # read photon ref_elev to get theta1
    # Does not account for curvature of Earth
    theta1 = np.pi/2 - photon_ref_elev
    
    # H = orbital altitude of IS2 (496km as mean)
    # H = 496km. we pass in the mean of the orbit from /geolocation/altitude_sc/
    # Diff from min to max of 100m over an orbit is 0.02% at 496km
    # More error probably introduced from Re (mean Earth radius) than intra-orbit changes in altitude
    # H = 496
    H = satellite_altitude/1000
    # Re = Radius of Earth (6371km mean)
    Re = 6371
    
    # remove as potentially more inaccurate of a correcttion
    #theta1 = np.arctan((H*np.tan(theta_1))/Re)
    
    # eq 1. Theta2
    theta2 = np.arcsin(((n1*np.sin(theta1))/n2))
    
    # eq 3. S
    # Approximate water Surface = 1.5
    # D  = raw uncorrected depth
    D = water_surface - photon_z
    
    # For Triangle DTS
    S = D/np.cos(theta1)
    
    # eq 2. R
    R = (S*n1)/n2
    Gamma = (np.pi/2)-theta1
    
    # For triangle RpS
    # phi is an angle needed
    phi = theta1-theta2
    
    # p is the difference between raw and corrected YZ location
    p = np.sqrt(R**2 + S**2 - 2*R*S*np.cos(phi))
    
    # alpha is an angle needed
    alpha = np.arcsin((R*np.sin(phi))/p)
    
    # Beta angle needed for Delta Y an d Delta Z
    Beta = Gamma - alpha
    
    # Delta Y
    DY = p*np.cos(Beta)
    
    # Delta Z
    DZ = p*np.sin(Beta)
    
    # Delta Easting
    DE = DY*np.sin(ph_ref_azimuth)
    
    # Delta Northing
    DN = DY*np.cos(ph_ref_azimuth)
    
    out_x = photon_x + DE
    out_y = photon_y + DN
    out_z = photon_z + DZ

    return(out_x, out_y, out_z, ph_conf, photon_x, photon_y, photon_z, ph_ref_azimuth, photon_ref_elev) # We are most interested in out_x, out_y, out_z

def get_bath_height(binned_data, percentile, WSHeight, height_resolution):
    '''Calculate bathymetry level per bin based on horizontal resolution'''
    # Create sea height list
    bath_height = []
    
    geo_photon_height = []
    geo_longitude = []
    geo_latitude = []
    ids = []

    # Group data by latitude
    # Filter out surface data that are two bins below median surface value calculated above
    binned_data_bath = binned_data[(binned_data['photon_height'] < WSHeight - (height_resolution * 2))]
        
    grouped_data = binned_data_bath.groupby(['lat_bins'], group_keys=True)
    data_groups = dict(list(grouped_data))
    
    # Create a percentile threshold of photon counts in each grid, grouped by both x and y axes.
    count_threshold = np.percentile(binned_data.groupby(['lat_bins', 'height_bins']).size().reset_index().groupby('lat_bins')[[0]].max(), percentile)
    #if count threshold is 0 or 1, set it to 2
    count_threshold = 2 if count_threshold <= 1 else count_threshold
    # Loop through groups and return average bathy height
    for k,v in data_groups.items():
        new_df = pd.DataFrame(v.groupby('height_bins').count())
        bath_bin = new_df['cor_latitude'].argmax()
        bath_bin_h = new_df.index[bath_bin]
        
        # Set threshold of photon counts per bin
        if new_df.iloc[bath_bin]['latitude'] >= count_threshold:
            ids.append(v.loc[v['height_bins']==bath_bin_h].index.values.astype(int))            
            geo_photon_height.append(v.loc[v['height_bins']==bath_bin_h, 'cor_photon_height'].values)
            geo_longitude.append(v.loc[v['height_bins']==bath_bin_h, 'cor_longitude'].values)
            geo_latitude.append(v.loc[v['height_bins']==bath_bin_h, 'cor_latitude'].values)
            
            bath_bin_median = v.loc[v['height_bins']==bath_bin_h, 'cor_photon_height'].median()
            bath_height.append(bath_bin_median)
            del new_df
            
        else:
            bath_height.append(np.nan)

            del new_df

    if len(geo_longitude) == 0 or len(geo_latitude) == 0:
        return(None, None)
    
    geo_longitude_list = np.concatenate(geo_longitude).ravel().tolist()
    geo_latitude_list = np.concatenate(geo_latitude).ravel().tolist()
    geo_photon_list = np.concatenate(geo_photon_height).ravel().tolist()
    ids_list = np.concatenate(ids).ravel().tolist()
    geo_depth = WSHeight - geo_photon_list
        
    geo_df = pd.DataFrame({'longitude': geo_longitude_list, 'latitude':geo_latitude_list, 'photon_height': geo_photon_list, 'depth':geo_depth, 'ids':ids_list})
    
    del geo_longitude_list, geo_latitude_list, geo_photon_list, ids_list
    
    return bath_height, geo_df

def get_bin_height(binned_data, percentile, surface_buffer=-0.5):
    '''Calculate mean sea height for easier calculation of depth and cleaner figures'''
    
    # Create sea height list
    sea_height = []
    bin_lat = []
    bin_lon = []
    
    # Group data by latitude
    binned_data_sea = binned_data[(binned_data['photon_height'] > surface_buffer)] # Filter out subsurface data
    grouped_data = binned_data_sea.groupby(['lat_bins'], group_keys=True)
    data_groups = dict(list(grouped_data))

    # Create a percentile threshold of photon counts in each grid, grouped by both x and y axes.
    count_threshold = np.percentile(binned_data.groupby(['lat_bins', 'height_bins']).size().reset_index().groupby('lat_bins')[[0]].max(), percentile)

    # Loop through groups and return average sea height
    for k,v in data_groups.items():
        # Create new dataframe based on occurance of photons per height bin
        new_df = pd.DataFrame(v.groupby('height_bins').count())

        # Return the bin with the highest count
        largest_h_bin = new_df['latitude'].argmax()

        # Select the index of the bin with the highest count
        largest_h = new_df.index[largest_h_bin]
        
        # Set threshold of photon counts per bin
        if new_df.iloc[largest_h_bin]['latitude'] >= count_threshold:        
            
            # Calculate the median value of all values within this bin
            lat_bin_sea_median = v.loc[v['height_bins']==largest_h, 'photon_height'].max()
            lat_bin_median = v.loc[v['height_bins']==largest_h, 'latitude'].median()
            lon_bin_median = v.loc[v['height_bins']==largest_h, 'longitude'].median()

            # Append to sea height list
            sea_height.append(lat_bin_sea_median)
            bin_lat.append(lat_bin_median)
            bin_lon.append(lon_bin_median)
            del new_df
        else:
            del new_df
        
    # Filter out sea height bin values outside 2 SD of mean.
    if np.all(np.isnan(sea_height)):
        return(None)
    
    mean = np.nanmean(sea_height, axis=0)
    sd = np.nanstd(sea_height, axis=0)
    sea_height_1 = np.where((sea_height > (mean + 2*sd)) | (sea_height < (mean - 2*sd)), np.nan, sea_height).tolist()
    
    return bin_lat, bin_lon, sea_height_1
