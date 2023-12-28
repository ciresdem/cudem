#!/usr/bin/env python
# coding: utf-8

'''
THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXpRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.'''

### cshelph.py
### Commentary:
## https://github.com/nmt28/C-SHELPh.git
##
## Changes were applied to make compatible with cudem
##
### Code:

import numpy as np
from pyproj import Proj
from pyproj import Transformer
import pandas as pd
import utm
# import xarray as xr
# import earthaccess # earthaccess removed to use fetches instead
    
def convert_wgs_to_utm(lat, lon):
	easting, northing, num, letter = utm.from_latlon(lat, lon)
	if letter >= 'N':
		epsg = 'epsg:326' + str(num)
	elif letter < 'N':
		epsg = 'epsg:327' + str(num)
	else:
		print('Error Finding UTM')
		
	return epsg

def orthometric_correction(lat, lon, Z, epsg):
    # transform ellipsod (WGS84) height to orthometric height
    transformerh = Transformer.from_crs("epsg:4326", "epsg:3855", always_xy=True)
    X_egm08, Y_egm08, Z_egm08 = transformerh.transform(lon, lat, Z)
    
    # transform WGS84 proj to local UTM
    myproj = Proj(epsg)
    X_utm, Y_utm = myproj(lon, lat)
    
    return Y_utm, X_utm, Z_egm08

    
def count_ph_per_seg(ph_index_beg, photon_h): # DEpRECATED
    
    ph_index_beg = ph_index_beg[ph_index_beg!=0]
    
    # add an extra val at the end of array for upper bounds
    ph_index_beg = np.hstack((ph_index_beg, len(photon_h)+1))-1
    
    photon_id = []
    #iterate over the photon indexes (ph_index_beg)
    for i, num in enumerate(np.arange(0, len(ph_index_beg)-1, 1)):
        photon_id.append(len(photon_h[ph_index_beg[i]:ph_index_beg[i+1]]))
    photon_id = np.array(photon_id)
    
    return photon_id

    
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
    
     # Duplicate dataframe
    dataset1 = dataset
    
    # Cut lat bins
    lat_bins = pd.cut(dataset['latitude'], lat_bin_number, labels = np.array(range(lat_bin_number)))
    
    # Add bins to dataframe
    dataset1['lat_bins'] = lat_bins
    
    # Cut height bins
    height_bins = pd.cut(dataset['photon_height'], height_bin_number, labels = np.round(np.linspace(dataset['photon_height'].min(), dataset['photon_height'].max(), num=height_bin_number), decimals = 1))
    
    # Add height bins to dataframe
    dataset1['height_bins'] = height_bins
    dataset1 = dataset1.reset_index(drop=True)

    return dataset1

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
    mean = np.nanmean(sea_height, axis=0)
    sd = np.nanstd(sea_height, axis=0)
    sea_height_1 = np.where((sea_height > (mean + 2*sd)) | (sea_height < (mean - 2*sd)), np.nan, sea_height).tolist()
    
    return sea_height_1

def get_water_temp(data_path, latitude, longitude):
    return(20)

# def get_water_temp(data_path, latitude, longitude):

#     try:
#         auth = earthaccess.login(strategy="netrc")
#     except:
#         print("Login credentials not found. Sign in manually, your netrc file will be created for next time")
#         auth = earthaccess.login(strategy="interactive", persist=True)
    
#     # Get date from data filename
#     file_date = data_path[-33:-25]
    
#     date_range = file_date[0:4] + '-' + file_date[4:6] + '-' + file_date[6:]
    
#     location_df = pd.DataFrame({'longitude':longitude,'latitude':latitude})
    
#     location_df = location_df.dropna(axis=0)
    
#     minx, miny, maxx, maxy = np.min(location_df['longitude']), np.min(location_df['latitude']), np.max(location_df['longitude']), np.max(location_df['latitude'])
    
#     lat_med = np.median(location_df['latitude'])
#     lon_med = np.median(location_df['longitude'])

#     Query = earthaccess.collection_query()

#     # Use chain methods to customize our query
#     collections = Query.keyword('GHRSST Level 4 CMC0.1deg Global Foundation Sea Surface Temperature Analysis').bounding_box(minx,miny,maxx,maxy).temporal(date_range,date_range).get(10)
    
#     short_name = collections[0]["umm"]["ShortName"]

#     Query = earthaccess.granule_query().short_name(short_name).version("3.0").bounding_box(minx,miny,maxx,maxy).temporal("2020-01-01","2020-01-01")
    
#     granules = Query.get(10)
    
#     ds_L3 = xr.open_mfdataset(earthaccess.open(granules),combine='nested',concat_dim='time',coords='minimal')
    
#     sea_temp = ds_L3['analysed_sst'].sel(lat=lat_med,lon=lon_med,method='nearest').load()
    
#     sst = round(np.nanmedian(sea_temp.values)-273,2)
    
#     return sst

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
    ph_ref_azimuth = ph_ref_azimuth[photon_z<=water_surface]
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
