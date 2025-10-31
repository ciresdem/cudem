### cdse.py
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
##
## fetches.py is part of CUDEM
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
##
### Code:
#import os
#import json
import requests
#import boto3
from tqdm import tqdm
#import time
import netrc
import datetime
from cudem import utils
from cudem import regions
from cudem.fetches import fetches
import certifi
try:
    from urllib.parse import urlparse
    from urllib.request import urlopen, Request, build_opener, HTTPCookieProcessor
    from urllib.error import HTTPError, URLError
except ImportError:
    from urlparse import urlparse
    from urllib2 import urlopen, Request, HTTPError, URLError, build_opener, HTTPCookieProcessor

import xml.etree.ElementTree as ET

## CDSE
class CDSE(fetches.FetchModule):
    """
    """
    
    def __init__(self, collection_name='SENTINEL-2', product_type='S2MSI1C',
                 max_cloud_cover=1, time_start='', time_end='', **kwargs):
        super().__init__(name='cdse', **kwargs)

        # base URL of the product catalogue
        self.catalogue_odata_url = "https://catalogue.dataspace.copernicus.eu/odata/v1"
        self._auth_openeo_url = 'https://identity.dataspace.copernicus.eu/auth/realms/CDSE/protocol/openid-connect/token'
        
        self.collection_name = collection_name
        self.product_type = product_type
        self.max_cloud_cover = utils.int_or(max_cloud_cover, 1)
        self.aoi = self.region.export_as_wkt().replace('POLYGON ', 'POLYGON')
        self.time_start = time_start
        self.time_end = time_end

        #print(self.time_start, self.time_end)
        if self.time_start != '' or self.time_end != '':
            self.time_start = datetime.datetime.fromisoformat(self.time_start).isoformat(timepsec='milliseconds') + 'Z' if self.time_start != '' else ''
            self.time_end = datetime.datetime.fromisoformat(self.time_end).isoformat(timepsec='milliseconds') + 'Z' if self.time_end != '' else ''
            
        self.access_token = self.get_auth_token()
        self.headers = {
            'Authorization': f'Bearer {self.access_token}',
        }


    def get_userpass(self):
        try:
            info = netrc.netrc()
            username, account, password \
                = info.authenticators(urlparse(self._auth_openeo_url).hostname)
            errprefix = 'netrc error: '
        except Exception as e:
            if (not ('No such file' in str(e))):
                print('netrc error: {0}'.format(str(e)))
            username = None
            password = None

        return(username, password)

    
    def get_auth_token(self):
        try:
            info = netrc.netrc()
            username, account, password \
                = info.authenticators(urlparse(self._auth_openeo_url).hostname)
            errprefix = 'netrc error: '
        except Exception as e:
            if (not ('No such file' in str(e))):
                print('netrc error: {0}'.format(str(e)))
            username = None
            password = None

        
        data = {
            "client_id": "cdse-public",
            "grant_type": "password",
            "username": username,
            "password": password,
        }

        try:
            # Make the POST request to the authentication server.
            response = requests.post(self._auth_openeo_url, data=data)
            response.raise_for_status() # Raise an exception for bad status codes (4xx or 5xx)

            # Extract the access token from the response.
            access_token = response.json()["access_token"]
            print("Successfully retrieved access token.")

        except requests.exceptions.RequestException as e:
            print(f"An error occurred: {e}")
            access_token = None

        return(access_token)

        
    def run(self):
        search_query = (f"{self.catalogue_odata_url}/Products?$filter=Collection/Name eq '{self.collection_name}'"
                        f" and Attributes/OData.CSC.StringAttribute/any(att:att/Name eq 'productType'"
                        f" and att/OData.CSC.StringAttribute/Value eq '{self.product_type}')"
                        f" and OData.CSC.Intersects(area=geography'SRID=4326;{self.aoi}')")
        
        if self.time_start != '':
            search_query += f' and ContentDate/Start gt {self.time_start}'

        if self.time_end != '':
            search_query += f' and ContentDate/Start lt {self.time_end}'
            
        response = requests.get(search_query).json()
        results = response['value']
        # Select identifier of the first product
        for result in results:
            product_identifier = result['Id']
            product_name = result['Name']

            url = f"{self.catalogue_odata_url}/Products({product_identifier})/Nodes({product_name})/Nodes(MTD_MSIL1C.xml)/$value"
            _req = fetches.Fetch(url, headers=self.headers, allow_redirects=False).fetch_req()
            while _req.status_code in (301, 302, 303, 307):
                url = _req.headers["Location"]
                _req = fetches.Fetch(url, headers=self.headers, allow_redirects=False).fetch_req()

            _req = fetches.Fetch(url, headers=self.headers, verify=True, allow_redirects=True).fetch_req()

            root = ET.fromstring(_req.text.encode('utf-8'))

            # Get the location of individual bands in Sentinel-2 granule
            band_location = []
            band_location.append(f"{product_name}/{root[0][0][12][0][0][1].text}.jp2".split("/"))
            band_location.append(f"{product_name}/{root[0][0][12][0][0][2].text}.jp2".split("/"))
            band_location.append(f"{product_name}/{root[0][0][12][0][0][3].text}.jp2".split("/"))

            bands = []
            for band_file in band_location:
                url = f"{self.catalogue_odata_url}/Products({product_identifier})/Nodes({product_name})/Nodes({band_file[1]})/Nodes({band_file[2]})/Nodes({band_file[3]})/Nodes({band_file[4]})/$value"
                #print(url)
                self.add_entry_to_results(
                    url,
                    band_file[4],
                    'SENTINEL 2LA'
                )

                
class Sentinel2(CDSE):
    def __init__(self, **kwargs):
        super().__init__(name='sentinel2', collection_name='SENTINEL-2', product_type='S2MSI1C', **kwargs)
        
            
# class CDSE_s3(fetches.FetchModule):
#     """
#     """
    
#     def __init__(self, eo_product_name=None, **kwargs):
#         super().__init__(name='cdse', **kwargs)
#         self._openeo_url = 'https://openeo.dataspace.copernicus.eu'
#         self._auth_openeo_url = 'https://identity.dataspace.copernicus.eu/auth/realms/CDSE/protocol/openid-connect/token'
#         self.eo_product_name = eo_product_name
#         #self._auth_openeo_url = 'https://identity.dataspace.copernicus.eu'
#         ## Set up the earthdata credentials, and add it to our headers
#         #credentials = fetches.get_credentials(None, authenticator_url=self._auth_openeo_url)
#         #utils.echo_msg(credentials)
#         # self.headers = {
#         #     'Authorization': 'Basic {0}'.format(credentials),
#         #     'User-Agent': ('Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:89.0) '
#         #                    'Gecko/20100101 Firefox/89.0')
#         # }

#         access_token = self.get_auth_token()
#         self.headers = {
#             'Authorization': f'Bearer {self.access_token}',
#         }
        
#         self.config = {
#             "auth_server_url": "https://identity.dataspace.copernicus.eu/auth/realms/CDSE/protocol/openid-connect/token",
#             "odata_base_url": "https://catalogue.dataspace.copernicus.eu/odata/v1/Products",
#             "s3_endpoint_url": "https://eodata.dataspace.copernicus.eu",
#         }

        
#     def get_userpass(self):
#         try:
#             info = netrc.netrc()
#             username, account, password \
#                 = info.authenticators(urlparse(self._auth_openeo_url).hostname)
#             errprefix = 'netrc error: '
#         except Exception as e:
#             if (not ('No such file' in str(e))):
#                 print('netrc error: {0}'.format(str(e)))
#             username = None
#             password = None

#         return(username, password)
        
#     def get_auth_token(self):
#         try:
#             info = netrc.netrc()
#             username, account, password \
#                 = info.authenticators(urlparse(self._auth_openeo_url).hostname)
#             errprefix = 'netrc error: '
#         except Exception as e:
#             if (not ('No such file' in str(e))):
#                 print('netrc error: {0}'.format(str(e)))
#             username = None
#             password = None

        
#         data = {
#             "client_id": "cdse-public",
#             "grant_type": "password",
#             "username": username,
#             "password": password,
#         }

#         try:
#             # Make the POST request to the authentication server.
#             response = requests.post(self._auth_openeo_url, data=data)
#             response.raise_for_status() # Raise an exception for bad status codes (4xx or 5xx)

#             # Extract the access token from the response.
#             access_token = response.json()["access_token"]
#             print("Successfully retrieved access token.")

#         except requests.exceptions.RequestException as e:
#             print(f"An error occurred: {e}")
#             access_token = None

#         return(access_token)


#     def get_eo_product_details(self, eo_product_name):
#         """
#         Retrieve EO product details using the OData API to determine the S3 path.
#         """
        
#         odata_url = f"{self.config['odata_base_url']}?$filter=Name eq '{eo_product_name}'"
#         response = requests.get(odata_url, headers=self.headers)
#         if response.status_code == 200:
#             eo_product_data = response.json()["value"][0]
#             return(eo_product_data["Id"], eo_product_data["S3Path"])
#         else:
#             print(f"Failed to retrieve EO product details. Status code: {response.status_code}")
#             exit(1)

            
#     def get_temporary_s3_credentials(self):
#         """
#         Create temporary S3 credentials by calling the S3 keys manager API.
#         """
        
#         credentials_response = requests.post("https://s3-keys-manager.cloudferro.com/api/user/credentials", headers=self.headers)
#         if credentials_response.status_code == 200:
#             s3_credentials = credentials_response.json()
#             print("Temporary S3 credentials created successfully.")
#             print(f"access: {s3_credentials['access_id']}")
#             print(f"secret: {s3_credentials['secret']}")
#             return(s3_credentials)        
#         else:
#             print(f"Failed to create temporary S3 credentials. Status code: {credentials_response.status_code}")
#             print("Product download aborted.")
#             exit(1)

#     def format_filename(self, filename, length=40):
#         """
#         Format a filename to a fixed length, truncating if necessary.
#         """
        
#         if len(filename) > length:
#             return(filename[:length - 3] + '...')
#         else:
#             return(filename.ljust(length))

        
#     def download_file_s3(self, s3, bucket_name, s3_key, local_path, failed_downloads):
#         """
#         Download a file from S3 with a progress bar.
#         Track failed downloads in a list.
#         """
        
#         try:
#             file_size = s3.head_object(Bucket=bucket_name, Key=s3_key)['ContentLength']
#             formatted_filename = self.format_filename(os.path.basename(local_path))
#             with tqdm(total=file_size, unit='B', unit_scale=True, desc=formatted_filename, ncols=80, bar_format='{desc:.40}|{bar:20}| {percentage:3.0f}% {n_fmt}/{total_fmt}B') as pbar:
#                 def progress_callback(bytes_transferred):
#                     pbar.update(bytes_transferred)

#                 s3.download_file(bucket_name, s3_key, local_path, Callback=progress_callback)
#         except Exception as e:
#             print(f"Failed to download {s3_key}. Error: {e}")
#             failed_downloads.append(s3_key)

            
#     def traverse_and_download_s3(self, s3_resource, bucket_name, base_s3_path, local_path, failed_downloads):
#         """
#         Traverse the S3 bucket and download all files under the specified prefix.
#         """
#         bucket = s3_resource.Bucket(bucket_name)
#         files = bucket.objects.filter(Prefix=base_s3_path)
        
#         for obj in files:
#             s3_key = obj.key
#             relative_path = os.path.relpath(s3_key, base_s3_path)
#             local_path_file = os.path.join(local_path, relative_path)
#             local_dir = os.path.dirname(local_path_file)
#             os.makedirs(local_dir, exist_ok=True)
#             self.download_file_s3(s3_resource.meta.client, bucket_name, s3_key, local_path_file, failed_downloads)

#     def run(self):
#         # # Step 1: Retrieve the access token
#         # access_token = get_access_token(config, args.username, args.password)

#         # # Step 2: Set up headers for API calls
#         # headers = {
#         #     "Authorization": f"Bearer {access_token}",
#         #     "Accept": "application/json"
#         # }

#         # Step 3: Get EO product details (including S3 path)
#         eo_product_id, s3_path = self.get_eo_product_details(self.eo_product_name)
#         bucket_name, base_s3_path = s3_path.lstrip('/').split('/', 1)

#         # Step 4: Get temporary S3 credentials
#         s3_credentials = self.get_temporary_s3_credentials()

#         # Step 5: Set up S3 client and resource with temporary credentials
#         time.sleep(5)  # Ensure the key pair is installed
#         s3_resource = boto3.resource(
#             's3',
#             endpoint_url=self.config["s3_endpoint_url"],
#             aws_access_key_id=s3_credentials["access_id"],
#             aws_secret_access_key=s3_credentials["secret"]
#         )

#         # Step 6: Create the top-level folder and start download
#         top_level_folder = self.eo_product_name
#         os.makedirs(top_level_folder, exist_ok=True)
#         failed_downloads = []
#         self.traverse_and_download_s3(s3_resource, bucket_name, base_s3_path, top_level_folder, failed_downloads)

#         # bucket = s3_resource.Bucket(bucket_name)
#         # files = bucket.objects.filter(Prefix=base_s3_path)
#         # #print(files)

#         # for obj in files:
#         #     s3_key = obj.key
#         #     relative_path = os.path.relpath(s3_key, base_s3_path)
#         #     #print(bucket_name, s3_key)
#         #     #local_path_file = os.path.join(local_path, relative_path)
#         #     #local_dir = os.path.dirname(local_path_file)
#         #     #os.makedirs(local_dir, exist_ok=True)
#         #     #self.download_file_s3(s3_resource.meta.client, bucket_name, s3_key, local_path_file, failed_downloads)

#         #     self.add_entry_to_results(
#         #         f'{self.config["s3_endpoint_url"]}/{s3_key}',
#         #         os.path.basename(s3_key),
#         #         'SENTINEL 2LA'
#         #     )

        
#         # # Step 7: Print final status
#         # if not failed_downloads:
#         #     print("Product download complete.")
#         # else:
#         #     print("Product download incomplete:")
#         #     for failed_file in failed_downloads:
#         #         print(f"- {failed_file}")

#         # # Step 7: Delete the temporary S3 credentials
#         # delete_response = requests.delete(f"https://s3-keys-manager.cloudferro.com/api/user/credentials/access_id/{s3_credentials['access_id']}", headers=self.headers)
#         # if delete_response.status_code == 204:
#         #     print("Temporary S3 credentials deleted successfully.")
#         # else:
#         #     print(f"Failed to delete temporary S3 credentials. Status code: {delete_response.status_code}")

        
#         return(self)

### End
