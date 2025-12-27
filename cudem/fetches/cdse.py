### cdse.py
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
##
## cdse.py is part of CUDEM
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
## Fetch data from the Copernicus Data Space Ecosystem (CDSE).
##
### Code:

import requests
import netrc
import datetime
import xml.etree.ElementTree as ET
from urllib.parse import urlparse
from typing import Optional, List, Dict, Tuple

from cudem import utils
from cudem.fetches import fetches

## ==============================================
## Constants
## ==============================================
CDSE_CATALOGUE_URL = "https://catalogue.dataspace.copernicus.eu/odata/v1"
CDSE_AUTH_URL = "https://identity.dataspace.copernicus.eu/auth/realms/CDSE/protocol/openid-connect/token"

## ==============================================
## CDSE Module
## ==============================================
class CDSE(fetches.FetchModule):
    """Copernicus Data Space Ecosystem Fetch Module

    Fetches satellite imagery (e.g., Sentinel-2) from the Copernicus Data Space Ecosystem.
    Requires a valid account in your .netrc file.


    Configuration Example:
    < cdse:collection_name=SENTINEL-2:product_type=S2MSI1C:time_start='':time_end='':max_cloud_cover=1 >
    """
    
    def __init__(self, collection_name: str = 'SENTINEL-2', product_type: str = 'S2MSI1C',
                 max_cloud_cover: Optional[float] = None, time_start: str = '', time_end: str = '', **kwargs):
        super().__init__(name='cdse', **kwargs)

        self.collection_name = collection_name
        self.product_type = product_type
        self.max_cloud_cover = utils.float_or(max_cloud_cover)
        
        ## Geometry for OData query
        if self.region:
            self.aoi = self.region.export_as_wkt().replace('POLYGON ', 'POLYGON')
        else:
            self.aoi = None

        ## Format Timestamps
        self.time_start = self._format_date(time_start)
        self.time_end = self._format_date(time_end)
            
        ## Authentication
        self.access_token = self.get_auth_token()
        if self.access_token:
            self.headers = {'Authorization': f'Bearer {self.access_token}'}
        else:
            self.headers = {}
            utils.echo_warning_msg("Could not acquire CDSE Access Token. Requests may fail.")

            
    def _format_date(self, date_str: str) -> str:
        """Formats an ISO date string for the OData filter."""
        
        if not date_str:
            return ''
        try:
            dt = datetime.datetime.fromisoformat(date_str)
            return dt.isoformat(timespec='milliseconds') + 'Z'
        except ValueError:
            return ''

        
    def get_credentials_from_netrc(self) -> Tuple[Optional[str], Optional[str]]:
        """Retrieve username and password from .netrc."""
        
        try:
            info = netrc.netrc()
            auth = info.authenticators(urlparse(CDSE_AUTH_URL).hostname)
            if auth:
                return auth[0], auth[2]
        except Exception as e:
            if 'No such file' not in str(e):
                utils.echo_error_msg(f"netrc error: {e}")
        return None, None

    
    def get_auth_token(self) -> Optional[str]:
        """Authenticate with CDSE and retrieve an access token."""
        
        username, password = self.get_credentials_from_netrc()
        
        if not username or not password:
            utils.echo_warning_msg("No credentials found in .netrc for CDSE.")
            return None

        data = {
            "client_id": "cdse-public",
            "grant_type": "password",
            "username": username,
            "password": password,
        }

        try:
            response = requests.post(CDSE_AUTH_URL, data=data)
            response.raise_for_status() # Raise an exception for bad status codes (4xx or 5xx)
            token = response.json().get("access_token")
            if self.verbose:
                utils.echo_msg("Successfully retrieved CDSE access token.")
            return token
        except requests.exceptions.RequestException as e:
            utils.echo_error_msg(f"CDSE Authentication failed: {e}")
            return None

        
    def _resolve_redirects(self, initial_url: str) -> str:
        """Manually resolve redirects to get the final download URL."""
        
        url = initial_url
        ## Loop limit to prevent infinite redirects
        for _ in range(10):
            req = fetches.Fetch(url, headers=self.headers, allow_redirects=False).fetch_req()
            if req and req.status_code in (301, 302, 303, 307):
                url = req.headers["Location"]
            else:
                break
        return url

    
    def run(self):
        """Execute the query and generate download links."""
        
        if not self.aoi or not self.access_token:
            return self

        ## Build OData Filter
        filters = [
            f"Collection/Name eq '{self.collection_name}'",
            f"Attributes/OData.CSC.StringAttribute/any(att:att/Name eq 'productType' and att/OData.CSC.StringAttribute/Value eq '{self.product_type}')",
            f"OData.CSC.Intersects(area=geography'SRID=4326;{self.aoi}')"
        ]

        if self.time_start:
            filters.append(f"ContentDate/Start gt {self.time_start}")
        if self.time_end:
            filters.append(f"ContentDate/Start lt {self.time_end}")
        if self.max_cloud_cover is not None:
            filters.append(f"Attributes/OData.CSC.DoubleAttribute/any(att:att/Name eq 'cloudCover' and att/OData.CSC.DoubleAttribute/Value le {self.max_cloud_cover})")

        query_url = f"{CDSE_CATALOGUE_URL}/Products?$filter={' and '.join(filters)}"
        utils.echo_msg(query_url)
        try:
            response = requests.get(query_url).json()
        except Exception as e:
            utils.echo_error_msg(f"Error querying CDSE Catalogue: {e}")
            return self

        results = response.get('value', [])
        
        with utils.ccp(total=len(results), desc='Parsing datasets...', leave=self.verbose) as pbar:
            for result in results:
                pbar.update()
                try:
                    product_id = result['Id']
                    product_name = result['Name']
                    
                    ## Fetch Metadata XML to determine file structure
                    ## Note: S2MSI1C usually has a specific structure node path
                    meta_url = f"{CDSE_CATALOGUE_URL}/Products({product_id})/Nodes({product_name})/Nodes(MTD_MSIL1C.xml)/$value"
                    
                    ## Resolve initial download URL redirect
                    final_meta_url = self._resolve_redirects(meta_url)
                    
                    ## Fetch content
                    meta_req = fetches.Fetch(final_meta_url, headers=self.headers, verify=True, allow_redirects=True).fetch_req()
                    if not meta_req:
                        continue

                    root = ET.fromstring(meta_req.text.encode('utf-8'))

                    ## Parse Band Locations (Specific to Sentinel-2 S2MSI1C structure)
                    ## Navigating XML: root -> Granule List -> Granule -> IMAGE_FILE (Band 2, 3, 4 typically)
                    ## Adjust indices based on known XML structure safety
                    granule_list = root.find(".//GranuleList") 
                    if granule_list is not None and len(granule_list) > 0:

                        ## We look for Bands 2, 3, 4 (Blue, Green, Red) usually
                        image_files = root.findall(".//IMAGE_FILE")
                        ## Filter for B02, B03, B04
                        target_bands = [img.text for img in image_files if img.text.endswith('B02') or img.text.endswith('B03') or img.text.endswith('B04')]

                        for band_path in target_bands:
                            ## Path in XML is like: GRANULE/L1C_T.../IMG_DATA/T..._B02
                            ## CDSE Node structure requires traversing the path components
                            parts = band_path.split('/')
                            ## Construct Node URL: Nodes(GRANULE)/Nodes(L1C...)/Nodes(IMG_DATA)/Nodes(File)
                            node_path = "/".join([f"Nodes({p})" for p in parts])
                             
                            file_url = f"{CDSE_CATALOGUE_URL}/Products({product_id})/Nodes({product_name})/{node_path}.jp2/$value"
                             
                            self.add_entry_to_results(
                                file_url,
                                f"{parts[-1]}.jp2",
                                'SENTINEL 2LA'
                            )

                except Exception as e:
                    if self.verbose:
                        utils.echo_warning_msg(f"Error parsing result {result.get('Name', 'unknown')}: {e}")
                    continue

        return self

class Sentinel2(CDSE):
    def __init__(self, **kwargs):
        super().__init__(collection_name='SENTINEL-2', product_type='S2MSI1C', **kwargs)

### End
