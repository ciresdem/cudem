### cdse.py
##
## Copyright (c) 2025 - 2026 Regents of the University of Colorado
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
import time
import re
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

    Automatically refreshes authentication tokens upon expiration.
    """

    def __init__(self, collection_name: str = 'SENTINEL-2', product_type: str = 'S2MSI1C',
                 max_cloud_cover: Optional[float] = None, time_start: str = '', time_end: str = '', **kwargs):
        # Initialize private headers attribute before super().__init__ calls the setter
        self._headers = {}
        self._token_expiry = 0.0

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
        ## Default to One year ago -> Today
        if not time_end:
            time_end = datetime.datetime.now().isoformat()
        if not time_start:
            time_start = (datetime.datetime.now() - datetime.timedelta(days=365)).isoformat()

        self.time_start = self._format_date(time_start)
        self.time_end = self._format_date(time_end)

        ## Initial Authentication
        self.access_token = self.refresh_token()


    @property
    def headers(self):
        """Dynamic headers property that checks for token expiry."""

        ## Buffer of 30 seconds to ensure we don't expire mid-request
        if time.time() > (self._token_expiry - 30):
            if self.verbose:
                utils.echo_msg("CDSE Token expired or expiring soon. Refreshing...")
            self.refresh_token()
        return self._headers


    @headers.setter
    def headers(self, value):
        self._headers = value


    def refresh_token(self):
        """Acquire a new token and update headers/expiry."""

        token, expires_in = self.get_auth_token_data()

        if token:
            self.access_token = token
            ## Set expiry time (current time + lifetime)
            self._token_expiry = time.time() + expires_in
            self._headers = {'Authorization': f'Bearer {token}'}
        else:
            self._headers = {}
            utils.echo_warning_msg("Could not acquire CDSE Access Token. Requests may fail.")
        return token


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


    def get_auth_token_data(self) -> Tuple[Optional[str], int]:
        """Authenticate with CDSE and retrieve access token AND expiration time."""

        username, password = self.get_credentials_from_netrc()

        if not username or not password:
            utils.echo_warning_msg("No credentials found in .netrc for CDSE.")
            return None, 0

        data = {
            "client_id": "cdse-public",
            "grant_type": "password",
            "username": username,
            "password": password,
        }

        try:
            response = requests.post(CDSE_AUTH_URL, data=data)
            response.raise_for_status()
            json_resp = response.json()
            token = json_resp.get("access_token")
            ## Default to 600s (10 mins) if not provided by server
            expires_in = json_resp.get("expires_in", 600)

            if self.verbose:
                utils.echo_msg(f"Successfully retrieved CDSE access token (expires in {expires_in}s).")
            return token, int(expires_in)
        except requests.exceptions.RequestException as e:
            utils.echo_error_msg(f"CDSE Authentication failed: {e}")
            return None, 0


    def _resolve_redirects(self, initial_url: str) -> str:
        """Manually resolve redirects to get the final download URL."""

        url = initial_url
        for _ in range(10):
            # Accessing self.headers property here ensures token is fresh for meta requests too
            req = fetches.Fetch(url, headers=self.headers, allow_redirects=False).fetch_req()
            if req and req.status_code in (301, 302, 303, 307):
                url = req.headers["Location"]
            else:
                break
        return url


    def _strip_ns(self, xml_text):
        """Strip namespaces from XML string to allow simple tag searching."""

        return re.sub(r'\sxmlns="[^"]+"', '', xml_text, count=1)


    def run(self):
        """Execute the query and generate download links."""

        if not self.aoi or not self.access_token:
            return self

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

        ## Initial query URL
        query_url = f"{CDSE_CATALOGUE_URL}/Products?$filter={' and '.join(filters)}&$top=100"

        if self.verbose:
            utils.echo_msg(f"Querying CDSE: {query_url}")

        page_count = 0

        ## Pagination loop: continues as long as a next link exists
        while query_url:
            page_count += 1
            if self.verbose:
                utils.echo_msg(f"Fetching page {page_count}...")

            try:
                response_req = requests.get(query_url)
                response_req.raise_for_status()
                response = response_req.json()
            except Exception as e:
                utils.echo_error_msg(f"Error querying CDSE Catalogue: {e}")
                break

            results = response.get('value', [])

            if not results:
                break

            ## Process current page of results
            with utils.ccp(total=len(results), desc=f'Parsing page {page_count}...', leave=self.verbose) as pbar:
                for result in results:
                    pbar.update()
                    try:
                        product_id = result['Id']
                        product_name = result['Name']

                        meta_url = f"{CDSE_CATALOGUE_URL}/Products({product_id})/Nodes({product_name})/Nodes(MTD_MSIL1C.xml)/$value"
                        final_meta_url = self._resolve_redirects(meta_url)

                        meta_req = fetches.Fetch(final_meta_url, headers=self.headers, verify=True, allow_redirects=True).fetch_req()
                        if not meta_req:
                            continue

                        xml_content = self._strip_ns(meta_req.text)
                        root = ET.fromstring(xml_content.encode('utf-8'))

                        image_files = root.findall(".//IMAGE_FILE")
                        target_bands = [img.text for img in image_files if img.text and (img.text.endswith('B02') or img.text.endswith('B03') or img.text.endswith('B04'))]

                        for band_path in target_bands:
                            parts = band_path.split('/')
                            parts[-1] = f"{parts[-1]}.jp2"
                            node_path = "/".join([f"Nodes({p})" for p in parts])
                            file_url = f"{CDSE_CATALOGUE_URL}/Products({product_id})/Nodes({product_name})/{node_path}/$value"

                            self.add_entry_to_results(
                                file_url,
                                parts[-1],
                                'SENTINEL 2LA'
                            )

                    except Exception as e:
                        if self.verbose:
                            utils.echo_warning_msg(f"Error parsing result {result.get('Name', 'unknown')}: {e}")
                        continue

            ## Check for the next page link
            query_url = response.get('@odata.nextLink', None)

        return self


class Sentinel2(CDSE):
    def __init__(self, **kwargs):
        super().__init__(collection_name='SENTINEL-2', product_type='S2MSI1C', **kwargs)

### End
