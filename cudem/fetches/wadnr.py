### wadnr.py
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
##
## wadnr.py is part of CUDEM
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
## Fetch Washington State Department of Natural Resources (WA DNR) LiDAR data.
##
### Code:

import json
from typing import Optional, List, Any, Dict
from cudem import regions
from cudem.fetches import fetches

## ==============================================
## Constants
## ==============================================
WA_DNR_BASE_URL = "https://lidarportal.dnr.wa.gov"
WA_DNR_DOWNLOAD_URL = f"{WA_DNR_BASE_URL}/download?"
WA_DNR_REST_URL = f"{WA_DNR_BASE_URL}/arcgis/rest/services/lidar/wadnr_hillshade/MapServer"
WA_DNR_LAYERS_URL = f"{WA_DNR_REST_URL}/layers?f=pjson"

## ==============================================
## WaDNR Module
## ==============================================
class WADNR(fetches.FetchModule):
    """Washington State Department of Natural Resources lidar data.
    
    Fetches LiDAR data from the WA DNR portal based on region intersection.
    """
    
    def __init__(self, ids: Optional[str] = None, **kwargs):
        super().__init__(name='waDNR', **kwargs) 

        ## Metadata
        self.data_format = -2 # zip file
        self.src_srs = None
        self.title = 'Washington DNR Lidar'
        self.source = 'Washington DNR'
        self.data_type = 'lidar'
        self.url = WA_DNR_DOWNLOAD_URL
        
        ## Headers (Browser emulation often required)
        self.headers = {
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:89.0) Gecko/20100101 Firefox/89.0'
        }

        ## ID filtering
        self.ids = [int(x) for x in ids.split('/')] if ids else []

    def find_ids(self) -> List[List[Any]]:
        """Query the map server to find layer IDs intersecting the region."""
        
        ## Fetch layer metadata
        req = fetches.Fetch(WA_DNR_LAYERS_URL, verbose=self.verbose).fetch_req(tries=10, timeout=5)
        if req is None:
            return []

        try:
            layers_json = req.json()
        except json.JSONDecodeError:
            return []

        layers_in_region = []
        
        ## WA DNR layers are typically in Web Mercator (EPSG:3857)
        ## We need to warp the extents to EPSG:4326 to check against self.region        
        for layer in layers_json.get('layers', []):
            layer_name = layer.get('name')
            sub_layers = layer.get('subLayers', [])
            
            ## Extract extent
            extent = layer.get('extent')
            if not extent:
                continue

            ## Construct region from extent and warp to WGS84
            try:
                layer_region = regions.Region(src_srs='epsg:3857').from_list([
                    extent['xmin'], extent['xmax'],
                    extent['ymin'], extent['ymax']
                ]).warp('epsg:4326')
            except Exception:
                continue

            ## Check Intersection and Extract IDs
            if sub_layers and regions.regions_intersect_ogr_p(layer_region, self.region):
                try:
                    ## WA DNR specific ID logic:
                    valid_ids = []
                    for sub in sub_layers:
                        try:
                            ## Attempt to parse ID from sublayer name
                            parsed_id = int(sub['name'][:-1]) - 1
                            valid_ids.append(parsed_id)
                        except (ValueError, TypeError):
                            continue
                    
                    if valid_ids:
                        layers_in_region.append([layer_name, min(valid_ids)])

                except Exception:
                    continue

        ## Filter by user-provided IDs if present
        if self.ids:
            layers_in_region = [x for x in layers_in_region if x[1] in self.ids]
            
        return layers_in_region

    
    def run(self):
        """Run the WA DNR fetching module."""

        if self.region is None:
            return []

        layers_in_region = self.find_ids()        
        
        for layer_name, layer_id in layers_in_region:
            
            ## Construct payload
            ## WA DNR download endpoint expects a polygon and a list of IDs
            payload = {
                'ids': [layer_id],
                'format': 'json',
                'geojson': json.dumps({
                    'type': 'Polygon',
                    'coordinates': [self.region.export_as_polygon()]
                })
            }
            
            ## Fetch the download link
            ## The endpoint usually returns a redirect or a JSON with the file link
            req = fetches.Fetch(
                WA_DNR_DOWNLOAD_URL, 
                verbose=self.verbose
            ).fetch_req(params=payload, tries=10, timeout=10)

            if req is not None and req.status_code == 200:
                out_fn = f"{layer_name}_{layer_id}.zip"
                self.add_entry_to_results(req.url, out_fn, 'wa_dnr')
                
        return self

    
### End
