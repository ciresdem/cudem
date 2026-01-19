### dav.py
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
##
## dav.py is part of CUDEM
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
## Fetch NOAA Lidar/Raster data via the Digital Coast Data Access Viewer (DAV) API v1.
##
### Code:

import os
import json
import requests
from typing import List, Dict, Optional, Any
from osgeo import ogr

from cudem import utils
from cudem import gdalfun
from cudem.fetches import fetches

## ==============================================
## Constants
## ==============================================
DAV_API_URL = 'https://coast.noaa.gov/dataviewer/api/v1/search/missions'
DAV_HEADERS = {'Content-Type': 'application/json'}


# @fetches.cli_opts(
#     help_text="NOAA's Digital Coast",
#     datatype="The desired data type, e.g. 'DEM', 'lidar', etc.",
#     want_footprints="Fetch dataset footprints",
#     footprints_only="Only fetch dataset footprints",
#     keep_footprints="Keep residual footprints",
#     title_filter="Filter the datasets 'title'"
# )

## ==============================================
## DAV Module
## ==============================================
class DAV(fetches.FetchModule):
    """Fetch NOAA lidar/raster data from DAV.

    Uses Digital Coast's Data Access Viewer API to discover
    dataset footprints and download tile indices.

    https://coast.noaa.gov

    Configuration Example:
    < digital_coast:datatype='Lidar':footprints_only=False >
    """
    
    def __init__(
            self,
            datatype: Optional[str] = None,
            want_footprints: bool = False,
            keep_footprints: bool = False,
            footprints_only: bool = False,
            title_filter: Optional[str] = None,
            name: Optional[str] = 'digital_coast',
            **kwargs
    ):
        super().__init__(name=name, **kwargs)
        self.datatype = datatype
        self.title_filter = title_filter
        
        self.want_footprints = want_footprints
        self.keep_footprints = keep_footprints
        self.footprints_only = footprints_only
        
        if self.footprints_only:
            self.want_footprints = True

            
    def _region_to_ewkt(self):
        """Convert the current region to NAD83 (SRID 4269) EWKT Polygon string."""
        
        if self.region is None:
            return None
            
        ## Ensure format is minx, maxx, miny, maxy
        ## Construct closed polygon: minx miny, maxx miny, maxx maxy, minx maxy, minx miny
        w = self.region.xmin
        e = self.region.xmax
        s = self.region.ymin
        n = self.region.ymax
        
        ## WKT format
        poly = f"POLYGON(({w} {s}, {e} {s}, {e} {n}, {w} {n}, {w} {s}))"
        return f"SRID=4269;{poly}"

    
    def _get_features(self) -> List[Dict]:
        """Query the DAV API for missions in the region."""
        
        if self.region is None:
            return []

        ## Map 'datatype' arg to API DataTypes
        ## API expects list: ["Lidar"], ["Imagery"], ["Elevation"], ["Land Cover"], ['DEM']
        dt_map = {
            'lidar': 'Lidar',
            'raster': 'DEM',
            'imagery': 'Imagery',
            'elevation': 'Elevation',
            'dem': 'DEM',
        }
        
        req_types = []
        if self.datatype:
            clean_dt = self.datatype.lower()
            if clean_dt in dt_map:
                req_types.append(dt_map[clean_dt])
            else:
                ## Fallback or pass through if user knows specific API key
                req_types.append(self.datatype)
        else:
            ## Fetching Lidar and Elevation by default for 'elevation' contexts
            req_types = ["Lidar", "Elevation"]

        payload = {
            "aoi": self._region_to_ewkt(),
            "published": "true",
            "dataTypes": req_types
        }

        ## Debug
        #utils.echo_msg(f"DAV Payload: {payload}")

        try:
            r = requests.post(DAV_API_URL, json=payload, headers=DAV_HEADERS)
            r.raise_for_status()
            response = r.json()
            ## Debug
            #utils.echo_msg(f"DAV Response: {response}")

            return response.get('data', {})
        except Exception as e:
            utils.echo_error_msg(f"DAV API Query Error: {e}")
            return []

        
    def _find_index_zip(self, bulk_url: str) -> Optional[str]:
        """Find the tile index zip file given the Bulk Download landing page URL."""
        
        ## Fetch the bulk download page (directory listing)
        try:
            page = fetches.Fetch(bulk_url, verbose=False).fetch_html()
        except:
            return None
            
        if page is None:
            return None

        ## Look for 'urllist' text file which contains direct links
        ## This is standard for NOAA DAV HTTP/S3 directories
        txt_links = page.xpath('//a[contains(@href, ".txt")]/@href')
        urllist_link = next((l for l in txt_links if 'urllist' in l), None)
        
        index_zip_url = None

        if urllist_link:
            ## Construct full URL
            if not urllist_link.startswith('http'):
                urllist_link = requests.compat.urljoin(bulk_url, urllist_link)

            ## Download urllist
            local_urllist = os.path.join(self._outdir, os.path.basename(urllist_link))
            if fetches.Fetch(urllist_link, verbose=self.verbose).fetch_file(local_urllist) == 0:
                # Scan file for the zip
                with open(local_urllist, 'r') as f:
                    for line in f:
                        if 'tileindex' in line and 'zip' in line:
                            index_zip_url = line.strip()
                            break
                utils.remove_glob(local_urllist)
        
        ## If no urllist, look for tileindex zip directly in HTML
        if not index_zip_url:
            zip_links = page.xpath('//a[contains(@href, ".zip")]/@href')
            tile_zip = next((l for l in zip_links if 'tileindex' in l), None)
            if tile_zip:
                if not tile_zip.startswith('http'):
                    index_zip_url = requests.compat.urljoin(bulk_url, tile_zip)
                else:
                    index_zip_url = tile_zip

        return index_zip_url

    
    def _process_index_shapefile(self, shp_path: str, dataset_id: str, data_type: str):
        """Parse the downloaded index shapefile and add intersecting tiles to results."""
        
        ## Read .prj if exists to handle projection
        prj_file = shp_path.replace('.shp', '.prj')
        warp_region = self.region.copy()
        
        if os.path.exists(prj_file):
            with open(prj_file, 'r') as f:
                prj_wkt = f.read()
            warp_region.warp(dst_crs=prj_wkt)
            
        ds = ogr.Open(shp_path)
        if not ds:
            return

        layer = ds.GetLayer(0)
        
        ## Common field names in DAV indices
        known_name_fields = ['Name', 'location', 'filename', 'tilename', 'NAME', 'TILE_NAME']
        known_url_fields = ['url', 'URL', 'path', 'link', 'HTTP_LINK', 'URL_Link']

        for feature in layer:
            geom = feature.GetGeometryRef()
            if geom and geom.Intersects(warp_region.export_as_geom()):
                
                tile_name = None
                tile_url = None
                
                ## Get field definitions
                feat_defn = layer.GetLayerDefn()
                field_names = [feat_defn.GetFieldDefn(i).GetName() for i in range(feat_defn.GetFieldCount())]

                ## Find Name
                for f in known_name_fields:
                    if f in field_names:
                        val = feature.GetField(f)
                        if val: tile_name = str(val).strip()
                    if tile_name: break

                ## Find URL
                for f in known_url_fields:
                    if f in field_names:
                        val = feature.GetField(f)
                        if val: tile_url = str(val).strip()
                    if tile_url: break

                if not tile_url or not tile_name:
                    continue

                ## URL Cleanup logic
                if not tile_url.endswith(tile_name):
                     if tile_url.endswith('/'):
                         tile_url += tile_name
                         ## Check if url ends with filename ignoring case/path
                     elif not tile_url.lower().endswith(os.path.basename(tile_name).lower()):
                         ## Often URL is just the dir, append filename
                         ## Use dirname of tile_url to be safe if it includes a partial file. 
                         ## Usually in DAV shapefiles, if it doesn't end in name, it's the dir.
                         tile_url = f"{tile_url.rstrip('/')}/{os.path.basename(tile_name)}"
                
                self.add_entry_to_results(
                    tile_url,
                    os.path.join(str(dataset_id), os.path.basename(tile_url)),
                    data_type
                )

        ds = None

        
    def run(self):
        """Run the DAV fetching module."""
        
        ## Query API
        data = self._get_features()

        #data = features.get('data', {})
        datasets = data.get('datasets', {})
        
        for dataset in datasets:
            ## DEBUG
            #utils.echo_msg(dataset)
            attrs = dataset.get('attributes', {})
            
            fid = attrs.get('id')
            name = attrs.get('title')
            f_datatype = attrs.get('dataType')
            links_list = attrs.get('links', [])

            if self.title_filter:
                if not self.title_filter.lower() in name.lower():
                    continue
            
            ## Find 'Bulk Download' Link (Service ID 46)
            bulk_url = None
            for link_obj in links_list:
                ## linkTypeId is a string in JSON ("46")
                if link_obj.get('linkTypeId') == "46":
                    bulk_url = link_obj.get('uri')
                    break
            
            if not bulk_url:
                if self.verbose:
                    utils.echo_msg(f"No bulk download found for {name} (ID: {fid})")
                continue

            ## Find Index ZIP from Bulk Page
            index_zip_url = self._find_index_zip(bulk_url)
            
            if not index_zip_url:
                utils.echo_warning_msg(f"Could not locate Tile Index ZIP for {name} at {bulk_url}")
                continue

            ## Handle Footprints Request
            if self.want_footprints:
                self.add_entry_to_results(
                    index_zip_url,
                    os.path.join(str(fid), os.path.basename(index_zip_url)),
                    'footprint'
                )
                if self.footprints_only:
                    continue

            ## Process Shapefile
            surv_name = f"dav_{fid}"
            local_zip = os.path.join(self._outdir, f'tileindex_{surv_name}.zip')
            
            try:
                ## Download
                if fetches.Fetch(index_zip_url, verbose=self.verbose).fetch_file(local_zip) == 0:
                    
                    ## Unzip
                    unzipped = utils.p_unzip(local_zip, ['shp', 'shx', 'dbf', 'prj'], outdir=self._outdir, verbose=self.verbose)
                    shp_file = next((f for f in unzipped if f.endswith('.shp')), None)
                    
                    if shp_file:
                        if self.verbose:
                            utils.echo_msg(f"Processing index: {shp_file}")
                        self._process_index_shapefile(shp_file, fid, f_datatype)
                    
                    ## Cleanup
                    if not self.keep_footprints:
                        utils.remove_glob(local_zip, *unzipped)
                else:
                    utils.echo_warning_msg(f"Failed to download index: {index_zip_url}")
                        
            except Exception as e:
                utils.echo_error_msg(f"Error processing DAV dataset {fid}: {e}")

        return self


## ==============================================
## Subclasses / Shortcuts
## ==============================================
class SLR(DAV):
    """Sea Level Rise DEMs via Digital Coast."""
    
    def __init__(self, **kwargs):
        super().__init__(name='SLR', title_filter='SLR', datatype='DEM', **kwargs)

        
class CoNED(DAV):
    """Coastal NED (CoNED) DEMs via Digital Coast."""
    
    def __init__(self, **kwargs):
        super().__init__(name='CoNED', title_filter='CoNED', datatype='DEM', **kwargs)

        
class CUDEM(DAV):
    """CUDEM Tiled DEMs via Digital Coast."""
    
    def __init__(self, datatype: Optional[str] = None, **kwargs):
        super().__init__(name='CUDEM', datatype='DEM', title_filter='CUDEM', **kwargs)

        
### End
