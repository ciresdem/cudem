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
## Fetch NOAA Lidar/Raster data via the Digital Coast Data Access Viewer (DAV).
##
### Code:

import os
import json
from typing import List, Dict, Optional, Any
from osgeo import ogr

from cudem import utils
from cudem import gdalfun
from cudem import vdatums
from cudem.fetches import fetches

## ==============================================
## Constants
## ==============================================
DAV_BASE_URL = 'https://maps.coast.noaa.gov/arcgis/rest/services/DAV/DAV_footprints/MapServer/'

## ==============================================
## DAV Module
## ==============================================
class DAV(fetches.FetchModule):
    """Fetch NOAA lidar/raster data from DAV.

    Uses Digital Coast's Data Access Viewer Mapserver to discover
    dataset footprints and download tile indices.

    https://coast.noaa.gov

    Configuration Example:
    < digital_coast:where='1=1':datatype=None:footprints_only=False >
    """
    
    def __init__(
            self,
            where: str = '1=1',
            index: bool = False,
            datatype: Optional[str] = None,
            layer: int = 0,
            name: str = 'digital_coast',
            want_footprints: bool = False,
            keep_footprints: bool = False,
            footprints_only: bool = False,
            **kwargs
    ):
        super().__init__(name=name, **kwargs)
        self.where = where
        self.index = index
        self.datatype = datatype
        self.layer = utils.int_or(layer)
        
        self.want_footprints = want_footprints
        self.keep_footprints = keep_footprints
        self.footprints_only = footprints_only
        
        if self.footprints_only:
            self.want_footprints = True

        self._dav_api_url = f"{DAV_BASE_URL}{self.layer}/query?"

        
    def _get_features(self) -> List[Dict]:
        """Query the ArcGIS REST endpoint for available datasets in the region."""
        if self.region is None:
            return []

        ## Construct Where Clause with Date Filters
        ## self.min_year/max_year are populated by the parent FetchModule
        ## if passed in kwargs
        current_where = self.where
        
        if self.min_year is not None:
            current_where += f" AND Year >= {self.min_year}"
            
        if self.max_year is not None:
            current_where += f" AND Year <= {self.max_year}"

        params = {
            'where': current_where,
            'outFields': '*',
            'geometry': self.region.format('bbox'),
            'inSR': 4326,
            'f': 'pjson',
            'returnGeometry': 'False',
        }
        
        req = fetches.Fetch(self._dav_api_url, verbose=self.verbose).fetch_req(params=params)
        if req is None:
            return []

        try:
            data = req.json()
            if 'error' in data:
                utils.echo_error_msg(f"DAV Query Error: {data['error']}, {req.url}")
                return []
            return data.get('features', [])
        except json.JSONDecodeError:
            utils.echo_error_msg("Failed to parse DAV response.")
            return []

        
    def _parse_metadata_for_index(self, metadata_url: str, feature_id: str) -> Optional[str]:
        """Scrape the metadata page to find the tile index zip file URL."""
        
        page = fetches.Fetch(metadata_url, verbose=False).fetch_html()
        if page is None:
            return None

        ## Look for index.html links (common in bulk download pages)
        index_urls = page.xpath('//a[contains(@href, "index.html")]/@href')
        
        ## Fallback: Look for ID-based directories if no index.html found
        if not index_urls:
            index_urls = page.xpath(f'//a[contains(@href, "_{feature_id}/")]/@href')

        if not index_urls:
            return None
        
        ## We assume the first found link is the correct directory
        base_url = index_urls[0]
        
        ## Fetch the directory page to find the text file list (urllist)
        dir_page = fetches.Fetch(base_url, verbose=False).fetch_html()
        if dir_page is None:
            return None
            
        txt_links = dir_page.xpath('//a[contains(@href, ".txt")]/@href')
        urllist_link = next((l for l in txt_links if 'urllist' in l), None)
        
        if not urllist_link:
            return None
            
        ## Construct full URL for the text list
        if 'http' in urllist_link:
            full_urllist_url = urllist_link
        else:
            ## Handle relative paths
            parent_url = os.path.dirname(base_url)
            full_urllist_url = f"{parent_url}/{urllist_link}" if parent_url.endswith('/') else f"{parent_url}/{urllist_link}"
            ## Sometimes base_url is the dir itself
            if not full_urllist_url.startswith('http'):
                #full_urllist_url = f"{base_url}/{urllist_link}".replace('//', '/')
                full_urllist_url =  requests.compat.urljoin(base_url, urllist_link)

        ## Download and read the urllist file to find the zip
        local_urllist = os.path.join(self._outdir, os.path.basename(full_urllist_url))
        if fetches.Fetch(full_urllist_url, verbose=self.verbose).fetch_file(local_urllist) != 0:
            return None
            
        index_zip_url = None
        if os.path.exists(local_urllist):
            with open(local_urllist, 'r') as f:
                for line in f:
                    if 'tileindex' in line and 'zip' in line:
                        index_zip_url = line.strip()
                        break
            utils.remove_glob(local_urllist)
            
        return index_zip_url

    
    def _process_index_shapefile(self, shp_path: str, dataset_id: str, data_type: str):
        """Parse the downloaded index shapefile and add intersecting tiles to results."""
        
        ## Handle Projection (Read .prj if exists)
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
        
        ## Standardize field names often found in these indices
        known_name_fields = ['Name', 'location', 'filename', 'tilename']
        known_url_fields = ['url', 'URL', 'path', 'link']

        for feature in layer:
            geom = feature.GetGeometryRef()
            if geom and geom.Intersects(warp_region.export_as_geom()):
                
                ## Extract Name and URL
                tile_name = None
                tile_url = None
                field_names = [f.name for f in layer.schema]

                for f in known_name_fields:
                    if f in field_names:
                        val = feature.GetField(f)
                        if val: tile_name = val.strip()
                    if tile_name: break

                for f in known_url_fields:
                    if f in field_names:
                        val = feature.GetField(f)
                        if val: tile_url = val.strip()
                    if tile_url: break

                if not tile_url or not tile_name:
                    continue

                ## Clean URL construction
                ## Often the URL field is just the base path, and we append the filename
                if not tile_url.endswith(tile_name):
                     if tile_url.endswith('/'):
                         tile_url += tile_name
                     elif not tile_url.lower().endswith(tile_name.lower().split('/')[-1]):
                         base = '/'.join(tile_url.split('/')[:-1])
                         tile_url = f"{base}/{os.path.basename(tile_name)}"
                
                self.add_entry_to_results(
                    tile_url,
                    os.path.join(str(dataset_id), os.path.basename(tile_url)),
                    data_type
                )

        ds = None

        
    def run(self):
        """Run the DAV fetching module."""
        
        features = self._get_features()
        
        for feature in features:
            attrs = feature.get('attributes', {})
            
            ## Filter by Datatype
            f_datatype = attrs.get('data_type', '')
            if self.datatype and self.datatype.lower() != f_datatype.lower():
                if self.datatype.lower() != 'sm':
                    continue

            fid = attrs.get('id')
            meta_link = attrs.get('metadata')
            
            if not meta_link or not fid:
                continue

            ## Find the bulk download index zip
            index_zip_url = self._parse_metadata_for_index(meta_link, fid)
            
            if not index_zip_url:
                continue

            ## Handle Footprints
            if self.want_footprints:
                self.add_entry_to_results(
                    index_zip_url,
                    os.path.join(str(fid), os.path.basename(index_zip_url)),
                    'footprint'
                )
                if self.footprints_only:
                    continue

            ## Download and Process Index Shapefile
            surv_name = f"dav_{fid}"
            local_zip = os.path.join(self._outdir, f'tileindex_{surv_name}.zip')
            
            try:
                if fetches.Fetch(index_zip_url, verbose=self.verbose).fetch_file(local_zip) == 0:
                    
                    unzipped = utils.p_unzip(local_zip, ['shp', 'shx', 'dbf', 'prj'], outdir=self._outdir)
                    shp_file = next((f for f in unzipped if f.endswith('.shp')), None)
                    
                    if shp_file:
                        self._process_index_shapefile(shp_file, fid, f_datatype)
                    
                    ## Cleanup
                    if not self.keep_footprints:
                        utils.remove_glob(local_zip, *unzipped)
                        
            except Exception as e:
                utils.echo_error_msg(f"Error processing DAV dataset {fid}: {e}")

        return self


## ==============================================
## Subclasses / Shortcuts
## ==============================================
class SLR(DAV):
    """Sea Level Rise DEMs via Digital Coast."""
    
    def __init__(self, **kwargs):
        super().__init__(name='SLR', where='ID=6230', **kwargs)

        
class CoNED(DAV):
    """Coastal NED (CoNED) DEMs via Digital Coast."""
    
    def __init__(self, **kwargs):
        super().__init__(name='CoNED', where="title LIKE '%CoNED%'", **kwargs)

        
class CUDEM(DAV):
    """CUDEM Tiled DEMs via Digital Coast."""
    
    def __init__(self, datatype: Optional[str] = None, **kwargs):
        datatype = utils.str_or(datatype, 'all')
        if datatype == '19' or datatype.lower() == 'ninth':
            where = "title LIKE '%CUDEM%Ninth%'"
        elif datatype == '13' or datatype.lower() == 'third':
            where = "title LIKE '%CUDEM%Third%'"
        else:
            where = "title LIKE '%CUDEM%'"
        
        super().__init__(name='CUDEM', where=where, **kwargs)

        
### End
