### osm.py
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
##
## osm.py is part of CUDEM
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
## Fetch and process OpenStreetMap (OSM) data.
## Supports fetching raw OSM XML via Overpass API and processing coastlines/water bodies into polygons.
##
### Code:

import os
import sys
import argparse
import numpy as np
from osgeo import ogr
from typing import Optional, List, Union

from cudem import utils
from cudem import regions
from cudem import gdalfun
from cudem import xyzfun
from cudem import regions
from cudem.fetches import fetches
from cudem.fetches import gmrt

## ==============================================
## Constants
## ==============================================
OSM_API_MAIN = 'https://lz4.overpass-api.de/api/interpreter'
OSM_API_MIRRORS = [
    'https://overpass.kumi.systems/api/interpreter',
    'https://overpass.openstreetmap.fr/api/interpreter'
]
OSM_PLANET_PBF = 'https://ftpmirror.your.org/pub/openstreetmap/pbf/planet-latest.osm.pbf'

import os
import numpy as np
from osgeo import ogr
from cudem import utils
from cudem import regions
from cudem import gdalfun
from cudem import xyzfun
from cudem.fetches import gmrt

## ==============================================
## Classes
## ==============================================
class OpenStreetMap(fetches.FetchModule):
    """OpenStreetMap data fetcher.
    
    Query the Overpass API for specific tags.
    
    Configuration Example:
    < osm:q=None:fmt=osm:planet=False:chunks=True:min_length=None >
    """
    
    def __init__(self, q=None, fmt='osm', planet=False, chunks=False,
                 min_length=None, **kwargs):
        super().__init__(name='osm', **kwargs)
        self.q = q
        self.q_fn = q if q else 'custom'
        self.fmt = fmt
        self.planet = planet
        self.chunks = chunks
        
        ## Define Standard Queries
        self.timeout_header = '[timeout:3600]'
        self.query_body = None

        if q == 'buildings':
            self.query_body = '''
            (way["building"]{};
            relation["building"]["type"="multipolygon"];
            );
            (._;>;);
            out meta;
            '''
            
        elif q == 'coastline':
            self.timeout_header = '[timeout:32]' # Shorter timeout for coastline chunks
            self.query_body = '''
            (way["natural"="coastline"]{};
            relation["type"="lines"];
            );
            (._;>;);
            out meta;
            '''

        elif q == 'water':
            self.query_body = '''
            (way["natural"="water"]{};
            rel["natural"="water"];
            );
            (._;>;);
            out meta;
            '''
        
        elif q == 'rivers':
            self.query_body = '''
            (way["natural"="water"]{};
            rel["natural"="river"];
            );
            (._;>;);
            out meta;
            '''
            
        elif q == 'lakes':
             self.query_body = '''
            (way["natural"="water"]{};
            rel["natural"="water"];
            rel["water"="lakes"];
            );
            (._;>;);
            out meta;
            '''
        
        ## Apply min_length filter if applicable
        filter_str = ''
        if min_length is not None:
            filter_str = f'(if: length() > {min_length})'
            
        if self.query_body:
            self.query_body = self.query_body.format(filter_str)
        else:
            self.query_body = q # Allow custom query string input

        self.headers = { 
            'User-Agent': 'CUDEM/OSM Fetcher',
            'referer': 'https://lz4.overpass-api.de/' 
        }

        
    def run(self):
        """Run the OSM fetches module."""
        
        if self.region is None:
            return []

        ## Case 1: Fetch Whole Planet
        if self.planet:
            self.add_entry_to_results(
                OSM_PLANET_PBF,
                os.path.join(self._outdir, 'planet-latest.osm.pbf'),
                'pbf'
            )
            return self

        ## Case 2: Fetch by Region (Chunked or Whole)
        if self.chunks:
            ## Divide region into manageable chunks for Overpass API
            x_delta = self.region.xmax - self.region.xmin
            y_delta = self.region.ymax - self.region.ymin
            
            ## Default to .25 degree chunks if region is large
            n_chunk = None
            if x_delta > 0.25 or y_delta > 0.25:
                ## Basic calculation to split into approx .25 deg chunks
                incs = self.region.increments(1000, 1000) # Dummy incs to init geotransform
                xcount, ycount, _ = self.region.geo_transform(x_inc=incs[0], y_inc=incs[1])
                
                pass 

            ## Create Chunks
            these_regions = self.region.chunk(0.1) # Hardcoded reasonable chunk size for Overpass
            utils.echo_msg(f'Chunking OSM request into {len(these_regions)} regions')
        else:
            these_regions = [self.region]

        for this_region in these_regions:
            c_bbox = this_region.format('osm_bbox')
            out_fn = f'osm_{self.q_fn}_{this_region.format("fn_full")}'
            
            ## Construct Overpass QL
            bbox_filter = f'[bbox:{c_bbox}]'
            output_format = f'[out:{self.fmt}]' if self.fmt != 'osm' else ''
            
            ## Combine headers and body
            full_query = f"{self.timeout_header}{output_format}{bbox_filter};\n{self.query_body}"
            
            ## Encode
            params = fetches.urlencode({'data': full_query})
            url = f"{OSM_API_MAIN}?{params}"
            
            self.add_entry_to_results(
                url,
                f'{out_fn}.{self.fmt}',
                'osm'
            )

        return self

## ==============================================
## OSM Coastline Polygonizer
## ==============================================
class OSMCoastlinePolygonizer:
    """
    A class to polygonize OSM coastline data (lines and polygons) into
    classified Land/Water polygons within a specific region.
    """

    def __init__(self, src_ogr, dst_ogr, region=None, include_landmask=True,
                 landmask_is_watermask=False, line_buffer=0.0000001, verbose=True):
        """Initialize the Polygonizer.

        Args:
            src_ogr (str): Path to source OGR dataset (OSM data).
            dst_ogr (str): Path to destination OGR dataset.
            region (regions.Region, optional): Region to clip/process. Defaults to source extent.
            include_landmask (bool): If True, output land polygons (watermask=0).
            landmask_is_watermask (bool): If True, flips the logic (Land=1, Water=0).
            line_buffer (float): Buffer distance for lines to ensure clean intersections.
            verbose (bool): Print progress messages.
        """
        
        self.src_ogr = src_ogr
        self.dst_ogr = dst_ogr
        self.region = region
        self.include_landmask = include_landmask
        self.landmask_is_watermask = landmask_is_watermask
        self.line_buffer = line_buffer
        self.verbose = verbose

        self.ds_src = None
        self.ds_dst = None
        self.layer_out = None
        self.region_geom = None
        self.has_feature = False

        
    def _determine_polygon_side(self, split_geom, line_geometries):
        """Determine if a split polygon lies on the 'land' or 'water' side of the coastline.
        Returns True if the polygon is likely water, False otherwise.
        """
        
        ss = []
        for line_geometry in line_geometries:
            ## Only check lines that actually touch this split polygon
            if not split_geom.Intersects(line_geometry.Buffer(self.line_buffer)):
                continue

            point_count = line_geometry.GetPointCount()
            for point_n in range(0, point_count - 1):
                ## Extract segment coordinates
                x_beg, y_beg = line_geometry.GetX(point_n), line_geometry.GetY(point_n)
                x_end, y_end = line_geometry.GetX(point_n + 1), line_geometry.GetY(point_n + 1)

                ## Create a test point perpendicular/offset from the segment
                y_ext, x_ext = y_end, x_beg
                
                xyz = xyzfun.XYZPoint(x=x_ext, y=y_ext)
                p_geom = ogr.CreateGeometryFromWkt(xyz.export_as_wkt())

                ## If test point is not within the split polygon, flip the test point
                if not p_geom.Within(split_geom):
                    y_ext, x_ext = y_beg, x_end
                    xyz = xyzfun.XYZPoint(x=x_ext, y=y_ext)
                    p_geom = ogr.CreateGeometryFromWkt(xyz.export_as_wkt())

                ## Check cross product/orientation
                ## s = (x2-x1)*(y_test-y1) > (y2-y1)*(x_test-x1)
                s = (x_end - x_beg) * (y_ext - y_beg) > (y_end - y_beg) * (x_ext - x_beg)
                ss.append(s)

        ## Voting mechanism to determine side
        if all(ss):
            return True
        elif not any(ss):
            return False
        else:
            return np.count_nonzero(ss) > len(ss) / 2

        
    def _write_feature(self, geometry, is_water_side):
        """Write a geometry feature to the output layer with the correct watermask attribute.
        
        - Default: Land=0, Water=1
        - If landmask_is_watermask: Land=1, Water=0
        """
        ## Determine the raw 'watermask' value based on the side calculation
        ## If is_water_side (s=True), normally watermask=1 (Water).
        ## If not is_water_side (s=False), normally watermask=0 (Land).
        
        ## Apply inversion if requested
        if self.landmask_is_watermask:
            is_water_side = not is_water_side

        ## Standard: 0=Land, 1=Water. 
        ## If is_water_side is True, we want 1. If False, we want 0.
        val = 0 if is_water_side else 1 

        ## Filter output based on settings
        ## We write if it is water (val=1) OR if we are including land (val=0)
        if val == 1 or (self.include_landmask and val == 0):
            feature = ogr.Feature(self.layer_out.GetLayerDefn())
            feature.SetGeometry(geometry)
            feature.SetField('watermask', val)
            self.layer_out.CreateFeature(feature)
            feature = None

            
    def _process_lines(self, line_layer):
        """Process LineString features (Coastlines)."""
        
        ## Union all lines into a single geometry collection for cleaner cutting
        line_geometries = gdalfun.ogr_union_geom(
            line_layer,
            ogr.wkbMultiLineString,
            verbose=True
        )
        
        if line_geometries.IsEmpty():
            return

        self.has_feature = True
        
        ## Buffer lines slightly to ensure clean topological operations
        poly_line = line_geometries.Buffer(self.line_buffer)
        
        ## Split the region polygon by the coastline
        split_geoms = self.region_geom.Difference(poly_line)

        ## Iterate over the resulting polygons (land chunks and water chunks)
        for split_geom in split_geoms:
            ## Determine if this chunk is Land or Water
            is_land_side = self._determine_polygon_side(split_geom, line_geometries)
            self._write_feature(split_geom, is_land_side)

            
    def _process_polygons(self, line_layer):
        """Process Polygon features (Islands)."""
        
        for feature in line_layer:
            geometry = feature.geometry()
            geometry = ogr.ForceTo(geometry, ogr.wkbLinearRing)
            
            if geometry.IsEmpty():
                continue

            self.has_feature = True

            ## Handle Overlaps: Subtract this island from existing features in output
            ## (e.g., removing an island hole from a water body)
            for out_feat in self.layer_out:
                out_geom = out_feat.geometry()
                if out_geom.Intersects(geometry):
                    diff_geom = out_geom.Difference(geometry)
                    out_feat.SetGeometry(diff_geom)
                    self.layer_out.SetFeature(out_feat)

            ## Write the island polygon (treated as Land, so is_land_side=True)
            self._write_feature(geometry, is_water_side=True)

            
    def _handle_no_features(self):
        """Handle case where no coastline features were found.
        Determines if the entire region is Land or Water using GMRT.
        """
        
        center_pnt = self.region.center()
        if center_pnt is None:
            return

        ## Fetch elevation at center to guess type
        center_z = utils.int_or(
            gmrt.gmrt_fetch_point(
                latitude=center_pnt[1],
                longitude=center_pnt[0]
            )
        )
        
        ## Positive Z = Land (True), Negative Z = Water (False)
        is_land = center_z >= 0
        self._write_feature(self.region_geom, is_land)

        
    def _setup(self):
        """Setup DataSources and Layers."""
        
        self.ds_src = ogr.Open(self.src_ogr)
        if self.ds_src is None:
            raise IOError(f"Could not open source OGR: {self.src_ogr}")

        ## Setup Region
        src_layer = self.ds_src.GetLayer()
        src_extent = src_layer.GetExtent()
        src_region_obj = regions.Region().from_list(src_extent)
        
        if self.region is not None and self.region.valid_p():
            self.region_geom = self.region.export_as_geom()
        else:
            self.region = src_region_obj
            self.region_geom = self.region.export_as_geom()

        ## Setup Output
        driver = ogr.GetDriverByName("GPKG")
        if os.path.exists(self.dst_ogr):
            self.ds_dst = driver.Open(self.dst_ogr, 1)
            self.layer_out = self.ds_dst.GetLayer()
        else:
            self.ds_dst = driver.CreateDataSource(self.dst_ogr)
            self.layer_out = self.ds_dst.CreateLayer(
                'split_polygons',
                src_layer.GetSpatialRef(),
                ogr.wkbMultiPolygon
            )
            self.layer_out.CreateField(ogr.FieldDefn('watermask', ogr.OFTInteger))

            
    def run(self):
        """Execute the polygonization process."""
        
        self._setup()

        with utils.ccp(total=len(self.ds_src), desc='Polygonizing OSM coastline', leave=self.verbose) as pbar:
            for layer in self.ds_src:
                geom_type = layer.GetGeomType()

                ## Process Lines (Coastlines)
                if geom_type == 2: # wkbLineString
                    self._process_lines(layer)

                ## Process Polygons (Islands)
                elif geom_type == 6: # wkbMultiPolygon / wkbPolygon
                    self._process_polygons(layer)
                
                pbar.update()

        ## Fallback if empty
        if not self.has_feature:
            self._handle_no_features()

        ## Cleanup
        self.ds_src = None
        self.ds_dst = None
        return self.dst_ogr


## ==============================================
## Helper Function
## ==============================================
def polygonize_osm_coastline(src_ogr, dst_ogr, region=None, include_landmask=True,
                             landmask_is_watermask=False, line_buffer=0.0000001, verbose=True):
    """Wrapper function for backward compatibility."""
    
    poly = OSMCoastlinePolygonizer(
        src_ogr, dst_ogr, region, include_landmask,
        landmask_is_watermask, line_buffer, verbose
    )
    return poly.run()


## ==============================================
## OSM Coastline/Water Generator
## ==============================================
class osmCoastline:
    """Wrapper to Fetch and Process OSM Coastline/Water data into Polygons."""
    
    def __init__(
            self,
            region=None,
            chunks=True,
            verbose=True,
            attempts=5,
            n_threads=1,
            cache_dir='.',
            landmask_is_watermask=False,
            line_buffer=0.0000001,
            include_landmask=False,
            q='coastline'
    ):
        self.region = region
        self.chunks = chunks
        self.verbose = verbose
        self.attempts = attempts
        self.n_threads = n_threads
        self.cache_dir = cache_dir
        self.landmask_is_watermask = landmask_is_watermask
        self.line_buffer = line_buffer
        self.include_landmask = include_landmask
        self.q = q
        self.osm_fetcher = None

        
    def init_fetch(self):
        """Initialize the OpenStreetMap fetcher."""
        
        self.osm_fetcher = OpenStreetMap(
            src_region=self.region,
            verbose=self.verbose,
            outdir=self.cache_dir,
            q=self.q,
            chunks=self.chunks,
        )

        
    def fetch(self):
        """Perform the download."""
        
        if self.osm_fetcher is None:
            self.init_fetch()
            
        self.osm_fetcher.run()
        
        ## Execute the fetch queue
        fr = fetches.fetch_results(
            self.osm_fetcher, 
            check_size=False, 
            attempts=self.attempts, 
            n_threads=self.n_threads
        )
        fr.daemon = True
        fr.start()
        fr.join()
        return fr

    
    def process(self, out_fn: str):
        """Process Coastlines into Polygons."""
        
        cst_geoms = []
        if self.osm_fetcher is None:
            self.osm_fetcher = self.fetch()

        if self.osm_fetcher.results:
            with utils.ccp(total=len(self.osm_fetcher.results), desc='Processing Coastlines', leave=self.verbose) as pbar:
                for entry in self.osm_fetcher.results:
                    ## Entry: [url, local_file, type, status]
                    if entry[-1] == 0: # Status 0 = Success
                        local_file = entry[1]
                        
                        ## Generate temporary output name if needed
                        temp_out = out_fn if out_fn else utils.make_temp_fn(
                            f'{utils.fn_basename2(local_file)}_coast.gpkg',
                            temp_dir=self.cache_dir
                        )

                        polygonize_osm_coastline(
                            local_file,
                            temp_out,
                            region=self.region,
                            include_landmask=self.include_landmask,
                            landmask_is_watermask=self.landmask_is_watermask,
                            line_buffer=self.line_buffer,
                            verbose=self.verbose,
                        )
                    pbar.update()

        ## If no output filename was provided, we might want to return geometry objects directly
        return out_fn, cst_geoms

    
    def process_water(self, out_fn: str):
        """Process Water bodies (Lakes/Rivers)."""
        
        if self.osm_fetcher is None:
            self.osm_fetcher = self.fetch()
            
        lk_geoms = []
        if self.osm_fetcher.results:
             with utils.ccp(total=len(self.osm_fetcher.results), desc=f'Processing {self.q}', leave=self.verbose) as pbar:
                for entry in self.osm_fetcher.results:
                    if entry[-1] == 0:
                        local_file = entry[1]
                        try:
                            ds = ogr.Open(local_file)
                            if ds:
                                ## OSM driver usually produces 'multipolygons' layer for area features
                                layer = ds.GetLayer('multipolygons')
                                if layer:
                                    if self.region:
                                        layer.SetSpatialFilter(self.region.export_as_geom())
                                    
                                    for feat in layer:
                                        ## Filter by natural=water if specific tagging wasn't pre-filtered enough
                                        ## (Though Overpass query usually handles this)
                                        geom = feat.GetGeometryRef()
                                        if geom and not geom.IsEmpty():
                                            lk_geoms.append(geom.ExportToWkt())
                        except Exception as e:
                            utils.echo_warning_msg(f"Error reading {local_file}: {e}")
                    pbar.update()

        if out_fn and lk_geoms:
            gdalfun.ogr_wktgeoms2ogr(lk_geoms, out_fn, ogr_format='GPKG')
            
        return out_fn, lk_geoms

    
    def __call__(self, out_fn: Optional[str] = None, return_geom: bool = True, overwrite: bool = False):
        """Main execution entry point."""
        
        if self.region is None or not self.region.valid_p():
            utils.echo_error_msg(f'{self.region} is an invalid region')
            return None

        ## Determine Output Filename
        if out_fn is None or not isinstance(out_fn, str):
            out_bn = utils.append_fn(f'osm_{self.q}', self.region, 1, high_res=True)
            out_fn = os.path.join(self.cache_dir, f'{out_bn}.gpkg')

        utils.echo_msg_bold(f"Output: {out_fn}")
        
        ## Check Overwrite
        if not overwrite and os.path.exists(out_fn):
            if not return_geom:
                return out_fn
            ## If returning geom, we might need to load it from existing file? 
            ## Skipping for brevity.

        if overwrite:
            utils.remove_glob(out_fn)

        ## Dispatch Processing
        if self.q == 'coastline':
            out_fn, geoms = self.process(out_fn)
        elif self.q in ['lakes', 'water', 'rivers']:
            out_fn, geoms = self.process_water(out_fn)
        else:
            utils.echo_warning_msg(f"Unknown query type '{self.q}', processing as water.")
            out_fn, geoms = self.process_water(out_fn)
            
        if return_geom:            
            return geoms
        else:
            return out_fn


## ==============================================
## Command-line Interface (CLI)
## ==============================================
def osm_coast_cli():
    parser = argparse.ArgumentParser(description="Fetch and Process OSM Data")
    
    parser.add_argument("-R", "--region", help="Region [xmin/xmax/ymin/ymax]", required=True)
    parser.add_argument("-q", "--query", help="Query type (coastline, water, lakes, rivers)", default="coastline")
    parser.add_argument("-O", "--output", help="Output filename")
    parser.add_argument("--chunks", action="store_true", help="Chunk request into smaller tiles")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing output")
    parser.add_argument("--landmask", action="store_true", help="Include landmask in output (for coastline)")
    parser.add_argument("--invert", action="store_true", help="Invert land/water mask")

    args = parser.parse_args()
    
    ## Parse Region
    try:
        region = regions.Region().from_string(args.region)
    except Exception:
        print("Error: Invalid region string.")
        return

    ## Initialize Processor
    osm_proc = osmCoastline(
        region=region,
        q=args.query,
        chunks=args.chunks,
        include_landmask=args.landmask,
        landmask_is_watermask=args.invert,
    )
    
    ## Run
    osm_proc(out_fn=args.output, return_geom=False, overwrite=args.overwrite)

    
if __name__ == "__main__":
    osm_coast_cli()

    
### End
