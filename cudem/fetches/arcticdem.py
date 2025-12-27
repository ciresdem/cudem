### arcticdem.py
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
##
## arcticdem.py is part of CUDEM
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
## Fetch data from the ArcticDEM project.
##
### Code:

import os
from typing import List, Optional
from osgeo import ogr
from cudem import utils
from cudem import regions
from cudem.fetches import fetches
from cudem.fetches import FRED

## ==============================================
## Constants
## ==============================================
ARCTIC_DEM_INDEX_URL = 'https://data.pgc.umn.edu/elev/dem/setsm/ArcticDEM/indexes/ArcticDEM_Tile_Index_Rel7.zip'

## ==============================================
## ArcticDEM Module
## ==============================================
class ArcticDEM(fetches.FetchModule):
    """Arctic DEM

    ArcticDEM is an NGA-NSF public-private initiative to automatically 
    produce a high-resolution, high quality, digital surface model (DSM) 
    of the Arctic using optical stereo imagery.

    https://www.pgc.umn.edu/data/arcticdem/

    Configuration Example:
    < arcticdem >
    """
    
    def __init__(self, where: str = '1=1', layer: int = 0, **kwargs):
        super().__init__(name='arcticdem', **kwargs)
        self.where = [where] if where else []
        
        ## Warp the input region to EPSG:3413 (Polar Stereographic North)
        ## for spatial filtering against the native index.
        if self.region:
            self.arctic_region = self.region.copy()
            self.arctic_region.warp('epsg:3413')
        else:
            self.arctic_region = None

        self.FRED = FRED.FRED(name=self.name, verbose=self.verbose)
        self.update_if_not_in_FRED()

        
    def update_if_not_in_FRED(self):
        """Check if ArcticDEM is in the local FRED database; if not, add it."""
        
        self.FRED._open_ds()
        try:
            self.FRED._attribute_filter([f"DataSource = '{self.name}'"])
            if len(self.FRED.layer) == 0:
                self.FRED._close_ds()
                self.update()
        finally:
             if self.FRED.ds is not None:
                self.FRED._close_ds()

            
    def _fetch_and_unzip_index(self, dst_dir: str) -> List[str]:
        """Helper to fetch and unzip the index shapefile."""
        
        v_zip = os.path.basename(ARCTIC_DEM_INDEX_URL)
        dst_zip = os.path.join(dst_dir, v_zip)
        
        status = fetches.Fetch(
            ARCTIC_DEM_INDEX_URL,
            verbose=self.verbose
        ).fetch_file(dst_zip)

        if status != 0:
            return []

        v_shps = utils.p_unzip(
            dst_zip, 
            ['shp', 'shx', 'dbf', 'prj'],
            outdir=dst_dir
        )
        
        ## Cleanup zip if downloaded locally
        if os.path.exists(dst_zip):
            utils.remove_glob(dst_zip)
            
        return v_shps

    
    def update(self):
        """Update the FRED reference vector with the ArcticDEM bounding box."""
        
        self.FRED._open_ds()
        
        try:        
            ## Fetch Index to current working directory or specific location
            v_shps = self._fetch_and_unzip_index(self._outdir)
            v_shp = next((v for v in v_shps if v.endswith('.shp')), None)
            
            if v_shp:
                ## Reproject to WGS84 for FRED storage
                temp_shp = os.path.join(self._outdir, 'arctic_tmp.shp')
                utils.run_cmd(
                    f'ogr2ogr {temp_shp} {v_shp} -t_srs epsg:4326',
                    verbose=self.verbose
                )
                
                temp_files = [
                    temp_shp, 
                    temp_shp.replace('.shp', '.dbf'), 
                    temp_shp.replace('.shp', '.shx'), 
                    temp_shp.replace('.shp', '.prj')
                ]

                ## Calculate merged bounding region from the reprojected shapefile
                shp_regions = regions.gdal_ogr_regions(temp_shp)
                combined_region = regions.Region()
                
                for this_region in shp_regions:
                    if combined_region.is_valid(check_xy=True):
                        combined_region = regions.regions_merge(combined_region, this_region)
                    else:
                        combined_region = this_region
                
                geom = combined_region.export_as_geom()

                ## Add to FRED
                self.FRED._attribute_filter(["ID = 'ARCTICDEM-1'"])
                if self.FRED.layer is None or len(self.FRED.layer) == 0:
                    self.FRED._add_survey(
                        Name='ArcticDEM',
                        ID='ARCTICDEM-1',
                        Agency='UMN',
                        Date=utils.this_year(),
                        MetadataLink='https://data.pgc.umn.edu/elev/dem/setsm/ArcticDEM/',
                        MetadataDate=str(utils.this_year()),
                        DataLink=ARCTIC_DEM_INDEX_URL,
                        IndexLink=ARCTIC_DEM_INDEX_URL,
                        DataType='raster',
                        DataSource='arcticdem',
                        Info='Arctic Only',
                        geom=geom
                    )

                ## Cleanup temp reprojected files
                utils.remove_glob(*temp_files)

        except Exception as e:
            utils.echo_error_msg(f"Failed to update ArcticDEM FRED entry: {e}")
        finally:                
            ## Cleanup downloaded index files
            utils.remove_glob(*v_shps)
            self.FRED._close_ds()

        
    def run(self):
        """Run the ArcticDEM fetches module."""
        
        if self.arctic_region is None:
            return self

        ## Get the Index Shapefile
        v_shps = self._fetch_and_unzip_index(self._outdir)
        v_shp = next((v for v in v_shps if v.endswith('.shp')), None)
        
        if not v_shp:
            return self

        v_ds = None
        try:
            v_ds = ogr.Open(v_shp)
            if v_ds is not None:
                layer = v_ds.GetLayer()
                
                ## Filter using the warped region (EPSG:3413)
                _bounds_geom = self.arctic_region.export_as_geom()
                layer.SetSpatialFilter(_bounds_geom)
                
                ## Apply additional attribute filters if specified
                for where_clause in self.where:
                    layer.SetAttributeFilter(where_clause)
                
                fcount = layer.GetFeatureCount()
                if self.verbose:
                    utils.echo_msg(f'Filtered {fcount} ArcticDEM features')
                
                ## Iterate results
                for feature in layer:
                    data_link = feature.GetField('fileurl')
                    if data_link:
                        self.add_entry_to_results(
                            data_link, 
                            os.path.basename(data_link), 
                            'raster'
                        )
        except Exception as e:
            utils.echo_error_msg(f"Error reading ArcticDEM index: {e}")
        finally:
            v_ds = None  # Close Datasource
            utils.remove_glob(*v_shps) 

        return self

### End
