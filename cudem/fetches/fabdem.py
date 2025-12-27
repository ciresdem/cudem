### fabdem.py
##
## Copyright (c) 2022 - 2025 Regents of the University of Colorado
##
## fabdem.py is part of CUDEM
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
## Fetch FABDEM elevation data.
##
### Code:

import os
from osgeo import ogr
from cudem import utils
from cudem.fetches import fetches

## ==============================================
## Constants
## ==============================================
FABDEM_FOOTPRINTS_URL = 'https://data.bris.ac.uk/datasets/s5hqmjcdj8yo2ibzi9b4ew3sn/FABDEM_v1-2_tiles.geojson'
FABDEM_INFO_URL = 'https://data.bris.ac.uk/data/dataset/s5hqmjcdj8yo2ibzi9b4ew3sn'
FABDEM_DATA_URL = 'https://data.bris.ac.uk/datasets/s5hqmjcdj8yo2ibzi9b4ew3sn'

## ==============================================
## FABDEM Module
## ==============================================
class FABDEM(fetches.FetchModule):
    """FABDEM elevation data
    
    FABDEM (Forest And Buildings removed Copernicus DEM) is a global elevation 
    map that removes building and tree height biases from the Copernicus GLO 30 
    Digital Elevation Model (DEM). The data is available at 1 arc second grid 
    spacing (approximately 30m at the equator).
    
    https://data.bris.ac.uk/data/dataset/s5hqmjcdj8yo2ibzi9b4ew3sn

    Configuration Example:    
    < fabdem >
    """
    
    def __init__(self, **kwargs):
        super().__init__(name='fabdem', **kwargs)

        self.data_format = -2 # zipfile
        self.src_srs = 'epsg:4326+3855' # WGS84 + EGM2008

        self.headers = {
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:89.0) Gecko/20100101 Firefox/89.0',
            'referer': FABDEM_INFO_URL
        }

        
    def run(self):
        """Run the FABDEM fetches module.
        
        Downloads the geojson footprint index, queries it spatially, 
        and adds intersecting tiles to the results.
        """
        
        if self.region is None:
            return []

        ## Define local path for the footprint file
        v_json = os.path.join(self._outdir, os.path.basename(FABDEM_FOOTPRINTS_URL))
        
        try:
            ## Fetch the footprint GeoJSON
            status = fetches.Fetch(
                FABDEM_FOOTPRINTS_URL, 
                verbose=self.verbose,
                headers=self.headers
            ).fetch_file(v_json)
            
            if status != 0:
                utils.echo_error_msg("Failed to download FABDEM footprints.")
                return self

            ## Open with OGR
            v_ds = ogr.Open(v_json)
            if v_ds is None:
                utils.echo_error_msg("Could not open FABDEM footprints file.")
                return self

            layer = v_ds.GetLayer()
            
            ## Spatial Filter
            bbox_geom = self.region.export_as_geom()
            layer.SetSpatialFilter(bbox_geom)
            
            for feature in layer:
                zipfile_name = feature.GetField('zipfile_name')
                if not zipfile_name:
                    continue
                    
                zipfile_url = f"{FABDEM_DATA_URL}/{zipfile_name}"
                
                ## Check for duplicates in current results
                if not any(entry['url'] == zipfile_url for entry in self.results):
                    self.add_entry_to_results(zipfile_url, zipfile_name, 'raster')
            
            v_ds = None # Close Datasource

        except Exception as e:
            utils.echo_error_msg(f"Error processing FABDEM: {e}")
            
        finally:
            ## Cleanup the index file
            if os.path.exists(v_json):
                utils.remove_glob(v_json)

        return self

### End
