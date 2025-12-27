### hrdem.py
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
##
## hrdem.py is part of CUDEM
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
## Fetch High-Resolution Digital Elevation Model (HRDEM) data from Canada (NRCAN).
##
### Code:

import os
from osgeo import ogr
from cudem import utils
from cudem.fetches import fetches

## ==============================================
## Constants
## ==============================================
HRDEM_FOOTPRINTS_URL = (
    'ftp://ftp.maps.canada.ca/pub/elevation/dem_mne/'
    'highresolution_hauteresolution/Datasets_Footprints.zip'
)
HRDEM_INFO_URL = (
    'https://open.canada.ca/data/en/dataset/'
    '957782bf-847c-4644-a757-e383c0057995#wb-auto-6'
)

## ==============================================
## HRDEM Module
## ==============================================
class HRDEM(fetches.FetchModule):
    """High-Resolution Digital Elevation Model data for Canada

    Fetch HRDEM data from Canada (NRCAN) via FTP footprints.
    
    https://open.canada.ca

    Configuration Example:
    < hrdem >
    """
    
    def __init__(self, **kwargs):
        super().__init__(name='hrdem', **kwargs)

    def run(self):
        """Run the HRDEM fetches module."""
        
        if self.region is None:
            return []

        ## Download Footprints Zip
        v_zip = os.path.join(self._outdir, 'Datasets_Footprints.zip')
        
        try:
            status = fetches.Fetch(
                HRDEM_FOOTPRINTS_URL,
                verbose=self.verbose
            ).fetch_ftp_file(v_zip)
            
            if status != 0:
                utils.echo_error_msg("Failed to download HRDEM footprints.")
                return self

            ## Unzip Shapefile components
            v_shps = utils.p_unzip(
                v_zip,
                ['shp', 'shx', 'dbf', 'prj'],
                outdir=self._outdir,
                verbose=self.verbose
            )
            
            ## Locate the .shp file
            v_shp = next((v for v in v_shps if v.endswith('.shp')), None)
            
            if not v_shp:
                utils.echo_error_msg("Could not find shapefile in HRDEM zip.")
                return self

            ## Open and Filter with OGR
            try:
                v_ds = ogr.Open(v_shp)
                if v_ds:
                    layer = v_ds.GetLayer()
                    
                    ## Spatial Filter
                    bbox_geom = self.region.export_as_geom()
                    layer.SetSpatialFilter(bbox_geom)
                    
                    for feature in layer:
                        data_link = feature.GetField('Ftp_dtm')
                        if data_link:
                            self.add_entry_to_results(
                                data_link,
                                data_link.split('/')[-1],
                                'raster'
                            )
                    
                    v_ds = None # Close Datasource
            except Exception as e:
                utils.echo_error_msg(f"Error reading HRDEM shapefile: {e}")

        except Exception as e:
            utils.echo_error_msg(f"Error running HRDEM fetch: {e}")
            
        finally:
            ## Cleanup
            if os.path.exists(v_zip):
                utils.remove_glob(v_zip)
            
            ## Cleanup extracted shapefiles if they exist
            if 'v_shps' in locals() and v_shps:
                utils.remove_glob(*v_shps)

        return self

### End
