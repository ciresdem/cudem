### bluetopo.py
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
##
## bluetopo.py is part of CUDEM
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
## Fetch data from NOAA's BlueTopo Bathymetric Source.
##
### Code:

import os
import boto3
from botocore import UNSIGNED
from botocore.client import Config
from osgeo import ogr
from cudem import utils
from cudem.fetches import fetches

## ==============================================
## Constants
## ==============================================
BLUETOPO_BUCKET = 'noaa-ocs-nationalbathymetry-pds'
BLUETOPO_PREFIX = 'BlueTopo'

## ==============================================
## BlueTopo Module
## ==============================================
class BlueTopo(fetches.FetchModule):
    """BlueTOPO DEM
    
    BlueTopo is a compilation of the nation's best available bathymetric data. 
    Created as part of the Office of Coast Survey nautical charting 
    mission and its National Bathymetric Source project.

    Output 'tiff' files are 3 bands:
    1 - Elevation (NAVD88)
    2 - Uncertainty
    3 - Data Source Table

    https://nauticalcharts.noaa.gov/data/bluetopo_specs.html
    https://noaa-ocs-nationalbathymetry-pds.s3.amazonaws.com/index.html

    Configuration Example:    
    < bluetopo:want_interpolation=False:unc_weights=False:keep_index=False >
    """
    
    def __init__(
            self,
            want_interpolation: bool = False,
            unc_weights: bool = False,
            keep_index: bool = False,
            **kwargs
    ):
        super().__init__(name='bluetopo', **kwargs)
        self.unc_weights = unc_weights
        self.want_interpolation = want_interpolation
        self.keep_index = keep_index
        
        ## DLIM variables
        self.data_format = 200
        
        ## State variables for index file
        self._bluetopo_index_url = None
        self._bluetopo_index_fn = None

        
    def _get_s3_client(self):
        """Return an anonymous S3 client."""
        
        return boto3.client('s3', config=Config(signature_version=UNSIGNED))

    
    def _get_index_url(self, s3_client) -> str:
        """Dynamically find the Tile Scheme index file URL from S3."""
        
        try:
            r = s3_client.list_objects(
                Bucket=BLUETOPO_BUCKET, 
                Prefix=f'{BLUETOPO_PREFIX}/_BlueTopo_Tile_Scheme'
            )
            
            if 'Contents' in r and len(r['Contents']) > 0:
                key = r['Contents'][0]['Key']
                return f'https://{BLUETOPO_BUCKET}.s3.amazonaws.com/{key}'
        except Exception as e:
            if self.verbose:
                utils.echo_error_msg(f"Error listing BlueTopo index: {e}")
        
        return None

    
    def run(self):
        """Run the BlueTopo fetch module."""
        
        s3 = self._get_s3_client()
        
        ## Locate and Fetch the Index File
        if self._bluetopo_index_url is None:
            self._bluetopo_index_url = self._get_index_url(s3)
            
        if not self._bluetopo_index_url:
            utils.echo_error_msg("Could not locate BlueTopo tile index.")
            return self

        self._bluetopo_index_fn = os.path.basename(self._bluetopo_index_url)
        
        try:
            ## Download index
            status = fetches.Fetch(
                self._bluetopo_index_url, 
                verbose=self.verbose
            ).fetch_file(self._bluetopo_index_fn)
            
            if status != 0:
                raise IOError("Failed to download BlueTopo index.")

            ## Open Index and Filter
            v_ds = ogr.Open(self._bluetopo_index_fn)
            if v_ds is None:
                raise IOError("Failed to open BlueTopo index OGR source.")

            layer = v_ds.GetLayer()
            if self.region is not None:
                _boundsGeom = self.region.export_as_geom()
                layer.SetSpatialFilter(_boundsGeom)            
            
            ## Iterate Features and Find Data
            for feature in layer:
                if feature is None:
                    continue
                
                tile_name = feature.GetField('tile')
                
                ## List objects for this specific tile to find the TIFF
                try:
                    r = s3.list_objects(
                        Bucket=BLUETOPO_BUCKET,
                        Prefix=f'{BLUETOPO_PREFIX}/{tile_name}'
                    )
                    
                    if 'Contents' in r:
                        for obj in r['Contents']:
                            key = obj['Key']
                            if key.endswith('.tiff'):
                                data_link = f'https://{BLUETOPO_BUCKET}.s3.amazonaws.com/{key}'
                                self.add_entry_to_results(
                                    data_link,
                                    os.path.basename(key),
                                    'raster'
                                )
                except Exception as e:
                    if self.verbose:
                        utils.echo_error_msg(f"Error searching tile {tile_name}: {e}")

            v_ds = None

        except Exception as e:
            utils.echo_error_msg(f"BlueTopo Run Error: {e}")
            
        finally:
            ## Cleanup
            if not self.keep_index and self._bluetopo_index_fn:
                utils.remove_glob(self._bluetopo_index_fn)
            
        return self

### End
