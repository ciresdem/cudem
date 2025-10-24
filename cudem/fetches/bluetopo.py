### bluetopo.py
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
##
## fetches.py is part of CUDEM
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
###############################################################################
### Commentary:
##
##
### Code:

import boto3 # boto3 for aws api
from osgeo import ogr
from cudem import utils
from cudem.fetches import fetches

## BlueTopo
class BlueTopo(fetches.FetchModule):
    """BlueTOPO DEM
    
    BlueTopo is a compilation of the nation's best available bathymetric data. 
    In the same way that topographic map details the height of land, BlueTopo 
    details the depth of lake beds and seafloor beneath navigationally significant 
    U.S. waters. Created as part of the Office of Coast Survey nautical charting 
    mission and its National Bathymetric Source project, BlueTopo is curated 
    bathymetric source data to provide a definitive nationwide model of the seafloor 
    and the Great Lakes.

    Output 'tiff' files are 3 bands
    1 - Elevation
    2 - Uncertainty
    3 - Data Source Table

    yield_xyz outputs elevation (band 1)
    elevation data is in NAVD88

    https://nauticalcharts.noaa.gov/data/bluetopo_specs.html
    https://noaa-ocs-nationalbathymetry-pds.s3.amazonaws.com/index.html#
    https://www.nauticalcharts.noaa.gov/data/bluetopo.html
    https://noaa-ocs-nationalbathymetry-pds.s3.amazonaws.com/index.html#BlueTopo/
    
    https://www.nauticalcharts.noaa.gov/data/bluetopo.html
    
    < bluetopo:want_interpolation=False:unc_weights=False:keep_index=False >
    """
    
    def __init__(
            self,
            want_interpolation=False,
            unc_weights=False,
            keep_index=False,
            **kwargs
    ):
        super().__init__(name='bluetopo', **kwargs)
        self.unc_weights = unc_weights
        self.want_interpolation = want_interpolation
        self.keep_index = keep_index
        
        ## BlueTopo uses AWS
        self._bt_bucket = 'noaa-ocs-nationalbathymetry-pds'
        s3 = boto3.client('s3', aws_access_key_id='', aws_secret_access_key='')
        s3._request_signer.sign = (lambda *args, **kwargs: None)
        r = s3.list_objects(Bucket = self._bt_bucket, Prefix='BlueTopo/_BlueTopo_Tile_Scheme')
        self._bluetopo_index_url = 'https://{}.s3.amazonaws.com/{}'.format(
            self._bt_bucket, r['Contents'][0]['Key']
        )
        self._bluetopo_index = self._bluetopo_index_url.split('/')[-1]

        
    def run(self):
        s3 = boto3.client('s3', aws_access_key_id='', aws_secret_access_key='')
        s3._request_signer.sign = (lambda *args, **kwargs: None)
        try:
            status = fetches.Fetch(
                self._bluetopo_index_url, verbose=self.verbose
            ).fetch_file(self._bluetopo_index)
            v_ds = ogr.Open(self._bluetopo_index)
        except:
            v_ds = None
            status = -1
            
        if v_ds is not None:
            layer = v_ds.GetLayer()
            _boundsGeom = self.region.export_as_geom()
            layer.SetSpatialFilter(_boundsGeom)            
            fcount = layer.GetFeatureCount()
            for feature in layer:
                if feature is None:
                    continue
                
                tile_name = feature.GetField('tile')
                r = s3.list_objects(
                    Bucket='noaa-ocs-nationalbathymetry-pds',
                    Prefix='BlueTopo/{}'.format(tile_name)
                )
                if 'Contents' in r:
                    for key in r['Contents']:
                        if key['Key'].split('.')[-1] == 'tiff':
                            data_link = ('https://noaa-ocs-nationalbathymetry-pds'
                                         f'.s3.amazonaws.com/{key["Key"]}')
                            self.add_entry_to_results(
                                data_link,
                                data_link.split('/')[-1],
                                'raster'
                            )
            v_ds = None

        if not self.keep_index:
            utils.remove_glob(self._bluetopo_index)
            
        return(self)

### End
