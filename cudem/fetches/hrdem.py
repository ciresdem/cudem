### hrdem.py
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

import os
from osgeo import ogr
from cudem import utils
from cudem.fetches import fetches

## HRDEM - Canada Topo
class HRDEM(fetches.FetchModule):
    """High-Resolution Digital Elevation Model data for Canada

    Fetch HRDEM data from Canada (NRCAN)
    
    https://open.canada.ca

    < hrdem >
    """
    
    def __init__(self, **kwargs):
        super().__init__(name='hrdem', **kwargs)

        ## hrdem URLs
        self._hrdem_footprints_url = ('ftp://ftp.maps.canada.ca/pub/elevation/'
                                      'dem_mne/highresolution_hauteresolution/'
                                      'Datasets_Footprints.zip')
        self._hrdem_info_url = ('https://open.canada.ca/data/en/dataset/'
                                '957782bf-847c-4644-a757-e383c0057995#wb-auto-6')

        
    def run(self):
        """Run the HRDEM fetches module"""

        # use the remote footprints to discover data.
        v_zip = os.path.join(self._outdir, 'Datasets_Footprints.zip') 
        status = fetches.Fetch(
            self._hrdem_footprints_url,
            verbose=self.verbose
        ).fetch_ftp_file(v_zip)
        v_shps = utils.p_unzip(
            v_zip,
            ['shp', 'shx', 'dbf', 'prj'],
            verbose=self.verbose
        )
        v_shp = None
        for v in v_shps:
            if v.split('.')[-1] == 'shp':
                v_shp = v
                break
        try:
            v_ds = ogr.Open(v_shp)
        except:
            v_ds = None
            status = -1
                
        if v_ds is not None:
            layer = v_ds.GetLayer()
            fcount = layer.GetFeatureCount()
            for f in range(0, fcount):
                feature = layer[f]
                geom = feature.GetGeometryRef()
                if geom.Intersects(self.region.export_as_geom()):
                    data_link = feature.GetField('Ftp_dtm')
                    self.add_entry_to_results(
                        data_link,
                        data_link.split('/')[-1],
                        'raster'
                    )
                    
            v_ds = None

        utils.remove_glob(v_zip, *v_shps)

### End
