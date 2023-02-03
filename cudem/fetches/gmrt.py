### gmrt.py - GMRT dataset
##
## Copyright (c) 2010 - 2023 Regents of the University of Colorado
##
## gmrt.py is part of CUDEM
##
## Permission is hereby granted, free of charge, to any person obtaining a copy 
## of this software and associated documentation files (the "Software"), to deal 
## in the Software without restriction, including without limitation the rights 
## to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies 
## of the Software, and to permit persons to whom the Software is furnished to do so, 
## subject to the following conditions:
##
## The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
##
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
## INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR 
## PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE 
## FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
## ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
##
### Commentary:
##
## GMRT Fetch
##
## fetch extracts of the GMRT. - Global Extents
## https://www.gmrt.org/index.php
## https://www.gmrt.org/services/gridserverinfo.php#!/services/getGMRTGridURLs
##
## The Global Multi-Resolution Topography (GMRT) synthesis is a multi-resolutional 
## compilation of edited multibeam sonar data collected by scientists and institutions worldwide, that is 
## reviewed, processed and gridded by the GMRT Team and merged into a single continuously updated compilation 
## of global elevation data. The synthesis began in 1992 as the Ridge Multibeam Synthesis (RMBS), was expanded 
## to include multibeam bathymetry data from the Southern Ocean, and now includes bathymetry from throughout 
## the global and coastal oceans.
##
## The GridServer service provides access to gridded data from the Global Multi-resolution Topography (GMRT) Synthesis.
## Requested data may be up to 2GB, or approximately 20 by 20 degrees at 100 meters per node (maximum available resolution).
## A variety of output formats are supported.
##
##    Data Formats
##        GMT v3 Compatible NetCDF (GMT id=cf)
##        COARDS/CF1.6 Compliant NetCDF (GMT id=nd)
##        ESRI ArcASCII
##        GeoTIFF
##    Metadata Formats
##        XML (metadata)
##        JSON (metadata)
##        Plain text (metadata)
##
## Layers:
##  topo - gmrt with gebco behind mb
##  topo-mask - gmrt mb surveys with gebco masked
##
## Data is assumed instantaneous MSL (5773?)
##
### Code:

import os

from osgeo import gdal

from cudem import utils
from cudem import regions
from cudem import datasets
from cudem import demfun

import cudem.fetches.utils as f_utils

class GMRT(f_utils.FetchModule):
    '''Fetch raster data from the GMRT'''
    
    def __init__(self, res='default', fmt='geotiff', bathy_only=False, upper_limit=None, lower_limit=None, layer='topo', **kwargs):
        super().__init__(name='gmrt', **kwargs) 

        self._gmrt_grid_url = "https://www.gmrt.org:443/services/GridServer?"
        self._gmrt_grid_urls_url = "https://www.gmrt.org:443/services/GridServer/urls?"
        self._gmrt_grid_metadata_url = "https://www.gmrt.org/services/GridServer/metadata?"        
        self.res = res
        self.fmt = fmt
        if layer != 'topo' and layer != 'topo-mask':
            self.layer = 'topo'
        else:
            self.layer = layer

        self.upper_limit = upper_limit
        self.lower_limit = lower_limit
        if bathy_only:
            self.upper_limit = 0
            
        self.gmrt_region = self.region.copy()
        self.gmrt_region.buffer(pct=2.33,x_inc=.0088,y_inc=.0088)
        self.gmrt_region._wgs_extremes(just_below=True)

        self.gmrt_p_region = self.region.copy()
        self.gmrt_p_region.zmax = utils.float_or(upper_limit)
        self.gmrt_p_region.zmin = utils.float_or(lower_limit)
        
    def run(self):
        '''Run the GMRT fetching module'''

        if self.region is None:
            return([])

        self.data = {
            'north':self.gmrt_region.ymax,
            'west':self.gmrt_region.xmin,
            'south':self.gmrt_region.ymin,
            'east':self.gmrt_region.xmax,
            'mformat':'json',
            'resolution':self.res,
            'format':self.fmt,
            'layer':self.layer,
        }

        req = f_utils.Fetch(self._gmrt_grid_url).fetch_req(
            params=self.data, tries=10, timeout=2
        )
        if req is not None:
            outf = 'gmrt_{}_{}_{}.{}'.format(self.layer, self.res, self.region.format('fn'), 'tif' if self.fmt == 'geotiff' else 'grd')
            self.results.append([req.url, os.path.join(self._outdir, outf), 'gmrt']) 
            
            # try:
            #     gmrt_urls = req.json()
            # except Exception as e:
            #     utils.echo_error_msg(e)
            #     gmrt_urls = []

            # for url in gmrt_urls:
            #     if self.layer == 'topo-mask':
            #         url = url.replace('topo', 'topo-mask')

            #     opts = {}
            #     for url_opt in url.split('?')[1].split('&'):
            #         opt_kp = url_opt.split('=')
            #         opts[opt_kp[0]] = opt_kp[1]
                    
            #     url_region = regions.Region().from_list([
            #         float(opts['west']),
            #         float(opts['east']),
            #         float(opts['south']),
            #         float(opts['north'])
            #     ])
            #     outf = 'gmrt_{}_{}.{}'.format(opts['layer'], url_region.format('fn'), 'tif' if self.fmt == 'geotiff' else 'grd')
            #     self.results.append([url, os.path.join(self._outdir, outf), 'gmrt'])
                
        return(self)

    def yield_ds(self, entry):
        src_data = 'gmrt_tmp.{}'.format('tif' if self.fmt == 'geotiff' else 'nc')
        if f_utils.Fetch(
                entry[0], callback=self.callback, verbose=self.verbose
        ).fetch_file(src_data,timeout=10, read_timeout=120) == 0:
            gmrt_ds = datasets.RasterFile(
                fn=src_data,
                data_format=200,
                #src_srs='epsg:4326+3855',
                src_srs='epsg:4326',
                dst_srs=self.dst_srs,
                x_inc=self.x_inc,
                y_inc=self.y_inc,
                weight=self.weight,
                src_region=self.gmrt_p_region,
                verbose=self.verbose
            )

            yield(gmrt_ds)

        else:
            utils.echo_error_msg('failed to fetch remote file, {}...'.format(src_data))
            
        utils.remove_glob('{}*'.format(src_data))
            
    def yield_xyz(self, entry):
        for ds in self.yield_ds(entry):
            for xyz in ds.yield_xyz():
                yield(xyz)

    def yield_array(self, entry):
        for ds in self.yield_ds(entry):
            for arr in ds.yield_array():
                yield(arr)
        
### End
