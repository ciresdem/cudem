### gmtfilters.py - GRId filTerS
##
## Copyright (c) 2024 - 2025 Regents of the University of Colorado
##
## Permission is hereby granted, free of charge, to any person obtaining a copy 
## of this software and associated documentation files (the "Software"), to deal 
## in the Software without restriction, including without limitation the rights 
## to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies 
## of the Software, and to permit persons to whom the Software is furnished to do so, 
## subject to the following conditions:
##
## The above copyright notice and this permission notice shall be included in all copies or
## substantial portions of the Software.
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

from cudem import utils
#from cudem import gdalfun
from cudem.grits import grits

class GMTgrdfilter(grits.Grits):
    """Filter a DEM through GMT's `grdfilter`; see `gmt grdfilter --help`

    -----------
    Parameters:

    filter_type(str) - The grdfilter filter type (grdfilter -F)
    dist(str) - The grdfilter distance value (grdfilter -D)
    node(str) - Either 'grid' or 'pixel'
    """
    
    def __init__(self, filter_type: str = 'c3s', dist: any = '1',
                 node: str = 'pixel', **kwargs: any):
        super().__init__(**kwargs)
        self.filter_type = filter_type
        self.dist = dist
        self.node = node

        
    def run(self):
        """filter `src_dem` using GMT grdfilter"""

        #if self.dst_dem is None:
        tmp_dst_dem = utils.make_temp_fn('{}_filtered.{}'.format(
            utils.fn_basename2(self.src_dem), utils.fn_ext(self.src_dem)
        ), self.cache_dir)
        #else:
        #    tmp_dst_dem = self.dst_dem

        ft_cmd1 = (
            'gmt grdfilter -V {} -G{}=gd:GTiff -F{} -D1{}'.format(
                self.src_dem,
                self.dst_dem,#tmp_dst_dem,
                self.filter_type,
                self.dist,
                ' -rp' if self.node == 'pixel' else ''
            )
        )

        out, status = utils.run_cmd(ft_cmd1, verbose=self.verbose)

        #out_array = gdalfun.gdal_get_array(tmp_dst_dem, 1)
        # with gdalfun.gdal_datasource(self.dst_dem, update=True) as dst_ds:
        #     if dst_ds is not None:
        #         dst_band = dst_ds.GetRasterBand(self.band)
        #         dst_band.Write(out_array)
        #     else:
        #         utils.echo_error_msg(f'coult not open output {self.dst_dem}')
        
        return(self.dst_dem, 0)

### End
