### mbs.py
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
##
## mbs.py is part of CUDEM
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
### Code:

import os

from cudem import utils
from cudem import gdalfun
from cudem.waffles.waffles import Waffle

class WafflesMBGrid(Waffle):
    """SPLINE DEM via MB-System's mbgrid.
    
    Generate a DEM using MB-System's mbgrid command.
    By default, will use MB-Systems datalist processes.
    set `use_datalist=True` to use CUDEM's dlim instead.
    see mbgrid --help for more info

    -----------
    Parameters:
    
    dist=[val] - the dist variable to use in mbgrid
    tension=[val] - the spline tension value (0-inf)
    use_stack=[True/False] - use built-in datalists rather than 
                             mbdatalist
    
    < mbgrid:dist='10/3':tension=35:use_datalists=False >
    """
    
    def __init__(
            self,
            dist='10/3',
            tension=0,
            use_stack=True,
            nc=False,
            **kwargs
    ):
        super().__init__(**kwargs)
        self.nc = nc
        self.dist = dist
        self.tension = tension
        self.use_stack = use_stack

        
    def _gmt_num_msk(self, num_grd, dst_msk):
        """generate a num-msk from a NUM grid using GMT grdmath

        Args:
          num_grd (str): pathname to a source `num` grid file
          dst_msk (str): pathname to a destination `msk` grid file

        Returns:
          list: [cmd-output, cmd-return-code]
        """

        if self.gc['MBGRID'] is None:
            utils.echo_error_msg(
                'MB-System must be installed to use the MBGRID module'
            )
            return(None, -1)

        if self.gc['GMT'] is None:
            utils.echo_error_msg(
                'GMT must be installed to use the MBGRID module'
            )
            return(None, -1)

        num_msk_cmd = 'gmt grdmath -V "{}" 0 MUL 1 ADD 0 AND = "{}"'.format(
            num_grd, dst_msk
        )
        return(
            utils.run_cmd(
                num_msk_cmd, verbose=self.verbose
            )
        )

    
    def _gmt_grd2gdal(self, src_grd, dst_fmt='GTiff'):
        """convert the grd file to tif using GMT

        Args:
          src_grd (str): a pathname to a grid file
          dst_fmt (str): the output GDAL format string

        Returns:
          str: the gdal file name or None
        """

        dst_gdal = '{}.{}'.format(
            os.path.basename(src_grd).split('.')[0], gdalfun.gdal_fext(dst_fmt)
        )        
        grd2gdal_cmd = 'gmt grdconvert "{}" "{}"=gd+n{}:{} -V'.format(
            src_grd, dst_gdal, self.ndv, dst_fmt
        )
        out, status = utils.run_cmd(
            grd2gdal_cmd, verbose=self.verbose
        )
        if status == 0:
            return(dst_gdal)
        else:
            return(None)

        
    def _gmt_grdsample(self, src_grd, dst_fmt='GTiff'):
        """convert the grd file to tif using GMT

        Args:
          src_grd (str): a pathname to a grid file
          dst_fmt (str): the output GDAL format string

        Returns:
          str: the gdal file name or None
        """

        dst_gdal = '{}.{}'.format(
            os.path.basename(src_grd).split('.')[0], galfun.gdal_fext(dst_fmt)
        )
        grdsample_cmd = 'gmt grdsample "{}" -T -G{}=gd+n{}:{} -V'.format(
            src_grd, dst_gdal, self.ndv, dst_fmt
        )        
        out, status = utils.run_cmd(
            grdsample_cmd, verbose=self.verbose
        )        
        if status == 0:
            return(dst_gdal)        
        else:
            return(None)

        
    def stack2mbdatalist(self):
        mb_datalist_fn = os.path.join(self.cache_dir, '_tmp_mb.datalist')
        mb_stack_xyz = os.path.join(
            self.cache_dir, '{}.xyz'.format(utils.fn_basename2(self.stack))
        )
        with open(mb_datalist_fn, 'w') as dl:
            dl.write('{} 168 1'.format(mb_stack_xyz))
                
        with open(mb_stack_xyz, 'w') as stack_xyz:
            self.dump_xyz(stack_xyz)

        return(mb_datalist_fn)

    
    def run(self):
        mb_datalist = self.stack2mbdatalist() if self.use_stack else self.data.fn
        #self.mb_region = self.ps_region.copy()
        self.mb_region = self.ps_region.copy()
        #out_name = os.path.join(self.cache_dir, self.name)
        out_name = self.name
        
        mb_xcount, mb_ycount, mb_gt \
            = self.ps_region.geo_transform(
                x_inc=self.xinc, y_inc=self.yinc, node=self.node
            )
        mbgrid_cmd = 'mbgrid -I{} {} -D{}/{} -O{} -A2 -F1 -N -C{} -S0 -X0 -T{}'.format(
            mb_datalist,
            self.mb_region.format('gmt'),
            mb_xcount,
            mb_ycount,
            out_name,
            self.dist,
            self.tension
        )
        out, status = utils.run_cmd(mbgrid_cmd, verbose=self.verbose)
        if status == 0:
            out = gdalfun.gdal2gdal(
                '{}.grd'.format(out_name),
                dst_dem=self.fn,
                dst_fmt=self.fmt,
                co=True
            )
                    
        return(self)


### End
