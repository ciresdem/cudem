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
### Commentary:
##
## Wraps the MB-System 'mbgrid' command for interpolation.
##
### Code:

import os

from cudem import utils
from cudem import gdalfun
from cudem.waffles.waffles import Waffle

class WafflesMBGrid(Waffle):
    """SPLINE DEM via MB-System's mbgrid.
    
    Generate a DEM using MB-System's mbgrid command.
    By default, will use built-in datalists. Set `use_stack=False` to use 
    an external MB-System datalist file.
    
    Parameters:
    -----------    
    dist (str) : the dist variable to use in mbgrid (default: '10/3')
    tension (float) : the spline tension value 0-inf (default: 0)
    use_stack (bool) : use built-in datalists rather than mbdatalist (default: True)
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
        self.gc = utils.config_check(chk_config_file=False)

        
    def _gmt_num_msk(self, num_grd, dst_msk):
        """Generate a num-msk from a NUM grid using GMT grdmath."""

        if self.gc['MBGRID'] is None:
            utils.echo_error_msg('MB-System must be installed to use the MBGRID module')
            return None, -1

        if self.gc['GMT'] is None:
            utils.echo_error_msg('GMT must be installed to use the MBGRID module')
            return None, -1

        num_msk_cmd = f'gmt grdmath -V "{num_grd}" 0 MUL 1 ADD 0 AND = "{dst_msk}"'
        return utils.run_cmd(num_msk_cmd, verbose=self.verbose)

    
    def _gmt_grd2gdal(self, src_grd, dst_fmt='GTiff'):
        """Convert the grd file to tif using GMT."""

        dst_gdal = f'{os.path.basename(src_grd).split('.')[0]}.{gdalfun.gdal_fext(dst_fmt)}'
        grd2gdal_cmd = f'gmt grdconvert "{src_grd}" "{dst_gdal}"=gd+n{self.ndv}:{dst_fmt} -V'
        
        out, status = utils.run_cmd(grd2gdal_cmd, verbose=self.verbose)
        if status == 0:
            return dst_gdal
        else:
            return None

        
    def _gmt_grdsample(self, src_grd, dst_fmt='GTiff'):
        """Resample the grd file to tif using GMT."""

        dst_gdal = f'{os.path.basename(src_grd).split('.')[0]}.{gdalfun.gdal_fext(dst_fmt)}'
        grdsample_cmd = f'gmt grdsample "{src_grd}" -T -G{dst_gdal}=gd+n{self.ndv}:{dst_fmt} -V'
        
        out, status = utils.run_cmd(grdsample_cmd, verbose=self.verbose)        
        if status == 0:
            return dst_gdal        
        else:
            return None

        
    def stack2mbdatalist(self):
        """Convert internal processing stack to MB-System Datalist."""
        
        mb_datalist_fn = os.path.join(self.cache_dir, '_tmp_mb.datalist')
        mb_stack_xyz = os.path.join(self.cache_dir, f'{utils.fn_basename2(self.stack)}.xyz')
        
        with open(mb_datalist_fn, 'w') as dl:
            # Format 168 is ASCII XYZ
            dl.write(f'{mb_stack_xyz} 168 1\n')
                
        with open(mb_stack_xyz, 'w') as stack_xyz:
            self.dump_xyz(stack_xyz)

        return mb_datalist_fn

    
    def run(self):
        """Execute the MBGRID process."""
        
        if self.gc['MBGRID'] is None:
            utils.echo_error_msg('MB-System must be installed to use the MBGRID module')
            return self

        mb_datalist = self.stack2mbdatalist() if self.use_stack else self.data.fn
        self.mb_region = self.ps_region.copy()
        out_name = self.name
        
        mb_xcount, mb_ycount, _ = self.ps_region.geo_transform(
            x_inc=self.xinc, y_inc=self.yinc, node=self.node
        )
        
        ## mbgrid arguments:
        # -I : Input Datalist
        # -R : Region (West/East/South/North)
        # -D : Dimensions (X/Y)
        # -O : Output Basename
        # -A2 : Thinning Algorithm (2=Mean)
        # -F1 : Interpolation Algorithm (1=Spline)
        # -N  : NoData value handling
        # -C  : Clipping Distance (dist)
        # -S0 : Slope calculation (0=Off)
        # -X0 : Extended output (0=Off)
        # -T  : Tension
        
        mbgrid_cmd = (
            f'mbgrid -I{mb_datalist} {self.mb_region.format("gmt")} '
            f'-D{mb_xcount}/{mb_ycount} -O{out_name} '
            f'-A2 -F1 -N -C{self.dist} -S0 -X0 -T{self.tension}'
        )
        
        out, status = utils.run_cmd(mbgrid_cmd, verbose=self.verbose)
        
        if status == 0:
            ## Convert resulting GRD to desired GDAL format
            gdalfun.gdal2gdal(
                f'{out_name}.grd',
                dst_dem=self.fn,
                dst_fmt=self.fmt,
                co=True
            )
            ## Cleanup GRD
            utils.remove_glob(f'{out_name}.grd')
            
        return self

### End
