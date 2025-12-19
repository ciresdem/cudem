### gmtgrid.py
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
##
## gmtgrid.py is part of CUDEM
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

from cudem import utils
from cudem import gdalfun
from cudem.waffles.waffles import Waffle

class GMTSurface(Waffle):
    """SPLINE DEM via GMT surface
    
    Generate a DEM using GMT's surface command
    see gmt surface --help for more info.

    -----------
    Parameters:
   
    tension=[0-1] - spline tension.
    relaxation=[val] - spline relaxation factor.
    aspect=[val/None] - gridding aspect
    breakline=[path/None] - use xyz dataset at `path` as a breakline
    convergence=[val/None] - gridding convergence
    blockmean=[True/False] - pipe the data through gmt blockmean before gridding
    geographic=[True/Faslse] - data/grid are geographic
    pixel_node=[True/False] - grid in pixel-node

    < gmt-surface:tension=.35:relaxation=1:max_radius=None:aspect=None:breakline=None:convergence=None:blockmean=False:geographic=True >
    """
    
    def __init__(self, tension=1, relaxation=1, max_radius=None,
                 aspect=None, breakline=None, convergence=None,
                 blockmean=False, geographic=True, pixel_node=False,
                 **kwargs):
        super().__init__(**kwargs)
        if utils.float_or(tension) is not None:
            if utils.float_or(tension) > 1:
                self.tension = 1
            else:
                self.tension = tension
        else:
            self.tension = .35

        self.convergence = utils.float_or(convergence)
        self.relaxation = relaxation
        self.breakline = breakline
        self.max_radius = max_radius
        self.aspect = aspect
        self.blockmean = blockmean
        self.geographic = geographic
        self.pixel_node = pixel_node

        self.gc = utils.config_check(chk_config_file=False)
        
    def run(self):
        if self.gc['GMT'] is None:
            utils.echo_error_msg(
                'GMT must be installed to use the SURFACE module'
            )
            return(None, -1)

        out, status = utils.run_cmd(
            'gmt gmtset IO_COL_SEPARATOR = SPACE',
            verbose = False,
            cwd = self.cache_dir
        )        

        #self.gmt_region = self.ps_region.copy()
        dem_surf_cmd = ('')
        if self.blockmean:
            dem_surf_cmd = (
                'gmt blockmean {} -I{:.16f}/{:.16f}+e{}{}{} -V |'.format(
                    self.ps_region.format('gmt') if not self.pixel_node else self.p_region.format('gmt'),
                    self.xinc, self.yinc,
                    ' -W' if self.want_weight else '',
                    ' -fg' if self.geographic else '',
                    ' -rp' if self.pixel_node else '',
                )
            )

        ## mrl: removed -rp and switched p_region to ps_region
        ## (pre 6.5.0 will shift the grid otherwise)
        dem_surf_cmd += (
            'gmt surface -V {} -I{:.16f}/{:.16f}+e -G"{}.tif=gd+n{}:GTiff" -T{} -Z{} {}{}{}{}{}{}{}'.format(
                self.ps_region.format('gmt') if not self.pixel_node else self.p_region.format('gmt'),
                self.xinc, self.yinc,
                self.name, self.ndv, self.tension, self.relaxation,
                ' -Qr' if self.gc['GMT'] >= '6.5.0' else '',
                ' -D{}'.format(self.breakline) if self.breakline is not None else '',
                ' -M{}'.format(self.max_radius) if self.max_radius is not None else '',
                ' -C{}'.format(self.convergence) if self.convergence is not None else '',
                ' -A{}'.format(self.aspect) if self.aspect is not None else '',
                ' -fg' if self.geographic else '',
                ' -rp' if self.pixel_node else '',
            )
        )

        # dem_surf_cmd += (
        #     'gmt surface -V {} -I{:.16f}/{:.16f} -G{}.tif=gd+n{}:GTiff -rp -T{} -Z{} {}{}{}{}{}{}'.format(
        #         self.p_region.format('gmt'), self.xinc, self.yinc,
        #         self.name, self.ndv, self.tension, self.relaxation,
        #         '' if self.gc['GMT'] >= '6.5.0' else '',
        #         ' -D{}'.format(self.breakline) if self.breakline is not None else '',
        #         ' -M{}'.format(self.max_radius) if self.max_radius is not None else '',
        #         ' -C{}'.format(self.convergence) if self.convergence is not None else '',
        #         ' -A{}'.format(self.aspect) if self.aspect is not None else '',
        #         ' -fg' if self.geographic else '',
        #     )
        # )

        out, status = utils.run_cmd(
            dem_surf_cmd,
            verbose=self.verbose,
            data_fun=lambda p: self.dump_xyz(
                dst_port=p, encode=True
            )
        )
            
        return(self)

    
class GMTTriangulate(Waffle):
    """TRIANGULATION DEM via GMT triangulate
    
    Generate a DEM using GMT's triangulate command.
    see gmt triangulate --help for more info.

    < gmt-triangulate >
    """
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)        
        
    def run(self):
        if self.gc['GMT'] is None:
            utils.echo_error_msg(
                'GMT must be installed to use the TRIANGULATE module'
            )
            return(None, -1)

        out, status = utils.run_cmd(
            'gmt gmtset IO_COL_SEPARATOR = SPACE',
            verbose = False
        )        
        
        dem_tri_cmd = 'gmt triangulate -V {} -I{:.14f}/{:.14f} -G{}.tif=gd:GTiff'.format(
            self.ps_region.format('gmt'), self.xinc, self.yinc, self.name
        )
        out, status = utils.run_cmd(
            dem_tri_cmd,
            verbose = self.verbose,
            data_fun = lambda p: self.dump_xyz(
                dst_port=p, encode=True
            )
        )
        
        return(self)

    
class GMTNearNeighbor(Waffle):
    """NEARNEIGHBOR DEM via GMT nearneighbor
    
    Generate a DEM using GMT's nearneighbor command.
    see gmt nearneighbor --help for more info.
    
    -----------
    Parameters:
    
    radius=[val] - search radius
    sectors=[val] - sector information
    
    < gmt-nearneighbor:radius=None:sectors=None >
    """
    
    def __init__(self, radius=None, sectors=None, **kwargs):
        super().__init__(**kwargs) 
        self.radius = radius
        self.sectors = sectors

        
    def run(self):
        if self.gc['GMT'] is None:
            utils.echo_error_msg(
                'GMT must be installed to use the NEARNEIGHBOR module'
            )
            return(None, -1)

        out, status = utils.run_cmd(
            'gmt gmtset IO_COL_SEPARATOR = SPACE',
            verbose = False
        )        

        dem_nn_cmd = 'gmt nearneighbor -V {} -I{:.14f}/{:.14f} -G{}.tif=gd+n{}:GTiff{}{}{}'.format(
            self.ps_region.format('gmt'),
            self.xinc,
            self.yinc,
            self.name,
            self.ndv,
            ' -W' \
            if self.want_weight \
            else '', ' -N{}'.format(self.sectors) \
            if self.sectors is not None \
            else '',
            ' -S{}'.format(self.radius) \
            if self.radius is not None \
            else ' -S{}'.format(self.xinc),
        )
        out, status = utils.run_cmd(
            dem_nn_cmd,
            verbose = self.verbose,
            data_fun = lambda p: self.dump_xyz(
                dst_port=p, encode=True
            )
        )
        
        return(self)


### End
