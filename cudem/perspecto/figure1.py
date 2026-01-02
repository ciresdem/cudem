### figure1.py 
##
## Copyright (c) 2023 - 2026 Regents of the University of Colorado
##
## figure1.py is part of CUDEM
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
## Generate standard figures (2D Map or 3D Perspective) using PyGMT.
##
### Code:

from cudem import utils
from . import gmtimage

try:
    import pygmt
    HAS_PYGMT = True
except ImportError:
    HAS_PYGMT = False


class Figure1(gmtimage.GMTImage):
    """
    Generate Figure 1 (Standard Visualization).
    
    Can produce a 2D contour map or a 3D perspective view.
    
    Configuration Example:
    < figure1:perspective=False:vertical_exaggeration=1.5:interval=100:azimuth=-130:elevation=30 >
    """
    
    def __init__(
        self,
        perspective=False,
        vertical_exaggeration=1.5,
        interval=100,
        azimuth=315,
        elevation=45,
        shade=True,
        want_colorbar=True,
        colorbar_text='Elevation',
        **kwargs
    ):
        super().__init__(**kwargs)
        self.perspective = perspective
        self.vertical_exaggeration = utils.float_or(vertical_exaggeration, 1)
        self.interval = interval
        self.azimuth = azimuth
        self.elevation = elevation
        self.shade = shade
        self.want_colorbar = want_colorbar
        self.colorbar_text = colorbar_text
        self.grad_grid = None

        
    def figure1_perspective(self):
        """Generates a 3D perspective view."""
        
        basename = utils.fn_basename2(self.src_dem)
        
        # Calculate Z-Scale based on latitude length
        if hasattr(utils, 'lll'):
            dem_lll = utils.lll(self.dem_region.ymin)
            deg_lat_len = dem_lll[0]
            
            zscale = (100 * (self.dem_infos['zr'][1] - self.dem_infos['zr'][0])) \
                / (deg_lat_len * (self.dem_region.ymax - self.dem_region.ymin))
        else:
            ## Fallback if lll is missing
            utils.echo_msg("utils.lll missing, using default zscale logic.")
            zscale = 0.1

        utils.echo_msg(f'zscale is {zscale}')
        
        outfile = f'{basename}_persp.png'
        fig = pygmt.Figure()
        
        fig.grdview(
            grid=self.grid,
            frame=['xaf', 'yaf'],
            perspective=[self.azimuth, self.elevation],
            zsize=f'{zscale * self.vertical_exaggeration}c',
            surftype='s',
            cmap=self.cpt,
            contourpen='0.1p',
        )
        
        fig.colorbar(
            perspective=True, 
            frame=["x+lElevation", "y+lm"]
        )
        
        fig.savefig(outfile)
        return outfile

    def figure1_standard(self):
        """Generates a standard 2D contour map."""
        
        basename = utils.fn_basename2(self.src_dem)
        
        ## Generate gradient for shading
        self.grad_grid = pygmt.grdgradient(
            grid=self.grid, 
            azimuth=[self.azimuth, self.elevation], 
            normalize='e0.9'
        )
        
        fig = pygmt.Figure()
        
        fig.grdimage(
            frame=['af', f'+t{basename}'],
            grid=self.grid,
            cmap=self.cpt,
            shading=self.grad_grid if self.shade else None,
        )
        
        fig.grdcontour(
            grid=self.grid,
            interval=self.interval,
            cut=60,
            pen="c0.1p",
        )
        
        if self.want_colorbar:
            fig.colorbar(frame=[f'x+l{self.colorbar_text}', 'y+1m'])

        outfile = f'{basename}_figure1.png'
        fig.savefig(outfile)
        return outfile

    def run(self):
        if not HAS_PYGMT:
            utils.echo_error_msg("PyGMT is not installed.")
            return None

        try:
            if self.perspective:
                result = self.figure1_perspective()
            else:
                result = self.figure1_standard()
            return result
            
        finally:
            ## Cleanup the temporary NetCDF file created in GMTImage.__init__
            ## TODO: update so as not to accidententally delete the acutal
            ## src_dem
            if hasattr(self, 'src_dem_netcdf') and self.src_dem_netcdf:
                utils.remove_glob(self.src_dem_netcdf)

### End
