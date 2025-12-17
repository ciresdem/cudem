### figure1.py 
##
## Copyright (c) 2023 - 2025 Regents of the University of Colorado
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

from cudem.perspecto import gmtimage
from cudem import utils
import pygmt

class figure1(gmtimage.GMTImage):
    """Generate Figure 1

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
        self.vertical_exaggeration=utils.float_or(vertical_exaggeration, 1)
        self.interval = interval
        self.azimuth=azimuth
        self.elevation=elevation
        self.shade=shade
        self.want_colorbar=want_colorbar
        self.colorbar_text=colorbar_text

        
    def figure1_perspective(self):
        dem_lll = lll(self.dem_region.ymin)
        deg_lat_len = dem_lll[0]
        zscale = (100 * (self.dem_infos['zr'][1] - self.dem_infos['zr'][0])) \
            / (deg_lat_len * (self.dem_region.ymax - self.dem_region.ymin))
        utils.echo_msg('zscale is {}'.format(zscale))
        fig = pygmt.Figure()
        fig.grdview(
            grid=self.grid,
            frame=['xaf', 'yaf'],
            perspective=[self.azimuth, self.elevation],
            zsize='{}c'.format(zscale*self.vertical_exaggeration),
            #zscale=zscale,
            surftype='s',
            #plane='{}+ggray'.format(self.dem_infos['zr'][1]*self.dem_infos['zr'][0]/2),
            cmap=self.cpt,
            #shading=self.grad_grid,
            contourpen='0.1p',
        )
        # fig.colorbar(
        #     perspective=True,
        #     frame=["a{}".format(self.interval),
        #            "x+lElevation", "y+lm"]
        # )
        fig.colorbar(perspective=True, frame=["x+lElevation", "y+lm"])
        fig.savefig('{}_persp.png'.format(utils.fn_basename2(self.src_dem)))
        return('{}_persp.png'.format(utils.fn_basename2(self.src_dem)))

    
    def figure1(self):
        self.grad_grid=pygmt.grdgradient(
            grid=self.grid, azimuth=[self.azimuth, self.elevation], normalize='e0.9'
        )
        fig = pygmt.Figure()
        # fig.basemap(
        #     region=self.dem_region.format('str'),
        #     frame=['a', '+t{}'.format(utils.fn_basename2(self.src_dem))],
        # )

        fig.grdimage(
            frame=['af', '+t{}'.format(utils.fn_basename2(self.src_dem))],
            grid=self.grid,
            cmap=self.cpt,
            #cmap=True,
            shading=self.grad_grid if self.shade else None,
            #dpi=100,
        )
        fig.grdcontour(
            grid=self.grid,
            #annotation=(self.interval, '+s2+pfaint'),
            interval=self.interval,
            cut=60,
            pen="c0.1p",
        )
        if self.want_colorbar:
            #fig.colorbar(frame=['a{}'.format(self.interval), 'x+lElevation', 'y+1m'])
            fig.colorbar(frame=['x+l{}'.format(self.colorbar_text), 'y+1m'])

        fig.savefig('{}_figure1.png'.format(utils.fn_basename2(self.src_dem)))
        return('{}_figure1.png'.format(utils.fn_basename2(self.src_dem)))

    
    def run(self):
        if self.perspective:
            return(self.figure1_perspective())
        else:
            return(self.figure1())

        ## this is removing the netcdf, self.src_dem is redefined in __init__
        utils.remove_glob(self.src_dem)

### End
