### gmtgrid.py
##
## Copyright (c) 2010 - 2026 Regents of the University of Colorado
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
### Commentary:
##
## Wrappers for GMT gridding modules (surface, triangulate, nearneighbor).
##
### Code:

from cudem import utils
from cudem import gdalfun
from cudem.waffles.waffles import Waffle

class GMTSurface(Waffle):
    """SPLINE DEM via GMT surface.

    Generate a DEM using GMT's surface command.
    See `gmt surface --help` for more info.

    Parameters:
    -----------
    tension (float) : Spline tension [0-1] (default: 0.35)
    relaxation (float) : Spline relaxation factor
    aspect (float) : Gridding aspect ratio
    breakline (str) : Path to XYZ breakline dataset
    convergence (float) : Gridding convergence limit
    blockmean (bool) : Pipe data through `gmt blockmean` before gridding
    geographic (bool) : Data/grid are geographic (adds -fg)
    pixel_node (bool) : Grid in pixel-node registration (adds -rp)
    """

    def __init__(self, tension=None, relaxation=1, max_radius=None,
                 aspect=None, breakline=None, convergence=None,
                 blockmean=False, geographic=True, pixel_node=False,
                 **kwargs):
        super().__init__(**kwargs)

        ## Validate Tension
        self.tension = utils.float_or(tension, 0.35)
        if self.tension > 1: self.tension = 1

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
            utils.echo_error_msg('GMT must be installed to use the SURFACE module')
            return None, -1

        out, status = utils.run_cmd(
            'gmt gmtset IO_COL_SEPARATOR = SPACE',
            verbose=False,
            cwd=self.cache_dir
        )

        ## Determine Region Format
        region_fmt = self.p_region.format('gmt') if self.pixel_node else self.ps_region.format('gmt')

        ## Build Blockmean Command
        dem_surf_cmd = ''
        if self.blockmean:
            dem_surf_cmd = (
                f'gmt blockmean {region_fmt} -I{self.xinc:.16f}/{self.yinc:.16f}+e '
                f'{"-W" if self.want_weight else ""} '
                f'{"-fg" if self.geographic else ""} '
                f'{"-rp" if self.pixel_node else ""} -V | '
            )

        ## Build Surface Command
        ## Note: -Qr is used for GMT >= 6.5.0
        ## Note: -D (breakline), -M (radius), -C (convergence), -A (aspect)
        dem_surf_cmd += (
            f'gmt surface -V {region_fmt} -I{self.xinc:.16f}/{self.yinc:.16f}+e '
            f'-G"{self.name}.tif=gd+n{self.ndv}:GTiff" '
            f'-T{self.tension} -Z{self.relaxation} '
            f'{"-Qr" if self.gc["GMT"] >= "6.5.0" else ""} '
            f'{f"-D{self.breakline}" if self.breakline else ""} '
            f'{f"-M{self.max_radius}" if self.max_radius else ""} '
            f'{f"-C{self.convergence}" if self.convergence else ""} '
            f'{f"-A{self.aspect}" if self.aspect else ""} '
            f'{"-fg" if self.geographic else ""} '
            f'{"-rp" if self.pixel_node else ""}'
        )

        out, status = utils.run_cmd(
            dem_surf_cmd,
            verbose=self.verbose,
            data_fun=lambda p: self.dump_xyz(dst_port=p, encode=True)
        )

        return self


class GMTTriangulate(Waffle):
    """TRIANGULATION DEM via GMT triangulate.

    Generate a DEM using GMT's triangulate command.
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.gc = utils.config_check(chk_config_file=False)


    def run(self):
        if self.gc['GMT'] is None:
            utils.echo_error_msg('GMT must be installed to use the TRIANGULATE module')
            return None, -1

        out, status = utils.run_cmd(
            'gmt gmtset IO_COL_SEPARATOR = SPACE',
            verbose=False
        )

        dem_tri_cmd = (
            f'gmt triangulate -V {self.ps_region.format("gmt")} '
            f'-I{self.xinc:.14f}/{self.yinc:.14f} '
            f'-G{self.name}.tif=gd:GTiff'
        )

        out, status = utils.run_cmd(
            dem_tri_cmd,
            verbose=self.verbose,
            data_fun=lambda p: self.dump_xyz(dst_port=p, encode=True)
        )

        return self


class GMTNearNeighbor(Waffle):
    """NEARNEIGHBOR DEM via GMT nearneighbor.

    Generate a DEM using GMT's nearneighbor command.

    Parameters:
    -----------
    radius (float) : search radius
    sectors (int) : sector information
    """

    def __init__(self, radius=None, sectors=None, **kwargs):
        super().__init__(**kwargs)
        self.radius = radius
        self.sectors = sectors
        self.gc = utils.config_check(chk_config_file=False)


    def run(self):
        if self.gc['GMT'] is None:
            utils.echo_error_msg('GMT must be installed to use the NEARNEIGHBOR module')
            return None, -1

        out, status = utils.run_cmd(
            'gmt gmtset IO_COL_SEPARATOR = SPACE',
            verbose=False
        )

        ## Determine Radius
        search_radius = f"-S{self.radius}" if self.radius is not None else f"-S{self.xinc}"

        dem_nn_cmd = (
            f'gmt nearneighbor -V {self.ps_region.format("gmt")} '
            f'-I{self.xinc:.14f}/{self.yinc:.14f} '
            f'-G{self.name}.tif=gd+n{self.ndv}:GTiff '
            f'{"-W" if self.want_weight else ""} '
            f'{f"-N{self.sectors}" if self.sectors else ""} '
            f'{search_radius}'
        )

        out, status = utils.run_cmd(
            dem_nn_cmd,
            verbose=self.verbose,
            data_fun=lambda p: self.dump_xyz(dst_port=p, encode=True)
        )

        return self

### End
