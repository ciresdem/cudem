### perspecto.py 
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
## Generate iMAGEs from a DEm
##
## uses:
##   povray
##   gdal
##   ImageMagick
##   GMT
##
### Code:

import os
import math
import sys

from osgeo import gdal
from osgeo_utils import gdal_calc
import numpy as np
from tqdm import tqdm
import colorsys

from cudem import utils
from cudem import gdalfun
from cudem import regions
from cudem import fetches
from cudem import factory

try:
    import pygmt
    has_pygmt = True
except ImportError as e:
    has_pygmt = False

    
## lll
def lll(src_lat):
    gds_equator = 111321.543
    gds_pi = 3.14159265358979323846
    degree_to_radian = lambda d: gds_pi * (d / 180)
    lonl_ = math.cos(degree_to_radian(src_lat)) * gds_equator
    latl_ = gds_equator

    return(lonl_, latl_)
                    

class Perspecto:
    def __init__(
            self,
            mod=None,
            src_dem=None,
            cpt=None,
            min_z=None,
            max_z=None,
            callback=lambda: False,
            outfile=None,
            outdir=None,
            verbose=True,
            split_cpt=False,
            params={}
    ):
        self.mod = mod
        self.mod_args = {}
        self.src_dem = src_dem
        self.outfile = outfile
        self.cpt = cpt
        self.min_z = utils.float_or(min_z)
        self.max_z = utils.float_or(max_z)
        self.split_cpt = split_cpt
        self.outdir = outdir
        self.callback = callback
        self.verbose = verbose
        self.params = params    
        self.dem_infos = gdalfun.gdal_infos(self.src_dem, scan=True)
        self.dem_region = regions.Region().from_geo_transform(
            self.dem_infos['geoT'], self.dem_infos['nx'], self.dem_infos['ny']
        )
        if self.outfile is None:
            self.outfile = '{}_{}.tif'.format(
                utils.fn_basename2(self.src_dem), self.mod
            )
        
        self.init_cpt(want_gdal=True)

        
    def __call__(self):
        return(self.run())

    
    def init_cpt(self, want_gdal=False):
        min_z = self.min_z if self.min_z is not None else self.dem_infos['zr'][0]
        max_z = self.max_z if self.max_z is not None else self.dem_infos['zr'][1]
        if self.cpt is None:
            # self.makecpt(
            #     'etopo1',
            #     output='{}_etopo1.cpt'.format(utils.fn_basename2(self.src_dem))
            # )
            self.cpt = generate_etopo_cpt(min_z, max_z)
        #else:
        #self.cpt = process_cpt(self.cpt, min_z, max_z)
        elif os.path.exists(self.cpt):
            #if has_pygmt:
            #    self.makecpt(cmap=self.cpt, color_model='r', output='{}.cpt'.format(utils.fn_basename2(self.src_dem)))
            #else:
            utils.echo_msg(
                'processing cpt {}, want_gdal is {}, split_cpt: {}'.format(
                    self.cpt, want_gdal, self.split_cpt
                )
            )
            self.cpt = process_cpt(
                self.cpt, min_z, max_z, gdal=want_gdal, split_cpt=self.split_cpt
            )
        else:
            self.cpt = process_cpt(
                fetch_cpt_city(q=self.cpt),
                min_z,
                max_z,
                gdal=want_gdal,
                split_cpt=self.split_cpt
            )
        
    
    def makecpt(self, cmap='etopo1', color_model='r', output=None):
        min_z = self.min_z if self.min_z is not None else self.dem_infos['zr'][0]
        max_z = self.max_z if self.max_z is not None else self.dem_infos['zr'][1]
        pygmt.makecpt(
            cmap=cmap,
            color_model='{}'.format(color_model),
            output=output,
            series=[min_z, max_z,50],
            no_bg=True,
            #continuous=True,
            #truncate=[min_z, max_z],
        )
        self.cpt = output

        
    def cpt_no_slash(self):
        # Read in the file
        with open(self.cpt, 'r') as file:
            filedata = file.read()
            
        # Replace the target string
        filedata = filedata.replace('\\', ' ')
            
        # Write the file out again
        with open(self.cpt, 'w') as file:
            file.write(filedata)

            
    def lll(self):
        gds_equator = 111321.543
        gds_pi = 3.14159265358979323846
        degree_to_radian = lambda d: gds_pi * (d / 180)
        lonl_ = math.cos(degree_to_radian(self.dem_region.ymin)) * gds_equator
        latl_ = gds_equator

        return(lonl_, latl_)

    
    def export_as_png(self, rgb=True, dem=True):
        if dem:
            utils.run_cmd(
                'gdal_translate -ot UInt16 -of PNG -scale {} {} 0 65535 {} _dem_temp.png'.format(
                    self.dem_infos['zr'][0],
                    self.dem_infos['zr'][1],
                    self.src_dem
                ),
                verbose=True
            )
            utils.run_cmd(
                'gdal_translate -srcwin 1 1 {} {} -of PNG _dem_temp.png {}_16bit.png'.format(
                    self.dem_infos['nx']-1,
                    self.dem_infos['ny']-1,
                    utils.fn_basename2(self.src_dem)
                ),
                verbose=True
            )
            utils.remove_glob('_dem_temp*')
            
        if rgb:
            self.init_cpt(want_gdal=True)
            #cr_fn = utils.make_temp_fn('gdaldem_cr.tif', self.outdir)
            # gdal.DEMProcessing(
            #     cr_fn, self.src_dem, 'color-relief', colorFilename=self.cpt,
            #     computeEdges=True, addAlpha=self.alpha
            # )
            
            utils.run_cmd(
                'gdaldem color-relief {} {} _rgb_temp.tif -alpha'.format(
                    self.src_dem, self.cpt
                ),
                verbose=True
            )
            utils.run_cmd(
                'gdal_translate -srcwin 1 1 {} {} -of PNG _rgb_temp.tif {}_rgb.png'.format(
                    self.dem_infos['nx']-1,
                    self.dem_infos['ny']-1,
                    utils.fn_basename2(self.src_dem)
                ),
                verbose=True
            )
            utils.remove_glob('_rgb_temp*')

            
    def run(self):
        raise(NotImplementedError)
             
### End
    
