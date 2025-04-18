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

## CPT          
def scale_el__(value, gmin, gmax, tr):
    if value > 0 and gmax > 0:
        return((gmax * tr) / 8000)
    elif value < 0 and gmin < 0:
        return((gmin * tr) / -11000)
    elif value == 0:
        return(0)
    else:
        print(value)
        return(None)

def scale_el_(value, gmin, gmax, tr, trs):
    if value > 0 and gmax > 0:
        return((gmax * tr) / max(trs))
    elif value < 0 and gmin < 0:
        if min(trs) == 0:
            return(gmin * tr)
        else:
            return((gmin * tr) / min(trs))
    elif value == 0:
        return(0)
    else:
        return(None)

def scale_el(value, gmin, gmax, tr, trs):
    p = (tr - min(trs)) / (max(trs) - min(trs))
    v = (1-p) * (gmin - gmax) + gmax
    #v = p * (gmax - gmin) + gmin
    
    # if value < 0:
    #     p = (tr - min(trs)) / (max(trs) - min(trs))
    #     v = (1-p) * (gmin - 0) + 0
    # elif value > 0:
    #     p = (tr - min(trs)) / (max(trs) - min(trs))
    #     v = (1-p) * (0 - gmax) + gmax
    # elif value == 0:
    #     v = 0

    return(v)
    
def generate_etopo_cpt(gmin, gmax):

    trs = [
        -11000, -10500, -10000, -9500, -9000, -8500,
        -8000, -7500, -7000, -6500, -6000, -5500,
        -5000, -4500, -4000, -3500, -3000, -2500,
        -2000, -1500, -1000, -500, -0.001, 0,
        100, 200, 500, 1000, 1500, 2000, 2500,
        3000, 3500, 4000, 4500, 5000, 5500,
        6000, 6500, 7000, 7500, 8000
    ]
    
    elevs = [
        -20, -19.09090909, -18.18181818, -17.27272727, -16.36363636,
        -15.45454545, -14.54545455, -13.63636364, -12.72727273, -11.81818182,
        -10.90909091, -10, -9.090909091, -8.181818182, -7.272727273,
        -6.363636364, -5.454545455, -4.545454545, -3.636363636, -2.727272727,
        -1.818181818, -0.909090909, -0.001, 0, 48.06865898,
        96.13731797, 240.3432949, 480.6865898, 721.0298848, 961.3731797,
        1201.716475, 1442.05977, 1682.403064, 1922.746359, 2163.089654,
        2403.432949, 2643.776244, 2884.119539, 3124.462834, 3364.806129,
        3605.149424,  3845.492719
    ]

    colors = (
        [10,0,121], [26,0,137], [38,0,152], [27,3,166], [16,6,180],
        [5,9,193], [0,14,203], [0,22,210], [0,30,216], [0,39,223],
        [12,68,231], [26,102,240], [19,117,244], [14,133,249], [21,158,252],
        [30,178,255], [43,186,255], [55,193,255], [65,200,255], [79,210,255],
        [94,223,255], [138,227,255], [138,227,255], [51,102,0], [51,204,102],
        [187,228,146], [255,220,185], [243,202,137], [230,184,88], [217,166,39],
        [168,154,31], [164,144,25], [162,134,19], [159,123,13], [156,113,7],
        [153,102,0], [162,89,89], [178,118,118], [183,147,147], [194,176,176],
        [204,204,204], [229,229,229], [138,227,255], [51,102,0]
    )

    with open('tmp.cpt', 'w') as cpt:
        for i,j in enumerate(elevs):
            if j != None and i + 1 != len(elevs):
                elev = scale_el_(j, gmin, gmax, trs[i], trs)
                elev1 = scale_el_(elevs[i + 1], gmin, gmax, trs[i + 1], trs)
                if elev is not None and elev1 is not None:
                    cpt.write(
                        '{} {} {} {} {} {} {} {}\n'.format(
                            elev, colors[i][0], colors[i][1], colors[i][2],
                            elev1, colors[i][0], colors[i][1], colors[i][2]
                        )
                    )
    return('tmp.cpt')

cpt_colors = {
    'black': [0, 0, 0],
    'white': [255, 255, 255],
}

def process_cpt(cpt, gmin, gmax, gdal = False, split_cpt = None):
    trs = []
    colors = []
    with open(cpt, 'r') as cpt:
        for l in cpt:
            ll = l.split()
            if utils.float_or(ll[0]) is not None:
                trs.append(float(ll[0]))

                if utils.int_or(ll[1]) is not None:
                    colors.append(
                        [int(float(ll[1])), int(float(ll[2])), int(float(ll[3]))]
                    )
                elif ll[1] in cpt_colors.keys():
                    colors.append(cpt_colors[ll[1]])
                elif len(ll[1].split('/')) > 1:
                    colors.append([int(float(x)) for x in ll[1].split('/')])

    if split_cpt is not None:
        _trs = np.array(trs)
        len_b = len(_trs[_trs<split_cpt])
        len_t = len(_trs[_trs>split_cpt])
        
        elevs_b = np.linspace(gmin, split_cpt, len_b)
        elevs_t = np.linspace(split_cpt, gmax, len_t)[1:]
        elevs = np.concatenate((elevs_b,elevs_t))
    else:
        elevs = np.linspace(gmin, gmax, len(trs))
        
    with open('tmp.cpt', 'w') as cpt:
        for i, j in enumerate(elevs):
            if j != None and i + 1 != len(elevs):
                elev = scale_el_(j, gmin, gmax, trs[i], trs)
                elev1 = scale_el_(elevs[i + 1], gmin, gmax, trs[i + 1], trs)
                if elev is not None and elev1 is not None:
                    if not gdal:
                        cpt.write(
                            '{} {} {} {} {} {} {} {}\n'.format(
                                elev, colors[i][0], colors[i][1], colors[i][2],
                                elev1, colors[i][0], colors[i][1], colors[i][2]
                            )
                        )
                    else:
                        cpt.write(
                            '{} {} {} {} 255\n'.format(
                                elev, colors[i][0], colors[i][1], colors[i][2]
                            )
                        )
        if gdal:
            cpt.write('nv 0 0 0 0\n')
    return('tmp.cpt')
    
def get_correctMap(path, luminosity, contrast):
        ds = gdal.Open(image_path)
        #To normalize
        band1 = ds.GetRasterBand(1)
        #Get the max value
        maxValue = int(2**16 -1)
        if band1.DataType == gdal.GDT_UInt16:
            maxValue = int(2**16 -1)
        elif band1.DataType == gdal.GDT_Byte:
            maxValue = int(2**8 -1)
        else:
            LOGGER.info(f"band type {band1.DataType} not handled: use default size of value (16 bits)")

        band1 = ds.ReadAsArray(0,0,ds.RasterXSize,ds.RasterYSize)[0]
        band2 = ds.ReadAsArray(0,0,ds.RasterXSize,ds.RasterYSize)[1]
        band3 = ds.ReadAsArray(0,0,ds.RasterXSize,ds.RasterYSize)[2] 

        for x in range(0,ds.RasterXSize):
                for y in range(0,ds.RasterXSize):

                    r = float(band1[x,y]) / maxValue
                    g = float(band2[x,y]) / maxValue
                    b = float(band3[x,y]) / maxValue

                    #Convert to HLS them apply luminosity and contrast
                    (h,l,s) = colorsys.rgb_to_hls(r, g, b)

                    l = min(max(0, l + (l - 0.5)*(luminosity - 0.5)) , 1)
                    s = min(max(0, s + (s - 0.5)*(contrast - 0.5)) , 1)

                    (r,g,b) = colorsys.hls_to_rgb(h, l, s)
                    
def fetch_cpt_city(q = 'grass/haxby', cache_dir = './'):
    utils.echo_msg('checking for `{}` at cpt-city'.format(q))
    this_fetches = fetches.FetchesFactory(
        mod='cpt_city', verbose=True, outdir=cache_dir, q=q
    )._acquire_module()
    this_fetches.run()
    utils.echo_msg(
        'found {} results for `{}` at cpt-city, using first entry `{}`'.format(
            len(this_fetches.results), q, this_fetches.results[0]['url']
        )
    )
    fetches.Fetch(
        this_fetches.results[0]['url']
    ).fetch_file(
        this_fetches.results[0]['dst_fn']
    )
    return(this_fetches.results[0]['dst_fn'])

class Perspecto:
    def __init__(
            self,
            mod = None,
            src_dem = None,
            cpt = None,
            min_z = None,
            max_z = None,
            callback = lambda: False,
            outfile = None,
            outdir = None,
            verbose = True,
            split_cpt = False,
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
            #self.makecpt('etopo1', output='{}_etopo1.cpt'.format(utils.fn_basename2(self.src_dem)))
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
                fetch_cpt_city(q=self.cpt), min_z, max_z, gdal=want_gdal, split_cpt=self.split_cpt
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
                    self.dem_infos['zr'][0], self.dem_infos['zr'][1], self.src_dem
                ),
                verbose=True
            )
            utils.run_cmd(
                'gdal_translate -srcwin 1 1 {} {} -of PNG _dem_temp.png {}_16bit.png'.format(
                    self.dem_infos['nx']-1, self.dem_infos['ny']-1, utils.fn_basename2(self.src_dem)
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
                'gdaldem color-relief {} {} _rgb_temp.tif -alpha'.format(self.src_dem, self.cpt),
                verbose=True
            )
            utils.run_cmd(
                'gdal_translate -srcwin 1 1 {} {} -of PNG _rgb_temp.tif {}_rgb.png'.format(
                    self.dem_infos['nx']-1, self.dem_infos['ny']-1, utils.fn_basename2(self.src_dem)
                ),
                verbose=True
            )
            utils.remove_glob('_rgb_temp*')
        
    def run(self):
        raise(NotImplementedError)

class Hillshade(Perspecto):
    """Generate a Hillshade Image

    https://en.wikipedia.org/wiki/Blend_modes#Overlay

< hillshade:vertical_exaggeration=1:projection=4326:azimuth=315:altitude=45 >
    """
    
    def __init__(
            self,
            vertical_exaggeration=1,
            projection=4326,
            azimuth=315,
            altitude=45,
            modulate=125,
            alpha=False,
            mode='multiply',
            gamma=.5,
            **kwargs,
    ):
        super().__init__(mod='hillshade', **kwargs)
        self.vertical_exaggeration = vertical_exaggeration
        self.projection = projection
        self.azimuth = azimuth
        self.altitude = altitude
        self.modulate = utils.float_or(modulate, 115)
        self.alpha = alpha
        self.mode = mode
        self.gamma = utils.float_or(gamma)
        self.cpt_no_slash()

    def _modulate(self, gdal_fn):
        utils.run_cmd(
            'mogrify -modulate {} -depth 8 {}'.format(self.modulate, gdal_fn),
            verbose=False
        )
        
    def gamma_correction_calc(self, hs_fn = None, gamma = .5, outfile = 'gdaldem_gamma.tif'):
        gdal_calc.Calc(
            "uint8(((A / 255.)**(1/0.5)) * 255)",
            A=hs_fn,
            outfile=outfile,
            quiet=True
        )
        return(outfile)

    def gamma_correction(self, arr):
        if self.gamma is not None:
            gc_arr = (((arr/255.)**(1/self.gamma)) * 255.).astype(np.uint8)
            return(gc_arr)
        else:
            return(arr)
    
    def blend(self, hs_file, rgb_file, mode = 'multiply', gamma = None, outfile = 'gdaldem_multiply.tif'):
        modes = ['multiply', 'screen', 'overlay', 'hard_light', 'soft_light']
        if mode not in modes:
            mode = 'multiply'
            
        with gdalfun.gdal_datasource(rgb_file, update=True) as rgb_ds:
            with gdalfun.gdal_datasource(hs_file) as cr_ds:
                ds_config = gdalfun.gdal_infos(cr_ds)
                n_chunk = int(ds_config['nx'] * .1)
                for srcwin in utils.yield_srcwin(
                        (ds_config['ny'], ds_config['nx']), n_chunk=n_chunk,
                        msg='Blending rgb and hillshade using {}'.format(mode),
                        end_msg='Generated color hillshade',
                        verbose=self.verbose
                ):
                    cr_band = cr_ds.GetRasterBand(1)
                    cr_arr = cr_band.ReadAsArray(*srcwin)
                    cr_arr = self.gamma_correction(cr_arr)
                    for band_no in [1,2,3]:
                        band = rgb_ds.GetRasterBand(band_no)
                        band_arr = band.ReadAsArray(*srcwin)
                        if mode == 'multiply':
                            band_arr = (((cr_arr/255.) * (band_arr/255.)) * 255.).astype(np.uint8)
                        elif mode == 'screen':
                            band_arr = ((1 - ((1 - cr_arr/255.) * (1 - (band_arr/255.)))) * 255.).astype(np.uint8)
                        elif mode == 'overlay':
                            out_arr = np.copy(band_arr)
                            out_arr[cr_arr<128] = (2 * (cr_arr[cr_arr<128]/255.) * (band_arr[cr_arr<128]/255.) *
                                                   255.).astype(np.uint8)
                            out_arr[cr_arr>=128] = ((1 - (2 * (1 - (cr_arr[cr_arr>=128]/255.)) *
                                                          (1 - (band_arr[cr_arr>=128]/255.))))
                                                    * 255.).astype(np.uint8)
                            band_arr[:] = out_arr
                        elif mode == 'hard_light':
                            out_arr = np.copy(band_arr)
                            out_arr[band_arr<128] = (2 * (cr_arr[band_arr<128]/255.) * (band_arr[band_arr<128]/255.)
                                                     * 255.).astype(np.uint8)
                            out_arr[band_arr>=128] = ((1 - (2 * (1 - (cr_arr[band_arr>=128]/255.)) *
                                                            (1 - (band_arr[band_arr>=128]/255.))))
                                                      * 255.).astype(np.uint8)
                            band_arr[:] = out_arr
                        elif mode == 'soft_light':
                            out_arr = np.copy(band_arr)
                            out_arr[band_arr<128] = (((2 * (cr_arr[band_arr<128]/255.) * (band_arr[band_arr<128]/255.)) +
                                                      (((cr_arr[band_arr<128]/255.)**2) *
                                                       (1 - (2 * band_arr[band_arr<128]/255.))))
                                                     * 255.).astype(np.uint8)
                            out_arr[band_arr>=128] = (((2 * (cr_arr[band_arr>=128]/255.) * (1 - band_arr[band_arr>=128]/255.)) +
                                                       (np.sqrt(cr_arr[band_arr>=128]/255.) *
                                                        ((2 * (band_arr[band_arr>=128]/255.)) - 1)))
                                                      * 255.).astype(np.uint8)
                            band_arr[:] = out_arr
                            
                        band.WriteArray(band_arr, srcwin[0], srcwin[1])
                        
        os.rename(rgb_file, outfile)
        return(outfile)
    
    def run(self):
        #self.init_cpt(want_gdal=True)
        hs_fn = utils.make_temp_fn('gdaldem_hs.tif', self.outdir)
        gdal.DEMProcessing(
            hs_fn, self.src_dem, 'hillshade', computeEdges=True, scale=111120,
            azimuth=self.azimuth, altitude=self.altitude, zFactor=self.vertical_exaggeration
        )

        cr_fn = utils.make_temp_fn('gdaldem_cr.tif', self.outdir)
        gdal.DEMProcessing(
            cr_fn, self.src_dem, 'color-relief', colorFilename=self.cpt,
            computeEdges=True, addAlpha=self.alpha
        )

        cr_hs_fn = utils.make_temp_fn('gdaldem_cr_hs.tif', self.outdir)
        self.blend(hs_fn, cr_fn, gamma=self.gamma, mode=self.mode, outfile=cr_hs_fn)

        utils.remove_glob(hs_fn, cr_fn)
        os.rename(cr_hs_fn, self.outfile)
        #self._modulate(self.outfile)
        
        return(self.outfile)

class Hillshade_cmd(Perspecto):
    """Generate a Hillshade Image

    uses gdal/ImageMagick

< hillshade:vertical_exaggeration=1:projection=4326:azimuth=315:altitude=45 >
    """
    
    def __init__(
            self,
            vertical_exaggeration=1,
            projection=None,
            azimuth=315,
            altitude=45,
            **kwargs,
    ):
        super().__init__(**kwargs)
        self.vertical_exaggeration = vertical_exaggeration
        self.projection = projection
        self.azimuth = azimuth
        self.altitude = altitude
        self.cpt_no_slash()
        
    def run(self):
        utils.run_cmd(
            'gdaldem hillshade -compute_edges -s 111120 -z {} -az {} -alt {} {} hillshade.tif'.format(
                self.vertical_exaggeration, self.azimuth, self.altitude, self.src_dem
            ), verbose=self.verbose
        )
        # Generate the color-releif
        utils.run_cmd(
            'gdaldem color-relief {} {} colors.tif'.format(self.src_dem, self.cpt),
            verbose=self.verbose
        )
        # Composite the hillshade and the color-releif
        utils.run_cmd(
            'composite -compose multiply -depth 8 colors.tif hillshade.tif output.tif',
            verbose=self.verbose
        )
        utils.run_cmd(
            'mogrify -modulate 115 -depth 8 output.tif',
            verbose=self.verbose
        )
        # Generate the combined georeferenced tif
        utils.run_cmd(
            'gdal_translate -co "TFW=YES" {} temp.tif'.format(self.src_dem),
            verbose=self.verbose
        )
        utils.run_cmd('mv temp.tfw output.tfw')
        utils.run_cmd(
            'gdal_translate {} output.tif temp2.tif'.format(
                '-a_srs epsg:{}'.format(self.projection) if self.projection is not None else ''
            )
        )
        # Cleanup
        utils.remove_glob(
            'output.tif*', 'temp.tif*', 'hillshade.tif*', 'colors.tif*', 'output.tfw*'
        )
        #subtract 2 cells from rows and columns 
        #gdal_translate -srcwin 1 1 $(gdal_rowcol.py $ingrd t) temp2.tif $outgrd
        utils.run_cmd(
            'gdal_translate temp2.tif {}_hs.tif'.format(utils.fn_basename2(self.src_dem))
        )
        utils.remove_glob('temp2.tif')

        return('{}_hs.tif'.format(utils.fn_basename2(self.src_dem)))
    
class HillShade_test(Perspecto):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def hillshade(self, array, azimuth, angle_altitude):
        azimuth = 360.0 - azimuth 
    
        x, y = np.gradient(array)
        slope = np.pi/2. - np.arctan(np.sqrt(x*x + y*y))
        aspect = np.arctan2(-x, y)
        azm_rad = azimuth*np.pi/180. #azimuth in radians
        alt_rad = angle_altitude*np.pi/180. #altitude in radians        
        shaded = np.sin(alt_rad) * np.sin(slope) + np.cos(alt_rad)*np.cos(slope)*np.cos((azm_rad - np.pi/2.) - aspect)
        
        return(255*(shaded + 1)/2)
    
class POVRay(Perspecto):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.cpt_no_slash()
        
        self.rgb_image = '{}_rgb.png'.format(utils.fn_basename2(self.src_dem))
        self.dem_image = '{}_16bit.png'.format(utils.fn_basename2(self.src_dem))

        #if not os.path.exists(self.rgb_image) or not os.path.exists(self.dem_image):
        self.export_as_png()            
        self.output_pov = '{}.pov'.format(utils.fn_basename2(self.src_dem))
        
    def run_povray(self, src_pov_template, pov_width=800, pov_height=800):
        utils.run_cmd(
            'povray {} +W{} +H{} -D'.format(src_pov_template, pov_width, pov_height),
            verbose=True
        )

class perspective(POVRay):
    """Generate a perspective image

< perspective:cam_azimuth=-130:cam_elevation=30:cam_distance=265:cam_view_angle=40:light_elevation=20:light_distance=10000:vertical_exaggeration=2 >
    """
    
    def __init__(
            self,
            cam_azimuth=-130,
            cam_elevation=30,
            cam_distance=265,
            cam_view_angle=40,
            light_elevation=20,
            light_distance=10000,
            vertical_exaggeration=2,
            **kwargs
    ):
        super().__init__(**kwargs)
        self.cam_azimuth = cam_azimuth
        self.cam_elevation = cam_elevation
        self.cam_distance = cam_distance
        self.cam_view_angle = cam_view_angle
        self.light_elevation = light_elevation
        self.light_distance = light_distance
        self.vertical_exaggeration = vertical_exaggeration
        self.dem_lll = lll(self.dem_region.ymin)
        self.output_pov = '{}_perspective.pov'.format(utils.fn_basename2(self.src_dem))
        
        self.template = """
// DEM

global_settings {{ assumed_gamma 1 }}
#include \"colors.inc\" 

//
// Custom parameters start here
//
#declare rgb_image = \"rgb.png\"
#declare dem_image = \"dem_16bit.png\"

#declare xdim = {xdim};  //number of pixels in X-direction
#declare ydim = {ydim};  //number of pixels in y-direction
#declare max_y = {ymax}; //maximum latitude extent
#declare min_y = {ymin}; //minimum latitude extent
#declare min_z = {zmin}; //minimum elevation
#declare max_z = {zmax}; //maximum elevation

// Obtained from http://www.csgnetwork.com/degreelenllavcalc.html  
#declare deg_lat_len = {deg_lat_len}; //length of a degree of latitude in meters  
#declare deg_lon_len = {deg_lon_len}; //length of a degree of longitude in meters

// Position of camera             
#declare cam_azimuth = {cam_azimuth};
#declare cam_elevation = {cam_elevation};
#declare cam_distance = {cam_distance}; 
#declare cam_view_angle = {cam_view_angle};
                     
// Position of the \"sun\"  
#declare light_azimuth = cam_azimuth+90;
#declare light_elevation = {light_elevation};
#declare light_distance = {light_distance};             

#declare vertical_exaggeration = {vertical_exaggeration};
//
// Custom parameters end here
//
             
#declare lon_scale = deg_lon_len / deg_lat_len;
#declare z_scale = (100 * (max_z - min_z)) / (deg_lat_len * (max_y - min_y));

#declare cam_x = cam_distance * cos(radians(cam_elevation)) * sin(radians(cam_azimuth));
#declare cam_y = cam_distance * sin(radians(cam_elevation));
#declare cam_z = cam_distance * cos(radians(cam_elevation)) * cos(radians(cam_azimuth));
#declare light_x = light_distance * cos(radians(light_elevation)) * sin(radians(light_azimuth));
#declare light_y = light_distance * sin(radians(light_elevation));
#declare light_z = light_distance * cos(radians(light_elevation)) * cos(radians(light_azimuth));

#declare Texture0 = // Color elevation image (24-bit RGB PNG)
texture {{
  pigment{{
    image_map {{ 
      png "{rgb_image}" map_type 0 once interpolate 2
    }} 
  }} 
  finish {{ ambient 0.4 diffuse 0.8 }}
  rotate x*90  
}}
                  
height_field {{ // Unsigned 16-bit PNG DEM
   png "{dem_image}"
   smooth
   clipped_by {{box {{ <0, 0, 0>, <0.999, 1, 0.999> }} }}
   texture {{ Texture0 }}
   translate <-0.5, 0, -0.5>
   scale <100*lon_scale*{xdim}/{ydim},
          vertical_exaggeration*z_scale,  //Vertical exaggeration
          100>
}} 

camera {{
   angle cam_view_angle
   location <cam_x, cam_y, cam_z>
   look_at <0, 0, 0> 
}}
                                    
light_source {{ <light_x, light_y, light_z> color White shadowless parallel }}
        
background {{ White }}""".format(
    rgb_image=self.rgb_image,
    dem_image=self.dem_image,
    xdim=self.dem_infos['nx'],
    ydim=self.dem_infos['ny'],
    ymax = self.dem_region.ymax,
    ymin = self.dem_region.ymin,
    zmax = self.dem_infos['zr'][1],
    zmin = self.dem_infos['zr'][0],
    deg_lat_len = self.dem_lll[0],
    deg_lon_len = self.dem_lll[1],
    cam_azimuth = self.cam_azimuth,
    cam_elevation = self.cam_elevation,
    cam_distance = self.cam_distance,
    cam_view_angle = self.cam_view_angle,
    light_elevation = self.light_elevation,
    light_distance = self.light_distance,
    vertical_exaggeration = self.vertical_exaggeration,
)

    def run(self):
        with open(self.output_pov, 'w') as pov_temp:
            pov_temp.write(self.template)
        
        self.run_povray(self.output_pov, self.dem_infos['nx'], self.dem_infos['ny'])

class sphere(POVRay):
    """Generate a sphere image

< sphere:cam_azimuth=310:cam_elevation=27:cam_distance=8:cam_view_angle=33:center_lat=None:center_long=None >
    """
    
    def __init__(
            self,
            cam_azimuth=310,
            cam_elevation=27,
            cam_distance=8,
            cam_view_angle=33,
            center_lat=None,
            center_long=None,
            **kwargs
    ):
        super().__init__(**kwargs)
        self.cam_azimuth = cam_azimuth
        self.cam_elevation = cam_elevation
        self.cam_distance = cam_distance
        self.cam_view_angle = cam_view_angle

        if center_lat is None or center_long is None:
            region_center = self.dem_region.center()
            self.center_long = region_center[0]
            self.center_lat = region_center[1]
        else:
            self.center_long = center_long
            self.center_lat = center_lat

        self.output_pov = '{}_sphere.pov'.format(utils.fn_basename2(self.src_dem))
            
        self.template = """
//Generates a hillshade color image on a sphere 
//CUDEM 2023MAR08

#include "colors.inc"  

#declare cam_azimuth = {cam_azimuth};
#declare cam_elevation = {cam_elevation};
#declare cam_distance = {cam_distance}; 
#declare cam_view_angle = {cam_view_angle};
        
#declare cam_x = cam_distance * cos(radians(cam_elevation)) * sin(radians(cam_azimuth));
#declare cam_y = cam_distance * sin(radians(cam_elevation));
#declare cam_z = cam_distance * cos(radians(cam_elevation)) * cos(radians(cam_azimuth));
        
// Colormap image  
#declare colorImage = 
pigment {{
  image_map {{
    png "{rgb_image}"     // Specify color map image
    map_type 1        // 0=planar, 1=spherical, 2=cylindrical, 5=torus
    interpolate 2     // 0=none, 1=linear, 2=bilinear, 4=normalized distance
  }} // image_map
}} 

// ETOPO1 topography used as a bump map
#declare hillShade = 
normal {{
  bump_map {{             // uses image color or index as bumpiness
    png "{dem_image}"  // the file to read (tiff/tga/gif/png/jpeg/tiff/sys)
    map_type 1           // 0=planar, 1=spherical, 2=cylindrical, 5=torus
    interpolate 2        // 0=none, 1=linear, 2=bilinear, 4=normalized distance
    bump_size 750        // Adjust the vertical exaggeration (typical values between 500 and 1000)
  }} // bump_map
}}

// Apply pigment to the texture 
#declare colorShade =
texture {{
  pigment {{
    colorImage
  }}
  finish {{ambient 0.3 diffuse 1}}
  normal {{ hillShade }}
}}

//Adjust values to desired cooridinates
#declare longitude = {center_long};
#declare latitude = {center_lat};

sphere {{
  <0,0,0>, 1
  texture {{ colorShade }}
  //rotate y*90 
  no_shadow
  //rotate y*longitude
  //rotate x*latitude
  //rotate x*20
}}

camera {{
   angle cam_view_angle
   location <cam_x, cam_y, cam_z>
   look_at <0, 0, 0> 
}}""".format(
    rgb_image=self.rgb_image,
    dem_image=self.dem_image,
    center_long=self.center_long,
    center_lat=self.center_lat,
    cam_azimuth = self.cam_azimuth,
    cam_elevation = self.cam_elevation,
    cam_distance = self.cam_distance,
    cam_view_angle = self.cam_view_angle,
)

    def run(self):

        with open(self.output_pov, 'w') as pov_temp:
            pov_temp.write(self.template)
        
        self.run_povray(
            self.output_pov, self.dem_infos['nx'], self.dem_infos['ny']
        )

class GMTImage(Perspecto):
    """Use GMTImage with pygmt"""
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        if self.dem_infos['fmt'] != 'NetCDF':
            utils.run_cmd(
                'gmt grdconvert {} {}.nc'.format(
                    self.src_dem, utils.fn_basename2(self.src_dem)
                )
            )
            self.src_dem = '{}.nc'.format(utils.fn_basename2(self.src_dem))
        
        self.grid = pygmt.load_dataarray(self.src_dem)
        self.makecpt(self.cpt, output=None)

class colorbar(GMTImage):
    """Generate a colorbar

< colorbar:colorbar_text='Elevation':widt=10:height=2 >
    """
    
    def __init__(self, colorbar_text = 'Elevation', width = 10, height = 2, **kwargs):
        super().__init__(**kwargs)
        self.colorbar_text=colorbar_text
        self.width = utils.int_or(width)
        self.height = utils.int_or(height)
        
    def run(self):
        fig = pygmt.Figure()
        fig.colorbar(
            region=[0,10,0,3],
            projection="X10c/3c",
            position='jTC+w{}c/{}c+h'.format(self.width, self.height),
            frame=['x+l{}'.format(self.colorbar_text), 'y+1m']
        )
        fig.savefig('{}_cbar.png'.format(utils.fn_basename2(self.src_dem)))
        return('{}_cbar.png'.format(utils.fn_basename2(self.src_dem)))
                
class figure1(GMTImage):
    """Generate Figure 1

< figure1:perspective=False:vertical_exaggeration=1.5:interval=100:azimuth=-130:elevation=30 >
    """
    
    def __init__(
            self,
            perspective=False,
            vertical_exaggeration=1.5,
            interval=100,
            azimuth=-130,
            elevation=30,
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
        #fig.colorbar(perspective=True, frame=["a{}".format(self.interval), "x+lElevation", "y+lm"])
        fig.colorbar(perspective=True, frame=["x+lElevation", "y+lm"])
        fig.savefig('{}_persp.png'.format(utils.fn_basename2(self.src_dem)))
        return('{}_persp.png'.format(utils.fn_basename2(self.src_dem)))
        
    def figure1(self):
        self.grad_grid=pygmt.grdgradient(
            grid=self.grid, azimuth=[self.azimuth, self.elevation], normalize='e0.6'
        )
        fig = pygmt.Figure()
        # fig.basemap(
        #     region=self.dem_region.format('str'),
        #     frame=['a', '+t{}'.format(utils.fn_basename2(self.src_dem))],
        # )

        fig.grdimage(
            frame=['af', '+t{}'.format(utils.fn_basename2(self.src_dem))],
            grid=self.grid,
            #cmap=self.cpt,
            cmap=True,
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
        
class PerspectoFactory(factory.CUDEMFactory):
    if has_pygmt:
        _modules = {
            'hillshade': {'call': Hillshade},
            'hillshade2': {'call': Hillshade_cmd},
            'perspective': {'call': perspective},
            'sphere': {'call': sphere},
            'figure1': {'call': figure1},
            'colorbar': {'call': colorbar},
        }
    else:
        _modules = {
            'hillshade': {'call': Hillshade},
            'hillshade2': {'call': Hillshade_cmd},
            'perspective': {'call': perspective},
            'sphere': {'call': sphere},
        }    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        
## ==============================================
## Command-line Interface (CLI)
## $ perspecto
##
## perspecto cli
## ==============================================
perspecto_cli_usage = lambda: """{cmd}

usage: {cmd} [ -hvCMZ [ args ] ] DEM ...

Options:

  -C, --cpt\t\tColor Pallette file (if not specified will auto-generate ETOPO CPT)
  -M, --module\t\tDesired perspecto MODULE and options. (see available Modules below)
\t\t\tWhere MODULE is module[:mod_opt=mod_val[:mod_opt1=mod_val1[:...]]]
  -Z, --split-cpt\tSplit the CPT values at zero

  --min_z\t\tMinimum z value to use in CPT
  --max_z\t\tMaximum z value to use in CPT

  --help\t\tPrint the usage text
  --modules\t\tDisplay the module descriptions and usage
  --version\t\tPrint the version information

Supported PERSPECTO modules (see perspecto --modules <module-name> for more info): 
  {d_formats}
""".format(cmd=os.path.basename(sys.argv[0]),
           d_formats=factory._cudem_module_short_desc(PerspectoFactory._modules))
        
#if __name__ == '__main__':
def perspecto_cli(argv = sys.argv):
    i = 1
    src_dem = None
    wg_user = None
    module = None
    src_cpt = None
    min_z = None
    max_z = None
    split_cpt = None
    
    while i < len(argv):
        arg = argv[i]
        if arg == '--module' or arg == '-M':
            module = str(argv[i + 1])
            i += 1
        elif arg[:2] == '-M':
            module = str(arg[2:])

        elif arg == '--cpt' or arg == '-C':
            src_cpt = str(argv[i + 1])
            i += 1
        elif arg[:2] == '-C':
            src_cpt = str(arg[2:])

        elif arg == '--min_z':
            min_z = utils.float_or(argv[i + 1])
            i += 1
            
        elif arg == '--max_z':
            max_z = utils.float_or(argv[i + 1])
            i += 1
        elif arg == '--split-cpt' or arg == '-Z':
            split_cpt = 0
        
        elif arg == '--modules' or arg == '-m':
            factory.echo_modules(PerspectoFactory._modules, None if i+1 >= len(argv) else sys.argv[i+1])
            sys.exit(0)            
        elif arg == '--help' or arg == '-h':
            sys.stderr.write(perspecto_cli_usage())
            sys.exit(0)
        elif arg == '--version' or arg == '-v':
            sys.stdout.write('{}\n'.format(cudem.__version__))
            sys.exit(0)
        elif arg[0] == '-':
            sys.stdout.write(perspecto_cli_usage())
            utils.echo_error_msg('{} is not a valid perspecto cli switch'.format(arg))
            sys.exit(0)
        else:
            wg_user = arg
        i += 1

    if module is None:
        sys.stderr.write(perspecto_cli_usage())
        utils.echo_error_msg(
            '''must specify a perspecto -M module.'''
        )
        sys.exit(-1)

    if module.split(':')[0] not in PerspectoFactory()._modules.keys():
        utils.echo_error_msg(
            '''{} is not a valid perspecto module, available modules are: {}'''.format(
                module.split(':')[0], factory._cudem_module_short_desc(PerspectoFactory._modules)
            )
        )
        sys.exit(-1)
        
    ## load the user wg json and run perspecto with that.
    if wg_user is not None:
        if os.path.exists(wg_user):
            try:
                with open(wg_user, 'r') as wgj:
                    wg = json.load(wgj)
                    if wg['src_region'] is not None:
                        wg['src_region'] = regions.Region().from_list(
                            wg['src_region']
                        )

                    this_waffle = waffles.WaffleFactory(**wg).acquire()
                    this_waffle.mask = True
                    this_waffle.clobber = False

                    if not this_waffle.valid_p():
                        this_waffle.generate()

                    src_dem = this_waffle.fn
            except:
                src_dem = wg_user
        else:
            sys.stderr.write(perspecto_cli_usage())
            utils.echo_error_msg(
                'specified waffles config file/DEM does not exist, {}'.format(wg_user)
            )
            sys.exit(-1)
    else:
        sys.stderr.write(perspecto_cli_usage())
        utils.echo_error_msg(
            'you must supply a waffles config file or an existing DEM; see waffles --help for more information.'
        )
        sys.exit(-1)

    this_perspecto = PerspectoFactory(mod=module, src_dem=src_dem, cpt=src_cpt, min_z=min_z, max_z=max_z, split_cpt=split_cpt)
    if this_perspecto is not None:
        this_perspecto_module = this_perspecto._acquire_module()
        if this_perspecto_module is not None:
            out_fig = this_perspecto_module()
            utils.echo_msg('generated {}'.format(out_fig))
            #this_perspecto_module.run()
        else:
            utils.echo_error_msg('could not acquire perspecto module {}'.format(module))

### End
    
