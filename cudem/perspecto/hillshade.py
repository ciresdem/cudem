### hillshade.py 
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

import os

from osgeo import gdal
from osgeo_utils import gdal_calc

from cudem.perspecto import perspecto
from cudem import utils
from cudem import gdalfun
import numpy as np

try:
    import pygmt
    has_pygmt = True
except ImportError as e:
    has_pygmt = False
    

class Hillshade(perspecto.Perspecto):
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
            'mogrify -modulate {} -depth 8 {}'.format(
                self.modulate, gdal_fn
            ),
            verbose=False
        )

        
    def gamma_correction_calc(self, hs_fn=None, gamma=.5,
                              outfile='gdaldem_gamma.tif'):
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

        
    def blend(self, hs_file, rgb_file, mode='multiply', gamma=None,
              outfile='gdaldem_multiply.tif'):
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

    
    def gmt_figure(self, colorbar_text='Elevation'):
        fig = pygmt.Figure()
        # fig.basemap(
        #     region=self.dem_region.format('str'),
        #     frame=[f'a', '+t{self.outfile}'],
        # )
        
        fig.grdimage(
            frame=[f'af', '+t{self.outfile}'],
            grid=self.outfile,
            cmap=False,
        )

        #if self.want_colorbar:
        #fig.colorbar(frame=['a{}'.format(self.interval), 'x+lElevation', 'y+1m'])
        fig.colorbar(frame=['x+l{}'.format(colorbar_text), 'y+1m'])

        fig.savefig('{}_gmt.png'.format(utils.fn_basename2(self.src_dem)))
            
        
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

        #self.gmt_figure()
        return(self.outfile)

    
class Hillshade_cmd(perspecto.Perspecto):
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

    
class HillShade_test(perspecto.Perspecto):
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


### End
