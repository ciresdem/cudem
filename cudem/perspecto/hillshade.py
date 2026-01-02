### hillshade.py
##
## Copyright (c) 2015 - 2026 Regents of the University of Colorado
##
## hillshade.py is part of CUDEM
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
## Generate Hillshades from a DEM using GDAL and Numpy blending.
##
### Code:

import os
import numpy as np
from osgeo import gdal
from osgeo_utils import gdal_calc

from cudem import utils
from cudem import gdalfun
from cudem.perspecto import perspecto

try:
    import pygmt
    HAS_PYGMT = True
except ImportError:
    HAS_PYGMT = False


class Hillshade(perspecto.Perspecto):
    """Generate a Hillshade Image using GDAL and Numpy blending.
    
    Supports various blending modes:
    https://en.wikipedia.org/wiki/Blend_modes#Overlay

    Configuration Example:
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
        gamma=0.5,
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
        """Adjusts brightness/saturation using ImageMagick (External utility)."""
        
        utils.run_cmd(
            f'mogrify -modulate {self.modulate} -depth 8 {gdal_fn}',
            verbose=False
        )

        
    def gamma_correction_calc(self, hs_fn=None, outfile='gdaldem_gamma.tif'):
        """Apply gamma correction using gdal_calc (External utility)."""
        
        gdal_calc.Calc(
            "uint8(((A / 255.)**(1/0.5)) * 255)",
            A=hs_fn,
            outfile=outfile,
            quiet=True
        )
        return outfile

    
    def gamma_correction(self, arr):
        """Apply gamma correction to a numpy array."""
        
        if self.gamma is not None:
            return (((arr / 255.)**(1 / self.gamma)) * 255.).astype(np.uint8)
        return arr

    
    def blend(self, hs_file, rgb_file, mode='multiply', gamma=None, outfile='gdaldem_multiply.tif'):
        """Blend a Hillshade (hs_file) and Color Relief (rgb_file) using the specified mode."""
        
        valid_modes = ['multiply', 'screen', 'overlay', 'hard_light', 'soft_light']
        if mode not in valid_modes:
            mode = 'multiply'

        with gdalfun.gdal_datasource(rgb_file, update=True) as rgb_ds:
            with gdalfun.gdal_datasource(hs_file) as cr_ds:
                ds_config = gdalfun.gdal_infos(cr_ds)
                n_chunk = int(ds_config['nx'] * 0.1)
                
                # Iterate over the image in chunks to manage memory
                for srcwin in utils.yield_srcwin(
                    (ds_config['ny'], ds_config['nx']), 
                    n_chunk=n_chunk,
                    msg=f'Blending rgb and hillshade using {mode}',
                    end_msg='Generated color hillshade',
                    verbose=self.verbose
                ):
                    cr_band = cr_ds.GetRasterBand(1)
                    cr_arr = cr_band.ReadAsArray(*srcwin)
                    cr_arr = self.gamma_correction(cr_arr)
                    
                    # Normalize hillshade for calculations
                    cr_norm = cr_arr / 255.0

                    for band_no in [1, 2, 3]:
                        band = rgb_ds.GetRasterBand(band_no)
                        band_arr = band.ReadAsArray(*srcwin)
                        band_norm = band_arr / 255.0
                        
                        out_arr = np.copy(band_arr)

                        if mode == 'multiply':
                            out_arr = (cr_norm * band_norm * 255.).astype(np.uint8)

                        elif mode == 'screen':
                            out_arr = ((1 - ((1 - cr_norm) * (1 - band_norm))) * 255.).astype(np.uint8)

                        elif mode == 'overlay':
                            # Contrast mode: Multiplies or Screens based on bottom layer
                            mask_low = cr_arr < 128
                            mask_high = cr_arr >= 128
                            
                            out_arr[mask_low] = (2 * cr_norm[mask_low] * band_norm[mask_low] * 255.).astype(np.uint8)
                            out_arr[mask_high] = ((1 - (2 * (1 - cr_norm[mask_high]) * (1 - band_norm[mask_high]))) * 255.).astype(np.uint8)

                        elif mode == 'hard_light':
                            # Like Overlay, but based on top layer
                            mask_low = band_arr < 128
                            mask_high = band_arr >= 128
                            
                            out_arr[mask_low] = (2 * cr_norm[mask_low] * band_norm[mask_low] * 255.).astype(np.uint8)
                            out_arr[mask_high] = ((1 - (2 * (1 - cr_norm[mask_high]) * (1 - band_norm[mask_high]))) * 255.).astype(np.uint8)

                        elif mode == 'soft_light':
                            mask_low = band_arr < 128
                            mask_high = band_arr >= 128
                            
                            # Soft light formula for darker pixels
                            term1_low = 2 * cr_norm[mask_low] * band_norm[mask_low]
                            term2_low = (cr_norm[mask_low]**2) * (1 - 2 * band_norm[mask_low])
                            out_arr[mask_low] = ((term1_low + term2_low) * 255.).astype(np.uint8)
                            
                            # Soft light formula for lighter pixels
                            term1_high = 2 * cr_norm[mask_high] * (1 - band_norm[mask_high])
                            term2_high = np.sqrt(cr_norm[mask_high]) * (2 * band_norm[mask_high] - 1)
                            out_arr[mask_high] = ((term1_high + term2_high) * 255.).astype(np.uint8)

                        band.WriteArray(out_arr, srcwin[0], srcwin[1])

        os.rename(rgb_file, outfile)
        return outfile

    
    def gmt_figure(self, colorbar_text='Elevation'):
        """Generate a figure using PyGMT (if available)."""
        
        if not HAS_PYGMT:
            utils.echo_msg("PyGMT is not installed.")
            return

        fig = pygmt.Figure()
        fig.grdimage(
            frame=[f'af', f'+t{self.outfile}'],
            grid=self.outfile,
            cmap=False,
        )
        fig.colorbar(frame=[f'x+l{colorbar_text}', 'y+1m'])
        fig.savefig(f'{utils.fn_basename2(self.src_dem)}_gmt.png')

        
    def run(self):
        hs_fn = utils.make_temp_fn('gdaldem_hs.tif', self.outdir)
        gdal.DEMProcessing(
            hs_fn, 
            self.src_dem, 
            'hillshade', 
            computeEdges=True, 
            scale=111120,
            azimuth=self.azimuth, 
            altitude=self.altitude, 
            zFactor=self.vertical_exaggeration
        )

        cr_fn = utils.make_temp_fn('gdaldem_cr.tif', self.outdir)
        gdal.DEMProcessing(
            cr_fn, 
            self.src_dem, 
            'color-relief', 
            colorFilename=self.cpt,
            computeEdges=True, 
            addAlpha=self.alpha
        )

        cr_hs_fn = utils.make_temp_fn('gdaldem_cr_hs.tif', self.outdir)
        self.blend(hs_fn, cr_fn, gamma=self.gamma, mode=self.mode, outfile=cr_hs_fn)

        utils.remove_glob(hs_fn, cr_fn)
        os.rename(cr_hs_fn, self.outfile)
        
        # self._modulate(self.outfile) # Optional ImageMagick post-processing
        # self.gmt_figure()            # Optional GMT figure generation
        
        return self.outfile


class Hillshade_cmd(perspecto.Perspecto):
    """
    Generate a Hillshade Image using command-line tools.
    Requires GDAL and ImageMagick installed in the system path.
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
        basename = utils.fn_basename2(self.src_dem)
        
        # 1. Generate Hillshade
        utils.run_cmd(
            f'gdaldem hillshade -compute_edges -s 111120 '
            f'-z {self.vertical_exaggeration} -az {self.azimuth} -alt {self.altitude} '
            f'{self.src_dem} hillshade.tif',
            verbose=self.verbose
        )

        # 2. Generate Color Relief
        utils.run_cmd(
            f'gdaldem color-relief {self.src_dem} {self.cpt} colors.tif',
            verbose=self.verbose
        )

        # 3. Composite (ImageMagick)
        utils.run_cmd(
            'composite -compose multiply -depth 8 colors.tif hillshade.tif output.tif',
            verbose=self.verbose
        )
        utils.run_cmd(
            'mogrify -modulate 115 -depth 8 output.tif',
            verbose=self.verbose
        )

        # 4. Restore Georeferencing
        # Create a TFW file from source
        utils.run_cmd(
            f'gdal_translate -co "TFW=YES" {self.src_dem} temp.tif',
            verbose=self.verbose
        )
        utils.run_cmd('mv temp.tfw output.tfw')
        
        # Translate back to GeoTiff
        srs_cmd = f'-a_srs epsg:{self.projection}' if self.projection else ''
        utils.run_cmd(
            f'gdal_translate {srs_cmd} output.tif temp2.tif'
        )

        # 5. Cleanup
        utils.remove_glob(
            'output.tif*', 'temp.tif*', 'hillshade.tif*', 'colors.tif*', 'output.tfw*'
        )
        
        # Final move
        final_output = f'{basename}_hs.tif'
        utils.run_cmd(f'gdal_translate temp2.tif {final_output}')
        utils.remove_glob('temp2.tif')

        return final_output

### End
