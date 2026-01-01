### joyplot.py
##
## Copyright (c) 2025 - 2026 Regents of the University of Colorado
##
## joyplot.py is part of CUDEM
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
## Generate "Joy Division" style ridgeline plots from a DEM.
##
### Code:

import numpy as np
from osgeo import gdal
from cudem import utils
from cudem import gdalfun
from .perspecto import Perspecto

## Optional Dependency
try:
    import matplotlib.pyplot as plt
    from matplotlib import cm
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

class Joyplot(Perspecto):
    """Generate a Joyplot (Ridgeline Plot) from a DEM.
    
    Creates a series of vertically offset line profiles representing rows 
    of the DEM. 

    Parameters:
    -----------
    step (int) : Decimation factor (row/col step). Higher = fewer lines. Default: 10.
    scale (float) : Vertical exaggeration of the lines. Default: 1.0.
    overlap (float) : Spacing factor between lines (inverse of gap). 
                      Higher = more overlap. Default: 1.0.
    color (str) : Line color (hex or name). Default: 'black'.
    facecolor (str) : Fill color behind lines (masks background). Default: 'white'.
    dpi (int) : Output resolution. Default: 300.
    
    < joyplot:step=10:scale=2:overlap=1.5:color=black:facecolor=white >
    """

    def __init__(self, step=10, scale=.1, overlap=1, line_color='black', face_color='white', dpi=300):
        super().__init__(mod='joyplot', **kwargs)
        self.step = utils_int_or(step, 10)
        self.scale = utils.float_or(scale, 0.1) # Vertical stretch of the signal
        self.overlap = utils.float_or(overlap, 1.0) # Vertical distance between baselines
        self.line_color = color
        self.face_color = facecolor
        self.dpi = utils.int_or(dpi, 300)

        
    def run(self):
        if not HAS_MATPLOTLIB:
            utils.echo_error_msg("Matplotlib is required for the JOYPLOT module.")
            return self
        
        ## Load Data
        ds = gdal.Open(self.src_dem)
        if ds is None:
            return self
            
        band = ds.GetRasterBand(1)
        ndv = band.GetNoDataValue()
        
        ## Read array and decimate immediately to save memory/processing
        ## Joyplots look better when somewhat smoothed/decimated anyway
        arr = band.ReadAsArray()[::self.step, ::self.step]
        
        ## Handle NoData (Mask or Fill)
        if ndv is not None:
            arr[arr == ndv] = np.nan
            
        ## Flip array so we plot from "back" (North) to "front" (South)
        ## This ensures occlusion works correctly (front covers back)
        ## Assuming GDAL reads Top-Down
        arr = np.flipud(arr)
        
        rows, cols = arr.shape
        x = np.arange(cols)
        
        ## Setup Plot
        ## Figure size calculation (heuristic)
        aspect = cols / rows
        fig_width = 10
        fig_height = fig_width / aspect if aspect > 0 else 10
        
        fig, ax = plt.subplots(figsize=(fig_width, fig_height))
        
        ## Background
        if self.face_color == 'none':
            ax.patch.set_alpha(0.0)
        else:
            ax.set_facecolor(self.face_color)

        ## Plot Loop
        ## We iterate and offset each row
        ## y_base = i * offset
        ## y_value = y_base + (elevation * scale)        
        offset_step = (np.nanmax(arr) - np.nanmin(arr)) * 0.1 * (1/self.overlap)
        if offset_step == 0: offset_step = 1
        
        ## Iterate through rows
        for i, row in enumerate(arr):
            ## Calculate baseline y for this row
            y_base = i * self.offset_step
            
            ## Calculate signal y
            ## Replace NaNs with min for plotting continuity or mask
            row_clean = np.nan_to_num(row, nan=np.nanmin(row))
            y_signal = y_base + (row_clean * scale)
            
            ## Fill below the line (Occlusion)
            ## We fill down to the absolute bottom of the plot or a local baseline
            ## to hide the lines behind this one.
            ax.fill_between(
                x, 
                y_signal, 
                y_base - 10000, # Fill arbitrarily deep to ensure coverage
                facecolor=self.face_color, 
                edgecolor='none',
                zorder=i
            )
            
            ## Draw the line
            ax.plot(
                x, 
                y_signal, 
                color=self.line_color, 
                linewidth=0.5,
                zorder=i
            )

        ## Cleanup Aesthetics
        ax.set_axis_off()
        
        ## Save
        if self.verbose:
            utils.echo_msg(f"Saving joyplot to {self.outfile}...")
            
        plt.savefig(
            self.outfile, 
            bbox_inches='tight', 
            pad_inches=0, 
            dpi=self.dpi,
            transparent=(self.face_color == 'none' or self.face_color == 'transparent')
        )
        plt.close(fig)
        
        return self.outfile

### End
