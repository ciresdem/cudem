### histogram.py
##
## Copyright (c) 2025 - 2026 Regents of the University of Colorado
##
## histogram.py is part of CUDEM
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
## Generate Hypsometric Histograms and CDF plots from a DEM.
##
### Code:

import os
import numpy as np
from osgeo import gdal
from cudem import utils
from .perspecto import Perspecto

## Optional Dependency
try:
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

class Histogram(Perspecto):
    """Generate a Hypsometric Histogram or CDF from a DEM.
    
    Visualizes the distribution of elevation values within the dataset.

    Parameters:
    -----------
    bins (int) : Number of histogram bins. Default: 100.
    type (str) : Plot type ['pdf', 'cdf']. Default: 'pdf'.
                 - pdf: Probability Density Function (Standard Histogram)
                 - cdf: Cumulative Distribution Function
    color (str) : Bar color. Default: 'gray'.
    stats (bool) : Overlay mean/median/std_dev lines. Default: False.
    title (str) : Custom title for the plot.
    dpi (int) : Output resolution. Default: 300.
    
    < histogram:bins=100:type=pdf:stats=True:color=steelblue >
    """

    def __init__(self, bins=100, type='pdf', color='gray', stats=False, title=None, dpi=300, **kwargs):
        super().__init__(mod='histogram', **kwargs)
        self.bins = utils.int_or(bins, 100)
        self.type = type
        self.color = color
        self.stats = stats
        self.title = title
        self.dpi = utils.int_or(dpi)

        
    def run(self):
        if not HAS_MATPLOTLIB:
            utils.echo_error_msg("Matplotlib is required for the HISTOGRAM module.")
            return self

        ## Parse Parameters
        bins = self.bins
        plot_type = self.type.lower()
        color = self.color
        show_stats = utils.str2bool(self.stats)
        title = self.title
        dpi = int(self.dpi)
        
        ## Load Data
        ds = gdal.Open(self.src_dem)
        if ds is None:
            return self
            
        band = ds.GetRasterBand(1)
        ndv = band.GetNoDataValue()
        
        ## Read array
        if self.verbose:
            utils.echo_msg("Reading raster data...")
        
        ## For very large files, this might need chunking, but for stats 
        ## usually decimation or full read is acceptable.
        arr = band.ReadAsArray()
        
        ## Flatten and Filter NoData
        if ndv is not None:
            data = arr[arr != ndv]
            ## Also filter NaNs just in case
            data = data[~np.isnan(data)]
        else:
            data = arr.flatten()
            data = data[~np.isnan(data)]
            
        if len(data) == 0:
            utils.echo_error_msg("No valid data found in DEM.")
            return self

        ## Calculate Statistics
        mean_val = np.mean(data)
        median_val = np.median(data)
        std_val = np.std(data)
        min_val = np.min(data)
        max_val = np.max(data)

        ## Plotting
        fig, ax = plt.subplots(figsize=(10, 6))
        
        cumulative = (plot_type == 'cdf')
        density = True # Normalize y-axis
        
        ## Histogram
        n, bins_out, patches = ax.hist(
            data, 
            bins=bins, 
            density=density, 
            cumulative=cumulative, 
            color=color, 
            alpha=0.75, 
            edgecolor='black',
            linewidth=0.5
        )

        ## Labels
        if title is None:
            title = f"Hypsometry: {utils.fn_basename2(self.src_dem)}"
            
        ax.set_title(title, fontsize=14, fontweight='bold')
        ax.set_xlabel("Elevation (z)", fontsize=12)
        ax.set_ylabel("Frequency / Probability", fontsize=12)
        ax.grid(True, linestyle='--', alpha=0.5)

        ## Overlay Stats
        if show_stats:
            stats_text = (
                f"Min: {min_val:.2f}\n"
                f"Max: {max_val:.2f}\n"
                f"Mean: {mean_val:.2f}\n"
                f"Median: {median_val:.2f}\n"
                f"Std Dev: {std_val:.2f}"
            )
            
            ## Vertical Lines
            ax.axvline(mean_val, color='red', linestyle='dashed', linewidth=1.5, label=f'Mean ({mean_val:.1f})')
            ax.axvline(median_val, color='green', linestyle='dashed', linewidth=1.5, label=f'Median ({median_val:.1f})')
            
            if not cumulative:
                ## Add Sigma lines for PDF
                ax.axvline(mean_val + std_val, color='orange', linestyle='dotted', label='1 Std Dev')
                ax.axvline(mean_val - std_val, color='orange', linestyle='dotted')

            ax.legend(loc='upper right')
            
            ## Add text box
            props = dict(boxstyle='round', facecolor='white', alpha=0.8)
            ax.text(0.02, 0.95, stats_text, transform=ax.transAxes, fontsize=10,
                    verticalalignment='top', bbox=props)

        ## Save
        if self.verbose:
            utils.echo_msg(f"Saving histogram to {self.outfile}...")
            
        plt.savefig(self.outfile, dpi=dpi, bbox_inches='tight')
        plt.close(fig)
        
        return self.outfile


### End
