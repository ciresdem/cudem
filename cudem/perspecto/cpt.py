### cpt.py
##
## Copyright (c) 2026 - 2025 Regents of the University of Colorado
##
## cpt.py is part of CUDEM
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
## CPT (Color Palette Table) generation and processing tools.
##
### Code:

import os
import colorsys
import numpy as np
from osgeo import gdal

try:
    import pygmt
except ImportError:
    pygmt = None

from cudem import utils
from cudem.fetches import fetches

## CPT Colors dictionary
CPT_COLORS = {
    'black': [0, 0, 0],
    'white': [255, 255, 255],
}


def scale_el_simple(value, gmin, gmax, tr):
    """Simple scaling of elevation based on predefined ranges."""
    
    if value > 0 and gmax > 0:
        return (gmax * tr) / 8000
    elif value < 0 and gmin < 0:
        return (gmin * tr) / -11000
    elif value == 0:
        return 0
    else:
        print(value)
        return None

    
def scale_el_relative(value, gmin, gmax, tr, trs):
    """Linearly scales 'tr' from the range [min(trs), max(trs)] to [gmin, gmax].
    Lowest input -> gmin
    Highest input -> gmax
    """
    
    input_min = min(trs)
    input_max = max(trs)
    
    input_range = input_max - input_min
    output_range = gmax - gmin
    
    if input_range == 0:
        return gmin

    percentage = (tr - input_min) / input_range    
    return gmin + (percentage * output_range)

    
def scale_el_relative_etopo(value, gmin, gmax, tr, trs):
    """Scaling relative to the max/min of the input ranges (trs)."""

    if value > 0 and gmax > 0:
        return (gmax * tr) / max(trs)
    elif value < 0 and gmin < 0:
        if min(trs) == 0:
            return gmin * tr
        else:
            return (gmin * tr) / min(trs)
    elif value == 0:
        return gmin
    else:
        return None


def scale_el_linear(value, gmin, gmax, tr, trs):
    """Linear scaling calculation."""
    
    p = (tr - min(trs)) / (max(trs) - min(trs))
    v = (1 - p) * (gmin - gmax) + gmax
    return v


def generate_etopo_cpt(gmin, gmax, output_file='tmp.cpt'):
    """Generates a CPT based on ETOPO1 color steps scaled to gmin/gmax."""
    
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
        3605.149424, 3845.492719
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

    with open(output_file, 'w') as cpt:
        for i, j in enumerate(elevs):
            if j is not None and i + 1 < len(elevs):
                elev_curr = scale_el_relative_etopo(j, gmin, gmax, trs[i], trs)
                elev_next = scale_el_relative_etopo(elevs[i + 1], gmin, gmax, trs[i + 1], trs)
                
                if elev_curr is not None and elev_next is not None:
                    c1 = colors[i]
                    cpt.write(
                        f'{elev_curr} {c1[0]} {c1[1]} {c1[2]} '
                        f'{elev_next} {c1[0]} {c1[1]} {c1[2]}\n'
                    )
    return output_file


def make_cpt(cmap='etopo1', color_model='r', min_z=None, max_z=None, output=None):
    """Wrapper around pygmt.makecpt."""
    
    if pygmt is None:
        utils.echo_error_msg("PyGMT is not installed.")
        return None

    pygmt.makecpt(
        cmap=cmap,
        color_model=color_model,
        output=output,
        series=[min_z, max_z, 50],
        no_bg=True,
    )
    return output


def process_cpt(cpt_file, gmin, gmax, gdal=False, split_cpt=None):
    """Parses an existing CPT and rescales it to gmin/gmax."""
    
    trs = []
    colors = []

    if cpt_file is None:
        return None
    
    with open(cpt_file, 'r') as f:
        for line in f:
            parts = line.split()
            if not parts:
                continue
            
            # Check if first part is a number (elevation)
            if utils.float_or(parts[0]) is not None:
                trs.append(float(parts[0]))

                # Parse color
                if utils.int_or(parts[1]) is not None:
                    # R G B format
                    colors.append([
                        int(float(parts[1])),
                        int(float(parts[2])),
                        int(float(parts[3]))
                    ])
                elif parts[1] in CPT_COLORS:
                    # Named color
                    colors.append(CPT_COLORS[parts[1]])
                elif '/' in parts[1]:
                    # Slash delimited
                    colors.append([int(float(x)) for x in parts[1].split('/')])

    # Determine elevation steps
    if split_cpt is not None:
        _trs = np.array(trs)
        len_b = len(_trs[_trs < split_cpt])
        len_t = len(_trs[_trs > split_cpt])
        
        elevs_b = np.linspace(gmin, split_cpt, len_b)
        # Skip first element of top to avoid duplicate/overlap if necessary
        elevs_t = np.linspace(split_cpt, gmax, len_t)[1:]
        elevs = np.concatenate((elevs_b, elevs_t))
    else:
        elevs = np.linspace(gmin, gmax, len(trs))

    output_fn = 'tmp.cpt'
    with open(output_fn, 'w') as f_out:
        for i, j in enumerate(elevs):
            if j is not None and i + 1 < len(elevs):
                # Choose scaling method based on range crossing zero
                if gmin < 0 and gmax > 0:
                    elev_curr = scale_el_relative(j, gmin, gmax, trs[i], trs)
                    elev_next = scale_el_relative(elevs[i + 1], gmin, gmax, trs[i + 1], trs)
                else:
                    elev_curr = scale_el_linear(j, gmin, gmax, trs[i], trs)
                    elev_next = scale_el_linear(elevs[i + 1], gmin, gmax, trs[i + 1], trs)
                    
                if elev_curr is not None and elev_next is not None:
                    c = colors[i]
                    if not gdal:
                        f_out.write(
                            f'{elev_curr} {c[0]} {c[1]} {c[2]} '
                            f'{elev_next} {c[0]} {c[1]} {c[2]}\n'
                        )
                    else:
                        # GDAL format typically requires index or specific structuring
                        f_out.write(f'{elev_curr} {c[0]} {c[1]} {c[2]} 255\n')
        
        if gdal:
            f_out.write('nv 0 0 0 0\n')

    return output_fn


def fetch_cpt_city(q='grass/haxby', cache_dir='./'):
    """Fetches a CPT file from the cpt-city archive."""
    
    utils.echo_msg(f'checking for `{q}` at cpt-city')
    
    this_fetches = fetches.FetchesFactory(
        mod='cpt_city', verbose=True, outdir=cache_dir, q=q
    )._acquire_module()
    
    this_fetches.run()
    
    if not this_fetches.results:
        utils.echo_error_msg(f"No results found for {q}")
        return None

    result = this_fetches.results[0]
    utils.echo_msg(
        f"found {len(this_fetches.results)} results for `{q}` at cpt-city, "
        f"using first entry `{result['url']}`"
    )

    result['dst_fn'] = os.path.join(this_fetches._outdir, result['dst_fn'])
    fetches.Fetch(result['url']).fetch_file(result['dst_fn'])
    
    return result['dst_fn']


def get_correct_map(path, luminosity, contrast):
    """Adjusts the luminosity and contrast of a GDAL raster.

    Note: This iterates pixel-by-pixel and is slow for large images.
    """
    
    if not os.path.exists(path):
        utils.echo_error_msg(f"File not found: {path}")
        return

    ds = gdal.Open(path)
    if ds is None:
        utils.echo_error_msg(f"Could not open {path}")
        return

    band1_obj = ds.GetRasterBand(1)
    
    # Determine Max Value based on DataType
    if band1_obj.DataType == gdal.GDT_UInt16:
        max_value = int(2**16 - 1)
    elif band1_obj.DataType == gdal.GDT_Byte:
        max_value = int(2**8 - 1)
    else:
        utils.echo_msg(
            f'band type {band1_obj.DataType} not handled: using default 16-bit'
        )
        max_value = int(2**16 - 1)

    # Read bands (Assuming RGB)
    # ReadAsArray returns [bands, y, x]
    data_array = ds.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
    if len(data_array) < 3:
        utils.echo_error_msg("Image requires at least 3 bands (RGB).")
        return

    band1 = data_array[0]
    band2 = data_array[1]
    band3 = data_array[2]

    # Process pixel by pixel (Slow)
    # TODO: Vectorize this loop using numpy for performance
    for x in range(ds.RasterXSize):
        for y in range(ds.RasterYSize): # Should check if dims match array shape (usually y,x)

            r = float(band1[y, x]) / max_value
            g = float(band2[y, x]) / max_value
            b = float(band3[y, x]) / max_value

            # Convert to HLS then apply luminosity and contrast
            h, l, s = colorsys.rgb_to_hls(r, g, b)

            l = min(max(0, l + (l - 0.5) * (luminosity - 0.5)), 1)
            s = min(max(0, s + (s - 0.5) * (contrast - 0.5)), 1)

            r_out, g_out, b_out = colorsys.hls_to_rgb(h, l, s)
            
### End
