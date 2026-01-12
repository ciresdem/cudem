### cpt.py
##
## Copyright (c) 2023 - 2026 Regents of the University of Colorado
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
    'red': [255, 0, 0],
    'green': [0, 255, 0],
    'blue': [0, 0, 255],
    'yellow': [255, 255, 0],
    'cyan': [0, 255, 255],
    'magenta': [255, 0, 255],
    'gray': [128, 128, 128],
    'lightgray': [211, 211, 211],
    'darkgray': [169, 169, 169],
    'orange': [255, 165, 0],
    'purple': [128, 0, 128],
    'brown': [165, 42, 42],
}


## The following scale_el_* functions are depreciated.
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
    """Generates a CPT based on ETOPO1 color steps scaled to gmin/gmax.
    
    The ETOPO1 color map assumes a split at 0 (Sea Level).
    This function scales the negative values to [gmin, 0] and 
    positive values to [0, gmax].
    """
    
    ## Original ETOPO1 Thresholds
    trs = [
        -11000, -10500, -10000, -9500, -9000, -8500,
        -8000, -7500, -7000, -6500, -6000, -5500,
        -5000, -4500, -4000, -3500, -3000, -2500,
        -2000, -1500, -1000, -500, -0.001, 0,
        100, 200, 500, 1000, 1500, 2000, 2500,
        3000, 3500, 4000, 4500, 5000, 5500,
        6000, 6500, 7000, 7500, 8000
    ]
    
    ## Original ETOPO1 Colors
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

    ## Use the process_cpt logic to calculate new elevations
    ## by treating the hardcoded trs as the input.
    ## ETOPO is hinged at 0.
    new_elevs = []
    split_val = 0
    t_min, t_max = min(trs), max(trs)

    for t in trs:
        if t <= split_val:
            ## Scale [t_min, 0] -> [gmin, 0]
            if t_min == split_val: pct = 0
            else: pct = (t - t_min) / (split_val - t_min)
            val = gmin + pct * (0 - gmin)
        else:
            ## Scale [0, t_max] -> [0, gmax]
            if t_max == split_val: pct = 0
            else: pct = (t - split_val) / (t_max - split_val)
            val = 0 + pct * (gmax - 0)
        new_elevs.append(val)

    with open(output_file, 'w') as cpt:
        for i in range(len(new_elevs) - 1):
            elev_curr = new_elevs[i]
            elev_next = new_elevs[i+1]
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
    """Parses an existing CPT and rescales it to [gmin, gmax].
    
    If split_cpt is provided (e.g., 0), the scaling is piece-wise:
    - Input values <= split_cpt are scaled to [gmin, split_cpt]
    - Input values >= split_cpt are scaled to [split_cpt, gmax]
    This preserves the hinge/break point in the color map (e.g. Sea Level).
    """
    
    if cpt_file is None:
        return None

    trs = []
    colors = []
    
    ## Parse the Input CPT
    with open(cpt_file, 'r') as f:
        for line in f:
            parts = line.split()
            if not parts:
                continue
            
            ## Use 'float_or' to safely check if line starts with a number
            if utils.float_or(parts[0]) is not None:
                trs.append(float(parts[0]))
                
                ## Parse RGB
                ## Assumes format: Z R G B ...
                if utils.int_or(parts[1]) is not None:
                    colors.append([
                        int(float(parts[1])),
                        int(float(parts[2])),
                        int(float(parts[3]))
                    ])
                elif parts[1] in CPT_COLORS:
                    colors.append(CPT_COLORS[parts[1]])
                elif '/' in parts[1]:
                    colors.append([int(float(x)) for x in parts[1].split('/')])

    if not trs:
        utils.echo_error_msg(f"No valid data found in CPT {cpt_file}")
        return None

    ## Calculate New Z-Values (Elevations)
    new_elevs = []
    t_min, t_max = min(trs), max(trs)

    if split_cpt is not None:
        split_val = float(split_cpt)
        ## Check if split_val is within input range to avoid extrapolation errors
        ## Though extrapolation might be desired if gmin/gmax force it        
        for t in trs:
            if t <= split_val:
                ## Bottom Segment: [t_min, split_val] -> [gmin, split_val]
                if split_val == t_min:
                    val = gmin 
                else:
                    ## Percentage of distance through input bottom segment
                    pct = (t - t_min) / (split_val - t_min)
                    ## Map to output bottom segment
                    val = gmin + pct * (split_val - gmin)
            else:
                ## Top Segment: [split_val, t_max] -> [split_val, gmax]
                if t_max == split_val:
                    val = gmax
                else:
                    ## Percentage of distance through input top segment
                    pct = (t - split_val) / (t_max - split_val)
                    ## Map to output top segment
                    val = split_val + pct * (gmax - split_val)
            
            new_elevs.append(val)
            
    else:
        ## Standard Linear Stretch: [t_min, t_max] -> [gmin, gmax]
        for t in trs:
            if t_max == t_min:
                val = gmin
            else:
                pct = (t - t_min) / (t_max - t_min)
                val = gmin + pct * (gmax - gmin)
            new_elevs.append(val)

    ## Write the Output CPT
    output_fn = 'tmp.cpt'
    with open(output_fn, 'w') as f_out:
        ## We iterate up to len-1 because a segment connects i to i+1
        for i in range(len(new_elevs) - 1):
            elev_curr = new_elevs[i]
            elev_next = new_elevs[i+1]
            c = colors[i]  # Use the color associated with the start of the segment
            
            if not gdal:
                ## Standard CPT Format
                f_out.write(
                    f'{elev_curr} {c[0]} {c[1]} {c[2]} '
                    f'{elev_next} {c[0]} {c[1]} {c[2]}\n'
                )
            else:
                ## GDAL Color Relief Format (one value per line)
                f_out.write(f'{elev_curr} {c[0]} {c[1]} {c[2]} 255\n')
        
        ## Ensure the last value is written for GDAL format
        if gdal and len(new_elevs) > 0:
            last_c = colors[-1]
            f_out.write(f'{new_elevs[-1]} {last_c[0]} {last_c[1]} {last_c[2]} 255\n')
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
    """Adjusts the luminosity and contrast of a GDAL raster using Vectorized NumPy.
    Saves the result back to the file.
    """
    
    if not os.path.exists(path):
        utils.echo_error_msg(f"File not found: {path}")
        return

    ## Open in Update mode to allow saving changes
    ds = gdal.Open(path, gdal.GA_Update)
    if ds is None:
        utils.echo_error_msg(f"Could not open {path}")
        return

    ## Check Band Count
    if ds.RasterCount < 3:
        utils.echo_error_msg("Image requires at least 3 bands (RGB).")
        return

    ## Check Data Type and Determine Max Value
    band_obj = ds.GetRasterBand(1)
    dtype = band_obj.DataType
    
    if dtype == gdal.GDT_UInt16:
        max_value = 65535
        np_dtype = np.uint16
    elif dtype == gdal.GDT_Byte:
        max_value = 255
        np_dtype = np.uint8
    else:
        # Default fallback
        max_value = 65535
        np_dtype = np.uint16

    ## Read All Bands into Memory (c, y, x)
    ## This reads R, G, B into indices 0, 1, 2
    try:
        data = ds.ReadAsArray()
    except Exception as e:
        utils.echo_error_msg(f"Failed to read raster array: {e}")
        return

    ## Normalize to Float [0, 1]
    ## We use a copy or cast to float for calculation
    img_norm = data[:3].astype(np.float32) / max_value
    
    r = img_norm[0]
    g = img_norm[1]
    b = img_norm[2]
    
    ## Convert to HLS (Vectorized)
    h, l, s = _rgb_to_hls_vectorized(r, g, b)
    
    ## Apply Luminosity and Contrast
    ## Formula: L_new = L + (L - 0.5) * (Factor - 0.5)
    ## Clamping is handled by clip at the end or logic
    if luminosity != 0.5:
        l = l + (l - 0.5) * (luminosity - 0.5)
        l = np.clip(l, 0, 1)
        
    if contrast != 0.5:
        s = s + (s - 0.5) * (contrast - 0.5)
        s = np.clip(s, 0, 1)

    ## Convert back to RGB (Vectorized)
    r_new, g_new, b_new = _hls_to_rgb_vectorized(h, l, s)
    
    ## Scale back to Integer Range
    r_out = (r_new * max_value).astype(np_dtype)
    g_out = (g_new * max_value).astype(np_dtype)
    b_out = (b_new * max_value).astype(np_dtype)
    
    ## Write Data Back to Bands
    ds.GetRasterBand(1).WriteArray(r_out)
    ds.GetRasterBand(2).WriteArray(g_out)
    ds.GetRasterBand(3).WriteArray(b_out)
    
    ## Force Flush/Close
    ds.FlushCache()
    ds = None
    
    utils.echo_msg(f"Updated color map for {path}")    

### End
