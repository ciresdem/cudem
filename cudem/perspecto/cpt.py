### cpt.py 
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

from cudem.perspecto import perspecto
from cudem import utils
import numpy as np
import pygmt


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


def makecpt(cmap='etopo1', color_model='r', min_z = None, max_z = None, output=None):
    pygmt.makecpt(
        cmap=cmap,
        color_model=color_model,
        output=output,
        series=[min_z, max_z, 50],
        no_bg=True,
        #continuous=True,
        #truncate=[min_z, max_z],
    )
    
    return(output)

def process_cpt(cpt, gmin, gmax, gdal=False, split_cpt=None):
    trs = []
    colors = []
    with open(cpt, 'r') as cpt:
        for l in cpt:
            ll = l.split()
            if len(ll) == 0:
                continue
            
            if utils.float_or(ll[0]) is not None:
                trs.append(float(ll[0]))

                if utils.int_or(ll[1]) is not None:
                    colors.append(
                        [int(float(ll[1])),
                         int(float(ll[2])),
                         int(float(ll[3]))]
                    )
                elif ll[1] in cpt_colors.keys():
                    colors.append(cpt_colors[ll[1]])
                elif len(ll[1].split('/')) > 1:
                    colors.append(
                        [int(float(x)) for x in ll[1].split('/')]
                    )

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
                if gmin < 0 and gmax > 0:
                    elev = scale_el_(j, gmin, gmax, trs[i], trs)
                    elev1 = scale_el_(elevs[i + 1], gmin, gmax, trs[i + 1], trs)
                else:
                    elev = scale_el(j, gmin, gmax, trs[i], trs)
                    elev1 = scale_el(elevs[i + 1], gmin, gmax, trs[i + 1], trs)
                    
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


def fetch_cpt_city(q='grass/haxby', cache_dir='./'):
    utils.echo_msg(f'checking for `{q}` at cpt-city')
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
            LOGGER.info(
                (f'band type {band1.DataType} not handled: '
                'use default size of value (16 bits)')
            )

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

                    
### End
