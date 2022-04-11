### vdatums.py
##
## Copyright (c) 2021, 2022 Regents of the University of Colorado
##
## vdatums.py is part of CUDEM
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
### Code:


import os
import sys
import json
import numpy as np

from osgeo import gdal
from osgeo import ogr

from cudem import regions
from cudem import demfun
from cudem import waffles
from cudem import utils
from cudem import htdpfun

import cudem.fetches.utils
import cudem.fetches.tides
import cudem.fetches.vdatum

_tidal_frames = {
    1089: {'name': 'mllw',
           'description': 'Mean Lower Low Water'},
    1091: {'name': 'mlw',
           'description': 'Mean Low Water'},
    1090: {'name': 'mhhw',
           'description': 'Mean Higher High Water'},
    5869: {'name': 'mhhw',
           'description': 'Mean Higher High Water'},
    5868: {'name': 'mhw',
           'description': 'Mean High Water'},
    5174: {'name': 'msl',
           'description': 'Mean Sea Level'},
}

_htdp_reference_frames = {
    4269: {'name': 'NAD_83(2011/CORS96/2007)',
           'description': '(North American plate fixed)',
           'htdp_id': 1},
    6781: {'name': 'NAD_83(2011/CORS96/2007)',
           'description': '(North American plate fixed)',
           'htdp_id': 1},
    6319: {'name': 'NAD_83(2011/CORS96/2007)',
           'description': '(North American plate fixed)',
           'htdp_id': 1},
    6321: {'name': 'NAD_83(PA11/PACP00)',
           'description': '(Pacific plate fixed)',
           'htdp_id': 2},
    6324: {'name': 'NAD_83(MA11/MARP00)',
           'description': '(Mariana plate fixed)',
           'htdp_id': 3},
    4979: {'name': 'WGS_84(original)',
           'description': '(NAD_83(2011) used)',
           'htdp_id': 4},
    7815: {'name': 'WGS_84(original)',
           'description': '(NAD_83(2011) used)',
           'htdp_id': 4},
    7816: {'name': 'WGS_84(original)',
           'description': '(NAD_83(2011) used)',
           'htdp_id': 4},
    7656: {'name': 'WGS_84(G730)',
           'description': '(ITRF91 used)',
           'htdp_id': 5},
    7657: {'name': 'WGS_84(G730)',
           'description': '(ITRF91 used)',
           'htdp_id': 5},
    7658: {'name': 'WGS_84(G873)',
           'description': '(ITRF94 used)',
           'htdp_id': 6},
    7659: {'name': 'WGS_84(G873)',
           'description': '(ITRF94 used)',
           'htdp_id': 6},
    7660: {'name': 'WGS_84(G1150)',
           'description': '(ITRF2000 used)',
           'htdp_id': 7},
    7661: {'name': 'WGS_84(G1150)',
           'description': '(ITRF2000 used)',
           'htdp_id': 7},
    7662: {'name': 'WGS_84(G1674)',
           'description': '(ITRF2008 used)',
           'htdp_id': 8},
    7663: {'name': 'WGS_84(G1674)',
           'description': '(ITRF2008 used)',
           'htdp_id': 8},
    7664: {'name': 'WGS_84(G1762)',
           'description': '(IGb08 used)',
           'htdp_id': 9},
    7665: {'name': 'WGS_84(G1762)',
           'description': '(IGb08 used)',
           'htdp_id': 9},
    7666: {'name': 'WGS_84(G2139)',
           'description': '(ITRF2014=IGS14=IGb14 used)',
           'htdp_id': 10},
    7667: {'name': 'WGS_84(G2139)',
           'description': '(ITRF2014=IGS14=IGb14 used)',
           'htdp_id': 10},
    4910: {'name': 'ITRF88',
           'description': '',
           'htdp_id': 11},
    4911: {'name': 'ITRF89',
           'description': '',
           'htdp_id': 12},
    7901: {'name': 'ITRF89',
           'description': '',
           'htdp_id': 12},
    7902: {'name': 'ITRF90',
           'description': '(PNEOS90/NEOS90)',
           'htdp_id': 13},
    7903: {'name': 'ITRF91',
           'description': '',
           'htdp_id': 14},
    7904: {'name': 'ITRF92',
           'description': '',
           'htdp_id': 15},
    7905: {'name': 'ITRF93',
           'description': '',
           'htdp_id': 16},
    7906: {'name': 'ITRF94',
           'description': '',
           'htdp_id': 17},
    7907: {'name': 'ITRF96',
           'description': '',
           'htdp_id': 18},
    7908: {'name': 'ITRF97',
           'description': 'IGS97',
           'htdp_id': 19},
    7909: {'name': 'ITRF2000',
           'description': 'IGS00/IGb00',
           'htdp_id': 20},
    7910: {'name': 'ITRF2005',
           'description': 'IGS05',
           'htdp_id': 21},
    7911: {'name': 'ITRF2008',
           'description': 'IGS08/IGb08',
           'htdp_id': 22},
    7912: {'name': 'ITRF2014',
           'description': 'IGS14/IGb14',
           'htdp_id': 23},
}

_cdn_reference_frames = {
    9245: {'name': 'CGVD2013(CGG2013a) height',},
    6647: {'name': 'CGVD2013(CGG2013) height',},
    3855: {'name': 'EGM2008 height',},
    5773: {'name': 'EGM96 height',},
    5703: {'name': 'NAVD88 height',},
    6644: {'name': 'GUVD04 height',},
    6641: {'name': 'PRVD02 height',},
    6643: {'name': 'ASVD02 height',},
    9279: {'name': 'SA LLD height',},
}

_geoids = [
    'g2018',
    'g2012b',
    'g1999',
    'geoid09',
    'geoid03',
]

class VerticalTransform:
    
    def __init__(self, src_region, src_x_inc, src_y_inc, epsg_in, epsg_out, verbose=True):
        self.src_region = src_region
        #print(self.src_region.full_region())
        #print(self.src_region.valid_p())
        #self.src_region.warp()
        #print(self.src_region.full_region())
        #print(self.src_region.valid_p())
        #print(self.src_region)
        #self.trans_region = src_region.copy()
        #self.trans_region.warp()
        self.src_x_inc = utils.str2inc(src_x_inc)
        self.src_y_inc = utils.str2inc(src_y_inc)
        self.epsg_in = self._datum_by_name(str(epsg_in))
        self.epsg_out = self._datum_by_name(str(epsg_out))

        self.verbose = verbose
        self.xcount, self.ycount, self.gt = self.src_region.geo_transform(x_inc=self.src_x_inc, y_inc=self.src_y_inc)

        self.ref_in, self.ref_out = self._frames(self.epsg_in, self.epsg_out)
        
    def _frames(self, epsg_in, epsg_out):
        
        if epsg_in in _tidal_frames.keys():
            ref_in = 'tidal'
        elif epsg_in in _htdp_reference_frames.keys():
            ref_in = 'htdp'
        elif epsg_in in _cdn_reference_frames.keys():
            ref_in = 'cdn'
        else: ref_in = None
            
        if epsg_out in _tidal_frames.keys():
            ref_out = 'tidal'
        elif epsg_out in _htdp_reference_frames.keys():
            ref_out = 'htdp'
        elif epsg_out in _cdn_reference_frames.keys():
            ref_out = 'cdn'
        else: ref_out = None

        return(ref_in, ref_out)

    def _datum_by_name(self, datum_name):
        ## tidal
        if utils.int_or(datum_name) not in _tidal_frames.keys():
            for t in _tidal_frames.keys():
                if datum_name.lower() in _tidal_frames[t]['name'].lower():
                    return(t)
        else:
            return(int(datum_name))
        ## htdp
        if utils.int_or(datum_name) not in _htdp_reference_frames.keys():
            for t in _htdp_reference_frames.keys():
                if datum_name.lower() in _htdp_reference_frames[t]['name'].lower():
                    return(t)
        else:
            return(int(datum_name))
        ## cdn
        if utils.int_or(datum_name) not in _cdn_reference_frames.keys():
            for t in _cdn_reference_frames.keys():
                if datum_name.lower() in _cdn_reference_frames[t]['name'].lower():
                    return(t)
        else:
            return(int(datum_name))
        
        return(None)
    
    def _tidal_transform(self, vdatum_tidal_in, vdatum_tidal_out):
        """generate tidal transformation grid"""

        v_in = cudem.fetches.vdatum.VDATUM(src_region=self.src_region, datatype=vdatum_tidal_in)
        v_out = cudem.fetches.vdatum.VDATUM(src_region=self.src_region, datatype=vdatum_tidal_out)

        v_in.run()
        v_out.run()

        if not v_in.results or not v_out.results:
            utils.echo_error_msg(
                'could not locate {} or {} in the region {}'.format(vdatum_tidal_in, vdatum_tidal_out, self.src_region)
            )
            return(np.zeros( (self.ycount, self.xcount) ), None)
        
        if vdatum_tidal_in != 5174 and vdatum_tidal_in != 'msl': 
            _trans_in = waffles.GMTSurface(
                data = ['vdatum:datatype={}'.format(vdatum_tidal_in)],
                src_region = self.src_region,
                xinc = self.src_x_inc,
                yinc = self.src_y_inc,
                name = '{}'.format(vdatum_tidal_in),
                verbose = True
            ).generate()
            utils.remove_glob('vdatum:datatype={}.inf'.format(vdatum_tidal_in))
        else:
            _trans_in = None

        if vdatum_tidal_out != 5174 and vdatum_tidal_out != 'msl':
            _trans_out = waffles.GMTSurface(
                data = ['vdatum:datatype={}'.format(vdatum_tidal_out)],
                src_region = self.src_region,
                xinc = self.src_x_inc,
                yinc = self.src_y_inc,
                name = '{}'.format(vdatum_tidal_out),
                verbose = True
            ).generate()
            utils.remove_glob('vdatum:datatype={}.inf'.format(vdatum_tidal_out))
                              
        else:
            _trans_out = None

        _trans_in_array, _trans_in_infos = demfun.get_array(_trans_in.fn)
        _trans_out_array, _trans_out_infos = demfun.get_array(_trans_out.fn)

        utils.remove_glob('{}*'.format(_trans_in.fn), '{}*'.format(_trans_out.fn))
        
        if _trans_in is None:
            if _trans_out is None:
                return(None)
            else:
                return(_trans_out_array*-1, self._datum_by_name(vdatum_tidal_out))
        else:
            if _trans_out is None:
                return(_trans_in_array, self._datum_by_name(vdatum_tidal_out))
            else:
                return(_trans_in_array - _trans_out_array, self._datum_by_name(vdatum_tidal_out))
            
    def _cdn_transform(self, epsg=None, name=None, invert=False):
        """create a cdn transofrmation grid"""

        if epsg is not None:
            cdn_results = cudem.fetches.vdatum.search_proj_cdn(self.src_region, epsg=epsg)
        else: cdn_results = cudem.fetches.vdatum.search_proj_cdn(self.src_region)

        for _result in cdn_results:
            for g in _geoids:
                if g in _result['name']:
                    cdn_results = [_result]
                    
        if len(cdn_results) > 0:
            for _result in cdn_results:
                src_code = int(_result['source_crs_code'].split(':')[-1])
                dst_code = int(_result['target_crs_code'].split(':')[-1])
                #if epsg == dst_code or epsg == src_code or np.any([g in _result['name'] for g in self._geoids]):
                if epsg == dst_code or np.any([g in _result['name'] for g in _geoids]):
                    if src_code in _htdp_reference_frames.keys():
                        _trans_grid = _result['name']
                        if cudem.fetches.utils.Fetch(_result['url'], verbose=self.verbose).fetch_file(_trans_grid) == 0:
                            tmp_infos = demfun.infos(_trans_grid)
                            tmp_region = regions.Region().from_geo_transform(tmp_infos['geoT'], tmp_infos['nx'], tmp_infos['ny'])
                            if os.path.exists('_{}'.format(_trans_grid)):
                                utils.remove_glob('_{}'.format(_trans_grid))
                            utils.run_cmd('gdalwarp {} {} -s_srs epsg:4326 -te {} -ts {} {} --config CENTER_LONG 0'.format(
                                _trans_grid, '_{}'.format(_trans_grid), self.src_region.format('te'), self.xcount, self.ycount
                            ), verbose=True)
                            
                            _tmp_array, _tmp_infos = demfun.get_array('_{}'.format(_trans_grid))
                            utils.remove_glob(_trans_grid, '_{}'.format(_trans_grid))
                            if invert:
                                _tmp_array = _tmp_array * -1
                            
                            return(_tmp_array, src_code)

        utils.echo_error_msg('failed to locate transformation for {}'.format(epsg))
        return(np.zeros( (self.ycount, self.xcount) ), epsg)
            
    def _htdp_transform(self, epsg_in, epsg_out):
        """create an htdp transformation grid"""

        htdp = htdpfun.HTDP()
        utils.echo_msg('{}: HTDP: {}->{}'.format(self.src_region, epsg_in, epsg_out))

        griddef = (self.src_region.xmin, self.src_region.ymax,
                   self.src_region.xmin, self.src_region.ymin,
                   self.xcount, self.ycount)
        grid = htdp._new_create_grid(griddef)

        htdp._write_grid(grid, '_tmp_input.xyz')
        htdp._write_control('_tmp_control.txt', '_tmp_output.xyz', '_tmp_input.xyz',
                            _htdp_reference_frames[epsg_in]['htdp_id'], 2012.0,
                            _htdp_reference_frames[epsg_out]['htdp_id'], 2012.0)
        htdp.run('_tmp_control.txt')

        out_grid = htdp._read_grid('_tmp_output.xyz', (griddef[5],griddef[4]))
        utils.remove_glob('_tmp_output.xyz', '_tmp_input.xyz', '_tmp_control.txt')
        return(out_grid, epsg_out)

    def _vertical_transform(self, epsg_in, epsg_out):

        trans_array = np.zeros( (self.ycount, self.xcount) )        
        while epsg_in != epsg_out and epsg_in is not None and epsg_out is not None:
            ref_in, ref_out = self._frames(epsg_in, epsg_out)
            #print(ref_in, ref_out, epsg_in, epsg_out)
            if ref_in == 'tidal':
                if ref_out == 'tidal':
                    tmp_trans, v = self._tidal_transform(_tidal_frames[epsg_in]['name'], _tidal_frames[epsg_out]['name'])
                    epsg_in = epsg_out
                else:
                    tg, tv = self._tidal_transform(_tidal_frames[self.epsg_in]['name'], 'tss')
                    cg, cv = self._cdn_transform(name='geoid', invert=False)
                    tmp_trans = tg + cg
                    epsg_in = cv
            elif ref_in == 'htdp':
                if ref_out == 'htdp':
                    tmp_trans, v = self._htdp_transform(epsg_in, epsg_out)
                    epsg_in = epsg_out
                else:
                    cg, cv = self._cdn_transform(epsg=epsg_out, invert=True)
                    hg, v = self._htdp_transform(epsg_in, cv)
                    tmp_trans = cg + hg
                    epsg_in = epsg_out
            elif ref_in == 'cdn':
                if ref_out == 'cdn':
                    tmp_trans, cv = self._cdn_transform(epsg=epsg_in, invert=False)
                    epsg_in = cv
                elif ref_out == 'tidal':
                    tmp_trans, tv = self._tidal_transform('tss', _tidal_frames[self.epsg_out]['name'])
                    epsg_in = tv
                else:
                    tmp_trans, cv = self._cdn_transform(epsg=epsg_in, invert=False)
                    epsg_in = cv
            else:
                utils.echo_error_msg('failed to locate transformation between {} and {}'.format(epsg_in, epsg_out))
                tmp_trans = np.zeros( (self.ycount, self.xcount) )
                epsg_in = epsg_out

            trans_array = trans_array + tmp_trans
            
        return(trans_array)
    
    def run(self, outfile='trans_grid.tif'):
        if self.epsg_in is None or self.epsg_out is None:
            utils.echo_error_msg('failed to parse vertical input or output, check inputs')
            return(None)
        else:
            trans_array = self._vertical_transform(self.epsg_in, self.epsg_out)
            trans_infos = demfun.set_infos(
                self.xcount, self.ycount, self.xcount*self.ycount, self.gt, None, gdal.GDT_Float32, -9999, 'GTiff'
            )

            if outfile is not None:
                utils.gdal_write(trans_array, outfile, trans_infos)
                return(outfile)
            else:
                return(trans_array, trans_infos)
        
### End
