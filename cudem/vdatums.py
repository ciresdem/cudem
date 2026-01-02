### vdatums.py
##
## Copyright (c) 2021 - 2026 Regents of the University of Colorado
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
### Commentary
##
##
### Code:

import os
import sys
import argparse
import numpy as np
from osgeo import gdal, osr

# CUDEM modules
from cudem import regions
from cudem import gdalfun
from cudem import utils
from cudem import htdpfun
from cudem import fetches
from cudem import __version__ as __cudem_version__

_vdatums_cache = utils.cudem_cache()
__version__ = '0.2.5'

## ==============================================
## Vertical Datum References
## ==============================================
_tidal_frames = {
    1089: {'name': 'mllw', 'description': 'Mean Lower Low Water', 'uncertainty': 0, 'epsg': 5866},
    5866: {'name': 'mllw', 'description': 'Mean Lower Low Water', 'uncertainty': 0, 'epsg': 5866},
    1091: {'name': 'mlw',  'description': 'Mean Low Water', 'uncertainty': 0, 'epsg': 1091},
    1090: {'name': 'mhhw', 'description': 'Mean Higher High Water', 'uncertainty': 0, 'epsg': 5869},
    5869: {'name': 'mhhw', 'description': 'Mean Higher High Water', 'uncertainty': 0, 'epsg': 5869},
    5868: {'name': 'mhw', 'description': 'Mean High Water', 'uncertainty': 0, 'epsg': 5868},
    5714: {'name': 'msl', 'description': 'Mean Sea Level', 'uncertainty': 0, 'epsg': 5714},
    5713: {'name': 'mtl', 'description': 'Mean Tide Level', 'uncertainty': 0, 'epsg': 5713},
    0000: {'name': 'crd', 'description': 'Columbia River Datum', 'uncertainty': 0, 'epsg': 0000},
    1:    {'name': 'xgeoid20b', 'description': 'xgeoid 2020 B', 'uncertainty': 0, 'epsg': 1},
    7968: {'name': 'NGVD', 'uncertainty': 0, 'epsg': 7968},
}

_htdp_reference_frames = {
    4269: {'name': 'NAD_83(2011/CORS96/2007)', 'description': '(North American plate fixed)', 'htdp_id': 1, 'uncertainty': .02, 'epoch': 1997.0},
    6781: {'name': 'NAD_83(2011/CORS96/2007)', 'description': '(North American plate fixed)', 'htdp_id': 1, 'uncertainty': .02, 'epoch': 1997.0},
    6319: {'name': 'NAD_83(2011/CORS96/2007)', 'description': '(North American plate fixed)', 'htdp_id': 1, 'uncertainty': .02, 'epoch': 1997.0},
    6321: {'name': 'NAD_83(PA11/PACP00)', 'description': '(Pacific plate fixed)', 'htdp_id': 2, 'uncertainty': .02, 'epoch': 1997.0},
    6324: {'name': 'NAD_83(MA11/MARP00)', 'description': '(Mariana plate fixed)', 'htdp_id': 3, 'uncertainty': .02, 'epoch': 1997.0},
    4979: {'name': 'WGS_84(original)', 'description': '(NAD_83(2011) used)', 'htdp_id': 4, 'uncertainty': 0, 'epoch': 1997.0},
    7815: {'name': 'WGS_84(original)', 'description': '(NAD_83(2011) used)', 'htdp_id': 4, 'uncertainty': 0, 'epoch': 1997.0},
    7816: {'name': 'WGS_84(original)', 'description': '(NAD_83(2011) used)', 'htdp_id': 4, 'uncertainty': 0, 'epoch': 1997.0},
    7656: {'name': 'WGS_84(G730)', 'description': '(ITRF91 used)', 'htdp_id': 5, 'uncertainty': 0, 'epoch': 1997.0},
    7657: {'name': 'WGS_84(G730)', 'description': '(ITRF91 used)', 'htdp_id': 5, 'uncertainty': 0, 'epoch': 1997.0},
    7658: {'name': 'WGS_84(G873)', 'description': '(ITRF94 used)', 'htdp_id': 6, 'uncertainty': 0, 'epoch': 1997.0},
    7659: {'name': 'WGS_84(G873)', 'description': '(ITRF94 used)', 'htdp_id': 6, 'uncertainty': 0, 'epoch': 1997.0},
    7660: {'name': 'WGS_84(G1150)', 'description': '(ITRF2000 used)', 'htdp_id': 7, 'uncertainty': 0, 'epoch': 1997.0},
    7661: {'name': 'WGS_84(G1150)', 'description': '(ITRF2000 used)', 'htdp_id': 7, 'uncertainty': 0, 'epoch': 1997.0},
    7662: {'name': 'WGS_84(G1674)', 'description': '(ITRF2008 used)', 'htdp_id': 8, 'uncertainty': 0, 'epoch': 2000.0},
    7663: {'name': 'WGS_84(G1674)', 'description': '(ITRF2008 used)', 'htdp_id': 8, 'uncertainty': 0, 'epoch': 2000.0},
    7664: {'name': 'WGS_84(G1762)', 'description': '(IGb08 used)', 'htdp_id': 9, 'uncertainty': 0, 'epoch': 2000.0},
    7665: {'name': 'WGS_84(G1762)', 'description': '(IGb08 used)', 'htdp_id': 9, 'uncertainty': 0, 'epoch': 2000.0},
    7666: {'name': 'WGS_84(G2139)', 'description': '(ITRF2014=IGS14=IGb14 used)', 'htdp_id': 10, 'uncertainty': 0, 'epoch': 1997.0},
    7667: {'name': 'WGS_84(G2139)', 'description': '(ITRF2014=IGS14=IGb14 used)', 'htdp_id': 10, 'uncertainty': 0, 'epoch': 1997.0},
    4910: {'name': 'ITRF88', 'description': '', 'htdp_id': 11, 'uncertainty': 0, 'epoch': 1988.0},
    4911: {'name': 'ITRF89', 'description': '', 'htdp_id': 12, 'uncertainty': 0, 'epoch': 1988.0},
    7901: {'name': 'ITRF89', 'description': '', 'htdp_id': 12, 'uncertainty': 0, 'epoch': 1988.0},
    7902: {'name': 'ITRF90', 'description': '(PNEOS90/NEOS90)', 'htdp_id': 13, 'uncertainty': 0, 'epoch': 1988.0},
    7903: {'name': 'ITRF91', 'description': '', 'htdp_id': 14, 'uncertainty': 0, 'epoch': 1988.0},
    7904: {'name': 'ITRF92', 'description': '', 'htdp_id': 15, 'uncertainty': 0, 'epoch': 1988.0},
    7905: {'name': 'ITRF93', 'description': '', 'htdp_id': 16, 'uncertainty': 0, 'epoch': 1988.0},
    7906: {'name': 'ITRF94', 'description': '', 'htdp_id': 17, 'uncertainty': 0, 'epoch': 1988.0},
    7907: {'name': 'ITRF96', 'description': '', 'htdp_id': 18, 'uncertainty': 0, 'epoch': 1996.0},
    7908: {'name': 'ITRF97', 'description': 'IGS97', 'htdp_id': 19, 'uncertainty': 0, 'epoch': 1997.0},
    7909: {'name': 'ITRF2000', 'description': 'IGS00/IGb00', 'htdp_id': 20, 'uncertainty': 0, 'epoch': 2000.0},
    7910: {'name': 'ITRF2005', 'description': 'IGS05', 'htdp_id': 21, 'uncertainty': 0, 'epoch': 2000.0},
    7911: {'name': 'ITRF2008', 'description': 'IGS08/IGb08', 'htdp_id': 22, 'uncertainty': 0, 'epoch': 2000.0},
    7912: {'name': 'ELLIPSOID', 'description': 'IGS14/IGb14/WGS84/ITRF2014 Ellipsoid', 'htdp_id': 23, 'uncertainty': 0, 'epoch': 2000.0},
    1322: {'name': 'ITRF2020', 'description': 'IGS20', 'htdp_id': 24, 'uncertainty': 0, 'epoch': 2000.0},
}

_cdn_reference_frames = {
    9245: {'name': 'CGVD2013(CGG2013a) height', 'uncertainty': 0},
    6647: {'name': 'CGVD2013(CGG2013) height', 'uncertainty': 0},
    3855: {'name': 'EGM2008 height', 'uncertainty': 0},
    5773: {'name': 'EGM96 height', 'uncertainty': 0},
    5703: {'name': 'NAVD88 height', 'uncertainty': .05},
    6360: {'name': 'NAVD88 height (usFt)', 'uncertainty': .05},
    8228: {'name': 'NAVD88 height (Ft)', 'uncertainty': .05},
    6644: {'name': 'GUVD04 height', 'uncertainty': 0},
    6641: {'name': 'PRVD02 height', 'uncertainty': 0},
    6643: {'name': 'ASVD02 height', 'uncertainty': 0},
    9279: {'name': 'SA LLD height', 'uncertainty': 0},
}

_geoids = {
    'g2018':   {'name': 'geoid 2018', 'uncertainty': .0127},
    'g2012b':  {'name': 'geoid 2012b', 'uncertainty': .017},
    'g2012a':  {'name': 'geoid 2012a', 'uncertainty': .017},
    'g1999':   {'name': 'geoid 1999', 'uncertainty': .046},
    'geoid09': {'name': 'geoid 2009', 'uncertainty': .05},
    'geoid03': {'name': 'geoid 2003', 'uncertainty': .046},
}

def get_vdatum_by_name(datum_name):
    """Return the vertical datum EPSG based on the vertical datum name."""
    
    if datum_name is None:
        return None
    
    datum_int = utils.int_or(datum_name)
    
    # Check Tidal
    if datum_int not in _tidal_frames:
        for t in _tidal_frames:
            if datum_name.lower() in _tidal_frames[t]['name'].lower():
                return t
    else:
        return int(datum_name)

    ## Check HTDP
    if datum_int not in _htdp_reference_frames:
        for t in _htdp_reference_frames:
            if datum_name.lower() in _htdp_reference_frames[t]['name'].lower():
                return t
    else:
        return int(datum_name)

    ## Check CDN
    if datum_int not in _cdn_reference_frames:
        for t in _cdn_reference_frames:
            if datum_name.lower() in _cdn_reference_frames[t]['name'].lower():
                return t
    else:
        return int(datum_name)

    return None

## ==============================================
## Vertical Transformation Grid Class
## ==============================================
class VerticalTransform:
    """Generate a vertical transformation grid based on input/output vertical EPSG."""
    
    def __init__(self, mode, src_region, src_x_inc, src_y_inc, epsg_in, epsg_out,
                 geoid_in=None, geoid_out='g2018', node='pixel', verbose=True,
                 wm=None, cache_dir=None):
        self.src_region = src_region
        self.src_x_inc = utils.str2inc(src_x_inc)
        self.src_y_inc = utils.str2inc(src_y_inc)
        self.epsg_in = get_vdatum_by_name(str(epsg_in))
        self.epsg_out = get_vdatum_by_name(str(epsg_out))
        self.geoid_in = geoid_in
        self.geoid_out = geoid_out
        self.cache_dir = _vdatums_cache if cache_dir is None else cache_dir
        self.verbose = verbose
        self.xcount, self.ycount, self.gt = self.src_region.geo_transform(
            x_inc=self.src_x_inc, y_inc=self.src_y_inc, node='grid'
        )
        self.ref_in, self.ref_out = self._frames(self.epsg_in, self.epsg_out)
        self.node = node
        self.mode = mode
        self.wm = wm

        
    def _frames(self, epsg_in, epsg_out):
        """Determine the input/output vertical datum frame."""
        
        ref_in = None
        ref_out = None
        
        if epsg_in in _tidal_frames:
            ref_in = 'tidal'
        elif epsg_in in _htdp_reference_frames:
            ref_in = 'htdp'
        elif epsg_in in _cdn_reference_frames:
            ref_in = 'cdn'
            
        if epsg_out in _tidal_frames:
            ref_out = 'tidal'
        elif epsg_out in _htdp_reference_frames:
            ref_out = 'htdp'
        elif epsg_out in _cdn_reference_frames:
            ref_out = 'cdn'

        return ref_in, ref_out

    
    def _tidal_transform(self, vdatum_tidal_in, vdatum_tidal_out):
        """Generate tidal transformation grid using NOAA VDatum."""

        from cudem.waffles.waffles import WaffleFactory
        
        ## Fetch the input tidal datum from VDatum
        v_in = fetches.vdatum.VDATUM(src_region=self.src_region, datatype=vdatum_tidal_in, verbose=self.verbose)
        v_in._outdir = self.cache_dir
        v_in.run()

        if vdatum_tidal_out is not None:
            ## Fetch the output tidal datum from VDatum if necessary
            v_out = fetches.vdatum.VDATUM(src_region=self.src_region, datatype=vdatum_tidal_out, verbose=self.verbose)
            v_out._outdir = self.cache_dir
            v_out.run()

        _trans_in = None
        _trans_in_array = np.zeros((self.ycount, self.xcount))
        
        if not v_in.results:
            utils.echo_error_msg(f'could not locate {vdatum_tidal_in} in the region {self.src_region}')
        else:
            if vdatum_tidal_in not in [5714, 'msl']: 
                _trans_in = WaffleFactory(
                    mod=self.mode, data=[f'vdatum:datatype={vdatum_tidal_in}'], 
                    src_region=self.src_region, xinc=self.src_x_inc, yinc=self.src_y_inc, 
                    name=str(vdatum_tidal_in), cache_dir=self.cache_dir, dst_srs='epsg:4326', 
                    node='pixel', verbose=self.verbose
                )._acquire_module()
                _trans_in.initialize()
                _trans_in.generate()
                
                _trans_in_array, _ = gdalfun.gdal_get_array(_trans_in.fn)
                utils.remove_glob(f'vdatum:datatype={vdatum_tidal_in}.inf')
                
                if _trans_in_array is None:
                    _trans_in = None
                elif _trans_in is not None:
                    for ofk in _trans_in.output_files:
                        utils.remove_glob(f'{_trans_in.output_files[ofk]}*')

        _trans_out = None
        _trans_out_array = np.zeros((self.ycount, self.xcount))

        if vdatum_tidal_out is None or not v_out.results:
            utils.echo_warning_msg(f'could not locate {vdatum_tidal_out} in the region {self.src_region}')
        else:
            if vdatum_tidal_out not in [5714, 'msl']:
                _trans_out = WaffleFactory(
                    mod=self.mode, data=[f'vdatum:datatype={vdatum_tidal_out}'], 
                    src_region=self.src_region, xinc=self.src_x_inc, yinc=self.src_y_inc, 
                    name=str(vdatum_tidal_out), dst_srs='epsg:4326', cache_dir=self.cache_dir, 
                    node='pixel', verbose=self.verbose
                )._acquire_module()
                _trans_out.initialize()
                _trans_out.generate()
                
                _trans_out_array, _ = gdalfun.gdal_get_array(_trans_out.fn)
                utils.remove_glob(f'vdatum:datatype={vdatum_tidal_out}.inf')

                if _trans_out_array is None:
                    _trans_out = None
                elif _trans_out is not None:
                    for ofk in _trans_out.output_files:
                        utils.remove_glob(f'{_trans_out.output_files[ofk]}*')

        if _trans_in is None:
            if _trans_out is None:
                return np.zeros((self.ycount, self.xcount)), None
            else:
                return _trans_out_array * -1, get_vdatum_by_name(vdatum_tidal_out)
        else:
            if _trans_out is None:
                return _trans_in_array, get_vdatum_by_name('msl')
            else:
                return _trans_in_array - _trans_out_array, get_vdatum_by_name(vdatum_tidal_out)

            
    def _cdn_transform(self, epsg=None, name=None, geoid='g2018', invert=False):
        """Create a CDN transformation grid (Proj)."""
        
        epsg = 5703 if (epsg == 6360 or epsg == 8228) else epsg
        geoid = 'g2018' if geoid is None else geoid
        
        ## Fetch the CDN transformation grids
        if epsg is not None:
            cdn_results = fetches.vdatum.search_proj_cdn(
                self.src_region, epsg=epsg, cache_dir=self.cache_dir, verbose=self.verbose
            )
        else:
            cdn_results = fetches.vdatum.search_proj_cdn(
                self.src_region, cache_dir=self.cache_dir, verbose=self.verbose
            )

        ## Filter by geoid name
        if len(cdn_results) > 0:
            for _result in cdn_results:
                if geoid in _result['name']:
                    cdn_results = [_result]
                    break
            
            for _result in cdn_results:
                src_code = _result.get('source_crs_code')
                dst_code = _result.get('target_crs_code')
                
                if src_code:
                    src_code = int(src_code.split(':')[-1])
                if dst_code:
                    dst_code = int(dst_code.split(':')[-1])
                    
                if epsg == dst_code or np.any([g in _result['name'] for g in _geoids]):
                    if src_code in _htdp_reference_frames:
                        _trans_grid = os.path.join(self.cache_dir, _result['name'])
                        
                        fetcher = fetches.fetches.Fetch(_result['url'], verbose=self.verbose)
                        if fetcher.fetch_file(_trans_grid) == 0:
                            
                            tmp_warp_name = f'_{os.path.basename(_trans_grid)}'
                            if os.path.exists(tmp_warp_name):
                                utils.remove_glob(tmp_warp_name)
                                
                            warp_cmd = (
                                f'gdalwarp "{_trans_grid}" "{tmp_warp_name}" '
                                f'-s_srs epsg:4326 -te {self.src_region.format("te")} '
                                f'-ts {self.xcount} {self.ycount} --config CENTER_LONG 0 '
                                f'-r cubicspline {"-wm " + str(self.wm) if self.wm else ""}'
                            )
                            utils.run_cmd(warp_cmd, verbose=self.verbose)
                            
                            _tmp_array, _ = gdalfun.gdal_get_array(tmp_warp_name)
                            utils.remove_glob(tmp_warp_name)
                            
                            if invert:
                                _tmp_array = _tmp_array * -1

                            return _tmp_array, src_code
                        
        utils.echo_error_msg(f'failed to locate transformation for {epsg}')
        return np.zeros((self.ycount, self.xcount)), epsg

    
    def _htdp_transform(self, epsg_in, epsg_out):
        """Create an HTDP transformation grid."""
        
        if utils.config_check()['HTDP'] is None:
            utils.echo_error_msg('you must have HTDP installed to perform HTDP vertical transformations')            
            return np.zeros((self.ycount, self.xcount)), epsg_out
        
        htdp = htdpfun.HTDP(verbose=self.verbose)
        utils.echo_msg(f'{self.src_region}: HTDP: {epsg_in}->{epsg_out}')

        griddef = (self.src_region.xmax, self.src_region.ymax,
                   self.src_region.xmin, self.src_region.ymin,
                   self.xcount, self.ycount)

        grid = htdp._new_create_grid(griddef)
        htdp._write_grid(grid, '_tmp_input.xyz')
        htdp._write_control('_tmp_control.txt', '_tmp_output.xyz', '_tmp_input.xyz',
                            _htdp_reference_frames[epsg_in]['htdp_id'], 1997., 
                            _htdp_reference_frames[epsg_out]['htdp_id'], 1997.)
        htdp.run('_tmp_control.txt')

        out_grid = htdp._read_grid('_tmp_output.xyz', (griddef[5], griddef[4]))
        utils.remove_glob('_tmp_output.xyz', '_tmp_input.xyz', '_tmp_control.txt')
        return out_grid, epsg_out

    
    def _vertical_transform(self, epsg_in, epsg_out):
        """Perform the combined vertical transformation."""
        
        trans_array = np.zeros((self.ycount, self.xcount))
        unc_array = np.zeros((self.ycount, self.xcount))
        
        if self.epsg_in == 0:
            self.geoid_in = 'geoid09'
            
        if self.geoid_in is not None:
            if self.epsg_in == 0:
                tg, _ = self._tidal_transform(_tidal_frames[self.epsg_in]['name'], None)
            else:
                tg = np.zeros((self.ycount, self.xcount))
                
            unc_array = np.sqrt(unc_array**2 + _geoids[self.geoid_in]['uncertainty']**2)
            tmp_trans_geoid, epsg_in = self._cdn_transform(name='geoid', geoid=self.geoid_in, invert=False)
            tmp_trans_geoid += tg
        else:
            tmp_trans_geoid = np.zeros((self.ycount, self.xcount))

        utils.echo_debug_msg(f'v-datum in: {epsg_in} v-datum_out: {epsg_out}')
        
        while epsg_in != epsg_out and epsg_in is not None and epsg_out is not None:
            ref_in, ref_out = self._frames(epsg_in, epsg_out)
            
            if ref_in == 'tidal':
                if ref_out == 'tidal':
                    tmp_trans, _ = self._tidal_transform(_tidal_frames[epsg_in]['name'], _tidal_frames[epsg_out]['name'])
                    epsg_in = epsg_out
                elif epsg_in == 7968:
                    tg, _ = self._tidal_transform(_tidal_frames[self.epsg_in]['name'], None)
                    cg, cv = self._cdn_transform(name='geoid', geoid=self.geoid_out, invert=False)
                    tmp_trans = tg + cg
                    epsg_in = cv                    
                else:
                    tg, _ = self._tidal_transform(_tidal_frames[self.epsg_in]['name'], 'tss')
                    cg, cv = self._cdn_transform(name='geoid', geoid=self.geoid_out, invert=False)
                    tmp_trans = tg + cg
                    epsg_in = cv
            elif ref_in == 'htdp':
                if ref_out == 'htdp':
                    unc_array = np.sqrt(unc_array**2 + float(_htdp_reference_frames[epsg_in]['uncertainty'] ** 2))
                    tmp_trans, _ = self._htdp_transform(epsg_in, epsg_out)
                    epsg_in = epsg_out
                else:
                    unc_array = np.sqrt(unc_array**2 + float(_htdp_reference_frames[epsg_in]['uncertainty'] ** 2))
                    unc_array = np.sqrt(unc_array**2 + float(_cdn_reference_frames[epsg_out]['uncertainty'] ** 2))
                    cg, cv = self._cdn_transform(epsg=epsg_out, invert=True)
                    unc_array = np.sqrt(unc_array**2 + float(_geoids[self.geoid_out]['uncertainty'] ** 2))
                    hg, _ = self._htdp_transform(epsg_in, cv)
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
                utils.echo_error_msg(f'failed to locate transformation between {epsg_in} and {epsg_out}')
                tmp_trans = np.zeros((self.ycount, self.xcount))
                epsg_in = epsg_out
            
            trans_array = trans_array + tmp_trans + tmp_trans_geoid
            
        return trans_array, unc_array

    
    def run(self, outfile=None):
        if outfile is None:
            outfile = utils.make_temp_fn('trans_grid.tif', self.cache_dir)

        if self.verbose:
            utils.echo_msg(outfile)
        
        if os.path.exists(outfile):
            if self.verbose:
                utils.echo_warning_msg(f'{outfile} exists, skipping...')
            return outfile

        unc_outfile = f'{utils.fn_basename2(outfile)}_unc.{utils.fn_ext(outfile)}'

        if self.verbose:
            utils.echo_msg(unc_outfile)
            utils.echo_msg(f'{self.epsg_in} {self.epsg_out}')
            
        if self.epsg_in is None or self.epsg_out is None:
            utils.echo_error_msg(f'failed to parse vertical input or output: {self.epsg_in}->{self.epsg_out}')
            return None
        
        trans_array, unc_array = self._vertical_transform(self.epsg_in, self.epsg_out)
        trans_infos = gdalfun.gdal_set_infos(
            self.xcount, self.ycount, self.xcount*self.ycount, self.gt, None, gdal.GDT_Float32, -9999, 'GTiff', None, None
        )

        if outfile is not None:
            gdalfun.gdal_write(trans_array, outfile, trans_infos)
            gdalfun.gdal_write(unc_array, unc_outfile, trans_infos)
            return outfile, unc_outfile
        else:
            return trans_array, unc_array, trans_infos

        
## ==============================================
## NOAA's VDATUM Wrapper
## ==============================================
class Vdatum:
    def __init__(self, jar=None, ivert='navd88:m:height', overt='mhw:m:height',
                 ihorz='NAD83_2011', ohorz='NAD83_2011', region='4', fmt='txt',
                 xyzl='0,1,2', skip=0, delim='space', result_dir='result',
                 verbose=False):
        self.jar = jar
        self.ivert = ivert
        self.overt = overt
        self.ihorz = ihorz
        self.ohorz = ohorz
        self.region = region
        self.fmt = fmt
        self.xyzl = xyzl
        self.skip = skip
        self.delim = delim
        self.result_dir = result_dir
        self.verbose = verbose
        self.epoch = None
        self.vdatum_set_horz()

        
    def vdatum_set_horz(self):
        if 'ITRF' in self.overt:
            self.ohorz = self.overt
            self.epoch = '1997.0:1997.0'

            
    def vdatum_locate_jar(self):
        """Find the VDatum executable on the local system."""
        
        results = []
        for root, dirs, files in os.walk('/'):
            if 'vdatum.jar' in files:
                results.append(os.path.abspath(os.path.join(root, 'vdatum.jar')))
                break
        if not results:
            return None
        else:
            self.jar = results[0]
            return results

        
    def vdatum_get_version(self):
        """Run vdatum and attempt to get its version."""
        
        if self.jar is None:
            self.vdatum_locate_jar()
        if self.jar is not None:
            out, _ = utils.run_cmd(f'java -jar {self.jar} -', verbose=self.verbose)
            for i in out.decode('utf-8').split('\n'):
                if '- v' in i.strip():
                    return i.strip().split('v')[-1]
        return None

    
    def vdatum_xyz(self, xyz):
        """Run vdatum on an xyz list [x, y, z]."""
        
        if self.jar is None:
            self.vdatum_locate_jar()
        if self.jar is not None:
            epoch_str = f'epoch:{self.epoch} ' if self.epoch is not None else ''
            vdc = (f'ihorz:{self.ihorz} ivert:{self.ivert} ohorz:{self.ohorz} overt:{self.overt} '
                   f'-nodata -pt:{xyz[0]},{xyz[1]},{xyz[2]} {epoch_str}region:{self.region}')
            
            out, _ = utils.run_cmd(
                f'java -Djava.awt.headless=false -jar {self.jar} {vdc}',
                verbose=False
            )
            z = xyz[2]
            for i in out.split('\n'):
                if 'Height/Z' in i:
                    try:
                        z = float(i.split()[2])
                        break
                    except ValueError:
                        pass
            return [xyz[0], xyz[1], z]
        else: 
            return xyz

        
    def vdatum_clean_result(self):
        """Clean the vdatum 'result' folder."""
        
        utils.remove_glob(f'{self.result_dir}/*')
        try:
            os.removedirs(self.result_dir)
        except OSError: 
            pass

        
    def run_vdatum(self, src_fn):
        """Run vdatum on src_fn which is an XYZ file."""
        
        if self.jar is None:
            self.vdatum_locate_jar()
        if self.jar is not None:
            epoch_str = f'epoch:{self.epoch} ' if self.epoch is not None else ''
            vdc = (f'ihorz:{self.ihorz} ivert:{self.ivert} ohorz:{self.ohorz} overt:{self.overt} '
                   f'-nodata -file:txt:{self.delim},{self.xyzl},skip{self.skip}:{src_fn}:{self.result_dir} '
                   f'{epoch_str}region:{self.region}')
            return utils.run_cmd(f'java -jar {self.jar} {vdc}', verbose=self.verbose)
        else: 
            return [], -1


## ==============================================
## Command-line Interface (CLI)
## $ vdatums
##
## vdatums cli
## ==============================================
def vdatums_cli():
    parser = argparse.ArgumentParser(
        description=f'%(prog)s ({__version__}): Transform a grid between vertical datums',
        epilog="CUDEM home page: <http://cudem.colorado.edu>"
    )
    
    parser.add_argument('input_grid', help='The input raster to transform')
    parser.add_argument('output_grid', nargs='?', help='The output transformed raster')
    parser.add_argument('-i', '--vdatum_in', default=5703, help='The input vertical datum as EPSG code or name')
    parser.add_argument('-o', '--vdatum_out', default=7662, help='The output vertical datum as EPSG code or name')
    parser.add_argument('-D', '--cache-dir', help='Cache Directory for storing temp data.')
    parser.add_argument('-k', '--keep-cache', action='store_true', help='Keep the cache data intact after run')
    parser.add_argument('-l', '--list-epsg', action='store_true', help='List the supported EPSG codes and their names')
    parser.add_argument('-q', '--quiet', action='store_true', help='Lower verbosity')
    parser.add_argument('-w', '--warp-mem', help='Warp memory for gdalwarp')
    parser.add_argument('--version', action='version', version=f'CUDEM {__cudem_version__} :: %(prog)s {__version__}')
    
    args = parser.parse_args()

    if args.list_epsg:
        def _print_epsg(title, data):
            print(f'{title}:')
            for key, val in data.items():
                print(f'  \033[1m{key}\033[0m\t{val["name"]}')
        
        _print_epsg('HTDP EPSG', _htdp_reference_frames)
        _print_epsg('CDN EPSG', _cdn_reference_frames)
        _print_epsg('Tidal EPSG', _tidal_frames)
        sys.exit(0)

    ## Set up cache directory
    cache_dir = args.cache_dir
    if not cache_dir:
        cache_dir = os.path.join(os.path.expanduser('~'), 'cudem_cache')

    src_grid = args.input_grid
    if not os.path.exists(src_grid):
        utils.echo_error_msg(f'Error: {src_grid} is not a valid file')
        sys.exit(1)

    ## Set default output grid name if not provided
    if not args.output_grid:
        base, ext = os.path.splitext(src_grid)
        v_out_str = str(args.vdatum_out).replace('(', '_').replace(')', '_')
        dst_grid = f"{base}_{v_out_str}{ext}"
    else:
        dst_grid = args.output_grid

    ## Process grid
    src_infos = gdalfun.gdal_infos(src_grid)
    src_region = regions.Region().from_geo_transform(src_infos['geoT'], src_infos['nx'], src_infos['ny'])
    src_horz, _ = gdalfun.split_srs(gdalfun.gdal_get_srs(src_grid))
    
    trans_region = src_region.copy()
    if src_horz is not None:
        src_region.src_srs = src_horz
        trans_region.src_srs = src_horz
        trans_region.warp()

    if not trans_region.valid_p(check_xy=True):
        utils.echo_warning_msg(f'failed to transform source region {src_region}!')
        trans_region = src_region.copy()
        
    trans_region.buffer(pct=2)
    trans_region._wgs_extremes()
    
    tmp_x_inc = 3/3600
    tmp_y_inc = 3/3600

    vt = VerticalTransform(
        'IDW', trans_region, tmp_x_inc, tmp_y_inc,
        args.vdatum_in, args.vdatum_out, wm=args.wm, cache_dir=cache_dir,
        verbose=not args.quiet
    )
    
    _trans_grid, _ = vt.run()
    out_trans_grid = utils.make_temp_fn('_trans_grid.tif', cache_dir)
    
    if os.path.exists(out_trans_grid):
        utils.remove_glob(out_trans_grid)
    
    if _trans_grid is not None:
        out_h, _ = gdalfun.epsg_from_input(gdalfun.gdal_get_srs(src_grid))
        
        warp_cmd = (
            f'gdalwarp "{_trans_grid}" "{out_trans_grid}" -te {src_region.format("te")} '
            f'-ts {src_infos["nx"]} {src_infos["ny"]} -s_srs epsg:4326 -t_srs {out_h} '
            f'{"-wm " + str(args.wm) if args.wm else ""}'
        )
        utils.run_cmd(warp_cmd, verbose=not args.quiet)

        gdc_cmd = (
            f'gdal_calc.py -A {src_grid.replace(" ", r"\ ")} -B {out_trans_grid.replace(" ", r"\ ")} '
            f'--calc "A+B" --outfile {dst_grid.replace(" ", r"\ ")} --co COMPRESS=LZW '
            f'--co TILED=YES --co PREDICTOR=3 --overwrite'
        )
        
        utils.run_cmd(gdc_cmd, verbose=not args.quiet)

        ## Set Metadata
        out_horz_srs = osr.SpatialReference()
        out_horz_srs.SetFromUserInput(src_horz)
        out_vert_srs = osr.SpatialReference()
        out_vert_srs.SetFromUserInput(f'epsg:{args.vdatum_out}')
        out_src_srs = osr.SpatialReference()
        out_src_srs.SetCompoundCS(
            'Combined', out_horz_srs, out_vert_srs
        )
        
        gdalfun.gdal_set_srs(dst_grid.replace(' ', r'\ '), out_src_srs.ExportToWkt())
    else:
        utils.echo_error_msg(
            f'could not parse input/output vertical datums: {args.vdatum_in} -> {args.vdatum_out}'
        )

    if not args.keep_cache:
        pass

if __name__ == '__main__':
    vdatums_cli()

    
### End
