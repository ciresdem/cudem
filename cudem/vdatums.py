### vdatums.py
##
## Copyright (c) 2021 - 2023 Regents of the University of Colorado
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
from cudem import gdalfun
from cudem import utils
from cudem import htdpfun
from cudem import fetches

_vdatums_cache = utils.cudem_cache()

_tidal_frames = {
    1089: {'name': 'mllw', 'description': 'Mean Lower Low Water'},
    5866: {'name': 'mllw', 'description': 'Mean Lower Low Water'},
    1091: {'name': 'mlw',  'description': 'Mean Low Water'},
    1090: {'name': 'mhhw', 'description': 'Mean Higher High Water'},
    5869: {'name': 'mhhw', 'description': 'Mean Higher High Water'},
    5868: {'name': 'mhw', 'description': 'Mean High Water'},
    5714: {'name': 'msl', 'description': 'Mean Sea Level'},
    5713: {'name': 'mtl', 'description': 'Mean Tide Level'},
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
    7912: {'name': 'ELLIPSOID',
           'description': 'IGS14/IGb14/WGS84 Ellipsoid',
           'htdp_id': 23},
}

_cdn_reference_frames = {
    9245: {'name': 'CGVD2013(CGG2013a) height',},
    6647: {'name': 'CGVD2013(CGG2013) height',},
    3855: {'name': 'EGM2008 height',},
    5773: {'name': 'EGM96 height',},
    5703: {'name': 'NAVD88 height',},
    6360: {'name': 'NAVD88 height (usFt)',},
    6644: {'name': 'GUVD04 height',},
    6641: {'name': 'PRVD02 height',},
    6643: {'name': 'ASVD02 height',},
    9279: {'name': 'SA LLD height',},
}

## todo: allow input/output geoids
_geoids = ['g2018', 'g2012b', 'g1999', 'geoid09', 'geoid03']
_geoids = ['g2018']

def get_vdatum_by_name(datum_name):
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

## ==============================================
## vertical transformation grid
## remove gmt dependence...maybe use IDW or cubic
## modules instead of gmt-surface
## ==============================================
class VerticalTransform:    
    def __init__(self, src_region, src_x_inc, src_y_inc, epsg_in, epsg_out, verbose=True, cache_dir=None):
        self.src_region = src_region
        self.src_x_inc = utils.str2inc(src_x_inc)
        self.src_y_inc = utils.str2inc(src_y_inc)
        self.epsg_in = self._datum_by_name(str(epsg_in))
        self.epsg_out = self._datum_by_name(str(epsg_out))
        self.cache_dir = _vdatum_cache if cache_dir is None else cache_dir
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

    def _feet_to_meters(self):
        c_array = np.zeros((self.ycount, self.xcount))
        c_array[:] = .3048
        return(c_array)

    def _meters_to_feet(self):
        c_array = np.zeros((self.ycount, self.xcount))
        c_array[:] = 3.28084
        return(c_array)
    
    def _tidal_transform(self, vdatum_tidal_in, vdatum_tidal_out):
        """generate tidal transformation grid

        This will fail over land or outside of US waters...
        """

        from cudem import waffles

        v_in = fetches.VDATUM(src_region=self.src_region, datatype=vdatum_tidal_in, verbose=self.verbose)
        v_out = fetches.VDATUM(src_region=self.src_region, datatype=vdatum_tidal_out, verbose=self.verbose)
        v_in._outdir = self.cache_dir
        v_in.run()
        v_out._outdir = self.cache_dir
        v_out.run()

        if not v_in.results:
            utils.echo_error_msg(
                'could not locate {} in the region {}'.format(
                    vdatum_tidal_in, self.src_region
                )
            )
            _trans_in_array, _trans_in_infos = [np.zeros( (self.ycount, self.xcount) ), None]
            _trans_in = None
            #return(np.zeros( (self.ycount, self.xcount) ), None)
        else:
            if vdatum_tidal_in != 5714 and vdatum_tidal_in != 'msl': 
                _trans_in = waffles.WaffleFactory(mod='IDW', data=['vdatum:datatype={}'.format(vdatum_tidal_in)], src_region=self.src_region,
                                                  xinc=self.src_x_inc, yinc=self.src_y_inc, name='{}'.format(vdatum_tidal_in),
                                                  cache_dir=self.cache_dir, dst_srs='epsg:4326', node='pixel', verbose=self.verbose)._acquire_module()
                _trans_in.initialize()
                _trans_in.generate()

                utils.remove_glob('vdatum:datatype={}.inf'.format(vdatum_tidal_in))
                _trans_in_array, _trans_in_infos = gdalfun.gdal_get_array(_trans_in.fn)
                if _trans_in_array is None:
                    _trans_in = None
                else:
                    utils.remove_glob('{}*'.format(_trans_in.fn))
            else:
                _trans_in = None
            
        if not v_out.results:
            utils.echo_error_msg(
                'could not locate {} in the region {}'.format(
                    vdatum_tidal_out, self.src_region
                )
            )
            _trans_out_array, _trans_out_infos = [np.zeros( (self.ycount, self.xcount) ), None]
            _trans_out = None
            #return(np.zeros( (self.ycount, self.xcount) ), None)
        else:
            if vdatum_tidal_out != 5714 and vdatum_tidal_out != 'msl':

                _trans_out = waffles.WaffleFactory(mod='IDW', data=['vdatum:datatype={}'.format(vdatum_tidal_out)], src_region=self.src_region,
                                                   xinc=self.src_x_inc, yinc=self.src_y_inc, name='{}'.format(vdatum_tidal_out), dst_srs='epsg:4326',
                                                   cache_dir=self.cache_dir, node='pixel', verbose=self.verbose)._acquire_module()
                _trans_out.initialize()
                _trans_out.generate()

                utils.remove_glob('vdatum:datatype={}.inf'.format(vdatum_tidal_out))
                _trans_out_array, _trans_out_infos = gdalfun.gdal_get_array(_trans_out.fn)
                if _trans_out_array is None:
                    _trans_out = None
                else:
                    utils.remove_glob('{}*'.format(_trans_out.fn))
            else:
                _trans_out = None
        
        if _trans_in is None:
            if _trans_out is None:
                return(np.zeros( (self.ycount, self.xcount) ), None)
            else:
                return(_trans_out_array*-1, self._datum_by_name(vdatum_tidal_out))
        else:
            if _trans_out is None:
                return(_trans_in_array, self._datum_by_name(vdatum_tidal_out))
            else:
                return(_trans_in_array - _trans_out_array, self._datum_by_name(vdatum_tidal_out))
            
    def _cdn_transform(self, epsg=None, name=None, invert=False):
        """create a cdn transofrmation grid"""

        epsg = 5703 if epsg == 6360 else epsg
        # c_array = None
        # if epsg == 6360:
        #     c_array = self._meters_to_feet()
        #     epsg = 5703
        
        if epsg is not None:
            cdn_results = fetches.search_proj_cdn(
                self.src_region, epsg=epsg, cache_dir=self.cache_dir, verbose=self.verbose
            )
        else:
            cdn_results = fetches.search_proj_cdn(
                self.src_region, cache_dir=self.cache_dir, verbose=self.verbose
            )

        for _result in cdn_results:
            for g in _geoids:
                if g in _result['name']:
                    print(g)
                    print(_result)
                    cdn_results = [_result]
                    break
                    
        if len(cdn_results) > 0:
            for _result in cdn_results:
                src_code = int(_result['source_crs_code'].split(':')[-1])
                dst_code = int(_result['target_crs_code'].split(':')[-1])
                #if epsg == dst_code or epsg == src_code or
                #np.any([g in _result['name'] for g in self._geoids]):
                
                if epsg == dst_code or np.any([g in _result['name'] for g in _geoids]):
                    if src_code in _htdp_reference_frames.keys():
                        _trans_grid = os.path.join(self.cache_dir, _result['name'])
                        if fetches.Fetch(
                                _result['url'], verbose=self.verbose
                        ).fetch_file(_trans_grid) == 0:
                            tmp_infos = gdalfun.gdal_infos(_trans_grid)
                            tmp_region = regions.Region().from_geo_transform(
                                tmp_infos['geoT'], tmp_infos['nx'], tmp_infos['ny']
                            )
                            if os.path.exists('_{}'.format(os.path.basename(_trans_grid))):
                                utils.remove_glob('_{}'.format(os.path.basename(_trans_grid)))
                            utils.run_cmd(
                                'gdalwarp {} {} -s_srs epsg:4326 -te {} -ts {} {} --config CENTER_LONG 0 -r cubicspline'.format(
                                    _trans_grid,
                                    '_{}'.format(os.path.basename(_trans_grid)),
                                    self.src_region.format('te'),
                                    self.xcount,
                                    self.ycount,
                                ), verbose=self.verbose
                            )
                            
                            _tmp_array, _tmp_infos = gdalfun.gdal_get_array(
                                '_{}'.format(os.path.basename(_trans_grid))
                            )
                            utils.remove_glob('_{}'.format(os.path.basename(_trans_grid)))
                            if invert:
                                _tmp_array = _tmp_array * -1

                            # if c_array is not None:
                            #     _tmp_array = _tmp_array * c_array
                                
                            return(_tmp_array, src_code)

        utils.echo_error_msg('failed to locate transformation for {}'.format(epsg))
            
        return(np.zeros( (self.ycount, self.xcount) ), epsg)
            
    def _htdp_transform(self, epsg_in, epsg_out):
        """create an htdp transformation grid"""

        htdp = htdpfun.HTDP(verbose=self.verbose)
        if self.verbose:
            utils.echo_msg('{}: HTDP: {}->{}'.format(self.src_region, epsg_in, epsg_out))

        griddef = (self.src_region.xmax, self.src_region.ymax,
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
        ## failure can cause endless loop
        while epsg_in != epsg_out and epsg_in is not None and epsg_out is not None:
            ref_in, ref_out = self._frames(epsg_in, epsg_out)
            #print(ref_in, ref_out, epsg_in, epsg_out)
            if ref_in == 'tidal':
                if ref_out == 'tidal':
                    #print(_tidal_frames)
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
            tmp_trans = None
            
        return(trans_array)
    
    def run(self, outfile=None):
        if outfile is None:
            outfile = utils.make_temp_fn('trans_grid.tif', self.cache_dir)
            
        if os.path.exists(outfile):
            if self.verbose:
                utils.echo_warning_msg('{} exists, skipping...'.format(outfile))
                
            return(outfile)
        
        if self.epsg_in is None or self.epsg_out is None:
            utils.echo_error_msg('failed to parse vertical input or output, check inputs')
                
            return(None)
        else:
            trans_array = self._vertical_transform(self.epsg_in, self.epsg_out)
            trans_infos = gdalfun.gdal_set_infos(
                self.xcount, self.ycount, self.xcount*self.ycount, self.gt, None, gdal.GDT_Float32, -9999, 'GTiff', None, None
            )

            if outfile is not None:
                gdalfun.gdal_write(trans_array, outfile, trans_infos)
                trans_array = None
                return(outfile)
            else:
                return(trans_array, trans_infos)

## ==============================================
## NOAA's VDATUM
## ==============================================
class Vdatum:
    def __init__(self, jar = None, ivert = 'navd88:m:height', overt = 'mhw:m:height',
                 ihorz = 'NAD83_2011', ohorz = 'NAD83_2011', region = '4', fmt = 'txt',
                 xyzl = '0,1,2', skip = 0, delim = 'space', result_dir = 'result',
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
        """Find the VDatum executable on the local system.
        
        returns a list of found vdatum.jar system paths
        """
        
        results = []
        for root, dirs, files in os.walk('/'):
            if 'vdatum.jar' in files:
                results.append(os.path.abspath(os.path.join(root, 'vdatum.jar')))
                break
        if len(results) == 0:
            return(None)
        else:
            self.jar = results[0]
            return(results)

    def vdatum_get_version(self):
        """run vdatum and attempt to get it's version

        return the vdatum version or None
        """

        if self.jar is None:
            self.vdatum_locate_jar()
        if self.jar is not None:
            out, status = utils.run_cmd('java -jar {} {}'.format(self.jar, '-'), verbose=self.verbose)
            for i in out.decode('utf-8').split('\n'):
                if '- v' in i.strip():
                    return(i.strip().split('v')[-1])
        return(None)

    def vdatum_xyz(self, xyz):
        """run vdatum on an xyz list [x, y, z]

        returns the transformed xyz list
        """

        if self.jar is None:
            self.vdatum_locate_jar()
        if self.jar is not None:
            vdc = 'ihorz:{} ivert:{} ohorz:{} overt:{} -nodata -pt:{},{},{} {}region:{}\
            '.format(self.ihorz, self.ivert, self.ohorz, self.overt, xyz[0], xyz[1], xyz[2], 'epoch:{} '.format(self.epoch) if self.epoch is not None else '', self.region)
            out, status = utils.run_cmd('java -Djava.awt.headless=false -jar {} {}'.format(self.jar, vdc), verbose = False)
            for i in out.split('\n'):
                if 'Height/Z' in i:
                    z = float(i.split()[2])
                    break
            return([xyz[0], xyz[1], z])
        else: return(xyz)

    def vdatum_clean_result(self):
        """clean the vdatum 'result' folder"""

        utils.remove_glob('{}/*'.format(self.result_dir))
        try:
            os.removedirs(self.result_dir)
        except: pass
    
    def run_vdatum(self, src_fn):
        """run vdatum on src_fn which is an XYZ file
        use vd_config to set vdatum parameters.

        returns [command-output, command-return-code]
        """

        if self.jar is None:
            self.vdatum_locate_jar()
        if self.jar is not None:
            vdc = 'ihorz:{} ivert:{} ohorz:{} overt:{} -nodata -file:txt:{},{},skip{}:{}:{} {}region:{}\
            '.format(self.ihorz, self.ivert, self.ohorz, self.overt, self.delim, self.xyzl, self.skip, src_fn, self.result_dir, 'epoch:{} '.format(self.epoch) if self.epoch is not None else '', self.region)
            #return(utils.run_cmd('java -Djava.awt.headless=true -jar {} {}'.format(self.jar, vdc), verbose=self.verbose))
            return(utils.run_cmd('java -jar {} {}'.format(self.jar, vdc), verbose=self.verbose))
        else: return([], -1)

## ==============================================
## Waffles VDATUM 'conversion grid' module
## U.S. Only
## ==============================================
# def waffles_vdatum(wg, ivert = 'navd88', overt = 'mhw',
#                    region = '4', jar = None):
#     """generate a 'conversion-grid' with vdatum.
    
#     output will be the differences (surfaced) between 
#     `ivert` and `overt` for the region

#     Args: 
#       wg (dict): a waffles config dictionary
#       ivert (str): input vertical datum string
#       overt (str): output vertical datum string
#       region (str): vdatum grid region
#       jar (path): path to vdatum .jar file
    
#     Returns:
#       list: [{'dem': ['dem-fn', 'raster']}, status]
#     """
    
#     vc = vdatumfun._vd_config
#     if jar is None:
#         vc['jar'] = vdatumfun.vdatum_locate_jar()[0]
#     else: vc['jar'] = jar
#     vc['ivert'] = ivert
#     vc['overt'] = overt
#     vc['region'] = region

#     gdalfun.gdal_null('empty.tif', waffles_proc_region(wg), 0.00083333, nodata = 0)
#     with open('empty.xyz', 'w') as mt_xyz:
#         for xyz in gdalfun.gdal_yield_entry(['empty.tif', 200, 1]):
#             xyzfun.xyz_line(xyz, mt_xyz, False)
    
#     vdatumfun.run_vdatum('empty.xyz', vc)
    
#     if os.path.exists('result/empty.xyz') and os.stat('result/empty.xyz').st_size != 0:
#         with open('result/empty.xyz') as infile:
#             empty_infos = xyzfun.xyz_inf(infile)
#         print(empty_infos)

#         ll = 'd' if empty_infos['minmax'][4] < 0 else '0'
#         lu = 'd' if empty_infos['minmax'][5] > 0 else '0'
#         wg['data'] = ['result/empty.xyz']
#         wg['spat'] = False
#         wg['unc'] = False
#         wg = waffles_config(**wg)
#         vd_out, status = waffles_gmt_surface(wg, tension = 0, upper_limit = lu, lower_limit = ll)
#     else:
#         utils.echo_error_msg('failed to generate VDatum grid, check settings')
#         vd_out = {}
#         status = -1
        
#     utils.remove_glob('empty.*', 'result/*', '.mjr.datalist', 'result')
#     return(vd_out, status)
        
## ==============================================
## Command-line Interface (CLI)
## $ vdatums
##
## vdatums cli
## ==============================================
_version = '0.2.0'
_epsg_desc = lambda t, x: '{}:\n '.format(t) + ' '.join(
    ['\033[1m{}\033[0m\t{}\n'.format(key, x[key]['name']) for key in x])

_usage = """{cmd} ({version}): transform a grid between vertical datums

usage: {cmd} [OPTIONS] input_grid output_grid

  input_grid\t\tThe input raster to transform
  output_grid\t\tThe output transformed raster

 Options:

  -i, --vdatum_in\tthe input vertical datum as EPSG code
  -o, --vdatum_out\tthe output vertical datum as EPSG code
  -D, --cache-dir\tCACHE Directory for storing temp data.
\t\t\tDefault Cache Directory is ~/.cudem_cache; cache will be cleared after a vdatums session
\t\t\tto retain the data, use the --keep-cache flag

  -k, --keep-cache\tKEEP the cache data intact after run
  -l, --list-epsg\tList the supported EPSG codes and their names
  -q, --quiet\t\tLower verbosity to a quiet.

  --help\t\tPrint the usage text
  --version\t\tPrint the version information

 Examples:
 % {cmd} my_dem_navd88.tif my_dem_wgs84_1674.tif --vdatum_in 5703 --vdatum_out 7662

CIRES DEM home page: <http://ciresgroups.colorado.edu/coastalDEM>
""".format(version=_version, cmd=os.path.basename(sys.argv[0]))
        
def vdatums_cli(argv = sys.argv):
    src_grid = None
    dst_grid = None
    vdatum_in = 5703
    vdatum_out = 7662
    verbose = True
    keep_cache = False
    cache_dir = _vdatums_cache
    i = 1

    while i < len(argv):
        arg = argv[i]
        if arg == '-i' or arg == '--vdatum_in':
            vdatum_in = argv[i + 1]
            i = i + 1
        elif arg == '-o' or arg == '--vdatum_out':
            vdatum_out = argv[i + 1]
            i = i + 1

        elif arg == '--cache-dir' or arg == '-D' or arg == '-cache-dir':
            cache_dir = os.path.join(utils.str_or(argv[i + 1], os.path.expanduser('~')), '.cudem_cache')
            i = i + 1
        elif arg[:2] == '-D':
            cache_dir = os.path.join(utils.str_or(argv[i + 1], os.path.expanduser('~')), '.cudem_cache')
        elif arg == '--list-epsg' or arg == '-l':
            print(_epsg_desc('htdp epsg', _htdp_reference_frames))
            print(_epsg_desc('cdn espg', _cdn_reference_frames))
            print(_epsg_desc('tidal epsg', _tidal_frames))
            sys.exit(1)
        elif arg == '-k' or arg == '--keep-cache':
            keep_cache = True
        elif arg == '--quiet':
            verbose = False
        elif arg == '-help' or arg == '--help' or arg == '-h':
            print(_usage)
            sys.exit(1)
        elif arg == '-version' or arg == '--version':
            print('vertical_datum_convert.py, version {}'.format(_version))
            sys.exit(1)
        elif src_grid is None:
            src_grid = arg
        elif dst_grid is None:
            dst_grid = arg
        else:
            sys.stderr.write(_usage)
            sys.exit(1)
        i = i + 1

    if src_grid is None or not os.path.exists(src_grid):
        sys.stderr.write(_usage)
        sys.exit(1)

    if dst_grid is None:
        dst_grid = '.'.join(src_grid.split('.')[:-1]) + '_' + str(vdatum_out.replace('(', '_').replace(')', '_')) + '.' + src_grid.split('.')[-1]

    if not os.path.exists(src_grid):
        utils.echo_error_msg('Error: {} is not a valid file'.format(src_grid))
    else:
        src_infos = gdalfun.gdal_infos(src_grid)
        
        src_region = regions.Region().from_geo_transform(src_infos['geoT'], src_infos['nx'], src_infos['ny'])
        src_region.src_srs = gdalfun.gdal_get_srs(src_grid)

        trans_region = src_region.copy()
        trans_region.warp()
        trans_region.buffer(pct=2)

        trans_region._wgs_extremes()
        
        #x_inc, y_inc = trans_region.increments(src_infos['nx']/3, src_infos['ny']/3)
        
        x_inc = src_infos['geoT'][1]
        y_inc = -src_infos['geoT'][5]
        tmp_x_inc = 3/3600
        tmp_y_inc = 3/3600
        
        vt = VerticalTransform(trans_region, tmp_x_inc, tmp_y_inc, vdatum_in, vdatum_out, cache_dir=cache_dir)
        _trans_grid = vt.run()
        out_trans_grid = utils.make_temp_fn('_trans_grid.tif')
        
        if os.path.exists(out_trans_grid):
            utils.remove_glob(out_trans_grid)
        
        if _trans_grid is not None:

            out_h, out_v = gdalfun.epsg_from_input(gdalfun.gdal_get_srs(src_grid))
            
            utils.run_cmd('gdalwarp {} {} -te {} -ts {} {} -s_srs epsg:4326 -t_srs {}'.format(
                _trans_grid, out_trans_grid,
                src_region.format('te'),
                src_infos['nx'],
                src_infos['ny'],
                out_h), verbose=verbose)
            
            # utils.run_cmd(
            #     'gdalwarp {} {} -te {} -tr {} {} -s_srs epsg:4326 -t_srs {} -co COMPRESS=LZW -co TILED=YES -co PREDICTOR=3'.format(
            #         _trans_grid,
            #         '_{}'.format(_trans_grid),
            #         src_region.format('te'),
            #         x_inc, y_inc,
            #         gdalfun.gdal_get_srs(src_grid)
            #     ), verbose=True
            # )

            # out, status = utils.run_cmd(
            #     'gdal_calc.py -A {} -B {} --calc "A+B" --outfile {} --co COMPRESS=LZW --co TILED=YES --co PREDICTOR=3 --overwrite'.format(
            #         src_grid.replace(' ', '\ '), '_{}'.format(_trans_grid).replace(' ', '\ '), dst_grid.replace(' ', '\ ')
            #     ),
            #     verbose=True
            # )
            # if status == 0:

            gdc_cmd = 'gdal_calc.py -A {} -B {} --calc "A+B" --outfile {} --co COMPRESS=LZW --co TILED=YES --co PREDICTOR=3 --overwrite'.format(
                src_grid.replace(' ', '\ '), out_trans_grid.replace(' ', '\ '), dst_grid.replace(' ', '\ '))
            os.system(gdc_cmd)
                
        else:
            utils.echo_error_msg('could not parse input/output vertical datums: {} -> {}; check spelling, etc'.format(vdatum_in, vdatum_out))

        #if not keep_cache:
        #    utils.remove_glob(cache_dir)
            
### End
