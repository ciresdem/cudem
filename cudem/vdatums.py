### vdatums.py
##
## Copyright (c) 2021 - 2024 Regents of the University of Colorado
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
from osgeo import osr
import pyproj
from tqdm import tqdm

from cudem import regions
from cudem import gdalfun
from cudem import utils
from cudem import htdpfun
from cudem import fetches

_vdatums_cache = utils.cudem_cache()

## ==============================================
## vertical datum references
## ==============================================
_tidal_frames = {
    1089: {'name': 'mllw', 'description': 'Mean Lower Low Water',
           'uncertainty': 0},
    5866: {'name': 'mllw', 'description': 'Mean Lower Low Water',
           'uncertainty': 0},
    1091: {'name': 'mlw',  'description': 'Mean Low Water',
           'uncertainty': 0},
    1090: {'name': 'mhhw', 'description': 'Mean Higher High Water',
           'uncertainty': 0},
    5869: {'name': 'mhhw', 'description': 'Mean Higher High Water',
           'uncertainty': 0},
    5868: {'name': 'mhw', 'description': 'Mean High Water',
           'uncertainty': 0},
    5714: {'name': 'msl', 'description': 'Mean Sea Level',
           'uncertainty': 0},
    5713: {'name': 'mtl', 'description': 'Mean Tide Level',
           'uncertainty': 0},
    0000: {'name': 'crd', 'description': 'Columbia River Datum',
           'uncertainty': 0},
}

_htdp_reference_frames = {
    4269: {'name': 'NAD_83(2011/CORS96/2007)',
           'description': '(North American plate fixed)',
           'htdp_id': 1,
           'uncertainty': .02},
    6781: {'name': 'NAD_83(2011/CORS96/2007)',
           'description': '(North American plate fixed)',
           'htdp_id': 1,
           'uncertainty': .02},
    6319: {'name': 'NAD_83(2011/CORS96/2007)',
           'description': '(North American plate fixed)',
           'htdp_id': 1,
           'uncertainty': .02},
    6321: {'name': 'NAD_83(PA11/PACP00)',
           'description': '(Pacific plate fixed)',
           'htdp_id': 2,
           'uncertainty': .02},
    6324: {'name': 'NAD_83(MA11/MARP00)',
           'description': '(Mariana plate fixed)',
           'htdp_id': 3,
           'uncertainty': .02},
    4979: {'name': 'WGS_84(original)',
           'description': '(NAD_83(2011) used)',
           'htdp_id': 4,
           'uncertainty': 0},
    7815: {'name': 'WGS_84(original)',
           'description': '(NAD_83(2011) used)',
           'htdp_id': 4,
           'uncertainty': 0},
    7816: {'name': 'WGS_84(original)',
           'description': '(NAD_83(2011) used)',
           'htdp_id': 4,
           'uncertainty': 0},
    7656: {'name': 'WGS_84(G730)',
           'description': '(ITRF91 used)',
           'htdp_id': 5,
           'uncertainty': 0},
    7657: {'name': 'WGS_84(G730)',
           'description': '(ITRF91 used)',
           'htdp_id': 5,
           'uncertainty': 0},
    7658: {'name': 'WGS_84(G873)',
           'description': '(ITRF94 used)',
           'htdp_id': 6,
           'uncertainty': 0},
    7659: {'name': 'WGS_84(G873)',
           'description': '(ITRF94 used)',
           'htdp_id': 6,
           'uncertainty': 0},
    7660: {'name': 'WGS_84(G1150)',
           'description': '(ITRF2000 used)',
           'htdp_id': 7,
           'uncertainty': 0},
    7661: {'name': 'WGS_84(G1150)',
           'description': '(ITRF2000 used)',
           'htdp_id': 7,
           'uncertainty': 0},
    7662: {'name': 'WGS_84(G1674)',
           'description': '(ITRF2008 used)',
           'htdp_id': 8,
           'uncertainty': 0},
    7663: {'name': 'WGS_84(G1674)',
           'description': '(ITRF2008 used)',
           'htdp_id': 8,
           'uncertainty': 0},
    7664: {'name': 'WGS_84(G1762)',
           'description': '(IGb08 used)',
           'htdp_id': 9,
           'uncertainty': 0},
    7665: {'name': 'WGS_84(G1762)',
           'description': '(IGb08 used)',
           'htdp_id': 9,
           'uncertainty': 0},
    7666: {'name': 'WGS_84(G2139)',
           'description': '(ITRF2014=IGS14=IGb14 used)',
           'htdp_id': 10,
           'uncertainty': 0},
    7667: {'name': 'WGS_84(G2139)',
           'description': '(ITRF2014=IGS14=IGb14 used)',
           'htdp_id': 10,
           'uncertainty': 0},
    4910: {'name': 'ITRF88',
           'description': '',
           'htdp_id': 11,
           'uncertainty': 0},
    4911: {'name': 'ITRF89',
           'description': '',
           'htdp_id': 12,
           'uncertainty': 0},
    7901: {'name': 'ITRF89',
           'description': '',
           'htdp_id': 12,
           'uncertainty': 0},
    7902: {'name': 'ITRF90',
           'description': '(PNEOS90/NEOS90)',
           'htdp_id': 13,
           'uncertainty': 0},
    7903: {'name': 'ITRF91',
           'description': '',
           'htdp_id': 14,
           'uncertainty': 0},
    7904: {'name': 'ITRF92',
           'description': '',
           'htdp_id': 15,
           'uncertainty': 0},
    7905: {'name': 'ITRF93',
           'description': '',
           'htdp_id': 16,
           'uncertainty': 0},
    7906: {'name': 'ITRF94',
           'description': '',
           'htdp_id': 17,
           'uncertainty': 0},
    7907: {'name': 'ITRF96',
           'description': '',
           'htdp_id': 18,
           'uncertainty': 0},
    7908: {'name': 'ITRF97',
           'description': 'IGS97',
           'htdp_id': 19,
           'uncertainty': 0},
    7909: {'name': 'ITRF2000',
           'description': 'IGS00/IGb00',
           'htdp_id': 20,
           'uncertainty': 0},
    7910: {'name': 'ITRF2005',
           'description': 'IGS05',
           'htdp_id': 21,
           'uncertainty': 0},
    7911: {'name': 'ITRF2008',
           'description': 'IGS08/IGb08',
           'htdp_id': 22,
           'uncertainty': 0},
    7912: {'name': 'ITRF2014',
           'description': 'IGS14/IGb14',
           'htdp_id': 23,
           'uncertainty': 0},
    1322: {'name': 'ITRF2020',
           'description': 'IGS20',
           'htdp_id': 24,
           'uncertainty': 0},
    7912: {'name': 'ELLIPSOID',
           'description': 'IGS14/IGb14/WGS84/ITRF2014 Ellipsoid',
           'htdp_id': 23,
           'uncertainty': 0},
}

_cdn_reference_frames = {
    9245: {'name': 'CGVD2013(CGG2013a) height',
           'uncertainty': 0},
    6647: {'name': 'CGVD2013(CGG2013) height',
           'uncertainty': 0},
    3855: {'name': 'EGM2008 height',
           'uncertainty': 0},
    5773: {'name': 'EGM96 height',
           'uncertainty': 0},
    5703: {'name': 'NAVD88 height',
           'uncertainty': .05},
    6360: {'name': 'NAVD88 height (usFt)',
           'uncertainty': .05},
    6644: {'name': 'GUVD04 height',
           'uncertainty': 0},
    6641: {'name': 'PRVD02 height',
           'uncertainty': 0},
    6643: {'name': 'ASVD02 height',
           'uncertainty': 0},
    9279: {'name': 'SA LLD height',
           'uncertainty': 0},
}

_geoids = {
    'g2018': {'name': 'geoid 2018',
              'uncertainty': .0127},
    'g2012b': {'name': 'geoid 2012b',
               'uncertainty': .017},
    'g2012a': {'name': 'geoid 2012a',
               'uncertainty': .017},
    'g1999': {'name': 'geoid 1999',
              'uncertainty': .046},
    'geoid09': {'name': 'geoid 2009',
                'uncertainty': .05},
    'geoid03': {'name': 'geoid 2003',
                'uncertainty': .046},
}

## todo: allow input/output geoids
#_geoids = ['g2018', 'g2012b', 'g1999', 'geoid09', 'geoid03']
#_geoids = ['g2018']

def set_transform(src_srs = None, dst_srs = None, region = None, verbose = True, cache_dir = None):
    """Set the pyproj horizontal and vertical transformations for the dataset"""

    want_vertical = True
    src_geoid = None
    dst_geoid = 'g2018'
    cache_dir = _vdatums_cache if cache_dir is None else cache_dir
    if region is None:
        utils.echo_error_msg('you must suuply a reigon')
        sys.exit(1)
    
    if src_srs is not None and dst_srs is not None:
        tmp_src_srs = src_srs.split('+geoid:')
        src_srs_ = tmp_src_srs[0]
        src_srs = src_srs_
        if len(tmp_src_srs) > 1:
            src_geoid = tmp_src_srs[1]

        tmp_dst_srs = dst_srs.split('+geoid:')
        dst_srs_ = tmp_dst_srs[0]
        dst_srs = dst_srs_
        if len(tmp_dst_srs) > 1:
            dst_geoid = tmp_dst_srs[1]

        is_esri = False
        in_vertical_epsg_esri = None
        if 'ESRI' in src_srs_.upper():
            is_esri = True
            srs_split = src_srs_.split('+')
            src_srs_ = srs_split[0]
            if len(srs_split) > 1:
                in_vertical_epsg_esri = srs_split[1]

        in_crs = pyproj.CRS.from_user_input(src_srs_)
        out_crs = pyproj.CRS.from_user_input(dst_srs_)

        if in_crs.is_compound:
            in_crs_list = in_crs.sub_crs_list
            in_horizontal_crs = in_crs_list[0]
            in_vertical_crs = in_crs_list[1]
            in_vertical_name = in_vertical_crs.name
            in_vertical_epsg = in_vertical_crs.to_epsg()
            if in_vertical_epsg is None:
                in_vertical_epsg = in_vertical_name
        else:
            in_horizontal_crs = in_crs
            want_vertical=False

        if out_crs.is_compound:            
            out_crs_list = out_crs.sub_crs_list
            out_horizontal_crs = out_crs_list[0]
            out_vertical_crs = out_crs_list[1]
            out_vertical_epsg = out_vertical_crs.to_epsg()
        else:
            out_horizontal_crs = out_crs
            want_vertical=False

        in_horizontal_epsg = in_horizontal_crs.to_epsg()
        out_horizontal_epsg = out_horizontal_crs.to_epsg()

        if (in_vertical_epsg_esri is not None and is_esri):
            in_vertical_epsg = in_vertical_epsg_esri
            want_vertical = True

        if want_vertical:
            if (in_vertical_epsg == out_vertical_epsg) and src_geoid is None:
                want_vertical = False

        ## horizontal Transformation
        try:
            transformer = pyproj.Transformer.from_crs(in_horizontal_crs, out_horizontal_crs, always_xy=True)
        except Exception as e:
            utils.echo_warning_msg('could not set transformation in: {}, out: {}, {}'.format(
                in_horizontal_crs.name, out_horizontal_crs.name, e
            ))
            transformer = None
            return

        #if region is not None:
        trans_region = region.copy()
        trans_region.src_srs = out_horizontal_crs.to_proj4() #'epsg:{}'.format(out_horizontal_epsg)
        trans_region.warp(in_horizontal_crs.to_proj4())

        src_proj4 = in_horizontal_crs.to_proj4()
        dst_proj4 = out_horizontal_crs.to_proj4()

        ## vertical Transformation
        if want_vertical:
            vd_region = region.copy()
            vd_region.src_srs = out_horizontal_crs.to_proj4()

            vd_region.warp('epsg:4326')
            if not vd_region.valid_p():
                utils.echo_warning_msg('failed to generate transformation')
                return

            vd_region.zmin = None
            vd_region.zmax = None
            vd_region.buffer(pct=5)

            ## trans_fn is the transformation grid, used in gdalwarp
            trans_fn = os.path.join(
                cache_dir, '_vdatum_trans_{}_{}_{}.tif'.format(
                    in_vertical_epsg, out_vertical_epsg, vd_region.format('fn')
                )
            )

            ## vertical transformation grid is generated in WGS84
            if not os.path.exists(trans_fn):
                with tqdm(
                        desc='generating vertical transformation grid {} from {} to {}'.format(
                            trans_fn, in_vertical_epsg, out_vertical_epsg
                        ),
                        leave=verbose
                ) as pbar:
                    vd_x_inc = vd_y_inc = utils.str2inc('3s')
                    xcount, ycount, dst_gt = vd_region.geo_transform(
                        x_inc=vd_x_inc, y_inc=vd_y_inc, node='grid'
                    )

                    while (xcount <=10 or ycount <=10):
                        vd_x_inc /= 2
                        vd_y_inc /= 2
                        xcount, ycount, dst_gt = vd_region.geo_transform(
                            x_inc=vd_x_inc, y_inc=vd_y_inc, node='grid'
                        )

                    trans_fn, trans_fn_unc = VerticalTransform(
                        'IDW', vd_region, vd_x_inc, vd_y_inc, in_vertical_epsg, out_vertical_epsg,
                        geoid_in=src_geoid, geoid_out=dst_geoid, cache_dir=cache_dir, verbose=True
                    ).run(outfile=trans_fn)                        

            if trans_fn is not None and os.path.exists(trans_fn):
                #utils.echo_msg('using vertical tranformation grid {} from {} to {}'.format(trans_fn, in_vertical_epsg, out_vertical_epsg))
                out_src_srs = '{} +geoidgrids={}'.format(in_horizontal_crs.to_proj4(), trans_fn)

                if utils.str_or(in_vertical_epsg) == '6360':# or 'us-ft' in utils.str_or(src_vert, ''):
                    out_src_srs = out_src_srs + ' +vto_meter=0.3048006096012192'
                    trans_to_meter = True

                # if utils.str_or(out_vertical_epsg) == '6360':# or 'us-ft' in utils.str_or(dst_vert, ''):
                #     out_dst_srs = out_dst_srs + ' +vto_meter=0.3048006096012192'
                #     trans_from_meter = True
            else:
                utils.echo_error_msg(
                    'failed to generate vertical transformation grid between {} and {} for this region!'.format(
                        in_vertical_epsg, out_vertical_epsg
                    )
                )

            in_vertical_crs = pyproj.CRS.from_user_input(out_src_srs)
            src_proj4 = in_vertical_crs.to_proj4()
            dst_proj4 = out_horizontal_crs.to_proj4()
            aux_src_proj4 = in_horizontal_crs.to_proj4()
            aux_dst_proj4 = out_horizontal_crs.to_proj4()
            #utils.echo_msg('{} {}'.format(in_vertical_crs, out_horizontal_crs))
            aoi = pyproj.aoi.AreaOfInterest(region.xmin, region.ymin, region.xmax, region.ymax)
            transformer = pyproj.Transformer.from_crs(in_vertical_crs, out_horizontal_crs, always_xy=True, area_of_interest=aoi)

    return(transformer)

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
##
## generate a vertical transformation grid based on input/output vertical epsg
## also generate an associated uncertainty grid
## ==============================================
class VerticalTransform:    
    def __init__(self, mode, src_region, src_x_inc, src_y_inc, epsg_in, epsg_out,
                 geoid_in=None, geoid_out='g2018', node='pixel', verbose=True, cache_dir=None):
        self.src_region = src_region
        self.src_x_inc = utils.str2inc(src_x_inc)
        self.src_y_inc = utils.str2inc(src_y_inc)
        self.epsg_in = self._datum_by_name(str(epsg_in))
        self.epsg_out = self._datum_by_name(str(epsg_out))
        self.geoid_in = geoid_in
        self.geoid_out = geoid_out
        self.cache_dir = _vdatums_cache if cache_dir is None else cache_dir
        self.verbose = verbose
        self.xcount, self.ycount, self.gt = self.src_region.geo_transform(x_inc=self.src_x_inc, y_inc=self.src_y_inc, node='grid')
        #utils.echo_msg('transform: {} {}'.format(self.xcount, self.ycount))
        self.ref_in, self.ref_out = self._frames(self.epsg_in, self.epsg_out)
        self.node = node
        self.mode = mode
        #self.geoid_in = 'geoid09'
        
    def _frames(self, epsg_in, epsg_out):
        if epsg_in in _tidal_frames.keys():
            ref_in = 'tidal'
        elif epsg_in in _htdp_reference_frames.keys():
            ref_in = 'htdp'
        elif epsg_in in _cdn_reference_frames.keys():
            ref_in = 'cdn'
        else:
            ref_in = None
            
        if epsg_out in _tidal_frames.keys():
            ref_out = 'tidal'
        elif epsg_out in _htdp_reference_frames.keys():
            ref_out = 'htdp'
        elif epsg_out in _cdn_reference_frames.keys():
            ref_out = 'cdn'
        else:
            ref_out = None

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

    ## ==============================================
    ## tidal transformation (VDatum)
    ## ==============================================
    def _tidal_transform(self, vdatum_tidal_in, vdatum_tidal_out):
        """generate tidal transformation grid

        This will fail over land or outside of US waters...
        """

        from cudem import waffles

        v_in = fetches.VDATUM(src_region=self.src_region, datatype=vdatum_tidal_in, verbose=self.verbose)
        v_in._outdir = self.cache_dir
        v_in.run()

        if vdatum_tidal_out is not None:
            v_out = fetches.VDATUM(src_region=self.src_region, datatype=vdatum_tidal_out, verbose=self.verbose)
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
            # if vdatum_tidal_in.lower() == 'crd':
            #     cg, cv = self._cdn_transform(name='geoid', geoid='geoid09', invert=False)
            # else:
            #     cg = np.zeros((self.ycount, self.xcount))
                
            if vdatum_tidal_in != 5714 and vdatum_tidal_in != 'msl': 
                _trans_in = waffles.WaffleFactory(mod=self.mode, data=['vdatum:datatype={}'.format(vdatum_tidal_in)], src_region=self.src_region,
                                                  xinc=self.src_x_inc, yinc=self.src_y_inc, name='{}'.format(vdatum_tidal_in),
                                                  cache_dir=self.cache_dir, dst_srs='epsg:4326', node='pixel', verbose=self.verbose)._acquire_module()
                _trans_in.initialize()
                _trans_in.generate()

                utils.remove_glob('vdatum:datatype={}.inf'.format(vdatum_tidal_in))
                _trans_in_array, _trans_in_infos = gdalfun.gdal_get_array(_trans_in.fn)
                #_trans_in_array += cg
                #print(_trans_in_array.shape)
                if _trans_in_array is None:
                    _trans_in = None
                else:
                    utils.remove_glob('{}*'.format(_trans_in.fn))
            else:
                _trans_in = None
            
        if vdatum_tidal_out is None or not v_out.results:
            utils.echo_warning_msg(
                'could not locate {} in the region {}'.format(
                    vdatum_tidal_out, self.src_region
                )
            )
            _trans_out_array, _trans_out_infos = [np.zeros( (self.ycount, self.xcount) ), None]
            _trans_out = None
            #return(np.zeros( (self.ycount, self.xcount) ), None)
        else:
            if vdatum_tidal_out != 5714 and vdatum_tidal_out != 'msl':
                _trans_out = waffles.WaffleFactory(mod=self.mode, data=['vdatum:datatype={}'.format(vdatum_tidal_out)], src_region=self.src_region,
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
                return(_trans_in_array, self._datum_by_name('msl'))#vdatum_tidal_out))
            else:
                return(_trans_in_array - _trans_out_array, self._datum_by_name(vdatum_tidal_out))

    ## ==============================================
    ## "CDN" transformation (proj)
    ## ==============================================
    def _cdn_transform(self, epsg=None, name=None, geoid='g2018', invert=False):
        """create a cdn transofrmation grid"""

        epsg = 5703 if epsg == 6360 else epsg
        geoid = 'g2018' if geoid is None else geoid
        ## fetch the cdn transformation grids
        if epsg is not None:
            cdn_results = fetches.search_proj_cdn(
                self.src_region, epsg=epsg, cache_dir=self.cache_dir, verbose=self.verbose
            )
        else:
            cdn_results = fetches.search_proj_cdn(
                self.src_region, cache_dir=self.cache_dir, verbose=self.verbose
            )

        ## get the proper cdn result based on the specified geoid
        if len(cdn_results) > 0:
            for _result in cdn_results:
                if geoid in _result['name']:
                    cdn_results = [_result]
                    break
                    
        if len(cdn_results) > 0:
            for _result in cdn_results:
                src_code = _result['source_crs_code']
                dst_code = _result['target_crs_code']
                if src_code is not None:
                    src_code = int(src_code.split(':')[-1])

                if dst_code is not None:
                    dst_code = int(dst_code.split(':')[-1])
                    
                if epsg == dst_code or np.any([g in _result['name'] for g in _geoids.keys()]):
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

                            return(_tmp_array, src_code)

        utils.echo_error_msg('failed to locate transformation for {}'.format(epsg))
        return(np.zeros( (self.ycount, self.xcount) ), epsg)

    ## ==============================================
    ## HTDP transformation
    ## ==============================================
    def _htdp_transform(self, epsg_in, epsg_out):
        """create an htdp transformation grid"""

        if utils.config_check()['HTDP'] is None:
            utils.echo_error_msg('you must have HTDP installed to perform HTDP vertical transformations')            
            return(np.zeros( (self.ycount, self.xcount) ), epsg_out)
        else:
            htdp = htdpfun.HTDP(verbose=self.verbose)
            # if self.verbose:
            #     utils.echo_msg('{}: HTDP: {}->{}'.format(self.src_region, epsg_in, epsg_out))

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

    ## ==============================================
    ## Combined vertical transformation from epsg_in -> epsg_out
    ## ==============================================
    def _vertical_transform(self, epsg_in, epsg_out):
        trans_array = np.zeros( (self.ycount, self.xcount) )
        unc_array = np.zeros( (self.ycount, self.xcount) )
        if self.epsg_in == 0:
            self.geoid_in = 'geoid09'
            
        if self.geoid_in is not None:# and self.geoid_out is not None:
            if self.epsg_in == 0:
                tg, tv = self._tidal_transform(_tidal_frames[self.epsg_in]['name'], None)
            else:
                tg = np.zeros((self.ycount, self.xcount))
                
            unc_array = np.sqrt(unc_array**2 + _geoids[self.geoid_in]['uncertainty']**2)
            tmp_trans_geoid, epsg_in = self._cdn_transform(name='geoid', geoid=self.geoid_in, invert=False)
            tmp_trans_geoid += tg
        else:
            tmp_trans_geoid = np.zeros((self.ycount, self.xcount))

        #utils.echo_msg('{} {}'.format(epsg_in, epsg_out))
        while epsg_in != epsg_out and epsg_in is not None and epsg_out is not None:
            utils.echo_msg('{} --> {}'.format(epsg_in, epsg_out))
            ref_in, ref_out = self._frames(epsg_in, epsg_out)
            utils.echo_msg('{} --> {}'.format(ref_in, ref_out))
            if ref_in == 'tidal':
                if ref_out == 'tidal':
                    tmp_trans, v = self._tidal_transform(_tidal_frames[epsg_in]['name'], _tidal_frames[epsg_out]['name'])
                    epsg_in = epsg_out
                else:
                    tg, tv = self._tidal_transform(_tidal_frames[self.epsg_in]['name'], 'tss')
                    ## crd here outputs navd88 geoid09
                    #cg, cv = self._cdn_transform(name='geoid', geoid='geoid09', invert=False)
                    cg, cv = self._cdn_transform(name='geoid', geoid=self.geoid_out, invert=False)
                    tmp_trans = tg + cg
                    epsg_in = cv
            elif ref_in == 'htdp':
                if ref_out == 'htdp':
                    unc_array = np.sqrt(unc_array**2 + float(_htdp_reference_frames[epsg_in]['uncertainty'] ** 2))
                    tmp_trans, v = self._htdp_transform(epsg_in, epsg_out)
                    epsg_in = epsg_out
                else:
                    unc_array = np.sqrt(unc_array**2 + float(_htdp_reference_frames[epsg_in]['uncertainty'] ** 2))
                    unc_array = np.sqrt(unc_array**2 + float(_cdn_reference_frames[epsg_out]['uncertainty'] ** 2))
                    cg, cv = self._cdn_transform(epsg=epsg_out, invert=True)
                    unc_array = np.sqrt(unc_array**2 + float(_geoids[self.geoid_out]['uncertainty'] ** 2))
                    #gg, gv = self._cdn_transform(name='geoid', geoid=self.geoid_out, invert=True)
                    hg, v = self._htdp_transform(epsg_in, cv)
                    tmp_trans = cg + hg
                    epsg_in = epsg_out
            elif ref_in == 'cdn':
                if ref_out == 'cdn':
                    ## todo: only do geoid transforms with 5703 (navd88)
                    tmp_trans, cv = self._cdn_transform(epsg=epsg_in, invert=False)
                    #if epsg_in == 5703:
                    #    cg, _ = self._cdn_transform(name='geoid', geoid=self.geoid_out, invert=False)
                    #    tmp_trans = tmp_trans + cg
                        
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

            #print(unc_array)
            #print(trans_array.shape)
            #print(tmp_trans.shape)
            #print(tmp_trans_geoid.shape)
            #print(trans_array)
            #print(tmp_trans)
            #print(tmp_trans_geoid)
            
            trans_array = trans_array + tmp_trans + tmp_trans_geoid
            tmp_trans = None
            
        return(trans_array, unc_array)
    
    def run(self, outfile=None):
        if outfile is None:
            outfile = utils.make_temp_fn('trans_grid.tif', self.cache_dir)

        if self.verbose:
            utils.echo_msg(outfile)
        
        if os.path.exists(outfile):
            if self.verbose:
                utils.echo_warning_msg('{} exists, skipping...'.format(outfile))
                
            return(outfile)

        unc_outfile = '{}_unc.{}'.format(utils.fn_basename2(outfile), utils.fn_ext(outfile))

        if self.verbose:
            utils.echo_msg(unc_outfile)

            utils.echo_msg('{} {}'.format(self.epsg_in, self.epsg_out))
        if self.epsg_in is None or self.epsg_out is None:
            utils.echo_error_msg('failed to parse vertical input or output, check inputs')
                
            return(None)
        else:
            if self.verbose:
                utils.echo_msg('{} {}'.format(self.epsg_in, self.epsg_out))
                
            trans_array, unc_array = self._vertical_transform(self.epsg_in, self.epsg_out)
            trans_infos = gdalfun.gdal_set_infos(
                self.xcount, self.ycount, self.xcount*self.ycount, self.gt, None, gdal.GDT_Float32, -9999, 'GTiff', None, None
            )

            if outfile is not None:
                gdalfun.gdal_write(trans_array, outfile, trans_infos)
                trans_array = None
                gdalfun.gdal_write(unc_array, unc_outfile, trans_infos)
                return(outfile, unc_outfile)
            else:
                return(trans_array, unc_array, trans_infos)

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
        src_horz, src_vert = gdalfun.split_srs(gdalfun.gdal_get_srs(src_grid))
        trans_region = src_region.copy()
        
        #if src_horz is not None:
        src_region.src_srs = src_horz
        
        trans_region.src_srs = src_horz
        trans_region.warp()
        #if trans_region is None:
        #    trans_region = src_region.copy()
        #trans_region.xmin = None
        #trans_region.xmax = None
        #trans_region.ymax = None
        #trans_region.ymin = None
        #utils.echo_msg(trans_region.valid_p(check_xy = True))
        if not trans_region.valid_p(check_xy = True):
            utils.echo_warning_msg('failed to transform source region {}!'.format(src_region))
            trans_region = src_region.copy()
            
        trans_region.buffer(pct=2)
        trans_region._wgs_extremes()

        #x_inc, y_inc = trans_region.increments(src_infos['nx']/3, src_infos['ny']/3)
        
        x_inc = src_infos['geoT'][1]
        y_inc = -src_infos['geoT'][5]
        tmp_x_inc = 3/3600
        tmp_y_inc = 3/3600

        # utils.echo_msg('input region: {}'.format(src_region))
        # utils.echo_msg('input horizontal proj: {}'.format(src_horz))
        # utils.echo_msg('input vertical proj: {}'.format(src_vert))
        # utils.echo_msg('trans region: {}'.format(trans_region))
        # utils.echo_msg('input inf: {}'.format(src_infos))
        
        vt = VerticalTransform('IDW', trans_region, tmp_x_inc, tmp_y_inc, vdatum_in, vdatum_out, cache_dir=cache_dir)
        _trans_grid, _trans_grid_unc = vt.run()
        out_trans_grid = utils.make_temp_fn('_trans_grid.tif', cache_dir)
        
        if os.path.exists(out_trans_grid):
            utils.remove_glob(out_trans_grid)
        
        if _trans_grid is not None:

            #print(src_grid)
            #print(gdalfun.gdal_get_srs(src_grid))
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

            out_horz_srs = osr.SpatialReference()
            out_horz_srs.SetFromUserInput(src_horz)
            out_vert_srs = osr.SpatialReference()
            out_vert_srs.SetFromUserInput('epsg:{}'.format(vdatum_out))
            out_src_srs = osr.SpatialReference()
            out_src_srs.SetCompoundCS('Combined'.format(out_horz_srs, out_vert_srs), out_horz_srs, out_vert_srs)
            
            gdalfun.gdal_set_srs(dst_grid.replace(' ', '\ '), out_src_srs.ExportToWkt())
            
        else:
            utils.echo_error_msg('could not parse input/output vertical datums: {} -> {}; check spelling, etc'.format(vdatum_in, vdatum_out))

        #if not keep_cache:
        #    utils.remove_glob(cache_dir)
            
### End
