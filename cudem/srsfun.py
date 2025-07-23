### srsfun.py - Projection functions
##
## Copyright (c) 2010 - 2024 Regents of the University of Colorado
##
## srsfun.py is part of CUDEM
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
## Functions, etc. for common pyproj/osr usage.
##
### Code:

import os

from osgeo import osr
from osgeo import ogr
import pyproj
from pyproj import CRS

import numpy as np

from cudem import utils
from cudem import vdatums

gc = utils.config_check()
ogr.DontUseExceptions()
osr.DontUseExceptions()


###############################################################################
## OSR/WKT/proj
###############################################################################
def split_srs(srs, as_epsg = False):
    """split an SRS into the horizontal and vertical elements.

    -----------
    Parameters:
    srs (str): an srs for input into SetFromUserInput
    as_epsg (bool): try to output EPSG codes if possible

    --------
    Returns:
    list: [horz_srs, vert_srs]
    """

    vert_epsg = None
    vert_wkt = None
    if srs is None:
        return(None, None)

    if np.any(['geoid' in x for x in srs.split('+')]):
        srs = '+'.join(srs.split('+')[:-1])
    
    if srs.split(':')[0] == 'ESRI':
        esri_split = srs.split('+')
        if len(esri_split) > 1:
            vert_epsg = srs.split('+')[1]
            vert_srs = osr.SpatialReference()
            vert_srs.SetFromUserInput('EPSG:{}'.format(vert_epsg))
            vert_wkt = vert_srs.ExportToWkt()
            
        srs = esri_split[0]
        
    src_srs = osr.SpatialReference()
    src_srs.SetFromUserInput(srs)
    #try:
    srs_wkt = src_srs.ExportToWkt()
    wkt_CRS = CRS.from_wkt(srs_wkt)
    #except:
    #    return(None, None)
    
    if wkt_CRS.is_compound:
        horz = wkt_CRS.sub_crs_list[0]
        horz_epsg = horz.to_epsg()
        horz_wkt = horz.to_wkt()
        vert = wkt_CRS.sub_crs_list[1]
        vert_epsg = vert.to_epsg()
        vert_wkt = vert.to_wkt()
    else:
        horz = wkt_CRS
        horz_epsg = horz.to_epsg()
        horz_wkt = horz.to_wkt()

    if as_epsg:
        return(horz_epsg if horz_epsg is not None else horz_wkt,
               vert_epsg if vert_epsg is not None else vert_wkt)
    else:
        return(horz_wkt, vert_epsg)

    
def combine_epsgs(src_horz, src_vert, name='Combined'):
    """combine src_horz and src_vert into a CompoundCS"""
    
    if src_horz is None or src_vert is None:
        return(None)
    
    horz_srs = osr.SpatialReference()
    horz_srs.SetFromUserInput(src_horz)
    vert_srs = osr.SpatialReference()
    vert_srs.SetFromUserInput('epsg:{}'.format(src_vert))
    src_srs = osr.SpatialReference()
    src_srs.SetCompoundCS('{}'.format(name, src_horz, src_vert),
                          horz_srs, vert_srs)
    return(src_srs.ExportToWkt())


def wkt2geom(wkt):
    """transform a wkt to an ogr geometry

    -----------
    Parameters:
    wkt (wkt): a wkt geometry

    --------
    Returns:
    ogr-geom: the ogr geometry
    """
    
    return(ogr.CreateGeometryFromWkt(wkt))


def osr_srs(src_srs):
    try:
        srs = osr.SpatialReference()
        srs.SetFromUserInput(src_srs)
        return(srs)
    except:
        return(None)

    
def osr_wkt(src_srs, esri=False):
    """convert a src_srs to wkt"""
    
    try:
        sr = osr.SpatialReference()
        sr.SetFromUserInput(src_srs)
        if esri:
            sr.MorphToESRI()
            
        return(sr.ExportToWkt())
    except:
        return(None)

    
def osr_prj_file(dst_fn, src_srs):
    """generate a .prj file given a src_srs

    -----------
    Parameters:
    dst_fn (str): the output filename
    src_srs (srs): the output srs
    """
    
    with open(dst_fn, 'w') as out:
        out.write(osr_wkt(src_srs, True))

    if os.path.exists(dst_fn):
        return(0)
    else:
        return(-1)

    
def srs_get_cstype(in_srs):
    src_srs = osr.SpatialReference()
    src_srs.SetFromUserInput(in_srs)
    
    ## HORZ
    if src_srs.IsGeographic() == 1:
        cstype = 'GEOGCS'
    else:
        cstype = 'PROJCS'

    src_srs = None

    return(cstype)


def epsg_from_input(in_srs):
    """get the epsg(s) from srs suitable as input to 
    SetFromUserInput

    -----------
    Parameters:
    in_srs (str): an srs as a string

    -----------
    Returns:
    list: [horz_epsg, vert_epsg]
    """

    src_vert = None
    if in_srs.split('+')[-1].startswith('geoid'):
        in_srs = '+'.join(in_srs.split('+')[:-1])
        
    if in_srs.split(':')[0] == 'ESRI':
        esri_split = in_srs.split('+')
        if len(esri_split) > 1:
            vert_epsg = in_srs.split('+')[1]
            vert_srs = osr.SpatialReference()
            vert_srs.SetFromUserInput('EPSG:{}'.format(vert_epsg))
            vert_wkt = vert_srs.ExportToWkt()

            src_vert = vert_epsg
            
        in_srs = esri_split[0]

    src_srs = osr.SpatialReference()
    src_srs.SetFromUserInput(in_srs)
    ## HORZ
    if src_srs.IsGeographic() == 1:
        cstype = 'GEOGCS'
    else:
        cstype = 'PROJCS'

    src_srs.AutoIdentifyEPSG()
    an = src_srs.GetAuthorityName(cstype)
    src_horz = src_srs.GetAuthorityCode(cstype)

    if src_horz is None:
        src_horz = src_srs.ExportToProj4()
        if src_horz == '':
            src_horz = src_srs.ExportToWkt()
    else:
        src_horz = '{}:{}'.format(an, src_horz)
        
    ## VERT
    if src_srs.IsVertical() == 1:
        csvtype = 'VERT_CS'
        src_vert = src_srs.GetAuthorityCode(csvtype)
    #else:
    #    src_vert = None

    return(src_horz, src_vert)


def osr_parse_srs(src_srs, return_vertcs=True):
    """parse an OSR SRS object and return a proj string"""
    
    if src_srs is not None:
        if src_srs.IsLocal() == 1:
            return(src_srs.ExportToWkt())
        if src_srs.IsGeographic() == 1:
            cstype = 'GEOGCS'
        else:
            cstype = 'PROJCS'
            
        src_srs.AutoIdentifyEPSG()
        an = src_srs.GetAuthorityName(cstype)
        ac = src_srs.GetAuthorityCode(cstype)

        #if return_vertcs:
        if src_srs.IsVertical() == 1:
            csvtype = 'VERT_CS'
            vn = src_srs.GetAuthorityName(csvtype)
            vc = src_srs.GetAuthorityCode(csvtype)
        else:
            csvtype = vc = vn = None

        if an is not None and ac is not None:
            if vn is not None and vc is not None:
                return('{}:{}+{}'.format(an, ac, vc))
            else:
                return('{}:{}'.format(an, ac))
        else:
            dst_srs = src_srs.ExportToProj4()
            if dst_srs:
                    return(dst_srs)
            else:
                return(None)
    else:
        return(None)

    
def parse_srs(src_srs=None, dst_srs=None):
    transform = {}
    if src_srs is not None and dst_srs is not None:
        want_vertical = True
        src_geoid = None
        dst_geoid = 'g2018' # dst_geoid is g2018 by default.
        in_vertical_epsg = in_vertical_crs = None
        out_vertical_epsg = out_vertical_crs = None

        ## parse out the source geoid, which is not standard in proj,
        ## reset `self.src_srs` without it.
        tmp_src_srs = src_srs.split('+geoid:')
        src_srs_ = tmp_src_srs[0]
        src_srs = src_srs_
        if len(tmp_src_srs) > 1:
            src_geoid = tmp_src_srs[1]

        ## parse out the destination geoid, which is not standard in proj
        ## reset `self.dst_srs` without it.
        tmp_dst_srs = dst_srs.split('+geoid:')
        dst_srs_ = tmp_dst_srs[0]
        dst_srs = dst_srs_
        if len(tmp_dst_srs) > 1:
            dst_geoid = tmp_dst_srs[1]

        ## check if this is an ESRI epsg code
        is_esri = False
        in_vertical_epsg_esri = None
        if 'ESRI' in src_srs.upper():
            is_esri = True
            srs_split = src_srs.split('+')
            src_srs = srs_split[0]
            if len(srs_split) > 1:
                in_vertical_epsg_esri = srs_split[1]

        ## check if the vertical epsg is tidal, fix the epsg code if so
        if utils.int_or(src_srs.split('+')[-1]) in vdatums._tidal_frames.keys():
            src_srs = '{}+{}'.format(
                src_srs.split('+')[0],
                vdatums._tidal_frames[utils.int_or(src_srs.split('+')[-1])]['epsg']
            )

        ## set the proj crs from the src and dst srs input
        #try:
        in_crs = pyproj.CRS.from_user_input(src_srs)
        out_crs = pyproj.CRS.from_user_input(dst_srs)
        #except:
        #    utils.echo_error_msg(src_srs)
        #    utils.echo_error_msg(src_srs)

        ## if the crs has vertical (compound), parse out the vertical crs
        ## and set the horizontal and vertical crs
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
            out_vertical_epsg=None

        ## check if esri vertical
        if (in_vertical_epsg_esri is not None and is_esri):
            in_vertical_epsg = in_vertical_epsg_esri
            if out_vertical_epsg is not None:
                want_vertical = True

        ## make sure the input and output vertical epsg is different
        if want_vertical:
            if (in_vertical_epsg == out_vertical_epsg) and src_geoid is None:
                want_vertical = False

        ## set `self.transform` with the parsed srs info
        transform['src_horz_crs'] = in_horizontal_crs
        transform['dst_horz_crs'] = out_horizontal_crs
        transform['src_vert_crs'] = in_vertical_crs
        transform['dst_vert_crs'] = out_vertical_crs
        transform['src_vert_epsg'] = in_vertical_epsg
        transform['dst_vert_epsg'] = out_vertical_epsg
        transform['src_geoid'] = src_geoid
        transform['dst_geoid'] = dst_geoid
        transform['want_vertical'] = want_vertical

        return(transform)

    
def set_vertical_transform(transform, region=None, infos=None,
                           cache_dir='./', verbose=True):
    ## set the region for the vdatum transformation grid.
    ## this is either from the input `self.region` or from the input
    ## data's bounding box. Transform it to WGS84.
    if region is None:
        vd_region = regions.Region().from_list(infos.minmax)
        vd_region.src_srs = transform['src_horz_crs'].to_proj4()
    else:
        vd_region = region.copy()
        vd_region.src_srs = transform['dst_horz_crs'].to_proj4()

    vd_region.zmin = None
    vd_region.zmax = None
    vd_region.warp('epsg:4326')
    vd_region.buffer(pct=10)
    if not vd_region.valid_p():
        utils.echo_warning_msg('failed to generate transformation')
        return
    # else:
    #     utils.echo_msg('generating vertical transformation to region {}'.format(vd_region))

    ## set `self.transform.trans_fn`, which is the transformation grid
    transform['trans_fn'] = os.path.join(
        cache_dir, '_vdatum_trans_{}_{}_{}.tif'.format(
            transform['src_vert_epsg'],
            transform['dst_vert_epsg'],
            vd_region.format('fn')
        )
    )

    ## if the transformation grid already exists, skip making a new one,
    ## otherwise, make the new one here with `vdatums.VerticalTransform()`
    if not os.path.exists(transform['trans_fn']):
        with tqdm(
                desc='generating vertical transformation grid {} from {} to {}'.format(
                    transform['trans_fn'],
                    transform['src_vert_epsg'],
                    transform['dst_vert_epsg']
                ),
                leave=verbose
        ) as pbar:
            ## set the vertical transformation grid to be 3 arc-seconds. This
            ## is pretty arbitrary, maybe it's too small...
            vd_x_inc = vd_y_inc = utils.str2inc('3s')
            xcount, ycount, dst_gt = vd_region.geo_transform(
                x_inc=vd_x_inc, y_inc=vd_y_inc, node='grid'
            )

            ## if the input region is so small it creates a tiny grid,
            ## keep increasing the increments until we are at least to
            ## a 10x10 grid.
            while (xcount <=10 or ycount <=10):
                vd_x_inc /= 2
                vd_y_inc /= 2
                xcount, ycount, dst_gt = vd_region.geo_transform(
                    x_inc=vd_x_inc, y_inc=vd_y_inc, node='grid'
                )

            ## run `vdatums.VerticalTransform()`, grid using `nearest`
            transform['trans_fn'], transform['trans_fn_unc'] = vdatums.VerticalTransform(
                'IDW',
                vd_region,
                vd_x_inc,
                vd_y_inc,
                transform['src_vert_epsg'],
                transform['dst_vert_epsg'],
                geoid_in=transform['src_geoid'],
                geoid_out=transform['dst_geoid'],
                cache_dir=cache_dir,
                verbose=False
            ).run(outfile=transform['trans_fn'])

    ## set the pyproj.Transformer for both horz+vert and vert only
    if transform['trans_fn'] is not None and os.path.exists(transform['trans_fn']):
        transform['pipeline'] = '+proj=pipeline +step {} +inv +step +proj=vgridshift +grids={} +inv +step {}'.format(
            transform['src_horz_crs'].to_proj4(),
            os.path.abspath(transform['trans_fn']),
            transform['dst_horz_crs'].to_proj4()
        )
        self.transform['vert_transformer'] = pyproj.Transformer.from_pipeline(
            '+proj=pipeline +step +proj=vgridshift +grids={} +inv'.format(os.path.abspath(transform['trans_fn']))
            )
    else:
        utils.echo_error_msg(
            ('failed to generate vertical transformation grid '
             f'between {transform["src_vert_epsg"]} '
             f'and {os.path.abspath(transform["dst_vert_epsg"])} for this region!')
        )

        
def set_transform(src_srs=None, dst_srs=None, region=None, infos=None):
    """Set the pyproj horizontal and vertical transformations for the dataset"""

    #in_horizontal_crs, out_horizontal_crs, in_vertical_crs, out_vertical_crs, src_geoid, dst_geoid, want_vertical = self.parse_srs()
    transform = parse_srs(src_srs=src_srs, dst_srs=dst_srs)
    if transform['src_horz_crs'] is not None \
       and transform['dst_horz_crs'] is not None:        
        ## horizontal Transformation
        transform['horz_pipeline'] = '+proj=pipeline +step {} +inv +step {}'.format(
            transform['src_horz_crs'].to_proj4(), transform['dst_horz_crs'].to_proj4()
        )
        if region is not None:
            transform['trans_region'] = region.copy()
            transform['trans_region'].src_srs = transform['dst_horz_crs'].to_proj4()
            transform['trans_region'].warp(transform['src_horz_crs'].to_proj4())
        else:
            transform['trans_region'] = regions.Region().from_list(infos.minmax)
            transform['trans_region'].src_srs = infos.src_srs
            transform['trans_region'].warp(transform['dst_horz_crs'].to_proj4())

        ## vertical Transformation
        if transform['want_vertical']:
            set_vertical_transform()
        else:
            transform['pipeline'] = transform['horz_pipeline']

        try:
            transform['transformer'] = pyproj.Transformer.from_pipeline(transform['pipeline'])
        except Exception as e:
            utils.echo_warning_msg('could not set transformation in: {}, out: {}, {}'.format(
                transform['src_horz_crs'].name, transform['dst_horz_crs'].name, e
            ))
            return

        ## dataset region
        if region is not None and region.valid_p():
            data_region = region.copy() \
                if transform['trans_region'] is None \
                   else transform['trans_region'].copy()
            inf_region = regions.Region().from_list(infos.minmax)
            data_region = regions.regions_reduce(data_region, inf_region)
            data_region.src_srs = infos.src_srs

            if not data_region.valid_p():
                data_region = self.region.copy() \
                    if transform['trans_region'] is None \
                       else transform['trans_region'].copy()
        else:
            data_region = regions.Region().from_list(infos.minmax)
            data_region.src_srs = infos.src_srs

    return(transform)

    
### End
