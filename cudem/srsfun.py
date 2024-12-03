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
## Functions, etc. for common pyproj/osr usage.
##
### Code:

import os

from osgeo import osr
from osgeo import ogr
from pyproj import CRS

import numpy as np

from cudem import utils

gc = utils.config_check()
ogr.DontUseExceptions()
osr.DontUseExceptions()

## OSR/WKT/proj
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
        return(horz_epsg if horz_epsg is not None else horz_wkt, vert_epsg if vert_epsg is not None else vert_wkt)
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
    src_srs.SetCompoundCS('{}'.format(name, src_horz, src_vert), horz_srs, vert_srs)
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
    """get the epsg(s) from srs suitable as input to SetFromUserInput

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

def osr_parse_srs(src_srs, return_vertcs = True):
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

### End
