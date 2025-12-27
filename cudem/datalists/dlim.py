### dlim.py - DataLists IMproved
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
##
## dlim.py is part of CUDEM
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
## A datalist is similar to an MBSystem datalist; it is a space-delineated file
## containing the following columns:
## data-path data-format data-weight data-uncertainty data-name data-source
## data-date data-resolution data-type data-horz data-vert data-url
## Minimally, data-path (column 1) is all that is needed.
##
## An associated inf and geojson file will be gerenated for each datalist
## only an associated inf file will be genereated for individual datasets
##
## Parse various dataset types by region/increments and yield data as xyz or
## array recursive data-structures which point to datasets (datalist, zip,
## fetches, etc) are negative format numbers, e.g. -1 for datalist
##
## supported datasets include: xyz, gdal, ogr, las/laz (laspy), mbs (MBSystem),
## fetches (see cudem.fetches)
##
## Initialize a datalist/dataset using init_data(list-of-datalist-entries)
## where the list of datalist entries can be any of the supported dataset
## formats. init_data will combine all the input datafiles into an internal
## scratch datalist and process that.
##
## If region, x_inc, and y_inc are set, all processing will go through
## Dataset.stacks() where the data will be combined
## either using the 'supercede', 'weighted-mean', 'min' or 'max' method.
## Dataset.stacks will output a multi-banded gdal file with the following
## bands: 1: z, 2: count, 3: weight, 4: uncerainty, 5: source-uncertainty
##
## If want_mask is set, stacks() will also generate a multi-band gdal raster
## file where each mask band contains the mask (0/1) of a specific dataset/datalist,
## depending on the input data. For a datalist, each band will contain a mask for
## each of the top-level datasets. If want_sm is also set, the multi-band mask
## will be processed to an OGR supported vector, with a feature for each band and
## the metadata items (cols > 4 in the datalist entry) as fields.
##
## Transform data between horizontal/vertical projections/datums by setting src_srs
## and dst_srs as 'EPSG:<>'
## if src_srs is not set, but dst_srs is, dlim will attempt to obtain the source srs
## from the data file itself or its respective inf file; otherwise, it will be
## assumed the source data file is in the same srs as dst_srs
##
## A dataset module should have at least an `__init__` and a `yield_points` method.
## `yield_points` should yield a numpy rec array with at least
## 'x', 'y' and 'z' fields, and optionally 'w' and 'u' fields.
## Other methods that can be re-written include `parse` which yields dlim
## dataset module(s).
##
### Examples:
##
## my_data_list = ['nos_data.xy', 'my_grid.tif', 'gmrt', 'other_data.datalist'] 
## my_processed_datalist = init_data(my_data_list) # initialize the data
## my_processed_datalist.dump_xyz() # dump all the xyz data
##
## from cudem import regions
## my_region = regions.Region().from_list([-160, -150, 30, 40])
## my_data_list = ['nos_data.xy', 'my_grid.tif', 'gmrt', 'other_data.datalist'] 
## my_processed_datalist = init_data(my_data_list, src_region=my_region, x_inc=1s, y_inc=1s)
## my_stack = my_processed_datalist.stacks() # stack the data to the region/increments
## my_processed_datalist.archive_xyz() # archive the data to the region/increments
## my_mask = my_processed_datalist._mask() # mask the data to the region/increments
##
### Code:

import os
import sys
import re
import copy
import json
import glob
import math
import traceback
import warnings
import numpy as np
import pyproj
from osgeo import gdal
from osgeo import ogr
import h5py as h5
import netCDF4 as nc

from cudem import utils
from cudem import regions
from cudem import xyzfun
from cudem import gdalfun
from cudem import factory
from cudem import vdatums
from cudem import fetches
from cudem import grits
from cudem import srsfun
from cudem import pointz
from . import __version__
from . import inf

## Config info and setup
gc = utils.config_check()
gdal.DontUseExceptions()
ogr.DontUseExceptions()
gdal.SetConfigOption(
    'CPL_LOG', 'NUL' if gc['platform'] == 'win32' else '/dev/null'
) 

## ==============================================
## dataset masking and spatial metadata functions
## ==============================================
def scan_mask_bands(
        src_ds, skip_band='Full Data Mask', mask_level=0, verbose=True
):
    mask_level = utils.int_or(mask_level, 0)
    src_infos = gdalfun.gdal_infos(src_ds)
    band_infos = {}    
    with utils.ccp(
            desc='scanning mask bands.',
            total=src_ds.RasterCount,
            leave=verbose
    ) as pbar:         
        for b in range(1, src_ds.RasterCount+1):
            pbar.update()
            this_band = src_ds.GetRasterBand(b)
            this_band_md = this_band.GetMetadata()
            this_band_name = this_band.GetDescription()
            if this_band_name == skip_band:
                continue

            if len(this_band_name.split('/')) > 1:
                this_band_name = this_band_name.split('/')[1]
                
            if not this_band_name in band_infos.keys():
                band_infos[this_band_name] = {'bands': [b],
                                              'metadata': this_band_md}
            else:
                band_infos[this_band_name]['bands'].append(b)
                
            band_md = band_infos[this_band_name]['metadata']
            for k in this_band_md.keys():
                if k not in band_md.keys() \
                   or band_md[k] is None \
                   or band_md[k] == 'None':
                    if this_band_md[k] != 'None':
                        band_md[k] = this_band_md[k]
                            
                if k == 'name':
                    band_md[k] = this_band_name                
                elif k == 'date' \
                     or k == 'weight' \
                     or k == 'uncertainty':
                    band_md[k] = utils.num_strings_to_range(band_md[k], this_band_md[k])
                elif band_md[k] not in this_band_md[k].split('/') \
                     and band_md[k] != this_band_md[k]:
                    if str(this_band_md[k]) != 'None':
                        band_md[k] = ','.join([band_md[k], this_band_md[k]])
                    
            band_infos[this_band_name]['metadata'] = band_md
            
    return(band_infos)


def ogr_mask_footprints(
        src_ds, ogr_format='GPKG', dst_srs='epsg:4326',
        mask_level=0, verbose=True
):
    
    src_infos = scan_mask_bands(
        src_ds, mask_level=mask_level, verbose=verbose
    )
    ## initialize output vector
    dst_ogr_fn = (f'{utils.fn_basename2(src_ds.GetDescription())}_sm'
                  '.{gdalfun.ogr_fext(ogr_format)}')
    driver = ogr.GetDriverByName(ogr_format)
    ds = driver.CreateDataSource(dst_ogr_fn)
    srs = srsfun.osr_srs(src_ds.GetProjectionRef())
    if ds is not None: 
        layer = ds.CreateLayer('footprint', srs, ogr.wkbMultiPolygon)
    else:
        layer = None

    #layer.CreateField(ogr.FieldDefn('location', ogr.OFTString))
    field_names = [field.name for field in layer.schema]
    for key in src_infos.keys():
        for md in src_infos[key]['metadata'].keys():
            #if md[:10] not in field_names and md not in field_names:
            if md not in field_names:
                layer.CreateField(ogr.FieldDefn(md, ogr.OFTString))
                field_names = [field.name for field in layer.schema]

    [layer.SetFeature(feature) for feature in layer]
    ds = layer = None
    ## generate a footprint of each group of bands from src_infos
    with utils.ccp(
            desc='generating mask footprints',
            total=len(src_infos.keys()),
            leave=verbose
    ) as pbar:
        for i, key in enumerate(src_infos.keys()):
            footprint_cmd = 'gdal_footprint {} {} -combine_bands union -max_points unlimited {} -no_location'.format(
                src_ds.GetDescription(),
                dst_ogr_fn,
                ' '.join(['-b {}'.format(x) for x in  src_infos[key]['bands']]),
            )
            utils.run_cmd(footprint_cmd, verbose=False)
            pbar.update()
            
    with utils.ccp(
            desc='setting mask metadata',
            total=len(src_infos.keys()),
            leave=verbose
    ) as pbar:
        for i, key in enumerate(src_infos.keys()):
            ## update db using ogrinfo
            key_val = ', '.join(
                ["'{}' = '{}'".format(
                    x, src_infos[key]['metadata'][x]
                ) for x in src_infos[key]['metadata'].keys()]
            )
            sql = "UPDATE {name} SET {key_val} WHERE rowid = {f}".format(
                name='footprint',
                key_val=key_val,
                value=src_infos[key]['metadata'][md],
                f=i+1
            )
            utils.run_cmd(
                'ogrinfo {} -dialect SQLite -sql "{}"'.format(dst_ogr_fn, sql),
                verbose=False
            )
            pbar.update()
            
    return(glob.glob('{}.*'.format(utils.fn_basename2(dst_ogr_fn))), ogr_format)


def merge_mask_bands_by_name(
        src_ds, verbose=True, skip_band='Full Data Mask', mask_level=0
):
    """merge data mask raster bands based on the top-level datalist name"""

    band_infos = scan_mask_bands(
        src_ds, mask_level=mask_level, verbose=verbose
    )
    mask_level = utils.int_or(mask_level, 0)
    src_infos = gdalfun.gdal_infos(src_ds)
    ## create new mem gdal ds to hold new raster
    gdt = gdal.GDT_Byte
    driver = gdal.GetDriverByName('MEM')
    m_ds = driver.Create(
        utils.make_temp_fn('tmp'),
        src_infos['nx'],
        src_infos['ny'],
        len(band_infos.keys()),
        gdt
    )
    m_ds.SetGeoTransform(src_infos['geoT'])
    m_ds.SetProjection(gdalfun.osr_wkt(src_infos['proj']))
    m_ds.SetDescription(src_ds.GetDescription())

    with utils.ccp(
            desc='generating merged mask dataset.',
            total=len(band_infos.keys()),
            leave=verbose
    ) as pbar:    
        for i, m in enumerate(band_infos.keys()):
            m_band = m_ds.GetRasterBand(i+1)
            m_band.SetDescription(m)
            m_band.SetNoDataValue(0)
            m_band_array = m_band.ReadAsArray()

            try:
                m_band.SetMetadata(band_infos[m]['metadata'])
            except:
                try:
                    band_md = m_band.GetMetadata()
                    for key in band_infos[m]['metadata'].keys():
                        if band_md[key] is None:
                            del band_md[key]

                    m_band.SetMetadata(band_md)
                except Exception as e:
                    utils.echo_error_msg(
                        'could not set band metadata: {}; {}'.format(
                            band_md, e
                        )
                    )

            with utils.ccp(
                    desc='merging {} mask bands.'.format(m),
                    total=len(band_infos[m]['bands']),
                    leave=verbose
            ) as merge_pbar:                        
                for b in band_infos[m]['bands']:
                    this_band = src_ds.GetRasterBand(b)
                    this_band_array = this_band.ReadAsArray()
                    m_band_array[this_band_array > 0] = 1
                    this_band = this_band_array = None
                    merge_pbar.update()
                    
            m_band.WriteArray(m_band_array)
            m_ds.FlushCache()

            pbar.update()

    return(m_ds)


def polygonize_mask_multibands(
        src_ds,
        output = None,
        ogr_format = 'GPKG',
        mask_level = 0,
        verbose = True
):
    """polygonze a multi-band mask raster. 

    -----------
    Parameters:
    src_ds (GDALDataset): the source multi-band raster as a gdal 
                          dataset object
    dst_srs (str): the output srs
    ogr_format (str): the output OGR format
    verbose (bool): be verbose

    -----------
    Returns:
    tuple: (dst_layer, ogr_format)
    """

    mask_level = utils.int_or(mask_level, 0)
    src_ds = merge_mask_bands_by_name(
        src_ds, mask_level=mask_level, verbose=verbose
    )
    dst_layer = '{}_sm'.format(
        utils.fn_basename2(src_ds.GetDescription()) \
        if output is None \
        else os.path.basename(output)
    )
    dst_vector = os.path.join(
        os.path.dirname(output) if output is not None else '',
        f'{dst_layer}.{gdalfun.ogr_fext(ogr_format)}'
    )
    #dst_vector = dst_layer + '.{}'.format(gdalfun.ogr_fext(ogr_format))
    #utils.remove_glob('{}.*'.format(dst_layer))
    utils.remove_glob('{}.*'.format(utils.fn_basename2(dst_vector)))
    srs = srsfun.osr_srs(src_ds.GetProjectionRef())
    driver = ogr.GetDriverByName(ogr_format)
    ds = driver.CreateDataSource(dst_vector)
    if ds is not None: 
        layer = ds.CreateLayer(
            dst_layer, srs, ogr.wkbMultiPolygon
        )
        [layer.SetFeature(feature) for feature in layer]
    else:
        layer = None

    layer.CreateField(ogr.FieldDefn('DN', ogr.OFTInteger))
    layer.StartTransaction()
    defn = None

    #field_names_to_delete = ['DN', 'Name', 'Uncertainty']
    
    with utils.ccp(
            desc='polygonizing mask bands.',
            total=src_ds.RasterCount,
            leave=verbose
    ) as pbar:
        for b in range(1, src_ds.RasterCount+1):
            pbar.update()
            this_band = src_ds.GetRasterBand(b)
            this_band_md = this_band.GetMetadata()
            this_band_name = this_band.GetDescription()
            b_infos = gdalfun.gdal_infos(src_ds, scan=True, band=b)
            field_names = [field.name for field in layer.schema]
            this_band_md = {k.title():v for k,v in this_band_md.items()}            
            for k in this_band_md.keys():
                if k not in field_names:
                    layer.CreateField(ogr.FieldDefn(k, ogr.OFTString))

            if 'Title' not in this_band_md.keys():
                if 'Title' not in field_names:
                    layer.CreateField(ogr.FieldDefn('Title', ogr.OFTString))

            if b_infos['zr'][1] == 1:
                tmp_ds = ogr.GetDriverByName('Memory').CreateDataSource(
                    '{}_poly'.format(this_band_name)
                )
                if tmp_ds is not None:
                    tmp_layer = tmp_ds.CreateLayer(
                        '{}_poly'.format(this_band_name), srs, ogr.wkbMultiPolygon
                    )
                    tmp_layer.CreateField(ogr.FieldDefn('DN', ogr.OFTInteger))
                    tmp_name = str(this_band_name)                                    
                    for k in this_band_md.keys():
                        tmp_layer.CreateField(ogr.FieldDefn(k, ogr.OFTString))

                    if 'Title' not in this_band_md.keys():
                        tmp_layer.CreateField(ogr.FieldDefn('Title', ogr.OFTString))

                    ## fill tmp_layer with all associated bands in band_infos...
                    status = gdal.Polygonize(
                        this_band,
                        None,
                        tmp_layer,
                        tmp_layer.GetLayerDefn().GetFieldIndex('DN'),
                        [],
                        #callback = gdal.TermProgress if verbose else None
                        callback = None
                    )

                    if len(tmp_layer) > 0:
                        if defn is None:
                            defn = layer.GetLayerDefn()

                        out_feat = gdalfun.ogr_mask_union(tmp_layer, 'DN', defn)
                        with utils.ccp(
                                desc='creating feature {}...'.format(this_band_name),
                                total=len(this_band_md.keys()),
                                leave=verbose
                        ) as feat_pbar:
                            for k in this_band_md.keys():
                                feat_pbar.update()
                                out_feat.SetField(k, this_band_md[k])

                            if 'Title' not in this_band_md.keys():
                                out_feat.SetField('Title', tmp_name)

                            out_feat.SetField('DN', b)
                            status = layer.CreateFeature(out_feat)
                            #layer.SetFeature(out_feat)
                
                # if verbose:
                #     utils.echo_msg('polygonized {}'.format(this_band_name))
                tmp_ds = tmp_layer = None
                
    layer.CommitTransaction()
    ds = src_ds = None
    
    ds = ogr.Open(dst_vector, 1)
    field_names_to_delete = ['DN', 'Name', 'Uncertainty']
    for fnd in field_names_to_delete:
        ds.ExecuteSQL(f'ALTER TABLE {dst_layer} DROP COLUMN {fnd}')

    ds = None

    return(dst_layer, ogr_format)


def polygonize_mask_multibands2(
        src_ds,
        output = None,
        ogr_format = 'GPKG',
        mask_level = 0,
        verbose = True
):
    """polygonze a multi-band mask raster. 

    -----------
    Parameters:
    src_ds (GDALDataset): the source multi-band raster as a gdal dataset object
    dst_srs (str): the output srs
    ogr_format (str): the output OGR format
    verbose (bool): be verbose

    -----------
    Returns:
    tuple: (dst_layer, ogr_format)
    """

    mask_level = utils.int_or(mask_level, 0)
    band_infos = scan_mask_bands(src_ds, mask_level=mask_level, verbose=verbose)
    dst_layer = '{}_sm'.format(
        utils.fn_basename2(src_ds.GetDescription()) if output is None else output
    )
    dst_vector = dst_layer + '.{}'.format(gdalfun.ogr_fext(ogr_format))
    utils.remove_glob('{}.*'.format(dst_layer))
    srs = srsfun.osr_srs(src_ds.GetProjectionRef())
    driver = ogr.GetDriverByName(ogr_format)
    ds = driver.CreateDataSource(dst_vector)
    if ds is not None: 
        layer = ds.CreateLayer(
            'footprints', srs, ogr.wkbMultiPolygon
        )
        [layer.SetFeature(feature) for feature in layer]
    else:
        layer = None

    layer.CreateField(ogr.FieldDefn('DN', ogr.OFTInteger))
    #layer.StartTransaction()
    defn = None
    
    with utils.ccp(
            desc='polygonizing mask bands.',
            total=len(band_infos.keys()),
            leave=verbose
    ) as pbar:
        for i, m in enumerate(band_infos.keys()):
            field_names = [field.name for field in layer.schema]
            this_md = {k.title():v for k,v in band_infos[m]['metadata'].items()}            
            for k in this_md.keys():
                if k not in field_names:
                    layer.CreateField(ogr.FieldDefn(k, ogr.OFTString))

            if 'Title' not in this_md.keys():
                if 'Title' not in field_names:
                    layer.CreateField(ogr.FieldDefn('Title', ogr.OFTString))

            tmp_ds = ogr.GetDriverByName('Memory').CreateDataSource(
                '{}_poly'.format(m)
            )
            if tmp_ds is not None:
                tmp_layer = tmp_ds.CreateLayer(
                    '{}_poly'.format(m), srs, ogr.wkbMultiPolygon
                )
                tmp_layer.CreateField(ogr.FieldDefn('DN', ogr.OFTInteger))
                tmp_name = str(m)                                    
                for k in this_md.keys():
                    tmp_layer.CreateField(ogr.FieldDefn(k, ogr.OFTString))

                if 'Title' not in this_md.keys():
                    tmp_layer.CreateField(ogr.FieldDefn('Title', ogr.OFTString))

                with utils.ccp(
                        desc='polygonizing {} mask bands from {}.'.format(
                            len(band_infos[m]['bands']), m,
                        ),
                        total=len(band_infos[m]['bands']),
                        leave=verbose
                ) as poly_pbar:                                            
                    for b in band_infos[m]['bands']:
                        this_band = src_ds.GetRasterBand(b)
                        status = gdal.Polygonize(
                            this_band,
                            None,
                            tmp_layer,
                            tmp_layer.GetLayerDefn().GetFieldIndex('DN'),
                            [],
                            #callback = gdal.TermProgress if verbose else None
                            callback = None
                        )
                        poly_pbar.update()

                if len(tmp_layer) > 0:
                    if defn is None:
                        defn = layer.GetLayerDefn()

                    out_feat = gdalfun.ogr_mask_union(tmp_layer, 'DN', defn)
                    with utils.ccp(
                            desc='creating feature {}...'.format(m),
                            total=len(this_md.keys()),
                            leave=verbose
                    ) as feat_pbar:
                        for k in this_md.keys():
                            feat_pbar.update()
                            out_feat.SetField(k, this_md[k])

                        if 'Title' not in this_md.keys():
                            out_feat.SetField('Title', tmp_name)

                        out_feat.SetField('DN', b)
                        status = layer.CreateFeature(out_feat)
                        layer.SetFeature(out_feat)
                                    
            tmp_ds = tmp_layer = None
            pbar.update()    
    #layer.CommitTransaction()            
    ds = src_ds = None
    return(dst_layer, ogr_format)


## ==============================================
## Datalist convenience functions
## data_list is a list of dlim supported datasets
## ==============================================
def get_factory_exts():
    fmts = []
    for key in DatasetFactory()._modules.keys():
        if key != '_factory' and int(key) > 0:
            for f in DatasetFactory()._modules[key]['fmts']:
                if f not in fmts:
                    fmts.append(f)
                    
    return(fmts)

def make_datalist(data_list, want_weight, want_uncertainty, region,
                  src_srs, dst_srs, x_inc, y_inc, sample_alg,
                  verbose):
    """Make a datalist object from a list of supported datasets"""

    ## Make the datalist object
    xdl = Datalist(
        fn='scratch',
        data_format=-1,
        weight=None if not want_weight else 1,
        uncertainty=None if not want_uncertainty else 0,
        src_region=region,
        verbose=verbose,
        parent=None,
        src_srs=src_srs,
        dst_srs=dst_srs,
        x_inc=x_inc,
        y_inc=y_inc,
        sample_alg=sample_alg
    )

    ## add the datasets from data_list to the datalist object
    xdl.data_entries = [DatasetFactory(
        fn=" ".join(['-' if x == "" else x for x in dl.split(",")]),
        weight=None if not want_weight else 1,
        uncertainty=None if not want_uncertainty else 0,
        src_region=region,
        verbose=verbose,
        src_srs=src_srs,
        dst_srs=dst_srs,
        x_inc=x_inc,
        y_inc=y_inc,
        sample_alg=sample_alg,
        parent=xdl,
        cache_dir=self.cache_dir
    )._acquire_module() for dl in data_list]
    return(xdl)


def write_datalist(data_list, outname=None):
    """write a datalist file from a datalist object"""
    
    if outname is None:
        outname = '{}_{}'.format(self.metadata['name'], utils.this_year())
    
    if os.path.exists('{}.datalist'.format(outname)):
        utils.remove_glob('{}.datalist*'.format(outname))
        
    with open('{}.datalist'.format(outname), 'w') as tmp_dl:
        [tmp_dl.write('{}\n'.format(x.format_entry())) for x in data_list]

    return('{}.datalist'.format(outname))


def h5_get_datasets_from_grp(grp, count = 0, out_keys = []):
    if isinstance(grp, h5.Dataset):
        count += 1
        out_keys.append(grp.name)

    elif isinstance(grp, h5.Group):
        for key in grp.keys():
            if isinstance(grp[key], h5.Dataset):
                count += 1
                out_keys.append(grp[key].name)
            elif isinstance(grp[key], h5.Group):
                for sub_key in grp[key].keys():
                    count, out_keys = h5_get_datasets_from_grp(
                        grp[key][sub_key], count, out_keys
                    )

    return(count, out_keys)


# ## initialize a list of datasets into a dataset object
# def init_data_(data_list,
#               region=None,
#               src_srs=None,
#               dst_srs=None,
#               src_geoid=None,
#               dst_geoid='g2018',
#               xy_inc=(None, None),
#               sample_alg='auto',
#               want_weight=False,
#               want_uncertainty=False,
#               want_verbose=True,
#               want_mask=False,
#               want_sm=False,
#               invert_region=False,
#               cache_dir=None,
#               dump_precision=4,
#               pnt_fltrs=None,
#               stack_fltrs=None,
#               stack_node=True,
#               stack_mode='mean',
#               upper_limit=None,
#               lower_limit=None,
#               mask=None):
#     """initialize a datalist object from a list of supported dataset entries"""

#     from . import datalists
#     from . import xyzfile
    
#     try:
#         #utils.echo_msg(data_list)
#         xdls = []
#         for dl in data_list:
#             #utils.echo_msg(dl)
#             if dl is sys.stdin:
#                 xdls.append(xyzfile.XYZFile(fn=dl, data_format=168))
#             elif isinstance(dl, str):
#                 xdls.append(
#                     DatasetFactory(
#                         mod=" ".join(['-' if x == "" else x for x in dl.split(",")]),
#                         data_format = None,
#                         weight=None if not want_weight else 1,
#                         uncertainty=None if not want_uncertainty else 0,
#                         src_srs=src_srs,
#                         dst_srs=dst_srs,
#                         src_geoid=src_geoid,
#                         dst_geoid=dst_geoid,
#                         x_inc=xy_inc[0],
#                         y_inc=xy_inc[1],
#                         sample_alg=sample_alg,
#                         parent=None,
#                         src_region=region,
#                         invert_region=invert_region,
#                         cache_dir=cache_dir,
#                         want_mask=want_mask,
#                         want_sm=want_sm,
#                         verbose=want_verbose,
#                         dump_precision=dump_precision,
#                         pnt_fltrs=pnt_fltrs,
#                         stack_fltrs=stack_fltrs,
#                         stack_node=stack_node,
#                         stack_mode=stack_mode,
#                         upper_limit=None,
#                         lower_limit=None,
#                         mask=mask
#                     )._acquire_module()
#                 )
#             elif isinstance(dl, dict):

#                 if dl['kwargs']['src_region'] is not None:
#                     dl['kwargs']['src_region'] = regions.Region().from_list(
#                         dl['kwargs']['src_region']
#                 )
#                 dl['kwargs']['xinc'] = xy_inc[0]
#                 dl['kwargs']['yinc'] = xy_inc[1]
#                 xdls.append(
#                     DatasetFactory(
#                         mod_name=dl['mod_name'], mod_args=dl['mod_args'], **dl['kwargs']
#                     )._acquire_module()
#                 )

#         #utils.echo_msg(region)
#         utils.echo_msg(xdls)
#         if len(xdls) > 1:
#             this_datalist = datalists.Scratch(
#                 fn=xdls,
#                 data_format=-3,
#                 weight=None if not want_weight else 1,
#                 uncertainty=None if not want_uncertainty else 0,
#                 src_srs=src_srs,
#                 dst_srs=dst_srs,
#                 src_geoid=src_geoid,
#                 dst_geoid=dst_geoid,
#                 x_inc=xy_inc[0],
#                 y_inc=xy_inc[1],
#                 sample_alg=sample_alg,
#                 parent=None,
#                 src_region=region,
#                 invert_region=invert_region,
#                 cache_dir=cache_dir,
#                 want_mask=want_mask,
#                 want_sm=want_sm,
#                 verbose=want_verbose,
#                 dump_precision=dump_precision,
#                 pnt_fltrs=pnt_fltrs,
#                 stack_fltrs=stack_fltrs,
#                 stack_node=stack_node,
#                 stack_mode=stack_mode,
#                 upper_limit=None,
#                 lower_limit=None,
#                 mask=mask
#             )

#         elif len(xdls) > 0:
#             this_datalist = xdls[0]
#         else:
#             this_datalist = None

#         return(this_datalist)
    
#     except Exception as e:
#         #utils.echo_warning_msg('failed to parse datalist obs')
#         utils.echo_error_msg(
#             f'could not initialize data, {data_list}: {e}'
#         )
#         return(None)

def init_data(data_list, **kwargs):
    """initialize a datalist object from a list of supported dataset entries"""

    from . import datalists
    from . import xyzfile
    
    try:
        xdls = []
        kwargs['weight'] = kwargs.get('want_weight')
        kwargs['uncertainty'] = kwargs.get('want_uncertainty')

        for dl in data_list:
            if dl is sys.stdin:
                xdls.append(xyzfile.XYZFile(fn=dl, data_format=168, **kwargs))

            elif isinstance(dl, str):
                kwargs['mod'] = " ".join(['-' if x == "" else x for x in dl.split(",")])
                utils.echo_msg(kwargs)
                xdls.append(DatasetFactory(**kwargs)._acquire_module())

            elif isinstance(dl, dict):
                if dl['kwargs']['src_region'] is not None:
                    dl['kwargs']['src_region'] = regions.Region().from_list(
                        dl['kwargs']['src_region']
                )
                #dl['kwargs']['xinc'] = xy_inc[0]
                #dl['kwargs']['yinc'] = xy_inc[1]
                xdls.append(
                    DatasetFactory(
                        mod_name=dl['mod_name'], mod_args=dl['mod_args'], **dl['kwargs']
                    )._acquire_module()
                )

        if len(xdls) > 1:
            this_datalist = datalists.Scratch(fn=xdls, data_format=-3, **kwargs)
        elif len(xdls) > 0:
            this_datalist = xdls[0]
        else:
            this_datalist = None

        return(this_datalist)
    
    except Exception as e:
        #utils.echo_warning_msg('failed to parse datalist obs')
        utils.echo_error_msg(
            f'could not initialize data, {data_list}: {e}'
        )
        return(None)
    
        
## ==============================================
## ElevationDataset and sub-modules
## ==============================================
class ElevationDataset:
    """representing an Elevation Dataset
    
    This is the super class for all datalist (dlim) datasets .    
    Each dataset sub-class should define a dataset-specific    
    data parser <self.parse> and a <self.generate_inf> function to generate
    inf files. 

    Specifically, each sub-dataset should minimally define the following 
    functions:

    sub_ds.__init__
    sub_ds.yield_points

    -----------
    Where:

    generate_inf generates a dlim compatible inf file,
    parse yields the dlim dataset module
    yield_points yields the data as a numpy rec-array with 'x', 'y', 'z' 
    and optionally 'w' and 'u'

    -----------
    Parameters:

    fn: dataset filename or fetches module
    data_format: dataset format
    weight: dataset weight
    uncertainty: dataset uncertainty
    mask: mask the dataset
    src_srs: dataset source srs
    dst_srs: dataset target srs
    src_geoid: dataset source geoid (if applicable)
    dst_geoid: dataset target geoid
    x_inc: target dataset x/lat increment
    y_inc: target dataset y/lon increment
    want_mask: generate a mask of the data
    want_sm: generate spatial metadata vector
    sample_alg: the gdal resample algorithm
    parent: dataset parent obj
    region: ROI
    invert_region: invert the region
    stack_node: output x/y data as mean of input x/y instead of cell center
    stack_mode: 'min', 'max', 'supercede', 'mean'
    cache_dir: cache_directory
    verbose: be verbose
    remote: dataset is remote
    dump_precision: specify the float precision of the dumped xyz data
    params: the factory parameters
    metadata: dataset metadata
    """

    gdal_sample_methods = [
        'near', 'bilinear', 'cubic', 'cubicspline', 'lanczos',
        'average', 'mode',  'max', 'min', 'med', 'Q1', 'Q3', 'sum',
        'auto'
    ]

    stack_modes = [
        'min', 'max', 'mean', 'supercede', 'mixed'
    ]

    ## todo: add transformation grid option
    ## (stacks += transformation_grid), geoids
    def __init__(self,
                 fn=None,
                 data_format=None,
                 weight=1,
                 uncertainty=0,
                 src_srs=None,
                 mask=None,
                 dst_srs='epsg:4326',
                 src_geoid=None,
                 dst_geoid='g2018',
                 x_inc=None,
                 y_inc=None,
                 want_mask=False,
                 want_sm=False,
                 sample_alg='auto',
                 parent=None,
                 src_region=None,
                 invert_region=False,
                 pnt_fltrs=[],
                 stack_fltrs=[],
                 stack_node=True,
                 stack_mode='mean',
                 cache_dir=None,
                 verbose=False,
                 remote=False,
                 dump_precision=6,
                 upper_limit=None,
                 lower_limit=None,
                 params={},
                 metadata={
                     'name':None, 'title':None, 'source':None, 'date':None,
                     'data_type':None, 'resolution':None, 'hdatum':None,
                     'vdatum':None, 'url':None
                 },
                 **kwargs):
        self.fn = fn # dataset filename or fetches module
        self.data_format = data_format # dataset format
        self.weight = weight # dataset weight
        self.uncertainty = uncertainty # dataset uncertainty
        self.mask = mask #{'mask': mask, 'invert': False}#mask # dataset mask
        self.src_srs = src_srs # dataset source srs
        self.dst_srs = dst_srs # dataset target srs
        self.src_geoid = src_geoid # dataset source geoid
        self.dst_geoid = dst_geoid # dataset target geoid
        self.x_inc = utils.str2inc(x_inc) # target dataset x/lat increment
        self.y_inc = utils.str2inc(y_inc) # target dataset y/lon increment
        self.stack_fltrs = stack_fltrs # pass the stack, if generated, through grits filters
        self.pnt_fltrs = pnt_fltrs # pass the points through filters
        self.want_mask = want_mask # mask the data
        self.want_sm = want_sm # generate spatial metadata vector
        if self.want_sm:
            self.want_mask = True
            
        self.sample_alg = sample_alg # the gdal resample algorithm
        self.metadata = copy.deepcopy(metadata) # dataset metadata
        self.parent = parent # dataset parent obj
        self.region = src_region # ROI
        if self.region is not None:
            self.region = regions.Region().from_user_input(self.region)
            
        self.invert_region = invert_region # invert the region
        self.cache_dir = cache_dir # cache_directory
        self.verbose = verbose # be verbose
        self.remote = remote # dataset is remote
        self.dump_precision = dump_precision # the precision of the dumped xyz data
        self.stack_node = stack_node # yield avg x/y data instead of center
        self.stack_mode = stack_mode # 'mean', 'min', 'max', 'supercede'
        self._init_stack_mode()

        self.upper_limit = utils.float_or(upper_limit)
        self.lower_limit = utils.float_or(lower_limit)

        self.infos = inf.INF(
            name=self.fn, file_hash='0', numpts=0, fmt=self.data_format
        ) # infos blob

        # for kwarg in kwargs:
        #     if kwarg not in self.__dict__.keys():
        #         utils.echo_warning_msg('{} is not a valid parameter'.format(kwarg))
        
        self.params = params # the factory parameters
        if not self.params:
            self.params['kwargs'] = self.__dict__.copy()
            self.params['mod'] = self.fn
            self.params['mod_name'] = self.data_format
            self.params['mod_args'] = {}

            
    def __str__(self):
        return(f'<Dataset: {self.metadata["name"]} - {self.fn}>')

    
    def __repr__(self):
        return(f'<Dataset: {self.metadata["name"]} - {self.fn}>')

    
    def __call__(self):
        self.initialize()

        
    def _init_stack_mode(self):
        opts, self.stack_mode_name, self.stack_mode_args \
            = factory.parse_fmod(self.stack_mode)

        if self.stack_mode_name not in self.stack_modes:
            utils.echo_warning_msg(
                f'{self.stack_mode_name} is not a valid stack mode'
            )
            self.stack_mode_name = 'mean'
        
        if 'mask_level' not in self.stack_mode_args.keys():
            self.stack_mode_args['mask_level'] = -1

            
    def _init_mask(self, mask):
        opts = factory.fmod2dict(mask, {})
        if 'mask_fn' not in opts.keys():
            if '_module' in opts.keys():
                opts['mask_fn'] = opts['_module']
            else:
                utils.echo_error_msg(f'could not parse mask {self.mask}')
                return(None)

        if 'invert' not in opts.keys():
            opts['invert'] = False
            
        if 'verbose' not in opts.keys():
            opts['verbose'] = False

        if 'min_z' not in opts.keys():
            opts['min_z'] = None
        else:
            opts['min_z'] = utils.float_or(opts['min_z'])

        if 'max_z' not in opts.keys():
            opts['max_z'] = None
        else:
            opts['max_z'] = utils.float_or(opts['max_z'])
            
        # mask is ogr, rasterize it
        opts['ogr_or_gdal'] = gdalfun.ogr_or_gdal(opts['mask_fn'])
        if opts['ogr_or_gdal'] == 1: 
            if self.region is not None \
               and self.x_inc is not None \
               and self.y_inc is not None:
                data_mask = gdalfun.ogr2gdal_mask(
                    opts['mask_fn'],
                    region=self.region,
                    x_inc=self.x_inc,
                    y_inc=self.y_inc,
                    dst_srs=self.dst_srs,
                    invert=True,
                    verbose=True,
                    temp_dir=self.cache_dir,
                )
                opts['ogr_or_gdal'] = 0
            else:    
                data_mask = opts['mask_fn']
        else:
            data_mask = opts['mask_fn']

        opts['data_mask'] = data_mask
        return(opts)
        
                    
    def _set_params(self, **kwargs):
        metadata = copy.deepcopy(self.metadata)
        #metadata['name'] = self.fn
        _params = {
            'parent': self,
            'weight': self.weight,
            'uncertainty': self.uncertainty,
            'src_region': self.region,
            'invert_region': self.invert_region,
            'x_inc': self.x_inc,
            'y_inc': self.y_inc,
            'mask': self.mask,            
            'src_srs': self.src_srs,
            'dst_srs': self.dst_srs,
            'stack_fltrs': self.stack_fltrs,
            'pnt_fltrs': self.pnt_fltrs,
            'cache_dir': self.cache_dir,
            'upper_limit': self.upper_limit,
            'lower_limit': self.lower_limit,
            'verbose': self.verbose,
            'metadata': metadata,
            'stack_mode': self.stack_mode,
        }
        for kw in kwargs.keys():
            _params[kw] = kwargs[kw]

        return(_params)

    
    def _copy_params(self, **kwargs):
        _params = {i:self.params['kwargs'][i] for i in self.params['kwargs'] if i!='params'}
        for kw in kwargs.keys():
            _params[kw] = kwargs[kw]

        return(_params)


    def _sub_init(self):
        pass
        
    
    def initialize(self):
        self._sub_init()
        self._fn = None # temp filename holder
        self.data_region = None # self.region and inf.region reduced
        self.inf_region = None # inf region
        self.archive_datalist = None # the datalist of the archived data
        self.data_entries = [] # 
        self.data_lists = {} #
        self.cache_dir = utils.cudem_cache() \
            if self.cache_dir is None \
               else self.cache_dir # cache directory
        if self.sample_alg not in self.gdal_sample_methods: # gdal_warp resmaple algorithm 
            utils.echo_warning_msg(
                (f'{self.sample_alg} is not a valid gdal warp resample algorithm, '
                 'falling back to `auto`')
            )
            self.sample_alg = 'auto'

        if utils.fn_url_p(self.fn):
            self.remote = True

        ## initialize transform
        self.transform = {
            'src_horz_crs': None,
            'dst_horz_crs': None,
            'src_vert_crs': None,
            'dst_vert_crs': None,
            'src_vert_epsg': None,
            'dst_vert_epsg': None,
            'pipeline': None,
            'trans_fn': None,
            'trans_fn_unc': None,
            'trans_region': None,
            'transformer': None,
            'vert_transformer': None,
            'want_vertical': False,
        } # pyproj transformation info


        if self.valid_p():
            #try:
            self.infos = self.inf(
                check_hash=True if self.data_format == -1 else False
            )
            # except:
            #    utils.echo_error_msg(
            #        f'could not parse dataset {self.fn}'
            #    )
            #    return(self)
             
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                try:
                    srsfun.set_transform(
                        src_srs=self.src_srs,
                        dst_srs=self.dst_srs,
                        region=self.region,
                        infos=self.infos,
                        cache_dir=self.cache_dir
                    )
                    #self.set_transform()
                except Exception as e:
                    utils.echo_error_msg(
                        f'could not set transformation on {self.fn}, {e}'
                    )

            self.set_yield(use_blocks=False)
            #self.set_yield()

        if self.pnt_fltrs is not None and isinstance(self.pnt_fltrs, str):
            self.pnt_fltrs = [self.pnt_fltrs]

        if self.mask is not None and isinstance(self.mask, str):
            self.mask = [self.mask]
            
        return(self)

    
    def fetch(self):
        """fetch remote data from self.data_entries"""
        
        for entry in self.data_entries:
            if entry.remote:
                if entry._fn is None:
                    entry._fn = os.path.basename(self.fn)
                    
                f = utils.Fetch(
                    url=entry.fn, verbose=entry.verbose
                )
                if f.fetch_file(entry._fn) == 0:
                    entry.fn = entry._fn
            else:
                utils.echo_warning_msg('nothing to fetch')

                
    def valid_p(self, fmts=['scratch']):
        """check if self appears to be a valid dataset entry"""

        if self.fn is None:
            return(False)

        if self.data_format is None:
            return(False)

        #check_path=True
        # if self.fn.startswith('http') \
        #    or self.fn.startswith('/vsicurl/') \
        #    or self.fn.startswith('BAG'):
        #     check_path = False

        if self.fn is sys.stdin or self.fn == '-':
            self.fn = sys.stdin
            return(True)
        
        #if check_path:
        if self.fn is not None:
            if not isinstance(self.fn, list) and \
               not isinstance(self.fn, np.ndarray) and \
               not isinstance(self.fn, np.core.records.recarray):
                if self.fn not in fmts:
                    if not self.fn.startswith('http') \
                       and not self.fn.startswith('/vsicurl/') \
                       and not self.fn.startswith('BAG') \
                       and not ':' in self.fn:
                        if not utils.fn_url_p(self.fn):
                            if self.data_format > -10:
                                if not os.path.exists(self.fn):
                                    utils.echo_warning_msg(f'{self.fn} does not exist')
                                    return (False)

                                if os.stat(self.fn).st_size == 0:
                                    utils.echo_warning_msg(f'{self.fn} is 0 bytes')
                                    return(False)
                        
        return(True)

    
    def echo_(self, sep=' ', **kwargs):
        """print self as a datalist entry string"""

        out = []
        for key in self.metadata.keys():
            if key != 'name':
                out.append(str(self.metadata[key]))

        return(sep.join(['"{}"'.format(str(x)) for x in out]))

    
    def echo(self, sep=' ', **kwargs):
        """print self.data_entries as a datalist entries."""

        out = []
        for entry in self.parse():
            entry_path = os.path.abspath(entry.fn) if not self.remote else entry.fn
            l = [entry_path, entry.data_format]
            if entry.weight is not None:
                l.append(entry.weight)

            entry_string = '{}'.format(sep.join([str(x) for x in l]))
            out.append(entry_string)
            print(entry_string)

        return(out)


    def format_data(self, sep=' '):
        out = []        
        for entry in self.parse():
            out.append(entry.format_entry(sep=sep))

        return(out)

    
    def format_entry(self, sep=' '):
        """format the dataset information as a `sep` separated string."""

        dl_entry = sep.join(
            [str(x) for x in [
                self.fn,
                f'{self.data_format}:{factory.dict2args(self.params["mod_args"])}',
                self.weight,
                self.uncertainty
            ]]
        )
        metadata = self.echo_()
        return(sep.join([dl_entry, metadata]))

            
    def format_metadata(self, **kwargs):
        """format metadata from self, for use as a datalist entry."""

        return(self.echo_())

    
    def set_yield(self, use_blocks=False):
        """set the yield strategy, either default (all points) or mask or stacks

        This sets both `self.array_yeild` and `self.xyz_yield`.
        
        if region, x_inc and y_inc are set, data will yield through stacks,
        otherwise, data will yield directly from `self.yield_xyz` and `self.yield_array`
        """

        #self.array_yield = self.yield_array()
        #self.xyz_yield = self.yield_xyz()    
        self.array_yield = self.mask_and_yield_array()
        self.xyz_yield = self.mask_and_yield_xyz()
        #if self.want_archive: # archive only works when yielding xyz data.
        #    self.xyz_yield = self.archive_xyz()
        # if (self.region is None \
        #     and self.transform['trans_region'] is None) \
        #     and self.x_inc is not None:
        if self.region is None and self.x_inc is not None:
            utils.echo_warning_msg(
                'must enter a region to output in increment blocks...'
            )

        if self.region is not None and self.x_inc is not None:
            self.x_inc = utils.str2inc(self.x_inc)
            if self.y_inc is None:
                self.y_inc = self.x_inc
            else:
                self.y_inc = utils.str2inc(self.y_inc)

            #out_name = utils.make_temp_fn('dlimstacks', temp_dir=self.cache_dir)
            # out_name = utils.make_temp_fn(
            #     utils.append_fn(
            #         'dlimstacks', self.region, self.x_inc
            #     ), temp_dir=self.cache_dir
            # )
            
            out_name = os.path.join(self.cache_dir, utils.append_fn('globato', self.region, self.x_inc))
            if not use_blocks:
                self.xyz_yield = self.stacks_yield_xyz(out_name=out_name)
            else:
                self.xyz_yield = self.blocks_yield_xyz(out_name=out_name)

            
    ## initialize transfom
    def parse_srs(self):
        if self.src_srs is not None and self.dst_srs is not None:
            want_vertical = True
            src_geoid = None
            dst_geoid = 'g2018' # dst_geoid is g2018 by default.
            in_vertical_epsg = in_vertical_crs = None
            out_vertical_epsg = out_vertical_crs = None
            
            ## parse out the source geoid, which is not standard in proj,
            ## reset `self.src_srs` without it.
            tmp_src_srs = self.src_srs.split('+geoid:')
            src_srs = tmp_src_srs[0]
            self.src_srs = src_srs
            if len(tmp_src_srs) > 1:
                src_geoid = tmp_src_srs[1]

            ## parse out the destination geoid, which is not standard in proj
            ## reset `self.dst_srs` without it.
            tmp_dst_srs = self.dst_srs.split('+geoid:')
            dst_srs = tmp_dst_srs[0]
            self.dst_srs = dst_srs
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
            try:
                in_crs = pyproj.CRS.from_user_input(src_srs)
                out_crs = pyproj.CRS.from_user_input(dst_srs)
            except:
                utils.echo_error_msg([src_srs, dst_srs])

            ## if the crs has vertical (compound), parse out the vertical crs
            ## and set the horizontal and vertical crs
            if in_crs.is_compound:
                in_crs_list = in_crs.sub_crs_list
                in_horizontal_crs = in_crs_list[0]
                in_vertical_crs = in_crs_list[1]
                in_vertical_name = in_vertical_crs.name
                in_vertical_epsg = in_vertical_crs.to_epsg()
                if in_vertical_epsg is None:
                    in_vertical_epsg = in_vertical_name.split(' ')[0]
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
                if (in_vertical_epsg == out_vertical_epsg) and self.src_geoid is None:
                    want_vertical = False

            ## set `self.transform` with the parsed srs info
            self.transform['src_horz_crs'] = in_horizontal_crs
            self.transform['dst_horz_crs'] = out_horizontal_crs
            self.transform['src_vert_crs'] = in_vertical_crs
            self.transform['dst_vert_crs'] = out_vertical_crs
            self.transform['src_vert_epsg'] = in_vertical_epsg
            self.transform['dst_vert_epsg'] = out_vertical_epsg
            self.transform['src_geoid'] = src_geoid
            self.transform['dst_geoid'] = dst_geoid
            self.transform['want_vertical'] = want_vertical        

            
    def set_vertical_transform(self):
        ## set the region for the vdatum transformation grid.
        ## this is either from the input `self.region` or from the input
        ## data's bounding box. Transform it to WGS84.
        if self.region is None:
            vd_region = regions.Region().from_list(self.infos.minmax)
            vd_region.src_srs = self.transform['src_horz_crs'].to_proj4()
        else:
            vd_region = self.region.copy()
            vd_region.src_srs = self.transform['dst_horz_crs'].to_proj4()

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
        self.transform['trans_fn'] = os.path.join(
            self.cache_dir, '_vdatum_trans_{}_{}_{}.tif'.format(
                self.transform['src_vert_epsg'],
                self.transform['dst_vert_epsg'],
                vd_region.format('fn')
            )
        )
        
        ## if the transformation grid already exists, skip making a new one,
        ## otherwise, make the new one here with `vdatums.VerticalTransform()`
        if not os.path.exists(self.transform['trans_fn']):
            with utils.ccp(
                    desc='generating vertical transformation grid {} from {} to {}'.format(
                        self.transform['trans_fn'], self.transform['src_vert_epsg'],
                        self.transform['dst_vert_epsg']
                    ),
                    leave=self.verbose
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
                self.transform['trans_fn'], self.transform['trans_fn_unc'] \
                    = vdatums.VerticalTransform(
                        'IDW',
                        vd_region,
                        vd_x_inc,
                        vd_y_inc,
                        self.transform['src_vert_epsg'],
                        self.transform['dst_vert_epsg'],
                        geoid_in=self.transform['src_geoid'],
                        geoid_out=self.transform['dst_geoid'],
                        cache_dir=self.cache_dir,
                        verbose=False
                    ).run(outfile=self.transform['trans_fn'])

        ## set the pyproj.Transformer for both horz+vert and vert only
        ## hack for navd88 datums in feet (6360 is us-feet, 8228 is international-feet
        if utils.str_or(self.transform['src_vert_epsg']) == '6360':
            # or 'us-ft' in utils.str_or(src_vert, ''):
            #out_src_srs = out_src_srs + ' +vto_meter=0.3048006096012192'
            uc = ' +step +proj=unitconvert +z_in=us-ft +z_out=m'
        elif utils.str_or(self.transform['src_vert_epsg']) == '8228':
            uc = ' +step +proj=unitconvert +z_in=ft +z_out=m'
        else:
            uc = ''
            
        if self.transform['trans_fn'] is not None \
           and os.path.exists(self.transform['trans_fn']):
            self.transform['pipeline'] \
                = (f'+proj=pipeline{uc} +step '
                   f'{self.transform["src_horz_crs"].to_proj4()} '
                   '+inv +step +proj=vgridshift '
                   f'+grids="{os.path.abspath(self.transform["trans_fn"])}" '
                   f'+inv +step {self.transform["dst_horz_crs"].to_proj4()}')
            self.transform['vert_transformer'] = pyproj.Transformer.from_pipeline(
                (f'+proj=pipeline{uc} +step +proj=vgridshift '
                 f'+grids="{os.path.abspath(self.transform["trans_fn"])}" +inv')
            )
            utils.echo_debug_msg(self.transform['pipeline'])
        else:
            utils.echo_error_msg(
                ('failed to generate vertical transformation grid between '
                 f'{self.transform["src_vert_epsg"]} and '
                 f'{os.path.abspath(self.transform["dst_vert_epsg"])} '
                 'for this region!')
            )

            
    def set_transform(self):
        """Set the pyproj horizontal and vertical transformations 
        for the dataset
        """

        self.parse_srs()
        if self.transform['src_horz_crs'] is not None \
           and self.transform['dst_horz_crs'] is not None:        
            ## horizontal Transformation
            self.transform['horz_pipeline'] \
                = ('+proj=pipeline +step '
                   f'{self.transform["src_horz_crs"].to_proj4()} '
                   f'+inv +step {self.transform["dst_horz_crs"].to_proj4()}')
            if self.region is not None:
                self.transform['trans_region'] = self.region.copy()
                self.transform['trans_region'].src_srs \
                    = self.transform['dst_horz_crs'].to_proj4()
                self.transform['trans_region'].warp(
                    self.transform['src_horz_crs'].to_proj4()
                )
            else:
                self.transform['trans_region'] \
                    = regions.Region().from_list(self.infos.minmax)
                self.transform['trans_region'].src_srs \
                    = self.infos.src_srs
                self.transform['trans_region'].warp(
                    self.transform['dst_horz_crs'].to_proj4()
                )

            ## vertical Transformation
            if self.transform['want_vertical']:
                self.set_vertical_transform()
            else:
                self.transform['pipeline'] = self.transform['horz_pipeline']

            try:
                self.transform['transformer'] \
                    = pyproj.Transformer.from_pipeline(
                        self.transform['pipeline']
                    )
            except Exception as e:
                utils.echo_warning_msg(
                    ('could not set transformation in: '
                     f'{self.transform["src_horz_crs"].name}, out: '
                     f'{self.transform["dst_horz_crs"].name}, {e}')
                    )

                return

        ## dataset region
        ## mrl: moved out of if block 
        if self.region is not None and self.region.valid_p():
            self.data_region = self.region.copy() \
                if self.transform['trans_region'] is None \
                   else self.transform['trans_region'].copy()
            #inf_region = regions.Region().from_list(self.infos.minmax)
            self.data_region = regions.regions_reduce(
                self.data_region, self.inf_region
            )
            self.data_region.src_srs = self.infos.src_srs

            if not self.data_region.valid_p():
                self.data_region = self.region.copy() \
                    if self.transform['trans_region'] is None \
                       else self.transform['trans_region'].copy()
        else:
            self.data_region = regions.Region().from_list(
                self.infos.minmax
            )
            self.data_region.src_srs = self.infos.src_srs
        #     self.region = self.data_region.copy()
        # self.region.zmax=self.upper_limit
        # self.region.zmin=self.lower_limit

        
    ## INF file reading/writing
    def generate_inf(self):
        """generate an inf file for the data source. this is generic and
        will parse through all the data via yield_xyz to get point count
        and region. if the datasource has a better way to do this, such as
        with las/gdal files, re-define this function in its respective
        dataset sub-class.

        todo: hull for wkt
        ## todo: put qhull in a yield and perform while scanning for min/max
        # try:
        #     out_hull = [pts[i] for i in ConvexHull(
        #         pts, qhull_options='Qt'
        #     ).vertices]
        #     out_hull.append(out_hull[0])
        #     self.infos['wkt'] = regions.create_wkt_polygon(out_hull, xpos=0, ypos=1)
        # except:
        #     ...
        """

        _region = None
        if self.region is not None:
            _region = self.region.copy()
            self.region = None

        _x_inc = None
        if self.x_inc is not None:
            _x_inc = self.x_inc
            self.x_inc = None
            
        this_region = regions.Region()
        point_count = 0

        for points in self.transform_and_yield_points():
            if point_count == 0:
                this_region.from_list(
                    [
                        points['x'].min(), points['x'].max(),
                        points['y'].min(), points['y'].max(),
                        points['z'].min(), points['z'].max()
                    ]
                )
            else:
                if points['x'].min() < this_region.xmin:
                    this_region.xmin = points['x'].min()
                elif points['x'].max() > this_region.xmax:
                    this_region.xmax = points['x'].max()
                    
                if points['y'].min() < this_region.ymin:
                    this_region.ymin = points['y'].min()
                elif points['y'].max() > this_region.ymax:
                    this_region.ymax = points['y'].max()
                    
                if points['z'].min() < this_region.zmin:
                    this_region.zmin = points['z'].min()
                elif points['z'].min() > this_region.zmax:
                    this_region.zmax = points['z'].min()
                
            point_count += len(points)

        self.infos.numpts = point_count
        if point_count > 0:
            self.infos.minmax = this_region.export_as_list(
                include_z=True
            )
            self.infos.wkt = this_region.export_as_wkt()
            
        self.infos.src_srs = self.src_srs
        if _region is not None:
            self.region = _region.copy()

        if _x_inc is not None:
            self.x_inc = _x_inc
            
        return(self.infos)

    
    def inf(self, check_hash=False, recursive_check=False,
            write_inf=True, **kwargs):
        """read/write an inf file

        If the inf file is not found, will attempt to generate one.
        The function `generate_inf` should be defined for each specific
        dataset sub-class.
        """

        inf_path = '{}.inf'.format(self.fn)
        generate_inf = False

        ## try to parse the existing inf file as either a native
        ## inf json file
        if os.path.exists(inf_path):
            try:
                self.infos.load_inf_file(inf_path)
                
            except ValueError as e:
                generate_inf = True
                utils.remove_glob(inf_path)
                if self.verbose:
                    utils.echo_warning_msg(
                        'failed to parse inf {}, {}'.format(inf_path, e)
                    )
        else:
            generate_inf = True        

        if self.fn is sys.stdin:
            generate_inf = False
        ## check hash from inf file vs generated hash,
        ## if hashes are different, then generate a new
        ## inf file...only do this if check_hash is set
        ## to True, as this can be time consuming and not
        ## always necessary...
        if check_hash:
            generate_inf = self.infos.generate_hash() != self.infos.file_hash

        ## this being set can break some modules (bags esp)
        # if self.remote:
        #     generate_inf = False

        if self.src_srs is not None:
            self.infos.src_srs = self.src_srs
        else:                            
            self.src_srs = self.infos.src_srs
        
        if generate_inf:
            self.infos = self.generate_inf()
            
            ## update this
            if self.data_format >= -2 and write_inf:
                self.infos.write_inf_file()

            if recursive_check and self.parent is not None:
                self.parent.inf(check_hash=True)

        try:
            self.inf_region = regions.Region().from_string(self.infos.wkt)
        except:
            try:
                self.inf_region = regions.Region().from_list(
                    self.infos.minmax
                )
            except:
                if self.region is not None:
                    self.inf_region = self.region.copy()
                else:
                    self.inf_region = None

        #self.transform['trans_region'] = self.inf_region.copy()
        if self.transform['transformer'] is not None \
           and self.transform['trans_region'] is None:
            ## trans_region is the inf_region reprojected to dst_srs
            ## self.region should be supplied in the same srs as dst_srs
            self.transform['trans_region'] = self.inf_region.copy()
            self.transform['trans_region'].src_srs = self.infos.src_srs
            self.transform['trans_region'].warp(self.dst_srs)

        # if self.region is None:
        #     if self.transform['trans_region'] is not None:
        #         self.region = self.transform['trans_region'].copy()
        #     else:
        #         self.region = self.inf_region.copy()

        return(self.infos)

    
    ## Datalist/Dataset parsing, reset in sub-dataset if needed
    def parse(self):
        """parse the datasets from the dataset.
        
        Re-define this method when defining a dataset sub-class
        that represents recursive data structures (datalists, zip, etc).

        This will fill self.data_entries

        -------
        Yields:

        dataset object
        """

        if self.region is not None and self.inf_region is not None:
            # try:
            #     inf_region = regions.Region().from_string(self.infos.wkt)
            # except:
            #     try:
            #         inf_region = regions.Region().from_list(self.infos.minmax)
            #     except:
            #         inf_region = self.region.copy()

            if regions.regions_intersect_p(
                    self.inf_region,
                    self.region \
                    if self.transform['trans_region'] is None \
                    else self.transform['trans_region']
            ):
                self.data_entries.append(self)
                yield(self)
                
        else:
            self.data_entries.append(self)
            yield(self)

            
    def parse_data_lists(self, gather_data=True):
        """parse the data into a datalist dictionary"""
        
        for e in self.parse():
            if e.parent is not None:
                if e.parent.metadata['name'] in self.data_lists.keys():
                    self.data_lists[e.parent.metadata['name']]['data'].append(e)
                else:
                    self.data_lists[e.parent.metadata['name']] = {
                        'data': [e], 'parent': e.parent
                    }
                    
            else:
                self.data_lists[e.metadata['name']] \
                    = {'data': [e], 'parent': e}
                
        return(self)
    
    
    ## Data Yield
    def yield_xyz_from_entries(self):
        """yield from self.data_entries, list of datasets"""
        
        for this_entry in self.data_entries:
            for xyz in this_entry.xyz_yield:
                yield(xyz)
                
            if this_entry.remote:
                utils.remove_glob('{}*'.format(this_entry.fn))

                
    def yield_entries(self):
        """yield from self.data_entries, list of datasets"""
        
        for this_entry in self.data_entries:
            yield(this_entry)

    ## ==============================================            
    ## yield points, where points are a numpy rec-array with 'xyzwu'
    ## ==============================================            
    def yield_points(self):
        """yield the numpy xyz rec array points from the dataset.

        reset in dataset if needed.
        """

        for ds in self.parse():
            for points in ds.yield_points():
                yield(points)
                

    def transform_and_yield_points(self):
        """points are an array of `points['x']`, `points['y']`, `points['z']`, 
        <`points['w']`, `points['u']`>`

        points will be transformed here, based on `self.transformer`, 
        which is set in `set_transform` after the points are transformed, the data 
        will pass through the region, if it exists and finally
        yield the transformed and reduced points.
        """            

        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            for points in self.yield_points():
                if self.transform['transformer'] is not None \
                   or self.transform['vert_transformer'] is not None:
                    if self.transform['transformer'] is not None:
                        points['x'], points['y'] \
                            = self.transform['transformer'].transform(
                                points['x'], points['y']
                            )
                    if self.transform['vert_transformer'] is not None:
                        _, _, points['z'] \
                            = self.transform['vert_transformer'].transform(
                                points['x'], points['y'], points['z']
                            )
                        
                    points = points[~np.isinf(points['z'])]

                if self.region is not None and self.region.valid_p():
                    xyz_region = self.region.copy()
                    #if self.transform['trans_region'] is None \
                        #else self.transform['trans_region'].copy()
                    if self.invert_region:
                        points = points[
                            ((points['x'] >= xyz_region.xmax) \
                             | (points['x'] <= xyz_region.xmin)) \
                            | ((points['y'] >= xyz_region.ymax) \
                               | (points['y'] <= xyz_region.ymin))
                        ]
                        if xyz_region.zmin is not None:
                            points = points[(points['z'] <= xyz_region.zmin)]

                        if xyz_region.zmax is not None:
                            points = points[(points['z'] >= xyz_region.zmax)]

                        if xyz_region.wmin is not None:
                            points = points[(points['w'] <= xyz_region.wmin)]
                            
                        if xyz_region.wmax is not None:
                            points = points[(points['w'] >= xyz_region.wmax)]
                            
                        if xyz_region.umin is not None:
                            points = points[(points['u'] <= xyz_region.umin)]
                            
                        if xyz_region.umax is not None:
                            points = points[(points['u'] >= xyz_region.umax)]
                            
                    else:
                        points = points[
                            ((points['x'] <= xyz_region.xmax) \
                             & (points['x'] >= xyz_region.xmin)) & \
                            ((points['y'] <= xyz_region.ymax) \
                             & (points['y'] >= xyz_region.ymin))
                        ]
                        if xyz_region.zmin is not None:
                            points = points[(points['z'] >= xyz_region.zmin)]

                        if xyz_region.zmax is not None:
                            points = points[(points['z'] <= xyz_region.zmax)]

                        if xyz_region.wmin is not None:
                            points = points[(points['w'] >= xyz_region.wmin)]
                            
                        if xyz_region.wmax is not None:
                            points = points[(points['w'] <= xyz_region.wmax)]
                            
                        if xyz_region.umin is not None:
                            points = points[(points['u'] >= xyz_region.umin)]
                            
                        if xyz_region.umax is not None:
                            points = points[(points['u'] <= xyz_region.umax)]

                if self.upper_limit is not None:
                    points = points[(points['z'] <= self.upper_limit)]
                    
                if self.lower_limit is not None:
                    points = points[(points['z'] >= self.lower_limit)]

                if len(points) > 0:
                    ## apply any dlim filters to the points
                    #utils.echo_debug_msg(points)
                    if isinstance(self.pnt_fltrs, list):
                        if self.pnt_fltrs is not None:
                            for f in self.pnt_fltrs:
                                point_filter = pointz.PointFilterFactory(
                                    mod=f,
                                    points=points,
                                    region=self.region,
                                    xyinc=[self.x_inc, self.y_inc],
                                    cache_dir=self.cache_dir,
                                )._acquire_module()
                                if point_filter is not None:
                                    points = point_filter()

                    if len(points) > 0:
                        yield(points)

        self.transform['transformer'] = self.transform['vert_transformer'] = None


    ## ==============================================            
    ## mask and yield, mask the incoming data to the masks specified in
    ## self.mask, which is a list of factory formatted strings, such as:
    ## ['masks/cudem_masks.shp:verbose=True:invert=False',
    ##  'mask_fn=usgs/OzetteLake.gpkg:verbose=True:invert=False']
    ## ==============================================            
    def mask_and_yield_array(self):
        """mask the incoming array from `self.yield_array` 
        and yield the results.
        
        the mask should be either an ogr supported vector or a 
        gdal supported raster.
        """
        
        mask_band = None
        mask_infos = None
        data_mask = None
        mask_count = 0
        data_masks = []
        if self.mask is not None:
            for mask in self.mask:
                opts = self._init_mask(mask)
                data_masks.append(opts)

        for out_arrays, this_srcwin, this_gt in self.yield_array():
            for data_mask in data_masks:
                if data_mask is None:
                    continue
                
                if os.path.exists(data_mask['data_mask']):
                    utils.echo_debug_msg(f'using mask dataset: {data_mask} to array')                        
                    src_mask = gdal.Open(data_mask['data_mask'])
                    if src_mask is not None:
                        mask_band = src_mask.GetRasterBand(1)
                        mask_infos = gdalfun.gdal_infos(src_mask)
                    else:
                        mask_band = None
                else:
                    utils.echo_warning_msg(f'could not load mask {data_mask["data_mask"]}')

                if data_mask['min_z'] is not None or data_mask['max_z'] is not None:
                    if data_mask['min_z'] is not None and data_mask['max_z'] is not None:
                        z_mask = ((out_arrays['z'] > data_mask['min_z']) & (out_arrays['z'] < data_mask['max_z']))
                    elif data_mask['min_z'] is not None:
                        z_mask = out_arrays['z'] > data_mask['min_z']
                    else:
                        z_mask = out_arrays['z'] < data_mask['max_z']
                else:
                    z_mask = np.full(out_arrays['z'].shape, True)

                if mask_band is not None:
                    ycount, xcount = out_arrays['z'].shape
                    this_region = regions.Region().from_geo_transform(
                        this_gt, xcount, ycount
                    )
                    mask_data = mask_band.ReadAsArray(*this_srcwin)
                    if mask_data is None or len(mask_data) == 0:
                        continue
                        
                    if not np.isnan(mask_infos['ndv']):
                        mask_data[mask_data==mask_infos['ndv']] = np.nan
                        
                    out_mask = ((~np.isnan(mask_data)) & (z_mask))                    
                    for arr in out_arrays.keys():
                        if out_arrays[arr] is not None:
                            if data_mask['invert']:
                                if arr == 'count':
                                    out_arrays[arr][~out_mask] = 0
                                else:                            
                                    out_arrays[arr][~out_mask] = np.nan

                            else:
                                if arr == 'count':
                                    out_arrays[arr][out_mask] = 0
                                else:
                                    out_arrays[arr][out_mask] = np.nan

                    mask_count += np.count_nonzero(~out_mask) \
                        if data_mask['invert'] \
                           else np.count_nonzero(out_mask)

                    mask_data = None
                    
                src_mask = mask_band = None
                
            yield(out_arrays, this_srcwin, this_gt)
                
        #if data_mask is not None:
        #utils.remove_glob('{}*'.format(data_mask))
        if self.mask is not None and self.verbose and mask_count > 0:
            utils.echo_msg_bold(
                'masked {} data records from {}'.format(
                    mask_count, self.fn,
                )
            )


    def mask_and_yield_xyz(self):
        """mask the incoming xyz data from `self.yield_xyz` and yield 
        the results.
        
        the mask should be either an ogr supported vector or a gdal 
        supported raster.
        """
        
        for this_entry in self.parse():
            if this_entry.mask is None:
                for this_xyz in this_entry.yield_xyz():
                    yield(this_xyz)
            else:
                data_masks = []
                mask_count = 0
                for mask in this_entry.mask:
                    opts = self._init_mask(mask)
                    data_masks.append(opts)

                utils.echo_debug_msg(
                    f'using mask dataset: {data_masks} to xyz'
                )                    
                for this_xyz in this_entry.yield_xyz():
                    masked = False
                    for data_mask in data_masks:                
                        if data_mask is None:
                            continue

                        z_masked = True
                        if data_mask['min_z'] is not None or data_mask['max_z'] is not None:
                            if data_mask['min_z'] is not None and data_mask['max_z'] is not None:
                                z_masked = (this_xyz.z > data_mask['min_z'] \
                                            & this_xyz.z < data_mask['mask_z'])
                            elif data_mask['min_z'] is not None:
                                z_masked = this_xyz.z > data_mask['min_z']
                            else:
                                z_masked = this_xyz.z < data_mask['max_z']

                        if data_mask['ogr_or_gdal'] == 0:
                            src_ds = gdal.Open(data_mask['data_mask'])
                            if src_ds is not None:
                                ds_config = gdalfun.gdal_infos(src_ds)
                                ds_band = src_ds.GetRasterBand(1)
                                ds_gt = ds_config['geoT']
                                ds_nd = ds_config['ndv']

                                xpos, ypos = utils._geo2pixel(
                                    this_xyz.x, this_xyz.y, ds_gt, node='pixel'
                                )
                                if xpos < ds_config['nx'] \
                                   and ypos < ds_config['ny'] \
                                   and xpos >=0 \
                                   and ypos >=0:
                                    tgrid = ds_band.ReadAsArray(xpos, ypos, 1, 1)
                                    if tgrid is not None:
                                        if not np.isnan(ds_config['ndv']):
                                            tgrid[tgrid==ds_config['ndv']] = np.nan

                                        if np.isnan(tgrid[0][0]):
                                            if data_mask['invert']:
                                                mask_count += 1
                                                masked = (True & z_masked)
                                            else:
                                                masked = False
                                        else:
                                            if not data_mask['invert']:
                                                mask_count += 1
                                                masked = (True & z_masked)
                                            else:
                                                masked = False
                        else:
                            src_ds = ogr.Open(data_mask['data_mask'])
                            if src_ds is not None:
                                layer = src_ds.GetLayer()
                                geomcol = gdalfun.ogr_union_geom(layer, verbose=False)
                                if not geomcol.IsValid():
                                    geomcol = geomcol.Buffer(0)

                                this_wkt = this_xyz.export_as_wkt()
                                this_point = ogr.CreateGeometryFromWkt(this_wkt)
                                if this_point.Within(geomcol):
                                    if not data_mask['invert']:
                                        masked = (True & z_masked)
                                        mask_count += 1
                                    else:
                                        masked = False
                                else:
                                    if not data_mask['invert']:
                                        masked = False
                                    else:
                                        masked = (True & z_masked)
                                        mask_count += 1

                                src_ds = layer = None

                    if not masked:
                        yield(this_xyz)

                    if mask_count > 0 and self.verbose:
                        utils.echo_msg_bold(
                            f'masked {mask_count} data records from {self.fn}'
                        )       
                    src_ds = ds_band = None

                    
    ## ==============================================
    ## yield the data either as xyz or as binned arrays, the latter of which
    ## depends on self.region, self.x_inc and self.y_inc to be set.
    ## ==============================================
    def yield_xyz(self):
        """Yield the data as xyz points

        incoming data are numpy rec-arrays of x,y,z<w,u> points.

        if 'w' (weight) is not set by the dataset, they will be given the
        default value from `self.weight`, otherwise the incoming weight 
        values will be multiplied by `self.weight`

        if 'u' (uncertainty) is not set by the dataset, they will be given 
        the default value from `self.uncertainty`, otherwise the incoming 
        uncertainty values will be caluculated by 
        `np.sqrt(u**2 + self.uncertainty**2)`
        """
        
        count = 0
        for points in self.transform_and_yield_points():
            try:
                points_w = points['w']
            except:
                points_w = np.ones(points['z'].shape).astype(float)

            points_w *= self.weight if self.weight is not None else 1.
            points_w[np.isnan(points_w)] = 1
                
            try:
                points_u = points['u']
            except:
                points_u = np.zeros(points['z'].shape)
                
            points_u = np.sqrt(
                points_u**2 \
                + (self.uncertainty if self.uncertainty is not None else 0)**2
            )
            points_u[np.isnan(points_u)] = 0
            dataset = np.vstack(
                (points['x'], points['y'], points['z'], points_w, points_u)
            ).transpose()
            count += len(dataset)
            points = None
            for point in dataset:
                this_xyz = xyzfun.XYZPoint(
                    x=point[0], y=point[1], z=point[2], w=point[3], u=point[4]
                )
                yield(this_xyz)

        if self.verbose:
            utils.echo_msg_bold(
                f'parsed {count} data records from {self.fn} @ a weight of {self.weight}'
            )

        
    def yield_array(self, want_sums=True):
        """Yield the data as an array which coincides with the desired 
        region, x_inc and y_inc

        incoming data are numpy rec-arrays of x,y,z<w,u> points

        if 'w' (weight) is not set by the dataset, they will be given the
        default value from `self.weight`, otherwise the incoming weight 
        values will be multiplied by `self.weight`

        if 'u' (uncertainty) is not set by the dataset, they will be given 
        the default value from `self.uncertainty`, otherwise the incoming 
        uncertainty values will be caluculated by 
        `np.sqrt(u**2 + self.uncertainty**2)`

        output arrays are weighted sums; get the actual values with:
        weight = weight / count
        z = (z / weight) / count
        x = (x / weight) / count
        y = (y / weight) / count

        """

        count = 0
        for points in self.transform_and_yield_points():
            count += points.size
            xcount, ycount, _ = self.region.geo_transform(
                x_inc=self.x_inc, y_inc=self.y_inc, node='grid'
            )
            
            point_array = pointz.PointPixels(
                src_region=self.region,
                x_size=xcount,
                y_size=ycount,
                verbose=self.verbose
            )
            yield(
                point_array(
                    points,
                    weight=self.weight,
                    uncertainty=self.uncertainty,
                    mode='sums' if want_sums else 'mean'
                )
            )

        if self.verbose:
            utils.echo_msg_bold(
                f'parsed {count} data records from {self.fn} @ a weight of {self.weight}'
            )    

    def stacks(self, out_name=None):
        from . import globato
        
        gbt = globato.GdalRasterStacker(
            region=self.region,
            x_inc=self.x_inc,
            y_inc=self.y_inc,
            dst_srs=self.dst_srs,
            cache_dir=self.cache_dir
        )

        stacked_fn = gbt.process_stack(self.parse(), out_name=out_name)

        return stacked_fn

            
    def stacks_yield_xyz(self, out_name=None, ndv=-9999, fmt='GTiff'):
        """yield the result of `stacks` as an xyz object"""
        
        #stacked_fn = self.stacks(out_name=out_name, ndv=ndv)
        from . import globato
        
        gbt = globato.GdalRasterStacker(
            region=self.region,
            x_inc=self.x_inc,
            y_inc=self.y_inc,
            dst_srs=self.dst_srs,
            cache_dir=self.cache_dir
        )

        stacked_fn = gbt.process_stack(self.parse())

        #, fmt=fmt)#, mode=mode)
        sds = gdal.Open(stacked_fn)
        sds_gt = sds.GetGeoTransform()
        sds_z_band = sds.GetRasterBand(1) # the z band from stacks
        sds_w_band = sds.GetRasterBand(3) # the weight band from stacks
        sds_u_band = sds.GetRasterBand(4) # the uncertainty band from stacks
        sds_x_band = sds.GetRasterBand(6) # the x band from stacks
        sds_y_band = sds.GetRasterBand(7) # the y band from stacks
        srcwin = (0, 0, sds.RasterXSize, sds.RasterYSize)
        for y in range(srcwin[1], srcwin[1] + srcwin[3], 1):
            sz = sds_z_band.ReadAsArray(srcwin[0], y, srcwin[2], 1)
            ## skip row if all values are ndv
            if np.all(sz == ndv):
                continue
            
            sw = sds_w_band.ReadAsArray(srcwin[0], y, srcwin[2], 1)
            su = sds_u_band.ReadAsArray(srcwin[0], y, srcwin[2], 1)
            sx = sds_x_band.ReadAsArray(srcwin[0], y, srcwin[2], 1)
            sy = sds_y_band.ReadAsArray(srcwin[0], y, srcwin[2], 1)
            for x in range(0, sds.RasterXSize):
                z = sz[0, x]
                if z != ndv:
                    if self.stack_node: # yield avg x/y data instead of center
                        out_xyz = xyzfun.XYZPoint(
                            x = sx[0, x], y = sy[0, x], z = z,
                            w = sw[0, x], u = su[0, x]
                        )
                    else: # yield center of pixel as x/y
                        geo_x, geo_y = utils._pixel2geo(x, y, sds_gt)
                        out_xyz = xyzfun.XYZPoint(
                            x=geo_x, y=geo_y, z=z, w=sw[0, x], u=su[0, x]
                        )
                        
                    yield(out_xyz)
                    
            sw = su = sx = sy = None
                    
        sds = sds_z_band = sds_w_band = sds_u_band \
            = sds_x_band = sds_y_band = None

        
    def blocks_yield_xyz(self, out_name=None):
        """yield the result of `stacks` as an xyz object"""

        #stacked_fn = self._blocks(out_name=out_name)
        from . import globato
        
        gbt = globato.GlobatoStacker(
            region=self.region,
            x_inc=self.x_inc,
            y_inc=self.y_inc,
            dst_srs=self.dst_srs,
            cache_dir=self.cache_dir
        )

        stacked_fn = gbt.process_blocks(self.parse())        
        sds = h5.File(stacked_fn, 'r')
        sds_gt = [float(x) for x in sds['crs'].attrs['GeoTransform'].split()]
        sds_stack = sds['block']
        #sds_stack = sds['sums']
        stack_shape = sds['block']['z'].shape
        srcwin = (0, 0, stack_shape[1], stack_shape[0])
        for y in range(srcwin[1], srcwin[1] + srcwin[3], 1):
            sz = sds_stack['z'][y:y+1, srcwin[0]:srcwin[0]+srcwin[2]]
            ## skip row if all values are ndv
            if np.all(np.isnan(sz)):
                continue

            sw = sds_stack['weight'][y:y+1, srcwin[0]:srcwin[0]+srcwin[2]]
            su = sds_stack['uncertainty'][y:y+1, srcwin[0]:srcwin[0]+srcwin[2]]
            sc = sds_stack['count'][y:y+1, srcwin[0]:srcwin[0]+srcwin[2]]

            sx = sds_stack['x'][y:y+1, srcwin[0]:srcwin[0]+srcwin[2]]
            sy = sds_stack['y'][y:y+1, srcwin[0]:srcwin[0]+srcwin[2]]
            for x in range(0, stack_shape[1]):
                z = sz[0, x]
                if not np.isnan(z):
                    if self.stack_node: # yield avg x/y data instead of center
                        out_xyz = xyzfun.XYZPoint(
                            x=sx[0, x], y=sy[0, x], z=z, w=sw[0, x], u=su[0, x]
                        )
                    else: # yield center of pixel as x/y
                        geo_x, geo_y = utils._pixel2geo(x, y, sds_gt)
                        out_xyz = xyzfun.XYZPoint(
                            x=geo_x, y=geo_y, z=z, w=sw[0, x], u=su[0, x]
                        )
                        
                    yield(out_xyz)
        sds.close()


    ## ==============================================
    ## Data Dump/Export/Archive
    ## ==============================================
    def _xyz_dump(self, this_xyz, dst_port=sys.stdout, encode=False):
        this_xyz.dump(
            include_w=True if self.weight is not None else False,
            include_u=True if self.uncertainty is not None else False,
            dst_port=dst_port,
            encode=encode,
            precision=self.dump_precision
        )

        
    def dump_xyz(self, dst_port=sys.stdout, encode=False):
        """dump the XYZ data from the dataset.

        data gets parsed through `self.xyz_yield`. 
        See `set_yield` for more info.
        """

        for this_xyz in self.xyz_yield:
            self._xyz_dump(this_xyz, dst_port=dst_port, encode=encode)

            
    def dump_xyz_direct(self, dst_port=sys.stdout, encode=False):
        """dump the XYZ data from the dataset

        data get dumped directly from `self.yield_xyz`, 
        by-passing `self.xyz_yield`.
        """
        
        for this_xyz in self.yield_xyz():
            self._xyz_dump(this_xyz, dst_port=dst_port, encode=encode)


    def dump_positions(self, dst_port=sys.stdout, encode=False):
        for arrs, srcwin, gt in this_entry.array_yield:            

            #np.vstack(arrs['pixel_x'], arrs['pixel_y'], arrs['z'])
            #l = f'{}'
            l = 'non'
            dst_port.write(l.encode('utf-8') if encode else l)
            # this_xyz.dump(
            #     include_w=True if self.weight is not None else False,
            #     include_u=True if self.uncertainty is not None else False,
            #     dst_port=dst_port,
            #     encode=encode,
            #     precision=self.dump_precision
            # )

            
    ## ==============================================            
    ## Data export
    ## ==============================================
    def export_xyz_as_list(self, z_only = False):
        """return the XYZ data from the dataset as python list

        This may get very large, depending on the input data.
        """

        return([xyz.z if z_only else xyz.copy() for xyz in self.xyz_yield])

    
    def export_xyz_as_ogr(self):
        """Make a point vector OGR DataSet Object from src_xyz

        for use in gdal gridding functions
        """

        dst_ogr = self.metadata['name']
        ogr_ds = gdal.GetDriverByName('Memory').Create(
            '', 0, 0, 0, gdal.GDT_Unknown
        )
        layer = ogr_ds.CreateLayer(
            dst_ogr,
            geom_type=ogr.wkbPoint25D
        )
        for x in ['long', 'lat', 'elev', 'weight']:
            fd = ogr.FieldDefn(x, ogr.OFTReal)
            fd.SetWidth(12)
            fd.SetPrecision(8)
            layer.CreateField(fd)
            
        f = ogr.Feature(feature_def=layer.GetLayerDefn())        
        for this_xyz in self.xyz_yield:
            f.SetField(0, this_xyz.x)
            f.SetField(1, this_xyz.y)
            f.SetField(2, this_xyz.z)
            f.SetField(3, this_xyz.w)
                
            wkt = this_xyz.export_as_wkt(include_z=True)
            g = ogr.CreateGeometryFromWkt(wkt)
            f.SetGeometryDirectly(g)
            layer.CreateFeature(f)
            
        return(ogr_ds)

    
    def export_xyz_as_pandas(self):
        """export the point data as a pandas dataframe.
        
        This may get very large, depending on the input data.
        """

        import pandas as pd
        
        dataset = pd.DataFrame(
            {}, columns=['x', 'y', 'z', 'weight', 'uncertainty']
        )
        for points in self.transform_and_yield_points():
            try:
                points_w = points['w']
            except:
                points_w = np.ones(points['z'].shape)

            points_w *= self.weight if self.weight is not None else 1
            points_w[np.isnan(points_w)] = 1
                
            try:
                points_u = points['u']
            except:
                points_u = np.zeros(points['z'].shape)

            points_u = np.sqrt(
                points_u**2 \
                + (self.uncertainty if self.uncertainty is not None else 0)**2
            )
            points_u[np.isnan(points_u)] = 0
            points_dataset = pd.DataFrame(
                {'x': points['x'],
                 'y': points['y'],
                 'z': points['z'],
                 'weight': points_w,
                 'uncertainty': points_u,
                },
                columns=['x', 'y', 'z', 'weight', 'uncertainty']
            )
            dataset = pd.concat(
                [dataset, points_dataset], ignore_index=True
            )

        return(dataset)

    
    ## TODO: fix/write these functions
    def export_stack_as_pandas(self):
        pass

    
    def export_stack_as_gdal(self, fmt='GTiff'):
        """export the stack h5 to gdal multi-band or separate file(s)
        """
        
        num_bands = len(stack_grp)
        driver = gdal.GetDriverByName(fmt)
        stack_dataset = driver.Create(
            stack_fn, xcount, ycount, num_bands, gdal.GDT_Float32,
            options=['COMPRESS=LZW',
                     'PREDICTOR=2',
                     'TILED=YES',
                     'BIGTIFF=YES']
        )
        stack_dataset.SetGeoTransform(dst_gt)
        if self.dst_srs is not None:
            stack_dataset.SetProjection(gdalfun.osr_wkt(self.dst_srs))

        with utils.ccp(
                total=len(stack_grp.keys()),
                desc='converting CSG stack to GeoTIFF',
                leave=self.verbose
        ) as pbar:
            for i, key in enumerate(reversed(stack_grp.keys())):
                pbar.update()
                stack_grp[key][np.isnan(stack_grp[key])] = ndv
                stack_band = stack_dataset.GetRasterBand(i + 1)
                stack_band.WriteArray(stack_grp[key][...])
                stack_band.SetDescription(key)
                stack_band.SetNoDataValue(ndv)
                stack_band.FlushCache()

        stack_dataset = None

        ## create the mask geotiff
        num_bands, mask_keys = h5_get_datasets_from_grp(mask_grp)
        driver = gdal.GetDriverByName(fmt)
        mask_dataset = driver.Create(
            mask_fn, xcount, ycount, num_bands, gdal.GDT_Byte,
            options=['COMPRESS=LZW', 'BIGTIFF=YES'] \
            if fmt == 'GTiff' \
            else ['COMPRESS=LZW']
        )
        mask_dataset.SetGeoTransform(dst_gt)
        if self.dst_srs is not None:
            mask_dataset.SetProjection(gdalfun.osr_wkt(self.dst_srs))

        with utils.ccp(
                total=len(mask_keys),
                desc='converting CSG mask to GeoTiff',
                leave=self.verbose
        ) as pbar:
            ii = 0
            for i, key in enumerate(mask_keys):
                pbar.update()
                srcwin = (0, 0, xcount, ycount)
                for y in range(
                        srcwin[1], srcwin[1] + srcwin[3], 1
                ):
                    mask_grp[key][
                        srcwin[1]:srcwin[1] + srcwin[3],
                        srcwin[0]:srcwin[0] + srcwin[2]
                    ][np.isnan(mask_grp[key][
                        srcwin[1]:srcwin[1] + srcwin[3],
                        srcwin[0]:srcwin[0] + srcwin[2]
                    ])] = 0
                    #if np.all(mask_grp[key][mask_grp[key] == 0]):
                    #    continue

                    #ii += 1
                    mask_band = mask_dataset.GetRasterBand(i+1)
                    mask_band.WriteArray(
                        mask_grp[key][srcwin[1]:srcwin[1] + srcwin[3],
                                      srcwin[0]:srcwin[0] + srcwin[2]]
                    )
                    md = {}
                    for attrs_key in mask_grp[key].attrs.keys():
                        md[attrs_key] = mask_grp[key].attrs[attrs_key]

                    if 'name' in md.keys():
                        mask_band.SetDescription(md['name'])
                    else:
                        mask_band.SetDescription(mask_grp[key].name)

                    mask_band.SetNoDataValue(0)
                    mask_band.SetMetadata(md)
                    mask_band.FlushCache()

        mask_dataset = None

    ## ==============================================
    ## Data archive
    ## ==============================================
    def _archive_xyz_test(self, **kwargs):
        srs_all = []
        a_name = None
        if 'dirname' in kwargs.keys() and kwargs['dirname'] is not None:
            a_name = kwargs['dirname']
        #else:
        aa_name = self.metadata['name'].split('/')[0]
        if self.region is None:
            aa_name = '{}_{}'.format(aa_name, utils.this_year())
        else:
            aa_name = '{}_{}_{}'.format(
                aa_name, self.region.format('fn'), utils.this_year())

        if a_name is None:
            a_name = aa_name
            self.archive_datalist = '{}.datalist'.format(a_name)
        else:
            self.archive_datalist = '{}.datalist'.format(a_name)
            a_name = os.path.join(a_name, aa_name)

        utils.echo_msg(self.archive_datalist)
        if not os.path.exists(os.path.dirname(self.archive_datalist)):
            try:
                os.makedirs(os.path.dirname(self.archive_datalist))
            except:
                pass

        archive_keys = []
        for this_entry in self.parse():
            #utils.echo_msg(this_entry)
            datalist_dirname = os.path.join(
                a_name, os.path.dirname(this_entry.metadata['name'])
            )
            this_key = datalist_dirname.split('/')[-1]
            if this_key not in archive_keys:
                archive_keys.append(this_key)
                utils.echo_msg_bold(this_key)
                metadata=this_entry.format_metadata()
                utils.echo_msg_bold(metadata)
            #if not os.path.exists(datalist_dirname):
            #    os.makedirs(datalist_dirname)

            
    def archive_xyz(self, **kwargs):
        """Archive data from the dataset to XYZ in the given dataset region.
        
        will convert all data to XYZ within the given region and will arrange
        the data as a datalist based on inputs.

        Data comes from `self.xyz_yield` set in `self.set_yield`. 
        So will process data through `stacks` before archival if region, 
        x_inc and y_inc are set.
        """
        
        srs_all = []
        a_name = None
        if 'dirname' in kwargs.keys() and kwargs['dirname'] is not None:
            a_name = kwargs['dirname']
        #else:
        aa_name = self.metadata['name'].split('/')[0]
        if self.region is None:
            aa_name = '{}_{}'.format(aa_name, utils.this_year())
        else:
            aa_name = '{}_{}_{}'.format(
                aa_name, self.region.format('fn'), utils.this_year())

        if a_name is None:
            a_name = aa_name
            self.archive_datalist = '{}.datalist'.format(a_name)
        else:
            self.archive_datalist = '{}.datalist'.format(a_name)
            a_name = os.path.join(a_name, aa_name)

        utils.echo_msg(self.archive_datalist)
        if not os.path.exists(os.path.dirname(self.archive_datalist)):
            try:
                os.makedirs(os.path.dirname(self.archive_datalist))
            except:
                pass

        archive_keys = []
        for this_entry in self.parse():
            srs_all.append(this_entry.dst_srs \
                           if this_entry.dst_srs is not None \
                           else this_entry.src_srs)
            datalist_dirname = os.path.join(
                a_name, os.path.dirname(this_entry.metadata['name'])
            )
            this_key = datalist_dirname.split('/')[-1]
            if not os.path.exists(datalist_dirname):
                os.makedirs(datalist_dirname)

            sub_datalist = os.path.join(datalist_dirname, f'{this_key}.datalist')
            if utils.fn_basename2(os.path.basename(this_entry.fn)) == '':
                sub_xyz_path = '.'.join(
                    [utils.fn_basename2(os.path.basename(this_entry.metadata['name'])),
                     'xyz']
                )
            else:
                sub_xyz_path = '.'.join(
                    [utils.fn_basename2(os.path.basename(this_entry.fn)),
                     'xyz']
                )

            this_xyz_path = os.path.join(datalist_dirname, sub_xyz_path)
            if os.path.exists(this_xyz_path):
                utils.echo_warning_msg(
                    f'{this_xyz_path} already exists, skipping...'
                )
                continue

            with open(sub_datalist, 'a') as sub_dlf:
                with open(this_xyz_path, 'w') as xp:
                    # data will be processed independently of each other
                    for this_xyz in this_entry.xyz_yield: 
                        this_xyz.dump(
                            include_w=True if self.weight is not None else False,
                            include_u=True if self.uncertainty is not None else False,
                            dst_port=xp, encode=False, precision=self.dump_precision
                        )

                if os.stat(this_xyz_path).st_size > 0:
                    sub_dlf.write(f'{sub_xyz_path} 168 1 0\n')
                else:
                    utils.remove_glob(f'{this_xyz_path}*')

            if os.stat(sub_datalist).st_size == 0:
                utils.remove_glob(sub_datalist)
            else:
                with open(self.archive_datalist, 'a') as dlf:
                    if not this_key in archive_keys:
                        archive_keys.append(this_key)                    
                        dlf.write(
                            '{name}.datalist -1 {weight} {uncertainty} {metadata}\n'.format(
                                name=os.path.join(datalist_dirname, this_key),
                                weight=utils.float_or(this_entry.weight, 1),
                                uncertainty=utils.float_or(this_entry.uncertainty, 0),
                                metadata=this_entry.format_metadata()
                            )
                        )

            if not os.listdir(datalist_dirname):
                os.rmdir(datalist_dirname)
                    
        ## generate datalist inf/json
        srs_set = set(srs_all)
        if len(srs_set) == 1:
            arch_srs = srs_set.pop()
        else:
            arch_srs = None

        utils.remove_glob(f'{self.archive_datalist}.*')
        this_archive = DatasetFactory(
            mod=self.archive_datalist,
            data_format=-1,
            parent=None,
            weight=1,
            uncertainty=0,
            src_srs=arch_srs,
            dst_srs=None,
            cache_dir=self.cache_dir,
        )._acquire_module().initialize()
        
        return(this_archive.inf())

    ## ==============================================
    ## Fetching
    ## ==============================================
    def fetch_data(self, fetches_module, check_size=True):
        this_fetches = fetches.FetchesFactory(
            mod=fetches_module,
            src_region=self.region,
            verbose=self.verbose,
            outdir=self.cache_dir,
            callback=fetches.fetches_callback
        )._acquire_module()        
        this_fetches.run()
        fr = fetches.fetch_results(this_fetches, check_size=check_size)
        fr.daemon = True
        fr.start()
        fr.join()

        return(fr)

                
class CUDEMFile(ElevationDataset):
    """CUDEM netcdf raster

    the cudem netcdf/h5 contains uninterpolated elevation data, 
    uncertainty weight and data mask...
    """
    
    def __init__(self, stack=False, **kwargs):
        super().__init__(**kwargs)
        self.stack = stack


    def yield_points_nc(self):

        ## extract z, unc and weight grids
        ## process through gdalfile
        
        nc_data = nc.Dataset(self.fn, 'r')
        stack_lat = nc_data['/stack/latitude'][...,]
        stack_lon = nc_data['/stack/longitude'][...,]
        stack_h = nc_data['/stack/stack_h'][...,]

        stack_u = nc_data['/stack/uncertainty'][...,]
        stack_w = nc_data['/stack/weight'][...,]
        counts = nc_data['/stack/count'][...,]

        nc_data.close()

        dataset = np.column_stack(
            (stack_lon.data, stack_lat.data,
             stack_h.data[0], stack_w.data[0],
             stack_u.data[0])
        )
        points = np.rec.fromrecords(
            dataset, names='x, y, z, w, u'
        )
        #points = points[points.z != -9999]
        
        yield(points)

    def yield_points(self):
        data = h5.File(self.fn)
        stack_lat = data['/stack/y'][...,]
        stack_lon = data['/stack/x'][...,]
        stack_h = data['/stack/z'][...,]

        stack_u = data['/stack/uncertainty'][...,]
        stack_w = data['/stack/weights'][...,]
        counts = data['/stack/count'][...,]

        data.close()

        nan_m = np.isnan(stack_h)
        stack_h = stack_h[~nan_m].flatten()                        
        stack_lat = stack_lat[~nan_m].flatten()                        
        stack_lon = stack_lon[~nan_m].flatten()                        
        stack_u = stack_u[~nan_m].flatten()                        
        stack_w = stack_w[~nan_m].flatten()                        
        
        dataset = np.column_stack(
            (stack_lon, stack_lat,
             stack_h, stack_w,
             stack_u)
        )
        points = np.rec.fromrecords(
            dataset, names='x, y, z, w, u'
        )
        #points = points[points.z != -9999]
        
        yield(points)
                

class FactoryDatalists(ElevationDataset):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        
class CRMDatalist(FactoryDatalists):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)    


## ==============================================
## Datasets Factory - module parser
## ==============================================
class DatasetFactory(factory.CUDEMFactory):
    """Dataset Factory Settings and Generator
    
    Parse a datalist entry and return the dataset object
    """

    from . import xyzfile
    from . import lasfile
    from . import gdalfile
    from . import mbsfile
    from . import ogrfile
    from . import icesat2file
    from . import swotfile
    from . import ziplists
    from . import fetchers
    from . import datalists
    
    _modules = {
        ## negative values are `datalists`, where they contain various other datasets.
        -1: {'name': 'datalist',
             'fmts': ['datalist', 'mb-1', 'dl'],
             'description': ('An extended MB-System style datalist containting '
                             'dlim-compatible datasets'),
             'call': datalists.Datalist},
        -2: {'name': 'zip',
             'fmts': ['zip', 'ZIP'],
             'description': 'A zipfile containing dlim-compatible datasets',
             'call': ziplists.ZIPlist}, # add other archive formats (gz, tar.gz, 7z, etc.)
        -3: {'name': 'scratch',
             'fmts': [''],
             'description': ('A scratch dataset, including a python list '
                             'of dlim-compatible datasets'),
             'call': datalists.Scratch},
        -4: {'name': 'points',
             'fmts': [''],
             'description': ('A points dataset, a numpy rec-array or a '
                             'pandas dataframe, defining x, y and z data columns'),
             'call': datalists.Points},
        ## data files
        167: {'name': 'yxz',
              'fmts': ['yxz'],
              'description': 'ascii DSV datafile formatted as y,x,z',
              'call': xyzfile.YXZFile},
        168: {'name': 'xyz',
              'fmts': ['xyz', 'csv', 'dat', 'ascii', 'txt', 'XYZ'],
              'description': 'An ascii DSV datafile formatted as x,y,z',
              'call': xyzfile.XYZFile},
        200: {'name': 'gdal',
              'fmts': ['tif', 'tiff', 'img', 'grd', 'nc', 'vrt'],
              'description': 'A gdal-compatible raster dataset',
              'call': gdalfile.GDALFile},
        201: {'name': 'bag',
              'fmts': ['bag'],
              'description': 'A BAG bathymetry dataset',
              'call': gdalfile.BAGFile},
        202: {'name': 'swot_pixc',
              'fmts': ['h5'],
              'description': 'An HDF5 SWOT PIXC datafile',
              'call': swotfile.SWOT_PIXC},
        203: {'name': 'swot_hr_raster',
              'fmts': ['nc'],
              'description': 'An HDF5 SWOT HR Raster datafile',
              'call': swotfile.SWOT_HR_Raster},
        300: {'name': 'las',
              'fmts': ['las', 'laz'],
              'description': 'An las or laz lidar datafile',
              'call': lasfile.LASFile},
        301: {'name': 'mbs',
              'fmts': ['fbt', 'mb'],
              'description': 'An MB-System-compatible multibeam datafile',
              'call': mbsfile.MBSParser},
        302: {'name': 'ogr',
              'fmts': ['000', 'shp', 'geojson', 'gpkg', 'gdb/'],
              'description': 'An ogr-compatible vector datafile',
              'call': ogrfile.OGRFile},
        303: {'name': 'icesat2_atl03',
              'fmts': ['h5'],
              'description': 'An HDF5 IceSat2 ATL03 datafile',
              'call': icesat2file.IceSat2_ATL03},
        304: {'name': 'icesat2_atl24',
              'fmts': ['h5'],
              'description': 'An HDF5 IceSat2 ATL24 datafile',
              'call': icesat2file.IceSat2_ATL24},
        310: {'name': 'cudem',
              'fmts': ['csg', 'nc', 'h5'],
              'description': 'A netCDF/h5 CUDEM file',
              'call': CUDEMFile},
        ## fetches modules
        -100: {'name': 'https',
               'fmts': ['https'],
               'description': 'A URL pointing to a dlim-compatible datafile',
               'call': fetchers.Fetcher},
        -101: {'name': 'gmrt',
               'fmts': ['gmrt'],
               'description': 'The GMRT fetches module',
               'call': fetchers.GMRTFetcher},
        -102: {'name': 'gebco',
               'fmts': ['gebco'],
               'description': 'The GEBCO fetches module',
               'call': fetchers.GEBCOFetcher},
        -103: {'name': 'copernicus',
               'fmts': ['copernicus'],
               'description': 'The Copernicus fetches module',
               'call': fetchers.CopernicusFetcher},
        -104: {'name': 'fabdem',
               'fmts': ['fabdem'],
               'description': 'The FABDEM fetches module',
               'call': fetchers.FABDEMFetcher},
        -105: {'name': 'nasadem',
               'fmts': ['nasadem'],
               'description': 'The NASADEM fetches module',
               'call': fetchers.Fetcher},
        -106: {'name': 'mar_grav',
               'fmts': ['mar_grav'],
               'description': 'The mar_grav fetches module',
               'call': fetchers.MarGravFetcher},
        -107: {'name': 'srtm_plus',
               'fmts': ['srtm_plus'],
               'description': 'The srtm_plus fetches module',
               'call': fetchers.Fetcher},
        -108: {'name': 'synbath',
               'fmts': ['synbath'],
               'description': 'The SynBath fetches module',
               'call': fetchers.Fetcher},
        -109: {'name': 'gedtm30',
               'fmts': ['gedtm30'],
               'description': 'Global DTM fetches module',
               'call': fetchers.GEDTM30Fetcher},
        -110: {'name': "swot",
               'fmts': ['swot'],
               'description': '	The SWOT fetches module',
               'call': fetchers.SWOTFetcher},
        -111: {'name': "icesat2",
               'fmts': ['icesat2'],
               'description': 'The IceSat2 fetches module',
               'call': fetchers.IceSat2Fetcher},
        -200: {'name': 'charts',
               'fmts': ['charts'],
               'description': 'The charts fetches module',
               'call': fetchers.ChartsFetcher},
        -201: {'name': 'multibeam',
               'fmts': ['multibeam'],
               'description': 'The multibeam fetches module',
               'call': fetchers.MBSFetcher},
        -202: {'name': 'hydronos',
               'fmts': ['hydronos'],
               'description': 'The hydronos fetches module',
               'call': fetchers.HydroNOSFetcher},
        -203: {'name': 'ehydro',
               'fmts': ['ehydro'],
               'description': 'The ehydro fetches module',
               'call': fetchers.eHydroFetcher},
        -204: {'name': 'bluetopo',
               'fmts': ['bluetopo'],
               'description': 'The bluetopo fetches module',
               'call': fetchers.BlueTopoFetcher},
        -205: {'name': 'ngs',
               'fmts': ['ngs'],
               'description': 'The ngs fetches module',
               'call': fetchers.NGSFetcher},
        -206: {'name': 'tides',
               'fmts': ['tides'],
               'description': 'The tides fetches module',
               'call': fetchers.TidesFetcher},
        -207: {'name': 'digital_coast',
               'fmts': ['digital_coast'],
               'description': 'The digital_coast fetches module',
               'call': fetchers.Fetcher},
        -208: {'name': 'ncei_thredds',
               'fmts': ['ncei_thredds'],
               'description': 'The ncei_thredds fetches module',
               'call': fetchers.Fetcher},
        -209: {'name': 'tnm',
               'fmts': ['tnm'],
               'description': 'The TNM fetches module',
               'call': fetchers.Fetcher},
        -210: {'name': "CUDEM",
               'fmts': ['CUDEM'],
               'description': 'The CUDEM fetches module',
               'call': fetchers.Fetcher},
        -211: {'name': "CoNED",
               'fmts': ['CoNED'],
               'description': 'The CoNED fetches module',
               'call': fetchers.DAVFetcher_CoNED},
        -212: {'name': "SLR",
               'fmts': ['SLR'],
               'description': '	The SLR fetches module',
               'call': fetchers.DAVFetcher_SLR},
        -213: {'name': 'waterservies',
               'fmts': ['waterservices'],
               'description': 'The waterservices fetches module',
               'call': fetchers.WaterServicesFetcher},
        -214: {'name': 'ned',
               'fmts': ['ned', 'ned1'],
               'description': 'The NED fetches module',
               'call': fetchers.NEDFetcher},        
        -215: {'name': "csb",
               'fmts': ['csb'],
               'description': 'The CSB fetches module',
               'call': fetchers.Fetcher},
        -216: {'name': 'wa_dnr',
               'fmts': ['wa_dnr'],
               'description': 'The Washington DNR lidar portal',
               'call': fetchers.DNRFetcher},
        -217: {'name': 'r2r',
               'fmts': ['r2r'],
               'description': 'The r2r fetches module',
               'call': fetchers.R2RFetcher},
        -300: {'name': 'emodnet',
               'fmts': ['emodnet'],
               'description': 'The emodnet fetches module',
               'call': fetchers.EMODNetFetcher},
        -301: {'name': 'chs',
               'fmts': ['chs'],
               'description': 'The chs fetches module',
               'call': fetchers.Fetcher}, 
        -302: {'name': 'hrdem',
               'fmts': ['hrdem'],
               'description': '	The hrdem fetches module',
               'call': fetchers.HRDEMFetcher},
        -303: {'name': 'arcticdem',
               'fmts': ['arcticdem'],
               'description': 'The arcticdem fetches module',
               'call': fetchers.Fetcher},
        -304: {'name': 'mrdem',
               'fmts': ['mrdem'],
               'description': '	The mrdem fetches module',
               'call': fetchers.Fetcher},
        -500: {'name': 'vdatum',
               'fmts': ['vdatum'],
               'description': 'The vdatum fetches module',
               'call': fetchers.VDatumFetcher},
        -600: {'name': 'hydrolakes',
               'fmts': ['hydrolakes'],
               'description': 'The hydrolakes fetches module',
               'call': fetchers.HydroLakesFetcher},        
    }
    _datalist_cols = ['path', 'format', 'weight', 'uncertainty', 'title', 'source',
                      'date', 'type', 'resolution', 'horz', 'vert',
                      'url']

    _metadata_keys = ['name', 'title', 'source', 'date', 'data_type', 'resolution',
                      'hdatum', 'vdatum', 'url']

    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    ## ==============================================
    ## redefine the factory default _parse_mod function for datasets
    ## the mod in this case is a datalist entry and the format key
    ## becomes the module
    ## TODO: use csv module to parse
    ## ==============================================
    def _parse_mod(self, mod=None):
        """parse the datalist entry line"""
        
        self.kwargs['fn'] = mod
        if self.kwargs['fn'] is None:
            return(self)

        ## mod exists as a file, no other entry items should occur, so
        ## guess the format and finish there...
        ## the format number becomes the mod_name
        ## check for specified data format as well
        ## if fn is not a path, parse it as a datalist entry
		## breaks on path with space, e.g. /meda/user/My\ Passport/etc
        this_entry = re.findall(r"(?:\".*?\"|\S)+", self.kwargs['fn'].rstrip())
        
        if os.path.exists(self.kwargs['fn']):
            if os.path.isfile(self.kwargs['fn']) \
               or self.kwargs['fn'].split('.')[-1] in get_factory_exts():
                this_entry = [self.kwargs['fn']]

        try:
            entry = [utils.str_or(x) if n == 0 \
                     else utils.str_or(x, replace_quote=False) if n < 2 \
                     else utils.float_or(x) if n < 3 \
                     else utils.float_or(x) if n < 4 \
                     else utils.str_or(x) \
                     for n, x in enumerate(this_entry)]

            #entry = [None if x == '-' else x for x in entry]
            utils.echo_debug_msg(f'initial parsed entry: {entry}')

        except Exception as e:
            utils.echo_error_msg(
                'could not parse entry {}, {}'.format(
                    self.kwargs['fn'], this_entry
                )
            )
            return(self)

        utils.echo_debug_msg(entry)        
        ## data format - entry[1]
        ## guess format based on fn if not specified otherwise
        ## parse the format for dataset specific opts.
        if len(entry) < 2:
            if self.kwargs.get('data_format') is not None:
                entry.append(self.kwargs.get('data_format'))
            else:
                entry = self.guess_and_insert_fmt(entry)

            if len(entry) < 2:
                utils.echo_error_msg(
                    'could not parse entry {}'.format(self.kwargs['fn'])
                )
                return(self)

        elif entry[1] is None or entry[1] == '-':
            if self.kwargs.get('data_format') is not None:
                entry[1] = self.kwargs['data_format']
            else:
                entry = self.guess_and_insert_fmt(entry)
                        
        ## parse the entry format options
        opts = factory.fmod2dict(str(entry[1]), {})
        utils.echo_debug_msg(opts)
        if '_module' in opts.keys():# and len(opts.keys()) > 1:
            entry[1] = utils.int_or((opts['_module']))
            self.mod_args = {i:opts[i] for i in opts if i!='_module'}

        try:
            assert isinstance(utils.int_or(entry[1]), int)
        except:
            utils.echo_error_msg(
                f'could not parse datalist entry {entry}'
            )
            return(self)

        ## entry is a COG to be read with gdal, append vsicurl to the fn
        if entry[0].startswith('http') and int(entry[1]) == 200:
            entry[0] = '/vsicurl/{}'.format(entry[0])
        
        self.kwargs['data_format'] = int(entry[1])
        self.mod_name = int(entry[1])
        
        ## file-name (or fetches module name) - entry[0]
        ## set to relative path from parent
        ## don't set relative path if 'fn' is a
        ## fetches module (entry[1] < -3)
        if 'parent' not in self.kwargs:
            self.kwargs['parent'] = None

        if self.kwargs['parent'] is None:
            self.kwargs['fn'] = entry[0]
        else:
            if self.mod_name >= -2 \
               and os.path.dirname(self.kwargs['parent'].fn) \
               != os.path.dirname(entry[0]) and \
                   ':' not in entry[0] and \
                   self.kwargs['parent'].data_format >= -2:
                self.kwargs['fn'] = os.path.join(
                    os.path.dirname(self.kwargs['parent'].fn), entry[0]
                )
                self.kwargs['fn'] = os.path.relpath(self.kwargs['fn'])
            else:
                self.kwargs['fn'] = entry[0]

        ## weight - entry[2]
        ## inherit weight of parent
        if len(entry) < 3:
            entry.append(self.set_default_weight())
        elif entry[2] is None:
            entry[2] = self.set_default_weight()

        if 'weight' not in self.kwargs:
            self.kwargs['weight'] = 1

        ## multiply the weight by the weight in the parent
        ## if it exists, otherwise weight is weight
        if self.kwargs['parent'] is not None:
            if self.kwargs['weight'] is not None:
                self.kwargs['weight'] *= entry[2]
        else:
            if self.kwargs['weight'] is not None:
                self.kwargs['weight'] = entry[2]

        ## uncertainty - entry[3]
        ## combine with partent using root sum of squares
        if len(entry) < 4:
            entry.append(self.set_default_uncertainty())
        elif entry[3] is None:
            entry[3] = self.set_default_uncertainty()

        if 'uncertainty' not in self.kwargs:
            self.kwargs['uncertainty'] = 0
            
        if self.kwargs['parent'] is not None:
            if self.kwargs['uncertainty'] is not None:
                self.kwargs['uncertainty'] = math.sqrt(
                    self.kwargs['uncertainty']**2 + entry[3]**2
                )
        else:
            if self.kwargs['uncertainty'] is not None:
                self.kwargs['uncertainty'] = entry[3]
                    
        ## Optional arguments follow, for metadata generation
        if 'metadata' not in self.kwargs:
            self.kwargs['metadata'] = {}

        for key in self._metadata_keys:
            if key not in self.kwargs['metadata'].keys():
                self.kwargs['metadata'][key] = None

        ## inherit metadata from parent, if available
        if self.kwargs['parent'] is not None:
            for key in self._metadata_keys:
                if key in self.kwargs['parent'].metadata:
                    if self.kwargs['metadata'][key] is None or key == 'name':
                        self.kwargs['metadata'][key] = self.kwargs['parent'].metadata[key]

        ## set or append metadata from entry
        for i, key in enumerate(self._metadata_keys):
            if key == 'name':
                if self.kwargs['metadata'][key] is None:
                    self.kwargs['metadata'][key] \
                        = utils.fn_basename2(
                            os.path.basename(self.kwargs['fn'].split(':')[0])
                        )
                else:
                    self.set_metadata_entry(
                        utils.fn_basename2(
                            os.path.basename(self.kwargs['fn'].split(':')[0])
                        ),
                        key, '/'
                    )

            else:
                if len(entry) < i+4:
                    entry.append(self.kwargs['metadata'][key])

                self.set_metadata_entry(entry[i+3], key, ', ')

                if key == 'date':
                    self.kwargs['metadata'][key] \
                        = utils.num_strings_to_range(
                            self.kwargs['metadata'][key], entry[i+3]
                        )
                
        return(self.mod_name, self.mod_args)

    
    def set_metadata_entry(self, entry, metadata_field, join_string = '/'):
        if entry != '-' and str(entry).lower() != 'none':
            if self.kwargs['metadata'][metadata_field] is not None:
                if str(self.kwargs['metadata'][metadata_field]).lower() != str(entry).lower():
                    self.kwargs['metadata'][metadata_field] \
                        = join_string.join(
                            [str(self.kwargs['metadata'][metadata_field]), str(entry)]
                        )
            else:
                self.kwargs['metadata'][metadata_field] = str(entry)

                
    def set_default_weight(self):
        return(1)

    
    def set_default_uncertainty(self):
        return(0)

    
    def guess_data_format(self, fn):
        """guess a data format based on the file-name"""
        
        for key in self._modules.keys():
            if fn.split('.')[-1] in self._modules[key]['fmts']:
                return(key)
            ## hack to accept .mb* mb-system files without having to record every one...
            elif fn.split('.')[-1][:2] in self._modules[key]['fmts']:
                return(key)


    def guess_and_insert_fmt(self, entry):
        if len(entry) >= 1:            
            for key in self._modules.keys():
                ## entry is stdin, set to 168 and break
                if entry[0] is None or entry[0] == '-':
                    if len(entry) > 1:
                        entry[1] = 168
                    else:
                        entry.append(168)
                        
                    break
                
                if entry[0].startswith('http') \
                   or entry[0].startswith('/vsicurl/'):
                    see = 'https'
                    
                else:
                    se = entry[0].split('.')
                    see = se[-1] \
                        if len(se) > 1 \
                           else entry[0].split(":")[0]

                if 'fmts' in self._modules[key].keys():
                    if see in self._modules[key]['fmts']:
                        if len(entry) > 1:
                            entry[1] = int(key)
                        else:
                            entry.append(int(key))
                            
                        break

        return(entry)

            
    def write_parameter_file(self, param_file: str):
        try:
            with open(param_file, 'w') as outfile:
                json.dump(self.__dict__, outfile)
                utils.echo_msg(
                    'New DatasetFactory file written to {}'.format(param_file)
                )
                
        except:
            raise ValueError(
                ('DatasetFactory: Unable to write new parameter '
                 f'file to {param_file}')
            )


## ==============================================
## Command-line Interface (CLI)
## $ dlim
##
## datalists cli
## ==============================================
datalists_usage = lambda: """{cmd} ({dl_version}): DataLists IMproved; 
Process and generate datalists

dlim is the elevation data processing tool using various dataset modules. 
dlim's native dataset format is a "datalist". 
A datalist is similar to an MBSystem datalist; 
it is a space-delineated file containing the following columns:

`data-path data-format data-weight data-uncertainty data-name data-source data-date data-resolution data-type data-horz data-vert data-url`

Minimally, `data-path` (column 1) is all that is needed.

An associated inf and geojson file will be gerenated for each datalist 
while only an associated inf file will be genereated for individual datasets

Parse various dataset types by region/increments and yield data as xyz or array. 
Recursive data-structures, which point to datasets (datalist, zip, fetches, etc), 
are negative format numbers, e.g. -1 for datalist. Fetches modules are <= -100.

usage: {cmd} [ -acdghijnquwEJPRT [ args ] ] DATALIST,FORMAT,WEIGHT,UNCERTAINTY ...

Options:
  -R, --region\t\t\tRestrict processing to the desired REGION 
\t\t\t\tWhere a REGION is xmin/xmax/ymin/ymax[/zmin/zmax[/wmin/wmax/umin/umax]]
\t\t\t\tUse '-' to indicate no bounding range; e.g. -R -/-/-/-/-10/10/1/-/-/-
\t\t\t\tOR an OGR-compatible vector file with regional polygons. 
\t\t\t\tWhere the REGION is /path/to/vector[:zmin/zmax[/wmin/wmax/umin/umax]].
\t\t\t\tIf a vector file is supplied, will use each region found therein.
\t\t\t\tOptionally, append `:pct_buffer=<value>` to buffer the region(s) by a percentage.
  -E, --increment\t\tBlock data to INCREMENT in native units.
\t\t\t\tWhere INCREMENT is x-inc[/y-inc]
  -X, --extend\t\t\tNumber of cells with which to EXTEND the output DEM REGION and a 
\t\t\t\tpercentage to extend the processing REGION.
\t\t\t\tWhere EXTEND is dem-extend(cell-count)[:processing-extend(percentage)]
\t\t\t\te.g. -X6:10 to extend the DEM REGION by 6 cells and the processing region by
\t\t\t\t10 percent of the input REGION.
  -J, --s_srs\t\t\tSet the SOURCE projection.
  -P, --t_srs\t\t\tSet the TARGET projection. (REGION should be in target projection) 
  -D, --cache-dir\t\tCACHE Directory for storing temp and output data.
  -Z, --z-precision\t\tSet the target precision of dumped z values. (default is 4)
  -A, --stack-mode\t\tSet the STACK MODE to 'mean', 'min', 'max', 'mixed' or 'supercede' 
\t\t\t\t(with -E and -R)
  -T, --stack_filter\t\tFILTER the data stack using one or multiple filters. 
\t\t\t\tWhere FILTER is filter-name[:opts] (see `grits --modules` for more information)
\t\t\t\tThe -T switch may be set multiple times to perform multiple filters.
\t\t\t\tAvailable FILTERS: {grits_modules}
  -F, --point_filter\t\tFILTER the POINT data using one or multiple filters. 
\t\t\t\tWhere FILTER is filter-name[:opts] (See {cmd} --point-filters for more information)
\t\t\t\tThe -F switch may be set multiple times to perform multiple filters.
\t\t\t\tAvailable FILTERS: {point_filter_modules}
  -V, --archive\t\t\tArchive the DATALIST to the given REGION[/INCREMENTs].
\t\t\t\tSpecify the name of the archive, if not specified an auto-generated 
\t\t\t\tname will be used.

  -m, --mask\t\t\tMASK the datalist to the given REGION/INCREMENTs
  -s, --spatial-metadata\tGenerate SPATIAL METADATA of the datalist to the given 
\t\t\t\tREGION/INCREMENTs
  -g, --glob\t\t\tGLOB the datasets in the current directory to stdout
  -i, --info\t\t\tGenerate and return an INFO dictionary of the dataset
  -l, --list\t\t\tList the assocated datasets from the datalist
  -w, --weights\t\t\tOutput WEIGHT values along with xyz
  -u, --uncertainties\t\tOutput UNCERTAINTY values along with xyz
  -n, --stack-node\t\tOutput stacked x/y data rather than pixel
  -q, --quiet\t\t\tLower the verbosity to a quiet

  --point-filters\t\tDisplay the POINT FILTER descriptions and usage
  --modules\t\t\tDisplay the datatype descriptions and usage
  --help\t\t\tPrint the usage text
  --version\t\t\tPrint the version information

Datalists and data formats:
  A datalist is a file that contains a number of datalist entries, 
  while an entry is a space-delineated line:
  `path [format weight uncertainty [name source date type resolution hdatum vdatum url]]`

  `path` can also be a supported fetches module (dataset IDs <= -100)

Supported datalist formats (see {cmd} --modules <dataset-key> for more information): 
  {dl_formats}

Examples:

  % {cmd} -R-90/-89/30/31/-100/100 *.tif -l -w > tifs_in_region.datalist
  % {cmd} tifs_in_region.datalist -R -90/-89/30/31 -E 1s > tifs_1s.xyz
  % {cmd} -R my_region.shp my_data.xyz -w -s_srs epsg:4326 -t_srs epsg:3565 > my_data_3565.xyz
  % {cmd} -R my_region.shp -w multibeam --archive
""".format(
    cmd=os.path.basename(sys.argv[0]), 
    dl_version=__version__,
    dl_formats=factory.get_module_name_short_desc(
        DatasetFactory._modules
    ),
    grits_modules=factory.get_module_short_desc(
        grits.grits.GritsFactory._modules
    ),
    point_filter_modules=factory.get_module_short_desc(
        pointz.PointFilterFactory._modules
    )
)

def datalists_cli(argv=sys.argv):
    """run datalists from command-line

See `datalists_usage` for full cli options.
    """

    dls = []
    src_srs = None
    dst_srs = None
    i_regions = []
    these_regions = []
    xy_inc = [None, None]
    extend = 0
    want_weights = False
    want_uncertainties = False
    want_mask = False
    want_inf = False
    want_list = False
    want_glob = False
    want_archive = False
    archive_dirname = None
    these_archives = []
    want_verbose = True
    want_region = False
    want_separate = False
    want_sm = False
    invert_region = False
    z_precision = 4
    stack_fltrs = []
    pnt_fltrs = []
    stack_node = False
    stack_mode = 'mean'
    cache_dir = utils.cudem_cache()
    
    ## parse command line arguments.
    i = 1
    while i < len(argv):
        arg = argv[i]
        if arg == '--region' or arg == '-R':
            i_regions.append(str(argv[i + 1]))
            i = i + 1
        elif arg[:2] == '-R':
            i_regions.append(str(arg[2:]))
        elif arg == '--increment' or arg == '-E':
            xy_inc = argv[i + 1].split('/')
            i = i + 1
        elif arg[:2] == '-E':
            xy_inc = arg[2:].split('/')
        elif arg == '--extend' or arg == '-X':
            extend = utils.int_or(argv[i + 1], 0)
            i += 1
        elif arg[:2] == '-X':
            extend = utils.int_or(argv[2:], 0)
        elif arg == '-s_srs' or arg == '--s_srs' or arg == '-J':
            src_srs = argv[i + 1]
            i = i + 1
        elif arg[:2] == '-J':
            srs_srs = arg[2:]
        elif arg == '-t_srs' or arg == '--t_srs' or arg == '-P':
            dst_srs = argv[i + 1]
            i = i + 1
        elif arg[:2] == '-P':
            dst_srs = arg[2:]
        elif arg == '-z_precision' \
             or arg == '--z-precision' \
             or arg == '-Z':
            z_precision = utils.int_or(argv[i + 1], 4)
            i = i + 1
        elif arg[:2] == '-Z':
            z_precision = utils.int_or(argv[2:], 4)

        elif arg == '-stack-mode' or arg == '--stack-mode' or arg == '-A':
            stack_mode = utils.str_or(argv[i + 1], 'mean')
            i = i + 1
        elif arg[:2] == '-A':
            stack_mode = utils.str_or(arg[2:], 'mean')
            
        elif arg == '-stack-filter' or arg == '--stack-filter' or arg == '-T':
            stack_fltrs.append(argv[i + 1])
            i = i + 1
        elif arg[:2] == '-T':
            stack_fltrs.append(argv[2:])
            
        elif arg == '-point-filter' or arg == '--point-filter' or arg == '-F':
            pnt_fltrs.append(argv[i + 1])
            i = i + 1
        elif arg[:2] == '-F':
            pnt_fltrs.append(argv[2:])

        elif arg == '--archive' or arg == '-V':
            want_archive = True
            dataexts = [
                xs for y in [
                    x['fmts'] for x in list(
                        DatasetFactory._modules.values()
                    )
                ] for xs in y
            ]
            if i+1 < len(argv) and not argv[i + 1].startswith('-'):
                archive_dirname = utils.str_or(argv[i + 1])
                i = i + 1
        elif arg[:2] == '-V':
            want_archive = True
            archive_dirname = utils.str_or(argv[2:])
            
        elif arg == '--cache-dir' or arg == '-D' or arg == '-cache-dir':
            cache_dir = utils.str_or(argv[i + 1], utils.cudem_cache())
            i = i + 1
        elif arg[:2] == '-D': cache_dir = utils.str_or(argv[i + 1], utils.cudem_cache)
            
        elif arg == '--mask' or arg == '-m' or arg == '--want-mask':
            want_mask = True
        elif arg == '--invert_region' or arg == '-v':
            invert_region = True
        # elif arg == '--archive' or arg == '-a':
        #     want_archive = True
        elif arg == '--weights' or arg == '-w':
            want_weights = True
        elif arg == '--uncertainties' or arg == '-u':
            want_uncertainties = True
        elif arg == '--info' or arg == '-i':
            want_inf = True
        elif arg == '--region_inf' or arg == '-r':
            want_region = True
        elif arg == '--list' or arg == '-l':
            want_list = True
        elif arg == '--glob' or arg == '-g':
            want_glob = True
        elif arg == '--spatial-metadata' or arg == '-sm':
            want_sm = True
            want_mask = True
        elif arg == '--stack-node' or arg == '-n':
            stack_node = True
            
        elif arg == '--separate' or arg == '-s':
            want_separate = True

        elif arg == '--modules':
            factory.echo_modules(
                DatasetFactory._modules,
                None if i+1 >= len(argv) else utils.int_or(sys.argv[i+1], str(sys.argv[i+1]))
            )
            sys.exit(0)
        elif arg == '--md_modules':
            factory.echo_modules(
                DatasetFactory._modules,
                None if i+1 >= len(argv) else int(sys.argv[i+1]), True
            )
            sys.exit(0)
        elif arg == '--point-filters':
            factory.echo_modules(
                PointFilterFactory._modules,
                None if i+1 >= len(argv) else utils.int_or(sys.argv[i+1], str(sys.argv[i+1]))
            )
            sys.exit(0)            
        elif arg == '--quiet' or arg == '-q':
            want_verbose = False
        elif arg == '--help' or arg == '-h':
            print(datalists_usage())
            sys.exit(1)
        elif arg == '--version' or arg == '-v':
            print('{}, version {}'.format(
                os.path.basename(sys.argv[0]), __version__)
                  )
            sys.exit(1)
        # elif arg[0] == '-':
        #     print(datalists_usage())
        #     sys.exit(0)
        else:
            dls.append(f'{arg}')#'"{}"'.format(arg)) # FIX THIS!!!
        
        i = i + 1

    if len(xy_inc) < 2:
        xy_inc.append(xy_inc[0])
    elif len(xy_inc) == 0:
        xy_inc = [None, None]

    if want_glob:
        import glob
        for key in DatasetFactory()._modules.keys():
            if key != -1 and key != '_factory':
                for f in DatasetFactory()._modules[key]['fmts']:
                    globs = glob.glob('*.{}'.format(f))
                    [sys.stdout.write(
                        '{}\n'.format(
                            ' '.join(
                                [x, str(key), '1', '0']
                            )
                        )
                    ) for x in globs]
                    
        sys.exit(0)

    #stack_fltrs = [':'.join(f.split('/')) for f in stack_fltrs]
    #pnt_fltrs = [':'.join(f.split('//')) for f in pnt_fltrs]
    if not i_regions: i_regions = [None]
    these_regions = regions.parse_cli_region(i_regions, want_verbose)
    if not these_regions: these_regions = [None]
    for rn, this_region in enumerate(these_regions):
        ## buffer the region by `extend` if xy_inc is set
        ## this effects the output naming of masks/stacks!
        ## do we want this like in waffles where the output name
        ## does not include the -X extend buffer?
        if xy_inc[0] is not None \
           and xy_inc[1] is not None \
           and this_region is not None:
            this_region.buffer(
                x_bv=(utils.str2inc(xy_inc[0])*extend),
                y_bv=(utils.str2inc(xy_inc[1])*extend)
            )


        if len(dls) == 0:
            dls = [sys.stdin]
            
        if len(dls) == 0:
            sys.stderr.write(datalists_usage())
            utils.echo_error_msg('you must specify some type of data')
        else:
            ## intiialze the input data. Treat data from CLI as a datalist.
            this_datalist = init_data(
                dls,
                src_region=this_region,
                src_srs=src_srs,
                dst_srs=dst_srs,
                x_inc=xy_inc[0],
                y_inc=xy_inc[1],
                sample_alg='auto',
                want_weight=want_weights,
                want_uncertainty=want_uncertainties,
                want_verbose=want_verbose,
                want_mask=want_mask,
                want_sm=want_sm,
                invert_region=invert_region,
                cache_dir=cache_dir,
                dump_precision=z_precision,
                pnt_fltrs=pnt_fltrs,
                stack_fltrs=stack_fltrs,
                stack_node=stack_node,
                stack_mode=stack_mode
            )

            if this_datalist is not None and this_datalist.valid_p(
                    fmts=DatasetFactory._modules[this_datalist.data_format]['fmts']
            ):
                this_datalist.initialize()                
                if not want_weights:
                    this_datalist.weight = None
                    
                if not want_uncertainties:
                    this_datalist.uncertainty = None
                    
                if want_inf:
                    # output the datalist inf blob
                    print(this_datalist.inf()) 
                elif want_list:
                    # output each dataset from the datalist
                    this_datalist.echo()
                
                elif want_region:
                    # get the region and warp it if necessary
                    this_inf = this_datalist.inf()
                    this_region = regions.Region().from_list(this_inf.minmax)
                    if dst_srs is not None:
                        if src_srs is not None:
                            this_region.src_srs = src_srs
                            this_region.warp(dst_srs)
                        elif this_inf.src_srs is not None:
                            this_region.src_srs = this_inf.src_srs
                            this_region.warp(dst_srs, include_z=False)
                    print(this_region.format('gmt'))
                elif want_archive:
                    # archive the datalist as xyz
                    this_archive = this_datalist.archive_xyz(dirname=archive_dirname) 
                    if this_archive.numpts == 0:
                        utils.remove_glob('{}*'.format(this_archive.name))
                # elif want_tables:
                #     this_datalist._archive_xyz_test()
                else:
                    try:
                        # process and dump each dataset independently
                        if want_separate: 
                            for this_entry in this_datalist.parse():
                                this_entry.dump_xyz()
                        else:
                            # process and dump the datalist as a whole
                            this_datalist.dump_xyz()
                    except KeyboardInterrupt:
                      utils.echo_error_msg('Killed by user')
                      break
                    except BrokenPipeError:
                      utils.echo_error_msg('Pipe Broken')
                      break
                    except Exception as e:
                      utils.echo_error_msg(e)
                      print(traceback.format_exc())

### End                      
