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
## A datalist is similar to an MBSystem datalist; it is a space-delineated file containing the following columns:
## data-path data-format data-weight data-uncertainty data-name data-source data-date data-resolution data-type data-horz data-vert data-url
## Minimally, data-path (column 1) is all that is needed.
##
## An associated inf and geojson file will be gerenated for each datalist
## only an associated inf file will be genereated for individual datasets
##
## Parse various dataset types by region/increments and yield data as xyz or array
## recursive data-structures which point to datasets (datalist, zip, fetches, etc) are negative format numbers, e.g. -1 for datalist
##
## supported datasets include: xyz, gdal, ogr, las/laz (laspy), mbs (MBSystem), fetches (see cudem.fetches)
##
## Initialize a datalist/dataset using init_data(list-of-datalist-entries) where the list of datalist entries can
## be any of the supported dataset formats. init_data will combine all the input datafiles into an internal scratch
## datalist and process that.
##
## If region, x_inc, and y_inc are set, all processing will go through Dataset._stacks() where the data will be combined
## either using the 'supercede', 'weighted-mean', 'min' or 'max' method. Dataset._stacks will output a multi-banded gdal file with the following
## bands: 1: z, 2: count, 3: weight, 4: uncerainty, 5: source-uncertainty
##
## If want_mask is set, _stacks() will also generate a multi-band gdal raster file where each
## mask band contains the mask (0/1) of a specific dataset/datalist, depending on the input data. For a datalist, each band will
## contain a mask for each of the top-level datasets. If want_sm is also set, the multi-band mask
## will be processed to an OGR supported vector, with a feature for each band and the metadata items (cols > 4 in the datalist entry) as fields.
##
## Transform data between horizontal/vertical projections/datums by setting src_srs and dst_srs as 'EPSG:<>'
## if src_srs is not set, but dst_srs is, dlim will attempt to obtain the source srs from the data file itself
## or its respective inf file; otherwise, it will be assumed the source data file is in the same srs as dst_srs
##
## A dataset module should have at least an `__init__` and a `yield_points` method. `yield_points` should yield a numpy rec array with at least
## 'x', 'y' and 'z' fields, and optionally 'w' and 'u' fields. Other methods that can be re-written include `parse` which yields
## dlim dataset module(s).
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
## my_processed_datalist = init_data(my_data_list, src_region = my_region, x_inc = 1s, y_inc = 1s)
## my_stack = my_processed_datalist._stacks() # stack the data to the region/increments
## my_processed_datalist.archive_xyz() # archive the data to the region/increments
## my_mask = my_processed_datalist._mask() # mask the data to the region/increments
##
### TODO:
## mask to stacks for supercede
## fetch results class
## temp files
### Code:

import os
import sys
import re
import copy
import json
import glob
import math
from tqdm import tqdm
import warnings
import traceback

# import threading
# import multiprocessing as mp
# mp.set_start_method('spawn')
# try:
#    import Queue as queue
# except: import queue as queue
        
import numpy as np
# from scipy.spatial import ConvexHull
# import lxml.etree

import pyproj
import utm
import laspy as lp
from osgeo import gdal
from osgeo import ogr
#from osgeo import osr
import h5py as h5
import netCDF4 as nc

import cudem
from cudem import utils
from cudem import regions
from cudem import xyzfun
from cudem import gdalfun
from cudem import factory
from cudem import vdatums
from cudem import fetches
from cudem import grits
from cudem import vrbag
from cudem import srsfun

# cshelph
import pandas as pd
from cudem import cshelph

## Config info and setup
gc = utils.config_check()
gdal.DontUseExceptions()
ogr.DontUseExceptions()
gdal.SetConfigOption('CPL_LOG', 'NUL' if gc['platform'] == 'win32' else '/dev/null') 

## dataset masking functions
def num_strings_to_range(*args):
    """parse args to a number range.

    e.g. if args == ['1934', '1920-1980', '2001'] will return '1920-2001'
    """
    
    dates = []
    for arg in args:
        if utils.str_or(arg) is not None:
            dd = re.findall(r'[-+]?\d*\.?\d+', str(arg))
            for d in dd:
                dates.append(abs(float(d)))
            
    if len(dates) > 0:
        if min(dates) != max(dates):
            return('-'.join([str(min(dates)), str(max(dates))]))
        else:
            return(str(min(dates)))
    else:
        return(None)

## spatial metadata
def scan_mask_bands(src_ds, skip_band = 'Full Data Mask', mask_level = 0, verbose = True):
    mask_level = utils.int_or(mask_level, 0)
    src_infos = gdalfun.gdal_infos(src_ds)
    band_infos = {}    
    with tqdm(
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
                band_infos[this_band_name] = {'bands': [b], 'metadata': this_band_md}
            else:
                band_infos[this_band_name]['bands'].append(b)
                #band_infos[this_band_name]['metadata'] = this_band_md
                
            band_md = band_infos[this_band_name]['metadata']
            for k in this_band_md.keys():
                if k not in band_md.keys() or band_md[k] is None or band_md[k] == 'None':
                    if this_band_md[k] != 'None':
                        band_md[k] = this_band_md[k]
                            
                if k == 'name':
                    band_md[k] = this_band_name                
                elif k == 'date' or k == 'weight' or k == 'uncertainty':
                    band_md[k] = num_strings_to_range(band_md[k], this_band_md[k])
                elif band_md[k] not in this_band_md[k].split('/') and band_md[k] != this_band_md[k]:
                    if str(this_band_md[k]) != 'None':
                        band_md[k] = ','.join([band_md[k], this_band_md[k]])
                    
            band_infos[this_band_name]['metadata'] = band_md
            
    return(band_infos)

## breaks on many bands
def ogr_mask_footprints(
        src_ds, ogr_format = 'GPKG', dst_srs = 'epsg:4326', mask_level = 0, verbose = True
):
    
    src_infos = scan_mask_bands(src_ds, mask_level=mask_level, verbose=verbose)
    #utils.echo_msg(src_infos)
    ## initialize output vector
    dst_ogr_fn = '{}_sm.{}'.format(utils.fn_basename2(src_ds.GetDescription()), gdalfun.ogr_fext(ogr_format))
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
    with tqdm(
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
            
    with tqdm(
            desc='setting mask metadata',
            total=len(src_infos.keys()),
            leave=verbose
    ) as pbar:
        for i, key in enumerate(src_infos.keys()):
            ## update db using ogrinfo
            key_val = ', '.join(["'{}' = '{}'".format(x, src_infos[key]['metadata'][x]) for x in src_infos[key]['metadata'].keys()])            
            sql = "UPDATE {name} SET {key_val} WHERE rowid = {f}".format(
                name='footprint',
                key_val=key_val,
                value=src_infos[key]['metadata'][md],
                f=i+1
            )
            utils.run_cmd('ogrinfo {} -dialect SQLite -sql "{}"'.format(dst_ogr_fn, sql), verbose=False)
            pbar.update()
            
    return(glob.glob('{}.*'.format(utils.fn_basename2(dst_ogr_fn))), ogr_format)
            
def merge_mask_bands_by_name(src_ds, verbose = True, skip_band = 'Full Data Mask', mask_level = 0):
    """merge data mask raster bands based on the top-level datalist name"""

    band_infos = scan_mask_bands(src_ds, mask_level=mask_level, verbose=verbose)
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

    with tqdm(
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

            with tqdm(
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
        
    # with tqdm(
    #         desc='merging mask bands.',
    #         total=src_ds.RasterCount,
    #         leave=verbose
    # ) as pbar:         
    #     for b in range(1, src_ds.RasterCount+1):
    #         pbar.update()
    #         this_band = src_ds.GetRasterBand(b)
    #         this_band_md = this_band.GetMetadata()
    #         this_band_name = this_band.GetDescription()
    #         if this_band_name == skip_band:
    #             continue

    #         # if mask_level > 0:
    #         #     if mask_level > len(this_band_name.split('/')):
    #         #         mask_level = len(this_band_name.split('/')) - 2

    #         #     this_band_name = '/'.join(this_band_name.split('/')[:-mask_level])

    #         if len(this_band_name.split('/')) > 1:
    #             this_band_name = this_band_name.split('/')[1]

    #         m_bands = {m_ds.GetRasterBand(i).GetDescription(): i for i in range(1, m_ds.RasterCount + 1)}
    #         if not this_band_name in m_bands.keys():
    #             m_ds.AddBand()
    #             m_band = m_ds.GetRasterBand(m_ds.RasterCount)
    #             m_band.SetNoDataValue(0)
    #             m_band.SetDescription(this_band_name)                
    #         else:
    #             m_band = m_ds.GetRasterBand(m_bands[this_band_name])                
                
    #         band_md = m_band.GetMetadata()
    #         for k in this_band_md.keys():
    #             if k not in band_md.keys() or band_md[k] is None:
    #                 band_md[k] = this_band_md[k]

    #             if k == 'name':
    #                 band_md[k] = this_band_name                
    #             elif k == 'date' or k == 'weight' or k == 'uncertainty':
    #                 band_md[k] = num_strings_to_range(band_md[k], this_band_md[k])
    #             elif band_md[k] != this_band_md[k]:
    #                 band_md[k] = ','.join([band_md[k], this_band_md[k]])

    #         try:
    #             m_band.SetMetadata(band_md)
    #         except:
    #             try:
    #                 for key in band_md.keys():
    #                     if band_md[key] is None:
    #                         del band_md[key]

    #                 m_band.SetMetadata(band_md)
    #             except Exception as e:
    #                 utils.echo_error_msg(
    #                     'could not set band metadata: {}; {}'.format(
    #                         band_md, e
    #                     )
    #                 )

    #         this_band_array = this_band.ReadAsArray()
    #         m_band_array = m_band.ReadAsArray()
    #         out_array = np.add(m_band_array, this_band_array)
    #         out_array[out_array > 1] = 1
    #         m_band.WriteArray(out_array)
    #         m_ds.FlushCache()
        
    # return(m_ds)
        
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
    src_ds (GDALDataset): the source multi-band raster as a gdal dataset object
    dst_srs (str): the output srs
    ogr_format (str): the output OGR format
    verbose (bool): be verbose

    -----------
    Returns:
    tuple: (dst_layer, ogr_format)
    """

    mask_level = utils.int_or(mask_level, 0)
    src_ds = merge_mask_bands_by_name(src_ds, mask_level=mask_level, verbose=verbose)
    dst_layer = '{}_sm'.format(utils.fn_basename2(src_ds.GetDescription()) if output is None else output)
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
    layer.StartTransaction()
    defn = None

    with tqdm(
            desc='polygonizing mask bands.',
            total=src_ds.RasterCount,
            leave=verbose
    ) as pbar:
        for b in range(1, src_ds.RasterCount+1):
            pbar.update()
            this_band = src_ds.GetRasterBand(b)
            this_band_md = this_band.GetMetadata()
            this_band_name = this_band.GetDescription()
            #if len(this_band_name.split('/')) > 1:
            #    this_band_name = this_band_name.split('/')[-2]

            b_infos = gdalfun.gdal_infos(src_ds, scan=True, band=b)
            field_names = [field.name for field in layer.schema]
            this_band_md = {k.title():v for k,v in this_band_md.items()}            
            for k in this_band_md.keys():
                if k[:9] not in field_names:
                    layer.CreateField(ogr.FieldDefn(k[:9], ogr.OFTString))

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
                        tmp_layer.CreateField(ogr.FieldDefn(k[:9], ogr.OFTString))

                    if 'Title' not in this_band_md.keys():
                        tmp_layer.CreateField(ogr.FieldDefn('Title', ogr.OFTString))

                    # if verbose:
                    #     utils.echo_msg('polygonizing {} mask...'.format(this_band_name))

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
                        with tqdm(
                                desc='creating feature {}...'.format(this_band_name),
                                total=len(this_band_md.keys()),
                                leave=verbose
                        ) as feat_pbar:
                            for k in this_band_md.keys():
                                feat_pbar.update()
                                out_feat.SetField(k[:9], this_band_md[k])

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
    dst_layer = '{}_sm'.format(utils.fn_basename2(src_ds.GetDescription()) if output is None else output)
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
    
    with tqdm(
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

                with tqdm(
                        desc='polygonizing {} mask bands from {}.'.format(len(band_infos[m]['bands']), m),
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
                    with tqdm(
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

## Datalist convenience functions
## data_list is a list of dlim supported datasets
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

## initialize a list of datasets into a dataset object
def init_data(data_list,
              region=None,
              src_srs=None,
              dst_srs=None,
              src_geoid=None,
              dst_geoid='g2018',
              xy_inc=(None, None),
              sample_alg='auto',
              want_weight=False,
              want_uncertainty=False,
              want_verbose=True,
              want_mask=False,
              want_sm=False,
              invert_region=False,
              cache_dir=None,
              dump_precision=4,
              pnt_fltrs=None,
              stack_fltrs=None,
              stack_node=True,
              stack_mode='mean',
              upper_limit=None,
              lower_limit=None,
              mask=None):
    """initialize a datalist object from a list of supported dataset entries"""

    try:
        xdls = [DatasetFactory(
            mod=" ".join(['-' if x == "" else x for x in dl.split(",")]),
            data_format = None,
            weight=None if not want_weight else 1,
            uncertainty=None if not want_uncertainty else 0,
            src_srs=src_srs,
            dst_srs=dst_srs,
            src_geoid=src_geoid,
            dst_geoid=dst_geoid,
            x_inc=xy_inc[0],
            y_inc=xy_inc[1],
            sample_alg=sample_alg,
            parent=None,
            src_region=region,
            invert_region=invert_region,
            cache_dir=cache_dir,
            want_mask=want_mask,
            want_sm=want_sm,
            verbose=want_verbose,
            dump_precision=dump_precision,
            pnt_fltrs=pnt_fltrs,
            stack_fltrs=stack_fltrs,
            stack_node=stack_node,
            stack_mode=stack_mode,
            upper_limit=None,
            lower_limit=None,
            mask=mask
        )._acquire_module() for dl in data_list]

        if len(xdls) > 1:
            this_datalist = Scratch(
                fn=xdls,
                data_format=-3,
                weight=None if not want_weight else 1,
                uncertainty=None if not want_uncertainty else 0,
                src_srs=src_srs,
                dst_srs=dst_srs,
                src_geoid=src_geoid,
                dst_geoid=dst_geoid,
                x_inc=xy_inc[0],
                y_inc=xy_inc[1],
                sample_alg=sample_alg,
                parent=None,
                src_region=region,
                invert_region=invert_region,
                cache_dir=cache_dir,
                want_mask=want_mask,
                want_sm=want_sm,
                verbose=want_verbose,
                dump_precision=dump_precision,
                pnt_fltrs=pnt_fltrs,
                stack_fltrs=stack_fltrs,
                stack_node=stack_node,
                stack_mode=stack_mode,
                upper_limit=None,
                lower_limit=None,
                mask=mask
            )
        elif len(xdls) > 0:
            this_datalist = xdls[0]
        else:
            this_datalist = None

        return(this_datalist)
    
    except Exception as e:
        utils.echo_error_msg('could not initialize data, {}: {}'.format(data_list, e))
        return(None)

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
    
## TODO: move pointfilter/binz/etc to own module
class PointFilter:
    def __init__(self, points = None, params = {}, verbose = True, **kwargs):
        self.points = points
        self.params = params
        self.verbose = verbose
        self.kwargs = kwargs

    def __call__(self):
        if self.verbose:
            utils.echo_msg('filtering points using {}'.format(self))
            
        return(self.run())

    def run(self):
        raise(NotImplementedError)
    
    def convert_wgs_to_utm(self, lat, lon):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            easting, northing, num, letter = utm.from_latlon(lat, lon)
            if letter >= 'N':
                epsg = 'epsg:326' + str(num)
            elif letter < 'N':
                epsg = 'epsg:327' + str(num)
            else:
                utils.echo_error_msg('Could not find UTM zone')

            return(epsg)

## Adapted from CShelph
class BinZ(PointFilter):
    def __init__(self, y_res = 1, z_res = .5, percentile = 50, z_min = None, z_max = None, **kwargs):
        super().__init__(**kwargs)
        self.y_res = utils.float_or(y_res, 1)
        self.z_res = utils.float_or(z_res, .5)
        self.percentile = utils.float_or(percentile, 50)
        self.z_min = utils.float_or(z_min)
        self.z_max = utils.float_or(z_max)
        
    def bin_points(self, points):
        '''Bin data along vertical and horizontal scales for later segmentation'''

        ## Calculate number of bins required both vertically and
        ## horizontally with resolution size
        y_bin_number = round(abs(points['y'].min() - points['y'].max())/self.y_res)
        if y_bin_number == 0:
            y_bin_number = 1
            
        #x_bin_number = round(abs(points['x'].min() - points['x'].max())/self.y_res)
        z_bin_number = round(abs(points['z'].min() - points['z'].max())/self.z_res)
        if z_bin_number == 0:
            z_bin_number = 1

        if (y_bin_number > 0 and z_bin_number > 0):    
            points1 = points
            y_bins = pd.cut(points['y'], y_bin_number, labels = np.array(range(y_bin_number)))
            #x_bins = pd.cut(points['x'], x_bin_number, labels = np.array(range(x_bin_number)))
            points1['y_bins'] = y_bins
            #points1['x_bins'] = x_bins
            #utils.echo_msg('{} {} {}'.format(points['z'].min(), points['z'].max(), z_bin_number))
            z_bins = pd.cut(
                points['z'], z_bin_number, labels = np.round(
                    np.linspace(points['z'].min(), points['z'].max(), num=z_bin_number),
                    decimals = 1
                )
            )
            points1['z_bins'] = z_bins
            points1 = points1.reset_index(drop=True)

            return(points1)

        return(None)

    def get_bin_height(self, binned_data):
        '''Calculate mean sea height for easier calculation of depth and cleaner figures'''

        # Create sea height list
        bin_height = []
        bin_lat = []
        bin_lon = []

        # Group data by latitude
        binned_data_sea = binned_data
        #grouped_data = binned_data_sea.groupby(['y_bins', 'x_bins'], group_keys=True)
        grouped_data = binned_data_sea.groupby(['y_bins'], group_keys=True)
        data_groups = dict(list(grouped_data))

        # Create a percentile threshold of photon counts in each grid, grouped by both x and y axes.
        # count_threshold = np.percentile(
        #     binned_data.groupby(['y_bins', 'x_bins', 'z_bins']).size().reset_index().groupby('y_bins')[[0]].max(),
        #     percentile
        # )
        count_threshold = np.percentile(
            binned_data.groupby(['y_bins', 'z_bins']).size().reset_index().groupby('y_bins')[[0]].max(),
            self.percentile
        )
        # Loop through groups and return average sea height
        for k,v in data_groups.items():
            # Create new dataframe based on occurance of photons per height bin
            new_df = pd.DataFrame(v.groupby('z_bins').count())

            # Return the bin with the highest count
            largest_h_bin = new_df['z'].argmax()

            # Select the index of the bin with the highest count
            largest_h = new_df.index[largest_h_bin]

            # Set threshold of photon counts per bin
            if new_df.iloc[largest_h_bin]['y'] >= count_threshold:
                [bin_lat.append(x) for x in v.loc[v['z_bins']==largest_h, 'y']]
                [bin_lon.append(x) for x in v.loc[v['z_bins']==largest_h, 'x']]
                [bin_height.append(x) for x in v.loc[v['z_bins']==largest_h, 'z']]                
                del new_df
            else:
                del new_df

        # Filter out sea height bin values outside 2 SD of mean.
        if np.all(np.isnan(bin_height)):
            return(None)

        mean = np.nanmean(bin_height, axis=0)
        sd = np.nanstd(bin_height, axis=0)
        bin_height_1 = np.where(
            (bin_height > (mean + 2*sd)) | (bin_height < (mean - 2*sd)), np.nan, bin_height
        ).tolist()
                
        return(bin_lat, bin_lon, bin_height_1)
    
    def run(self):
        try:
            epsg_code = self.convert_wgs_to_utm(self.points['y'][0], self.points['x'][0])
            epsg_num = int(epsg_code.split(':')[-1])
            utm_proj = pyproj.Proj(epsg_code)
            x_utm, y_utm = utm_proj(self.points['x'], self.points['y'])
        except:
            utils.echo_warning_msg('could not transform to utm')
            x_utm, y_utm = self.points['x'], self.points['y']

        points_1 = pd.DataFrame(
            {'y': y_utm,
             'x': x_utm,
             'z': self.points['z']},
            columns=['y', 'x', 'z']
        )

        if utils.float_or(self.z_min) is not None:
            points_1 = points_1[(points_1['z'] > self.z_min)]

        if utils.float_or(self.z_max) is not None:
            points_1 = points_1[(points_1['z'] < self.z_max)]

        #print(points_1)
        #print(len(points_1))
        if len(points_1) > 0:
            binned_points = self.bin_points(points_1)
            points_1 = None
            
            if binned_points is not None:
                ys, xs, zs = self.get_bin_height(binned_points)
                binned_points = None
                
                bin_ds = np.column_stack((xs, ys, zs))
                bin_ds = np.rec.fromrecords(bin_ds, names='x, y, z')
                
                xs = ys = zs = None
                bin_ds = bin_ds[~np.isnan(bin_ds['z'])]
                med_surface_h = np.nanmedian(bin_ds['z'])
                #bin_ds = bin_ds[bin_ds['z'] < med_surface_h + (z_res * 2)]
                #bin_ds = bin_ds[bin_ds['z'] > med_surface_h - (z_res * 2)]

                transformer = pyproj.Transformer.from_crs(
                    "EPSG:"+str(epsg_num), "EPSG:4326", always_xy=True
                )
                lon_wgs84, lat_wgs84 = transformer.transform(bin_ds['x'], bin_ds['y'])
                bin_points = np.column_stack((lon_wgs84, lat_wgs84, bin_ds['z']))
                lon_wgs84 = lat_wgs84 = bin_ds = None
                bin_points = np.rec.fromrecords(bin_points, names='x,y,z')
                #utils.echo_msg('bin_points: {}'.format(bin_points))
                return(bin_points)
            # return(bin_ds)
            
        return(points_1)

## Adapted from https://medium.com/@anthony.klemm/dealing-with-outliers-in-crowdsourced-bathymetry-data-using-predictive-mean-matching-pmm-bc490c158c7e
class PMMOutliers(PointFilter):
    def __init__(
            self, percentile = 98, max_percentile = 99.9,
            multipass = 1, want_iqr = False, return_outliers = False,
            **kwargs
    ):
        super().__init__(**kwargs)
        self.percentile = utils.float_or(percentile, 98)
        self.max_percentile = utils.float_or(max_percentile, 99)
        self.multipass = utils.int_or(multipass, 1)
        self.want_iqr = want_iqr
        self.return_outliers = return_outliers

    def pmm(self, points, percentile = 98, want_iqr=False, return_outliers=False):
        p = np.array([self.points['y'], self.points['x'], self.points['z']]).T
        try:
            from sklearn.experimental import enable_iterative_imputer
            from sklearn.impute import IterativeImputer
            from sklearn.preprocessing import StandardScaler

            scaler = StandardScaler(with_mean=False)
            #data_scaled = scaler.fit_transform(data[[‘lat’, ‘lon’, ‘depth’]])
            #print(np.array([self.points['y'], self.points['x'], self.points['z']]))
            #p = np.array([self.points['y'], self.points['x'], self.points['z']]).T
            #data_scaled = scaler.fit_transform(np.array([self.points['y'], self.points['x'], self.points['z']]))
            data_scaled = scaler.fit_transform(p)
            imputer = IterativeImputer(random_state=36, sample_posterior=True)
            data_imputed = imputer.fit_transform(data_scaled)

        except:
            utils.echo_warning_msg('you must install sklearn for imputing')
            data_imputed = p

        from scipy.ndimage import uniform_filter1d
        smoothed_depth = uniform_filter1d(data_imputed[:, 2], size=50)
        residuals = np.abs(data_imputed[:, 2] - smoothed_depth)
        
        data_imputed = smoothed_depth = imputer = None
        #smoothed_depth = uniform_filter1d(points['z'], size=50)
        #residuals = np.abs(points['z'] - smoothed_depth)
        if want_iqr:
            perc_max = np.nanpercentile(residuals, percentile)            
            iqr_p = (perc_max - (100-percentile)) * 1.5
            outlier_threshold = perc_max + iqr_p
        else:
            outlier_threshold = np.percentile(residuals, percentile)

        outliers = residuals >= outlier_threshold

        if self.verbose:
            utils.echo_msg_bold('removed {} outliers @ {}'.format(np.count_nonzero(outliers), percentile))

        if return_outliers:
            return(points[outliers])
        else:
            return(points[~outliers])
            #return(data_imputed[~outliers])
        
    def run(self):

        #if self.want_iqr:
        #    percs_it = np.linspace(self.max_percentile, self.percentile, self.multipass)
        #else:
        percs_it = np.linspace(self.percentile, self.max_percentile, self.multipass)            
        for mpass in range(0, self.multipass):
            self.points = self.pmm(self.points, percentile=percs_it[mpass], want_iqr=self.want_iqr, return_outliers=self.return_outliers)
            
        return(self.points)

class RQOutliers(PointFilter):
    def __init__(
            self, percentile = 98, max_percentile = 99.9,
            multipass = 2, want_iqr = False, return_outliers = False,
            src_raster = None, **kwargs
    ):
        super().__init__(**kwargs)
        self.percentile = utils.float_or(percentile, 98)
        self.max_percentile = utils.float_or(max_percentile, 99)
        self.multipass = utils.int_or(multipass, 1)
        self.want_iqr = want_iqr
        self.return_outliers = return_outliers
        self.src_raster = src_raster

    def _fetch_gmrt(self, gmrt_region = None, dst_srs = None):
        """GMRT - Global low-res.
        """
        
        this_gmrt = fetches.GMRT(
            src_region=gmrt_region, verbose=self.verbose, layer='topo', outdir='./'#outdir=self.cache_dir
        )
        this_gmrt.run()

        fr = fetches.fetch_results(this_gmrt)
        fr.daemon = True
        fr.start()
        fr.join()

        # dst_srs = osr.SpatialReference()
        # dst_srs.SetFromUserInput(self.dst_srs)
        
        gmrt_tif = os.path.join(this_gmrt._outdir, this_gmrt.results[0]['dst_fn'])
        # gmrt_ds = gdalfun.gdal_mem_ds(self.ds_config, name='gmrt', co=self.co)
        # gdal.Warp(gmrt_ds, gmrt_tif, dstSRS=dst_srs, resampleAlg=self.sample)
        # return(gmrt_ds)
        return(gmrt_tif)
    
    def rq(self, points, raster, percentile = 98, want_iqr=False, return_outliers=False):
        p = np.array([points['y'], points['x'], points['z']]).T
        raster_z = gdalfun.gdal_query(points, raster, 'g').flatten()            
        residuals = np.abs(p[:,2] - raster_z)
        #smoothed_depth = uniform_filter1d(points['z'], size=50)
        #residuals = np.abs(points['z'] - smoothed_depth)
        if want_iqr:
            perc_max = np.nanpercentile(residuals, percentile)            
            iqr_p = (perc_max - (100-percentile)) * 1.5
            outlier_threshold = perc_max + iqr_p
        else:
            outlier_threshold = np.percentile(residuals, percentile)

        outliers = residuals >= outlier_threshold

        if self.verbose:
            utils.echo_msg_bold('removed {} outliers @ {}'.format(np.count_nonzero(outliers), percentile))

        if return_outliers:
            return(points[outliers])
        else:
            return(points[~outliers])
            #return(data_imputed[~outliers])
        
    def run(self):

        if self.src_raster is None:
            region = regions.Region().from_list([np.min(self.points['x']), np.max(self.points['x']), np.min(self.points['y']), np.max(self.points['y'])])            
            self.src_raster = self._fetch_gmrt()
        #if self.want_iqr:
        #    percs_it = np.linspace(self.max_percentile, self.percentile, self.multipass)
        #else:
        percs_it = np.linspace(self.percentile, self.max_percentile, self.multipass)
            
        for mpass in range(0, self.multipass):
            self.points = self.rq(
                self.points,
                self.src_raster,
                percentile=percs_it[mpass],
                want_iqr=self.want_iqr,
                return_outliers=self.return_outliers
            )
            
        return(self.points)
    
class PointFilterFactory(factory.CUDEMFactory):
    _modules = {
        'bin_z': {'name': 'bin_z', 'call': BinZ},
        'pmm': {'name': 'pmm', 'call': PMMOutliers},
        'rq': {'name': 'rq', 'call': RQOutliers},
    }
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
    
class INF:
    """INF Files contain information about datasets"""
    
    def __init__(self,
                 name = None,
                 file_hash = None,
                 numpts = 0,
                 minmax = [],
                 wkt = None,
                 fmt = None,
                 src_srs = None):
        self.name = name
        self.file_hash = file_hash
        self.hash = file_hash
        self.numpts = numpts
        self.minmax = minmax
        self.wkt = wkt
        self.fmt = fmt
        self.format = fmt
        self.src_srs = src_srs

    def __str__(self):
        return('<Dataset Info: {}>'.format(self.__dict__))
    
    def __repr__(self):
        return('<Dataset Info: {}>'.format(self.__dict__))

    def generate_hash(self, fn = None, sha1 = False):
        """generate a hash of the xyz-dataset source file"""

        import hashlib
        fn = self.name if fn is None else fn
        BUF_SIZE = 65536
        if sha1:
            this_hash = hashlib.sha1()
        else:
            this_hash = hashlib.md5()
            
        try:
            with open(fn, 'rb') as f:
                while True:
                    data = f.read(BUF_SIZE)
                    if not data:
                        break

                    this_hash.update(data)

            self.file_hash = this_hash.hexdigest()
        
        except:
            self.file_hash = '0'

        return(self.file_hash)

    def generate_mini_grid(self, x_size = 10, y_size = 10):
        """generate a 'mini-grid' of the data in about a 10x10 grid.
        this will help us determine the location of data before processing.
        """
        
        raise(NotImplementedError)
    
    def generate(self):
        if self.name is None:
            return(self)

        this_ds = DatasetFactory(mod=self.name)
        this_ds_mod = this_ds._acquire_module()        
        self.fmt = this_ds.mod_name            
        
    def load_inf_file(self, inf_fn = None):
        valid_data = False
        data = {}
        if inf_fn is None:
            return(self)
        
        if os.path.exists(inf_fn):
            try:
                with open(inf_fn, 'r') as inf_ob:
                    data = json.load(inf_ob)
            except ValueError:
                try:
                    data = MBSParser(fn=inf_fn).inf_parse().infos.__dict__
                    #self.check_hash = False
                    mb_inf = True
                except Exception as e:
                    raise ValueError(
                        'CUDEMFactory: Unable to read data from {} as mb-system inf, {}'.format(inf_fn, e)
                    )
            except:
                raise ValueError(
                    'CUDEMFactory: Unable to read data from {} as json'.format(inf_fn)
                )

        for ky, val in data.items():
            if ky in self.__dict__:
                self.__setattr__(ky, val)
                valid_data = True

        return(self)
            
    def write_inf_file(self, inf_fn = None):
        if inf_fn is None:
            if self.name is not None:
                inf_fn = '{}.inf'.format(self.name)
        try:
            with open(inf_fn, 'w') as outfile:
               json.dump(self.__dict__, outfile)
        except:
            pass
        
class ElevationDataset:
    """representing an Elevation Dataset
    
    This is the super class for all datalist (dlim) datasets .    
    Each dataset sub-class should define a dataset-specific    
    data parser <self.parse> and a <self.generate_inf> function to generate
    inf files. 

    Specifically, each sub-dataset should minimally define the following functions:

    sub_ds.__init__
    sub_ds.yield_points

    -----------
    Where:

    generate_inf generates a dlim compatible inf file,
    parse yields the dlim dataset module
    yield_points yields the data as a numpy rec-array with 'x', 'y', 'z' and optionally 'w' and 'u'

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

    ## todo: add transformation grid option (stacks += transformation_grid), geoids
    def __init__(self,
                 fn = None,
                 data_format = None,
                 weight = 1,
                 uncertainty = 0,
                 src_srs = None,
                 mask = None,
                 dst_srs = 'epsg:4326',
                 src_geoid = None,
                 dst_geoid = 'g2018',
                 x_inc = None,
                 y_inc = None,
                 want_mask = False,
                 want_sm = False,
                 sample_alg = 'auto',
                 parent = None,
                 src_region = None,
                 invert_region = False,
                 pnt_fltrs=[],
                 stack_fltrs=[],
                 stack_node = True,
                 stack_mode = 'mean',
                 cache_dir = None,
                 verbose = False,
                 remote = False,
                 dump_precision = 6,
                 upper_limit = None,
                 lower_limit = None,
                 params = {},
                 metadata = {'name':None, 'title':None, 'source':None, 'date':None,
                             'data_type':None, 'resolution':None, 'hdatum':None,
                             'vdatum':None, 'url':None}, **kwargs):
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
        self.invert_region = invert_region # invert the region
        self.cache_dir = cache_dir # cache_directory
        self.verbose = verbose # be verbose
        self.remote = remote # dataset is remote
        self.dump_precision = dump_precision # the precision of the dumped xyz data
        self.stack_node = stack_node # yield avg x/y data instead of center
        self.stack_mode = stack_mode # 'mean', 'min', 'max', 'supercede'
        self._init_stack_mode()
        self.mask_keys = ['mask', 'invert_mask', 'ogr_or_gdal'] # options for input data mask
        if self.mask is not None:
            if isinstance(self.mask, str):
                self.mask = {'mask': self.mask}
                
            for kpam, kval in kwargs.items():
                if kpam not in self.__dict__:
                    self.mask[kpam] = kval

            for kpam, kval in self.mask.items():
                if kpam in kwargs:
                    del kwargs[kpam]

            self.mask['ogr_or_gdal'] = gdalfun.ogr_or_gdal(self.mask['mask'])
            for key in self.mask_keys:
                if key not in self.mask.keys():
                    self.mask[key] = None
                    
        self.upper_limit = utils.float_or(upper_limit)
        self.lower_limit = utils.float_or(lower_limit)

        self.infos = INF(name=self.fn, file_hash='0', numpts=0, fmt=self.data_format) # infos blob

        self.params = params # the factory parameters
        if not self.params:
            self.params['kwargs'] = self.__dict__.copy()
            self.params['mod'] = self.fn
            self.params['mod_name'] = self.data_format
            self.params['mod_args'] = {}
        
    def __str__(self):
        return('<Dataset: {} - {}>'.format(self.metadata['name'], self.fn))
    
    def __repr__(self):
        return('<Dataset: {} - {}>'.format(self.metadata['name'], self.fn))

    def __call__(self):
        self.initialize()

    def _init_stack_mode(self):
        opts = self.stack_mode.split(':')
        self.stack_mode_name = opts[0]
        self.stack_mode_args = {}
        if self.stack_mode_name not in self.stack_modes:
            self.stack_mode_name = 'mean'
            
        elif len(opts) > 1:
            self.stack_mode_args = factory.args2dict(list(opts[1:]), self.stack_mode_args)

        if 'mask_level' not in self.stack_mode_args.keys():
            self.stack_mode_args['mask_level'] = 0
        
    def _init_mask(self):
        if self.mask is not None:
            if isinstance(self.mask, str):
                self.mask = {'mask': self.mask}
                
            for kpam, kval in kwargs.items():
                if kpam not in self.__dict__:
                    self.mask[kpam] = kval

            for kpam, kval in self.mask.items():
                if kpam in kwargs:
                    del kwargs[kpam]

            self.mask['ogr_or_gdal'] = gdalfun.ogr_or_gdal(self.mask['mask'])
            for key in self.mask_keys:
                if key not in self.mask.keys():
                    self.mask[key] = None
        
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
        
    def initialize(self):
        self._fn = None # temp filename holder
        self.data_region = None # self.region and inf.region reduced
        self.archive_datalist = None # the datalist of the archived data
        self.data_entries = [] # 
        self.data_lists = {} #
        self.cache_dir = utils.cudem_cache() if self.cache_dir is None else self.cache_dir # cache directory
        if self.sample_alg not in self.gdal_sample_methods: # gdal_warp resmaple algorithm 
            utils.echo_warning_msg(
                '{} is not a valid gdal warp resample algorithm, falling back to `auto`'.format(
                    self.sample_alg
                )
            )
            self.sample_alg = 'auto'

        if utils.fn_url_p(self.fn):
            self.remote = True

        ## Dataset Mask, if the input mask is a vector, rasterize it if possible.
        # if self.mask is not None:
        #     if self.mask['ogr_or_gdal'] == 1: # mask is ogr, rasterize it
        #         if self.region is not None and self.x_inc is not None and self.y_inc is not None:
        #             self.mask['mask'] = gdalfun.ogr2gdal_mask(
        #                 self.mask['mask'], region=self.region, x_inc=self.x_inc, y_inc=self.y_inc, dst_srs=self.dst_srs,
        #                 invert=True, verbose=False, temp_dir=self.cache_dir
        #             )
        #             self.mask['ogr_or_gdal'] = 0

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
            self.infos = self.inf(check_hash=True if self.data_format == -1 else False)
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                try:
                    self.set_transform()
                except Exception as e:
                    utils.echo_error_msg('could not set transformation on {}, {}'.format(self.fn, e))

            self.set_yield()

        ## initialize filters
        if isinstance(self.stack_fltrs, str):
            self.stack_fltrs = [':'.join(self.stack_fltrs.split('/'))]
            #self.stack_fltrs = [self.stack_fltrs]

        if isinstance(self.pnt_fltrs, str):
            self.pnt_fltrs = [':'.join(self.pnt_fltrs.split('//'))]

        #self._init_stack_mode()
            
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
        #if self.fn.startswith('http') or self.fn.startswith('/vsicurl/') or self.fn.startswith('BAG'):
        #    check_path = False
            
        #if check_path:
        if self.fn is not None:
            if not isinstance(self.fn, list) and \
               not isinstance(self.fn, np.ndarray) and \
               not isinstance(self.fn, np.core.records.recarray):
                if self.fn not in fmts:
                    if self.fn.startswith('http') or self.fn.startswith('/vsicurl/') or self.fn.startswith('BAG'):
                        if not utils.fn_url_p(self.fn):
                            if self.data_format > -10:
                                if not os.path.exists(self.fn):
                                    return (False)

                                if os.stat(self.fn).st_size == 0:
                                    return(False)
                        
        return(True)
                
    def echo_(self, sep=' ', **kwargs):
        """print self as a datalist entry string"""

        out = []
        for key in self.metadata.keys():
            if key != 'name':
                out.append(str(self.metadata[key]))

        return(sep.join(['"{}"'.format(str(x)) for x in out]))
    
    def echo(self, **kwargs):
        """print self.data_entries as a datalist entries."""

        for entry in self.parse():
            entry_path = os.path.abspath(entry.fn) if not self.remote else entry.fn
            l = [entry_path, entry.data_format]
            if entry.weight is not None:
                l.append(entry.weight)
                
            print('{}'.format(" ".join([str(x) for x in l])))

    def format_entry(self, sep = ' '):
        """format the dataset information as a `sep` separated string."""
        
        dl_entry = sep.join(
            [str(x) for x in [
                self.fn,
                '{}:{}'.format(self.data_format, factory.dict2args(self.params['mod_args'])),
                self.weight,
                self.uncertainty
            ]]
        )
        metadata = self.echo_()
        return(sep.join([dl_entry, metadata]))

            
    def format_metadata(self, **kwargs):
        """format metadata from self, for use as a datalist entry."""

        return(self.echo_())

    def set_yield(self):
        """set the yield strategy, either default (all points) or mask or stacks

        This sets both `self.array_yeild` and `self.xyz_yield`.
        
        if region, x_inc and y_inc are set, data will yield through _stacks,
        otherwise, data will yield directly from `self.yield_xyz` and `self.yield_array`
        """

        self.array_yield = self.mask_and_yield_array()
        self.xyz_yield = self.mask_and_yield_xyz()
        #if self.want_archive: # archive only works when yielding xyz data.
        #    self.xyz_yield = self.archive_xyz()
        #if (self.region is None and self.transform['trans_region'] is None) and self.x_inc is not None:
        if self.region is None and self.x_inc is not None:
            utils.echo_msg(self.transform['trans_region'])
            utils.echo_warning_msg('must enter a region to output in increment blocks...')

        if self.region is not None and self.x_inc is not None:
            self.x_inc = utils.str2inc(self.x_inc)
            if self.y_inc is None:
                self.y_inc = self.x_inc
            else:
                self.y_inc = utils.str2inc(self.y_inc)

                #out_name = utils.make_temp_fn('dlim_stacks', temp_dir=self.cache_dir)
            out_name = utils.make_temp_fn(utils.append_fn('dlim_stacks', self.region, self.x_inc), temp_dir=self.cache_dir)
            self.xyz_yield = self.stacks_yield_xyz(out_name=out_name)#, fmt='GTiff')#, mode=self.stack_mode)

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
                utils.echo_error_msg(src_srs)
                utils.echo_error_msg(src_srs)

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
                self.transform['src_vert_epsg'], self.transform['dst_vert_epsg'], vd_region.format('fn')
            )
        )
        
        ## if the transformation grid already exists, skip making a new one,
        ## otherwise, make the new one here with `vdatums.VerticalTransform()`
        if not os.path.exists(self.transform['trans_fn']):
            with tqdm(
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
                self.transform['trans_fn'], self.transform['trans_fn_unc'] = vdatums.VerticalTransform(
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

        #utils.echo_msg('using transformation grid {} for {}'.format(self.transform['trans_fn'], self.fn))
        ## set the pyproj.Transformer for both horz+vert and vert only
        ## hack for navd88 datums in feet (6360 is us-feet, 8228 is international-feet
        if utils.str_or(self.transform['src_vert_epsg']) == '6360':# or 'us-ft' in utils.str_or(src_vert, ''):
            #out_src_srs = out_src_srs + ' +vto_meter=0.3048006096012192'
            uc = ' +step +proj=unitconvert +z_in=us-ft +z_out=m'
        elif utils.str_or(self.transform['src_vert_epsg']) == '8228':
            uc = ' +step +proj=unitconvert +z_in=ft +z_out=m'
        else:
            uc = ''
            
        if self.transform['trans_fn'] is not None and os.path.exists(self.transform['trans_fn']):
            self.transform['pipeline'] = '+proj=pipeline{} +step {} +inv +step +proj=vgridshift +grids={} +inv +step {}'.format(
                uc,
                self.transform['src_horz_crs'].to_proj4(),
                os.path.abspath(self.transform['trans_fn']),
                self.transform['dst_horz_crs'].to_proj4()
            )
            self.transform['vert_transformer'] = pyproj.Transformer.from_pipeline(
                '+proj=pipeline{} +step +proj=vgridshift +grids={} +inv'.format(uc, os.path.abspath(self.transform['trans_fn']))
                )
        else:
            utils.echo_error_msg(
                'failed to generate vertical transformation grid between {} and {} for this region!'.format(
                    self.transform['src_vert_epsg'], os.path.abspath(self.transform['dst_vert_epsg'])
                )
            )
            
    def set_transform(self):
        """Set the pyproj horizontal and vertical transformations for the dataset"""

        #in_horizontal_crs, out_horizontal_crs, in_vertical_crs, out_vertical_crs, src_geoid, dst_geoid, want_vertical = self.parse_srs()
        self.parse_srs()
        if self.transform['src_horz_crs'] is not None and self.transform['dst_horz_crs'] is not None:        
            ## horizontal Transformation
            self.transform['horz_pipeline'] = '+proj=pipeline +step {} +inv +step {}'.format(
                self.transform['src_horz_crs'].to_proj4(), self.transform['dst_horz_crs'].to_proj4()
            )
            if self.region is not None:
                self.transform['trans_region'] = self.region.copy()
                self.transform['trans_region'].src_srs = self.transform['dst_horz_crs'].to_proj4()
                self.transform['trans_region'].warp(self.transform['src_horz_crs'].to_proj4())
            else:
                self.transform['trans_region'] = regions.Region().from_list(self.infos.minmax)
                self.transform['trans_region'].src_srs = self.infos.src_srs
                self.transform['trans_region'].warp(self.transform['dst_horz_crs'].to_proj4())

            ## vertical Transformation
            if self.transform['want_vertical']:
                self.set_vertical_transform()
            else:
                self.transform['pipeline'] = self.transform['horz_pipeline']

            try:
                self.transform['transformer'] = pyproj.Transformer.from_pipeline(self.transform['pipeline'])
            except Exception as e:
                utils.echo_warning_msg('could not set transformation in: {}, out: {}, {}'.format(
                    self.transform['src_horz_crs'].name, self.transform['dst_horz_crs'].name, e
                ))
                return

        ## dataset region
        ## mrl: moved out of if block 
        if self.region is not None and self.region.valid_p():
            self.data_region = self.region.copy() \
                if self.transform['trans_region'] is None \
                   else self.transform['trans_region'].copy()
            inf_region = regions.Region().from_list(self.infos.minmax)
            self.data_region = regions.regions_reduce(self.data_region, inf_region)
            self.data_region.src_srs = self.infos.src_srs

            if not self.data_region.valid_p():
                self.data_region = self.region.copy() \
                    if self.transform['trans_region'] is None \
                       else self.transform['trans_region'].copy()
        else:
            self.data_region = regions.Region().from_list(self.infos.minmax)
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
            self.infos.minmax = this_region.export_as_list(include_z=True)
            self.infos.wkt = this_region.export_as_wkt()
            
        self.infos.src_srs = self.src_srs
        if _region is not None:
            self.region = _region.copy()

        if _x_inc is not None:
            self.x_inc = _x_inc
            
        return(self.infos)
    
    def inf(self, check_hash = False, recursive_check = False, write_inf = True, **kwargs):
        """read/write an inf file

        If the inf file is not found, will attempt to generate one.
        The function `generate_inf` should be defined for each specific
        dataset sub-class.
        """

        inf_path = '{}.inf'.format(self.fn)
        generate_inf = False

        ## try to parse the existing inf file as either a native inf json file
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

        if generate_inf:
            self.infos = self.generate_inf()
            
            ## update this
            if self.data_format >= -2 and write_inf:
                self.infos.write_inf_file()

            if recursive_check and self.parent is not None:
                self.parent.inf(check_hash=True)

        if self.infos.src_srs is None:
            self.infos.src_srs = self.src_srs

        else:
            #if self.src_srs is None:
            self.src_srs = self.infos.src_srs

            if self.transform['transformer'] is not None and self.transform['trans_region'] is None:
                self.transform['trans_region'] = regions.Region().from_list(self.infos.minmax)
                self.transform['trans_region'].src_srs = self.infos.src_srs
                self.transform['trans_region'].warp(self.dst_srs)

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

        if self.region is not None:
            try:
                inf_region = regions.Region().from_string(self.infos.wkt)
            except:
                try:
                    inf_region = regions.Region().from_list(self.infos.minmax)
                except:
                    inf_region = self.region.copy()

            if regions.regions_intersect_p(
                    inf_region,
                    self.region if self.transform['trans_region'] is None else self.transform['trans_region']
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
                    self.data_lists[e.parent.metadata['name']] = {'data': [e], 'parent': e.parent}
                    
            else:
                self.data_lists[e.metadata['name']] = {'data': [e], 'parent': e}
                
        return(self)

    ## todo: 'separate mode': multi-band z?
    ## todo: cleanup masking
    def _stacks(
            self,
            out_name = None,
            ndv = -9999,
            fmt = 'GTiff',
            #mask_only = False,
            #mask_level = 0
    ):
        """stack and mask incoming arrays (from `array_yield`) together

        -----------
        Parameters:
        out_name (str): the output stacked raster basename
        ndv (float): the desired no data value
        fmt (str): the output GDAL file format
        mask_level (int): the granularity of the mask, 0 is every file

        --------
        Returns:
        output-file-name{_msk} of a multi-band raster with a band for each dataset.
        output-file-name of a multi-band raster with the following bands:
          x
          y
          z
          weights
          count
          uncertainty
          src uncertainty
        """
        
        def add_mask_band(m_ds, this_entry, mask_level = 0):
            # if mask_level < 0:
            #     m_band = None
            #else:
            m_bands = {m_ds.GetRasterBand(i).GetDescription(): i for i in range(2, m_ds.RasterCount + 1)}
            entry_name = this_entry.metadata['name']
            if mask_level > 0:
                if mask_level > len(entry_name.split('/')):
                    mask_level = len(entry_name.split('/')) - 2

                entry_name = '/'.join(entry_name.split('/')[:-mask_level])
            elif mask_level < 0:
                entry_name = '/'.join(entry_name.split('/')[:2])

            if not entry_name in m_bands.keys():
                m_ds.AddBand()
                m_band = m_ds.GetRasterBand(m_ds.RasterCount)
                m_band.SetNoDataValue(0)
                m_band.SetDescription(entry_name)
            else:
                m_band = m_ds.GetRasterBand(m_bands[entry_name])

            band_md = m_band.GetMetadata()
            for k in this_entry.metadata.keys():
                if k not in band_md.keys() or band_md[k] is None:
                    #utils.echo_msg(this_entry.metadata[k])
                    band_md[k] = this_entry.metadata[k]
                ## mrl commented out for no_wa.datalist
                #else:
                #    band_md[k] = None

            band_md['weight'] = this_entry.weight if this_entry.weight is not None else 1
            band_md['uncertainty'] = this_entry.uncertainty if this_entry.uncertainty is not None else 0
            try:
                m_band.SetMetadata(band_md)
            except:
                try:
                    for key in band_md.keys():
                        if band_md[key] is None:
                            del band_md[key]

                    m_band.SetMetadata(band_md)
                except Exception as e:
                    utils.echo_error_msg(
                        'could not set band metadata: {}; {}'.format(
                            band_md, e
                        )
                    )

            m_band.FlushCache()
                        
            return(m_band)

        def reset_mask_bands(m_ds, srcwin, except_band_name = None, mask = None):
            for b in range(1, m_ds.RasterCount+1):
                this_band = m_ds.GetRasterBand(b)
                this_band_name = this_band.GetDescription()
                
                if this_band_name != except_band_name:
                    this_band_array = this_band.ReadAsArray(srcwin[0], srcwin[1], srcwin[2], srcwin[3])
                    if mask is None:
                        mask = ~np.isnan(this_band_array)

                    masked = this_band_array[(mask)] == 1
                    if np.any(masked):
                        this_band_array[(mask)][(masked)] = 0                        
                        #this_band_array[(mask)] = 0                        
                        this_band.WriteArray(this_band_array, srcwin[0], srcwin[1])
                        m_ds.FlushCache()
                        
                    this_band_array = None
                    
                this_band = None
        
        ## this is a bit convoluted, we're writing the mem mask to disk, then
        ## creating a VRT with the `good_bands` and re-writing that to a final
        ## mask raster...This is all to remove bands that have no data...other methods
        ## are either too slow or take up too much memory...
        ## this is only used if we add a mask band for every parsed entry
        def remove_empty_mask_bands(m_ds):
            if m_ds.RasterCount > 0:
                good_bands = []
                with tqdm(
                        total=m_ds.RasterCount,
                        desc='checking mask bands',
                        leave=self.verbose
                ) as pbar:
                    for band_num in range(1, m_ds.RasterCount+1):
                        pbar.update()
                        band_infos = gdalfun.gdal_infos(m_ds, scan=True, band=band_num)
                        if not np.isnan(band_infos['zr'][0]) and not np.isnan(band_infos['zr'][1]):
                            good_bands.append(band_num)

                pbar = tqdm(desc='writing mask to vrt', total=100, leave=self.verbose)
                pbar_update = lambda a,b,c: pbar.update((a*100)-pbar.n)        
                band_count = len(good_bands)
                msk_ds = gdal.GetDriverByName(fmt).CreateCopy(mask_tmp_fn, m_ds, 0, options=msk_opts, callback=pbar_update)
                pbar.close()
                msk_ds = None
                vrt_options_specific_bands = gdal.BuildVRTOptions(bandList=good_bands)
                vrt_ds = gdal.BuildVRT(mask_vrt_fn, mask_tmp_fn, options=vrt_options_specific_bands)
                for i, b in enumerate(good_bands):                
                    v_band = vrt_ds.GetRasterBand(i+1)
                    m_band = m_ds.GetRasterBand(b)
                    v_band.SetDescription(m_band.GetDescription())
                    v_band.SetMetadata(m_band.GetMetadata())

                pbar = tqdm(desc='writing vrt to disk', total=100, leave=self.verbose)
                pbar_update = lambda a,b,c: pbar.update((a*100)-pbar.n)        
                msk_ds = gdal.GetDriverByName(fmt).CreateCopy(mask_fn, vrt_ds, 0, options=msk_opts, callback=pbar_update)
                pbar.close()
                vrt_ds = None
                utils.remove_glob(mask_tmp_fn, mask_vrt_fn)
            else:
                msk_ds = m_ds
                if self.verbose:
                    utils.echo_msg('no bands found for {}'.format(mask_fn))

            return(msk_ds)
        
        utils.set_cache(self.cache_dir)
        #mask_level = utils.int_or(mask_level, 0)
        # if self.stack_mode not in ['mean', 'min', 'max', 'supercede', 'mixed']:
        #     mode = 'mean'
        # else:
        #     mode = self.stack_mode
        #mode = self.stack_mode.split(':')[0]
        mode = self.stack_mode_name
        if 'mask_level' in self.stack_mode_args.keys():
            mask_level = utils.int_or(self.stack_mode_args['mask_level'], 0)
        else:
            mask_level = 0

        if self.verbose:
            utils.echo_msg('stacking using {} mode with mask level of {}'.format(mode, mask_level))
        
        ## initialize the output rasters
        if out_name is None:
            out_name = os.path.join(self.cache_dir, '{}'.format(
                utils.append_fn('dlim_stacks', self.region, self.x_inc)
            ))

        xcount, ycount, dst_gt = self.region.geo_transform(
            x_inc=self.x_inc, y_inc=self.y_inc, node='grid'
        )
        if xcount <= 0 or ycount <=0:
            utils.echo_error_msg(
                'could not create grid of {}x{} cells with {}/{} increments on region: {}'.format(
                    xcount, ycount, self.x_inc, self.y_inc, self.region
                )
            )
            sys.exit(-1)

        ## initialize stack grid
        out_file = '{}.{}'.format(out_name, gdalfun.gdal_fext(fmt))
        driver = gdal.GetDriverByName(fmt)
        if os.path.exists(out_file):
            status = driver.Delete(out_file)
            if status != 0:
                utils.remove_glob('{}*'.format(out_file))

        if not os.path.exists(os.path.dirname(out_file)):
            os.makedirs(os.path.dirname(out_file))
                
        gdt = gdal.GDT_Float32
        c_gdt = gdal.GDT_Int32
        stack_opts = ['COMPRESS=LZW',
                     'PREDICTOR=2',
                     'TILED=YES']
        if fmt == 'GTiff':
            stack_opts.append('BIGTIFF=YES')
        elif fmt == 'MEM':
            stack_opts = []
        
        dst_ds = driver.Create(
            out_file,
            xcount,
            ycount,
            7,
            gdt,
            options=stack_opts
        )
        if dst_ds is None:
            utils.echo_error_msg(
                'failed to create stack grid {} {} {} {} {}...'.format(
                    out_file, xcount, ycount, gdt, fmt
                )
            )
            sys.exit(-1)

        dst_ds.SetGeoTransform(dst_gt)
        if self.dst_srs is not None:
            dst_ds.SetProjection(gdalfun.osr_wkt(self.dst_srs))

        ## the stack keys are the output bands in the stack grid
        ## stacked_data holds tempory data from incoming entries
        ## stacked_bands holds the gdal_ds band for each of the stakc keys
        stack_keys = ['z', 'count', 'weights', 'uncertainty', 'src_uncertainty', 'x', 'y']
        stacked_data = {x: None for x in stack_keys}            
        stacked_bands = {x: dst_ds.GetRasterBand(i+1) for i, x in enumerate(stack_keys)}
        for key in stack_keys:
            stacked_bands[key].SetNoDataValue(np.nan)
            stacked_bands[key].SetDescription(key)

        ## initialize data mask grid
        if self.want_mask:
            mask_fn = '{}_msk.{}'.format(out_name, gdalfun.gdal_fext(fmt))
            mask_tmp_fn = '{}_tmpmsk.{}'.format(out_name, gdalfun.gdal_fext(fmt))
            mask_vrt_fn = '{}_msk.vrt'.format(out_name)
            if os.path.exists(mask_fn):
                status = driver.Delete(mask_fn)
                if status != 0:
                    utils.remove_glob('{}*'.format(mask_fn))                

            if not os.path.exists(os.path.dirname(mask_fn)):
                os.makedirs(os.path.dirname(mask_fn))

            m_gdt = gdal.GDT_Byte
            msk_opts = ['COMPRESS=LZW']
            if fmt == 'GTiff':
                msk_opts.append('BIGTIFF=YES')
            elif fmt == 'MEM':
                msk_opts = []

            driver = gdal.GetDriverByName('MEM')
            m_ds = driver.Create('', xcount, ycount, 1, m_gdt)
            m_ds.SetDescription('_msk'.format(out_name))
            m_ds.SetGeoTransform(dst_gt)
            if self.dst_srs is not None:
                m_ds.SetProjection(gdalfun.osr_wkt(self.dst_srs))

            ## set first band to hold the data mask for all the data
            m_band_all = m_ds.GetRasterBand(1)
            m_band_all.SetNoDataValue(0)
            m_band_all.SetDescription('Full Data Mask')
        
        ## parse each entry and process it
        for this_entry in self.parse():
            ## only check to add a new band per entry
            #m_band = add_mask_band(m_ds, this_entry, mask_level=mask_level)                
            ## yield entry arrays for stacks
            ##
            ## incoming arrays arrs['z'], arrs['weight'] arrs['uncertainty'], and arrs['count']
            ##
            ## srcwin is the srcwin of the waffle, represented by the gt, relative to the incoming arrays
            for arrs, srcwin, gt in this_entry.array_yield:
                ## update the mask
                if np.all(arrs['count'] == 0):
                    continue

                if self.want_mask:
                    ## check to add a new band if the arrs has data
                    m_all_array = m_band_all.ReadAsArray(srcwin[0], srcwin[1], srcwin[2], srcwin[3])
                    m_band = add_mask_band(m_ds, this_entry, mask_level=mask_level)
                    if m_band is not None:
                        m_array = m_band.ReadAsArray(srcwin[0], srcwin[1], srcwin[2], srcwin[3])
                        #m_array[arrs['count'] != 0] = 1
                        #m_band.WriteArray(m_array, srcwin[0], srcwin[1])

                    #m_all_array[arrs['count'] != 0] = 1
                    #m_band_all.WriteArray(m_all_array, srcwin[0], srcwin[1])
                    #m_ds.FlushCache()
                    #m_array = m_all_array = None
                    # if mask_only:
                    #     continue
                    ## mrl moved masking to mode sections
                    
                ## Read the saved accumulated rasters at the incoming srcwin and set ndv to zero
                for key in stack_keys:
                    stacked_data[key] = stacked_bands[key].ReadAsArray(
                        srcwin[0], srcwin[1], srcwin[2], srcwin[3]
                    )
                    if mode != 'min' and mode != 'max':
                        stacked_data[key][np.isnan(stacked_data[key])] = 0
                    #else:
                    if key == 'count':
                        stacked_data[key][np.isnan(stacked_data[key])] = 0
                    
                ## set incoming np.nans to zero and mask to non-nan count
                arrs['count'][np.isnan(arrs['count'])] = 0
                arrs['weight'][np.isnan(arrs['z'])] = 0
                arrs['uncertainty'][np.isnan(arrs['z'])] = 0

                #arrs['weight'][arrs['count']==0] = 0
                #arrs['uncertainty'][arrs['count']==0] = 0

                if mode != 'min' and mode != 'max':
                    arrs['x'][np.isnan(arrs['x'])] = 0
                    arrs['y'][np.isnan(arrs['y'])] = 0
                    arrs['z'][np.isnan(arrs['z'])] = 0
                    for arr_key in arrs:
                        if arrs[arr_key] is not None:
                            arrs[arr_key][np.isnan(arrs[arr_key])] = 0

                ## add the count to the accumulated rasters
                if mode != 'mixed':
                    stacked_data['count'] += arrs['count']

                ## supercede based on weights, else do weighted mean
                ## todo: do (weighted) mean on cells with same weight
                if mode == 'supercede':
                    sup_mask = arrs['weight'] > (stacked_data['weights'] / stacked_data['count'])
                    if self.want_mask:
                        # remove the mask from other bands
                        reset_mask_bands(m_ds, srcwin, except_band_name=m_band.GetDescription(), mask=sup_mask)
                        m_array[(sup_mask) & (arrs['count'] != 0)] = 1
                        m_all_array[(sup_mask) & (arrs['count'] != 0)] = 1
                    
                    ## higher weight supercedes lower weight (first come first served atm)
                    stacked_data['z'][(sup_mask)] = arrs['z'][(sup_mask)]
                    stacked_data['x'][(sup_mask)] = arrs['x'][(sup_mask)]
                    stacked_data['y'][(sup_mask)] = arrs['y'][(sup_mask)]
                    stacked_data['src_uncertainty'][(sup_mask)] = arrs['uncertainty'][(sup_mask)]
                    stacked_data['weights'][(sup_mask)] = arrs[(sup_mask)]
                    ## uncertainty is src_uncertainty, as only one point goes into a cell
                    stacked_data['uncertainty'][(sup_mask)] = np.array(stacked_data['src_uncertainty'][(sup_mask)])

                elif mode == 'mixed':
                    ## weights above threshold supercede weights below threshold, otherwise meaned...
                    wt = 1
                    if 'weight_threshold' in self.stack_mode_args.keys():
                        wt = utils.float_or(self.stack_mode_args['weight_threshold'], 1)

                    # above
                    tmp_stacked_weight = (stacked_data['weights'] / stacked_data['count'])
                    tmp_stacked_weight[np.isnan(tmp_stacked_weight)] = 0

                    # supercede existing data below weight_threshold
                    weight_above_sup = (arrs['weight'] >= wt) & (tmp_stacked_weight < wt)
                    #if np.any(weight_above_sup):
                    if self.want_mask:
                        # remove the mask from superceded bands
                        reset_mask_bands(m_ds, srcwin, except_band_name=m_band.GetDescription(), mask=weight_above_sup)
                        # for b in range(1, m_ds.RasterCount+1):
                        #     this_band = m_ds.GetRasterBand(b)
                        #     this_band_name = this_band.GetDescription()
                            
                        #     if this_band_name != m_band.GetDescription():
                        #         this_band_array = this_band.ReadAsArray(srcwin[0], srcwin[1], srcwin[2], srcwin[3])
                        #         masked = this_band_array[(weight_above_sup)] == 1
                        #         if np.any(masked):
                        #             this_band_array[(weight_above_sup)][(masked)] = 0                        
                        #             this_band.WriteArray(this_band_array, srcwin[0], srcwin[1])
                        #             m_ds.FlushCache()
                        #             this_band_array = None

                        #     this_band = None
                        
                        m_array[(weight_above_sup) & (arrs['count'] != 0)] = 1
                        m_all_array[(weight_above_sup) & (arrs['count'] != 0)] = 1

                    stacked_data['count'][weight_above_sup] = arrs['count'][weight_above_sup]
                    stacked_data['z'][weight_above_sup] = (arrs['z'][weight_above_sup] * arrs['weight'][weight_above_sup])
                    stacked_data['x'][weight_above_sup] = (arrs['x'][weight_above_sup] * arrs['weight'][weight_above_sup])
                    stacked_data['y'][weight_above_sup] = (arrs['y'][weight_above_sup] * arrs['weight'][weight_above_sup])
                    stacked_data['src_uncertainty'][weight_above_sup] = arrs['uncertainty'][weight_above_sup]
                    stacked_data['weights'][weight_above_sup] = arrs['weight'][weight_above_sup]
                    stacked_data['uncertainty'][weight_above_sup] = np.array(stacked_data['src_uncertainty'][weight_above_sup])
                    
                    tmp_stacked_weight = (stacked_data['weights'] / stacked_data['count'])
                    tmp_stacked_weight[np.isnan(tmp_stacked_weight)] = 0
                        
                    # average of incoming data with existing data above weight_threshold
                    weight_above = (arrs['weight'] >= wt) & (arrs['weight'] >= tmp_stacked_weight) & (~weight_above_sup)
                    if self.want_mask:
                        m_array[(weight_above) & (arrs['count'] != 0)] = 1
                        m_all_array[(weight_above) & (arrs['count'] != 0)] = 1
                        
                    stacked_data['count'][weight_above] += arrs['count'][weight_above]
                    stacked_data['z'][weight_above] += (arrs['z'][weight_above] * arrs['weight'][weight_above])
                    stacked_data['x'][weight_above] += (arrs['x'][weight_above] * arrs['weight'][weight_above])
                    stacked_data['y'][weight_above] += (arrs['y'][weight_above] * arrs['weight'][weight_above])
                    stacked_data['src_uncertainty'][weight_above] = np.sqrt(np.power(stacked_data['src_uncertainty'][weight_above], 2) \
                                                                            + np.power(arrs['uncertainty'][weight_above], 2))
                    stacked_data['weights'][weight_above] += arrs['weight'][weight_above]
                    ## accumulate variance * weight
                    stacked_data['uncertainty'][weight_above] += arrs['weight'][weight_above] \
                        * np.power(
                            (arrs['z'][weight_above] - (stacked_data['z'][weight_above] / stacked_data['weights'][weight_above])),
                            2
                        )
                    
                    # below
                    tmp_stacked_weight = (stacked_data['weights'] / stacked_data['count'])
                    tmp_stacked_weight[np.isnan(tmp_stacked_weight)] = 0
                    weight_below = (arrs['weight'] < wt) & (tmp_stacked_weight < wt)#& (arrs['weight'] >= tmp_stacked_weight)
                    if self.want_mask:
                        m_array[(weight_below) & (arrs['count'] != 0)] = 1
                        m_all_array[(weight_below) & (arrs['count'] != 0)] = 1
                        
                    stacked_data['count'][weight_below] += arrs['count'][weight_below]
                    stacked_data['z'][weight_below] += (arrs['z'][weight_below] * arrs['weight'][weight_below])
                    stacked_data['x'][weight_below] += (arrs['x'][weight_below] * arrs['weight'][weight_below])
                    stacked_data['y'][weight_below] += (arrs['y'][weight_below] * arrs['weight'][weight_below])
                    stacked_data['src_uncertainty'][weight_below] = np.sqrt(np.power(stacked_data['src_uncertainty'][weight_below], 2) \
                                                                            + np.power(arrs['uncertainty'][weight_below], 2))
                    stacked_data['weights'][weight_below] += arrs['weight'][weight_below]
                    ## accumulate variance * weight
                    stacked_data['uncertainty'][weight_below] += arrs['weight'][weight_below] \
                        * np.power(
                            (arrs['z'][weight_below] - (stacked_data['z'][weight_below] / stacked_data['weights'][weight_below])),
                            2
                        )               

                elif mode == 'min' or mode == 'max':
                    ## set nodata values in stacked_data to whatever the value is in arrs
                    if self.want_mask:
                        m_array[arrs['count'] != 0] = 1
                        m_all_array[arrs['count'] != 0] = 1

                    mask = np.isnan(stacked_data['z'])
                    stacked_data['x'][mask] = arrs['x'][mask]
                    stacked_data['y'][mask] = arrs['y'][mask]
                    stacked_data['src_uncertainty'][mask] = arrs['uncertainty'][mask]
                    stacked_data['uncertainty'][mask] = arrs['uncertainty'][mask]
                    stacked_data['weights'][mask] = arrs['weight'][mask]
                    stacked_data['z'][mask] = arrs['z'][mask]

                    ## mask the min or max and apply it to stacked_data 
                    if mode == 'min':
                        mask = arrs['z'] <= stacked_data['z']
                    else:
                        mask = arrs['z'] >= stacked_data['z']
                        
                    stacked_data['x'][mask] = arrs['x'][mask]
                    stacked_data['y'][mask] = arrs['y'][mask]

                    ## accumulate uncertainty and weight
                    stacked_data['src_uncertainty'][mask] += (arrs['uncertainty'][mask] * arrs['weight'][mask])

                    ## accumulate incoming weights (weight*weight?) and set results to np.nan for calcs
                    stacked_data['weights'][mask] += arrs['weight'][mask]
                    stacked_data['weights'][mask][stacked_data['weights'][mask] == 0] = np.nan
                    
                    ## accumulate variance * weight
                    stacked_data['uncertainty'][mask] += arrs['weight'][mask] \
                        * np.power((arrs['z'][mask] - (stacked_data['z'][mask] / stacked_data['weights'][mask])), 2)
                    stacked_data['z'][mask] = arrs['z'][mask]
                    
                elif mode == 'mean':
                    if self.want_mask:
                        m_array[arrs['count'] != 0] = 1
                        m_all_array[arrs['count'] != 0] = 1
                    
                    ## accumulate incoming z*weight and uu*weight
                    stacked_data['z'] += (arrs['z'] * arrs['weight'])
                    stacked_data['x'] += (arrs['x'] * arrs['weight'])
                    stacked_data['y'] += (arrs['y'] * arrs['weight'])
                    #stacked_data['src_uncertainty'] += (arrs['uncertainty'] * arrs['weight'])
                    stacked_data['src_uncertainty'] = np.sqrt(np.power(stacked_data['src_uncertainty'], 2) \
                                                              + np.power(arrs['uncertainty'], 2))
                    
                    ## accumulate incoming weights (weight*weight?) and set results to np.nan for calcs
                    stacked_data['weights'] += arrs['weight']
                    stacked_data['weights'][stacked_data['weights'] == 0] = np.nan
                    
                    ## accumulate variance * weight
                    stacked_data['uncertainty'] += arrs['weight'] \
                        * np.power((arrs['z'] - (stacked_data['z'] / stacked_data['weights'])), 2)

                if self.want_mask:
                    m_band.WriteArray(m_array, srcwin[0], srcwin[1])            
                    m_band_all.WriteArray(m_all_array, srcwin[0], srcwin[1])
                    m_ds.FlushCache()
                    m_array = m_all_array = None
                    
                ## write out results to accumulated rasters
                stacked_data['count'][stacked_data['count'] == 0] = np.nan

                for key in stack_keys:
                    stacked_data[key][np.isnan(stacked_data['count'])] = np.nan
                    stacked_bands[key].WriteArray(stacked_data[key], srcwin[0], srcwin[1])
                                    
        ## Finalize weighted mean rasters and close datasets
        ##
        ## incoming arrays have all been processed, if weighted mean the
        ## "z" is the sum of z*weight, "weights" is the sum of weights
        ## "uncertainty" is the sum of variance*weight
        if self.verbose:
            utils.echo_msg('finalizing stacked raster bands...')

        if self.want_mask:
            ## output the mask raster to disk
            pbar = tqdm(desc='writing mask to disk', total=100, leave=self.verbose)
            pbar_update = lambda a,b,c: pbar.update((a*100)-pbar.n)        
            msk_ds = gdal.GetDriverByName(fmt).CreateCopy(mask_fn, m_ds, 0, options=msk_opts, callback=pbar_update)
            pbar.close()

            ## if we added a band per entry, we have to go through and remove empty bands
            ## mask data will be checked for bands with nodata and those
            ## will be omitted from the final mask raster
            #msk_ds = remove_empty_mask_bands(m_ds)
            m_ds = None
            
        #if not mask_only:
        ## by scan-line
        srcwin = (0, 0, dst_ds.RasterXSize, dst_ds.RasterYSize)
        for y in range(
                srcwin[1], srcwin[1] + srcwin[3], 1
        ):
            for key in stack_keys:
                stacked_data[key] = stacked_bands[key].ReadAsArray(srcwin[0], y, srcwin[2], 1)
                stacked_data[key][stacked_data[key] == ndv] = np.nan

            stacked_data['weights'] = stacked_data['weights'] / stacked_data['count']
            #if mode == 'mean' or mode == 'min' or mode == 'max' or mode == 'mixed':
            if mode != 'supercede':
                if mode == 'mean' or mode == 'mixed':
                    ## average the accumulated arrays for finalization
                    ## x, y, z and u are weighted sums, so divide by weights
                    stacked_data['x'] = (stacked_data['x'] / stacked_data['weights']) / stacked_data['count']
                    stacked_data['y'] = (stacked_data['y'] / stacked_data['weights']) / stacked_data['count']
                    stacked_data['z'] = (stacked_data['z'] / stacked_data['weights']) / stacked_data['count']

                ## apply the source uncertainty with the sub-cell variance uncertainty
                ## caclulate the standard error (sqrt( uncertainty / count))
                stacked_data['uncertainty'] = np.sqrt((stacked_data['uncertainty'] / stacked_data['weights']) / stacked_data['count'])
                stacked_data['uncertainty'] = np.sqrt(np.power(stacked_data['src_uncertainty'], 2) + np.power(stacked_data['uncertainty'], 2))

            ## write out final rasters
            for key in stack_keys:
                stacked_data[key][np.isnan(stacked_data[key])] = ndv
                stacked_bands[key].WriteArray(stacked_data[key], srcwin[0], y)

        ## set the final output nodatavalue
        for key in stack_keys:
            stacked_bands[key].DeleteNoDataValue()

        for key in stack_keys:
            stacked_bands[key].SetNoDataValue(ndv)
                
        ## create a vector of the masks (spatial-metadata)
        if self.want_mask and self.want_sm:
            polygonize_mask_multibands(
                msk_ds, output=os.path.basename(out_name), verbose=False, mask_level=0#mask_level
            )
            # ogr_mask_footprints(                
            #     msk_ds, output=os.path.basename(out_name), verbose=False, mask_level=mask_level
            # )
            
        msk_ds = dst_ds = None
        if self.want_mask:
            os.replace(mask_fn, os.path.basename(mask_fn))
        
        ## apply any grits filters to the stack
        for f in self.stack_fltrs:
            utils.echo_msg('filtering stacks module with {}'.format(f))
            grits_filter = grits.GritsFactory(
                mod=f, src_dem=out_file, uncertainty_mask=4, weight_mask=3, count_mask=2
            )._acquire_module()
            if grits_filter is not None:
                grits_filter()
                os.replace(grits_filter.dst_dem, out_file)

        if self.verbose:
            utils.echo_msg('generated stack: {}'.format(os.path.basename(out_file)))
            if self.want_mask:
                utils.echo_msg('generated mask: {}'.format(os.path.basename(mask_fn)))
                           
        return(out_file)
    
    def _stacks_h5(self, out_name = None, ndv = -9999, mask_level = 0):
        """stack and mask incoming arrays (from `self.array_yield`) together

        -----------
        Parameters:
        out_name (str): the output stacked raster basename
        ndv (float): the desired no data value
        mask_level (int): the granularity of the mask, 0 is every file

        --------
        Returns:
        out_name.csg - an hdf5 file containing the stack and mask, with datasets:
        stack:
          x
          y
          z
          weights
          count
          uncertainty
          src uncertainty
        mask:
          full_data_mask
          ...<dataset mask>...
        """

        utils.set_cache(self.cache_dir)
        # if self.stack_mode not in ['mean', 'min', 'max', 'supercede']:
        #     mode = 'mean'
        # else:
        #     mode = self.stack_mode

        if self.verbose:
            utils.echo_msg('stacking using {} mode'.format(self.stack_mode))

        mask_level = utils.int_or(mask_level, 0)
        ## initialize the output rasters
        if out_name is None:
            out_name = os.path.join(self.cache_dir, '{}'.format(
                utils.append_fn('dlim_stacks', self.region, self.x_inc)
            ))

        ## .csg is h5 Cudem Stack Grid
        out_file = '{}.csg'.format(out_name)
        stack_fn = '{}.{}'.format(out_name, gdalfun.gdal_fext(fmt))
        mask_fn = '{}_msk.{}'.format(out_name, gdalfun.gdal_fext(fmt))        
        if os.path.exists(out_file):
            utils.remove_glob('{}.*'.format(out_file))
            
        if not os.path.exists(os.path.dirname(out_file)):
            os.makedirs(os.path.dirname(out_file))
            
        xcount, ycount, dst_gt = self.region.geo_transform(
            x_inc=self.x_inc, y_inc=self.y_inc, node='grid'
        )
        if xcount <= 0 or ycount <=0:
            utils.echo_error_msg(
                'could not create grid of {}x{} cells with {}/{} increments on region: {}'.format(
                    xcount, ycount, self.x_inc, self.y_inc, self.region
                )
            )
            sys.exit(-1)

        lon_start = dst_gt[0] + (dst_gt[1] / 2)
        lat_start = dst_gt[3] + (dst_gt[5] / 2)
        lon_end = dst_gt[0] + (dst_gt[1] * xcount)
        lat_end = dst_gt[3] + (dst_gt[5] * ycount)
        lon_inc = dst_gt[1]
        lat_inc = dst_gt[5]
        stack_ds = h5.File(out_file, 'w', rdcc_nbytes=(1024**2)*12000)#, rdcc_nslots=1e7)
        
        # ************  Grid Mapping  ********
        crs_dset = stack_ds.create_dataset('crs', dtype=h5.string_dtype())
        crs_dset.attrs['GeoTransform'] = ' '.join([str(x) for x in dst_gt])
        if self.dst_srs is not None:
            crs_dset.attrs['crs_wkt'] = self.dst_srs
        
        # ************  LATITUDE  ************
        lat_array = np.arange(lat_start, lat_end, lat_inc)
        lat_dset = stack_ds.create_dataset ('lat', data=lat_array)
        lat_dset.make_scale('latitude')
        lat_dset.attrs["long_name"] = "latitude"
        lat_dset.attrs["units"] = "degrees_north"
        lat_dset.attrs["standard_name"] = "latitude"
        lat_dset.attrs["actual_range"] = [lat_end, lat_start]
        #lat_dset.attrs[''] = ''

        # ************  LONGITUDE  ***********
        lon_array = np.arange(lon_start, lon_end, lon_inc)
        lon_dset = stack_ds.create_dataset ('lon', data=lon_array)
        lon_dset.make_scale('longitude')
        lon_dset.attrs["long_name"] = "longitude"
        lon_dset.attrs["units"] = "degrees_east" 
        lon_dset.attrs["standard_name"]= "longitude"
        lon_dset.attrs["actual_range"] = [lon_start, lon_end]

        # ************  STACK  ***************
        stack_grp = stack_ds.create_group('stack')
        stacked_data = {
            'z': None,
            'count': None,
            'weights': None,
            'uncertainty': None,
            'src_uncertainty': None,
            'x': None,
            'y': None
        }

        for key in stacked_data.keys():
            stack_dset = stack_grp.create_dataset(
                key, data=np.full((ycount, xcount), np.nan),
                compression='lzf', maxshape=(ycount, xcount),
                chunks=(100, xcount), rdcc_nbytes=1024*xcount*400,
                dtype=np.float32
            )
            stack_dset.dims[0].attach_scale(lat_dset)
            stack_dset.dims[1].attach_scale(lon_dset)
            stack_dset.attrs['grid_mapping'] = "crs"
        
        # ************  MASK  ****************   
        mask_grp = stack_ds.create_group('mask')
        
        ## set first dataset to hold the data mask for all the data
        mask_all_dset = mask_grp.create_dataset(
            'full_data_mask', data=np.zeros((ycount,xcount)),
            compression='lzf', maxshape=(ycount, xcount),
            chunks=(100, xcount), rdcc_nbytes=1024*xcount*100,
            dtype=np.uint8
        )
        mask_all_dset.dims[0].attach_scale(lat_dset)
        mask_all_dset.dims[1].attach_scale(lon_dset)
        mask_all_dset.attrs['grid_mapping'] = "crs"
        
        ## parse each entry and process it
        ## todo: mask here instead of in each dataset module
        for this_entry in self.parse():
            ## MASK
            if mask_level < 0:
                mask_dset = None
            else:
                entry_name = this_entry.metadata['name']
                if mask_level > 0:
                    if mask_level > len(entry_name.split('/')):
                        mask_level = len(entry_name.split('/')) - 2

                    entry_name = '/'.join(entry_name.split('/')[:-mask_level])

                if not entry_name in mask_grp.keys():
                    mask_dset = mask_grp.create_dataset(
                        entry_name, data=np.zeros((ycount,xcount)),
                        compression='lzf', maxshape=(ycount, xcount),
                        chunks=(1, xcount), rdcc_nbytes=1024*xcount,
                        dtype=np.uint8
                    )
                    
                    mask_dset.dims[0].attach_scale(lat_dset)
                    mask_dset.dims[1].attach_scale(lon_dset)
                    mask_dset.attrs['grid_mapping'] = "crs"
                else:
                    mask_dset = mask_grp[key]
                    
                for k in this_entry.metadata.keys():
                    if k not in mask_dset.attrs.keys() or mask_dst.attrs[k] is None:
                        try:
                            mask_dset.attrs[k] = str(this_entry.metadata[k])
                        except:
                            utils.echo_warning_msg(k)
                            utils.echo_warning_msg(this_entry.metadata[k])
                            mask_dset.attrs[k] = None
                    else:
                        mask_dset.attrs[k] = None

                mask_dset.attrs['weight'] = this_entry.weight if this_entry.weight is not None else 1
                mask_dset.attrs['uncertainty'] = this_entry.uncertainty if this_entry.uncertainty is not None else 0
                
            ## yield entry arrays for stacks
            ## incoming arrays arrs['z'], arrs['weight'] arrs['uncertainty'], and arrs['count']
            ## srcwin is the srcwin of the waffle relative to the incoming arrays
            ## gt is the geotransform of the incoming arrays
            for arrs, srcwin, gt in this_entry.array_yield:
                ## update the mask
                mask_all_dset[srcwin[1]:srcwin[1]+srcwin[3], srcwin[0]:srcwin[0]+srcwin[2]][arrs['count'] != 0] = 1
                if mask_dset is not None:
                    mask_dset[srcwin[1]:srcwin[1]+srcwin[3], srcwin[0]:srcwin[0]+srcwin[2]][arrs['count'] != 0] = 1

                # if mask_only:
                #     continue
                
                ## Read the saved accumulated rasters at the incoming srcwin and set ndv to zero
                for key in stack_grp.keys():
                    stacked_data[key] = stack_grp[key][srcwin[1]:srcwin[1]+srcwin[3], srcwin[0]:srcwin[0]+srcwin[2]]
                    if self.stack_mode != 'min' and self.stack_mode != 'max':
                        stacked_data[key][np.isnan(stacked_data[key])] = 0

                    if key == 'count':
                        stacked_data[key][np.isnan(stacked_data[key])] = 0
                    
                ## set incoming np.nans to zero and mask to non-nan count
                arrs['count'][np.isnan(arrs['count'])] = 0
                arrs['weight'][np.isnan(arrs['z'])] = 0
                arrs['uncertainty'][np.isnan(arrs['z'])] = 0
                
                if self.stack_mode != 'min' and self.stack_mode != 'max':
                    arrs['x'][np.isnan(arrs['x'])] = 0
                    arrs['y'][np.isnan(arrs['y'])] = 0
                    arrs['z'][np.isnan(arrs['z'])] = 0
                    for arr_key in arrs:
                        if arrs[arr_key] is not None:
                            arrs[arr_key][np.isnan(arrs[arr_key])] = 0

                ## add the count to the accumulated rasters
                stacked_data['count'] += arrs['count']

                ## supercede based on weights, else do weighted mean
                ## todo: do (weighted) mean on cells with same weight
                if self.stack_mode == 'supercede':
                    ## higher weight supercedes lower weight (first come first served atm)
                    stacked_data['z'][arrs['weight'] > stacked_data['weights']] = arrs['z'][arrs['weight'] > stacked_data['weights']]
                    stacked_data['x'][arrs['weight'] > stacked_data['weights']] = arrs['x'][arrs['weight'] > stacked_data['weights']]
                    stacked_data['y'][arrs['weight'] > stacked_data['weights']] = arrs['y'][arrs['weight'] > stacked_data['weights']]
                    stacked_data['src_uncertainty'][arrs['weight'] > stacked_data['weights']] = arrs['uncertainty'][arrs['weight'] > stacked_data['weights']]
                    stacked_data['weights'][arrs['weight'] > stacked_data['weights']] = arrs['weight'][arrs['weight'] > stacked_data['weights']]
                    ## uncertainty is src_uncertainty, as only one point goes into a cell
                    stacked_data['uncertainty'][:] = np.array(stacked_data['src_uncertainty'])

                elif self.stack_mode == 'min' or self.stack_mode == 'max':
                    ## set nodata values in stacked_data to whatever the value is in arrs
                    mask = np.isnan(stacked_data['z'])
                    stacked_data['x'][mask] = arrs['x'][mask]
                    stacked_data['y'][mask] = arrs['y'][mask]
                    stacked_data['src_uncertainty'][mask] = arrs['uncertainty'][mask]
                    stacked_data['uncertainty'][mask] = arrs['uncertainty'][mask]
                    stacked_data['weights'][mask] = arrs['weight'][mask]
                    stacked_data['z'][mask] = arrs['z'][mask]

                    ## mask the min or max and apply it to stacked_data 
                    if self.stack_mode == 'min':
                        mask = arrs['z'] <= stacked_data['z']
                    else:
                        mask = arrs['z'] >= stacked_data['z']
                        
                    stacked_data['x'][mask] = arrs['x'][mask]
                    stacked_data['y'][mask] = arrs['y'][mask]

                    ## accumulate uncertainty and weight
                    stacked_data['src_uncertainty'][mask] += (arrs['uncertainty'][mask] * arrs['weight'][mask])

                    ## accumulate incoming weights (weight*weight?) and set results to np.nan for calcs
                    stacked_data['weights'][mask] += arrs['weight'][mask]
                    stacked_data['weights'][mask][stacked_data['weights'][mask] == 0] = np.nan
                    
                    ## accumulate variance * weight
                    stacked_data['uncertainty'][mask] += arrs['weight'][mask] \
                        * np.power((arrs['z'][mask] - (stacked_data['z'][mask] / stacked_data['weights'][mask])), 2)
                    stacked_data['z'][mask] = arrs['z'][mask]
                    
                elif self.stack_mode == 'mean':
                    ## accumulate incoming z*weight and uu*weight
                    stacked_data['z'] += (arrs['z'] * arrs['weight'])
                    stacked_data['x'] += (arrs['x'] * arrs['weight'])
                    stacked_data['y'] += (arrs['y'] * arrs['weight'])
                    stacked_data['src_uncertainty'] = np.sqrt(np.power(stacked_data['src_uncertainty'], 2) + np.power(arrs['uncertainty'], 2))
                    
                    ## accumulate incoming weights (weight*weight?) and set results to np.nan for calcs
                    stacked_data['weights'] += arrs['weight']
                    stacked_data['weights'][stacked_data['weights'] == 0] = np.nan
                    
                    ## accumulate variance * weight
                    stacked_data['uncertainty'] += arrs['weight'] * np.power((arrs['z'] - (stacked_data['z'] / stacked_data['weights'])), 2)

                ## write out results to accumulated rasters
                stacked_data['count'][stacked_data['count'] == 0] = np.nan
                for key in stack_grp.keys():
                    stacked_data[key][np.isnan(stacked_data['count'])] = np.nan
                    stack_grp[key][srcwin[1]:srcwin[1]+srcwin[3], srcwin[0]:srcwin[0]+srcwin[2]] = stacked_data[key]
                                    
        ## Finalize weighted mean rasters and close datasets
        ## incoming arrays have all been processed, if weighted mean the
        ## "z" is the sum of z*weight, "weights" is the sum of weights
        ## "uncertainty" is the sum of variance*weight
        if self.verbose:
            utils.echo_msg('finalizing stacked raster bands...')

        for key in mask_grp.keys():
            key_dset = mask_grp[key]
            if np.all(key_dset == 0):
                del key_dset
                
        #if not mask_only:
        ## by scan-line
        # srcwin = (0, 0, xcount, ycount)
        # for y in range(
        #         srcwin[1], srcwin[1] + srcwin[3], 1
        # ):
        #     utils.echo_msg(y)
        for key in stack_grp.keys():
            stacked_data[key] = stack_grp[key][...]
            stacked_data[key][stacked_data[key] == ndv] = np.nan

        #utils.echo_msg(stacked_data['count'])
        if self.stack_mode == 'mean' or self.stack_mode == 'min' or self.stack_mode == 'max':
            stacked_data['weights'] = stacked_data['weights'] / stacked_data['count']
            if self.stack_mode == 'mean':
                ## average the accumulated arrays for finalization
                ## x, y, z and u are weighted sums, so divide by weights
                stacked_data['x'] = (stacked_data['x'] / stacked_data['weights']) / stacked_data['count']
                stacked_data['y'] = (stacked_data['y'] / stacked_data['weights']) / stacked_data['count']
                stacked_data['z'] = (stacked_data['z'] / stacked_data['weights']) / stacked_data['count']

            ## apply the source uncertainty with the sub-cell variance uncertainty
            ## caclulate the standard error (sqrt( uncertainty / count))
            stacked_data['uncertainty'] = np.sqrt((stacked_data['uncertainty'] / stacked_data['weights']) / stacked_data['count'])
            stacked_data['uncertainty'] = np.sqrt(np.power(stacked_data['src_uncertainty'], 2) + np.power(stacked_data['uncertainty'], 2))

        ## write out final rasters
        for key in stack_grp.keys():
            stacked_data[key][np.isnan(stacked_data[key])] = np.nan
            stack_grp[key][:] = stacked_data[key]

        stack_ds.close()                        
        return(stack_fn)

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
    
    def yield_points(self):
        """yield the numpy xyz rec array points from the dataset.

        reset in dataset if needed.
        """

        for ds in self.parse():
            for points in ds.yield_points():
                yield(points)

    def transform_and_yield_points(self):
        """points are an array of `points['x']`, `points['y']`, `points['z']`, <`points['w']`, `points['u']`>`

        points will be transformed here, based on `self.transformer`, which is set in `set_transform`
        after the points are transformed, the data will pass through the region, if it exists and finally
        yield the transformed and reduced points.
        """            
        
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            for points in self.yield_points():

                if self.transform['transformer'] is not None or self.transform['vert_transformer'] is not None:
                    #transformer = self.transform['transformer'] if self.transform['transformer'] is not None else self.transform['vert_transformer']
                    if self.transform['transformer'] is not None:
                        points['x'], points['y'] = self.transform['transformer'].transform(
                            points['x'], points['y']
                        )
                    if self.transform['vert_transformer'] is not None:
                        _, _, points['z'] = self.transform['vert_transformer'].transform(
                            points['x'], points['y'], points['z']
                        )
                        
                    points = points[~np.isinf(points['z'])]

                if self.region is not None and self.region.valid_p():
                    xyz_region = self.region.copy() #if self.transform['trans_region'] is None else self.transform['trans_region'].copy()
                    if self.invert_region:
                        points = points[((points['x'] >= xyz_region.xmax) | (points['x'] <= xyz_region.xmin)) | \
                                        ((points['y'] >= xyz_region.ymax) | (points['y'] <= xyz_region.ymin))]
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
                        points = points[((points['x'] <= xyz_region.xmax) & (points['x'] >= xyz_region.xmin)) & \
                                        ((points['y'] <= xyz_region.ymax) & (points['y'] >= xyz_region.ymin))]
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
                    if self.pnt_fltrs is not None:
                        for f in self.pnt_fltrs:
                            point_filter = PointFilterFactory(
                                mod=f, points=points, verbose=True
                            )._acquire_module()
                            if point_filter is not None:
                                points = point_filter()

                    if len(points) > 0:
                        yield(points)
        
        self.transform['transformer'] = self.transform['vert_transformer'] = None
            
    def mask_and_yield_array(self):
        """mask the incoming array from `self.yield_array` and yield the results.
        
        the mask should be either an ogr supported vector or a gdal supported raster.
        """
        
        mask_band = None
        mask_infos = None
        data_mask = None
        if self.mask is not None:
            if self.mask['ogr_or_gdal'] == 1: # mask is ogr, rasterize it
                if self.region is not None and self.x_inc is not None and self.y_inc is not None:
                    data_mask = gdalfun.ogr2gdal_mask(
                        self.mask['mask'], region=self.region, x_inc=self.x_inc, y_inc=self.y_inc, dst_srs=self.dst_srs,
                        invert=True, verbose=False, temp_dir=self.cache_dir
                    )
                    #self.mask['ogr_or_gdal'] = 0
            else:
                data_mask = self.mask['mask']
            
            if os.path.exists(data_mask):            
                utils.echo_msg('using mask dataset: {} to array'.format(data_mask))
                if self.region is not None and self.x_inc is not None and self.y_inc is not None:
                    src_mask = gdalfun.sample_warp(
                        data_mask, None, self.x_inc, self.y_inc,
                        src_region=self.region,
                        sample_alg='nearest',
                        dst_srs=self.transform['dst_horz_crs'] if self.transform['dst_horz_crs'] is not None else None,
                        #dst_srs=self.transform['dst_horz_crs'].to_proj4() if self.transform['dst_horz_crs'] is not None else None,
                        #dst_srs=self.dst_srs.split('+')[0],
                        ndv=gdalfun.gdal_get_ndv(self.mask['mask']),
                        verbose=True,
                        co=["COMPRESS=DEFLATE", "TILED=YES"]
                    )[0]
                else:
                    src_mask = gdal.Open(data_mask)
                    
                mask_band = src_mask.GetRasterBand(1)
                mask_infos = gdalfun.gdal_infos(src_mask)
            else:
                utils.echo_warning_msg('could not load mask {}'.format(data_mask))

        for out_arrays, this_srcwin, this_gt in self.yield_array():
            #utils.echo_msg(mask_band)
            if mask_band is not None:
                ycount, xcount = out_arrays['z'].shape
                this_region = regions.Region().from_geo_transform(this_gt, xcount, ycount)
                mask_data = mask_band.ReadAsArray(*this_srcwin)
                if not np.isnan(mask_infos['ndv']):
                    mask_data[mask_data==mask_infos['ndv']] = np.nan

                for arr in out_arrays.keys():
                    if out_arrays[arr] is not None:
                        if self.mask['invert_mask']:
                            if arr == 'count':
                                out_arrays[arr][~np.isnan(mask_data)] = 0
                            else:                            
                                out_arrays[arr][~np.isnan(mask_data)] = np.nan                            
                        else:
                            if arr == 'count':
                                out_arrays[arr][np.isnan(mask_data)] = 0
                            else:
                                out_arrays[arr][np.isnan(mask_data)] = np.nan

                mask_data = None
                
            yield(out_arrays, this_srcwin, this_gt)
        
        src_mask = mask_band = None
        if data_mask is not None:
            utils.remove_glob('{}*'.format(data_mask))

    def mask_and_yield_xyz(self):
        """mask the incoming xyz data from `self.yield_xyz` and yield the results.
        
        the mask should be either an ogr supported vector or a gdal supported raster.
        """

        for this_entry in self.parse():
            if this_entry.mask is None:
                for this_xyz in this_entry.yield_xyz():
                    yield(this_xyz)
            else:
                utils.echo_msg('using mask dataset: {} to xyz'.format(this_entry.mask['mask']))
                if this_entry.mask['ogr_or_gdal'] == 0:
                    src_ds = gdal.Open(this_entry.mask['mask'])
                    if src_ds is not None:
                        ds_config = gdalfun.gdal_infos(src_ds)
                        ds_band = src_ds.GetRasterBand(1)
                        ds_gt = ds_config['geoT']
                        ds_nd = ds_config['ndv']
                    
                        for this_xyz in this_entry.yield_xyz():
                            xpos, ypos = utils._geo2pixel(this_xyz.x, this_xyz.y, ds_gt, node='pixel')
                            if xpos < ds_config['nx'] and ypos < ds_config['ny'] and xpos >=0 and ypos >=0:
                                tgrid = ds_band.ReadAsArray(xpos, ypos, 1, 1)
                                if tgrid is not None:
                                    if tgrid[0][0] == ds_nd:
                                        if this_entry.mask['invert_mask']:
                                            yield(this_xyz)
                                    else:
                                        if not this_entry.mask['invert_mask']:
                                            yield(this_xyz)
                                            
                        src_ds = ds_band = None
                else:
                    ## this is very slow! find another way.
                    src_ds = ogr.Open(this_entry.mask['mask'])
                    if src_ds is not None:
                        layer = src_ds.GetLayer()
                        #utils.echo_msg(len(layer))

                        geomcol = gdalfun.ogr_union_geom(layer)
                        if not geomcol.IsValid():
                            geomcol = geomcol.Buffer(0)

                        for this_xyz in this_entry.yield_xyz():
                            this_wkt = this_xyz.export_as_wkt()
                            this_point = ogr.CreateGeometryFromWkt(this_wkt)
                            if this_point.Within(geomcol):
                                if not this_entry.mask['invert_mask']:
                                    yield(this_xyz)
                            else:
                                if not this_entry.mask['invert_mask']:
                                    yield(this_xyz)

                        src_ds = ds_band = None
                    else:
                        yield(this_xyz)
                                
    def yield_xyz(self):
        """Yield the data as xyz points

        incoming data are numpy rec-arrays of x,y,z<w,u> points.

        if 'w' (weight) is not set by the dataset, they will be given the
        default value from `self.weight`, otherwise the incoming weight values will be
        multiplied by `self.weight`

        if 'u' (uncertainty) is not set by the dataset, they will be given the
        default value from `self.uncertainty`, otherwise the incoming uncertainty values
        will be caluculated by `np.sqrt(u**2 + self.uncertainty**2)`
        """
        
        count = 0
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

            points_u = np.sqrt(points_u**2 + (self.uncertainty if self.uncertainty is not None else 0)**2)
            points_u[np.isnan(points_u)] = 0
            dataset = np.vstack((points['x'], points['y'], points['z'], points_w, points_u)).transpose()
            count += len(dataset)
            points = None
            for point in dataset:
                this_xyz = xyzfun.XYZPoint(x=point[0], y=point[1], z=point[2], w=point[3], u=point[4])
                yield(this_xyz)

        if self.verbose:
            utils.echo_msg_bold(
                'parsed {} data records from {}{}'.format(
                    count, self.fn, ' @{}'.format(self.weight) if self.weight is not None else ''
                )
            )
                     
    def yield_array(self):
        """Yield the data as an array which coincides with the desired region, x_inc and y_inc

        incoming data are numpy rec-arrays of x,y,z<w,u> points

        if 'w' (weight) is not set by the dataset, they will be given the
        default value from `self.weight`, otherwise the incoming weight values will be
        multiplied by `self.weight`

        if 'u' (uncertainty) is not set by the dataset, they will be given the
        default value from `self.uncertainty`, otherwise the incoming uncertainty values
        will be caluculated by `np.sqrt(u**2 + self.uncertainty**2)`
        """
        
        out_arrays = {
            'z':None,
            'count':None,
            'weight':None,
            'uncertainty': None,
            'mask':None,
            'x': None,
            'y': None
        }
        count = 0
        for points in self.transform_and_yield_points():
            xcount, ycount, dst_gt = self.region.geo_transform(
                x_inc=self.x_inc, y_inc=self.y_inc, node='grid'
            )
            ## convert the points to pixels based on the geotransform
            ## and calculate the local srcwin of the points
            pixel_x = np.floor((points['x'] - dst_gt[0]) / dst_gt[1]).astype(int)
            pixel_y = np.floor((points['y'] - dst_gt[3]) / dst_gt[5]).astype(int)
            # points_x = dst_gt[0] + (pixel_x+0.5) * dst_gt[1] + (pixel_y+0.5) * dst_gt[2]
            # points_y = dst_gt[3] + (pixel_x+0.5) * dst_gt[4] + (pixel_y+0.5) * dst_gt[5]
            points_x = np.array(points['x'])
            points_y = np.array(points['y'])
            
            pixel_z = np.array(points['z'])
            try:
                pixel_w = np.array(points['w'])
            except:
                pixel_w = np.ones(pixel_z.shape)

            try:
                pixel_u = np.array(points['u'])
            except:
                pixel_u = np.zeros(pixel_z.shape)
                
            ## remove pixels that will break the srcwin
            out_idx = np.nonzero((pixel_x >= xcount) | (pixel_x < 0) | (pixel_y >= ycount) | (pixel_y < 0))
            pixel_x = np.delete(pixel_x, out_idx)
            pixel_y = np.delete(pixel_y, out_idx)
            pixel_z = np.delete(pixel_z, out_idx)
            pixel_w = np.delete(pixel_w, out_idx)
            pixel_u = np.delete(pixel_u, out_idx)
            points_x = np.delete(points_x, out_idx)
            points_y = np.delete(points_y, out_idx)
            if len(pixel_x) == 0 or len(pixel_y) == 0:
                continue
            
            points = None
            pixel_w[np.isnan(pixel_w)] = 1
            pixel_u[np.isnan(pixel_u)] = 0
            
            ## set the srcwin of the incoming points
            this_srcwin = (int(min(pixel_x)), int(min(pixel_y)),
                           int(max(pixel_x) - min(pixel_x))+1,
                           int(max(pixel_y) - min(pixel_y))+1)
            count += len(pixel_x)

            ## adjust the pixels to the srcwin and stack together
            pixel_x = pixel_x - this_srcwin[0]
            pixel_y = pixel_y - this_srcwin[1]
            pixel_xy = np.vstack((pixel_y, pixel_x)).T
            
            ## find the non-unique x/y points and mean/min/max their z values together
            ## while calculating the std for uncertainty
            unq, unq_idx, unq_inv, unq_cnt = np.unique(
                pixel_xy, axis=0, return_inverse=True, return_index=True, return_counts=True
            )
            cnt_msk = unq_cnt > 1
            cnt_idx, = np.nonzero(cnt_msk)
            idx_msk = np.in1d(unq_inv, cnt_idx)
            idx_idx, = np.nonzero(idx_msk)
            srt_idx = np.argsort(unq_inv[idx_msk])
            dup_idx = np.split(idx_idx[srt_idx], np.cumsum(unq_cnt[cnt_msk])[:-1])
            zz = pixel_z[unq_idx]
            ww = pixel_w[unq_idx]
            uu = pixel_u[unq_idx]
            xx = points_x[unq_idx]
            yy = points_y[unq_idx]
            #u = np.zeros(zz.shape)
            if np.any([len(dup) for dup in dup_idx]):
                if self.stack_mode == 'min':
                    dup_stack = [np.min(pixel_z[dup]) for dup in dup_idx]
                    dup_stds = np.zeros(dup_stack.shape)
                elif self.stack_mode == 'max':
                    dup_stack = [np.max(pixel_z[dup]) for dup in dup_idx]
                    dup_stds = np.zeros(dup_stack.shape)
                else:
                    dup_stack = [np.mean(pixel_z[dup]) for dup in dup_idx]
                    dup_stack_x = [np.mean(points_x[dup]) for dup in dup_idx]
                    dup_stack_y = [np.mean(points_y[dup]) for dup in dup_idx]
                    dup_stds = [np.std(pixel_z[dup]) for dup in dup_idx]
                    
                zz[cnt_msk] = dup_stack
                #uu[cnt_msk] = np.sqrt(dup_stds)
                uu[cnt_msk] = np.sqrt(np.power(uu[cnt_msk],2) + np.power(dup_stds,2))
                
            ## make the output arrays to yield
            out_x = np.zeros((this_srcwin[3], this_srcwin[2]))
            out_x[unq[:,0], unq[:,1]] = xx
            out_x[out_x == 0] = np.nan
            out_arrays['x'] = out_x

            out_y = np.zeros((this_srcwin[3], this_srcwin[2]))
            out_y[unq[:,0], unq[:,1]] = yy
            out_y[out_y == 0] = np.nan
            out_arrays['y'] = out_y
            
            out_z = np.zeros((this_srcwin[3], this_srcwin[2]))
            out_z[unq[:,0], unq[:,1]] = zz
            out_z[out_z == 0] = np.nan
            out_arrays['z'] = out_z
            
            out_arrays['count'] = np.zeros((this_srcwin[3], this_srcwin[2]))
            out_arrays['count'][unq[:,0], unq[:,1]] = unq_cnt
            
            out_arrays['weight'] = np.ones((this_srcwin[3], this_srcwin[2]))
            out_arrays['weight'][:] = self.weight if self.weight is not None else 1
            out_arrays['weight'][unq[:,0], unq[:,1]] *= (ww * unq_cnt)
            #out_arrays['weight'][unq[:,0], unq[:,1]] *= unq_cnt
            
            out_arrays['uncertainty'] = np.zeros((this_srcwin[3], this_srcwin[2]))
            #out_arrays['uncertainty'][:] = self.uncertainty if self.uncertainty is not None else 0
            out_arrays['uncertainty'][unq[:,0], unq[:,1]] = np.sqrt(uu**2 + (self.uncertainty if self.uncertainty is not None else 0)**2)                

            ## by scan-line for h5
            # for y in range(0, this_srcwin[3], 1):
            #     scanline_arrays = {}
            #     scanline_srcwin = (this_srcwin[0], this_srcwin[1]+y, this_srcwin[2], 1)
            #     for key in out_arrays.keys():
            #         if out_arrays[key] is not None:                        
            #             scanline_arrays[key] = out_arrays[key][y:y+1, 0:this_srcwin[2]]
            #         else:
            #             scanline_arrays[key] = out_arrays[key]

            #     yield(scanline_arrays, scanline_srcwin, dst_gt)
                        
            ## for gdal
            yield(out_arrays, this_srcwin, dst_gt)

        if self.verbose:
            utils.echo_msg_bold(
                'parsed {} data records from {}{}'.format(
                    count, self.fn, ' @{}'.format(self.weight) if self.weight is not None else ''
                )
            )    

    def stacks_yield_xyz(self, out_name = None, ndv = -9999, fmt = 'GTiff'):#, mode = 'mean'):
        """yield the result of `_stacks` as an xyz object"""

        stacked_fn = self._stacks(out_name=out_name, ndv=ndv, fmt=fmt)#, mode=mode)
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
                            x = sx[0, x], y = sy[0, x], z = z, w = sw[0, x], u = su[0, x]
                        )
                    else: # yield center of pixel as x/y
                        geo_x, geo_y = utils._pixel2geo(x, y, sds_gt)
                        out_xyz = xyzfun.XYZPoint(
                            x=geo_x, y=geo_y, z=z, w=sw[0, x], u=su[0, x]
                        )
                        
                    yield(out_xyz)
                    
            sw = su = sx = sy = None
                    
        sds = sds_z_band = sds_w_band = sds_u_band = sds_x_band = sds_y_band = None

    def stacks_yield_xyz_h5(self, out_name = None, ndv = -9999):#, mode = 'mean'):
        """yield the result of `_stacks` as an xyz object"""

        stacked_fn = self._stacks(out_name=out_name, ndv=ndv)#, mode=mode)
        sds = h5.File(stacked_fn, 'r')
        sds_gt = [float(x) for x in sds['crs'].attrs['GeoTransform'].split()]
        sds_stack = sds['stack']
        stack_shape = sds['stack']['z'].shape
        srcwin = (0, 0, stack_shape[1], stack_shape[0])
        for y in range(srcwin[1], srcwin[1] + srcwin[3], 1):
            sz = sds_stack['z'][y:y+1, srcwin[0]:srcwin[0]+srcwin[2]]
            ## skip row if all values are ndv
            if np.all(sz == ndv):
                continue

            sw = sds_stack['weights'][y:y+1, srcwin[0]:srcwin[0]+srcwin[2]]
            su = sds_stack['uncertainty'][y:y+1, srcwin[0]:srcwin[0]+srcwin[2]]

            sx = sds_stack['x'][y:y+1, srcwin[0]:srcwin[0]+srcwin[2]]
            sy = sds_stack['y'][y:y+1, srcwin[0]:srcwin[0]+srcwin[2]]
            for x in range(0, stack_shape[1]):
                z = sz[0, x]
                if z != ndv:
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
            
    ## Data Dump/Export/Archive
    def dump_xyz(self, dst_port=sys.stdout, encode=False):
        """dump the XYZ data from the dataset.

        data gets parsed through `self.xyz_yield`. See `set_yield` for more info.
        """

        for this_xyz in self.xyz_yield:
            this_xyz.dump(
                include_w=True if self.weight is not None else False,
                include_u=True if self.uncertainty is not None else False,
                dst_port=dst_port,
                encode=encode,
                precision=self.dump_precision
            )

    def dump_xyz_direct(self, dst_port=sys.stdout, encode=False):
        """dump the XYZ data from the dataset

        data get dumped directly from `self.yield_xyz`, by-passing `self.xyz_yield`.
        """
        
        for this_xyz in self.yield_xyz():
            this_xyz.dump(
                include_w=True if self.weight is not None else False,
                include_u=True if self.uncertainty is not None else False,
                dst_port=dst_port,
                encode=encode,
                precision=self.dump_precision
            )
            
    def export_xyz_as_list(self, z_only = False):
        """return the XYZ data from the dataset as python list

        This may get very large, depending on the input data.
        """

        return([xyz.z if z_only else xyz.copy() for xyz in self.xyz_yield])
    
    def export_xyz_as_ogr(self):
        """Make a point vector OGR DataSet Object from src_xyz

        for use in gdal gridding functions
        """

        dst_ogr = '{}'.format(self.metadata['name'])
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
        
        dataset = pd.DataFrame({}, columns=['x', 'y', 'z', 'weight', 'uncertainty'])
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
                points_u**2 + (self.uncertainty if self.uncertainty is not None else 0)**2
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
            dataset = pd.concat([dataset, points_dataset], ignore_index=True)

        return(dataset)

    ## TODO: fix/write these functions
    def export_stack_as_pandas(self):
        pass
    
    def export_stack_as_gdal(self, fmt='GTiff'):
        """export the stack h5 to gdal multi-band or seperate file(s)
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

        with tqdm(
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
            options=['COMPRESS=LZW', 'BIGTIFF=YES'] if fmt == 'GTiff' else ['COMPRESS=LZW']
        )
        mask_dataset.SetGeoTransform(dst_gt)
        if self.dst_srs is not None:
            mask_dataset.SetProjection(gdalfun.osr_wkt(self.dst_srs))

        with tqdm(
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
                    mask_grp[key][srcwin[1]:srcwin[1]+srcwin[3], srcwin[0]:srcwin[0]+srcwin[2]][np.isnan(mask_grp[key][srcwin[1]:srcwin[1]+srcwin[3], srcwin[0]:srcwin[0]+srcwin[2]])] = 0
                    #if np.all(mask_grp[key][mask_grp[key] == 0]):
                    #    continue

                    #ii += 1
                    mask_band = mask_dataset.GetRasterBand(i+1)
                    mask_band.WriteArray(mask_grp[key][srcwin[1]:srcwin[1]+srcwin[3], srcwin[0]:srcwin[0]+srcwin[2]])
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

    def archive_xyz(self, **kwargs):
        """Archive data from the dataset to XYZ in the given dataset region.
        
        will convert all data to XYZ within the given region and will arrange
        the data as a datalist based on inputs.

        Data comes from `self.xyz_yield` set in `self.set_yield`. So will process data through
        `_stacks` before archival if region, x_inc and y_inc are set.
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
            srs_all.append(this_entry.dst_srs if this_entry.dst_srs is not None else this_entry.src_srs)
            datalist_dirname = os.path.join(a_name, os.path.dirname(this_entry.metadata['name']))
            this_key = datalist_dirname.split('/')[-1]
            if not os.path.exists(datalist_dirname):
                os.makedirs(datalist_dirname)

            sub_datalist = os.path.join(datalist_dirname, '{}.datalist'.format(this_key))
            with open(sub_datalist, 'a') as sub_dlf:
                if utils.fn_basename2(os.path.basename(this_entry.fn)) == '':
                    sub_xyz_path = '.'.join([utils.fn_basename2(os.path.basename(this_entry.metadata['name'])), 'xyz'])
                else:
                    sub_xyz_path = '.'.join([utils.fn_basename2(os.path.basename(this_entry.fn)), 'xyz'])                        
                this_xyz_path = os.path.join(datalist_dirname, sub_xyz_path)
                if os.path.exists(this_xyz_path):
                    utils.echo_warning_msg('{} already exists, skipping...'.format(this_xyz_path))
                    continue

                with open(this_xyz_path, 'w') as xp:
                    for this_xyz in this_entry.xyz_yield: # data will be processed independently of each other
                        this_xyz.dump(include_w=True if self.weight is not None else False,
                                      include_u=True if self.uncertainty is not None else False,
                                      dst_port=xp, encode=False, precision=self.dump_precision)

                if os.stat(this_xyz_path).st_size > 0:
                    sub_dlf.write('{} 168 1 0\n'.format(sub_xyz_path))
                else:
                    utils.remove_glob('{}*'.format(this_xyz_path))

            if os.stat(sub_datalist).st_size == 0:
                utils.remove_glob(sub_datalist)
            else:
                with open(self.archive_datalist, 'a') as dlf:
                    if not this_key in archive_keys:
                        archive_keys.append(this_key)                    
                        dlf.write(
                            '{name}.datalist -1 {weight} {uncertainty} {metadata}\n'.format(
                                name=os.path.join(datalist_dirname, this_key),
                                weight = utils.float_or(this_entry.weight, 1),
                                uncertainty = utils.float_or(this_entry.uncertainty, 0),
                                metadata = this_entry.format_metadata()
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

        utils.remove_glob('{}.*'.format(self.archive_datalist))
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
        
    ## Fetching
    ## TODO: move these out of this class
    def fetch_water_temp(self):
        """fetch and return the average water temperature over the region"""
        
        if self.region is not None:
            this_region = self.region.copy()
        else:
            this_region = regions.Region().from_list(self.infos.minmax)
            
        time_start = self.f.attrs['time_coverage_start'].decode('utf-8')
        time_end = self.f.attrs['time_coverage_end'].decode('utf-8')
        this_sst = fetches.MUR_SST(
            src_region=this_region,
            verbose=self.verbose,
            outdir=self.cache_dir,
            time_start=time_start,
            time_end=time_end
        )
        this_sst.run()
        
        # sea_temp = ds.analysed_sst.sel(lat=lats, lon=lons)
        # sst = round(np.nanmedian(sea_temp.values)-273,2)
        
        return(20)
        # return(sst)
            
    def fetch_buildings(self, verbose = True):
        """fetch building footprints from BING"""
        
        if self.region is not None:
            this_region = self.region.copy()
        else:
            this_region = regions.Region().from_list(self.infos.minmax)

        if this_region.valid_p():
            if verbose:
                utils.echo_msg('fetching buildings for region {}'.format(this_region))
                
            this_bldg = fetches.BingBFP(
                src_region=this_region, verbose=self.verbose, outdir=self.cache_dir
            )
            this_bldg.run()
            fr = fetches.fetch_results(this_bldg)#, check_size=False)
            fr.daemon=True
            fr.start()
            fr.join()
            return(fr)
        
        return(None)

    def process_buildings(self, this_bing, verbose = True):
        bldg_geoms = []
        if this_bing is not None:
            with tqdm(
                    total=len(this_bing.results),
                    desc='processing buildings',
                    leave=verbose
            ) as pbar:
                for n, bing_result in enumerate(this_bing.results):                
                    if bing_result[-1] == 0:
                        pbar.update()
                        bing_gz = bing_result[1]
                        try:
                            bing_gj = utils.gunzip(bing_gz, self.cache_dir)
                            os.rename(bing_gj, bing_gj + '.geojson')
                            bing_gj = bing_gj + '.geojson'
                            bldg_ds = ogr.Open(bing_gj, 0)
                            bldg_layer = bldg_ds.GetLayer()
                            bldg_geom = gdalfun.ogr_union_geom(bldg_layer, verbose=False)
                            bldg_geoms.append(bldg_geom)
                            bldg_ds = None
                            utils.remove_glob(bing_gj)
                        except Exception as e:
                            utils.echo_error_msg('could not process bing bfp, {}'.format(e))

                        utils.remove_glob(bing_gz)

        return(bldg_geoms)

    def fetch_coastline(self, chunks = True, verbose = True):
        if self.region is not None:
            this_region = self.region.copy()
        else:
            this_region = regions.Region().from_list(self.infos.minmax)

        if this_region.valid_p():
            if verbose:
                utils.echo_msg('fetching coastline for region {}'.format(this_region))
                
            this_cst = fetches.OpenStreetMap(
                src_region=this_region, verbose=verbose, outdir=self.cache_dir,
                q='coastline', chunks=chunks,
            )
            this_cst.run()
            fr = fetches.fetch_results(this_cst)
            fr.daemon=True
            fr.start()
            fr.join()
            return(fr)
        
        return(None)

    def process_coastline(
            self, this_cst, return_geom = True, landmask_is_watermask = False,
            line_buffer = 0.0000001, include_landmask = False, verbose = True
    ):
        if self.region is not None:
            this_region = self.region.copy()
        else:
            this_region = regions.Region().from_list(self.infos.minmax)
            
        cst_geoms = []
        if this_cst is not None:
            with tqdm(
                    total=len(this_cst.results),
                    desc='processing coastline',
                    leave=verbose
            ) as pbar:
                for n, cst_result in enumerate(this_cst.results):                
                    if cst_result[-1] == 0:
                        pbar.update()
                        cst_osm = cst_result[1]
                        out = fetches.polygonize_osm_coastline(
                            cst_osm,
                            utils.make_temp_fn(
                                utils.fn_basename2(cst_osm) + '_coast.shp',
                                temp_dir=self.cache_dir
                            ),
                            region=this_region, include_landmask=include_landmask,
                            landmask_is_watermask=landmask_is_watermask,
                            line_buffer=line_buffer
                        )

                        if out is not None:
                            cst_ds = ogr.Open(out, 0)
                            cst_layer = cst_ds.GetLayer()
                            cst_geom = gdalfun.ogr_union_geom(cst_layer, verbose=False)
                            cst_geoms.append(cst_geom)
                            cst_ds = None
                            utils.remove_glob(cst_osm)
                        
            if return_geom:
                utils.remove_glob('{}.*'.format(utils.fn_basename2(out)))
                        
        if return_geom:            
            return(cst_geoms)
        else:
            return(out)
    
class XYZFile(ElevationDataset):
    """representing an ASCII xyz dataset stream.

    Parse data from an xyz file/stdin

    generate_inf - generate an inf file for the xyz data
    yield_xyz - yield the xyz data as xyz
    yield_array - yield the xyz data as an array must set the x_inc/y_inc in the super class
    
    -----------
    Parameters:
    
    delim: the delimiter of the xyz data (str)
    xpos: the position (int) of the x value
    ypos: the position (int) of the y value
    zpos: the position (int) of the z value
    wpos: the position (int) of the w value (weight)
    upos: the position (int) of the u value (uncertainty)
    skip: number of lines to skip
    x_scale: scale the x value
    y_scale: scale the y value
    z_scale: scale the z value
    x_offset: offset the x value
    y_offset: offset the y value    
    """
            
    def __init__(self, delim = None, xpos = 0, ypos = 1, zpos = 2,
                 wpos = None, upos = None, skip = 0, x_scale = 1, y_scale = 1,
                 z_scale = 1, x_offset = 0, y_offset = 0, use_numpy = True,
                 iter_rows = 1000000, **kwargs):
        
        super().__init__(**kwargs)
        self.delim = delim # the file delimiter
        self.xpos = utils.int_or(xpos, 0) # the position of the x/lat value
        self.ypos = utils.int_or(ypos, 1) # the position of the y/lon value
        self.zpos = utils.int_or(zpos, 2) # the position of the z value
        self.wpos = utils.int_or(wpos) # the position of the weight value
        self.upos = utils.int_or(upos) # the position of the uncertainty value
        self.skip = utils.int_or(skip, 0) # number of header lines to skip
        self.x_scale = utils.float_or(x_scale, 1) # multiply x by x_scale
        self.y_scale = utils.float_or(y_scale, 1) # multiply y by y_scale
        self.z_scale = utils.float_or(z_scale, 1) # multiply z by z_scale
        self.x_offset = x_offset # offset x by x_offset
        self.y_offset = utils.int_or(y_offset, 0) # offset y by y_offset
        self.rem = False # x is 360 instead of 180
        self.use_numpy = use_numpy # use numpy.loadtxt to load the xyz points
        self.iter_rows = iter_rows # max rows to process at a time

    def yield_points(self):
        if self.delim is not None:
            xyzfun._known_delims.insert(0, self.delim)
                
        if self.x_offset == 'REM':
            self.x_offset = 0
            self.rem = True
            
        self.scoff = True if self.x_scale != 1 or self.y_scale != 1 or self.z_scale != 1 \
           or self.x_offset != 0 or self.y_offset != 0 else False
        self.field_names  = [x for x in ['x' if self.xpos is not None else None,
                                         'y' if self.ypos is not None else None,
                                         'z' if self.zpos is not None else None,
                                         'w' if self.wpos is not None else None,
                                         'u' if self.upos is not None else None]
                             if x is not None]
        self.field_formats = [float for x in [self.xpos, self.ypos,
                                              self.zpos, self.wpos,
                                              self.upos] if x is not None]
        #if self.use_numpy:
        try:
            if self.delim is None:
                self.guess_delim()

            skip_ = self.skip            
            with open(self.fn, 'r') as src_data:
                while True:
                    points = np.loadtxt(
                        src_data,
                        delimiter=self.delim,
                        comments='#',
                        ndmin = 1,
                        skiprows=skip_,
                        usecols=[x for x in [self.xpos, self.ypos,
                                             self.zpos, self.wpos,
                                             self.upos] if x is not None],
                        dtype={'names': self.field_names, 'formats': self.field_formats},
                        max_rows=self.iter_rows
                    )
                    skip_ = 0

                    if self.scoff:
                        points['x'] = (points['x'] + self.x_offset) * self.x_scale
                        points['y'] = (points['y'] + self.y_offset) * self.y_scale
                        points['z'] *= self.z_scale

                    if self.rem:
                        points['x'] = np.fmod(points['x'] + 180, 360) - 180 
                    
                    yield(points)
                    
                    if self.iter_rows is None or len(points) < self.iter_rows:
                        break

        ## old processing function used as a fallback for when numpy.loadtxt fails
        except Exception as e:
            utils.echo_warning_msg(
                'could not load xyz data from {}, {}, falling back'.format(self.fn, e)
            )
            if self.fn is not None:
                if os.path.exists(str(self.fn)):
                    self.src_data = open(self.fn, "r")
                else:
                    self.src_data = self.fn
            else:
                self.src_data = sys.stdin

            points_x = []
            points_y = []
            points_z = []
            points_w = []
            points_u = []
            count = 0
            skip = self.skip
            for xyz_line in self.src_data:
                if count >= skip:
                    this_xyz = self.line_delim(xyz_line)
                    if this_xyz is None:
                        continue

                    if len(this_xyz) < 3:
                        continue
                    
                    x = utils.float_or(this_xyz[self.xpos])
                    y = utils.float_or(this_xyz[self.ypos])
                    z = utils.float_or(this_xyz[self.zpos])
                    w = utils.float_or(this_xyz[self.wpos]) if self.wpos is not None else 1
                    u = utils.float_or(this_xyz[self.upos]) if self.upos is not None else 0

                    if x is None or y is None or z is None:
                        continue
                    
                    if self.scoff:
                        x = (x+self.x_offset) * self.x_scale
                        y = (y+self.y_offset) * self.y_scale
                        z *= self.z_scale

                    if self.rem:
                        x = math.fmod(x+180,360)-180 

                    if self.data_region is not None and self.data_region.valid_p():
                        try:
                            this_xyz = xyzfun.XYZPoint(
                                x=this_xyz[self.xpos],
                                y=this_xyz[self.ypos],
                                z=this_xyz[self.zpos]
                            )
                        except Exception as e:
                            utils.echo_error_msg('{} ; {}'.format(e, this_xyz))
                            this_xyz = xyzfun.XYZPoint()
                            
                    points_x.append(x)
                    points_y.append(y)
                    points_z.append(z)
                    points_w.append(w)
                    points_u.append(u)
                    count += 1
                else:
                    skip -= 1

            self.src_data.close()
            dataset = np.column_stack((points_x, points_y, points_z, points_w, points_u))
            points = np.rec.fromrecords(dataset, names='x, y, z, w, u')
            
            yield(points)
        
    def guess_delim(self):
        """guess the xyz delimiter"""

        if self.fn is not None:
            if os.path.exists(str(self.fn)):
                self.src_data = open(self.fn, "r")
            else:
                self.src_data = self.fn
        else:
            self.src_data = sys.stdin

        for xyz_line in self.src_data:
            for delim in xyzfun._known_delims:
                try:
                    this_xyz = xyz_line.split(delim)
                    if len(this_xyz) > 1:
                        self.delim = delim
                        break
                except:
                    pass
            break
            
        self.src_data.close()

    def line_delim(self, xyz_line):
        """guess a line delimiter and return the split line."""

        for delim in xyzfun._known_delims:
            try:
                this_xyz = xyz_line.split(delim)
                if len(this_xyz) > 1:
                    return(this_xyz)
            except:
                pass

class YXZFile(XYZFile):
    """yxz file shortcut (167 mbdatalist)
    """
    
    def __init__(self, **kwargs):
        super().__init__(xpos = 1, ypos = 0, zpos = 2, **kwargs)
            
class LASFile(ElevationDataset):
    """representing an LAS/LAZ dataset.

    Process LAS/LAZ lidar files using pylas.
    
    get_epsg - attempt to parse the EPSG from the LAS file header
    generate_inf - generate an inf file for the LAS data
    yield_xyz - yield the LAS data as xyz
    yield_array - yield the LAS data as an array must set the x_inc/y_inc in the super class
    
    -----------
    Parameters:
    
    classes (str): a list of classes to parse, being a string with `/` seperator 
    """

    def __init__(self, classes='2/29/40', **kwargs):
        super().__init__(**kwargs)
        self.classes = [int(x) for x in classes.split('/')] # list of lidar classes to retain
        if self.src_srs is None:
            self.src_srs = self.get_epsg()
            if self.src_srs is None:
                self.src_srs = self.infos.src_srs
                                
    def valid_p(self, fmts = ['scratch']):
        """check if self appears to be a valid dataset entry"""

        if self.fn is None: # and not self.fn.startswith('http'):
            return(False)
        else:
            if os.path.exists(self.fn) :
                if os.stat(self.fn).st_size == 0:
                    return(False)
            else:
                return(False)

            try:
                lp.open(self.fn)
            except:
                utils.echo_warning_msg('{} could not be opened by the lasreader'.format(self.fn))
                return(False)
                        
        return(True)
        
    def get_epsg(self):
        with lp.open(self.fn) as lasf:
            lasf_vlrs = lasf.header.vlrs
            for vlr in lasf_vlrs:
                if vlr.record_id == 2112:
                    src_srs = vlr.string
                    return(src_srs)
            return(None)

    def generate_inf(self):
        """generate an inf file for a lidar dataset."""
        
        with lp.open(self.fn) as lasf:
            self.infos.numpts = lasf.header.point_count
            this_region = regions.Region(xmin=lasf.header.x_min, xmax=lasf.header.x_max,
                                         ymin=lasf.header.y_min, ymax=lasf.header.y_max,
                                         zmin=lasf.header.z_min, zmax=lasf.header.z_max)
            self.infos.minmax = this_region.export_as_list(include_z=True)
            self.infos.wkt = this_region.export_as_wkt()

        utils.echo_msg(self.get_epsg())
        self.infos.src_srs = self.src_srs if self.src_srs is not None else self.get_epsg()            
        return(self.infos)

    def yield_points(self):
        with lp.open(self.fn) as lasf:
            try:
                for points in lasf.chunk_iterator(2_000_000):
                    points = points[(np.isin(points.classification, self.classes))]
                    dataset = np.column_stack((points.x, points.y, points.z))
                    points = np.rec.fromrecords(dataset, names='x, y, z')                    
                    yield(points)
                    
            except Exception as e:
                utils.echo_warning_msg(
                    'could not read points from lasfile {}, {}'.format(self.fn, e)
                )
                
class GDALFile(ElevationDataset):
    """providing a GDAL raster dataset parser.

    Process/Parse GDAL supported raster files.
    See GDAL for more information regarding supported formats.
    
    generate_inf - generate an inf file for the RASTER data
    yield_xyz - yield the RASTER data as xyz
    yield_array - yield the RASTER data as an array

    -----------
    Parameters:

    mask: raster dataset to use as MASK dataset OR a band number
    weight_mask: raster dataset to use as weights (per-cell) OR a band number
    otherwise will use single value (weight) from superclass.
    uncertainty_mask: raster dataset to use as uncertainty (per-cell) OR a band number
    uncertainty_mask_to_meter: conversion value for uncertainty data to meter
    otherwise will use a single value (uncertainty) from superclass.
    open_options: GDAL open_options for raster dataset
    sample: sample method to use in resamplinig
    check_path: check to make sure path exists
    super_grid: Force processing of a supergrid (BAG files) (True/False)
    band_no: the band number of the elevation data
    remove_flat: remove flattened data from the input
    node: force node registration of either 'pixel' or 'grid'
    """
    
    def __init__(self,
                 weight_mask = None,
                 uncertainty_mask = None,
                 x_band = None,
                 y_band = None,
                 uncertainty_mask_to_meter = 1,
                 open_options = None,
                 sample = None,
                 check_path = True,
                 super_grid = False,
                 band_no = 1,
                 remove_flat = False,
                 node = None,
                 resample_and_warp = True,
                 **kwargs):
        super().__init__(**kwargs)
        self.weight_mask = weight_mask # associated raster file/band holding weight data
        self.uncertainty_mask = uncertainty_mask # associated raster file/band holding uncertainty data
        self.uncertainty_mask_to_meter = uncertainty_mask_to_meter # uncertainty data factor to meters
        self.open_options = open_options # GDAL open-options
        self.sample = sample # GDAL resampling method
        self.check_path = check_path # check for the path to the input data file
        self.super_grid = super_grid # input has super_grids (force their use)
        self.band_no = band_no # band number holding elevation data
        self.tmp_elev_band = None # temporary elevation band
        self.tmp_unc_band = None # temporary uncertainty band
        self.tmp_weight_band = None # temporary weight band
        self.src_ds = None # the GDAL dataset object
        self.remove_flat = remove_flat # remove flattened data from input
        self.flat_removed = False
        self.x_band = x_band # band holding x values
        self.y_band = y_band # band holding y values
        self.node = node # input is 'pixel' or 'grid' registered (force)
        self.resample_and_warp = resample_and_warp

        if self.fn.startswith('http') \
           or self.fn.startswith('/vsicurl/') \
           or self.fn.startswith('BAG'):
            self.check_path = False

        if self.valid_p() and self.src_srs is None:
            if self.infos.src_srs is None:
                self.src_srs = self.init_srs(self.fn)
            else:
                self.src_srs = self.infos.src_srs

    def destroy_ds(self):
        self.src_ds = None

    def init_srs(self, src_ds):
        """initialize the srs from the gdal file.

        try to split the horizontal and vertical and them combine them...
        """
        
        if self.src_srs is None:
            return(gdalfun.gdal_get_srs(src_ds))
        else:
            return(self.src_srs)

    def generate_inf(self):
        with gdalfun.gdal_datasource(self.fn) as src_ds:
            if src_ds is not None:

                if self.src_srs is None:
                    self.infos.src_srs = self.init_srs(src_ds)
                else:
                    self.infos.src_srs = self.src_srs

                ds_infos = gdalfun.gdal_infos(src_ds)
                this_region = regions.Region(src_srs=self.src_srs).from_geo_transform(
                    geo_transform=ds_infos['geoT'],
                    x_count=ds_infos['nx'],
                    y_count=ds_infos['ny']
                )

                #if scan:
                #    zr = src_ds.GetRasterBand(utils.int_or(self.band_no, 1)).ComputeRasterMinMax()

                ##gdalfun.gdal_polygonize(src_ds, dst_layer, verbose=True)
                
                #this_region.zmin, this_region.zmax = zr[0], zr[1]
                #self.infos.minmax = this_region.export_as_list(include_z=True)
                self.infos.minmax = this_region.export_as_list()
                self.infos.wkt = this_region.export_as_wkt()
                self.infos.numpts = ds_infos['nb']

        return(self.infos)
        
    def get_srcwin(self, gt, x_size, y_size, node = 'grid'):
        if self.region is not None:
            if self.invert_region:
                srcwin = (
                    0, 0, x_size, y_size, node
                )
            else:
                if self.transform['transformer'] is not None:
                    if self.transform['trans_region'] is not None \
                       and self.transform['trans_region'].valid_p(
                            check_xy = True
                    ):
                        srcwin = self.transform['trans_region'].srcwin(
                            gt, x_size, y_size, node
                        )
                    else:
                        srcwin = (
                            0, 0, x_size, y_size, node
                        )

                else:
                    srcwin = self.region.srcwin(
                        gt, x_size, y_size, node
                    )

        else:
            srcwin = (
                0, 0, x_size, y_size
            )
            
        return(srcwin)
        
    def yield_points(self):
        """initialize the raster dataset

        if x/y incs are set, will warp raster to that resolution.
        """

        if self.check_path and not os.path.exists(self.fn):
            return(None)

        ## apply any open_options that are specified
        try:
            self.open_options = self.open_options.split('/')
        except AttributeError:
            self.open_options = self.open_options
        except:
            self.open_options = None
            
        ## set up any transformations and other options
        inf_region = regions.Region().from_string(self.infos.wkt)
        self.sample_alg = self.sample if self.sample is not None else self.sample_alg
        self.dem_infos = gdalfun.gdal_infos(self.fn)
        if self.node is None:
            self.node = gdalfun.gdal_get_node(self.fn, 'pixel')
            
        #self.resample_and_warp = True        
        if (self.x_inc is None and self.y_inc is None) or self.region is None:
            self.resample_and_warp = False
        else:
            ## only perform 'upsamples' if the input raster is of higher resolution than the
            ## target self.x_inc, etc. then treat the data as points and do any 'downsampling'
            ## in stacks.
            inf_trans_region = inf_region.copy()
            inf_trans_region.src_srs = self.src_srs

            if self.dst_srs is not None:
                try:
                    inf_trans_region.warp(self.dst_srs)
                except Exception as e:
                    utils.echo_error_msg(
                        'could not perform transformation of region {} from {} to {}, {}'.format(
                            inf_trans_region, self.src_srs, self.dst_srs, e
                        )
                    )

            raster_is_higher_res = np.prod(
                inf_region.geo_transform(x_inc=self.dem_infos['geoT'][1])[:2]
            ) > np.prod(
                inf_trans_region.geo_transform(x_inc=self.x_inc)[:2]
            )
            if raster_is_higher_res:
                self.resample_and_warp = False
            
        if self.node == 'grid':
            self.resample_and_warp = False

        #self.resample_and_warp = False
        ndv = utils.float_or(gdalfun.gdal_get_ndv(self.fn), -9999)
        if self.region is not None:
            self.warp_region = self.region.copy()
        else:
            self.warp_region = inf_region.copy() #regions.Region().from_list(self.infos.minmax)
            if self.transform['transformer'] is not None:
                self.warp_region.src_srs = self.src_srs
                self.warp_region.warp(self.dst_srs)
  
        tmp_elev_fn = utils.make_temp_fn('{}'.format(self.fn), temp_dir=self.cache_dir)
        tmp_unc_fn = utils.make_temp_fn('{}'.format(self.fn), temp_dir=self.cache_dir)
        tmp_weight_fn = utils.make_temp_fn('{}'.format(self.fn), temp_dir=self.cache_dir)

        if self.remove_flat:
            grits_filter = grits.GritsFactory(
                mod='flats', src_dem=self.fn, cache_dir=self.cache_dir, verbose=True
            )._acquire_module()
            grits_filter()
            self.fn = grits_filter.dst_dem
            self.flat_removed = True
            
        ## resample/warp src gdal file to specified x/y inc/transformer respectively
        if self.resample_and_warp:
            if self.transform['transformer'] is not None:
                self.transform['transformer'] = None

            if self.sample_alg == 'auto':
                if self.stack_mode == 'min':
                    self.sample_alg = 'min'
                elif self.stack_mode == 'max':
                    self.sample_alg = 'max'
                elif not raster_is_higher_res:
                    self.sample_alg = 'bilinear'
                else:
                    self.sample_alg = 'average'

            tmp_ds = self.fn
            if self.open_options is not None:
                self.src_ds = gdal.OpenEx(self.fn, open_options=self.open_options)
                if self.src_ds is None:
                    self.src_ds = gdal.Open(self.fn)
            else:
                self.src_ds = gdal.Open(self.fn)

            if self.src_ds is None:
                return(None)

            ## Sample/Warp
            ## resmaple and/or warp dataset based on target srs and x_inc/y_inc
            ## doing this in MEM has a bug, fix if able
            ## extract necessary bands before warping, gadlwarp does not work with multiband rasters!
            tmp_warp = utils.make_temp_fn('{}'.format(tmp_ds), temp_dir=self.cache_dir)
            in_bands = self.src_ds.RasterCount
            src_ds_config = gdalfun.gdal_infos(self.src_ds)
            src_gt = src_ds_config['geoT']

            if in_bands > 1:
                ## the srcwin for to extract data
                if self.transform['trans_region'] is not None:
                    srcwin_region = self.transform['trans_region'].copy()
                elif self.region is not None:
                    srcwin_region = self.region.copy()
                else:
                    srcwin_region = None

                if srcwin_region is not None:
                    srcwin = srcwin_region.srcwin(
                        src_gt, src_ds_config['nx'], src_ds_config['ny'], node='pixel'
                    )
                else:
                    srcwin = None

                self.tmp_elev_band = tmp_elev_fn
                if self.verbose:
                    utils.echo_msg(
                        'extracting elevation data from {} to {}'.format(
                            self.fn, self.tmp_elev_band
                        )
                    )
                    
                self.tmp_elev_band, status = gdalfun.gdal_extract_band(
                    self.src_ds, self.tmp_elev_band,
                    band=self.band_no,
                    exclude=[],
                    srcwin=srcwin,
                    inverse=False
                )
                tmp_ds = self.tmp_elev_band
                if tmp_ds is None:
                    return(None)
                
                ## band is now 1!
                self.band_no = 1
                if utils.int_or(self.uncertainty_mask) is not None:
                    self.tmp_unc_band = tmp_unc_fn
                    if self.verbose:
                        utils.echo_msg(
                            'extracting uncertainty mask from {} to {}'.format(
                                self.fn, self.tmp_unc_band
                            )
                        )
                        
                    self.tmp_unc_band, status = gdalfun.gdal_extract_band(
                        self.src_ds, self.tmp_unc_band,
                        band=self.uncertainty_mask,
                        exclude=[],
                        srcwin=srcwin,
                        inverse=False
                    )
                    self.uncertainty_mask = self.tmp_unc_band

                if utils.int_or(self.weight_mask) is not None:                    
                    self.tmp_weight_band = tmp_weight_fn
                    if self.verbose:
                        utils.echo_msg(
                            'extracting weight mask from {} to {}'.format(
                                self.fn, self.tmp_weight_band
                            )
                        )
                        
                    self.tmp_weight_band, status = gdalfun.gdal_extract_band(
                        self.src_ds, self.tmp_weight_band,
                        band=self.weight_mask,
                        exclude=[],
                        srcwin=srcwin,
                        inverse=False
                    )
                    self.weight_mask = self.tmp_weight_band

            # if self.remove_flat:
            #     grits_filter = grits.GritsFactory(
            #         mod='flats', src_dem=tmp_ds, cache_dir=self.cache_dir, verbose=True
            #     )._acquire_module()
            #     grits_filter()
            #     tmp_ds = grits_filter.dst_dem
                
            warp_ = gdalfun.sample_warp(
                tmp_ds, tmp_warp, self.x_inc, self.y_inc,
                src_srs=self.transform['src_horz_crs'].to_proj4() if self.transform['src_horz_crs'] is not None else None,
                dst_srs=self.transform['dst_horz_crs'].to_proj4() if self.transform['dst_horz_crs'] is not None else None,                
                src_region=self.warp_region,
                sample_alg=self.sample_alg,
                ndv=ndv,
                verbose=self.verbose,
                co=["COMPRESS=DEFLATE", "TILED=YES"]
            )[0]
            tmp_ds = warp_ = None
            
            ## the following seems to be redundant...
            warp_ds = gdal.Open(tmp_warp)
            if warp_ds is not None:
                ## clip wapred ds to warped srcwin
                warp_ds_config = gdalfun.gdal_infos(warp_ds)
                gt = warp_ds_config['geoT']
                srcwin = self.warp_region.srcwin(gt, warp_ds.RasterXSize, warp_ds.RasterYSize, node='grid')
                dst_gt = (gt[0] + (srcwin[0] * gt[1]), gt[1], 0., gt[3] + (srcwin[1] * gt[5]), 0., gt[5])
                out_ds_config = gdalfun.gdal_set_infos(
                    srcwin[2], srcwin[3], srcwin[2] * srcwin[3], dst_gt, warp_ds_config['proj'],
                    warp_ds_config['dt'], warp_ds_config['ndv'], warp_ds_config['fmt'], None, None
                )

                in_bands = warp_ds.RasterCount
                self.src_ds = gdalfun.gdal_mem_ds(out_ds_config, bands=in_bands)
                if self.src_ds is not None:
                    for band in range(1, in_bands+1):
                        this_band = self.src_ds.GetRasterBand(band)
                        this_band.WriteArray(warp_ds.GetRasterBand(band).ReadAsArray(*srcwin))
                        self.src_ds.FlushCache()

                warp_ds = None
                utils.remove_glob(tmp_warp)
        else:
            # if self.remove_flat:
            #     grits_filter = grits.GritsFactory(
            #         mod='flats', src_dem=self.fn, cache_dir=self.cache_dir
            #     )._acquire_module()
            #     grits_filter()
            #     self.fn = grits_filter.dst_dem
            
            if self.open_options:
                self.src_ds = gdal.OpenEx(self.fn, open_options=self.open_options)
                if self.src_ds is None:
                    if self.verbose:
                        utils.echo_warning_msg(
                            'could not open file using open options {}'.format(
                                self.open_options
                            )
                        )
                        
                    self.src_ds = gdal.Open(self.fn)
            else:
                self.src_ds = gdal.Open(self.fn)
                
        self.src_dem_infos = gdalfun.gdal_infos(self.src_ds)
        if self.src_ds is None:
            if self.verbose:
                utils.echo_error_msg('could not load raster file {}'.format(self.fn))

        band = self.src_ds.GetRasterBand(utils.int_or(self.band_no, 1))
        gt = self.src_ds.GetGeoTransform()
        ndv = utils.float_or(band.GetNoDataValue())
        mask_band = None
        weight_band = None
        uncertainty_band = None
        nodata = ['{:g}'.format(-9999), 'nan', float('nan')]            
        if ndv is not None:
            nodata.append('{:g}'.format(ndv))

        band_data = count_data = weight_data = mask_data = None
                
        src_dem_x_inc = self.src_dem_infos['geoT'][1]
        src_dem_y_inc = -1*self.src_dem_infos['geoT'][5]
        src_dem_region = regions.Region().from_geo_transform(
            self.src_dem_infos['geoT'], self.src_dem_infos['nx'], self.src_dem_infos['ny']
        )
            
        ## todo: always warp these to src_ds
        ## weight mask, each cell should have the corresponding weight
        ## weight_mask can either be a seperate gdal file or a band number
        ## corresponding to the appropriate band in src_ds
        if self.weight_mask is not None:
            if utils.int_or(self.weight_mask) is not None:
                weight_band = self.src_ds.GetRasterBand(int(self.weight_mask))
            elif os.path.exists(self.weight_mask): # some numbers now return true here (file-descriptors), check for int first!
                src_weight = gdalfun.sample_warp(
                    self.weight_mask, None, src_dem_x_inc, src_dem_y_inc,
                    src_srs=self.transform['src_horz_crs'].to_proj4() if self.transform['src_horz_crs'] is not None else None,
                    dst_srs=self.transform['dst_horz_crs'].to_proj4() if self.transform['dst_horz_crs'] is not None else None,
                    src_region=src_dem_region,
                    sample_alg=self.sample_alg,
                    ndv=ndv,
                    verbose=self.verbose,
                    co=["COMPRESS=DEFLATE", "TILED=YES"]
                )[0]
                weight_band = src_weight.GetRasterBand(1)

            else:
                utils.echo_warning_msg('could not load weight mask {}'.format(self.weight_mask))
                weight_band = None

        ## uncertainty mask, each cell should have the corresponding uncertainty
        ## uncertainty_mask can either be a seperate gdal file or a band number
        ## corresponding to the appropriate band in src_ds
        if self.uncertainty_mask is not None:
            if utils.int_or(self.uncertainty_mask):
                uncertainty_band = self.src_ds.GetRasterBand(int(self.uncertainty_mask))
            elif os.path.exists(self.uncertainty_mask):
                src_uncertainty = gdalfun.sample_warp(
                    self.uncertainty_mask, None, src_dem_x_inc, src_dem_y_inc,
                    src_srs=self.transform['src_horz_crs'].to_proj4() if self.transform['src_horz_crs'] is not None else None,
                    dst_srs=self.transform['dst_horz_crs'].to_proj4() if self.transform['dst_horz_crs'] is not None else None,
                    src_region=src_dem_region,
                    sample_alg='bilinear',
                    ndv=ndv,
                    verbose=self.verbose,
                    co=["COMPRESS=DEFLATE", "TILED=YES"]
                )[0]
                uncertainty_band = src_uncertainty.GetRasterBand(1)
            else:
                utils.echo_warning_msg('could not load uncertainty mask {}'.format(self.uncertainty_mask))
                uncertainty_band = None

        ## uncertainty from the vertical transformation
        if self.transform['trans_fn_unc'] is not None:
            trans_uncertainty = gdalfun.sample_warp(
                self.transform['trans_fn_unc'], None, src_dem_x_inc, src_dem_y_inc,
                src_srs='+proj=longlat +datum=WGS84 +ellps=WGS84',
                dst_srs=self.transform['dst_horz_crs'].to_proj4() if self.transform['dst_horz_crs'] is not None else None,
                src_region=src_dem_region,
                sample_alg='bilinear',
                ndv=ndv,
                verbose=self.verbose,
                co=["COMPRESS=DEFLATE", "TILED=YES"]
            )[0]

            if uncertainty_band is not None:
                trans_uncertainty_band = trans_uncertainty.GetRasterBand(1)
                trans_uncertainty_arr = trans_uncertainty_band.ReadAsArray()
                uncertainty_arr = uncertainty_band.ReadAsArray()
                uncertainty_arr *= self.uncertainty_mask_to_meter
                uncertainty_arr = np.sqrt(uncertainty_arr**2 + trans_uncertainty_arr**2)
                uncertainty_band.WriteArray(uncertainty_arr)
                trans_uncertainty_band = None
            else:
                uncertainty_band = trans_uncertainty.GetRasterBand(1)

        ## parse through the data
        srcwin = self.get_srcwin(gt, self.src_ds.RasterXSize, self.src_ds.RasterYSize, node=self.node)
        for y in range(srcwin[1], (srcwin[1] + srcwin[3]), 1):
            band_data = band.ReadAsArray(srcwin[0], y, srcwin[2], 1).astype(float)
            if ndv is not None and not np.isnan(ndv):
                band_data[band_data == ndv] = np.nan

            if np.all(np.isnan(band_data)):
                continue

            ## weights
            if weight_band is not None:
                weight_data = weight_band.ReadAsArray(srcwin[0], y, srcwin[2], 1)
                weight_ndv = float(weight_band.GetNoDataValue())
                if not np.isnan(weight_ndv):
                    weight_data[weight_data==weight_ndv] = np.nan
            else:
                weight_data = np.ones(band_data.shape)

            ## uncertainty
            if uncertainty_band is not None:
                uncertainty_data = uncertainty_band.ReadAsArray(srcwin[0], y, srcwin[2], 1)
                uncertainty_ndv = float(uncertainty_band.GetNoDataValue())
                if not np.isnan(uncertainty_ndv):
                    uncertainty_data[uncertainty_data==uncertainty_ndv] = np.nan

                if self.transform['trans_fn_unc'] is None:
                    uncertainty_data *= self.uncertainty_mask_to_meter
                    
            else:
                uncertainty_data = np.zeros(band_data.shape)

            ## convert grid array to points
            if self.x_band is None and self.y_band is None:
                x_precision = 12#len(str(gt[0]).split('.')[-1])
                y_precision = 12#len(str(gt[3]).split('.')[-1])
                while True:
                    geo_x_origin, geo_y_origin = utils._pixel2geo(
                        srcwin[0], y, gt, node='pixel', x_precision=x_precision, y_precision=y_precision
                    ) # 'self.node to 'pixel', this breaks if set to 'grid' even if 'grid-node'
                    geo_x_end, geo_y_end = utils._pixel2geo(
                        srcwin[0] + srcwin[2], y, gt, node='grid', x_precision=x_precision, y_precision=y_precision
                    )
                    lon_array = np.arange(geo_x_origin, geo_x_end, gt[1])

                    if lon_array.shape == band_data[0].shape:
                        break
                    else:
                        if x_precision < 0 or y_precision < 0:
                            break

                        x_precision -= 1
                        y_precision -= 1

                lat_array = np.zeros((lon_array.shape))
                lat_array[:] = geo_y_origin
                dataset = np.column_stack((lon_array, lat_array, band_data[0], weight_data[0], uncertainty_data[0]))

                try:
                    assert lon_array.shape == lat_array.shape
                    assert lon_array.shape == band_data[0].shape
                    dataset = np.column_stack((lon_array, lat_array, band_data[0], weight_data[0], uncertainty_data[0]))
                except Exception as e:
                    utils.echo_error_msg(e)
                    pass
            else:
                lon_band = self.src_ds.GetRasterBand(self.x_band)
                lon_array = lon_band.ReadAsArray(srcwin[0], y, srcwin[2], 1).astype(float)
                lon_array[np.isnan(band_data)] = np.nan

                lat_band = self.src_ds.GetRasterBand(self.y_band)
                lat_array = lat_band.ReadAsArray(srcwin[0], y, srcwin[2], 1).astype(float)
                lat_array[np.isnan(band_data)] = np.nan
                dataset = np.column_stack((lon_array[0], lat_array[0], band_data[0], weight_data[0], uncertainty_data[0]))

            points = np.rec.fromrecords(dataset, names='x, y, z, w, u')
            points =  points[~np.isnan(points['z'])]
            dataset = band_data = weight_data = uncertainty_data = lat_array = lon_array = None
            utils.remove_glob(tmp_elev_fn, tmp_unc_fn, tmp_weight_fn)
            yield(points)
            
        src_uncertainty = src_weight = trans_uncertainty = self.src_ds = None
        # delete the filtered dem...
        if self.remove_flat and self.flat_removed:
            utils.remove_glob(self.fn)

## todo: update to h5py
class BAGFile(ElevationDataset):
    """providing a BAG raster dataset parser.

    Process supergrids at native resolution if they
    exist, otherwise process as normal grid.

    generate_inf - generate an inf file for the BAG data
    init_srs - discover the srs of the BAG file
    parse - parse the gdal datasets out of the BAG file

    -----------
    Parameters:

    explode (bool): Explode the BAG and process each super grid seperately.
    force_vr (bool): Force VR processing (if BAG file has bad header info)
    vr_strategy (str): VR strategy to use (MIN/MAX/AUTO)
    """

    def __init__(self,
                 explode = False,
                 force_vr = False,
                 vr_resampled_grid = True,
                 vr_interpolate = False,
                 vr_strategy = 'MIN',
                 **kwargs):
        super().__init__(**kwargs)
        self.explode = explode
        self.force_vr = force_vr
        self.vr_resampled_grid = vr_resampled_grid
        self.vr_interpolate = vr_interpolate
        if self.vr_interpolate:
            self.vr_resampled_grid = False
            
        self.vr_strategy = vr_strategy
        if self.src_srs is None:
            self.src_srs = self.init_srs()

    def init_srs(self):
        if self.src_srs is None:
            src_horz, src_vert = gdalfun.split_srs(gdalfun.gdal_get_srs(self.fn))
            if src_horz is None and src_vert is None:
                return(None)
            
            if src_vert is None:
                src_vert = '5866'

            self.src_srs = gdalfun.combine_epsgs(
                src_horz, src_vert, name='BAG Combined'
            )
            return(self.src_srs)
        else:
            return(self.src_srs)
                    
    def generate_inf(self):
        if self.src_srs is None:
            self.infos.src_srs = self.init_srs()
        else:
            self.infos.src_srs = self.src_srs

        with gdalfun.gdal_datasource(self.fn) as src_ds:
            if src_ds is not None:
                ds_infos = gdalfun.gdal_infos(src_ds)
                this_region = regions.Region(src_srs=self.src_srs).from_geo_transform(
                    geo_transform=ds_infos['geoT'],
                    x_count=ds_infos['nx'],
                    y_count=ds_infos['ny']
                )
                zr = src_ds.GetRasterBand(1).ComputeRasterMinMax() # bag band 1 is elevation                
                this_region.zmin, this_region.zmax = zr[0], zr[1] 
                self.infos.minmax = this_region.export_as_list(include_z=True)
                self.infos.wkt = this_region.export_as_wkt()
                self.infos.numpts = ds_infos['nb']
                
        return(self.infos)

    def parse(self, resample=True):
        mt = gdal.Info(self.fn, format='json')['metadata']['']
        oo = []
        if self.data_region is not None and self.data_region.valid_p():
            oo.append('MINX={}'.format(self.data_region.xmin))
            oo.append('MAXX={}'.format(self.data_region.xmax))
            oo.append('MINY={}'.format(self.data_region.ymin))
            oo.append('MAXY={}'.format(self.data_region.ymax))
                        
        if ('HAS_SUPERGRIDS' in mt.keys() and mt['HAS_SUPERGRIDS'] == 'TRUE') \
           or self.force_vr \
           or 'MAX_RESOLUTION_X' in mt.keys() \
           or 'MAX_RESOLUTION_Y' in mt.keys():
            if self.explode:
                oo.append("MODE=LIST_SUPERGRIDS")
                src_ds = gdal.OpenEx(self.fn, open_options=oo)
                sub_datasets = src_ds.GetSubDatasets()
                src_ds = None

                with tqdm(
                        total=len(sub_datasets),
                        desc='parsing {} supergrids from BAG file {}'.format(len(sub_datasets), self.fn),
                        leave=self.verbose
                ) as pbar:                
                    for sub_dataset in sub_datasets:
                        pbar.update()
                        sub_ds = DatasetFactory(
                            **self._set_params(
                                mod=sub_dataset[0],
                                data_format=200,
                                band_no=1,
                                uncertainty_mask_to_meter=0.01,
                                check_path=False,
                                super_grid=True,
                                uncertainty_mask=2)
                        )._acquire_module()
                        self.data_entries.append(sub_ds)
                        sub_ds.initialize()
                        for gdal_ds in sub_ds.parse():
                            yield(gdal_ds)

            elif self.vr_resampled_grid or self.vr_interpolate:
                if self.vr_resampled_grid:
                    oo.append("MODE=RESAMPLED_GRID")
                elif self.vr_interpolate:
                    oo.append("MODE=INTERPOLATED")
                    
                oo.append("RES_STRATEGY={}".format(self.vr_strategy))
                sub_ds = DatasetFactory(
                    **self._set_params(
                        mod=self.fn,
                        data_format=200,
                        band_no=1,
                        open_options=oo,
                        uncertainty_mask=2,
                        uncertainty_mask_to_meter=0.01
                    )
                )._acquire_module()
                self.data_entries.append(sub_ds)
                sub_ds.initialize()
                for gdal_ds in sub_ds.parse():
                    yield(gdal_ds)
            else: # use vrbag.py
                tmp_bag_as_tif = utils.make_temp_fn('{}_tmp.tif'.format(utils.fn_basename2(self.fn)))
                #self.x_inc * 111120 # scale cellsize to meters, todo: check if input is degress/meters/feet
                sr_cell_size = None
                vrbag.interpolate_vr_bag(
                    self.fn, tmp_bag_as_tif, self.cache_dir, sr_cell_size=sr_cell_size,
                    use_blocks=True, method='linear', nodata=3.4028234663852886e+38
                )
                sub_ds = DatasetFactory(
                    **self._set_params(
                        mod=tmp_bag_as_tif,
                        data_format=200,
                        band_no=1,
                        uncertainty_mask=2,
                        uncertainty_mask_to_meter=0.01
                    )
                )._acquire_module()
                self.data_entries.append(sub_ds)
                sub_ds.initialize()
                for gdal_ds in sub_ds.parse():
                    yield(gdal_ds)
                
        else:
            sub_ds = DatasetFactory(
                **self._set_params(
                    mod=self.fn,
                    data_format=200,
                    band_no=1,
                    uncertainty_mask=2,
                    uncertainty_mask_to_meter=0.01,
                )
            )._acquire_module()
            self.data_entries.append(sub_ds)
            sub_ds.initialize()
            for gdal_ds in sub_ds.parse():
                yield(gdal_ds)

class CUDEMFile(ElevationDataset):
    """CUDEM netcdf raster

    the cudem netcdf contains uninterpolated elevation data, uncertainty weight and data mask...
    """
    
    def __init__(self, stack=False, **kwargs):
        super().__init__(**kwargs)
        self.stack = stack


    def yield_points(self):

        ## extract z, unc and weight grids
        ## process through gdalfile
        
        nc_data = nc.Dataset(self.fn, 'r')
        stack_lat = nc_data['/stack/lat'][...,]
        stack_lon = nc_data['/stack/lon'][...,]
        stack_h = nc_data['/stack/stack_h'][...,]

        stack_u = nc_data['/stack/uncertainty'][...,]
        stack_w = nc_data['/stack/weight'][...,]
        counts = nc_data['/stack/count'][...,]

        nc_data.close()

        dataset = np.column_stack((stack_lon.data, stack_lat.data, stack_h.data[0], stack_w.data[0], stack_u.data[0]))
        points = np.rec.fromrecords(dataset, names='x, y, z, w, u')
        #points = points[points.z != -9999]
        
        yield(points)
        
## NASA SWOT Data class (hdf5)
## uses h5py
class SWOTFile(ElevationDataset):
    """NASA SWOT Data super class

    Uses h5py to parse data. Make a subclass of this to 
    process various types of SWOT data
    """
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        
    def _init_h5File(self, short_name='L2_HR_PIXC'):
        src_h5 = None
        try:
            src_h5 = h5.File(self.fn, 'r')
            if src_h5 is not None:
                if 'short_name' in src_h5.attrs.keys():
                    if src_h5.attrs['short_name'] != short_name.encode('utf-8'):
                        utils.echo_error_msg(
                            '{} does not appear to be a SWOT {} file'.format(
                                self.fn, short_name
                            )
                        )
                        self._close_h5File(src_h5)
                else:
                    utils.echo_error_msg(
                        '{} does not appear to be a SWOT file'.format(
                            self.fn
                        )
                    )
                    self._close_h5File(src_h5)
                    
        except Exception as e:
            utils.echo_error_msg(e)

        return(src_h5)
            
    def _close_h5File(self, src_h5):
        if src_h5 is not None:
            src_h5.close()

    def _get_var_arr(self, src_h5, var_path):
        return(src_h5['/{}'.format(var_path)][...,])
       
class SWOT_PIXC(SWOTFile):
    """NASA SWOT PIXC data file.

    Extract data from a SWOT PIXC file.

    classes: 1UB, 2UB, 3UB, 4UB, 5UB, 6UB, 7UB
    "land, land_near_water, water_near_land, open_water, dark_water, low_coh_water_near_land, open_low_coh_water"

    classes_qual: 1U, 2U, 4U, 8U, 16U, 2048U, 8192U, 16384U, 32768U, 262144U, 524288U, 134217728U, 536870912U, 1073741824U, 2147483648U
    "no_coherent_gain power_close_to_noise_floor detected_water_but_no_prior_water detected_water_but_bright_land 
    water_false_detection_rate_suspect coherent_power_suspect tvp_suspect sc_event_suspect small_karin_gap in_air_pixel_degraded 
    specular_ringing_degraded coherent_power_bad tvp_bad sc_event_bad large_karin_gap"

    anc_classes: 0UB, 1UB, 2UB, 3UB, 4UB, 5UB, 6UB 
    "open_ocean land continental_water aquatic_vegetation continental_ice_snow floating_ice salted_basin"  		
    """
    
    def __init__(self,
                 group = 'pixel_cloud',
                 var = 'height',
                 apply_geoid = True,
                 classes = None,
                 classes_qual = None,
                 anc_classes = None,
                 remove_class_flags = False,
                 **kwargs):
        super().__init__(**kwargs)
        self.group = group
        self.var = var
        self.apply_geoid = apply_geoid
        self.classes = [int(x) for x in classes.split('/')] if classes is not None else []
        self.classes_qual = [int(x) for x in classes_qual.split('/')] if classes_qual is not None else []
        self.anc_classes = [int(x) for x in anc_classes.split('/')] if anc_classes is not None else []
        self.remove_class_flags = remove_class_flags
        # if self.remove_class_flags:
        #     self.classes_qual = [1, 2, 4, 8, 16, 2048, 8192, 16384, 32768,
        #                          262144, 524288, 134217728, 536870912,
        #                          1073741824, 2147483648]

    def yield_points(self):
        src_h5 = self._init_h5File(short_name='L2_HR_PIXC')
        src_h5_vec = None
        
        #if self.pixc_vec is not None:
        #    src_h5_vec = self._init_h5File(short_name='L2_HR_PIXCVec')
        
        if src_h5 is not None:
            latitude = self._get_var_arr(src_h5, '{}/latitude'.format(self.group))
            longitude = self._get_var_arr(src_h5, '{}/longitude'.format(self.group))
            var_data = self._get_var_arr(src_h5, '{}/{}'.format(self.group, self.var))
            if self.apply_geoid:
                geoid_data = self._get_var_arr(src_h5, '{}/geoid'.format(self.group))
                out_data = var_data - geoid_data
            else:
                out_data = var_data
                
            dataset = np.column_stack((longitude, latitude, out_data))
            points = np.rec.fromrecords(dataset, names='x, y, z')
            #points = points[points['z'] != 9.96921e+36]

            ## Classification Filter
            if len(self.classes) > 0:
                class_data = self._get_var_arr(src_h5, '{}/classification'.format(self.group))
                points = points[(np.isin(class_data, self.classes))]

                ## Classification Quality Filter
                if self.remove_class_flags:
                    class_qual_data = self._get_var_arr(src_h5, '{}/classification_qual'.format(self.group))
                    class_qual_data = class_qual_data[(np.isin(class_data, self.classes))]
                    points = points[class_qual_data == 0]
                                   
                elif len(self.classes_qual) > 0:
                    class_qual_data = self._get_var_arr(src_h5, '{}/classification_qual'.format(self.group))
                    class_qual_data = class_qual_data[(np.isin(class_data, self.classes))]
                    points = points[(~np.isin(class_qual_data, self.classes_qual))]
                
            ## Ancilliary Classification Filter
            if len(self.anc_classes) > 0:
                anc_class_data = self._get_var_arr(
                    src_h5, '{}/ancillary_surface_classification_flag'.format(self.group)
                )
                points = points[(np.isin(anc_class_data, self.anc_classes))]
                
            points = points[points['z'] != -9.969209968386869e+36]
            self._close_h5File(src_h5)
            self._close_h5File(src_h5_vec)

            #tmp_points = points
            #while len(points) > 0:
            #tmp_points = points[:1000000]
            #points = points[1000000:]
            #utils.echo_msg(len(points))
            #yield(tmp_points)
            yield(points)

## todo: update to h5
class SWOT_HR_Raster(ElevationDataset):
    """NASA SWOT HR_Raster data file.

    Extract data from a SWOT HR_Raster file.
    """
        
    def __init__(self, data_set='wse', **kwargs):
        super().__init__(**kwargs)
        self.data_set = data_set

    def parse(self):
        src_ds = gdal.Open(self.fn)
        if src_ds is not None:
            sub_datasets = src_ds.GetSubDatasets()

            idx = 2
            if utils.int_or(self.data_set) is not None:
                idx = utils.int_or(self.data_set)
            else:
                for j, sd in enumerate(sub_datasets):
                    _name = sd[0].split(':')[-1]
                    if self.data_set == _name:
                        idx = j
                        break

            src_ds = None
            src_srs = gdalfun.gdal_get_srs(sub_datasets[idx][0])
            if self.data_set == 'wse':
                src_srs = gdalfun.combine_epsgs(
                    src_srs, '3855', name='SWOT Combined'
                )
            sub_ds = DatasetFactory(
                **self._set_params(
                    mod=sub_datasets[idx][0],
                    data_format=200,
                    node='grid',
                    check_path=False
                )
            )._acquire_module()
            self.data_entries.append(sub_ds)
            sub_ds.initialize()
            for gdal_ds in sub_ds.parse():
                yield(gdal_ds)

class IceSat2File(ElevationDataset):
    """IceSat2 data from NASA

    Parameters:
    
    water_surface: 'geoid' # this is the vertical datum, can be 'geoid', 'ellipsoid' or 'mean_tide'
    classes: None # return only data with the specified classes, e.g. '2/3/4'
    confidence_levels: None # return only data with the specified confidence levels, e.g. '2/3/4'
    columns: {} # the additional columns to export in transform_and_yield_points {'/h5/atl/path', column_name}
    classify_bathymetry: True # extract bathymetry with CShelph
    classify_buildings: True # classify buildings using BING BFP
    classify_water: True # classify water using OSM
    reject_failed_qa: True # skip granules that failed QA

    Classes:
    -1 - no classification (ATL08)
    0 - noise / atmosphere (ATL08)
    1 - ground surface (ATL08)
    2 - canopy (ATL08)
    3 - canopy top (ATL08)
    40 - bathymetry floor surface (CShelph, ATL24)
    41 - bathymetry water surface (OSM coastline, ATL24)
    6 - ice surface (ATL06) (unused for now, just planning ahead for possible future ATL06 integration)
    7 - built structure (OSM or Bing)
    8 - "urban" (WSF, if used)

    Confidence Levels:
    0, 1, 2, 3, 4
    """
    
    def __init__(self,
                 water_surface = 'geoid',
                 classes = None,
                 confidence_levels = '2/3/4',
                 columns = {},
                 classify_bathymetry = True,
                 classify_buildings = True,
                 classify_inland_water = True,
                 reject_failed_qa = True,
                 classify_water = True,
                 **kwargs):
        super().__init__(**kwargs)
        self.data_format = 303
        self.water_surface = water_surface
        if self.water_surface not in ['mean_tide', 'geoid', 'ellipsoid']:
            self.water_surface = 'mean_tide'

        self.classes = [int(x) for x in classes.split('/')] if classes is not None else []
        self.confidence_levels = [int(x) for x in confidence_levels.split('/')] if confidence_levels is not None else []
        self.columns = columns
        self.atl_fn = None
        self.want_bathymetry = classify_bathymetry
        self.want_buildings = classify_buildings
        self.want_watermask = classify_water
        self.want_inland_water = classify_inland_water
        self.reject_failed_qa = reject_failed_qa
        
    def init_atl_h5(self):
        """initialize all the relevant ATL h5 files"""

        self.atl03_f = None
        self.atl08_f = None
        self.atl24_f = None
        self.atl13_f = None

        if self.atl03_fn is not None and os.path.exists(self.atl03_fn):
            self.atl03_f = h5.File(self.atl03_fn, 'r')            
            if 'short_name' not in self.atl03_f.attrs.keys():
                raise UnboundLocalError('this file does not appear to be an ATL03 file')
        
            if self.reject_failed_qa:
                if self.atl03_f['/quality_assessment/qa_granule_pass_fail'][...,][0] != 0:
                    raise UnboundLocalError(
                        'this granule has failed qa {}'.format(
                            self.atl03_f['/quality_assessment/qa_granule_pass_fail'][...,][0]
                        )
                    )
        
        if self.atl08_fn is not None and os.path.exists(self.atl08_fn):
            self.atl08_f = h5.File(self.atl08_fn, 'r')
            if 'short_name' not in self.atl08_f.attrs.keys():
                utils.echo_warning_msg('this file does not appear to be an ATL file')
                self.atl08_f.close()

        if self.atl24_fn is not None and os.path.exists(self.atl24_fn):
            self.atl24_f = h5.File(self.atl24_fn, 'r')
            if 'short_name' not in self.atl24_f.attrs.keys():
                utils.echo_warning_msg('this file does not appear to be an ATL file')
                self.atl24_f.close()
                
        if self.atl13_fn is not None and os.path.exists(self.atl13_fn):
            self.atl13_f = h5.File(self.atl13_fn, 'r')
            if 'short_name' not in self.atl13_f.attrs.keys():
                utils.echo_warning_msg('this file does not appear to be an ATL file')
                self.atl13_f.close()                

    def close_atl_h5(self):
        """close all open atl files"""
        
        self.atl03_fn = self.atl08_fn = self.atl24_fn = None
        if self.atl03_f is not None:
            self.atl03_f.close()
            
        if self.atl08_f is not None:
            self.atl08_f.close()

        if self.atl24_f is not None:
            self.atl24_f.close()

        if self.atl13_f is not None:
            self.atl13_f.close()
            
    def yield_points(self):
        """yield the points from the dataset.

        In this case, it will yield a pandas dataframe
        """

        dataset = None
        self.atl03_fn = self.fn
        self.atl08_fn = None
        self.atl24_fn = None
        self.atl13_fn = None
        if len(self.classes) > 0:
            atl08_result = self.fetch_atlxx(short_name='ATL08')
            self.atl08_fn = atl08_result
            atl13_result = self.fetch_atlxx(short_name='ATL13')
            self.atl13_fn = atl13_result
            atl24_result = self.fetch_atlxx(short_name='ATL24')
            self.atl24_fn = atl24_result

        try:
            self.init_atl_h5()
        except Exception as e:
            utils.echo_error_msg('could not initialize data {}'.format(e))
            self.close_atl_h5()
            return

        ## fetch and process buildings, if wanted
        this_bing = None
        if self.want_buildings:
            if isinstance(self.want_buildings, bool):
                this_bing = self.process_buildings(
                    self.fetch_buildings(verbose=False),
                    verbose=False
                )
            elif isinstance(self.want_buildings, list):
                this_bing = self.want_buildings         

        ## fetch and process watermask, if wanted
        this_wm = None
        if self.want_watermask:
            if isinstance(self.want_watermask, bool):
                this_wm = self.process_coastline(
                    self.fetch_coastline(chunks=False, verbose=False),
                    return_geom=True,
                    verbose=False
                )
            elif isinstance(self.want_watermask, list):
                this_wm = self.want_watermask                

        ## parse through the icesat2 file by laser number
        for i in range(1, 4):
            #try:
            dataset = self.read_atl_data('{}'.format(i))
            if dataset is None or len(dataset) == 0:
                continue

            ## keep only photons with confidence levels mentioned in `self.confidence_levels`
            if len(self.confidence_levels) > 0:
                dataset = dataset[(np.isin(dataset['confidence'], self.confidence_levels))]

            if dataset is None or len(dataset) == 0:
                continue

            ## re-classify photons based on buildings/watermask/bathymetry
            if self.want_buildings and this_bing is not None:
                dataset = self.classify_buildings(dataset, this_bing)

            if self.want_watermask and this_wm is not None:
                dataset = self.classify_water(dataset, this_wm)

            # if self.want_bathymetry:
            #     dataset = self.classify_bathymetry(dataset)
                
            if dataset is None or len(dataset) == 0:
                continue

            ## keep only photons with a classification mentioned in `self.classes`
            if len(self.classes) > 0:
                dataset = dataset[(np.isin(dataset['ph_h_classed'], self.classes))]
                
            if dataset is None or len(dataset) == 0:
                continue

            ## rename the x,y,z columns for `transform_and_yield_points`
            dataset.rename(
                columns={'longitude': 'x', 'latitude': 'y', 'photon_height': 'z'},
                inplace=True
            )
            yield(dataset)

        self.close_atl_h5()

    # def fetch_atl08(self):
    #     """fetch the associated ATL08 file"""
        
    #     atl08_filter = utils.fn_basename2(self.atl03_fn).split('ATL03_')[1]
    #     this_atl08 = fetches.IceSat2(
    #         src_region=None, verbose=self.verbose, outdir=self.cache_dir, short_name='ATL08',
    #         filename_filter=atl08_filter
    #     )
    #     this_atl08.run()
    #     if len(this_atl08.results) == 0:
    #         utils.echo_warning_msg('could not locate associated atl08 file for {}'.format(atl08_filter))
    #         return(None)
    #     else:
    #         if this_atl08.fetch(this_atl08.results[0], check_size=True) == 0:
    #             return(os.path.join(this_atl08._outdir, this_atl08.results[0][1]))

    # def fetch_atl24(self):
    #     """fetch the associated ATL24 file"""
        
    #     atl24_filter = utils.fn_basename2(self.atl03_fn).split('ATL03_')[1]
    #     this_atl24 = fetches.IceSat2(
    #         src_region=None, verbose=self.verbose, outdir=self.cache_dir, short_name='ATL24',
    #         filename_filter=atl24_filter
    #     )
    #     this_atl24.run()
    #     if len(this_atl24.results) == 0:
    #         utils.echo_warning_msg('could not locate associated atl24 file for {}'.format(atl24_filter))
    #         return(None)
    #     else:
    #         if this_atl24.fetch(this_atl24.results[0], check_size=True) == 0:
    #             return(os.path.join(this_atl24._outdir, this_atl24.results[0][1]))            

    def fetch_atlxx(self, short_name='ATL08'):
        """fetch an associated ATLxx file"""

        #utils.echo_msg(self.atl03_fn)
        atlxx_filter = utils.fn_basename2(self.atl03_fn).split('ATL03_')[1]
        this_atlxx = fetches.IceSat2(
            src_region=None, verbose=self.verbose, outdir=self.cache_dir, short_name=short_name,
            filename_filter=atlxx_filter
        )
        this_atlxx.run()
        if len(this_atlxx.results) == 0:
            atlxx_filter = '_'.join(utils.fn_basename2(self.atl03_fn).split('ATL03_')[1].split('_')[:-1])
            this_atlxx = fetches.IceSat2(
                src_region=None, verbose=self.verbose, outdir=self.cache_dir, short_name=short_name,
                filename_filter=atlxx_filter
            )
            this_atlxx.run()
            if len(this_atlxx.results) == 0:
                utils.echo_warning_msg('could not locate associated {} file for {}'.format(short_name, atlxx_filter))
                return(None)
        else:
            if this_atlxx.fetch(this_atlxx.results[0], check_size=True) == 0:
                return(os.path.join(this_atlxx._outdir, this_atlxx.results[0][1]))            
            
    def read_atl_data(self, laser_num, orientation = None):
        """Read data from an ATL03 file

        Adapted from 'cshelph' https://github.com/nmt28/C-SHELPh.git 
        and 'iVert' https://github.com/ciresdem/ivert.git

        laser_num is 1, 2 or 3
        surface is 'mean_tide', 'geoid' or 'ellipsoid'
        """

        if orientation is None:
            orientation = self.atl03_f['/orbit_info/sc_orient'][0]
            
        ## selects the strong beams only [we can include weak beams later on]
        orientDict = {0:'l', 1:'r', 21:'error'}
        laser = 'gt' + laser_num + orientDict[orientation]

        ## for 'subsets', where heights don't come through
        if 'heights' not in self.atl03_f['/{}'.format(laser)].keys():
            return(None)
        
        ## Read in the required photon level data
        photon_h = self.atl03_f['/' + laser + '/heights/h_ph'][...,]
        latitude = self.atl03_f['/' + laser + '/heights/lat_ph'][...,]
        longitude = self.atl03_f['/' + laser + '/heights/lon_ph'][...,]
        ph_count = self.atl03_f['/' + laser + '/heights/ph_id_count'][...,]
        conf = self.atl03_f['/' + laser + '/heights/signal_conf_ph/'][...,0]
        qual = self.atl03_f['/' + laser + '/heights/quality_ph/'][...,0]
        dist_ph_along = self.atl03_f['/' + laser + '/heights/dist_ph_along'][...,]
        this_N = latitude.shape[0]

        ## Read in the geolocation level data
        segment_ph_cnt = self.atl03_f['/' + laser + '/geolocation/segment_ph_cnt'][...,]
        segment_id = self.atl03_f['/' + laser + '/geolocation/segment_id'][...,]
        segment_dist_x = self.atl03_f['/' + laser + '/geolocation/segment_dist_x'][...,]
        ref_elev = self.atl03_f['/' + laser + '/geolocation/ref_elev'][...,]
        ref_azimuth = self.atl03_f['/' + laser + '/geolocation/ref_azimuth'][...,]
        ref_photon_index = self.atl03_f['/' + laser + '/geolocation/reference_photon_index'][...,]
        ph_index_beg = self.atl03_f['/' + laser + '/geolocation/ph_index_beg'][...,]
        altitude_sc = self.atl03_f['/' + laser + '/geolocation/altitude_sc'][...,]

        ## Read in the geoid data
        photon_geoid = self.atl03_f['/' + laser + '/geophys_corr/geoid'][...,]
        photon_geoid_f2m = self.atl03_f['/' + laser + '/geophys_corr/geoid_free2mean'][...,]

        ## Create a dictionary with (segment_id --> index into ATL03 photons)
        ## lookup pairs, for the starting photon of each segment
        #segment_indices = ph_index_beg+(ref_photon_index-1)
        #segment_index_dict = dict(zip(segment_id, ph_index_beg+ref_photon_index-1))
        segment_indices = np.concatenate(([0], np.cumsum(segment_ph_cnt)[:-1]))
        segment_index_dict = dict(zip(segment_id, segment_indices))
        ph_segment_ids = segment_id[np.searchsorted(segment_indices, np.arange(0.5, len(photon_h), 1))-1]
        
        ## Compute the total along-track distances.
        #segment_dist_dict = dict(zip(segment_id, segment_dist_x))

        ## Determine where in the array each segment index needs to look.
        #ph_segment_dist_x = np.array(list(map((lambda pid: segment_dist_dict[pid]), ph_segment_ids)))
        #dist_x = ph_segment_dist_x + dist_ph_along

        ## meantide/geoid heights
        h_geoid_dict = dict(zip(segment_id, photon_geoid))
        ph_h_geoid = np.array(list(map((lambda pid: h_geoid_dict[pid]), ph_segment_ids)))        
        h_meantide_dict = dict(zip(segment_id, photon_geoid_f2m))
        ph_h_meantide = np.array(list(map((lambda pid: h_meantide_dict[pid]), ph_segment_ids)))
        photon_h_geoid = photon_h - ph_h_geoid
        photon_h_meantide = photon_h - (ph_h_geoid + ph_h_meantide)

        ## setup classifications
        ph_h_classed = np.zeros(photon_h.shape)
        ph_h_classed[:] = -1

        ## append the laser to each record
        laser_arr = np.empty(photon_h.shape, dtype='object')
        laser_arr[:] = laser

        ## append the filename to each record
        fn_arr = np.empty(photon_h.shape, dtype='object')
        fn_arr[:] = self.fn
        
        ## ref values
        h_ref_elev_dict = dict(zip(segment_id, ref_elev))
        ph_ref_elev = np.array(list(map((lambda pid: h_ref_elev_dict[pid]), ph_segment_ids)))#.astype(float)        
        h_ref_azimuth_dict = dict(zip(segment_id, ref_azimuth))
        ph_ref_azimuth = np.array(list(map((lambda pid: h_ref_azimuth_dict[pid]), ph_segment_ids)))#.astype(float)
        h_altitude_sc_dict = dict(zip(segment_id, altitude_sc))
        ph_altitude_sc = np.array(list(map((lambda pid: h_altitude_sc_dict[pid]), ph_segment_ids)))#.astype(float)

        ## Read in the ATL08 data
        if self.atl08_f is not None:
            ## classed flag (signal_photons)
            atl08_classed_pc_flag  = self.atl08_f['/' + laser + '/signal_photons/classed_pc_flag'][...,]
            atl08_ph_segment_id = self.atl08_f['/' + laser + '/signal_photons/ph_segment_id'][...,] # photon src 20 m seg id
            atl08_classed_pc_indx = self.atl08_f['/' + laser + '/signal_photons/classed_pc_indx'][...,]

            ## set the classifications from ATL08
            atl08_segment_id_msk = [True if x in segment_id else False for x in atl08_ph_segment_id]
            atl08_ph_segment_indx = np.array(list(map((lambda pid: segment_index_dict[pid]), atl08_ph_segment_id[atl08_segment_id_msk])))
            atl08_ph_index = np.array(atl08_ph_segment_indx + (atl08_classed_pc_indx[atl08_segment_id_msk] - 1), dtype=int)
            class_mask = atl08_ph_index < len(ph_segment_ids)
            ph_h_classed[atl08_ph_index[class_mask]] = atl08_classed_pc_flag[atl08_segment_id_msk][class_mask]

            ## old from IVERT
            # dict_success = False
            # while not dict_success:
            #     try:
            #         atl_08_ph_segment_indx = np.array(list(map((lambda pid: segment_index_dict[pid]), atl_08_ph_segment_id)))
            #     except KeyError as e:
            #         # One of the atl08_ph_segment_id entries does not exist in the atl03 granule, which
            #         # causes problems here. Eliminate it from the list and try again.
            #         problematic_id = e.args[0]
            #         good_atl08_mask = (atl_08_ph_segment_id != problematic_id)
            #         atl_08_classed_pc_flag = atl_08_classed_pc_flag[good_atl08_mask]
            #         atl_08_ph_segment_id = atl_08_ph_segment_id[good_atl08_mask]
            #         atl_08_classed_pc_indx = atl_08_classed_pc_indx[good_atl08_mask]
            #         # Then, try the loop again.
            #         continue
            #     dict_success = True
            # atl_08_ph_index = np.array(atl_08_ph_segment_indx + atl_08_classed_pc_indx - 1 , dtype=int)
            # ph_h_classed[atl_08_ph_index] = atl_08_classed_pc_flag[atl_08_segment_id_msk]            

        ## Read in the ATL24 data
        ## todo: check re subsets
        if self.atl24_f is not None:
            atl24_classed_pc_flag  = self.atl24_f['/' + laser + '/class_ph'][...,]
            atl24_classed_pc_indx = self.atl24_f['/' + laser + '/index_ph'][...,]
            atl24_longitude = self.atl24_f['/' + laser + '/lon_ph'][...,]
            atl24_latitude = self.atl24_f['/' + laser + '/lat_ph'][...,]
            atl24_surface_h = self.atl24_f['/' + laser + '/surface_h'][...,]
            atl24_ellipse_h = self.atl24_f['/' + laser + '/ellipse_h'][...,]
            atl24_ortho_h = self.atl24_f['/' + laser + '/ortho_h'][...,]

            class_40_mask = atl24_classed_pc_flag == 40
            ph_h_classed[atl24_classed_pc_indx] = atl24_classed_pc_flag

            # we also need to change the lon/lat/height values to the updated bathymetry values (we'll just do it to class 40)
            longitude[atl24_classed_pc_indx[class_40_mask]] = atl24_longitude[class_40_mask]
            latitude[atl24_classed_pc_indx[class_40_mask]] = atl24_latitude[class_40_mask]
            photon_h[atl24_classed_pc_indx[class_40_mask]] = atl24_ellipse_h[class_40_mask]
            photon_h_geoid[atl24_classed_pc_indx[class_40_mask]] = atl24_ortho_h[class_40_mask]

        if self.atl13_f is not None:
            atl13_refid = self.atl13_f['/' + laser + '/segment_id_beg'][...,]
            ph_h_classed[atl13_refid] = 44
            
        ## set the photon height, either 'mean_tide' or 'geoid', else ellipsoid
        if self.water_surface == 'mean_tide':
            ph_height = photon_h_meantide
        elif self.water_surface == 'geoid':
            ph_height = photon_h_geoid
        else:
            ph_height = photon_h
            
        ## create the pandas dataframe            
        dataset = pd.DataFrame(
            {'latitude': latitude,
             'longitude': longitude,
             'photon_height': ph_height,
             'laser': laser_arr,
             'fn': fn_arr,
             'confidence': conf,
             'ref_elevation':ph_ref_elev,
             'ref_azimuth':ph_ref_azimuth,
             'ref_sat_alt':ph_altitude_sc,
             'ph_h_classed': ph_h_classed},
            columns=['latitude', 'longitude', 'photon_height', 'laser', 'fn',
                     'confidence', 'ref_elevation', 'ref_azimuth', 'ref_sat_alt',
                     'ph_h_classed']
        )

        ## Process extra columns specified in `self.columns`
        for col in self.columns.keys():
            try:
                if 'gtx' in col:
                    _col = col.replace('/gtx', '/' + laser)
                    col_arr = self.atl03_f[_col][...,]                
                    if 'heights' not in _col: 
                        col_dict = dict(zip(segment_id, col_arr))
                        col_arr = np.array(list(map((lambda pid: col_dict[pid]), ph_segment_ids)))
                else:
                    col_arr = np.empty(photon_h.shape, dtype='object')
                    col_arr[:] = self.atl03_f[col][...,]

                extra_dataset = pd.DataFrame(
                    {self.columns[col]: col_arr},
                    columns = [self.columns[col]]
                )
                dataset = dataset.join(extra_dataset)
            except:
                utils.echo_warning_msg('could not find and/or process {}'.format(col))

        return(dataset)        
            
    def classify_bathymetry(
            self, dataset, thresh = 95, min_buffer = -40, max_buffer = 5,
            start_lat = False, end_lat = False, lat_res = 10, h_res = .5,
            surface_buffer = -.5, water_temp = None
    ):
        """Classify bathymetry in an ATL03 file. 

        This uses C-Shelph to locate, extract and process bathymetric photons.
        This function is adapted from the C-Shelph CLI
        Depreciated once ATL24 are online
        """
        
        water_temp = utils.float_or(water_temp)
        epsg_code = cshelph.convert_wgs_to_utm(dataset.latitude.iloc[0], dataset.longitude.iloc[0])
        epsg_num = int(epsg_code.split(':')[-1])
        utm_proj = pyproj.Proj(epsg_code)
        lon_utm, lat_utm = utm_proj(dataset.longitude, dataset.latitude)
                       
        ## Aggregate data into dataframe
        dataset_sea = pd.DataFrame(
            {'latitude': lat_utm,
             'longitude': lon_utm,
             'photon_height': dataset.photon_height,
             'laser': dataset.laser,
             'fn': dataset.fn,
             'confidence': dataset.confidence,
             'ref_elevation': dataset.ref_elevation,
             'ref_azimuth': dataset.ref_azimuth,
             'ref_sat_alt': dataset.ref_sat_alt,
             'ph_h_classed': dataset.ph_h_classed},
            columns=['latitude', 'longitude', 'photon_height',
                     'laser', 'fn', 'confidence', 'ref_elevation',
                     'ref_azimuth', 'ref_sat_alt', 'ph_h_classed']
        )
        dataset_sea1 = dataset_sea[(dataset_sea['photon_height'] > min_buffer) & (dataset_sea['photon_height'] < max_buffer)]
        dataset_sea1 = dataset_sea1[(dataset_sea1['ph_h_classed'] == 5)]
        
        if self.region is not None:
            xyz_region = self.region.copy()
            xyz_region.epsg = 'epsg:4326'
            xyz_region.warp('epsg:{}'.format(epsg_num))
            dataset_sea1 = dataset_sea1[(dataset_sea1['latitude'] > xyz_region.ymin) & (dataset_sea1['latitude'] < xyz_region.ymax)]
            dataset_sea1 = dataset_sea1[(dataset_sea1['longitude'] > xyz_region.xmin) & (dataset_sea1['longitude'] < xyz_region.xmax)]

        if len(dataset_sea1) > 0:
            binned_data_sea = cshelph.bin_data(dataset_sea1, lat_res, h_res)
            if binned_data_sea is not None:
                if water_temp is None:
                    try:
                        #water_temp = cshelph.get_water_temp(self.fn, latitude, longitude)
                        ## water_temp via fetches instead of earthaccess
                        water_temp = self.get_water_temp()
                    except Exception as e:
                        #utils.echo_warning_msg('NO SST PROVIDED OR RETRIEVED: 20 degrees C assigned')
                        water_temp = 20

                #utils.echo_msg('water temp is {}'.format(water_temp))
                sea_height = cshelph.get_sea_height(binned_data_sea, surface_buffer)
                sea_height1 = cshelph.get_bin_height(binned_data_sea, 60, surface_buffer)
                
                if sea_height is not None:
                    med_water_surface_h = np.nanmedian(sea_height) #* -1
                    
                    ## Correct for refraction
                    ref_x, ref_y, ref_z, ref_conf, raw_x, raw_y, raw_z, ph_ref_azi, ph_ref_elev = cshelph.refraction_correction(
                        water_temp, med_water_surface_h, 532, dataset_sea1.ref_elevation, dataset_sea1.ref_azimuth,
                        dataset_sea1.photon_height, dataset_sea1.longitude, dataset_sea1.latitude, dataset_sea1.confidence,
                        dataset_sea1.ref_sat_alt
                    )
                    depth = med_water_surface_h - ref_z

                    # Create new dataframe with refraction corrected data
                    dataset_bath = pd.DataFrame({'latitude': raw_y, 'longitude': raw_x, 'cor_latitude':ref_y, 'cor_longitude':ref_x,
                                                 'cor_photon_height':ref_z, 'photon_height': raw_z, 'confidence':ref_conf, 'depth':depth},
                                                columns=['latitude', 'longitude', 'photon_height', 'cor_latitude','cor_longitude',
                                                         'cor_photon_height', 'confidence', 'depth'])

                    # Bin dataset again for bathymetry
                    if len(dataset_bath) > 0:
                        binned_data = cshelph.bin_data(dataset_bath, lat_res, h_res)
                        if binned_data is not None:
                            # Find bathymetry
                            bath_height, geo_df = cshelph.get_bath_height(binned_data, thresh, med_water_surface_h, h_res)
                            if bath_height is not None:
                                transformer = pyproj.Transformer.from_crs("EPSG:"+str(epsg_num), "EPSG:4326", always_xy=True)
                                lon_wgs84, lat_wgs84 = transformer.transform(geo_df.longitude.values, geo_df.latitude.values)
                                bath_height = [x for x in bath_height if ~np.isnan(x)]
                                
                                for n, id in enumerate(geo_df.ids.values):
                                    dataset.at[id, 'ph_h_classed'] = 4
                                    #dataset.at[id, 'latitude'] = lat_wgs84[n]
                                    #dataset.at[id, 'longitude'] = lon_wgs84[n]
                                    dataset.at[id, 'photon_height'] = geo_df.depth.values[n] * -1
                            
                                return(dataset)
            
        return(dataset)

    def _vectorize_df(self, dataset):
        """Make a point vector OGR DataSet Object from a pandas dataframe

        This is to allow spatial filtering for watermask and buildings.
        """

        dst_ogr = '{}'.format('icesat_dataframe')
        ogr_ds = gdal.GetDriverByName('Memory').Create(
            '', 0, 0, 0, gdal.GDT_Unknown
        )
        layer = ogr_ds.CreateLayer(
            dst_ogr, geom_type=ogr.wkbPoint
        )
        fd = ogr.FieldDefn('index', ogr.OFTInteger)
        layer.CreateField(fd)
        # for x in dataset.columns.to_list():
        #     fd = ogr.FieldDefn(x, ogr.OFTReal)
        #     fd.SetWidth(12)
        #     fd.SetPrecision(8)
        #     layer.CreateField(fd)
            
        f = ogr.Feature(feature_def=layer.GetLayerDefn())        
        with tqdm(desc='vectorizing dataframe', leave=False) as pbar:
            for index, this_row in dataset.iterrows():
                pbar.update()
                f.SetField(0, index)
                #[f.SetField(n+1, this_row[x]) for n,x in enumerate(dataset.columns.to_list())]
                g = ogr.CreateGeometryFromWkt('POINT ({} {})'.format(this_row['longitude'], this_row['latitude']))
                f.SetGeometryDirectly(g)
                layer.CreateFeature(f)
            
        return(ogr_ds)
    
    def classify_buildings(self, dataset, these_bldgs):
        """classify building photons using BING building footprints 
        """

        ## vectorize the icesate2 photons
        ogr_df = self._vectorize_df(dataset)
        
        os.environ["OGR_OSM_OPTIONS"] = "INTERLEAVED_READING=YES"
        os.environ["OGR_OSM_OPTIONS"] = "OGR_INTERLEAVED_READING=YES"
        with tqdm(
                total=len(these_bldgs),
                desc='classifying building photons',
                leave=False
        ) as pbar:
            for n, bldg_geom in enumerate(these_bldgs):                
                pbar.update()
                icesat_layer = ogr_df.GetLayer()
                icesat_layer.SetSpatialFilter(bldg_geom)
                #[dataset.at[f.GetField('index'), 'ph_h_classed'] = 7 for f in icesat_layer]
                for f in icesat_layer:
                    idx = f.GetField('index')
                    dataset.at[idx, 'ph_h_classed'] = 7

                icesat_layer.SetSpatialFilter(None)

        ogr_df = None
        return(dataset)

    def classify_water(self, dataset, these_wms):
        """classify water photons using OSM coastline 
        """

        ## vectorize the photons
        ogr_df = self._vectorize_df(dataset)
        
        os.environ["OGR_OSM_OPTIONS"] = "INTERLEAVED_READING=YES"
        os.environ["OGR_OSM_OPTIONS"] = "OGR_INTERLEAVED_READING=YES"
        with tqdm(
                total=len(these_wms),
                desc='classifying water photons',
                leave=False
        ) as pbar:
            for n, wm_geom in enumerate(these_wms):                
                pbar.update()
                icesat_layer = ogr_df.GetLayer()
                icesat_layer.SetSpatialFilter(wm_geom)
                #[dataset.at[f.GetField('index'), 'ph_h_classed'] = 7 for f in icesat_layer]
                for f in icesat_layer:
                    idx = f.GetField('index')
                    dataset.at[idx, 'ph_h_classed'] = 5

                icesat_layer.SetSpatialFilter(None)

        ogr_df = None
        return(dataset)
    
class MBSParser(ElevationDataset):
    """providing an mbsystem parser

    Process MB-System supported multibeam data files.
    See MB-System for more information regarding supported
    file formats, etc.
    
    generate_inf - generate an inf file for the MBS data
    yield_xyz - yield the MBS data as xyz
    yield_array - yield the MBS data as an array
    
    -----------
    Parameters:
    
    mb_fmt=[]
    mb_exclude=[]
    """

    def __init__(self,
                 mb_fmt = None,
                 mb_exclude = 'A',
                 want_mbgrid = False,
                 want_binned = False,
                 min_year = None,
                 auto_weight = True,
                 auto_uncertainty = True,
                 **kwargs):
        super().__init__(**kwargs)
        self.mb_fmt = mb_fmt
        self.mb_exclude = mb_exclude
        self.want_mbgrid = want_mbgrid
        self.want_binned = want_binned
        self.min_year = min_year
        self.auto_weight = auto_weight
        self.auto_uncertainty = auto_uncertainty

    def inf_parse(self):
        self.infos.minmax = [0,0,0,0,0,0]
        this_row = 0
        xinc = 0
        yinc = 0
        dims = []
        inf_fn = '{}.inf'.format(self.fn) if self.fn.split('.')[-1] != 'inf' else self.fn
        with open(inf_fn) as iob:
            for il in iob:
                til = il.split()
                if len(til) > 1:

                    if til[0] == 'Swath':
                        if til[2] == 'File:':
                            self.infos.name = til[3]

                    if ' '.join(til[:-1]) == 'MBIO Data Format ID:':
                        self.mb_fmt = til[-1]
                            
                    if til[0] == 'Number':
                        if til[2] == 'Records:':
                            self.infos.numpts = utils.int_or(til[3])
                            
                    if til[0] == 'Minimum':
                        if til[1] == 'Longitude:':
                            self.infos.minmax[0] = utils.float_or(til[2])
                            self.infos.minmax[1] = utils.float_or(til[5])
                        elif til[1] == 'Latitude:':
                            self.infos.minmax[2] = utils.float_or(til[2])
                            self.infos.minmax[3] = utils.float_or(til[5])
                        elif til[1] == 'Depth:':
                            self.infos.minmax[4] = utils.float_or(til[5]) * -1
                            self.infos.minmax[5] = utils.float_or(til[2]) * -1
                            
                    if til[0] == 'CM':
                        if til[1] == 'dimensions:':
                            dims = [utils.int_or(til[2]), utils.int_or(til[3])]
                            cm_array = np.zeros((dims[0], dims[1]))
                            
                    if til[0] == 'CM:':
                        for j in range(0, dims[0]):
                            cm_array[this_row][j] = utils.int_or(til[j+1])
                        this_row += 1

        mbs_region = regions.Region().from_list(self.infos.minmax)
        if len(dims) > 0:
            xinc = (mbs_region.xmax - mbs_region.xmin) / dims[0]
            yinc = (mbs_region.ymin - mbs_region.ymax) / dims[1]
        
        if abs(xinc) > 0 and abs(yinc) > 0:
            xcount, ycount, dst_gt = mbs_region.geo_transform(
                x_inc=xinc, y_inc=yinc
            )
            ds_config = {'nx': dims[0], 'ny': dims[1], 'nb': dims[1] * dims[0],
                         'geoT': dst_gt, 'proj': gdalfun.osr_wkt(self.src_srs),
                         'dt': gdal.GDT_Float32, 'ndv': 0, 'fmt': 'GTiff'}
            driver = gdal.GetDriverByName('MEM')
            ds = driver.Create(
                'tmp', ds_config['nx'], ds_config['ny'], 1, ds_config['dt']
            )
            ds.SetGeoTransform(ds_config['geoT'])
            if ds_config['proj'] is not None:
                ds.SetProjection(ds_config['proj'])
                
            ds_band = ds.GetRasterBand(1)
            ds_band.SetNoDataValue(ds_config['ndv'])
            ds_band.WriteArray(cm_array)
            tmp_ds = ogr.GetDriverByName('Memory').CreateDataSource('tmp_poly')
            tmp_layer = tmp_ds.CreateLayer('tmp_poly', None, ogr.wkbMultiPolygon)
            tmp_layer.CreateField(ogr.FieldDefn('DN', ogr.OFTInteger))
            gdal.Polygonize(ds_band, ds_band, tmp_layer, 0)
            multi = ogr.Geometry(ogr.wkbMultiPolygon)
            for feat in tmp_layer:
                feat.geometry().CloseRings()
                wkt = feat.geometry().ExportToWkt()
                multi.AddGeometryDirectly(ogr.CreateGeometryFromWkt(wkt))
                
            wkt = multi.ExportToWkt()
            tmp_ds = ds = None
        else:
            wkt = mbs_region.export_as_wkt()

        self.infos.wkt = wkt
        return(self)

    def parse_(self):
        """use mblist to convert data to xyz then set the dataset as xyz and 
        use that to process...
        """
        
        mb_fn = os.path.join(self.fn)
        if self.region is not None:
            mb_region = self.region.copy()
            mb_region.buffer(pct=25)
        else:
            mb_region = None

        out_mb = utils.make_temp_fn('{}_.xyz'.format(utils.fn_basename2(mb_fn)), self.cache_dir)
        out, status = utils.run_cmd('mblist -M{}{} -OXYZ -I{} > {}'.format(
            self.mb_exclude, ' {}'.format(
                mb_region.format('gmt') if mb_region is not None else ''
            ), mb_fn, out_mb
        ), verbose=True)

        if status == 0:
            data_set = DatasetFactory(
                **self._set_params(
                    mod=out_mb, data_foramt=168
                )
            )._acquire_module.initialize()
            for ds in data_set.parse(): # fill self.data_entries with each dataset for use outside the yield.
                self.data_entries.append(ds) 
                yield(ds)

        utils.remove_glob('{}*'.format(out_mb))

    def mb_inf_format(self, src_inf):
        """extract the format from the mbsystem inf file."""

        with open(src_inf, errors='ignore') as iob:
            for il in iob:
                til = il.split()
                if len(til) > 1:
                    #utils.echo_msg(til[0].strip())
                    if til[0].strip() == 'MBIO':
                        return(til[4])
        return(None)
        
    def mb_inf_data_date(self, src_inf):
        """extract the date from the mbsystem inf file."""

        with open(src_inf, errors='ignore') as iob:
            for il in iob:
                til = il.split()
                if len(til) > 1:
                    if til[0] == 'Time:':
                        return(til[3])
        return(None)

    def mb_inf_perc_good(self, src_inf):
        """extract the data format from the mbsystem inf file."""

        with open(src_inf, errors='ignore') as iob:
            for il in iob:
                til = il.split(':')
                if len(til) > 1:
                    if til[0].strip() == 'Number of Good Beams':
                        return(til[1].split()[-1].split('%')[0])
        return(None)
        
    def yield_mbgrid_ds(self):
        """process the data through mbgrid and use GDALFile to further process the gridded data"""
        
        mb_fn = os.path.join(self.fn)
        with open('_mb_grid_tmp.datalist', 'w') as tmp_dl:
            tmp_dl.write(
                '{} {} {}\n'.format(
                    self.fn,
                    self.mb_fmt if self.mb_fmt is not None else '',
                    self.weight if self.mb_fmt is not None else ''
                )
            )

        ofn = '_'.join(os.path.basename(mb_fn).split('.')[:-1])
        try:
            utils.run_cmd(
                'mbgrid -I_mb_grid_tmp.datalist {} -E{}/{}/degrees! -O{} -A2 -F1 -C2/1 -S0 -T35 -W.1'.format(
                    self.data_region.format('gmt'), self.x_inc, self.y_inc, ofn
                ), verbose=True
            )

            gdalfun.gdal2gdal('{}.grd'.format(ofn))
            utils.remove_glob(
                '_mb_grid_tmp.datalist',
                '{}.cmd'.format(ofn),
                '{}.mb-1'.format(ofn),
                '{}.grd*'.format(ofn)
            )
            
            grits_filter = grits.GritsFactory(
                mod='outliers:multipass=2', src_dem='{}.tif'.format(ofn))._acquire_module()
            
            if grits_filter is not None:
                try:
                    grits_filter()
                    os.replace(grits_filter.dst_dem, '{}'.format(ofn))
                except:
                    pass

            mbs_ds = DatasetFactory(
                **self._set_params(mod='{}.tif'.format(ofn), data_format=200)
            )._acquire_module()
            mbs_ds.initialize()
            for gdal_ds in mbs_ds.parse():
                for pts in gdal_ds.transform_and_yield_points():
                    yield(pts)
            
            #yield(mbs_ds)
            utils.remove_glob('{}*.tif*'.format(ofn))
        except Exception as e:
            #pass
            raise(e)

    def yield_mblist_ds(self):
        """use mblist to process the multibeam data"""
        
        mb_fn = os.path.join(self.fn)
        xs = []
        ys = []
        zs = []
        ws = []
        us = []
        if self.region is not None:
            mb_region = self.region.copy()
            mb_region.buffer(pct=25)
        else:
            mb_region = None

        src_inf = '{}.inf'.format(self.fn)
        try:
            #mb_date = self.mb_inf_data_date(src_inf)
            #mb_perc = self.mb_inf_perc_good(src_inf)
            mb_format = self.mb_inf_format(src_inf)
            if mb_fn.split('.')[-1] == 'fbt':
                mb_format = None
        except:
            mb_format = None

        #'mblist -M{}{} -OXYZDAGgFPpRrSCc -I{}{}'.format(
        for line in utils.yield_cmd(
                'mblist -M{}{} -OXYZDSc -I{}{}'.format(
                    self.mb_exclude, ' {}'.format(
                        mb_region.format('gmt') if mb_region is not None else ''
                    ), mb_fn, ' -F{}'.format(mb_format) if mb_format is not None else ''
                ),
                verbose=False,
        ):
            this_line = [float(x) for x in line.strip().split('\t')]
            if self.auto_weight or self.auto_uncertainty:
                x = this_line[0]
                y = this_line[1]
                z = this_line[2]
                crosstrack_distance = this_line[3]
                #crosstrack_slope = this_line[4]
                #flat_bottom_grazing_angle = this_line[5]
                #seafloor_grazing_angle = this_line[6]
                #beamflag = this_line[7]
                #pitch = this_line[8]
                #draft = this_line[9]
                #roll = this_line[10]
                #heave = this_line[11]
                speed = this_line[4]
                #sonar_alt = this_line[13]
                sonar_depth = this_line[5]

                #utils.echo_msg(this_line)

                #if int(beamflag) == 0:# and abs(crosstrack_distance) < 1000:#abs(this_line[4]) < .15:

                xs.append(x)
                ys.append(y)
                zs.append(z)
                ## uncertainty
                #u_depth = 0
                u_depth = ((2+(0.02*(z*-1)))*0.51)
                u_s_depth = ((2+(0.02*(sonar_depth*-1)))*0.51)
                #u_depth = math.sqrt(1 + ((.023 * (z * -1))**2))
                u_cd = math.sqrt(1 + ((.023 * abs(crosstrack_distance))**2))
                #u_cd = ((2+(0.02*abs(crosstrack_distance)))*0.51) ## find better alg.
                #u_cd = 0
                u_s = math.sqrt(1 + ((.023 * abs(speed))**2))
                u = math.sqrt(u_depth**2 + u_cd**2 + u_s**2 + u_s_depth**2)
                us.append(u)
                if self.auto_weight:
                    ## weight

                    # if mb_perc is not None and mb_date is not None:
                    #     this_year = int(utils.this_year()) if self.min_year is None else self.min_year
                    #     w = float(mb_perc) * ((int(mb_date)-2000)/(this_year-2000))/100.            
                    #     w *= self.weight if self.weight is not None else 1
                    # else:
                    w = (1/u)
                    w *= self.weight if self.weight is not None else 1

                    ws.append(w)
            else:
                x = this_line[0]
                y = this_line[1]
                z = this_line[2]
                xs.append(x)
                ys.append(y)
                zs.append(z)
                        
        if len(xs) > 0:
            mb_points = np.column_stack((xs, ys, zs, ws, us))
            xs = ys = zs = ws = us = None
            mb_points = np.rec.fromrecords(mb_points, names='x, y, z, w, u')
            if self.want_binned:
                point_filter = PointFilterFactory(
                    mod='bin_z:percentile=25:y_res=3:z_res=20', points=mb_points
                )._acquire_module()
                if point_filter is not None:
                    mb_points = point_filter()
                    
            if mb_points is not None:
                yield(mb_points)
                
            mb_points = None

    def yield_mblist2_ds(self):
        """use mblist to process the multibeam data"""
        
        mb_fn = os.path.join(self.fn)
        if self.region is None and self.data_region is None:
            self.want_mbgrid = False

        if self.region is not None:
            mb_region = self.region.copy()
            mb_region.buffer(pct=25)
        else:
            mb_region = None

        mb_points = [[float(x) for x in l.strip().split('\t')] for l in utils.yield_cmd(
                'mblist -M{}{} -OXYZ -I{}'.format(
                    self.mb_exclude, ' {}'.format(
                        mb_region.format('gmt') if mb_region is not None else ''
                    ), mb_fn#, '-F{}'.format(self.mb_fmt) if self.mb_fmt is not None else ''
                ),
                verbose=False,
        )]

        if len(mb_points) > 0:
            mb_points = np.rec.fromrecords(mb_points, names='x, y, z')
        else:
            mb_points = None

        if self.want_binned:
            #mb_points = self.bin_z_points(mb_points)
            point_filter = PointFilterFactory(mod='bin_z', points=mb_points)._acquire_module()
            if point_filter is not None:
                mb_points = point_filter()

        if mb_points is not None:
            yield(mb_points)
                
    def yield_points(self):        
        mb_fn = os.path.join(self.fn)
        if self.region is None or self.data_region is None:
            self.want_mbgrid = False

        ## update want_mbgrid to output points!
        if self.want_mbgrid and (self.x_inc is not None and self.y_inc is not None):
            for pts in self.yield_mbgrid_ds():
                yield(pts)
        else:
            for pts in self.yield_mblist_ds():
                yield(pts)
                        
class OGRFile(ElevationDataset):
    """providing an OGR 3D point dataset parser.

    Useful for data such as S-57, ENC, E-Hydro, Etc.

    -----------
    Parameters:

    ogr_layer (str/int): the OGR layer containing elevation data
    elev_field (str): the field containing the z values
    weight_field (str): the field containing weight values
    uncertainty_field (str): the field containing uncertainty_values
    z_scale (float): scale the output z values    
    """

    _known_layer_names = ['SOUNDG', 'SurveyPoint_HD', 'SurveyPoint']
    _known_elev_fields = ['Elevation', 'elev', 'z', 'height', 'depth',
                          'topography', 'surveyPointElev', 'Z_depth',
                          'Z_height']
    
    def __init__(self,
                 ogr_layer = None,
                 elev_field = None,
                 weight_field = None,
                 uncertainty_field = None,
                 z_scale = None,
                 elevation_value = None,
                 **kwargs):
        super().__init__(**kwargs)
        self.ogr_layer = ogr_layer
        self.elev_field = elev_field
        self.weight_field = weight_field
        self.uncertainty_field = uncertainty_field
        self.z_scale = utils.float_or(z_scale)
        self.elevation_value = utils.float_or(elevation_value)

    def find_elevation_layer(self, ds_ogr):
        for l in self._known_layer_names:
            test_layer = ds_ogr.GetLayerByName(l)
            if test_layer is not None:
                return((l, test_layer))
        return(None, None)

    def find_elevation_field(self, ogr_feature):
        for f in self._known_elev_fields:            
            test_field = ogr_feature.GetField(f)
            if test_field is not None:
                return((l, test_field))
            
        return(None, None)
    
    def yield_points(self):
        ds_ogr = ogr.Open(self.fn)
        count = 0
        #utils.echo_msg(self.ogr_layer)
        if ds_ogr is not None:
            layer_name = None
            if self.ogr_layer is None:
                layer_name, layer_s = self.find_elevation_layer(ds_ogr)

                if layer_name is None:
                    layer_s = ds_ogr.GetLayer()
                    
            elif utils.str_or(self.ogr_layer) is not None:
                layer_name = self.ogr_layer
                layer_s = ds_ogr.GetLayer(str(self.ogr_layer))
                if layer_s is None:
                    layer_name, layer_s = self.find_elevation_layer(ds_ogr)
                    
            elif utils.int_or(self.ogr_layer) is not None:
                layer_s = ds_ogr.GetLayer(self.ogr_layer)
            else:
                layer_s = ds_ogr.GetLayer()

            if layer_s is not None:
                field_names = [field.name for field in layer_s.schema]
                if self.elev_field not in field_names:                    
                    for field_name in field_names:
                        if field_name in self._known_elev_fields:
                            self.elev_field = field_name
                            break
                    
                if self.region is not None:
                    layer_s.SetSpatialFilter(
                        self.region.export_as_geom() if self.transform['transformer'] is None else self.transform['trans_region'].export_as_geom()
                    )

                for f in layer_s:
                    geom = f.GetGeometryRef()
                    g = json.loads(geom.ExportToJson())
                    #utils.echo_msg(g)
                    xyzs = g['coordinates']
                    if geom.GetGeometryName() != 'MULTIPOINT':# and 'POLYGON' not in geom.GetGeometryName():
                        xyzs = [xyzs]

                    out_xyzs = []
                    for xyz in xyzs:
                        if not geom.Is3D():
                            # if self.elev_field is None:
                            #     self.elev_field = self.find_elevation_field(f)

                            if self.elev_field is None:
                                if self.elevation_value is None:
                                    elev = 0
                                else:
                                    elev = self.elevation_value
                            else:
                                elev = utils.float_or(f.GetField(self.elev_field))
                                
                            if elev is not None:                                
                                if isinstance(xyz[0], list):
                                    for x in xyz:
                                        for xx in x:
                                            #utils.echo_msg(xx)
                                            xx.append(elev)
                                else:
                                    xyz.append(elev)
                            else:
                                continue

                        else:
                            if self.elev_field is None:
                                if self.elevation_value is None:
                                    elev = 0
                                else:
                                    elev = self.elevation_value
                                    
                            else:
                                elev = utils.float_or(f.GetField(self.elev_field))
                                
                            if isinstance(xyz[0], list):
                                for x in xyz:
                                    for xx in x:
                                        xx[2] = elev

                    if isinstance(xyzs[0], list):
                        for x in xyzs:
                            if isinstance(x[0], list):
                                for xx in x:
                                    points = np.rec.fromrecords(xx, names='x, y, z')
                                    if self.z_scale is not None:
                                        points['z'] *= self.z_scale

                                    yield(points)
                            else:
                                points = np.rec.fromrecords([x], names='x, y, z')
                                if self.z_scale is not None:
                                    points['z'] *= self.z_scale

                                yield(points)
                                
                    # #points = np.rec.fromrecords(out_xyzs, names='x, y, z')
                    # if self.z_scale is not None:
                    #     points['z'] *= self.z_scale

                    # yield(points)
                            
            ds_ogr = layer_s = None

class Points(ElevationDataset):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.metadata['name'] = 'points'
        self.data_format = -4

    def yield_points(self):
        points = self.fn
        if isinstance(points, np.ndarray):
            points = np.rec.fromrecords(points, names='x, y, z')
        elif isinstance(points, np.core.records.recarray):
            points = points
        elif isinstance(points, pd.DataFrame):
            points = points
        
        yield(points)
            
class Scratch(ElevationDataset):
    """Scratch Dataset

    Process a python list of valid dataset entries
    """
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.metadata['name'] = 'scratch'
        self.data_format = -3

    def generate_inf(self):
        """generate and return the infos blob of the datalist.

        datasets parsed from the datalist may have variable srs, so
        we want to transform those to whatever self.dst_srs is, if it
        exists, in order to properly fill the datalist/json srs and regions...
        """

        _region = self.region
        self.region = None
        out_regions = []
        out_srs = []
        for entry in self.parse():
            entry_minmax = entry.infos.minmax
            ## entry has an srs and dst_srs is set, so lets transform the region to suit
            if entry.src_srs is not None:
                out_srs.append(entry.src_srs)
                if self.dst_srs is not None:
                    entry_region = regions.Region().from_list(entry_minmax)
                    if entry_region.valid_p():
                        entry_region.src_srs = entry.src_srs
                        entry_region.warp(self.dst_srs)
                        entry_minmax = entry_region.export_as_list(include_z=True)
                        
            entry_region = regions.Region().from_list(entry_minmax)
            if entry_region.valid_p():
                out_regions.append(entry_region)
                self.infos.numpts += entry.infos.numpts

        ## merge all the gathered regions
        region_count = 0
        out_region = None
        for this_region in out_regions:
            if this_region.valid_p():
                if region_count == 0:
                    out_region = this_region
                    region_count += 1
                else:
                    out_region = regions.regions_merge(out_region, this_region)
                    
        if out_region is not None:
            self.infos.minmax = out_region.export_as_list(include_z=True)
            self.infos.wkt = out_region.export_as_wkt()

        ## set the epsg for the datalist
        if self.infos.src_srs is None:
            if self.src_srs is not None:
                self.infos.src_srs = self.src_srs
            elif out_srs:
                if all(x == out_srs[0] for x in out_srs):
                    self.infos.src_srs = out_srs[0]

                self.src_srs = self.infos.src_srs
                
        self.region = _region
        return(self.infos)
    
    def parse(self):
        if isinstance(self.fn, list):
            for this_ds in self.fn:
                if this_ds is not None and this_ds.valid_p(
                        fmts=DatasetFactory._modules[this_ds.data_format]['fmts']
                ):
                    this_ds.initialize()
                    for ds in this_ds.parse():
                        # fill self.data_entries with each dataset for use outside the yield.
                        self.data_entries.append(ds)
                        yield(ds)
                        
class Datalist(ElevationDataset):
    """representing a datalist parser
    
    A datalist is an extended MB-System style datalist.
    
    Each datalist consists of datalist-entries, where a datalist entry has the following columns:
    path format weight uncertainty title source date data_type resolution hdatum vdatum url

    the datalist can contain datalist-entries to other datalist files, distributed across a
    file-system.

    see `cudem.dlim.datasets` for superclass ElevationDataset
    """

    _datalist_json_cols = ['path', 'format', 'weight', 'uncertainty', 'name', 'title', 'source',
                           'date', 'data_type', 'resolution', 'hdatum', 'vdatum',
                           'url', 'mod_args']

    def __init__(self, fmt=None, **kwargs):
        super().__init__(**kwargs)
        
    def _init_datalist_vector(self):
        """initialize the datalist geojson vector.
        
        this vector is used to quickly parse data from distributed files
        across a file-system. Each data file will have it's own entry in
        the geojson datalist vector, containing gathered metadata based on the
        source datalist.
        """

        #self.set_transform()
        self.dst_layer = '{}'.format(self.fn)
        self.dst_vector = self.dst_layer + '.json'

        utils.remove_glob('{}.json'.format(self.dst_layer))
        #if self.src_srs is not None:
        #    gdalfun.osr_prj_file('{}.prj'.format(self.dst_layer), self.src_srs)
            
        self.ds = ogr.GetDriverByName('GeoJSON').CreateDataSource(self.dst_vector)
        srs = srsfun.osr_srs(self.src_srs)
        
        if self.ds is not None: 
            self.layer = self.ds.CreateLayer(
                '{}'.format(self.dst_layer), srs, ogr.wkbMultiPolygon
            )
            [self.layer.CreateField(
                ogr.FieldDefn('{}'.format(f), ogr.OFTString)
            ) for f in self._datalist_json_cols]
            
            [self.layer.SetFeature(feature) for feature in self.layer]
            return(0)
        else:
            self.layer = None
            return(-1)

    def _create_entry_feature(self, entry, entry_region):
        """create a datalist entry feature and insert it into the datalist-vector geojson

        -----------
        Parameters:
        
        entry - the datalist entry object
        entry_region - the region of the datalist entry
        """

        #utils.echo_msg(entry)
        #utils.echo_msg(entry.params['kwargs'])
        entry_path = os.path.abspath(entry.fn) if not entry.remote else entry.fn
        entry_fields = [entry_path,
                        entry.data_format,
                        entry.weight,
                        entry.uncertainty,
                        entry.metadata['name'],
                        entry.metadata['title'],
                        entry.metadata['source'],
                        entry.metadata['date'],
                        entry.metadata['data_type'],
                        entry.metadata['resolution'],
                        entry.metadata['hdatum'],
                        entry.metadata['vdatum'],
                        entry.metadata['url'],
                        utils.dict2args(entry.params['mod_args'])]
        dst_defn = self.layer.GetLayerDefn()
        entry_geom = ogr.CreateGeometryFromWkt(entry_region.export_as_wkt())
        out_feat = ogr.Feature(dst_defn)
        out_feat.SetGeometry(entry_geom)
        for i, f in enumerate(self._datalist_json_cols):
            out_feat.SetField(f, entry_fields[i])
            
        self.layer.CreateFeature(out_feat)

    def generate_inf(self):
        """generate and return the infos blob of the datalist.

        datasets parsed from the datalist may have variable srs, so
        we want to transform those to whatever self.dst_srs is, if it
        exists, in order to properly fill the datalist/json srs and regions...
        """

        self.region = None
        self.infos.file_hash = self.infos.generate_hash()
        _region = self.region
        out_regions = []
        out_srs = []

        ## attempt to generate a datalist-vector geojson and
        ## if successful, fill it wil the datalist entries, using `parse`
        if self._init_datalist_vector() == 0:
            for entry in self.parse():
                entry_minmax = entry.infos.minmax
                if entry.mask is not None: ## add all duplicate params???
                    entry.params['mod_args']['mask'] = entry.mask['mask']
                    for key in entry.mask.keys():
                        entry.params['mod_args'][key] = entry.mask[key]
                        #entry.params['mod_args'] = {'mask': entry.mask}

                ## entry has an srs and dst_srs is set, so lets transform the region to suit
                if entry.src_srs is not None:
                    out_srs.append(entry.src_srs)
                    if self.dst_srs is not None:
                        entry_region = regions.Region(
                            src_srs=entry.src_srs
                        ).from_list(entry_minmax).warp(self.dst_srs)
                        entry_minmax = entry_region.export_as_list(include_z=True)

                ## create the feature for the geojson
                entry_region = regions.Region().from_list(entry_minmax)
                if entry_region.valid_p():
                    out_regions.append(entry_region)
                    self._create_entry_feature(entry, entry_region)
                    self.infos.numpts += entry.infos.numpts

            self.ds = self.layer = None # close the geojson ogr dataset

        else:
            utils.echo_warning_msg('could not initialize datalist vector')
            return(self.infos)

        ## merge all the gathered regions
        region_count = 0
        out_region = None
        for this_region in out_regions:
            if this_region.valid_p():
                if region_count == 0:
                    out_region = this_region
                    region_count += 1
                else:
                    out_region = regions.regions_merge(out_region, this_region)
                    
        if out_region is not None:
            self.infos.minmax = out_region.export_as_list(include_z=True)
            self.infos.wkt = out_region.export_as_wkt()

        ## set the epsg for the datalist
        if self.infos.src_srs is None:
            if self.src_srs is not None:
                self.infos.src_srs = self.src_srs
            elif out_srs:
                if all(x == out_srs[0] for x in out_srs):
                    self.infos.src_srs = out_srs[0]

                self.src_srs = self.infos.src_srs

        self.region = _region
        return(self.infos)
        
    def parse_json(self):
        """parse the datalist using the datalist-vector geojson.

        Quickly find data in a given region using the datalist-vector. The datalist-vector
        must have been previously generated using `parse`. If the datlist-vector is not found
        will fall-back to `parse` and generate a new datalist-vector geojson.

        -------
        Yields:
        dataset object of each dataset found in the datalist-vector
        """

        ## check for the datalist-vector geojson
        status = 0
        count = 0
        
        ## user input to re-gerenate json...?
        if os.path.exists('{}.json'.format(self.fn)):
            driver = ogr.GetDriverByName('GeoJSON')
            dl_ds = driver.Open('{}.json'.format(self.fn))
            if dl_ds is None:
                #utils.echo_warning_msg('could not open {}.json'.format(self.fn))
                status = -1
        else:
            status = -1

        ## parse the datalist-vector geojson and yield the results
        if status != -1:
            dl_layer = dl_ds.GetLayer()
            ldefn = dl_layer.GetLayerDefn()
            _boundsGeom = None
            if self.region is not None:
                _boundsGeom = self.region.export_as_geom() \
                    if self.transform['transformer'] is None \
                       else self.transform['trans_region'].export_as_geom()

            dl_layer.SetSpatialFilter(_boundsGeom)
            count = len(dl_layer)
            with tqdm(
                    total=len(dl_layer),
                    desc='parsing {} datasets from datalist json {}.json @ {}'.format(
                        count, self.fn, self.weight
                    ),
                    leave=self.verbose
            ) as pbar:
                for l,feat in enumerate(dl_layer):
                    pbar.update()
                    ## filter by input source region extras (weight/uncertainty)
                    if self.region is not None:
                        w_region = self.region.w_region()
                        if w_region[0] is not None:
                            if float(feat.GetField('weight')) < w_region[0]:
                                continue

                        if w_region[1] is not None:
                            if float(feat.GetField('weight')) > w_region[1]:
                                continue

                        u_region = self.region.u_region()
                        if u_region[0] is not None:
                            if float(feat.GetField('uncertainty')) < u_region[0]:
                                continue

                        if u_region[1] is not None:
                            if float(feat.GetField('uncertainty')) > u_region[1]:
                                continue

                    ## extract the module arguments from the datalist-vector
                    try:
                        ds_args = feat.GetField('mod_args')
                        data_set_args = utils.args2dict(list(ds_args.split(':')), {})
                        # for kpam in list(data_set_args.keys()):
                        #     if kpam in self.__dict__:
                        #         kval = data_set_args[kpam]
                        #         self.__dict__[kpam] = kval
                        #         del data_set_args[kpam]

                        ds_args = utils.dict2args(data_set_args)
                        # for kpam, kval in data_set_args.items():
                        #     if kpam in self.__dict__:
                        #         utils.echo_msg(kpam)
                        #         self.__dict__[kpam] = kval
                        #         #del data_set_args[kpam]
                    except:
                        ds_args = None
                        data_set_args = {}

                    ## update existing metadata
                    md = copy.deepcopy(self.metadata)
                    for key in self.metadata.keys():
                        md[key] = feat.GetField(key)

                    ## generate the dataset object to yield
                    data_mod = '"{}" {}{} {} {}'.format(
                        feat.GetField('path'),
                        feat.GetField('format'),
                        ':{}'.format(ds_args) if ds_args is not None else '',
                        feat.GetField('weight'),
                        feat.GetField('uncertainty')
                    )
                    data_set = DatasetFactory(
                        **self._set_params(mod=data_mod, metadata=md)#, **data_set_args)
                    )._acquire_module()
                    if data_set is not None and data_set.valid_p(
                            fmts=DatasetFactory._modules[data_set.data_format]['fmts']
                    ):
                        data_set.initialize()
                        #utils.echo_msg(data_set.params)
                        ## fill self.data_entries with each dataset for use outside the yield.
                        for ds in data_set.parse(): 
                            self.data_entries.append(ds)
                            yield(ds)

            dl_ds = dl_layer = None
                
        else:
            ## failed to find/open the datalist-vector geojson, so run `parse` instead and
            ## generate one for future use...
            # utils.echo_warning_msg(
            #     'could not load datalist-vector json {}.json, falling back to parse, generate a json file for the datalist using `dlim -i`'.format(self.fn)
            # )
            for ds in self.parse_no_json():
                yield(ds)
                                        
    def parse(self):
        """parse the datalist file.

        -------
        Yields:
        dataset object of each dataset found in the datalist        
        """

        status = 0
        if os.path.exists(self.fn):
            with open(self.fn, 'r') as f:
                count = sum(1 for _ in f)

            with open(self.fn, 'r') as op:
                with tqdm(
                        total=count,
                        desc='parsing datalist {}...'.format(self.fn),
                        leave=self.verbose
                ) as pbar:
                    for l, this_line in enumerate(op):
                        pbar.update()
                        ## parse the datalist entry line
                        if this_line[0] != '#' and this_line[0] != '\n' and this_line[0].rstrip() != '':
                            md = copy.deepcopy(self.metadata)
                            md['name'] = utils.fn_basename2(os.path.basename(self.fn))
                            
                            ## generate the dataset object to yield
                            data_set = DatasetFactory(
                               **self._set_params(
                                   mod=this_line,
                                   metadata=md,
                                   src_srs=self.src_srs,
                                   parent=self,
                                   fn=None
                               )
                            )._acquire_module()
                            if data_set is not None and data_set.valid_p(
                                    fmts=DatasetFactory._modules[data_set.data_format]['fmts']
                            ):
                                data_set.initialize()
                                ## filter with input source region, if necessary
                                ## check input source region against the dataset region found
                                ## in its inf file.
                                if self.region is not None and self.region.valid_p(check_xy=True):
                                    inf_region = regions.Region().from_list(data_set.infos.minmax)
                                    if inf_region.valid_p():
                                        inf_region.wmin = data_set.weight
                                        inf_region.wmax = data_set.weight
                                        inf_region.umin = data_set.uncertainty
                                        inf_region.umax = data_set.uncertainty

                                        if not regions.regions_intersect_p(
                                                inf_region,
                                                self.region if data_set.transform['transformer'] is None else data_set.transform['trans_region']
                                        ):
                                            continue
                                        
                                ## fill self.data_entries with each dataset for use outside the yield and
                                ## yield the dataset object.
                                for ds in data_set.parse(): 
                                    self.data_entries.append(ds)
                                    yield(ds)

        ## self.fn is not a file-name, so check if self.data_entries not empty
        ## and return the dataset objects found there.
        elif len(self.data_entries) > 0:
            for data_set in self.data_entries:
                for ds in data_set.parse():
                    yield(ds)
        else:
            if self.verbose:
                utils.echo_warning_msg(
                    'could not open datalist/entry {}'.format(self.fn)
                )
                
class ZIPlist(ElevationDataset):
    """Zip file parser.

    Parse supported datasets from a zipfile.
    """
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def generate_inf(self):
        """generate and return the infos blob of the datalist.

        datasets parsed from the datalist may have variable srs, so
        we want to transform those to whatever self.dst_srs is, if it
        exists, in order to properly fill the datalist/json srs and regions...
        """

        _region = self.region
        self.region = None
        out_regions = []
        out_srs = []
        self.infos.file_hash = self.infos.generate_hash()
        for entry in self.parse():            
            entry_minmax = entry.infos.minmax

            ## entry has an srs and dst_srs is set, so lets transform the region to suit
            if entry.src_srs is not None:
                out_srs.append(entry.src_srs)
                if self.dst_srs is not None:
                    entry_region = regions.Region().from_list(entry_minmax)
                    entry_region.src_srs = entry.src_srs
                    entry_region.warp(self.dst_srs)
                    entry_minmax = entry_region.export_as_list(include_z=True)

            entry_region = regions.Region().from_list(entry_minmax)
            if entry_region.valid_p():
                out_regions.append(entry_region)
                self.infos.numpts += entry.infos.numpts

        ## merge all the gathered regions
        region_count = 0
        out_region = None
        for this_region in out_regions:
            if this_region.valid_p():
                if region_count == 0:
                    out_region = this_region
                    region_count += 1
                else:
                    out_region = regions.regions_merge(out_region, this_region)
                    
        if out_region is not None:
            self.infos.minmax = out_region.export_as_list(include_z=True)
            self.infos.wkt = out_region.export_as_wkt()

        ## set the epsg for the datalist
        if self.infos.src_srs is None:
            if self.src_srs is not None:
                self.infos.src_srs = self.src_srs
            elif out_srs:
                if all(x == out_srs[0] for x in out_srs):
                    self.infos.src_srs = out_srs[0]

                self.src_srs = self.infos.src_srs
                
        self.region = _region
        return(self.infos)
                
    def parse(self):
        import zipfile
        exts = [DatasetFactory()._modules[x]['fmts'] for x in DatasetFactory()._modules.keys()]
        exts = [x for y in exts for x in y]
        datalist = []
        if self.fn.split('.')[-1].lower() == 'zip':
            try:
                with zipfile.ZipFile(self.fn) as z:
                    zfs = z.namelist()
                    for ext in exts:
                        for zf in zfs:
                            if ext == zf.split('.')[-1]:
                                datalist.append(os.path.basename(zf))
            except Exception as e:
                utils.echo_error_msg('could not unzip {}, {}'.format(self.fn, e))
                            
        for this_data in datalist:
            this_line = utils.p_f_unzip(self.fn, fns=[this_data], outdir=os.path.dirname(self.fn), tmp_fn=True)[0]
            #utils.echo_msg(this_line)
            data_set = DatasetFactory(
                **self._set_params(
                    mod=os.path.basename(this_line),
                    data_format=None,
                    src_srs=self.src_srs,
                    parent=self
                )
            )._acquire_module()
            #utils.echo_msg(data_set)
            #sys.exit()
            if data_set is not None and data_set.valid_p(
                    fmts=DatasetFactory._modules[data_set.data_format]['fmts']
            ):
                data_set.initialize()
                if self.region is not None and self.region.valid_p(check_xy=True):
                    try:
                        inf_region = regions.Region().from_string(
                            data_set.infos.wkt
                        )
                    except:
                        try:
                            inf_region = regions.Region().from_list(
                                data_set.infos.minmax
                            )
                        except:
                            inf_region = self.region.copy()
                            
                    inf_region.wmin = data_set.weight
                    inf_region.wmax = data_set.weight
                    inf_region.umin = data_set.uncertainty
                    inf_region.umax = data_set.uncertainty
                    if inf_region.valid_p(check_xy=True):
                        if regions.regions_intersect_p(
                                inf_region,
                                self.region if data_set.transform['transformer'] is None else data_set.transform['trans_region']
                        ):
                            for ds in data_set.parse():
                                self.data_entries.append(ds)
                                yield(ds)
                    else:
                        if self.verbose:
                            utils.echo_warning_msg(
                                'invalid inf file: {}.inf, skipping'.format(
                                    data_set.fn
                                )
                            )
                            utils.remove_glob('{}*'.format(data_set.fn))
                            
                else:
                    for ds in data_set.parse():
                        self.data_entries.append(ds)
                        yield(ds)
                        
            utils.remove_glob('{}*'.format(data_set.fn))
            #utils.remove_glob('{}*'.format(this_line))
        
class Fetcher(ElevationDataset):
    """The default fetches dataset type; dlim Fetcher dataset class

    This is used in waffles/dlim for on-the-fly remote data
    parsing and processing.
    
    See `fetches`, `cudem.fetches` for more information on usage and the
    various fetches modules supported.

    If a fetch module needs special processing define a sub-class
    of Fetcher and redefine the yield_ds(self, result) function which yields a
    list of dlim dataset objects, where result is an item from the fetch result list.
    Otherwise, this Fetcher class can be used as default if the fetched data comes
    in a normal sort of way.

    Generally, though not always, if the fetched data is a raster then
    there is no need to redefine yield_ds, though if the raster has insufficient
    information, such as with Copernicus, whose nodata value is not
    specified in the geotiff files, it may be best to create a simple
    sub-class for it.
    """

    def __init__(self,
                 keep_fetched_data = True,
                 outdir = None,
                 check_size = True,
                 want_single_metadata_name = False,
                 callback = fetches.fetches_callback,
                 **kwargs
    ):
        super().__init__(**kwargs)
        ## cache directory to store fetched data
        #self.outdir = outdir if outdir is not None else self.cache_dir 
        self.fetch_module = fetches.FetchesFactory(
            mod=self.fn, src_region=self.region, callback=callback,
            verbose=self.verbose, outdir=outdir
        )._acquire_module() # the fetches module
        if self.fetch_module is None:
            utils.echo_warning_msg('fetches modules {} returned None'.format(self.fn))
            pass
        
        self.want_single_metadata_name = want_single_metadata_name
        self.check_size = check_size # check the size of the fetched data
        self.keep_fetched_data = keep_fetched_data # retain fetched data after processing
        self.outdir = self.fetch_module._outdir
        self.cache_dir = self.outdir                
        try:
            self.fetch_module.run()
        except:
            self.fetch_module.results = []

        ## breaks when things not set...
        # src_horz, src_vert = gdalfun.epsg_from_input(self.fetch_module.src_srs)

        ## set the metadata from the fetches module
        md = copy.deepcopy(self.metadata)
        md['name'] = self.metadata['name']
        md['title'] = self.fetch_module.title
        md['source'] = self.fetch_module.source
        md['date'] = self.fetch_module.date
        md['data_type'] = self.data_format
        md['resolution'] = self.fetch_module.resolution
        md['hdatum'] = self.fetch_module.hdatum
        md['vdatum'] = self.fetch_module.vdatum
        md['url'] = self.fetch_module.url

        self.fetches_params = self._set_params(
            data_format=self.fetch_module.data_format,
            src_srs=self.fetch_module.src_srs,
            cache_dir=self.fetch_module._outdir,
            remote=True,
            metadata=md,
            parent=self,
        )

    def _reset_params(self):
        self.fetches_params = self._set_params(**self.fetches_params)
        
    def generate_inf(self):
        """generate a infos dictionary from the Fetches dataset"""

        tmp_region = self.fetch_module.region if self.region is None else self.region.copy()
        self.infos.minmax = tmp_region.export_as_list()    
        self.infos.wkt = tmp_region.export_as_wkt()
        self.infos.src_srs = self.fetch_module.src_srs
        return(self.infos)

    def parse(self):
        with tqdm(
                total=len(self.fetch_module.results),
                desc='parsing datasets from datalist fetches {} @ {}'.format(
                    self.fetch_module, self.weight
                ),
            leave=self.verbose
        ) as pbar:
            for result in self.fetch_module.results:
                status = self.fetch_module.fetch(result, check_size=self.check_size)
                if status == 0:
                    self.fetches_params['mod'] = os.path.join(self.fetch_module._outdir, result['dst_fn'])
                    for this_ds in self.yield_ds(result):
                        if this_ds is not None:
                            f_name = os.path.relpath(this_ds.fn.split(':')[0], self.fetch_module._outdir)
                            if f_name == '.':
                                f_name = this_ds.fn.split(':')[0]

                            mod_name = os.path.dirname(utils.fn_basename2(f_name))
                            if mod_name == '':
                                mod_name = self.fetch_module.name

                            #if self.want_single_metadata_name:
                            #    this_ds.metadata['name'] = mod_name
                            #else:
                            this_ds.metadata['name'] = '/'.join(['/'.join(this_ds.metadata['name'].split('/')[:-1]), f_name])
                            this_ds.remote = True
                            this_ds.initialize()
                            for ds in this_ds.parse():
                                yield(ds)
                        else:
                            utils.echo_warning_msg(
                                'could not set fetches datasource {}'.format(result)
                            )
                else:
                    utils.echo_warning_msg(
                        'data not fetched {}:{}'.format(status, result)
                    )
                        
                pbar.update()
                
        if not self.keep_fetched_data:
            utils.remove_glob('{}*'.format(self.fn))
            
    def yield_ds(self, result):
        ## try to get the SRS info from the result if it's a gdal file
        ## fix this.
        try:
            vdatum = self.fetch_module.vdatum
            src_srs = gdalfun.gdal_get_srs(os.path.join(self.fetch_module._outdir, result['dst_fn']))
            horz_epsg, vert_epsg = gdalfun.epsg_from_input(src_srs)
            if vert_epsg is None:
                vert_epsg = vdatum

            if vert_epsg is not None:
                self.fetch_module.src_srs = '{}+{}'.format(horz_epsg, vert_epsg)
            else:
                self.fetch_module.src_srs = src_srs
        except:
            pass

        yield(DatasetFactory(**self.fetches_params)._acquire_module())

class NEDFetcher(Fetcher):
    """National Elevation Dataset from USGS

    This is a wrapper shortcut for fetching NED DEMs from USGS' The National Map (TNM)
    """

    __doc__ = '''{}    
    Fetches Module: <ned> - {}'''.format(__doc__, fetches.NED.__doc__)
    
    def __init__(self, coast_buffer=0.00001, **kwargs):
        super().__init__(**kwargs)
        self.coast_buffer = utils.float_or(coast_buffer, 0.00001)
        
    def yield_ds(self, result):
        ## coastline generation
        # from cudem import waffles
        # coast_mask_fn = utils.make_temp_fn('ned_coast')
        # _coast_mask = waffles.WaffleFactory(mod='coastline', src_region=self.region, xinc=self.x_inc, yinc=self.y_inc,
        #                                     want_stack=False, name=coast_mask_fn, node='pixel', verbose=True)._acquire_module()()        
        # self.mask = '{}.shp'.format(_coast_mask.name)
        
        ## todo: merge the coast mask with user input self.mask
        coast_mask = None
        ned_mask = self.mask
        if self.mask is None:
            coast_mask = self.process_coastline(
                self.fetch_coastline(
                    chunks=False, verbose=False
                ),
                return_geom=False,
                landmask_is_watermask=True,
                include_landmask=False,
                line_buffer=self.coast_buffer,
                verbose=False
            )
            if coast_mask is not None:
                ned_mask = {'mask': coast_mask, 'invert_mask': True}
        
        src_dem = os.path.join(self.fetch_module._outdir, result['dst_fn'])
        #ned_metadata = copy.deepcopy(self.metadata)
        #ned_metadata['name'] = src_dem
        #self.fetches_params['metadata'] = ned_metadata
        self.fetches_params['mod'] = src_dem
        self.fetches_params['mask'] = ned_mask
        self.fetches_params['remove_flat'] = True
        yield(DatasetFactory(**self.fetches_params)._acquire_module())

        if coast_mask is not None:
            utils.remove_glob('{}.*'.format(utils.fn_basename2(coast_mask)))

class DNRFetcher(Fetcher):
    """
    """

    __doc__ = '''{}    
    Fetches Module: <wadnr> - {}'''.format(__doc__, fetches.waDNR.__doc__)
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        self.data_format = -2
        
    def yield_ds(self, result):
        src_dnr_dems = utils.p_unzip(
            os.path.join(self.fetch_module._outdir, result['dst_fn']),
            exts=['tif'],
            outdir=self.fetch_module._outdir,
            verbose=self.verbose
        )
        for src_dnr_dem in src_dnr_dems:
            self.fetches_params['mod'] = src_dnr_dem
            #self.fetches_params['data_format'] = 200
            #self.fetches_params['node'] = 'pixel'
            self.fetches_params['remove_flat'] = True
            yield(DatasetFactory(**self.fetches_params)._acquire_module())
            
class DAVFetcher_CoNED(Fetcher):
    """CoNED from the digital coast 

    This is a wrapper shortcut for fetching CoNED DEMs from the Digital Coast,
    mainly so we can pull the vertical datum info from the DAV metadata since the
    CoNED doesn't assign one to their DEMs.
    """

    __doc__ = '''{}    
    Fetches Module: <CoNED> - {}'''.format(__doc__, fetches.CoNED.__doc__)
    
    
    def __init__(self, keep_fetched_data = True, cog = True, **kwargs):
        super().__init__(**kwargs)
        self.keep_fetched_data = keep_fetched_data
        self.cog = cog

    def parse(self):
        #self.fetch_module.run()
        for result in self.fetch_module.results:
            if not self.cog:
                status = self.fetch_module.fetch(result, check_size=self.check_size)
                if status != 0:
                    break
                
            for this_ds in self.yield_ds(result):
                if this_ds is not None:
                    #this_ds.metadata['name'] = utils.fn_basename2(this_ds.fn)
                    this_ds.initialize()
                    yield(this_ds)
            
    def yield_ds(self, result):
        ## try to get the SRS info from the result
        try:
            vdatum = self.fetch_module.vdatum
            if not self.cog:
                src_srs = gdalfun.gdal_get_srs(os.path.join(self.fetch_module._outdir, result['dst_fn']))
            else:
                src_srs = gdalfun.gdal_get_srs(result['url'])
                
            horz_epsg, vert_epsg = gdalfun.epsg_from_input(src_srs)
            if vert_epsg is None:
                vert_epsg = vdatum

            #utils.echo_msg('srs: {}+{}'.format(horz_epsg, vert_epsg))
            if vert_epsg is not None:
                self.fetch_module.src_srs = '{}+{}'.format(horz_epsg, vert_epsg)
            else:
                self.fetch_module.src_srs = src_srs
        except:
            pass

        self.fetches_params['mod'] = os.path.join(self.fetch_module._outdir, result['dst_fn']) \
            if not self.cog \
               else result['url']
        self.fetches_params['check_path'] = False if self.cog else True
        self.fetches_params['src_srs'] = self.fetch_module.src_srs
        self.fetches_params['data_format'] = 200
        yield(DatasetFactory(**self.fetches_params)._acquire_module())

class DAVFetcher_SLR(Fetcher):
    """SLR DEM from the digital coast 

    This is a wrapper shortcut for fetching SLR DEMs from the Digital Coast,
    mainly so we can pull and remove the flattened ring around the actual data.
    """

    __doc__ = '''{}    
    Fetches Module: <SLR> - {}'''.format(__doc__, fetches.SLR.__doc__)
    
    def __init__(self, keep_fetched_data = True, **kwargs):
        super().__init__(**kwargs)

    def yield_ds(self, result):
        ## this doesn't work in all cases, update to find and remove flattened areas
        #gdalfun.gdal_set_ndv(os.path.join(self.fetch_module._outdir, result[1]), ndv=-99.0000, convert_array=True)
        self.fetches_params['remove_flat'] = True
        yield(DatasetFactory(**self.fetches_params)._acquire_module())    
        
class SWOTFetcher(Fetcher):
    """SWOT L2_HR_Raster data from NASA

    set `data_set` to one of (for L2_HR_Raster):
    1 - longitude (64-bit floating-point)
    2 - latitude (64-bit floating-point)
    3 - wse (32-bit floating-point)
    4 - status_flag (8-bit unsigned integer)
    5 - status_flag (32-bit unsigned integer)
    6 - wse_uncert (32-bit floating-point)
    7 - water_area (32-bit floating-point)
    8 - status_flag (8-bit unsigned integer)
    9 - status_flag (32-bit unsigned integer)
    10 - water_area_uncert (32-bit floating-point)
    11 - water_frac (32-bit floating-point)
    12 - water_frac_uncert (32-bit floating-point)
    13 - sig0 (32-bit floating-point)
    14 - status_flag (8-bit unsigned integer)
    15 - status_flag (32-bit unsigned integer)
    16 - sig0_uncert (32-bit floating-point)
    17 - inc (32-bit floating-point)
    18 - cross_track (32-bit floating-point)
    19 - time (64-bit floating-point)
    20 - time (64-bit floating-point)
    21 - n_wse_pix (32-bit unsigned integer)
    22 - n_water_area_pix (32-bit unsigned integer)
    23 - n_sig0_pix (32-bit unsigned integer)
    24 - n_other_pix (32-bit unsigned integer)
    25 - dark_frac (32-bit floating-point)
    26 - status_flag (8-bit unsigned integer)
    27 - status_flag (8-bit unsigned integer)
    28 - layover_impact (32-bit floating-point)
    29 - sig0_cor_atmos_model (32-bit floating-point)
    30 - height_cor_xover (32-bit floating-point)
    31 - geoid_height_above_reference_ellipsoid (32-bit floating-point)
    32 - solid_earth_tide (32-bit floating-point)
    33 - load_tide_fes (32-bit floating-point)
    34 - load_tide_got (32-bit floating-point)
    35 - pole_tide (32-bit floating-point)
    36 - model_dry_tropo_cor (32-bit floating-point)
    37 - model_wet_tropo_cor (32-bit floating-point)
    38 - iono_cor_gim_ka (32-bit floating-point)

    currently supported SWOT products and product options:
    - L2_HR_Raster
      apply_geoid
      data_set

    - L2_HR_PIXC
      apply_geoid
      
    """

    __doc__ = '''{}    
    Fetches Module: <swot> - {}'''.format(__doc__, fetches.SWOT.__doc__)
    
    def __init__(self,
                 data_set = 'wse',
                 apply_geoid = True,
                 classes = None,
                 classes_qual = None,
                 anc_classes = None,
                 remove_class_flags = False,
                 **kwargs):
        super().__init__(**kwargs)
        self.fetches_params['var'] = data_set
        self.fetches_params['apply_geoid'] = apply_geoid
        self.fetches_params['classes'] = classes
        self.fetches_params['anc_classes'] = anc_classes
        self.fetches_params['classes_qual'] = clases_qual
        self.fetches_params['remove_class_flags'] = remove_class_flags
        
    # def fetch_pixc_vec(self, swot_fn):
    #     pixc_vec_filter = utils.fn_basename2(swot_fn).split('PIXC_')[1]
    #     this_pixc_vec = fetches.SWOT(
    #         src_region=None, verbose=self.verbose, outdir=self.fetch_module._outdir, product='L2_HR_PIXCVec', filename_filter=pixc_vec_filter
    #     )
    #     this_pixc_vec.run()
    #     if len(this_pixc_vec.results) == 0:
    #         utils.echo_warning_msg('could not locate associated PIXCVec file for {}'.format(pixc_vec_filter))
    #         return(None)
    #     else:
    #         if this_pixc_vec.fetch(this_pixc_vec.results[0], check_size=self.check_size) == 0:
    #             return(os.path.join(this_pixc_vec._outdir, this_pixc_vec.results[0][1]))
        
    def yield_ds(self, result):
        #swot_fn = os.path.join(self.fetch_module._outdir, result[1])
        if 'L2_HR_PIXC_' in result['data_type']:
            #pixc_vec_result = self.fetch_pixc_vec(swot_fn)
            #swot_pixc_vec_fn = pixc_vec_result
            self.fetches_params['data_format'] = 202
            yield(DatasetFactory(**self.fetches_params).acquire_module())
        elif 'L2_HR_Raster' in result['data_type']:
            self.fetch_module['data_format'] = 203
            yield(DatasetFactory(**self.fetches_params).acquire_module())
        else:
            utils.echo_warning_msg('{} is not a currently supported dlim dataset'.format(result['data_type']))

class IceSat2Fetcher(Fetcher):
    """IceSat2 data from NASA

    See `fetches --modules icesat2` for fetching parameters

    Parameters:
    
    water_surface: 'geoid' # this is the vertical datum, can be 'geoid', 'ellipsoid' or 'mean_tide'
    classes: None # return only data with the specified classes, e.g. '2/3/4'
    confidence_levels: None # return only data with the specified confidence levels, e.g. '2/3/4'
    columns: {} # the additional columns to export in transform_and_yield_points
    classify_bathymetry: True # extract bathymetry with CShelph
    classify_buildings: True # classify buildings with BING BFP
    classify_water: True # classify water using OSM
    reject_failed_qa: True # skip granules that failed QA

    Classes:
    -1 - no classification (ATL08)
    0 - noise / atmosphere (ATL08)
    1 - ground surface (ATL08)
    2 - canopy (ATL08)
    3 - canopy top (ATL08)
    4 - bathymetry floor surface (CShelph, future ATL24)
    5 - bathymetry water surface (CShelph, future ATL24)
    6 - ice surface (ATL06) (unused for now, just planning ahead for possible future ATL06 integration)
    7 - built structure (OSM or Bing)
    8 - "urban" (WSF, if used)

    Confidence Levels:
    0, 1, 2, 3, 4
    """

    __doc__ = '''{}    
    Fetches Module: <icesat2> - {}'''.format(__doc__, fetches.IceSat2.__doc__)
    
    def __init__(self,
                 water_surface = 'geoid',
                 classes = None,
                 confidence_levels = None,
                 columns = {},
                 classify_bathymetry = True,
                 classify_buildings = True,
                 classify_water = True,
                 reject_failed_qa = True,
                 **kwargs):
        super().__init__(**kwargs)
        self.fetches_params['water_suface'] = water_surface
        self.fetches_params['classes'] = classes
        self.fetches_params['confidence_levels'] = confidence_levels
        self.fetches_params['columns'] = columns
        self.fetches_params['classify_bathymetry'] = classify_bathymetry
        self.fetches_params['classify_buildings'] = classify_buildings
        self.fetches_params['classify_water'] = classify_water
        self.fetches_params['reject_failed_qa'] = reject_failed_qa        
    
    def yield_ds(self, result):
        icesat2_fn= os.path.join(self.fetch_module._outdir, result['dst_fn'])
        if self.fetches_params['classify_buildings']:
            self.fetches_params['classify_buildings'] = self.process_buildings(self.fetch_buildings(verbose=True))

        if self.fetches_params['classify_water']:
            self.fetches_params['classify_water'] = self.process_coastline(
                self.fetch_coastline(chunks=False, verbose=False), return_geom=True, verbose=False
            )
            
        if 'processed_zip' in result['data_type']:
            icesat2_h5s = utils.p_unzip(
                icesat2_fn,
                exts=['h5'],
                outdir=self.cache_dir,
                verbose=self.verbose
            )
            for icesat2_h5 in icesat2_h5s:
                self.fetches_params['fn'] = icesat2_h5
                yield(IceSat2File(**self.fetches_params)._acquire_module())
                yield(sub_ds)
            
        else:
            self.fetches_params['fn'] = icesat2_fn
            yield(IceSat2File(**self.fetches_params)._acquire_module())        
                
class GMRTFetcher(Fetcher):
    """GMRT Gridded data.

    -----------
    Parameters:
    
    swath_only: only return MB swath data
    """
    
    __doc__ = '''{}    
    Fetches Module: <gmrt> - {}'''.format(__doc__, fetches.GMRT.__doc__)
    
    def __init__(self, swath_only = False, **kwargs):
        super().__init__(**kwargs)
        self.swath_only = swath_only

    def yield_ds(self, result):
        swath_mask=None
        gmrt_fn = os.path.join(self.fetch_module._outdir, result['dst_fn'])
        with gdalfun.gdal_datasource(gmrt_fn, update = 1) as src_ds:
            md = src_ds.GetMetadata()
            md['AREA_OR_POINT'] = 'Point'
            src_ds.SetMetadata(md)
            gdalfun.gdal_set_srs(src_ds)
            gdalfun.gdal_set_ndv(src_ds, verbose=False)

        if self.swath_only:
            if fetches.Fetch(
                    self.fetch_module._gmrt_swath_poly_url, verbose=self.verbose
            ).fetch_file(os.path.join(self.fetch_module._outdir, 'gmrt_swath_polygons.zip')) == 0:
                swath_shps = utils.p_unzip(
                    os.path.join(self.fetch_module._outdir, 'gmrt_swath_polygons.zip'),
                    exts=['shp', 'shx', 'prj', 'dbf'],
                    outdir=self.cache_dir,
                    verbose=self.verbose
                )
                for v in swath_shps:
                    if '.shp' in v:
                        #swath_shp = v
                        swath_mask = {'mask': v, 'invert_mask': True}
                        break
                    
                if not os.path.exists(swath_mask['mask']):
                    utils.echo_error_msg('could not find gmrt swath polygons...')
                    self.swath_only = False
                    swath_mask = None
                else:
                    self.fetches_params['mask'] = swath_mask

        yield(DatasetFactory(**self.fetches_params)._acquire_module())
        
class GEBCOFetcher(Fetcher):
    """GEBCO Gridded data

    Note: Fetches entire zip file.
    """
    
    __doc__ = '''{}    
    Fetches Module: <gebco> - {}'''.format(__doc__, fetches.GEBCO.__doc__)
        
    def __init__(self, exclude_tid = None, **kwargs):
        super().__init__(**kwargs)
        self.exclude_tid = []
        if utils.str_or(exclude_tid) is not None:
            for tid_key in exclude_tid.split('/'):
                self.exclude_tid.append(utils.int_or(tid_key))
        
    def yield_ds(self, result):
        wanted_gebco_fns = []
        gebco_fns = utils.p_unzip(
            os.path.join(self.fetch_module._outdir, result['dst_fn']),
            exts=['tif'],
            outdir=self.fetch_module._outdir,
            verbose=self.verbose
        )
        
        ## fetch the TID zip if needed
        if self.exclude_tid:
            if fetches.Fetch(
                    self.fetch_module._gebco_urls['gebco_tid']['geotiff'], verbose=self.verbose
            ).fetch_file(os.path.join(self.fetch_module._outdir, 'gebco_tid.zip')) == 0:
                tid_fns = utils.p_unzip(
                    os.path.join(self.fetch_module._outdir, 'gebco_tid.zip'),
                    exts=['tif'],
                    outdir=self.cache_dir,
                    verbose=self.verbose
                )
                for tid_fn in tid_fns:
                    ds_config = gdalfun.gdal_infos(tid_fn)                    
                    if self.region is not None and self.region.valid_p(check_xy=True):
                        inf_region = regions.Region().from_geo_transform(
                            ds_config['geoT'], ds_config['nx'], ds_config['ny']
                        )
                        inf_region.wmin = self.weight
                        inf_region.wmax = self.weight
                        inf_region.umin = self.uncertainty
                        inf_region.umax = self.uncertainty
                        if regions.regions_intersect_p(inf_region, self.region):
                            wanted_gebco_fns.append(tid_fn)
                    else:
                        wanted_gebco_fns.append(tid_fn)

                for tid_fn in wanted_gebco_fns:
                    tmp_tid = utils.make_temp_fn(
                        'tmp_tid.tif', temp_dir=self.fetch_module._outdir
                    )
                    with gdalfun.gdal_datasource(tid_fn) as tid_ds:
                        tid_config = gdalfun.gdal_infos(tid_ds)
                        tid_band = tid_ds.GetRasterBand(1)
                        tid_array = tid_band.ReadAsArray().astype(float)

                    tid_config['ndv'] = -9999
                    tid_config['dt'] = gdal.GDT_Float32                        
                    for tid_key in self.exclude_tid:
                        tid_array[tid_array == tid_key] = tid_config['ndv']
                        
                    for tid_key in self.fetch_module.tid_dic.keys():
                        tid_array[tid_array == tid_key] = self.fetch_module.tid_dic[tid_key][1]
                            
                    gdalfun.gdal_write(tid_array, tmp_tid, tid_config)
                    #utils.echo_msg(tmp_tid)
                    #utils.echo_warning_msg('mask: {}'.format(self.mask))
                    if self.mask is not None:
                        new_mask = utils.make_temp_fn('test_tmp_mask')
                        gdalfun.gdal_mask(
                            tmp_tid, self.mask['mask'], new_mask,
                            msk_value=1, verbose=True
                        )
                        os.replace(new_mask, tmp_tid)

                    self.fetches_params['mod'] = tid_fn.replace('tid_', '')
                    self.fetches_params['data_format'] = 200
                    self.fetches_params['mask'] = tmp_tid
                    self.fetches_params['weight_mask'] = tmp_tid
                    yield(DatasetFactory(**self.fetches_params)._acquire_module())
                    utils.remove_glob(tmp_tid)
        else:
            for gebco_fn in gebco_fns:
                ds_config = gdalfun.gdal_infos(gebco_fn)
                inf_region = regions.Region().from_geo_transform(
                    ds_config['geoT'], ds_config['nx'], ds_config['ny']
                )
                if self.region is not None and self.region.valid_p(check_xy=True):                
                    inf_region.wmin = self.weight
                    inf_region.wmax = self.weight
                    inf_region.umin = self.uncertainty
                    inf_region.umax = self.uncertainty

                    if regions.regions_intersect_p(inf_region, self.region):
                        wanted_gebco_fns.append(gebco_fn)
                else:
                    wanted_gebco_fns.append(gebco_fn)    
            
            for gebco_fn in wanted_gebco_fns:
                self.fetches_params['mod'] = gebco_fn
                self.fetches_params['data_format'] = 200
                yield(DatasetFactory(**self.fetches_params)._acquire_module())
                
class CopernicusFetcher(Fetcher):
    """Gridded Copernicus sattelite data.
    """
    
    __doc__ = '''{}    
    Fetches Module: <copernicus> - {}'''.format(__doc__, fetches.CopernicusDEM.__doc__)
    
    def __init__(self, datatype=None, **kwargs):
        super().__init__(**kwargs)
        self.check_size=False
        self.datatype=datatype
        
    def yield_ds(self, result):
        if self.datatype is None or result['data_type'] == self.datatype:
            src_cop_dems = utils.p_unzip(
                os.path.join(self.fetch_module._outdir, result['dst_fn']),
                exts=['tif'],
                outdir=self.fetch_module._outdir,
                verbose=self.verbose
            )
            for src_cop_dem in src_cop_dems:
                gdalfun.gdal_set_ndv(src_cop_dem, ndv=0, verbose=False)
                self.fetches_params['mod'] = src_cop_dem
                self.fetches_params['data_format'] = 200
                self.fetches_params['node'] = 'pixel'
                yield(DatasetFactory(**self.fetches_params)._acquire_module())

class FABDEMFetcher(Fetcher):
    """FABDEM Gridded data
    """
    
    __doc__ = '''{}
    Fetches Module: <fabdem> - {}'''.format(__doc__, fetches.FABDEM.__doc__)

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def yield_ds(self, result):
        src_fab_dems = utils.p_unzip(
            os.path.join(self.fetch_module._outdir, result['dst_fn']),
            exts=['tif'],
            outdir=self.fetch_module._outdir,
            verbose=self.verbose
        )
        for src_fab_dem in src_fab_dems:
            gdalfun.gdal_set_ndv(src_fab_dem, ndv=0, verbose=False)
            self.fetches_params['mod'] = src_fab_dem
            self.fetches_params['data_format'] = 200
            yield(DatasetFactory(**self.fetches_params)._acquire_module())

class MarGravFetcher(Fetcher):
    """Marine Gravity Bathymetry
    """
    
    __doc__ = '''{}    
    Fetches Module: <mar_grav> - {}'''.format(__doc__, fetches.MarGrav.__doc__)
    
    def __init__(self,
                 rasterize = False,
                 bathy_only = False,
                 upper_limit = None,
                 lower_limit = None,
                 **kwargs):
        super().__init__(**kwargs)
        self.rasterize = rasterize
        self.bathy_only = bathy_only
        self.upper_limit = utils.float_or(upper_limit)
        self.lower_limit = utils.float_or(lower_limit)
        if self.bathy_only:
            self.upper_limit = 0

        self.region.zmax=self.upper_limit
        self.region.zmin=self.lower_limit
            
    def yield_ds(self, result):
        if result['data_type'] == 'mar_grav_img':
            nc_fn = utils.make_temp_fn(
                '{}.nc'.format(utils.fn_basename2(result['dst_fn'])),
                temp_dir=self.fetch_module._outdir)
            img2grd_cmd = 'gmt img2grd {} {} -G{} -D -T1 -I1m -E'.format(
                os.path.join(
                    self.fetch_module._outdir,
                    result['dst_fn']
                ), self.region.format('gmt'), nc_fn
            )
            out, status = utils.run_cmd(img2grd_cmd, verbose=True)
            out, status = utils.run_cmd('gmt grdedit {} -T'.format(nc_fn))
            self.fetches_params['mod'] = nc_fn
            self.fetches_params['data_format'] = 200
            self.fetches_params['resample_and_warp'] = False
            self.fetches_params['node'] = 'grid'
            
        elif self.rasterize:
            from cudem import waffles
            mg_region = self.region.copy()
            mg_region.zmax = self.upper_limit
            mg_region.zmin = self.lower_limit        
            mar_grav_fn = utils.make_temp_fn('mar_grav')
            _raster = waffles.WaffleFactory(
                mod='IDW:min_points=16',
                data=['{},168:x_offset=REM,1'.format(
                    os.path.join(self.fetch_module._outdir, result['dst_fn'])
                )],
                src_region=mg_region,
                xinc=utils.str2inc('30s'),
                yinc=utils.str2inc('30s'),
                upper_limit = self.upper_limit,
                name=mar_grav_fn,
                node='pixel',
                verbose=True
            )._acquire_module()()            
            if self.upper_limit is not None or self.lower_limit is not None:
                ds = gdal.Open(_raster.fn)
                ds_band = ds.GetRasterBand(1)
                ds_arr = ds_band.ReadAsArray()
                if self.upper_limit is not None:
                    ds_arr[ds_arr >= self.upper_limit] = ds_band.GetNoDataValue()

                if self.lower_limit is not None:
                    ds_arr[ds_arr <= self.lower_limit] = ds_band.GetNoDataValue()
                    
                ds = None

            self.fetches_params['mod'] = _raster.fn
            self.fetches_params['data_format'] = 200
            
        yield(DatasetFactory(**self.fetches_params)._acquire_module())

class ChartsFetcher(Fetcher):
    """NOAA ENC Charts Fetcher

    Digital Soundings
    """
    
    __doc__ = '''{}
    Fetches Module: <charts> - {}'''.format(__doc__, fetches.Charts.__doc__)

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        
    def yield_ds(self, result):
        src_000s = utils.p_unzip(
            os.path.join(self.fetch_module._outdir, result['dst_fn']),
            exts=['000'],
            outdir=self.fetch_module._outdir,
            verbose=self.verbose            
        )
        for src_000 in src_000s:
            self.fetches_params['mod'] = src_000
            self.fetches_params['data_format'] = 302
            self.fetches_params['ogr_layer'] = 'SOUNDG'
            self.fetches_params['z_scale'] = -1
            yield(DatasetFactory(**self.fetches_params)._acquire_module())

class MBSFetcher(Fetcher):
    """NOAA Multibeam Fetcher
    """

    __doc__ = '''{}
    Fetches Module: <multibeam> - {}'''.format(__doc__, fetches.Multibeam.__doc__)

    def __init__(self, mb_exclude = 'A', want_binned = False, want_mbgrid = False, **kwargs):
        super().__init__(**kwargs)
        self.fetches_params['mb_exclude'] = mb_exclude
        self.fetches_params['want_binned'] = want_binned
        self.fetches_params['want_mbgrid'] = want_mbgrid
        self.fetches_parasm['want_inf'] = False

    def yield_ds(self, result):            
        mb_infos = self.fetch_module.parse_entry_inf(result, keep_inf=True)
        yield(DatasetFactory(**self.fetches_params)._acquire_module())

class HydroNOSParser(ElevationDataset):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
                
class HydroNOSFetcher(Fetcher):
    """NOAA HydroNOS Data Fetcher
    """
    
    __doc__ = '''{}
    Fetches Module: <hydronos> - {}'''.format(__doc__, fetches.HydroNOS.__doc__)

    def __init__(self, explode=False, **kwargs):
        super().__init__(**kwargs)
        self.explode=explode

    def yield_ds(self, result):
        if result['data_type'] == 'xyz':
            nos_fns = utils.p_unzip(
                os.path.join(self.fetch_module._outdir, result['dst_fn']),
                exts=['xyz', 'dat'],
                outdir=os.path.dirname(os.path.join(self.fetch_module._outdir, result['dst_fn'])),
                verbose=self.verbose
            )
            for nos_fn in nos_fns:
                self.fetches_params['mod'] = nos_fn
                self.fetches_params['data_format'] = '168:skip=1:xpos=2:ypos=1:zpos=3:z_scale=-1:delimiter=,'
                self.fetches_params['src_srs'] = 'epsg:4326+5866'
                yield(DatasetFactory(**self.fetches_params)._acquire_module())

        elif result['data_type'] == 'bag':
            bag_fns = utils.p_unzip(
                os.path.join(self.fetch_module._outdir, result['dst_fn']),
                exts=['bag'],
                outdir=os.path.dirname(os.path.join(self.fetch_module._outdir, result['dst_fn'])),
                verbose=self.verbose
            )
            for bag_fn in bag_fns:
                # bag_fn = os.path.join(self.fetch_module._outdir, result[1])
                if 'ellipsoid' not in bag_fn.lower():
                    self.fetches_params['mod'] = bag_fn
                    self.fetches_params['data_format'] = 201
                    self.fetches_params['src_srs'] = None
                    yield(DatasetFactory(**self.fetches_params)._acquire_module())

class CSBFetcher(Fetcher):
    """Crowd Sourced Bathymetry data fetcher
    """
    
    __doc__ = '''{}
    Fetches Module: <csb> - {}'''.format(__doc__, fetches.CSB.__doc__)

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def yield_ds(self, result):
        yield(DatasetFactory(**self.fetches_params)._acquire_module())

class EMODNetFetcher(Fetcher):
    """EMODNet Data Fetcher
    """
    
    __doc__ = '''{}
    Fetches Module: <emodnet> - {}'''.format(__doc__, fetches.EMODNet.__doc__)

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def yield_ds(self, result):
        if result['data_type'] == 'csv':
            self.fetches_params['data_format'] = '168:skip=1:xpos=2:ypos=1:zpos=3:delimiter=,'
        elif result['data_type'] == 'nc':
            self.fetches_params['data_format'] = 200
            
        yield(DatasetFactory(**self.fetches_params)._acquire_module())
        
class GEDTM30Fetcher(Fetcher):
    """GEDTM30 Data Fetcher
    """
    
    __doc__ = '''{}
    Fetches Module: <gedtm30> - {}'''.format(__doc__, fetches.GEDTM30.__doc__)

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def parse(self):
        for result in self.fetch_module.results:
            for this_ds in self.yield_ds(result):
                if this_ds is not None:
                    this_ds.remote = True
                    this_ds.initialize()
                    for ds in this_ds.parse():
                        yield(ds)
                                        
    def yield_ds(self, result):
        self.fetches_params['mod'] = '/vsicurl/{}'.format(result['url'])
        self.fetches_params['data_format'] = 200            
        yield(DatasetFactory(**self.fetches_params)._acquire_module())

class eHydroFetcher(Fetcher):
    """USACE eHydro soundings
    """

    __doc__ = '''{}
    Fetches Module: <ehydro> - {}'''.format(__doc__, fetches.eHydro.__doc__)
        
    def __init__(self, want_contours = False, **kwargs):
        super().__init__(**kwargs)
        self.want_contours = want_contours
        
    def yield_ds(self, result):
        try:
            src_gdb = utils.gdb_unzip(
                os.path.join(self.fetch_module._outdir, result['dst_fn']),
                outdir=self.fetch_module._outdir,
                verbose=False
            )
        except Exception as e:
            utils.echo_error_msg('{}: {}'.format(self.fn, e))
            src_gdb = None

        if src_gdb is not None:
            tmp_gdb = ogr.Open(src_gdb)
            tmp_layer = tmp_gdb.GetLayer('SurveyPoint')
            src_srs = tmp_layer.GetSpatialRef()
            src_epsg = gdalfun.osr_parse_srs(src_srs)
            elev_datum = tmp_layer[1].GetField('elevationDatum')
            tmp_gdb = tmp_layer = None
            v = vdatums.get_vdatum_by_name(elev_datum)
            self.fetches_params['mod'] = src_gdb
            self.fetches_params['src_srs'] = '{}+{}'.format(
                src_epsg, v if v is not None else '5866'
            ) if src_epsg is not None else None
            self.fetches_params['data_format'] = '302:ogr_layer=SurveyPoint_HD:elev_field=Z_label:z_scale=-0.3048006096012192'
            self.metadata['name'] = self.fn
            yield(DatasetFactory(**self.fetches_params)._acquire_module())            

            if self.want_contours:
                self.metadata['name'] = '{}_contours'.format(utils.fn_basename2(self.fn))
                #self.fetches_params['data_format'] = '302:ogr_layer=ElevationContour_ALL:elev_field=contourElevation:z_scale=-0.3048'
                self.fetches_params['data_format'] = '302:ogr_layer=ElevationContour_ALL:elev_field=contourElevation:z_scale=-0.3048006096012192'
                yield(DatasetFactory(**self.fetches_params)._acquire_module())            

    def yield_ds_XYZ(self, result):
        src_gdb = utils.gdb_unzip(
            os.path.join(self.fetch_module._outdir, result['dst_fn']),
            outdir=self.fetch_module._outdir,
            verbose=False
        )
        if src_gdb is not None:
            tmp_gdb = ogr.Open(src_gdb)
            tmp_layer = tmp_gdb.GetLayer('SurveyPoint')
            src_srs = tmp_layer.GetSpatialRef()
            src_epsg = gdalfun.osr_parse_srs(src_srs)
            tmp_gdb = None
            src_usaces = utils.p_unzip(
                os.path.join(self.fetch_module._outdir, result['dst_fn']),
                ['XYZ', 'xyz', 'dat'],
                outdir=self.fetch_module._outdir,
                verbose=self.verbose
            )
            for src_usace in src_usaces:
                self.fetches_params['mod'] = src_usace
                self.fetches_params['data_format'] = '168:z_scale=.3048'
                self.fetches_params['src_srs'] = '{}+{}'.format(src_epsg, v if v is not None else '5866') \
                    if src_epsg is not None \
                       else None
                yield(DatasetFactory(**self.fetches_params)._acquire_module())        
        
class BlueTopoFetcher(Fetcher):
    """BlueTopo Gridded bathymetric data Fetcher

    -----------
    Parameters:
    
    want_interpolation: True/False to include interpolated cells
    unc_weights: use the uncertainty mask as weights
    """

    __doc__ = '''{}
    Fetches Module: <bluetopo> - {}'''.format(__doc__, fetches.BlueTopo.__doc__)

    def __init__(self, want_interpolation=False, unc_weights=False, **kwargs):
        super().__init__(**kwargs)
        self.want_interpolation = want_interpolation
        self.unc_weights = unc_weights

    def yield_ds(self, result):
        sid = None
        if not self.want_interpolation:
            sid = gdalfun.gdal_extract_band(
                os.path.join(self.fetch_module._outdir, result['dst_fn']),
                utils.make_temp_fn('tmp_bt_tid.tif', self.fetch_module._outdir),
                band=3,
                exclude=[0]
            )[0]

        if self.mask is not None:
            new_mask = utils.make_temp_fn('test_tmp_mask')
            gdalfun.gdal_mask(sid, self.mask['mask'], new_mask, msk_value=1, verbose=True)
            os.replace(new_mask, sid)

        self.fetches_params['data_format'] = '200:band_no=1:mask={}:uncertainty_mask=2{}'.format(
            sid, ':weight_mask=2' if self.unc_weights else ''
        )
        yield(DatasetFactory(**self.fetches_params)._acquire_module())        

class NGSFetcher(Fetcher):
    """NGS Monument data
    """

    __doc__ = '''{}
    Fetches Module: <ngs> - {}'''.format(__doc__, fetches.NGS.__doc__)
    
    def __init__(self, datum = 'geoidHt', **kwargs):
        super().__init__(**kwargs)
        self.datum = datum
        if self.datum not in ['orthoHt', 'geoidHt', 'z', 'ellipHeight']:
            utils.echo_warning_msg('could not parse {}, falling back to geoidHt'.format(datum))
            self.datum = 'geoidHt'

    def yield_ds(self, result):
        with open(os.path.join(self.fetch_module._outdir, result['dst_fn']), 'r') as json_file:
            r = json.load(json_file)
            if len(r) > 0:
                with open(os.path.join(self.fetch_module._outdir, '_tmp_ngs.xyz'), 'w') as tmp_ngs:
                    for row in r:
                        z = utils.float_or(row[self.datum])
                        if z is not None:
                            xyz = xyzfun.XYZPoint(src_srs='epsg:4326').from_list(
                                [float(row['lon']), float(row['lat']), z]
                            )
                            xyz.dump(dst_port=tmp_ngs)

        self.fetches_params['mod'] = os.path.join(self.fetch_module._outdir, '_tmp_ngs.xyz')
        self.fetches_params['data_format'] = 168
        yield(DatasetFactory(**self.fetches_params)._acquire_module())

        utils.remove_glob(os.path.join(self.fetch_module._outdir, '_tmp_ngs.xyz'))

class TidesFetcher(Fetcher):
    """NOS Tide Station data
    """

    __doc__ = '''{}
    Fetches Module: <tides> - {}'''.format(__doc__, fetches.Tides.__doc__)
    
    def __init__(self, s_datum='mllw', t_datum='msl', units='m', **kwargs):
        super().__init__(**kwargs)
        self.s_datum = s_datum
        self.t_datum = t_datum
        self.units = units
        
    def yield_ds(self, result):
        with open(os.path.join(self.fetch_module._outdir, result['dst_fn']), 'r') as json_file:
            r = json.load(json_file)
            if len(r) > 0:
                with open(os.path.join(self.fetch_module._outdir, '_tmp_tides.xyz'), 'w') as tmp_ngs:
                    for feature in r['features']:
                        # if self.fetch_module.station_id is not None:
                        #     if self.fetch_module.station_id != feature['attributes']['id']:
                        #         continue
                            
                        lon = feature['attributes']['longitude']
                        lat = feature['attributes']['latitude']
                        if feature['attributes'][self.s_datum] != -99999.99 and feature['attributes'][self.t_datum] != -99999.99:
                            z = feature['attributes'][self.s_datum] - feature['attributes'][self.t_datum]
                            if self.units == 'm':
                                z = z * 0.3048

                            xyz = xyzfun.XYZPoint(src_srs='epsg:4326').from_list([lon, lat, z])
                            xyz.dump(dst_port=tmp_ngs)

        self.fetches_params['mod'] = os.path.join(self.fetch_module._outdir, '_tmp_tides.xyz')
        self.fetches_params['data_format'] = 168
        yield(DatasetFactory(**self.fetches_params)._acquire_module())

        utils.remove_glob(os.path.join(self.fetch_module._outdir, '_tmp_tides.xyz'))
        
class WaterServicesFetcher(Fetcher):
    """USGS Water Services

    -----------
    Parameters:
    
    site_code: the site code to fetch
    units: 'm' for meters

    site_codes:
    00065 - Gate Height
    00060 - StreamFlow
    63160 - Stream water level elevation above NAVD 1988
    62611 - Groundwater level above NAVD 1988
    72019 - Depth to water level, units below land surface
    """

    __doc__ = '''{}
    Fetches Module: <waterservices> - {}'''.format(__doc__, fetches.WaterServices.__doc__)
    
    def __init__(self, site_code='00065', units='m', **kwargs):
        super().__init__(**kwargs)
        self.units = units
        self.site_code = site_code
        
    def yield_ds(self, result):
        with open(os.path.join(self.fetch_module._outdir, result['dst_fn']), 'r') as json_file:
            r = json.load(json_file)
            if len(r) > 0:
                with open(os.path.join(self.fetch_module._outdir, '_tmp_ws.xyz'), 'w') as tmp_ws:
                    features = r['value']['timeSeries']
                    for feature in features:
                        if feature['variable']['variableCode'][0]['value'] == self.site_code:
                            lon = float(feature['sourceInfo']['geoLocation']['geogLocation']['longitude'])
                            lat = float(feature['sourceInfo']['geoLocation']['geogLocation']['latitude'])
                            z = float(feature['values'][0]['value'][0]['value'])

                            if self.units == 'm':
                                z = z * 0.3048
                            
                            xyz = xyzfun.XYZPoint(src_srs='epsg:4326').from_list([lon, lat, z])
                            xyz.dump(dst_port=tmp_ws)

        self.fetches_params['mod'] = os.path.join(self.fetch_module._outdir, '_tmp_ws.xyz')
        self.fetches_params['data_format'] = 168
        yield(DatasetFactory(**self.fetches_params)._acquire_module())

        utils.remove_glob(os.path.join(self.fetch_module._outdir, '_tmp_ws.xyz'))
                            
class VDatumFetcher(Fetcher):
    """VDatum transformation grids.
    """

    __doc__ = '''{}
    Fetches Module: <vdatum> - {}'''.format(__doc__, fetches.VDATUM.__doc__)
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def yield_ds(self, result):
        src_tif = os.path.join(
            self.fetch_module._outdir, '{}.tif'.format(
                utils.fn_basename2(os.path.basename(result['dst_fn']))
            )
        )
        if result['dst_fn'].endswith('.zip'):
            v_gtx = utils.p_f_unzip(
                os.path.join(self.fetch_module._outdir, result['dst_fn']), [result['data_type']], outdir=self.fetch_module._outdir
            )[0]
            utils.run_cmd(
                'gdalwarp {} {} -t_srs epsg:4269 --config CENTER_LONG 0'.format(v_gtx, src_tif),
                verbose=self.verbose
            )

        self.fetches_params['mod'] = src_tif
        self.fetches_params['data_format'] = 200
        self.fetches_params['node'] = 'pixel'
        yield(DatasetFactory(**self.fetches_params)._acquire_module())
        
## todo: allow lakes bathymetry
## as well as lakes breaklines (shape nodes)
## see: https://www.esri.com/arcgis-blog/products/arcgis-pro/3d-gis/hydro-flattening-of-river-shorelines-in-lidar-based-dem-production/
class HydroLakesFetcher(Fetcher):
    """HydroLakes lake bathymetric data
    """
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def yield_ds(self, result):
        pass
    
class DatasetFactory(factory.CUDEMFactory):
    """Dataset Factory Settings and Generator
    
    Parse a datalist entry and return the dataset object
    """
    
    _modules = {
        ## negative values are `datalists`, where they contain various other datasets.
        -1: {'name': 'datalist',
             'fmts': ['datalist', 'mb-1', 'dl'],
             'description': 'An extended MB-System style datalist containting dlim-compatible datasets',
             'call': Datalist},
        -2: {'name': 'zip',
             'fmts': ['zip', 'ZIP'],
             'description': 'A zipfile containing dlim-compatible datasets',
             'call': ZIPlist}, # add other archive formats (gz, tar.gz, 7z, etc.)
        -3: {'name': 'scratch',
             'fmts': [''],
             'description': 'A scratch dataset, including a python list of dlim-compatible datasets',
             'call': Scratch},
        -4: {'name': 'points',
             'fmts': [''],
             'description': 'A points dataset, a numpy rec-array or a pandas dataframe, defining x, y and z data columns',
             'call': Points},
        ## data files
        167: {'name': 'yxz',
              'fmts': ['yxz'],
              'description': 'ascii DSV datafile formatted as y,x,z',
              'call': YXZFile},
        168: {'name': 'xyz',
              'fmts': ['xyz', 'csv', 'dat', 'ascii', 'txt', 'XYZ'],
              'description': 'An ascii DSV datafile formatted as x,y,z',
              'call': XYZFile},
        200: {'name': 'gdal',
              'fmts': ['tif', 'tiff', 'img', 'grd', 'nc', 'vrt'],
              'description': 'A gdal-compatible raster dataset',
              'call': GDALFile},
        201: {'name': 'bag',
              'fmts': ['bag'],
              'description': 'A BAG bathymetry dataset',
              'call': BAGFile},
        202: {'name': 'swot_pixc',
              'fmts': ['h5'],
              'description': 'An HDF5 SWOT PIXC datafile',
              'call': SWOT_PIXC},
        203: {'name': 'swot_hr_raster',
              'fmts': ['nc'],
              'description': 'An HDF5 SWOT HR Raster datafile',
              'call': SWOT_HR_Raster},
        300: {'name': 'las',
              'fmts': ['las', 'laz'],
              'description': 'An las or laz lidar datafile',
              'call': LASFile},
        301: {'name': 'mbs',
              'fmts': ['fbt', 'mb'],
              'description': 'An MB-System-compatible multibeam datafile',
              'call': MBSParser},
        302: {'name': 'ogr',
              'fmts': ['000', 'shp', 'geojson', 'gpkg', 'gdb/'],
              'description': 'An ogr-compatible vector datafile',
              'call': OGRFile},
        303: {'name': 'icesat2_atl',
              'fmts': ['h5'],
              'description': 'An HDF5 IceSat2 ATL03 datafile',
              'call': IceSat2File},
        304: {'name': 'cudem',
              'fmts': ['csg', 'nc', 'h5'],
              'description': 'A netCDF/h5 CUDEM file',
              'call': CUDEMFile},
        ## fetches modules
        -100: {'name': 'https',
               'fmts': ['https'],
               'description': 'A URL pointing to a dlim-compatible datafile',
               'call': Fetcher},
        -101: {'name': 'gmrt',
               'fmts': ['gmrt'],
               'description': 'The GMRT fetches module',
               'call': GMRTFetcher},
        -102: {'name': 'gebco',
               'fmts': ['gebco'],
               'description': 'The GEBCO fetches module',
               'call': GEBCOFetcher},
        -103: {'name': 'copernicus',
               'fmts': ['copernicus'],
               'description': 'The Copernicus fetches module',
               'call': CopernicusFetcher},
        -104: {'name': 'fabdem',
               'fmts': ['fabdem'],
               'description': 'The FABDEM fetches module',
               'call': FABDEMFetcher},
        -105: {'name': 'nasadem',
               'fmts': ['nasadem'],
               'description': 'The NASADEM fetches module',
               'call': Fetcher},
        -106: {'name': 'mar_grav',
               'fmts': ['mar_grav'],
               'description': 'The mar_grav fetches module',
               'call': MarGravFetcher},
        -107: {'name': 'srtm_plus',
               'fmts': ['srtm_plus'],
               'description': 'The srtm_plus fetches module',
               'call': Fetcher},
        -108: {'name': 'synbath',
               'fmts': ['synbath'],
               'description': 'The SynBath fetches module',
               'call': Fetcher},
        -109: {'name': 'gedtm30',
               'fmts': ['gedtm30'],
               'description': 'Global DTM fetches module',
               'call': GEDTM30Fetcher},
        -110: {'name': "swot",
               'fmts': ['swot'],
               'description': '	The SWOT fetches module',
               'call': SWOTFetcher},
        -111: {'name': "icesat2",
               'fmts': ['icesat2'],
               'description': 'The IceSat2 fetches module',
               'call': IceSat2Fetcher},
        -200: {'name': 'charts',
               'fmts': ['charts'],
               'description': 'The charts fetches module',
               'call': ChartsFetcher},
        -201: {'name': 'multibeam',
               'fmts': ['multibeam'],
               'description': 'The multibeam fetches module',
               'call': MBSFetcher},
        -202: {'name': 'hydronos',
               'fmts': ['hydronos'],
               'description': 'The hydronos fetches module',
               'call': HydroNOSFetcher},
        -203: {'name': 'ehydro',
               'fmts': ['ehydro'],
               'description': 'The ehydro fetches module',
               'call': eHydroFetcher},
        -204: {'name': 'bluetopo',
               'fmts': ['bluetopo'],
               'description': 'The bluetopo fetches module',
               'call': BlueTopoFetcher},
        -205: {'name': 'ngs',
               'fmts': ['ngs'],
               'description': 'The ngs fetches module',
               'call': NGSFetcher},
        -206: {'name': 'tides',
               'fmts': ['tides'],
               'description': 'The tides fetches module',
               'call': TidesFetcher},
        -207: {'name': 'digital_coast',
               'fmts': ['digital_coast'],
               'description': 'The digital_coast fetches module',
               'call': Fetcher},
        -208: {'name': 'ncei_thredds',
               'fmts': ['ncei_thredds'],
               'description': 'The ncei_thredds fetches module',
               'call': Fetcher},
        -209: {'name': 'tnm',
               'fmts': ['tnm'],
               'description': 'The TNM fetches module',
               'call': Fetcher},
        -210: {'name': "CUDEM",
               'fmts': ['CUDEM'],
               'description': 'The CUDEM fetches module',
               'call': Fetcher},
        -211: {'name': "CoNED",
               'fmts': ['CoNED'],
               'description': 'The CoNED fetches module',
               'call': DAVFetcher_CoNED},
        -212: {'name': "SLR",
               'fmts': ['SLR'],
               'description': '	The SLR fetches module',
               'call': DAVFetcher_SLR},
        -213: {'name': 'waterservies',
               'fmts': ['waterservices'],
               'description': 'The waterservices fetches module',
               'call': WaterServicesFetcher},
        -214: {'name': 'ned',
               'fmts': ['ned', 'ned1'],
               'description': 'The NED fetches module',
               'call': NEDFetcher},        
        -215: {'name': "csb",
               'fmts': ['csb'],
               'description': 'The CSB fetches module',
               'call': Fetcher},
        -216: {'name': 'wa_dnr',
               'fmts': ['wa_dnr'],
               'description': 'The Washington DNR lidar portal',
               'call': DNRFetcher},
        -217: {'name': 'r2r',
               'fmts': ['r2r'],
               'description': 'The r2r fetches module',
               'call': MBSFetcher},
        -300: {'name': 'emodnet',
               'fmts': ['emodnet'],
               'description': 'The emodnet fetches module',
               'call': EMODNetFetcher},
        -301: {'name': 'chs',
               'fmts': ['chs'],
               'description': 'The chs fetches module',
               'call': Fetcher}, # chs is broken
        -302: {'name': 'hrdem',
               'fmts': ['hrdem'],
               'description': '	The hrdem fetches module',
               'call': Fetcher},
        -303: {'name': 'arcticdem',
               'fmts': ['arcticdem'],
               'description': 'The arcticdem fetches module',
               'call': Fetcher},
        -500: {'name': 'vdatum',
               'fmts': ['vdatum'],
               'description': 'The vdatum fetches module',
               'call': VDatumFetcher},
        -600: {'name': 'hydrolakes',
               'fmts': ['hydrolakes'],
               'description': 'The hydrolakes fetches module',
               'call': HydroLakesFetcher},        
    }
    _datalist_cols = ['path', 'format', 'weight', 'uncertainty', 'title', 'source',
                      'date', 'type', 'resolution', 'horz', 'vert',
                      'url']

    _metadata_keys = ['name', 'title', 'source', 'date', 'data_type', 'resolution',
                      'hdatum', 'vdatum', 'url']
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    ## redefine the factory default _parse_mod function for datasets
    ## the mod in this case is a datalist entry and the format key
    ## becomes the module
    ## TODO: use csv module to parse
    def _parse_mod(self, mod=None):
        """parse the datalist entry line"""
        
        self.kwargs['fn'] = mod
        if self.kwargs['fn'] is None:
            return(self)

        ## mod exists as a file, no other entry items should occur, so
        ## guess the format and finish there...
        ## the format number becomes the mod_name
        ## check for specified data format as well
        # if os.path.exists(self.kwargs['fn']):
        #     if 'data_format' not in self.kwargs.keys() or self.kwargs['data_format'] is None:
        #         self.mod_name = self.guess_data_format(self.kwargs['fn'])
        #         self.mod_args = {}
        #         self.kwargs['data_format'] = self.mod_name
        #     else:
        #         opts = str(self.kwargs['data_format']).split(':')
        #         if len(opts) > 1:
        #             self.mod_name = int(opts[0])
        #             self.mod_args = utils.args2dict(list(opts[1:]), {})
        #         else:
        #             self.mod_name = int(self.kwargs['data_format'])
        #             self.mod_args = {}
                    
        #         self.kwargs['data_format'] = self.mod_name

        #     # inherit metadata from parent if available
        #     # something something!
        #     if 'metadata' not in self.kwargs.keys():
        #         self.kwargs['metadata'] = {}
                
        #     self.kwargs['metadata']['name'] = utils.fn_basename2(os.path.basename(self.kwargs['fn']))
            
        #     for key in self._metadata_keys:
        #         if key not in self.kwargs['metadata'].keys():
        #             self.kwargs['metadata'][key] = None
                    
        #     return(self.mod_name, self.mod_args)

        ## if fn is not a path, parse it as a datalist entry
		## breaks on path with space, e.g. /meda/user/My\ Passport/etc
        #this_entry = re.findall(r'[^"\s]\S*|".+?"', self.kwargs['fn'].rstrip())
        #this_entry = shlex.split(shlex.quote(self.kwargs['fn'].rstrip()))#.replace("'", "\'")))
        #this_entry = [p for p in re.split("( |\\\".*?\\\"|'.*?')", self.kwargs['fn']) if p.strip()]
        #this_entry = [t.strip('"') for t in re.findall(r'[^\s"]+|"[^"]*"', self.kwargs['fn'].rstrip())]
        this_entry = re.findall(r"(?:\".*?\"|\S)+", self.kwargs['fn'].rstrip())
        try:
            entry = [utils.str_or(x) if n == 0 \
                     else utils.str_or(x) if n < 2 \
                     else utils.float_or(x) if n < 3 \
                     else utils.float_or(x) if n < 4 \
                     else utils.str_or(x) \
                     for n, x in enumerate(this_entry)]

        except Exception as e:
            utils.echo_error_msg('could not parse entry {}, {}'.format(self.kwargs['fn'], this_entry))
            return(self)

        ## data format - entry[1]
        ## guess format based on fn if not specified otherwise
        ## parse the format for dataset specific opts.
        if len(entry) < 2:
            if 'data_format' in self.kwargs.keys() and self.kwargs['data_format'] is not None:
                entry.append(self.kwargs['data_format'])
            else:
                for key in self._modules.keys():
                    if entry[0].startswith('http') or entry[0].startswith('/vsicurl/'):
                        see = 'https'
                    else:
                        se = entry[0].split('.')
                        see = se[-1] if len(se) > 1 else entry[0].split(":")[0]

                    if 'fmts' in self._modules[key].keys():
                        if see in self._modules[key]['fmts']:
                            entry.append(int(key))
                            break

            if len(entry) < 2:
                utils.echo_error_msg('could not parse entry {}'.format(self.kwargs['fn']))
                return(self)
            
        ## parse the entry format options
        opts = str(entry[1]).split(':')
        if len(opts) > 1:
            self.mod_args = utils.args2dict(list(opts[1:]), {})
            entry[1] = int(opts[0])
        else:
            self.mod_args = {}
                                                       
        try:
            assert isinstance(utils.int_or(entry[1]), int)
        except:
            utils.echo_error_msg('could not parse datalist entry {}'.format(entry))
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
            if self.mod_name >= -2:
                self.kwargs['fn'] = os.path.join(
                    os.path.dirname(self.kwargs['parent'].fn), entry[0]
                )
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
                self.kwargs['uncertainty'] = math.sqrt(self.kwargs['uncertainty']**2 + entry[3]**2)
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
                    self.kwargs['metadata'][key] = utils.fn_basename2(os.path.basename(self.kwargs['fn'].split(':')[0]))
                else:
                    self.set_metadata_entry(utils.fn_basename2(os.path.basename(self.kwargs['fn'].split(':')[0])), key, '/')
                    #self.kwargs['metadata'][key] = '/'.join([self.kwargs['metadata'][key], utils.fn_basename2(os.path.basename(self.kwargs['fn']))])

            else:
                if len(entry) < i+4:
                    entry.append(self.kwargs['metadata'][key])

                self.set_metadata_entry(entry[i+3], key, ', ')

                if key == 'date':
                    self.kwargs['metadata'][key] = num_strings_to_range(self.kwargs['metadata'][key], entry[i+3])
                
        return(self.mod_name, self.mod_args)

    def set_metadata_entry(self, entry, metadata_field, join_string = '/'):
        if entry != '-' and str(entry).lower() != 'none':
            if self.kwargs['metadata'][metadata_field] is not None:
                if str(self.kwargs['metadata'][metadata_field]).lower() != str(entry).lower():
                    self.kwargs['metadata'][metadata_field] = join_string.join([str(self.kwargs['metadata'][metadata_field]), str(entry)])
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
            
    def write_parameter_file(self, param_file: str):
        try:
            with open(param_file, 'w') as outfile:
                json.dump(self.__dict__, outfile)
                utils.echo_msg('New DatasetFactory file written to {}'.format(param_file))
                
        except:
            raise ValueError('DatasetFactory: Unable to write new parameter file to {}'.format(param_file))            

## ==============================================
## Command-line Interface (CLI)
## $ dlim
##
## datalists cli
## ==============================================
datalists_usage = lambda: """{cmd} ({dl_version}): DataLists IMproved; Process and generate datalists

dlim is the elevation data processing tool using various dataset modules. dlim's native dataset format is a "datalist". 
A datalist is similar to an MBSystem datalist; it is a space-delineated file containing the following columns:

`data-path data-format data-weight data-uncertainty data-name data-source data-date data-resolution data-type data-horz data-vert data-url`

Minimally, `data-path` (column 1) is all that is needed.

An associated inf and geojson file will be gerenated for each datalist while only an associated inf file will be genereated 
for individual datasets

Parse various dataset types by region/increments and yield data as xyz or array recursive data-structures which point to 
datasets (datalist, zip, fetches, etc) are negative format numbers, e.g. -1 for datalist

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
\t\t\t\te.g. -X6:10 to extend the DEM REGION by 6 cells and the processing region by 10 
\t\t\t\tpercent of the input REGION.
  -J, --s_srs\t\t\tSet the SOURCE projection.
  -P, --t_srs\t\t\tSet the TARGET projection. (REGION should be in target projection) 
  -D, --cache-dir\t\tCACHE Directory for storing temp and output data.
  -Z, --z-precision\t\tSet the target precision of dumped z values. (default is 4)
  -A, --stack-mode\t\tSet the STACK MODE to 'mean', 'min', 'max', 'mixed' or 'supercede' (with -E and -R)
  -T, --stack_filter\t\tFILTER the data stack using one or multiple filters. 
\t\t\t\tWhere FILTER is fltr_name[:opts] (see `grits --modules` for more information)
\t\t\t\tThe -T switch may be set multiple times to perform multiple filters.
\t\t\t\tAvailable FILTERS: {grits_modules}
  -F, --point_filter\t\tFILTER the POINT data using one or multiple filters. 
\t\t\t\tWhere FILTER is fltr_name[:opts] 
\t\t\t\tThe -F switch may be set multiple times to perform multiple filters.
\t\t\t\tAvailable FILTERS: {point_filter_modules}
  -V, --archive\t\t\tArchive the DATALIST to the given REGION[/INCREMENTs].
\t\t\t\tSpecify the name of the archive, if not specified an auto-generated name will be used.

  -m, --mask\t\t\tMASK the datalist to the given REGION/INCREMENTs
  -s, --spatial-metadata\tGenerate SPATIAL METADATA of the datalist to the given REGION/INCREMENTs
  -g, --glob\t\t\tGLOB the datasets in the current directory to stdout
  -i, --info\t\t\tGenerate and return an INFO dictionary of the dataset
  -l, --list\t\t\tList the assocated datasets from the datalist
  -w, --weights\t\t\tOutput WEIGHT values along with xyz
  -u, --uncertainties\t\tOutput UNCERTAINTY values along with xyz
  -n, --stack-node\t\tOutput stacked x/y data rather than pixel
  -q, --quiet\t\t\tLower the verbosity to a quiet

  --modules\t\t\tDisplay the datatype descriptions and usage
  --help\t\t\tPrint the usage text
  --version\t\t\tPrint the version information

Supported datalist formats (see {cmd} --modules <dataset-key> for more info): 
  {dl_formats}

Examples:

  % {cmd} -R-90/-89/30/31/-100/100 *.tif -l -w > tifs_in_region.datalist
  % {cmd} tifs_in_region.datalist -R -90/-89/30/31 -E 1s > tifs_1s.xyz
  % {cmd} -R my_region.shp my_data.xyz -w -s_srs epsg:4326 -t_srs epsg:3565 > my_data_3565.xyz
""".format(cmd=os.path.basename(sys.argv[0]), 
           dl_version=cudem.__version__,
           dl_formats=factory._cudem_module_name_short_desc(DatasetFactory._modules),
           grits_modules=factory._cudem_module_short_desc(grits.GritsFactory._modules),
           point_filter_modules=factory._cudem_module_short_desc(PointFilterFactory._modules))

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
        elif arg == '-t_srs' or arg == '--t_srs' or arg == '-P':
            dst_srs = argv[i + 1]
            i = i + 1
        elif arg == '-z_precision' or arg == '--z-precision' or arg == '-Z':
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
            dataexts = [xs for y in [x['fmts'] for x in list(DatasetFactory._modules.values())] for xs in y]
            #utils.echo_msg_bold(i+1 < len(argv))
            #utils.echo_msg_bold('{} {}'.format(len(argv), i+1))
            if i+1 < len(argv) and not argv[i + 1].startswith('-'):# and argv[i+1].split(':')[0] not in dataexts and argv[i+1].split('.')[-1] not in dataexts:
                archive_dirname = utils.str_or(argv[i + 1])
                #utils.echo_msg_bold(archive_dirname)
                i = i + 1
        elif arg[:2] == '-V':
            want_archive = True
            archive_dirname = utils.str_or(argv[2:])
            
        elif arg == '--cache-dir' or arg == '-D' or arg == '-cache-dir':
            cache_dir = utils.str_or(argv[i + 1], utils.cudem_cache)
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
        elif arg == '--quiet' or arg == '-q':
            want_verbose = False
        elif arg == '--help' or arg == '-h':
            print(datalists_usage())
            sys.exit(1)
        elif arg == '--version' or arg == '-v':
            print('{}, version {}'.format(
                os.path.basename(sys.argv[0]), cudem.__version__)
                  )
            sys.exit(1)
        elif arg[0] == '-':
            print(datalists_usage())
            sys.exit(0)
        else: dls.append(arg)#'"{}"'.format(arg)) # FIX THIS!!!
        
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

    stack_fltrs = [':'.join(f.split('/')) for f in stack_fltrs]
    pnt_fltrs = [':'.join(f.split('//')) for f in pnt_fltrs]
    if not i_regions: i_regions = [None]
    these_regions = regions.parse_cli_region(i_regions, want_verbose)
    for rn, this_region in enumerate(these_regions):
        ## buffer the region by `extend` if xy_inc is set
        ## this effects the output naming of masks/stacks!
        ## do we want this like in waffles where the output name
        ## does not include the -X extend buffer?
        if xy_inc[0] is not None and xy_inc[1] is not None and this_region is not None:
            this_region.buffer(
                x_bv=(utils.str2inc(xy_inc[0])*extend),
                y_bv=(utils.str2inc(xy_inc[1])*extend)
            )

        if len(dls) == 0:
            sys.stderr.write(datalists_usage())
            utils.echo_error_msg('you must specify some type of data')
        else:
            ## intiialze the input data. Treat data from CLI as a datalist.
            this_datalist = init_data(
                dls,
                region=this_region,
                src_srs=src_srs,
                dst_srs=dst_srs,
                xy_inc=xy_inc,
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
                    print(this_datalist.inf()) # output the datalist inf blob
                elif want_list:
                    this_datalist.echo() # output each dataset from the datalist
                elif want_region: # get the region and warp it if necessary
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
                    this_archive = this_datalist.archive_xyz(dirname=archive_dirname) # archive the datalist as xyz
                    if this_archive.numpts == 0:
                        utils.remove_glob('{}*'.format(this_archive.name))
                else:
                    try:
                        if want_separate: # process and dump each dataset independently
                            for this_entry in this_datalist.parse():
                                this_entry.dump_xyz()
                        else: # process and dump the datalist as a whole
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
