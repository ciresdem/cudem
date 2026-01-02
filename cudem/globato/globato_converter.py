### globato_converter.py
##
## Copyright (c) 2025 - 2026 Regents of the University of Colorado
##
## globato_converter.py is part of CUDEM
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
## Utilities to convert between GlobatoStack (HDF5) and GdalStack (GeoTIFF) formats.
## Preserves data, metadata, and masks where possible.
##
### Code:

import os
import sys
import argparse
import numpy as np
import h5py as h5
from osgeo import gdal, osr

from cudem import utils
from cudem import srsfun
from cudem import gdalfun
from cudem.globato.globato import GlobatoStack, GdalStack

def globato_to_gdal(h5_path, tif_path=None, co=None, verbose=False):
    """
    Convert a Globato HDF5 stack to a GDAL GeoTIFF stack.
    
    This flattens the '/stack' group of the HDF5 file into a multi-band GeoTIFF.
    It also attempts to convert the '/mask' group into a sidecar mask file (_msk.tif).
    
    Args:
        h5_path (str): Path to the source .h5 (gbt) file.
        tif_path (str): Path to the destination .tif file. (Defaults to h5_path with .tif extension).
        co (list): GDAL Creation Options (e.g. ["COMPRESS=LZW"]).
        verbose (bool): Verbose output.
        
    Returns:
        str: Path to the generated GeoTIFF.
    """

    #from cudem.globato import GlobatoStack
    
    if not os.path.exists(h5_path):
        utils.echo_error_msg(f"Source file {h5_path} does not exist.")
        return None

    if tif_path is None:
        tif_path = os.path.splitext(h5_path)[0] + ".tif"
        
    if co is None:
        co = ["COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES", "PREDICTOR=2"]

    with GlobatoStack(h5_path) as gs:
        if verbose: utils.echo_msg(f"Converting {h5_path} -> {tif_path}...")
        
        ## Metadata extraction
        gt = gs.geotransform
        wkt = gs.crs_wkt
        
        ## Get dimensions from the z array
        z_ds = gs.stack['z']
        y_size, x_size = z_ds.shape
        
        ## Initialize GeoTIFF
        ## GdalStack expects bands: 1:z, 2:count, 3:weights, 4:uncertainty, 5:src_uncertainty, 6:x, 7:y
        band_map = ['z', 'count', 'weights', 'uncertainty', 'src_uncertainty', 'x', 'y']
        
        driver = gdal.GetDriverByName("GTiff")
        ds = driver.Create(tif_path, x_size, y_size, len(band_map), gdal.GDT_Float32, options=co)
        
        if ds is None:
            utils.echo_error_msg(f"Failed to create output file {tif_path}")
            return None
            
        ds.SetGeoTransform(gt)
        if wkt:
            ds.SetProjection(wkt)
            
        ## Transfer Data Bands
        ## Process block-by-block to save memory
        block_size = 256
        
        for i, key in enumerate(band_map):
            band_idx = i + 1
            dst_band = ds.GetRasterBand(band_idx)
            dst_band.SetDescription(key)
            dst_band.SetNoDataValue(np.nan)
            
            if key not in gs.stack:
                if verbose: utils.echo_warning_msg(f"Dataset '{key}' missing in HDF5, skipping band.")
                continue
                
            src_dset = gs.stack[key]
            
            ## Simple chunked copy
            for y in range(0, y_size, block_size):
                if y + block_size < y_size:
                    rows = block_size
                else:
                    rows = y_size - y
                    
                data_chunk = src_dset[y:y+rows, :]
                dst_band.WriteArray(data_chunk, 0, y)
                
        ds.FlushCache()
        ds = None # Close
        
        ## Transfer Masks (Sidecar)
        ## Globato stores masks in /mask group. GdalStack expects {tif_path}_msk.tif
        mask_keys = list(gs.masks.keys())
        mask_keys = None # TODO: need to fix the processing, as the h5 groups can go deep.
        if mask_keys:
            mask_path = os.path.splitext(tif_path)[0] + "_msk.tif"
            if verbose: utils.echo_msg(f"Generating mask sidecar: {mask_path}")
            
            ## Prioritize 'full_dataset_mask' as Band 1
            ordered_keys = []
            if 'full_dataset_mask' in mask_keys:
                ordered_keys.append('full_dataset_mask')
                mask_keys.remove('full_dataset_mask')
            ordered_keys.extend(mask_keys)
            
            msk_ds = driver.Create(mask_path, x_size, y_size, len(ordered_keys), gdal.GDT_Byte, options=co)
            msk_ds.SetGeoTransform(gt)
            if wkt: msk_ds.SetProjection(wkt)
            
            for i, key in enumerate(ordered_keys):
                band_idx = i + 1
                dst_band = msk_ds.GetRasterBand(band_idx)
                dst_band.SetDescription(key)
                dst_band.SetNoDataValue(0)
                
                src_dset = gs.masks[key]
                
                ## Copy attributes as metadata
                md = {}
                for attr_key, attr_val in src_dset.attrs.items():
                    ## Skip internal netcdf/hdf attributes
                    if attr_key not in ['DIMENSION_LIST', 'grid_mapping']:
                        md[attr_key] = str(attr_val)
                dst_band.SetMetadata(md)
                
                for y in range(0, y_size, block_size):
                    if y + block_size < y_size: rows = block_size
                    else: rows = y_size - y
                    dst_band.WriteArray(src_dset[y:y+rows, :], 0, y)

                    
            msk_ds = None

    return tif_path


def gdal_to_globato(tif_path, h5_path=None, verbose=False):
    """Convert a GDAL GeoTIFF stack to a Globato HDF5 stack.
    
    This reconstructs the hierarchical structure. 
    Note: The '/sums' group will be populated with a copy of '/stack' since 
    raw sums are usually not preserved in the GeoTIFF output format.
    
    Args:
        tif_path (str): Path to source .tif file.
        h5_path (str): Path to dest .h5 file.
        verbose (bool): Verbose output.
        
    Returns:
        str: Path to the generated HDF5 file.
    """

    #from cudem.globato import GdalStack
    
    if not os.path.exists(tif_path):
        utils.echo_error_msg(f"Source file {tif_path} does not exist.")
        return None
        
    if h5_path is None:
        h5_path = os.path.splitext(tif_path)[0] + ".h5"
        
    ## Use existing GdalStack interface to read safely
    gs_in = GdalStack(tif_path)
    gs_in.open()
    
    gt = gs_in.geotransform
    wkt = gs_in.projection
    y_size, x_size = gs_in.shape
    
    ## Initialize HDF5
    with h5.File(h5_path, 'w') as h5_out:
        if verbose: utils.echo_msg(f"Converting {tif_path} -> {h5_path}...")
        
        h5_out.attrs['short_name'] = 'globato'
        
        ## CRS
        crs_dset = h5_out.create_dataset('crs', dtype=h5.string_dtype())
        crs_dset.attrs['GeoTransform'] = ' '.join([str(x) for x in gt])
        if wkt:
            crs_dset.attrs['crs_wkt'] = wkt
            
        ## Coordinates
        ## GeoTransform: [top_left_x, w_e_res, 0, top_left_y, 0, n_s_res]
        ## Pixel Center convention for coordinates
        lon_start = gt[0] + (gt[1] / 2)
        lon_inc = gt[1]
        lon_end = gt[0] + (gt[1] * x_size)
        
        lat_start = gt[3] + (gt[5] / 2)
        lat_inc = gt[5]
        lat_end = gt[3] + (gt[5] * y_size)
        
        lat_array = np.arange(lat_start, lat_end, lat_inc)
        ## Ensure exact size match (floating point errors)
        if len(lat_array) > y_size: lat_array = lat_array[:y_size]
        
        lon_array = np.arange(lon_start, lon_end, lon_inc)
        if len(lon_array) > x_size: lon_array = lon_array[:x_size]
        
        lat_dset = h5_out.create_dataset('lat', data=lat_array)
        lat_dset.make_scale('latitude')
        
        lon_dset = h5_out.create_dataset('lon', data=lon_array)
        lon_dset.make_scale('longitude')
        
        ##  Create Groups
        stack_grp = h5_out.create_group('stack')
        sums_grp = h5_out.create_group('sums')
        mask_grp = h5_out.create_group('mask')
        h5_out.create_group('datasets')
        
        ## Transfer Data
        ## Keys match GdalStack.BAND_MAP
        band_keys = ['z', 'count', 'weights', 'uncertainty', 'src_uncertainty', 'x', 'y']
        
        block_size = 256
        chunks = (min(100, y_size), x_size)
        
        for key in band_keys:
            ## Create HDF5 dataset
            dset = stack_grp.create_dataset(
                key, shape=(y_size, x_size), dtype=np.float32,
                chunks=chunks, compression='lzf'
            )
            dset.dims[0].attach_scale(lat_dset)
            dset.dims[1].attach_scale(lon_dset)
            
            ## Also create in sums (placeholder)
            sums_dset = sums_grp.create_dataset(
                key, shape=(y_size, x_size), dtype=np.float32,
                chunks=chunks, compression='lzf'
            )
            
            ## Copy Data from TIFF
            ## GdalStack accesses bands by name via property (e.g. gs_in.z)
            idx = GdalStack.BAND_MAP.get(key)
            if idx and idx <= gs_in._ds.RasterCount:
                band = gs_in._ds.GetRasterBand(idx)
                ndv = band.GetNoDataValue()
                
                for y in range(0, y_size, block_size):
                    if y + block_size < y_size: rows = block_size
                    else: rows = y_size - y
                    
                    data = band.ReadAsArray(0, y, x_size, rows)
                    if ndv is not None:
                        data[data == ndv] = np.nan
                        
                    dset[y:y+rows, :] = data
                    sums_dset[y:y+rows, :] = data # Mirror to sums
                    
        ## Transfer Masks (if sidecar exists)
        if gs_in.has_mask:
            mask_ds = gs_in._mask_ds
            for b in range(1, mask_ds.RasterCount + 1):
                m_band = mask_ds.GetRasterBand(b)
                name = m_band.GetDescription()
                if not name: name = f"mask_{b}"
                
                m_dset = mask_grp.create_dataset(
                    name, shape=(y_size, x_size), dtype=np.uint8,
                    chunks=chunks, compression='lzf'
                )
                m_dset.dims[0].attach_scale(lat_dset)
                m_dset.dims[1].attach_scale(lon_dset)
                
                ## Copy Metadata
                md = m_band.GetMetadata()
                for k, v in md.items():
                    m_dset.attrs[k] = v
                    
                ## Copy Data
                for y in range(0, y_size, block_size):
                    if y + block_size < y_size: rows = block_size
                    else: rows = y_size - y
                    m_dset[y:y+rows, :] = m_band.ReadAsArray(0, y, x_size, rows)

    gs_in.close()
    return h5_path


## ==============================================
## Command-line Interface (CLI)
## $ globato_converter
##
## globato_converter cli
## ==============================================
def globato_converter_cli():
    parser = argparse.ArgumentParser(
        description="Convert between Globato (HDF5) and GDAL (GeoTIFF) stack formats.",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="CUDEM home page: <http://cudem.colorado.edu>"
    )
    
    parser.add_argument(
        "src", 
        help="Input file path (.h5 or .tif)"
    )
    parser.add_argument(
        "dst", 
        nargs="?", 
        help="Output file path (optional)"
    )
    parser.add_argument(
        "-co", "--creation-option", 
        action="append", 
        help="GDAL creation option (e.g. -co COMPRESS=LZW). Only for HDF5 -> TIFF."
    )
    parser.add_argument(
        "-q", "--quiet", 
        action="store_true", 
        help="Suppress output."
    )
    parser.add_argument(
        "--to-gdal", 
        action="store_true", 
        help="Force conversion to GDAL GeoTIFF."
    )
    parser.add_argument(
        "--to-globato", 
        action="store_true", 
        help="Force conversion to Globato HDF5."
    )
    
    args = parser.parse_args()
    verbose = not args.quiet
    
    ## Detect conversion direction
    src_ext = os.path.splitext(args.src)[1].lower()
    
    if args.to_gdal or src_ext in ['.h5', '.hdf5', '.he5', 'gbt']:
        globato_to_gdal(
            args.src, 
            args.dst, 
            co=args.creation_option, 
            verbose=verbose
        )
        
    elif args.to_globato or src_ext in ['.tif', '.tiff', '.gtiff']:
        gdal_to_globato(
            args.src, 
            args.dst, 
            verbose=verbose
        )
        
    else:
        utils.echo_error_msg(
            "Could not determine conversion direction from file extension. "
            "Please use --to-gdal or --to-globato."
        )

        
if __name__ == "__main__":
    globato_converter_cli()

    
### End
