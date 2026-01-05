### spatial_metadata.py 
##
## Copyright (c) 2018 - 2026 Regents of the University of Colorado
##
## spatial_metadata.py is part of CUDEM
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
## Tools for generating spatial metadata (footprints) from raster masks.
## Supports both standard GDAL Datasets and Globato HDF5 Stacks.
##
### Code:

import os
import sys
import argparse
import numpy as np
from osgeo import gdal, ogr, osr
from typing import Optional, Dict, List, Union

from cudem import utils
from cudem import gdalfun
from cudem import srsfun

try:
    from cudem.globato import GlobatoStack
    HAS_GLOBATO = True
except ImportError:
    HAS_GLOBATO = False

## ==============================================
##  The Spatial Metadata Generator (from globato/stacks)
## ==============================================
class SpatialMetadataGenerator:
    """Generates vector spatial metadata (footprints) from raster masks.
    Can accept a GDAL Dataset or a GlobatoStack.
    """

    def __init__(self, input_source: Union[str, gdal.Dataset, 'GlobatoStack'], 
                 output_vector: Optional[str] = None, 
                 ogr_format: str = 'GPKG', 
                 verbose: bool = True):
        
        self.input_source = input_source
        self.output_vector = output_vector
        self.ogr_format = ogr_format
        self.verbose = verbose
        
        ## Determine source type
        self.is_globato = False
        self.src_ds = None
        self.globato_stack = None
        
        if HAS_GLOBATO and isinstance(input_source, GlobatoStack):
            self.is_globato = True
            self.globato_stack = input_source
            if not self.globato_stack.is_open:
                self.globato_stack.open()
                
        elif isinstance(input_source, str):
            # Check file extension or try opening
            if input_source.endswith('.h5'):
                if not HAS_GLOBATO:
                    raise ImportError("Globato module needed for .h5 files.")
                self.is_globato = True
                self.globato_stack = GlobatoStack(input_source)
                self.globato_stack.open()
            else:
                self.src_ds = gdal.Open(input_source)
                if not self.src_ds:
                    raise IOError(f"Could not open {input_source}")
                    
        elif isinstance(input_source, gdal.Dataset):
            self.src_ds = input_source
        else:
            raise ValueError("Invalid input source type.")

        ## Set default output name
        if self.output_vector is None:
            base = ""
            if self.is_globato:
                base = os.path.splitext(os.path.basename(self.globato_stack.filepath))[0]
            else:
                base = os.path.splitext(os.path.basename(self.src_ds.GetDescription()))[0]
            
            ext = gdalfun.ogr_fext(self.ogr_format)
            self.output_vector = f"{base}_sm.{ext}"

            
    def _get_globato_geo_info(self):
        """Extract geo-info from GlobatoStack."""
        
        gt = self.globato_stack.geotransform
        wkt = self.globato_stack.crs_wkt
        return gt, wkt

    
    def _get_h5_mask_as_mem_ds(self, mask_name: str) -> gdal.Dataset:
        """Convert an HDF5 mask array into an in-memory GDAL Dataset 
        so gdal.Polygonize can process it.
        """
        
        if not self.is_globato:
            raise ValueError("Not a Globato source.")
            
        mask_grp = self.globato_stack.masks
        if mask_name not in mask_grp:
            raise KeyError(f"Mask {mask_name} not found in stack.")
            
        data = mask_grp[mask_name][:] # Read into numpy array
        rows, cols = data.shape
        
        gt, wkt = self._get_globato_geo_info()
        
        drv = gdal.GetDriverByName('MEM')
        ds = drv.Create('', cols, rows, 1, gdal.GDT_Byte)
        ds.SetGeoTransform(gt)
        if wkt:
            ds.SetProjection(wkt)
            
        band = ds.GetRasterBand(1)
        band.WriteArray(data)
        band.SetNoDataValue(0)
        
        ## Copy attributes as metadata
        attrs = dict(mask_grp[mask_name].attrs)
        ## Convert non-string attrs
        str_attrs = {k: str(v) for k, v in attrs.items()}
        band.SetMetadata(str_attrs)
        
        return ds

    
    def scan_masks(self, skip_band: str = 'Full Data Mask') -> Dict:
        """Scan the source for mask layers/bands and return their metadata.
        """
        
        infos = {}
        
        if self.is_globato:
            ## Globato Mode: Iterate HDF5 keys
            mask_names = self.globato_stack.get_mask_names()
            with utils.ccp(desc='Scanning Globato masks', total=len(mask_names), leave=self.verbose) as pbar:
                for name in mask_names:
                    # Retrieve attributes directly from HDF5 dataset
                    attrs = dict(self.globato_stack.masks[name].attrs)
                    clean_md = {k: str(v) for k, v in attrs.items()}
                    
                    ## Normalize metadata structure
                    infos[name] = {
                        'source': 'hdf5',
                        'metadata': clean_md
                    }
                    pbar.update()
        else:
            ## GDAL Mode: Iterate Raster Bands
            count = self.src_ds.RasterCount
            with utils.ccp(desc='Scanning Raster bands', total=count, leave=self.verbose) as pbar:
                for b in range(1, count + 1):
                    band = self.src_ds.GetRasterBand(b)
                    name = band.GetDescription()
                    if name == skip_band: 
                        pbar.update()
                        continue
                        
                    md = band.GetMetadata()
                    infos[name] = {
                        'band_index': b,
                        'metadata': md
                    }
                    pbar.update()
                    
        return infos

    
    def polygonize(self) -> str:
        """Main execution method. Generates the vector file.
        """
        
        ## Cleanup existing
        utils.remove_glob(f"{os.path.splitext(self.output_vector)[0]}.*")

        ## Create Output Vector
        drv = ogr.GetDriverByName(self.ogr_format)
        ds = drv.CreateDataSource(self.output_vector)
        if ds is None:
            raise IOError(f"Could not create {self.output_vector}")

        ## Determine SRS
        srs = osr.SpatialReference()
        if self.is_globato:
            wkt = self.globato_stack.crs_wkt
            if wkt: srs.ImportFromWkt(wkt)
        else:
            wkt = self.src_ds.GetProjectionRef()
            if wkt: srs.ImportFromWkt(wkt)

        layer = ds.CreateLayer('footprints', srs, ogr.wkbMultiPolygon)
        
        ## Standard Fields
        layer.CreateField(ogr.FieldDefn('Name', ogr.OFTString))
        layer.CreateField(ogr.FieldDefn('Title', ogr.OFTString))
        
        ## Scan to find all potential metadata fields
        mask_infos = self.scan_masks()
        all_keys = set()
        for info in mask_infos.values():
            all_keys.update(info['metadata'].keys())
            
        ## Add dynamic fields
        current_fields = ['Name', 'Title']
        for k in all_keys:
            if k not in current_fields:
                layer.CreateField(ogr.FieldDefn(k, ogr.OFTString))

        ## --- Polygonize Loop ---
        with utils.ccp(desc='Polygonizing masks', total=len(mask_infos), leave=self.verbose) as pbar:
            for name, info in mask_infos.items():
                
                ## Get Source Band/Dataset
                tmp_src_band = None
                tmp_ds_ref = None # Keep reference to prevent GC during operation
                
                if self.is_globato:
                    ## Create MEM ds from HDF5
                    tmp_ds_ref = self._get_h5_mask_as_mem_ds(name)
                    tmp_src_band = tmp_ds_ref.GetRasterBand(1)
                else:
                    tmp_src_band = self.src_ds.GetRasterBand(info['band_index'])

                ## Create Temp Layer for this mask
                mem_drv = ogr.GetDriverByName('Memory')
                tmp_ogr_ds = mem_drv.CreateDataSource('mem')
                tmp_layer = tmp_ogr_ds.CreateLayer('poly', srs, ogr.wkbMultiPolygon)
                tmp_layer.CreateField(ogr.FieldDefn('DN', ogr.OFTInteger))
                
                ## Run Polygonize
                gdal.Polygonize(tmp_src_band, None, tmp_layer, 0, [], callback=None)
                
                ## Merge Polygons (Union) to get single footprint for this source
                ## Note: cudem.gdalfun.ogr_mask_union logic assumed here or simplified
                if tmp_layer.GetFeatureCount() > 0:
                    geom_union = ogr.Geometry(ogr.wkbMultiPolygon)
                    for feat in tmp_layer:
                        geom = feat.GetGeometryRef()
                        if geom:
                            geom_union.AddGeometry(geom)
                    
                    ## Clean Union
                    geom_union = geom_union.UnionCascaded()
                    
                    ## Create Final Feature
                    out_feat = ogr.Feature(layer.GetLayerDefn())
                    out_feat.SetGeometry(geom_union)
                    out_feat.SetField('Name', name)
                    
                    ## Set Metadata Attributes
                    md = info['metadata']
                    for k, v in md.items():
                        ## Truncate if needed based on driver limits? (GPKG is lenient)
                        out_feat.SetField(k, str(v))
                        
                    layer.CreateFeature(out_feat)
                
                tmp_ds_ref = None
                pbar.update()

        ds = None
        if self.verbose:
            utils.echo_msg(f"Created spatial metadata: {self.output_vector}")
            
        return self.output_vector

    
## ==============================================
## Command-line Interface (CLI)
## $ spatial_metadata
##
## spatial_metadata cli
## ==============================================
def spatial_metadata_cli():
    parser = argparse.ArgumentParser(
        description="Generate Spatial Metadata (Footprints) from Raster/Globato Masks",
        epilog="CUDEM home page: <http://cudem.colorado.edu>",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("input", help="Input file (GeoTIFF or .h5 Globato Stack)")
    parser.add_argument("-o", "--output", help="Output vector file (default: input_sm.gpkg)")
    parser.add_argument("-f", "--format", default="GPKG", help="Output OGR Format (default: GPKG)")
    parser.add_argument("--quiet", action="store_true", help="Suppress progress output")

    args = parser.parse_args()
    
    try:
        gen = SpatialMetadataGenerator(
            input_source=args.input,
            output_vector=args.output,
            ogr_format=args.format,
            verbose=not args.quiet
        )
        gen.polygonize()
        
    except Exception as e:
        sys.stderr.write(f"Error: {e}\n")
        sys.exit(1)

if __name__ == "__main__":
    spatial_metadata_cli()

### End
