### visus.py - OpenVisus integration for CUDEM
##
## Copyright (c) 2023 - 2026 Regents of the University of Colorado
##
## This module requires the OpenVisus library:
## pip install OpenVisus
##

import os
import sys
import numpy as np
from osgeo import gdal
from cudem import utils
from cudem import gdalfun
from cudem.perspecto.perspecto import Perspecto

class OpenVisusViz(Perspecto):
    """Convert DEM to OpenVisus (.idx) format for interactive streaming/visualization.
    
    Convert DEMs into the OpenVisus Hierarchical Z-Order indexing format, 
    enabling real-time interactive viewing of datasets.
    
    Args:
        view (bool): If True, launch the OpenVisus viewer after conversion.
        compression (str): Compression type (zip, lz4, etc.). Default: zip.
    """
    
    def __init__(self, view=False, compression="zip", **kwargs):
        super().__init__(**kwargs)
        self.view = view
        self.compression = compression
        
    def run(self):
        try:
            import OpenVisus as ov
        except ImportError:
            utils.echo_error_msg("OpenVisus not installed. Run: pip install OpenVisus")
            return None

        ds = gdal.Open(self.src_dem)
        if ds is None: return None
        
        band = ds.GetRasterBand(1)
        width = ds.RasterXSize
        height = ds.RasterYSize
        
        ov_dtype = "float32" 
        
        idx_name = self.outfile if self.outfile.endswith(".idx") else f"{os.path.splitext(self.outfile)[0]}.idx"
        
        db = ov.CreateIdx(
            url=idx_name,
            dims=[width, height],
            fields=[ov.Field("elevation", ov_dtype)],
            compression=self.compression
        )

        if self.verbose:
            utils.echo_msg(f"Generating OpenVisus Index: {idx_name} [{width}x{height}]")

        chunk_size = 2048        
        for srcwin in utils.yield_srcwin((height, width), chunk_size, verbose=self.verbose):
            x, y, w, h = map(np.int64, srcwin)
                
            data = band.ReadAsArray(*srcwin).astype(np.float32)
           
            # OpenVisus expects [x_start, x_end) and [y_start, y_end)
            #p1 = ov.PointNi(x, y)
            #p2 = ov.PointNi(x + w, y + h)
            #box = ov.BoxNi(p1, p2)
            chunk_box = db.getLogicBox(x=[x, x+w], y=[y, y+h])
            
            # Write using the box argument
            #access = db.createAccess()
            #db.write(data)#, access, p1, p2)#, access=access, x=x, y=y)
            #db.write(data, access=access, box=box)
            db.write(data, logic_box=chunk_box)
                    
        if self.verbose:
            utils.echo_msg(f"Successfully created {idx_name}")
        
        # 4. Interactive Viewer
        if self.view:
            if self.verbose:
                utils.echo_msg("Launching OpenVisus Viewer...")
            try:
                dataset = ov.LoadDataset(idx_name)
                data = dataset.read(max_resolution=20) 
                if self.verbose:
                    utils.echo_msg(f"Dataset Loaded. Max Resolution: {dataset.getMaxResolution()}")
                
                # To launch a GUI here:
                from OpenVisus.gui import PyVisus
                PyVisus.Run()
                
            except Exception as e:
                utils.echo_error_msg(f"Viewer error: {e}")

        return idx_name
    

### End
