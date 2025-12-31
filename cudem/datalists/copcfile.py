import json
from cudem import utils
from cudem import regions
from cudem.datalists.dlim import ElevationDataset

try:
    import pdal
    HAS_PDAL = True
except ImportError:
    HAS_PDAL = False

class COPCFile(ElevationDataset):
    """
    Cloud Optimized Point Cloud (COPC) Parser.
    Streams points from a remote or local COPC file using PDAL.
    
    Parameters:
    -----------
    remote : bool
        If True, treats fn as a URL (s3://, http://, etc.)
    """
    
    def __init__(self, remote=False, **kwargs):
        super().__init__(**kwargs)
        self.remote = remote
        self.data_format = 300 # Assign a unique ID for Factory

    def generate_inf(self):
        """
        Generate metadata using PDAL's quickinfo (hexbin/stats).
        COPC headers contain bounds, so this is fast.
        """
        try:
            # Simple PDAL pipeline to get metadata
            pipeline_json = [
                {
                    "type": "readers.copc",
                    "filename": self.fn,
                    "count": 0  # Metadata only
                }
            ]
            pipeline = pdal.Pipeline(json.dumps(pipeline_json))
            pipeline.execute()
            meta = pipeline.metadata['metadata']['readers.copc']
            
            # Extract Bounds
            # PDAL metadata structure varies slightly, checking standard keys
            bounds = meta.get('bounds')
            if bounds:
                # bounds are usually {'minx': ..., 'maxx': ..., ...}
                self.infos.minmax = [
                    bounds['minx'], bounds['maxx'],
                    bounds['miny'], bounds['maxy'],
                    bounds['minz'], bounds['maxz']
                ]
                self.infos.numpts = meta.get('count', 0)
                
            # Extract SRS
            srs_wkt = meta.get('srs', {}).get('horizontal')
            if srs_wkt:
                self.infos.src_srs = srs_wkt
                if self.src_srs is None:
                    self.src_srs = srs_wkt

            r = regions.Region().from_list(self.infos.minmax)
            self.infos.wkt = r.export_as_wkt()

        except Exception as e:
            utils.echo_warning_msg(f"Failed to generate INF for COPC {self.fn}: {e}")
            
        return self.infos

    def yield_points(self):
        """
        Stream points using PDAL.
        PDAL handles the chunking and requests for us.
        """
        # Build PDAL Pipeline
        # We can add filters here (e.g. range filter for self.region)
        reader_args = {
            "type": "readers.copc",
            "filename": self.fn
        }
        
        # If we have a region, PDAL can request only those tiles!
        # This is the "Cloud Optimized" power.
        if self.region is not None:
            # PDAL expects bounds in format: "([minx, maxx], [miny, maxy])"
            # or just a polygon WKT
            reader_args["bounds"] = (
                f"([{self.region.xmin}, {self.region.xmax}], "
                f"[{self.region.ymin}, {self.region.ymax}])"
            )

        pipeline_json = [reader_args]
        
        try:
            # Stream via PDAL iterator if supported or standard execute
            # For massive files, we might want to implement a custom iterator
            # that fetches standard chunks if PDAL python bindings load all into RAM.
            # Standard pdal.Pipeline loads results into numpy arrays.
            
            pipeline = pdal.Pipeline(json.dumps(pipeline_json))
            
            # Execute (Note: For massive regions, this might still be too much RAM.
            # Ideally, we'd use PDAL's stream mode, but python bindings support is limited.
            # A workaround is to break self.region into smaller tiles here.)
            count = pipeline.execute()
            
            arrays = pipeline.arrays
            for arr in arrays:
                # Convert PDAL structured array to CUDEM convention
                # PDAL fields: X, Y, Z, Intensity, etc.
                # CUDEM expects: x, y, z
                
                # Create view/copy with correct names
                # PDAL uses Capital X, Y, Z usually
                import numpy as np
                out_ds = np.rec.fromarrays(
                    [arr['X'], arr['Y'], arr['Z']],
                    names=['x', 'y', 'z']
                )
                yield out_ds

        except Exception as e:
            utils.echo_error_msg(f"PDAL Pipeline failed: {e}")

### End
