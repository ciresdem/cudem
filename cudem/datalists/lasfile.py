### lasfile.py 
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
##
## lasfile.py is part of CUDEM
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
## LAS/LAZ Lidar Data Parser
##
### Code:

import os
import numpy as np
import laspy as lp

from cudem import utils
from cudem import regions
from cudem.datalists.dlim import ElevationDataset

class LASFile(ElevationDataset):
    """
    Representing an LAS/LAZ dataset.
    Process LAS/LAZ lidar files using laspy.
    """

    def __init__(self, classes='2/29/40', **kwargs):
        super().__init__(**kwargs)
        
        ## List of lidar classes to retain (Default: Ground(2), LowNoise(29), BathymetricPoint(40))
        try:
            if isinstance(classes, str):
                self.classes = [int(x) for x in classes.split('/')]
            elif isinstance(classes, (list, tuple)):
                self.classes = [int(x) for x in classes]
            else:
                self.classes = []
        except Exception:
            self.classes = []

        ## Attempt to get SRS from header if not provided
        if self.src_srs is None:
            self.src_srs = self.get_epsg()
            if self.src_srs is None:
                self.src_srs = self.infos.src_srs

                
    def valid_p(self, fmts=None):
        """Check if self appears to be a valid dataset entry."""
        
        if self.fn is None: 
            return False
        
        if not os.path.exists(self.fn) or os.stat(self.fn).st_size == 0:
            return False

        try:
            ## Quick check if laspy can read the header
            with lp.open(self.fn) as lasf:
                pass
        except Exception as e:
            utils.echo_warning_msg(f'{self.fn} could not be opened by lasreader: {e}')
            return False
                        
        return True

    
    def get_epsg(self):
        """Attempt to parse EPSG/WKT from LAS VLRs."""
        
        try:
            with lp.open(self.fn) as lasf:
                for vlr in lasf.header.vlrs:
                    ## Record ID 2112 is "OGC Coordinate System WKT"
                    if vlr.record_id == 2112:
                        try:
                            ## Decode bytes if necessary
                            srs = vlr.string
                            if isinstance(srs, bytes):
                                return srs.decode('utf-8').strip('\0')
                            return srs
                        except:
                            pass
                    ## Record ID 34735 is "GeoKeyDirectoryTag" (GeoTIFF keys) - harder to parse directly here
                    ## without external libs, but laspy handles some of this internally in newer versions.
        except Exception:
            pass
            
        return None

    
    def generate_inf(self, make_grid=True, make_block_mean=False, block_inc=None):
        """Generate metadata (INF) for the LAS data.
        
        Uses the LAS header for fast bounds extraction.
        Falls back to full scan via parent class if grids/block-means are requested.
        """
        
        ## Quick Metadata from Header (Fast)
        try:
            with lp.open(self.fn) as lasf:
                self.infos.numpts = lasf.header.point_count
                
                ## Check for empty header bounds
                if lasf.header.x_min == 0 and lasf.header.x_max == 0:
                    ## Header might be empty, force scan
                    raise ValueError("Empty Header Bounds")

                this_region = regions.Region(
                    xmin=lasf.header.x_min, xmax=lasf.header.x_max,
                    ymin=lasf.header.y_min, ymax=lasf.header.y_max,
                    zmin=lasf.header.z_min, zmax=lasf.header.z_max
                )
                self.infos.minmax = this_region.export_as_list(include_z=True)
                self.infos.wkt = this_region.export_as_wkt()
                
                if self.infos.src_srs is None:
                    self.infos.src_srs = self.src_srs if self.src_srs else self.get_epsg()

        except Exception:
            ## Fallback to full scan if header read fails or bounds are bad
            pass

        ## If Grids are requested, we must scan the data.
        if make_grid or make_block_mean:
            return super().generate_inf(
                make_grid=make_grid, 
                make_block_mean=make_block_mean, 
                block_inc=block_inc
            )
        
        return self.infos

    
    def yield_points(self):
        """Yield points from the LAS file using chunked reading.
        Applies class filtering.
        """
        
        try:
            with lp.open(self.fn) as lasf:
                ## Iterate in chunks to handle large files efficiently
                for chunk in lasf.chunk_iterator(2_000_000):
                    ## Filter by classification if classes are set
                    if self.classes:
                        mask = np.isin(chunk.classification, self.classes)
                        points_x = chunk.x[mask]
                        points_y = chunk.y[mask]
                        points_z = chunk.z[mask]
                    else:
                        points_x = chunk.x
                        points_y = chunk.y
                        points_z = chunk.z
                    
                    if len(points_x) == 0:
                        continue

                    ## Create structured array
                    dataset = np.column_stack((points_x, points_y, points_z))
                    points = np.rec.fromrecords(dataset, names='x, y, z')
                    
                    yield points
                    
        except Exception as e:
            utils.echo_warning_msg(f'Could not read points from lasfile {self.fn}: {e}')

### End
