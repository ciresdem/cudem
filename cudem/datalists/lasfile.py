### lasfile.py 
##
## Copyright (c) 2015 - 2026 Regents of the University of Colorado
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
## LAS/LAZ Lidar Data Parser.
## Handles both:
## * Local LAS/LAZ files via laspy
## * Remote COPC (Cloud Optimized Point Cloud) files via HTTP Range Requests
## * Remote standard LAS/LAZ files (Fallback -> Download -> Process)
##
### Code:

import os
import io
import struct
import numpy as np
import laspy as lp

from cudem import utils
from cudem import regions
from cudem.fetches import fetches
from cudem.datalists.dlim import ElevationDataset

class LASFile(ElevationDataset):
    """
    Representing an LAS/LAZ dataset.
    Process LAS/LAZ lidar files using laspy.
    
    Supports:
    - Local Files: Standard processing via laspy.
    - Remote URLs (COPC): Intelligent chunk-based streaming via HTTP Range Requests.
    - Remote URLs (Standard): Automatic download and local processing.
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

        ## Check for Remote/COPC status
        self.is_remote = False
        if self.fn is not None and self.fn.startswith(('http', 'https', 'ftp')):
            self.is_remote = True
        
        ## COPC specific state
        self.copc_info = None
        self.chunk_table = []
        self.scale = None
        self.offset = None

        ## Attempt to get SRS from header if not provided
        if self.src_srs is None:
            if self.is_remote:
                pass 
            else:
                self.src_srs = self.get_epsg()
            
            if self.src_srs is None:
                self.src_srs = self.infos.src_srs

                
    def valid_p(self, fmts=None):
        """Check if self appears to be a valid dataset entry."""
        
        if self.fn is None: 
            return False
        
        if self.is_remote:
            return True
        else:
            if not os.path.exists(self.fn) or os.stat(self.fn).st_size == 0:
                return False
            try:
                with lp.open(self.fn) as lasf:
                    pass
            except Exception as e:
                utils.echo_warning_msg(f'{self.fn} could not be opened by lasreader: {e}')
                return False
        return True


    def get_epsg(self):
        """Attempt to parse EPSG/WKT from LAS Header using laspy."""
        
        if self.is_remote: return None
        
        try:
            with lp.open(self.fn) as lasf:
                ## Try Native Laspy CRS parsing laspy >= 2.0
                ## This handles both OGC WKT and GeoTIFF VLRs automatically.
                try:
                    crs = lasf.header.parse_crs()
                    if crs is not None:
                        return crs.to_wkt()
                except Exception:
                    pass

                ## Fallback: Manual VLR Check (Legacy)
                for vlr in lasf.header.vlrs:
                    ## Record ID 2112 is "OGC Coordinate System WKT"
                    if vlr.record_id == 2112:
                        try:
                            srs = vlr.string
                            if isinstance(srs, bytes):
                                return srs.decode('utf-8').strip('\0')
                            return srs
                        except:
                            pass
        except Exception:
            pass
            
        return None    

    
    ## ==============================================
    ## COPC / Remote Helper Methods
    ## ==============================================
    def _fetch_bytes(self, start, length):
        """Fetch a byte range using cudem.fetches"""
        
        fetcher = fetches.Fetch(url=self.fn, verbose=self.verbose)
        fetcher.headers['Range'] = f'bytes={start}-{start + length - 1}'
        req = fetcher.fetch_req()
        if req is not None and req.status_code in [200, 206]:
            return req.content
        return None

    
    def _fetch_and_switch_to_local(self):
        """Fallback: Download the remote file to cache and switch 
        this instance to treat it as a local file.
        """
        
        if not self.cache_dir:
            self.cache_dir = '.'
            
        ## Clean filename (handle query strings)
        clean_name = self.fn.split('?')[0].split('/')[-1]
        local_path = os.path.join(self.cache_dir, clean_name)
        
        if self.verbose:
            utils.echo_msg(f"Downloading non-COPC file to {local_path}...")
            
        fetcher = fetches.Fetch(url=self.fn, verbose=self.verbose)
        if fetcher.fetch_file(local_path) == 0:
            self.fn = local_path
            self.is_remote = False
            return True
        return False


    def _parse_copc_metadata(self):
        """Parse LAS Header and COPC VLRs from remote source."""
        
        if self.copc_info is not None: return True

        ## Fetch LAS Header (375 bytes)
        header_bytes = self._fetch_bytes(0, 375)
        if not header_bytes: return False

        ## Parse Scale/Offset
        self.scale = struct.unpack('<3d', header_bytes[131:155])
        self.offset = struct.unpack('<3d', header_bytes[155:179])

        vlr_offset = struct.unpack('<I', header_bytes[96:100])[0]
        num_vlrs = struct.unpack('<I', header_bytes[100:104])[0]
        
        ## Fetch VLR block (heuristic 2000 bytes)
        vlr_chunk = self._fetch_bytes(vlr_offset, 2000)
        curr_pos = 0
        
        for _ in range(num_vlrs):
            if curr_pos + 54 > len(vlr_chunk): break
            
            ## VLR Header: Reserved(2), UserID(16), RecordID(2), Length(2)
            user_id_bytes = vlr_chunk[curr_pos+2 : curr_pos+18]
            record_id = struct.unpack('<H', vlr_chunk[curr_pos+18 : curr_pos+20])[0]
            rec_len = struct.unpack('<H', vlr_chunk[curr_pos+20 : curr_pos+22])[0]
            
            data_pos = curr_pos + 54
            
            ## Check for COPC Info (User: "entwine", Record: 1)
            ## Compare bytes directly to avoid UnicodeDecodeError
            if b'entwine' in user_id_bytes and record_id == 1:
                copc_data = vlr_chunk[data_pos : data_pos + rec_len]
                vals = struct.unpack('<4d d Q Q', copc_data)
                self.copc_info = {
                    'center': vals[0:3],
                    'half_size': vals[3],
                    'spacing': vals[4],
                    'root_hier_offset': vals[5],
                    'root_hier_len': vals[6]
                }
            
            ## Check for WKT SRS (Record 2112)
            if record_id == 2112:
                try:
                    wkt_bytes = vlr_chunk[data_pos : data_pos + rec_len]
                    ## Use 'ignore' or 'replace' to handle bad characters safely
                    self.src_srs = wkt_bytes.decode('utf-8', errors='ignore').strip('\0')
                except: pass

            curr_pos += 54 + rec_len

        return self.copc_info is not None
    
    
    def _traverse_copc_hierarchy(self, offset, length):
        """Recursively fetch and parse hierarchy pages."""
        
        page_bytes = self._fetch_bytes(offset, length)
        if not page_bytes: return

        num_nodes = len(page_bytes) // 32
        for i in range(num_nodes):
            base = i * 32
            d, x, y, z = struct.unpack('<4i', page_bytes[base : base+16])
            part2 = page_bytes[base+16 : base+32]
            offset_val = int.from_bytes(part2[0:6], 'little')
            byte_size = int.from_bytes(part2[6:10], 'little')
            point_count = int.from_bytes(part2[10:14], 'little')

            ## Calculate Node Bounds
            node_size = (self.copc_info['half_size'] * 2) / (2**d)
            node_min_x = (self.copc_info['center'][0] - self.copc_info['half_size']) + (x * node_size)
            node_max_x = node_min_x + node_size
            node_min_y = (self.copc_info['center'][1] - self.copc_info['half_size']) + (y * node_size)
            node_max_y = node_min_y + node_size
            
            node_region = regions.Region().from_list([node_min_x, node_max_x, node_min_y, node_max_y])
            
            ## Intersection Check
            if self.region is not None:
                if not regions.regions_intersect_p(self.region, node_region):
                    continue

            ## Leaf Node with data
            if byte_size > 0 and point_count > 0:
                self.chunk_table.append({
                    'offset': offset_val,
                    'length': byte_size,
                    'count': point_count
                })

                
    ## ==============================================
    ## Main Methods (Generate INF & Yield)
    ## ==============================================
    def generate_inf(self, make_grid=True, make_block_mean=False, block_inc=None):
        """Generate metadata (INF). Dispatches to Local or Remote logic."""
        
        if self.is_remote:
            ## Try COPC Parsing
            if self._parse_copc_metadata():
                c = self.copc_info['center']
                hs = self.copc_info['half_size']
                self.infos.minmax = [
                    c[0] - hs, c[0] + hs,
                    c[1] - hs, c[1] + hs,
                    c[2] - hs, c[2] + hs
                ]
                self.infos.wkt = regions.Region().from_list(self.infos.minmax).export_as_wkt()
                
                if self.infos.src_srs is None and self.src_srs is not None:
                    self.infos.src_srs = self.src_srs
            else:
                ## Fallback: Download and recurse
                if self.verbose:
                    utils.echo_warning_msg(f"COPC parse failed for {self.fn}. Downloading...")
                
                if self._fetch_and_switch_to_local():
                    ## Now self.is_remote is False, recurse to use local logic
                    return self.generate_inf(make_grid, make_block_mean, block_inc)
                else:
                    utils.echo_error_msg(f"Failed to download remote file {self.fn}")

        else:
            ## Local
            try:
                with lp.open(self.fn) as lasf:
                    self.infos.numpts = lasf.header.point_count
                    if lasf.header.x_min == 0 and lasf.header.x_max == 0:
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
                pass

        if make_grid or make_block_mean:
            return super().generate_inf(
                make_grid=make_grid, 
                make_block_mean=make_block_mean, 
                block_inc=block_inc
            )
        
        return self.infos

    
    def yield_points(self):
        """Yield points. Dispatches to Local or Remote logic."""
        
        if self.is_remote:
            ## COPC check (if not already done)
            if self._parse_copc_metadata():
                yield from self._yield_points_copc()
            else:
                ## Fallback
                if self.verbose:
                    utils.echo_warning_msg(f"File not COPC compatible: {self.fn}. Downloading...")
                
                if self._fetch_and_switch_to_local():
                    yield from self._yield_points_local()
        else:
            yield from self._yield_points_local()

            
    def _yield_points_local(self):
        """Yield points from local file using standard laspy."""
        
        try:
            with lp.open(self.fn) as lasf:
                for chunk in lasf.chunk_iterator(2_000_000):
                    if self.classes:
                        mask = np.isin(chunk.classification, self.classes)
                        points_x = chunk.x[mask]
                        points_y = chunk.y[mask]
                        points_z = chunk.z[mask]
                    else:
                        points_x = chunk.x
                        points_y = chunk.y
                        points_z = chunk.z
                    
                    if len(points_x) == 0: continue

                    dataset = np.column_stack((points_x, points_y, points_z))
                    points = np.rec.fromrecords(dataset, names='x, y, z')
                    yield points
        except Exception as e:
            utils.echo_warning_msg(f'Could not read points from local lasfile {self.fn}: {e}')

            
    def _yield_points_copc(self):
        """Yield points from remote COPC using chunk streaming."""
        
        if not self._parse_copc_metadata():
            return

        if self.verbose:
            utils.echo_msg(f"Traversing COPC hierarchy for {self.fn}...")

        self.chunk_table = []
        self._traverse_copc_hierarchy(
            self.copc_info['root_hier_offset'], 
            self.copc_info['root_hier_len']
        )
        
        if self.verbose:
            utils.echo_msg(f"Fetching {len(self.chunk_table)} chunks...")

        for chunk in self.chunk_table:
            raw_bytes = self._fetch_bytes(chunk['offset'], chunk['length'])
            
            if raw_bytes:
                with io.BytesIO(raw_bytes) as f:
                    try:
                        las = lp.open(f)
                        las_data = las.read()
                        
                        if self.classes:
                            mask = np.isin(las_data.classification, self.classes)
                            points_x = np.array(las_data.x[mask])
                            points_y = np.array(las_data.y[mask])
                            points_z = np.array(las_data.z[mask])
                        else:
                            points_x = np.array(las_data.x)
                            points_y = np.array(las_data.y)
                            points_z = np.array(las_data.z)
                        
                        if len(points_x) == 0: continue
                        
                        dataset = np.column_stack((points_x, points_y, points_z))
                        points = np.rec.fromrecords(dataset, names='x, y, z')
                        yield points
                        
                    except Exception as e:
                        utils.echo_warning_msg(f"Failed to decompress COPC chunk: {e}")

### End
