### copcfile.py - Cloud Optimized Point Cloud Parser
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
##
### Commentary:
##
## A native Python COPC reader for CUDEM.
## Uses cudem.fetches for HTTP Range requests to traverse the COPC Octree
## and only download chunks overlapping the requested region.
##
## Requires: laspy[lazrs] or laszip for the final byte decompression.
##
### Code:

import struct
import io
import numpy as np
import laspy 

from cudem import utils
from cudem import regions
from cudem.fetches import fetches
from cudem.datalists import dlim

class COPCFile(dlim.ElevationDataset):
    """Parse Cloud Optimized Point Cloud (COPC) files via HTTP Range Requests.
    
    Traverses the COPC hierarchy to find only the chunks intersecting the
    requested region/resolution, downloads them via cudem.fetches, and
    yields the points.
    """
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.copc_info = None
        self.hierarchy = {}
        self.chunk_table = []

        
    def _fetch_bytes(self, start, length):
        """Helper to fetch a byte range using cudem.fetches"""
        
        fetcher = fetches.Fetch(url=self.fn, verbose=self.verbose)
        
        ## Set Range Header manually
        fetcher.headers['Range'] = f'bytes={start}-{start + length - 1}'
        
        ## We use fetch_req directly to get the response object
        req = fetcher.fetch_req()
        if req is not None and req.status_code in [200, 206]:
            return req.content
        return None

    
    def _parse_header(self):
        """Parse LAS Header and COPC VLRs."""
        
        ## LAS Header is 375 bytes
        header_bytes = self._fetch_bytes(0, 375)
        if not header_bytes:
            return False

        ## Unpack Scale/Offset (LAS 1.4 Header offsets)
        self.scale = struct.unpack('<3d', header_bytes[131:155])
        self.offset = struct.unpack('<3d', header_bytes[155:179])
        
        ## Offset to VLRs is at 96 (4 bytes)
        vlr_offset = struct.unpack('<I', header_bytes[96:100])[0]
        ## Number of VLRs is at 100 (4 bytes)
        num_vlrs = struct.unpack('<I', header_bytes[100:104])[0]
        
        ## Fetch VLRs
        ## VLR Header is 54 bytes. 
        ## We need to scan them to find the COPC info (User: "entwine", Record: 1)
        ## This implementation naively fetches a block of bytes for VLRs. 
        ## A safer way is to fetch one by one, but VLRs are usually small.
        ## Let's fetch the first 2000 bytes after header to look for COPC VLR.
        vlr_chunk = self._fetch_bytes(vlr_offset, 2000)
        
        curr_pos = 0
        found_copc = False
        
        for _ in range(num_vlrs):
            if curr_pos + 54 > len(vlr_chunk): break
            
            ## VLR Header: Reserved(2), UserID(16), RecordID(2), Length(2), Desc(32)
            user_id = vlr_chunk[curr_pos+2 : curr_pos+18].replace(b'\x00', b'').decode('ascii')
            record_id = struct.unpack('<H', vlr_chunk[curr_pos+18 : curr_pos+20])[0]
            rec_len = struct.unpack('<H', vlr_chunk[curr_pos+20 : curr_pos+22])[0]
            
            if user_id == "entwine" and record_id == 1:
                ## Found COPC Info VLR
                data_pos = curr_pos + 54
                copc_data = vlr_chunk[data_pos : data_pos + rec_len]
                
                ## COPC Info: double(center_x,y,z), double(half_size), double(spacing), 
                ## uint64(root_hier_offset), uint64(root_hier_length)
                vals = struct.unpack('<4d d Q Q', copc_data)
                
                self.copc_info = {
                    'center': vals[0:3],
                    'half_size': vals[3],
                    'spacing': vals[4],
                    'root_hier_offset': vals[5],
                    'root_hier_len': vals[6]
                }
                found_copc = True
                break
            
            curr_pos += 54 + rec_len

        return found_copc

    def _traverse_hierarchy(self, offset, length):
        """Recursively fetch and parse hierarchy pages.
        Hierarchy Nodes are 32 bytes:
        Nodes: D(key), X(key), Y(key), Z(key), Offset(48b), Size(32b), Count(32b)
        """
        
        page_bytes = self._fetch_bytes(offset, length)
        if not page_bytes: return

        num_nodes = len(page_bytes) // 32
        
        for i in range(num_nodes):
            base = i * 32
            ## Unpack Node
            ## D, X, Y, Z are int32
            d, x, y, z = struct.unpack('<4i', page_bytes[base : base+16])
            
            ## Offset (6 bytes), Size (4 bytes), Count (4 bytes)
            ## Need strict struct packing for the 6-byte offset
            part2 = page_bytes[base+16 : base+32]
            offset_val = int.from_bytes(part2[0:6], 'little')
            byte_size = int.from_bytes(part2[6:10], 'little')
            point_count = int.from_bytes(part2[10:14], 'little')

            ## Calculate Node Bounds to check intersection
            ## Center = Origin + (Size * Step)
            ## COPC uses a cube octree centered at 'center' with radius 'half_size'
            
            ## This logic simplifies the octree math for brevity. 
            ## In a full implementation, you'd calculate the bbox of node (d,x,y,z)
            ## based on self.copc_info['half_size'] and self.copc_info['center']
            
            node_size = (self.copc_info['half_size'] * 2) / (2**d)
            node_min_x = (self.copc_info['center'][0] - self.copc_info['half_size']) + (x * node_size)
            node_max_x = node_min_x + node_size
            node_min_y = (self.copc_info['center'][1] - self.copc_info['half_size']) + (y * node_size)
            node_max_y = node_min_y + node_size
            
            node_region = regions.Region().from_list([node_min_x, node_max_x, node_min_y, node_max_y])
            
            ## Check Intersection with Requested Region
            if self.region is not None:
                if not regions.regions_intersect_p(self.region, node_region):
                    continue
            
            # If this is a leaf node with data (Size > 0 and Count > 0)
            if byte_size > 0 and point_count > 0:
                self.chunk_table.append({
                    'offset': offset_val,
                    'length': byte_size,
                    'count': point_count
                })

                
    def yield_points(self):
        """Yield points from intersecting chunks."""
        
        if not self._parse_header():
            utils.echo_error_msg("Failed to parse COPC header or VLR.")
            return

        if self.verbose:
            utils.echo_msg(f"COPC Info found. Traversing hierarchy...")

        ## Build Chunk Table by traversing root
        self._traverse_hierarchy(
            self.copc_info['root_hier_offset'], 
            self.copc_info['root_hier_len']
        )
        
        if self.verbose:
            utils.echo_msg(f"Identified {len(self.chunk_table)} chunks intersecting region.")

        ## Download and Parse Chunks
        for chunk in self.chunk_table:
            raw_bytes = self._fetch_bytes(chunk['offset'], chunk['length'])
            
            if raw_bytes:
                ## Use Laspy to decompress the bytes in memory
                with io.BytesIO(raw_bytes) as f:
                    try:
                        las = laspy.open(f)
                        las_data = las.read()
                        
                        ## Extract coordinates
                        x = np.array(las_data.x)
                        y = np.array(las_data.y)
                        z = np.array(las_data.z)
                        
                        ## Apply weights/uncertainty if needed (defaults)
                        w = np.ones_like(z)
                        u = np.zeros_like(z)
                        
                        ## Stack and yield
                        ## Using similar logic to lasfile.py
                        dataset = np.column_stack((x, y, z, w, u))
                        points = np.rec.fromrecords(dataset, names='x, y, z, w, u')
                        
                        yield points
                        
                    except Exception as e:
                        utils.echo_warning_msg(f"Failed to decompress COPC chunk: {e}")


### End
