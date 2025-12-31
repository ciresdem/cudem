### inf.py - DataLists IMproved
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
##
## inf.py is part of CUDEM
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
## Manage INF metadata files for datasets.
##
### Code:

import os
import json
import hashlib
import numpy as np
from typing import Optional, List, Dict, Any

from cudem import utils
from cudem import regions

class INF:
    """INF Files contain information about datasets."""
    
    def __init__(self,
                 name: Optional[str] = None,
                 file_hash: Optional[str] = None,
                 numpts: int = 0,
                 minmax: Optional[List[float]] = None,
                 wkt: Optional[str] = None,
                 fmt: Optional[str] = None,
                 src_srs: Optional[str] = None):
        
        self.name = name
        self.file_hash = file_hash
        self.hash = file_hash # Alias
        self.numpts = numpts
        self.minmax = minmax if minmax is not None else []
        self.wkt = wkt
        self.fmt = fmt
        self.format = fmt # Alias
        self.src_srs = src_srs
        self.mini_grid = None # Placeholder for 2D grid list

    def __str__(self):
        return f'<Dataset Info: {self.__dict__}>'

    def __repr__(self):
        return f'<Dataset Info: {self.__dict__}>'

    def generate_hash(self, fn: Optional[str] = None, sha1: bool = False) -> str:
        """Generate a hash of the source file."""
        
        target_file = fn if fn is not None else self.name
        
        if not target_file or not os.path.exists(target_file):
            self.file_hash = '0'
            return self.file_hash

        buf_size = 65536
        hasher = hashlib.sha1() if sha1 else hashlib.md5()
            
        try:
            with open(target_file, 'rb') as f:
                while True:
                    data = f.read(buf_size)
                    if not data:
                        break
                    hasher.update(data)

            self.file_hash = hasher.hexdigest()
        except Exception:
            self.file_hash = '0'

        return self.file_hash

    def generate_mini_grid(self, x_size: int = 10, y_size: int = 10):
        """
        Generate a 'mini-grid' of the data.
        
        Uses DatasetFactory to load the data and PointPixels to bin it into
        a coarse grid (default 10x10). This stored in the INF file for 
        quick spatial checking and visualization.
        """
        if self.name is None: return None

        # Local imports to avoid circular dependency
        from cudem import pointz
        from cudem.datalists.dlim import DatasetFactory

        # Ensure we have bounds
        if not self.minmax:
            return None

        # Create Region
        region = regions.Region().from_list(self.minmax)
        
        # Initialize Factory
        ds = DatasetFactory(mod=self.name, src_region=region)._acquire_module().initialize()
        if ds is None: return None

        # Accumulators
        running_sum = np.zeros((y_size, x_size))
        running_count = np.zeros((y_size, x_size))
        
        # PointPixels Processor
        pp = pointz.PointPixels(src_region=region, x_size=x_size, y_size=y_size, verbose=False)

        try:
            for points in ds.transform_and_yield_points():
                grid_arrays, _, _ = pp(points, mode='sums')
                if grid_arrays['z'] is not None:
                    running_sum += np.nan_to_num(grid_arrays['z'])
                    running_count += np.nan_to_num(grid_arrays['count'])

            # Finalize Mean
            with np.errstate(divide='ignore', invalid='ignore'):
                final_grid = running_sum / running_count
                
            # Convert to list for JSON serialization (NaN -> None)
            self.mini_grid = np.where(np.isnan(final_grid), None, final_grid).tolist()
            return self.mini_grid

        except Exception as e:
            utils.echo_warning_msg(f"Failed to generate mini-grid for {self.name}: {e}")
            return None

    def generate_block_mean(self, x_inc: float, y_inc: float = None, output_name: str = None):
        """
        Generate a block-mean XYZ file from the source data.
        
        Calculates the mean value of points within grid cells defined by 
        x_inc and y_inc and saves the result as an XYZ file.
        """
        if self.name is None: return None

        from cudem import pointz
        from cudem.datalists.dlim import DatasetFactory
        
        if y_inc is None: y_inc = x_inc
        
        if output_name is None:
            base = os.path.splitext(self.name)[0]
            output_name = f"{base}_blockmean.xyz"

        if not self.minmax: return None

        region = regions.Region().from_list(self.minmax)
        
        try:
            x_count, y_count, _ = region.geo_transform(x_inc=x_inc, y_inc=y_inc)
        except Exception:
            # Fallback if geotransform fails
            return None

        ds = DatasetFactory(mod=self.name, src_region=region)._acquire_module().initialize()
        if ds is None: return None

        pp = pointz.PointPixels(src_region=region, x_size=x_count, y_size=y_count, verbose=False)
        
        running_sum = np.zeros((y_count, x_count))
        running_count = np.zeros((y_count, x_count))
        running_x = np.zeros((y_count, x_count))
        running_y = np.zeros((y_count, x_count))

        try:
            for points in ds.transform_and_yield_points():
                grid_arrays, _, _ = pp(points, mode='sums')
                if grid_arrays['z'] is not None:
                    running_sum += np.nan_to_num(grid_arrays['z'])
                    running_count += np.nan_to_num(grid_arrays['count'])
                    running_x += np.nan_to_num(grid_arrays['x'])
                    running_y += np.nan_to_num(grid_arrays['y'])

            with np.errstate(divide='ignore', invalid='ignore'):
                final_z = running_sum / running_count
                final_x = running_x / running_count
                final_y = running_y / running_count

            valid_mask = (running_count > 0) & np.isfinite(final_z)
            
            out_x = final_x[valid_mask]
            out_y = final_y[valid_mask]
            out_z = final_z[valid_mask]
            
            with open(output_name, 'w') as f:
                for x, y, z in zip(out_x, out_y, out_z):
                    f.write(f"{x:.6f} {y:.6f} {z:.6f}\n")
            
            utils.echo_msg(f"Generated block-mean file: {output_name}")
            return output_name

        except Exception as e:
            utils.echo_warning_msg(f"Failed to generate block-mean for {self.name}: {e}")
            return None

    def generate(self, make_grid: bool = True, make_block_mean: bool = False, block_inc: float = None):
        """Generate metadata using DatasetFactory.
        
        Args:
            make_grid (bool): If True, generate the 10x10 mini-grid in the INF.
            make_block_mean (bool): If True, generate a companion block-mean XYZ file.
            block_inc (float): Resolution for block mean. If None, auto-calculates ~500px width.
        """
        
        if self.name is None:
            return self

        ## Local import to avoid circular dependency
        from cudem.datasets import DatasetFactory
        
        try:
            this_ds = DatasetFactory(mod=self.name)._acquire_module()
            if this_ds:
                # Populate basic info
                generated_info = this_ds.generate_inf()
                
                self.minmax = generated_info.minmax
                self.numpts = generated_info.numpts
                self.wkt = generated_info.wkt
                self.src_srs = generated_info.src_srs
                self.fmt = this_ds.data_format
                
                # Generate Mini Grid (Embedded in INF)
                if make_grid:
                    self.generate_mini_grid()
                    
                # Generate Block Mean (Sidecar File)
                if make_block_mean and self.minmax:
                    # Calculate default increment if not provided (Aim for ~500 pixels wide)
                    if block_inc is None:
                        width = self.minmax[1] - self.minmax[0]
                        block_inc = width / 500.0 if width > 0 else 0.001
                        
                    self.generate_block_mean(x_inc=block_inc)

        except Exception as e:
            utils.echo_warning_msg(f"Failed to generate INF metadata: {e}")
            pass
            
        return self

    
    def load_inf_file(self, inf_path: Optional[str] = None):
        """Load metadata from an existing INF file (JSON or MBSystem)."""
        
        if inf_path is None:
            return self
        
        if not os.path.exists(inf_path):
            return self

        data = {}
        
        ## Try JSON first
        try:
            with open(inf_path, 'r') as f:
                data = json.load(f)
        except (ValueError, json.JSONDecodeError):
            # Fallback to MBSystem parsing
            try:
                from cudem.datalists.mbsfile import MBSParser
                data = MBSParser(fn=inf_path).inf_parse().infos.__dict__
            except Exception as e:
                raise ValueError(f'Unable to read data from {inf_path} as JSON or MBSystem INF: {e}')

        ## Apply Loaded Data
        for key, val in data.items():
            if hasattr(self, key):
                setattr(self, key, val)

        return self

    def write_inf_file(self, inf_path: Optional[str] = None):
        """Write current metadata to a JSON INF file."""
        
        if inf_path is None:
            if self.name:
                inf_path = f'{self.name}.inf'
            else:
                return # Cannot write without a filename

        try:
            with open(inf_path, 'w') as outfile:
               json.dump(self.__dict__, outfile, indent=4)
        except Exception:
            pass 

### End

# ### inf.py - DataLists IMproved
# ##
# ## Copyright (c) 2010 - 2025 Regents of the University of Colorado
# ##
# ## inf.py is part of CUDEM
# ##
# ## Permission is hereby granted, free of charge, to any person obtaining a copy 
# ## of this software and associated documentation files (the "Software"), to deal 
# ## in the Software without restriction, including without limitation the rights 
# ## to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies 
# ## of the Software, and to permit persons to whom the Software is furnished to do so, 
# ## subject to the following conditions:
# ##
# ## The above copyright notice and this permission notice shall be included in all
# ## copies or substantial portions of the Software.
# ##
# ## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
# ## INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR 
# ## PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE 
# ## FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
# ## ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# ## SOFTWARE.
# ##
# ### Commentary:
# ##
# ## Manage INF metadata files for datasets.
# ##
# ### Code:

# import os
# import json
# import hashlib
# import numpy as np
# from typing import Optional, List, Dict, Any

# from cudem import utils
# from cudem import regions

# class INF:
#     """INF Files contain information about datasets."""
    
#     def __init__(self,
#                  name: Optional[str] = None,
#                  file_hash: Optional[str] = None,
#                  numpts: int = 0,
#                  minmax: Optional[List[float]] = None,
#                  wkt: Optional[str] = None,
#                  fmt: Optional[str] = None,
#                  src_srs: Optional[str] = None):
        
#         self.name = name
#         self.file_hash = file_hash
#         self.hash = file_hash # Alias
#         self.numpts = numpts
#         self.minmax = minmax if minmax is not None else []
#         self.wkt = wkt
#         self.fmt = fmt
#         self.format = fmt # Alias
#         self.src_srs = src_srs
#         self.mini_grid = None # Placeholder for the grid

        
#     def __str__(self):
#         return f'<Dataset Info: {self.__dict__}>'

    
#     def __repr__(self):
#         return f'<Dataset Info: {self.__dict__}>'

    
#     def generate_hash(self, fn: Optional[str] = None, sha1: bool = False) -> str:
#         """Generate a hash of the source file."""
        
#         target_file = fn if fn is not None else self.name
        
#         if not target_file or not os.path.exists(target_file):
#             self.file_hash = '0'
#             return self.file_hash

#         buf_size = 65536
#         hasher = hashlib.sha1() if sha1 else hashlib.md5()
            
#         try:
#             with open(target_file, 'rb') as f:
#                 while True:
#                     data = f.read(buf_size)
#                     if not data:
#                         break
#                     hasher.update(data)

#             self.file_hash = hasher.hexdigest()
#         except Exception:
#             self.file_hash = '0'

#         return self.file_hash

    
#     def generate_mini_grid(self, x_size: int = 10, y_size: int = 10):
#         """Generate a 'mini-grid' of the data.
        
#         Uses DatasetFactory to load the data and PointPixels to bin it into
#         a coarse grid (default 10x10). This is useful for quick spatial 
#         checking and visualization of data distribution.
        
#         Args:
#             x_size (int): Number of columns in the mini grid.
#             y_size (int): Number of rows in the mini grid.
            
#         Returns:
#             list: A 2D list representing the grid (Z values).
#         """
        
#         if self.name is None:
#             return None

#         ## Local imports to avoid circular dependency
#         from cudem import pointz
#         from cudem.datasets import DatasetFactory

#         ## Ensure we have bounds to define the grid
#         if not self.minmax:
#             self.generate() # Try to generate basic info first
#             if not self.minmax:
#                 return None

#         ## Create Region from minmax [xmin, xmax, ymin, ymax, ...]
#         region = regions.Region().from_list(self.minmax)
        
#         ## Initialize Factory
#         ds = DatasetFactory(mod=self.name, src_region=region)._acquire_module()
#         if ds is None:
#             return None

#         ## Initialize Accumulators for Mean Calculation
#         ## PointPixels with mode='sums' returns weighted sums and counts
#         running_sum = np.zeros((y_size, x_size))
#         running_count = np.zeros((y_size, x_size))
        
#         ## Initialize PointPixels Processor
#         ## We calculate the grid interactively chunk-by-chunk
#         pp = pointz.PointPixels(src_region=region, x_size=x_size, y_size=y_size, verbose=False)

#         try:
#             ## Iterate through data chunks
#             for points in ds.transform_and_yield_points():
#                 ## Process chunk into grid (getting Sums and Counts)
#                 ## mode='sums' returns Z as sum(z * w) and weight as sum(w)
#                 ## We assume weight=1 for simple mini-grid or use dataset weights
#                 grid_arrays, _, _ = pp(points, mode='sums')
                
#                 if grid_arrays['z'] is not None:
#                     ## Accumulate valid data
#                     ## PointPixels returns NaNs for empty cells, convert to 0 for accumulation
#                     chunk_z = np.nan_to_num(grid_arrays['z'])
#                     chunk_count = np.nan_to_num(grid_arrays['count'])
                    
#                     running_sum += chunk_z
#                     running_count += chunk_count

#             ## Finalize Mean: Sum / Count
#             with np.errstate(divide='ignore', invalid='ignore'):
#                 final_grid = running_sum / running_count
                
#             ## Replace Infs/NaNs (empty cells) with None for JSON compatibility
#             final_grid[~np.isfinite(final_grid)] = np.nan
            
#             ## Convert to list for JSON serialization
#             ## Using nan_to_num to None conversion logic if needed, 
#             ## but JSON dump handles None (null), not NaN usually without allow_nan=True
#             ## Here we convert numpy array to python list, replacing nans with None
#             self.mini_grid = np.where(np.isnan(final_grid), None, final_grid).tolist()
            
#             return self.mini_grid

#         except Exception as e:
#             utils.echo_warning_msg(f"Failed to generate mini-grid for {self.name}: {e}")
#             return None


#     def generate_block_mean(self, x_inc: float, y_inc: float = None, output_name: str = None):
#         """Generate a block-mean XYZ file from the source data.
        
#         This calculates the mean value of points within grid cells defined by 
#         x_inc and y_inc and saves the result as an XYZ file. This serves as 
#         a decimated or 'thinned' version of the original dataset.
        
#         Args:
#             x_inc (float): The X resolution (increment) for blocking.
#             y_inc (float): The Y resolution. Defaults to x_inc if None.
#             output_name (str): The output filename. Defaults to 
#                                '{self.name}_blockmean.xyz'.
        
#         Returns:
#             str: The path to the generated block-mean file, or None on failure.
#         """
        
#         if self.name is None:
#             return None

#         # Local imports
#         from cudem import pointz
#         from cudem.datasets import DatasetFactory

#         # Default Y increment
#         if y_inc is None:
#             y_inc = x_inc

#         # Determine Output Filename
#         if output_name is None:
#             # Strip existing extension and append _blockmean.xyz
#             base = os.path.splitext(self.name)[0]
#             output_name = f"{base}_blockmean.xyz"

#         # Ensure we have bounds
#         if not self.minmax:
#             self.generate()
#             if not self.minmax:
#                 return None

#         # Create Region
#         region = regions.Region().from_list(self.minmax)
        
#         # Calculate Grid Dimensions
#         try:
#             x_count, y_count, _ = region.geo_transform(x_inc=x_inc, y_inc=y_inc)
#         except Exception:
#             # Fallback if geotransform calc fails (e.g. region size < inc)
#             x_count = 10
#             y_count = 10

#         # Initialize Factory
#         ds = DatasetFactory(mod=self.name, src_region=region)._acquire_module()
#         if ds is None:
#             return None

#         # Initialize PointPixels Processor
#         pp = pointz.PointPixels(src_region=region, x_size=x_count, y_size=y_count, verbose=False)
        
#         # Accumulators
#         running_sum = np.zeros((y_count, x_count))
#         running_count = np.zeros((y_count, x_count))
#         running_x = np.zeros((y_count, x_count))
#         running_y = np.zeros((y_count, x_count))

#         try:
#             # Iterate and Accumulate
#             for points in ds.transform_and_yield_points():
#                 # We need 'sums' mode to aggregate properly across chunks
#                 grid_arrays, _, _ = pp(points, mode='sums')
                
#                 if grid_arrays['z'] is not None:
#                     # Accumulate Z sums
#                     running_sum += np.nan_to_num(grid_arrays['z'])
#                     running_count += np.nan_to_num(grid_arrays['count'])
                    
#                     # Accumulate X/Y sums for true block center positions
#                     # (Note: PointPixels returns aggregated X/Y sums in 'sums' mode)
#                     running_x += np.nan_to_num(grid_arrays['x'])
#                     running_y += np.nan_to_num(grid_arrays['y'])

#             # Finalize Means
#             with np.errstate(divide='ignore', invalid='ignore'):
#                 final_z = running_sum / running_count
#                 final_x = running_x / running_count
#                 final_y = running_y / running_count

#             # Mask invalid cells
#             valid_mask = (running_count > 0) & np.isfinite(final_z)
            
#             # Extract Valid Points
#             out_x = final_x[valid_mask]
#             out_y = final_y[valid_mask]
#             out_z = final_z[valid_mask]
            
#             # Write to File
#             with open(output_name, 'w') as f:
#                 for x, y, z in zip(out_x, out_y, out_z):
#                     f.write(f"{x:.6f} {y:.6f} {z:.6f}\n")
            
#             utils.echo_msg(f"Generated block-mean file: {output_name}")
#             return output_name

#         except Exception as e:
#             utils.echo_warning_msg(f"Failed to generate block-mean for {self.name}: {e}")
#             return None
        
        
#     # def generate(self):
#     #     """Generate metadata using DatasetFactory (stub)."""
#     #     if self.name is None:
#     #         return self

#     #     # Local import to avoid circular dependency
#     #     from cudem.datasets import DatasetFactory
        
#     #     try:
#     #         this_ds = DatasetFactory(mod=self.name)._acquire_module()
#     #         if this_ds:
#     #             ## In datasets.py, generate_inf() populates self.infos
#     #             ## So we update this object with the result.
#     #             generated_info = this_ds.generate_inf()
                
#     #             ## Update self attributes
#     #             self.minmax = generated_info.minmax
#     #             self.numpts = generated_info.numpts
#     #             self.wkt = generated_info.wkt
#     #             self.src_srs = generated_info.src_srs
#     #             self.fmt = this_ds.data_format
#     #     except Exception:
#     #         pass
            
#     #     return self

    
#     def load_inf_file(self, inf_path: Optional[str] = None):
#         """Load metadata from an existing INF file (JSON or MBSystem)."""
        
#         if inf_path is None:
#             return self
        
#         if not os.path.exists(inf_path):
#             return self

#         data = {}
        
#         ## Try JSON first
#         try:
#             with open(inf_path, 'r') as f:
#                 data = json.load(f)
#         except (ValueError, json.JSONDecodeError):
#             ## Fallback to MBSystem parsing
#             try:
#                 from cudem.datalists.mbsfile import MBSParser
#                 data = MBSParser(fn=inf_path).inf_parse().infos.__dict__
#             except Exception as e:
#                 raise ValueError(f'Unable to read data from {inf_path} as JSON or MBSystem INF: {e}')

#         ## Apply Loaded Data
#         for key, val in data.items():
#             if hasattr(self, key):
#                 setattr(self, key, val)

#         return self

    
#     def write_inf_file(self, inf_path: Optional[str] = None):
#         """Write current metadata to a JSON INF file."""
        
#         if inf_path is None:
#             if self.name:
#                 inf_path = f'{self.name}.inf'
#             else:
#                 return # Cannot write without a filename

#         try:
#             with open(inf_path, 'w') as outfile:
#                json.dump(self.__dict__, outfile, indent=4)
#         except Exception:
#             pass # Silent fail

# ### End
