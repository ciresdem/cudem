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
from typing import Optional, List, Dict, Any

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
        """Generate a 'mini-grid' of the data (placeholder)."""
        
        raise NotImplementedError

    
    def generate(self):
        """Generate metadata using DatasetFactory (stub)."""
        
        if self.name is None:
            return self

        from cudem.datasets import DatasetFactory
        
        this_ds = DatasetFactory(mod=self.name)
        this_ds_mod = this_ds._acquire_module()        
        self.fmt = this_ds.mod_name
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
            ## Fallback to MBSystem parsing
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
            pass # Silent fail, logging would be better.

### End
