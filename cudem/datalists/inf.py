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
###############################################################################
### Commentary:
##
### Examples:
##
### TODO:
##
### Code:

import os
import json
import hashlib

###############################################################################
## INF files and processing
###############################################################################
class INF:
    """INF Files contain information about datasets"""
    
    def __init__(self,
                 name=None,
                 file_hash=None,
                 numpts=0,
                 minmax=[],
                 wkt=None,
                 fmt=None,
                 src_srs=None):
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

    
    def generate_hash(self, fn=None, sha1=False):
        """generate a hash of the xyz-dataset source file"""

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

    
    def generate_mini_grid(self, x_size=10, y_size=10):
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
                        (f'CUDEMFactory: Unable to read data from {inf_fn} '
                         f'as mb-system inf, {e}')
                    )
            except:
                raise ValueError(
                    f'CUDEMFactory: Unable to read data from {inf_fn} as json'
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
        

### End
