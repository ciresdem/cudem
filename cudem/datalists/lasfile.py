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
###############################################################################
### Commentary:
##
### Examples:
##
### TODO:
##
### Code:

import os
import numpy as np
import laspy as lp

from cudem import utils
from cudem import regions
from cudem.datalists.dlim import ElevationDataset

class LASFile(ElevationDataset):
    """representing an LAS/LAZ dataset.

    Process LAS/LAZ lidar files using pylas.
    
    get_epsg - attempt to parse the EPSG from the LAS file header
    generate_inf - generate an inf file for the LAS data
    yield_xyz - yield the LAS data as xyz
    yield_array - yield the LAS data as an array must set the 
                  x_inc/y_inc in the super class
    
    -----------
    Parameters:
    
    classes (str): a list of classes to parse, being a string 
    with `/` seperator 
    """

    def __init__(self, classes='2/29/40', **kwargs):
        super().__init__(**kwargs)
        # list of lidar classes to retain
        self.classes = [int(x) for x in classes.split('/')] 
        if self.src_srs is None:
            self.src_srs = self.get_epsg()
            if self.src_srs is None:
                self.src_srs = self.infos.src_srs

                
    def valid_p(self, fmts=['scratch']):
        """check if self appears to be a valid dataset entry"""

        if self.fn is None: # and not self.fn.startswith('http'):
            return(False)
        else:
            if os.path.exists(self.fn) :
                if os.stat(self.fn).st_size == 0:
                    return(False)
            else:
                return(False)

            try:
                lp.open(self.fn)
            except:
                utils.echo_warning_msg(
                    f'{self.fn} could not be opened by the lasreader'
                )
                return(False)
                        
        return(True)

    
    def get_epsg(self):
        with lp.open(self.fn) as lasf:
            lasf_vlrs = lasf.header.vlrs
            for vlr in lasf_vlrs:
                if vlr.record_id == 2112:
                    src_srs = vlr.string
                    return(src_srs)
                
            return(None)

        
    def generate_inf(self):
        """generate an inf file for a lidar dataset."""
        
        with lp.open(self.fn) as lasf:
            self.infos.numpts = lasf.header.point_count
            this_region = regions.Region(
                xmin=lasf.header.x_min, xmax=lasf.header.x_max,
                ymin=lasf.header.y_min, ymax=lasf.header.y_max,
                zmin=lasf.header.z_min, zmax=lasf.header.z_max
            )
            self.infos.minmax = this_region.export_as_list(include_z=True)
            self.infos.wkt = this_region.export_as_wkt()

        #utils.echo_msg(self.get_epsg())
        self.infos.src_srs = self.src_srs \
            if self.src_srs is not None \
               else self.get_epsg()
        
        return(self.infos)

    
    def yield_points(self):
        with lp.open(self.fn) as lasf:
            try:
                for points in lasf.chunk_iterator(2_000_000):
                    points = points[(np.isin(points.classification, self.classes))]
                    dataset = np.column_stack((points.x, points.y, points.z))
                    points = np.rec.fromrecords(dataset, names='x, y, z')                    
                    yield(points)
                    
            except Exception as e:
                utils.echo_warning_msg(
                    f'could not read points from lasfile {self.fn}, {e}'
                )


### End
