### htdp.py
##
## Copyright (c) 2022 - 2024  Regents of the University of Colorado
##
## htdp.py is part of CUDEM
##
## Permission is hereby granted, free of charge, to any person obtaining a copy 
## of this software and associated documentation files (the "Software"), to deal 
## in the Software without restriction, including without limitation the rights 
## to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies 
## of the Software, and to permit persons to whom the Software is furnished to do so, 
## subject to the following conditions:
##
## The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
##
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
## INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR 
## PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE 
## FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
## ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
##
### Commentary:
##
## wrapper and functions for using the htdp program
## https://geodesy.noaa.gov/TOOLS/Htdp/Htdp.shtml
##
## version 3.3.0
##  1...NAD_83(2011/CORS96/2007)  (North American plate fixed) 
##  2...NAD_83(PA11/PACP00)       (Pacific plate fixed) 
##  3...NAD_83(MA11/MARP00)       (Mariana plate fixed)                                                  
##  5...WGS_84(original) (NAD_83(2011) used) 15...ITRF91 
##  6...WGS_84(G730) (ITRF91 used)           16...ITRF92 
##  7...WGS_84(G873) (ITRF94 used)           17...ITRF93 
##  8...WGS_84(G1150) (ITRF2000 used)        18...ITRF94 
##  9...WGS_84(G1674) (ITRF2008 used)        19...ITRF96 
## 10...WGS_84(G1762) (IGb08 used)           20...ITRF97 or IGS97
## 11...SIO/MIT_92 (ITRF91 used)    21...ITRF2000 or IGS00/IGb00
## 12...ITRF88                      22...ITRF2005 or IGS05 
## 13...ITRF89                      23...ITRF2008 or IGS08/IGb08
## 14...ITRF90 or (PNEOS90/NEOS90)  24...ITRF2014 or IGS14/IGb14
##
## version 3.4.0
##  1...NAD_83(2011/CORS96/2007)  North America plate fixed     
##  2...NAD_83(PA11/PACP00)       Pacific plate fixed           
##  3...NAD_83(MA11/MARP00)       Mariana plate fixed           
                                                             
##  4...WGS84 original (Transit)                                
##  5...WGS84(G730)   ITRF91 used                               
##  6...WGS84(G873)   ITRF94=ITRF96=ITRF97 used                 
##  7...WGS84(G1150)  ITRF2000=IGS00=IGb00 used                 
##  8...WGS84(G1674)  ITRF2008=IGS08=IGb08 used                 
##  9...WGS84(G1762)  ITRF2008=IGS08=IGb08 used                 
## 10...WGS84(G2139)  ITRF2014=IGS14=IGb14 used                                                                             
## 11...ITRF88                      18...ITRF96 (=ITRF94=ITRF97)
## 12...ITRF89                      19...ITRF97 (=ITRF94=ITRF96)
## 13...ITRF90 (or PNEOS90/NEOS90)  20...ITRF2000 or IGS00/IGb00
## 14...ITRF91 (or SIO/MIT_92)      21...ITRF2005 or IGS05      
## 15...ITRF92                      22...ITRF2008 or IGS08/IGb08
## 16...ITRF93                      23...ITRF2014 or IGS14/IGb14
## 17...ITRF94 (=ITRF96=ITRF97) 
### Code:

import sys
import numpy as np
from cudem import utils

## =============================================================================
## HTDP
## =============================================================================

class HTDP:
    def __init__(self, htdp_bin='htdp', verbose=True):
        self.htdp_bin = htdp_bin
        self.verbose=verbose

        if utils.config_check()['HTDP'] is None:
            utils.echo_error_msg('you must have HTDP installed to perform vertical transformations')
        
    def _next_point(self, fd):
        line = fd.readline().strip()
        while line != '' and line.find('PNT_') == -1:
            line = fd.readline()

        if line == '':
            return None

        line = line.strip()
        name_tokens = line.split('_')
        tokens = line.split()
        lat_dst = float(tokens[0])
        lon_dst = float(tokens[1])
        eht_dst = float(tokens[2])
        return((int(name_tokens[1].strip('"')),
                int(name_tokens[2].strip('"')),
                lat_dst, lon_dst, eht_dst))

    def _read_grid(self, filename, shape):
        """read grid created by `_create_grid`"""
        
        fd = open(filename)
        grid = np.zeros(shape)

        if self.verbose:
            for i in range(5):
                utils.echo_msg(fd.readline().rstrip())

        points_found = 0
        ptuple = self._next_point(fd)
        while ptuple != None:
            grid[ptuple[1],ptuple[0]] = ptuple[4]
            points_found += 1
            ptuple = self._next_point(fd)

        if points_found < shape[0] * shape[1]:
            print('points found:   ', points_found)
            print('points expected:', shape[0] * shape[1])
            sys.exit(1)

        return(grid)
        
    def _new_create_grid(self, griddef):
        """This function creates a regular grid of lat/long values with one
        "band" for latitude, and one for longitude.
        """

        lon_start = -1 * griddef[0]
        lon_end = -1 * griddef[2]

        #lon_start = griddef[0]
        #lon_end = griddef[2]

        lat_start = griddef[1]
        lat_end = griddef[3]

        lon_steps = int(griddef[4])
        lat_steps = int(griddef[5])

        lon_axis = np.linspace(lon_start,lon_end,lon_steps)
        lat_axis = np.linspace(lat_start,lat_end,lat_steps)

        lon_list = []
        for i in range(lat_steps):
            lon_list.append(lon_axis)

        lon_band = np.vstack(lon_list)
        
        lat_list = []
        for i in range(lon_steps):
            lat_list.append(lat_axis)

        lat_band = np.column_stack(lat_list)

        t = np.array([lon_band, lat_band])
        return(t)

    def _write_grid(self, grid, out_filename):
        """This function writes a grid out in form suitable to use as input to the
        # htdp program.
        """
        
        fd_out = open(out_filename, 'w')

        for i in range(grid.shape[2]):
            for j in range(grid.shape[1]):
                fd_out.write('{} {} 0 "PNT_{}_{}"\n'.format(grid[1,j,i], grid[0,j,i], i, j))
                
        fd_out.close()
        
    def _write_control(self, control_fn, out_grid_fn,
                      in_grid_fn, src_crs_id, src_crs_date,
                      dst_crs_id, dst_crs_date):
        
        # start_date, end_date should be something like "2011.0"
        control_template = """
4
{out_grid}
{src_id}
{dst_id}
2
{src_date}
2
{dst_date}
3
{in_grid}
0
0
""".format(out_grid=out_grid_fn, src_id=src_crs_id,
           dst_id=dst_crs_id, src_date=src_crs_date,
           dst_date=dst_crs_date, in_grid=in_grid_fn)
        
        open(control_fn,'w').write(control_template)
        
    def run(self, htdp_control):

        if utils.config_check()['HTDP'] is not None:
            if utils.config_check()['platform'] == 'win32':
                os.system('cat {} | {}'.format(htdp_control, self.htdp_bin))
            else:
                utils.run_cmd('{} < {}'.format(self.htdp_bin, htdp_control), verbose=self.verbose)
        
### End
