### htdp.py
##
## Copyright (c) 2022  Regents of the University of Colorado
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

import os
import sys
import numpy as np
from cudem import utils

## =============================================================================
## HTDP
## =============================================================================

class HTDP:

    _reference_frames = {
        6781: {'name': 'NAD_83(2011/CORS96/2007)',
               'description': '(North American plate fixed)',
               'htdp_id': 1},
        6321: {'name': 'NAD_83(PA11/PACP00)',
               'description': '(Pacific plate fixed)',
               'htdp_id': 2},
        6324: {'name': 'NAD_83(MA11/MARP00)',
               'description': '(Mariana plate fixed)',
               'htdp_id': 3},
        # 4979: {'name': 'WGS_84(original)',
        #        'description': '(NAD_83(2011) used)',
        #        'htdp_id': 4},
        7815: {'name': 'WGS_84(original)',
               'description': '(NAD_83(2011) used)',
               'htdp_id': 4},
        7816: {'name': 'WGS_84(original)',
               'description': '(NAD_83(2011) used)',
               'htdp_id': 4},
        6319: {'name': 'WGS_84(original)',
               'description': '(NAD_83(2011) used)',
               'htdp_id': 4},
        7656: {'name': 'WGS_84(G730)',
               'description': '(ITRF91 used)',
               'htdp_id': 5},
        7657: {'name': 'WGS_84(G730)',
               'description': '(ITRF91 used)',
               'htdp_id': 5},
        7658: {'name': 'WGS_84(G873)',
               'description': '(ITRF94 used)',
               'htdp_id': 6},
        7659: {'name': 'WGS_84(G873)',
               'description': '(ITRF94 used)',
               'htdp_id': 6},
        7660: {'name': 'WGS_84(G1150)',
               'description': '(ITRF2000 used)',
               'htdp_id': 7},
        7661: {'name': 'WGS_84(G1150)',
               'description': '(ITRF2000 used)',
               'htdp_id': 7},
        4979: {'name': 'WGS_84(G1674)',
               'description': '(ITRF2008 used)',
               'htdp_id': 8},
        7662: {'name': 'WGS_84(G1674)',
               'description': '(ITRF2008 used)',
               'htdp_id': 8},
        7663: {'name': 'WGS_84(G1674)',
               'description': '(ITRF2008 used)',
               'htdp_id': 8},
        7664: {'name': 'WGS_84(G1762)',
               'description': '(IGb08 used)',
               'htdp_id': 9},
        7665: {'name': 'WGS_84(G1762)',
               'description': '(IGb08 used)',
               'htdp_id': 9},
        7666: {'name': 'WGS_84(G2139)',
               'description': '(ITRF2014=IGS14=IGb14 used)',
               'htdp_id': 10},
        7667: {'name': 'WGS_84(G2139)',
               'description': '(ITRF2014=IGS14=IGb14 used)',
               'htdp_id': 10},
        4910: {'name': 'ITRF88',
               'description': '',
               'htdp_id': 11},
        4911: {'name': 'ITRF89',
               'description': '',
               'htdp_id': 12},
        7901: {'name': 'ITRF89',
               'description': '',
               'htdp_id': 12},
        7902: {'name': 'ITRF90',
               'description': '(PNEOS90/NEOS90)',
               'htdp_id': 13},
        7903: {'name': 'ITRF91',
               'description': '',
               'htdp_id': 14},
        7904: {'name': 'ITRF92',
               'description': '',
               'htdp_id': 15},
        7905: {'name': 'ITRF93',
               'description': '',
               'htdp_id': 16},
        7906: {'name': 'ITRF94',
               'description': '',
               'htdp_id': 17},
        7907: {'name': 'ITRF96',
               'description': '',
               'htdp_id': 18},
        7908: {'name': 'ITRF97',
               'description': 'IGS97',
               'htdp_id': 19},
        7909: {'name': 'ITRF2000',
               'description': 'IGS00/IGb00',
               'htdp_id': 20},
        7910: {'name': 'ITRF2005',
               'description': 'IGS05',
               'htdp_id': 21},
        7911: {'name': 'ITRF2008',
               'description': 'IGS08/IGb08',
               'htdp_id': 22},
        7912: {'name': 'ITRF2014',
               'description': 'IGS14/IGb14',
               'htdp_id': 23},
    }
    
    def __init__(self, htdp_bin='htdp', verbose=True):
        self.htdp_bin = htdp_bin
        self.verbose=verbose

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
        fd = open(filename)
        grid = np.zeros(shape)

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

        lat_axis = np.linspace(lat_start,lat_end,lat_steps)
        lon_axis = np.linspace(lon_start,lon_end,lon_steps)

        lon_list = []
        for i in range(lat_steps):
            lon_list.append(lon_axis)

        lon_band = np.vstack(lon_list)

        lat_list = []
        for i in range(lon_steps):
            lat_list.append(lat_axis)

        lat_band = np.column_stack(lat_list)
            
        return(np.array([lon_band,lat_band]))

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
        utils.run_cmd('{} < {}'.format(self.htdp_bin, htdp_control), verbose=self.verbose)


if __name__ == '__main__':

    in_epsg = 6781
    out_epsg = 7912

    htdp = HTDP()
    griddef = (144.5, 13.75, 144.75, 13.5, 212, 212)
    grid = htdp._new_create_grid(griddef)
    htdp._write_grid(grid, '_tmp_input.xyz')
    htdp._write_control('_tmp_control.txt', '_tmp_output.xyz', '_tmp_input.xyz', htdp._reference_frames[in_epsg]['htdp_id'], 2012.0, htdp._reference_frames[out_epsg]['htdp_id'], 2012.0)
    htdp.run('_tmp_control.txt')
    out_grid = htdp._read_grid('_tmp_output.xyz', (griddef[4],griddef[5]))
    print(out_grid)
        
### End
