#!/usr/bin/env python
##
import sys
from geomods import waffles

if __name__ == '__main__':
    '''adjust the x value in an xyz file by 360'''
    which_way = None
    try:
        infile = sys.argv[1]
    except:
        waffles.echo_error_msg('you must enter an infile; x360.py input')
        sys.exit(0)
        
    if len(sys.argv) > 2: which_way = sys.argv[2]
    with open(infile, 'r') as iob:
        for xyz in waffles.xyz_parse(iob):
            if which_way is None: x = xyz[0] - 360
            else: x = xyz[0] + 360
            waffles.xyz_line([x] + xyz[1:])

### End

