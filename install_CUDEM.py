#!/usr/bin/env python

import os
import sys

try:
    from osgeo import gdal
except ImportError:
    raise (""" ERROR: Could not find the GDAL/OGR Python library bindings. 
               On Debian based systems you can install it with this command:
               apt install python-gdal
               On Redhat based systems (fedora) you can install it with this command:
               dnf install gdal gdal-devel python3-gdal 
               On windows, either use `OSGEO4W` or install GDAL via ...""") 


cmd_exists = lambda x: any(os.access(os.path.join(path, x), os.X_OK) for path in os.environ['PATH'].split(os.pathsep))

if not cmd_exists('gmt'):
    print(""" WARNING: Could not find GMT on the system.
              Some functionality of CUDEM will be unavailable without it.
              On Debian based systems you can install it with this command:
              apt install GMT
              On Redhat based systems (fedora/centos) you can install it with this command:
              dnf install GMT""")

if not cmd_exists('mblist'):
    print(""" WARNING: Could not find MB-System on the system.
              Some functionality of CUDEM will be unavailable without it.
              Get the latest from ...""")

if not cmd_exists('laszip'):
    print(""" WARNING: Could not find laszip on the system.
              Some functionality of CUDEM will be unavailable without it.
              Get the latest from ...""")

if __name__ == '__main__':
    clone = False
    pull = False
    install = True
    i = 1
    while i < len(sys.argv):
        arg = sys.argv[i]
        if arg == '--pull' or arg == '-p':
            pull = True
        i = i + 1

    if clone:
        os.system('git clone https://github.com/ciresdem/cudem.git')
    if pull:
        os.system('git pull')
    if install:
        os.system('pip install --user --upgrade ./')
