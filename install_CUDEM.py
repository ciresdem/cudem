#!/usr/bin/env python

import os

try:
    from osgeo import gdal
except ImportError:
    raise (""" ERROR: Could not find the GDAL/OGR Python library bindings. 
               On Debian based systems you can install it with this command:
               apt install python-gdal
               On Redhat based systems (fedora) you can install it with this command:
               dnf install gdal gdal-devel python3-gdal 
               On windows, either use `OSGEO4W` or install GDAL via ...""") 

clone_and_install = False
pull_and_install = False
pip_install = True

cmd_exists = lambda x: any(os.access(os.path.join(path, x), os.X_OK) for path in os.environ['PATH'].split(os.pathsep))

if not cmd_exists('gmt'):
    raise(""" WARNING: Could not find GMT on the system.
              Some functionality of CUDEM will be unavailable without it.
              On Debian based systems you can install it with this command:
              apt install GMT
              On Redhat based systems (fedora/centos) you can install it with this command:
              dnf install GMT""")

if not cmd_exists('mblist'):
    raise(""" WARNING: Could not find MB-System on the system.
              Some functionality of CUDEM will be unavailable without it.
              Get the latest from ...""")

if not cmd_exists('las2txt'):
    raise(""" WARNING: Could not find LASTools on the system.
              Some functionality of CUDEM will be unavailable without it.
              Get the latest from ...""")

if clone_and_install:
    os.system('git clone https://github.com/ciresdem/cudem.git')
    os.system('pip install --user --upgrade ./cudem')    
elif pull_and_install:
    os.system('git pull')
    os.system('pip install --user --upgrade ./')
elif pip_install:
    os.system('pip install --user --upgrade ./')
