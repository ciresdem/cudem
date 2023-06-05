#!/usr/bin/env python

import os
import sys

cmd_exists = lambda x: any(os.access(os.path.join(path, x), os.X_OK) for path in os.environ['PATH'].split(os.pathsep))

use_dnf = False
use_apt = False
want_mbsystem = False

if cmd_exists('dnf'):
    use_dnf = True
elif cmd_exists('apt'):
    use_apt = True

def install_htdp():
    htdp_url = 'https://geodesy.noaa.gov/TOOLS/Htdp/HTDP-download.zip'
    os.system('wget {}'.format(htdp_url))
    os.system('unzip HTDP-download.zip -d htdp_tmp')
    os.system('gfortran htdp_tmp/htdp.f')
    os.system('mv a.out ~/.local/bin/htdp')
    os.system('rm -rf htdp_tmp HTDP-download.zip')

def install_mbsystem():
    mbs_url = ''
    os.system('wget {}'.format(mbs_url))

def check_dependencies(install_missing=False):

    try:
        from osgeo import gdal
    except ImportError:
        if install_missing:
            if use_dnf:
                os.system('sudo dnf install gdal gdal-devel python3-gdal gdal-python-tools')
            elif use_apt:
                os.system('sudo apt install python-is-python3 python3-pip python3-gdal gdal-bin gdal-data')
        else:
            print(""" ERROR: Could not find the GDAL/OGR Python library bindings. 
            On Debian based systems you can install it with this command:
            apt install python-gdal
            On Redhat based systems (fedora) you can install it with this command:
            dnf install gdal gdal-devel python3-gdal 
            On windows, either use `OSGEO4W` or install GDAL via ...""")

    if not cmd_exists('gmt'):
        if install_missing:
            if use_dnf:
                os.system('sudo dnf install GMT GMT-devel GMT-doc')
            elif use_apt:
                os.system('sudo apt install gmt gmt-common gmt-dcw gmt-gshhg gmt-gshhg-full gmt-gshhg-low gmt-gshhg-high')
        else:
            print(""" WARNING: Could not find GMT on the system.
            Some functionality of CUDEM will be unavailable without it.
            On Debian based systems you can install it with this command:
            apt install GMT
            On Redhat based systems (fedora/centos) you can install it with this command:
            dnf install GMT""")

    if not cmd_exists('mblist') or want_mbsystem:
        if install_missing and use_dnf:
            if use_dnf:
                os.system('sudo dnf install g++ motif motif-devel fftw fftw-devel netcdf netcdf-devel proj proj-devel gdal-devel GMT GMT-devel boost boost-python3 glibc-devel mesa* xorg-x11-fonts* gcc-c++ libtirpc-devel')
            elif use_apt:
                os.system('sudo apt install g++ libgdal-dev libgmt-dev libmotif-common libmotif-dev libfftw3-dev fftw-dev libnetcdf-dev libproj-dev libboost-python-dev mesa* libtirpc-dev libglu1-mesa-dev freeglut3-dev mesa-common-dev')
                
            os.system('git clone https://github.com/dwcaress/MB-System.git')
            os.chdir('MB-System')
            os.system('./configure --prefix ~/.local')
            os.system('make')
            os.system('make install')
            os.system('cd ..')
            os.chdir('..')
            os.system('rm -rf MB-System')
        else:
            print(""" WARNING: Could not find MB-System on the system.
            Some functionality of CUDEM will be unavailable without it.
            Get the latest from ...""")

    if not cmd_exists('htdp'):
        if install_missing:
            if cmd_exists('wget') and cmd_exists('unzip'):
                
                if not cmd_exists('gfortran'):
                    if use_dnf:
                        os.system('sudo dnf install gfortran')
                    elif use_apt:
                        os.system('sudo apt install gfortran')

                if cmd_exists('gfortran'):
                    install_htdp()
        else:
            print(""" WARNING: Could not find HTDP on the system.
            Some functionality of CUDEM will be unavailable without it.
            Get the latest version: https://geodesy.noaa.gov/TOOLS/Htdp/HTDP-download.zip...""")

install_usage = """{cmd}: install CUDEM and dependencies

Options:
  --pull\t\tpull updates from git
  --dependencies\tattempt to install dependencies
Note: installing dependencies only supports RedHat based systems.

  --help\t\tPrint the usage text
""".format(cmd=os.path.basename(sys.argv[0]))
            
if __name__ == '__main__':
    clone = False
    pull = False
    install = True
    install_dep = False
    
    i = 1
    while i < len(sys.argv):
        arg = sys.argv[i]
        if arg == '--pull' or arg == '-p':
            pull = True
        elif arg == '--dependencies' or arg == '-d':
            install_dep = True
        elif arg == '--mbsystem' or arg == '-m':
            want_mbsystem = True

        elif arg == '--help' or arg == '-h':
            sys.stderr.write(install_usage)
            sys.exit(0)
            
        i = i + 1

    check_dependencies(install_missing = install_dep)
        
    if clone:
        os.system('git clone https://github.com/ciresdem/cudem.git')
    if pull:
        os.system('git pull')
    if install:
        if not cmd_exists('pip') and install_dep:
            os.system('sudo dnf install python3-pip')
            
        os.system('pip install --user --upgrade ./')
