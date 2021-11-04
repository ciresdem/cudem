#!/usr/bin/env python

# $Id: grd2mesh.py,v 1.1 \$

vcid = "$Id: grd2mesh.py,v 1.2"

import getopt, os, sys, shutil, glob
import numpy as np
from osgeo import gdal

def cleanup(files, keep):
    for i in files:
        if keep:
            for j in glob.glob(i):
                shutil.move(j, "tmp_grd2mesh")
        else:
            for j in glob.glob(i):
                os.remove(j)

def proc_cmds(cmds):
    for i in cmds:
        print(i)
        #os.system(i)
        
##---------------------------##
## Output Command Generation ##
##---------------------------##

def triangle_commands(grid,a):
    tmp0 = ('orig_len=$(wc -l clipped_bm.xyz | awk \'{print $1}\') \n\
    wc -l clipped_bm.xyz | awk \'{print $1,2,1,0}\' > header.txt \n\
    cat header.txt clipped_bm.xyz | awk \'{if (NR!=1) {print NR-1,$1,$2,$3} else {print}}\' > clipped_bm.node \n\
    triangle -q%d clipped_bm \n\
    node_len=$(wc -l clipped_bm.1.node | awk \'{print $1}\') \n\
    len_diff=$(echo ${node_len}-${orig_len}-2 | bc) \n\
    tail -${len_diff} clipped_bm.1.node | sed \'$d\' | awk \'{print $2,$3,$4}\' | gdal_query.py -d_format xyg %s >> clipped_bm.xyz \n\
    ' % (a,grid))
    proc_cmds([tmp0])
#grdhisteq tmp_curv_sub_elev.grd -Gtmp_nd.grd -C%d -Q \n\
def elevation_commands(grid,elevels,glevels):
    elev_cmds = ('grdmath -M %s ABS = tmp_curv_elev.grd \n\
    grdhisteq -C%d tmp_curv_elev.grd -Gtmp_curv_elnd.grd \n\
    grdmath -M tmp_curv_nd.grd tmp_curv_elnd.grd SUB = tmp_curv_sub_elev.grd \n\
    cp tmp_curv_sub_elev.grd tmp_nd.grd \n\
    ' % (grid,elevels))
    proc_cmds([elev_cmds])

def gradient_commands(grid,gradlevels,glevels):
    elev_cmds = ('grdgradient -M -D %s -Gtmp_grad.grd \n\
    grdhisteq tmp_grad.grd -Gtmp_grad_l%d.grd -C%d -Q \n\
    grdmath -M tmp_curv_l%d.grd tmp_curv_el%d.grd SUB = tmp_curv_l.grd \n\
    grdhisteq tmp_curv_l.grd -Gtmp_l%d.grd -C%d -Q \n\
    ' % (grid, glevels, gradlevels, glevels, glevels, glevels, glevels))
    proc_cmds([elev_cmds])
    
def mk_area_bounds(grid,inc,glevels):
    mk_area = ('grdsample %s -I%fc -Ggrid_sample.grd \n\
ymin=$(grdinfo grid_sample.grd | grep y_min | awk \'{print $3}\') \n\
ymax=$(grdinfo grid_sample.grd | grep y_min | awk \'{print $5}\') \n\
xmin=$(grdinfo grid_sample.grd | grep x_min | awk \'{print $3}\') \n\
xmax=$(grdinfo grid_sample.grd | grep x_min | awk \'{print $5}\') \n\
grd2xyz grid_sample.grd | grep -w -- $ymin > area_bounds.xyz \n\
grd2xyz grid_sample.grd | grep -w -- $ymax >> area_bounds.xyz \n\
grd2xyz grid_sample.grd | grep -w -- $xmin >> area_bounds.xyz \n\
grd2xyz grid_sample.grd | grep -w -- $xmax >> area_bounds.xyz \n\
cat area_bounds.xyz >> clipped_bm.xyz \n\
' % (grid,(inc*2)))
    proc_cmds([mk_area])

def cat_datalist(src_xyz):
    new_src = ("%s_src.xyz" % src_xyz.split(".")[0])
    mk_cat_srcs = ('\n\
cat %s | sed \'s/.*#.*//g\' | grep -v \'^$\' | awk \'{print $1}\' | rev | cut -f2- -d"/" | rev | sed \'s/$/\/*.xyz/g\' | awk \'{print $1,"\\\\"}\' | sed \'1i cat \\\\\' > "cat_%s" \n\
echo "| awk -F, \'{print \\$1,\\$2,\\$3}\' \\\\" >> cat_%s \n\
echo "> "%s >> cat_%s \n\
chmod +x $(echo cat_%s)\n\
./cat_%s \n\
' % (src_xyz, new_src, new_src, new_src, new_src, new_src, new_src))
    proc_cmds([mk_cat_srcs])
    return new_src

def generate_commands(grid,output,src_xyz,inc,minc,\
                          region,glevels,elevations_p,elevels,\
                          triangle_p,triangle_a,use_native_res_p,\
                          use_sph,keep_tmp,use_area_bounds):
    "Generate the Commands to make the mesh, given the user parameters"

    ## Obtain the NoDataValue of the input grid using GDAL
    ## We should change this to use GMT
    gd = gdal.Open(grid)
    nd = gd.GetRasterBand(1).GetNoDataValue()
    gd = None

    if nd is not None:
        if np.isnan(nd):
            nodata="NaN"
        else:
            nodata=str(nd)

    ## Discover the input xyz data
    ## This can be the xyz points from the input grid,
    ## or a given xyz file or an mbSystem style datalist
    if src_xyz == 'grd_xyz.xyz':
        mk_srcxyz = ('grd2xyz %s -V | grep -v %s > %s' % (grid,nodata,src_xyz))
        proc_cmds([mk_srcxyz])

    if "datalist" in src_xyz:
        src_xyz = cat_datalist(src_xyz)
        print(src_xyz)

    ## sample input grid to lowest resolution
    mk_sample_grd = ('grdsample %s -Gtmp_sample.grd -I%sc' %(grid,inc))
    proc_cmds([mk_sample_grd])
    grid="tmp_sample.grd"
    ## Generate the Curvature Grid, ABS and Normal Dist
    mk_curv_grd = ('grdmath -M %s CURV = tmp_curv.grd' % grid)
    mk_abs_curv = ('grdmath -M tmp_curv.grd ABS = tmp_curv_abs.grd')
    ## Assuming Normal Distribution
    mk_histeq = ('grdhisteq tmp_curv_abs.grd -Gtmp_curv_nd.grd -V -C%d' %(glevels))
    # mk_ceil_his = ('grdmath -M tmp_curv_nd.grd -V CEIL = tmp_curv_nd_c.grd')
    # mk_mv_his = ('mv tmp_curv_nd_c.grd tmp_curv_nd.grd')

    proc_cmds([mk_curv_grd, mk_abs_curv, mk_histeq])
    cleanup(["tmp_curv.grd", "tmp_curv_abs.grd"], keep_tmp)

    ## Process the Elevation grid, if user wants to.
    if elevations_p:
        elevation_commands(grid,elevels,glevels)
        cleanup(["tmp_curv_sub_elev.grd","tmp_curv_elev.grd","tmp_curv_nd.grd","tmp_curv_elnd.grd"], keep_tmp)
    else:
        mk_hist_grd = ('mv tmp_curv_nd.grd tmp_nd.grd')
        proc_cmds([mk_hist_grd])

    ## Query the regional domain grid with the input xyz data,
    ## returning an xyz data file which includes the normal dist
    ## value as the 4th column.
    #grdrs = ('grdsample tmp_nd.grd -I%dc -Gtmp_nd_rs.grd' % (inc*2))
    #proc_cmds([grdrs])
    #mk_grdquery = ('grdquery -g tmp_nd_rs.grd -i %s -r "xyzg" -d " " | grep -v nan | gmtselect -R%s > tmp_nd.xyzg' % (src_xyz, region))
    mk_gdalquery = ('gdal_query.py tmp_nd.grd %s -d_format xyzg -delimiter " " | grep -v nan | gmtselect -R%s > tmp_nd.xyzg' % (src_xyz, region))
    proc_cmds([mk_gdalquery])
    #cleanup(["tmp_nd.grd"], keep_tmp)

    ## First remove user-defined priority elevation values
    ## to be added to high-resolution xyz data

    mk_awkb = ('mv tmp_nd.xyzg tmp_rest%d.xyzg' % ((glevels*-1)+1))
    proc_cmds([mk_awkb])

    ## With the remainder, create new files for each of
    ## the desired normalized scores, with decreasing resolution.

    distr = [-8,-6,-4,-2,0,2,4,6,8]

    for i in range((glevels*-1)+1,glevels,3):  
        mk_awka = ('awk \'{if ($4 <= %d) {print $1,$2,$3 > \"grid_points%d.xyz\" } else {print $1,$2,$3,$4 > \"tmp_rest%d.xyzg\" }}\' tmp_rest%d.xyzg ' % (i,i,i+3,i))
        mk_awkb = ('awk \'{if ($4 != %d) {print $1,$2,$3,$4}}\' tmp_rest%d.xyzg > tmp_rest%d.xyzg' % (i,i,i+1))
        proc_cmds([mk_awka])

#     ## Make an xyz file of the combined regional domain grid
#     mk_hist_xyz = ('grd2xyz -V tmp_nd.grd > tc.xyz')
#     proc_cmds([mk_hist_xyz])
#     cleanup(["tmp_l%d.grd" % glevels], keep_tmp)

#     ## Parse the combined regional domain xyz file into it's various level peices
#     for i in range(0,glevels):
#         mk_awk = ('awk \'{if ($3==%d) {print}}\' tc.xyz > mag%d.xyz' % (i,(glevels-(i+1))))
#         proc_cmds([mk_awk])
#     cleanup(["tc.xyz"], keep_tmp)

#     ## Make the level xyz files into GMT mask grids
#     for i in range(0,glevels):
#         mk_grdmask = ('grdmask mag%d.xyz -R%s -I%fc -Gres_msk%d.grd -S%fc -V -N2/1/1' % (i,region,minc,i,minc))
#         proc_cmds([mk_grdmask])
#         cleanup(["mag%d.xyz" % i], keep_tmp)

#     ## Query the mask grids with the input xyz data
#     for i in range(0,glevels):
#         if i == 0:
#             mk_grdquery = ('grdquery -g res_msk%d.grd -i %s -d " " -r "xyzg" | grep -v nan > res_msk%d.xyzg' % (i,src_xyz,i))
#             mk_qawk0 = ('awk \'{if ($4==1) {print $1,$2,$3}}\' res_msk%d.xyzg > grid_points%d.xyz' % (i,i))
#         else:
#             mk_grdquery = ('grdquery -g res_msk%d.grd -i testarea-%d.xyz -d " " -r "xyzg" | grep -v nan > res_msk%d.xyzg' % (i,(i-1),i))
#             mk_qawk0 = ('awk \'{if ($4==1) {print $1,$2,$3}}\' res_msk%d.xyzg > grid_points%d.xyz' % (i,i))

#         mk_qawk1 = ('awk \'{if ($4!=1) {print $1,$2,$3}}\' res_msk%d.xyzg > testarea-%d.xyz' % (i,i))

#         proc_cmds([mk_grdquery, mk_qawk0, mk_qawk1])
#     cleanup(["res_msk*.grd", "res_msk*.xyzg", "testarea-*.xyz"], keep_tmp)

    ## Blockmedian the results of the mask query
    # zones=0
    # for i in range(glevels*-1,glevels,3):
    #     mk_zone = ('cat grid_points%d.xyz')
    #     zones+=1
    incs=inc
    gr = range((glevels*-1)+1,glevels)
    for i in gr[::-3]:
        if i == glevels*-1:
            if use_native_res_p:
                mk_blk1 = ('mv grid_points%d.xyz grid_block%d.xyz' % (i,i))
            else:
                mk_blk1 = ('blockmedian grid_points%d.xyz -I%fc -R%s > grid_block%d.xyz' % (i,inc,region,i))
        else:
            mk_blk1 = ('blockmedian grid_points%d.xyz -I%fc -R%s > grid_block%d.xyz' % (i,incs,region,i))
        proc_cmds([mk_blk1])
        incs=incs*3
    #cleanup(["grid_points*.xyz"], keep_tmp)

    ## Combine the blockmedian'd xyz data
    cat_blks = ('cat grid_block*.xyz > clipped_all.xyz')
    proc_cmds([cat_blks])
    cleanup(["grid_block*.xyz"], keep_tmp)

    ## Blockmedian the combined xyz data to the highest resolution
    if use_native_res_p:
        mk_clipped_bm = ('mv clipped_all.xyz clipped_bm.xyz')
        proc_cmds([mk_clipped_bm])
    else:
        mk_clipped_bm = ('blockmedian -R%s -I%fc clipped_all.xyz > clipped_bm.xyz' % (region,inc))
        proc_cmds([mk_clipped_bm])
        cleanup(["clipped_all.xyz"], keep_tmp)

    ## Generate the area bounding points if wanted
    if use_area_bounds:
        mk_area_bounds(grid,inc,glevels)
        cleanup(["grid_sample.grd", "area_bounds.xyz"], keep_tmp)

    ## If wanted, pre-triangulate the data with Schewchuck's Triangle, thus constraining the
    ## Delaunay Triangulation
    if triangle_p:
        triangle_commands(grid,triangle_a)
        cleanup(["header.txt","clipped_bm.node","clipped_bm.1.*"], keep_tmp)

    ##--------##
    ## Output ##
    ##--------##

    ## Triangulate the final xyz data with GMT
    if use_sph:
        mk_triangulate = ('sphtriangulate clipped_bm.xyz -M -Qd > %s_triangles.gmt' % output)
    else:
        mk_triangulate = ('triangulate clipped_bm.xyz -M -R%s -I%fc -G%s.grd > %s_triangles.gmt' % (region,inc,output,output))

    ## Translate the GMT output into Shapefiles
    mk_ogr2ogr = ('ogr2ogr %s_triangles.shp %s_triangles.gmt' % (output,output))
    mk_node_shp = ('xyz2shp.py -delimiter " " %s_nodes.shp clipped_bm.xyz ' % output)
    mk_node_xyz = ('mv clipped_bm.xyz %s_nodes.xyz' % output)
    proc_cmds([mk_triangulate, mk_ogr2ogr, mk_node_shp, mk_node_xyz])

##----------------##
## User Interface ##
##----------------##

def usage():
    print('usage: grd2mesh.py [options]')
    print('')
    print('Required Options:')
    print('-R, --region          The output region')
    print('-G, --grid            The Input GRD file')
    print('-O, --output          The Output basename')
    print('')
    print('Voluntary Options:')
    print('-X, --xyz             The xyz dataset to use in triangulation')
    print('-I, --increment       The output minimum spacing')
    print('-M, --mask-increment  The masking increment, used to classify values')
    print('-L, --levels          The number of levels')
    print('-E, --elevations      Account for elevation value in making levels,')
    print('                      use the given number of elevation levels')
    print('-T, --triangle        Use \'triangle\' to constrain triangle angle')
    print('-N, --native-res      Use \'native\' resolution for first class triangles')
    print('-B, --bounds          Generate and use area region xyz data for final triangulation')
    print('-S, --sphtriangulate  Use \'sphtriangulate\' instead of \'triangulate\'')
    print('-K, --keep-tmp        Do not remove temporary files, if set, will keep all')
    print('                      tmp files in temporary directory.')
    print('')
    print('Example:')
    print('--Generate triangulation from crm_s1.grd using triangles with')
    print('a 20 degree angles at 3 arc seconds minimum resolution using 22 levels:')
    print('grd2mesh.py -R-156.675/-155.5/20.07/20.725 -I3 -Gcrm_s1.grd -L22 -Ocrm1 -T20')
    print('')

def main():

    grid=None
    region=None
    inc=3
    minc=3
    glevels=3
    elevels=int(glevels/2)
    triangle_p=False
    triangle_a=20
    elevations_p=False
    src_xyz = 'grd_xyz.xyz'
    use_native_res_p=False
    use_sph=False
    keep_tmp=False
    use_area_bounds=False

    try:
        opts, args = getopt.getopt(sys.argv[1:], \
                                       "hO:R:G:I:M:L:X:T:E:NSKBv", \
                                       ["help", "output=", "region=", "grid=", \
                                            "increment=", "mask-increment=", "levels=", \
                                            "triangle=", "native-res", "sphtriangle", "keep", "area-bounds"])
    except getopt.GetoptError:
        # print help information and exit:
        #print(str(err)) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)
    output = None
    verbose = False
    for o, a in opts:
        if o == "-v":
            verbose = True
        elif o in ("-h", "--help"):
            usage()
            sys.exit()
        elif o in ("-R", "--region"):
            region = a
        elif o in ("-G", "--grid"):
            grid = a
        elif o in ("-I", "--increment"):
            inc = float(a)
        elif o in ("-M", "--mask-increment"):
            minc = float(a)
        elif o in ("-L", "--levels"):
            glevels = int(a)
        elif o in ("-E", "--elevations"):
            elevations_p = True
            elevels = int(a)
        elif o in ("-N", "--native-res"):
            use_native_res_p = True
        elif o in ("-O", "--output"):
            output = a
        elif o in ("-X", "--xyz"):
            src_xyz = a
        elif o in ("-T", "--triangle"):
            triangle_p = True
            triangle_a = float(a)
        elif o in ("-S", "--sphtriangulate"):
            use_sph = True
        elif o in ("-K", "--keep-tmp"):
            keep_tmp = True
        elif o in ("-B", "--area-bounds"):
            use_area_bounds = True
        else:
            assert False, "unhandled option"
    
    if grid is None or region is None:
        usage()
        sys.exit()
    else:
        if keep_tmp:
            os.mkdir("tmp_grd2mesh")
        generate_commands(grid,output,src_xyz,inc,minc,region,\
                          glevels,elevations_p,elevels,triangle_p,\
                          triangle_a,use_native_res_p,use_sph,keep_tmp,\
                          use_area_bounds)
    # ...

if __name__ == "__main__":
    main()
