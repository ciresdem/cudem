#!/usr/bin/env python
### Code:
import os
import sys
import tempfile
import traceback
from xml.etree import ElementTree as et

import h5py
import numpy as np
from numpy import ma
import scipy.interpolate
#MRL: import scipy.spatial
from scipy import spatial
from osgeo import gdal
from tqdm import tqdm

DEBUGGING = False
if DEBUGGING:
    import matplotlib.pyplot as plt

count_file_ext = ".temp.cnt"

# plt.ion()

gmd = '{http://www.isotc211.org/2005/gmd}'
gco = '{http://www.isotc211.org/2005/gco}'

class BAGError(Exception):
    """ BAG class for exceptions"""

    def __init__(self, message, *args):
        if isinstance(message, Exception):
            msg = message.args[0] if len(message.args) > 0 else ''
        else:
            msg = message

        self.message = msg
        # allow users initialize misc. arguments as any other builtin Error
        Exception.__init__(self, message, *args)

class VR_Bag:
    def __init__(self, input_file_full_path):
        self.f = h5py.File(input_file_full_path)
        self.bag_root = self.f['/BAG_root']
        try:
            self.varres_metadata = self.bag_root['varres_metadata']
        except KeyError:
            raise BAGError("Could not find VR metadata, probably an SR BAG")
        self.varres_refinements = self.bag_root['varres_refinements']
        self.depth = self.varres_refinements[self.varres_refinements.dtype.names[0]][0]
        self.fill_value = self.varres_refinements.fillvalue[0]  # the depth fill value, not the stddev fill

        metadata = ''.join(map(bytes.decode, self.bag_root['metadata']))
        self.xml_root = et.fromstring(metadata)

    def get_valid_refinements(self):
        return self.varres_metadata[:, :, "dimensions_x"] > 0

    def get_res_x(self):
        return self.varres_metadata[:, :, "resolution_x"]

    def get_res_y(self):
        return self.varres_metadata[:, :, "resolution_y"]

    def get_max_depth(self):
        return self.depth[self.depth < self.fill_value].max()

    def get_min_depth(self):
        return self.depth.min()

    def read_meta_row_col(self):
        """
        Read meta_cols, meta_dx, meta_dx_uom info

        """
        meta_cols_arr = {}
        md_dimension = self.xml_root.findall(
            gmd + 'spatialRepresentationInfo/' + gmd + 'MD_Georectified/' + gmd + 'axisDimensionProperties/' + gmd + 'MD_Dimension')
        for elem in md_dimension:
            for sub_elem in elem:
                for sub_sub_elem in sub_elem:
                    # print (sub_sub_elem.tag, sub_sub_elem.text, sub_sub_elem.attrib)
                    if sub_sub_elem.text in ('column', "row"):
                        if sub_sub_elem.text == 'column':
                            meta_rc = 'meta_cols'
                            meta_dxy = 'meta_dx'
                            meta_uom = 'meta_dx_uom'
                        else:
                            meta_rc = 'meta_rows'
                            meta_dxy = 'meta_dy'
                            meta_uom = 'meta_dy_uom'

                        for ss_elem in elem:
                            if (ss_elem.tag == gmd + 'dimensionSize'):
                                for sss_elem in ss_elem:
                                    meta_cols = sss_elem.text
                                    meta_cols_arr[meta_rc] = meta_cols
                            if (ss_elem.tag == gmd + 'resolution'):
                                for sss_elem in ss_elem:
                                    meta_dx = sss_elem.text
                                    meta_cols_arr[meta_dxy] = meta_dx
                                    meta_dx_uom = sss_elem.get('uom')
                                    # print ('meta_dx_uom: ' + meta_dx_uom)
                                    meta_cols_arr[meta_uom] = meta_dx_uom
        return meta_cols_arr


    def read_meta_llx_lly(self):
        """
        Read meta_llx, meta_lly info

        """
        meta_llx_lly = {}
        corner_points = self.xml_root.findall(gmd + 'spatialRepresentationInfo/' + gmd + 'MD_Georectified/' + gmd + 'cornerPoints')
        for elem in corner_points:
            for sub_elem in elem:
                for sub_sub_elem in sub_elem:
                    # print (sub_sub_elem.tag, sub_sub_elem.text, sub_sub_elem.attrib)
                    corner_pt = str(sub_sub_elem.text).split()
                    # print (corner_pt)
                    meta_sw_pt = corner_pt[0].split(',')
                    # print (meta_sw_pt)
                    meta_llx = float(meta_sw_pt[0])
                    meta_lly = float(meta_sw_pt[1])
                    # print ('meta_llx: ' + str(meta_llx))
                    # print ('meta_lly: ' + str(meta_lly))
                    meta_llx_lly['meta_llx'] = meta_llx
                    meta_llx_lly['meta_lly'] = meta_lly
                    return meta_llx_lly


    def read_meta_horizontal_proj(self):
        """
        Read meta_horizontal_proj info

        """
        md_reference_system = self.xml_root.findall(
            gmd + 'referenceSystemInfo/' + gmd + 'MD_ReferenceSystem/' + gmd + 'referenceSystemIdentifier/' + gmd + 'RS_Identifier')
        for elem in md_reference_system:
            for sub_elem in elem:
                if (sub_elem.tag == gmd + 'code'):
                    for sub_sub_elem in sub_elem:
                        # print (sub_sub_elem.tag, sub_sub_elem.text, sub_sub_elem.attrib)
                        meta_horizontal_proj = sub_sub_elem.text.replace('\n', '')
                        # print ("meta_horizontal_proj:" + meta_horizontal_proj)
                        return meta_horizontal_proj


class Grid(object):
    '''Class to handle an ESRI grid
        NOTE: MUST UPDATE OLDER SCRIPTS TO INCLUDE BUFFER DISTANCE (4 nm = ~7500 m)'''

    def __init__(self, ll, nxy, cellSize, buffer_dist=0, allocate=True):
        """cellSize is either a single number for a square grid or a 2-tuple for x_size, ysize
        !!!buffer_dist is the size to expand a grid -- not sure this should be here!!!
        """
        self.allocate = allocate
        try:
            self.cell_size_x, self.cell_size_y = cellSize
        except:
            self.cell_size_x = cellSize
            self.cell_size_y = cellSize
        # Identify spatial characteristics of input grid raster
        self.minx = ll[0] - buffer_dist  # Desired grid cell size = 500 m ## BUFFER
        self.maxx = self.minx + nxy[0] * self.cell_size_x + buffer_dist
        self.miny = ll[1] - buffer_dist  # Desired grid cell size = 500 m ## BUFFER
        self.maxy = self.miny + nxy[1] * self.cell_size_y + buffer_dist

    @property
    def cell_size_x(self):
        return self._cell_size_y

    @cell_size_x.setter
    def cell_size_x(self, value):
        self._cell_size_y = value
        self.reset_bounds()

    @property
    def cell_size_y(self):
        return self._cell_size_y

    @cell_size_y.setter
    def cell_size_y(self, value):
        self._cell_size_y = value
        self.reset_bounds()

    @property
    def minx(self):
        return self._minx

    @minx.setter
    def minx(self, value):
        self._minx = value
        self.reset_bounds()

    @property
    def miny(self):
        return self._miny

    @miny.setter
    def miny(self, value):
        self._miny = value
        self.reset_bounds()

    @property
    def maxx(self):
        return self._maxx

    @maxx.setter
    def maxx(self, value):
        self._maxx = value
        self.reset_bounds()

    @property
    def maxy(self):
        return self._maxy

    @maxy.setter
    def maxy(self, value):
        self._maxy = value
        self.reset_bounds()

    @property
    def orig_x(self):
        return self.minx

    @property
    def orig_y(self):
        return self.miny

    @property
    def origin(self):
        return (self.orig_x, self.orig_y)

    def reset_bounds(self):
        try:
            self.x_edges = np.arange(self.orig_x, self.maxx + self.cell_size_x, self.cell_size_x)
            self.y_edges = np.arange(self.orig_y, self.maxy + self.cell_size_y, self.cell_size_y)
            if self.allocate:
                self.grid_val = self.zeros()
                self.grid_val.fill(np.nan)
                self.grid_count = self.zeros()
        except AttributeError:
            pass  # may not be set up yet

    # Define extents of Groundings Histogram based on extents of grid raster
    # Note: If Must keep origin at same grid cell node, set orig_x and orig_y as min coords of grid (bottom left)

    @property
    def orig(self):
        return np.array([self.orig_x, self.orig_y])  # Min x and y coordinate of grid raster; used to generate raster of Groundings

    @property
    def numx(self):
        return len(self.x_edges)-1

    @property
    def numy(self):
        return len(self.y_edges)-1

    @property
    def numcol(self):
        return self.numx

    @property
    def numrow(self):
        return self.numy

    def RowColFromXY(self, x, y):
        return self.row_index(y), self.col_index(x)
        raise Exception("This is backwards -- use ArrayIndicesFromXY")

    def row_index(self, y):
        return (y - self.orig_y) / self.cell_size_y

    def col_index(self, x):
        return (x - self.orig_x) / self.cell_size_x

    def ArrayIndicesFromXY(self, inputarray):
        '''Pass in a numpy array of XY values and returns the row,column indices as a numpy array'''
        output = np.array((inputarray - np.array(self.origin)) / np.array((self.cell_size_x, self.cell_size_y)),
                             dtype=np.int32)  # Translating polygon enveloppe indices to grid indices
        return output

    def zeros(self, dtype=np.float64):
        return np.zeros([self.numx, self.numy], dtype=dtype)

def plot(a, title=''):
    if DEBUGGING:
        plt.figure()
        plt.imshow(a, aspect='auto')
        plt.title(title)
        plt.draw()
        plt.show()

#MRL: add idw class
class Invdisttree():
    """ inverse-distance-weighted interpolation using KDTree:
    @Denis via https://stackoverflow.com/questions/3104781/inverse-distance-weighted-idw-interpolation-with-python
    https://creativecommons.org/licenses/by-nc-sa/3.0/
    """
    
    def __init__( self, X, z, leafsize=10, stat=0 ):
        assert len(X) == len(z), "len(X) %d != len(z) %d" % (len(X), len(z))
        self.tree = spatial.cKDTree( X, leafsize=leafsize )  # build the tree
        self.z = z
        self.stat = stat
        self.wn = 0
        self.wsum = None;

    def __call__(self, q, nnear=6, eps=0, p=1, dub=np.inf, weights=None):
        # nnear nearest neighbours of each query point --
        q = np.asarray(q)
        qdim = q.ndim
        if qdim == 1:
            q = np.array([q])
        if self.wsum is None:
            self.wsum = np.zeros(nnear)

        eps = float(eps)
        if eps is None: eps = .1
        self.distances, self.ix = self.tree.query(q, k=nnear, eps=eps, distance_upper_bound=dub)
        interpol = np.zeros((len(self.distances),) + np.shape(self.z[0]))
        jinterpol = 0
        for dist, ix in zip( self.distances, self.ix ):
            if np.any(np.isinf(dist)):
                wz = np.nan
            elif nnear == 1:
                wz = self.z[ix]
            elif dist[0] < 1e-10:
                wz = self.z[ix[0]]
            else:  # weight z s by 1/dist --
                w = 1 / dist**p
                if weights is not None:
                    w *= weights[ix]  # >= 0
                    
                w /= np.sum(w)
                wz = np.dot( w, self.z[ix] )
                if self.stat:
                    self.wn += 1
                    self.wsum += w
                    
            interpol[jinterpol] = wz
            jinterpol += 1
        return interpol if qdim > 1  else interpol[0]
        
# treat data as points, get the mask or treat as cell areas (combine with min/max/mean)
POINT, MASK, DATA_LOC_MASK, MIN, MAX, MEAN = range(6)
# when using cell areas determine how to treat the potential gaps from the individual cells to the supergrid boundary,
STRETCH, FILL, NOTHING = range(3)


def VRBag_to_TIF(input_file_full_path, dst_filename, sr_cell_size=None, mode=MIN, edge_option=STRETCH, use_blocks=True, nodata=np.nan):
    """

    Parameters
    ----------
    input_file_full_path : str
        Full path to BAG file
    dst_filename : str
        Full path to output TIFF file
    sr_cell_size : float
        cell size for resulting TIF.
        None will use the highest resolution cell size from the BAG file
    mode : int
        One of the defined modes of operation from the enumeration of POINT, DATA_LOC_MASK, DATA, MIN, MAX, MEAN.
        If mode is MIN, MAX or MEAN the BAG cells are treated as pixels that have area.
        If mode is POINT or MASK then the cells are treated as points that fall at the center of the cell.
        The difference between MASK and DATA_LOC_MASK is that the MASK shows all cells that could have had data
        while DATA_LOC_MASK only returns cells that will would data (Which doesn Glen want?).
    edge_option : int
        Supply a value from the enumeration STRETCH, FILL, NOTHING.
        Used if `mode` is one of the cell areas options, otherwise ignored.
    use_blocks : bool
        Boolean value to determine if the TIFF should use blocks in storage and if the
        program should load everything into memory (faster) or load the data piecewise (more robust)
    nodata : float
        The value to use in the resulting TIFF as no data.

    Returns
    -------

    """

    vr = VR_Bag(input_file_full_path)

    meta_rc = vr.read_meta_row_col()
    meta_llx_lly = vr.read_meta_llx_lly()
    meta_horizontal_proj = vr.read_meta_horizontal_proj()
    bag_supergrid_dx = float(meta_rc['meta_dx'])
    bag_supergrid_nx = int(meta_rc['meta_cols'])
    bag_supergrid_dy = float(meta_rc['meta_dy'])
    bag_supergrid_ny = int(meta_rc['meta_rows'])
    bag_llx = float(meta_llx_lly['meta_llx']) - bag_supergrid_dx/2.0  # @todo seems the llx is center of the supergridd cel?????
    bag_lly = float(meta_llx_lly['meta_lly']) - bag_supergrid_dy/2.0

    good_refinements = vr.get_valid_refinements()  # bool matrix of which refinements have data
    index2d = np.argwhere(good_refinements)
    # sort indices based on resolution areas
    res_area = vr.get_res_x()[good_refinements] * vr.get_res_y()[good_refinements]
    index2d = index2d[np.argsort(res_area)]

    # Adjust coordinates for the bag being able to extend outside the individual supergrids
    # Really it would be better to find the largest overlap and then fit to that exactly but this works for now.

    # So start from 1/2 supergrid below and run n+1 supergrids which ends up 1/2 above the end of the supergrids
    # sr_grid = Grid([bag_llx - bag_supergrid_dx/2.0, bag_lly - bag_supergrid_dy/2.0], [(bag_supergrid_nx+1) * bag_supergrid_dx, (bag_supergrid_ny+1) * bag_supergrid_dy], sr_cell_size, allocate=False)
    # Revise this to be the largest resolution refinement and buffer on either side for that amount.  Better than 1/2 supergrid at least
    possible_overlap_x = vr.get_res_x()[good_refinements].max()
    possible_overlap_y = vr.get_res_y()[good_refinements].max()
    if sr_cell_size is None or sr_cell_size == 0:
        sr_cell_size = min(vr.get_res_x()[good_refinements].min(), vr.get_res_y()[good_refinements].min())  # / 2.0
    # For deomstration purposes, make our grid line up with the GDAL grid.
    # GDAL uses the bottom left (BAG origin) and doesn't consider the possible overrun that BAG has
    # so extend by the number of cells necessary to have the BAG origin fall on a cell corner
    num_cells_to_extend_x = int(possible_overlap_x / sr_cell_size) + 1
    possible_overlap_x = num_cells_to_extend_x * sr_cell_size
    num_cells_to_extend_y = int(possible_overlap_y / sr_cell_size) + 1
    possible_overlap_y = num_cells_to_extend_y * sr_cell_size

    xcells = int((2 * possible_overlap_x + bag_supergrid_nx * bag_supergrid_dx) / sr_cell_size) + 1
    ycells = int((2 * possible_overlap_y + bag_supergrid_ny * bag_supergrid_dy)/ sr_cell_size) + 1
    sr_grid = Grid([bag_llx - possible_overlap_x, bag_lly - possible_overlap_y], [xcells, ycells], sr_cell_size, allocate=False)

    bagds = gdal.Open(input_file_full_path)
    if os.path.exists(dst_filename):
        os.remove(dst_filename)
    if os.path.exists(dst_filename + count_file_ext):
        os.remove(dst_filename + count_file_ext)
    fileformat = "GTiff"
    driver = gdal.GetDriverByName(fileformat)
    options = ["BLOCKXSIZE=256", "BLOCKYSIZE=256", "TILED=YES", "COMPRESS=LZW", "BIGTIFF=YES"]
    ds_val = driver.Create(dst_filename, sr_grid.numx, sr_grid.numy, bands=1, eType=gdal.GDT_Float32, options=options)
    ds_cnt = driver.Create(dst_filename+count_file_ext, sr_grid.numx, sr_grid.numy, bands=1, eType=gdal.GDT_Float32, options=options)

    ds_val.SetProjection(bagds.GetProjection())
    # bag_xform = bagds.GetGeoTransform()
    xform = (sr_grid.orig_x, sr_grid.cell_size_x, 0, sr_grid.maxy, 0, -sr_grid.cell_size_y)
    ds_val.SetGeoTransform(xform)
    r_val = ds_val.GetRasterBand(1)
    r_val.SetNoDataValue(0)
    ds_cnt.SetProjection(bagds.GetProjection())
    ds_cnt.SetGeoTransform(xform)
    r_cnt = ds_cnt.GetRasterBand(1)
    r_cnt.SetNoDataValue(0)
    bagds = None  # close the bag file

    if use_blocks is None:
        use_blocks = True if int(bag_supergrid_ny * bag_supergrid_dy) * int(bag_supergrid_nx * bag_supergrid_dx) * 4 > 1000000000 else False

    if not use_blocks:
        grid_val = r_val.ReadAsArray()
        grid_count = r_cnt.ReadAsArray()
        col_offset = 0
        row_offset = 0

    for i, j in tqdm(np.flipud(index2d), mininterval=.7):  # [:25]:
        # index_start,dimensions_x,dimensions_y,resolution_x,resolution_y,sw_corner_x,sw_corner_y = varres_metadata_np[i, j]
        index_start = vr.varres_metadata[i, j, "index"]
        dimensions_x = vr.varres_metadata[i, j, "dimensions_x"]
        dimensions_y = vr.varres_metadata[i, j, "dimensions_y"]
        resolution_x = vr.varres_metadata[i, j, "resolution_x"]
        resolution_y = vr.varres_metadata[i, j, "resolution_y"]
        sw_corner_x = vr.varres_metadata[i, j, "sw_corner_x"]
        sw_corner_y = vr.varres_metadata[i, j, "sw_corner_y"]

        # print('Processing Tile_{0}_{1}  row x col: {2}x{3}   res:{4:.2f}  cornerx:{5:.2f}  size:x={6:.2f} y={7:.2f}'.format(i, j, dimensions_x, dimensions_y, resolution_x, sw_corner_x, dimensions_x*resolution_x, dimensions_y*resolution_y))
        # @TODO check this for accuracy
        supergrid_x = j * bag_supergrid_dx
        supergrid_y = i * bag_supergrid_dy
        refinement_llx = bag_llx + supergrid_x + sw_corner_x - resolution_x / 2.0  # @TODO implies swcorner is to the center and not the exterior
        refinement_lly = bag_lly + supergrid_y + sw_corner_y - resolution_y / 2.0
        index_end = index_start + int(dimensions_x * dimensions_y)
        # big speedup to read the index range then get the depth array, so switching order of read.
        # tile = vr.depth[index_start:index_end].reshape(dimensions_y, dimensions_x)
        tile = vr.varres_refinements[:, index_start:index_end][vr.varres_refinements.dtype.names[0]].reshape(dimensions_y, dimensions_x)
        tile[tile == vr.fill_value] = np.nan

        # @TODO remember to pad the first/last row+col to remove the dead space "pane" between the refinement adn the supergrid edge
        # compute the edges of the grid rows/columns and their related indices in the final SR grid
        if mode in (POINT, MASK, DATA_LOC_MASK):
            # use the cell centers which means the start/end index will be the same
            yends = ystarts = refinement_lly + np.arange(dimensions_y) * resolution_y + resolution_y / 2.0
            xends = xstarts = refinement_llx + np.arange(dimensions_x) * resolution_x + resolution_x / 2.0
        else:
            # determine the edges of the cells so they can be mapped into the new geotiff overlapped coordinate
            ystarts = refinement_lly + np.arange(dimensions_y) * resolution_y
            yends = refinement_lly + (np.arange(dimensions_y) + 1) * resolution_y - .000001
            xstarts = refinement_llx + np.arange(dimensions_x) * resolution_x
            xends = refinement_llx + (np.arange(dimensions_x) + 1) * resolution_x - .000001

            if edge_option in (STRETCH, FILL):
                # FILL always stretches the fist+last row/col to the edge of the supercell
                # STRETCH only fills the gap if the refinement is close the the edge
                # -- thinking the user didn't mean to leave the gap inherent in the BAG format
                apply_stretch = edge_option == STRETCH and sw_corner_x <= resolution_x and sw_corner_y <= resolution_y
                if edge_option == FILL or apply_stretch:
                    # print(xstarts[0], xends[-1], ystarts[0], yends[-1])
                    xstarts[0] = min(xstarts[0], bag_llx + supergrid_x)
                    xends[-1] = max(xends[-1], bag_llx + supergrid_x + bag_supergrid_dx - .000001)
                    ystarts[0] = min(ystarts[0], bag_lly + supergrid_y)
                    yends[-1] = max(yends[-1], bag_lly + supergrid_y + bag_supergrid_dy - .000001)
                    # print(xstarts[0], xends[-1], ystarts[0], yends[-1])
                if edge_option == STRETCH and not apply_stretch:
                    print("not stretched??")

        # convert the BAG coordinates into geotiff pixel indices
        # also reverse the rows since bag is lowerleft and tif is upper left
        row_end_indices = (sr_grid.numrow - 1) - np.array(sr_grid.row_index(ystarts), int)
        row_start_indices = (sr_grid.numrow - 1) - np.array(sr_grid.row_index(yends), int)
        col_start_indices = np.array(sr_grid.col_index(xstarts), int)
        col_end_indices = np.array(sr_grid.col_index(xends), int)

        if DEBUGGING:
            row_min = int(min(row_start_indices.min(), row_end_indices.min()))
            col_min = int(min(col_start_indices.min(), col_end_indices.min()))
            row_max = int(max(row_start_indices.max(), row_end_indices.max())) + 1
            col_max = int(max(col_start_indices.max(), col_end_indices.max())) + 1
            if row_min<0 or row_max<0 or col_min<0 or col_max<0:
                print("something is less than zero")
        if use_blocks:  # read the part of the geotiff that we need and modify it then write it back after applying the refinement grid
            row_offset = int(min(row_start_indices.min(), row_end_indices.min()))
            col_offset = int(min(col_start_indices.min(), col_end_indices.min()))
            row_max = int(max(row_start_indices.max(), row_end_indices.max())) + 1
            col_max = int(max(col_start_indices.max(), col_end_indices.max())) + 1
            # grid_val = r_val.ReadAsArray(row_offset, col_offset, row_max - row_offset, col_max - col_offset)
            # grid_count = r_cnt.ReadAsArray(row_offset, col_offset, row_max - row_offset, col_max - col_offset)  # col_offset, row_offset, col_max-col_offset, row_max-row_offset
            grid_val = r_val.ReadAsArray(col_offset, row_offset, col_max - col_offset, row_max - row_offset)
            grid_count = r_cnt.ReadAsArray(col_offset, row_offset, col_max - col_offset, row_max - row_offset)  # col_offset, row_offset, col_max-col_offset, row_max-row_offset


        # iterate the refinement cells and place each into the final SR grid
        for tile_i in range(tile.shape[0]):
            for tile_j in range(tile.shape[1]):
                cur_val = tile[tile_i, tile_j]
                # If there is real data or we are creating the MASK (which is jsut a map of cell centers) then write a value into output
                if not np.isnan(cur_val) or mode is MASK:
                    col_row_slice = np.s_[row_start_indices[tile_i] - row_offset:row_end_indices[tile_i] + 1 - row_offset,
                                    col_start_indices[tile_j] - col_offset:col_end_indices[tile_j] + 1 - col_offset]
                    if mode == MEAN:
                        grid_val[col_row_slice] += cur_val
                    elif mode == MIN:
                        grid_val[col_row_slice][grid_count[col_row_slice] == 0] = cur_val
                        grid_val[col_row_slice] = np.minimum(grid_val[col_row_slice], cur_val)
                    elif mode == MAX:
                        grid_val[col_row_slice][grid_count[col_row_slice] == 0] = cur_val
                        grid_val[col_row_slice] = np.maximum(grid_val[col_row_slice], cur_val)
                    elif mode in (POINT, DATA_LOC_MASK):
                        grid_val[col_row_slice] = cur_val
                    grid_count[col_row_slice] += 1
        if DEBUGGING:
            row_min = int(min(row_start_indices.min(), row_end_indices.min()))
            col_min = int(min(col_start_indices.min(), col_end_indices.min()))
            row_max = int(max(row_start_indices.max(), row_end_indices.max())) + 1
            col_max = int(max(col_start_indices.max(), col_end_indices.max())) + 1
            plot(grid_val[row_min-row_offset:row_max-row_offset, col_min-col_offset:col_max-col_offset], "tile %d %d"%(i,j))

        # write out the current block to the TIFF file before we load/operate on the next one
        if use_blocks:
            r_val.WriteArray(grid_val, col_offset, row_offset)
            r_cnt.WriteArray(grid_count, col_offset, row_offset)
            if DEBUGGING:
                plot(r_val.ReadAsArray(), 'grid before averaging')

    # normalize the TIF if needed (divide value by count for MEAN)
    # then write the array into the raster band
    if use_blocks:
        # plot(r_val.ReadAsArray(), 'grid before averaging')

        block_sizes = r_val.GetBlockSize()
        row_block_size = block_sizes[0]
        col_block_size = block_sizes[1]
        row_size = r_val.XSize
        col_size = r_val.YSize
        r_val.SetNoDataValue(nodata)
        for ic in tqdm(range(0, col_size, col_block_size), mininterval=.7):
            if ic + col_block_size < col_size:
                cols = col_block_size
            else:
                cols = col_size - ic
            #JDV for ir in tqdm(range(0, row_size, row_block_size), mininterval=.7):
            for ir in range(0, row_size, row_block_size):
                if ir + row_block_size < row_size:
                    rows = row_block_size
                else:
                    rows = row_size - ir
                data = r_val.ReadAsArray(ir, ic, rows, cols)
                cnt = r_cnt.ReadAsArray(ir, ic, rows, cols)
                if mode == MEAN:
                    res = data/cnt
                    res[cnt == 0] = nodata
                elif mode in (MIN, MAX, POINT):
                    res = data
                    res[cnt == 0] = nodata
                elif mode in (MASK, DATA_LOC_MASK):
                    res = data
                    res[cnt == 0] = nodata
                    res[cnt != 0] = 1
                r_val.WriteArray(res, ir, ic)
    else:
        if mode == MEAN:
            grid_val = grid_val / grid_count
            grid_val[grid_count == 0] = nodata
        elif mode in (MIN, MAX, POINT):
            grid_val[grid_count == 0] = nodata
        elif mode == MASK:
            grid_val[grid_count == 0] = nodata
            grid_val[grid_count != 0] = 1
        # print("min {}, max {}".format(np.nanmin(grid_val), np.nanmax(grid_val)))
        # grid_val = np.flipud(grid_val)
        r_val.SetNoDataValue(nodata)
        r_val.WriteArray(grid_val)
        # plot(r_val.ReadAsArray(), 'read back from tiff')

    # if use_blocks:
    #     grid_val = r_val.ReadAsArray()  # for testing read the whole thing -- this breaks on a super large file of course
    #     grid_count = r_cnt.ReadAsArray()  # for testing read the whole thing -- this breaks on a super large file of course
    #
    # plt.gca().invert_yaxis()
    # plot(grid_val, 'SR grid at res {}'.format(sr_cell_size, {MIN: "Min", MAX: "Max", MEAN: "Mean"}[mode]))
    #
    # plot(grid_count, "count")

    # a = np.fromfunction(lambda i, j: i + j, a.shape, dtype=np.float32)  # add a gradient

    ds_cnt = None  #close the count file so it can be deleted
    os.remove(dst_filename+count_file_ext)
    return bag_supergrid_dx, bag_supergrid_dy, sr_cell_size

def interpolate_vr_bag(input_file_full_path, dst_filename, tmp_dir, sr_cell_size=None, method='linear', use_blocks=True, nodata=np.nan):
    """ Interpolation scheme
    Create the POINT version of the TIFF with only data at precise points of VR BAG
    Load in blocks with enough buffer around the outside (nominally 3x3 supergrids with 1 supergrid buffer)
        run scipy.interpolate.griddata on the block (use linear as cubic causes odd peaks and valleys)
        copy the interior (3x3 supergrids without the buffer area) into the output TIFF

    Create the MIN (or MAX) version of the TIFF
    Load blocks of data and copy any NaNs from the MIN (cell based coverage) into the INTERP grid to remove erroneous interpolations,
    this essentially limits coverage to VR cells that were filled
    """
    fobj, point_filename = tempfile.mkstemp(".point.tif", dir=tmp_dir)
    os.close(fobj)
    fobj, min_filename = tempfile.mkstemp(".min.tif", dir=tmp_dir)
    os.close(fobj)
    # if not DEBUGGING:
    if 1:
        print('Creating point output...') #JDV
        dx, dy, cell_sz = VRBag_to_TIF(input_file_full_path, point_filename, sr_cell_size=sr_cell_size, mode=POINT, use_blocks=use_blocks, nodata=nodata)
        print (dx, dy, cell_sz)
        print('Creating min output...') #JDV
        VRBag_to_TIF(input_file_full_path, min_filename, sr_cell_size=sr_cell_size, mode=MIN, use_blocks=use_blocks, nodata=nodata)
    else:
        dx, dy, cell_sz = 46.879258, 46.879258, 1.0
        point_filename=r"/data/BAG/tmp/H13000.min.tif"
        min_filename=r"/data/BAG/tmp/H13000.point.tif"

    points_ds = gdal.Open(point_filename)
    points_band = points_ds.GetRasterBand(1)
    points_no_data = points_band.GetNoDataValue()
    coverage_ds = gdal.Open(min_filename)
    coverage_band = coverage_ds.GetRasterBand(1)
    coverage_no_data = coverage_band.GetNoDataValue()
    interp_ds = points_ds.GetDriver().Create(dst_filename, points_ds.RasterXSize, points_ds.RasterYSize, bands=1, eType=points_band.DataType,
                                             options=["BLOCKXSIZE=256", "BLOCKYSIZE=256", "TILED=YES", "COMPRESS=LZW", "BIGTIFF=YES"])
    interp_ds.SetProjection(points_ds.GetProjection())
    interp_ds.SetGeoTransform(points_ds.GetGeoTransform())
    interp_band = interp_ds.GetRasterBand(1)
    interp_band.SetNoDataValue(nodata)

    print ('Interpolating...') #JDV
    if use_blocks:
        pixels_per_supergrid = int(max(dx / cell_sz, dy / cell_sz)) + 1
        row_block_size = col_block_size = 3 * pixels_per_supergrid
        row_buffer_size = col_buffer_size = 1 * pixels_per_supergrid
        row_size = interp_band.XSize
        col_size = interp_band.YSize
        for ic in tqdm(range(0, col_size, col_block_size), mininterval=.7):
            # print(str(ic) + '/' + str(col_size))
            cols = col_block_size
            if ic + col_block_size > col_size:  # read a full set of data by offsetting the column index back a bit
                ic = col_size - cols
            col_buffer_lower = col_buffer_size if ic >= col_buffer_size else ic
            col_buffer_upper = col_buffer_size if col_size - (ic + col_block_size) >= col_buffer_size else col_size - (ic + col_block_size)
            read_cols = col_buffer_lower + cols + col_buffer_upper
            # JDV for ir in tqdm(range(0, row_size, row_block_size), mininterval=.7):
            for ir in range(0, row_size, row_block_size):
                
                rows = row_block_size
                if ir + row_block_size > row_size:
                    ir = row_size - rows
                row_buffer_lower = row_buffer_size if ir >= row_buffer_size else ir
                row_buffer_upper = row_buffer_size if row_size - (ir + row_block_size) >= row_buffer_size else row_size - (ir + row_block_size)
                read_rows = row_buffer_lower + rows + row_buffer_upper
                points_array = points_band.ReadAsArray(ir-row_buffer_lower, ic-col_buffer_lower, read_rows, read_cols)
                
                # Find the points that actually have data as N,2 array shape that can index the data arrays
                if np.isnan(points_no_data):
                    if np.all(~np.isnan(points_array)):
                        interp_band.WriteArray(points_array, ir, ic)
                        continue
                    else:
                        point_indices = np.nonzero(~np.isnan(points_array))
                else:
                    if np.all(points_array != points_no_data):
                        interp_band.WriteArray(points_array, ir, ic)
                        continue
                    point_indices = np.nonzero(points_array != points_no_data)

                # if there were any data points then do interpolation -- could be all empty space too which raises Exception in griddata
                if len(point_indices[0]):                    
                    # get the associated data values
                    point_values = points_array[point_indices]
                    # interpolate all the other points in the array
                    # (actually it's interpolating everywhere which is a waste of time where there is already data)
                    xi, yi = np.mgrid[row_buffer_lower:row_buffer_lower+row_block_size,
                                      col_buffer_lower:col_buffer_lower+col_block_size]
                    
                    #JDV: put the interpolation step into a try-except block. 
                    #Sometimes an exception gets raised if there aren't enough points to do the interpolation. In that case, we should skip that cell and continue.
                    try:
                        #MRL: add IDW interpolation method
                        if method == 'idw':                            
                            invdisttree = Invdisttree(np.transpose(point_indices), point_values, leafsize=10, stat=1)
                            interp_data = invdisttree(
                                np.vstack((xi.flatten(), yi.flatten())).T,
                                nnear=3,
                                eps=.1,
                                p=2,
                                dub=12,
                                weights=None
                            )
                            interp_data = np.reshape(interp_data, (col_block_size, row_block_size))
                        elif method == 'linear' or method == 'cubic':
                            interp_data = scipy.interpolate.griddata(
                                np.transpose(point_indices), point_values,
                                (xi, yi), method=method
                            )
                            
                        # mask based on the cell coverage found using the MIN mode                        
                        coverage_data = coverage_band.ReadAsArray(ir, ic, row_block_size, col_block_size)
                        interp_data[coverage_data == coverage_no_data] = nodata

                        # Write the data into the TIF on disk
                        interp_band.WriteArray(interp_data, ir, ic)
                    except:
                        print('An error occurred in interpolation, continuing...')
                        traceback.print_exc()
                        continue
        if DEBUGGING:
            points_array = points_band.ReadAsArray()
            interp_array = interp_band.ReadAsArray()
            plot(points_array)
            plot(interp_array)
    else:
        points_array = points_band.ReadAsArray()
        coverage_data = coverage_band.ReadAsArray()

        # Find the points that actually have data
        if np.isnan(points_no_data):
            point_indices = np.nonzero(~np.isnan(points_array))
        else:
            point_indices = np.nonzero(points_array!=points_no_data)
        # get the associated data values
        point_values = points_array[point_indices]
        # interpolate all the other points in the array (actually interpolating everywhere which is a waste of time where there is already data)
        xi, yi = np.mgrid[0:points_array.shape[0], 0:points_array.shape[1]]

        #MRL: add IDW interpolation method
        if method == 'idw':
            invdisttree = Invdisttree(np.transpose(point_indices), point_values, leafsize=10, stat=1)
            interp_data = invdisttree(
                np.vstack((xi.flatten(), yi.flatten())).T,
                nnear=12,
                eps=.1,
                p=1,
                dub=64,
                weights=None
            )
            interp_data = np.reshape(interp_data, points_array.shape[0], points_array.shape[1])
        elif method == 'linear' or method == 'cubic':
            interp_data = scipy.interpolate.griddata(np.transpose(point_indices), point_values,
                                                     (xi, yi), method=method)
        
        plot(interp_data)
        # mask based on the cell coverage found using the MIN mode
        interp_data[coverage_data==coverage_no_data]=nodata
        plot(interp_data)
        # Write the data into the TIF on disk
        interp_band.WriteArray(interp_data)


    # release the temporary tif files and delete them
    #MRL: point_ -> points_
    points_band = None
    points_ds = None
    coverage_band=None
    coverage_ds = None
    if not DEBUGGING:
        os.remove(min_filename)
        os.remove(point_filename)
    print('Finished interpolating.') #JDV


if __name__ == "__main__":
    # input_file = r"C:\Data\Survey_Outline_Examples\Variable_resolution\H13070_VR\H13070_MB_VR_MLLW_1of2.bag"
    # input_file = r"C:\Data\BAG\GDAL_VR\H-10771\ExampleForEven\H-10771.bag"
    # input_file = r"C:\Data\BAG\GDAL_VR\F00788\F00788_MB_VR_MLLW_Final.bag"
    # input_file = r"C:\Data\BAG\BAG_to_XYZ_Problem\H13080_MB_1m_MLLW_1of1.bag"
    # input_file = r"C:\Data\Survey_Outline_Examples\Variable_resolution\H12993_MB_VR_MLLW.bag"
    
    #JDV: read arguments for input, output, and temp dir from the command line
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    #MRL: add default tmp_dir and method
    #MRL: add interpolation method option to CLI and make sure the method is 'idw', 'linear' or 'cubic', use 'linear as default
    if len(sys.argv) > 3:
        method = sys.argv[3]
    else:
        method = 'linear'

    if len(sys.argv) > 4:
        tmp_dir = sys.argv[4]
    else:
        tmp_dir = './'
        
    if method != 'idw' and method != 'linear' and method != 'cubic':
        print('invalid interpolation method {}, falling back to linear'.format(method))
        method = 'linear'
    
    if 1:
        interpolate_vr_bag(input_file, output_file, tmp_dir, use_blocks=True, method=method, nodata=3.4028234663852886e+38)
    if 0:
        dst_filename = input_file+".1m_mean.tif"
        print("\n\nmean")
        VRBag_to_TIF(input_file, dst_filename, sr_cell_size=1.0, mode=MEAN, use_blocks=True, nodata=3.4028234663852886e+38)
        dst_filename = input_file+".1m_min.tif"
        print("\n\nmin")
        VRBag_to_TIF(input_file, dst_filename, sr_cell_size=1.0, mode=MIN, use_blocks=True, nodata=3.4028234663852886e+38)
        dst_filename = input_file+".1m_max.tif"
        print("\n\nmax")
        VRBag_to_TIF(input_file, dst_filename, sr_cell_size=1.0, mode=MAX, use_blocks=True, nodata=3.4028234663852886e+38)
        dst_filename = input_file+".1m_mask.tif"
        print("\n\nmask")
        VRBag_to_TIF(input_file, dst_filename, sr_cell_size=1.0, mode=MASK, use_blocks=True, nodata=3.4028234663852886e+38)
        dst_filename = input_file+".1m_point.tif"
        print("\n\npoint")
        VRBag_to_TIF(input_file, dst_filename, sr_cell_size=1.0, mode=POINT, use_blocks=True, nodata=3.4028234663852886e+38)
        dst_filename = input_file + ".1m_data_mask.tif"
        print("\n\nmask")
        VRBag_to_TIF(input_file, dst_filename, sr_cell_size=1.0, mode=DATA_LOC_MASK, use_blocks=True, nodata=3.4028234663852886e+38)

