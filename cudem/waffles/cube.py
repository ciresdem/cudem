### flatten.py
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
##
## flatten.py is part of CUDEM
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
### Code:

from cudem import gdalfun
from cudem.waffles.waffles import Waffle

class WafflesCUBE(Waffle):
    """
    BathyCUBE - doesn't seem to work as expected, 
    likely doing something wrong here...

    https://github.com/noaa-ocs-hydrography/bathygrid
    
    """
    
    def __init__(
            self,
            chunk_size=None,
            chunk_buffer=40,
            **kwargs):
        """generate a `CUBE` dem"""
        
        super().__init__(**kwargs)
        self.chunk_size = utils.int_or(chunk_size)
        self.chunk_step = None
        self.chunk_buffer = chunk_buffer

        
    def run(self):

        try:
            import bathycube.cube as cube
        except:
            utils.echo_error_msg(
                'could not import bathycube, it may not be installed'
            )
            return(self)

        from scipy import optimize
        
        # _x = []
        # _y = []
        # _z = []
        # for this_xyz in self.yield_xyz():
        #     _x.append(this_xyz.x)
        #     _y.append(this_xyz.y)
        #     _z.append(this_xyz.z)

        # _x = np.array(_x)
        # _y = np.array(_y)
        # _z = np.array(_z)

        # print(_z.size)
        # print(_z.shape)
        # print(len(_z))
        
        if self.chunk_size is None:
            n_chunk = int(self.ds_config['nx'] * .1)
            n_chunk = 10 if n_chunk < 10 else n_chunk
        else: n_chunk = self.chunk_size
        if self.chunk_step is None:
            n_step = int(n_chunk/2)
        else: n_step = self.chunk_step

        points_ds = gdal.Open(self.stacked_rasters['z'])
        points_band = points_ds.GetRasterBand(1)
        points_no_data = points_band.GetNoDataValue()
        
        #try:
        interp_ds = points_ds.GetDriver().Create(
            self.fn,
            points_ds.RasterXSize,
            points_ds.RasterYSize,
            bands=1,
            eType=points_band.DataType,
            options=["BLOCKXSIZE=256",
                     "BLOCKYSIZE=256",
                     "TILED=YES",
                     "COMPRESS=LZW",
                     "BIGTIFF=YES"]
        )
        interp_ds.SetProjection(points_ds.GetProjection())
        interp_ds.SetGeoTransform(points_ds.GetGeoTransform())
        #interp_ds.SetGeoTransform(ds_config['geoT'])
        interp_band = interp_ds.GetRasterBand(1)
        interp_band.SetNoDataValue(np.nan)

        uncert_ds = points_ds.GetDriver().Create(
            '{}_unc.tif'.format(self.name),
            points_ds.RasterXSize,
            points_ds.RasterYSize,
            bands=1,
            eType=points_band.DataType,
            options=["BLOCKXSIZE=256",
                     "BLOCKYSIZE=256",
                     "TILED=YES",
                     "COMPRESS=LZW",
                     "BIGTIFF=YES"]
        )
        uncert_ds.SetProjection(points_ds.GetProjection())
        uncert_ds.SetGeoTransform(points_ds.GetGeoTransform())
        #uncert_ds.SetGeoTransform(ds_config['geoT'])
        uncert_band = uncert_ds.GetRasterBand(1)
        uncert_band.SetNoDataValue(np.nan)
            
        # except:
        #     return(self)
        
        if self.verbose:
            utils.echo_msg(
                'buffering srcwin by {} pixels'.format(self.chunk_buffer)
            )

        _x, _y = np.mgrid[0:ds_config['nx'], 0:ds_config['ny']]
        _x = _x.ravel()
        _y = _y.ravel()

        _z = points_band.ReadAsArray()
        _z = _z.T
        _z = _z.ravel()
        point_indices = np.nonzero(_z != ds_config['ndv'])
        point_values = _z[point_indices]        
        xi = _x[point_indices]
        yi = _y[point_indices]

        thu = np.ones(len(point_values))
        thu[:] = 2
        #thu = thu[point_indices]
        a = .25
        b = .0075
        tvu = np.sqrt((a**2 + (b * abs(point_values))**2))
        numrows, numcols = (ds_config['nx'], ds_config['ny'])
        res_x, res_y = ds_config['geoT'][1], ds_config['geoT'][5]*-1
        depth_grid, uncertainty_grid, ratio_grid, numhyp_grid = cube.run_cube_gridding(
            point_values, thu, tvu, xi, yi, self.ds_config['nx'], self.ds_config['ny'],
            min(_x), max(_y), 'local', 'order1a', 1, 1)
        print(depth_grid)
        depth_grid = np.flip(depth_grid)
        depth_grid = np.fliplr(depth_grid)
        interp_band.WriteArray(depth_grid)
        uncertainty_grid = np.flip(uncertainty_grid)
        uncertainty_grid = np.fliplr(uncertainty_grid)
        uncert_band.WriteArray(uncertainty_grid)
            
        return(self)


class WafflesBGrid(Waffle):
    """
    BathyCUBE - doesn't seem to work as expected, 
    likely doing something wrong here...

    https://github.com/noaa-ocs-hydrography/bathygrid
    
    """
    
    def __init__(
            self,
            chunk_size=None,
            chunk_buffer=40,
            **kwargs):
        """generate a `CUBE` dem"""
        
        super().__init__(**kwargs)
        self.chunk_size = utils.int_or(chunk_size)
        self.chunk_step = None
        self.chunk_buffer = chunk_buffer

        
    def run(self):

        try:
            import bathycube.cube as cube
        except:
            utils.echo_error_msg(
                'could not import bathycube, it may not be installed'
            )
            return(self)

        from scipy import optimize
        
        # self._stacks_array(
        #     out_name='{}_cube_stack'.format(self.name),
        #     supercede=self.supercede
        # )
        # n = '{}_cube_stack_s.tif'.format(self.name)
        # w = '{}_cube_stack_w.tif'.format(self.name)
        # c = '{}_cube_stack_c.tif'.format(self.name)
        _x = []
        _y = []
        _z = []
        for this_xyz in self.yield_xyz():
            _x.append(this_xyz.x)
            _y.append(this_xyz.y)
            _z.append(this_xyz.z)

        _x = np.array(_x)
        _y = np.array(_y)
        _z = np.array(_z)

        dtyp = [('x', np.float64),
                ('y', np.float64),
                ('z', np.float32),
                ('tvu', np.float32),
                ('thu', np.float32)]
        data = np.empty(len(_x), dtype=dtyp)
        data['x'] = _x
        data['y'] = _y
        data['z'] = _z
        
        # print(_z.size)
        # print(_z.shape)
        # print(len(_z))
        
        if self.chunk_size is None:
            n_chunk = int(ds_config['nx'] * .1)
            n_chunk = 10 if n_chunk < 10 else n_chunk
        else: n_chunk = self.chunk_size
        if self.chunk_step is None:
            n_step = int(n_chunk/2)
        else: n_step = self.chunk_step

        if self.verbose:
            utils.echo_msg(
                'buffering srcwin by {} pixels'.format(self.chunk_buffer)
            )

        thu = np.ones(len(_z))
        thu[:] = 2
        #thu = thu[point_indices]
        a = .25
        b = .0075
        tvu = np.sqrt((a**2 + (b * abs(_z))**2))

        data['tvu'] = tvu
        data['thu'] = thu
        numrows, numcols = (ds_config['nx'], ds_config['ny'])
        res_x, res_y = ds_config['geoT'][1], ds_config['geoT'][5]*-1


        from bathygrid.convenience import create_grid
        # a single resolution grid that is entirely within computer memory
        bg = create_grid(folder_path='./', grid_type='single_resolution')
        # add points from two multibeam lines, EPSG:26917 with vertical reference 'waterline'
        bg.add_points(data, 'test1', ['line1', 'line2'], 26965, 'waterline')
        print(bg.points_count)
        assert not bg.is_empty

        # grid by looking up the mean depth of each tile to determine resolution
        bg.grid()

        print(bg.resolutions)
        out_tif = os.path.join(bg.output_folder, self.fn)
        bg.export(out_tif, export_format='geotiff')
        
        # # check to see if the new bags are written
        # new_bag = os.path.join(bg.output_folder, 'outtiff_0.5.bag')
        # new_bag_two = os.path.join(bg.output_folder, 'outtiff_1.0.bag')
        # assert os.path.exists(new_bag)
        # assert os.path.exists(new_bag_two)

        # # Get the total number of cells in the variable resolution grid for each resolution
        # bg.cell_count
        # bg.coverage_area
        
        # depth_grid, uncertainty_grid, ratio_grid, numhyp_grid = cube.run_cube_gridding(
        #     point_values,
        #     thu,
        #     tvu,
        #     xi,
        #     yi,
        #     ds_config['nx'],
        #     ds_config['ny'],
        #     min(_x),
        #     max(_y),
        #     'local',
        #     'order1a',
        #     1,
        #     1,
        # )
        # print(depth_grid)
        # depth_grid = np.flip(depth_grid)
        # depth_grid = np.fliplr(depth_grid)
        # interp_band.WriteArray(depth_grid)
        # uncertainty_grid = np.flip(uncertainty_grid)
        # uncertainty_grid = np.fliplr(uncertainty_grid)
        # uncert_band.WriteArray(uncertainty_grid)
            
        return(self)

    
### End
