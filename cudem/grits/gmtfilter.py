### gmtfilter.py - GRId filTerS
##
## Copyright (c) 2024 - 2025 Regents of the University of Colorado
##
## gmtfilter.py is part of cudem
##
## Permission is hereby granted, free of charge, to any person obtaining a copy 
## of this software and associated documentation files (the "Software"), to deal 
## in the Software without restriction, including without limitation the rights 
## to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies 
## of the Software, and to permit persons to whom the Software is furnished to do so, 
## subject to the following conditions:
##
## The above copyright notice and this permission notice shall be included in all copies or
## substantial portions of the Software.
##
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
## INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR 
## PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE 
## FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
## ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
## SOFTWARE.
##
### Commentary:
##
## Wraps GMT's grdfilter module.
## Supports both subprocess (CLI) and PyGMT execution modes.
##
### Code:

from cudem import utils
from cudem.grits import grits

try:
    import pygmt
    HAS_PYGMT = True
except ImportError:
    HAS_PYGMT = False

class GMTgrdfilter(grits.Grits):
    """Filter a DEM through GMT's `grdfilter`.
    See `gmt grdfilter --help` for details on filter types and distance modes.

    Parameters:
    -----------
    filter_type : str
        The grdfilter filter definition (e.g., 'c100', 'g500').
        Maps to the `-F` argument.
    dist : str
        The distance calculation mode (e.g., '0' for Cartesian, '4' for Geodetic).
        Maps to the `-D` argument. Default is '0'.
    node : str
        Registration type: 'pixel' or 'grid'. Default is 'pixel'.
    use_pygmt : bool
        If True and PyGMT is installed, use the python bindings instead of subprocess.
    """
    
    def __init__(self, filter_type: str = 'c3s', dist: str = '0',
                 node: str = 'pixel', use_pygmt: bool = False, **kwargs: any):
        super().__init__(**kwargs)
        self.filter_type = filter_type
        self.dist = str(dist)
        self.node = node
        self.use_pygmt = use_pygmt

        
    def run(self):
        """Execute the GMT grdfilter."""
        
        ## Ensure destination is set (Grits base does this, but being safe)
        if self.dst_dem is None:
            self.dst_dem = utils.make_temp_fn(
                f'{utils.fn_basename2(self.src_dem)}_filtered.tif', 
                temp_dir=self.cache_dir
            )

        if self.use_pygmt and HAS_PYGMT:
            return self._run_pygmt()
        else:
            return self._run_subprocess()

        
    def _run_pygmt(self):
        """Run using PyGMT bindings."""
        
        if self.verbose:
            utils.echo_msg(f"Running PyGMT grdfilter: -F{self.filter_type} -D{self.dist}")

        try:
            ## Registration: 'p' for pixel (gridline), 'g' for grid
            reg = 'p' if self.node == 'pixel' else 'g'
            
            pygmt.grdfilter(
                grid=self.src_dem,
                filter=self.filter_type,
                distance=self.dist,
                registration=reg,
                outgrid=self.dst_dem,
                verbose='q' if not self.verbose else None
            )
            return self.dst_dem, 0
            
        except Exception as e:
            utils.echo_error_msg(f"PyGMT grdfilter failed: {e}")
            # Fallback to subprocess? Or just fail.
            # return self._run_subprocess() 
            return self.src_dem, -1

        
    def _run_subprocess(self):
        """Run using system GMT command."""
        
        ## Registration flag
        reg_flag = '-rp' if self.node == 'pixel' else ''
        
        ## Construct Command
        ## -G<out>=gd:GTiff forces GeoTIFF output via GDAL driver in GMT
        cmd = (
            f"gmt grdfilter -V {self.src_dem} "
            f"-G{self.dst_dem}=gd:GTiff "
            f"-F{self.filter_type} "
            f"-D{self.dist} {reg_flag}"
        )

        if self.verbose:
            utils.echo_msg(f"Running GMT: {cmd}")

        out, status = utils.run_cmd(cmd, verbose=self.verbose)

        if status != 0:
            utils.echo_error_msg(f"GMT grdfilter failed (Status {status})")
            return self.src_dem, -1
            
        return self.dst_dem, 0

    
### End
