### gdalfile.py - DataLists IMproved
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
##
## gdalfile.py is part of CUDEM
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
### Commentary:
##
## GDAL Raster Parsers (GeoTiff, BAG, etc.)
##
### Code:

import os
import re
import numpy as np
from osgeo import gdal

from cudem import utils
from cudem import regions
from cudem import gdalfun
from cudem.grits import grits
from cudem.datalists.dlim import ElevationDataset

## ==============================================
## dlim GDALFile Module
## ==============================================
class GDALFile(ElevationDataset):
    """
    Providing a GDAL raster dataset parser.
    Process/Parse GDAL supported raster files.
    """
    
    def __init__(self,
                 weight_mask=None,
                 uncertainty_mask=None,
                 x_band=None,
                 y_band=None,
                 uncertainty_mask_to_meter=1,
                 open_options=None,
                 sample=None,
                 check_path=True,
                 super_grid=False,
                 band_no=1,
                 remove_flat=False,
                 node=None,
                 resample_and_warp=True,
                 yield_chunk=True,
                 chunk_step=None,
                 **kwargs):
        
        super().__init__(**kwargs)
        
        self.weight_mask = weight_mask 
        self.uncertainty_mask = uncertainty_mask 
        self.uncertainty_mask_to_meter = uncertainty_mask_to_meter
        self.open_options = open_options
        self.sample = sample
        self.check_path = check_path
        self.super_grid = super_grid
        self.band_no = band_no
        self.remove_flat = remove_flat 
        self.flat_removed = False
        self.x_band = x_band
        self.y_band = y_band
        self.node = node 
        self.resample_and_warp = resample_and_warp
        self.yield_chunk = yield_chunk
        self.chunk_step = utils.int_or(chunk_step)
        
        self.tmp_elev_band = None
        self.tmp_unc_band = None
        self.tmp_weight_band = None
        self.src_ds = None

        ## Network paths check skip
        if self.fn.startswith(('http', '/vsicurl/', 'BAG')) or self.fn.lower().endswith('.vrt'):
            self.check_path = False

            if self.fn.lower().endswith('.vrt'):
                self.chunk_size = 4096
                
            ## Set GDAL Config options for network performance
            ## Enable caching for seek operations (crucial for TIFFs)
            gdal.SetConfigOption('VSI_CACHE', 'TRUE')
            gdal.SetConfigOption('VSI_CACHE_SIZE', '50000000') # 50MB cache
            
            ## Prevent GDAL from scanning the parent directory of the URL
            gdal.SetConfigOption('GDAL_DISABLE_READDIR_ON_OPEN', 'EMPTY_DIR')
            
            ## Increase chunk size to reduce HTTP requests (128KB is default, 16MB is better for broadband)
            gdal.SetConfigOption('CPL_VSIL_CURL_CHUNK_SIZE', '16384000')
            gdal.SetConfigOption('GDAL_CACHEMAX', '1024')

        if self.valid_p() and self.src_srs is None:
            if self.infos.src_srs is None:
                self.src_srs = self.init_srs(self.fn)
            else:
                self.src_srs = self.infos.src_srs

                
    def destroy_ds(self):
        self.src_ds = None

        
    def init_srs(self, src_ds):
        """Initialize the srs from the gdal file."""
        
        if self.src_srs is None:
            return gdalfun.gdal_get_srs(src_ds)
        else:
            return self.src_srs

        
    def generate_inf(self, make_grid=True, make_block_mean=False, block_inc=None):
        """Generate metadata (INF) for the data source.
        
        Overrides ElevationDataset.generate_inf to leverage GDAL for basic info,
        but falls back to the parent class implementation to perform the scan
        for Mini-Grids and Block-Means if requested.
        """
        
        ## Quick Metadata from GDAL Header (Fast)
        with gdalfun.gdal_datasource(self.fn) as src_ds:
            if src_ds is not None:
                if self.src_srs is None:
                    self.infos.src_srs = self.init_srs(src_ds)
                else:
                    self.infos.src_srs = self.src_srs

                #is_remote = self.fn.startswith('/vsicurl/')
                is_remote = self.fn.startswith(('/vsicurl/', 'http')) or self.fn.lower().endswith('.vrt')
                    
                ds_infos = gdalfun.gdal_infos(src_ds)
                
                ## Get bounds from header
                this_region = regions.Region(src_srs=self.src_srs).from_geo_transform(
                    geo_transform=ds_infos['geoT'],
                    x_count=ds_infos['nx'],
                    y_count=ds_infos['ny']
                )
                
                ## Try to get Min/Max from statistics if available (Fast)
                force_scan = 0 if is_remote else 1
                
                band = src_ds.GetRasterBand(utils.int_or(self.band_no, 1))
                stats = band.GetStatistics(1, force_scan)
                if stats and (stats[0] != stats[1]):
                    this_region.zmin, this_region.zmax = stats[0], stats[1]
                else:
                    if self.verbose:
                         utils.echo_warning_msg(f"Skipping stats scan for remote file: {self.fn}")

                self.infos.minmax = this_region.export_as_list(include_z=True)
                self.infos.wkt = this_region.export_as_wkt()
                self.infos.numpts = ds_infos['nb']

        ## If Mini-Grid or Block-Mean is requested, we MUST scan the data.
        if make_grid or make_block_mean:
            ## We call the parent method. It will re-calculate numpts/minmax via scan,
            ## which is safer anyway, and generate the grids.
            return super().generate_inf(
                make_grid=make_grid, 
                make_block_mean=make_block_mean, 
                block_inc=block_inc
            )

        return self.infos

    
    def get_region(self):
        with gdalfun.gdal_datasource(self.fn) as src_ds:
            if src_ds is not None:
                ds_infos = gdalfun.gdal_infos(src_ds)
                return regions.Region(src_srs=self.src_srs).from_geo_transform(
                    geo_transform=ds_infos['geoT'],
                    x_count=ds_infos['nx'],
                    y_count=ds_infos['ny']
                )
        return None

    
    def get_srcwin(self, gt, x_size, y_size, node='grid'):
        """Calculate source window based on region and transform."""
        
        if self.region is not None:
            if self.invert_region:
                srcwin = (0, 0, x_size, y_size, node)
            else:
                if self.transform['transformer'] is not None:
                    if self.transform['trans_region'] is not None and self.transform['trans_region'].valid_p(check_xy=True):
                        srcwin = self.transform['trans_region'].srcwin(gt, x_size, y_size, node)
                    else:
                        srcwin = (0, 0, x_size, y_size, node)
                else:
                    srcwin = self.region.srcwin(gt, x_size, y_size, node)
        else:
            srcwin = (0, 0, x_size, y_size, node)
            
        return srcwin

    
    def yield_points(self):
        """Yield points from the raster.
        Handles resampling, chunking, and region subsetting.
        """

        if self.fn is None or (self.check_path and not os.path.exists(self.fn)):
            utils.echo_warning_msg(f'{self.fn} doesn\'t exist!')
            return None

        ## Attempt to open the file to verify it is a valid/supported raster.
        ## This prevents crashes in gdalfun.gdal_infos() or later processing.
        try:
            test_ds = gdal.Open(self.fn, gdal.GA_ReadOnly)
            if test_ds is None:
                utils.echo_warning_msg(f"Skipping {self.fn}: Not a valid GDAL raster source.")
                return None
            test_ds = None # Close dataset
        except Exception as e:
            utils.echo_warning_msg(f"Skipping {self.fn}: GDAL Open failed ({e})")
            return None
        finally:
            test_ds = None
        
        ## Process Open Options
        try:
            if isinstance(self.open_options, str):
                self.open_options = self.open_options.split('/')
        except Exception:
            self.open_options = None
            
        ## Initialize Infos
        self.dem_infos = gdalfun.gdal_infos(self.fn)
        if self.node is None:
            self.node = gdalfun.gdal_get_node(self.fn, 'pixel')
            
        ## Determine Resample/Warp Necessity
        if (self.x_inc is None and self.y_inc is None) or self.region is None:
            self.resample_and_warp = False
        else:
            ## Only upsample if target res is finer than source
            if self.transform['trans_region'] is not None:
                src_res = np.prod(self.transform['trans_region'].geo_transform(x_inc=self.dem_infos['geoT'][1])[:2])
                dst_res = np.prod(self.transform['trans_region'].geo_transform(x_inc=self.x_inc)[:2])
                if src_res > dst_res:
                    self.resample_and_warp = False
                
        if self.node == 'grid':
            self.resample_and_warp = False

        ndv = utils.float_or(gdalfun.gdal_get_ndv(self.fn), -9999)
        
        ## Setup Region
        if self.region is not None:
            self.warp_region = self.region.copy()
        else:
            if self.transform['transformer'] is not None:
                self.warp_region = self.transform['trans_region'].copy()
            else:
                self.warp_region = self.inf_region
            
        ## Temp Files
        tmp_elev_fn = utils.make_temp_fn(self.fn, temp_dir=self.cache_dir)
        tmp_unc_fn = utils.make_temp_fn(self.fn, temp_dir=self.cache_dir)
        tmp_weight_fn = utils.make_temp_fn(self.fn, temp_dir=self.cache_dir)

        # Handle Flats Removal (Grits)
        if self.remove_flat:
            grits_filter = grits.GritsFactory(
                mod='flats', src_dem=self.fn, cache_dir=self.cache_dir, verbose=False
            )._acquire_module()
            if grits_filter is not None:
                grits_filter()                
                self.fn = grits_filter.dst_dem
                self.flat_removed = True
            
        ## --- Resample / Warp ---
        if self.resample_and_warp:
            if self.transform['transformer'] is not None:
                self.transform['transformer'] = None

            ## Determine Resample Algorithm
            if self.sample_alg == 'auto':
                if self.stack_mode in ['min', 'max']:
                    self.sample_alg = self.stack_mode
                else:
                    self.sample_alg = 'bilinear' # Default fallback

            ## Open Source DataSet
            if self.open_options:
                self.src_ds = gdal.OpenEx(self.fn, open_options=self.open_options)
            else:
                self.src_ds = gdal.Open(self.fn)

            if self.src_ds is None: return None
            
            ## Prepare for Warp
            tmp_warp = utils.make_temp_fn(f'{self.fn}', temp_dir=self.cache_dir)
            in_bands = self.src_ds.RasterCount
            src_ds_config = gdalfun.gdal_infos(self.src_ds)
            
            ## If multi-band, extract bands first (gdalwarp limitation handling)
            ## This block handles extraction of Z, U, W bands to temp files
            tmp_ds_to_warp = self.fn
            if in_bands > 1:
                ## the srcwin for to extract data
                if self.transform['trans_region'] is not None:
                    srcwin_region = self.transform['trans_region'].copy()
                elif self.region is not None:
                    srcwin_region = self.region.copy()
                else:
                    srcwin_region = None

                if srcwin_region is not None:
                    srcwin = srcwin_region.srcwin(
                        src_ds_config['geoT'],
                        src_ds_config['nx'],
                        src_ds_config['ny'],
                        node='pixel'
                    )
                else:
                    srcwin = None

                self.tmp_elev_band = tmp_elev_fn
                if self.verbose:
                    utils.echo_msg(
                        (f'extracting elevation data from {self.fn} '
                        f'to {self.tmp_elev_band}')
                    )
                    
                self.tmp_elev_band, status = gdalfun.gdal_extract_band(
                    self.src_ds, self.tmp_elev_band,
                    band=self.band_no,
                    exclude=[],
                    srcwin=srcwin,
                    inverse=False
                )
                tmp_ds = self.tmp_elev_band
                if tmp_ds is None:
                    return(None)
                
                ## band is now 1!
                self.band_no = 1
                if utils.int_or(self.uncertainty_mask) is not None:
                    self.tmp_unc_band = tmp_unc_fn
                    if self.verbose:
                        utils.echo_msg(
                            (f'extracting uncertainty mask from {self.fn} '
                             f'to {self.tmp_unc_band}')
                        )
                        
                    self.tmp_unc_band, status = gdalfun.gdal_extract_band(
                        self.src_ds, self.tmp_unc_band,
                        band=self.uncertainty_mask,
                        exclude=[],
                        srcwin=srcwin,
                        inverse=False
                    )
                    self.uncertainty_mask = self.tmp_unc_band

                if utils.int_or(self.weight_mask) is not None:                    
                    self.tmp_weight_band = tmp_weight_fn
                    if self.verbose:
                        utils.echo_msg(
                            (f'extracting weight mask from {self.fn} '
                             f'to {self.tmp_weight_band}')
                        )
                        
                    self.tmp_weight_band, status = gdalfun.gdal_extract_band(
                        self.src_ds, self.tmp_weight_band,
                        band=self.weight_mask,
                        exclude=[],
                        srcwin=srcwin,
                        inverse=False
                    )
                    self.weight_mask = self.tmp_weight_band


            ## Perform Warp
            warped_ds = gdalfun.sample_warp(
                tmp_ds_to_warp, tmp_warp, self.x_inc, self.y_inc,
                src_region=self.warp_region,
                src_srs=self.transform['src_horz_crs'].to_proj4()
                if self.transform['src_horz_crs'] is not None \
                else None,
                dst_srs=self.transform['dst_horz_crs'].to_proj4() \
                if self.transform['dst_horz_crs'] is not None \
                else None,                
                sample_alg=self.sample_alg,
                ndv=ndv,
                verbose=False,
                co=["COMPRESS=DEFLATE", "TILED=YES"]
            )[0]
            
            ## Load Warped DS into Memory
            self.src_ds = gdal.Open(warped_ds) if warped_ds else None
            utils.remove_glob(tmp_warp)

        else:
            ## No Warp
            if self.open_options:
                self.src_ds = gdal.OpenEx(self.fn, open_options=self.open_options)
            else:
                self.src_ds = gdal.Open(self.fn)

        if self.src_ds is None:
            if self.verbose: utils.echo_error_msg(f'could not load raster file {self.fn}')
            return None

        self.src_dem_infos = gdalfun.gdal_infos(self.src_ds)
        
        ## --- Read Data ---
        band = self.src_ds.GetRasterBand(utils.int_or(self.band_no, 1))
        gt = self.src_ds.GetGeoTransform()
        ndv = utils.float_or(band.GetNoDataValue())
        
        ## Load Masks (Weights/Uncertainty) if present
        mask_band = None
        weight_band = None
        uncertainty_band = None
        nodata = ['{:g}'.format(-9999), 'nan', float('nan')]            
        if ndv is not None:
            nodata.append('{:g}'.format(ndv))

        band_data = count_data = weight_data = mask_data = None
                
        src_dem_x_inc = self.src_dem_infos['geoT'][1]
        src_dem_y_inc = -1*self.src_dem_infos['geoT'][5]
        src_dem_region = regions.Region().from_geo_transform(
            self.src_dem_infos['geoT'],
            self.src_dem_infos['nx'],
            self.src_dem_infos['ny']
        )

        ## ==============================================
        ## todo: always warp these to src_ds
        ## weight mask, each cell should have the corresponding weight
        ## weight_mask can either be a seperate gdal file or a band number
        ## corresponding to the appropriate band in src_ds
        ## ==============================================
        if self.weight_mask is not None:
            if utils.int_or(self.weight_mask) is not None:
                weight_band = self.src_ds.GetRasterBand(int(self.weight_mask))
            elif os.path.exists(self.weight_mask):
                # some numbers now return true here (file-descriptors),
                # check for int first!
                src_weight = gdalfun.sample_warp(
                    self.weight_mask, None, src_dem_x_inc, src_dem_y_inc,
                    src_srs=self.transform['src_horz_crs'].to_proj4() \
                    if self.transform['src_horz_crs'] is not None \
                    else None,
                    dst_srs=self.transform['dst_horz_crs'].to_proj4() \
                    if self.transform['dst_horz_crs'] is not None \
                    else None,
                    src_region=src_dem_region,
                    sample_alg=self.sample_alg,
                    ndv=ndv,
                    verbose=self.verbose,
                    co=["COMPRESS=DEFLATE", "TILED=YES"]
                )[0]
                weight_band = src_weight.GetRasterBand(1)

            else:
                utils.echo_warning_msg(
                    f'could not load weight mask {self.weight_mask}'
                )
                weight_band = None

        ## ==============================================
        ## uncertainty mask, each cell should have the corresponding uncertainty
        ## uncertainty_mask can either be a seperate gdal file or a band number
        ## corresponding to the appropriate band in src_ds
        ## ==============================================
        if self.uncertainty_mask is not None:
            if utils.int_or(self.uncertainty_mask):
                uncertainty_band = self.src_ds.GetRasterBand(
                    int(self.uncertainty_mask)
                )
            elif os.path.exists(self.uncertainty_mask):
                src_uncertainty = gdalfun.sample_warp(
                    self.uncertainty_mask, None, src_dem_x_inc, src_dem_y_inc,
                    src_srs=self.transform['src_horz_crs'].to_proj4() \
                    if self.transform['src_horz_crs'] is not None \
                    else None,
                    dst_srs=self.transform['dst_horz_crs'].to_proj4() \
                    if self.transform['dst_horz_crs'] is not None \
                    else None,
                    src_region=src_dem_region,
                    sample_alg='bilinear',
                    ndv=ndv,
                    verbose=self.verbose,
                    co=["COMPRESS=DEFLATE", "TILED=YES"]
                )[0]
                uncertainty_band = src_uncertainty.GetRasterBand(1)
            else:
                utils.echo_warning_msg(
                    f'could not load uncertainty mask {self.uncertainty_mask}'
                )
                uncertainty_band = None

        ## ==============================================        
        ## uncertainty from the vertical transformation
        ## ==============================================
        if self.transform['trans_fn_unc'] is not None:
            trans_uncertainty = gdalfun.sample_warp(
                self.transform['trans_fn_unc'], None, src_dem_x_inc, src_dem_y_inc,
                src_srs='+proj=longlat +datum=WGS84 +ellps=WGS84',
                dst_srs=self.transform['dst_horz_crs'].to_proj4() \
                if self.transform['dst_horz_crs'] is not None \
                else None,
                src_region=src_dem_region,
                sample_alg='bilinear',
                ndv=ndv,
                verbose=self.verbose,
                co=["COMPRESS=DEFLATE", "TILED=YES"]
            )[0]

            if uncertainty_band is not None:
                trans_uncertainty_band = trans_uncertainty.GetRasterBand(1)
                trans_uncertainty_arr = trans_uncertainty_band.ReadAsArray()
                uncertainty_arr = uncertainty_band.ReadAsArray()
                uncertainty_arr *= self.uncertainty_mask_to_meter
                uncertainty_arr = np.sqrt(uncertainty_arr**2 + trans_uncertainty_arr**2)
                uncertainty_band.WriteArray(uncertainty_arr)
                trans_uncertainty_band = None
            else:
                uncertainty_band = trans_uncertainty.GetRasterBand(1)

        ## ==============================================
        ## --- Chunked Yield ---
        ## ==============================================
        if self.yield_chunk:
            
            ## Determine Optimal Chunk Size
            if self.chunk_step is not None:
                ## User specified size
                read_chunk_size = self.chunk_step
            else:
                ## Auto-detect best size based on internal tiling
                try:
                    block_size = self.src_ds.GetRasterBand(1).GetBlockSize()
                    bx, by = block_size[0], block_size[1]
                    
                    if bx > 1 and by > 1:
                        target_size = 4096
                        read_chunk_size = (target_size // bx) * bx
                        if read_chunk_size == 0: read_chunk_size = bx
                    else:
                        read_chunk_size = 4096
                except:
                    ## Fallback
                    read_chunk_size = 4096

            for srcwin in utils.yield_srcwin(
                    n_size=(self.src_ds.RasterYSize, self.src_ds.RasterXSize),
                    n_chunk=read_chunk_size,
                    verbose=False
            ):                
                band_data = band.ReadAsArray(*srcwin).astype(float)
                if ndv is not None:
                    band_data[band_data == ndv] = np.nan

                if np.all(np.isnan(band_data)): continue

                ## Weights
                weight_data = np.ones(band_data.shape)
                if weight_band:
                    w_d = weight_band.ReadAsArray(*srcwin)
                    weight_data = w_d 

                ## Uncertainty
                uncertainty_data = np.zeros(band_data.shape)
                if uncertainty_band:
                    u_d = uncertainty_band.ReadAsArray(*srcwin)
                    uncertainty_data = u_d * self.uncertainty_mask_to_meter

                ## Coordinates
                if self.x_band is None and self.y_band is None:
                    ## Generate Coords from GT
                    x0, y0 = utils._pixel2geo(srcwin[0], srcwin[1], gt, node='pixel')
                    x1, y1 = utils._pixel2geo(srcwin[0]+srcwin[2], srcwin[1]+srcwin[3], gt, node='grid')
                    
                    lon_array = np.linspace(x0, x1, srcwin[2], endpoint=False)
                    lat_array = np.linspace(y0, y1, srcwin[3], endpoint=False)
                    
                    lon_data, lat_data = np.meshgrid(lon_array, lat_array)
                    
                    dataset = np.column_stack((
                        lon_data.flatten(), lat_data.flatten(), band_data.flatten(),
                        weight_data.flatten(), uncertainty_data.flatten()
                    ))
                else:
                    ## Use Band Coords
                    lon_band = self.src_ds.GetRasterBand(self.x_band)
                    lon_array = lon_band.ReadAsArray(*srcwin).astype(float)
                    lon_array[np.isnan(band_data)] = np.nan

                    lat_band = self.src_ds.GetRasterBand(self.y_band)
                    lat_array = lat_band.ReadAsArray(*srcwin).astype(float)
                    lat_array[np.isnan(band_data)] = np.nan
                    dataset = np.column_stack(
                        (lon_array[0], lat_array[0], band_data[0],
                         weight_data[0], uncertainty_data[0])
                    )

                ## Filter NaNs
                points = np.rec.fromrecords(dataset, names='x, y, z, w, u')
                points = points[~np.isnan(points['z'])]
                dataset = band_data = weight_data = None
                uncertainty_data = lat_array = lon_array = None
                utils.remove_glob(tmp_elev_fn, tmp_unc_fn, tmp_weight_fn)
                
                yield points
                
        else:
            ## Scanline Yield (Legacy)
            srcwin = self.get_srcwin(
                gt, self.src_ds.RasterXSize, self.src_ds.RasterYSize,
                node=self.node
            )
            with utils.ccp(total=(srcwin[1] + srcwin[3]), desc='Parsing Scanlines') as pbar:
                for y in range(srcwin[1], (srcwin[1] + srcwin[3]), 1):
                    pbar.update()
                    band_data = band.ReadAsArray(
                        srcwin[0], y, srcwin[2], 1
                    ).astype(float)
                    if ndv is not None and not np.isnan(ndv):
                        band_data[band_data == ndv] = np.nan

                    if np.all(np.isnan(band_data)):
                        continue

                    ## weights
                    if weight_band is not None:
                        weight_data = weight_band.ReadAsArray(
                            srcwin[0], y, srcwin[2], 1
                        )
                        weight_ndv = float(weight_band.GetNoDataValue())
                        if not np.isnan(weight_ndv):
                            weight_data[weight_data==weight_ndv] = np.nan
                    else:
                        weight_data = np.ones(band_data.shape)

                    ## uncertainty
                    if uncertainty_band is not None:
                        uncertainty_data = uncertainty_band.ReadAsArray(
                            srcwin[0], y, srcwin[2], 1
                        )
                        uncertainty_ndv = float(uncertainty_band.GetNoDataValue())
                        if not np.isnan(uncertainty_ndv):
                            uncertainty_data[uncertainty_data==uncertainty_ndv] = np.nan

                        if self.transform['trans_fn_unc'] is None:
                            uncertainty_data *= self.uncertainty_mask_to_meter

                    else:
                        uncertainty_data = np.zeros(band_data.shape)

                    ## convert grid array to points
                    if self.x_band is None and self.y_band is None:
                        x_precision = 12#len(str(gt[0]).split('.')[-1])
                        y_precision = 12#len(str(gt[3]).split('.')[-1])
                        while True:
                            ## 'self.node to 'pixel', this breaks if set to
                            ## 'grid' even if 'grid-node'
                            geo_x_origin, geo_y_origin = utils._pixel2geo(
                                srcwin[0], y, gt,
                                node='pixel',
                                x_precision=x_precision,
                                y_precision=y_precision
                            ) 
                            geo_x_end, geo_y_end = utils._pixel2geo(
                                srcwin[0] + srcwin[2], y, gt,
                                node='grid',
                                x_precision=x_precision,
                                y_precision=y_precision
                            )
                            lon_array = np.arange(geo_x_origin, geo_x_end, gt[1])

                            if lon_array.shape == band_data[0].shape:
                                break
                            else:
                                if x_precision < 0 or y_precision < 0:
                                    break

                                x_precision -= 1
                                y_precision -= 1

                        lat_array = np.zeros((lon_array.shape))
                        lat_array[:] = geo_y_origin

                        try:
                            assert lon_array.shape == lat_array.shape
                            assert lon_array.shape == band_data[0].shape
                            dataset = np.column_stack(
                                (lon_array, lat_array, band_data[0],
                                 weight_data[0], uncertainty_data[0])
                            )
                        except Exception as e:
                            utils.echo_error_msg(e)
                            pass
                    else:
                        lon_band = self.src_ds.GetRasterBand(self.x_band)
                        lon_array = lon_band.ReadAsArray(
                            srcwin[0], y, srcwin[2], 1
                        ).astype(float)
                        lon_array[np.isnan(band_data)] = np.nan

                        lat_band = self.src_ds.GetRasterBand(self.y_band)
                        lat_array = lat_band.ReadAsArray(
                            srcwin[0], y, srcwin[2], 1
                        ).astype(float)
                        lat_array[np.isnan(band_data)] = np.nan
                        dataset = np.column_stack(
                            (lon_array[0], lat_array[0], band_data[0],
                             weight_data[0], uncertainty_data[0])
                        )

                    points = np.rec.fromrecords(dataset, names='x, y, z, w, u')
                    points =  points[~np.isnan(points['z'])]
                    dataset = band_data = weight_data = None
                    uncertainty_data = lat_array = lon_array = None
                    utils.remove_glob(tmp_elev_fn, tmp_unc_fn, tmp_weight_fn)
                    yield points

        self.src_ds = None
        if self.remove_flat and self.flat_removed:
            utils.remove_glob(self.fn)

            
## ==============================================
## BAGFile Module
## ==============================================
class BAGFile(ElevationDataset):
    """Providing a BAG raster dataset parser.
    Process supergrids at native resolution if they exist.
    """

    def __init__(self,
                 explode = False,
                 force_vr = False,
                 vr_resampled_grid = False,
                 vr_interpolate = True,
                 vr_strategy = 'MIN',
                 min_weight = 0,
                 **kwargs):
        super().__init__(**kwargs)
        self.explode = explode
        self.min_weight = utils.float_or(min_weight, 0)
        self.force_vr = force_vr
        self.vr_resampled_grid = vr_resampled_grid
        self.vr_interpolate = vr_interpolate
        if self.vr_interpolate:
            self.vr_resampled_grid = False
            
        self.vr_strategy = vr_strategy
        if self.src_srs is None:
            self.src_srs = self.init_srs()

            
    def init_srs(self):
        if self.src_srs is None:
            src_horz, src_vert = gdalfun.split_srs(gdalfun.gdal_get_srs(self.fn), as_epsg=False)
            if src_horz is None and src_vert is None: return None
            
            if src_vert is None: src_vert = 5866 # Default to MLLW if unknown.
            self.src_srs = gdalfun.combine_epsgs(src_horz, src_vert, name='BAG Combined')
            return self.src_srs
        else:
            return self.src_srs

        
    def generate_inf(self, make_grid=True, make_block_mean=False, block_inc=None):
        """Generate metadata (INF) for the BAG.
        """
        
        ## BAG specific basic info
        if self.src_srs is None: self.infos.src_srs = self.init_srs()
        else: self.infos.src_srs = self.src_srs

        with gdalfun.gdal_datasource(self.fn) as src_ds:
            if src_ds is not None:
                ds_infos = gdalfun.gdal_infos(src_ds)
                this_region = regions.Region(src_srs=self.src_srs).from_geo_transform(
                    geo_transform=ds_infos['geoT'],
                    x_count=ds_infos['nx'],
                    y_count=ds_infos['ny']
                )
                
                ## BAG Band 1 is Elevation
                zr = src_ds.GetRasterBand(1).ComputeRasterMinMax() 
                this_region.zmin, this_region.zmax = zr[0], zr[1] 
                
                self.infos.minmax = this_region.export_as_list(include_z=True)
                self.infos.wkt = this_region.export_as_wkt()
                self.infos.numpts = ds_infos['nb']

        ## Delegate to parent for Grids if requested
        if make_grid or make_block_mean:
            return super().generate_inf(make_grid, make_block_mean, block_inc)
            
        return self.infos

    
    def parse(self):
        from .dlim import DatasetFactory
        
        ## Check intersection first
        if self.region is not None:
            check_region = self.transform['trans_region'] if self.transform['trans_region'] else self.region
            if not regions.regions_intersect_p(self.inf_region, check_region):
                return

        mt = gdal.Info(self.fn, format='json')['metadata']['']
        ds_infos = gdalfun.gdal_infos(self.fn)
        x_res = ds_infos['geoT'][1]
        
        ## Heuristic for weight based on resolution
        sub_weight = max((3 * (10 if x_res <=3 else 1))/x_res, self.min_weight)
        
        oo = []
        if self.data_region is not None and self.data_region.valid_p():
            oo.append(f'MINX={self.data_region.xmin}')
            oo.append(f'MAXX={self.data_region.xmax}')
            oo.append(f'MINY={self.data_region.ymin}')
            oo.append(f'MAXY={self.data_region.ymax}')

        ## VR Processing
        is_vr = ('HAS_SUPERGRIDS' in mt and mt['HAS_SUPERGRIDS'] == 'TRUE') or self.force_vr
        
        if is_vr:
            if self.explode:
                # Process Supergrids separately
                utils.echo_debug_msg(f'Initializing BAG EXPLODE {self.fn}')
                oo.append("MODE=LIST_SUPERGRIDS")
                src_ds = gdal.OpenEx(self.fn, open_options=oo)
                sub_datasets = src_ds.GetSubDatasets()
                src_ds = None

                with utils.ccp(total=len(sub_datasets), desc=f'Parsing supergrids: {self.fn}') as pbar:                
                    for sub_dataset in sub_datasets:
                        pbar.update()
                        # Parse resolution from desc to set weight
                        res_match = re.findall(r'\d+\.\d+|\d+', sub_dataset[-1].split(',')[-2])
                        res_val = float(res_match[0]) if res_match else x_res
                        s_w = max((3 * (10 if res_val <=3 else 1))/res_val, self.min_weight)

                        sub_ds = DatasetFactory(
                            **self._set_params(
                                mod=sub_dataset[0],
                                data_format=200,
                                src_srs=self.src_srs,
                                band_no=1,
                                uncertainty_mask_to_meter=0.01,
                                check_path=False,
                                super_grid=True,
                                #node='pixel',
                                #weight_multiplier=sub_weight,
                                weight=sub_weight*(self.weight if self.weight else 1),
                                uncertainty_mask=2,
                            )
                        )._acquire_module()
                                                
                        self.data_entries.append(sub_ds)
                        sub_ds.initialize()
                        for gdal_ds in sub_ds.parse():
                            yield gdal_ds
                            utils.remove_glob(f'{gdal_ds.fn}.inf')

            elif self.vr_resampled_grid or self.vr_interpolate:
                ## Use GDAL internal VR modes
                utils.echo_debug_msg(f'Initializing BAG VR {self.fn}')
                mode = "RESAMPLED_GRID" if self.vr_resampled_grid else "INTERPOLATED"
                oo.append(f"MODE={mode}")
                oo.append(f"RES_STRATEGY={self.vr_strategy}")

                sub_ds = DatasetFactory(
                    **self._set_params(
                        mod=self.fn,
                        data_format=200,
                        band_no=1,
                        open_options=oo,
                        weight=sub_weight*(self.weight if self.weight else 1),
                        uncertainty_mask=2,
                        uncertainty_mask_to_meter=0.01,
                    )
                )._acquire_module()
                                
                self.data_entries.append(sub_ds)
                sub_ds.initialize()
                for gdal_ds in sub_ds.parse():
                    yield gdal_ds
            else:
                utils.echo_debug_msg(f'Initializing BAG VRBAG {self.fn}')
                tmp_bag_as_tif = utils.make_temp_fn(
                    '{}_tmp.tif'.format(utils.fn_basename2(self.fn))
                )
                ## scale cellsize to meters,
                ## todo: check if input is degress/meters/feet
                #self.x_inc * 111120 
                sr_cell_size = None
                vrbag.interpolate_vr_bag(
                    self.fn, tmp_bag_as_tif,
                    self.cache_dir, sr_cell_size=sr_cell_size,
                    use_blocks=True,
                    method='linear',
                    nodata=3.4028234663852886e+38
                )
                sub_ds = DatasetFactory(
                    **self._set_params(
                        mod=tmp_bag_as_tif,
                        data_format=200,
                        band_no=1,
                        weight=sub_weight*(self.weight if self.weight is not None else 1),
                        uncertainty_mask=2,
                        uncertainty_mask_to_meter=0.01,
                    )
                )._acquire_module()
                self.data_entries.append(sub_ds)
                sub_ds.initialize()
                for gdal_ds in sub_ds.parse():
                    yield(gdal_ds)

        else:
            ## Standard BAG
            utils.echo_debug_msg(f'Initializing BAG STANDARD {self.fn}')
            sub_ds = DatasetFactory(
                **self._set_params(
                    mod=self.fn,
                    data_format=200,
                    band_no=1,
                    uncertainty_mask=2,
                    uncertainty_mask_to_meter=0.01,
                    weight=sub_weight*(self.weight if self.weight is not None else 1),
                )
            )._acquire_module()

            self.data_entries.append(sub_ds)
            sub_ds.initialize()
            yield sub_ds

### End
