### sdb.py
##
## Copyright (c) 2025 - 2026 Regents of the University of Colorado
##
## sdb.py is part of CUDEM
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
## Waffles SDB: Satellite Derived Bathymetry Module
##
## Integrates Sentinel-2 fetching (CDSE) with Physics-based RF Inversion.
##
### Code:

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Waffles SDB: Satellite Derived Bathymetry Module
Integrates Sentinel-2 fetching (CDSE) with Physics-based RF Inversion.
"""

import os
import sys
import numpy as np
from osgeo import gdal
from cudem import utils
from cudem import gdalfun
from cudem.waffles.waffles import Waffle
from cudem.fetches import fetches

## Attempt to import the Sentinel-2 fetcher (handles installed or local file)
try:
    from cudem.fetches.cdse import Sentinel2
    HAS_CDSE = True
except ImportError:
    try:
        ## Fallback for local cdse.py file
        import importlib.util
        spec = importlib.util.spec_from_file_location("cdse", "cdse.py")
        if spec and spec.loader:
            cdse = importlib.util.module_from_spec(spec)
            sys.modules["cdse"] = cdse
            spec.loader.exec_module(cdse)
            Sentinel2 = cdse.Sentinel2
            HAS_CDSE = True
        else:
            HAS_CDSE = False
    except Exception:
        HAS_CDSE = False

## Optional ML dependencies
try:
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.linear_model import LinearRegression
    HAS_ML = True
except ImportError:
    HAS_ML = False

class WafflesSDB(Waffle):
    """Satellite Derived Bathymetry (SDB) via Random Forest and Physics-based inversion.
    
    Fetches Sentinel-2 data on-the-fly if time parameters are provided.
    
    Args:
        max_depth (float): Maximum depth to resolve (cutoff for training/prediction).
        n_trees (int): Number of trees in Random Forest.
        smooth (bool): Apply median filter to features before training.
        physics (bool): Enforce physics-based constraints (Stumpf ratio).
        time (str): Time range for Sentinel-2 fetch (e.g. '2023-01-01/2023-03-01').
        cloud (float): Max cloud cover percentage for Sentinel-2 fetch.
    """

    def __init__(self, max_depth=30.0, n_trees=100, smooth=False, physics=True, 
                 time=None, cloud=10.0, **kwargs):
        super().__init__(**kwargs)
        self.max_depth = utils.float_or(max_depth)
        self.n_trees = utils.int_or(n_trees)
        self.smooth = smooth
        self.physics = physics
        
        ## Fetch parameters
        self.time = time
        self.cloud = utils.float_or(cloud)
        
        ## Internal state
        self.predictor_path = None
        
        if not HAS_ML:
            utils.echo_error_msg("SDB module requires scikit-learn.")

            
    def _init(self):
        """Initialize the SDB module.
        If time window is provided, fetch, download, mosaic, and stack Sentinel-2 data.
        """
        
        if self.time:
            if not HAS_CDSE:
                utils.echo_error_msg("cudem.fetches.cdse (or local cdse.py) not found. Cannot fetch Sentinel-2 data.")
                return -1
                
            self._fetch_and_process_s2()
        
        return 0

    
    def _fetch_and_process_s2(self):
        """Fetches Sentinel-2 bands, downloads them, and creates a VRT stack.
        """
        
        utils.echo_msg(f"Fetching Sentinel-2 data for {self.region.format('str')} [{self.time}]...")
        
        ## Parse Time
        try:
            t_start, t_end = self.time.split('/')
        except ValueError:
            utils.echo_error_msg("Time must be in format YYYY-MM-DD/YYYY-MM-DD")
            return

        ## Setup Cache Directory
        if self.cache_dir is None:
            self.cache_dir = os.getcwd()
        
        self.cache_dir = os.path.abspath(self.cache_dir)
        s2_cache = os.path.join(self.cache_dir, 's2_cache')
        
        try:
            if not os.path.exists(s2_cache):
                os.makedirs(s2_cache, exist_ok=True)
        except OSError as e:
            utils.echo_error_msg(f"Failed to create cache directory {s2_cache}: {e}")
            return

        # self.predictor_path = os.path.join(s2_cache, "s2_stack.vrt")
        # if os.path.exists(self.predictor_path):
        #     return
        
        ## Run Fetcher
        s2_fetcher = Sentinel2(
            time_start=t_start,
            time_end=t_end,
            max_cloud_cover=self.cloud,
            src_region=self.src_region,
            verbose=self.verbose,
            outdir=s2_cache
        )
        s2_fetcher.run()
        
        if not s2_fetcher.results:
            utils.echo_warning_msg("No Sentinel-2 scenes found for the given criteria.")
            return

        ## Download Data using fetches
        band_files = {'B02': [], 'B03': [], 'B04': []}
        utils.echo_msg(f"Downloading {len(s2_fetcher.results)} band files to {s2_cache}...")

        s2_daemon = fetches.fetch_results(s2_fetcher, check_size=True)
        s2_daemon.daemon = True
        s2_daemon.start()
        s2_daemon.join()
        
        s2_results = s2_daemon.results
        
        ## Get Headers (Authorization) from the fetcher instance
        #headers = s2_fetcher.headers 

        #utils.echo_msg(s2_results)
        for entry in s2_results:
            fname = os.path.join(s2_fetcher._outdir, entry[1])
            
            ## Identify band
            band_key = None
            if 'B02' in fname: band_key = 'B02'
            elif 'B03' in fname: band_key = 'B03'
            elif 'B04' in fname: band_key = 'B04'
            
            if band_key:
                local_path = os.path.join(s2_cache, fname)

                ## Below will fetch each file one by one, use fetches.fetch_results
                ## to fetch data in threads to get them quicker...
                ## Initialize Fetch object with headers (for auth) and verbosity
                #f = Fetch(url=url, headers=headers, verbose=self.verbose)
                
                ## fetch_file handles retries, resumption, and size checks automatically
                #status = f.fetch_file(local_path)
                
                #if status == 0:
                band_files[band_key].append(local_path)
                #else:
                #    utils.echo_warning_msg(f"Failed to download {fname}")

        ## Create VRT Stack of Sentinel data
        mosaic_vrts = []
        valid_bands = [b for b in ['B02', 'B03', 'B04'] if band_files[b]]
        
        if len(valid_bands) < 2:
            utils.echo_error_msg("Insufficient bands fetched (need at least B02 and B03).")
            return

        for band in valid_bands:
            files = band_files[band]
            if not files: continue
            
            mosaic_path = os.path.join(s2_cache, f"mosaic_{band}.vrt")

            if self.dst_srs:
                gdal.Warp(
                    mosaic_path, 
                    files, 
                    format="VRT", 
                    dstSRS=self.dst_srs, 
                    srcNodata=0, 
                    dstNodata=0
                )
            else:
                gdal.BuildVRT(
                    mosaic_path, 
                    files, 
                    options=gdal.BuildVRTOptions(srcNodata=0, VRTNodata=0)
                )
            
            mosaic_vrts.append(mosaic_path)
            
            #gdal.BuildVRT(mosaic_path, files, options=gdal.BuildVRTOptions(srcNodata=0, VRTNodata=0))
            #mosaic_vrts.append(mosaic_path)

        self.predictor_path = os.path.join(s2_cache, "s2_stack.vrt")
        gdal.BuildVRT(self.predictor_path, mosaic_vrts, options=gdal.BuildVRTOptions(separate=True))
        
        utils.echo_msg(f"Created Sentinel-2 Stack: {self.predictor_path}")

        
    def _extract_features(self, band_arrays):
        """Feature Engineering"""
        
        eps = 1e-6
        b_blue = np.maximum(band_arrays[0], eps)
        b_green = np.maximum(band_arrays[1], eps)
        b_red = band_arrays[2] if len(band_arrays) > 2 else np.zeros_like(b_blue)

        log_b = np.log(b_blue)
        log_g = np.log(b_green)
        stumpf = np.log(1000 * b_blue) / np.log(1000 * b_green)
        brightness = (b_blue + b_green + b_red) / 3.0
        
        features_stack = np.stack([
            b_blue, b_green, b_red,
            log_b, log_g,
            stumpf,
            brightness
        ], axis=-1)
        
        rows, cols, n_feats = features_stack.shape
        X = features_stack.reshape(rows * cols, n_feats)
        return X, (rows, cols)

    
    def run(self):
        ## Get Predictor
        predictor_fn = self.predictor_path
        if not predictor_fn:
            for ds in self.data.datasets:
                if ds.type == 'gdal':
                    predictor_fn = ds.fn
                    break
        
        if not predictor_fn:
            utils.echo_error_msg("No Sentinel-2 data found. Provide --time argument or input raster.")
            return

        ## Gather Training Points
        training_xyz = []
        for xyz in self.stack_ds.yield_xyz():
            if xyz.valid_p():
                if xyz.z < 0 and abs(xyz.z) <= self.max_depth:
                    training_xyz.append([xyz.x, xyz.y, xyz.z])

        if not training_xyz:
            utils.echo_error_msg("No valid training points found in input data.")
            return

        training_data = np.array(training_xyz)
        
        ## Load Predictors (Sentinel-2 data)
        ds = gdal.Open(predictor_fn)
        gt = ds.GetGeoTransform()
        srcwin = self.src_region.srcwin(gt, ds.RasterXSize, ds.RasterYSize)
        utils.echo_msg(f'Sentinel-2 coverage: {srcwin}')
        # if srcwin[2] <= 0 or srcwin[3] <= 0:
        #     utils.echo_warning_msg("Region outside Sentinel-2 coverage.")
        #     return 

        band_data = []
        for b in range(1, ds.RasterCount + 1):
            band = ds.GetRasterBand(b)
            arr = band.ReadAsArray(srcwin[0], srcwin[1], srcwin[2], srcwin[3])
            arr = np.nan_to_num(arr, nan=0.0)
            band_data.append(arr.astype(np.float32))
            
        X_img, shape = self._extract_features(band_data)
        
        ## Sample
        train_feats = []
        train_z = []
        inv_gt = gdal.InvGeoTransform(gt)
        
        for x, y, z in training_data:
            px = int(inv_gt[0] + inv_gt[1] * x + inv_gt[2] * y) - srcwin[0]
            py = int(inv_gt[3] + inv_gt[4] * x + inv_gt[5] * y) - srcwin[1]
            
            if 0 <= px < shape[1] and 0 <= py < shape[0]:
                idx = py * shape[1] + px
                feats = X_img[idx]
                if np.all(np.isfinite(feats)) and feats[0] > 0:
                    train_feats.append(feats)
                    train_z.append(z)

        if len(train_feats) < 20:
            utils.echo_error_msg(f"Insufficient training samples ({len(train_feats)}).")
            return

        X_train = np.array(train_feats)
        y_train = np.array(train_z)

        ## Train
        utils.echo_msg(f"Training RF on {len(X_train)} points...")
        if self.physics:
            stumpf_vals = X_train[:, 5].reshape(-1, 1)
            lr = LinearRegression()
            lr.fit(stumpf_vals, y_train)
            stumpf_pred = lr.predict(stumpf_vals)
            X_train = np.column_stack([X_train, stumpf_pred])
            
            stumpf_img = X_img[:, 5].reshape(-1, 1)
            stumpf_img_pred = lr.predict(stumpf_img)
            X_img = np.column_stack([X_img, stumpf_img_pred])

        rf = RandomForestRegressor(n_estimators=self.n_trees, n_jobs=-1, random_state=42)
        rf.fit(X_train, y_train)
        
        ## Predict & Write
        valid_mask = np.all(np.isfinite(X_img), axis=1) & (X_img[:, 0] > 0)
        out_z_flat = np.full(shape[0] * shape[1], self.ndv, dtype=np.float32)
        
        if np.any(valid_mask):
            pred_z = rf.predict(X_img[valid_mask])
            out_z_flat[valid_mask] = pred_z
            
        out_grid = out_z_flat.reshape(shape)

        driver = gdal.GetDriverByName(self.fmt)
        ds_out = driver.Create(self.fn, self.xcount, self.ycount, 1, gdal.GDT_Float32, options=self.co)
        ds_out.SetGeoTransform(self.dst_gt)
        ds_out.SetProjection(self.dst_srs)
        band_out = ds_out.GetRasterBand(1)
        band_out.SetNoDataValue(self.ndv)
        band_out.WriteArray(out_grid)
        ds_out = None
        
        utils.echo_msg(f"SDB Grid generated: {self.fn}")
        return self

### End
