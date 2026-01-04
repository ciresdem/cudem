### ml_interp.py
##
## Copyright (c) 2025 - 2026 Regents of the University of Colorado
##
## ml_interp.py is part of CUDEM
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
## Machine Learning Interpolation module using Scikit-Learn.
## Treats interpolation as a regression problem: Z = f(X, Y).
##
### Code:

import numpy as np
from osgeo import gdal

from cudem import utils
from cudem import gdalfun
from cudem.waffles.waffles import Waffle

## Optional Dependency
try:
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.neighbors import KNeighborsRegressor
    from sklearn.neural_network import MLPRegressor
    from sklearn.preprocessing import StandardScaler
    from sklearn.pipeline import make_pipeline
    HAS_SKLEARN = True
except ImportError:
    HAS_SKLEARN = False

    
class WafflesML(Waffle):
    """MACHINE LEARNING DEM via Scikit-Learn Regression.
    
    Interpolates elevation data by training a regression model on 
    spatial coordinates (X, Y) -> Z.

    Algorithms:
      - 'rf': Random Forest Regressor (Good for complex, non-linear terrains)
      - 'knn': K-Nearest Neighbors (Similar to IDW but often faster/adaptive)
      - 'mlp': Multi-Layer Perceptron (Neural Network)

    Parameters:
    -----------
    algo (str) : Algorithm to use ['rf', 'knn', 'mlp']. Default: 'rf'.
    max_points (int) : Max training points (decimated if exceeded). Default: 50000.
    trees (int) : Number of trees (Random Forest only). Default: 100.
    neighbors (int) : Number of neighbors (KNN only). Default: 5.
    hidden_layers (tuple) : Neurons per layer (MLP only). Default: (100, 50).
    chunk_size (int) : Processing chunk size in pixels.

    < ml_interp:algo=rf:trees=100:max_points=50000 >
    """
    
    def __init__(
            self,
            algo='rf',
            max_points=50000,
            trees=100,
            neighbors=5,
            hidden_layers='100/50',
            chunk_size=None,
            **kwargs
    ):
        super().__init__(**kwargs)
        self.algo = algo.lower()
        self.max_points = utils.int_or(max_points, 50000)
        self.trees = utils.int_or(trees, 100)
        self.neighbors = utils.int_or(neighbors, 5)
        
        ## Parse hidden layers string "100/50" -> (100, 50)
        if isinstance(hidden_layers, str):
            self.hidden_layers = tuple([int(x) for x in hidden_layers.split('/')])
        else:
            self.hidden_layers = (100, 50)
            
        self.chunk_size = chunk_size
        self.chunk_step = None

    def _get_model(self):
        """Initialize the requested Scikit-Learn model."""
        
        if self.algo == 'rf':
            ## Random Forest: Robust, handles non-linearities well, no scaling needed usually
            return RandomForestRegressor(
                n_estimators=self.trees,
                n_jobs=1, #set to -1 for all cores, waffles already starts mp
                random_state=42
            )
            
        elif self.algo == 'knn':
            ## KNN: Distance-based, requires scaling for Lat/Lon
            ## We wrap it in a pipeline to auto-scale inputs
            return make_pipeline(
                StandardScaler(),
                KNeighborsRegressor(
                    n_neighbors=self.neighbors,
                    weights='distance',
                    n_jobs=1 # set to -1 for all-cores, waffles already starts mp
                )
            )
            
        elif self.algo == 'mlp':
            ## MLP: Neural Net, definitely requires scaling
            return make_pipeline(
                StandardScaler(),
                MLPRegressor(
                    hidden_layer_sizes=self.hidden_layers,
                    activation='relu',
                    solver='adam',
                    max_iter=500,
                    random_state=42
                )
            )
        else:
            utils.echo_warning_msg(f"Unknown algorithm '{self.algo}', falling back to 'knn'")
            return KNeighborsRegressor(n_neighbors=5, weights='distance')

        
    def run(self):
        if not HAS_SKLEARN:
            utils.echo_error_msg("Scikit-Learn must be installed to use the ML_INTERP module.")
            return self

        if self.verbose:
            utils.echo_msg(
                f'Generating ML Grid ({self.algo}) @ {self.ycount}/{self.xcount} '
                f'using max {self.max_points} training points...'
            )

        if self.chunk_size is None:
            n_chunk = 512
        else:
            n_chunk = self.chunk_size

        ## Load Data from Stack
        ## Read the stack to get training samples (X, Y, Z)
        with gdalfun.gdal_datasource(self.stack) as stack_ds:
            points_band = stack_ds.GetRasterBand(1)
            points_no_data = points_band.GetNoDataValue()
            
            ## Read into memory (Caution with massive files)
            points_array = points_band.ReadAsArray()
            
            ## Extract valid data mask
            valid_mask = points_array != points_no_data
            if np.count_nonzero(valid_mask) == 0:
                utils.echo_warning_msg("No valid data points found in stack.")
                return self

            ## Convert to coordinate arrays
            gt = stack_ds.GetGeoTransform()
            py, px = np.nonzero(valid_mask)
            
            z_train = points_array[py, px]
            x_train, y_train = utils._pixel2geo(px, py, gt, node='pixel')
            
            ## Decimate training data if too large
            if len(z_train) > self.max_points:
                if self.verbose:
                    utils.echo_msg(f"Decimating input from {len(z_train)} to {self.max_points} points.")
                indices = np.random.choice(len(z_train), self.max_points, replace=False)
                x_train = x_train[indices]
                y_train = y_train[indices]
                z_train = z_train[indices]

            ## Feature Matrix (X) and Target Vector (y)
            ## X features: [Longitude, Latitude]
            ## Note: For complex terrain, you could add features here (e.g. distance to center)
            X_features = np.column_stack((x_train, y_train))
            y_target = z_train

            ## Clean memory
            points_array = None

            ## Train Model
            model = self._get_model()
            if self.verbose:
                utils.echo_msg(f"Training {self.algo.upper()} model...")
            
            try:
                model.fit(X_features, y_target)
            except Exception as e:
                utils.echo_error_msg(f"Model training failed: {e}")
                return self

            ## Predict / Interpolate (Chunked)
            ## We generate a grid of coordinates for each chunk and predict Z
            
            ## Prepare Output Dataset
            driver = gdal.GetDriverByName(self.fmt)
            dst_ds = driver.Create(
                self.fn,
                self.xcount,
                self.ycount,
                1,
                gdal.GDT_Float32,
                options=self.co
            )
            
            dst_ds.SetGeoTransform(self.dst_gt)
            dst_ds.SetProjection(gdalfun.osr_wkt(self.dst_srs))
            
            elev_band = dst_ds.GetRasterBand(1)
            elev_band.SetNoDataValue(self.ndv)
            elev_band.SetDescription(f"ML Interpolation ({self.algo})")

            if self.verbose:
                utils.echo_msg("Predicting surface...")

            for srcwin in utils.yield_srcwin(
                    (self.ycount, self.xcount),
                    n_chunk=n_chunk,
                    msg='ML Prediction',
                    verbose=self.verbose
            ):
                ## Calculate coordinate vectors for this chunk
                x_off, y_off, x_size, y_size = srcwin
                
                ## Get Pixel Centers
                x_start, y_start = utils._pixel2geo(x_off, y_off, self.dst_gt, node='pixel')
                
                ## Create Grid Mesh for Chunk
                ## X coordinates (1D)
                grid_x = np.linspace(
                    x_start, 
                    x_start + (x_size * self.dst_gt[1]), 
                    x_size, endpoint=False
                )
                
                ## Y coordinates (1D)
                grid_y = np.linspace(
                    y_start, 
                    y_start + (y_size * self.dst_gt[5]), 
                    y_size, endpoint=False
                )
                
                ## Meshgrid -> Flatten -> Feature Matrix
                ## We use 'xy' indexing for coordinate consistency
                GX, GY = np.meshgrid(grid_x, grid_y)
                
                ## Flatten to (N_pixels, 2) features
                X_pred = np.column_stack((GX.ravel(), GY.ravel()))
                
                ## Predict
                try:
                    z_pred = model.predict(X_pred)
                    
                    ## Reshape back to chunk dimensions (Rows, Cols)
                    z_chunk = z_pred.reshape((y_size, x_size))
                    
                    ## Write to disk
                    elev_band.WriteArray(z_chunk, x_off, y_off)
                    
                except Exception as e:
                    utils.echo_warning_msg(f"Prediction failed @ {srcwin}: {e}")

            dst_ds = None
            
        return self

### End
