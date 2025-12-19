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

class WafflesLakes(Waffle):
    """Estimate lake bathymetry.
    
    By default, will return lake bathymetry as depth values (positive down), 
    to get elevations (positive up), set apply_elevations=True.

    -----------
    Parameters:
    
    apply_elevations=[True/False] - use COPERNICUS to apply lake level elevations to output
    min_area=[val] - minimum lake area to consider
    max_area=[val] - maximum lake area to consider
    min_id=[val] - minimum lake ID to consider
    max_id=[val] - maximum lake ID to consider
    depth=[globathy/hydrolakes/val] - obtain the depth value from GloBathy, 
                                      HydroLakes or constant value
    
    < lakes:apply_elevations=False:min_area=None:max_area=None:min_id=None:max_id=None:depth=globathy:elevations=copernicus >
    """
    
    def __init__(
            self,
            apply_elevations=False,
            min_area=None,
            max_area=None,
            min_id=None,
            max_id=None,
            depth='globathy',
            elevations='copernicus',
            **kwargs
    ):
        super().__init__(**kwargs)
        self._mask = None
        self.apply_elevations = apply_elevations
        self.min_area = min_area
        self.max_area = max_area
        self.min_id = min_id
        self.max_id = max_id
        self.depth = depth
        self.elevations = elevations
        self.ds_config = None

        #self.initialize()
        #self.init_lakes()

        
    def _fetch_lakes(self):
        """fetch hydrolakes polygons"""

        this_lakes = fetches.hydrolakes.HydroLakes(
            src_region=self.p_region,
            verbose=self.verbose,
            outdir=self.cache_dir
        )
        this_lakes.run()
        fr = fetches.fetches.fetch_results(this_lakes, check_size=False)
        fr.daemon = True
        fr.start()
        fr.join()

        lakes_shp = None
        lakes_zip = os.path.join(this_lakes._outdir, this_lakes.results[0]['dst_fn'])
        lakes_shps = utils.unzip(lakes_zip, self.cache_dir)
        for i in lakes_shps:
            if i.split('.')[-1] == 'shp':
                lakes_shp = i

        return(lakes_shp)

    
    def _fetch_globathy(self, ids=[]):
        """fetch globathy csv data and process into dict"""
        
        import csv
        
        _globathy_url = 'https://springernature.figshare.com/ndownloader/files/28919991'
        globathy_zip = os.path.join(self.cache_dir, 'globathy_parameters.zip')
        fetches.fetches.Fetch(
            _globathy_url, verbose=self.verbose
        ).fetch_file(globathy_zip, check_size=False)
        globathy_csvs = utils.unzip(globathy_zip, self.cache_dir)        
        globathy_csv = os.path.join(
            self.cache_dir,
            'GLOBathy_basic_parameters/GLOBathy_basic_parameters(ALL_LAKES).csv'
        )
        with open(globathy_csv, mode='r') as globc:
            reader = csv.reader(globc)
            next(reader)
            if len(ids) > 0:
                globd = {}
                for row in reader:
                    if int(row[0]) in ids:
                        globd[int(row[0])] = float(row[-1])
                        ids.remove(int(row[0]))
                        
                    if len(ids) == 0:
                        break
            else:
                globd = {int(row[0]):float(row[-1]) for row in reader}

        return(globd)

    
    def _fetch_gmrt(self, gmrt_region=None):
        """GMRT - Global low-res.
        """

        if gmrt_region is None:
            gmrt_region = self.p_region
        
        this_gmrt = fetches.gmrt.GMRT(
            src_region=gmrt_region,
            verbose=self.verbose,
            layer='topo',
            outdir=self.cache_dir
        )
        this_gmrt.run()

        fr = fetches.fetches.fetch_results(this_gmrt)
        fr.daemon = True
        fr.start()
        fr.join()

        dst_srs = osr.SpatialReference()
        dst_srs.SetFromUserInput(self.dst_srs)
        
        gmrt_tif = os.path.join(this_gmrt._outdir, this_gmrt.results[0]['dst_fn'])
        gmrt_ds = gdalfun.gdal_mem_ds(self.ds_config, name='gmrt', co=self.co)
        gdal.Warp(gmrt_ds, gmrt_tif, dstSRS=dst_srs, resampleAlg=self.sample)
        return(gmrt_ds)

    
    def _fetch_copernicus(self, cop_region=None):
        """copernicus"""

        if cop_region is None:
            cop_region = self.p_region
            
        this_cop = fetches.copernicus.CopernicusDEM(
            src_region=cop_region,
            verbose=self.verbose,
            datatype='1',
            outdir=self.cache_dir
        )
        this_cop.run()

        fr = fetches.fetches.fetch_results(this_cop, check_size=False)
        fr.daemon = True
        fr.start()
        fr.join()

        dst_srs = osr.SpatialReference()
        dst_srs.SetFromUserInput(self.dst_srs)
        
        cop_ds = gdalfun.gdal_mem_ds(self.ds_config, name='copernicus', co=self.co)
        [gdal.Warp(
            cop_ds,
            os.path.join(this_cop._outdir, cop_result['dst_fn']),
            dstSRS=dst_srs,
            resampleAlg=self.sample
        ) for cop_result in this_cop.results]
        
        return(cop_ds)

    
    def generate_mem_ogr(self, geom, srs):
        """Create temporary polygon vector layer from feature geometry 
        so that we can rasterize it (Rasterize needs a layer)
        """
        
        ds = ogr.GetDriverByName('MEMORY').CreateDataSource('')
        Layer = ds.CreateLayer('', geom_type=ogr.wkbPolygon, srs=srs)
        outfeature = ogr.Feature(Layer.GetLayerDefn())
        outfeature.SetGeometry(geom)
        Layer.SetFeature(outfeature)

        return(ds)


    ## Adapted from GLOBathy
    def apply_calculation(self, shore_distance_arr, lake_depths_arr, shore_arr=None):
        """
        Apply distance calculation, which is each pixel's distance to shore, multiplied
        by the maximum depth, all divided by the maximum distance to shore. This provides
        a smooth slope from shore to lake max depth.

        shore_distance - Input numpy array containing array of distances to shoreline.
            Must contain positive values for distance away from shore and 0 elsewhere.
        max_depth - Input value with maximum lake depth.
        NoDataVal - value to assign to all non-lake pixels (zero distance to shore).
        """

        labels, nfeatures = ndimage.label(shore_distance_arr)
        nfeatures = np.arange(1, nfeatures +1)
        maxes = ndimage.maximum(shore_distance_arr, labels, nfeatures)
        max_dist_arr = np.zeros(np.shape(shore_distance_arr))
        with utils.ccp(
                total=len(nfeatures),
                desc='applying labels.',
                leave=self.verbose
        ) as pbar:
            for n, x in enumerate(nfeatures):
                pbar.update()
                max_dist_arr[labels==x] = maxes[x-1]
        
        max_dist_arr[max_dist_arr == 0] = np.nan
        bathy_arr = (shore_distance_arr * lake_depths_arr) / max_dist_arr
        bathy_arr[bathy_arr == 0] = np.nan

        if shore_arr is None \
           or shore_arr.size == 0 \
           or shore_arr[~np.isnan(bathy_arr)].size == 0 \
           or shore_arr[~np.isnan(bathy_arr)].max() == 0:
            utils.echo_warning_msg(
                'invalid shore array, using default shore value of zero'
            )
            bathy_arr = 0 - bathy_arr
        else:
            bathy_arr = shore_arr - bathy_arr

        utils.echo_msg('applied shore elevations to lake depths')
        bathy_arr[np.isnan(bathy_arr)] = 0
        return(bathy_arr)    

    
    def _init(self):
        self.wgs_region = self.p_region.copy()
        if self.dst_srs is not None:
            self.wgs_region.warp('epsg:4326')            
        else:
            self.dst_srs = 'epsg:4326'
        
        #self.p_region.buffer(pct=2)
        
        lakes_shp = self._fetch_lakes()
        self.lk_ds = ogr.Open(lakes_shp, 1)
        self.lk_layer = self.lk_ds.GetLayer()

        ## filter layer to region
        filter_region = self.p_region.copy()
        self.lk_layer.SetSpatialFilter(filter_region.export_as_geom())

        ## filter by ID
        if self.max_id is not None:
            self.lk_layer.SetAttributeFilter(
                f'Hylak_id < {self.max_id}'
            )
            
        if self.min_id is not None:
            self.lk_layer.SetAttributeFilter(
                f'Hylak_id > {self.min_id}'
            )

        ## filter by Area
        if self.max_area is not None:
            self.lk_layer.SetAttributeFilter(
                f'Lake_area < {self.max_area}'
            )
            
        if self.min_area is not None:
            self.lk_layer.SetAttributeFilter(
                f'Lake_area > {self.min_area}'
            )
            
        lk_features = self.lk_layer.GetFeatureCount()
        if lk_features == 0:
            utils.echo_error_msg('no lakes found in region')
            return(-1)
            #return(self)

        ## get lake ids and globathy depths
        self.lk_ids = []
        [self.lk_ids.append(feat.GetField('Hylak_id')) for feat in self.lk_layer]
        utils.echo_msg('using Lake IDS: {}'.format(self.lk_ids))
        
        lk_regions = self.p_region.copy()
        with utils.ccp(
                total=len(self.lk_layer),
                desc=f'processing {lk_features} lakes',
                leave=self.verbose
        ) as pbar:            
            for lk_f in self.lk_layer:
                pbar.update()
                this_region = regions.Region()
                lk_geom = lk_f.GetGeometryRef()
                lk_wkt = lk_geom.ExportToWkt()
                this_region.from_list(
                    ogr.CreateGeometryFromWkt(lk_wkt).GetEnvelope()
                )
                lk_regions = regions.regions_merge(lk_regions, this_region)

        while not regions.regions_within_ogr_p(self.p_region, lk_regions):
            utils.echo_msg(
                ('buffering region by 2 percent to gather all '
                 f'lake boundaries...{self.p_region}')
            )
            self.p_region.buffer(pct=2, x_inc=self.xinc, y_inc=self.yinc)

        return(0)

    
    def run(self):
        ## fetch and initialize the shoreline data
        cop_band = None
        cop_arr = None
        ds_config = gdalfun.gdal_infos(self.stack)
        if self.elevations == 'copernicus':
            cop_ds = self._fetch_copernicus(cop_region=self.p_region)
            cop_band = cop_ds.GetRasterBand(1)
        elif self.elevations == 'gmrt':
            cop_ds = self._fetch_gmrt(gmrt_region=self.p_region)
            cop_band = cop_ds.GetRasterBand(1)
        elif utils.float_or(self.elevations) is not None: # single value
            cop_band = None
            cop_arr = np.zeros((ds_config['nx'], ds_config['ny']))
            cop_arr[:] = self.elevations
        elif self.elevations == 'self': # from stacks, interpolated
            cop_band = None
            tmp_arr = gdalfun.gdal_get_array(self.stack, band=1)[0]
            cop_arr = gdalfun.generate_mem_ds(
                gdalfun.gdal_infos(self.stack),
                band_data=tmp_arr,
                srcwin=None,
                return_array=True,
                interpolation='flats'
            )
        elif os.path.exists(self.elevations): # from input raster
            elev_ds = gdal.Open(self.elevations)
            if elev_ds is not None:
                dst_srs = osr.SpatialReference()
                dst_srs.SetFromUserInput(self.dst_srs)
                cop_ds = gdalfun.gdal_mem_ds(
                    ds_config,
                    name='cop',
                    co=self.co
                )
                gdal.Warp(
                    cop_ds,
                    elev_ds,
                    dstSRS=dst_srs,
                    resampleAlg=self.sample
                )
                cop_band = cop_ds.GetRasterBand(1)

        if cop_band is not None:
            cop_arr = cop_band.ReadAsArray()
                
        ## initialize the tmp datasources
        prox_ds = gdalfun.gdal_mem_ds(ds_config, name='prox', co=self.co)
        msk_ds = gdalfun.gdal_mem_ds(ds_config, name='msk', co=self.co)
        msk_band = None
        globd = None
        
        if len(self.lk_ids) == 0:
            return(self)
        
        if self.depth == 'globathy':
            globd = self._fetch_globathy(ids=self.lk_ids[:])
            ## rasterize hydrolakes using id
            gdal.RasterizeLayer(
                msk_ds,
                [1],
                self.lk_layer,
                options=["ATTRIBUTE=Hylak_id"],
                callback=gdal.TermProgress
            )
            msk_ds.FlushCache()
            msk_band = msk_ds.GetRasterBand(1)
            msk_band.SetNoDataValue(ds_config['ndv'])

            ## assign max depth from globathy
            msk_arr = msk_band.ReadAsArray()
            with utils.ccp(
                    total=len(self.lk_ids),
                    desc='Assigning Globathy Depths to rasterized lakes...',
                    leave=self.verbose
            ) as pbar:
                
                for n, this_id in enumerate(self.lk_ids):
                    depth = globd[this_id]
                    msk_arr[msk_arr == this_id] = depth
                    pbar.update()
            
        elif self.depth == 'hydrolakes':
            gdal.RasterizeLayer(
                msk_ds,
                [1],
                self.lk_layer,
                options=["ATTRIBUTE=Depth_avg"],
                callback=gdal.TermProgress
            )
            msk_ds.FlushCache()
            msk_band = msk_ds.GetRasterBand(1)
            msk_band.SetNoDataValue(ds_config['ndv'])
            msk_arr = msk_band.ReadAsArray()
            
        elif utils.float_or(self.depth) is not None:
            msk_arr = np.zeros((ds_config['nx'], ds_config['ny']))
            msk_arr[:] = self.depth            
        else:            
            msk_arr = np.ones((ds_config['nx'], ds_config['ny']))

        ## calculate proximity of lake cells to shore
        if msk_band is None:
            gdal.RasterizeLayer(
                msk_ds,
                [1],
                self.lk_layer,
                options=["ATTRIBUTE=Hylak_id"],
                callback=gdal.TermProgress
            )
            msk_ds.FlushCache()
            msk_band = msk_ds.GetRasterBand(1)
            msk_band.SetNoDataValue(ds_config['ndv'])
            
        self.lk_ds = None
        prox_band = prox_ds.GetRasterBand(1)
        proximity_options = ["VALUES=0", "DISTUNITS=PIXEL"]
        gdal.ComputeProximity(
            msk_band,
            prox_band,
            options=proximity_options,
            callback=gdal.TermProgress
        )
        prox_arr = prox_band.ReadAsArray()

        ## apply calculation from globathy
        utils.echo_msg('Calculating simulated lake depths...')

        if self.depth == 'flatten':
            cop_arr[prox_arr == 0] = 0
            bathy_arr = cop_arr
        else:
            bathy_arr = self.apply_calculation(
                prox_arr,
                msk_arr,
                shore_arr=cop_arr,
            )
            
        bathy_arr[bathy_arr == 0] = self.ndv            
        gdalfun.gdal_write(
            bathy_arr, '{}.tif'.format(self.name), ds_config,
        )            

        prox_ds = msk_ds = cop_ds = None
        return(self)


### End
