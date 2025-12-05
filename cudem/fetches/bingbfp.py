### bingbfp.py
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
##
## fetches.py is part of CUDEM
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
##
### Code:

import os
import mercantile
import csv
from tqdm import tqdm
from osgeo import ogr
from cudem import utils
from cudem import gdalfun
from cudem.fetches import fetches


class bingBuildings:
    def __init__(self, region=None, verbose=True, attempts=5, n_threads=5, cache_dir='.'):
        self.region = region
        self.verbose = verbose
        self.attempts = attempts
        self.n_threads = n_threads
        self.cache_dir = cache_dir
        
        self.this_bldg = None

        
    def __call__(self, out_fn=None, return_geom=True, overwrite=False):
        if self.region is None or not self.region.valid_p():
            utils.echo_error_msg(f'{self.region} is an invalid region')
            return(None)

        if not return_geom:
            if out_fn is None or not isinstance(out_fn, str):
                out_fn = '{}.gpkg'.format(utils.append_fn('bing_buildings', self.region, 1))
                
            if not overwrite:
                if os.path.exists(out_fn):
                    return(out_fn)
            else:
                utils.remove_glob(out_fn)
        else:
            out_fn = None

        out_fn, bing_geoms = self.process(out_fn)
                
        if return_geom:            
            return(bing_geoms)
        else:
            return(out_fn)
                
                
    def init_fetch(self):
        self.this_bldg = BingBFP(
            src_region=self.region,
            verbose=self.verbose,
            outdir=self.cache_dir
        )

        
    def fetch(self):
        if self.this_bldg is None:
            self.init_fetch()
            
        self.this_bldg.run()
        fr = fetches.fetch_results(
            self.this_bldg, check_size=False, attempts=self.attempts, n_threads=self.n_threads
        )
        fr.daemon=True
        fr.start()
        fr.join()
        return(fr)


    def process(self, out_fn):
        bldg_geoms = []
        if self.this_bldg is None:
            self.this_bldg = self.fetch()

        with tqdm(
                total=len(self.this_bldg.results),
                desc='processing buildings',
                leave=self.verbose
        ) as pbar:
            for n, bing_result in enumerate(self.this_bldg.results):
                if bing_result[-1] == 0:
                    bing_gz = bing_result[1]
                    try:
                        bing_gj = utils.gunzip(bing_gz, self.cache_dir)
                        os.rename(bing_gj, bing_gj + '.geojson')
                        bing_gj = bing_gj + '.geojson'
                        bldg_ds = ogr.Open(bing_gj, 0)
                        bldg_layer = bldg_ds.GetLayer()
                        bldg_geom = gdalfun.ogr_union_geom(
                            bldg_layer, verbose=False
                        )
                        bldg_geoms.append(bldg_geom)
                        bldg_ds = None
                        #utils.remove_glob(bing_gj)
                    except Exception as e:
                        utils.echo_error_msg(f'could not process bing bfp, {e}')

                    #utils.remove_glob(bing_gz)

                pbar.update()

        if out_fn is not None:
            gdalfun.ogr_geoms2ogr(bldg_geoms, out_fn, ogr_format='GPKG')
        
        return(out_fn, bldg_geoms)        


def fetch_buildings(region=None, verbose=True, cache_dir='./'):
    """fetch building footprints from BING"""

    if region.valid_p():
        if verbose:
            utils.echo_msg(
                f'fetching buildings for region {region}'
            )

        this_bldg = BingBFP(
            src_region=region,
            verbose=verbose,
            outdir=cache_dir
        )
        this_bldg.run()
        fr = fetches.fetch_results(this_bldg)
        #, check_size=False)
        fr.daemon=True
        fr.start()
        fr.join()
        return(fr)

    return(None)


def process_buildings(this_bing, region=None, return_geom=False, verbose=True, cache_dir='./', overwrite=False):
    if not return_geom:
        out_fn = utils.append_fn('bing_buildings.gpkg', region, 1)
        if overwrite:
            if os.path.exists(out_fn):
                return(out_fn)
            else:
                utils.remove_glob(out_fn)
            
    bldg_geoms = []
    if this_bing is not None:
        with tqdm(
                total=len(this_bing.results),
                desc='processing buildings',
                leave=verbose
        ) as pbar:
            for n, bing_result in enumerate(this_bing.results):
                if bing_result[-1] == 0:
                    bing_gz = bing_result[1]
                    try:
                        bing_gj = utils.gunzip(bing_gz, cache_dir)
                        os.rename(bing_gj, bing_gj + '.geojson')
                        bing_gj = bing_gj + '.geojson'
                        bldg_ds = ogr.Open(bing_gj, 0)
                        bldg_layer = bldg_ds.GetLayer()
                        bldg_geom = gdalfun.ogr_union_geom(
                            bldg_layer, verbose=False
                        )
                        bldg_geoms.append(bldg_geom)
                        bldg_ds = None
                        utils.remove_glob(bing_gj)
                    except Exception as e:
                        utils.echo_error_msg(f'could not process bing bfp, {e}')

                    utils.remove_glob(bing_gz)

                pbar.update()

    if return_geom:
        return(bldg_geoms)
    else:
        return(gdalfun.geoms2ogr(bldg_geoms, out_fn, ogr_format='GPKG'))

    
## BING Building Footprints
class BingBFP(fetches.FetchModule):
    """Bing Building Footprints

    https://github.com/microsoft/GlobalMLBuildingFootprints
    """
    
    def __init__(self, **kwargs):
        super().__init__(name='bingbfp', **kwargs)

        ## The various BING-BFP URLs
        self._bing_bfp_csv = ('https://minedbuildings.z5.web.core.windows.net/'
                              'global-buildings/dataset-links.csv')

        ## Set the user-agent and referer
        self.headers = {
            'User-Agent': ('Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:89.0) '
                           'Gecko/20100101 Firefox/89.0'),
            'referer': 'https://lz4.overpass-api.de/'
        }

        
    def run(self):
        """Run the Bing BFP fetches module"""
        
        if self.region is None:
            return([])

        bbox = self.region.format('bbox')
        quad_keys = set()
        for tile in list(
                mercantile.tiles(
                    self.region.xmin, self.region.ymin,
                    self.region.xmax, self.region.ymax, zooms=9
                )
        ):
            quad_keys.add(int(mercantile.quadkey(tile)))
            
        quad_keys = list(quad_keys)
        #utils.echo_msg('The input area spans {} tiles: {}'.format(len(quad_keys), quad_keys))
        #utils.echo_msg('The input area spans {} tiles.'.format(len(quad_keys)))
        bing_csv = os.path.join(self._outdir, os.path.basename(self._bing_bfp_csv))
        try:
            status = fetches.Fetch(self._bing_bfp_csv, verbose=self.verbose).fetch_file(bing_csv)
        except:
            status = -1
            
        if status == 0 and os.path.exists(bing_csv):
            with open(bing_csv, mode='r') as bc:
                reader = csv.reader(bc)
                next(reader)
                bd = [[row[2], row[1], row[0]] for row in reader if int(row[1]) in quad_keys]

            [self.add_entry_to_results(
                line[0], '{}_{}_{}'.format(line[2], line[1], os.path.basename(line[0])), 'bing'
            ) for line in bd]
        else:
            utils.echo_error_msg('could not fetch BING dataset-links.csv')

### End
