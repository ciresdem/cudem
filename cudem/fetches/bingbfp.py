### bingbfp.py
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
##
## bingbfp.py is part of CUDEM
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
## Fetch and process Microsoft Bing Building Footprints.
##
### Code:

import os
import sys
import csv
import mercantile
import argparse
from typing import List, Optional, Union, Tuple, Any
from osgeo import ogr
from cudem import utils
from cudem import gdalfun
from cudem import regions
from cudem import factory
from cudem.fetches import fetches

## ==============================================
## Constants
## ==============================================
BING_CSV_URL = 'https://minedbuildings.z5.web.core.windows.net/global-buildings/dataset-links.csv'
BING_HEADERS = {
    'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:89.0) Gecko/20100101 Firefox/89.0',
    'referer': 'https://lz4.overpass-api.de/'
}

## ==============================================
## Bing Building Footprints Fetch Module
## ==============================================
class BingBFP(fetches.FetchModule):
    """Bing Building Footprints Fetch Module
    
    Fetches building footprint data from Microsoft's Global ML Building Footprints.
    https://github.com/microsoft/GlobalMLBuildingFootprints

    Configuration Example:
    < bingbfp >
    """
    
    def __init__(self, **kwargs):
        super().__init__(name='bingbfp', **kwargs)
        self.headers = BING_HEADERS

    def run(self):
        """Run the Bing BFP fetches module.
        
        Calculates QuadKeys for the region and parses the Microsoft CSV to find
        relevant download links.
        """
        
        if self.region is None:
            return []

        ## Calculate QuadKeys for Zoom Level 9 covering the region
        quad_keys = set()
        for tile in list(mercantile.tiles(*self.region.xy_extent, zooms=9)):
            quad_keys.add(int(mercantile.quadkey(tile)))
            
        bing_csv_path = os.path.join(self._outdir, os.path.basename(BING_CSV_URL))
        
        ## Download the dataset links CSV
        try:
            status = fetches.Fetch(BING_CSV_URL, verbose=self.verbose).fetch_file(bing_csv_path)
        except Exception as e:
            utils.echo_error_msg(f"Error fetching Bing CSV: {e}")
            status = -1
            
        if status == 0 and os.path.exists(bing_csv_path):
            with open(bing_csv_path, mode='r', newline='', encoding='utf-8') as bc:
                reader = csv.reader(bc)
                next(reader) # Skip header
                
                ## Filter rows: [Location, QuadKey, Url]
                ## CSV format expected: Location, QuadKey, Url, Size, etc.
                for row in reader:
                    try:
                        qk = int(row[1])
                        if qk in quad_keys:
                            url = row[2]
                            location = row[0]
                            filename = f"{url.split('/')[-1]}"
                            ## Prefixing filename with QuadKey and Location to avoid collisions
                            dst_fn = f"{location}_{qk}_{filename}"
                            
                            self.add_entry_to_results(url, dst_fn, 'bing')
                    except ValueError:
                        continue
        else:
            utils.echo_error_msg('Could not fetch BING dataset-links.csv')

        return self

## ==============================================
## Workflow Class
## ==============================================
class BingBuildings:
    """High-level class to manage fetching and processing of Bing Buildings."""

    def __init__(self, region=None, verbose: bool = True, attempts: int = 5, n_threads: int = 5, cache_dir: str = '.'):
        self.region = region
        self.verbose = verbose
        self.attempts = attempts
        self.n_threads = n_threads
        self.cache_dir = cache_dir
        self.fetch_module: Optional[BingBFP] = None

        
    def __call__(self, out_fn: str = None, return_geom: bool = True, overwrite: bool = False) -> Union[str, List[ogr.Geometry], None]:
        """Main entry point to run the fetch and process pipeline."""
        
        if self.region is None or not self.region.is_valid():
            utils.echo_error_msg(f'{self.region} is an invalid region')
            return None

        ## Determine output filename
        if not return_geom:
            if out_fn is None or not isinstance(out_fn, str):
                out_fn = utils.append_fn('bing_buildings', self.region, 1) + '.gpkg'
                
            if not overwrite and os.path.exists(out_fn):
                if self.verbose:
                    utils.echo_msg(f'Output file {out_fn} exists, skipping...')
                return out_fn
            elif overwrite:
                utils.remove_glob(out_fn)
        else:
            out_fn = None

        #utils.echo_msg(out_fn)
        ## Run Processing
        out_fn, bing_geoms = self.process(out_fn)
                
        if return_geom:            
            return bing_geoms
        else:
            return out_fn

        
    def init_fetch(self):
        """Initialize the FetchModule."""
        
        self.fetch_module = BingBFP(
            src_region=self.region,
            verbose=self.verbose,
            outdir=self.cache_dir
        )

        
    def fetch(self) -> fetches.fetch_results:
        """Execute the fetching process."""
        
        if self.fetch_module is None:
            self.init_fetch()
            
        self.fetch_module.run()
        
        ## Threaded downloader
        fr = fetches.fetch_results(
            self.fetch_module, 
            check_size=False, 
            attempts=self.attempts, 
            n_threads=self.n_threads
        )
        fr.daemon = True
        fr.start()
        fr.join()
        return fr

    
    def process(self, out_fn: Optional[str]) -> Tuple[Optional[str], List[ogr.Geometry]]:
        """Process downloaded files into a unified geometry list or file."""
        
        bldg_geoms = []
        
        ## Ensure data is fetched
        results_handler = self.fetch()
        ## Note: self.fetch() runs the fetcher. access results from self.fetch_module.results        
        ## We need to access the results that were populated in the module
        results = self.fetch_module.results

        with utils.ccp(total=len(results), desc='Processing Buildings', leave=self.verbose) as pbar:
            for entry in results:
                ## entry is a dict: {'url':..., 'dst_fn':..., 'data_type':...}
                ## Check actual file existence since fetch_results updates status in place or we check local path
                
                local_path = os.path.join(self.fetch_module._outdir, entry['dst_fn'])
                if os.path.exists(local_path):
                    try:
                        ## Unzip
                        # utils.gunzip returns the path to the unzipped file
                        bing_gj = utils.gunzip(local_path, self.cache_dir)
                        
                        ## Rename/Prepare for OGR
                        ## Ensure it has a recognized extension if gunzip didn't provide one
                        if not bing_gj.endswith('.geojson') and not bing_gj.endswith('.json'):
                            new_name = bing_gj + '.geojson'
                            os.rename(bing_gj, new_name)
                            bing_gj = new_name

                        ## Open with OGR
                        bldg_ds = ogr.Open(bing_gj, 0) # 0 = ReadOnly
                        if bldg_ds:
                            bldg_layer = bldg_ds.GetLayer()

                            ## Spatial Filter
                            if self.region is not None:
                                _boundsGeom = self.region.export_as_geom()                            
                                bldg_layer.SetSpatialFilter(_boundsGeom)
                            
                            ## Union Geometries
                            bldg_geom = gdalfun.ogr_union_geom(bldg_layer, verbose=False)
                            if bldg_geom:
                                bldg_geoms.append(bldg_geom)
                            
                            bldg_ds = None # Close DS

                        ## Cleanup unzipped file to save space
                        utils.remove_glob(bing_gj) 

                    except Exception as e:
                        utils.echo_error_msg(f'Could not process Bing BFP file {local_path}: {e}')

                    ## Cleanup downloaded gz file
                    utils.remove_glob(local_path)

                pbar.update()

        ## Output to file if requested
        if out_fn is not None and bldg_geoms:
            gdalfun.ogr_geoms2ogr(bldg_geoms, out_fn, ogr_format='GPKG')
        
        return out_fn, bldg_geoms        

## Alias
bingBuildings = BingBuildings


## ==============================================
## Command-line Interface (CLI)
## $ bing_bfp
##
## bing_bfp cli
## ==============================================
def bing_bfp_cli():
    parser = argparse.ArgumentParser(
        description="Fetch and process Microsoft Bing Building Footprints.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="CUDEM home page: <http://cudem.colorado.edu>"
    )

    parser.add_argument(
        '-R', '--region',
        required=True,
        help=("Restrict processing to the desired REGION \n"
              "Where a REGION is xmin/xmax/ymin/ymax[/zmin/zmax[/wmin/wmax]]\n"
              "OR an OGR-compatible vector file with regional polygons.")
    )
    parser.add_argument(
        '-o', '--output',
        help="Output filename (GPKG). If not provided, a default name is generated."
    )
    parser.add_argument(
        '-c', '--cache-dir',
        default='.',
        help="Directory to cache downloaded files."
    )
    parser.add_argument(
        '-t', '--threads',
        type=int,
        default=5,
        help="Number of download threads."
    )
    parser.add_argument(
        '-a', '--attempts',
        type=int,
        default=5,
        help="Number of retry attempts for downloads."
    )
    parser.add_argument(
        '-q', '--quiet',
        action='store_true',
        help="Suppress verbose output."
    )
    parser.add_argument(
        '--overwrite',
        action='store_true',
        help="Overwrite existing output file."
    )

    fixed_argv = factory.fix_argparse_region(sys.argv[1:])
    
    args = parser.parse_args(fixed_argv)

    ## Create Region
    try:
        region = regions.parse_cli_region([args.region], not args.quiet)[0]
        if not region.is_valid():
            raise ValueError("Region is invalid.")
    except Exception as e:
        utils.echo_error_msg(f"Error parsing region: {e}")
        sys.exit(1)

    ## Initialize BingBuildings Workflow
    bing = BingBuildings(
        region=region,
        verbose=not args.quiet,
        attempts=args.attempts,
        n_threads=args.threads,
        cache_dir=args.cache_dir
    )

    ## Run
    try:
        out_file = bing(
            out_fn=args.output,
            return_geom=False,
            overwrite=args.overwrite
        )

        if out_file and os.path.exists(out_file):
            utils.echo_msg(f"Successfully generated: {out_file}")
        else:
            utils.echo_error_msg("No output generated (data might not exist in this region).")

    except Exception as e:
        utils.echo_error_msg(f"Process failed: {e}")
        sys.exit(1)

        
if __name__ == '__main__':
    bing_bfp_cli()
    
### End
