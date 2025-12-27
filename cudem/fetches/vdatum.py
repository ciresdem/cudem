### vdatum.py
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
##
## vdatum.py is part of CUDEM
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
## Fetch vertical datum conversion grids from NOAA and PROJ CDN.
##
### Code:

import os
import json
from typing import Optional, Dict, List, Any
from osgeo import ogr
from cudem import utils
from cudem import regions
from cudem.fetches import fetches
from cudem.fetches import FRED

## ==============================================
## Constants
## ==============================================
VDATUM_DATA_URL = 'https://vdatum.noaa.gov/download/data/'
PROJ_CDN_INDEX_URL = 'https://cdn.proj.org/files.geojson'

## Known Tidal References for EPSG lookups
TIDAL_REFERENCES = {
    1089: {'name': 'mllw', 'description': 'Mean Lower Low Water Height', 'grid': 'mllw.gtx'},
    1091: {'name': 'mlw', 'description': 'Mean Low Water Height', 'grid': 'mlw.gtx'},
    5868: {'name': 'mhw', 'description': 'Mean High Water', 'grid': 'mhw.gtx'},
    5869: {'name': 'mhhw', 'description': 'Mean Higher High Water', 'grid': 'mhhw.gtx'},
    5703: {'name': 'tss', 'description': 'NAVD88 tss geoid', 'grid': 'tss.gtx'},
    6641: {'name': 'tss', 'description': 'PRVD02 tss geoid', 'grid': 'tss.gtx'},
    6642: {'name': 'tss', 'description': 'VIVD09 tss geoid', 'grid': 'tss.gtx'},
    5714: {'name': 'tss', 'description': 'to MSL tss geoid', 'grid': 'tss.gtx'},
}

VDATUM_LIST = [
    'TIDAL', 'CRD', 'IGLD85',
    'XGEOID16B', 'XGEOID17B', 'XGEOID18B',
    'XGEOID19B', 'XGEOID20B', 'VERTCON'
]

TIDAL_DATUMS = ['mhw', 'mhhw', 'mlw', 'mllw', 'tss', 'mtl']

CDN_HEADERS = {
    'Host': 'cdn.proj.org',
    'User-Agent': 'Mozilla/5.0 (X11; Linux x86_64; rv:45.0) Gecko/20100101 Firefox/45.0',
    'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8',
    'Accept-Language': 'pt-BR,pt;q=0.8,en-US;q=0.5,en;q=0.3',
    'Accept-Encoding': 'gzip, deflate',
    'Connection': 'keep-alive',
    'Pragma': 'no-cache',
    'Cache-Control': 'no-cache'
}

## ==============================================
## Helper Functions
## ==============================================
def proc_vdatum_inf(vdatum_inf: str, name: Optional[str] = 'vdatum') -> Dict[str, Any]:
    """Process a VDatum INF file to extract region and grid info."""
    
    inf_areas = {}
    
    try:
        with open(vdatum_inf, 'r') as f:
            for line in f:
                line_list = line.split('.')
                if len(line_list) > 1:
                    ## Find where the assignment happens (e.g., key=val)
                    try:
                        eq_val = [i for i, x in enumerate(['=' in x for x in line_list]) if x][0]
                    except IndexError:
                        continue
                        
                    if eq_val != 0:
                        line_key = '.'.join(line_list[:eq_val])
                        if line_key not in inf_areas:
                            inf_areas[line_key] = {}

                        line_val = '.'.join(line_list[eq_val:]).strip()
                        utils.args2dict([line_val], inf_areas[line_key])
    except IOError as e:
        utils.echo_error_msg(f"Error reading VDatum INF file: {e}")
        return {}

    inf_areas_fmt = {}
    for key, data in inf_areas.items():
        if 'minlon' in data:
            out_key = f"{name}_{key}" if name else key
            
            if out_key not in inf_areas_fmt:
                inf_areas_fmt[out_key] = {}

            xmin = utils.x360(float(data.get('minlon', 0)))
            xmax = utils.x360(float(data.get('maxlon', 0)))
            ymin = float(data.get('minlat', 0))
            ymax = float(data.get('maxlat', 0))
            
            inf_areas_fmt[out_key]['region'] = [xmin, xmax, ymin, ymax]
            ## Handle path separators for cross-platform compatibility
            inf_areas_fmt[out_key]['grid'] = data.get('source', '').replace('\\', '/').split('/')[-1]
            
    return inf_areas_fmt


def search_proj_cdn(
        region: Optional['regions.Region'] = None, 
        epsg: Optional[int] = None, 
        crs_name: Optional[str] = None, 
        name: Optional[str] = None,
        verbose: bool = True, 
        cache_dir: str = './'
) -> List[Dict]:
    """Search PROJ CDN for transformation grids."""
    
    cdn_index = utils.make_temp_fn('proj_cdn_files.geojson', cache_dir)
    
    ## Download Index
    try:
        status = fetches.Fetch(
            PROJ_CDN_INDEX_URL, headers=CDN_HEADERS, verbose=verbose
        ).fetch_file(
            cdn_index, timeout=10, read_timeout=10, check_size=False
        )
    except Exception:
        status = -1

    results = []
    
    if status == 0:
        try:
            cdn_driver = ogr.GetDriverByName('GeoJSON')
            cdn_ds = cdn_driver.Open(cdn_index, 0)
            if cdn_ds:
                cdn_layer = cdn_ds.GetLayer()
                
                # Apply Filters
                filter_str = "type != 'HORIZONTAL_OFFSET'"
                if crs_name:
                    filter_str += f" AND (target_crs_name LIKE '%{name.upper()}%' OR source_crs_name LIKE '%{name.upper()}%')"
                elif epsg:
                    filter_str += f" AND (target_crs_code LIKE '%{epsg}%' OR source_crs_code LIKE '%{epsg}%')"
                elif name:
                    filter_str += f" AND name LIKE '%{name}%'"
                
                cdn_layer.SetAttributeFilter(filter_str)

                bounds_geom = region.export_as_geom() if region else None

                for feat in cdn_layer:
                    if bounds_geom:
                        geom = feat.GetGeometryRef()
                        if geom and not bounds_geom.Intersects(geom):
                            continue
                    
                    ## Extract Properties
                    props = json.loads(feat.ExportToJson()).get('properties', {})
                    ## Ensure all fields are captured even if not in properties dict
                    feat_dict = props.copy()
                    results.append(feat_dict)

                cdn_ds = None # Close
        except Exception as e:
            if verbose:
                utils.echo_error_msg(f"Error reading PROJ CDN index: {e}")
        finally:
             if os.path.exists(cdn_index):
                utils.remove_glob(cdn_index)
                
    return results


## ==============================================
## VDATUM Module
## ==============================================
class VDATUM(fetches.FetchModule):
    """NOAA's VDATUM transformation grids

    Fetch vertical datum conversion grids from NOAA and PROJ CDN.
    
    VDatum is designed to vertically transform geospatial data among a 
    variety of tidal, orthometric and ellipsoidal vertical datums.

    https://vdatum.noaa.gov
    https://cdn.proj.org

    Configuration Example:
    < vdatum:datatype=None:gtx=False >
    """
    
    def __init__(self, where: str = '', datatype: Optional[str] = None, gtx: bool = False, epsg: Optional[int] = None, **kwargs):
        super().__init__(name='vdatum', **kwargs)
        
        self.where = [where] if where else []
        self.datatype = datatype
        self.epsg = utils.int_or(epsg)
        self.gtx = gtx
        self.src_srs = 'epsg:4326'

        ## Initialize FRED
        self.FRED = FRED.FRED(name=self.name, verbose=self.verbose)
        self.update_if_not_in_FRED()

        
    def update_if_not_in_FRED(self):
        """Update the fetches module in FRED if it's not already in there."""
        
        self.FRED._open_ds()
        try:
            self.FRED._attribute_filter([f"DataSource = '{self.name}'"])
            if len(self.FRED.layer) == 0:
                self.FRED._close_ds()
                self.update()
        finally:
             if self.FRED.ds is not None:
                self.FRED._close_ds()

                
    def update(self):
        """Update or create the reference vector file."""
        
        self.FRED._open_ds(1)
        
        try:
            for vd in VDATUM_LIST:
                surveys = []
                
                ## Determine URL and INF file
                if vd in ['TIDAL', 'IGLD85', 'CRD']:
                    vd_dl_name = 'DEVAemb12_8301' if vd == 'TIDAL' else vd
                    v_inf = 'tidal_area.inf' if vd == 'TIDAL' else f'{vd}.inf'
                    vd_zip_url = f'{VDATUM_DATA_URL}{vd_dl_name}.zip'
                    
                elif 'XGEOID' in vd:
                    vd_zip_url = f'{VDATUM_DATA_URL}vdatum_{vd}.zip'
                    v_inf = f'{vd.lower()}.inf'
                    
                elif vd == 'VERTCON':
                    vd_zip_url = f'{VDATUM_DATA_URL}vdatum_{vd}.zip'
                    v_inf = 'vcn.inf'
                else:
                    vd_zip_url = f'{VDATUM_DATA_URL}vdatum_{vd}.zip'
                    v_inf = f'{vd.lower()}.inf'

                ## Download Zip
                local_zip = f'{vd}.zip'
                try:
                    status = fetches.Fetch(vd_zip_url, verbose=self.verbose).fetch_file(local_zip)
                except Exception:
                    status = -1
                    
                if status == 0:
                    v_infs = utils.p_unzip(local_zip, ['inf'])
                    v_dict = {}
                    
                    for inf in v_infs:
                        try:
                            # Parse the INF file
                            v_dict = proc_vdatum_inf(inf, name=vd if vd != 'TIDAL' else None)
                            if v_dict: break
                        except Exception:
                            pass

                    ## Augment Dictionary
                    for key in v_dict:
                        v_dict[key]['vdatum'] = vd
                        v_dict[key]['remote'] = vd_zip_url
                    
                    ## Handle Tidal Special Case
                    if vd == 'TIDAL':
                        v_dict_expanded = {}
                        for tidal_key, data in v_dict.items():                        
                            for t in TIDAL_DATUMS:
                                key_ = f'{t}_{tidal_key}'
                                v_dict_expanded[key_] = {
                                    'region': data['region'],
                                    'vdatum': t,
                                    'grid': f'{t}.gtx',
                                    'remote': f'{VDATUM_DATA_URL}{tidal_key}.zip'
                                }
                        v_dict = v_dict_expanded

                    ## Add to Surveys
                    for key, data in v_dict.items():
                        ## Check existence
                        self.FRED._attribute_filter([f"ID = '{key}'"])
                        if self.FRED.layer is None or len(self.FRED.layer) == 0:
                            geom = regions.Region().from_list(data['region']).export_as_geom()
                            if geom:
                                surveys.append({
                                    'Name': data['grid'],
                                    'ID': key,
                                    'Agency': 'NOAA',
                                    'Date': utils.this_date(),
                                    'MetadataLink': "",
                                    'MetadataDate': utils.this_date(),
                                    'DataLink': data['remote'],
                                    'Link': VDATUM_DATA_URL,
                                    'DataType': data['vdatum'],
                                    'DataSource': 'vdatum',
                                    'HorizontalDatum': 4326,
                                    'VerticalDatum': data['vdatum'],
                                    'Info': "",
                                    'geom': geom
                                })
                    
                    ## Cleanup extracted INFs and Zip
                    if v_infs:
                        utils.remove_glob(*v_infs)
                    utils.remove_glob(local_zip)

                self.FRED._add_surveys(surveys)
        except Exception as e:
            utils.echo_error_msg(f"Error updating VDATUM FRED: {e}")                
        finally:
            self.FRED._close_ds()

            
    def run(self):
        """Run the VDatum fetching module."""
        
        ## Search FRED for NOAA VDatum Grids
        where_clauses = []
        if self.datatype:
            where_clauses.append(f"DataType = '{self.datatype}'")
        elif self.epsg and self.epsg in TIDAL_REFERENCES:
             where_clauses.append(f"DataType = '{TIDAL_REFERENCES[self.epsg]['name']}'")

        fred_results = self.FRED._filter(self.region, where_clauses, [self.name])
        
        for surv in fred_results:
            if self.gtx:
                ## Download and Extract specific GTX
                dst_zip = f"{surv['ID']}.zip"
                try:
                    status = fetches.Fetch(
                        surv['DataLink'],
                        callback=self.callback,
                        verbose=self.verbose
                    ).fetch_file(dst_zip)
                except Exception:
                    status = -1
                    
                if status == 0:
                    v_gtxs = utils.p_f_unzip(dst_zip, [surv['Name']])
                    for v_gtx in v_gtxs:
                        if os.path.exists(v_gtx):
                            os.replace(v_gtx, f"{surv['ID']}.gtx")
                    utils.remove_glob(dst_zip)
            else:
                self.add_entry_to_results(
                    surv['DataLink'], 
                    f"{surv['ID']}.zip", 
                    surv['Name'].lower()
                )

        ## Search PROJ CDN for other Grids
        proj_results = search_proj_cdn(
            region=self.region,
            epsg=self.epsg,
            name=self.datatype,
            verbose=self.verbose,
            cache_dir=self._outdir
        )
        
        for result in proj_results:
             self.add_entry_to_results(
                result.get('url'), 
                result.get('name'), 
                result.get('source_crs_code', 'proj_cdn')
            )

        return self

### End
