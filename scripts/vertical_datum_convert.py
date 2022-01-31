#!/usr/bin/env python
### vertical_datum_convert.py
##
## Copyright (c) 2022  Regents of the University of Colorado
##
## Permission is hereby granted, free of charge, to any person obtaining a copy 
## of this software and associated documentation files (the "Software"), to deal 
## in the Software without restriction, including without limitation the rights 
## to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies 
## of the Software, and to permit persons to whom the Software is furnished to do so, 
## subject to the following conditions:
##
## The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
##
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
## INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR 
## PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE 
## FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
## ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
##
### Commentary:
##
## tranform a grid between vertical datums
##
## uses htdp and proj cdn
##
### Code:

import os
import sys
import json
from osgeo import gdal
from osgeo import ogr

from cudem import utils
from cudem import regions
from cudem import demfun
from cudem import htdpfun

import cudem.fetches.utils as f_utils

_version = '0.0.2'
_htdp_epsg_desc = lambda x: 'htdp epsg:\n  ' + ' '.join(
    ['\033[1m{}\033[0m\t{}\n'.format(key, x[key]['name']) for key in x]) + '\n'

_usage = """{cmd} ({version}): transform a grid between vertical datums

usage: {cmd} [OPTIONS] input_grid output_grid

  input_grid\t\tThe input raster to transform
  output_grid\t\tThe output transformed raster

 Options:

  -i, --vdatum_in\t\tthe input vertical datum as EPSG code
  -o, --vdatum_out\t\tthe output vertical datum as EPSG code

  --help\tPrint the usage text
  --version\tPrint the version information

 Examples:
 % {cmd} my_dem_navd88.tif my_dem_wgs84_1674.tif --vdatum_in 5703 --vdatum_out 7662

CIRES DEM home page: <http://ciresgroups.colorado.edu/coastalDEM>
""".format(version=_version, cmd=os.path.basename(sys.argv[0]))
#{htdp_epsg}
#, htdp_epsg=_htdp_epsg_desc(htdpfun.HTDP()._reference_frames))

geoids = {'geoid18': 'us_noaa_g2018u0.tif',
          'geoid12b': 'us_noaa_g2012bu0.tif',
          'geoid99': 'us_noaa_g1999u01.tif',
          'geoid03': 'us_noaa_geoid03_conus.tif',
          'geoid09': 'us_noaa_geoid09_conus.tif'}

def list_proj_cdn_epsg(verbose=False):
    _proj_vdatum_index = 'https://cdn.proj.org/files.geojson'
    cdn_index = 'proj_cdn_files.geojson'
    if f_utils.Fetch(_proj_vdatum_index, verbose=verbose).fetch_file(cdn_index) == 0:
        cdn_driver = ogr.GetDriverByName('GeoJSON')
        cdn_ds = cdn_driver.Open(cdn_index, 0)
        cdn_layer = cdn_ds.GetLayer()

        cdn_layer.SetAttributeFilter("type != 'HORIZONTAL_OFFSET'")

        for feat in cdn_layer:
            try:
                print('{} {}-{}'.format(feat.GetField('target_crs_code').split(':')[-1], feat.GetField('target_crs_name'), feat.GetField('name')))
            except:
                pass
                #print(feat.GetField('target_crs_code'))
                #print(feat.GetField('source_crs_code'))
                  
        cdn_ds = None
            
def search_proj_cdn(region, epsg=None, crs_name=None, geoid=None, verbose=True):
    """Search PROJ CDN for transformation grids:
    the PROJ CDN holds transformation grids from around the
    world, including global transformations such as EGM
    """
    
    _proj_vdatum_index = 'https://cdn.proj.org/files.geojson'
    cdn_index = 'proj_cdn_files.geojson'
    if f_utils.Fetch(_proj_vdatum_index, verbose=verbose).fetch_file(cdn_index) == 0:
        cdn_driver = ogr.GetDriverByName('GeoJSON')
        cdn_ds = cdn_driver.Open(cdn_index, 0)
        cdn_layer = cdn_ds.GetLayer()
        _boundsGeom = region.export_as_geom()
        _results = []

        if crs_name is not None:
            cdn_layer.SetAttributeFilter("type != 'HORIZONTAL_OFFSET' AND (target_crs_name LIKE '%{}%' OR source_crs_name LIKE '%{}%')".format(name.upper(), name.upper()))
        elif epsg is not None:
            cdn_layer.SetAttributeFilter("type != 'HORIZONTAL_OFFSET' AND (target_crs_code LIKE '%{}%' OR source_crs_code LIKE '%{}%')".format(epsg, epsg))
        else:
            cdn_layer.SetAttributeFilter("type != 'HORIZONTAL_OFFSET'")

        for feat in cdn_layer:
            if _boundsGeom is not None:
                geom = feat.GetGeometryRef()
                if geom is not None:
                    if _boundsGeom.Intersects(geom):
                        _results.append({})
                        f_j = json.loads(feat.ExportToJson())
                        for key in f_j['properties'].keys():
                            _results[-1][key] = feat.GetField(key)
            else:
                _results.append({})
                f_j = json.loads(feat.ExportToJson())
                for key in f_j['properties'].keys():
                    _results[-1][key] = feat.GetField(key)

        cdn_ds = None
        utils.remove_glob(cdn_index)
        return(_results)

def transform_grid_vertical_htdp(src_grid, dst_grid, epsg_in, epsg_out):
    """Run HTDP to transform src_grd with a vertical datum of epsg_in
    to dst_grd with a vertical datum of epsg_out
    """

    htdp = htdpfun.HTDP()
    utils.echo_msg('{}: HTDP: {}->{}'.format(src_grid, epsg_in, epsg_out))
    
    src_infos = demfun.infos(src_grid)
    src_infos['proj'] = None
    src_region = regions.Region().from_geo_transform(src_infos['geoT'], src_infos['nx'], src_infos['ny'])
    griddef = (src_region.xmin, src_region.ymax, src_region.xmin, src_region.ymin, src_infos['nx'], src_infos['ny'])
    grid = htdp._new_create_grid(griddef)
    
    htdp._write_grid(grid, '_tmp_input.xyz')
    htdp._write_control('_tmp_control.txt', '_tmp_output.xyz', '_tmp_input.xyz', htdp._reference_frames[epsg_in]['htdp_id'], 2012.0, htdp._reference_frames[epsg_out]['htdp_id'], 2012.0)
    htdp.run('_tmp_control.txt')

    out_grid = htdp._read_grid('_tmp_output.xyz', (griddef[5],griddef[4]))
    src_ds = gdal.Open(src_grid)
    src_band = src_ds.GetRasterBand(1)
    src_array = src_band.ReadAsArray()
    final_array = src_array + out_grid

    src_ds = None
    utils.gdal_write(final_array, dst_grid, src_infos)
    utils.remove_glob('_tmp_control.txt', '_tmp_input.xyz', '_tmp_output.xyz')

def transform_grid_vertical_datum(src_grid, dst_grid, vdatum_in, vdatum_out, verbose=True):
    """Transform src_grid with a vertical datum of vdatum_in to dst_grid
    with a vertical datum of vdatum_out
    """
    
    htdp_rf = htdpfun.HTDP()._reference_frames

    utils.echo_msg('{} ({}) -> {} ({})'.format(src_grid, vdatum_in, dst_grid, vdatum_out))
    src_infos = demfun.infos(src_grid)
    #utils.echo_msg(src_infos)
    src_region = regions.Region().from_geo_transform(src_infos['geoT'], src_infos['nx'], src_infos['ny'])
    #utils.echo_msg(src_region)
    
    if vdatum_in in htdp_rf.keys():
        #utils.echo_msg('IN - HTDP')
        ## if output vdatum is also HTDP, do that transformation and done.
        ## this can be sped up by making a lower res transformation grid first
        if vdatum_out in htdp_rf.keys():
            transform_grid_vertical_htdp(src_grid, dst_grid, vdatum_in, vdatum_out)
            return(dst_grid, vdatum_out)
        else:
            cdn_results = search_proj_cdn(src_region, epsg=vdatum_out)
            if len(cdn_results) > 0:
                for _result in cdn_results:
                    ## check results better !! different geoids especially.
                    src_code = int(_result['source_crs_code'].split(':')[-1])
                    dst_code = int(_result['target_crs_code'].split(':')[-1])

                    if src_code in htdp_rf.keys():
                        _trans_grid = _result['name']
                        if f_utils.Fetch(_result['url'], verbose=verbose).fetch_file(_trans_grid) == 0:
                            tmp_infos = demfun.infos(_trans_grid)
                            tmp_region = regions.Region().from_geo_transform(tmp_infos['geoT'], tmp_infos['nx'], tmp_infos['ny'])
                            #rr = regions.regions_reduce(tmp_region, regions.Region().from_list([-180,180,-90,90]))
                            #if not rr.valid_p(check_xy=True):
                            utils.run_cmd('gdalwarp {} {} -s_srs epsg:4326 --config CENTER_LONG 0'.format(_trans_grid, '_{}'.format(_trans_grid)))
                            os.rename('_{}'.format(_trans_grid), _trans_grid)
                                
                            demfun.cut(_trans_grid, src_region, '_{}'.format(_trans_grid))
                            os.rename('_{}'.format(_trans_grid), _trans_grid)
                            
                            _tmp_trans = '_tmp_trans.tif'
                            transform_grid_vertical_htdp(_trans_grid, _tmp_trans, src_code, vdatum_in)
                            utils.run_cmd('gdalwarp {} {} -ts {} {}'.format(_tmp_trans, '_{}'.format(_tmp_trans), src_infos['nx'], src_infos['ny']), verbose=True)
                            os.rename('_{}'.format(_tmp_trans), _tmp_trans)
                            utils.run_cmd('gdal_calc.py -A {} -B {} --calc "A-B" --outfile {}'.format(src_grid, _tmp_trans, dst_grid), verbose=True)
                            utils.remove_glob(_trans_grid, _tmp_trans)
                            return(dst_grid, vdatum_out)
                        break
    else:        
        cdn_results = search_proj_cdn(src_region, epsg=vdatum_in)
        print(cdn_results)
        if len(cdn_results) > 0:
            for _result in cdn_results:
                src_code = int(_result['source_crs_code'].split(':')[-1])
                dst_code = int(_result['target_crs_code'].split(':')[-1])
                if vdatum_in == src_code:
                    if dst_code in htdp_rf.keys():
                        htdp_epsg = dst_code
                        _result_final = _result
                        break
                elif vdatum_in == dst_code:
                    if src_code in htdp_rf.keys():
                        htdp_epsg = src_code
                        _result_final = _result
                        break
                    
            _trans_grid = _result_final['name']
            if f_utils.Fetch(_result_final['url'], verbose=verbose).fetch_file(_trans_grid) == 0:
                tmp_infos = demfun.infos(_trans_grid)
                tmp_region = regions.Region().from_geo_transform(tmp_infos['geoT'], tmp_infos['nx'], tmp_infos['ny'])
                rr = regions.regions_reduce(tmp_region, regions.Region().from_list([-180,180,-90,90]))
                #print(tmp_region)
                #if not rr.valid_p(check_xy=True):
                utils.run_cmd('gdalwarp {} {} -s_srs epsg:4326 --config CENTER_LONG 0'.format(_trans_grid, '_{}'.format(_trans_grid)), verbose=True)
                os.rename('_{}'.format(_trans_grid), _trans_grid)

                demfun.cut(_trans_grid, src_region, '_{}'.format(_trans_grid))
                os.rename('_{}'.format(_trans_grid), _trans_grid)

        if vdatum_out in htdp_rf.keys():
            transform_grid_vertical_htdp(_trans_grid, '_tmp_output.tif', htdp_epsg, vdatum_out)
            utils.run_cmd('gdalwarp {} {} -ts {} {}'.format('_tmp_output.tif', '_tmp_output_sampled.tif', src_infos['nx'], src_infos['ny']), verbose=True)
            utils.run_cmd('gdal_calc.py -A {} -B {} --calc "A+B" --outfile {}'.format(src_grid, '_tmp_output_sampled.tif', dst_grid), verbose=True)
            utils.remove_glob(_trans_grid, '_tmp_output.tif', '_tmp_output_sampled.tif')
            return(dst_grid, vdatum_out)
        else:
            _tmp_dst = '_{}'.format(dst_grid)
            utils.run_cmd('gdalwarp {} {} -ts {} {}'.format(_trans_grid, '_tmp_output_sampled.tif', src_infos['nx'], src_infos['ny']), verbose=True)
            utils.run_cmd('gdal_calc.py -A {} -B {} --calc "A+B" --outfile {}'.format(src_grid, '_tmp_output_sampled.tif', _tmp_dst), verbose=True)
            out_grid, out_vdatum = transform_grid_vertical_datum(_tmp_dst, dst_grid, htdp_epsg, vdatum_out, verbose=True)
            utils.remove_glob(_trans_grid, '_tmp_output_sampled.tif', _tmp_dst)
            return(out_grid, out_vdatum)
        
def main():
    src_grid = None
    dst_grid = None
    vdatum_in = 5703
    vdatum_out = 7662
    verbose = False
    #list_epsg = false

    i = 1

    argv = sys.argv
    while i < len(argv):
        arg = argv[i]

        if arg == '-i' or arg == '--vdatum_in':
            vdatum_in = int(argv[i + 1])
            i = i + 1
        elif arg == '-o' or arg == '--vdatum_out':
            vdatum_out = int(argv[i + 1])
            i = i + 1
        elif arg == '--list-epsg':
            print(_htdp_epsg_desc(htdpfun.HTDP()._reference_frames))
            list_proj_cdn_epsg()
            sys.exit(1)
        elif arg == '--verbose': verbose = True
        elif arg == '-help' or arg == '--help' or arg == '-h':
            print(_usage)
            sys.exit(1)
        elif arg == '-version' or arg == '--version':
            print('vertical_datum_convert.py, version {}'.format(_version))
            sys.exit(1)
        elif src_grid is None:
            src_grid = arg
        elif dst_grid is None:
            dst_grid = arg
        else:
            print(_usage)
            sys.exit(1)
        i = i + 1

    if src_grid is None:
        print(_usage)
        sys.exit(1)

    if dst_grid is None:
        dst_grid = '.'.join(src_grid.split('.')[:-1]) + '_' + str(vdatum_out) + '.' + src_grid.split('.')[-1]

    if not os.path.exists(src_grid):
        print('Error: {} is not a valid file'.format(src_grid))
    else:
        transform_grid_vertical_datum(src_grid, dst_grid, vdatum_in, vdatum_out, verbose=True)
        
if __name__ == '__main__':
    main()

### End
