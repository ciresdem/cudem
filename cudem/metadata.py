### metadata.py
##
## Copyright (c) 2019 - 2021 CIRES Coastal DEM Team
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
## generate spatial metadata from a datalist
##
### Code:

import os
import sys
from osgeo import ogr
from osgeo import gdal
from cudem import utils
from cudem import dlim
from cudem import regions
from cudem import demfun

_version = '0.0.2'

def gdal_ogr_mask_union(src_layer, src_field, dst_defn = None):
    '''`union` a `src_layer`'s features based on `src_field` where
    `src_field` holds a value of 0 or 1. optionally, specify
    an output layer defn for the unioned feature.

    returns the output feature class'''
    
    if dst_defn is None: dst_defn = src_layer.GetLayerDefn()
    multi = ogr.Geometry(ogr.wkbMultiPolygon)
    feats = len(src_layer)
    utils.echo_msg('unioning {} features'.format(feats))
    for n, f in enumerate(src_layer):
        gdal.TermProgress_nocb((n+1 / feats) * 100)
        if f.GetField(src_field) == 0:
            src_layer.DeleteFeature(f.GetFID())
        elif f.GetField(src_field) == 1:
            f.geometry().CloseRings()
            wkt = f.geometry().ExportToWkt()
            multi.AddGeometryDirectly(ogr.CreateGeometryFromWkt(wkt))
            src_layer.DeleteFeature(f.GetFID())
    #union = multi.UnionCascaded() ## slow on large multi...
    out_feat = ogr.Feature(dst_defn)
    out_feat.SetGeometryDirectly(multi)
    #union = multi = None
    return(out_feat)

def ogr_clip(src_ogr, dst_ogr, clip_region = None, dn = "ESRI Shapefile"):
    driver = ogr.GetDriverByName(dn)
    ds = driver.Open(src_ogr, 0)
    layer = ds.GetLayer()

    clip_region.export_as_ogr('tmp_clip.shp')
    c_ds = driver.Open('tmp_clip.shp', 0)
    c_layer = c_ds.GetLayer()
    
    dst_ds = driver.CreateDataSource(dst_ogr)
    dst_layer = dst_ds.CreateLayer(dst_ogr.split('.')[0], geom_type=ogr.wkbMultiPolygon)

    layer.Clip(c_layer, dst_layer)

    ds = c_ds = dst_ds = None

def ogr_empty_p(src_ogr):
    driver = ogr.GetDriverByName('ESRI Shapefile')
    ds = driver.Open(src_ogr, 0)

    if ds is not None:
        layer = ds.GetLayer()
        fc = layer.GetFeatureCount()
        if fc == 0:
            return(True)
        else: return(False)
    else: return(True)

## ==============================================
## Waffles Spatial Metadata
## Polygonize each datalist entry into vector data source
## ==============================================
class SpatialMetadata: #(waffles.Waffle):

    def __init__(self, data=[], src_region=None, inc=None, name='waffles_sm', epsg=4326,
                 warp=None, extend=0, node='pixel', verbose=False):
        """generate spatial-metadata

        Args:
          geojson(bool): generate a geojson output
        """
        
        #super().__init__(**kwargs)
        self.data = data
        self.inc = utils.float_or(inc)
        self.epsg = utils.int_or(epsg)
        self.warp = utils.int_or(warp)
        self.extend = extend
        self.node = node
        
        self.region = src_region
        self.d_region = self.dist_region()
        
        self.name = name
        self.verbose = verbose

        self._init_data()
        self._init_vector()

    def dist_region(self):
            
        dr = regions.Region().from_region(self.region)
        return(dr.buffer((self.inc*self.extend)))

    def _init_data(self):

        self.data = [dlim.DatasetFactory(
            fn=" ".join(['-' if x == "" else x for x in dl.split(":")]),
            src_region=self.d_region, verbose=self.verbose,
            epsg=self.epsg).acquire_dataset() for dl in self.data]

        self.data = [d for d in self.data if d is not None]
        
        for d in self.data:
            d.parse()        
    
    def _init_vector(self):
        self.dst_layer = '{}_sm'.format(self.name)
        self.dst_vector = self.dst_layer + '.shp'
        self.v_fields = ['Name', 'Agency', 'Date', 'Type', 'Resolution', 'HDatum', 'VDatum', 'URL']
        self.t_fields = [ogr.OFTString, ogr.OFTString, ogr.OFTString, ogr.OFTString,
                         ogr.OFTString, ogr.OFTString, ogr.OFTString, ogr.OFTString]
        utils.remove_glob('{}.*'.format(self.dst_layer))
        utils.gdal_prj_file('{}.prj'.format(self.dst_layer), self.epsg)
    
        self.ds = ogr.GetDriverByName('ESRI Shapefile').CreateDataSource(self.dst_vector)
        if self.ds is not None: 
            self.layer = self.ds.CreateLayer('{}'.format(self.dst_layer), None, ogr.wkbMultiPolygon)
            [self.layer.CreateField(ogr.FieldDefn('{}'.format(f), self.t_fields[i])) for i, f in enumerate(self.v_fields)]
            [self.layer.SetFeature(feature) for feature in self.layer]
        else: self.layer = None

    def run(self):
        for xdl in self.data:
            for x in xdl.data_lists.keys():
                xdl.data_entries = xdl.data_lists[x]
                dl_name = x
                o_v_fields = [dl_name, 'Unknown', '0', 'xyz_elevation', 'Unknown', 'WGS84', 'NAVD88', 'URL']
                defn = None if self.layer is None else self.layer.GetLayerDefn()
                [x for x in xdl.mask_xyz('{}.tif'.format(dl_name), self.inc)]

                if demfun.infos('{}.tif'.format(dl_name), scan=True)['zr'][1] == 1:
                    tmp_ds = ogr.GetDriverByName('ESRI Shapefile').CreateDataSource('{}_poly.shp'.format(dl_name))
                    if tmp_ds is not None:
                        tmp_layer = tmp_ds.CreateLayer('{}_poly'.format(dl_name), None, ogr.wkbMultiPolygon)
                        tmp_layer.CreateField(ogr.FieldDefn('DN', ogr.OFTInteger))
                        demfun.polygonize('{}.tif'.format(dl_name), tmp_layer, verbose=self.verbose)

                        if len(tmp_layer) > 1:
                            if defn is None: defn = tmp_layer.GetLayerDefn()
                            out_feat = gdal_ogr_mask_union(tmp_layer, 'DN', defn)
                            [out_feat.SetField(f, o_v_fields[i]) for i, f in enumerate(self.v_fields)]
                            self.layer.CreateFeature(out_feat)
                    tmp_ds = None
                    utils.remove_glob('{}_poly.*'.format(dl_name), 'tmp.tif')
        self.ds = None

        utils.run_cmd('ogrinfo -dialect SQLITE -sql "UPDATE {} SET geometry = ST_MakeValid(geometry)" {}\
        '.format(self.dst_layer, self.dst_vector), verbose=True)

    def yield_xyz(self):
        for xdl in self.data:
            for x in xdl.data_lists.keys():
                xdl.data_entries = xdl.data_lists[x]
                dl_name = x
                o_v_fields = [dl_name, 'Unknown', '0', 'xyz_elevation', 'Unknown', 'WGS84', 'NAVD88', 'URL']
                defn = None if self.layer is None else self.layer.GetLayerDefn()
                for xyz in xdl.mask_xyz('{}.tif'.format(dl_name), self.inc):
                    yield(xyz)

                if demfun.infos('{}.tif'.format(dl_name), scan=True)['zr'][1] == 1:
                    tmp_ds = ogr.GetDriverByName('ESRI Shapefile').CreateDataSource('{}_poly.shp'.format(dl_name))
                    if tmp_ds is not None:
                        tmp_layer = tmp_ds.CreateLayer('{}_poly'.format(dl_name), None, ogr.wkbMultiPolygon)
                        tmp_layer.CreateField(ogr.FieldDefn('DN', ogr.OFTInteger))
                        demfun.polygonize('{}.tif'.format(dl_name), tmp_layer, verbose=self.verbose)

                        if len(tmp_layer) > 1:
                            if defn is None: defn = tmp_layer.GetLayerDefn()
                            out_feat = gdal_ogr_mask_union(tmp_layer, 'DN', defn)
                            [out_feat.SetField(f, o_v_fields[i]) for i, f in enumerate(self.v_fields)]
                            self.layer.CreateFeature(out_feat)
                    tmp_ds = None
                    utils.remove_glob('{}_poly.*'.format(dl_name), 'tmp.tif')
        self.ds = None

        utils.run_cmd('ogrinfo -spat {} -dialect SQLITE -sql "UPDATE {} SET geometry = ST_MakeValid(geometry)" {}\
        '.format(dr.format('ul_lr'), self.dst_layer, self.dst_vector))
        
_usage = '''spatial_metadata.py ({}): generate spatial metadata from a datalist

usage: spatial_metadata.py [ datalist [ OPTIONS ] ]

  datalist\t\tThe input datalist/entry

 Options:

  -R, --region\t\tSpecifies the desired REGION;
\t\t\tThis can either be a GMT-style region ( -R xmin/xmax/ymin/ymax )
\t\t\tor an OGR-compatible vector file with regional polygons. 
\t\t\tIf a vector file is supplied it will search each region found therein.
\t\t\tIf omitted, use the region gathered from the data in DATALIST.
  -E, --increment\tGridding CELL-SIZE in native units or GMT-style increments.
\t\t\tappend :<inc> to resample the output to the given <inc>: -E.3333333s:.1111111s
  -O, --output-name\tBASENAME for all outputs.
  -P, --epsg\t\tHorizontal projection of data as EPSG code [4326]
  -X, --extend\t\tNumber of cells with which to EXTEND the REGION.
\t\t\tappend :<num> to extend the processing region: -X6:12

  -p, --prefix\t\tSet BASENAME (-O) to PREFIX (append inc/region/year info to output BASENAME).
  -r, --grid-node\tUse grid-node registration, default is pixel-node
  -q, --quiet\t\tLower verbosity to a quiet. (overrides --verbose)
  --help\tPrint the usage text
  --version\tPrint the version information

 Examples:
  % spatial_metadata.py my_data.datalist

CIRES DEM home page: <http://ciresgroups.colorado.edu/coastalDEM>
'''.format(_version)


def spat_meta_cli(argv = sys.argv):
    i = 1
    dls = []
    i_regions = []
    these_regions = []
    epsg = 4326
    inc = utils.str2inc('1s')
    node = 'pixel'
    name = 'waffles_spat'
    extend = 0
    want_verbose = True
    want_prefix = False

    argv = sys.argv
    while i < len(argv):
        arg = sys.argv[i]
        if arg == '--region' or arg == '-R':
            i_regions.append(str(argv[i + 1]))
            i = i + 1
        elif arg[:2] == '-R':
            i_regions.append(str(arg[2:]))
        elif arg == '--outname' or arg == '-O':
            name = argv[i + 1]
            i += 1
        elif arg[:2] == '-O': name = arg[2:]
        elif arg == '-s_epsg' or arg == '--s_epsg' or arg == '-P':
            epsg = argv[i + 1]
            i = i + 1
        elif arg == '--increment' or arg == '-E':
            inc = utils.str2inc(argv[i + 1])
            i = i + 1
        elif arg[:2] == '-E':
            inc = utils.str2inc(arg[2:])
        elif arg == '--extend' or arg == '-X':
            extend = utils.int_or(argv[i + 1], 0)
            i = i + 1
        elif arg[:2] == '-X':
            extend = utils.int_or(arg[2:], 0)
        elif arg == '-r' or arg == '--grid-node': node = 'grid'
        elif arg == '-p' or arg == '--prefix': want_prefix = True
        elif arg == '--quiet' or arg == '-q': want_verbose = False
        elif arg == '-help' or arg == '--help' or arg == '-h':
            sys.stderr.write(_usage)
            sys.exit(1)
        elif arg == '-version' or arg == '--version':
            sys.stderr.write('{}\n'.format(_version))
            sys.exit(1)
        else: dls.append(arg)
        i = i + 1

    for i_region in i_regions:
        tmp_region = regions.Region().from_string(i_region)
        if tmp_region.valid_p():
            these_regions.append(tmp_region)
        else:
            tmp_region = regions.ogr_wkts(i_region)
            for i in tmp_region:
                if i.valid_p():
                    these_regions.append(i)
                    
    if len(these_regions) == 0:
        these_regions = [None]
    else:
        if want_verbose: utils.echo_msg('parsed {} region(s)'.format(len(these_regions)))

    for rn, this_region in enumerate(these_regions):
        if len(dls) == 0:
            sys.stderr.write(_usage)
            utils.echo_error_msg('you must specify some type of data')
        else:
            if want_prefix or len(these_regions) > 1:
                name = utils.append_fn(name, this_region, inc)
            #[x for x in waffles.Waffle(data=dls, src_region=this_region, inc=inc, extend=extend, epsg=epsg, node=node, name=name, verbose=want_verbose).spat_meta(yxyz=False)]
            SpatialMetadata(data=dls, src_region=this_region, inc=inc, extend=extend, epsg=epsg,
                            node=node, name=name, verbose=want_verbose).run()
        
### End
