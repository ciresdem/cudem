### metadata.py
##
## Copyright (c) 2019 - 2021 CIRES Coastal DEM Team
##
## metadata.py is part of CUDEM
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

__version__ = '0.0.8'

def gdal_ogr_mask_union(src_layer, src_field, dst_defn=None):
    '''`union` a `src_layer`'s features based on `src_field` where
    `src_field` holds a value of 0 or 1. optionally, specify
    an output layer defn for the unioned feature.

    returns the output feature class'''
    
    if dst_defn is None:
        dst_defn = src_layer.GetLayerDefn()
        
    multi = ogr.Geometry(ogr.wkbMultiPolygon)
    src_layer.SetAttributeFilter("{} = 1".format(src_field))
    feats = len(src_layer)
    
    _prog = utils.CliProgress('unioning {} features...'.format(feats))
    if feats > 0:
        for n, f in enumerate(src_layer):
            _prog.update_perc((n, feats))
            f_geom = f.geometry()
            #f_geom.CloseRings()
            #try:
            #    f_geom_valid = f_geom.MakeValid()
            #except:
            f_geom_valid = f_geom
                
            #wkt = f_geom_valid.ExportToWkt()
            #wkt_geom = ogr.CreateGeometryFromWkt(wkt)
            #multi.AddGeometryDirectly(wkt_geom)
            multi.AddGeometry(f_geom_valid)
            
    #union = multi.UnionCascaded() ## slow on large multi...
    _prog.end(0, 'unioned {} features'.format(feats))
    utils.echo_msg('setting geometry to unioned feature...')
    out_feat = ogr.Feature(dst_defn)
    #out_feat.SetGeometryDirectly(multi)
    out_feat.SetGeometry(multi)
    union = multi = None
    
    return(out_feat)

def ogr_clip(src_ogr, dst_ogr, clip_region=None, dn="ESRI Shapefile"):
    driver = ogr.GetDriverByName(dn)
    ds = driver.Open(src_ogr, 0)
    layer = ds.GetLayer()

    clip_region.export_as_ogr('tmp_clip.{}'.format(utils.ogr_fext(dn)))
    c_ds = driver.Open('tmp_clip.{}'.format(utils.ogr_fext(dn)), 0)
    c_layer = c_ds.GetLayer()
    
    dst_ds = driver.CreateDataSource(dst_ogr)
    dst_layer = dst_ds.CreateLayer(
        dst_ogr.split('.')[0], geom_type=ogr.wkbMultiPolygon
    )

    layer.Clip(c_layer, dst_layer)
    ds = c_ds = dst_ds = None

def ogr_empty_p(src_ogr, dn='ESRI Shapefile'):
    driver = ogr.GetDriverByName(dn)
    ds = driver.Open(src_ogr, 0)

    if ds is not None:
        layer = ds.GetLayer()
        fc = layer.GetFeatureCount()
        if fc == 0:
            return(True)
        else: return(False)
        
    else: return(True)

class SpatialMetadata:
    """Waffles Spatial Metadata
    Polygonize each datalist entry into vector data source.
    """

    def __init__(
            self,
            data=[],
            src_region=None,
            #inc=None,
            xinc=None,
            yinc=None,
            name='waffles_sm',
            src_srs='epsg:4326',
            dst_srs=None,
            extend=0,
            node='pixel',
            ogr_format='ESRI Shapefile',
            verbose=False
    ):
        """generate spatial-metadata"""
        
        self.data = data
        #self.inc = utils.float_or(inc)
        self.xinc = utils.float_or(xinc)
        self.yinc = utils.float_or(yinc)
        self.src_srs = utils.str_or(src_srs, 'epsg:4326')
        self.dst_srs = utils.str_or(dst_srs, 'epsg:4326')
        self.extend = extend
        self.node = node
        
        self.region = src_region
        self.d_region = self.dist_region()
        
        self.name = name
        self.ogr_format = ogr_format
        self.verbose = verbose

        self._init_data()
        self._init_vector()

    def dist_region(self):
            
        dr = regions.Region().from_region(self.region)
        return(
            dr.buffer(
                x_bv=(self.xinc*self.extend),
                y_bv=(self.yinc*self.extend)
            )
        )

    def _init_data(self):

        self.data = [dlim.DatasetFactory(
            fn=" ".join(['-' if x == "" else x for x in dl.split(",")]),
            src_region=self.d_region,
            verbose=self.verbose,
            src_srs=self.src_srs
        ).acquire() for dl in self.data]

        self.data = [d for d in self.data if d is not None]
        
    def _init_vector(self):
        self.dst_layer = '{}_sm'.format(self.name)
        self.dst_vector = self.dst_layer + '.{}'.format(utils.ogr_fext(self.ogr_format))

        self.v_fields = [
            'Title',
            'Agency',
            'Date',
            'Type',
            'Resolution',
            'HDatum',
            'VDatum',
            'URL'
        ]

        self.t_fields = [
            ogr.OFTString,
            ogr.OFTString,
            ogr.OFTString,
            ogr.OFTString,
            ogr.OFTString,
            ogr.OFTString,
            ogr.OFTString,
            ogr.OFTString
        ]
        utils.remove_glob('{}.*'.format(self.dst_layer))
        utils.gdal_prj_file('{}.prj'.format(self.dst_layer), self.src_srs)
    
        self.ds = ogr.GetDriverByName(self.ogr_format).CreateDataSource(self.dst_vector)
        if self.ds is not None: 
            self.layer = self.ds.CreateLayer(
                '{}'.format(self.dst_layer), None, ogr.wkbMultiPolygon
            )
            [self.layer.CreateField(
                ogr.FieldDefn('{}'.format(f), self.t_fields[i])
            ) for i, f in enumerate(self.v_fields)]
            
            [self.layer.SetFeature(feature) for feature in self.layer]
        else: self.layer = None

    def run(self):
        for xdl in self.data:
            dls = {}
            #xdl.parse_data_lists()
            for e in xdl.parse():
                while e.parent != xdl:
                    e = e.parent
                if e.metadata['name'] not in dls.keys():
                    dls[e.metadata['name']] = {'data': [e] ,'dl': e}

            for x in dls.keys():
                utils.echo_msg('Working on {}'.format(x))
                xdl.data_entries = dls[x]['data']
                p = dls[x]['dl']
                o_v_fields = [
                    p.metadata['title'] if p.metadata['title'] is not None else x,
                    p.metadata['source'],
                    p.metadata['date'],
                    p.metadata['data_type'],
                    p.metadata['resolution'],
                    p.metadata['hdatum'],
                    p.metadata['vdatum'],
                    p.metadata['url']
                ]

                defn = None if self.layer is None else self.layer.GetLayerDefn()
                dl_name = x

                mask_ds, mask_config = xdl.mask_xyz(self.xinc, self.yinc)
                if demfun.gather_infos(mask_ds, scan=True)['zr'][1] == 1:
                    tmp_ds = ogr.GetDriverByName('Memory').CreateDataSource(
                        '{}_poly'.format(dl_name)
                    )
                    
                    if tmp_ds is not None:
                        tmp_layer = tmp_ds.CreateLayer(
                            '{}_poly'.format(dl_name), None, ogr.wkbMultiPolygon
                        )
                        
                        tmp_layer.CreateField(ogr.FieldDefn('DN', ogr.OFTInteger))

                        if self.verbose:
                            utils.echo_msg('polygonizing {} mask...'.format(dl_name))
                            
                        mask_band = mask_ds.GetRasterBand(1)
                        status = gdal.Polygonize(
                            mask_band,
                            None,
                            tmp_layer,
                            0,
                            callback = gdal.TermProgress if self.verbose else None
                        )
                        
                        if len(tmp_layer) > 0:
                            if defn is None:
                                defn = tmp_layer.GetLayerDefn()
                                
                            out_feat = gdal_ogr_mask_union(tmp_layer, 'DN', defn)

                            utils.echo_msg('creating feature {}...'.format(dl_name))
                            for i, f in enumerate(self.v_fields):
                                out_feat.SetField(f, o_v_fields[i])

                            self.layer.CreateFeature(out_feat)

                        if self.verbose:
                            utils.echo_msg('polygonized {}'.format(dl_name))

                            
                    tmp_ds = tmp_layer = out_feat = None
                mask_ds = mask_band = None
                
        utils.echo_msg('Generated SPATIAL METADATA {}'.format(self.name))

## ==============================================
## Command-line Interface (CLI)
## $ spatial_metadata
##
## spatial_metadata cli
## ==============================================
_usage = '''spatial_metadata ({}): generate spatial metadata from a datalist

usage: spatial_metadata [ datalist [ OPTIONS ] ]

  datalist\t\tThe input datalist/entry

 Options:

  -R, --region\t\tSpecifies the desired REGION;
\t\t\tThis can either be a GMT-style region ( -R xmin/xmax/ymin/ymax )
\t\t\tor an OGR-compatible vector file with regional polygons. 
\t\t\tIf a vector file is supplied it will search each region found therein.
\t\t\tIf omitted, use the region gathered from the data in DATALIST.
  -E, --increment\tGridding INCREMENT in native units or GMT-style increments.
\t\t\tWhere INCREMENT is x-inc[/y-inc]
  -O, --output-name\tBASENAME for all outputs.
  -P, --src_srs\t\tPROJECTION of the data and output.
  -X, --extend\t\tNumber of cells with which to EXTEND the REGION.
  -F, --format\t\tOutput OGR format. (ESRI Shapefile)

  -p, --prefix\t\tSet BASENAME (-O) to PREFIX (append inc/region/year info to output BASENAME).
  -r, --grid-node\tUse grid-node registration, default is pixel-node
  -q, --quiet\t\tLower verbosity to a quiet. (overrides --verbose)
  --help\t\tPrint the usage text
  --version\t\tPrint the version information

 Examples:
  % spatial_metadata my_data.datalist

CIRES DEM home page: <http://ciresgroups.colorado.edu/coastalDEM>
'''.format(__version__)

def spat_meta_cli(argv = sys.argv):
    i = 1
    dls = []
    i_regions = []
    these_regions = []
    src_srs = 'epsg:4326'
    #inc = utils.str2inc('1s')
    xinc = utils.str2inc('1s')
    yinc = utils.str2inc('1s')
    node = 'pixel'
    name = 'waffles_spat'
    ogr_format = 'ESRI Shapefile'
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
        elif arg == '-s_srs' or arg == '--s_srs' or arg == '-P':
            src_srs = argv[i + 1]
            i = i + 1
        elif arg == '--increment' or arg == '-E':
            incs = argv[i + 1].split(':')
            xy_inc = incs[0].split('/')
            xinc = utils.str2inc(xy_inc[0])
            if len(xy_inc) > 1:
                yinc = utils.str2inc(xy_inc[1])
            else:
                yinc = utils.str2inc(xy_inc[0])
            i = i + 1
        elif arg[:2] == '-E':
            incs = arg[2:].split(':')
            xy_inc = incs[0].split('/')
            xinc = utils.str2inc(xy_inc[0])
            if len(xy_inc) > 1:
                yinc = utils.str2inc(xy_inc[1])
            else:
                yinc = utils.str2inc(xy_inc[0])
        elif arg == '--extend' or arg == '-X':
            exts = argv[i + 1].split(':')
            extend = utils.int_or(exts[0], 0)
            i += 1
        elif arg[:2] == '-X':
            exts = arg[2:].split(':')
            extend = utils.int_or(exts[0], 0)
        elif arg == '--format' or arg == '-F':
            ogr_format = argv[i + 1]
            i += 1
        elif arg[:2] == '-F':
            ogr_format = argv[2:]
        elif arg == '-r' or arg == '--grid-node': node = 'grid'
        elif arg == '-p' or arg == '--prefix': want_prefix = True
        elif arg == '--quiet' or arg == '-q': want_verbose = False
        elif arg == '-help' or arg == '--help' or arg == '-h':
            sys.stderr.write(_usage)
            sys.exit(1)
        elif arg == '-version' or arg == '--version':
            sys.stdout.write('{}\n'.format(__version__))

            sys.exit(1)
        else: dls.append(arg)
        i = i + 1

    for i_region in i_regions:
        tmp_region = regions.Region().from_string(i_region)
        if tmp_region.valid_p(check_xy=True):
            these_regions.append(tmp_region)
        else:
            i_region_s = i_region.split(':')
            tmp_region = regions.ogr_wkts(i_region_s[0])
            for i in tmp_region:
                if i.valid_p():
                    if len(i_region_s) > 1:
                        these_regions.append(
                            regions.Region().from_string(
                                '/'.join([i.format('str'), i_region_s[1]])
                            )
                        )
                        
                    else:
                        these_regions.append(i)

    if len(these_regions) == 0:
        these_regions = [None]
        utils.echo_error_msg('Could not parse region {}'.format(these_regions))
        sys.stderr.write('{}\n'.format(_usage))
        sys.exit(1)
    else:
        if want_verbose:
            utils.echo_msg(
                'parsed {} region(s)'.format(len(these_regions))
            )

    name_ = name
    for rn, this_region in enumerate(these_regions):
        utils.echo_msg('using region {}'.format(this_region.format('gmt')))
        if len(dls) == 0:
            sys.stderr.write(_usage)
            utils.echo_error_msg('you must specify some type of data')
        else:
            if want_prefix or len(these_regions) > 1:
                name_ = utils.append_fn(name, this_region, inc)

            if os.path.exists('{}_sm.{}'.format(name_, utils.ogr_fext(ogr_format))):
                utils.echo_msg(
                'SPATIAL METADATA {} already exists, skipping...'.format('{}_sm.{}'.format(name_, utils.ogr_fext(ogr_format)))
                    )
            else:
                
                SpatialMetadata(
                    data=dls,
                    src_region=this_region,
                    xinc=xinc,
                    yinc=yinc,
                    extend=extend,
                    src_srs=src_srs,
                    node=node,
                    name=name_,
                    verbose=want_verbose,
                    ogr_format=ogr_format
                ).run()
### End
