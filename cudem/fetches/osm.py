### osm.py
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
import numpy as np
from tqdm import tqdm
from osgeo import ogr
from cudem import utils
from cudem import regions
from cudem import gdalfun
from cudem import xyzfun
from cudem.fetches import fetches
from cudem.fetches import gmrt

## OSM - Open Street Map
## todo: make wrapper modules for 'buildings' and 'coastline' and whaterver else...perhaps
def polygonize_osm_coastline(
        src_ogr, dst_ogr, region=None, include_landmask=True,
        landmask_is_watermask=False, line_buffer=0.0000001, verbose=True
):
    """Polygonize an OSM coastline LineString to the given region

    if include_landmask is True, polygon(s) will be returned for the land areas 
    with a `watermask` valule of 0, otherwise, only the watermask polygon will
    be returned with a `watermask` value of 1.
    """

    def get_ss(point_n, line_geometry, split_geom):
        x_beg = line_geometry.GetX(point_n)
        y_beg = line_geometry.GetY(point_n)
        x_end = line_geometry.GetX(point_n+1)
        y_end = line_geometry.GetY(point_n+1)
        y_ext = y_end
        x_ext = x_beg
        xyz = xyzfun.XYZPoint(x=x_ext, y=y_ext)
        xyz_wkt = xyz.export_as_wkt()
        p_geom = ogr.CreateGeometryFromWkt(xyz_wkt)
        if not p_geom.Within(split_geom):
            y_ext = y_beg
            x_ext = x_end
            xyz = xyzfun.XYZPoint(x=x_ext, y=y_ext)
            xyz_wkt = xyz.export_as_wkt()
            p_geom = ogr.CreateGeometryFromWkt(xyz_wkt)
        else:
            s = (x_end - x_beg)*(y_ext - y_beg) \
                > (y_end - y_beg)*(x_ext - x_beg)
            # ss.append(s)
            return(s)

    
    # Open the input line layer
    line_ds = ogr.Open(src_ogr)
    if line_ds is None:
        return(None)
    
    line_layer = line_ds.GetLayer()
    line_region = regions.Region().from_list(line_layer.GetExtent())
    region_geom = line_region.export_as_geom()
    #region_geom = regions.Region().from_list(line_layer.GetExtent()).export_as_geom()
    ## todo: check if input region is larger than the line region, if so,
    ##       reduce the region to the size of the line region...
    if region is not None and region.valid_p():
        region_geom = region.export_as_geom()
    else:
        region = line_region.copy()

    # Create the output layer
    driver = ogr.GetDriverByName("GPKG")
    if os.path.exists(dst_ogr):
        output_ds = driver.Open(dst_ogr, 1)
        output_layer = output_ds.GetLayer()
    else:
        output_ds = driver.CreateDataSource(dst_ogr)
        output_layer = output_ds.CreateLayer(
            'split_polygons',
            line_layer.GetSpatialRef(),
            ogr.wkbMultiPolygon
        )
        output_layer.CreateField(
            ogr.FieldDefn('watermask', ogr.OFTInteger)
        )
        
    has_feature = False
    with tqdm(
            total=len(line_ds),
            desc='polygonizing osm coastline',
            leave=verbose
    ) as pbar:        
        for line_layer in line_ds:
            utils.echo_msg(line_layer)
            line_type = line_layer.GetGeomType()
            ## feature is a line, polygonize water/land based on which side of
            ## the line each polygon falls...
            if line_type == 2:
                line_geometries = gdalfun.ogr_union_geom(
                    line_layer,
                    ogr.wkbMultiLineString if line_type == 2 else ogr.wkbMultiPolygon,
                    verbose=verbose
                )
                if line_geometries.IsEmpty():
                    continue

                has_feature = True
                poly_line = line_geometries.Buffer(line_buffer)
                split_geoms = region_geom.Difference(poly_line)
                for split_geom in split_geoms:
                    #ss = []
                    #[get_ss(point_n, split_geom) if split_geom.Intersects(line_geometry.Buffer(line_buffer)) else None for point_n in range(0, line_geometry.GetPointCount()) for line_geometry in line_geometries]
                    ss = [get_ss(point_n, line_geometry, split_geom) \
                          for line_geometry in line_geometries \
                          for point_n in range(0, line_geometry.GetPointCount()) \
                          if split_geom.Intersects(line_geometry.Buffer(line_buffer))]
                    # for line_geometry in line_geometries:
                    #     if split_geom.Intersects(line_geometry.Buffer(line_buffer)):
                    #         point_count = line_geometry.GetPointCount()
                    #         utils.echo_msg(point_count)
                    #         for point_n in range(0, point_count-1):
                    #             ss.append(get_ss(point_n, split_geom))
                    #             # x_beg = line_geometry.GetX(point_n)
                    #             # y_beg = line_geometry.GetY(point_n)
                    #             # x_end = line_geometry.GetX(point_n+1)
                    #             # y_end = line_geometry.GetY(point_n+1)
                    #             # y_ext = y_end
                    #             # x_ext = x_beg
                    #             # xyz = xyzfun.XYZPoint(x=x_ext, y=y_ext)
                    #             # xyz_wkt = xyz.export_as_wkt()
                    #             # p_geom = ogr.CreateGeometryFromWkt(xyz_wkt)
                    #             # if not p_geom.Within(split_geom):
                    #             #     y_ext = y_beg
                    #             #     x_ext = x_end
                    #             #     xyz = xyzfun.XYZPoint(x=x_ext, y=y_ext)
                    #             #     xyz_wkt = xyz.export_as_wkt()
                    #             #     p_geom = ogr.CreateGeometryFromWkt(xyz_wkt)

                    #             # if p_geom.Within(split_geom):
                    #             #     s = (x_end - x_beg)*(y_ext - y_beg) \
                    #             #         > (y_end - y_beg)*(x_ext - x_beg)
                    #             #     ss.append(s)
                    #ss = ss[ss is not None]
                    utils.echo_msg(ss)
                    if all(ss):
                        s = True
                    elif not any(ss):
                        s = False
                    else:
                        if np.count_nonzero(ss) > len(ss) / 2:
                            s = True
                        else:
                            s = False

                    out_feature = ogr.Feature(output_layer.GetLayerDefn())
                    out_feature.SetGeometry(split_geom)
                    if landmask_is_watermask:
                        s = False if s else True

                    if s == 0:
                        out_feature.SetField('watermask', 1)
                        output_layer.CreateFeature(out_feature)

                    if include_landmask:
                        if s == 1:
                            out_feature.SetField('watermask', 0)
                            output_layer.CreateFeature(out_feature)                            
            
            ## feature is a polygon, which in osm means an island.
            if line_type == 6:
                for line_feature in line_layer:
                    line_geometry = line_feature.geometry()
                    line_geometry = ogr.ForceTo(line_geometry, ogr.wkbLinearRing)
                    if line_geometry.IsEmpty():
                        continue

                    has_feature = 1
                    for feature in output_layer:
                        feature_geom = feature.geometry()
                        if feature_geom.Intersects(line_geometry):
                            feature_geoms = feature_geom.Difference(line_geometry)
                            feature.SetGeometry(feature_geoms)
                            output_layer.SetFeature(feature)

                    out_feature = ogr.Feature(output_layer.GetLayerDefn())
                    out_feature.SetGeometry(line_geometry)
                    s = True
                    if landmask_is_watermask:
                        s = False if s else True

                    if s == 0:
                        out_feature.SetField('watermask', 1)
                        output_layer.CreateFeature(out_feature)

                    if include_landmask:
                        if s == 1:
                            out_feature.SetField('watermask', 0)
                            output_layer.CreateFeature(out_feature)

            pbar.update()
                            
    ## no features in the input osm coastline, so the entire region is either land or water.
    ## find the center point of the region and check the z value from gmrt.
    if not has_feature:
        center_pnt = region.center()
        if center_pnt is not None:
            center_z = utils.int_or(
                gmrt.gmrt_fetch_point(latitude=center_pnt[1],
                                      longitude=center_pnt[0])
            )
            out_feature = ogr.Feature(output_layer.GetLayerDefn())
            out_feature.SetGeometry(region_geom)
            if center_z >= 0:
                s = True
            else:
                s = False

            if landmask_is_watermask:
                s = False if s else True

            if s == 0:
                out_feature.SetField('watermask', 1)
                output_layer.CreateFeature(out_feature)

            if include_landmask:
                if s == 1:
                    out_feature.SetField('watermask', 0)
                    output_layer.CreateFeature(out_feature)
                    
    line_ds = output_ds = None
    return(dst_ogr)


class OpenStreetMap(fetches.FetchModule):
    """OpenStreetMap data.
    
    OpenStreetMap is a free, editable map of the whole world that is 
    being built by volunteers largely from scratch and released with an 
    open-content license.
    
    https://wiki.openstreetmap.org/

    coastline: https://wiki.openstreetmap.org/wiki/Tag:natural=coastline

    `q` should be the query to send to the osm interpreter, such as:
    q = '''
    (node;
    <;
    >;
    );
    out meta;
    '''
    or `q` can be a keyword to fetch specific things, currently we accept:
    'buildings', 'coastline'

    < osm:q=None:fmt=osm:planet=False:chunks=True:min_length=None >
    """
    
    def __init__(self, q=None, h='', fmt='osm', planet=False, chunks=False,
                 min_length = None, **kwargs):
        super().__init__(name='osm', **kwargs)
        self.q = q
        self.fmt = fmt
        self.planet = planet
        self.chunks = chunks
        self.h = h
        if self.q == 'buildings':
            #self.h = '[maxsize:2000000000]'
            self.h = '[timeout:3600]'
            self.q = '''
            (way["building"]{};
            relation["building"]["type"="multipolygon"];
            );
            (._;>;);
            out meta;
            '''.format(
                '(if: length() > {})'.format(
                    min_length
                ) if min_length is not None else ''
            )
            
        if self.q == 'coastline':
            #self.h = '[maxsize:2000000000]'
            self.h = '[timeout:32]'
            self.q = '''
            (way["natural"="coastline"]{};
            relation["type"="lines"];
            );
            (._;>;);
            out meta;
            '''.format(
                '(if: length() > {})'.format(
                    min_length
                ) if min_length is not None else ''
            )

        if self.q == 'rivers':
            #self.h = '[maxsize:2000000000]'
            self.h = '[timeout:3600]'
            self.q = '''
            //(way["natural"="water"]{};
            //relation["type"="lines"];
            nwr["water"="river"];
            //);
            (._;>;);
            out meta;
            '''.format(
                '(if: length() > {})'.format(
                    min_length
                ) if min_length is not None else ''
            )

        if self.q == 'place2':
            #self.h = '[maxsize:2000000000]'
            self.h = '[timeout:3600]'
            self.q = '''
            nwr["place"="city"];
            nwr["name"~"Paris",1]["type"="multipolygon"];
            (._;>;);
            out meta;
            '''.format()

        if self.q == 'place':
            #self.h = '[maxsize:2000000000]'
            self.h = '[timeout:3600]'
            self.q = '''
            area[boundary=administrative][name="Paris"];
            wr[place~"^(sub|town|city|count|state)"](area);
            (._;>;);
            out meta;
            '''.format()
            
            
        ## various OSM URLs
        self._osm_api = 'https://lz4.overpass-api.de/api/interpreter'
        self._osm_api2 = 'https://overpass.kumi.systems/api/interpreter'
        self._osm_api3 = 'https://overpass.openstreetmap.fr/api/interpreter'
        self._osm_planet_bz2 = ('https://ftpmirror.your.org/pub/openstreetmap/'
                                'planet/planet-latest.osm.bz2')
        self._osm_planet = ('https://ftpmirror.your.org/pub/openstreetmap/'
                            'pbf/planet-latest.osm.pbf')
        self._osm_continents = 'https://download.geofabrik.de/'

        ## Set user-agent and referer
        self.headers = { 'User-Agent': ('Mozilla/5.0 (Windows NT 10.0; '
                                        'Win64; x64; rv:89.0) Gecko/20100101 '
                                        'Firefox/89.0'),
                         'referer': 'https://lz4.overpass-api.de/' }

        self.check_size=False

        
    def run(self):
        """Run the OSM fetches module"""
        
        if self.region is None:
            return([])

        ## fetch whole planet
        if self.planet:
            self.add_entry_to_results(
                self._osm_planet,
                os.path.join(self._outdir, 'planet-latest.osm.pbf'),
                'pbf'
            )

        ## fetch in chunks
        elif self.chunks:
            x_delta = self.region.xmax - self.region.xmin
            y_delta = self.region.ymax - self.region.ymin
            incs = self.region.increments(1000,1000)

            ## break up the requests into .05 degree chunks for
            ## better usage of the OSM API
            if x_delta > .25 or y_delta > .25:
                xcount, ycount, gt = self.region.geo_transform(x_inc=incs[0], y_inc=incs[1])
                if x_delta >= y_delta:
                    n_chunk = int(xcount*(.1/x_delta))
                elif y_delta > x_delta:
                    n_chunk = int(ycount*(.1/y_delta))
            else:
                n_chunk = None

            these_regions = self.region.chunk(incs[0], n_chunk=n_chunk)
            utils.echo_msg(
                'chunking OSM request into {} regions'.format(
                    len(these_regions)
                )
            )
            for this_region in these_regions:
                c_bbox = this_region.format('osm_bbox')
                out_fn = 'osm_{}'.format(this_region.format('fn_full'))
                osm_q_bbox  = '''
                {1}{2}[bbox:{0}];'''.format(
                    c_bbox,
                    '[out:{}]'.format(self.fmt) if self.fmt != 'osm' else '',
                    self.h
                )
                osm_q = '''
                (node;
                <;
                >;
                );
                out meta;
                '''                
                osm_q_ = osm_q_bbox + (osm_q if self.q is None else self.q)
                osm_data = fetches.urlencode({'data': osm_q_})
                osm_data_url = self._osm_api + '?' + osm_data
                self.add_entry_to_results(
                    osm_data_url,
                    '{}.{}'.format(out_fn, self.fmt),
                    'osm'
                )
        else:
            c_bbox = self.region.format('osm_bbox')
            out_fn = 'osm_{}'.format(self.region.format('fn_full'))
            osm_q_bbox  = '''
            {1}[bbox:{0}];'''.format(
                c_bbox, '[out:{}]'.format(
                    self.fmt
                ) if self.fmt != 'osm' else ''
            )

            osm_q = '''
            (node;
            <;
            >;
            );
            out meta;
            '''

            osm_q_ = osm_q_bbox + (osm_q if self.q is None else self.q)
            osm_data = fetches.urlencode({'data': osm_q_})
            osm_data_url = self._osm_api + '?' + osm_data            
            self.add_entry_to_results(
                osm_data_url,
                '{}.{}'.format(out_fn, self.fmt),
                'osm'
            )

### End
