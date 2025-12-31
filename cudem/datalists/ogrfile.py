### ogrfile.py 
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
##
## ogrfile.py is part of CUDEM
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
### Examples:
##
### TODO:
##
### Code:

import json
import numpy as np
from osgeo import ogr

from cudem import utils
from cudem.datalists.dlim import ElevationDataset

class OGRFile(ElevationDataset):
    """providing an OGR 3D point dataset parser.

    Useful for data such as S-57, ENC, E-Hydro, Etc.

    -----------
    Parameters:

    ogr_layer (str/int): the OGR layer containing elevation data
    elev_field (str): the field containing the z values
    weight_field (str): the field containing weight values
    uncertainty_field (str): the field containing uncertainty_values
    z_scale (float): scale the output z values    
    """

    _known_layer_names = ['SOUNDG', 'SurveyPoint_HD', 'SurveyPoint']
    _known_elev_fields = ['Elevation', 'elev', 'z', 'height', 'depth',
                          'topography', 'surveyPointElev', 'Z_depth',
                          'Z_height']
    
    def __init__(self,
                 ogr_layer=None,
                 elev_field=None,
                 weight_field=None,
                 uncertainty_field=None,
                 z_scale=None,
                 elevation_value=None,
                 **kwargs):
        super().__init__(**kwargs)
        self.ogr_layer = ogr_layer
        self.elev_field = elev_field
        self.weight_field = weight_field
        self.uncertainty_field = uncertainty_field
        self.z_scale = utils.float_or(z_scale)
        self.elevation_value = utils.float_or(elevation_value)

        
    def find_elevation_layer(self, ds_ogr):
        for l in self._known_layer_names:
            test_layer = ds_ogr.GetLayerByName(l)
            if test_layer is not None:
                return((l, test_layer))
        return(None, None)

    
    def find_elevation_field(self, ogr_feature):
        for f in self._known_elev_fields:            
            test_field = ogr_feature.GetField(f)
            if test_field is not None:
                return((l, test_field))
            
        return(None, None)


    def generate_inf(self, make_grid=True, make_block_mean=False, block_inc=None):
        """Generate metadata (INF) for the OGR file."""
        
        # 1. Open Datasource
        ds_ogr = ogr.Open(self.fn)
        if ds_ogr is None: 
            return self.infos

        # 2. Find the correct Layer
        layer = None
        if self.ogr_layer is not None:
            layer = ds_ogr.GetLayer(self.ogr_layer)
        
        if layer is None:
            # Auto-detect layer using your existing logic logic
            layer_name, layer = self.find_elevation_layer(ds_ogr)
            
        if layer is None:
            # Fallback to first layer
            layer = ds_ogr.GetLayer(0)

        if layer is None:
            return self.infos

        # 3. Extract SRS from Layer
        if self.src_srs is None:
            spatial_ref = layer.GetSpatialRef()
            if spatial_ref:
                # Export to Proj4 or WKT for CUDEM
                self.infos.src_srs = spatial_ref.ExportToProj4()
                # Also set the instance variable so it persists
                self.src_srs = self.infos.src_srs
        else:
            self.infos.src_srs = self.src_srs

        # 4. Optional: Get Extent/Count (Fast Metadata)
        try:
            extent = layer.GetExtent()
            # OGR Extent is (minx, maxx, miny, maxy)
            self.infos.minmax = [extent[0], extent[1], extent[2], extent[3], None, None]
            self.infos.numpts = layer.GetFeatureCount()
            
            # Create WKT for the region
            r = regions.Region().from_list(self.infos.minmax)
            self.infos.wkt = r.export_as_wkt()
        except Exception:
            pass

        # 5. Handle Grids (if requested via flags)
        if make_grid or make_block_mean:
            return super().generate_inf(make_grid, make_block_mean, block_inc)

        ds_ogr = None
        return self.infos
    
    
    def yield_points(self):
        ds_ogr = ogr.Open(self.fn)
        count = 0
        utils.echo_debug_msg(f'ogrfile layer name is: {self.ogr_layer}')
        if ds_ogr is not None:
            layer_name = None
            if self.ogr_layer is None:
                layer_name, layer_s = self.find_elevation_layer(ds_ogr)

                if layer_name is None:
                    layer_s = ds_ogr.GetLayer()
                    
            elif utils.str_or(self.ogr_layer) is not None:
                layer_name = self.ogr_layer
                layer_s = ds_ogr.GetLayer(str(self.ogr_layer))
                if layer_s is None:
                    layer_name, layer_s = self.find_elevation_layer(ds_ogr)
                    
            elif utils.int_or(self.ogr_layer) is not None:
                layer_s = ds_ogr.GetLayer(self.ogr_layer)
            else:
                layer_s = ds_ogr.GetLayer()

            if layer_s is not None:
                field_names = [field.name for field in layer_s.schema]
                if self.elev_field not in field_names:                    
                    for field_name in field_names:
                        if field_name in self._known_elev_fields:
                            self.elev_field = field_name
                            break
                    
                if self.region is not None:
                    layer_s.SetSpatialFilter(
                        self.region.export_as_geom() \
                        if self.transform['transformer'] is None \
                        else self.transform['trans_region'].export_as_geom()
                    )

                for f in layer_s:
                    geom = f.GetGeometryRef()
                    g = json.loads(geom.ExportToJson())
                    #xyzs = g['coordinates']
                    xyzs = []
                    if g['type'] == 'Point':
                        xyzs = [g['coordinates']]
                    else:
                        for i in g['coordinates']:
                            if isinstance(i[0], list):
                                xyzs += i
                            else:
                                xyzs.append(i)

                    if self.uncertainty_field is not None:
                        unc = utils.float_or(f.GetField(self.uncertainty_field))
                    else:
                        unc = None

                    if self.weight_field is not None:
                        weight = utils.float_or(f.GetField(self.weight_field))
                    else:
                        weight = None
                                
                    for xyz in xyzs:                            
                        if not geom.Is3D():
                            # if self.elev_field is None:
                            #     self.elev_field = self.find_elevation_field(f)

                            if self.elev_field is None:
                                if self.elevation_value is None:
                                    elev = 0
                                else:
                                    elev = self.elevation_value
                            else:
                                elev = utils.float_or(f.GetField(self.elev_field))

                            if elev is not None:
                                xyz.append(elev)
                            else:
                                continue

                        else:
                            if self.elev_field is None:
                                if self.elevation_value is None:
                                    elev = 0
                                else:
                                    elev = self.elevation_value
                                    
                            else:
                                elev = utils.float_or(f.GetField(self.elev_field))

                            if elev is not None:
                                xyz.append(elev)
                            else:
                                continue
                                
                            if isinstance(xyz[0], list):
                                for x in xyz:
                                    #for xx in x:
                                    if utils.float_or(x) is not None:
                                        x[2] = elev

                    if isinstance(xyzs[0], list):
                        #[x.append(weight if weight is not None else 1) for x in xyzs]
                        #[x.append(unc if unc is not None else 0) for x in xyzs]
                        
                        for x in xyzs:
                            x.append(weight if weight is not None else 1.)
                            x.append(unc if unc is not None else 0.)
                            
                        points = np.rec.fromrecords(
                            xyzs, names='x, y, z, w, u'
                        )
                        #points = points[points['z'] is not None]
                        if self.z_scale is not None:
                            points['z'] *= self.z_scale

                        yield(points)

                        # for x in xyzs:
                        #     x.append(weight if weight is not None else 1)
                        #     x.append(unc if unc is not None else 0)
                        #     points = np.rec.fromrecords(
                        #         [x], names='x, y, z, w, u'
                        #     )
                        #     if self.z_scale is not None:
                        #         points['z'] *= self.z_scale

                        #     yield(points)
                                
            ds_ogr = layer_s = None


### End

# ### ogrfile.py 
# ##
# ## Copyright (c) 2010 - 2025 Regents of the University of Colorado
# ##
# ## ogrfile.py is part of CUDEM
# ##
# ## Permission is hereby granted, free of charge, to any person obtaining a copy 
# ## of this software and associated documentation files (the "Software"), to deal 
# ## in the Software without restriction, including without limitation the rights 
# ## to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies 
# ## of the Software, and to permit persons to whom the Software is furnished to do so, 
# ## subject to the following conditions:
# ##
# ## The above copyright notice and this permission notice shall be included in all
# ## copies or substantial portions of the Software.
# ##
# ## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
# ## INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR 
# ## PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE 
# ## FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
# ## ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# ## SOFTWARE.
# ##
# ### Commentary:
# ##
# ## OGR Vector Data Parser (S-57, Shapefiles, GeoJSON, etc.)
# ##
# ### Code:

# import os
# import sys
# import numpy as np
# from osgeo import ogr, osr

# from cudem import utils
# from cudem import regions
# from cudem import srsfun 
# from cudem.datalists.dlim import ElevationDataset

# class OGRFile(ElevationDataset):
#     """Providing an OGR 3D point dataset parser.
#     Useful for data such as S-57, ENC, E-Hydro, Shapefiles, etc.
#     """

#     _known_layer_names = ['SOUNDG', 'SurveyPoint_HD', 'SurveyPoint', 'Mass_Point', 'Spot_Elevation']
#     _known_elev_fields = ['Elevation', 'elev', 'z', 'height', 'depth', 'val', 'value',
#                           'topography', 'surveyPointElev', 'Z_depth', 'Z_height']
    
#     def __init__(self,
#                  ogr_layer=None,
#                  elev_field=None,
#                  weight_field=None,
#                  uncertainty_field=None,
#                  z_scale=None,
#                  elevation_value=None,
#                  **kwargs):
        
#         super().__init__(**kwargs)
        
#         self.ogr_layer = ogr_layer
#         self.elev_field = elev_field
#         self.weight_field = weight_field
#         self.uncertainty_field = uncertainty_field
        
#         self.z_scale = utils.float_or(z_scale)
#         self.elevation_value = utils.float_or(elevation_value)

        
#     def generate_inf(self, make_grid=True, make_block_mean=False, block_inc=None):
#         """Generate metadata (INF) for the OGR file.
        
#         Attempts to use OGR's fast GetExtent/GetFeatureCount first.
#         If grids are requested, falls back to full scan via parent.
#         """
        
#         ## If grids are requested, we MUST do a full scan anyway.
#         if make_grid or make_block_mean:
#             return super().generate_inf(make_grid, make_block_mean, block_inc)

#         ## Fast Metadata Extraction
#         ds_ogr = ogr.Open(self.fn)
#         if ds_ogr is None: return self.infos

#         layer = self._get_layer(ds_ogr)
#         if layer is None: return self.infos

#         ## Extent
#         try:
#             extent = layer.GetExtent() # (minx, maxx, miny, maxy)
#             ## OGR Extent is X/Y only. We don't know Z min/max without scanning.
#             ## We set a placeholder or skip Z bounds here.
#             ## (minx, maxx, miny, maxy, minz, maxz)
#             self.infos.minmax = [extent[0], extent[1], extent[2], extent[3], None, None]
            
#             ## Create Region for WKT
#             r = regions.Region().from_list(self.infos.minmax)
#             self.infos.wkt = r.export_as_wkt()
#         except Exception:
#             pass

#         ## Count
#         try:
#             self.infos.numpts = layer.GetFeatureCount()
#         except Exception:
#             pass
            
#         ## SRS
#         if self.src_srs is None:
#             spatial_ref = layer.GetSpatialRef()
#             if spatial_ref:
#                 self.infos.src_srs = spatial_ref.ExportToProj4()
#         else:
#             self.infos.src_srs = self.src_srs

#         return self.infos

    
#     def _get_layer(self, ds):
#         """Internal helper to resolve the OGR Layer to process."""
        
#         layer = None
        
#         ## By Index/Name from Config
#         if self.ogr_layer is not None:
#             layer = ds.GetLayer(self.ogr_layer)
            
#         ## Auto-detect Known Names
#         if layer is None:
#             for lname in self._known_layer_names:
#                 layer = ds.GetLayerByName(lname)
#                 if layer: break
        
#         ## Default to first layer
#         if layer is None:
#             layer = ds.GetLayer(0)
            
#         return layer

    
#     def _resolve_fields(self, layer_defn):
#         """Internal helper to auto-detect field names if not provided."""
        
#         field_count = layer_defn.GetFieldCount()
#         field_names = [layer_defn.GetFieldDefn(i).GetName() for i in range(field_count)]
        
#         ## Elevation
#         if self.elev_field is None:
#             for f in self._known_elev_fields:
#                 if f in field_names:
#                     self.elev_field = f
#                     break
                    
#         ## Weight (No auto-detect)
#         ## Uncertainty (No auto-detect)


#     def yield_points(self):
#         """Yield points from the OGR datasource."""
        
#         if self.fn is None: return

#         ds_ogr = ogr.Open(self.fn)
#         if ds_ogr is None:
#             utils.echo_warning_msg(f"Could not open OGR file: {self.fn}")
#             return

#         layer = self._get_layer(ds_ogr)
#         if layer is None: return

#         ## Auto-detect fields
#         self._resolve_fields(layer.GetLayerDefn())

#         ## --- Spatial Filter Logic ---
#         if self.region is not None:
#             ## Determine the Check Region (User ROI)
#             ## Use trans_region if available (often handles pre-calc transforms), else user region
#             check_region = self.transform['trans_region'] if self.transform['trans_region'] else self.region
            
#             if check_region:
#                 ## Get Layer Native SRS
#                 layer_srs = layer.GetSpatialRef()
                
#                 ## Create Filter Geometry
#                 filter_geom = ogr.CreateGeometryFromWkt(check_region.export_as_wkt())
                
#                 ## Reproject Filter to Layer SRS if needed
#                 if layer_srs:
#                     ## Determine SRS of the Check Region
#                     filter_srs_str = check_region.src_srs if check_region.src_srs else 'epsg:4326'                    

#                     filter_srs = srsfun.osr_srs(filter_srs_str)

#                     if filter_srs and not layer_srs.IsSame(filter_srs):
#                         try:
#                             ## Create Transform: Region SRS -> Layer SRS
#                             transform = osr.CoordinateTransformation(filter_srs, layer_srs)
#                             filter_geom.Transform(transform)
#                         except Exception as e:
#                             ## If transform fails, warn but proceed (might filter incorrectly or empty)
#                             if self.verbose:
#                                 utils.echo_warning_msg(f"Failed to warp spatial filter to layer SRS: {e}")

#                 ## Apply Filter
#                 layer.SetSpatialFilter(filter_geom)

#         ## --- Buffer Setup for Chunking ---
#         chunk_x = []
#         chunk_y = []
#         chunk_z = []
#         chunk_w = []
#         chunk_u = []
#         chunk_size = 0
#         max_chunk = 100000

#         def flush_chunk():
#             if chunk_size == 0: return None
            
#             dataset = np.column_stack((chunk_x, chunk_y, chunk_z, chunk_w, chunk_u))
#             rec_arr = np.rec.fromrecords(dataset, names='x, y, z, w, u')
            
#             ## Apply Z Scale
#             if self.z_scale is not None:
#                 rec_arr['z'] *= self.z_scale
            
#             return rec_arr

#         ## --- Feature Iteration ---
#         for feature in layer:
#             geom = feature.GetGeometryRef()
#             if geom is None: continue

#             ## Get Weight/Uncertainty (Per Feature)
#             w_val = 1.0
#             if self.weight_field:
#                 w_val = utils.float_or(feature.GetField(self.weight_field), 1.0)
            
#             u_val = 0.0
#             if self.uncertainty_field:
#                 u_val = utils.float_or(feature.GetField(self.uncertainty_field), 0.0)

#             ## Get Coordinates
#             pts = []
            
#             ## Simple Geometry
#             if geom.GetGeometryCount() == 0:
#                 pts = geom.GetPoints()
#             else:
#                 ## Flatten Multi-Geometries
#                 for i in range(geom.GetGeometryCount()):
#                     sub_geom = geom.GetGeometryRef(i)
#                     pts.extend(sub_geom.GetPoints())

#             if not pts: continue

#             ## Extract Z
#             is_3d = geom.Is3D()
            
#             for pt in pts:
#                 x, y = pt[0], pt[1]
#                 z = None
                
#                 ## 3D Geometry Z
#                 if is_3d and len(pt) > 2:
#                     z = pt[2]
                
#                 ## Attribute Field
#                 if z is None and self.elev_field:
#                     z = utils.float_or(feature.GetField(self.elev_field))
                
#                 ## Constant Value
#                 if z is None and self.elevation_value is not None:
#                     z = self.elevation_value
                
#                 if z is None: continue
                
#                 chunk_x.append(x)
#                 chunk_y.append(y)
#                 chunk_z.append(z)
#                 chunk_w.append(w_val)
#                 chunk_u.append(u_val)
#                 chunk_size += 1

#                 if chunk_size >= max_chunk:
#                     yield flush_chunk()
#                     chunk_x, chunk_y, chunk_z, chunk_w, chunk_u = [], [], [], [], []
#                     chunk_size = 0

#         ## Final flush
#         if chunk_size > 0:
#             yield flush_chunk()

#         ds_ogr = None
        

# ### End
