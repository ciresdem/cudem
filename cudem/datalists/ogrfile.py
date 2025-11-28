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

    
    def yield_points(self):
        ds_ogr = ogr.Open(self.fn)
        count = 0
        #utils.echo_msg(self.ogr_layer)
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
                                
                            if isinstance(xyz[0], list):
                                for x in xyz:
                                    #for xx in x:
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
