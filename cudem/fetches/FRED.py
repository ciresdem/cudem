### FRED.py
##
## Copyright (c) 2012 - 2026 CIRES Coastal DEM Team
##
## FRED.py is part of CUDEM
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
## Fetches Remote Elevation Datalist (FRED)
##
## The fetches reference vector location and related functions
## for generating and parsing FRED.
##
### Code:

import os
import json
from typing import List, Dict, Optional, Any
from osgeo import ogr
from cudem import utils
from cudem import regions

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
FETCH_DATA_DIR = os.path.join(THIS_DIR, 'data')

class FRED:
    def __init__(self, name: str = 'FRED', verbose: bool = False, local: bool = False):
        self.verbose = verbose
        self.name = name
        self.driver = ogr.GetDriverByName('GeoJSON')
        self.filename = f'{name}.geojson'
        
        ## Determine file location
        if local:
            self.path = self.filename
        elif os.path.exists(self.filename):
            self.path = self.filename
        elif os.path.exists(os.path.join(FETCH_DATA_DIR, self.filename)):
            self.path = os.path.join(FETCH_DATA_DIR, self.filename)
        else:
            self.path = self.filename
            
        if self.verbose:
            utils.echo_msg(f'using {self.path}')
            
        self.ds = None
        self.layer = None
        self.is_open = False
        
        self._fields = {
            'Name': ogr.OFTString,
            'ID': ogr.OFTString,
            'Date': ogr.OFTInteger,
            'Agency': ogr.OFTString,
            'MetadataLink': ogr.OFTString,
            'MetadataDate': ogr.OFTString,
            'DataLink': ogr.OFTString,
            'IndexLink': ogr.OFTString,
            'Link': ogr.OFTString,
            'DataType': ogr.OFTString,
            'DataSource': ogr.OFTString,
            'Resolution': ogr.OFTString,
            'HorizontalDatum': ogr.OFTString,
            'VerticalDatum': ogr.OFTString,
            'LastUpdate': ogr.OFTString,
            'Etcetra': ogr.OFTString,
            'Info': ogr.OFTString
        }

        
    def __enter__(self):
        self._open_ds(mode=1)  # Open in update mode by default for context manager
        return self

    
    def __exit__(self, exc_type, exc_value, traceback):
        self._close_ds()

        
    def _create_ds(self):
        """Creates a new OGR DataSource and Layer."""
        
        utils.remove_glob(self.path)
        self.ds = self.driver.CreateDataSource(self.path)        
        self.layer = self.ds.CreateLayer(self.name, None, ogr.wkbMultiPolygon)
        
        ## Create fields dynamically from the schema
        for field_name, field_type in self._fields.items():
            self.layer.CreateField(ogr.FieldDefn(field_name, field_type))

            
    def _open_ds(self, mode: int = 0) -> int:
        """Opens the Data Source.
        Mode 0: Read-only
        Mode 1: Read/Write
        """
        
        if not self.is_open:
            try:
                self.ds = self.driver.Open(self.path, mode)
            except Exception:
                self.ds = None

            if self.ds is None or self.ds.GetLayerCount() == 0:
                self.ds = None
                self._create_ds()
                
            self.layer = self.ds.GetLayer()
            self.is_open = True
            return 0
        return -1

    
    def _close_ds(self):
        """Closes the Data Source and clears memory."""
        
        self.layer = None
        self.ds = None
        self.is_open = False
        return 0

    
    def _get_fields(self) -> List[str]:
        """Returns a list of field names in the current layer."""
        
        if self.is_open:
            schema = []
            ldefn = self.layer.GetLayerDefn()
            for n in range(ldefn.GetFieldCount()):
                fdefn = ldefn.GetFieldDefn(n)
                schema.append(fdefn.name)
            return schema
        return []

    
    def _add_feature(self, attributes: Dict[str, Any], geometry_json: str) -> int:
        """Add a feature (attributes + geometry) to the reference vector layer."""
        
        if not self.is_open:
            return -1

        layer_defn = self.layer.GetLayerDefn()
        feat = ogr.Feature(layer_defn)
        
        ## Set Geometry
        geom = ogr.CreateGeometryFromJson(geometry_json)
        if geom:
            geom_valid = geom.MakeValid()
            feat.SetGeometry(geom_valid)
        
        ## Set Fields
        for field in self._fields.keys():
            val = attributes.get(field)

            if val is not None:
                feat.SetField(field, val)
            else:
                feat.SetField(field, -1) 

        self.layer.CreateFeature(feat)
        feat = None
        return 0

    
    def _add_survey(self, geom=None, **kwargs):
        """Add a single survey entry to the FRED database.
        Accepts keywords matching the schema (Name, ID, Agency, etc.)
        """
        
        if not self.is_open:
            return None
            
        if geom is not None and self.layer is not None:
            ## Prepare attributes
            attrs = kwargs.copy()
            attrs['LastUpdate'] = utils.this_date()
            
            # # Ensure all schema keys exist, default to None if missing
            # for key in self._fields.keys():
            #     if key not in attrs:
            #         attrs[key] = None

            self._add_feature(attrs, geom.ExportToJson())
            for feature in self.layer:
                self.layer.SetFeature(feature)
        else:
            return None

        
    def _add_surveys(self, surveys: List[Dict]):
        """Update or create a reference vector using a list of surveys."""
        
        if self.is_open and self.layer is not None:
            for survey in surveys:
                ## Unpack survey dict into kwargs
                self._add_survey(**survey)
            return 0
        return -1

    
    def _edit_feature(self, feature, survey_data: tuple):
        """Edit an existing feature."""
        
        if self.is_open:
            attributes, geom_json = survey_data
            geom = ogr.CreateGeometryFromJson(geom_json)
            feature.SetGeometry(geom)
            
            for field in self._fields.keys():
                try:
                    feature.SetField(field, attributes[field])
                except KeyError:
                    feature.SetField(field, -1)
            
            self.layer.SetFeature(feature)
            return 0
        return -1

    
    def _attribute_filter(self, where: List[str] = []) -> int:
        """Apply an Attribute Filter (SQL WHERE clause style)."""
        
        if self.is_open:
            if where:
                wf = ' AND '.join(where)
                self.layer.SetAttributeFilter(wf)
            else:
                self.layer.SetAttributeFilter(None)
            return 0
        return -1

    
    def _get_region(self, where: List[str] = [], layers: List[str] = []) -> Any:
        """Calculate the bounding region for specific filters."""
        
        out_regions = []
        for i, layer_name in enumerate(layers):
            self.layer.SetAttributeFilter(f"DataSource = '{layer_name}'")
            
            ## Apply additional where filters
            current_filter = self.layer.GetLayerDefn() 
            if where:
                existing_filter = self.layer.GetAttributeFilter()
                combined = f"({existing_filter}) AND ({' AND '.join(where)})" if existing_filter else ' AND '.join(where)
                self.layer.SetAttributeFilter(combined)

            for feat in self.layer:
                geom = feat.GetGeometryRef()
                if geom:
                    wkt = geom.ExportToWkt()
                    this_region = ogr.CreateGeometryFromWkt(wkt).GetEnvelope()
                    if len(out_regions) > 0:
                        out_regions = regions.regions_merge(out_regions, this_region)
                    else:
                        out_regions = this_region
        return out_regions

    
    def _filter(self, region=None, where: List[str] = [], layers: List[str] = []) -> List[Dict]:
        """Search for data in the reference vector file."""
        
        _results = []
        
        ## Handle geometry filter
        _boundsGeom = region.export_as_geom() if region is not None else None

        ## Context management for opening/closing
        close_on_exit = False
        if not self.is_open:
            self._open_ds()
            close_on_exit = True

        try:
            ## Iterate through requested layers (DataSources)
            with utils.ccp(total=len(layers), desc=f'filtering {self.path}', leave=self.verbose) as pbar:
                for layer_name in layers:
                    pbar.update(1)
                    
                    ## Construct filter
                    layer_where = where.copy()
                    layer_where.append(f"DataSource = '{layer_name}'")
                    
                    if self.verbose:
                        utils.echo_msg(f'FRED region: {region}')
                        utils.echo_msg(f'FRED filter: {layer_where}')

                    self._attribute_filter(where=layer_where)
                    
                    for feat in self.layer:
                        ## Spatial Filter
                        if _boundsGeom is not None:
                            geom = feat.GetGeometryRef()
                            if geom and not _boundsGeom.Intersects(geom):
                                continue 
                        
                        ## Process matching feature
                        f_json = json.loads(feat.ExportToJson())
                        properties = f_json.get('properties', {})
                        
                        ## Ensure we get actual field values from the feature object
                        ## (properties dict from ExportToJson is usually sufficient, 
                        ## but GetField ensures strict type casting if OGR layer is defined)
                        result_dict = {}
                        for key in properties.keys():
                            result_dict[key] = feat.GetField(key)
                        _results.append(result_dict)

        finally:
            if close_on_exit:
                self._close_ds()
            
        return _results

## ==============================================
## Lambdas for the FRED using the module object `mod`
## These lambdas assume `mod` is an object with attributes: FRED, region, where, name.
## ==============================================
_filter_FRED = lambda mod: mod.FRED._filter(region=mod.region, where=mod.where, layers=[mod.name])
_update_FRED = lambda mod, s: mod.FRED._add_surveys(s)
_filter_FRED_index = lambda mod: [utils.echo_msg(json.dumps(f, indent=2)) for f in _filter_FRED(mod)]

### End
