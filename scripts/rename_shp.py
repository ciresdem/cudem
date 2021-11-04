#!/usr/bin/env python

import sys
from osgeo import ogr

inshp = sys.argv[1]
#newfld = sys.argv[2]

# Open a Shapefile, and get field names
source = ogr.Open(inshp, 1)
layer = source.GetLayer()
layer_defn = layer.GetLayerDefn()
field_names = [layer_defn.GetFieldDefn(i).GetName() for i in range(layer_defn.GetFieldCount())]
#print len(field_names)

# Add a new field
new_field = ogr.FieldDefn('Name', ogr.OFTString)
#f = ogr.Feature(layer_defn)
layer.CreateField(new_field)
#new_field.SetField(1, "test")
#layer.SetField(1, "test")
#layer.SetFeature(

for i in range(layer.GetFeatureCount()):
    feat = layer.GetFeature(i)
    feat.SetField('Name', inshp[:-4])
    #print i
    layer.SetFeature(feat)
    feat = None

# Close the Shapefile
source = None
