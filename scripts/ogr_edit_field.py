#!/usr/bin/env python

from osgeo import ogr
import sys

driver = ogr.GetDriverByName('ESRI Shapefile')
fn = sys.argv[1]
#print(fn)
dataSource = driver.Open(fn, 1)

layer = dataSource.GetLayer()
ldefn = layer.GetLayerDefn()
feature = layer.GetNextFeature()

dist = 233

while feature:
    for n in range(ldefn.GetFieldCount()):
        fdefn = ldefn.GetFieldDefn(n)
        fval = feature.GetField(fdefn.name)
        print(fdefn.name, fval)
        try:
            feature.SetField(fdefn.name, fval.replace('"', ''))
        except: pass
    layer.SetFeature(feature)
    feature = layer.GetNextFeature()

dataSource.Destroy()
