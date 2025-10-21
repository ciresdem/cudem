#!/usr/bin/env python

from osgeo import ogr
import sys

driver = ogr.GetDriverByName('GPKG')
fn = sys.argv[1]
#f_name = sys.argv[2]
f_val = sys.argv[2]
f_val_rep = sys.argv[3]

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
        #try:
        if fval == f_val:
            print(fdefn.name, fval)
            feature.SetField(fdefn.name, f_val_rep)#fval.replace('"', ''))
        #except: pass
    layer.SetFeature(feature)
    feature = layer.GetNextFeature()

dataSource.Destroy()
