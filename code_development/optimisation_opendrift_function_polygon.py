# -*- coding: utf-8 -*-
"""
Created on Wed May 26 10:29:19 2021

@author: Romain
"""

#!/usr/bin/python

# Import necessary modules
import shapefile
from shapely.geometry import Polygon, Point, MultiPolygon
import time

# Import release polygons
shp = shapefile.Reader('./habitat/rock_lobster_polygons.shp')
bins = shp.shapes()
records = shp.records()

# Define two points to test the algorithm:
#p1 = Point(174.71091, -38.07374) # Centroid of one of the polygons
#p2 = Point(175.71091, -37.07374) # Pt not inside of any polygon of the shapefile
lons = [174.71091, 174.71091, 173.7666]
lats = [-38.07374, -37.07374, -39.1275]

## Current code:
# Loop to check within polygon or not
#def within_habitat(shp, current_lat, current_lon)
start = time.time()
for i in range(len(lons)):
    pt = Point(lons[i], lats[i])
    for index, item in enumerate(shp):
        pts = bins[index].points # Access one polygon
        #pts = [list(elem) for elem in pts]
        poly = Polygon(pts)
        in_area = pt.within(poly)
        if in_area == True:
            #print("The habitat id is", index)
            print("In habitat")
        #else:
            #print("No habitat")
# get time taken to run the for loop code         
elapsed_time_fl = (time.time() - start)
print(elapsed_time_fl) # 0.2 sec for 3 particles


## Remove loop
import fiona
polyShp = fiona.open('./habitat/rock_lobster_polygons_fixed.shp')

polyList = []
polyProperties = []
for poly in polyShp:
    polyGeom = Polygon(poly['geometry']['coordinates'][0])
    polyList.append(polyGeom)
    polyProperties.append(poly['properties'])
#print(polyList[10])
#print(polyProperties[10])
multiShp = MultiPolygon(polyList)
multiShp = multiShp.buffer(0)
#print(multiShp.is_valid)
#print(type(multiShp))

start = time.time()
for i in range(len(lons)):
    pt = Point(lons[i], lats[i])
    in_area = pt.within(multiShp)
    if in_area == True:
        print("In habitat")
    else:
        print("No habitat")
# get time taken to run       
elapsed_time_fl = (time.time() - start)
print(elapsed_time_fl) # 0.01 sec for the 3 particles

