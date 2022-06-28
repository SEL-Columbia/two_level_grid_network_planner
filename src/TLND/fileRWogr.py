#!/usr/bin/env python
"""
Functions for reading coordinate and segment files from ArcGIS, and for
converting generated networks into an ArcGIS coverage or shapefile. Also
handles writing network statistics fields to a given shapefile. Uses OGR
to process shapefiles.

Alex Zvoleff, aiz2101@columbia.edu
"""

import sys
import os

import network

try:
    from osgeo import ogr
except ImportError:
    import ogr

try:
    from osgeo import osr
except ImportError:
    import osr

def convertToShpUnits(shapefile, input):
    'Converts input units to match linear units of datasource.'
    #TODO: Fix this function so it actually handles conversions, rather
    # than just throwing an error
    inputLength, inputUnits = input.split()
    inputUnits = inputUnits.lower()

    projFile = shapefile[:-3] + 'prj'
    projText = []
    inFile = open(projFile, 'r')
    for line in inFile:
        projText.append(line)
    inFile.close()

    spatialRef = osr.SpatialReference()
    spatialRef.ImportFromESRI(projText)
    shpUnits = spatialRef.GetLinearUnitsName().lower()
    if shpUnits != inputUnits and shpUnits + 's' != inputUnits:
        print("ERROR: shapefile units do not match input units")
        sys.exit(1)
    else:
        return float(inputLength)

def getFieldType(fieldValue):
    'Returns OGR field type appropriate for the given value'
    if type(fieldValue) == float:
        return ogr.OFTReal
    elif type(fieldValue) == int:
        return ogr.OFTInteger
    elif type(fieldValue) == str:
        return ogr.OFTString

def readNodesFromShp(shapefile):
    'Reads nodes and node weights from a point shapefile.'
    ds = ogr.Open(shapefile)
    ptLayer = ds.GetLayer(0)
    nodes = {}
    feat = ptLayer.GetNextFeature()
    while feat is not None:
        try:
            weightField = feat.GetFieldIndex("Weight")
            nodeWeight = feat.GetField(weightField)
        except:
            nodeWeight = 1
        geomRef = feat.GetGeometryRef()
        x = geomRef.GetX()
        y = geomRef.GetY()
        FID = feat.GetFID()
        nodes[FID] = network.Node(FID, x, y, nodeWeight)
        feat = ptLayer.GetNextFeature()
    ds.Destroy()
    return nodes

def readNetFromShp(shapefile):
    'Reads segs and nodes from the given shapefile'
    ds = ogr.Open(shapefile)
    layer = ds.GetLayer(0)
    net = network.Network()
    feat = layer.GetNextFeature()
    lengthField = "Length"
    nodeWeightFields = ["pt1Weight", "pt2Weight"]
    nodeIDFields = ["pt1", "pt2"]
    while feat is not None:
        geomRef = feat.GetGeometryRef()
        length = feat.GetField(lengthField)

        endPts = []
        for n in xrange(2):
            x, y = geomRef.GetX(n), geomRef.GetY(n)
            try:
                nodeWeight = feat.GetField(nodeWeightFields[n])
                print("try nodefield", nodeWeight)
            except ValueError as msg:
                print(msg)
                print("ERROR: field \""  + nodeWeightFields[n] + "\" doesn't exist")
            nodeID = feat.GetField(nodeIDFields[n])
            endPts.append(network.Node(nodeID, x, y, nodeWeight))
            print("node weight", nodeWeight)
        newSeg = network.Seg(feat.GetFID(), endPts[0], endPts[1], length)
        net.addSeg(newSeg)
        feat = layer.GetNextFeature()
    ds.Destroy()
    return net

def writeFieldToShp(shapefile, fieldValues, field):
    'Writes a field (provided as a dictionary by FID) to a shapefile.'
    #TODO: fix to allow creation of a new field in a shapefile that already has
    #features... if OGR will allow it. Can copy into a new shapefle as in the
    #soundg examples.
    ds = ogr.Open(shapefile)
    layer = ds.GetLayer(0)
    feat = layer.GetNextFeature()
    fieldIndex = feat.GetFieldIndex(field)
    if fieldIndex == -1:
        fieldType = getFieldType(fieldValues.values()[0])
        fieldDefn = ogr.FieldDefn(field, fieldType)
        layer.CreateField(fieldDefn)
    while feat is not None:
        FID = feat.GetFID()
        try:
            fieldValue = fieldValues[FID]
        except KeyError:
            pass
        else:
            feat.SetField(fieldIndex, fieldValue)
            feat = layer.GetNextFeature()
    ds.Destroy()
    return 0

def readFieldFromShp(shapefile, field):
    'Reads field values from a shapefile and returns them as a dictionary by FID.'
    ds = ogr.Open(shapefile)
    layer = ds.GetLayer(0)
    feat = layer.GetNextFeature()
    fieldValues = {}
    while feat is not None:
        try:
            fieldValue = feat.GetField(field)
        except ValueError as msg:
            print(msg)
            print("ERROR: field \""  + field + "\" doesn't exist")
        FID = feat.GetFID()
        fieldValues[FID] = fieldValue
        feat = layer.GetNextFeature()
    ds.Destroy()
    return fieldValues

def genShapefile(network, projFile, outputShape):
    'Generates a shapefile from a network class instance.'
    if not os.path.exists(projFile): 
        print("ERROR:", projFile, "does not exist")
        return 1
    elif os.path.exists(outputShape): 
        print("ERROR:", outputShape, "already exists")
        return 1
    shpDriver = ogr.GetDriverByName('ESRI Shapefile')
    ds = shpDriver.CreateDataSource(outputShape)
    
    projText = []
    inFile = open(projFile, 'r')
    for line in inFile:
        projText.append(line)
    inFile.close()

    spatialRef = osr.SpatialReference()
    spatialRef.ImportFromESRI(projText)

    layer = ds.CreateLayer('network', srs=spatialRef, geom_type=ogr.wkbLineString)

    # setup fields
    lengthFD = ogr.FieldDefn("Length", ogr.OFTReal)
    layer.CreateField(lengthFD)
    netidFD = ogr.FieldDefn("netID", ogr.OFTString)
    layer.CreateField(netidFD)
    pt1idFD = ogr.FieldDefn("pt1", ogr.OFTInteger)
    layer.CreateField(pt1idFD)
    pt2idFD = ogr.FieldDefn("pt2", ogr.OFTInteger)
    layer.CreateField(pt2idFD)
    pt1WeightFD = ogr.FieldDefn("pt1Weight", ogr.OFTReal)
    layer.CreateField(pt1WeightFD)
    pt2WeightFD = ogr.FieldDefn("pt2Weight", ogr.OFTReal)
    layer.CreateField(pt2WeightFD)

    # write features to shapefile

    for segment in network.getEdges():
        node1, node2 = segment.getNodes()
        newFeature = ogr.Feature(feature_def = layer.GetLayerDefn())

        newFeature.SetField("netID", network.getNetID(node1))
        newFeature.SetField("Length", segment.getWeight())
        newFeature.SetField("pt1", node1.getID())
        newFeature.SetField("pt2", node2.getID())
        newFeature.SetField("pt1Weight", node1.getWeight())
        newFeature.SetField("pt2Weight", node2.getWeight())
        
        x1, y1 = node1.getX(), node1.getY()
        x2, y2 = node2.getX(), node2.getY()
        wkt = 'LINESTRING (%f %f, %f %f)' %(x1, y1, x2, y2)
        geom = ogr.CreateGeometryFromWkt(wkt)

        newFeature.SetGeometry(geom)

        layer.CreateFeature(newFeature)
    ds.Destroy()
