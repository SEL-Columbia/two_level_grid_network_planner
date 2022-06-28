#!/usr/bin/env python
"""
Minimum Spanning Tree Alg
Ayse Selin Kocaman
ask2170@columbia.edu
"""

import os
import sys
import getopt
import time
import copy
import gc
import collections
#import scipy
#import pylab
#import numpy
#import batchPrimsforTransformers
from heapq import heappush, heappop
from osgeo import ogr
import network
import fileRW


class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

class Error(Exception):
    def __init__(self, msg):
        self.msg = msg

def generateDictsFromShp(shapeFile,outputPath):
    'Reads nodes and node weights from a point shapefile.'
    rootDir, fc = os.path.split(shapeFile)
    file, ext = os.path.splitext(fc)

    if not os.path.exists(outputPath):
        try:
            os.mkdir(outputPath)
        except:
            print "ERROR: could not create new directory", outputPath
    ds = ogr.Open(shapeFile)
    ptLayer = ds.GetLayer(0)

    nodes={}
    feat = ptLayer.GetNextFeature()
    while feat is not None:
        nodeWeight = 1
        geomRef = feat.GetGeometryRef()
        x = geomRef.GetX()
        y = geomRef.GetY()
        FID = feat.GetFID()
        nodes[FID] = network.Node(FID, x, y, nodeWeight) #Households
        feat = ptLayer.GetNextFeature()
    ds.Destroy()
    return nodes

def generateSegments(centers):
    segments=[]
    nodeCopy = centers.copy()

    segID=0
    for startNode in centers.values():
        del nodeCopy[startNode.getID()]
        for endNode in nodeCopy.values():
            dist = ((startNode.getX() - endNode.getX())**2 +
                    (startNode.getY() - endNode.getY())**2)**(.5)
            #if dist < searchRadius:

            segments.append(network.Seg(segID, startNode, endNode, dist))
            segID+=1
    return segments

def kruskalsAlg(segments,nodes):
    'Kruskal\'s algorithm for finding a minimum spanning tree'
    segments.sort(key=lambda obj:obj.getWeight())
    tree = network.Network()
    numNodes=len(nodes)

    for segment in segments:
        node1 = segment.getNode1()
        node2 = segment.getNode2()
        node1InNet = tree.inNet(node1)
        node2InNet = tree.inNet(node2)

        if (not node1InNet and not node2InNet) or (node1InNet != node2InNet):
            tree.addSeg(segment)
        else:
             if node1InNet and node2InNet and \
                    (tree.getNetID(node1) != tree.getNetID(node2)):
                        tree.addSeg(segment)
        if tree.numNodes() > numNodes:
            break
    return tree

def primsAlg(segments, numNodes, firstNodeID, excludedNodes = [],
        nodeDict = None):
    'Prim\'s Algorithm for finding a minimum spanning tree'

    # How to start from old tree
        # make parameter with existing tree (instance of network.Network)
        # can use old buildAssocDict (even though it's not optimal)
        # when getting the best first segment (algorithm starter segment),
            # add a check to make sure the first segment is not in the old tree

    tree = network.Network()
    if nodeDict == None:
        nodeDict = buildAssocDict(segments)
    segHeap = []

    # Find the shortest segment emanating from the node with the firstNodeID
    try:
        segs = nodeDict[firstNodeID]
    except KeyError:
        return tree

    leastWeight = None
    for seg in segs:
        if (seg.getWeight() < leastWeight) or (leastWeight == None):
            leastWeight = seg.getWeight()
            firstSeg = seg
    tree.addSeg(firstSeg)

    # Starter to algorithm
    # Add the segs emanating from the first two endpoints to the heap
    for endNode in [firstSeg.getNode1(), firstSeg.getNode2()]:
        addToHeap(segHeap, nodeDict[endNode.getID()])

    # Pick best from heap and repeat
    while tree.numNodes() < numNodes:
        try:
            # Get best segment from heap
            seg = heappop(segHeap)
        except:
            # Tree is finished (not all nodes contained).
            break
        node1, node2 = seg.getNodes()
        node1InNet = tree.inNet(node1)
        node2InNet = tree.inNet(node2)
        # Add the seg if it's terminal node isn't already in the cluster.
        if (not node1InNet) or (not node2InNet):
            if not node1InNet:
                endNode = node1
            else:
                endNode = node2
            tree.addSeg(seg)
            # Add all emanating segs to the heap:
            # nodeDict returns all segments coming out from the endNode
            # endNode is the node that is outside of the tree
            addToHeap(segHeap, nodeDict[endNode.getID()])
            # And we are sure that everything in the heap is adjacent to the tree because
            # we only add the adjacent segments in the first place using nodeDict
    return tree

def addToHeap(heap, newSegs):
    'Adds new segments to the segHeap.'
    for seg in newSegs:
        heappush(heap, seg)
    return heap

def buildAssocDict(segments):
    'Builds dictionary with nodeID key where values are segs from/to that node'
    segList = {}
    for seg in segments:
        node1, node2 = seg.getNodes()
        for nodeID in [node1.getID(), node2.getID()]:
            if segList.has_key(nodeID):
                segList[nodeID].append(seg)
            else:
                segList[nodeID] = [seg]
    return segList



def main(argv=None):
    if argv is None:
        argv = sys.argv
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "h", ["help"])
        except getopt.GetoptError, msg:
            raise Usage(msg)
        try:
            #inputShapeFile = argv[1]

            inputShapeFile=r"C:\Users\Selin\Documents\test\random100.shp"

            #outputDir = argv[2]
            outputDir=r"C:\Users\Selin\Documents\test\random100"
            #algorithm = argv[3]
            searchRadius = 6000 # in meters

        except IndexError:
            raise Error("Not enough arguments provided to script.")

        nodes=generateDictsFromShp(inputShapeFile,outputDir)
        numNodes=len(nodes)
        segments= generateSegments(nodes)
        tree= primsAlg(segments, numNodes, firstNodeID=0, excludedNodes = [],nodeDict = None)
        fileRW.genShapefile(tree, outputDir + ".prj", outputDir + os.sep + "MST.shp")


    except Usage, err:
        print >>sys.stderr, err.msg

        return 2
    except Error, err:
        print >>sys.stderr, "ERROR:", err.msg
        return 1

if __name__ == "__main__":
    sys.exit(main())
