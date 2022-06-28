#!/usr/bin/env python
"""
Chooses either fileRW-arcgis.py or fileRW-ogr.py to handle shapefile i/o.
Returns an error if neither is available.

Alex Zvoleff, aiz2101@columbia.edu
"""


try:
    from fileRWogr import * 
except:
    try:
        from fileRWarcgis import *   
    except ImportError:
        print("ERROR: failed to load ArcGIS or OGR. Cannot process shapefiles.")

def readNodesFromTxt(filename):
    'Reads nodes in from a text file and returns them as a dictionary by FID.'
    try:
        fid = open(filename,"r")
    except:
        print("ERROR: Could not load " + filename)
        return 1
    line = fid.readline()
    nodes = {}
    while "END" not in line:
        nodeID, x, y, weight = line.split()
        nodeID = int(nodeID)
        x = float(x)
        y = float(y)
        weight = float(weight)
        nodes[nodeID] = network.Node(nodeID, x, y, weight)
        line = fid.readline()
    fid.close()
    return nodes
    
def readSegsFromTxt(filename, nodes, cutoff=None):
    'Reads segments from a text file and them as a list of Seg instances.'
    try:
        fid = open(filename,"r")
    except:
        print("ERROR: Could not load " + filename)
        return 1
    segments = []
    line = fid.readline()
    segID = 0
    while line:
        try:
            node1ID, node2ID, weight = line.split()
        except ValueError:
            # Handles the case of node lists without stored seg weights
            node1ID, node2ID = line.split()
            weight = -1
        node1ID = int(node1ID)
        node2ID = int(node2ID)
        if nodes != None:
            node1 = nodes[node1ID]
            node2 = nodes[node2ID]
        else:
            node1 = node1ID
            node2 = node2ID
        weight = float(weight)
        if cutoff ==  None or weight <= cutoff:
            segments.append(network.Seg(segID, node1, node2, weight))
            segID += 1
        line = fid.readline()
    fid.close()
    return segments

def writeNodesToTxt(nodes, outputFile):
    'Writes input node instances to a text file.'
    ofile = open(outputFile, "w")
    for node in nodes.values():
        id = node.getID()
        x, y = node.getX(), node.getY()
        weight = node.getWeight()
        ofile.write("%(id)i %(x).11e %(y).11e %(weight)f\n" %vars())
    ofile.write("END\n")
    ofile.close()
    return 0
 
def writeSegsToTxt(dists, outputFile):
    'Writes distances to a text file.'
    ofile = open(outputFile, "w")
    for dist in dists:
        p1 = dist[0]
        p2 = dist[1]
        dist = float(dist[2])
        ofile.write("%(p1)i %(p2)i %(dist)f\n" %vars())
    ofile.close()
    return 0

def writeNodesToPP(nodes, outputFile):
    'Writes input node instances to a text file in R point-process format.'
    ofile = open(outputFile, "w")
    # Calculate bounding box:
    xvals = []
    yvals = []
    for node in nodes.values():
        xvals.append(node.getX())
        yvals.append(node.getY())
    xl = min(xvals)
    xu = max(xvals)
    yl = min(yvals)
    yu = max(yvals)
    ofile.write(str(len(nodes)))
    ofile.write("\n***R spatial point process***\n" %vars())
    ofile.write("%(xl).11e %(xu).11e %(yl).11e %(yu).11e 1\n" %vars())
    for x, y in zip(xvals, yvals):
        ofile.write("%(x).11e %(y).11e\n" %vars())
    ofile.close()
    return 0
