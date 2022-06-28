#!/usr/bin/env python
"""
A heuristic approach for two-level network design - rural electrification
Ayse Selin Kocaman
ask2170@columbia.edu
"""

import os
import sys
import getopt
import time
import copy
import CMST_dfs
import gc
import collections
#import scipy
#import pylab
import numpy as np
#import batchPrimsforTransformers
from heapq import heappush, heappop
from osgeo import ogr
import network
import fileRW
#import prims5

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

class Error(Exception):
    def __init__(self, msg):
        self.msg = msg


    
def mergeCluster(ClusterByNode,NodesByClusterID,Centers,segment): 
    center1,center2=segment.getNodes()
   
    centerX=(center1.getWeight()*center1.getCenterX()
        +center2.getWeight()*center2.getCenterX())/(center2.getWeight()+center1.getWeight())
    centerY=(center1.getWeight()*center1.getCenterY()
        +center2.getWeight()*center2.getCenterY())/(center2.getWeight()+center1.getWeight())

    
    weight=center2.getWeight()+center1.getWeight()
    baseClusterID=min(ClusterByNode[center1],ClusterByNode[center2])
    mergingClusterID=max(ClusterByNode[center1],ClusterByNode[center2])
   
    NodesByClusterID[baseClusterID].extend(NodesByClusterID.pop(mergingClusterID))
    
    Centers[baseClusterID].setXY(centerX,centerY)
    Centers[baseClusterID].setWeight(weight)
    
    del Centers[mergingClusterID] 
   
    for node in NodesByClusterID[baseClusterID]:
        ClusterByNode[node]=baseClusterID


        
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
    
    nodesByClusterID=collections.defaultdict(list)
    clusterByNode={}
    nodes={}
    centers={}
    LVCostDict={}
    
    feat = ptLayer.GetNextFeature()
    while feat is not None:
        nodeWeight = 1
        geomRef = feat.GetGeometryRef()
        x = geomRef.GetX()
        y = geomRef.GetY()
        FID = feat.GetFID()
        nodes[FID] = network.Node(FID, x, y, nodeWeight) #Households
        centers[FID]=network.Node(FID, x, y, nodeWeight) #Transformers (center of mass of the cluster)
        
        clusterByNode[nodes[FID]]=FID 
        nodesByClusterID[FID].append(nodes[FID])
        LVCostDict[FID]=0
        feat = ptLayer.GetNextFeature()
    ds.Destroy()
    return nodesByClusterID,clusterByNode,nodes,centers,LVCostDict

def generateSegments(centers,searchRadius): 
    segments=[]
    nodeCopy = centers.copy()
   
    segID=0
    for startNode in centers.values():
        del nodeCopy[startNode.getID()]
        for endNode in nodeCopy.values():
            dist = ((startNode.getX() - endNode.getX())**2 + 
                    (startNode.getY() - endNode.getY())**2)**(.5)
            if dist < searchRadius:
                
                segments.append(network.Seg(segID, startNode, endNode, dist))
                segID+=1
    return segments





def maxInClusterDist(centerNode,nodesByClusterID): #Returns maxDist within the cluster
    maxdist=0
    for node in nodesByClusterID[centerNode.getID()]: #uses the fact that centerID and ClusterID are same
        dist=((centerNode.getX()-node.getX())**2+
                (centerNode.getY()-node.getY())**2)**(.5)
        if dist>=maxdist:
            maxdist=dist
    return maxdist


def maxTempInClusterDist(segment,ClusterByNode,nodesByClusterID):
    maxDist=0

    tempCenter1,tempCenter2=segment.getNodes()
       
    tempCenterX=(tempCenter1.getWeight()*tempCenter1.getX()
        +tempCenter2.getWeight()*tempCenter2.getX())/(tempCenter2.getWeight()+tempCenter1.getWeight())
    tempCenterY=(tempCenter1.getWeight()*tempCenter1.getY()
        +tempCenter2.getWeight()*tempCenter2.getY())/(tempCenter2.getWeight()+tempCenter1.getWeight())
    
    for node in nodesByClusterID[ClusterByNode[segment.getNode1()]]:
        dist=((tempCenterX-node.getX())**2+(tempCenterY-node.getY())**2)**(.5)
        if dist>=maxDist:
            maxDist=dist

    
    for node in nodesByClusterID[ClusterByNode[segment.getNode2()]]:
        dist=((tempCenterX-node.getX())**2+(tempCenterY-node.getY())**2)**(.5)
        if dist>=maxDist:
            maxDist=dist

    return maxDist,tempCenterX,tempCenterY



    
def totalInClusterCost(nodesByClusterID,centers):
    totalCost=0
    for centerID in centers.keys():
        for node in nodesByClusterID[centerID]:
            totalCost+=((node.getX()-centers[centerID].getX())**2+
                        (node.getY()-centers[centerID].getY())**2)**(.5)
    return totalCost
    
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
    return tree,segments

def primsAlg(segments, numNodes, firstNodeID, nodeDict):
    'Prim\'s Algorithm for finding a minimum spanning tree'

    
    tree = network.Network()
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


   
def run(centers,nodesByClusterID,clusterByNode,LVCostDict,sr,MV,LV,TCost,distFromT,maxLVLenghtInCluster,outputDir):
    
    #print "First Stage starts without MST"
    sumLVCostAtEachStep={}
    minCenters=copy.deepcopy(centers)
    st=time.time()
    segments=generateSegments(minCenters,sr)
    
    
    # To write total cost to a text file
 #   statFile= outputDir + os.sep + "TotalCost_FirstStage2.txt"
#    outFile = open(statFile,"w")
    
    #minTree,segments=kruskalsAlg(segments,centers) # can use either Kruskal or Prims
    m1=time.time()
    #nodeDict=buildAssocDict(segments)
    
    #minTree=primsAlg(segments, len(minCenters), 0, nodeDict) # 0 is the starting node of Prims algorithm
   
    minTotalCost=len(centers)*TCost
    #minTotalCost=minTree.getTotalEdgeWeight()*MV+len(centers)*TCost
 #   outFile.write("%(minTotalCost)f\n" %vars())
    minLVCostDict=copy.deepcopy(LVCostDict)
    
    minNodesByClusterID=copy.deepcopy(nodesByClusterID)
    minClusterByNode=copy.deepcopy(clusterByNode)
    minSeg=min(segments,key=lambda obj:obj.getWeight())
    #fixes tempCenterX , tempCenterY and maxDist of the first merge before going into while
    if minSeg.getWeight()<=distFromT*2:
        maxDist=0 # can be anything less than 500
    else:
        maxDist=distFromT+10 # can be anything greater than 500
        print "NO CLUSTER POSSIBLE"
    
    tempCenter1,tempCenter2=minSeg.getNodes()
       
    tempCenterX=(tempCenter1.getWeight()*tempCenter1.getX()
        +tempCenter2.getWeight()*tempCenter2.getX())/(tempCenter2.getWeight()+tempCenter1.getWeight())
    tempCenterY=(tempCenter1.getWeight()*tempCenter1.getY()
        +tempCenter2.getWeight()*tempCenter2.getY())/(tempCenter2.getWeight()+tempCenter1.getWeight())
    i=len(centers)
    
    while(maxDist<=distFromT):
        i-=1
        
        
        #if i%20==0:
            #print i
        center1,center2=minSeg.getNodes()

	
        weight=center2.getWeight()+center1.getWeight()
        baseClusterID=min(clusterByNode[center1],clusterByNode[center2])
        mergingClusterID=max(clusterByNode[center1],clusterByNode[center2])
   
        nodesByClusterID[baseClusterID].extend(nodesByClusterID.pop(mergingClusterID))
    
        centers[baseClusterID].setXY(tempCenterX,tempCenterY)
        centers[baseClusterID].setWeight(weight)
    
        del centers[mergingClusterID] 
        
        for node in nodesByClusterID[baseClusterID]:
            clusterByNode[node]=baseClusterID
        
        segments=generateSegments(centers,sr) # generate segments for new graph
        #nodeDict=buildAssocDict(segments)
        #newTree,segments=kruskalsAlg(segments,centers) #returns sorted segments
        #newTree=primsAlg(segments, len(centers), 0, nodeDict) # 0 is the starting node of prims.
        #TotalMVCost=newTree.getTotalEdgeWeight()*MV
        #TotalTransformerCost=len(centers)*TCost
        #del LVCostDict[mergingClusterID]
        gc.collect()
        #segmentsCMST, LVCostDict[baseClusterID] =CMST_dfs.CMST(nodesByClusterID[baseClusterID],maxLVLenghtInCluster,centers[baseClusterID])
        #newTotalCost=TotalMVCost+TotalTransformerCost+(sum(LVCostDict.values()))*LV
        #sumLVCostAtEachStep[len(centers)]=sum(LVCostDict.values())*LV
        #newTotalCost=TotalTransformerCost+(sum(LVCostDict.values()))*LV

        #outFile.write("%(newTotalCost)f\n" %vars())
        #outFile.write("%i %f\n" %(i,sumLVCostAtEachStep[len(centers)]))
        '''                      
        if(newTotalCost<=minTotalCost):
            minNodesByClusterID=copy.deepcopy(nodesByClusterID)
            #minTree=copy.deepcopy(newTree)
            minCenters=copy.deepcopy(centers)
            minLVCostDict=LVCostDict.copy()
            minTotalCost=newTotalCost
            minClusterByNode=copy.deepcopy(clusterByNode)
            minMaxDist=maxDist
            '''
        # Calculate maxDist below for next graph and continue if it is less than 500
        
        try: # to check if there is a segment on the graph or there is only one cluster  # bir tane break eden varsa bile devamini check ediyor!!!!!
            minSeg=min(segments,key=lambda obj:obj.getWeight()) 
	    maxDist,tempCenterX,tempCenterY=maxTempInClusterDist(minSeg,clusterByNode,nodesByClusterID)
	    if maxDist>distFromT:
		segments.sort(key=lambda obj:obj.getWeight())
	
		for seg in segments:
		    if seg.getWeight()>distFromT*2:
			break #break from for loop
		    else:
		        maxDist,tempCenterX,tempCenterY=maxTempInClusterDist(seg,clusterByNode,nodesByClusterID)    
                        if maxDist <=distFromT:
			    minSeg=seg
			    break #break from for loop
		
	
        except:
            break
        
    #outFile.close()
    '''
    print "Second Stage starts with MST"
       # Second stage with MST starts

    if len(minCenters)==len(centers) or len(minCenters)==1:
          segments_ST=generateSegments(centers,sr)
          nodeDict=buildAssocDict(segments_ST)
          minTree=primsAlg(segments_ST, len(centers), 0, nodeDict) # 0 is the starting node of Prims algorithm
          minTotalCost_ST=minTree.getTotalEdgeWeight()*MV+minTotalCost
          return minTotalCost_ST,minTree,centers,nodesByClusterID,sum(LVCostDict.values())*LV
    
    centers_ST=copy.deepcopy(minCenters)
    
    segments_ST=generateSegments(centers_ST,sr)
    
    # To write total cost to a text file
    statFile= outputDir + os.sep + "TotalCost_SecondStage2.txt"
    outFile = open(statFile,"w")
    
    #minTree,segments=kruskalsAlg(segments,centers) # can use either Kruskal or Prims
    m1=time.time()
    nodeDict=buildAssocDict(segments_ST)
    
    minTree=primsAlg(segments_ST, len(centers_ST), 0, nodeDict) # 0 is the starting node of Prims algorithm
   
    i=len(centers_ST)
    minTotalCost_ST=minTree.getTotalEdgeWeight()*MV+len(centers_ST)*TCost+(sum(minLVCostDict.values()))*LV
 #   outFile.write("%(minTotalCost_ST)f\n" %vars())
    outFile.write("%i %f %f %f\n" %(i,(sum(minLVCostDict.values()))*LV, minTree.getTotalEdgeWeight()*MV, minTotalCost_ST))
    minLVCostSum_ST=9999999999999999 #a big number
    nodesByClusterID_ST=copy.deepcopy(minNodesByClusterID)
    clusterByNode_ST=copy.deepcopy(minClusterByNode)

    try: # to check if there is a segment on the graph or there is only one cluster
            minSeg_ST=min(segments_ST,key=lambda obj:obj.getWeight()) 
	    maxDist,tempCenterX,tempCenterY=maxTempInClusterDist(minSeg_ST,clusterByNode_ST,nodesByClusterID_ST)
	    if maxDist>distFromT:
		segments_ST.sort(key=lambda obj:obj.getWeight())
	
		for seg in segments_ST:
		    if seg.getWeight()>distFromT*2:
			break #break from for loop
		    else:
		        maxDist,tempCenterX,tempCenterY=maxTempInClusterDist(seg,clusterByNode_ST,nodesByClusterID_ST)    
                        if maxDist <=distFromT:
			    minSeg_ST=seg
			    break #break from for loop
			    
    except:
            return minTotalCost_ST,minTree,centers_ST,nodesByClusterID_ST,sum(LVCostDict_ST.values())*LV

    
  #  centers_ST=copy.deepcopy(minCenters)
  
    #try:
#        minSeg_ST=min(segments_ST,key=lambda obj:obj.getWeight())
    #except:
#        return minTotalCost_ST,minTree,centers_ST,nodesByClusterID_ST,sum(LVCostDict_ST.values())*LV


     #### bu asagida ne yaptigimi bilmiyor olabilirim!!!

    minNodesByClusterID_ST=copy.deepcopy(minNodesByClusterID)
    minCenters_ST=copy.deepcopy(minCenters)
    minLVCostDict_ST=minLVCostDict.copy()

        
    
    #fixes tempCenterX , tempCenterY and maxDist of the first merge before going into while
    if minSeg_ST.getWeight()<=distFromT*2:
       maxDist=0 # can be anything less than 500
    else:
       maxDist=distFromT+10 # can be anything greater than 500
       print "NO CLUSTER POSSIBLE"

    tempCenter1,tempCenter2=minSeg_ST.getNodes()
       
    tempCenterX=(tempCenter1.getWeight()*tempCenter1.getX()
        +tempCenter2.getWeight()*tempCenter2.getX())/(tempCenter2.getWeight()+tempCenter1.getWeight())
    tempCenterY=(tempCenter1.getWeight()*tempCenter1.getY()
        +tempCenter2.getWeight()*tempCenter2.getY())/(tempCenter2.getWeight()+tempCenter1.getWeight())

    
    i=len(minCenters)
    while(maxDist<=distFromT):
        i-=1
    
        if i%20==0:
            print i
        center1,center2=minSeg_ST.getNodes()

	
        weight=center2.getWeight()+center1.getWeight()
        baseClusterID=min(clusterByNode_ST[center1],clusterByNode_ST[center2])
        
        mergingClusterID=max(clusterByNode_ST[center1],clusterByNode_ST[center2])
        
        nodesByClusterID_ST[baseClusterID].extend(nodesByClusterID_ST.pop(mergingClusterID))
    
        centers_ST[baseClusterID].setXY(tempCenterX,tempCenterY)
        centers_ST[baseClusterID].setWeight(weight)
    
        del centers_ST[mergingClusterID] 
   
        for node in nodesByClusterID_ST[baseClusterID]:
            clusterByNode_ST[node]=baseClusterID

        segments_ST=generateSegments(centers_ST,sr) # generate segments for new graph
        nodeDict=buildAssocDict(segments_ST)
 #       newTree,segments=kruskalsAlg(segments,centers) #returns sorted segments
        newTree=primsAlg(segments_ST, len(centers_ST), 0, nodeDict) # 0 is the starting node of prims.
        TotalMVCost_ST=newTree.getTotalEdgeWeight()*MV
        TotalTransformerCost_ST=len(centers_ST)*TCost
        #del LVCostDict_ST[mergingClusterID]
        gc.collect()
        #segmentsCMST, LVCostDict_ST[baseClusterID] =CMST_dfs.CMST(nodesByClusterID_ST[baseClusterID],maxLVLenghtInCluster,centers_ST[baseClusterID])

  #      newTotalCost_ST=TotalMVCost_ST+TotalTransformerCost_ST+(sum(LVCostDict_ST.values()))*LV
 #       newTotalCost_ST=TotalMVCost_ST+TotalTransformerCost_ST+sumLVCostAtEachStep[i]

       
        newTotalCost_ST=TotalMVCost_ST+TotalTransformerCost_ST+sumLVCostAtEachStep[len(centers_ST)]
            

 

        
        #outFile.write("%(newTotalCost_ST)f\n" %vars())
        #outFile.write("%i %f %f %f\n" %(i,sumLVCostAtEachStep[i], TotalMVCost_ST, minTotalCost_ST))
        if(newTotalCost_ST<=minTotalCost_ST):
            minNodesByClusterID_ST=copy.deepcopy(nodesByClusterID_ST)
            minTree=copy.deepcopy(newTree)
            minCenters_ST=copy.deepcopy(centers_ST)
            #minLVCostDict_ST=LVCostDict_ST.copy()
            minLVCostSum_ST=sumLVCostAtEachStep[len(centers_ST)]
            minTotalCost_ST=newTotalCost_ST
       
        # Calculate maxDist below for next graph and continue if it is less than 500

        try: # to check if there is a segment on the graph or there is only one cluster
            minSeg_ST=min(segments_ST,key=lambda obj:obj.getWeight()) 
	    maxDist,tempCenterX,tempCenterY=maxTempInClusterDist(minSeg_ST,clusterByNode_ST,nodesByClusterID_ST)
	    if maxDist>distFromT:
		segments_ST.sort(key=lambda obj:obj.getWeight())
	
		for seg in segments_ST:
		    if seg.getWeight()>distFromT*2:
			break #break from for loop
		    else:
		        maxDist,tempCenterX,tempCenterY=maxTempInClusterDist(seg,clusterByNode_ST,nodesByClusterID_ST)    
                        if maxDist <=distFromT:
			    minSeg_ST=seg
			    break #break from for loop
			    
        except:
            break
        
    outFile.close()
    '''
    #return minTotalCost_ST,minTree,minCenters_ST,minNodesByClusterID_ST,minLVCostSum_ST
    return centers,nodesByClusterID


def addLVSeg(tree,centers,nodesByClusterID):#single points line from the root
    
    SegID=1000000

    for centerID in centers.keys():
        try:
            netID=tree.getNetID(centers[centerID])
        except:
            netID=0
            tree._nodesByNetID[0]=[]
            tree._network[netID]=[]
        for node in nodesByClusterID[centerID]:
            length=((node.getX()-centers[centerID].getX())**2+
                        (node.getY()-centers[centerID].getY())**2)**(.5)
            newSeg=network.Seg(SegID, node, centers[centerID], length)
            tree._netIDByNode[node] = netID
            tree._nodesByNetID[netID].append(node)
            tree._network[netID].append(newSeg)
   
    return tree


def writeLVDictToText(statsFile, Dict):
    'Writes LVCostDict to a text file for batchPrimsForTransformers.py.'
    outFile = open(statsFile,"w")
    for key in Dict.keys():
            LVCost=Dict[key]*10
            outFile.write("%(key)i %(LVCost)f\n" %vars())
    outFile.close()
    return 0

def writeNumNodesinClusterToText(statsFile, Dict):
    'Writes LVCostDict to a text file for batchPrimsForTransformers.py.'
    outFile = open(statsFile,"w")
    for key in Dict.keys():
            NumNodes=len(Dict[key])
            outFile.write("%(key)i %(NumNodes)f\n" %vars())
    outFile.close()
    return 0


def writeCenterSizeToText(statsFile, Dict):
   
    outFile = open(statsFile,"w")
    for key in Dict.keys():
            size=Dict[key].getWeight()
            outFile.write("%(size)i \n" %vars())
    outFile.close()
    return 0


def get_debug_subgrids(txtpath, batch_num):
    #grid_files = []
    ward_files = os.listdir(txtpath)
    grid_files = sorted(ward_files)
    print("Number of grid files found", len(grid_files))
    mygrids= grid_files[int(batch_num):int(batch_num+1.0)]
    mygrids = [os.path.join(txtpath,m) for m in mygrids]
    return mygrids


def get_scale_and_subgrids(grid,start,stop=None):
    #allowed_subgrid = np.arange(start,stop+1)
    subgrids = []
    for root, dirs, files in os.walk(grid):
        for name in files:
            if 'MV' not in name and 'FinalGrid' not in name:
                if 'scale' in root and '.shp' in name:
                    valid = int(name[:-4].split('_')[-1])
                    if valid ==start:
                        subgrids.append(os.path.join(root, name))
                        return subgrids
   
def main(txtpath, batch_number, start_sub, stop_sub=None):
    my_grids = get_debug_subgrids(txtpath,batch_number)
    print(my_grids)
    #Cost parameters:
    MV =25# Cost of MV per meter
    LV = 10 # Cost of LV per meter
    TCost=2000 # Transformer Cost
    distFromT=500 #Dmax, direct distance from transformers
    maxLVLenghtInCluster=distFromT #Lmax
    #if argv is None:
    #argv = sys.argv
    for grid in my_grids:
        files = get_scale_and_subgrids(grid,start_sub,stop_sub)
        for file in files:
            if file.endswith('.shp'):
                inputShapeFile = file
                print(file)
                outputDir = file[:-4]
                searchRadius = 1000 # meters

                startTime = time.time()
                print "Generating Dictionaries"
       
                nodesByClusterID,clusterByNode,nodes,centers,LVCostDict=generateDictsFromShp(inputShapeFile,outputDir)
       
                print "Run function starts..."
                timeBeforeRun=time.time()
                centers,nodesByClusterID=run(centers,nodesByClusterID,clusterByNode,LVCostDict,searchRadius,MV,LV,TCost,distFromT,maxLVLenghtInCluster,outputDir)
        
                statsFile1= outputDir + os.sep + "NumberofNodesinEachCluster.txt" 
	        #statsFile2= outputDir + os.sep + "CenterSize.txt" 
                writeNumNodesinClusterToText(statsFile1, nodesByClusterID)
                #writeCenterSizeToText(statsFile2, centers)
                #batchPrimsforTransformers.batchPrims(tree,centers,LVCostDict,outputDir)
       
                #print "LVLength_run", sum(LVCostDict.values())
                #MVLength=tree.getTotalEdgeWeight()
        
                #MVCost=MVLength*MV
                numTransformers=len(centers)
                print "Number of Transformers:", numTransformers

                print "Number of Nodes in Each Cluster:"
                clusters_in_nodes = []
                for id in nodesByClusterID.keys():
                    clusters_in_nodes.append(len(nodesByClusterID[id]))
                    print len(nodesByClusterID[id])
                clusters_in_nodes = np.array(clusters_in_nodes)

                with open(outputDir + 'clusteringOutput.txt', 'w') as dst:
                    dst.write("NumberofClusters:" + str(len(clusters_in_nodes)) + "\n")
                    dst.write("AverageNumberofStructsperCluster:"+ str(np.mean(clusters_in_nodes))+ "\n")
                    dst.write("StdofNumberofStructsperCluster:" + str(np.std(clusters_in_nodes)) + "\n")
                    dst.write("MaxNumberofStructsperCluster:" + str(np.max(clusters_in_nodes)) + "\n")
                    dst.write("TotalRuntime:" + str(time.time()-startTime) + "\n")
                print "Total Running Time:",time.time()-startTime
        
         

if __name__ == "__main__":
    current_batch = float(sys.argv[1])
    start_subgrid =float(sys.argv[2])
    #stop_subgrid = float(sys.argv[3])
    #grid_paths = '../../viable_locs/4000band/'
    #grid_paths = '../../selected_wards_shapefiles/'
    grid_paths = '../../selected_wards_gypsum/'
    main(grid_paths, current_batch,start_subgrid)
    #sys.exit(main())
