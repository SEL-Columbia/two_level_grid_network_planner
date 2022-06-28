#!/usr/bin/env python
"""
A heuristic for Capaciated Minimum Spanning Tree

Ayse Selin Kocaman
"""

import os
import sys
import copy
import gc
import collections

import network

class T:
    def __init__(self,value1=None,value2=None,value3=None,value4=None):
        self._removedSegmentNode= value1
        self._addedSegmentNode1=value2
        self._addedSegmentNode2=value3
        self._value=value4

    def __repr__(self):
        removedSegment = self.getRemovedSegmentNode()
        addedSegmentNode1 = self.getAddedSegmentNode1()
        addedSegmentNode2 = self.getAddedSegmentNode2()
        value=self.getValue()
        return "T ( %(removedSegment)s), %(addedSegmentNode1)s, %(addedSegmentNode2)s, %(value)s)" %vars()


    def __str__(self):
        removedSegment = self.getRemovedSegmentNode()
        addedSegmentNode1 = self.getAddedSegmentNode1()
        addedSegmentNode2 = self.getAddedSegmentNode2()
        value=self.getValue()
        return "T ( %(removedSegment)s), %(addedSegmentNode1)s, %(addedSegmentNode2)s, %(value)s)" %vars()

    def getRemovedSegmentNode(self):
        return self._removedSegmentNode
    def getAddedSegmentNode1(self):
        return self._addedSegmentNode1
    def getAddedSegmentNode2(self):
        return self._addedSegmentNode2
    def getValue(self):
        return self._value

def buildAssocDict(segments):
    'Builds dictionary with nodeID key where values are segs from/to that node'
    segList = {}
    for seg in segments:
        node1, node2 = seg.getNodes()
        for nodeID in [node1.getID(), node2.getID()]:
            #if segList.has_key(nodeID):
            if nodeID in segList.keys():
                segList[nodeID].append(seg)
            else:
                segList[nodeID] = [seg]
    return segList    


def dfs(root, visited = None, 
        preorder_process  = lambda x: None):
    """
    Given a starting vertex, root, do a depth-first search.
    """
    to_visit = []  # a list can be used as a stack in Python
    if visited is None: visited = set()
 
    to_visit.append(root) # Start with root
    while len(to_visit) != 0:
        v = to_visit.pop()
        if v not in visited:
            visited.add(v)
            preorder_process(v)
            to_visit.extend(v.neighbors)


def depthFirstSet(node1,node2,root,segList,distDict):
    """
    Given a starting vertex, root, do a depth-first search and set the weights.
    """
    to_visit = []  # a list can be used as a stack in Python
    visited=[]
    node1.setWeight(node2.getWeight()+distDict[(node1.getID(),node2.getID())])
    to_visit.append(node1) # Start with root
    while len(to_visit)!= 0:
        v = to_visit.pop()
        #print v
        if v not in visited:
            visited.append(v)
            vNeighbors=[]
            #print segList[v.getID()]
            for seg in segList[v.getID()]:
                firstNode,secondNode=seg.getNodes()
                if firstNode==root or secondNode==root :
                    continue
                if firstNode==v and secondNode not in visited:
                    vNeighbors.append(secondNode)
                    secondNode.setWeight(v.getWeight()+distDict[(v.getID(),secondNode.getID())])
                if secondNode==v and firstNode not in visited:
                    vNeighbors.append(firstNode)
                    firstNode.setWeight(v.getWeight()+distDict[(v.getID(),firstNode.getID())])
            to_visit.extend(vNeighbors)

def maxDistFromNode(node1,root,treeSegments,distDict):
    """
    Given a starting vertex, do a depth-first search and find the furthest node to that vertex
    """
    segList=buildAssocDict(treeSegments.values())
    to_visit = []  # a list can be used as a stack in Python
    visited=[]
    
    tempWeightByNode={}
    to_visit.append(node1)
    tempWeightByNode[node1.getID()]=0
    while len(to_visit)!= 0:
        v = to_visit.pop()
        if v not in visited:
            visited.append(v)
            vNeighbors=[]
            for seg in segList[v.getID()]:
                firstNode,secondNode=seg.getNodes()
                if firstNode==root or secondNode==root :
                    continue
                if firstNode==v and secondNode not in visited:
                    vNeighbors.append(secondNode)
                    tempWeightByNode[secondNode.getID()]=tempWeightByNode[v.getID()]+distDict[(v.getID(),secondNode.getID())]
                if secondNode==v and firstNode not in visited:
                    vNeighbors.append(firstNode)
                    tempWeightByNode[firstNode.getID()]=tempWeightByNode[v.getID()]+distDict[(v.getID(),firstNode.getID())]
            to_visit.extend(vNeighbors)
    maxDist=max(tempWeightByNode.values())
    return maxDist


            
def CMST(households,capacity,root):
    #Connects hhs directly to the root first
    treeSegments={}
    distDict={}
    households_Copy=copy.deepcopy(households)
    SegID=10000000
    root_Copy=copy.deepcopy(root)
    newRootID=root_Copy.getID()*(-1)-100  #### not to be confused with the same nodeID
    root_Copy.setID(newRootID)
    maxTvalue=0
    branchNodeByNode={}# which branch is the node on?
    nodesByBranchNode=collections.defaultdict(list )# what are the nodes on a spesific branch?
    for node in households_Copy:
        length=((node.getX()-root_Copy.getX())**2+(node.getY()-root_Copy.getY())**2)**(.5)
        treeSegments[(node.getID(),newRootID)]=network.Seg(SegID, node, root_Copy, length)
        SegID+=1
        node.setWeight(length)
        branchNodeByNode[node]=node
        nodesByBranchNode[node].append(node)
        distDict[(newRootID,node.getID())]=length
        distDict[(node.getID(),newRootID)]=length
       
    for node1 in households_Copy:
        for node2 in households_Copy:
            if node1==node2:
                continue
            else:
                distance=((node1.getX()-node2.getX())**2+(node1.getY()-node2.getY())**2)**(.5)
                distDict[(node1.getID(), node2.getID())]=distance
                Tvalue=treeSegments[(node1.getID(),newRootID)].getWeight()-distance

                if ((node2.getWeight()+distance)<=capacity):
                    newT=T(node1,node1,node2,Tvalue)

                    
                    if (newT.getValue()>=maxTvalue):
                        maxTObject=newT
                        maxTvalue=newT.getValue()
    totalLVCost=0
    for segment in treeSegments.values():
        totalLVCost=totalLVCost+segment.getWeight() 
    while(maxTvalue>0):

        maxTvalue=0
        SegID+=1
        node1=maxTObject.getAddedSegmentNode1() # node1 and node2 of new segment
        #print "node1", node1
        node1Weigth=node1.getWeight()
        node2=maxTObject.getAddedSegmentNode2()
        #print "node2", node2
        
        #delete root segment of branch
        del treeSegments[(maxTObject.getRemovedSegmentNode().getID(),newRootID)]
       
        
        #if node1 is the first node in the branch
        if node1==branchNodeByNode[node1]:
            tempWeight=node1.getWeight() # I need this becouse node1 is updated first and it effects the others
            for node in nodesByBranchNode[branchNodeByNode[node1]]:
                node.setWeight(node.getWeight()-tempWeight+ node2.getWeight()+distDict[(node1.getID(),node2.getID())])# digerlerinin de set olmasi lazim

        #if not things get complicated and we need dfs
        else:
            segList=buildAssocDict(treeSegments.values()) # daha efficient yapmak icin sadece update edebilirim
            depthFirstSet(node1,node2,root_Copy,segList,distDict) # root  (node1) hala icinde unutma


        #tree updated after weights are set
        treeSegments[(node1.getID(), node2.getID())]=network.Seg(SegID,node1, node2,distDict[(node1.getID(),node2.getID())])
      
        
        # Update dictionaries
        nodesByBranchNode[branchNodeByNode[node2]].extend(nodesByBranchNode.pop(branchNodeByNode[node1]))
        for node in nodesByBranchNode[branchNodeByNode[node2]]:
           branchNodeByNode[node]=branchNodeByNode[node2]

        # Rebuilt TStore & select maxT object
        for node1 in households_Copy:  #
            #print "node1", node1
            for node2 in households_Copy:
                if node1==node2 or branchNodeByNode[node1]==branchNodeByNode[node2]:
                    continue
                else:
                    maxDistFromNode1=maxDistFromNode(node1,root_Copy,treeSegments,distDict)
                    #print "maaxx",maxDistFromNode1
                    if (node2.getWeight()+distDict[node1.getID(),node2.getID()]+maxDistFromNode1<=capacity): #1 2ye baslansa ne olur?
                        #print "TTTTT", Tvalue
                        Tvalue=treeSegments[(branchNodeByNode[node1].getID(),newRootID)].getWeight()-distDict[(node1.getID(),node2.getID())]
                        if Tvalue>=maxTvalue:
                            maxTvalue=Tvalue
                            maxTObject=T(branchNodeByNode[node1],node1,node2,Tvalue)
        #print maxTvalue
        #print maxTObject
    totalLVCost=0    
    for segment in treeSegments.values():
        totalLVCost=totalLVCost+segment.getWeight()
    del root_Copy
    gc.collect()
    #print treeSegments.keys()
    #print households_Copy
    return treeSegments, totalLVCost





#capacity=10
#root=network.Node(0, 452520.229406, 9445325.65803, 10)
'''root=network.Node(0,0,0,5)
households=[]
households.append(network.Node(1,0,2,1))
households.append(network.Node(2,1,2,1))
households.append(network.Node(3,-1,2,1))
households.append(network.Node(4,1,1,1))
households.append(network.Node(5,-1,1,1))
households.append(network.Node(6,2,1,1))
households.append(network.Node(7,1,4,1))
households.append(network.Node(8,-2,3,1))


households=[]
households.append(network.Node(0, 453143.895802, 9441304.74666, 1))
households.append(network.Node(1, 453141.369002, 9441342.64865, 1))
households.append(network.Node(2, 453176.112492, 9441319.90746, 1))
households.append(network.Node(3, 453062.406526, 9441304.74666, 1))
households.append(network.Node(4, 4453088.937918, 9441567.53378, 1))
households.append(network.Node(5, 452770.75615539739, 9445716.9034971781, 1))
households.append(network.Node(6, 452751.25842492335, 9445688.1985050924, 1))
households.append(network.Node(7, 452713.88777484809, 9445680.6160543524, 1))
households.append(network.Node(8, 452759.38247928751, 9445604.249943329, 1))
households.append(network.Node(9, 452783.75464238005, 9445863.1364757344, 1))

'''

#CMST(households,capacity,root)



