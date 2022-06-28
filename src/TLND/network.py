#!/use/bin/env python
"""
Defines classes to allow creation of undirected networks. Separate classes
are used for the full network, the constituent line segments, and their
constituent nodes.

Alex Zvoleff, aiz2101@columbia.edu
"""

class Network():
    'Implements an undirected network made up of segments and nodes.'
    def __init__(self):
        self._network = {}
        self._netIDByNode = {}
        self._nodesByNetID = {}

    def __repr__(self):
        #TODO: Finish this
        return "__repr__ UNDEFINED"
        
    def __str__(self):
        #TODO: Finish this
        return "__str__ UNDEFINED"

    def addSeg(self, newSeg):
        'Adds a new segment to a network, taking care of netID assignment.'
        node1 = newSeg.getNode1()
        node2 = newSeg.getNode2()
        # Neither node already in network:
        if (not self.inNet(node1)) and (not self.inNet(node2)):
            try:
                newID = max(self.listNetIDs()) + 1 # fails for first seg
            except ValueError:
                # only occurs when first two nodes are added to the network
                newID = 0;
            self._netIDByNode[node1] = newID
            self._netIDByNode[node2] = newID
            self._nodesByNetID[newID] = [node1, node2]
            self._network[newID] = [newSeg]
        # Both or only one node already in network:
        else:
            # Both nodes already in network:
            if self.inNet(node1) and self.inNet(node2):
                baseNetID = min(self.getNetID(node1), self.getNetID(node2))
                mergingNetID = max(self.getNetID(node1), self.getNetID(node2))
                mergingNodes = self._nodesByNetID[mergingNetID]
                for node in mergingNodes:
                    self._netIDByNode[node] = baseNetID
                self._network[baseNetID].extend(
                        self._network.pop(mergingNetID))
                self._nodesByNetID[baseNetID].extend(
                        self._nodesByNetID.pop(mergingNetID))
            # Only one node already in network:
            else:
                if self.inNet(node1):
                    baseNode = node1
                    newNode = node2
                else:
                    baseNode = node2
                    newNode = node1
                baseNetID = self.getNetID(baseNode)
                self._nodesByNetID[baseNetID].append(newNode)
                self._netIDByNode[newNode] = baseNetID
            self._network[baseNetID].append(newSeg)
        return 0

    def inNet(self, value):
        if value in self._netIDByNode.keys():
            return True
        else:
            return False

    # def setInNetValue(self, node, value):
        # self._inNet(node)=value

    def getNetID(self, node):
        if node in self._netIDByNode.keys():
            return self._netIDByNode[node]
        else:
            raise ValueError("Node not in network")
        
    def getEdges(self):
        netIDs = self.listNetIDs()
        edgeList = []
        for netID in netIDs:
            for edge in self._network[netID]:
                edgeList.append(edge)
        return edgeList
    
    def getEdgesDict(self):
        netIDs = self.listNetIDs()
        edgesDict = {}
        for netID in netIDs:
            for edge in self._network[netID]:
                edgesDict[edge.getID()] = edge
        return edgesDict

    def getNodes(self):
        return self._netIDByNode.keys()

    def numEdges(self):
        return len(self.getEdges())

    def numNodes(self):
        return len(self._netIDByNode.keys())

    def numSubnets(self):
        return len(self._nodesByNetID.keys())


    def listNetIDs(self):
        return self._nodesByNetID.keys()

    #TODO: Below is deprecated. Use getTotalEdgeWeight() instead
    def getTotalWeight(self):
        print("WARNING: getTotalWeight() is deprecated")
        totalWeight = 0
        for edge in self.getEdges():
            totalWeight += edge.getWeight()
        return totalWeight

    def getTotalEdgeWeight(self):
        totalWeight = 0
        for edge in self.getEdges():
            totalWeight += edge.getWeight()
        return totalWeight

    def getTotalNodeWeight(self):
        totalWeight = 0
        for node in self.getNodes():
            totalWeight += node.getWeight()
        return totalWeight

    def getSubnetNodeWeight(self, netID):
        weight = 0
        for node in self._nodesByNetID[netID]:
            weight += node.getWeight()
        return weight

    def getSubnetEdgeWeight(self, netID):
        weight = 0
        for seg in self._network[netID]:
            weight += seg.getWeight()
        return weight

    def getSubnetEdges(self, netID):
        return self._network[netID]

    def getSubnetNodes(self, netID):
        return self._nodesByNetID[netID]

    def getNetworkArea(self, bufferDist):
        from shapely.geometry import asShape
        segs = self.getEdges()
        segUnion = asShape(segs[0]).buffer(bufferDist)
        for n in xrange(len(segs) - 1):
            print(n)
            segUnion = segUnion.union(asShape(segs[n + 1]).buffer(bufferDist))
        return segUnion.area

    def writeARC(self, filename, coords):
        'Write network to text file in ArcGIS generate format.'
        try:
            fid = open(filename,"w")
        except:
            print("ERROR: could not access " + filename)
            return 1
        segNum = 0 # required as an identifier for ArcGIS
        for segment in self.getEdges():
            segNum += 1
            x1, y1 = coords[segment.getNode1()]
            x2, y2 = coords[segment.getNode2()]
            fid.write(str(segNum) + "\n")
            fid.write("%(x1)E,%(y1)E\n" %vars())
            fid.write("%(x2)E,%(y2)E\n" %vars())
            fid.write("END\n")
        fid.close()

    def writeCSV(self, filename):
        'Write a network to a text file in CSV format.'
        try:
            fid = open(filename,"w")
        except:
            print("ERROR: could not access " + filename)
            return 1
        for segment in self.getEdges():
            ID = segment.getID()
            node1 = segment.getNode1()
            node2 = segment.getNode2()
            fid.write("%(node1)r,%(node2)r\n" %vars())
        fid.close()

    def readCSV(self, filename):
        'Read network from CSV format textfile.'
        net = Network()
        try:
            fid = open(filename,"r")
        except:
            print("ERROR: could not open " + filename)
            return 1
        line = fid.readline
        while line:
            line = fid.readline()
            ID, node1, node2, weight = line.split()
            net.addSeg(Seg(node1, node2, weight))
        fid.close()
        return net

class Seg:
    'A class representing undirected segs.'
    def __init__(self, ID=None, node1=None, node2=None, weight=None,demand=None):
        self._ID = ID
        self._node1 = node1
        self._node2 = node2
        self._weight = weight
        self._demand = demand

    def __lt__(self, other):
        if self.getWeight() < other.getWeight():
            return True
        else:
            return False

    def __le__(self, other):
        if self.getWeight() <= other.getWeight():
            return True
        else:
            return False

    def __eq__(self, other):
        'Two segments are defined as equal iff they have identical IDs.'
        if self.getID() == other.getID():
            return True
        else:
            return False

    def __ne__(self, other):
        if self.getID() != other.getID():
            return True
        else:
            return False
    
    def __gt__(self, other):
        if self.getWeight() > other.getWeight():
            return True
        else:
            return False

    def __ge__(self, other):
        if self.getWeight() >= other.getWeight():
            return True
        else:
            return False

    def __repr__(self):
        id = self.getID()
        node1 = self.getNode1()
        node2 = self.getNode2()
        weight = self.getWeight()
        demand = self.getDemand()
        return "Seg(%(id)r, %(node1)r, %(node2)r, %(weight)r, %(demand)r)" %vars()

    def __str__(self):
        id = self.getID()
        node1 = self.getNode1()
        node2 = self.getNode2()
        weight = self.getWeight()
        demand = self.getDemand()
        return "(%(id)s, %(node1)s, %(node2)s, %(weight)s, %(demand)s)" %vars()

    def __hash__(self):
        return hash(str(self.getID()))
    
    @property
    def __geo_interface__(self):
        'Provides an interface to allow conversion for use with Shapely.'
        x1 = self.getNode1().getX()
        y1 = self.getNode1().getY()
        x2 = self.getNode2().getX()
        y2 = self.getNode2().getY()
        return {'type': 'LineString', 'coordinates': ((x1, y1), (x2, y2))}

    def getWeight(self):
        return self._weight
    
    def getDemand(self):
        return self._demand

    def getID(self):
        return self._ID

    def getNode1(self):
        return self._node1

    def getNode2(self):
        return self._node2

    def getNodes(self):
        return self._node1, self._node2
    
    def lensegment(self):
        lenght=((self.getNode1().getX()-self.getNode2().getX())**2+(self.getNode1().getY()-self.getNode2().getY())**2)**0.5
        return lenght
    
class Node:
    'Defines a node class, with ID, x, y, weight and demand attributes'
    def __init__(self, value1=None, value2=None, value3=None, value4=None, value5=None):
        self._id = value1
        self._x = value2
        self._y = value3
        self._weight = value4
        self._demand = value5

    def __hash__(self):
        return hash(str(self.getID()))

    def __repr__(self):
        id = self.getID()
        x = self.getX()
        y = self.getY()
        weight = self.getWeight()
        demand = self.getDemand()
        return "Node(%(id)r, %(x)r, %(y)r, %(weight)r ,%(demand)r)" %vars()

    def __str__(self):
        id = self.getID()
        x = self.getX()
        y = self.getY()
        weight = self.getWeight()
        demand = self.getDemand()
        return "(%(id)s, %(x)s, %(y)s, %(weight)s, %(demand)s)" %vars()

    def __lt__(self, other):
        if self.getWeight() < other.getWeight():
            return True
        else:
            return False

    def __le__(self, other):
        if self.getWeight() <= other.getWeight():
            return True
        else:
            return False

    def __eq__(self, other):
        if self.getID() == other.getID():
            return True
        else:
            return False

    def __ne__(self, other):
        if self.getID() != other.getID():
            return True
        else:
            return False
    
    def __gt__(self, other):
        if self.getWeight() > other.getWeight():
            return True
        else:
            return False

    def __ge__(self, other):
        if self.getWeight() >= other.getWeight():
            return True
        else:
            return False

    @property
    def __geo_interface__(self):
        'Provides an interface to allow conversion for use with Shapely.'
        return {'type': 'Point', 'coordinates': (self._x, self._y)}

    def getID(self):
        return self._id

    def setID(self,newID):
        self._id=newID

    def getX(self):
        return self._x

    def setXY(self,X,Y):
        self._x=X
        self._y=Y

    def getY(self):
        return self._y

    def getWeight(self):
        return self._weight
    
    def setWeight(self,newWeight): # newWeight is NewMVMax or new distance in CMST
        self._weight=newWeight

    def getDemand(self):
        return self._demand
    
    def setDemand(self,newDemand): # newWeight is NewMVMax or new distance in CMST
        self._demand = newDemand
    
    def getCoords(self):
        return x, y
