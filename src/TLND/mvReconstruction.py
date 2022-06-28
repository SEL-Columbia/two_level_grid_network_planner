import os, fiona
import numpy as np
import network
from heapq import heappush, heappop
import pickle


# def make_args():
#     parser = argparse.ArgumentParser()
#     parser.add_Argument('search-radius', default=searchRadius)
#     flags = parser.parse_args()
#
#     return flags
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


def generateSegments(centers):
    segments = []
    all_nodes = {}
    FID = 0
    nodeWeight = 1

    for tx in centers:
        all_nodes[FID] = network.Node(FID, tx[0], tx[1], nodeWeight)
        FID += 1

    nodeCopy = all_nodes.copy()
    segID = 0
    for startNode in all_nodes.values():
        del nodeCopy[startNode.getID()]
        for endNode in nodeCopy.values():
            dist = ((startNode.getX() - endNode.getX()) ** 2 +
                    (startNode.getY() - endNode.getY()) ** 2) ** (.5)
            segments.append(network.Seg(segID, startNode, endNode, dist))
            segID += 1
    return segments


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


def retrieve_tx_locations(folder):
    tx_locations = []
    for root, dirs, files in os.walk(folder):
        for name in files:
            if 'MV.shp' in name:
                src = fiona.open(os.path.join(root, name))
                if len(src) > 0:
                    for elem in src:
                        line = elem['geometry']['coordinates']
                        for pt in line:
                            if pt not in tx_locations:
                                tx_locations.append(pt)
                    src.close()

    ## if we didn't have any tx cause too small compute the centroid of points and return
    if len(tx_locations) == 0:
        src_points = fiona.open(os.path.join(folder, folder.split('/')[-1] + '.shp'))
        x, y, counter = 0, 0, 0
        for elem in src_points:
            x += elem['geometry']['coordinates'][0]
            y += elem['geometry']['coordinates'][1]
            counter += 1
        src_points.close()
        centroid_x, centroid_y = x / float(counter), y / float(counter)
        tx_locations.append([centroid_x, centroid_y])

    return tx_locations


def retrieve_tx_locations_v2(folder):
    tx_locations = []
    for root, dirs, files in os.walk(folder):
        for name in files:
            if 'MV.shp' in name:
                src = fiona.open(os.path.join(root, name))
                if len(src) > 0:
                    for elem in src:
                        line = elem['geometry']['coordinates']
                        for pt in line:
                            if pt not in tx_locations:
                                tx_locations.append(pt)

                else:
                    shp = root + '.shp'
                    src_points = fiona.open(shp)
                    x, y, counter = 0, 0, 0
                    for elem in src_points:
                        x += elem['geometry']['coordinates'][0]
                        y += elem['geometry']['coordinates'][1]
                        counter += 1
                    src_points.close()
                    centroid_x, centroid_y = x / float(counter), y / float(counter)
                    tx_locations.append([centroid_x, centroid_y])

                src.close()
    return tx_locations


def compute_distance(point, x, y):
    return ((y - point[:, 1]) ** 2 + (x - point[:, 0]) ** 2) ** 0.5


def tx_to_edge(tx_points, bbox):
    tx_points = np.array(tx_points)
    centroid_x, centroid_y = bbox[0] + (bbox[2] - bbox[0]) / 2.0, bbox[1] + (bbox[-1] - bbox[1]) / 2.0

    closest_output = 1000.0 * np.ones((len(tx_points), len(bbox)))

    for i in range(4):
        pos = int(np.remainder(i, 2))
        closest_output[:, i] = np.abs(tx_points[:, pos] - bbox[i])
    furthest_output = compute_distance(tx_points, centroid_x, centroid_y)

    closest_to_edge_tx = tx_points[np.unravel_index(closest_output.argmin(), closest_output.shape)[0]]
    furthest_from_edge_tx = tx_points[furthest_output.argmin()]

    return closest_to_edge_tx, furthest_from_edge_tx


def is_ward_complete(folder):
    num_subgrids = len(os.listdir(folder))
    output_counts = 0
    for root, dirs, files in os.walk(folder):
        for name in files:
            if 'modelOutput.txt' in name:
                output_counts += 1
    if num_subgrids == output_counts:
        return True
    else:
        # return False
        return True


def get_num_structs(folder):
    num_structs = 0
    for root, dirs, files in os.walk(folder):
        for name in files:
            if 'subgrid_' in name and '.shp' in name:
                src = fiona.open(os.path.join(root, name))
                num_structs += len(src)
                src.close()
    return num_structs


def mv_correction_option_b(cur_paths):
    recomputed_mv = {}
    MV_per_meter_cost = 25
    for cur_path in cur_paths:
        for cur_folder in os.listdir(cur_path):
            if os.path.isdir(os.path.join(cur_path, cur_folder)) and (cur_folder not in recomputed_mv):
                # check if ward computation is done
                ward_status = is_ward_complete(os.path.join(cur_path, cur_folder))
                if ward_status:
                    num_structs = get_num_structs(os.path.join(cur_path, cur_folder))
                    tx_locations = retrieve_tx_locations_v2(os.path.join(cur_path, cur_folder))
                    if len(tx_locations) > 0:
                        segments_ST = generateSegments(tx_locations)
                        nodeDict = buildAssocDict(segments_ST)
                        tree = primsAlg(segments_ST, len(tx_locations), 0, nodeDict)
                        MVLength = tree.getTotalEdgeWeight()
                        MVCost = MVLength * MV_per_meter_cost
                        MVperstruct = MVLength / float(num_structs)
                        recomputed_mv[cur_folder] = MVperstruct

    with open('../../recomputed_mv.pck', 'w') as hdl:
        pickle.dump(recomputed_mv, hdl, protocol=pickle.HIGHEST_PROTOCOL)

if __name__ == "__main__":
    grid_paths = ['../../structures_in_missing_wards_20m/','../../wards_for_swarm/', '../../wards_for_habanero']
    mv_correction_option_b(grid_paths)
