'''
THIS SCRIPT DOES THE FOLLOWING:
    1) Picks a specific ward and gets all the structures in ward
    2) Saves them to pickle
'''

import pickle
import re
import os
import sys
import fiona
import pandas as pd
from shapely.geometry import Polygon, Point, mapping
from rtree import index
from fiona.crs import from_epsg
import numpy as np
import simone_agg_clustering as ac


def cleanNames(name):
    name = name.title()  # converts to CamelCase
    valids = re.findall(r"[\w']+", name)  # Removes all unwanted xters
    valids = ''.join(valids)
    return (valids)


def createIndex(structures):
    # initialize index
    idx = index.Index()
    for pos, row in structures.iterrows():
        try:
            left, bottom, right, top = pointBBox(row)
            idx.insert(int(pos), (left, bottom, right, top))
        except:
            print(row)
    return (idx)


def pointBBox(row):
    x, y = row['gps_x'], row['gps_y']
    return (x, y, x, y)


def jointBBox(poly):
    left, bottom, right, top = [], [], [], []
    for geom in poly.geoms:
        left.append(geom.bounds[0])
        bottom.append(geom.bounds[1])
        right.append(geom.bounds[2])
        top.append(geom.bounds[3])
    return ((min(left), min(bottom), max(right), max(top)))


def structuresInPoly(poly, idx):
    largestBB = jointBBox(poly)
    structureIdx_in_polygon = list(idx.intersection(largestBB))  ### returns index of structures in each polygon
    return (structureIdx_in_polygon)


def getPointsInLevel(structures, poly):
    df = pd.DataFrame(columns=['utm_x', 'utm_y', 'gps_x', 'gps_y'])
    for ind, row in structures.iterrows():
        bldgPoint = Point(row[['gps_x', 'gps_y']])
        for geom in poly.geoms:
            if geom.contains(bldgPoint):
                df = df.append(
                    {'utm_x': row['utm_x'], 'utm_y': row['utm_y'], 'gps_x': row['gps_x'], 'gps_y': row['gps_y']},
                    ignore_index=True)
    return (df)


def getSaveToken(levelDir, ward, county):
    return (levelDir + '/' + ward + '_' + county + '.pck')


def loggers(log_filename, data):
    with open(log_filename, 'w') as hdl:
        pickle.dump(data, hdl, protocol=pickle.HIGHEST_PROTOCOL)


def getAllStructs(structureIdx_in_polygon, structures, poly, ward, county, levelDir):
    IDXlist = [int(k) for k in structureIdx_in_polygon]  ## verify on cluster might not need list version
    try:
        curDF = structures.loc[IDXlist]
        ## checkif DF is not empty
        if not curDF.empty:
            inShapeStructs = getPointsInLevel(curDF, poly)
            print('Number of structures in ', ward, county, ':', inShapeStructs.shape[0])
            savetoken = getSaveToken(levelDir, ward, county)
            loggers(savetoken, inShapeStructs)
            print('Saved structures in ', ward)
    except Exception as LErr:
        print('Looping Error for index:', LErr)
    return inShapeStructs


def is_ward_merged(ward, county, data):
    valid = False
    for f in data:
        if ward in f and county in f:
            valid = True
    return valid


# def make_save_token(ward,county):
#     token = ward + '_' + county
#     outputDir = '../../all_wards_20m/'+token
#     if not os.path.exists(outputDir):
#         os.mkdir(outputDir)
#     savepath = outputDir+ '/'+token + '.shp'
#     return savepath


def make_savetoken(srcpath, idx):
    subgrid = 'subgrid' + '_' + str(idx)

    if not os.path.exists(os.path.join(srcpath, subgrid)):
        os.mkdir(os.path.join(srcpath, subgrid))

    savetoken = os.path.join(srcpath, subgrid) + '/' + subgrid + '.shp'
    return savetoken


def make_save_dir(ward, county):
    token = ward + '_' + county
    outputDir = '../../structures_in_missing_wards_20m/' + token
    if not os.path.exists(outputDir):
        os.mkdir(outputDir)
    savepath = outputDir + '/'
    return savepath


def make_shp(savepath, structs):
    myschema = {'geometry': 'Point',
                'properties': {
                    'utm_x': 'float',
                    'utm_y': 'float'}}

    with fiona.open(savepath, 'w', crs=from_epsg(4326), driver='ESRI Shapefile', schema=myschema) as output:
        for idx, row in structs.iterrows():
            point = Point(float(row['utm_x']), float(row['utm_y']))
            prop = {'utm_x': row['utm_x'], 'utm_y': row['utm_y']}
            output.write({'geometry': mapping(point), 'properties': prop})


def get_extent(structures):
    '''
    Get the structure extent and add a 10 m buffer
    :param structures:
    :return:
    '''
    minx = min(structures['utm_x']) - 10
    maxx = max(structures['utm_x']) + 10
    miny = min(structures['utm_y']) - 10
    maxy = max(structures['utm_y']) + 10
    return minx, miny, maxx, maxy


def make_df(src, bbox):
    df = []
    count = 0
    for idx, elem in src.iterrows():
        if int(np.floor(elem['utm_x'])) >= bbox[0] and int(np.floor(elem['utm_x'] < bbox[2])):
            if int(np.floor(elem['utm_y'])) >= bbox[1] and int(np.floor(elem['utm_y'] < bbox[-1])):
                df.append([elem['utm_x'], elem['utm_y']])
                count += 1
    df = pd.DataFrame(df, columns=['utm_x', 'utm_y'])
    return df


def recursive_split(start_bbox, src,numStructures,clusterSize,radius):
    stack, valid = [start_bbox], []
    while stack:
        bbox = stack.pop()
        # print(bbox)
        cur_df = make_df(src, bbox)
        # max_nodes = ac.get_max_nodes_in_cluster(cur_df)
        print('Number of Structures: ', cur_df.shape[0])
        if cur_df.shape[0] < numStructures:
            if cur_df.shape[0] <= clusterSize:
                valid.append(bbox)
            else:
                max_nodes = ac.get_max_nodes_in_cluster(cur_df)
                if max_nodes <= clusterSize:
                    valid.append(bbox)
                else:
                    y_step = (bbox[-1] - bbox[1]) / 2.0
                    x_step = (bbox[2] - bbox[0]) / 2.0
                    cur_radius = int(((y_step ** 2) + (x_step ** 2)) ** 0.5)
                    if cur_radius <= radius:
                        valid.append(bbox)
                    else:
                        for x in np.arange(bbox[0], bbox[2], x_step):
                            for y in np.arange(bbox[1], bbox[-1], y_step):
                                sub_bbox = x, y, x + x_step, y + y_step
                                sub_bbox = [int(np.floor(b)) for b in sub_bbox]
                                stack.append(sub_bbox)
        else:
            y_step = (bbox[-1] - bbox[1]) / 2.0
            x_step = (bbox[2] - bbox[0]) / 2.0
            cur_radius = int(((y_step ** 2) + (x_step ** 2)) ** 0.5)
            if cur_radius <= 700:
                valid.append(bbox)
            else:
                for x in np.arange(bbox[0], bbox[2], x_step):
                    for y in np.arange(bbox[1], bbox[-1], y_step):
                        sub_bbox = x, y, x + x_step, y + y_step
                        sub_bbox = [int(np.floor(b)) for b in sub_bbox]
                        stack.append(sub_bbox)
        print('Stack size: ', len(stack))
    return valid


def split_by_max_nodes(start_bbox, structures):
    stack, valid = [start_bbox], []
    while stack:
        bbox = stack.pop()
        print(bbox)
        cur_df = make_df(structures, bbox)
        max_nodes = ac.get_max_nodes_in_cluster(cur_df)
        print('Max Nodes: ', max_nodes)
        if max_nodes < 300:
            valid.append(bbox)
        else:
            y_step = (bbox[-1] - bbox[1]) / 2.0
            x_step = (bbox[2] - bbox[0]) / 2.0
            if x_step < 500 or y_step < 500:
                valid.append(bbox)
            else:
                for x in np.arange(bbox[0], bbox[2], x_step):
                    for y in np.arange(bbox[1], bbox[-1], y_step):
                        sub_bbox = x, y, x + x_step, y + y_step
                        sub_bbox = [int(np.floor(b)) for b in sub_bbox]
                        stack.append(sub_bbox)
        print('Stack size: ', len(stack))
    return valid


def save_splits(splits, structures, ward_dir):
    print('Number of Acceptable extents:', len(splits))
    rand_ext = np.random.randint(low=0, high=len(splits))
    print('Sample acceptable extent:', splits[rand_ext])
    print('Starting the make shapefile process')
    idx = 0
    for pos, split in enumerate(splits):
        cur_df = make_df(structures, split)
        if not cur_df.empty:
            savetoken = make_savetoken(ward_dir, idx)
        make_shp(savetoken, cur_df)
        idx += 1


def main(ward, county,numStructures,clusterSize,radius):
    # levelDir = '../../structures_in_wards/'
    merged_csv_folder = '../../ward_merged_points_ratio2_20meters/'
    # merged_csv_folder = '../../missingward_merged_points_ratio2_20meters/'
    merged_csv_files = os.listdir(merged_csv_folder)
    print('Retrieving structures in bounding rectangles for:', ward, county)
    is_merged = is_ward_merged(ward, county, merged_csv_files)
    if is_merged:
        correct_file = [f for f in merged_csv_files if (ward in f) and (county in f)][0]
        all_structs = pd.read_csv(os.path.join(merged_csv_folder, correct_file), index_col=False)
        all_structs.columns = ['utm_x', 'utm_y', 'num_structures_in_cluster', 'area', 'population']
    else:
        # correct_file = "../../wards_csv_20m/" + ward + "_" + county + "_.pck"
        unmerged_files = os.listdir("../../missingwards_csv_20m/")
        for cur_file in unmerged_files:
            tokens = cur_file.split('/')[-1].split('_')
            if (tokens[0].lower() == ward.lower()) and (tokens[1].lower() == county.lower()):
                correct_file = os.path.join("../../missingwards_csv_20m/",cur_file)

        all_structs = pd.read_pickle(correct_file)
    if not all_structs.empty:
        try:
            save_dir = make_save_dir(ward, county)
            bbox = get_extent(all_structs)
            bbox = [int(np.floor(b)) for b in bbox]
            valid_splits = recursive_split(bbox, all_structs,numStructures,clusterSize,radius)
            save_splits(valid_splits, all_structs, save_dir)
        except Exception as err:
            print(err)
    else:
        print('No structures found for ', ward, county)


if __name__ == '__main__':
    ward = "YimboWest"
    county = "Siaya"
    numStructures = 3000
    clusterSize = 300
    radius = 700
    main(ward, county)