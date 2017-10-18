# -*- coding: utf-8 -*-
###################################################################################################
# Calculate a criticality score for various linear features arranged in a network configuration
# Yan Deng, Charles Fox and Ben Stewart, September 2017
# Purpose: determine the criticality score for each linear feature in the network dataset
###################################################################################################
import os, sys, inspect, logging
import networkx as nx
import pandas as pd
import geopandas as gpd
import numpy as np
import shapely.geometry.base
import shapely.wkt
import copy
import time
import datetime
from datetime import datetime
import warnings

warnings.filterwarnings("ignore")

module_path = os.path.join(os.path.realpath(os.path.abspath(os.path.split(inspect.getfile(inspect.currentframe() ))[0])))
module_path = os.path.join(os.path.split(module_path)[0], "PCS/Criticality")
if module_path not in sys.path:
    sys.path.append(module_path)

# Modules developed by TU Delft team
from network_lib import network_prep as net_p
from network_lib import od_prep as od_p

verbose = 1
dump = 1

def main(adminIsPoint = False):
    path = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile( inspect.currentframe() ))[0]))
    path = os.path.split(path)[0]
    dash = os.path.join(path,r'dashboard.xlsm')
    ctrl = pd.read_excel(dash, sheetname = "AGGREGATE", index_col = 0)
    district = ctrl['Weight'].loc['DISTRICT']

    logging.basicConfig(filename = os.path.join(path, 'runtime', district, "PCS_Criticality_log.log"), level=logging.INFO, format="%(asctime)s-%(levelname)s: %(message)s")
    logging.info("Starting Criticality Process")
    print "Running: Criticality Analysis on %s. Do not interrupt" % district
    # Path Settings
    outpath = os.path.join(path, 'Outputs', '%s' % district)
    runtime = os.path.join(path, r'PCS\Criticality\runtime\%s\\' % district)
    for d in [outpath, runtime]:
        if not os.path.isdir(d):
            os.mkdir(d)
    NETWORK_IN = os.path.join(path, r'runtime\%s\\' % district)
    OD_IN = os.path.join(path, 'PCS\Criticality\input', '%s' % district)
    DATA_IN = os.path.join(path, 'PCS\Criticality\Vietnam_Data_Layers')
    inAdmin = os.path.join(DATA_IN,'Poverty_Communes_2009.shp')
    inNetworkFile = os.path.join(NETWORK_IN, 'Network.csv')

    crs_in = {'init': 'epsg:4326'}   #WGS 84

    #Create folders for analysis
    for d in [outpath, runtime, OD_IN]:
        if not os.path.isdir(d):
            os.mkdir(d)
    #Error checking - Check input data
    for curFile in [dash, inNetworkFile, inAdmin,DATA_IN,OD_IN,NETWORK_IN]:
        if not os.path.exists(curFile):
            logging.error("No input found: %s" % curFile)
            raise ValueError("No input found: %s" % curFile)


    inNetwork = pd.read_csv(inNetworkFile)
    ctrldf = pd.read_excel(dash, sheetname = "CRITICALITY", index_col = 'COL_ID')
    #Inputs
    network = os.path.join(runtime,'Network.shp')

    #Network Prep
    fillvalue = inNetwork['iri_med'].mean()
    inNetwork['TC_iri_med'] = inNetwork['iri_med'].fillna(fillvalue)
    inNetwork['total_cost'] = inNetwork['length'] * (ctrldf['Base_cost_km'][0] + (ctrldf['IRI_Coeff'][0] * inNetwork['TC_iri_med']))
    ginNetwork = gpd.GeoDataFrame(inNetwork,crs = crs_in, geometry = inNetwork['Line_Geometry'].map(shapely.wkt.loads))
    ginNetwork.to_file(network, driver = 'ESRI Shapefile')
    logging.info("Successfully loaded data")
    if not adminIsPoint:
        prepareAdminCentroids(ginNetwork, inAdmin, crs_in, os.path.join(OD_IN, 'adm_centroids.shp'))
        logging.info("Created admin centroids")
    def makeOrigin(n, ctrldf):
        origindict = {
            'name': ctrldf['OName'][n],
            'file': os.path.join(path,'PCS','Criticality','input',district,'%s.shp'% ctrldf['OName'][n]),
            'scalar_column': ctrldf['OScalar'][n]
            }
        return origindict
    def makeDestination(n, ctrldf):
        destdict = {
            'name': ctrldf['DName'][n],
            'file': os.path.join(path,'PCS','Criticality','input',district,'%s.shp'% ctrldf['DName'][n]),
            'penalty': ctrldf['DPenalty'][n],
            'importance':ctrldf['DImportance'][n],
            'annual':ctrldf['DAnnual'][n],
            'scalar_column': ctrldf['DScalar'][n]
            }
        return destdict

    origin_1, origin_2, origin_3, origin_4, origin_5 = makeOrigin(0, ctrldf), makeOrigin(1, ctrldf), makeOrigin(2, ctrldf), makeOrigin(3, ctrldf), makeOrigin(4, ctrldf)
    originlist = {
        '%s' % ctrldf['OName'][0]: origin_1,
        '%s' % ctrldf['OName'][1]: origin_2,
        '%s' % ctrldf['OName'][2]: origin_3,
        '%s' % ctrldf['OName'][3]: origin_4,
        '%s' % ctrldf['OName'][4]: origin_5,
        }
    destination_1, destination_2, destination_3, destination_4, destination_5 = makeDestination(0, ctrldf), makeDestination(1, ctrldf), makeDestination(2, ctrldf), makeDestination(3, ctrldf), makeDestination(4, ctrldf)
    destinationlist = {
        '%s' % ctrldf['DName'][0]: destination_1,
        '%s' % ctrldf['DName'][1]: destination_2,
        '%s' % ctrldf['DName'][2]: destination_3,
        '%s' % ctrldf['DName'][3]: destination_4,
        '%s' % ctrldf['DName'][4]: destination_5,
        }
    logging.debug("Opened origins and destinations")
    # Prepation of network
    gdf_points, gdf_node_pos, gdf = net_p.prepare_centroids_network(origin_1['file'], network)
    # Create Networkx MultiGraph object from the GeoDataFrame
    G = net_p.gdf_to_simplified_multidigraph(gdf_node_pos, gdf, simplify=False)
    # Change the MultiGraph object to Graph object to reduce computation cost
    G_tograph = net_p.multigraph_to_graph(G)
    logging.debug('Loaded road network: number of disconnected components is: %d' % nx.number_connected_components(G_tograph))
    # Observe the properties of the Graph object
    nx.info(G_tograph)
    # Take only the largest subgraph with all connected links
    len_old = 0
    for g in nx.connected_component_subgraphs(G_tograph):
        if len(list(g.edges())) > len_old:
            G1 = g
            len_old = len(list(g.edges()))
    G_sub = G1.copy()

    nx.info(G_sub)

    # Save the simplified transport network into a GeoDataFrame
    gdf_sub = net_p.graph_to_df(G_sub)
    blank, gdf_node_pos2, gdf_new = net_p.prepare_newOD(origin_1['file'], gdf_sub)

    #Road Network Graph prep
    G2_multi = net_p.gdf_to_simplified_multidigraph(gdf_node_pos2, gdf_new, simplify=False)
    Filedump(gdf_new, 'Road_Lines', runtime)
    Filedump(gdf_node_pos2,'Road_Nodes', runtime)
    G2 = net_p.multigraph_to_graph(G2_multi)
    gdf2 = net_p.graph_to_df(G2)
    nLink = len(G2.edges())

    Outputs, cost_list, iso_list = [], [], []

    for z in ctrldf.index:
        if (((ctrldf['ComboO'][z]) != 0) & ((ctrldf['ComboD'][z]) != 0) & (pd.notnull(ctrldf['ComboO'][z])) & (pd.notnull(ctrldf['ComboO'][z]))):
            Q = int(ctrldf['ComboNumber'][z])
            logging.info('Computing | combination %s as origin and %s as destination ' % (ctrldf['ComboO'][z],ctrldf['ComboD'][z]))
            xx = calculateOD(originlist['%s' % ctrldf['ComboO'][z]], destinationlist['%s' % ctrldf['ComboD'][z]], Q, gdf_sub, G2, nLink, gdf2, runtime, ctrldf)
            Outputs.append(xx)
            cost_list.append("Social_Cost_%s" % Q)
            iso_list.append("Isolated_Trips_%s" % Q)

    Output = inNetwork.drop(["geometry",'TC_iri_med','total_cost'],axis =1)
    for o_d_calc in range(0,len(Outputs)):
        Output = Output.merge(Outputs[o_d_calc]['summary'],how = 'left', on = 'ID')

    Output['Cost_total'] = Output[cost_list].sum(axis = 1)
    Output['Iso_total']  = Output[iso_list].sum(axis = 1)
    Output['CRIT_SCORE'] = (
        ctrldf['Disrupt_Weight'][0] * Output['Cost_total'] +
        ctrldf['Isolate_Weight'][0] * Output['Iso_total']
    )
    Output['CRIT_SCORE'] = ((Output['CRIT_SCORE'] - Output['CRIT_SCORE'].min()) / (Output['CRIT_SCORE'].max() - Output['CRIT_SCORE'].min()))
    logging.info("Calculated PCS Criticality")
    FileOut(Output,'criticality_output', outpath)

def prepareAdminCentroids(ginNetwork, inAdmin, crs_in, outputFile):
    #Centroid Prep
    SingleNetworkObj = pd.DataFrame([str(ginNetwork.unary_union)], columns = ['geom'])
    gSingleNetworkObj = gpd.GeoDataFrame(SingleNetworkObj,crs = crs_in, geometry = SingleNetworkObj['geom'].map(shapely.wkt.loads))
    ginAdmin = gpd.read_file(inAdmin)
    ginAdmin = ginAdmin.to_crs(crs_in)
    #TODO - This may need to be removed
    ginAdmin = ginAdmin[['P_EName','D_EName','EN_name','Admin_EN','COMPOPULA','PROCODE04','DISTCODE04','CCode04New','geometry']]
    ginAdmin = ginAdmin.rename(columns = {'COMPOPULA': 'Pop'})
    ginAdmin['ID'] = ginAdmin.index
    SelectedAdmins = gpd.sjoin(gSingleNetworkObj, ginAdmin, how="inner",op='intersects')
    SelectedAdmins = SelectedAdmins.drop(['geom','geometry'], axis = 1)
    ginAdmin = ginAdmin.loc[ginAdmin['ID'].isin(SelectedAdmins['ID']) == True]
    ginAdmin['geometry'] = ginAdmin['geometry'].centroid
    ginAdmin.to_file(outputFile, driver = 'ESRI Shapefile')

#Utility Functions
def Filedump(df, name, runtime):
    # Dumps file to runtime folder
    if dump == 1:
        df.to_csv(os.path.join(runtime, r'%s.csv' % name), encoding = 'utf-8')

def FileOut(df, name, outpath):
    # Dumps file to outpath folder
    df.to_csv(os.path.join(outpath, r'%s.csv' % name), index = False, encoding = 'utf-8')
    try:
        df = df.drop('Line_Geometry', axis = 1)
        df.to_excel(os.path.join(outpath, r'%s.xlsx' % name), index = False, encoding = 'utf-8')
    except:
        pass
def Vprint(s):
    if verbose == 1:
        ts = time.time()
        st = datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        print '\n%s -- %s' % (st, s)

##########Criticality Script #########

def PrepSet(point_set, gdf_sub):
    '''
    Prepares a small df of a given origin / destination set, expressed as 'item : Nearest Node ID'
    '''
    Prepared_point_set, gdf_node_pos2, gdf_new = net_p.prepare_newOD(point_set['file'], gdf_sub)
    Prepared_point_set = Prepared_point_set['Node']
    return Prepared_point_set

def Pathfinder(origin, destination, G2):
    '''
    find the shortest path for each od to minimize the total travel cost;
    output:
    1) baseCost ($): total travel cost between all OD pairs;
    2) basePath the shortest path between all OD pairs;
    3) baseLength: total travel distance between all OD pairs.
    '''
    basePath = [[ [] for d in range(len(destination))] for d in range(len(origin))]
    baseCost = np.zeros((len(origin),len(destination)))
    baseLength = np.zeros((len(origin),len(destination)))
    for o in range(len(origin)):
        for d in range(len(destination)):
            basePath[o][d] = nx.dijkstra_path(G2,origin[o],destination[d], weight = 'total_cost')
            baseCost[o][d] = nx.dijkstra_path_length(G2,origin[o],destination[d], weight = 'total_cost')
            baseLength[o][d] = nx.dijkstra_path_length(G2,origin[o],destination[d], weight = 'length')
    return basePath, baseCost, baseLength

def BreakEdge(origin, destination, penalty, baseCost, name, G2, runtime, nLink, gdf2):
    '''
    fundamental criticaility analysis: remove each link by changing it to a extremely large value,
    calculate the shortest path, then re-set total cost of traversing link.
    output:
    1.) diff: a 3D matrix describing the change in cost for each 2D O-Dmatrix when a given link is broken
    2.) iso: a 3D matrix describing the cost of isolated journeys for each 2D O-D matrix when a given link is broken
    '''
    to_allNode = []
    G = copy.deepcopy(G2)
    cost_disrupt = np.zeros((nLink, len(origin),len(destination)))
    ts = time.time()
    for road in range(len(gdf2)):
        ts1 = time.time()
        w = G[gdf2['FNODE_'][road]][gdf2['TNODE_'][road]]['total_cost']
        G[gdf2['FNODE_'][road]][gdf2['TNODE_'][road]]['total_cost'] = 1e10
        for o in range(len(origin)):
            for d in range(len(destination)):
                cost_disrupt[road][o][d] = nx.dijkstra_path_length(G,origin[o],destination[d], weight = 'total_cost')
        G[gdf2['FNODE_'][road]][gdf2['TNODE_'][road]]['total_cost'] = w # Resetting traverse cost to original (w)
        logging.info("Computation completed for road: %s of %s. Compute time: %s seconds" % (road+1, len(gdf2)+1, (time.time() - ts1)))
        if ((road+1) % 10) == 0:
            logging.info("estimated time remaining for this O-D combination (hours): %s" % (((1 / ((road+1) / (len(gdf2)+1))) * (time.time() - ts))/3600))
    diff = (cost_disrupt - baseCost)
    Filedump(pd.DataFrame(cost_disrupt[1]),'cost_disrupt_%s' % name, runtime)
    # change cost of the isolated OD to penalty value
    iso = np.zeros((nLink, len(origin),len(destination)))
    for index,item in np.ndenumerate(cost_disrupt):
        if item >= 1e9:
            diff[index] = penalty
            iso[index] = 1
    return diff, iso

def GenerateDemandFunction(origin,destination,baseLength):
    '''
    1.) Output: a demand function of form:
    Trips[o,d] = (Pop[o] * Annual_Number_of_trips * x)
        where sum(x) is 1 and  x = f(norm(Scalar[d]), (exp -(Distance[o,d] * importance))
    '''
    demand = np.zeros((len(origin['P']),len(destination['P'])))
    O_DF = gpd.read_file(origin['file'])
    D_DF = gpd.read_file(destination['file'])
    for o in O_DF.index:
        destscalar = pd.Series(D_DF[destination['scalar_column']] / sum(D_DF[destination['scalar_column']]))
        distscalar = pd.Series(np.exp(-1*(baseLength[o] * (destination['importance']/10))))
        distscalar[distscalar >= 1] = 0
        weight = pd.Series(destscalar * distscalar)
        normweight = pd.Series(weight / sum(weight))
        X = pd.DataFrame({'x': normweight})
        for d in D_DF.index:
            if (o == d) & (origin['name'] == destination['name']):
                pass
            else:
                demand[o][d] = ((O_DF[origin['scalar_column']][o] * destination['annual']) * X['x'][d])
    return demand

def summarise(diff, iso, demand, origin, destination, Q, nLink, gdf2, ctrldf):
    '''
    For each link disrupted, we find the number of trips that cannot be completed due to link disruption: isolate_sumTrip
    Total social cost is defined as a link disrupted under traffic: disrupt_sumCost
    Total number of trips being disrupted: disrupt_sumTrip
    '''
    disrupt_sumCost, disrupt_sumTrip, isolate_sumTrip = np.zeros(nLink),np.zeros(nLink),np.zeros(nLink)
    for i in range((nLink)):
        disrupt_sumCost[i] = np.nansum(np.multiply(diff[i,:,:],demand))
        isolate_sumTrip[i] = np.nansum(np.multiply(iso[i,:,:],demand))
    # criticality analysis with traffic for total user cost
    criticality = np.column_stack((gdf2['ID'],disrupt_sumCost,isolate_sumTrip))
    criticality[:,0] = criticality[:,0].astype(int)
    out = pd.DataFrame({'ID':criticality[:,0],'Social_Cost_%s' % Q:criticality[:,1],'Isolated_Trips_%s' % Q:criticality[:,2]})
    a = out['Social_Cost_%s' % Q]
    b = out['Isolated_Trips_%s' % Q]
    out['Social_Cost_score_%s' % Q] = ((a - a.min()) / (a.max()- a.min()))
    out['Isolated_score_%s' % Q] = ((b - b.min()) / (b.max()- b.min()))
    out['CRIT_SCORE_%s' % Q] = (
        ctrldf['Disrupt_Weight'][0] * out['Social_Cost_score_%s' % Q] +
        ctrldf['Isolate_Weight'][0] * out['Isolated_score_%s' % Q]
        )
    out = out[['ID','Social_Cost_%s' % Q,'Isolated_Trips_%s' % Q,'CRIT_SCORE_%s' % Q]]
    return out

def calculateOD(origin, destination, Q, gdf_sub, G2, nLink, gdf2, runtime, ctrldf):
    origin['P'], destination['P'] = PrepSet(origin, gdf_sub), PrepSet(destination, gdf_sub)
    basePath, baseCost, baseLength = Pathfinder(origin['P'], destination['P'], G2)
    diff, iso = BreakEdge(origin['P'], destination['P'], destination['penalty'], baseCost, destination['name'], G2, runtime, nLink, gdf2)
    demand = GenerateDemandFunction(origin, destination, baseLength)
    summary = summarise(diff, iso, demand, origin, destination, Q, nLink, gdf2, ctrldf)
    return({
        'origin': origin['name'],
        'destination': destination['name'],
        'basePath': basePath,
        'baseCost': baseCost,
        'baseLength': baseLength,
        'diff': diff,
        'iso': iso,
        'summary': summary})

main()
