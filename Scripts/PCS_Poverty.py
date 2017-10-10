# -*- coding: utf-8 -*-
###################################################################################################
# Calculate a poverty score for various linear features
# Charles Fox (and a little by Ben), September 2017
# Purpose: determine the poverty index for each linear feature in the network dataset
###################################################################################################
import os, sys, inspect, logging
import pandas as pd
import geopandas as gpd
import shapely.wkt

roadID = ''

def main(district="test", admin="Poverty_Communes_2009.shp", curRoadID="ID"):
    path = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile( inspect.currentframe() ))[0])).split("\scripts")[0]
    dash = os.path.join(path, r'dashboard.xlsm')
    ctrl = pd.read_excel(dash, sheetname = "AGGREGATE", index_col = 0)
    district = ctrl['Weight'].loc['DISTRICT']

    admin = r'Poverty_Communes_2009.shp'
    roadID = curRoadID

    #Set up logging
    logging.basicConfig(filename = os.path.join(path, 'runtime',district,"PCS_Poverty_log.log"), level=logging.INFO, format="%(asctime)s-%(levelname)s: %(message)s")
    logging.info("Starting PCS Poverty Process")

    povpath = os.path.join(path,'PCS','Poverty')
    povAdmin = os.path.join(povpath, admin)
    roadpath = os.path.join(path,'Runtime', '%s' % district, 'Network.csv')

    crs = {'init': 'epsg:4326'}

    for curFile in [povpath, povAdmin, roadpath, dash]:
        if not os.path.exists(curFile):
            logging.error("No input found: %s" % curFile)
            raise ValueError("No input found: %s" % curFile)

    #read CMS output (roads) in as GDF
    df = pd.read_csv(roadpath)
    if df.index.max() < 2:
        logging.error("roadfile contains fewer than 2 roads")
        raise ValueError("roadfile contains fewer than 2 roads")
    try:
        geometry = df['Line_Geometry'].map(shapely.wkt.loads)
    except:
        logging.error("column 'Line Geometry' may not exist, or does not contain WKT formatted string")
        raise ValueError("column 'Line Geometry' may not exist, or does not contain WKT formatted string")
    gdf_road = gpd.GeoDataFrame(df, crs=crs, geometry=geometry)
    #Read in Shapefile as GDF
    gdf_adm = gpd.read_file(povAdmin)
    gdf_adm = gdf_adm.to_crs(crs)
    #read in weights as dictionary
    try:
        weightsdf = pd.read_excel(dash, sheetname = "POV")
    except:
        logging.error("sheetname 'POV' could not be found in dashboard (typically hidden)")
        raise ValueError("sheetname 'POV' could not be found in dashboard (typically hidden)")
    w8s = weightsdf.to_dict(orient='records')

    for x in range(0,len(w8s)):
        logging.debug("Doing calcs for: %s" % w8s[x]['B'])
        try:
            plc = ((gdf_adm[w8s[x]['B']]-gdf_adm[w8s[x]['B']].min())/(gdf_adm[w8s[x]['B']].max()-gdf_adm[w8s[x]['B']].min()))
            if w8s[x]['D'] == False:
                gdf_adm['Povcomp_%d'% (x+1)] = (1 - plc)*float(w8s[x]['C'])
            else:
                gdf_adm['Povcomp_%d'% (x+1)] = plc*float(w8s[x]['C'])
        except:
            gdf_adm['Povcomp_%d'% (x+1)] = 0
            logging.info('factor %s either weighted at 0 percent or not present')
    gdf_adm['POV_SCORE'] = gdf_adm[[
        'Povcomp_1','Povcomp_2','Povcomp_3','Povcomp_4','Povcomp_5','Povcomp_6','Povcomp_7','Povcomp_8','Povcomp_9','Povcomp_10']].sum(axis = 1, skipna = True)

    #gdf_adm = gdf_adm[['CCode04New','Ccode02_05','DISTCODE02','P_EName','D_EName','DISTCODE04','__Commun_1','COMNAME','COMPOPULA','__Commune','POV_SCORE','geometry']]
    gdf_adm.to_file(os.path.join(path, 'PCS','Poverty','pov_layer.shp'), driver = 'ESRI Shapefile')
    gdf_adm_join = gdf_adm[['geometry','POV_SCORE']]
    gdf_adm = gdf_adm.drop('geometry', axis=1)
    gdf_adm.to_excel(os.path.join(povpath, 'Pov_all_communes.xlsx'))

    #create mini DF with just the average poverty score and VPROMMS_ID
    gdf_road_join = gdf_road[['geometry',roadID]]
    join_pov = gpd.sjoin(gdf_road_join, gdf_adm_join, how = "inner", op = 'intersects')
    join_pov = join_pov.groupby([roadID]).apply(lambda x: LINEGROUPER(x))
    gdf_road = gdf_road.merge(join_pov, how = "inner", on = roadID)
    gdf_road = gdf_road.drop('geometry',axis = 1)
    gdf_road = pd.DataFrame(gdf_road)
    gdf_road['POV_SCORE'] = ((gdf_road['POV_SCORE'] - gdf_road['POV_SCORE'].min()) / (gdf_road['POV_SCORE'].max() - gdf_road['POV_SCORE'].min()))
    outpath = os.path.join(path, 'Outputs','%s' % district)
    if not os.path.isdir(outpath):
        os.mkdir(outpath)
    gdf_road.to_csv(os.path.join(outpath,'poverty_output.csv'), index = False, encoding = 'utf-8')
    gdf_road = gdf_road.drop('Line_Geometry', axis = 1)
    gdf_road.to_excel(os.path.join(outpath,'poverty_output.xlsx'), index = False, encoding = 'utf-8')
    logging.info("Finished PCS Poverty Calculations")

#Spatial Join
def LINEGROUPER(x2):    #EDIT
    y = pd.DataFrame()
    y['POV_SCORE'] = [x2.POV_SCORE.mean()]
    y['ID'] = [x2.ID.iloc[0]]
    return y

main()
