import os, sys, logging, inspect
import pandas as pd
import numpy as np

def Main(district = 'YD'):
    path = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile( inspect.currentframe() ))[0])).split("\Scripts")[0]
    dash = os.path.join(path,r'dashboard.xlsm')
    ctrl = pd.read_excel(dash, sheetname = "AGGREGATE", index_col = 0)
    district = ctrl['Weight'].loc['DISTRICT']

    #Set up logging
    logging.basicConfig(filename = os.path.join(path, 'runtime', district,"RONET_log.log"), level=logging.INFO, format="%(asctime)s-%(levelname)s: %(message)s")
    logging.info("Starting RONET input preparation process")

    # Picks up the post-manipulation road df - after PCS has been added
    roadpath = os.path.join(path, 'Outputs', district, 'PCS.csv')
    if not os.path.exists(roadpath):
        logging.error("No input found: %s" % roadpath)
        raise ValueError("No input found: %s" % roadpath)
    df = pd.read_csv(roadpath)
    df['VPROMMS_type'] = df['VPROMMS_ID'].astype(str).str[2]
    df['VPROMMS_type'].loc[df.VPROMMS_type.isin(['4','3','2']) == False] = 'unrecognised'

    # Pick up dashboard
    if not os.path.exists(dash):
        logging.error("No input found: %s" % dash)
        raise ValueError("No input found: %s" % dash)
    ctrldf = pd.read_excel(dash, sheetname = 'RONET_mirror', index_col = 0)
    ctrldf.index = ctrldf.index.map(unicode)

    # Define defaults
    for col in ctrldf.columns[0:10]:
        if col in df.columns:
            df[col] = df[col].fillna(df['VPROMMS_type'].replace(ctrldf[col]))
        else:
            df[col] = df['VPROMMS_type'].replace(ctrldf[col])

    # Define Condition Class based on iri:

    df['ConditionClass'] = ctrldf['noIRI'].iloc[3]
    df['ConditionClass'].loc[df['iri_med'] > 0] = 1
    df['ConditionClass'].loc[df['iri_med'] > ctrldf['CC_vgood'].iloc[3]] = 2
    df['ConditionClass'].loc[df['iri_med'] > ctrldf['CC_good'].iloc[3]] = 3
    df['ConditionClass'].loc[df['iri_med'] > ctrldf['CC_fair'].iloc[3]] = 4
    df['ConditionClass'].loc[df['iri_med'] > ctrldf['CC_poor'].iloc[3]] = 5
    Outpath = os.path.join(path, 'RONET', district)
    if not os.path.exists(Outpath):
        try:
            os.mkdir(Outpath)
        except:
            logging.error("not possible to construct output directory")
            sys.exit
    df.to_excel(os.path.join(Outpath,'prepared_network.xlsx'),sheet_name = 'RONET_Input', index = False, engine = "openpyxl")
    logging.info("RONET input preparation complete")
Main()
