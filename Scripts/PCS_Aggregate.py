# -*- coding: utf-8 -*-
###################################################################################################
# Calculate the aggregate PCS score for various linear features
# Charles Fox, September 2017
# Purpose: determine the overall spatial PCS score for each linear feature in the network dataset
###################################################################################################
import os, sys, inspect
import shapely
import rtree
import pandas as pd
import numpy as np
import geopandas as gdp

path = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile( inspect.currentframe() ))[0]))
path = os.path.split(path)[0]
dash = os.path.join(path, r'dashboard.xlsm')
ctrl = pd.read_excel(dash, sheetname = "AGGREGATE", index_col = 0)
district = ctrl['Weight'].loc['DISTRICT']
Outpath = os.path.join(path,'Outputs','%s' % district)
if not os.path.isdir(Outpath):
    os.mkdir(Outpath)
roadpath = os.path.join(path, 'runtime', '%s' % district,'Network.csv')
weightsdf = pd.read_excel(dash, sheetname = 'AGGREGATE', index_col = 0)
df = pd.read_csv(roadpath)

#Define component values
#Component1: Poverty
#Poverty - Vietnam poverty map. Tati has his name. African. (simple extract attributes)
poverty = pd.read_csv(os.path.join(Outpath,'Poverty_output.csv'))
poverty = poverty[['ID','POV_SCORE']]

#Component2: Disaster Risk
risk = pd.read_csv(os.path.join(Outpath,'risk_output.csv'))
risk = risk[['ID','RISK_SCORE']]

#Component3: Criticality
criticality = pd.read_csv(os.path.join(Outpath,'criticality_output.csv'))
criticality = criticality[['ID','CRIT_SCORE']]

#Component4: Access Index
#access index - 'keith megatool' KG
component4 = 1
df['ACCESS_SCORE'] = 0

#Component5: Roughness
roughness = pd.read_csv(os.path.join(Outpath,'roughness_output.csv'))
roughness = roughness[['ID','ROUGHNESS_SCORE']]

#Calculate PCS
for component in [poverty, risk, criticality, roughness]:
    df = df.merge(component, how = 'left', on = 'ID')

df['Access_weight'] = weightsdf['Weight']['PCS_ACCESS']
df['Pov_weight'] = weightsdf['Weight']['PCS_POV']
df['Risk_weight'] = weightsdf['Weight']['PCS_RISK']
df['Crit_weight'] = weightsdf['Weight']['PCS_CRIT']
df['Rough_weight'] = weightsdf['Weight']['PCS_ROUGH']
df['a'] = df['POV_SCORE'] * df['Pov_weight']
df['b'] = df['RISK_SCORE'] * df['Risk_weight']
df['c'] = df['CRIT_SCORE'] * df['Crit_weight']
df['d'] = df['ACCESS_SCORE'] * df['Access_weight']
df['e'] = df['ROUGHNESS_SCORE'] * df['Rough_weight']
df['PCS'] = df[['a','b','c','d','e']].sum(axis = 1, skipna = True)
#df cleanup
df = df.drop(['a','b','c','d','e'],axis = 1)
df.to_csv(os.path.join(Outpath,'PCS.csv'), index = False, encoding = 'utf-8')
df = df.drop('Line_Geometry', axis = 1)
df.to_excel(os.path.join(Outpath,'PCS.xlsx'), index = False, encoding = 'utf-8')
