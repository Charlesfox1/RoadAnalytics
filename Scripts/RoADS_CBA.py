# -*- coding: utf-8 -*-
"""
Created on Wed 20 Sep
@author: Charles Fox
"""
import pandas as pd
import numpy as np
import os
from copy import copy
# Inputs
path = r'C:\Users\charl\Documents\Vietnam\Ronet\CFOX\\'
f = 'RONET.xlsx'
tables = os.path.join(path, f)

# Construct DFs
print 'reading in'
Control = pd.read_excel(tables, sheetname = 'Control', header = 0, index_col = 0)
VehicleFleet = pd.read_excel(tables, sheetname = 'VehicleFleet', header = 0, index_col = 0)
TrafficGrowth = pd.read_excel(tables, sheetname = 'TrafficGrowth', header = 0, index_col = 0)
RoadWorks = pd.read_excel(tables, sheetname = 'RoadWorks', header = [0,1], index_col = 0)
RecurrentMaintenance = pd.read_excel(tables, sheetname = 'RecurrentMaintenance', header = [0,1], index_col = 0)
Width = pd.read_excel(tables, sheetname = 'Width', header = 0, index_col = 0)
MCoeff = pd.read_excel(tables, sheetname = 'MCoeff', header = 0, index_col = 0)
SurfaceType = pd.read_excel(tables, sheetname = 'SurfaceType', header = 0, index_col = 0).squeeze()
ConditionData = pd.read_excel(tables, sheetname = 'ConditionData', header = 0, index_col = [0,1])
TrafficLevels = pd.read_excel(tables, sheetname = 'TrafficLevels', header = [0,1], index_col = 0)
VOC = pd.read_excel(tables, sheetname = 'VOC', header = 0, index_col = [0,1,2])
Speed = pd.read_excel(tables, sheetname = 'Speed', header = 0, index_col = [0,1,2])
RoadDeterioration = pd.read_excel(tables, sheetname = 'RoadDeterioration', header = [0,1], index_col = [0])
WorkEvaluated = pd.read_excel(tables, sheetname = 'WorkEvaluated', header = [0,1], index_col = [0,1,2])
WorkEvaluated.index.set_names(['SurfaceType','RoadClass','ConditionClass'],inplace = True)
VehicleTypes = ['Motor_cycle','Small_Car','Medium_Car','Delivery_Vehicle','4_Wheel_Drive','Light_Truck','Medium_Truck','Heavy_Truck','Articulated_Truck','Small_Bus','Medium_Bus','Large_Bus']

# Import & Name Constants
EcoFactor = Control['Value'].loc['Economic_Cost_Factor_(#)']
DiscountRate = Control['Value'].loc['Discount_Rate_(%)']

# Add defaults
Roads = pd.read_excel(tables, sheetname = 'Roads', header = 0, index_col = 0)
Roads['ID'] = Roads.index
Roads['SurfaceType'] = Roads['SurfaceType'].fillna(Roads.RoadType.replace(SurfaceType))    # Fills Surface type
Roads['Width'] = Roads['Width'].fillna(Roads.Lanes.replace(Width['Default Carriageway Width (m)']))   # Fills Width
Roads = Roads.set_index(['SurfaceType','ConditionClass']).groupby(level=[0,1]).fillna(ConditionData).reset_index() # Fills Roughness AND Pavement Age
foo = pd.DataFrame(TrafficLevels['Pavement Condition Class'].stack(), columns = ['StructuralNo']) # Fills StructuralNo for Surface Type < 4
foo.index.set_names(['TrafficLevel','ConditionClass'],inplace = True)
Roads.loc[Roads.SurfaceType < 4] = Roads.loc[Roads.SurfaceType < 4].set_index(['TrafficLevel','ConditionClass']).groupby(level=[0,1]).fillna(foo).reset_index()
Roads = Roads.sort_values(by = 'ID', ascending = 1)

AADTVehicleTypes = []
for vehicle in VehicleTypes:
    Roads['AADT_%s' % vehicle] = Roads['AADT_%s' % vehicle].fillna(Roads.TrafficLevel.replace((TrafficLevels['Default AADT (%)','%s' % vehicle]) * TrafficLevels['Default Total Traffic','(veh/day)']))
    AADTVehicleTypes.append('AADT_%s' % vehicle)

Roads['AADT_Total'] = Roads[AADTVehicleTypes].sum(axis = 1)
RoadsDict = Roads.to_dict(orient = 'index')

# Alternatives Frame
def makeAlternativeFrame(road):
    a = WorkEvaluated['BaseCase','RoadWorkNumber'].loc[road['SurfaceType'],road['RoadClass'],road['ConditionClass']]
    b = WorkEvaluated['AlternativesRoadWorks','Alternative1'].loc[road['SurfaceType'],road['RoadClass'],road['ConditionClass']]
    c = WorkEvaluated['AlternativesRoadWorks','Alternative2'].loc[road['SurfaceType'],road['RoadClass'],road['ConditionClass']]
    Work = pd.Series([a,b,b,b,b,b,b,c,c,c,c,c,c], index=[range(0,13)])
    Year = pd.Series(np.nan, index=[range(0,13)])
    Alternatives = pd.DataFrame.from_dict({'Work':Work,'Year':Year}, orient = 'columns')
    Alternatives['Year'].loc[0] = WorkEvaluated['BaseCase','RoadWorkYear'].loc[road['SurfaceType'],road['RoadClass'],road['ConditionClass']]
    if pd.notnull(Alternatives['Work'].loc[1]):
        Alternatives['Year'].loc[1:6] = range(1,7)
    if pd.notnull(Alternatives['Work'].loc[7]):
        Alternatives['Year'].loc[7:13] = range(1,7)
    Alternatives.index.set_names('Alternative')
    if Alternatives['Work'].count() == 1:
        Alternatives.loc[1] = Alternatives.loc[0] # Question for Rodrigo - when only 2, we are evaluating the same thing?  why?
    return Alternatives

# Annual Traffic
def makeAADTGrowthFrame(road):
    AADTFrame = pd.DataFrame(index = range(1,21))
    for vehicle in AADTVehicleTypes:
        g = float(TrafficGrowth[vehicle].loc[road['TrafficGrowth']])
        AADTFrame['time'] = range(1,21)
        AADTFrame['%s' % vehicle] = road[vehicle] * (1+g) ** (AADTFrame['time']-1)
    AADTFrame['Total'] = AADTFrame[AADTVehicleTypes].sum(axis = 1)
    return AADTFrame

# ESA loading
def makeESAframe(AADT, VehicleFleet, road):
    ESA = pd.DataFrame(AADT / 1000000 * 365 / road['Width'])
    Axles = pd.DataFrame({"Axles": VehicleFleet['Equivalent_Standard_Axles'], "VehicleTypes":AADTVehicleTypes}).set_index("VehicleTypes")
    ESA = ESA.drop(['time','Total'], axis = 1).transpose()
    ESA['Axles'] = Axles['Axles']
    ESA = ESA.multiply(ESA['Axles'], axis = 0).drop(['Axles'], axis = 1).transpose()
    ESAtotal = ESA.sum(axis = 1)
    return ESA, ESAtotal

# Roughness Progression Function
def makeRoughness(road_a_,a,y,dYearRoughness,iYearAge,ESAtotal):
    if road_a_['SurfaceType'] == 1 or road_a_['SurfaceType'] == 6 or road_a_['SurfaceType'] == 7:
        dYearRoughness = min((dYearRoughness * (1 + RoadDeterioration['Annual Percentage Increase','%'].loc[road_a_['SurfaceType']])),16)
    elif road_a_['SurfaceType'] == 2 or road_a_['SurfaceType'] == 3:
        case = RoadDeterioration['Roughness Progression Type','(1: Percentage, 2: Traffic & Climate: 3: Climate)'].loc[road_a_['SurfaceType']]
        m = MCoeff[road_a_['Moisture']].loc[road_a_['Temperature']]
        if case == 1: # Constant Increase
            dYearRoughness = min((dYearRoughness * (1 + RoadDeterioration['Annual Percentage Increase','%'].loc[road_a_['SurfaceType']])),16)
        elif case == 2: # Simplified HDM 4: RIb = RIa + Kgp * (a0 * Exp (Kgm * m * AGE3) * [(1 + SNC * a1)] ^ -5 * YE4 + a2 * AGE3) + (Kgm *m * RIa)
            kgp = RoadDeterioration['Roughness Progression Coefficients','kgp'].loc[road_a_['SurfaceType']]
            kgm = RoadDeterioration['Roughness Progression Coefficients','kgm'].loc[road_a_['SurfaceType']]
            a0 = RoadDeterioration['Roughness Progression Coefficients','a0'].loc[road_a_['SurfaceType']]
            a1 = RoadDeterioration['Roughness Progression Coefficients','a1'].loc[road_a_['SurfaceType']]
            a2 = RoadDeterioration['Roughness Progression Coefficients','a2'].loc[road_a_['SurfaceType']]
            YE4 = ESAtotal.loc[y]
            dYearRoughness = dYearRoughness + kgp * (a0 * np.exp(kgm * m * iYearAge) * ((1 + road_a_['StructuralNo'] + a1)) ** -5 * YE4 + a2 * iYearAge) + (kgm * m * dYearRoughness) #HDM4
        else: # Climate Related Only
            dYearRoughness = min(dYearRoughness * (1 + m),16)
    elif road_a_['SurfaceType'] == 4 or road_a_['SurfaceType'] == 5:
        dYearRoughness = min((dYearRoughness * (1 + RoadDeterioration['Annual Percentage Increase','%'].loc[road_a_['SurfaceType']])),25)
    iYearAge += 1
    return dYearRoughness, iYearAge

# Main Call
Outputs, AADTlist, IRIlist, NetBenefitslist = [], [], [], []
for x in Roads.index:
    road = RoadsDict[x]
    Alternatives = makeAlternativeFrame(road)
    iNoAlternatives = Alternatives['Work'].count()
    print '\nNumber of alternatives: %s' % iNoAlternatives
    AADT = makeAADTGrowthFrame(road)
    Load = AADT.drop(['time','Total'], axis = 1)
    Load.columns = VehicleTypes
    ESA, ESAtotal = makeESAframe(AADT, VehicleFleet, road)
    trucks = ['AADT_Light_Truck','AADT_Medium_Truck','AADT_Heavy_Truck','AADT_Articulated_Truck']
    truck_pct = (AADT[trucks].loc[1].sum() / AADT['Total'].loc[1])
    utilization = (AADT['Total'].loc[1] * road['Length'] * 365) / 1000000
    dCostFactor = WorkEvaluated['UnitCostMultiplier','-'].loc[road['SurfaceType'],road['RoadClass'],road['ConditionClass']]
    alternatives = range(0,iNoAlternatives)
    years = range(1,21)
    #zf = pd.DataFrame(0, index = alternatives, columns = years)
    dCostCapitalFin = pd.DataFrame(0, index = alternatives, columns = years)
    dCostCapitalEco = pd.DataFrame(0, index = alternatives, columns = years)
    dCostRepairFin = pd.DataFrame(0, index = alternatives, columns = years)
    dCostRepairEco = pd.DataFrame(0, index = alternatives, columns = years)
    dCostRecurrentFin = pd.DataFrame(0, index = alternatives, columns = years)
    dCostRecurrentEco = pd.DataFrame(0, index = alternatives, columns = years)
    dCostVOCs = pd.DataFrame(0, index = alternatives, columns = years)
    sRoadCode = pd.DataFrame(0, index = alternatives, columns = years)
    Sol = pd.DataFrame(0, index = alternatives, columns = ['Name','Code','Class','Year','Interval','Cost','Costkm'])
    dCondIRI = pd.DataFrame(0, index = alternatives, columns = years)
    iCondAge = pd.DataFrame(0, index = alternatives, columns = years)
    dCostVOC = pd.DataFrame(0, index = alternatives, columns = years)
    dSpeed = pd.DataFrame(0, index = alternatives, columns = years)
    dCostTime = pd.DataFrame(0, index = alternatives, columns = years)
    print 'Starting road %s of %s. Starting surface type = %s' % (x, len(Roads.index),road['SurfaceType'])
    for a in alternatives:
        road_a_ = copy(road)
        iTheWork, iTheYear = Alternatives['Work'].loc[a], Alternatives['Year'].loc[a]
        iTheRepair = RoadWorks['Future Periodic Maintenace','Road Work Number'].loc[iTheWork]
        interval = RoadWorks['Future Periodic Maintenace','Interval (years)'].loc[iTheWork]
        iTheRepairY1, iTheRepairY2, iTheRepairY3, iTheRepairY4 = iTheYear + interval,iTheYear + interval*2, iTheYear + interval*3, iTheYear + interval*4
        RepYears = [iTheRepairY1, iTheRepairY2, iTheRepairY3, iTheRepairY4]
        print "work, year: %s, %s; Repair: %s, Repair interval: %s)" % (iTheWork, iTheYear, iTheRepair, interval)
        dYearRoughness = road_a_['Roughness']
        iYearAge = road_a_['PavementAge']
        for y in years:
            # Capital Road Works
            if y == iTheYear:
                if RoadWorks[('After Road Works Condition','Lanes Class (#)')].loc[iTheWork] > 0:
                    road_a_['Lanes'] = RoadWorks[('After Road Works Condition','Lanes Class (#)')].loc[iTheWork]
                    road_a_['Width'] = RoadWorks[('After Road Works Condition','Width (m)')].loc[iTheWork]
                    road_a_['SurfaceType'] = RoadWorks[('After Road Works Condition','Surface (1 to 7)')].loc[iTheWork]
                if RoadWorks[('Bituminous Works Characteristics','Periodic M. Thickness (mm)')].loc[iTheWork] > 0:
                    road_a_['StructuralNo'] = (
                        RoadWorks[('Bituminous Works Characteristics','Periodic M. Thickness (mm)')].loc[iTheWork] *
                        RoadWorks[('Bituminous Works Characteristics','Periodic M. Strength (mm)')].loc[iTheWork] * 0.0393701)
                if RoadWorks[('Bituminous Works Characteristics','Rehab/Upgra Structural No (#)')].loc[iTheWork] > 0:
                    road_a_['StructuralNo'] = RoadWorks[('Bituminous Works Characteristics','Rehab/Upgra Structural No (#)')].loc[iTheWork]
                CapCostVar = RoadWorks.iloc[int(iTheWork)-1,road_a_['Terrain']+4] * int(road_a_['Length']) * int(road_a_['Width']) * dCostFactor / 1000
                dCostCapitalFin[y].loc[a] = CapCostVar
                dCostCapitalEco[y].loc[a] = CapCostVar * EcoFactor
                sRoadCode[y].loc[a] = RoadWorks['Work Code','-'].loc[iTheWork]
                Sol['Name'].loc[a] = RoadWorks['Work Name','-'].loc[iTheWork]
                Sol['Code'].loc[a] = RoadWorks['Work Code','-'].loc[iTheWork]
                Sol['Class'].loc[a] = RoadWorks['Work Class','-'].loc[iTheWork]
                Sol['Year'].loc[a] = y
                Sol['Interval'].loc[a] = interval
                Sol['Cost'].loc[a] = CapCostVar
                Sol['Costkm'].loc[a] = CapCostVar / int(road_a_['Length'])
            # Repair Road Works
            if y in RepYears:
                sRoadCode[y].loc[a] = RoadWorks['Work Code','-'].loc[iTheRepair]
                if RoadWorks[('After Road Works Condition','Lanes Class (#)')].loc[iTheRepair] > 0:
                    road_a_['Lanes'] = RoadWorks[('After Road Works Condition','Lanes Class (#)')].loc[iTheRepair]
                    road_a_['Width'] = RoadWorks[('After Road Works Condition','Width (m)')].loc[iTheRepair]
                    road_a_['SurfaceType'] = RoadWorks[('After Road Works Condition','Surface (1 to 7)')].loc[iTheRepair]
                if RoadWorks[('Bituminous Works Characteristics','Rehab/Upgra Structural No (#)')].loc[iTheRepair] > 0:
                    road_a_['StructuralNo'] = RoadWorks[('Bituminous Works Characteristics','Rehab/Upgra Structural No (#)')].loc[iTheRepair]
                RepCostVar = RoadWorks.iloc[int(iTheRepair)-1,road_a_['Terrain']+4] * int(road_a_['Length']) * int(road_a_['Width']) / 1000      # WHY NO COST FACTOR HERE?     # * dCostFactor
                dCostRepairFin[y].loc[a] = RepCostVar
                dCostRepairEco[y].loc[a] = RepCostVar * EcoFactor
            dCostRecurrentFin[y].loc[a] = float(RecurrentMaintenance['Number of Lanes Class', road_a_['Lanes']].loc[road_a_['SurfaceType']] * int(road_a_['Length'])) / 1000000
            dCostRecurrentEco[y].loc[a] = float(RecurrentMaintenance['Number of Lanes Class', road_a_['Lanes']].loc[road_a_['SurfaceType']] * int(road_a_['Length']) * EcoFactor) / 1000000
            # Roughness
            if y > 1:
                dYearRoughness, iYearAge = makeRoughness(road_a_,a,y,dYearRoughness,iYearAge,ESAtotal)
            if y == iTheYear:
                dYearRoughness, iYearAge = RoadWorks[('After Road Works Condition','Roughness (IRI)')].loc[iTheWork], 1
            elif y in RepYears:
                dYearRoughness, iYearAge = RoadWorks[('After Road Works Condition','Roughness (IRI)')].loc[iTheRepair], 1
            dCondIRI[y].loc[a], iCondAge[y].loc[a] = dYearRoughness, iYearAge
            road_a_['Roughness'], road_a_['PavementAge'] = dYearRoughness, iYearAge
            # Vehicle Operating Costs
            a0 = (VOC.loc[road_a_['Lanes'],road_a_['Terrain'],'a0'])
            a1 = (VOC.loc[road_a_['Lanes'],road_a_['Terrain'],'a1'] * road_a_['Roughness'])
            a2 = (VOC.loc[road_a_['Lanes'],road_a_['Terrain'],'a2'] * (road_a_['Roughness'] ** 2))
            a3 = (VOC.loc[road_a_['Lanes'],road_a_['Terrain'],'a3'] * (road_a_['Roughness'] ** 3))
            CostVOC = Load.loc[y] * (a0 + a1 + a2 +a3)
            dCostVOC[y].loc[a] = CostVOC[VehicleTypes].sum() * road_a_['Length'] * 365 / 1000000
            # Speed
            b0 = (Speed.loc[road_a_['Lanes'],road_a_['Terrain'],'b0'])
            b1 = (Speed.loc[road_a_['Lanes'],road_a_['Terrain'],'b1'] * road_a_['Roughness'])
            b2 = (Speed.loc[road_a_['Lanes'],road_a_['Terrain'],'b2'] * (road_a_['Roughness'] ** 2))
            b3 = (Speed.loc[road_a_['Lanes'],road_a_['Terrain'],'b3'] * (road_a_['Roughness'] ** 3))
            Speed_by_vehicle_type = (b0 + b1 + b2 + b3)
            dSpeed[y].loc[a] = (Speed_by_vehicle_type[VehicleTypes].sum() / 12)
            #Time
            time = 1 / Speed_by_vehicle_type * road_a_['Length'] * VehicleFleet['Number_of_Passengers'] * VehicleFleet['Passengers_Time_Cost'] * Load.loc[y] * 365 / 1000000
            dCostTime[y].loc[a] = time[VehicleTypes].sum()
    dCostAgencyFin = dCostCapitalFin + dCostRepairFin + dCostRecurrentFin
    dCostAgencyEco = dCostCapitalEco + dCostRepairEco + dCostRecurrentEco
    # Users and Total
    dCostTotal = dCostAgencyEco + dCostTime + dCostVOC
    # Net Benefits + NPV
    dNetTotal = dCostTotal.loc[0] - dCostTotal
    dNetTotal['list'] = dNetTotal[years].values.tolist()
    NPV = dNetTotal['list'].apply(lambda x: np.npv(DiscountRate, x))
    NPV_perkm = NPV / road['Length']
    IRR = dNetTotal['list'].apply(lambda x: np.irr(x))
    # Select best road work schedule alternative
    Selection = NPV.idxmax()
    Result = {
        'Selected_case': Selection,
        'Class': Sol['Class'].loc[Selection],
        'Code': Sol['Code'].loc[Selection],
        'Name': Sol['Name'].loc[Selection],
        'Cost': Sol['Cost'].loc[Selection],
        'Costkm': Sol['Costkm'].loc[Selection],
        'Year': Sol['Year'].loc[Selection],
        'Interval':Sol['Interval'].loc[Selection],
        'NPV': NPV.loc[Selection],
        'NPVperKm': NPV_perkm.loc[Selection],
        'IRR': IRR.loc[Selection],
        'AADTTotal': AADT['Total'].loc[1],
        'truckpct': truck_pct,
        'utilization': utilization,
        'ESA total': ESAtotal.loc[1],
        'Width' : road['Width'],   #ERROR - will take latest, not optimum
        'SurfaceType': road['SurfaceType'],    #ERROR - will take latest, not optimum
        'Roughness': dCondIRI[1].loc[Selection],
        'Structural Number' : road['StructuralNo'],   #ERROR - will take latest, not optimum
        'Pavement Age' : iCondAge[1].loc[Selection]
        }
    Outputs.append(Result)
    AADTlist.append(AADT['Total'])
    IRIlist.append(dCondIRI.loc[Selection])
    NetBenefitslist.append(dNetTotal.loc[Selection].drop('list',axis = 0))

#Output
col_order = ['Selected_case','Class','Code','Name','Cost','Costkm','Year','Interval','NPV','NPVperKm','IRR','AADTTotal','truckpct','utilization','ESA total','Width','SurfaceType','Roughness','Structural Number','Pavement Age']
Output = pd.DataFrame(Outputs, columns = col_order)
IRIdf = pd.DataFrame(IRIlist, columns = years).reset_index()
AADTdf = pd.DataFrame(AADTlist, columns = years).reset_index()
NetBenefitsdf = pd.DataFrame(NetBenefitslist, columns = years).reset_index()
Output = pd.concat([Output, IRIdf, AADTdf, NetBenefitsdf],axis = 1)
Output = Output.drop(['index','index','index'],axis = 1)
Output.to_csv(os.path.join(path, 'output.csv'))

"""
#Rodrigo

- your selection mechanism only finds local maxima, and even then badly
- some roadwork profiles are well off....what are you doing

"""
