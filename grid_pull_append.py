# Name:         grid_pull_USBRIrrAcres.py
# Purpose:      Pulls Daily GRIDMET data and calculates ETr for ETDemands model input
# Author:       Chris Pearson
# Created       2017-01-10
# Python:       2.7
#--------------------------------
import pandas as pd
import datetime as dt
import ee
import numpy as np
from time import sleep
import timeit
import sys
import os

ee.Initialize()
    
if len(sys.argv) < 2:
    print('Please Specify Username Input Argument. Example Command: python grid_ee_pull_USBRIrrAcres.py Chris')
    sys.exit()

user=str(sys.argv[1])  
if user == 'Chris':
    start = 0
    end = 1
elif user == 'Chris2':
    start = 401
    end = 700
elif user == 'Chris3':
    start = 701
    end = 1000
elif user == 'Britta':
    start = 1001
    end = 1500
elif user == 'Charles':
    start = 1501
    end = 2000 
elif user == 'Dan':
    start = 2001
    end = 2500
elif user == 'Justin':
    start = 2501
    end = 3000
elif user == 'Matt':
    start = 3000
    end = 3693
elif user == 'test':
    start = 0
    end = 0
elif user == 'all':
    start = 0
    end = 4000
else:
    print('Username Not Found. Check Spelling. Possible Options: Britta, Charles, Chris, Dan, Justin, Matt, test')
    sys.exit()
    
#Run RefET on imported dataset
RefET_flag = True
if RefET_flag: 
    import refet
    
#Input .csv containing GRIDMET_ID, LAT, LON,ELEV_M
#Create input .csv from: GIS-DATA\Shapefiles\gridmet_4km\gridmet_4km_dd_pts_full.shp 
input_df = pd.read_csv('gridmet_4km_USBRIrrAcresIntersect_UPDATE.txt')

#Specify column order for output .csv Variables: 'Date','Year','Month','Day','Tmax','Tmin','Srad','Ws_10m','q','Prcp','ETo'
output_order = ['Date','Year','Month','Day','Tmax_C','Tmin_C','Srad','Ws_2m','q','Prcp','ETo']

# Date Range to Retrieve Values (start=inclusive, end=exclusive)
# start_date = '1979-01-01'
# end_date = '2018-01-01'

start_date = '1979-01-01'
end_date = '2018-01-01'

if user == 'test':
    start_date = '1979-01-01'
    end_date = '1984-01-01'
        
#List of ee GRIDMET varibles to retrieve
met_bands = ['tmmx', 'tmmn', 'srad', 'vs', 'sph', 'pr','eto']
#The GRIDMET dataset contains the following bands:
#pr: precipitation amount (mm, daily total)
#rmax: maximum relative humidity (%)
#rmin: minimum relative humidity (%)
#sph: specific humidity (kg/kg)
#srad: surface downward shortwave radiation (W/m^2)
#th: wind direction (degrees clockwise from north)
#tmmn: minimum temperature (K)
#tmmx: maximum temperature (K)
#vs: wind velocity at 10 m (m/s)
#erc: energy release component (NFDRS fire danger index)
#eto: Daily reference evapotranspiration (grass, mm)
#bi: burning index (NFDRS fire danger index)
#fm100: 100-hour dead fuel moisture (%)
#fm1000: 1000-hour dead fuel moisture (%)
#etr: Daily reference evapotranspiration (alfalfa, mm)
#vpd: Mean vapor pressure deficit (kPa)

#Rename GRIDMET variables during ee export
met_names= ['Tmax', 'Tmin', 'Srad', 'Ws_10m', 'q', 'Prcp','ETo']

#Exponential getinfo call from ee-tools/utils.py
def ee_getinfo(ee_obj, n=30):
    """Make an exponential backoff getInfo call on the Earth Engine object"""
    output = None
    for i in range(1, n):
        try:
            output = ee_obj.getInfo()
        except Exception as e:
            print('    Resending query ({}/10)'.format(i))
            print('    {}'.format(e))
            sleep(i ** 2)
        if output:
            break
    return output
#%%
#Loop through dataframe row by row and grab desired met data from GRID collections based on Lat/Lon and Start/End dates#
for index, row in input_df.iterrows():
    start_time = timeit.default_timer()
    #Limit iteration during development#
    if index < start:
        continue
    if index > end:
        break
    lat = row.LAT
    lon = row.LON

    GRIDMET_ID_str=str(row.GRIDMET_ID)
    print(['Cell Loop Counter',index+1])
    print(['GRIDMET ID:',row.GRIDMET_ID])

    output_name = ['gridmet_historical_' + GRIDMET_ID_str + '.csv']

    if os.path.isfile(output_name[0]):
        print('{} Exists. Skipping.').format(output_name)
        continue

    point = ee.Geometry.Point(lon,lat);
    
    #Loop through ee pull by year (max 5000 records for getInfo())
    #Append each new year on end of dataframe
    #Process dataframe units and output after                        
    start_date_yr = dt.datetime.strptime(start_date,'%Y-%m-%d').year
    end_date_yr=dt.datetime.strptime(end_date,'%Y-%m-%d').year                             
    for iter_year in range(start_date_yr, end_date_yr+1, 1):
        print(iter_year)                               
    #Filter Collection by Start/End Date and Lat/Lon Point#
        gridmet_coll = ee.ImageCollection('IDAHO_EPSCOR/GRIDMET') \
                .filterDate(start_date, end_date) \
                .filter(ee.Filter.calendarRange(iter_year, iter_year, 'year')) \
                .select(met_bands,met_names)
       
        def get_values(image):
            #Pull out Date from Image
            dateStr = image.date();
            dateNum = ee.Image.constant(ee.Number.parse(dateStr.format("YYYYMMdd"))).rename(['Date'])
            #Add DateNum Band to Image
            image = image.addBands([dateNum])      
            #Reduce image taking mean of all pixels in geometry (4km resolution)
            input_mean = ee.Image(image) \
                .reduceRegion(
                        reducer = ee.Reducer.mean(), geometry=point, scale = 4000)
            return ee.Feature(None, input_mean)   
        #Run get_values fuction over all images in gridmet collection
        data=gridmet_coll.map(get_values)
    
    #Export dictionary to pandas dataframe using expornential getInfo fnc
        if iter_year == start_date_yr:
                export_df = pd.DataFrame([
                        ftr['properties']
                        for ftr in ee_getinfo(data)['features']])
        else:
                export_df2 = pd.DataFrame([
                ftr['properties']
                for ftr in ee_getinfo(data)['features']])
                export_df = export_df.append(export_df2)    
    #Reset Index
    export_df = export_df.reset_index(drop=False)

    #Convert dateNum to datetime and create Year, Month, Day, DOY variables
    export_df.Date = pd.to_datetime(export_df.Date.astype(str), format='%Y%m%d')
    export_df['Year'] = export_df['Date'].dt.year
    export_df['Month'] = export_df['Date'].dt.month
    export_df['Day'] = export_df['Date'].dt.day       
    export_df['DOY'] = export_df['Date'].dt.dayofyear 
             
    #Remove all negative Prcp values (GRIDMET Bug)
    export_df.Prcp = export_df.Prcp.clip(lower=0)       

    #Convert 10m windspeed to 2m (ASCE Eqn. 33)
    zw = 10
    export_df['Ws_2m'] = export_df.Ws_10m*((4.87/np.log(67.8*zw-5.42)))

    #Unit Conversions
    export_df.Tmax = export_df.Tmax-273.15; #K to C
    export_df.Tmin = export_df.Tmin-273.15; #K to C
    export_df.rename(columns = {'Tmax': 'Tmax_C','Tmin':'Tmin_C'}, inplace=True)
      
    if RefET_flag: 
        #Create RefET Specific Dataframe
        ref_df = pd.DataFrame()
        ref_df['Date'] = export_df.Date
        ref_df['DOY'] = export_df.DOY      
        #Temperature Min/Max (C)
        ref_df['Tmin'] = export_df['Tmin_C']
        ref_df['Tmax'] = export_df['Tmax_C']
        
        #Pressure from elevation (ASCE Eqn. 34) kPa
        P= 101.3*((293-0.0065*row.ELEV_M)/293)**5.26
            
        #Specific Humidity and elevation to Vapor Pressure (kPa)
        ref_df['ea'] = export_df.q*P/(0.622+0.378*export_df.q)      
              
        #Srad to MJ m-2 day-1
        ref_df['Srad_MJ'] = 0.0864*export_df.Srad
              
        #Windspeed Sensor Ht (m)      
        ref_df['Ws_10m'] = export_df.Ws_10m
#        ref_df['Ws_2m'] = export_df.Ws_2m
              
        #Windspeed Sensor Ht (m)      
        ref_df['Sensor_Ht'] = zw     
        
        #Elevation (m)                
        ref_df['ELEV_M'] = row.ELEV_M
                
        #Latitude (Radians)
        ref_df['LAT_rad'] = lat*np.pi/180
               
        #Run RefET code on ref_df dataframe
        ref_df['ETr_ASCE'] = refet.daily(ref_df.Tmin, ref_df.Tmax, ref_df.ea, ref_df.Srad_MJ, ref_df.Ws_10m, ref_df.Sensor_Ht, ref_df.ELEV_M, ref_df.LAT_rad, ref_df.DOY, ref_type='etr',
              rso_type='full', rso=None, asce_flag=True)
#        ref_df['ETo_ASCE']= RefET.daily(ref_df.Tmin, ref_df.Tmax, ref_df.ea, ref_df.Srad_MJ, ref_df.Ws_10m, ref_df.Sensor_Ht, ref_df.ELEV_M, ref_df.LAT_rad, ref_df.DOY, ref_type='eto',
#              rso_type='full', rso=None, asce_flag=True)

        #Remove all negative evaporation estimates
        ref_df.ETr_ASCE = ref_df.ETr_ASCE.clip(lower=0)
        #Copy ETr_ASCE data to export_df
        export_df['ETr_ASCE'] = ref_df.ETr_ASCE
        #RefET Flag Output Order
        output_order = ['Date','Year','Month','Day','Tmax_C','Tmin_C','Srad','Ws_2m','q','Prcp','ETo','ETr_ASCE']
    
    #Output Filename with GRIDMET ID
    # output_name = ['gridmet_historical_'+GRIDMET_ID_str+'ee_test'+'.csv']

    # # Update files with new data???
    # original_df = pd.read_csv(output_name)
    # update_df = original_df.update(export_df)

    #Write csv files to working directory
    export_df.to_csv(output_name[0],columns=output_order)
    elapsed = timeit.default_timer() - start_time
    print(elapsed)   
   
