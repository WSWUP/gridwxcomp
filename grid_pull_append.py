import pandas as pd
import datetime as dt
import ee
import numpy as np
from time import sleep
import timeit
#import sys
import os

ee.Initialize()

<<<<<<< HEAD
os.chdir('D:\etr-biascorrect')

#Input .csv containing GRIDMET_ID, LAT, LON,ELEV_M
#Create input .csv from: GIS-DATA\Shapefiles\gridmet_4km\gridmet_4km_dd_pts_full.shp 
input_df = pd.read_csv('D:\etr-biascorrect\gridmet_testlist.txt')



=======
os.chdir('C:\etr-biascorrect')

#Input .csv containing GRIDMET_ID, LAT, LON,ELEV_M
#Create input .csv from: GIS-DATA\Shapefiles\gridmet_4km\gridmet_4km_dd_pts_full.shp 
input_df = pd.read_csv('C:\etr-biascorrect\gridmet_testlist.txt')

#Specify column order for output .csv Variables: 'Date','Year','Month','Day','Tmax','Tmin','Srad','Ws_10m','q','Prcp','ETo'
output_order = ['Date', 'Year', 'Month', 'Day', 'Tmax_C', 'Tmin_C', 'Srad',
                'Ws_2m', 'q', 'Prcp', 'ETo']

>>>>>>> master
#List of ee GRIDMET varibles to retrieve
met_bands = ['tmmx', 'tmmn', 'srad', 'vs', 'sph', 'pr','etr']
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
met_names= ['Tmax', 'Tmin', 'Srad', 'Ws_10m', 'q', 'Prcp','ETr']
<<<<<<< HEAD

#Specify column order for output .csv Variables: 'Date','Year','Month','Day','Tmax','Tmin','Srad','Ws_10m','q','Prcp','ETo'
output_order = ['Date', 'Year', 'Month', 'Day', 'Tmax_C', 'Tmin_C', 'Srad',
                'Ws_2m', 'q', 'Prcp', 'ETr']
=======
>>>>>>> master

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
<<<<<<< HEAD

=======
#%%
>>>>>>> master
# Loop through dataframe row by row and grab desired met data from
# GRID collections based on Lat/Lon and Start/End dates#
for index, row in input_df.iterrows():
    start_time = timeit.default_timer()
<<<<<<< HEAD
    # Reset original_df
    original_df = None
    
    GRIDMET_ID_str=str(row.GRIDMET_ID)
    # print(['Cell Loop Counter',index+1])
    print('Processing GRIDMET ID: {}',format(GRIDMET_ID_str))
=======
    # #Limit iteration during development#
    # if index < start:
    #     continue
    # if index > end:
    #     break
    lat = row.LAT
    lon = row.LON
    GRIDMET_ID_str=str(row.GRIDMET_ID)
    print(['Cell Loop Counter',index+1])
    print(['GRIDMET ID:',row.GRIDMET_ID])
>>>>>>> master

    # determine end date of data collection
    current_date = dt.datetime.today()
 
<<<<<<< HEAD
    end_date = dt.date(current_date.year, current_date.month,
                                 current_date.day-1)
    # Create List of all dates
    full_date_list = pd.date_range(dt.datetime.strptime('1979-01-01','%Y-%m-%d'), end_date)    
    output_name = 'gridmet_historical_' + GRIDMET_ID_str + '.csv'
    output_file = os.path.join('D:\etr-biascorrect', output_name)
    
    if os.path.isfile(output_file):
        print('{} Exists. Checking for missing data.'.format(output_name))
        original_df = pd.read_csv(output_file, parse_dates=True)
        missing_dates = list(set(full_date_list) - set(pd.to_datetime(original_df['Date'])))
    else:
        print('{} does not exists. Creating file.'.format(output_name))
        missing_dates = list(set(full_date_list))
    if not missing_dates:
        print('No missing data found. Skipping')
        continue
    

    # Min and Max of Missing Dates (Start: Inclusive; End: Exclusive)
    start_date = min(missing_dates)
    end_date = max(missing_dates)

    # Create ee point from lat and lon
    point = ee.Geometry.Point(row.LON, row.LAT);
=======
    provisional_flag = False

    # gridMET data is provision for approx. 2 months
    # remove provisional data unless provisional flag is set to TRUE
    if provisional_flag:
        end_date = current_date
    else:
        end_date = dt.date(current_date.year, current_date.month-2,
                                 current_date.day)

    output_name = 'gridmet_historical_' + GRIDMET_ID_str + '.csv'
    output_file = os.path.join('C:\etr-biascorrect', output_name)
    if os.path.isfile(output_file):
        print('{} Exists. Updating file.').format(output_name)
        original_df = pd.read_csv(output_file)
    else:
        print('{} Does Not Exists. Creating file.').format(output_name)
        # Inclusive
        start_date = dt.datetime.strptime('1979-01-01','%Y-%m-%d')
    # Create List of all dates
    def daterange(date1, date2):
        for n in range(int ((date2 - date1).days)+1):
            yield dt.datetime.strptime(date1 + dt.timedelta(n))
    full_date_list = daterange(start_date, end_date)

    # Find missing dates in DF
    missing_dates = full_date_list - original_df['Date']

    # Min and Max of Missing Dates
    start_date = missing_dates.min()

    # Create ee point from lat and lon
    point = ee.Geometry.Point(lon, lat);
>>>>>>> master
    
    #Loop through ee pull by year (max 5000 records for getInfo())
    #Append each new year on end of dataframe
    #Process dataframe units and output after                        
    start_date_yr = start_date.year
    end_date_yr = end_date.year                             
    for iter_year in range(start_date_yr, end_date_yr+1, 1):
        print(iter_year)                               
    #Filter Collection by Start/End Date and Lat/Lon Point
    #Should image status be an input argument ('early', 'provisional', or 'permanent')
        gridmet_coll = ee.ImageCollection('IDAHO_EPSCOR/GRIDMET') \
                .filterDate(start_date, end_date+1) \
                .filter(ee.Filter.calendarRange(iter_year, iter_year, 'year')) \
                .filter(ee.Filter.eq('status', 'permanent')) \
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
        data = gridmet_coll.map(get_values)
    
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
    # Reset Index
    export_df = export_df.reset_index(drop=False)

    # Convert dateNum to datetime and create Year, Month, Day, DOY variables
    export_df.Date = pd.to_datetime(export_df.Date.astype(str), format='%Y%m%d')
    export_df['Year'] = export_df['Date'].dt.year
    export_df['Month'] = export_df['Date'].dt.month
    export_df['Day'] = export_df['Date'].dt.day       
    # export_df['DOY'] = export_df['Date'].dt.dayofyear
    # Format Date for export
    export_df['Date'] = export_df.Date.apply(lambda x: x.strftime('%Y-%m-%d'))
             
    # Remove all negative Prcp values (GRIDMET Bug)
    export_df.Prcp = export_df.Prcp.clip(lower=0)       

    # Convert 10m windspeed to 2m (ASCE Eqn. 33)
    zw = 10
    export_df['Ws_2m'] = export_df.Ws_10m*((4.87/np.log(67.8*zw-5.42)))

    # Unit Conversions
    export_df.Tmax = export_df.Tmax-273.15; #K to C
    export_df.Tmin = export_df.Tmin-273.15; #K to C
    export_df.rename(columns = {'Tmax': 'Tmax_C','Tmin':'Tmin_C'}, inplace=True)
      
<<<<<<< HEAD

    # Add new data to original dataframe, remove duplicates
    export_df = pd.concat([original_df, export_df], sort = True)
    export_df = export_df[output_order].drop_duplicates().reset_index(drop=False)
    # Write csv files to working directory
    export_df.to_csv(output_name, columns=output_order)
=======

    # Update files with new data???
    export_df = original_df.update(export_df).drop_duplicates().reset_index(drop=1).to_frame()

    # Write csv files to working directory
    export_df.to_csv(output_name[0], columns=output_order)
>>>>>>> master
    elapsed = timeit.default_timer() - start_time
    print(elapsed)   
   
