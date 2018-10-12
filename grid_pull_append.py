import pandas as pd
import datetime as dt
import ee
import numpy as np
from time import sleep
import timeit
import os

ee.Initialize()

os.chdir('D:\etr-biascorrect')

# Input .csv containing GRIDMET_ID, LAT, LON,ELEV_M
input_df = pd.read_csv('D:\etr-biascorrect\gridmet_testlist.txt')

# List of ee GRIDMET varibles to retrieve
# https://explorer.earthengine.google.com/#detail/IDAHO_EPSCOR%2FGRIDMET
met_bands = ['tmmx', 'tmmn', 'srad', 'vs', 'sph', 'pr', 'etr']

# Rename GRIDMET variables during ee export
met_names= ['Tmax', 'Tmin', 'Srad', 'Ws_10m', 'q', 'Prcp', 'ETr']

# Specify column order for output .csv Variables: 'Date','Year','Month','Day',
# 'Tmax','Tmin','Srad','Ws_10m','q','Prcp','ETo'
output_order = ['Date', 'Year', 'Month', 'Day', 'Tmax_C', 'Tmin_C', 'Srad',
                'Ws_2m', 'q', 'Prcp', 'ETr']

# Exponential getinfo call from ee-tools/utils.py
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


# Loop through dataframe row by row and grab desired met data from
# GRID collections based on Lat/Lon and Start/End dates
for index, row in input_df.iterrows():
    start_time = timeit.default_timer()
    # Reset original_df
    original_df = None
    export_df = None
    
    GRIDMET_ID_str=str(row.GRIDMET_ID)
    # print(['Cell Loop Counter',index+1])
    print('Processing GRIDMET ID: {}', format(GRIDMET_ID_str))

    # determine end date of data collection
    current_date = dt.datetime.today()
    end_date = dt.date(current_date.year, current_date.month,
                       current_date.day-1)
    # Create List of all dates
    full_date_list = pd.date_range(dt.datetime.strptime('1979-01-01',
                                                        '%Y-%m-%d'), end_date)
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
    point = ee.Geometry.Point(row.LON, row.LAT)
    
    # Loop through ee pull by year (max 5000 records for getInfo())
    # Append each new year on end of dataframe
    # Process dataframe units and output after                        
    start_date_yr = start_date.year
    end_date_yr = end_date.year                             
    for iter_year in range(start_date_yr, end_date_yr+1, 1):
        print(iter_year)                               
    # Filter Collection by Start/End Date and Lat/Lon Point
    # Should image status be an input argument: 'early', 'provisional', 'permanent'
        gridmet_coll = ee.ImageCollection('IDAHO_EPSCOR/GRIDMET') \
            .filterDate(start_date, end_date+1) \
            .filter(ee.Filter.calendarRange(iter_year, iter_year, 'year')) \
            .filter(ee.Filter.eq('status', 'permanent')) \
            .select(met_bands, met_names)
        # Check if collection is empty        
        image_count = ee.Number(gridmet_coll.limit(1).
                                reduceColumns('count', ['system:index']).get('count'))
        empty = ee.Algorithms.If(image_count.neq(1), 'True', 'False')
        if empty:       
            print('No new "permanent" data found. Skipping.')
            export_df = None
            continue

        def get_values(image):
            # Pull out Date from Image
            datestr = image.date()
            datenum = ee.Image.constant(ee.Number.parse(datestr.format("YYYYMMdd"))).rename(['Date'])
            # Add DateNum Band to Image
            image = image.addBands([datenum])
            # Reduce image taking mean of all pixels in geometry (4km resolution)
            input_mean = ee.Image(image) \
                .reduceRegion(
                        reducer=ee.Reducer.mean(), geometry=point, scale=4000)
            return ee.Feature(None, input_mean)   
        # Run get_values function over all images in gridmet collection
        data = gridmet_coll.map(get_values)
    
    # Export dictionary to pandas dataframe using exponential getInfo fnc
        if iter_year == start_date_yr:
                export_df = pd.DataFrame([
                        ftr['properties']
                        for ftr in ee_getinfo(data)['features']])
        else:
                export_df2 = pd.DataFrame([
                ftr['properties']
                for ftr in ee_getinfo(data)['features']])
                export_df = export_df.append(export_df2)    
    # If export_df is None (skip to next ID)
    if export_df is None:
        continue
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
    export_df.Tmax = export_df.Tmax-273.15  #K to C
    export_df.Tmin = export_df.Tmin-273.15  #K to C
    export_df.rename(columns = {'Tmax': 'Tmax_C','Tmin':'Tmin_C'}, inplace=True)

    # Add new data to original dataframe, remove duplicates
    export_df = pd.concat([original_df, export_df], sort=True)
    export_df = export_df[output_order].drop_duplicates().reset_index(drop=False)
    
    # Write csv files to working directory
    export_df.to_csv(output_name, columns=output_order)
    elapsed = timeit.default_timer() - start_time
    print(elapsed)   
