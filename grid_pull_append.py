import datetime as dt
import ee
import numpy as np
import os
import pandas as pd
import refet
import sys
from time import sleep
import timeit

ee.Initialize()

os.chdir('D:\etr-biascorrect')

# Input .csv containing GRIDMET_ID, LAT, LON
input_df = pd.read_csv('D:\etr-biascorrect\gridmet_testlist.txt')

# List of ee GRIDMET varibles to retrieve
# https://explorer.earthengine.google.com/#detail/IDAHO_EPSCOR%2FGRIDMET
met_bands = ['tmmx', 'tmmn', 'srad', 'vs', 'sph', 'rmin', 'rmax',
             'pr', 'etr', 'eto']

# Rename GRIDMET variables during ee export
met_names= ['tmax', 'tmin', 'srad_wm2', 'u10_ms', 'q_kgkg', 'rh_min', 'rh_max',
            'prcp_mm', 'etr_mm', 'eto_mm']

# Specify column order for output .csv Variables:
output_order = ['date', 'year', 'month', 'day', 'centroid_lat', 'centroid_lon',
                'elev_m', 'u2_ms', 'tmin_c', 'tmax_c', 'srad_wm2', 'ea_kpa',
                'prcp_mm', 'etr_mm', 'eto_mm']

# Exponential getinfo call from ee-tools/utils.py
def ee_getinfo(ee_obj, n=30):
    """Make an exponential backoff getInfo call on the Earth Engine object"""
    output = None
    for i in range(1, n):
        try:
            output = ee_obj.getInfo()
        except Exception as e:
            print('    Resending query ({}/{})'.format(i, n))
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
    print('Processing GRIDMET ID: {}'.format(GRIDMET_ID_str))

    # determine end date of data collection
    current_date = dt.datetime.today()
    end_date = dt.date(current_date.year, current_date.month,
                       current_date.day-1)
    # Create List of all dates
    full_date_list = pd.date_range(dt.datetime.strptime('2018-01-01',
                                                        '%Y-%m-%d'), end_date)
    output_name = 'gridmet_historical_' + GRIDMET_ID_str + '.csv'
    output_file = os.path.join('C:\etr-biascorrect', output_name)
    
    if os.path.isfile(output_file):
        print('{} Exists. Checking for missing data.'.format(output_name))
        original_df = pd.read_csv(output_file, parse_dates=True)
        missing_dates = list(set(full_date_list) - set(pd.to_datetime(
            original_df['date'])))
        original_df.date = pd.to_datetime(original_df.date.astype(str),
                                          format='%Y-%m-%d')
        original_df = original_df.Date.apply(lambda x: x.strftime('%Y-%m-%d'))
    else:
        print('{} does not exists. Creating file.'.format(output_name))
        missing_dates = list(set(full_date_list))
    if not missing_dates:
        print('No missing data found. Skipping')
        continue
    
    # Min and Max of Missing Dates (Start: Inclusive; End: Exclusive)
    start_date = min(missing_dates)
    end_date = max(missing_dates)

    missing_years = []
    for date in missing_dates:
        missing_years = missing_years + [date.year]
    missing_years = sorted(list(set(missing_years)))

    # Add check to verify lat/lon fall within the gridmet extent
    # -124.78749996666667 25.04583333333334 -67.03749996666667 49.42083333333334
    # Create ee point from lat and lon
    point = ee.Geometry.Point(row.LON, row.LAT)

    # gridmet elevation image:
    # ee.Image('projects/climate-engine/gridmet/elevation')
    elev= ee.Image('projects/climate-engine/gridmet/elevation') \
        .reduceRegion(reducer=ee.Reducer.mean(), geometry=point,
                      scale=4000)
    elev = ee_getinfo(elev)['b1']

    # Calculate out grid cell centroid
    # gridMET elevation asset lower left corner coordinates
    gridmet_lon = -124.78749996666667
    gridmet_lat = 25.04583333333334
    gridmet_cs = 0.041666666666666664
    gridcell_lat = int(abs(row.LAT - gridmet_lat) / gridmet_cs) * gridmet_cs +\
                   gridmet_lat + gridmet_cs/2
    gridcell_lon = int(abs(row.LON - gridmet_lon) / gridmet_cs) * gridmet_cs +\
                   gridmet_lon + gridmet_cs/2

    # Loop through ee pull by year (max 5000 records for getInfo())
    # Append each new year on end of dataframe
    # Process dataframe units and output after                        
    start_date_yr = start_date.year
    end_date_yr = end_date.year                             
    # for iter_year in range(start_date_yr, end_date_yr+1, 1):
    for iter_year in missing_years:
        print(iter_year)                               
    # Filter Collection by start/end date and lat/lon Point
    # Only include 'permanent' data
        gridmet_coll = ee.ImageCollection('IDAHO_EPSCOR/GRIDMET') \
            .filterDate(start_date, end_date+1) \
            .filter(ee.Filter.calendarRange(iter_year, iter_year, 'year')) \
            .filter(ee.Filter.eq('status', 'permanent')) \
            .select(met_bands, met_names)

        # Check if collection is empty        
        image_count = ee.Number(gridmet_coll.limit(1) \
                                .reduceColumns('count', ['system:index']) \
                                .get('count'))
        empty = ee.Algorithms.If(image_count.eq(1), False, True)

        if ee_getinfo(empty):
            print('No new "permanent" data found. Skipping.')
            export_df = None
            continue

        def get_values(image):
            # Pull out date from Image
            datestr = image.date()
            datenum = ee.Image.constant(ee.Number.parse(
                datestr.format("YYYYMMdd"))).rename(['date'])
            # Add dateNum Band to Image
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
    export_df.date = pd.to_datetime(export_df.date.astype(str), format='%Y%m%d')
    export_df['year'] = export_df['date'].dt.year
    export_df['month'] = export_df['date'].dt.month
    export_df['day'] = export_df['date'].dt.day
    # export_df['DOY'] = export_df['Date'].dt.dayofyear
    # Format Date for export
    export_df['date'] = export_df.date.apply(lambda x: x.strftime('%Y-%m-%d'))

    # Remove all negative Prcp values (GRIDMET Bug)
    export_df.prcp_mm = export_df.prcp_mm.clip(lower=0)

    # Convert 10m windspeed to 2m (ASCE Eqn. 33)
    zw = 10
    export_df['u2_ms'] = refet.calcs._wind_height_adjust(export_df.u10_ms, zw)
    # elevation from gridMET elevation layer
    export_df['elev_m'] = elev
    export_df['centroid_lat'] = gridcell_lat
    export_df['centroid_lon'] = gridcell_lon

    # air pressure from gridmet elevation using refet module
    export_df['pair_kpa'] = refet.calcs._air_pressure(export_df.elev_m,
                                                      method='asce')

    # actual vapor pressure (kg/kg) using refet module
    export_df['ea_kpa'] = refet.calcs._actual_vapor_pressure(export_df.q_kgkg,
                                                             export_df.pair_kpa)

    # Unit Conversions
    export_df.tmax = export_df.tmax-273.15  #K to C
    export_df.tmin = export_df.tmin-273.15  #K to C
    export_df.rename(columns={'tmax': 'tmax_c', 'tmin': 'tmin_c'}, inplace=True)
    # export_df['Tavg_C'] = (export_df.Tmax_C + export_df.Tmin_C)/2

    # Relative Humidity from gridMET min and max
    # export_df['RH_avg'] = (export_df.RH_max + export_df.RH_min)/2

    # Add new data to original dataframe, remove duplicates
    export_df = pd.concat([original_df, export_df], ignore_index=True, sort=True)
    export_df = export_df[output_order].drop_duplicates('date')
    export_df = export_df.sort_values(by=['year', 'month', 'day'])
    
    # Write csv files to working directory
    export_df.to_csv(output_name, columns=output_order, index=False)
    elapsed = timeit.default_timer() - start_time
    print(elapsed)   


# Saturated vapor pressure
# export_df['esat_min_kPa'] = refet.calcs._sat_vapor_pressure(export_df.Tmin_C)
# export_df['esat_max_kPa'] = refet.calcs._sat_vapor_pressure(export_df.Tmax_C)
# export_df['esat_avg_kPa'] = (export_df.esat_min_kPa + export_df.esat_max_kPa)/2