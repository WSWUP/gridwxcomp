# -*- coding: utf-8 -*-
"""
Download gridMET climatic time series for multiple variables, e.g. ETr, temp,
wind speed 2m, swrad, etc. 

"""
import argparse
import logging
import os
import sys
import timeit
import datetime as dt
from time import sleep

import ee
import refet
import pandas as pd

# for connection to Earth Engine
ee.Initialize() 

def download_gridmet_ee(input_csv, out_folder, year_filter='', year_update=''): 
    """
    Download gridMET time series data for multiple climate variables for 
    select gridMET cells as listed in ``input_csv``.

    Arguments:
        input_csv (str): file path of input CSV produced by :mod:`gridwxcomp.prep_input`
        out_folder (str): directory path to save gridmet timeseries CSV files

    Keyword Arguments:
        year_filter (str): default ''. Single year YYYY or range YYYY-YYYY 
            to download.
        year_update (str): default ''. Re-download existing data for year or
            range, YYYY or YYYY-YYYY.

    Returns:
        None

    Examples:
        Say we wanted to download data for 2016 through 2018, from the command 
        line,

        .. code-block:: sh

            $ download_gridmet_ee.py -i merged_input.csv -o gridmet_data -y 2016-2018

        note, "merged_input.csv" should have been created by first running 
        :mod:`gridwxcomp.prep_input`. If the ``[-y, --years]`` option is not 
        given the default behaviour is to download gridMET data from 1979 up 
        through yesterday.

        If the data for 2018 has changed since the last run or for debugging
        purposes you can re-download data for all or select years with the
        ``[-u, --update-data]`` option

        .. code-block:: sh

            $ python download_gridmet_ee.py -i merged_input.csv -o gridmet_data -y 2018 -u

        To download the same gridMET data within Python

        >>> from gridwxcomp import download_gridmet_ee
        >>> download_gridmet_ee('merged_input.csv', 'gridmet_data', '2016-2018')

        Running :func:`download_gridmet_ee` also updates the CSV file
        produced from :mod:`gridwxcomp.prep_input` to include file paths to 
        gridMET time series files that are paired with climate stations. 
    """
    if not os.path.exists(out_folder):
        logging.info('\nCreating output folder: {}'.format(out_folder))
        os.makedirs(out_folder)

    # Input .csv containing GRIDMET_ID, LAT, LON
    input_df = pd.read_csv(input_csv)

    # List of ee GRIDMET varibles to retrieve
    # https://explorer.earthengine.google.com/#detail/IDAHO_EPSCOR%2FGRIDMET
    met_bands = ['tmmx', 'tmmn', 'srad', 'vs', 'sph', 'rmin', 'rmax',
                 'pr', 'etr', 'eto']

    # Rename GRIDMET variables during ee export
    met_names = ['tmax', 'tmin', 'srad_wm2', 'u10_ms', 'q_kgkg', 'rh_min',
                 'rh_max', 'prcp_mm', 'etr_mm', 'eto_mm']

    # Specify column order for output .csv Variables:
    output_order = ['date', 'year', 'month', 'day', 'centroid_lat',
                    'centroid_lon', 'elev_m', 'u2_ms', 'tmin_c', 'tmax_c',
                    'srad_wm2', 'ea_kpa', 'prcp_mm', 'etr_mm', 'eto_mm']

    # Year Filter
    if year_filter:
            year_list = sorted(list(_parse_int_set(year_filter)))
            logging.info('\nDownloading Years: {0}-{1}'.format(min(year_list),
                                                               max(year_list)))
            date_list = pd.date_range(
                dt.datetime.strptime('{}-01-01'.format(min(year_list)),
                                     '%Y-%m-%d'),
                dt.datetime.strptime('{}-12-31'.format(max(year_list)),
                                     '%Y-%m-%d'))
    else:
        logging.info('\nDownloading full historical record (1979-present).')
        # Create List of all dates
        # determine end date of data collection
        current_date = dt.datetime.today()
        end_date = dt.date(current_date.year, current_date.month,
                           current_date.day - 1)
        date_list = pd.date_range(dt.datetime.strptime('1979-01-01',
                                                       '%Y-%m-%d'), end_date)

    # Year Update List
    if year_update:
        update_list = sorted(list(_parse_int_set(year_update)))
        logging.info('\nUpdating Years: {0}-{1}'.format(min(update_list),
                                                           max(update_list)))
        # update_list = pd.date_range(
        #     dt.datetime.strptime('{}-01-01'.format(min(year_list)),
        #                          '%Y-%m-%d'),
        #     dt.datetime.strptime('{}-12-31'.format(max(year_list)),
        #                          '%Y-%m-%d'))

    # Exponential getinfo call from ee-tools/utils.py
    def ee_getinfo(ee_obj, n=30):
        """Make an exponential backoff getInfo call on the EarthEngine object"""
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

        GRIDMET_ID_str = str(row.GRIDMET_ID)
        logging.info('\nProcessing GRIDMET ID: {}'.format(GRIDMET_ID_str))

        output_name = 'gridmet_historical_' + GRIDMET_ID_str + '.csv'
        output_file = os.path.join(out_folder, output_name)

        if os.path.isfile(output_file):
            logging.info('{} exists. Checking for missing data.'.format(
                output_name))
            original_df = pd.read_csv(output_file, parse_dates=True)
            # Apply update filter (remove original data based on year)
            if year_update:
                original_df = original_df[~original_df['year']
                    .isin(update_list)]

            missing_dates = list(set(date_list) - set(pd.to_datetime(
                original_df['date'])))
            original_df.date = pd.to_datetime(original_df.date.astype(str),
                                              format='%Y-%m-%d')
            original_df['date'] = original_df.date.apply(lambda x: x.strftime(
                '%Y-%m-%d'))
        else:
            logging.info('{} does not exists. Creating file.'.format(
                output_name))
            missing_dates = list(set(date_list))
        if not missing_dates:
            logging.info('No missing data found. Skipping')
            # Add gridMET file path to input table if not already there
            input_df.loc[input_df.GRIDMET_ID == row.GRIDMET_ID,\
                'GRIDMET_FILE_PATH'] = os.path.abspath(output_file)
            input_df.to_csv(input_csv, index=False)
            continue

        # Min and Max of Missing Dates (Start: Inclusive; End: Exclusive)
        start_date = min(missing_dates)
        end_date = max(missing_dates)

        missing_years = []
        for date in missing_dates:
            missing_years = missing_years + [date.year]
        missing_years = sorted(list(set(missing_years)))

        # Add check to verify lat/lon fall within the gridmet extent
        # -124.78749996666667 25.04583333333334
        # -67.03749996666667 49.42083333333334
        # Create ee point from lat and lon
        point = ee.Geometry.Point(row.LON, row.LAT)

        # gridmet elevation image:
        # ee.Image('projects/climate-engine/gridmet/elevation')
        elev = ee.Image('projects/climate-engine/gridmet/elevation') \
            .reduceRegion(reducer=ee.Reducer.mean(), geometry=point,
                          scale=4000)
        elev = ee_getinfo(elev)['b1']

        # Calculate out grid cell centroid
        # gridMET elevation asset lower left corner coordinates
        gridmet_lon = -124.78749996666667
        gridmet_lat = 25.04583333333334
        gridmet_cs = 0.041666666666666664
        gridcell_lat = int(
            abs(row.LAT - gridmet_lat) / gridmet_cs) * gridmet_cs +\
            gridmet_lat + gridmet_cs/2
        gridcell_lon = int(
            abs(row.LON - gridmet_lon) / gridmet_cs) * gridmet_cs +\
            gridmet_lon + gridmet_cs/2

        # Loop through ee pull by year (max 5000 records for getInfo())
        # Append each new year on end of dataframe
        # Process dataframe units and output after
        start_date_yr = start_date.year

        for iter_year in missing_years:
            logging.info(iter_year)
        # Filter Collection by start/end date and lat/lon Point
        # Only include 'permanent' data
            gridmet_coll = ee.ImageCollection('IDAHO_EPSCOR/GRIDMET') \
                .filterDate(start_date, end_date+1) \
                .filter(ee.Filter.calendarRange(iter_year, iter_year, 'year')) \
                .filter(ee.Filter.eq('status', 'permanent')) \
                .select(met_bands, met_names)

            # Check if collection is empty
            image_count = ee.Number(gridmet_coll.limit(1)
                                    .reduceColumns('count', ['system:index'])
                                    .get('count'))
            empty = ee.Algorithms.If(image_count.eq(1), False, True)

            if ee_getinfo(empty):
                logging.info('No new "permanent" data found. Skipping.')
                export_df = None
                continue

            def get_values(image):
                # Pull out date from Image
                datestr = image.date()
                datenum = ee.Image.constant(ee.Number.parse(
                    datestr.format("YYYYMMdd"))).rename(['date'])
                # Add dateNum Band to Image
                image = image.addBands([datenum])
                # Reduce image taking mean of all pixels in geometry (4km res)
                input_mean = ee.Image(image) \
                    .reduceRegion(
                            reducer=ee.Reducer.mean(), geometry=point,
                            scale=4000)
                return ee.Feature(None, input_mean)
            # Run get_values function over all images in gridmet collection
            data = gridmet_coll.map(get_values)

            # Export dictionary to pandas dataframe using exponential getInfo
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
        export_df.date = pd.to_datetime(export_df.date.astype(str),
                                        format='%Y%m%d')
        export_df['year'] = export_df['date'].dt.year
        export_df['month'] = export_df['date'].dt.month
        export_df['day'] = export_df['date'].dt.day
        # export_df['DOY'] = export_df['Date'].dt.dayofyear
        # Format Date for export
        export_df['date'] = export_df.date.apply(lambda x: x.strftime(
            '%Y-%m-%d'))

        # Remove all negative Prcp values (GRIDMET Bug)
        export_df.prcp_mm = export_df.prcp_mm.clip(lower=0)

        # Convert 10m windspeed to 2m (ASCE Eqn. 33)
        zw = 10
        export_df['u2_ms'] = refet.calcs._wind_height_adjust(
            export_df.u10_ms, zw)
        # elevation from gridMET elevation layer
        export_df['elev_m'] = elev
        export_df['centroid_lat'] = gridcell_lat
        export_df['centroid_lon'] = gridcell_lon

        # air pressure from gridmet elevation using refet module
        export_df['pair_kpa'] = refet.calcs._air_pressure(export_df.elev_m,
                                                          method='asce')

        # actual vapor pressure (kg/kg) using refet module
        export_df['ea_kpa'] = refet.calcs._actual_vapor_pressure(
            export_df.q_kgkg, export_df.pair_kpa)

        # Unit Conversions
        export_df.tmax = export_df.tmax-273.15  # K to C
        export_df.tmin = export_df.tmin-273.15  # K to C
        export_df.rename(columns={'tmax': 'tmax_c', 'tmin': 'tmin_c'},
                         inplace=True)
        # export_df['Tavg_C'] = (export_df.Tmax_C + export_df.Tmin_C)/2

        # Relative Humidity from gridMET min and max
        # export_df['RH_avg'] = (export_df.RH_max + export_df.RH_min)/2

        # Add new data to original dataframe, remove duplicates
        export_df = pd.concat([original_df, export_df], ignore_index=True,
                              sort=True)
        export_df = export_df[output_order].drop_duplicates('date')
        export_df = export_df.sort_values(by=['year', 'month', 'day'])
        export_df = export_df.dropna()

        # Add gridMET file path to input table
        input_df.loc[input_df.GRIDMET_ID == row.GRIDMET_ID,\
                'GRIDMET_FILE_PATH'] = os.path.abspath(output_file)
        input_df.to_csv(input_csv, index=False)

        # Write csv files to working directory
        export_df.to_csv(output_file, columns=output_order, index=False)
        elapsed = timeit.default_timer() - start_time
        logging.info('\nDownload Time: {}'.format(elapsed))


def _parse_int_set(nputstr=""):
    """Return list of numbers given a string of ranges

    http://thoughtsbyclayg.blogspot.com/2008/10/
    parsing-list-of-numbers-in-python.html
    """
    selection = set()
    invalid = set()
    # tokens are comma separated values
    tokens = [x.strip() for x in nputstr.split(',')]
    for i in tokens:
        try:
            # typically tokens are plain old integers
            selection.add(int(i))
        except:
            # if not, then it might be a range
            try:
                token = [int(k.strip()) for k in i.split('-')]
                if len(token) > 1:
                    token.sort()
                    # we have items seperated by a dash
                    # try to build a valid range
                    first = token[0]
                    last = token[len(token) - 1]
                    for x in range(first, last + 1):
                        selection.add(x)
            except:
                # not an int and not a range...
                invalid.add(i)
    # Report invalid tokens before returning valid selection
    # print "Invalid set: " + str(invalid)
    return selection


def arg_parse():
    """
    Command line usage of download_gridmet_ee.py for downloading
    gridMET time series of several climatic variables using google 
    earth engine API.
    """
    parser = argparse.ArgumentParser(
        description=arg_parse.__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    optional = parser._action_groups.pop() # optionals listed second
    required = parser.add_argument_group('required arguments')
    required.add_argument(
        '-i', '--input', metavar='', required=True,
        help='Input file containing station and gridMET IDs created by '+\
            'prep_input.py')
    required.add_argument(
        '-o', '--out-dir', metavar='', required=True,
        help='Output directory to save time series CSVs of gridMET data')
    optional.add_argument(
        '-y', '--years', metavar='', default=None, type=str,
        help='Year(s) to download, single year (YYYY) or range (YYYY-YYYY)')
    optional.add_argument(
        '-u', '--update', metavar='', default=None, type=str,
        help='Year(s) to update, single year (YYYY) or range (YYYY-YYYY)')
    optional.add_argument(
        '--debug', default=logging.INFO, const=logging.DEBUG,
        help='Debug level logging', action="store_const", dest="loglevel")
    parser._action_groups.append(optional)# to avoid optionals listed first
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = arg_parse()

    logging.basicConfig(level=args.loglevel, format='%(message)s')
    logging.info('\n{}'.format('#' * 80))
    logging.info('{0:<20s} {1}'.format(
        'Run Time Stamp:', dt.datetime.now().isoformat(' ')))
    logging.info('{0:<20s} {1}'.format('Current Directory:', os.getcwd()))
    logging.info('{0:<20s} {1}'.format(
        'Script:', os.path.basename(sys.argv[0])))

    download_gridmet_ee(input_csv=args.input, out_folder=args.out_dir,
         year_filter=args.years, year_update=args.update)

    # Saturated vapor pressure
    # export_df['esat_min_kPa'] =
    # refet.calcs._sat_vapor_pressure(export_df.Tmin_C)
    # export_df['esat_max_kPa'] =
    # refet.calcs._sat_vapor_pressure(export_df.Tmax_C)
    # export_df['esat_avg_kPa'] =
    # (export_df.esat_min_kPa + export_df.esat_max_kPa)/2
