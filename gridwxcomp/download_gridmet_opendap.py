# -*- coding: utf-8 -*-
"""
Download gridMET climatic time series for multiple variables, e.g. ETr, temp,
wind speed 2m, swrad, etc. 

"""
import argparse
import datetime as dt
import logging
import os
import sys
import timeit
from time import sleep

import pandas as pd
pd.options.display.float_format = '{:,.10f}'.format
import refet
import xarray

# for connection to Earth Engine
# ee.Initialize() 

def download_gridmet_opendap(input_csv, out_folder, year_filter='', 
        update_data=False):

    """
    Download gridMET time series data for multiple climate variables for 
    select gridMET cells as listed in ``input_csv``.

    Arguments:
        input_csv (str): file path of input CSV produced by 
            :mod:`prep_input.py`
        out_folder (str): directory path to save gridmet timeseries CSV files
        year_filter (list): default ''. Single year or range to download
        update_data (bool): default False. Re-download existing data

    Returns:
        None

    Examples:
        Say we wanted to download data for 2016 through 2018, from the command 
        line,

        .. code::
            $ download_gridmet_opendap.py -i merged_input.csv -o gridmet_data -y 2016-2018

        note, "merged_input.csv" should have been created by first running 
        :mod:`prep_input.py`. If the ``[-y, --years]`` option is note given 
        the default behaviour is to download gridMET data from 1979 up through
        yesterday.

        If the data for 2018 has changed since the last run or for debugging
        purposes you can re-download data for all or select years with the
        ``[-u, --update-data]`` option

        .. code::
            $ download_gridmet_opendap.py -i merged_input.csv -o gridmet_data -y 2018 -u

        To download the same gridMET data within Python

        >>> from gridwxcomp import download_gridmet_opendap
        >>> download_gridmet_opendap('merged_input.csv',
                'gridmet_data',
                '2016-2018'
            )

        Running :func:`download_gridmet_opendap.py` also updates the CSV file
        produced from :mod:`prep_input.py` to include file paths to gridMET
        time series files that are paired with climate stations. 

    """
    if not os.path.exists(out_folder):
        logging.info('\nCreating output folder: {}'.format(out_folder))
        os.makedirs(out_folder)

    # Input .csv containing GRIDMET_ID, LAT, LON
    input_df = pd.read_csv(input_csv)

    # Specify column order for output .csv Variables:
    output_order = ['date', 'year', 'month', 'day', 'centroid_lat',
                    'centroid_lon', 'elev_m', 'u2_ms', 'tmin_c', 'tmax_c',
                    'srad_wm2', 'ea_kpa', 'prcp_mm', 'etr_mm', 'eto_mm']

    opendap_url = 'http://thredds.northwestknowledge.net:8080/thredds/dodsC'
    elev_nc = '{}/{}'.format(opendap_url, '/MET/elev/metdata_elevationdata.nc')
    params = {
        'etr': {
            'nc': 'agg_met_etr_1979_CurrentYear_CONUS',
            'var': 'daily_mean_reference_evapotranspiration_alfalfa',
            'col': 'etr_mm'},
        'pet': {
            'nc': 'agg_met_pet_1979_CurrentYear_CONUS',
            'var': 'daily_mean_reference_evapotranspiration_grass',
            'col': 'eto_mm'},
        'pr': {
            'nc': 'agg_met_pr_1979_CurrentYear_CONUS',
            'var': 'precipitation_amount',
            'col': 'prcp_mm'},
        'sph': {
            'nc': 'agg_met_sph_1979_CurrentYear_CONUS',
            'var': 'daily_mean_specific_humidity',
            'col': 'q_kgkg'},
        'srad': {
            'nc': 'agg_met_srad_1979_CurrentYear_CONUS',
            'var': 'daily_mean_shortwave_radiation_at_surface',
            'col': 'srad_wm2'},
        'vs': {
            'nc': 'agg_met_vs_1979_CurrentYear_CONUS',
            'var': 'daily_mean_wind_speed',
            'col': 'u10_ms'},
        'tmmx': {
            'nc': 'agg_met_tmmx_1979_CurrentYear_CONUS',
            'var': 'daily_maximum_temperature',
            'col': 'tmax_k'},
        'tmmn': {
            'nc': 'agg_met_tmmn_1979_CurrentYear_CONUS',
            'var': 'daily_minimum_temperature',
            'col': 'tmin_k'},
    }
    
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

    # Loop through dataframe row by row and grab desired met data from
    # GRID collections based on Lat/Lon and Start/End dates
    for index, row in input_df.iterrows():
        start_time = timeit.default_timer()
        # Reset original_df
        original_df = None
        export_df = None

        GRIDMET_ID_str = str(row.GRIDMET_ID)
        logging.info('Processing GRIDMET ID: {}'.format(GRIDMET_ID_str))

        output_name = 'gridmet_historical_' + GRIDMET_ID_str + '.csv'
        output_file = os.path.join(out_folder, output_name)

        if os.path.isfile(output_file):
            logging.info('{} exists. Checking for missing data.'.format(
                output_name))
            original_df = pd.read_csv(output_file, parse_dates=True)
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
        if not missing_dates and not update_data:
            logging.info('No missing data found. Skipping')
            # Add gridMET file path to input table if not already there
            input_df.loc[input_df.GRIDMET_ID == row.GRIDMET_ID,\
                'GRIDMET_FILE_PATH'] = os.path.abspath(output_file)
            input_df.to_csv(input_csv, index=False)
            continue
        # for option to redownload all data in given year range
        elif not missing_dates and update_data:
            if min(date_list.year) == max(date_list.year):
                yr_rng = min(date_list.year)
            else:
                yr_rng = '{}-{}'.format(
                        min(date_list.year), max(date_list.year))
            logging.info('Updating data for years: {}'.format(yr_rng))
            start_date = min(date_list)
            end_date = max(date_list)
            # missing years are all years for updating
            missing_years = sorted(list(set(date_list.year)))

        # only download years with missing data
        else:
            # Min and Max of Missing Dates (Start: Inclusive; End: Exclusive)
            start_date = min(missing_dates)
            end_date = max(missing_dates)
    
            missing_years = []
            for date in missing_dates:
                missing_years = missing_years + [date.year]
            missing_years = sorted(list(set(missing_years)))
        logging.debug('  Missing Years: {}'.format(
            ', '.join(map(str, missing_years))))

        # Loop through ee pull by year (max 5000 records for getInfo())
        # Append each new year on end of dataframe
        # Process dataframe units and output after
        start_date_yr = start_date.year    
            
        # Calculate out grid cell centroid
        # gridMET elevation asset lower left corner coordinates
        gridmet_lon = -124.78749996666667
        gridmet_lat = 25.04583333333334
        gridmet_cs = 0.041666666666666664
        gridcell_lat = (int(abs(row.LAT - gridmet_lat) / gridmet_cs) * gridmet_cs + 
                        gridmet_lat + 0.5 * gridmet_cs)
        gridcell_lon = (int(abs(row.LON - gridmet_lon) / gridmet_cs) * gridmet_cs + 
                        gridmet_lon + 0.5 * gridmet_cs)
        gridmet_crs = 'EPSG:4326'
        gridmet_geo = [gridmet_cs, 0, gridmet_lon, 
                       0, -gridmet_cs, gridmet_lat + gridmet_cs * 585]
        logging.debug('  Latitude:  {}'.format(gridcell_lat))
        logging.debug('  Longitude: {}'.format(gridcell_lon))

        # OpenDAP call for elevation
        elev_ds = xarray.open_dataset(elev_nc)\
            .sel(lon=gridcell_lon, lat=gridcell_lat, method='nearest')
        elev = elev_ds['elevation'].values[0]
        logging.debug('  Elevation: {}'.format(elev))
        
        # OpenDAP call for each variable
        met_df_list = []
        for met_name in params.keys():
            logging.debug('  Variable: {}'.format(met_name))
            # Pulling the full time series then filtering later seems faster than selecting here
            # # day=pd.date_range(start=start_date, end=end_date), 
            met_nc = '{}/{}.nc'.format(opendap_url, params[met_name]['nc'])
            met_ds = xarray.open_dataset(met_nc)\
                .sel(lon=gridcell_lon, lat=gridcell_lat, method='nearest')\
                .drop(['crs', 'lat', 'lon'])\
                .rename({params[met_name]['var']: params[met_name]['col'],
                         'day': 'date'})
            met_df = met_ds.to_dataframe()
            # logging.debug(met_df.head())
            met_df_list.append(met_df)
            # print(met_df.head())

        # This might need to be a merge call if the indices don't match
        export_df = pd.concat(met_df_list, axis=1)
        logging.debug(export_df.head())
        # End OpenDAP stuff

        # If export_df is None (skip to next ID)
        if export_df is None:
            continue
        # Reset Index
        export_df = export_df.reset_index(drop=False)

        # Convert dateNum to datetime and create Year, Month, Day, DOY variables
        export_df.date = pd.to_datetime(export_df.date.astype(str),
                                        format='%Y-%m-%d')
        export_df['year'] = export_df['date'].dt.year
        export_df['month'] = export_df['date'].dt.month
        export_df['day'] = export_df['date'].dt.day
        # export_df['DOY'] = export_df['Date'].dt.dayofyear
        # Format Date for export
        export_df['date'] = export_df.date.apply(lambda x: x.strftime(
            '%Y-%m-%d'))

        # CGM - This is needed here since filtering by date on the OpenDAP call was really slow
        # Should the dataframe be filtered to missing years or based on start/end date?
        # export_df = export_df.loc[(export_df['date'] >= start_date.strftime('%Y-%m-%d')) & 
        #                           (export_df['date'] <= end_date.strftime('%Y-%m-%d'))]
        export_df = export_df.loc[export_df['year'].isin(missing_years)]
        # export_df = export_df.reset_index(drop=False)
        
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
        export_df['pair_kpa'] = refet.calcs._air_pressure(
            export_df.elev_m, method='asce')

        # actual vapor pressure (kg/kg) using refet module
        export_df['ea_kpa'] = refet.calcs._actual_vapor_pressure(
            export_df.q_kgkg, export_df.pair_kpa)

        # Unit Conversions
        export_df['tmax_c'] = export_df.tmax_k - 273.15  # K to C
        export_df['tmin_c'] = export_df.tmin_k - 273.15  # K to C
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
        export_df.to_csv(output_file, columns=output_order, index=False, 
                         float_format='%.10f')
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
    Command line usage of download_gridmet_opendap.py for downloading
    gridMET time series of several climatic variables using OpenDAP .
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
        '-u','--update-data', required=False, default=False, 
        action='store_true', help='Flag to re-download existing data')
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

    download_gridmet_opendap(input_csv=args.input, out_folder=args.out_dir,
         year_filter=args.years, update_data=args.update_data)

    # Saturated vapor pressure
    # export_df['esat_min_kPa'] =
    # refet.calcs._sat_vapor_pressure(export_df.Tmin_C)
    # export_df['esat_max_kPa'] =
    # refet.calcs._sat_vapor_pressure(export_df.Tmax_C)
    # export_df['esat_avg_kPa'] =
    # (export_df.esat_min_kPa + export_df.esat_max_kPa)/2
