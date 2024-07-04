# -*- coding: utf-8 -*-
"""
This module has tools to read a CSV of climate station metadata information and
verify it has the contents necessary to proceed with the later steps. The
output from this module will be the in a standardized format that is used by
the :mod:`gridwxcomp.ee_download` and :mod:`gridwxcomp.calc_bias_ratios`
modules for the main bias correction workflows.
"""
import ee
import os

import numpy as np
import pandas as pd                                                             
from pathlib import Path
from gridwxcomp.util import read_config, reproject_crs_for_point


def _read_station_list(station_path):
    """
    Helper function that reads station list CSV file and return modified 
    version as a :obj:`Pandas.DataFrame` that includes file paths to each 
    station time series file. Renames some columns for consistency with other 
    ``gridwxcomp`` functions and scripts.

    Arguments:
        station_path (str): path to CSV file containing list of climate
            stations that will later be used to calculate monthly
            bias rations to gridded data.

    Returns:
        station_list (:class:`pandas.DataFrame`): ``Pandas.DataFrame`` that
            contains station name, latitude, longitude, and others for each
            climate station.

    """

    station_list = pd.read_csv(station_path)
    # mandatory columns 
    need_cols = ['Latitude', 'Longitude', 'Filename', 'Station']

    # make sure mandatory columns exist else abort
    station_cols = station_list.columns
    if not set(need_cols).issubset(set(station_cols)):
        err_msg = ('One or more of the mandatory columns is missing from the station input file, it must contain:',
                   ', '.join(c for c in need_cols))
        raise ValueError(err_msg)

    station_list.rename(
            columns={
                'Latitude': 'STATION_LAT',
                'Longitude': 'STATION_LON',
                'Elev_m': 'STATION_ELEV_M',
                'Elev_FT': 'STATION_ELEV_FT',
                'Station': 'STATION_ID',
                'Filename': 'STATION_FILE_PATH'}, inplace=True)

    # get station name only for matching to file name without extension
    station_list.STATION_FILE_PATH = station_list.STATION_FILE_PATH.str.split('.').str.get(0)

    # look at path for station CSV, look for time series files in same directory
    station_path_tuple = os.path.split(station_path)
    path_root = station_path_tuple[0]
    file_name = station_path_tuple[1]

    # look in parent directory that contains station CSV file
    if path_root != '' and file_name != '':
        file_names = os.listdir(path_root)
    # if station CSV file is in cwd look there
    else:
        file_names = os.listdir(os.getcwd())
    # match station name with time series Excel files full path,
    # assumes no other files in the directory have station names in their name
    # will accept files of any extension, e.g. xlx, csv, txt
    for i, station in enumerate(station_list.STATION_FILE_PATH):
        try:
            match = [s for s in file_names if station in s][0]
        except:
            match = None
        if match:
            station_list.loc[station_list.STATION_FILE_PATH == station, 'STATION_FILE_PATH'] = \
                os.path.abspath(os.path.join(path_root, match))
        else:
            missing_station = station_list.iloc[i]['STATION_ID']
            print('WARNING: no file was found that matches station: ', missing_station, '\nin directory: ',
                  os.path.abspath(path_root), '\nskipping.\n')
            continue

    return station_list


def prep_metadata(station_path, config_path, grid_name, 
        out_path='formatted_input.csv'):
    """
    Read list of climate stations in metadata and verify all needed parameters
    exist. An output CSV file is saved that will be the formatted in a way that
    is standardized for the variables that are needed by the subsequent 
    Earth Engine download and bias calculation modules. 

    Station time series files must be in the same directory as the main input
    to this function, i.e., the `station_path` metadata file.

    Arguments:
        station_path (str): path to CSV file containing metadata of climate
            stations that will later be used to calculate bias ratios to 
            the gridded dataset.
        config_path (str): path to config file containing projection info
        grid_name (str): name of the gridded dataset that is being used
            for comparison against observed data.
        out_path (str): path to save output CSV, default is to save as 
            'merged_input.csv' to current working directory.


    Returns:
        None

    Example:

        >>> from gridwxcomp import prep_metadata
        >>> prep_metadata('example_metadata.txt','outfile.csv')
        
        outfile.csv will be created containing station and corresponding
        gridded data. This file is later used as input for
        :mod:`gridwxcomp.ee_download` and
        :mod:`gridwxcomp.calc_bias_ratios`.

    Important:
        Make sure the following column headers exist in your input station 
        metadata file (``station_path``) and are spelled exactly:

          * Latitude
          * Longitude
          * Station
          * Filename

        Also, the "Filename" column should match the names of the climate time
        series files that should be in the same directory as the station
        metadata file. For example, if one of the time series files is named
        "Bluebell_daily_data.csv" then the following are permissiable entries
        as the "Filename": "Bluebell_daily_data" or "Bluebell_daily_data.csv".
        
    Raises:
        ValueError: if one or more of the following mandatory columns are 
            missing from the input CSV file (``station_path`` parameter): 
            'Longitude', 'Latitude', 'Station', or 'Filename'.   
    """

    # Create parent directories if necessary
    path_root = Path(out_path).parent
    if not path_root.is_dir():
        print('The directory: ', path_root.absolute(), ' does not exist, creating directory')
        os.makedirs(path_root)

    print('station list CSV: ', os.path.abspath(station_path))
    print('merged CSV will be saved to: ', os.path.abspath(out_path))

    config = read_config(config_path)

    stations = _read_station_list(station_path)
    stations[f'GRID_ID'] = f'{grid_name}_' + stations['STATION_ID']

    if 'ELEV_M' in stations.columns:
        stations['ELEV_FT'] = stations.ELEV_M * 3.28084  # m to ft


    # Add WGS84 projection columns for earth engine requests
    temp_proj_df = stations[['STATION_LAT', 'STATION_LON']].copy(deep=True)
    temp_proj_df['STATION_LAT_WGS84'] = np.nan
    temp_proj_df['STATION_LON_WGS84'] = np.nan

    for index, row in temp_proj_df.iterrows():
        (temp_proj_df.loc[index, 'STATION_LON_WGS84'],
         temp_proj_df.loc[index, 'STATION_LAT_WGS84']) =\
            reproject_crs_for_point(
            row['STATION_LON'], row['STATION_LAT'],
            config['input_data_projection'], 'EPSG:4326')

    stations['STATION_LAT_WGS84'] = temp_proj_df['STATION_LAT_WGS84']
    stations['STATION_LON_WGS84'] = temp_proj_df['STATION_LON_WGS84']

    # save CSV
    stations.to_csv(out_path, index=False)


