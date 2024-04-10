# -*- coding: utf-8 -*-
"""
Read a CSV of climate station information and verifies it has the contents necessary to proceed
with the later steps. The output from this module will be the input for main bias correction workflow.

TODO: add logging

"""

import os
import pandas as pd                                                             
from pathlib import Path


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


def prep_metadata(station_path, grid_name, out_path='formatted_input.csv'):
    """
    Read list of climate stations in metadata and verify all needed parameters exist

    Station time series files must be in the same directory as `station_path`
    metadata file.

    Arguments:
        station_path (str): path to CSV file containing metadata of climate
            stations that will later be used to calculate bias ratios to 
            the gridded dataset.
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
        :mod:`gridwxcomp.download_gridmet_opendap` and
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

    stations = _read_station_list(station_path)
    stations[f'GRID_ID'] = f'{grid_name}_' + stations['STATION_ID']

    if 'ELEV_M' in stations.columns:
        stations['ELEV_FT'] = stations.ELEV_M * 3.28084  # m to ft

    # save CSV 
    stations.to_csv(out_path, index=False)


def get_subgrid_bounds(in_path, grid_res, rounding_decimals, buffer=0):
    """
    Calculate bounding box for spatial interpolation grid from
    output of prep_metadata.prep_metadata()

    Arguments:
        in_path (str): path to comprehensive summary file containing
            monthly bias ratios, created by :func:`gridwxcomp.calc_bias_ratios`.
        grid_res (float): the resolution of the fishnet grid, in the units projection the input data is in
        rounding_decimals (int): number of decimal places to round the subgrid bounds to
        buffer (int): number of grid cells to expand the rectangular extent
            of the subgrid fishnet.

    Returns:
        bounds (dict): dictionary with coordinates that
            define the outer bounds of the subgrid fishnet in the format
            (min long, max long, min lat, max lat)

    Raises:
        FileNotFoundError: if input summary CSV file is not found.

    Note:
        By expanding the grid to a larger area encompassing the climate
        stations of interest :func:`interpolate` can be used to extrapolate
        passed the bounds of the outer station locations.

    """
    if not os.path.isfile(in_path):
        raise FileNotFoundError('Input metadata file given was invalid or not found')


    in_df = pd.read_csv(in_path)

    lons = in_df.sort_values('STATION_LON')['STATION_LON'].values
    lats = in_df.sort_values('STATION_LAT')['STATION_LAT'].values

    if rounding_decimals == 0:
        # likely projection is in meters, use modulo to ensure grid alignment
        lon_min = (lons[0] - (lons[0] % grid_res)) - grid_res
        lon_max = (lons[0] + (lons[0] % grid_res)) + grid_res
        lat_min = (lons[0] - (lats[0] % grid_res)) - grid_res
        lat_max = (lons[0] + (lats[0] % grid_res)) + grid_res
    else:
        lon_min = lons[0] - grid_res
        lon_max = lons[-1] + grid_res
        lat_min = lats[0] - grid_res
        lat_max = lats[-1] + grid_res

    # expand bounding extent based on buffer cells
    lon_min -= grid_res * buffer
    lat_min -= grid_res * buffer
    lon_max += grid_res * buffer
    lat_max += grid_res * buffer

    return {'xmin': round(lon_min, rounding_decimals), 'xmax': round(lon_max, rounding_decimals),
            'ymin': round(lat_min, rounding_decimals), 'ymax': round(lat_max, rounding_decimals)}
