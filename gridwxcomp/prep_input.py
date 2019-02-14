# -*- coding: utf-8 -*-
"""
Read a CSV of climate station information and match each with nearest gridMET
cell information. Produce a CSV file that will be used for input for main
bias correction workflow.

Todo:
    *  add logging 
"""

import os                                                                       
import argparse                                                                 
import logging
import pkg_resources
                                                    
import pandas as pd                                                             
import numpy as np        
from scipy import spatial



def main(station_file, out_path, gridmet_meta_file):
    """
    Take list of climate stations and merge each with overlapping gridMET cell
    information, write new CSV for next step in bias correction workflow.

    Arguments:
        station_file (str): path to CSV file containing list of climate
            stations that will later be used to calculate monthly
            bias ratios to GridMET reference ET.
        gridmet_meta_file (str): path to metadata CSV file that contains
            all gridMET cells for the contiguous United States. Can be
            found at ``etr-biascorrect/gridmet_cell_data.csv``.
        out_path (str): path to save output CSV, default is to save as 
            "merged_input.csv" to current working directory if not passed
            at command line to script.

    Example:
        From the command line interface,

        .. code::
            $ python prep_input.py -i ETrBias_DataPackage/Station_Data.txt -g gridmet_cell_data.csv

        The result is "merged_input.csv" being created in the working 
        directory which contains metadata from climate staions as well as the 
        lat, long, and gridMET ID of the nearest gridMET cell centroid.
        This file is used as input to :mod:`download_gridmet_ee.py`
        followed by :mod:`calc_bias_ratios.py`.

    """

    # station info with overlapping gridMET and save CSV
    prep_input(
        station_file, 
        out_path,
        gridmet_meta_path=gridmet_meta_file 
    )

def gridMET_centroid(lat,lon):
    """
    Calculate the nearest centroid lattitude and longitude for an arbitrary
    coordinate. Used for finding closest neighboring grdiMET cell to a climate 
    station.
    
    Arguments:
        lat (float): decimal degree latitude of location
        lon (float): decimal degree longitude of location 
            
    Returns:
        gridcell_lat,gridcell_lon (tuple): tuple of latitude and 
            longitude of nearest gridMET cell centroid location.
    """
    gridmet_lon = -124.78749996666667
    gridmet_lat = 25.04583333333334
    gridmet_cs = 0.041666666666666664
    gridcell_lat = int(
        abs(lat - gridmet_lat) / gridmet_cs) * gridmet_cs +\
        gridmet_lat + gridmet_cs/2
    gridcell_lon = int(
        abs(lon - gridmet_lon) / gridmet_cs) * gridmet_cs +\
        gridmet_lon + gridmet_cs/2
    
    return gridcell_lat, gridcell_lon

def _read_station_list(station_path):
    """
    Helper function that reads station list CSV file and return modified 
    version as a :obj:`Pandas.DataFrame` that includes file paths to each 
    station time series file. Renames some columns for consistency with other 
    ``etr-biascorrect`` functions and scripts. 

    Arguments:
        station_path (str): path to CSV file containing list of climate
            stations that will later be used to calculate monthly
            bias rations to GridMET reference ET.

    Returns:
        station_list (:class:`pandas.DataFrame`): ``Pandas.DataFrame`` that
            contains station ID, lattitude, longitude, elevation, and others 
            for each climate station.

    """

    station_list = pd.read_csv(station_path)
    cols = [
            'FID',
            'LATDECDEG',
            'LONGDECDEG',
            'Elev_m',
            'FileName',
            'Station_ID',
            'Elev_FT',
            'State',
            'Source',
            'Station',
            'Website',
            'Comments',
            'Irrigation'
            ]
    station_list = station_list[cols]
    station_list.rename(
            columns={
                'LATDECDEG':'STATION_LAT',
                'LONGDECDEG':'STATION_LON',
                'Elev_m':'STATION_ELEV_M',
                'Elev_FT':'STATION_ELEV_FT',
                'Station_ID':'STATION_ID',
                'FileName':'STATION_FILE_PATH'},
            inplace=True
            )
    # get station name only for matching to file name
    station_list.STATION_FILE_PATH =\
            station_list.STATION_FILE_PATH.str.split('_').str.get(0)
    # look at path for station CSV, look for time series files in same directory
    station_path_tuple = os.path.split(station_path)
    path_root = station_path_tuple[0]
    file_name = station_path_tuple[1]
    # look in parent directory that contains station CSV file
    if path_root != '' and file_name != '':
        file_names = os.listdir(path_root)   
    # if station CSV file is in same directory look there
    else:
        file_names = os.listdir(os.getcwd())
    # match station name with time series excel files full path,
    # assumes no other files in the directory have station names in their name
    # will accept files of any extension, e.g. xlx, csv, txt
    for station in station_list.STATION_FILE_PATH:
        try:
            match = [s for s in file_names if station in s][0]
        except:
            match = None
        if match:
            station_list.loc[station_list.STATION_FILE_PATH == station,\
                'STATION_FILE_PATH'] = os.path.abspath(
                                       os.path.join(path_root,match))
        else:
            print('No file was found that matches station: ', station,
                    '\nin directory: ', os.path.abspath(path_root),
                    '\nskipping.\n')
            continue

    return station_list

def prep_input(station_path, out_path='merged_input.csv', 
                                                gridmet_meta_path=None):
    """
    Read list of climate stations and match each with its
    closest GridMET cell, save CSV with information from both.

    Arguments:
        station_path (str): path to CSV file containing list of climate
            stations that will later be used to calculate monthly
            bias ratios to GridMET reference ET.

    Keyword Arguments:
        out_path (str): path to save output CSV, default is to save as 
            "merged_input.csv" to current working directory.
        gridmet_meta_path (str): path to metadata CSV file that contains
            all gridMET cells for the contiguous United States. If None
            it is looked for at ``etr-biascorrect/gridmet_cell_data.csv``.

    Returns:
        None

    Example:

        >>> from prep_input import prep_input 
        >>> prep_input(
                'ETrBias_DataPackage/Station_Data.txt',
                'gridmet_cell_data.csv',
                'outfile.csv'
            )
        
        outfile.csv will be created containing station and corresponding
        gridMET cell data. This file is later used as input for 
        :mod:`download_gridmet_ee.py` and :mod:`calc_bias_ratios.py`.
        
    Raises:
        FileNotFoundError: if the ``gridmet_meta_path`` is not passed as a 
        command line argument and it is not in the current working directory
        and named "gridmet_cell_data.csv".

    Note:
        The CSV file that is saved contains latitude, longitude, and elevation
        fields for both the station and nearest gridMET cell centroid. Fields 
        that may refer to both gridMET and station data have prefixes to 
        distinguish, the climate station data are prefixed with "STATION_" and 
        those refering to gridMET have no prefix. Other fields without a 
        prefix are not in all capital letters and refer to the climate station, 
        e.g. "Website". 

    """

    # look for pacakged gridmet_cell_data.csv if path not given
    if not gridmet_meta_path:
        try:
            if pkg_resources.resource_exists('gridwxcomp', 
                    "gridmet_cell_data.csv"):
                gridmet_meta_path = pkg_resources.resource_filename(
                    'gridwxcomp', 
                    "gridmet_cell_data.csv"
                    )
        except:
            gridmet_meta_path = 'gridmet_cell_data.csv'
    if not os.path.exists(gridmet_meta_path):
        raise FileNotFoundError('GridMET file path was not given and '+\
                'gridmet_cell_data.csv was not found in the gridwxcomp '+\
                'install directory. Please assign the path or put '+\
                '"gridmet_cell_data.csv" in the current working directory.\n')

    path_root = os.path.split(os.path.abspath(out_path))[0]
    if not os.path.exists(path_root):
        print(
            'The directory: ', 
            os.path.abspath(path_root),
            '\ndoes not exist, creating directory'
        )
        os.makedirs(path_root)

    print(
          'station list CSV: ',
          os.path.abspath(station_path),
          '\ngridMET cell info CSV: ',
          os.path.abspath(gridmet_meta_path),
          '\nmerged CSV will be saved to: ',
          os.path.abspath(out_path)
    )

    stations = _read_station_list(station_path)
    gridmet_meta = pd.read_csv(gridmet_meta_path)
    gridmet_pts = list(zip(gridmet_meta.LAT,gridmet_meta.LON))
    # scipy KDTree to find nearest neighbor between station and centroids
    tree = spatial.KDTree(gridmet_pts)
    # loop through each station find closest GridMET
    for index, row in stations.iterrows():
        try:
            station_lat = row.STATION_LAT
            station_lon = row.STATION_LON
            pt = np.array([station_lat,station_lon])
            # index of nearest GridMET point, same as GRIDMET_ID
            ind = tree.query(pt)[1]
            stations.loc[index,'GRIDMET_ID'] = ind
        except:
            print('Failed to find matching gridMET info for climate '\
                    +'station with STATION_ID = ', row.STATION_ID,'\n')  
    stations.GRIDMET_ID = stations.GRIDMET_ID.astype(int)
    out_df = stations.merge(gridmet_meta,on='GRIDMET_ID')
    out_df['ELEV_FT'] = out_df.ELEV_M * 3.28084 # m to ft
    out_df = out_df.reindex(
             columns=['GRIDMET_ID',
                      'LAT',
                      'LON',
                      'ELEV_M',
                      'ELEV_FT',
                      'STATION_ID',
                      'FID',
                      'STATION_LAT',
                      'STATION_LON',
                      'STATION_ELEV_M',
                      'STATION_ELEV_FT',
                      'STATION_FILE_PATH',
                      'State',
                      'Source',
                      'Station',
                      'Website',
                      'Comments',
                      'Irrigation'
            ])
    # save CSV 
    out_df.to_csv(out_path, index=False)

def arg_parse():
    """
    Command line usage for prep_input.py for merging climate station 
    and corresponding gridMET data into a single table (CSV file).
    The CSV file produced is used as input to download_gridmet_ee.py 
    and calc_bias_ratios.py.
    """
    parser = argparse.ArgumentParser(
        description=arg_parse.__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    optional = parser._action_groups.pop() # optionals listed second
    required = parser.add_argument_group('required arguments')
    required.add_argument(
        '-i', '--input', metavar='PATH', required=True,
        help='Input CSV file of climate stations')
    optional.add_argument(
        '-g', '--gridmet', metavar='PATH', required=False,
        default='gridmet_cell_data.csv',
        help='GridMET master CSV file with cell data, packaged with '+\
             'etr-biascorrect at etr-biascorrect/gridmet_cell_data.csv '+\
             'if not given it needs to be located in the currect directory')
    optional.add_argument(
        '-o', '--out', metavar='PATH', required=False, 
        default='merged_input.csv',
        help='Optional output path for CSV with merged climate/gridMET data')
    parser._action_groups.append(optional)# to avoid optionals listed first
#    parser.add_argument(
#        '--debug', default=logging.INFO, const=logging.DEBUG,
#        help='Debug level logging', action="store_const", dest="loglevel")
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = arg_parse()

    main(station_file=args.input, 
         out_path=args.out,
         gridmet_meta_file=args.gridmet
         )
