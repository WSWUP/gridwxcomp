# -*- coding: utf-8 -*-
"""
Read a CSV of climate station information and match each with nearest gridMET
(or other gridded product) cell information. Has routines to build metadata file from arbitrary uniform gridded dataset and ultimately produce a CSV file that 
will be used for input for main bias correction workflow.

Todo:
    *  add logging 
"""

import argparse         
import fiona
import logging
import os                                                                       
from pathlib import Path
                                                    
import pandas as pd                                                             
import numpy as np        
from scipy import spatial
from shapely.geometry import Polygon

# allows for CL script usage if gridwxcomp not installed
try:
    from .util import get_gridmet_meta_csv
except:
    from util import get_gridmet_meta_csv

def main(station_file, out_path, grid_meta_file, grid_path, grid_id_name,
        grid_data_dir):
    """
    Take list of climate stations and merge each with overlapping gridMET cell
    information, write new CSV for next step in bias correction workflow.

    Arguments:
        station_file (str): path to CSV file containing list of climate
            stations that will later be used to calculate monthly
            bias ratios to gridMET reference ET.
        grid_meta_file (str): path to metadata CSV file that contains
            all gridMET cells for the contiguous United States. Can be
            found at ``gridwxcomp/gridmet_cell_data.csv``.
        out_path (str or None): path to save output CSV, default is to save 
            as "merged_input.csv" to current working directory if not passed
            at command line to script.
        grid_path (str): path to grid vector file if not using gridMET.
        grid_id_name (str): name of gridcell identifier present in grid,
            ID data values should be integers, only if using custom grid.
        grid_data_dir (str): directory that contains grid time series files,
            each file should have the integer grid ID value in its name and 
            should be in CSV format. Only used when gridded time series data 
            already exists on disk, i.e. when not using gridMET as the gridded
            data.

    Example:
        From the command line interface within the ``gridwxcomp/gridwxcomp``
        directory (or replace input path with correct path),

        .. code-block:: sh

            $ gridwxcomp prep-input <station_metadata>  

        where ``station_metadata`` is a file containing metadata of climate 
        stations built from `PyWeatherQAQC <https://github.com/WSWUP/pyWeatherQAQC>`_.
        The file should be in CSV format and contain at least these four columns:
        * Latitude
        * Longitude
        * Station
        * Filename 

        The result is "merged_input.csv" being created in the working 
        directory which contains metadata from climate staions as well as the 
        lat, long, and gridMET ID of the nearest gridMET cell centroid.
        This file is used as input to :mod:`gridwxcomp.download_gridmet_opendap`
        followed by :mod:`gridwxcomp.calc_bias_ratios`.

    """

    # station info with overlapping gridMET and save CSV
    prep_input(
        station_file, 
        out_path,
        grid_meta_path=gridmet_meta_file,
        grid_path=grid_path,
        grid_id_name=grid_id_name,
        grid_data_dir=grid_data_dir
    )


def _read_station_list(station_path):
    """
    Helper function that reads station list CSV file and return modified 
    version as a :obj:`Pandas.DataFrame` that includes file paths to each 
    station time series file. Renames some columns for consistency with other 
    ``gridwxcomp`` functions and scripts. 

    Arguments:
        station_path (str): path to CSV file containing list of climate
            stations that will later be used to calculate monthly
            bias rations to GridMET reference ET.

    Returns:
        station_list (:class:`pandas.DataFrame`): ``Pandas.DataFrame`` that
            contains station name, lattitude, longitude, and others for each 
            climate station.

    """

    station_list = pd.read_csv(station_path)
    # mandatory columns 
    need_cols = [
            'Latitude',
            'Longitude',
            'Filename',
            'Station',
            ]

    # make sure mandatory columns exist else abort
    station_cols = station_list.columns
    if not set(need_cols).issubset(set(station_cols)):
        err_msg = ('One or more of the mandatory columns is missing',
            'from the station input file, it must contain:',
            ', '.join(c for c in need_cols))
        raise ValueError(err_msg)

    station_list.rename(
            columns={
                'Latitude':'STATION_LAT',
                'Longitude':'STATION_LON',
                'Elev_m':'STATION_ELEV_M',
                'Elev_FT':'STATION_ELEV_FT',
                'Station':'STATION_ID',
                'Filename':'STATION_FILE_PATH'},
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
    # if station CSV file is in cwd look there
    else:
        file_names = os.listdir(os.getcwd())
    # match station name with time series excel files full path,
    # assumes no other files in the directory have station names in their name
    # will accept files of any extension, e.g. xlx, csv, txt
    for i, station in enumerate(station_list.STATION_FILE_PATH):
        try:
            match = [s for s in file_names if station in s][0]
        except:
            match = None
        if match:
            station_list.loc[station_list.STATION_FILE_PATH == station,\
                'STATION_FILE_PATH'] = os.path.abspath(
                                       os.path.join(path_root,match))
        else:
            missing_station = station_list.iloc[i]['STATION_ID']
            print('WARNING: no file was found that matches station: ', 
                missing_station, '\nin directory: ', 
                os.path.abspath(path_root), '\nskipping.\n'
            )
            continue

    return station_list

def _get_cell_centroid(coords, x_cell_size, y_cell_size):
    """get centroid of gridcell (Shapely Polygon)"""
    poly = Polygon(coords)
    # bounds gives a (minx, miny, maxx, maxy) tuple
    lon_c = poly.bounds[0] + x_cell_size / 2
    lat_c = poly.bounds[1] + y_cell_size / 2
    
    return lat_c, lon_c

def build_grid_meta(grid_path, grid_id_name, out_path=None):
    """
    Build a metadata (CSV) file for an arbitrary georeferenced grid
    (vector of polygons) that represents a master (full extent) grid for
    a corresponding gridded meterological dataset. The output CSV will
    include for each cell at least: an ID and gridcell centroid latitude
    and longitude.

    The grid that is passed should include integer IDs for each cell which
    should be named as `grid_id_name`. Any other attributes will also be saved
    to the output CSV meta file. The grid file should have a coordinate
    reference system in decimal degrees, e.g. WGS 84 - Geographic system.

    If the file already exists at `out_path` this function will **not**
    overwrite it.

    Arguments:
        grid_path (str): path to grid vector file
        grid_id_name (str): name of gridcell identifier present in grid,
            ID data values should be integers.

    Keyword Arguments:
        out_path (str or None): default None. Path to save output metadata
            CSV, if None save to "grid_cell_data.csv" in current directory.

    Returns:
        out_path (:obj:`pathlib.Path`): absolute path to saved gridcell metadata
            CSV file

    """
    # make sure grid file exists
    if grid_path and not Path(grid_path).is_file():
        raise FileNotFoundError('ERROR: Grid file was not found')

    # check output directory
    if out_path is None:
        out_path = Path.cwd() / 'grid_cell_data.csv'
    else:
        out_path = Path(out_path)

    # exit if meta file already exists, do not overwrite
    if Path(out_path).is_file():
        print('{} already exists, it will not be overwritten, skipping\n'.\
            format(out_path)
        )
        out_path = out_path.absolute()
        return out_path

    # create any sub-directories in out_path if they do not exist
    if not out_path.parent.is_dir():
        print(
            '\nOutput directory: {}\ndoes not exist, creating it now.\n'.format(
                out_path.parent.absolute()
            )
        )
        out_path.parent.mkdir(parents=True, exist_ok=True)

    # read grid file, read attributes for each gridcell, write
    print('Extracting attributes from grid shapefile:\n {}\n'
        '\nAnd saving to: {}\n'.format(
            Path(grid_path).absolute(), out_path.absolute()
        )
    )
    grid_meta_df = pd.DataFrame()
    # attributes to NOT write to metadata file, bnds from Shapely
    exclude_attrs = ['left', 'top', 'right', 'bottom', grid_id_name]

    with fiona.open(grid_path, 'r') as source:
        n_cells = len([f for f in source])
        print(
            'Looking up and assigning cell data for', n_cells,
            'gridcells.\n'
        )
        if n_cells >= 10000:
            time_est_min = round(((n_cells // 1000) * 2.5) / 60)
            print(
                'This will take "roughly"', time_est_min, 'minutes.\n'
            )
        for i, feature in enumerate(source):
            coords = feature['geometry']['coordinates'][0]

            if i == 0:
                # read names of any extra cell attributes to save
                extra_attrs =\
                    set(feature['properties'].keys()) - set(exclude_attrs)
                # calculate the X and Y cell size of grid
                X_CS = abs(coords[0][0] - coords[1][0])
                Y_CS = abs(coords[1][1] - coords[2][1])

            lat,lon = _get_cell_centroid(coords, X_CS, Y_CS)
            grid_id = int(feature['properties'][grid_id_name])
            grid_meta_df.loc[grid_id, 'LAT'] = lat
            grid_meta_df.loc[grid_id, 'LON'] = lon
            # add any extra attributes if they exist (e.g. elevation)
            for attr in extra_attrs:
                grid_meta_df.loc[grid_id, attr] = feature['properties'][attr]

    grid_meta_df.sort_index().to_csv(out_path)
    print(
        'Successfully saved gridcell metadata for grid at:\n {}\n'
        '\nto: {}'.format(
            Path(grid_path).absolute(), out_path.absolute()
        )
    )
    out_path = out_path.absolute()

    return out_path


def prep_input(station_path, out_path='merged_input.csv', grid_meta_path=None, 
        grid_path=None, grid_id_name=None, grid_data_dir=None):
    """
    Read list of climate stations and match each with its closest gridcell, 
    save CSV with information from both.

    Station time series files must be in the same directory as `station_path`
    metadata file and each file must begin with the name specifed as 'Station'
    in the station metadata file.

    If using gridded data other than gridMET this function may be used to 
    create a metadata CSV file of cell data for any arbitrary rectangular grid.
    The grid must be passed to ``grid_path`` and it must contain a cell 
    identifier attribute (name must be passed in as ``grid_id_name``) that is 
    an integer which increases monotonically by steps of 1 without gaps, 
    e.g. 1,2,3,4,... although the cell order does not have to follow any rule. 
    For example, the first cell may be in any location and the next may be 
    anywhere and so forth. Also if using a different grid, the time series 
    files associated with the cell IDs that match your station locations 
    should be in the directory given as ``grid_data_dir``.

    Arguments:
        station_path (str): path to CSV file containing metadata of climate
            stations that will later be used to calculate bias ratios to 
            GridMET.

    Keyword Arguments:
        out_path (str): path to save output CSV, default is to save as 
            'merged_input.csv' to current working directory.
        grid_meta_path (str or None): default None. Path to save grid metadata
            CSV, if None save to "grid_cell_data.csv" in current directory. This
            is only used if working with a user provided gridded dataset, i.e. 
            if `grid_path` and `grid_id_name` (if creating a new grid meta file 
            are given. 
        grid_path (str): path to grid vector file if not using gridMET.
        grid_id_name (str): name of gridcell identifier present in grid,
            ID data values should be integers.
        grid_data_dir (str): directory that contains grid time series files,
            each file should have the integer grid ID value in its name and 
            should be in CSV format. Only used when gridded time series data 
            already exists on disk.

    Returns:
        None

    Example:

        >>> from gridwxcomp import prep_input 
        >>> prep_input('gridwxcomp/example_data/Station_Data.txt','outfile.csv')
        
        outfile.csv will be created containing station and corresponding
        gridMET cell data. This file is later used as input for
        :mod:`gridwxcomp.download_gridmet_opendap` and
        :mod:`gridwxcomp.calc_bias_ratios`.

    Important:
        Make sure the following column headers exist in your input station 
        metadata file (``station_path``) and are spelled exactly:

          * Latitude
          * Longitude
          * Station
          * Filename

        Also, the "Filename" column should have names that begin with the names
        of climate time series files that should be in the same directory as
        the station metadata file. For example, if one of the time series files
        is named "Bluebell_daily_data.csv" then all of the following are
        permissiable entries as the "Filename": "Bluebell", "Bluebell_daily",
        or "Bluebell_daily_data.csv".
        
    Raises:
        FileNotFoundError: if the ``grid_meta_path`` is not passed as a 
            command line argument and it is not in the current working directory
            and named "gridmet_cell_data.csv" (i.e. if other grid data is not 
            given) and if ``gridwxcomp`` was not installed to the user's PATH, 
            i.e. pip or python setup.py install.
        ValueError: if one or more of the following mandatory columns are 
            missing from the input CSV file (``station_path`` parameter): 
            'Longitude', 'Latitude', 'Station', or 'Filename'.   

    Note:
        If climate station time series files do **NOT** follow the format 
        created by `pyWeatherQAQC <https://github.com/WSWUP/pyWeatherQAQC>`_ 
        i.e. microsoft excel files with data stored in a tab named 'Corrected Data'.
        Then station files should be in text (CSV) format with a column
        containing datetime strings e.g. '12/01/2018', that are 
        able to be parsed by Pandas. The CSV file produced by :func:`prep_input`
        contains latitude, longitude, and other fields for both the station 
        and nearest gridcell centroid coordinates. Fields that may refer 
        to both grid and station data have prefixes to distinguish, the 
        climate station data are prefixed with 'STATION' and those refering 
        to grid data have no prefix. 

    """
    # for building from user's grid (not gridMET)
    if grid_path:
        grid_meta_path = build_grid_meta(
            grid_path, grid_id_name, out_path=grid_meta_path)
    # otherwise assume gridMET data
    else:
        # look for pacakged gridmet_cell_data.csv if path not given
        grid_meta_path = get_gridmet_meta_csv(
                gridmet_meta_path=grid_meta_path)
        grid_id_name = 'GRIDMET_ID'

    path_root = Path(out_path).parent
    if not path_root.is_dir():
        print(
            'The directory: ', 
            path_root.absolute(),
            '\ndoes not exist, creating directory'
        )
        os.makedirs(path_root)

    print(
          'station list CSV: ',
          os.path.abspath(station_path),
          '\ngridcell meta info CSV: ',
          os.path.abspath(grid_meta_path)
    )

    if grid_data_dir is not None:
        print('grid data files in dir: ', Path(grid_data_dir).absolute())

    print(
          'merged CSV will be saved to: ',
          os.path.abspath(out_path)
    )

    stations = _read_station_list(station_path)
    grid_meta = pd.read_csv(grid_meta_path, index_col=grid_id_name)
    # make sure gridcell integer index ID is sorted ascending
    grid_meta.sort_index(inplace=True)
    # array of grid lat long for searching with KD tree
    gridmet_pts = list(zip(grid_meta.LAT,grid_meta.LON))
    grid_id_start_int = grid_meta.index[0]
    # scipy KDTree to find nearest neighbor between station and centroids
    tree = spatial.KDTree(gridmet_pts)
    # loop through each station find closest gridcell
    for index, row in stations.iterrows():
        try:
            station_lat = row.STATION_LAT
            station_lon = row.STATION_LON
            pt = np.array([station_lat,station_lon])
            # index of nearest gridcell, same as grid_id because starts at 
            # first integer grid ID which should be presorted
            ind = tree.query(pt)[1] + grid_id_start_int
            stations.loc[index, grid_id_name] = ind 
        except:
            print('Failed to find matching gridcell info for climate '\
                    +'station with STATION_ID = ', row.STATION_ID,'\n')  
    stations[grid_id_name] = stations[grid_id_name].astype(int)
    out_df = stations.merge(grid_meta, on=grid_id_name)
    if 'ELEV_M' in out_df.columns:
        out_df['ELEV_FT'] = out_df.ELEV_M * 3.28084 # m to ft

    # if grid_data_dir is given look their to add grid_file_paths
    if grid_data_dir is not None:
        out_df = find_grid_files(out_df, grid_data_dir, grid_id_name)
    # save CSV 
    out_df.to_csv(out_path, index=False)

def find_grid_files(merged_df, grid_data_dir, grid_id_name):
    """
    Given a directory that contains time series climate files from a gridded
    dataset, find all file paths that match the grid ID for each 
    station-gridcell pair.

    This is used when the gridded data has already been downloaded or otherwise
    exists on the file system. For example, when using gridded products other 
    than gridMET or when creating a new gridMET merged input file that reuses
    predownloaded gridMET data with a different set of climate station data.

    Arguments:
        merged_df (:obj:`pandas.DataFrame`): dataframe created by 
            :func:`prep_input` that contains the integer grid ID or identifier 
            for each grid cell that needs to be matched to overlapping stations.
        grid_data_dir (str): directory that contains grid time series files,
            each file should have the integer grid ID value in its name and 
            should be in CSV format. 
        grid_id_name (str): name of gridcell integer identifier present in grid.

    Returns:
        merged_df (:obj:`pandas.DataFrame`): dataframe passed in as an argument
            modified to include the column 'GRID_FILE_PATH' which contains full
            paths to each grid time series file that corresponds with 
            overlapping station locations. 
    """
    if grid_data_dir is not None:
        grid_data_dir = Path(grid_data_dir).absolute()

    if not Path(grid_data_dir).is_dir():
        print('ERROR: the directory given for gridded time series was not' 
            'found at:\n{}\nThe merged input file will not include paths to'
            'paired gridded time series, fix and rerun this routine.'\
                .format(grid_data_dir)
        )
        return merged_df

    if grid_id_name is None:
        print('ERROR: please specify the name of the integer grid ID used'
            'in your grid and gridded time series names\nThe merged input'
            'file will not include paths to paired gridded time series,\n'
            'fix and rerun this routine.'
        )
        return merged_df

    grid_data_files = [str(f) for f in grid_data_dir.glob('*') if f.is_file()]

    for index, row in merged_df.iterrows():
        grid_id = row[grid_id_name]
        try:
            match = [f for f in grid_data_files if str(grid_id) in f][0]
        except:
            match = None
        if match:
            merged_df.loc[index, 'GRID_FILE_PATH'] = match
        else:
            print('WARNING: no grid file was found for {} = {}'
                '\nin directory: {} \nskipping.\n'.format(
                    grid_id_name, grid_id, grid_data_dir
                )
            )
            continue
            
    return merged_df

def arg_parse():
    """
    Command line usage for prep_input.py for merging climate station 
    and corresponding gridMET data into a single table (CSV file).
    The CSV file produced is used as input to download_gridmet_opendap.py 
    and calc_bias_ratios.py.
    """
    parser = argparse.ArgumentParser(
        description=arg_parse.__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    optional = parser._action_groups.pop() # optionals listed second
    required = parser.add_argument_group('required arguments')
    required.add_argument(
        '-i', '--input', metavar='PATH', required=True,
        help='Climate station metadata CSV file')
    optional.add_argument(
        '-g', '--grid-meta', metavar='PATH', required=False,
        default=None,
        help='GridMET metadata CSV file with cell data, packaged with '+\
             'gridwxcomp and automatically found if pip was used to install '+\
             'if not given it needs to be located in the current directory')
    optional.add_argument(
        '-o', '--out', metavar='PATH', required=False, 
        default='merged_input.csv',
        help='Optional output path for CSV with merged climate/gridMET data')
    optional.add_argument(
        '--grid-path', metavar='PATH', required=False, default=None,
        help='Path to grid shapefile if not using gridMET')
    optional.add_argument(
        '--grid-id-name', metavar='STR', required=False, default=None,
        help='Name of gridcell integer ID used in grid if not using gridMET')
    optional.add_argument(
        '--grid-data-dir', metavar='PATH', required=False, default=None,
        help='Path to gridded time series files if not using gridMET')
    parser._action_groups.append(optional)# to avoid optionals listed first
#    parser.add_argument(
#        '--debug', default=logging.INFO, const=logging.DEBUG,
#        help='Debug level logging', action="store_const", dest="loglevel")
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = arg_parse()

    main(
        station_file=args.input, 
        out_path=args.out,
        grid_meta_file=args.grid_meta,
        grid_path=args.grid_path,
        grid_id_name=args.grid_id_name,
        grid_data_dir=args.grid_data_dir
    )
