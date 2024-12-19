# -*- coding: utf-8 -*-
"""
Utility functions or classes for ``gridwxcomp`` package
"""
import configparser as cp
import ee
import numpy as np
import os
import pandas as pd
import pathlib as pl
import pkg_resources
import pyproj


def affine_transform(img):
    """
    Get the affine transform of the image as an EE object

    Arguments:
        img: ee.Image object

    Returns
        ee.List object

    """
    return ee.List(ee.Dictionary(
        ee.Algorithms.Describe(img.projection())).get('transform'))


def parse_yr_filter(dt_df, years, label):
    """
    Parse string year filter and apply it to datetime-indexed
    DataFrame.

    Arguments:
        dt_df (:obj:`pandas.DataFrame`): datetime-indexed DataFrame
        years (str or int): years to select, e.g. 2015 or 2000-2010
        label (str): identifier to print warning message if ``years``
            filter partially overlaps with actual date index

    Returns:

        ret (tuple of (:obj:`pandas.DataFrame`, str)): first element is
            input DataFrame ``dt_df`` indexed to ``years`` filter,
            second element is string of year range, e.g. '2001_2011'

    Example:

        >>> df = pd.DataFrame(index=pd.date_range('2000', '2015'))
        >>> df, yr_str = parse_yr_filter(df, '1998-2002', 'station1')
        WARNING: data for station1 starts in 2000 but you gave 1998
        Years used will only include 2000 to 2002

        Now df will only contain indices with dates between 2000 and
        2002 and

        >>> yr_str
        '1998_2002'

    Raises:
        ValueError: if ``years`` is invalid or not found
            in time series index of DataFrame.

    """
    err_msg = ('{} is not a valid years option,\n'.format(years), 'use single or range e.g. 2015 or 2000-2010')
    if years == 'all':
        year_str = 'all_yrs'
    else:
        try:
            if years and isinstance(years, str) and '-' in years:
                start, end = years.strip().split('-')
                year_str = '{}_{}'.format(start, end)
                data_start = start
                data_end = end
                # the assignment on the next line will not raise an
                # exception even if the full date range is missing
                dt_df = dt_df.loc[start:end]
                if start not in dt_df.index:
                    data_start = dt_df.index.year.min()
                    print('WARNING: data for {l} starts in {d}'.format(l=label, d=data_start) +
                          ' but you gave {s}'.format(s=start))
                if end not in dt_df.index:
                    data_end = dt_df.index.year.max()
                    print('WARNING: data for {l} ends in {d}'\
                              .format(l=label, d=data_end) +\
                         ' but you gave {e}'.format(e=end))
                if data_start != start or data_end != end:
                    print('Years used will only include {} to {}'\
                              .format(data_start, data_end))
            else:
                year_str = str(int(years))
                if not len(year_str) == 4:
                    raise ValueError(err_msg)
                if not years in dt_df.index:
                    print('WARNING:', label, 'is missing data',
                        'for year:', years)
                    data_start = dt_df.index.year.min()
                    data_end = dt_df.index.year.max()
                    print('Years used will only include {} to {}'\
                              .format(data_start, data_end))
                else:
                    dt_df = dt_df.loc[years]
        except:
            raise ValueError(err_msg)

    ret = dt_df, year_str
    return ret


def validate_file(file_path, expected_extensions):
    """
    Checks to see if provided path is valid, while also checking to see if file is of expected type.
    Raises exceptions if either of those fail.

    Args:
        file_path: string of path to file
        expected_extensions: list of strings of expected file types

    Returns:
        None
    """
    # Check to see if provided config file path actually points to a file.
    if pl.Path(file_path).is_file():

        # Next check to see if provided file is of the appropriate type.
        file_extension = pl.PurePath(file_path).suffix
        file_extension = file_extension.split('.', 1)[1]  # Remove period
        file_extension = file_extension.lower()  # Make it lowercase

        if file_extension not in expected_extensions:
            raise IOError('\n\nProvided file was of type \'{}\' but script was expecting type \'{}\'.'
                          .format(file_extension, expected_extensions))
        else:
            pass
    else:
        raise IOError('\n\nUnable to find the file at path \'{}\'.'.format(file_path))


def read_config(config_file_path):
    """
    Opens config file at provided path and stores all required values in a python dictionary. This dictionary will be
    used both to import data and elsewhere in the code to refer to what type of data was passed in

    Args:
        config_file_path: string of path to config file

    Returns:
        config_dict: a dictionary of all required config file parameters

    """

    # Check to see if provided file exists and also that it is the correct type
    validate_file(config_file_path, 'ini')

    # Open ConfigParser and point it to file.
    config_reader = cp.ConfigParser()
    config_reader.read(config_file_path)

    # Create config file dictionary and start adding entries to it
    # The DATA and UNITS sections are all strings, so just import the config_reader dictionaries
    config_dict = {**config_reader._sections['DATA'], **config_reader._sections['UNITS']}

    # METADATA Section
    # Projection information
    config_dict['input_data_projection'] =\
        config_reader['METADATA']['input_data_projection'] 
    config_dict['grid_resolution'] =\
        config_reader.getfloat('METADATA','grid_resolution', fallback=0.1)
    config_dict['interpolation_projection'] =\
        config_reader.get('METADATA','interpolation_projection', fallback='ESRI:102004')
    config_dict['interpolation_resolution'] =\
        config_reader.getfloat('METADATA','interpolation_resolution', fallback=1000)
    config_dict['output_data_projection'] =\
        config_reader.get('METADATA','output_data_projection', fallback='ESRI:4326')
    config_dict['output_data_resolution'] =\
        config_reader.getfloat('METADATA','output_data_resolution', fallback=0.1)

    # Below variables are for obtaining decimal places on resolution if it's a float
    # might be useful in developing eventual way to force snapping to grid
    for res in ['grid_resolution', 'interpolation_resolution', 'output_data_resolution']:
        if '.' in str(config_dict[res]):
            config_dict[f'{res}_decimals'] = \
                len(str(config_dict[res]).split('.')[1])
        else:
            config_dict[f'{res}_decimals'] = 0

    # Bounding information
    config_dict['input_bounds'] = {}
    config_dict['input_bounds']['xmin'] = config_reader['METADATA'].getfloat('xmin')
    config_dict['input_bounds']['xmax'] = config_reader['METADATA'].getfloat('xmax')
    config_dict['input_bounds']['ymin'] = config_reader['METADATA'].getfloat('ymin')
    config_dict['input_bounds']['ymax'] = config_reader['METADATA'].getfloat('ymax')

    # Gridded dataset information
    config_dict['collection_info'] = {}
    config_dict['collection_info']['name'] = config_reader['METADATA']['collection_name']
    config_dict['collection_info']['path'] = config_reader['METADATA']['collection_path']
    config_dict['collection_info']['start_date'] = config_reader['METADATA']['start_date']
    config_dict['collection_info']['end_date'] = config_reader['METADATA']['end_date']

    # File structure information
    config_dict['station_anemometer_height'] = config_reader['METADATA'].getfloat('station_anemometer_height')
    config_dict['station_lines_of_header'] = config_reader['METADATA'].getint('station_lines_of_header')
    config_dict['station_missing_data_value'] = config_reader['METADATA']['station_missing_data_value']
    config_dict['gridded_anemometer_height'] = config_reader['METADATA'].getfloat('gridded_anemometer_height')
    config_dict['gridded_lines_of_header'] = config_reader['METADATA'].getint('gridded_lines_of_header')
    config_dict['gridded_missing_data_value'] = config_reader['METADATA']['gridded_missing_data_value']

    # Check to see that all expected variables are provided, for now just print out a warning letting the user know
    # but also change all empty strings into None
    if '' in config_dict.values():
        missing_keys = [key for (key, value) in config_dict.items() if value == '']

        for key in missing_keys:
            config_dict[key] = None

        print('\n\nThe following parameters were unspecified in the config file: {}.'.format(missing_keys))
    else:
        pass

    return config_dict


def read_data(config_dictionary, version, filepath):
    """
    Uses config_dict parameters to read in the data and rename it to standard parameters

    Args:
        config_dictionary: dictionary of everything
        version: a string that will be either 'station' or 'gridded'

    Returns:
        filtered_df: a dataframe containing only the variable we want to plot, with a standardized naming convention

    """

    # Generate vars corresponding to config_dict keys
    version = version + '_'
    loh = version + 'lines_of_header'
    missing_val = version + 'missing_data_value'
    date_column = version + 'date_col'

    # Open file, or station, data file
    (_file_name, file_extension) = os.path.splitext(filepath)
    validate_file(filepath, ['csv', 'xls', 'xlsx'])
    if file_extension == '.csv':  # csv file provided
        raw_file_data = pd.read_csv(
            filepath, delimiter=',', header=config_dictionary[loh]-1,
            index_col=config_dictionary[date_column], parse_dates=True, engine='python',
            na_values=config_dictionary[missing_val], keep_default_na=True,
            na_filter=True, skip_blank_lines=True)

    elif file_extension in ['.xls', '.xlsx']:
        raw_file_data = pd.read_excel(
            filepath, sheet_name=0, header=config_dictionary[loh]-1,
            index_col=config_dictionary[date_column], parse_dates=True, engine='openpyxl',
            na_values=config_dictionary[missing_val], keep_default_na=True,
            na_filter=True)

    else:
        # This script is only handles csv and excel files. Validate_file() already catches this case
        raise IOError('\n\nProvided file was of type \'{}\' but script was expecting type \'{}\'.'
                      .format(file_extension, ['csv', 'xls', 'xlsx']))

    # Create handling for 'unnamed:0' and 'datetime' column in station data files
    if 'Unnamed: 0' in raw_file_data.columns:
        raw_file_data.rename(columns={'Unnamed: 0': 'date'}, inplace=True)
        raw_file_data.set_index('date', drop=True, inplace=True)
    elif 'datetime' in raw_file_data.columns:
        raw_file_data.rename(columns={'datetime': 'date'}, inplace=True)
        raw_file_data.set_index('date', drop=True, inplace=True)

    # iterate through an expected list of vars and append a column should one be missing, to prevent a key error later
    var_list = ['tmax', 'tmin', 'tdew', 'rs', 'wind', 'rhmax', 'rhmin', 'rhavg', 'ea', 'eto', 'etr', 'prcp']
    for var in var_list:
        var_col = version + var + '_col'

        if config_dictionary[var_col] is None:  # var wasn't provided, create empty column
            empty_col = np.empty(len(raw_file_data))
            empty_col[:] = np.nan
            raw_file_data[var] = empty_col

        elif config_dictionary[var_col] is not None and config_dictionary[var_col] not in list(raw_file_data.columns):
            # var is provided but doesn't match any column in the data file
            raise ValueError(
                '\n\n\'{}\' was specified in the config file as \'{}\' but that '
                'column was not found in the data file \'{}\'.'
                .format(var_col, config_dictionary[var_col], filepath))

        else:  # var was provided, so just rename it to the standard naming convention
            raw_file_data.rename(columns={config_dictionary[var_col]: var}, inplace=True)

    filtered_df = pd.DataFrame(data=raw_file_data[var_list])

    return filtered_df


def convert_units(config_dictionary, version, df):
    """
        Uses config_dict parameters to check what units provided variables are in and convert them if needed

        Args:
            config_dictionary: dictionary of everything contained within config file
            version: a string that will be either 'station' or 'gridded'
            df: pandas dataframe of input data, at this point naming of dataframe columns has been standardized

        Returns:
            converted_df: a dataframe containing data in the correct units

        """
    version = version + '_'

    converted_df = df.copy(deep=True)
    # iterate through list of vars to convert each
    # todo make these lists into a dict, and allow for column order parameters in the config file instead of names
    var_list = ['tmax', 'tmin', 'tdew', 'rs', 'wind', 'ea', 'rhmax', 'rhmin', 'rhavg', 'eto', 'etr', 'prcp']
    units_list = ['temp', 'temp', 'temp', 'solar', 'wind', 'ea', 'rh', 'rh', 'rh', 'et', 'et', 'prcp']
    for i in range(len(var_list)):
        var_col = version + var_list[i] + '_col'
        var_units_key = version + units_list[i] + '_units'
        var_units = str(config_dictionary[var_units_key]).lower()

        if config_dictionary[var_col] is None:
            # var is not provided, so just pass through empty column
            converted_data = np.array(df.shape[0] * np.nan)
        elif config_dictionary[var_col] is not None and config_dictionary[var_units_key] is None:
            # var is provided but units aren't specified, raise an error
            raise ValueError('\n\n\'{}\' was specified in the config file but the parameter \'{}\' was unspecified.'
                             .format(var_col, var_units_key))
        else:
            # everything is provided, convert units if necessary

            if units_list[i] == 'temp':
                if var_units == 'c':
                    converted_data = np.array(df[var_list[i]])
                elif var_units == 'f':
                    converted_data = np.array(((df[var_list[i]] - 32.0) * (5.0 / 9.0)))
                elif var_units == 'k':
                    converted_data = np.array(df[var_list[i]] - 273.15)
                else:
                    raise ValueError(
                        '\n\n\'{}\' was specified in the config file as having units \'{}\' which is not a valid option.'
                        .format(var_units_key, config_dictionary[var_units_key]))

            elif units_list[i] == 'solar':
                if var_units == 'w/m2':
                    converted_data = np.array(df[var_list[i]])
                elif var_units == 'j/m2':
                    converted_data = np.array((df[var_list[i]] / 1000000) * 11.574)  # j/m2 to w/m2
                elif var_units == 'mj/m2':
                    converted_data = np.array(df[var_list[i]] * 11.574)  # mj/m2 to w/m2
                elif var_units == 'langleys' or var_units == 'lang':
                    converted_data = np.array((df[var_list[i]] * 0.484583)) # langleys to w/m2
                elif var_units == 'kw-hr/m2':
                    converted_data = np.array((df[var_list[i]] * 1000) / 24)  # kw-hr/m2 to w/m2
                else:
                    raise ValueError(
                        '\n\n\'{}\' was specified in the config file as having units \'{}\' which is not a valid option.'
                        .format(var_units_key, config_dictionary[var_units_key]))

            elif units_list[i] == 'wind':
                if var_units == 'm/s':
                    converted_data = np.array(df[var_list[i]])
                elif var_units == 'mph':
                    converted_data = np.array(df[var_list[i]] * 0.44704)  # mph to m/s
                elif var_units == 'kmhr':
                    converted_data = np.array(df[var_list[i]] / 3.6)  # Km/hr to m/s
                else:
                    raise ValueError(
                        '\n\n\'{}\' was specified in the config file as having units \'{}\' which is not a valid option.'
                        .format(var_units_key, config_dictionary[var_units_key]))

            elif units_list[i] == 'ea':
                if var_units == 'kpa':
                    converted_data = np.array(df[var_list[i]])
                elif var_units == 'hpa':
                    converted_data = np.array(df[var_list[i]] * 0.1)  # hPa to kPa
                elif var_units == 'torr':
                    converted_data = np.array(df[var_list[i]] * 0.133322)  # Torr to kPa
                elif var_units == 'mbar':
                    converted_data = np.array(df[var_list[i]] * 0.1)  # Mbar to kPa
                else:
                    raise ValueError(
                        '\n\n\'{}\' was specified in the config file as having units \'{}\' which is not a valid option.'
                        .format(var_units_key, config_dictionary[var_units_key]))

            elif units_list[i] == 'rh':
                if var_units == 'percent':
                    converted_data = np.array(df[var_list[i]])
                elif var_units == 'fraction':
                    converted_data = np.array(df[var_list[i]] * 100.0)  # fraction to %
                else:
                    raise ValueError(
                        '\n\n\'{}\' was specified in the config file as having units \'{}\' which is not a valid option.'
                        .format(var_units_key, config_dictionary[var_units_key]))

            elif units_list[i] == 'et' or units_list[i] == 'prcp':
                if var_units == 'mm':
                    converted_data = np.array(df[var_list[i]])
                elif var_units == 'inches' or var_units == 'in':
                    converted_data = np.array(df[var_list[i]] * 25.4)  # inches to mm
                else:
                    raise ValueError(
                        '\n\n\'{}\' was specified in the config file as having units \'{}\' which is not a valid option.'
                        .format(var_units_key, config_dictionary[var_units_key]))
            else:
                raise ValueError(units_list[i] + ' is not a valid var unit code')

        # add converted var into dataframe of converted data
        converted_df[var_list[i]] = converted_data

    return converted_df


def reproject_crs_for_point(orig_lon, orig_lat, orig_crs, requested_crs):
    """
        Uses the pyproj library to reproject point data from one CRS to another
            ex. will be used to make input coords wgs84 for earth engine

        Will return original data without any reprojection if orig_crs
            and requested_crs are the same

        Args:
            orig_lon: float of original longitude
            orig_lat: float of original latitude
            orig_crs: string of EPSG code for orig_lat and orig_lon
            requested_crs: string of EPSG code to reproject into
        Returns:
            Reprojected latitude and longitude for point
    """
    if orig_crs == requested_crs:
        return orig_lon, orig_lat

    proj_transformer = pyproj.Transformer.from_crs(orig_crs, requested_crs,
                                                   always_xy=True)
    return proj_transformer.transform(orig_lon, orig_lat)


def reproject_crs_for_bounds(bounds, resolution, orig_crs, requested_crs,
                             requested_decimals):
    """
        Uses the pyproj library to reproject dictionary of bounds for
            interpolation extent. This is done in more than just two calls
            (ex. NW and SE corners) as some projections may have curvature

        Afterwords it rounds the coordinates to the requested decimals

        If orig_crs and requested_crs are the same it will just round the coords
            without reprojecting

        Args:
            bounds: dictionary of bounds, containing the following keys:
                xmin, xmax, ymin, ymax
            resolution: resolution used for interpolation, coordinates will
                be rounded in an attempt to snap to grid
            orig_crs: string of EPSG code for original bounds
            requested_crs: string of EPSG code to reprojected bounds
            requested_decimals: int of number of decimals to round coords to
        Returns:
            Reprojected bounds into new CRS
    """

    if orig_crs == requested_crs:
        projection_dict = {key: value for key, value in bounds.items()}
    else:
        projection_dict = {}
        proj_transformer = (
            pyproj.Transformer.from_crs(orig_crs, requested_crs,
                                        always_xy=True))

        # Calculate xmin at the SW corner
        projection_dict['xmin'], _ignore =\
            proj_transformer.transform(bounds['xmin'], bounds['ymin'])
        # Calculate xmax at the SE corner
        projection_dict['xmax'], _ignore =\
            proj_transformer.transform(bounds['xmax'], bounds['ymin'])
        # Calculate ymax at the NE corner, could've also been NW corner
        _ignore, projection_dict['ymax'] =\
            proj_transformer.transform(bounds['xmax'], bounds['ymax'])
        # Calculate ymin as at the average between the east and west extent
        _ignore, projection_dict['ymin'] =\
            proj_transformer.transform(
                (bounds['xmax'] + bounds['xmin']) / 2, bounds['ymin']
            )

    # Round the entries in the reproj dict to the resolution if above 1
    # Mainly used to cut off decimal places on projections defined in meters
    for key in projection_dict.keys():
        if requested_decimals > 0:
            projection_dict[key] = round(projection_dict[key], requested_decimals)
        else:
            # if resolution is above 1, turn it into an int and subtract modulo
            int_res = int(projection_dict[key])
            remainder = int_res % int(resolution)
            projection_dict[key] = int_res - remainder

    return projection_dict
