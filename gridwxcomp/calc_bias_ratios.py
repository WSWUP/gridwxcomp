# -*- coding: utf-8 -*-
"""
Calculate monthly bias ratios of variables from climate station 
to overlapping gridded dataset cells.

Input file for this module must first be created by running 
:mod:`gridwxcomp.prep_metadata` followed by :mod:`gridwxcomp.ee_download`.

Note: The module is reliant on values within the specified config file
    in order to interpret the station and gridded data values successfully.
    If you are experiencing errors or bad values please check it is set up
    correctly.
    
"""


import os
import calendar
import numpy as np
import pandas as pd
import warnings
from refet.calcs import _wind_height_adjust
from .util import parse_yr_filter, read_config, read_data, convert_units


VAR_LIST = ['tmax', 'tmin', 'tdew', 'rs', 'wind', 'ea', 'rhmax', 'rhmin', 'rhavg', 'eto', 'etr', 'prcp']
GROW_THRESH = 65
SUM_THRESH = 35
ANN_THRESH = 125


def _save_output(out_df, comp_out_df, out_dir, grid_id, var_name, yrs):
    """
    Save short summary file or overwrite existing data for a single
    climate station.
    
    Arguments:
        out_df (:class:`pandas.DataFrame`): data containing short
            summary info, mainly mean monthly bias ratios for a 
            single climate station to save.
        comp_out_df (:class:`pandas.DataFrame`, bool): either a 
            single row dataframe with comprehensive summary data
            or False (default). depends on ``comp`` argument to 
            :func:`calc_bias_ratios`. If :class:`pandas.DataFrame` 
            is passed then save or update existing file.
        out_dir (str): path to directory to save or update summary data
            for monthly bias ratios.
        grid_id (int, optional): depends on ``grid_id`` argument
            passed to :func:`calc_bias_ratios`. If not None (default)
            then save summary files for stations that correspond with
            the given gridMET ID with the suffix "_X" where X is the
            gridMET ID value.
        var_name (str): name of gridMET variable that is being processed.
        yrs (str): years used to calc ratios, save to out files as suffix.
    
    Returns:
        None       
        
    """

    def __save_update(out_df, out_file):

        """ 
        Helper function that is reused to save or update both short and
        long summary files one station or row at a time. Saves station ratio
        data by appending to existing file or overwriting data for a station if
        it was previously calculated. `out_df` is a single row from the ratio
        results table representing data for a single climate station-gridcell
        pair. 
        """

        # if short file exists add/overwrite row for station
        if os.path.isfile(out_file):
            existing_df = pd.read_csv(out_file, index_col='STATION_ID')
            if not out_df.index.values[0] in existing_df.index.values:
                out_df = pd.concat([existing_df, out_df], sort=False)
                out_df.to_csv(out_file, na_rep=-999, index=True)
            # overwrite if station is in existing, could change to
            # allow for duplicates if values are different
            else:
                existing_df.loc[out_df.index.values[0], :] =\
                    out_df.loc[out_df.index.values[0]]
                existing_df.to_csv(out_file, na_rep=-999, index=True)
        else:
            out_df.to_csv(out_file, na_rep=-999, index=True)
            
    # save or update short and comprehensive summary files
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)   
        
    # save/update short summary file or update existing with new station
    if not grid_id:
        out_file = os.path.join(
            out_dir, 
            '{v}_summary_{y}.csv'.format(v=var_name, y=yrs)
        )
    else: 
        out_file = os.path.join(
            out_dir,
            '{v}_summary_grid_{g}_{y}.csv'.format(
                v=var_name, g=grid_id, y=yrs)
        )
    __save_update(out_df, out_file)

    # if comprehensive summary is requested save/update
    if isinstance(comp_out_df, pd.DataFrame):
        if not grid_id:
            comp_out_file = os.path.join(
                out_dir, 
                '{v}_summary_comp_{y}.csv'.format(v=var_name, y=yrs)
            )
        else:
            comp_out_file = os.path.join(
                out_dir,
                '{v}_summary_comp_{g}_{y}.csv'.format(
                    v=var_name, g=grid_id, y=yrs)
            )
        __save_update(comp_out_df, comp_out_file)


def calc_bias_ratios(input_path, config_path, out_dir, method='long_term_mean', grid_id_name='GRID_ID',
                     comparison_var='etr', grid_id=None, day_limit=10, years='all', comp=True):
    """
    Read input metadata CSV file and config file, use them to calculate mean
    monthly bias ratios between station to corresponding grid cells for all
    station and grid pairs, optionally calculate ratios for a single gridcell.
    
    Arguments:
        input_path (str): path to input CSV file with matching station climate and grid metadata.
            This file is created by running :func:`gridwxcomp.prep_metadata` followed by
            :func:`gridwxcomp.ee_download`.
        config_path (str): path to the configuration file that has the 
            parameters used to interpret the station and gridded data files.
        out_dir (str): path to directory to save CSV files with monthly bias ratios of etr.
        method (str): default 'long_term_mean'. How to calculate mean station to grid ratios,
            currently two options 'long_term_mean' takes the mean of all dates for the station variable
            that fall in a time periods, e.g. the month of January, to the mean of all paired
            January dates in the gridded product. The other option is  'mean_of_annual' which calculates ratios,
            for each time period if enough paired days exist, the ratio of sums for each year in the
            record and then takes the mean of the annual ratios. This method is always used to calculate
            standard deviation and coefficient of variation of ratios which describe interannual variation of ratios.
        grid_id_name (str): default 'GRID_ID'. Name of column containing index/cell identifiers for gridded dataset.
        comparison_var (str): default 'etr'. Grid climate variable to calculate bias ratios. This value must be found
            within the module parameter VAR_LIST.
        grid_id (int or str or None): default None. Grid ID (int cell identifier) to only
            calculate bias ratios for a single gridcell.
        day_limit (int): default 10. Threshold number of days in month of missing data,
            if less, exclude month from calculations. Ignored when ``method='long_term_mean'``.
        years (int or str): default 'all'. Years to use for calculations e.g. 2000-2005 or 2011.
        comp (bool): default True. Flag to save a "comprehensive" 
            summary output CSV file that contains station metadata and 
            statistics in addition to the mean monthly ratios.
                    
    Returns:
        None
        
    Example:
        To use within Python for observed ETo
        
        >>> from gridwxcomp import calc_bias_ratios
        >>> input_file = 'prepped_metadata.csv'
        >>> config_file = 'gridwxcomp_config.ini'
        >>> out_directory = 'monthly_ratios'
        >>> grid_variable = 'eto'
        >>> calc_bias_ratios(input_file, config_file, out_directory, comparison_var=grid_variable, comp=True)

        This results in two CSV files in ``out_directory`` named
        "eto_summary_all_yrs.csv" and "eto_summary_comp_all_yrs.csv".

    Raises:
        FileNotFoundError: if input file or config file are invalid or not found.
        KeyError: if the input file does not contain file paths to the climate station and 
            grid time series files. This occurs if, for example, the 
            :mod:`gridwxcomp.prep_metadata` and/or
            :mod:`gridwxcomp.ee_download` scripts have not been run first.
            Also raised if the given the values specified in the config file are not found
            within the station and gridded data files.
        ValueError: if the ``method`` kwarg is invalid.

    Note:
        Growing season and summer periods over which ratios are calculated are
        defined as April through October and June through August respectively. 
    
    Note:
        If an existing summary file contains a climate station that is being 
        reprocessed its monthly bias ratios and other data will be overwritten. 
        Also, to proceed with spatial analysis scripts, the comprehensive 
        summary file must be produced using this function first.
    """

    # ignore np runtime warnings due to calcs with nans, div by 0
    np.seterr(divide='ignore', invalid='ignore')
    # specific for standard deviation of nans
    std_warning = "Degrees of freedom <= 0 for slice"
    warnings.filterwarnings("ignore", message=std_warning)

    # Define months belonging to growing season and annual
    grow_months = list(range(4, 11))
    ann_months = list(range(1, 13))

    method_options = ('long_term_mean', 'mean_of_annual')
    if method not in method_options:
        raise ValueError('{} is not a valid method, use one of: {}'.format(
            method, method_options))

    if not os.path.isdir(out_dir):
        print('{} does not exist, creating directory'.format(out_dir))
        os.mkdir(out_dir)
    if not os.path.isfile(input_path):
        raise FileNotFoundError('Input CSV file given was invalid or not found')
    # If only calculating ratios for a single cell, change console message
    if grid_id:
        single_grid_cell_msg = f'For grid cell ID: {grid_id}.'
    else:
        single_grid_cell_msg = ''

    if comparison_var not in VAR_LIST:
        raise ValueError('{} is not a valid option, use one of: {}'.format(
            comparison_var, VAR_LIST))

    station_var = f'station_{comparison_var}'
    grid_var = f'gridded_{comparison_var}'

    print(
        f'Calculating ratios between climate station variable: {station_var}'
        f'\nand grid variable: {grid_var} using the "{method.replace("_"," ")}"'
        f' method. {single_grid_cell_msg}'
    )

    # Read in metadata file to iterate over and config file to interpret data files
    input_df = pd.read_csv(input_path)
    config_dict = read_config(config_path)

    # loop through each station and calculate monthly ratio
    for index, row in input_df.iterrows():
        if 'STATION_FILE_PATH' not in row or 'GRID_FILE_PATH' not in row:
            raise KeyError('Missing station and/or grid file paths in ' +
                'input file. Run prep_metadata.py followed by ee_download.py first.')

        # if only doing a single grid cell check for matching ID
        if grid_id and grid_id != row[grid_id_name]:
            continue

        # load station and grid time series files
        try:
            # open file, convert units, then adjust wind measurement height
            raw_station_df = read_data(config_dict, 'station', row.STATION_FILE_PATH)
            converted_station_df = convert_units(config_dict, 'station', raw_station_df)
            converted_station_df['wind'] = _wind_height_adjust(
                uz=converted_station_df['wind'], zw=config_dict['station_anemometer_height'])
            prepped_station_df = converted_station_df.add_prefix('station_', axis='columns')
        except IOError:
            print(
                'Time series file for station: ', row.STATION_ID, ' was not found, skipping.')
            continue

        try:
            # open file, convert units, then adjust wind measurement height
            raw_gridded_df = read_data(config_dict, 'gridded', row.GRID_FILE_PATH)
            converted_gridded_df = convert_units(config_dict, 'gridded', raw_gridded_df)
            converted_gridded_df['wind'] = _wind_height_adjust(
                uz=converted_gridded_df['wind'], zw=config_dict['gridded_anemometer_height'])
            prepped_gridded_df = converted_gridded_df.add_prefix('gridded_', axis='columns')
        except IOError:
            print(
                'Time series file for gridded: ', row.STATION_ID, 'was not found, skipping.')
            continue

        if station_var not in prepped_station_df.columns:
            err_msg = '{v} not found in the station file: {p}'.format(
                v=station_var, p=row.STATION_FILE_PATH)
            raise KeyError(err_msg)
        print('\nIndex {u}: Calculating {v} bias ratios for station:'.format(
            u=index, v=grid_var), row.STATION_ID)

        # merge both datasets drop missing days
        result = pd.concat(
            [prepped_station_df[station_var], prepped_gridded_df[grid_var]],
            axis=1 
        )
        result = result.reindex(prepped_gridded_df.index)
        result.dropna(inplace=True)
        # make datetime index
        result.index = pd.to_datetime(result.index)
        # apply year filter
        result, years_str = parse_yr_filter(result, years, row.STATION_ID)
        # for calculating ratios with long-term means later
        orig = result.copy()
        # monthly sums and day counts for each year
        result = result.groupby([result.index.year, result.index.month]).agg(['sum', 'mean', 'count'])
        result.index.set_names(['year', 'month'], inplace=True)
        # remove totals with less than XX days
        result = result[result[grid_var, 'count'] >= day_limit]

        # calc mean growing season and June to August ratios with month sums
        # temperature deltas
        if grid_var in ('gridded_tmin', 'gridded_tmax', 'gridded_tdew'):
            grow_season = \
                result.loc[result.index.get_level_values('month').isin(grow_months), station_var]['mean'].mean() - \
                result.loc[result.index.get_level_values('month').isin(grow_months), grid_var]['mean'].mean()
            june_to_aug = \
                result.loc[result.index.get_level_values('month').isin([6, 7, 8]), station_var]['mean'].mean() - \
                result.loc[result.index.get_level_values('month').isin([6, 7, 8]), grid_var]['mean'].mean()
            annual = \
                result.loc[result.index.get_level_values('month').isin(ann_months), station_var]['mean'].mean() - \
                result.loc[result.index.get_level_values('month').isin(ann_months), grid_var]['mean'].mean()
        # ratios for other variables
        else:
            grow_season = \
                result.loc[result.index.get_level_values('month').isin(grow_months), station_var]['mean'].mean() / \
                result.loc[result.index.get_level_values('month').isin(grow_months), grid_var]['mean'].mean()
            june_to_aug = \
                result.loc[result.index.get_level_values('month').isin([6, 7, 8]), station_var]['mean'].mean() / \
                result.loc[result.index.get_level_values('month').isin([6, 7, 8]), grid_var]['mean'].mean()
            annual = \
                result.loc[result.index.get_level_values('month').isin(ann_months), station_var]['mean'].mean() / \
                result.loc[result.index.get_level_values('month').isin(ann_months), grid_var]['mean'].mean()

        ratio = pd.DataFrame(columns=['ratio', 'count'])
        # ratio/deltas of monthly sums for each year
        if grid_var in ('gridded_tmin', 'gridded_tmax', 'gridded_tdew'):
            ratio['ratio'] = (result[station_var, 'mean']) - (result[grid_var, 'mean'])
        else:
            ratio['ratio'] = (result[station_var, 'sum']) / (result[grid_var, 'sum'])

        # monthly counts and stddev
        ratio['count'] = result.loc[:, (grid_var, 'count')]
        if result.empty:
            print(f'WARNING: no data for site: {row.STATION_ID}, skipping')
            continue

        # rebuild Index DateTime
        ratio['year'] = ratio.index.get_level_values('year').values.astype(int)
        ratio['month'] = ratio.index.get_level_values('month').values.astype(int)
        ratio.index = pd.to_datetime(ratio.year*10000+ratio.month*100+15, format='%Y%m%d')

        # useful to know how many years were used in addition to day counts
        start_year = ratio.year.min()
        end_year = ratio.year.max()
        counts = ratio.groupby(ratio.index.month).sum()['count']

        # get standard deviation of each years' monthly mean ratio
        stdev = {month: np.std(ratio.loc[ratio.month.isin([month]), 'ratio'].values) for month in ann_months}
        stdev = pd.Series(stdev, name='stdev')

        # mean of monthly means of all years, can change to median or other meth
        final_ratio = ratio.groupby(ratio.index.month).mean()
        final_ratio.drop(['year', 'month'], axis=1, inplace=True)
        final_ratio['count'] = counts
        final_ratio['stdev'] = stdev
        final_ratio['cv'] = stdev / final_ratio['ratio']
        # calc mean growing season, June through August, ann stdev
        grow_season_std = np.std(ratio.loc[ratio.month.isin(grow_months), 'ratio'].values)
        june_to_aug_std = np.std(ratio.loc[ratio.month.isin([6, 7, 8]), 'ratio'].values)
        annual_std = np.std(ratio.loc[ratio.month.isin(ann_months), 'ratio'].values)
        # get month abbreviations in a column and drop index values
        for m in final_ratio.index:
            final_ratio.loc[m, 'month'] = calendar.month_abbr[m]
        # restructure as a row with station index
        months = final_ratio.month.values
        final_ratio = final_ratio.T
        final_ratio.columns = months
        final_ratio.drop('month', inplace=True)
        # add monthly means and counts into single row dataframe
        ratio_cols = [c + '_mean' for c in final_ratio.columns]
        count_cols = [c + '_count' for c in final_ratio.columns]
        stddev_cols = [c + '_stdev' for c in final_ratio.columns]
        coef_var_cols = [c + '_cv' for c in final_ratio.columns]
        # combine all monthly stats
        out_cols = ratio_cols + count_cols + stddev_cols + coef_var_cols
        final_ratio = pd.concat([
            final_ratio.loc['ratio'], 
            final_ratio.loc['count'],
            final_ratio.loc['stdev'],
            final_ratio.loc['cv']
        ])
        final_ratio.index = out_cols
        # transpose so that each station is one row in final output
        final_ratio = final_ratio.to_frame().T
        # assign non-monthly stats, growing season, annual, june-aug
        final_ratio['grow_mean'] = grow_season
        final_ratio['summer_mean'] = june_to_aug
        final_ratio['annual_mean'] = annual
        # day counts for all years in non-monthly periods
        final_ratio['grow_count'] =\
            counts.loc[counts.index.isin(grow_months)].sum()
        final_ratio['summer_count'] =\
            counts.loc[counts.index.isin([6, 7, 8])].sum()
        final_ratio['annual_count'] =\
            counts.loc[counts.index.isin(ann_months)].sum()
        # assign stdev, coef. var. 
        final_ratio['grow_stdev'] = grow_season_std
        final_ratio['summer_stdev'] = june_to_aug_std
        final_ratio['annual_stdev'] = annual_std
        # coefficient of variation
        final_ratio['grow_cv'] = grow_season_std / grow_season
        final_ratio['summer_cv'] = june_to_aug_std / june_to_aug
        final_ratio['annual_cv'] = annual_std / annual
        # start and end years for interpreting annual CV, stdev...
        final_ratio['start_year'] = start_year
        final_ratio['end_year'] = end_year

        # round numerical data before adding string metadata
        for v in final_ratio:
            if '_mean' or '_stdev' or '_cv' in v:
                final_ratio[v] = final_ratio[v].astype(float).round(3)
            else:
                final_ratio[v] = final_ratio[v].astype(float).round(0)

        # set station ID as index
        final_ratio['STATION_ID'] = row.STATION_ID
        final_ratio.set_index('STATION_ID', inplace=True)

        out = final_ratio.copy()
        out.drop(count_cols+stddev_cols+coef_var_cols, axis=1, inplace=True)

        # save grid ID for merging with input table, merge other metadata
        final_ratio[grid_id_name] = row[grid_id_name]    
        final_ratio = final_ratio.merge(input_df, on=grid_id_name)

        # if more than one site in same gridcell- will have multiple rows 
        # after merge, select the one for the current station 
        if final_ratio.shape[0] > 1:
            final_ratio = final_ratio.loc[final_ratio.STATION_ID == row.STATION_ID]
            final_ratio.reset_index(inplace=True)  # for slicing with .at[0]
 
        # long term mean station to mean grid ratio calc as opposed to mean of
        # annual ratios - default less bias potential
        if method == 'long_term_mean':
            month_means = orig.groupby(orig.index.month).mean()
            # cast as str to avoid deprecation warning assigning str to int
            month_means['month'] = month_means.index.astype(str) 
            for m in month_means.index:
                month_means.loc[m, 'month'] = f'{calendar.month_abbr[int(m)]}_mean'
            month_means.set_index('month', inplace=True)
            if grid_var in ('gridded_tmin', 'gridded_tmax', 'gridded_tdew'):
                month_means['ratios'] = month_means[station_var] - month_means[grid_var]
            else:
                month_means['ratios'] = month_means[station_var] / month_means[grid_var]

            long_term = month_means.drop(labels=[station_var, grid_var], axis=1).T

            # non-monthly periods long-term mean to mean ratios
            grow_season = orig.loc[orig.index.month.isin(grow_months)]
            summer_season = orig.loc[orig.index.month.isin([6, 7, 8])]
            if grid_var in ('gridded_tmin', 'gridded_tmax', 'gridded_tdew'):
                long_term['grow_mean'] = grow_season[station_var].mean() - grow_season[grid_var].mean()
                long_term['summer_mean'] = summer_season[station_var].mean() - summer_season[grid_var].mean()
                long_term['annual_mean'] = orig[station_var].mean() - orig[grid_var].mean()

            else:
                long_term['grow_mean'] = grow_season[station_var].mean() / grow_season[grid_var].mean()
                long_term['summer_mean'] = summer_season[station_var].mean() / summer_season[grid_var].mean()
                long_term['annual_mean'] = orig[station_var].mean() / orig[grid_var].mean()

            # overwrite only mean ratios (keep stats from mean of annual ratios)
            overwrite = long_term.columns.intersection(final_ratio.columns)
            # return long_term, overwrite, final_ratio
            final_ratio[overwrite] = long_term[overwrite].values
            out[overwrite] = long_term[overwrite].values

        final_ratio['ratio_method'] = method

        # round numeric columns
        final_ratio = final_ratio.round({
            'LAT': 10,
            'LON': 10,
            'ELEV_M': 0,
            'ELEV_FT': 0,
            'STATION_LAT': 10,
            'STATION_LON': 10,
            'STATION_ELEV_M': 0
        })

        # check if day counts for non-monthly periods are too low, if assign na
        if final_ratio.at[0, 'summer_count'] < SUM_THRESH:
            print('WARNING: less than:', SUM_THRESH, 'days in summer period', '\nfor station:',
                  row.STATION_ID, 'assigning -999 for all stats')
            cols = [col for col in final_ratio.columns if 'summer_' in col and '_count' not in col]
            final_ratio.loc[:, cols] = np.nan
        
        if final_ratio.at[0, 'grow_count'] < GROW_THRESH:
            print('WARNING: less than:', GROW_THRESH, 'days in growing season', '\nfor station:',
                  row.STATION_ID, 'assigning -999 for all stats')
            cols = [col for col in final_ratio.columns if 'grow_' in col and '_count' not in col]
            final_ratio.loc[:, cols] = np.nan
            
        if final_ratio.at[0, 'annual_count'] < ANN_THRESH:
            print('WARNING: less than:', ANN_THRESH, 'days in annual period', '\nfor station:',
                  row.STATION_ID, 'assigning -999 for all stats')
            cols = [col for col in final_ratio.columns if 
                    'annual_' in col and '_count' not in col]
            final_ratio.loc[:, cols] = np.nan

        if comp:
            out[grid_id_name] = row[grid_id_name]
            out[grid_id_name] = final_ratio[grid_id_name].unique()
            # build comprehensive output summary 
            comp_out = final_ratio
            comp_out.set_index('STATION_ID', inplace=True)

            # no longer need grid ID in short summary 
            out.drop(columns=grid_id_name, inplace=True)
        # if comp False
        else:
            comp_out = comp

        # save output depending on options
        _save_output(
            out, comp_out, out_dir, grid_id, 
            grid_var.replace('gridded_',''), years_str)

    print('\nSummary file(s) for bias ratios saved to: \n', os.path.abspath(out_dir))    


if __name__ == '__main__':
    print('\n--------------------------------------------------------'
          ' Functionality for running this library from the terminal'
          ' was removed. Please refer to the documentation on how to'
          ' make calls to these functions. \n\n')
