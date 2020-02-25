# -*- coding: utf-8 -*-
"""
Calculate monthly bias ratios of variables from climate station 
to overlapping gridMET (or other gridded dataset) cells. 

Input file for this module must first be created by running 
:mod:`gridwxcomp.prep_input` followed by :mod:`gridwxcomp.download_gridmet_opendap`. 

Attributes:
    GRIDMET_STATION_VARS (:obj:`dict`): mapping dictionary with gridMET
        variable names as keys and station variable names as values.
        Used to determine which station variable to calculate bias
        ratios according to the given gridMET variable. 

Default values::

    GRIDMET_STATION_VARS = {
        'u2_ms' : 'ws_2m (m/s)',
        'tmin_c' : 'TMin (C)',
        'tmax_c' : 'TMax (C)',
        'srad_wm2' : 'Rs (w/m2)',
        'ea_kpa' : 'Vapor Pres (kPa)',
        'prcp_mm' : 'Precip (mm)',
        'etr_mm' : 'ETr (mm)',
        'eto_mm' : 'ETo (mm)'
    }

Note: The module attribute ``GRIDMET_STATION_VARS`` can be manually adjusted,
    if ``gridwxcomp`` is installed in editable mode or used as scripts from the
    root directory. New pairs of station-to-grid variables can then be made
    or removed to efficiently use :mod:`gridwxcomp` on station data that was
    **not** created by `PyWeatherQAQC
    <https://github.com/WSWUP/pyWeatherQAQC>`_.  Otherwise, the same can be
    achieved by the ``var_dict`` or ``grid_var`` and ``station_var`` arguments 
    to :func:`calc_bias_ratios`.
    
"""

import os
import calendar
import argparse
import warnings

import pandas as pd
import numpy as np

# allows for CL script usage if gridwxcomp not installed
try:
    from .util import parse_yr_filter
except:
    from util import parse_yr_filter

# keys = gridMET variable name
# values = climate station variable name
GRIDMET_STATION_VARS = {
    'u2_ms' : 'ws_2m (m/s)',
    'tmin_c' : 'TMin (C)',
    'tmax_c' : 'TMax (C)',
    'srad_wm2' : 'Rs (w/m2)',
    'ea_kpa' : 'Vapor Pres (kPa)',
    'prcp_mm' : 'Precip (mm)',
    'etr_mm' : 'ETr (mm)',
    'eto_mm' : 'ETo (mm)'
}

OPJ = os.path.join

def main(input_file_path, out_dir, method='long_term_mean', 
        grid_id_name='GRIDMET_ID', grid_var='etr_mm', station_var=None, 
        station_date_name='date', grid_date_name='date', grid_ID=None, 
        day_limit=10, years='all', comp=True):
    """
    Calculate monthly bias ratios between station climate and gridMET
    cells that correspond with each other geographically. Saves data
    to CSV files in the given output directory. If run later with
    new station data, bias ratios for new stations will be appended
    to existing output summary CSV.
    
    Arguments:
        input_file_path (str): path to input CSV file containing
            paired station/gridMET metadata. This file is 
            created by running :mod:`gridwxcomp.prep_input` followed by 
            :mod:`gridwxcomp.download_gridmet_opendap`.
        out_dir (str): path to directory to save CSV files with
            monthly bias ratios of etr.
            
    Keyword Arguments:
        method (str): default 'long_term_mean'. How to calculate mean station to
            grid ratios, currently two options 'long_term_mean' takes the 
            mean of all dates for the station variable that fall in a time 
            periods, e.g. the month of January, to the mean of all paired 
            January dates in the gridded product. The other option is 
            'mean_of_annual' which calculates ratios, for each time period if
            enough paired days exist, the ratio of sums for each year in the 
            record and then takes the mean of the annual ratios. This method 
            is always used to calculate standard deviation and coefficient of 
            variation of ratios which describe interannual variation of ratios. 
        grid_var (str): default 'etr_mm'. Grid climate variable
            to calculate bias ratios.
        station_var (str): default None. Climate station variable to use
            to calculate bias ratios. If None, look up using ``grid_var`` 
            as a key to :attr:`GRIDMET_STATION_VARS` dictionary found as a 
            module attribute to :mod:`gridwxcomp.calc_bias_ratios`.
        grid_ID (int): default None. Grid ID (int cell identifier) to only 
            calculate bias ratios for a single gridcell.
        day_limit (int): default 10. Threshold number of days in month
            of missing data, if less exclude month from calculations.
        years (int or str): default 'all'. Years to use for calculations
            e.g. 2000-2005 or 2011.
        comp (bool): default True. Flag to save a "comprehensive" 
            summary output CSV file that contains station metadata and 
            statistics in addition to the mean monthly ratios.

    Returns:
        None

    Examples:
        From the command line interface,

        .. code-block:: sh

            $ # for all gridMET cells in input file for gridMET var "etr_mm" (default)
            $ python calc_bias_ratios.py -i merged_input.csv -o monthly_ratios
            $ # for all gridMET cells in input file for gridMET var "eto_mm"
            $ python calc_bias_ratios.py -i merged_input.csv -o monthly_ratios -gv eto_mm
            $ # for a specific gridMET cell ID for "etr_mm"
            $ python calc_bias_ratios.py -i merged_input.csv -o monthly_ratios -id 509011
            $ # to exclude any months with less than 15 days of data
            $ python calc_bias_ratios.py -i merged_input.csv -o monthly_ratios -d 15
            
        It is also possible for the user to define their own station 
        variable name if, for example, they are using station data that was
        **not** created by `PyWeatherQAQC <https://github.com/WSWUP/pyWeatherQAQC>`_.
        Let's say our station time series has ETo named as 'EO' then 
        use the ``[-sv, --station-var]`` and ``[-gv, --grid-var]`` options
        
        .. code-block:: sh

            $ python calc_bias_ratios.py -i merged_input.csv -o monthly_ratios -sv EO -gv eto_mm

        This will produce two CSV files in ``out_dir`` named 
        "eto_mm_summary_all_yrs.csv" and "eto_mm_summary_comp_all_yrs.csv". If 
        the ``[-y, --years]`` option is assigned, e.g. as '2010', then the 
        output CSVs will have '2010' suffix, i.e. 'eto_mm_summary_comp_2010.csv'
        
        For use within Python see :func:`calc_bias_ratios`.

    Note:
        If ``[-gv, --grid-var]`` command line option or ``grid_var`` 
        keyword argument is given but the station variable is left as default 
        (None), the corresponding station variable is looked up from the mapping
        dictionary in :mod:`gridwxcomp.calc_bias_ratios` named 
        ``GRIDMET_STATION_VARS``. To efficiently use climate data that was  
        **not** created by `PyWeatherQAQC <https://github.com/WSWUP/pyWeatherQAQC>`_ 
        which is where the default names are derived you can manually adjust
        ``GRIDMET_STATION_VARS`` near the top of the :mod:`gridwxcomp.calc_bias_ratios`
        submodule file. Alternatively, the gridMET and station variable names
        may be explicitly passed as command line or function arguments. 
        
    """

    # calculate monthly bias ratios and save to CSV files
    calc_bias_ratios(
        input_file_path, 
        out_dir, 
        method=method,
        grid_id_name=grid_id_name,
        grid_var=grid_var,
        station_var=station_var, 
        station_date_name=station_date_name,
        grid_date_name=grid_date_name,
        grid_ID=grid_ID, 
        day_limit=day_limit,
        comp=comp
    )

def _save_output(out_df, comp_out_df, out_dir, grid_ID, var_name, yrs):
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
        grid_ID (int, optional): depends on ``grid_ID`` argument
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
    if not grid_ID:
        out_file = OPJ(
            out_dir, 
            '{v}_summary_{y}.csv'.format(v=var_name, y=yrs)
        )
    else: 
        out_file = OPJ(
            out_dir,
            '{v}_summary_grid_{g}_{y}.csv'.format(
                v=var_name, g=grid_ID, y=yrs)
        )
    __save_update(out_df, out_file)

    # if comprehensive summary is requested save/update
    if isinstance(comp_out_df, pd.DataFrame):
        if not grid_ID:
            comp_out_file = OPJ(
                out_dir, 
                '{v}_summary_comp_{y}.csv'.format(v=var_name, y=yrs)
            )
        else:
            comp_out_file = OPJ(
                out_dir,
                '{v}_summary_comp_{g}_{y}.csv'.format(
                    v=var_name, g=grid_ID, y=yrs)
            )
        __save_update(comp_out_df, comp_out_file)
    
def calc_bias_ratios(input_path, out_dir, method='long_term_mean',
        grid_id_name='GRIDMET_ID', grid_var='etr_mm', station_var=None, 
        var_dict=None, station_date_name='date', grid_date_name='date', 
        grid_ID=None, day_limit=10, years='all', comp=True):
    """
    Read input CSV file and calculate mean monthly bias ratios between
    station to corresponding grid cells for all station and grid 
    pairs, optionally calculate ratios for a single gridcell.
    
    Arguments:
        input_path (str): path to input CSV file with matching
            station climate and grid metadata. This file is 
            created by running :func:`gridwxcomp.prep_input` followed by 
            :func:`gridwxcomp.download_gridmet_opendap`.
        out_dir (str): path to directory to save CSV files with
            monthly bias ratios of etr.
            
    Keyword Arguments:
        method (str): default 'long_term_mean'. How to calculate mean station to
            grid ratios, currently two options 'long_term_mean' takes the 
            mean of all dates for the station variable that fall in a time 
            periods, e.g. the month of January, to the mean of all paired 
            January dates in the gridded product. The other option is 
            'mean_of_annual' which calculates ratios, for each time period if
            enough paired days exist, the ratio of sums for each year in the 
            record and then takes the mean of the annual ratios. This method 
            is always used to calculate standard deviation and coefficient of 
            variation of ratios which describe interannual variation of ratios. 
        grid_id_name (str): default 'GRIDMET_ID'. Name of index/cell identifier
            for gridded dataset, only change if supplying user grid data.
        grid_var (str): default 'etr_mm'. Grid climate variable
            to calculate bias ratios.
        station_var (str): default None. Climate station variable to use
            to calculate bias ratios. If None, look up using ``grid_var`` 
            as a key to :attr:`GRIDMET_STATION_VARS` dictionary found as a 
            module attribute to :mod:`gridwxcomp.calc_bias_ratios`.
        var_dict (dict): default None. Dictionary that maps grid variable names
            to station variable names to overide gridMET and PyWeatherQaQc
            defaules used by :attr:`GRIDMET_STATION_VARS`.
        grid_ID (int): default None. Grid ID (int cell identifier) to only 
            calculate bias ratios for a single gridcell.
        day_limit (int): default 10. Threshold number of days in month
            of missing data, if less exclude month from calculations. Ignored 
            when ``method='long_term_mean'``.
        years (int or str): default 'all'. Years to use for calculations
            e.g. 2000-2005 or 2011.
        comp (bool): default True. Flag to save a "comprehensive" 
            summary output CSV file that contains station metadata and 
            statistics in addition to the mean monthly ratios.
                    
    Returns:
        None
        
    Examples:
        To use within Python for observed ET,
        
        >>> from gridwxcomp import calc_bias_ratios
        >>> input_path = 'merged_input.csv'
        >>> out_dir = 'monthly_ratios'
        >>> grid_variable = 'eto_mm'
        >>> calc_bias_ratios(input_path, out_dir, grid_var=grid_variable)
        
        To use custom station data, give the keyword argument ``station_var``, 
        e.g. if we had climate daily time series data for precipitation 
        with the column named "p" then,
        
        >>> calc_bias_ratios(input_path, out_dir, grid_var='prcp_mm', 
        >>>     station_var='p')
                
        This results in two CSV files in ``out_dir`` named 
        "prcp_mm_summary_all_yrs.csv" and "prcp_mm_summary_comp_all_yrs.csv". 

    Raises:
        FileNotFoundError: if input file is invalid or not found.
        KeyError: if the input file does not contain file paths to
            the climate station and grid time series files. This
            occurs if, for example, the :mod:`gridwxcomp.prep_input` and/or 
            :mod:`gridwxcomp.download_gridmet_opendap` scripts have not been 
            run first (if using gridMET data). Also raised if the given 
            ``grid_var``, ``station_var``, or values of ``var_dict`` kwargs 
            are invalid.
        ValueError: if the ``method`` kwarg is invalid.
    
    Note:
        If an existing summary file contains a climate station that is being 
        reprocessed its monthly bias ratios and other data will be overwritten. 
        Also, to proceed with spatial analysis scripts, the comprehensive 
        summary file must be produced using this function first. If 
        ``grid_var`` keyword argument is given but the ``station_var`` is 
        left as default (None), the corresponding station variable is looked 
        up from the mapping dictionary in :mod:`calc_bias_ratios.py` 
        named :attr:`GRIDMET_STATION_VARS`. To use climate data 
        that was  **not** created by `pyWeatherQAQC <https://github.com/WSWUP/pyWeatherQAQC>`_ 
        for station data and/or gridded data other than gridMET, which is 
        where the default names are derived, the grid and station 
        variable names need to be explicitly passed as function arguments. 
        
    """
    # ignore np runtime warnings due to calcs with nans, div by 0
    np.seterr(divide='ignore', invalid='ignore')
    # specific for standard deviation of nans
    std_warning = "Degrees of freedom <= 0 for slice"
    warnings.filterwarnings("ignore", message=std_warning)

    method_options = ('long_term_mean','mean_of_annual')
    if method not in method_options:
        raise ValueError('{} is not a valid method, use one of: {}'.format(
            method, method_options)
        )

    if var_dict is None:
        var_dict = GRIDMET_STATION_VARS

    if not var_dict.get(grid_var, None):
        print(
            'Valid grid variable names:\n',
            '\n'.join([i for i in var_dict.keys()]),
            '\n'
        )
        err_msg = 'Invalid grid variable name {}'.format(grid_var)
        raise KeyError(err_msg)

    if not os.path.isdir(out_dir):
        print('{} does not exist, creating directory'.format(out_dir))
        os.mkdir(out_dir)
    if not os.path.isfile(input_path):
        raise FileNotFoundError('Input CSV file given was invalid or not found')

    input_df = pd.read_csv(input_path)
    # get matching station variable name
    if not station_var:
        station_var = var_dict.get(grid_var)
    # If only calculating ratios for a single cell, change console message
    if grid_ID:
        single_grid_cell_msg = \
            'For grid cell ID: {g}\n'.format(g=grid_ID)
    else:
        single_grid_cell_msg = ''
    print(
        'Calculating bias ratios between climate station variable: ',
        station_var,
        '\n',
        'and grid climate variable: ',
        grid_var,
        '\n{g}'.format(g=single_grid_cell_msg)
    )
    # loop through each station and calculate monthly ratio
    for index, row in input_df.iterrows():
        if not 'STATION_FILE_PATH' in row or not 'GRID_FILE_PATH' in row:
            raise KeyError('Missing station and/or grid file paths in '+\
                           'input file. Run prep_input.py followed '+\
                           'by download_gridmet_opendap.py first.')
        # if only doing a single grid cell check for matching ID
        if grid_ID and int(grid_ID) != row[grid_id_name]:
            continue

        # load station and grid time series files
        try:
            # if time series not from PyWeatherQaQc, CSV with 'date' column
            if not row.STATION_FILE_PATH.endswith('.xlsx'):
                station_df = pd.read_csv(
                    row.STATION_FILE_PATH, parse_dates=True,
                    index_col=station_date_name
                )
                station_df.index = station_df.index.date # for joining
            # if excel file, assume PyWeatherQaQc format
            else:
                station_df = pd.read_excel(
                    row.STATION_FILE_PATH,
                    sheet_name='Corrected Data', parse_dates=True,
                    index_col=0
                )
        except:
            print('Time series file for station: ', row.STATION_ID, 
                  'was not found, skipping.')
            continue

        if not station_var in station_df.columns:
            err_msg = '{v} not found in the station file: {p}'.\
                format(v=station_var, p=row.STATION_FILE_PATH)
            raise KeyError(err_msg)
        print(
            '\nCalculating {v} bias ratios for station:'.format(v=grid_var),
            row.STATION_ID
        )
        grid_df = pd.read_csv(row.GRID_FILE_PATH, parse_dates=True, 
                                 index_col=grid_date_name)
        # merge both datasets drop missing days
        result = pd.concat(
            [
                station_df[station_var], 
                grid_df[grid_var]
            ], 
            axis=1 
        )
        result = result.reindex(grid_df.index)
        result.dropna(inplace=True)
        # make datetime index
        result.index = pd.to_datetime(result.index)
        # apply year filter
        result, years_str = parse_yr_filter(result, years, row.STATION_ID)
        # for calculating ratios with long-term means later
        orig = result.copy()
        # monthly sums and day counts for each year
        result = result.groupby([result.index.year, result.index.month])\
                .agg(['sum','count'])
        result.index.set_names(['year', 'month'], inplace=True)
        # remove totals with less than XX days
        result = result[result[grid_var,'count']>=day_limit]
        # calc mean growing season and June to August ratios with month sums
        grow_season = result.loc[
            result.index.get_level_values('month').isin([4,5,6,7,8,9]),\
                    (station_var)]['sum'].sum() / result.loc[
                result.index.get_level_values('month').isin([4,5,6,7,8,9]),\
                        (grid_var)]['sum'].sum()
        june_to_aug = result.loc[
            result.index.get_level_values('month').isin([6,7,8]), (station_var)
            ]['sum'].sum() / result.loc[result.index.get_level_values('month')\
                    .isin([6,7,8]), (grid_var)]['sum'].sum()
        ann_months = list(range(1,13))
        annual = result.loc[
            result.index.get_level_values('month').isin(ann_months),\
                    (station_var)]['sum'].sum() / result.loc[
                result.index.get_level_values('month').isin(ann_months),\
                        (grid_var)]['sum'].sum()
        ratio = pd.DataFrame(columns = ['ratio', 'count'])
        # ratio of monthly sums for each year
        ratio['ratio'] = (result[station_var,'sum'])/(result[grid_var,'sum'])
        # monthly counts and stddev
        ratio['count'] = result.loc[:,(grid_var,'count')]
        if result.empty:
            print(f'WARNING: no data for site: {row.STATION_ID}, skipping')
            continue

        # rebuild Index DateTime
        ratio['year'] = ratio.index.get_level_values('year').values.astype(int)
        ratio['month']=ratio.index.get_level_values('month').values.astype(int)
        ratio.index = pd.to_datetime(
            ratio.year*10000+ratio.month*100+15,format='%Y%m%d'
        )
        # useful to know how many years were used in addition to day counts
        start_year = ratio.year.min()
        end_year = ratio.year.max()
        counts = ratio.groupby(ratio.index.month).sum()['count']
        # get standard deviation of each years' monthly mean ratio
        stdev = {
            month: np.std(
                ratio.loc[ratio.month.isin([month]), 'ratio'].values
            ) for month in ann_months
        }
        stdev = pd.Series(stdev, name='stdev')

        # mean of monthly means of all years, can change to median or other meth
        final_ratio = ratio.groupby(ratio.index.month).mean()
        final_ratio.drop(['year', 'month'], axis=1, inplace=True)
        final_ratio['count'] = counts
        final_ratio['stdev'] = stdev
        final_ratio['cv'] = stdev / final_ratio['ratio']
        # calc mean growing season, June through August, ann stdev
        grow_season_std = np.std(
            ratio.loc[ratio.month.isin([4,5,6,7,8,9]), 'ratio'].values
        )
        june_to_aug_std = np.std(
            ratio.loc[ratio.month.isin([6,7,8]), 'ratio'].values
        )
        annual_std = np.std(
            ratio.loc[ratio.month.isin(ann_months), 'ratio'].values
        )
        # get month abbreviations in a column and drop index values
        for m in final_ratio.index:
            final_ratio.loc[m,'month'] = calendar.month_abbr[m]
        # restructure as a row with station index
        months = final_ratio.month.values
        final_ratio = final_ratio.T
        final_ratio.columns = months
        final_ratio.drop('month', inplace=True)
        # add monthy means and counts into single row dataframe
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
        final_ratio['growseason_mean'] = grow_season
        final_ratio['summer_mean'] = june_to_aug
        final_ratio['annual_mean'] = annual
        # day counts for all years in non monthly periods
        final_ratio['growseason_count'] =\
            counts.loc[counts.index.isin([4,5,6,7,8,9])].sum()
        final_ratio['summer_count'] =\
            counts.loc[counts.index.isin([6,7,8])].sum()
        final_ratio['annual_count'] =\
            counts.loc[counts.index.isin(ann_months)].sum()
        # assign stdev, coef. var. 
        final_ratio['growseason_stdev'] = grow_season_std
        final_ratio['summer_stdev'] = june_to_aug_std
        final_ratio['annual_stdev'] = annual_std
        # coefficient of variation
        final_ratio['growseason_cv'] = grow_season_std / grow_season
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

        # save grid ID for merging with input table
        final_ratio[grid_id_name] = row[grid_id_name]    
        final_ratio = final_ratio.merge(input_df, on=grid_id_name)

        # long term mean station to mean grid ratio calc as opposed to mean of
        # annual ratios- default less bias potential
        if method == 'long_term_mean':
            month_means = orig.groupby(orig.index.month).mean()
            month_means['month'] = month_means.index
            for m in month_means.index:
                month_means.loc[m,'month'] = f'{calendar.month_abbr[m]}_mean'
            month_means.set_index('month', inplace=True)
            month_means['ratios'] =\
                month_means[station_var] / month_means[grid_var]

            long_term = month_means.drop([station_var, grid_var],1).T
            # non-monthly periods long-term mean to mean ratios
            grow_season = orig.loc[orig.index.month.isin([4,5,6,7,8,9])]
            long_term['growseason_mean'] =\
                grow_season[station_var].mean() / grow_season[grid_var].mean()
            summer_season = orig.loc[orig.index.month.isin([6,7,8])]
            long_term['summer_mean'] =\
                summer_season[station_var].mean()/summer_season[grid_var].mean()
            long_term['annual_mean'] =\
                orig[station_var].mean() / orig[grid_var].mean()
            # overwrite only mean ratios (keep stats from mean of annual ratios)
            overwrite = long_term.columns.intersection(final_ratio.columns)
            final_ratio[overwrite] = long_term[overwrite].values


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
        grow_thresh = 65
        sum_thresh = 35
        ann_thresh = 125
        
        if final_ratio.at[0,'summer_count'] < sum_thresh:
            print('WARNING: less than:', sum_thresh, 'days in summer period',
                 '\nfor station:',row.STATION_ID,'assigning -999 for all stats')
            cols = [col for col in final_ratio.columns if 
                    'summer_' in col and '_count' not in col]
            final_ratio.loc[:,cols] = np.nan
        
        if final_ratio.at[0,'growseason_count'] < grow_thresh:
            print('WARNING: less than:',grow_thresh,'days in growing season',
                 '\nfor station:',row.STATION_ID,'assigning -999 for all stats')
            cols = [col for col in final_ratio.columns if 
                    'growseason_' in col and '_count' not in col]
            final_ratio.loc[:,cols] = np.nan
            
        if final_ratio.at[0,'annual_count'] < ann_thresh:
            print('WARNING: less than:',ann_thresh,'days in annual period',
                 '\nfor station:',row.STATION_ID,'assigning -999 for all stats')
            cols = [col for col in final_ratio.columns if 
                    'annual_' in col and '_count' not in col]
            final_ratio.loc[:,cols] = np.nan

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
        _save_output(out, comp_out, out_dir, grid_ID, grid_var, years_str)

    print(
        '\nSummary file(s) for bias ratios saved to: \n', 
         os.path.abspath(out_dir)
         )    
    
def arg_parse():
    """
    Command line usage of calc_bias_ratios.py which calculates monthly bias 
    ratios between station climate and grid cells that correspond with 
    each other geographically. Saves data to CSV files in the given output 
    directory. If run later with new station data, bias ratios for new 
    stations will be appended to existing output summary CSV.
    """
    parser = argparse.ArgumentParser(
        description=arg_parse.__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    optional = parser._action_groups.pop() # optionals listed second
    required = parser.add_argument_group('required arguments')
    required.add_argument(
        '-i', '--input', metavar='PATH', required=True,
        help='Input CSV file of merged climate/grid data that '+\
            'was created by running prep_input.py and '+\
            'download_gridmet_opendap.py')
    required.add_argument(
        '-o', '--out', metavar='PATH', required=True,
        help='Output directory to save CSV files containing bias ratios')
    optional.add_argument('-meth', '--method', metavar='', required=False, 
        default='long_term_mean', help='ratio calc method "long_term_mean" or'+\
            '"mean_of_annual"')
    optional.add_argument('-gin', '--grid-id-name', metavar='', required=False, 
        default='GRIDMET_ID', help='Name of gridcell identifier if not using '+\
            'gridMET grid')
    optional.add_argument(
        '-y', '--years', metavar='', required=False, default='all',
        help='Years to use, single or range e.g. 2018 or 1995-2010')
    optional.add_argument(
        '-gv', '--grid-var', metavar='', required=False, default='etr_mm',
        help='Grid variable name for bias ratio calculation')
    optional.add_argument(
        '-sv', '--station-var', metavar='', required=False, default=None,
        help='Station variable name for bias ratio calculation')
    optional.add_argument(
        '-sdn', '--station-date-name', metavar='',required=False,default='date',
        help='Date column name in station time series files if not using '+\
            'gridMET.')
    optional.add_argument(
        '-gdn', '--grid-date-name', metavar='', required=False, default='date',
        help='Date column name in grid time series files if not using gridMET.')
    optional.add_argument(
        '-id', '--grid-id', metavar='', required=False, default=None,
        help='Optional grid ID to calculate bias ratios for a single '+\
            'gridcell')
    optional.add_argument('-d', '--day-limit', metavar='', required=False, 
        default=10, help='Number of days of valid data per month to '+\
            'include it in bias correction calculation.')
    optional.add_argument('-c', '--comprehensive', required=False, 
        default=True, action='store_false', dest='comprehensive', 
        help='Flag, if given, to NOT save comprehensive summary file with '+\
            'extra metadata and statistics with the suffix "_comp"')
#    parser.add_argument(
#        '--debug', default=logging.INFO, const=logging.DEBUG,
#        help='Debug level logging', action="store_const", dest="loglevel")
    parser._action_groups.append(optional)# to avoid optionals listed first
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = arg_parse()

    main(
        input_file_path=args.input, 
        out_dir=args.out,
        method=args.method,
        grid_id_name=args.grid_id_name,
        grid_var=args.grid_var, 
        station_var=args.station_var,
        station_date_name=args.station_date_name,
        grid_date_name=args.grid_date_name,
        grid_ID=args.grid_id, 
        day_limit=args.day_limit,
        years=args.years, 
        comp=args.comprehensive
    )
