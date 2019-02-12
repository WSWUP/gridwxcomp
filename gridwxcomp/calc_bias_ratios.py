# -*- coding: utf-8 -*-
"""
Calculate monthly bias ratios of variables from climate station 
to overlapping gridMET cells. 

Input file for this module must first be created by running 
:mod:`prep_input.py` followed by :mod:`download_gridmet_ee.py`. 

Attributes:
    GRIDMET_STATION_VARS (:obj:`dict`): mapping dictionary with gridMET
        variable names as keys and station variable names as values.
        Used to determine which station variable to calculate bias
        ratios according to the given gridMET variable. 
        
Note:
    ``GRIDMET_STATION_VARS`` can be manually adjusted, e.g. new pairs
    can be made or removed to efficiently use :mod:`gridwxcomp` on custom
    station data that was **not** created by 
    `pyWeatherQAQC <https://github.com/DRI-WSWUP/pyWeatherQAQC>`_.
    
"""

import os
import calendar
import argparse

import pandas as pd
import numpy as np

# keys = gridMET variable name
# values = climate station variable name
GRIDMET_STATION_VARS = {
    'u2_ms' : 'Windspeed (m/s)',
    'tmin_c' : 'TMin (C)',
    'tmax_c' : 'TMax (C)',
    'srad_wm2' : 'Rs (w/m2)',
    'ea_kpa' : 'Vapor Pres (kPa)',
    'prcp_mm' : 'Precip (mm)',
    'etr_mm' : 'Calc_ETr (mm)',
    'eto_mm' : 'Calc_ETo (mm)'
}

OPJ = os.path.join

def main(input_file_path, out_dir, gridmet_var='etr_mm', station_var=None,
         gridmet_id=None, comp=True):
    """
    Calculate monthly bias ratios between station climate and gridMET
    cells that correspond with each other geographically. Saves data
    to CSV files in the given output directory. If run later with
    new station data, bias ratios for new stations will be appended
    to existing output summary CSV.
    
    Arguments:
        input_file_path (str): path to input CSV file containing
            paired station/gridMET metadata. This file is 
            created by running :mod:`prep_input.py` followed by 
            :mod:`download_gridmet_ee.py`.
        out_dir (str): path to directory to save CSV files with
            monthly bias ratios of etr.
            
    Keyword Arguments:
        gridmet_var (str): default 'etr_mm'. GridMET climate variable
            to calculate bias ratios.
        station_var (str): default None. Climate station variable to use
            to calculate bias ratios. If None, look up using ``gridmet_var`` 
            as a key to ``GRIDMET_STATION_VARS`` dictionary found as a 
            module attribute to :mod:`calc_bias_ratios.py`.
        gridmet_ID (str): optional gridMET ID number if user wants
            to only calculate bias ratios for a single gridMET cell.
        comp (bool): default True. Save a "comprehensive" summary
            output CSV file that contains additional station metadata
            and statistics in addition to the mean monthly ratios.

    Returns:
        None

    Examples:
        From the command line interface,

        .. code::
            $ # for all gridMET cells in input file for gridMET var "etr_mm" (default)
            $ python calc_bias_ratios.py -i merged_input.csv -o monthly_ratios
            $ # for all gridMET cells in input file for gridMET var "eto_mm"
            $ python calc_bias_ratios.py -i merged_input.csv -o monthly_ratios -gv eto_mm
            $ # for a specific gridMET cell for "etr_mm"
            $ python calc_bias_ratios.py -i merged_input.csv -o monthly_ratios -id 509011
            
        It is also possible for the user to define their own station 
        variable name if, for example, they are using station data that was
        **not** created by `pyWeatherQAQC <https://github.com/DRI-WSWUP/pyWeatherQAQC>`_.
        Let's say our station time series has ETo named as 'EO' then 
        use the ``[-sv, --station-var]`` and ``[-gv, --gridmet-var]`` options
        
        .. code::
            $ python calc_bias_ratios.py -i merged_input.csv -o monthly_ratios -sv EO -gv eto_mm

        This will produce two CSV files in ``out_dir`` named "eto_mm_summary.csv"
        and "eto_mm_summary_comp.csv". 
        
        For use within Python see :func:`calc_bias_ratios`.

    Note:
        If ``[-gv, --gridmet-var]`` command line option or ``gridmet_var`` 
        keyword argument is given but the station variable is left as default 
        (None), the corresponding station variable is looked up from the mapping 
        dictionary in :mod:`calc_bias_ratios.py` named ``GRIDMET_STATION_VARS``.
        To efficiently use climate data that was  **not** created by 
        `pyWeatherQAQC <https://github.com/DRI-WSWUP/pyWeatherQAQC>`_ which
        is where the default names are derived we recommend manually adjusting
        ``GRIDMET_STATION_VARS`` near the top of the :mod:`calc_bias_ratios.py`
        submodule file. Alternatively, the gridMET and station variable names
        need to be explicitly passed as command line or function arguments. 
        
    """

    # calculate monthly bias ratios and save to CSV files
    calc_bias_ratios(
        input_file_path, 
        out_dir, 
        gridmet_var=gridmet_var,
        station_var=station_var, 
        gridmet_ID=gridmet_id, 
        comp=comp
    )

def _save_output(out_df, comp_out_df, out_dir, gridmet_ID, var_name):
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
        gridmet_ID (int, optional): depends on ``gridmet_ID`` argument
            passed to :func:`calc_bias_ratios`. If not None (default)
            then save summary files for stations that correspond with
            the given gridMET ID with the suffix "_X" where X is the
            gridMET ID value.
        var_name (str): name of gridMET variable that is being processed.
    
    Returns:
        None       
        
    """
        
    def __save_update(out_df, out_file):
        """
        Helper function that is reused for both short and long summary
        files.        
        """
        # if short file exists add/overwrite row for station
        if os.path.isfile(out_file):
            existing_df = pd.read_csv(out_file, index_col='STATION_ID')
            if not out_df.index.values[0] in existing_df.index.values:
                out_df = pd.concat([existing_df, out_df])
                out_df.to_csv(out_file, index=True)
            # overwrite if station is in existing, could change to
            # allow for duplicates if values are different
            else:
                existing_df.loc[out_df.index.values[0], :] =\
                    out_df.loc[out_df.index.values[0]]
                existing_df.to_csv(out_file, index=True)
        else:
            out_df.to_csv(out_file, index=True)
            
    # save or update short and comprehensive summary files
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)   
        
    # save/update short summary file or update existing with new station
    if not gridmet_ID:
        out_file = OPJ(
            out_dir, 
            '{v}_summary.csv'.format(v=var_name)
        )
    else: 
        out_file = OPJ(
            out_dir,
            '{v}_summary_grid_{g}.csv'.format(v=var_name, g=gridmet_ID)
        )
    __save_update(out_df, out_file)

    # if comprehensive summary is requested save/update
    if isinstance(comp_out_df, pd.DataFrame):
        if not gridmet_ID:
            comp_out_file = OPJ(
                out_dir, 
                '{v}_summary_comp.csv'.format(v=var_name)
            )
        else:
            comp_out_file = OPJ(
                out_dir,
                '{v}_summary_comp_{g}.csv'.format(v=var_name, g=gridmet_ID)
            )
        __save_update(comp_out_df, comp_out_file)
    
def calc_bias_ratios(input_path, out_dir, gridmet_var='etr_mm', 
                     station_var=None, gridmet_ID=None, comp=True):
    """
    Read input CSV file and calculate mean monthly bias ratios between
    station to corresponding gridMET cells for all station and gridMET 
    pairs, optionally calculate ratios for a single gridMET cell.
    
    Arguments:
        input_path (str): path to input CSV file with matching
            station climate and gridMET metadata. This file is 
            created by running prep_input.py followed by 
            download_gridmet_ee.py.
        out_dir (str): path to directory to save CSV files with
            monthly bias ratios of etr.
            
    Keyword Arguments:
        gridmet_var (str): default 'etr_mm'. GridMET climate variable
            to calculate bias ratios.
        station_var (str): default None. Climate station variable to use
            to calculate bias ratios. If None, look up using ``gridmet_var`` 
            as a key to ``GRIDMET_STATION_VARS`` dictionary found as a 
            module attribute to :mod:`calc_bias_ratios.py`.
        gridmet_ID (int): default None. GridMET ID number if user wants
            to only calculate bias ratios for a single gridMET cell.
        comp (bool): default True. Flag to save a "comprehensive" 
            summary output CSV file that contains additional station 
            metadata and statistics in addition to the mean monthly ratios.
                    
    Returns:
        None
        
    Examples:
        To use within Python for observed ET,
        
        >>> from gridwxcomp import calc_bias_ratios
        >>> input_path = 'merged_input.csv'
        >>> out_dir = 'monthly_ratios'
        >>> gridmet_variable = 'eto_mm'
        >>> calc_bias_ratios(input_path, out_dir, gridmet_var=gridmet_variable)
        
        To use custom station data, give the keyword argument ``station_var``, 
        e.g. if we had climate daily time series data for precipitation 
        with the column named "p" then,
        
        >>> calc_bias_ratios(input_path, out_dir, 
                gridmet_var='prcp_mm', station_var='p')
                
        This results in two CSV files in ``out_dir`` named "prcp_mm_summary.csv"
        and "prcp_mm_summary_comp.csv". 
        
        For command line examples see :func:`main`.
        
    Raises:
        FileNotFoundError: if input file is invalid or not found.
        KeyError: if the input file does not contain file paths to
            the climate station and gridMET time series files. This
            occurs if, for example, the :mod:`prep_input.py` and/or 
            the :mod:`download_gridmet_ee.py` scripts have not been 
            run first. Also raised if the given ``gridmet_var`` or 
            ``station_var`` kwargs are invalid.
    
    Note:
        If an existing summary file contains a climate station that
        is being reprocessed its monthly bias ratios and other data
        will be overwritten. Also, to proceed with spatial analysis
        scripts, the comprehensive summary file must be produced 
        using this function first (default). If ``gridmet_var`` 
        keyword argument is given but the ``station_var`` is left as 
        default (None), the corresponding station variable is looked 
        up from the mapping dictionary in :mod:`calc_bias_ratios.py` 
        named ``GRIDMET_STATION_VARS``. To efficiently use climate data 
        that was  **not** created by `pyWeatherQAQC <https://github.com/DRI-WSWUP/pyWeatherQAQC>`_ 
        which is where the default names are derived we recommend 
        manually adjusting ``GRIDMET_STATION_VARS`` near the top of 
        the :mod:`calc_bias_ratios.py` submodule file. Alternatively, 
        the gridMET and station variable names need to be explicitly 
        passed as function arguments. 
        
    """
    if not GRIDMET_STATION_VARS.get(gridmet_var, None):
        print(
            'Valid gridMET variable names:\n',
            '\n'.join([i for i in GRIDMET_STATION_VARS.keys()]),
            '\n'
        )
        raise KeyError('Invalid gridMET variable name {}'.format(gridmet_var))
        
    if not os.path.isdir(out_dir):
        print('{} does not exist, creating directory'.format(out_dir))
        os.mkdir(out_dir)
    if not os.path.isfile(input_path):
        raise FileNotFoundError('Input CSV file given was invalid or not found')

    input_df = pd.read_csv(input_path)
    # get matching station variable name
    if not station_var:
        station_var = GRIDMET_STATION_VARS.get(gridmet_var)
    # If only calculating ratios for a single cell, change console message
    if gridmet_ID:
        single_gridmet_cell_msg = \
            'For gridmet cell ID: {g}\n'.format(g=gridmet_ID)
    else:
        single_gridmet_cell_msg = ''
    print(
        'Calculating bias ratios between climate station variable: ',
        station_var,
        '\n',
        'and gridMET climate variable: ',
        gridmet_var,
        '\n{g}'.format(g=single_gridmet_cell_msg)
    )
    # loop through each station and calculate monthly ratio
    for index, row in input_df.iterrows():
        if not 'STATION_FILE_PATH' in row or not 'GRIDMET_FILE_PATH' in row:
            raise KeyError('Missing station and/or gridMET file paths in '+\
                           'input file. Run prep_input.py followed '+\
                           'by download_gridmet_ee.py first.')
        # if only doing a single gridMET cell check for matching ID
        if gridmet_ID and int(gridmet_ID) != row.GRIDMET_ID:
            continue

        # load station and gridMET time series files
        try:
            # if station files not from pyWeatherQAQC this needs changed
            station_df = pd.read_excel(row.STATION_FILE_PATH,
                                       sheet_name='Corrected Data')
        except:
            print('Station time series file: ', row.STATION_FILE_PATH, 
                  '\nwas not found, skipping.')
            continue
            
        if not station_var in station_df.columns:
            raise KeyError('{v} not found in the station file: \n{p}'.\
                           format(v=station_var, p=row.STATION_FILE_PATH))
        print(
             'Calculating {v} bias ratios for station: '.format(v=gridmet_var), 
             row.STATION_ID
             )
        gridmet_df = pd.read_csv(row.GRIDMET_FILE_PATH, parse_dates=True, 
                                 index_col='date')
        # merge both datasets drop missing days
        result = pd.concat([station_df[station_var], 
                            gridmet_df[gridmet_var]], axis=1, 
                           join_axes=[station_df.index])
        result.dropna(inplace=True)
        # calculate monthy mean, median, and counts of both datasets
        result = result.groupby(result.index.month).\
                agg({
                    station_var:['sum','median','mean','count'],
                     gridmet_var:['sum','median','mean','count']
                     })
        # calculate ratios of monthly mean and median values and day counts
        result[('ratio','mean_ratio')] =\
               result.loc[:,(station_var,'sum')]\
               / result.loc[:,(gridmet_var,'sum')]
        result[('ratio','median_ratio')] =\
               result.loc[:,(station_var,'median')]\
               / result.loc[:,(gridmet_var,'median')]
        result[('ratio','count_days')] = result.loc[:,(gridmet_var,'count')]
        # drop columns and multiindex
        result.drop([gridmet_var,station_var], axis=1, inplace=True)
        result.columns = result.columns.droplevel(0)
        # calc mean growing season and June to August ratios
        result['April_to_oct_mean'] = result.loc[3:10,'mean_ratio'].mean()
        result['June_to_aug_mean'] = result.loc[5:8,'mean_ratio'].mean()
        result['Annual_mean'] = result.loc[:,'mean_ratio'].mean()
        # get month abbreviations in a column and drop index values
        for m in result.index:
            result.loc[m,'month'] = calendar.month_abbr[m]
        # save GRIDMET_ID for merging with input table
        result['GRIDMET_ID'] = row.GRIDMET_ID    
        result = result.merge(input_df, on='GRIDMET_ID')

        # round numeric columns
        result = result.round({
            'month': 0,
            'mean_ratio': 3,
            'median_ratio': 3,
            'count_days': 0,
            'April_to_oct_mean': 2,
            'June_to_aug_mean': 2,
            'Annual_mean': 2,
            'LAT': 10,
            'LON': 10,
            'ELEV_M': 0,
            'ELEV_FT': 0,
            'STATION_LAT': 10,
            'STATION_LON': 10,
            'STATION_ELEV_M': 0
        })

        # pivot table on monthly statistics, concatenate, merge for comp
        df_mean = result.pivot(
            index='STATION_ID', 
            columns='month', 
            values='mean_ratio'
        )
        df_mean.columns = [c+'_mean' for c in df_mean.columns]

        df_median = result.pivot(
            index='STATION_ID', 
            columns='month', 
            values='median_ratio'
        )
        df_median.columns = [c+'_median' for c in df_median.columns]

        df_count = result.pivot(
            index='STATION_ID', 
            columns='month', 
            values='count_days'
        )
        df_count.columns = [c+'_count' for c in df_count.columns]

        out = df_mean    
        # concat only growing season and annual ratios to short summary
        tmp = pd.DataFrame(result.loc[0,['April_to_oct_mean','Annual_mean']]).T
        tmp.index = result.STATION_ID.unique()
        out = pd.concat([out,tmp], axis=1)

        if comp:
            out['GRIDMET_ID'] = result.GRIDMET_ID.unique()
            # build comprehensive output summary 
            comp_out = pd.concat([df_mean,df_median,df_count], axis=1)
            comp_out['GRIDMET_ID'] = result.GRIDMET_ID[0]
            # add other non-monthly data by merge
            cols = [
                c for c in result.columns if c not in
                [
                    'month', 
                    'mean_ratio',
                    'median_ratio',
                    'count_days'
            ]]

            comp_out = comp_out.merge(
                    result[cols], 
                    on='GRIDMET_ID'
                ).drop_duplicates()
            comp_out.set_index('STATION_ID', inplace=True)

            # no longer need GRIDMET_ID in short summary 
            out.drop(columns='GRIDMET_ID', inplace=True)
        # if comp False
        else:
            comp_out = comp
        # save output depending on options
        _save_output(out, comp_out, out_dir, gridmet_ID, gridmet_var)
        
    print(
        '\nSummary file(s) for bias ratios saved to: \n', 
         os.path.abspath(out_dir)
         )    

def arg_parse():
    """
    Command line usage of calc_bias_ratios.py which calculates monthly bias 
    ratios between station climate and gridMET cells that correspond with 
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
        help='Input CSV file of merged climate/gridMET data that '+\
             'was created by running prep_input.py and download_gridmet_ee.py')
    required.add_argument(
        '-o', '--out', metavar='PATH', required=True,
        help='Output directory to save CSV files containing bias ratios')
    optional.add_argument(
        '-gv', '--gridmet-var', metavar='', required=False, default='etr_mm',
        help='GridMET variable name for bias ratio calculation')
    optional.add_argument(
        '-sv', '--station-var', metavar='', required=False, default=None,
        help='Station variable name for bias ratio calculation')
    optional.add_argument(
        '-id', '--gridmet-id', metavar='', required=False, default=None,
        help='Optional gridMET ID to calculate bias ratios for a single '+\
             'gridMET cell')
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

    main(input_file_path=args.input, out_dir=args.out,
         gridmet_var=args.gridmet_var, station_var=args.station_var,
         gridmet_id=args.gridmet_id, comp=args.comprehensive)
