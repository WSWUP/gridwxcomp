# -*- coding: utf-8 -*-
"""
Calculate monthly bias ratios of climate station to gridMET cell 
for each station given or a select gridMET cell, write 
summary statistics to CSV file(s).
"""

import pandas as pd
import numpy as np
import os
import calendar
import argparse

def main(input_file_path, out_dir, gridmet_id=None, comp=False):
    """
    Calculate monthly bias ratios between station climate and gridMET
    cells that correspond with each other geographically. Save data
    to CSV files in the given output directory. Does not overwrite 
    previously calculated data.
    
    Arguments:
        input_file_path (str): path to input CSV file containing
            paired station/gridMET metadata. This file is 
            created by running :mod:`prep_input.py` followed by 
            :mod:`download_gridmet_ee.py`.
        out_dir (str): path to directory to save CSV files with
            monthly bias ratios of etr.
            
    Keyword Arguments:
        gridmet_ID (str): optional gridMET ID number if user wants
            to only calculate bias ratios for a single gridMET cell.
        comp (bool): optional flag to save a "comprehensive" summary
            output CSV file that contains additional station metadata
            and statistics in addition to the mean monthly ratios.

    Returns:
        None

    Example:
        From the command line interface,

        .. code::
            $ # for all gridMET cells in input file
            $ python calc_bias_ratios.py -i merged_input.csv -o monthly_ratios
            $ # for a specific gridMET cell
            $ python calc_bias_ratios.py -i merged_input.csv -o monthly_ratios -id 509011
            $ # for all gridMET cells and output comprehensive summary
            $ python calc_bias_ratios.py -i merged_input.csv -o monthly_ratios -c 

        To use within Python for all station data and comprehensive output,
        
        >>> from calc_bias_ratios import calc_bias_ratios
        >>> input_path = 'merged_input.csv'
        >>> out_dir = 'test_out_ratios'
        >>> comp = True
        >>> calc_bias_ratios(input_path, out_dir, comp=comp)
        
        This will produce two CSV files in ``out_dir`` named summary.csv
        and summary_comp.csv.
        
    """

    if not os.path.isdir(out_dir):
        print('{} does not exist, creating directory'.format(out_dir))
        os.mkdir(out_dir)
    # calculate monthly bias ratios and save to CSV files
    calc_bias_ratios(input_file_path, out_dir, gridmet_ID=gridmet_id,
            comp=comp)

def _save_output(out_df, comp_out_df, out_dir, gridmet_ID):
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
        out_file = os.path.join(out_dir, 'summary.csv')
    else: 
        out_file = os.path.join(
                out_dir, 
                'summary_grid_{}.csv'.format(gridmet_ID)
                )
    __save_update(out_df, out_file)

    # if comprehensive summary is requested save/update
    if isinstance(comp_out_df, pd.DataFrame):
        if not gridmet_ID:
            comp_out_file = os.path.join(out_dir, 'summary_comp.csv')
        else:
            comp_out_file = os.path.join(
                    out_dir, 
                    'summary_comp_{}.csv'.format(gridmet_ID)
                    )
        __save_update(comp_out_df, comp_out_file)
    
    
def calc_bias_ratios(input_path, out_dir, gridmet_ID=None, comp=False):
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
        gridmet_ID (str): optional gridMET ID number if user wants
            to only calculate bias ratios for a single gridMET cell.
        comp (bool): optional flag to save a "comprehensive" summary
            output CSV file that contains additional station metadata
            and statistics in addition to the mean monthly ratios.
                    
    Returns:
        None
        
    Raises:
        FileNotFoundError: if input file is invalid or not found.
        KeyError: if the input file does not contain file paths to
            the climate station and gridMET time series files. This
            occurs if, for example, the :mod:`prep_input.py` and/or 
            the :mod:`download_gridmet_ee.py` scripts have not been 
            run first.
    
    Note:
        If an existing summary file contains a climate station that
        is being processed its monthly bias ratios and other data
        will be overwritten. Also, to proceed with spatial analysis
        scripts, the "-c" or comprehensive keyword argument must be 
        True.
        
    """
    if not os.path.isfile(input_path):
        raise FileNotFoundError('Input CSV file given was invalid or not found')

    input_df = pd.read_csv(input_path)
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
            station_df = pd.read_excel(row.STATION_FILE_PATH,
                                       sheet_name='Corrected Data')
        except:
            print('Station time series file: ', row.STATION_FILE_PATH, 
                  '\nwas not found, skipping.')
            continue
        print(
             'Calculating bias ratios for station: ', 
             row.STATION_ID
             )
        gridmet_df = pd.read_csv(row.GRIDMET_FILE_PATH, parse_dates=True, 
                                 index_col='date')
        # merge both datasets drop missing days
        result = pd.concat([station_df['Calc_ETr (mm)'], 
                            gridmet_df['etr_mm']], axis=1, 
                           join_axes=[station_df.index])
        result.dropna(inplace=True)
        # calculate monthy mean, median, and counts of both datasets
        result = result.groupby(result.index.month).\
                agg({
                    'Calc_ETr (mm)':['sum','median','mean','count'],
                     'etr_mm':['sum','median','mean','count']
                     })
        # calculate ratios of monthly mean and median values and day counts
        result[('ratio','mean_ratio')] =\
               result.loc[:,('Calc_ETr (mm)','sum')]\
               / result.loc[:,('etr_mm','sum')]
        result[('ratio','median_ratio')] =\
               result.loc[:,('Calc_ETr (mm)','median')]\
               / result.loc[:,('etr_mm','median')]
        result[('ratio','count_days')] = result.loc[:,('etr_mm','count')]
        # drop columns and multiindex
        result.drop(['etr_mm','Calc_ETr (mm)'],axis=1, inplace=True)
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
        _save_output(out, comp_out, out_dir, gridmet_ID)
        
    print(
        '\nSummary file(s) for bias ratios saved to: ', 
         os.path.abspath(out_dir)
         )

def arg_parse():
    """
    Parse command line arguments for calculating monthly bias ratios 
    between climate station and gridMET etr estimates.
    """
    parser = argparse.ArgumentParser(
        description='Calculate monthly bias ratios for gridMET cells.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '-i', '--input', metavar='PATH', required=True,
        help='Input CSV file of merged climate/gridMET data that '+\
             'was created by running prep_input.py and download_gridmet_ee.py')
    parser.add_argument(
        '-o', '--out', metavar='PATH', required=True,
        help='Output directory to save CSV files containing bias ratios '+\
             'for gridMET cells')
    parser.add_argument(
        '-id', '--gridmet_id', metavar='INTEGER', required=False,
        help='Optional gridMET ID to calculate bias ratios for a single '+\
             'gridMET cell')
    parser.add_argument('-c', '--comprehensive', required=False, default=False, 
        action='store_true', dest='comprehensive', help='Optional flag '+\
             'to save a summary file with bias ratios and extra metadata '+\
             'and statistics with the suffix "_comp" (default=False when '+\
             'flag ommitted)')
#    parser.add_argument(
#        '--debug', default=logging.INFO, const=logging.DEBUG,
#        help='Debug level logging', action="store_const", dest="loglevel")
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = arg_parse()

    main(input_file_path=args.input, out_dir=args.out,
         gridmet_id=args.gridmet_id, comp=args.comprehensive)
