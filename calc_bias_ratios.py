"""
"""

import pandas as pd
import numpy as np
import os
import calendar
import argparse

def main(input_file_path, out_dir, gridmet_id=None):
    """
    Calculate monthly bias ratios between station climate and gridMET
    cells that correspond with each other geographically. Save data
    to CSV files in the given output directory. Does not overwrite 
    previously calculated data.

    Example:
        From the command line interface,

        .. code::
            $ # for all gridMET cells in input file
            $ python calc_bias_ratios.py -i merged_input.csv -o monthly_ratios
            $ # for a specific gridMET cell
            $ python calc_bias_ratios.py -i merged_input.csv -o monthly_ratios -id 509011

    Raises:
        FileNotFoundError: if input file is invalid or not found.
    """
    if not os.path.isfile(input_file_path):
        raise FileNotFoundError('Input CSV file given was invalid or not found')
    if not os.path.isdir(out_dir):
        print('{} does not exist, creating directory'.format(out_dir))
        os.mkdir(out_dir)
    # calculate monthly bias ratios and save to CSV files
    calc_bias_ratios(input_file_path,out_dir,gridmet_ID=gridmet_id)

def calc_bias_ratios(input_path, out_dir, gridmet_ID=None):
    """
    Read input CSV file and calculate monthly bias ratios between
    station and corresponding gridMET cells for all station
    and gridMET time series, optionally calculate ratios for a 
    single gridMET cell.
    
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
            
    Raises:
        KeyError: if the input file does not contain file paths to
            the climate station and gridMET time series files. This
            occurs if, for example, the prep_input.py and/or the
            download_gridmet_ee.py scripts have not been run first.
    """
    input_df = pd.read_csv(input_path)
    # loop through each station and calculate monthly ratio
    for index, row in input_df.iterrows():
        if not 'STATION_FILE_PATH' in row or not 'GRIDMET_FILE_PATH' in row:
            raise KeyError('Missing station and/or gridMET file paths in '+\
                          'input file. Please run prep_input.py followed '+\
                          'by download_gridmet_ee.py first.')
        # if only doing a single gridMET cell check for matching ID
        if gridmet_ID and int(gridmet_ID) != row.GRIDMET_ID:
            continue
            
        # load station and gridMET time series files
        station_df = pd.read_excel(row.STATION_FILE_PATH,
                                     sheet_name='Corrected Data')
        gridmet_df = pd.read_csv(row.GRIDMET_FILE_PATH, parse_dates=True, 
                                 index_col='date')
        # merge both datasets drop missing days
        result = pd.concat([station_df['Calc_ETr (mm)'], 
                            gridmet_df['etr_mm']], axis=1, 
                           join_axes=[station_df.index])
        result.dropna(inplace=True)
        # calculate monthy mean, median, and counts of both datasets
        result = result.groupby(result.index.month).\
                agg({'Calc_ETr (mm)':['sum','median','mean','count'],
                                      'etr_mm':['sum','median','mean','count']})
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
        # get month abbreviations in a column and drop index values
        for m in result.index:
            result.loc[m,'month'] = calendar.month_abbr[m]
        # save GRIDMET_ID for merging with input table
        result['GRIDMET_ID'] = row.GRIDMET_ID    
        result = result.merge(input_df, on='GRIDMET_ID')
        # reorder columns and save
        result = result.reindex(columns=['month','mean_ratio',
               'median_ratio','count_days','April_to_oct_mean',
               'June_to_aug_mean','GRIDMET_ID','LAT','LON', 'ELEV_M',
               'FID','STATION_LAT','STATION_LON','STATION_ELEV_M',
               'STATION_FILE_PATH','GRIDMET_FILE_PATH'])
        # round numeric columns
        result = result.round({'month': 0,
                               'mean_ratio': 3,
                               'median_ratio': 3,
                               'count_days': 0,
                               'April_to_oct_mean': 2,
                               'June_to_aug_mean': 2,
                               'GRIDMET_ID': 0,
                               'LAT': 10,
                               'LON': 10,
                               'ELEV_M': 0,
                               'FID': 0,
                               'STATION_LAT': 10,
                               'STATION_LON': 10,
                               'STATION_ELEV_M': 0

        })
        out_file = os.path.join(out_dir, '{grid_id}_ratios.csv'.\
                                   format(grid_id=row.GRIDMET_ID))
        # if CSV ratio file already exists for gridMET cell
        # and station id (FID) is not already there, append data
        if os.path.isfile(out_file):
            existing_df = pd.read_csv(out_file)
            if not row.FID in set(existing_df.FID):
                result = pd.concat([existing_df, result])
        # if file does not exist or specific station is not there, save
        else:
            result.to_csv(out_file, index=False)

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
        '-id', '--gridmet_id', metavar='PATH', required=False,
        help='Optional gridMET ID to calculate bias ratios for a single '+\
             'gridMET cell')
#    parser.add_argument(
#        '--debug', default=logging.INFO, const=logging.DEBUG,
#        help='Debug level logging', action="store_const", dest="loglevel")
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = arg_parse()

    main(input_file_path=args.input, out_dir=args.out,
         gridmet_id=args.gridmet_id)
