# -*- coding: utf-8 -*-
"""
Create shapefile containing monthly mean bias ratios for all stations 
found in an input comprehensive summary file. The input file is first 
created by :mod:`calc_bias_ratios`.

TODO:
    snap to gridmet cells, add spatial interpolation, docs, logging 
"""

import os
import argparse
import pandas as pd
from shapely.geometry import Point, mapping
from fiona import collection
from fiona.crs import from_epsg

OPJ = os.path.join

def main(input_file_path):
    """
    Create shapefile of monthly mean bias ratios from comprehensive
    CSV file created by :mod:`calc_bias_ratios`.
    
    Arguments:
        in_path (str): path to summary_comp.CSV file containing monthly
            bias ratios, lat, long, and other data. Shapefile "summary.shp"
            is saved to parent directory of ``in_path`` under a subdirectory
            named "spatial".
            
    Returns:
        None

    Example:
        From the command line interface,

        .. code::
            $ python make_spatial.py -i monthly_ratios/summary_comp.csv

        Or from within Python,

        >>> from make_spatial import make_spatial
        >>> input_path = 'monthly_ratios/summary_comp.csv'
        >>> make_spatial(input_path)
        
        This will produce a subdirectory under "monthly_ratios" named
        "spatial" where the shapefile will be saved as "summary.shp".

    """

    make_spatial(input_file_path)

def make_spatial(in_path):
    """
    Create shapefile of monthly mean bias ratios from comprehensive
    CSV file created by :mod:`calc_bias_ratios`.
    
    Arguments:
        in_path (str): path to summary_comp.CSV file containing monthly
            bias ratios, lat, long, and other data. Shapefile "summary.shp"
            is saved to parent directory of ``in_path`` under a subdirectory
            named "spatial".
            
    Returns:
        None
        
    Raises:
        FileNotFoundError: if input summary CSV file is not found.
        
    Note:
        In order to conduct spatial analysis the "comprehensive" output
        summary file needs to produced from the previous step using the
        "-c" command line flag or comp=True keyword argument to
        :mod:`calc_bias_ratios`. The summary_comp.csv file that is produced
        from that step is the input to this function.
        
    """
    if not os.path.isfile(in_path):
        raise FileNotFoundError('Input summary CSV file given'+\
                                ' was invalid or not found')
        
    in_df = pd.read_csv(in_path, index_col='STATION_ID')
    # save shapefile to "spatial" subdirectory of in_path
    path_tuple = os.path.split(in_path)
    path_root = path_tuple[0]
    file_name = path_tuple[1]

    out_dir = OPJ(path_root, 'spatial')
    out_file = OPJ(out_dir, 'summary.shp')
    # create output directory if does not exist
    if not os.path.isdir(out_dir):
        print(out_dir, ' does not exist, creating directory.')
        os.mkdir(out_dir)

    crs = from_epsg(4326) # WGS 84 projection
    # attributes of shapefile
    schema = { 
        'geometry': 'Point', 
        'properties': { 
            'Jan': 'float',
            'Feb': 'float',
            'Mar': 'float',
            'Apr': 'float',
            'May': 'float',
            'Jun': 'float',
            'Jul': 'float',
            'Aug': 'float',
            'Sep': 'float',
            'Oct': 'float',
            'Nov': 'float',
            'Dec': 'float',
            'grow_season': 'float',
            'annual': 'float',
            'STATION_ID': 'str',
            'GRIDMET_ID': 'int'

        }}
    
    # create shapefile from points, overwrite if exists
    with collection(
        out_file, 'w', 
        driver='ESRI Shapefile', 
        crs=crs, 
        schema=schema) as output:
        print(
            'Creating shapefile of bias ratios, saving to: \n',
             os.path.abspath(out_file)
        )
        for index, row in in_df.iterrows():
            point = Point(float(row.STATION_LON), float(row.STATION_LAT))
            output.write({
                'properties': {
                    'Jan': row['Jan_mean'],
                    'Feb': row['Feb_mean'],
                    'Mar': row['Mar_mean'],
                    'Apr': row['Apr_mean'],
                    'May': row['May_mean'],
                    'Jun': row['Jun_mean'],
                    'Jul': row['Jul_mean'],
                    'Aug': row['Aug_mean'],
                    'Sep': row['Sep_mean'],
                    'Oct': row['Oct_mean'],
                    'Nov': row['Nov_mean'],
                    'Dec': row['Dec_mean'],
                    'grow_season': row['April_to_oct_mean'],
                    'annual': row['Annual_mean'],
                    'STATION_ID': index,
                    'GRIDMET_ID': row['GRIDMET_ID']
                },
                'geometry': mapping(point)
            }
        )

def arg_parse():
    """
    Parse command line arguments for creating shapefile of all stations
    found in comprehensive summary CSV file.
    """
    parser = argparse.ArgumentParser(
        description='Calculate monthly bias ratios for gridMET cells.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '-i', '--input', metavar='PATH', required=True,
        help='Input CSV file of merged climate/gridMET data that '+\
             'was created by running prep_input.py, download_gridmet_ee.py'+\
             ', and calc_bias_ratios.py')
#    parser.add_argument(
#        '--debug', default=logging.INFO, const=logging.DEBUG,
#        help='Debug level logging', action="store_const", dest="loglevel")
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = arg_parse()

    main(input_file_path=args.input)
