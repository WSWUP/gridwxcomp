# -*- coding: utf-8 -*-
"""
Utility functions or classes for ``gridwxcomp`` package
"""

import pkg_resources
from pathlib import Path

def get_gridmet_meta_csv(gridmet_meta_path=None):
    """Find path to 'gridmet_cell_data.csv' packaged with gridwxcomp"""
    # look for pacakged gridmet_cell_data.csv if path not given
    if not gridmet_meta_path:
        try:
            if pkg_resources.resource_exists('gridwxcomp',
                    'gridmet_cell_data.csv'):
                gridmet_meta_path = pkg_resources.resource_filename(
                    'gridwxcomp',
                    'gridmet_cell_data.csv'
                    )
        except:
            gridmet_meta_path = 'gridmet_cell_data.csv'
    if not Path(gridmet_meta_path).is_file():
        raise FileNotFoundError('GridMET file path was not given and '+\
                'gridmet_cell_data.csv was not found in the gridwxcomp '+\
                'install directory. Please assign the path or put '+\
                '"gridmet_cell_data.csv" in the current working directory.\n')
    return gridmet_meta_path


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
    err_msg = ('{} is not a valid years option,\n'.format(years),
                    'use single or range e.g. 2015 or 2000-2010')
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
                if not start in dt_df.index:
                    data_start = dt_df.index.year.min()
                    print('WARNING: data for {l} starts in {d}'\
                              .format(l=label, d=data_start) +\
                         ' but you gave {s}'.format(s=start))
                if not end in dt_df.index:
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



