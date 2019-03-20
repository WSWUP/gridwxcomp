# -*- coding: utf-8 -*-
"""
Utility functions or classes for ``gridwxcomp`` package
"""

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
    if years == 'all':
        dt_df = dt_df
        year_str = 'all_yrs'
    else:
        try:
            if years and '-' in years:
                err_msg = '{} is not a valid years option,\n'.format(years) +\
                    'use single or range e.g. 2015 or 2000-2010'
                start, end = years.split('-')
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
                year_str = str(years)
                if not years in dt_df.index:
                    print('WARNING: station:', label, 'is missing data',
                        'for year:', years)
                else:
                    
                    dt_df = dt_df.loc[years]
        except:
            raise ValueError(err_msg)

    ret = dt_df, year_str
    return ret



