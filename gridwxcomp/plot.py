# -*- coding: utf-8 -*-
"""
Create interactive HTML comparison plots between paired station and
gridded climatic variables or bar comparison plots between interpolated and
station point data.

"""
import logging
import os
from math import pi
from pathlib import Path
import pandas as pd
import numpy as np
from bokeh import models
from bokeh.plotting import ColumnDataSource, figure, output_file, save
from bokeh.layouts import gridplot
from refet.calcs import _wind_height_adjust
from .util import parse_yr_filter, read_config, read_data, convert_units

VAR_LIST = [
    'tmax',
    'tmin',
    'tdew',
    'rs',
    'wind',
    'ea',
    'rhmax',
    'rhmin',
    'rhavg',
    'eto',
    'etr',
    'prcp']
UNITS_DICT = {'tmax': '(C)', 'tmin': '(C)', 'tdew': '(C)',
	'rs': '(w/m2)', 'wind': '(m/s)', 'ea': '(kpa)',
	'rhmax': '(%)', 'rhmin': '(%)', 'rhavg': '(%)',
	'eto': '(mm)', 'etr': '(mm)', 'prcp': '(mm)'}

# TITLES_LIST and MONTHLY_TITLES_LIST are formatted as LaTeX
TITLES_LIST = [
    '$$Maximum\\:Temperature$$',
    '$$Minimum\\:Temperature$$',
    '$$Dewpoint\\:Temperature$$',
    '$$Solar\\:Radiation$$',
    '$$Wind\\:Speed\\:_{2m}$$',
    '$$Vapor\\:Pressure$$',
    '$$RH\\:Maximum$$',
    '$$RH\\:Minimum$$',
    '$$RH\\:Average$$',
    '$$ET_{O}$$',
    '$$ET_{r}$$',
    '$$Precipitation$$']

MONTHLY_TITLES_LIST = [
    '$$Maximum\\:Temperature\\:Monthly\\:Averages$$',
    '$$Minimum\\:Temperature\\:Monthly\\:Averages$$',
    '$$Dewpoint\\:Temperature\\:Monthly\\:Averages$$',
    '$$Solar\\:Radiation\\:Monthly\\:Averages$$',
    '$$Wind\\:Speed\\:_{2m}\\:Monthly\\:Averages$$',
    '$$Vapor\\:Pressure\\:Monthly\\:Averages$$',
    '$$RH\\:Maximum\\:Monthly\\:Averages$$',
    '$$RH\\:Minimum\\:Monthly\\:Averages$$',
    '$$RH\\:Average\\:Monthly\\:Averages$$',
    '$$ET_{O}\\:Monthly\\:Averages$$',
    '$$ET_{r}\\:Monthly\\:Averages$$',
    '$$Precipitation\\:Monthly\\:Averages$$']


# list of x (station), y (gridded) variables
x_var_list = [f'station_{var}' for var in VAR_LIST]
# list of y variables
y_var_list = [f'gridded_{var}' for var in VAR_LIST]

# timeseries y label list
ts_ylabel_list = [f'{var.upper()} {UNITS_DICT[var]}' for var in VAR_LIST]
# scatter x label, y label lists
xlabel_list = [f'Station {var.upper()} {UNITS_DICT[var]}' for var in VAR_LIST]
ylabel_list = [f'Gridded {var.upper()} {UNITS_DICT[var]}' for var in VAR_LIST]

# legend x lists
legendx_list = ['Station'] * len(TITLES_LIST)


def daily_comparison(
        input_csv,
        config_path,
        dataset_name='gridded',
        out_dir=None,
        year_filter=None):
    """
    Compare daily weather station data from
    `PyWeatherQAQC <https://github.com/WSWUP/pyWeatherQAQC>`_ with gridded data
    for each month in year specified.

    The :func:`daily_comparison` function produces HTML files with time series
    and scatter plots of station versus gridded climate variables. It uses the
    `bokeh <https://bokeh.pydata.org/en/latest/>`_ module to create interactive
    plots, e.g. they can be zoomed in/out and panned. Separate plot files are
    created for each month of a single year.

    The scatterplots for each month will allow you to visualize the overall
    correlation and relationship between the gridded and station variables.

    The timeseries plots will show all daily observations of a single month
    together, which will highlight any months that differ from their
    neighboring years, ex: a year that had a significantly colder June than
    other Junes

    Arguments:
        input_csv (str): path to input CSV file containing paired station/
            gridded metadata. This file is created by running
            :mod:`gridwxcomp.prep_metadata` followed by 
            :mod:`gridwxcomp.ee_download`.
        config_path (str): path to the config file that has the parameters used
            to interpret the station and gridded data files
        dataset_name (str): Name of gridded dataset to be used in plots
        out_dir (str or None): default None. Directory to save comparison
            plots, if None save to "daily_comp_plots" in currect directory.
        year_filter (str or None): default None. Single year YYYY or range
            YYYY-YYYY

    Returns:
        None

    Example:
        The :func:`daily_comparison` function will generate HTML files with
        bokeh plots for all paired climate variables within the config file


        or within Python,

        >>> from gridwxcomp.plot import daily_comparison
        >>> daily_comparison('merged_input.csv', 'config_file.ini', 'comp_plots_2016', '2016')

        Both methods result in monthly HTML `bokeh
        <https://bokeh.pydata.org/en/latest/>`_ plots being saved to
        "comp_plots_2016/STATION_ID/" where "STATION_ID" is the station ID as
        found in the input CSV file. A file is saved for each month with the
        station ID, month, and year in the file name.  If ``out_dir`` keyword
        argument is not given the plots will be saved to a directory named
        "daily_comp_plots".

    Note:
        If there are less than five days of data in a month the plot for that
        month will not be created.

    """

    if not out_dir:
        out_dir = os.getcwd()

    if not os.path.isdir(out_dir):
        print('{} does not exist, creating directory'.format(out_dir))
        os.makedirs(out_dir)

    year = year_filter
    logging.info('\nProcessing Year: {}'.format(year))

    # Import Station/gridded meta data shapefile
    input_df = pd.read_csv(input_csv, sep=',')
    # Import config file
    config_dict = read_config(config_path)

    # loop through each station and calculate monthly ratio
    for index, row in input_df.iterrows():
        if 'STATION_FILE_PATH' not in row or 'GRID_FILE_PATH' not in row:
            raise KeyError(
                'Missing station and/or grid file paths in ' +
                'input file. Run prep_metadata.py followed by ee_download.py '
                'first.')

        # load station and grid time series files
        try:
            # open file, convert units, then adjust wind measurement height
            raw_station_df = read_data(
                config_dict, 'station', row.STATION_FILE_PATH)
            converted_station_df = convert_units(
                config_dict, 'station', raw_station_df)
            converted_station_df['wind'] = _wind_height_adjust(
                uz=converted_station_df['wind'], 
                zw=config_dict['station_anemometer_height'])
            prepped_station_df = converted_station_df.add_prefix(
                'station_', axis='columns')
        except IOError:
            print(
                'Time series file for station: ',
                row.STATION_ID,
                ' was not found, skipping.')
            continue

        try:
            # open file, convert units, then adjust wind measurement height
            raw_gridded_df = read_data(
                config_dict, 'gridded', row.GRID_FILE_PATH)
            converted_gridded_df = convert_units(
                config_dict, 'gridded', raw_gridded_df)
            converted_gridded_df['wind'] = _wind_height_adjust(
                uz=converted_gridded_df['wind'], 
                zw=config_dict['gridded_anemometer_height'])
            prepped_gridded_df = converted_gridded_df.add_prefix(
                'gridded_', axis='columns')
        except IOError:
            print(
                'Time series file for gridded: ',
                row.STATION_ID,
                'was not found, skipping.')
            continue

        # Filter years
        if year:
            prepped_station_df, year_str = parse_yr_filter(
                prepped_station_df, year, label=row.STATION_ID)
        else:
            start_yr = int(prepped_station_df.index.year.min())
            end_yr = int(prepped_station_df.index.year.max())
            year_str = '{}_{}'.format(start_yr, end_yr)

        # Combine station and gridded dataframes
        merged = pd.concat([prepped_station_df, prepped_gridded_df], axis=1)

        # Check to see if any overlapping data exists
        if merged.shape[0] == 0:
            print(
                'No overlapping data between gridded and station source for '
                f'{row.STATION_ID}')
            continue

        for month in range(1, 13):
            logging.info('Month: {}'.format(month))
            monthly_data = merged[merged.index.month == month]

            if len(monthly_data.index) <= 5:
                logging.info('Skipping. Less than 5 observations in month.')
                continue
            # Output Folder
            out_folder = os.path.join(
                out_dir, 'daily_comp_plots', '{}'.format(
                    row.STATION_ID.replace(
                        " ", "")))

            # Create path if it doesn't exist
            if not os.path.exists(out_folder):
                os.makedirs(out_folder)

            # Output to HTML file
            out_file_path = os.path.join(out_folder, '{}_{:02}_{}.html')\
                .format(row.STATION_ID.replace(" ", ""), month, year_str)

            output_file(out_file_path)

            # empty list to append figures to
            figure_list = []
            legendy_list = [dataset_name] * len(TITLES_LIST)

            # loop through and create figures for each variable using vars
            # and plot labels from lists above
            first_plot = True
            for i, (x_var, y_var, title, ts_ylabel, xlabel, ylabel, legendx,
                    legendy) in enumerate(zip(x_var_list, y_var_list,
                                              TITLES_LIST, ts_ylabel_list,
                                              xlabel_list, ylabel_list,
                                              legendx_list, legendy_list)):

                # least squares cannot have nans (drop nas for each var
                # separately)
                monthly_data_subset = monthly_data[[x_var, y_var]]
                monthly_data_subset = monthly_data_subset.dropna()

                monthly_data_subset['date'] = monthly_data_subset.index
                monthly_data_subset.index.name = ''
                monthly_data_subset.reset_index(inplace=True)

                if monthly_data_subset.empty:
                    logging.info("Skipping {}. No Data.".format(x_var))
                    continue

                if first_plot:
                    # Initial timeseries plot to establish xrange for link axes
                    p1 = figure(width=800, height=400,
                                title=title, x_axis_type="datetime",
                                y_axis_label=ts_ylabel)
                    p1.line(monthly_data_subset.index,
                            monthly_data_subset[x_var], color="navy",
                            alpha=0.5, legend_label=legendx, line_width=2)
                    p1.line(monthly_data_subset.index,
                            monthly_data_subset[y_var], color="red",
                            alpha=0.5, legend_label=legendy, line_width=2)
                    p1.xaxis.major_label_overrides = {
                        i: date.strftime('%Y %b %d') for i, date in enumerate(
                            pd.to_datetime(monthly_data_subset.date)
                        )
                    }
                    first_plot = False

                else:
                    # Timeseries plots after first pass
                    p1 = figure(width=800, height=400,
                                title=title, x_axis_type="datetime",
                                y_axis_label=ts_ylabel,
                                x_range=p1.x_range)
                    p1.line(
                        monthly_data_subset.index,
                        monthly_data_subset[x_var],
                        color="navy",
                        alpha=0.5,
                        legend_label=legendx,
                        line_width=2)
                    p1.line(monthly_data_subset.index,
                            monthly_data_subset[y_var], color="red", alpha=0.5,
                            legend_label=legendy, line_width=2)

                p1.xaxis.major_label_overrides = {
                    i: date.strftime('%Y %b %d') for i, date in enumerate(
                        pd.to_datetime(monthly_data_subset.date)
                    )
                }

                # 1 to 1 Plot
                # Regression through Zero
                # https://stackoverflow.com/questions/9990789/how-to-force-
                #   zero-interception-in-linear-regression/9994484#9994484
                m = np.linalg.lstsq(
                    monthly_data_subset[x_var].values.reshape(
                        -1, 1), monthly_data_subset[y_var], rcond=None)[0][0]
                r_x, r_y = zip(*((i, i * m) for i in range(
                    int(np.min([monthly_data_subset[y_var], 
                        monthly_data_subset[x_var]]) - 2),
                    int(np.max([monthly_data_subset[y_var], 
                        monthly_data_subset[x_var]]) + 3), 1)))
                # Plots
                p2 = figure(width=400, height=400,
                            x_axis_label=xlabel, y_axis_label=ylabel,
                            title='Slope Through Zero: m = {}'.format(
                                round(m, 4)))
                p2.scatter(
                    monthly_data_subset[x_var],
                    monthly_data_subset[y_var],
                    size=15,
                    color="navy",
                    alpha=0.5)
                p2.line(
                    [int(np.min([monthly_data_subset[y_var], 
                        monthly_data_subset[x_var]]) - 2),
                     int(np.max([monthly_data_subset[y_var], 
                         monthly_data_subset[x_var]]) + 2)],
                    [int(np.min([monthly_data_subset[y_var], 
                        monthly_data_subset[x_var]]) - 2),
                     int(np.max([monthly_data_subset[y_var], 
                         monthly_data_subset[x_var]]) + 2)],
                    color="black", legend_label='1 to 1 line')

                p2.line(r_x, r_y, color="red", legend_label='Reg thru zero')
                p2.legend.location = "top_left"

                # Append [p1, p2] to figure_list (create list of lists)
                figure_list.append([p1, p2])

            # Plot all figures in list
            fig = gridplot(figure_list, toolbar_location="left")
            # Save the figure
            save(fig)


def monthly_comparison(
        input_csv,
        config_path,
        dataset_name='gridded',
        out_dir=None,
        day_limit=10):
    """
    Compare monthly average weather station data from `PyWeatherQAQC
    <https://github.com/WSWUP/pyWeatherQAQC>`_ with the gridded dataset.

    The :func:`monthly_comparison` function produces HTML files with time series
    and scatter plots of station versus gridded climate variables of monthly
    mean data. It uses the `bokeh <https://bokeh.pydata.org/en/latest/>`_ module
    to create interactive plots, e.g. they can be zoomed in/out and panned.

    Arguments:
        input_csv (str): path to input CSV file containing paired 
            station:gridded metadata. This file is created by running 
            :mod:`gridwxcomp.prep_metadata` followed by
            :mod:`gridwxcomp.ee_download`.
        config_path (str): path to the config file that has the parameters used
            to interpret the station and gridded data files'
        dataset_name (str): Name of gridded dataset to be used in plots
        out_dir (str): default None. Directory to save comparison plots.
        day_limit (int): default 10. Number of paired days per month that must
            exist for variable to be plotted.

    Returns:
        None

    Example:
        The :func:`monthly_comparison` function will generate HTML files with
        bokeh plots for paired climate variable, e.g. etr_mm, tmax_c

        >>> from gridwxcomp.plot import monthly_comparison
        >>> monthly_comparison('merged_input.csv', 'monthly_plots')

        Both methods result in monthly HTML bokeh plots being saved to
        "monthly_plots/" which contains a plot file for each station as found
        in the input CSV file. If ``out_dir`` keyword argument is not given the
        plots will be saved to a directory named "monthly_comp_plots".

    Note:
        If there are less than 2 months of data the plot for that station 
        will not be created.
    """
    if not out_dir:
        out_dir = os.getcwd()

    if not os.path.isdir(out_dir):
        print('{} does not exist, creating directory'.format(out_dir))
        os.makedirs(out_dir)

    # Import Station/gridded meta data shapefile
    input_df = pd.read_csv(input_csv, sep=',')
    # Import config file
    config_dict = read_config(config_path)

    # loop through each station and calculate monthly ratio
    for index, row in input_df.iterrows():
        if 'STATION_FILE_PATH' not in row or 'GRID_FILE_PATH' not in row:
            raise KeyError(
                'Missing station and/or grid file paths in ' +
                'input file. Run prep_metadata followed by download_grid_data '
                'first.')

        # load station and grid time series files
        try:
            # open file, convert units, then adjust wind measurement height
            raw_station_df = read_data(
                config_dict, 'station', row.STATION_FILE_PATH)
            converted_station_df = convert_units(
                config_dict, 'station', raw_station_df)
            converted_station_df['wind'] = _wind_height_adjust(
                uz=converted_station_df['wind'], 
                zw=config_dict['station_anemometer_height'])
            prepped_station_df = converted_station_df.add_prefix(
                'station_', axis='columns')
        except IOError:
            print(
                'Time series file for station: ',
                row.STATION_ID,
                ' was not found, skipping.')
            continue

        try:
            # open file, convert units, then adjust wind measurement height
            raw_gridded_df = read_data(
                config_dict, 'gridded', row.GRID_FILE_PATH)
            converted_gridded_df = convert_units(
                config_dict, 'gridded', raw_gridded_df)
            converted_gridded_df['wind'] = _wind_height_adjust(
                uz=converted_gridded_df['wind'], 
                zw=config_dict['gridded_anemometer_height'])
            prepped_gridded_df = converted_gridded_df.add_prefix(
                'gridded_', axis='columns')
        except IOError:
            print(
                'Time series file for gridded: ',
                row.STATION_ID,
                'was not found, skipping.')
            continue

        # Combine station and gridded dataframes and crop to just the shared
        # data
        start_date = max(
            prepped_station_df.index.date.min(),
            prepped_gridded_df.index.date.min())
        end_date = min(
            prepped_station_df.index.date.max(),
            prepped_gridded_df.index.date.max())
        merged = pd.concat([prepped_station_df, prepped_gridded_df], axis=1)
        merged = merged.loc[start_date:end_date]

        # Check to see if any overlapping data exists
        if merged.shape[0] == 0:
            print(
                'No overlapping data between gridded and station source '
                f'for {row.STATION_ID}')
            continue

        # remove all pairs where one var missing
        for (x_var, y_var) in zip(x_var_list, y_var_list):
            merged[[x_var, y_var]] = merged[[x_var, y_var]].dropna()

        # Monthly averages including count
        monthly = merged.groupby([lambda x: x.year, lambda x: x.month]).agg(
            ['mean', 'sum', 'count'])

        # Remove months with Less Than XX Days in average
        var_names = list(monthly.columns.levels)[0]
        for v in var_names:
            mask = monthly.loc[:, (v, 'count')] < day_limit
            monthly.loc[mask, ('sum', 'mean')] = np.nan

        # Rebuild Index DateTime
        monthly['year'] = monthly.index.get_level_values(0).values
        monthly['month'] = monthly.index.get_level_values(1).values
        monthly.index = pd.to_datetime(
            monthly.year * 10000 + monthly.month * 100 + 15,
            format='%Y%m%d')

        if len(monthly.index) < 2:
            logging.info('Skipping. Less than 2 months of observations.')
            continue

        # Output Folder
        out_folder = os.path.join(out_dir, 'monthly_comp_plots')

        # Create path if it doesn't exist
        if not os.path.exists(out_folder):
            os.makedirs(out_folder)

        # Output to HTML file
        out_file_path = os.path.join(out_folder, '{}.html')\
            .format(row.STATION_ID.replace(" ", ""))
        output_file(out_file_path)

        # empty list to append figures to
        figure_list = []
        legendy_list = [dataset_name] * len(TITLES_LIST)

        # loop through and create figures for each variable using vars
        # and plot labels from lists above
        first_plot = True
        for i, (x_var, y_var, title, ts_ylabel, xlabel, 
                ylabel, legendx, legendy) in enumerate(zip(
            x_var_list, y_var_list, MONTHLY_TITLES_LIST, ts_ylabel_list,
                xlabel_list, ylabel_list, legendx_list, legendy_list)):

            # todo put an if statement here if precip sums are needed
            stat = 'mean'
            # lstsq cannot have nans (drop nas for each var separately)
            monthly2 = monthly[[x_var, y_var]]
            monthly2 = monthly2.dropna()

            if monthly2.empty:
                logging.info("Skipping {}. No Data.".format(x_var))
                continue

            if first_plot:
                # Initial timeseries plot to establish xrange for link axes
                p1 = figure(width=800, height=400,
                            x_axis_type="datetime", title=title,
                            y_axis_label=ts_ylabel)
                p1.line(monthly.index.to_pydatetime(),
                        monthly[x_var, stat], color="navy",
                        alpha=0.5, legend_label=legendx, line_width=2)
                p1.line(monthly.index.to_pydatetime(),
                        monthly[y_var, stat], color="red",
                        alpha=0.5, legend_label=legendy, line_width=2)

                first_plot = False
            else:
                # Timeseries plots after first pass
                p1 = figure(width=800, height=400,
                            x_axis_type="datetime", title=title,
                            y_axis_label=ts_ylabel,
                            x_range=p1.x_range)
                p1.line(monthly.index.to_pydatetime(),
                        monthly[x_var, stat], color="navy", alpha=0.5,
                        legend_label=legendx, line_width=2)
                p1.line(monthly.index.to_pydatetime(),
                        monthly[y_var, stat], color="red", alpha=0.5,
                        legend_label=legendy, line_width=2)

            # 1 to 1 Plot
            # Regression through Zero
            # https://stackoverflow.com/questions/9990789/how-to-force-
            #   zero-interception-in-linear-regression/9994484#9994484
            m = np.linalg.lstsq(monthly2[x_var, stat].values.reshape(-1, 1),
                                monthly2[y_var, stat], rcond=None)[0][0]
            r_x, r_y = zip(*((i, i * m) for i in range(
                int(np.min([monthly2[y_var, stat],
                            monthly2[x_var, stat]]) - 2),
                int(np.max([monthly2[y_var, stat],
                            monthly2[x_var, stat]]) + 3), 1)))
            # Plots
            p2 = figure(width=400, height=400,
                        x_axis_label=xlabel, y_axis_label=ylabel,
                        title='Slope Through Zero: m = {}'.format(
                            round(m, 4)))
            p2.scatter(monthly2[x_var, stat], monthly2[y_var, stat],
                       size=15, color="navy", alpha=0.5)
            p2.line(
                [int(np.min([monthly2[y_var, stat], 
                    monthly2[x_var, stat]]) - 2),
                 int(np.max([monthly2[y_var, stat], 
                     monthly2[x_var, stat]]) + 2)],
                [int(np.min([monthly2[y_var, stat], 
                    monthly2[x_var, stat]]) - 2),
                 int(np.max([monthly2[y_var, stat], 
                     monthly2[x_var, stat]]) + 2)],
                color="black", legend_label='1 to 1 line')
            p2.line(r_x, r_y, color="red", legend_label='Reg thru zero')
            p2.legend.location = "top_left"

            # Append [p1, p2] to figure_list (create list of lists)
            figure_list.append([p1, p2])

        # Plot all figures in list
        fig = gridplot(figure_list, toolbar_location="left")
        
        # Save the figure
        save(fig)


def station_bar_plot(
        summary_csv,
        bar_plot_layer,
        out_dir=None,
        y_label=None,
        title=None,
        subtitle=None,
        year_subtitle=True):
    """
    Produce an interactive bar chart comparing multiple climate stations to 
    each other for a particular variable, e.g. bias ratios or interpolated 
    residuals.

    This function may also be used for any numerical data in the summary CSV
    files that are created by :func:`gridwxcomp.interpolate` in addition to
    those created by :func:`gridwxcomp.calc_bias_ratios`.  The main requirement
    is that ``summary_csv`` must contain the column 'STATION_ID' and the
    ``bar_plot_layer`` keyword argument.

    Arguments:
        summary_csv (str, Path): path to summary CSV produced by either 
            :func:`gridwxcomp.calc_bias_ratios` or by 
            :func:`gridwxcomp.interpolate`. Should contain ``bar_plot_layer`` 
            data for plot.
        bar_plot_layer (str): name of variable to plot.
        out_dir (str or None): default None. Output directory path, default is
            'station_bar_plots' in parent directory of ``summary_csv``.
        y_label (str or None): default None. Label for y-axis, defaults 
            to ``bar_plot_layer``.
        title (str or None): default None. Title of plot.
        subtitle (str, list, or None): default None. Additional subtitle(s) 
            for plot.
        year_subtitle (bool): default True. If true print subtitle on plot with
            the max year range used for station data, e.g. 'years: 1995-2005'

    Example:
        Let's say we want to compare the mean growing seasion bias ratios of
        reference evapotranspiration (ETr) for the selection of stations we
        used to calculate bias ratios.

        The summary CSV file containing the ratios should be first
        created using :func:`gridwxcomp.calc_bias_ratios`.

            >>> from gridwxcomp.plot import station_bar_plot
            >>> # path to summary CSV with station data
            >>> in_file = 'monthly_ratios/etr_mm_summary_all_yrs.csv'
            >>> example_layer = 'grow_mean'
            >>> station_bar_plot(in_file, example_layer)
        
        The resulting file will be saved using the bar_plot_layer name as a
        file name::
        
        'monthly_ratios/station_bar_plots/grow_mean.html'
        
        The plot file will contain the mean growing season bias ratios
        of ETr for each station, sorted from smallest to largest values.

    Raises:
        FileNotFoundError: if ``summary_csv`` is not found.
        KeyError: if ``bar_plot_layer`` does not exist as a column name in ``summary_csv``.
    """

    if not Path(summary_csv).is_file():
        err_msg = '\n{} is not a valid path to a summary CSV file!'.\
            format(summary_csv)
        raise FileNotFoundError(err_msg)

    df = pd.read_csv(summary_csv, na_values=[-999])

    if bar_plot_layer not in df.columns:
        err_msg = '\nColumn {} was not found in {}'.format(
            bar_plot_layer, summary_csv)
        raise KeyError(err_msg)

    df.sort_values(bar_plot_layer, inplace=True)
    df.index.name = 'dummy_name'  # fix internal to bokeh- reset_index
    source = ColumnDataSource(df)
    # hover tooltip with station and value
    tooltips = [
        ("station", "@STATION_ID"),
        ("value", "@{}".format(bar_plot_layer)),
    ]
    hover = models.HoverTool(tooltips=tooltips)

    if not y_label:
        y_label = bar_plot_layer
    # save to working directory in 'station_bar_plots' if not specified
    if not out_dir:
        out_dir = Path(summary_csv).parent / 'station_bar_plots'
    else:
        out_dir = Path(out_dir)
    if not out_dir.is_dir():
        print('\n{}\nDoes not exist, making directory'.format(
            out_dir.absolute()))
        out_dir.mkdir(parents=True, exist_ok=True)

    out_file = out_dir / '{}.html'.format(bar_plot_layer)
    print('\nCreating station bar plot for variable: ', bar_plot_layer,
          '\nUsing data from file: ', Path(summary_csv).absolute())

    output_file(out_file)

    p = figure(x_range=df.STATION_ID, y_axis_label=y_label, title=title)
    p.vbar(x='STATION_ID', top=bar_plot_layer, width=0.8, source=source)
    p.xaxis.major_label_orientation = pi / 2
    p.add_tools(hover, models.BoxSelectTool())

    if year_subtitle:
        # add data range (years start to end) as subtitle
        min_yr = int(df.start_year.min())
        max_yr = int(df.end_year.max())
        if min_yr == max_yr:
            year_str = 'year: {}'.format(min_yr)
        else:
            year_str = 'years: {}-{}'.format(min_yr, max_yr)
        # caution note if not all stations use full year range
        if not (
                df.end_year == max_yr).all() or not (
                df.start_year == min_yr).all():
            year_str = '{} (less years exist for some stations)'.format(
                year_str)

        p.add_layout(
            models.Title(
                text=year_str,
                text_font_style="italic"),
            'above')

    # add arbitrary number of custom subtitles as lines above plot
    if isinstance(subtitle, (list, tuple)):
        for st in subtitle:
            p.add_layout(
                models.Title(
                    text=st,
                    text_font_style="italic"),
                'above')
    elif subtitle:
        p.add_layout(
            models.Title(
                text=subtitle,
                text_font_style="italic"),
            'above')

    save(p)
    print('\nPlot saved to: ', out_file.absolute())
