# -*- coding: utf-8 -*-
"""
Create interactive HTML comparison plots between paired station and 
gridMET climatic variables or bar comparison plots between interpolated and 
station point data. 

"""
import argparse
import logging
import os
import sys
import datetime as dt
from math import pi
from pathlib import Path

import pandas as pd
import numpy as np
from bokeh import models
from bokeh.plotting import ColumnDataSource, figure, output_file, save
from bokeh.layouts import gridplot

# allows for CL script usage if gridwxcomp not installed
try:
    from .util import parse_yr_filter
except:
    from util import parse_yr_filter

def daily_comparison(input_csv, out_dir=None, year_filter=None):
    """
    Compare daily weather station data from 
    `PyWeatherQAQC <https://github.com/WSWUP/pyWeatherQAQC>`_ with gridMET 
    for each month in year specified.

    The :func:`daily_comparison` function produces HTML files with time series 
    and scatter plots of station versus gridMET climate variables. It uses the 
    `bokeh <https://bokeh.pydata.org/en/latest/>`_ module to create interactive 
    plots, e.g. they can be zoomed in/out and panned. Separate plot files are 
    created for each month of a single year. 

    Arguments:
        input_csv (str): path to input CSV file containing paired station/
            gridMET metadata. This file is created by running 
            :mod:`gridwxcomp.prep_input` followed by :mod:`gridwxcomp.download_gridmet_ee`.

    Keyword Arguments:
        out_dir (str or None): default None. Directory to save comparison 
            plots, if None save to "daily_comp_plots" in currect directory. 
        year_filter (str or None): default None. Single year YYYY or range 
            YYYY-YYYY

    Returns:
        None

    Example:
        The :func:`daily_comparison` function will generate HTML files with 
        bokeh plots for paired climate variables, e.g. etr_mm, eto_mm, 
        u2_ms, tmin_c, tmax_c, srad_wm2, ea_kpa, and Ko (dew point depression). 
        Monthly plots are created for a single year.
        
        From the command line, use the "plot" command with the 
        ``[-t, --plot-type]`` option set to station-grid-comp and 
        the ``[-f, --freq]`` option left as default ("daily"),

        .. code-block:: sh

            $ gridwxcomp plot merged_input.csv -t station-grid-comp -o comp_plots_2016 -y 2016

        or within Python,

        >>> from gridwxcomp.plot import daily_comparison
        >>> daily_comparison('merged_input.csv', 'comp_plots_2016', '2016')

        Both methods result in monthly HTML `bokeh <https://bokeh.pydata.org/en/latest/>`_ 
        plots being saved to "comp_plots_2016/STATION_ID/" where "STATION_ID" 
        is the station ID as found in the input CSV file. A file is saved for 
        each month with the station ID, month, and year in the file name. 
        If ``out_dir`` keyword argument or ``[-o, --out-dir]`` command line 
        option is not given the plots will be saved to a directory named 
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

    # # Import Station/GRIDMET meta data shapefile
    paired_data = pd.read_csv(input_csv, sep=',')

    # List of variables to compare (STATION/gridMET ORDER SHOULD MATCH)
    station_vars = ['TMin (C)', 'TMax (C)', 'wx_Ko_c', 'Rs (w/m2)',
                    'ws_2m (m/s)', 'Vapor Pres (kPa)', 'RHAvg (%)',
                    'Precip (mm)', 'Calc_ETo (mm)', 'Calc_ETr (mm)']

    gridmet_vars = ['tmin_c', 'tmax_c', 'grid_Ko_c', 'srad_wm2', 'u2_ms',
                    'ea_kpa', 'rh_avg', 'prcp_mm', 'eto_mm', 'etr_mm']

    # # Limit row processing range (testing)
    # start = 0
    # end = 1
    #Loop through each station/gridmet pair
    for index, row in paired_data.iterrows():
    # #    Limit iteration during development
    #     if index < start:
    #         continue
    #     if index >= end:
    #       break

        # clear previous datasets
        grid_data = []
        station_data = []

        station_path = row.STATION_FILE_PATH
        logging.info('\nStation: {}'.format(row.STATION_ID))

        # Check is station path is given in input
        if pd.isnull(station_path):
            logging.info('Station path is not given. Skipping.')
            continue

        # Skip If FILE DOES NOT EXIST
        if not os.path.exists(station_path):
            logging.info('SKIPPING {}. NO STATION FILE FOUND.'.format(
                station_path))
            continue
        else:
            station_data = pd.read_excel(station_path,
                                         sheet_name='Corrected Data')
            # Filter years
            if year:
                station_data, year_str = parse_yr_filter(
                    station_data, year, label=row.STATION_ID)
            else:
                start_yr = station_data.year.min().astype(int) 
                end_yr = station_data.year.max().astype(int) 
                year_str = '{}_{}'.format(start_yr, end_yr)

        # Import GRIDMET Data
        grid_path = row.GRIDMET_FILE_PATH
        # Skip if GRIDMET FILE DOES NOT EXIST
        if not os.path.exists(grid_path):
            print('SKIPPING {}. NO GRIDMET FILE FOUND.'.format(grid_path))
            continue
        else:
            grid_data = pd.read_csv(grid_path, sep=',',parse_dates=True,
                                    index_col='date')
            # Filter to specific year
            # grid_data = grid_data[grid_data['year'] == year]

            # Add Tdew to gridmet dataset Teten's equation ASCE REF-ET
            #  supporting equations Appendix 2-1

            grid_data['tdew_c'] = (116.91 + 237.3 * np.log(grid_data.ea_kpa)) /\
                                  (16.78 - np.log(grid_data.ea_kpa))

            # Calculate Tmin - Tdew = Ko for both Station and GridMET
            # Dew Point Depression
            grid_data['grid_Ko_c'] = grid_data.tmin_c - grid_data.tdew_c

            station_data['wx_Ko_c'] = station_data['TMin (C)'] - \
                                      station_data['TDew (C)']

            # grid RH Avg calc
            # Saturated Vapor Pressure
            grid_data['tavg_c'] = (grid_data.tmin_c + grid_data.tmax_c) / 2
            grid_data['e_sat_kpa'] = 0.6108 * np.exp(
                (17.27 * grid_data.tavg_c) /
                (grid_data.tavg_c + 237.3))
            # Average RH (%)
            grid_data['rh_avg'] = (grid_data.ea_kpa / grid_data.e_sat_kpa) * 100

            # Combine station and gridMET dataframes (only plotting variables)
            merged = pd.concat([station_data[station_vars],
                                grid_data[gridmet_vars]], axis=1,
                               join_axes=[station_data.index])
            # Remove results with na
            # merged = merged.dropna()

            for month in range(1,13):
                logging.info('Month: {}'.format(month))
                monthly_data = merged[merged.index.month==month]

                if len(monthly_data.index)<= 5:
                     logging.info('Skipping. Less than 5 observations in '
                                  'month.')
                     continue
                # Output Folder
                out_folder =  os.path.join(out_dir, 'daily_comp_plots',
                                           '{}'.format(
                                               row.STATION_ID.replace(" ","")))

                # Create path if it doesn't exist
                if not os.path.exists(out_folder):
                    os.makedirs(out_folder)

                # Output to HTML file
                out_file_path = os.path.join(out_folder, '{}_{:02}_{}.html')\
                    .format(row.STATION_ID.replace(" ", ""), month, year_str)

                output_file(out_file_path)

                station_vars = ['TMin (C)', 'TMax (C)', 'wx_Ko_c', 'Rs (w/m2)',
                                'ws_2m (m/s)', 'Vapor Pres (kPa)', 'RHAvg (%)',
                                'Precip (mm)', 'Calc_ETo (mm)', 'Calc_ETr (mm)']

                gridmet_vars = ['tmin_c', 'tmax_c', 'grid_Ko_c', 'srad_wm2',
                                'u2_ms',
                                'ea_kpa', 'rh_avg', 'prcp_mm', 'eto_mm',
                                'etr_mm']

                # list of x variables
                x_var_list= station_vars
                # list of y variables
                y_var_list= gridmet_vars
                # title list
                title_list= ['TMin', 'TMax', 'Ko' , 'Rs', 'WS 2m',
                               'ea', 'RH', 'Prcp', 'ETo', 'ETr']
                # timeseries y label list
                ts_ylabel_list = ['TMin (C)', 'TMax (C)', 'Ko (C)', 'Rs (w/m2)',
                                  'WS 2m (m/s)', 'ea (kPa)', 'Avg RH (%)',
                                  'Prcp (mm)',
                                  'ETo (mm)', 'ETr (mm)']
                # scatter xlabel list
                xlabel_list = ['Station TMin (C)', 'Station TMax (C)',
                               'Station Ko (C)', 'Station Rs (w/m2)',
                               'Station WS 2m (m/s)', 'Station ea (kPa)',
                               'Station RH (%)', 'Station Prcp (mm)',
                               'Station ETo (mm)', 'Station ETr (mm)']
                # scatter ylabel list
                ylabel_list = ['gridMET TMin (C)', 'gridMET TMax (C)',
                               'gridMET Ko (C)', 'gridMET Rs (w/m2)',
                               'gridMET WS 2m (m/s)', 'gridMET ea (kPa)',
                               'gridMET RH (%)', 'gridMET Prcp (mm)',
                               'gridMET ETo (mm)', 'gridMET ETr (mm)']
                # legendx list
                legendx_list = ['Station'] * len(title_list)
                # legend y list
                legendy_list = ['gridMET'] * len(title_list)

                # empty list to append figures to
                figure_list = []

                # loop through and create figures for each variable using vars
                # and plot labels from lists above
                for i, (x_var, y_var, title, ts_ylabel, xlabel, ylabel, legendx,
                        legendy) in enumerate(zip(x_var_list, y_var_list,
                                                  title_list, ts_ylabel_list,
                                                  xlabel_list, ylabel_list,
                                                  legendx_list, legendy_list)):

                    # lstsq cannot have nans (drop nas for each var separately)
                    monthly_data2 = monthly_data[[x_var, y_var]]
                    monthly_data2 = monthly_data2.dropna()

                    if monthly_data2.empty:
                        logging.info("Skipping {}. No Data.".format(x_var))
                        continue

                    if i == 0:
                        # Initial timeseries plot to establish xrange for link axes
                        p1 = figure(plot_width=800, plot_height=400,
                                    x_axis_type="datetime",title = title,
                                    y_axis_label = ts_ylabel)
                        p1.line(monthly_data2.index.to_pydatetime(),
                                monthly_data2[x_var],  color="navy",
                                alpha=0.5, legend=legendx,line_width=2)
                        p1.line(monthly_data2.index.to_pydatetime(),
                                monthly_data2[y_var],  color="red",
                                alpha=0.5, legend=legendy,line_width=2)
                    else:
                        # Timeseries plots after first pass
                        p1 = figure(plot_width=800, plot_height=400,
                                    x_axis_type="datetime",title = title,
                                    y_axis_label = ts_ylabel,
                                    x_range=p1.x_range)
                        p1.line(monthly_data2.index.to_pydatetime(),
                                monthly_data2[x_var],  color="navy", alpha=0.5,
                                legend=legendx,line_width=2)
                        p1.line(monthly_data2.index.to_pydatetime(),
                                monthly_data2[y_var],  color="red", alpha=0.5,
                                legend=legendy,line_width=2)

                    # 1 to 1 Plot
                    # Regression through Zero
                    # https://stackoverflow.com/questions/9990789/how-to-force-
                    # zero-interception-in-linear-regression/9994484#9994484

                    m = np.linalg.lstsq(monthly_data2[x_var].values.reshape(-1,1),
                                        monthly_data2[y_var], rcond=None)[0][0]
                    r_x, r_y = zip(*((i, i*m ) for i in range(
                        int(np.min([monthly_data2[y_var],monthly_data2[x_var]])-2),
                                     int(np.max([monthly_data2[y_var],
                                                 monthly_data2[x_var]])+3),1)))
                    # Plots
                    p2 = figure(plot_width=400, plot_height=400,
                                x_axis_label = xlabel, y_axis_label = ylabel,
                                title = 'Slope Through Zero: m = {}'.format(
                                    round(m,4)))
                    p2.circle(monthly_data2[x_var], monthly_data2[y_var],
                              size=15, color="navy", alpha=0.5)
                    p2.line([int(np.min([monthly_data2[y_var],
                                         monthly_data2[x_var]])-2),int(np.max(
                        [monthly_data2[y_var],monthly_data2[x_var]])+2)],
                             [int(np.min([monthly_data2[y_var],
                                          monthly_data2[x_var]])-2),int(np.max(
                                 [monthly_data2[y_var],monthly_data2[x_var]])+2)],
                              color = "black", legend = '1 to 1 line')
                    p2.line(r_x, r_y, color="red", legend = 'Reg thru zero')
                    p2.legend.location = "top_left"

                    # Append [p1, p2] to figure_list (create list of lists)
                    figure_list.append([p1, p2])

                # Plot all figures in list
                fig = gridplot(figure_list, toolbar_location="left")

                # Save the figure
                save(fig)

def monthly_comparison(input_csv, out_dir=None):
    """
    Compare monthly average weather station data from
    `PyWeatherQAQC <https://github.com/WSWUP/pyWeatherQAQC>`_ with gridMET.

    The :func:`monthly_comparison` function produces HTML files with time series
    and scatter plots of station versus gridMET climate variables of monthly 
    mean data. It uses the `bokeh <https://bokeh.pydata.org/en/latest/>`_ module
    to create interactive plots, e.g. they can be zoomed in/out and panned. 

    Arguments:
        input_csv (str): path to input CSV file containing
            paired station/gridMET metadata. This file is
            created by running :mod:`gridwxcomp.prep_input` followed by
            :mod:`gridwxcomp.download_gridmet_ee`.

    Keyword Arguments:
        out_dir (str): default None. Directory to save comparison plots.

    Returns:
        None

    Example:
        The :func:`monthly_comparison` function will generate HTML files with
        bokeh plots for paired climate variable, e.g. etr_mm,
        eto_mm, u2_ms, tmin_c, tmax_c, srad_wm2, ea_kpa, and Ko (dew point
        depression).
        
        From the command line, use the "plot" command with the 
        ``[-t, --plot-type]`` option set to station-grid-comp and
        the ``[-f, --freq]`` option set to "monthly",

        .. code-block:: sh

            $ gridwxcomp plot merged_input.csv -t station-grid-comp -freq monthly -o monthly_plots

        or within Python,

        >>> from gridwxcomp.plot import monthly_comparison
        >>> monthly_comparison('merged_input.csv', 'monthly_plots')

        Both methods result in monthly HTML bokeh plots being saved
        to "monthly_plots/" which contains a plot file for each station
        as found in the input CSV file. If ``out_dir`` keyword argument or
        ``[-o, --out-dir]`` command line option is not given the plots will
        be saved to a directory named "monthly_comp_plots".

    Note:
        If there are less than 2 months of data the plot for that
        station will not be created.

    """
    if not out_dir:
        out_dir = os.getcwd()

    if not os.path.isdir(out_dir):
        print('{} does not exist, creating directory'.format(out_dir))
        os.makedirs(out_dir)

    # # Import Station/GRIDMET meta data shapefile
    paired_data = pd.read_csv(input_csv, sep=',')

    # List of variables to compare (STATION/gridMET ORDER SHOULD MATCH)
    station_vars = ['TMin (C)', 'TMax (C)', 'wx_Ko_c', 'Rs (w/m2)',
                    'ws_2m (m/s)', 'Vapor Pres (kPa)', 'RHAvg (%)',
                    'Precip (mm)', 'Calc_ETo (mm)', 'Calc_ETr (mm)']

    gridmet_vars = ['tmin_c', 'tmax_c', 'grid_Ko_c', 'srad_wm2', 'u2_ms',
                    'ea_kpa', 'rh_avg', 'prcp_mm', 'eto_mm', 'etr_mm']

    # # Limit row processing range (testing)
    # start = 0
    # end = 1
    #Loop through each station/gridmet pair
    for index, row in paired_data.iterrows():
    # #    Limit iteration during development
    #     if index < start:
    #         continue
    #     if index >= end:
    #       break

        # clear previous datasets
        grid_data = []
        station_data = []


        logging.info('\nStation: {}'.format(row.STATION_ID))
        station_path = row.STATION_FILE_PATH

        # Check is station path is given in input
        if pd.isnull(station_path):
            logging.info('Station path is not given. Skipping.')
            continue

        # Skip If FILE DOES NOT EXIST
        if not os.path.exists(station_path):
            logging.info('SKIPPING {}. NO STATION FILE FOUND.'.format(
                station_path))
            continue
        else:
            station_data = pd.read_excel(station_path,
                                         sheet_name='Corrected Data')

        # Import GRIDMET Data
        grid_path = row.GRIDMET_FILE_PATH
        # Skip if GRIDMET FILE DOES NOT EXIST
        if not os.path.exists(grid_path):
            print('SKIPPING {}. NO GRIDMET FILE FOUND.'.format(grid_path))
            continue
        else:
            grid_data = pd.read_csv(grid_path, sep=',',parse_dates=True,
                                    index_col='date')
            # Filter to specific year
            # grid_data = grid_data[grid_data['year'] == year]

            # Add Tdew to gridmet dataset Teten's equation ASCE REF-ET
            #  supporting equations Appendix 2-1

            grid_data['tdew_c'] = (116.91 + 237.3 * np.log(grid_data.ea_kpa)) /\
                                  (16.78 - np.log(grid_data.ea_kpa))

            # Calculate Tmin - Tdew = Ko for both Station and GridMET
            # Dew Point Depression
            grid_data['grid_Ko_c'] = grid_data.tmin_c - grid_data.tdew_c

            station_data['wx_Ko_c'] = station_data['TMin (C)'] - \
                                      station_data['TDew (C)']

            # grid RH Avg calc
            # Saturated Vapor Pressure
            grid_data['tavg_c'] = (grid_data.tmin_c + grid_data.tmax_c )/2
            grid_data['e_sat_kpa'] = 0.6108*np.exp((17.27*grid_data.tavg_c)/
                                                   (grid_data.tavg_c+237.3))
            # Average RH (%)
            grid_data['rh_avg'] = (grid_data.ea_kpa/grid_data.e_sat_kpa)*100

            # Combine station and gridMET dataframes (only plotting variables)
            merged = pd.concat([station_data[station_vars],
                                grid_data[gridmet_vars]], axis=1,
                               join_axes=[station_data.index])

            # Remove results with na
            # merged = merged.dropna()


            # Monthly averages including count
            monthly = merged.groupby([lambda x: x.year, lambda x: x.month]).agg(
                ['mean', 'sum' ,'count'])

            # Remove months with Less Than XX Days in average
            day_limit = 10
            monthly = monthly[monthly['Calc_ETr (mm)', 'count'] >= day_limit]

            # Rebuild Index DateTime
            monthly['year'] = monthly.index.get_level_values(0).values
            monthly['month'] = monthly.index.get_level_values(1).values
            monthly.index = pd.to_datetime(
                monthly.year * 10000 + monthly.month * 100 + 15,
                format='%Y%m%d')

            if len(monthly.index) <= 2:
                logging.info('Skipping. Less than 2 months of observations.')
                continue

            # Output Folder
            out_folder =  os.path.join(out_dir, 'monthly_comp_plots')
                                       # '{}'.format(
                                       #     row.STATION_ID.replace(" ","")))

            # Create path if it doesn't exist
            if not os.path.exists(out_folder):
                os.makedirs(out_folder)

            # Output to HTML file
            out_file_path = os.path.join(out_folder, '{}.html')\
                .format(row.STATION_ID.replace(" ", ""))
            output_file(out_file_path)

            station_vars = ['TMin (C)', 'TMax (C)', 'wx_Ko_c', 'Rs (w/m2)',
                            'ws_2m (m/s)', 'Vapor Pres (kPa)', 'RHAvg (%)',
                                                               'Precip (mm)',
                            'Calc_ETo (mm)', 'Calc_ETr (mm)']

            gridmet_vars = ['tmin_c', 'tmax_c', 'grid_Ko_c', 'srad_wm2',
                            'u2_ms',
                            'ea_kpa', 'rh_avg', 'prcp_mm', 'eto_mm', 'etr_mm']

            # list of x variables
            x_var_list= station_vars
            # list of y variables
            y_var_list= gridmet_vars
            # title list
            title_list= ['TMin: Monthly Average', 'TMax: Monthly Average',
                         'Ko: Monthly Average' ,
                         'Rs: Monthly Average: Monthly Average',
                         'WS 2m: Monthly Average',
                         'ea: Monthly Average',
                         'RH Avg: Monthly Average', 'Prcp: Monthly Total',
                         'ETo: Monthly Average', 'ETr: Monthly Average']
            # timeseries y label list
            ts_ylabel_list = ['TMin (C)', 'TMax (C)', 'Ko (C)', 'Rs (w/m2)',
                              'WS 2m (m/s)', 'ea (kPa)', 'Avg RH (%)',
                              'Prcp (mm)',
                              'ETo (mm)', 'ETr (mm)']
            # scatter xlabel list
            xlabel_list= ['Station TMin (C)', 'Station TMax (C)',
                          'Station Ko (C)','Station Rs (w/m2)',
                          'Station WS 2m (m/s)', 'Station ea (kPa)',
                          'Station RH (%)', 'Station Prcp (mm)',
                          'Station ETo (mm)', 'Station ETr (mm)']
            # scatter ylabel list
            ylabel_list=['gridMET TMin (C)', 'gridMET TMax (C)',
                         'gridMET Ko (C)','gridMET Rs (w/m2)',
                         'gridMET WS 2m (m/s)', 'gridMET ea (kPa)',
                         'gridMET RH (%)', 'gridMET Prcp (mm)',
                         'gridMET ETo (mm)', 'gridMET ETr (mm)']
            stat_list = ['mean','mean','mean','mean',
                         'mean','mean','mean', 'sum',
                         'mean','mean']
            # legendx list
            legendx_list = ['Station'] * len(title_list)
            # legend y list
            legendy_list = ['gridMET'] * len(title_list)

            # empty list to append figures to
            figure_list = []

            # loop through and create figures for each variable using vars
            # and plot labels from lists above
            for i, (x_var, y_var, title, ts_ylabel, xlabel, ylabel, legendx,
                    legendy, stat) in enumerate(zip(x_var_list, y_var_list,
                                              title_list, ts_ylabel_list,
                                              xlabel_list, ylabel_list,
                                              legendx_list, legendy_list,
                                                    stat_list)):

                # lstsq cannot have nans (drop nas for each var separately)
                monthly2 = monthly[[x_var, y_var]]
                monthly2 = monthly2.dropna()

                if monthly2.empty:
                    logging.info("Skipping {}. No Data.".format(x_var))
                    continue


                if i == 0:
                    # Initial timeseries plot to establish xrange for link axes
                    p1 = figure(plot_width=800, plot_height=400,
                                x_axis_type="datetime",title = title,
                                y_axis_label = ts_ylabel)
                    p1.line(monthly2.index.to_pydatetime(),
                            monthly2[x_var, stat],  color="navy",
                            alpha=0.5, legend=legendx,line_width=2)
                    p1.line(monthly2.index.to_pydatetime(),
                            monthly2[y_var, stat],  color="red",
                            alpha=0.5, legend=legendy,line_width=2)
                else:
                    # Timeseries plots after first pass
                    p1 = figure(plot_width=800, plot_height=400,
                                x_axis_type="datetime",title = title,
                                y_axis_label = ts_ylabel,
                                x_range=p1.x_range)
                    p1.line(monthly2.index.to_pydatetime(),
                            monthly2[x_var, stat],  color="navy", alpha=0.5,
                            legend=legendx,line_width=2)
                    p1.line(monthly2.index.to_pydatetime(),
                            monthly2[y_var, stat],  color="red", alpha=0.5,
                            legend=legendy,line_width=2)

                # 1 to 1 Plot
                # Regression through Zero
                # https://stackoverflow.com/questions/9990789/how-to-force-
                # zero-interception-in-linear-regression/9994484#9994484
                m = np.linalg.lstsq(monthly2[x_var, stat].values.reshape(-1,1),
                                    monthly2[y_var, stat], rcond=None)[0][0]
                r_x, r_y = zip(*((i, i*m ) for i in range(
                    int(np.min([monthly2[y_var, stat],
                                monthly2[x_var, stat]])-2),
                    int(np.max([monthly2[y_var, stat],
                                monthly2[x_var, stat]])+3),1)))
                # Plots
                p2 = figure(plot_width=400, plot_height=400,
                            x_axis_label = xlabel, y_axis_label = ylabel,
                            title = 'Slope Through Zero: m = {}'.format(
                                round(m,4)))
                p2.circle(monthly2[x_var, stat], monthly2[y_var, stat],
                          size=15, color="navy", alpha=0.5)
                p2.line([int(np.min([monthly2[y_var, stat],
                                     monthly2[x_var, stat]])-2),int(np.max(
                    [monthly2[y_var, stat],monthly2[x_var, stat]])+2)],
                         [int(np.min([monthly2[y_var, stat],
                                      monthly2[x_var, stat]])-2),int(np.max(
                             [monthly2[y_var, stat],monthly2[x_var, stat]])+2)],
                          color = "black", legend = '1 to 1 line')
                p2.line(r_x, r_y, color="red", legend = 'Reg thru zero')
                p2.legend.location = "top_left"

                # Append [p1, p2] to figure_list (create list of lists)
                figure_list.append([p1, p2])

            # Plot all figures in list
            fig = gridplot(figure_list, toolbar_location="left")

            # Save the figure
            save(fig)

def station_bar_plot(summary_csv, layer, out_dir=None, x_label=None, 
        y_label=None, title=None, subtitle=None, year_subtitle=True):
    """
    Produce an interactive bar chart comparing multiple climate stations to each
    other for a particular variable, e.g. bias ratios or interpolated residuals.
    
    Arguments:
        summary_csv (str): path to summary CSV produced by either :func:`gridwxcomp.calc_bias_ratios`
            or by :func:`gridwxcomp.interpolate`. Should contain ``layer`` 
            data for plot.
        layer (str): name of variable to plot.
    
    Keyword Arguments:
        out_dir (str or None): default None. Output directory path, default is 
            'station_bar_plots' in parent directory of ``summary_csv``.
        x_label (str or None): default None. Label for x-axis.
        y_label (str or None): default None. Label for y-axis, defaults to 
            ``layer``.
        title (str or None): default None. Title of plot.
        subtitle (str, list, or None): default None. Additional subtitle(s) 
            for plot.
        year_subtitle (bool): default True. If true print subtitle on plot with
            the max year range used for station data, e.g. 'years: 1995-2005'
            
    Example:
        Let's say we want to compare the mean growing seasion bias ratios of 
        reference evapotranspiration (ETr) for the selection of stations we 
        used to calculate bias ratios. The summary CSV file containing the 
        ratios should be first created using :func:`gridwxcomp.calc_bias_ratios`.
        
        >>> from gridwxcomp.plot import station_bar_plot
        >>> # path to summary CSV with station data
        >>> in_file = 'monthly_ratios/etr_mm_summary_all_yrs.csv'
        >>> layer = 'growseason_mean' 
        >>> station_bar_plot(in_file, layer)
        
        The resulting file will be saved using the layer name as a file name::
        
            'monthly_ratios/station_bar_plots/growseason_mean.html'
            
        The plot file will contain the mean growing season bias ratios 
        of ETr for each station, sorted from smallest to largest values. 
        
        This function may also be used for any numerical data in the summary CSV
        files that are created by :func:`gridwxcomp.interpolate` in addition to
        those created by :func:`gridwxcomp.calc_bias_ratios`. The main  
        requirement is that ``summary_csv`` must contain the column 'STATION_ID'
        and the ``layer`` keyword argument.
        
    Raises:
        FileNotFoundError: if ``summary_csv`` is not found.
        KeyError: if ``layer`` does not exist as a column name in ``summary_csv``.
        
    """
    if not Path(summary_csv).is_file():
        err_msg = '\n{} is not a valid path to a summary CSV file!'.\
                format(in_path)
        raise FileNotFoundError(err_msg)
        
    df = pd.read_csv(summary_csv)
    
    if not layer in df.columns:
        err_msg = '\nColumn {} was not found in {}'.format(layer, summary_csv)
        raise KeyError(err_msg)
    
    df.sort_values(layer, inplace=True)
    source = ColumnDataSource(df)
    # hover tooltip with station and value
    tooltips = [
        ("station", "@STATION_ID"),
        ("value", "@{}".format(layer)), 
    ]
    hover = models.HoverTool(tooltips=tooltips)
    
    if not y_label:
        y_label = layer
    # save to working directory in 'station_bar_plots' if not specified
    if not out_dir:
        out_dir = Path(summary_csv).parent/'station_bar_plots'
        if not out_dir.is_dir():
            print('\n{} does not exist, making directory'.format(out_dir))
            out_dir.mkdir(parents=True, exist_ok=True)
    
    out_file = out_dir/'{}.html'.format(layer)
    print('\nCreating station bar plot for variable: ', layer,
         '\nUsing data from file: ', Path(summary_csv).absolute())
    
    output_file(out_file)
    
    p = figure(x_range=df.STATION_ID, y_axis_label=y_label, title=title)
    p.vbar(x='STATION_ID', top=layer, width=0.8, source=source)
    p.xaxis.major_label_orientation = pi/2
    p.add_tools(hover, models.BoxSelectTool())
    
    if year_subtitle:
        # add data range (years start to end) as subtitle 
        min_yr = df.start_year.min().astype(int)
        max_yr = df.end_year.max().astype(int)
        if min_yr == max_yr:
            year_str = 'year: {}'.format(min_yr)
        else:
            year_str = 'years: {}-{}'.format(min_yr, max_yr)
        # caution note if not all stations use full year range
        if not (df.end_year==max_yr).all() or not (df.start_year==min_yr).all():
            year_str = '{} (less years exist for some stations)'.\
                    format(year_str)
            
        p.add_layout(models.Title(text=year_str, text_font_style="italic"), 
                'above')
    # add arbitrary number of custom subtitles as lines above plot 
    if isinstance(subtitle, (list, tuple)):
        for st in subtitle:
            p.add_layout(models.Title(text=st, text_font_style="italic"), 
                    'above')
    elif subtitle:
        p.add_layout(models.Title(text=subtitle, text_font_style="italic"), 
                'above')
        
    save(p)
    print('\nPlot saved to: ', out_file.absolute())


def arg_parse():
    """
    Create time series and  scatter comparison plots between paired station 
    and gridMET climatic variables or bar plots comparing multiple stations 
    for a single variable. 
    """
    parser = argparse.ArgumentParser(
        description=arg_parse.__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    optional = parser._action_groups.pop()  # optionals listed second
    required = parser.add_argument_group('required arguments')
    required.add_argument(
        '-i', '--input', metavar=None, required=True,
        help='Input CSV, if type="station-grid-comp" use file created by '+\
        'download_gridmet_ee.py, if type="station-bar" use file created '+\
        'by calc_bias_ratios.py or spatial.py')
    optional.add_argument(
        '-t', '--plot-type', required=False, type=str, 
        default='station-grid-comp',
        help='Plot type, station comparison bar plot "station-bar" for a '+\
        'single variable, or station to gridMET time series comparison '+\
        'for multiple variables "station-grid-comp"')
    optional.add_argument(
        '-v', '--variable', required=False, type=str, default='annual_mean',
        help='Variable to plot for --type="station-bar" plots, default '+\
                '"annual_mean", can also pass a comma separated list: e.g. '+\
                '"Jan_mean,Jan_res,Jan_est"')
    optional.add_argument(
        '-o', '--out-dir', required=False, type=str, default=None,
        help='Output directory to save comparison plots')
    optional.add_argument(
        '-f', '--freq', required=False, type=str, default='daily',
        help='Time frequency for station-grid-comp plots, "daily" or "monthly"')
    optional.add_argument(
        '--x-label', required=False, type=str, default=None,
        help='X-axis label for station bar plot, when --type="station-bar"')
    optional.add_argument(
        '--y-label', required=False, type=str, default=None,
        help='Y-axis label for station bar plot, when --type="station-bar"')
    optional.add_argument(
        '--title', required=False, type=str, default=None,
        help='Plot title for station bar plot, when --type="station-bar"')
    optional.add_argument(
        '--subtitle', required=False, type=str, default=None,
        help='Plot subtitle for station bar plot, when --type="station-bar" '+\
        ', supports multiple lines as comma separated list')
    optional.add_argument(
        '--year-subtitle', required=False, default=True, action='store_false', 
        help='Add years used for station bar plot as subtitle')
    optional.add_argument(
        '-y', '--year', required=False, default='', type=str,
        help='Year to plot when freq="daily" type="station-grid-comp", '+\
        'single year or range YYYY or YYYY-YYYY')
    optional.add_argument(
        '--debug', default=logging.INFO, const=logging.DEBUG,
        help='Debug level logging', action="store_const", dest="loglevel")
    parser._action_groups.append(optional)  # to avoid optionals listed first
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = arg_parse()

    logging.basicConfig(level=args.loglevel, format='%(message)s')
    logging.info('\n{}'.format('#' * 80))
    logging.info('{0:<20s} {1}'.format(
        'Run Time Stamp:', dt.datetime.now().isoformat(' ')))
    logging.info('{0:<20s} {1}'.format('Current Directory:', os.getcwd()))
    logging.info('{0:<20s} {1}'.format(
        'Script:', os.path.basename(sys.argv[0])))

    if args.plot_type == 'station-grid-comp':
        # all options for time aggregation of plotting data
        time_freqs = ['daily', 'monthly']
        # check plot frequency option
        if not args.freq in time_freqs:
            print(
                '\n{} is not a valid time frequency, available options: {}'.\
                format(args.freq, ', '.join([t for t in time_freqs]))
            )
        elif args.freq == 'daily':
            # call gridwxcomp.daily_comparison
            daily_comparison(input_csv=args.input, out_dir=args.out_dir,
                year_filter=args.year)
        elif args.freq == 'monthly':
            if args.year:
                print('\nWarning: the --year, -y option is not used for'+\
                    ' creating monthly avg. plots, all years will be used.')
            monthly_comparison(input_csv=args.input, out_dir=args.out_dir)

    elif args.plot_type == 'station-bar':
        if args.year:
            print('\nWarning: the --year, -y option is not used for'+\
                ' creating monthly avg. plots, all years will be used.')
        if args.subtitle:
            args.subtitle = args.subtitle.split(',')
        # run one or multiple variables for station bar plots
        variables = args.variable.split(',')
        if len (variables) == 1:
            station_bar_plot(summary_csv=args.input, layer=args.variable,
                out_dir=args.out_dir, x_label=args.x_label, 
                y_label=args.y_label, title=args.title,subtitle=args.subtitle,
                year_subtitle=args.year_subtitle)
        else:
            for v in variables:
                station_bar_plot(summary_csv=args.input, layer=v, 
                    out_dir=args.out_dir, x_label=args.x_label, 
                    y_label=args.y_label, title=args.title, 
                    subtitle=args.subtitle, year_subtitle=args.year_subtitle)
    else:
        print('{} is an invalid option for plot type'.format(args.plot_type),
            '\navailable options include: ', 
            '\n'.join(['station-gridmet','station-station']))

