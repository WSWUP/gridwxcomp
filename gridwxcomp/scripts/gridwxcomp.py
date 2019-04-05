# -*- coding: utf-8 -*-
"""
Command line interface for gridwxcomp.
"""
import click
import os
import logging

from gridwxcomp.prep_input import prep_input as prep
from gridwxcomp.download_gridmet_ee import download_gridmet_ee as download
from gridwxcomp.calc_bias_ratios import calc_bias_ratios as calc_ratios
from gridwxcomp.spatial import main as interp 
from gridwxcomp.plot import daily_comparison, monthly_comparison, station_bar_plot

logging.basicConfig(level=logging.INFO, format='%(message)s')


@click.group()
def gridwxcomp():
    """
    Access gridwxcomp functionality from the command line
    
    The ``gridwxcomp`` tools are accesible via a command line interface using
    the console command "gridwxcomp". To print a list of all available commands
    from the command line type,::

        gridwxcomp --help

    Or to print information on a specific command,::

        gridwxcomp COMMAND --help

    Basic usage to run a command: 

    """
    click.echo('\n*** Welcome to gridwxcomp ***\n')


@gridwxcomp.command()
@click.argument('station_meta_path', nargs=1)
@click.option('--out-path', '-o', nargs=1, type=str, default='merged_input.csv',
              help='File path to save output CSV')
@click.option('--gridmet-meta', '-g', nargs=1, type=str, default=None,
              help='File path to gridmet_cell_data.csv metadata')
@click.option('--quiet', default=False, is_flag=True, 
        help='Supress command line output')
def prep_input(station_meta_path, out_path, gridmet_meta, quiet):
    """
    Pairs climate metadata with gridMET

    Read climate station metadata file, e.g. from `PyWeatherQAQC <https://github.com/WSWUP/pyWeatherQAQC>`_
    and match each station with gridMET cell locations, save CSV. The required 
    argument ``STATION_META_PATH`` is a file containing metadata of climate 
    stations, if it was not built from `PyWeatherQAQC <https://github.com/WSWUP/pyWeatherQAQC>`_, 
    the file should be in CSV format and contain at least these four columns:

    * Latitude
    * Longitude
    * Station
    * Filename 

    For example the "Station_Data.txt" file is included as an example at
    `gridwxcomp/gridwxcomp/example_data/ <https://raw.githubusercontent.com/WSWUP/gridwxcomp/master/gridwxcomp/example_data/Station_Data.txt>`_, contains

        ======== ========== ======================= ==================== 
        Latitude Longitude  Station                 Filename             
        ======== ========== ======================= ==================== 
        40.37232 -110.20918 Bluebell (Neola Area)   BluebellUT_Daily     
        38.38346 -111.63583 Loa                     LoaUT_Daily          
        38.32829 -108.85549 Bedrock                 BedrockCO_Daily      
        38.64294 -109.39880 Castle Valley near Moab CastleValleyUT_Daily 
        ======== ========== ======================= ====================
    
    Also, if not created by ``PyWeatherQaQc``, the climate station time series
    files should be in CSV format with a "date" column with datetime formats,
    e.g. '12/01/2018'. By default ``gridwxcomp prep-input`` creates the file 
    "merged_input.csv" in the current working directory which contains 
    metadata from climate stations as well as the lat, long, and gridMET ID 
    (cell index) of the nearest gridMET cell's centroid.
    
    """
    if quiet:
        logging.getLogger().setLevel(logging.ERROR)
    else:
        logging.getLogger().setLevel(logging.INFO)
    # calling gridwxcomp.prep_input 
    prep(
            station_meta_path, 
            out_path=out_path, 
            gridmet_meta_path=gridmet_meta
        )

@gridwxcomp.command()
@click.argument('input_csv', nargs=1)
@click.option('--out-dir', '-o', nargs=1, type=str, default='gridmet_data',
        help='Folder to save downloaded gridMET time series')
@click.option('--years', '-y', nargs=1, type=str, default=None,
        help='Year(s) to download, single year (YYYY) or range (YYYY-YYYY)')
@click.option('--update-years', '-u', nargs=1, type=str, default=None,
        help='Year(s) to redownload or update, YYYY or YYYY-YYYY')
@click.option('--quiet', default=False, is_flag=True, 
        help='Supress command line output')
def download_gridmet_ee(input_csv, out_dir, years, update_years, quiet):
    """
    Download gridMET climate time series

    Download gridMET time series for cells that are paired to climate stations
    in a CSV file that is first created by ``gridwxcomp prep-input``. Options 
    allow for downloading all years available (default) or select years, it is 
    also possible to redownload data for specified year(s). Uses the Google 
    Earth Engine Python API. If ``--out-dir`` is not specified, gridMET time 
    series CSVs are saved to a new directory named "gridmet_data" within the
    current working directory.
    """
    if quiet:
        logging.getLogger().setLevel(logging.ERROR)
    else:
        logging.getLogger().setLevel(logging.INFO)
    # call gridwxcomp.download_gridmet_ee
    download(input_csv, out_dir, year_filter=years, year_update=update_years)


@gridwxcomp.command()
@click.argument('input_csv', nargs=1)
@click.option('--out-dir', '-o', nargs=1, type=str, default='monthly_ratios',
        help='Folder to save correction ratio summary CSV files')
@click.option('--gridmet-var', '-gv', nargs=1, type=str, default='etr_mm',
        help='Name of gridMET climatic variable, e.g. etr_mm')
@click.option('--station-var', '-sv', nargs=1, type=str, default=None,
        help='Name of station climatic variable')
@click.option('--gridmet-id', '-id', nargs=1, type=str, default=None,
        help='gridMET ID for calculating ratios only for stations in that cell')
@click.option('--day-limit', '-d', nargs=1, type=str, default=10,
        help='Monthly day threshold, if missing more exclude month from calc')
@click.option('--years', '-y', nargs=1, type=str, default='all',
        help='Year(s) to use, e.g. 2010 or 2000-2010')
@click.option('--comp', '-c', default=True, is_flag=True,
        help='Flag to NOT save comprehensive output CSV')
@click.option('--quiet', default=False, is_flag=True, 
        help='Supress command line output')
def calc_bias_ratios(input_csv, out_dir, gridmet_var, station_var, 
        gridmet_id, day_limit, years, comp, quiet):
    """
    Bias ratio statistics of station-to-gridMET 

    Calculates mean bias ratios between climate station to gridMET variables 
    at multiple time periods: monthly, annual, growing season. Other statistics
    include day counts, standard deviation, and coefficient of variation. 
    Options allow for excluding statistics on months with less than a certain 
    number of days (default 10). gridMET variable names available and their 
    default paired climate variable names include:

        ==========  ==================
        gridMET     station
        ==========  ==================
        'u2_ms'     'ws_2m (m/s)'
        'tmin_c'    'TMin (C)'
        'tmax_c'    'TMax (C)'
        'srad_wm2'  'Rs (w/m2)'
        'ea_kpa'    'Vapor Pres (kPa)'
        'prcp_mm'   'Precip (mm)'
        'etr_mm'    'Calc_ETr (mm)'
        'eto_mm'    'Calc_ETo (mm)'
        ==========  ==================

    If ``--station-var`` is given from the default options, the corresponding
    default gridMET variable will be used, otherwise you may assign these based
    on the variable names in your station or gridMET time series files. 
    """
    if quiet:
        logging.getLogger().setLevel(logging.ERROR)
    else:
        logging.getLogger().setLevel(logging.INFO)
    # call gridwxcomp.calc_bias_ratios
    calc_ratios(
        input_csv, 
        out_dir, 
        gridmet_var=gridmet_var,
        station_var=station_var, 
        gridmet_ID=gridmet_id, 
        day_limit=day_limit,
        years=years,
        comp=comp
    ) 


@gridwxcomp.command()
@click.argument('summary_comp_csv', nargs=1)
@click.option('--out', '-o', nargs=1, type=str, default=None,
        help='Subdirectory for saving interpolated rasters')
@click.option('--layer', '-l', nargs=1, type=str, default='all',
        help='Layers to interpolate comma separated, e.g. Jan_mean,Aug_mean')
@click.option('--buffer', '-b', nargs=1, type=int, default=25,
        help='Buffer for expanding interpolation region (gridMET cells)')
@click.option('--scale', '-s', nargs=1, type=float, default=0.1,
        help='Scale factor for gridMET interpolation, applied to 4 km res')
@click.option('--function', '-f', nargs=1, type=str, default='invdist',
        help='Algorithm name for spatial interpolation')
@click.option('--smooth', nargs=1, type=float, default=0, is_flag=False,
        help='Smoothing parameter for radial basis funciton interpolation')
@click.option('--params', '-p', nargs=1, type=str, default=None, is_flag=False,
        help='Parameters for gdal_grid interpolation e.g. :power=2:smooth=0')
@click.option('--no-zonal-stats', '-z', default=True, is_flag=True,
        help='Flag to NOT extract zonal means of interpolated results')
@click.option('--overwrite-grid', default=False, is_flag=True,
        help='Flag to overwrite grid for zonal stats if already exists')
@click.option('--options', nargs=1, type=str, default=None, is_flag=False,
        help='Extra command line arguments for gdal_grid interpolation')
@click.option('--gridmet-meta', '-g', nargs=1, type=str, default=None,
              help='file path to gridmet_cell_data.csv metadata')
@click.option('--quiet', default=False, is_flag=True, 
        help='Supress command line output')
def spatial(summary_comp_csv, layer, out, buffer, scale, function, smooth, 
        params, no_zonal_stats, overwrite_grid, options, gridmet_meta, quiet):
    """
    Spatially interpolate ratio statistics 

    2-D interpolation of mean bias ratios or other statistics that exist in the
    ``SUMMARY_COMP_CSV`` file created by ``gridwxcomp calc-bias-ratios``. 
    Interpolation algorithms include those provided by `gdal_grid <https://www.gdal.org/gdal_grid.html>`_:
    'invdist' (default), 'invdistnn', 'linear', 'average', and 'nearest' 
    and :class:`scipy.interpolate.Rbf`: 'multiquadric', 'inverse', 'gaussian', 
    'linear_rbf', 'cubic', 'quintic', and 'thin_plate'. Parameters for gdal 
    algorithms are set using ``--params`` whereas the ``--smooth`` parameter 
    (default 0) is used by scipy radial basis functions. Default parameters 
    for gdal are stored in :attr:`gridwxcomp.InterpGdal.default_params`. 
    
    Interpolation resampling resolution can be down- or up-scaled using 
    ``--scale-factor`` (default 0.1) which is applied to gridMET resolution of 
    4 km. 
    
    On the first run this command will produce a fishnet of gridMET cells
    that bound the climate stations in ``SUMMARY_COMP_CSV`` with optional 
    buffer of gridMET cells which can later be overwritten. 
    
    All mean bias ratios found within :attr:`gridwxcomp.InterpGdal.default_layers` 
    are interpolated by default unless ``--layers`` are specified. Interpolated 
    residuals to station point ratios are also estimated. Zonal means for
    gridMET cells in the fishnet grid are extracted by default and assigned to
    each gridMET cell by its index (GRIDMET_ID) in the full gridMET fishnet 
    that covers the CONUS. In total interpolated rasters, station point 
    shapefiles, fishnet grid for zonal stats, and CSVs of bias ratios and zonal
    statistics are all created and stored in a file structure that is explained
    in :func:`gridwxcomp.spatial.make_grid`. and :func:`gridwxcomp.spatial.interpolate`.
    """
    if quiet:
        logging.getLogger().setLevel(logging.ERROR)
    else:
        logging.getLogger().setLevel(logging.INFO)

    # parse multiple layers option from comma separated string
    layer = layer.split(',')
    if len(layer) == 1:
        layer = layer[0] # get single layer as string

    # call gridwxcomp.spatial.main
    interp(
        summary_comp_csv, 
        layer=layer,
        out=out,
        buffer=buffer,
        scale_factor=scale,
        function=function,
        smooth=smooth,
        params=params,
        zonal_stats=no_zonal_stats,
        overwrite=overwrite_grid,
        options=options,
        gridmet_meta_path=gridmet_meta
    )

@gridwxcomp.command()
@click.argument('input_csv', nargs=1)
@click.option('--plot-type', '-t', nargs=1,type=str,default='station-grid-comp',
        help='Plot type, station comparison bar plot "station-bar" for a '+\
        'single variable, or station to gridMET time series comparison '+\
        'for multiple variables "station-grid-comp"')
@click.option('--variable', '-v', nargs=1, type=str, default='annual_mean',
        help='Variable to plot for --type="station-bar" plots, default '+\
        '"annual_mean", can also pass a comma separated list: e.g. '+\
        '"Jan_mean,Jan_res,Jan_est"')
@click.option('--freq', '-f', nargs=1, type=str, default='daily',
        help='Time frequency for station-grid-comp plots "daily" or "monthly"')
@click.option('--out-dir', '-o', nargs=1, type=str, default=None,
        help='Folder to save plots')
@click.option('--year', '-y', nargs=1, type=int, default=None,
        help='Years for time series, single year YYYY or range YYYY-YYYY')
@click.option('--x-label', nargs=1, type=str, default=None,
        help='X-axis label for station-bar plot')
@click.option('--y-label', nargs=1, type=str, default=None,
        help='Y-axis label for station-bar plot')
@click.option('--title', nargs=1, type=str, default=None,
        help='Title for station-bar plot')
@click.option('--subtitle', nargs=1, type=str, default=None,
        help='Subtitle for station-bar plot, comma separated list for multiple')
@click.option('--year-subtitle', default=True, is_flag=True, 
        help='Print subtitle with years used for station-bar plots')
@click.option('--quiet', default=False, is_flag=True, 
        help='Supress command line output')
def plot(input_csv, plot_type, freq, variable, out_dir, year, x_label,
        y_label, title, subtitle, year_subtitle, quiet):
    """
    Create comparison or diagnostic graphics

    Plot time series and scatter comparison plots between station and gridMET
    paired variables at daily or monthly average time periods. Also allows for
    plotting a single variable as a bar plot that compares multiple climate
    stations. Plots are interactive (can be panned, zoomed, etc.) HTML files.

    Currently there are two major plot types provided which are specified by
    ``[-t, --plot-type]``: 1) "station-grid-comp" which creates time series 
    comparison plots of all gridMET variables to station values; and 2)
    the second type "station-bar" which compares stations for a single variable.
    For time series comparisons a plot is created for each station if time 
    frequency of plots, i.e. ``--freq`` is "monthly" and for each month for 
    each station if ``--freq`` is "daily" (default). For producing a station
    bar plot you must assign the variable to plot using ``[-v, --variable]``.

    By default, if ``--out-dir`` is not specified time series plots are saved 
    in "[freq]_comp_plots" in the current directory where [freq] is set to 
    ``--freq``. Station bar plots are saved by default to "station_bar_plots" 
    directory created in the same directory containing the summary CSV assigned
    to ``INPUT_CSV``.

    **Note:**
    Both daily and monthly comparison plots between station and gridMET
    i.e. when ``plot-type`` = "station-grid-comp" require station data to
    be within excel files that follow the format produced by 
    `PyWeatherQAQC <https://github.com/WSWUP/pyWeatherQAQC>`_. 
    """
    if quiet:
        logging.getLogger().setLevel(logging.ERROR)
    else:
        logging.getLogger().setLevel(logging.INFO)

    if plot_type == 'station-grid-comp':
        # all options for time aggregation of plotting data
        time_freqs = ['daily', 'monthly']
        # check plot frequency option
        if not freq in time_freqs:
            click.echo(
                '\n{} is not a valid time frequency, available options: {}'.\
                format(freq, ', '.join([t for t in time_freqs]))
            )
        elif freq == 'daily':
            # call gridwxcomp.daily_comparison
            daily_comparison(input_csv=input_csv, out_dir=out_dir,
                year_filter=year)
        elif freq == 'monthly':
            if year:
                click.echo('\nWarning: the --year, -y option is not used for'+\
                    ' creating monthly avg. plots, all years will be used.')
            monthly_comparison(input_csv=input_csv, out_dir=out_dir)

    elif plot_type == 'station-bar':
        if year:
            click.echo('\nWarning: the --year, -y option is not used for'+\
                ' creating monthly avg. plots, all years will be used.')
        if subtitle:
            subtitle = subtitle.split(',')
        # run one or multiple variables for station bar plots
        variables = variable.split(',')
        if len (variables) == 1:
            station_bar_plot(summary_csv=input_csv, layer=variable,
                out_dir=out_dir, x_label=x_label, y_label=y_label,
                title=title, subtitle=subtitle, 
                year_subtitle=year_subtitle)
        else:
            for v in variables:
                station_bar_plot(summary_csv=input_csv, layer=v, 
                    out_dir=out_dir, x_label=x_label, 
                    y_label=y_label, title=title, 
                    subtitle=subtitle, year_subtitle=year_subtitle)
    else:
        click.echo('{} is an invalid option for plot type'.format(plot_type),
            '\navailable options include: ', 
            '\n'.join(['station-gridmet','station-station']))


