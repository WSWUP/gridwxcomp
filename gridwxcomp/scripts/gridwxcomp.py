# -*- coding: utf-8 -*-
"""
Command line interface for gridwxcomp.
"""
import click
import os
import logging

from ..calc_bias_ratios import calc_bias_ratios as calc_ratios
from ..daily_comparison import daily_comparison as daily_comp
from ..download_gridmet_ee import download_gridmet_ee as download
from ..prep_input import prep_input as prep
from ..spatial import main as interp #make_points_file, make_grid, interpolate

logging.basicConfig(level=logging.INFO, format='%(message)s')


@click.group()
def gridwxcomp():
    "access gridwxcomp functionality from the command line"
    click.echo('\n*** Welcome to gridwxcomp! ***\n')


@gridwxcomp.command()
@click.argument('station_meta_path', nargs=1)
@click.option('--out-path', '-o', nargs=1, type=str, default='merged_input.csv',
              help='file path to save output CSV')
@click.option('--gridmet-meta', '-g', nargs=1, type=str, default=None,
              help='file path to gridmet_cell_data.csv metadata')
@click.option('--quiet', default=False, is_flag=True, 
        help='supress command line output')
def prep_input(station_meta_path, out_path, gridmet_meta, quiet):
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
        help='folder to save downloaded gridMET time series')
@click.option('--years', '-y', nargs=1, type=str, default=None,
        help='Year(s) to download, single year (YYYY) or range (YYYY-YYYY)')
@click.option('--quiet', default=False, is_flag=True, 
        help='supress command line output')
def download_gridmet_ee(input_csv, out_dir, years, quiet):
    if quiet:
        logging.getLogger().setLevel(logging.ERROR)
    else:
        logging.getLogger().setLevel(logging.INFO)
    # call gridwxcomp.download_gridmet_ee
    download(input_csv, out_dir, year_filter=years)


@gridwxcomp.command()
@click.argument('input_csv', nargs=1)
@click.option('--out-dir', '-o', nargs=1, type=str, default='monthly_ratios',
        help='folder to save correction ratio summary CSV files')
@click.option('--gridmet-var', '-gv', nargs=1, type=str, default='etr_mm',
        help='name of gridMET climatic variable, e.g. etr_mm')
@click.option('--station-var', '-sv', nargs=1, type=str, default=None,
        help='name of station climatic variable')
@click.option('--gridmet-id', '-id', nargs=1, type=str, default=None,
        help='gridMET ID for calculating ratios only for stations in that cell')
@click.option('--day-limit', '-d', nargs=1, type=str, default=10,
        help='day threshold per month, if missing more exclude month from calc')
@click.option('--comp', '-c', default=True, is_flag=True,
        help='flag to NOT save comprehensive output CSV')
@click.option('--quiet', default=False, is_flag=True, 
        help='supress command line output')
def calc_bias_ratios(input_csv, out_dir, gridmet_var, station_var, 
        gridmet_id, day_limit, comp, quiet):
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
        comp=comp
    ) 


@gridwxcomp.command()
@click.argument('summary_comp_csv', nargs=1)
@click.option('--out', '-o', nargs=1, type=str, default=None,
        help='subdirectory for saving interpolated rasters')
@click.option('--layer', '-l', nargs=1, type=str, default='all',
        help='layers to interpolate comma separated, e.g. Jan_mean,Aug_mean')
@click.option('--buffer', '-b', nargs=1, type=int, default=25,
        help='buffer for expanding interpolation region (gridMET cells)')
@click.option('--scale', '-s', nargs=1, type=float, default=0.1,
        help='scale factor for gridMET interpolation, applied to 4 km res')
@click.option('--function', '-f', nargs=1, type=str, default='invdist',
        help='algorithm name for spatial interpolation')
@click.option('--smooth', nargs=1, type=float, default=0, is_flag=False,
        help='smoothing parameter for radial basis funciton interpolation')
@click.option('--params', '-p', nargs=1, type=str, default=None, is_flag=False,
        help='parameters for gdal_grid interpolation e.g. :power=2:smooth=0')
@click.option('--zonal-stats', default=True, is_flag=True,
        help='flag to NOT extract zonal means of interpolated results')
@click.option('--overwrite-grid', default=False, is_flag=True,
        help='flag to overwrite grid for zonal stats if already exists')
@click.option('--options', nargs=1, type=str, default=None, is_flag=False,
        help='extra command line arguments for gdal_grid interpolation')
@click.option('--gridmet-meta', '-g', nargs=1, type=str, default=None,
              help='file path to gridmet_cell_data.csv metadata')
@click.option('--quiet', default=False, is_flag=True, 
        help='supress command line output')
def spatial(summary_comp_csv, layer, out, buffer, scale, function, smooth, 
        params, zonal_stats, overwrite_grid, options, gridmet_meta, quiet):
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
        zonal_stats=zonal_stats,
        overwrite=overwrite_grid,
        options=options,
        gridmet_meta_path=gridmet_meta
    )

@gridwxcomp.command()
@click.argument('input_csv', nargs=1)
@click.option('--out-dir', '-o', nargs=1, type=str, default=os.getcwd(),
        help='folder to save time series comparison plots')
@click.option('--year', '-y', nargs=1, type=int, default=None,
        help='Year to plot, single year (YYYY)')
@click.option('--quiet', default=False, is_flag=True, 
        help='supress command line output')
def daily_comparison(input_csv, out_dir, year, quiet):
    if quiet:
        logging.getLogger().setLevel(logging.ERROR)
    else:
        logging.getLogger().setLevel(logging.INFO)
    # call gridwxcomp.daily_comparison
    daily_comp(input_csv, out_dir, year_filter=year)


