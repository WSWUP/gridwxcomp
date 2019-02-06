# -*- coding: utf-8 -*-
"""
Perform multiple workflows needed to estimate the spatial surface of monthly
bias ratios between climate station and gridMET estimated of reference 
evapotranspiration. The input file is first created by :mod:`calc_bias_ratios`.

Attributes:
    CELL_SIZE (float): constant gridMET cell size in decimal degrees.

Note:
    All spatial files, i.e. vector and raster files, utilize the
    *ESRI Shapefile* file format. 

Todo:
    * add extraction of zonal stats to gridMET cells
    * logging 

"""

import os
import re
import argparse
import copy
from math import ceil
from shutil import move

import fiona
import numpy as np
import pandas as pd
from scipy.interpolate import Rbf
from shapely.geometry import Point, Polygon, mapping
from fiona import collection
from fiona.crs import from_epsg
from osgeo import gdal, osr, ogr
from rasterstats import zonal_stats

CELL_SIZE = 0.041666666666666664

OPJ = os.path.join

def main(input_file_path, gridmet_meta_path=None, 
         buffer=25, scale_factor=0.1, function='inverse'):
    """
    Create point shapefile of monthly mean bias ratios from comprehensive
    CSV file created by :mod:`calc_bias_ratios`. Build fishnet grid around
    the climate stations that matches the gridMET dataset. Perform spatial
    interpolation to estimate 2-dimensional surface of bias ratios and extract
    interpolated values to gridMET cells. 
    
    Arguments:
        input_file_path (str): path to summary_comp CSV file containing monthly
            bias ratios, lat, long, and other data. Shapefile "summary_pts.shp"
            is saved to parent directory of ``input_file_path`` under a 
            subdirectory named "spatial".

    Keyword Arguments:
        gridmet_meta_path (str): default None. Path to metadata CSV file 
            that contains all gridMET cells for the contiguous United 
            States. If None it is looked for at 
            ``etr-biascorrect/gridmet_cell_data.csv``.
        buffer (int): number of gridMET cells to expand the rectangular extent
            of the subgrid fishnet (default=25). 
        scale_factor (float, int): scaling factor to apply to original
            gridMET fishnet to create resampling fishnet. For example,
            if scale_factor = 0.1, the resolution will be one tenth of 
            the original girdMET resolution resulting in a 400 m resolution.
        function (str): default 'inverse', radial basis function to use 
            for interpolation. Options include 
            * 'multiquadric'
            * 'inverse'
            * 'gaussian'
            * 'linear'
            * 'cubic'
            * 'quintic'
            * 'thin_plate'
            
    Returns:
        None

    Examples:
        From the command line interface,

        .. code::
            $ python spatial.py -i monthly_ratios/summary_comp.csv -g gridmet_cell_data.csv

        This will produce a subdirectory under "monthly_ratios" named
        "spatial" where the shapefile will be saved as "summary_pts.shp".
        It will also create a spatially referenced fishnet of gridMET
        cells (at "monthly_ratios/spatial/grid.shp") with a buffer of 25 
        gridMET cells added around the encompassed climate stations. Next
        a 2-dimensional surface is interpolated from the point data of each
        mean bias ratio in summary_comp.csv which includes monthly as well 
        as growing season and annual means. The surfaces are saved to
        "monthly_ratios/spatial" with file names of the form::
        
            [time_period]_[resolution]m.tiff
            
        where [time_period] refers to the bias ratio time aggregate e.g.
        Apr_mean or Annual_mean. [resolution] is the pixel size in meters
        of the spatially interpolated surface. Zonal mean statistics are also 
        calculated for each gridMET cell in the fishnet that is created.
        The mean gridMET cell values for monthly mean, growing season, and
        annual mean ratios are saved to the file::
        
            'monthly_ratios/gridmet_summary_inverse_400m.csv'
            
        The CSV name contains information on the radial basis interpolation
        function and the raster resolution used to calculate zonal statistics.
        The contents of the file will look something like::
        
            ==========  ========  ========  ========  ========
            GRIDMET_ID  Apr_mean  Aug_mean  Dec_mean  ...
            ==========  ========  ========  ========  ========
            515902      1.028707  0.856440  1.058291  ...
            514516      1.035543  0.862066  1.065963  ...
            ...         ...       ...       ...       ...
            ==========  ========  ========  ========  ========
        
        Alternatively, to use a smaller buffer zone on the fishnet grid 
        used for interpolation if, for example, extrapolation is not needed, 
        use the -b or --buffer option

        .. code::
            $ python spatial.py -i monthly_ratios/summary_comp.csv -g gridmet_cell_data.csv -b 30

        Or if we wanted to interpolate to a 200 m resolution 
        (i.e. scale_factor = 0.05, 0.05 x 4 km = 200 m) using the 'linear'
        radial basis function, assuming "gridmet_cell_data.csv" is in the 
        current working directory
        
        .. code::
            $ python spatial.py -i monthly_ratios/summary_comp.csv -s 0.05 -f 'linear'        

        In this case the final zonal statistics of spatial mean bias ratios 
        will be saved to::
        
            'monthly_ratios/gridmet_summary_linear_200m.csv'

        For use within Python see examples of :func:`make_points_file`, 
        :func:`build_subgrid`, and :func:`interpolate`.

    """
    # build point shapefile 
    make_points_file(input_file_path)
    # build fishnet for interpolation
    build_subgrid(
        input_file_path, 
        gridmet_meta_path=gridmet_meta_path,
        buffer=buffer
    )
    # get column names from summary csv for interpolation (mean ratios)
    in_df = pd.read_csv(input_file_path)
    cols = [c for c in in_df.columns if '_mean' in c and not 'June_to' in c]
    # 2D interpolate using RBF, save geoTIFFs for each mean ratio
    for v in cols:
        interpolate(
            input_file_path,                 
            var_name=v, 
            scale_factor=scale_factor, 
            function=function, 
            gridmet_meta_path=gridmet_meta_path,
            buffer=buffer
        ) 

def make_points_file(in_path):
    """
    Create vector shapefile of points with monthly mean bias ratios 
    for climate stations using all stations found in a comprehensive
    CSV file created by :mod:`calc_bias_ratios`.
    
    Arguments:
        in_path (str): path to summary_comp.CSV file containing monthly
            bias ratios, lat, long, and other data. Shapefile 
            "summary_pts.shp" is saved to parent directory of ``in_path`` 
            under a subdirectory named "spatial".
            
    Returns:
        None
    
    Example:
        Create shapefile containing point data for all climate stations
        in input summary file created by :mod:`etr_biascorrect`
        
        >>> from etr_biascorrect import spatial
        >>> # path to comprehensive summary CSV
        >>> summary_file = 'monthly_ratios/summary_comp.csv'
        >>> spatial.make_points_file(summary_file)
        
        The result is the file "summary_pts.shp" being saved to a
        subdirectory created in the directory containing ``in_path``
        named "spatial", i.e. "monthly_ratios/spatial/summary_pts.shp". 
        

    Raises:
        FileNotFoundError: if input summary CSV file is not found.
        
    Note:
        In order to create the point shapefile the "comprehensive" output
        summary file needs to produced from the previous step using the
        "-c" command line flag or comp=True keyword argument to
        :mod:`etr_biascorrect.calc_bias_ratios`. The summary_comp.csv file 
        that is produced from that step is the input to this function.
        
    """
    if not os.path.isfile(in_path):
        raise FileNotFoundError('Input summary CSV file given'+\
                                ' was invalid or not found')
    print(
        '\nMapping point data for climate stations in: \n',
        in_path, '\n'
    )
    in_df = pd.read_csv(in_path, index_col='STATION_ID')
    # save shapefile to "spatial" subdirectory of in_path
    path_root = os.path.split(in_path)[0]
    out_dir = OPJ(path_root, 'spatial')
    out_file = OPJ(out_dir, 'summary_pts.shp')
    print(            
        'Creating point shapefile of station bias ratios, saving to: \n',
        os.path.abspath(out_file),
        '\n'
    )
    # create output directory if does not exist
    if not os.path.isdir(out_dir):
        print(
            out_dir, 
            ' does not exist, creating directory.\n'
        )
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
        # loop through stations and add point data to shapefile
        for index, row in in_df.iterrows():
            print(
                'Saving point data for station: ',
                index, 
            )
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

def get_subgrid_bounds(in_path, buffer):
    """
    Calculate bounding box for spatial interpolation grid from 
    comprehensive summary CSV file containing monthly bias ratios
    and station locations. 
    
    Arguments:
        in_path (str): path to comprehensive summary file containing
            monthly bias ratios, created by :mod:`etr_biascorrect.calc_bias_ratios.py`.
        buffer (int): number of gridMET cells to expand the rectangular extent
            of the subgrid fishnet. 
        
    Returns:
        bounds (tuple): tuple with coordinates in decimal degrees that
            define the outer bounds of the subgrid fishnet in the format
            (min long, max long, min lat, max lat)

    Raises:
        FileNotFoundError: if input summary CSV file is not found.

    Note:
        By expanding the grid to a larger area encompassing the climate 
        stations of interest :func:`interpolate` can be used to extrapolate
        passed the bounds of the outer station locations.
        
    """
    if not os.path.isfile(in_path):
        raise FileNotFoundError('Input summary CSV file given'+\
                                ' was invalid or not found')

    in_df = pd.read_csv(in_path)
    
    # values are centroids of gridmet, get corners
    lons = in_df.sort_values('LON')['LON'].values
    lon_min = lons[0] - CELL_SIZE / 2
    lon_max = lons[-1] + CELL_SIZE / 2

    lats = in_df.sort_values('LAT')['LAT'].values
    lat_min = lats[0] - CELL_SIZE / 2
    lat_max = lats[-1] + CELL_SIZE / 2
    
    # expand bounding extent based on buffer cells
    lon_min -= CELL_SIZE*buffer
    lat_min -= CELL_SIZE*buffer
    lon_max += CELL_SIZE*buffer
    lat_max += CELL_SIZE*buffer    
    
    bounds = lon_min, lon_max, lat_min, lat_max
    
    return bounds

def make_grid(out_path, bounds):
    """
    Make fishnet grid (vector file with polygon geometry) for 
    select gridMET cells based on bounding coordinates. The grid
    is later used to spatially interpolate monthly bias ratios. 
    
    Slightly modified from
    https://pcjericks.github.io/py-gdalogr-cookbook/vector_layers.html#create-fishnet-grid
    
    Arguments:
        out_path (str): file path to save subgrid fishnet shapefile. 
        bounds (tuple): tuple of bounding coordinates in the following order
            (min long, max long, min lat, max lat) which need to be in 
            decimal degrees.
            
    Returns:
        None
        
    Examples:
        Build a fishnet of all gridMET cells that contain climate stations
        found in a summary CSV file of monthly bias ratios that was created 
        using :mod:`etr_biascorrect`.
        
        >>> from etr_biascorrect import spatial
        >>> grid_file = 'grid.shp'    
        >>> # path to summary CSV, most be comprehensive version
        >>> summary_file = 'monthly_ratios/summary_comp.csv'
        >>> # calculate bounding box with extent limited (no buffer cells)
        >>> bnds = get_subgrid_bounds(summary_file, buffer=0)
        >>> # build grid and save
        >>> spatial.make_subgrid(grid_file, bnds)
    
        Build a 5 by 5 cell grid with the lower left corner matching 
        the full gridMET fishnet. 
        
        >>> grid_file = 'grid_5_by_5.shp'
        >>> # get gridMET cell length from spatial module
        >>> CELL_SIZE = spatial.CELL_SIZE
        >>> lon_min = -124.78749996666667
        >>> lon_max = lon_min + CELL_SIZE * 5
        >>> lat_min = 25.04583333333334
        >>> lat_max = lat_min + CELL_SIZE * 5
        >>> bnds = (lon_min, lon_max, lat_min, lat_max)
        >>> # build grid and save 
        >>> spatial.make_subgrid(grid_file, bnds)
            
    Note:
        The fishnet created will not contain gridMET IDs or any other
        attributes for grid cells nor will it be assigned a spatial
        reference system. These are added by subsequently calling 
        :func:`update_subgrid` or by running :mod:`spatial` from
        the command line, see :func:`main`.
        
    """
    
    xmin, xmax, ymin, ymax = bounds
    # read path and make parent directories if they don't exist
    path_root = os.path.split(out_path)[0]
    
    if not os.path.isdir(path_root):
        print(
            os.path.abspath(path_root), 
            ' does not exist, creating directory.\n'
        )
        os.makedirs(path_root)
        
    print(
        '\nCreating fishnet grid for subset of gridMET cells, \n',
        '\nSouthwest corner (lat, long): {:9.4f}, {:9.4f}'.format(ymin, xmin),
        '\nNortheast corner (lat, long): {:9.4f}, {:9.4f}'.format(ymax, xmax),
    )
    
    # get rows
    rows = ceil((ymax-ymin) / CELL_SIZE)
    # get columns
    cols = ceil((xmax-xmin) / CELL_SIZE)

    # start grid cell envelope
    ringXleftOrigin = xmin
    ringXrightOrigin = xmin + CELL_SIZE
    ringYtopOrigin = ymax
    ringYbottomOrigin = ymax - CELL_SIZE
    
    # add spatial reference system WGS 84, add argument to CreateLayer
    #dest_srs = ogr.osr.SpatialReference()
    #dest_srs.ImportFromEPSG(4326)
    
    # create output file
    outDriver = ogr.GetDriverByName('ESRI Shapefile')
    if os.path.exists(out_path):
        os.remove(out_path)
    outDataSource = outDriver.CreateDataSource(out_path)
    outLayer = outDataSource.CreateLayer(out_path,geom_type=ogr.wkbPolygon )
    featureDefn = outLayer.GetLayerDefn()

    # create grid cells
    countcols = 0
    while countcols < cols:
        countcols += 1

        # reset envelope for rows
        ringYtop = ringYtopOrigin
        ringYbottom = ringYbottomOrigin
        countrows = 0

        while countrows < rows:
            countrows += 1
            ring = ogr.Geometry(ogr.wkbLinearRing)
            ring.AddPoint(ringXleftOrigin, ringYtop)
            ring.AddPoint(ringXrightOrigin, ringYtop)
            ring.AddPoint(ringXrightOrigin, ringYbottom)
            ring.AddPoint(ringXleftOrigin, ringYbottom)
            ring.AddPoint(ringXleftOrigin, ringYtop)
            poly = ogr.Geometry(ogr.wkbPolygon)
            poly.AddGeometry(ring)

            # add new geom to layer
            outFeature = ogr.Feature(featureDefn)
            outFeature.SetGeometry(poly)
            outLayer.CreateFeature(outFeature)
            outFeature = None

            # new envelope for next poly
            ringYtop = ringYtop - CELL_SIZE
            ringYbottom = ringYbottom - CELL_SIZE

        # new envelope for next poly
        ringXleftOrigin = ringXleftOrigin + CELL_SIZE
        ringXrightOrigin = ringXrightOrigin + CELL_SIZE

    # Save and close DataSources
    outDataSource = None
    
    print(
        '\nFishnet shapefile successfully saved to: \n',
        os.path.abspath(out_path),
        '\n'
    )

def get_cell_ID(coords, cell_data):
    """
    Helper function that calculates the gridMET ID for gridMET cells using
    the coordinates of a cell in the fishnet grid. Uses the
    :class:`Fiona.geometry.Polygon` object to calculate centroids for
    gridMET cells, then lookup gridMET ID from gridmet_cells_data.csv
    
    Arguments:
        coords (list): list of coordinates that constitute a polygon for a 
            gridMET cell. Coordinates must be acceptable as arguments to
            :class:`fiona.geometry.Polygon`.
        cell_data (:obj:`pandas.DataFrame`): pandas dataframe of the 
            gridMET metadata CSV file, i.e. "gridmet_cell_data.csv".
        
    Returns:
        gridmet_id (int): gridMET ID value for gridMET cell that is defined
            by the bounding polygon defined by ``coords``.
            
    """
    poly = Polygon(coords)
    lon_c = poly.bounds[0] + CELL_SIZE / 2
    lat_c = poly.bounds[1] + CELL_SIZE / 2
    
    # use centroid of cell and centroid coords in cell_data
    row = cell_data.loc[
            (np.isclose(cell_data.LON, lon_c)) & (np.isclose(cell_data.LAT, lat_c))
            ]
    
    # if cell falls outside of master gridMET fishnet assign -999 id
    if len(row.GRIDMET_ID.values) == 1:
        gridmet_id = int(row.GRIDMET_ID.values[0])
    else:
        print(
            'Cell centroid (lat, long): {:9.4f}, {:9.4f}'.format(lat_c, lon_c),
            '\nfalls outside of the gridMET dataset, ',
            'assigning GRIDMET_ID attribute -999'
        )
        gridmet_id = -999
    
    return gridmet_id

def update_subgrid(grid_path, gridmet_meta_path=None):
    """
    Assign gridMET ID values and EPSG WGS 84 projection to fishnet grid.
    
    Arguments:
        grid_path (str): path to fishnet grid shapefile for subset of
            gridMET cells which contain climate stations.
            
    Keyword Arguments:
        gridmet_meta_path (str): default None. Path to metadata CSV file 
            that contains all gridMET cells for the contiguous United 
            States. If None it is looked for at 
            ``etr-biascorrect/gridmet_cell_data.csv``.
            
    Returns:
        None
        
    Example:
        First build the subgrid fishnet of gridMET cells of interest 
        by calling :func:`make_subgrid`. Next run this function to
        add gridMET ID values and spatial projection WGS 84.
        
        >>> from etr_biascorrect import spatial
        >>> # set path to existing fishnet shapefile
        >>> grid_file = 'grid.shp' 
        >>> spatial.update_subgrid(grid_file)
        
        Note that ``gridmet_meta_path`` keyword argument was not 
        given, therefore it must be in the current working directory
        and named "gridmet_cell_data.csv".
        
    Raises:
        FileNotFoundError: if ``grid_path`` or ``gridmet_meta_path`` 
        are not found. If ``gridmet_meta_path`` is not passed as a 
        command line argument it is not in the current working directory
        and named "gridmet_cell_data.csv".    
        
    Note:
        If cells in the existing fishnet grid lie outside of the
        gridMET master fishnet, the GRIDMET_ID will be assigned
        the value of -999. Also, if the existing fishnet polygons
        do not exactly align with the gridMET dataset cells this function 
        will assign -999 for GRIDMET_ID.
    
    """
    
    # check if paths to input files were given
    if not os.path.isfile(grid_path):
        raise FileNotFoundError('The file path for the gridMET fishnet '\
                               +'was invalid or does not exist. ')
    tmp_out = grid_path.replace('.shp', '_tmp.shp')
    if not gridmet_meta_path:
        gridmet_meta_path = 'gridmet_cell_data.csv'
    if not os.path.isfile(gridmet_meta_path):
        raise FileNotFoundError('GridMET file path was not given and '\
                +'gridmet_cell_data.csv was not found in the current '\
                +'directory. Please assign the correct path or put '\
                +'gridmet_cell_data.csv in the current directory.\n')
    
    # load gridMET metadata file for looking up gridMET IDs
    gridmet_meta_df = pd.read_csv(gridmet_meta_path)
    # WGS 84 projection
    crs = from_epsg(4326) 

    # overwrite fishnet grid with updated GRIDMET_ID field
    with fiona.open(grid_path, 'r') as source:
        print(
            'Adding gridMET IDs to fishnet grid, saving to: \n',
             os.path.abspath(grid_path), '\n'
        )
        
        n_cells = len([f for f in source])
        print(
            'Looking up and assigning values for ', n_cells, 
            ' gridcells.\n'
        )        
        
        # Copy the source schema and add GRIDMET_ID property.
        sink_schema = source.schema
        sink_schema['properties']['GRIDMET_ID'] = 'int'
        # overwrite file add spatial reference
        with fiona.open(
                tmp_out, 
                'w', 
                crs=crs, 
                driver=source.driver, 
                schema=sink_schema
            ) as sink:
            # add GRIDMET_ID feature to outfile
            for feature in source:
                coords = feature['geometry']['coordinates'][0]
                gridmet_id = get_cell_ID(coords, gridmet_meta_df)
                feature['properties']['GRIDMET_ID'] = gridmet_id
                sink.write(feature)
    # cannot open same file and write to it on Windows, overwrite temp
    root_dir = os.path.split(grid_path)[0]
    for f in os.listdir(root_dir):
        if '_tmp' in f:
            move(OPJ(root_dir, f), OPJ(root_dir, f.replace('_tmp', '')))
    print(
        'Completed assigning gridMET IDs to fishnet. \n'
    )
	
def build_subgrid(in_path, gridmet_meta_path=None, buffer=25):
    """
    Condenses workflow to create fishnet grid of gridMET cells 
    around climate stations utilizing multiple functions in 
    :mod:`spatial`. Namely, :func:`get_subgrid_bounds`, :func:`make_grid`, 
    and :func:`update_subgrid`.
    
    Arguments:
        in_path (str): path to summary_comp.csv file containing monthly
            bias ratios, lat, long, and other data. Created by 
            :mod:`etr_biascorrect.calc_bias_ratios.py`.
            
    Keyword Arguments:
        gridmet_meta_path (str): default None. Path to metadata CSV file 
            that contains all gridMET cells for the contiguous United 
            States. If None it is looked for at 
            ``etr-biascorrect/gridmet_cell_data.csv``.
        buffer (int): default 25. Number of gridMET cells to expand 
            the rectangular extent of the subgrid fishnet and interpolated
            output raster.
            
    Returns:
        None

    Examples:
        To build the fishnet with a 25 cell outer buffer and update with
        gridMET ID values and projection in one step,
        
        >>> from etr_biascorrect import spatial
        >>> # assign input paths
        >>> summary_file = 'monthly_ratios/summary_comp.csv'
        >>> gridmet_metadata_path = 'gridmet_cell_data.csv'
        >>> # make fishnet of gridMET cells for interpolation
        >>> spatial.build_subgrid(
                summary_file,
                gridmet_meta_path=gridmet_metadata_path
            )
            
        For more examples and details on functions that ``build_subgrid``
        utilizes see :func:`get_subgrid_bounds`, :func:`make_grid`, 
        and :func:`update_subgrid`.
            
    Raises:
        FileNotFoundError: if input summary CSV file is not found. 
        see also :func:`update_subgrid` for error with ``gridmet_meta_path``.
        
    Note:
        The fishnet grid is saved to a subdirectory named "spatial" in
        the directory where the input summary file is found. Also, this 
        function performs a workflow comprised of multiple functions which
        can be used independently if for example the the fishnet was created
        already and the user only want's to update it to add gridMET ID values
        using :func:`update_subgrid`.
        
    """
    if not os.path.isfile(in_path):
        raise FileNotFoundError('Input summary CSV file given'+\
                                ' was invalid or not found')
    # save grid to "spatial" subdirectory of in_path
    path_root = os.path.split(in_path)[0]
    out_dir = OPJ(path_root, 'spatial')
    grid_file = OPJ(out_dir, 'grid.shp')
    
    # create output directory if does not exist
    if not os.path.isdir(out_dir):
        print(
            os.path.abspath(out_dir), 
            ' does not exist, creating directory.\n'
        )
        os.mkdir(out_dir)
        
    # get bounds based on station locations in CSV
    bnds = get_subgrid_bounds(in_path, buffer=buffer)
    # make grid
    make_grid(grid_file, bnds) 
    # add gridMET_IDs and WGS 84 projection
    update_subgrid(grid_file, gridmet_meta_path)
    
# add to argparse "method" for rbf, scale_factor

def interpolate(in_path, var_name, scale_factor=0.1, 
                function='inverse', gridmet_meta_path=None, buffer=25,
                zonal_stats=True):
    """
    Use a radial basis function to interpolate 2-dimensional surface of
    calculated bias ratios for climate stations in input summary CSV. 
    Options allow for modifying (down/up scaling) the resolution of the 
    resampling grid and to select from multiple radial basis functions. 
    Interploated surface is saved as a GeoTIFF raster.
    
    Arguments:
        in_path (str): path to summary_comp.csv file containing monthly
            bias ratios, lat, long, and other data. Created by 
            :mod:`etr_biascorrect.calc_bias_ratios.py`.
        var_name (str): name of variable in ``in_path`` to conduct 2-D
            interpolation. e.g. "Annual_mean" that is found in summary
            input CSV (``in_path``).
    
    Keyword Arguments:
        scale_factor (float, int): default 0.1. Scaling factor to apply to 
            original gridMET fishnet to create resampling fishnet. For example,
            if scale_factor = 0.1, the resolution will be one tenth of 
            the original gridMET resolution resulting in a 400 m resolution.
        function (str): default 'inverse', radial basis function to use 
            for interpolation. Options include 
            * 'multiquadric'
            * 'inverse'
            * 'gaussian'
            * 'linear'
            * 'cubic'
            * 'quintic'
            * 'thin_plate'
        gridmet_meta_path (str): default None. Path to metadata CSV file 
            that contains all gridMET cells for the contiguous United 
            States. If None it is looked for at 
            ``etr-biascorrect/gridmet_cell_data.csv``. 
        buffer (int): default 25. Number of gridMET cells to expand 
            the rectangular extent of the subgrid fishnet and interpolated
            output raster.
        zonal_stats (bool): default True. Calculate zonal means of interpolated
            surface to gridMET cells in fishnet and save to a CSV file. 
            The CSV file will be saved to the same directory of ``in_path``
            and will be named as "gridmet_summary_[function]_[res]m.csv"
            where [function] is the radial basis interpolation function and
            [res] is the resolution of the interpolated surface in meters.
            Also see :func:`gridmet_zonal_stats`.
            
    Returns:
        None
        
    Examples:
        Let's say we wanted to interpolate the "Annual_mean" bias 
        ratio in an input CSV first created by :mod:`calc_bias_ratios.py`.
        This example uses the gaussian radial basis function to interpolate 
        to a 400 m resolution surface. The result is a GeoTIFF raster that 
        has an extent that encompasses station locations in the input file 
        plus an additional optional buffer of outer gridMET cells.
        
        >>> from etr_biascorrect import spatial
        >>> # assign input paths
        >>> summary_file = 'monthly_ratios/summary_comp.csv'
        >>> gridmet_metadata_path = 'gridmet_cell_data.csv'
        >>> buffer = 10
        >>> var_name = 'Annual_mean'
        >>> func = 'gaussian'
        >>> interpolate(summary_file, var_name, scale_factor=0.1, 
                    function=func, gridmet_meta_path=gridmet_metadata_path, 
                    buffer=buffer)
                    
        The raster will be saved to "monthly_ratios/spatial/Annual_mean_400m.tiff"
        where the file name is based on the variable being interpolated 
        followed by the raster resolution. In the case the original gridMET
        resolution is 4 km therefore the scale facter of 0.1 results in 
        a 400 m resolution. 
        
        By using this function within Python it allows for interpolation 
        and calculation of zonal statistics of bias ratios that
        are not part of the command line workflow. For example if we wanted to 
        use the spatial surface of the *median* July bias ratio as opposed to 
        the *mean* (the command line usage calculates utilizes only select 
        *mean* values of bias ratios),
        
        >>> var_name = 'Jul_median'
        >>> func = 'gaussian'
        >>> interpolate(summary_file, var_name, scale_factor=0.1, 
                    function=func, gridmet_meta_path=gridmet_metadata_path, 
                    buffer=buffer)
                    
        This will create the GeoTIFF raster::
        
            'monthly_ratios/spatial/Jul_median_400m.tiff'
        
        Now to calculate zonal statistics for gridMET cell,
        
        >>> raster_file = 'monthly_ratios/spatial/Jul_median_400m.tiff'
        >>> spatial.gridmet_zonal_stats(
                summary_file, 
                raster_file, 
                function=func,
                res=400
            )
        
        The final result will be the creation of::
            
            'monthly_ratios/gridmet_summary_gaussian_400m.csv'
            
        In "gridmet_summary_gaussian_400m.csv" the zonal mean for median 
        ratios of July station to gridMET etr will be stored along with
        gridMET IDs, e.g.::
        
            ==========  =================
            GRIDMET_ID  Jul_median
            ==========  =================
            515902      0.87439453125
            514516      0.888170013427734
            513130      0.90002197265625
            ...         ...
            ==========  =================
        
    
    Raises:
        FileNotFoundError: if input summary CSV file is not found. 
        see also :func:`update_subgrid` for error with ``gridmet_meta_path``.
        
    Note:
        This function can be used independently of :func:`build_subgrid`
        however, if the buffer and input summary_comp.csv files arguments
        differ from those used for :func:`interpolate` the raster may not
        fully cover the fishnet which may later cause errors or gaps
        in the final zonal statistics to gridMET cells.
        
    """
    if not os.path.isfile(in_path):
        raise FileNotFoundError('Input summary CSV file given'+\
                                ' was invalid or not found')
    if not gridmet_meta_path:
        gridmet_meta_path = 'gridmet_cell_data.csv'
    if not os.path.isfile(gridmet_meta_path):
        raise FileNotFoundError('GridMET file path was not given and '\
                +'gridmet_cell_data.csv was not found in the current '\
                +'directory. Please assign the correct path or put '\
                +'gridmet_cell_data.csv in the current directory.\n')
    # calc raster resolution in meters (as frac of 4 km)
    res = int(4 * scale_factor * 1000)
    print(
        'Interpolating point bias ratios for: {}\n'.format(var_name),
        'Using the "{}" radial basis function\n'.format(function),
        'Resolution (pixel size) of output raster: {} m'.format(res)
    )
    # path to save raster of interpolated grid scaled by scale_factor
    path_root = os.path.split(in_path)[0]

    out_dir = OPJ(path_root, 'spatial')
    out_file = OPJ(out_dir, '{name}_{res}m.tiff'.format(name=var_name, res=res))
    # create output directory if does not exist
    if not os.path.isdir(out_dir):
        print(
            os.path.abspath(out_dir), 
            ' does not exist, creating directory.\n'
        )
        os.mkdir(out_dir)
    
    print(            
        'GeoTIFF raster will be saved to: \n',
        os.path.abspath(out_file),
        '\n'
    )
    
    in_df = pd.read_csv(in_path)

    lon_pts, lat_pts = in_df.STATION_LON.values, in_df.STATION_LAT.values
    values = in_df[var_name].values
    

    # get bounding area from summary CSV
    bnds = get_subgrid_bounds(in_path, buffer=buffer)
    lon_min, lon_max, lat_min, lat_max = bnds
    
    nx_cells = int(np.round(np.abs((lon_min - lon_max) / CELL_SIZE)))
    ny_cells = int(np.round(np.abs((lat_min - lat_max) / CELL_SIZE)))
    # rbf requires uniform grid (n X n) so 
    # extend short dimension and clip later 
    nx_cells_out = copy.copy(nx_cells)
    ny_cells_out = copy.copy(ny_cells)
    # gdal requires "upper left" corner coordinates
    lat_max_out = copy.copy(lat_max)
    lon_max_out = copy.copy(lon_max)
    
    if not nx_cells == ny_cells:
        diff = np.abs(nx_cells - ny_cells)
        if nx_cells > ny_cells:
            lat_max += diff * CELL_SIZE
            ny_cells += diff
        else:
            lon_max += diff * CELL_SIZE
            nx_cells += diff

            
    # make finer/coarse grid by scale factor
    lons = np.linspace(lon_min, lon_max, int(np.round(nx_cells/scale_factor))+1)
    lats = np.linspace(lat_min, lat_max, int(np.round(ny_cells/scale_factor))+1)
    # extent for original created by spatial.build_subgrid
    # add one to make sure raster covers full extent
    lons_out = np.linspace(lon_min, lon_max_out, int(np.round(nx_cells_out/scale_factor))+1)
    lats_out = np.linspace(lat_min, lat_max_out, int(np.round(ny_cells_out/scale_factor))+1)

    XI, YI = np.meshgrid(lons, lats)
    # apply rbf interpolation
    rbf = Rbf(lon_pts, lat_pts, values, function=function)
    ZI = rbf(XI, YI)
    # clip to original extent, rbf array flips axes, and row order... 
    ZI_out = ZI[0:len(lats_out),0:len(lons_out)]
    ZI_out = np.flip(ZI_out,axis=0)
    #### save interpolated data as raster 
    pixel_size = CELL_SIZE * scale_factor
    # number of pixels in each direction
    x_size = len(lons_out)
    y_size = len(lats_out)
    # set geotransform info
    gt = [
            lon_min,
            pixel_size,
            0,
            lat_max_out,
            0,
            -pixel_size
    ]
    # make geotiff raster
    driver = gdal.GetDriverByName('GTiff')
    ds = driver.Create(
        out_file,
        x_size, 
        y_size, 
        1, 
        gdal.GDT_Float32, 
    )

    # set projection geographic lat/lon WGS 84
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    ds.SetProjection(srs.ExportToWkt())
    # assign spatial dimensions 
    ds.SetGeoTransform(gt)
    outband = ds.GetRasterBand(1)
    # save rbf interpolated array as geotiff raster close
    outband.WriteArray(ZI_out)
    ds = None
    
    # calculate zonal statistics save means for each gridMET cell
    if zonal_stats:
        gridmet_zonal_stats(in_path, out_file, function=function, res=res)

def gridmet_zonal_stats(in_path, raster, function=None, res=None):
    """
    Calculate zonal means from interpolated surface of etr bias ratios
    created by :func:`interpolate` using the fishnet grid created by 
    :func:`build_subgrid`. Save mean values for each gridMET cell to
    a CSV file joined to gridMET IDs. 
    
    Arguments:
        in_path (str): path to summary_comp.csv file containing monthly
            bias ratios, lat, long, and other data. Created by 
            :mod:`etr_biascorrect.calc_bias_ratios.py`. The path is used
            to locate the fishnet shapefile that should be under "spatial"
            subdirectory. 
        raster (str): path to interpolated raster of bias ratios to
            be used for zonal stats. First created by :func:`interpolate`.
            
    Keyword Arguments:
        function (str): default None. Used only to rename output summary
            file.
        res (int): default None. Used only to rename output summary file.
        
    Example:
        Although it is prefered to use this function as part of 
        :func:`interpolate` or indirectly using the :mod:`spatial.py`
        command line usage. However if the grid shapefile and spatial
        interpolated raster(s) have already been created,
        
        >>> from etr_biascorrect import spatial
        >>> # assign input paths
        >>> summary_file = 'monthly_ratios/summary_comp.csv'        
        >>> raster_file = 'monthly_ratios/spatial/Jul_med_400m.tiff'
        >>> spatial.gridmet_zonal_stats(
                summary_file, 
                raster_file, 
                function=func,
                res=400
            )
        
        The final result will be the creation of::
            
            'monthly_ratios/gridmet_summary_gaussian_400m.csv'
            
        The resulting CSV contains the gridMET IDS and zonal means
        for each gridMET cell in the fishnet which must exist at::
        
            'monthly_ratios/spatial/grid.shp'
            
        also see :func:`interpolate`
        
    Raises:
        FileNotFoundError: if the input summary CSV file or the 
            fishnet for extracting zonal statistics do not exist.
            The fishnet should be in the subdirectory of ``in_path``
            at "/spatial/grid.shp".

    Note:
        The final CSV file with gridMET zonal mean etr bias ratios
        will *not* be overwritten after it is first created if the
        same interpolation method and grid resolution is used subsequently.

    """
    if not os.path.isfile(in_path):
        raise FileNotFoundError('Input summary CSV file given'+\
                                ' was invalid or not found')
    # look for fishnet created in 'in_path/spatial'
    path_root = os.path.split(in_path)[0]
    out_dir = OPJ(path_root, 'spatial')
    grid_file = OPJ(out_dir, 'grid.shp')
    out_file = OPJ(path_root, 'gridmet_summary.csv')
    if function and res:
        out_file = out_file.replace(
            '.csv', '_{func}_{res}m.csv'.format(func=function,res=res)
        )
    # this error would only occur when using within Python
    if not os.path.isfile(grid_file):
        raise FileNotFoundError(
            os.path.abspath(grid_file),
            '\ndoes not exist, create it using spatial.build_subgrid first'
        )
    # get var name from raster file (always before resolution from interpolate)
    var_match = re.compile(r'(.*)_(\d+)m.tiff')
    raster_name = os.path.split(raster)[1]
    var_name = var_match.match(raster_name).group(1)
    res = var_match.match(raster_name).group(2)
    print(
        'calculating gridMET zonal means for',
        var_name, 'from', res, 'm resolution raster'
    )
    # calc zonal stats and open for grid IDs
    with fiona.open(grid_file, 'r') as source:
        zs = zonal_stats(source, raster)
        gridmet_ids = [f['properties'].get('GRIDMET_ID') for f in source]

    # get just mean values, zonal_stats also calcs max, min, pixel count
    means = [z['mean'] for z in zs]
    out_df = pd.DataFrame(
        data={
            'GRIDMET_ID': gridmet_ids, 
            var_name: means
        }
    )
    out_df.GRIDMET_ID = out_df.GRIDMET_ID.astype(int)
    # drop rows for cells outside of gridMET master grid
    out_df = out_df.drop(out_df[out_df.GRIDMET_ID == -999].index)

    # save or update existing csv file
    if not os.path.isfile(out_file):
        print(
            os.path.abspath(out_file),
            '\ndoes not exist, creating file'
        )
        out_df.to_csv(out_file, index=False)
    else:
        existing_df = pd.read_csv(out_file)
        existing_df.GRIDMET_ID = existing_df.GRIDMET_ID.astype(int)
        # avoid adding same column twice and duplicates
        if var_name in existing_df.columns:
            pass
        else:
            existing_df = existing_df.merge(out_df, on='GRIDMET_ID')
            #existing_df = pd.concat([existing_df, out_df], axis=1).drop_duplicates()
            existing_df.to_csv(out_file, index=False)   
        
def arg_parse():
    """
    Command line usage of spatial.py for creating shapefile of all stations
    found in comprehensive summary CSV file, build fishnet around stations
    and perform spatial interpolation for gridMET cells.
    """
    parser = argparse.ArgumentParser(
        description=arg_parse.__doc__,
        #formatter_class=argparse.RawDescriptionHelpFormatter)
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    optional = parser._action_groups.pop() # optionals listed second
    required = parser.add_argument_group('required arguments')
    required.add_argument(
        '-i', '--input', metavar='PATH', required=True,
        help='Input summary_comp CSV file of climate/gridMET bias ratios '+\
             'created by first running calc_bias_ratios.py')
    optional.add_argument(
        '-g', '--gridmet', metavar='PATH', required=False,
        help='GridMET master CSV file with cell data, packaged with '+\
             'etr-biascorrect at etr-biascorrect/gridmet_cell_data.csv '+\
             'if not given it needs to be located in the currect directory')
    optional.add_argument(
        '-b', '--buffer', required=False, default=25, type=int, metavar='',
        help='Number of gridMET cells to expand outer bounds of fishnet '+\
             'which can be used for extrapolation')
    optional.add_argument(
        '-s', '--scale', required=False, default=0.1, type=float, metavar='',
        help='Scale facter used on gridMET resolution to down/upscale '+\
             'interpolation output which is applied to the gridMET cell '+\
             'size of 4 km, default 400 m')
    optional.add_argument(
        '-f', '--function', required=False, default='inverse', type=str,
        metavar='', help='Radial basis function to use for interpolation '+\
                'Options include: multiquadric, inverse, gaussian, linear '+\
                'cubic, quintic, and thin_plate')
#    optional.add_argument(
#        '--debug', default=logging.INFO, const=logging.DEBUG,
#        help='Debug level logging', action="store_const", dest="loglevel")
    parser._action_groups.append(optional)# to avoid optionals listed first
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = arg_parse()

    main(
        input_file_path=args.input, 
        gridmet_meta_path=args.gridmet, 
        buffer=args.buffer,
        scale_factor=args.scale,
        function=args.function
    )
