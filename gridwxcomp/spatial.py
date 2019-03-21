# -*- coding: utf-8 -*-
"""
Perform multiple workflows needed to estimate the spatial surface of monthly 
and annual bias ratios for multiple climatic variables. The input file is first 
created by :mod:`gridwxcomp.calc_bias_ratios`.

Attributes:
    CELL_SIZE (float): constant gridMET cell size in decimal degrees,
        value = 0.041666666666666664.

Note:
    All spatial files, i.e. vector and raster files, utilize the
    *ESRI Shapefile* file or GeoTiff format. 
    
Todo:
    * logging 

"""

import os
import re
import argparse
import copy
import pkg_resources
from math import ceil, pow, sqrt
from pathlib import Path
from shutil import move

import fiona
import numpy as np
import pandas as pd
import rasterio
from scipy.interpolate import Rbf
from shapely.geometry import Point, Polygon, mapping
from fiona import collection
from fiona.crs import from_epsg
from osgeo import gdal, osr, ogr
from rasterstats import zonal_stats

# constant gridmet resolution in decimal degrees
CELL_SIZE = 0.041666666666666664

OPJ = os.path.join
   
def main(input_file_path, layer='all', out=None, buffer=25, scale_factor=0.1, 
         function='invdist', smooth=0, params=None, zonal_stats=True,
         overwrite=False, options=None, gridmet_meta_path=None):
    """
    Create point shapefile of monthly mean bias ratios from comprehensive
    CSV file created by :mod:`gridwxcomp.calc_bias_ratios`. Build fishnet grid 
    around the climate stations that matches the gridMET dataset. Perform 
    spatial interpolation to estimate 2-dimensional surface of bias ratios and 
    extract interpolated values to gridMET cells. 
    
    Arguments:
        input_file_path (str): path to [var]_summary_comp CSV file containing 
            monthly bias ratios, lat, long, and other data. Shapefile 
            "[var]_summary_pts.shp" is saved to parent directory of 
            ``input_file_path`` under a subdirectory named "spatial".

    Keyword Arguments:
        out (str): default None. Subdirectory to save GeoTIFF rasters.
        gridmet_meta_path (str): default None. Path to metadata CSV file 
            that contains all gridMET cell info. By default it is looked
            for in the package installation directory or current directory
            as "gridmet_cell_data.csv".
        buffer (int): number of gridMET cells to expand the rectangular extent
            of the subgrid fishnet (default=25). 
        scale_factor (float, int): scaling factor to apply to original
            gridMET fishnet to create resampling fishnet. For example,
            if scale_factor = 0.1, the resolution will be one tenth of 
            the original girdMET resolution resulting in a 400 m resolution.
        function (str): default 'invdist'. Interpolation method, gdal 
            methods include 
            * 'invdist' 
            * 'indistnn' 
            * 'linear' 
            * 'average'
            * 'nearest'
            see `gdal_grid <https://www.gdal.org/gdal_grid.html>`_.
            Radial basis functions, see :class:`scipy.interpolate.Rbf`, 
            include
            * 'multiquadric'
            * 'inverse'
            * 'gaussian'
            * 'linear_rbf'
            * 'cubic'
            * 'quintic'
            * 'thin_plate'
        smooth (float): default 0. Smooth parameter for Rbf functions.
        params (dict, str, or None): default None. Parameters for interpolation
            using gdal, see defaults in :class:`gridwxcomp.InterpGdal`.
        overwrite (bool): default False. If True overwrite the grid 
            shapefile that already exists.
        zonal_stats (bool): default True. Calculate zonal means of interpolated
            surface to gridMET cells in fishnet and save to a CSV file. 
            The CSV file will be saved to the same directory as the interpolated
            raster file(s).
        options (str or None): default None. Extra command line arguments
            for gdal interpolation.
        gridmet_meta_path (str): default None. Path to metadata CSV file 
            that contains all gridMET cells for the contiguous United 
            States. If None it is looked for at the install directory of 
            gridwxcomp (i.e. with pip install) or within the current directory
            as "gridmet_cell_data.csv".

    Returns:
        None

    Examples:
        From the command line,

        .. code-block:: sh

            $ python spatial.py -i monthly_ratios/etr_mm_summary_comp.csv 

        This will produce a subdirectory under "monthly_ratios" named
        "spatial" where a point shapefile will be saved as 
        "etr_mm_summary_pts.shp". It will also create a spatially referenced 
        fishnet of gridMET cells (at "monthly_ratios/spatial/grid.shp") with 
        a buffer of 25 gridMET cells added around the encompassed climate 
        stations. Next a 2-dimensional surface is interpolated from the point 
        data of each mean bias ratio in etr_mm_summary_comp.csv which includes 
        monthly as well as growing season and annual means using the inverse 
        distance weighting method. The interpolated rasters are saved to::
        
            'monthly_ratios/spatial/etr_mm_inverse_dist_400m/'' 
            
        with file names of the form::
        
            [time_period].tiff
            
        where [time_period] refers to the bias ratio time aggregate e.g. 
        Apr_mean or Annual_mean. Resolution in meters of the spatially 
        interpolated raster, gridMET variable name, and interpolation function 
        are saved to the output directory which can have additional parent 
        subdirectories if the ``out`` argument is given. Zonal mean statistics 
        are calculated for gridMET cells in the fishnet that was created. They 
        are saved to::
        
            'monthly_ratios/spatial/etr_mm_inverse_dist_400m/gridMET_stats.csv'
            
        The contents of the CSV file will look something like::
        
            ==========  ========  ========  ========  ========
            GRIDMET_ID  Apr_mean  Aug_mean  Dec_mean  ...
            ==========  ========  ========  ========  ========
            515902      1.028707  0.856440  1.058291  ...
            514516      1.035543  0.862066  1.065963  ...
            ...         ...       ...       ...       ...
            ==========  ========  ========  ========  ========
        
        Alternatively, to use a smaller buffer zone on the fishnet grid 
        used for interpolation if, for example, extrapolation is not needed, 
        use the ``[-b, --buffer]`` option

        .. code-block:: sh

            $ python spatial.py -i monthly_ratios/etr_mm_summary_comp.csv -b 5

        Or, if we wanted to interpolate to a 200 m resolution 
        (i.e. scale_factor = 0.05, 0.05 x 4 km = 200 m) using the 'inverse'
        radial basis function,
        
        .. code-block:: sh

            $ python spatial.py -i monthly_ratios/etr_mm_summary_comp.csv -s 0.05 -f 'inverse'        

        To calculate zonal means for a different climate variable, e.g. 
        observed ET ("eto_mm"), as opposed to reference ET (default) use the 
        corresponding input file
        
        .. code-block:: sh

            $ python spatial.py -i monthly_ratios/eto_mm_summary_comp.csv -s 0.05 -f 'inverse'                

        In this case the final zonal statistics of zonal mean bias ratios 
        will be saved to::
        
            'monthly_ratios/spatial/etr_mm_inverse_200m/gridMET_stats.csv'

        To run interpolation of a single layer that is not part of the default
        mean layers, e.g. the coefficient of variation of the growing season
        bias ratio, assign the ``[-l, --layer]`` option,
        
        .. code-block:: sh

            $ python spatial.py -i monthly_ratios/etr_mm_summary_comp.csv -l April_to_oct_cv
            
        Also see examples of :func:`make_points_file`, :func:`make_grid`, 
        :func:`interpolate`, and the :class:`gridwxcomp.InterGdal` class.

    """
    # build fishnet for interpolation
    make_grid(input_file_path, 
              gridmet_meta_path=gridmet_meta_path, 
              buffer=buffer, 
              overwrite=overwrite)
    
    # run spatial interpolation depending on options
    interpolate(
        input_file_path, 
        layer=layer, 
        out=out,
        scale_factor=scale_factor, 
        function=function, 
        smooth=smooth,
        params=params,
        buffer=buffer,
        zonal_stats=zonal_stats,
        options=options,
        gridmet_meta_path=gridmet_meta_path) 

def make_points_file(in_path):
    """
    Create vector shapefile of points with monthly mean bias ratios 
    for climate stations using all stations found in a comprehensive
    CSV file created by :mod:`gridwxcomp.calc_bias_ratios`.
    
    Arguments:
        in_path (str): path to [var]_summary_comp.CSV file containing 
            monthly bias ratios, lat, long, and other data. Shapefile 
            "[var]_summary_pts.shp" is saved to parent directory of 
            ``in_path`` under "spatial" subdirectory.
            
    Returns:
        None
    
    Example:
        Create shapefile containing point data for all climate stations
        in input summary file created by :mod:`gridwxcomp.calc_bias_ratios`
        
        >>> from gridwxcomp import spatial
        >>> # path to comprehensive summary CSV
        >>> summary_file = 'monthly_ratios/etr_mm_summary_comp.csv'
        >>> spatial.make_points_file(summary_file)
        
        The result is the file "etr_mm_summary_pts.shp" being saved to 
        a subdirectory created in the directory containing ``in_path``
        named "spatial", i.e.::

            "monthly_ratios/spatial/etr_mm_summary_pts.shp". 
        
    Raises:
        FileNotFoundError: if input summary CSV file is not found.
        
    Note:
        :func:`make_points_file` will overwrite any existing point
        shapefile of the same climate variable. 
        
    """
    if not os.path.isfile(in_path):
        raise FileNotFoundError('Input summary CSV file: '+\
                                '{}\nwas not found.'.format(in_path))
    print(
        '\nMapping point data for climate stations in: \n',
        in_path, '\n'
    )
    in_df = pd.read_csv(in_path, index_col='STATION_ID', na_values=[-999])
    # save shapefile to "spatial" subdirectory of in_path
    path_root = os.path.split(in_path)[0]
    file_name = os.path.split(in_path)[1]
    # get variable name from input file prefix
    var_name = file_name.split('_summ')[0]
    out_dir = OPJ(path_root, 'spatial')
    out_file = OPJ(out_dir, '{v}_summary_pts.shp'.format(v=var_name))
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
            'summer': 'float',
            'growseason': 'float',
            'annual': 'float',
            'Jan_cnt': 'float',
            'Feb_cnt': 'float',
            'Mar_cnt': 'float',
            'Apr_cnt': 'float',
            'May_cnt': 'float',
            'Jun_cnt': 'float',
            'Jul_cnt': 'float',
            'Aug_cnt': 'float',
            'Sep_cnt': 'float',
            'Oct_cnt': 'float',
            'Nov_cnt': 'float',
            'Dec_cnt': 'float',
            'summer_cnt': 'float',
            'grow_cnt': 'float',
            'annual_cnt': 'float',
            'Jan_std': 'float',
            'Feb_std': 'float',
            'Mar_std': 'float',
            'Apr_std': 'float',
            'May_std': 'float',
            'Jun_std': 'float',
            'Jul_std': 'float',
            'Aug_std': 'float',
            'Sep_std': 'float',
            'Oct_std': 'float',
            'Nov_std': 'float',
            'Dec_std': 'float',
            'summer_std': 'float',
            'grow_std': 'float',
            'annual_std': 'float',
            'Jan_cv': 'float',
            'Feb_cv': 'float',
            'Mar_cv': 'float',
            'Apr_cv': 'float',
            'May_cv': 'float',
            'Jun_cv': 'float',
            'Jul_cv': 'float',
            'Aug_cv': 'float',
            'Sep_cv': 'float',
            'Oct_cv': 'float',
            'Nov_cv': 'float',
            'Dec_cv': 'float',
            'summer_cv': 'float',
            'grow_cv': 'float',
            'annual_cv': 'float',
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
                    'summer': row['summer_mean'],
                    'growseason': row['growseason_mean'],
                    'annual': row['annual_mean'],
                    'Jan_cnt': row['Jan_count'],
                    'Feb_cnt': row['Feb_count'],
                    'Mar_cnt': row['Mar_count'],
                    'Apr_cnt': row['Apr_count'],
                    'May_cnt': row['May_count'],
                    'Jun_cnt': row['Jun_count'],
                    'Jul_cnt': row['Jul_count'],
                    'Aug_cnt': row['Aug_count'],
                    'Sep_cnt': row['Sep_count'],
                    'Oct_cnt': row['Oct_count'],
                    'Nov_cnt': row['Nov_count'],
                    'Dec_cnt': row['Dec_count'],
                    'summer_cnt': row['summer_count'],
                    'grow_cnt': row['growseason_count'],
                    'annual_cnt': row['annual_count'],
                    'Jan_std': row['Jan_stdev'],
                    'Feb_std': row['Feb_stdev'],
                    'Mar_std': row['Mar_stdev'],
                    'Apr_std': row['Apr_stdev'],
                    'May_std': row['May_stdev'],
                    'Jun_std': row['Jun_stdev'],
                    'Jul_std': row['Jul_stdev'],
                    'Aug_std': row['Aug_stdev'],
                    'Sep_std': row['Sep_stdev'],
                    'Oct_std': row['Oct_stdev'],
                    'Nov_std': row['Nov_stdev'],
                    'Dec_std': row['Dec_stdev'],
                    'summer_std': row['summer_stdev'],
                    'grow_std': row['growseason_stdev'],
                    'annual_std': row['annual_stdev'],
                    'Jan_cv': row['Jan_cv'],
                    'Feb_cv': row['Feb_cv'],
                    'Mar_cv': row['Mar_cv'],
                    'Apr_cv': row['Apr_cv'],
                    'May_cv': row['May_cv'],
                    'Jun_cv': row['Jun_cv'],
                    'Jul_cv': row['Jul_cv'],
                    'Aug_cv': row['Aug_cv'],
                    'Sep_cv': row['Sep_cv'],
                    'Oct_cv': row['Oct_cv'],
                    'Nov_cv': row['Nov_cv'],
                    'Dec_cv': row['Dec_cv'],
                    'summer_cv': row['summer_cv'],
                    'grow_cv': row['growseason_cv'],
                    'annual_cv': row['annual_cv'],
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
            monthly bias ratios, created by :func:`gridwxcomp.calc_bias_ratios`.
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

def make_grid(in_path, bounds=None, buffer=25, overwrite=False, 
        gridmet_meta_path=None):
    """
    Make fishnet grid (vector file of polygon geometry) for 
    select gridMET cells based on bounding coordinates. 
    
    Add gridMET ID values to each cell based on their centroid 
    lookup in ``gridwxcomp/gridmet_cell_data.csv``. Assigns the 
    WGS84 reference coordinate system. The grid is later used to 
    spatially interpolate point data. Modified from the 
    `Python GDAL/OGR Cookbook <https://pcjericks.github.io/py-gdalogr-cookbook/vector_layers.html#create-fishnet-grid>`_.
    
    Arguments:
        in_path (str): path to [var]_summary_comp.csv file containing monthly
            bias ratios, lat, long, and other data. Created by 
            :func:`gridwxcomp.calc_bias_ratios`.
    
    Keyword Arguments:
        bounds (tuple or None): default None. Tuple of bounding coordinates 
            in the following order (min long, max long, min lat, max lat) 
            which need to be in decimal degrees. Need to align with gridMET
            resolution outer corners. If None, get extent from centoid
            locations of climate stations in ``in_path`` summary CSV. 
        buffer (int): default 25. Number of gridMET cells to expand 
            the rectangular extent of the subgrid fishnet and interpolated
            output raster.
        overwrite (bool): default False. If True overwrite the grid 
            shapefile at ``out_path`` if it already exists.
        gridmet_meta_path (str): default None. Path to metadata CSV file 
            that contains all gridMET cells for the contiguous United 
            States. If None it is looked for at the install directory of 
            gridwxcomp (i.e. with pip install) or within the current directory
            as "gridmet_cell_data.csv".

    Returns:
        None
        
    Examples:
        Build a fishnet uniform grid that matches gridMET cells around
        climate stations found in ``in_path`` summary CSV file with a 
        25 cell buffer. 
        
        >>> from gridwxcomp import spatial
        >>> # assign input paths
        >>> summary_file = 'monthly_ratios/etr_mm_summary_comp.csv'
        >>> # make fishnet of gridMET cells for interpolation
        >>> spatial.make_grid(summary_file)
            
        The file will be saved as "grid.shp" to a newly created subdirectory
        "spatial" in the same directory as the input summary CSV file. 
        To use a smaller buffer to the extent of the grid assign the 
        ``buffer`` keyword argument
        
        >>> spatial.make_grid(summary_file, buffer=5)        
            
        If the grid file already exists the default action is to not 
        overwrite. To overwrite an existing grid if, for example, you 
        are working with a new set of climate stations as input, then 
        set the ``overwrite`` keyword argument to True. 
        
        >>> spatial.make_grid(summary_file, overwrite=True, buffer=5)       
            
    Raises:
        FileNotFoundError: if input summary CSV file is not found. 
        
    Note:
        If cells in the fishnet grid lie outside of the gridMET master 
        fishnet for the contiguous U.S., the GRIDMET_ID will be assigned 
        to -999. 
        
    """
        
    if not os.path.isfile(in_path):
        raise FileNotFoundError('Input summary CSV file given'+\
                                ' was invalid or not found')
    # save grid to "spatial" subdirectory of in_path
    path_root = os.path.split(in_path)[0]
    out_dir = OPJ(path_root, 'spatial')
    out_path = OPJ(out_dir, 'grid.shp')
    
    # skip building grid if already exists
    if os.path.isfile(out_path) and not overwrite:
        print(
            '\nFishnet grid already exists at: \n',
            out_path,
            '\nnot overwriting.\n'
        )
        return
    # print message if overwriting existing grid
    elif os.path.isfile(out_path) and overwrite:
        print(
            '\nOverwriting fishnet grid at: \n',
            out_path,
            '\n'
        )
    # create output directory if does not exist
    if not os.path.isdir(out_dir):
        print(
            os.path.abspath(out_dir), 
            ' does not exist, creating directory.\n'
        )
        os.mkdir(out_dir)
        
    # get grid extent based on station locations in CSV
    if not bounds:
        bounds = get_subgrid_bounds(in_path, buffer=buffer)    
    xmin, xmax, ymin, ymax = bounds
    # read path and make parent directories if they don't exist    
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
    
    # get n rows
    rows = ceil((ymax-ymin) / CELL_SIZE)
    # get n columns
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
    # reopen grid and assign gridMET attribute and coord. ref.
    _update_subgrid(out_path, gridmet_meta_path=gridmet_meta_path)

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
            gridMET metadata CSV file "gridmet_cell_data.csv".
        
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

def _update_subgrid(grid_path, gridmet_meta_path=None):
    """
    Helper function to assign gridMET ID values and EPSG WGS 84 
    projection to fishnet grid.
    
    Arguments:
        grid_path (str): path to fishnet grid shapefile for subset of
            gridMET cells which contain climate stations.
            
    Keyword Arguments:
        gridmet_meta_path (str): default None. Path to metadata CSV file 
            that contains all gridMET cells for the contiguous United 
            States. If None it is looked for at the install directory of 
            gridwxcomp (i.e. with pip install) or within the current directory
            as "gridmet_cell_data.csv".

    Returns:
        None
        
    Raises:
        FileNotFoundError: if ``grid_path`` or ``gridmet_meta_path`` 
        are not found. If ``gridmet_meta_path`` is not passed as a 
        command line argument and is not found in the gridwxcomp install
        directory and it is not in the current working directory
        and named "gridmet_cell_data.csv".    
        
    """

    if not os.path.isfile(grid_path):
        raise FileNotFoundError('The file path for the gridMET fishnet '\
                               +'was invalid or does not exist. ')
    # look for pacakged gridmet_cell_data.csv if path not given
    if not gridmet_meta_path:
        try:
            if pkg_resources.resource_exists('gridwxcomp', 
                    "gridmet_cell_data.csv"):
                gridmet_meta_path = pkg_resources.resource_filename(
                    'gridwxcomp', 
                    "gridmet_cell_data.csv"
                    )
        except:
            gridmet_meta_path = 'gridmet_cell_data.csv'
    if not os.path.exists(gridmet_meta_path):
        raise FileNotFoundError('GridMET file path was not given and '+\
                'gridmet_cell_data.csv was not found in the gridwxcomp '+\
                'install directory. Please assign the path or put '+\
                '"gridmet_cell_data.csv" in the current working directory.\n')
      
    tmp_out = grid_path.replace('.shp', '_tmp.shp')

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

    
def interpolate(in_path, layer='all', out=None, scale_factor=0.1, 
                function='invdist', smooth=0, params=None, bounds=None, 
                buffer=25, zonal_stats=True, options=None, 
                gridmet_meta_path=None):
    """
    Use various methods to interpolate a 2-dimensional surface of
    calculated bias ratios or other statistics for station/gridMET
    pairs found in input summary CSV. 
    
    Options allow for modifying (down/up scaling) the resolution of the 
    resampling grid and to select from multiple interpolation methods.
    Interploated surfaces are saved as GeoTIFF rasters.
    
    Arguments:
        in_path (str): path to [var]_summary_comp.csv file containing 
            monthly bias ratios, lat, long, and other data. Created by 
            :func:`gridwxcomp.calc_bias_ratios`.

    Keyword Arguments:
        layer (str or list): default 'all'. Name of variable in ``in_path`` 
            to conduct 2-D interpolation. e.g. 'Annual_mean'.
        out (str): default None. Subdirectory to save GeoTIFF raster(s).
        scale_factor (float, int): default 0.1. Scaling factor to apply to 
            original gridMET fishnet to create resampling resolution. If 
            scale_factor = 0.1, the resolution will be one tenth gridMET 
            resolution or 400 m.
        function (str): default 'invdist'. Interpolation method, gdal 
            methods include: 'invdist', 'indistnn', 'linear', 'average',
            and 'nearest' see `gdal_grid <https://www.gdal.org/gdal_grid.html>`_.
            Radial basis functions, see :class:`scipy.interpolate.Rbf`, 
            include: 'multiquadric', 'inverse', 'gaussian', 'linear_rbf',
            'cubic', 'quintic', and 'thin_plate'.
        smooth (float): default 0. Smooth parameter for Rbf functions.
        params (dict, str, or None): default None. Parameters for interpolation
            using gdal, see defaults in :class:`gridwxcomp.InterpGdal`.
        bounds (tuple or None): default None. Tuple of bounding coordinates 
            in the following order (min long, max long, min lat, max lat) 
            which need to be in decimal degrees. Need to align with gridMET
            resolution outer corners. If None, get extent from centoid
            locations of climate stations in ``in_path`` summary CSV. 
        buffer (int): default 25. Number of gridMET cells to expand 
            the rectangular extent of the subgrid fishnet and interpolated
            output raster.
        zonal_stats (bool): default True. Calculate zonal means of interpolated
            surface to gridMET cells in fishnet and save to a CSV file. 
            The CSV file will be saved to the same directory as the interpolated
            raster file(s).
        options (str or None): default None. Extra command line arguments
            for gdal interpolation.
        gridmet_meta_path (str or None): default None. Path to metadata CSV 
            file that contains all gridMET cells for the contiguous United 
            States. If None it is looked for at the install directory of 
            gridwxcomp (e.g. after pip install gridwxcomp) or within the 
            current directory as 'gridmet_cell_data.csv'.

    Returns:
        None
        
    Examples:
        Let's say we wanted to interpolate the "Annual_mean" bias 
        ratio in an input CSV first created by :func:`gridwxcomp.calc_bias_ratios`.
        This example uses the "invdist" method (default) to interpolate 
        to a 400 m resolution surface. The result is a GeoTIFF raster that 
        has an extent that encompasses station locations in the input file 
        plus an additional optional buffer of outer gridMET cells.
        
        >>> from gridwxcomp import spatial
        >>> summary_file = 'monthly_ratios/etr_mm_summary_comp.csv'
        >>> buffer = 10
        >>> layer = 'Annual_mean'
        >>> params = {'power':1, 'smooth':20}
        >>> out_dir = 's20_p1' # optional subdir name for saving rasters
        >>> interpolate(summary_file, layer=layer, out=out_dir, 
        >>>     scale_factor=0.1, params=params, buffer=buffer)
                    
        The interpolated raster will be saved to::
        
            'monthly_ratios/spatial/etr_mm_invdist_400m/s20_p1/Annual_mean.tiff'
            
        where the file name and directory is based on the variable being 
        interpolated, methods, and the raster resolution. The ``out`` 
        keyword argument lets us add any number of subdirectories to the final 
        output directory, in this case the 's20_p1' dir contains info on params. 
        
        In this case the original gridMET resolution is 4 km therefore the 
        scale facter of 0.1 results in a 400 m resolution. 
        
        The final result will be the creation of the CSV::
            
            'monthly_ratios/spatial/etr_mm_invdist_400m/s20_p1/gridMET_stats.csv'
            
        In "gridMET_stats.csv" the zonal mean for ratios of annual station to 
        gridMET ETr will be stored along with gridMET IDs, e.g.
        
            ==========  =================
            GRIDMET_ID  Annual_mean
            ==========  =================
            515902      0.87439453125
            514516      0.888170013427734
            513130      0.90002197265625
            ...         ...
            ==========  =================      
            
        To calculate zonal statistics of bias ratios that are not part of 
        the default or command line workflow we can assign any numeric layer
        in the input summary CSV to be interpolations. For example if 
        we wanted to interpolate the coefficient of variation of the growing
        season bias ratio "April_to_oct_cv", then we could 
        estimate the surface of this variable straightforwardly,
        
        >>> layer = 'April_to_oct_cv'
        >>> func = 'invdistnn'
        >>> # we can also 'upscale' the interpolation resolution
        >>> interpolate(summary_file, layer=layer, scale_factor=2, 
        >>>     function=func, buffer=buffer)
                    
        This will create the GeoTIFF raster::
        
            'monthly_ratios/spatial/etr_mm_invdistnn_400m/April_to_oct_cv.tiff'
               
        And the zonal means will be added as a column named "April_to_oct_cv"
        to:: 
        
            'monthly_ratios/spatial/etr_mm_invdistnn_400m/gridMET_stats.csv'

        As with other components of ``gridwxcomp``, any other climatic
        variables that exist in the gridMET dataset can be used along
        with any corresponding station time series data from the user.
        The input (``in_path``) to all routines in :mod:`gridwxcomp.spatial` 
        is the summary CSV created by :func:`gridwxcomp.calc_bias_ratios`, the 
        prefix to this file is what determines the climatic variable 
        that spatial analysis is conducted. See :func:`gridwxcomp.calc_bias_ratios` 
        for examples of how to use different climatic variables, e.g. 
        TMax or ETo.
        
    Raises:
        FileNotFoundError: if the input summary CSV file or the 
            fishnet for extracting zonal statistics do not exist.
            The fishnet should be in the subdirectory of ``in_path``
            i.e. "<in_path>/spatial/grid.shp".
    Note:
        This function can be used independently of :func:`make_grid`
        however, if the buffer and input [var]_summary_comp.csv files 
        arguments differ from those used for :func:`interpolate` the 
        raster may not fully cover the fishnet which may result in 
        gaps in the zonal statistics.
        
    """
    # avoid circular import for InterpGdal for gdal interpolation methods
    from gridwxcomp.interpgdal import InterpGdal

    if not os.path.isfile(in_path):
        raise FileNotFoundError('Input summary CSV file given'+\
                                ' was invalid or not found')
    # look for packaged gridmet_cell_data.csv if path not given
    if not gridmet_meta_path:
        try:
            if pkg_resources.resource_exists('gridwxcomp', 
                    "gridmet_cell_data.csv"):
                gridmet_meta_path = pkg_resources.resource_filename(
                    'gridwxcomp', 
                    "gridmet_cell_data.csv"
                    )
        except:
            gridmet_meta_path = 'gridmet_cell_data.csv'
    if not os.path.exists(gridmet_meta_path):
        raise FileNotFoundError('GridMET file path was not given and '+\
                'gridmet_cell_data.csv was not found in the gridwxcomp '+\
                'install directory. Please assign the path or put '+\
                '"gridmet_cell_data.csv" in the current working directory.\n')
    # calc raster resolution in meters (as frac of 4 km)
    res = int(4 * scale_factor * 1000)
    # path to save raster of interpolated grid scaled by scale_factor
    path_root = os.path.split(in_path)[0]
    file_name = os.path.split(in_path)[1]
    # get variable name from input file prefix
    grid_var = file_name.split('_summ')[0]
    
    if not out: 
        out_dir = OPJ(path_root,
                      'spatial', 
                      '{}_{}_{}m'.format(grid_var, function, res))
    else:
        out_dir = OPJ(path_root,
                      'spatial',
                      '{}_{}_{}m'.format(grid_var, function, res),
                      out)
        
    # create output directory if does not exist
    if not os.path.isdir(out_dir):
        print(
            os.path.abspath(out_dir), 
            ' does not exist, creating directory.\n'
        )
        Path(out_dir).mkdir(parents=True, exist_ok=True)
    
    
    def _run_rbf_interpolation(layer, bounds, function, smooth):
        """Workflow for running scipy Rbf interpolation of scatter points"""
        out_file = OPJ(
            out_dir, 
            '{time_agg}.tiff'.format(time_agg=layer)
        )
        print(
            '\nInterpolating {g} point bias ratios for: {t}\n'.\
                format(g=grid_var, t=layer),
            'Using the "{}" method\n'.format(function),
            'Resolution (pixel size) of output raster: {} m'.format(res)
        )
        print(            
            'GeoTIFF raster will be saved to: \n',
            os.path.abspath(out_file)
        )
        # get grid extent based on station locations in CSV
        if not bounds:
            bounds = get_subgrid_bounds(in_path, buffer=buffer) 
        lon_min, lon_max, lat_min, lat_max = bounds
        # check if layer is in summary CSV 
        existing_layers = pd.read_csv(in_path).columns
        if not layer in existing_layers:
            print('column {} does not exist in input CSV:\n {}'.format(
               layer, in_path),
                 '\nSkipping interpolation.'
            )
            return
        
        # get point station data from summary CSV
        in_df = pd.read_csv(in_path, na_values=[-999])
        lon_pts, lat_pts = in_df.STATION_LON.values, in_df.STATION_LAT.values
        values = in_df[layer].values
    
        # mask out stations with missing data
        if in_df[layer].isnull().sum() > 0:
            mask = in_df[layer].notnull()
            n_missing = in_df[layer].isna().sum()
            # if one point or less data points exists exit
            if len(mask) == n_missing or len(mask) == 1:
                print('Missing sufficient bias ratios for variable: {} {}'.\
                        format(grid_var, layer),
                        '\nNeed at least two stations with data, skipping.')
                return
            print('Warning:\n',
                    'Data missing for {} of {} stations for variable: {} {}'.\
                    format(n_missing, len(values), grid_var, layer),
                    '\nproceeding with interpolation.')
            # get locations where ratio is not nan
            values = values[mask]
            lon_pts = lon_pts[mask]
            lat_pts = lat_pts[mask]

        nx_cells = int(np.round(np.abs((lon_min - lon_max) / CELL_SIZE)))
        ny_cells = int(np.round(np.abs((lat_min - lat_max) / CELL_SIZE)))
        # rbf requires uniform grid (n X n) so 
        # extend short dimension and clip later 
        nx_cells_out = copy.copy(nx_cells)
        ny_cells_out = copy.copy(ny_cells)
        # gdal requires "upper left" corner coordinates
        lat_max_out = copy.copy(lat_max)
        lon_max_out = copy.copy(lon_max)
        # extend short dimension to make square grid
        if not nx_cells == ny_cells:
            diff = np.abs(nx_cells - ny_cells)
            if nx_cells > ny_cells:
                lat_max += diff * CELL_SIZE
                ny_cells += diff
            else:
                lon_max += diff * CELL_SIZE
                nx_cells += diff

        # make finer/coarse grid by scale factor
        lons = np.linspace(lon_min, lon_max, 
                int(np.round(nx_cells/scale_factor))+1)
        lats = np.linspace(lat_min, lat_max, 
                int(np.round(ny_cells/scale_factor))+1)
        # extent for original created by spatial.build_subgrid
        # add one to make sure raster covers full extent
        lons_out = np.linspace(lon_min, lon_max_out, 
                int(np.round(nx_cells_out/scale_factor))+1)
        lats_out = np.linspace(lat_min, lat_max_out, 
                int(np.round(ny_cells_out/scale_factor))+1)

        # if function was 'linear_rbf' 
        function = function.replace('_rbf', '')
        # make sampling square grid
        XI, YI = np.meshgrid(lons, lats)
        # apply rbf interpolation
        rbf = Rbf(lon_pts, lat_pts, values, function=function, smooth=smooth)
        ZI = rbf(XI, YI)
        # clip to original extent, rbf array flips axes, and row order... 
        ZI_out = ZI[0:len(lats_out),0:len(lons_out)]
        ZI_out = np.flip(ZI_out,axis=0)

        #### save scipy interpolated data as raster 
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

        # calc residuals add to shapefile and in_path CSV, move shape to out_dir
        calc_pt_error(in_path, out_dir, layer, grid_var)
        # calculate zonal statistics save means for each gridMET cell
        if zonal_stats:
            gridmet_zonal_stats(in_path, out_file)

        
    # run gdal_grid interpolation 
    if function in InterpGdal.interp_methods:
        if not bounds:
            bounds = get_subgrid_bounds(in_path, buffer=buffer) 
        lon_min, lon_max, lat_min, lat_max = bounds
        gg = InterpGdal(in_path)
        gg.gdal_grid(layer=layer, out_dir=out_dir, interp_meth=function,
                    params=params, bounds=bounds, scale_factor=scale_factor,
                    zonal_stats=zonal_stats, options=options)
        
    # scipy radial basis function interpolation for now
    # run interpolation and zonal statistics depending on layer kwarg
    else: 
        if layer == 'all': # potential for multiprocessing
            for l in InterpGdal.default_layers:
                _run_rbf_interpolation(l, bounds, function, smooth)
        # single layer option
        elif isinstance(layer, str):
            _run_rbf_interpolation(layer, bounds, function, smooth)
        # run select list or tuple of layers
        elif isinstance(layer, (list, tuple)):
            for l in layer:
                _run_rbf_interpolation(l, bounds, function, smooth)

def calc_pt_error(in_path, out_dir, layer, grid_var):
    """
    Calculate point ratio estimates from interpolated raster, residuals,
    and add to output summary CSV and point shapefile. Make copies of
    updated files and saves to directory with interpolated rasters.
    
    Arguments:
        in_path (str): path to comprehensive summary CSV created by 
            :mod:`gridwxcomp.calc_bias_ratios`
        out_dir (str): path to dir that contains interpolated raster
        layer (str): layer to calculate error e.g. "annual_mean"
        grid_var (str): name of gridMET variable e.g. "etr_mm"

    Returns:
        None

    Note:
        This function should be run **after** :func:`make_points_file`
        because it copies data from the shapefile it created.
    """
    raster = str(Path(out_dir)/'{}.tiff'.format(layer))
    pt_shp = '{}_summary_pts.shp'.format(grid_var)
    pt_shp = str(Path(in_path).parent/'spatial'/pt_shp)

    if not Path(pt_shp).is_file():
        make_points_file(in_path)

    pt_shp_out = str(Path(out_dir)/'{}_summary_pts.shp'.format(grid_var))
    # mean fields in point shapefile does not include '_mean'
    pt_layer = layer.replace('_mean', '')
    if pt_layer == 'growseason':
        pt_layer = 'grow'
    # names of new fields for estimated and residual e.g. Jan_est, Jan_res
    pt_est = '{}_est'.format(pt_layer)
    pt_res = '{}_res'.format(pt_layer)
    
    print('\nExtracting interpolated data at station locations and \n',
        'calculating residuals for layer:', layer)
    pt_err = pd.DataFrame(columns=[pt_est, pt_res])
    # read raster for layer and get interpolated data for each point
    with fiona.open(pt_shp) as shp:
        for feature in shp:
            STATION_ID = feature['properties']['STATION_ID']
            coords = feature['geometry']['coordinates']
            # Read pixel value at the given coordinates using Rasterio
            # sample() returns an iterable of ndarrays.
            with rasterio.open(raster) as src:
                value = [v for v in src.sample([coords])][0][0]
            # store interpolated point estimates of ratios 
            pt_err.loc[STATION_ID, pt_est] = value

    # merge estimated point data with observed to calc residual
    pt_err['STATION_ID'] = pt_err.index
    # read summary CSV with observed ratios
    in_df = pd.read_csv(in_path, index_col='STATION_ID', na_values=[-999])
    in_df.loc[pt_err.index, pt_est] = pt_err.loc[:, pt_est]
    # calculate residual estimated minus observed
    in_df.loc[:,pt_res] = in_df.loc[:,pt_est] - in_df.loc[:,layer]
    # save/overwrite error to input CSV for future interpolation 
    in_df.to_csv(in_path, index=True, na_rep=-999)

    # save copy of CSV with updated error info to out_dir with rasters
    out_summary_csv = Path(out_dir)/Path(in_path).name
    if not out_summary_csv.is_file():
        in_df.to_csv(str(out_summary_csv), index=True, na_rep=-999)
    else:
        out_df = pd.read_csv(str(out_summary_csv), index_col='STATION_ID')
        out_df.loc[pt_err.index, pt_est] = pt_err.loc[:, pt_est]
        out_df.loc[pt_err.index, pt_res] = in_df.loc[pt_err.index, pt_res]
        out_df.to_csv(out_summary_csv, index=True, na_rep=-999)
    
    # error info to new point shapefile in raster directory
    if not Path(pt_shp_out).is_file():
        with fiona.open(pt_shp, 'r') as inf:
            schema = inf.schema.copy()
            input_crs = inf.crs
            # add attributes for point estimate and residual to output points
            schema['properties'][pt_est] = 'float'
            schema['properties'][pt_res] = 'float'
            with fiona.open(pt_shp_out, 'w', 'ESRI Shapefile', 
                    schema, input_crs) as outf:
                for feat in inf:
                    STATION_ID = feat['properties']['STATION_ID']
                    feat['properties'][pt_est] =\
                            in_df.loc[STATION_ID, pt_est].astype(float)
                    feat['properties'][pt_res] =\
                            in_df.loc[STATION_ID, pt_res].astype(float)
                    outf.write(feat)
    # if already exists update point shapefile
    else:
        tmp_out = pt_shp_out.replace('.shp', '_tmp.shp')
        with fiona.open(pt_shp_out, 'r') as inf:
            schema = inf.schema.copy()
            input_crs = inf.crs
            # add attributes for point estimate and residual to output points
            schema['properties'][pt_est] = 'float'
            schema['properties'][pt_res] = 'float'
            with fiona.open(tmp_out, 'w', 
                    'ESRI Shapefile', schema, input_crs) as outf:
                for feat in inf:
                    STATION_ID = feat['properties']['STATION_ID']
                    feat['properties'][pt_est] =\
                            in_df.loc[STATION_ID, pt_est].astype(float)
                    feat['properties'][pt_res] =\
                            in_df.loc[STATION_ID, pt_res].astype(float)
                    outf.write(feat)
        # keep tmp point file with new data and remove old version
        for f in os.listdir(out_dir):
            if '_tmp.' in f:
                move(OPJ(out_dir, f), OPJ(out_dir, f.replace('_tmp', '')))

    # remove point shapefile from "spatial" directory and tmp files
    spatial_dir = Path(in_path).parent.joinpath('spatial')
    for f in os.listdir(spatial_dir):
        if Path(f).stem == '{}_summary_pts'.format(grid_var):
            # delete temp point shapefile
            (Path(in_path).parent/'spatial'/f).resolve().unlink()

    # delete tmp summary csv used in interpgdal _make_pt_vrt method 
    tmp_csv = str(in_path).replace('.csv','_tmp.csv')
    if Path(tmp_csv).resolve().is_file():
        Path(tmp_csv).resolve().unlink()


def gridmet_zonal_stats(in_path, raster):
    """
    Calculate zonal means from interpolated surface of etr bias ratios
    created by :func:`interpolate` using the fishnet grid created by 
    :func:`make_grid`. Save mean values for each gridMET cell to
    a CSV file joined to gridMET IDs. 
    
    Arguments:
        in_path (str): path to [var]_summary_comp.csv file containing 
            monthly bias ratios, lat, long, and other data. Created by 
            :mod:`gridwxcomp.calc_bias_ratios`. 
        raster (str): path to interpolated raster of bias ratios to
            be used for zonal stats. First created by :func:`interpolate`.
        
    Example:
        Although it is prefered to use this function as part of 
        :func:`interpolate` or indirectly using the :mod:`gridwxcomp.spatial`
        command line usage. However if the grid shapefile and spatial
        interpolated raster(s) have already been created without zonal
        stats then,
        
        >>> from gridwxcomp import spatial
        >>> # assign input paths
        >>> summary_file = 'monthly_ratios/etr_mm_summary_comp.csv'  
        >>> raster_file = 'monthly_ratios/spatial/etr_mm_invdist_400m/Jan_mean.tiff'
        >>> spatial.gridmet_zonal_stats(summary_file, raster_file)
        
        The final result will be the creation of::
            
            'monthly_ratios/spatial/etr_mm_invdist_400m/gridMET_stats.csv'
            
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
        If zonal statistics are estimated for the same variable on the
        same raster more than once, the contents of that column in the 
        gridMET_stats.csv file will be overwritten. 
        
    """
    if not os.path.isfile(in_path):
        raise FileNotFoundError('Input summary CSV file given'+\
                                ' was invalid or not found')
    # look for fishnet created in 'in_path/spatial'
    path_root = os.path.split(in_path)[0]
    file_name = os.path.split(in_path)[1]
    # get variable names from input file prefix
    grid_var = file_name.split('_summ')[0]
    var_name = Path(raster).name.split('.')[0]
    # grid is always in the "spatial" subdir of in_path
    grid_file = OPJ(path_root, 'spatial', 'grid.shp')
    # save zonal stats to summary CSV in same dir as raster as of version 0.3
    raster_root = os.path.split(raster)[0]
    out_file = OPJ(raster_root, 'gridMET_stats.csv')

    # this error would only occur when using within Python 
    if not os.path.isfile(grid_file):
        raise FileNotFoundError(
            os.path.abspath(grid_file),
            '\ndoes not exist, create it using spatial.make_grid first'
        )
    print(
        'Calculating', grid_var, 'zonal means for', var_name
    )

    # calc zonal stats and open for grid IDs
    with fiona.open(grid_file, 'r') as source:
        zs = zonal_stats(source, raster, all_touched=True)
        gridmet_ids = [f['properties'].get('GRIDMET_ID') for f in source]

    # get just mean values, zonal_stats can do other stats...
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
        # overwrite column values if exists, else append
        existing_df = pd.read_csv(out_file)
        existing_df.GRIDMET_ID = existing_df.GRIDMET_ID.astype(int)
        if var_name in existing_df.columns:
            # may throw error if not same size as original grid
            try:
                existing_df.update(out_df)
                existing_df.to_csv(out_file, index=False)   
            except:
                print('Zonal stats for this variable already exist but they',
                      'appear to have been calculated with a different grid',
                      'overwriting existing file at:\n',
                      os.path.abspath(out_file)
                )
                out_df.to_csv(out_file, index=False)
        else:
            existing_df = existing_df.merge(out_df, on='GRIDMET_ID')
            #existing_df = pd.concat([existing_df, out_df], axis=1).drop_duplicates()
            existing_df.to_csv(out_file, index=False)   
    
    
def arg_parse():
    """
    Command line usage of grdwxcomp spatial.py for creating shapefiles of 
    climate station point data of bias ratios, build fishnet around stations,
    perform spatial interpolation of ratios, and extract zonal means for gridMET    
    cells.
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
        '-l', '--layer', metavar='', required=False, default='all', nargs='*',
        help='Layer(s) to perform spatial interpolation in input summary '+\
             'CSV, if multiple use comma separation e.g. Jan_mean,Aug_mean')
    optional.add_argument(
        '-o', '--out-dir', metavar='PATH', required=False,
        help='Subdirectory to save interpolated rasters')
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
        '-f', '--function', required=False, default='invdist', type=str,
        metavar='', help='Function to use for spatial interpolation, gdal'+\
                ' methods include invdist, invdistnn, average, nearest,'+\
                ' and linear. Radial basis functions include: multiquadric,'+\
                ' inverse, gaussian, linear_rbf, cubic, quintic, thin_plate.') 
    optional.add_argument(
        '--smooth', required=False, default=0, type=float, metavar='',
        help='Smooth parameter for radial basis func interpolation methods')
    optional.add_argument(
        '-p', '--params', required=False, default=None, type=str, metavar='',
        help='Parameter string e.g. ":power=2:smoothing=.1:nodata=-999" for'+\
            ' gdal_grid spatial interpolation methods')
    optional.add_argument(
        'z', '--zonal-stats', required=False, default=True, 
        action='store_false', help='Flag to NOT extract zonal means of'+\
            ' interpolated rasters to gridMET cells.')
    optional.add_argument(
        '--overwrite-grid', required=False, default=False, 
        action='store_true', help='Flag to overwrite existing fishnet grid')
    optional.add_argument(
        '--options', required=False, default=None, type=str, metavar='',
        help='Additional command line options for gdal_grid interpolation')
    optional.add_argument(
        '-g', '--gridmet-meta', metavar='', required=False, default=None,
        help='GridMET metadata CSV file with cell data, packaged with '+\
             'gridwxcomp and automatically found if pip was used to install '+\
             'if not given it needs to be located in the currect directory')
#    optional.add_argument(
#        '--debug', default=logging.INFO, const=logging.DEBUG,
#        help='Debug level logging', action="store_const", dest="loglevel")
    parser._action_groups.append(optional)# to avoid optionals listed first
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = arg_parse()
        
    # parse layers if multiple
    layer = args.layer
    layer = layer.split(',')
    if len(layer) == 1:
        layer = layer[0] # get single layer as string

    main(
        input_file_path=args.input,
        layer=args.layer,
        out=args.out_dir,
        buffer=args.buffer,
        scale_factor=args.scale,
        function=args.function,
        smooth=args.smooth,
        params=args.params,
        zonal_stats=args.zonal_stats,
        overwrite=args.overwrite_grid,
        options=args.options,
        gridmet_meta_path=args.gridmet_meta
    )
