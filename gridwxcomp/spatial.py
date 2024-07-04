# -*- coding: utf-8 -*-
""" 
Perform multiple workflows needed to estimate the spatial surface (and
other related outputs) of monthly and annual station-to-grid bias ratios for
meteorological variables. The input file is first created by
:mod:`gridwxcomp.calc_bias_ratios`.

Attributes:
    PT_ATTRS (tuple): all attributes expected to be in point shapefile
        created for stations except station and Grid IDs.

Note:
    All spatial files, i.e. vector and raster files, utilize the
    *ESRI Shapefile* or GeoTiff format. 
    
Todo:
    * logging 

"""

import os
import copy
from math import ceil
from pathlib import Path
from shutil import move

import fiona
import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
from shapely.geometry import Point, mapping
from fiona import collection
from fiona.crs import from_epsg
from osgeo import gdal, osr, ogr
from rasterstats import zonal_stats as zstats


# point shapefile attributes except for station and grid IDs
PT_ATTRS = (
   'Jan_mean', 'Feb_mean', 'Mar_mean', 'Apr_mean', 'May_mean', 
   'Jun_mean', 'Jul_mean', 'Aug_mean', 'Sep_mean', 'Oct_mean', 
   'Nov_mean', 'Dec_mean', 'Jan_count', 'Feb_count', 'Mar_count', 
   'Apr_count', 'May_count', 'Jun_count', 'Jul_count', 'Aug_count', 
   'Sep_count', 'Oct_count', 'Nov_count', 'Dec_count', 'Jan_stdev', 
   'Feb_stdev', 'Mar_stdev', 'Apr_stdev', 'May_stdev', 'Jun_stdev', 
   'Jul_stdev', 'Aug_stdev', 'Sep_stdev', 'Oct_stdev', 'Nov_stdev', 
   'Dec_stdev', 'Jan_cv', 'Feb_cv', 'Mar_cv', 'Apr_cv', 'May_cv', 
   'Jun_cv', 'Jul_cv', 'Aug_cv', 'Sep_cv', 'Oct_cv', 'Nov_cv', 'Dec_cv',
   'growseason_mean', 'summer_mean', 'annual_mean', 'growseason_count', 
   'summer_count', 'annual_count', 'growseason_stdev', 'summer_stdev', 
   'annual_stdev', 'growseason_cv', 'summer_cv', 'annual_cv'
)

OPJ = os.path.join


def make_points_file(in_path, grid_id_name='GRID_ID'):
    """
    Create vector shapefile of points with monthly mean bias ratios 
    for climate stations using all stations found in a comprehensive
    CSV file created by :mod:`gridwxcomp.calc_bias_ratios`.
    
    Arguments:
        in_path (str): path to [var]_summary_comp.CSV file containing 
            monthly bias ratios, lat, long, and other data. Shapefile 
            "[var]_summary_pts.shp" is saved to parent directory of 
            ``in_path`` under "spatial" subdirectory.
        grid_id_name (str): name of the column containing grid ID's
            
    Returns:
        None
    
    Example:
        Create shapefile containing point data for all climate stations
        in input summary file created by :mod:`gridwxcomp.calc_bias_ratios`
        
        >>> from gridwxcomp import spatial
        >>> # path to comprehensive summary CSV
        >>> summary_file = 'monthly_ratios/etr_mm_summary_comp_all_yrs.csv'
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
        raise FileNotFoundError('Input summary CSV file: {}\nwas not found.'.format(in_path))
    print('\nMapping point data for climate stations in: \n', in_path, '\n')

    in_df = pd.read_csv(in_path, index_col='STATION_ID', na_values=[-999])
    # add in potentially missing columns to avoid errors when no ratios exist
    # in input that are expected by schema/attribute table
    missing_vars = list(set(PT_ATTRS).difference(in_df.columns))
    in_df = in_df.reindex(columns=list(in_df.columns) + missing_vars)
    # save shapefile to "spatial" subdirectory of in_path
    path_root = os.path.split(in_path)[0]
    file_name = os.path.split(in_path)[1]
    # get variable name from input file prefix
    var_name = file_name.split('_summ')[0]
    out_dir = OPJ(path_root, 'spatial')
    out_file = OPJ(out_dir, '{v}_summary_pts_wgs84.shp'.format(v=var_name))
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

    crs = from_epsg(4326)  # WGS 84 projection
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
            grid_id_name: 'int'
        }}

    # remove nans- gdal will not recognize  
    in_df = in_df.where(pd.notnull(in_df), None)

    # create shapefile from points, overwrite if exists
    with collection(out_file, 'w', driver='ESRI Shapefile', crs=crs, schema=schema) as output:
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
                    grid_id_name: row[grid_id_name]
                },
                'geometry': mapping(point)
            })

    # Now that the shapefile has been created, open it and
    # reproject it from wsg84 to LCC before saving
    # ESRI:102004 is what I have inferred the correct projection is from John Volk's workflow
    lcc_out_file = OPJ(out_dir, '{v}_summary_pts_lcc.shp'.format(v=var_name))
    wgs84_shapefile = gpd.read_file(out_file)
    lcc_shapefile = wgs84_shapefile.to_crs(crs='ESRI:102004')
    lcc_shapefile.to_file(lcc_out_file)


def make_grid(in_path, grid_res, bounds, overwrite=False, grid_id_name='GRID_ID'):
    """
    Make fishnet grid (vector file of polygon geometry) for 
    select gridcells based on bounding coordinates. Each cell in the grid will be assigned
    a unique numerical identifier (property specified by grid_id_name) and the grid will be in
    the WGS84 coordinate system.
    
    Modified from the
    `Python GDAL/OGR Cookbook <https://pcjericks.github.io/py-gdalogr-cookbook/vector_layers.html#create-fishnet-grid>`.
    
    Arguments:
        in_path (str): path to [var]_summary_comp_[years].csv file containing 
            monthly bias ratios, lat, long, and other data. Created by 
            :func:`gridwxcomp.calc_bias_ratios`.
        grid_res (float):  Cell size of grid in decimal degrees
        bounds (dict): dict of bounding coordinates
            the rectangular extent of the subgrid fishnet and interpolated
            output raster.
        overwrite (bool): default False. If True, overwrite the grid
            shapefile at ``out_path`` if it already exists.
        grid_id_name (str): default "GRID_ID". Name of gridcell identifier column


    Returns:
        None
        
    Examples:
        Build a fishnet uniform grid that matches gridded data cells around
        climate stations found in ``in_path`` summary CSV file with a 
        25 cell buffer. The gridded dataset has a resolution of 0.1 decimal degrees.
        
        >>> from gridwxcomp import spatial
        >>> # assign input paths
        >>> summary_file = 'monthly_ratios/etr_mm_summary_comp_all_yrs.csv'
        >>> # make fishnet of grid cells for interpolation
        >>> spatial.make_grid(summary_file, 0.1)
            
        The file will be saved as "grid.shp" to a newly created subdirectory
        "spatial" in the same directory as the input summary CSV file. i.e.::
        
            monthly_ratios/
            ├── etr_mm_summary_all_yrs.csv
            ├── etr_mm_summary_comp_all_yrs.csv
            └── spatial/
                ├── grid.cpg
                ├── grid.dbf
                ├── grid.prj
                ├── grid.shp
                └── grid.shx
        
        To use a smaller buffer to the extent of the grid assign the 
        ``buffer`` keyword argument
        
        >>> spatial.make_grid(summary_file, 0.1, buffer=5)
            
        If the grid file already exists the default action is to not 
        overwrite. To overwrite an existing grid if, for example, you 
        are working with a new set of climate stations as input, then 
        set the ``overwrite`` keyword argument to True. 
        
        >>> spatial.make_grid(summary_file, 0.1, overwrite=True, buffer=5)
            
    Raises:
        FileNotFoundError: if input summary CSV file is not found. 

        
    """
        
    if not os.path.isfile(in_path):
        raise FileNotFoundError('Input summary CSV file given was invalid or not found')
    # save grid to "spatial" subdirectory of in_path
    path_root = os.path.split(in_path)[0]
    out_dir = OPJ(path_root, 'spatial')
    out_path = OPJ(out_dir, 'grid.shp')
    
    # skip building grid if already exists
    if os.path.isfile(out_path) and not overwrite:
        print('\nFishnet grid already exists at: \n', out_path, '\nnot overwriting.\n')
        return
    # print message if overwriting existing grid
    elif os.path.isfile(out_path) and overwrite:
        print('\nOverwriting fishnet grid at: \n', out_path, '\n')

    # create output directory if it does not exist
    if not os.path.isdir(out_dir):
        print(os.path.abspath(out_dir), ' does not exist, creating directory.\n')
        os.mkdir(out_dir)

    # extract values from dict for readability
    xmin = bounds['xmin']
    xmax = bounds['xmax']
    ymin = bounds['ymin']
    ymax = bounds['ymax']

    # read path and make parent directories if they don't exist    
    if not os.path.isdir(path_root):
        print(
            os.path.abspath(path_root), 
            ' does not exist, creating directory.\n'
        )
        os.makedirs(path_root)
        
    print(
        '\nCreating fishnet grid for subset of gridcells, \n',
        '\nSouthwest corner (lat, long): {:9.4f}, {:9.4f}'.format(ymin, xmin),
        '\nNortheast corner (lat, long): {:9.4f}, {:9.4f}'.format(ymax, xmax),
    )
    
    # get n rows
    rows = ceil((ymax-ymin) / grid_res)
    # get n columns
    cols = ceil((xmax-xmin) / grid_res)

    # start grid cell envelope
    ringXleftOrigin = xmin
    ringXrightOrigin = xmin + grid_res
    ringYtopOrigin = ymax
    ringYbottomOrigin = ymax - grid_res

    # create output file
    outDriver = ogr.GetDriverByName('ESRI Shapefile')
    if os.path.exists(out_path):
        os.remove(out_path)
    outDataSource = outDriver.CreateDataSource(out_path)
    outLayer = outDataSource.CreateLayer(out_path, geom_type=ogr.wkbPolygon)
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
            ringYtop = ringYtop - grid_res
            ringYbottom = ringYbottom - grid_res

        # new envelope for next poly
        ringXleftOrigin = ringXleftOrigin + grid_res
        ringXrightOrigin = ringXrightOrigin + grid_res

    # Save and close DataSources
    outDataSource = None
    
    print('\nFishnet shapefile successfully saved to: \n', os.path.abspath(out_path), '\n')
    # reopen grid and assign grid attributes and coord. ref.
    _update_subgrid(out_path, grid_id_name=grid_id_name)


def _update_subgrid(grid_path, grid_id_name='GRID_ID'):
    """
    Helper function to assign grid ID values and EPSG WGS 84 
    projection to fishnet grid.
    
    Arguments:
        grid_path (str): path to fishnet grid shapefile for subset of
            gridcells which contain climate stations being analyzed.
        grid_id_name (str): the name of the property that will be added to the grid shapefile

    Returns:
        None
        
    Raises:
        FileNotFoundError: if ``grid_path``  are not found.
        
    """

    if not os.path.isfile(grid_path):
        raise FileNotFoundError('The file path for the grid fishnet was invalid or does not exist. ')

    tmp_out = grid_path.replace('.shp', '_tmp.shp')

    # WGS 84 projection
    crs = from_epsg(4326) 

    # overwrite fishnet grid with updated grid_id field
    with fiona.open(grid_path, 'r') as source:
        print('Adding grid IDs ({}) to fishnet grid, saving to: \n'.format(grid_id_name),
              os.path.abspath(grid_path), '\n')
        
        n_cells = len([f for f in source])
        print('Assigning values for ', n_cells, ' gridcells.\n')
        
        # Copy the source schema and add grid name property
        sink_schema = source.schema
        sink_schema['properties'][grid_id_name] = 'int'
        # overwrite file add spatial reference
        with fiona.open(tmp_out, 'w', crs=crs, driver=source.driver, schema=sink_schema) as sink:
            # add grid_id feature to outfile
            grid_id = 1000000
            for feature in source:
                feature['properties'][grid_id_name] = grid_id
                sink.write(feature)
                grid_id += 1
    # cannot open same file and write to it on Windows, overwrite temp
    root_dir = os.path.split(grid_path)[0]
    for f in os.listdir(root_dir):
        if '_tmp' in f:
            move(OPJ(root_dir, f), OPJ(root_dir, f.replace('_tmp', '')))
    print(
        'Completed assigning grid IDs to fishnet. \n'
    )

    
def interpolate(in_path, proj_dict, proj_name, layer='all', out=None, scale_factor=1,
                function='invdist', params=None, buffer=25, z_stats=False, res_plot=True,
                grid_id_name='GRID_ID', options=None):
    """
    Use various methods to interpolate a 2-dimensional surface of
    calculated bias ratios or other statistics for station/gridded data
    pairs found in input comprehensive summary CSV. 
    
    Options allow for modifying down- or up-scaling the resolution of the 
    resampling grid and to select from multiple interpolation methods.
    Interploated surfaces are saved as GeoTIFF rasters. Zonal statistics 
    using :func:`zonal_stats` are also extracted to grid cells in
    the fishnet grid built first by :func:`make_grid`. 
    
    Arguments:
        in_path (str): path to [var]_summary_comp_[years].csv file containing 
            monthly bias ratios, lat, long, and other data. Created by 
            :func:`gridwxcomp.calc_bias_ratios`.
        proj_dict (dict): Dictionary containing CRS information, grid resolution, and bounds.
        proj_name (str): entry in proj_dict to pull projection info from
        layer (str or list): default 'all'. Name of variable(s) in ``in_path``
            to conduct 2-D interpolation. e.g. 'Annual_mean'.
        out (str): default None. Subdirectory to save GeoTIFF raster(s).
        scale_factor (float, int): default 0.1. Scaling factor to apply to 
            original grid resolution to create resampling resolution. If 
            scale_factor = 0.1, the interpolation resolution will be
            one tenth of the grid resolution
        function (str): default 'invdist'. Interpolation method, gdal 
            methods include: 'invdist', 'indistnn', 'linear', 'average',
            and 'nearest' see `gdal_grid <https://www.gdal.org/gdal_grid.html>`.
        params (dict, str, or None): default None. Parameters for interpolation
            using gdal, see defaults in :class:`gridwxcomp.InterpGdal`.
        buffer (int): default 25. Number of grid cells to expand 
            the rectangular extent of the subgrid fishnet and interpolated
            output raster.
        z_stats (bool): default True. Calculate zonal means of interpolated
            surface to grid cells in fishnet and save to a CSV file.
            The CSV file will be saved to the same directory as the interpolated
            raster file(s).            
        res_plot (bool): default True. Make bar plot for residual (error)
            between interpolated and station value for ``layer``.
        grid_id_name (str): default "GRID_ID". Name of gridcell identifier
        options (str or None): default None. Extra command line arguments
            for gdal interpolation.

    Returns:
        None
        
    Examples:
        Let's say we wanted to interpolate the "Annual_mean" bias 
        ratio in an input CSV first created by :func:`gridwxcomp.calc_bias_ratios` and a fishnet
        grid was first created by :func:`make_grid`. This example uses the 
        "invdist" method (default) to interpolate to a 0.1 decimal degree grid scaled down
        to a 0.01 decimal degree surface. The result is a GeoTIFF raster that has an extent that
        encompasses station locations in the input file plus an additional 
        optional buffer of outer grid cells. Additionally, point residuals
        of bias ratios are added to CSV and newly created point shapefiles, 
        zonal (grid cell) means are also extracted and stored in a CSV.
        
        >>> from gridwxcomp import spatial
        >>> summary_file = 'monthly_ratios/etr_mm_summary_comp_all_yrs.csv'
        >>> buffer = 10
        >>> layer = 'annual_mean'
        >>> params = {'power':1, 'smooth':20}
        >>> proj_dict = {
        >>>     'wgs84': {'bounds': {'xmin': -115.0, 'xmax': -101.0, 'ymin': 35.5, 'ymax': 42.5},
        >>>               'resolution': 0.1, 'crs_id': 'EPSG:4326'}}
        >>> out_dir = 's20_p1' # optional subdir name for saving rasters
        >>> interpolate(summary_file, proj_dict, 'wgs84', layer=layer, out=out_dir,
        >>>     scale_factor=0.1, params=params, buffer=buffer)

        The resulting file structure that is created by the above command is::

            monthly_ratios/
            ├── etr_mm_summary_all_yrs.csv
            ├── etr_mm_summary_comp_all_yrs.csv
            └── spatial/
                ├── etr_mm_invdist_400m/
                │ └── s20_p1/
                │     ├── annual_mean.tiff
                │     ├── annual_mean.vrt
                │     ├── etr_mm_summary_comp_all_yrs.csv
                │     ├── etr_mm_summary_pts_wgs84.cpg
                │     ├── etr_mm_summary_pts_wgs84.dbf
                │     ├── etr_mm_summary_pts_wgs84.prj
                │     ├── etr_mm_summary_pts_wgs84.shp
                │     ├── etr_mm_summary_pts_wgs84.shx
                │     ├── etr_mm_summary_pts_lcc.cpg
                │     ├── etr_mm_summary_pts_lcc.dbf
                │     ├── etr_mm_summary_pts_lcc.prj
                │     ├── etr_mm_summary_pts_lcc.shp
                │     ├── etr_mm_summary_pts_lcc.shx
                │     ├── zonal_stats.csv
                │     └── residual_plots
                │         └── annual_res.html
                ├── grid.cpg
                ├── grid.dbf
                ├── grid.prj
                ├── grid.shp
                └── grid.shx
            
        Specifically, the interpolated raster is saved to::
        
            'monthly_ratios/spatial/etr_mm_invdist_400m/s20_p1/Annual_mean.tiff'
            
        where the file name and directory is based on the variable being 
        interpolated, methods, and the raster resolution. The ``out`` 
        keyword argument lets us add any number of subdirectories to the final 
        output directory, in this case the 's20_p1' dir contains info on params.
        
        The final result will be the creation of the CSV::
            
            'monthly_ratios/spatial/etr_mm_invdist_400m/s20_p1/zonal_stats.csv'
            
        In "zonal_stats.csv" the zonal mean for ratios of annual station to 
        grid ETr will be stored along with grid IDs, e.g.
        
            ==========  =================
            GRID_ID  Annual_mean
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
        >>> interpolate(summary_file, dataset_resolution, layer=layer, scale_factor=2,
        >>>     function=func, buffer=buffer)
                    
        This will create the GeoTIFF raster::
        
            'monthly_ratios/spatial/etr_mm_invdistnn_400m/April_to_oct_cv.tiff'
               
        And the zonal means will be added as a column named "April_to_oct_cv"
        to:: 
        
            'monthly_ratios/spatial/etr_mm_invdistnn_400m/zonal_stats.csv'

        As with other components of ``gridwxcomp``, any other climatic
        variables that exist in the grid dataset can be used along
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
        however, if the buffer and input [var]_summary_comp_[years].csv files 
        arguments differ from those used for :func:`interpolate` the 
        raster may not fully cover the fishnet which may result in 
        gaps in the zonal statistics.
        
    """
    # avoid circular import for InterpGdal for gdal interpolation methods
    try:
        from gridwxcomp.interpgdal import InterpGdal
    except: # for use as script, i.e. $ python spatial ...
        from interpgdal import InterpGdal
    if not os.path.isfile(in_path):
        raise FileNotFoundError('Input summary CSV file given was invalid or not found')

    # calc raster resolution
    interpolation_res = scale_factor * proj_dict[proj_name]['resolution']
    # path to save raster of interpolated grid scaled by scale_factor
    path_root = os.path.split(in_path)[0]
    file_name = os.path.split(in_path)[1]
    # get variable name from input file prefix
    grid_var = file_name.split('_summary')[0]
    
    if not out: 
        out_dir = OPJ('spatial', '{}_{}_{:.1f}_meters'.format(grid_var, function, interpolation_res))
    elif out == str(Path(in_path).parent):
        out_dir = OPJ('spatial', '{}_{}_{:.1f}_meters'.format(grid_var, function, interpolation_res))
        print(
            'WARNING: output subdirectory for rasters cannot be named '
            'the same as the parent directory holding the input '
            'summary CSV file. Output will be saved to:\n{}'.format(out_dir)
        )
    else:
        out_dir = OPJ('spatial', '{}_{}_{:.1f}_meters'.format(grid_var, function, interpolation_res), out)

    # run gdal_grid interpolation
    if function in InterpGdal.interp_methods:
        gg = InterpGdal(in_path)
        gg.gdal_grid(proj_dict, proj_name, layer=layer, out_dir=out_dir, interp_meth=function,
                     params=params, scale_factor=scale_factor, z_stats=z_stats, res_plot=res_plot,
                     grid_id_name=grid_id_name, options=options)
    else:
        # Invalid method, raise error
        raise ValueError(f'Interpolation method "{function}" not included in list of GDAL methods of interpolation\n'
                         f'must be one of: {InterpGdal.interp_methods}')


def calc_pt_error(in_path, out_dir, layer, grid_var, grid_id_name='GRID_ID'):
    """
    Calculate point ratio estimates from interpolated raster, residuals,
    and add to output summary CSV and point shapefile. Make copies of
    updated files and saves to directory with interpolated rasters.

    This function will search for the LCC projection of the points shapefile
    
    Arguments:
        in_path (str): path to comprehensive summary CSV created by 
            :mod:`gridwxcomp.calc_bias_ratios`
        out_dir (str): path to dir that contains interpolated raster
        layer (str): layer to calculate error e.g. "annual_mean"
        grid_var (str): name of grid variable e.g. "etr_mm"

    Returns:
        None

    Note:
        This function should be run **after** :func:`make_points_file`
        because it copies data from the shapefile it created.
    """
    raster = str(Path(out_dir)/'{}_lcc.tiff'.format(layer))
    pt_shp = '{}_summary_pts_lcc.shp'.format(grid_var)
    pt_shp = str(Path(in_path).parent/'spatial'/pt_shp)

    if not Path(pt_shp).is_file():
        make_points_file(in_path, grid_id_name=grid_id_name)

    pt_shp_out = str(Path(out_dir)/'{}_summary_pts_lcc.shp'.format(grid_var))
    # mean fields in point shapefile does not include '_mean'
    pt_layer = layer.replace('_mean', '')
    if pt_layer == 'growseason':
        pt_layer = 'grow'
    # names of new fields for estimated and residual e.g. Jan_est, Jan_res
    pt_est = '{}_est'.format(pt_layer)
    pt_res = '{}_res'.format(pt_layer)
    
    print('\nExtracting interpolated data at station locations and calculating residuals for layer: ', layer)
    pt_err = pd.DataFrame(columns=[pt_est, pt_res])
    # read raster for layer and get interpolated data for each point
    with fiona.open(pt_shp) as shp:
        for feature in shp:
            STATION_ID = feature['properties']['STATION_ID']
            coords = feature['geometry']['coordinates']
            # Read pixel value at the given coordinates using Rasterio
            # sample() returns an iterable of ndarrays.
            try:
                with rasterio.open(raster) as src:
                    value = [v for v in src.sample([coords])][0][0]
            except:
                raise Exception('ERROR: at least one station location does not overlap with the interpolated raster.')
            # store interpolated point estimates of ratios 
            pt_err.loc[STATION_ID, pt_est] = value

    # merge estimated point data with observed to calc residual
    pt_err['STATION_ID'] = pt_err.index
    # read summary CSV with observed ratios
    in_df = pd.read_csv(in_path, index_col='STATION_ID', na_values=[-999])
    in_df.loc[pt_err.index, pt_est] = pt_err.loc[:, pt_est]
    # calculate residual estimated minus observed
    in_df.loc[:, pt_res] = in_df.loc[:, pt_est] - in_df.loc[:, f'{layer}_mean']
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
                            float(in_df.loc[STATION_ID, pt_est])
                    feat['properties'][pt_res] =\
                            float(in_df.loc[STATION_ID, pt_res])
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
                            float(in_df.loc[STATION_ID, pt_est])
                    feat['properties'][pt_res] =\
                            float(in_df.loc[STATION_ID, pt_res])
                    outf.write(feat)
        # keep tmp point file with new data and remove old version
        for f in os.listdir(out_dir):
            if '_tmp.' in f:
                move(OPJ(out_dir, f), OPJ(out_dir, f.replace('_tmp', '')))

    # # remove point shapefile from "spatial" directory and tmp files
    # spatial_dir = Path(in_path).parent.joinpath('spatial')
    # for f in os.listdir(spatial_dir):
    #     if Path(f).stem == '{}_summary_pts_lcc'.format(grid_var):
    #         # delete temp point shapefile
    #         (Path(in_path).parent/'spatial'/f).resolve().unlink()
    #
    # # delete tmp summary csv used in interpgdal _make_pt_vrt method
    # tmp_csv = str(in_path).replace('.csv', '_tmp.csv')
    # if Path(tmp_csv).resolve().is_file():
    #     Path(tmp_csv).resolve().unlink()


def zonal_stats(in_path, raster, grid_id_name='GRID_ID'):
    """
    Calculate zonal means from interpolated surface of etr bias ratios
    created by :func:`interpolate` using the fishnet grid created by 
    :func:`make_grid`. Save mean values for each gridcell to
    a CSV file joined to grid IDs. 
    
    Arguments:
        in_path (str): path to [var]_summary_comp_[years].csv file containing 
            monthly bias ratios, lat, long, and other data. Created by 
            :mod:`gridwxcomp.calc_bias_ratios`.
        raster (str): path to interpolated raster of bias ratios to
            be used for zonal stats. First created by :func:`interpolate`.
        
    Example:
        Although it is prefered to use this function as part of 
        :func:`interpolate` or indirectly using the :mod:`gridwxcomp.spatial`
        command line usage. However, if the grid shapefile and spatial
        interpolated raster(s) have already been created without zonal
        stats then,
        
        >>> from gridwxcomp import spatial
        >>> # assign input paths
        >>> summary_file = 'monthly_ratios/etr_mm_summary_comp_[years].csv'  
        >>> raster_file = 'monthly_ratios/spatial/etr_mm_invdist_400m/Jan_mean.tiff'
        >>> spatial.zonal_stats(summary_file, raster_file)
        
        The final result will be the creation of::
            
            'monthly_ratios/spatial/etr_mm_invdist_400m/grid_stats.csv'
            
        The resulting CSV contains the grid IDS and zonal means
        for each grid cell in the fishnet which must exist at::
        
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
        zonal_stats.csv file will be overwritten. 
        
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
    # grid is in the "spatial" subdir of in_path
    grid_file = OPJ(path_root, 'spatial', 'grid.shp')
    # save zonal stats to summary CSV in same dir as raster as of version 0.3
    raster_root = os.path.split(raster)[0]
    out_file = OPJ(raster_root, 'zonal_stats.csv')

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
        zs = zstats(source, raster, all_touched=True)
        grid_ids = [f['properties'].get(grid_id_name) for f in source]

    # get just mean values, zonal_stats can do other stats...
    means = [z['mean'] for z in zs]
    out_df = pd.DataFrame(
        data={
            grid_id_name: grid_ids, 
            var_name: means
        }
    )
    out_df[grid_id_name] = out_df[grid_id_name].astype(int)
    # drop rows for cells outside of grid master grid
    out_df = out_df.drop(out_df[out_df[grid_id_name] == -999].index)

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
        existing_df[grid_id_name] = existing_df[grid_id_name].astype(int)
        if var_name in existing_df.columns:
            # may throw error if not same size as original grid
            try:
                existing_df.update(out_df)
                existing_df.to_csv(out_file, index=False)   
            except:
                print('Zonal stats for this variable already exist but they',
                      'appear to have been calculated with a different grid',
                      'overwriting existing file at:\n',os.path.abspath(out_file))
                out_df.to_csv(out_file, index=False)
        else:
            existing_df = existing_df.merge(out_df, on=grid_id_name)
            existing_df.to_csv(out_file, index=False)   


def get_subgrid_bounds(in_path, config_path, buffer=0):
    """
    Calculate bounding box for spatial interpolation grid from
    output of prep_metadata.prep_metadata(), will attempt to make it snap
    to grid by querying catalog on earth engine

    Arguments:
        in_path (str): path to metadata file created by
            :func:`gridwxcomp.prep_metadata.prep_metadata()`.
        config_path (str): path to config file containing projection/catalog
        info buffer (int): number of grid cells to expand the rectangular
            extent of the subgrid fishnet.

    Returns:
        bounds (dict): dictionary with coordinates that
            define the outer bounds of the subgrid fishnet in the format
            (min long, max long, min lat, max lat)

    Raises:
        FileNotFoundError: if input summary CSV file is not found.

    Notes:
        By expanding the grid to a larger area encompassing the climate
        stations of interest :func:`interpolate` can be used to extrapolate
        passed the bounds of the outer station locations.

        You must authenticate with Google Earth Engine before using
        this function.

    """

    # TODO - a potential oversight of this function is that it assumes
    #   the CRS of the metadata file and the CRS of the EE catalog are the same
    #   right now we're testing with the CONUS404 collection and its in meters
    #   Idea:
    #       As a stopgap raise an error if the resolution of the ee dataset is
    #       Value that wouldn't make sense for WGS84

    def _find_nearest(val, array):
        """Element in numpy array `array` closest to the scalar value `val`"""
        idx = (np.abs(array - val)).argmin()
        return array[idx]

    if not os.path.isfile(in_path):
        raise FileNotFoundError('Input metadata file given was invalid or not found')

    # read metadata from prep_metadata
    in_df = pd.read_csv(in_path)

    # read in config information from file

if __name__ == '__main__':
    print('\n--------------------------------------------------------'
          ' Functionality for running this library from the terminal'
          ' was removed. Please refer to the documentation on how to'
          ' make calls to these functions. \n\n')
