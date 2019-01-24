# -*- coding: utf-8 -*-
"""
Tools to estimate the spatial surface of monthly ratios between climate station
and gridMET estimated reference evapotranspiration. 

Attributes:
    CELL_SIZE (float): constant gridMET cell size in decimal degrees.

Note:
    All spatial files, i.e. vector and raster files, utilize the
    *ESRI Shapefile* file format. 

Todo:
    * add spatial interpolation, logging 

"""

import os
import fiona
import ogr
import argparse
import pandas as pd
import numpy as np
from math import ceil
from shapely.geometry import Point, Polygon, mapping
from fiona import collection
from fiona.crs import from_epsg

# decimal degrees gridMET cell size
CELL_SIZE = 0.041666666666666664

OPJ = os.path.join

def main(input_file_path, gridmet_meta_path=None, buffer=25):
    """
    Create point shapefile of monthly mean bias ratios from comprehensive
    CSV file created by :mod:`calc_bias_ratios`. Build fishnet grid around
    the climate stations that matches the gridMET dataset. Perform spatial
    interpolation to estimate 2-dimensional surface of bias ratios and extract
    interpolated values to gridMET cells. 
    
    Arguments:
        input_file_path (str): path to summary_comp.CSV file containing monthly
            bias ratios, lat, long, and other data. Shapefile "summary_pts.shp"
            is saved to parent directory of ``in_path`` under a subdirectory
            named "spatial".

    Keyword Arguments:
        gridmet_meta_path (str): path to metadata CSV file that contains
            all gridMET cells for the contiguous United States. If None
            it is looked for at ``etr-biascorrect/gridmet_cell_data.csv``.
        buffer (int): number of gridMET cells to expand the rectangular extent
            of the subgrid fishnet (default=25). 

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
        gridMET cells added around the encompassed climate stations. 

        Alternatively, to make a larger buffer zone on the fishnet grid 
        used for interpolation, use the -b or --buffer option

        .. code::
            $ python spatial.py -i monthly_ratios/summary_comp.csv -g gridmet_cell_data.csv -b 30

        For use within Python see examples of :func:`make_points_file` 
        and :func:`build_subgrid`.

    Note:
        The input file is first created by :mod:`calc_bias_ratios`.

    """
    # build point shapefile 
    make_points_file(input_file_path)
    # build fishnet for interpolation
    build_subgrid(
        input_file_path, 
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
            of the subgrid fishnet (default=25). 
        
    Returns:
        bounds (tuple): tuple with coordinates in decimal degrees that
            define the outer bounds of the subgrid fishnet in the format
            (min long, max long, min lat, max lat)

    Raises:
        FileNotFoundError: if input summary CSV file is not found.

    Note:
        By expanding the grid to a larger area encompassing the climate 
        stations of interest spatial interpolation is less skewed by the 
        boundary values.
        
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
        gridmet_meta_path (str): path to metadata CSV file that contains
            all gridMET cells for the contiguous United States. If None
            it is looked for at ``etr-biascorrect/gridmet_cell_data.csv``.
            
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
                grid_path, 
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
            bias ratios, lat, long, and other data.
            
    Keyword Arguments:
        gridmet_meta_path (str): path to metadata CSV file that contains
            all gridMET cells for the contiguous United States. If None
            it is looked for at ``etr-biascorrect/gridmet_cell_data.csv``.
        buffer (int): number of gridMET cells to expand the rectangular extent
            of the subgrid fishnet (default=25). 
            
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

def arg_parse():
    """
    Parse command line arguments for creating shapefile of all stations
    found in comprehensive summary CSV file, build fishnet around stations
    and perform spatial interpolation for gridMET cells.
    """
    parser = argparse.ArgumentParser(
        description='Perform spatial analysis to estimate spatial'+\
                ' distribution of etr bias ratios.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '-i', '--input', metavar='PATH', required=True,
        help='Input CSV file of merged climate/gridMET data that '+\
             'was created by running prep_input.py, download_gridmet_ee.py'+\
             ', and calc_bias_ratios.py')
    parser.add_argument(
        '-g', '--gridmet', metavar='PATH', required=False,
        help='GridMET master CSV file with cell data, packaged with '+\
             'etr-biascorrect at etr-biascorrect/gridmet_cell_data.csv '+\
             'if not given it needs to be located in the currect directory')
    parser.add_argument(
        '-b', '--buffer', required=False, default=25, type=int,
        help='Number of gridMET cells to expand outer bounds of fishnet '+\
             'in order to avoid skewing of interpolation values by boundary')
#    parser.add_argument(
#        '--debug', default=logging.INFO, const=logging.DEBUG,
#        help='Debug level logging', action="store_const", dest="loglevel")
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = arg_parse()

    main(
        input_file_path=args.input, 
        gridmet_meta_path=args.gridmet, 
        buffer=args.buffer
    )
