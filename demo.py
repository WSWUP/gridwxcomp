import pyproj
from gridwxcomp.prep_metadata import prep_metadata, get_subgrid_bounds
from gridwxcomp.ee_download import download_grid_data
from gridwxcomp.plot import daily_comparison, monthly_comparison, station_bar_plot
from gridwxcomp.calc_bias_ratios import calc_bias_ratios
from gridwxcomp.spatial import make_grid, interpolate


gridded_dataset_name = 'conus404'  # name of the dataset comparison is being made with
station_meta_path = './gridwxcomp/example_data/Station_Data.txt'  # Path to station metadata file with lat/long coords
config_path = f'./gridwxcomp/example_data/gridwxcomp_config_{gridded_dataset_name}.ini'  # local path for config file

export_bucket = 'openet'  # name of bucket data will be exported to once the earth engine calls complete
export_path = f'bias_correction_gridwxcomp_testing/gridwxcomp_{gridded_dataset_name}/'  # path in bucket to export to

gridwxcomp_input = f'{gridded_dataset_name}_gridwxcomp_metadata.csv'  # local path for prep_metadata output
output_dir = f'bias_outputs_{gridded_dataset_name}'  # Directory that bias ratio/interpolation outputs will be saved to

# Defining the parameters for the interpolation
search_radius = 15*100*1000
params = {'power': 4, 'smoothing': 0.0, 'max_points': 25, 'min_points': 0, 'radius': search_radius, 'nodata': -999}
lcc_interpolation_resolution = 1000  # in meters, since interpolation will happen in Lambert Conformal Conic

'''
prep_metadata
   The purpose of this module is to open the file provided at 'station_meta_path' and verify
        it has everything needed in order to proceed with acquiring gridded data and comparing it against
        the observed station data
'''
prep_metadata(station_meta_path, gridded_dataset_name, out_path=gridwxcomp_input)

'''
define_projection_parameters
    Functions in gridwxcomp.spatial and gridwxcomp.interpgdal are expecting projection information in a dictionary
    
    This might seem cumbersome but if this is fleshed out a little more we could have gridwxcomp function with any two
    EPSG projections, although the first is required to be the same projection as the input data coordinates
    
    This should eventually be made its own function under prep_metadata.py, where it asks for:
        - two projection names (arbitrary, user will refer to them later in other functions)
        - two projection resolutions
        - optionally, bounds for the first projection, if none provided uses spatial.get_subgrid_bounds(),
            which should also be reclassified under prep_metadata.py
'''

# TODO: this should be made a function
projection_dict = {
    'wgs84': {  # WGS84 is EPSG:4326 and is in decimal degrees
        'bounds': {},
        'resolution': 0.1,
        'crs_id': 'EPSG:4326'},
    'lcc': {  # Lambert Conformal Conic is ESRI:102004 and is in meters
        'bounds': {},
        'resolution': 1000,
        'crs_id': 'ESRI:102004'}
}

# User has option to specify bounds
# projection_dict['wgs84']['bounds'] = {'xmin': -115.0, 'xmax': -101.0, 'ymin': 35.5, 'ymax': 42.5}

if len(projection_dict['wgs84']['bounds']) == 0:
    # Calculate interpolation grid bounds if not provided
    # TODO improve alignment of get_subgrid_bounds, currently wont snap on NLDAS resolution (0.125) but will for (0.1)
    projection_dict['wgs84']['bounds'] = get_subgrid_bounds(gridwxcomp_input, projection_dict['wgs84']['resolution'],
                                                            1, buffer=25)

# Convert WGS84 coordinates into LCC coordinates and add them to projection_dict
lcc_transformer = pyproj.Transformer.from_crs(projection_dict['wgs84']['crs_id'],
                                              projection_dict['lcc']['crs_id'], always_xy=True)
# Calculate xmin at the SW corner
projection_dict['lcc']['bounds']['xmin'], projection_dict['lcc']['bounds']['_ignore'] = (
    lcc_transformer.transform(projection_dict['wgs84']['bounds']['xmin'], projection_dict['wgs84']['bounds']['ymin']))
# Calculate xmax at the SE corner
projection_dict['lcc']['bounds']['xmax'], projection_dict['lcc']['bounds']['_ignore'] = (
    lcc_transformer.transform(projection_dict['wgs84']['bounds']['xmax'], projection_dict['wgs84']['bounds']['ymin']))
# Calculate ymax at the NE corner, could've also been NW corner
projection_dict['lcc']['bounds']['_ignore'], projection_dict['lcc']['bounds']['ymax'] = (
    lcc_transformer.transform(projection_dict['wgs84']['bounds']['xmax'], projection_dict['wgs84']['bounds']['ymax']))
# Calculate lcc's ymin as at the average between the east and west extent
projection_dict['lcc']['bounds']['_ignore'], projection_dict['lcc']['bounds']['ymin'] = (
    lcc_transformer.transform((projection_dict['wgs84']['bounds']['xmax']+projection_dict['wgs84']['bounds']['xmin'])/2, projection_dict['wgs84']['bounds']['ymin']))

# Round the entries in the LCC dict to the nearest meter
for key in projection_dict['lcc']['bounds'].keys():
    projection_dict['lcc']['bounds'][key] = round(projection_dict['lcc']['bounds'][key], 0)

'''
download_grid_data
   The purpose of this module is to make calls to the earth engine API to export gridded timeseries at the station
        points contained within the station metadata file generated by prep_metadata. Exported data will be saved to
        a cloud storage bucket that the user specifies.

   The user will need to manually download the data via gsutil before proceeding further. The module amends the
        station metadata file to have paths for the gridded dataset as they would exist locally once the user has
        downloaded them.

   Example: the example run is with the dataset conus404, the gridded path will be filled in such that if you
        download the 'conus404' folder from the bucket to the local directory of trial_runs.py it will work without
        further modification.
'''
download_grid_data(gridwxcomp_input, dataset=gridded_dataset_name, export_bucket=export_bucket,
                   export_path=export_path, force_download=False, authorize=False)

'''
plotting
    The purpose of this module is to generate bokeh plots at both the daily and monthly timesteps to visualize
        differences/similarities between the station and gridded datasets. Both methods will generate a subplot
        for every variable which is shared between the two datasets. Each timeseries plot is accompanied by
        a scatterplot featuring a LSR line that is forced through zero, giving a sense of the overall bias and
        correlation between the two datasets

    The 'monthly_comparison()' method will generate monthly averages to get a better sense of the overall bias
        by reducing the noise of daily observations (Ex. observations from Jan 1st - 31st 2020 become one average).
        This generates one plot per station included in the input file.

    The 'daily_comparison()' will generate 12 plots per station, one for each month of the year. These months are
        selected and plotted together, which allows the user to visualize if one particular year has had a substantial
        deviation for other years in the record

    These plots are mainly diagnostic and are not used as part of the later steps.

    The plotting module also contains a function to visualize the results from 'calc_bias_ratios' or 'interpolate'
        in the form of bar plots.

'''
monthly_comparison(gridwxcomp_input, config_path, dataset_name=gridded_dataset_name)
daily_comparison(gridwxcomp_input, config_path, dataset_name=gridded_dataset_name)

for var in ['tmax', 'tmin', 'eto']:  # Iterate over vars in list. Valid entries found in calc_bias_ratios.py VAR_LIST
    ratio_filepath = f'{output_dir}/gridded_{var}_summary_comp_1980_2020.csv'  # path to bias ratios output file
    interpolation_out_path = (f'{var}_invdistnn_p{params["power"]}_'  # directory for interpolation outputs
                              f's{params["smoothing"]}_maxpoints{params["max_points"]}_radius{params["radius"]}')

    '''
    calc_bias_ratios
       The purpose of this module is to calculate the monthly bias
            (either the long-term-mean [default] or mean-of-annual)
            between the gridded dataset and the station observations.
    
       IMPORTANT: The bias factors are either
            (station - gridded) [Temperature] or (station / gridded) [Everything else].
            This is the reverse relationship from how we normally think of modeled / observed, so a temperature bias
            of -2 means that the gridded dataset is hotter than observed, not colder. This reversal was done to make
            applying bias a matter of addition/multiplication.
    
       This module requires the configuration .ini file (stored under the example_data folder)
            to be set up so that it can correct interpret the columns within the station and gridded data files.
            The acceptable parameters for 'comparison_var' are coded in ACCEPTABLE_VAR_LIST within calc_bias_ratios.py
    '''
    calc_bias_ratios(gridwxcomp_input, config_path, output_dir, method='long_term_mean', grid_id_name='GRID_ID',
                     comparison_var=var, grid_id=None, day_limit=10, years='1980-2020', comp=True)

    '''
    station_bar_plot
        Generates boxplots of the bias ratios to visualize overall performance.
        Requires the outputs of calc_bias_ratios as an input.
    '''
    station_bar_plot(ratio_filepath, bar_plot_layer='growseason_mean')

    '''
    make grid
        This will generate a fishnet grid to Make fishnet grid (vector file of polygon geometry)
            based on bounding coordinates. Each cell in the grid will be assigned
            a unique numerical identifier (property specified by grid_id_name) and the grid will be in
            the WGS84 coordinate system.

        This step is not required if the desired output is just the interpolation surfaces. It becomes necessary
            zonal statistics are going to be generated by spatial.interpolate.
    '''
    make_grid(ratio_filepath, projection_dict['wgs84']['resolution'],
              bounds=projection_dict['wgs84']['bounds'], overwrite=False, grid_id_name='GRID_ID')

    '''
    interpolate
        Contains several methods to create a 2-D interpolated surface based on the lat/long points and values
            generated by calc_bias_ratios.py. Provides options allow for up/down-scaling the resolution of the
            resampling grid and to select from multiple interpolation methods.

            Interpolation will occur in lambert conformal conic, however the end product
            Interploated surfaces are saved as GeoTIFF rasters and will be in the WGS84 CRS.

        The default behavior will also calculate zonal statistics (which requires having run spatial.make_grid).
            To disable this you must specifically pass z_stats=False
    '''
    interpolate(ratio_filepath, projection_dict, proj_name='lcc', layer='all', out=interpolation_out_path, scale_factor=1,
                function='invdistnn', params=params, buffer=5, z_stats=False, res_plot=False,
                grid_id_name='GRID_ID', options=None)