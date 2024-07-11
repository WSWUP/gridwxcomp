Tutorial
========

This tutorial describes in detail how to use the ``gridwxcomp`` Python
package including preparation of input data, downloading gridded data
and pairing it with station weather data, calculating monthly bias
ratios between station and gridded data, spatial interpolation of point
results, and generating interactive graphics files.

.. code:: python
    
    # module and function imports
    import pandas as pd
    import numpy as np
        
    from gridwxcomp import prep_metadata
    from gridwxcomp import calc_bias_ratios
    from gridwxcomp import download_grid_data
    from gridwxcomp import prep_metadata
    from gridwxcomp import spatial
    from gridwxcomp import plot

Input data files and formatting requirements
--------------------------------------------

There are three input files used by ``gridwxcomp``, they include:

1. configuration text file [.INI]
2. text or Microsoft Excel files containing time series of station
   weather data [.CSV, .XLS, XLSX]
3. text file containing metadata regarding the weather stations [.CSV]

An example of each of these files is included in with this package, they
can be found in the “example_data” folder in the package after
installation, if you installed using PiP you can find the path to the
example data station metadata file where the other example data also
exists (in the same directory) by running:

.. code:: python

   import pkg_resources; print(pkg_resources.resource_filename('gridwxcomp', 'example_data/Station_Data.txt'))

Or you can directly download the example files from GitHub
`here <https://github.com/WSWUP/gridwxcomp/tree/master/gridwxcomp/example_data>`__.

**Note:** To follow this tutorial it is recommended to start a Python 
script or Jupyter Notebook from within the provided "example_data" folder.

The configuration file
~~~~~~~~~~~~~~~~~~~~~~

The configuration file is a text file with the “.INI” extension and is
read using the Python :obj:`configparser.ConfigParser` class to read
in metadata about the control parameters that the user specifies for
``gridwxcomp``.

The configuration file has three sections, **METADATA**, **DATA**, and
**UNITS**, and these sections are defined by there names inside of
brackets, e.g., the units section comes after the line that contains the
text “[UNITS]”.

The **METADATA** section includes information on the spatial domain,
resolution, and projection for spatial interpolation of results, as well
as the name of the gridded dataset to download from Google Earth
including the dataset’s collection path and start and end date of data
to download. Most of these options are self explanatory, here is an
example **METADATA** snippet from the provided example configuration
file that was built to download gridded data from the CONUS404 dataset:

::

   [METADATA]
   # Projection information
   #   gridwxcomp will reproject data for interpolation and output if requested. 
   #   Fill out the lines below with EPSG codes and desired resolutions
   #   If no reprojection is desired then set input/interpolation/output to the same values.
   input_data_projection = EPSG:4326
   input_data_resolution = 0.1
   interpolation_projection = ESRI:102004
   interpolation_resolution = 1000
   output_data_projection = EPSG:4326
   output_data_resolution = 0.1


   # Bounding information
   #   You can manually specify the bounds/extents for the interpolation area below.
   #   If any of xmin, xmax, ymin, ymax are left blank then gridwxcomp will 
   #   generate these values automatically
   xmin = -111.8708332996666428
   xmax = -108.6208332996662733
   ymin = 38.0874999999668162
   ymax = 40.5874999999666741


   # Gridded dataset information
   #   Specify the Earth Engine image collection you'd like to use for comparison.
   #   collection_name will be used in the generation of filenames
   #   You may also specify the start and end dates (Format: YYYY-MM-DD) of the data to download.
   #   If the dates are left blank then gridwxcomp will generate these values automatically.
   collection_name = conus404
   collection_path = projects/openet/assets/meteorology/conus404/daily
   start_date = 
   end_date = 


   # File structure information
   #   These values are necessary for gridwxcomp to parse the data files.
   #   Use 'station' for observed and 'gridded' for any model data
   #   Anemometer height required in meters
   station_anemometer_height = 2
   station_lines_of_header = 1
   station_missing_data_value = nan

   gridded_anemometer_height = 10
   gridded_lines_of_header = 1
   gridded_missing_data_value = nan

**Note:** The station and gridded data wind speed height are needed so
that the wind speed variables can both be scaled to 2 m using the
logarithmic vertical velocity profile, see equation 33 in [Allen2005]_.

The second section of the configuration file is called **DATA**; this
section is exclusivly for the user to specify the names of the station
and gridded weather data as they are found in the station weather data
CSV files (in the headers) and as they are named for the specified
Google Earth Engine data collection. Here is an example for the CONUS404
dataset and the provided weather data:

::

   [DATA]
   # For the below parameters, enter the name of the column containing the following values
   #   If a column is not provided, leave the parameter blank.

   station_date_col = date
   station_tmax_col = TMax (C)
   station_tmin_col = TMin (C)
   station_rs_col = Rs (w/m2)
   station_wind_col = ws_2m (m/s)
   station_ea_col =
   station_tdew_col = TDew (C)
   station_rhmax_col = RHMax (%)
   station_rhmin_col = RHMin (%)
   station_rhavg_col = RHAvg (%)
   station_eto_col = ETo (mm)
   station_etr_col = ETr (mm)

   gridded_date_col = date
   gridded_tmax_col = T2_MAX
   gridded_tmin_col = T2_MIN
   gridded_rs_col = ACSWDNB
   gridded_wind_col = WIND10
   gridded_ea_col = 
   gridded_tdew_col = TD2
   gridded_rhmax_col =
   gridded_rhmin_col =
   gridded_rhavg_col =
   gridded_eto_col = ETO_ASCE
   gridded_etr_col = ETR_ASCE

The final and third section of the ``gridwxcomp`` configuration input
file is the **UNITS** section, which as the name implies, allows the
user to specify the units of the station and gridded weather data that
the software will parse. This is critical so that the software can
convert units is necessary so that they match before computing
station:gridded monthy bias ratios. The unit conversion is done by the
:func:`gridwxcomp.calc_bias_ratios` function. Here is an example of
this section from the provided example data:

::

   [UNITS]
   # For the parameters in this section, enter the corresponding units from the options commented above.

   # K, F, C
   station_temp_units = C
   gridded_temp_units = K

   # kw-hr/m2, j/m2, mj/m2, langleys, w/m2
   station_solar_units = w/m2
   gridded_solar_units = j/m2

   # m/s, mph, kmph
   station_wind_units = m/s
   gridded_wind_units = m/s

   # kPa, torr, mbar
   station_ea_units = kpa
   gridded_ea_units =

   # percent, fraction
   station_rh_units = percent
   gridded_rh_units =

   # inches, mm
   station_et_units = mm
   gridded_et_units = mm

.. _variable_list:
Weather variables processed by ``gridwxcomp``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The available input options for weather variables and their units
currently allowed by ``gridwxcomp`` are as follows:

   =================== ================================================================== ======================================
   Variable             Description                                                        Allowable Unit(s)         
   =================== ================================================================== ======================================
   tmax, tmin, tdew     maximum, minimum and dew point air temperature                     c, f, k
   rs                   solar radiation                                                    kw-hr/m2, j/m2, mj/m2, langleys, w/m2
   wind                 wind speed                                                         m/s, mph, kmph                                
   ea                   vapor pressure                                                     kPa, torr, mbar
   rhmax, rhmin, rhavg  maximum, minimum and average relative humidity                     percent, fraction                              
   eto, etr             short (grass) and tall (alfalfa) ASCE standardized reference ET    inches, mm
   =================== ================================================================== ======================================

The converted weather variables will not be written to files, they are
converted so that the pairing of station:gridded data can be done before
computing and saving average bias ratios or temperature differences.

The weather station’s data files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Files containing daily time series of weather station data are the key
input to ``gridwxcomp``. These files should be formatted as comma
separated variable [.CSV] text files or Microsoft Excel files [.XLS or .
XLSX]. The names of variables that can be used by ``gridwxcomp`` should
be listed in the configuration file and they should match the data as
they are found in the weather station and gridded data file headers.
Here is an example of the first three rows and first seven columns of an
example weather station data:


+---------------------+----------+----------+----------+----------+------------------+-----------+
| date                | TAvg (C) | TMax (C) | TMin (C) | TDew (C) | Vapor Pres (kPa) | RHAvg (%) |
+=====================+==========+==========+==========+==========+==================+===========+
| 2013-11-07 00:00:00 | 4.382    | 15.83    | -4.331   | -4.7     | 0.431            | 55.25     |
+---------------------+----------+----------+----------+----------+------------------+-----------+
| 2013-11-08 00:00:00 | 4.005    | 19.3     | -7.252   | -5.65    | 0.401            | 55.65     |
+---------------------+----------+----------+----------+----------+------------------+-----------+
| 2013-11-09 00:00:00 | 3.019    | 19.1     | -6.842   | -4.98    | 0.422            | 54.95     |
+---------------------+----------+----------+----------+----------+------------------+-----------+


**Note:** The “date” column in the provided weather data will be parsed
by :mod:`Pandas` and should be in a format that is able to
automatically converted to a :obj:`Pandas.datetime` object. For
example, “YYYY/MM/DD” or “YYYY-MM-DD HH:MM:SS”

The weather station’s metadata file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Within the same folder of the station weather data files the user must
provide a text file [.CSV] that lists all the weather stations that are
to be included in the ``gridwxcomp`` routines and for each station, this
file lists some key metadata. There are four columns that are required
by ``gridwxcomp`` to be provided in this file: ‘Latitude’, ‘Longitude’,
‘Filename’, and ‘Station’. Filename refers to the name of the weather
station data file, e.g., “BedrockCO_Daily_output.xlsx”. The “Station”
column should contain the ID that the user wants to use for that station
and this will be used for output file names that apply to that station
and in different outputs, e.g., the ID given to to the stations in the
bias ratio files and point shapefiles. Here is an example of a station
metadata file with the four required columns:

+-------------------------+------------------+-------------------+----------------------------------+---------+
| Station                 | Latitude         | Longitude         | Filename                         | Elev_FT |
+=========================+==================+===================+==================================+=========+
| Bluebell (Neola Area)   | 40.3723213601075 | -110.209184085302 | BluebellUT_Daily_output.xlsx     | 6186    |
+-------------------------+------------------+-------------------+----------------------------------+---------+
| Loa                     | 38.3834675639262 | -111.635832870077 | LoaUT_Daily_output.xlsx          | 7116    |
+-------------------------+------------------+-------------------+----------------------------------+---------+
| Bedrock                 | 38.328297440752  | -108.855494308994 | BedrockCO_Daily_output.xlsx      | 4973    |
+-------------------------+------------------+-------------------+----------------------------------+---------+
| Castle Valley near Moab | 38.6429447999517 | -109.398808843297 | CastleValleyUT_Daily_output.xlsx | 4687    |
+-------------------------+------------------+-------------------+----------------------------------+---------+



**Note:** Any additional columns that exist in the weather station
metadata file will be retained and added to the formatted output CSV
file that is produced by the :func:`gridwxcomp.prep_metadata`
function. However they will not be used by any of the following
procedures, only the four required columns’ values are used (‘Latitude’,
‘Longitude’, ‘Filename’, and ‘Station’). In the exampe above, the extra
columns that were provided are “Elev_FT” and “Location”.


Step 1: Parse input data
------------------------

The first step to running ``gridwxcomp`` after preparing the required
input data as specified in
:ref:`Input data files and formatting requirements` is to run the
:func:`gridwxcomp.prep_metadata` function which reads the station
metadata file and prepares for downloading gridded data. This step is 
straightforward with minimal options involved:

.. code:: python

    # specify the paths to input data files, in this case using the provided example data:
    station_meta_path = '/home/john/gridwxcomp/gridwxcomp/example_data/Station_Data.txt'
    conus404_config = '/home/john/gridwxcomp/gridwxcomp/example_data/gridwxcomp_config_conus404.ini'
    gridded_dataset_name = 'conus404'
    
    # run the function 
    prep_metadata(station_meta_path, conus404_config, gridded_dataset_name)


The file that was produced from running
:func:`gridwxcomp.prep_metadata` is named “formatted_input.csv” by
default and it will be saved to the workspace where the function is
called from unless otherwise stated in the . It has updated the paths to the station weather data and
reformatted the station metadata file. This will be the input file used
for the next two steps in the ``gridwxcomp`` workflow which are
:func:`gridwxcomp.ee_download` and
:func:`gridwxcomp.calc_bias_ratios`.

Step 2: Download gridded timeseries data from Google Earth Engine
-----------------------------------------------------------------

After running :func:`gridwxcomp.prep_metadata` the next step is to use
the formatted CSV file that was created alongwith the configuration
input file as input to download the specified gridded data that
corresponds with the locations and variables of the weather stations.
Some of that required information is in the configuration file, such as
the dataset collection path on Google Earth Engine and its name. Some
data required to download Earth Engine gridded climate data needs to be
specified as arguments to the :func:`gridwxcomp.ee_download.download_grid_data`
function, such as the bucket to export the extracted point time series
data to and the local folder to download the same data to.

**Important:** Before downloading data using the Earth Engine Python
API, the use must initialize Earth Engine locally and have permissions
to access the requested data as well as to export data on the Google
Cloud. After setting up Google Earth Engine locally following the
`online
instructions <https://developers.google.com/earth-engine/guides/python_install>`__,
one can initialize Earth Engine in Python using the following line:

.. code:: python

   import ee
   ee.Authenticate()
   ee.Initialize(project='my-project')

Now we can download gridded data:

.. code:: ipython

    # Specify the path to the file created by running prep_metadata
    formatted_input_file = '/home/john/gridwxcomp/gridwxcomp/example_data/formatted_input.csv'
    
    import ee
    ee.Initialize()
    # download the gridded data
    download_grid_data(
        formatted_input_file, 
        conus404_config, 
        export_bucket='openet', # bucket root to export to
        export_path=f'bias_correction_gridwxcomp_testing/gridwxcomp_conus404/', # path to export data to
        local_folder=None, # If not specified then the gridded data will be downloaded to a new folder
        force_download=False, # if False check if data already exists locally, if True overwrite
    )


**Note:** If the start and end dates for downloading gridded weather
data are not specified in the configuration file, the entire period of
record of gridded data will be downloaded for each station (at the
overlapping grid cell).

After running :func:`gridwxcomp.ee_download.download_grid_data` time series of the
weather data will be saved to a folder that is named using the gridded
data collection name as specified in the configuration file. This folder
will be created where the download function is called, in this case in
the “example_data” folder. The individual files containing the gridded
time series at the station locations will be named using the gridded
dataset name, the station name, and the start and end dates that were
used for downloading, for example:
``"[collection_name]_[station]_[start_date]_[end_date]_all_vars.csv"``

Here is the file structure that should have been produced after up to
this stage assuming that the “example_data” folder was used as the
working space for running this tutorial:

::

   example_data/
   ├── BedrockCO_Daily_output.xlsx
   ├── BluebellUT_Daily_output.xlsx
   ├── CastleValleyUT_Daily_output.xlsx
   ├── conus404
   │   ├── conus404_bedrock_19791001_20220928_all_vars.csv
   │   ├── conus404_bluebell_neola_area_19791001_20220928_all_vars.csv
   │   ├── conus404_castle_valley_near_moab_19791001_20220928_all_vars.csv
   │   └── conus404_loa_19791001_20220928_all_vars.csv
   ├── formatted_input.csv
   ├── gridwxcomp_config_conus404.ini
   ├── LoaUT_Daily_output.xlsx
   └── Station_Data.txt

At this step in the normal workflow of ``gridwxcomp`` the output file created by :func:`gridwxcomp.ee_download.download_grid_data` can be used for making interactive daily and monthly time series and scatter plots of paired station and gridded weather data using the :mod:`gridwxcomp.plot` module.

Step 3: Calculate monthly, seasonal, and annual station:gridded biases and statistics
-------------------------------------------------------------------------------------

After parsing the input station weather data and configuration options,
and downloading the corresponding gridded weather data of choice, the
next step in the ``gridwxcomp`` workflow is computing station:gridded
biases. This process involved pairing the station and gridded time
series together for overlapping time periods, making necessary unit
conversions, and computing monthly, seasonal, and annual average bias
ratios (or differences for air temperature) between the station and
gridded data for each variable that is available or specified.
Additional metrics are calculated that are helpful to evaluate the
variability in the station:gridded ratios such as the annual standard
deviation and coefficients of variation for the bias ratios or
differences, as well as the number of paired data points used to compute
the bias ratios or differences. In addition to calculating long-term
average monthly bias ratios or differences between station:gridded data,
summer periods (JJA), growing season (AMJJASO), and annual periods are
also used for computing the metrics.

To run the bias corrections, the :func:`gridwxcomp.calc_bias_ratios`
reads the formatted metadata file created by
:func:`gridwxcomp.prep_metadata` and the configuration file. The user
should also specify the folder to save the output file, which variable
to use for the calculations from the list of available variables: see
:ref:`variable_list`, the maximum number of gaps days per month
allowed for computations (``day_limit`` kwarg to
:func:`gridwxcomp.calc_bias_ratios`, default is ten days maximum of
gap days), and the year range to use for the calculations in case one is
not interested in using the full data record.

There are two methods for calculating the bias ratios or differences,
the “long_term_mean” and the “mean_of_annual”. The default method
(``method='long_term_mean'``) first groups the paired station and
gridded data for each time period (monthly, etc.) and then takes the
average of station and gridded data respectively before taking the ratio
or difference, for example,

.. math::  \frac{ \frac{\sum_{i=1}^{n} station_i}{n}} {\frac{\sum_{i=1}^{n} grid_i}{n}} 

where :math:`station_i` and :math:`grid_i` are the :math:`i^{th}` paired
daily weather data in the full record for a given temporal period, such
as all the summer days or all the days that fall within the month of
May. For air temperature variables, as opposed to taking the ratio the
calculation is

.. math::   \frac{\sum_{i=1}^{n} station_i}{n} - \frac{\sum_{i=1}^{n} grid_i}{n}. 

The other option for calculating the bias ratios or temperature
differences between station and gridded data
(``method='mean_of_annual'``) is similar except it makes the calculation
as shown above for each year in the paired data record separately, and
then it takes the average of those annual ratios or differences. This
approach is always used for calculating the statndard deviation and
coefficient of variation variables that are also computed by the
:func:`gridwxcomp.calc_bias_ratios` function.

Here is an example code using the default method, next we will examine
the output:

.. code:: python3

    # directory to save results of point calculations
    output_dir = 'test_data_bias_results'
    
    calc_bias_ratios(
        input_path=formatted_input_file,
        config_path=conus404_config,
        out_dir=output_dir,
        method='long_term_mean',
        comparison_var='wind'
    )


Here is a selection of the results for the month of January from the
output CSV file that was created which was named
“wind_summary_comp_all_yrs.csv”:

+-----------------------+-------------------+-----+-----------+-----+-----------+-----+--------+
| STATION_ID            | Jan_mean          | ... | Jan_count | ... | Jan_stdev | ... | Jan_cv |
+=======================+===================+=====+===========+=====+===========+=====+========+
| Bluebell (Neola Area) | 0.648012105097692 | ... | 62        | ... | 0.108     | ... | 0.163  |
+-----------------------+-------------------+-----+-----------+-----+-----------+-----+--------+
| Loa                   | 1.15758848442987  | ... | 62        | ... | 0.001     | ... | 0      |
+-----------------------+-------------------+-----+-----------+-----+-----------+-----+--------+


The file retains the structure of the station metadata that was
previously reformmated by the :func:`gridwxcomp.prep_metadata` and :func:`gridwxcomp.ee_download.download_grid_data`, 
in that it each rows refers to a distinct weather station and any metadata 
that was in the original station metadata file created by the user 
is retained. There are four major variables calculated by :func:`gridwxcomp.calc_bias_ratios` 
that were added to this file, they are the long-term mean bias ratios (suffix
“\_mean”), the count of paired days used in those calculations (suffix
“\_count”), the standard deviation of the annual bias ratios or
differences (suffix “stdev”), and the coefficient of variation (suffix
“\_cv”).

At this step in the normal workflow of ``gridwxcomp`` the output file created by :func:`gridwxcomp.calc_bias_ratios` can be used for spatial mapping of point data and interpolation of the results using the :mod:`gridwxcomp.spatial` module.

    
    

References
----------

   .. [Allen2005] R.G. Allen, I.A. Walter, R.L. Elliott, T.A. Howell, D. Itenfisu, M.E. Jensen, and R.L. Snyder, The ASCE Standardized Reference Evapotranspiration Equation, American Society of Civil Engineers, 2005. https://doi.org/10.1061/9780784408056.
