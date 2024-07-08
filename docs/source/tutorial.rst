Tutorial
========

This tutorial describes in detail how to use the ``gridwxcomp`` Python
package including preparation of input data, downloading gridded data
and pairing it with station weather data, calculating monthly bias
ratios between station and gridded data, spatial interpolation of point
results, and generating interactive graphics files.

.. code:: ipython3
    
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

.. code:: verbatim

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

.. code:: verbatim

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

.. code:: verbatim

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
‘Filename’, and ‘Station’.





References
----------

   .. [Allen2005] R.G. Allen, I.A. Walter, R.L. Elliott, T.A. Howell, D. Itenfisu, M.E. Jensen, and R.L. Snyder, The ASCE Standardized Reference Evapotranspiration Equation, American Society of Civil Engineers, 2005. https://doi.org/10.1061/9780784408056.
