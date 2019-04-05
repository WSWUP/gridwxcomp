Change Log
**********

Version 0.0.6
=============

Update all ``gridwxcomp`` modules (with exception of ``interpgdal.py``) to ensure that they can be used as standalone scripts without installing ``gridwxcomp``. For example, the spatial interpolation routines in ``spatial.py`` can be used from the command line

.. code-block:: sh

    $ python spatial.py <options>

The only requirement is that all the files within ``gridwxcomp/gridwxcomp`` are within the same directory. 

Add ``gridwxcomp.plot`` module to consolidate current and future graphics tools. Current tools include the following functions: ``gridwxcomp.plot.daily_comparison``, ``gridwxcomp.plot.monthly_comparison``, and the newly added ``gridwxcomp.plot.station_bar_plot``. Changed the ``gridwxcomp`` command line interface plot command to handle the three options using the new option ``[-t, --plot-type]``. 

Add year range option for ``gridwxcomp.plot.daily_comparison``, useful for adding additional data to scatter plot comparisons using multiple years data for a particular month. 

Changed docs hosting to `GitHub <https://wswup.github.io/gridwxcomp/>`_

Version 0.0.5
=============

Functionality for climate station data that was **NOT** created by `PyWeatherQaQc <https://github.com/WSWUP/pyWeatherQAQC>`_ after ``gridwxcomp >= 0.0.55``. Climate station time series files should be in CSV format and need a "date" column with date strings that can be parsed as datetime objects, e.g. '12/01/2018' or '12-01-2018'. ``daily_comparison.py`` and ``monthly_comparison.py`` plotting modules however still require climate station input data in the format of ``PyWeatherQaQc``. 

Add monthly plotting to command line interface, change command line command "daily-comparison" to
"plot" with "daily" and "monthly" options. 

Add documentation page at `ReadTheDocs <http://gridwxcomp.readthedocs.io/>`_

Add option to re-download gridMET time series data using ``download_gridmet_ee`` for specified year range.

Version 0.0.4
=============

Improve handling of missing data, if ratio data is missing it is generally represented by ``-999`` in text files (i.e. CSVs) and by ``nan`` in geospatial files, e.g. within point shapefiles of bias ratios. Importantly, fixed bug where gdal interpolation methods used missing data in interpolation as zeros.

Add ``util.py`` module to hold utility functions or classes which may be useful to multiple ``gridwxcomp`` modules. Added function to index a pandas DataFrame or Series that has a datetime index based on a user input year filter.

Add year filter option to ``calc_bias_ratios.py``, so that certain years or ranges of years are only used to calculate bias ratios and statistics, the file names of the summary CSV files are also modified with the year or range added as a suffix so that they can be distinguished and used for spatial interpolation. 

New function: ``spatial.calc_pt_error`` which
* calculates interpolated point ratios and residuals betwen station data
* updates summary CSV and point shapefile, copies to directory with rasters

For example, now after building point shapefile, making the extraction grid, and interpolating point bias ratios using the ``spatial`` module with default options but only interpolating two layers, the following file structure is created from the root directory holding the ratio sumary CSVs::

    .
    ├── etr_mm_summary_comp.csv
    ├── etr_mm_summary.csv
    └── spatial
        ├── etr_mm_invdist_400m
        │   ├── annual_mean.tiff
        │   ├── annual_mean.vrt
        │   ├── etr_mm_summary_comp.csv
        │   ├── etr_mm_summary_pts.cpg
        │   ├── etr_mm_summary_pts.dbf
        │   ├── etr_mm_summary_pts.prj
        │   ├── etr_mm_summary_pts.shp
        │   ├── etr_mm_summary_pts.shx
        │   ├── gridMET_stats.csv
        │   ├── growseason_mean.tiff
        │   └── growseason_mean.vrt
        ├── grid.cpg
        ├── grid.dbf
        ├── grid.prj
        ├── grid.shp
        └── grid.shx

Note, now there is a copy of the summary_comp.csv file in the directory containing the interpolated rasters, and the point shapefile is also saved there as opposed to the "spatial" dir in previous versions. The CSV in the root directory is needed for running additional interpolations, the copy also contains newly added interpolation estimates at points and error residuals which are unique to a specific interpolation run.

Change calculations of annual, growing season, and summer bias ratios to use period sum of data as opposed to mean of monthly ratios. Same for standard deviation calculations and coefficient of variation. Results in slightly more accurate values. Also add total day accounts for these time periods, add all of these fields to georeferenced point shapefile as opposed to only bias ratios in previous versions.

Version 0.0.3
=============

First version available on `PyPI <https://pypi.org/project/gridwxcomp/>`_.

Add class ``gridwxcomp.interpgdal.InterpGdal`` for interpolation methods provided by the `gdal_grid <https://www.gdal.org/gdal_grid.html>`_ command, the most useful being inverse distance weighting to a power and inverse distance weighting to a power with n nearest neighbors. The ``InterpGdal`` object can be used on its own within Python to efficiently produce interpolated rasters of arbitrary variables from point data that is calculated by ``gridwxcomp.calc_bias_ratios``, it is also used in the main spatial interpolation workflow, e.g. the command line usage of ``gridwxcomp.spatial``, by providing additional interpolation routines in addition to the radial basis functions. Instance attributes allow for managing metadata of different interpolation outcomes such as parameter values and paths to output files.  

Added calculation of standard deviation and coefficient of variation for bias ratios to the ``gridwxcomp.calc_bias_ratios`` function.

Update file structure format for spatial interpolation and calculation of zonal statistics to gridMET cells. In previous versions a CSV file containing zonal statistics for gridMET cells was created based on the interpolation method, gridMET variable name, and interpolated raster resolution, e.g.::

        'etr_mm_gridmet_summary_linear_400m.csv'

which was saved to the output directory of ``calc_bias_ratios``, i.e. where the CSV file containing station point ratios and other statistics exists. This was problematic for tracking results created by multiple interpolation parameters such as changing the power parameter of the inverse distance weighting algorithm. So the new structure is saving a file named 'gridMET_stats.csv' to the output directory where interpolated rasters are saved for any interpolation routine, which can now be modified when conducting any interpolation. The columns in the CSV are updatedwhen layers are interpolated and zonal stats are extracted with the same out directory specified. 


Version 0.0.2
=============

Add more robust and intuitive command line interface ``gridwxcomp`` which interfaces with all major workflows of the module as opposed to needing to access multiple submodules of ``gridwxcomp``, e.g. ``gridwxcomp.prep_input``. Also add changelog. Example use of new CLI

.. code-block:: bash

    $ gridwxcomp prep-input <station_metadata_file>

old method (still possible if ``prep_input.py`` in working directory),

.. code-block:: bash

    $ python prep-input.py -i <station_metadata_file>

Added dependencies:

* `click >= 7.0 <https://click.palletsprojects.com/en/7.x/>`_

Version 0.0.1
=============

First numbered version. Many changes occured for initial development under this version which were not released or registered to PyPI. Main workflow has beed tested on Linux and Windows including: 

* pairing climate stations with gridMET cells
* calculation of bias correction ratios of climatic variables 
* created georeferenced point shapefiles, fishnet grid 
* perform 2-D interpolation of bias ratio surface with multiple options
* exctract zonal statistics to gridMET cells of bias ratio surface
* produce interactive plots comparing time series of station and gridMET data

Package not yet hosted on PyPI however it is packaged and can be installed to the Python and system env PATHs with 

.. code-block:: bash

    $ pip install --editable .

