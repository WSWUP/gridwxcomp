Change log
**********

Version 0.0.4
=============

new function: ``spatial.calc_pt_error`` which
* calculates interpolated point ratios and residuals betwen station data
  * updates summary CSV and point shapefile, copies to directory with rasters

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

