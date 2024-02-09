gridwxcomp
==========

|Build| |Coverage| |Documentation Status| |Downloads per month| |PyPI version|

-----------

A package for comparing weather station data to gridded weather data including `gridMET <http://www.climatologylab.org/gridmet.html>`_ and other gridded datasets that are hosted on Google Earth Engine. Major functionality includes: 

* pairing of station locations with overlapping grid cells 
* downloading point data from gridded datasets on Google Earth Engine or gridMET data from OpenDap server 
* calculation of mean bias ratios between station and gridded data and other statistics 
* performing spatial interpolation of bias ratios with multiple options 
* building geo-referenced vector and raster data of spatially interpolated bias ratios and statistics
* extraction of zonal means from spatially interpolated bias ratios using the gridded dataset resolution 
* production of interactive time series and scatter plot comparisons between station and gridded data
* calculation and plotting of the residuals between spatially interpolated bias ratios and those computed at station locations 

Bias ratios calculated by ``gridwxcomp`` can be used to correct bias of grid to station data based on the properties of the stations. For example, monthly humidity ratios between station and grid for stations within agricultural settings can be used to estimate grid bias relative to agricultural locations. ``gridwxcomp`` includes an intuitive command line interface and a Python API.

``gridwxcomp`` was used to create monthly bias ratios of `gridMET <http://www.climatologylab.org/gridmet.html>`_ reference evapotranspiration (ETo) data relative to ETo calculated at irrigated weather stations. The bias ratios were subsequently interpolated and used to correct gridMET ETo which is a key scaling flux for most of the remote sensing models that are part of the `OpenET <http://www.openetdata.org>`_ platform. 

Documentation
-------------
`Online documentation <https://wswup.github.io/gridwxcomp/>`_

Installation
------------

Currently we recommend using the provided conda environment file to install ``gridwxcomp`` and its dependencies in a virtual environment. Download the `environment.yml <https://raw.githubusercontent.com/WSWUP/gridwxcomp/master/gridwxcomp/env/environment.yml>`_ file and then install and activate it. If you don't have conda `get it here <https://conda.io/projects/conda/en/latest/user-guide/install/index.html>`_. To install dependencies in a virtual environment run 

.. code-block:: bash

    $ conda env create -f environment.yml

To activate the environment before using ``gridwxcomp`` run

.. code-block:: bash

    $ conda activate gridwxcomp

Optionally, install using `pip <https://pip.pypa.io/en/stable/installing/>`_,

.. code-block:: bash

    $ pip install gridwxcomp

Due to dependency conflicts you may have issues directly installing with pip before activating the conda environment.

Alternatively, or if there are installation issues, you can manually install. First activate the ``gridwxcomp`` conda environment (above). Next, clone or download the package from `GitHub <https://github.com/WSWUP/gridwxcomp>`_ or `PyPI <https://pypi.org/project/gridwxcomp/>`_ and then install locally with pip in "editable" mode. For example with cloning,

.. code-block:: bash

    $ git clone https://github.com/WSWUP/gridwxcomp.git
    $ cd gridwxcomp

If you are experiencing errors on installing the ``gridwxcomp`` conda environment above with dependencies. For example, if the Shapely package is not installing from the enironment.yml file, remove it or modify it from the "setup.py" file in the install requirements section before you install gridwxcomp from source with:

.. code-block:: bash

    $ pip install -e .

More help with installation issues related to dependency conflicts can be found in the ``gridwxcomp`` `issues <https://github.com/WSWUP/gridwxcomp/issues>`_ on GitHub, be sure to check the closed issues as well.


Quick start from command line
-----------------------------

This example uses data provided with ``gridwxcomp`` including climate variable time series data for four climate stations, it uses the gridMET as the gridded dataset however any uniform gridded data product can be used with ``gridwxcomp`` if extra information including a vector grid file is provided. 

After installation you can find the location of the data needed for the example by typing the following at the command line,

.. code-block:: bash

    $ python -c "import pkg_resources; print(pkg_resources.resource_filename('gridwxcomp', 'example_data/Station_Data.txt'))"

Once complete, this example will calculate bias ratios between station and gridMET ETr (reference evapotranspiration), spatially interpolate GeoTIFF rasters of bias ratios at 400 meter resolution, and calculate zonal statistics of mean bias ratios for each gridMET cell in the region of the stations, similar to what is shown in the figure below.

.. image:: https://raw.githubusercontent.com/WSWUP/gridwxcomp/master/docs/source/_static/test_case.png
   :align: center

The same procedure can be done for climate variables other than ETr, e.g. observed evapotranspiration, temperature, precipitation, wind speed, short wave radiation, etc.

After installing with pip the ``gridwxcomp`` command line interface can be used from any directory, the first step pairs climate station data with their nearest gridMET cell and produces a CSV file used in the following steps,

.. code-block:: bash

    $ gridwxcomp prep-input <PATH_TO example_data/Station_Data.txt>  

This will result in the file "merged_input.csv". Next download matching gridMET climate time series with `OpeNDAP <https://www.opendap.org>`_ by running

.. code-block:: bash

    $ gridwxcomp download-gridmet-opendap merged_input.csv -y 2016-2017

The time series of gridMET data that correpond with the stations in "merged_input.csv" will be saved to a new folder called "gridmet_data" by default. In this case only the years 2016-2017 are used. 

Next, to calculate mean monthly and annual bias ratios for each station/gridMET pair along with other statistics and metadata and save to CSV files, 

.. code-block:: bash

    $ gridwxcomp calc-bias-ratios merged_input.csv -o monthly_ratios 

Last, to calculate interpolated surfaces of mean bias ratios and extract zonal means to gridMET cells using the default interpolation method (inverse distance weighting):

.. code-block:: bash

    $ gridwxcomp spatial monthly_ratios/etr_mm_summary_comp_all_yrs.csv -b 5

The ``[-b 5]`` option indicates that we want to expand the rectangular bounding area for interpolation by five gridMET cells (extrapolation in the outer regions).

GeoTIFF rasters of interpolated ratios will be saved to "monthly_ratios/spatial/etr_mm_invdist_400m/". Note, the gridMET variable name (etr_mm), the interpolation method (invdist), and the raster resolution (400m) are specified in the output directory. A fishnet grid with gridMET id values and a point shapefile of station ratios should all be created and saved in the "monthly_ratios/spatial/" directory.

The output file "monthly_ratios/spatial/etr_mm_invdist_400m/gridMET_stats.csv" contains monthly bias ratios for each gridMET cell in the interpolation region, similar to what is shown below. 

    ========== ======== ======== ======== ===
    GRIDMET_ID Jan_mean Feb_mean Mar_mean ...
    ========== ======== ======== ======== ===
    515902     0.66     0.76     0.96     ...
    514516     0.66     0.77     0.96     ...
    513130     0.67     0.77     0.97     ...
    511744     0.67     0.78     0.97     ...
    510358     0.68     0.79     0.97     ...
    ...        ...      ...      ...      ...
    ========== ======== ======== ======== ===

Note ``GRIDMET_ID`` is the index of the master gridMET dataset 4 km fishnet grid starting at 0 in the upper left corner and moving across rows and down columns. This value can be joined with previously created data, e.g. the ID values can be joined to centroid coordinates of gridMET cells. 

Bar plots that show the residual between station mean ratios and interpolated estimates are saved to "monthly_ratios/spatial/etr_mm_invdist_400m/residual_plots/".

To get abbreviated descriptions for any of the above ``gridwxcomp`` commands use the ``[--help]`` option, e.g.

.. code-block:: bash

    $ gridwxcomp spatial --help



.. |Coverage| image:: https://coveralls.io/repos/github/WSWUP/gridwxcomp/badge.svg?branch=master&kill_cache=1
   :target: https://coveralls.io/github/WSWUP/gridwxcomp?branch=master&kill_cache=1

.. |Build| image:: https://travis-ci.org/WSWUP/gridwxcomp.svg?branch=master
   :target: https://travis-ci.org/WSWUP/gridwxcomp

.. |Downloads per month| image:: https://img.shields.io/pypi/dm/gridwxcomp.svg
   :target: https://pypi.python.org/pypi/gridwxcomp/

.. |Documentation Status| image:: https://img.shields.io/website-up-down-green-red/http/shields.io.svg
   :target: https://wswup.github.io/gridwxcomp/

.. |PyPI version| image:: https://img.shields.io/pypi/v/gridwxcomp.svg
   :target: https://pypi.python.org/pypi/gridwxcomp/
