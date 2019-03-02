gridwxcomp
==========

Station-based bias correction of gridded weather for agricultural applications

A Python package for calculating bias correction factors between climate stations and gridMET variables. Correction ratios can be used to correct gridMET estimated reference evapotranspiration (ETr) or other climatic variables to weather station observed data e.g. stations within agriculture settings. The package includes tools to pair stations locations with gridMET cells, download gridMET data using Google Earth Engine API, calculate point ratios, build geo-referenced files of point data,  conduct spatial interpolation of correction ratios with various interpolation options, and perform zonal extraction of mean correction ratios to a subset of gridMET cells around the stations. ``gridwxcomp`` includes an intuitive command line interface and a set of Python functions.

Documentation
-------------
A full documentation website is under development.

Installation
------------

We recommend installing ``gridwxcomp`` with `pip <https://pip.pypa.io/en/stable/installing/>`_, 

.. code-block:: bash

    $ pip install gridwxcomp 

Alternatively you can install manually. First navigate to the directory where you would like to use the software and then clone the repository using git from the command line:

.. code-block:: bash

    $ git clone https://github.com/WSWUP/gridwxcomp.git


If you do not have git, `get it here <https://git-scm.com/book/en/v2/Getting-Started-Installing-Git>`_. Alternatively, you can download the repository on GitHub or on `PyPI <https://pypi.org/project/gridwxcomp/>`_. Next, to install dependencies manually use `conda <https://conda.io/projects/conda/en/latest/user-guide/install/index.html>`_ and the provided virtual environment for ``gridwxcomp``. 

Once conda is installed move to the ``env`` directory and create the ``gridwxcomp`` virtual environment using the appropriate environment file included, e.g. on Windows:

.. code-block:: bash

    $ cd gridwxcomp/env
    $ conda env create -f env_windows.yml

To activate the environment whenever you need to use ``gridwxcomp`` just run

.. code-block:: bash

    $ activate gridwxcomp

or, on Linux

.. code-block:: bash

    $ source activate gridwxcomp

Lastly, ``gridwxcomp`` uses the Google Earth Engine API to download gridMET data, therefore you will need a Google account and before the first use on a machine you will need to verify your account. From the command line type:

.. code-block:: bash

    $ python -c "import ee; ee.Initialize()"

and follow the instructions.

If you installed manually, you can make all of ``gridwxcomp`` available on your environment variable PATH with `pip <https://pip.pypa.io/en/stable/installing/>`_. From the root directory of the cloned or downloaded copy of ``gridwxcomp`` run

.. code-block:: bash

    $ pip install --editable .

This will allow for running the ``gridwxcomp`` command line interface from any directory and aslo for importing ``gridwxcomp`` modules and functions in any Python environment on your system. However, for the example below, you will need to know where the example data is residing on your computer. 

Quick start from command line
-----------------------------

This workflow will use the example data given in "gridwxcomp/gridwxcomp/example_data" which includes four climate stations. To find the location of this data from the command line type

.. code-block:: bash

    $ python -c "import pkg_resources; print(pkg_resources.resource_filename('gridwxcomp', 'example_data/Station_Data.txt'))"

Once complete, this example workflow will calculate bias ratios between station and gridMET ETr, spatially interpolate GeoTIFF rasters of bias ratios at 400m resolution, and calculate zonal statistics of mean bias ratios for each gridMET cell in the region of the stations, similar to what is shown below.

.. image:: https://raw.githubusercontent.com/DRI-WSWUP/gridwxcomp/master/docs/source/_static/test_case.png?sanitize=true
   :align: center

The same workflow can be done on climate variables other than ETr using ``gridwxcomp``, e.g. observed ET, temperature, precipitation, wind speed, short wave radiation, etc.

After installing with pip the ``gridwxcomp`` command line interface can be used from any directory,

.. code-block:: bash

    $ gridwxcomp prep-input -i <PATH_TO example_data/Station_Data.txt>  

This will result in the file "merged_input.csv". Next download matching gridMET climate time series by running

.. code-block:: bash

    $ gridwxcomp download-gridmet-ee merged_input.csv -y 2016-2017

The time series of gridMET data that correpond with the stations in "merged_input.csv" will be saved to a new folder called "gridmet_data" by defualt. In this case only the years 2016-2017 are used because the test station data time coverage only includes recent years plus it saves time. 

Next, this command calculates monthly (and annual) bias ratios for each station/gridMET pair and saves them to CSV files 

.. code-block:: bash

    $ gridwxcomp calc-bias-ratios merged_input.csv -o monthly_ratios 

Last, to calculate interpolated surfaces of mean bias ratios and extract zonal means to gridMET cells using the default interpolation method (inverse distance weighting):

.. code-block:: bash

    $ gridwxcomp spatial monthly_ratios/etr_mm_summary_comp.csv -b 5

The ``[-b 5]`` option indicates that we want to expand the rectangular bounding area for interpolation by five gridMET cells (extrapolation in the outer regions).

The final output file "monthly_ratios/spatial/etr_mm_invdist_400m/gridMET_stats.csv" contains monthly bias ratios for each gridMET cell in the interpolation region, similar to what is shown below. 

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

Note ``GRIDMET_ID`` is the index of the master gridMET dataset 4 km fishnet grid starting at 0 in the upper left corner and moving across rows and down columns. This value can be joined with previously created data to relate the ID values to centroid locations of cells. 

GeoTIFF rasters of interpolated ratios will be saved to "monthly_ratios/spatial/etr_mm_invdist_400m/". Note, the gridMET variable name (etr_mm), the interpolation method (invdist), and the raster resolution (400m) are specified in the output directory. A fishnet grid with gridMET id values and a point shapefile of station ratios should all be created and saved in the "monthly_ratios/spatial/" directory.

To get help with any of the above ``gridwxcomp`` commands use the ``[--help]`` option, e.g.

.. code-block:: bash

    $ gridwxcomp spatial --help
