gridwxcomp
==========

Station-based bias correction of gridded weather for agricultural applications

A package for comparing climate station time series data to `gridMET <http://www.climatologylab.org/gridmet.html>`_ data. Functionality includes tools to: pair station locations with overlapping gridMET cells; download gridMET data using the Google Earth Engine API; calculate mean bias ratios between station and gridMET data and other statistics; build geo-referenced files of point data; conduct spatial interpolation of correction ratios with multiple interpolation options; extract zonal mean data to a subset of gridMET cells; and produce graphical plots of time series comparisons. Bias ratios calculated by ``gridwxcomp`` can be used to correct bias of gridMET to station data based on the properties of the stations. For example, monthly humidity ratios between station and gridMET for stations within agricultural settings can be used to estimate gridMET bias relative to agricultural locations. gridwxcomp includes an intuitive command line interface and a Python API.

Documentation
-------------
At `ReadTheDocs <http://gridwxcomp.readthedocs.io/>`_

Installation
------------

Currently we recommend using the provided conda environment file to install dependencies in a virtual environment. Download the `environment.yml <https://raw.githubusercontent.com/WSWUP/gridwxcomp/master/gridwxcomp/env/environment.yml>`_ file and then install and activate it. If you don't have conda `get it here <https://conda.io/projects/conda/en/latest/user-guide/install/index.html>`_. To install dependencies in a virtual environment run 

.. code-block:: bash

    $ conda env create -f environment.yml

To activate the environment before using ``gridwxcomp`` run

.. code-block:: bash

    $ activate gridwxcomp

on windows, or on Linux, Mac

.. code-block:: bash

    $ source activate gridwxcomp

Next install using `pip <https://pip.pypa.io/en/stable/installing/>`_,

.. code-block:: bash

    $ pip install gridwxcomp

Due to dependency conflicts you may have issues directly installing with pip before activating the conda environment.

Alternatively, or if there are installation issues, you can manually install. First activate the ``gridwxcomp`` conda environment (above). Next, clone or download the package from `GitHub <https://github.com/WSWUP/gridwxcomp>`_ or `PyPI <https://pypi.org/project/gridwxcomp/>`_ and then install locally with pip in "editable" mode. For example with cloning,

.. code-block:: bash

    $ git clone https://github.com/WSWUP/gridwxcomp.git
    $ cd gridwxcomp
    $ pip install -e .

If you downloaded the source distribution then run ``pip install -e .`` in the root directory where the setup.py file is located. This installation method is ideal if you want to be able to modify the source code.

Lastly, ``gridwxcomp`` uses the Google Earth Engine API to download gridMET data, therefore you will need a Google account and before the first use on a machine you will need to verify your account. From the command line type:

.. code-block:: bash

    $ python -c "import ee; ee.Initialize()"

and follow the instructions.


Quick start from command line
-----------------------------

This example workflow uses data provided with ``gridwxcomp`` including climate variable time series data for four climate stations. After installation you can find the location of the data needed for the example by typing the following at the command line,

.. code-block:: bash

    $ python -c "import pkg_resources; print(pkg_resources.resource_filename('gridwxcomp', 'example_data/Station_Data.txt'))"

Once complete, this example will calculate bias ratios between station and gridMET ETr (reference evapotranspiration), spatially interpolate GeoTIFF rasters of bias ratios at 400m resolution, and calculate zonal statistics of mean bias ratios for each gridMET cell in the region of the stations, similar to what is shown in the figure below.

.. image:: https://raw.githubusercontent.com/WSWUP/gridwxcomp/master/docs/source/_static/test_case.png
   :align: center

The same procedure can be done for climate variables other than ETr, e.g. observed evapotranspiration, temperature, precipitation, wind speed, short wave radiation, etc.

After installing with pip the ``gridwxcomp`` command line interface can be used from any directory, the first step pairs climate station data with their nearest gridMET cell and produces a CSV file used in the following steps,

.. code-block:: bash

    $ gridwxcomp prep-input <PATH_TO example_data/Station_Data.txt>  

This will result in the file "merged_input.csv". Next download matching gridMET climate time series with Google Earth Engine by running

.. code-block:: bash

    $ gridwxcomp download-gridmet-ee merged_input.csv -y 2016-2017

The time series of gridMET data that correpond with the stations in "merged_input.csv" will be saved to a new folder called "gridmet_data" by defualt. In this case only the years 2016-2017 are used. 

Next, to calculate mean monthly and annual bias ratios for each station/gridMET pair along with other statistics and metadata and save to CSV files, 

.. code-block:: bash

    $ gridwxcomp calc-bias-ratios merged_input.csv -o monthly_ratios 

Last, to calculate interpolated surfaces of mean bias ratios and extract zonal means to gridMET cells using the default interpolation method (inverse distance weighting):

.. code-block:: bash

    $ gridwxcomp spatial monthly_ratios/etr_mm_summary_comp_all_yrs.csv -b 5

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

Note ``GRIDMET_ID`` is the index of the master gridMET dataset 4 km fishnet grid starting at 0 in the upper left corner and moving across rows and down columns. This value can be joined with previously created data, e.g. the ID values can be joined to centroid coordinates of gridMET cells. 

GeoTIFF rasters of interpolated ratios will be saved to "monthly_ratios/spatial/etr_mm_invdist_400m/". Note, the gridMET variable name (etr_mm), the interpolation method (invdist), and the raster resolution (400m) are specified in the output directory. A fishnet grid with gridMET id values and a point shapefile of station ratios should all be created and saved in the "monthly_ratios/spatial/" directory.

To get help with any of the above ``gridwxcomp`` commands use the ``[--help]`` option, e.g.

.. code-block:: bash

    $ gridwxcomp spatial --help
