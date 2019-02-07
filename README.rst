gridwxcomp
==========

Station-based bias correction of gridded weather for agricultural applications

A Python module that calculates monthly bias correction factors that can be utilized to bias correct gridMET estimated ETr or other climatic variables to weather station observed data (e.g. stations within agriculture settings). The process includes spatial interpolation of correction ratios with various interpolation options and zonal extraction of mean correction ratios to gridMET cells. ``gridwxcomp`` includes a command line interface as well as a set of Python submodules.

Documentation
-------------
A full documentation website is under development.

Installation
------------

Currently ``gridwxcomp`` requires manual download for installation. However it will be available on on the Python Package Index `PyPI <https://pypi.org/>`_ for installation with pip soon. 

To download first navigate to the directory where you would like to use the software and then clone the repository using git from the command line:

.. code-block:: bash

    git clone https://github.com/DRI-WSWUP/gridwxcomp.git


If you do not have git, `get it here <https://git-scm.com/book/en/v2/Getting-Started-Installing-Git>`_. Alternatively, you can download the repository using the download option near the top this page on GitHub.

Next, to install dependencies you will need `conda <https://conda.io/projects/conda/en/latest/user-guide/install/index.html>`_. 

Once conda is installed move to the ``env`` directory and create the ``gridwxcomp`` environment using the appropriate environment file included, e.g. on Windows:

.. code-block:: bash

    $ cd env
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

Quick start from command line
-----------------------------

This workflow will use the example station time series data given in "gridwxcomp/gridwxcomp/ETrBias_DataPackage". It will calculate bias ratios between station and gridMET ETr, estimate spatial GeoTIFF rasters at 400m resolution, and calculate zonal statistics of mean bias ratios for each gridMET cell in the region of the stations as shown below.

.. image:: https://raw.githubusercontent.com/DRI-WSWUP/gridwxcomp/master/docs/source/_static/test_case.png?sanitize=true
   :align: center

From ``gridwxcomp`` root directory run

.. code-block:: bash

    $ cd gridwxcomp
    $ python prep_input.py -i ETrBias_DataPackage/Station_Data.txt  

This will result in the file "merged_input.csv", next to download matching gridMET climate time series run

.. code-block:: bash

    $ python download_gridmet_ee.py -i merged_input.csv -o test_gridmet_data -y 2016-2017

In this case the years 2016-2017 are used because the test data is short plus it saves time by downloading a single year. Next to calculate monthly bias ratios and save to CSV files run

.. code-block:: bash

    $ python calc_bias_ratios.py -i merged_input.csv -o test_ratios -c

Last, to calculate interpolated spatial surfaces of bias ratios and extract zonal means:

.. code-block:: bash

    $ python spatial.py -i test_ratios/summary_comp.csv -b 5

The ``[-b 5]`` option indicates to exapnd the rectangular bounding area for interpolation by five gridMET cells.

The final output files including monthly bias ratios for each gridMET cell, GeoTIFF rasters of interpolated ratios, a fishnet grid with gridMET id values, and a point shapefile of station ratios should all be created within the "test_ratios" directory.
