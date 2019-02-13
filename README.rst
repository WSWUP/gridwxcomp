gridwxcomp
==========

Station-based bias correction of gridded weather for agricultural applications

A Python module that calculates monthly bias correction factors that can be utilized to bias correct gridMET estimated reference evapotranspiration (ETr) or other climatic variables to weather station observed data (e.g. stations within agriculture settings). The process includes spatial interpolation of correction ratios with various interpolation options and zonal extraction of mean correction ratios to gridMET cells. ``gridwxcomp`` includes a command line interface as well as a set of Python submodules.

Documentation
-------------
A full documentation website is under development.

Installation
------------

Currently ``gridwxcomp`` requires manual download for installation. However it will be available on the Python Package Index `PyPI <https://pypi.org/>`_ for installation with pip soon. 

To download first navigate to the directory where you would like to use the software and then clone the repository using git from the command line:

.. code-block:: bash

    $ git clone https://github.com/DRI-WSWUP/gridwxcomp.git


If you do not have git, `get it here <https://git-scm.com/book/en/v2/Getting-Started-Installing-Git>`_. Alternatively, you can download the repository using the download option near the top this page on GitHub.

Next, to install dependencies we recommend using `conda <https://conda.io/projects/conda/en/latest/user-guide/install/index.html>`_ and the provided virtual environment for ``gridwxcomp``. 

Once conda is installed move to the ``env`` directory and create the ``gridwxcomp`` virtual environment using the appropriate environment file included, e.g. on Windows:

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

To make all of ``gridwxcomp`` available on your environment variable PATH install it with `pip <https://pip.pypa.io/en/stable/installing/>`_. From the root directory of the cloned or downloaded copy of ``gridwxcomp`` run

.. code-block:: bash

    $ pip install --editable .

This will allow for importing ``gridwxcomp`` modules and functions in any Python environment on your system and also running the command line scripts, as shown below, from any directory. However you will need to know where the example data is residing on your computer if running from a different directory. 

Quick start from command line
-----------------------------

This workflow will use the example data given in "gridwxcomp/gridwxcomp/example_data" which includes four climate stations. It will calculate bias ratios between station and gridMET ETr, spatially interpolate GeoTIFF rasters of bias ratios at 400m resolution, and calculate zonal statistics of mean bias ratios for each gridMET cell in the region of the stations as shown below.

.. image:: https://raw.githubusercontent.com/DRI-WSWUP/gridwxcomp/master/docs/source/_static/test_case.png?sanitize=true
   :align: center

The same workflow can be done on climate variables other than ETr using ``gridwxcomp``, e.g. observed ET, temperature, precipitation, wind speed, short wave radiation, etc.

From ``gridwxcomp`` root directory run

.. code-block:: bash

    $ cd gridwxcomp
    $ python prep_input.py -i example_data/Station_Data.txt  

This will result in the file "merged_input.csv". Next download matching gridMET climate time series by running

.. code-block:: bash

    $ python download_gridmet_ee.py -i merged_input.csv -o test_gridmet_data -y 2016-2017

In this case the years 2016-2017 are used because the test station data time coverage only includes recent years plus it saves time as an example run by downloading a single year. Next to calculate monthly bias ratios and save to CSV files run

.. code-block:: bash

    $ python calc_bias_ratios.py -i merged_input.csv -o test_ratios 

Last, to calculate interpolated spatial surfaces of bias ratios and extract zonal means use the file produced from the previous step as input:

.. code-block:: bash

    $ python spatial.py -i test_ratios/etr_mm_summary_comp.csv -b 5

The ``[-b 5]`` option indicates that we want to expand the rectangular bounding area for interpolation by five gridMET cells (extrapolation in the outer regions).

The final output file "test_ratios/etr_mm_gridmet_summary_inverse_400m.csv" contains monthly bias ratios for each gridMET cell in the interpolation region, similar to what is shown below. 

    ========== ======== ======== ======== 
    GRIDMET_ID Apr_mean Aug_mean Dec_mean 
    ========== ======== ======== ======== 
    515902     0.75     0.68     0.63     
    514516     0.76     0.69     0.64     
    513130     0.77     0.69     0.65     
    511744     0.77     0.70     0.65     
    510358     0.78     0.70     0.66     
    ...        ...      ...      ...
    ========== ======== ======== ========

GeoTIFF rasters of interpolated ratios, a fishnet grid with gridMET id values, and a point shapefile of station ratios should all be created within the "test_ratios" directory.


