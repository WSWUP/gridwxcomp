gridwxcomp
==========

|Build| |Coverage| |Documentation Status| |Downloads per month| |PyPI version|

-----------

A package for comparing weather station data to gridded weather data that are hosted on Google Earth Engine. Major functionality includes: 

* pairing of station locations with overlapping grid cells 
* downloading point data from gridded datasets on Google Earth Engine 
* calculation of mean bias ratios between station and gridded data and other statistics 
* performing spatial interpolation of bias ratios with multiple options 
* building geo-referenced vector and raster data of spatially interpolated bias ratios and statistics
* extraction of zonal means from spatially interpolated bias ratios using the gridded dataset resolution 
* production of interactive time series and scatter plot comparisons between station and gridded data
* calculation and plotting of the residuals between spatially interpolated bias ratios and those computed at station locations 

Bias ratios calculated by ``gridwxcomp`` can be used to correct bias of grid to station data based on the properties of the stations. For example, monthly humidity ratios between station and grid for stations within agricultural settings can be used to estimate grid bias relative to agricultural locations. 

``gridwxcomp`` was used to create monthly bias ratios of `gridMET <http://www.climatologylab.org/gridmet.html>`_ reference evapotranspiration (ETo) data relative to ETo calculated at irrigated weather stations. The bias ratios were subsequently interpolated and used to correct gridMET ETo which is a key scaling flux for most of the remote sensing models that are part of the `OpenET <http://www.openetdata.org>`_ platform. 

Documentation
-------------
`Online documentation <https://gridwxcomp.readthedocs.io/en/latest/>`_

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




.. |Coverage| image:: https://coveralls.io/repos/github/WSWUP/gridwxcomp/badge.svg?branch=master&kill_cache=1
   :target: https://coveralls.io/github/WSWUP/gridwxcomp?branch=master&kill_cache=1

.. |Build| image:: https://github.com/WSWUP/gridwxcomp/actions/workflows/gridwxcomp_tests.yml/badge.svg
   :target: https://github.com/WSWUP/gridwxcomp/actions

.. |Downloads per month| image:: https://img.shields.io/pypi/dm/gridwxcomp.svg
   :target: https://pypi.python.org/pypi/gridwxcomp/

.. |Documentation Status| image:: https://img.shields.io/website-up-down-green-red/http/shields.io.svg
   :target: https://wswup.github.io/gridwxcomp/

.. |PyPI version| image:: https://img.shields.io/pypi/v/gridwxcomp.svg
   :target: https://pypi.python.org/pypi/gridwxcomp/
