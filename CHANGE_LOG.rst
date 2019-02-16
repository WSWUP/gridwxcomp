Change log
**********

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

