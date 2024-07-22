
Automated testing 
=================

Software tests are automatically run each time a change to ``gridwxcomp`` is made on the master 
branch in GitHub using a `GitHub Actions 
workflow <https://github.com/WSWUP/gridwxcomp/actions/workflows/gridwxcomp_tests.yml>`__ 
that uses `pytest <https://docs.pytest.org/en/8.2.x/>`__.  Automated tests help spot potential 
bugs early so that they be identified and corrected efficiently.  

Running tests manually
^^^^^^^^^^^^^^^^^^^^^^

``pytest`` is required to run software tests that are provided. You can install ``pytest`` with PIP:

.. code-block:: bash

    $ pip install pytest

The tests utilize the example weather station data and configuration file provided with the software whether installed from PyPI or GitHub. Alternatively, these files can be found `here <https://github.com/Open-ET/flux-data-qaqc/tree/master/examples>`__.

.. note::
   The example files utilize the `CONUS404 gridded data <https://support.climateengine.org/article/117-conus404>`__ 
   which is hosted by OpenET on Google Earth Engine, 
   it is a public asset and as long as you have access to Google Earth Engine you should have no issues 
   accessing the data. The export path that is specified in the tests will be automatically created and 
   you must have authenticated and initialized Google Earth Engine for Python before running the tests. 

To run the tests, navivgate to the root directory of the source code (from the command line or shell) and run pytest:

.. code-block:: bash

    $ pytest
    
This will print out basic test results, usage of ``pytest`` plugins and command line options can be used for getting more information out of the tests.

