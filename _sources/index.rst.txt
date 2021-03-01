gridwxcomp - Grid Weather Station Comparison
===============================================

|Build| |Coverage| |Downloads per month| |PyPI version| 

-----------

A package for comparing climate station time series data to `gridMET <http://www.climatologylab.org/gridmet.html>`_ (or other gridded) data. Major functionality includes tools to: 

* pair station locations with overlapping grid cells 
* download gridMET data from OpenDap server 
* calculate mean bias ratios between station and gridded data and other statistics 
* spatial interpolation of bias ratios with multiple interpolation options 
* build geo-referenced files of point data and interpolated rasters
* extract zonal means for a subset of gridcells  
* graphics, e.g. time series comparisons and interpolated residuals at stations 
* tools for doing the same analysis with arbitrary gridded and point data

Bias ratios calculated by ``gridwxcomp`` can be used to correct bias of grid to station data based on the properties of the stations. For example, monthly humidity ratios between station and grid for stations within agricultural settings can be used to estimate grid bias relative to agricultural locations. ``gridwxcomp`` includes an intuitive command line interface and a Python API.

.. toctree::
   :maxdepth: 2
   :caption: Contents:
   
   getting_started
   faq
   api
   changelog


.. |Coverage| image:: https://coveralls.io/repos/github/WSWUP/gridwxcomp/badge.svg?branch=master
   :target: https://coveralls.io/github/WSWUP/gridwxcomp?branch=master

.. |Build| image:: https://travis-ci.org/WSWUP/gridwxcomp.svg?branch=master
   :target: https://travis-ci.org/WSWUP/gridwxcomp

.. |Downloads per month| image:: https://img.shields.io/pypi/dm/gridwxcomp.svg
   :target: https://pypi.python.org/pypi/gridwxcomp/

.. |PyPI version| image:: https://img.shields.io/pypi/v/gridwxcomp.svg
   :target: https://pypi.python.org/pypi/gridwxcomp/

