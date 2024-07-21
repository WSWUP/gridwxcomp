gridwxcomp - Grid Weather Station Comparison
===============================================

|Build| |Downloads per month| |PyPI version| 

-----------

A package for comparing weather station data to gridded weather data that are hosted on Google Earth Engine. Major functionality includes: 

* parsing of multiple weather stations and weather variables and metadata
* downloading point data from gridded datasets on Google Earth Engine at weather station locations 
* temporal pairing of station and gridded data
* unit handling and automated conversions
* calculation of mean bias ratios between station and gridded data and related statistics 
* performing spatial mapping and interpolation of bias ratios with multiple options 
* calculation of residuals between spatially interpolated bias ratios and those computed at station locations 
* building geo-referenced vector and raster data of spatially interpolated and point data
* zonal averaging of spatially interpolated bias results using a fishnet grid  
* interactive graphics (time series, scatter, and bar charts) comparing station and gridded data


Bias ratios calculated by ``gridwxcomp`` can be used to correct bias of grid to station data based on the properties of the stations. For example, monthly humidity ratios between station and grid for stations within agricultural settings can be used to estimate grid bias relative to agricultural locations. 

.. toctree::
   :maxdepth: 3
   :caption: Contents:
   
   install
   tutorial
   faq
   api
   changelog


.. |Build| image:: https://github.com/WSWUP/gridwxcomp/actions/workflows/gridwxcomp_tests.yml/badge.svg
   :target: https://github.com/WSWUP/gridwxcomp/actions

.. |Downloads per month| image:: https://img.shields.io/pypi/dm/gridwxcomp.svg
   :target: https://pypi.python.org/pypi/gridwxcomp/

.. |PyPI version| image:: https://img.shields.io/pypi/v/gridwxcomp.svg
   :target: https://pypi.python.org/pypi/gridwxcomp/

