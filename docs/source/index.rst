gridwxcomp - Grid Weather Station Comparison
===============================================

|Build| |Downloads per month| |PyPI version| 

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

