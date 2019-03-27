===============================================
gridwxcomp - gridMET Weather Station Comparison
===============================================

A package for comparing climate station time series data to `gridMET <http://www.climatologylab.org/gridmet.html>`_ data. Major functionality includes tools to: 

* pair station locations with overlapping gridMET cells 
* download gridMET data using Google Earth Engine 
* calculate mean bias ratios between station and gridMET data and other statistics 
* spatial interpolation of bias ratios with multiple interpolation options 
* build geo-referenced files of point data and interpolated rasters
* extract zonal means for a subset of gridMET cells  
* graphics, e.g. time series comparisons and interpolated residuals at stations 

Bias ratios calculated by ``gridwxcomp`` can be used to correct bias of gridMET to station data based on the properties of the stations. For example, monthly humidity ratios between station and gridMET for stations within agricultural settings can be used to estimate gridMET bias relative to agricultural locations. ``gridwxcomp`` includes an intuitive command line interface and a Python API.

.. toctree::
   :maxdepth: 2
   :caption: Contents:
   
   getting_started
   api
   changelog

Indices and tables
==================

* :ref:`genindex`
* :ref:`search`
