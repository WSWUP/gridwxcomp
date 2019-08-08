===============================================
gridwxcomp - gridMET Weather Station Comparison
===============================================

A package for comparing climate station time series data to `gridMET <http://www.climatologylab.org/gridmet.html>`_ (or other gridded) data. Major functionality includes tools to: 

* pair station locations with overlapping grid cells 
* download gridMET data from OpenDap or Google Earth Engine servers 
* calculate mean bias ratios between station and gridded data and other statistics 
* spatial interpolation of bias ratios with multiple interpolation options 
* build geo-referenced files of point data and interpolated rasters
* extract zonal means for a subset of gridcells  
* graphics, e.g. time series comparisons and interpolated residuals at stations 

Bias ratios calculated by ``gridwxcomp`` can be used to correct bias of grid to station data based on the properties of the stations. For example, monthly humidity ratios between station and grid for stations within agricultural settings can be used to estimate grid bias relative to agricultural locations. ``gridwxcomp`` includes an intuitive command line interface and a Python API.

.. toctree::
   :maxdepth: 2
   :caption: Contents:
   
   getting_started
   api
   changelog


