gridwxcomp
==========

A Python package for calculating bias correction factors between climate stations and gridMET variables. Correction ratios can be used to adjust gridMET estimated reference evapotranspiration (ETr) or other climatic variables to weather station data e.g. stations within agriculture settings. The package includes tools to pair stations locations with gridMET cells, download gridMET data using Google Earth Engine API, calculate point ratios, build geo-referenced files of point data, conduct spatial interpolation of correction ratios with various interpolation options, and perform zonal extraction of mean correction ratios to a subset of gridMET cells around the stations. gridwxcomp includes an intuitive command line interface and a set of Python functions.

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
