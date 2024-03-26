---
title: 'gridwxcomp: A Python package to evaluate and interpolate biases between station and gridded weather data.'
tags:
  - Python
  - hydrology
  - interpolation
  - weather station
  - meteorology
  - Google Earth Engine
authors:
  - name: John M. Volk^[corresponding author]
    orcid: 0000-0001-9994-1545
    affiliation: 1
  - name: Christian Dunkerly
    orcid: 0000-0003-3592-4118
    affiliation: 1
  - name: Sayantan Majumdar 
    orcid: 0000-0002-3539-0147
    affiliation: 1
  - name: Christopher Pearson
    affiliation: 1 
  - name: Charles G. Morton
    affiliation: 1
  - name: Matt Bromley
    affiliation: 1
  - name: Justin L. Huntington
    affiliation: 1 
affiliations:
  - name: Desert Research Institute
    index: 1
date: January 2024
bibliography: paper.bib
---

# Introduction

Gridded weather data has become increasingly accessible and accurate over the recent decades and such data enables a variety of applications and research that require spatially continuous and high spatio-temporal resolution data [@Thornton2021]. Gridded weather data are often developed from an assimilation of data production and measurement techniques including the incorporation of in situ observational data networks, land surface modeling techniques, and remote sensing techniques. Although well-curated in situ measurements of weather variables are typically more accurate than gridded products, they provide data at a single point in space and they involve difficulties and expense related to deployment and sensor calibration and maintenance, resulting in incomplete spatial and temporal coverage. Gaps in spatial and temporal coverage in in situ weather data are often filled by gridded data products, however, the increased coverage provided by gridded data comes with the tradeoff of increased uncertainty and potential for bias [@Blankenau2020] that is introduced in the modeling and statistical data assimilation techniques used for gridded data production that can be difficult to characterize and quantify. 

# Statement of Need

Commonly, in situ measurements of weather data are used to validate and assess the uncertainty in their gridded counterparts. Direct comparisons between point measurements of in situ measured weather data and gridded data are useful to assess biases and uncertainties in gridded data. Point comparisons can be interpolated to investigate spatial biases given sufficient density of measurement stations. Finally, maps of spatial bias can subsequently be used to adjust or correct the gridded weather data for the observed bias. ``gridwxcomp`` was developed precisely to streamline these objectives in a reproducible Python framework. 

This package has functionality to download point data from a variety of gridded meteorological datasets that are hosted on [Google Earth Engine](https://developers.google.com/earth-engine/datasets/) (e.g., NLDAS, ERA5, gridMET) and pair those with station data, make comparison plots, monthly bias ratios and metrics, and interpolate those data to make spatially complete georeferenced raster images of bias between the gridded and station data using multiple interpolation techniques such as inverse distance weighting and linear interpolation. As far as the authors know, this is the only open-source software that accomplishes these tasks.

# Design and Features

``gridwxcomp`` is a Python 3 package that consists of five core submodules and two utility submodules (\autoref{fig:fig1}). The software can be run in Python environment or using the command line interface provided. Currently, ``gridwxcomp`` can process the following meteorological variables: air temperature (minimum and maximum), dew point temperature, precipitation, shortwave radiation, wind speed, vapor pressure, relative humidity (minimum, maximum, and average), and grass (short) and alfalfa (tall) reference evapotranspiration (ET). The following gridded weather datasets are available to be compared to station data using this package: CONUS404, ERA5-Land, gridMET, NLDAS, RTMA, and spatial CIMIS.

![Flowchart diagram of submodules and data processing pipeline of ``gridwxcomp``.\label{fig:fig1}](figure1.pdf)

Starting with ``prep_input``, this module parses metadata regarding one or more meteorological stations’ daily weather data and regarding the gridded dataset to compare with, this script outputs a CSV text file with coordinates used for pairing each meteorological station with the gridded data. 

The output file from ``prep_input`` is used by the next submodule ``ee_download`` which queries the Google Earth Engine for the previously specified weather data variables and gridded product and then downloads those data at gridcell centroids that are nearest to the corresponding station coordinates for the period that matches the station data. This process also updates the output file from ``prep_input`` to include the file paths of the downloaded gridded timeseries data. 

The ``calc_bias_ratios`` submodule reads in the output from the previous submodules and pairs the station and gridded daily data, aggregated the data to monthly periods, and computes average monthly, seasonal, and annual bias ratios or differences (for temperature variables) between the station and gridded data for the specified variable(s). In addition to ratios and differences, the interannual variability (standard deviation and coefficient of variation) of the results are also calculated as well as the number of paired data used in each calculated ratio/difference.  There are currently two options to calculate the long-term average bias ratios between station gridded weather data, the results are saved to CSV text files.

After computing station-gridded bias metrics using the ``calc_bias_ratios`` submodule, the ``spatial`` submodule offers automated tools to conduct spatial interpolation of the point biases and variability metrics and other spatial mapping functions. The ``spatial`` submodule has tools for developing a regular grid or fishnet polygon file that matches the domain of the area containing the weather stations of interest and that is at the same resolution of the gridded dataset, this grid can later be used for conducting zonal statistics from the interpolated bias surfaces. Point bias results that were previously computed in the processing pipeline are assigned a spatial reference and saved as a point polygon (shapefile). The main function of the ``spatial`` submodule is the interpolation of the point bias results, options include radial basis functions from the [SciPy Python library](https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.RBFInterpolator.html#scipy.interpolate.RBFInterpolator) and interpolation algorithms from the GDAL’s [gdal_grid](https://www.gdal.org/gdal_grid.html). The GDAL methods are contained in the ``interpgdal`` submodule for improved code organization and structure. Lastly, methods for computing point residuals between points containing measurements and the interpolated surfaces are provided in the ``spatial`` submodule, as well as methods for extracting zonal statistics from the grid that covers the interpolation domain. 

The ``plot`` submodule in ``gridwxcomp`` provides tools for making interactive diagnostic plots of both daily and monthly data, pairing the station data alongside its gridded counterpart. The plots are created using the [Bokeh](https://docs.bokeh.org/en/latest/index.html) Python package and are saved as HTML files that give the viewer the ability to pan and zoom into individual subplots. 

# Research Enabled by gridwxcomp

The most significant application of  ``gridwxcomp`` was the development of bias correction surfaces that are applied to gridded reference evapotransipiration (ET) data which are key inputs to some of the remote sensing ET models that comprise the OpenET platform [@Melton2021]. Daily from approximately 800 weather stations located in irrigated agricultural plots were curated and the American Society of Civil Engineers (ASCE) standardized Penman-Monteith reference ET equation [@allen2005] was used to estimate reference ET at the stations, then ``gridwxcomp`` was used to pair these data with the equivalent reference ET data from the gridMET [@Abatzoglou2013] dataset at their nearest gridcells and for temporally consistent periods. The long-term average monthly ratios for station reference ET to gridded reference ET were then calculated for each point and saved as georeferenced data by ``gridwxcomp`` and were subsequently spatially interpolated using a kriging approach. The final interpolated monthly surfaces are used within the OpenET platform to correct daily and monthly reference ET before it is used by the remote sensing ET models. 

# Acknowledgments

We would like to thank the Bureau of Reclamation, NASA Applied Sciences Program, and the Western States Water Use Program at the Desert Research Institute for providing funding for the development of this software.

# References

