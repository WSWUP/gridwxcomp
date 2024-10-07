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
  - name: Christopher Pearson
    affiliation: 1 
  - name: Charles G. Morton
    affiliation: 1
  - name: Justin L. Huntington
    affiliation: 1 
affiliations:
  - name: Desert Research Institute
    index: 1
date: April 2024
bibliography: paper.bib
---

# Introduction

Gridded weather data has become increasingly accessible and accurate over the recent decades, and such data enables a variety of applications and research that require spatially continuous and high spatio-temporal resolution data [@Thornton2021;@Rasmussen2023;@MuozSabater2021]. Gridded weather data are often developed from an assimilation of data production and measurement techniques, including the incorporation of in situ observational data networks, land surface modeling techniques, and remote sensing techniques. Although well-curated in situ measurements of weather variables are typically more accurate than gridded products, they provide data at a single point in space, and they involve difficulties and expenses related to deployment and sensor calibration and maintenance, resulting in incomplete spatial and temporal coverage. Gaps in spatial and temporal coverage in in situ weather data are often filled by gridded data products, however, the increased coverage provided by gridded data comes with the tradeoff of increased uncertainty and potential for bias [@Blankenau2020] that is introduced in the modeling and statistical data assimilation techniques used for gridded data production that can be difficult to characterize and quantify. 

# Statement of Need

Commonly, in situ measurements of weather data are used to validate and assess the bias and uncertainty in their gridded counterparts. Point biases can be interpolated to investigate spatial biases given sufficient density of measurement stations. Maps of spatial bias can subsequently be used to adjust the gridded weather data for the observed bias. ``gridwxcomp`` was developed to streamline these objectives in a reproducible Python framework. 

This package has the functionality to download point data from a variety of gridded meteorological datasets that are hosted on [Google Earth Engine](https://developers.google.com/earth-engine/datasets/) (e.g., NLDAS, ERA5, gridMET) and pair those with station data. It also has functionality to make comparison plots, calculate monthly bias ratios and metrics, and interpolate those data to make spatially complete georeferenced raster images of bias between the gridded and station data using multiple interpolation techniques such as inverse distance weighting and linear interpolation. As far as the authors know, this is the only open-source software that accomplishes these tasks.

# Design and Features

``gridwxcomp`` is a Python 3 package that consists of five core submodules and two utility submodules (\autoref{fig:fig1}). ``gridwxcomp`` can process the following meteorological variables: air temperature (minimum and maximum), dew point temperature, shortwave radiation, wind speed, vapor pressure, relative humidity (minimum, maximum, and average), and grass (short) and alfalfa (tall) reference evapotranspiration (ET). Daily gridded weather datasets hosted on Google Earth Engine can be accessed and compared against station data using ``gridwxcomp`` as long as the user has access to the dataset collection. Example public datasets include: CONUS404 [@Rasmussen2023], ERA5-Land [@MuozSabater2021], gridMET [@Abatzoglou2013], NLDAS [@Mitchell2004], RTMA [@DePondeca2011], and spatial CIMIS [@Hart2009]. 

![Flowchart diagram of submodules and data processing pipeline of ``gridwxcomp``.\label{fig:fig1}](figure1.pdf)

The ``prep_metadata`` submodule parses metadata of meteorological stations, including reprojection of station coordinates. The output file from ``prep_metadata`` is used by the ``ee_download`` submodule which queries Google Earth Engine for the gridded weather data specifed by the user and downloads time series data at the corresponding station coordinates. 

The ``calc_bias_ratios`` submodule pairs the station and gridded daily data, performs unit conversions, and computes average monthly, seasonal, and annual bias ratios or differences (for temperature variables) between the station and gridded data for a specified variable. In addition to computing station-to-gridded biases, the ``calc_bias_ratios`` routine calculates the interannual variability (standard deviation and coefficient of variation) of those metrics and the number of paired data used in each metric. 

After computing station-gridded bias metrics, the ``spatial`` submodule offers tools to conduct spatial interpolation of the point biases and variability metrics and other spatial mapping functions. The main function of the ``spatial`` submodule is the interpolation of the point bias results; options include the algorithms from GDAL’s [gdal_grid](https://www.gdal.org/gdal_grid.html) tool. The coordinate reference system and resolution used for spatial interpolation can be specified by the user. In addition, methods for computing zonal statistics using a fishnet grid, generation of a point shapefile, and calculation of point residuals between the bias ratio results at station locations and the overlapping interpolated surfaces are provided in the ``spatial`` submodule. 

The ``plot`` submodule in ``gridwxcomp`` provides tools for generating interactive diagnostic plots of both daily and monthly data, pairing the station data alongside its gridded counterpart. 

# Research Enabled by gridwxcomp

The most significant application of ``gridwxcomp`` was the development of bias correction surfaces that are applied to gridded reference evapotransipiration (ETo) data which are key inputs to some of the remote sensing ET models that comprise the OpenET platform [@Volk2024;@Melton2021]. Daily data from approximately 800 weather stations located in irrigated agricultural sites were curated, and the American Society of Civil Engineers (ASCE) standardized Penman-Monteith reference ET equation [@allen2005] was used to estimate ETo at the stations. Then ``gridwxcomp`` was used to pair these data with the nearest ETo data from the gridMET [@Abatzoglou2013] dataset over temporally consistent periods. The long-term average monthly ratios for station ETo relative to the gridded ETo were calculated for each point and saved as georeferenced data by ``gridwxcomp`` and were subsequently spatially interpolated using a kriging approach. The interpolated monthly surfaces are used within the OpenET platform to correct gridMET ETo data before it is used by most of the remote sensing ET models as a major scaling flux. 

# Co-author Roles

``gridwxcomp`` was developed through the following efforts:

* John M. Volk: Conceptualization, Software, Validation, Writing & Editting
* Christian Dunkerly: Conceptualization, Software, Validation
* Christopher Pearson: Conceptualization, Software, Validation
* Charles Morton: Conceptualization, Software
* Justin L. Huntington: Conceptualization, Funding acquisition, Resources 


# Acknowledgments

We would like to thank the Bureau of Reclamation, NASA Applied Sciences Program, and the Western States Water Use Program at the Desert Research Institute for providing funding for the development of this software. We also thank Sayantan Majumdar for his edits and suggestions.

# References

