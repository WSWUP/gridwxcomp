"""
A Python package for calculating bias correction factors between climate stations and gridMET variables. Correction ratios can be used to adjust gridMET estimated reference evapotranspiration (ETr) or other climatic variables to weather station data e.g. stations within agriculture settings. The package includes tools to pair stations locations with gridMET cells, download gridMET data using Google Earth Engine API, calculate point ratios, build geo-referenced files of point data, conduct spatial interpolation of correction ratios with various interpolation options, and perform zonal extraction of mean correction ratios to a subset of gridMET cells around the stations. gridwxcomp includes an intuitive command line interface and a set of Python functions.
"""

__name__ = 'gridwxcomp'
__author__ = 'John Volk and Chris Pearson'
__version__ = '0.0.51'


from gridwxcomp.prep_input import prep_input
from gridwxcomp.download_gridmet_ee import download_gridmet_ee
from gridwxcomp.calc_bias_ratios import calc_bias_ratios
from gridwxcomp.interpgdal import InterpGdal
from gridwxcomp.spatial import make_points_file, make_grid, interpolate
from gridwxcomp.daily_comparison import daily_comparison
