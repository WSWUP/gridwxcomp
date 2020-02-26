"""
A package for comparing climate station time series data to gridded datasets, with tools for downloading `gridMET <http://www.climatologylab.org/gridmet.html>`_ data. Major functionality includes tools to: pair station locations with overlapping gridcells; download gridMET data from OpenDap or Google EarthEngine servers; calculate mean bias ratios between station and grid data and other statistics; spatial interpolation of bias ratios with multiple interpolation options; build geo-referenced files of point data and interpolated rasters; extract zonal means to gridcells; and graphics, e.g. time series comparisons and interpolated residuals at stations. Bias ratios calculated by ``gridwxcomp`` can be used to correct bias of grid to station data based on the properties of the stations. For example, monthly humidity ratios between station and grid for stations within agricultural settings can be used to estimate grid bias relative to agricultural locations. gridwxcomp includes an intuitive command line interface and a Python API.
"""

__name__ = 'gridwxcomp'
__author__ = 'John Volk and Chris Pearson'
__version__ = '0.1.3.post2'


from gridwxcomp.prep_input import prep_input
from gridwxcomp.download_gridmet_opendap import download_gridmet_opendap
from gridwxcomp.calc_bias_ratios import calc_bias_ratios
from gridwxcomp.interpgdal import InterpGdal
from gridwxcomp.spatial import make_points_file, make_grid, interpolate
from gridwxcomp.plot import daily_comparison, monthly_comparison
