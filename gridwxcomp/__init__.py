"""
A package for comparing climate station time series data to `gridMET <http://www.climatologylab.org/gridmet.html>`_ data. Functionality includes tools to: pair station locations with overlapping gridMET cells; download gridMET data using the Google Earth Engine API; calculate mean bias ratios between station and gridMET data and other statistics; build geo-referenced files of point data; conduct spatial interpolation of correction ratios with multiple interpolation options; extract zonal mean data to a subset of gridMET cells; and produce graphical plots of time series comparisons. Bias ratios calculated by ``gridwxcomp`` can be used to correct bias of gridMET to station data based on the properties of the stations. For example, monthly humidity ratios between station and gridMET for stations within agricultural settings can be used to estimate gridMET bias relative to agricultural locations. gridwxcomp includes an intuitive command line interface and a Python API.
"""

__name__ = 'gridwxcomp'
__author__ = 'John Volk and Chris Pearson'
__version__ = '0.0.54'


from gridwxcomp.prep_input import prep_input
from gridwxcomp.download_gridmet_ee import download_gridmet_ee
from gridwxcomp.calc_bias_ratios import calc_bias_ratios
from gridwxcomp.interpgdal import InterpGdal
from gridwxcomp.spatial import make_points_file, make_grid, interpolate
from gridwxcomp.daily_comparison import daily_comparison
from gridwxcomp.monthly_comparison import monthly_comparison
