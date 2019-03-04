"""
gridwxcomp is a Python package for comparing station to gridMET climatic variables such as reference ET or air temperature. Its primary function is to calculate ratios between station and gridMET data for bias correction of gridMET values, e.g. in agricultural regions. Multiple options for spatial interpolation of bias ratios are available as well as tools for viewing time series comparisons of station and gridMET data. The final output includes zonal statistics of bias ratios for gridMET cells at monthly and annual time periods. 
"""

__name__ = 'gridwxcomp'
__author__ = 'John Volk and Chris Pearson'
__version__ = '0.0.33.dev3'


from gridwxcomp.prep_input import prep_input
from gridwxcomp.download_gridmet_ee import download_gridmet_ee
from gridwxcomp.calc_bias_ratios import calc_bias_ratios
from gridwxcomp.interpgdal import InterpGdal
from gridwxcomp.spatial import make_points_file, make_grid, interpolate
from gridwxcomp.daily_comparison import daily_comparison
