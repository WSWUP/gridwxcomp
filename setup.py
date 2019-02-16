from gridwxcomp import __version__
from setuptools import setup

requires = [
    'bokeh >= 1.0.4',
    'click >= 7.0',
    'cryptography == 2.3.1',
    'fiona == 1.7.13',
    'earthengine-api >= 0.1.164',
    'gdal == 2.2.4',
    'google-api-python-client >= 1.7.7',
    'numpy >= 1.15.4',
    'oauth2client >= 4.1.2', 
    'pandas == 0.23.4',
    'rasterstats == 0.13.0',
    'refet >= 0.3.7',
    'scipy == 1.1.0',
    'shapely == 1.6.4',
    'xlrd == 1.2.0'
]

tests_require = []

classifiers = [
    'License :: OSI Approved :: Apache Software License',
    'Programming Language :: Python :: 3.7',
    'Environment :: Console',
    'Development Status :: 4 - Beta',
    'Topic :: Scientific/Engineering',
    'Intended Audience :: Science/Research'
]

setup(
    name='gridwxcomp',
    description='Tools for comparing climate station data to gridMET data.',
    long_description='''gridwxcomp is a Python package for comparing station to gridMET climatic variables such as reference ET or air temperature. Its primary function is to calculate ratios between station and gridMET data for bias correction of gridMET values, e.g. in agricultural regions. Multiple options for spatial interpolation of bias ratios are available as well as tools for viewing time series comparisons of station and gridMET data. The final output includes zonal statistics of bias ratios for gridMET cells at monthly and annual time periods. 
    ''',
    author='Chris Pearson and John Volk',
    author_email='jmvolk@unr.edu',
    license='Apache',
    version=__version__,
    url='https://github.com/DRI-WSWUP/gridwxcomp',
    platforms=['Windows','Linux','Mac OS X'],
    classifiers=classifiers,
    packages=['gridwxcomp'],
    install_requires=requires,
    tests_require=tests_require,
    package_data={'gridwxcomp': ['example_data/*'],
        'gridwxcomp': ['gridmet_cell_data.csv']},
    include_package_data=True,
    entry_points='''
        [console_scripts]
        gridwxcomp=gridwxcomp.scripts.gridwxcomp:gridwxcomp
    '''
)
