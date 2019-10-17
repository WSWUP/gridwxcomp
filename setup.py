import io, re
from setuptools import setup

with io.open("README.rst", "rt", encoding="utf8") as f:
    readme = f.read()

with io.open("gridwxcomp/__init__.py", "rt", encoding="utf8") as f:
    version = re.search(r"__version__ = \'(.*?)\'", f.read()).group(1)

requires = [
    'bokeh>=1.0.4',
    'click>=7.0',
    'fiona>=1.7.13',
    'gdal',
    'netCDF4',
    'numpy>=1.15',
    'pandas>=0.24',
    'rasterstats>=0.13',
    'refet>=0.3.7',
    'scipy>=1.1.0',
    'shapely==1.6.4',
    'xlrd==1.2.0'
]

tests_require = ['pytest']

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
    version=version,
    description='Compare meterological station data to gridded data',
    long_description=readme,
    author='John Volk and Chris Pearson',
    author_email='jmvolk@unr.edu',
    license='Apache',
    url='https://github.com/WSWUP/gridwxcomp',
    download_url='https://github.com/WSWUP/gridwxcomp/archive/v0.1.0.tar.gz',
    platforms=['Windows','Linux','Mac OS X'],
    classifiers=classifiers,
    packages=['gridwxcomp', 'gridwxcomp.scripts'],
    install_requires=requires,
    tests_require=tests_require,
    package_data={'gridwxcomp': ['example_data/*'],
        'gridwxcomp': ['env/*.yml'],
        'gridwxcomp': ['gridmet_cell_data.csv']},
    include_package_data=True,
    entry_points='''
        [console_scripts]
        gridwxcomp=gridwxcomp.scripts.gridwxcomp:gridwxcomp
    '''
)
