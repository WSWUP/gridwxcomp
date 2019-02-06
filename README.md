# gridwxcomp
Station-based Bias Correction of Gridded Weather for Agricultural Applications

A Python module that calculates monthly bias correction factors that can be utilized to bias correct gridMET estimated ETr or other climatic variables to weather station observed data (e.g. stations within agriculture settings). The process includes spatial interpolation of correction ratios with various interpolation options and zonal extraction of mean correction ratios to gridMET cells. ``gridwxcomp`` includes a command line interface as well as a set of Python submodules.

A full documentation website is under development.

## Quick start

### Installation

Currently ``gridwxcomp`` needs requires manual download for installation but will be available on on the Python Package Index [PyPI](https://pypi.org/) for installation with pip soon. 

To download first navigate to the directory where you would like to use the software and then clone the repository using git from the command line:

```
git clone https://github.com/DRI-WSWUP/gridwxcomp.git
```

If you do not have git, [here are installation instructions](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git).

Next to install dependencies you will need [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html). 

Once conda is installed move to the ``env`` directory and create the ``gridwxcomp`` environment using the appropriate environment file included:

```
cd env
conda env create -f env_windows.yml
```

To activate the environment whenever you need to use ``gridwxcomp`` just run

```
activate gridwxcomp
```

or, on Linux

```
source activate gridwxcomp
```

Lastly, ``gridwxcomp`` uses the Google Earth Engine API to download gridMET data, therefore you will need a Google account and before the first use on a machine you will need to verify your account. From the command line type:

```
python -c "import ee; ee.Initialize()"
```

and follow the instructions.

### Basic use from command line

This workflow will use the example station data given in "gridwxcomp/gridwxcomp/ETrBias_DataPackage". It will calculate bias ratios between station and gridMET ETr, estimate spatial GeoTIFF rasters at 400m resolution, and calculate zonal statistics of mean bias ratios for each gridMET cell in the region of the stations. 

From ``gridwxcomp`` root directory run

```
cd gridwxcomp
python prep_input.py -i ETrBias_DataPackage/Station_Data.txt  
```

This will result in the file "merged_input.csv", next to download matching gridMET climate time series run

```
python download_gridmet_ee.py -i merged_input.csv -o test_gridmet_data -y 2016-2017
```

In this case the years 2016-2017 are used because the test data is short plus it saves time by downloading a single year. Next to calculate monthly bias ratios and save to CSV files run

```
python calc_bias_ratios.py -i merged_input.csv -o test_ratios -c
```

Last, to calculate interpolated spatial surfaces of bias ratios and extract zonal means:

```
python spatial.py -i test_ratios/summary_comp.csv -b 5
```

The ``[-b 5]`` option indicates to exapnd the rectangular bounding area for interpolation by five gridMET cells.

The final output files including monthly bias ratios for each gridMET cell, GeoTIFF rasters of interpolated ratios, a fishnet grid with gridMET id values, and a point shapefile of station ratios should all be created within the "test_ratios" directory.
