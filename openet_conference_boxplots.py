import ee
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
import seaborn as sb
from gridwxcomp.prep_metadata import prep_metadata
from gridwxcomp.ee_download import download_grid_data
from gridwxcomp.plot import monthly_comparison
from gridwxcomp.calc_bias_ratios import calc_bias_ratios


# Config settings, must provide a project
ee_project = ""

# Optional flags to save time on rerunning data
force_download_of_data = False
generate_monthly_timeseries_plots = False
generate_bias_ratios = False

# Path to station metadata file with lat/long coords
station_meta_path = Path(Path.cwd() / "standardized_data/metadata_for_publication.csv")

# Path for output box plots
boxplot_output_path = Path(Path.cwd() / "bias_boxplots/")

# Below dictionaries hold all the information for making bias ratios and box plots
gridded_dataset_dict = {"urma": "URMA", 'gridmet': 'gridMET', 'nldas2': 'NLDAS',
                        'era5land': 'ERA5-Land', 'era5ag': 'ERA5-Ag'}

column_dict = {'Jan': 'January', 'Feb': 'February', 'Mar': 'March', 'Apr': 'April', 'May': 'May', 'Jun': 'June',
               'Jul': 'July', 'Aug': 'August', 'Sep': 'September', 'Oct': 'October', 'Nov': 'November',
               'Dec': 'December', 'grow': 'Growing Season', 'summer': 'Summer', 'annual': 'Annual'}

var_dict = {'eto': '$Daily\\:ET_{o}$', 'etr': '$Daily\\:ET_{r}$', 'prcp': 'Precipitation',
            'tmax': 'Temperature Max', 'tmin': 'Temperature Min', 'rs': 'Shortwave Radiation',
            'wind': '$Windspeed_{2m}$', 'tdew': 'Dewpoint Temperature', 'ea': 'Vapor Pressure'}


# run code
if __name__ == "__main__":

    # initialize earth engine
    ee.Authenticate()
    ee.Initialize(project=ee_project, opt_url='https://earthengine-highvolume.googleapis.com')

    # name of the dataset comparison is being made with
    for gridded_dataset_name in gridded_dataset_dict:

        # local path for config file
        config_path = f'gridwxcomp_config_{gridded_dataset_name}.ini'


        # local path for prep_metadata output
        gridwxcomp_input = f'gridwxcomp_metadata_{gridded_dataset_name}.csv'

        prep_metadata(
            station_meta_path,
            config_path,
            gridded_dataset_name,
            out_path=gridwxcomp_input)
        download_grid_data(
            gridwxcomp_input,
            config_path,
            local_folder=None,
            force_download=force_download_of_data)

        if generate_monthly_timeseries_plots:
            monthly_comparison(
                gridwxcomp_input,
                config_path,
                dataset_name=gridded_dataset_name)

        for var in var_dict:  # Iterate over vars in list. Valid entries found in calc_bias_ratios.py VAR_LIST
            '''
            calc_bias_ratios
               The purpose of this module is to calculate the monthly bias
                    (either the long-term-mean [default] or mean-of-annual)
                    between the gridded dataset and the station observations.
        
               IMPORTANT: The bias factors are either
                    (station - gridded) [Temperature] or (station / gridded) [Everything else].
                    This is the reverse relationship from how we normally think of modeled / observed, so a temperature bias
                    of -2 means that the gridded dataset is hotter than observed, not colder. This reversal was done to make
                    applying bias a matter of addition/multiplication.
        
               This module requires the configuration .ini file (stored under the example_data folder)
                    to be set up so that it can correct interpret the columns within the station and gridded data files.
                    The acceptable parameters for 'comparison_var' are coded in ACCEPTABLE_VAR_LIST within calc_bias_ratios.py
            '''
            if not generate_bias_ratios:
                continue

            calc_bias_ratios(
                gridwxcomp_input,
                config_path,
                f'gridwxcomp_bias_ratios_{gridded_dataset_name}',
                method='long_term_mean',
                grid_id_name='GRID_ID',
                comparison_var=var,
                grid_id=None,
                day_limit=10,
                years='all',
                comp=True)

    # generate boxplots now that all bias ratios have been generated
    sb.set_palette('Set2')
    plt.rc('xtick', labelsize=12)
    plt.rc('ytick', labelsize=12)
    plt.rc('axes', labelsize=12)
    plt.rc('axes', labelsize=12)
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Arial']
    for var in var_dict:
        for column in column_dict:
            boxplot_df = pd.DataFrame()

            for gridded_dataset_name in gridded_dataset_dict:
                try:
                    grid_df = pd.read_csv(
                        f'gridwxcomp_bias_ratios_{gridded_dataset_name}/{var}_summary_comp_all_yrs.csv',
                        keep_default_na=True, na_values=['-999', '-9999', -999, 'inf'])
                except FileNotFoundError:
                    continue
                grid_sel = grid_df[['STATION_ID', f'{column}_mean']].copy(deep=True)

                if var in ['tmax', 'tmin', 'tdew']:
                    grid_sel[f'{column}_mean_reciprocal'] = grid_sel[f'{column}_mean'] * -1
                    col_name = 'Mean Bias Ratio (Gridded Dataset - Observed)'
                else:
                    grid_sel[f'{column}_mean_reciprocal'] = 1 / grid_sel[f'{column}_mean']
                    col_name = 'Mean Bias Ratio (Gridded Dataset / Observed)'

                grid_sel.rename(columns={f'{column}_mean_reciprocal': col_name}, inplace=True)
                grid_sel['Gridded Dataset'] = gridded_dataset_dict[gridded_dataset_name]
                boxplot_df = pd.concat([boxplot_df, grid_sel], axis=0, ignore_index=True)

            # make figure
            plt.clf()  # reset current figure
            fig, ax = plt.subplots(figsize=(9, 4.5))
            sea_box = sb.boxplot(data=boxplot_df, x='Gridded Dataset', y=col_name, hue='Gridded Dataset', ax=ax)

            plt.title(f'{var_dict[var]} Mean Bias Ratios for {column_dict[column]}', fontsize=14)
            plt.tight_layout()

            if column in  ['grow', 'summer', 'annual']:
                boxplot_folder = boxplot_output_path / "seasonal_statistics"
            else:
                boxplot_folder = boxplot_output_path / "individual_months"

            Path.mkdir(boxplot_folder, exist_ok=True, parents=True)
            sea_box.figure.savefig(f'{boxplot_folder}/{var}_dataset_bias_for_{column}.png')

            plt.close()
