# -*- coding: utf-8 -*-

import os
import pkg_resources
import pytest
from datetime import datetime
from pathlib import Path
from shutil import move, copy, rmtree

import numpy as np
import pandas as pd

from gridwxcomp.util import read_config, parse_yr_filter
from gridwxcomp.prep_metadata import prep_metadata, _read_station_list
#from gridwxcomp.download_gridmet_opendap import download_gridmet_opendap
from gridwxcomp.ee_download import download_grid_data
from gridwxcomp.calc_bias_ratios import calc_bias_ratios
from gridwxcomp import spatial
from gridwxcomp import plot

@pytest.fixture(scope="session")
def data(request):
    """Prepare input data used for tests"""
    d = {}
    d['station_meta_path'] = Path(
            pkg_resources.resource_filename(
                'gridwxcomp', 'example_data/Station_Data.txt'
                )
            )
    # make copy of all example data for use in tests in tests/example_data
    example_files = [f for f in d.get('station_meta_path').parent.glob('*')]
    if not (Path('tests')/'example_data').is_dir():
        (Path('tests')/'example_data').mkdir(parents=True, exist_ok=True)
    for f in example_files:
        copy(f, Path('tests')/'example_data')
    d['station_meta_path'] = Path('tests')/'example_data'/'Station_Data.txt'
    bad_station_meta_path = Path('tests')/'example_data'/'Bad_Stations.txt'
    copy(d.get('station_meta_path'), bad_station_meta_path)
    df = pd.read_csv(bad_station_meta_path)
    df.drop('Station', axis=1, inplace=True)
    df.to_csv(bad_station_meta_path, index=False)
    d['bad_station_meta_path'] = bad_station_meta_path

    # mandatory renamed columns for station metadata
    d['need_cols'] = [
        'STATION_LAT',
        'STATION_LON',
        'STATION_FILE_PATH',
        'STATION_ID',
        ]
    
    d['prep_metadata_outpath'] =\
        d.get('station_meta_path').parent.parent/'merged_input.csv'
    d['prep_metadata_outpath_copy'] =\
        d.get('station_meta_path').parent.parent/'merged_input_cp.csv'
    d['conus404_config_path'] =\
        d.get('station_meta_path').parent/'gridwxcomp_config_conus404.ini'
    d['conus404_download_dir'] =\
        d.get('prep_metadata_outpath').parent/'conus404'
    d['conus404_output_dir'] =\
        d.get('prep_metadata_outpath').parent/'bias_outputs_conus404'

    def teardown():
        rmtree(Path('tests')/'example_data')
        if d['prep_metadata_outpath'].is_file():
            d['prep_metadata_outpath'].unlink()
        if d['prep_metadata_outpath_copy'].is_file():
            d['prep_metadata_outpath_copy'].unlink()
        if d['conus404_download_dir'].is_dir():
            rmtree(d['conus404_download_dir'])
        if d['conus404_output_dir'].is_dir():
            rmtree(d['conus404_output_dir'])
        if (Path('tests')/'__pycache__').is_dir():
            rmtree(Path('tests')/'__pycache__')
        if (Path('tests')/'daily_comp_plots').is_dir():
            rmtree(Path('tests')/'daily_comp_plots')
        if (Path('tests')/'monthly_comp_plots').is_dir():
            rmtree(Path('tests')/'monthly_comp_plots')

    request.addfinalizer(teardown)

    return d

class TestUtil(object):
    def setup_method(self):
        # for util.parse_yr_filter
        df = pd.DataFrame(index=pd.date_range('2000', '2015'))
        self.yr_df, self.yr_str = parse_yr_filter(df, '1998-2002', 'station1')

    
    #### tests for func util.parse_yr_filter
    def test_parse_yr_filter_data_type(self):
        assert type(self.yr_df.index) is pd.core.indexes.datetimes.DatetimeIndex

    def test_parse_yr_filter_data(self):
        # correct start date
        start = datetime(year=2000, month=1, day=1)
        assert datetime.date(self.yr_df.index.min()) == datetime.date(start)

    def test_parse_yr_filter_string(self):
        assert self.yr_str == '1998_2002'    

    def test_parse_yr_filter_bad_year_range(self):
        df = pd.DataFrame(index=pd.date_range('2000', '2015'))
        with pytest.raises(Exception) as e_info:
            d,y = parse_yr_filter(df, 'not-a_year', 'station1')

    def test_parse_yr_filter_bad_year(self):
        df = pd.DataFrame(index=pd.date_range('2000', '2015'))
        with pytest.raises(Exception) as e_info:
            d,y = parse_yr_filter(df, 25, 'station1')



class TestPrepMetadata(object):

    def test_read_station_list_inpath(self, data):
        # get example station metadata from install dir
        assert data['station_meta_path'].is_file(),\
                'example station metadata file not found'

    def test_read_station_list_ret_data(self, data):
        station_meta_df = _read_station_list(data['station_meta_path'])
        assert isinstance(station_meta_df, pd.DataFrame)
        # has mandatory columns (renamed)
        assert not set(station_meta_df).isdisjoint(data['need_cols'])
        # no missing data
        assert station_meta_df.loc[:,data['need_cols']].notna().all().all()
        # specific data value
        test_lon = station_meta_df.loc[
                station_meta_df['STATION_ID']=='Loa'
                ]['STATION_LON'].values[0]
        assert np.isclose(test_lon, -111.635832870077)

    def test_read_station_list_missing_input_col(self, data):
        with pytest.raises(Exception) as e_info:
            _read_station_list(data['bad_station_meta_path'])

    def test_prep_metadata(self, data):
        prep_metadata(
            data['station_meta_path'], 
            data['conus404_config_path'],
            'conus404',
            out_path=data['prep_metadata_outpath']
        )
        assert data['prep_metadata_outpath'].is_file()

    def test_prep_metadata_missing_input_col(self, data):
        with pytest.raises(Exception) as e_info:
            prep_metadata(data['bad_station_meta_path'])

    def test_prep_metadata_outfile_data(self, data):
        df = pd.read_csv(data['prep_metadata_outpath'])
        assert isinstance(df, pd.DataFrame)
        assert not set(df.columns).isdisjoint(data['need_cols'])
        test_lon = df.loc[df['STATION_ID']=='Loa']['STATION_LON'].values[0]
        assert np.isclose(test_lon, -111.635832870077)

class TestEEDownload(object):

    def setup_method(self):
        self.gridmet_cols = ('date', 'year', 'month', 'day', 'centroid_lat',
                    'centroid_lon', 'elev_m', 'u2_ms', 'tmin_c', 'tmax_c',
                    'srad_wm2', 'ea_kpa', 'prcp_mm', 'etr_mm', 'eto_mm')
        self.conus404_cols = ('ACSWDNB', 'ETO_ASCE', 'ETR_ASCE', 'PREC_ACC_NC',
            'PSFC', 'T2_MAX', 'T2_MIN', 'TD2', 'WIND10', 'date', 'station_name')
        self.export_path =\
            'bias_correction_gridwxcomp_testing/gridwxcomp_conus404'
 
                
    def test_download_grid_data_conus404(self, data):
        download_grid_data(
            data['prep_metadata_outpath'],
            config_path=data['conus404_config_path'],
            export_bucket='openet',
            export_path=self.export_path,
            local_folder='tests',
            force_download=False)

        # download dir
        assert Path('tests/conus404').is_dir()

        # check if output metadata file has been updated
        meta_df = pd.read_csv(data['prep_metadata_outpath'])

        # check all four station time series were downloaded
        for grid_file_path in meta_df.GRID_FILE_PATH:
            assert Path(grid_file_path).is_file()

        # open a downloaded grid data file and make checks
        df = pd.read_csv(
            meta_df.loc[
                meta_df.Station_ID == 'loau', 'GRID_FILE_PATH'
            ].values[0])
            
        assert not set(df.columns).isdisjoint(self.conus404_cols)
        assert df.notna().all().all()
        assert df.T2_MAX.dtype == 'float'

    @pytest.mark.slow
    def test_download_grid_data_conus404_overwrite(self, data):
        download_grid_data(
            data['prep_metadata_outpath'],
            config_path=data['conus404_config_path'],
            export_bucket='openet',
            export_path=self.export_path,
            local_folder='tests',
            force_download=True)

        # download dir
        assert Path('tests/conus404').is_dir()

        # check if output metadata file has been updated
        meta_df = pd.read_csv(data['prep_metadata_outpath'])

        # check all four station time series were downloaded
        for grid_file_path in meta_df.GRID_FILE_PATH:
            assert Path(grid_file_path).is_file()

        # open a downloaded grid data file and make checks
        df = pd.read_csv(
            meta_df.loc[
                meta_df.Station_ID == 'loau', 'GRID_FILE_PATH'
            ].values[0])
            
        assert not set(df.columns).isdisjoint(self.conus404_cols)
        assert df.notna().all().all()
        assert df.T2_MAX.dtype == 'float'
    
        


class TestCalcBiasRatios(object):

    def setup_method(self):
        # header of station metadata file should carry through plus ratio stats 
        self.comp_header = (
           'Jan_mean', 'Feb_mean', 'Mar_mean', 'Apr_mean', 'May_mean', 
           'Jun_mean', 'Jul_mean', 'Aug_mean', 'Sep_mean', 'Oct_mean', 
           'Nov_mean', 'Dec_mean', 'Jan_count', 'Feb_count', 'Mar_count', 
           'Apr_count', 'May_count', 'Jun_count', 'Jul_count', 'Aug_count', 
           'Sep_count', 'Oct_count', 'Nov_count', 'Dec_count', 'Jan_stdev', 
           'Feb_stdev', 'Mar_stdev', 'Apr_stdev', 'May_stdev', 'Jun_stdev', 
           'Jul_stdev', 'Aug_stdev', 'Sep_stdev', 'Oct_stdev', 'Nov_stdev', 
           'Dec_stdev', 'Jan_cv', 'Feb_cv', 'Mar_cv', 'Apr_cv', 'May_cv', 
           'Jun_cv', 'Jul_cv', 'Aug_cv', 'Sep_cv', 'Oct_cv', 'Nov_cv', 'Dec_cv',
           'growseason_mean', 'summer_mean', 'annual_mean', 'growseason_count', 
           'summer_count', 'annual_count', 'growseason_stdev', 'summer_stdev', 
           'annual_stdev', 'growseason_cv', 'summer_cv', 'annual_cv', 
           'start_year', 'end_year', 'GRIDMET_ID', 'FID', 'OBJECTID', 'Id', 
           'State', 'Source', 'Status', 'STATION_LAT', 'STATION_LON', 'Date', 
           'Station_ID', 'STATION_ELEV_FT', 'Comments', 'Location', 
           'STATION_FILE_PATH', 'Irrigation', 'Website', 'STATION_ELEV_M', 
           'LAT', 'LON', 'ELEV_M', 'ELEV_FT', 'FIPS_C', 'STPO', 'COUNTYNAME', 
           'CNTYCATEGO', 'STATENAME', 'HUC8', 'GRID_FILE_PATH'
        )

    def test_calc_bias_ratios_default_options(self, data):
        comp_out_file = data.get(
            'conus404_output_dir')/'etr_summary_comp_all_yrs.csv'
        short_out_file = data.get(
            'conus404_output_dir')/'etr_summary_all_yrs.csv'

        calc_bias_ratios(
            input_path=data.get('prep_metadata_outpath'), 
            config_path=data.get('conus404_config_path'),
            out_dir=data.get('conus404_output_dir')
        )

        assert comp_out_file.is_file()
        assert short_out_file.is_file()
        df = pd.read_csv(
            comp_out_file, 
            index_col='STATION_ID', 
            na_values=[-999]
        )
        assert not set(df.columns).isdisjoint(self.comp_header)

        assert df.loc[df.Station_ID == 'loau', 'ratio_method'].values[0]\
            == 'long_term_mean'


    def test_calc_bias_ratios_temp_delta_mean_of_annual_meth(self, data):
        comp_out_file = data.get(
            'conus404_output_dir')/'tmax_summary_comp_all_yrs.csv'
        short_out_file = data.get(
            'conus404_output_dir')/'tmax_summary_all_yrs.csv'

        calc_bias_ratios(
            input_path=data.get('prep_metadata_outpath'), 
            config_path=data.get('conus404_config_path'),
            out_dir=data.get('conus404_output_dir'),
            comparison_var='tmax',
            method='mean_of_annual'
        )

        assert comp_out_file.is_file()
        assert short_out_file.is_file()
        df = pd.read_csv(
            comp_out_file, 
            index_col='STATION_ID', 
            na_values=[-999]
        )
        assert not set(df.columns).isdisjoint(self.comp_header)

        assert df.loc[df.Station_ID == 'loau', 'ratio_method'].values[0]\
            == 'mean_of_annual'

class TestPlot(object):

    # note station bar plot gets tested in spatial default settings
    def setup_method(self):
        self.daily_plot_dir = Path('tests')/'daily_comp_plots'
        self.monthly_plot_dir = Path('tests')/'monthly_comp_plots'
    
    @pytest.mark.slow
    def test_daily_comparison(self, data):
        input_path=data.get('prep_metadata_outpath') 
        plot.daily_comparison(
            input_path, 
            data.get('conus404_config_path'),
            out_dir='tests')
        assert len(list(self.daily_plot_dir.glob('*'))) == 4
        assert Path(self.daily_plot_dir/'Bedrock').is_dir()
        assert len(
            list(Path(self.daily_plot_dir/'Bedrock').glob('*.html'))) == 12

    def test_monthly_comparison(self, data):
        input_path=data.get('prep_metadata_outpath') 
        plot.monthly_comparison(
            input_path, 
            data.get('conus404_config_path'),
            out_dir='tests')
        assert len(list(self.monthly_plot_dir.glob('*.html'))) == 4

    def test_station_bar_plot(self, data):
        input_path = data.get(
            'conus404_output_dir')/'tmax_summary_all_yrs.csv'
        plot.station_bar_plot(
            input_path, 
            bar_plot_layer='growseason_mean')
        outpath = data.get('conus404_output_dir')\
                /'station_bar_plots'/'growseason_mean.html'
        assert outpath.is_file()

        # same but using other file and coef. of variation
        input_path = data.get(
            'conus404_output_dir')/'etr_summary_comp_all_yrs.csv'
        plot.station_bar_plot(
            input_path, 
            bar_plot_layer='Jun_cv')
        outpath = data.get('conus404_output_dir')\
                /'station_bar_plots'/'Jun_cv.html'
        assert outpath.is_file()

class TestSpatial(object):

    def test_make_points_file(self, data):
        input_path = data.get(
            'conus404_output_dir')/'etr_summary_comp_all_yrs.csv'

        spatial.make_points_file(input_path)
        
        point_file = data.get(
            'conus404_output_dir')/'spatial'/'etr_summary_pts_wgs84.shp'
        assert point_file.is_file()

    def test_make_grid(self, data):
        config = read_config(data['conus404_config_path'])

        input_path = data.get(
            'conus404_output_dir')/'etr_summary_comp_all_yrs.csv'

        spatial.make_grid(
            input_path, 
            grid_res=config['output_data_resolution'],
            bounds=config['input_bounds'],
        )

        grid_path = data.get('conus404_output_dir')/'spatial'/'grid.shp'
        assert grid_path.is_file()
