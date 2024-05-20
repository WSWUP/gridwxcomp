# -*- coding: utf-8 -*-

import os
import pkg_resources
import pytest
from datetime import datetime
from pathlib import Path
from shutil import move, copy, rmtree

import numpy as np
import pandas as pd

from gridwxcomp.util import get_gridmet_meta_csv, parse_yr_filter
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
    d['ratio_dir'] = d.get('prep_metadata_outpath').parent/'test_ratios'

    def teardown():
        rmtree(Path('tests')/'example_data')
        if d['prep_metadata_outpath'].is_file():
            d['prep_metadata_outpath'].unlink()
        if d['prep_metadata_outpath_copy'].is_file():
            d['prep_metadata_outpath_copy'].unlink()
        if d['conus404_download_dir'].is_dir():
            rmtree(d['conus404_download_dir'])
        if d['ratio_dir'].is_dir():
            rmtree(d['ratio_dir'])
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

                
    def test_download_grid_data_conus404(self, data):
    	download_grid_data(
            data['prep_metadata_outpath'],
            dataset='conus404',
            export_bucket='openet',
            export_path=f'bias_correction_gridwxcomp_testing/gridwxcomp_conus404/',
            local_folder='tests',
            force_download=False)
    
        


#class TestDownloadGridmetOpenDap(object):
#
#    def setup_method(self):
#        self.gridmet_cols = ('date', 'year', 'month', 'day', 'centroid_lat',
#                    'centroid_lon', 'elev_m', 'u2_ms', 'tmin_c', 'tmax_c',
#                    'srad_wm2', 'ea_kpa', 'prcp_mm', 'etr_mm', 'eto_mm')
#
#            year_filter='2016',
#            update_data=False
#        )
#        assert download_dir.is_dir()
#        # check all four gridMET cells time series were downloaded
#        for id in self.gridmet_ids:
#            gridmet_file = download_dir.joinpath(
#                'gridmet_historical_{}.csv'.format(id)
#            )
#            assert gridmet_file.is_file()
#        df = pd.read_csv(download_dir/'gridmet_historical_441130.csv')
#        assert not set(df.columns).isdisjoint(self.gridmet_cols)
#        assert df.notna().all().all()
#        assert df.etr_mm.dtype == 'float'
#        #### test is input file from prep_metadata is updated with gridmet ts path
#        df = pd.read_csv(data['prep_metadata_outpath'])
#        assert 'GRID_FILE_PATH' in df.columns
#        # check correct path to gridmet time series file
#        assert not set(
#            df.loc[0,'GRID_FILE_PATH'].split(os.sep)).isdisjoint(
#                ['tests','od_gridmet_data','gridmet_historical_509011.csv'])
#
#    @pytest.mark.slow
#    def test_download_gridmet_od_multiple_years(self, data):
#        download_dir = data['prep_metadata_outpath'].parent/'od_gridmet_data'
#        download_gridmet_opendap(
#            data['prep_metadata_outpath'],
#            out_folder=download_dir,
#            year_filter='2016-2018',
#            update_data=False
#        )
#        for id in self.gridmet_ids:
#            gridmet_file = download_dir.joinpath(
#                'gridmet_historical_{}.csv'.format(id)
#            )
#            assert gridmet_file.is_file()
#        df = pd.read_csv(download_dir/'gridmet_historical_441130.csv',
#            parse_dates=True, index_col='date')
#        # make sure new years data exists and full range
#        assert df.index.year.min() == 2016
#        assert df.index.year.max() == 2018
#        assert not set(df.columns).isdisjoint(self.gridmet_cols)
#        assert df.notna().all().all()
#        assert df.etr_mm.dtype == 'float'
#        #### test is input file from prep_metadata is updated with gridmet ts path
#        df = pd.read_csv(data['prep_metadata_outpath'])
#        assert 'GRID_FILE_PATH' in df.columns
#        # check correct path to gridmet time series file
#        assert not set(
#            df.loc[0,'GRID_FILE_PATH'].split(os.sep)).isdisjoint(
#                ['tests','od_gridmet_data','gridmet_historical_509011.csv'])
#
#    @pytest.mark.slow
#    def test_download_gridmet_od_update_years(self, data):
#        """
#        Open one gridmet time series file after downloading and change
#        etr values to -1 for April 2016, redownload with update option and 
#        check that they are not -1 (i.e. they should have correct values)
#        """
#        download_dir = data['prep_metadata_outpath'].parent/'od_gridmet_data'
#        # first download only if this test is run first
#        if not download_dir.joinpath('gridmet_historical_441130.csv').is_file():
#            download_gridmet_opendap(
#                data['prep_metadata_outpath'],
#                out_folder=download_dir,
#                year_filter='2016',
#                update_data=False
#            )
#        df = pd.read_csv(download_dir/'gridmet_historical_441130.csv',
#            parse_dates=True, index_col='date')
#        df.loc['04/2016', 'etr_mm'] = -1
#        assert (df.loc['04/2016', 'etr_mm'] == -1).all()
#        # second download should update 2016
#        download_gridmet_opendap(
#            data['prep_metadata_outpath'],
#            out_folder=download_dir,
#            year_filter='2016',
#            update_data=True
#        )
#        # should be updated and no values of -1
#        df = pd.read_csv(download_dir/'gridmet_historical_441130.csv',
#            parse_dates=True, index_col='date')
#        assert not (df.loc['04/2016', 'etr_mm'] == -1).any()
#        # make sure other tests also pass
#        for id in self.gridmet_ids:
#            gridmet_file = download_dir.joinpath(
#                'gridmet_historical_{}.csv'.format(id)
#            )
#            assert gridmet_file.is_file()
#
#        assert not set(df.columns).isdisjoint(self.gridmet_cols)
#        assert df.notna().all().all()
#        assert df.etr_mm.dtype == 'float'
#        #### test is input file from prep_metadata is updated with gridmet ts path
#        df = pd.read_csv(data['prep_metadata_outpath'])
#        assert 'GRID_FILE_PATH' in df.columns
#        # check correct path to gridmet time series file
#        assert not set(
#            df.loc[0,'GRID_FILE_PATH'].split(os.sep)
#        ).isdisjoint(['tests','gridmet_data','gridmet_historical_509011.csv'])
#
#
#class TestCalcBiasRatios(object):
#
#    def setup_method(self):
#        # header of station metadata file should carry through plus ratio stats 
#        self.comp_header = (
#           'Jan_mean', 'Feb_mean', 'Mar_mean', 'Apr_mean', 'May_mean', 
#           'Jun_mean', 'Jul_mean', 'Aug_mean', 'Sep_mean', 'Oct_mean', 
#           'Nov_mean', 'Dec_mean', 'Jan_count', 'Feb_count', 'Mar_count', 
#           'Apr_count', 'May_count', 'Jun_count', 'Jul_count', 'Aug_count', 
#           'Sep_count', 'Oct_count', 'Nov_count', 'Dec_count', 'Jan_stdev', 
#           'Feb_stdev', 'Mar_stdev', 'Apr_stdev', 'May_stdev', 'Jun_stdev', 
#           'Jul_stdev', 'Aug_stdev', 'Sep_stdev', 'Oct_stdev', 'Nov_stdev', 
#           'Dec_stdev', 'Jan_cv', 'Feb_cv', 'Mar_cv', 'Apr_cv', 'May_cv', 
#           'Jun_cv', 'Jul_cv', 'Aug_cv', 'Sep_cv', 'Oct_cv', 'Nov_cv', 'Dec_cv',
#           'growseason_mean', 'summer_mean', 'annual_mean', 'growseason_count', 
#           'summer_count', 'annual_count', 'growseason_stdev', 'summer_stdev', 
#           'annual_stdev', 'growseason_cv', 'summer_cv', 'annual_cv', 
#           'start_year', 'end_year', 'GRIDMET_ID', 'FID', 'OBJECTID', 'Id', 
#           'State', 'Source', 'Status', 'STATION_LAT', 'STATION_LON', 'Date', 
#           'Station_ID', 'STATION_ELEV_FT', 'Comments', 'Location', 
#           'STATION_FILE_PATH', 'Irrigation', 'Website', 'STATION_ELEV_M', 
#           'LAT', 'LON', 'ELEV_M', 'ELEV_FT', 'FIPS_C', 'STPO', 'COUNTYNAME', 
#           'CNTYCATEGO', 'STATENAME', 'HUC8', 'GRID_FILE_PATH'
#        )
#
#    def test_calc_bias_ratios_default_options(self, data):
#        comp_out_file = data.get('ratio_dir')/'etr_mm_summary_comp_all_yrs.csv'
#        short_out_file = data.get('ratio_dir')/'etr_mm_summary_all_yrs.csv'
#
#        calc_bias_ratios(
#            input_path=data.get('prep_metadata_outpath'), 
#            out_dir=data.get('ratio_dir')
#        )
#
#        assert comp_out_file.is_file()
#        assert short_out_file.is_file()
#        df = pd.read_csv(
#            comp_out_file, 
#            index_col='STATION_ID', 
#            na_values=[-999]
#        )
#        assert not set(df.columns).isdisjoint(self.comp_header)
#
#    #TODO:
#    # make a long version to check year range mark it
#    # along with multiple year download, mark everything else as short 
#    # add tests for spatial and plot modules
#
#
#class TestSpatial(object):
#
#    def setup_method(self):
#        pass
#
#    def test_spatial_basic_options(self, data):
#        comp_out_file = data.get('ratio_dir')/'etr_mm_summary_comp_all_yrs.csv'
#
#        spatial.main(
#            input_file_path=comp_out_file, 
#            buffer=1,
#            scale_factor=1
#        )
#
#class TestPlot(object):
#
#    # note station bar plot gets tested in spatial default settings
#    def setup_method(self):
#        self.daily_plot_dir =  Path('tests')/'daily_comp_plots'
#        self.monthly_plot_dir =  Path('tests')/'monthly_comp_plots'
#
#    def test_daily_comparison(self, data):
#        input_path=data.get('prep_metadata_outpath') 
#        plot.daily_comparison(input_path, out_dir=self.daily_plot_dir)
#
#    def test_monthly_comparison(self, data):
#        input_path=data.get('prep_metadata_outpath') 
#        plot.monthly_comparison(input_path, out_dir=self.monthly_plot_dir)
# -*- coding: utf-8 -*-

import os
import pkg_resources
import pytest
from datetime import datetime
from pathlib import Path
