# -*- coding: utf-8 -*-
"""
Contains objects for interfacing ``gridwxcomp`` with gdal command line tools
used for spatial interpolation of scattered point data.
"""
import numpy as np
import os
import pandas as pd
from pathlib import Path
import subprocess
import xml.dom.minidom
import xml.etree.cElementTree

from .spatial import get_subgrid_bounds, zonal_stats, calc_pt_error, make_points_file
from .plot import station_bar_plot


class InterpGdal(object):
    """
    Usage of gdal command line tools within ``gridwxcomp``, currently
    utilizes the `gdal_grid <https://www.gdal.org/gdal_grid.html>` command line tool.
    
    Arguments:
        summary_csv_path (str): path to [var]_summary_comp CSV file created 
            by :mod:`gridwxcomp.calc_bias_ratios` containing point bias ratios,
            lat, and long. 
            
    Attributes:
        summary_csv_path (:obj:`pathlib.Path`): absolute path object to 
            input ``summary_csv_path``.
        layers (list): list of layers in ``summary_csv_path`` to interpolate
            e.g. when using :meth:`InterpGdal.gdal_grid` with ``layer="all"``
            defaults to :attr:`InterpGdal.default_layers`.
        grid_bounds (tuple or None): default None. Extent for interpolation
            raster in order (min long, max long, min lat, max lat).
        interp_meth (str): default 'invdist'. gdal_grid interpolation method.
        interped_rasters (list): empty list that is appended with 
            :obj:`pathlib.Path` objects to interpolated rasters after using
            :meth:`InterpGdal.gdal_grid`.
        params (dict or None): default None. After :meth:`InterpGdal.gdal_grid`
            ``params`` is updated with the last used interpolation parameters
            in the form of a dictionary with parameter names as keys.

    Example:
        The :class:`InterpGdal` class is useful for experimenting with multiple 
        interpolation algorithms provided by gdal which are optimized and
        often sensitive to multiple parameters. We can use the object
        to loop over a range of parameter combinations to test how algorithms
        perform, we might pick a single layer to test, in this case the
        growing season bias ratios between station and gridded reference
        evapotranspiration (etr_mm). Below is a routine to conduct a
        basic sensitivity analysis of the power and smooth parameters of the 
        inverse distance to a power interpolation method
        
            >>> from gridwxcomp.interpgdal import InterpGdal
            >>> import os
            >>> # create instance variable from summary csv
            >>> summary_file = 'PATH/TO/etr_mm_summary_comp.csv'
            >>> # define grid resolution in meters
            >>> grid_resolution = 1000  # will be modified by argument scale_factor if it's provided
            >>> # create a InterpGdal instance
            >>> test = InterpGdal(summary_file)
            >>> layer = 'growseason_mean' # could be a list
            >>> # run inverse distance interpolation with different params
            >>> for p in range(1,10):
            >>>     for s in [0, 0.5, 1, 1.5, 2]:
            >>>         # build output directory based on parameter sets
            >>>         out_dir = os.path.join('spatial', 'invdist',
            >>>             'power={}_smooth={}'.format(p, s))
            >>>         params = {'power': p, 'smooth': s}
            >>>         test.gdal_grid(grid_resolution, out_dir=out_dir, layer=layer, params=params)

        Note, we did not assign the interpolation method 'invdist' because it 
        is the default. To use another interpolation method we would
        assign the ``interp_meth`` kwarg to :meth:`InterpGdal.gdal_grid`.
        Similarly, we could experiment with other parameters which all can be
        found in the class attribute :attr:`InterpGdal.default_params`. The
        instance variable :attr:`InterpGdal.params` can also be used to save
        metadata on parameters used for each run.
        
    """
    # gdal_grid interpolation methods
    interp_methods = ('average',
                      'invdist',
                      'invdistnn',
                      'linear',
                      'nearest')
    """
        interp_methods (tuple): gdal_grid interpolation algorithms.
    """

    default_layers = ('Jan',
                      'Feb',
                      'Mar',
                      'Apr',
                      'May',
                      'Jun',
                      'Jul',
                      'Aug',
                      'Sep',
                      'Oct',
                      'Nov',
                      'Dec',
                      'growseason',
                      'annual',
                      'summer')
    """
    default_layers (tuple): Layers to interpolate created by 
        :mod:`gridwxcomp.calc_bias_ratios` and then :mod: `gridwxcomp.make_points_file, 
        e.g. "Jan" in the shapefile which is "Jan_mean" found in ``summary_csv_path``.
    """

    # interp params, method key, param dic as values
    default_params = {
        'invdist': {
            'power': 2,
            'smoothing': 0,
            'radius1': 0,
            'radius2': 0,
            'angle': 0,
            'max_points': 0,
            'min_points': 0,
            'nodata': -999
        },
        'invdistnn': {
            'power': 2,
            'smoothing': 0,
            'radius': 10,
            'max_points': 12,
            'min_points': 0,
            'nodata': -999
        },
        'average': {
            'radius1': 0,
            'radius2': 0,
            'angle': 0,
            'min_points': 0,
            'nodata': -999
        },
        'linear': {
            'radius': -1,
            'nodata': -999
        },
        'nearest': {
            'radius1': 0,
            'radius2': 0,
            'angle': 0,
            'nodata': -999
        }
    }
    """
    default_params (dict): Dictionary with default parameters for each
        interpolation algorithm, slightly modified from GDAL defaults.
        Keys are interpolation method names, keys are dictionaries with
        parameter names keys and corresponding values.
    """

    var_residual_names = {
        'annual_mean': 'annual_res',
        'growseason_mean': 'grow_res',
        'summer_mean': 'summer_res'
    }
    """
    var_residual_names (dict): Dictionary that maps names of bias ratios 
        to the name of their intepolated residual if they are too long for
        storing as a field name in an ESRI shapefile (i.e. > 10 chars).
    """

    def __init__(self, summary_csv_path):
        
        if not Path(summary_csv_path).is_file():
            raise FileNotFoundError(
                'Summary CSV file: {} not found!'.format(
                    Path(summary_csv_path).absolute())
            )
            
        self.summary_csv_path = Path(summary_csv_path).absolute()
        self.layers = list(self.default_layers)  # mutable instance attr.
        self.grid_bounds = None
        self.interp_meth = 'invdist'
        self.interped_rasters = []  # appended with Path objects
        self.params = None  # to hold last used interp. parameters as dict

    def gdal_grid(self, proj_dict, proj_name, layer='all', out_dir='', interp_meth='invdist',
                  params=None, scale_factor=1, z_stats=True, res_plot=True,
                  grid_id_name='GRID_ID', options=None):
        """
        Run gdal_grid command line tool to interpolate point ratios.
        
        For further information on theinterpolation algorithms including 
        their function, parameters, and options see 
        `gdal_grid <https://www.gdal.org/gdal_grid.html>`
        
        Arguments:
            proj_dict (dict): Dictionary containing CRS information, grid resolution, and bounds.
            proj_name (str): entry in proj_dict to pull projection info from
            layer (str, list): default 'all'. Name of summary file column
                to interpolate, e.g. 'Jan_mean', or list of names. If 'all'
                use all variables in mutable instance attribute "layers".
            out_dir (str): default ''. Output directory to save rasters and
                zonal stats CSV, always appended to the root dir of the
                input summary CSV parent path that contains point ratios.
            interp_meth (str): default 'invdist'. gdal interpolation algorithm
            params (dict, str, or None): default None. Parameters for 
                interpolation algorithm. See examples for format rules.
            bounds (tuple or None): default None. Tuple of bounding coordinates
                in the following order (min long, max long, min lat, max lat) 
                which need to be in decimal degrees or meters. If None, get extent from
                locations of climate stations in input summary CSV with 25 cell
                buffer.
            scale_factor (float, int): default 1. Scaling factor to apply to 
                original grid_res fishnet to create resampling resolution.
            z_stats (bool): default True. Calculate zonal means of 
                interpolated surface to gridded cells in fishnet and save to a
                CSV file. The CSV file will be saved to the same directory as 
                the interpolated raster file(s).
            res_plot (bool): default True. Make bar plot for residual (error)
                between interpolated and station value for ``layer``.
            options (str or None): default None. Extra command line options for
                gdal_grid spatial interpolation.
                
        Returns:
            None
                
        Examples:
            The default interpolation algorithm 'invdist' or inverse distance 
            weighting to a power to interpolate bias ratios in a summary CSV 
            file that was first created by :mod:`gridwxcomp.calc_bias_ratios`.
            The default option will interpolate all layers in :obj:`InterpGdal.default_layers` 
            and calculate zonal statistics for all layers. It will also assume 
            the boundaries of the rasters are defined by the centroid locations 
            of the *outer* station locations in the input summary CSV plus a 25 
            gridded cell buffer, the pixel size will be 400m (scale_factor=0.1).
            We also limit the interpolation to two layers, growing season and
            annual mean bias ratios,
            
            >>> from gridwxcomp.interpgdal import InterpGdal
            >>> summary_file = 'PATH/TO/[var]_summary_comp.csv'
            >>> out_dir = 'default_params'
            >>> layers = ['growseason_mean', 'annual_mean']
            >>> proj_dict = {
            >>>     'wgs84': {'bounds': {'xmin': -115.0, 'xmax': -101.0, 'ymin': 35.5, 'ymax': 42.5},
            >>>               'resolution': 0.1, 'crs_id': 'EPSG:4326'}}
            >>> # create a InterpGdal instance
            >>> test = InterpGdal(summary_file)
            >>> # run inverse distance interpolation
            >>> test.gdal_grid(proj_dict, 'wgs84', out_dir=out_dir, layer=layers)
            
            Note, zonal statistics to gridded cells and interpolated residual
            plots are computed by default. A gridded fishnet must have been
            previously created using :func:`gridwxcomp.spatial.make_grid`.
            
            After running the code above the following files will be created 
            in the 'default_params' directory which will be built in the same
            location as the input summary CSV::

                default_params/           
                ├── annual_mean.tiff
                ├── annual_mean.vrt
                ├── etr_mm_summary_comp_all_yrs.csv
                ├── etr_mm_summary_pts_wgs84.cpg
                ├── etr_mm_summary_pts_wgs84.dbf
                ├── etr_mm_summary_pts_wgs84.prj
                ├── etr_mm_summary_pts_wgs84.shp
                ├── etr_mm_summary_pts_wgs84.shx
                ├── etr_mm_summary_pts_lcc.cpg
                ├── etr_mm_summary_pts_lcc.dbf
                ├── etr_mm_summary_pts_lcc.prj
                ├── etr_mm_summary_pts_lcc.shp
                ├── etr_mm_summary_pts_lcc.shx
                ├── zonal_stats.csv
                ├── growseason_mean.tiff
                ├── growseason_mean.vrt
                └── residual_plots/
                    ├── annual_res.html
                    └── grow_res.html

            GeoTiff interpolated raster files are now created for select layers
            as well as VRT (virtual raster) meta files that store info
            on each raster's data source. The file "gridded_stats.csv" contains
            grid_id as an index and each layer zonal mean as columns. For
            example,
            
                ========== ================== ================== 
                grid_id      growseason_mean    annual_mean
                ========== ================== ================== 
                511747     0.9650671287940088 0.9078723876243023 
                510361     0.9658465063428492 0.9097255715561022 
                508975     0.9667075970344162 0.9117676407214926 
                ========== ================== ================== 
                
            There are several :class:`InterpGdal` instance attributes that 
            may be useful, for example to see the parameters that were used 
            for the last call to :meth:`InterpGdal.gdal_grid`
            
            >>> test.params
            {'power': 2,
             'smoothing': 0,
             'radius1': 0,
             'radius2': 0,
             'angle': 0,
             'max_points': 0,
             'min_points': 0,
             'nodata': -999}
                 
            Or to find the paths to the interpolated raster files that
            have been created by the instance (all), the "interped_rasters"
            instance attribute is a list of all :obj:`pathlib.Path` objects
            of absolute paths of raster files. To get them as strings,
            
            >>> list(map(str, test.interped_rasters))
            ['PATH/TO/growseason_mean.tiff',
             'PATH/TO/annual_mean.tiff']
             
            Similary, the raster extent that was used and will be used again for
            any subsequent calls of :meth:`InterpGdal.gdal_grid` can be 
            retrieved by
             
            >>> test.grid_bounds 
            (-111.74583329966664,
             -108.74583330033335,
             38.21250000003333,
             40.462499999966674)
            
        Raises:
            KeyError: if interp_meth is not a valid gdal_grid interpolation
                algorithm name. 
        """
        
        cwd = os.getcwd()

        # out_dir as subdir from dir with summary CSV 
        out_dir = (Path(self.summary_csv_path).parent/Path(out_dir)).resolve()

        if not out_dir.is_dir():
            print('{}\nDoes not exist, creating directory'.format(str(out_dir)))
            out_dir.mkdir(parents=True, exist_ok=True)
    
        source_file = Path(self.summary_csv_path).name
        source = source_file.replace('.csv', '_tmp')
        
        if interp_meth not in InterpGdal.interp_methods:
            err_msg = '{} not a valid interpolation method'.format(interp_meth)
            raise KeyError(err_msg)
        self.interp_meth = interp_meth
            
        # look up default parameters for interpolation method
        if not params:
            params = InterpGdal.default_params.get(self.interp_meth)
        # avoid zero NA values for zonal_stats, set default NA rep to -999
        if not params.get('nodata'):
            params['nodata'] = -999
        # update instance parameters for later reference
        self.params = params
        # parameters to command line input str :name=value:name=value ...
        param_str = ':'+':'.join('{}={}'.format(key, val) for (key, val) in params.items())

        # get boundary and projection info
        self.grid_bounds = proj_dict[proj_name]['bounds']
        xmin = self.grid_bounds['xmin']
        xmax = self.grid_bounds['xmax']
        ymin = self.grid_bounds['ymin']
        ymax = self.grid_bounds['ymax']
        grid_res = proj_dict[proj_name]['resolution']

        # create the points shapefile and construct path to the correct one in LCC projection
        make_points_file(self.summary_csv_path, grid_id_name='GRID_ID')
        path_root = os.path.split(self.summary_csv_path)[0]
        file_name = os.path.split(self.summary_csv_path)[1]
        var_name = file_name.split('_summ')[0]
        spatial_dir = os.path.join(path_root, 'spatial')
        lcc_shapefile_path = os.path.join(spatial_dir, '{v}_summary_pts_{p}.shp'.format(v=var_name, p=proj_name))

        # to parse options, like --config GDAL_NUM_THREADS update here
        if not options:
            options = ''
                
        def _run_gdal_grid(layer):
            """reuse if running multiple layers"""
            existing_layers = pd.read_csv(self.summary_csv_path).columns
            if f'{layer}_mean' not in existing_layers:
                print('\nError: {} does not exist in input CSV:\n {}'.format(
                    f'{layer}_mean', self.summary_csv_path), '\nSkipping interpolation.')
                return
            
            # move to out_dir to run gdal command
            os.chdir(out_dir)

            lcc_tiff_file = '{}_lcc.tiff'.format(layer)
            resampled_tiff_file = '{}.tiff'.format(layer)

            # TODO out_file is currently the LCC version but is being passed to zonal stats which will expect WGS84
            #   fix this and update the output status message
            # add raster path to instance if not already there (overwritten)
            out_file = out_dir.joinpath(lcc_tiff_file)
            # print message to console/logging about interpolation
            grid_var = Path(self.summary_csv_path).name.split('_summ')[0]
            _interp_msg(grid_var, layer, self.interp_meth, (scale_factor * grid_res), out_file)
            # build command line arguments
            cmd = (f'gdal_grid -l {var_name}_summary_pts_lcc -zfield {layer} '
                   f'-a {interp_meth}{param_str} '
                   f'-txe {xmin} {xmax} -tye {ymax} {ymin} '
                   f'-tr {grid_res * scale_factor} {grid_res * scale_factor} '
                   f'-ot Float32 -of GTiff {lcc_shapefile_path} '
                   f'{lcc_tiff_file}')

            # run gdal_grid with arguments, x-platform
            p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE)
            out, err = p.communicate()
            if err:
                print(err)
            else:
                if out_file not in self.interped_rasters:
                    self.interped_rasters.append(out_file)

            p.stdout.close()
            p.stderr.close()

            grid_extent = proj_dict['wgs84']['bounds']
            resample_res = 0.1  # in decimal degrees
            # Resample using gdal warp, to the gridded extent, and exact cell size, using bilinear method
            # instead of "bilinear" you could use "near" for nearest neighbor, documentation says the fastest and worst
            cmd = (f'gdalwarp -overwrite -r bilinear '
                   f'-te {grid_extent["xmin"]} {grid_extent["ymin"]} {grid_extent["xmax"]} {grid_extent["ymax"]} '
                   f'-t_srs EPSG:4326 -tr {resample_res} {resample_res} {lcc_tiff_file} {resampled_tiff_file}')
            p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = p.communicate()
            if err:
                print(err)
            p.stdout.close()
            p.stderr.close()

            # Change directory back to CWD
            os.chdir(cwd)

            # calculate interpolated values and error at stations
            # only calc residuals for mean bias ratios, i.e. not std dev, etc.
            if layer in InterpGdal.default_layers:
                calc_pt_error(self.summary_csv_path, out_dir, layer, grid_var, grid_id_name=grid_id_name)
            else:
                # delete tmp summary csv used in interpgdal _make_pt_vrt method 
                # normally deleted in calc_pt_error
                tmp_csv = str(self.summary_csv_path).replace('.csv', '_tmp.csv')
                if Path(tmp_csv).resolve().is_file():
                    Path(tmp_csv).resolve().unlink()

            # zonal means extracted to gridded dataset cell index
            if z_stats:
                zonal_stats(self.summary_csv_path, out_file, grid_id_name=grid_id_name)
            # residual (error) bar plot, only for mean bias ratios
            if res_plot and layer in InterpGdal.default_layers:
                layer = self.var_residual_names.get(
                    layer, 
                    layer.replace('mean', 'res')
                )
                y_label = 'residual (interpolated minus station value)'
                title = 'layer: {} algorithm: {} (gdal_grid) res: {} deg.'.\
                        format(layer, self.interp_meth, grid_res)
                res_plot_dir = Path(out_dir)/'residual_plots'
                subtitle = 'parameters: {}'.format(params)
                source_file = Path(out_dir)/Path(self.summary_csv_path).name

                station_bar_plot(source_file, layer, out_dir=res_plot_dir,
                                 y_label=y_label, title=title, subtitle=subtitle)

        # run interpolation and zonal statistics depending on layer kwarg
        if layer == 'all':  # potential for multiprocessing
            for item in self.layers:
                _run_gdal_grid(item)
        # run single field
        elif isinstance(layer, str):
            _run_gdal_grid(layer)
        # run select list or tuple of layers
        elif isinstance(layer, (list, tuple)):
            for item in layer:
                _run_gdal_grid(item)
                
                
def _prettify(elem):
    """
    Return an indented, multiline XML string for a XML element tree.
    """
    rough_string = xml.etree.cElementTree.tostring(elem, 'utf-8', method='xml')
    reparsed = xml.dom.minidom.parseString(rough_string)
    return reparsed.toprettyxml(indent='  ')


def _interp_msg(grid_var, layer, function, res, out_file):
    print(
        '\nInterpolating {g} point bias ratios for: {t}\n'.format(g=grid_var, t=layer),
        'Using the "{}" method\n'.format(function),
        'Resolution (pixel size) of output raster: {} meters'.format(res)
    )
    print('GeoTIFF raster will be saved to: \n', os.path.abspath(out_file))
