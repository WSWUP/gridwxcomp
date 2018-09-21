""" 
Extract monthly ETr correction ratios from IDW correction grids for each
cell in model run
"""


import pandas as pd
#import numpy as np
import os
import arcpy

##%%
## Set Working Directory
#os.chdir('D:\USBR_UC')
#out_path = "D:\USBR_UC\GIS"

#%% Interpolate Bias Correction Maps for Each Month (IDW 4km cell size?)
#Checkout spatial analyst license
arcpy.CheckOutExtension('Spatial')

# Set file overwrite
overwrite_flag = True
arcpy.env.overwriteOutput = overwrite_flag

# Paths should be taken from .ini eventually
gis_ws = 'D:\upper_co\gis'
et_cells_path = 'D:\upper_co\gis\ETCells.shp'
cells_dd_path = os.path.join(gis_ws, 'ETCells_dd.shp')
cells_ras_path = os.path.join(gis_ws, 'ETCells_ras.img')

#Create dd et_cells
arcpy.Project_management(et_cells_path, cells_dd_path,
                         arcpy.SpatialReference('WGS 1984'))

# Set arcpy environmental parameters
arcpy.env.extent = cells_dd_path
arcpy.env.outputCoordinateSystem = cells_dd_path
station_id_field = 'GRIDMET_ID'  
 
# Convert cells_dd to cells_ras (0.041666667 from GEE GRIDMET tiff) HARDCODED
arcpy.FeatureToRaster_conversion(cells_dd_path, station_id_field,
                                 cells_ras_path, 0.041666667)

#Month List
month_list= ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep',
             'Oct', 'Nov', 'Dec']

#create list of bias correction rasters
ras_list =[]
for month in month_list:
    out_file_path = os.path.join('D:\USBR_UC\GIS',
                                 '{}_Median_IDW_Final.img').format(month)
    ras_list.append(out_file_path)


pt_path = os.path.join('D:\upper_co\gis', 'ETCells_pt.shp')
arcpy.RasterToPoint_conversion(cells_ras_path, pt_path, 'VALUE')   

# Extract all idw raster values to point .shp
arcpy.sa.ExtractMultiValuesToPoints(pt_path, ras_list, 'NONE')

#Write Values to ETrRatiosMon.txt
out_header = ['Met Node ID', 'Met Node Name', 'Jan', 'Feb', 'Mar', 'Apr', 'May',
              'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']

fields =['GRID_CODE', 'Jan_MEDIAN', 'Feb_MEDIAN', 'Mar_MEDIAN', 'Apr_MEDIAN',
         'May_MEDIAN', 'Jun_MEDIAN', 'Jul_MEDIAN', 'Aug_MEDIAN', 'Sep_MEDIAN',
         'Oct_MEDIAN', 'Nov_MEDIAN', 'Dec_MEDIAN']

pt_array = arcpy.da.FeatureClassToNumPyArray(pt_path, fields)



#Write pandas dataframe to .csv for et-demands input
out_df = pd.DataFrame(pt_array, columns = fields)
#Add empty 'Met Name' Column (OPTIONAL IN ET DEMANDS)
out_df['Met Node Name'] =""

out_order =['GRID_CODE', 'Met Node Name', 'Jan_MEDIAN', 'Feb_MEDIAN',
            'Mar_MEDIAN', 'Apr_MEDIAN', 'May_MEDIAN', 'Jun_MEDIAN',
            'Jul_MEDIAN', 'Aug_MEDIAN', 'Sep_MEDIAN', 'Oct_MEDIAN',
            'Nov_MEDIAN', 'Dec_MEDIAN']

out_file = os.path.join('D:\upper_co\static', 'ETrRatiosMon.txt')
out_df.to_csv(out_file,  columns=out_order, sep='\t',
          header =out_header, index=False)

#
##np array export option
#fmt_str = ['%10.0f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f' 
#    '%10.5f %10.5f %10.5f %10.5f %10.5f']
#out_file_np = os.path.join(out_path, 'ETrRatiosMon.txt')
#np.savetxt(fname=out_file_np, X=pt_array, fmt=fmt_str, delimiter = '\t',
#            header = '\t'.join(out_header), comments='')

#%% Testing
# refet_ratios_df = pd.read_table('D:\upper_co\static\ETrRatiosMon.txt',
#                                 delimiter='\t', header=None,
#                                 skiprows=1 - 1, na_values=['NaN'])