""" 
Desert Research Institute
Calculate Monthly ETr Bias Correction shapefile from raw station and GRIDMET datasets
"""
	
import pandas as pd
import numpy as np
import os

#%%
# Set Working Directory
os.chdir('D:\USBR_UC')

#Specify data paths
corrected_data_path = 'D:\USBR_UC\ETr_Bias_Correction\Corrected_Output_Data'
gridmet_data_path = 'D:\USBR_UC\ETr_Bias_Correction\GRIDMET_Data'

#Import Station/GRIDMET meta data shapefile           
station_attribute = pd.read_csv('ETr_Bias_Correction\station_gridmetID_joinFINAL.csv', sep=',')

##Limit row processing range
#start = 5
#end = 6

# Speicfy desired output fields

out = pd.DataFrame(columns =['GRIDMET_ID', 'GRID_Lat','GRID_Lon',
                             'Jan',	'Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'])


out_comp = pd.DataFrame(columns =['GRIDMET_ID', 'GRID_Lat','GRID_Lon','GRID_Elev_ft','State','FileName', 'Station_ID', 'Station_Lat','Station_Lon', 'Station_Elev_ft',
                                  'Source', 'Website',  'Date First Ob.',  'Irrigation', 'Comments', 'Elev_Diff_ft',
                             'Jan',	'Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec', 'Ann_avg', 'JunthruAug_avg', 'GS_avg',
                             'Jan_count',	'Feb_count','Mar_count','Apr_count','May_count','Jun_count','Jul_count','Aug_count',
                             'Sep_count','Oct_count','Nov_count','Dec_count'])


#Station filename skip_list based on Corrected_Output_Data Station in station_attribute .shp
skip_list = ['Cortez','Center 1','Afton','Pleasant Valley','Ferron','Olathe 2',
             'Buena Vista','Montrose','Orchard Mesa','Bedrock','Encampment','Birch Springs Draw',
             'Safford','Elk Mountain']
#skip_list = ['AftonDaily_Raw', 'BigPineyWY_Daily', 'BirchSpringsUT_Daily', 'LymanWY_Daily',
#             'ElmoUT_Daily', 'FerronUT_Daily', 'SiltCO_Daily', 'ElkMountainDaily_Raw',
#             'YellowJacketCO_Daily', 'TowaocCO_Daily', 'FarmingtonNM_Daily', 
#             'Center1Daily_Raw', 'CortezCO_Daily'] 

#Loop through each station/gridmet pair
for index, row in station_attribute.iterrows(): 
    #Limit iteration during development#
#    if index < start:
#        continue
#    if index >= end:
#       break
    print('Shapefile Row Number: {}').format(index)
    #Skip if station in skip_list
    if row.Station in skip_list:
        print('Station in Skip List. Skipping {}.').format(row.Station)
        continue
    else:
        #Import Corrected Station Data
        station_path = os.path.join(corrected_data_path, '{}_output.xlsx').format(row.FileName) 

        #Skip If FILE DOES NOT EXIST
        if not os.path.exists(station_path):
            print('SKIPPING {}. NO STATION FILE FOUND.').format(station_path)
            continue
        else:
            station_data = pd.read_excel(station_path, sheet_name = 'Corrected Data')
            
            #Import GRIDMET Data
            grid_path = os.path.join(gridmet_data_path, 'gridmet_historical_{}.csv').format(row.GRIDMET_ID)
            #Skip if GRIDMET FILE DOES NOT EXIST
            if not os.path.exists(grid_path):
                print('SKIPPING {}. NO FILE GRIDMET FOUND.').format(grid_path)
                continue
            else:
                grid_data = pd.read_csv(grid_path, sep=',',parse_dates=True, index_col='Date')
                result = pd.concat([station_data['Calc_ETr (mm)'], grid_data['ETr_ASCE']], axis=1, join_axes=[station_data.index])
                #Remove results with 
                result = result.dropna()
                
                # Monthly ETr Sums including Count
                monthly = result.groupby([lambda x: x.year, lambda x: x.month]).agg(['sum','count'])
                  
                # Remove Totals with Less Than XX Days
                day_limit = 10
                monthly = monthly[monthly['ETr_ASCE','count']>=day_limit]
                
                #Calculated Station/GRIDMET ratio from monthly sums
                ratio = pd.DataFrame(columns = ['Ratio'])
                ratio['Ratio'] =(monthly['Calc_ETr (mm)','sum'])/(monthly['ETr_ASCE','sum'])
                
                #Rebuild Index DateTime
                ratio['year'] = ratio.index.get_level_values(0).values 
                ratio['month'] = ratio.index.get_level_values(1).values
                ratio.index = pd.to_datetime(ratio.year*10000+ratio.month*100+15,format='%Y%m%d')       
                
                #Create final monthly ratios df; Fill with monthly CHECK STAT BELOW
                final_ratios = pd.DataFrame(np.nan, index=range(1,13,1), columns=['Ratio'])                              
                final_ratios['Ratio'] = ratio.groupby(ratio.index.month).median()
                
                final_counts = pd.DataFrame(np.nan, index=range(1,13,1), columns=['Count'])
                final_counts['Count'] = ratio.groupby(ratio.index.month).count()
                
                ann_ratio = final_ratios.mean()
                
                #Growing Season Apr (inclusive) through Oct (exlusive); Zero-based
                gs_ratio = final_ratios.iloc[3:10].mean()
                #Summer Season Jun (inclusive) through Aug (exlusive); Zero-based
                juntoaug_ratio = final_ratios.iloc[5:8].mean()
               
                #Write to output dataframe
                out = out.append ({'GRIDMET_ID': row.GRIDMET_ID,
                                   'GRID_Lat': row.LAT,
                                   'GRID_Lon': row.LON,                               
    #                               'FileName': row.FileName,
    #                               'Station_Lat': row.LATDECDEG,
    #                               'Station_Lon': row.LONGDECDEG,
    #                               'Elev_m': row.Elev_m, 
                                   'Jan': final_ratios.Ratio[1],
                                   'Feb': final_ratios.Ratio[2],
                                   'Mar': final_ratios.Ratio[3],
                                   'Apr': final_ratios.Ratio[4],
                                   'May': final_ratios.Ratio[5],
                                   'Jun': final_ratios.Ratio[6],
                                   'Jul': final_ratios.Ratio[7],
                                   'Aug': final_ratios.Ratio[8],
                                   'Sep': final_ratios.Ratio[9],
                                   'Oct': final_ratios.Ratio[10],
                                   'Nov': final_ratios.Ratio[11],
                                   'Dec': final_ratios.Ratio[12]}, ignore_index=True)
    
    
    
                #Write to output dataframe
                out_comp = out_comp.append ({'GRIDMET_ID': row.GRIDMET_ID,
                                   'GRID_Lat': row.LAT,
                                   'GRID_Lon': row.LON,
                                   'GRID_Elev_ft': row.ELEV_FT_GR,
                                   'State': row.State,
                                   'FileName': row.FileName,
                                   'Station_ID': row.Station_ID,
                                   'Station_Lat': row.LATDECDEG,
                                   'Station_Lon': row.LONGDECDEG,
                                   'Station_Elev_ft': row.Elev_FT,
                                   'Source': row.Source,
                                   'Website': row.Website,
                                   'Date First Ob.': row.Date,
                                   'Irrigation': row.Irrigation,
                                   'Comments': row.Comments,                                   
                                   'Elev_Diff_ft': row.Elev_DIFF,
                                   'Jan': final_ratios.Ratio[1],
                                   'Feb': final_ratios.Ratio[2],
                                   'Mar': final_ratios.Ratio[3],
                                   'Apr': final_ratios.Ratio[4],
                                   'May': final_ratios.Ratio[5],
                                   'Jun': final_ratios.Ratio[6],
                                   'Jul': final_ratios.Ratio[7],
                                   'Aug': final_ratios.Ratio[8],
                                   'Sep': final_ratios.Ratio[9],
                                   'Oct': final_ratios.Ratio[10],
                                   'Nov': final_ratios.Ratio[11],
                                   'Dec': final_ratios.Ratio[12],
                                   'Ann_avg': ann_ratio.Ratio,
                                   'JunthruAug_avg':  juntoaug_ratio.Ratio,
                                   'GS_avg':  gs_ratio.Ratio,
                                   'Jan_count': final_counts.Count[1],
                                   'Feb_count': final_counts.Count[2],
                                   'Mar_count': final_counts.Count[3],
                                   'Apr_count': final_counts.Count[4],
                                   'May_count': final_counts.Count[5],
                                   'Jun_count': final_counts.Count[6],
                                   'Jul_count': final_counts.Count[7],
                                   'Aug_count': final_counts.Count[8],
                                   'Sep_count': final_counts.Count[9],
                                   'Oct_count': final_counts.Count[10],
                                   'Nov_count': final_counts.Count[11],
                                   'Dec_count': final_counts.Count[12],}, ignore_index=True)
    
###%% Export Full Dataset .txt file
out_path = "D:\USBR_UC\GIS"
out_file = os.path.join(out_path, 'BiasCorrectionData_Median_Final.txt')
out_comp.to_csv(out_file, sep ='\t', index = False)    
#    
       
#%% Create Point Shapefile with Monthly Ratio Fields
#USES GRIDMET GRIDCELL CENTER FOR LAT/LON
# Import arcpy module
import arcpy
#Set file overwrite 
overwrite_flag = True
arcpy.env.overwriteOutput = overwrite_flag

#Set arcpy workspace
arcpy.env.workspace = 'D:\USBR_UC\GIS'

#Set file paths
out_path = "D:\USBR_UC\GIS"
out_name = "BiasCorrection_point.shp"
geometry_type = "POINT"
file_path = os.path.join(out_path, out_name)

#Delete files if they exist (overwrite?)
if arcpy.Exists(file_path):
    arcpy.Delete_management(file_path)

#Set Geographic Coordinate system "WGS 1984" (factory code=4326)
sr = arcpy.SpatialReference(4326)

#Execute CreateFeatureclass
arcpy.CreateFeatureclass_management(out_path, out_name, geometry_type, spatial_reference = sr)


#Add Fields to .shp
fields =  ['GRIDMET_ID', 'Jan', 'Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
for field in fields:
    arcpy.AddField_management(file_path, field, "DOUBLE")

#Insert cursor fieldnames
cur_fields = ['GRIDMET_ID', 'SHAPE@Y', 'SHAPE@X', 'Jan', 'Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

#cur_fields = ['GRIDMET_ID', 'SHAPE@Y', 'SHAPE@X', 'Jan', 'Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

#Set-up cursor
cursor = arcpy.da.InsertCursor(file_path, cur_fields)

#Insert new rows (ORDER MUST MATCH CURSOR FIELD LIST ABOVE; pandas df to np.array)
for row in np.array(out):
    cursor.insertRow(row)
# Delete cursor object
del cursor

#%% Interpolate Bias Correction Maps for Each Month (IDW 4km cell size?)
#Checkout spatial analyst license
arcpy.CheckOutExtension('Spatial')
#Month List
month_list= ['Jan', 'Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
#month_list= ['Jan']

# Set Processing Extent
#Original Extent .shp file
extent_file = os.path.join(out_path,'UpperCO_HUC2_Boundary.shp')

#Extent + buffer file path
extent_buffer_file = os.path.join(out_path,'UpperCO_HUC2_Boundary_Buffer.shp')

#Add Buffer (4km ~ '0.041666667 DecimalDegrees'based on ee GRIDMET Grid; Apply 10x Buffer)
arcpy.Buffer_analysis(extent_file, extent_buffer_file, '3 DecimalDegrees')

#Set Processing Extent to extent_buffer file
arcpy.env.extent = extent_buffer_file

#Set Snap Raster to UC ET Cells File (Optional)
arcpy.env.snapRaster = 'D:\et-demands\upper_co_full\gis\ETCells_ras.img'

#Loop through month list creating IDW correction grids
for month in month_list:
    #IDW Inputs
    inPointFeatures = file_path
    zField = month
    #(0.041666667 decimal degrees taken from GEE GRIDMET tiff) HARDCODED FOR NOW
    cellSize = 0.041666667
    power = 2
    
    # Execute IDW
    outIDW = arcpy.sa.Idw(inPointFeatures, zField, cellSize, power)
    
    # Save the output (Delete/Overwrite OLD FILES)
    out_file_path = os.path.join(out_path,'{}_Median_IDW_Final.img').format(month)
    if arcpy.Exists(out_file_path):
        arcpy.Delete_management(out_file_path)
    outIDW.save(out_file_path)
    del outIDW

#%% Create Growing Season Avg Raster   
growing_season= ['Apr','May','Jun','Jul','Aug','Sep','Oct']
ras_list =[]
for month in growing_season:
    out_file_path = os.path.join(out_path,'{}_Median_IDW_Final.img').format(month)
    ras_list.append(out_file_path)

calc = arcpy.sa.CellStatistics(ras_list, statistics_type = "MEAN")

calc.save(os.path.join(out_path,'GrowingSeasonMean_AprOct_Median_IDW_Final.img'))    
del calc  

#%% Create Annual Avg Raster   
ras_list =[]
for month in month_list:
    out_file_path = os.path.join(out_path,'{}_Median_IDW_Final.img').format(month)
    ras_list.append(out_file_path)

calc = arcpy.sa.CellStatistics(ras_list, statistics_type = "MEAN")

calc.save(os.path.join(out_path,'Annual_Median_IDW_Final.img'))    
del calc  

#%% Create Growing Season Avg Raster   
summer_season= ['Jun','Jul','Aug']
ras_list =[]
for month in summer_season:
    out_file_path = os.path.join(out_path,'{}_Median_IDW_Final.img').format(month)
    ras_list.append(out_file_path)

calc = arcpy.sa.CellStatistics(ras_list, statistics_type = "MEAN")

calc.save(os.path.join(out_path,'Summer_JunAug_Median_IDW_Final.img'))    
del calc  



#%% Test Code
##%% Example Scatter Plot using Seaborn package
#import seaborn
#plot1=(seaborn.jointplot(out.GRID_Lat, out.Elev_m, kind="reg",annot_kws=dict(stat="r"))
#    .set_axis_labels("Elevation (m)", "Ratio (Station/GRIDMET)"));                                
#plot1
#    
#
##%%Multiple Linear Regression (Testing)
#import statsmodels.formula.api as smf
#formula = 'Jun~GRID_Lat+Elev_m'
##formula = 'Insitu_Chl_a~green_sur'
#
#est=smf.ols(formula=formula, data=out).fit()
#est.summary()

