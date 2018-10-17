""" 
C.Pearson 3/28/2018
Desert Research Institute
Compare Monthly gridMET and Weather Station Values (windspeed, solar rad, Tmin, Tmax, (Tdew, RH, vp)?) Always RH for Now
"""
	
import pandas as pd
import numpy as np
import os

from bokeh.plotting import figure, output_file, show, save
from bokeh.layouts import gridplot
#from bokeh.models import Legend

#%%
# Set Working Directory
os.chdir('D:\USBR_UC')

#Specify data paths
corrected_data_path = 'D:\USBR_UC\ETr_Bias_Correction\Corrected_Output_Data'
gridmet_data_path = 'D:\USBR_UC\ETr_Bias_Correction\GRIDMET_Data'

#Import Station/GRIDMET meta data shapefile           
station_attribute = pd.read_csv('ETr_Bias_Correction\station_gridmetID_joinFINAL.csv', sep=',')

# List of variables to compare (STATION/gridMET ORDER SHOULD MATCH)
station_vars = ['Windspeed (m/s)','Rs (w/m2)','RHAvg (%)','TMax (C)','TMin (C)','Calc_ETr (mm)']

gridmet_vars = ['Ws_2m','Srad','RH','Tmax_C','Tmin_C','ETr_ASCE']

##Limit row processing range
#start = 0
#end = 1
#Loop through each station/gridmet pair
for index, row in station_attribute.iterrows(): 
#    Limit iteration during development
#    if index < start:
#        continue
#    if index >= end:
#       break
    print('Shapefile Row Number: {}').format(index)
    
    grid_elev = row.ELEV_M_GRI
    
    station_elev = row.Elev_FT

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
            
            #Pressure from elevation (ASCE Eqn. 34) kPa
            grid_data['P_kPa'] = 101.3*((293-0.0065*grid_elev)/293)**5.26
                     
            #Specific Humidity and elevation to Vapor Pressure (kPa)
            grid_data['ea'] = grid_data.q*grid_data.P_kPa/(0.622+0.378*grid_data.q)
            
            #Saturated Vapor Pressure
            grid_data['Tavg_C'] = (grid_data.Tmin_C + grid_data.Tmax_C )/2
            grid_data['e_sat'] = 0.6108*np.exp((17.27*grid_data.Tavg_C)/(grid_data.Tavg_C+237.3))
            #RH (%)
            grid_data['RH'] = (grid_data.ea/grid_data.e_sat)*100

            #Combine station and gridMET dataframes (only plotting variables)                                      
            result = pd.concat([station_data[station_vars], grid_data[gridmet_vars]], axis=1, join_axes=[station_data.index])
            #Remove results with na
            result = result.dropna()
            
            # Monthly averages including count
            monthly = result.groupby([lambda x: x.year, lambda x: x.month]).agg(['mean','count'])
              
            # Remove Totals with Less Than XX Days
            day_limit = 10
            monthly = monthly[monthly['Windspeed (m/s)','count']>=day_limit]
            
            #Rebuild Index DateTime
            monthly['year'] = monthly.index.get_level_values(0).values 
            monthly['month'] = monthly.index.get_level_values(1).values
            monthly.index = pd.to_datetime(monthly.year*10000+monthly.month*100+15,format='%Y%m%d')
          
            #Output to HTML file
            out_file_path = os.path.join('D:\USBR_UC','Comparison_Plots','{}_MonthComp.html').format(row.Station.replace(" ", "").replace("(","").replace(")",""))
            output_file(out_file_path)
                       
            #list of x variables
            x_var_list= ['Calc_ETr (mm)','Rs (w/m2)','Windspeed (m/s)','RHAvg (%)','TMax (C)','TMin (C)']
            #list of y variables
            y_var_list= ['ETr_ASCE','Srad','Ws_2m','RH','Tmax_C','Tmin_C']
            #title list
            title_list= ['ETr: Monthly Average','Solar Radiation: Monthly Average','Windspeed: Monthly Average','Relative Humidity: Monthly Average',
                    'Tmax: Monthly Average', 'Tmin: Monthly Average']
            #timeseries y label list
            ts_ylabel_list = ['ETr (mm)','Solar Radiation (w/m2)','Windspeed (m/s)','Relative Humidity (%)','Temperature (C)','Temperature (C)']
            #scatter xlabel list
            xlabel_list= ['Station ETr (mm)','Station SR (w/m2)','Station Windspeed (m/s)','Station RH (%)','Station Temperature (C)','Station Temperature (C)']
            #scatter ylabel list
            ylabel_list=['gridMET ETr (mm)','gridMET SR (w/m2)','gridMET Windspeed (m/s)','gridMET RH (%)','gridMET Temperature (C)','gridMET Temperature (C)']
            #legendx list
            legendx_list = ['Station','Station','Station','Station','Station','Station']
            #legend y list
            legendy_list = ['gridMET','gridMET','gridMET','gridMET','gridMET','gridMET']
            #empty list to append figures to
            figure_list = [] 
                       
            #loop through and create figures for each variable using variables and plot labels from lists above
            for i, (x_var, y_var, title, ts_ylabel, xlabel, ylabel, legendx, legendy) in enumerate(zip(x_var_list, y_var_list, title_list, ts_ylabel_list,
                                                                                        xlabel_list, ylabel_list, legendx_list, legendy_list)):
                if i == 0:
                    #Initial Timeseries Plot to establish xrange for link axes
                    p1 = figure(plot_width=1000, plot_height=400, x_axis_type="datetime",title = title, y_axis_label = ts_ylabel)
                    p1.line(monthly.index.to_pydatetime(), monthly[x_var,'mean'],  color="navy", alpha=0.5, legend=legendx,line_width=2)
                    p1.line(monthly.index.to_pydatetime(), monthly[y_var,'mean'],  color="red", alpha=0.5, legend=legendy,line_width=2)
                else:
                    #Timeseries Plots after first pass using x_range from first pass
                    p1 = figure(plot_width=1000, plot_height=400, x_axis_type="datetime",title = title, y_axis_label = ts_ylabel, x_range=p1.x_range)
                    p1.line(monthly.index.to_pydatetime(), monthly[x_var,'mean'],  color="navy", alpha=0.5, legend=legendx,line_width=2)
                    p1.line(monthly.index.to_pydatetime(), monthly[y_var,'mean'],  color="red", alpha=0.5, legend=legendy,line_width=2)                      
                    
                #1 to 1 Plot
                #Regression through Zero (https://stackoverflow.com/questions/9990789/how-to-force-zero-interception-in-linear-regression/9994484#9994484)                               
                m = np.linalg.lstsq(monthly[x_var,'mean'].values.reshape(-1,1), monthly[y_var,'mean'])[0][0]                              
                r_x, r_y = zip(*((i, i*m ) for i in range(int(np.min([monthly[y_var,'mean'],monthly[x_var,'mean']])-2),
                                 int(np.max([monthly[y_var,'mean'],monthly[x_var,'mean']])+3),1)))
                #Plots
                p2 = figure(plot_width=400, plot_height=400,x_axis_label = xlabel, y_axis_label = ylabel, title = 'Slope Through Zero: m = {}'.format(round(m,4)))
                p2.circle(monthly[x_var,'mean'], monthly[y_var,'mean'], size=20, color="navy", alpha=0.5)
                p2.line([int(np.min([monthly[y_var,'mean'],monthly[x_var,'mean']])-2),int(np.max([monthly[y_var,'mean'],monthly[x_var,'mean']])+2)],
                         [int(np.min([monthly[y_var,'mean'],monthly[x_var,'mean']])-2),int(np.max([monthly[y_var,'mean'],monthly[x_var,'mean']])+2)],
                          color = "black", legend = '1 to 1 line')
                p2.line(r_x, r_y, color="red", legend = 'Reg thru zero')
                p2.legend.location = "top_left"
                
                #Append [p1, p2] to figure_list (create list of lists)
                figure_list.append([p1, p2])
                
            #Plot all figures in list         
            fig = gridplot(figure_list, toolbar_location="left")
            #Show the figure
#            show(fig)
            #Save the figure
            save(fig)
                 
