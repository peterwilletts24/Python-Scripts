#Monthly

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.mlab as ml

import datetime

from dateutil.relativedelta import relativedelta

import re

import numpy as np

from math import sin, cos, atan2, radians, sqrt

import scipy.interpolate

import gc

import pdb

Cross_Section_Title = 'Vizag_to_Afghanistan'
station_list_cs=[43150, 42867, 43014, 42339, 40990, 40948]

first_station=43150

date_min=datetime.datetime(2011,5,1,0,0,0)
date_max=datetime.datetime(2011,10,1,0,0,0)

delta = relativedelta(weeks=+1)


    
def station_info_search(stat):
    for line in station_metadata:
     if "%s" % stat in line: 
      st = line[2].lower().title().replace('_',' ')
      lo = float(line[3])
      la = float(line[4])
    return st,la,lo

def calculate_distance_from_first_station(stat, first_station_lon, first_station_lat, station_lat, station_lon):

 fslat_rad = radians(first_station_lat)
 fslon_rad = radians(first_station_lon)
 lat_rad = radians(station_lat)
 lon_rad = radians(station_lon)

 #Haversine Formula
    
 a = sin((lat_rad-fslat_rad)/2)**2 + cos(lat_rad) * cos(fslat_rad) * sin((lon_rad-fslon_rad)/2)**2
 c = 2 * atan2(sqrt(a), sqrt(1-a))
 d = 6371 * c

 return d

def variable_name_index_match(variable, variable_list):
 for key, value in variable_list.iteritems():   # iter on both keys and values
        if key.startswith('%s' % variable):
                arr_index_var=value 
 return arr_index_var
            
def variable_cat(var_index, station_list_cs):
    var_cat=[]
    distances=[]
    date_min_max=[]
    
    for stat in station_list_cs:
     load_file = np.load('/nfs/a90/eepdw/Data/Observations/Radiosonde_Numpy/Radiosonde_Cross_Section_'
                        'IND_SOUNDING_INTERP_MEAN_%s_%s_%s_%s_%s.npz' 
                        % (Cross_Section_Title, date_min.strftime('%Y%m%d'), date_max.strftime('%Y%m%d'), delta, stat))
  
     print load_file['date_bin_mean_all_dates_one_station'].shape

     if date_min_max ==[]:
      date_min_max=np.empty(load_file['min_max_date_bin'].shape)
        
     station_title, station_lon, station_lat = station_info_search(stat)
    
     dist_from_first_station = calculate_distance_from_first_station(stat, first_station_lon, first_station_lat, station_lat, station_lon)
    
     print dist_from_first_station
     #print load_file['date_bin_mean_all_dates_one_station'][:,var_index,:].shape
     var_cat.append(load_file['date_bin_mean_all_dates_one_station'][:,var_index,:])
     distances.append(dist_from_first_station)
     #pdb.set_trace()
     #if load_file['min_max_date_bin'].any() != np.NAN:
     #date_min_max=np.ma.masked_outside(load_file['min_max_date_bin'], date_min, date_max ).data
     date_min_max = np.where((load_file['min_max_date_bin']>date_min) & (load_file['min_max_date_bin']<date_max), load_file['min_max_date_bin'], date_min_max )
    print np.array(var_cat).shape
    print date_min_max
    return np.array(var_cat), np.array(distances, dtype=float), date_min_max

def station_name_plot(station_list_cs, first_station, yi):    
 y_offset_text=0
 first_station_title, first_station_lon, first_station_lat = station_info_search(first_station)
    
 
 for stat in station_list_cs: 
     
   station_title, station_lon, station_lat = station_info_search(stat)
    
   dist_from_first_station = calculate_distance_from_first_station(stat, first_station_lon, first_station_lat, station_lat, station_lon)
    
   plt.axvline(x=dist_from_first_station, ymin=0, ymax=1, label=station_title, color='k')
   plt.text(dist_from_first_station+0.1,max(yi)/100+20,station_title,rotation=-45)

   y_offset_text=+1   
    
def u_v_winds(wind_direction, wind_speed):
 wind_rad = math.rad(wind_direction)
 u_wind=-(wind_speed/10)*sin(wind_rad)
 v_wind=-(wind_speed/10)*cos(wind_rad)
 return u_wind,v_wind

def grid_data_cs(pressure, distance, param):
 xi=np.linspace(0, max(distance), 200) 
 #yi=np.linspace(np.nanmin(pressure), np.nanmax(pressure), 500) 
 yi=np.linspace(5000, 100000, 50) # Points for pressure interpolation

 #yi=np.array([1000, 925, 850, 700, 500, 400, 300, 250, 200, 150, 100, 70, 50, 30, 20,10], dtype=float)
 #yi=np.array([10, 20, 30, 50, 70, 100, 150, 200, 250, 300, 400, 500, 700, 850, 925, 1000]*100, dtype=float)
 

 try:
    zi = ml.griddata(distance, pressure,param,xi, yi, interp='nn')
    #zi = scipy.interpolate.griddata((distance, pressure),  param, (xi[None,:],yi[:,None]), method='linear')
 except Exception, e:
  print e
 return xi,yi,zi 
 #return xi,yi
# def plot_rad_cs(xi,yi,zi, min_contour, max_contour):
#  clevs = np.linspace(min_contour, max_contour,256)
#  ticks = (np.arange(min_contour, max_contour+tick_interval,tick_interval))

#  plt.figure(figsize=(14,8))
#  cmap=plt.cm.jet
#  cont = plt.contourf(xi,yi/100, zi, clevs, cmap=cmap, extend='both')

 

#  cbar = plt.colorbar(cont, orientation='vertical', pad=0.05, extend='both', format = '$%d$')
#  #cbar.set_label('$W m^{-2}$') 
#  cbar.set_ticks(np.arange(min_contour, max_contour+tick_interval,tick_interval))
#  cbar.set_ticklabels(['${%d}$' % i for i in ticks])
    
#  plt.gca().invert_yaxis()
#  plt.ylabel('Pressure (hPa)')
#  plt.xlabel('km from first station')
    
#  return cont,cbar

def plot_rad_cs_winds(xi,yi,zi, min_contour, max_contour, wind_gridded):
 clevs = np.linspace(min_contour, max_contour,256)
 ticks = (np.arange(min_contour, max_contour+tick_interval,tick_interval))

 plt.figure(figsize=(14,8))
 cmap=plt.cm.jet
 cont = plt.contourf(xi,yi/100, zi, clevs, cmap=cmap, extend='both')
 plt.contour(xi,yi/100, zi, clevs, cmap=cmap, extend='both')
 

 cbar = plt.colorbar(cont, orientation='vertical', pad=0.05, extend='both', format = '$%d$')
 #cbar.set_label('$W m^{-2}$') 
 cbar.set_ticks(np.arange(min_contour, max_contour+tick_interval,tick_interval))
 cbar.set_ticklabels(['${%d}$' % i for i in ticks])
    
 plt.gca().invert_yaxis()
 plt.ylabel('Pressure (hPa)')
 plt.xlabel('km from first station')
    
 return cont,cbar

# def date_bin_plot(i, date_bin, concat_plot_variable, pressures, distances, min_contour, max_contour):
#   nan_mask = np.ma.masked_array(np.array(concat_plot_variable[:,i,:], dtype=float).flatten(), np.isnan(np.array(concat_plot_variable[:,i,:], dtype=float).flatten()))
#  #print nan_mask
#   print concat_plot_variable.shape
 
#   try:    
#    if nan_mask.mask.all() == False:
#     print nan_mask
#     xi,yi, zi = grid_data_cs(np.array(pressures[:,i,:], dtype=float).flatten(), np.repeat(distances, concat_plot_variable[:,i,:].shape[1]).flatten(), nan_mask) 
    
#     cont,cbar = plot_rad_cs(xi, yi, zi, min_contour, max_contour)
#     station_name_plot(station_list_cs, first_station, yi)
#   except Exception, e:
#    print e
#   return cont,cbar

def date_bin_plot_winds(i, date_bin, concat_plot_variable, pressures, distances, min_contour, max_contour, wind_to_plot):
  nan_mask = np.ma.masked_array(np.array(concat_plot_variable[:,i,:], dtype=float).flatten(), np.isnan(np.array(concat_plot_variable[:,i,:], dtype=float).flatten()))
 #print nan_mask
  print concat_plot_variable.shape
 
  try:    
   if nan_mask.mask.all() == False:
    print nan_mask
    xi,yi, zi = grid_data_cs(np.array(pressures[:,i,:], dtype=float).flatten(), np.repeat(distances, concat_plot_variable[:,i,:].shape[1]).flatten(), nan_mask) 
    xiw,yiw, ziw = grid_data_cs(np.array(pressures[:,i,:], dtype=float).flatten(), np.repeat(distances, concat_plot_variable[:,i,:].shape[1]).flatten(), wind_to_plot[nan_mask.mask])
    cont,cbar = plot_rad_cs_winds(xi, yi, zi, min_contour, max_contour, ziw)
    station_name_plot(station_list_cs, first_station, yi)
  except Exception, e:
   print e
  return cont,cbar

        
station_list_search='/nfs/a90/eepdw/Data/Observations/Radiosonde_downloaded_from_NOAA_GUAN/igra-stations.txt'
station_metadata=[]
f = open(station_list_search,'r')
for line in f:
     line = line.strip()
     line=re.sub(r'([A-Z])\s([A-Z])', r'\1_\2',line)
     line=re.sub(r'([A-Z])\s\s([A-Z])', r'\1_\2',line)
     station_metadata.append(line.split())
f.close()

first_station_title, first_station_lon, first_station_lat = station_info_search(first_station)

variable_list={'pressures': 0, 'temps':1, 'dewpoints':2, 'winddirs':3, 'windspeeds':4, 'pot_temp':5, 
               'sat_vap_pres':6, 'vap_press':7, 'rel_hum':8, 'wvmr':9, 'sp_hum':10, 'sat_temp':11, 'theta_e':12, 'theta_e_sat':13}

variable='pressures'

var_index = variable_name_index_match(variable, variable_list)
pressures, distances, date_min_max = variable_cat(var_index, station_list_cs)

variable='rel_hum'

var_index = variable_name_index_match(variable, variable_list)
concat_plot_variable, distances, date_min_max = variable_cat(var_index, station_list_cs)

variable='windspeeds'

var_index = variable_name_index_match(variable, variable_list)
wind_direction, distances, date_min_max = variable_cat(var_index, station_list_cs)

variable='winddirs'

var_index = variable_name_index_match(variable, variable_list)
wind_speed, distances, date_min_max = variable_cat(var_index, station_list_cs)

u_wind,v_wind = u_v_winds(wind_direction, wind_speed)

max_contour=100
min_contour=0
tick_interval=10
    
for i, date_bin in enumerate(date_min_max[:,0]):
        
 try:
  cont,cbar = date_bin_plot_wind(i, date_bin, concat_plot_variable, pressures, distances, min_contour, max_contour, v_wind)
        
  cbar.set_label('\%', rotation=90)

  print date_bin
  
  plt.title('%s %s Cross-Section of Relative Humidity from Radiosonde Soundings' % (date_bin.strftime("%d %B"), Cross_Section_Title.replace('_',' ') ))

  plt.show()
  
  #plt.savefig('/nfs/a90/eepdw/Figures/Radiosonde/Cross_Sections/%s_%s_%s_Relative_Humidity.png' % (Cross_Section_Title, date_bin.strftime("%y"), date_bin.strftime("%d_%B")),  format='png', bbox_inches='tight')
  plt.close()
  plt.clf()
  gc.collect()
 except Exception, e:
   print e
