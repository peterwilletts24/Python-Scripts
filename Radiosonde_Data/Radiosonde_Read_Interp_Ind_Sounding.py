# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import datetime
import scipy.interpolate
import re
import numpy as np

# <codecell>

station_list_search='/nfs/a90/eepdw/Data/Observations/Radiosonde_downloaded_from_NOAA_GUAN/igra-stations.txt'
station_metadata=[]
f = open(station_list_search,'r')
for line in f:
     line = line.strip()
     line=re.sub(r'([A-Z])\s([A-Z])', r'\1_\2',line)
     line=re.sub(r'([A-Z])\s\s([A-Z])', r'\1_\2',line)
     station_metadata.append(line.split())
f.close()

# <codecell>

station_list_cs=[43150, 42867 43014, 42339, 40990, 40948]

# <codecell>

#Raw data read
# /nfs/a90/eepdw/Data/Observations/Radiosonde_downloaded_from_NOAA_GUAN/igra-stations.txt
# /nfs/a90/eepdw/Data/Observations/Radiosonde_downloaded_from_NOAA_GUAN/readme.txt

# Station numbers for Bombay - 43003; New Delhi - 42182, Lucknow - 42369, Cochin, Bangalore, Madras, Minicoy, Nagpur Sonegaon)
import datetime
import scipy.interpolate
import re
import numpy as np

date_min=datetime.datetime(2011,5,1,0,0,0)
date_max=datetime.datetime(2011,7,1,0,0,0)

import datetime


test_count=0

match_header = re.compile(r'#.....20110[56789]')

match_bad = re.compile(r'9999')
    
lev_count=False

y_points=np.linspace(5000, 100000, 200) # Points for pressure interpolation

def interp_sounding(variable, pressures_s,y_points):
    #print pressures_s
    #print variable
    
    #pressures_s,variable=np.sort(np.array((pressures_s, variable), dtype=float), axis=1)
    #a=np.array((pressures_s, variable))
   # a=a[a[0].argsort()]
    #print a
    #print pressures_s
    #print variable
    nan_mask = np.ma.masked_array(np.array(variable, dtype=float), np.isnan(np.array(variable, dtype=float)))
    nan_mask_p = np.ma.masked_array(np.array(pressures_s, dtype=float), np.isnan(np.array(variable, dtype=float)))
    variable = [x for (y,x) in sorted(zip(pressures_s, variable), key=lambda pair: pair[0])]
    pressures_s = [y for (y,x) in sorted(zip(pressures_s, variable), key=lambda pair: pair[0])]
    #print nan_mask_p
    #print nan_mask
    #print sorted_pre_interp
    interp = scipy.interpolate.interp1d(pressures_s, variable, bounds_error=False, fill_value=np.nan)
    y_interp = interp(y_points)
    return y_interp

for stat in station_list_cs:
  single_values=[]
  pressures=[]
  temps=[]
  dewpoints=[]
  winddirs=[]
  windspeeds=[]
    
  pot_temp=[]
  sat_vap_pres=[]
  vap_press=[]
  rel_hum=[]
  wvmr=[]
  sp_hum=[]
  sat_temp=[]
  theta_e=[]
  theta_e_sat=[]
    
  pressures_s=[]
  temps_s=[]
  dewpoints_s=[]
  winddirs_s=[]
  windspeeds_s=[]
    
  #pressures_for_plotting = np.array((dates_for_plotting, pressures, temps, dewpoints, winddirs, windspeeds))
    
  dates_for_plotting_single=[]
 
  levs=-1
  numlevs=0

  #print stat

  i = '/nfs/a90/eepdw/Data/Observations/Radiosonde_downloaded_from_NOAA_GUAN/%s.dat' % stat
  print i
  f = open(i,'r')
  #print f
  for line in f:
   #line = line.strip().replace('-99999',' %s ' %float('NaN'))
   line = line.strip()
   #print line
   #columns=line.split()
            
   if levs >= numlevs and pressures_s and temps_s and dewpoints_s and winddirs_s and windspeeds_s:
     lev_count=False
     #print pressures_s
     #print temps_s
    
     temps_s_mask = np.ma.masked_array(np.array(temps_s, dtype=float), np.isnan(np.array(temps_s, dtype=float)))
     pressures_s_mask =  np.ma.masked_array(np.array(pressures_s, dtype=float), np.isnan(np.array(pressures_s, dtype=float)))   
     dewpoints_s_mask = np.ma.masked_array(np.array(dewpoints_s, dtype=float), np.isnan(np.array(dewpoints_s, dtype=float)))
        
     temps_scent = np.array(temps_s_mask, dtype=float)/10
     temps_s = (np.array(temps_s_mask, dtype=float)/10)+273.15
     pressures_s_hPa = np.array(pressures_s_mask, dtype=float)/100 # Pressures in raw files as mb*100 - which is Pa - convert to hPa for calcs
     dewpoints_s = np.array(dewpoints_s_mask, dtype=float)/10
        
     dewpoint_temp = temps_scent-dewpoints_s
    
     pot_temp_s = temps_s*((1000/pressures_s_hPa)**(2/7))

     sat_vap_pres_s = 611.2*np.exp(17.67*(temps_scent/(temps_scent+243.5)))
     
     vap_press_s = 611.2*np.exp(17.67*(dewpoint_temp/(dewpoint_temp+243.5)))
        
     rel_hum_s = 100*vap_press_s/sat_vap_pres_s

     wvmr_s = 0.622*vap_press_s/(pressures_s-vap_press_s) # kg/kg
        
     sp_hum_s = wvmr_s/(1+wvmr_s)
        
     sat_temp_s = 55+2840/(3.5*np.log(temps_s)-np.log(vap_press_s/100)-4.805)
    
     theta_e_s = pot_temp_s*((1000/pressures_s_hPa)**(0.2854*(1-0.28*wvmr_s)))*np.exp(((3376/sat_temp_s)-2.54)*wvmr_s*(1+0.81*wvmr_s))
        

     wvmr_sat_s = 0.622*sat_vap_pres_s/(pressures_s-sat_vap_pres_s) # kg/kg
    
     sat_temp_sat_s = 55+2840/(3.5*np.log(temps_s)-np.log(sat_vap_pres_s/100)-4.805)
        
     theta_e_sat_s = pot_temp_s*((1000/pressures_s_hPa)**(0.2854*(1-0.28*wvmr_sat_s)))*np.exp(((3376/sat_temp_sat_s)-2.54)*wvmr_sat_s*(1+0.81*wvmr_sat_s))
     
     try:
      temp_interp = interp_sounding(temps_s, pressures_s,y_points)
            
      #print pressures_s
      dewpoints_interp = interp_sounding(dewpoints_s, pressures_s,y_points)
      #print dewpoints_interp
      winddirs_interp= interp_sounding(np.array(winddirs_s, dtype=float), pressures_s, y_points)
      windspeeds_interp= interp_sounding(np.array(windspeeds_s, dtype=float)/10, pressures_s, y_points)
            
      pot_temp_interp = interp_sounding(pot_temp_s, pressures_s,y_points)
      #print pot_temp_interp
      sat_vap_pres_interp = interp_sounding(sat_vap_pres_s, pressures_s,y_points)
      vap_press_interp = interp_sounding(vap_press_s, pressures_s,y_points)
      rel_hum_interp = interp_sounding(rel_hum_s, pressures_s,y_points)
      #print rel_hum_interp
      wvmr_interp = interp_sounding(wvmr_s, pressures_s,y_points)
      sp_hum_interp = interp_sounding(sp_hum_s, pressures_s,y_points)
      sat_temp_interp = interp_sounding(sat_temp_s, pressures_s,y_points)
      theta_e_interp = interp_sounding(theta_e_s, pressures_s,y_points)
            
      theta_e_sat_interp = interp_sounding(theta_e_sat_s, pressures_s, y_points)
            
      pressures.append(list(y_points))
      temps.append(list(temp_interp))
      dewpoints.append(list(dewpoints_interp))
      winddirs.append(list(winddirs_interp))
      windspeeds.append(list(windspeeds_interp))
            
      pot_temp.append(list(pot_temp_interp))
      sat_vap_pres.append(list(sat_vap_pres_interp))
      vap_press.append(list(vap_press_interp))
      rel_hum.append(list(rel_hum_interp))
      wvmr.append(list(wvmr_interp))
      sp_hum.append(list(sp_hum_interp))
      sat_temp.append(list(sat_temp_interp))
      theta_e.append(list(theta_e_interp))
            
      theta_e_sat.append(list(theta_e_sat_interp))
    
      dates_for_plotting_single.append(date)
        
     # dates_for_plotting.append(list([dates_for_plotting_single for i in range(len(y_points))]))
      
    
      pp=pressures_s
      tt=temps_s
      ppp=pressures
      ttt=temps
      ti=temp_interp
     
      del pressures_s
      del temps_s
      del dewpoints_s
      del winddirs_s
      del windspeeds_s
      del temps_s_mask
      del pressures_s_mask 
      del dewpoints_s_mask 
      del temps_scent
      del pressures_s_hPa
      del dewpoint_temp
      del pot_temp_s
      del sat_vap_pres_s
      del vap_press_s 
      del rel_hum_s
      del wvmr_s
      del sp_hum_s
      del sat_temp_s
      del theta_e_s
        
      
      pressures_s=[]
      temps_s=[]
      dewpoints_s=[]
      winddirs_s=[]
      windspeeds_s=[]
        
      temps_s_mask=[]
      pressures_s_mask =[]
      dewpoints_s_mask =[]
        
      temps_scent = []
      temps_s = []
      pressures_s_hPa = []
      dewpoints_s = []
        
      dewpoint_temp = []
    
      pot_temp_s = []

      sat_vap_pres_s = []
     
      vap_press_s = []
        
      rel_hum_s = []

      wvmr_s =[]
        
      sp_hum_s = []
        
      sat_temp_s =[]
    
      theta_e_s = []
        
      theta_e_sat_s = []
        
     except ValueError, e:
      print e
      pass
    # np.append(pressures_for_plotting, (y_points.tolist(), temp_interp.tolist(), dewpoints_interp.tolist(), winddirs_interp.tolist(), windspeeds_interp.tolist()), axis=1)
  # levs=0
   if lev_count is True: 
    
    # Line length varies depending on length of first variable (pressure)
 #   li=len(line)
  # ld=142-len(line)
                      
    pressures_s.append(float(line[2:8].replace('-9999','%s' %float('NaN')).replace('-8888','%s' %float('NaN'))))     
    temps_s.append(float(line[15:20].replace('-9999','%s' %float('NaN')).replace('-8888','%s' %float('NaN')))) 
    dewpoints_s.append(float(line[21:26].replace('-9999','%s' %float('NaN')).replace('-8888','%s' %float('NaN')))) 
    winddirs_s.append(float(line[26:31].replace('-9999','%s' %float('NaN')).replace('-8888','%s' %float('NaN')))) 
    windspeeds_s.append(float(line[31:36].replace('-9999','%s' %float('NaN')).replace('-8888','%s' %float('NaN')))) 
    
    
    levs+=1
    
    #print pressures_s
    
# Extract independent headers 

   #print uwind
   if match_header.search(line) is not None:
    #print line
    
    year=int(line[6:10])
    month=int(line[10:12])
    day=int(line[12:14])
    hour=int(line[14:16])
    
    date = datetime.datetime(year,month,day,hour)
    
    numlevs=int(line[20:24])
    
    #print numlevs
    
    lev_count=True
    levs=0
    
    pressures_s=[]
    temps_s=[]
    dewpoints_s=[]
    winddirs_s=[]
    windspeeds_s=[]
    
  #print pressures.shape  
         
      
  pressures_for_plotting = np.array((pressures, temps, dewpoints, winddirs, windspeeds, pot_temp, sat_vap_pres, vap_press, rel_hum, wvmr, sp_hum, sat_temp,theta_e, theta_e_sat))
  print pressures_for_plotting.shape
  np.savez('/nfs/a90/eepdw/Data/Observations/Radiosonde_Numpy/Radiosonde_Single_Station_PRESSURES__IND_INTERP_SOUNDING_%s_%s_to_%s'\
          % (stat, date_min.strftime('%Y%m%d'), date_max.strftime('%Y%m%d') ), pressures_for_plotting=pressures_for_plotting, dates_for_plotting_single=dates_for_plotting_single)
 

# <codecell>


# Loading individual stations numpy files, sorting into datetime bins
# AND then interpolating for each station 
# BEFORE appending so that all 1-D interpolated values for stations are in one array with distances for gridding and plotting
# Plot variable index must be changed in script (plot_variable)

#Indexes -  pressures = 1, ptemps = 2, vptemps = 3, rh = 4, uwind = 5, vwind = 6, calculated geopotential heights = 7 , 
#temps=8, actual vapour pressure=9, wvmr = 10, theta_e = 11, theta_es = 12. . . may change
from dateutil.relativedelta import relativedelta

date_min=datetime.datetime(2011,5,1,0,0,0)
date_max=datetime.datetime(2011,10,1,0,0,0)

#delta = relativedelta(months=+1)
delta = relativedelta(weeks=+1)

Cross_Section_Title = 'Vizag_to_Afghanistan'

from dateutil.relativedelta import relativedelta
from math import sin, cos, atan2, radians, sqrt
from collections import defaultdict
import collections
import bisect

#from scipy.interpolate import spline

#import scipy.interpolate.interp1d

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

all_stat_date_range_pressure_interp=[]
grid=[]

d = date_min
while (d <= date_max): 
 grid.append(d)
 d += delta

#print grid
                 
bins=collections.defaultdict(list)

for stat in station_list_cs: 
  date_bin_mean_all_dates_one_station=[]
  min_max_date_bin=[]
  bins=collections.defaultdict(list)
        
# Pressure Levels

 #pressure_file = np.load('/nfs/a90/eepdw/Data/Observations/Radiosonde_Numpy/Radiosonde_Single_Station_PRESSURES_%s_%s.npy' % (stat, date_title))
  npz_file = np.load('/nfs/a90/eepdw/Data/Observations/Radiosonde_Numpy/Radiosonde_Single_Station_PRESSURES__IND_INTERP_SOUNDING_%s_%s_to_%s.npz'\
          % (stat, date_min.strftime('%Y%m%d'), date_max.strftime('%Y%m%d') ))
  pressure_file=npz_file['pressures_for_plotting']
  date_file=npz_file['dates_for_plotting_single']
  #print '/nfs/a90/eepdw/Data/Observations/Radiosonde_Numpy/Radiosonde_Single_Station_PRESSURES_%s_%s_to_%s.npy'\
  #        % (stat, date_min.strftime('%Y%m%d'), date_max.strftime('%Y%m%d') )

  #print pressure_file
  for di, dates in enumerate(date_file):
   #print idx

   idx=bisect.bisect(grid,dates)
   bins[idx].append((pressure_file[:,di],dates))
   #print idx
  #print len(bins)
  #print len(grid)
        
  for i in range(len(grid)):
    #if bins[l] is None:
   #bins[l].append(([],[]))
        
   if np.array(bins[i]).size != 0:
            
    empty_array =  np.empty((np.nanmean(np.dstack(np.array(bins[i])[:,0]), axis=2).shape))
    empty_array[:] = np.NAN
    
    empty_array_list =  (np.NAN, np.NAN)
    #empty_array_list[:] = np.NAN
    
    
  for i in range(len(grid)):
    #if bins[l] is None:
   #bins[l].append(([],[]))
        
   if np.array(bins[i]).size != 0:
  #for i in bins:
    print i
   #print np.array(bins[i])[:,0].shape
  #bins = np.sort(bins, axis=1)
   #print grid[i]   
   #date_bin_mean = np.nanmean(np.dstack(np.array(bins[i])[:,0]), axis=2)
   #print date_bin_mean
   #nan_mask = np.ma.masked_array(np.array(date_bin_mean, dtype=float), np.isnan(np.array(date_bin_mean, dtype=float)))
            
  # if np.all(np.isnan(np.dstack(np.array(bins[i])[:,0]))) == False:
        

    date_bin_mean_all_dates_one_station.append(np.nanmean(np.dstack(np.array(bins[i])[:,0]), axis=2))
   #print date_bin_mean
    min_max_date_bin.append((min(np.array(bins[i])[:,1]), max(np.array(bins[i])[:,1])))
    
   elif np.array(bins[i]).size == 0:
     date_bin_mean_all_dates_one_station.append(empty_array)
     min_max_date_bin.append(empty_array_list)
    
  #print date_bin_mean_all_dates_one_station

  print np.array(date_bin_mean_all_dates_one_station).shape
  np.savez('/nfs/a90/eepdw/Data/Observations/Radiosonde_Numpy/Radiosonde_Cross_Section_IND_SOUNDING_INTERP_MEAN_%s_%s_%s_%s_%s' \
       % (Cross_Section_Title, date_min.strftime('%Y%m%d'), date_max.strftime('%Y%m%d'), delta, stat)\
    , date_bin_mean_all_dates_one_station=date_bin_mean_all_dates_one_station, min_max_date_bin=min_max_date_bin)

# <codecell>

nan_mask.mask.all()

# <codecell>

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


Cross_Section_Title = 'Vizag_to_Afghanistan'
station_list_cs=[43150, 43014, 42867, 42339, 40948, 40990]

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
    for stat in station_list_cs:
     load_file = np.load('/nfs/a90/eepdw/Data/Observations/Radiosonde_Numpy/Radiosonde_Cross_Section_'
                        'IND_SOUNDING_INTERP_MEAN_%s_%s_%s_%s_%s.npz' 
                        % (Cross_Section_Title, date_min.strftime('%Y%m%d'), date_max.strftime('%Y%m%d'), delta, stat))
    #print stat
     print load_file['date_bin_mean_all_dates_one_station'].shape
        
     station_title, station_lon, station_lat = station_info_search(stat)
    
     dist_from_first_station = calculate_distance_from_first_station(stat, first_station_lon, first_station_lat, station_lat, station_lon)
    
     print dist_from_first_station
     #print load_file['date_bin_mean_all_dates_one_station'][:,var_index,:].shape
     var_cat.append(load_file['date_bin_mean_all_dates_one_station'][:,var_index,:])
     distances.append(dist_from_first_station)
    print np.array(var_cat).shape           
    return np.array(var_cat), np.array(distances, dtype=float), load_file['min_max_date_bin']

def station_name_plot(station_list_cs, first_station, yi):    
 y_offset_text=0
 first_station_title, first_station_lon, first_station_lat = station_info_search(first_station)
    
 
 for stat in station_list_cs: 
     
   station_title, station_lon, station_lat = station_info_search(stat)
    
   dist_from_first_station = calculate_distance_from_first_station(stat, first_station_lon, first_station_lat, station_lat, station_lon)
    
   plt.axvline(x=dist_from_first_station, ymin=0, ymax=1, label=station_title, color='k')
   plt.text(dist_from_first_station+0.1,max(yi)/100+20,station_title,rotation=-45)

   y_offset_text=+1   
    
def grid_data_cs(pressure, distance, param):
 xi=np.linspace(0, max(distance), 1000) 
 #yi=np.linspace(np.nanmin(pressure), np.nanmax(pressure), 500) 
 yi=np.linspace(5000, 100000, 50) # Points for pressure interpolation

 #yi=np.array([1000, 925, 850, 700, 500, 400, 300, 250, 200, 150, 100, 70, 50, 30, 20,10], dtype=float)
 #yi=np.array([10, 20, 30, 50, 70, 100, 150, 200, 250, 300, 400, 500, 700, 850, 925, 1000]*100, dtype=float)
 
 #zi = ml.griddata(distance, pressure,param,xi, yi, interp='nn')
 try:
    zi = scipy.interpolate.griddata((distance, pressure),  param, (xi[None,:],yi[:,None]), method='linear')
 except Exception, e:
  print e
 return xi,yi,zi 
 #return xi,yi
def plot_rad_cs(xi,yi,zi, min_contour, max_contour):
 clevs = np.linspace(min_contour, max_contour,256)
 ticks = (np.arange(min_contour, max_contour+tick_interval,tick_interval))

 plt.figure(figsize=(14,8))
 cmap=plt.cm.jet
 cont = plt.contourf(xi,yi/100, zi, clevs, cmap=cmap, extend='both')

 cbar = plt.colorbar(cont, orientation='vertical', pad=0.05, extend='both', format = '$%d$')
 #cbar.set_label('$W m^{-2}$') 
 cbar.set_ticks(np.arange(min_contour, max_contour+tick_interval,tick_interval))
 cbar.set_ticklabels(['${%d}$' % i for i in ticks])
    
 plt.gca().invert_yaxis()
 plt.ylabel('Pressure (hPa)')
 plt.xlabel('km from first station')
    
 return cont,cbar

def date_bin_plot(i, date_bin, concat_plot_variable, pressures, distances, min_contour, max_contour):
  nan_mask = np.ma.masked_array(np.array(concat_plot_variable[:,i,:], dtype=float).flatten(), np.isnan(np.array(concat_plot_variable[:,i,:], dtype=float).flatten()))
 #print nan_mask
  print concat_plot_variable.shape
 
  try:    
   if nan_mask.mask.all() == False:
    xi,yi, zi = grid_data_cs(np.array(pressures[:,i,:], dtype=float).flatten(), np.repeat(distances, concat_plot_variable[:,i,:].shape[1]).flatten(), nan_mask) 
    
    conts,cbars = plot_rad_cs(xi, yi, zi, min_contour, max_contour)
    station_name_plot(station_list_cs, first_station, yi)
  except Exception, e:
   print e
  return conts,cbars
        
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

max_contour=100
min_contour=0
tick_interval=10
    
for i, date_bin in enumerate(date_min_max[:,0]):
        
 try:
  conts,cbars = date_bin_plot(i, date_bin, concat_plot_variable, pressures, distances, min_contour, max_contour)
        
  cbars.set_label('\%', rotation=90)
  plt.title('%s %s Cross-Section of Relative Humidity from Radiosonde Soundings' % (date_bin.strftime("%d %B"), Cross_Section_Title.replace('_',' ') ))

  #plt.show()
  plt.savefig('/nfs/a90/eepdw/Figures/Radiosonde/Cross_Sections/%s_%s_%s_Relative_Humidity.png' % (Cross_Section_Title, date_bin.strftime("%y"), date_bin.strftime("%d_%B")),  format='png', bbox_inches='tight')
 
 except Exception, e:
   print e

#Potential Temperature

variable='pot_temp'

var_index = variable_name_index_match(variable, variable_list)
concat_plot_variable, distances, date_min_max = variable_cat(var_index, station_list_cs)

max_contour=300
min_contour=250
tick_interval=10
    
for i, date_bin in enumerate(date_min_max[:,0]):
        
 try:
  conts,cbars = date_bin_plot(i, date_bin, concat_plot_variable, pressures, distances, min_contour, max_contour)
        
  cbars.set_label('\%', rotation=90)
  plt.title('%s %s Cross-Section of from Radiosonde Soundings' % (date_bin.strftime("%d %B"), Cross_Section_Title.replace('_',' ') ))

  #plt.show()
  plt.savefig('/nfs/a90/eepdw/Figures/Radiosonde/Cross_Sections/%s_%s_%_.png' % (Cross_Section_Title, date_bin.strftime("%y"), date_bin.strftime("%d_%B")),  format='png', bbox_inches='tight')
 
 except Exception, e:
   print e


#Theta E

variable='theta_e'

var_index = variable_name_index_match(variable, variable_list)
concat_plot_variable, distances, date_min_max = variable_cat(var_index, station_list_cs)

max_contour=380
min_contour=330
tick_interval=10
    
for i, date_bin in enumerate(date_min_max[:,0]):
        
 try:
  conts,cbars = date_bin_plot(i, date_bin, concat_plot_variable, pressures, distances, min_contour, max_contour)
        
  cbars.set_label('K', rotation=90)
  plt.title('%s %s Cross-Section of from Radiosonde Soundings' % (date_bin.strftime("%d %B"), Cross_Section_Title.replace('_',' ') ))

  #plt.show()
  plt.savefig('/nfs/a90/eepdw/Figures/Radiosonde/Cross_Sections/%s_%s_%_.png' % (Cross_Section_Title, date_bin.strftime("%y"), date_bin.strftime("%d_%B")),  format='png', bbox_inches='tight')
 
 except Exception, e:
   print e

#Theta Es

variable='theta_e_sat'

var_index = variable_name_index_match(variable, variable_list)
concat_plot_variable, distances, date_min_max = variable_cat(var_index, station_list_cs)

max_contour=380
min_contour=330
tick_interval=10
    
for i, date_bin in enumerate(date_min_max[:,0]):
        
 try:
  conts,cbars = date_bin_plot(i, date_bin, concat_plot_variable, pressures, distances, min_contour, max_contour)
        
  cbars.set_label('K', rotation=90)
  plt.title('%s %s Cross-Section of from Radiosonde Soundings' % (date_bin.strftime("%d %B"), Cross_Section_Title.replace('_',' ') ))

  #plt.show()
  plt.savefig('/nfs/a90/eepdw/Figures/Radiosonde/Cross_Sections/%s_%s_%_.png' % (Cross_Section_Title, date_bin.strftime("%y"), date_bin.strftime("%d_%B")),  format='png', bbox_inches='tight')
 
 except Exception, e:
   print e
        
        
# WVMR

variable = 'wvmr'

var_index = variable_name_index_match(variable, variable_list)
concat_plot_variable, distances, date_min_max = variable_cat(var_index, station_list_cs)

max_contour=5
min_contour=0
tick_interval=1
    



#Theta Es - Theta E - Buoyancy

variable='theta_e'

var_index = variable_name_index_match(variable, variable_list)
concat_plot_variable_theta_e, distances, date_min_max = variable_cat(var_index, station_list_cs)

variable='theta_e_sat'

var_index = variable_name_index_match(variable, variable_list)
concat_plot_variable_theta_e_sat, distances, date_min_max = variable_cat(var_index, station_list_cs)

concat_plot_variable = concat_plot_variable_theta_e_sat - concat_plot_variable_theta_e

max_contour=50
min_contour=0
tick_interval=10
    
for i, date_bin in enumerate(date_min_max[:,0]):
        
 try:
  conts,cbars = date_bin_plot(i, date_bin, concat_plot_variable, pressures, distances, min_contour, max_contour)
        
  cbars.set_label('K', rotation=90)
  plt.title('%s %s Cross-Section of from Radiosonde Soundings' % (date_bin.strftime("%d %B"), Cross_Section_Title.replace('_',' ') ))

  #plt.show()
  plt.savefig('/nfs/a90/eepdw/Figures/Radiosonde/Cross_Sections/%s_%s_%_.png' % (Cross_Section_Title, date_bin.strftime("%y"), date_bin.strftime("%d_%B")),  format='png', bbox_inches='tight')
 
 except Exception, e:
   print e

# <codecell>

import numpy as np
if np.var(nan_mask)

# <codecell>

#Weekly

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.mlab as ml

import datetime

from dateutil.relativedelta import relativedelta

import re

import numpy as np

from math import sin, cos, atan2, radians, sqrt

import scipy.interpolate


Cross_Section_Title = 'Vizag_to_Afghanistan'
station_list_cs=[43150, 43014, 42867, 42339, 40948, 40990]

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
    for stat in station_list_cs:
     load_file = np.load('/nfs/a90/eepdw/Data/Observations/Radiosonde_Numpy/Radiosonde_Cross_Section_'
                        'IND_SOUNDING_INTERP_MEAN_%s_%s_%s_%s_%s.npz' 
                        % (Cross_Section_Title, date_min.strftime('%Y%m%d'), date_max.strftime('%Y%m%d'), delta, stat))
    #print stat
     print load_file['date_bin_mean_all_dates_one_station'].shape
        
     station_title, station_lon, station_lat = station_info_search(stat)
    
     dist_from_first_station = calculate_distance_from_first_station(stat, first_station_lon, first_station_lat, station_lat, station_lon)
    
     print dist_from_first_station
        
     #print load_file['date_bin_mean_all_dates_one_station'][0]
     #if load_file['date_bin_mean_all_dates_one_station'][:,var_index,:].size!=0:
            
      #var_cat.append(load_file['date_bin_mean_all_dates_one_station'][:,var_index,:])
     
     var_cat.append(load_file['date_bin_mean_all_dates_one_station'][:,var_index,:])
     distances.append(dist_from_first_station)
    print np.array(var_cat).shape           
    return np.array(var_cat), np.array(distances, dtype=float), load_file['min_max_date_bin']

def station_name_plot(station_list_cs, first_station, yi):    
 y_offset_text=0
 first_station_title, first_station_lon, first_station_lat = station_info_search(first_station)
    
 
 for stat in station_list_cs: 
     
   station_title, station_lon, station_lat = station_info_search(stat)
    
   dist_from_first_station = calculate_distance_from_first_station(stat, first_station_lon, first_station_lat, station_lat, station_lon)
    
   plt.axvline(x=dist_from_first_station, ymin=0, ymax=1, label=station_title, color='k')
   plt.text(dist_from_first_station+0.1,max(yi)/100+20,station_title,rotation=-45)

   y_offset_text=+1   
    
def grid_data_cs(pressure, distance, param):
 xi=np.linspace(0, max(distance), 1000) 
 #yi=np.linspace(np.nanmin(pressure), np.nanmax(pressure), 500) 
 yi=np.linspace(5000, 100000, 50) # Points for pressure interpolation

 #yi=np.array([1000, 925, 850, 700, 500, 400, 300, 250, 200, 150, 100, 70, 50, 30, 20,10], dtype=float)
 #yi=np.array([10, 20, 30, 50, 70, 100, 150, 200, 250, 300, 400, 500, 700, 850, 925, 1000]*100, dtype=float)
 
 #zi = ml.griddata(distance, pressure,param,xi, yi, interp='nn')
 zi = scipy.interpolate.griddata((distance, pressure),  param, (xi[None,:],yi[:,None]), method='linear')
 return xi,yi,zi 
 #return xi,yi
def plot_rad_cs(xi,yi,zi, min_contour, max_contour):
 clevs = np.linspace(min_contour, max_contour,256)
 ticks = (np.arange(min_contour, max_contour+tick_interval,tick_interval))

 plt.figure(figsize=(14,8))
 cmap=plt.cm.jet
 cont = plt.contourf(xi,yi/100, zi, clevs, cmap=cmap, extend='both')

 cbar = plt.colorbar(cont, orientation='vertical', pad=0.05, extend='both', format = '$%d$')
 #cbar.set_label('$W m^{-2}$') 
 cbar.set_ticks(np.arange(min_contour, max_contour+tick_interval,tick_interval))
 cbar.set_ticklabels(['${%d}$' % i for i in ticks])
    
 plt.gca().invert_yaxis()
 plt.ylabel('Pressure (hPa)')
 plt.xlabel('km from first station')
    
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

for i, date_bin in enumerate(date_min_max):
 #nan_mask = np.ma.masked_array(np.array(concat_plot_variable[:,i,:], dtype=float).flatten(), np.isnan(np.array(concat_plot_variable[:,i,:], dtype=float).flatten()))
 #print nan_mask
 #print concat_plot_variable.shape
 #try:    
  nan_mask = np.ma.masked_array(np.array(concat_plot_variable[:,i,:], dtype=float).flatten(), np.isnan(np.array(concat_plot_variable[:,i,:], dtype=float).flatten()))
  nan_mask_pr = np.ma.masked_array(np.array(pressures[:,i,:], dtype=float).flatten(), np.isnan(np.array(pressures[:,i,:], dtype=float).flatten()))
  
  if nan_mask.mask.all() == False:
   print date_bin
   print np.max(nan_mask)
   xi,yi, zi = grid_data_cs(nan_mask_pr, np.repeat(distances, concat_plot_variable[:,i,:].shape[1]).flatten(), nan_mask) 

   max_contour=100
   min_contour=0
   tick_interval=10
   cont,cbar = plot_rad_cs(xi, yi, zi, min_contour, max_contour)
   station_name_plot(station_list_cs, first_station, yi)
        
   cbar.set_label('\%', rotation=90)
   plt.title('%s %s Cross-Section of Relative Humidity from Radiosonde Soundings' % (date_bin[:,0].strftime("%d %B"), Cross_Section_Title.replace('_',' ') ))

  #plt.show()
   plt.savefig('/nfs/a90/eepdw/Figures/Radiosonde/Cross_Sections/%s_%s_%s_Relative_Humidity.png' % (Cross_Section_Title, year, date_bin.strftime("%d_%B")),  format='png', bbox_inches='tight')
 
# except Exception, e:
 # print e
 # print 
 # r=nan_mask
  








# <codecell>

nan_mask_pr.mask

# <codecell>

np.array(pressures[:,i,:], dtype=float).flatten()

# <codecell>

np.repeat(distances, concat_plot_variable[:,i,:].shape[1]).flatten()

# <codecell>

xi,yi, zi = grid_data_cs(np.array(pressures[:,i,:], dtype=float).flatten(), np.repeat(distances, concat_plot_variable[:,i,:].shape[1]).flatten(), nan_mask) 

# <codecell>

np.isnan(np.array(pressures[:,i,:], dtype=float).flatten())

# <markdowncell>

# * Madras.  43279  
# * Machilipatnam. 43185
# * Visakhpatnam. 43150
# * Hyderabad. 43128
# * Aurangabad. 43014
# * Nagpur. 42867
# * Jodhpur. 42339
# * (Karachi). 41780
# * Srinagar. 42027

# <headingcell level=2>

# Visakhpatnam to Afghanistan cross-section

# <codecell>

# Vizag to Afghanistan
station_list_cs=[43150, 43014, 42867, 42339, 40948, 40990]
#station_list=[42027]
#Indexes -  pressures = 1, ptemps = 2, vptemps = 3, rh = 4, uwind = 5, vwind = 6, calculated geopotential heights = 7 . . . may change
# Single level indexes - dates_for_plotting_single = 0, pwater = 1, lclpress = 2, lfcpress = 3, lnbpress = 4, cape = 5,cin = 6

# <headingcell level=2>

# Madras to Srinagar cross-section

# <codecell>

#Madras to Srinagar
station_list_cs=[43279, 43185, 43150, 43128, 43014, 42867, 42339, 41780, 42027]
#station_list=[42027]

# Single level indexes - dates_for_plotting_single = 0, pwater = 1, lclpress = 2, lfcpress = 3, lnbpress = 4, cape = 5,cin = 6

# <headingcell level=2>

# Save in cross-section numpy files

# <markdowncell>

# Change 'Cross_section_title' for different cross-section files

# <codecell>

def grid_data_cs(pressure, distance, param):
 xi=np.linspace(0, max(distance), 1000) 
 yi=np.linspace(min(pressure), max(pressure), 25) 
 #yi=np.array([1000, 925, 850, 700, 500, 400, 300, 250, 200, 150, 100, 70, 50, 30, 20,10], dtype=float)
 #yi=np.array([10, 20, 30, 50, 70, 100, 150, 200, 250, 300, 400, 500, 700, 850, 925, 1000]*100, dtype=float)
 zi = ml.griddata(distance, pressure,param,xi, yi, interp='nn')
 #zi = scipy.interpolate.griddata((matplotlib.dates.date2num(q[0]),q[1]), param, (xi[None,:],yi[:,None]), method='linear')
 return xi,yi,zi        

# <headingcell level=1>

# Weekly Radiosonde Cross-Sections for June

