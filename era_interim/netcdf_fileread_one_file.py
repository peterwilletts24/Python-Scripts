########################################
# Read TRMM netcdf files.py
#
# Created by: Peter Willetts
# Created on: 12/11/2013
#
########################################
#
#
###################################################
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
import glob
import re
import os
import pickle

import datetime
 
#first_month=8
#first_day_of_month=
#last_month=9
#last_day=

time_min=datetime.datetime(2011,8,18,0,0,0,0)
time_max=datetime.datetime(2011,9,8,0,0,0,0)

lon_max = 116 
lon_min = 34

lat_max= 40.
lat_min=-11.25

nc = Dataset('/nfs/a90/eepdw/Data/ERA_Iinterim_Heat_Rad_Fluxes/netcdf-web237-20140610105812-24189-67930.nc')

hours_since=datetime.datetime(1900,1,1,0,0,0,0)
# Get min and max  index positions for latitude and longitude

datetimes = np.array([datetime.timedelta(hours=int(i))+hours_since for i in nc.variables['time'][:]])

time_index= np.where((datetimes<=time_max) & (datetimes >= time_min))

la_index = np.where((nc.variables['latitude'][:]<=lat_max) & (nc.variables['latitude'][:] >= lat_min))
lo_index = np.where((nc.variables['longitude'][:]<=lon_max) & (nc.variables['longitude'][:] >= lon_min))

la_i_max = np.max(la_index)
la_i_min = np.min(la_index)
lo_i_max = np.max(lo_index)
lo_i_min = np.min(lo_index)

t_i_max = np.max(time_index)
t_i_min = np.min(time_index) 

lat_amounts=la_i_max-la_i_min
lon_amounts=lo_i_max-lo_i_min

print nc


latent_in = nc.variables['slhf'][t_i_min:t_i_max,la_i_min:la_i_max, lo_i_min:lo_i_max]
sensible_in = nc.variables['sshf'][t_i_min:t_i_max,la_i_min:la_i_max, lo_i_min:lo_i_max]
swave_in = nc.variables['ssrd'][t_i_min:t_i_max,la_i_min:la_i_max, lo_i_min:lo_i_max]
lwave_in = nc.variables['str'][t_i_min:t_i_max,la_i_min:la_i_max, lo_i_min:lo_i_max]

latitude_in = nc.variables['latitude'][la_index]
longitude_in = nc.variables['longitude'][lo_index]
time_in = datetimes[time_index]



# pickle.dump([sphum_dom, longitude_dom, latitude_dom, time_dom, time_hour], open('/nfs/a90/eepdw/Data/Saved_data/era_i/era_i_emb_time_update_large_sphum.p', 'wb'))
# pickle.dump([geopotential_dom, longitude_dom, latitude_dom, time_dom, time_hour], open('/nfs/a90/eepdw/Data/Saved_data/era_i/era_i_emb_time_update_large_geopotential.p', 'wb'))
# pickle.dump([u_wind_dom, longitude_dom, latitude_dom, time_dom, time_hour], open('/nfs/a90/eepdw/Data/Saved_data/era_i/era_i_emb_time_update_large_u_wind.p', 'wb'))
# pickle.dump([v_wind_dom, longitude_dom, latitude_dom, time_dom, time_hour], open('/nfs/a90/eepdw/Data/Saved_data/era_i/era_i_emb_time_update_large_v_wind.p', 'wb'))
# pickle.dump([cloud_cover_dom, longitude_dom, latitude_dom, time_dom, time_hour], open('/nfs/a90/eepdw/Data/Saved_data/era_i/era_i_emb_time_update_large_cloud_cover.p', 'wb'))
# pickle.dump([temperature_dom, longitude_dom, latitude_dom, time_dom, time_hour], open('/nfs/a90/eepdw/Data/Saved_data/era_i/era_i_emb_time_update_large_temperature.p', 'wb'))
        
if '__name__' == '__netcdf_fileread___':
  TRMM_fileread()                                    
