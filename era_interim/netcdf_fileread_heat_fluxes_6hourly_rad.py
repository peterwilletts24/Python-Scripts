########################################
# Read ECMWF netcdf files for heat fluxes
#
# Created by: Peter Willetts
# Created on: 12/06/2014
#
# ECMWF heat and radiation flux - Read from netcdf 
#                                 filter by date, latitude and longitude
#                                 calculate mean and total heat flux
# BEWARE!!! ECMWF flux descrpitions may be the wrong way round, as well as upwards/downwards signs
# This script is deisgned to work with the total accumulated time-integrated fluxes at 0 timesteps - every 12 hours
# So average of 12 hourly accumulations in J/m2, divided by seconds, minutes and 12 hours gives Wm^-2
#

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

nc = Dataset('/nfs/a90/eepdw/Data/ERA_Iinterim_Heat_Rad_Fluxes/era_interim_netcdf_heat_rad_flux_evap_precip_6hr_timestep.nc')

hours_since=datetime.datetime(1900,1,1,0,0,0,0)
# Get min and max  index positions for latitude and longitude

datetimes = np.array([datetime.timedelta(hours=float(i))+hours_since for i in nc.variables['time'][:]])

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


#latent_in = nc.variables['slhf'][t_i_min:t_i_max+1,la_i_min:la_i_max+1, lo_i_min:lo_i_max+1]/(60*60*12)
#sensible_in = nc.variables['sshf'][t_i_min:t_i_max+1,la_i_min:la_i_max+1, lo_i_min:lo_i_max+1]/(60*60*6)
lwave_in = nc.variables['str'][t_i_min:t_i_max+1,la_i_min:la_i_max+1, lo_i_min:lo_i_max+1]/(60*60*6)
swave_in = nc.variables['ssr'][t_i_min:t_i_max+1,la_i_min:la_i_max+1, lo_i_min:lo_i_max+1]/(60*60*6)

latitude_in = nc.variables['latitude'][la_index]
longitude_in = nc.variables['longitude'][lo_index]
time_in = datetimes[time_index]

##

#latent_mean = -np.mean(latent_in, axis=0, dtype=np.float64)
#sensible_mean = -np.mean(sensible_in, axis=0, dtype=np.float64)

swave_mean = np.mean(swave_in, axis=0, dtype=np.float64)
lwave_mean = np.mean(lwave_in, axis=0, dtype=np.float64)

# I don't think the ECMWF data is very well documented
#  According to ECMWF descriptions 'swave_in' is solar (longwave - sounds wrong to me) downward radiation
#                                 'lwave_in is thermal (shortwave - again sounds wrong) upward radiation
#                                 'latent_mean' is upward 
#                                 'sensible_mean' is upward

# From UM calc - pcubetotal = Downward shortwave + Downward longwave flux - Upward sensible - Upward latent heat

# Even though the ECMRWF latent/sensible heat flux is 'upwards', the sign's are opposite to those in the EMBRACE data etc

#total_mean = swave_mean + lwave_mean - sensible_mean - latent_mean

#np.save('/nfs/a90/eepdw/Data/Saved_data/era_i/Fluxes/era_i_emb_latent_mean', latent_mean)
#np.save('/nfs/a90/eepdw/Data/Saved_data/era_i/Fluxes/era_i_emb_sensible_mean', sensible_mean)
np.save('/nfs/a90/eepdw/Data/Saved_data/era_i/Fluxes/era_i_emb_swave_mean', swave_mean)
np.save('/nfs/a90/eepdw/Data/Saved_data/era_i/Fluxes/era_i_emb_lwave_mean', lwave_mean)
#np.save('/nfs/a90/eepdw/Data/Saved_data/era_i/Fluxes/era_i_emb_total_mean', total_mean)
#np.save('/nfs/a90/eepdw/Data/Saved_data/era_i/Fluxes/era_i_emb_lats', latitude_in)
#np.save('/nfs/a90/eepdw/Data/Saved_data/era_i/Fluxes/era_i_emb_longs', longitude_in)



#np.save('/nfs/a90/eepdw/Data/Saved_data/era_i/Fluxes/era_i_emb_times_6h', time_in)
      
if '__name__' == '__netcdf_fileread___':
  TRMM_fileread()                                    
