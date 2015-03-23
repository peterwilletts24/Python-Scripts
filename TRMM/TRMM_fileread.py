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
#from mpl_toolkits.basemap import Basemap
#import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
import glob
import re
import os
import pickle

import pdb
 
#first_month=8
#first_day_of_month=
#last_month=9
#last_day=

lon_max = 116 
lon_min = 30.5

lat_max= 40
lat_min=-11.25

#Get number of lines 
no_lines=0
with open('/nfs/see-fs-01_users/eepdw/python_scripts/filenamelist/trmm_netcdffilelist') as f:
 for i, l in enumerate(f):

    # if . . .
   no_lines=no_lines+1


# Get min and max  index positions for latitude and longitude

linestr=l.rstrip()
nc = Dataset(linestr)

la_index = np.where((nc.variables['latitude'][:]<=lat_max) & (nc.variables['latitude'][:] >= lat_min))
lo_index = np.where((nc.variables['longitude'][:]<=lon_max) & (nc.variables['longitude'][:] >= lon_min))

la_i_max = np.max(la_index)
la_i_min = np.min(la_index)
lo_i_max = np.max(lo_index)
lo_i_min = np.min(lo_index)
   
lat_amounts=la_i_max-la_i_min
lon_amounts=lo_i_max-lo_i_min

# Pre-declare arrays

pcp_in = np.zeros((lat_amounts,lon_amounts),dtype=float)
latitude_in = np.zeros((lat_amounts),dtype=float)
longitude_in = np.zeros((lon_amounts),dtype=float)

pcp_dom = np.zeros((no_lines,lat_amounts,lon_amounts),dtype=float)
latitude_dom = np.zeros((no_lines,lat_amounts),dtype=float)
longitude_dom =  np.zeros((no_lines,lon_amounts),dtype=float)
time_dom = np.empty((no_lines), dtype=(str,12))
time_hour = np.empty((no_lines), dtype=(str,2))

 #fs=glob.iglob('/nfs/a80/earceb/model_data/observations/satellite/TRMM/2011/*20110[89]*.nc')    

# Go through TRMM files and get specified time and lat/lon range
with open('/nfs/see-fs-01_users/eepdw/python_scripts/filenamelist/trmm_netcdffilelist') as f: 
  for i,line in enumerate(f):
   linestr=line.rstrip()
   #print linestr
   nc = Dataset(linestr)
   

   pcp_in = nc.variables['pcp']
   latitude_in = nc.variables['latitude']
   longitude_in = nc.variables['longitude']
   time_in = nc.variables['time']
   #print longitude_in
   #print pcp_in.shape
   #print time_in
   pcp_dom[i,:,:] = pcp_in[:, la_i_min:la_i_max, lo_i_min:lo_i_max]  
   latitude_dom[i,:] = latitude_in[la_i_min:la_i_max]
   longitude_dom[i,:] = longitude_in[lo_i_min:lo_i_max]
   time_dom[i]=nc.variables['time'].units[12:23]
   time_hour[i]=nc.variables['time'].units[-2:]
   #print time_hour[i]
   
   pdb.set_trace()

   nc.close()

pickle.dump([pcp_dom, longitude_dom, latitude_dom, time_dom, time_hour], open('/nfs/a90/eepdw/Data/Saved_data/TRMM/trmm_emb_time_update_large.p', 'wb'))

        
if '__name__' == '__TRMM_fileread___':
  TRMM_fileread()                                    
