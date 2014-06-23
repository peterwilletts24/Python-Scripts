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

lon_max = 116 
lon_min = 34

lat_max= 40.
lat_min=-11.25

#dys = 21
#hrs = [0,6,12,18]


#Get number of lines 
no_lines=0
with open('/nfs/see-fs-01_users/eepdw/python_scripts/filenamelist/erai_embrace_netcdffilelist') as f:
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

print nc

pressure_level_amounts=nc.variables['p'].shape[0]


vort_dom = np.zeros((no_lines, pressure_level_amounts, lat_amounts,lon_amounts),dtype=float)
pt_vort_dom = np.zeros((no_lines, pressure_level_amounts, lat_amounts,lon_amounts),dtype=float)
div_dom = np.zeros((no_lines, pressure_level_amounts, lat_amounts,lon_amounts),dtype=float)

latitude_dom = np.zeros((no_lines,lat_amounts),dtype=float)
longitude_dom =  np.zeros((no_lines,lon_amounts),dtype=float)
time_dom = np.empty((no_lines), dtype=(str,12))
time_hour = np.empty((no_lines), dtype=(str,2))

 #fs=glob.iglob('/nfs/a80/earceb/model_data/observations/satellite/TRMM/2011/*20110[89]*.nc')    

# Go through era-interim files and get specified time and lat/lon range
with open('/nfs/see-fs-01_users/eepdw/python_scripts/filenamelist/erai_embrace_netcdffilelist') as f: 
  for i,line in enumerate(f):
   linestr=line.rstrip()
   #print linestr
   nc = Dataset(linestr)

   vort_in = nc.variables['VO']
   div_in = nc.variables['D']
   pt_vort_in = nc.variables['PV']
  

   #print sphum_in[0].shape
   #print  sphum_dom[i,:,:,:].shape

   latitude_in = nc.variables['latitude']
   longitude_in = nc.variables['longitude']
   time_in = nc.variables['t']
   #print longitude_in
   #print pcp_in.shape
   #print time_in
   vort_dom[i,:,:,:] = vort_in[:, :, la_i_min:la_i_max, lo_i_min:lo_i_max]  
   div_dom[i,:,:,:] = div_in[:, :, la_i_min:la_i_max, lo_i_min:lo_i_max]  
   pt_vort_dom[i,:,:,:] = pt_vort_in[:, :, la_i_min:la_i_max, lo_i_min:lo_i_max] 
   

   latitude_dom[i,:] = latitude_in[la_i_min:la_i_max]
   longitude_dom[i,:] = longitude_in[lo_i_min:lo_i_max]
   time_dom[i]=nc.variables['t'].units[11:21]
   time_hour[i]=nc.variables['t'].units[-8:-6]
   #print time_hour[i]
   nc.close()

#variables_list=['sphum_dom, 'geopotential', 'u_wind', 'v_wind', 'cloud_cover', 'temperature']

#for i in variables_list:

#save_variable=

pickle.dump([vort_dom, longitude_dom, latitude_dom, time_dom, time_hour], open('/nfs/a90/eepdw/Data/Saved_data/era_i/era_i_emb_time_update_large_vort.p', 'wb'))
pickle.dump([div_dom, longitude_dom, latitude_dom, time_dom, time_hour], open('/nfs/a90/eepdw/Data/Saved_data/era_i/era_i_emb_time_update_large_div.p', 'wb'))
pickle.dump([pt_vort_dom, longitude_dom, latitude_dom, time_dom, time_hour], open('/nfs/a90/eepdw/Data/Saved_data/era_i/era_i_emb_time_update_large_ptvort.p', 'wb'))
        
if '__name__' == '__netcdf_fileread_vort_div___':
  TRMM_fileread()                                    


#for i in (nc.variables):
 #    print nc.variables[str(i)]
