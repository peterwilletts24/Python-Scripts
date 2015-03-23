########################################
# Read TRMM netcdf files.py
#
# Created by: Peter Willetts
# Created on: 22/09/2014
#
########################################
#
#
###################################################

from netCDF4 import Dataset

import numpy as np
import re
import os
import pickle

from datetime import datetime
from time import mktime

import pdb

lon_max = 116 
lon_min = 30.5

lat_max= 40
lat_min=-11.25

date_max=datetime(2011, 9, 8, 0, 0)
date_min=datetime(2011, 8, 18, 0, 0)

date_max_unix=mktime(date_max.timetuple())
date_min_unix=mktime(date_min.timetuple())

#Get number of lines 
no_lines=0
with open('/nfs/see-fs-01_users/eepdw/python_scripts/filenamelist/cmorph_netcdffilelist') as f:
 for i, l in enumerate(f):

   #print l
   linestr=l.rstrip()
   #print linestr
   if linestr.endswith('.nc'):
    nc = Dataset(linestr)   

    date_index = np.where((nc.variables['time'][:]<=date_max_unix) & (nc.variables['time'][:] >= date_min_unix))
 
    if len(date_index[0])>0:
     #print date_index
     date_i_max = np.max(date_index)
     date_i_min = np.min(date_index)
   
     date_amounts=date_i_max-date_i_min 

     no_lines=no_lines+1
     #print no_lines
 
# Get min and max  index positions for latitude and longitude

#nc = Dataset(linestr)

la_index = np.where((nc.variables['lat'][:]<=lat_max) & (nc.variables['lat'][:] >= lat_min))
lo_index = np.where((nc.variables['lon'][:]<=lon_max) & (nc.variables['lon'][:] >= lon_min))

la_i_max = np.max(la_index)
la_i_min = np.min(la_index)
lo_i_max = np.max(lo_index)
lo_i_min = np.min(lo_index)
   
lat_amounts=la_i_max-la_i_min
lon_amounts=lo_i_max-lo_i_min

# Pre-declare arrays

pcp_in = np.zeros((date_amounts,lat_amounts,lon_amounts),dtype=float)
latitude_in = np.zeros((lat_amounts),dtype=float)
longitude_in = np.zeros((lon_amounts),dtype=float)
time_in =  np.zeros((date_amounts),dtype=float)

pcp_dom = np.zeros((no_lines,date_amounts,lat_amounts,lon_amounts),dtype=float)
print pcp_dom.shape
print '(Days, Hours, Latitudes, Longitudes)'
latitude_dom = np.zeros((no_lines,lat_amounts),dtype=float)
longitude_dom =  np.zeros((no_lines,lon_amounts),dtype=float)
time_dom = np.empty((no_lines, date_amounts), dtype=float)


# Go through CMORPH files and get specified time and lat/lon range
l_n=0
with open('/nfs/see-fs-01_users/eepdw/python_scripts/filenamelist/cmorph_netcdffilelist') as f: 
  for i,line in enumerate(f):

   linestr=line.rstrip()
   #print linestr   
   if linestr.endswith('.nc'):
    nc = Dataset(linestr)

    date_index = np.where((nc.variables['time'][:]<=date_max_unix) & (nc.variables['time'][:] >= date_min_unix))
    #print date_index
 
    if len(date_index[0])>0:

     #print 'write'
     #print linestr
     date_i_max = np.max(date_index)
     date_i_min = np.min(date_index)
   
     date_amounts=date_i_max-date_i_min

     pcp_in = nc.variables['cmorph_precip']
     latitude_in = nc.variables['lat']
     longitude_in = nc.variables['lon']
     time_in = nc.variables['time']

     pcp_dom[l_n,:,:,:] = pcp_in[date_i_min:date_i_max,la_i_min:la_i_max, lo_i_min:lo_i_max]  
     latitude_dom[l_n,:] = latitude_in[la_i_min:la_i_max]
     longitude_dom[l_n,:] = longitude_in[lo_i_min:lo_i_max]
     time_dom[l_n,:]=time_in[date_i_min:date_i_max]
  
     pdb.set_trace()
     nc.close()

     l_n=l_n+1

pickle.dump([pcp_dom, longitude_dom, latitude_dom, time_dom], open('/nfs/a90/eepdw/Data/Saved_data/CMORPH/cmorph_emb_time_update_large.p', 'wb'))

        
if '__name__' == '__CMORPH_fileread___':
  CMORPH_fileread()                  
