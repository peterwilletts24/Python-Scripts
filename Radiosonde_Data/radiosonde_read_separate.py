########################################
# Read monthly radiosonde files in, separate soundings and save in python
#
# Created by: Peter Willetts
# Created on: 12/11/2013
#
########################################
#
#
###################################################import numpy as np

import glob
import numpy as np
import re
import os
import cPickle as pickle

# Create list of monthly radiosonde files to be used as input

rad_flist = glob.glob ('/nfs/a80/eepdw/Observations/Radiosonde/*/*/2011/*_20110[89].na')
match_header = re.compile(r'2011\d\d\d\d \d\d:\d\d')

latmin=-10.
latmax=30.
lonmin=60.
lonmax=105.
#
# Initialise variables

num_soundings=0

date=[]
time=[]
nafile=[]
no_levels=[]
st_height=[]
st_lat=[]
st_lon=[]
st_wmo_no=[]
station=[]
country=[]
lines=[]


#
# Read file searching for start of each indpendent sounding header

for i in rad_flist:
    print i
    f = open(i,'r')
    for line in f:
     line = line.strip()
     columns = line.split()

# Extract independent header variable

     if match_header.search(line) is not None:
      date.append(columns[0])
      time.append(columns[1])
      nafile.append(i)
      
    # Next line
      in_ht= next(f)
      ind_head = in_ht.split()
      no_levels.append(ind_head[0]) 
      st_height.append(ind_head[1]) 
      st_lat.append(ind_head[2])
      st_lon.append(ind_head[3])
      
    # Next line
      wmono = next(f)
      wmono = wmono.rstrip()
      st_wmo_no.append(wmono)
      
    # Next line
      sta = next(f)
      sta = sta.rstrip()
      station.append(sta)

    # Next line
      co = next(f)
      co=co.rstrip()
      country.append(co)
      
      num_soundings=num_soundings+1
    
      rc=int(ind_head[0])

ind_header=date,time,nafile,no_levels,st_height,st_lat,st_lon,st_wmo_no,station,country
in_header=np.array(ind_header) 
                                                                       

f.close()

# Loop through again extracting sounding data      
dsc=0
ml=map(int, no_levels)
max_no_levels=max(ml) 
no_columns=7
np_soundings=num_soundings
data=np.zeros((no_columns, np_soundings, max_no_levels))

for i in rad_flist:
  
    f = open(i,'r')
    for line in f:

# Initialize numpy array for sounding data

     if match_header.search(line) is not None:
    
 # Skip independent header lines
      in_ht= next(f)
      ind_head = in_ht.split()
      next(f)
      next(f)
      next(f)
 # Get sounding data     
      rc=int(ind_head[0])
      if ind_head[2] > latmin and ind_head[2] < latmax and ind_head[3] > lonmin and ind_head[3] < lonmax:
       dlc=0
       for x in range(rc):
        if line is not None:
         da = next(f) 
         das = da.rstrip()
         data[:,dsc,dlc] = das.split()
        dlc = dlc+1
       dsc = dsc+1
      
      #print i
      #print dlc
f.close()

pickle.dump([data, in_header, np_soundings], open('india_radiosonde_aug_sep2011.p', 'wb'))
pickle.dump([date, time, station], open('date_time_station_rsondes.p', 'wb'))
