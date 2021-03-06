########################################
# Read marine observations .txt file and list of variables
#
# Created by: Peter Willetts
# Created on: 26/11/2013
#
########################################
#
#
###################################################

import glob
import numpy as np
import re
import os
import cPickle as pickle
from collections import defaultdict
import csv
import datetime
from operator import itemgetter
import operator
import calendar


# Set domain boundaries

latmin=-10.
latmax=30.
lonmin=60.
lonmax=105.

# Set date boundaries

date_min = datetime.datetime(2011, 8, 1)
date_max = datetime.datetime(2011, 9, 30) 

# Initialize array

full_data=[]

i='/nfs/a80/eepdw/Observations/Marine observations MIDAS/MO_Column_Headers.txt'
f = open(i,'r')
for line in (f):
     line = line.strip()
     line = line.translate(None, ',')
     column_headers = line.split()
     
#for i in column_headers:
#      vars()[i]
f.close()


ifile = open('/nfs/a80/eepdw/Observations/Marine observations MIDAS/midas_marine-obs-lon-band-h_201101-201112.txt', "rb")
reader = csv.reader(ifile)

for row in reader:
 if float(row[1]) > latmin and float(row[1]) < latmax and float(row[2]) > lonmin and float(row[2]) < lonmax:
   date_python = datetime.datetime.strptime(row[0], "%Y-%m-%d %H:%M")
   if date_python >= date_min and date_python <= date_max:
       ap = date_python, row
       full_data.append(ap)

 # Count for types of id
d = defaultdict(int)

for location in [row[1][4] for row in full_data]:
    d[location] += 1

 # Count for different I.D's e.g individual ship I.D
id_count = defaultdict(int)

for location in [row[1][3] for row in full_data]:
    id_count[location] += 1

# Sort into lists grouped by station_id

id_data_sort = defaultdict(list)
for k, v in zip((row[1][3] for row in full_data), (row for row in full_data)):
    id_data_sort[k].append(v)


# Sort into list grouped by datetime

date_data_sort = defaultdict(list)

for k, v in zip((row[1][0] for row in full_data), (row for row in full_data)):
    date_data_sort[k].append(v)

plt_lat=[[] for i in range (len(date_data_sort))]
plt_lon=[[] for i in range (len(date_data_sort))]
plt_tim=[0] * len(date_data_sort)
plt_nam=[[] for i in range (len(date_data_sort))]

plt_com = list(plt_lat), list(plt_lon), list(plt_tim), list(plt_nam)



for d, date in enumerate(date_data_sort.items()):
# try:
     plt_tim[d] =  datetime.datetime.strptime(date[0], "%Y-%m-%d %H:%M")
     for ind in date[1]: 
         plt_lat[d].append(float(ind[1][1]))
         plt_lon[d].append(float(ind[1][2]))
         plt_nam[d].append(ind[1][3])

plt_com = plt_lat, plt_lon, plt_tim, plt_nam

#Transpose so it's the 'right' way round 

plt_com = zip(*plt_com)

plt_com = sorted(plt_com, key = lambda entry: entry[2])    

#Transpose so it's the 'right' way round, again !!! List of list indexing baffling 

plt_com = zip(*plt_com)

plt_com_np=np.array(plt_com)

plt_lat = list(plt_com_np[0])  
plt_lon = list(plt_com_np[1])
plt_tim = list(plt_com_np[2])
plt_nam = list(plt_com_np[3])

pickle.dump([plt_nam, plt_lat, plt_lon, plt_tim], open('marine_obs_times_and_positions.p', 'wb'))   
#pickle.dump([date_data_sort], open('marine_obs_times_and_positions.p', 'wb'))          

 #except IndexError:
 # pass
# print 'Index error'
 #continue
