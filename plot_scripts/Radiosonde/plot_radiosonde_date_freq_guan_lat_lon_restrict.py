#################
# Plot frequency of radiosonde soundings by station vs date/time

##########

import os
import cPickle as pickle 
import datetime

import matplotlib
#matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from collections import defaultdict
import numpy as np
import matplotlib.cm as cm 
# Load file from radiosonde_read_separate.py

import re

from matplotlib import rc
from matplotlib.font_manager import FontProperties
from matplotlib import rcParams

from itertools import chain

rc('font', family = 'serif', serif = 'cmr10')
rc('text', usetex=True)

rcParams['text.usetex']=True
rcParams['text.latex.unicode']=True
rcParams['font.family']='serif'
rcParams['font.serif']='cmr10'

station_data=np.load('/nfs/a90/eepdw/Data/Observations/Radiosonde_downloaded_from_NOAA_GUAN/Embrace_Period_India_Station_and_sounding_Info_measured_derived.npy')
#station_data=np.load('/nfs/a90/eepdw/Data/Observations/Radiosonde_downloaded_from_NOAA_GUAN/Embrace_Period_India_Station_and_sounding_Info_derived_parameters.npy')
date_min=datetime.datetime(2011,5,1,0,0,0)
date_max=datetime.datetime(2011,9,30,0,0,0)

lat_min=15
lat_max=30

lon_min = 70
lon_max = 90

match_bad = re.compile(r'9999')
# Combine date and time lists
dt=[]
x=[]
pdt=[]
dts=[]


plt_st=[0.] * len(station_data)
#plt_tim=np.zeros((len(station_data)))
plt_tim=[]


s= sorted(station_data, key=lambda station: station[0])
#plt_st,plt_tim=zip(*s)
# Convert to python date time

for l,q in enumerate(s):
 print q
 print q[2]
 if q[1] > lat_min and float(q[1]) < lat_max and float(q[2]) > lon_min and float(q[2]) < lon_max:
    times=[]
    plt_st[l]=q[0].lower().title()
    for n,m in enumerate(q[3]):
        if match_bad.search(m[0]) is None:
         if match_bad.search(m[1]) is None:
             p = datetime.datetime.strptime('%s%s' % (m[0][5:16], m[1]), "%Y%m%d%H%M")
   
             #print p
             if p>=date_min and p<=date_max:
                 times.append(p)
         #if match_bad.search(m[1]) is not None: 
            #p = datetime.datetime.strptime('%s%s' % (m[0][5:16], m[2]), "%Y%m%d%H")
            #print p
            #if p>=date_min and p<=date_max:
             #   times.append(p)
    plt_tim.append(times)     
    #print plt_tim
#print len(plt_tim)
#print min(list(chain(*plt_tim)))

    b = np.arange(10)
    ys = [i+b+(i*b)**2 for i in range(l+1)]
    colors = cm.rainbow(np.linspace(0, 1, len(ys)))

time_min = min(list(chain(*plt_tim)))
time_max = max(list(chain(*plt_tim)))

plt.figure(figsize=(18,12))
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%m/%d/%H:%M'))
plt.gca().xaxis.set_major_locator(mdates.DayLocator(interval=2))
plt.gca().xaxis.set_minor_locator(mdates.DayLocator(interval=1))
plt.gca().yaxis.set_ticks([])
#gridlines=plt.gca().get_xgridlines()

plt.gca().xaxis.grid(True)

#for line in gridlines:
 #   line.set_linestyle('-')
pos_count=0
for t, x in enumerate(plt_tim):
 if x!=[]:
     y = [1*pos_count] * len(x)
     pos_count+= 1
 
     c = colors[t]
     l=plt_st[t].replace('_',' ')
 
     plt.scatter(x, y, label=l, color = c)
     plt.gcf().autofmt_xdate()

     #print x
     pos = [min(x), (y[0])]
     plt.text(time_min-datetime.timedelta(hours=8), pos[1], l, size=12, rotation=0, ha="right", va="center")

#plt.legend(loc='upper left', prop={'size':4})


plt.ylim([-0.5,y[0]+0.5])
plt.xlim([time_min-datetime.timedelta(hours=5), time_max+datetime.timedelta(hours=5)])
plt.title('Sounding date and time for each station in model domain (Measured)')
plt.show() 
#plt.savefig('/nfs/a90/eepdw/Figures/Observation_Plots/sounding_date_time_each_station_embrace_guan_measured_9999_times_not_included.png',  format='png', bbox_inches='tight')
