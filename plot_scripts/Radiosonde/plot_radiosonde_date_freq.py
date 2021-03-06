#################
# Plot frequency of radiosonde soundings by station vs date/time

##########

import os
import cPickle as pickle 
import datetime

import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from collections import defaultdict
import numpy as np
import matplotlib.cm as cm 
# Load file from radiosonde_read_separate.py


from matplotlib import rc
from matplotlib.font_manager import FontProperties
from matplotlib import rcParams

rc('font', family = 'serif', serif = 'cmr10')
rc('text', usetex=True)

rcParams['text.usetex']=True
rcParams['text.latex.unicode']=True
rcParams['font.family']='serif'
rcParams['font.serif']='cmr10'


date, time, station = pickle.load(open('/nfs/a90/eepdw/Saved_data/date_time_station_rsondes.p', 'rb'))

# Combine date and time lists
dt=[]
x=[]
pdt=[]
dts=[]

# Combine date and time variable into one variable

for i,j in zip(date,time):
 dt.append(i+j)

# Create lists with station name and date and time
for w,e in zip(station,dt):
 d =w,e
 dts.append(d)

# Create sorted list of times at each station i.e at this station, these times
d = defaultdict(list)
for k, v in dts:
    d[k].append(v)

plt_st=[0] * len(d)
plt_tim=[0] * len(d)

for e,p in enumerate(d):
    plt_st[e]=(d.items()[e][0]).title()
    plt_tim[e]=(d.items()[e][1])

s= sorted(zip(plt_st, plt_tim), key=lambda station: station[0], reverse=True)
plt_st,plt_tim=zip(*s)
# Convert to python date time

for l,q in enumerate(plt_tim):
    for n,m in enumerate(q):
  # Some times with 99 value eg 00:99 - replace
     lr=m.replace("99", "00")
    
     p = datetime.datetime.strptime(lr, "%Y%m%d%H:%M")
     plt_tim[l][n]=p

b = np.arange(10)
ys = [i+b+(i*b)**2 for i in range(l+1)]
colors = cm.rainbow(np.linspace(0, 1, len(ys)))

time_min = min(min(plt_tim))
time_max = max(max(plt_tim))

plt.figure(figsize=(18,12))
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%m/%d/%H:%M'))
plt.gca().xaxis.set_major_locator(mdates.DayLocator(interval=2))
plt.gca().xaxis.set_minor_locator(mdates.DayLocator(interval=1))
plt.gca().yaxis.set_ticks([])
#gridlines=plt.gca().get_xgridlines()

plt.gca().xaxis.grid(True)

#for line in gridlines:
 #   line.set_linestyle('-')

for t, x in enumerate(plt_tim):
 y = [1*t] * len(x)

 
 c = colors[t]
 l=plt_st[t].replace('_',' ')
 
 plt.scatter(x, y, label=l, color = c)
 plt.gcf().autofmt_xdate()

 min(min(plt_tim))
 pos = [min(plt_tim[t]), (y[0])]
 plt.text(time_min-datetime.timedelta(hours=8), pos[1], l, size=12, rotation=0, ha="right", va="center")

#plt.legend(loc='upper left', prop={'size':4})


plt.ylim([-0.5,y[0]+0.5])
plt.xlim([time_min-datetime.timedelta(hours=5), time_max+datetime.timedelta(hours=5)])
plt.title('Sounding date and time for each station in model domain')
#plt.show() 
plt.savefig('/nfs/see-fs-01_users/eepdw/Figures/sounding_date_time_each_station_embrace.png',  format='png', bbox_inches='tight')
