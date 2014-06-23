#############
# Count number of strings in list that are the same
#
#
#############

from collections import defaultdict
import os
import numpy as np
import cPickle as pickle 

import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!

import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

station_data=np.load('/nfs/a90/eepdw/Data/Observations/Radiosonde_downloaded_from_NOAA_GUAN/Embrace_Period_India_Station_and_sounding_Info_measured.npy')

#match_bad = re.compile(r'9999')

lc=len(station_data)
st_lat_c=[0] * lc
st_lon_c=[0] * lc
st_nam_c=[0] * lc
st_cnt_c=[0] * lc
st_namcnt_c=[0] * lc

for i, c in enumerate(station_data):
    
       st_lat_c[i]=c[1]
       st_lon_c[i]=c[2]
       st_nam_c[i]=c[0].lower().title()
       #st_cnt_c[i]=str()
       #st_namcnt_c[i]=st_nam_c[i] + ': ' + st_cnt_c[i]
       st_namcnt_c[i]=st_nam_c[i]
#PLOT TIME AND DATE vS FREQ ALL STATIONS

# create figure and axes instances
fig = plt.figure(figsize=(16,16))
ax = fig.add_axes([0.1,0.1,0.8,0.8])
# create polar stereographic Basemap instance.
#m = Basemap(projection='ste',lon_0=lon_mid,lat_0=lat_mid,lat_ts=lat_mid,\

m = Basemap(projection='cyl',\
            llcrnrlat=-10.,urcrnrlat=30.,\
            llcrnrlon=60.,urcrnrlon=105.,\
            rsphere=6371200.,resolution='l',area_thresh=10000)
# draw coastlines, state and country boundaries, edge of map.
m.drawcoastlines()
# draw parallels.
parallels = np.arange(0.,90,10.)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
# draw meridians
meridians = np.arange(0.,360.,10.)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)


x, y = m(st_lon_c,st_lat_c)
m.scatter(x,y,3,marker='o',color='red')


for i,j,s in zip(x, y, st_namcnt_c):
 plt.text(i, j, s, fontsize=10)

plt.title('Position and number of soundings in August and September 2011')
plt.savefig('/nfs/a90/eepdw/Figures/Observation_Plots/sounding_station_map_igra_guan.png',  format='png', bbox_inches='tight')
#plt.show()
