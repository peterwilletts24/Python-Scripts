#############
# Count number of strings in list that are the same
#
#
#############

from collections import defaultdict
import os
import numpy as np
import cPickle as pickle 
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

data, in_header, np_soundings = pickle.load(open('india_radiosonde_aug_sep2011.p', 'rb'))

d = defaultdict(int)

for q, location in enumerate(in_header[8]):
 d[location] += 1

lc=0
for i2, c2 in enumerate(d):
 lc=lc+1
st_lat_c=[0] * lc
st_lon_c=[0] * lc
st_nam_c=[0] * lc
st_cnt_c=[0] * lc
st_namcnt_c=[0] * lc
for i, c in enumerate(d):
    
    for p, l in enumerate(in_header[8]):
     if c in l:
       st_lat_c[i]=float(in_header[5][p])
       st_lon_c[i]=float(in_header[6][p])
       st_nam_c[i]=(d.items()[i][0]).title()
       st_cnt_c[i]=str(d.items()[i][1])
       st_namcnt_c[i]=st_nam_c[i] + ': ' + st_cnt_c[i]


stt=st_nam_c, st_lat_c, st_lon_c, st_cnt_c
st_name_count=np.array(stt)


 #d.items()[location][0], d.items()[location][1]


#PLOT TIME AND DATE vS FREQ ALL STATIONS

# create figure and axes instances
fig = plt.figure(figsize=(8,8))
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
 plt.text(i, j, s, fontsize=6)

    #ax.annotate( st_nam_c, xy = (x, y), xycoords='data', horizontalalignment='right', verticalalignment='center' )
plt.title('Position and number of soundings in August and September 2011')
plt.show()
