import cPickle as pickle
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, cm
import numpy as np
import matplotlib.colors as colors

import os

import cartopy.crs as ccrs

sum_dom, latitude_domsingle, longitude_domsingle= pickle.load(open('/nfs/see-fs-01_users/eepdw/Saved_data/TRMM/trmm_emb_pcpsum.p', 'rb'))

# Calculate total at each lat,lon position

#mean_dom = np.mean(pcp_dom, axis=0)

#sum_dom = np.sum(pcp_dom, axis=0)

lon_mid=longitude_domsingle[90]
lat_mid=latitude_domsingle[80]
lons= longitude_domsingle[:]
lats = latitude_domsingle[:]

lons,lats = np.meshgrid(lons, lats)
#date_range_check = 

#lon_0 = -nc.variables['true_lon'].getValue()
#lat_0 = nc.variables['true_lat'].getValue()

# create figure and axes instances
fig = plt.figure(figsize=(10,12))
ax = fig.add_axes([0.05,0.05,0.9,0.85])
# create polar stereographic Basemap instance.
#m = Basemap(projection='ste',lon_0=lon_mid,lat_0=lat_mid,lat_ts=lat_mid,\

m = Basemap(projection='mill',\
            llcrnrlat=-10.,urcrnrlat=30.,\
            llcrnrlon=60.,urcrnrlon=105.,\
            rsphere=6371200.,resolution='h',area_thresh=10000)
# draw coastlines, state and country boundaries, edge of map.
m.drawcoastlines(linewidth=0.5,color='#262626')
#m.drawstates()
m.drawcountries(linewidth=0.5,color='#262626')
# draw parallels.
parallels = np.arange(0.,90,10.)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=16)
# draw meridians
meridians = np.arange(0.,360.,10.)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=16)
#ny = mean_dom.shape[0]; nx = mean_dom.shape[1]
#lons, lats = m.makegrid(longitude_dom[1,:], latitude_dom[1,:]) # get lat/lons of ny by nx evenly space grid.
x, y = m(lons, lats) # compute map proj coordinates.
# draw filled contours.

#clevs = [0,1,2.5,5,7.5,10,15,20,30,40,50,70,100,150,200,250,300,400,500,600,750]
clevs = np.linspace(0, 360,256)
#clevs = np.logspace(0, 8.9657842847,num=50, endpoint=True, base=2.0)
#clevs = [0.10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 350,  400]
cs = m.contourf(x,y,sum_dom,clevs,cmap=cm.s3pcpn_l, extend='both')
# add colorbar.
cbar = m.colorbar(cs,location='bottom',pad="5%", ticks=[0,100,200,300,360])
cbar.set_label('mm')
cbar.ax.tick_params(labelsize=16) 

# add title
plt.title('TRMM rainfall Estimate Total for EMBRACE Period')
for item in([ax.title]):
  item.set_fontsize(20)
#plt.xlabel('Longitude')
#plt.ylabel('Latitude')
#plt.show()

save_path='/nfs/see-fs-01_users/eepdw/Figures/TRMM/'
if not os.path.exists('%s' % save_path): os.makedirs('%s'  % save_path)
plt.savefig('%s/TRMM_total_Embrace_360dpi.png' % save_path, format='png', bbox_inches='tight', dpi=360)
