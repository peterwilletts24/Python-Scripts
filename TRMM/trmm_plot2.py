import pickle
import
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

pickle.load(open('trmm_emb.p', 'rb'))

# Calculate mean at each lat,lon position

mean_dom = np.mean(pcp_dom, axis=0)

#date_range_check = 

#lon_0 = -nc.variables['true_lon'].getValue()
#lat_0 = nc.variables['true_lat'].getValue()

# create figure and axes instances
fig = plt.figure(figsize=(8,8))
ax = fig.add_axes([0.1,0.1,0.8,0.8])
# create polar stereographic Basemap instance.
m = Basemap(projection='stere',lon_0=lon_0,lat_0=90.,lat_ts=lat_0,\
            llcrnrlat=-10.,urcrnrlat=30.,\
            llcrnrlon=60.,urcrnrlon=105.,\
            rsphere=6371200.,resolution='l',area_thresh=10000)
# draw coastlines, state and country boundaries, edge of map.
m.drawcoastlines()
m.drawstates()
m.drawcountries()
# draw parallels.
parallels = np.arange(0.,90,10.)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
# draw meridians
meridians = np.arange(180.,360.,10.)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
ny = mean_dom.shape[0]; nx = mean_dom.shape[1]
lons, lats = m.makegrid[longitude_dom[1,:], latitude_dom[1,:]] # get lat/lons of ny by nx evenly space grid.
x, y = m(lons, lats) # compute map proj coordinates.
# draw filled contours.
clevs = [0,1,2.5,5,7.5,10,15,20,30,40,50,70,100,150,200,250,300,400,500,600,750]
cs = m.contourf(x,y,data,clevs,cmap=cm.s3pcpn)
# add colorbar.
cbar = m.colorbar(cs,location='bottom',pad="5%")
cbar.set_label('mm')
# add title
plt.title(prcpvar.long_name+' for period ending '+prcpvar.dateofdata)
plt.show()
