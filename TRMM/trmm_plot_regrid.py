import matplotlib

matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
from matplotlib import rc
from matplotlib.font_manager import FontProperties
from matplotlib import rcParams
from matplotlib import cm

from mpl_toolkits.basemap import Basemap
from mpl_toolkits.basemap import cm as cm_base

import cPickle as pickle
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as colors
from matplotlib.patches import Polygon

rc('font', family = 'serif', serif = 'cmr10')
rc('text', usetex=True)

rcParams['text.usetex']=True
rcParams['text.latex.unicode']=True
rcParams['font.family']='serif'
rcParams['font.serif']='cmr10'



min_contour = -1
max_contour = 1
tick_interval=0.2

numbers_after_point=len(str(tick_interval)[str(tick_interval).find('.')+1])

lon_high = 102
lon_low = 64
lat_high= 30.
lat_low=-10.5

divisor=10  # for lat/lon rounding

experiment_ids = ['djznw', 'djznq', 'djzny', 'djzns', 'dkmbq', 'dklyu', 'dklwu', 'dklzq', 'dkbhu', 'djznu', 'dkhgu', 'dkjxq' ]

sum_dom_trmm, latitude_domsingle, longitude_domsingle= pickle.load(open('/nfs/a90/eepdw/Data/Saved_data/TRMM/trmm_emb_pcpmean.p', 'rb'))

lons= longitude_domsingle[:]
lats = latitude_domsingle[:]

lons,lats = np.meshgrid(lons, lats)

print lons
print lats
for experiment_id in experiment_ids:

  #expmin1 = experiment_id[:-1]

  sum_dom=np.load("/nfs/a90/eepdw/Data/Rain_TRMM_regrid/rain_mean_regrid_onto_trmm_%s.npy" %  experiment_id)

  print sum_dom.shape
# create figure and axes instances
  fig = plt.figure(figsize=(8,8))
  ax = fig.add_axes([0.1,0.1,0.8,0.8])

  m = Basemap(projection='mill',\
            llcrnrlat=lat_low,urcrnrlat=lat_high,\
            llcrnrlon=lon_low,urcrnrlon=lon_high,\
           rsphere=6371229.,resolution='h',area_thresh=10000)


# draw coastlines, state and country boundaries, edge of map.
  m.drawcoastlines(linewidth=0.5,color='#262626')
#m.drawstates()
  m.drawcountries(linewidth=0.5,color='#262626')
# draw parallels.
  parallels = np.arange(0.,90,divisor)
  m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
# draw meridians
  meridians = np.arange(0.,360., divisor)
  m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)


  x, y = m(lons, lats) # compute map proj coordinates.
# draw filled contours.

  clevs = np.arange(min_contour, max_contour+tick_interval,tick_interval)
  ticks = (np.arange(min_contour, max_contour+tick_interval,tick_interval))

  cs = m.contourf(x,y,((sum_dom*3600)-sum_dom_trmm), clevs, cmap=cm.RdBu, extend='both')
# add colorbar.
#cbar = m.colorbar(cs,location='bottom',pad="5%")
  cbar = m.colorbar(cs,location='bottom',pad="5%")
  
  cbar.set_ticks(ticks)
  if numbers_after_point==1:
    cbar.set_ticklabels(["${%.1f}$" % i for i in ticks])
  if numbers_after_point==2:
    cbar.set_ticklabels(["${%.2f}$" % i for i in ticks])
  if numbers_after_point==0:
    cbar.set_ticklabels(["${%d}$" % i for i in ticks])
  cbar.set_label('mm/h')

  plt.savefig('/nfs/a90/eepdw/Figures/TRMM/Rain_regridded_onto_TRMM_and_differenced_%s_notitle.png' % experiment_id, format='png', bbox_inches='tight')

  plt.title('EMBRACE mean rain differenced from TRMM Rainfall Retrieval mean for EMBRACE Period - 21 days from 21st August 2011 (regridded to TRMM first)' , fontsize=16, color='#262626')

  plt.savefig('/nfs/a90/eepdw/Figures/TRMM/Rain_regridded_onto_TRMM_and_differenced_%s.png' % experiment_id, format='png', bbox_inches='tight')

  #plt.show()
