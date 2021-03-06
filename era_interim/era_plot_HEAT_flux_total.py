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

from matplotlib.colors import from_levels_and_colors

import scipy.interpolate

rc('font', family = 'serif', serif = 'cmr10')
rc('text', usetex=True)

rcParams['text.usetex']=True
rcParams['text.latex.unicode']=True
rcParams['font.family']='serif'
rcParams['font.serif']='cmr10'
rcParams['font.weight']='normal'
rcParams['text.color']='#262626'

plot_type='mean'
plot_diag=''

cmap= cmap=plt.cm.RdBu_r

min_contour = 0
max_contour = 200
tick_gap=20.  

lon_high_plot = 102
lon_low_plot = 64
lat_high_plot= 30.
lat_low_plot=-10

divisor=10  # for lat/lon rounding

latent_mean = np.load('/nfs/a90/eepdw/Data/Saved_data/era_i/Fluxes/era_i_emb_latent_mean.npy')
sensible_mean = np.load('/nfs/a90/eepdw/Data/Saved_data/era_i/Fluxes/era_i_emb_sensible_mean.npy')
#swave_mean = np.load('/nfs/a90/eepdw/Data/Saved_data/era_i/Fluxes/era_i_emb_swave_mean.npy')
#lwave_mean = np.load('/nfs/a90/eepdw/Data/Saved_data/era_i/Fluxes/era_i_emb_lwave_mean.npy')

total_mean = sensible_mean + latent_mean
 
lats = np.load('/nfs/a90/eepdw/Data/Saved_data/era_i/Fluxes/era_i_emb_lats.npy')
lons = np.load('/nfs/a90/eepdw/Data/Saved_data/era_i/Fluxes/era_i_emb_longs.npy')

lon_low= np.min(lons)
lon_high = np.max(lons)
lat_low = np.min(lats)
lat_high = np.max(lats)

lons,lats = np.meshgrid(lons, lats)

m =\
Basemap(llcrnrlon=lon_low_plot,llcrnrlat=lat_low_plot,urcrnrlon=lon_high_plot,urcrnrlat=lat_high_plot,projection='mill', rsphere=6371229)

x, y = m(lons, lats)
fig=plt.figure(figsize=(8,8))
ax = fig.add_axes([0.05,0.05,0.9,0.85])

clevs = np.linspace(min_contour, max_contour,256)
midpoint=0
midp = np.mean(np.c_[clevs[:-1], clevs[1:]], axis=1)

vals = np.interp(midp, [min_contour, midpoint, max_contour], [0, 0.5, 1])
cols = plt.cm.RdBu_r(vals)

clevs_extend = np.linspace(min_contour, max_contour,254)
cmap, norm = from_levels_and_colors(clevs_extend, cols, extend='both')  

# draw coastlines, state and country boundaries, edge of map.
m.drawcoastlines(linewidth=0.5,color='#262626')
#m.drawstates()
m.drawcountries(linewidth=0.5,color='#262626')
# draw parallels.
parallels = np.arange(0.,90,divisor)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10, color='#262626' )
# draw meridians
meridians = np.arange(0.,360., divisor)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10, color='#262626')

cs_col = m.contourf(x,y, total_mean,  clevs, cmap=cmap, extend='both')
cbar = m.colorbar(cs_col,location='bottom',pad="5%")
                               #cbar.ax.tick_params(labelsize=12,  colors='#262626')

ticks= np.arange(int(min_contour),int(max_contour)+tick_gap,tick_gap)
cbar.set_ticks(ticks, update_ticks=True)
cbar.set_ticklabels(([r"${%s}$" % x for x in ticks]))

cbar.set_label('$W m^{-2}$', fontsize=12, color='#262626')  
plt.suptitle('ERA-Interim Reanalysis Mean Total Energy Flux for EMBRACE period', fontsize=16, color='#262626')  

#plt.show()

plt.savefig('/nfs/a90/eepdw/Figures/ERA_Interim/Era_Interim_mean_EMBRACE_period_notitle_total_HEAT_mean.png', format='png', bbox_inches='tight')

plt.suptitle('ERA-Interim Reanalysis Mean Total Heat (Latent and Sensible) Flux for EMBRACE period', fontsize=16, color='#262626')
plt.savefig('/nfs/a90/eepdw/Figures/ERA_Interim/Era_Interim_mean_EMBRACE_period_total_HEAT_mean.png', format='png', bbox_inches='tight')

