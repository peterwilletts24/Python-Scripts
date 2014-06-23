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



min_contour = 0
max_contour = 2
tick_interval=0.2                  

lon_high = 102
lon_low = 64
lat_high= 30.
lat_low=-10.5

divisor=10  # for lat/lon rounding

sum_dom, latitude_domsingle, longitude_domsingle= pickle.load(open('/nfs/a90/eepdw/Data/Saved_data/TRMM/trmm_emb_pcpmean.p', 'rb'))

def draw_screen_poly( lats, lons, m):
    x, y = m( lons, lats )
    xy = zip(x,y)
    poly = Polygon( xy, edgecolor='#262626', facecolor='none', alpha=1, linewidth=2 )
    if (plot_coords[2]=='Southern Indian Ocean'):
       poly = Polygon( xy, edgecolor='red', facecolor='none', alpha=0.4, linewidth=5, label=l+1 ) 

    plt.gca().add_patch(poly)

    legendEntries.append(l+1)
    legendtext.append(plot_coords[2])

def label(lats, lons,  text):
    #y = xy[1] - 0.15 # shift y-value for label so that it's below the artist   
    lons_label = (np.max(lons)+np.min(lons)) / 2
    lats_label = (np.max(lats) + np.min(lats)) / 2
    x, y = m( lons_label, lats_label ) 
    #plt.text(x, y, text, color='#262626', ha="center", va="center", size=32, backgroundcolor='white', alpha=0.4 )
    font0 = FontProperties()
    font0.set_family('sans-serif')
   
    plt.text(x, y, text, color='black', ha="center", va="center", size=64 , fontweight='bold', fontproperties=font0)
    if (plot_coords[2]=='Southern Indian Ocean'):
       plt.text(x, y, text, color='red', ha="center", va="center", size=64, fontproperties=font0)

# Bit above Western Ghats
lats_1 = [20, 28, 28, 20]
lons_1 = [67, 67, 71, 71]
label_1 = 'Bit above Western Ghats'

# Western Ghats
#lats_2 = [8, 21, 21, 8]
#lons_2 = [72, 72, 77, 77]
#label_2 = 'Western Ghats'

# Western Ghats
lats_2 = [8.75, 22., 22., 8.75]
lons_2 = [73.75, 70., 73.75, 77.75]
label_2 = 'Western Ghats'

# Bay of Bengal
lats_3 = [10, 25, 25, 10]
lons_3 = [80, 80, 100, 100]
label_3 = 'Bay of Bengal'

# Southern , western Indian Ocean
lats_4 = [-10, 5, 5, -10]
lons_4 = [64.12, 64.12, 80, 80]
label_4 = 'Southern, western Indian Ocean'

# Southern , western Indian Ocean
lats_5 = [-10, 5, 5, -10]
lons_5 = [80, 80, 101.87, 101.87]
label_5 = 'Southern, eastern Indian Ocean'

# Southern Indian Ocean
lats_6 = [-10, 5, 5, -10]
lons_6 = [64.12, 64.12, 101.87, 101.87]
label_6 = 'Southern Indian Ocean'

# Monsoon Trough
lats_7 = [21., 16., 22., 27]
lons_7 = [73., 83., 87., 75]
label_7 = 'Monsoon Trough'

# Himalayas
lats_8 = [25.8, 26.3, 30., 30., 28.5, 27.8, 27.8, 25.8]
lons_8 = [90., 83., 76.3, 82.7, 86.3, 90., 95., 95.]
label_8 = 'Himalayas'

#Ganga Basin
lats_9 = [22, 27., 30., 26.2, 25.8, 22]
lons_9 = [87, 75, 76.3, 83, 90., 90.]
label_9 = 'Ganga Basin'

lats_poly = lats_1, lats_2, lats_3, lats_4, lats_5, lats_6, lats_7, lats_8, lats_9
lons_poly = lons_1, lons_2, lons_3, lons_4, lons_5, lons_6, lons_7, lons_8, lons_9
labels = label_1, label_2, label_3, label_4, label_5, label_6, label_7, label_8, label_9

NUM_COLOURS = len(labels)
cmap=cm.get_cmap(cm.Set1, NUM_COLOURS*2)

legendEntries=[]
legendtext=[]


# Calculate total at each lat,lon position

#mean_dom = np.mean(pcp_dom, axis=0)

#sum_dom = np.sum(pcp_dom, axis=0)

lon_mid=longitude_domsingle[90]
lat_mid=latitude_domsingle[80]
lons= longitude_domsingle[:]
lats = latitude_domsingle[:]

lons,lats = np.meshgrid(lons, lats)

#lon_0 = -nc.variables['true_lon'].getValue()
#lat_0 = nc.variables['true_lat'].getValue()

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
#ny = mean_dom.shape[0]; nx = mean_dom.shape[1]
#lons, lats = m.makegrid(longitude_dom[1,:], latitude_dom[1,:]) # get lat/lons of ny by nx evenly space grid.
x, y = m(lons, lats) # compute map proj coordinates.
# draw filled contours.

clevs = np.linspace(min_contour, max_contour,256)
ticks = (np.arange(min_contour, max_contour+tick_interval,tick_interval))

cs = m.contourf(x,y,sum_dom, clevs, cmap=cm_base.s3pcpn_l, extend='both')
# add colorbar.
#cbar = m.colorbar(cs,location='bottom',pad="5%")
cbar = m.colorbar(cs,location='bottom',pad="5%")

cbar.set_ticklabels(['%.1f' % i for i in ticks])
cbar.set_label('mm/h')

for l,plot_coords in enumerate(zip(lats_poly,lons_poly, labels)):
    colour = cmap(1.*(l*2)/(NUM_COLOURS*2))
    print labels
    draw_screen_poly( plot_coords[0], plot_coords[1], m)
    label(plot_coords[0], plot_coords[1], l+1)
    
#leg = plt.legend(legendEntries, legendtext, bbox_to_anchor=(1.05, 1), frameon=False, loc=2, borderaxespad=0.)

 # Change the legend label colors to almost black
#texts = leg.texts
#for t in texts:
#    t.set_color('#262626')
# add title
#plt.title('TRMM Rainfall Retrieval Total for EMBRACE Period - 21 days from 21st August 2011', , fontsize=16, color='#262626')
#for item in([ax.title]):
 # item.set_fontsize(8)
#plt.xlabel('Longitude')
#plt.ylabel('Latitude')
#plt.show()
plt.savefig('/nfs/a90/eepdw/Figures/TRMM/TRMM_mean_EMBRACE_period_with_boxes_drawn_on_notitle.png', format='png', bbox_inches='tight')

plt.title('TRMM Rainfall Retrieval Total for EMBRACE Period - 21 days from 21st August 2011' , fontsize=16, color='#262626')

plt.savefig('/nfs/a90/eepdw/Figures/TRMM/TRMM_mean_EMBRACE_period_with_boxes_drawn_on.png', format='png', bbox_inches='tight')
