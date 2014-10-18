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

import scipy.interpolate

rc('font', family = 'serif', serif = 'cmr10')
rc('text', usetex=True)

rcParams['text.usetex']=True
rcParams['text.latex.unicode']=True
rcParams['font.family']='serif'
rcParams['font.serif']='cmr10'
rcParams['font.weight']='normal'
rcParams['text.color']='#262626'

plot_levels = [925, 850, 700, 500] 
#plot_levels = [925] 
plot_type='mean'
plot_diag='temperature'

min_contour = 0
max_contour = 2
tick_interval=0.2                  

lon_high_plot = 102
lon_low_plot = 64
lat_high_plot= 30.
lat_low_plot=-10



divisor=10  # for lat/lon rounding

geopotential, latitude_domsingle, longitude_domsingle= pickle.load(open('/nfs/a90/eepdw/Data/Saved_data/era_i/era_i_emb_geopotential_mean.p', 'rb'))

variable, latitude_domsinglet, longitude_domsinglet= pickle.load(open('/nfs/a90/eepdw/Data/Saved_data/era_i/era_i_emb_%s_mean.p' % plot_diag, 'rb'))

u_wind, latitude_domsingleu, longitude_domsingleu= pickle.load(open('/nfs/a90/eepdw/Data/Saved_data/era_i/era_i_emb_u_wind_mean.p', 'rb'))
v_wind, latitude_domsinglev, longitude_domsinglev= pickle.load(open('/nfs/a90/eepdw/Data/Saved_data/era_i/era_i_emb_v_wind_mean.p', 'rb'))

pressure_levels =  pickle.load(open('/nfs/a90/eepdw/Data/Saved_data/era_i/era_i_emb_pressure_levels.p', 'rb'))

# convert geopotential to geopotential height

# Calculate total at each lat,lon position

#mean_dom = np.mean(pcp_dom, axis=0)

#sum_dom = np.sum(pcp_dom, axis=0)

lons= longitude_domsingle[:]
lats = latitude_domsingle[:]

lon_low= np.min(lons)
lon_high = np.max(lons)
lat_low = np.min(lats)
lat_high = np.max(lats)

lons,lats = np.meshgrid(lons, lats)

for p in plot_levels:

# Find index of plot level in datasets

    s = np.searchsorted(pressure_levels[::-1], p)

# Get plot grids on pressure level
   # /9.81 to convert geopotential to geopotential height
    plt_h = geopotential[-(s+1),:,:]/9.81
    plt_u_wind = u_wind[-(s+1),:,:]
    plt_v_wind = v_wind[-(s+1),:,:]
 
    plt_v = variable[-(s+1),:,:]
   

 ## Regrid winds onto 2 degree grid for clarity when plotting

    # 2 degree lats lon lists for wind regridding
    lat_wind_1deg = np.arange(lat_low,lat_high, 2)
    lon_wind_1deg = np.arange(lon_low,lon_high, 2)

    lons_wi, lats_wi = np.meshgrid(lon_wind_1deg, lat_wind_1deg)
    fl_la_lo = (lats.flatten(),lons.flatten())

    u = scipy.interpolate.griddata(fl_la_lo, plt_u_wind.flatten(), (lats_wi, lons_wi), method='linear')
    v = scipy.interpolate.griddata(fl_la_lo, plt_v_wind.flatten(), (lats_wi, lons_wi), method='linear')


    m_title = 'Height of %s-hPa level (m)' % (p)

# Set pressure height contour min/max
    if p == 925:
            clev_min = 680.
            clev_max = 810.
    elif p == 850:
            clev_min = 1435.
            clev_max = 1530.
    elif p == 700:
            clev_min = 3090.
            clev_max = 3155.
    elif p == 500:
            clev_min = 5800.
            clev_max = 5890.
    else:
            print 'Contour min/max not set for this pressure level'

# Set potential temperature min/max       
    if p == 925:
            clevpt_min = 298.
            clevpt_max = 310.
    elif p == 850:
            clevpt_min = 302.
            clevpt_max = 312.
    elif p == 700:
            clevpt_min = 312.
            clevpt_max = 320.
    elif p == 500:
            clevpt_min = 325.
            clevpt_max = 332.
    else:
            print 'Potential temperature min/max not set for this pressure level'


  # Set specific humidity min/max       
    if p == 925:
            clevsh_min = 0.012
            clevsh_max = 0.020
    elif p == 850:
            clevsh_min = 0.007
            clevsh_max = 0.017
    elif p == 700:
            clevsh_min = 0.002
            clevsh_max = 0.010
    elif p == 500:
            clevsh_min = 0.001
            clevsh_max = 0.005
    else:
            print 'Specific humidity min/max not set for this pressure level'

    clevvort_min = -5
    clevvort_max = 5
    #clevs_col = np.arange(clev_min, clev_max,256)
    clevs_lin = np.linspace(clev_min, clev_max, num=20)

    m =\
Basemap(llcrnrlon=lon_low_plot,llcrnrlat=lat_low_plot,urcrnrlon=lon_high_plot,urcrnrlat=lat_high_plot,projection='mill', rsphere=6371229)

    x, y = m(lons, lats)
    x_w, y_w = m(lons_wi, lats_wi)
    fig=plt.figure(figsize=(8,8))
    ax = fig.add_axes([0.05,0.05,0.9,0.85])


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
#ny = mean_dom.shape[0]; nx = mean_dom.shape[1]
#lons, lats = m.makegrid(longitude_dom[1,:], latitude_dom[1,:]) # get lat/lons of ny by nx evenly space grid.

# draw geopotential contour lines 

    cs_lin = m.contour(x,y, plt_h, clevs_lin,colors='#262626',linewidths=0.5)
        
    if plot_diag=='temperature':
             clevspt_nums=clevpt_max-clevpt_min+1
             plt_v = np.ma.masked_outside(plt_v, clevpt_max+20,  clevpt_min-20)
             tick_gap=2
             cs_col = m.contourf(x,y, plt_v,  np.linspace(clevpt_min, clevpt_max, clevspt_nums), cmap=plt.cm.jet, extend='both')
             cbar = m.colorbar(cs_col,location='bottom',pad="5%")
                               #cbar.ax.tick_params(labelsize=12,  colors='#262626')
             tick_gap=2
        
             ticks= np.arange(int(clevpt_min),int(clevpt_max)+tick_gap,tick_gap)
             cbar.set_ticks(ticks, update_ticks=True)
             cbar.set_ticklabels(([r"${%s}$" % x for x in ticks]))

             cbar.set_label('Potential Temperature ${\\theta}$(K)', fontsize=12, color='#262626')  
             plt.suptitle('Height, Potential Temperature and Wind Vectors at %s hPa'% (p), fontsize=16, color='#262626')  

    elif plot_diag=='sphum':
             clevssh_nums=clevpt_max-clevpt_min+1
             plt_v = np.ma.masked_outside(plt_v, clevsh_max+20,  clevsh_min-20)

             cs_col = m.contourf(x,y, plt_v,  np.linspace(clevsh_min, clevsh_max, clevssh_nums), cmap=plt.cm.jet_r, extend='both')
             cbar = m.colorbar(cs_col,location='bottom',pad="5%", format = '%.3f') 
             
             tick_gap=0.002
             ticks= np.arange(clevsh_min,clevsh_max+tick_gap,tick_gap)

             cbar.set_ticks(ticks)
             cbar.set_ticklabels((["${%.3f}$" % x for x in ticks]) )

             cbar.set_label('Specific Humidity ${\\phi}$(kg/kg)', fontsize=12, color='#262626')
             plt.suptitle('Height, Specific Humidity and Wind Vectors at %s hPa'% (p),  fontsize=16, color='#262626')

    elif plot_diag=='vort':
        clevvort_min = -5    
        clevvort_max = 5

        cs_col = m.contourf(x,y, plt_v*(10**5), np.linspace(clevvort_min, clevvort_max), cmap=plt.cm.RdBu_r, extend='both')
        cbar = m.colorbar(cs_col,location='bottom',pad="5%", format = '%i') 
        
        tick_gap=1
        
     
        ticks= np.arange(clevvort_min,clevvort_max+tick_gap,tick_gap)

        cbar.set_ticks(ticks)
        cbar.set_ticklabels((["${%i}$" % x for x in ticks]) )

        cbar.set_label('Relative Vorticity (${10^{-5}\ s^{-1}}$)', fontsize=12, color='#262626')
        plt.suptitle('Height, Relative Vorticity and Wind Vectors at %s hPa'% (p),  fontsize=16, color='#262626')
    elif plot_diag=='ptvort':
        clevvort_min = -0.1
        clevvort_max = 0.5

        cs_col = m.contourf(x,y, plt_v*(10**6), np.linspace(clevvort_min, clevvort_max), cmap=plt.cm.RdBu_r, extend='both')
        cbar = m.colorbar(cs_col,location='bottom',pad="5%") 
        #K m**2 kg**-1 s**-1
       
        tick_gap=0.1
        ticks= np.arange(clevvort_min,clevvort_max+tick_gap,tick_gap)

        cbar.set_ticks(ticks)
        cbar.set_ticklabels((["${%.1f}$" % x for x in ticks]) )

        cbar.set_label('Potential Vorticity (${K\ m^{2}\ kg^{-1}\ s^{-1}}$)', fontsize=12, color='#262626')
        plt.suptitle('Height, Potential Vorticity and Wind Vectors at %s hPa'% (p),  fontsize=16, color='#262626')
    elif plot_diag=='div':
        clevvort_min = -1.5
        clevvort_max = 1.5

        cs_col = m.contourf(x,y, plt_v*(10**5), np.linspace(clevvort_min, clevvort_max), cmap=plt.cm.RdBu_r, extend='both')
        cbar = m.colorbar(cs_col,location='bottom',pad="5%") 
      
        tick_gap=0.3
        ticks= np.arange(clevvort_min,clevvort_max+tick_gap,tick_gap)

        cbar.set_ticks(ticks)
        cbar.set_ticklabels((["${%.1f}$" % x for x in ticks]) )

        cbar.set_label('Divergence of Wind (${s^{-1}}$)', fontsize=12, color='#262626')
        plt.suptitle('Height, Divergence and Wind Vectors at %s hPa'% (p),  fontsize=16, color='#262626')

# Scale 150 for diff plots, scale 400 for mean state plots
 # wind = m.quiver(x_w,y_w, u, v, scale=150,color='#262626' )
    wind = m.quiver(x_w,y_w, u, v, scale=400, color='#262626' )
    qk = plt.quiverkey(wind, 0.14, 0.072, 5, '${5\ ms^{-1}}$', labelpos='W', fontproperties={'weight':'heavy', 'size':'14'}, labelcolor='#262626', color='#262626' )
                
    plt.clabel(cs_lin, fontsize=10, fmt="${%i}$", color='#262626')
   # cbar.ax.tick_params(labelsize=10,  color='#262626', ')

    #plt.show()

    plt.savefig('/nfs/a90/eepdw/Figures/ERA_Interim/Era_Interim_mean_EMBRACE_period_%shPa_%s_shorttitle.png' % (p,plot_diag), format='png', bbox_inches='tight')

    #plt.title('TRMM Ra for EMBRACE Period ' , fontsize=16, color='#262626')

    plt.savefig('/nfs/a90/eepdw/Figures/ERA_Interim/Era_Interim_mean_EMBRACE_period_%shPa_%s.png' % (p,plot_diag), format='png', bbox_inches='tight')

    plt.suptitle('', visible=False)

    plt.savefig('/nfs/a90/eepdw/Figures/ERA_Interim/Era_Interim_mean_EMBRACE_period_%shPa_%s_notitle.png' % (p,plot_diag), format='png', bbox_inches='tight')
