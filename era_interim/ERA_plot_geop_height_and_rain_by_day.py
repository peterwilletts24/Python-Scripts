import matplotlib

#matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
from matplotlib import rc
from matplotlib.font_manager import FontProperties
from matplotlib import rcParams
from matplotlib import cm

from mpl_toolkits.basemap import Basemap
from mpl_toolkits.basemap import cm

import cPickle as pickle
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as colors
from matplotlib.patches import Polygon

import scipy.interpolate

import datetime

import pdb

import gc

rc('font', family = 'serif', serif = 'cmr10')
rc('text', usetex=True)

rcParams['text.usetex']=True
rcParams['text.latex.unicode']=True
rcParams['font.family']='serif'
rcParams['font.serif']='cmr10'
rcParams['font.weight']='normal'
rcParams['text.color']='#262626'

#plot_levels = [925, 850, 700, 500] 
plot_levels = [925] 
plot_type='mean'
plot_diag='precip'

   
clevpr_min = 0.
clevpr_max = 3.

lon_high_plot = 102
lon_low_plot = 64
lat_high_plot= 30.
lat_low_plot=-10



divisor=10  # for lat/lon rounding

geop_file = np.load('/nfs/a90/eepdw/Data/Era_Interim/Era_Interim_Daily_Geopotential_Height_EMBRACE_Period.npz')
precip_file = np.load('/nfs/a90/eepdw/Data/Era_Interim/Era_Interim_Daily_Total_Precip_EMBRACE_Period.npz')
# convert geopotential to geopotential height

# Calculate total at each lat,lon position

#mean_dom = np.mean(pcp_dom, axis=0)

#sum_dom = np.sum(pcp_dom, axis=0)

#lons= longitude_domsingle[:]
#lats = latitude_domsingle[:]

#lon_low= np.min(lons)
#lon_high = np.max(lons)
#lat_low = np.min(lats)
#lat_high = np.max(lats)

#lons,lats = np.meshgrid(lons, lats)

geopotential = geop_file['data']
variable = precip_file['data']

lons_geop = geop_file['longitudes']
lats_geop = geop_file['latitudes']

lons_var = precip_file['longitudes']
lats_var = precip_file['latitudes']

pressure_levels = geop_file['pressures']

geop_dates=[datetime.datetime.strptime(g, '%Y-%m-%d') for g in geop_file['time_coords']]
var_dates=[datetime.datetime.strptime(g, '%Y-%m-%d') for g in precip_file['time_coords']]

#pdb.set_trace()

for p in plot_levels:

# Find index of plot level in datasets

    s = np.searchsorted(pressure_levels[::-1], p)

    for date in geop_dates:

        dg = np.searchsorted(geop_dates, date)
        dv = np.searchsorted(var_dates, date)
        # Get plot grids on pressure level
        # /9.81 to convert geopotential to geopotential height
        

        plt_h = geopotential[dg, -(s+1),:,:]/9.81
        
        #plt_v = variable[-(s+1),:,:]
        plt_v = variable[dv,:,:]


        # Set pressure height contour min/max
        if p == 925:
            clev_min = 660.
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


        m_title = 'Height of %s-hPa level (m)' % (p)
        # Set precip min/max
        
        clevs_lin = np.arange(clev_min, clev_max, 5)

        m =\
            Basemap(llcrnrlon=lon_low_plot,llcrnrlat=lat_low_plot,urcrnrlon=lon_high_plot,urcrnrlat=lat_high_plot,projection='mill', rsphere=6371229)

        #pdb.set_trace()

        x, y = m(lons_geop, lats_geop)
        x_v, y_v = m(lons_var, lats_var)
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
        
        if plot_diag=='precip':
            clevspt_nums=64
            plt_v = np.ma.masked_outside(plt_v, clevpr_max+0.5,  clevpr_min-0.5)
            tick_gap=0.2
            cs_col = m.contourf(x_v,y_v, plt_v*1000/6,  np.linspace(clevpr_min, clevpr_max, clevspt_nums), cmap=cm.s3pcpn_l, extend='both')
            #cbar = m.colorbar(cs_col,location='bottom',pad="5%")
            #cbar.ax.tick_params(labelsize=12,  colors='#262626')
            tick_gap=0.5
        
            #ticks= np.arange(int(clevpr_min),int(clevpr_max)+tick_gap,tick_gap)
            #cbar.set_ticks(ticks, update_ticks=True)
            #cbar.set_ticklabels(([r"${%s}$" % x for x in ticks]))
            
            #cbar.set_label('Precipitation mm h$^{-1}$', fontsize=12, color='#262626')  
            #plt.suptitle('Height of %s hPa, and Total Precipitation' % (p), fontsize=16, color='#262626')  

        elif plot_diag=='temperature':
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


                
        plt.clabel(cs_lin, fontsize=10, fmt="${%i}$", color='#262626')
        # cbar.ax.tick_params(labelsize=10,  color='#262626', ')
        plt.title('%s UTC' % (datetime.datetime.strftime(date, '%d%b')))

        #plt.show()

        plt.savefig('/nfs/a90/eepdw/Figures/ERA_Interim/Era_Interim_mean_EMBRACE_period_%shPa_height_and_rain_by_day_%s_%s_shorttitle.png' % (p,plot_diag, datetime.datetime.strftime(date, '%d%b')), format='png', bbox_inches='tight')

        plt.clf()
        plt.close()
        
        gc.collect()

        #plt.title('TRMM Ra for EMBRACE Period ' , fontsize=16, color='#262626')

                #plt.savefig('/nfs/a90/eepdw/Figures/ERA_Interim/Era_Interim_mean_EMBRACE_period_%shPa_%s.png' % (p,plot_diag), format='png', bbox_inches='tight')

                #plt.suptitle('', visible=False)

                #plt.savefig('/nfs/a90/eepdw/Figures/ERA_Interim/Era_Interim_mean_EMBRACE_period_%shPa_%s_notitle.png' % (p,plot_diag), format='png', bbox_inches='tight')
