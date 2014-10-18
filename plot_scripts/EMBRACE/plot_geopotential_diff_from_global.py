"""

Load mean geopotential heights and plot in colour

"""
import os, sys
import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt

from matplotlib import rc

import matplotlib.cm as mpl_cm
from mpl_toolkits.basemap import Basemap
import iris
import iris.analysis.cartography
import iris.unit as unit

import datetime

import numpy as np
import imp
import h5py
import cartopy.crs as ccrs

import scipy.interpolate

from textwrap import wrap


model_name_convert_title = imp.load_source('util', '/nfs/see-fs-01_users/eepdw/python_scripts/Modules/model_name_convert_title.py')


def main():

# Plot diagnostics, model and pressure levels etc. to plot on for looping through

 plot_type='mean_masked'
 plot_type_h5py_var='mean'
 plot_diags=['temp', 'sp_hum']
 plot_levels = [925, 850, 700, 500] 
 #plot_levels = [925]
 experiment_ids = ['djzny', 'djznq', 'dklyu', 'dkmbq', 'dklwu', 'dklzq', 'dkjxq', 'dklyu', 'dkmbq', 'dklwu', 'dklzq']
 #experiment_id = 'djzny'

 p_levels = [1000, 950, 925, 850, 700, 500, 400, 300, 250, 200, 150, 100, 70, 50, 30, 20, 10]

 dtmindt = datetime.datetime(2011,8,19,0,0,0)
 dtmaxdt = datetime.datetime(2011,9,7,23,0,0)
 dtmin = unit.date2num(dtmindt, 'hours since 1970-01-01 00:00:00', unit.CALENDAR_STANDARD)
 dtmax = unit.date2num(dtmaxdt, 'hours since 1970-01-01 00:00:00', unit.CALENDAR_STANDARD)
 time_constraint = iris.Constraint(time= lambda t: dtmin <= t.point <= dtmax)

###### Unrotate global model data ##############################

######### Regrid to global, and difference  #######
############################################################################

## Load global wind and orography
 
 fw_global = '/nfs/a90/eepdw/Data/EMBRACE/Mean_State/pp_files/djzn/djznw/30201_mean.pp'
 fo_global = '/nfs/a90/eepdw/Data/EMBRACE/Mean_State/pp_files/djzn/djznw/33.pp'
    
 u_global,v_global = iris.load(fw_global)

 print u_global
 print v_global
 oro_global = iris.load_cube(fo_global)
    
# Unrotate global coordinates

 cs_glob = u_global.coord_system('CoordSystem')
 cs_glob_v = v_global.coord_system('CoordSystem')

 cs_glob_oro = oro_global.coord_system('CoordSystem')

 lat_g = u_global.coord('grid_latitude').points
 lon_g = u_global.coord('grid_longitude').points

 lat_g_oro = oro_global.coord('grid_latitude').points
 lon_g_oro = oro_global.coord('grid_longitude').points
    
 if cs_glob!=cs_glob_v:
    print 'Global model u and v winds have different poles of rotation'

# Unrotate global winds

 if isinstance(cs_glob, iris.coord_systems.RotatedGeogCS):
        print ' Global Model - Winds - djznw - Unrotate pole %s' % cs_glob
        lons_g, lats_g = np.meshgrid(lon_g, lat_g)
        lons_g,lats_g = iris.analysis.cartography.unrotate_pole(lons_g,lats_g, cs_glob.grid_north_pole_longitude, cs_glob.grid_north_pole_latitude)
        
        lon_g=lons_g[0]
        lat_g=lats_g[:,0]

 for i, coord in enumerate (u_global.coords()):
            if coord.standard_name=='grid_latitude':
                lat_dim_coord_uglobal = i
            if coord.standard_name=='grid_longitude':
                lon_dim_coord_uglobal = i

 csur_glob=cs_glob.ellipsoid
 u_global.remove_coord('grid_latitude')
 u_global.remove_coord('grid_longitude')
 u_global.add_dim_coord(iris.coords.DimCoord(points=lat_g, standard_name='grid_latitude', units='degrees', coord_system=csur_glob), lat_dim_coord_uglobal)
 u_global.add_dim_coord(iris.coords.DimCoord(points=lon_g, standard_name='grid_longitude', units='degrees', coord_system=csur_glob), lon_dim_coord_uglobal)

#print u_global
    
 v_global.remove_coord('grid_latitude')
 v_global.remove_coord('grid_longitude')
 v_global.add_dim_coord(iris.coords.DimCoord(points=lat_g, standard_name='grid_latitude', units='degrees', coord_system=csur_glob),  lat_dim_coord_uglobal)
 v_global.add_dim_coord(iris.coords.DimCoord(points=lon_g, standard_name='grid_longitude', units='degrees', coord_system=csur_glob),  lon_dim_coord_uglobal)
    
 #print v_global

# Unrotate global model

 if isinstance(cs_glob_oro, iris.coord_systems.RotatedGeogCS):
        print ' Global Model - Orography - djznw - Unrotate pole %s - Winds and other diagnostics may have different number of grid points' % cs_glob_oro
        lons_go, lats_go = np.meshgrid(lon_g_oro, lat_g_oro)
        lons_go,lats_go = iris.analysis.cartography.unrotate_pole(lons_go,lats_go, cs_glob_oro.grid_north_pole_longitude, cs_glob_oro.grid_north_pole_latitude)
        
        lon_g_oro=lons_go[0]
        lat_g_oro=lats_go[:,0]

 for i, coord in enumerate (oro_global.coords()):
            if coord.standard_name=='grid_latitude':
                lat_dim_coord_og = i
            if coord.standard_name=='grid_longitude':
                lon_dim_coord_og = i 

 csur_glob_oro=cs_glob_oro.ellipsoid
 oro_global.remove_coord('grid_latitude')
 oro_global.remove_coord('grid_longitude')
 oro_global.add_dim_coord(iris.coords.DimCoord(points=lat_g_oro, standard_name='grid_latitude', units='degrees', coord_system=csur_glob_oro), lat_dim_coord_og)
 oro_global.add_dim_coord(iris.coords.DimCoord(points=lon_g_oro, standard_name='grid_longitude', units='degrees', coord_system=csur_glob_oro), lon_dim_coord_og)

   ###############################################################################
####################  Load global heights and temp/sp_hum #####################

 
 f_glob_h = '/nfs/a90/eepdw/Data/EMBRACE/Pressure_level_means/408_pressure_levels_interp_pressure_djznw_%s' % (plot_type)
 
######################################################################################
 with h5py.File(f_glob_h, 'r') as i:
     mh = i['%s' % plot_type_h5py_var]
     mean_heights_global = mh[. . .]

######################################################################################
## Loop through experiment id's ######################################################
 
    
 for  pl in plot_diags:
  plot_diag=pl

  f_glob_d = '/nfs/a90/eepdw/Data/EMBRACE/Pressure_level_means/%s_pressure_levels_interp_djznw_%s' % (plot_diag, plot_type)
  
  with h5py.File(f_glob_d, 'r') as i:
   mg = i['%s' % plot_type_h5py_var]
   mean_var_global = mg[. . .]

  for experiment_id in experiment_ids:
    expmin1 = experiment_id[:-1]

 

    ###############################################################################
####################  Load model heights and temp/sp_hum #####################


    fname_h = '/nfs/a90/eepdw/Data/EMBRACE/Pressure_level_means/408_pressure_levels_interp_pressure_%s_%s' % (experiment_id, plot_type)
    fname_d = '/nfs/a90/eepdw/Data/EMBRACE/Pressure_level_means/%s_pressure_levels_interp_%s_%s' % (plot_diag, experiment_id, plot_type)
    # print fname_h
    # print fname_d
#  Height data file
    with h5py.File(fname_h, 'r') as i:
        mh = i['%s' % plot_type_h5py_var]
        mean_heights = mh[. . .]
    # print mean_heights.shape
    with h5py.File(fname_d, 'r') as i:
        mh = i['%s' % plot_type_h5py_var]
        mean_var = mh[. . .]
    # print mean_var.shape

    f_oro =  '/nfs/a90/eepdw/Data/EMBRACE/Mean_State/pp_files/%s/%s/33.pp' % (expmin1, experiment_id)
    oro = iris.load_cube(f_oro)

    #print oro
    
    for i, coord in enumerate (oro.coords()):
        if coord.standard_name=='grid_latitude':
            lat_dim_coord_oro = i
        if coord.standard_name=='grid_longitude':
            lon_dim_coord_oro = i

    fu = '/nfs/a90/eepdw/Data/EMBRACE/Mean_State/pp_files/%s/%s/30201_mean.pp' % (expmin1, experiment_id)

    time_list = u_global.coord('time').points
    glob_tc = iris.Constraint(time=time_list)
    
    u_wind,v_wind = iris.load(fu, glob_tc)
    
# Wind may have different number of grid points so need to do unrotate again for wind grid points

    lat_w = u_wind.coord('grid_latitude').points
    lon_w = u_wind.coord('grid_longitude').points
    p_levs = u_wind.coord('pressure').points

    lat = oro.coord('grid_latitude').points
    lon = oro.coord('grid_longitude').points
    

    cs_w = u_wind.coord_system('CoordSystem')
    cs = oro.coord_system('CoordSystem')

    if isinstance(cs_w, iris.coord_systems.RotatedGeogCS):
        print ' Wind - %s - Unrotate pole %s' % (experiment_id, cs_w)
        lons_w, lats_w = np.meshgrid(lon_w, lat_w)
        lons_w,lats_w = iris.analysis.cartography.unrotate_pole(lons_w,lats_w, cs_w.grid_north_pole_longitude, cs_w.grid_north_pole_latitude)
        
        lon_w=lons_w[0]
        lat_w=lats_w[:,0]

        csur_w=cs_w.ellipsoid

        for i, coord in enumerate (u_wind.coords()):
            if coord.standard_name=='grid_latitude':
                lat_dim_coord_uwind = i
            if coord.standard_name=='grid_longitude':
                lon_dim_coord_uwind = i
       
        u_wind.remove_coord('grid_latitude')
        u_wind.remove_coord('grid_longitude')
        u_wind.add_dim_coord(iris.coords.DimCoord(points=lat_w, standard_name='grid_latitude', units='degrees', coord_system=csur_w),lat_dim_coord_uwind )
        u_wind.add_dim_coord(iris.coords.DimCoord(points=lon_w, standard_name='grid_longitude', units='degrees', coord_system=csur_w), lon_dim_coord_uwind)

        v_wind.remove_coord('grid_latitude')
        v_wind.remove_coord('grid_longitude')
        v_wind.add_dim_coord(iris.coords.DimCoord(points=lat_w, standard_name='grid_latitude', units='degrees', coord_system=csur_w), lat_dim_coord_uwind)
        v_wind.add_dim_coord(iris.coords.DimCoord(points=lon_w, standard_name='grid_longitude', units='degrees', coord_system=csur_w),lon_dim_coord_uwind )
        
    if isinstance(cs, iris.coord_systems.RotatedGeogCS):
        print ' 33.pp  - %s - Unrotate pole %s' % (experiment_id, cs)
        lons, lats = np.meshgrid(lon, lat) 

        lon_low= np.min(lons)
        lon_high = np.max(lons)
        lat_low = np.min(lats)
        lat_high = np.max(lats)

        lon_corners, lat_corners = np.meshgrid((lon_low, lon_high), (lat_low, lat_high))

        lons,lats = iris.analysis.cartography.unrotate_pole(lons,lats, cs.grid_north_pole_longitude, cs.grid_north_pole_latitude)
        lon_corner_u,lat_corner_u = iris.analysis.cartography.unrotate_pole(lon_corners, lat_corners, cs.grid_north_pole_longitude, cs.grid_north_pole_latitude)
        #lon_highu,lat_highu = iris.analysis.cartography.unrotate_pole(lon_high, lat_high, cs.grid_north_pole_longitude, cs.grid_north_pole_latitude)

        lon=lons[0]
        lat=lats[:,0]

        lon_low = lon_corner_u[0,0]
        lon_high = lon_corner_u[0,1]
        lat_low = lat_corner_u[0,0]
        lat_high = lat_corner_u[1,0]
                
        for i, coord in enumerate (oro.coords()):
            if coord.standard_name=='grid_latitude':
                lat_dim_coord_oro = i
            if coord.standard_name=='grid_longitude':
                lon_dim_coord_oro = i

        csur=cs.ellipsoid  
     
        oro.remove_coord('grid_latitude')
        oro.remove_coord('grid_longitude')
        oro.add_dim_coord(iris.coords.DimCoord(points=lat, standard_name='grid_latitude', units='degrees', coord_system=csur), lat_dim_coord_oro)
        oro.add_dim_coord(iris.coords.DimCoord(points=lon, standard_name='grid_longitude', units='degrees', coord_system=csur), lon_dim_coord_oro)

    else:

        lons, lats = np.meshgrid(lon, lat)
        lons_w, lats_w = np.meshgrid(lon_w, lat_w)

        lon_low= np.min(lons)
        lon_high = np.max(lons)
        lat_low = np.min(lats)
        lat_high = np.max(lats)
        
    print u_wind
        

############## Regrid and Difference #################################

  # Regrid Height and Temp/Specific humidity to global grid
    h_regrid = np.empty((len(lat_g_oro), len(lon_g_oro)))
    v_regrid = np.empty((len(lat_g_oro), len(lon_g_oro)))

    for p in plot_levels:
  
        ### Search for pressure level match
    
        s = np.searchsorted(p_levels[::-1], p)
        h_regrid = scipy.interpolate.griddata((lats.flatten(),lons.flatten()),mean_heights[:,:,-(s+1)].flatten() , (lats_go,lons_go),method='linear')


        v_regrid = scipy.interpolate.griddata((lats.flatten(),lons.flatten()),mean_var[:,:,-(s+1)].flatten() , (lats_go,lons_go),method='linear')
   
# Difference heights

        plt_h = np.where(np.isnan(h_regrid), np.nan, h_regrid - mean_heights_global[:,:,-(s+1)])

#Difference temperature/specific humidity
   
        plt_v = np.where(np.isnan(v_regrid), np.nan, v_regrid - mean_var_global[:,:,-(s+1)])
    
    # Set u,v for winds, linear interpolate to approx. 2 degree grid
        sc =  np.searchsorted(p_levs, p)

                        
    ##### Does not work on iris1.0 as on Leeds computers. Does work on later versions

        #u_interp = u_wind[sc,:,:]
        #v_interp = v_wind[sc,:,:].
        #sample_points = [('grid_latitude', np.arange(lat_low,lat_high,2)), ('grid_longitude', np.arange(lon_low,lon_high,2))]

        #u = iris.analysis.interpolate.linear(u_interp, sample_points, extrapolation_mode='linear')
        #v = iris.analysis.interpolate.linear(v_interp, sample_points).data


    ##### Does work on Iris 1.0

        # 2 degree lats lon lists for wind regridding on plot

        lat_wind_1deg = np.arange(lat_low,lat_high, 2)
        lon_wind_1deg = np.arange(lon_low,lon_high, 2)

    ### Regrid winds to global, difference, and then regrid to 2 degree spacing

        fl_la_lo = (lats_w.flatten(),lons_w.flatten())

        # print u_wind[sc,:,:].data.flatten().shape

        u_wind_rg_to_glob = scipy.interpolate.griddata(fl_la_lo, u_wind[sc,:,:].data.flatten(), (lats_g, lons_g), method='linear')
        v_wind_rg_to_glob = scipy.interpolate.griddata(fl_la_lo, v_wind[sc,:,:].data.flatten(), (lats_g, lons_g), method='linear')

        u_w=u_wind_rg_to_glob-u_global[sc,:,:].data
        v_w=v_wind_rg_to_glob-v_global[sc,:,:].data

        #u_interp = u_wind[sc,:,:].data
        #v_interp = v_wind[sc,:,:].data
       
        lons_wi, lats_wi = np.meshgrid(lon_wind_1deg, lat_wind_1deg)

        fl_la_lo = (lats_g.flatten(),lons_g.flatten())

        u = scipy.interpolate.griddata(fl_la_lo, u_w.flatten(), (lats_wi, lons_wi), method='linear')
        v = scipy.interpolate.griddata(fl_la_lo, v_w.flatten(), (lats_wi, lons_wi), method='linear')
       


        
#######################################################################################

### Plotting #########################################################################

        rc('font', family = 'serif', serif = 'cmr10')
        rc('text', usetex=True)

        #m_title = 'Height of %s-hPa level (m)' % (p)

# Set pressure height contour min/max
        if p == 925:
            clev_min = -24.
            clev_max = 24.
        elif p == 850:
            clev_min = -24.
            clev_max = 24.
        elif p == 700:
            clev_min = -24.
            clev_max = 24.
        elif p == 500:
            clev_min = -24.
            clev_max = 24.
        else:
            print 'Contour min/max not set for this pressure level'

# Set potential temperature min/max       
        if p == 925:
            clevpt_min = -3.
            clevpt_max = 3.
        elif p == 850:
            clevpt_min = -3.
            clevpt_max = 3.
        elif p == 700:
            clevpt_min = -3.
            clevpt_max = 3.
        elif p == 500:
            clevpt_min = -3.
            clevpt_max = 3.
        else:
            print 'Potential temperature min/max not set for this pressure level'

 # Set specific humidity min/max       
        if p == 925:
            clevsh_min = -0.0025
            clevsh_max = 0.0025
        elif p == 850:
            clevsh_min = -0.0025
            clevsh_max = 0.0025
        elif p == 700:
            clevsh_min = -0.0025
            clevsh_max = 0.0025
        elif p == 500:
            clevsh_min = -0.0025
            clevsh_max = 0.0025
        else:
            print 'Specific humidity min/max not set for this pressure level'

       #clevs_col = np.arange(clev_min, clev_max)
        clevs_lin = np.linspace(clev_min, clev_max, num=24)


        m =\
Basemap(llcrnrlon=lon_low,llcrnrlat=lat_low,urcrnrlon=lon_high,urcrnrlat=lat_high, rsphere = 6371229)

        #x, y = m(lons, lats)
        x,y = m(lons_go,lats_go)
        print x.shape
        x_w,y_w = m(lons_wi, lats_wi)
        fig=plt.figure(figsize=(8,8))
        ax = fig.add_axes([0.05,0.05,0.9,0.85], axisbg='#262626')

        m.drawcoastlines(color='#262626')  
        m.drawcountries(color='#262626')  
        m.drawcoastlines(linewidth=0.5)
        #m.fillcontinents(color='#CCFF99')
        m.drawparallels(np.arange(-80,81,10),labels=[1,1,0,0])
        m.drawmeridians(np.arange(0,360,10),labels=[0,0,0,1])
    
        cs_lin = m.contour(x,y, plt_h, clevs_lin,colors='#262626',linewidths=0.8)
       
        cmap=plt.cm.RdBu_r
        #cmap.set_bad('#262626', 1.)
        #cmap.set_over('#262626')
        #cmap.set_under('#262626')
        

        if plot_diag=='temp':
            
             plt_v = np.ma.masked_outside(plt_v, clevpt_max+20,  clevpt_min-20)

             cs_col = m.contourf(x,y, plt_v,  np.linspace(clevpt_min, clevpt_max), cmap=cmap, extend='both')
             
             cbar = m.colorbar(cs_col,location='bottom',pad="5%", format = '%d')  
             cbar.set_ticks(np.arange(clevpt_min,clevpt_max+2,2.))
             cbar.set_ticklabels(np.arange(clevpt_min,clevpt_max+2,2.))
             cbar.set_label('K')  
             #plt.suptitle('Difference from global model (Model - global ) of Height, Potential Temperature and Wind Vectors  at %s hPa'% (p), fontsize=10)  

        elif plot_diag=='sp_hum':
             plt_v = np.ma.masked_outside(plt_v, clevsh_max+20,  clevsh_min-20)

             cs_col = m.contourf(x,y, plt_v,  np.linspace(clevsh_min, clevsh_max), cmap=cmap, extend='both')
             cbar = m.colorbar(cs_col,location='bottom',pad="5%", format = '%.3f') 
             cbar.set_label('kg/kg')
             
             #plt.suptitle('Difference from global model (Model - Global Model ) of Height, Specific Humidity and Wind Vectors  at %s hPa'% (p), fontsize=10) 

        wind = m.quiver(x_w,y_w, u, v, scale=150,color='#262626' )
        qk = plt.quiverkey(wind, 0.1, 0.1, 5, '5 m/s', labelpos='W')
                
        plt.clabel(cs_lin, fontsize=10, fmt='%d', color='black')

        #plt.title('%s\n%s' % (m_title, model_name_convert_title.main(experiment_id)), fontsize=10)
        plt.savefig('/nfs/a90/eepdw/Figures/EMBRACE/%s/%s/geop_height_difference_globalLAM_%shPa_%s_%s_no_title.png' % (experiment_id, plot_diag, p, experiment_id, plot_diag), format='png', bbox_inches='tight')

        plt.title('\n'.join(wrap('%s' % (model_name_convert_title.main(experiment_id)), 75, replace_whitespace=False)), fontsize=16)
        
        #plt.show()  
        if not os.path.exists('/nfs/a90/eepdw/Figures/EMBRACE/%s/%s' % (experiment_id, plot_diag)): os.makedirs('/nfs/a90/eepdw/Figures/EMBRACE/%s/%s' % (experiment_id, plot_diag))
        plt.savefig('/nfs/a90/eepdw/Figures/EMBRACE/%s/%s/geop_height_difference_globalLAM_%shPa_%s_%s.png' % (experiment_id, plot_diag, p, experiment_id, plot_diag), format='png', bbox_inches='tight')


        #  Save fig - update dpi if need for printing
if __name__ == '__main__':
    main()
