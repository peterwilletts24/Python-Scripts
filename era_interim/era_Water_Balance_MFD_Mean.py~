import os, sys

import matplotlib

matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
from matplotlib import rc
from matplotlib.font_manager import FontProperties
from matplotlib import rcParams
from matplotlib.colors import from_levels_and_colors

from mpl_toolkits.basemap import Basemap

rc('font', family = 'serif', serif = 'cmr10')
rc('text', usetex=True)

rcParams['text.usetex']=True
rcParams['text.latex.unicode']=True
rcParams['font.family']='serif'
rcParams['font.serif']='cmr10'

import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as mpl_cm
import numpy as np

import iris
import iris.coords as coords
import iris.quickplot as qplt
import iris.plot as iplt
import iris.coord_categorisation

import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

import datetime
from mpl_toolkits.basemap import cm

import imp
from textwrap import wrap

import re

import iris.analysis.cartography

import math

from netCDF4 import Dataset

#nc = Dataset('/nfs/a90/eepdw/Data/ERA_Iinterim_Heat_Rad_Fluxes/era_interim_netcdf_heat_rad_flux_evap_precip_6hr_timestep.nc')
nc12 = Dataset('/nfs/a90/eepdw/Data/ERA_Iinterim_Heat_Rad_Fluxes/era_interim_netcdf_heat_rad_flux_evap_precip_00_timestep.nc')

time_min=datetime.datetime(2011,8,18,0,0,0,0)
time_max=datetime.datetime(2011,9,8,0,0,0,0)

min_contour = -15
max_contour = 10
tick_interval=5

lon_max = 116 
lon_min = 34

lat_max= 40.
lat_min=-11.25

lon_high_plot = 102
lon_low_plot = 64
lat_high_plot= 30.
lat_low_plot=-10

divisor=10  # for lat/lon rounding

#latent_mean = np.load('/nfs/a90/eepdw/Data/Saved_data/era_i/Fluxes/era_i_emb_latent_mean.npy')
# sensible_mean = np.load('/nfs/a90/eepdw/Data/Saved_data/era_i/Fluxes/era_i_emb_sensible_mean.npy')
# swave_mean = np.load('/nfs/a90/eepdw/Data/Saved_data/era_i/Fluxes/era_i_emb_swave_mean.npy')
# lwave_mean = np.load('/nfs/a90/eepdw/Data/Saved_data/era_i/Fluxes/era_i_emb_lwave_mean.npy')
# total_mean = np.load('/nfs/a90/eepdw/Data/Saved_data/era_i/Fluxes/era_i_emb_total_mean.npy')
lat = np.load('/nfs/a90/eepdw/Data/Saved_data/era_i/Fluxes/era_i_emb_lats.npy')
lon = np.load('/nfs/a90/eepdw/Data/Saved_data/era_i/Fluxes/era_i_emb_longs.npy')

lons,lats = np.meshgrid(lon, lat)

hours_since=datetime.datetime(1900,1,1,0,0,0,0)
# Get min and max  index positions for latitude and longitude

la_index = np.where((nc12.variables['latitude'][:]<=lat_max) & (nc12.variables['latitude'][:] >= lat_min))
lo_index = np.where((nc12.variables['longitude'][:]<=lon_max) & (nc12.variables['longitude'][:] >= lon_min))

la_i_max = np.max(la_index)
la_i_min = np.min(la_index)
lo_i_max = np.max(lo_index)
lo_i_min = np.min(lo_index)

lat_amounts=la_i_max-la_i_min
lon_amounts=lo_i_max-lo_i_min

print nc12

# Load evaporation and precipitation accumulations (in metres)
# datetimes = np.array([datetime.timedelta(hours=float(i))+hours_since for i in nc.variables['time'][:]])
# time_index= np.where((datetimes<=time_max) & (datetimes >= time_min))

# t_i_max = np.max(time_index)
# t_i_min = np.min(time_index) 

#evap_in = nc.variables['e'][t_i_min:t_i_max+1,la_i_min:la_i_max+1, lo_i_min:lo_i_max+1]
#precip_in = nc.variables['tp'][t_i_min:t_i_max+1,la_i_min:la_i_max+1, lo_i_min:lo_i_max+1]

datetimes = np.array([datetime.timedelta(hours=float(i))+hours_since for i in nc12.variables['time'][:]])
time_index= np.where((datetimes<=time_max) & (datetimes >= time_min))

t_i_max = np.max(time_index)
t_i_min = np.min(time_index) 

precip_in = nc12.variables['tp'][t_i_min:t_i_max+1,la_i_min:la_i_max+1, lo_i_min:lo_i_max+1]
#evap_in = nc12.variables['e'][t_i_min:t_i_max+1,la_i_min:la_i_max+1, lo_i_min:lo_i_max+1]
# Evaporation and precipitaiton fields from ERA-Interim are in 12 hourly accumulations, and in metres, E-P for the EMBRACE models have been calculated in kg/m^2/day) - kg/m^2/day = mm/day - (m * 1000)*(24/12) = kg/m^2/day - based on how much accumulates in that 12 hours, from which the average will be calculated
# By looking at data, it appears that precipitation accumulates for 12 hours -ECMWF site link to documentation broken currently

#evap_mean = np.mean((evap_in * 1000)*(24/12), axis=0, dtype=np.float64)
precip_mean = np.mean((precip_in * 1000)*(24/12), axis=0, dtype=np.float64)

# Evaporation same method as in EMBRACE calculation

latent_mean = np.load('/nfs/a90/eepdw/Data/Saved_data/era_i/Fluxes/era_i_emb_latent_mean.npy')
latent_heat_of_vapourisation = 2.5E06
convert_to_kg_m2_day = 86400


evap_rate = (latent_mean/latent_heat_of_vapourisation)*86400
# rain_daily=pcuberain*convert_to_kg_m2_day


waterbalance = evap_rate - precip_mean
m =\
Basemap(llcrnrlon=lon_low_plot,llcrnrlat=lat_low_plot,urcrnrlon=lon_high_plot,urcrnrlat=lat_high_plot,projection='mill', rsphere=6371229)

x, y = m(lons, lats)
fig=plt.figure(figsize=(8,8))
ax = fig.add_axes([0.05,0.05,0.9,0.85])

cmap=plt.cm.RdBu_r

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

cs_col = m.contourf(x,y, waterbalance, clevs, cmap=cmap, norm=norm, extend='both')
cbar = m.colorbar(cs_col,location='bottom',pad="5%")
                               #cbar.ax.tick_params(labelsize=12,  colors='#262626')

ticks= np.arange(int(min_contour),int(max_contour)+tick_interval,tick_interval)
cbar.set_ticks(ticks, update_ticks=True)
cbar.set_ticklabels(([r"${%s}$" % x for x in ticks]))

cbar.set_label('$kgm^{-2}day^{-1}$', fontsize=12, color='#262626')  

#plt.show()

plt.savefig('/nfs/a90/eepdw/Figures/ERA_Interim/Era_Interim_EMBRACE_period_evap_minus_precip_mean_notitle.png', format='png', bbox_inches='tight')

plt.suptitle('ERA-Interim Reanalysis Mean Evaporation - Precipitation for EMBRACE period', fontsize=16, color='#262626')

plt.savefig('/nfs/a90/eepdw/Figures/ERA_Interim/Era_Interim_EMBRACE_period_evap_minus_precip_mean.png', format='png', bbox_inches='tight')
