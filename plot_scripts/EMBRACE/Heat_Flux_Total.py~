"""

Load pp, plot and save


"""

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

save_path='/nfs/a90/eepdw/Figures/EMBRACE/'

model_name_convert_title = imp.load_source('util', '/nfs/see-fs-01_users/eepdw/python_scripts/model_name_convert_title.py')
unrotate = imp.load_source('util', '/nfs/see-fs-01_users/eepdw/python_scripts/unrotate_pole.py')
pp_file1201 = '1201_mean'
pp_file2201 = '2201_mean'
pp_file3217 = '3217_mean'
pp_file3234 = '3234_mean'

pp_file='total_downward_durface_heat_flux_mean'

degs_crop_top = 1.7
degs_crop_bottom = 2.5

min_contour = -30
max_contour = 160
tick_interval=20
#
# cmap= cm.s3pcpn_l

divisor=10  # for lat/lon rounding

def main():

 experiment_ids = ['djzny','djznq', 'djzns', 'djznw', 'dkjxq', 'dklyu', 'dkmbq', 'dklwu', 'dklzq' ] 
 
 #experiment_ids = ['djzny' ] 
 for experiment_id in experiment_ids:

  expmin1 = experiment_id[:-1]
  pfile1201 = '/nfs/a90/eepdw/Data/EMBRACE/Mean_State/pp_files/%s/%s/%s.pp' % (expmin1, experiment_id, pp_file1201)
  pfile2201 = '/nfs/a90/eepdw/Data/EMBRACE/Mean_State/pp_files/%s/%s/%s.pp' % (expmin1, experiment_id, pp_file2201)
  pfile3217 = '/nfs/a90/eepdw/Data/EMBRACE/Mean_State/pp_files/%s/%s/%s.pp' % (expmin1, experiment_id, pp_file3217)
  pfile3234 = '/nfs/a90/eepdw/Data/EMBRACE/Mean_State/pp_files/%s/%s/%s.pp' % (expmin1, experiment_id, pp_file3234)
     #pc =  iris(pfile)
  pcube1201 = iris.load_cube(pfile1201)
  pcube2201 = iris.load_cube(pfile2201)
  pcube3217 = iris.load_cube(pfile3217)
  pcube3234 = iris.load_cube(pfile3234)

  pcubetotal = pcube1201 + pcube2201 - pcube3217 - pcube3234
 
  print pcubetotal
     #print pc
 
 # Get min and max latitude/longitude and unrotate  to get min/max corners to crop plot automatically - otherwise end with blank bits on the edges 
  lats = pcubetotal.coord('grid_latitude').points
  lons = pcubetotal.coord('grid_longitude').points
  
  cs = pcubetotal.coord_system('CoordSystem')
  if isinstance(cs, iris.coord_systems.RotatedGeogCS):

      print 'Rotated CS %s' % cs
     
      lon_low= np.min(lons)
      lon_high = np.max(lons)
      lat_low = np.min(lats)
      lat_high = np.max(lats)

      lon_corners, lat_corners = np.meshgrid((lon_low, lon_high), (lat_low, lat_high))
      
      lon_corner_u,lat_corner_u = unrotate.unrotate_pole(lon_corners, lat_corners, cs.grid_north_pole_longitude, cs.grid_north_pole_latitude)
      lon_low = lon_corner_u[0,0]
      lon_high = lon_corner_u[0,1]
      lat_low = lat_corner_u[0,0]
      lat_high = lat_corner_u[1,0]

  else: 
      lon_low= np.min(lons)
      lon_high = np.max(lons)
      lat_low = np.min(lats)
      lat_high = np.max(lats)

  lon_low_tick=lon_low -(lon_low%divisor)
  lon_high_tick=math.ceil(lon_high/divisor)*divisor

  lat_low_tick=lat_low - (lat_low%divisor)
  lat_high_tick=math.ceil(lat_high/divisor)*divisor
 
  print lat_high_tick
  print lat_low_tick
  plt.figure(figsize=(8,8))
         
  cmap=plt.cm.RdBu_r

  ax = plt.axes(projection=ccrs.PlateCarree(), extent=(lon_low,lon_high,lat_low+degs_crop_bottom,lat_high-degs_crop_top))

  clevs = np.linspace(min_contour, max_contour,256)

  midpoint=0
  midp = np.mean(np.c_[clevs[:-1], clevs[1:]], axis=1)

  vals = np.interp(midp, [min_contour, midpoint, max_contour], [0, 0.5, 1])
  cols = plt.cm.RdBu_r(vals)

  clevs_extend = np.linspace(min_contour, max_contour,254)
  cmap, norm = from_levels_and_colors(clevs_extend, cols, extend='both')  
 
  cont = iplt.contourf(pcubetotal, clevs, cmap=cmap, extend='both', norm=norm)
  #cont = iplt.contourf(pcubetotal, cmap=cmap, extend='both')                  
  #plt.clabel(cont, fmt='%d')
  #ax.stock_img()
  ax.coastlines(resolution='110m', color='#262626') 
                     
  gl = ax.gridlines(draw_labels=True,linewidth=0.5, color='#262626', alpha=0.5, linestyle='--')
  gl.xlabels_top = False
  gl.ylabels_right = False
            #gl.xlines = False

  gl.xlocator = mticker.FixedLocator(range(int(lon_low_tick),int(lon_high_tick)+divisor,divisor))
  gl.ylocator = mticker.FixedLocator(range(int(lat_low_tick),int(lat_high_tick)+divisor,divisor))
  gl.xformatter = LONGITUDE_FORMATTER
  gl.yformatter = LATITUDE_FORMATTER
  
  gl.xlabel_style = {'size': 12, 'color':'#262626'}
  #gl.xlabel_style = {'color': '#262626', 'weight': 'bold'}
  gl.ylabel_style = {'size': 12, 'color':'#262626'}         

  cbar = plt.colorbar(cont, orientation='horizontal', pad=0.05, extend='both', format = '$%d$')
  cbar.set_label('$W m^{-2}$') 
  #cbar.set_label(pcubetotal.units, fontsize=10)
  cbar.set_ticks(np.arange(min_contour, max_contour+tick_interval,tick_interval))
  ticks = (np.arange(min_contour, max_contour+tick_interval,tick_interval))
  cbar.set_ticklabels(['${%d}$' % i for i in ticks])
  #main_title=pcubetotal.standard_name.title().replace('_',' ')
  #main_title=('Total 

  #model_info=re.sub('(.{68} )', '\\1\n', str(model_name_convert_title.main(experiment_id)), 0, re.DOTALL)
  #model_info = re.sub(r'[(\']', ' ', model_info)
  #model_info = re.sub(r'[\',)]', ' ', model_info)
  #print model_info
  
  if not os.path.exists('%s%s/%s' % (save_path, experiment_id, pp_file)): os.makedirs('%s%s/%s' % (save_path, experiment_id, pp_file))

  plt.savefig('%s%s/%s/%s_%s_notitle.png' % (save_path, experiment_id, pp_file, experiment_id, pp_file), format='png', bbox_inches='tight')

  #plt.title('\n'.join(wrap('%s\n%s' % (main_title, model_info), 1000,replace_whitespace=False)), fontsize=16)
 
  #plt.show()
 
  #plt.savefig('%s%s/%s/%s_%s.png' % (save_path, experiment_id, pp_file, experiment_id, pp_file), format='png', bbox_inches='tight')
  
  plt.close()
 

if __name__ == '__main__':
   main()
