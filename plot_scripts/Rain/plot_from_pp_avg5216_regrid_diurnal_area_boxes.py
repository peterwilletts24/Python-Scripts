"""

Load pp, plot and save


"""

import os, sys

import matplotlib

matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
from matplotlib import rc
from matplotlib.font_manager import FontProperties
from matplotlib import rcParams

from mpl_toolkits.basemap import Basemap

from matplotlib.patches import Polygon

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
pp_file = 'rain_mean_regrid'

degs_crop_top = 1.7
degs_crop_bottom = 2.5

min_contour = 0
max_contour = 3.5
tick_interval=0.5
#
# cmap= cm.s3pcpn_l

divisor=10  # for lat/lon rounding



def main():

 def draw_screen_poly_iris( lats, lons):
    #x, y = m( lons, lats )
    xy = zip(lons,lats)
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
    x, y = (lons_label, lats_label ) 
    #plt.text(x, y, text, color='#262626', ha="center", va="center", size=32, backgroundcolor='white', alpha=0.4 )
    font0 = FontProperties()
    font0.set_family('sans-serif')
   
    plt.text(x, y, text, color='black', ha="center", va="center", size=64 , fontweight='bold', fontproperties=font0)
    if (plot_coords[2]=='Southern Indian Ocean'):
       plt.text(x, y, text, color='red', ha="center", va="center", size=64, fontproperties=font0)

#experiment_ids = ['djzny', 'djzns', 'djznw', 'dkjxq', 'dklyu', 'dkmbq', 'dklwu', 'dklzq' ] 
 #experiment_ids = ['djznq', 'djzny', 'djzns', 'djznu', 'dkbhu', 'dkjxq', 'dklyu', 'dkmbq', 'dklwu', 'dklzq', 'dkhgu'] 
 experiment_ids = ['djznq' ] 
 #experiment_ids = ['dkhgu','dkjxq']
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
 for experiment_id in experiment_ids:

  expmin1 = experiment_id[:-1]
  pfile = '/nfs/a90/eepdw/Data/EMBRACE/Mean_State/pp_files/%s/%s/%s.pp' % (expmin1, experiment_id, pp_file)
  if (experiment_id=='djznq'):
     pfile = '/nfs/a90/eepdw/Data/EMBRACE/Mean_State/pp_files/%s/%s/rain_mean.pp' % (expmin1, experiment_id)

  legendEntries=[]
  legendtext=[]

  #pc =  iris(pfile)
  pcube = iris.load_cube(pfile)
  print pcube
     #print pc
 
 # Get min and max latitude/longitude and unrotate  to get min/max corners to crop plot automatically - otherwise end with blank bits on the edges 
  if (experiment_id=='djznq'):
   lats = pcube.coord('grid_latitude').points
   lons = pcube.coord('grid_longitude').points
  else:
   lats = pcube.coord('latitude').points
   lons = pcube.coord('longitude').points
  
  cs = pcube.coord_system('CoordSystem')
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

  #lon_low= 62
  #lon_high = 102
  #lat_low = -7
  #lat_high = 33

  #lon_high_box = 101.866 
  #lon_low_box = 64.115
  #lat_high_box = 33.
  #lat_low_box =-6.79

  #lon_high = 101.866 
  #lon_low = 64.115
  #lat_high = 33.
  #lat_low =-6.79

  lon_low_tick=lon_low -(lon_low%divisor)
  lon_high_tick=math.ceil(lon_high/divisor)*divisor

  lat_low_tick=lat_low - (lat_low%divisor)
  lat_high_tick=math.ceil(lat_high/divisor)*divisor
 
  print lat_high_tick
  print lat_low_tick
  plt.figure(figsize=(8,8))
         
  cmap=cm.s3pcpn_l
    
  ax = plt.axes(projection=ccrs.PlateCarree(), extent=(lon_low,lon_high,lat_low+degs_crop_bottom,lat_high-degs_crop_top))
  
  #ax = plt.axes(projection=ccrs.PlateCarree(), extent=(lon_low,lon_high,lat_low,lat_high))

  #ax = plt.axes(projection=ccrs.PlateCarree())

  clevs = np.linspace(min_contour, max_contour,256)

  pcubeplot=iris.analysis.maths.multiply(pcube,3600)

  cont = iplt.contourf(pcubeplot, clevs, cmap=cmap, extend='both')
                     
  #plt.clabel(cont, fmt='%d')
  #ax.stock_img()
  ax.coastlines(resolution='110m', color='#262626') 
                     
  gl = ax.gridlines(draw_labels=True,linewidth=0.5, color='#262626', alpha=0.5, linestyle='--')
  gl.xlabels_top = False
  gl.ylabels_right = False
            #gl.xlines = False
  dx, dy = 10, 10

  gl.xlocator = mticker.FixedLocator(range(int(lon_low_tick),int(lon_high_tick)+dx,dx))
  gl.ylocator = mticker.FixedLocator(range(int(lat_low_tick),int(lat_high_tick)+dy,dy))
  gl.xformatter = LONGITUDE_FORMATTER
  gl.yformatter = LATITUDE_FORMATTER
  
  gl.xlabel_style = {'size': 12, 'color':'#262626'}
  #gl.xlabel_style = {'color': '#262626', 'weight': 'bold'}
  gl.ylabel_style = {'size': 12, 'color':'#262626'}         

  cbar = plt.colorbar(cont, orientation='horizontal', pad=0.05, extend='both')
  cbar.set_label('mm/h', fontsize=10, color='#262626') 
  #cbar.set_label(pcube.units, fontsize=10, color='#262626')
  cbar.set_ticks(np.arange(min_contour, max_contour+tick_interval,tick_interval))
  ticks = (np.arange(min_contour, max_contour+tick_interval,tick_interval))
  cbar.set_ticklabels(['%.1f' % i for i in ticks])
  cbar.ax.tick_params(labelsize=10, color='#262626')
  
  main_title='Mean Rainfall for EMBRACE Period (smoothed to 24km)'
  #main_title=pcube.standard_name.title().replace('_',' ')
  model_info=re.sub('(.{68} )', '\\1\n', str(model_name_convert_title.main(experiment_id)), 0, re.DOTALL)
  #model_info = re.sub(r'[(\']', ' ', model_info)
  #model_info = re.sub(r'[\',)]', ' ', model_info)
  #print model_info
  for l,plot_coords in enumerate(zip(lats_poly,lons_poly, labels)):
    #colour = cmap(1.*(l*2)/(NUM_COLOURS*2))
    draw_screen_poly_iris( plot_coords[0], plot_coords[1])
    label(plot_coords[0], plot_coords[1], l+1)

  if not os.path.exists('%s%s/%s' % (save_path, experiment_id, pp_file)): os.makedirs('%s%s/%s' % (save_path, experiment_id, pp_file))

  plt.savefig('%s%s/%s/%s_%s_area_boxes_notitle.png' % (save_path, experiment_id, pp_file, experiment_id, pp_file), format='png', bbox_inches='tight')

  plt.title('\n'.join(wrap('%s\n%s' % (main_title, model_info), 1000,replace_whitespace=False)), fontsize=16)
 
  #plt.show()
 
  plt.savefig('%s%s/%s/%s_%s_area_boxes.png' % (save_path, experiment_id, pp_file, experiment_id, pp_file), format='png', bbox_inches='tight')
  
  plt.close()
 

if __name__ == '__main__':
   main()
