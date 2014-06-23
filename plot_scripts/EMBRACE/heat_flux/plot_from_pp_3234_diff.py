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
import iris.analysis.cartography

import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

import scipy.interpolate

import datetime
from mpl_toolkits.basemap import cm

import imp
from textwrap import wrap

import re

import iris.analysis.cartography

import math

experiment_ids = ['djzny', 'djznq', 'djzns', 'dkjxq', 'dklyu', 'dkmbq', 'dklwu', 'dklzq' ]
#experiment_ids = ['dkjxq', 'dklyu', 'dkmbq', 'dklwu', 'dklzq' ]
save_path='/nfs/a90/eepdw/Mean_State_Plot_Data/Figures/'

model_name_convert_title = imp.load_source('util', '/nfs/see-fs-01_users/eepdw/python_scripts/model_name_convert_title.py')
unrotate = imp.load_source('util', '/nfs/see-fs-01_users/eepdw/python_scripts/unrotate_pole.py')
pp_file = '3234_mean'

degs_crop_top = 1.7
degs_crop_bottom = 2.5

min_contour = -50
max_contour = 50
tick_interval=20
#
# cmap= cm.s3pcpn_l

divisor=10  # for lat/lon rounding

def main():
 # Min and max lats lons from smallest model domain (dkbhu) - see spreadsheet

 latmin=-6.79
 latmax=29.721
 lonmin=340.
 lonmax=379.98

 lat_constraint=iris.Constraint(grid_latitude= lambda la: latmin <= la.point <= latmax)
 lon_constraint=iris.Constraint(grid_longitude= lambda lo: lonmin <= lo.point <= lonmax)

# Global LAM not rotated - so different coord constraints
 lonmin_g=64.1153327
 lonmax_g=101.865817

 lon_constraint_g = iris.Constraint(grid_longitude= lambda lo: lonmin_g <= lo.point <= lonmax_g)
# Load global cube
                                
 gl = '/nfs/a90/eepdw/Data/EMBRACE/Mean_State/pp_files/djzn/djznw/%s.pp' % pp_file
 glob = iris.load_cube(gl, lat_constraint & lon_constraint_g)
 #glob = iris.load_cube(gl)
 cs_glob = glob.coord_system('CoordSystem')
 # Unrotate global cube

 lat_g = glob.coord('grid_latitude').points
 lon_g = glob.coord('grid_longitude').points

 #print lat_g
 if isinstance(cs_glob, iris.coord_systems.RotatedGeogCS):
        print ' Global Model - djznw - Unrotate pole %s' % cs_glob
        lons_g, lats_g = np.meshgrid(lon_g, lat_g)
        lons_g,lats_g = iris.analysis.cartography.unrotate_pole(lons_g,lats_g, cs_glob.grid_north_pole_longitude, cs_glob.grid_north_pole_latitude)
        
        lon_g=lons_g[0]
        lat_g=lats_g[:,0]

        #print lats_g

 for i, coord in enumerate (glob.coords()):
            if coord.standard_name=='grid_latitude':
                lat_dim_coord_glob = i
            if coord.standard_name=='grid_longitude':
                lon_dim_coord_glob = i

 csur_glob=cs_glob.ellipsoid
 glob.remove_coord('grid_latitude')
 glob.remove_coord('grid_longitude')
 glob.add_dim_coord(iris.coords.DimCoord(points=lat_g, standard_name='grid_latitude', units='degrees', coord_system=csur_glob), lat_dim_coord_glob)
 glob.add_dim_coord(iris.coords.DimCoord(points=lon_g, standard_name='grid_longitude', units='degrees', coord_system=csur_glob), lon_dim_coord_glob)


 experiment_ids = ['djzny', 'djznq', 'djzns', 'djznw', 'dkjxq', 'dklyu', 'dkmbq', 'dklwu', 'dklzq' ] 
 
 #experiment_ids = ['djzny' ] 
 for experiment_id in experiment_ids:

  expmin1 = experiment_id[:-1]
  pfile = '/nfs/a90/eepdw/Data/EMBRACE/Mean_State/pp_files/%s/%s/%s.pp' % (expmin1, experiment_id, pp_file)

     #pc =  iris(pfile)
  #pcube = iris.load_cube(pfile, lat_constraint & lon_constraint)
  pcube = iris.load_cube(pfile)
  #print pcube
     #print pc
 
 # Get min and max latitude/longitude and unrotate  to get min/max corners to crop plot automatically - otherwise end with blank bits on the edges 

 # Unrotate cube

  lat = pcube.coord('grid_latitude').points
  lon = pcube.coord('grid_longitude').points

  #print lat
  #print 'lat'
  #print lon
  cs = pcube.coord_system('CoordSystem')

  if isinstance(cs, iris.coord_systems.RotatedGeogCS):
        print ' %s - Unrotate pole %s' % (experiment_id,cs)

        lons, lats = np.meshgrid(lon, lat)

        lons,lats = iris.analysis.cartography.unrotate_pole(lons,lats, cs.grid_north_pole_longitude, cs.grid_north_pole_latitude)
       
        lon=lons[0]
        lat=lats[:,0]
    
        for i, coord in enumerate (pcube.coords()):
            if coord.standard_name=='grid_latitude':
                lat_dim_coord = i
            if coord.standard_name=='grid_longitude':
                lon_dim_coord = i

        csur=cs.ellipsoid
        
        pcube.remove_coord('grid_latitude')
        pcube.remove_coord('grid_longitude')
        pcube.add_dim_coord(iris.coords.DimCoord(points=lat, standard_name='grid_latitude', units='degrees', coord_system=csur), lat_dim_coord)
        pcube.add_dim_coord(iris.coords.DimCoord(points=lon, standard_name='grid_longitude', units='degrees', coord_system=csur), lon_dim_coord)

  lon_min=np.min(lons_g)
  lon_max=np.max(lons_g)
  
  lon_low_tick=lon_min -(lon_min%divisor)
  lon_high_tick=math.ceil(lon_max/divisor)*divisor

  lat_min=np.min(lats_g)
  lat_max=np.max(lats_g)
  lat_low_tick=lat_min - (lat_min%divisor)
  lat_high_tick=math.ceil(lat_max/divisor)*divisor
 
  print lon_high_tick
  print lon_low_tick

  pcube_regrid_data = scipy.interpolate.griddata((lats.flatten(), lons.flatten()),pcube.data.flatten(), (lats_g, lons_g), method='linear')

  
  #pcube_regrid = iris.analysis.interpolate.linear(pcube, sample_points)
  #print pcube.data.flatten()
  pcube_regrid = glob.copy(data=pcube_regrid_data)
  pcubediff=pcube_regrid-glob
  #print pcube.data[0,0]
  #print pcube_regrid_data[0,0]
  #print pcubediff.data
  #print glob.data[0,0]
  plt.figure(figsize=(8,8))
         
  cmap= cmap=plt.cm.RdBu_r
    
  ax = plt.axes(projection=ccrs.PlateCarree(), extent=(lonmin_g+2,lonmax_g-2,latmin+degs_crop_bottom,latmax-degs_crop_top))
  
  clevs = np.linspace(min_contour, max_contour,256)
  cont = iplt.contourf(pcubediff, clevs, cmap=cmap, extend='both')
                     
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
  
  gl.xlabel_style = {'size': 12, 'color':'black'}
  #gl.xlabel_style = {'color': '#262626', 'weight': 'bold'}
  gl.ylabel_style = {'size': 12, 'color':'black'}         

  cbar = plt.colorbar(cont, orientation='horizontal', pad=0.05, extend='both', format = '%d')
  #cbar.set_label('') 
  cbar.set_label(pcube.units, fontsize=10)
  cbar.set_ticks(np.arange(min_contour, max_contour+tick_interval,tick_interval))
  ticks = (np.arange(min_contour, max_contour+tick_interval,tick_interval))
  cbar.set_ticklabels(['%d' % i for i in ticks])
  main_title='%s - Difference' % pcube.standard_name.title().replace('_',' ')
  model_info=re.sub('(.{68} )', '\\1\n', str(model_name_convert_title.main(experiment_id)), 0, re.DOTALL)
  model_info = re.sub(r'[(\']', ' ', model_info)
  model_info = re.sub(r'[\',)]', ' ', model_info)
  print model_info
  
  if not os.path.exists('%s%s/%s' % (save_path, experiment_id, pp_file)): os.makedirs('%s%s/%s' % (save_path, experiment_id, pp_file))

  plt.savefig('%s%s/%s/%s_%s_notitle_diff.png' % (save_path, experiment_id, pp_file, experiment_id, pp_file), format='png', bbox_inches='tight')

  plt.title('\n'.join(wrap('%s\n%s' % (main_title, model_info), 1000,replace_whitespace=False)), fontsize=16)
 
  #plt.show()
 
  plt.savefig('%s%s/%s/%s_%s_diff.png' % (save_path, experiment_id, pp_file, experiment_id, pp_file), format='png', bbox_inches='tight')
  
  plt.close()
 

if __name__ == '__main__':
   main()
