"""

Load geopotential heights and orography cube to get lat/lon cross-section

13/08/2014
"""

import os, sys

import pdb

import iris
import iris.analysis.cartography

from update_pp_cube_coords import update_coords

#import h5py

import numpy as np

#c_section_lon=74.

c_lon_min=75.
c_lon_max=85.
gap=1.

c_section_lat=0

diagnostic=4
#experiment_ids = ['djznw', 'djzny', 'djznq', 'dklyu', 'dkmbq', 'dklwu', 'dklzq', 'djzns'  ]
#experiment_ids = ['djznq', 'dklyu', 'dkmbq', 'dklzq', 'djzns'  ] #djznw and dklwu missing
#experiment_ids = ['dklyu', 'dkmbq', 'dklwu', 'dklzq' ]
#experiment_ids = ['dkbhu']
experiment_ids = ['djznw']

for experiment_id in experiment_ids:

 expmin1 = experiment_id[:-1]

 diag = iris.load_cube('/nfs/a90/eepdw/Data/EMBRACE/%s/%s/%s.pp' % (expmin1,experiment_id, diagnostic))
 
 #diag = iris.load(f_diag)

 print diag

 cs = diag.coord_system('CoordSystem')
 print cs
 csur=cs.ellipsoid  

 lat = diag.coord('grid_latitude').points
 lon = diag.coord('grid_longitude').points

 lons, lats = np.meshgrid(lon, lat)  

 lons,lats = iris.analysis.cartography.unrotate_pole(lons,lats, cs.grid_north_pole_longitude, cs.grid_north_pole_latitude)

 lon=lons[0]
 lat=lats[:,0]

 for i, coord in enumerate (diag.coords()):
     if coord.standard_name=='grid_latitude':
         lat_dim_coord_diag = i
     if coord.standard_name=='grid_longitude':
         lon_dim_coord_diag = i

 diag.remove_coord('grid_latitude')
 diag.remove_coord('grid_longitude')
 diag.add_dim_coord(iris.coords.DimCoord(points=lat, standard_name='grid_latitude', units='degrees', coord_system=csur), lat_dim_coord_diag)
 diag.add_dim_coord(iris.coords.DimCoord(points=lon, standard_name='grid_longitude', units='degrees', coord_system=csur), lon_dim_coord_diag)

 for c_section_lon in np.arange(c_lon_min,c_lon_max+1, gap):

  if (c_section_lon != 0):
   
    print c_section_lon
    l=diag.coord('grid_longitude').nearest_neighbour_index(c_section_lon) 

    if lon_dim_coord_diag==0:
       xc=diag[l,:]  
    if lon_dim_coord_diag==1:
       xc=diag[:,l,:]
    if lon_dim_coord_diag==2:
       xc=diag[:,:,l,:]
    if lon_dim_coord_diag==3:
       xc=diag[:,:,:,l,:]
  
    iris.save(xc, '/nfs/a90/eepdw/Figures/EMBRACE/Cross_Sections/%s_%s_height_XC_Longitude_%s.pp' % (experiment_id, diagnostic, str(c_section_lon).replace(".", "")))


   #THESE METHODS MIGHT WORK BUT TAKE A LONG TIME - I THINK BECAUSE THEY LOAD THE WHOLD CUBE IN TO INDEX
    #xc = iris.analysis.interpolate.extract_nearest_neighbour(diag, [('grid_longitude', c_section_lon)]).data
    #lon_slice = iris.analysis.interpolate.linear(diag, [('grid_longitude', l), ('grid_latitude', np.linspace(20, 30, 50))])
    #print lon_slice
    #pdb.set_trace
    #iris.save(lon_slice, '/nfs/a90/eepdwCross_Sections/%s_%s_height_XC_Longitude_%s.pp' % (experiment_id, diagnostic, str(c_section_lon).replace(".", "")))
    #iris.save(iris.analysis.interpolate.extract_nearest_neighbour(diag, [('grid_longitude', c_section_lon)]), 
     #         '/nfs/a90/eepdwCross_Sections/%s_%s_height_XC_Longitude_%s.pp' 
      #        % (experiment_id, diagnostic, str(c_section_lon).replace(".", "")))
   
    #xc=lon_slice.data
    #np.savez('/nfs/a90/eepdwCross_Sections/%s_%s_height_XC_Longitude_%s' % (experiment_id, diag, c_section_lon), xc=xc, coord=diag.coord('grid_latitude').points)
