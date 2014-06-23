"""
Regrid pp file to TRMM and difference

"""

import os, sys
import datetime

import iris
import iris.unit as unit

import scipy.interpolate

import cPickle as pickle

import numpy as np

diag = 'rain_mean'

pp_file_path='/nfs/a90/eepdw/Data/EMBRACE/Mean_State/pp_files/'

regrid_model='trmm'
fg = '%sdjzn/djznw/%s.pp' % (pp_file_path, diag)

sum_dom, latitude_domsingle, longitude_domsingle= pickle.load(open('/nfs/a90/eepdw/Data/Saved_data/TRMM/trmm_emb_pcpmean.p', 'rb'))

lons_trmm= longitude_domsingle[:]
lats_trmm = latitude_domsingle[:]
lons_trmm,lats_trmm = np.meshgrid(lons_trmm, lats_trmm)

#experiment_ids = ['djznw', 'djznq', 'djzny', 'djzns', 'dkmbq', 'dklyu', 'dklwu', 'dklzq', 'dkbhu', 'djznu', 'dkhgu' ]
experiment_ids = ['dkjxq' ]
def unrotate_pole(rotated_lons, rotated_lats, pole_lon, pole_lat):
     import cartopy.crs as ccrs

     """
      Convert rotated-pole lons and lats to unrotated ones.

      Example::

      lons, lats = unrotate_pole(grid_lons, grid_lats, pole_lon, pole_lat)

      .. note:: Uses proj.4 to perform the conversion.

      """
     src_proj = ccrs.RotatedGeodetic(pole_longitude=pole_lon,
                                    pole_latitude=pole_lat)
     target_proj = ccrs.Geodetic()
     res = target_proj.transform_points(x=rotated_lons, y=rotated_lats,
                                       src_crs=src_proj)
     unrotated_lon = res[..., 0]
     unrotated_lat = res[..., 1]

     return unrotated_lon, unrotated_lat

def unrotate_and_update_cube(rot_cube):
    
    import iris
    import numpy as np

    latr = rot_cube.coord('grid_latitude').points
    lonr = rot_cube.coord('grid_longitude').points
    #p_levs = rot_cube.coord('pressure').points
    
    cs = rot_cube.coord_system('CoordSystem')

    if isinstance(cs, iris.coord_systems.RotatedGeogCS):

        print '%s Unrotate cube %s' % (experiment_id, cs)

        lons, lats = np.meshgrid(lonr, latr)
        lons ,lats = unrotate_pole(lons,lats, cs.grid_north_pole_longitude, cs.grid_north_pole_latitude)

        lon=lons[0]
        lat=lats[:,0]

        csur=cs.ellipsoid

        for i, coord in enumerate (rot_cube.coords()):
            if coord.standard_name=='grid_latitude':
                lat_dim_coord_uwind = i
            if coord.standard_name=='grid_longitude':
                lon_dim_coord_uwind = i

        rot_cube.remove_coord('grid_latitude')
        rot_cube.remove_coord('grid_longitude')
        rot_cube.add_dim_coord(iris.coords.DimCoord(points=lat, standard_name='grid_latitude', units='degrees', coord_system=csur),lat_dim_coord_uwind )
        rot_cube.add_dim_coord(iris.coords.DimCoord(points=lon, standard_name='grid_longitude', units='degrees', coord_system=csur), lon_dim_coord_uwind)

    return rot_cube

for experiment_id in experiment_ids:

  expmin1 = experiment_id[:-1]
 
  fu = '/nfs/a90/eepdw/Data/EMBRACE/Mean_State/pp_files/%s/%s/%s.pp' % (expmin1, experiment_id, diag)

  print experiment_id
  sys.stdout.flush()

  try:
        #cube_names = ['%s' % cube_name_param, '%s' % cube_name_explicit]
        cube  = iris.load_cube(fu)
  except Exception,e:
      print e

  cube = unrotate_and_update_cube(cube)

  lat = cube.coord('grid_latitude').points
  lon = cube.coord('grid_longitude').points

  lons, lats = np.meshgrid(lon,lat)
  

  cube_regrid_data = scipy.interpolate.griddata((lats.flatten(), lons.flatten()),cube.data.flatten(), (lats_trmm, lons_trmm), method='linear')

  np.save("/nfs/a90/eepdw/Data/Rain_TRMM_regrid/rain_mean_regrid_onto_trmm_%s" %  experiment_id, cube_regrid_data)
  
