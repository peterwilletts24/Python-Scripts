import os, sys
import datetime

import iris
import iris.unit as unit
import iris.analysis.cartography

import numpy as np

import iris.analysis.geometry
from shapely.geometry import Polygon

from iris.coord_categorisation import add_categorised_coord

import imp
imp.load_source('UnrotateUpdateCube', '/nfs/see-fs-01_users/eepdw/python_scripts/Monsoon_Python_Scripts/modules/unrotate_and_update_pole.py')

from UnrotateUpdateCube import *

import pdb

diag = 'avg.5216'
cube_name_explicit='stratiform_rainfall_rate'
cube_name_param='convective_rainfall_rate'

pp_file_path='/nfs/a90/eepdw/Data/EMBRACE/'

experiment_ids = ['djzny', 'djznw', 'djznq', 'djzns', 'dkjxq', 'dklyu', 'dkmbq', 'dklwu', 'dklzq', 'dkbhu', 'djznu', 'dkhgu' ] # All 12

# Min and max lats lons from smallest model domain (dkbhu) - see spreadsheet

latmin=-6.79
latmax=33.038
lonmin=340
lonmax=379.98


lonmin_g=64.115
lonmax_g=101.866

lat_constraint=iris.Constraint(grid_latitude= lambda la: latmin <= la.point <= latmax)


#experiment_ids = ['dkhgu']
#experiment_ids = ['djzns', 'dklyu', 'dkmbq', 'dklwu', 'dklzq', 'dkbhu', 'djznu', 'dkhgu' ]
#experiment_ids = [ 'dklwu', 'dklzq', 'dklyu', 'dkmbq', 'dkbhu', 'djznu', 'dkhgu', 'djzns' ]
#experiment_ids = ['djznu', 'dkhgu' ] # High Res
#experiment_ids = ['djznw', 'djzny', 'djznq', 'dkjxq']
#experiment_ids = ['djznw', 'djzny', 'djznq', 'dkmbq', 'dklzq', 'dkjxq' ] # Params
# Load global LAM
# dtmindt = datetime.datetime(2011,8,19,0,0,0)
# dtmaxdt = datetime.datetime(2011,9,7,23,0,0)
# dtmin = unit.date2num(dtmindt, 'hours since 1970-01-01 00:00:00', unit.CALENDAR_STANDARD)
# dtmax = unit.date2num(dtmaxdt, 'hours since 1970-01-01 00:00:00', unit.CALENDAR_STANDARD)
# time_constraint = iris.Constraint(time= lambda t: dtmin <= t.point <= dtmax)

# Min and max lats lons from smallest model domain (dkbhu) - see spreadsheet


for experiment_id in experiment_ids:


    if (experiment_id=='djznw'):
        lon_constraint=iris.Constraint(grid_longitude= lambda lo: lonmin_g <= lo.point <= lonmax_g)
    else:
        lon_constraint=iris.Constraint(grid_longitude= lambda lo: lonmin <= lo.point <= lonmax) 


    expmin1 = experiment_id[:-1]

    fu = '%s%s/%s/%s.pp' % (pp_file_path, expmin1, experiment_id,  diag)

    flsm = '%s%s/%s/30.pp' % (pp_file_path, expmin1, experiment_id)
 
    print experiment_id
    sys.stdout.flush()

    try:
        #cube_names = ['%s' % cube_name_param, '%s' % cube_name_explicit]
        cubeconv  = iris.load_cube(fu,'%s' % cube_name_param  &  lat_constraint & lon_constraint)
        cubestrat  = iris.load_cube(fu,'%s' % cube_name_explicit &  lat_constraint & lon_constraint)
        cube=cubeconv+cubestrat
        cube.rename('total_precipitation_rate')
    except iris.exceptions.ConstraintMismatchError:        
        cube = iris.load_cube(fu, ('%s' % cube_name_explicit)   &  lat_constraint & lon_constraint)

     
        #cube.coord('grid_longitude').guess_bounds()
        #cube.coord('grid_latitude').guess_bounds()

        #cube= unrotate_pole_update_cube(cube)
   
    cube.coord('grid_longitude').guess_bounds()
    cube.coord('grid_latitude').guess_bounds()

    sys.stdout.flush()

    #pdb.set_trace()

# Calculate weights

        #l=iris.analysis.geometry.geometry_area_weights(cube, polygon) # Polygon weights

    # Load land/sea mask 

    lsm = iris.load_cube(flsm, ('land_binary_mask' ) &  lat_constraint & lon_constraint)

    print lsm
   
    sys.stdout.flush()

# For Sea and Land, mask area and calculate mean of each hour for sea/land and SAVE as numpy array

    for s in (['sea','land']):
      
        #pdb.set_trace()
        if s=='sea':
            w=1-lsm.data
        if s=='land':
            w=lsm.data

        w=np.repeat(w,cube.coord('time').points.shape[0]).reshape(cube.shape)

        coords = ('grid_latitude', 'grid_longitude')

        collapsed_cube = cube.collapsed(coords,
                                                   iris.analysis.MEAN,
                                                   weights=w)

     
        np.savez('%s%s/%s/%s_%s_rainfall_TimeVar_np_domain_constrain' % (pp_file_path, expmin1, experiment_id, diag, s), \
                 data=collapsed_cube.data.data, time_coords=collapsed_cube.coord('time').points)
     
    del lsm
  



 



