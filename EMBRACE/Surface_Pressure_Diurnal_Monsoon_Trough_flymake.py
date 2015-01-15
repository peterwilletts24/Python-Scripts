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

import linecache

def PrintException():
    exc_type, exc_obj, tb = sys.exc_info()
    f = tb.tb_frame
    lineno = tb.tb_lineno
    filename = f.f_code.co_filename
    linecache.checkcache(filename)
    line = linecache.getline(filename, lineno, f.f_globals)
    print 'EXCEPTION IN ({}, LINE {} "{}"): {}'.format(filename, lineno, line.strip(), exc_obj)


diag = '_408_on_p_levs_mean_by_hour'

pp_file_path='/nfs/a90/eepdw/Data/EMBRACE/'

experiment_ids = ['djznw', 'djzny', 'djznq', 'djzns', 'dkjxq', 'dklyu', 'dkmbq', 'dklwu', 'dklzq', 'dkbhu', 'djznu', 'dkhgu' ] # All 12
#experiment_ids = ['dklyu', 'dkmbq', 'dklwu', 'dklzq'] 

#experiment_ids = ['djzns']

# Min and max lats lons from smallest model domain (dkbhu) - see spreadsheet

#latmin=-10
#latmax=5
#lonmin=64.115
#lonmax=80

# Monsoon Trough

polygon = Polygon(((73., 21.), (83., 16.), (87., 22.), (75., 27.))) 
        
#lat_constraint=iris.Constraint(grid_latitude= lambda la: latmin <= la.point <= latmax)
#lon_constraint=iris.Constraint(grid_longitude= lambda lo: lonmin <= lo.point <= lonmax)


for experiment_id in experiment_ids:

    expmin1 = experiment_id[:-1]

    fu = '%s%s/%s/%s%s.pp' % (pp_file_path, expmin1, experiment_id, experiment_id, diag)

    flsm = '%s%s/%s/30.pp' % (pp_file_path, expmin1, experiment_id)
 
    print experiment_id
    sys.stdout.flush()

    try:
        #cube_names = ['%s' % cube_name_param, '%s' % cube_name_explicit]
        cube  = iris.load_cube(fu)

        cube= unrotate_pole_update_cube(cube)

        cube.coord('grid_longitude').guess_bounds()
        cube.coord('grid_latitude').guess_bounds()
     
# Calculate weights

        l=iris.analysis.geometry.geometry_area_weights(cube, polygon)

# For Sea and Land, mask area and calculate mean of each hour for sea/land and SAVE as numpy array
 
        coords = ('grid_latitude', 'grid_longitude')

   
        collapsed_cube = cube.collapsed(coords,
                                            iris.analysis.MEAN,
                                            weights=l) 
    
        np.savez('%s%s/%s/%s%s_diurnal_monsoon_trough' % (pp_file_path, expmin1, experiment_id, experiment_id, diag), \
                data=collapsed_cube.data.data, time=collapsed_cube.coord('time').points, pressure=collapsed_cube.coord('pressure').points)
     
    except iris.exceptions.ConstraintMismatchError:  
        PrintException()






