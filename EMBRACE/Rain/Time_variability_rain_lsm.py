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

#experiment_ids = ['djzny', 'djznw', 'djznq', 'djzns', 'dkjxq', 'dklyu', 'dkmbq', 'dklwu', 'dklzq', 'dkbhu', 'djznu', 'dkhgu' ] # All 12
experiment_ids = ['dkmgw']
#experiment_ids = ['dkhgu']
#experiment_ids = ['djzns', 'dklyu', 'dkmbq', 'dklwu', 'dklzq', 'dkbhu', 'djznu', 'dkhgu' ]
#experiment_ids = [ 'dklwu', 'dklzq', 'dklyu', 'dkmbq', 'dkbhu', 'djznu', 'dkhgu', 'djzns' ]
#experiment_ids = ['djznu', 'dkhgu' ] # High Res
#experiment_ids = ['djznw', 'djzny', 'djznq', 'dkjxq']
#experiment_ids = ['djznw', 'djzny', 'djznq', 'dkmbq', 'dklzq', 'dkjxq' ] # Params
# Load global LAM
dtmindt = datetime.datetime(2011,8,19,0,0,0)
dtmaxdt = datetime.datetime(2011,9,7,23,0,0)
dtmin = unit.date2num(dtmindt, 'hours since 1970-01-01 00:00:00', unit.CALENDAR_STANDARD)
dtmax = unit.date2num(dtmaxdt, 'hours since 1970-01-01 00:00:00', unit.CALENDAR_STANDARD)
time_constraint = iris.Constraint(time= lambda t: dtmin <= t.point <= dtmax)

# Min and max lats lons from smallest model domain (dkbhu) - see spreadsheet

latmin=-20
latmax=19.98
lonmin=340
lonmax=379.98

latmin_g=-6.79043923611
latmax_g=33.037779928
lonmin_g=64.115
lonmax_g=101.866

lat_constraint=iris.Constraint(grid_latitude= lambda la: latmin <= la.point <= latmax)
#lon_constraint=iris.Constraint(grid_longitude= lambda lo: lonmin <= lo.point <= lonmax)

fg = '%sdjzn/djznw/%s.pp' % (pp_file_path, diag)
glob_load = iris.load_cube(fg, ('%s' % cube_name_param)  & time_constraint)   

## Get time points from global LAM to use as time constraint when loading other runs
time_list = glob_load.coord('time').points
glob_tc = iris.Constraint(time=time_list)

del glob_load

#polygon = Polygon(((73., 21.), (83., 16.), (87., 22.), (75., 27.)))

for experiment_id in experiment_ids:

    if (experiment_id=='djznw'):
        lon_constraint=iris.Constraint(grid_longitude= lambda lo: lonmin_g <= lo.point <= lonmax_g)
        lat_constraint=iris.Constraint(grid_latitude= lambda la: latmin_g <= la.point <= latmax_g)
    else:
        lat_constraint=iris.Constraint(grid_latitude= lambda la: latmin <= la.point <= latmax)
        lon_constraint=iris.Constraint(grid_longitude= lambda lo: lonmin <= lo.point <= lonmax) 

    expmin1 = experiment_id[:-1]

    fu = '%s%s/%s/%s.pp' % (pp_file_path, expmin1, experiment_id,  diag)

    flsm = '%s%s/%s/30.pp' % (pp_file_path, expmin1, experiment_id)
 
    print experiment_id
    sys.stdout.flush()

    #pdb.set_trace()

    try:
        #cube_names = ['%s' % cube_name_param, '%s' % cube_name_explicit]
        cubeconv  = iris.load_cube(fu,('%s' % cube_name_param) & glob_tc & lat_constraint & lon_constraint)
        cubestrat  = iris.load_cube(fu,('%s' % cube_name_explicit) & glob_tc & lat_constraint & lon_constraint)
        cube=cubeconv+cubestrat
        cube.rename('total_precipitation_rate')
    except iris.exceptions.ConstraintMismatchError:        
        cube = iris.load_cube(fu, ('%s' % cube_name_explicit)  & glob_tc & lat_constraint & lon_constraint)
    
        #except iris.exceptions.ConstraintMismatchError:        
     
        cube.coord('grid_longitude').guess_bounds()
        cube.coord('grid_latitude').guess_bounds()

        cube= unrotate_pole_update_cube(cube)
   
        cube.coord('grid_longitude').guess_bounds()
        cube.coord('grid_latitude').guess_bounds()

        sys.stdout.flush()

# Calculate weights

        #l=iris.analysis.geometry.geometry_area_weights(cube, polygon) # Polygon weights
# Load land/sea mask 

    lsm = iris.load_cube(flsm, ('land_binary_mask' ) &  lat_constraint & lon_constraint)

# For Sea and Land, mask area and calculate mean of each hour for sea/land and SAVE as numpy array
    time_coords = cube.coord('time')
    iris.coord_categorisation.add_day_of_year(cube, time_coords, name='day_of_year')

    for s in ([0,1]):
        #pdb.set_trace()
        sm = np.where(lsm.data==s, cube.data, np.nan)
        maskedcube = np.ma.masked_array(sm, sm==np.nan)
        masked_iris_cube = cube.copy(data=maskedcube)

        mean = masked_iris_cube.aggregated_by('day_of_year', iris.analysis.MEAN)
        mean2=mean.collapsed(('grid_latitude', 'grid_longitude'), iris.analysis.MEAN)

        #trnp =[mean2.data.data, masked_iris_cube.coord('day_of_year').points]
        if s == 0:
             # Areas of ocean
             #print total_rainfall
            np.savez('%s%s/%s/%s_sea_rainfall_np_domain_constrain_daily' % (pp_file_path, expmin1, experiment_id, diag), 
                     data=mean2.data.data, time_coords=mean2.coord('time').points)
        
        if s == 1:
            # Areas of land
            np.savez('%s%s/%s/%s_land_rainfall_np_domain_constrain_daily' % (pp_file_path, expmin1, experiment_id, diag), 
                    data=mean2.data.data, time_coords=mean2.coord('time').points)      
     

   
