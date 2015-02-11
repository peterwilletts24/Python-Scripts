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

diag = '408_on_p_levs_mean_by_day'

pp_file_dir='/nfs/a90/eepdw/Data/EMBRACE/'

#experiment_ids = ['djzny', 'djznw', 'djznq', 'djzns', 'dkjxq', 'dklyu', 'dkmbq', 'dklwu', 'dklzq', 'dkbhu', 'djznu', 'dkhgu' ] # All 12

experiment_ids = ['dkmgw']

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

#latmin=-10
#latmax=5
#lonmin=64.115
#lonmax=80

# Monsoon Trough and Ganga Basin combined
polygon = Polygon(((73., 21.), (83., 16.), (87., 22.), (90.,22.), (90.,23.8), (83., 24.2), (76.3, 28.))) 

# Ganga Basin
#polygon = Polygon(((87., 22), (75., 27), (76.3, 30.), (83, 26.2), (90, 25.8), (90., 22)))

#lat_constraint=iris.Constraint(grid_latitude= lambda la: latmin <= la.point <= latmax)
#lon_constraint=iris.Constraint(grid_longitude= lambda lo: lonmin <= lo.point <= lonmax)

# fg = '%sdjzn/djznw/%s.pp' % (pp_file_path, diag)
# glob_load = iris.load_cube(fg, ('%s' % cube_name_param)  & time_constraint)   

## Get time points from global LAM to use as time constraint when loading other runs
#time_list = glob_load.coord('time').points
#glob_tc = iris.Constraint(time=time_list)

#del glob_load

for experiment_id in experiment_ids:

    expmin1 = experiment_id[:-1]

    pp_file_path = '%s_%s.pp' % (experiment_id, diag)

    floro = '%s%s/%s/33.pp' % (pp_file_dir, expmin1, experiment_id)
    fu = '%s%s/%s/%s' % (pp_file_dir, expmin1, experiment_id,  pp_file_path)

    #flsm = '%s%s/%s/30.pp' % (pp_file_path, expmin1, experiment_id)
 
    print experiment_id
    sys.stdout.flush()

    try:
      
        cube  = iris.load_cube(fu)
     
    
        #except iris.exceptions.ConstraintMismatchError:        
     
        cube.coord('grid_longitude').guess_bounds()
        cube.coord('grid_latitude').guess_bounds()

        cube= unrotate_pole_update_cube(cube)
   
        cube.coord('grid_longitude').guess_bounds()
        cube.coord('grid_latitude').guess_bounds()

        sys.stdout.flush()

# Calculate area weights to define whether in or out area polygon

        l=iris.analysis.geometry.geometry_area_weights(cube, polygon) # Polygon weights



        oro = iris.load_cube(floro, ('surface_altitude' ))
        oro = unrotate_pole_update_cube(oro)
        #oro=oro.extract(lat_constraint & lon_constraint)
        heights= oro.slices(['grid_latitude', 'grid_longitude']).next().data

 # Update cube data with masked array for area polygon

        #cube.data =  np.ma.masked_array(cube.data, l<=0.)
        cube.data = np.ma.masked_array(cube.data, (l<=0.) | (heights>=cube.data))
        #pdb.set_trace()

        for s in ([1]):
   
            coords = ('grid_latitude', 'grid_longitude')

            collapsed_cube = cube.collapsed(coords,
                                                   iris.analysis.MIN)
            #pdb.set_trace()

     #trnp =[collapsed_cube.data.data, collapsed_cube.coord('time').points]
     #if s == 0:
         # Areas of ocean
     #    print total_rainfall  
         #np.save('%s%s/%s/%s_sea_rainfall_diurnal_np_domain_constrain_lat_%s-%s_lon-%s-%s_monsoon_trough' % (pp_file_path, expmin1, experiment_id, diag, latmin, latmax, lonmin, lonmax), trnp)
         #np.save('%s%s/%s/%s_sea_rainfall_diurnal_np_domain_constrain_lat_%s-%s_lon-%s-%s_MASKED_ARRAY' % (pp_file_path, expmin1, experiment_id, diag, latmin, latmax, lonmin, lonmax), maskedcube)
            if s == 1:
          # Areas of land
                np.savez('%s%s/%s/TimeVar_%s_land_domain_constrain__and_oro_not_greater_than_data_monsoon_trough' \
                        % (pp_file_dir, expmin1, experiment_id, diag), \
                         data=collapsed_cube.data.data, time_coords=collapsed_cube.coord('time').points, pressures=collapsed_cube.coord('pressure').points) 
                 #iris.save(collapsed_cube, '%s%s/%s/TimeVar_%s_land_diurnal_domain_constrain_monsoon_trough.pp'  \
                #           % (pp_file_path, expmin1, experiment_id, diag)) 

    except Exception, e:
         print e
         pass



 



