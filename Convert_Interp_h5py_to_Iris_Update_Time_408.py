from iris.coord_systems import RotatedGeogCS, GeogCS
from iris.analysis.cartography import unrotate_pole
from iris.coords import DimCoord
import iris

import h5py

import numpy as np

import pdb

#experiment_ids = ['dklyu', 'dkmbq']
experiment_ids = ['djzny', 'djznq', 'djzns', 'dkjxq', 'dklyu', 'dkmbq', 'dklwu', 'dklzq', 'dkbhu', 'djznu', 'dkhgu' ] # All 12
#experiment_ids = ['dkjxq','dklwu', 'dklzq', 'dkbhu', 'djznu', 'dkhgu' ] 
#experiment_ids = [ 'dklyu'] 

diag = '_408_on_p_levs'

for experiment_id in experiment_ids:

    expmin1 = experiment_id[:-1]

    fname = '/nfs/a90/eepdw/Data/EMBRACE/%s/%s/30201.pp' %(expmin1, experiment_id)

    time_coord_cube = iris.load(fname)

    time_cube_loop_points=[]

    for t, time_cube in enumerate(time_coord_cube[0].slices(['pressure', 'grid_latitude', 'grid_longitude'])):
      
        time_cube_loop_points.append(time_cube.coord('time').points[0])

    time_coord_argsort = np.argsort(np.array(time_cube_loop_points))

    fu = '/nfs/a90/eepdw/Data/EMBRACE/%s/%s/%s%s.pp' % (expmin1, experiment_id, experiment_id, diag)

    cube = iris.load_cube(fu)

    cube.remove_coord('time')

    #pdb.set_trace()

    cube.add_dim_coord(DimCoord(points=np.array(time_cube_loop_points), 
                    long_name='time', standard_name='time', units=time_coord_cube[0].coord('time').units),0)
                        
    print experiment_id
    fu_2 = '/nfs/a90/eepdw/Data/EMBRACE/%s/%s/%s%s_updated.pp' % (expmin1, experiment_id, experiment_id, diag)    
    iris.save(cube, '%s' % fu_2, field_coords=('grid_latitude','grid_longitude')) 


    #pdb.set_trace

    

  

