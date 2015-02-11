import iris 

import imp

imp.load_source('UnrotateUpdateCube', '/nfs/see-fs-01_users/eepdw/python_scripts/Monsoon_Python_Scripts/modules/unrotate_and_update_pole.py')
from UnrotateUpdateCube import *



experiment_ids = ['djzny', 'djznq', 'djzns', 'dkjxq', 'dklyu', 'dkmbq', 'dklwu', 'dklzq', 'dkbhu', 'djznu', 'dkhgu' ] # All 12
pp_file_path = '/nfs/a90/eepdw/Data/EMBRACE/'

diag='30'
for experiment_id in experiment_ids:

    expmin1 = experiment_id[:-1]
    fu = '%s%s/%s/%s.pp' % (pp_file_path, expmin1, experiment_id,  diag)

    cube=iris.load(fu)
    
    print '%s, %s, %s, %s, %s' % (experiment_id, cube[0].coord('grid_latitude').points.min(), cube[0].coord('grid_latitude').points.max(), 
                                                                            cube[0].coord('grid_longitude').points.min(), cube[0].coord('grid_longitude').points.max())

for experiment_id in experiment_ids:

    expmin1 = experiment_id[:-1]
    fu = '%s%s/%s/%s.pp' % (pp_file_path, expmin1, experiment_id,  diag)

    cube=iris.load(fu)    

    cube[0].coord('grid_longitude').guess_bounds()
    cube[0].coord('grid_latitude').guess_bounds()

    cube= unrotate_pole_update_cube(cube[0])

    print '%s, %s, %s, %s, %s' % (experiment_id, cube.coord('grid_latitude').points.min(), cube.coord('grid_latitude').points.max(), 
                                                                            cube.coord('grid_longitude').points.min(), cube.coord('grid_longitude').points.max())

