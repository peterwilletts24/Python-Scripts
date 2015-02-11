from iris.coord_systems import RotatedGeogCS, GeogCS
from iris.analysis.cartography import unrotate_pole
from iris.coords import DimCoord
import iris

import h5py

import numpy as np

import pdb

#experiment_ids = ['dklyu', 'dkmbq']
#experiment_ids = ['djznw', 'djzny', 'djznq', 'djzns', 'dkjxq', 'dklyu', 'dkmbq', 'dklwu', 'dklzq', 'dkbhu', 'djznu', 'dkhgu' ] # All 12
#experiment_ids = ['djznw', 'djzny', 'djznq', 'djzns', 'dklyu', 'dkmbq', 'dklwu', 'dklzq', 'dkbhu' ] # All small-ish res minus dkjxq
#experiment_ids = ['dkjxq','dklwu', 'dklzq', 'dkbhu', 'djznu', 'dkhgu' ] 
experiment_ids = [ 'dkbhu'] 
#experiment_ids = ['djznw']
#experiment_ids = ['dklyu', 'dkmbq', 'dklwu', 'dklzq']
#experiment_ids = ['dkmgw']

data_to_save = ['408']
dset = ['interps']

p_levels = [1000, 950, 925, 850, 700, 500, 400, 300, 250, 200, 150, 100, 70, 50, 30, 20, 10]

import linecache, sys

def PrintException():
    exc_type, exc_obj, tb = sys.exc_info()
    f = tb.tb_frame
    lineno = tb.tb_lineno
    filename = f.f_code.co_filename
    linecache.checkcache(filename)
    line = linecache.getline(filename, lineno, f.f_globals)
    print 'EXCEPTION IN ({}, LINE {} "{}"): {}'.format(filename, lineno, line.strip(), exc_obj)

for experiment_id in experiment_ids:
 try:
    #pdb.set_trace()

    expmin1 = experiment_id[:-1]

    info_cube = iris.load('/nfs/a90/eepdw/Data/EMBRACE/%s/%s/30201.pp' % (expmin1, experiment_id))[0]
    info_cube2 = iris.load('/nfs/a90/eepdw/Data/EMBRACE/%s/%s/33.pp' % (expmin1, experiment_id))[0]
    for a, dm in enumerate(data_to_save):
    
        fname = '/nfs/a90/eepdw/Data/EMBRACE/On_Heights_Interpolation_Data/%s_pressure_levels_interp_pressure_%s' % (dm,experiment_id)
        ds = dset[a]

        time_cube_loop_points=[]

        for t, time_cube in enumerate(info_cube.slices(['pressure', 'grid_latitude', 'grid_longitude'])):
      
            time_cube_loop_points.append(time_cube.coord('time').points[0])

        time_coord_argsort = np.argsort(np.array(time_cube_loop_points))

        time_cube_loop_points = np.array(time_cube_loop_points)[time_coord_argsort]

        #pdb.set_trace()

        with h5py.File(fname, 'r') as f:

            save_as_cube=iris.cube.Cube(f['%s' % ds][:][time_coord_argsort], standard_name=None, long_name=None, var_name=None, units=None, attributes=None, cell_methods=None, dim_coords_and_dims=None, aux_coords_and_dims=None, aux_factories=None)

# for i, coord in enumerate (info_cube.coords()):
#             if coord.standard_name=='grid_latitude':
#                 lat_dim_coord_cube = i
#             if coord.standard_name=='grid_longitude':
#                 lon_dim_coord_cube = i
#             if coord.standard_name=='time':
#                 time_dim_coord_cube = i
#             if (coord.standard_name=='pressure') | (coord.long_name=='pressure'):
#                 pres_dim_coord_cube = i
        #pdb.set_trace()

        save_as_cube.add_dim_coord(DimCoord(points=p_levels, long_name='pressure', units='hPa'),3)
        save_as_cube.add_dim_coord((info_cube2.coord('grid_latitude')),1)
        save_as_cube.add_dim_coord((info_cube2.coord('grid_longitude')),2)
        save_as_cube.add_dim_coord(DimCoord(points=np.array(time_cube_loop_points), long_name='time', standard_name='time', units=info_cube.coord('time').units),0)

        if dm=='temp':
            save_as_cube.rename('potential_temperature')
            save_as_cube.units=iris.unit.Unit('K')
        if dm=='sp_hum':
            save_as_cube.rename('specific_humidity')
            save_as_cube.units=iris.unit.Unit('kg kg-1')
        if dm=='408':
            save_as_cube.rename('height')
            save_as_cube.units=iris.unit.Unit('m')
        iris.save(save_as_cube, '/nfs/a90/eepdw/Data/EMBRACE/%s/%s/%s_%s_on_p_levs.pp' % (expmin1, experiment_id, experiment_id, dm), field_coords=('grid_latitude','grid_longitude')) 

        
 except Exception:
  PrintException()
