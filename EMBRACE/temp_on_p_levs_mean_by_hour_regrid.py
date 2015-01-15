"""
Load multiple pp diagnostic files, aggregate by year, day etc, calcualte mean, sum etc and save
"""

import os, sys
import datetime

import iris
import iris.unit as unit
from iris.coord_categorisation import add_categorised_coord
from iris.coords import DimCoord

import numpy as np

from iris.analysis.cartography import unrotate_pole

import pdb

diag = '_temp_on_p_levs'
#cube_name='surface_upward_latent_heat_flux'
dm='temp'

pp_file_path='/nfs/a90/eepdw/Data/EMBRACE/'

#experiment_ids = ['djznw', 'djzny', 'djznq', 'djzns', 'dkjxq', 'dklyu', 'dkmbq', 'dklwu', 'dklzq'] # All minus large 3
#experiment_ids = ['djzns', 'dkjxq', 'dklyu', 'dkmbq', 'dklwu', 'dklzq'] # All minus large 3
#experiment_ids = ['dklyu', 'dkmbq', 'dklwu', 'dklzq', 'dkjxq'] # All minus large 3
experiment_ids = ['djzny', 'djznq', 'dklyu', 'dkmbq', 'dklwu', 'dklzq', 'dkjxq'] # All minus large 3

#experiment_ids = ['djznw', 'djzny', 'djznq', 'dkjxq', 'dkmbq', 'dklzq']
#experiment_ids = ['dklyu', 'dkmbq']
regrid_model='djznw'
regrid_model_min1=regrid_model[:-1]

def add_hour_of_day(cube, coord, name='hour'):
    add_categorised_coord(cube, name, coord,
          lambda coord, x: coord.units.num2date(x).hour)

dtmindt = datetime.datetime(2011,8,19,0,0,0)
dtmaxdt = datetime.datetime(2011,9,7,23,0,0)
dtmin = unit.date2num(dtmindt, 'hours since 1970-01-01 00:00:00', unit.CALENDAR_STANDARD)
dtmax = unit.date2num(dtmaxdt, 'hours since 1970-01-01 00:00:00', unit.CALENDAR_STANDARD)
time_constraint = iris.Constraint(time= lambda t: dtmin <= t.point <= dtmax)

fr = '%s%s/%s/%s%s.pp' % (pp_file_path, regrid_model_min1, regrid_model, regrid_model, diag)
fg = '%sdjzn/djznw/%s%s.pp' % (pp_file_path, regrid_model, diag)
#glob_load = iris.load_cube(fg, ('%s' % cube_name)  & time_constraint)
glob_load = iris.load_cube(fg,  time_constraint)
## Get time points from global LAM to use as time constraint when loading other runs
time_list = glob_load.coord('time').points
del glob_load

# Some models have radiation diagnostics that are 10s offset from others so checking int values of time 
#glob_tc = iris.Constraint(time= lambda t: int(t.point) in time_list.astype(int))
glob_tc = iris.Constraint(time=time_list)

regrid_cube_init = iris.load_cube(fr, glob_tc)



for experiment_id in experiment_ids:

 expmin1 = experiment_id[:-1]

 fu = '%s%s/%s/%s%s.pp' % (pp_file_path, expmin1, experiment_id, experiment_id, diag)

 print experiment_id
 sys.stdout.flush()

      
 cube = iris.load_cube(fu,   glob_tc)
 #cube = iris.load_cube(fu)       
 cs = cube.coord_system('CoordSystem')

 lons, lats = np.meshgrid(cube.coord('grid_longitude').points, cube.coord('grid_latitude').points)

 unrot_lons, unrot_lats = unrotate_pole(lons, 
                                        lats,
                                        cs.grid_north_pole_longitude,
                                        cs.grid_north_pole_latitude)
                                         
 latmin = unrot_lats.min()
 latmax = unrot_lats.max()
 lonmin = unrot_lons.min()
 lonmax = unrot_lons.max()

 lat_constraint=iris.Constraint(grid_latitude= lambda la: latmin <= la.point <= latmax)
 lon_constraint=iris.Constraint(grid_longitude= lambda lo: lonmin <= lo.point <= lonmax)  
#time_constraint=
 regrid_cube = regrid_cube_init.extract(lat_constraint & lon_constraint)   

 #pdb.set_trace()

 mean_list=[]
 for pressure_cube in (cube.slices(['time', 'grid_latitude', 'grid_longitude'])):

     pressure_cube = iris.analysis.interpolate.regrid(pressure_cube, regrid_cube, mode='bilinear')

     time_coords = pressure_cube.coord('time')
     add_hour_of_day(pressure_cube, time_coords)

     pc_time_merge = pressure_cube.aggregated_by('hour', iris.analysis.MEAN)
     pc_time_merge.add_dim_coord(DimCoord(points=pc_time_merge.coords('time')[0].bounds[:,0].flatten(),\
                                        long_name='time', standard_name='time',units=pc_time_merge.coords('time')[0].units),0)


     mean_list.extend(iris.cube.CubeList([pc_time_merge]))
     #mean_list.sort()

 #pdb.set_trace()

 mean = iris.cube.CubeList(mean_list).merge_cube() 

 mean.coord('grid_latitude').guess_bounds()
 mean.coord('grid_longitude').guess_bounds()
 #mean.coords('time')[0].guess_bounds()
 #mean.add_dim_coord(DimCoord(points=mean.coords('time')[0].bounds[:,0].flatten(), long_name='time', standard_name='time',\
                        #    units=mean.coords('time')[0].units),0)

 if dm=='temp':
        mean.rename('potential_temperature')
        mean.units=iris.unit.Unit('K')
 if dm=='sp_hum':
        mean.rename('specific_humidity')
        mean.units=iris.unit.Unit('kg kg-1')
 if dm=='408':
        mean.rename('height')
        mean.units=iris.unit.Unit('m')

 iris.save((mean),'%s%s/%s/%s%s_mean_by_hour_regrid.pp' % (pp_file_path, expmin1, experiment_id, experiment_id, diag),\
                            field_coords=('grid_latitude','grid_longitude'))

 try:
     iris.load_cube('%s%s/%s/%s%s_mean_by_hour_regrid.pp' % (pp_file_path, expmin1, experiment_id, experiment_id, diag))
 except iris.exceptions.ConstraintMismatchError:

     #pdb.set_trace()

     save_as_cube=iris.cube.Cube(mean.data)
     save_as_cube.add_dim_coord(mean.coord('pressure'),0)
     save_as_cube.add_dim_coord(mean.coord('grid_latitude'),2)
     save_as_cube.add_dim_coord(mean.coord('grid_longitude'),3)      
     save_as_cube.add_dim_coord(mean.coords('time')[0],1)

     if dm=='temp':
        save_as_cube.rename('potential_temperature')
        save_as_cube.units=iris.unit.Unit('K')
     if dm=='sp_hum':
        save_as_cube.rename('specific_humidity')
        save_as_cube.units=iris.unit.Unit('kg kg-1')
     if dm=='408':
        save_as_cube.rename('height')
        save_as_cube.units=iris.unit.Unit('m')

     iris.save((save_as_cube),'%s%s/%s/%s%s_mean_by_hour_regrid.pp'  % (pp_file_path, expmin1, experiment_id, experiment_id, diag)\
                                             ,field_coords=('grid_latitude','grid_longitude'))
     del save_as_cube











