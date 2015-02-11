"""
Load multiple pp diagnostic files, aggregate by year, day etc, calcualte mean, sum etc and save
"""

import os, sys
import datetime

import iris
import iris.unit as unit
from iris.coord_categorisation import add_categorised_coord
from iris.analysis.cartography import unrotate_pole
from iris.coords import DimCoord

import numpy as np

import re

import pdb

diag = '30201'
cube_names=['eastward_wind', 'northward_wind']
cube_names2=['x_wind','y_wind']

#pp_file_path='/projects/cascade/pwille/moose_retrievals/'
pp_file_path='/nfs/a90/eepdw/Data/EMBRACE/'

#experiment_ids = ['djzny', 'djznw', 'djznq', 'djzns', 'dkjxq', 'dklyu', 'dkmbq', 'dklwu', 'dklzq'] # All minus large 3
#experiment_ids = ['djzns', 'dkjxq', 'dklyu', 'dkmbq', 'dklwu', 'dklzq'] # All minus large 3

date_ranges=[[datetime.datetime(2011,8,18,0,0,0), datetime.datetime(2011,8,27,0,0,0)], 
                        [datetime.datetime(2011,8,28,0,0,0),datetime.datetime(2011,9,5,0,0,0)]]

experiment_ids = ['dklyu', 'dkmbq', 'dkmgw']

regrid_model='djznw'
regrid_model_min1=regrid_model[:-1]

#def add_hour_of_day(cube, coord, name='hour'):
 #   add_categorised_coord(cube, name, coord,
  #        lambda coord, x: coord.units.num2date(x).hour)

dtmindt = datetime.datetime(2011,8,19,0,0,0)
dtmaxdt = datetime.datetime(2011,9,7,23,0,0)
dtmin = unit.date2num(dtmindt, 'hours since 1970-01-01 00:00:00', unit.CALENDAR_STANDARD)
dtmax = unit.date2num(dtmaxdt, 'hours since 1970-01-01 00:00:00', unit.CALENDAR_STANDARD)
time_constraint = iris.Constraint(time= lambda t: dtmin <= t.point <= dtmax)

fr = '%s%s/%s/%s.pp' % (pp_file_path, regrid_model_min1, regrid_model, diag)
fg = '%sdjzn/djznw/%s.pp' % (pp_file_path, diag)
try:
    glob_load = iris.load_cube(fg, ('%s' % cube_names[0])  & time_constraint)
except iris.exceptions.ConstraintMismatchError:
    glob_load = iris.load_cube(fg, ('%s' % cube_names2[0])  & time_constraint)

## Get time points from global LAM to use as time constraint when loading other runs
time_list = glob_load.coord('time').points
# Some models have radiation diagnostics that are 10s offset from others so checking int values of time 
glob_tc = iris.Constraint(time= lambda t: int(t.point) in time_list.astype(int))
#glob_tc = iris.Constraint(time=time_list)

del glob_load

for experiment_id in experiment_ids:

 expmin1 = experiment_id[:-1]

 fu = '%s%s/%s/%s.pp' % (pp_file_path, expmin1, experiment_id, diag)
 save_name = '%s%s/%s/%s_%s_mean_by_date_range.pp' % (pp_file_path, expmin1, experiment_id, experiment_id, diag)

 print experiment_id
 sys.stdout.flush()

 try:
     os.remove(save_name)
 except OSError:
     print '%s NOT REMOVED' % save_name
     pass 



 
 try:
    for cube_name in cube_names:
     
         cube = iris.load_cube(fu, ('%s' % cube_name)  & glob_tc)

         #pdb.set_trace()

         time_coords = cube.coord('time')
 #add_hour_of_day(cube, time_coords)
         iris.coord_categorisation.add_day_of_year(cube, time_coords, name='day_of_year')

         day_bin = np.zeros(cube.coord('day_of_year').points.shape)
         for dr in date_ranges:
             day_of_year_range = [d.timetuple().tm_yday for d in dr]
             day_bin = np.where((cube.coord('day_of_year').points>min(day_of_year_range)) 
                                       & (cube.coord('day_of_year').points<max(day_of_year_range)), 
                                       str(day_of_year_range), day_bin)

         cube.add_aux_coord(iris.coords.AuxCoord(day_bin, var_name='day_bin_agg'), 0)

         for height_slice in cube.slices(['time', 'grid_latitude', 'grid_longitude']):
            mean = height_slice.aggregated_by('day_bin_agg', iris.analysis.MEAN)

            float_or_string = [float(x) if re.match(r'[-+]?(\d*\.?\d+|\d+\.)$', x) else x for x in mean.coord('day_bin_agg').points]
            idx_to_keep = np.where([type(fs)==np.string_ for fs in float_or_string])
            
            iris.save(mean[idx_to_keep], save_name, append=True)

 
                 
 except iris.exceptions.ConstraintMismatchError:        
         for cube_name in cube_names2:
                      cube = iris.load_cube(fu, ('%s' % cube_name)  & glob_tc)

                      #pdb.set_trace()

                      time_coords = cube.coord('time')
                      #add_hour_of_day(cube, time_coords)
                      iris.coord_categorisation.add_day_of_year(cube, time_coords, name='day_of_year')

                      day_bin = np.zeros(cube.coord('day_of_year').points.shape)
                      for dr in date_ranges:
                          day_of_year_range = [d.timetuple().tm_yday for d in dr]
                          day_bin = np.where((cube.coord('day_of_year').points>min(day_of_year_range)) 
                                             & (cube.coord('day_of_year').points<max(day_of_year_range)), 
                                             str(day_of_year_range), day_bin)

                      cube.add_aux_coord(iris.coords.AuxCoord(day_bin, var_name='day_bin_agg'), 0)

                      for height_slice in cube.slices(['time', 'grid_latitude', 'grid_longitude']):
                          mean = height_slice.aggregated_by('day_bin_agg', iris.analysis.MEAN)

                          float_or_string = [float(x) if re.match(r'[-+]?(\d*\.?\d+|\d+\.)$', x) else x for x in mean.coord('day_bin_agg').points]
                          idx_to_keep = np.where([type(fs)==np.string_ for fs in float_or_string])
            
                          iris.save(mean[idx_to_keep], save_name, append=True)

 np.save( '%s%s/%s/%s_mean_by_date_range_anc_file' \
               % (pp_file_path, expmin1, experiment_id, experiment_id), date_ranges)



