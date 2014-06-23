"""

Load pp which has same diagnostic split into smaller cubes ( append save for alrge datasets) and saves as one


"""

import os, sys

import iris
import iris.experimental
import iris.coords

import numpy as np

def main():

 pp_file = '1201_mean'

 pp_file_dir ='/nfs/a90/eepdw/Data/EMBRACE/Mean_State/pp_files/'
 experiment_ids = ['dkhgu', 'dkbhu', 'djznu'] 
 
 #experiment_ids = ['djznu' ] 
 for experiment_id in experiment_ids:
 
  expmin1 = experiment_id[:-1]
  pfile = '%s%s/%s/%s.pp' % (pp_file_dir, expmin1, experiment_id, pp_file)
  sfile = '%s%s/%s/%s_cat.pp' % (pp_file_dir, expmin1, experiment_id, pp_file)
  ofile = '%s%s/%s/33.pp' % (pp_file_dir, expmin1, experiment_id)

  oc = iris.load_cube(ofile)
  
  plist = iris.load(pfile)
  
  #pcubef=np.empty(oc.shape, np.float32)
  
  latitude = iris.coords.DimCoord(oc.coord('grid_latitude').points, standard_name='grid_latitude',
                    units='degrees')
  longitude = iris.coords.DimCoord(oc.coord('grid_longitude').points, standard_name='grid_longitude',
                     units='degrees')
  

  for pcube in plist:

      armin = np.searchsorted(oc.coord('grid_longitude').points, min(pcube.coord('grid_longitude').points))
      armax = np.searchsorted(oc.coord('grid_longitude').points, max(pcube.coord('grid_longitude').points))

      #pcubef[:,armin:armax+1] = pcube.data
      pcubef[:,armin:armax+1] = pcube.data               
  
      pc = iris.cube.Cube(pcubef, standard_name=pcube.standard_name, units =pcube.units,dim_coords_and_dims=[(latitude, 0), (longitude, 1)])
      print pc

      iris.save((pcube), '%s%s/%s/%s_cat.pp' % (pp_file_dir, expmin1, experiment_id, pp_file), append=True) 
if __name__ == '__main__':
   main()    
