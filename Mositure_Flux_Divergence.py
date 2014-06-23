import os, sys
sys.path.append('/nfs/see-fs-01_users/eepdw/python_scripts/modules')

from update_pp_cube_coords import update_coordsimport iris

import iris.analysis.cartography

import h5py

import numpy as np

import pdb

import scipy

#Load specific humidity and wind

#/nfs/a90/eepdw/Data/EMBRACE/Pressure_level\_means/sp_hum_pressure_levels_interp_djzns_mean_masked/

#experiment_ids = ['djznw', 'djznq', 'djzny', 'djzns', 'dkmbq', 'dklyu', 'dklwu', 'dklzq']
experiment_ids = ['dkmbq', 'dklyu']
p_levels = [1000, 950, 925, 850, 700, 500, 400, 300, 250, 200, 150, 100, 70, 50, 30, 20, 10]

for experiment_id in experiment_ids:
    
   expmin1 = experiment_id[:-1]

   fname_h = '/nfs/a90/eepdw/Data/EMBRACE/On_Heights_Interpolation_Data/sp_hum_pressure_levels_interp_%s' % (experiment_id)

   with h5py.File(fname_h, 'r') as i:
        
#  q = i['%s' % 'mean'][. . .]
    q = i['%s' % 'sh_on_p'][. . .]
   pdb.set_trace()


   f_oro =  '/nfs/a90/eepdw/Data/EMBRACE/Mean_State/pp_files/%s/%s/33.pp' % (expmin1, experiment_id)
   oro = iris.load_cube(f_oro)
   oro,lats,lons = update_coords(oro)

   fu = '/nfs/a90/eepdw/Data/EMBRACE/Mean_State/pp_files/%s/%s/30201.pp' % (expmin1, experiment_id)
    
   u_wind,v_wind = iris.load(fu)

   u_wind,lats_w,lons_w = update_coords(u_wind)
   v_wind,lats_w,lons_w = update_coords(v_wind)

   print u_wind
   print u_wind.coord('pressure')

   qu_div = np.empty((u_wind.shape[1], u_wind.shape[2], u_wind.shape[0]))
   qv_div = np.empty((u_wind.shape[1], u_wind.shape[2], u_wind.shape[0]))
 
   fl_la_lo = (lats.flatten(),lons.flatten())

   p_lev_delta = np.diff( np.append(u_wind.coord('pressure').points, u_wind.coord('pressure').points[-1]))
   
   for p, pressure_cube in enumerate(u_wind.slices(['grid_latitude', 'grid_longitude'])):
    
    print p
    s = np.searchsorted(p_levels[::-1], p)
    #sc =  np.searchsorted(p_levs, p)

    q_slice = q[:,:,-(s+1)]
    q_interp = scipy.interpolate.griddata(fl_la_lo, q_slice.flatten(), (lats_w, lons_w), method='linear')

    #pdb.set_trace()
    qu_div[:,:,pn] = (u_wind.coord('pressure').points[pn]*q_interp)/9.81
    qv_div[:,:,pn] = (v_wind.coord('pressure').points[pn]*q_interp)/9.81

   Qu_div = np.sum(qu_div*p_lev_delta, axis=-1)
   Qv_div = np.sum(qv_div*p_lev_delta, axis=-1)

   


    
    

