import pdb
import numpy as np

import iris

experiment_ids = ['djzny', 'dklyu', 'dkmbq', 'djznq']

pp_file_path='/nfs/a90/eepdw/Data/EMBRACE/'

pp_file = 'tcwv'
#for ls in ['land', 'sea']:

for experiment_id in experiment_ids:
        
        print experiment_id
    
        expmin1 = experiment_id[:-1]

        #data = np.load('%s%s/%s/%s_%s_%s_hourly_mean.npz' % (pp_file_path, expmin1, experiment_id, pp_file, experiment_id, ls))

        data = iris.load_cube('%s%s/%s/%s_%s.pp' % (pp_file_path, expmin1, experiment_id, experiment_id, pp_file))
        #pdb.set_trace()

        data_fix = data/100000

        iris.save(data_fix, '%s%s/%s/%s_%s.pp' % (pp_file_path, expmin1, experiment_id, experiment_id, pp_file))

         #np.savez('%s%s/%s/%s_%s_%s_hourly_mean' 
          #       % (pp_file_path, expmin1, experiment_id, pp_file, experiment_id, ls), 
         #    means_hourly=data_fix, hours=data['hours'])

