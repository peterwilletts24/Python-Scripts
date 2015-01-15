import cPickle as pickle
import numpy as np

import pdb

pcp_dom, longitude_dom, latitude_dom, time_dom, time_hour = pickle.load(open('/nfs/a90/eepdw/Data/Saved_data/TRMM/trmm_emb_time_update.p', 'rb'))

# Calculate mean at each lat,lon position

longitude_domsingle = longitude_dom[1,:]
latitude_domsingle = latitude_dom[1,:]

# lat_min=-6.79
# lat_max=33.038
# lon_min=64.115
# lon_max=101.866

# la_index = np.where((latitude_dom<=lat_max) & (latitude_dom >= lat_min))
# lo_index = np.where((longitude_dom<=lon_max) & (longitude_dom >= lon_min))

# pdb.set_trace()

# la_i_max = np.max(la_index)
# la_i_min = np.min(la_index)
# lo_i_max = np.max(lo_index)
# lo_i_min = np.min(lo_index)
   
# lat_amounts=la_i_max-la_i_min
# lon_amounts=lo_i_max-lo_i_min



mean_dom = np.mean(pcp_dom,axis=tuple(range(1, pcp_dom.ndim)))

#sum_dom = np.sum(pcp_dom,axis=tuple(range(1, pcp_dom.ndim)))

pickle.dump([mean_dom, latitude_domsingle, longitude_domsingle, time_dom], open('/nfs/a90/eepdw/Data/Saved_data/TRMM/trmm_emb_pcpmean_time_va_dkbhu_domain.p', 'wb'))

#pickle.dump([sum_dom, latitude_domsingle, longitude_domsingle], open('/nfs/see-fs-01_users/eepdw/Saved_data/TRMM/trmm_emb_pcpsum_by_hour.p', 'wb'))
