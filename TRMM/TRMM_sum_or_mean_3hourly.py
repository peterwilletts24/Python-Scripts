import cPickle as pickle
import numpy as np

pcp_dom, longitude_dom, latitude_dom, time_dom, time_hour = pickle.load(open('/nfs/a90/eepdw/Data/Saved_data/TRMM/trmm_emb_time_update_large.p', 'rb'))

# Calculate mean at each lat,lon position

longitude_domsingle = longitude_dom[1,:]
latitude_domsingle = latitude_dom[1,:]

hours=[' 0', ' 3', ' 6', ' 9', '12', '15', '18', '21']
mean_dom=np.empty((len(hours), pcp_dom.shape[1],  pcp_dom.shape[2]))
for i,h in enumerate(hours):

    mean_dom[i]=np.mean(pcp_dom[time_hour==h],axis=0)

###############################

## Need to check time interval - TRMM data in 3 hourly intervals, given in mm/hr
############################

time_interval=3

print mean_dom

#sum_dom = np.sum(pcp_dom*time_interval, axis=0)

pickle.dump([mean_dom, latitude_domsingle, longitude_domsingle, time_dom, hours], open('/nfs/a90/eepdw/Data/Saved_data/TRMM/trmm_emb_pcpmean_by_hour.p', 'wb'))

#pickle.dump([sum_dom, latitude_domsingle, longitude_domsingle], open('/nfs/see-fs-01_users/eepdw/Saved_data/TRMM/trmm_emb_pcpsum_by_hour.p', 'wb'))
