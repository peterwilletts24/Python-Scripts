import cPickle as pickle
import numpy as np

import pdb

from datetime import datetime
from datetime import timedelta

def daterange( start_date, end_date ):
    if start_date <= end_date:
        for n in range( ( end_date - start_date ).days + 1 ):
            yield start_date + timedelta( n )
    else:
        for n in range( ( start_date - end_date ).days + 1 ):
            yield start_date - timedelta( n )
 
 


pcp_dom, longitude_dom, latitude_dom, time_dom, time_hour = pickle.load(open('/nfs/a90/eepdw/Data/Saved_data/TRMM/trmm_emb_time_update_large.p', 'rb'))

Time_Dom_Datetime = np.array([datetime.strptime('%s' % t, '%Y-%m-%d %H') for t in time_dom])
mo=np.array([t.month for t in Time_Dom_Datetime])
da=np.array([t.day for t in Time_Dom_Datetime])
# Calculate mean at each lat,lon position

longitude_domsingle = longitude_dom[1,:]
latitude_domsingle = latitude_dom[1,:]

hours=[' 0', ' 3', ' 6', ' 9', '12', '15', '18', '21']

dates=daterange(datetime.date(min(Time_Dom_Datetime)), datetime.date(max(Time_Dom_Datetime)))

mean_dom=np.empty((len(list(dates)), pcp_dom.shape[1],  pcp_dom.shape[2]))
#mean_dom=[]

#for i,h in enumerate(Time_Dom_Datetime():
dates=daterange(datetime.date(min(Time_Dom_Datetime)), datetime.date(max(Time_Dom_Datetime)))
for i, date in enumerate(dates):
    
    #Time_Dom_Datetime
    print date
    mean_dom[i]=np.mean(pcp_dom[np.where((mo==date.month) & (da==date.day))[0]],axis=0)

###############################

## Need to check time interval - TRMM data in 3 hourly intervals, given in mm/hr
############################

#time_interval=3

print mean_dom

dates_list=list(daterange(datetime.date(min(Time_Dom_Datetime)), datetime.date(max(Time_Dom_Datetime))))

#sum_dom = np.sum(pcp_dom*time_interval, axis=0)

pickle.dump([mean_dom, latitude_domsingle, longitude_domsingle, dates_list, hours], open('/nfs/a90/eepdw/Data/Saved_data/TRMM/trmm_emb_pcpmean_by_day.p', 'wb'))

#pickle.dump([sum_dom, latitude_domsingle, longitude_domsingle], open('/nfs/see-fs-01_users/eepdw/Saved_data/TRMM/trmm_emb_pcpsum_by_hour.p', 'wb'))
