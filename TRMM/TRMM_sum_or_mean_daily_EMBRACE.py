import cPickle as pickle
import numpy as np

import pdb

from collections import defaultdict
from netCDF4 import Dataset
from scipy.interpolate import griddata

from datetime import datetime
from datetime import timedelta

def daterange( start_date, end_date ):
    if start_date <= end_date:
        for n in range( ( end_date - start_date ).days + 1 ):
            yield start_date + timedelta( n )
    else:
        for n in range( ( start_date - end_date ).days + 1 ):
            yield start_date - timedelta( n )
 
latmin=-6.79
latmax=33.04
lonmin=64.12
lonmax=101.87


pcp_dom, longitude_dom, latitude_dom, time_dom, time_hour = pickle.load(open('/nfs/a90/eepdw/Data/Saved_data/TRMM/trmm_emb_time_update_large.p', 'rb'))

Time_Dom_Datetime = np.array([datetime.strptime('%s' % t, '%Y-%m-%d %H') for t in time_dom])
mo=np.array([t.month for t in Time_Dom_Datetime])
da=np.array([t.day for t in Time_Dom_Datetime])
# Calculate mean at each lat,lon position

longitude_domsingle = longitude_dom[1,:]
latitude_domsingle = latitude_dom[1,:]

hours=[' 0', ' 3', ' 6', ' 9', '12', '15', '18', '21']

# Load land sea mask.  TRMM land sea mask is in % of water coverage so 100% is all water

nc = Dataset('/nfs/a90/eepdw/Data/Observations/Satellite/TRMM/TMPA_mask.nc')

#  Regrid lsm to data grid (offset b 0.125 degrees

lsm_lons, lsm_lats = np.meshgrid(nc.variables['lon'][:],nc.variables['lat'][:])
lons_data, lats_data = np.meshgrid(longitude_dom[0], latitude_dom[0])
lsm_regrid = griddata((lsm_lats.flatten(), lsm_lons.flatten()), nc.variables['landseamask'][:].flatten(), (lats_data,lons_data), method='linear')

#  Find points that are within defined area

lat_max_idx = np.max(np.where((latmin<=latitude_dom[0]) & (latitude_dom[0]<=latmax)))
lat_min_idx = np.min(np.where((latmin<=latitude_dom[0]) & (latitude_dom[0]<=latmax)))

lon_max_idx = np.max(np.where((lonmin<=longitude_dom[0]) & (longitude_dom[0]<=lonmax)))
lon_min_idx = np.min(np.where((lonmin<=longitude_dom[0]) & (longitude_dom[0]<=lonmax)))

lsm = lsm_regrid[lat_min_idx:lat_max_idx, lon_min_idx:lon_max_idx]


dates=daterange(datetime.date(min(Time_Dom_Datetime)), datetime.date(max(Time_Dom_Datetime)))
mean_dom=np.empty((len(list(dates)), pcp_dom.shape[1],  pcp_dom.shape[2]))
#mean_dom=[]
lsm=np.resize(lsm, mean_dom.shape)[0]

bad_values=np.ma.masked_array(pcp_dom,pcp_dom<0.)

for i, date in enumerate(dates):
    
    #Time_Dom_Datetime
    print date
    mean_dom[i]=np.mean(bad_values[np.where((mo==date.month) & (da==date.day))[0]],axis=0)

###############################

## Need to check time interval - TRMM data in 3 hourly intervals, given in mm/hr
############################

#time_interval=3

print mean_dom

dates_list=list(daterange(datetime.date(min(Time_Dom_Datetime)), datetime.date(max(Time_Dom_Datetime))))

#sum_dom = np.sum(pcp_dom*time_interval, axis=0)

mean=np.mean(mean_dom,axis=(1,2))

np.savez('/nfs/a90/eepdw/Data/Observations/Satellite/TRMM/trmm_daily_mean_EMBRACE_total.npz', dates=dates_list, means=mean)

#pdb.set_trace()

for ls in ['land', 'sea']:

    if ls=='land':

        lsm_weights = (100. - lsm)/100.

    if ls=='sea':

        lsm_weights = (lsm)/100.


    mean=np.average(mean_dom.reshape(mean_dom.shape[0],mean_dom.shape[1]*mean_dom.shape[2]), axis=1, weights=lsm_weights.flatten())

    print mean

    np.savez('/nfs/a90/eepdw/Data/Observations/Satellite/TRMM/trmm_daily_mean_EMBRACE_%s.npz' % ls, dates=dates_list, means=mean)

#pickle.dump([mean_dom, latitude_domsingle, longitude_domsingle, dates_list, hours], open('/nfs/a90/eepdw/Data/Saved_data/TRMM/trmm_emb_pcpmean_by_day.p', 'wb'))

#pickle.dump([sum_dom, latitude_domsingle, longitude_domsingle], open('/nfs/see-fs-01_users/eepdw/Saved_data/TRMM/trmm_emb_pcpsum_by_hour.p', 'wb'))
