import cPickle as pickle
import numpy as np
from collections import defaultdict
from netCDF4 import Dataset
from scipy.interpolate import griddata

lon_high = 71 
lon_low = 67

lat_high= 28
lat_low=20

#lon_high = 116 
#lon_low = 30.5

#lat_high= 40
#lat_low=-11.25

pcp_dom, longitude_dom, latitude_dom, time_dom, time_hour = pickle.load(open('/nfs/a90/eepdw/Data/Saved_data/TRMM/trmm_emb_time_update_large.p', 'rb'))

# Load land sea mask.  TRMM land sea mask is in % of water coverage so 100% is all water

nc = Dataset('/nfs/a90/eepdw/Data/Observations/Satellite/TRMM/TMPA_mask.nc')

#  Regrid lsm to data grid (offset b 0.125 degrees

lsm_lons, lsm_lats = np.meshgrid(nc.variables['lon'][:],nc.variables['lat'][:])
lons_data, lats_data = np.meshgrid(longitude_dom[0], latitude_dom[0])
lsm_regrid = griddata((lsm_lats.flatten(), lsm_lons.flatten()), nc.variables['landseamask'][:].flatten(), (lats_data,lons_data), method='linear')


# Get min and max  index positions for latitude and longitude - FOR LSM



# Get min and max  index positions for latitude and longitude - FOR PCP

la_index_pcp = np.where((latitude_dom[0]<=lat_high) & (latitude_dom[0] >= lat_low))
lo_index_pcp = np.where((longitude_dom[0]<=lon_high) & (longitude_dom[0] >= lon_low))

la_i_max_pcp = np.max(la_index_pcp)
la_i_min_pcp = np.min(la_index_pcp)
lo_i_max_pcp = np.max(lo_index_pcp)
lo_i_min_pcp = np.min(lo_index_pcp)

print la_i_min_pcp,la_i_max_pcp, lo_i_min_pcp,lo_i_max_pcp

pcp_dom_2 = pcp_dom[:,la_i_min_pcp:la_i_max_pcp, lo_i_min_pcp:lo_i_max_pcp]
lsm= lsm_regrid[la_i_min_pcp:la_i_max_pcp, lo_i_min_pcp:lo_i_max_pcp]

print pcp_dom.shape
print pcp_dom_2.shape
print lsm.shape
####################################################

# Calculate mean for every time in the date range for entire area

mean_of_each_time = pcp_dom_2.reshape((pcp_dom_2.shape[0], -1)).mean(axis=1)
print pcp_dom_2.reshape(pcp_dom_2.shape[0], -1).shape
mean_and_hour=zip(mean_of_each_time,time_hour)

# OCEAN - Calculate mean for every time in the date range

lsm_weights=lsm/100
print pcp_dom_2.reshape(pcp_dom_2.shape[0], -1).shape
print lsm_weights.flatten().shape

mean_oc = np.ma.average(pcp_dom_2.reshape(pcp_dom_2.shape[0], -1), axis=1, weights=lsm_weights.flatten())
oc_mean_and_hour=zip(mean_oc,time_hour)

# LAND - Calculate mean for every time in the date range

lsm_weights=1-(lsm/100)
mean_la = np.ma.average(pcp_dom_2.reshape(pcp_dom_2.shape[0], -1), weights=lsm_weights.flatten(), axis=1)
la_mean_and_hour=zip(mean_la,time_hour)

#####################################################


# Sort into time of day #################

# Total

i = defaultdict(list)

for v,k in mean_and_hour:
    i[k.strip()].append(v)

mean=[]
hour=[]

# Average for each time of day

for q,a in enumerate(i.items()):
    #print a[1]
    if a[1]:
        mean.append(np.mean(a[1]))
        hour.append(a[0])

print mean
print hour

# Land

i = defaultdict(list)

for v,k in la_mean_and_hour:
    i[k.strip()].append(v)

mean_l=[]
hour_l=[]

# Average for each time of day

for q,a in enumerate(i.items()):
    #print a[1]
    if a[1]:
        mean_l.append(np.mean(a[1]))
        hour_l.append(a[0])

print mean_l
print hour_l

# Ocean

i = defaultdict(list)

for v,k in oc_mean_and_hour:
    i[k.strip()].append(v)

mean_o=[]
hour_o=[]

      # Average for each time of day

for q,a in enumerate(i.items()):
    #print a[1]
    if a[1]:
        mean_o.append(np.mean(a[1]))
        hour_o.append(a[0])

print mean_o
print hour_o

# Save

np.savez("/nfs/a90/eepdw/Data/Observations/Satellite/TRMM/Diurnal/total_trmm_diurnal_average_lat_%s_%s_lon_%s_%s_bit_above_western_ghats" % (lat_low,lat_high, lon_low, lon_high), mean=mean, hour=hour )
np.savez("/nfs/a90/eepdw/Data/Observations/Satellite/TRMM/Diurnal/sea_trmm_diurnal_average_lat_%s_%s_lon_%s_%s_bit_above_western_ghats" % (lat_low,lat_high, lon_low, lon_high), mean=mean_o, hour=hour_o )
np.savez("/nfs/a90/eepdw/Data/Observations/Satellite/TRMM/Diurnal/land_trmm_diurnal_average_lat_%s_%s_lon_%s_%s_bit_above_western_ghats" % (lat_low,lat_high, lon_low, lon_high), mean=mean_l, hour=hour_l )
