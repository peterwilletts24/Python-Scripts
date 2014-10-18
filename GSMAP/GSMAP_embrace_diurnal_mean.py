import numpy as np

from netCDF4 import Dataset
from scipy.interpolate import griddata

from collections import defaultdict

lon_high = 101.866 
lon_low = 64.115

lat_high= 33.
lat_low=-6.79


numpy_cube=np.load('/nfs/a90/eepdw/Data/Observations/Satellite/GSMAP_Aug_Sep_2011/GSMAP_EMBRACE.npz')

# Load land sea mask.  TRMM land sea mask is in % of water coverage so 100% is all water

nc = Dataset('/nfs/a90/eepdw/Data/Observations/Satellite/TRMM/TMPA_mask.nc')

#  Regrid lsm to data grid (offset b 0.125 degrees

lsm_lons, lsm_lats = np.meshgrid(nc.variables['lon'][:],nc.variables['lat'][:])
lons_data, lats_data = np.meshgrid(numpy_cube['lons'], numpy_cube['lats'])
#lsm_regrid = griddata((lsm_lats.flatten(), lsm_lons.flatten()), nc.variables['landseamask'][:].flatten(), (lats_data,lons_data), method='linear')
lsm_regrid = griddata((lsm_lats.flatten(), lsm_lons.flatten()), nc.variables['landseamask'][:].flatten(), (lats_data,lons_data), method='linear')
####################################################
# Get min and max  index positions for latitude and longitude - FOR PCP

la_index_pcp = np.where((numpy_cube['lats']<=lat_high) & (numpy_cube['lats'] >= lat_low))
lo_index_pcp = np.where((numpy_cube['lons']<=lon_high) & (numpy_cube['lons'] >= lon_low))

la_i_max_pcp = np.max(la_index_pcp)
la_i_min_pcp = np.min(la_index_pcp)
lo_i_max_pcp = np.max(lo_index_pcp)
lo_i_min_pcp = np.min(lo_index_pcp)

print la_i_min_pcp,la_i_max_pcp, lo_i_min_pcp,lo_i_max_pcp

pcp_dom_2 = numpy_cube['rain_rates'][la_i_min_pcp:la_i_max_pcp, lo_i_min_pcp:lo_i_max_pcp, :]
lsm= lsm_regrid[la_i_min_pcp:la_i_max_pcp, lo_i_min_pcp:lo_i_max_pcp]
# Calculate mean for every time in the date range for entire area
#pdb.set_trace()

#nan_bad_values=np.where(numpy_cube['rain_rates']>=0, numpy_cube['rain_rates'], np.NaN)
bad_values=np.ma.masked_array(pcp_dom_2,pcp_dom_2<0.)

mean_of_each_time = bad_values.reshape((-1, bad_values.shape[2])).mean(axis=0)
#print pcp_dom_2.reshape(pcp_dom_2.shape[0], 1).shape

#time_hour = [datetime.utcfromtimestamp(v).hour for v in time_dom.flatten()]
time_hour = numpy_cube['time_list'][3]

mean_and_hour=zip(mean_of_each_time,time_hour)

# OCEAN - Calculate mean for every time in the date range

lsm_weights=lsm/100
#print pcp_dom_2.reshape(pcp_dom_2.shape[0], -1).shape
#print lsm_weights.flatten().shape

mean_oc = np.ma.average(bad_values.reshape((-1, bad_values.shape[2])), axis=0, weights=lsm_weights.flatten())
oc_mean_and_hour=zip(mean_oc,time_hour)

# LAND - Calculate mean for every time in the date range

lsm_weights=1-(lsm/100)
mean_la = np.ma.average(bad_values.reshape((-1, bad_values.shape[2])), axis=0, weights=lsm_weights.flatten())
la_mean_and_hour=zip(mean_la,time_hour)

#####################################################

#####################################################


# Sort into time of day #################

# Total

i = defaultdict(list)

for v,k in mean_and_hour:
    i[k].append(v)

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
    i[k].append(v)

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
    i[k].append(v)

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

np.savez("/nfs/a90/eepdw/Data/Observations/Satellite/GSMAP_Aug_Sep_2011/Diurnal/total_GSMAP_diurnal_average_rainfall_lat_%s_%s_lon_%s_%s" % (lat_low,lat_high, lon_low, lon_high), mean=mean, hour=hour )
np.savez("/nfs/a90/eepdw/Data/Observations/Satellite/GSMAP_Aug_Sep_2011/Diurnal/sea_GSMAP_diurnal_average_rainfall_lat_%s_%s_lon_%s_%s" % (lat_low,lat_high, lon_low, lon_high), mean=mean_o, hour=hour_o )
np.savez("/nfs/a90/eepdw/Data/Observations/Satellite/GSMAP_Aug_Sep_2011/Diurnal/land_GSMAP_diurnal_average_rainfall_lat_%s_%s_lon_%s_%s" % (lat_low,lat_high, lon_low, lon_high) , mean=mean_l, hour=hour_l )
